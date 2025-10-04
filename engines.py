# conda install -c conda-forge slim=4.0
import inspect
import shlex
import tempfile
from typing import Dict, Iterable, Optional, Sequence
import stdpopsim as sps
import warnings

# Python < 3.10: swallow ignore_cleanup_errors so stdpopsim's SLiM engine works.
if "ignore_cleanup_errors" not in inspect.signature(tempfile.TemporaryDirectory.__init__).parameters:
    _orig_tempdir_init = tempfile.TemporaryDirectory.__init__

    def _patched_init(self, *args, ignore_cleanup_errors=False, **kwargs):
        _orig_tempdir_init(self, *args, **kwargs)

    tempfile.TemporaryDirectory.__init__ = _patched_init


warnings.filterwarnings("ignore")


def extended_events_define(contig, sweep_pos, sweep_population,sweep_time,fixation_time,selection_coeff):
    id = f"hard_sweep_{sweep_population}"
    contig.add_single_site(id=id, coordinate=sweep_pos)
    return sps.selective_sweep(
        single_site_id=id,
        population=sweep_population,
        selection_coeff=selection_coeff,
        mutation_generation_ago=sweep_time,
        end_generation_ago = fixation_time,
        min_freq_at_end=1
    )

def slim_simulate(species_std,model_std,chromosome,length,population_dict,slim_scaling_factor,slim_burn_in,sweep_population=None,sweep_pos=None,sweep_time=None,fixation_time=None,selection_coeff=None):
    contig = species_std.get_contig(chromosome,mutation_rate=model_std.mutation_rate,right=length)
    engine_std = sps.get_engine("slim")
    if sweep_population is not None and sweep_pos is not None and sweep_time is not None and fixation_time is not None:
        extended_events = extended_events_define(contig,sweep_pos,sweep_population,sweep_time,fixation_time,selection_coeff)
    else:
        extended_events = None
    ts =  engine_std.simulate(model_std,
                        contig,
                        population_dict,
                        extended_events=extended_events,
                        slim_scaling_factor=slim_scaling_factor,
                        slim_burn_in=slim_burn_in)
    t = ts.tables
    t.sequence_length = float(length)
    ts = t.tree_sequence()
    return ts


def msprime_simulation(species_std,model_std,chromosome,length,population_dict):
    contig = species_std.get_contig(chromosome,mutation_rate=model_std.mutation_rate,right=length)
    engine_std = sps.get_engine("msprime")
    ts =  engine_std.simulate(
            model_std,
            contig,
            population_dict)
    t = ts.tables
    t.sequence_length = float(length)
    ts = t.tree_sequence()
    return ts


def msms_command(
    species_std,
    model_std,
    chromosome,
    length,
    population_dict: Dict[str, int],
    sweep_population=None,
    sweep_pos=None,
    sweep_time=None,
    selection_coeff=None,
    *,
    replicates: int = 1,
    contig=None,
    pop_models: Optional[Sequence] = None,
    demo_dict: Optional[Dict] = None,
):
    def pop_idx(name):
        if isinstance(name, int):
            return max(0, min(len(pops) - 1, int(name)))
        for i, p in enumerate(pops):
            if p.name == name:
                return i
        return 0

    if contig is None:
        contig = species_std.get_contig(
            chromosome, mutation_rate=model_std.mutation_rate, right=length
        )

    if demo_dict is None:
        demo_dict = model_std.model.asdict()

    pops: Sequence = pop_models or tuple(model_std.populations)

    ploidy = int(getattr(species_std, "ploidy", 2))
    L = int(length or getattr(contig, "length", 1)) or 1

    counts = [int(population_dict.get(p.name, 0)) * ploidy for p in pops]
    if not any(counts) and pops:
        counts = [ploidy] + [0] * (len(pops) - 1)
    n = sum(counts) if counts else ploidy

    ref = next((i for i, c in enumerate(counts) if c > 0), 0)
    Ne0 = float(getattr(pops[ref], "initial_size", 1.0) or 1.0)
    fourNe = 2.0 * ploidy * Ne0

    mu = getattr(contig, "mutation_rate", getattr(model_std, "mutation_rate", 0.0)) or 0.0
    r = getattr(getattr(contig, "recombination_map", None), "mean_rate", None)
    if r is None:
        r = getattr(model_std, "recombination_rate", 0.0) or 0.0

    theta = fourNe * mu * L
    rho = fourNe * r * L

    replicates = max(1, int(replicates or 1))

    cmd = [
        "msms",
        str(n),
        str(replicates),
        "-t",
        f"{theta}",
        "-r",
        f"{rho}",
        str(L),
        "-I",
        str(len(pops) or 1),
    ]
    cmd += [str(c) for c in (counts if counts else [ploidy])]

    time_units = demo_dict.get(
        "time_units", getattr(model_std.model, "time_units", "generations")
    )
    if time_units == "years":
        gt = getattr(species_std, "generation_time", 1.0)
        upg = float(gt) if gt else 1.0
    else:
        upg = 1.0

    for i, p in enumerate(pops):
        if getattr(p, "initial_size", Ne0) != Ne0:
            cmd += ["-n", str(i + 1), f"{float(p.initial_size) / Ne0}"]
        g = float(getattr(p, "growth_rate", 0.0) or 0.0) * upg
        if g:
            cmd += ["-g", str(i + 1), f"{g * fourNe}"]

    migration_matrix = demo_dict.get("migration_matrix")
    if migration_matrix is not None:
        for i, row in enumerate(migration_matrix):
            for j, rate in enumerate(row):
                if i != j and rate:
                    cmd += [
                        "-m",
                        str(i + 1),
                        str(j + 1),
                        f"{float(rate) * upg * fourNe}",
                    ]

    events: Iterable[Dict] = list(
        sorted(demo_dict.get("events", []), key=lambda e: e.get("time", 0.0))
    )
    for e in events:
        t = (float(e.get("time", 0.0)) / upg) / max(fourNe, 1e-12)
        if "proportion" in e:
            cmd += [
                "-ej",
                f"{t}",
                str(pop_idx(e.get("source")) + 1),
                str(pop_idx(e.get("dest")) + 1),
            ]
        elif "matrix_index" in e:
            ij = e.get("matrix_index")
            rate = float(e.get("rate", 0.0)) * upg * fourNe
            if ij is None:
                cmd += ["-eM", f"{t}", f"{rate}"]
            else:
                cmd += [
                    "-em",
                    f"{t}",
                    str(pop_idx(ij[0]) + 1),
                    str(pop_idx(ij[1]) + 1),
                    f"{rate}",
                ]
        elif "initial_size" in e or "growth_rate" in e:
            i = pop_idx(e.get("population", e.get("population_id")))
            if e.get("initial_size") is not None:
                cmd += ["-en", f"{t}", str(i + 1), f"{float(e['initial_size']) / Ne0}"]
            if e.get("growth_rate") is not None:
                cmd += ["-eg", f"{t}", str(i + 1), f"{float(e['growth_rate']) * upg * fourNe}"]

    cmd += ["-oFP", "0.000000000000000000000000000000000"] # increase precision of floating point output

    if sweep_pos is not None and sweep_population is not None and selection_coeff is not None:
        loc = float(sweep_pos) / float(L)
        loc = 1e-6 if loc <= 0.0 else (1.0 - 1e-6 if loc >= 1.0 else loc)

        s = float(selection_coeff)
        h = 0.5
        S_AA = fourNe * s
        S_aA = fourNe * h * s

        cmd += [
            "-Sp",
            f"{loc}",
            "-SAA",
            f"{S_AA}",
            "-SAa",
            f"{S_aA}",
            "-Saa",
            "0",
            "-N",
            f"{Ne0}",
        ]

        t_sel = (
            0.0
            if sweep_time in (None, "beginning")
            else float(sweep_time) / max(fourNe, 1e-12)
        )
        denom = max(float(ploidy) * Ne0, 1.0)
        init_f = max(1.0 / denom, 1e-6)
        vec = ["0"] * max(len(pops), 1)
        vec[pop_idx(sweep_population)] = f"{init_f}"
        cmd += ["-SI", f"{t_sel}", str(len(vec)), *vec, "-Smark"]

        cmd += ["-SFC"]

    return " ".join(shlex.quote(str(x)) for x in cmd)



def discoal_command(
    species_std,
    model_std,
    chromosome,
    length,
    population_dict: Dict[str, int],
    sweep_population=None,
    sweep_pos=None,
    sweep_time=None,
    selection_coeff=None,
    replicates: int = 1,
    contig=None,
    pop_models: Optional[Sequence] = None,
    demo_dict: Optional[Dict] = None,
    growth_steps: int = 12,
):
    """
    Build a discoal command string from a stdpopsim demographic model.
    Notes (per discoal docs):
      - Multiple populations use `-p numPops size1 size2 ...` (zero-indexed IDs)
      - Times are in 4N0 units
      - Recombination: -r rho, where rho = 4*N0*r*(L-1)
      - Mutation: -t theta = 4*N0*u*L
      - Population splits: -ed time pop1 pop2 (merge backward in time)
      - Size changes: -en time popID sizeRel
      - Migration: only *constant* rates supported (`-m` and `-M`), time-varying migration not supported
      - We currently do *not* enable selection flags via discoal in this pipeline.
    """
    # If caller requested a selective sweep, delegate to the msms pipeline
    # which supports selection; we return the msms command string so the
    # same downstream pipeline can be used.
    if any(x is not None for x in (sweep_population, sweep_pos, sweep_time, selection_coeff)):
        return msms_command(
            species_std,
            model_std,
            chromosome,
            length,
            population_dict,
            sweep_population=sweep_population,
            sweep_pos=sweep_pos,
            sweep_time=sweep_time,
            selection_coeff=selection_coeff,
            replicates=replicates,
            contig=contig,
            pop_models=pop_models,
            demo_dict=demo_dict,
        )

    import math

    def pop_idx(name):
        if isinstance(name, int):
            return max(0, min(len(pops) - 1, int(name)))
        for i, p in enumerate(pops):
            if p.name == name:
                return i
        return 0

    if contig is None:
        contig = species_std.get_contig(
            chromosome, mutation_rate=model_std.mutation_rate, right=length
        )

    if demo_dict is None:
        demo_dict = model_std.model.asdict()

    pops: Sequence = pop_models or tuple(model_std.populations)

    ploidy = int(getattr(species_std, "ploidy", 2))
    L = int(length or getattr(contig, "length", 1)) or 1

    # sample sizes per population are HAPLOID counts
    counts = [int(population_dict.get(p.name, 0)) * ploidy for p in pops]
    if not any(counts) and pops:
        counts = [ploidy] + [0] * (len(pops) - 1)
    n = sum(counts) if counts else ploidy

    # Reference pop with first nonzero sample count
    ref = next((i for i, c in enumerate(counts) if c > 0), 0)
    Ne0 = float(getattr(pops[ref], "initial_size", 1.0) or 1.0)
    fourNe = 2.0 * ploidy * Ne0

    mu = getattr(contig, "mutation_rate", getattr(model_std, "mutation_rate", 0.0)) or 0.0
    r = getattr(getattr(contig, "recombination_map", None), "mean_rate", None)
    if r is None:
        r = getattr(model_std, "recombination_rate", 0.0) or 0.0

    # Discoal scaling (per msprime docs): theta = 4*N*u*L ; rho = 4*N*r*(L-1)
    theta = fourNe * mu * L
    rho = fourNe * r * max(L - 1, 1)

    replicates = max(1, int(replicates or 1))

    cmd = [
        "discoal",
        str(n),
        str(replicates),
        str(L),
        "-t",
        f"{theta}",
        "-r",
        f"{rho}",
    ]

    # Multiple populations declaration
    if len(pops) > 1:
        cmd += ["-p", str(len(pops))]
        cmd += [str(c) for c in counts]

    # Units: if model time units are years, convert to generations via species generation_time
    time_units = demo_dict.get(
        "time_units", getattr(model_std.model, "time_units", "generations")
    )
    if time_units == "years":
        gt = getattr(species_std, "generation_time", 1.0)
        upg = float(gt) if gt else 1.0
    else:
        upg = 1.0

    # Initial constant migration matrix (only present-time supported in discoal; no time-varying rates)
    migration_matrix = demo_dict.get("migration_matrix")
    if migration_matrix is not None:
        # If any later migration change events exist, error out (not supported by discoal)
        events: Iterable[Dict] = demo_dict.get("events", [])
        mig_change = any("matrix_index" in e and float(e.get("time", 0.0)) > 0.0 for e in events)
        if mig_change:
            raise ValueError("discoal backend does not support time-varying migration rates; adjust demography or use msms/msprime.")
        for i, row in enumerate(migration_matrix):
            for j, rate in enumerate(row):
                if i != j and rate:
                    cmd += ["-m", str(i), str(j), f"{float(rate) * upg * fourNe}"]

    # Build growth approximation timeline and apply -en size changes
    events: Iterable[Dict] = list(
        sorted(demo_dict.get("events", []), key=lambda e: e.get("time", 0.0))
    )

    # helper: add -en at time t (in generations) for population i with absolute size Ne_abs
    def add_en(t_generations: float, pop_index: int, Ne_abs: float):
        t = (float(t_generations) / upg) / max(fourNe, 1e-12)
        cmd.extend(["-en", f"{t}", str(pop_index), f"{Ne_abs / Ne0}"])

    # Track per-pop sizes and exponential growth rates (per generation)
    current_sizes = {i: float(getattr(p, "initial_size", Ne0) or Ne0) for i, p in enumerate(pops)}
    current_growth = {i: float(getattr(p, "growth_rate", 0.0) or 0.0) * upg for i, p in enumerate(pops)}
    last_t = 0.0

    # Create a timeline: include all event times and add a grid between them for growth approximation
    change_times = sorted(set([float(e.get("time", 0.0)) for e in events] + [0.0]))
    grid = []
    for a, b in zip(change_times, change_times[1:] + [max(change_times[-1], 0.0)]):
        if b <= a:
            continue
        if growth_steps > 1:
            step = (b - a) / float(growth_steps)
            grid.extend([a + k * step for k in range(1, growth_steps)])
    timeline = sorted(set(change_times + grid))

    # If any populations have initial sizes != Ne0, emit an -en at time 0
    for i, p in enumerate(pops):
        if getattr(p, "initial_size", Ne0) != Ne0:
            add_en(0.0, i, float(getattr(p, "initial_size", Ne0)))

    # Walk the timeline, apply growth and events
    for t in timeline:
        dt = t - last_t
        if dt > 0.0:
            for i in range(len(pops)):
                g = current_growth[i]
                if g:
                    current_sizes[i] = current_sizes[i] * (2.718281828459045 ** (g * dt))
                    add_en(t, i, current_sizes[i])
        # apply discrete events at exact t
        for e in [e for e in events if float(e.get("time", 0.0)) == t]:
            if "proportion" in e:
                # stdpopsim 'mass migration' with fraction!=1 is NOT supported here (no direct -ea mapping)
                fraction = float(e.get("proportion", 0.0))
                if not math.isclose(fraction, 1.0, rel_tol=1e-9, abs_tol=0.0):
                    raise ValueError(
                        "discoal backend cannot reproduce partial mass-migration "
                        f"(proportion={fraction}); use msprime or adjust the demography."
                    )
                src = pop_idx(e.get("source"))
                dst = pop_idx(e.get("dest"))
                # population split backward (merge src into dst)
                t_coal = (t / upg) / max(fourNe, 1e-12)
                cmd += ["-ed", f"{t_coal}", str(src), str(dst)]
            elif "matrix_index" in e:
                # Already handled above: time-varying migration not supported
                if float(e.get("time", 0.0)) > 0.0:
                    raise ValueError("discoal backend does not support time-varying migration rates (`-em`-like events); use msms/msprime.")
            elif "initial_size" in e or "growth_rate" in e:
                i = pop_idx(e.get("population", e.get("population_id")))
                if e.get("initial_size") is not None:
                    current_sizes[i] = float(e["initial_size"])
                    add_en(t, i, current_sizes[i])
                if e.get("growth_rate") is not None:
                    current_growth[i] = float(e["growth_rate"]) * upg
        last_t = t

    return " ".join(shlex.quote(str(x)) for x in cmd)
