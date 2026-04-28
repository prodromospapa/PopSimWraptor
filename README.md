# PopSimWraptor

PopSimWraptor is a small simulation wrapper for generating population-genetics data from `stdpopsim` demographic models.
It provides a single CLI entry point that can run simulations with four backends:

- `msprime` — exact coalescent simulation
- `slim` — forward simulation with optional selective sweeps
- `msms` — command-based coalescent simulator with sweep support
- `discoal` — command-based coalescent simulator with sweep support

The repository is structured around the files that are actually used at runtime:

- `simulator.py` — main command-line interface
- `engines.py` — simulation backends and command builders
- `export.py` — output conversion helpers
- `install_dependencies.sh` — environment/setup helper
- `raisd-ai.yml` — conda environment definition

## What it does

The tool wraps `stdpopsim` species and demographic models and can export simulated data as:

- `ms`
- `ms.gz`
- `vcf`
- `vcf.gz`
- `bcf`
- site frequency spectrum (`.sfs.csv`)

It is intended for training and preparing data for RAiSD-AI-style workflows, but it is also useful as a general-purpose population simulation wrapper.

## Requirements

The recommended setup is the `raisd-ai` conda environment defined in `raisd-ai.yml`.
That environment includes the Python packages and external tools used by the repository, including:

- `stdpopsim`
- `msprime`
- `pyslim`
- `tskit`
- `bcftools`
- `vcftools`
- `openjdk`
- build tools needed for `discoal`

For `msms` and `discoal`, the corresponding command-line tools must be available on `PATH`.
The installer script can set these up for you.

## Installation

### Recommended: conda environment + helper script

1. Create the environment from `raisd-ai.yml`.
2. Run the installer script:

```bash
./install_dependencies.sh
```

The installer will:

- create or update the `raisd-ai` conda environment
- verify Java
- install the `msms` wrapper
- build and install `discoal` when a compiler toolchain is available
- copy the Python wrapper files into the environment prefix
- optionally install or update RAiSD-AI itself when the toolchain is available

After installation, the wrapper can be run as:

```bash
simulator
```

or directly:

```bash
simulator.py
```

### Manual use

If you already have the dependencies installed, you can run the script directly from the repository root:

```bash
python simulator.py --help
```

## Usage

`simulator.py` requires a species, engine, chromosome, demographic model, simulation populations, sample counts, and an output destination.

### Required arguments

- `--species` — stdpopsim species id or full name, for example `HomSap` or `Homo sapiens`
- `--engine` — one of `msprime`, `slim`, `msms`, `discoal`
- `--chromosome` — chromosome id from the selected species
- `--demography` — demographic model id from the selected species
- `--sim-population` — comma-separated simulation populations
- `--sample-counts` — comma-separated sample counts matching `--sim-population`
- `--output-file` — output file prefix or directory prefix depending on format

### Common optional arguments

- `--length` — chromosome length to simulate
- `--target-snps` — choose a length that aims for a target number of SNPs
- `--ref-population` — population used as the reference for metadata and some calculations
- `--simulations` — number of replicate simulations
- `--parallel` — number of worker processes
- `--output-format` — one of `ms`, `ms.gz`, `vcf`, `vcf.gz`, `bcf`
- `--sfs` — write the site frequency spectrum to `.sfs.csv`
- `--sfs-mean` — average SFS across replicates
- `--folded` — compute folded SFS
- `--sfs-normalized` — normalize the SFS
- `--get-commands` — print the external `msms`/`discoal` command that would be executed

### Sweep-related arguments

These are supported for `slim`, `msms`, and `discoal`:

- `--sweep-population`
- `--sweep-pos` — position as a proportion of chromosome length, between 0 and 1
- `--sweep-time` — `beginning` or a numeric time
- `--fixation-time`
- `--selection-coeff`

## Examples

### 1) Run an exact coalescent simulation with msprime

```bash
python simulator.py \
  --species HomSap \
  --engine msprime \
  --chromosome chr22 \
  --demography OutOfAfrica_2T12 \
  --sim-population AFR,EUR \
  --sample-counts 10,10 \
  --length 200000 \
  --simulations 10 \
  --output-format ms \
  --output-file results/homsap_msprime
```

### 2) Run SLiM with a selective sweep and export VCF

```bash
python simulator.py \
  --species HomSap \
  --engine slim \
  --chromosome chr22 \
  --demography OutOfAfrica_2T12 \
  --sim-population AFR,EUR \
  --sample-counts 10,10 \
  --length 200000 \
  --simulations 1 \
  --output-format vcf \
  --output-file results/homsap_slim \
  --sweep-population AFR \
  --sweep-pos 0.5 \
  --sweep-time beginning \
  --selection-coeff 0.1
```

### 3) Generate an msms command without executing downstream export

```bash
python simulator.py \
  --species HomSap \
  --engine msms \
  --chromosome chr22 \
  --demography OutOfAfrica_2T12 \
  --sim-population AFR,EUR \
  --sample-counts 10,10 \
  --length 200000 \
  --simulations 5 \
  --get-commands \
  --output-file results/homsap_msms
```

### 4) Compute the folded SFS

```bash
python simulator.py \
  --species HomSap \
  --engine msprime \
  --chromosome chr22 \
  --demography OutOfAfrica_2T12 \
  --sim-population AFR,EUR \
  --sample-counts 10,10 \
  --length 200000 \
  --simulations 20 \
  --sfs \
  --folded \
  --sfs-normalized \
  --output-file results/homsap_sfs
```

## Output files

Depending on the options you choose, the script will write:

- `OUTPUT.ms` or `OUTPUT.ms.gz`
- `OUTPUT.vcf`, `OUTPUT.vcf.gz`, or `OUTPUT.bcf`
- `OUTPUT.sfs.csv`

When writing VCF/BCF output for multiple simulations, the script creates an output directory and stores one file per replicate inside it.
Metadata describing the run is embedded in the output headers.

## Notes on engines

### `msprime`

Uses `stdpopsim`'s msprime engine and returns a tree sequence.
This is the most faithful backend for standard demographic models.

### `slim`

Uses `stdpopsim`'s SLiM engine and can simulate selective sweeps.

### `msms`

Builds an external `msms` command from the selected demographic model and can batch multiple replicates per job.

### `discoal`

Builds an external `discoal` command from the selected demographic model.
For time-varying migration, the implementation approximates the history with a constant rate derived from the model timeline.

## Supported model assumptions

The wrapper validates the selected species, chromosome, demographic model, and populations against `stdpopsim`.
It also checks that the requested sampling populations are valid for the chosen model.

Sweep simulations have extra restrictions enforced by the script, including valid sweep time, fixation time, and selection coefficient values.

## Tips

- Use comma-separated values for `--sim-population` and `--sample-counts`.
- `--target-snps` and `--length` are mutually exclusive.
- For VCF/BCF output, `bcftools` must be available.
- For `msms`/`discoal`, having the simulator binaries in the active environment is required.
- If you want a quick command preview, use `--get-commands` with `msms` or `discoal`.

## License

This repository is licensed under the MIT License. See `LICENSE` for the full text.
