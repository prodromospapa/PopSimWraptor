#!/bin/sh
# -----------------------------------------------------------------------------
# install_dependencies.sh — robust POSIX installer with spinner UI (raisd-ai)
# -----------------------------------------------------------------------------
# Steps
#  1) Detect conda/mamba/micromamba
#  2) Create/Update env from ./raisd-ai.yml (override: YML_PATH=...)
#  3) Resolve env prefix & verify Java
#  4) Install msms wrapper
#  5) Install simulator wrappers (+ copy engines/ export/ for imports)
#  6) Install/Update RAiSD-AI (skipped if toolchain absent)
# -----------------------------------------------------------------------------

set -eu

START_TS="$(date +%s)"
YML="${YML_PATH-./raisd-ai.yml}"
MSMS_SHA256="${MSMS_SHA256-}"
RAISD_AI_REF="${RAISD_AI_REF-}"

# ---------- UI helpers ---------------------------------------------------------
is_tty() { [ -t 1 ] && [ -t 2 ]; }
if is_tty; then
  B="$(printf '\033[1m')"; RED="$(printf '\033[31m')"; GRN="$(printf '\033[32m')"; YEL="$(printf '\033[33m')"; CYA="$(printf '\033[36m')"; RST="$(printf '\033[0m')"
else
  B=""; RED=""; GRN=""; YEL=""; CYA=""; RST=""
fi
ok()   { printf "%s%s%s\n" "$GRN" "$1" "$RST"; }
warn() { printf "%s%s%s\n" "$YEL" "$1" "$RST"; }
err()  { printf "%s%s%s\n" "$RED" "$1" "$RST" >&2; }
info() { printf "%s%s%s\n" "$CYA" "$1" "$RST"; }

# step: supports either (n total title) or (title)
step() {
  if [ "$#" -ge 3 ]; then
    n="$1"; total="$2"; title="$3"
    printf "\n%s▶ %s [%s/%s]%s\n" "$B" "$title" "$n" "$total" "$RST"
  else
    title="$1"
    printf "\n%s▶ %s%s\n" "$B" "$title" "$RST"
  fi
}

elapsed() {
  end="$(date +%s)"
  mins=$(( (end - START_TS) / 60 )); secs=$(( (end - START_TS) % 60 ))
  printf "%02dm%02ds" "$mins" "$secs"
}

# Cursor control (best-effort)
_show_cursor() { command -v tput >/dev/null 2>&1 && tput cnorm 2>/dev/null || true; }
_hide_cursor() { command -v tput >/dev/null 2>&1 && tput civis 2>/dev/null || true; }
trap '_show_cursor' EXIT INT

# runq: spinner with last-40-lines-on-failure
runq() {
  title="$1"; shift
  tmpf="$(mktemp)"
  spinner_pid=""
  if is_tty; then
    _hide_cursor || true
    (
      i=0
      while :; do
        case $((i % 4)) in 0) f='|';; 1) f='/';; 2) f='-';; 3) f='\\';; esac
        printf "\r\033[K⏳ %s %s" "$title" "$f"
        i=$((i+1))
        sleep 0.1
      done
    ) &
    spinner_pid="$!"
  else
    printf "⏳ %s..." "$title"
  fi

  if ( "$@" >/dev/null 2>>"$tmpf" ); then
    if [ -n "${spinner_pid:-}" ]; then
      kill "$spinner_pid" 2>/dev/null || true
      wait "$spinner_pid" 2>/dev/null || true
      _show_cursor || true
      printf "\r\033[K✅ %s\n" "$title"
    else
      printf "\r✅ %s\n" "$title"
    fi
    rm -f "$tmpf"
    return 0
  else
    if [ -n "${spinner_pid:-}" ]; then
      kill "$spinner_pid" 2>/dev/null || true
      wait "$spinner_pid" 2>/dev/null || true
      _show_cursor || true
      printf "\r\033[K❌ %s\n" "$title"
    else
      printf "\r❌ %s\n" "$title"
    fi
    err "   Failed."
    echo
    warn "Last 40 stderr lines:"
    tail -n 40 "$tmpf" | sed 's/^/   /'
    rm -f "$tmpf"
    exit 1
  fi
}

# ---------- Status flags for summary ------------------------------------------
MSMS_STATUS="skipped"
SIM_STATUS="skipped"
RAISD_STATUS="skipped"

# ---------- Plan ---------------------------------------------------------------
TOTAL_STEPS=6
CUR=0

# 1) Detect conda
CUR=$((CUR+1)); step "$CUR" "$TOTAL_STEPS" "Checking for Conda/Mamba"
if ! command -v conda >/dev/null 2>&1 && ! command -v mamba >/dev/null 2>&1 && ! command -v micromamba >/dev/null 2>&1; then
  err "Conda/Mamba not installed."
  exit 1
fi
if command -v mamba >/dev/null 2>&1; then CONDA_EXE="$(command -v mamba)"
elif command -v micromamba >/dev/null 2>&1; then CONDA_EXE="$(command -v micromamba)"
else CONDA_EXE="$(command -v conda)"; fi
ok "Conda front-end: $CONDA_EXE"

# 2) Create/Update env
CUR=$((CUR+1)); step "$CUR" "$TOTAL_STEPS" "Updating conda env 'raisd-ai'"
[ -f "$YML" ] || { err "Environment file not found: $YML"; exit 1; }
if "$CONDA_EXE" env list | grep -Eq '^[^#]*\braisd-ai\b'; then
  runq "conda env update" "$CONDA_EXE" env update -f "$YML" --prune --quiet
else
  runq "conda env create" "$CONDA_EXE" env create -f "$YML" --quiet
fi

# 3) Resolve prefix & Java
CUR=$((CUR+1)); step "$CUR" "$TOTAL_STEPS" "Verifying Java (openjdk) in env"
PREFIX="$("$CONDA_EXE" run -n raisd-ai sh -lc 'printf %s "$CONDA_PREFIX"')"
[ -n "${PREFIX:-}" ] && [ -d "$PREFIX" ] || { err "Could not resolve CONDA_PREFIX."; exit 1; }
run_in_env() { "$CONDA_EXE" run -n raisd-ai "$@"; }
mkdir -p "$PREFIX/bin" "$PREFIX/share/raisd-ai-tools" "$PREFIX/share/msms"
TOOLS_DIR="$PREFIX/share/raisd-ai-tools"
BIN_DIR="$PREFIX/bin"
if run_in_env sh -lc 'command -v java >/dev/null 2>&1'; then ok "Java OK."
else err "'java' not available in env. Add 'openjdk' to your YAML."; exit 1; fi

# 4) Install msms
CUR=$((CUR+1)); step "$CUR" "$TOTAL_STEPS" "Installing msms"
install_msms() {
  jar="$PREFIX/share/msms/msms.jar"
  bin="$BIN_DIR/msms"
  if [ ! -f "$jar" ]; then
    runq "Downloading msms.jar" sh -c '
      set -eu
      url="https://www.mabs.at/fileadmin/user_upload/p_mabs/msms3.2rc-b163.jar"
      tmp="$(mktemp)"
      if command -v curl >/dev/null 2>&1; then curl -fsSL "$url" -o "$tmp"; else wget -q -O "$tmp" "$url"; fi
      mv "$tmp" "$1"; chmod 0644 "$1"
      # PK signature check
      head -c4 "$1" | od -An -t x1 | tr -d " \n" | grep -qi 504b0304
    ' sh "$jar"
    if [ -n "$MSMS_SHA256" ]; then
      runq "Verifying msms.jar SHA256" sh -c '
        set -eu
        if command -v sha256sum >/dev/null 2>&1; then got="$(sha256sum "$1" | awk "{print \$1}")"; else got="$(shasum -a 256 "$1" | awk "{print \$1}")"; fi
        [ "$got" = "$2" ]
      ' sh "$jar" "$MSMS_SHA256"
    fi
  else
    info "msms jar already present."
  fi
  # Wrapper
  runq "Installing msms wrapper" sh -c '
    set -eu
    cat > "$1" <<EOF
#!/bin/sh
set -eu
exec java -jar "$2" "\$@"
EOF
    chmod 0755 "$1"
  ' sh "$bin" "$jar"
  MSMS_STATUS="ok"
  ok "msms ready."
}
install_msms

# 5) Install simulator wrappers (+ copy engines/export)
CUR=$((CUR+1)); step "$CUR" "$TOTAL_STEPS" "Installing simulator wrapper"
install_simulator() {
  src_dir="$(cd "$(dirname "$0")" && pwd)"
  # Require simulator.py at source
  if [ ! -f "$src_dir/simulator.py" ]; then
    warn "simulator.py not found beside installer; skipping."
    return 0
  fi
  # Copy top-level .py files
  runq "Copying simulator module files" sh -c '
    set -eu
    src="$1"; dst="$2"
    mkdir -p "$dst"
    copied=0
    for f in "$src"/*.py; do
      [ -f "$f" ] || continue
      cp "$f" "$dst/"
      copied=1
    done
    [ "$copied" -eq 1 ] || cp "$src/simulator.py" "$dst/"
    chmod 755 "$dst/simulator.py" 2>/dev/null || true
  ' sh "$src_dir" "$TOOLS_DIR"
  # Copy engines/ and export/ recursively if present
  runq "Copying engines/ and export/ (if present)" sh -c '
    set -eu
    src="$1"; dst="$2"
    copy_pkg() {
      name="$1"
      if [ -d "$src/$name" ]; then
        mkdir -p "$dst/$name"
        (cd "$src/$name" && find . -type d -print | sed "s|^\./||" | while read -r d; do mkdir -p "$dst/$name/$d"; done)
        (cd "$src/$name" && find . -type f -print | sed "s|^\./||" | while read -r f; do cp "$src/$name/$f" "$dst/$name/$f"; done)
        [ -f "$dst/$name/__init__.py" ] || : > "$dst/$name/__init__.py"
      fi
    }
    copy_pkg engines
    copy_pkg export
  ' sh "$src_dir" "$TOOLS_DIR"
  # Create wrappers in bin (both names)
  runq "Installing simulator entry points" sh -c '
    set -eu
    bin="$1"; bin2="$2"
    cat > "$bin" << "EOF"
#!/bin/sh
set -eu
PREFIX="$(cd "$(dirname "$0")/.." && pwd)"
exec "$PREFIX/bin/python" "$PREFIX/share/raisd-ai-tools/simulator.py" "$@"
EOF
    chmod 0755 "$bin"
    cp "$bin" "$bin2"
    chmod 0755 "$bin2"
  ' sh "$BIN_DIR/simulator" "$BIN_DIR/simulator.py"
  SIM_STATUS="ok"
  ok "simulator installed to PATH."
}
install_simulator

# 6) RAiSD-AI (non-blocking)
CUR=$((CUR+1)); step "$CUR" "$TOTAL_STEPS" "Installing / updating RAiSD-AI"
install_raisd() {
  repo="$TOOLS_DIR/RAiSD-AI"
  if ! run_in_env sh -lc 'command -v git >/dev/null && command -v make >/dev/null && command -v gcc >/dev/null'; then
    info "Toolchain not in env; skipping RAiSD-AI."
    return 0
  fi
  runq "Preparing RAiSD-AI repo" sh -c '
    set -eu
    repo="$1"; url="https://github.com/alachins/RAiSD-AI.git"; conda="$2"; ref="$3"
    if [ ! -d "$repo/.git" ]; then
      rm -rf "$repo" || true
      "$conda" run -n raisd-ai git clone --depth 1 "$url" "$repo" >/dev/null 2>&1 || true
      if [ -n "$ref" ]; then
        "$conda" run -n raisd-ai git -C "$repo" fetch --tags >/dev/null 2>&1 || true
        "$conda" run -n raisd-ai git -C "$repo" checkout --quiet "$ref" >/dev/null 2>&1 || true
      fi
    else
      "$conda" run -n raisd-ai git -C "$repo" fetch --all --tags >/dev/null 2>&1 || true
      "$conda" run -n raisd-ai git -C "$repo" pull --ff-only >/dev/null 2>&1 || true
    fi
  ' sh "$repo" "$CONDA_EXE" "$RAISD_AI_REF"
  # Try to compile but NEVER fail the installer
  runq "Compiling RAiSD-AI" sh -c '
    set +e
    repo="$1"; conda="$2"
    [ -x "$repo/compile-RAiSD-AI.sh" ] || [ -x "$repo/compile-RAiSD-AI-ZLIB.sh" ] || exit 0
    "$conda" run -n raisd-ai sh -lc "cd \"$repo\" && [ -x ./compile-RAiSD-AI.sh ] && ./compile-RAiSD-AI.sh || true"
    "$conda" run -n raisd-ai sh -lc "cd \"$repo\" && [ -x ./compile-RAiSD-AI-ZLIB.sh ] && ./compile-RAiSD-AI-ZLIB.sh || true"
    exit 0
  ' sh "$repo" "$CONDA_EXE"
  # Install binaries if present
  runq "Installing RAiSD-AI binaries (if built)" sh -c '
    set -eu
    repo="$1"; bindir="$2"
    installed=0
    if [ -e "$repo/bin/release/RAiSD-AI" ]; then cp -L "$repo/bin/release/RAiSD-AI" "$bindir/RAiSD-AI"; chmod 755 "$bindir/RAiSD-AI"; installed=1; fi
    if [ -e "$repo/bin/release/RAiSD-AI-ZLIB" ]; then cp -L "$repo/bin/release/RAiSD-AI-ZLIB" "$bindir/RAiSD-AI-ZLIB"; chmod 755 "$bindir/RAiSD-AI-ZLIB"; installed=1; fi
    [ "$installed" -eq 0 ] && exit 0 || exit 0
  ' sh "$repo" "$BIN_DIR"
  if [ -x "$BIN_DIR/RAiSD-AI" ] || [ -x "$BIN_DIR/RAiSD-AI-ZLIB" ]; then RAISD_STATUS="ok"; fi
  info "RAiSD-AI step complete."
}
install_raisd || true

# ---------- Summary ------------------------------------------------------------
printf "\n%s %s (elapsed %s)\n" "🎉" "All tools installed into env." "$(elapsed)"
printf "\n%s" "┌──────────────── Installation Summary ────────────────┐"
printf "\n│ %-22s : %s" "msms" "$MSMS_STATUS"
printf "\n│ %-22s : %s" "simulator wrappers" "$SIM_STATUS"
printf "\n│ %-22s : %s" "RAiSD-AI" "$RAISD_STATUS"
printf "\n%s\n\n" "└──────────────────────────────────────────────────────┘"
info "Activate with: ${B}conda activate raisd-ai${RST}"
info "Run: ${B}simulator${RST} or ${B}simulator.py${RST} (avoid './simulator.py' to prefer env)"
