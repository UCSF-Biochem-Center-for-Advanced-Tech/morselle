#!/usr/bin/env bash
#
# sbx10x installer for Ubuntu (also works on most Debian-family systems).
#
# Installs build prerequisites, the Rust toolchain (if missing), builds the
# release binary, and installs it (plus the ARC whitelist, if available) under
# a prefix. Idempotent and safe to re-run.
#
# Usage:
#   ./install.sh [--prefix DIR] [--whitelist PATH] [--skip-apt] [--skip-rust]
#
# Options:
#   --prefix DIR      install root (default: /usr/local). Binary -> DIR/bin,
#                     whitelist -> DIR/share/sbx10x.
#   --whitelist PATH  path to 737K-arc-v1.txt.gz to install. If omitted, the
#                     copy bundled in resources/cellranger_barcodes/ is used
#                     when present.
#   --skip-apt        do not run apt-get (use when packages are already present
#                     or you lack sudo).
#   --skip-rust       do not install rustup (use the cargo already on PATH).
#   -h, --help        show this help.
#
set -euo pipefail

PREFIX="/usr/local"
WHITELIST=""
SKIP_APT=0
SKIP_RUST=0

log()  { printf '\033[1;36m[sbx10x]\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[sbx10x] WARN:\033[0m %s\n' "$*" >&2; }
die()  { printf '\033[1;31m[sbx10x] ERROR:\033[0m %s\n' "$*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix)    PREFIX="${2:?--prefix needs a value}"; shift 2 ;;
    --whitelist) WHITELIST="${2:?--whitelist needs a value}"; shift 2 ;;
    --skip-apt)  SKIP_APT=1; shift ;;
    --skip-rust) SKIP_RUST=1; shift ;;
    -h|--help)   sed -n '2,30p' "$0"; exit 0 ;;
    *)           die "unknown option: $1 (try --help)" ;;
  esac
done

# Repo root = directory of this script, so it can be run from anywhere.
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_DIR"
[[ -f Cargo.toml ]] || die "Cargo.toml not found in $REPO_DIR — run this from the sbx10x source tree."

# Pick sudo only when the install prefix is not writable.
writable_prefix() {
  local d="$PREFIX"
  while [[ ! -e "$d" && "$d" != "/" ]]; do d="$(dirname "$d")"; done
  [[ -w "$d" ]]
}
SUDO=""
if [[ $EUID -ne 0 ]] && ! writable_prefix; then
  command -v sudo >/dev/null 2>&1 || die "no write access to $PREFIX and sudo not found; use --prefix \$HOME/.local"
  SUDO="sudo"
fi

# 1. System build prerequisites -------------------------------------------------
if [[ $SKIP_APT -eq 0 ]] && command -v apt-get >/dev/null 2>&1; then
  APT_SUDO=""
  if [[ $EUID -ne 0 ]]; then
    command -v sudo >/dev/null 2>&1 || die "apt-get needs root; install prereqs manually or re-run with --skip-apt"
    APT_SUDO="sudo"
  fi
  log "installing build prerequisites via apt-get..."
  $APT_SUDO apt-get update -y
  $APT_SUDO apt-get install -y build-essential curl pkg-config git ca-certificates
else
  log "skipping apt-get (--skip-apt or apt-get not found)."
  command -v cc >/dev/null 2>&1 || warn "no C compiler found; install build-essential."
fi

# 2. Rust toolchain -------------------------------------------------------------
if ! command -v cargo >/dev/null 2>&1; then
  if [[ $SKIP_RUST -eq 1 ]]; then
    die "cargo not found and --skip-rust set. Install Rust (https://rustup.rs) first."
  fi
  log "installing Rust toolchain via rustup..."
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal
  # shellcheck disable=SC1091
  source "$HOME/.cargo/env"
fi
command -v cargo >/dev/null 2>&1 || die "cargo still not on PATH; open a new shell or 'source \$HOME/.cargo/env'."
log "using $(cargo --version)"

# 3. Build ----------------------------------------------------------------------
log "building release binary (this may take a few minutes the first time)..."
if [[ -f Cargo.lock ]]; then
  cargo build --release --locked
else
  cargo build --release
fi
[[ -x target/release/sbx10x ]] || die "build did not produce target/release/sbx10x"

# 4. Install binary -------------------------------------------------------------
BIN_DIR="$PREFIX/bin"
log "installing binary -> $BIN_DIR/sbx10x"
$SUDO install -d "$BIN_DIR"
$SUDO install -m 0755 target/release/sbx10x "$BIN_DIR/sbx10x"

# 5. Install whitelist (optional) ----------------------------------------------
SHARE_DIR="$PREFIX/share/sbx10x"
BUNDLED_WL="resources/cellranger_barcodes/737K-arc-v1.txt.gz"
WL_SRC=""
if [[ -n "$WHITELIST" ]]; then
  [[ -f "$WHITELIST" ]] || die "whitelist not found: $WHITELIST"
  WL_SRC="$WHITELIST"
elif [[ -f "$BUNDLED_WL" ]]; then
  WL_SRC="$BUNDLED_WL"
fi
if [[ -n "$WL_SRC" ]]; then
  log "installing whitelist -> $SHARE_DIR/737K-arc-v1.txt.gz"
  $SUDO install -d "$SHARE_DIR"
  $SUDO install -m 0644 "$WL_SRC" "$SHARE_DIR/737K-arc-v1.txt.gz"
else
  warn "no whitelist installed. The ARC barcode list (737K-arc-v1.txt.gz) is NOT"
  warn "in the repo. Obtain it from cellranger-arc-<ver>/lib/python/cellranger/barcodes/"
  warn "and pass it at runtime with --whitelist, or re-run with --whitelist PATH."
fi

# 6. Verify + smoke test --------------------------------------------------------
log "installed: $("$BIN_DIR/sbx10x" --version 2>/dev/null || echo '(version check failed)')"
SMOKE_FQ="tests/fixtures/smoke.fastq.gz"
SMOKE_SS="tests/fixtures/CE15204_bclconvert_input.csv"
if [[ -f "$SMOKE_FQ" && -f "$SMOKE_SS" ]]; then
  SMOKE_OUT="$(mktemp -d)"
  log "running smoke test on $SMOKE_FQ ..."
  if "$BIN_DIR/sbx10x" --fastqs "$SMOKE_FQ" -o "$SMOKE_OUT" --sample-sheet "$SMOKE_SS" \
       ${WL_SRC:+--whitelist "$WL_SRC"} >/dev/null 2>&1 \
     && [[ -f "$SMOKE_OUT/metrics_summary.json" ]]; then
    log "smoke test OK (metrics_summary.json written)."
  else
    warn "smoke test did not complete cleanly — check 'sbx10x --help' manually."
  fi
  rm -rf "$SMOKE_OUT"
fi

log "done. Ensure $BIN_DIR is on your PATH, then run: sbx10x --help"
