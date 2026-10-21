#!/usr/bin/env bash
#
# morselle installer for Ubuntu/Debian and macOS.
#
# Installs build prerequisites, minimap2 companions, the Rust toolchain (if
# missing), builds the release binaries, and installs them under a prefix.
# Idempotent and safe to re-run.
#
# Usage:
#   ./install.sh [--prefix DIR] [--skip-packages] [--skip-rust]
#
# Options:
#   --prefix DIR      install root (default: /usr/local). Binaries, k8, and
#                     paftools.js -> DIR/bin.
#   --skip-packages   do not run apt-get/Homebrew (packages must already exist).
#   --skip-apt        backward-compatible alias for --skip-packages.
#   --skip-rust       do not install rustup (use the cargo already on PATH).
#   -h, --help        show this help.
#
set -euo pipefail

PREFIX="/usr/local"
SKIP_PACKAGES=0
SKIP_RUST=0

log()  { printf '\033[1;36m[morselle]\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[morselle] WARN:\033[0m %s\n' "$*" >&2; }
die()  { printf '\033[1;31m[morselle] ERROR:\033[0m %s\n' "$*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix)    PREFIX="${2:?--prefix needs a value}"; shift 2 ;;
    --skip-packages|--skip-apt) SKIP_PACKAGES=1; shift ;;
    --skip-rust) SKIP_RUST=1; shift ;;
    -h|--help)   sed -n '2,21p' "$0"; exit 0 ;;
    *)           die "unknown option: $1 (try --help)" ;;
  esac
done

# Repo root = directory of this script, so it can be run from anywhere.
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_DIR"
[[ -f Cargo.toml ]] || die "Cargo.toml not found in $REPO_DIR — run this from the morselle source tree."

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

# 1. System packages ------------------------------------------------------------
OS="$(uname -s)"
if [[ $SKIP_PACKAGES -eq 1 ]]; then
  log "skipping system packages (--skip-packages)."
elif [[ "$OS" == "Linux" ]] && command -v apt-get >/dev/null 2>&1; then
  APT_SUDO=""
  if [[ $EUID -ne 0 ]]; then
    command -v sudo >/dev/null 2>&1 || die "apt-get needs root; install packages manually or re-run with --skip-packages"
    APT_SUDO="sudo"
  fi
  log "installing build prerequisites, minimap2, and samtools via apt-get..."
  $APT_SUDO apt-get update -y
  $APT_SUDO apt-get install -y \
    build-essential curl pkg-config git ca-certificates cmake clang libclang-dev \
    minimap2 samtools
elif [[ "$OS" == "Darwin" ]]; then
  command -v brew >/dev/null 2>&1 || die "Homebrew is required on macOS; install it from https://brew.sh or use --skip-packages"
  log "installing build prerequisites, minimap2, and samtools via Homebrew..."
  brew install cmake pkg-config minimap2 samtools
else
  die "unsupported package manager; install a C compiler, cmake, clang, minimap2, and samtools, then use --skip-packages"
fi

command -v cc >/dev/null 2>&1 || die "no C compiler found"
command -v minimap2 >/dev/null 2>&1 || die "minimap2 not found after package setup"
command -v samtools >/dev/null 2>&1 || die "samtools not found after package setup"
log "using minimap2 $(minimap2 --version)"
log "using samtools $(samtools --version-only)"

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
log "building release binaries (this may take a few minutes the first time)..."
if [[ -f Cargo.lock ]]; then
  cargo build --release --locked
else
  cargo build --release
fi
# Every private engine `morselle` dispatches to (engine() in src/bin/morselle.rs).
# 0.5.2-0.6.0 omitted the four Flex/CRISPR helpers here, so an installed
# morselle could not run probeset-align, probe-demux, guides or
# protospacer-call; tests/pipeline_scripts.rs now checks this list against
# the dispatch table.
HELPERS=(morselle-extract morselle-illumina morselle-cells morselle-annotate morselle-molecules bam-audit cr-bam-to-mtx
  morselle-probe-demux morselle-guides morselle-protospacer-call morselle-probeset-align)
for binary in morselle "${HELPERS[@]}"; do
  [[ -x target/release/$binary ]] || die "build did not produce target/release/$binary"
done

# 4. Install public CLI and private engines ------------------------------------
BIN_DIR="$PREFIX/bin"
LIBEXEC_DIR="$PREFIX/libexec/morselle"
log "installing public command -> $BIN_DIR/morselle"
$SUDO install -d "$BIN_DIR"
$SUDO install -m 0755 target/release/morselle "$BIN_DIR/morselle"
$SUDO install -d "$LIBEXEC_DIR"
for helper in "${HELPERS[@]}"; do
  log "installing private engine -> $LIBEXEC_DIR/$helper"
  $SUDO install -m 0755 "target/release/$helper" "$LIBEXEC_DIR/$helper"
done

# 5. Install minimap2 companion tools -------------------------------------------
# Package-manager minimap2 distributions generally install only minimap2 itself.
# paftools.js requires k8 (not Node), so install the official small upstream
# runtime and the paftools script matching the installed minimap2 release.
K8_VERSION="0.2.4"
K8_TMP="$(mktemp -d)"
trap 'rm -rf "$K8_TMP"' EXIT
K8_OS="$(uname -s)"
case "$K8_OS" in
  Linux|Darwin) ;;
  *) die "no official k8 binary is available for $K8_OS" ;;
esac
log "installing k8 $K8_VERSION for paftools.js..."
curl --fail --location --silent --show-error \
  "https://github.com/attractivechaos/k8/releases/download/v${K8_VERSION}/k8-${K8_VERSION}.tar.bz2" \
  | tar -xj -C "$K8_TMP"
K8_BIN="$K8_TMP/k8-${K8_VERSION}/k8-${K8_OS}"
[[ -x "$K8_BIN" ]] || die "k8 archive did not contain k8-${K8_OS}"
$SUDO install -m 0755 "$K8_BIN" "$BIN_DIR/k8"

MM2_VERSION="$(minimap2 --version 2>&1 | sed -E 's/^v?([0-9]+\.[0-9]+).*/\1/')"
[[ "$MM2_VERSION" =~ ^[0-9]+\.[0-9]+$ ]] || die "could not determine minimap2 version"
PAFTOOLS_URL="https://raw.githubusercontent.com/lh3/minimap2/v${MM2_VERSION}/misc/paftools.js"
log "installing paftools.js for minimap2 $MM2_VERSION..."
curl --fail --location --silent --show-error "$PAFTOOLS_URL" > "$K8_TMP/paftools.js"
[[ -s "$K8_TMP/paftools.js" ]] || die "downloaded paftools.js is empty"
$SUDO install -m 0755 "$K8_TMP/paftools.js" "$BIN_DIR/paftools.js"
# paftools intentionally exits nonzero when invoked without a command after
# printing its help. Capture that expected status explicitly: with `pipefail`,
# piping the invocation directly to grep incorrectly reports verification failure.
PAFTOOLS_HELP="$("$BIN_DIR/k8" "$BIN_DIR/paftools.js" 2>&1 || true)"
grep -qi '^Usage:' <<< "$PAFTOOLS_HELP" || die "paftools.js verification failed"
log "using k8 $K8_VERSION and paftools.js (minimap2 $MM2_VERSION)"

# 6. Install staged barcode resources (optional) --------------------------------
SHARE_DIR="$PREFIX/share/morselle"
BARCODE_SHARE_DIR="$SHARE_DIR/barcodes"
if compgen -G "resources/10x_barcodes/*" >/dev/null; then
  log "installing 10x barcode resources -> $BARCODE_SHARE_DIR/"
  $SUDO install -d "$BARCODE_SHARE_DIR"
  for barcode in resources/10x_barcodes/*; do
    [[ -f "$barcode" ]] || continue
    $SUDO install -m 0644 "$barcode" "$BARCODE_SHARE_DIR/$(basename "$barcode")"
  done
fi

# 7. Verify + smoke test --------------------------------------------------------
log "installed: $("$BIN_DIR/morselle" --version 2>/dev/null || echo '(version check failed)')"
log "installed: $(sed -n '1p' <<< "$PAFTOOLS_HELP")"
SMOKE_FQ="tests/fixtures/smoke.fastq.gz"
SMOKE_SS="tests/fixtures/CE15204_bclconvert_input.csv"
if [[ -f "$SMOKE_FQ" && -f "$SMOKE_SS" ]]; then
  SMOKE_OUT="$(mktemp -d)"
  log "running smoke test on $SMOKE_FQ ..."
  if "$BIN_DIR/morselle" extract --fastqs "$SMOKE_FQ" -o "$SMOKE_OUT" --sample-sheet "$SMOKE_SS" \
       >/dev/null 2>&1 \
     && [[ -f "$SMOKE_OUT/metrics_summary.json" ]]; then
    log "smoke test OK (metrics_summary.json written)."
  else
    warn "smoke test did not complete cleanly — check 'morselle extract --help' manually."
  fi
  rm -rf "$SMOKE_OUT"
fi

log "done. Ensure $BIN_DIR is on your PATH, then run: morselle --help"
