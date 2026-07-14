#!/usr/bin/env bash
#
# Remove an sbxpress installation created by install.sh.
#
# Usage:
#   ./uninstall.sh [--prefix DIR] [--dry-run] [--keep-data] [--remove-companions] [--remove-legacy]
#
# Options:
#   --prefix DIR          installation root used with install.sh (default: /usr/local)
#   --dry-run             print what would be removed without changing anything
#   --keep-data           preserve DIR/share/sbxpress (including the whitelist)
#   --remove-companions   also remove k8 and paftools.js from DIR/bin
#   --remove-legacy       remove old sbx10x-era public app binaries
#   -h, --help            show this help
#
set -euo pipefail

PREFIX="/usr/local"
DRY_RUN=0
KEEP_DATA=0
REMOVE_COMPANIONS=0
REMOVE_LEGACY=0

log() { printf '\033[1;36m[sbxpress]\033[0m %s\n' "$*"; }
die() { printf '\033[1;31m[sbxpress] ERROR:\033[0m %s\n' "$*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix) PREFIX="${2:?--prefix needs a value}"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    --keep-data) KEEP_DATA=1; shift ;;
    --remove-companions) REMOVE_COMPANIONS=1; shift ;;
    --remove-legacy) REMOVE_LEGACY=1; shift ;;
    -h|--help) sed -n '2,15p' "$0"; exit 0 ;;
    *) die "unknown option: $1 (try --help)" ;;
  esac
done

[[ -n "$PREFIX" ]] || die "installation prefix cannot be empty"
[[ "$PREFIX" = /* ]] || die "installation prefix must be an absolute path: $PREFIX"
[[ "$PREFIX" != "/" ]] || die "refusing to use / as the installation prefix"

SUDO=""
if [[ $EUID -ne 0 && -e "$PREFIX" && ! -w "$PREFIX" ]]; then
  command -v sudo >/dev/null 2>&1 || die "no write access to $PREFIX and sudo not found"
  SUDO="sudo"
fi

remove_file() {
  local path=$1
  if [[ -e "$path" || -L "$path" ]]; then
    if [[ $DRY_RUN -eq 1 ]]; then
      log "would remove $path"
    else
      log "removing $path"
      $SUDO rm -f -- "$path"
    fi
  else
    log "not installed: $path"
  fi
}

remove_tree() {
  local path=$1
  if [[ -d "$path" ]]; then
    if [[ $DRY_RUN -eq 1 ]]; then
      log "would remove $path/"
    else
      log "removing $path/"
      $SUDO rm -rf -- "$path"
    fi
  else
    log "not installed: $path/"
  fi
}

remove_file "$PREFIX/bin/sbxpress"
remove_tree "$PREFIX/libexec/sbxpress"

if [[ $KEEP_DATA -eq 0 ]]; then
  remove_tree "$PREFIX/share/sbxpress"
else
  log "preserving $PREFIX/share/sbxpress/ (--keep-data)"
fi

if [[ $REMOVE_COMPANIONS -eq 1 ]]; then
  remove_file "$PREFIX/bin/k8"
  remove_file "$PREFIX/bin/paftools.js"
else
  log "preserving shared k8 and paftools.js (use --remove-companions to remove)"
fi

if [[ $REMOVE_LEGACY -eq 1 ]]; then
  for legacy in sbx10x sbx10x-annotate sbx10x-molecules bam-audit cr-bam-to-mtx; do
    remove_file "$PREFIX/bin/$legacy"
  done
fi

if [[ $DRY_RUN -eq 1 ]]; then
  log "dry run complete; nothing was removed"
else
  log "uninstall complete"
fi
