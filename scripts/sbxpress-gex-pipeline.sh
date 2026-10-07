#!/usr/bin/env bash
# sbxpress-gex-pipeline.sh — Forwarder wrapper for sbxpress-pipeline.sh (GEX/SBX mode)
# Preserved for backward compatibility.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$SCRIPT_DIR/sbxpress-pipeline.sh" --input-mode sbx --analysis-mode gex "$@"
