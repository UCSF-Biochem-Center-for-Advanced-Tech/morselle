#!/usr/bin/env bash
# Helpers shared by the unified and standalone GEX pipeline drivers.
# The caller defines SCRIPT_DIR, die(), and (optionally) dry-run mode.

bind_pipeline_state() {
    [[ ${MORSELLE_PIPELINE_DRY_RUN:-0} == 1 ]] && return 0
    command -v python3 >/dev/null 2>&1 || die "python3 is required for safe pipeline restarts"
    python3 "$SCRIPT_DIR/pipeline_state.py" "$@"
}

prepare_matrix_stage() {
    [[ ${MORSELLE_PIPELINE_DRY_RUN:-0} == 1 ]] && return 0
    rm -rf "$1"
    mkdir -p "$1"
}

# Publish a successfully completed matrix/cells stage. A .done marker is
# written last so an interrupted publication will be rebuilt on restart.
publish_matrix_stage() {
    local stage=$1 destination=$2 marker=$3 item name
    [[ ${MORSELLE_PIPELINE_DRY_RUN:-0} == 1 ]] && return 0
    for item in "$stage"/*; do
        [[ -e $item ]] || continue
        name=$(basename "$item")
        # Only replace the artifacts produced by this completed stage.
        rm -rf "$destination/$name"
        mv "$item" "$destination/$name"
    done
    rmdir "$stage"
    touch "$destination/$marker"
}
