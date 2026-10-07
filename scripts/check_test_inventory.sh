#!/usr/bin/env bash
# Fail if any test listed in tests/test_inventory.txt no longer exists.
#
# The 0.5.7 modularization silently deleted two molecules tests (one of
# them the only guard on --tagseq semantics) while "refactoring". Removing
# or renaming a test must now be an explicit, reviewable edit of
# tests/test_inventory.txt in the same commit.
#
#   scripts/check_test_inventory.sh            # check (CI, scripts/check.sh)
#   scripts/check_test_inventory.sh --update   # rewrite after adding/renaming tests
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
INVENTORY=tests/test_inventory.txt

# "<test binary> <test name>" for every test in every target.
current=$(cargo test --release -- --list --format terse 2>&1 \
    | awk '/Running .*\(target/ { split($2, p, "/"); bin = p[length(p)]; sub(/\.rs$/, "", bin) }
           / Running unittests / { bin = $3 }
           /: test$/ { sub(/: test$/, ""); print bin " " $0 }' \
    | sort)

if [[ ${1:-} == --update ]]; then
    printf '%s\n' "$current" > "$INVENTORY"
    echo "wrote $(wc -l < "$INVENTORY") tests to $INVENTORY"
    exit 0
fi

[[ -s $INVENTORY ]] || { echo "missing $INVENTORY; run $0 --update" >&2; exit 1; }
missing=$(comm -23 <(sort "$INVENTORY") <(printf '%s\n' "$current"))
added=$(comm -13 <(sort "$INVENTORY") <(printf '%s\n' "$current"))
if [[ -n $missing ]]; then
    echo "ERROR: tests removed or renamed without updating $INVENTORY:" >&2
    echo "$missing" | sed 's/^/  /' >&2
    echo "If intentional, run $0 --update and explain the removal in the commit." >&2
    exit 1
fi
if [[ -n $added ]]; then
    echo "note: $(wc -l <<< "$added") new test(s) not yet in $INVENTORY (run $0 --update):"
    echo "$added" | sed 's/^/  /'
fi
echo "test inventory OK ($(wc -l <<< "$current") tests)"
