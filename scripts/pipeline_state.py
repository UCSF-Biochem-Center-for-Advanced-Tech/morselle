#!/usr/bin/env python3
"""Bind resumable pipeline outputs to their inputs and analysis configuration.

Large FASTQs/reference indexes are identified by absolute path, size and nanosecond
mtime rather than rereading them. This is a restart check, not a content checksum.
"""
import json
import os
from pathlib import Path
import sys


def identity(value):
    if not value:
        return value
    path = Path(value)
    if not path.exists():
        return value
    path = path.resolve()
    if path.is_dir():
        files = sorted(p for p in path.rglob("*") if p.is_file())
    else:
        files = [path]
    return {
        "path": str(path),
        "files": [[str(p), p.stat().st_size, p.stat().st_mtime_ns] for p in files],
    }


def bind(directory, values):
    root = Path(directory)
    manifest = root / ".pipeline-state.json"
    state = {"schema": 1, "configuration": [identity(v) for v in values]}
    if manifest.exists():
        if json.loads(manifest.read_text()) != state:
            raise ValueError(f"inputs or analysis settings changed for {root}; use a new output directory")
        return
    if root.exists() and any(p.is_file() for p in root.rglob("*")):
        raise ValueError(f"existing outputs in {root} have no restart configuration; use a new output directory")
    root.mkdir(parents=True, exist_ok=True)
    temporary = manifest.with_suffix(".partial")
    temporary.write_text(json.dumps(state, indent=2) + "\n")
    os.replace(temporary, manifest)


if __name__ == "__main__":
    try:
        bind(sys.argv[1], sys.argv[2:])
    except (ValueError, OSError) as error:
        sys.exit(f"ERROR: {error}")
