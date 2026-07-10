#!/usr/bin/env python3
"""Assert that a sliced-BAM reconstruction is a non-overcounting matrix subset."""

import gzip
import sys


def read(path):
    with gzip.open(path, "rt") as handle:
        line = next(handle)
        while line.startswith("%"):
            line = next(handle)
        shape = tuple(map(int, line.split()))
        values = {}
        for line in handle:
            row, column, value = map(int, line.split())
            values[row, column] = value
    return shape, values


expected_shape, expected = read(sys.argv[1])
actual_shape, actual = read(sys.argv[2])
assert actual_shape[:2] == expected_shape[:2], (actual_shape, expected_shape)
overcounts = {key: value for key, value in actual.items() if value > expected.get(key, 0)}
assert not overcounts, f"unexpected/over-counted entries: {list(overcounts.items())[:10]}"
recovered = sum(actual.values())
total = sum(expected.values())
print(f"recovered {recovered}/{total} expected molecules ({recovered / total:.3%}); no overcounts")
