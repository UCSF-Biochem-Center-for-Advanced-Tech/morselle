#!/usr/bin/env python3
"""Pipeline restart regressions using only the standard library."""
import ast
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))
from pipeline_state import bind, identity


class ScriptRegressions(unittest.TestCase):
    def test_empty_configuration_values_are_not_directory_inputs(self):
        self.assertEqual(identity(""), "")

    def test_restart_requires_identical_inputs_and_analysis_settings(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.fastq"
            source.write_text("original")
            output = root / "output"
            bind(output, [str(source), "shared"])
            bind(output, [str(source), "shared"])
            with self.assertRaisesRegex(ValueError, "settings changed"):
                bind(output, [str(source), "fractional"])
            source.write_text("changed input length")
            with self.assertRaisesRegex(ValueError, "settings changed"):
                bind(output, [str(source), "shared"])

    def test_existing_outputs_without_manifest_are_not_silently_adopted(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "raw.h5").write_text("unfinished")
            with self.assertRaisesRegex(ValueError, "no restart configuration"):
                bind(root, ["shared"])


if __name__ == "__main__":
    unittest.main()
