#!/usr/bin/env python3
"""Pipeline restart and reporting regressions using only the standard library."""
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
from presentation_reporting import load_metric, metric_rows, export_insert_sizes


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

    def test_report_discovers_actual_pipeline_names_and_uses_molecule_cell_values(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sample = root / "S1"
            (root / "extract").mkdir()
            (sample / "molecules").mkdir(parents=True)
            (sample / "outs").mkdir()
            mol = {"molecules_counted": 37, "representatives_tagged": 38}
            cells = {"filtered_cells": 4, "filtered_molecules": 29}
            illumina = {"insert_length_histogram": {"120": 3}}
            (sample / "molecules/S1.molecules.metrics.json").write_text(json.dumps(mol))
            (sample / "outs/cell_calling_metrics.json").write_text(json.dumps(cells))
            (root / "extract/S1.illumina.metrics.json").write_text(json.dumps(illumina))
            self.assertEqual(load_metric(None, str(sample), "illumina"), illumina)
            rows = metric_rows(load_metric(None, str(sample), "molecules"), load_metric(None, str(sample), "cell"))
            self.assertIn(["Counted molecules", "37"], rows)
            self.assertIn(["Called cells/barcodes", "4"], rows)
            self.assertIn(["Filtered molecules", "29"], rows)
            path = root / "insert.tsv"
            export_insert_sizes(illumina, path)
            self.assertEqual(path.read_text(), "insert_length\tcount\n120\t3\n")

    def test_report_rejects_missing_or_ambiguous_metrics(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaises(FileNotFoundError):
                load_metric(str(root / "missing.json"), None, "cell")
            for sample in ["A", "B"]:
                (root / sample).mkdir()
                (root / sample / "cell_calling_metrics.json").write_text("{}")
            with self.assertRaisesRegex(ValueError, "multiple cell metrics"):
                load_metric(None, root, "cell")

    def test_reporting_sources_parse_and_label_historical_examples(self):
        for name in ["build_presentation", "presentation_layout", "presentation_history", "presentation_reporting"]:
            ast.parse((SCRIPTS / f"{name}.py").read_text())
        self.assertIn("Historical mapping example", (SCRIPTS / "build_presentation.py").read_text())


if __name__ == "__main__":
    unittest.main()
