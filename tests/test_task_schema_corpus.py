"""Corpus safety gate for strict task-card validation."""

from __future__ import annotations

from pathlib import Path
import unittest

from jarvishep2.task_config import load_task_yaml
from jarvishep2.task_validation import validate_task_config


class TaskSchemaCorpusTests(unittest.TestCase):
    def test_example_cards_have_no_schema_vocabulary_rejections(self) -> None:
        cards = sorted(Path(__file__).parents[2].glob("Jarvis-Examples/*/bin/*.yaml"))
        self.assertGreaterEqual(len(cards), 65)
        bad: list[str] = []
        for card in cards:
            report = validate_task_config(load_task_yaml(str(card)))
            schema_errors = [item.format_line() for item in report.errors() if item.code.startswith("JV2-SCH")]
            if schema_errors:
                bad.append(f"{card}:\n" + "\n".join(schema_errors))
        self.assertEqual(bad, [])

    def test_unmigrated_example_methods_are_runtime_errors_not_schema_errors(self) -> None:
        cards = sorted(Path(__file__).parents[2].glob("Jarvis-Examples/*/bin/*.yaml"))
        expected = {"DNN", "DREAM", "ESS", "HMC", "MALA", "NUTS", "RLTPMCMC", "RobustAM", "SliceMCMC"}
        seen: set[str] = set()
        for card in cards:
            report = validate_task_config(load_task_yaml(str(card)))
            for item in report.errors():
                if item.code == "JV2-MTH-003":
                    seen.add(str(load_task_yaml(str(card))["Sampling"]["Method"]))
        self.assertEqual(seen, expected)


if __name__ == "__main__":
    unittest.main()
