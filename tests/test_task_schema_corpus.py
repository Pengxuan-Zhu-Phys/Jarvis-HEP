"""Corpus safety gate for strict task-card validation."""

from __future__ import annotations

from pathlib import Path
import unittest

from jarvishep2.task_config import load_task_yaml
from jarvishep2.task_validation import validate_task_card_encoding, validate_task_config


EXAMPLES_ROOT = Path(__file__).parents[2] / "Jarvis-Examples"


class TaskSchemaCorpusTests(unittest.TestCase):
    @unittest.skipUnless(EXAMPLES_ROOT.is_dir(), "requires sibling Jarvis-Examples checkout")
    def test_example_cards_only_need_the_documented_utils_removal(self) -> None:
        cards = sorted(EXAMPLES_ROOT.glob("*/bin/*.yaml"))
        self.assertGreaterEqual(len(cards), 65)
        bad: list[str] = []
        legacy_utils_cards: list[str] = []
        for card in cards:
            report = validate_task_config(load_task_yaml(str(card)))
            schema_errors = [item for item in report.errors() if item.code.startswith("JV2-SCH")]
            if not schema_errors:
                continue
            if (
                len(schema_errors) == 1
                and schema_errors[0].path == "$"
                and "'Utils' was unexpected" in schema_errors[0].message
                and "moved to Jarvis-Operas" in (schema_errors[0].suggestion or "")
            ):
                legacy_utils_cards.append(str(card))
                continue
            bad.append(f"{card}:\n" + "\n".join(item.format_line() for item in schema_errors))
        self.assertEqual(bad, [])
        self.assertGreater(legacy_utils_cards, [])

    def test_unmigrated_example_methods_are_runtime_errors_not_schema_errors(self) -> None:
        if not EXAMPLES_ROOT.is_dir():
            self.skipTest("requires sibling Jarvis-Examples checkout")
        cards = sorted(EXAMPLES_ROOT.glob("*/bin/*.yaml"))
        expected = {"DNN", "DREAM", "ESS", "HMC", "MALA", "NUTS", "RLTPMCMC", "RobustAM", "SliceMCMC"}
        seen: set[str] = set()
        for card in cards:
            report = validate_task_config(load_task_yaml(str(card)))
            for item in report.errors():
                if item.code == "JV2-MTH-003":
                    seen.add(str(load_task_yaml(str(card))["Sampling"]["Method"]))
        self.assertEqual(seen, expected)

    @unittest.skipUnless(EXAMPLES_ROOT.is_dir(), "requires sibling Jarvis-Examples checkout")
    def test_example_cards_have_no_multimode_dependency_regressions(self) -> None:
        cards = sorted(EXAMPLES_ROOT.glob("*/bin/*.yaml"))
        self.assertGreaterEqual(len(cards), 65)
        bad: list[str] = []
        for card in cards:
            report = validate_task_config(load_task_yaml(str(card)))
            mode_errors = [
                item for item in report.errors() if item.code.startswith("JV2-MOD-")
            ]
            if mode_errors:
                bad.append(
                    f"{card}:\n"
                    + "\n".join(item.format_line() for item in mode_errors)
                )
        self.assertEqual(bad, [])

    def test_shipped_yaml_corpus_is_ascii_only(self) -> None:
        template_root = Path(__file__).parents[1] / "jarvishep2" / "project_template"
        parity_root = Path(__file__).parent / "parity_project"
        cards = [
            *EXAMPLES_ROOT.glob("*/bin/*.yaml"),
            *template_root.rglob("*.yaml"),
            *parity_root.glob("*.yaml"),
        ]
        self.assertGreaterEqual(len(cards), 78 if EXAMPLES_ROOT.is_dir() else 13)
        bad = []
        for card in cards:
            import yaml
            document = yaml.safe_load(card.read_text(encoding="utf-8"))
            issues = validate_task_card_encoding(document)
            if issues:
                bad.append(f"{card}: {[issue.path for issue in issues]}")
        self.assertEqual(bad, [])


if __name__ == "__main__":
    unittest.main()
