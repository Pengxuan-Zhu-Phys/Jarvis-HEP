"""Corpus safety gate for strict task-card validation."""

from __future__ import annotations

from pathlib import Path
import unittest

from jarvishep2.task_config import TaskCardLoadError, load_task_yaml
from jarvishep2.task_validation import validate_task_card_encoding, validate_task_config


EXAMPLES_ROOT = Path(__file__).parents[2] / "Jarvis-Examples"

# Existing external V1 card repeats Bounds.seed (lines 37/40 in the current
# checkout). It must now fail loading rather than silently overwrite a value.
# Do not waive duplicate-key diagnostics for any other corpus file.
_KNOWN_DUPLICATE_EXAMPLE = "Eggbox-V1/bin/Example_MultiNest_Operas.yaml"

_LEGACY_PTMCMC_BOUNDS_MARKERS = (
    "'proposal_scale' is a required property",
    "'temperature_ladder' is a required property",
    "'proposal_scales' was unexpected",
)


def _is_legacy_utils_error(item) -> bool:
    return (
        item.path == "$"
        and "'Utils' was unexpected" in item.message
        and "moved to Jarvis-Operas" in (item.suggestion or "")
    )


def _is_legacy_ptmcmc_bounds_error(item) -> bool:
    """Eggbox-V1 PTMCMC cards still use Bounds.proposal_scales and omit the ladder.

    Runtime accepts ``proposal_scales`` as an alias; the closed V2 schema key is
    ``proposal_scale`` plus ``temperature_ladder``.
    """
    if item.path != "$.Sampling.Bounds":
        return False
    return any(marker in item.message for marker in _LEGACY_PTMCMC_BOUNDS_MARKERS)


class TaskSchemaCorpusTests(unittest.TestCase):
    def _load_external_card(self, card: Path):
        try:
            return load_task_yaml(str(card))
        except TaskCardLoadError as exc:
            if exc.code != "JV2-YAML-003" or card.relative_to(EXAMPLES_ROOT).as_posix() != _KNOWN_DUPLICATE_EXAMPLE:
                raise
            self.assertEqual(exc.path, "$.Sampling.Bounds.seed")
            self.assertEqual(exc.location["file"], str(card))
            self.assertEqual(exc.related_locations[0]["file"], str(card))
            self.assertLess(exc.related_locations[0]["line"], exc.location["line"])
            self.assertIn("duplicate YAML key 'seed'", str(exc))
            return None

    @unittest.skipUnless(EXAMPLES_ROOT.is_dir(), "requires sibling Jarvis-Examples checkout")
    def test_example_cards_only_need_the_documented_utils_removal(self) -> None:
        cards = sorted(EXAMPLES_ROOT.glob("*/bin/*.yaml"))
        self.assertGreaterEqual(len(cards), 65)
        bad: list[str] = []
        legacy_utils_cards: list[str] = []
        for card in cards:
            config = self._load_external_card(card)
            if config is None:
                continue  # Its structured loading failure was asserted above.
            report = validate_task_config(config)
            schema_errors = [item for item in report.errors() if item.code.startswith("JV2-SCH")]
            if not schema_errors:
                continue
            if any(_is_legacy_utils_error(item) for item in schema_errors):
                legacy_utils_cards.append(str(card))
            leftover = [
                item
                for item in schema_errors
                if not _is_legacy_utils_error(item)
                and not _is_legacy_ptmcmc_bounds_error(item)
            ]
            if leftover:
                bad.append(
                    f"{card}:\n" + "\n".join(item.format_line() for item in leftover)
                )
        self.assertEqual(bad, [])
        self.assertGreater(legacy_utils_cards, [])

    def test_unmigrated_example_methods_are_runtime_errors_not_schema_errors(self) -> None:
        if not EXAMPLES_ROOT.is_dir():
            self.skipTest("requires sibling Jarvis-Examples checkout")
        cards = sorted(EXAMPLES_ROOT.glob("*/bin/*.yaml"))
        expected = {"DNN", "DREAM", "ESS", "HMC", "MALA", "NUTS", "RLTPMCMC", "RobustAM", "SliceMCMC"}
        seen: set[str] = set()
        for card in cards:
            config = self._load_external_card(card)
            if config is None:
                continue
            report = validate_task_config(config)
            for item in report.errors():
                if item.code == "JV2-MTH-003":
                    seen.add(str(config["Sampling"]["Method"]))
        self.assertEqual(seen, expected)

    @unittest.skipUnless(EXAMPLES_ROOT.is_dir(), "requires sibling Jarvis-Examples checkout")
    def test_example_cards_have_no_multimode_dependency_regressions(self) -> None:
        cards = sorted(EXAMPLES_ROOT.glob("*/bin/*.yaml"))
        self.assertGreaterEqual(len(cards), 65)
        bad: list[str] = []
        for card in cards:
            config = self._load_external_card(card)
            if config is None:
                continue
            report = validate_task_config(config)
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
