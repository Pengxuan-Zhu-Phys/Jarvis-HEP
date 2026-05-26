# Jarvis-HEP Standalone Project

This directory is a self-contained standalone project scaffold created by:

```bash
Jarvis project create <PROJECT_NAME>
```

It is intended to serve as the project-local workspace for runnable cards, input data, bundled defaults, and runtime-generated outputs.

## Quick Start

Run one of the built-in template examples:

```bash
Jarvis bin/quickstart_mcmc_operas.yaml
```

Run the CSV replay example:

```bash
Jarvis bin/quickstart_csv_operas.yaml
```

## Layout

- `.jarvis-project.json` / `jarvis.project.yaml`: markers that identify the standalone project root
- `bin/`: runnable YAML entry cards
- `data/`: project input datasets
- `deps/`: project-local bundled defaults, including `environment_default.yaml`

Runtime directories are created automatically on first run:

- `outputs/<scan>/DATABASE`: HDF5, CSV, schema, and run metadata
- `outputs/<scan>/SAMPLE`: per-sample artifacts and sample-local logs
  (`Jarvis ... --check-modules` uses `outputs/<scan>/SAMPLE/tests`)
  (`DATABASE` remains `outputs/<scan>/DATABASE`, not under `SAMPLE/`)
- `logs/<scan>/`: Jarvis / sampler / factory logs
- `images/<scan>/`: plots, generated plotting YAML, semantic `flowchart.json`,
  and JarvisPLOT-rendered `flowchart.png` when JarvisPLOT is installed

Add optional project directories such as `calculators/`, `configs/`, `scripts/`, `assets/`, and `docs/` only when your workflow actually needs them.

## Path Rules

- `&J/...` resolves against the standalone project root
- Prefer `deps/` for project-local bundled defaults such as `&J/deps/environment_default.yaml`

## LogLikelihood Notes

`Sampling.LogLikelihood` is a list of named expressions. Jarvis-HEP evaluates
every item and writes each value into the sample outputs.

If no item is named `LogL`, the total `LogL` is the sum of all listed items. If
an item is named `LogL`, that expression defines the total likelihood instead,
and it can reference earlier named items:

```yaml
Sampling:
  LogLikelihood:
    - {name: "LogL_Z", expression: "LogGauss(z, 100, 10)"}
    - {name: "LogL_H", expression: "LogGauss(h, 125, 2)"}
    - {name: "LogL", expression: "LogL_Z + 0.5 * LogL_H"}
```
