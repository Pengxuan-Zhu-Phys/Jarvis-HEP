# Jarvis-HEP Standalone Project

This folder was created by:

```bash
Jarvis --mkproject <PROJECT_NAME>
```

## Quick Start

Run the built-in toy workflow:

```bash
Jarvis bin/quickstart_mcmc_operas.yaml
```

Run CSV replay example:

```bash
Jarvis bin/quickstart_csv_operas.yaml
```

## Layout

- `.jarvis-project.json` / `jarvis.project.yaml`: project-root markers and metadata
- `bin/`: runnable YAML entry cards
- `data/`: project input datasets
- `outputs/`: default runtime outputs (`Scan.save_dir: "&J/outputs"`)
- `calculators/`: calculator runtime/shadow workspaces
- `deps/`: external program/library dependencies and project-local env baseline (`environment_default.yaml`)
- `images/`: project-local plots and exported figures
- `logs/`: project-level runtime logs
- `checkpoints/`: reserved for future resume/restart support

Path rules:
- `&J/...` resolves against project root.
- `&SRC/...` resolves against installed `jarvishep` package resources.
