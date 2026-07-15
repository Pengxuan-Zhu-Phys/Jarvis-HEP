# Jarvis-HEP V2 Project

Created with:

```bash
Jarvis2 project create <name>
```

## Layout

| Path | Role |
|------|------|
| `bin/` | Task YAML cards |
| `data/` | Small input tables / fixtures |
| `deps/` | Default environment policy (`environment_default.yaml`) |
| `jarvis.project.yaml` | Project descriptor (`&J` root) |
| `.jarvis-project.json` | Machine-readable layout marker |

Runtime directories (`outputs/`, `logs/`, `images/`, `checkpoints/`) appear on first run.

## Quick start

```bash
cd <project>
Jarvis2 run bin/quickstart_bridson_operas.yaml
# or
Jarvis2 bin/quickstart_bridson_operas.yaml
```

CSV operas smoke:

```bash
Jarvis2 run bin/quickstart_csv_operas.yaml
```

## Scheduling defaults

`deps/environment_default.yaml` supplies `EnvReqs.V2` (workers, SAMPLE buckets, archiver).
Override in the task YAML or edit the defaults file. There is **no** top-level `Runtime`
block on V2 — use `EnvReqs.V2` instead.

## Package for sharing

```bash
Jarvis2 project pack . --share
Jarvis2 project pack . --repro
Jarvis2 project pack . --full
Jarvis2 project pack . --man    # write a pack manifest only
```
