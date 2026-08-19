# Jarvis-HEP V2 Project

Created with:

```bash
Jarvis project create <name>
```

## Layout

| Path | Role |
|------|------|
| `bin/` | Task YAML cards |
| `bin/sampling/` | Nested sampling `Sampling:` templates (Dynesty / MultiNest) |
| `data/` | Small input tables / fixtures |
| `deps/` | Default environment policy (`environment_default.yaml`) |
| `jarvis.project.yaml` | Project descriptor (`&J` root) |
| `.jarvis-project.json` | Machine-readable layout marker |

Runtime directories (`outputs/`, `logs/`, `images/`, `checkpoints/`) appear on first run.

## Quick start

```bash
cd <project>
Jarvis run bin/quickstart_bridson_operas.yaml
# or
Jarvis bin/quickstart_bridson_operas.yaml
```

CSV operas smoke:

```bash
Jarvis run bin/quickstart_csv_operas.yaml
```

### Look up YAML while editing

`Jarvis man` documents **task-card YAML** (keys, paths, copy-paste examples that should validate).
List-valued YAML fields are queried by field name alone; the Keys table reports `list`.
For Portal adapter runtime or Operas operator catalogs, use their CLIs:

```bash
Jarvis man yaml.Calculators.Modules.execution
Jarvis man calculator.execution.output --type JSON
Jarvis man operas
Jarvis man example.calculator
Jarvis man example.random --json
Jarvis man example.random-operas --json
Jarvis man example.random-calculator --json
Jarvis man yaml.EnvReqs.V2
Jarvis portal man JSON          # runtime adapter behaviour
Jarvis operas info helper.eggbox2d
Jarvis validate bin/quickstart_bridson_operas.yaml
```

Random migration examples:

```bash
Jarvis validate bin/quickstart_random.yaml
Jarvis validate bin/quickstart_random_operas.yaml
Jarvis check bin/quickstart_random_calculator.yaml
```

### Nested sampling templates

Under `bin/sampling/` (copy the `Sampling:` block into your card):

| Template | Notes |
|----------|--------|
| `Sampling_Dynesty_Simple.yaml` | Dynesty = DynamicNestedSampler (everyday) |
| `Sampling_Dynesty_Full.yaml` | Dynesty full dynamic `run_nested` + constructor |
| `Sampling_MultiNest_Simple.yaml` | MultiNest = static NestedSampler (everyday) |
| `Sampling_MultiNest_Full.yaml` | MultiNest full static surface only |

See `bin/sampling/README.md`.

## Scheduling defaults

`deps/environment_default.yaml` supplies `EnvReqs.V2` (workers, SAMPLE buckets, archiver,
and `checkpoint_heartbeat_sec`).  The heartbeat defaults to 30 seconds; lower it for
expensive samples when you want a tighter resume replay bound.
Override in the task YAML or edit the defaults file. There is **no** top-level `Runtime`
block on V2 — use `EnvReqs.V2` instead.

## Package for sharing

```bash
Jarvis project pack . --share
Jarvis project pack . --repro
Jarvis project pack . --full
Jarvis project pack . --man    # write a pack manifest only
```

`--share` includes the project markers, README, `bin/`, `data/`, `deps/`,
`outputs/`, and `images/`. It deliberately excludes runtime `calculators/`,
`logs/`, and `checkpoints/`. To omit a project-specific dependency subtree,
declare it explicitly in `jarvis.project.yaml`:

```yaml
pack:
  exclude:
    - deps/private-models/
```

### Restricted (encrypted) release — CLI only

Do **not** run `openssl` by hand. Use:

```bash
# Pack + encrypt → *.tar.gz.jenc
Jarvis project pack . --repro --encrypt --key 'YOUR_KEY'

# Or encrypt an existing tarball
Jarvis project encrypt path/to/archive.tar.gz --key 'YOUR_KEY'
```

Collaborators fetch with:

```bash
Jarvis project list
Jarvis project fetch YourProjectName --key 'YOUR_KEY'
# or: export JARVIS_PROJECT_FETCH_KEY='YOUR_KEY'
```

Official catalog (public list + restricted entries) lives in
**Jarvis-Examples** `catalog/official_project_library.json` — not a PyPI package.
See `Jarvis-Books/Jarvis-HEP V2/components/project_tools.md` and `INSTALL.md`.
