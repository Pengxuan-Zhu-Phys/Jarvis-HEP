# Jarvis-HEP

Jarvis-HEP (Just a Really Viable Interface to Suites for High Energy Physics) is an open-source, modular Python framework for
likelihood-based parameter scanning and global fits in
high-energy physics phenomenology.

The project focuses on **robust sampling strategies**, **nuisance-parameter handling**, and **scalable asynchronous workflows**, with an emphasis on *profile likelihood*–oriented studies commonly encountered in modern collider phenomenology.

---

## Motivation

Modern HEP analyses often suffer from:
- Extremely sparse viable regions in high-dimensional parameter spaces
- Expensive external likelihood evaluations (collider simulations, Higgs constraints, flavour physics, etc.)
- Non-trivial treatment of nuisance parameters
- Poor reproducibility and opaque analysis pipelines

Jarvis-HEP is designed to address these issues by providing:
- A **unified orchestration layer** for sampling, likelihoods, and external tools
- Explicit separation between **parameters of interest** and **nuisance parameters**
- Flexible, engineering-oriented solutions such as **profile likelihood–based inference**
- Transparent data management and diagnostic logging

---

## Key Features

### Sampling & Exploration
- Multiple sampling strategies (random, grid-like, adaptive, nested-style, and custom algorithms)
- Designed for **highly constrained and fine-tuned** parameter spaces
- Iterator-style, point-level sampling (not generation-locked)

### Likelihood & Nuisance Parameters
- Native support for **profile likelihood construction**
- Dedicated nuisance-parameter samplers decoupled from main exploration
- Fast evaluation paths for nuisance optimization

### Architecture
- Pure Python implementation
- Asynchronous execution of expensive external programs
- Modular factory-based design (samplers, likelihoods, IO, monitors)
- YAML-driven configuration for reproducibility

### Data Handling & Diagnostics
- HDF5-based structured output
- Schema-driven CSV flattening for structured observables (`samples.schema.json`)
- Explicit resource and file-handle monitoring
- Detailed logging via `loguru`
- Designed to survive partial failures (no silent data loss)

### Visualization
- Dynamic sampling visualizations (under active development)
- Designed to interface with the standalone plotting package **JarvisPLOT**

Example visualization of an adaptive sampling procedure:  
👉 https://github.com/Pengxuan-Zhu-Phys/Jarvis-HEP/blob/master/simu/sample_dynamic_viz.gif

---
📘 **Full documentation and tutorials are hosted on a dedicated documentation site:**

👉 https://pengxuan-zhu-phys.github.io/Jarvis-Docs/
---

## Installation

Jarvis-HEP is a pure Python project.

### Install from PyPI

```bash
python3 -m pip install Jarvis-HEP
```

Dependency installation is handled by `pip` from package metadata at install time.
Jarvis runtime no longer provides a dependency-installer CLI mode.

After installation, run directly:

```bash
Jarvis --help
Jarvis --version
```

### Install from source (developer mode)

```bash
python3 -m pip install -e .
```

### Release helper (maintainer)

```bash
./jhrel 1.0.1 --dry
./jhrel 1.0.1 --testpypi
./jhrel 1.0.1 --reinstall
```

Maintainer release workflow:
- [docs/RELEASE_PLAYBOOK.md](docs/RELEASE_PLAYBOOK.md)
- Version/banner metadata source (inside package, includes resources and citation references):
  - `jarvishep/card/document_links.json`

### Create a standalone project workspace

Use `--mkproject` to create a fresh Jarvis project directory:

```bash
Jarvis --mkproject PROJECT_NAME
```

This creates a standalone project workspace with:

- project markers: `.jarvis-project.json`, `jarvis.project.yaml`
- runnable YAML: `bin/`
- project data: `data/`
- run outputs: `outputs/`
- toolchain folders: `calculators/`, `deps/`
- project artifacts: `images/`, `logs/`, `checkpoints/`

Quick start after scaffold:

```bash
cd PROJECT_NAME
Jarvis bin/quickstart_mcmc_operas.yaml
```

### Pack a standalone project

After a run finishes, package the project with:

```bash
Jarvis --packproject [path-to-project-root] --profile share|repro|full
```

Rules:

- `path-to-project-root` is optional; default is current directory.
- if no Jarvis project marker is found in target path, CLI reports:
  - "Jarvis project not found in target path."
- `--profile` default is `repro`.

Profile intent:

- `share`: lightweight sharing package (results-focused)
- `repro` (default): reproducible rerun package
- `full`: full archival package

Recommended default output path for generated archive:

- parent directory of project root (`dirname(&J)`)

## Running

Jarvis-HEP can now be launched without changing into the source root:

```bash
Jarvis /path/to/project/bin/task.yaml
```

Or via module mode:

```bash
python -m jarvishep /path/to/project/bin/task.yaml
```

### Path Markers

- `&J/...` resolves to the runtime task root (project marker first, then YAML-location inference)
- `&SRC/...` resolves to the installed `jarvishep` package resources (internal cards/schema/logo)

### Observable Schema And CSV Flattening

Jarvis-HEP stores scan data in HDF5 and exports CSV through a user-editable schema file:

- raw records: `.../DATABASE/samples.N.hdf5`
- schema rules: `.../DATABASE/samples.schema.json`
- CSV export: `.../DATABASE/samples.N.csv`

`samples.schema.json` records each column's type metadata and flatten rules.  
You can edit flatten rules and regenerate CSV without rerunning sampling:

```bash
Jarvis /path/to/task.yaml --convert
```

Supported column flatten modes:

- `scalar`: scalar-first export (structured values fall back to JSON cell text)
- `json`: write structured value as one JSON string cell
- `split`: expand structured value to multiple CSV columns (supports `name_map` rename mapping)
- `drop`: omit column from CSV export

Detailed schema fields and examples:

- [docs/OBSERVABLE_SCHEMA.md](docs/OBSERVABLE_SCHEMA.md)

---

## Contributing

Contributions are welcome in all forms, including feature proposals, documentation improvements, bug reports, and bug fixes.  
Please refer to `CONTRIBUTING.md` for guidelines on how to get started.

---

## License

Jarvis-HEP is released under the **MIT License**.  
See the [`LICENSE`](LICENSE) file for details.

---

## Acknowledgements

The author thanks  
**Yang Zhang** and **Liangliang Shang**  
for helpful discussions during the development of this project.

---

## References

- **Exploring supersymmetry with machine learning**  
  Jie Ren, Lei Wu, Jin-Ming Yang, Jun Zhao  
  *Nuclear Physics B*, 2019  
  arXiv:1708.06615
