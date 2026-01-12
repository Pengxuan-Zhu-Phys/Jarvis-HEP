# Jarvis-HEP

**Jarvis-HEP** (Just a Really Viable Interface to Suites for High Energy Physics) is a modular, Python-based framework for **high-energy physics phenomenology**, designed to streamline **parameter-space exploration, likelihood evaluation, and post-processing** in complex BSM and collider analyses.

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

Jarvis-HEP is a **pure Python** project. No need to compile, but need to install all python Dependencies.

### Dependencies

Install required packages via `pip`:

```bash
python3 -m pip install \
  numpy scipy pandas pyyaml h5py matplotlib \
  dynesty pyhf networkx shapely sympy \
  xslha pyslha xmltodict \
  sqlalchemy aiofiles dill emoji prettytable loguru
```

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
