<nav>
  
**[Previous](intro2yaml_format.md)** · **[Home](index.md)** · **[Next]()**

</nav>

# Jarvis-HEP YAML Configurations

**Author:** Pengxuan Zhu

**Date:** 11 Feb 2026

This document describes **how to use Jarvis-HEP by writing YAML configuration files**.

Jarvis-HEP is configured entirely through YAML.  
Users are not expected to interact with internal Python APIs.

> **If you understand the YAML structure and conventions, you understand how to use Jarvis-HEP.**

---

## Scope of This Document

This document explains:

- The **top-level YAML blocks** recognised by Jarvis-HEP
- The **responsibility** of each block
- The **conventions** governing how blocks interact
- How user intent is expressed *declaratively* through YAML

This document does **not** cover:
- Physics models or domain-specific assumptions
- Internal class or function implementations
- Auto-generated API references

---

## Design Philosophy and Workflow

This section explains the **conceptual design** of Jarvis-HEP and the rationale behind its configuration model.
It is intended to help users understand *why* the YAML structure is designed in this way, before learning *how* to use it.

---

### Separation of Sampling and Calculation

Jarvis-HEP enforces a strict separation between:

- **Sampling**: deciding *which points* in parameter space to evaluate
- **Calculation**: deciding *how a given point is evaluated*

Sampling algorithms are responsible only for proposing points.
They do not perform physics calculations, manage external programs, or interpret results.

Conversely, physics calculations are treated as **black-box evaluations**:
a set of inputs is provided, and a set of outputs is produced.

This separation avoids tight coupling between numerical exploration strategies and domain-specific computation logic.

---

### Sampler and Factory

In Jarvis-HEP, the **Sampler** and the **Factory** have clearly distinct roles:

- The **Sampler** generates candidate points according to a defined strategy.
- The **Factory** is responsible for constructing executable workflows from the YAML configuration.

The Factory:
- Instantiates calculators and utilities
- Connects inputs and outputs
- Manages execution order and data flow

The Sampler never interacts directly with external programs.
Instead, it delegates all execution-related responsibilities to the Factory.

This design allows sampling strategies to be changed without modifying calculation logic, and vice versa.

---

### Calculators as Black Boxes

Jarvis-HEP treats physics programs as **opaque calculators or libraries**.

From the framework’s perspective, a calculator is defined by:
- Its required inputs
- Its produced outputs
- The commands needed to execute it

Internal implementation details are intentionally ignored.

Users are encouraged **not** to repeatedly modify calculation scripts during a project.
Instead, external programs should be treated as stable black boxes once validated.

This approach:
- Reduces debugging risk
- Improves reproducibility
- Prevents accidental divergence between runs

---

### YAML as the Single Source of Truth

In Jarvis-HEP, the YAML configuration file is the **authoritative description** of a workflow.

A complete analysis should be reproducible using only:
- The YAML configuration file
- The original external programs
- The generated data

No auxiliary scripts or ad-hoc glue code are required.

This design ensures that:
- Workflow logic is explicit
- Configuration changes are version-controlled
- Long-term projects remain maintainable

---

### Long-Term Reproducibility

The primary goal of this design is **long-term reproducibility**.

After completing a project, users should only need to archive:
- Data outputs
- The YAML configuration
- The external program packages

To reproduce results in the future, the same YAML file can be executed again using the original programs.

This is significantly more robust than maintaining collections of fragmented scripts whose intent and dependencies may be unclear over time.

---

### Where to Find Concrete Tutorials

This document focuses on **concepts and structure**.

Detailed tutorials for each module and YAML block are provided in dedicated documents:
- The overall workflow tutorial is covered in *Guide to YAML*
- Concrete configuration examples are provided in *Examples*
- Complete key-by-key definitions are available in *YAML Reference*

Users are encouraged to read this section once to understand the design intent, then proceed to the practical tutorials.

---

## Top-Level YAML Structure

A Jarvis-HEP configuration file is composed of the following **top-level blocks**:

| YAML Block | Purpose |
|-----------|--------|
| `Scan` | Scan identity, bookkeeping, and global run metadata |
| `Sampling` | Definition of variables and sampling strategy |
| `EnvReqs` | Environment and dependency requirements |
| `Calculators` | External programs and execution workflow |
| `Utils` | Auxiliary utilities (interpolations, helpers, etc.) |

Blocks are **activated only when explicitly defined**.  
Undefined blocks are ignored.

---

## Block Overview

### `Scan`

Defines global metadata and output location for a scan.

| Key | Type | Description |
|----|------|-------------|
| `name` | string | Human-readable scan identifier |
| `save_dir` | string | Output directory (supports path macros) |

---

### `Sampling`

Defines **what is sampled** and **how sampling is performed**.

#### `Sampling.Method`

| Key | Type | Description |
|----|------|-------------|
| `Method` | string | Sampling method identifier (e.g. grid, random) |

#### `Sampling.Variables`

A list of variable definitions.

Each variable entry follows the pattern:

| Key | Type | Description |
|----|------|-------------|
| `name` | string | Variable identifier |
| `description` | string | Optional human-readable description |
| `distribution` | mapping | Sampling distribution definition |

#### `distribution`

| Key | Type | Description |
|----|------|-------------|
| `type` | string | Distribution type |
| `parameters` | mapping | Distribution-specific parameters |

---

### `EnvReqs`

Defines runtime environment requirements.

| Key | Type | Description |
|----|------|-------------|
| `OS` | list | Supported operating systems |
| `Check_default_dependences` | mapping | Dependency checking behaviour |

---

### `Calculators`

Defines **external programs**, their setup, execution, and data exchange.

| Key | Type | Description |
|----|------|-------------|
| `make_paraller` | integer | Parallel execution level |
| `path` | string | Base working directory |
| `Modules` | list | External module definitions |

Each module defines:
- Installation steps
- Initialisation steps
- Execution commands
- Input/output data mapping

---

### `Utils`

Defines auxiliary utilities available to the workflow.

Typical use cases include:
- Interpolation tables
- Lookup functions
- Reusable helper definitions

---

## General Conventions

Jarvis-HEP follows a small number of **strict conventions**:

| Convention | Meaning |
|-----------|--------|
| Declarative configuration | YAML declares *what*, not *how* |
| Explicit activation | Only defined blocks are active |
| Separation of concerns | Sampling, execution, and utilities are isolated |
| No implicit defaults | Behaviour is not assumed without YAML |

Violating these conventions may result in undefined behaviour.

---

## How to Learn Jarvis-HEP

Recommended learning order:

1. Understand the top-level YAML blocks (this document)
2. Study example YAML files
3. Modify examples incrementally to match your use case
4. Introduce additional blocks only when required

Jarvis-HEP is designed so that **correct YAML structure leads to correct execution behaviour**.

---
