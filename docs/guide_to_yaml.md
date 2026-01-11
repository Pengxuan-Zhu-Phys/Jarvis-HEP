# Jarvis-HEP YAML Documentation

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
