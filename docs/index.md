# Jarvis-HEP Documentation

This documentation is designed to help users **learn how to use Jarvis-HEP by writing YAML configuration files**, without needing to read the source code.

The guiding principle is simple:

> **If you understand the YAML, you understand how to use Jarvis-HEP.**

Jarvis-HEP does not expose a large user-facing Python API.  
Instead, functionality is composed by **declaring components and conventions in YAML**, which are then orchestrated by the framework.

---

## What This Documentation Covers

This documentation focuses on:

- What sections a Jarvis-HEP YAML file can contain
- What each section is responsible for
- What conventions Jarvis-HEP assumes
- How different choices in YAML activate different functionality

It intentionally avoids:
- Auto-generated API references
- Internal class-level details
- Exhaustive parameter lists without context

---

## Core YAML Structure (Overview)

A typical Jarvis-HEP configuration file is composed of the following logical blocks:

| YAML Section | Purpose |
|-------------|---------|
| `run` | Global run control (modes, bookkeeping, output behaviour) |
| `parameters` | Definition of parameters of interest |
| `nuisance_parameters` | Definition and handling of nuisance parameters |
| `sampler` | Choice and configuration of the sampling algorithm |
| `likelihood` | Likelihood construction and evaluation logic |
| `pass_condition` | Hard acceptance / rejection criteria |
| `output` | Data storage and output formats |

Not all sections are mandatory.  
Jarvis-HEP activates functionality **only when a section is explicitly defined**.

---

## Minimal Working Example

The following is a **minimal but valid** Jarvis-HEP YAML file:

```yaml
run:
  mode: scan

parameters:
  M1:
    range: [100, 1000]
  M2:
    range: [-2000, 2000]

sampler:
  type: random
  n_points: 1000

likelihood:
  type: external
```

This configuration instructs Jarvis-HEP to:
- Scan over parameters `M1` and `M2`
- Use a random sampler
- Evaluate a user-defined external likelihood

---

## Design Conventions You Must Know

Jarvis-HEP relies on a small number of **strict conventions**:

| Convention | Meaning |
|-----------|--------|
| Explicit is better than implicit | Undefined YAML blocks are ignored |
| Parameters are continuous by default | Discrete variables must be handled explicitly |
| Samplers operate on parameters of interest | Nuisance parameters are treated separately |
| Likelihoods do not control sampling | Sampling strategy and likelihood evaluation are decoupled |

Understanding these conventions is essential for correct usage.

---

## Learning Path (Recommended Order)

To learn Jarvis-HEP efficiently, read this documentation in the following order:

1. YAML structure and conventions (this page)
2. Parameter and nuisance-parameter definitions
3. Sampler configuration patterns
4. Likelihood and profiling behaviour
5. Output and data interpretation

Each topic is documented using **tables, short examples, and explicit rules** rather than prose.

---

## Philosophy of the Documentation

This documentation is intentionally:
- Compact
- Declarative
- Convention-driven

Jarvis-HEP is designed for users who prefer **understanding behaviour through configuration**, rather than through large APIs.

If you can write a correct YAML file, you can use Jarvis-HEP.

---
