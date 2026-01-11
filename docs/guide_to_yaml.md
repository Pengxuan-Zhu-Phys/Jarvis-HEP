<nav>
  
**[Previous](intro2yaml_format.md)** · **[Home](index.md)** · **[Next]()**

</nav>

# Jarvis-HEP YAML Configurations

**Author:** Pengxuan Zhu

**Date:** 11 Feb 2026

Jarvis-HEP is configured entirely through **YAML configuration files**.  
This documentation explains how to use Jarvis-HEP by *declaring workflows, sampling strategies, and external program calls in YAML*.

> **The YAML file is the user interface of Jarvis-HEP.**

This documentation is intentionally **configuration-driven** rather than API-driven.

---

## Scope of This Document

This document introduces the **conceptual structure and design principles**
behind Jarvis-HEP YAML configurations.

It explains:
- How a complete computational workflow is expressed in YAML
- Why sampling and calculation are treated as separate concerns
- How external programs are integrated as black-box calculators
- Why the YAML file serves as the single source of truth for a workflow

This document does **not**:
- Provide a key-by-key YAML reference
- Describe individual samplers or calculators in full detail
- Replace module-specific tutorials or examples

Detailed documentation for each top-level section is provided in
dedicated documents referenced throughout this guide.

---

## Design Philosophy and Workflow

This section explains the **conceptual design** of Jarvis-HEP and the rationale behind its configuration model.
It is intended to help users understand *why* the YAML structure is designed in this way, before learning *how* to use it.

---

## How a Jarvis-HEP Workflow Is Described

A Jarvis-HEP YAML file describes a complete computational workflow.
Rather than focusing on syntax, this section explains **what each major part of the configuration is responsible for**.

The YAML file is organised into a small number of **top-level sections**, each corresponding to a dedicated document.
These sections together define *what is executed*, *how it is executed*, and *how results are managed*.

**[`Scan`](scan.md)**  
Defines the identity of a scan, bookkeeping information, and where results are stored.
This section records *what this run is* and *how it should be archived*.

**[`Sampling`](sampling.md)**  
Declares what is sampled and which strategy is used to explore the parameter space.
This section controls *which points are proposed* without performing any calculations.

**[`EnvReqs`](envreqs.md)**  
Specifies environment assumptions, platform constraints, and dependency checks.
This section documents *what the execution environment must satisfy* before a workflow is run.

**[`Calculators`](calculators.md)**  
Describes external programs, how they are executed, and how inputs and outputs are mapped.
This section defines *how physics or numerical evaluations are performed* while treating programs as black boxes.

**[`Utils`](utils.md)**  
Defines auxiliary utilities and shared resources used across the workflow.
This section provides *supporting functionality* that is reused by calculators and other components.

Only sections that are explicitly defined are active.
If a section is omitted, the corresponding functionality is not used.

---

### A Concrete Example

For a complete, working example of how these sections are composed, see:

- [Example_Grid.yaml](../bin/EggBox/Example_Grid.yaml)

This example illustrates how a full workflow can be described using only YAML,
without modifying external programs or writing auxiliary scripts.

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
