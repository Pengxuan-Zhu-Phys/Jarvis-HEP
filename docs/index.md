<nav>
  
**Home** · _Previous_ · **[Next](guide_to_yaml.md)**

</nav>

---
# Jarvis-HEP Documentation

Jarvis-HEP is configured entirely through **YAML configuration files**.  
This documentation explains how to use Jarvis-HEP by *declaring workflows, sampling strategies, and external program calls in YAML*.

> **The YAML file is the user interface of Jarvis-HEP.**

This documentation is intentionally **configuration-driven** rather than API-driven.

---

## Tutorial Overview

The **Tutorial** explains how to construct a complete Jarvis-HEP YAML file step by step.

The goal of the tutorial is:
- To teach users how to *express their intent* in YAML
- To understand how different YAML blocks activate different functionality
- To enable users to design scans without reading internal source code

The tutorial does **not**:
- Teach physics models or domain-specific assumptions
- Document internal Python classes or functions
- Replace reference tables or example configurations

📘 **Start here:**  
→ [Guide to YAML](guide_to_yaml.md)

---

## Documentation Overview

The documentation is organised around **how users interact with Jarvis-HEP**, rather than internal implementation details.

**[Guide to YAML](guide_to_yaml.md)**  
A step-by-step tutorial that walks through the construction of a complete YAML configuration file.  
This document introduces the overall workflow and explains how different YAML blocks work together.

**[YAML Reference](yaml_reference.md)**  
A compact reference of supported YAML blocks and keys.  
This document is intended as a lookup table when editing or extending existing configurations.

**[Examples](examples.md)**  
A collection of curated example YAML files illustrating common usage patterns.  
These examples are designed to be copied and adapted to new use cases.

**[FAQ](faq.md)**  
A list of common pitfalls, silent failure modes, and debugging hints.  
This document focuses on practical issues encountered when running scans.

---

## Examples and External Package Interfaces

This section provides a collection of **standard interface patterns** for working with external programs.

The focus is on:
- Declaring external tools in YAML
- Defining input and output mappings
- Composing multiple tools within a single workflow

These examples are intended to be generic and reusable, illustrating common interaction patterns rather than project-specific integrations.

📂 **Examples index:**  
→ [Examples](examples.md)

📄 **Repository example YAML:**  
→ [Example_Grid.yaml](../bin/EggBox/Example_Grid.yaml)

---

## How to Use This Documentation

Recommended reading order:

1. **Guide to YAML** — understand the full workflow once
2. **Examples** — start from a working configuration
3. **YAML Reference** — look up block/key meanings when needed
4. **FAQ** — diagnose common mistakes and silent failures

Jarvis-HEP is designed so that **a well-structured YAML file fully defines the behaviour of a scan**.

---
<nav>
  
**Home** · _Previous_ · **[Next](guide_to_yaml.md)**

</nav>
