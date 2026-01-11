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

## Documentation Structure

The documentation is organised around **usage patterns**, not internal modules.

```
docs/
├── index.md              # This page: overview and navigation
├── guide_to_yaml.md      # Step-by-step tutorial for writing YAML
├── yaml_reference.md     # Reference tables for YAML blocks and keys
├── examples.md           # Curated example configurations
├── faq.md                # Common issues and debugging hints
```

Each document has a distinct role.  
Users are encouraged to **read the tutorial once**, then rely on references and examples.

---

## Examples and External Package Interfaces

In addition to the tutorial, Jarvis-HEP provides **standardised example patterns** for interfacing with external programs.

These examples summarise:
- How external codes are declared in YAML
- How inputs are generated and passed
- How outputs are collected and mapped back into the workflow
- How multiple external tools can be composed

This approach is similar in spirit to *interface examples* used in projects such as BSMArt, but adapted to Jarvis-HEP’s declarative YAML model.

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
