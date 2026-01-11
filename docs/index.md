<nav>
  
**[Home](index.md)** · **[Next](intro2yaml_format.md)**

</nav>

---
# Jarvis-HEP Documentation
**Author:** Pengxuan Zhu

**Date:** 11 Feb 2026

---

## Tutorial Overview

📘 **Start here:**  

**[Eggbox example](tutorial_eggbox.md)**  

The tutorial begins with a **concrete working example** based on the Eggbox function. Starting from this simple example, you are guided through how a numerical or physical problem is translated into a Jarvis-HEP scan.

Using the Eggbox example as a reference, the tutorial then introduces the
**YAML configuration** used by Jarvis-HEP, including:
- The overall structure of a YAML configuration file
- The role of different top-level YAML blocks
- Basic YAML syntax and common patterns used in Jarvis-HEP

The goal of the tutorial is:
- To teach users how to *express their intent* in YAML
- To explain how different YAML blocks activate different functionality
- To help users understand the structure of a complete Jarvis-HEP configuration file

The tutorial does **not**:
- Teach physics models or domain-specific assumptions
- Document internal Python classes or implementation details
- Replace reference tables or example configurations

→ [Guide to YAML](guide_to_yaml.md)

---

## Documentation Overview

**[Introduction to YAML Format](intro2yaml_format.md)**  
A concise introduction to YAML syntax and rules, including correct and incorrect examples.  
This document is intended for users who are not yet familiar with YAML.

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

1. **Introduction to YAML Format** — learn basic YAML syntax and rules
2. **Guide to YAML** — understand the full Jarvis-HEP workflow once
3. **Examples** — start from a working configuration
4. **YAML Reference** — look up block/key meanings when needed
5. **FAQ** — diagnose common mistakes and silent failures

Jarvis-HEP is designed so that **a well-structured YAML file fully defines the behaviour of a scan**.

---
<nav>
  
**[Home](index.md)** · **[Next](intro2yaml_format.md)**

</nav>
