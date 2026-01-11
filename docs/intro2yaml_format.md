# Introduction to YAML Format

**Author:** Pengxuan Zhu

This document provides a **short, practical introduction to YAML** for users who are not familiar with the format.
It is intended to help you read and write YAML files correctly before using them to configure Jarvis-HEP.

YAML is a human-readable data-serialization format commonly used for configuration files.

---

## What Is YAML?

YAML stands for **“YAML Ain’t Markup Language”**.
It is designed to be:
- Easy to read and write
- Indentation-based
- Friendly to version control systems

YAML is often used as an alternative to JSON or XML for configuration purposes.

---

## Core Syntax Rules

### 1. Indentation Matters

YAML uses **indentation to represent structure**.
Spaces define hierarchy.

- Use **spaces**, not tabs
- Be consistent (commonly 2 or 4 spaces)

**Correct**
```yaml
parent:
  child:
    key: value
```

**Incorrect**
```yaml
parent:
 child:
    key: value
```

---

### 2. Key–Value Pairs

A basic YAML element is a key–value pair:

```yaml
key: value
```

Keys are strings.
Values can be strings, numbers, booleans, lists, or mappings.

---

### 3. Lists

Lists are written using hyphens (`-`).

```yaml
items:
  - apple
  - orange
  - banana
```

Lists of mappings are also common:

```yaml
users:
  - name: Alice
    role: admin
  - name: Bob
    role: user
```

---

### 4. Strings and Quoting

Most strings do **not** require quotes.

```yaml
name: example
path: /tmp/data
```

Use quotes when:
- The string contains special characters
- The value could be misinterpreted

```yaml
expression: "a > b"
version: "1.0"
```

---

### 5. Numbers and Booleans

YAML automatically interprets numbers and booleans.

```yaml
count: 10
ratio: 0.5
enabled: true
disabled: false
```

Avoid quoting numbers unless they must be treated as strings.

---

### 6. Comments

Comments start with `#`.

```yaml
# This is a comment
key: value  # Inline comment
```

Comments are ignored by the parser.

---

## Common Mistakes

### ❌ Using Tabs

```yaml
key:
\tchild: value
```

Tabs are not allowed. Always use spaces.

---

### ❌ Inconsistent Indentation

```yaml
key:
  child1: value
   child2: value
```

Even a single extra space changes the structure.

---

### ❌ Missing Spaces After Colon

```yaml
key:value
```

Correct form:

```yaml
key: value
```

---

### ❌ Ambiguous Values

```yaml
yes: no
on: off
```

These may be interpreted as booleans in YAML 1.1.
Use explicit values or quotes if unsure.

---

## Useful Features (Optional)

### Anchors and Aliases

YAML supports reusable blocks using anchors (`&`) and aliases (`*`).

```yaml
defaults: &defaults
  timeout: 10
  retries: 3

task1:
  <<: *defaults
  name: task_one
```

This reduces duplication in large configuration files.

---

## How to Validate YAML

Before using a YAML file, it is strongly recommended to validate it.

Options include:
- Online validators (e.g. yamllint)
- Command-line tools
- Parsing the file in Python using `yaml.safe_load`

---

## Recommended References

- YAML Official Specification  
  https://yaml.org/spec/

- YAML Tutorial (Learn X in Y Minutes)  
  https://learnxinyminutes.com/docs/yaml/

- PyYAML Documentation  
  https://pyyaml.org/wiki/PyYAMLDocumentation

These resources provide deeper explanations and advanced usage examples.

---

## Summary

- YAML is indentation-based and human-readable
- Consistent spacing is critical
- Start simple and validate often
- When in doubt, prefer clarity over brevity

Once you are comfortable with YAML syntax, you can proceed to the Jarvis-HEP tutorial.

→ **Next:** [Guide to YAML](guide_to_yaml.md)
