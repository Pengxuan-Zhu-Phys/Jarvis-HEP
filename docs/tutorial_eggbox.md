<nav>

**[Previous](guide_to_yaml.md)** · **[Home](index.md)** · **[Next](examples.md)**

</nav>

# Tutorial: Declarative Scanning with a 2D Eggbox Example

This tutorial uses `Example_Bridson.yaml` as a **working example**.

The point is not the Eggbox model.
The point is **how to describe a scan in Jarvis-HEP**.

If you can follow this example, you can adapt Jarvis-HEP to other problems.

---

## The Problem Pattern

Most scan problems look like this:

- Pick some parameters
- Explore their values with a sampler
- Run an external program at each point
- Save the results

Jarvis-HEP is built for this pattern. 

### Eggbox Function as an Simplifed Package

The Eggbox function is a simple two-dimensional test function.
It takes two input parameters, `xx` and `yy`, and returns a single numerical value.

In this example, the function is defined as

``z(x, y) = (\sin x \cdot \cos y + 2)^5``.

The exact formula is not important for Jarvis-HEP.
What matters is that:
- One number is produced for each input point
- The output changes smoothly with the inputs
- The function has many local peaks and valleys

These features make it useful as a test case for scanning and workflow logic.

In Jarvis-HEP, this function is wrapped as a **minimal Python program**.
The program behaves like a very small external physics or analysis code.

The program does the following:
- Reads parameters from `input.json`
- Computes the Eggbox value
- Adds a dummy timing value to mimic runtime information
- Writes all results to `output.json`

The file `input.json` contains the sampled parameters provided by Jarvis-HEP.
The file `output.json` contains the results, for example:
- `z`: the computed Eggbox value
- `Time`: a random number used as a stand-in for execution time

From Jarvis-HEP’s point of view, this Python script is just an external executable.
Jarvis-HEP does not inspect the code or the mathematical expression.

Jarvis-HEP only needs to know:
- How to write `input.json`
- How to run the program
- How to read values from `output.json`

All of this is declared in the YAML configuration.
Jarvis-HEP then runs this program automatically for each sampled point.

---

## How Jarvis-HEP Models the Problem

Jarvis-HEP splits the task into parts:

- **Sampling**: which points to try
- **Calculation**: what to run at each point
- **Workflow**: how everything is connected and executed

You do not write scripts to glue these together.
You **declare** them in YAML.

See:
- [Guide to YAML](guide_to_yaml.md)

---

## Step 1: Name the Scan and Outputs

Open `Example_Bridson.yaml`.

The `Scan` section does basic setup:
- A scan name
- Output paths
- Bookkeeping info

This section does not change physics or algorithms.

See:
- [Scan](scan.md)

---

## Step 2: Choose a Sampling Method

The `Sampling` section answers two questions:
- Which sampler is used?
- Which variables are sampled?

In this example:
- The sampler is `Bridson`
- Two variables are defined
- Each variable has bounds

Jarvis-HEP uses this to generate points.
Nothing is executed yet.

See:
- [Sampling](sampling.md)

---

## Step 3: Set Sampler Controls

Some samplers need extra settings.

In `Example_Bridson.yaml` these include:
- A radius
- A limit on attempts

These settings affect *how points are proposed*.
They do not affect calculations.

---

## Step 4: Declare the Calculation

The `Calculators` section describes an external program.

Jarvis-HEP treats it as a **black box**.

You declare:
- How inputs are prepared
- How the program is run
- Where outputs are read from

Jarvis-HEP handles:
- Working directories
- Execution order
- Parallel runs

No glue scripts are needed.

See:
- [Calculators](calculators.md)

---

## Step 5: Connect Inputs and Outputs

Inside `Calculators`, mappings are defined.

These mappings say:
- Which sampled variables go to which inputs
- Which files or values are read as outputs

This is the only place where sampling meets calculation.

---

## Step 6: Check the Environment

The `EnvReqs` section lists requirements:
- Supported systems
- Required tools or libraries

Jarvis-HEP checks these before running.

See:
- [EnvReqs](envreqs.md)

---

## Step 7: Optional Utilities

The `Utils` section defines shared helpers.

Examples:
- Tables
- Interpolations
- Reusable data

These helpers do not control sampling or execution.

See:
- [Utils](utils.md)

---

## What Jarvis-HEP Does Automatically

From one YAML file, Jarvis-HEP will:
- Build the workflow
- Generate sampling points
- Run the external program
- Collect outputs
- Save results in a fixed layout

You do not manage loops or job control.

---

## Adapting This to Your Own Problem

To use Jarvis-HEP for a new task:

1. List the parameters you want to scan
2. Pick a sampling method
3. Treat your code as a black box
4. Declare inputs and outputs in YAML
5. Run Jarvis-HEP

If your task fits this list, it fits Jarvis-HEP.

---

## What to Read Next

- YAML keys and options:
  - [YAML Reference](yaml_reference.md)

- More examples:
  - [Examples](examples.md)

Read this tutorial once.
Then copy an example and modify it.
````