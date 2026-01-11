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

To study a Beyond Standard Model numerically, one typically needs to scan its parameter space. Most scan problems look like this:

- Pick some parameters
- Explore their values with a sampler
- Run an external program at each point
- Save the results

Jarvis-HEP is built for this pattern. 

### Eggbox Function as an Simplifed Package

We begin with the Eggbox function, a purely mathematical function.

The Eggbox function here is defined as a two-dimensional function of two real variables,
denoted by `x` and `y`.

`z(x, y) = (sin(x) * cos(y) + 2)^5`

For each input point `(x, y)`, the function returns a single scalar value `z`. In a realistic physics study, `z` could represent any theoretical prediction or observable, for example the Higgs boson mass in the MSSM. Without a numerical scan, it is hard to see how such an observable behaves across the parameter space.

![Eggbox surface](includes/eggbox_vis.png)

As shown in the above figure, a three-dimensional surface plot of `z(xx, yy)` clearly illustrates this structure. This function has several characteristic features:
- It is smooth and continuous
- It exhibits a highly structured, multi-peak landscape
- Many local maxima and minima appear across the parameter space

These features make the Eggbox function a standard test case
for sampling algorithms and scanning workflows.

---

Suppose an experiment measures this observable as

`z = 100 ± 10`.

Our task is then to identify which regions in the `(x, y)` parameter space
are compatible with this measurement.

To make this numerical study practical, we implement a small package
that computes `z` for any given `(x, y)`.

---
#### EggBox package （Code see [eggbox](../External/Inertial/EggBox)）


We now move from the mathematical description
to the computational abstraction used by Jarvis-HEP.

In Jarvis-HEP, an *input point* refers to one complete set of parameter values.
For the Eggbox example, one input point is the pair `(xx, yy)`.

The Eggbox calculation is not treated as a direct function call.
Instead, it is packaged as a small external program,
referred to as a **calculator**.

The calculator follows a strict input–output contract:

1. Jarvis-HEP writes the current input point to `input.json`
2. Jarvis-HEP executes the calculator program
3. The calculator reads `input.json` and computes the Eggbox value
4. The calculator writes results to `output.json`
5. Jarvis-HEP reads `output.json` and records the requested outputs

In this example, the calculator writes:
- `z`: the value of the Eggbox function
- `Time`: a dummy field used to mimic runtime information

From Jarvis-HEP’s perspective, the calculator is a complete black box.
Jarvis-HEP does not inspect the internal code or the mathematical expression.

This example demonstrates the core abstraction used by Jarvis-HEP:
any function or model can be used,
as long as it can be packaged as a program
with well-defined inputs and outputs.
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