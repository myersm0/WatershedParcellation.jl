# WatershedParcellation
This package introduces a high performance Julia implementation of the parcellation method presented originally in Tim Laumann and Evan Gordon's 2016 paper ["Generation and Evaluation of a Cortical Area Parcellation from Resting-State Correlations."](https://pubmed.ncbi.nlm.nih.gov/25316338/), and based on their original MATLAB code. It builds upon recently registered packages [CorticalSurfaces.jl](https://github.com/myersm0/CorticalSurfaces.jl) and [CorticalParcels.jl](https://github.com/myersm0/CorticalParcels.jl).

Our new, soon-to-be published results from this method, a neonatal parcellation generated from a dataset of 261 subjects, can be found [here](https://github.com/myersm0/Myers-Labonte_parcellation). A link to the paper on BioArxiv will be available shortly.

The goals of this package (and its related dependencies, menioned above) are:
- Extend the Julia language's ecosystem of fMRI-related packages
- Build a framework of operations for more easily and performantly working with spatial and functional neuroimaging data in conjunction
- Improve accessibility of this parcellation method to other researchers by:
	- Implementing it with open source tools
	- Substantially improving the execution speed
	- Bringing RAM usage under control (<= 32 GB)
	- Breaking the method into modular, easy-to-digest pieces for better comprehensibility and customizability

## System requirements
The method requires handling several large matrices. To generate and evaluate a bilateral parcellational of the cortical surface in a space of 64,000 vertices:
- 32 GB RAM (required)
- 4 to 8 CPU cores (recommended)
- Julia >= v1.9

For working with a single-hemisphere only, you should be able to reduce the RAM usage to 12 GB.

## Development status
In terms of the basic stages of operations, the below chart summarizes what's currently ready for use:

| |Functionality|
|-|:----------------------------------------------------------------|
|☑|The core watershed algorithm to generate edge maps from gradients|
|☑|Flooding of edge map basins to generate parcels|
|☑|Homogeneity evaluation and null-model testing via rotation|

Missing in this implementation is code to generate gradient maps from connectivity, which is a preliminary step for the above. However, doing that step mostly just involves a straightforward application of Connectome Workbench's [cifti-gradient](https://humanconnectome.org/software/workbench-command/-cifti-gradient) command, and the details of running that will be specific to your dataset and analysis strategy. The file `examples/make_gradients.jl` demonstrates how you could do this for the single-subject case.

As for the functionality that currently exists, there are a few gaps still to be filled and that will be remedied in coming days (by approximately the end of October 2023).
| |Goals|
|-|----------------------------------------------------|
|☑|Single-hemisphere functionality|
|☐|Bilateral functionality|
|☑|Exclusion of low-signal regions|
|☐|Support for loading of spatial data from GIFTI files|
|☐|A full demo|

## Installation
From within a Julia session:
```
using Pkg
Pkg.add("WatershedParcellation")
```

## Performance
Approximate execution times noted below were achieved in running the methods on a single hemisphere in 32k resolution, on a Macbook Pro with 8 Apple M2 cores.
|Processing stage|Benchmark|
|-----------------------------------------|---------|
|Finding local minima in the gradient maps|4 minutes|
|Generation of an edge map from gradients|3 minutes|
|Generation of 1000 rotated parcellations|7 seconds|
|Homogeneity testing of 1000 parcellations|90 seconds|

## Usage
Before loading Julia, you need to inform it of the number of available processing cores by setting an environment variable, for example `export JULIA_NUM_THREADS=8` in bash. Then, within Julia:
```
using CorticalSurfaces
using CorticalParcels
using WatershedParcellation
```

### Running the watershed algorithm to produce an edge map
Details coming soon.

### Flooding basins in the edge map to produce parcels
Not yet implemented. Check back soon.

### Generating a rotation-based null model
Details coming soon.

### Homogeneity evalaution
Details coming soon.

[![Build Status](https://github.com/myersm0/WatershedParcellation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/myersm0/WatershedParcellation.jl/actions/workflows/CI.yml?query=branch%3Amain)
