# WatershedParcellation

Under development. I aim to release my code here over the next several days and weeks, both in the form of a Julia module and as a standalone executable, as adapted from the original MATLAB code by Tim Laumann and Evan Gordon from their 2016 paper ["Generation and Evaluation of a Cortical Area Parcellation from Resting-State Correlations."](https://pubmed.ncbi.nlm.nih.gov/25316338/).

Our soon-to-be published results, a neonatal parcellation generated from a dataset of 262 subjects, can be found ["here"](https://github.com/myersm0/Myers-Labonte_parcellation).

## Resource requirements
The whole process requires passing around several large matrices. Care has been taken to reduce the RAM footprint by using sparse representations and small element types where possible. 32 GB RAM should be sufficient to run all the code.

Availability of 8 CPU cores is recommended. In the current state of implementation, the benefits of multithreaded parallelism are saturated at about 8 cores so I don't recommend going higher than that.

## Roadmap
1. The core watershed algorithm to generate edge maps from gradients (status: almost ready to use)
2. Homogeneity evaluation and null-model testing via rotation (ETA: week of Sept 17th)
3. Code to get from connectivity data to gradients (ETA: week of Sept 24th)

The basic sequence of operations is connectivity -> gradients -> watershed -> homogeneity testing. However, the creation of gradients has several external dependencies and process complexities so I'm postponing the release of that step.

## Installation
I aim to add this code to the Julia general repository. Until then, installation is as follows:
```
using Pkg
Pkg.add(url = "https://github.com/myersm0/WatershedParcellation.jl")
```

Eventually I aim to add the option to run a standalone executable so that you don't have to know Julia or manage a Julia session yourself in order to generate results.

## Usage
### Running the watershed algorithm to produce an edge map
Before loading Julia, you need to inform it of the number of available processing cores by setting an environment variable, for example `export JULIA_NUM_THREADS=8` in bash. Then, within Julia:

```
using WatershedParcellation

grads = load_gradients(filename) # not yet implemented
neigh = load_neighbors() # 5 seconds
adjmat = make_adjmat(neigh) # 5 seconds
minima = find_minima(grads, adjmat) # 4 minutes
edgemap = run_watershed(grads, minima, neigh) # 30 minutes
```

Approximate execution times noted above were achieved on a CentOS 7 Linux machine using 8 cores (Intel(R) Xeon(R) Gold 6154 CPU @ 3.00GHz).

[![Build Status](https://github.com/myersm0/WatershedParcellation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/myersm0/WatershedParcellation.jl/actions/workflows/CI.yml?query=branch%3Amain)
