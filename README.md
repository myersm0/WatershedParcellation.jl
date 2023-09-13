# WatershedParcellation

Under development. I aim to release my code here over the next several days and weeks, both in the form of a Julia module and as a standalone executable, as adapted from the original MATLAB code by Tim Laumann and Evan Gordon from their 2016 paper ["Generation and Evaluation of a Cortical Area Parcellation from Resting-State Correlations."](https://pubmed.ncbi.nlm.nih.gov/25316338/). 

## Resource requirements
The whole process requires passing around several large matrices. Care has been taken to reduce RAM overhead, by using sparse representations and low-storage element types where possible. About 40 GB should be sufficient to run the whole thing.

Availability of 8 CPU cores is recommended.

## Roadmap
1. The core watershed algorithm to generate edge maps from gradients (status: almost ready to use)
2. Homogeneity evaluation and null-model testing via rotation (ETA: week of Sept 17th)
3. Code to get from connectivity data to gradients (ETA: week of Sept 24th)

The basic sequence of operations is connectivity -> gradients -> watershed -> homogeneity testing. However, the creation of gradients has several external dependencies and process complexities so I'm postponing the release of that step.

## Usage
### Running the watershed algorithm to produce an edge map
Before loading Julia, you need to inform it of the number of available processing cores by setting an environment variable, for example `export JULIA_NUM_THREADS=8` in bash. Then, within Julia:

```
grads = load_gradients(filename) # not yet implemented
neigh = load_neighbors()
adjmat = make_adjmat(neigh)
minima = find_minima(grads, adjmat)
edgemap = run_watershed(grads, minima, neigh)
```

[![Build Status](https://github.com/myersm0/WatershedParcellation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/myersm0/WatershedParcellation.jl/actions/workflows/CI.yml?query=branch%3Amain)
