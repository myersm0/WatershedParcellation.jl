# WatershedParcellation

Under development. I aim to release my code here over the next several days and weeks, both in the form of a Julia module and as a standalone executable, as adapted from the original MATLAB code by Tim Laumann and Evan Gordon from their 2016 paper ["Generation and Evaluation of a Cortical Area Parcellation from Resting-State Correlations."](https://pubmed.ncbi.nlm.nih.gov/25316338/).

Our soon-to-be published results, a neonatal parcellation generated from a dataset of 262 subjects, can be found [here](https://github.com/myersm0/Myers-Labonte_parcellation).

## System requirements
32 GB RAM (required), 8 CPU cores (recommended).

## Roadmap
1. The core watershed algorithm to generate edge maps from gradients (status: almost ready to use)
2. Homogeneity evaluation and null-model testing via rotation (ETA: week of Sept 17th)
3. Code to get from connectivity data to gradients (ETA: week of Sept 24th)

The basic order of operations is gradients -> watershed -> homogeneity testing. But the creation of gradients has several external dependencies and process complexities so I'm postponing the release of that step.

## Installation
I aim to submit this package to the Julia general registry. Until then, to install run the following from within a Julia session:
```
using Pkg
Pkg.add(url = "https://github.com/myersm0/WatershedParcellation.jl")
```

## Usage
### Running the watershed algorithm to produce an edge map
Before loading Julia, you need to inform it of the number of available processing cores by setting an environment variable, for example `export JULIA_NUM_THREADS=8` in bash. Then, within Julia:

```
using WatershedParcellation

grads = load_gradients(filename) # not yet implemented
neigh = load_neighbors() # 5 seconds
adjmat = make_adjmat(neigh) # 30 seconds
minima = find_minima(grads, adjmat) # 4 minutes
edgemap = run_watershed(grads, minima, neigh) # 30 minutes
```

Approximate execution times noted above were achieved on a CentOS 7 Linux machine using 8 cores (Intel(R) Xeon(R) Gold 6154 CPU @ 3.00GHz).

[![Build Status](https://github.com/myersm0/WatershedParcellation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/myersm0/WatershedParcellation.jl/actions/workflows/CI.yml?query=branch%3Amain)
