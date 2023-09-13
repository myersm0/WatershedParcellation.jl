# WatershedParcellation

Under development. I aim to release my code here over the next several days and weeks, both in the form of a Julia module and as a standalone executable, as adapted from the original MATLAB code by Tim Laumann and Evan Gordon from their 2016 paper ["Generation and Evaluation of a Cortical Area Parcellation from Resting-State Correlations."](https://pubmed.ncbi.nlm.nih.gov/25316338/)

Roadmap for coming days and weeks:
1. The core watershed algorithm to generate edgemaps from gradients (status: almost ready to use)
2. Homogeneity evaluation and null-model testing via rotation (ETA: week of Sept 17th)
3. Code to get from connectivity data to gradients (ETA: week of Sept 24th)
  - mostly this involves wrapping some algorithms from ["Connectome Workbench"](https://humanconnectome.org/software/workbench-command))

The basic sequence of operations is connectivity -> gradients -> watershed -> homogeneity testing. However, the creation of gradients has several external dependencies and process complexities so I'm postponing the release of the step.

[![Build Status](https://github.com/myersm0/WatershedParcellation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/myersm0/WatershedParcellation.jl/actions/workflows/CI.yml?query=branch%3Amain)
