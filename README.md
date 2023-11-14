# WatershedParcellation
This package introduces a high performance Julia implementation of the parcellation method presented originally in Tim Laumann and Evan Gordon's 2016 paper ["Generation and Evaluation of a Cortical Area Parcellation from Resting-State Correlations."](https://pubmed.ncbi.nlm.nih.gov/25316338/), and based on their original MATLAB code. It builds upon recently registered packages [CIFTI.jl](https://github.com/myersm0/CIFTI.jl), [CorticalSurfaces.jl](https://github.com/myersm0/CorticalSurfaces.jl), and [CorticalParcels.jl](https://github.com/myersm0/CorticalParcels.jl).

Our results from this method, a neonatal parcellation generated from a dataset of 261 subjects, can be found [here](https://github.com/myersm0/Myers-Labonte_parcellation). The paper pre-print is available [here](https://www.biorxiv.org/content/10.1101/2023.11.10.566629v1).

The goals of this package are:

- Improve accessibility of this parcellation method to other researchers by:
	- Implementing it with open source tools
	- Substantially improving the execution speed
	- Bringing RAM usage under control (<= 32 GB)
	- Breaking the method into modular, customizable pieces that encourage further experimentation and improvement
- Extend the Julia language's ecosystem of fMRI-related packages
- Demonstrate utility of a framework of operations (laid out in [CorticalSurfaces.jl](https://github.com/myersm0/CorticalSurfaces.jl) and [CorticalParcels.jl](https://github.com/myersm0/CorticalParcels.jl)) for more easily and efficiently working with spatial and functional neuroimaging data in conjunction

## System requirements
The method requires handling several large matrices. To generate and evaluate a bilateral parcellational of the cortical surface in a space of 64,000 vertices:
- 32 GB RAM
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

Missing in this implementation is code to generate gradient maps from connectivity, which is a preliminary step for the above. However, that step mostly just involves a straightforward application of Connectome Workbench's [cifti-gradient](https://humanconnectome.org/software/workbench-command/-cifti-gradient) command, and the details of running that will be specific to your dataset and analysis strategy. The file `examples/make_gradients.jl` shows how you could do this for the single-subject case.

| |Goals|
|-|----------------------------------------------------|
|☑|Single-hemisphere functionality|
|☐|Bilateral functionality|
|☑|Exclusion of low-signal regions|
|☐|Support for loading of spatial data from GIFTI files|
|☑|A full demo|

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

See `examples/demo.jl` for a run-through of the major steps. You would need to edit the `config.json` in the same folder to point to paths of all the necessary input objects, to adjust parameters, etc. Unfortunately, due to the sizes of some of the required inputs, I'm unable to include these artifacts in this repo but I aim to find a way to get around that soon. Some excerpts below from `demo.jl` demonstrate the main functions available in this package. You may need to manually trigger the garbage collector with `GC.gc()` at times, if you're operating under RAM constraints.

### Edgemap creation
Supposing you have a matrix `grads` of size #vertices x #vertices, such as would be computed in `examples/make_gradients.jl`, and a `CorticalSurface` struct `c` to provide the vertex space and related spatial information (see [CorticalSurfaces.jl](https://github.com/myersm0/CorticalSurfaces.jl) for details):
```
minima = find_minima(grads, c)
edgemap = run_watershed(grads, minima, c)
```

### Parcel creation and cleanup
The edgemap then goes into a second pass of watershed, and a much faster one because this time it's operating on a #vertices vector rather than a large square matrix. This is where the initial parcellation is created:
```
labels = run_watershed(edgemap, c)
```

Then, because CorticalParcels.jl does not yet support working with a bilateral `Parcellation` struct, just pull out the left hemisphere for now and make a `Parcellation{Int}` struct from that:
```
hem = L
verts = @collapse vertices(c[hem])
px = Parcellation{Int}(c[hem], labels[verts])
edges = pad(edgemap[hem][:], c[hem])
```

The following are the cleanup operations then applied, to do things like merging small parcels and thresholding the parcel map at a certain edgemap height value. All of these are following the original methods from Gordon et al, but the intention of this package is to allow you a lot of freedom in omitting and/or modifying these steps or supplying your own:
```
remove_weak_boundaries!(px, edges; threshold = 0.3, radius = 30)
merge_small_parcels!(px, edges; minsize = 30, radius = 30)
threshold!(px, edges; threshold = 0.9)
merge_small_parcels!(px, edges; minsize = 30, radius = 30)
threshold!(px, edges; threshold = 0.9)
remove_articulation_points!(px; minsize = 4)
remove_small_parcels!(px; minsize = 10)
```

### Homogeneity evaluation
Rotate the parcellation around the surface by random x, y, z rotational parameters 1000 times to generate a null distribution of parcellations:
```
using NearestNeighbors

# make a KDTree that will assist in nearest neighbors search in the rotation process
tree = KDTree(coordinates(c[hem]); leafsize = 10)

# generate 1000 random 3x3x3 arrays that will be used to rotate the parcels
rotational_params = make_rotations(1000)

# now with the help of those two items, create a vector of 1000 rotated parcellations
# which will be compared against the real parcellation for homogeneity testing below
pxθ = rotation_wrapper(px, rotational_params, tree)
```

Now, load in a dconn (dense connectivity matrix file) with which you want to evaluate parcel homogeneity, and make a covariance of correlations matrix:
```
using CIFTI
using StatsBase: cov
dconn = CIFTI.load(config["dconn"])
cov_corr = make_cov_corr(dconn[L, L], c[L])
```

Now use it to test the parcel homogeneity of the real parcellation and of the 1000 rotations:
```
real_result = homogeneity_test(px, cov_corr; criteria = p -> default_criteria(p))
rotated_result = homogeneity_test(pxθ, cov_corr; criteria = p -> default_criteria(p))
```

Notice the keyword arg for which we're passing `p -> default_criteria(p)`. This is one of the strenghts of the current implementation: instead of the supplied `default_criteria` function, you could pass in a function specifying your own parcel inclusion criteria, possibly involving other objects in your workspace such as a map of low-signal regions that you want to exclude. Then, in any individual parcel undergoing a homogeneity test, upon receiving a `false` return value from this criteria function the result will be `NaN` and a homogeneity score will not be calculated. This is in order to exclude calculation of homogeneity for parcels that overlap with the medial wall, at a minimum, since functional measures (and therefore homogeneity) will not be defined within the medial wall. See the full demo for more details.

Finally, compare the real and the rotated results (null model) by taking the mean homogoneity of the real parcellation minus the mean of the same from all rotated parcellations, divided by the standard deviation of homogeneity from the rotated parcellations, to get a z-score:
```
zscore = summarize_homogeneity(real_result, rot_result)
```


[![Build Status](https://github.com/myersm0/WatershedParcellation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/myersm0/WatershedParcellation.jl/actions/workflows/CI.yml?query=branch%3Amain)
