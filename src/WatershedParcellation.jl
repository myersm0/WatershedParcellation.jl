
module WatershedParcellation

using Chain
using CIFTI
using CorticalSurfaces
using CorticalParcels
using JLD
using HDF5
using LinearAlgebra
using NamedArrays
using NearestNeighbors
using SparseArrays
using StatsBase: cov, mean, median, std, quantile, sample
using ThreadsX

const assets_dir = joinpath(@__DIR__, "..", "data")

import CorticalSurfaces: AdjacencyList, AdjacencyMatrix, DistanceMatrix

include("minima.jl")
export find_minima

include("watershed.jl")
export run_watershed

include("cleanup.jl")
export remove_weak_boundaries!, threshold!, merge_small_parcels!
export remove_articulation_points!, remove_small_parcels!

include("rotation.jl")
export make_rotations, rotation_wrapper

include("homogeneity.jl")
export make_cov_corr, homogeneity_test, default_criteria, summarize_homogeneity

end

