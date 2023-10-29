
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
using StatsBase: cov, mean, median, quantile, sample
using ThreadsX

const assets_dir = joinpath(@__DIR__, "..", "data")

import CorticalSurfaces: AdjacencyList, AdjacencyMatrix, DistanceMatrix

include("minima.jl")
include("watershed.jl")
include("fill.jl")
include("rotation.jl")
include("homogeneity.jl")

end

