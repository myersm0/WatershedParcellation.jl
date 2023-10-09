
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
using Statistics: cov, mean
using StatsBase: sample
using ThreadsX

const assets_dir = joinpath(@__DIR__, "..", "data")

include("minima.jl")
include("watershed.jl")
include("rotation.jl")
include("homogeneity.jl")

end

