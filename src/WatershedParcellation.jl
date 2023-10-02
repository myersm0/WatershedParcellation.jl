
module WatershedParcellation

using Chain
using CIFTI
using CorticalSurfaces
using CorticalParcels
using JLD
using LinearAlgebra
using NamedArrays
using NearestNeighbors
using SparseArrays
using Statistics
using StatsBase: sample
using ThreadsX

const assets_dir = joinpath(@__DIR__, "../data/")

include("generation/minima.jl")
include("generation/watershed.jl")
include("evaluation/rotation.jl")
include("evaluation/homogeneity.jl")

end

