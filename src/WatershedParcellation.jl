
module WatershedParcellation

const nverts = 59412
const VertexList = Vector{Vector{UInt16}}
const assets_dir = joinpath(@__DIR__, "../data/32k_tools/")

include("generation/minima.jl")
include("generation/watershed.jl")

end

