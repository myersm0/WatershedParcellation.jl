
module WatershedParcellation

const nverts = 59412
const VertexList = Vector{Vector{UInt16}}

include("minima.jl")
include("watershed.jl")

end
