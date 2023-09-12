
using Cifti
using Statistics
using HDF5
using JLD
using SparseArrays
using SharedArrays

data_dir = "./test/data/"
grads = # load smoothed gradient here
grads = SharedArray(Matrix(grads'))
GC.gc()

nverts = size(grads, 1)
neigh = load("./32k_tools/neighbors.jld", "neigh")

function make_adjmat(neigh::Vector{Vector{Int}})::SpareMatrixCSC
	adjmat = spzeros(Bool, nverts, nverts)
	for i in 1:nverts
		adjmat[i, i] = true
		adjmat[i, neigh[i]] .= true
	end
	return adjmat
end

adjmat = make_adjmat(neigh)

# let reachability be a nverts x nverts matrix where, in each row, 
# the non-zero elements will be the set of vertices reachable within 3 steps
reachability = adjmat ^ 3

function find_minima(metric::AbstractMatrix, reachability::AbstractMatrix, v::Int)::BitMatrix
	neighbors = setdiff(findall(reachability[v, :] .!= 0), v)
	a = repeat(metric[v, :]', outer = [length(neighbors), 1])
	b = @view metric[neighbors, :]
	all(b - a .> 0; dims = 1)
end

minima = BitMatrix(undef, nverts, nverts)
Threads.@threads for v in 1:nverts
	minima[v, :] = find_minima(grads, reachability, v)
end



