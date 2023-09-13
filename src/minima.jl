
using Cifti
using SparseArrays

function make_adjmat(neigh::Vector{Vector{Int}})::SpareMatrixCSC
	adjmat = spzeros(Bool, nverts, nverts)
	for i in 1:nverts
		adjmat[i, i] = true
		adjmat[i, neigh[i]] .= true
	end
	return adjmat
end

function find_minima(metric::AbstractMatrix, reachability::AbstractMatrix, v::Int)::BitMatrix
	neighbors = setdiff(findall(reachability[v, :] .!= 0), v)
	a = repeat(metric[v, :]', outer = [length(neighbors), 1])
	b = @view metric[neighbors, :]
	all(b - a .> 0; dims = 1)
end

function load_gradients(filename::String)
	grads = # TODO: load smoothed gradient here
	return SharedArray(Matrix(grads))
end

function load_neighbors()
	load("./32k_tools/neighbors.jld", "neigh")
end

function make_minima_metrics(grads::AbstractMatrix, adjmat::SparseMatrixCSC)
	# let reachability be a nverts x nverts matrix where, in each row, 
	# the non-zero elements will be the set of vertices reachable within 3 steps
	reachability = adjmat ^ 3

	minima = BitMatrix(undef, nverts, nverts)
	Threads.@threads for v in 1:nverts
		minima[v, :] = find_minima(grads, reachability, v)
	end
	return minima
end


