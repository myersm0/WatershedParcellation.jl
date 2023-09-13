
using Cifti
using SparseArrays

function make_adjmat(neigh::VertexList)
	adjmat = spzeros(Bool, nverts, nverts)
	for v in 1:nverts
		adjmat[v, v] = true
		adjmat[v, neigh[v]] .= true
	end
	return adjmat
end

function load_gradients(filename::String)
	grads = # TODO: load smoothed gradient here
	return SharedArray(Matrix(grads))
end

function load_neighbors()
	return load("$assets_dir/neighbors.jld", "neigh")
end

function find_minima(metric::AbstractMatrix, reachability::AbstractMatrix, v::Int)::BitMatrix
	neighbors = setdiff(findall(reachability[v, :] .!= 0), v)
	a = repeat(metric[v, :]', outer = [length(neighbors), 1])
	b = @view metric[neighbors, :]
	return all(b - a .> 0; dims = 1)
end

function find_minima(grads::AbstractMatrix, adjmat::SparseMatrixCSC)
	reachability = adjmat ^ 3 # for each vert, identify 3-step-reachable vertices
	minima = BitMatrix(undef, nverts, nverts)
	Threads.@threads for v in 1:nverts
		minima[v, :] = find_minima(grads, reachability, v)
	end
	return minima
end



