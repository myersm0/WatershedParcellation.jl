
using SparseArrays
using JLD

export find_minima

function find_minima(metric::AbstractMatrix, Aᵖ::AbstractMatrix, v::Int)
	neighbors = setdiff(findall(Aᵖ[v, :] .!= 0), v)
	a = repeat(metric[v, :]', outer = [length(neighbors), 1])
	b = @view metric[neighbors, :]
	return all(b - a .> 0; dims = 1)
end

function find_minima(metric::AbstractMatrix, A::SparseMatrixCSC; power::Int = 3)
	nverts = size(grads, 1)
	minima = BitMatrix(undef, nverts, nverts)
	Threads.@threads for v in 1:nverts
		minima[v, :] = find_minima(metric, A^power, v)
	end
	return minima
end



