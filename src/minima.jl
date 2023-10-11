
export find_minima

function find_minima(metric::AbstractMatrix, Aᵖ::AbstractMatrix, v::Int)
	neighbors = setdiff(findall(Aᵖ[v, :] .!= 0), v)
	a = repeat(metric[v, :]', outer = [length(neighbors), 1])
	b = @view metric[neighbors, :]
	return all(b - a .> 0; dims = 1)
end

"""
    find_minima(metric, A; power = 3)

Given a matrix `metric`, such as a gradient matrix, find local minima within a
window of size `radius` steps, guided by the topology given in adjaceny matrix `A`
"""
function find_minima(metric::AbstractMatrix, A::SparseMatrixCSC; power::Int = 3)
	nverts = size(metric, 1)
	minima = BitMatrix(undef, nverts, nverts)
	Threads.@threads for v in 1:nverts
		minima[v, :] = find_minima(metric, A^power, v)
	end
	return minima
end

function find_minima(metric::AbstractMatrix, hem::Hemisphere; power::Int = 3)
	haskey(hem, :A) || error("Hemisphere must contain adjacency matrix :A")
	m, n = size(metric)
	allequal([m, n]) || error("`metric` must be a square Matrix")
	if m == size(hem, Inclusive())
		return find_minima(metric, hem[:A]; power = power)
	elseif m == size(hem, Exclusive())
		return find_minima(metric, hem[:A, Exclusive()]; power = power)
	end
end


