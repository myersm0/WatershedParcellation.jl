
export find_minima

function find_minima(metric::AbstractMatrix, Aᵖ::AbstractMatrix, v::Int)::BitMatrix
	neighbors = setdiff(findall(Aᵖ[:, v] .!= 0), v)
	a = repeat(metric[v, :]', outer = [length(neighbors), 1])
	b = @view metric[neighbors, :]
	return all(b - a .> 0; dims = 1)
end

function find_minima(
		metric::AbstractMatrix, A::AbstractMatrix; radius::Int = 3
	)::BitMatrix
	nverts = size(metric, 1)
	Aᵖ = A ^ radius
	minima = BitMatrix(undef, nverts, nverts)
	Threads.@threads :dynamic for v in 1:nverts
		minima[:, v] = find_minima(metric, Aᵖ, v)
	end
	return minima'
end

"""
    find_minima(metric, surface; radius = 3)

Given a matrix `metric`, such as a gradient matrix, find local minima within a
window of size `radius` steps. Argument `surface <: SurfaceSpace` must contain 
symbol `:A` (an adjacency matrix)
"""
function find_minima(metric::AbstractMatrix, surface::SurfaceSpace; radius::Int = 3)
	haskey(surface, :A) || error("surface must contain adjacency matrix :A")
	m, n = size(metric)
	m == n || error("`metric` must be a square Matrix")
	if m == size(surface, Inclusive())
		A = surface[:A]
	elseif m == size(surface, Exclusive())
		A = surface[:A, Exclusive()]
	else
		error(DimensionMismatch)
	end
	return find_minima(metric, A; radius = radius)
end


