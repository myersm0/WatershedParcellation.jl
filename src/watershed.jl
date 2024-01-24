
# doesn't return a val, but modifies labels and watershed_zones vectors in place
function eval_at_height!(h, labels, edgemetric, watershed_zones, neighbors)
	nodes_at_threshold = @chain begin
		(edgemetric .< h) .& (labels .== 0) .& .!watershed_zones
		findall
		sample(_, length(_); replace = false)
	end
	for node in nodes_at_threshold
		nodeneighlab = setdiff(labels[neighbors[node]], 0)
		n = length(nodeneighlab)
		n > 0 || continue
		if n == 1
			labels[node] = nodeneighlab[1]
		else
			watershed_zones[node] = true
			labels[node] = 0
		end
	end
end

function watershed_iter(edgemetric, minima, neighbors, heights, nverts)::BitVector
	watershed_zones = falses(nverts)
	labels = zeros(UInt16, nverts)
	nbasins = sum(minima)
	labels[minima] .= sample(0x0001:nbasins, nbasins; replace = false)
	for h in heights
		eval_at_height!(h, labels, edgemetric, watershed_zones, neighbors)
	end
	return labels .== 0
end

"""
    run_watershed(metric, minima, neighbors)

Run the watershed algorithm on `metric`, using its `minima` (as returned from 
`find_minima()`) as initialization points and guided by the topology given in
the adjacency list `neighbors`.

Some optional parameters can be tuned:
- `nsteps`: the number of bins into which to discretize the heights (default 400)
- `fracmaxh`: a scaling factor to determine a ceiling on height values (default 1.0)
"""
function run_watershed(
		metric::Matrix, minima::BitMatrix, neighbors::AdjacencyList;
		nsteps::Int = 400, fracmaxh::Float64 = 1.0
	)
	all([nsteps, fracmaxh] .> 0) || error(DomainError)
	minheight = minimum(metric)
	maxheight = maximum(metric) * fracmaxh
	heights = range(minheight, maxheight, length = nsteps)
	nverts = length(neighbors)
	edges = falses(nverts, nverts)
	Threads.@threads :dynamic for v in 1:nverts
		edges[:, v] .= watershed_iter(metric[:, v], minima[:, v], neighbors, heights, nverts) 
	end
	return mean(edges; dims = 2)[:]
end

function run_watershed(
		metric::Matrix, minima::BitMatrix, surface::SurfaceSpace; kwargs...
	)
	haskey(surface, :neighbors) || error("surface must have adjacency list :neighbors")
	nverts = size(metric, 1)
	if nverts == size(surface)
		mw_indexing = Inclusive()
	elseif nverts == size(surface, Exclusive())
		mw_indexing = Exclusive()
	else
		error(DimensionMismatch)
	end
	neighbors = surface[:neighbors, mw_indexing]
	return run_watershed(metric, minima, neighbors; kwargs...)
end

function make_basins!(metric::Vector, surface::SurfaceSpace; thresh_quantile = 0.75)
	haskey(surface, :neighbors) || error(KeyError)
	 nelem = length(metric)
	if nelem == size(surface)
		mw_indexing = Inclusive()
	elseif nelem == size(surface, Exclusive())
		mw_indexing = Exclusive()
	else
		error(DimensionMismatch)
	end
	neighbors = surface[:neighbors, mw_indexing]
	nverts = size(neighbors, 1)
	threshold = quantile(metric, [thresh_quantile])[1]
	out = falses(nverts)
	for v in 1:nverts
		val = metric[v]
		nodeneigh = neighbors[v]
		min_neigh = minimum(metric[nodeneigh])
		if val <= min_neigh && val < threshold
			out[v] = true
			temp = findall(metric[nodeneigh] .== val)
			if length(temp) > 0
				metric[nodeneigh[temp]] .+= 0.00001
			end
		end
	end
	return out
end

function run_watershed(metric::Vector, surface::SurfaceSpace; thresh_quantile = 0.75)
	nelem = length(metric)
	if nelem == size(surface)
		mw_indexing = Inclusive()
	elseif nelem == size(surface, Exclusive())
		mw_indexing = Exclusive()
	else
		error(DimensionMismatch)
	end
	neighbors = surface[:neighbors, mw_indexing]
	metric2 = deepcopy(metric)
	basins = make_basins!(metric2, surface; thresh_quantile = 0.75)
	labels = zeros(Int, nelem)
	labelpos = findall(basins)
	nlabels = length(labelpos)
	sorti = sortperm(metric[labelpos])
	labels[labelpos[sorti]] .= 1:nlabels
	watershed_zones = falses(nelem)
	hiter = sort(unique(metric))
	[eval_at_height!(h, labels, metric2, watershed_zones, neighbors) for h in hiter]
	return labels
end

