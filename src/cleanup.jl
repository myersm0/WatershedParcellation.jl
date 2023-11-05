
export remove_weak_boundaries!, cap_at_height!, merge_small_parcels!
export remove_articulation_points!, remove_small_parcels!

function edge_strength(margin::Parcel, metric::Vector; radius::Number)
	verts = vertices(margin)
	median_edgeval = median(metric[verts])
	verts_in_radius =
		any(margin.surface[:distances][:, verts] .< radius; dims = 2)[:]
	localarea = metric[verts_in_radius]
	return mean(localarea .< median_edgeval)
end

"""
    remove_weak_boundaries!(px, metric; threshold)

For each marginal or interstitial region separating two parcels in `px::Parcellation`,
determine a merge priority for that region based on its edge strength: the median 
value from `metric` (such as an edgemap) relative to the values in a neighborhood of 
size `radius`. Iteratively merge pairs of parcels having such regions, starting with 
the weakest edge and proceeding until the edge strength reaches or exceeds `threshold`.

Returns the number of boundaries that have been merged.
"""
function remove_weak_boundaries!(
		px::Parcellation, metric::Vector; threshold = 0.38, radius = 30
	)
	n = 0
	while true
		margins = interstices(px)
		pairs = collect(keys(margins))
		vals = zeros(length(pairs))
		Threads.@threads :dynamic for i in 1:length(pairs)
			pair = pairs[i]
			vals[i] = edge_strength(
				Parcel(px.surface, margins[pair]), metric; radius = radius
			)
		end
		i = argmin(vals)
		vals[i] < threshold || break
		merge!(px, pairs[i]...)
		n += 1
	end
	return n
end

"""
    merge_small_parcels!(px, metric; threshold)

For each `Parcel` in `px::Parcellation` smaller than `minsize` vertices, find its
neighboring parcels (if any) and merge it with the one that has the weakest edge,
i.e. the one for which the median value of `metric` (such as an edgemap) in the 
interstitial region has the lowest value relative to values in a neighborhood
of size `radius`.

Returns the number of merge operations that have occurred. 
"""
function merge_small_parcels!(
		px::Parcellation, metric::Vector; minsize = 30, radius = 30
	)
	n = 0
	for k in keys(px)
		# skip if this key has already been removed (merged) in a previous iter:
		haskey(px, k) || continue
		p = px[k]
		size(p) < minsize || continue
		neigh_parcels = filter(
			k -> any(interstices(p, px[k])),
			setdiff(collect(keys(px)), k)
		)
		length(neigh_parcels) > 0 || continue
		vals = [
			edge_strength(Parcel(px.surface, interstices(p, px[k])), metric; radius = radius) 
			for k in neigh_parcels
		]
		merge!(px, k, neigh_parcels[argmin(vals)])
		n += 1
	end 
	return n
end

"""
    cap_at_height!(px, metric; threshold)

Remove vertices from `px::Parcellation` where `metric` (such as an edgemap) exceeds
a maximum value, defined as the `threshold` quantile of values from `metric`. If any
parcels are separated or disconnected in this process, then split their new connected
components that emerged into new parcels.

Returns the number of high vertices removed in this process.
"""
function cap_at_height!(px::Parcellation, metric::Vector; threshold = 0.9)
	0.0 < threshold <= 1.0 || error(DomainError)
	threshold = quantile(metric, [threshold])[1]
	n = 0
	for k in keys(px)
		p = px[k]
		verts = vertices(p)
		high_edge_verts = filter(x -> metric[x] > threshold, verts)
		any(high_edge_verts) || continue
		new_parcels = split(p, high_edge_verts)
		if length(new_parcels) == 1
			intersect!(p, new_parcels[1])
		else
			delete!(px, k)
			for i in 1:length(new_parcels)
				verts = vertices(new_parcels[i])
				new_key = maximum(keys(px)) + 1
				px[new_key] = new_parcels[i]
			end
		end
		n += length(high_edge_verts)
	end
	return n
end

"""
    remove_articulation_points!(px; threshold)

If any parcels in `px::Parcellation` have an articulation point or cut vertex,
the removal of which will disconnect the parcel, then remove that vertex and 
separate the parcel into its resulting connected components, as long as there
are at least two "reasonably sized" components of size >= `minsize` vertices.

Returns the number of articulation points that were handled.
"""
function remove_articulation_points!(px::Parcellation; minsize::Int = 4)
	n = 0
	for k in keys(px)
		new_parcels = cut(px[k])
		if sum(size.(new_parcels) .> minsize) > 1
			println("Handing parcel $k ...")
			delete!(px, k)
			for p in new_parcels
				new_key = maximum(keys(px)) + 1
				px[new_key] = p
			end
			n += 1
		end
	end
	return n
end

"""
    remove_articulation_points!(px; threshold)

Remove any parcels in `px::Parcellation` smaller than `minsize` vertices.

Returns the number of small parcels that were removed.
"""
function remove_small_parcels!(px::Parcellation; minsize = 10)
	n = 0
	for k in keys(px)
		size(px[k]) < minsize || continue
		delete!(px, k)
		n += 1
	end
	return n
end

