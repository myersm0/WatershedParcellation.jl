
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
		px::AbstractParcellation, metric::Vector; threshold = 0.3, radius = 30
	)
end

function remove_weak_boundaries!(
		px::HemisphericParcellation, metric::Vector; threshold = 0.3, radius = 30
	)
	n = 0
	margins = interstices(px)

	while true
		# find boundaries to merge based on edge_strength
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

		# delete margins that involved this pair because they may no longer be accurate
		for pair in pairs
			if pairs[i][1] in pair || pairs[i][2] in pair
				delete!(margins, pair)
			end
		end

		# add new margins for the retained parcel because it has grown
		k = pairs[i][1]
		ua = unassigned(px)
		for k2 in setdiff(collect(keys(px)), k)
			i = interstices(px[k], px[k2]) .& ua
			if any(i)
				a = minimum([k, k2])
				b = maximum([k, k2])
				margins[(a, b)] = i
			end
		end
		n += 1
	end

	return n
end

function remove_weak_boundaries!(
		px::BilateralParcellation, metric::Vector; threshold = 0.3, radius = 30
	)
	return remove_weak_boundaries!(px[L]) + remove_weak_boundaries!(px[R])
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
		px::AbstractParcellation, metric::Vector; minsize = 30, radius = 30
	)
end

function merge_small_parcels!(
		px::HemisphericParcellation, metric::Vector; minsize = 30, radius = 30
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

function merge_small_parcels!(
		px::BilateralParcellation, metric::Vector; minsize = 30, radius = 30
	)
	return merge_small_parcels!(px[L]) + merge_small_parcels!(px[R])
end

"""
    threshold!(px, metric; threshold)

Remove vertices from `px::Parcellation` where `metric` (such as an edgemap) exceeds
a maximum value, defined as the `threshold` quantile of values from `metric`. If any
parcels are separated or disconnected in this process, then split their new connected
components that emerged into new parcels.

Returns the number of high vertices removed in this process.
"""
function threshold!(px::AbstractParcellation, metric::Vector; threshold = 0.9) end

function threshold!(px::HemisphericParcellation, metric::Vector; threshold = 0.9)
	0.0 < threshold <= 1.0 || error(DomainError)
	threshold = quantile(metric[isfinite.(metric)], [threshold])[1]
	n = 0
	for k in keys(px)
		p = px[k]
		verts = vertices(p)
		high_edge_verts = filter(x -> metric[x] > threshold, verts)
		length(high_edge_verts) > 0 || continue
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

function threshold!(px::BilateralParcellation, metric::Vector; threshold = 0.9)
	return threshold!(px[L]) + threshold!(px[R])
end

"""
    remove_articulation_points!(px; threshold)

If any parcels in `px::Parcellation` have an articulation point or cut vertex,
the removal of which will disconnect the parcel, then remove that vertex and 
separate the parcel into its resulting connected components, as long as there
are at least two "reasonably sized" components of size >= `minsize` vertices.

Returns the number of articulation points that were handled.
"""
function remove_articulation_points!(px::AbstractParcellation; minsize::Int = 4) end

function remove_articulation_points!(px::HemisphericParcellation; minsize::Int = 4)
	n = 0
	for k in keys(px)
		new_parcels = cut(px[k])
		if sum(size.(new_parcels) .> minsize) > 1
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

function remove_articulation_points!(px::BilateralParcellation; minsize::Int = 4)
	return remove_articulation_points!(px[L]) + remove_articulation_points!(px[R])
end

"""
    remove_articulation_points!(px; threshold)

Remove any parcels in parcellation `px` smaller than `minsize` vertices.

Returns the number of small parcels that were removed.
"""
function remove_small_parcels!(px::AbstractParcellation; minsize = 10) end

function remove_small_parcels!(px::HemisphericParcellation; minsize = 10)
	n = 0
	for k in keys(px)
		size(px[k]) < minsize || continue
		delete!(px, k)
		n += 1
	end
	return n
end

function remove_small_parcels!(px::BilateralParcellation; minsize = 10)
	return remove_small_parcels!(px[L]) + remove_small_parcels!(px[R])
end

