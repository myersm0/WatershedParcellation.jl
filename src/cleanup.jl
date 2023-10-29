
export remove_weak_boundaries!, remove_high_edges!, merge_small_parcels!

function localarea_test(margin::Parcel, metric::Vector; radius = 30)
	verts = vertices(margin)
	median_edgeval = median(metric[verts])
	verts_in_radius =
		any(margin.surface[:distances][:, verts] .< radius; dims = 2)[:]
	localarea = metric[verts_in_radius]
	return mean(localarea .< median_edgeval)
end

function remove_weak_boundaries!(px::Parcellation, metric::Vector; threshold = 0.38)
	n = 0
	while true
		margins = interstices(px)
		pairs = collect(keys(margins))
		vals = [localarea_test(Parcel(px.surface, margins[p]), metric) for p in pairs]
		i = argmin(vals)
		vals[i] < threshold || break
		merge!(px, pairs[i]...)
		n += 1
	end
	return n
end

function merge_small_parcels!(px::Parcellation, metric::Vector; threshold = 30)
	n = 0
	for k in keys(px)
		k in keys(px) || continue
		p = px[k]
		size(p) < threshold || continue
		neigh_parcels = filter(k -> sum(interstices(p, px[k])) > 0, collect(keys(px)))
		length(neigh_parcels) > 0 || continue
		vals = [localarea_test(Parcel(px.surface, interstices(p, px[k])), metric) for k in neigh_parcels]
		merge!(px, k, neigh_parcels[argmin(vals)])
		n += 1
	end 
	return n
end

function remove_high_edges!(px::Parcellation, metric::Vector; threshold = 0.336085307)
	for k in keys(px)
		p = px[k]
		verts = vertices(p)
		high_edge_verts = filter(x -> metric[x] > threshold, verts)
		new_parcels = split(p, high_edge_verts)
		if length(new_parcels) == 1
			intersect!(p, new_parcels[1])
		else
			delete!(px, k)
			for i in 1:length(new_parcels)
				verts = vertices(new_parcels[i])
				new_k = verts[1]
				px[new_k] = new_parcels[i]
			end
		end
	end
end


