
using Statistics: mean
using Chain
using StatsBase: sample
using ThreadsX

export run_watershed

# doesn't return a val, but modifies label and watershed_zones vectors in place
function eval_at_height!(h, label, edgemetric, watershed_zones, neighbors)
	nodes_at_threshold = @chain begin
		(edgemetric .< h) .& (label .== 0) .& (watershed_zones .== 0)
		findall
		sample(_, length(_); replace = false)
	end
	for node in nodes_at_threshold
		nodeneighlab = setdiff(label[neighbors[node]], 0)
		length(nodeneighlab) > 0 || continue
		minnodeneighlab = minimum(nodeneighlab)
		maxnodeneighlab = maximum(nodeneighlab)
		if minnodeneighlab != maxnodeneighlab
			watershed_zones[node] = true
			label[node] = 0
		else
			label[node] = minnodeneighlab
		end
	end
end

# edges arg will be modified in place
function watershed_chunk!(
		edges::AbstractMatrix, 
		metric::AbstractMatrix, 
		minima::BitMatrix, 
		neigh::Vector{Vector{Int}},
		chunk::UnitRange, 
		heights::StepRangeLen
	)
	for i in chunk
		nverts = size(minima, 1)
		label = @view edges[:, i] 
		edgemetric = @view metric[:, i] 
		labelpos = findall(minima[:, i])
		randval = randn(length(labelpos))
		labelnums = sortperm(randval)
		label[labelpos] = labelnums
		watershed_zones = zeros(nverts)
		[eval_at_height!(h, label, edgemetric, watershed_zones, neigh) for h in heights]
	end
end

"""
    run_watershed(metric, minima, neighbors)

Run the watershed algorithm on `metric`, using its `minima` (as returned from 
`find_minima()`) as initialization points and guided by the topology given in
the adjacency list `neighbors`.

Some optional parameters can be tuned:
- `nsteps2`: the number of bins into which to discretize the heights (default `400`)
- `fracmaxh`: a scaling factor to determine a ceiling on height values (default `1.0`)
- `nchunks`: split up the work into this many chunks for multithreading (default `64`)
"""
function run_watershed(
		metric::Matrix, minima::BitMatrix, neigh::Vector{Vector{Int}};
		nsteps2::Int = 400, fracmaxh::Float64 = 1.0, nchunks::Int = 64
	)
	@assert all([nsteps2, nchunks, fracmaxh] .> 0)
	minheight = minimum(grads)
	maxheight = maximum(grads) * fracmaxh
	heights = range(minheight, maxheight, length = nsteps2)
	nverts = length(neigh)
	chunk_size = Int(floor(nverts / nchunks))
	edges = zeros(UInt16, nverts, nverts)
	chunks = [((c - 1) * chunk_size + 1):min(nverts, c * chunk_size) for c in 1:nchunks]
	ThreadsX.foreach(
		chunk -> @time watershed_chunk!(edges, metric, minima, neigh, chunk, heights), 
		chunks
	)
	return mean(edges .== 0; dims = 2)[:]
end

