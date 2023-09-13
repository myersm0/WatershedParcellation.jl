
using Statistics: mean
using Chain
using StatsBase: sample
using ThreadsX

const nverts = 59412

# doesn't return a val, but modifies label and watershed_zones vectors in place
function eval_at_height(h, label, edgemetric, watershed_zones, neighbors)
	nodes_at_threshold = @chain begin
		(edgemetric .< h) .&& (label .== 0) .&& (watershed_zones .== 0)
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

function watershed_chunk(
		edges::AbstractMatrix, 
		grads::AbstractMatrix, 
		minima::BitMatrix, 
		neigh::VertexList,
		chunk::UnitRange, 
		heights::StepRangeLen
	)
	for i in chunk
		label = @view edges[:, i] 
		edgemetric = @view grads[:, i] 
		labelpos = findall(minima[:, i])
		randval = randn(length(labelpos))
		labelnums = sortperm(randval)
		label[labelpos] = labelnums
		watershed_zones = zeros(nverts)
		[eval_at_height(h, label, edgemetric, watershed_zones, neigh) for h in heights]
	end
end

function run_watershed(
		grads::Matrix, minima::BitMatrix, neigh::VertexList;
		nsteps::Int = 400, fracmaxh::Float64 = 1.0, nchunks::Int = 64
	)
	minheight = minimum(grads)
	maxheight = maximum(grads) * fracmaxh
	heights = range(minheight, maxheight, length = nsteps)
	chunk_size = Int(floor(nverts / nchunks))
	edges = zeros(UInt16, nverts, nverts)
	chunks = [((c - 1) * chunk_size + 1):min(nverts, c * chunk_size) for c in 1:nchunks]
	ThreadsX.foreach(
		chunk -> watershed_chunk(edges, grads, minima, neigh, chunk, heights), 
		chunks
	)
	return mean(edges .== 0; dims = 2)[:]
end

