
using Distributed
addprocs(8)

@everywhere using Chain
@everywhere using SharedArrays
@everywhere using ParallelDataTransfer
@everywhere using StatsBase: sample

const nverts = 59412

# doesn't return a val, but modifies label and watershed_zones vectors in place
@everywhere function eval_at_height(h, label, edgemetric, watershed_zones, neighbors)
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

@everywhere function watershed_chunk(edges::AbstractMatrix, grads::AbstractMatrix, minima::BitMatrix, neigh::Vector{Vector{UInt16}}, chunk::Int, chunk_size::Int, hiter::StepRangeLen)
	start = (chunk - 1) * chunk_size + 1
	stop = min(nverts, chunk * chunk_size)
	for i in start:stop
		label = @view edges[:, i] 
		edgemetric = @view grads[:, i] 
		labelpos = findall(minima[:, i])
		randval = randn(length(labelpos))
		labelnums = sortperm(randval)
		label[labelpos] = labelnums
		watershed_zones = zeros(nverts)
		[eval_at_height(h, label, edgemetric, watershed_zones, neigh) for h in hiter]
	end
end

sendto(workers(); minima = minima)
sendto(workers(); neigh = neigh)

function run_watershed(
		grads::Matrix, minima::BitMatrix, neigh::Vector{Vector{UInt16}};
		nsteps::Int = 400, fracmaxh::Float64 = 1.0, nchunks::Int = 64
	)
	minheight = minimum(grads)
	maxheight = maximum(grads) * fracmaxh
	heights = range(minheight, maxheight, length = nsteps)
	chunk_size = Int(floor(nverts / nchunks))
	edges = SharedArray(zeros(UInt16, nverts, nverts))
	pmap(
		chunk -> watershed_chunk(edges, grads, minima, neigh, chunk, chunk_size, heights),
		1:nchunks
	)
	return mean(edges .== 0; dims = 2)[:]
end









