
using CIFTI
using NearestNeighbors
using HDF5
using OrderedCollections

const sphere = h5read(sphere_file, "sphereLR")
const trees = Dict(
	CORTEX_LEFT => KDTree(sphere[1:nverts_L_full, :]'; leafsize = 10),
	CORTEX_RIGHT => KDTree(sphere[(nverts_L_full + 1):end, :]'; leafsize = 10),
)

# take specified 3d rotational params, generate three 3x3 rotational matrices
# (so result is a 3x3x3 array)
function compute_rotation_mats(xrot::Number, yrot::Number, zrot::Number)
	rot = zeros(Float64, 3, 3, 3)
	rot[1, :, :] = [1 0 0; 0 cos(xrot) -sin(xrot); 0 sin(xrot) cos(xrot)]
	rot[2, :, :] = [cos(yrot) 0 sin(yrot); 0 1 0; -sin(yrot) 0 cos(yrot)]
	rot[3, :, :] = [cos(zrot) -sin(zrot) 0; sin(zrot) cos(zrot) 0; 0 0 1]
	return rot
end

# compute rotations from a pre-defined file of rotational params
function make_rotations(filename::String)
	return @chain begin
		[h5read(filename, x) for x in ["xrot", "yrot", "zrot"]]
		vcat(_...)
		transpose
		[compute_rotation_mats(x...) for x in eachrow(_)]
	end
end

# compute rotations from a randomly initialized set of rotational params
function make_rotations(nrot::Int)
	return @chain begin
		randn(nrot, 3) * 2π
		[compute_rotation_mats(x...) for x in eachrow(_)]
	end
end

# get rotated spherical coords from a trio of rot matrices generated from the above fn
function rotate_on_sphere(rotmat::Array, sphere_coords::Matrix)
	xrot_coords = rotmat[1, :, :] * sphere_coords'
	xyrot_coords = rotmat[2, :, :] * xrot_coords
	return rotmat[3, :, :] * xyrot_coords
end

# take a set of rotated spherical coords from the above fn, and map them to 
# the set of sphere-projected cortical vertices via nearest neighbor calcuation
function get_rotated_parcel(xyzrot_coords::Matrix, tree::KDTree, hem::BrainStructure)
	rotated_inds, dists = knn(tree, xyzrot_coords, 1)
	rotated_inds = [x[1] for x in rotated_inds] # flatten result to just a vector
	rotated_inds .+= (hem == CORTEX_LEFT ? 0 : nverts_L_full)
	rotated_inds = full2trunc[rotated_inds] # reduce to surface verts only
	rotated_parcel = BitVector(undef, nverts) * false
	# if any verts ended up in the medial wall, we'll throw it out
	if !any(rotated_inds .== 0) # a val of 0 would indicate medial wall
		rotated_parcel[rotated_inds] .= true
	end
	return rotated_parcel
end

function fill_in_gaps!(rotated_parcel::BitVector, neigh::VertexList)
	while true
		temp = [sum(rotated_parcel[x]) for x in neigh]
		# to adjust for the first 11 entries of neigh, per hem, having one less elem:
		temp[[1:11; 29697:29707]] .+= 1 
		add_inds = findall(.!rotated_parcel .&& temp .> 4)
		length(add_inds) > 0 || break
		rotated_parcel[add_inds] .= true
	end
end

function dilate_or_contract_parcel!(
		rotated_parcel::BitVector, desired_size::Int, adjmat::AbstractMatrix; maxiter::Int = 20
	)
	curr_size = sum(rotated_parcel)
	Δ = curr_size - desired_size
	i = 1
	# if rotated parcel is smaller than the real one, grow it:
	while Δ < 0 && i < maxiter 
		# find wherever there's no rot parcel assignment but at least one of the neighbors
		# belongs to the parcel; fill up the rot parcel with as many of these border verts
		# as you can, without making the rot parcel bigger than the real one
		temp = findall(rotated_parcel .&& .!baddata)
		border_verts = setdiff(adjmat[:, temp].rowval, temp)
		border_verts = border_verts[1:min(abs(Δ), length(border_verts))]
		rotated_parcel[border_verts] .= true
		Δ += length(border_verts)
		i += 1
		i == maxiter && println("Warning: maxiter condition reached")
	end
	# if rotated parcel is bigger than the real one, shrink it:
	while Δ > 0 && i < maxiter 
		# find wherever there's a rotated parcel vert but at least one of its neighbors 
		# doesn't belong to the parcel; get rid of as many of these border verts as you 
		# can, without making the rot parcel smaller than the real one
		temp = rotated_parcel .&& .!baddata
		temp2 = findall(temp)
		border_verts = temp2[adjmat[temp, .!temp].rowval]
		border_verts = border_verts[1:min(Δ, length(border_verts))]
		rotated_parcel[border_verts] .= false
		Δ -= length(border_verts)
		i += 1
		i == maxiter && println("Warning: maxiter condition reached")
	end
	temp = findall(rotated_parcel .&& .!baddata)
	border_verts = setdiff(adjmat[:, temp].rowval, temp)
	if length(border_verts) / length(temp) >= 3
		rotated_parcel .= false
	end
end

@kwdef struct Parcel
	id::UInt16
	vertices::Vector
	size::Int = length(vertices)
	hem::BrainStructure = vertices[1] <= nverts_L_trunc ? L : R
	spherical_projection::Matrix = sphere[trunc2full[vertices], :]
end

function make_parcels(x)::Parcellation
	Dict(
		[p => Parcel(id = p, vertices = findall(x .== p)) for p in setdiff(x, 0)]
	)
end

function read_parcels(filename::String)
	@chain begin
		CIFTI.load(filename)[LR][:]
		convert.(UInt16, _)
		make_parcels
	end
end

const Parcellation = Dict{T, Parcel} where T <: Integer

# perform a single rotation for a single parcel; modifies rot_verts in place
function process_rotation!(
		rot_verts::Vector{UInt16}, 
		parcel::Parcel,
		rotmats::Array, 
		tree::KDTree, 
		neigh::VertexList, 
		adjmat::AbstractMatrix
	)
	xyzrot_coords = rotate_on_sphere(rotmats, parcel.spherical_projection)
	rotated_parcel = get_rotated_parcel(xyzrot_coords, tree, parcel.hem)
	sum(rotated_parcel) > 0 || return
	fill_in_gaps!(rotated_parcel, neigh)
	sum(rotated_parcel .&& baddata) < minsize || return
	dilate_or_contract_parcel!(rotated_parcel, parcel.size, adjmat)
	rot_verts[rotated_parcel] .= parcel.id
end

function rotation_wrapper(parcels::Parcellation, rotations::Vector{Array{Float64, 3}})
	parcel_ids = collect(keys(parcels))
	result = Vector{Parcellation}(undef, nrot)
	Threads.@threads for r in 1:nrot
		rotmats = rotations[r]
		rot_verts = zeros(UInt16, nverts)
		for id in parcel_ids
			tree = trees[parcels[id].hem]
			process_rotation!(rot_verts, parcels[id], rotmats, tree, neigh, adjmat)
		end
		result[r] = make_parcels(rot_verts)
	end
	return result
end


