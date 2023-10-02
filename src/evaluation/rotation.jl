
"Given the specified 3D rotational parameters, generate three 3x3 rotational matrices"
function compute_rotation_mats(x::T, y::T, z::T)::Array{T, 3} where T <: Real
	rot = zeros(Float64, 3, 3, 3)
	rot[1, :, :] = [1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)]
	rot[2, :, :] = [cos(y) 0 sin(y); 0 1 0; -sin(y) 0 cos(y)]
	rot[3, :, :] = [cos(z) -sin(z) 0; sin(z) cos(z) 0; 0 0 1]
	return rot
end

"Compute rotations from a pre-defined file of rotational parameters"
function make_rotations(filename::String)
	return @chain begin
		[h5read(filename, x) for x in ["xrot", "yrot", "zrot"]]
		vcat(_...)
		transpose
		[compute_rotation_mats(x...) for x in eachrow(_)]
	end
end

"Compute rotations from a randomly initialized set of rotational parameters"
function make_rotations(nrot::Int)
	return @chain begin
		randn(nrot, 3) * 2π
		[compute_rotation_mats(x...) for x in eachrow(_)]
	end
end

"Get rotated spherical coordinates from a set of three rotation matrices"
function rotate_on_sphere(rotmat::Array, sphere_coords::Matrix)
	xrot_coords = rotmat[1, :, :] * sphere_coords'
	xyrot_coords = rotmat[2, :, :] * xrot_coords
	return rotmat[3, :, :] * xyrot_coords
end

# perform a single rotation for a single parcel
function process_rotation(
		px::Parcellation{T},
		id::T,
		surf::Hemisphere,
		rotmats::Array, 
		tree::KDTree, 
		A::AbstractMatrix,
		neigh::Vector{Vector{Int}}
	) where T
	θcoords = rotate_on_sphere(rotmats, coordinates(surf)[vertices(px[id]), :])
	θp = Parcel(θcoords, tree)
	size(θp) > 0 || return
	close!(θp, neigh)
	overlap(θp, baddata) < 15 || return Parcel(size(surf))
	resize!(θp, size(px[id]); A = A, neighbors = neigh)
	return θp
end

function rotation_wrapper(
		px::Parcellation{T}, rotations::Vector{Array{Float64, 3}}
	) where T
	parcel_ids = collect(keys(px))
	nrot = size(rotations, 1)
	surf = px.surface
	result = Vector{Parcellation{T}}(undef, nrot)
	Threads.@threads for r in 1:nrot
		rotmats = rotations[r]
		temp = Dict{T, Parcel}()
		for id in parcel_ids
			temp[id] = process_rotation(px, id, surf, rotmats, tree, A, neigh)
		end
		result[r] = Parcellation(surf, temp)
	end
	return result
end


