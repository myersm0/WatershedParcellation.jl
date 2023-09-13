
using NearestNeighbors

# take specified 3d rotational params, generate three 3x3 rotational matrices
# (so result is a 3x3x3 array)
function make_rotmats(xrot::Number, yrot::Number, zrot::Number)
	rotmat = zeros(Float64, 3, 3, 3)
	rotmat[1, :, :] = [1 0 0; 0 cos(xrot) -sin(xrot); 0 sin(xrot) cos(xrot)]
	rotmat[2, :, :] = [cos(yrot) 0 sin(yrot); 0 1 0; -sin(yrot) 0 cos(yrot)]
	rotmat[3, :, :] = [cos(zrot) -sin(zrot) 0; sin(zrot) cos(zrot) 0; 0 0 1]
	return rotmat
end

# take a trio of rotation matrices generated from the above fn, map it to the sphere
function rotate_on_sphere(rotmat::Array{Number, 3}, sphere_coords::Matrix{Number})
	xrot_coords = rotmat[1, :, :] * sphere_coords'
	xyrot_coords = rotmat[2, :, :] * xrot_coords
	rotmat[3, :, :] * xyrot_coords
end

# take a set of rotated spherical coords from the above fn, and map them to 
# the set of sphere-projected cortical vertices via nearest neighbor calcuation
function get_rotated_parcel(xyzrot_coords::Matrix{Number}, tree::KDTree, hem::Char)
	rotated_inds, dists = knn(tree, xyzrot_coords, 1)
	rotated_inds = [x[1] for x in rotated_inds] # simplify to just a vector
	rotated_inds .+= (hem == 'L' ? 0 : nverts_L)
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
	return rotated_parcel
end

function dilate_or_contract_parcel!(rotated_parcel::BitVector, desired_size::Int64)
	curr_size = sum(rotated_parcel)
	delta = curr_size - desired_size
	maxiter = 20
	i = 1
	# if rotated parcel is smaller than the real one, grow it
	while delta < 0 && i < maxiter 
		# find wherever there's no rot parcel assignment, 
		# but at least one of the neighbors belongs to the parcel;
		# fill up the rot parcel with as many of these border verts as you can, 
		# without making the rot parcel bigger than the real one
		temp = findall(rotated_parcel .&& .!baddata)
		border_verts = setdiff(adjmat[:, temp].rowval, temp)
		border_verts = border_verts[1:min(abs(delta), length(border_verts))]
		rotated_parcel[border_verts] .= true
		delta += length(border_verts)
		i += 1
		if i == maxiter
			println("Warning: maxiter condition reached")
		end
	end
	# if rotated parcel is bigger than the real one, shrink it
	while delta > 0 && i < maxiter 
		# find wherever there's a rotated parcel vert 
		# but at least one of its neighbors doesn't belong to the parcel;
		# get rid of as many of these border verts as you can, 
		# without making the rot parc smaller than the real one
		temp = rotated_parcel .&& .!baddata
		temp2 = findall(temp)
		border_verts = temp2[adjmat[temp, .!temp].rowval]
		border_verts = border_verts[1:min(delta, length(border_verts))]
		rotated_parcel[border_verts] .= false
		delta -= length(border_verts)
		i += 1
		if i == maxiter
			println("Warning: maxiter condition reached")
		end
	end
	temp = findall(rotated_parcel .&& .!baddata)
	border_verts = setdiff(adjmat[:, temp].rowval, temp)
	if length(border_verts) / length(temp) >= 3
		rotated_parcel .= false
	end
end
    



