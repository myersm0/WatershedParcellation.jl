
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

