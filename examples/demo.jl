
using CIFTI
using CorticalSurfaces
using CorticalParcels
using WatershedParcellation

using JLD
using NearestNeighbors
using JSON

config = open("config.json", "r") do fid
	JSON.parse(fid)
end

include("surface_setup.jl")

template = config["template"]


# ===== make edgemap =====
grads = CIFTI.load(config["outputs"]["smoothed gradients"])
minima = find_minima(grads.data, c)
edgemap = run_watershed(grads[LR], minima, c; nsteps = config["nsteps"])
outname = config["outputs"]["edgemap"]
CIFTI.save(outname, Float32.(edgemap); template = template)

# gradients no longer needed; free memory
grads = nothing
GC.gc()


# ===== make parcels =====
edgemap = CIFTI.load(config["outputs"]["edgemap"])
labels = run_watershed(edgemap[LR][:], c; thresh_quantile = 0.75)

hem = L
verts = @collapse vertices(c[hem])
px = Parcellation{Int}(c[hem], labels[verts])
edges = pad(edgemap[hem][:], c[hem])


# ==== cleanup (enforce maximum height values, etc) =====
remove_weak_boundaries!(px, edges; threshold = 0.3, radius = 30)
merge_small_parcels!(px, edges; minsize = 30, radius = 30)
cap_at_height!(px, edges; threshold = 0.9)
merge_small_parcels!(px, edges; minsize = 30, radius = 30)
cap_at_height!(px, edges; threshold = 0.9)
remove_articulation_points!(px; minsize = 4)
remove_small_parcels!(px; minsize = 10)

temp = Float32[trim(vec(px), c[hem]); zeros(29716)]
CIFTI.save(config["outputs"]["parcels"], temp; template = template)


# ===== make rotations =====
tree = KDTree(coordinates(c[hem]); leafsize = 10)
nrot = config["nrotations"]
rotations = make_rotations(nrot)
pxθ = rotation_wrapper(px, rotations, tree)


# ===== homogeneity testing =====
dconn = CIFTI.load(config["dconn"])
GC.gc()
cov_corr = make_cov_corr(dconn[L, L], c[hem])
GC.gc()

# test a parcel from the parcellation
p = px[37]
homogeneity_test(p, cov_corr; criteria = x -> default_criteria(x))

# test a parcel from one of the rotated parcellations
pθ = pxθ[1][1479]
homogeneity_test(pθ, cov_corr; criteria = x -> default_criteria(x))

# run a homogeneity test on the whole real parcellation
homogeneity_test(px, cov_corr; criteria = x -> default_criteria(x))

# load in the "baddata" mask (map of low-signal regions):
temp = CIFTI.load(config["low-signal map"])[L][:]
baddata = Parcel(c[hem])
baddata[findall(pad(temp, c[hem]) .== 1)] .= true

# define custom criteria function involving low-signal region exclusion:
function criteria(p::Parcel, baddata::Parcel)
	overlap(p, medial_wall(p.surface)) == 0 || return false
	complement(p, baddata) >= 15 || return false
	return true
end

homogeneity_test(px, cov_corr; criteria = x -> criteria(x, baddata)

# now save the resulting homogeneity to a NamedVector (indexed by parcel keys)
ks = collect(keys(px))
real_result = homogeneity_test(px, cov_corr; criteria = x -> criteria(x, baddata))

# do the same for the vector of 1000 rotated parcellations 
# (also to be indexed by the same set of keys as in the above)
rot_result = homogeneity_test(pxθ, cov_corr; criteria = x -> criteria(x, baddata), ks = ks)

# get the mean homogeneity of parcels from the real parcellation
summarize_homogeneity(real_result)

# get the mean homogeneity and std-dev of parcels from the rotated parcellations
summarize_homogeneity(rot_result)

# calculate a z-score of the real parcellation versus the rotated ones
summarize_homogeneity(real_result, rot_result)

