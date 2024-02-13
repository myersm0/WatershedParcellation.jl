
using CIFTI
using CorticalSurfaces
using CorticalParcels
using WatershedParcellation

using JLD
using NearestNeighbors
using JSON

# you could manually set up a config.json file yourself, like this example;
# the repo provides some of the needed params and input materials, but
# unfortunately you'll need some things to exist elsewhere, in your own filesystem,
# like a dconn, because they're too big for me to put on github
config = open("config.json", "r") do fid
	JSON.parse(fid)
end

# this loads in a surface definition that's been prepared in advance from a GIFTI file
# into a .jld data serialization. The main thing it creates is the
# CorticalSurface struct we'll call c, which contains two Hemisphere structs L and R
# that will be needed for spatial- and adjacency-related info throughout this demo
include("surface_setup.jl")

# the path of a CIFTI file in the same space as the parcellation to be created;
# probably a dscalar or a dtseries with one timepoint and no subcortical structures
template = config["template"]


# ===== make edgemap =====

# load in the gradients that you've computed and saved previously, 
# as demonstrated in make_gradients.jl
grads = CIFTI.load(config["outputs"]["smoothed gradients"])

# generate a BitMatrix of local minima which will serve as initialization points
# for the watershed algrotihm
minima = find_minima(grads[LR], c)

# now generate an edgemap by iterating through the gradient map one vertex at a time
edgemap = run_watershed(grads[LR], minima, c; nsteps = config["nsteps"])

# save out the results
outname = config["outputs"]["edgemap"]
CIFTI.save(outname, Float32.(edgemap); template = template)

# gradients no longer needed; free memory
grads = nothing
GC.gc()


# ===== make parcels =====
edgemap = CIFTI.load(config["outputs"]["edgemap"])

# this will be our initial parcellation, to be tuned in the next section
labels = run_watershed(edgemap[LR][:], c; thresh_quantile = 0.75)

# For now, we'll just work with a single hemisphere (L) from the above output.
# You could run the results for each hem separately, then concatenate
# them later for homogeneity testing. Eventually I plan to add support
# for doing a full bilateral parcellation all at once, but for now it's more
# efficient and probably conceptually better (though a little more complicated)
# to do it this way.
hem = L
verts = @collapse vertices(c[hem])
px = HemisphericParcellation{Int}(c[hem], labels[verts])
edges = pad(edgemap[hem][:], c[hem])


# ==== cleanup (enforce maximum height values, etc) =====
remove_weak_boundaries!(px, edges; threshold = 0.3, radius = 30)
merge_small_parcels!(px, edges; minsize = 30, radius = 30)
threshold!(px, edges; threshold = 0.9)
merge_small_parcels!(px, edges; minsize = 30, radius = 30)
threshold!(px, edges; threshold = 0.9)
remove_articulation_points!(px; minsize = 4)
remove_small_parcels!(px; minsize = 10)

temp = Float32[trim(vec(px), c[hem]); zeros(29716)]
CIFTI.save(config["outputs"]["parcels"], temp; template = template)


# ===== make rotations =====
# make a KDTree that will assist in nearest neighbors search in the rotation process
tree = KDTree(coordinates(c[hem]); leafsize = 10)

# generate some random 3x3x3 arrays that will be used to rotate the parcels
rotation_matrices = make_rotations(config["nrotations"])

# now with the help of those two items, create a vector of rotated parcellations
# which will be compared against the real parcellation for homogeneity testing below
pxθ = rotation_wrapper(px, rotation_matrices, tree)


# ===== homogeneity testing =====
dconn = CIFTI.load(config["dconn"])
GC.gc()
cov_corr = make_cov_corr(dconn[L, L], c[hem])
GC.gc()

# run a homogeneity test on a parcel from the parcellation; this simply:
# checks whether the parcel meets inclusion criteria for analysis
# (in this case, the provided function `default_criteria()` is used which 
# simply checks that the parcel doesn't overlap with the medial wall, because
# it's expected that there's no connectivity information defined for the 
# medial wall), and if so returns the homogeneity of the matrix cov_corr when
# subset by the vertices of p, or NaN otherwise
p = px[37]
homogeneity_test(p, cov_corr; criteria = p -> default_criteria(p))

# test the same parcel from one of the rotated parcellations; it will 
# probably have lower homogeneity
pθ = pxθ[1][37]
homogeneity_test(pθ, cov_corr; criteria = p -> default_criteria(p))

# run a homogeneity test on the whole real parcellation
homogeneity_test(px, cov_corr; criteria = p -> default_criteria(p))

# load in the "baddata" mask (map of low-signal regions)
temp = CIFTI.load(config["low-signal map"])[L][:]
baddata = Parcel(c[hem])
baddata[findall(pad(temp, c[hem]) .== 1)] .= true

# define a custom criteria function that will reject parcels overlapping with
# the medial wall (the same as in default_criteria() that we used before),
# and also requiring the parcel to have at least 15 vertices outside of the
# low-signal mask regions
function criteria(p::Parcel, baddata::Parcel)
	overlap(p, medial_wall(p.surface)) == 0 || return false
	complement(p, baddata) >= 15 || return false
	return true
end

# run a homogeneity test on the whole real parcellation again, but this time
# using our custom criteria function that will reject parcels on the basis of
# overlap with the "baddata" mask
homogeneity_test(px, cov_corr; criteria = p -> criteria(p, baddata)

# now save the resulting homogeneity to a NamedVector (indexed by parcel keys)
ks = collect(keys(px))
real_result = homogeneity_test(px, cov_corr; criteria = p -> criteria(p, baddata))

# do the same for the vector of 1000 rotated parcellations (also to be indexed by 
# the same set of keys as in the above, so that the results are directly comparable)
rot_result = homogeneity_test(pxθ, cov_corr; criteria = p -> criteria(p, baddata), ks = ks)

# get the mean homogeneity of parcels from the real parcellation
summarize_homogeneity(real_result)

# get the mean homogeneity and std-dev of parcels from the rotated parcellations
summarize_homogeneity(rot_result)

# calculate a z-score of the real parcellation versus the rotated ones
zscore = summarize_homogeneity(real_result, rot_result)

