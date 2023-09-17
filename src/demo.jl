
# edgemap creation
grads = load_gradients(filename)
neigh = load_neighbors()
adjmat = make_adjmat(neigh)
minima = find_minima(grads, adjmat)
edgemap = run_watershed(grads, minima, neigh)


# rotations
# TODO: I need to exclude baddata real parcels before evaluating
include("evaluation/rotation.jl")
rotations = make_rotations(rotations_file)
parcels = read_parcels(parcel_file)
rotated_parcels = rotation_wrapper(parcels, rotations) # 14 min on 6 cores

#make_rotmats(x, y, z) -> compute_rotation_mats(...)
#make_rotmats(filename) -> make_rotations(...)

include("evaluation/homogeneity.jl")
dts = CIFTI.load("/Users/myersm/5000_all_sessions.dtseries.nii")[L]
dconn = cor(dts')
cov_corr = cov(dconn)
real_homog = [test_parcel(p, cov_corr) for (k, p) in parcels]


# parcel creation
edgemap = CIFTI.load(
	"/Users/myersm/.julia/dev/WatershedParcellation/test/data/5003_post-massage_wateredge_avg.dtseries.nii"
)[LR][:]










