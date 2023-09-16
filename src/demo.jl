
# edgemap creation
grads = load_gradients(filename)
neigh = load_neighbors()
adjmat = make_adjmat(neigh)
minima = find_minima(grads, adjmat)
edgemap = run_watershed(grads, minima, neigh)

# TODO: I need to exclude baddata real parcels before evaluating

# rotations
include("evaluation/rotation.jl")
rotations = make_rotations(rotations_file)
parcels = read_parcels(parcel_file)
rotated_parcels = rotation_wrapper(parcels, rotations) # 14 min on 6 cores

#make_rotmats(x, y, z) -> compute_rotation_mats(...)
#make_rotmats(filename) -> make_rotations(...)






