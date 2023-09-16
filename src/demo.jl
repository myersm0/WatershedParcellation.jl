
# edgemap creation
grads = load_gradients(filename)
neigh = load_neighbors()
adjmat = make_adjmat(neigh)
minima = find_minima(grads, adjmat)
edgemap = run_watershed(grads, minima, neigh)

# TODO: I need to exclude baddata real parcels before evaluating

# rotations
include("evaluation/rotation.jl")
rotations = make_rotmats(rotations_file)
test = rotation_wrapper(parcel_file, rotations) # 14 min on 6 cores




