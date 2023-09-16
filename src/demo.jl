
# edgemap creation
grads = load_gradients(filename)
neigh = load_neighbors()
adjmat = make_adjmat(neigh)
minima = find_minima(grads, adjmat)
edgemap = run_watershed(grads, minima, neigh)

# rotations
rotations = make_rotmats(rotations_file)




