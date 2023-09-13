
grads = load_gradients(filename)
neigh = load_neighbors()
adjmat = make_adjmat(neigh)
minima = make_minima_metrics(grads, adjmat)
edgemap = run_watershed(grads, minima, neigh)



