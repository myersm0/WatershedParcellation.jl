
temp = load(config["32k objects"])

surf = temp["pointsets"]["sphere"][L]
mid = temp["pointsets"]["midthickness"][L]
mw = temp["medial wall"][L]
triangle = temp["triangle"][L]
sphere = Hemisphere(surf, mw; triangles = triangle)
hemL = sphere
hemL[:neighbors] = temp["adjacency list"]
hemL[:A] = make_adjacency_matrix(sphere)

surf = temp["pointsets"]["sphere"][R]
mid = temp["pointsets"]["midthickness"][R]
mw = temp["medial wall"][R]
triangle = temp["triangle"][R]
sphere = Hemisphere(surf, mw; triangles = triangle)
hemR = sphere
hemR[:neighbors] = temp["adjacency list"]
hemR[:A] = make_adjacency_matrix(sphere)

c = CorticalSurface(hemL, hemR)

dmat = load(config["distance matrix"])

c[L][:distances] = dmat["geodesic"][L]
c[R][:distances] = dmat["geodesic"][R]



