
temp = load(config["32k objects"])

surf = temp["pointsets"]["sphere"][L]
mid = temp["pointsets"]["midthickness"][L]
mw = temp["medial wall"][L]
triangle = temp["triangle"][L]
sphere = Hemisphere(L, surf, mw; triangles = triangle)
hemL = sphere

surf = temp["pointsets"]["sphere"][R]
mid = temp["pointsets"]["midthickness"][R]
mw = temp["medial wall"][R]
triangle = temp["triangle"][R]
sphere = Hemisphere(R, surf, mw; triangles = triangle)
hemR = sphere

c = CorticalSurface(hemL, hemR)
initialize_adjacencies!(c)

dmat = load(config["distance matrix"])

c[L][:distances] = dmat["geodesic"][L]
c[R][:distances] = dmat["geodesic"][R]



