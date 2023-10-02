
using CorticalSurfaces
using CorticalParcels
using JLD
using CIFTI
using OrderedCollections
using Chain
using NamedArrays
using SparseArrays
using NearestNeighbors
using StatsBase
using BenchmarkTools
using GLMakie
using Colors

data_dir = "../../ParcelOps/data/"
conte_file = "$data_dir/conte69.32k_fs_LR.jld"
conte = load(conte_file)
surfL = conte["pointsets"]["sphere"][L]
surfR = conte["pointsets"]["sphere"][R]
midL = conte["pointsets"]["midthickness"][L]
midR = conte["pointsets"]["midthickness"][R]
mwL = conte["medial wall"][L]
mwR = conte["medial wall"][R]
neigh = conte["adjacency list"]
A = make_adjacency_matrix(neigh)
c = CorticalSurface(Hemisphere(surfL, mwL), Hemisphere(surfR, mwR))

midthickness = CorticalSurface(Hemisphere(midL, mwL), Hemisphere(midR, mwR))

baddata_file = "$data_dir/elabe_baddata.dtseries.nii"
temp = CIFTI.load(baddata_file)
baddata = Dict(
	hem => 
		@chain begin
			findall(temp[hem] .!= 0)
			expand(_, c[hem])
			Parcel(_; n = size(c[hem]))
		end
	for hem in LR
)

parcel_file = "$data_dir/test_parcels.dtseries.nii"
temp = CIFTI.load(parcel_file)
labeltype = Int
px = Dict(hem => Parcellation{labeltype}(c[hem], temp[hem]) for hem in LR)

trees = Dict(
	L => KDTree(coordinates(c[L])'; leafsize = 10),
	R => KDTree(coordinates(c[R])'; leafsize = 10),
)

const nrot = 1000
px = px[L]
tree = trees[L]
hem = c[L]

append!(hem, :A, A)
append!(hem, :neighbors, neigh)

include("generation/minima.jl")
include("generation/watershed.jl")
include("evaluation/rotation.jl")
include("evaluation/homogeneity.jl")

# edgemap creation
grads = load_gradients(filename)
minima = find_minima(grads, A)
edgemap = run_watershed(grads, minima, neigh)

# rotations
rotations = make_rotations(nrot)
p = deepcopy(px[11])
θp = Parcel(coordinates(hem)[vertices(p), :], tree)
θpx = rotation_wrapper(parcels, rotations) # 14 min on 6 cores

fig = Figure()
ax = Axis3(fig[1, 1])
scatter!(ax, coordinates(midthickness[L]), color = RGBA(0, 0, 0, 0.2), markersize = 1)
scatter!(ax, coordinates(hem), color = RGBA(0, 0, 0, 0.2), markersize = 1)

scatter!(ax, coordinates(midthickness[L])[vertices(p), :], color = "blue")
scatter!(ax, coordinates(surf)[vertices(p), :], color = "blue")
scatter!(ax, coordinates(surf)[vertices(θp), :], color = "purple")

close!(θp, neigh)
scatter!(ax, coordinates(surf)[vertices(θp), :], color = "red")

dconn = load("$data_dir/5003_pre-massage.dconn_L.jld", "dconn")
cov_corr = make_cov_corr(dconn, hem)
dconn = zeros(Float32, 2, 2)
GC.gc()

include("evaluation/homogeneity.jl")
dts = CIFTI.load("/Users/myersm/5000_all_sessions.dtseries.nii")[L]
dconn = cor(dts')
cov_corr = cov(dconn)
real_homog = [test_parcel(p, cov_corr) for (k, p) in parcels]

# parcel creation
edgemap = CIFTI.load(
	"/Users/myersm/.julia/dev/WatershedParcellation/test/data/5003_post-massage_wateredge_avg.dtseries.nii"
)[LR][:]

grads = load("../../ParcelOps/data/5003_all_sessions_fullgrads_L.jld", "grads")

verts = vertices(hem, Exclusive())
minima = find_minima(grads, A[verts, verts])

neigh2 = [collapse(x, hem) for x in neigh[verts]]

edgemap = run_watershed(grads, minima, neigh2; nsteps = 50) # 2.5 minutes
GC.gc()
edgemap = run_watershed(grads, minima, neigh2; nsteps = 400)




