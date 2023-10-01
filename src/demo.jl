
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
pθ = Parcel(coordinates(hem)[vertices(p), :], tree)
pxθ = rotation_wrapper(parcels, rotations) # 14 min on 6 cores

fig = Figure()
ax = Axis3(fig[1, 1])
scatter!(ax, coordinates(midthickness[L]), color = RGBA(0, 0, 0, 0.2), markersize = 1)
scatter!(ax, coordinates(hem), color = RGBA(0, 0, 0, 0.2), markersize = 1)

scatter!(ax, coordinates(midthickness[L])[vertices(p), :], color = "blue")
scatter!(ax, coordinates(surf)[vertices(p), :], color = "blue")
scatter!(ax, coordinates(surf)[vertices(θp), :], color = "purple")

close!(θp, neigh)
scatter!(ax, coordinates(surf)[vertices(θp), :], color = "red")


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






fig = Figure()
ax = Axis3(fig[1, 1])

p1 = deepcopy(px[4984])
p2 = deepcopy(px[4984])

scatter!(ax, coordinates(hem)[vertices(p), :], color = "black")

dilate!(p2, A)
temp = setdiff(p, p_orig)
scatter!(ax, coordinates(hem)[temp, :], color = "red")

erode!(p2, neigh)
p2.membership == p1.membership

erode!(p2, neigh)
scatter!(ax, coordinates(hem)[vertices(p2), :], color = "purple")

# make some holes in the parcel
exclude_verts = [6702, 6866, 6974]
p2[exclude_verts] .= false
# equivalent to `setdiff!(p2, exclude_verts)`

fig = Figure()
ax = Axis3(fig[1, 1])
scatter!(ax, coordinates(hem)[vertices(p2), :], color = "black")

p3 = deepcopy(p2)
close!(p3, neigh)
setdiff!(p3, p2)
scatter!(ax, coordinates(hem)[vertices(p3), :]', color = "blue")

union!(p3, p2)
resize!(p3, 400; A = A, neigh = neigh)
scatter!(ax, coordinates(hem)[vertices(p3), :], color = "purple")

resize!(p3, size(p); A = A, neigh = neigh)
scatter!(ax, coordinates(hem)[vertices(p3), :], color = "red")

complement(p, p3)
complement(p3, p)

using BenchmarkTools
using StatsBase: sample
test_ids = sample(collect(keys(px)), 10)
@benchmark [size(px[p]) for p in test_ids]

@kwdef struct mstest
	size::Int = rand(1:100, 1)[1]
end

test_structs = [mstest() for _ in 1:10]

@benchmark [x.size for x in test_structs]

test = falses(59412)
test[[9, 99, 999, 9999]] .= true
@benchmark findall(test)
@benchmark sum(test)

test = spzeros(Bool, 59412)
test[[9, 99, 999, 9999]] .= true
@benchmark vertices(px[5051])
@benchmark size(px[5051])

@benchmark vec(px)

# getting vertices or non-zero elements:
# - SparseVector:    150 ns
# - BitVector:       368 ns
# - Vector{T}:         1 ns

# getting intersection of two parcels:
# - SparseVector:    637 ns
# - BitVector:       753 ns
# - Vector{T}:      5083 ns

# adding or removing ~300 elements:
# - SparseVector:   3047 ns
# - BitVector:        85 ns (with union!(a, b))
# - Vector{T}:      7692 ns

# getting size of parcel:
# - SparseVector:     83 ns
# - BitVector:       104 ns
# - Vector{T}:         9 ns

# checking unassigned or nnz:
# - SparseVector:  39000 ns
# - BitVector:     22084 ns
# - BitMatrix:   4517000 ns

# checking overlap of two parcels:
# - SparseVector:    812 ns
# - BitVector:       108 ns
# - Vector{T}:     49110 ns

a = px[6779]
b = px[5051]
@benchmark intersect(a, b).nzval

a = vertices(px[6779])
b = vertices(px[5051])
@benchmark intersect(a, b)

a = falses(length(px))
a[vertices(px[6779])] .= true

b = falses(length(px))
b[vertices(px[5051])] .= true
@benchmark a .& b

temp = vertices(px[1722])
a = deepcopy(px[6779])
b = falses(length(px))
b[vertices(px[6779])] .= true
c = Parcel(temp; n = length(px))

@benchmark union!(a, c)
a = deepcopy(px[6779])
@benchmark a[temp] .= true
@benchmark b[temp] .= true

@benchmark a[temp] .= false
@benchmark b[temp] .= false

c = vertices(a)
d = vertices(b)
@benchmark union(c, temp)

@benchmark setdiff(c, temp)

c = spzeros(Bool, 32492)
d = spzeros(Bool, 32492)
c[vertices(a)] .= true
d[vertices(b)] .= true

a = vertices(a)

temp = vertices(px[1722])
a = deepcopy(px[6779])
b = falses(length(px))
b[vertices(px[6779])] .= true

@benchmark a[temp] .= true
@benchmark b[temp] .= true

@benchmark size(a)
@benchmark sum(b)
@benchmark length(c)


@benchmark vertices(a)
@benchmark findall(b)
@benchmark c

@benchmark unassigned(px[L])
@benchmark nnz(px[L])




