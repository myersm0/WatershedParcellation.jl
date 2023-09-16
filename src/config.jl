
using JLD
using CIFTI
using Chain

const VertexList = Vector{Vector{UInt16}}
const assets_dir = joinpath(@__DIR__, "../data/32k_tools/")

# params that should be configurable but hard-coding for now
const mw_verts = load("$assets_dir/mw_verts.jld", "mw_verts")
const refcifti = CIFTI.load("$assets_dir/elabe_baddata.dtseries.nii")
const baddata = @chain refcifti[LR] vec convert(BitVector, _)
const nverts = length(baddata)
const nverts_L_trunc = length(refcifti.brainstructure[L])
const nverts_L_full = nverts_L_trunc + sum(mw_verts .<= nverts_L_trunc)
const full2trunc = zeros(Int, nverts_L_full * 2)
const trunc2full = setdiff(1:(nverts + length(mw_verts)), mw_verts)
full2trunc[setdiff(1:length(full2trunc), mw_verts)] .= 1:length(baddata)
const minsize = 15
const nrot = 1000

parcellation_name = "test"

# hard-coded paths; later give people the option to specify themselves
sphere_file = "$assets_dir/sphereLR.h5"
rotations_file = "$assets_dir/../rotations.h5"
parcel_file = "$assets_dir/../test_parcels.dtseries.nii"



