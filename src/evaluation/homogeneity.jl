
using LinearAlgebra
using Statistics
using NamedArrays

function make_cov_corr(dconn::Matrix, nverts_L::Int)
	nelem = size(dconn, 1)
	covL = cov(dconn[:, 1:nverts_L])
	covR = cov(dconn[:, (nverts_L + 1):end])
	sizeL = size(covL, 1)
	sizeR = size(covR, 1)
	cov_corr = zeros(nelem, nelem)
	cov_corr[1:sizeL, 1:sizeL] .= covL
	cov_corr[(sizeL + 1):end, (sizeL + 1):end] .= covR
	return cov_corr
end

"""
    homogeneity(mat)

Compute the homogeneity of a `Matrix` as the percent variance explained by its
first prinipal component in a PCA
"""
function homogeneity(mat::Matrix)
	eigenvals = eigvals(mat)
	return maximum(eigenvals) / sum(eigenvals)
end

"""
	 homogeneity_test(p, cov_corr; minsize = 15)

Test homogeneity of a `Parcel` with respect to a covariance of correlations matrix,
or return `NaN` if the parcel is smaller than `minsize`
"""
function homogeneity_test(p::Parcel, cov_corr::Matrix; minsize::Int = 15)
	verts = vertices(p)
	return length(verts) < minsize ? NaN : homogeneity(cov_corr[verts, verts])
end

"""
	 homogeneity_test(px, cov_corr; minsize = 15)

Iterate over `Parcel`s within a `Parcellation`, testing the homogeneity of each.

Optionally supply a set of parcel keys `ks`, which are not necessarily those of `px`!
(The idea is to allow the possibility of using a fixed set of keys across many calls
to this function, and the `Parcellation`s passed each time may not exactly share
all of them. See the method below that accepts a `Vector{Parcellation{T}}`.)
"""
function homogeneity_test(
		px::Parcellation{T}, cov_corr::Matrix; 
		minsize::Int = 15, ks::Vector{T} = collect(keys(px))
	) where T
	result = NamedArray(zeros(length(ks)) * NaN, (ks,))
	for k in intersect(ks, keys(p))
		result[Name(k)] = homogeneity_test(px[k], cov_corr; minsize = minsize)
	end
	return result
end

"""
	 homogeneity_test(px, cov_corr; minsize = 15)

Iterate over a `Vector` of `Parcellation`s, testing the homogeneity of each. Returns
a `NamedMatrix{T}` where parcel keys are given along the rows.
"""
function homogeneity_test(
		vp::Vector{Parcellation{T}}, cov_corr::Matrix; minsize::Int = 15, ks::Vector{T}
	) where T
	ntests = length(vp)
	result = Vector{NamedVector}(undef, ntests)
	Threads.@threads for i in 1:ntests
		result[i] = homogeneity_test(vp[i], cov_corr; minsize = minsize, ks = ks)
	end
	return hcat(result...)
end


