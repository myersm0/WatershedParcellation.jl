
using LinearAlgebra
using Statistics
using NamedArrays

function calc_homog(mat::Matrix)
	eigenvals = eigvals(mat)
	return maximum(eigenvals) / sum(eigenvals)
end

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

# 0.005165 - 0.0269 sec
function homogeneity_test(parcel::Parcel, cov_corr::Matrix; minsize::Int = 15)
	verts = parcel.vertices
	return parcel.size < minsize ? NaN : calc_homog(cov_corr[verts, verts])
end

function homogeneity_test(p::Parcellation, cov_corr::Matrix, ks::Vector; minsize::Int = 15)
	result = NamedArray(zeros(length(ks)) * NaN, (ks,))
	for k in intersect(ks, keys(p))
		result[Name(k)] = homogeneity_test(p[k], cov_corr; minsize = minsize)
	end
	return result
end

function homogeneity_test(
		vp::Vector{Parcellation}, cov_corr::Matrix, ks::Vector; minsize::Int = 15
	)
	ntests = length(vp)
	result = Vector{NamedVector}(undef, ntests)
	Threads.@threads for i in 1:ntests
		result[i] = homogeneity_test(vp[i], cov_corr, ks; minsize = minsize)
	end
	return hcat(result...)
end







