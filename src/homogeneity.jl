
using StatsBase: std

export make_cov_corr, homogeneity_test, default_criteria, summarize_homogeneity

"""
	 make_cov_corr(mat, hem)

Make a covariance matrix with indices aligned with those of a Hemisphere struct `hem`
"""
function make_cov_corr(mat::Matrix, hem::Hemisphere)
	nelem = size(mat, 1)
	if nelem == size(hem, Inclusive())
		inds = vertices(hem)
	else
		inds = collapse(vertices(hem, Exclusive()), hem)
	end
	return cov(mat[:, inds])
end

function default_criteria(p::Parcel)
	return overlap(p, medial_wall(p.surface)) == 0
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
or return `NaN` if the parcel fails to satisfy inclusion `criteria`
"""
function homogeneity_test(
		p::Parcel, cov_corr::Matrix; criteria::Function
	)
	verts = collapse(vertices(p), p.surface)
	return criteria(p) ? homogeneity(cov_corr[verts, verts]) : NaN
end

"""
	 homogeneity_test(px, cov_corr; ks, criteria)

Iterate over `Parcel`s within a `Parcellation`, testing the homogeneity of each.

Optionally supply a set of parcel keys `ks`, which are not necessarily those of `px`!
(The idea is to allow the possibility of using a fixed set of keys across many calls
to this function, and the `Parcellation`s passed each time may not exactly share
all of them. See the method below that accepts a `Vector{Parcellation{T}}`.)
"""
function homogeneity_test(
		px::Parcellation{T}, cov_corr::Matrix; 
		ks::Vector{T} = collect(keys(px)), criteria::Function = x -> default_criteria
	) where T
	result = NamedArray(zeros(length(ks)) * NaN, (ks,))
	for k in intersect(ks, keys(px))
		result[Name(k)] = homogeneity_test(px[k], cov_corr; criteria = criteria)
	end
	return result
end

"""
	 homogeneity_test(px, cov_corr; ks, criteria)

Iterate over a `Vector` of `Parcellation`s, testing the homogeneity of each. Returns
a `NamedMatrix{T}` where parcel keys are given along the rows.
"""
function homogeneity_test(
		vp::Vector{Parcellation{T}}, cov_corr::Matrix;
		ks::Vector{T}, criteria::Function = x -> default_criteria
	) where T
	ntests = length(vp)
	names = (ks, collect(1:ntests))
	dimnames = ("parcel", "iteration")
	result = NamedArray(zeros(length(ks), ntests), names, dimnames)
	for i in 1:ntests
		result[:, i] = homogeneity_test(vp[i], cov_corr; ks = ks, criteria = criteria)
	end
	return result
end

function summarize_homogeneity(m::NamedMatrix)
	mat = rot_result
	nparc, nrot = size(mat)
	rowmeans = mapslices(x -> mean(x[isfinite.(x)]), rot_result.array; dims = 2)[:]
	for i in 1:nparc
		for j in 1:nrot
			if isnan(mat[i, j])
				mat[i, j] = rowmeans[i]
			end
		end
	end
	result = mapslices(x -> mean(x[isfinite.(x)]), mat.array; dims = 1)
	return (mean = mean(result), std = std(result))
end

function summarize_homogeneity(x::NamedVector)
	return mean(x.array[isfinite.(x.array)])
end

function summarize_homogeneity(real_result::NamedVector, rot_result::NamedMatrix)
	real_mean = summarize_homogeneity(real_result)
	rot_mean, rot_sd = summarize_homogeneity(rot_result)
	return (real_mean - rot_mean) / rot_sd
end

