
export make_cov_corr, homogeneity_test

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
		p::Parcel, cov_corr::Matrix; criteria::Function = x -> true
	)
	verts = @view collapse(vertices(p), p.surface)
	return criteria(p, verts) ? homogeneity(cov_corr[verts, verts]) : NaN
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
		ks::Vector{T} = collect(keys(px), criteria::Function = x -> true)
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
		ks::Vector{T}, criteria::Function = x -> true
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


