

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





