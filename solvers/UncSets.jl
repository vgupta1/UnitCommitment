####
# Convenience Fcns for creating Polyhedral Sets
####
using JuMPeR, Gurobi

#data stored row-wise
function createBertSimU(m, resid_data, Gamma; Gamma_bound=2.0, force_pos=false)
	N, d = size(resid_data)
	mu = mean(resid_data, 1); sigma = std(resid_data, 1)

	#compute bounding box.  Needed in other models.
	#Computation assumes that 0 in interior(U), so check
	bContainsZero = true 
	for ix = 1:d
		(mu[ix] / sigma[ix] > Gamma_bound) && (bContainsZero = false)
	end
	(sum(abs(mu ./ sigma)) > Gamma) && (bContainsZero = false)
	! bContainsZero && warn("BertSim Set Not Contain zero")
	minG = min(Gamma, Gamma_bound)
	alphas = [(mu[ix] + minG * sigma[ix])::Float64 for ix = 1:d]
	@defUnc(m, u[1:d])
	@defUnc(m, -Gamma_bound <= z[1:d] <= Gamma_bound)
	for ix = 1:d
		force_pos && addConstraint(u[ix] >= 0)
		addConstraint(m, (u[ix] - mu[ix])/sigma[ix] <= z[ix])
		addConstraint(m, (u[ix] - mu[ix])/sigma[ix] >= -z[ix])
	end
	addConstraint(m, sum(z[:]) <= Gamma)
	return alphas, u
end

function eigenProjMatrix(Sigma, numEigs)
	D, V = eig(Sigma);
	projMatrix = V[:, (end-numEigs):end]'  #should thsi be scaled?
	return (ihr->projMatrix)
end
eigenProjMatrixData(data, numEigs) = eigenProjMatrix(cov(data), numEigs)


#create a simple polyhedral outter approximation to UCS
# VG Refactor the bertsimas/sim norm function
#data stored row-wise
createPolyUCS(rm, resid_data::Matrix{Float64}, Gamma1::Real, Gamma2::Real, kappa::Real, isOuter=true) =
	createPolyUCS(rm, mean(resid_data, 1), cov(resid_data), Gamma1, Gamma2, kappa, isOuter)

function createPolyUCS(rm, mu::Matrix{Float64}, Sigma::Matrix{Float64}, Gamma1::Real, Gamma2::Real, kappa::Real, isOuter)
	! isOuter && error("Inner Approx Not Yet Implemented")
	d = length(mu)
	Covbar = Sigma + Gamma2 * eye(d)
	C = chol(Covbar)::Array{Float64, 2}
	@defUnc(rm, pp[1:d] >= 0)
	@defUnc(rm, pm[1:d] >= 0)
	@defUnc(rm, qp[1:d] >= 0)
	@defUnc(rm, qm[1:d] >= 0)
	@defUnc(rm, w1[1:d] >= 0)
	@defUnc(rm, w2[1:d] >= 0)
	@defUnc(rm, theta[1:2] >= 0)
	for ix = 1:d
		addConstraint(rm, pp[ix] + pm[ix] == w1[ix] + theta[1] )
		addConstraint(rm, qp[ix] + qm[ix] == w2[ix] + theta[2] )
	end		
	addConstraint(rm, sum(w1[:])/sqrt(d) + theta[1] == Gamma1)
	addConstraint(rm, sum(w2[:])/sqrt(d) + theta[2] == kappa)
	@defUnc(rm, us[1:d])
	for ix = 1:d
		cq_ix = sum([C[ix, j] * (qp[j] - qm[j]) for j = 1:d])
		addConstraint(rm, us[ix] == mu[ix] + pp[ix] - pm[ix] + cq_ix)
	end
	#next compute the bounding box
	alphas = Float64[]
	sqrt_d = sqrt(d)
	for ix = 1:d
		alphai = Gamma1 * sqrt_d + sqrt_d * kappa * max(norm(C[:, ix], Inf), norm(C[:, ix], 1) / sqrt_d)
		mu[ix] > alphai && warn("UCS Set may not contain 0")
		push!(alphas, alphai + mu[ix])
	end
	return alphas, us
end
