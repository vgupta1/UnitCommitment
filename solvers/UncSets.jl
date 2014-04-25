####
# Convenience Fcns for creating Polyhedral Sets
####
using JuMPeR, Gurobi
using Distributions


function calcUMBounds(data, s::Int)
	N = size(data, 1)
	(N - s + 1 >= s) && error("Condition on s not met, $s, $N")

	#First sort each column
	data_s = copy(data)
	for jx = 1:size(data_s, 2)
		data_s[:, jx] = sort!(data_s[:, jx])  #Julia doesn't allow sort in place
	end
	lbounds = transpose(data_s[N-s+1, :])
	ubounds = transpose(data_s[s, :])
	return lbounds, ubounds
end
calcUMBounds(data, r::Float64) = calcUMBounds(data, int(r * size(data, 1)))

function calcUMBounds(data, epsilon, delta)
	#compute s
	N, d     = size(data_s)
	eps_d = epsilon / d
	if (1-eps_d)^N > delta / 2d
		warn("Desired epsilon $epsilon, delta $delta  exceed bounds")
		s = N
	else
		binom_dist = Binomial(N, 1 - eps_d )
		s = quantile(binom_dist, 1-delta /2d )
	end
	return calcUMBounds(data, s)
end

function createUM(m, lbounds, ubounds)
	d = length(lbounds)
	@defUnc(m, u[1:d])
	for i = 1:d
		addConstraint(m, lbounds[i] <= u[i])
		addConstraint(m, ubounds[i] >= u[i])
	end
	return u
end
createUM(m, data, epsilon, delta) = createUM(m, calcUMBounds(data, epsilon, delta)...)
createUM(m, data, r::Float64) = createUM(m, calcUMBounds(data, r)...)

##############################

function createBertSimU(m, mu, sigma, Gamma, Gamma_bound, force_pos)
	#compute bounding box.  Needed in other models.
	#Computation assumes that 0 in interior(U), so check
	d = length(mu)
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

#data stored row-wise
createBertSimU(m, resid_data, Gamma; Gamma_bound=2.0, force_pos=false) =
	createBertSimU(m, mean(resid_data, 1), std(resid_data, 1), Gamma, Gamma_bound, force_pos)

##############################

function eigenProjMatrix(Sigma, numEigs)
	D, V = eig(Sigma);
	projMatrix = V[:, (end-numEigs + 1):end]'  #should thsi be scaled?
	return (ihr->projMatrix)
end
eigenProjMatrixData(data, numEigs) = eigenProjMatrix(cov(data), numEigs)

##############################

function identProjMatrixData(resids, numEigs)
	stds = std(resids, 1)
	perm = sortperm(vec(stds))
	diagelem = zeros(size(resids, 2))
	diagelem[ perm[(end-numEigs+1):end] ] = 1
	return(ihr->diagm(diagelem))
end



###############################

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
	@defUnc(rm, us[1:d])
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
	for ix = 1:d
		cq_ix = sum([C[j, ix] * (qp[j] - qm[j]) for j = 1:d])
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

##############################
