###
# By Hand jackknife for Affine Methods
###
using Iterators

include("readGenerator.jl")
include("readLoads.jl")
include("nomsolver.jl")
include("robustsolver.jl")
include("UncSets.jl")
include("adaptivesolver.jl")

gens, scaling       = loadISO("../Data/AndysGenInstance", 1)
dts, vals           = readLoads("../Data/ISO-NE Load Data/PredTest.csv")
dts_true, vals_true = readLoads("../Data/ISO-NE Load Data/LoadTest.csv")
vals               *= scaling
vals_true          *= scaling
resids              = map(float, vals_true - vals)
kappa(eps)          = sqrt(1/eps - 1)
penalty             = 5e3
numEigs             = 2
ofile               = open(ARGS[1], "a")

#VG Revisit these....
eps_grid = [.7, .8, .9, .95]

#Corresponds to delta/2 = 
#               95%       90%       85%       80%       75% 
g1_grid = scaling * [0.5993808 0.5171897 0.4588753 0.4105277 0.3695186] 
g2_grid = scaling * scaling * [4.974080  4.190187  3.713862  3.349935  3.035816]

#limit life to the indices that matter
INDXSET = [72, 88, 107, 160, 1, 2, 5, 204, 78]
resids    = resids[INDXSET, :]
vals      = vals[INDXSET, :]
vals_true = vals_true[INDXSET, :]
dts       = dts[INDXSET]
dts_true  = dts_true[INDXSET]

tic()
for (eps, g1, g2) in product(eps_grid, g1_grid, g2_grid)
	#now iterate over everyone
	costs = Float64[]
	for ix = 1:size(resids, 1)
		#Solve one robust model that will be used to warmstart everyone
		m = RobustModel(solver=GurobiSolver(OutputFlag=0))
		alphas, uncs = createPolyUCS(m, resids, g1, g2, kappa(eps))
		rob = UCRob(m, gens, penalty, uncs)
		solve(rob, vals[1, :], report=false)
		w = WarmStartInfo()
		copyWarmStart(rob, w)


		#train and solve an affine model
		train_indx = [1:ix-1, ix+1:size(resids,1)]
		rm2 = RobustModel(solver=GurobiSolver(MIPGap=1e-3, OutputFlag=0, Method=-1, TimeLimit=60*15), 
						 cutsolver=GurobiSolver(OutputFlag=0))
		alphas, uncs = createPolyUCS(rm2, resids[train_indx, :], g1, g2, kappa(eps), true)
		aff = UCAff(rm2, gens, penalty, uncs)
		aff.proj_fcn = eigenProjMatrixData(resids[train_indx, :], numEigs)

		#Add some random cuts for good luck. ;)
		if length(ARGS) >= 2
			samples = readdlm(open(ARGS[2], "r"), '\t')
			aff.sample_uncs = samples
		end

		try
			solve(aff, vals[ix, :], report=false, usebox=false, prefer_cuts=true)
			aff2     = secondSolve(aff, vals_true[ix, :], report=false)

			#solve a nominal model for the variance reduction
	        m = RobustModel(solver=GurobiSolver(OutputFlag=0))
	        nom = UCNom(m, gens, penalty)
	        solve(nom, vals[ix, :])
	        nom2 = secondSolve(nom, vals_true[ix, :], report=false)

			push!(costs, getObjectiveValue(nom2.m) - getObjectiveValue(aff2.m) )
		catch e
			show(e)
		end
	end
	#tally up what you've got and write to file
	writedlm(ofile, [eps g1 g2 mean(costs) std(costs) length(costs) ])
	flush(ofile)
	println(eps, "  ", g1, "  ", g2)
	toc()
	tic()
end
