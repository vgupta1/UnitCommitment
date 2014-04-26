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
penalty             = 5e3
ofile               = open(ARGS[1], "a")

G1_grid = sqrt(24) * linspace(1, 3, 5)
G2_grid = linspace(1.75, 3.25, 4)

INDXSET = [116, 51, 239, 118, 73, 59, 218, 220, 99, 227]  #obtained by 10 k means
resids    = resids[INDXSET, :]
vals      = vals[INDXSET, :]
vals_true = vals_true[INDXSET, :]
dts       = dts[INDXSET]
dts_true  = dts_true[INDXSET]

tic()
for (g1, g2) in product(G1_grid, G2_grid)
	#Solve one robust model that will be used to warmstart everyone
	m = RobustModel(solver=GurobiSolver(OutputFlag=0))
	alphas, uncs = createBertSimU(m, mean(resids, 1), std(resids, 1), g1, g2, false)

	rob = UCRob(m, gens, penalty, uncs)
	solve(rob, vals[1, :], report=false)
	w = WarmStartInfo()
	copyWarmStart(rob, w)

	#now iterate over everyone
	costs = Float64[]
	for ix = 1:size(resids, 1)
		#train and solve an affine model
		train_indx = [1:ix-1, ix+1:size(resids,1)]
		rm2 = RobustModel(solver=GurobiSolver(MIPGap=1e-3, OutputFlag=0, Method=0, TimeLimit=60*15), 
						 cutsolver=GurobiSolver(OutputFlag=0))
		alphas, uncs = createBertSimU(rm2, mean(resids[train_indx, :], 1), cov(resids[train_indx, :]), g1, g2, false)

		aff = UCAff(rm2, gens, penalty, uncs)
		aff.proj_fcn = hr -> eye(HRS)

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
	writedlm(ofile, [g1 g2 mean(costs) std(costs) length(costs) ])
	flush(ofile)
	println(g1, "  ", g2, "  ", toq())
	tic()
end
