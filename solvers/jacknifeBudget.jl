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

gens, scaling       = loadISO("../Data/AndysGenInstance", .1)
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

## Solve the nominal problems for variance reduction
nomVals = Dict{Int, Float64}()
for ix in INDXSET
    m = RobustModel(solver=GurobiSolver(OutputFlag=0, MIPGap=1e-3, TimeLimit=15*60))
    nom = UCNom(m, gens, penalty)
    solve(nom, vals[ix, :])
    nom2 = secondSolve(nom, vals_true[ix, :], report=false)
    nomVals[ix] = getObjectiveValue(nom2.m)
end

writedlm(ofile, ["Gamma1" "GammaBound" "Indx" "TotCost" "StartCost" "VarCost" "Shed" "NomVals"])
tic()
for (g1, g2) in product(G1_grid, G2_grid)
	#now iterate over everyone
	for ix in INDXSET
		#Robust model that will be used to warmstart everyone
		m = RobustModel(solver=GurobiSolver(OutputFlag=0, MIPGap=1e-2, TimeLimit=15*60))
		alphas, uncs = createBertSimU(m, mean(resids, 1), std(resids, 1), g1, g2, false)
		rob = UCRob(m, gens, penalty, uncs)
		solve(rob, vals[ix, :], report=false)
		w = copyWarmStart(rob, WarmStartInfo())

		#train and solve an affine model
		train_indx = [1:ix-1, ix+1:size(resids,1)]
		rm2 = RobustModel(solver=GurobiSolver(MIPGap=1e-3, OutputFlag=0, Method=3, TimeLimit=60*30), 
						 cutsolver=GurobiSolver(OutputFlag=0))
		alphas, uncs = createBertSimU(rm2, mean(resids[train_indx, :], 1), cov(resids[train_indx, :]), g1, g2, false)
		aff = UCAff(rm2, gens, penalty, uncs)
		aff.proj_fcn = hr -> eye(HRS)
		aff.warmstart = w


		try
			solve(aff, vals[ix, :], report=false, usebox=false, prefer_cuts=true)
			aff2     = secondSolve(aff, vals_true[ix, :], report=false)
			writedlm(ofile, [g1 g2 ix getObjectiveValue(aff2.m) getStartCost(aff) getVarCost(aff2) totShed(aff2) nomVals[ix] ])
			flush(ofile)
		catch e
			show(e)
		end
	end
	println(g1, "  ", g2, "  ", toq())
	tic()
end
