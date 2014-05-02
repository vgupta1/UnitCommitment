###
# Proj Affine policies
###
# Runs an experiment to compute the incremental cost differences for affine projections 
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
resids              = map(float, vals_true - vals);
kappa(eps)          = sqrt(1/eps - 1)
penalty             = 5e3
ofile               = open(ARGS[1], "w")


Gamma1              = .14  #chosen via rough jacknife
Gamma2              = .27 
eps                 = .1
GammaBS             = 1.282 * sqrt(HRS)  #These are guestimated.
GammaBound          = 3.0

INDXSET = [116, 51, 239, 118, 73, 59, 218, 220, 99, 227]  #obtained by 10 k means

#solve the nominal problems and the robust problems for variance reduction
tic()
nomVals = Dict{Int, Float64}()
warmStartsUCS = Dict{Int, WarmStartInfo}()
warmStartsBudget = Dict{Int, WarmStartInfo}()
for ix in INDXSET
    m = RobustModel(solver=GurobiSolver(OutputFlag=0, MIPGap=1e-3, TimeLimit=15*60))
    nom = UCNom(m, gens, penalty)
    solve(nom, vals[ix, :])
    nom2 = secondSolve(nom, vals_true[ix, :], report=false)
    nomVals[ix] = getObjectiveValue(nom2.m)

	#solve a robust problem for a UCS mipstart
	m = RobustModel(solver=GurobiSolver(OutputFlag=0, MIPGap=1e-3, TimeLimit=15*60))
	alphas, uncs = createPolyUCS(m, resids, Gamma1, Gamma2, kappa(eps))
	rob = UCRob(m, gens, penalty, uncs)
	solve(rob, vals[ix, :], report=false)
	warmStartsUCS[ix] = copyWarmStart(rob, WarmStartInfo())

	#solve a robust problem for the budget mipstart
	m = RobustModel(solver=GurobiSolver(OutputFlag=0, MIPGap=1e-3, TimeLimit=15*60))
	alphas, uncs = createBertSimU(m, mean(resids, 1), cov(resids), GammaBS, GammaBound, false)
	rob = UCRob(m, gens, penalty, uncs)
	solve(rob, vals[ix, :], report=false)
	warmStartsBudget[ix] = copyWarmStart(rob, WarmStartInfo())
end
println("Setting up MipStart Stuff", toc() )

writedlm(ofile, ["Set" "NumEigs" "Ind" "Status" "TotCost" "StartCost" "VarCost" "Shed" "NomTotCost" "SolveTime" "Gap"])
# iterate over directions
for (numDirs, ix) in product([1 2 3 5 7 10],INDXSET)
	#solve a UC 
	rm2 = RobustModel(solver=GurobiSolver(MIPGap=1e-3, OutputFlag=0, Method=3, TimeLimit=60*60), 
					  cutsolver=GurobiSolver(OutputFlag=0))
	alphas, uncs = createPolyUCS(rm2, resids, Gamma1, Gamma2, kappa(eps));
	aff = UCAff(rm2, gens, penalty, uncs);
	aff.proj_fcn = eigenProjMatrixData(resids, numDirs)
	aff.warmstart = warmStartsUCS[ix]
	ucstime = 0; tic()
	status = solve(aff, vals[ix, :], report=false, usebox=false, prefer_cuts=true) 
	ucstime += toq()
	aff2 = secondSolve(aff, vals_true[ix, :], report=false)

	#write the adjusted values
	writedlm(ofile, ["UCS" numDirs ix status getObjectiveValue(aff2.m) getStartCost(aff) getVarCost(aff2) totShed(aff2) nomVals[ix] ucstime getgap(aff)])
	flush(ofile)

	#solve the budget
	rm2 = RobustModel(solver=GurobiSolver(MIPGap=1e-3, OutputFlag=0, Method=3, TimeLimit=60*60), 
					  cutsolver=GurobiSolver(OutputFlag=0))
	alphas, uncs = createBertSimU(rm2, mean(resids, 1), cov(resids), GammaBS, GammaBound, false)
	aff = UCAff(rm2, gens, penalty, uncs);
	aff.proj_fcn = identProjMatrixData(resids, numDirs)
	aff.warmstart = warmStartsBudget[ix]
	budtime = 0; tic()
	status = solve(aff, vals[ix, :], report=false, usebox=false, prefer_cuts=true) 
	budtime += toq()
	aff2 = secondSolve(aff, vals_true[ix, :], report=false)

	#write the adjusted values
	writedlm(ofile, ["Budget" numDirs ix status getObjectiveValue(aff2.m) getStartCost(aff) getVarCost(aff2) totShed(aff2) nomVals[ix] budtime getgap(aff) ])
	flush(ofile)
	println(ix, "  ", numDirs, "  ", ucstime, "  ", budtime)
end

