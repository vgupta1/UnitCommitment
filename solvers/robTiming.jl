###
# Profiling the robust Set
###
include("readGenerator.jl")
include("readLoads.jl")
include("nomsolver.jl")
include("robustsolver.jl")
include("UncSets.jl")

gens, scaling       = loadISO("../Data/AndysGenInstance", 1-1e-8)
dts, vals           = readLoads("../Data/ISO-NE Load Data/PredTest.csv")
dts_true, vals_true = readLoads("../Data/ISO-NE Load Data/LoadTest.csv")
vals               *= scaling
vals_true          *= scaling
resids              = map(float, vals_true - vals);

##VG DEBUG
resid_mu = mean(resids, 1)
for i = 1:size(resids, 1)
	resids[i, :] -= resid_mu
end

kappa(eps)          = sqrt(1/eps - 1)
penalty             = 5e3
Gamma1              = .5 * scaling
Gamma2              = 4  * scaling * scaling
eps                 = .1
indx                = int(ARGS[1])
##################################
#solve the robust problem for a warmstart
m = RobustModel(solver=GurobiSolver(Presolve=0), cutsolver=GurobiSolver(OutputFlag=0))
alphas, uncs = createPolyUCS(m, resids, Gamma1, Gamma2, kappa(eps))
#
alphas, uncs = createBertSimU(m, resids, Gamma1, Gamma_bound=Gamma2)
rob = UCRob(m, gens, penalty, uncs)
println(indx)
println( @elapsed solve(rob, vals[indx, :], report=true, 
			prefer_cuts= (ARGS[2]=="true")) )


println("\n \nNominal Value:\n")
nom = UCNom(RobustModel(solver=GurobiSolver(OutputFlag=0)), gens, penalty)
solve(nom, vals[indx, :])
println(getObjectiveValue(nom.m))
