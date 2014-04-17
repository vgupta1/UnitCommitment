###
# Profiling the adaptive Set
###
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
Gamma1              = .5 * scaling
Gamma2              = 4  * scaling * scaling
eps                 = .1
##################################
#solve the robust problem for a warmstart
m = RobustModel(solver=GurobiSolver())
alphas, uncs = createPolyUCS(m, resids, Gamma1, Gamma2, kappa(eps))
rob = UCRob(m, gens, penalty, uncs)
solve(rob, vals[int(ARGS[1]), :], report=true)

w = WarmStartInfo()
copyWarmStart(rob, w)

##################################
## a test function
function testRun()
	rm2 = RobustModel(solver=GurobiSolver(MIPGap=5e-3, Method=2))  #MIPGap=1e-2, OptimalityTol=1e-4)
	alphas, uncs = createPolyUCS(rm2, resids, Gamma1, Gamma2, kappa(eps));
	aff = UCAff(rm2, gens, penalty, uncs);
	aff.proj_fcn = eigenProjMatrixData(resids, 1)
	aff.warmstart = w
	solve(aff, vals[int(ARGS[1]), :], report=true, usebox=false) #, lprelax=true
end

println( @elapsed testRun() )

# .3 small, eigs = 1, pen = 1
#  LP: 1451.48  total:1497

# .3 small, eigs = 1, pen = 5e3
#  LP: 686.71  total:???
#650272.2295467586
#127.50386707066602

#solving the mip versin of above:
# Root Relaxation: 6.5049432e+05  (Seems to be some noise)

