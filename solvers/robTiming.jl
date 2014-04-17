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
kappa(eps)          = sqrt(1/eps - 1)
penalty             = 5e3
Gamma1              = .5 * scaling
Gamma2              = 4  * scaling * scaling
eps                 = .1
indx                = int(ARGS[1])
##################################
#solve the robust problem for a warmstart
m = RobustModel(solver=GurobiSolver())
#alphas, uncs = createPolyUCS(m, resids, Gamma1, Gamma2, kappa(eps))
alphas, uncs = createBertSimU(m, resids, 10, Gamma_bound=3)
rob = UCRob(m, gens, penalty, uncs)

println( @elapsed solve(rob, vals[indx, :], report=true) )
