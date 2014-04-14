###
# Profiling the adaptive Set
###

include("readGenerator.jl")
include("readLoads.jl")
include("nomsolver.jl")
include("UncSets.jl")
include("adaptivesolver.jl")

gens = loadISO("../Data/AndysGenInstance")
dts, vals = readLoads("../Data/ISO-NE Load Data/PredTest.csv")
dts_true, vals_true = readLoads("../Data/ISO-NE Load Data/LoadTest.csv")
resids = map(float, vals_true - vals);


D, V = eig(cov(resids));
numEigs = 5
projMatrix = V[:, (end-numEigs):end]'  #should thsi be scaled?
kappa(eps) = sqrt(1/eps - 1)

function testRun()
	rm2 = RobustModel(solver=GurobiSolver())
	alphas, uncs = createPolyUCS(rm2, resids, .1, .1, kappa(.1));
	aff = UCAff(rm2, gens, 5e3, uncs);
	aff.proj_fcn = (hr)-> projMatrix
	solve(aff, vals[1, :], report=true)
end

println( @elapsed testRun() )


