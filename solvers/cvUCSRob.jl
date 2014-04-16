###
# Cross Validation for the UCS Robust set
###
using DataFrames, Resampling, Iterators

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

########################
mydf = DataFrame(1:size(resids, 1))
trainUCS(df) = (mean(resids[df[1], :], 1), cov(resids[df[1], :]))
function testUCSRob(df, params, epsilon, Gamma1, Gamma2)
    mu, Sigma = params
    avg = 0
    for i in df[1]
        rm = RobustModel(solver=GurobiSolver(MIPGap=5e-3, OutputFlag=0))
        alphas, uncs = createPolyUCS(rm, mu, Sigma, Gamma1, Gamma2, kappa(epsilon), true)
        rob = UCRob(rm, gens, 5e3, uncs)
        solve(rob, vals[i, :], usebox=false, report=false)

        #now solve a second model to get real costs
        nom2 = secondSolve(rob, vals_true[i, :], report=false);        
        avg += getObjectiveValue(nom2.m)
    end
    avg/size(df, 1)
end    
##########################
#VG Revisit these....
eps_grid = linspace(1e-3, 1-1e-3, 10)
g1_grid  = linspace(0, 1, 10)
g2_grid  = linspace(0, 2, 10)

##DEBUG
eps_grid = [.1, .2]
g1_grid  = [.5, ]
g2_grid  = [.5, ]

ofile = open(ARGS[1], "a")

for (eps, g1, g2) in product(eps_grid, g1_grid, g2_grid)
	testUCSRob_(df, params) = testUCSRob(df, params, eps, g1, g2)
	try
		dummy, results = kfold_crossvalidate(DataFrame(mydf), trainUCS, testUCSRob_, 10)
		#write a value and flush it
		writedlm(ofile, [eps g1 g2 mean(results) std(results)])
		flush(ofile)
	catch
	end
end