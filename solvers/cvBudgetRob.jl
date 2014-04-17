###
# Cross Validation for the Budget Robust set
###
using DataFrames, Resampling, Iterators

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

########################
mydf = DataFrame(1:size(resids, 1))
trainUCS(df) = (mean(resids[df[1], :], 1), cov(resids[df[1], :]))
function testUCSRob(df, params, epsilon, Gamma1, Gamma2)
    mu, Sigma = params
    avg = 0
    for i in df[1]
        rm = RobustModel(solver=GurobiSolver(MIPGap=5e-3, OutputFlag=0))
        alphas, uncs = createBertSimU(rm, mu, Sigma, Gamma1, Gamma2, false)
        rob = UCRob(rm, gens, penalty, uncs)
        solve(rob, vals[i, :], usebox=false, report=false)

        #now solve a second model to get real costs
        nom2 = secondSolve(rob, vals_true[i, :], report=false);        
        avg += getObjectiveValue(nom2.m)
    end
    avg/size(df, 1)
end    
##########################
#VG Revisit these....
#Corresponds to delta/2 = 
#               95%       90%       85%       80%       75% 
g1_grid = sqrt(24) * linspace(.25, 3, 5) 
g2_grid = [3,]
eps_grid = [.05 .1 .15 .2 .25]

ofile = open(ARGS[1], "a")
for eps in eps_grid
    for g1 in g1_grid
        g2 = g2_grid[1]  #VG Hacky
    	testUCSRob_(df, params) = testUCSRob(df, params, eps, g1, g2)
    	#try
    		dummy, results = kfold_crossvalidate(DataFrame(mydf), trainUCS, testUCSRob_, 5)
    		#write a value and flush it
    		writedlm(ofile, [eps g1 g2 mean(results) std(results)])
    		flush(ofile)
            println(eps, "  ", g1, "  ", g2)
    	#catch
    	#end
    end
end
