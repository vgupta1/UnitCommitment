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
        alphas, uncs = createPolyUCS(rm, mu, Sigma, Gamma1, Gamma2, kappa(epsilon), true)
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
g1_grid = scaling * [0.5993808 0.5171897 0.4588753 0.4105277 0.3695186] 
g2_grid = scaling * scaling * [4.974080  4.190187  3.713862  3.349935  3.035816]
eps_grid = [.05 .1 .15 .2 .25]

ofile = open(ARGS[1], "a")
for eps in eps_grid
    for iG = 1:length(g1_grid)
        g1 = g1_grid[iG]; g2 = g2_grid[iG]
    	testUCSRob_(df, params) = testUCSRob(df, params, eps, g1, g2)
    	try
    		dummy, results = kfold_crossvalidate(DataFrame(mydf), trainUCS, testUCSRob_, 5)
    		#write a value and flush it
    		writedlm(ofile, [eps g1 g2 mean(results) std(results)])
    		flush(ofile)
            println(eps, "  ", g1, "  ", g2)
    	catch
    	end
    end
end
