###
# Cross Validation for the UCS Aff set
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
resids              = map(float, vals_true - vals)
kappa(eps)          = sqrt(1/eps - 1)
penalty             = 5e3
numEigs             = int(ARGS[2])

########################
mydf = DataFrame(1:size(resids, 1))
trainUCS(df) = (mean(resids[df[1], :], 1), cov(resids[df[1], :]), eigenProjMatrixData(resids[df[1], :], numEigs))
function testUCSAff(df, params, epsilon, Gamma1, Gamma2)
    mu, Sigma, proj_fcn = params
    rob_avg   = 0.
    aff_avg   = 0.
    rob_gap   = 0.
    aff_gap   = 0.
    for i in df[1]
    	#First solve a Robust Model for the warm start
        rm = RobustModel(solver=GurobiSolver(MIPGap=5e-3, OutputFlag=0, TimeLimit=60*20))
        alphas, uncs = createPolyUCS(rm, mu, Sigma, Gamma1, Gamma2, kappa(epsilon), true)
        rob = UCRob(rm, gens, penalty, uncs)
        solve(rob, vals[i, :], usebox=false, report=false)

    	#Copy over the warmstart info
		w = WarmStartInfo()
		copyWarmStart(rob, w)

    	#solve it for real cost
    	nom2     = secondSolve(rob, vals_true[i, :], report=false)
    	rob_avg += getObjectiveValue(nom2.m)
    	# rob_gap += getgap(rob)

    	#solve an affine model
    	rm2 = RobustModel(solver=GurobiSolver(MIPGap=5e-3, OutputFlag=0, TimeLimit=60*20))
		alphas, uncs = createPolyUCS(rm2, mu, Sigma, Gamma1, Gamma2, kappa(epsilon), true)
		aff = UCAff(rm2, gens, penalty, uncs);
		aff.proj_fcn = proj_fcn
		aff.warmstart = w
		solve(aff, vals[i, :], report=false, usebox=false)

		#solve once more for the real cost
		nom2     = secondSolve(aff, vals_true[i, :], report=false)
		aff_avg += getObjectiveValue(nom2.m)
    	# aff_gap += getgap(aff)
    end

    #return the average cost of the rob and the avg cost of the aff
    n = size(df, 1)
    rob_avg/n, aff_avg/n #, rob_gap/n, aff_gap/n
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
	testUCSAff_(df, params) = testUCSAff(df, params, eps, g1, g2)
	# try
		dummy, results = kfold_crossvalidate(DataFrame(mydf), trainUCS, testUCSAff_, 5)
		#write a value and flush it
		rob_res, aff_res = zip(results...)
		rob_res = [rob_res...]
		aff_res = [aff_res...]
		writedlm(ofile, [eps g1 g2 mean(rob_res) std(rob_res) mean(aff_res) std(aff_res)])
		flush(ofile)
		println(eps, g1, g2)
	# catch
	# end
end

##VG fix the grids.

