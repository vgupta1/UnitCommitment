###
# Cross Validation for the Robust Set
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
penalty             = 5e3

########################
function testM(df, params)
    avg = 0
    for i in df[1]
        rm = RobustModel(solver=GurobiSolver(OutputFlag=0))
        uncs = createUM(rm, params...)
        rob = UCRob(rm, gens, penalty, uncs)
        solve(rob, vals[i, :], usebox=false, report=false, prefer_cuts=true)

        #resolve to get real costs
        rob2 = secondSolve(rob, vals_true[i, :], report=false)        
    
        #solve a nominal problem as a variance reduction
        m = RobustModel(solver=GurobiSolver(OutputFlag=0))
        nom = UCNom(m, gens, penalty)
        solve(nom, vals[i, :])
        nom2 = secondSolve(nom, vals_true[i, :], report=false)
        avg += getObjectiveValue(rob2.m) - getObjectiveValue(nom2.m)

    end
    avg/size(df, 1)
end    
##########################
ofile = open(ARGS[1], "a")

mydf  = DataFrame(1:size(resids, 1))
N, d = size(resids)
s_grid = [int(N * p) for p in linspace(.7, 1, 10)]

for r in linspace(.7, 1., 10)
	trainUM(df) = calcUMBounds(resids[df[1], :], r)
	try
		dummy, results = kfold_crossvalidate(mydf, trainUM, testM, 5)
		#write a value and flush it
		writedlm(ofile, [r  mean(results) std(results)])
		flush(ofile)
        println(r, "  ", mean(results))
	catch e
        show(e)
	end
end
