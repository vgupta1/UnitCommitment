###
# Back Test for the Robust Set
###
### VG Currently bug in the way residuals handled

using DataFrames, Resampling, Iterators

include("readGenerator.jl")
include("readLoads.jl")
include("nomsolver.jl")
include("robustsolver.jl")
include("UncSets.jl")

gens, scaling       = loadISO("../Data/AndysGenInstance", .1)
dts, vals           = readLoads("../Data/ISO-NE Load Data/PredTest.csv")
dts_true, vals_true = readLoads("../Data/ISO-NE Load Data/LoadTest.csv")
vals               *= scaling
vals_true          *= scaling
resids              = map(float, vals_true - vals);

dts, vals           = readLoads("../Data/ISO-NE Load Data/PredValidate.csv")
dts_true, vals_true = readLoads("../Data/ISO-NE Load Data/LoadValidate.csv")
vals               *= scaling
vals_true          *= scaling


penalty             = 5e3
box_ratio           = .933
file_out             = fileout = open(ARGS[1], "w")

for iRun = 1:length(dts)
    rm = RobustModel(solver=GurobiSolver(OutputFlag=0, MIPGap=1e-3, TimeLimit=60*30))
    uncs = createUM(rm, resids, box_ratio)
    rob    = UCRob(rm, gens, penalty, uncs)

    tic()
    status = solve(rob, vals[iRun, :], usebox=false, report=false, prefer_cuts=true) 
    solvetime   = toq()
    rob2   = secondSolve(rob, vals_true[iRun, :], report=false)        

	if iRun == 1
		writedlm(file_out, ["Date" "status" "Time" "FixedCost" "PredVarCost" "VarCost" "PredShed" "Shed"])
        println(["Date" "status" "Time" "FixedCost" "PredVarCost" "VarCost" "PredShed" "Shed"])
        flush(file_out)
	end
	writedlm(file_out, [dts[iRun] status solvetime getStartCost(rob) getVarCost(rob) getVarCost(rob2) totShed(rob) totShed(rob2)])
    println([dts[iRun] status solvetime getStartCost(rob) getVarCost(rob) getVarCost(rob2) totShed(rob) totShed(rob2)])
    flush(file_out)
end

