###
# Back Test for the Nominal Approach
###
using DataFrames, Resampling, Iterators
include("readGenerator.jl")
include("readLoads.jl")
include("nomsolver.jl")
include("robustsolver.jl")
include("UncSets.jl")

gens, scaling       = loadISO("../Data/AndysGenInstance", 1-1e-8)
dts, vals           = readLoads("../Data/ISO-NE Load Data/PredValidate.csv")
dts_true, vals_true = readLoads("../Data/ISO-NE Load Data/LoadValidate.csv")
vals               *= scaling
vals_true          *= scaling
resids              = map(float, vals_true - vals);
penalty             = 5e3
box_ratio           = .933
file_out             = fileout = open(ARGS[1], "w")

for iRun = 1:length(dts)
    #solve a nominal problem as a variance reduction
    m = RobustModel(solver=GurobiSolver(OutputFlag=0))
    nom = UCNom(m, gens, penalty)
    tic()
    status = solve(nom, vals[iRun, :])
    solvetime = toq()
    nom2 = secondSolve(nom, vals_true[iRun, :], report=false)

	if iRun == 1
		writedlm(file_out, ["Date" "status" "Time" "FixedCost" "PredVarCost" "VarCost" "PredShed" "Shed"])
        println(file_out,  ["Date" "status" "Time" "FixedCost" "PredVarCost" "VarCost" "PredShed" "Shed"])
	end
	writedlm(file_out, [dts[iRun] status solvetime getStartCost(nom) getVarCost(nom) getVarCost(nom2) totShed(nom) totShed(nom2)])
    println( file_out, [dts[iRun] status solvetime getStartCost(nom) getVarCost(nom) getVarCost(nom2) totShed(nom) totShed(nom2)])
end

