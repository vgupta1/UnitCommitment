####
# A Hindsight test
####
# For each day in the test period, compute the nominal vs. the forecast, and what it would receive
# Also compute the perfect foresight value for the difference...
include("readGenerator.jl")
include("readLoads.jl")
include("nomsolver.jl")

genstub  = "../Data/AndysGenInstance"
loadstub = "../Data/ISO-NE\ Load\ Data"

gens, scaling       = loadISO(genstub, .1)
dts, vals           = readLoads("$loadstub/PredTest.csv")
dts_true, vals_true = readLoads("$loadstub/LoadTest.csv")
vals               *= scaling
vals_true          *= scaling
resids              = map(float, vals_true - vals);

fileout = open(ARGS[1], "w");
numruns = length(ARGS) > 1 ?  int(ARGS[2]) : size(vals, 1)

for iRun = 1:numruns
	forecast = vals[iRun, :]; loads = vals_true[iRun, :]
	m = RobustModel(solver=GurobiSolver())
	nom = UCNom(m, gens, 5000)
	status1 = solve(nom, forecast)
	nom2 = secondSolve(nom, loads)

	nomhind = UCNom(RobustModel(solver=GurobiSolver()), gens, 5000)
	status2 = solve(nomhind, loads)

	#outputSummary
	if iRun == 1
		writedlm(fileout, ["Date" "status" "Gap" "FixedCost" "VarCost" "Shed" "HindFixed" "HindVar" "HindShed"])
	end
	writedlm(fileout, [dts[iRun] status1 getgap(nom) getStartCost(nom) getVarCost(nom2) totShed(nom2) getStartCost(nomhind) getVarCost(nomhind) totShed(nomhind)])
	flush(fileout)
end


