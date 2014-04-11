####
# A Hindsight test
####
# For each day in the test period, compute the nominal vs. the forecast, and what it would receive
# Also compute the perfect foresight value for the difference...
include("readGenerator.jl")
include("readLoads.jl")
include("nomsolver.jl")

genstub = "../Data/AndysGenInstance"
gens = loadISO(genstub)
loadstub = "../Data/ISO-NE\ Load\ Data"

dts, vals = readLoads("$loadstub/predtest.csv");
dtstrue, valstrue = readLoads("$loadstub/LoadTest.csv")    
assert(size(vals, 1) == size(valstrue, 1))
fileout = open(ARGS[1], "w");
numruns = length(ARGS) > 1 ?  int(ARGS[2]) : size(vals, 1)

for iRun = 1:numruns
	forecast = vals[iRun, :]; loads = valstrue[iRun, :]
	m = RobustModel(solver=GurobiSolver(TimeLimit=5*60))
	nom = UCNom(m, gens, 5000)
	status1 = solve(nom, forecast)
	nom2 = secondSolve(nom, loads)
	totShed = sum(map(getValue, nom2.sheds))

	nomhind = UCNom(RobustModel(solver=GurobiSolver(TimeLimit=5*60)), gens, 5000)
	status2 = solve(nomhind, loads)
	totShed2 = sum(map(getValue, nomhind.sheds))

	#outputSummary
	if iRun == 1
		writedlm(fileout, ["Date" "status" "Gap" "FixedCost" "VarCost" "Shed" "HindFixed" "HindVar" "HindShed"])
	end
	writedlm(fileout, [dts[iRun] status1 getgap(nom) getStartCost(nom) getVarCost(nom2) totShed getStartCost(nomhind) getVarCost(nomhind) totShed2])
end


