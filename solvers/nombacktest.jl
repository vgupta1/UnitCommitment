include("readGenerator.jl")
include("readLoads.jl")
include("nomsolver.jl")

genstub = "../Data/AndysGenInstance"
gens = loadISO(genstub)
loadstub = "../Data/ISO-NE\ Load\ Data"

dts, vals = readLoads("$loadstub/predtest.csv");
dtstrue, valstrue = readLoads("$loadstub/LoadTest.csv")    
assert(size(vals, 1) == size(valstrue, 1))
filesummary = open(ARGS[1], "w"); filehrs = open(ARGS[2], "w")
numruns = length(ARGS) > 2 ?  int(ARGS[3]) : size(vals, 1)

for iRun = 1:numruns
	forecast = vals[iRun, :]; loads = valstrue[iRun, :]
	m = RobustModel(solver=GurobiSolver(TimeLimit=1*60))
	nom = UCNom(m, gens, 5000)
	status1 = solve(nom, forecast)

	capbyhr = Float64[]
	for ihr = 1:24
		cap = 0.
		for g in values(nom.gendata)
			if getValue(nom.starts[g.name][ihr]) > .5
				cap += getCap(g, ihr)
			end
		end
		push!(capbyhr, cap)
	end
	predvarcost = getVarCost(nom)

	nom2 = secondSolve(nom, loads)
	totShed = sum(map(getValue, nom2.sheds))

	#outputSummary
	if iRun == 1
		writedlm(filesummary, ["Date" "status" "Gap" "FixedCost" "Capacity" "PredVarCost" "VarCost" "Shed"])
	end

	writedlm(filesummary, [dts[iRun] status1 getgap(nom) getStartCost(nom) sum(capbyhr) predvarcost getVarCost(nom2) totShed])

	#dt  DispatchedCapacityByHr   PredVarCostByHr   ShedsByHr 
	if iRun ==1 
		writedlm(filehrs, hcat(dts[iRun], 
						  ["Cap$ihr" for ihr = 1:HRS]', 
						  ["PredVar$ihr" for ihr = 1:HRS]', 
						  ["Shed$ihr" for ihr = 1:HRS]' ) )
	end
	writedlm(filehrs, [dts[iRun] capbyhr' [getVarCost(nom, ihr) for ihr = 1:HRS]' map(getValue, nom.sheds)'])
end

