###
# Generate LP Cluster cuts
###
using Iterators, JuMP
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
penalty             = 5e3

Gamma1              = .14  #chosen via rough jacknife
Gamma2              = .27 
eps                 = .1
numEigs             = 1

GammaBS             = 1.282 * sqrt(HRS)  #These are guestimated.
GammaBound          = 3.0

clusters = [72, 88, 221, 160]

for ix in clusters
	# rm2 = RobustModel(solver=GurobiSolver(MIPGap=1e-3, OutputFlag=1, Method=3, TimeLimit=30*60), 
	# 				  cutsolver=GurobiSolver(OutputFlag=0))
	# alphas, uncs = createPolyUCS(rm2, resids, Gamma1, Gamma2, kappa(eps))
	# aff = UCAff(rm2, gens, penalty, uncs)
	# aff.proj_fcn = eigenProjMatrixData(resids, numEigs)


	rm2 = RobustModel(solver=GurobiSolver(MIPGap=1e-3, OutputFlag=1, Method=3, TimeLimit=30*60), 
					  cutsolver=GurobiSolver(OutputFlag=0))
	alphas, uncs = createBertSimU(rm2, mean(resids, 1), std(resids, 1), GammaBS, GammaBound, false)
	aff = UCAff(rm2, gens, penalty, uncs);
	aff.proj_fcn = identProjMatrixData(resids, numEigs)

	solve(aff, vals[ix, :], report=true, usebox=false, 
				prefer_cuts=true,  lprelax=true, active_cuts=true) 

	println( "Cluster: $ix")
	println( "Active Cuts:")
	ac = rm2.ext[:Robust].activecuts
	println(" Total: $(length(ac))")
	ac = unique(ac)
	println(" Unique: $(length(ac))")

	outpath = "lpbudgetcuts$ix.txt"
	ofile = open(outpath, "w")
	for ix = 1:length(ac)
		writedlm(ofile, transpose(ac[ix][1:HRS]) )  #this is a dangeorus hack
	end
	close(ofile)
end