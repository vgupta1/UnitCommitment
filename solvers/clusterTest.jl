###
# Clustering Test
###
using Iterators
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
ofile               = open(ARGS[1], "a")

#VG Revisit these
Gamma1              = .5 * scaling
Gamma2              = 4  * scaling * scaling
eps                 = .1
GammaBS             = sqrt(HRS) * 2.5
GammaBound          = 3


#map the above indices to their appropriate clusters (4)
clusters = [72, 88, 221, 160]
INDXSET = [116, 51, 239, 118, 73, 59, 218, 220, 99, 227]  #obtained by 10 k means
clustermap = [72  72  88  88  221 221 160 160 160 160]  #only pertains ot indx set
clustermap = Dict(INDXSET, clustermap)

#solve the nominal problems and the robust problems for variance reduction
tic()
nomVals = Dict{Int, Float64}()
warmStartsUCS = Dict{Int, WarmStartInfo}()
warmStartsBudget = Dict{Int, WarmStartInfo}()
for ix in INDXSET
    m = RobustModel(solver=GurobiSolver(OutputFlag=0, MIPGap=1e-3, TimeLimit=15*60))
    nom = UCNom(m, gens, penalty)
    solve(nom, vals[ix, :])
    nom2 = secondSolve(nom, vals_true[ix, :], report=false)
    nomVals[ix] = getObjectiveValue(nom2.m)

	#solve a robust problem for a UCS mipstart
	m = RobustModel(solver=GurobiSolver(OutputFlag=0, MIPGap=1e-3, TimeLimit=15*60))
	alphas, uncs = createPolyUCS(m, resids, Gamma1, Gamma2, kappa(eps))
	rob = UCRob(m, gens, penalty, uncs)
	solve(rob, vals[ix, :], report=false)
	warmStartsUCS[ix] = copyWarmStart(rob, WarmStartInfo())

	# #solve a robust problem for the budget mipstart
	# m = RobustModel(solver=GurobiSolver(OutputFlag=0, 2*60))
	# alphas, uncs = createBertSimU(m, mean(resids, 1), cov(resids), GammaBS, GammaBound, false)
	# rob = UCRob(m, gens, penalty, uncs)
	# solve(rob, vals[ix, :], report=false)
	# warmStartsBudget[ix] = copyWarmStart(rob, WarmStartInfo())
end
println("Setting up MipStart Stuff", toc() )

writedlm(ofile, ["Set" "Mode" "Indx" "Status" "TotCost" "StartCost" "VarCost" "Shed" "NomTotCost" "SolveTime"])
for (mode, ix) in product([:true :random :none], INDXSET)
	cluster = clustermap[ix]

	#solve a UC 
	rm2 = RobustModel(solver=GurobiSolver(MIPGap=1e-3, OutputFlag=1, Method=3, TimeLimit=45*60), 
					  cutsolver=GurobiSolver(OutputFlag=0))
	alphas, uncs = createPolyUCS(rm2, resids, Gamma1, Gamma2, kappa(eps));
	aff = UCAff(rm2, gens, penalty, uncs);
	aff.proj_fcn = eigenProjMatrixData(resids, 1)
	aff.warmstart = warmStartsUCS[ix]
	
	if mode == :true
		aff.sample_uncs = readdlm(open("../results/Size_10/cuts$cluster.txt", "r"), '\t')
	elseif mode ==:random
		#choose a random number that's not cluster the given cluster
		indx = ceil(rand() * 3)
		cluster = filter(x-> x != cluster, clusters)[indx]
		aff.sample_uncs = readdlm(open("../results/Size_10/cuts$cluster.txt", "r"), '\t')
	end		
	ucstime = 0
	tic()
	status = solve(aff, vals[ix, :], report=false, usebox=false, prefer_cuts=true) 
	ucstime += toq()

	#do a second solve
	aff2 = secondSolve(aff, vals_true[ix, :], report=false)

	#write the adjusted values
	writedlm(ofile, ["UCS" mode ix status getObjectiveValue(aff2.m) getStartCost(aff) getVarCost(aff2) totShed(aff2) nomVals[ix] ucstime])
	flush(ofile)
	println(["UCS" mode ix status getObjectiveValue(aff2.m) getStartCost(aff) getVarCost(aff2) totShed(aff2) nomVals[ix] ucstime])


	# #solve the budget
	# rm2 = RobustModel(solver=GurobiSolver(MIPGap=1e-3, OutputFlag=1, Method=3, TimeLimit=2*60), 
	# 				  cutsolver=GurobiSolver(OutputFlag=0))
	# alphas, uncs = createBertSimU(rm2, mean(resids, 1), cov(resids), GammaBS, GammaBound, false)
	# aff = UCAff(rm2, gens, penalty, uncs);
	# aff.proj_fcn = identProjMatrixData(resids, numEigs)

	# aff.warmstart = warmStartsBudget[ix]
	# aff.sample_uncs = readdlm(open("../results/Size_10/budgetcut$cluster.txt", "r"), '\t')
	# budtime = 0
	# tic()
	# status = solve(aff, vals[ix, :], report=false, usebox=false, prefer_cuts=true) 
	# budtime += toq()

	# #write the adjusted values
	# writedlm(ofile, [numDirs ix status getObjectiveValue(rm2) nomVals[ix] budtime ])
	# flush(ofile)
	# println(ix, "  ", numDirs, "  ", ucstime, "  ", budtime)
end







