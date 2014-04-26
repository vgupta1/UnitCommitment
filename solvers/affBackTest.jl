###
# Back Test for the Affine Sets
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
resids              = map(float, vals_true - vals);

dts, vals           = readLoads("../Data/ISO-NE Load Data/PredValidate.csv")
dts_true, vals_true = readLoads("../Data/ISO-NE Load Data/LoadValidate.csv")
vals               *= scaling
vals_true          *= scaling

penalty             = 5e3
file_out            = fileout = open(ARGS[1], "w")

numEigs             = 1
Gamma1              = 0.5171897 * scaling
Gamma2              = 4.190187   * scaling * scaling
kappa                = sqrt(1/.1 - 1)
GammaBS             = 1.282 * sqrt(24)
GammaBound          = 3.0

clustermap = readdlm("../results/All/ValidateSet/validateClusterMap2.csv", '\r')
println(size(clustermap))
clustermap = [int(c) for c in clustermap]
println(size(clustermap))

for iRun = 1:length(dts)
    if iRun == 1
        writedlm(file_out, ["Date" "Method" "status" "Time" "FixedCost" "PredVarCost" "VarCost" "PredShed" "Shed"])
        # println(file_out,  ["Date" "Method" "status" "Time" "FixedCost" "PredVarCost" "VarCost" "PredShed" "Shed"])
    end

    cluster = clustermap[iRun]

    #solve the nominal problem to use as variance reduction
    m = RobustModel(solver=GurobiSolver(OutputFlag=0, MIPGap=1e-3, TimeLimit=60*15))
    nom = UCNom(m, gens, penalty)
    solvetime = 0; tic()
    status = solve(nom, vals[iRun, :])
    solvetime = toq()
    nom2 = secondSolve(nom, vals_true[iRun, :], report=false)
    writedlm(file_out, [dts[iRun] "Nominal" status solvetime getStartCost(nom) getVarCost(nom) getVarCost(nom2) totShed(nom) totShed(nom2)])
    println( file_out, [dts[iRun] "Nominal" status solvetime getStartCost(nom) getVarCost(nom) getVarCost(nom2) totShed(nom) totShed(nom2)])
    flush(file_out)

    #Solve a robust problem for warm start info
    #VG This should be change to properly use a warm start file
    # This biases the timings  a little.
    rm = RobustModel(solver=GurobiSolver(OutputFlag=0, MIPGap=1e-2, TimeLimit=60*15))
    uncs = createUM(rm, resids, .9333)
    rob    = UCRob(rm, gens, penalty, uncs)
    solve(rob, vals[iRun, :], usebox=false, report=false, prefer_cuts=true)
    w = copyWarmStart(rob, WarmStartInfo())

    println("Got through Robust")


    #Solve against UCS
    rm2 = RobustModel(solver=GurobiSolver(MIPGap=1e-3, OutputFlag=0, Method=3, TimeLimit=60*30), 
                     cutsolver=GurobiSolver(OutputFlag=0))
    alphas, uncs = createPolyUCS(rm2, resids, Gamma1, Gamma2, kappa, true)
    aff = UCAff(rm2, gens, penalty, uncs)
    aff.proj_fcn = eigenProjMatrixData(resids, numEigs)

    ##VG Change this to use proper cluster map
    aff.sample_uncs = readdlm(open("../results/Size_10/cuts$cluster.txt", "r"), '\t')
    aff.warmstart = w

    try
        solvetime = 0; tic()
        solve(aff, vals[iRun, :], report=false, usebox=false, prefer_cuts=true)
        solvetime = toq()
        aff2     = secondSolve(aff, vals_true[iRun, :], report=false)
        writedlm(file_out, [dts[iRun] "UCS" status solvetime getStartCost(aff) getVarCost(aff) getVarCost(aff2) totShed(aff) totShed(aff2)])
        # println( file_out, [dts[iRun] "UCS" status solvetime getStartCost(aff) getVarCost(aff) getVarCost(aff2) totShed(aff) totShed(aff2)])
        flush(file_out)
    catch e
            show(e)
    end

    #Solve against UB
    rm2 = RobustModel(solver=GurobiSolver(MIPGap=1e-3, OutputFlag=1, Method=3, TimeLimit=30*60), 
                      cutsolver=GurobiSolver(OutputFlag=0))
    alphas, uncs = createBertSimU(rm2, mean(resids, 1), cov(resids), GammaBS, GammaBound, false)
    aff = UCAff(rm2, gens, penalty, uncs);
    aff.proj_fcn = identProjMatrixData(resids, numEigs)

    aff.sample_uncs = readdlm(open("../results/Size_10/budgetcut$cluster.txt", "r"), '\t')    
    aff.warmstart = w


    try
        solvetime = 0; tic()
        status = solve(aff, vals[iRun, :], report=false, usebox=false, prefer_cuts=true) 
        solvetime = toq()
        aff2     = secondSolve(aff, vals_true[iRun, :], report=false)
        writedlm(file_out, [dts[iRun] "UBudget" status solvetime getStartCost(aff) getVarCost(aff) getVarCost(aff2) totShed(aff) totShed(aff2)])
        # println( file_out, [dts[iRun] "UBudget" status solvetime getStartCost(aff) getVarCost(aff) getVarCost(aff2) totShed(aff) totShed(aff2)])
        flush(file_out)
    catch e
            show(e)
    end

end
close(file_out)
