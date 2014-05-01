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
file_out            = fileout = open(ARGS[1], "a")

numEigs             = 1
Gamma1              = .14   #From rough jacknife analysis
Gamma2              = .27
kappa                = sqrt(1/.1 - 1)

GammaBS             = 1.282 * sqrt(24)  #These are guestimated.
GammaBound          = 3.0

box_size            =.9333  #From 5-fold crossvalidation

for iRun = 1:length(dts)
    if iRun == 1
        writedlm(file_out, ["Date" "Method" "status" "Time" "FixedCost" "PredVarCost" "VarCost" "PredShed" "Shed"])
        println(["Date" "Method" "status" "Time" "FixedCost" "PredVarCost" "VarCost" "PredShed" "Shed"])
    end

    #Nominal problem
    m = RobustModel(solver=GurobiSolver(OutputFlag=0, MIPGap=1e-3, TimeLimit=60*15))
    nom = UCNom(m, gens, penalty)
    solvetime = 0; tic()
    status = solve(nom, vals[iRun, :])
    solvetime = toq()
    nom2 = secondSolve(nom, vals_true[iRun, :], report=false)
    writedlm(file_out, [dts[iRun] "Nominal" status solvetime getStartCost(nom) getVarCost(nom) getVarCost(nom2) totShed(nom) totShed(nom2)])
    println([dts[iRun] "Nominal" status solvetime getStartCost(nom) getVarCost(nom) getVarCost(nom2) totShed(nom) totShed(nom2)])
    flush(file_out)

    #Robust problem
    rm = RobustModel(solver=GurobiSolver(OutputFlag=0, MIPGap=1e-3, TimeLimit=60*15))
    uncs = createUM(rm, resids, box_size)
    rob    = UCRob(rm, gens, penalty, uncs)
    solvetime = 0; tic()
    status = solve(rob, vals[iRun, :], usebox=false, report=false, prefer_cuts=true)
    solvetime = toq()
    w = copyWarmStart(rob, WarmStartInfo())
    rob2 = secondSolve(rob, vals_true[iRun, :], report=false)
    writedlm(file_out, [dts[iRun] "Robust" status solvetime getStartCost(rob) getVarCost(rob) getVarCost(rob2) totShed(rob) totShed(rob2)])
    println([dts[iRun] "Robust" status solvetime getStartCost(rob) getVarCost(rob) getVarCost(rob2) totShed(rob) totShed(rob2)])
    flush(file_out)

    #UCS
    rm2 = RobustModel(solver=GurobiSolver(MIPGap=1e-3, OutputFlag=0, Method=3, TimeLimit=60*30), 
                     cutsolver=GurobiSolver(OutputFlag=0))
    alphas, uncs = createPolyUCS(rm2, resids, Gamma1, Gamma2, kappa, true)
    aff = UCAff(rm2, gens, penalty, uncs)
    aff.proj_fcn = eigenProjMatrixData(resids, numEigs)
    aff.warmstart = w

    try
        solvetime = 0; tic()
        solve(aff, vals[iRun, :], report=false, usebox=false, prefer_cuts=true)
        solvetime = toq()
        aff2     = secondSolve(aff, vals_true[iRun, :], report=false)
        writedlm(file_out, [dts[iRun] "UCS" status solvetime getStartCost(aff) getVarCost(aff) getVarCost(aff2) totShed(aff) totShed(aff2)])
        println([dts[iRun] "UCS" status solvetime getStartCost(aff) getVarCost(aff) getVarCost(aff2) totShed(aff) totShed(aff2)])
        flush(file_out)
    catch e
            show(e)
    end

    #Solve against UB
    rm2 = RobustModel(solver=GurobiSolver(MIPGap=1e-3, OutputFlag=0, Method=3, TimeLimit=30*60), 
                      cutsolver=GurobiSolver(OutputFlag=0))
    alphas, uncs = createBertSimU(rm2, mean(resids, 1), cov(resids), GammaBS, GammaBound, false)
    aff = UCAff(rm2, gens, penalty, uncs);
    aff.proj_fcn = identProjMatrixData(resids, numEigs)
    aff.warmstart = w

    try
        solvetime = 0; tic()
        status = solve(aff, vals[iRun, :], report=false, usebox=false, prefer_cuts=true) 
        solvetime = toq()
        aff2     = secondSolve(aff, vals_true[iRun, :], report=false)
        writedlm(file_out, [dts[iRun] "UBudget" status solvetime getStartCost(aff) getVarCost(aff) getVarCost(aff2) totShed(aff) totShed(aff2)])
        println( [dts[iRun] "UBudget" status solvetime getStartCost(aff) getVarCost(aff) getVarCost(aff2) totShed(aff) totShed(aff2)])
        flush(file_out)
    catch e
            show(e)
    end

end
close(file_out)
