###
# Profiling the adaptive Set
###

include("readGenerator.jl")
include("readLoads.jl")
include("nomsolver.jl")
include("UncSets.jl")
include("adaptivesolver.jl")

gens = loadISO("../Data/AndysGenInstance")

println("num Gens:\t", length(gens))
hr13Cap = sum([getCap(g, 13) for g in values(gens)])
println("HR 13 Cap:\t", hr13Cap)
function filtGens!(gens)
    for g in values(gens)
        if g.startcost < .1
            g.startcost = 0.
        end
    end
    #These generators seemingly have bad data...
    delete!(gens, "VLCO-I_ASCUTNEY_26_22")
    delete!(gens, "VLCO-I_N_RUTLND_44_22")

    #drop the wind and hydro for now
    for g in values(gens)
        if g.fueltype == "Wind"
            delete!(gens, g.name)
        elseif g.fueltype=="Hydro"
            delete!(gens, g.name)
        end
    end
    
    #halve the steam, CT and Diesel
    srand(8675309)
    p = .5
    for g in values(gens)
        if g.fueltype =="Steam"
            rand() > p && delete!(gens, g.name)
        elseif g.fueltype =="CT"
            rand() > p && delete!(gens, g.name)
        elseif g.fueltype == "Diesel"
            rand() > p && delete!(gens, g.name)
        end
    end
    #Generators which cause needless numerical instability
    for g in values(gens)
        blocks, prices = getCurve(g)
        if minimum(blocks) < 1e-3
            println(g.name)
            delete!(gens, g.name)
        end
    end
end

newhr13Cap = sum([getCap(g, 13) for g in values(gens)])
scaling = newhr13Cap / hr13Cap
dts, vals = readLoads("../Data/ISO-NE Load Data/PredTest.csv")
dts_true, vals_true = readLoads("../Data/ISO-NE Load Data/LoadTest.csv")
vals *= scaling
vals_true *= scaling
resids = map(float, vals_true - vals);


D, V = eig(cov(resids));
numEigs = 1
projMatrix = V[:, (end-numEigs):end]'  #should thsi be scaled?
kappa(eps) = sqrt(1/eps - 1)

function testRun()
	rm2 = RobustModel(solver=GurobiSolver())
	alphas, uncs = createPolyUCS(rm2, resids, .1, .1, kappa(.1));
	aff = UCAff(rm2, gens, 1, uncs);
	aff.proj_fcn = (hr)-> projMatrix
	solve(aff, vals[1, :], report=true)
end

println( @elapsed testRun() )


