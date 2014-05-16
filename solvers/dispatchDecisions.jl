###
# Dispatch Decisions by fueltype
###
include("readGenerator.jl")
include("readLoads.jl")
include("nomsolver.jl")
include("robustsolver.jl")
include("UncSets.jl")
include("adaptivesolver.jl")

gens, scaling       = loadISO("../Data/AndysGenInstance", .1, true)
dts, vals           = readLoads("../Data/ISO-NE Load Data/PredTest.csv")
dts_true, vals_true = readLoads("../Data/ISO-NE Load Data/LoadTest.csv")
vals               *= scaling
vals_true          *= scaling
resids              = map(float, vals_true - vals);

dts, vals           = readLoads("../Data/ISO-NE Load Data/PredValidate.csv")
dts_true, vals_true = readLoads("../Data/ISO-NE Load Data/LoadValidate.csv")
vals               *= scaling
vals_true          *= scaling

kappa          = sqrt(1/.1 - 1)
penalty             = 5e3
Gamma1              = .5 * scaling
Gamma2              = 4  * scaling * scaling
numEigs             = 1

box_size            =.9333  #From 5-fold crossvalidation

indx                = 1

function cap_hr_fuel(ondict, gens, hr)
	fuel_dict = Dict{String, Float64}()
	for gname in keys(ondict)
		g = gens[gname]
		if !haskey(fuel_dict, g.fueltype)
			fuel_dict[g.fueltype] = 0
		end
		fuel_dict[g.fueltype] += getValue(ondict[gname][hr]) * getCap(g, hr)
	end
	fuel_dict
end

#Solve a Robust box problem
rm = RobustModel(solver=GurobiSolver(OutputFlag=0, MIPGap=1e-3, TimeLimit=60*15))
uncs = createUM(rm, resids, box_size)
rob    = UCRob(rm, gens, penalty, uncs)
status = solve(rob, vals[indx, :], usebox=false, report=false, prefer_cuts=true)

caprob12am = cap_hr_fuel(rob.ons, gens, 1)
println("\n Loadhr 12 am")
println(vals[indx, 1])
for k in keys(caprob12am)
	println(k, "  ", caprob12am[k])
end

caprob1pm = cap_hr_fuel(rob.ons, gens, 13)
println("\n Loadhr 1pm")
println(vals[indx, 13])
for k in keys(caprob1pm)
	println(k, "  ", caprob1pm[k])
end


#Solve the Affine UCS Problem
rm2 = RobustModel(solver=GurobiSolver(MIPGap=1e-3, OutputFlag=0, Method=3, TimeLimit=60*30), 
                 cutsolver=GurobiSolver(OutputFlag=0))
alphas, uncs = createPolyUCS(rm2, resids, Gamma1, Gamma2, kappa, true)
aff = UCAff(rm2, gens, penalty, uncs)
aff.proj_fcn = eigenProjMatrixData(resids, numEigs)
solve(aff, vals[indx, :], report=false, usebox=false, prefer_cuts=true)

caprob12am = cap_hr_fuel(aff.ons, gens, 1)
println("\n Loadhr 12 am")
println(vals[indx, 1])
for k in keys(caprob12am)
	println(k, "  ", caprob12am[k])
end

caprob1pm = cap_hr_fuel(aff.ons, gens, 13)
println("\n Loadhr 1pm")
println(vals[indx, 13])
for k in keys(caprob1pm)
	println(k, "  ", caprob1pm[k])
end

