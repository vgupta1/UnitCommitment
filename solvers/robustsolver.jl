###
# Helper Functions for building the Nominal UC Model
###
using JuMPeR, Gurobi

#VG 
#There are a non-trivial number where 
#ecomin = ecomax by hour for all hours (or some hours)
#Some even where Ecomax = 0 for all hours (or some hours)
#create this subclass within secondStage fcns
#for gen-hrs where eco-min == eco-max, 
#fix the production level (constant, non-adaptive)

#Move this to the ucbase structure instead
type UCRob
    m  #Robust Model
    gendata #handle to the generator dictionary
    penalty

    #::Dict{String, Array{Variable, 1}}
    ons 
    starts 
    stops 
    prod      #(nonadaptive) robust model
    varcosts  
    maxprod   #represents maximal production in that hour by that generator
    
    load_diffs##::Array{AffExpr, 1}  #doesn't include unc
    sheds##Array{Variable, 1}

    #Ask Iain about this... seems poor form
    UCRob(model, gens, penalty) =new(model, gens, penalty, Dict{String, Any}(), 
            Dict{String, Any}(), Dict{String, Any}(), Dict{String, Any}(), Dict{String, Any}(),
            Dict{String, Any}(), AffExpr[], Variable[])
end
UCRob(gens; PENALTY=5000.) = UCRob(RobustModel(solver=GurobiSolver()), gens, PENALTY)

#add uncertainty set constraints to the model
#returns the bounding box parameters and the uncertainties
#data stored row-wise
function createBertSimU(m, resid_data, Gamma; Gamma_bound=2.0)
	N, d = size(resid_data)
	mu = mean(resid_data, 1); sigma = std(resid_data, 1)
	@defUnc(m, u[1:d] >= 0)
	@defUnc(m, z[1:d])
	for ix = 1:d
		addConstraint(m, u[ix] <= mu[ix] + Gamma_bound*sigma[ix])
		addConstraint(m, u[ix] >= mu[ix] - Gamma_bound*sigma[ix])
		addConstraint(m, (u[ix] - mu[ix])/sigma[ix] <= z[ix])
		addConstraint(m, (u[ix] - mu[ix])/sigma[ix] >= -z[ix])
	end
	addConstraint(m, sum(z[:]) <= Gamma)

	#now compute the bounding box
	minG = min(Gamma, Gamma_bound)
	return [(mu[ix] + minG * sigma[ix])::Float64 for ix = 1:d], u
end

## addSecondStage! of the nominal model pertains here bc non-adaptive

function addLoadBalance!(rob::UCRob, forecasts, uncs)
    for ihr = 1:length(forecasts)
        total_prod = sum([rob.prod[gname][ihr] for gname in keys(rob.gendata)])
        @defVar(rob.m, shed >= 0)
        addConstraint(rob.m, total_prod - forecasts[ihr] - uncs[ihr]  <= shed)
        addConstraint(rob.m, total_prod - forecasts[ihr] - uncs[ihr] >= -shed)

        push!(rob.load_diffs, total_prod - forecasts[ihr])
        push!(rob.sheds, shed)
    end
end

#VG Try to factor this some more to reduce overhead
function solve(rob::UCRob, forecasts, uncs; forceserve=false)
   addFirstStage!(rob)
   addSecondStage!(rob)
   addLoadBalance!(rob, forecasts, uncs)
   @setObjective(rob.m, Min, sum{g.startcost*rob.starts[g.name][ihr], g in values(rob.gendata), ihr=1:HRS} + 
                        + sum{rob.varcosts[g.name][ihr], g in values(rob.gendata), ihr=1:HRS} + 
                        rob.penalty * sum{rob.sheds[ihr], ihr=1:HRS})

   print

    if forceserve
        for ix =1:HRS
            @addConstraint(rob.m, rob.sheds[ix] == 0)
        end
    end

   #Load the model only for timing
   status= JuMPeR.solveRobust(rob.m, report=true)

   println("Start Costs\t $(getStartCost(rob))")
   println("Objective\t $(getObjectiveValue(rob.m))")
   return status
end



