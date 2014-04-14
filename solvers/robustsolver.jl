###
# Helper Functions for building the Nominal UC Model
###
using JuMPeR, Gurobi

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
    
    load_diffs##::Array{AffExpr, 1}  #doesn't include unc
    sheds##Array{Variable, 1}

    uncs  #the uncertainties

    #Ask Iain about this... seems poor form
    UCRob(robmodel, gens, penalty, uncs) =new(robmodel, gens, penalty, Dict{String, Any}(), 
            Dict{String, Any}(), Dict{String, Any}(), Dict{String, Any}(), Dict{String, Any}(),
            AffExpr[], Variable[], uncs)
end
#VG Add a better convenience constructor
#Current interace requires user to understand the uncertainty set on their own

function addLoadBalance!(rob::UCRob, forecasts)
	uncs = rob.uncs
    for ihr = 1:length(forecasts)
        total_prod = sum([rob.prod[gname][ihr] for gname in keys(rob.gendata)])
        @defVar(rob.m, shed >= 0)
        addConstraint(rob.m, total_prod - forecasts[ihr] - uncs[ihr]  <= shed)
        addConstraint(rob.m, total_prod - forecasts[ihr] - uncs[ihr] >= -shed)

        push!(rob.load_diffs, total_prod - forecasts[ihr])
        push!(rob.sheds, shed)
    end
end

## addSecondStage! of the nominal model pertains here bc non-adaptive
solve_(rob::UCRob, report) = JuMPeR.solveRobust(rob.m, report=report)
