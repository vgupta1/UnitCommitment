###
# Helper Functions for building the Nominal UC Model
###
using JuMPeR, Gurobi, MathProgSolverInterface
#include("generators.jl")  #Speak to iain about file heirarchy
include("UCBase.jl")

##########################
function addLoadBalance!(nom::UCNom, forecasts)
    for ihr = 1:length(forecasts)
        total_prod = sum([nom.prod[gname][ihr] for gname in keys(nom.gendata)])
        @defVar(nom.m, shed >= 0)
        @addConstraint(nom.m, total_prod - forecasts[ihr] <= shed)
        @addConstraint(nom.m, total_prod - forecasts[ihr] >= -shed)

        push!(nom.load_diffs, total_prod - forecasts[ihr] )
        push!(nom.sheds, shed)
    end
end

solve_(nom::UCNom, report; args...) = JuMPeR.solve(nom.m)  #VG should probably extend somehow
