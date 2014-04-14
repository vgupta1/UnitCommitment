####
# The UC Nom Object
####
# VG: Refactor this into a data object with its own access methods.  Belongs in UCBase
## Then place all UCNom functionalty in the nomsolver file.
type UCNom
    m  #model
    gendata #handle to the generator dictionary
    penalty

    #::Dict{String, Array{Variable, 1}}
    ons 
    starts 
    stops 
    prod 
    varcosts  #cost of generation
    
    ##::Array{AffExpr, 1}
    load_diffs
    sheds

    #Ask Iain about this... seems poor form
    UCNom(model, gens, penalty) =new(model, gens, penalty, Dict{String, Any}(), 
            Dict{String, Any}(), Dict{String, Any}(), Dict{String, Any}(), 
            Dict{String, Any}(), AffExpr[], Variable[])
end
UCNom(gens; PENALTY=5000.) = UCNom(RobustModel(solver=GurobiSolver(OutputFlag=0)), gens, PENALTY)
