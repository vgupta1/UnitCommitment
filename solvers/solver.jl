###
# Helper Functions for building the UC Models
###
using JuMP, Gurobi, MathProgSolverInterface
#include("generators.jl")

#VG 
#There are a non-trivial number where 
#ecomin = ecomax by hour for all hours (or some hours)
#Some even where Ecomax = 0 for all hours (or some hours)
#create this subclass within secondStage fcns
#for gen-hrs where eco-min == eco-max, 
#fix the production level (constant, non-adaptive)

######
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
UCNom(gens; PENALTY=5000.) = UCNom(Model(solver=GurobiSolver()), gens, PENALTY)

#Assumes already solved
getStartCost(nom, gen, hr) = gen.startcost * getValue(nom.starts[gen.name][hr])
getStartCost(nom, hr) = sum([getStartCost(nom, g, hr) for g in values(nom.gendata)])
getStartCost(nom) = sum([getStartCost(nom, hr) for hr in 1:HRS])
getVarCost(nom, gen, hr) = getValue(nom.varcosts[gen.name][hr])
getVarCost(nom, hr) = sum([ getVarCost(nom, g, hr) for g in values(nom.gendata)])
getVarCost(nom) = sum([getVarCost(nom, hr) for hr in 1:HRS])
getgap(nom) = (t=getInternalModel(nom.m); (getobjval(t) - getobjbound(t))/getobjval(t))




#The whole shebang
function solve(nom::UCNom, forecasts; forceserve=false)
   tic();
   addFirstStage!(nom)
   addSecondStage!(nom)
   addLoadBalance!(nom, forecasts)
   @setObjective(nom.m, Min, sum{g.startcost*nom.starts[g.name][ihr], g in values(nom.gendata), ihr=1:HRS} + 
                        + sum{nom.varcosts[g.name][ihr], g in values(nom.gendata), ihr=1:HRS} + 
                        nom.penalty * sum{nom.sheds[ihr], ihr=1:HRS})

    if forceserve
        for ix =1:HRS
            @addConstraint(nom.m, nom.sheds[ix] == 0)
        end
    end

   #Load the model only for timing
   JuMP.solve(nom.m, load_model_only=true)
   buildTime = toc();
   tic(); status = JuMP.solve(nom.m); solveTime = toc();

   println("Build Time \t $buildTime")
   println("Solve Time \t $solveTime")
   println("Start Costs\t $(getStartCost(nom))")
   println("Objective\t $(getObjectiveValue(nom.m))")
   return status
end

function secondSolve(nom, loads)
    tic()
    nom2 = UCNom(nom.gendata, PENALTY = nom.penalty)

    #copy over first stage solution
    for g in values(nom2.gendata)
        nom2.ons[g.name] = map(v->getValue(v), nom.ons[g.name][:])
        nom2.starts[g.name] = map(v->getValue(v), nom.starts[g.name][:])
        nom2.stops[g.name] = map(v->getValue(v), nom.stops[g.name][:])
    end
    addSecondStage!(nom2)  
    addLoadBalance!(nom2, loads)
    @setObjective(nom2.m, Min, sum{g.startcost*nom2.starts[g.name][ihr], g in values(nom2.gendata), ihr=1:HRS} + 
                        + sum{nom2.varcosts[g.name][ihr], g in values(nom2.gendata), ihr=1:HRS} + 
                        nom2.penalty * sum{nom2.sheds[ihr], ihr=1:HRS})

   JuMP.solve(nom2.m, load_model_only=true)
   buildTime = toc()
   tic(); status = JuMP.solve(nom2.m); solveTime = toc()

   println("Build Time \t $buildTime")
   println("Solve Time \t $solveTime")
   println("Start Costs\t $(getStartCost(nom))")
   println("Objective\t $(getObjectiveValue(nom.m))")
   return nom2
end

#adds first stage variables and constraints for all generators 
function addFirstStage!(nom; TOL=1e-8)   
    for g in values(nom.gendata)
        #omit useless things from the formulation
        if maximum(g.ecomax) <= TOL
            delete!(nom.gendata, g.name)
            warn("Deleting: $(g.name)")
            continue
        end

        ons, starts, stops = addFirstStage(nom.m)
        addMinUp(nom.m, ons, starts, g.minup)
        addMinDown(nom.m, ons, stops, g.mindown)
        if g.maxstart <= HRS
            @addConstraint(nom.m, sum{starts[ix], ix=1:length(starts)} <= g.maxstart)  
        end          

        #Constrain it to be off in meaningless hours
        for ihr = 1:HRS
            if g.ecomax[ihr] <= TOL
                @addConstraint(nom.m, ons[ihr] == 0)
            end
        end
        nom.ons[g.name] = ons; nom.starts[g.name] = starts; nom.stops[g.name] = stops
    end
end

function addSecondStage!(nom; TOL=1e-8)
    for g in values(nom.gendata)
        nom.prod[g.name] = [addProdVar!(nom.m, nom.ons[g.name][ihr], g.ecomin[ihr], g.ecomax[ihr], TOL) for ihr=1:HRS]
        blocks, relcosts = getCurve(g)
        nom.varcosts[g.name] = [addVarCost!(nom.m, nom.prod[g.name][ihr], blocks, relcosts) for ihr = 1:HRS]
    end
end

function addLoadBalance!(nom, forecasts)
    for ihr = 1:length(forecasts)
        total_prod = sum([nom.prod[gname][ihr] for gname in keys(nom.gendata)])
        @defVar(nom.m, shed >= 0)
        @addConstraint(nom.m, total_prod - forecasts[ihr] <= shed)
        @addConstraint(nom.m, total_prod - forecasts[ihr] >= -shed)

        push!(nom.load_diffs, total_prod - forecasts[ihr] )
        push!(nom.sheds, shed)
    end
end

#Adds on, starting, stopping variables for a single generator
function addFirstStage(m; HRS=24)
    @defVar(m, xs[1:HRS], Bin)  #on
    @defVar(m, zs[1:HRS], Bin)  #starting
    @defVar(m, ys[1:HRS], Bin)  #stopping
    
    for ix = 1:HRS
        @addConstraint(m, zs[ix] <= xs[ix])  #Start now -> On
        @addConstraint(m, ys[ix] <= 1-xs[ix])  #Stop now -> Off
    end
    
    #assume everyone starts off
    @addConstraint(m, xs[1] <= zs[1])
    @addConstraint(m, ys[1] == 0)
    
    for ix = 2:HRS
        @addConstraint(m, zs[ix] <= 1-xs[ix-1]) #Start now -> Was off
        @addConstraint(m, ys[ix] <= xs[ix-1])  #Stop now -> Was on
    
        @addConstraint(m, xs[ix] - xs[ix-1] <= zs[ix])  #Was off, now on -> Start
        @addConstraint(m, xs[ix-1] - xs[ix] <= ys[ix])  #Was on, now off -> Start
    end
    return xs, zs, ys
end

function addMinUp(m, xs, zs, minUp::Int)
    #Neglects ending effects
    #Interpret minUp <= 1 to be redundant.
    if minUp <= 1
        return
    end
    hrs = length(xs)
    for ix = 1:hrs
        for jx = ix+1:min(hrs, ix + minUp - 1)  #ix counts as 1 hour by itself
            @addConstraint(m, zs[ix] <= xs[jx])
        end
    end
end

function addMinDown(m, xs, ys, minDown::Int)
    if minDown <= 1
        return
    end
    hrs = length(xs)
    for ix = 1:hrs
        for jx = ix+1:min(hrs, ix + minDown - 1)  #ix counts as 1 hour by itself
            @addConstraint(m, ys[ix] <= 1 - xs[jx])
        end
    end
end

#creates the variable and adds teh ecomin/ecomax constraint
function addProdVar!(m, x, min, max, TOL)
    assert(max >= min)
    if max - min <= TOL  #If trivial, return expression
        return max * x
    else
        @defVar(m, p)
        @addConstraint(m, p <= max * x)
        @addConstraint(m, p >= min * x)
        return p
    end
end

#could be smarter to take advantage of the fixed case
function addVarCost!(m, prod, blocks, relcosts)
    numBlocks = length(blocks)
    #formulation exploits the fact that relcosts are non-decreasing
    @defVar(m, 0 <= deltas[1:numBlocks] <= 1)
    @addConstraint(m, sum{blocks[ix] * deltas[ix], ix=1:numBlocks} == prod)
    @defVar(m, varCost >= 0)
    @addConstraint(m, sum{blocks[ix] * relcosts[ix] * deltas[ix], ix=1:numBlocks} == varCost)
    return varCost
end

