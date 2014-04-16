### 
# UC Base
###
# Contains generic functionality which should pertain to all UC Models (nominal robust, affine)
# VG Introduce a subdata type to make this more obvious/consistent
include("UCNom.jl")  #VG Fix Heirarchy

#Assumes already solved  VG: Add a check
getStartCost(ucbase, gen, hr) = gen.startcost * getValue(ucbase.starts[gen.name][hr])
getStartCost(ucbase, hr) = sum([getStartCost(ucbase, g, hr) for g in values(ucbase.gendata)])
getStartCost(ucbase) = sum([getStartCost(ucbase, hr) for hr in 1:HRS])
getVarCost(ucbase, gen, hr) = getValue(ucbase.varcosts[gen.name][hr])
getVarCost(ucbase, hr) = sum([ getVarCost(ucbase, g, hr) for g in values(ucbase.gendata)])
getVarCost(ucbase) = sum([getVarCost(ucbase, hr) for hr in 1:HRS])
getgap(ucbase) = (t=getInternalModel(ucbase.m); (getobjval(t) - getobjbound(t))/getobjval(t))
getCap(ucbase, hr) = sum([getCap(g, hr) * getValue(ucbase.ons[g.name][hr]) for g in values(ucbase.gendata)])
getCap(ucbase) = sum([getCap(ucbase, hr) for hr = 1:HRS])
totShed(ucbase) = sum(map(getValue, ucbase.sheds))

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

function secondSolve(ucbase, loads; report=true)
	tic()
	nom2 = UCNom(ucbase.gendata, PENALTY = ucbase.penalty)

	#copy over first stage solution
	for g in values(nom2.gendata)
	    nom2.ons[g.name] = map(v->getValue(v), ucbase.ons[g.name][:])
	    nom2.starts[g.name] = map(v->getValue(v), ucbase.starts[g.name][:])
	    nom2.stops[g.name] = map(v->getValue(v), ucbase.stops[g.name][:])
	end
	addSecondStage!(nom2)  
	addLoadBalance!(nom2, loads)
	@setObjective(nom2.m, Min, sum{g.startcost*nom2.starts[g.name][ihr], g in values(nom2.gendata), ihr=1:HRS} + 
	                    + sum{nom2.varcosts[g.name][ihr], g in values(nom2.gendata), ihr=1:HRS} + 
	                    nom2.penalty * sum{nom2.sheds[ihr], ihr=1:HRS})
	JuMPeR.solve(nom2.m, load_model_only=true)
	buildTime = toq()
	tic(); status = JuMPeR.solve(nom2.m); solveTime = toq()

	if report
		println("Build Time \t $buildTime")
		println("Solve Time \t $solveTime")
		println("Start Costs\t $(getStartCost(ucbase))")
		println("Capacity \t $(getCap(ucbase))")
		println("Second Stage\t $(getObjectiveValue(nom2.m) - getStartCost(ucbase))")
		println("Tot Shed\t $(totShed(nom2))")
        println("Objective\t $(getObjectiveValue(nom2.m))")
	end
	return nom2
end

#adds first stage variables and constraints for all generators 
function addFirstStage!(ucbase; TOL=1e-8)   
    for g in values(ucbase.gendata)
        #omit useless things from the formulation
        if maximum(g.ecomax) <= TOL
            delete!(ucbase.gendata, g.name)
            warn("Deleting: $(g.name)")
            continue
        end

        ons, starts, stops = addFirstStage(ucbase.m)
        addMinUp(ucbase.m, ons, starts, g.minup)
        addMinDown(ucbase.m, ons, stops, g.mindown)
        if g.maxstart <= HRS
            @addConstraint(ucbase.m, sum{starts[ix], ix=1:length(starts)} <= g.maxstart)  
        end          

        #Constrain it to be off in meaningless hours
        for ihr = 1:HRS
            if g.ecomax[ihr] <= TOL
                @addConstraint(ucbase.m, ons[ihr] == 0)
            end
        end
        ucbase.ons[g.name] = ons; ucbase.starts[g.name] = starts; ucbase.stops[g.name] = stops
    end
end

#The whole shebang
function solve(ucbase, forecasts; forceserve=false, report=true, args...)
	tic();
	addFirstStage!(ucbase)
	addSecondStage!(ucbase)
	addLoadBalance!(ucbase, forecasts)
	@setObjective(ucbase.m, Min, sum{g.startcost*ucbase.starts[g.name][ihr], g in values(ucbase.gendata), ihr=1:HRS} + 
                        + sum{ucbase.varcosts[g.name][ihr], g in values(ucbase.gendata), ihr=1:HRS} + 
                        ucbase.penalty * sum{ucbase.sheds[ihr], ihr=1:HRS})

    if forceserve
        for ix =1:HRS
            @addConstraint(ucbase.m, ucbase.sheds[ix] == 0)
        end
    end
    buildTime = toq();

	tic(); status = solve_(ucbase, report; args...); solveTime = toq();

	if report 
		println("Build Time \t $buildTime")
		println("Solve Time \t $solveTime")
		println("Start Costs\t $(getStartCost(ucbase))")
		println("Capacity \t $(getCap(ucbase))")
		println("Second Stage\t $(getObjectiveValue(ucbase.m) - getStartCost(ucbase))")
		println("Tot Shed \t $(totShed(ucbase))")
		println("Objective\t $(getObjectiveValue(ucbase.m))")
	end
	return status
end

#used both by nominal and robust
function addSecondStage!(nom; TOL=1e-8)
    for g in values(nom.gendata)
        nom.prod[g.name] = [addProdVar!(nom.m, nom.ons[g.name][ihr], g.ecomin[ihr], g.ecomax[ihr], TOL) for ihr=1:HRS]
        blocks, relcosts = getCurve(g)
        nom.varcosts[g.name] = [addVarCost!(nom.m, nom.prod[g.name][ihr], blocks, relcosts) for ihr = 1:HRS]
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

#could be smarter to take advantage of the fixed production case
function addVarCost!(m, prod, blocks, relcosts)
    numBlocks = length(blocks)
    #formulation exploits the fact that relcosts are non-decreasing
    @defVar(m, 0 <= deltas[1:numBlocks] <= 1)
    @addConstraint(m, sum{blocks[ix] * deltas[ix], ix=1:numBlocks} == prod)
    @defVar(m, varCost >= 0)
    @addConstraint(m, sum{blocks[ix] * relcosts[ix] * deltas[ix], ix=1:numBlocks} == varCost)
    return varCost
end

