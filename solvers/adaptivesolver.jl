####
# Adaptive Solver
####
type UCAff
    m  #Robust Model
    gendata #handle to the generator dictionary
    penalty

    #::Dict{String, Array{Variable, 1}}
    ons 
    starts 
    stops 

    prod      #Dict{String, Array{FullAffExpr, 1}}   i.e. Aff Fcn of Unc.
    varcosts  #simple const for now since splitting everything 
    maxprod   #represents maximal production in that hour by that generator
    
    load_diffs##::Array{AffExpr, 1}  #VG onsider dropping
    sheds##Array{Variable, 1}  Simply variables for now, revisit

    uncs  #the uncertainties
    proj_fcn  # hr -> projMatrix

    #Ask Iain about this... seems poor form
    UCAff(robmodel, gens, penalty, uncs) =new(robmodel, gens, penalty, 
    		Dict{String, Any}(), Dict{String, Any}(), Dict{String, Any}(), 
    		Dict{String, Any}(), Dict{String, Any}(), Dict{String, Any}(),
            AffExpr[], Variable[], uncs, nothing)
end
#VG Add a better convenience constructor

solve_(aff::UCAff, report) = JuMPeR.solveRobust(aff.m, report=report)

function addSecondStage!(aff::UCAff; TOL=1e-8)
    for g in values(aff.gendata)   	
		#add the affine expressions for the production.  also ads ecomin-max constraints
        aff.prod[g.name] = [addAdaptProdVar!(aff.m, aff.proj_fcn(ihr), aff.uncs, aff.ons[g.name][ihr], 
        						g.ecomin[ihr], g.ecomax[ihr], TOL) for ihr=1:HRS]
		#populate maxprod
		#VG: Revisit if this should be affine
		@defVar(aff.m, maxprods[1:HRS])
		for ihr = 1:HRS
			addConstraint(aff.m, maxprods[ihr] >= aff.prod[g.name][ihr])
		end
		aff.maxprod[g.name] = maxprods[:]

		#This works because everything is disaggregated for now
        blocks, relcosts = getCurve(g)
        aff.varcosts[g.name] = [addVarCost!(aff.m, aff.maxprod[g.name][ihr], blocks, relcosts) for ihr = 1:HRS]
    end
end

# using JuMPeR
function addAdaptProdVar!(m, projMatrix, uncs, x, min, max, TOL)
	if max - min < TOL
		#A trivial case where adaptivity non-helpful
		return max * x
	end
	k = size(projMatrix, 1)
	@defVar(m, prods[1:k])
	@defVar(m, prods0)
	proj_uncs = projMatrix * uncs[:]
#	proj_uncs = [sum([projMatrix[ix, jx] * uncs[jx] for jx = 1:HRS]) for ix=1:k]
	p = sum([proj_uncs[ix] * prods[ix] for ix = 1:k]) + prods0 
    addConstraint(m, p <= max * x)
	addConstraint(m, p >= min * x)
	return p
end

function addLoadBalance!(aff::UCAff, forecasts, split_fcn = nothing)
	(split_fcn != nothing) && error("Splitting Loads not yet Implemented")
	uncs = aff.uncs
    for ihr = 1:length(forecasts)
        total_prod = sum([aff.prod[gname][ihr] for gname in keys(aff.gendata)])
        @defVar(aff.m, shed >= 0)
        addConstraint(aff.m, total_prod - forecasts[ihr] - uncs[ihr]  <= shed)
        addConstraint(aff.m, total_prod - forecasts[ihr] - uncs[ihr] >= -shed)
#        push!(aff.load_diffs, total_prod - forecasts[ihr])
        push!(aff.sheds, shed)
    end
end

