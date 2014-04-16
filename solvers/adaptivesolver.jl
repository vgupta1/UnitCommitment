####
# Adaptive Solver
####
using JuMPeR

type WarmStartInfo
	#::Dict{String, Vector{Int}}
	ons
	starts
	stops

	#::Dict{String, Vector{Float64}}
	prod  
	varcosts 

	#::Vector{Float64}
	sheds
	WarmStartInfo() = new(Dict{String, Vector{Float64}}(), Dict{String, Vector{Float64}}(), Dict{String, Vector{Float64}}(), 
						  Dict{String, Vector{Float64}}(), Dict{String, Vector{Float64}}(), Float64[])
end

function copyWarmStart(ucbase, w)
	for gn in keys(ucbase.gendata)
		w.ons[gn]      = Int[] 
		w.starts[gn]   = Int[] 
		w.stops[gn]    = Int[]
		w.prod[gn]     = Float64[] 
		w.varcosts[gn] = Float64[]
		append!(w.ons[gn],      map(getValue, ucbase.ons[gn])[:])
		append!(w.starts[gn],   map(getValue, ucbase.starts[gn])[:])
		append!(w.stops[gn],    map(getValue, ucbase.stops[gn])[:])
		append!(w.prod[gn],     map(getValue, ucbase.prod[gn])[:])
		append!(w.varcosts[gn], map(getValue, ucbase.varcosts[gn])[:])
	end
	append!(w.sheds, map(getValue, ucbase.sheds))
	w.sheds += 1e-4  #to deal with numeical stability stuff
end

type UCAff
    m  #::RobustModel
    gendata 
    penalty

    #::Dict{String, Array{Variable, 1}}
    ons 
    starts 
    stops 

    prod      #Dict{String, Vector{FullAffExpr}}   i.e. Aff Fcn of Unc.
    varcosts  #Vector{Variable} since we are disaggregating everything
    maxprod   #Vector{Variable} maximal production in that hour by that generator
    
    #actual affine policy
    prodcoeff
    prodint

    sheds ##Vector{Variable}  Simply variables for now, revisit
    uncs  
    proj_fcn  # hr -> projMatrix

    warmstart

    UCAff(robmodel, gens, penalty, uncs) =new(robmodel, gens, penalty, 
    		Dict{String, Any}(), Dict{String, Any}(), Dict{String, Any}(), #ons starts stops 
    		Dict{String, Any}(), Dict{String, Any}(), Dict{String, Any}(),  #prod, varcost, maxprod
            Dict{String, Any}(), Dict{String, Any}(), #prodcoeff, prodintercept
            Variable[], uncs, nothing, nothing)  #sheds, uncs, proj_fcn, warmstart
end
#VG Add a better convenience constructor

solve_(aff::UCAff, report; args...) = JuMPeR.solveRobust(aff.m, report=report; args...)
isBaseLoad(g) = g.fueltype == "Nuclear"

function addSecondStage!(aff::UCAff; TOL=1e-8)
    for g in values(aff.gendata)
		gn = g.name
		if isBaseLoad(g)
			ps = [addProdVar!(aff.m, aff.ons[gn][ihr], g.ecomin[ihr], g.ecomax[ihr], TOL) for ihr=1:HRS]
			aff.prodcoeff[gn]    = [ [] for ihr = 1:HRS ]
			aff.prodint[gn] = ps
			aff.prod[gn]          = ps
			aff.maxprod[gn]       = ps
		else
	    	jumble = [addAdaptProdVar!(aff.m, aff.proj_fcn(ihr), aff.uncs, aff.ons[gn][ihr], 
        								g.ecomin[ihr], g.ecomax[ihr], TOL) for ihr=1:HRS]
    		aff.prodcoeff[gn], aff.prodint[gn], aff.prod[gn], aff.maxprod[gn] = zip(jumble...)
    	end

		#This works because everything is disaggregated for now
        blocks, relcosts = getCurve(g)
        aff.varcosts[gn] = [addVarCost!(aff.m, aff.maxprod[gn][ihr], blocks, relcosts) for ihr = 1:HRS]
    end
end

function addWarmStart!(aff::UCAff; force=false)
	if force
		return forceWarmStart!(aff)
	end
	w = aff.warmstart
	if w == nothing
		return
	end
	for gn in keys(aff.gendata)
		for ihr = 1:HRS
			# ons, starts and stops
			setValue(aff.ons[gn][ihr],    w.ons[gn][ihr])
			setValue(aff.starts[gn][ihr], w.starts[gn][ihr])
			setValue(aff.stops[gn][ihr],  w.stops[gn][ihr])
			map(x ->setValue(x, 0.0), aff.prodcoeff[gn][ihr])
			if typeof(aff.prodint[gn][ihr]) == Variable
				setValue(aff.prodint[gn][ihr], w.prod[gn][ihr])
			end
			setValue(aff.varcosts[gn][ihr], w.varcosts[gn][ihr])
		end
	end
	for ihr = 1:HRS
		setValue(aff.sheds[ihr], w.sheds[ihr])
	end
end

function forceWarmStart!(aff::UCAff)
	w = aff.warmstart
	if w == nothing
		return
	end
	for gname in keys(aff.gendata)
		for ihr = 1:HRS
			# ons, starts and stops
			@addConstraint(aff.m, aff.ons[gname][ihr]    == int(w.ons[gname][ihr]))
			@addConstraint(aff.m, aff.starts[gname][ihr] == int(w.starts[gname][ihr]))
			@addConstraint(aff.m, aff.stops[gname][ihr]  == int(w.stops[gname][ihr]))

			# prods and varcosts
			for ix = 1:length(aff.prodcoeff[gname][ihr])
				@addConstraint(aff.m, aff.prodcoeff[gname][ihr][ix] == 0.0 )
			end
			if aff.prodint[gname][ihr] != nothing
				@addConstraint(aff.m, aff.prodint[gname][ihr] == w.prod[gname][ihr])
			end
			@addConstraint(aff.m, aff.varcosts[gname][ihr] == w.varcosts[gname][ihr])
		end
	end
	for ihr = 1:HRS
		@addConstraint(aff.m, aff.sheds[ihr] == w.sheds[ihr])
	end
end

function addAdaptProdVar!(m, projMatrix, uncs, x, min, max, TOL)
	if max - min < TOL
		#A trivial case where adaptivity non-helpful
		return [], nothing, max * x, max * x
	end
	k = size(projMatrix, 1)
	@defVar(m, prods[1:k])
	@defVar(m, prods0 >= 0)
	@defVar(m, maxprod >=0)
	proj_uncs = projMatrix * uncs[:]
	p = sum([proj_uncs[ix] * prods[ix] for ix = 1:k]) + prods0 
    
    addConstraint(m, p <= maxprod)
    @addConstraint(m, maxprod <= max * x)
	addConstraint(m, p >= min * x)
	return prods, prods0, p, maxprod
end

function addLoadBalance!(aff::UCAff, forecasts, split_fcn = nothing)
	(split_fcn != nothing) && error("Splitting Loads not yet Implemented")
	uncs = aff.uncs
    for ihr = 1:length(forecasts)
        total_prod = sum([aff.prod[gname][ihr] for gname in keys(aff.gendata)])
        @defVar(aff.m, shed >= 0)
        addConstraint(aff.m, total_prod - forecasts[ihr] - uncs[ihr]  <= shed)
        addConstraint(aff.m, total_prod - forecasts[ihr] - uncs[ihr] >= -shed)
        push!(aff.sheds, shed)
    end
    addWarmStart!(aff)  ##VG This should live somewhere else...
end

