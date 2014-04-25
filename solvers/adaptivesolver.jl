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

function copyWarmStart(ucbase, w; TOL=1e-4)
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
	w.sheds += TOL  #to deal with numeical stability stuff
	return w
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
    proj_fcn  #hr -> projMatrix

	#may or may not be set by user
    warmstart  
    sample_uncs

    UCAff(robmodel, gens, penalty, uncs) =new(robmodel, gens, penalty, 
    		Dict{String, Any}(), Dict{String, Any}(), Dict{String, Any}(), #ons starts stops 
    		Dict{String, Any}(), Dict{String, Any}(), Dict{String, Any}(),  #prod, varcost, maxprod
            Dict{String, Any}(), Dict{String, Any}(), #prodcoeff, prodintercept
            Variable[], uncs, nothing, nothing, nothing)  #sheds, uncs, proj_fcn, warmstart
end
#VG Add a better convenience constructor

solve_(aff::UCAff, report; args...) = JuMPeR.solveRobust(aff.m, report=report; args...)
isBaseLoad(g) = g.fueltype == "Nuclear"
addSampleUncs!(aff::UCAff, samples) = aff.sample_uncs = copy(samples)

function addWarmStart!(aff::UCAff; force=false)
	w = aff.warmstart
	if w == nothing
		return
	end
	if force
		return forceWarmStart!(aff)
	end

	for gn in keys(aff.gendata)
		for ihr = 1:HRS
			# # ons, starts and stops
			setValue(aff.ons[gn][ihr],    w.ons[gn][ihr])
			setValue(aff.starts[gn][ihr], w.starts[gn][ihr])
			setValue(aff.stops[gn][ihr],  w.stops[gn][ihr])
			map(x ->setValue(x, 0.0), aff.prodcoeff[gn][ihr])  #if empty, do nothing
			setValue(aff.varcosts[gn][ihr], w.varcosts[gn][ihr])

			#VG This sort of embedded logic about baseload and trivial elements is yucky
			if isa(aff.prodint[gn][ihr], Variable)
				setValue(aff.prodint[gn][ihr], w.prod[gn][ihr])
			elseif isa(aff.prod[gn][ihr], Variable)
				setValue(aff.prod[gn][ihr], w.prod[gn][ihr])
			end
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
	for gn in keys(aff.gendata)
		for ihr = 1:HRS
			# ons, starts and stops
			@addConstraint(aff.m, aff.ons[gn][ihr]    == int(w.ons[gn][ihr]))
			@addConstraint(aff.m, aff.starts[gn][ihr] == int(w.starts[gn][ihr]))
			@addConstraint(aff.m, aff.stops[gn][ihr]  == int(w.stops[gn][ihr]))
			@addConstraint(aff.m, aff.varcosts[gn][ihr] == w.varcosts[gn][ihr])

			# prods
			for ix = 1:length(aff.prodcoeff[gn][ihr])
				@addConstraint(aff.m, aff.prodcoeff[gn][ihr][ix] == 0.0 )
			end

			#VG This sort of embedded logic about baseload and trivial elements is yucky
			if isa(aff.prodint[gn][ihr], Variable)
				@addConstraint(aff.m, aff.prodint[gn][ihr] ==  w.prod[gn][ihr])
			elseif isa(aff.prod[gn][ihr], Variable)
				@addConstraint(aff.m, aff.prod[gn][ihr]==  w.prod[gn][ihr])
			end
		end
	end
	for ihr = 1:HRS
		@addConstraint(aff.m, aff.sheds[ihr] == w.sheds[ihr])
	end
end

function addSecondStage!(aff::UCAff; TOL=1e-8)
    for g in values(aff.gendata)
		gn = g.name
		if isBaseLoad(g)
			ps = [addProdVar!(aff.m, aff.ons[gn][ihr], g.ecomin[ihr], g.ecomax[ihr], TOL) for ihr=1:HRS]
			aff.prodcoeff[gn]     = [ [] for ihr = 1:HRS ]
			aff.prodint[gn]       = [nothing for ihr = 1:HRS]
			aff.prod[gn]          = [1 * ps[ihr] for ihr = 1:HRS] #make it an affine expression

			#for type consistency, make maxprods a variable
			@defVar(aff.m, maxprods[1:HRS] >= 0)
			for ihr= 1:HRS
				@addConstraint(aff.m, ps[ihr] <= maxprods[ihr])
			end
			aff.maxprod[gn] = maxprods
		else
	    	jumble = [addAdaptProdVar!(aff.m, aff.proj_fcn(ihr), aff.uncs, aff.ons[gn][ihr], 
        								g.ecomin[ihr], g.ecomax[ihr], TOL) for ihr=1:HRS]
    		aff.prodcoeff[gn], aff.prodint[gn], aff.prod[gn], aff.maxprod[gn] = zip(jumble...)
    		
    		#if there are samples, add them now
    		if aff.sample_uncs != nothing
	    		for indx = 1:size(aff.sample_uncs, 1)
	    			for ihr = 1:HRS
	    				addEcoMinMax!(aff.m, aff.proj_fcn(ihr), aff.sample_uncs[indx, :], 
	    								aff.ons[gn][ihr], g.ecomin[ihr], aff.maxprod[gn][ihr], 
	    								aff.prodcoeff[gn][ihr], aff.prodint[gn][ihr])
	    			end
	    		end
	    	end
    	end

		#This works because everything is disaggregated for now
        blocks, relcosts = getCurve(g)
        aff.varcosts[gn] = [addVarCost!(aff.m, aff.maxprod[gn][ihr], blocks, relcosts) for ihr = 1:HRS]
    end
end

function addEcoMinMax!(m, projMatrix, uncs, x, min, maxprod, prodcoeff, prods0)
	k = length(prodcoeff)  
	if k == 0  #Not affine adaptive
		return prods0
	else
		proj_uncs = projMatrix * uncs[:]
		# println("Sizes: ", size(projMatrix), " ", size(uncs[:]), " ", size(proj_uncs[:]))
		prod = sum([proj_uncs[ix] * prodcoeff[ix] for ix = 1:k]) + prods0 
	    addConstraint(m, prod - maxprod <= 0)
		addConstraint(m, prod - min * x >= 0)
		return prod
	end
end

function addAdaptProdVar!(m, projMatrix, uncs, x, min, max, TOL)
	@defVar(m, maxprod >=0)
    @addConstraint(m, maxprod <= max * x)

	if max - min < TOL
		#A trivial case where adaptivity non-helpful
		return [], nothing, max * x, maxprod
	end
	k = size(projMatrix, 1)
	@defVar(m, prodcoeff[1:k])
	@defVar(m, prodint >= 0)
	prod = addEcoMinMax!(m, projMatrix, uncs, x, min, maxprod, prodcoeff, prodint)
	return prodcoeff, prodint, prod, maxprod
end

function addLoadBalance!(aff::UCAff, forecasts)
	uncs = aff.uncs
    for ihr = 1:length(forecasts)
        total_prod = sum([aff.prod[gn][ihr] for gn in keys(aff.gendata)])
        @defVar(aff.m, shed >= 0)
        addConstraint(aff.m, total_prod - shed - uncs[ihr]  <= forecasts[ihr])
        addConstraint(aff.m, total_prod + shed - uncs[ihr]  >= forecasts[ihr])
        push!(aff.sheds, shed)

        # Add the sampled points, if they exist
        # VG Factor this out somehow
        if aff.sample_uncs != nothing
	        projMatrix = aff.proj_fcn(ihr)
	    	for ix = 1:size(aff.sample_uncs, 1)
        		# form the value of the affine value here
        		total_prod = AffExpr()
        		proj_uncs = projMatrix * aff.sample_uncs[ix, :]'
        		for gn in keys(aff.gendata)
	    			for jx = 1:length(aff.prodcoeff[gn][ihr])  #if non-affine, empty loop
	    				push!(total_prod, proj_uncs[jx], aff.prodcoeff[gn][ihr][jx])
	    			end
	    			#VG Yucky logic
	    			if aff.prodint[gn][ihr] != nothing 
	    				push!(total_prod, 1., aff.prodint[gn][ihr])
	    			elseif isa(aff.prod[gn][ihr], Variable) #robust policy
	    				push!(total_prod, 1., aff.prod[gn][ihr])
	    			else  #a trivial min=max case
	    				append!(total_prod, aff.prod[gn][ihr]) 
	    			end
        		end
        		# add the sampled cnsts
		        @addConstraint(aff.m, total_prod - shed - aff.sample_uncs[ix, ihr]  <= forecasts[ihr])
		        @addConstraint(aff.m, total_prod + shed - aff.sample_uncs[ix, ihr]  >= forecasts[ihr])
        	end
        end
    end
    addWarmStart!(aff)  ##VG Move this to main with a try-catch
end

