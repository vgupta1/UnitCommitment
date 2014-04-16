#####
#  Does all the reading in from the amply style files from andy's gen instance
#####
include("generators.jl")

#Notes:
# Ignoring NoLoadCost  - assume only pertains to Spinning Reserve
# Ignoring MaxEnergy  -- don't understand what this means
# Ignoring RampRates (for now)
# Ignoring ReserveCap -- only pertains to reseres
# Ignoring "mustRun" file (for now)
# Ignoring inihour - Using warm starts for everyone for now

#"../Data/AndysGenInstance"
function loadISO( dir, filt_ratio=.3 )
	##Elements of the "others" collection include INC Decs, and loads
	## the file pResType indexes type of everything (not currently used)
	gens, others = createGens("$dir/SetRes.txt");
	addEcoMax!("$dir/pEcomax.txt", gens)
	addEcoMin!("$dir/pEcomin.txt", gens)
	addVal!("$dir/pMaxStart.txt", gens, :maxstart, 
	        test_condition=x->1 <= x <= HRS)
	addVal!("$dir/pIniHour.txt", gens, :initonline)
	addVal!("$dir/pMinUpTime.txt", gens, :minup, 
	        test_condition=x->0 <= x <= HRS)
	addVal!("$dir/pMinDnTime.txt", gens, :mindown, 
	        test_condition=x->0 <= x <= HRS)

	fixCycles!(gens)  

	addOffers!("$dir/pEnergyBidPrice.txt", "$dir/pEnergyBidSize.txt", gens)
	addStartCosts!("$dir/pStartUpCost.txt", gens)
	
	scaling = 1.
	if filt_ratio !=nothing
		scaling = filtGens!(gens, filt_ratio)
	end
	return gens, scaling
end

function filtGens!(gens, p)
	println("Initial Num Gens:\t", length(gens))
	hr13Cap = sum([getCap(g, 13) for g in values(gens)])
    for g in values(gens)
        if g.startcost < .1
            g.startcost = 0.
        end
    end
    #These generators seemingly have bad data...
    delete!(gens, "VLCO-I_ASCUTNEY_26_22")
    delete!(gens, "VLCO-I_N_RUTLND_44_22")

    #drop the wind and hydro for now
    for g in values(gens)
        if g.fueltype == "Wind"
            delete!(gens, g.name)
        elseif g.fueltype=="Hydro"
            delete!(gens, g.name)
        end
    end
    
    #halve the steam, CT and Diesel
    srand(8675309)
    for g in values(gens)
        if g.fueltype =="Steam"
            rand() > p && delete!(gens, g.name)
        elseif g.fueltype =="CT"
            rand() > p && delete!(gens, g.name)
        elseif g.fueltype == "Diesel"
            rand() > p && delete!(gens, g.name)
        end
    end
    #Generators which cause needless numerical instability
    for g in values(gens)
        blocks, prices = getCurve(g)
        if minimum(blocks) < 1e-3
            println(g.name)
            delete!(gens, g.name)
        end
    end
	newhr13Cap = sum([getCap(g, 13) for g in values(gens)])
	println("Final Num Gens:\t", length(gens))
	println("Scaling:\t", newhr13Cap / hr13Cap)
	return newhr13Cap / hr13Cap
end


#returns incBlocks, prices for the curve
function getCurve( g::Generator )
    nblocks = length(g.offerblocks)
    maxval = maximum(g.ecomax)
    blocks = [g.offerblocks[ix] for ix =1:nblocks]
    
    #clip the last block at ecomax
    amounts = [min(b.amount, maxval) for b in blocks]
    
    #check increasing costs
    amounts != sort(amounts) && error("Amounts not ordered: $(g.name)")
    prices = [b.price for b in blocks]
    prices != sort(prices) && error("Prices not ordered: $(g.name)")
    
    if nblocks == 1
        return amounts, prices
    else
        amt_diffs = amounts[2:nblocks] - amounts[1:nblocks-1]
        #drop the tail as necessary
        nblocks = sum(amt_diffs .> 0) + 1
        return [amounts[1]; amt_diffs[1:nblocks-1]], prices[1:nblocks]
    end
end    


function splitName( sname )
    #strip surrounding single quotes
    sname2 = strip(sname)
    sname2 = strip(sname2, '\'')
    for iFuel in FUEL_TYPES
        indx = searchindex(sname2, iFuel)
        if 0 != indx
            name = sname2[1:indx-1]
            return strip(name), iFuel
        end
    end
    #Probably not a fuel
    return sname2, ""
end

function parseVal( line ) 
    namefuel, val = split(line, '\t')
    name, fuel = splitName(namefuel)
    return name, fuel, float(val)
end

function parseHourlyVal(line)
    trial = split(line, ".'H")
    if length(trial)> 1
        namefuel, hrval = trial
    else
        namefuel, hrval = split(line, ".H")
    end
    hr, val = split(hrval)
    hr = strip(hr, '\'')
    name, fuel = splitName(namefuel)
    return name, int(hr), float(val)
end

#assumes that it's a H1*H36 type
function parseH24Block( line )
    #Try to split by as if a total
    namefuelblock, val = split(line, ".H1*H36")
    namefuel, block = split(namefuelblock, ".")
    name, fuel = splitName(namefuel)
    block = strip(block, '\'')
    block = block[4:]  #drops the "BLK" prefix
    return name, int(block), float(val)
end    

####
#Read the file SetRes.txt and add a dictionary of generators
function createGens(path)
	fp = open(path, "r")
	readline(fp) #burn the first line
	gens = Dict{String, Generator}()
	others = String[]

	while !eof(fp)
		line = readline(fp)
	    if contains(line, ";")
	        break
	    end
		name, fuel = splitName(line)
		if fuel != ""
			gens[name] = Generator(name, fuel)
		else
			push!(others, name)
		end
	end
	close(fp)
	return gens, others
end	

#uses the file ../Data/AndysGenInstance/pEcomax.txt
function addEcoMax!(path, gens; log=false, scaling=1e-3)
	fp = open(path)
	readline(fp)
	while !eof(fp)
	    line = readline(fp)
	    if contains(line, ";")
	        break
	    end
	    #parse the line
	    name, hr, val = parseHourlyVal(line)
	    if !haskey(gens, name)
	        if log 
	        	println("NonGen EcoMax: $name $hr $val")
	        end
	    else
			#Take advantage of the fact that data listed in order
	        assert(length(gens[name].ecomax) == hr-1)
	        push!(gens[name].ecomax, val * scaling)
	    end
	end
	close(fp)
end

#uses the file ../Data/AndysGenInstance/pEcomin.txt
function addEcoMin!(path, gens; log=false, scaling=1e-3)
	fp = open(path)
	readline(fp)
	while !eof(fp)
	    line = readline(fp)
	    if contains(line, ";")
	        break
	    end
	    #parse the line
	    name, hr, val = parseHourlyVal(line)
	    if !haskey(gens, name)
	        if log 
	        	println("NonGen EcoMin: $name $hr $val")
	        end
	    else
			#Take advantage of the fact that data listed in order
	        assert(length(gens[name].ecomin) == hr-1)
	        push!(gens[name].ecomin, val * scaling)
	    end
	end
	close(fp)
end

#Scalar values like: "../Data/AndysGenInstance/pMaxStart.txt"
function addVal!(path, gens, field; log=false, test_condition=x->true, convert=x->int(x))
	fp = open(path); readline(fp)
	while !eof(fp)
	    line = readline(fp)
	    if contains(line, ";")
	        break
	    end
	    name, fuel, val = parseVal(line)
	    if !haskey(gens, name)
	        if log 
	        	println("NonGen: $field $name $fuel $val")
	        end
	    else
	    	if test_condition(val)
	    		setfield(gens[name], field, convert(val))
			end
	    end
	end
	close(fp)
end

# fixes maxstarts to the minimal number feasible given up/dwn constraints
# The remove_dn functionality assumes that everything starts off
function fixCycles!(gens; remove_dn = false)  
	maxstarts(minup, mindown) = floor((HRS-1) / (minup + mindown)) + 1
	for g in values(gens)
		numstarts = maxstarts(g.minup, g.mindown)
		if numstarts < g.maxstart
			g.maxstart = numstarts
		end

		#ASSUMES ALL THINGS START OFF
		if remove_dn
			g.mindown = 0
		end
	end
end

#"../Data/AndysGenInstance/pEnergyBidPrice.txt"
#"../Data//AndysGenInstance/pEnergyBidSize.txt"
function addOffers!(price_path, size_path, gens; log=false, scaling=1e-3)
	#Right now ignoring bid info from everything in the "other" pile
	fprices = open(price_path) ; fsizes = open(size_path)
	readline(fprices) ; readline(fsizes)

	#assume that both files of same size
	while !eof(fprices)
	    pline = readline(fprices); sline = readline(fsizes)
	    if contains(pline, ";")
	        continue
	    end
	    if contains(pline, ".H1*H36")
	        #parse both.  Assert they're equal
	        pname, pblock, price = parseH24Block(pline)
	        sname, sblock, size = parseH24Block(sline)
	        assert( pname == sname ) ; assert(pblock == sblock)

	        #just add the blocks directly for now
	        #assumes they are added in order
	        if ! haskey(gens, sname)
	        	if log
		            println("Bad Offer: \t", sname)
		        end
	            continue
	        end
	        gens[sname].offerblocks[sblock] = Block(size * scaling, price)
	    else
	        #For now, toss things that aren't real generators
	    end
	end
	close(fprices); close(fsizes)
end

#"../Data/AndysGenInstance/pStartupCost.txt"
function addStartCosts!(path, gens; scaling=1e-3)
	fp = open(path); readline(fp)
	while !eof(fp)
	    line = readline(fp)
	    if contains(line, ";")
	        break
	    end

	    #parse the line
	    namefuel, blockval = split(line, ".BLK")
	    name, fuel = splitName(namefuel)
	    block, val = split(blockval)

	    if int(block) == 2  #only take warm starts for now.  Data must contain blocks for everyone
	        gens[name].startcost = float(val) * scaling
	    end
	end
	close(fp)
end



