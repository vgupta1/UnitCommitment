###
# Reading in Gen Info
###
# File defines the relevant generator types
FUEL_TYPES = ("Steam", "CT", "CC", "Diesel", "Hydro", "Nuclear", "Wind")
HRS = 24
type Block
	amount::Float64
	price::Float64
end

type Generator
	name::String
	fueltype::String
	startcost::Float64   #Assumes warm starts for everyone for now
	maxstart::Int
	initonline::Int
	mindown::Int
	minup::Int
	offerblocks::Dict{Int, Block}  #assumes that all true generators have H1*H36 offer specs. Int is the block num

	#of Length 24
	ecomin::Vector{Float64}
	ecomax::Vector{Float64}
end
Generator(name, fueltype) = Generator(name, fueltype, 0., HRS, 0, 0, 0, Dict{Int, Block}(), Float64[], Float64[])

getCap(gen::Generator, hr) = gen.ecomax[hr]
getCap(gen::Generator) = maximum(gen.ecomax)