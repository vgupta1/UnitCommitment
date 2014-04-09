""" Generator

Defines a class for generators with some helper functions for initializing from .txt files
"""
import pdb, csv, math
from config import *

def getFuelType( name_fuel ):
    """Translates: AES LOND_GR_RIDGE_9_99Steam to ("AES LOND_GR_RIDGE_9_99", "Steam") """
    for type in FUEL_TYPES:
        if type in name_fuel:
            fuel_type = type
            indx = name_fuel.find(type)
            name = name_fuel[:indx]
            name = name.strip("'")
            break
    else: # Might be a load or something else.
        name = name_fuel
        fuel_type = "Other"
    return (squish(name), fuel_type)    

def squish( string ):
    return "".join(string.split())

class StartUpBlock:
    """ Simple Struct"""
    def __init__(self, offHr, cost, time, noti, blockNum):
        self.offHr, self.cost, self.time, self.noti, self.num = (
            offHr, cost, time, noti, blockNum )

class OfferBlock:
    """ Simple Struct"""
    def __init__(self, blockNum, size, price):
        self.size, self.price, self.num = size, price, blockNum

class GenUnit:
    """ A generation unit is defined by various pieces of info
        This struct collects all of this information
    """
    #Class Constant
    RESERVE_PRODUCTS = ("TMSR_Cap", "TMNSR_Cap", "TMOR_Cap")
    RESERVE_REQUIREMENTS = ("TMSR_Cap", "T10_Cap", "T30_Cap")
        
    def __init__(self, name, fuel_type):
        self.name = name.strip()

        #dictionaries indexed by hour
        self.eco_max, self.eco_min, self.fixed_energy={}, {}, {}

        #List of hours must run
        self.mustRun=[]

        #blocks are sorted in increasing order for knot value
        self.startupblocks = []

        #indexed by time.  each time has a list ordered by knot value
        self.offerBlocks = {}

        #constants
        self.reserveBidSize = 0.0
        self.fuel_type, self.res_type = fuel_type, None
        self.max_energy = None
        self.max_start = None
        self.init_online = None
        self.min_down, self.min_up = None, None
        self.no_load = None
        self.ramp_rate = None
        self.isFlexible = False
        self.isFlexibleReserve = False

###########################################
#aimed at the setRes.txt file
def createGenerators( resources_path ):
    gen_dict = {}
    gen_file = open(resources_path, 'rU')
    gen_file.readline()    #burn a line.
    for line in gen_file:
        if "/;" in line:
            break         #EOF Character
        gen_name, fuel_type = getFuelType( line.strip() )
        if gen_name in gen_dict:
            raise RunTimeError("Generator %s listed twice." % gen_name)
        else:
            gen_dict[gen_name] = GenUnit(gen_name, fuel_type)
    return gen_dict

def parseGen( line, hasTime, isFloatParam, scaling ):
        t = line.split()
        s_param = t[-1]
        indx = line.rfind(s_param)
        gen_name = line[:indx]
        if isFloatParam:
            param = float(s_param) * scaling
        else:
            param = s_param
        
        if hasTime:
            gen_name, time = gen_name.split(".")
            time = time.strip()
            time = time.strip("'")        
            if "H" in time:
                time = int(time.lstrip("H")) - 1
            else:
                raise RuntimeError("Time Indx without 'H' ")
        else:
            time = None
        gen_name, fuel_type = getFuelType( gen_name.strip() )
        gen_name = gen_name.strip("'")
        return squish(gen_name), fuel_type, param, time
        
def addGenParameter( param_path, param_name, gen_dict, hasTime=True, isFloatParam=True, scaling=1. ):
    """Assumes a listing by generator and time and a parameter.  Useful for Ecomax, econmin"""
    gen_file = open(param_path, 'rU')
    gen_file.readline()    #burn a line.
    for line in gen_file:
        if "/;" in line:
            break       #EOF character
        gen_name, fuel_type, param, time = parseGen( line, hasTime, isFloatParam, scaling )
        if gen_name not in gen_dict:
            raise RuntimeError("Generator %s not found" % gen_name)
        else:
            g = gen_dict[ gen_name ]
    
        if hasTime:
            param_dict = getattr(g, param_name)
            if time in param_dict:
                raise RuntimeError("Parameter set twice.")
            param_dict[time] = param
        else:
            if getattr(g, param_name) is not None:
                raise RuntimeError("Parameter set Twice.")
            setattr(g, param_name, param)
    return gen_dict

#VG Ideally, this should be lumped into the functionality above.
def addReserveCap(param_path, gen_dict, scaling=1.):
    gen_file = open(param_path, 'rU')
    gen_file.readline() #burn a line.
    for line in gen_file:
        if "/;" in line:
            break   #EOF Character
        t = line.split()
        s_value = t[-1]
        indx = line.rfind(s_value)
        gen_name = line[:indx].strip()
        gen_name, cap_type, time = gen_name.split(".")
        if cap_type not in GenUnit.RESERVE_REQUIREMENTS:
            raise RuntimeError("Cap type %s not found" % cap_type)
        gen_name = gen_name.strip("'")
        gen_name, fuel_type = getFuelType(gen_name)
        gen_name = squish(gen_name)
        if gen_name not in gen_dict:
            raise RuntimeError("Generator %s not found" % gen_name)
        else:
            g = gen_dict[ gen_name ]
        
        #May already be set because data is by hour and we ignore hourly
        if hasattr(g, cap_type):
            assert abs( getattr(g, cap_type) - float(s_value) * scaling) < 1e-8
        setattr(g, cap_type, float(s_value) * scaling)

    return gen_dict

def addFlexResources( param_path, gen_dict):
    gen_file = open(param_path, 'rU')
    gen_file.readline()     #burn a line.
    for line in gen_file:
        if "/;" in line:
            break         #EOF character
        gen_name, fuel_type = getFuelType( line.strip() )
        gen_name = gen_name.strip("'")
        if gen_name not in gen_dict:
            raise RunTimeError("Flex Resourece %s not found." % gen_name)
        else:
            gen_dict[ gen_name ].isFlexible = True
    return gen_dict

def addFlexReserves(param_path, gen_dict):
    """ This is the set under SetFlexReserve. Appears to be a subest of 
    SetFlexEgr.  Everything is one of the Type ((TMSR,TMNSR,TMOR)"""
    gen_file = open(param_path, 'rU')
    gen_file.readline()    #burn a line.
    for line in gen_file:
        if "/;" in line:
            break       #EOF character
        gen_name, t = line.split(".")
        gen_name, fuel_type = getFuelType( gen_name.strip("'") )
        if gen_name not in gen_dict:
            raise RunTimeError("Flex Reserve %s not found." % gen_name)
        else:
            gen_dict[ gen_name ].isFlexibleReserve = True
    return gen_dict

def addMustRunHours(param_path, gen_dict):
    gen_file = open(param_path, 'rU')
    gen_file.readline()     #burn a line
    for line in gen_file:
        if "/;" in line:
            break       #EOF Characer
        gen_name, hour = line.split(".")
        gen_name, fuel_type = getFuelType( gen_name.strip() )
        hour = hour.strip()
        gen_name = gen_name.strip("'")

        if gen_name not in gen_dict:
            raise RunTimeError("Generator %s not found." % gen_name)
        else:
            gen_dict[gen_name].mustRun.append( hour )
    return gen_dict

def addReserveBidSizes(param_path, gen_dict):
    gen_file = open(param_path, 'rU')
    gen_file.readline()    #burn a line.
    for line in gen_file:
        if "/;" in line:
            break   #EOF CHARACTER
        t = line.split()
        s_value = t[-1]
        indx = line.rfind(s_value)
        gen_name = line[:indx].strip()
        gen_name = gen_name.split(".")[0]
        gen_name, fuel_type = getFuelType( gen_name.strip("'") )
        if gen_name not in gen_dict:
            raise RuntimeError("Generator %s not found" % gen_name)
        else:
            gen_dict[ gen_name ].reserveBidSize = float( s_value )
    return gen_dict

def addStartUpBlocks( offHoursCurve_path, costs_path, time_path, 
                                            notification_path, gen_dict, scaling=1.):
    file_offHrs = open( offHoursCurve_path, 'rU')
    file_offHrs.readline() #burn a line
    file_costs = open( costs_path, 'rU')
    file_costs.readline() #burn a line
    file_times = open( time_path, 'rU')
    file_times.readline() #burn a line
    file_noti = open( notification_path, 'rU')
    file_noti.readline() #burn a line
    
    for line_offHrs, line_costs, line_times, line_noti in zip(
                file_offHrs, file_costs, file_times, file_noti):
        if "/;" in line_offHrs:
            break   #EOF Character
        t = line_offHrs.split()
        s_value = t[-1]
        indx_offHrs = line_offHrs.rfind(s_value)
        t = line_costs.split()
        s_value = t[-1]
        indx_costs = line_costs.rfind(s_value)
        t = line_times.split()
        s_value = t[-1]
        indx_times = line_times.rfind(s_value)
        t = line_noti.split()
        s_value = t[-1]
        indx_noti = line_noti.rfind(s_value)

        #logic relies on the fact that everything is listed in the same order.
        assert ( line_offHrs[:indx_offHrs].strip() == 
                      line_costs[:indx_costs].strip() == 
                      line_times[:indx_times].strip() == 
                      line_noti[:indx_noti].strip() )

        gen_name, blocknum = line_offHrs[:indx_offHrs].split(".")
        gen_name, dummy = getFuelType( gen_name.strip("'") )
        blocknum = blocknum.strip().lstrip("BLK")
        
        if gen_name not in gen_dict:
            raise RuntimeError("Generator %s not found" % gen_name)
        else:
            offHr = float( line_offHrs[indx_offHrs:] )
            cost = float( line_costs[indx_costs:] ) * scaling
            time = float( line_times[indx_times:] )
            noti = float( line_noti[indx_noti:] )
            num = int( blocknum )
            gen_dict[ gen_name ].startupblocks.append( StartUpBlock(offHr, cost, time, noti, num) )

    #go back through and reorder everybody
    for g in gen_dict.values():
        g.startupblocks.sort( key = lambda block : block.num )

    return gen_dict

def addOfferBlocks( price_path, size_path, gen_dict, mw_scaling=1., price_scaling=1.):
    """ Adds the hourly offer curves """
    file_price = open( price_path, 'rU')
    file_price.readline() #burn a line
    file_size = open( size_path, 'rU')
    file_size.readline() # burn a line
    for line_price, line_size in zip(file_price, file_size):
        if "/;" in line_price:
            break       #EOF CHARACTER
        t = line_price.split()
        s_value = t[-1]
        indx_price = line_price.rfind(s_value)
        t = line_size.split()
        s_value = t[-1]
        indx_size = line_size.rfind(s_value)

        assert line_price[:indx_price].strip() == line_size[:indx_size].strip()
        gen_name, block, time = line_price[:indx_price].split(".")
        gen_name, fueltype = getFuelType( gen_name.strip("'") )
        block_num = int( block.strip("'").lstrip("BLK") )
        time = time.strip()
        
        if gen_name not in gen_dict:
            raise RuntimeError("Gen %s not found" % gen_name )
        if time not in gen_dict[gen_name].offerBlocks:
            gen_dict[gen_name].offerBlocks[time] = []

        block = OfferBlock(block_num, float(line_size[indx_size:]) * mw_scaling, float(line_price[indx_price:] ) * price_scaling )
        gen_dict[gen_name].offerBlocks[time].append( block ) 

    for g in gen_dict.values():
        for blocks in g.offerBlocks.values():
            blocks.sort( key= lambda block : block.num )

def fixRampRates(true_gens):
    for name, gen in true_gens.items():
        ramp_rate = gen.ramp_rate
        if ramp_rate is None:
            print name, "No Ramp"
            continue
        max_diff = ramp_rate
        for iHr in xrange(HORIZON_LENGTH):
            eco_min = gen.eco_min[iHr] 
            eco_max_p = gen.eco_max[iHr + 1]
            #check the ramping down constraints
            max_diff = max(eco_min - eco_max_p, max_diff)

            eco_max = gen.eco_max[iHr] 
            eco_min_p = gen.eco_min[iHr + 1]
            #check the ramping up constraints
            max_diff = max(eco_min_p - eco_max, max_diff)

        gen.ramp_rate = max_diff
        
def doEverything():
    gen_dict = createGenerators("../AndysGenInstance/SetRes.txt")
    addGenParameter("../AndysGenInstance/pResTyp.txt", "res_type", gen_dict, hasTime=False, isFloatParam=False)
    true_gens = dict( filter( lambda (name, gen): gen.res_type == "GEN" and gen.fuel_type <> "FixedImport", 
                                        gen_dict.items() ) )

    addGenParameter( "../AndysGenInstance/pEcomax.txt", "eco_max", gen_dict, scaling=1e-3 )
    addGenParameter( "../AndysGenInstance/pEcomin.txt", "eco_min", gen_dict, scaling=1e-3 )
    addGenParameter( "../AndysGenInstance/pMaxEnergy.txt", "max_energy", gen_dict, hasTime=False, scaling=1e-3 )
    addGenParameter("../AndysGenInstance/pMaxStart.txt", "max_start", gen_dict, hasTime=False)
    addGenParameter("../AndysGenInstance/pIniHour.txt", "init_online", gen_dict, hasTime=False)
    addGenParameter("../AndysGenInstance/pMinDnTime.txt", "min_down", gen_dict, hasTime=False)
    addGenParameter("../AndysGenInstance/pMinUpTime.txt", "min_up", gen_dict, hasTime=False)
    addGenParameter("../AndysGenInstance/pNoloadcost.txt", "no_load", gen_dict, hasTime=False, scaling=1e-3)
    addGenParameter("../AndysGenInstance/pRampRate.txt", "ramp_rate", gen_dict, hasTime=False, scaling=1e-3)
    addReserveCap("../AndysGenInstance/pReserveCapacity.txt", gen_dict, scaling=1e-3)
    
    ##VG Treat these quickly
    addGenParameter("../AndysGenInstance/pFixEnergy.txt", "fixed_energy", gen_dict, hasTime=True, scaling=1e-3)
    addFlexResources( "../AndysGenInstance/SetFlxEgr.txt", gen_dict)
    addMustRunHours("../AndysGenInstance/SetMustRun.txt", gen_dict)
    addFlexReserves("../AndysGenInstance/SetFlexReserve.txt", gen_dict)

    addStartUpBlocks( "../AndysGenInstance/pOfflineHr.txt", 
                                        "../AndysGenInstance/pStartUpCost.txt", 
                                        "../AndysGenInstance/pStartupTime.txt", 
                                        "../AndysGenInstance/pNotTime.txt", gen_dict, scaling=1e-3)
    addOfferBlocks( "../AndysGenInstance/pEnergyBidPrice.txt", "../AndysGenInstance/pEnergyBidSize.txt", gen_dict)


    fixRampRates(true_gens)

    for g in true_gens.values():
        if max(g.eco_max.values() ) < .001:
            del gen_dict[ g.name ] 

    return gen_dict

def smallTestCase( gen_dict, filt_percent=.1):
    """Generates a smaller test case to test algorithms
    Takes given percentage of each fuel-type and the flex loads and fixed loads
    """
    true_gens = filter(lambda (n, g): g.res_type == "GEN" and g.fuel_type <> "FixedImport", 
                                    gen_dict.items() )
    gens_filt = {}
    tot_cap, tot_new_cap = 0., 0.
    #For each fuel_type, filter down
    for fuel in FUEL_TYPES:
        gen_fuels = filter(lambda(n, g): g.fuel_type == fuel, true_gens )
        if fuel == "Nuclear":
            gen_fuels_filt  = gen_fuels[:3 ]
        else:
            gen_fuels_filt  = gen_fuels[:int(math.ceil(len(gen_fuels) * filt_percent)) ]

        orig_cap = sum( max(g.eco_max.values() ) for n, g in gen_fuels)
        new_cap = sum( max(g.eco_max.values() ) for n, g in gen_fuels_filt)
        gens_filt.update( gen_fuels_filt )
        tot_cap += orig_cap
        tot_new_cap += new_cap
        print fuel, orig_cap, new_cap, new_cap / (orig_cap + .0001 )

    flex_loads = filter(lambda (n, g): g.res_type =="LOAD" and g.isFlexible, gen_dict.items() )
    flex_loads_filt = flex_loads[:int(math.ceil(len(flex_loads) * filt_percent))] 
    gens_filt.update( flex_loads_filt )
    orig_cap = sum( max(g.eco_max.values() ) for n, g in flex_loads )
    new_cap = sum( max(g.eco_max.values() )  for n, g in flex_loads_filt )
    tot_cap -= orig_cap
    tot_new_cap -= new_cap
    print "Flex Loads", orig_cap, new_cap, new_cap / orig_cap
    print "total:", tot_cap, tot_new_cap, tot_new_cap / tot_cap
    load_ratio = tot_new_cap / tot_cap
    
    return gens_filt, load_ratio

###########################################
if __name__ == "__main__":
#    pdb.set_trace()
    gen_dict = doEverything()
    print len(filter(lambda g: g.init_online is None, gen_dict.values() ) ), len(gen_dict)

    fixRampRates(gen_dict)


    print "\n \n Take 2 \n \n"
    fixRampRates(gen_dict)

    sys.exit()

    #Notes
    #VG Weirdly there are things with no eco-max or eco-min
    # There are also 36 hours worth of data?
    #Is the fixed energy set the complement of flexible energy set?
    # What are the loads/fixed imports?
    
    #let's do some sanity checks on the curves


    # many items don't have max energy either...
    # other items dont' have max_start, or have max_start = 0
    # other items don't have init_oinline

####################################    
    file_out = csv.writer(open("resource_details.csv", "w"))
    
    headers = ["Name", "FuelType", "ResType", 
                        "NoLoadCost", "InitOnline",
                        "IsFlexible", "IsFlexibleReserve", 
                        "MaxStart", "MaxEnergy", 
                        "MinDown", "MinUp", "RampRate"]
    
    ## Make a serpate file for: EcoMin, EcoMax, "FixedEnergy",  "MustRun"
    
    file_out.writerow( headers )
    for ix, g in enumerate(gen_dict.values()):
        line = [ix, g.name, g.fuel_type, g.res_type, 
                    g.no_load, g.init_online, 
                    g.isFlexible, g.isFlexibleReserve, 
                    g.max_start, g.max_energy, 
                    g.min_down, g.min_up, g.ramp_rate]
        file_out.writerow(line)
        print line
    
#######################################
    #Organize a dump of the time dependent values
    file_out = csv.writer( open("resource_details_time.csv", "w") )
    headers = ["Name", "ResType", "FuelType", "Parameter", "Time", "Value"]
    file_out.writerow( headers )
    
    for g in gen_dict.values():
        stub = [g.name, g.res_type, g.fuel_type]
        for time, value in g.eco_max.items():
            file_out.writerow( stub + ["EcoMax", time, value] )
    
        for time, value in g.eco_min.items():
            file_out.writerow( stub + ["EcoMin", time, value] )
    
        for time, value in g.fixed_energy.items():
            file_out.writerow( stub + ["FixedEnergy", time, value] )
    
        for time in g.mustRun:
            file_out.writerow( stub + ["MustRun", time, 1] )
    
    
    
