"""
    Collects some helper functions to summarize output of a UC model
"""
from config import *
from collections import  Counter
def writeHourlyGens(file_out, label, prod_by_type):
    #These items are definitely here
    file_out.writerow([label, "Total"] + 
        [ prod_by_type["TOTAL", hr] for hr in xrange(HORIZON_LENGTH) ] )
    file_out.writerow([label, "Slack"] + 
        [ prod_by_type["Slack", hr] for hr in xrange(HORIZON_LENGTH) ] )
    file_out.writerow([label, "FLEX"] + 
        [ prod_by_type["FLEX", hr] for hr in xrange(HORIZON_LENGTH) ] )
    
    for fuel_type in FUEL_TYPES:
        if (fuel_type, 0) in prod_by_type:
            file_out.writerow([label, fuel_type] + 
                    [prod_by_type[fuel_type, hr] for hr in xrange(HORIZON_LENGTH) ] )           
    if ("INC", 0) in prod_by_type:
        file_out.writerow([label, "INC"] + 
                [ prod_by_type["INC", hr] for hr in xrange(HORIZON_LENGTH) ] )
    if ("DEC", 0) in prod_by_type:
        file_out.writerow([label, "DEC"] + 
            [ prod_by_type["DEC", hr] for hr in xrange(HORIZON_LENGTH) ] )
    file_out.writerow([label, "FLEX"] + 
        [ prod_by_type["FLEX", hr] for hr in xrange(HORIZON_LENGTH) ] )

def writeHourlySchedCap(file_out, label, on_vals, gen_dict):
    #tally the amount of scheduled capacity by fuel_type
    cap_tally = Counter()
    for (name, iHr), val in on_vals.items():
        if round(val) > .5:
            cap_tally[iHr, gen_dict[name].fuel_type ] += gen_dict[name].eco_max[iHr ]

    fuel_types = set( fuel for (hr, fuel) in cap_tally )
    for fuel in fuel_types:
        file_out.writerow([label, fuel] + [cap_tally[hr, fuel] for hr in xrange(HORIZON_LENGTH) ] )        
        
def calcResReqs(on_vals, gen_dict):
    cap_producers= {}
    #not quite right, but it's a reasonable to assume the biggest prods will be nukes that are always on
    for name, hr in on_vals:
        if name not in cap_producers:
            cap_producers[name] = gen_dict[name].eco_max[hr]
        else:
            cap_producers[name] = max( gen_dict[name].eco_max[hr], cap_producers[name] )
            
    t = sorted(cap_producers.items(), key = lambda (n, cap): cap, reverse=True)
    big1, big2 = (cap for n, cap in t[:2] )
    return .5 * big1, big1, big1 + .5 * big2

def calcOpLimits(on_vals, gen_dict):
    """Heuristically compute the min/max loads per hour for which this dispatch still optimal."""
    maxCap = Counter()
    for (name, iHr) in on_vals:
        if on_vals[name, iHr] > .9:
            maxCap[iHr] += gen_dict[name].eco_max[iHr]
    return maxCap
            
    