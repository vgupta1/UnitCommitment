"""BuildUC2

Encapsulates the logic to build the various UC models.  Nominal logic only goes in this file
New file for affine logic.
New file for solving and experiments
"""

import csv, pdb, numpy, sys
import gurobipy as grb
from collections import  Counter
import generator, readData 
from AffineHelpers import *

HORIZON_LENGTH = 24
EPSILON_ZERO = 0.
PENALTY = 5
TOL = 1e-5

## To Do
#Redo the read in data to only use hrs 0 to 23 to index everything.  No "Hr24" stuff.
# Cold/warm/hot starts
# Incorporate startup times correctly
# Revisit the eco-min max constraints.  Are these locally ideal?  Use union of polyhedra results.

#To speed up
#Do a single filter for true_gens for each of the functions
# rampingConstsAffine does multiple variable adds... can reduce.
# Variable elimination:
  # Gens with identical eco-min eco-max
  # Gens with zero eco_max in particular hours
  # a couple palce where we block add variables we can remove them after adding if they're not used...
#Rescale the affine policies so that they are of order 1 instead of order large....

## Notes
# The "GEN" with "FixedImport" have no generator characteristics, but have fixedEnergy.  We ignore
# Amounts in the hours after 24....  hence we just drop them from the data set.
# We are ignoring network structure for now
# For the incremental bid, actually vary by hour, but single block cost structure

class UCNomModel():
    def __init__(self):
        self.model = None
        self.on_vars, self.start_vars, self.stop_vars = None, None, None
        self.prod_vars = None
        self.reserve_vars = None
        self.variable_cost_vars = None
        self.flex_loads = None
        self.fixed_cost_var = None
        self.dec_vars_amt = None
        self.inc_vars = None

    def rampingCallback(self, gen_dict, where):
        """Callback to be used when using sparse ramping constraints"""
        model = self.model
        if where == grb.GRB.callback.MIPSOL:
            for name, gen in gen_dict.items():
                if gen.res_type <> "GEN" or gen.fuel_type in ("FixedImport", "Steam", "CT"):
                    continue
                for hr in xrange(1, HORIZON_LENGTH):
                    #check to see if current solution violates either ramping constraint
                    prod_, prod, start, stop = model.cbGetSolution(
                        [self.prod_vars[name, hr-1], self.prod_vars[name, hr], self.start_vars[name, hr], self.stop_vars[name, hr]] )        
                    if prod - prod_ > gen.ramp_rate and not start > .9:
                        eco_max = gen.eco_max["H" + str(hr + 1)]
                        model.addConstr( self.prod_vars[name, hr] - self.prod_vars[name, hr - 1] <= 
                                                            gen.ramp_rate + eco_max * self.start_vars[name, hr], 
                                                            name= "RampUp" + name + "H" + str(hr) )
                    if prod_ - prod > gen.ramp_rate and not stop > .9:
                        eco_max_m = gen.eco_max["H" + str(hr)] 
                        model.addConstr( self.prod_vars[name, hr -1] - self.prod_vars[name, hr] <= 
                                                            gen.ramp_rate + eco_max_m * self.stop_vars[name, hr], 
                                                            name="RampDown" + name + "H" + str(hr) )

def hrToInt( sHr ):
    """ Convert "H12" to 12 """
    return int( sHr.lstrip("H") )

######### 
# Helper functions for creating the nominal model
#########

def genStage1Vars( model, gen_dict, init_on={}, init_start={} ):
    """ Creates all stage 1 (committment) variables
    returns dicts:  on[gen, time], start[gen, time], stop[gen, time], cost_var
    cost_var is an auxiliary variable that allows us to access the first stage costs directly
    """
    on_vars, start_vars, stop_vars = {}, {}, {}
    cost_var_expr = grb.LinExpr()
    for name, gen in gen_dict.items():
        if gen.res_type <> "GEN" or gen.fuel_type == "FixedImport":
            continue

        for iHr in xrange(HORIZON_LENGTH ):
                sHr = "H" + str(iHr)
                on_vars[name, iHr] = model.addVar( vtype=grb.GRB.BINARY, name="On" + name + sHr )
                stop_vars[name, iHr] = model.addVar( vtype=grb.GRB.BINARY, name="Stop" + name + sHr)  
                assert gen.no_load >= 0
    
                #for now do a simple warm start for everyone
                #VG expand to cold, warm, Hot
                warm_start_block = gen.startupblocks[1]
                assert warm_start_block.cost >= 0
                start_vars[name, iHr] = model.addVar( vtype=grb.GRB.BINARY, name="Start" + name + sHr)

                cost_var_expr += (on_vars[name, iHr] * gen.no_load + 
                                                start_vars[name, iHr] * warm_start_block.cost + 
                                                stop_vars[name, iHr] * EPSILON_ZERO) #dataset has no shutdown costs

    cost_var = model.addVar(obj = 1.0)
    model.update()  
    model.addConstr( cost_var >= cost_var_expr )

    #only bother if they've been set
    if init_on:
        for name, gen in gen_dict.items():
            if gen.res_type <> "GEN" or gen.fuel_type == "FixedImport":
                continue
            for iHr in xrange(HORIZON_LENGTH ):
                on_vars[name, iHr].start = init_on[name, iHr]
                start_vars[name, iHr].start = init_start[name, iHr]                    

    return on_vars, start_vars, stop_vars, cost_var

def addPiecewiseCosts(model, gen_dict, prod_vars):
    """Creates variables and the defining constraints for the piecewise cost functions
    for the variable cost """
    ys_all, cost_vars = {}, {}
    for name, gen in gen_dict.items():
        if gen.res_type <> "GEN" or gen.fuel_type == "FixedImport":
            continue
        #There is a nicety that all the true generators pertain H1*H36
        assert gen.offerBlocks.keys() == ['H1*H36']
        offerBlocks = gen.offerBlocks['H1*H36'] 
        #For numerical stability, makes sense to clip these at eco_max.
        eco_max = max( gen.eco_max.values() )
        for numKnots, b in enumerate(offerBlocks):
            if b.size >= eco_max:
                break
#         sizes = [b.size for b in offerBlocks[:ix] ] + [ eco_max ] 
#         prices = [b.price for b in offerBlocks[:ix] ] + [ offerBlocks[ix].price ] 
#         size_diffs = [s  - sm for s, sm in zip(sizes, [0] + sizes) ] 
#        numKnots = len(size_diffs)
        numKnots += 1

        #simpler formulation relies on the fact that size_diffs are all > 0
        for iHr in xrange(HORIZON_LENGTH):
            ys_all[name, iHr] = [model.addVar(ub=1.0) for ix in xrange(numKnots ) ] 
            cost_vars[name, iHr] = model.addVar(obj=1.0)

    #single call to save time at memory
    model.update()
    for name, gen in gen_dict.items():
        if gen.res_type <> "GEN" or gen.fuel_type == "FixedImport":
            continue
        offerBlocks = gen.offerBlocks['H1*H36'] 
        #VG for numerical stability, makes sense to clip these at eco_max.
        eco_max = max( gen.eco_max.values() )
        for numKnots, b in enumerate(offerBlocks):
            if b.size >= eco_max:
                break
        sizes = [b.size for b in offerBlocks[:numKnots] ] + [ eco_max ] 
        prices = [b.price for b in offerBlocks[:numKnots] ] + [ offerBlocks[numKnots].price ] 
        size_diffs = [s  - sm for s, sm in zip(sizes, [0] + sizes) ] 
    
        for iHr in xrange(HORIZON_LENGTH):     
            model.addConstr( cost_vars[name, iHr] == 
                    grb.quicksum( y * p  * s for y, p, s  in zip(ys_all[name, iHr], prices, size_diffs) ) )    
            model.addConstr( prod_vars[name, iHr] == 
                    grb.quicksum( y * s for y, s in zip(ys_all[name, iHr], size_diffs) ) )
    return cost_vars

def startStopConstraints(m, gen_dict, on_vars, start_vars, stop_vars):
    for gen in gen_dict.values():
        if gen.res_type <> "GEN" or gen.fuel_type == "FixedImport":
            continue
        name = gen.name

        #Treat the first time slice specially
        if gen.init_online <= 0 :
            #if its running in the first hour, needs a start
            m.addConstr( on_vars[name, 0] <= start_vars[name, 0] )
        else:
            #if its on and stopping, needs a stop
            m.addConstr( 1 - on_vars[name, 0] <= stop_vars[name, 0] )
            m.addConstr( stop_vars[name, 0] <= 1 - on_vars[name, 0] )

        #add constraints for all other time slices.
        for iHr in xrange(1, HORIZON_LENGTH):
            m.addConstr( on_vars[name, iHr - 1] - on_vars[name, iHr] + start_vars[name, iHr] >= 0 )
            m.addConstr( on_vars[name, iHr] - on_vars[name, iHr - 1] + stop_vars[name, iHr] >= 0 )
            m.addConstr( stop_vars[name, iHr] <= 1 - on_vars[name, iHr] )
            #check for eco_max values
            #the on_vars is handled in the eco-min_max constraints
            if gen.eco_max["H" + str(iHr + 1 )] < TOL:
                m.addConstr( start_vars[name, iHr] == 0 )
    return m            

def minUpConstraints(model, gen_dict, on_vars):
    for gen in gen_dict.values():
        if gen.res_type <> "GEN" or gen.fuel_type == "FixedImport":
            continue
        name = gen.name
        
        if gen.min_up <= 1:
            continue
        #Treat the first time slice specially
        #if its been on for less than min-run, keep it on.
        if 0 <= gen.init_online < gen.min_up: #VG Check the off by one here
            for tau in xrange(int(gen.min_up - gen.init_online)):
                model.addConstr( on_vars[name, tau] ==1 )                
        elif 0 >= gen.init_online:
           #if its starting up, force it to stay on
            for tau in xrange( int(gen.min_up) ):
                model.addConstr( on_vars[name, 0] <= on_vars[name, tau] )
    
        #add constraints for all other time slices
        for iHr in xrange(1, HORIZON_LENGTH):
            for tau in xrange( iHr + 1, min( iHr + int(gen.min_up), HORIZON_LENGTH) ):
                model.addConstr( on_vars[name, iHr] - on_vars[name, iHr - 1] <=
                                        on_vars[name, tau] )

def minDownConstraints(model, gen_dict, on_vars):
    for gen in gen_dict.values():
        if gen.res_type <> "GEN" or gen.fuel_type == "FixedImport":
            continue
        name = gen.name
        
        if gen.min_down <= 1:
            continue
        #Treat the first time slice specially
        #if its been off for less than min-down, keep it off.
        if gen.init_online < 0 and -gen.init_online < gen.min_up: #VG Check the off by one here
            for tau in xrange(int(gen.min_down + gen.init_online)):
                model.addConstr( on_vars[name, tau] == 0 )                
        elif 0 <= gen.init_online:
           #if its shutting down, force it to stay down
            for tau in xrange( int(min(gen.min_down, HORIZON_LENGTH)) ):
                model.addConstr( on_vars[name, 0] >= on_vars[name, tau] )
    
        #add constraints for all other time slices
        for iHr in xrange(1, HORIZON_LENGTH):
            for tau in xrange( iHr + 1, min( iHr + int(gen.min_up), HORIZON_LENGTH) ):
                model.addConstr( on_vars[name, iHr - 1] - on_vars[name, iHr] <=
                                        1 - on_vars[name, tau] )

def genIncVarsNom( model, gen_dict):
    """Creates the variables for the incremental bids"""
    inc_vars = {}
    for name, gen in gen_dict.items():
        if gen.res_type <> "INC":
            continue
        #VG Think about changing these to be unit size...
        for iHr, blocks in gen.offerBlocks.items():
            inc_vars[name, iHr] = model.addVar( ub = blocks[0].size, 
                                                                                      obj = blocks[0].price, 
                                                                                name = "Inc" + name + "H" + str(iHr) )
    return inc_vars

def genDecVarsNom(model, gen_dict):
    """Creates the dec bids.  notice that we assume the cost structure here
    is in terms of blocks, different form bid curve of true generators."""
    #generate a big set of dictionaries at expense of memory to do only one model update
    dec_vars_amt, dec_vars_cost = {}, {}
    ys_all = {}
    #iterate twice so only one model update
    for name, gen in gen_dict.items():
        if gen.res_type <> "DEC":
            continue
            
        #VG Think about changing these to be unit size...
        for iHr, blocks in gen.offerBlocks.items():
            dec_vars_amt[name, iHr] = model.addVar()
            dec_vars_cost[name, iHr] = model.addVar(obj=-1.0)
            ys_all[name, iHr] = [ model.addVar(ub = b.size)  for b in blocks]

    model.update()
    for name, gen in gen_dict.items():
        if gen.res_type <> "DEC":
            continue
        for iHr, blocks in gen.offerBlocks.items():
            model.addConstr( dec_vars_amt[name, iHr] == 
                    grb.quicksum( ys_all[name, iHr] ) )
            model.addConstr( dec_vars_cost[name, iHr] == 
                    grb.quicksum( y * b.price for y, b in zip(ys_all[name, iHr], blocks) ) )

    return dec_vars_amt, dec_vars_cost    

def genStage2VarsNom( model, gen_dict, useReserve):
    """Creates the continuous variables for nominal model.
    prod[gen, time], reserve[gen, time, type]
    """
    prod, reserve = {}, {}
    for name, gen in gen_dict.items():
        if gen.res_type <> "GEN" or gen.fuel_type == "FixedImport":
            continue

        for iHr in range(HORIZON_LENGTH):
            prod[name, iHr] = model.addVar(name="Prod" + name + "H" + str(iHr) )

            if not useReserve:
                continue
            for cap_type in generator.GenUnit.RESERVE_PRODUCTS:
                    reserve[name, iHr, cap_type] = model.addVar(name="Res" + name + "H" + str(iHr) + cap_type)
    model.update()
    return prod, reserve

def genFlexibleDemandVarsNom( model, gen_dict ):
    """Creates 2nd stage varibles for the flexible demands... Only 7 of them currently.
    Interpret these as revenue earned, and extra load you must satisfy.  
    """
    flex_loads = {}
    for name, gen in gen_dict.items():
        if gen.res_type <> "LOAD" or not gen.isFlexible:
            continue
        
        #these flex loads have a single block bid structure
        for iHr in range(HORIZON_LENGTH):
            flex_loads[name, iHr] = model.addVar(lb=gen.eco_min["H" + str(iHr + 1)], 
                                                                                 ub = gen.eco_max["H" + str(iHr + 1)], 
                                                                                 obj=-gen.offerBlocks.values()[0][0].price, 
                                                                                 name="FlexLoad" + name + "H" + str(iHr) )
    return flex_loads

def reserveRequirements(model, gen_dit, reserve_vars, TMSR_REQ, T10_REQ, T30_REQ):
    for iHr in xrange(HORIZON_LENGTH):
        TMSR_vars_by_hr = filter( lambda (name, hr, type): hr == iHr and type == "TMSR_Cap", 
                                                        reserve_vars.keys() )
        model.addConstr( grb.quicksum(reserve_vars[k] for k in TMSR_vars_by_hr ) >= TMSR_REQ, 
                                         name="TMSR_CAP " + "H" + str(iHr) )
        T10_vars_by_hr = filter( lambda (name, hr, type): hr == iHr and 
                                                    type in ("TMSR_CAP", "TMNSR_Cap"), 
                                                    reserve_vars.keys() )
        model.addConstr( grb.quicksum(reserve_vars[k] for k in T10_vars_by_hr ) >= T10_REQ, 
                name="T10_Req" + "H" + str(iHr) )
        T30_vars_by_hr = filter( lambda (name, hr, type): hr == iHr and 
                                                    type in ("TMSR_CAP", "TMNSR_Cap", "TMOR_CAP"), 
                                                    reserve_vars.keys())
        model.addConstr( grb.quicksum(reserve_vars[k] for k in T30_vars_by_hr ) >= T30_REQ, 
                                name="T10_Req" + "H" + str(iHr) )
    return []

def reserveCapacity(model, gen_dict, reserve_vars, res_cnsts):
    """ No generator can offer more reserve of any type than its capacity allows."""
    for name, gen in gen_dict.items():
        if gen.res_type <> "GEN" or gen.fuel_type == "FixedImport":
            continue
        #Notice that the TMSR cap is already handled by the on-off
        for iHr in xrange(HORIZON_LENGTH):
            model.addConstr( reserve_vars[name, iHr, "TMSR_Cap"] + 
                                                reserve_vars[name, iHr, "TMNSR_Cap"] <= gen.T10_Cap )
        for iHr in xrange(HORIZON_LENGTH):
            model.addConstr( reserve_vars[name, iHr, "TMSR_Cap"] + 
                                                reserve_vars[name, iHr, "TMNSR_Cap"] +
                                                reserve_vars[name, iHr, "TMOR_Cap"] <= gen.T30_Cap )
    return res_cnsts

def rampingConsts(model, gen_dict, prod_vars, start_vars, stop_vars, sparse):
    for name, gen in gen_dict.items():
        if gen.res_type <> "GEN" or gen.fuel_type == "FixedImport":
            continue
        if sparse and gen.fuel_type not in ("Steam", "CT"): #only add ramping for marginal turbines
            continue
        #ignore ramping constraints for the first slice.
        for hr in xrange(1, HORIZON_LENGTH):
            ##VG This logic is partially duplicated in "fixRampRates of generator.py"
            eco_min_m = gen.eco_min["H" + str(hr)]
            eco_max_m = gen.eco_max["H" + str(hr)] 
            eco_min = gen.eco_min["H" + str(hr + 1)]
            eco_max = gen.eco_max["H" + str(hr + 1)]
            if eco_max <= TOL: #won't run anyway
                continue
            if gen.ramp_rate < eco_max_m - eco_min:
                model.addConstr( prod_vars[name, hr -1] - prod_vars[name, hr] <= 
                                                    gen.ramp_rate + eco_max_m * stop_vars[name, hr], 
                                                    name="RampDown" + name + "H" + str(hr) )
            if gen.ramp_rate < eco_max - eco_min_m:        
                model.addConstr( prod_vars[name, hr] - prod_vars[name, hr - 1] <= 
                                                    gen.ramp_rate + eco_max * start_vars[name, hr], 
                                                    name= "RampUp" + name + "H" + str(hr) )
    return model    

def ecoMinMaxConsts(model, gen_dict, prod_vars, on_vars, reserve_vars):
    """
    on_var * eco_min <= prod_vars + reserve <= on_vars * eco_max
    """
    #VG Right now no on_vars for nonstandard generation.  Later update
    for name, iHr in on_vars.keys():
        sHr = "H" + str(iHr + 1)
        reserve = grb.quicksum(reserve_vars[name, iHr, type] 
                            for type in generator.GenUnit.RESERVE_PRODUCTS)
        if gen_dict[name].eco_max[sHr] < TOL:
            model.addConstr(on_vars[name, iHr] == 0 )
            model.addConstr(prod_vars[name, iHr] == 0 )
            model.addConstr(reserve == 0)
            continue
        model.addConstr( on_vars[name, iHr] * gen_dict[name].eco_min[sHr] <=
                                          prod_vars[name, iHr] + reserve, 
                                          name="EcoMin" + name + "H" + str(iHr) )  
        model.addConstr( prod_vars[name, iHr]  + reserve <= 
                                         on_vars[name, iHr] * gen_dict[name].eco_max[sHr], 
                                         name="Ecomax" + name + "H" + str(iHr) )
        #you need to be on in order to offer spinning reserve
        if (name, iHr, "TMSR_Cap") in reserve_vars:
            model.addConstr( reserve_vars[name, iHr, "TMSR_Cap"] <= 
                on_vars[name, iHr] * gen_dict[name].TMSR_Cap, 
                name="SpinningReserve" + "name" + sHr )

def __buildNomNoLoad(gen_dict, TMSR_REQ, T10_REQ, T30_REQ, includeIncDecs, useReserve, sparseRamps):
    """ Builds the nominal version of the capacity model
    Excludes the load balance constraints
    Returns UCNomModel obj"""
    m = grb.Model("UCNominal")
    on_vars, start_vars, stop_vars, cost_var = genStage1Vars( m, gen_dict )
    prod_vars, reserve_vars = genStage2VarsNom(m, gen_dict, useReserve )
    variable_cost_vars = addPiecewiseCosts(m, gen_dict, prod_vars )
    startStopConstraints(m, gen_dict, on_vars, start_vars, stop_vars)
    ecoMinMaxConsts(m, gen_dict, prod_vars, on_vars, reserve_vars)
    minUpConstraints(m, gen_dict, on_vars)
    minDownConstraints(m, gen_dict, on_vars)
    flex_loads = genFlexibleDemandVarsNom( m, gen_dict )

    if useReserve:
        reserve_cnsts = reserveRequirements(m, gen_dict, reserve_vars, TMSR_REQ, T10_REQ, T30_REQ)
        reserveCapacity(m, gen_dict, reserve_vars, reserve_cnsts)
    else:
        raise NotImplementedError()
    if includeIncDecs:
       dec_vars_amt, dec_vars_price = genDecVarsNom(m, gen_dict)
       inc_vars = genIncVarsNom( m, gen_dict)
    else:
        dec_vars_amt = dec_vars_price = inc_vars = {}

    #VG comment out now for consistency to affine.  
    rampingConsts(m, gen_dict, prod_vars, start_vars, stop_vars, sparseRamps)
    out = UCNomModel()
    out.model = m
    out.on_vars = on_vars
    out.start_vars = start_vars
    out.stop_vars = stop_vars
    out.fixed_cost_var = cost_var    
    out.prod_vars = prod_vars
    out.reserve_vars = reserve_vars
    out.variable_cost_vars = variable_cost_vars
    out.flex_loads = flex_loads
    out.dec_vars_amt = dec_vars_amt
    out.inc_vars = inc_vars
    return out    

def __addLoadBalanceCnst( nomUCObj, gen_dict, load_by_hr ):
    """ Adds a constraint to minimize the L1 deviation from the load."""
    #add a slack variable for the amount missed.    
    model = nomUCObj.model
    slack = [model.addVar(obj=PENALTY)  for ix in xrange(HORIZON_LENGTH) ]
    model.update()
    prod_by_hr = {}
    for name, hr in nomUCObj.prod_vars.keys():
        if hr not in prod_by_hr:
            prod_by_hr[hr] = nomUCObj.prod_vars[name, hr]
        else:
            prod_by_hr[hr] += nomUCObj.prod_vars[name, hr]

    for name, hr in nomUCObj.flex_loads.keys():
        prod_by_hr[hr] -= nomUCObj.flex_loads[name, hr]

    for name, hr in nomUCObj.inc_vars.keys():
        prod_by_hr[int(hr.lstrip("H")) - 1 ] += nomUCObj.inc_vars[name, hr]
    
    for name, hr in nomUCObj.dec_vars_amt.keys():
        prod_by_hr[int(hr.lstrip("H")) - 1] -= nomUCObj.dec_vars_amt[name, hr]

    balance_cnsts = []        
    for hr in xrange(HORIZON_LENGTH):
        balance_cnsts.append ( model.addConstr( prod_by_hr[hr] - load_by_hr[hr] <= slack[hr] ) )
        balance_cnsts.append ( model.addConstr( load_by_hr[hr] - prod_by_hr[hr] <= slack[hr] ) )
    return slack, balance_cnsts

def __summarizeSolution(UCObj, gen_dict, slacks):
    start_vars = UCObj.start_vars
    start_vals = {}
    for name, iHr in start_vars:
        start_vals[name, iHr] = start_vars[name, iHr].x
    on_vars = UCObj.on_vars
    on_vals = {}
    for name, iHr in start_vars:
        on_vals[name, iHr] = on_vars[name, iHr].x

    prod_by_hour = Counter()
    for name, hr in UCObj.prod_vars.keys():
        prod_by_hour[ "TOTAL", hr ] += UCObj.prod_vars[name, hr].x

    for name, hr in UCObj.prod_vars.keys():
        prod_by_hour[ gen_dict[name].fuel_type, hr ] += UCObj.prod_vars[name, hr].x
        
    for name, hr in UCObj.flex_loads:
        prod_by_hour["FLEX", hr] += UCObj.flex_loads[name, hr].x
        
    for name, hr in UCObj.inc_vars.keys():
        prod_by_hour["INC", int(hr.lstrip("H")) - 1] += UCObj.inc_vars[name, hr].x
    
    for name, hr, in UCObj.dec_vars_amt.keys():
        prod_by_hour["DEC", int(hr.lstrip("H")) -1 ] += UCObj.dec_vars_amt[name, hr].x    

    variable_costs = Counter()
    for ((name, hr), v) in UCObj.variable_cost_vars.items():
        variable_costs[hr] += v.x

    for hour in xrange(len(slacks)):
        prod_by_hour["Slack", hour] += slacks[hour].x

    #identify the three largest producers at hour 15 and just print to screen
    cap_producers= {}
    for name, hr in UCObj.on_vars.keys():
        if hr <> 15:
            continue
        cap_producers[name] = gen_dict[name].eco_max["H16"]
    t = sorted(cap_producers.items(), key = lambda (n, cap): cap, reverse=True)
    print "Biggest Produers:\t", t[:3]                

    return on_vals, start_vals, UCObj.fixed_cost_var.x, UCObj.model.objVal, prod_by_hour, variable_costs

def writeHourlyGens(file_out, label, prod_by_type):
    file_out.writerow([label, "Total"] + 
        [ prod_by_type["TOTAL", hr] for hr in xrange(HORIZON_LENGTH) ] )
    file_out.writerow([label, "Nuclear"] + 
        [ prod_by_type["Nuclear", hr] for hr in xrange(HORIZON_LENGTH) ] )
    file_out.writerow([label, "Hydro"] + 
        [ prod_by_type["Hydro", hr] for hr in xrange(HORIZON_LENGTH) ] )
    file_out.writerow([label, "Steam"] + 
        [ prod_by_type["Steam", hr] for hr in xrange(HORIZON_LENGTH) ] )
    file_out.writerow([label, "CT"] + 
        [ prod_by_type["CT", hr] for hr in xrange(HORIZON_LENGTH) ] )
    file_out.writerow([label, "Diesel"] + 
        [ prod_by_type["Diesel", hr] for hr in xrange(HORIZON_LENGTH) ] )
    file_out.writerow([label, "FixedImport"] + 
        [ prod_by_type["FixedImport", hr] for hr in xrange(HORIZON_LENGTH) ] )
    file_out.writerow([label, "Other"] + 
        [ prod_by_type["Other", hr] for hr in xrange(HORIZON_LENGTH) ] )
    file_out.writerow([label, "INC"] + 
        [ prod_by_type["INC", hr] for hr in xrange(HORIZON_LENGTH) ] )
    file_out.writerow([label, "DEC"] + 
        [ prod_by_type["DEC", hr] for hr in xrange(HORIZON_LENGTH) ] )
    file_out.writerow([label, "FLEX"] + 
        [ prod_by_type["FLEX", hr] for hr in xrange(HORIZON_LENGTH) ] )
    file_out.writerow([label, "FLEX"] + 
        [ prod_by_type["FLEX", hr] for hr in xrange(HORIZON_LENGTH) ] )
    file_out.writerow([label, "LOAD"] + 
        [ prod_by_type["LOAD", hr] for hr in xrange(HORIZON_LENGTH) ] )
    file_out.writerow([label, "Slack"] + 
        [ prod_by_type["Slack", hr] for hr in xrange(HORIZON_LENGTH) ] )


def writeHourlySchedCap(file_out, label, on_vals, gen_dict):
    #tally the amount of scheduled capacity by fuel_type
    cap_tally = Counter()
    for (name, iHr), val in on_vals.items():
        if round(val) > .5:
            cap_tally[iHr, gen_dict[name].fuel_type ] += gen_dict[name].eco_max["H%d" % (iHr + 1) ]

    fuel_types = set( fuel for (hr, fuel) in cap_tally )
    for fuel in fuel_types:
        file_out.writerow([label, fuel] + [cap_tally[hr, fuel] for hr in xrange(HORIZON_LENGTH) ] )        

############  the real workers
def buildSolveNom(gen_dict, TMSR_REQ, T10_REQ, T30_REQ, load_by_hr, 
                                    includeIncDecs = False, useReserve = True, sparseRamps = False):
    """Builds and solves nominal model from scratch.  Returns
    on_vals[name, hr], start_vals[name, hr], fixed cost, tot_cost, prod_by_hour[fuel_type, hr], variable_costs[hr]
    """
    nomUCObj = __buildNomNoLoad(gen_dict, TMSR_REQ, T10_REQ, T30_REQ, includeIncDecs, useReserve, sparseRamps)                               
    slacks, cnsts = __addLoadBalanceCnst( nomUCObj, gen_dict, load_by_hr)                            
    
    #add the callback if necessary
    if sparseRamps:
        nomUCObj.model.optimize( lambda model, where: nomUCObj.rampingCallback(gen_dict, where) )
    else:
        nomUCObj.model.optimize()
    
    return __summarizeSolution(nomUCObj, gen_dict, slacks)

def calcResReqs(gen_dict, load_by_hr):
    return buildSolveNom(gen_dict, 0, 0., 0., load_by_hr)    

#right now build second stage as an ip and fix everything
#check to see how much slower this is than a direct solution
def updateSolveSecondStage( UCObj, new_load_by_hr, old_objs, start_vals, on_vals, gen_dict ):
    model = UCObj.model
    for o in old_objs:
        model.remove( o )

    old_objs = []
    on_vars = UCObj.on_vars
    for k, v in on_vars.items():
        old_objs.append( model.addConstr( v == round(on_vals[k] )) )
    start_vars = UCObj.start_vars
    for k, v in start_vars.items():
        old_objs.append( model.addConstr( v == round(start_vals[k] )) )

    slacks, balance_cnsts = __addLoadBalanceCnst( UCObj, gen_dict, new_load_by_hr )
    old_objs += slacks
    old_objs += balance_cnsts
    model.optimize()

    return __summarizeSolution(UCObj, gen_dict, slacks), old_objs

def buildAffineModel(gen_dict, TMSR_REQ, T10_REQ, T30_REQ, resids, eps, delta, k, 
                                    includeIncDecs = False, sparseRamps = False, omitRamping=False):
    """builds and solves an affine modeling object.  Resulting object is reuseable iteration to iteration"""
    #build the uncertainty set and teh three functions
    uncSet = SparseAffineCutGen(resids, eps, 5, delta * .5, delta * .5)

    #now construct the model
    print "Beginning Affine Model:"
    m = grb.Model("Affine")
    on_vars, start_vars, stop_vars, cost_var = genStage1Vars( m, gen_dict )
    startStopConstraints(m, gen_dict, on_vars, start_vars, stop_vars)
    minUpConstraints(m, gen_dict, on_vars)
    minDownConstraints(m, gen_dict, on_vars)

    #Stage 2
    #build the norm variables for production because they're shared in many places
    prod_vars, reserve_vars = genStage2VarsAffine( m, gen_dict, k )

    #reserves
    reserveRequirementsAffine(m, reserve_vars, TMSR_REQ, T10_REQ, T30_REQ, uncSet)
    reserveCapacityAffine(m, gen_dict, reserve_vars, uncSet)

    #ramping
    if not omitRamping:
        raise ValueError()
#         rampingConstsAffine(m , gen_dict, prod_vars, start_vars, stop_vars,
#             addVars, addLess, addGreater, sparse=sparseRamps)

    #flex_loads
    flex_loads = genFlexibleDemandVarsAffine( m, gen_dict, k )
    flex_loads_norm = {}
    temp_vars = uncSet.addVecVars( m, len(flex_loads) )
    for ix, k in enumerate(flex_loads):
        flex_loads_norm[k] = uncSet.createNormVars(m, flex_loads[k][0], temp_vars[ix] )
    boundFlexDemandAffine( m, gen_dict, flex_loads, flex_loads_norm, uncSet)

    #eco Min Max
    ecoMinMaxConstsAffine(m, gen_dict, prod_vars, on_vars, reserve_vars, uncSet, M=1e4)
    #PiecewiseCosts
    variable_cost_vars = addPiecewiseCostsAffine(m, gen_dict, prod_vars, uncSet)
    #prep load balance
    slack, fvecs, fvecs_norm, gprod_sys = prepForLoadBalance(m, prod_vars, flex_loads, uncSet)

    return ( m, gen_dict, prod_vars, flex_loads, on_vars, start_vars, stop_vars, 
            cost_var, reserve_vars, variable_cost_vars, slack, fvecs, fvecs_norm, gprod_sys, uncSet)

def addSolveAffineLoadBalanceNaive(m, gen_dict, prod_vars, flex_loads, on_vars, start_vars, stop_vars, 
                                                                    fixed_cost_var, reserve_vars, avg_load_by_hr, variable_cost_vars, 
                                                                    old_objs, slack, fvecs, fvecs_norm, gprod_sys, uncSet):
    for o in old_objs:
        m.remove(o)
    slacks, balance_objs = affineLoadBalanceNaive(m, gen_dict, fvecs, fvecs_norm, gprod_sys, slack, avg_load_by_hr, uncSet) 

    #DEBUG ONLY
    m.update()
    m.printStats()

    m.optimize()
    return summarizeAffine(m, on_vars, start_vars, stop_vars, prod_vars, reserve_vars, variable_cost_vars, flex_loads, 
                                        fixed_cost_var, slacks), balance_objs    