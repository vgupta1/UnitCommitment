""" Nominal UC Model

Contains all the logic for constructing a nominal UC Model
Top level convenience functions for other tests
"""

import csv, pdb, numpy, sys
import gurobipy as grb
import generator, readData 
from AffineHelpers import *
from config import *

class UCNomModel():
    def __init__(self, model, on_vars, start_vars, stop_vars, prod_vars, reserve_vars, variable_cost_vars, flex_loads, 
                                        fixed_cost_var, dec_vars_amt, inc_vars ):
        self.model, self.on_vars, self.start_vars, self.stop_vars, self.prod_vars, self.reserve_vars = (
                model, on_vars, start_vars, stop_vars, prod_vars, reserve_vars )
        self.variable_cost_vars, self.flex_loads, self.fixed_cost_var, self.dec_vars_amt, self.inc_vars = (
                variable_cost_vars, flex_loads, fixed_cost_var, dec_vars_amt, inc_vars )
        self.slacks, self.balance_cnsts, self.fixed_var_cnsts = [], [], []

    def removeOldCnsts(self):
        for o in self.slacks:
            self.model.remove( o )
        for o in self.balance_cnsts:
            self.model.remove( o )
        for o in self.fixed_var_cnsts:
            self.model.remove(  o )

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

    def computeProdByType(self, gen_dict):
        prod_by_hour = Counter()
        for name, hr in self.prod_vars.keys():
            prod_by_hour[ "TOTAL", hr ] += self.prod_vars[name, hr].x
        for name, hr in self.prod_vars.keys():
            prod_by_hour[ gen_dict[name].fuel_type, hr ] += self.prod_vars[name, hr].x
        for name, hr in self.flex_loads:
            prod_by_hour["FLEX", hr] += self.flex_loads[name, hr].x
        for name, hr in self.inc_vars:
            prod_by_hour["INC", hr] += self.inc_vars[name, hr].x
        for name, hr, in self.dec_vars_amt:
            prod_by_hour["DEC", hr ] += self.dec_vars_amt[name, hr].x    
        for hour in xrange(HORIZON_LENGTH):
            prod_by_hour["Slack", hour] += self.slacks[hour].x
        return prod_by_hour

    def summarizeSolution(self, gen_dict):
        start_vals = {}
        for (name, iHr) in self.start_vars:
            start_vals[name, iHr] = self.start_vars[name, iHr].x
        on_vals = {}
        for name, iHr in self.on_vars:
            on_vals[name, iHr] = self.on_vars[name, iHr].x
        variable_costs = Counter()
        for ((name, hr), v) in self.variable_cost_vars.items():
            variable_costs[hr] += v.x
        prod_by_hour = self.computeProdByType(gen_dict)
        return on_vals, start_vals, self.fixed_cost_var.x, self.model.objVal, prod_by_hour, variable_costs

def __buildNomNoLoad(gen_dict, TMSR_REQ, T10_REQ, T30_REQ, includeIncDecs, sparseRamps):
    """ Build everything except load-balance constraints
    Returns UCNomModel obj with various properties """
    m = grb.Model("UCNominal")
    true_gens = dict( filter( lambda (n, g): g.res_type == "GEN" and g.fuel_type <> "FixedImport", gen_dict.items() ) )
    on_vars, start_vars, stop_vars, cost_var = genStage1Vars( m, true_gens )
    prod_vars, reserve_vars = genStage2VarsNom(m, true_gens )
    variable_cost_vars = addPiecewiseCosts(m, true_gens, prod_vars )
    startStopConstraints(m, true_gens, on_vars, start_vars, stop_vars)
    ecoMinMaxConsts(m, gen_dict, prod_vars, on_vars, reserve_vars)
    minUpConstraints(m, true_gens, on_vars)
    minDownConstraints(m, true_gens, on_vars)
    flex_loads = genFlexibleDemandVarsNom( m, gen_dict )

    reserveRequirements(m, reserve_vars, TMSR_REQ, T10_REQ, T30_REQ)
    reserveCapacity(m, true_gens, reserve_vars)
    if includeIncDecs:
       dec_vars_amt, dec_vars_price = genDecVarsNom(m, gen_dict)
       inc_vars = genIncVarsNom( m, gen_dict)
    else:
        dec_vars_amt = dec_vars_price = inc_vars = {}
    rampingConsts(m, true_gens, prod_vars, start_vars, stop_vars, sparseRamps)

    return UCNomModel( m, on_vars, start_vars, stop_vars, prod_vars, reserve_vars, 
                                            variable_cost_vars, flex_loads, cost_var, dec_vars_amt, inc_vars )

def buildSolveNom(gen_dict, TMSR_REQ, T10_REQ, T30_REQ, load_by_hr, 
                                    includeIncDecs=False, sparseRamps=False):
    """Builds and solves nominal model from scratch.  Returns
    on_vals[name, hr], start_vals[name, hr], fixed cost, tot_cost, prod_by_hour[fuel_type, hr], variable_costs[hr]
    """
    nomUCObj = __buildNomNoLoad(gen_dict, TMSR_REQ, T10_REQ, T30_REQ, includeIncDecs, sparseRamps)                               
    __addLoadBalanceCnst( nomUCObj, load_by_hr)                            

    #VG DEBUG
    nomUCObj.model.params.mipgap = 5e-3

    if sparseRamps:
        nomUCObj.model.optimize( lambda model, where: nomUCObj.rampingCallback(gen_dict, where) )
    else:
        nomUCObj.model.optimize()
    return nomUCObj.summarizeSolution(gen_dict)

def updateSolveSecondStage( UCObj, new_load_by_hr, gen_dict, on_vals, start_vals, forceBalance=False ):
    UCObj.removeOldCnsts()
    model = UCObj.model
    old_objs = []
    for key, val in on_vals.items() :
        old_objs.append( UCObj.model.addConstr( UCObj.on_vars[key] == round( val ), 
                                                                                        name="FixOnVar%s%d" % key) )
    for key, val in start_vals.items() :
        old_objs.append( model.addConstr( UCObj.start_vars[key] == round( val ), 
                                                                                        name="FixStartVal%s%d" % key) )

    UCObj.fixed_var_cnsts = old_objs
    __addLoadBalanceCnst( UCObj, new_load_by_hr, forceBalance )
    model.optimize()
    return UCObj.summarizeSolution(gen_dict)

    

#############
#  The following functions are all helpers that do various tasks
##############
def genStage1Vars( model, true_gens, init_on={}, init_start={} ):
    """ Creates all stage 1 (committment) variables
    returns dicts:  on[gen, time], start[gen, time], stop[gen, time], cost_var
    cost_var is an auxiliary variable that allows us to access the first stage costs directly
    """
    on_vars, start_vars, stop_vars = {}, {}, {}
    cost_var_expr = grb.LinExpr()
    for name, gen in true_gens.items():
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

    #MIP Start Info
    if init_on:
        for name, gen in true_gens.items():
            for iHr in xrange(HORIZON_LENGTH ):
                on_vars[name, iHr].start = init_on[name, iHr]
                start_vars[name, iHr].start = init_start[name, iHr]                    

    return on_vars, start_vars, stop_vars, cost_var
    
def addPiecewiseCosts(model, true_gens, prod_vars):
    """Creates variables and the defining constraints for the piecewise cost functions
    for the variable cost """
    #First generate all the variables so only need one call to model.update()
    ys_all, cost_vars = {}, {}
    for name, gen in true_gens.items():
        #There is a nicety that all the true generators pertain H1*H36
        assert gen.offerBlocks.keys() == ['H1*H36']
        offerBlocks = gen.offerBlocks['H1*H36'] 
        #For numerical stability, makes sense to clip these at eco_max.
        eco_max = max( gen.eco_max.values() )
        for numKnots, b in enumerate(offerBlocks):
            if b.size >= eco_max:
                break
        numKnots += 1
        for iHr in xrange(HORIZON_LENGTH):
            ys_all[name, iHr] = [model.addVar(ub=1.0) for ix in xrange(numKnots ) ] 
            cost_vars[name, iHr] = model.addVar(obj=1.0)
    model.update()

    #now go back and add relevant constraints against the variables
    for name, gen in true_gens.items():
        offerBlocks = gen.offerBlocks['H1*H36'] 
        #clip these at eco-max.
        #VG in reality, could clip them also at eco min and do this by hour.
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
    
def startStopConstraints(m, true_gens, on_vars, start_vars, stop_vars):
    for name, gen in true_gens.items():
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
            if gen.eco_max[iHr] < TOL:
                m.addConstr( start_vars[name, iHr] == 0 )
    return m            

def minUpConstraints(model, true_gens, on_vars):
    for name, gen in true_gens.items():
        if gen.min_up <= 1:
            continue
        #Treat the first time slice specially
        #if its been on for less than min-run, keep it on.
        if 0 <= gen.init_online < gen.min_up: 
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

def minDownConstraints(model, true_gens, on_vars):
    for name, gen in true_gens.items():
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
    dec_vars_amt, dec_vars_cost, ys_all = {}, {}, {}
    #iterate twice so only one model update
    for name, gen in gen_dict.items():
        if gen.res_type <> "DEC":
            continue
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
    
def genStage2VarsNom(model, true_gens):
    """Creates the continuous variables for nominal model.
    prod[gen, time], reserve[gen, time, type]
    """
    prod, reserve = {}, {}
    for name, gen in true_gens.items():
        for iHr in range(HORIZON_LENGTH):
            prod[name, iHr] = model.addVar(name="Prod%sH%d" % (name, iHr))
            for cap_type in generator.GenUnit.RESERVE_PRODUCTS:
                    reserve[name, iHr, cap_type] = model.addVar(name="Res%sH%d%s" % (name iHr, cap_type))
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
            flex_loads[name, iHr] = model.addVar(lb=gen.eco_min[iHr ], 
                                                                                 ub = gen.eco_max[iHr], 
                                                                                 obj=-gen.offerBlocks.values()[0][0].price, 
                                                                                 name="FlexLoad" + name + "H" + str(iHr) )
    return flex_loads

def reserveRequirements(model, reserve_vars, TMSR_REQ, T10_REQ, T30_REQ):
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

def reserveCapacity(model, true_gens, reserve_vars):
    """ No generator can offer more reserve of any type than its capacity allows."""
    for name, gen in true_gens.items():
        #TMSR cap already handled by the on-off
        for iHr in xrange(HORIZON_LENGTH):
            model.addConstr( reserve_vars[name, iHr, "TMSR_Cap"] + 
                                                reserve_vars[name, iHr, "TMNSR_Cap"] <= gen.T10_Cap )
        for iHr in xrange(HORIZON_LENGTH):
            model.addConstr( reserve_vars[name, iHr, "TMSR_Cap"] + 
                                                reserve_vars[name, iHr, "TMNSR_Cap"] +
                                                reserve_vars[name, iHr, "TMOR_Cap"] <= gen.T30_Cap )

def rampingConsts(model, true_gens, prod_vars, start_vars, stop_vars, sparse):
    for name, gen in true_gens.items():
        if sparse and gen.fuel_type not in ("Steam", "CT", "CC", "Diesel"): #only add ramping for marginal turbines
            continue
        #ignore ramping constraints for the first slice.
        for hr in xrange(1, HORIZON_LENGTH):
            ##VG This logic is partially duplicated in "fixRampRates of generator.py"
            eco_min_m = gen.eco_min[hr - 1]
            eco_max_m = gen.eco_max[hr - 1] 
            eco_min = gen.eco_min[hr]
            eco_max = gen.eco_max[hr]
            if eco_max <= TOL: #won't run anyway
                continue
            if gen.ramp_rate < eco_max_m - eco_min:
                model.addConstr( prod_vars[name, hr - 1] - prod_vars[name, hr] <= 
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
    #No on_vars for nonstandard generation.
    for name, iHr in on_vars.keys():
        reserve = grb.quicksum(reserve_vars[name, iHr, type] 
                            for type in generator.GenUnit.RESERVE_PRODUCTS)
        if gen_dict[name].eco_max[iHr] < TOL:
            model.addConstr(on_vars[name, iHr] == 0 )
            model.addConstr(prod_vars[name, iHr] == 0 )
            model.addConstr(reserve == 0)
            continue
        model.addConstr( on_vars[name, iHr] * gen_dict[name].eco_min[iHr] <=
                                          prod_vars[name, iHr] + reserve, 
                                          name="EcoMin" + name + "H" + str(iHr) )  
        model.addConstr( prod_vars[name, iHr]  + reserve <= 
                                         on_vars[name, iHr] * gen_dict[name].eco_max[iHr], 
                                         name="Ecomax" + name + "H" + str(iHr) )
        #you need to be on in order to offer spinning reserve
#        if (name, iHr, "TMSR_Cap") in reserve_vars:  VG eliminate
        model.addConstr( reserve_vars[name, iHr, "TMSR_Cap"] <= 
            on_vars[name, iHr] * gen_dict[name].TMSR_Cap, 
            name="SpinningReserve%sH%d" % (name, iHr))
            
def __addLoadBalanceCnst(nomUCObj, load_by_hr, forceBalance=False):
    """ Adds a constraint to minimize the L1 deviation from the load."""
    #Assume that old objects have already been removed
    model = nomUCObj.model

    #add a slack variable for the amount missed.    
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
        prod_by_hr[hr] += nomUCObj.inc_vars[name, hr]
    
    for name, hr in nomUCObj.dec_vars_amt.keys():
        prod_by_hr[hr] -= nomUCObj.dec_vars_amt[name, hr]

    balance_cnsts = []        
    for hr in xrange(HORIZON_LENGTH):
        balance_cnsts.append ( model.addConstr( prod_by_hr[hr] - load_by_hr[hr] <= slack[hr], 
                                                                                     name = "Balance2%d" % hr ) )
        balance_cnsts.append ( model.addConstr( load_by_hr[hr] - prod_by_hr[hr] <= slack[hr], 
                                                                                     name="Balance1%d" % hr ) )

    nomUCObj.slacks = slack
    nomUCObj.balance_cnsts = balance_cnsts


