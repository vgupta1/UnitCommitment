""" Build the Unit Committment Problem

Sequence of helper functions used in creating a unit committment instance
"""
import gurobipy as grb
import csv, pdb, numpy, sys
import generator, readData
from collections import *
        
HORIZON_LENGTH = 24
EPSILON_ZERO = 1
PENALTY = 5000
TOL = 1e-2

TMSR_REQ = 622.5
T10_REQ  = 1245.
T30_REQ = 1883.

## To DO:
# Add cold/warm/hot starts
# Only add ramping constraints for "interesting" Steam turbies.  Treat rest as cut generation + lazy
# incorporate startup times correctly
# Go through and find various places we can do some variable elimination:
  # gens with identical eco-min eco-max
  # gens with zero eco_max in particular hours


#notes:
#The "GEN" with "FixedImport" have no generator characteristics, but have fixedEnergy
#amounts in the hours after 24....  hence we just drop them from the data set.
#We are ignoring network structure for now
#For the incremental bid, actually vary by hour, but single block cost structure

def hrToInt( sHr ):
    """ Convert "H12" to 12 """
    return int( sHr.lstrip("H") )

def genStage1Vars( model, gen_dict ):
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

    return on_vars, start_vars, stop_vars, cost_var

def genStage2VarsNom( model, gen_dict, useReserve=False ):
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

def genStage2VarsAffine( model, gen_dict ):
    """Includes reserve varibles by default
    prod[gen,time] = (f vec, g consant), reserve[gen, time, type] = (fvec, gconstant)
    """
    #VG Experiment computationally with value of adding explicit upper bounds to these variables
    prod, reserve = {}, {}
    for name, gen in gen_dict.items():
        if gen.res_type <> "GEN" or gen.fuel_type == "FixedImport":
            continue
        for iHr in xrange(HORIZON_LENGTH):
            fvec = [model.addVar(lb = -grb.GRB.INFINITY) for ix in xrange(HORIZON_LENGTH) ] 
            prod[name, iHr] = (fvec, model.addVar(lb=-grb.GRB.INFINITY) )    
        
            for cap_type in generator.GenUnit.RESERVE_PRODUCTS:
                fvec = [model.addVar(lb = -grb.GRB.INFINITY) for ix in xrange(HORIZON_LENGTH) ] 
                reserve[name, iHr, cap_type] = (fvec, model.addVar(lb=-grb.GRB.INFINITY))    
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

def genFlexibleDemandVarsAffine( model, gen_dict ):
    """Creates Affine variables for flex demands.  Interpret these as revenue earned and extra load
    to satisfy 
    flex_loads[name, iHr] = (fvec, gconstant) 
    This formulation requires we call a separate function to add the lower and upper bounds and
    associated cost
    """
    flex_loads = {}
    for name, gen in gen_dict.items():
        if gen.res_type <> "LOAD" or not gen.isFlexible:
            continue
        
        #these flex loads have a single block bid structure
        for iHr in range(HORIZON_LENGTH):
            fvec = [model.addVar(lb = -grb.GRB.INFINITY) for ix in xrange(HORIZON_LENGTH) ]             
            flex_loads[name, iHr] = (fvec, model.addVar(lb= -grb.GRB.INFINITY) )
    return flex_loads

def boundFlexDemandAffine( model, gen_dict, flex_loads, avg_load, kappa, gam1, C):
    """Adds upper and lower bounds and the cost variable for the affine formulation for flex."""
    #requisition everything upfront
    revenue_vars = {}
    for (name, iHr), (fvec, g0) in flex_loads.items():
        revenue_vars[name, iHr] = model.addVar(obj = -gen_dict[name].offerBlocks.values()[0][0].price) #Notice neg obj 
    t1s, t2s, t3s, t4s = __addVecVars(model, 3 * len(flex_loads) )
    ix = 0
    for (name, iHr), (fvec, g0) in flex_loads.items():
        #lower bound
        __addGreaterEqual( model, fvec, g0, gen_dict[name].eco_min["H" + str(iHr + 1)], 
                                            avg_load, kappa, gam1, C, "FlexDemandLB" + name + str(iHr), 
                                            t1s[ix], t2s[ix], t3s[ix], t4s[ix])        

        #upper bound
        __addLessEqual( model, fvec, g0, gen_dict[name].eco_max["H" + str(iHr + 1)], 
                                        avg_load, kappa, gam1, C, "FlexDemandUB" + name + str(iHr), 
                                        t1s[ix + 1], t2s[ix + 1], t3s[ix + 1], t4s[ix + 1] )
        #Revenue
        __addGreaterEqual( model, fvec, g0, revenue_vars[name, iHr], 
                                                avg_load, kappa, gam1, C, "RevenueFlex" + name + str(iHr), 
                                                t1s[ix + 2], t2s[ix+2], t3s[ix + 2], t4s[ix + 2] )
        ix +=3

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

def __addVecVars( model, numConstrs):
    """ Adds all the variables for the vec version of Less Equal.  a tuple of variable handles"""
    t1s, t2s, t3s, t4s = [], [], [], []
    for ix in xrange(numConstrs):
        t1s.append(model.addVar() )
        t2s.append(model.addVar() )
        t3s.append( [model.addVar(lb = -grb.GRB.INFINITY) for ix in xrange(HORIZON_LENGTH)] )
        t4s.append( [model.addVar(lb = -grb.GRB.INFINITY) for ix in xrange(HORIZON_LENGTH)] )

    model.update()
    return t1s, t2s, t3s, t4s

def __addLessEqual( model, fvec, g0, rhs, avg_load, kappa, gam1, C, sname, 
                                        t1 = None, t2 = None, t3 = None, t4 = None):
    """Adds the equivalent to the robust affine constraint f0' d + g <= rhs for all d in U """
    if t1 is None:
        t1 = model.addVar()
        t2 = model.addVar()
        t3 = [ model.addVar( lb = -grb.GRB.INFINITY) for ix in xrange(HORIZON_LENGTH) ] 
        t4 = [ model.addVar(lb = -grb.GRB.INFINITY) for ix in xrange(HORIZON_LENGTH) ]
        model.update()
    
    objs = [ t1, t2]
    objs += t3
    objs += t4

    objs.append( model.addConstr(grb.quicksum(f * m for f, m in zip(fvec, avg_load) ) + 
                                    gam1 * t1 + kappa * t2 + g0 <= rhs, name=sname )  )
    for row, t in zip(C, t3):
        objs.append( model.addConstr( grb.quicksum( r * f for r,f in zip(row, fvec) ) == t ) )

    for f, t in zip(fvec, t4):
        objs.append( model.addConstr( f == t ) )
    objs.append( model.addQConstr( grb.quicksum( t * t for t in t4 ) <= t1*t1, name=sname + "Q1" ) )
    objs.append( model.addQConstr( grb.quicksum( t * t for t in t3 ) <= t2 * t2, name=sname + "Q2" ) )

    return objs

def __addGreaterEqual( model, fvec, g0, rhs, avg_load, kappa, gam1, C, sname,
                                        t1 = None, t2 = None, t3 = None, t4 = None):
    if t1 is None:
        t1 = model.addVar()
        t2 = model.addVar()
        t3 = [ model.addVar( lb = -grb.GRB.INFINITY) for ix in xrange(HORIZON_LENGTH) ] 
        t4 = [ model.addVar(lb = -grb.GRB.INFINITY) for ix in xrange(HORIZON_LENGTH) ]
        model.update()
    objs = [ t1, t2]
    objs += t3
    objs += t4
    objs.append( model.addConstr(grb.quicksum(f * m for f, m in zip(fvec, avg_load) ) - 
                                    gam1 * t1 - kappa * t2 + g0 >= rhs, name=sname ) )
    for row, t in zip(C, t3):
        objs.append( model.addConstr( numpy.dot( row, fvec ) == t ) )              
    for f, t in zip(fvec, t4):
        objs.append( model.addConstr( f == t ) )

    objs.append( model.addQConstr( numpy.dot(t4, t4) <= t1*t1, name = sname + "Q1" ) )
    objs.append( model.addQConstr( numpy.dot(t3, t3) <= t2 * t2, name=sname + "Q2" ) )
    return objs
   
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
def addPiecewiseCosts(model, gen_dict, prod_vars, useAffine=False, Udict = {} ):
    """Creates variables and the defining constraints for the piecewise cost functions
    for the variable cost """
    ys_all = {}
    cost_vars = {}
    for name, gen in gen_dict.items():
        if gen.res_type <> "GEN" or gen.fuel_type == "FixedImport":
            continue
        #There is a nicety that all the true generators pertain H1*H36
        assert gen.offerBlocks.keys() == ['H1*H36']
        offerBlocks = gen.offerBlocks['H1*H36'] 
        #For numerical stability, makes sense to clip these at eco_max.
        eco_max = max( gen.eco_max.values() ) + TOL
        for ix, b in enumerate(offerBlocks):
            if b.size > eco_max:
                break
        sizes = [b.size for b in offerBlocks[:ix] ] + [ eco_max ] 
        prices = [b.price for b in offerBlocks[:ix] ] + [ offerBlocks[ix].price ] 
        size_diffs = [s  - sm for s, sm in zip(sizes, [0] + sizes) ] 
        numKnots = len(size_diffs)
        assert numKnots >= 1

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
        eco_max = max( gen.eco_max.values() ) + TOL
        for ix, b in enumerate(offerBlocks):
            if b.size > eco_max:
                break
        sizes = [b.size for b in offerBlocks[:ix] ] + [ eco_max ] 
        prices = [b.price for b in offerBlocks[:ix] ] + [ offerBlocks[ix].price ] 
        size_diffs = [s  - sm for s, sm in zip(sizes, [0] + sizes) ] 
        numKnots = len(sizes)
        assert numKnots >= 1
    
        for iHr in xrange(HORIZON_LENGTH):     
            model.addConstr( cost_vars[name, iHr] == 
                                grb.quicksum( y * p  * s for y, p, s  in zip(ys_all[name, iHr], prices, size_diffs) ) )    
            if useAffine:
                f, g = prod_vars[name, iHr]
                __addLessEqual(model, f, g, grb.quicksum( y * s for y, s in zip(ys_all[name, iHr], size_diffs) ), 
                                                Udict["avg_load"], Udict["kappa"], Udict["gam1"], Udict["C"], 
                                                "PieceWiseCost" + name + str(iHr))
            else:
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
            m.addConstr( stop_vars[name, iHr] <= 1 - on_vars[name, iHr] )

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

def nominalLoadBalance(model, gen_dict, prod_vars, load_by_hr, 
                                              flex_loads, dec_vars_amt, inc_vars):
    """ Adds a constraint to minimize the L1 deviation from the load."""
    #add a slack variable for the amount missed.    
    slack = [model.addVar(obj=PENALTY)  for ix in xrange(HORIZON_LENGTH) ]
    model.update()
    prod_by_hr = {}
    for name, hr in prod_vars.keys():
        if hr not in prod_by_hr:
            prod_by_hr[hr] = prod_vars[name, hr]
        else:
            prod_by_hr[hr] += prod_vars[name, hr]

    for name, hr in flex_loads.keys():
        prod_by_hr[hr] -= flex_loads[name, hr]

    for name, hr in inc_vars.keys():
        prod_by_hr[int(hr.lstrip("H")) - 1 ] += inc_vars[name, hr]
    
    for name, hr in dec_vars_amt.keys():
        prod_by_hr[int(hr.lstrip("H")) - 1] -= dec_vars_amt[name, hr]

    balance_cnsts = []        
    for hr in xrange(HORIZON_LENGTH):
        balance_cnsts.append ( model.addConstr( prod_by_hr[hr] - load_by_hr[hr] <= slack[hr] ) )
        balance_cnsts.append ( model.addConstr( load_by_hr[hr] - prod_by_hr[hr] <= slack[hr] ) )
    return slack, balance_cnsts

def affineLoadBalanceNaive(model, gen_dict, fprod_sys, gprod_sys, 
                                                    avg_load, kappa, gam1, C):
    """ Adds a constraint to minimize the L1 deviation from the load."""
    #add a slack variable for the amount missed.    
    slack = [model.addVar(obj=PENALTY)  for ix in xrange(HORIZON_LENGTH)]
    aux_fs = [model.addVar(lb = -grb.GRB.INFINITY ) for ix in xrange(HORIZON_LENGTH)]
    t1s, t2s, t3s, t4s = __addVecVars(model, 2 * HORIZON_LENGTH)

    balance_cnsts = []        
    for hr in xrange(HORIZON_LENGTH):
        balance_cnsts.append( model.addConstr( aux_fs[hr] == fprod_sys[hr] - 1 ) )
        fvec = [f for f in fprod_sys[hr] ]
        fvec[hr] = aux_fs[hr]
        balance_cnsts += __addLessEqual( model, fvec, gprod_sys[hr] , slack[hr], avg_load, kappa, gam1, C, 
            "LoadBalanceGB" + str(hr), 
            t1s[2 * hr], t2s[2 * hr], t3s[2 * hr], t4s[2 *hr] )
        balance_cnsts += __addGreaterEqual( model, fvec, gprod_sys[hr], -slack[hr], avg_load, kappa, gam1, C, 
            "LoadBalanceLB" + str(hr), 
            t1s[2* hr + 1], t2s[2 * hr + 1], t3s[2*hr + 1], t4s[2 * hr + 1])
    return slack, balance_cnsts

def formSysLevelAffine( prod_vars, flex_vars ):
    f_sys, g_sys = {}, {}
    for name, hr in prod_vars.keys():
        if hr not in f_sys:
            f_sys[hr] = numpy.zeros(HORIZON_LENGTH)
            g_sys[hr] = 0.
        f, g = prod_vars[name, hr]
        f_sys[hr] =  f_sys[hr] + f
        g_sys[hr] += g

    for name, hr in flex_loads.keys():
        f, g = flex_loads[ name, hr]
        f_sys[hr] = f_sys[hr] - f
        g_sys[hr] -= g

    return f_sys, g_sys

def rampingConsts(model, gen_dict, prod_vars, start_vars, stop_vars, M=5e3):
    for name, gen in gen_dict.items():
        if gen.res_type <> "GEN" or gen.fuel_type == "FixedImport":
            continue

        if gen.ramp_rate is None:
            print "No Ramp Rate:\t", name
            continue
        
        #ignore ramping constraints for the first slice.
        for hr in xrange(1, HORIZON_LENGTH):
            model.addConstr( prod_vars[name, hr] - prod_vars[name, hr - 1] <= 
                                                gen.ramp_rate + M * start_vars[name, hr], 
                                                name= "RampUP" + name + "H" + str(hr) )
            model.addConstr( prod_vars[name, hr -1] - prod_vars[name, hr] <= 
                                                gen.ramp_rate + M * stop_vars[name, hr], 
                                                name="RampDown" + name + "H" + str(hr) )
    return model    

def rampingConstsAffine(model, gen_dict, prod_vars, start_vars, stop_vars, 
                                                avg_load, kappa, gam1, C, M=5e3, sparse=False):
    for name, gen in gen_dict.items():
        if gen.res_type <> "GEN" or gen.fuel_type == "FixedImport":
            continue
        if sparse and gen.fuel_type not in ("Steam", "CT"): #only add ramping for marginal turbines
            continue

        #ignore ramping constraints for the first slice.
        for hr in xrange(1, HORIZON_LENGTH):
            ##VG This logic is partially duplicated in "fixRampRates of generator.py"
            # Combine these two.
                eco_min_m = gen.eco_min["H" + str(hr)]
                eco_max_m = gen.eco_max["H" + str(hr)] 
                eco_min = gen.eco_min["H" + str(hr + 1)]
                eco_max = gen.eco_max["H" + str(hr + 1)]
                if eco_max <= TOL: #won't run anyway
                    continue

                f, g = prod_vars[name, hr]
                fm, gm = prod_vars[name, hr - 1]

                #check the ramping down constraint could be tight
                if gen.ramp_rate < eco_max_m - eco_min:
                    __addLessEqual(model, numpy.array(fm) - numpy.array(f), gm - g, 
                                                gen.ramp_rate + eco_max_m * stop_vars[name, hr], 
                                                avg_load, kappa, gam1, C, 
                                                "rampingUB" + name + str(hr) )
                                
                #check the ramping up constraint could be tight
                if gen.ramp_rate < eco_max - eco_min_m:        
                    __addLessEqual(model, numpy.array(f) - numpy.array(fm), g - gm, 
                                                gen.ramp_rate + eco_max * start_vars[name, hr], 
                                                avg_load, kappa, gam1, C, 
                                                "rampingLB" + name + str(hr))
            ## End duplication
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
        eco_max_dict = gen_dict[name].eco_max
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
        reserve = grb.quicksum(reserve_vars[name, iHr, type] 
                            for type in generator.GenUnit.RESERVE_PRODUCTS)

def ecoMinMaxConstsAffine(model, gen_dict, prod_vars, on_vars, reserve_vars, avg_load, 
kappa, gam1, C):
    """
    on_var * eco_min <= prod_vars + reserve <= on_vars * eco_max
    """
    #iterate through setting the eco_max zero elements, and computing how many constraints we need
    numConstr = 0
    for name, iHr in on_vars.keys():
        #special case for eco_max == 0
        sHr = "H" + str(iHr + 1)
        if gen_dict[name].eco_max[sHr] < TOL:
            fprod, gprod = prod_vars[name, iHr]
            model.addConstr(on_vars[name, iHr] == 0 )
            model.addConstr(gprod == 0)
            for f in fprod:
                model.addConstr(f == 0)
            for type in generator.GenUnit.RESERVE_PRODUCTS:
                model.addConstr( reserve_vars[name, iHr, type][1] == 0 )
                for f in reserve_vars[name, iHr, type][0]:
                    model.addConstr(f == 0)
        else:
            numConstr +=3            

    #requisition all the variables in one go
    t1s, t2s, t3s, t4s =  __addVecVars(model, numConstr) 
    ix = 0
    for name, iHr in on_vars.keys():
        sHr = "H" + str(iHr + 1)
        if gen_dict[name].eco_max[sHr] >= TOL:        
            fprod, gprod = prod_vars[name, iHr]
            fres = numpy.zeros( HORIZON_LENGTH )
            gres = 0.
            for type in generator.GenUnit.RESERVE_PRODUCTS:
                fres = fres + reserve_vars[name, iHr, type][0]
                gres += reserve_vars[name, iHr, type][1]

            __addGreaterEqual(model, numpy.array(fprod) + fres, gprod + gres, 
                    on_vars[name, iHr] * gen_dict[name].eco_min[sHr], 
                    avg_load, kappa, gam1, C, 
                    "Ecomin" + name + str(iHr), 
                    t1s[ix], t2s[ix], t3s[ix], t4s[ix] )
            __addLessEqual(model, numpy.array(fprod) + fres, gprod + gres,  
                    on_vars[name, iHr] * gen_dict[name].eco_max[sHr], 
                   avg_load, kappa, gam1, C, 
                   "Ecomax" + name + str(iHr), 
                   t1s[ix + 1], t2s[ix + 1], t3s[ix + 1], t4s[ix + 1])
    
            #you need to be on in order to offer spinning reserve
            if (name, iHr, "TMSR_Cap") in reserve_vars:
                __addLessEqual( model, reserve_vars[name, iHr, "TMSR_Cap"][0], 
                        reserve_vars[name, iHr, "TMSR_Cap"][1], 
                        on_vars[name, iHr] * gen_dict[name].TMSR_Cap, 
                        avg_load, kappa, gam1, C, 
                        "OnForSpin" + name + str(iHr), 
                        t1s[ix + 2], t2s[ix + 2], t3s[ix + 2], t4s[ix + 2])
            else:
                print "No TMSR_Cap:\t", name, iHr
            ix += 3     

def reserveRequirements(model, gen_dit, reserve_vars):
    cnsts = []
    for iHr in xrange(HORIZON_LENGTH):
        TMSR_vars_by_hr = filter( lambda (name, hr, type): hr == iHr and type == "TMSR_Cap", 
                                                        reserve_vars.keys() )
        cnsts.append ( model.addConstr( 
                grb.quicksum(reserve_vars[k] for k in TMSR_vars_by_hr ) >= TMSR_REQ ) )
        T10_vars_by_hr = filter( lambda (name, hr, type): hr == iHr and 
                                                    type in ("TMSR_CAP", "TMNSR_Cap"), 
                                                    reserve_vars.keys() )
        cnsts.append( model.addConstr( 
                grb.quicksum(reserve_vars[k] for k in T10_vars_by_hr ) >= T10_REQ, 
                name="T10_Req" + "H" + str(iHr) ) )
        T30_vars_by_hr = filter( lambda (name, hr, type): hr == iHr and 
                                                    type in ("TMSR_CAP", "TMNSR_Cap", "TMOR_CAP"), 
                                                    reserve_vars.keys())
        cnsts.append ( model.addConstr( 
                grb.quicksum(reserve_vars[k] for k in T30_vars_by_hr ) >= T30_REQ, 
                                name="T10_Req" + "H" + str(iHr) ) )
    return cnsts

def reserveRequirementsAffine(model, reserve_vars, Udict):
    cnsts = []
    t1s, t2s, t3s, t4s = __addVecVars(model, 3 * HORIZON_LENGTH)
    for iHr in xrange(HORIZON_LENGTH):
        TMSR_vars_by_hr = ( v for ((name, hr, type), v) in reserve_vars.items() if 
                                                        hr == iHr and type == "TMSR_Cap" )
        TMSR_fs, TMSR_gs = zip( * TMSR_vars_by_hr )
        f = reduce( lambda x, y: numpy.add(x, y), TMSR_fs )
        g = reduce( lambda x, y: x + y, TMSR_gs )
        cnsts += __addGreaterEqual(model, f, g, TMSR_REQ, Udict["avg_load"], Udict["kappa"], 
                Udict["gam1"], Udict["C"], "TMSRReq" + str(iHr), 
                t1s[3 * iHr], t2s[3 *iHr], t3s[ 3* iHr], t4s[ 3* iHr] )

        T10_vars_by_hr = ( v for ((name, hr, type), v) in reserve_vars.items() if 
                                                        hr == iHr and type in ("TMSR_CAP", "TMNSR_Cap") )
        T10_fs, T10_gs = zip(* T10_vars_by_hr )            
        f = reduce( lambda x, y: numpy.add(numpy.array(x), numpy.array(y) ), T10_fs )
        g = reduce( lambda x, y: x + y, T10_gs )
        cnsts += __addGreaterEqual(model, f, g, T10_REQ,Udict["avg_load"], Udict["kappa"], 
                Udict["gam1"], Udict["C"], "T10Req" + str(iHr), 
                t1s[3 * iHr + 1], t2s[3 *iHr + 1], t3s[3* iHr + 1], t4s[3* iHr + 1] )

        T30_vars_by_hr = ( v for ((name, hr, type), v) in reserve_vars.items() if 
                                                        hr == iHr and type in ("TMSR_CAP", "TMNSR_Cap", "TMOR_CAP") )
        T30_fs, T30_gs = zip(* T30_vars_by_hr )            
        f = reduce( lambda x, y: numpy.add(numpy.array(x), numpy.array(y) ), T30_fs )
        g = reduce( lambda x, y: x + y, T30_gs )
        cnsts += __addGreaterEqual(model, f, g, T30_REQ, Udict["avg_load"], Udict["kappa"], 
                Udict["gam1"], Udict["C"], "T30Req" + str(iHr), 
                t1s[3 * iHr + 2], t2s[3 *iHr + 2], t3s[3* iHr + 2], t4s[3* iHr + 2] )
    return cnsts

def reserveCapacityAffine(model, gen_dict, reserve_vars, res_cnsts, avg_load, kappa, gam1, C):
    """ No generator can offer more reserve of any type than its capacity allows."""
    numConstrs = sum(1 for g in gen_dict.values() if g.res_type == "GEN" and g.fuel_type <> "FixedImport")
    t1s, t2s, t3s, t4s= __addVecVars(model, numConstrs)
    ix = 0
    for name, gen in gen_dict.items():
        if gen.res_type <> "GEN" or gen.fuel_type == "FixedImport":
            continue

        #Notice that the TMSR cap is already handled by the on-off
        if gen.T10_Cap is not None:
            for iHr in xrange(HORIZON_LENGTH):
                f, g = reserve_vars[name, iHr, "TMSR_Cap"]
                f2, g2 = reserve_vars[name, iHr, "TMNSR_Cap"]
                res_cnsts += __addLessEqual( model, numpy.array(f) + numpy.array(f2), 
                                                                    g + g2, gen.T10_Cap, avg_load, kappa, gam1, C, 
                                                                    "T10_Cap" + name + "H" + str(iHr), 
                                                                    t1s[ix], t2s[ix], t3s[ix], t4s[ix])
        if gen.T30_Cap is not None:
            for iHr in xrange(HORIZON_LENGTH):
                f, g = reserve_vars[name, iHr, "TMSR_Cap"]
                f2, g2 = reserve_vars[name, iHr, "TMNSR_Cap"]
                f3, g3 = reserve_vars[name, iHr, "TMOR_Cap"]
                res_cnsts += __addLessEqual(model, numpy.array(f) + numpy.array(f2) + numpy.array(f3), 
                                                                    g + g2 + g3, gen.T30_Cap ,avg_load, kappa, gam1, C, 
                                                                    "T30_Cap" + name + str(iHr), 
                                                                    t1s[ix + 1], t2s[ix + 1], t3s[ix + 1], t4s[ix + 1])
    return res_cnsts

def reserveCapacity(model, gen_dict, reserve_vars, res_cnsts):
    """ No generator can offer more reserve of any type than its capacity allows."""
    for name, gen in gen_dict.items():
        if gen.res_type <> "GEN" or gen.fuel_type == "FixedImport":
            continue
        #Notice that the TMSR cap is already handled by the on-off
        if gen.T10_Cap is not None:
            for iHr in xrange(HORIZON_LENGTH):
                model.addConstr( reserve_vars[name, iHr, "TMSR_Cap"] + 
                                                    reserve_vars[name, iHr, "TMNSR_Cap"] <= gen.T10_Cap )
        if gen.T30_Cap is not None:
            for iHr in xrange(HORIZON_LENGTH):
                model.addConstr( reserve_vars[name, iHr, "TMSR_Cap"] + 
                                                    reserve_vars[name, iHr, "TMNSR_Cap"] +
                                                    reserve_vars[name, iHr, "TMOR_Cap"] <= gen.T30_Cap )
    return res_cnsts

def solveSecondStage(model, balance_cnsts, gen_dict, prod_vars, new_load_by_hr, 
                                              flex_loads, dec_vars_amt, inc_vars, on_vars, start_vars, stop_vars, 
                                              reserve_cnsts, reserve_vars):
    for cnst in balance_cnsts:
        model.remove( cnst )

    for cnst in reserve_cnsts:
        model.remove( cnst )

    #this seems risky but worth it
    for var in reserve_vars.values():
        model.remove(var)

    #"Fixing" the model seems to destroy our handles
    # try just manually constraining all binaries
    for v in on_vars.values():
        model.addConstr( v.x == v )
    for v in start_vars.values():
        model.addConstr( v.x == v)
    for v in stop_vars.values():
        model.addConstr( v.x == v)

    slack, balance_cnsts2 = nominalLoadBalance(model, gen_dict, prod_vars, new_load_by_hr, 
                                              flex_loads, dec_vars_amt, inc_vars)
    return slack, balance_cnsts2, model

def computeGenByHr(gen_dict, prod_vars, flex_loads, dec_vars_amt, inc_vars):
    """Use the solved model to create a nice grid of generation type by day"""
    prod_by_hour = Counter()
    for name, hr in prod_vars.keys():
        prod_by_hour[ "TOTAL", hr ] += prod_vars[name, hr].x

    for name, hr in prod_vars.keys():
        prod_by_hour[ gen_dict[name].fuel_type, hr ] += prod_vars[name, hr].x
        
    for name, hr in flex_loads:
        prod_by_hour["FLEX", hr] += flex_loads[name, hr].x
        
    for name, hr in inc_vars.keys():
        prod_by_hour["INC", int(hr.lstrip("H")) - 1] += inc_vars[name, hr].x
    
    for name, hr, in dec_vars_amt.keys():
        prod_by_hour["DEC", int(hr.lstrip("H")) -1 ] += dec_vars_amt[name, hr].x    

    return prod_by_hour

#########################
if __name__ == "__main__":
    gen_dict = generator.doEverything(adjustRamp=True)
    load_by_hr = readData.readLoads("ISO-data/pRandomAverage_v2.txt")    

#     m = grb.Model("UCNominal")
#     on_vars, start_vars, stop_vars, cost_var = genStage1Vars( m, gen_dict )
#     prod_vars, reserve_vars = genStage2VarsNom(m, gen_dict, True )
#     variable_cost_vars = addPiecewiseCosts(m, gen_dict, prod_vars )
#     startStopConstraints(m, gen_dict, on_vars, start_vars, stop_vars)
#     ecoMinMaxConsts(m, gen_dict, prod_vars, on_vars, reserve_vars)
#     minUpConstraints(m, gen_dict, on_vars)
#     minDownConstraints(m, gen_dict, on_vars)
#     flex_loads = genFlexibleDemandVarsNom( m, gen_dict )
#     reserve_cnsts = reserveRequirements(m, gen_dict, reserve_vars)
#     reserveCapacity(m, gen_dict, reserve_vars, reserve_cnsts)
#     rampingConsts(m, gen_dict, prod_vars, start_vars, stop_vars)
# #    dec_vars_amt, dec_vars_price = genDecVarsNom(m, gen_dict)
# #    inc_vars = genIncVarsNom( m, gen_dict)
#     dec_vars_amt = dec_vars_price = inc_vars = {}
#     slack, balance_cnsts = nominalLoadBalance(m, gen_dict, prod_vars, load_by_hr, 
#                                                         flex_loads, dec_vars_amt, inc_vars)
#     m.update()
#     m.printStats()
#     m.params.mipgap = 1e-2
#     m.optimize()

#make a much sparser generator list.
#     gens_2 = filter(lambda (n, g): g.res_type == "GEN" and g.fuel_type <> "FixedImport", 
#                                     gen_dict.items() )
#     gens_2 = dict(gens_2[:10])
#     flex_loads = filter(lambda (n, g): g.res_type == "LOAD" and g.isFlexible, gen_dict.items() )
#     gens_2.update(dict(flex_loads) )
#     
#     inc_vars = filter(lambda (n, g): g.res_type == "INC", gen_dict.items() )
#     gens_2.update( dict(inc_vars[1:10] ) )
#     
#     dec_vars = filter( lambda (n, g): g.res_type == "DEC", gen_dict.items() )
#     gens_2.update( dict( dec_vars[1:10] ) )
#     len(gens_2)
#     gen_dict = gens_2
# 
#     TMSR_REQ = 62.25
#     T10_REQ  = 124.5
#     T30_REQ = 188.3
# 
#     load_by_hr = readData.readLoads("ISO-data/pRandomAverage_v2.txt")    
#     dummy, load_by_hr = zip(*load_by_hr.items())
#     load_by_hr = [l / 15. for l in load_by_hr ]
# 
###  Now build an affine model
    eps = .05
    kappa = numpy.sqrt( 1/eps  -1 )
    gam1 = .0001
    C = numpy.eye(24) * 1.001
    load_by_hr = readData.readLoads("ISO-data/pRandomAverage_v2.txt")    
    dummy, load_by_hr = zip(*load_by_hr.items())
    Udict = {"kappa":kappa, "gam1":gam1, "C":C, "avg_load": load_by_hr}

    print "Building Affine Model"
    m = grb.Model("UCAffine")
    on_vars, start_vars, stop_vars, cost_var = genStage1Vars( m, gen_dict )
    prod_vars, reserve_vars = genStage2VarsAffine(m, gen_dict)
  
    print "Reserves"
    res_cnsts = reserveRequirementsAffine(m, reserve_vars, Udict)
    reserveCapacityAffine(m, gen_dict, reserve_vars, res_cnsts, load_by_hr, kappa, gam1, C)

    print "Ramping"
    rampingConstsAffine(m , gen_dict, prod_vars, start_vars, stop_vars, 
           load_by_hr, kappa, gam1, C)

    print "Flex Loads"
    flex_loads = genFlexibleDemandVarsAffine( m, gen_dict )
    boundFlexDemandAffine( m, gen_dict, flex_loads, load_by_hr, kappa, gam1, C)

    print "Rest"
    variable_cost_vars = addPiecewiseCosts(m, gen_dict, prod_vars, useAffine=True, Udict=Udict )
    ecoMinMaxConstsAffine(m, gen_dict, prod_vars, on_vars, reserve_vars, load_by_hr, kappa, gam1, C)

    f_sys, g_sys  = formSysLevelAffine( prod_vars, flex_loads )
    startStopConstraints(m, gen_dict, on_vars, start_vars, stop_vars)
    minUpConstraints(m, gen_dict, on_vars)
    minDownConstraints(m, gen_dict, on_vars)

    slack, balance_cnsts = affineLoadBalanceNaive(m, gen_dict, f_sys, g_sys, 
                                                    load_by_hr, kappa, gam1, C)

    m.update()
    m.printStats()
    m.params.mipgap = 1e-2
    m.optimize()
    
    sys.exit()    

    print "Load by Hour:\t"
    for k, v in load_by_hr.items():
        print k, v

    ########
    ### Fix the models and reset the lodads
    for k, v in load_by_hr.items():
        load_by_hr[k] = v

    sys.exit()

#     print "Second Stage Problem"
# 
#     slack, balance_cnsts, m = solveSecondStage(m, balance_cnsts, gen_dict, prod_vars, load_by_hr, 
#                                               flex_loads, dec_vars_amt, inc_vars, on_vars, start_vars, stop_vars, 
#                                               reserve_cnsts, reserve_vars)
# 
#     m.optimize()

    #compute how much is produced each hour.
    prod_by_hour = Counter()
    for name, hr in prod_vars.keys():
        prod_by_hour[ int(hr) ] += prod_vars[name, hr].x

    flex_loads_by_hr = Counter()
    for name, hr in flex_loads.keys():
        flex_loads_by_hr[hr] += flex_loads[name, hr].x

    print "Costs:\t", m.objVal, cost_var.x, m.objVal - cost_var.x
    print "\t Production\t Load \t Slack:\t"
    for ix in xrange(HORIZON_LENGTH):
        print ix, prod_by_hour[ix], load_by_hr[ix], slack[ix].x

    prod_by_type = computeGenByHr(gen_dict, prod_vars, flex_loads, dec_vars_amt, inc_vars)

    print "\t", 
    print "Nuclear \t Hydro \t Steam \t CT \t Diesel \t FixedImport \t Other \t INC \t  Flex \t DEC \t Load \t Slack"
    for hr in xrange(HORIZON_LENGTH):
        print hr,
        print prod_by_type["Nuclear", hr], 
        print prod_by_type["Hydro", hr], 
        print prod_by_type["Steam", hr], 
        print prod_by_type["CT", hr], 
        print prod_by_type["Diesel", hr], 
        print prod_by_type["FixedImport", hr],
        print prod_by_type["Other", hr],
        print prod_by_type["INC", hr],
        print prod_by_type["FLEX", hr],
        print prod_by_type["DEC", hr],
        print load_by_hr[hr], 
        print slack[hr].x

    print "Initial State:", sum( 1 for g in gen_dict.values() if g.init_online > 0 )
    
    print "Starts:\t"
    starts_by_hr, stops_by_hr = Counter(), Counter()
    for name, hr in start_vars.keys():
        starts_by_hr[ hr ] += start_vars[name, hr].x
    
    for name, hr in stop_vars.keys():
        stops_by_hr[ hr ] += stop_vars[name, hr].x
        
    for hr in xrange(HORIZON_LENGTH):
        print hr, starts_by_hr[ hr], stops_by_hr[ hr ]

    sys.exit()
    
## Determine the 2 largest generators by hour and their associated megawatts
    #First identify true generators
    true_gens = filter( lambda (name, g) : g.res_type == "GEN" and 
                                                                        g.fuel_type <> "FixedImport", 
                                     gen_dict.items() )
    #now identify by hour
    for hr in xrange(HORIZON_LENGTH):
        gens_by_hr = sorted( true_gens, 
                                    key= lambda(name, g): prod_vars[name, hr].x )
        gen1 = gens_by_hr[-1][0]
        gen2 = gens_by_hr[-2][0]
        gen3 = gens_by_hr[-3][0]
        print hr, gen1, prod_vars[gen1, hr].x, 
        print gen2, prod_vars[gen2, hr].x, 
        print gen3, prod_vars[gen3, hr].x