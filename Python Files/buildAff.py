""" 
    Functions to build the Affine model
"""
import csv, pdb, numpy, sys, math
import gurobipy as grb
from config import *
import sparseAffineCG as cg
from buildNom import genStage1Vars, startStopConstraints, minUpConstraints, minDownConstraints
from collections import Counter
import generator

class AffModel:
    def __init__(self, model, onVars, startVars, stopVars, prodVars, reserveVars, variableCostVars, 
                                    flexLoads, fixedCostVar, affCG, genDict):
        self.model, self.onVars, self.startVars, self.stopVars, self.prodVars, self.reserveVars = \
                model, onVars, startVars, stopVars, prodVars, reserveVars 
        self.variableCostVars, self.flexLoads, self.fixedCostVar= variableCostVars, flexLoads, fixedCostVar
        self.slacks, self.balanceObjs, self.fixed_var_cnsts = [], [], [] #VG Combine these two to one object
        self.affCG = affCG
        self.genDict = genDict
        self.trueGens = dict(filter(lambda (n, g): g.res_type == "GEN" and g.fuel_type <> "FixedImport", 
                                genDict.items() ) )
        self.hasCallBack = False

    def removeOldObjs(self):
        for o in self.balanceObjs:
            self.model.remove( o )
        for o in self.slacks:
            self.model.remove( o )
        for o in self.fixed_var_cnsts:
            self.model.remove(  o )

    def computeProdByType(self, genDict):
        #populate the production variables by the g amount, i.e. nominal behavior
        prod_by_hour = Counter()
        for name, hr in self.prodVars:
            prod_by_hour[ "TOTAL", hr ] += self.prodVars[name, hr][1].x
        for name, hr in self.prodVars:
            prod_by_hour[ genDict[name].fuel_type, hr ] += self.prodVars[name, hr][1].x
        for name, hr in self.flexLoads:
            prod_by_hour["FLEX", hr] += self.flexLoads[name, hr][1].x
        for hour in xrange(HORIZON_LENGTH):
            prod_by_hour["Slack", hour] += self.slacks[hour].x
        return prod_by_hour

    def summarizeSolution(self, genDict):
        startVals = {}
        for (name, iHr) in self.startVars:
            startVals[name, iHr] = self.startVars[name, iHr].x
        onVals = {}
        for name, iHr in self.onVars:
            onVals[name, iHr] = self.onVars[name, iHr].x
        variable_costs = Counter()
        for ((name, hr), v) in self.variableCostVars.items():
            variable_costs[hr] += v.x
        prod_by_hour = self.computeProdByType(genDict)
        return onVals, startVals, self.fixedCostVar.x, self.model.objVal, prod_by_hour, variable_costs

    #ramping call-back
    def setRampCallback(self):
        self.hasCallBack = True


    def rampingCallback(self, where):
        """Callback to be used when using sparse ramping constraints"""
        model = self.model
        if where == grb.GRB.callback.MIPSOL:
            for name, gen in self.trueGens.items():
                if gen.fuel_type in ("Steam", "CT"):
                    continue
                for hr in xrange(1, HORIZON_LENGTH):
                    f, g = self.prodVars[name, hr]
                    fm, gm = self.prodVars[name, hr - 1]
                    fval = numpy.array( model.cbGetSolution(f) )
                    fmval = numpy.array( model.cbGetSolution(fm) )
                    gval, gmval = model.cbGetSolution([g, gm])

                    #check to see if its violated
                    eco_max_m = gen.eco_max[hr - 1] 
                    if model.cbGetSolution(self.stopVars[name, hr]) <= .9:
                        resid_star, value = self.affCG.suppFcn(fmval - fval)
                        value += gmval - gval
                        if value > gen.ramp_rate + 1e-5: #VG return to this
                            self.affCG.model.cbLazy(grb.quicksum((fmi - fi) * ri for (fmi, fi, ri) in 
                                                                zip(fm, f, resid_star) ) + gm - g <= 
                                    gen.ramp_rate + eco_max_m * self.stopVars[name, hr])
                    eco_max = gen.eco_max[hr]
                    if model.cbGetSolution(self.startVars[name, hr]) <= .9:
                        resid_star, value = self.affCG.suppFcn(fval - fmval)
                        value += gval - gmval
                        if value > gen.ramp_rate + 1e-5: #VG return to this
                            self.affCG.model.cbLazy(grb.quicksum((fi - fmi) * ri for (fi, fmi, ri) in 
                                                                zip(f, fm, resid_star) ) + g - gm <=
                                    gen.ramp_rate + eco_max * self.startVars[name, hr])
    
def __buildAffNoLoad(affCG, genDict, TMSR_REQ, T10_REQ, T30_REQ, sparseRamps):
    m = affCG.model
    trueGens = dict(filter(lambda (n, g): g.res_type == "GEN" and g.fuel_type <> "FixedImport", 
                                genDict.items() ) )
    onVars, startVars, stopVars, costVar = genStage1Vars(m, trueGens)
    startStopConstraints(m, trueGens, onVars, startVars, stopVars)
    minUpConstraints(m, trueGens, onVars)
    minDownConstraints(m, trueGens, onVars)

    prodVars, reserveVars = genStage2VarsAff(m, trueGens, affCG.k)
#    flexLoads = genFlexibleDemandVarsAff(m, genDict, affCG.k)
    flexLoads = {}
    variableCostVars = addPiecewiseCostsAff(m, trueGens, prodVars, affCG)
    ecoMinMaxConstsAff(m, genDict, prodVars, onVars, reserveVars, affCG)
    reserveRequirementsAff(m, reserveVars, TMSR_REQ, T10_REQ, T30_REQ, affCG)    
    reserveCapacityAff(m, trueGens, reserveVars, affCG )

    o = AffModel(m, onVars, startVars, stopVars, prodVars, reserveVars, variableCostVars, 
                                    flexLoads, costVar, affCG, genDict)
    rampingConstsAff(o, sparseRamps, trueGens)
    return o

def __addLoadBalanceCnst(affModel, predLoads):
    """ Adds a constraint to minimize the L1 deviation from the load."""
    #Assume that old objects have already been removed
    model = affModel.model

    #add a slack variable for the amount missed.    
    slacks = [model.addVar(obj=PENALTY, name="SlackH%d" % ix)  for ix in xrange(HORIZON_LENGTH) ]
    model.update()
    f_sys, g_sys = {}, {}     #remember, we've ignored Inc/Decs
    for name, hr in affModel.prodVars:
        f, g = affModel.prodVars[name, hr]
        if hr not in f_sys:
            f_sys[hr], g_sys[hr] =  numpy.array(f), g
        else:
            f_sys[hr] +=  f 
            g_sys[hr] += g

    for name, hr in affModel.flexLoads:
        f, g = affModel.flexLoads[ name, hr]
        f_sys[hr] -= f
        g_sys[hr] -= g
   
    affModel.affCG.createVars( HORIZON_LENGTH, odd=True )
    balanceObjs = []        
    for hr in xrange(HORIZON_LENGTH):
        balanceObjs += affModel.affCG.addBoth(f_sys[hr], g_sys[hr], 
                                                -slacks[hr] + predLoads[hr], slacks[hr] + predLoads[hr], 
                                                "BalanceH%d" % hr, ixD=hr)         
    affModel.slacks = slacks
    affModel.balanceObjs = balanceObjs
    print "Length of New Objs:\t", len(balanceObjs), len(slacks)
    return

def resolve( affModel, predLoads, genDict, onValsHint={}, startValsHint={}):
    affModel.removeOldObjs()
    model = affModel.model
    #add a hint for the second stage if its available
    for key, val in onValsHint.items():
        affModel.onVars[key].start = val
    for key, val in startValsHint.items():
        affModel.startVars[key].start = val
    
    __addLoadBalanceCnst(affModel, predLoads )    
    
    model.update()
    model.printStats()
    model.params.mipgap = 5e-3
    if affModel.hasCallBack:
        def mycallback(model, where):
            return affModel.rampingCallback(where)
        model.optimize( mycallback )
    else:
#         model.write("tuning.lp")
#         sys.exit()
        model.optimize()
    return affModel.summarizeSolution(genDict)

#############
#  The following functions are all helpers that do various tasks
##############
def genStage2VarsAff(model, trueGens, k):
    """Includes reserve varibles by default
    prod[gen,time] = (fvec, gconstant), 
    reserve[gen, time, type] = gconstant
    len(fvec) = k.
    """
    #VG Experiment computationally with value of adding explicit upper bounds to these variables
    prod, reserve = {}, {}
    for name, gen in trueGens.items():
        for iHr in xrange(HORIZON_LENGTH):
            fvec = numpy.array([model.addVar(lb = -grb.GRB.INFINITY) for ix in xrange(k)]) 
            prod[name, iHr] = (fvec, model.addVar())   #Assumes 0 in U
        
            for cap_type in generator.GenUnit.RESERVE_PRODUCTS:
                reserve[name, iHr, cap_type] = model.addVar(name="%sH%d%s" % (name, iHr, cap_type)) #uses lowebound zero    
    model.update()
    return prod, reserve

def genFlexibleDemandVarsAff( model, genDict, k ):
    """Creates Affine variables for flex demands.  Interpret these as revenue earned and extra load
    to satisfy 
    flex_loads[name, iHr] = (fvec, gconstant) 
    This formulation requires we call a separate function to add the lower and upper bounds and
    associated cost
    """
    raise NotImplementedError()
    #right now the costs of flex loads are not being correctly accounted for....
    flexLoads = {}
    for name, gen in genDict.items():
        if gen.res_type <> "LOAD" or not gen.isFlexible:
            continue
        
        #these flex loads have a single block bid structure
        for iHr in range(HORIZON_LENGTH):
            fvec = [model.addVar(lb = -grb.GRB.INFINITY) for ix in xrange(k) ]             
            flexLoads[name, iHr] = (fvec, model.addVar() )  # assumes u in 0
    return flexLoads

def addPiecewiseCostsAff(model, trueGens, prodVars, affCG):
    """Creates variables and the defining constraints for the piecewise cost functions
    for the variable cost """
    ysAll, costVars = {}, {}
    for name, gen in trueGens.items():
        assert gen.offerBlocks.keys() == ['H1*H36'] #nicety in data set
        offerBlocks = gen.offerBlocks['H1*H36'] 
        #Clip curves at ecoMax.
        ecoMax = max( gen.eco_max.values() )
        for numKnots, b in enumerate(offerBlocks):
            if b.size >= ecoMax:
                break
        numKnots += 1

        #Relies on the fact that size_diffs are all > 0
        for iHr in xrange(HORIZON_LENGTH):
            ysAll[name, iHr] = [model.addVar(ub=1.) for ix in xrange(numKnots) ] 
            costVars[name, iHr] = model.addVar(obj=1.)
    #single call to save time at memory
    model.update()

    #requisition the variables for the affineness
    affCG.createVars(len(trueGens) * HORIZON_LENGTH)
    for ixGen, (name, gen) in enumerate(trueGens.items()):
        offerBlocks = gen.offerBlocks['H1*H36'] 
        ecoMax = max( gen.eco_max.values() ) 
        for numKnots, b in enumerate(offerBlocks):
            if b.size >= ecoMax:
                break
        sizes = [b.size for b in offerBlocks[:numKnots] ] + [ ecoMax ] 
        prices = [b.price for b in offerBlocks[:numKnots] ] + [ offerBlocks[numKnots].price ] 
        size_diffs = [s  - sm for s, sm in zip(sizes, [0] + sizes) ] 
        numKnots += 1

        for iHr in xrange(HORIZON_LENGTH):     
            model.addConstr(costVars[name, iHr] == 
                                grb.quicksum( y * p  * s for y, p, s  in zip(ysAll[name, iHr], prices, size_diffs) ) )    
            ybar, gbar = prodVars[name, iHr]
            affCG.addLessEqual(ybar, gbar, 
                            grb.quicksum(y * s for y, s in zip(ysAll[name, iHr], size_diffs)), 
                            "PieceWiseCost%s%d" % (name, iHr), ixD=None)
    return costVars
#VG Test the value of M is sufficient?
def ecoMinMaxConstsAff(model, genDict, prodVars, onVars, reserveVars, affCG, M = 5e1):
    """
    on_var * eco_min <= prodVars + reserve <= onVars * eco_max
    actually implemented slightly smarter than this to get a locally ideal formulation
    """
    #iterate through setting the eco_max zero elements, and computing how many constraints we need
    numConstr = 0
    for name, iHr in onVars:
        #special case for eco_max == 0
        if genDict[name].eco_max[iHr] < TOL:
            model.addConstr(onVars[name, iHr] == 0 )
            model.addConstr( reserveVars[name, iHr,  "TMSR_Cap"] == 0 )
            model.addConstr( reserveVars[name, iHr,  "TMNSR_Cap"] == 0 )
            model.addConstr( reserveVars[name, iHr,  "TMOR_Cap"] == 0 )

            fprod, gprod = prodVars[name, iHr]
            model.addConstr(gprod == 0 )
            for f in fprod:
                model.addConstr(f == 0 )
        else:
            numConstr += 2 # one for the eco-min/max and one for the spinning reserve           

    #requisition all the variables in one go
    affCG.createVars(numConstr) 
    for name, iHr in onVars:
        if genDict[name].eco_max[iHr]  < TOL:
            continue
        #if you're off, no production and no spinning reserve
        eco_max = genDict[name].eco_max[iHr]
        fprod, gprod = prodVars[name, iHr]
        model.addConstr( gprod <= onVars[name, iHr] * eco_max )  #assumes 0 in U
        bounds = affCG.boundElem( eco_max )
        for f, (l,u) in zip(fprod, bounds):
            model.addConstr( f <= onVars[name, iHr] * u )
            model.addConstr( f >= onVars[name, iHr] * l )  # minus sign is subsumed within l
        
        #spinning reserve bounded by TMSR_Cap
        #Other bounds on reserve handled in the reserve capacity function
        gspin = reserveVars[name, iHr,  "TMSR_Cap"]
        model.addConstr( gspin <= onVars[name, iHr] * genDict[name].TMSR_Cap )

        gres = grb.quicksum(reserveVars[name, iHr, type] for type in genDict[name].RESERVE_PRODUCTS)

        #first constrain the production and reserves to be between eco min-max
        affCG.addBoth(fprod, gprod + gres, 
                onVars[name, iHr] * genDict[name].eco_min[iHr], 
                onVars[name, iHr] * genDict[name].eco_max[iHr], 
                "EcoMinMax%s%d" % (name, iHr) )

        #next production itself must be at least eco_min
        #VG make sure this correctly matches what we're doing in the nominal
        affCG.addGreaterEqual(fprod, gprod,
                onVars[name, iHr] * genDict[name].eco_min[iHr], 
                "ProdAtLeastEcoMin%s%d" % (name, iHr) )  
    return

def reserveRequirementsAff(model, reserveVars, TMSR_REQ, T10_REQ, T30_REQ, affCG):
    reserveSlacks = [ (model.addVar(obj=PENALTY, name="SlackTMSR_REQ%d" % ix), 
                       model.addVar(obj=PENALTY, name="SlackT10_REQ%d" % ix), 
                       model.addVar(obj=PENALTY, name="SlackT30_REQ%d" % ix) ) 
                                    for ix in xrange(HORIZON_LENGTH)]
    model.update()
    for iHr in xrange(HORIZON_LENGTH):
        TMSR_vars_by_hr = ( v for ((name, hr, type), v) in reserveVars.items() if 
                                                        hr == iHr and type == "TMSR_Cap" )
        g = reduce( lambda x, y: x + y, TMSR_vars_by_hr )
        affCG.model.addConstr(g + reserveSlacks[iHr][0] >= TMSR_REQ, name="TMSRReqH%d" % iHr)

        T10_vars_by_hr = ( v for ((name, hr, type), v) in reserveVars.items() if 
                                                        hr == iHr and type in ("TMSR_CAP", "TMNSR_Cap") )
        g = reduce( lambda x, y: x + y, T10_vars_by_hr )
        affCG.model.addConstr( g + reserveSlacks[iHr][1] >= T10_REQ, name="T10ReqH%d" % iHr)

        T30_vars_by_hr = ( v for ((name, hr, type), v) in reserveVars.items() if 
                                                        hr == iHr and type in ("TMSR_CAP", "TMNSR_Cap", "TMOR_CAP") )
        g = reduce( lambda x, y: x + y, T30_vars_by_hr )
        affCG.model.addConstr(g + reserveSlacks[iHr][2] >= T30_REQ, name="T30ReqH%d" % iHr)

def reserveCapacityAff(model, trueGens, reserveVars, affCG):
    """ No generator can offer more reserve of any type than its capacity allows."""
    ix = 0
    for name, gen in trueGens.items():
        #TMSR cap is already handled by the on-off
        for iHr in xrange(HORIZON_LENGTH):
            g = reserveVars[name, iHr, "TMSR_Cap"]
            g2 = reserveVars[name, iHr, "TMNSR_Cap"]
            g3 = reserveVars[name, iHr, "TMOR_Cap"]
            #eco_max 0 case is handled by ecoMinMaxConsts
            #VG Check nominal behavior: total reserves + prod bound by ecoMax.  no spin unless on.
            if gen.eco_max[iHr] <= TOL:
                continue   

            T10_Cap = min(gen.T10_Cap, gen.eco_max[iHr])
            affCG.model.addConstr(g + g2 <= T10_Cap, "T10_Cap%sH%d" % (name, iHr))

            T30_Cap = min(gen.T30_Cap, gen.eco_max[iHr])
            affCG.model.addConstr(g + g2 + g3 <= T30_Cap, "T30_Cap%sH%d" % (name, iHr))

def rampingConstsAff(affModelObj, sparse, trueGens):
    affCG = affModelObj.affCG
    prodVars, startVars, stopVars = affModelObj.prodVars, affModelObj.startVars, affModelObj.stopVars
    for name, gen in trueGens.items():
        if sparse and gen.fuel_type not in ("Steam", "CT"): #only add ramping for marginal turbines
            continue

        #ignore ramping constraints for the first slice.
        for hr in xrange(1, HORIZON_LENGTH):
            ##VG This logic is partially duplicated in "fixRampRates of generator.py"
            # Combine these two.
                eco_min_m = gen.eco_min[hr - 1]
                eco_max_m = gen.eco_max[hr - 1] 
                eco_min = gen.eco_min[hr]
                eco_max = gen.eco_max[hr]
                if eco_max <= TOL: #won't run anyway
                    continue

                f, g = prodVars[name, hr]
                fm, gm = prodVars[name, hr - 1]

                #check the ramping down constraint could be tight
                if gen.ramp_rate < eco_max_m - eco_min:
                        affCG.addLessEqual(fm - f, gm - g, 
                                            gen.ramp_rate + eco_max_m * stopVars[name, hr], 
                                            "rampingUB%sH%d" % (name, hr))
                                
                #check the ramping up constraint could be tight
                #VG There is a possibility to combine logic with the above for shared constraint
                if gen.ramp_rate < eco_max - eco_min_m:        
                    affCG.addLessEqual(f - fm, g - gm, 
                                        gen.ramp_rate + eco_max * startVars[name, hr], 
                                        "rampingLB%sH%d" % (name, hr) )
    
    if sparse:
        affCG.model.params.lazyconstraints = True
        affModelObj.setRampCallback()
        
    return
