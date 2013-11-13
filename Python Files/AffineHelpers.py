"""" Build UC2 Affine

Helper functions for creating the nominal models

"""
import csv, pdb, numpy, sys, math
import gurobipy as grb
import generator, readData
from scipy import stats

#VG THIS IS BAD
HORIZON_LENGTH = 24
EPSILON_ZERO = 0.
PENALTY = 5
TOL = 1e-5

def bootstrap2(sample, statfunc, delta, nsamples = 1000):
    """
    Arguments:
       sample - input sample of values, assummed to be a nparray with obs in rows
       nsamples - number of samples to generate
       statfunc- statistical function to apply to each generated sample.
 
    Returns a 1-delta confidence interval for the statistic
    """
    n = sample.shape[0]
    X = numpy.zeros(nsamples)
    for i in xrange(nsamples):
        indxs = stats.randint.rvs(0, n-1, size=n)
        X[i] = statfunc(sample[indxs])

    X = numpy.sort(X) 
    return X[ [int( nsamples* .5 * delta), int(nsamples * (1 - .5 * delta) ) ] ]

class CSCutGen:
    """A class for the classical ellipsoid heuristic"""
    def __init__(self, data, eps, delta1=.025, delta2=.025, useBootStrap=True):
        """
        Data is N x d... each row is an observation
        """
        self.N = data.shape[0]
        self.d = data.shape[1]
        #calc sample mean and cov
        self.mean = numpy.average(data, axis=0)
        self.cov = numpy.cov(data, rowvar = False)
        self.kappa = math.sqrt(1/eps -1)
        self.reset_delta_boot(data, delta1, delta2)

    def reset_delta_boot(self, data, delta1, delta2):
        """Defines the thresholds for the test in terms of a bootstrap."""
        def muDist(sample):
            return numpy.linalg.norm(numpy.mean(sample, axis=0) - self.mean)
            
        def sigmaDist(sample):
            sigma_samp = numpy.cov(sample, rowvar = False)
            return numpy.linalg.norm(sigma_samp - self.cov)

        dummy, self.gamma1 = bootstrap2(data, muDist, 2 * delta1)
        dummy, self.gamma2 = bootstrap2(data, sigmaDist, 2 * delta2)

        self.shrunk_cov = self.cov + self.gamma2 * numpy.eye( self.d ) 
        self.chol = numpy.linalg.cholesky( self.shrunk_cov )
        self.chol.transpose()

    def addLessEqual(self, model, lhs_terms, xs, rhs, sname, aux_vars = None):
        """Adds the unwrapped version of
                lhs_terms + xs * u <= rhs
                if aux_vars is not None, will use the vars in there"""
        if aux_vars is None:
            t1 = model.addVar()
            t2 = model.addVar()
            t3 = [ model.addVar( lb = -grb.GRB.INFINITY) for ix in xrange(HORIZON_LENGTH) ] 
            t4 = [ model.addVar(lb = -grb.GRB.INFINITY) for ix in xrange(HORIZON_LENGTH) ]
            model.update()
        else:
            (t1, t2, t3, t4) = aux_vars

        objs = [t1, t2]
        objs += t3
        objs += t4
    
        #assume that x represents an expression, add extra variables to unpack it
        #VG: If you want, test on <class 'gurobipy.LinExpr'> or <class 'gurobipy.Var'>
        for x, t in zip(xs, t4):
            objs.append( model.addConstr( x == t ) )

        for ix, (row, t) in enumerate(zip(self.chol, t3)):
            objs.append( model.addConstr( grb.quicksum( r * f for r,f in zip(row[:(ix + 1)], xs) ) == t ) )
    
        objs.append( model.addQConstr( grb.quicksum( t * t for t in t4 ) <= t1*t1, name=sname + "Q1" ) )
        objs.append( model.addQConstr( grb.quicksum( t * t for t in t3 ) <= t2 * t2, name=sname + "Q2" ) )
        objs.append( model.addConstr(lhs_terms +  grb.quicksum( mu * x for mu, x in zip(self.mean, xs) ) + 
                self.gamma1 * t1 + self.kappa * t2 <= rhs, name=sname )  )
    
        return objs
    
    def addGreaterEqual(self, model, lhs_terms, xs, rhs, sname, aux_vars = None):
        """Adds the unwrapped version of
            lhs_terms + xs *u >= rhs """
        if aux_vars is None:
            t1 = model.addVar()
            t2 = model.addVar()
            t3 = [ model.addVar( lb = -grb.GRB.INFINITY) for ix in xrange(HORIZON_LENGTH) ] 
            t4 = [ model.addVar(lb = -grb.GRB.INFINITY) for ix in xrange(HORIZON_LENGTH) ]
            model.update()
        else:
            (t1, t2, t3, t4) = aux_vars

        objs = [ t1, t2]
        objs += t3
        objs += t4
    
        #assume that x represents an expression, add extra variables to unpack it
        #VG test type of x Can test if you want to
        for x, t in zip(xs, t4):
            objs.append( model.addConstr( x == t ) )

        for ix, (row, t) in enumerate(zip(self.chol, t3)):
            objs.append( model.addConstr( grb.quicksum( r * f for r,f in zip(row[:(ix+1)], xs) ) == t ) )
    
        objs.append( model.addQConstr( grb.quicksum( t * t for t in t4 ) <= t1*t1, name=sname + "Q1" ) )
        objs.append( model.addQConstr( grb.quicksum( t * t for t in t3 ) <= t2 * t2, name=sname + "Q2" ) )
        objs.append( model.addConstr(lhs_terms + grb.quicksum( mu * x for mu, x in zip(self.mean, xs) )
                -  self.gamma1 * t1 - self.kappa * t2 >= rhs, name=sname )  )
        return objs
    
    def addVecVars(self, model, numConstrs):
        """Adds all the auxiliary variables needed in one go, and then returns them as a list"""
        out = []
        for ix in xrange(numConstrs):
            t1 = model.addVar()
            t2 = model.addVar()
            t3 = [model.addVar(lb = -grb.GRB.INFINITY) for ix in xrange(HORIZON_LENGTH)] 
            t4 = [model.addVar(lb = -grb.GRB.INFINITY) for ix in xrange(HORIZON_LENGTH)] 
            out.append( (t1, t2, t3, t4) )
        model.update()    
        return out
##################################
class BertSimCutGen:
    """A class for an outter approximation to the CS Cut Generation"""
    def __init__(self, data, eps, delta1=.025, delta2=.025, useBootStrap=True):
        """
        Data is N x d... each row is an observation
        """
        self.N, self.d = data.shape
        #calc sample mean and cov
        self.mean = numpy.average(data, axis=0)
        self.cov = numpy.cov(data, rowvar = False)
        self.kappa = math.sqrt(1/eps -1)
        self.reset_delta_boot(data, delta1, delta2)

    def reset_delta_boot(self, data, delta1, delta2):
        """Defines the thresholds for the test in terms of a bootstrap."""
        def muDist(sample):
            return numpy.linalg.norm(numpy.mean(sample, axis=0) - self.mean)
            
        def sigmaDist(sample):
            sigma_samp = numpy.cov(sample, rowvar = False)
            return numpy.linalg.norm(sigma_samp - self.cov)

        dummy, self.gamma1 = bootstrap2(data, muDist, 2 * delta1)
        dummy, self.gamma2 = bootstrap2(data, sigmaDist, 2 * delta2)

        print "Gammas:\t", self.gamma1, self.gamma2

        self.shrunk_cov = self.cov + self.gamma2 * numpy.eye( self.d ) 
        self.chol = numpy.linalg.cholesky( self.shrunk_cov )
        self.chol.transpose()
        print "MaxChol\t", self.chol.max()

    def addQConstr(self, model, xs, rhs, theta, w1, w2, s1, s2):
        """Adds the simple unwrapped version of the bertsimas sim approx to norm(xs, 2) <= rhs """
        model.addConstr( theta * math.sqrt( len(xs) ) + 
                    grb.quicksum( w1_elem + w2_elem for w1_elem, w2_elem in zip(w1, w2) ) <= rhs)
        for ix in xrange(len(xs) ):
            model.addConstr( w1[ix] - w2[ix] + s1[ix] - s2[ix] == xs[ix] )
            model.addConstr( -s1[ix] -s2[ix] + theta == 0 )

    def addLessEqual(self, model, lhs_terms, xs, rhs, sname, aux_vars = None):
        """Adds the unwrapped version of
                lhs_terms + xs * u <= rhs
                if aux_vars is not None, will use the vars in there"""
        #VG Note that in this new structure we don't return all new objects... this coudl ause errors in updates
        if aux_vars is None:
            aux_vars = self.addVecVars(model, 1)
        t1, t2 = self.createNormVars(model, xs, aux_vars)
        return model.addConstr(lhs_terms +  grb.quicksum( mu * x for mu, x in zip(self.mean, xs) ) + 
                self.gamma1 * t1 + self.kappa * t2 <= rhs, name=sname ) 

    def addGreaterEqual(self, model, lhs_terms, xs, rhs, sname, aux_vars = None):
        """Adds the unwrapped version of
                lhs_terms + xs * u <= rhs
                if aux_vars is not None, will use the vars in there"""
        if aux_vars is None:
            aux_vars = self.addVecVars(model, 1)
        t1, t2 = self.createNormVars(model, xs, aux_vars)
        return model.addConstr(lhs_terms +  grb.quicksum( mu * x for mu, x in zip(self.mean, xs) ) - 
                self.gamma1 * t1 - self.kappa * t2 >= rhs, name=sname ) 
    
    def addLessEqualFast(self, model, lhs_terms, xs, rhs, sname, norm_x, std_x):
        """Adds the unwrapped version of
            lhs_terms + xs *u >= rhs 
            norm_x and std_x are variables corrsponding to norm(x, 2) and sqrt(x'Sigma x)
        """
        return model.addConstr(lhs_terms + grb.quicksum( mu * x for mu, x in zip(self.mean, xs) ) +
                self.gamma1 * norm_x + self.kappa * std_x <= rhs, name=sname ) 

    def addGreaterEqualFast(self, model, lhs_terms, xs, rhs, sname, norm_x, std_x):
        """Adds the unwrapped version of
            lhs_terms + xs *u >= rhs 
            norm_x and std_x are variables corrsponding to norm(x, 2) and sqrt(x'Sigma x)
        """
        return model.addConstr(lhs_terms + grb.quicksum( mu * x for mu, x in zip(self.mean, xs) )
                -  self.gamma1 * norm_x - self.kappa * std_x >= rhs, name=sname )
    
    def addVecVars(self, model,  numVars):
        """Adds all the auxiliary variables needed to define the appropriate upperbounds to norms
        in one go, and then returns them as a list"""
        out = []
        for ix in xrange(numVars):
            t1 = model.addVar()
            t2 = model.addVar()
            t3 = [model.addVar(lb = -grb.GRB.INFINITY) for ix in xrange(HORIZON_LENGTH)] 
            t4 = [model.addVar(lb = -grb.GRB.INFINITY) for ix in xrange(HORIZON_LENGTH)] 
            theta = [model.addVar() for jx in xrange(2)]
            w1 = [model.addVar() for ix in xrange(2 * HORIZON_LENGTH) ] 
            w2 = [model.addVar() for ix in xrange(2 * HORIZON_LENGTH)]
            s1 = [model.addVar() for ix in xrange(2 * HORIZON_LENGTH)]
            s2 = [model.addVar() for ix in xrange(2 * HORIZON_LENGTH)]
            out.append( (t1, t2, t3, t4, theta, w1, w2, s1, s2) )
        model.update()    
        return out

    def createNormVars(self, model, xs, aux_vars=None):
        """Adds only the constraints to define norm(xs, 2) and sqrt(x' Sigmabar x)
        returns variables corresponding to these two"""
        if aux_vars is None:
            print "VG Why is this happening?"
            (t1, t2, t3, t4, theta, w1, w2, s1, s2) = self.addVecVars(model, 1)
        else:
            (t1, t2, t3, t4, theta, w1, w2, s1, s2) = aux_vars

        #assume that x represents an expression, add extra variables to unpack it
        #VG test type of x Can test if you want to
        for x, t in zip(xs, t4):
            model.addConstr( x == t )

        for ix, (row, t) in enumerate(zip(self.chol, t3)):
            model.addConstr( grb.quicksum( r * f for r,f in zip(row[:(ix+1)], xs) ) == t )
        d = len(xs)
        self.addQConstr(model, t4, t1, theta[0], w1[:d], w2[:d], s1[:d], s2[:d])
#        objs.append( model.addQConstr( grb.quicksum( t * t for t in t4 ) <= t1*t1, name=sname + "Q1" ) )
        self.addQConstr(model, t3, t2, theta[1], w1[d:], w2[d:], s1[d:], s2[d:])
#        objs.append( model.addQConstr( grb.quicksum( t * t for t in t3 ) <= t2 * t2, name=sname + "Q2" ) )
        return t1, t2
                        
################################
##################################
# class SparseAffineCutGen:
#     """A class for adding cuts for the sparse affine cut generator
#     Uses a bertsimas -sim style outer approximation for the quadratic constraints"""
#     def __init__(self, data, eps, k, delta1=.025, delta2=.025, useBootStrap=True):
#         """
#         Data is N x d... each row is an observation
#         k is the number of dimensions we allow the affine solution to span
#         """
#         self.N, self.d = data.shape
#         self.k = k
#         #calc sample mean and cov
#         self.mean = numpy.average(data, axis=0)
#         self.cov = numpy.cov(data, rowvar = False)
#         self.kappa = math.sqrt(1/eps -1)
#         self.reset_delta_boot(data, delta1, delta2)
# 
#     def reset_delta_boot(self, data, delta1, delta2):
#         """Defines the thresholds for the test in terms of a bootstrap."""
#         def muDist(sample):
#             return numpy.linalg.norm(numpy.mean(sample, axis=0) - self.mean)
#             
#         def sigmaDist(sample):
#             sigma_samp = numpy.cov(sample, rowvar = False)
#             return numpy.linalg.norm(sigma_samp - self.cov)
# 
#         dummy, self.gamma1 = bootstrap2(data, muDist, 2 * delta1)
#         dummy, self.gamma2 = bootstrap2(data, sigmaDist, 2 * delta2)
# 
#         print "Gammas:\t", self.gamma1, self.gamma2
#         eig_vals, eig_vecs = numpy.linalg.eigh(self.cov + self.gamma2 * numpy.eye(self.d) )
#         indx = numpy.argsort(eig_vals)
#         self.Lambda = numpy.diag( eig_vals[indx[-self.k:]] )
#         self.Qbar = eig_vecs[:, indx[-self.k:] ]
#         self.muQ = numpy.dot(self.mean.T, self.Qbar)
# 
#     def addVecVars(self, model, numConstrs):
#         """Adds all the auxiliary variables needed in one go, and then returns them as a list"""
#         out = []
#         for ix in xrange(numConstrs):
#             t1 = model.addVar()
#             t2 = model.addVar()
#             t3 = [model.addVar(lb = -grb.GRB.INFINITY) for ix in xrange(HORIZON_LENGTH)] 
#             t4 = [model.addVar(lb = -grb.GRB.INFINITY) for ix in xrange(HORIZON_LENGTH)] 
#             theta = [model.addVar() for jx in xrange(2)]
#             w1 = [model.addVar() for ix in xrange(2 * HORIZON_LENGTH) ] 
#             w2 = [model.addVar() for ix in xrange(2 * HORIZON_LENGTH)]
#             s1 = [model.addVar() for ix in xrange(2 * HORIZON_LENGTH)]
#             s2 = [model.addVar() for ix in xrange(2 * HORIZON_LENGTH)]
#             out.append( (t1, t2, t3, t4, theta, w1, w2, s1, s2) )
#         model.update()    
#         return out
# 
#     def addLessEqual(self, model, lhs_terms, fbar, rhs, sname, aux_vars = None):
#         """Adds the unwrapped version of
#                 lhs_terms + Qbar *fbar * u <= rhs
#                 if aux_vars is not None, will use the vars in there"""
#         if aux_vars is None:
#             t1 = model.addVar()
#             t2 = model.addVar()
#             t3 = [ model.addVar( lb = -grb.GRB.INFINITY) for ix in xrange(HORIZON_LENGTH) ] 
#             t4 = [ model.addVar(lb = -grb.GRB.INFINITY) for ix in xrange(HORIZON_LENGTH) ]
#             theta = [model.addVar() for jx in xrange(2)]
#             w1 = [model.addVar() for ix in xrange(2 * HORIZON_LENGTH) ] 
#             w2 = [model.addVar() for ix in xrange(2 * HORIZON_LENGTH)]
#             s1 = [model.addVar() for ix in xrange(2 * HORIZON_LENGTH)]
#             s2 = [model.addVar() for ix in xrange(2 * HORIZON_LENGTH)]
#             model.update()
#         else:
#             (t1, t2, t3, t4, theta, w1, w2, s1, s2) = aux_vars
# 
#         objs = [t1, t2]
#         objs += theta
#         objs += t3 
#         objs += t4 
#         objs += w1 
#         objs += w2 
#         objs += s1 
#         objs += s2
# 
#         objs.append( model.addConstr(lhs_terms +  grb.quicksum( mu * f for mu, f in zip(self.muQ, fbar) ) + 
#                 self.gamma1 * t1 + self.kappa * t2 <= rhs, name=sname )  )
# 
# 
#         #assume that x represents an expression, add extra variables to unpack it
#         #VG: If you want, test on <class 'gurobipy.LinExpr'> or <class 'gurobipy.Var'>
#         for x, t in zip(xs, t4):
#             objs.append( model.addConstr( x == t ) )
# 
#         for ix, (row, t) in enumerate(zip(self.chol, t3)):
#             objs.append( model.addConstr( grb.quicksum( r * f for r,f in zip(row[:(ix + 1)], xs) ) == t ) )
#     
#         d = len(xs)
#         self.addQConstr(model, t4, t1, theta[0], w1[:d], w2[:d], s1[:d], s2[:d], objs)
#         self.addQConstr(model, t3, t2, theta[1], w1[d:], w2[d:], s1[d:], s2[d:], objs)
#     
#         return objs
#     
#     def addGreaterEqual(self, model, lhs_terms, xs, rhs, sname, aux_vars = None):
#         """Adds the unwrapped version of
#             lhs_terms + xs *u >= rhs """
#         if aux_vars is None:
#             t1 = model.addVar()
#             t2 = model.addVar()
#             t3 = [ model.addVar( lb = -grb.GRB.INFINITY) for ix in xrange(HORIZON_LENGTH) ] 
#             t4 = [ model.addVar(lb = -grb.GRB.INFINITY) for ix in xrange(HORIZON_LENGTH) ]
#             theta = [model.addVar() for jx in xrange(2)]
#             w1 = [model.addVar() for ix in xrange(2 * HORIZON_LENGTH) ] 
#             w2 = [model.addVar() for ix in xrange(2 * HORIZON_LENGTH)]
#             s1 = [model.addVar() for ix in xrange(2 * HORIZON_LENGTH)]
#             s2 = [model.addVar() for ix in xrange(2 * HORIZON_LENGTH)]
#             model.update()
#         else:
#             (t1, t2, t3, t4, theta, w1, w2, s1, s2) = aux_vars
# 
#         objs = [t1, t2]
#         objs += theta
#         objs += t4 
#         objs += w1 
#         objs += w2 
#         objs += s1 
#         objs += s2
#         #assume that x represents an expression, add extra variables to unpack it
#         #VG test type of x Can test if you want to
#         for x, t in zip(xs, t4):
#             objs.append( model.addConstr( x == t ) )
# 
#         for ix, (row, t) in enumerate(zip(self.chol, t3)):
#             objs.append( model.addConstr( grb.quicksum( r * f for r,f in zip(row[:(ix+1)], xs) ) == t ) )
#     
#             objs.append( model.addConstr(lhs_terms + grb.quicksum( mu * x for mu, x in zip(self.mean, xs) )
#                 -  self.gamma1 * t1 - self.kappa * t2 >= rhs, name=sname )  )
# 
#         d = len(xs)
#         self.addQConstr(model, t4, t1, theta[0], w1[:d], w2[:d], s1[:d], s2[:d], objs)
#         self.addQConstr(model, t3, t2, theta[1], w1[d:], w2[d:], s1[d:], s2[d:], objs)
#         return objs
#     
#     def addQConstr(self, model, xs, rhs, theta, w1, w2, s1, s2, objs=[]):
#         """Adds the simple unwrapped version of the bertsimas sim approx to norm(xs, 2) <= rhs """
#         objs.append( model.addConstr( theta * math.sqrt( len(xs) ) + 
#                     grb.quicksum( w1_elem + w2_elem for w1_elem, w2_elem in zip(w1, w2) ) <= rhs) )
#         for ix in xrange(len(xs) ):
#             objs.append( model.addConstr( w1[ix] - w2[ix] + s1[ix] - s2[ix] == xs[ix] ) )
#             objs.append( model.addConstr( -s1[ix] -s2[ix] + theta == 0 ) )


##################################


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

def boundFlexDemandAffine( model, gen_dict, flex_loads, flex_loads_norm, 
        funAddLessEqualFast, funAddGreaterEqualFast):
    """Adds upper and lower bounds and the cost variable for the affine formulation for flex.
    funAddLessEqualFast, funAddGreaterFast should be pointers to the appropriate funcitons
    Notice that affine function of residuals, and hence mean performance does not occur here...
    """
    #requisition everything upfront
    revenue_vars = {}
    for (name, iHr), (fvec, g0) in flex_loads.items():
        revenue_vars[name, iHr] = model.addVar(obj = -gen_dict[name].offerBlocks.values()[0][0].price) #Notice neg obj 
    model.update()

    for ix, (k, v) in enumerate(flex_loads.items() ):
        name, iHr = k
        fvec, g0 = v
        #lower bound
        funAddGreaterEqualFast(model, g0, fvec, gen_dict[name].eco_min["H" + str(iHr + 1)], 
                                        "FlexDemandLB" + name + str(iHr), flex_loads_norm[k][0], flex_loads_norm[k][1] )
#         funAddGreaterEqual(model, g0, fvec, gen_dict[name].eco_min["H" + str(iHr + 1)], 
#                                         "FlexDemandLB" + name + str(iHr), aux_vars[3 * ix])

        #upper bound
        funAddLessEqualFast(model, g0, fvec, gen_dict[name].eco_max["H" + str(iHr + 1)], 
                                            "FlexDemandUB" + name + str(iHr), flex_loads_norm[k][0], flex_loads_norm[k][1] )
#         funAddLessEqual(model, g0, fvec, gen_dict[name].eco_max["H" + str(iHr + 1)], 
#                                             "FlexDemandUB" + name + str(iHr), aux_vars[3 * ix + 1] )

        #Revenue
        funAddGreaterEqualFast(model, g0, fvec, revenue_vars[name, iHr], 
                                                "RevenueFlex" + name + str(iHr), flex_loads_norm[k][0], flex_loads_norm[k][1] )
#         funAddGreaterEqual(model, g0, fvec, revenue_vars[name, iHr], 
#                                                 "RevenueFlex" + name + str(iHr), aux_vars[3 * ix + 2] )

def affineLoadBalanceNaive(model, gen_dict, fvecs, fvecs_norm, gprod_sys, slacks, avg_load, USet):
    """After forming the relevant variables, this adds just the constraints for load balance"""
    balance_cnsts = []
    for hr in xrange(HORIZON_LENGTH):
        t1, t2 = fvecs_norm[hr]
        balance_cnsts.append( USet.addLessEqualFast(model, gprod_sys[hr] - avg_load[hr], fvecs[hr], slacks[hr], 
                "LoadBalanceGB" + str(hr), t1, t2) )
        balance_cnsts.append( USet.addGreaterEqualFast(model, gprod_sys[hr] - avg_load[hr], fvecs[hr], -slacks[hr], 
                "LoadBalanceLB" + str(hr), t1, t2) )

    return slacks, balance_cnsts

def formSysLevelAffine( prod_vars, flex_vars ):
    f_sys, g_sys = {}, {}
    for name, hr in prod_vars.keys():
        if hr not in f_sys:
            f_sys[hr] = numpy.zeros(HORIZON_LENGTH)
            g_sys[hr] = 0.
        f, g = prod_vars[name, hr]
        f_sys[hr] =  f_sys[hr] + f
        g_sys[hr] += g

    for name, hr in flex_vars.keys():
        f, g = flex_vars[ name, hr]
        f_sys[hr] = f_sys[hr] - f
        g_sys[hr] -= g

    return f_sys, g_sys

def prepForLoadBalance(model, prod_vars, flex_vars, USet):
    """Forms the variables needed to define load balance L1"""    
    fprod_sys, gprod_sys = formSysLevelAffine( prod_vars, flex_vars )
    #add a slack variable for the amount missed.    
    slack = [model.addVar(obj=PENALTY)  for ix in xrange(HORIZON_LENGTH)]
    aux_fs = [model.addVar(lb = -grb.GRB.INFINITY ) for ix in xrange(HORIZON_LENGTH)]
    aux_vars = USet.addVecVars(model, HORIZON_LENGTH)
    fvecs, fvecs_norm = {}, {}
    for hr in xrange(HORIZON_LENGTH):
        model.addConstr( aux_fs[hr] == fprod_sys[hr] - 1 )
        fvec = list(fprod_sys[hr])
        fvec[hr] = aux_fs[hr]
        fvecs[hr] = fvec
        fvecs_norm[hr] = USet.createNormVars(model, fvec, aux_vars[hr] )

    return slack, fvecs, fvecs_norm, gprod_sys
# def affineLoadBalanceNaive(model, gen_dict, prod_vars, flex_vars, avg_load, 
#                                                     USet):
#     """ Adds a constraint to minimize the L1 deviation from the load."""
#     #Form the system level policy
#     fprod_sys, gprod_sys = formSysLevelAffine( prod_vars, flex_vars )
# 
#     #add a slack variable for the amount missed.    
#     slack = [model.addVar(obj=PENALTY)  for ix in xrange(HORIZON_LENGTH)]
#     aux_fs = [model.addVar(lb = -grb.GRB.INFINITY ) for ix in xrange(HORIZON_LENGTH)]
#     aux_vars = USet.addVecVars(model, HORIZON_LENGTH)
# 
#     balance_cnsts = []        
#     for hr in xrange(HORIZON_LENGTH):
#         balance_cnsts.append( model.addConstr( aux_fs[hr] == fprod_sys[hr] - 1 ) )
#         fvec = [f for f in fprod_sys[hr] ]
#         fvec[hr] = aux_fs[hr]
#         t1, t2 = USet.createNormVars(model, fvec, aux_vars[hr] )
#         USet.addLessEqualFast(model, gprod_sys[hr] - avg_load[hr], fvec, slack[hr], 
#                 "LoadBalanceGB" + str(hr), t1, t2)
# #         balance_cnsts += funAddLessEqual(model, gprod_sys[hr] - avg_load[hr], fvec, slack[hr], 
# #                                 "LoadBalanceGB" + str(hr), aux_vars[2 * hr] )
# 
#         USet.addGreaterEqualFast(model, gprod_sys[hr] - avg_load[hr], fvec, -slack[hr], 
#                 "LoadBalanceLB" + str(hr), t1, t2)
# #         balance_cnsts += funAddGreaterEqual(model, gprod_sys[hr] - avg_load[hr], fvec, -slack[hr], 
# #                                 "LoadBalanceLB" + str(hr), aux_vars[2 * hr + 1] )
# 
#     return slack, balance_cnsts

def rampingConstsAffine(model, gen_dict, prod_vars, start_vars, stop_vars, 
                                                funAddVars, funAddLessEqual, funAddGreaterEqual, sparse=False):
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
#                     __addLessEqual(model, numpy.array(fm) - numpy.array(f), gm - g, 
#                                                 gen.ramp_rate + eco_max_m * stop_vars[name, hr], 
#                                                 avg_load, kappa, gam1, C, 
#                                                 "rampingUB" + name + str(hr) )
                    funAddLessEqual(model, gm -g, numpy.array(fm) - numpy.array(f), 
                                                    gen.ramp_rate + eco_max_m * stop_vars[name, hr], 
                                                    "rampingUB" + name + str(hr) )
                                
                #check the ramping up constraint could be tight
                if gen.ramp_rate < eco_max - eco_min_m:        
#                     __addLessEqual(model, numpy.array(f) - numpy.array(fm), g - gm, 
#                                                 gen.ramp_rate + eco_max * start_vars[name, hr], 
#                                                 avg_load, kappa, gam1, C, 
#                                                 "rampingLB" + name + str(hr))
                    funAddLessEqual(model, g - gm, numpy.array(f) - numpy.array(fm), 
                                                gen.ramp_rate + eco_max * start_vars[name, hr], 
                                                "rampingLB" + name + str(hr) )
    return

def ecoMinMaxConstsAffine(model, gen_dict, prod_vars, on_vars, reserve_vars, USet, M = 1e4):
    """
    on_var * eco_min <= prod_vars + reserve <= on_vars * eco_max
    actually implemented slightly smarter than this to get a locally ideal formulation
    """
    #iterate through setting the eco_max zero elements, and computing how many constraints we need
    numVars = 0
    for name, iHr in on_vars.keys():
        #special case for eco_max == 0
        sHr = "H" + str(iHr + 1)
        if gen_dict[name].eco_max[sHr] < TOL:
            model.addConstr(on_vars[name, iHr] == 0 )
        else:
            numVars +=2 # one for the eco-min/max and one for the spinning reserve           

    #requisition all the variables in one go
    aux_vars =  USet.addVecVars(model, numVars) 
    ix = 0
    for name, iHr in on_vars.keys():
        sHr = "H" + str(iHr + 1)
        #if you're off, no production and no spinning reserve
        fprod, gprod = prod_vars[name, iHr]
        model.addConstr( gprod <= on_vars[name, iHr] * M )
        model.addConstr( gprod >= -on_vars[name, iHr] * M)
        for f in fprod:
            model.addConstr( f <= on_vars[name, iHr] * M )
            model.addConstr( f >= -on_vars[name, iHr] * M )
        fspin, gspin = reserve_vars[name, iHr,  "TMSR_Cap"]
        model.addConstr( gspin <= on_vars[name, iHr] * M )
        model.addConstr( gspin >= -on_vars[name, iHr] * M )
        for f in fspin:
            model.addConstr ( f <= on_vars[name, iHr] * M)
            model.addConstr ( f >= -on_vars[name, iHr] * M )

        #only bother adding the robust constraints if eco-Max is non-zero
        if gen_dict[name].eco_max[sHr] >= TOL:        
            fres = numpy.zeros( HORIZON_LENGTH )
            gres = 0.
            for type in generator.GenUnit.RESERVE_PRODUCTS:
                fres = fres + reserve_vars[name, iHr, type][0]
                gres += reserve_vars[name, iHr, type][1]

            #VG Check this.  Don't we need to account for the start time?
            norm_f, std_f = USet.createNormVars(model, fres  + fprod, aux_vars[ix] )
            USet.addGreaterEqualFast( model, gres + gprod, fres + fprod, 
                    on_vars[name, iHr] * gen_dict[name].eco_min[sHr], 
                    "Ecomin" + name + str(iHr), norm_f, std_f )
#             funAddGreaterEqual(model, gprod + gres, numpy.array(fprod) + fres, 
#                     on_vars[name, iHr] * gen_dict[name].eco_min[sHr], 
#                     "Ecomin" + name + str(iHr),  aux_vars[ix] )
            USet.addLessEqualFast( model, gres + gprod, fres + fprod, 
                    on_vars[name, iHr] * gen_dict[name].eco_max[sHr], 
                    "Ecomax" + name + str(iHr), norm_f, std_f )
            ix += 1
    
            #spinning reserve bounded by TMSR_Cap
            if (name, iHr, "TMSR_Cap") in reserve_vars:
                USet.addLessEqual(model, gspin, fspin, 
                        on_vars[name, iHr] * gen_dict[name].TMSR_Cap, 
                        "OnForSpin" + name + str(iHr), aux_vars[ix] )
#                 funAddLessEqual(model, reserve_vars[name, iHr, "TMSR_Cap"][1], 
#                                                 reserve_vars[name, iHr, "TMSR_Cap"][0], 
#                                                 on_vars[name, iHr] * gen_dict[name].TMSR_Cap, 
#                                                 "OnForSpin" + name + str(iHr), aux_vars[ix + 2] )
                ix += 1     
            else:
                raise ValueError()

def reserveRequirementsAffine(model, reserve_vars, TMSR_REQ, T10_REQ, T30_REQ,
                                                        funAddVars, funAddLessEqual, funAddGreaterEqual):
    cnsts = []
    aux_vars = funAddVars(model, 3 * HORIZON_LENGTH)
    for iHr in xrange(HORIZON_LENGTH):
        TMSR_vars_by_hr = ( v for ((name, hr, type), v) in reserve_vars.items() if 
                                                        hr == iHr and type == "TMSR_Cap" )
        TMSR_fs, TMSR_gs = zip( * TMSR_vars_by_hr )
        f = reduce( lambda x, y: numpy.add(x, y), TMSR_fs )
        g = reduce( lambda x, y: x + y, TMSR_gs )
        cnsts.append ( funAddGreaterEqual(model, g, f, TMSR_REQ, "TMSRReq" + str(iHr), aux_vars[3 * iHr] ) )

        T10_vars_by_hr = ( v for ((name, hr, type), v) in reserve_vars.items() if 
                                                        hr == iHr and type in ("TMSR_CAP", "TMNSR_Cap") )
        T10_fs, T10_gs = zip(* T10_vars_by_hr )            
        f = reduce( lambda x, y: numpy.add(numpy.array(x), numpy.array(y) ), T10_fs )
        g = reduce( lambda x, y: x + y, T10_gs )
        cnsts.append( funAddGreaterEqual(model, g, f, T10_REQ, "T10Req" + str(iHr), aux_vars[3 * iHr + 1] ) )

        T30_vars_by_hr = ( v for ((name, hr, type), v) in reserve_vars.items() if 
                                                        hr == iHr and type in ("TMSR_CAP", "TMNSR_Cap", "TMOR_CAP") )
        T30_fs, T30_gs = zip(* T30_vars_by_hr )            
        f = reduce( lambda x, y: numpy.add(numpy.array(x), numpy.array(y) ), T30_fs )
        g = reduce( lambda x, y: x + y, T30_gs )
        cnsts.append ( funAddGreaterEqual(model, g, f, T30_REQ, "T30Req" + str(iHr), aux_vars[3 * iHr + 2] ) )
    return cnsts

def reserveCapacityAffine(model, gen_dict, reserve_vars, res_cnsts, 
                                                funAddVars, funAddLessEqual, funAddGreaterEqual):
    """ No generator can offer more reserve of any type than its capacity allows."""
    numConstrs = sum(1 for g in gen_dict.values() if g.res_type == "GEN" and g.fuel_type <> "FixedImport")
    aux_vars = funAddVars(model, 2 * numConstrs)
    ix = 0
    for name, gen in gen_dict.items():
        if gen.res_type <> "GEN" or gen.fuel_type == "FixedImport":
            continue

        #Notice that the TMSR cap is already handled by the on-off
        if gen.T10_Cap is not None:
            for iHr in xrange(HORIZON_LENGTH):
                f, g = reserve_vars[name, iHr, "TMSR_Cap"]
                f2, g2 = reserve_vars[name, iHr, "TMNSR_Cap"]
#                 res_cnsts += __addLessEqual( model, numpy.array(f) + numpy.array(f2), 
#                                                                     g + g2, gen.T10_Cap, avg_load, kappa, gam1, C, 
#                                                                     "T10_Cap" + name + "H" + str(iHr), 
#                                                                     t1s[ix], t2s[ix], t3s[ix], t4s[ix])
                res_cnsts.append ( funAddLessEqual(model, g + g2, numpy.array(f) + numpy.array(f2), gen.T10_Cap, 
                                                                    "T10_Cap" + name + "H" + str(iHr), aux_vars[ix] ) )

        if gen.T30_Cap is not None:
            for iHr in xrange(HORIZON_LENGTH):
                f, g = reserve_vars[name, iHr, "TMSR_Cap"]
                f2, g2 = reserve_vars[name, iHr, "TMNSR_Cap"]
                f3, g3 = reserve_vars[name, iHr, "TMOR_Cap"]
#                 res_cnsts += __addLessEqual(model, numpy.array(f) + numpy.array(f2) + numpy.array(f3), 
#                                                                     g + g2 + g3, gen.T30_Cap ,avg_load, kappa, gam1, C, 
#                                                                     "T30_Cap" + name + str(iHr), 
#                                                                     t1s[ix + 1], t2s[ix + 1], t3s[ix + 1], t4s[ix + 1])
                res_cnsts.append( funAddLessEqual(model, g + g2 + g3, numpy.array(f) + numpy.array(f2) + numpy.array(f3), 
                                                                            gen.T30_Cap, "T30_Cap" + name + str(iHr), aux_vars[ix + 1] ) )
        ix +=2

    return res_cnsts

def addPiecewiseCostsAffine(model, gen_dict, prod_vars, funAddVars, funAddLessEqual, funAddGreaterEqual):
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
        sizes = [b.size for b in offerBlocks[:numKnots] ] + [ eco_max ] 
        prices = [b.price for b in offerBlocks[:numKnots] ] + [ offerBlocks[numKnots].price ] 
        size_diffs = [s  - sm for s, sm in zip(sizes, [0] + sizes) ] 
        numKnots += 1

        #simpler formulation relies on the fact that size_diffs are all > 0
        for iHr in xrange(HORIZON_LENGTH):
            ys_all[name, iHr] = [model.addVar(ub=1.0) for ix in xrange(numKnots ) ] 
            cost_vars[name, iHr] = model.addVar(obj=1.0)

    #single call to save time at memory
    model.update()

    #requisition the variables for the affineness
    true_gens = filter( lambda g: g.res_type == "GEN" and gen.fuel_type <> "FixedImport", 
                                     gen_dict.values() )
    aux_vars = funAddVars(model, len(true_gens) * HORIZON_LENGTH)
    for ixGen, gen in enumerate(true_gens):
        if gen.res_type <> "GEN" or gen.fuel_type == "FixedImport":
            continue
        name = gen.name
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
            f, g = prod_vars[name, iHr]
            funAddLessEqual(model, g, f, grb.quicksum( y * s for y, s in zip(ys_all[name, iHr], size_diffs) ), 
                                            "PieceWiseCost" + name + str(iHr), aux_vars[HORIZON_LENGTH * ixGen + iHr] )
    return cost_vars

def summarizeAffine(model, on_vars, start_vars, stop_vars, prod_vars, reserve_vars, variable_cost_vars, flex_loads, 
                                        fixed_cost_var, slacks):
    start_vals = {}
    for name, iHr in start_vars:
        start_vals[name, iHr] = start_vars[name, iHr].x
    on_vals = {}
    for name, iHr in start_vars:
        on_vals[name, iHr] = on_vars[name, iHr].x

    variable_costs = Counter()
    for ((name, hr), v) in variable_cost_vars.items():
        variable_costs[hr] += v.x

    prod_by_hour = {}    
    for hour in xrange(len(slacks)):
        prod_by_hour["Slack", hour] += slacks[hour].x

    return on_vals, start_vals, UCObj.fixed_cost_var.x, UCObj.model.objVal, prod_by_hour, variable_costs

