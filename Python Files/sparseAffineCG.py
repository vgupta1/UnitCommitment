""" 
    An uncertainty set for affine sets where policy specified in eigen basis.
"""
import math, numpy, pdb
from scipy import stats
import gurobipy as grb

#get a better version of bootstrapping for speed.  use the one in scikit?
def bootstrap(sample, statfunc, delta, nsamples = 1000):
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

class SparseAffineCutGen:
    """A class for an outter approximation to the CS Cut Generation"""
    def __init__(self, data, eps, k, model, delta1=.1, delta2=.1):
        """
        Data is N x d... each row is an observation
        k is the number of eigenvectors considered
        """
        self.N, self.d = data.shape
        self.k, self.kappa = k, math.sqrt(1/eps -1)
        self.auxVars = []
        self.model = model

        #part of the logic for approximating 2 norm depends on 
        #particular choice of approximation. This will change based on approx
        #scheme
        self.numAuxVars2Norm = self.k + self.d

        mean = numpy.average(data, axis=0)
        cov = numpy.cov(data, rowvar = False)
        def muDist(sample):
            return numpy.linalg.norm(numpy.mean(sample, axis=0) - mean)
        def sigmaDist(sample):
            sigma_samp = numpy.cov(sample, rowvar = False)
            return numpy.linalg.norm(sigma_samp - cov)

        dummy, self.gamma1 = bootstrap(data, muDist, 2 * delta1)
        dummy, self.gamma2 = bootstrap(data, sigmaDist, 2 * delta2)

        #Cheap Sanity check that zero is in the U set
        print "Gammas:\t", self.gamma1, self.gamma2
        assert numpy.dot(numpy.dot(mean, numpy.linalg.inv(cov)), mean) <= self.kappa**2
        cov += self.gamma2 * numpy.eye( self.d ) 
        eig_vals, eig_vecs = numpy.linalg.eigh(cov)
        indx = numpy.argsort(eig_vals)
        self.lambdas = eig_vals[indx[-self.k:] ]
        self.M = numpy.dot( eig_vecs[:, indx[-self.k:] ], numpy.diag( numpy.sqrt( self.lambdas ) ) )

        self.mu = mean
        self.invChol = numpy.linalg.cholesky( numpy.linalg.inv(cov) ).T
        
    def createVars(self, numConstrs):
        """create and cache the vars for a numConstrs distinct robust constraints"""
        for ix in xrange(numConstrs):
            t1 = self.model.addVar()
            t2 = self.model.addVar()
            w1 = [self.model.addVar(lb = -grb.GRB.INFINITY) for ix in xrange(self.d)] 
            w2 = [self.model.addVar(lb = -grb.GRB.INFINITY) for ix in xrange(self.k)] 
            norm_vars = [self.model.addVar(lb=-grb.GRB.INFINITY) for ix in xrange(self.numAuxVars2Norm) ]
            self.auxVars.append( (t1, t2, w1, w2, norm_vars) )
        self.model.update()    

    #this should be implemented on children
    def addNormConstr(self, w, t, v, tag):
        """Adds an appropriate approximation to || w||_2 \leq t """
        model = self.model
        assert len(w) == len(v)
        for ix, (wi, vi) in enumerate(zip(w, v)):
            model.addConstr(wi <= t, name =  "%sInfNorm%d" % (tag, ix ))
            model.addConstr( wi <= vi, name= "%sOneNormA%d" % (tag, ix) )
            model.addConstr( -wi <= vi, name= "%sOneNormA%d" % (tag, ix) )
        model.addConstr( grb.quicksum(v) <= math.sqrt(self.d) * t )

    def addLessEqual(self, ybar, gbar, b, tag, ixD=None):
        """Adds constraints for ybar ^T M_k d + gbar <= d_{ix} + b 
        If ix is None, drop the d term
        assert len(ybar) == k"""
        if not self.auxVars : #we didn't precache
            print "No Caching?"
            self.createVars(1)

        t1, t2 ,w1, w2, norm_vars = self.auxVars.pop()
        v1 = norm_vars[:len(w1)]
        v2 = norm_vars[len(w1):]
        self.addNormConstr( w1, t1, v1, tag )
        self.addNormConstr( w2, t2, v2, tag )                
        model = self.model

        #VG Check: We are going to remove w1, w2, t1, t2, v1, v2, but not defining constraints
        #check to make sure gurobi knows what it's doing....
        out_objs = [t1, t2, w2, w2] + norm_vars
                
        #VG since no longer using 2nd order constraint, can likely eliminate the variables w1, w2
        for i, (mi, wi) in enumerate(zip(self.M, w1)):
            if ixD <> i:
                out_objs.append( model.addConstr( grb.quicksum( mi[j] * ybar[j] for j in xrange(self.k) ) == wi ) )
            else:
                out_objs.append( model.addConstr( grb.quicksum( mi[j] * ybar[j] for j in xrange(self.k) ) == wi + 1 ) )
            
        if ixD is None:
            for i in xrange(self.k):
                out_objs.append( model.addConstr( ybar[i] * self.lambdas[i] == w2[i] ) )
            #the master equation
            out_objs.append( model.addConstr( numpy.dot(self.mu, numpy.dot(self.M, ybar) ) 
                                                + self.gamma1 * t1 + self.kappa * t2 <= b - gbar, 
                                                name="%sRob" % tag ) )
        else:
            for i in xrange(self.k):
                out_objs.append( model.addConstr(  ybar[i] * self.lambdas[i] - self.M[ixD, i] == w2[i] ) )
            #master equation
            out_objs.append( model.addConstr( numpy.dot(self.mu, numpy.dot(self.M, ybar ) ) - self.mu[ixD] 
                                                + self.gamma1 * t1 + self.kappa * t2 <= b - gbar,  
                                                name = "%sRob" % tag ) )
        return out_objs

    def addGreaterEqual(self, ybar, gbar, b, tag, ixD=None):
        """Adds constraints for ybar ^T M_k d + gbar >= d_{ix} + b 
                If ix is None, drop the d term"""
        if not self.auxVars : #we didn't precache
            print "No Caching?"
            self.createVars(1)

        t1, t2 ,w1, w2, norm_vars = self.auxVars.pop()
        v1 = norm_vars[:len(w1)]
        v2 = norm_vars[len(w1):]
        self.addNormConstr( w1, t1, v1, tag )
        self.addNormConstr( w2, t2, v2, tag )                
        model = self.model
        #VG Check: We are going to remove w1, w2, t1, t2, v1, v2, but not defining constraints
        #check to make sure gurobi knows what it's doing....
        out_objs = [t1, t2, w2, w2] + norm_vars

        for i, (mi, wi) in enumerate(zip(self.M, w1)):
            if ixD <> i:
                temp = model.addConstr( grb.quicksum( mi[j] * ybar[j] for j in xrange(self.k) ) == wi )
            else:
                temp = model.addConstr( grb.quicksum( mi[j] * ybar[j] for j in xrange(self.k) ) == wi + 1 )
            out_objs.append(temp) 
            
        if ixD is None:
            for i in xrange(self.k):
                out_objs.append( model.addConstr( ybar[i] * self.lambdas[i] == w2[i] ) )
            #the master equation
            temp = model.addConstr( numpy.dot(self.mu, numpy.dot(self.M, ybar) ) 
                                                - self.gamma1 * t1 - self.kappa * t2 >= b - gbar, 
                                                name="%sRob" % tag )
            out_objs.append( temp )
        else:
            for i in xrange(self.k):
                out_objs.append( model.addConstr(  ybar[i] * self.lambdas[i] - self.M[ixD, i] == w2[i] ) )
            #master equation
            temp = model.addConstr( numpy.dot(self.mu, numpy.dot(self.M, ybar ) ) - self.mu[ixD] 
                                                - self.gamma1 * t1 - self.kappa * t2 >= b - gbar, 
                                                name = "%sRob" % tag )
            out_objs.append( temp )
        return out_objs
        
    def addBoth(self, ybar, gbar, bg, bl, tag, ixD=None):
        """ d_{ix} + bg <= ybar ^T M_k d + gbar<= d_{ix} + bl 
                If ix is None, drop the d term"""
        if not self.auxVars : #we didn't pre-cache
            print "No Caching?"
            self.createVars(1)

        #check if this logic can be moved to within the variable creation
        #without significantly affecting solve times.
        t1, t2 ,w1, w2, norm_vars = self.auxVars.pop()
        v1 = norm_vars[:len(w1)]
        v2 = norm_vars[len(w1):]
        self.addNormConstr( w1, t1, v1, tag )
        self.addNormConstr( w2, t2, v2, tag )                
        model = self.model
        #VG Check: We are going to remove w1, w2, t1, t2, v1, v2, but not defining constraints
        #check to make sure gurobi knows what it's doing....
        out_objs = [t1, t2, w2, w2] + norm_vars
        
        for i, (mi, wi) in enumerate(zip(self.M, w1)):
            if ixD <> i:
                temp = model.addConstr( grb.quicksum( mi[j] * ybar[j] for j in xrange(self.k) ) == wi )
            else:
                temp = model.addConstr( grb.quicksum( mi[j] * ybar[j] for j in xrange(self.k) ) == wi + 1 )
            out_objs.append(temp)

        if ixD is None:
            for i in xrange(self.k):
                out_objs.append( model.addConstr( ybar[i] * self.lambdas[i] == w2[i] ) )
            #the master equations
            temp = model.addConstr( numpy.dot(self.mu, numpy.dot(self.M, ybar) ) 
                                                - self.gamma1 * t1 - self.kappa * t2 >= bg - gbar, 
                                                name="%sRobGreater" % tag )
            out_objs.append( temp )
            temp = model.addConstr( numpy.dot(self.mu, numpy.dot(self.M, ybar) ) 
                                                + self.gamma1 * t1 + self.kappa * t2 <= bl - gbar, 
                                                name="%sRobLess" % tag )
            out_objs.append( temp )
        else:
            for i in xrange(self.k):
                out_objs.append( model.addConstr(  ybar[i] * self.lambdas[i] - self.M[ixD, i] == w2[i] ) )
            #master equation
            temp = model.addConstr( numpy.dot(self.mu, numpy.dot(self.M, ybar ) ) - self.mu[ixD] 
                                                - self.gamma1 * t1 - self.kappa * t2 >= bg - gbar, 
                                                name = "%sRobGreater" % tag )
            out_objs.append( temp )
            temp = model.addConstr( numpy.dot(self.mu, numpy.dot(self.M, ybar ) ) - self.mu[ixD] 
                                                + self.gamma1 * t1 + self.kappa * t2 <= bl - gbar, 
                                                name = "%sRobLessr" % tag )
            out_objs.append( temp )
        return out_objs

    def sampleConstarint(self, ybar, gbar, b, tag, numPts, isLessEqual, ixD=None):
        """Samples numPts versions of the constraint 
                ybar ^T M_k d + gbar <= d_{ix} + b """
        #assume the sample mean is exact for the sampling
        #simply use a gaussian approximation around the sample mean with shrunken covar
        zetas = numpy.random.randn( self.d, numPts )
        zetas /= map(numpy.linalg.norm, zetas.T)
        zetas *= self.kappa
        us = self.mu + self.invChol * zetas
        Mus = numpy.dot(self.M.T, self.u)
        
        #VG Check the dimensions of Mus and us and stuff
        pdb.set_trace()
        if ixD is None:
            for ix in xrange(numPts):
                self.model.addConstr( 
                        grb.quicksum( ybar_i * mu_i for ybar_i, mu_i in zip(ybar, Mus[:, ix] ) )
                        + gbar <= b, name = tag + "Sample%d" % ix )
            return
        else:
            for ix in xrange(numPts):
                self.model.addConstr( 
                        grb.quicksum( ybar_i * mu_i for ybar_i, mu_i in zip(ybar, Mus[:, ix] ) )
                        + gbar <= b + us[ixD, ix], name = tag +"Sample%d" % ix  )
            return
            
        

