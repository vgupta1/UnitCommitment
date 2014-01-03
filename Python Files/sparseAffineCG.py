""" 
    An uncertainty set for affine sets where policy specified in eigen basis.
"""
import math, numpy, pdb
from scipy import stats
import gurobipy as grb
from scikits import bootstrap as boot

#get a better version of bootstrapping for speed.  use the one in scikit?
def bootstrap(sample, statfunc, delta, nsamples=5e6):
    return boot.ci(sample, statfunction=statfunc, alpha = delta, n_samples=int(nsamples) )


class SparseAffineCutGen:
    """A class for an outter approximation to the CS Cut Generation"""
    def __init__(self, data, eps, k, model, delta1=.1, delta2=.1, gamma1 = None, gamma2 = None):
        """
        Data is N x d... each row is an observation
        k is the number of eigenvectors considered
        """
        self.N, self.d = data.shape
        self.k, self.kappa = k, math.sqrt(1/eps -1)
        self.auxVars = []
        self.model = model
        mean = numpy.average(data, axis=0)
        cov = numpy.cov(data, rowvar = False)

        if gamma1 is None:
            def muDist(sample):
                return numpy.linalg.norm(numpy.mean(sample, axis=0) - mean)
            dummy, self.gamma1 = bootstrap(data, muDist, 2 * delta1)
        else:
            self.gamma1 = gamma1
            
        if gamma2 is None:
            def sigmaDist(sample):
                sigma_samp = numpy.cov(sample, rowvar = False)
                return numpy.linalg.norm(sigma_samp - cov)
            dummy, self.gamma2 = bootstrap(data, sigmaDist, 2 * delta2)
        else:
            self.gamma2 = gamma2
            
        #Cheap Sanity check that zero is in the U set
        print "Gammas:\t", self.gamma1, self.gamma2
        assert numpy.dot(numpy.dot(mean, numpy.linalg.inv(cov)), mean) <= self.kappa**2

        cov += self.gamma2 * numpy.eye( self.d ) 
        eig_vals, eig_vecs = numpy.linalg.eigh(cov)
        indx = numpy.argsort(eig_vals)
        self.lambdas = eig_vals[indx[-self.k:] ]
        self.Mk = numpy.dot( eig_vecs[:, indx[-self.k:] ], numpy.diag( numpy.sqrt( self.lambdas ) ) )
        self.M = numpy.dot(eig_vecs[:, indx], numpy.diag( numpy.sqrt( eig_vals[indx] ) ) )
        self.mu = mean
        self.invChol = numpy.linalg.cholesky( numpy.linalg.inv(cov) ).T
        
    def createVars(self, numConstrs, odd=False):
        """create and cache the vars for a numConstrs distinct robust constraints"""
        for ix in xrange(numConstrs):
            t1 = self.model.addVar()
            t2 = self.model.addVar()
            if odd:
                norm_vars = [self.model.addVar(lb=-grb.GRB.INFINITY) for ix in xrange( 2 * self.d ) ]
            else:
                norm_vars = [self.model.addVar(lb=-grb.GRB.INFINITY) for ix in xrange(self.d + self. k) ]
            self.auxVars.append( (t1, t2, norm_vars) )
        self.model.update()    

    #this should be implemented on children
    def addNormConstr(self, w, t, v, tag):
        """Adds an appropriate approximation to || w||_2 \leq t """
        model = self.model
        assert len(w) == len(v)
        sqrt_d  = math.sqrt(len(v))
        out_objs = []
        for ix, (wi, vi) in enumerate(zip(w, v)):
            out_objs += [model.addConstr( grb.LinExpr(wi) <= vi, name= "%sOneNormA%d" % (tag, ix) ),
                                    model.addConstr( grb.LinExpr(-wi) <= vi, name= "%sOneNormB%d" % (tag, ix) ),
                                    model.addConstr(vi * sqrt_d <= t, name =  "%sInfNorm%d" % (tag, ix )) ]
        out_objs.append( model.addConstr( grb.quicksum(v) <= sqrt_d * t ) )
        return out_objs
        
    def addLessEqual(self, ybar, gbar, b, tag, ixD=None):
        """Adds constraints for ybar ^T M_k^T d + gbar <= d_{ix} + b 
        If ix is None, drop the d term
        assert len(ybar) == k"""
        if not self.auxVars : #we didn't precache
#            raise NotImplementedError()
#           VG go back through and fix this....
            print "No Caching?"
            self.createVars(1)

        t1, t2, norm_vars = self.auxVars.pop()
        v1 = norm_vars[:self.d]
        v2 = norm_vars[self.d:]
        model = self.model

        #VG Check: We are going to remove w1, w2, t1, t2, v1, v2, but not defining constraints
        #check to make sure gurobi knows what it's doing....
        out_objs = [t1, t2] + norm_vars
                
        #Create expressions of requisite size for w1 and w2
        w1 = numpy.dot(self.Mk, ybar)
        if ixD is None:
            w2 = [ ybar[i] * self.lambdas[i] for i in xrange(self.k) ]
        else:
            w1[ixD] -= 1 
            w2 = [ybar[i] * self.lambdas[i] - self.M[ixD, i] for i in xrange(self.k) ]
            w2 += [ -self.M[ixD, i] for i in xrange(self.k , self.d ) ]

        self.addNormConstr( w1, t1, v1, "%sw1" % tag )
        self.addNormConstr( w2, t2, v2, "%sw2" % tag )                
        #master equation
        out_objs.append( model.addConstr( numpy.dot(self.mu, w1)  #VG Double check this in all3 palces
                                            + self.gamma1 * t1 + self.kappa * t2 <= b - gbar,  
                                            name = "%sRob" % tag ) )
        return out_objs

    def addGreaterEqual(self, ybar, gbar, b, tag, ixD=None):
        """Adds constraints for ybar ^T M_k d + gbar >= d_{ix} + b 
                If ix is None, drop the d term"""
        if not self.auxVars : #we didn't precache
#            raise NotImplementedError()
            #VG fix this
            print "No Caching?"
            self.createVars(1)

        t1, t2, norm_vars = self.auxVars.pop()
        v1 = norm_vars[:self.d]
        v2 = norm_vars[self.d:]
        model = self.model
        #VG Check: We are going to remove w1, w2, t1, t2, v1, v2, but not defining constraints
        #check to make sure gurobi knows what it's doing....
        out_objs = [t1, t2] + norm_vars
        #Create expressions of requisite size for w1 and w2
        w1 = numpy.dot(self.Mk, ybar)
        if ixD is None:
            w2 = [ ybar[i] * self.lambdas[i] for i in xrange(self.k) ]
        else:
            w1[ixD] -= 1 
            w2 = [ybar[i] * self.lambdas[i] - self.M[ixD, i] for i in xrange(self.k) ]
            w2 += [ -self.M[ixD, i] for i in xrange(self.k, self.d ) ]

        self.addNormConstr( w1, t1, v1, "%sw1" % tag )
        self.addNormConstr( w2, t2, v2, "%sw2" % tag )                
        #the master equation
        temp = model.addConstr( numpy.dot(self.mu, w1 ) 
                                            - self.gamma1 * t1 - self.kappa * t2 >= b - gbar, 
                                            name="%sRob" % tag )
        out_objs.append( temp )
        return out_objs
        
    def addBoth(self, ybar, gbar, bg, bl, tag, ixD=None):
        """ d_{ix} + bg <= ybar ^T M_k d + gbar<= d_{ix} + bl 
                If ix is None, drop the d term"""
        if not self.auxVars : #we didn't pre-cache
            raise NotImplementedError()

        #check if this logic can be moved to within the variable creation
        #without significantly affecting solve times.
        t1, t2, norm_vars = self.auxVars.pop()
        v1 = norm_vars[:self.d]
        v2 = norm_vars[self.d:]
        model = self.model
        #VG Check: We are going to remove t1, t2, v1, v2, but not defining constraints
        #check to make sure gurobi knows what it's doing....
        out_objs = [t1, t2] + norm_vars

        #Create expressions of requisite size for w1 and w2
        w1 = numpy.dot(self.Mk, ybar)
        if ixD is None:
            w2 = [ ybar[i] * self.lambdas[i] for i in xrange(self.k) ]
        else:
            w1[ixD] -= 1 
            w2 = [ybar[i] * self.lambdas[i] - self.M[ixD, i] for i in xrange(self.k) ]
            w2 += [ -self.M[ixD, i] for i in xrange(self.k, self.d ) ]

        out_objs += self.addNormConstr( w1, t1, v1, "%sw1" % tag )
        out_objs += self.addNormConstr( w2, t2, v2, "%sw2" % tag )                
        #the master equations
        temp = model.addConstr( numpy.dot(self.mu, w1 ) 
                                            - self.gamma1 * t1 - self.kappa * t2 >= bg - gbar, 
                                            name="%sRobGreater" % tag )
        out_objs.append( temp )
        temp = model.addConstr( numpy.dot(self.mu, w1 ) 
                                            + self.gamma1 * t1 + self.kappa * t2 <= bl - gbar, 
                                            name="%sRobLess" % tag )
        out_objs.append( temp )
        return out_objs

    def sampleConstraint(self, ybar, gbar, b, tag, numPts, isLessEqual, ixD=None):
        """Samples numPts versions of the constraint 
                ybar ^T M_k d + gbar <= d_{ix} + b """
        #assume the sample mean is exact for the sampling
        #simply use a gaussian approximation around the sample mean with shrunken covar
        zetas = numpy.random.randn( self.d, numPts )
        zetas /= map(numpy.linalg.norm, zetas.T)
        zetas *= self.kappa
        us = self.mu + self.invChol * zetas
        Mus = numpy.dot(self.Mk.T, self.u)
        
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
            
    def subGradNorm(self,  x ):
        """Return a subgradient of the scaled approx norm. 
        x should be a d dimensional vector in cannonical basis (not eigen) """
        sqrt_d = math.sqrt( self.d )
        one_norm = numpy.linalg.norm(x, 1)
        inf_norm = numpy.linalg.norm(x, numpy.inf)
        def sgn(t):
            if t > 0:
                return 1.
            elif t < 0:
                return -1.
            else:
                return 0.
        if one_norm > self.d * inf_norm:
            return numpy.array([sgn(xi) for xi in x]) / sqrt_d
        else:
            indx = numpy.argmax(numpy.absolute(x) )
            out = numpy.zeros( len(x) )
            out[indx] = sgn( x[ indx ]) * sqrt_d
            return out

    def suppFcn(self, y ):
        """ Computes max_u in u y^T M^T u.  returns both ustar and value """
        f = numpy.dot(self.Mk, y)
        #compute the subgradient explicitly
        ustar = self.mu + self.gamma1 * self.subGradNorm( f ) + \
                self.kappa * numpy.dot(self.M, self.subGradNorm( numpy.dot(self.M.T , f ))) ##VG Double check
        val = numpy.dot(ustar, f)

        #VG Debug double check:
        sqrt_d = math.sqrt( self.d )
        one_norm = numpy.linalg.norm(f, 1)
        inf_norm = numpy.linalg.norm(f, numpy.inf)
        mf = numpy.dot(self.M.T, f)
        one_normMf = numpy.linalg.norm(mf, 1)
        inf_normMf = numpy.linalg.norm(mf, numpy.inf)
        val2 = numpy.dot(self.mu, f) + self.gamma1 * max(one_norm / sqrt_d, inf_norm * sqrt_d) + \
            self.kappa * max(one_normMf / sqrt_d, inf_normMf * sqrt_d )
        if abs(val - val2) > 1e-7:
            #Something wonky
            pdb.set_trace()
        return numpy.dot(self.Mk.T, ustar), val