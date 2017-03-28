###########################################################################
#                                                                         #
#    physical_validation,                                                 #
#    a python package to test the physical validity of MD results         #
#                                                                         #
#    Written by Michael R. Shirts <michael.shirts@colorado.edu>           #
#               Pascal T. Merz <pascal.merz@colorado.edu>                 #
#                                                                         #
#    Copyright (C) 2012 University of Virginia                            #
#              (C) 2017 University of Colorado Boulder                    #
#                                                                         #
#    This library is free software; you can redistribute it and/or        #
#    modify it under the terms of the GNU Lesser General Public           #
#    License as published by the Free Software Foundation; either         #
#    version 2.1 of the License, or (at your option) any later version.   #
#                                                                         #
#    This library is distributed in the hope that it will be useful,      #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU    #
#    Lesser General Public License for more details.                      #
#                                                                         #
#    You should have received a copy of the GNU Lesser General Public     #
#    License along with this library; if not, write to the                #
#    Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,     #
#    Boston, MA 02110-1301 USA                                            #
#                                                                         #
###########################################################################
r"""
This file largely corresponds to the checkensemble.py code originally
published on https://github.com/shirtsgroup/checkensemble. It now serves
as the low-level functionality of the high-level module
:mod:`physical_validation.ensemble`.
"""
from __future__ import division
from __future__ import print_function
import numpy
import numpy.random
import scipy
import scipy.optimize
import scipy.stats

from . import plot

#==========================
# HELPER FUNCTIONS
#=========================

def check_twodtype(type):  # check if it's a valid type
    if (type=='dbeta-dpressure') or (type=='dbeta-dmu') or (type=='dbeta-ddmu'):
        #print 'Warning: can\'t do 3D fits currently' 
    # throw an exception?
        return False
    else:
        return True

def PrepConversionFactors(eunits='kJ/mol',punits='bar',vunits='nm^3'):

    if (vunits == 'nm^3') and (punits == 'bar'):
    # default conversion is gromacs nm3*bar to kJ/mol
    # 1 nm3.bar = 0.00060221415 m3.bar / mol, 0.01 m3.bar/mol = 1 kJ/mol --> 6.0221415x10^-2 kJ/mol / nm3/bar
        pvconvert = 0.06221415
    elif (vunits == 'kT' and punits == 'kT'):
        pvconvert = 1
    else:    
        print("I don't know the conversion factor for %s volume units and %s pressure units" % (vunits,punits))
    if (eunits == 'kJ/mol') or (eunits == 'kT'):    
        econvert = 1;
    elif (eunits == 'kcal/mol'):
        pvconvert /= 4.184  # if eunits are in kcal/mol, then we need to change pconvert to kcal/mol / nm3/bar
    else:
        print("I don't know those energy units")

    muconvert = -1*econvert    
    return econvert,pvconvert,muconvert    


def prepare_conversion_factors(units):
    # GROMACS standard units are
    #   energy: kJ/mol
    #   volume: nm^3
    #   pressure: bar
    #   => pV-term: bar * nm^3 == 1e-25 kJ == 6.022140857e-2 kJ/mol
    #   => pvconvert = 6.022140857e-2

    econvert = 1
    pvconvert = 6.022140857e-2

    # UnitData stores conversion factors relative to GROMACS units
    #   energy: energy_conversion * kJ/mol
    #   volume: volume_conversion * nm^3
    #   pressure: pressure_conversion * bar
    #   => pV-term: [p]*[V] == pressure_conversion * volume_conversion bar * nm^3
    pvconvert *= units.pressure_conversion * units.volume_conversion
    pvconvert /= units.energy_conversion

    # TODO: check this once muVT is implemented
    muconvert = -econvert

    return econvert, pvconvert, muconvert


def PrepStrings(type,vunits='kT'):
    
    if (type == 'dbeta-constV'):
        vt = 'E'
        plinfit = r'$-(\beta_2-\beta_1)E$'
        pnlfit = r'$\exp(-(\beta_2-\beta_1)E)$'
        varstring = r'$E (kT)$'
        legend_location = 'upper left'

    elif (type == 'dbeta-constP'):    
        vt = 'H'
        plinfit = r'$-(\beta_2-\beta_1)H$'
        pnlfit = r'$\exp(-(\beta_2-\beta_1)H)$'
        varstring = r'$H (kT)$'
        legend_location = 'upper left'

    elif (type == 'dbeta-constmu'):    
        #averages are A = <E>-mu<N>. 
        vt = 'A'
        plinfit = r'$-(\beta_2-\beta_1)A$'
        pnlfit = r'$\exp(-(\beta_2-\beta_1)A)$'
        varstring = r'$A (kT)$'
        legend_location = 'upper left'

    elif (type == 'dbeta-constdmu'):    
        #averages are <E>-dmu<N>. 
        vt = '\langle E \rangle - \Delta \mu \langle N_2\rangle'
        plinfit = r'$-(\beta_2-\beta_1)A$'
        pnlfit = r'$\exp(-(\beta_2-\beta_1)A)$'
        varstring = r'$A (kT)$'
        legend_location = 'upper left'

    elif (type == 'dpressure-constB'):    
        vt = 'V'
        plinfit = r'$-\beta(P_2-P_1)V$'
        pnlfit = r'$\exp(-\beta(P_2-P_1)V)$'
        varstring = r'$V (' + vunits + r')$'
        legend_location = 'upper right'

    elif (type == 'dmu-constB'):    
        vt = 'N'
        plinfit = r'$\beta(\mu_2-\mu_1)N$'
        pnlfit = r'$\exp(\beta(\mu_2-\mu_1)N)$'
        varstring = r'$N (' + 'number' + r')$'
        legend_location = 'upper left'

    elif (type == 'ddmu-constB'):    
        vt = 'N_2'
        plinfit = r'$\beta(\Delta \mu_2-\ Delta \mu_1)N_2$'
        pnlfit = r'$\exp(\beta(\Delta \mu_2-\Delta \mu_1)N)_2$'
        varstring = r'$N (' + 'number' + r')$'
        legend_location = 'upper left'

    elif (type == 'dbeta-dpressure') or (type == 'dbeta-dmu') or (type == 'dbeta-ddmu'):    
        vt = ''
        plinfit = ''
        pnlfit = ''
        varstring = ''
        legend_location = ''
    else:
        print("Type is not defined for plotting!")

    pstring = 'ln(P_2(' + vt + ')/P_1(' + vt + '))'
    
    return vt,pstring,plinfit,pnlfit,varstring,legend_location       
        

def PrepInputs(N_k,conversions,type='dbeta-constV',beta=None,beta_ave=None,P=None,P_ave=None,mu=None,mu_ave=None,U_kn=None,V_kn=None,N_kn=None):

    """
    useg can be either "scale", where uncertainties are scaled, or "subsample" resulting in subsampled data.
    """

    pvconvert = conversions[1]
    # convenience variables 
    N0 = N_k[0]
    N1 = N_k[1]
    maxN = numpy.max(N_k);

    # Currently seven types; fitting parameters are: 
    # NVT
    # 1) free energy, dbeta, no P, N  - constants are beta_ave, variables (vectors) are E 
    # NPT
    # 2) free energy, dpressure - constants are p_ave, variables (vectors) are V  
    # 3) free energy, dbeta, constP  - constants are beta_ave, variables (vectors) are H 
    # 4) free energy, dbeta, dpressure - constants are beta_ave, p_ave, variables (vectors) are E and V
    # muPT
    # 5) free energy, dmu = constants are mu_ave, variables (vectors) are N
    # 6) free energy, dbeta, constmu  - constants are beta_ave, variables (vectors) are A 
    # 7) free energy, dbeta, dmu - constants are beta_ave, mu_ave, variables (vectors) are E and N
    #
    # 8) free energy, ddmu = constants are dmu_ave, variables (vectors) are N_2
    # 9) free energy, dbeta, constdmu = constants are beta_ave, variables (vectors) are E-\delta\mu\N
    # 10) free energy, dbeta, ddmu - constants are beta_ave, dmu_ave, variables (vectors) are E and N_2

    # NVT types
    if (type == 'dbeta-constV'):
        # allocate space 
        v = numpy.zeros([1,2,maxN],float)        # the variables  
        vr = numpy.zeros([1,2,maxN],float)
        const = numpy.zeros(1,float)             # parameter constants
        dp = numpy.zeros(1,float)                # "true" change in constants

        v[0,0,0:N0] = U_kn[0,0:N0]
        v[0,1,0:N1] = U_kn[1,0:N1]
        const[0] = 0.5*(beta[0] + beta[1])
        dp[0] = beta[0] - beta[1]

    # NPT types    
    elif (type == 'dbeta-constP'):
        # allocate space 
        v = numpy.zeros([1,2,maxN],float)        # the variables  
        vr = numpy.zeros([1,2,maxN],float)
        const = numpy.zeros(1,float)             # parameter constants
        dp = numpy.zeros(1,float)                # "true" change in constants

        v[0,0,0:N0] = U_kn[0,0:N0] + conversions[1]*P_ave*V_kn[0,0:N0]  # everything goes into energy units
        v[0,1,0:N1] = U_kn[1,0:N1] + conversions[1]*P_ave*V_kn[1,0:N1]
        const[0] = 0.5*(beta[0] + beta[1])
        dp[0] = beta[0] - beta[1]

    elif (type == 'dpressure-constB'):    
        # allocate space 
        v = numpy.zeros([1,2,maxN],float)
        vr = numpy.zeros([1,2,maxN],float)
        const = numpy.zeros(1,float)
        dp = numpy.zeros(1,float)

        v[0,0,0:N0] = V_kn[0,0:N0]
        v[0,1,0:N1] = V_kn[1,0:N1]
        const[0] = 0.5*pvconvert*beta_ave*(P[0] + P[1])  # units of 1/volume
        dp[0] = pvconvert*beta_ave*(P[0] - P[1])   # units of 1/volume

    elif (type == 'dbeta-dpressure'):    
        # allocate space 
        v = numpy.zeros([2,2,maxN],float)
        vr = numpy.zeros([2,2,maxN],float)
        const = numpy.zeros(2,float)
        dp = numpy.zeros(2,float)
        v[0,0,0:N0] = U_kn[0,0:N0]
        v[0,1,0:N1] = U_kn[1,0:N1]
        v[1,0,0:N0] = V_kn[0,0:N0]
        v[1,1,0:N1] = V_kn[1,0:N1]
        const[0] = 0.5*(beta[0] + beta[1]) # units of 1/E
        const[1] = 0.5*pvconvert*(P[0] + P[1])  # units of E/V?
        dp[0] = beta[0] - beta[1]   # units of 1/Energy
        dp[1] = pvconvert*(beta[0]*P[0] - beta[1]*P[1])  # units of 1/Volume

# mu V T types         

    elif (type == 'dbeta-constmu'):
        # allocate space 
        v = numpy.zeros([1,2,maxN],float)        # the variables  
        vr = numpy.zeros([1,2,maxN],float)
        const = numpy.zeros(1,float)             # parameter constants
        dp = numpy.zeros(1,float)                # "true" change in constants

        v[0,0,0:N0] = U_kn[0,0:N0] - mu_ave*N_kn[0,0:N0]  # everything goes into energy units
        v[0,1,0:N1] = U_kn[1,0:N1] - mu_ave*N_kn[1,0:N1]  # everything is in energy units. 
        const[0] = 0.5*(beta[0] + beta[1])
        dp[0] = beta[0] - beta[1]

    elif (type == 'dmu-constB'):    
        # allocate space 
        v = numpy.zeros([1,2,maxN],float)
        vr = numpy.zeros([1,2,maxN],float)
        const = numpy.zeros(1,float)
        dp = numpy.zeros(1,float)

        v[0,0,0:N0] = N_kn[0,0:N0]
        v[0,1,0:N1] = N_kn[1,0:N1]
        const[0] = -1*(beta_ave*(mu[0] + mu[1]))  # units of 1/energy?
        dp[0] = -1*(beta_ave*(mu[0] - mu[1]))   # units of 1/energy?
        
    elif (type == 'dbeta-dmu'):    
        # allocate space 
        v = numpy.zeros([2,2,maxN],float)
        vr = numpy.zeros([2,2,maxN],float)
        const = numpy.zeros(2,float)
        dp = numpy.zeros(2,float)
        v[0,0,0:N0] = U_kn[0,0:N0]
        v[0,1,0:N1] = U_kn[1,0:N1]
        v[1,0,0:N0] = N_kn[0,0:N0]
        v[1,1,0:N1] = N_kn[1,0:N1]
        const[0] = 0.5*(beta[0] + beta[1]) # units of 1/E
        const[1] = -0.5*(mu[0] + mu[1])  # units of E/V?
        dp[0] = beta[0] - beta[1]   # units of 1/Energy
        dp[1] = -(beta[0]*mu[0] - beta[1]*mu[1])  # units of 1/number?

    else:
        print("Warning:  Type of analysis %s is not defined!" % (type))

    return dp,const,v,vr

def LogLikelihood(x,N_k,const,v):

    L = len(x)

    N0 = N_k[0]
    N1 = N_k[1]
    N  = N0+N1

    M = numpy.log((1.0*N1)/(1.0*N0))

    #D0 = M + beta_ave*x[0] + U0*x[1]
    #D1 = M + beta_ave*x[0] + U1*x[1]

    D0 = D1 = M + const[0]*x[0]
    for i in range(L-1):
        D0 = D0 + v[i,0,0:N0]*x[i+1]
        D1 = D1 + v[i,1,0:N1]*x[i+1]

    # old version:
    # E0 = 1 + numpy.exp(D0)
    # E1 = 1 + numpy.exp(-D1)
    #
    # # this is the negative of the log likelihood, since we want to maximize it using fmin
    #
    # of = ((numpy.sum(numpy.log(E0)) + numpy.sum(numpy.log(E1)))) / N

    def log_1_plus_exp(y):
        def f(yy):
            with numpy.errstate(over='raise'):
                try:
                    xx = numpy.log(1 + numpy.exp(yy))
                except FloatingPointError:
                    xx = yy + numpy.log(1 + numpy.exp(-yy))
                return xx
        return numpy.vectorize(f)(y)

    E0 = log_1_plus_exp(D0)
    E1 = log_1_plus_exp(-D1)

    # this is the negative of the log likelihood, since we want to maximize it using fmin

    of = (numpy.sum(E0) + numpy.sum(E1)) / N

    return of

def dLogLikelihood(x,N_k,const,v):
    """
    Derivative with respect to the parameters, to aid the minimization.
    
    """

    L = len(x)

    N0 = N_k[0]
    N1 = N_k[1]
    N  = N0+N1

    M = numpy.log((1.0*N1)/(1.0*N0))

    D0 = D1 = M + const[0]*x[0]
    for i in range(L-1):
        D0 = D0 + v[i,0,0:N0]*x[i+1]
        D1 = D1 + v[i,1,0:N1]*x[i+1]

    # old version:
    # E0 = 1/(1 + numpy.exp(-D0))
    # E1 = 1/(1 + numpy.exp(D1))

    def inv_1_plus_exp(y):
        def f(yy):
            with numpy.errstate(over='raise'):
                try:
                    xx = 1. / (1 + numpy.exp(yy))
                except FloatingPointError:
                    xx = 0.
                return xx
        return numpy.vectorize(f)(y)

    E0 = inv_1_plus_exp(-D0)
    E1 = inv_1_plus_exp(D1)

    g = numpy.zeros(L,dtype=numpy.float64)

    #this is the gradient of -log likelihood
    #g[0] = (1.0/N)*(numpy.sum(beta*E0) - numpy.sum(beta*E1))
    #g[1] = (1.0/N)*(numpy.sum(U0*E0) - numpy.sum(U1*E1))
    #g[2] = (1.0/N)*(numpy.sum(V0*E0) - numpy.sum(V1*E1))

    g[0] = const[0]*(numpy.sum(E0) - numpy.sum(E1))
    for i in range(L-1):
        g[i+1] = numpy.sum(v[i,0,0:N0]*E0) - numpy.sum(v[i,1,0:N1]*E1)
    return (1.0/N)*g

def d2LogLikelihood(x,N_k,const,v):

    """

    beta = const[0]
    pave = const[1]

    if D = M + beta*x[0] + x[1]*U
    I = \sum_{i=1}^N [[-beta^2/S,-beta*U/S],[-beta*U/S,-U^2/S]] where S = [(1+exp(-D))*(1+exp(D))]

    if D = M + beta*x[0] + x[1]*U + x[2]*V
    I = \sum_{i=1}^N [[-beta^2/S,-beta*U/S,-beta*V/S],[-beta*U/S,-U^2/S,-U*V/S],[-beta*V/S,-U*V^2/S,-V^2/S]] where S = [(1+exp(-D))*(1+exp(D))]
    """

    L = len(x)

    N0 = N_k[0]
    N1 = N_k[1]
    N  = N0+N1
    M = numpy.log((1.0*N1)/(1.0*N0))

    vall = numpy.zeros([L-1,N],dtype=numpy.float64)
    for i in range(L-1):
        vall[i,0:N0] = v[i,0,0:N0]
        vall[i,N0:N] = v[i,1,0:N1]
    
    D = M + const[0]*x[0]
    for i in range(L-1):
        D = D + vall[i,:]*x[i+1]

    hf = numpy.zeros([L,L,N],dtype=numpy.float64)

    cones = const[0] * numpy.ones(N,dtype=numpy.float64)

    # fix this to match the     
    for i in range(L):
        if (i == 0):
            a = cones
        else: 
            a = vall[i-1,:]
        for j in range(L):
            if (j == 0): 
                b = cones
            else:
                b = vall[j-1,:]    
            hf[i,j,:] = a*b

    def log_E(y):
        def f(yy):
            with numpy.errstate(over='raise'):
                try:
                    numpy.exp(yy)
                except FloatingPointError:
                    return yy + numpy.log(1 + 2*numpy.exp(-yy) + numpy.exp(-2*yy))
                try:
                    numpy.exp(-yy)
                except FloatingPointError:
                    return -yy + numpy.log(1 + 2*numpy.exp(yy) + numpy.exp(2*yy))
                return numpy.log(2 + numpy.exp(yy) + numpy.exp(-yy))
        return numpy.vectorize(f)(y)

    #E = (1 + numpy.exp(-D)) * (1 + numpy.exp(D))
    hf_E = hf*numpy.exp(-log_E(D))

    # this is the hessian of the minimum function (not the max)
    h = -numpy.sum(hf_E,axis=2)/N
                            
    return h

def SolveMaxLike(x, N_k, const, v, tol = 1e-10, maxiter=20, quiet=False):
    
    converge = False
    itol = 1e-2
    rtol = 1e-2
    lasttol = 100000;

    for i in range(maxiter):
        lastx = x
        gx = numpy.transpose(dLogLikelihood(x,N_k,const,v))
        nh = d2LogLikelihood(x,N_k,const,v)
        dx = numpy.linalg.solve(nh,gx)
        if numpy.any(numpy.isinf(dx)):
            checktol = numpy.inf
            checkrtol = numpy.inf
        else:
            x += dx  # adding, not subtracting because of the handling of negatives
            rx = dx/x
            checktol = numpy.sqrt(numpy.dot(dx,dx))
            checkrtol = numpy.sqrt(numpy.dot(rx,rx))
        if (checkrtol < tol):
            break
            converge = True
        if (checkrtol > 1.0) and (checktol > lasttol):  # we are possibly diverging. Switch to cg for a bit.
            x = scipy.optimize.fmin_cg(LogLikelihood,lastx,fprime=dLogLikelihood,gtol=itol,args=(N_k,const,v),disp=(not quiet))
            itol *= rtol
            dx = x - lastx
            checktol = numpy.sqrt(numpy.dot(dx,dx))
        lasttol = checktol

    if (i == maxiter) and (converge == False):
        print("Too many iterations, convergence failing")

    return x    

def MaxLikeUncertain(x,N_k,const,v,vave):

    L = len(x)
    d = numpy.zeros(L,float)

    # multiply back by N, since we were dealing with smaller numbers for numerical robustness.
    fi = -(N_k[0] + N_k[1])*d2LogLikelihood(x,N_k,const,v)

    d2 = numpy.linalg.inv(fi)

    # We have a fit to the line y = m(x-Uave) + b, so to get the uncertainty in the free energy back, we need 
    # to add M*Uave back to f.  The uncertainty in cov(b + m*Uave,b+m*Uave) = var(b) + Uave**2*var(m) + Uave*cov(v,m)

    # For two dimensioms, we have the line y = m1(x1-vave1) + m2(x2-vave2) + b
    #  Uncertainty will be cov(b + m1vave1 + m2vave2) = var(b) + vave1^2 var(m1) + vave2^2 var(m2)
    #                                                          + 2vave1 cov(m1,b) + 2vave2 cov(m2,b)
    #                                                          + 2vave1 cov(m1,m2)
    d[0] = const[0]**2*d2[0,0] 
    for i in range(1,L):
        d[0] += vave[i-1]**2*d2[i,i] + 2*vave[i-1]*d2[0,i]   # should this last one be plus or minus
        d[i] = d2[i,i]
        for j in range(i+1,L-1):
            d[0] += 2*vave[i-1]*vave[j-1]*d2[i,j]
    d = numpy.sqrt(d)
    return d

def MaxLikeParams(N_k,dp,const,v,df=0,analytic_uncertainty=False,g=1,quiet=False):
    
    L = len(const)
    optimum = numpy.zeros(L+1,float)
    vave = numpy.zeros(L,dtype=numpy.float64)
    vmod = numpy.zeros([L,2,numpy.max(N_k)],dtype=numpy.float64)
    # for numerical stability, we need to translate the curve
    for i in range(L):
        vave[i] = (numpy.sum(v[i,0,0:N_k[0]]) + numpy.sum(v[i,1,0:N_k[1]]))/numpy.sum(N_k)
        vmod[i,0,0:N_k[0]] = v[i,0,0:N_k[0]] - vave[i]
        vmod[i,1,0:N_k[1]] = v[i,1,0:N_k[1]] - vave[i]

    xstart = numpy.zeros(L+1,float)
    for i in range(L):
        xstart[0] += vave[i]*dp[i]
        xstart[i+1] = dp[i]
    xstart[0] += df
    xstart[0] /= const[0]

    ofit = SolveMaxLike(xstart,N_k,const,vmod,tol=1e-10,quiet=quiet)

    optimum[0] = ofit[0]*const[0]
    for i in range(L):
        optimum[i+1] = ofit[i+1]
        optimum[0] -= (vave[i]*ofit[i+1])

    results = []    
    results.append(optimum)
    if (analytic_uncertainty):
        doptimum = MaxLikeUncertain(ofit,N_k,const,vmod,vave)*numpy.sqrt(numpy.average(g))
        results.append(doptimum)

    return results

#========================================================================================
# Functions for computing Bennett acceptance ratio
#==========================================================================================
def logsum(a_n):
  """
  Compute the log of a sum of exponentiated terms exp(a_n) in a numerically-stable manner:

    logsum a_n = max_arg + \log \sum_{n=1}^N \exp[a_n - max_arg]

  where max_arg = max_n a_n.  This is mathematically (but not numerically) equivalent to

    logsum a_n = \log \sum_{n=1}^N \exp[a_n]

  ARGUMENTS
    a_n (numpy array) - a_n[n] is the nth exponential argument
  
  RETURNS
    log_sum (float) - the log of the sum of exponentiated a_n, log (\sum_n exp(a_n))

  EXAMPLE  
    
  """

  # Compute the maximum argument.
  max_log_term = numpy.max(a_n)

  # Compute the reduced terms.
  terms = numpy.exp(a_n - max_log_term)

  # Compute the log sum.
  log_sum = numpy.log(sum(terms)) + max_log_term
        
  return log_sum

#=============================================================================================
# Bennett acceptance ratio function to be zeroed to solve for BAR.
#=============================================================================================
def BARzero(w_F,w_R,DeltaF):
    """
    ARGUMENTS
      w_F (numpy.array) - w_F[t] is the forward work value from snapshot t.
                        t = 0...(T_F-1)  Length T_F is deduced from vector.
      w_R (numpy.array) - w_R[t] is the reverse work value from snapshot t.
                        t = 0...(T_R-1)  Length T_R is deduced from vector.

      DeltaF (float) - Our current guess

    RETURNS

      fzero - a variable that is zeroed when DeltaF satisfies BAR.
    """

    # Recommended stable implementation of BAR.

    # Determine number of forward and reverse work values provided.
    T_F = float(w_F.size) # number of forward work values
    T_R = float(w_R.size) # number of reverse work values

    # Compute log ratio of forward and reverse counts.
    M = numpy.log(T_F / T_R)
    
    # Compute log numerator.
    # log f(W) = - log [1 + exp((M + W - DeltaF))]
    #          = - log ( exp[+maxarg] [exp[-maxarg] + exp[(M + W - DeltaF) - maxarg]] )
    #          = - maxarg - log[exp[-maxarg] + (T_F/T_R) exp[(M + W - DeltaF) - maxarg]]
    # where maxarg = max( (M + W - DeltaF) )
    exp_arg_F = (M + w_F - DeltaF)
    max_arg_F = numpy.choose(numpy.greater(0.0, exp_arg_F), (0.0, exp_arg_F))
    log_f_F = - max_arg_F - numpy.log( numpy.exp(-max_arg_F) + numpy.exp(exp_arg_F - max_arg_F) )
    log_numer = logsum(log_f_F) - numpy.log(T_F)
    
    # Compute log_denominator.
    # log_denom = log < f(-W) exp[-W] >_R
    # NOTE: log [f(-W) exp(-W)] = log f(-W) - W
    exp_arg_R = (M - w_R - DeltaF)
    max_arg_R = numpy.choose(numpy.greater(0.0, exp_arg_R), (0.0, exp_arg_R))
    log_f_R = - max_arg_R - numpy.log( numpy.exp(-max_arg_R) + numpy.exp(exp_arg_R - max_arg_R) ) - w_R 
    log_denom = logsum(log_f_R) - numpy.log(T_R)

    # This function must be zeroed to find a root
    fzero  = DeltaF - (log_denom - log_numer)

    return fzero

def BAR(w_F, w_R, DeltaF=0.0, maximum_iterations=500, relative_tolerance=1.0e-10, verbose=False):
  """
  Compute free energy difference using the Bennett acceptance ratio (BAR) method using false position

  ARGUMENTS
    w_F (numpy.array) - w_F[t] is the forward work value from snapshot t.
                        t = 0...(T_F-1)  Length T_F is deduced from vector.
    w_R (numpy.array) - w_R[t] is the reverse work value from snapshot t.
                        t = 0...(T_R-1)  Length T_R is deduced from vector.

  OPTIONAL ARGUMENTS

    DeltaF (float) - DeltaF can be set to initialize the free energy difference with a guess (default 0.0)
    maximum_iterations (int) - can be set to limit the maximum number of iterations performed (default 500)
    relative_tolerance (float) - can be set to determine the relative tolerance convergence criteria (defailt 1.0e-5)
    verbose (boolean) - should be set to True if verbse debug output is desired (default False)

  RETURNS

    [DeltaF, dDeltaF] where dDeltaF is the estimated std dev uncertainty

  REFERENCE

    [1] Shirts MR, Bair E, Hooker G, and Pande VS. Equilibrium free energies from nonequilibrium
    measurements using maximum-likelihood methods. PRL 91(14):140601, 2003.

  EXAMPLES

  Compute free energy difference between two specified samples of work values.

  """

  UpperB = numpy.average(w_F)
  LowerB = -numpy.average(w_R)

  FUpperB = BARzero(w_F,w_R,UpperB)
  FLowerB = BARzero(w_F,w_R,LowerB)
  nfunc = 2;
    
  if (numpy.isnan(FUpperB) or numpy.isnan(FLowerB)):
      # this data set is returning NAN -- will likely not work.  Return 0, print a warning:
      print("Warning: BAR is likely to be inaccurate because of poor sampling. Guessing 0.")
      return [0.0, 0.0]
      
  while FUpperB*FLowerB > 0:
      # if they have the same sign, they do not bracket.  Widen the bracket until they have opposite signs.
      # There may be a better way to do this, and the above bracket should rarely fail.
      if verbose:
          print('Initial brackets did not actually bracket, widening them')
      FAve = (UpperB+LowerB)/2
      UpperB = UpperB - max(abs(UpperB-FAve),0.1)
      LowerB = LowerB + max(abs(LowerB-FAve),0.1)
      FUpperB = BARzero(w_F,w_R,UpperB)
      FLowerB = BARzero(w_F,w_R,LowerB)
      nfunc += 2

  # Iterate to convergence or until maximum number of iterations has been exceeded.

  for iteration in range(maximum_iterations):

      DeltaF_old = DeltaF
    
      # Predict the new value
      if (LowerB==0.0) and (UpperB==0.0):
        DeltaF = 0.0
        FNew = 0.0
      else:
        DeltaF = UpperB - FUpperB*(UpperB-LowerB)/(FUpperB-FLowerB)
        FNew = BARzero(w_F,w_R,DeltaF)
      nfunc += 1
     
      if FNew == 0: 
        # Convergence is achieved.
        if verbose: 
          print("Convergence achieved.")
        relative_change = 10^(-15)
        break

      # Check for convergence.
      if (DeltaF == 0.0):
          # The free energy difference appears to be zero -- return.
          if verbose: print("The free energy difference appears to be zero.")
          return [0.0, 0.0]
        
      relative_change = abs((DeltaF - DeltaF_old)/DeltaF)
      if verbose:
          print("relative_change = %12.3f" % relative_change)
          
      if ((iteration > 0) and (relative_change < relative_tolerance)):
          # Convergence is achieved.
          if verbose: 
              print("Convergence achieved.")
          break

      if FUpperB*FNew < 0:
          # these two now bracket the root
          LowerB = DeltaF
          FLowerB = FNew
      elif FLowerB*FNew <= 0:
          # these two now bracket the root
          UpperB = DeltaF
          FUpperB = FNew
      else:
          message = 'WARNING: Cannot determine bound on free energy'
          raise BoundsError(message)        

      if verbose:
        print("iteration %5d : DeltaF = %16.3f" % (iteration, DeltaF))

  # Report convergence, or warn user if not achieved.
  if iteration < maximum_iterations:
      if verbose: 
          print('Converged to tolerance of %e in %d iterations (%d function evaluations)' % (relative_change, iteration,nfunc))
  else:
      message = 'WARNING: Did not converge to within specified tolerance. max_delta = %f, TOLERANCE = %f, MAX_ITS = %d' % (relative_change, tolerance, maximum_iterations)
      raise ConvergenceException(message)

  # Compute asymptotic variance estimate using Eq. 10a of Bennett, 1976 (except with n_1<f>_1^2 in 
  # the second denominator, it is an error in the original
  # NOTE: The numerical stability of this computation may need to be improved.
      
  # Determine number of forward and reverse work values provided.

  T_F = float(w_F.size) # number of forward work values
  T_R = float(w_R.size) # number of reverse work values
  
  # Compute log ratio of forward and reverse counts.
  M = numpy.log(T_F / T_R)

  T_tot = T_F + T_R

  C = M-DeltaF
  
  fF =  1/(1+numpy.exp(w_F + C))
  fR =  1/(1+numpy.exp(w_R - C))
  
  afF2 = (numpy.average(fF))**2
  afR2 = (numpy.average(fR))**2
  
  vfF = numpy.var(fF)/T_F
  vfR = numpy.var(fR)/T_R
  
  variance = vfF/afF2 + vfR/afR2
  
  dDeltaF = numpy.sqrt(variance)
  if verbose: 
      print("DeltaF = %8.3f +- %8.3f" % (DeltaF, dDeltaF))
  return (DeltaF, dDeltaF)

def Print1DStats(title,type,fitvals,convert,trueslope,const,dfitvals='N/A'):

    # if dB, 'convert' is kB
    # if dP, 'convert' is beta*PV_convert
    # first element in fitvals is free energies df
    dfs = fitvals[0]
    # second element in fitvals is the slope
    slopes = fitvals[1]

    # Need to fix this so that uncertainties aren't printed when ddf is 'N/A'

    df = numpy.average(dfs) # true even if there is only one per slope
    if (numpy.size(dfs) > 1):
        ddf  = numpy.std(dfs)
    else:
        ddf = dfitvals[0]

    slope = numpy.average(slopes) # true even if there is only one
    if (numpy.size(slopes) > 1):    
        dslope = numpy.std(slopes)
    else:
        dslope = dfitvals[1]

    print("")
    print("---------------------------------------------")
    print("     %20s        " % (title))
    print("---------------------------------------------")
    print("     df = %.5f +/- %.5f " % (df,ddf))
    print("---------------------------------------------")
    print("     Estimated slope       vs.   True slope")
    print("---------------------------------------------")
    print("%11.6f +/- %11.6f  |  %11.6f" % (slope, dslope, trueslope))
    print("---------------------------------------------")

    quant = numpy.abs((slope-trueslope)/dslope)
    print("")
    print("(That's %.2f quantiles from true slope=%5f, FYI.)" % (quant,trueslope))
    if (quant > 5):
        print(" (Ouch!)")
    else:
        print("")

    if (type[0:5] == 'dbeta'):    
        #trueslope = B1 - B0, const = (B1 + B0)/2, B = 1/(k_B T)
        # so B0 = (const-trueslope/2), T0 = 1/(k_B*B0)
        # so B1 = (const+trueslope/2), T1 = 1/(k_B*B1)
        T0 = (convert*(const-trueslope/2))**(-1)
        T1 = (convert*(const+trueslope/2))**(-1)

        print("---------------------------------------------")
        print(" True dT = %7.3f, Eff. dT = %7.3f+/-%.3f" % (T0-T1, convert*T0*T1*slope,convert*dslope*T0*T1))
        print("---------------------------------------------")

    elif (type == 'dpressure-constB'):
        # trueslope = B*PV_conv*(P1-P0), const = B*PV_conv*(P1+P0)/2, 
        # we need to convert this slope to a pressure.  This should just be dividing by pvconvert*beta
        #
        print("---------------------------------------------")
        print(" True dP = %7.3f, Eff. dP = %7.3f+/-%.3f" % (-trueslope/convert, -slope/convert, numpy.abs(dslope/convert)))
        print("---------------------------------------------")

    elif (type == 'dmu-constB'):
        # trueslope = B*(mu1-mu0), const = B*(mu1+mu0)/2, 
        # we need to convert this slope to a chemical potential.  This should just be dividing by beta
        #
        print("---------------------------------------------")
        print(" True dmu = %7.3f, Eff. dmu = %7.3f+/-%.3f" % (-trueslope/convert, -slope/convert, numpy.abs(dslope/convert)))
        print("---------------------------------------------")

    elif (type == 'ddmu-constB'):
        # trueslope = B*(dmu1-dmu0), const = B*(dmu1+dmu0)/2, 
        # we need to convert this slope to a chemical potential.  This should just be dividing by beta
        #
        print("---------------------------------------------")
        print(" True ddmu = %7.3f, Eff. ddmu = %7.3f+/-%.3f" % (-trueslope/convert, -slope/convert, numpy.abs(dslope/convert)))
        print("---------------------------------------------")

def Print2DStats(title,type,fitvals,kB,convertback,trueslope,const,dfitvals='N/A'):

    # first element in fitvals is free energies df
    dfs = fitvals[0]
    # Need to fix this so that uncertainties aren't printed when ddf is 'N/A'
    df = numpy.average(dfs) # true even if there is only one per slope
    if (numpy.size(dfs) > 1):
        ddf  = numpy.std(dfs)
    else:
        ddf = dfitvals[0]

    slopes = []
    # second element in fitvals is the energy slope
    # third element in fitvals is the PV slope 
    for i in range(2):
        slopes.append(fitvals[i+1])

    slope = numpy.zeros(2,float)
    dslope = numpy.zeros(2,float)
    for i in range(2):
        slope[i] = numpy.average(slopes[i]) # true even if there is only one
        if (numpy.size(slopes[i]) > 1):    
            dslope[i] = numpy.std(slopes[i])
        else:
            dslope[i] = dfitvals[i+1]

    print("")
    print("---------------------------------------------------")
    print("     %20s        " % (title))
    print("---------------------------------------------------")
    print("     df = %.5f +/- %.5f " % (df,ddf))
    for i in range(2):    
        print("---------------------------------------------------")
        print("     Estimated slope[%d]       vs.   True slope[%d]" % (i,i))
        print("---------------------------------------------------")
        print("%11.6f +/- %11.6f  |  %11.6f" % (slope[i], dslope[i], trueslope[i]))
        
        quant = numpy.abs((slope[i]-trueslope[i])/dslope[i])
        print("")
        print("(That's %.2f quantiles from true slope=%5f, FYI.)" % (quant,trueslope[i]))
        if (quant > 5):
            print(" (Ouch!)")
        else:
            print("")

    #dp = B1 - B0, const = (B1 + B0)/2, B = 1/kbT
    # so B0 = (const[0]-trueslope[0]/2), T0 = 1/(kB*B0)
    # so B1 = (const[0]+trueslope[0]/2), T1 = 1/(kB*B1)
    T0 = (kB*(const[0]-trueslope[0]/2))**(-1)
    T1 = (kB*(const[0]+trueslope[0]/2))**(-1)

    print("---------------------------------------------")
    print(" True dT = %7.3f, Eff. dT = %7.3f+/-%.3f" % (T0-T1, kB*T0*T1*slope[0],kB*dslope[0]*T0*T1))
    print("---------------------------------------------")
    
    if (type == 'dbeta-dpressure'):
        text = 'dP'
    elif (type == 'dbeta-dmu'):
        text = 'dmu' 
    elif (type == 'dbeta-ddmu'):
        text = 'ddmu' 
        
    print("---------------------------------------------")
    print(" True %s = %7.3f, Eff. %s = %7.3f+/-%.3f" % (text,-trueslope[1]/convertback,text, -slope[1]/convertback, dslope[1]/convertback))


def PrintPicture(xaxis,true,y,dy,fit,type,name,figname,fittype,vunits='kT',show=False):

    try:
        import matplotlib
    except:
        print('*************')
        print('Note: Figures not generated because matplotlib not found. ', end=' ')
        print('Please install matplotlib to allow generation of pictures')
        return

    import matplotlib.pyplot as plt

    matplotlib.rc('lines',lw=2)
    font = {'family' : 'serif',
            'weight' : 'bold',
            'size'   : '14'}
    matplotlib.rc('font',**font)

    [vt,pstring,plinfit,pnlfit,varstring,legend_location] = PrepStrings(type,vunits=vunits)

    pstringtex = r'$\frac{P_2(' + vt + r')}{P_1(' + vt + r')}$' 
    pstringlntex = r'$\ln\frac{P_2(' + vt + r')}{P_1(' + vt + r')}$' 

    print("Now printing figure %s" % (figname))
    plt.clf()
    plt.xlabel(varstring)
    if (fittype == 'linear'):
        plt.title(vt + ' vs. log probability ratio \n for ' + name)
        plt.errorbar(xaxis,y,fmt='b-',yerr=dy,label = pstringlntex)  # make this general!
        plt.errorbar(xaxis,true,fmt='k-',label = plinfit)
        plt.errorbar(xaxis,fit,fmt='r-',label = 'Fit to $y = b+aB$')
        plt.ylabel(pstringlntex)
    elif (fittype == 'nonlinear'):
        plt.title(vt + ' vs. probability ratio \n for ' + name)
        plt.errorbar(xaxis,y,fmt='b-',yerr=dy,label = pstringtex) 
        plt.errorbar(xaxis,true,fmt='k-',label = pnlfit)
        plt.errorbar(xaxis,fit,fmt='r-',label = 'Fit to $y = \exp(b+aE)$')
        plt.ylabel(pstringtex)
    elif (fittype == 'maxwell'):
        # only valid for kinetic energy
        plt.title('E_kin vs. probability \n for' + name)
        plt.errorbar(xaxis,y,fmt='b-',yerr=dy,label = r'$P(E_{\mathrm{kin}})$')
        if (true is not None):  # sometimes, true will be none.
            plt.errorbar(xaxis,true,fmt='k-',label = 'Fit to Analytical')
        plt.errorbar(xaxis,fit,fmt='r-',label = 'Fit to Normal')
        plt.ylabel(r'$P(E_{\mathrm{kin}})$')
    else:
        print("I'm crying foul!  %s is not an allowed chart type!" % (fittype))

    plt.legend(loc=legend_location)
    if show:
        plt.show()
    plt.savefig(figname + '.pdf')

def GenHistogramProbs(N_k,bins,v,g):

    K = len(N_k)
    
    hlist = []
    dhlist = []

    # check for zero bins
    left = 0
    right = len(bins)
    center = int(len(bins) / 2)
    for k in range(0,K):
        hstat = numpy.histogram(v[0, k, 0:N_k[k]], bins=bins)
        empty = numpy.where(hstat[0] == 0)[0]
        for e in empty:
            if e < center:
                if e >= left:
                    left = e + 1
            else:
                if e < right:
                    right = e
    bins = bins[left:right]
    for k in range(0,K):
        hstat = numpy.histogram(v[0,k,0:N_k[k]], bins = bins)
        h = (1.0*hstat[0])/N_k[k] 
        hlist.append(h)
        dh = numpy.sqrt(g[k]*h*(1.0-h)/N_k[k])
        dhlist.append(dh)

    return hlist,dhlist,bins

def LinFit(bins, N_k, dp, const, v, df=0,
           analytic_uncertainty=False, name="", g=(1, 1),
           screen=False, filename=None, eunits=None):

    [hlist,dhlist,bins] = GenHistogramProbs(N_k,bins,v,g)

    ratio = numpy.log(hlist[1]/hlist[0]) # this should have the proper exponential distribution 
    dratio = numpy.sqrt((dhlist[0]/hlist[0])**2 + (dhlist[1]/hlist[1])**2)

    usedat = numpy.isfinite(ratio)
    y = ratio[usedat]
    nuse = len(y)
    weights = 1.0/dratio[usedat]

    xaxis = (bins[0:len(bins)-1] + bins[1:len(bins)])/2    
    x = xaxis[usedat]

    X = numpy.ones([nuse,2],float)
    X[:,1] = x

    w = numpy.diag(weights) 
    WX = numpy.dot(w,X)
    WY = numpy.dot(w,y)
    WXT = numpy.transpose(WX)
    Z = numpy.dot(WXT,WX)
    WXY = numpy.dot(WXT,WY)

    a = numpy.linalg.solve(Z,WXY)
    da_matrix = numpy.transpose(numpy.linalg.inv(Z))
    da = numpy.zeros(2,float)
    da[0] = numpy.sqrt(da_matrix[0,0])
    da[1] = numpy.sqrt(da_matrix[1,1])

    # the true line is y = df + dp*x, where y is ln P_1(X)/P_2(X)

    do_plot = screen or filename is not None
    if do_plot:
        trueslope = dp
        true = df+trueslope*xaxis 
        fit = a[0] + a[1]*xaxis

        #PrintData(xaxis,true,fit,ratio,dratio,'linear')

        name = name + ' (linear)'
        #PrintPicture(xaxis,true,ratio,dratio,fit,type,name,figname,'linear',vunits)

        data = [{'x': xaxis,
                 'y': ratio,
                 'y_err': dratio,
                 'name': 'Simulation'},
                {'x': xaxis,
                 'y': fit,
                 'name': 'Fit to simulation'},
                {'x': xaxis,
                 'y': true,
                 'name': 'Analytical ratio'}]

        units = ''
        if eunits is not None:
            units = ' [' + eunits + ']'

        annot = None
        if analytic_uncertainty:
            annot = ('{:.1f}'.format(abs((a[1] - dp[0]) / da[1])) +
                     ' quantiles')

        plot.plot(data,
                  legend='best',
                  title='Log probability ratio',
                  xlabel='Energy' + units,
                  ylabel=r'$\log\frac{P_2(E)}{P_1(E)}$',
                  filename=filename,
                  screen=screen,
                  axtext=annot)



    results = []    
    results.append(a)
    if (analytic_uncertainty):
        results.append(da)

    return results

def SolveNonLin(f,df,a,data,ddata,xaxis,maxiter=20,tol=1e-10):

    K = numpy.size(a)
    usedat = numpy.isfinite(data)
    y = data[usedat]
    nuse = len(y)
    weights = 1.0/ddata[usedat]
    w = numpy.diag(weights) 
    x = xaxis[usedat]
    J = numpy.zeros([nuse,K],dtype = numpy.float64)

    # do the newton-raphson solution
    endnext = False
    for n in range(maxiter):
        
        expt = f(a,x)
        
        J = numpy.transpose(df(a,x))
        WJ = numpy.dot(w,J)
        JTW = numpy.transpose(WJ)
        dy = y - expt
        Z = numpy.dot(JTW,WJ)
        incr_a = numpy.linalg.solve(Z,numpy.dot(JTW,dy))
        a += incr_a
        ra = incr_a/a
        chtol = numpy.sqrt(numpy.dot(ra,ra))
        if (chtol < tol):
            if (endnext == True) or (analytical_estimate == False):
                    endnow == True   # we put in this logic so that we calculate the matrix at the minimum
                                     # if we want the analytic uncertainty 
            endnext = True
            if (endnow == True):
                break

        if (n == maxiter):
             print("Too many iterations for nonlinear least squares")

    da_matrix = numpy.linalg.inv(Z)
    da = numpy.zeros(K,float)
    for k in range(K):
        da[k] = numpy.sqrt(da_matrix[k,k])

    return a,da    

def ExpFit(a,x):  # assume only 2D, since we are not generating histograms
    return numpy.exp(a[0] + a[1]*x)

def dExpFit(a,x):
    s = a[0] + a[1]*x
    e = numpy.exp(s)
    return numpy.array([e,x*e])

def NonLinFit(bins,N_k,dp,const,v,df=0,analytic_uncertainty=False,bGraph=False,name="",
              figname='nonlin_figure', tol=1e-10,g=[1,1], type = 'dbeta-constV',vunits='kT'):

    # nonlinear model is exp(A + B*E_i), where the i are the bin energies.
    # residuals are y_i - exp(A + B*E_i)
    # dS/dbeta_j = 2\sum_i r_i dr_i/dB_j = 0 
    # 
    # dr_i/dA = exp(A + B*E_i)
    # dr_i/dB = E_i*exp(A + B*E_i)

    [hlist,dhlist,bins] = GenHistogramProbs(N_k,bins,v,g)

    ratio = (hlist[1]/hlist[0]) # this should have the proper exponential distribution 
    dratio = ratio*(numpy.sqrt((dhlist[0]/hlist[0])**2 + (dhlist[1]/hlist[1])**2))

    xaxis = (bins[0:len(bins)-1] + bins[1:len(bins)])/2    

    # starting point for nonlinear fit
    L = numpy.size(dp)+1
    a = numpy.zeros(L)
    a[0] = df
    a[1:L] = dp[:]

    (a,da) = SolveNonLin(ExpFit,dExpFit,a,ratio,dratio,xaxis,tol=tol)

    if (bGraph):
        trueslope = dp
        true = numpy.exp(df+trueslope*xaxis) 
        fit = ExpFit(a,xaxis)
        
        PrintData(xaxis,true,fit,ratio,dratio,'nonlinear')

        name = name + ' (nonlinear)'
        PrintPicture(xaxis,true,ratio,dratio,fit,type,name,figname,'nonlinear',vunits=vunits)

    results = []    
    results.append(a)
    if (analytic_uncertainty):
        results.append(da)

    return results

def MaxwellBoltzFit(bins,U,N,kT,figname,name="",ndof=None,g=1):

    # generate histogram
    hstat = numpy.histogram(U, bins = bins)
    # normalize the histogram
    h = (1.0*hstat[0])/N 
    # compute the error bars
    dh = numpy.sqrt(g*h*(1.0-h)/N)
    xaxis = (bins[0:len(bins)-1] + bins[1:len(bins)])/2    

    # we assume we have the correct mean for now, since presumably the temperature works
    mean = numpy.mean(U) 
    std_fit = numpy.std(U)
    std_true = numpy.sqrt(mean*kT)
    if (mean > 50*kT):  #if too big, we use a normal distribution -- we'll use limit of 50 DOF as suggest (by Wikipedia!)
        # note that for a normal distribution, the sample mean and standard deviation give the maximum likelihood information.
        fit = numpy.exp(-(xaxis-mean)**2/(2*std_fit**2))/(numpy.sqrt(2*numpy.pi*std_fit**2))
        true = numpy.exp(-(xaxis-mean)**2/(2*std_true**2))/(numpy.sqrt(2*numpy.pi*std_true**2))
        # check this with paper?
    else:
        # should be a gamma distribution; no std fit
        fit = (xaxis/kT)**(mean/kT-1)*numpy.exp(-xaxis/kT)/(kT*scipy.special.gamma(mean/kT))
        if (ndof is not None):
            mean_true = 0.5*ndof*kT
            true = (xaxis/kT)**(mean_true/kT-1)*numpy.exp(-xaxis/kT)/(kT*scipy.special.gamma(mean/kT))
        else:
            true = None # unless we know the number of DOF, we don't know the true distribution:

    print("--- Kinetic energy analysis ---")
    print("")
    print("kT = %10.4f" % (kT))
    if (ndof is None):
        print("Effective # of DOF = %10.4f" % (2*mean/kT))
    else:     
        print("Reported # of DOF = %10.4f" % ndof) 
    if (mean > 25*kT):
        "Approximating the Maxwell-Boltzmann with a normal distribution, as # DOF > 50"
    print("Direct Std = %10.4f, Std from sqrt(U*kT) = %10.4f" % (std_fit,std_true))
    print("")

    # name = name + str
    # normalize histogram and error bars
    width = bins[1]-bins[0]  # assumes equal spacing (currently true)
    h /= width
    dh /= width
    PrintPicture(xaxis,true,h,dh,fit,'dbeta-constV',name,figname,'maxwell')
    
def PrintData(xaxis,true,fit,collected,dcollected,type):

    if (type == 'linear'):
        print("----  Linear Fit  ----")
    elif (type == 'nonlinear'):
        print("----  Nonlinear Fit  ----")
    elif (type == 'maxwell'):
        print("----  fit to Maxwell-Boltzmann ----")
    else:
        quit("Sorry, no allowed type of fit (%s) specified!" % (type))

    print("      X         True     Observed     Error   d(true/obs) sig(true/obs)  Fit   ")
    print("---------------------------------------------------------------------------------------")
    for i in range(len(collected)):
        diff = collected[i]-true[i]
        sig = numpy.abs(collected[i]-true[i])/dcollected[i]
        print("%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f" % (xaxis[i],true[i],collected[i],dcollected[i],diff,sig,fit[i]))


def ProbabilityAnalysis(N_k, type='dbeta-constV',
                        T_k=None, P_k=None, mu_k=None, U_kn=None, V_kn=None, N_kn=None,
                        kB=0.0083144624, nbins=40,
                        bMaxLikelihood=True, bLinearFit=True, bNonLinearFit=True, reptype=None, nboots=200,
                        g=[1,1], reps=None, cuttails=0.001, bMaxwell=False,
                        units=None, seed=None,
                        screen=False, filename=None, quiet=False):

    K = len(N_k)  # should be 2 pretty much always . . . 

    # decide if we are printing figures:
    do_plot = screen or filename is not None

    # get correct conversion terms between different units.
    conversions = 1, 1, 1
    if units is not None:
        conversions = prepare_conversion_factors(units)

    if (seed):
        numpy.random.seed(seed)  # so there is the ability to make the RNG repeatable
        if not quiet:
            print("setting random number seed for bootstrapping %d" % (seed))
    # initialize constant terms
    beta_ave = None
    P_ave = None
    mu_ave = None

    if (T_k is None):
        T_k = numpy.zeros(2,float)
    else:
        beta_k = (1.0/(kB*T_k))
        beta_ave = numpy.average(beta_k)

    if (P_k is None):    
        P_k = numpy.zeros(2,float)    
    else:
        P_ave = numpy.average(P_k)

    if (mu_k is None):    
        mu_k = numpy.zeros(2,float)    
    else:
        mu_ave = numpy.average(mu_k)
 
   # prepare the variables we are going to work with    
    [dp,const,v,vr] = PrepInputs(N_k,conversions,type,beta_k,beta_ave,P_k,P_ave,mu_k,mu_ave,U_kn,V_kn,N_kn)
    [vt,pstring,plinfit,pnlfit,varstring,legend_location] = PrepStrings(type)
    
    if (check_twodtype(type)):  # if it's 2D, we can graph, otherwise, there is too much histogram error 

        # determine the bin widths
        maxk = numpy.zeros(K,float)
        mink = numpy.zeros(K,float)

        # cuttails indicates how many we leave out on each tail
        # for now, we choose the range that cuts 0.1% from the tails of the smallest distribution.
        prange = cuttails

        for k in range(K):
            maxk[k] = scipy.stats.scoreatpercentile(v[0,k,0:N_k[k]],100*(1-prange))
            mink[k] = scipy.stats.scoreatpercentile(v[0,k,0:N_k[k]],100*(prange))

        binmax = numpy.min(maxk)
        binmin = numpy.max(mink)

        if (binmax < binmin):
            quit("\nThere is no overlap between the two ensembles; high distribution min is %f and low distribution max is %f; quitting!" %(binmin,binmax));
            
        if (type == 'dmu-constB'):   # special code for N, since it's discrete
                                     # if the spread is greater than
            if ((binmax-binmin) < nbins): 
                bins = numpy.arange(binmin-0.5,binmax+0.6,1)
                nbins = len(bins)
            else:
                print("Warning: since the range of particle number is greater than the ", end=' ')
                print("number of bins specified, particle number is not discrete for ", end=' ')
                print("methods using histograms.  Set nbins larger (using --nbins) to ")
                print("obtain discrete N distributions")
            bins = numpy.zeros(nbins+1,float)
            for i in range(nbins+1):
                bins[i] = binmin + (binmax-binmin)*(i/(1.0*nbins))    
        else:
            bins = numpy.zeros(nbins+1,float)
            for i in range(nbins+1):
                bins[i] = binmin + (binmax-binmin)*(i/(1.0*nbins))    

    #===================================================================================================
    # Calculate free energies with different methods
    #===================================================================================================    

    if (type == 'dbeta-dpressure') or (type == 'dbeta-dmu'):
        if (dp[0] == 0):
            quit("Warning: two input temperatures are equal, can't do joint variable fit!");
        if (dp[1] == 0):
            if (type == 'dbeta-dpressure'):
                quit("Warning: two input pressures are equal, can't do joint E,V fit!");            
            elif (type == 'dbeta-dmu'):
                quit("Warning: two input chemical potentials are equal, can't do joint E,N fit!");
    else:
        trueslope = dp
        if not quiet:
            print("True slope of %s should be %.8f" % (pstring,trueslope))
  
    if type == 'dpressure-constB' or type == 'dbeta-dpressure':
        convertback = beta_ave*conversions[1] # this is the pv component
    elif type == 'dmu-constB' or type == 'dbeta-dmu':
        convertback = -1*kB
    else:
        convertback = kB
        
    w_F = (beta_k[1]-beta_k[0])*U_kn[0,0:N_k[0]]
    w_R = -((beta_k[1]-beta_k[0])*U_kn[1,0:N_k[1]])

    if (type == 'dbeta-constP') or (type == 'dpressure-constB') or (type == 'dbeta-dpressure'):
        w_F += conversions[1]*(beta_k[1]*P_k[1]-beta_k[0]*P_k[0])*V_kn[0,0:N_k[0]]       
        w_R += -conversions[1]*(beta_k[1]*P_k[1]-beta_k[0]*P_k[0])*V_kn[1,0:N_k[1]]       

    if (type == 'dbeta-constmu') or (type == 'dmu-constB') or (type == 'dbeta-dmu'):  
        w_F += -(beta_k[1]*mu_k[1]-beta_k[0]*mu_k[0])*N_kn[0,0:N_k[0]]       
        w_R += (beta_k[1]*mu_k[1]-beta_k[0]*mu_k[0])*N_kn[1,0:N_k[1]]       
        # it's not entirely clear if this right because of lack of overlap in phase space between different N's! 
        
    if not quiet:
        print("Now computing log of partition functions using BAR")
    
    (df,ddf) = BAR(w_F,w_R)

    if not quiet:
        print("using %.5f for log of partition functions computed from BAR" % (df))
        print("Uncertainty in quantity is %.5f" % (ddf))
        print("Assuming this is negligible compared to sampling error at individual points")

    if (bMaxwell):  # only applies for kinetic energies
        if not quiet:
            print("Now fitting to a Maxwell-Boltzmann distribution")
        for k in range(2):
            fn = filename + '_maxboltz' + str(T_k[k])
            MaxwellBoltzFit(bins,U_kn[k,0:N_k[k]],N_k[k],kB*T_k[k],fn)

    quant = {}

    if (bLinearFit and check_twodtype(type)):
        if not quiet:
            print("Now computing the linear fit parameters")
        if filename is not None:
            fn = filename + '_linear'
        else:
            fn = None
        (fitvals,dfitvals) = LinFit(bins, N_k, dp, const, v, df=df,
                                    analytic_uncertainty=True, name='', g=g,
                                    filename=fn, screen=screen)
        slope = fitvals[1]
        dslope = dfitvals[1]
        quant['linear'] = [abs((slope - dp[0])/dslope)]
        if not quiet:
            Print1DStats('Linear Fit Analysis (analytical error)',type,fitvals,convertback,dp,const,dfitvals=dfitvals)

    if (bNonLinearFit and check_twodtype(type)):
        if not quiet:
            print("Now computing the nonlinear fit parameters")
        fn = filename + '_nonlinear'
        (fitvals,dfitvals) = NonLinFit(bins,N_k,dp,const,v,df=df,name='',figname=fn,bGraph=do_plot,analytic_uncertainty=True,g=g,type=type,vunits=vunits)
        slope = fitvals[1]
        dslope = dfitvals[1]
        quant['nonlinear'] = [abs((slope - dp[0])/dslope)]
        if not quiet:
            Print1DStats('Nonlinear Fit Analysis (analytical error)',type,fitvals,convertback,dp,const,dfitvals=dfitvals)

    if (bMaxLikelihood):
        if not quiet:
            print("Now computing the maximum likelihood parameters")
        (fitvals,dfitvals) = MaxLikeParams(N_k,dp,const,v,df=df,analytic_uncertainty=True,g=numpy.average(g),quiet=quiet)
        if (check_twodtype(type)):
            if not quiet:
                Print1DStats('Maximum Likelihood Analysis (analytical error)',type,fitvals,convertback,dp,const,dfitvals=dfitvals)
            slope = fitvals[1]
            dslope = dfitvals[1]
            quant['maxLikelihood'] = [abs((slope - dp[0])/dslope)]
        else:
            if not quiet:
                Print2DStats('2D-Maximum Likelihood Analysis (analytical error)',type,fitvals,kB,convertback,dp,const,dfitvals=dfitvals)
            eneslope = fitvals[1]
            pvslope = fitvals[2]
            deneslope = dfitvals[1]
            dpvslope = dfitvals[2]
            quant['maxLikelihood'] = [abs((eneslope - dp[0])/deneslope),
                                      abs((pvslope - dp[1]) / dpvslope)]

    returnvalue = None
    if bMaxLikelihood:
        returnvalue = quant['maxLikelihood']
    elif bLinearFit:
        returnvalue = quant['linear']
    elif bNonLinearFit:
        returnvalue = quant['nonlinear']

    if (reptype is None):
        return returnvalue

    if (reptype == 'bootstrap'):
        nreps = nboots
        if (nreps < 50):
            if (nreps > 1):
                print("Warning, less than 50 bootstraps (%d requested) is likely not a good statistical idea" % (nreps))
            else:
                print("Cannot provide bootstrap statisics, only %d requested" % (nreps))
                return returnvalue

        print("Now bootstrapping (n=%d) for uncertainties . . . could take a bit of time!" % (nreps))
    elif (reptype == 'independent'):
        nreps = len(reps)
        print("Now analyzing %d independent samples . . . could take a bit of time!" % (nreps))
    else:
        quit("Don't understand reptype = %s; quitting" % (reptype))

    if check_twodtype(type):   # how many values do we have to deal with?
        rval = 2
    else:
        rval = 3

    linvals = numpy.zeros([rval,nreps],float)
    nlvals = numpy.zeros([rval,nreps],float)
    mlvals = numpy.zeros([rval,nreps],float)

    for n in range(nreps):
        if n % 10 == 0 and not quiet:
            print("Finished %d samples . . ." % (n))
        
        if (reptype == 'bootstrap'):    
            for k in range(K):
                if ((g is None) or (g[0]==1 and g[1]==1)):
                    #do normal bootstrap
                    rindex = numpy.random.randint(0,high=N_k[k],size=N_k[k]);  # bootstrap it 
                    for i in range(len(const)):
                        vr[i,k,0:N_k[k]] = v[i,k,rindex]
                else:
                    # we are given correlation times.  Do block bootstrapping.
                    gk = int(numpy.ceil(g[k]))
                    nblocks = int(numpy.floor(N_k[k]/gk))
                    # moving blocks bootstrap; all contiguous segments of length gk

                    rindex = numpy.random.randint(0,high=gk*(nblocks-1),size=nblocks); 
                    for nb in range(nblocks):
                        for i in range(len(const)):
                            vr[i,k,nb*gk:(nb+1)*gk] = v[i,k,rindex[nb]:rindex[nb]+gk]
                    N_k[k] = nblocks*gk  # we could have a few samples less now
    
        if (reptype == 'independent'):
            for k in range(K):
                for i in range(len(const)):
                    vr[i,k,0:N_k[k]] = reps[n][i][k,0:N_k[k]] 
            
        if (bLinearFit and check_twodtype(type)):    
            fitvals = LinFit(bins,N_k,dp,const,vr) 
            for i in range(rval):
                linvals[i,n] = fitvals[0][i]

        if (bNonLinearFit and check_twodtype(type)):
            fitvals = NonLinFit(bins,N_k,dp,const,vr,df=df)
            for i in range(rval):
                nlvals[i,n] = fitvals[0][i]

        if (bMaxLikelihood):
            fitvals = MaxLikeParams(N_k,dp,const,vr,df=df,quiet=quiet)
            for i in range(rval):
                mlvals[i,n] = fitvals[0][i]

    if (bLinearFit and check_twodtype(type)):
        if not quiet:
            Print1DStats('Linear Fit Analysis',type,[linvals[0],linvals[1]],convertback,dp,const)
        slopes = linvals[1]
        slope = numpy.average(slopes)
        dslope = numpy.std(slopes)
        quant['linear'] = [abs((slope - dp[0])/dslope)]
        
    if (bNonLinearFit and check_twodtype(type)):
        if not quiet:
            Print1DStats('Nonlinear Fit Analysis',type,[nlvals[0],nlvals[1]],convertback,dp,const)
        slopes = nlvals[1]
        slope = numpy.average(slopes)
        dslope = numpy.std(slopes)
        quant['nonlinear'] = [abs((slope - dp[0]) / dslope)]

    if (bMaxLikelihood):
        if check_twodtype(type):
            if not quiet:
                Print1DStats('Maximum Likelihood Analysis',type,[mlvals[0],mlvals[1]],convertback,dp,const)
            slopes = nlvals[1]
            slope = numpy.average(slopes)
            dslope = numpy.std(slopes)
            quant['maxLikelihood'] = [abs((slope - dp[0]) / dslope)]
        else:
            if not quiet:
                Print2DStats('2D-Maximum Likelihood Analysis',type,[mlvals[0],mlvals[1],mlvals[2]],kB,convertback,dp,const)
            eneslopes = mlvals[1]
            pvslopes = mlvals[2]
            eneslope = numpy.average(eneslopes)
            deneslope = numpy.average(deneslopes)
            pvslope = numpy.average(pvslopes)
            dpvslope = numpy.average(dpvslopes)
            quant['maxLikelihood'] = [abs((eneslope - dp[0]) / deneslope),
                                      abs((pvslope - dp[1]) / dpvslope)]

    returnvalue = None
    if bMaxLikelihood:
        returnvalue = quant['maxLikelihood']
    elif bLinearFit:
        returnvalue = quant['linear']
    elif bNonLinearFit:
        returnvalue = quant['nonlinear']

    return returnvalue
    
# to do: fix the drawing directions so that correct data has the legend in the right place.
