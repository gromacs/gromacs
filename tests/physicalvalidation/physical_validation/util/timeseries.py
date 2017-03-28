#!/usr/local/bin/env python

"""
A module for extracting uncorrelated samples from correlated timeseries data.

This module provides various tools that allow one to examine the correlation functions and
integrated autocorrelation times in correlated timeseries data, compute statistical inefficiencies,
and automatically extract uncorrelated samples for data analysis.

REFERENCES

[1] Shirts MR and Chodera JD. Statistically optimal analysis of samples from multiple equilibrium states.
J. Chem. Phys. 129:124105, 2008
http://dx.doi.org/10.1063/1.2978177

[2] J. D. Chodera, W. C. Swope, J. W. Pitera, C. Seok, and K. A. Dill. Use of the weighted
histogram analysis method for the analysis of simulated and parallel tempering simulations.
JCTC 3(1):26-41, 2007.

"""
from __future__ import print_function
from __future__ import division

#=============================================================================================
# COPYRIGHT NOTICE
#
# Written by John D. Chodera <jchodera@gmail.com> and Michael R. Shirts <mrshirts@gmail.com>.
#
# Copyright (c) 2007 The Regents of the University of California.  All Rights Reserved.
# Portions of this software are Copyright (c) 2007 Stanford University and Columbia University.
#
# This program is free software; you can redistribute it and/or modify it under the terms of
# the GNU General Public License as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with this program;
# if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.
#=============================================================================================

#=============================================================================================
# TODO
# * Implement unit tests that generate timeseries with various levels of Gaussian correlation to test all methods.
# * Add Zwanzig procedure for estimating statistical uncertainties in correlation functions
# (by making Gaussian process assumptions).
#=============================================================================================

#=============================================================================================
# VERSION CONTROL INFORMATION
#=============================================================================================


__version__ = "2.0beta"
__authors__ = "Michael R. Shirts and John D. Chodera."
__licanse__ = "GPL 2.0"

#=============================================================================================
# IMPORTS
#=============================================================================================

import math
import numpy 
import numpy.linalg
import sys

#=============================================================================================
# Exception class.
#=============================================================================================

class ParameterError(Exception):
  """An error in the input parameters has been detected.

  """

#=============================================================================================
# Issue warning on import.
#=============================================================================================

LongWarning = "Warning: If the inherent timescales of the system are long compared to those being analyzed, this statistical inefficiency may be an underestimate.  The estimate presumes the use of many statistically independent samples.  Tests should be performed to assess whether this condition is satisfied.   Be cautious in the interpretation of the data."

# sys.stderr.write(LongWarning + '\n')

#=============================================================================================
# METHODS
#=============================================================================================

#=============================================================================================
def statisticalInefficiency(A_n, B_n=None, fast=False, mintime=3):
  """
  Compute the (cross) statistical inefficiency of (two) timeseries.

  REQUIRED ARGUMENTS  
    A_n (numpy array) - A_n[n] is nth value of timeseries A.  Length is deduced from vector.

  OPTIONAL ARGUMENTS
    B_n (numpy array) - B_n[n] is nth value of timeseries B.  Length is deduced from vector.
       If supplied, the cross-correlation of timeseries A and B will be estimated instead of the
       autocorrelation of timeseries A.  
    fast (boolean) - if True, will use faster (but less accurate) method to estimate correlation
       time, described in Ref. [1] (default: False)
    mintime (int) - minimum amount of correlation function to compute (default: 3)
       The algorithm terminates after computing the correlation time out to mintime when the
       correlation function furst goes negative.  Note that this time may need to be increased
       if there is a strong initial negative peak in the correlation function.

  RETURNS
    g is the estimated statistical inefficiency (equal to 1 + 2 tau, where tau is the correlation time).
       We enforce g >= 1.0.

  NOTES 
    The same timeseries can be used for both A_n and B_n to get the autocorrelation statistical inefficiency.
    The fast method described in Ref [1] is used to compute g.

  REFERENCES  
    [1] J. D. Chodera, W. C. Swope, J. W. Pitera, C. Seok, and K. A. Dill. Use of the weighted
    histogram analysis method for the analysis of simulated and parallel tempering simulations.
    JCTC 3(1):26-41, 2007.

  EXAMPLES

  Compute statistical inefficiency of timeseries data with known correlation time.  

  >>> import testsystems
  >>> A_n = testsystems.generateCorrelatedTimeseries(N=100000, tau=5.0)
  >>> g = statisticalInefficiency(A_n, fast=True)
  
  """

  # Create numpy copies of input arguments.
  A_n = numpy.array(A_n)
  if B_n is not None:  
    B_n = numpy.array(B_n)
  else:
    B_n = numpy.array(A_n) 
  
  # Get the length of the timeseries.
  N = A_n.size

  # Be sure A_n and B_n have the same dimensions.
  if(A_n.shape != B_n.shape):
    raise ParameterError('A_n and B_n must have same dimensions.')

  # Initialize statistical inefficiency estimate with uncorrelated value.
  g = 1.0
    
  # Compute mean of each timeseries.
  mu_A = A_n.mean()
  mu_B = B_n.mean()

  # Make temporary copies of fluctuation from mean.
  dA_n = A_n.astype(numpy.float64) - mu_A
  dB_n = B_n.astype(numpy.float64) - mu_B

  # Compute estimator of covariance of (A,B) using estimator that will ensure C(0) = 1.
  sigma2_AB = (dA_n * dB_n).mean() # standard estimator to ensure C(0) = 1

  # Trap the case where this covariance is zero, and we cannot proceed.
  if(sigma2_AB == 0):
    raise ParameterException('Sample covariance sigma_AB^2 = 0 -- cannot compute statistical inefficiency')

  # Accumulate the integrated correlation time by computing the normalized correlation time at
  # increasing values of t.  Stop accumulating if the correlation function goes negative, since
  # this is unlikely to occur unless the correlation function has decayed to the point where it
  # is dominated by noise and indistinguishable from zero.
  t = 1
  increment = 1
  while (t < N-1):
    # compute normalized fluctuation correlation function at time t
    C = numpy.sum( dA_n[0:(N-t)]*dB_n[t:N] + dB_n[0:(N-t)]*dA_n[t:N] ) / (2.0 * float(N-t) * sigma2_AB)
    # Terminate if the correlation function has crossed zero and we've computed the correlation
    # function at least out to 'mintime'.
    if (C <= 0.0) and (t > mintime):
      break
    
    # Accumulate contribution to the statistical inefficiency.
    g += 2.0 * C * (1.0 - float(t)/float(N)) * float(increment)
    # Increment t and the amount by which we increment t.
    t += increment

    # Increase the interval if "fast mode" is on.
    if fast: increment += 1

  # g must be at least unity
  if (g < 1.0): g = 1.0
   
  # Return the computed statistical inefficiency.
  return g
#=============================================================================================
def statisticalInefficiencyMultiple(A_kn, fast=False, return_correlation_function=False):
  """
  Estimate the statistical inefficiency from multiple stationary timeseries (of potentially differing lengths).

  REQUIRED ARGUMENTS  
    A_kn (Python list of numpy arrays) - A_kn[k] is the kth timeseries, and A_kn[k][n] is nth value of timeseries k.  Length is deduced from arrays.

  OPTIONAL ARGUMENTS  
    fast can be set to True to give a less accurate but very quick estimate (default False)
    return_correlation_function - if True, will also return estimates of normalized fluctuation correlation function that were computed (default: False)

  RETURNS
    g is the statistical inefficiency (equal to 1 + 2 tau, where tau is the integrated autocorrelation time).
    Ct (list of tuples) - Ct[n] = (t, C) with time t and normalized correlation function estimate C is returned as well if return_correlation_function is set to True
    
  NOTES 
    The autocorrelation of the timeseries is used to compute the statistical inefficiency.
    The normalized fluctuation autocorrelation function is computed by averaging the unnormalized raw correlation functions.
    The fast method described in Ref [1] is used to compute g.

  REFERENCES  
    [1] J. D. Chodera, W. C. Swope, J. W. Pitera, C. Seok, and K. A. Dill. Use of the weighted
    histogram analysis method for the analysis of simulated and parallel tempering simulations.
    JCTC 3(1):26-41, 2007.

  EXAMPLES

  Estimate statistical efficiency from multiple timeseries of different lengths.

  >>> import testsystems
  >>> N_k = [1000, 2000, 3000, 4000, 5000]
  >>> tau = 5.0 # exponential relaxation time
  >>> A_kn = [ testsystems.generateCorrelatedTimeseries(N=N, tau=tau) for N in N_k ]
  >>> g = statisticalInefficiencyMultiple(A_kn)

  Also return the values of the normalized fluctuation autocorrelation function that were computed.

  >>> [g, Ct] = statisticalInefficiencyMultiple(A_kn, return_correlation_function=True)

  """

  # Convert A_kn into a list of arrays if it is not in this form already.
  if (type(A_kn) == numpy.ndarray):
    A_kn_list = list()
    if A_kn.ndim == 1:      
      A_kn_list.append(A_kn.copy())
    else:
      [K,N] = A_kn.shape
      for k in range(K):
        A_kn_list.append(A_kn[k,:].copy())
    A_kn = A_kn_list
      
  # Determine number of timeseries.
  K = len(A_kn)

  # Get the length of each timeseries.
  N_k = numpy.zeros([K], numpy.int32)
  for k in range(K):
    N_k[k] = A_kn[k].size

  # Compute average timeseries length.
  Navg = numpy.array(N_k, numpy.float64).mean()

  # Determine total number of samples.
  N = numpy.sum(N_k)

  # Initialize statistical inefficiency estimate with uncorrelated value.
  g = 1.0
    
  # Compute sample mean.
  mu = 0.0
  for k in range(K):
    mu += A_kn[k].sum()
  mu /= float(N)
  
  # Construct and store fluctuation timeseries.
  dA_kn = list()
  for k in range(K):
    dA_n = A_kn[k] - mu
    dA_kn.append(dA_n.copy())

  # Compute sample variance from mean of squared fluctuations, to ensure that C(0) = 1.
  sigma2 = 0.0
  for k in range(K):
    sigma2 += (dA_kn[k]**2).sum()
  sigma2 /= float(N)
  
  # Initialize statistical inefficiency estimate with uncorrelated value.
  g = 1.0

  # Initialize storage for correlation function.
  Ct = list() # Ct[n] is a tuple (t, C) of the time lag t and estimate of normalized fluctuation correlation function C
    
  # Accumulate the integrated correlation time by computing the normalized correlation time at
  # increasing values of t.  Stop accumulating if the correlation function goes negative, since
  # this is unlikely to occur unless the correlation function has decayed to the point where it
  # is dominated by noise and indistinguishable from zero.
  t = 1
  increment = 1  
  while (t < N_k.max()-1):
    # compute unnormalized correlation function
    numerator = 0.0
    denominator = 0.0
    for k in range(K):
      if (t >= N_k[k]): continue # skip trajectory if lag time t is greater than its length
      dA_n = dA_kn[k] # retrieve trajectory
      x = dA_n[0:(N_k[k]-t)] * dA_n[t:N_k[k]]
      numerator += x.sum() # accumulate contribution from trajectory k
      denominator += float(x.size) # count how many overlapping time segments we've included

    C = numerator / denominator

    # compute normalized fluctuation correlation function at time t
    C = C / sigma2
    #print "C[%5d] = %16f (%16f / %16f)" % (t, C, numerator, denominator)    

    # Store estimate of correlation function.
    Ct.append( (t,C) )

    # Terminate if the correlation function has crossed zero.
    # Note that we've added a hack (t > 10) condition to avoid terminating too early in correlation functions that have a strong negative peak at 
    if (C <= 0.0) and (t > 10):
      break
  
    # Accumulate contribution to the statistical inefficiency.
    g += 2.0 * C * (1.0 - float(t)/Navg) * float(increment)

    # Increment t and the amount by which we increment t.
    t += increment

    # Increase the interval if "fast mode" is on.
    if fast: increment += 1

  # g must be at least unity
  if (g < 1.0): g = 1.0

  # Return statistical inefficency and correlation function estimate, if requested.
  if return_correlation_function:
    return (g, Ct)

  # Return the computed statistical inefficiency.
  return g
#=============================================================================================
def integratedAutocorrelationTime(A_n, B_n=None, fast=False, mintime=3):
  """
  Estimate the integrated autocorrelation time.

  """
  
  g = statisticalInefficiency(A_n, B_n, fast, mintime)
  tau = (g-1.0)/2.0
  return tau
#=============================================================================================
def integratedAutocorrelationTimeMultiple(A_kn, fast=False):
  """
  Estimate the integrated autocorrelation time from multiple timeseries.

  """
  
  g = statisticalInefficiencyMultiple(A_kn, fast, False)
  tau = (g-1.0)/2.0
  return tau
#=============================================================================================
def normalizedFluctuationCorrelationFunction(A_n, B_n=None, N_max=None):
  """
  Compute the normalized fluctuation (cross) correlation function of (two) stationary timeseries.

  C(t) = (<A(t) B(t)> - <A><B>) / (<AB> - <A><B>)

  This may be useful in diagnosing odd time-correlations in timeseries data.

  REQUIRED ARGUMENTS  
    A_n[n] is nth value of timeseries A.  Length is deduced from vector.
    B_n[n] is nth value of timeseries B.  Length is deduced from vector.

  OPTIONAL ARGUMENTS
    N_max - if specified, will only compute correlation function out to time lag of N_max

  RETURNS
    C_n[n] is the normalized fluctuation auto- or cross-correlation function for timeseries A(t) and B(t).

  NOTES 
    The same timeseries can be used for both A_n and B_n to get the autocorrelation statistical inefficiency.
    This procedure may be slow.
    The statistical error in C_n[n] will grow with increasing n.  No effort is made here to estimate the uncertainty.

  REFERENCES  
    [1] J. D. Chodera, W. C. Swope, J. W. Pitera, C. Seok, and K. A. Dill. Use of the weighted
    histogram analysis method for the analysis of simulated and parallel tempering simulations.
    JCTC 3(1):26-41, 2007.

  EXAMPLES

  Estimate normalized fluctuation correlation function.

  >>> import testsystems
  >>> A_t = testsystems.generateCorrelatedTimeseries(N=10000, tau=5.0)
  >>> C_t = normalizedFluctuationCorrelationFunction(A_t, N_max=25)

  """

  # If B_n is not specified, set it to be identical to A_n.
  if B_n is None:
    B_n = A_n

  # Create numpy copies of input arguments.
  A_n = numpy.array(A_n)
  B_n = numpy.array(B_n)

  # Get the length of the timeseries.
  N = A_n.size

  # Set maximum time to compute correlation functon for.
  if (not N_max) or (N_max > N-1):
    N_max = N-1

  # Be sure A_n and B_n have the same dimensions.
  if(A_n.shape != B_n.shape):
    raise ParameterError('A_n and B_n must have same dimensions.')

  # Initialize statistical inefficiency estimate with uncorrelated value.
  g = 1.0
    
  # Compute means and variance.
  mu_A = A_n.mean()
  mu_B = B_n.mean()

  # Make temporary copies at high precision with means subtracted off.
  dA_n = A_n.astype(numpy.float64) - mu_A
  dB_n = B_n.astype(numpy.float64) - mu_B

  # sigma2_AB = sum((A_n-mu_A) * (B_n-mu_B)) / (float(N)-1.0) # unbiased estimator
  sigma2_AB = (dA_n * dB_n).mean() # standard estimator to ensure C(0) = 1
  if(sigma2_AB == 0):
    raise ParameterException('Sample covariance sigma_AB^2 = 0 -- cannot compute statistical inefficiency')

  # allocate storage for normalized fluctuation correlation function
  C_n = numpy.zeros([N_max+1], numpy.float64)  

  # Compute normalized correlation funtion.
  t = 0
  for t in range(0,N_max+1):
    # compute normalized fluctuation correlation function at time t
    C_n[t] = sum( dA_n[0:(N-t)]*dB_n[t:N] + dB_n[0:(N-t)]*dA_n[t:N] ) / (2.0 * float(N-t) * sigma2_AB)

  # Return the computed correlation function
  return C_n
#=============================================================================================
def normalizedFluctuationCorrelationFunctionMultiple(A_kn, B_kn=None, N_max=None, suppress_warning=False):
  """
  Compute the normalized fluctuation (cross) correlation function of (two) timeseries from multiple timeseries samples.

  C(t) = (<A(t) B(t)> - <A><B>) / (<AB> - <A><B>)

  This may be useful in diagnosing odd time-correlations in timeseries data.

  REQUIRED ARGUMENTS  
    A_kn (Python list of numpy arrays) - A_kn[k] is the kth timeseries, and A_kn[k][n] is nth value of timeseries k.  Length is deduced from arrays.
    B_kn (Python list of numpy arrays) - B_kn[k] is the kth timeseries, and B_kn[k][n] is nth value of timeseries k.  B_kn[k] must have same length as A_kn[k]

  OPTIONAL ARGUMENTS
    N_max - if specified, will only compute correlation function out to time lag of N_max
    suppress_warning - if we are calculating a lot of these, the warning could get a little annoying. Make it possible to suppress it, but don't make that the default.

  RETURNS
    C_n[n] is the normalized fluctuation auto- or cross-correlation function for timeseries A(t) and B(t).

  NOTES
    The same timeseries can be used for both A_n and B_n to get the autocorrelation statistical inefficiency.
    This procedure may be slow.
    The statistical error in C_n[n] will grow with increasing n.  No effort is made here to estimate the uncertainty.

  REFERENCES  
    [1] J. D. Chodera, W. C. Swope, J. W. Pitera, C. Seok, and K. A. Dill. Use of the weighted
    histogram analysis method for the analysis of simulated and parallel tempering simulations.
    JCTC 3(1):26-41, 2007.

  EXAMPLES

  Estimate a portion of the normalized fluctuation autocorrelation function from multiple timeseries of different length.

  >>> import testsystems
  >>> N_k = [1000, 2000, 3000, 4000, 5000]
  >>> tau = 5.0 # exponential relaxation time
  >>> A_kn = [ testsystems.generateCorrelatedTimeseries(N=N, tau=tau) for N in N_k ]
  >>> C_n = normalizedFluctuationCorrelationFunctionMultiple(A_kn, N_max=25)

  """

  # If B_kn is not specified, define it to be identical with A_kn.
  if B_kn is None:
    B_kn = A_kn  

  # TODO: Change this to support other iterable types, like sets.
  # Make sure A_kn and B_kn are both lists
  if (type(A_kn) is not list) or (type(B_kn) is not list):
    raise ParameterError("A_kn and B_kn must each be a list of numpy arrays.")
  
  # Ensure the same number of timeseries are stored in A_kn and B_kn.
  if (len(A_kn) != len(B_kn)):
    raise ParameterError("A_kn and B_kn must contain corresponding timeseries -- different numbers of timeseries detected in each.")

  # Determine number of timeseries stored.
  K = len(A_kn)  

  # Ensure both observable trajectories in each timeseries are of the same length.
  for k in range(K):
    A_n = A_kn[k]
    B_n = B_kn[k]
    if A_n.size != B_n.size:
      raise Exception("A_kn and B_kn must contain corresponding timeseries -- lack of correspondence in timeseries lenghts detected.")
  
  # Get the length of each timeseries.
  N_k = numpy.zeros([K], numpy.int32)
  for k in range(K):
    N_k[k] = A_kn[k].size

  # Determine total number of samples.
  N = sum(N_k)

  # Set maximum time to compute correlation functon for.
  if (not N_max) or (N_max > max(N_k) - 1):
    N_max = max(N_k) - 1

  # Compute means.
  mu_A = 0.0
  mu_B = 0.0
  for k in range(K):
    mu_A += A_kn[k].sum()
    mu_B += B_kn[k].sum()
  mu_A /= float(N)
  mu_B /= float(N)
  
  # Compute fluctuation timeseries.
  dA_kn = list()
  dB_kn = list()
  for k in range(K):
    dA_n = A_kn[k] - mu_A
    dB_n = B_kn[k] - mu_B
    dA_kn.append(dA_n)
    dB_kn.append(dB_n)    

  # Compute covariance.
  sigma2_AB = 0.0
  for k in range(K):
    sigma2_AB += (dA_kn[k] * dB_kn[k]).sum()
  sigma2_AB /= float(N)

  # allocate storage for normalized fluctuation correlation function
  C_n = numpy.zeros([N_max+1], numpy.float64)  

  # Accumulate the integrated correlation time by computing the normalized correlation time at
  # increasing values of t.  Stop accumulating if the correlation function goes negative, since
  # this is unlikely to occur unless the correlation function has decayed to the point where it
  # is dominated by noise and indistinguishable from zero.
  t = 0
  for t in range(0,N_max+1):
    # compute unnormalized correlation function
    numerator = 0.0
    denominator = 0.0
    for k in range(K):
      if (t >= N_k[k]): continue # skip this trajectory if t is longer than the timeseries
      numerator += (dA_kn[k][0:(N_k[k]-t)] * dB_kn[k][t:N_k[k]]).sum()
      denominator += float(N_k[k]-t)
    C = numerator / denominator
    
    # compute normalized fluctuation correlation function at time t
    C /= sigma2_AB

    # Store correlation function.
    C_n[t] = C

  # Return the computed fluctuation correlation function.
  return C_n
#=============================================================================================
def subsampleCorrelatedData(A_t, g=None, fast=False, conservative=False, verbose=False):
  """Determine the indices of an uncorrelated subsample of the data.

  REQUIRED ARGUMENTS  
    A_t (T array) - A_t[t] is the t-th value of timeseries A(t).  Length is deduced from vector.

  OPTIONAL ARGUMENTS
    g (float) - if provided, the statistical inefficiency g is used to subsample the timeseries -- otherwise it will be computed (default: None)
    fast (logical) - fast can be set to True to give a less accurate but very quick estimate (default: False)
    conservative (logical) - if set to True, uniformly-spaced indices are chosen with interval ceil(g), where g is the statistical inefficiency.  Otherwise, indices are chosen non-uniformly with interval of approximately g in order to end up with approximately T/g total indices
    verbose (logical) - if True, some output is printed

  RETURNS  
    indices (list of int) - the indices of an uncorrelated subsample of the data

  NOTES
    The statistical inefficiency is computed with the function computeStatisticalInefficiency().

  TODO
    Instead of using regular stride, use irregular stride so more data can be fit in when g is non-integral.

  EXAMPLES

  Subsample a correlated timeseries to extract an effectively uncorrelated dataset.

  >>> import testsystems
  >>> A_t = testsystems.generateCorrelatedTimeseries(N=10000, tau=5.0) # generate a test correlated timeseries
  >>> indices = subsampleCorrelatedData(A_t) # compute indices of uncorrelated timeseries
  >>> A_n = A_t[indices] # extract uncorrelated samples

  Extract uncorrelated samples from multiple timeseries data from the same process.

  >>> # Generate multiple correlated timeseries data of different lengths.
  >>> T_k = [1000, 2000, 3000, 4000, 5000]
  >>> K = len(T_k) # number of timeseries
  >>> tau = 5.0 # exponential relaxation time
  >>> A_kt = [ testsystems.generateCorrelatedTimeseries(N=T, tau=tau) for T in T_k ] # A_kt[k] is correlated timeseries k
  >>> # Estimate statistical inefficiency from all timeseries data.
  >>> g = statisticalInefficiencyMultiple(A_kt)
  >>> # Count number of uncorrelated samples in each timeseries.
  >>> N_k = numpy.array([ len(subsampleCorrelatedData(A_t, g=g)) for A_t in A_kt ]) # N_k[k] is the number of uncorrelated samples in timeseries k
  >>> N = N_k.sum() # total number of uncorrelated samples
  >>> # Subsample all trajectories to produce uncorrelated samples
  >>> A_kn = [ A_t[subsampleCorrelatedData(A_t, g=g)] for A_t in A_kt ] # A_kn[k] is uncorrelated subset of trajectory A_kt[t]
  >>> # Concatenate data into one timeseries.
  >>> A_n = numpy.zeros([N], numpy.float32) # A_n[n] is nth sample in concatenated set of uncorrelated samples
  >>> A_n[0:N_k[0]] = A_kn[0]
  >>> for k in range(1,K): A_n[N_k[0:k].sum():N_k[0:k+1].sum()] = A_kn[k]

  """

  # Create numpy copy of arrays.
  A_t = numpy.array(A_t)

  # Get the length of the timeseries.
  T = A_t.size

  # Compute the statistical inefficiency for the timeseries.
  if not g:
    if verbose: print("Computing statistical inefficiency...")
    g = statisticalInefficiency(A_t, A_t, fast = fast)
    if verbose: print("g = %f" % g)

  if conservative:
    # Round g up to determine the stride we can use to pick out regularly-spaced uncorrelated samples.
    import math
    stride = int(math.ceil(g))
    if verbose: print("conservative subsampling: using stride of %d" % stride)
    
    # Assemble list of indices of uncorrelated snapshots.
    indices = list(range(0, T, stride))
  else:
    # Choose indices as floor(n*g), with n = 0,1,2,..., until we run out of data.
    import math
    indices = []
    n = 0
    while int(round(n*g)) < T:
      t = int(round(n*g))
      # ensure we don't sample the same point twice
      if (n == 0) or (t != indices[n-1]):
        indices.append(t)
      n += 1
    if verbose: print("standard subsampling: using average stride of %f" % g)

  # Number of samples in subsampled timeseries.
  N = len(indices)
  
  if verbose: print("The resulting subsampled set has %d samples (original timeseries had %d)." % (N, T))

  # Return the list of indices of uncorrelated snapshots.
  return indices

#=============================================================================================
# MAIN AND TESTS
#=============================================================================================

if __name__ == "__main__":
  import doctest
  doctest.testmod()

