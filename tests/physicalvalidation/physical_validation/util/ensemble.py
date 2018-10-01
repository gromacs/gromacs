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
This file reimplements most functionality of the checkensemble.py code
originally published on https://github.com/shirtsgroup/checkensemble. It
serves as the low-level functionality of the high-level module
:mod:`physical_validation.ensemble`.
"""
from __future__ import division
import numpy as np
import scipy.optimize

import pymbar

from . import trajectory
from . import error as pv_error
from . import plot


def generate_histograms(traj1, traj2, g1, g2, bins):

    n1 = np.size(traj1)
    n2 = np.size(traj2)

    h1 = np.histogram(traj1, bins=bins)[0]/n1
    h2 = np.histogram(traj2, bins=bins)[0]/n2
    dh1 = np.sqrt(g1 * h1 * (1 - h1) / n1)
    dh2 = np.sqrt(g2 * h2 * (1 - h2) / n2)

    return h1, h2, dh1, dh2


def do_linear_fit(traj1, traj2, g1, g2, bins,
                  screen=False, filename=None,
                  trueslope=0.0, trueoffset=0.0,
                  units=None):

    h1, h2, dh1, dh2 = generate_histograms(traj1, traj2, g1, g2, bins)

    #  v  copied from checkensemble.py  v
    ratio = np.log(h2 / h1)
    dratio = np.sqrt((dh1/h1)**2 + (dh2/h2)**2)

    usedat = np.isfinite(ratio)
    y = ratio[usedat]
    nuse = len(y)
    weights = 1.0/dratio[usedat]

    xaxis = (bins[:-1] + bins[1:])/2
    x = xaxis[usedat]

    x_mat = np.ones([nuse, 2])
    x_mat[:, 1] = x

    w = np.diag(weights) 
    wx = np.dot(w, x_mat)
    wy = np.dot(w, y)
    wx_t = np.transpose(wx)
    z = np.dot(wx_t, wx)
    wxy = np.dot(wx_t, wy)

    a = np.linalg.solve(z, wxy)
    da_matrix = np.transpose(np.linalg.inv(z))
    da = np.zeros(2)
    da[0] = np.sqrt(da_matrix[0, 0])
    da[1] = np.sqrt(da_matrix[1, 1])

    # the true line is y = df + dp*x, where y is ln P_1(X)/P_2(X)
    #  ^  end copied from checkensemble.py  ^

    do_plot = screen or filename is not None
    if do_plot:
        true = trueoffset+trueslope*xaxis
        fit = a[0] + a[1]*xaxis

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

        if units is not None:
            units = ' [' + units + ']'
        else:
            units = ''

        annot = ('{:.1f}'.format(abs((a[1] - trueslope) / da[1])) +
                 ' quantiles')

        plot.plot(data,
                  legend='best',
                  title='Log probability ratio',
                  xlabel='Energy' + units,
                  ylabel=r'$\log\frac{P_2(E)}{P_1(E)}$',
                  filename=filename,
                  screen=screen,
                  axtext=annot)

    return a, da


def do_max_likelihood_fit(traj1, traj2, g1, g2,
                          init_params=None,
                          verbose=False):

    # ============================================================= #
    # Define (negative) log-likelihood function and its derivatives #
    # ============================================================= #
    def log_likelihood(a, ene1, ene2):
        # Returns negative of eq (8) of check_ensemble paper
        #
        # Uses log (1/f(x)) == -log(f(x))
        # and log(1 + e^x) == log(e^x (e^-x + 1)) == x + log(1 + e^-x)
        #     ^(a)                                   ^(b)
        # form (a) -> 0 for x->-inf, -> inf for x->inf
        # form (b) -> NaN for x->-inf, -> x for x->inf
        # combined: -> 0 for x-> -inf, -> x for x-> inf
        def log_1_plus_exp(y):
            def f(yy):
                with np.errstate(over='raise'):
                    try:
                        xx = np.log(1 + np.exp(yy))
                    except FloatingPointError:
                        xx = yy + np.log(1 + np.exp(-yy))
                    return xx
            return np.vectorize(f)(y)

        if a.size == 2:
            return (np.sum(log_1_plus_exp(a[0] + a[1]*ene1)) +
                    np.sum(log_1_plus_exp(-a[0] - a[1]*ene2)))
        else:
            return (np.sum(log_1_plus_exp(a[0] + a[1]*ene1[0] + a[2]*ene1[1])) +
                    np.sum(log_1_plus_exp(-a[0] - a[1]*ene2[0] - a[2]*ene2[1])))

    def da_log_likelihood(a, ene1, ene2):
        # Returns the first derivative wrt the parameters a of log_likelihood
        #
        # d/da0 log(1 + exp(a0 + a1*E)) == exp(a0 + a1*E) / (1 + exp(a0 + a1*E))
        #                               == 1 / (1 + exp(-a0 - a1*E))
        # d/da1 log(1 + exp(a0 + a1*E)) == E * exp(a0 + a1*E) / (1 + exp(a0 + a1*E))
        #                               == E / (1 + exp(-a0 - a1*E))
        def inv_1_plus_exp(y):
            def f(yy):
                with np.errstate(over='raise'):
                    try:
                        xx = 1. / (1 + np.exp(yy))
                    except FloatingPointError:
                        xx = 0.
                    return xx
            return np.vectorize(f)(y)

        if a.size == 2:
            d = np.zeros(2)
            d[0] = (np.sum(inv_1_plus_exp(-a[0] - a[1]*ene1)) -
                    np.sum(inv_1_plus_exp(a[0] + a[1]*ene2)))
            d[1] = (np.sum(inv_1_plus_exp(-a[0] - a[1]*ene1) * ene1) -
                    np.sum(inv_1_plus_exp(a[0] + a[1]*ene2) * ene2))
        else:
            d = np.zeros(3)
            d[0] = (np.sum(inv_1_plus_exp(-a[0] - a[1]*ene1[0] - a[2]*ene1[1])) -
                    np.sum(inv_1_plus_exp(a[0] + a[1]*ene2[0] + a[2]*ene2[1])))
            d[1] = (np.sum(inv_1_plus_exp(-a[0] - a[1]*ene1[0] - a[2]*ene1[1]) * ene1[0]) -
                    np.sum(inv_1_plus_exp(a[0] + a[1]*ene2[0] + a[2]*ene2[1]) * ene2[0]))
            d[2] = (np.sum(inv_1_plus_exp(-a[0] - a[1]*ene1[0] - a[2]*ene1[1]) * ene1[1]) -
                    np.sum(inv_1_plus_exp(a[0] + a[1]*ene2[0] + a[2]*ene2[1]) * ene2[1]))

        return d

    def hess_log_likelihood(a, ene1, ene2):
        # Returns the hessian wrt the parameters a of log_likelihood
        # fac1 = 1 / (2 + 2*cosh(a0 + a1*ene1))
        # h1 = [[ fac1,      ene1*fac1    ],
        #       [ ene1*fac1, ene1**2*fac1 ]]
        # fac2 = 1 / (2 + 2*cosh(a0 + a1*ene2))
        # h2 = [[ fac2,      ene2*fac2    ],
        #       [ ene2*fac2, ene2**2*fac2 ]]
        # h = h1 + h2

        if a.size == 2:
            fac1 = 1 / (2 + 2*np.cosh(a[0] + a[1]*ene1))
            fac2 = 1 / (2 + 2*np.cosh(a[0] + a[1]*ene2))

            h = np.zeros((2, 2))

            h[0, 0] = np.sum(fac1) + np.sum(fac2)
            h[0, 1] = h[1, 0] = np.sum(ene1 * fac1) + np.sum(ene2 * fac2)
            h[1, 1] = np.sum(ene1 * ene1 * fac1) + np.sum(ene2 * ene2 * fac2)

        else:
            fac1 = 1 / (2 + 2*np.cosh(a[0] + a[1]*ene1[0] + a[2]*ene1[1]))
            fac2 = 1 / (2 + 2*np.cosh(a[0] + a[1]*ene2[0] + a[2]*ene2[1]))

            h = np.zeros((3, 3))

            h[0, 0] = np.sum(fac1) + np.sum(fac2)
            h[1, 1] = np.sum(ene1[0] * ene1[0] * fac1) + np.sum(ene2[0] * ene2[0] * fac2)
            h[2, 2] = np.sum(ene1[1] * ene1[1] * fac1) + np.sum(ene2[1] * ene2[1] * fac2)

            h[0, 1] = h[1, 0] = np.sum(ene1[0] * fac1) + np.sum(ene2[0] * fac2)
            h[0, 2] = h[2, 0] = np.sum(ene1[1] * fac1) + np.sum(ene2[1] * fac2)
            h[1, 2] = h[2, 1] = (np.sum(ene1[0] * ene1[1] * fac1) +
                                 np.sum(ene2[0] * ene2[1] * fac2))

        return h

    # ==================================================== #
    # Minimize the negative of the log likelihood function #
    # ==================================================== #
    if init_params is None:
        init_params = np.zeros(traj1.ndim + 1)
    else:
        init_params = np.array(init_params)

    min_res = scipy.optimize.minimize(
        log_likelihood,
        x0=init_params,
        args=(traj1, traj2),
        method='dogleg',
        jac=da_log_likelihood,
        hess=hess_log_likelihood
    )

    # fallback options
    if not min_res.success:
        if verbose:
            print('Note: Max-Likelihood minimization failed using \'dogleg\' method. '
                  'Trying to vary initial parameters.')
        min_res_1 = scipy.optimize.minimize(
            log_likelihood,
            x0=init_params * 0.9,
            args=(traj1, traj2),
            method='dogleg',
            jac=da_log_likelihood,
            hess=hess_log_likelihood
        )
        min_res_2 = scipy.optimize.minimize(
            log_likelihood,
            x0=init_params * 1.1,
            args=(traj1, traj2),
            method='dogleg',
            jac=da_log_likelihood,
            hess=hess_log_likelihood
        )
        if min_res_1.success and min_res_2.success and np.allclose(min_res_1.x, min_res_2.x):
            min_res = min_res_1

    if not min_res.success:
        # dogleg was unsuccessful using alternative starting point
        if verbose:
            print('Note: Max-Likelihood minimization failed using \'dogleg\' method. '
                  'Trying method \'nelder-mead\'.')
        min_res = scipy.optimize.minimize(
            log_likelihood,
            x0=init_params * 0.9,
            args=(traj1, traj2),
            method='nelder-mead'
        )

    if not min_res.success:
        raise RuntimeError('MaxLikelihood: Unable to minimize function.')

    final_params = min_res.x

    # ======================= #
    # Calculate uncertainties #
    # ======================= #
    cov = np.linalg.inv(hess_log_likelihood(final_params, traj1, traj2))
    final_error = np.sqrt(np.diag(cov))*np.sqrt(np.average([g1, g2]))

    return final_params, final_error


def check_bins(traj1, traj2, bins):
    # check for empty bins
    h1, _ = np.histogram(traj1, bins=bins)
    h2, _ = np.histogram(traj2, bins=bins)
    empty = np.where((h1 == 0) | (h2 == 0))[0]

    if np.size(empty) == 0:
        return bins
    elif np.size(empty) == 1:
        empty = empty[0]
        if empty > np.size(bins) / 2:
            return bins[:empty]
        else:
            return bins[empty+1:]
    else:
        # find longest non-empty interval
        empty = np.insert(np.append(empty, [40]), 0, [-1])
        max_interval = np.argmax(empty[1:] - empty[:-1])
        left = empty[max_interval] + 1
        right = empty[max_interval + 1]
        return bins[left:right]


def print_stats(title,
                fitvals, dfitvals,
                kb, param1, param2, trueslope,
                temp=None, pvconvert=None,
                dtemp=False, dpress=False, dmu=False,
                dtempdpress=False, dtempdmu=False):

    # if simple 1d:
    #     fitvals = [df, slope]
    #     dfitvals = [ddf, dslope]
    # if simple 2d:
    #     fitvals = [df, slope0, slope1]
    #     dfitvals = [ddf, dslope0, dslope1]
    # if bootstrapped 1d:
    #     fitvals = [[df, slope], [df, slope], ...]
    #     dfitvals = None
    # if bootstrapped 2d:
    #     fitvals = [[df, slope0, slope1], [df, slope0, slope1], ...]
    #     dfitvals = None
    if fitvals.ndim > 1:
        dfitvals = np.std(fitvals, axis=0)
        fitvals = np.average(fitvals, axis=0)

    if np.ndim(trueslope) == 0:
        trueslopes = np.array([trueslope])
    else:
        trueslopes = trueslope

    free_energy = fitvals[0]
    slopes = fitvals[1:]
    dfree_energy = dfitvals[0]
    dslopes = dfitvals[1:]

    print('='*50)
    print(title)
    print('='*50)
    print('Free energy')
    print('    {:.5f} +/- {:.5f}'.format(free_energy, dfree_energy))
    print('{:27s}      |  {:s}'.format('Estimated slope', 'True slope'))
    for slope, dslope, trueslope in zip(slopes, dslopes, trueslopes):
        print('    {:<9.6f} +/- {:<9.6f}      |  {:<9.6f}'.format(slope, dslope, trueslope))
        quant = np.abs((slope-trueslope)/dslope)
        print('    ({:.2f} quantiles from true slope)'.format(quant))

    if dtemp or dtempdpress or dtempdmu:
        # slope is estimated beta2 - beta1
        # kb * slope == 1/T1' - 1/T2' == (T2' - T1')/(T1'*T2')
        # So we'll assume dT' == T1' - T2' ~= kb * slope * T1*T2
        slope = slopes[0]
        dslope = dslopes[0]
        if dtemp:
            t1 = param1
            t2 = param2
        else:
            t1 = param1[0]
            t2 = param2[0]
        print('{:27s}      |  {:s}'.format('Estimated dT', 'True dT'))
        print('    {:<6.1f} +/- {:<6.1f}            |  {:<6.1f}'.format(
            kb * slope * t1 * t2,
            kb * dslope * t1 * t2,
            t2 - t1
        ))
    if dpress or dtempdpress:
        # slope is estimated (P1 - P2)/beta*pvconvert (1d), or
        #                    (P1/b1 - P2/b2)*pvconvert (2d)
        if temp is None and dtempdpress:
            temp = .5*(param1[0] + param2[0])
        if dpress:
            press = -slopes[0] * (kb*temp) / pvconvert
            ddpress = -dslopes[0] * (kb*temp) / pvconvert
            truepress = -trueslopes[0] * (kb*temp) / pvconvert
        else:
            press = -slopes[1] * (kb*temp) / pvconvert
            ddpress = -dslopes[1] * (kb*temp) / pvconvert
            truepress = -trueslopes[1] * (kb*temp) / pvconvert
        print('{:27s}      |  {:s}'.format('Estimated dP', 'True dP'))
        print('    {:<6.1f} +/- {:<6.1f}            |  {:<6.1f}'.format(
            press, np.abs(ddpress), truepress
        ))
    if dmu or dtempdmu:
        pass
    print('='*50)


def estimate_interval(ens_string, ens_temp,
                      energy, kb,
                      ens_press=None,
                      volume=None, pvconvert=None,
                      verbosity=1, cutoff=0.001,
                      tunit='', punit=''):
    result = {}
    if ens_string == 'NVT':
        # Discard burn-in period and time-correlated frames
        energy = trajectory.equilibrate(energy, verbose=(verbosity > 1), name='Energy')
        energy = trajectory.decorrelate(energy, verbose=(verbosity > 1), name='Energy')
        energy = trajectory.cut_tails(energy, cut=cutoff, verbose=(verbosity > 2), name='Energy')

        # dT
        sig = np.std(energy)
        result['dT'] = 2*kb*ens_temp*ens_temp/sig
    elif ens_string == 'NPT':
        enthalpy = energy + pvconvert * ens_press * volume
        traj_2d = np.array([energy, volume])
        # Discard burn-in period and time-correlated frames
        enthalpy = trajectory.equilibrate(enthalpy, verbose=(verbosity > 1), name='Enthalpy')
        enthalpy = trajectory.decorrelate(enthalpy, verbose=(verbosity > 1), name='Enthalpy')
        enthalpy = trajectory.cut_tails(enthalpy, cut=cutoff, verbose=(verbosity > 2), name='Enthalpy')
        volume_1d = trajectory.equilibrate(volume, verbose=(verbosity > 1), name='Volume')
        volume_1d = trajectory.decorrelate(volume_1d, verbose=(verbosity > 1), name='Volume')
        volume_1d = trajectory.cut_tails(volume_1d, cut=cutoff, verbose=(verbosity > 2), name='Volume')
        traj_2d = trajectory.equilibrate(traj_2d, verbose=(verbosity > 1), name='2D-Trajectory')
        traj_2d = trajectory.decorrelate(traj_2d, facs=[1, pvconvert * ens_press], verbose=(verbosity > 1), name='2D-Trajectory')
        traj_2d = trajectory.cut_tails(traj_2d, cut=cutoff, verbose=(verbosity > 2), name='2D-Trajectory')

        # dT
        sig = np.std(enthalpy)
        result['dT'] = 2*kb*ens_temp*ens_temp/sig
        # dP
        sig = np.std(volume_1d)*pvconvert
        result['dP'] = 2*kb*ens_temp/sig
        # dTdP
        cov = np.cov(traj_2d)
        sig = np.sqrt(np.diag(cov))
        sig[1] *= pvconvert
        result['dTdP'] = [2*kb*ens_temp*ens_temp/sig[0],
                          2*kb*ens_temp/sig[1]]
    else:
        raise pv_error.InputError('ens_str', 'Unrecognized ensemble string.')

    if verbosity > 0:
        print('A rule of thumb states that good error recognition can be expected when\n'
              'spacing the tip of the distributions by about two standard deviations.\n'
              'Based on this rule, and the assumption that the standard deviation of the\n'
              'distributions is largely independent of the state point, here\'s an estimate\n'
              'for the interval given the current simulation:')
        if ens_string == 'NVT':
            print('Current trajectory: NVT, T = {:.2f} {:s}'.format(ens_temp, tunit))
            print('Suggested interval: dT = {:.1f} {:s}'.format(result['dT'], tunit))
        if ens_string == 'NPT':
            print('Current trajectory: NPT, T = {:.2f} {:s}, P = {:.2f} {:s}'.format(
                ens_temp, tunit, ens_press, punit))
            print('Suggested interval:')
            print('  Temperature-only: dT = {:.1f} {:s}'.format(result['dT'], tunit))
            print('  Pressure-only: dP = {:.1f} {:s}'.format(result['dP'], punit))
            print('  Combined: dT = {:.1f} {:s}, dP = {:.1f} {:s}'.format(
                result['dTdP'][0], tunit, result['dTdP'][1], punit))


def check_1d(traj1, traj2, param1, param2, kb,
             quantity, dtemp=False, dpress=False, dmu=False,
             temp=None, pvconvert=None,
             nbins=40, cutoff=0.001, seed=None,
             verbosity=1, screen=False, filename=None):
    r"""
    Checks whether the energy trajectories of two simulation performed at
    different temperatures have sampled distributions at the analytically
    expected ratio.

    Parameters
    ----------
    traj1 : array-like
        Trajectory of the first simulation
        If dtemp:

            * NVT: Potential energy U or total energy E = U + K
            * NPT: Enthalpy H = U + pV or total energy E = H + K

        If dpress:

            * NPT: Volume V

    traj2 : array-like
        Trajectory of the second simulation
        If dtemp:

            * NVT: Potential energy U or total energy E = U + K
            * NPT: Enthalpy H = U + pV or total energy E = H + K

        If dpress:

            * NPT: Volume V

    param1 : float
        Target temperature or pressure of the first simulation
    param2 : float
        Target temperature or pressure of the second simulation
    kb : float
        Boltzmann constant in same units as the energy trajectories
    quantity : str
        Name of quantity analyzed (used for printing only)
    dtemp : bool, optional
        Set to True if trajectories were simulated at different temperature
        Default: False.
    dpress : bool, optional
        Set to True if trajectories were simulated at different pressure
        Default: False.
    temp : float, optional
        The temperature in equal temperature, differring pressure NPT simulations.
        Needed to print optimal dP.
    pvconvert : float, optional
        Conversion from pressure * volume to energy units.
        Needed to print optimal dP.
    dmu : bool, optional
        Set to True if trajectories were simulated at different chemical potential
        Default: False.
    nbins : int, optional
        Number of bins used to assess distributions of the trajectories
        Default: 40
    cutoff : float, optional
        Tail cutoff of distributions.
        Default: 0.001 (0.1%)
    seed : int, optional
        If set, bootstrapping will be reproducible.
        Default: None, bootstrapping non-reproducible.
    verbosity : int, optional
        Verbosity level.
        Default: 1 (only most important output)
    screen : bool, optional
        Plot distributions on screen.
        Default: False.
    filename : string, optional
        Plot distributions to `filename`.pdf.
        Default: None.

    Returns
    -------

    """

    if (not (dtemp or dpress or dmu) or
       (dtemp and dpress) or
       (dtemp and dmu) or
       (dpress and dmu)):
        raise pv_error.InputError(['dtemp', 'dpress', 'dmu'],
                                  'Need to specify exactly one of `dtemp`, `dpress` and `dmu`.')

    if dmu:
        raise NotImplementedError('check_1d: Testing of `dmu` not implemented.')

    if seed is not None:
        raise NotImplementedError('check_1d: Bootstrapping not implemented.')

    if dpress and (temp is None or pvconvert is None):
        raise pv_error.InputError(['dpress', 'temp', 'pvconvert'],
                                  '`ensemble.check_1d` with `dpress=True` requires `temp` and `pvconvert`.')

    # =============================== #
    # prepare constants, strings etc. #
    # =============================== #
    pstring = 'ln(P_2(' + quantity + ')/P_1(' + quantity + '))'
    trueslope = 0
    if dtemp:
        trueslope = 1/(kb * param1) - 1/(kb * param2)
    elif dpress:
        trueslope = (param1 - param2) / (kb * temp) * pvconvert

    if verbosity > 1:
        print('Analytical slope of {:s}: {:.8f}'.format(pstring, trueslope))

    quant = {}

    # ==================== #
    # prepare trajectories #
    # ==================== #
    # Discard burn-in period and time-correlated frames
    traj1 = trajectory.prepare(traj1, cut=cutoff, verbosity=verbosity, name='Trajectory 1')
    traj2 = trajectory.prepare(traj2, cut=cutoff, verbosity=verbosity, name='Trajectory 2')

    # calculate inefficiency
    g1 = pymbar.timeseries.statisticalInefficiency(traj1)
    g2 = pymbar.timeseries.statisticalInefficiency(traj2)

    # calculate overlap
    traj1_full = traj1
    traj2_full = traj2
    traj1, traj2, min_ene, max_ene = trajectory.overlap(
        traj1=traj1_full, traj2=traj2_full,
    )
    if verbosity > 0:
        print('Overlap is {:.1%} of trajectory 1 and {:.1%} of trajectory 2.'.format(
            traj1.shape[0] / traj1_full.shape[0],
            traj2.shape[0] / traj2_full.shape[0]
        ))
    if verbosity > 0 and dtemp:
        sig1 = np.std(traj1_full)
        sig2 = np.std(traj2_full)
        dt1 = 2*kb*param1*param1/sig1
        dt2 = 2*kb*param2*param2/sig2
        if verbosity > 1:
            print('A rule of thumb states that a good overlap is found when dT/T = (2*kB*T)/(sig),\n'
                  'where sig is the standard deviation of the energy distribution.\n'
                  'For the current trajectories, dT = {:.1f}, sig1 = {:.1f} and sig2 = {:.1f}.\n'
                  'According to the rule of thumb, given T1, a good dT is dT = {:.1f}, and\n'
                  '                                given T2, a good dT is dT = {:.1f}.'.format(
                      param2-param1, sig1, sig2, dt1, dt2)
                  )
        print('Rule of thumb estimates that dT = {:.1f} would be optimal '
              '(currently, dT = {:.1f})'.format(.5*(dt1+dt2), param2-param1))
    if verbosity > 0 and dpress:
        sig1 = np.std(traj1_full)*pvconvert
        sig2 = np.std(traj2_full)*pvconvert
        dp1 = 2*kb*temp/sig1
        dp2 = 2*kb*temp/sig2
        if verbosity > 1:
            print('A rule of thumb states that a good overlap is found when dP = (2*kB*T)/(sig),\n'
                  'where sig is the standard deviation of the volume distribution.\n'
                  'For the current trajectories, dP = {:.1f}, sig1 = {:.1g} and sig2 = {:.1g}.\n'
                  'According to the rule of thumb, given P1, a good dP is dP = {:.1f}, and\n'
                  '                                given P2, a good dP is dP = {:.1f}.'.format(
                      param2-param1, sig1, sig2, dp1, dp2)
                  )
        print('Rule of thumb estimates that dP = {:.1f} would be optimal '
              '(currently, dP = {:.1f})'.format(.5*(dp1+dp2), param2-param1))
    if not min_ene:
        raise pv_error.InputError(['traj1', 'traj2'],
                                  'No overlap between trajectories.')
    # calculate bins
    bins = np.linspace(min_ene, max_ene, nbins+1)
    bins = check_bins(traj1, traj2, bins)
    if np.size(bins) < 3:
        raise pv_error.InputError(['traj1', 'traj2', 'nbins', 'cutoff'],
                                  'Less than 3 bins were filled in the overlap region.\n'
                                  'Ensure sufficient overlap between the trajectories, and '
                                  'consider increasing `cutoff` or `nbins` if there is '
                                  'sufficient overlap but unusually long tails.')

    w_f = -trueslope * traj1_full
    w_r = trueslope * traj2_full

    if verbosity > 2:
        print('Computing log of partition functions using pymbar.BAR...')
    df, ddf = pymbar.BAR(w_f, w_r)
    if verbosity > 2:
        print('Using {:.5f} for log of partition functions as computed from BAR.'.format(df))
        print('Uncertainty in quantity is {:.5f}.'.format(ddf))
        print('Assuming this is negligible compared to sampling error at individual points.')

    # ========== #
    # linear fit #
    # ========== #
    if verbosity > 2:
        print('Computing linear fit parameters (for plotting / comparison)')

    fitvals, dfitvals = do_linear_fit(
        traj1=traj1_full, traj2=traj2_full, g1=g1, g2=g2, bins=bins,
        screen=screen, filename=filename,
        trueslope=trueslope, trueoffset=df,
        units=None
    )

    slope = fitvals[1]
    dslope = dfitvals[1]
    quant['linear'] = [abs((slope - trueslope)/dslope)]
    if verbosity > 1:
        print_stats(
            title='Linear Fit Analysis (analytical error)',
            fitvals=fitvals,
            dfitvals=dfitvals,
            kb=kb,
            param1=param1,
            param2=param2,
            trueslope=trueslope,
            temp=temp, pvconvert=pvconvert,
            dtemp=dtemp, dpress=dpress, dmu=dmu
        )

    # ================== #
    # max-likelihood fit #
    # ================== #
    if verbosity > 2:
        print('Computing the maximum likelihood parameters')

    fitvals, dfitvals = do_max_likelihood_fit(traj1_full, traj2_full, g1, g2,
                                              init_params=[df, trueslope],
                                              verbose=(verbosity > 1))

    slope = fitvals[1]
    dslope = dfitvals[1]
    quant['maxLikelihood'] = [abs((slope - trueslope)/dslope)]
    if verbosity > 0:
        print_stats(
            title='Maximum Likelihood Analysis (analytical error)',
            fitvals=fitvals,
            dfitvals=dfitvals,
            kb=kb,
            param1=param1,
            param2=param2,
            trueslope=trueslope,
            temp=temp, pvconvert=pvconvert,
            dtemp=dtemp, dpress=dpress, dmu=dmu
        )

    return quant['maxLikelihood']


def check_2d(traj1, traj2, param1, param2, kb, pvconvert,
             quantity, dtempdpress=False, dtempdmu=False,
             cutoff=0.001, seed=None,
             verbosity=1, screen=False, filename=None):
    r"""
    Checks whether the energy trajectories of two simulation performed at
    different temperatures have sampled distributions at the analytically
    expected ratio.

    Parameters
    ----------
    traj1 : array-like, 2d
        Trajectory of the first simulation
        If dtempdpress:

            * traj[0,:]: Potential energy U or total energy E = U + K
            * traj[1,:]: Volume V
    traj2 : array-like, 2d
        Trajectory of the second simulation
        If dtempdpress:

            * traj[0,:]: Potential energy U or total energy E = U + K
            * traj[1,:]: Volume V
    param1 : array-like
        If dtempdpress:
            Target temperature and pressure of the first simulation
    param2 : array-like
        If dtempdpress:
            Target temperature and pressure of the first simulation
    kb : float
        Boltzmann constant in same units as the energy trajectories
    pvconvert : float
        Conversion from pressure * volume to energy units
    quantity : List[str]
        Names of quantities analyzed (used for printing only)
    dtempdpress : bool, optional
        Set to True if trajectories were simulated at different
        temperature and pressure
        Default: False.
    dtempdmu : bool, optional
        Set to True if trajectories were simulated at different
        temperature and chemical potential
        Default: False.
    cutoff : float
        Tail cutoff of distributions.
        Default: 0.001 (0.1%)
    seed : int
        If set, bootstrapping will be reproducible.
        Default: None, bootstrapping non-reproducible.
    verbosity : int
        Verbosity level.
        Default: 1 (only most important output)
    screen : bool, optional
        Plot distributions on screen.
        Default: False.
    filename : string, optional
        Plot distributions to `filename`.pdf.
        Default: None.

    Returns
    -------

    """

    if not (dtempdpress or dtempdmu) or (dtempdpress and dtempdmu):
        raise pv_error.InputError(['dtempdpress', 'dtempdmu'],
                                  'Need to specify exactly one of `dtempdpress` and `dtempdmu`.')

    if dtempdmu:
        raise NotImplementedError('check_2d: Testing of `dtempdmu` not implemented.')

    if seed is not None:
        raise NotImplementedError('check_2d: Bootstrapping not implemented.')

    if screen or filename is not None:
        raise NotImplementedError('check_2d: Plotting not implemented.')

    # =============================== #
    # prepare constants, strings etc. #
    # =============================== #
    pstring = ('ln(P_2(' + quantity[0] + ', ' + quantity[1] + ')/' +
               'P_1(' + quantity[0] + ', ' + quantity[1] + '))')
    trueslope = np.zeros(2)
    facs = [None, None]
    if dtempdpress:
        trueslope = np.array([
            1/(kb * param1[0]) - 1/(kb * param2[0]),
            pvconvert*(1/(kb * param1[0]) * param1[1] - 1/(kb * param2[0]) * param2[1])
        ])
        facs = [[1, param1[1]], [1, param2[1]]]

    if verbosity > 1:
        print('Analytical slope of {:s}: {:.8f}, {:.8f}'.format(
            pstring, trueslope[0], trueslope[1]
        ))

    quant = {}

    # ==================== #
    # prepare trajectories #
    # ==================== #
    # Discard burn-in period and time-correlated frames
    traj1 = trajectory.prepare(traj1, cut=cutoff, facs=facs[0],
                               verbosity=verbosity, name='Trajectory 1')
    traj2 = trajectory.prepare(traj2, cut=cutoff, facs=facs[1],
                               verbosity=verbosity, name='Trajectory 2')

    # calculate inefficiency
    g1 = np.array([
            pymbar.timeseries.statisticalInefficiency(traj1[0]),
            pymbar.timeseries.statisticalInefficiency(traj1[1])
        ])
    g2 = np.array([
            pymbar.timeseries.statisticalInefficiency(traj2[0]),
            pymbar.timeseries.statisticalInefficiency(traj2[1])
        ])

    # calculate overlap
    traj1_full = traj1
    traj2_full = traj2
    traj1, traj2, min_ene, max_ene = trajectory.overlap(
        traj1=traj1_full, traj2=traj2_full,
    )
    if verbosity > 0:
        print('Overlap is {:.1%} of trajectory 1 and {:.1%} of trajectory 2.'.format(
            traj1.shape[1] / traj1_full.shape[1],
            traj2.shape[1] / traj2_full.shape[1]
        ))
    if verbosity > 0 and dtempdpress:
        cov1 = np.cov(traj1_full)
        sig1 = np.sqrt(np.diag(cov1))
        sig1[1] *= pvconvert
        cov2 = np.cov(traj2_full)
        sig2 = np.sqrt(np.diag(cov2))
        sig2[1] *= pvconvert
        dt1 = 2*kb*param1[0]*param1[0]/sig1[0]
        dt2 = 2*kb*param2[0]*param2[0]/sig2[0]
        dp1 = 2*kb*param1[0]/sig1[1]
        dp2 = 2*kb*param2[0]/sig2[1]
        if verbosity > 1:
            print('A rule of thumb states that a good overlap can be expected when choosing state\n'
                  'points separated by about 2 standard deviations.\n'
                  'For the current trajectories, dT = {:.1f}, and dP = {:.1f},\n'
                  'with standard deviations sig1 = [{:.1f}, {:.1g}], and sig2 = [{:.1f}, {:.1g}].\n'
                  'According to the rule of thumb, given point 1, the estimate is dT = {:.1f}, dP = {:.1f}, and\n'
                  '                                given point 2, the estimate is dT = {:.1f}, dP = {:.1f}.'.format(
                      param2[0]-param1[0], param2[1]-param1[1],
                      sig1[0], sig1[1], sig2[0], sig2[1],
                      dt1, dt2, dp1, dp2)
                  )
        print('Rule of thumb estimates that (dT,dP) = ({:.1f},{:.1f}) would be optimal '
              '(currently, (dT,dP) = ({:.1f},{:.1f}))'.format(.5*(dt1+dt2), .5*(dp1+dp2),
                                                              param2[0]-param1[0], param2[1]-param1[1]))
    if min_ene is None:
        raise pv_error.InputError(['traj1', 'traj2'],
                                  'No overlap between trajectories.')

    w_f = -trueslope[0] * traj1_full[0] - trueslope[1] * traj1_full[1]
    w_r = trueslope[0] * traj2_full[0] + trueslope[1] * traj2_full[1]

    if verbosity > 2:
        print('Computing log of partition functions using pymbar.BAR...')
    df, ddf = pymbar.BAR(w_f, w_r)
    if verbosity > 2:
        print('Using {:.5f} for log of partition functions as computed from BAR.'.format(df))
        print('Uncertainty in quantity is {:.5f}.'.format(ddf))
        print('Assuming this is negligible compared to sampling error at individual points.')

    # ================== #
    # max-likelihood fit #
    # ================== #
    if verbosity > 2:
        print('Computing the maximum likelihood parameters')

    fitvals, dfitvals = do_max_likelihood_fit(traj1_full, traj2_full, g1, g2,
                                              init_params=[df, trueslope[0], trueslope[1]],
                                              verbose=(verbosity > 1))

    slope = fitvals[1:]
    dslope = dfitvals[1:]
    quant['maxLikelihood'] = np.abs((slope - trueslope)/dslope)
    if verbosity > 0:
        print_stats(
            title='Maximum Likelihood Analysis (analytical error)',
            fitvals=fitvals,
            dfitvals=dfitvals,
            kb=kb,
            param1=param1,
            param2=param2,
            trueslope=trueslope,
            pvconvert=pvconvert,
            dtempdpress=dtempdpress, dtempdmu=dtempdmu
        )

    return quant['maxLikelihood']
