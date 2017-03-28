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
This module contains low-level functionality of the
`physical_validation.kinetic_energy` module. The functions in this module
should generally not be called directly. Please use the high-level
functions from `physical_validation.kinetic energy`.
"""
from __future__ import print_function
from __future__ import division

import scipy.stats as stats
import numpy as np

from . import plot


def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def temperature(kin, ndof, kb=8.314e-3):
    r"""
    Calculates the temperature acccording to the equipartition theorem.

    .. math::
        T(K) = \frac{2K}{N k_B}

    Parameters
    ----------
    kin : float
        Kinetic energy.
    ndof : float
        Number of degrees of freedom.
    kb : float
        Boltzmann constant :math:`k_B`. Default: 8.314e-3 (kJ/mol).

    Returns
    -------
    temperature : float
        Calculated temperature.
    """
    # ndof * kb * T = 2 * kin
    if isclose(ndof, 0):
        return 0
    return 2 * float(kin) / (float(ndof) * float(kb))


def check_mb_ensemble(kin, temp, ndof, alpha, kb=8.314e-3, verbose=False,
                      screen=False, filename=None, ene_unit=None):
    r"""
    Checks if a kinetic energy trajectory is Maxwell-Boltzmann distributed.

    .. warning: This is a low-level function. Additionally to being less
       user-friendly, there is a higher probability of erroneous and / or
       badly documented behavior due to unexpected inputs. Consider using
       the high-level version based on the SimulationData object. See
       physical_validation.kinetic_energy.check_mb_ensemble for more
       information and full documentation.

    Parameters
    ----------
    kin : array-like
        Kinetic energy snapshots of the system.
    temp : float
        Target temperature of the system. Used to construct the
        Maxwell-Boltzmann distribution.
    ndof : float
        Number of degrees of freedom in the system. Used to construct the
        Maxwell-Boltzmann distribution.
    alpha : float
        Confidence. TODO: Check proper statistical definition.
    kb : float
        Boltzmann constant :math:`k_B`. Default: 8.314e-3 (kJ/mol).
    verbose : bool
        Print result details. Default: False.
    screen : bool
        Plot distributions on screen. Default: False.
    filename : string
        Plot distributions to `filename`.pdf. Default: None.
    ene_unit : string
        Energy unit - used for plotting only.

    Returns
    -------
    result : float
        The p value of the test.

    See Also
    --------
    physical_validation.kinetic_energy.check_mb_ensemble : High-level version
    """

    kt = kb * temp
    d, p = stats.kstest(kin, 'chi2', (ndof, 0, kt/2))

    do_plot = screen or filename is not None

    if do_plot:
        hist_sim, k_sim = np.histogram(kin, bins=25, normed=True)
        k_sim = (k_sim[1:] + k_sim[:-1])/2

        ana_dist = stats.chi2(df=ndof, scale=kt/2)
        k_ana = np.linspace(ana_dist.ppf(0.0001),
                            ana_dist.ppf(0.9999), 100)
        hist_ana = ana_dist.pdf(k_ana)

        data = [{'x': k_sim,
                 'y': hist_sim * 100,
                 'name': 'Simulation result'},
                {'x': k_ana,
                 'y': hist_ana * 100,
                 'name': 'Maxwell-Boltzmann'}]

        unit = ''
        if ene_unit is not None:
            unit = ' [' + ene_unit + ']'

        plot.plot(data,
                  legend='best',
                  title='Simulation vs. Maxwell-Boltzmann',
                  xlabel='Kinetic energy' + unit,
                  ylabel='Probability [%]',
                  filename=filename,
                  screen=screen)

    if verbose:
        message = ('Kolmogorov-Smirnov test result: p = {:g}\n'
                   'Null hypothesis: Kinetic energy is Maxwell-Boltzmann distributed'.format(p))
        if alpha is not None:
            if p >= alpha:
                message += ('\nConfidence alpha = {:f}\n'
                            'Result: Hypothesis stands'.format(alpha))
            elif p < alpha:
                message += ('\nConfidence alpha = {:f}\n'
                            'Result: Hypothesis rejected'.format(alpha))
        print(message)

    return p


def check_equipartition(positions, velocities, masses,
                        molec_idx, molec_nbonds,
                        natoms, nmolecs,
                        ndof_reduction_tra=0, ndof_reduction_rot=0,
                        dtemp=0.1, temp=None, alpha=0.05,
                        molec_groups=None,
                        random_divisions=0, random_groups=2,
                        ndof_molec=None, kin_molec=None,
                        verbosity=2,
                        screen=False, filename=None):
    r"""
    Checks the equipartition of a simulation trajectory.

    .. warning: This is a low-level function. Additionally to being less
       user-friendly, there is a higher probability of erroneous and / or
       badly documented behavior due to unexpected inputs. Consider using
       the high-level version based on the SimulationData object. See
       physical_validation.kinetic_energy.check_equipartition for more
       information and full documentation.

    Parameters
    ----------
    positions : array-like (nframes x natoms x 3)
        3d array containing the positions of all atoms for all frames
    velocities : array-like (nframes x natoms x 3)
        3d array containing the velocities of all atoms for all frames
    masses : array-like (natoms x 1)
        1d array containing the masses of all atoms
    molec_idx : array-like (nmolecs x 1)
        Index of first atom for every molecule
    molec_nbonds : array-like (nmolecs x 1)
        Number of bonds for every molecule
    natoms : int
        Number of atoms in the system
    nmolecs : int
        Number of molecules in the system
    ndof_reduction_tra : int, optional
        Number of center-of-mass translational degrees of freedom to
        remove. Default: 0.
    ndof_reduction_rot : int, optional
        Number of center-of-mass rotational degrees of freedom to remove.
        Default: 0.
    dtemp : float, optional
        Fraction of temperature deviation tolerated between groups.
        Default: 0.05 (5%).
    temp : float, optional
        Target temperature of the simulation. If None, the kinetic
        energies will not be tested for Maxwell-Boltzmann distribution,
        but only compared amongst each others. Default: None.
    alpha : float, optional
        Confidence for Maxwell-Boltzmann test. Default: 0.05 (5%).
    molec_groups : List[array-like] (ngroups x ?), optional
        List of 1d arrays containing molecule indeces defining groups.
        Useful to pre-define groups of molecules (e.g. solute / solvent,
        liquid mixture species, ...). If None, no pre-defined molecule
        groups will be tested. Default: None.
        Note: If an empty 1d array is found as last element in the list, the remaining
              molecules are collected in this array. This allows, for example, to only
              specify the solute, and indicate the solvent by giving an empty array.
    random_divisions : int, optional
        Number of random division tests attempted. Default: 0 (random
        division tests off).
    random_groups : int, optional
        Number of groups the system is randomly divided in. Default: 2.
    ndof_molec : List[dict], optional
        Pass in the degrees of freedom per molecule. Slightly increases speed of repeated
        analysis of the same simulation run.
    kin_molec : List[List[dict]], optional
        Pass in the kinetic energy per molecule. Greatly increases speed of repeated
        analysis of the same simulation run.
    verbosity : int, optional
        Verbosity level, where 0 is quiet and 3 very chatty. Default: 2.
    screen : bool
        Plot distributions on screen. Default: False.
    filename : string
        Plot distributions to `filename`.pdf. Default: None.

    Returns
    -------
    result : int
        Number of equipartition violations. Tune up verbosity for details.
    ndof_molec : List[dict]
        List of the degrees of freedom per molecule. Can be saved to increase speed of
        repeated analysis of the same simulation run.
    kin_molec : List[List[dict]]
        List of the kinetic energy per molecule per frame. Can be saved to increase speed
        of repeated analysis of the same simulation run.

    See Also
    --------
    physical_validation.kinetic_energy.check_equipartition : High-level version

    """

    dict_keys = ['tot', 'tra', 'rni', 'rot', 'int']

    # for each molecule, calculate total / translational / rotational & internal /
    #   rotational / internal degrees of freedom
    #   returns: List[dict] of floats (shape: nmolecs x 5 x 1)
    if ndof_molec is None:
        ndof_molec = calc_ndof(natoms, nmolecs, molec_idx, molec_nbonds,
                               ndof_reduction_tra, ndof_reduction_rot)

    # for each frame, calculate total / translational / rotational & internal /
    #   rotational / internal kinetic energy for each molecule
    if kin_molec is None:
        kin_molec = []
        for r, v in zip(positions, velocities):
            kin_molec.append(calc_molec_kinetic_energy(r, v, masses,
                                                       molec_idx, natoms, nmolecs))

    result = []

    # test system-wide tot, tra, rni, rot, int
    if temp is not None:
        # check for Maxwell-Boltzmann distribution of
        # partitioned kinetic energy trajectories
        result.extend(test_mb_dist(kin_molec=kin_molec,
                                   ndof_molec=ndof_molec,
                                   nmolecs=nmolecs,
                                   temp=temp,
                                   alpha=alpha,
                                   dict_keys=dict_keys,
                                   verbosity=verbosity,
                                   screen=screen,
                                   filename=filename))
    else:
        # compare partitioned temperatures to total temperature
        result.extend(test_temp_diff(kin_molec=kin_molec,
                                     ndof_molec=ndof_molec,
                                     nmolecs=nmolecs,
                                     dtemp=dtemp,
                                     dict_keys=dict_keys,
                                     verbosity=verbosity,
                                     screen=screen,
                                     filename=filename))
    # divide in random groups
    for i in range(random_divisions):
        # randomly assign a group index to each molecule
        group_idx = np.random.randint(random_groups, size=nmolecs)
        # create molecule index for each group
        groups = []
        for rg in range(random_groups):
            groups.append(np.arange(nmolecs)[group_idx == rg])
        # test each group separately
        for rg, group in enumerate(groups):
            if verbosity > 0:
                print('Testing randomly divided group {:d}'.format(rg))
            if verbosity > 3:
                print(group)
            if temp is not None:
                result.extend(test_mb_dist(kin_molec, ndof_molec, nmolecs, temp,
                                           alpha, dict_keys, group, verbosity))
            else:
                result.extend(test_temp_diff(kin_molec, ndof_molec, nmolecs,
                                             dtemp, dict_keys, group, verbosity))
        # test groups against each others
        for rg1, group1 in enumerate(groups):
            for group2 in groups[rg1 + 1:]:
                result.extend(test_temp_diff_groups(kin_molec, ndof_molec, nmolecs,
                                                    group1, group2,
                                                    dtemp, dict_keys, verbosity))

    # use predefined group division?
    # if no groups, return
    if not molec_groups:
        return result, ndof_molec, kin_molec
    # is last group empty?
    last_empty = not molec_groups[-1]
    # get rid of empty groups
    molec_groups = [g for g in molec_groups if g]
    # if no groups now (all were empty), return now
    if not molec_groups:
        return result, ndof_molec, kin_molec

    if last_empty:
        # last group is [] -> insert remaining molecules
        combined = []
        for group in molec_groups:
            combined.extend(group)
        molec_groups[-1] = [m for m in range(nmolecs) if m not in combined]

    for mg, group in enumerate(molec_groups):
        if verbosity > 0:
            print('Testing predifined divided group {:d}'.format(mg))
        if verbosity > 3:
            print(group)
        if temp is not None:
            result.extend(test_mb_dist(kin_molec, ndof_molec, nmolecs, temp,
                                       alpha, dict_keys, group, verbosity))
        else:
            result.extend(test_temp_diff(kin_molec, ndof_molec, nmolecs,
                                         dtemp, dict_keys, group, verbosity))
    # test groups against each others
    if len(molec_groups) > 1:
        for rg1, group1 in enumerate(molec_groups):
            for group2 in molec_groups[rg1 + 1:]:
                result.extend(test_temp_diff_groups(kin_molec, ndof_molec, nmolecs,
                                                    group1, group2,
                                                    dtemp, dict_keys, verbosity))

    return result, ndof_molec, kin_molec


def calc_system_ndof(natoms, nmolecs, nbonds,
                     stop_com_tra, stop_com_rot):
    r"""
    Calculates the total / translational / rotational & internal /
    rotational / internal degrees of freedom of the system.

    Parameters
    ----------
    natoms : int
        Total number of atoms in the system
    nmolecs : int
        Total number of molecules in the system
    nbonds : int
        Total number of bonds in the system
    stop_com_tra : bool
        Was the center-of-mass translation removed during the simulation?
    stop_com_rot : bool
        Was the center-of-mass translation removed during the simulation?

    Returns
    -------
    ndof : dict
        Dictionary containing the degrees of freedom.
        Keys: ['tot', 'tra', 'rni', 'rot', 'int']
    """
    # total ndof
    ndof_tot = 3*natoms - nbonds

    # ndof reduction due to COM motion constraining
    if stop_com_tra:
        ndof_tot -= 3
    if stop_com_rot:
        ndof_tot -= 3

    # translational ndof
    ndof_tot_tra = 3*nmolecs
    if stop_com_tra:
        ndof_tot -= 3

    # rotational & internal ndof
    ndof_tot_rni = ndof_tot - ndof_tot_tra

    # rotational ndof
    ndof_tot_rot = 3*nmolecs
    if stop_com_tra:
        ndof_tot -= 3

    # internal ndof
    ndof_tot_int = ndof_tot_rni - ndof_tot_rot

    # return dict
    ndof = {'tot': ndof_tot,
            'tra': ndof_tot_tra,
            'rni': ndof_tot_rni,
            'rot': ndof_tot_rot,
            'int': ndof_tot_int}
    return ndof


def calc_ndof(natoms, nmolecs,
              molec_idx, molec_nbonds,
              ndof_reduction_tra, ndof_reduction_rot):
    r"""
    Calculates the total / translational / rotational & internal /
    rotational / internal degrees of freedom per molecule.

    Parameters
    ----------
    natoms : int
        Total number of atoms in the system
    nmolecs : int
        Total number of molecules in the system
    molec_idx : List[int]
        Index of first atom for every molecule
    molec_nbonds : List[int]
        Number of bonds for every molecule
    ndof_reduction_tra : int
        Number of center-of-mass translational degrees of freedom to
        remove. Default: 0.
    ndof_reduction_rot : int
        Number of center-of-mass rotational degrees of freedom to remove.
        Default: 0.

    Returns
    -------
    ndof_molec : List[dict]
        List of dictionaries containing the degrees of freedom for each molecule
        Keys: ['tot', 'tra', 'rni', 'rot', 'int']
    """
    # check whether there are monoatomic molecules:
    nmono = (np.array(molec_idx[1:]) - np.array(molec_idx[:-1]) == 1).sum()

    # ndof to be deducted per molecule
    # ndof reduction due to COM motion constraining
    ndof_com_tra_pm = ndof_reduction_tra / nmolecs
    ndof_com_rot_pm = ndof_reduction_rot / (nmolecs - nmono)

    ndof_molec = []
    # add last idx to molec_idx to ease looping
    molec_idx = np.append(molec_idx, [natoms])
    # loop over molecules
    for idx_molec, (idx_atm_init, idx_atm_end) in enumerate(zip(molec_idx[:-1], molec_idx[1:])):
        natoms = idx_atm_end - idx_atm_init
        nbonds = molec_nbonds[idx_molec]
        ndof_tot = 3*natoms - nbonds - ndof_com_tra_pm - ndof_com_rot_pm
        ndof_tra = 3 - ndof_com_tra_pm
        ndof_rni = ndof_tot - ndof_tra
        ndof_rot = 3 - ndof_com_rot_pm
        ndof_int = ndof_tot - ndof_tra - ndof_rot
        if isclose(ndof_int, 0, abs_tol=1e-09):
            ndof_int = 0
        if natoms == 1:
            ndof_tot = 3 - ndof_com_tra_pm
            ndof_tra = 3 - ndof_com_tra_pm
            ndof_rni = 0
            ndof_rot = 0
            ndof_int = 0
        ndof_molec.append({'tot': ndof_tot,
                           'tra': ndof_tra,
                           'rni': ndof_rni,
                           'rot': ndof_rot,
                           'int': ndof_int})

    return ndof_molec


def calc_molec_kinetic_energy(pos, vel, masses,
                              molec_idx, natoms, nmolecs):
    r"""
    Calculates the total / translational / rotational & internal /
    rotational / internal kinetic energy per molecule.

    Parameters
    ----------
    pos : nd-array (natoms x 3)
        2d array containing the positions of all atoms
    vel : nd-array (natoms x 3)
        2d array containing the velocities of all atoms
    masses : nd-array (natoms x 1)
        1d array containing the masses of all atoms
    molec_idx : nd-array (nmolecs x 1)
        Index of first atom for every molecule
    natoms : int
        Total number of atoms in the system
    nmolecs : int
        Total number of molecules in the system

    Returns
    -------
    kin : List[dict]
        List of dictionaries containing the kinetic energies for each molecule
        Keys: ['tot', 'tra', 'rni', 'rot', 'int']
    """
    # add last idx to molec_idx to ease looping
    molec_idx = np.append(molec_idx, [natoms])

    # calculate kinetic energy
    kin_tot = np.zeros(nmolecs)
    kin_tra = np.zeros(nmolecs)
    kin_rni = np.zeros(nmolecs)
    kin_rot = np.zeros(nmolecs)
    kin_int = np.zeros(nmolecs)
    # loop over molecules
    for idx_molec, (idx_atm_init, idx_atm_end) in enumerate(zip(molec_idx[:-1], molec_idx[1:])):
        # if monoatomic molecule
        if idx_atm_end == idx_atm_init + 1:
            v = vel[idx_atm_init]
            m = masses[idx_atm_init]
            kin_tot[idx_molec] = .5 * m * np.dot(v, v)
            kin_tra[idx_molec] = kin_tot[idx_molec]
            kin_rni[idx_molec] = 0
            kin_rot[idx_molec] = 0
            kin_int[idx_molec] = 0
            continue
        # compute center of mass position, velocity and total mass
        com_r = np.zeros(3)
        com_v = np.zeros(3)
        com_m = 0
        # loop over atoms in molecule
        for r, v, m in zip(pos[idx_atm_init:idx_atm_end],
                           vel[idx_atm_init:idx_atm_end],
                           masses[idx_atm_init:idx_atm_end]):
            com_r += m*r
            com_v += m*v
            com_m += m

            # total kinetic energy is straightforward
            kin_tot[idx_molec] += .5 * m * np.dot(v, v)

        com_r /= com_m
        com_v /= com_m

        # translational kinetic energy
        kin_tra[idx_molec] = .5 * com_m * np.dot(com_v, com_v)
        # combined rotational and internal kinetic energy
        kin_rni[idx_molec] = kin_tot[idx_molec] - kin_tra[idx_molec]

        # compute tensor of inertia and angular momentum
        inertia = np.zeros((3, 3))
        angular_mom = np.zeros(3)
        # loop over atoms in molecule
        for r, v, m in zip(pos[idx_atm_init:idx_atm_end],
                           vel[idx_atm_init:idx_atm_end],
                           masses[idx_atm_init:idx_atm_end]):
            # relative positions and velocities
            rr = r - com_r
            rv = v - com_v
            rr2 = np.dot(rr, rr)
            # inertia tensor:
            #   (i,i) = m*(r*r - r(i)*r(i))
            #   (i,j) = m*r(i)*r(j) (i != j)
            atm_inertia = -m*np.tensordot(rr, rr, axes=0)
            for i in range(3):
                atm_inertia[i][i] += m*rr2
            inertia += atm_inertia
            # angular momentum: r x p
            angular_mom += m * np.cross(rr, rv)

        # angular velocity of the molecule: inertia^{-1} * angular_mom
        angular_v = np.dot(np.linalg.inv(inertia), angular_mom)

        # test_kin = 0
        # for r, v, m in zip(pos[idx_atm_init:idx_atm_end],
        #                    vel[idx_atm_init:idx_atm_end],
        #                    masses[idx_atm_init:idx_atm_end]):
        #     # relative positions and velocities
        #     rr = r - com_r
        #     rv = v - com_v - np.cross(angular_v, rr)
        #     test_kin += .5 * m * np.dot(rv, rv)
        #
        # print(test_kin)

        kin_rot[idx_molec] = .5 * np.dot(angular_v, angular_mom)
        kin_int[idx_molec] = kin_rni[idx_molec] - kin_rot[idx_molec]

        # end loop over molecules

    return {'tot': kin_tot,
            'tra': kin_tra,
            'rni': kin_rni,
            'rot': kin_rot,
            'int': kin_int}


def group_kinetic_energy(kin_molec, nmolecs, molec_group=None):
    r"""
    Sums up the partitioned kinetic energy for a
    given group or the entire system.

    Parameters
    ----------
    kin_molec : List[dict]
        Partitioned kinetic energies per molecule.
    nmolecs : int
        Total number of molecules in the system.
    molec_group : iterable
        Indeces of the group to be summed up. None defaults to all molecules
        in the system. Default: None.

    Returns
    -------
    kin : dict
        Dictionary of partitioned kinetic energy for the group.
    """
    #
    kin = {'tot': 0, 'tra': 0, 'rni': 0, 'rot': 0, 'int': 0}
    #
    if molec_group is None:
        molec_group = range(nmolecs)
    # loop over molecules
    for idx_molec in molec_group:
        for key in kin.keys():
            kin[key] += kin_molec[key][idx_molec]

    return kin


def group_ndof(ndof_molec, nmolecs, molec_group=None):
    r"""
    Sums up the partitioned degrees of freedom for a
    given group or the entire system.

    Parameters
    ----------
    ndof_molec : List[dict]
        Partitioned degrees of freedom per molecule.
    nmolecs : int
        Total number of molecules in the system.
    molec_group : iterable
        Indeces of the group to be summed up. None defaults to all molecules
        in the system. Default: None.

    Returns
    -------
    ndof : dict
        Dictionary of partitioned degrees of freedom for the group.
    """
    #
    ndof = {'tot': 0, 'tra': 0, 'rni': 0, 'rot': 0, 'int': 0}
    #
    if molec_group is None:
        molec_group = range(nmolecs)
    # loop over molecules
    for idx_molec in molec_group:
        for key in ndof.keys():
            ndof[key] += ndof_molec[idx_molec][key]

    return ndof


def calc_temperatures(kin_molec, ndof_molec, nmolecs, molec_group=None):
    r"""
    Calculates the partitioned temperature for a
    given group or the entire system.

    Parameters
    ----------
    kin_molec : List[dict]
        Partitioned kinetic energies per molecule.
    ndof_molec : List[dict]
        Partitioned degrees of freedom per molecule.
    nmolecs : int
        Total number of molecules in the system.
    molec_group : iterable
        Indeces of the group to be summed up. None defaults to all molecules
        in the system. Default: None.

    Returns
    -------
    temp : dict
        Dictionary of partitioned temperatures for the group.
    """

    kin = group_kinetic_energy(kin_molec, nmolecs, molec_group)
    ndof = group_ndof(ndof_molec, nmolecs, molec_group)

    temp = {}
    for key in kin:
        temp[key] = temperature(kin[key], ndof[key])

    return temp


def test_mb_dist(kin_molec, ndof_molec, nmolecs,
                 temp, alpha, dict_keys, group=None,
                 verbosity=0, screen=False, filename=None,
                 ene_unit=None):
    r"""
    Tests if the partitioned kinetic energy trajectory of a group (or,
    if group is None, of the entire system) are separately Maxwell-Boltzmann
    distributed.

    Parameters
    ----------
    kin_molec : List[List[dict]]
        Partitioned kinetic energies per molecule for every frame.
    ndof_molec : List[dict]
        Partitioned degrees of freedom per molecule.
    nmolecs : int
        Total number of molecules in the system.
    temp : float
        Target temperature of the simulation.
    alpha : float
        Confidence for Maxwell-Boltzmann test.
    dict_keys : List[str]
        List of dictionary keys representing the partitions of the degrees
        of freedom.
    group : iterable
        Indeces of the group to be tested. None defaults to all molecules
        in the system. Default: None.
    verbosity : int
        Verbosity level, where 0 is quiet and 3 very chatty. Default: 0.
    screen : bool
        Plot distributions on screen. Default: False.
    filename : string
        Plot distributions to `filename`.pdf. Default: None.
    ene_unit : string
        Energy unit - used for plotting only.

    Returns
    -------
    result : List[float]
        p value for every partition
    """
    # save the partitioned kinetic energy trajectories
    group_kin = {key: [] for key in dict_keys}
    ndof = group_ndof(ndof_molec, nmolecs, group)
    # loop over frames
    for k in kin_molec:
        frame_kin = group_kinetic_energy(k, nmolecs, group)
        for key in dict_keys:
            group_kin[key].append(frame_kin[key])

    result = []
    failed = 0
    # test tot, tra, rni, rot, int
    if verbosity > 1:
        print('Testing whether\n'
              '* total,\n'
              '* translational,\n'
              '* rotational & internal,\n'
              '* rotational, and\n'
              '* internal\n'
              'kinetic energies are Maxwell-Boltzmann distributed.')
    elif verbosity > 0:
        print('Testing whether kinetic energies are Maxwell-Boltzmann distributed.')

    for key in dict_keys:
        p = check_mb_ensemble(kin=group_kin[key], temp=temp, ndof=ndof[key],
                              alpha=alpha, verbose=(verbosity > 2),
                              screen=screen, filename=filename+'_'+key,
                              ene_unit=ene_unit)
        result.append(p)
        if alpha is not None and p < alpha:
            failed += 1
        if verbosity > 1 and alpha is not None:
            if p >= alpha:
                print('* {}: passed'.format(key))
            else:
                print('* {}: failed'.format(key))

    if verbosity > 0 and alpha is not None:
        if failed == 0:
            print('-> Passed')
        else:
            print('-> Failed')

    return result


def test_temp_diff(kin_molec, ndof_molec, nmolecs,
                   dtemp, dict_keys, group=None,
                   verbosity=0, screen=False, filename=None,
                   ene_unit=None):
    r"""
    Tests if the partitioned temperatures (averaged over a trajectory)
    of a group (or, if group is None, of the entire system) are within a
    range `dtemp` of the total temperature.

    Parameters
    ----------
    kin_molec : List[List[dict]]
        Partitioned kinetic energies per molecule for every frame.
    ndof_molec : List[dict]
        Partitioned degrees of freedom per molecule.
    nmolecs : int
        Total number of molecules in the system.
    dtemp : float
        Target temperature of the simulation.
    dict_keys : List[str]
        List of dictionary keys representing the partitions of the degrees
        of freedom.
    group : iterable
        Indeces of the group to be tested. None defaults to all molecules
        in the system. Default: None.
    verbosity : int
        Verbosity level, where 0 is quiet and 3 very chatty. Default: 0.
    screen : bool
        Plot distributions on screen. Default: False.
    filename : string
        Plot distributions to `filename`.pdf. Default: None.
    ene_unit : string
        Energy unit - used for plotting only.

    Returns
    -------
    result : List[float]
        Temperature ratio to the total temperature for every partition.
    """
    # save the partitioned temperature trajectories
    group_temp = {key: [] for key in dict_keys}
    # loop over frames
    for k in kin_molec:
        frame_temp = calc_temperatures(k, ndof_molec, nmolecs, group)
        for key in dict_keys:
            group_temp[key].append(frame_temp[key])
    # average temperature
    group_temp_avg = {}
    for key in dict_keys:
        group_temp_avg[key] = np.mean(group_temp[key])

    failed = 0
    result = []
    if verbosity > 0 and dtemp is not None:
        print('Testing whether temperatures')
        print('  ' + str(dict_keys[1:]))
        print('are within {:.1f}% of temperature'.format(dtemp*100))
        print('  ' + dict_keys[0] + ' (' + str(group_temp_avg[dict_keys[0]]) + ')')
    elif verbosity > 0:
        print('Testing difference between temperatures')
        print('  ' + str(dict_keys[1:]))
        print('and')
        print('  ' + dict_keys[0] + ' (' + str(group_temp_avg[dict_keys[0]]) + ')')

    temp0 = group_temp_avg[dict_keys[0]]
    for key in dict_keys[1:]:
        result.append(group_temp_avg[key] / temp0)
        if dtemp is not None:
            if (1 - dtemp) * temp0 <= group_temp_avg[key] <= (1 + dtemp) * temp0:
                if verbosity > 1:
                    print('* {} ({:f}): passed'.format(key, group_temp_avg[key]))
            else:
                failed += 1
                if verbosity > 1:
                    print('* {} ({:f}): failed'.format(key, group_temp_avg[key]))

    if verbosity > 0 and dtemp is not None:
        if failed == 0:
            print('-> Passed')
        else:
            print('-> Failed')

    do_plot = screen or filename is not None
    if do_plot:
        data = []
        for key in dict_keys:
            t = group_temp[key]
            data.append({'x': range(0, len(t)),
                         'y': t,
                         'name': 'T(' + key + ')'})

        unit = ''
        if ene_unit is not None:
            unit = ' [' + ene_unit + ']'
        plot.plot(data,
                  legend='best',
                  title='Temperature trajectories',
                  xlabel='Frames',
                  ylabel='Temperature' + unit,
                  filename=filename,
                  screen=screen)

    return result


def test_temp_diff_groups(kin_molec, ndof_molec, nmolecs,
                          group1, group2,
                          dtemp, dict_keys, verbosity=0):
    r"""
    Tests if the partitioned temperatures (averaged over a trajectory)
    of two groups are (individually) within a range `dtemp` of each others.

    Parameters
    ----------
    kin_molec : List[List[dict]]
        Partitioned kinetic energies per molecule for every frame.
    ndof_molec : List[dict]
        Partitioned degrees of freedom per molecule.
    nmolecs : int
        Total number of molecules in the system.
    group1 : iterable
        Indeces of the first group to be compared.
    group2 : iterable
        Indeces of the second group to be compared.
    dtemp : float
        Target temperature of the simulation.
    dict_keys : List[str]
        List of dictionary keys representing the partitions of the degrees
        of freedom.
    verbosity : int
        Verbosity level, where 0 is quiet and 3 very chatty. Default: 0.

    Returns
    -------
    result : List[float]
        Temperature ratio (first group / second group) for every partition.
    """
    # save the partitioned temperature trajectories (group1)
    group1_temp = {key: [] for key in dict_keys}
    # loop over frames
    for k in kin_molec:
        frame_temp = calc_temperatures(k, ndof_molec, nmolecs, group1)
        for key in dict_keys:
            group1_temp[key].append(frame_temp[key])
    # average temperature
    for key in dict_keys:
        group1_temp[key] = np.mean(group1_temp[key])

    # save the partitioned temperature trajectories (group2)
    group2_temp = {key: [] for key in dict_keys}
    # loop over frames
    for k in kin_molec:
        frame_temp = calc_temperatures(k, ndof_molec, nmolecs, group2)
        for key in dict_keys:
            group2_temp[key].append(frame_temp[key])
    # average temperature
    for key in dict_keys:
        group2_temp[key] = np.mean(group2_temp[key])

    result = []
    failed = 0
    if verbosity > 0 and dtemp is not None:
        print('Testing whether temperatures of both groups are within {:.1f}%'.
              format(dtemp*100))

    for key in dict_keys:
        result.append(group1_temp[key]/group2_temp[key])
        if dtemp is not None:
            if (1. - dtemp) <= group1_temp[key]/group2_temp[key] <= (1 + dtemp):
                if verbosity > 1:
                    print('* {}: passed'.format(key))
            else:
                failed += 1
                if verbosity > 1:
                    print('* {}: failed'.format(key))

    if verbosity > 0 and dtemp is not None:
        if failed == 0:
            print('-> Passed')
        else:
            print('-> Failed')

    return result
