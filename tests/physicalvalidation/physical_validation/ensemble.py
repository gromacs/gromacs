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
This software allows users to perform statistical test to determine if a
given molecular simulation is consistent with the thermodynamic ensemble
it is performed in.

Users should cite the JCTC paper: Shirts, M. R. "Simple Quantitative
Tests to Validate Sampling from Thermodynamic Ensembles",
J. Chem. Theory Comput., 2013, 9 (2), pp 909-926,
https://doi.org/10.1021/ct300688p
"""

import numpy as np

from .util import ensemble
from .data import SimulationData
from .util import error as pv_error


def check(data_sim_one, data_sim_two,
          total_energy=False,
          screen=False, filename=None,
          verbosity=1):
    r"""
    Check the ensemble. The correct check is inferred from the
    simulation data given.

    Parameters
    ----------
    data_sim_one : SimulationData
    data_sim_two : SimulationData
    total_energy : bool
    screen : bool
        Plot distributions on screen. Default: False.
    filename : string
        Plot distributions to `filename`.pdf. Default: None.
    verbosity : int
        Level of verbosity, from 0 (quiet) to 3 (very verbose).
        Default: 1

    Returns
    -------
    quantiles : List[float]
        The number of quantiles the computed result is off the analytical one.

    """
    if not SimulationData.compatible(data_sim_one,
                                     data_sim_two):
        raise pv_error.InputError(['data_sim_one', 'data_sim_two'],
                                  'Simulation data not compatible.')

    if data_sim_one.ensemble.ensemble != data_sim_two.ensemble.ensemble:
        raise pv_error.InputError(['data_sim_one', 'data_sim_two'],
                                  'The two simulations were sampling different ensembles. '
                                  'The simulations are expected to differ in state point '
                                  '(e.g. target temperature, target pressure), but not '
                                  'in their sampled ensemble (e.g. NVT, NPT).')

    sampled_ensemble = data_sim_one.ensemble.ensemble

    if sampled_ensemble == 'NVE' or sampled_ensemble == 'muVE':
        raise pv_error.InputError(['data_sim_one', 'data_sim_two'],
                                  'Test of ensemble ' + sampled_ensemble + ' is not implemented '
                                  '(yet).')

    if total_energy:
        eneq = 'E'
        e1 = data_sim_one.observables.total_energy
        e2 = data_sim_two.observables.total_energy
    else:
        eneq = 'U'
        e1 = data_sim_one.observables.potential_energy
        e2 = data_sim_two.observables.potential_energy

    quantiles = None

    if sampled_ensemble == 'NVT':
        quantiles = ensemble.check_1d(
            traj1=e1, traj2=e2,
            param1=data_sim_one.ensemble.temperature,
            param2=data_sim_two.ensemble.temperature,
            kb=data_sim_one.units.kb,
            quantity=eneq,
            dtemp=True, dpress=False,
            verbosity=verbosity,
            filename=filename, screen=screen
        )

    elif sampled_ensemble == 'NPT':
        temperatures = np.array([data_sim_one.ensemble.temperature,
                                 data_sim_two.ensemble.temperature])
        pressures = np.array([data_sim_one.ensemble.pressure,
                              data_sim_two.ensemble.pressure])
        equal_temps = temperatures[0] == temperatures[1]
        equal_press = pressures[0] == pressures[1]

        v1 = data_sim_one.observables.volume
        v2 = data_sim_two.observables.volume

        # Calculate conversion from p*V to energy units
        #
        # GROMACS standard units are
        #   energy: kJ/mol
        #   volume: nm^3
        #   pressure: bar
        #   => pV-term: bar * nm^3 == 1e-25 kJ == 6.022140857e-2 kJ/mol
        #   => pvconvert = 6.022140857e-2
        # UnitData stores conversion factors relative to GROMACS units
        #   energy: energy_conversion * kJ/mol
        #   volume: volume_conversion * nm^3
        #   pressure: pressure_conversion * bar
        #   => pV-term: [p]*[V] == pressure_conversion * volume_conversion bar * nm^3
        # Units were checked earlier, so we can use either simulation data structure
        pvconvert = 6.022140857e-2
        pvconvert *= (data_sim_one.units.pressure_conversion *
                      data_sim_one.units.volume_conversion)
        pvconvert /= data_sim_one.units.energy_conversion

        if equal_press and not equal_temps:
            e1 = e1 + pvconvert * pressures[0] * v1
            e2 = e2 + pvconvert * pressures[1] * v2
            if eneq == 'U':
                eneq = 'H'
            quantiles = ensemble.check_1d(
                traj1=e1, traj2=e2,
                param1=temperatures[0],
                param2=temperatures[1],
                kb=data_sim_one.units.kb,
                quantity=eneq,
                dtemp=True, dpress=False,
                verbosity=verbosity,
                filename=filename, screen=screen
            )
        elif equal_temps and not equal_press:
            quantiles = ensemble.check_1d(
                traj1=v1, traj2=v2,
                param1=pressures[0],
                param2=pressures[1],
                kb=data_sim_one.units.kb,
                quantity='V',
                dtemp=False, dpress=True,
                temp=temperatures[0],
                pvconvert=pvconvert,
                verbosity=verbosity,
                filename=filename, screen=screen
            )
        else:
            traj1 = np.array([e1, v1])
            traj2 = np.array([e2, v2])
            param1 = np.array([temperatures[0], pressures[0]])
            param2 = np.array([temperatures[1], pressures[1]])
            quantiles = ensemble.check_2d(
                traj1=traj1, traj2=traj2,
                param1=param1, param2=param2,
                kb=data_sim_one.units.kb,
                pvconvert=pvconvert,
                quantity=[eneq, 'V'],
                dtempdpress=True,
                verbosity=verbosity,
                filename=filename, screen=screen
            )

    return quantiles


def estimate_interval(data, verbosity=1, total_energy=False):
    r"""
    In order to perform an ensemble check, two simulations at distinct state
    point are needed. Choosing two state points too far apart will result
    in poor or zero overlap between the distributions, leading to very noisy
    results (due to sample errors in the tails) or a breakdown of the method,
    respectively. Choosing two state points very close to each others, on the
    other hand, makes it difficult to distinguish the slope from statistical
    error in the samples.

    This function implements a rule of thumb based on the standard deviations
    of distributions. It takes a single simulation and suggests appropriate
    intervals for a second simulation to be used for ensemble checking.

    Parameters
    ----------
    data : SimulationData
        The performed simulation.
    verbosity : int, optional
        If 0, no output is printed on screen. If 1, estimated intervals are
        printed. If larger, additional information during calculation are
        printed.
        Default: 1
    total_energy : bool, optional
        Use total energy instead of potential energy only.
        Default: False

    Returns
    -------
    intervals : Dict
        If `data` was performed under NVT conditions, `intervals` contains only
        one entry:

            * `'dT'`, containing the suggested temperature interval.

        If `data` was performed under NPT conditions, `intervals` contains three
        entries:

            * `'dT'`: Suggested temperature interval at constant pressure
            * `'dP'`: Suggested pressure interval at constant temperature
            * `'dTdP'`: Suggested combined temperature and pressure interval

    """

    if total_energy:
        ene = data.observables.total_energy
    else:
        ene = data.observables.potential_energy

    if data.ensemble.ensemble == 'NVT':
        result = ensemble.estimate_interval(
            ens_string='NVT',
            ens_temp=data.ensemble.temperature,
            energy=ene,
            kb=data.units.kb,
            verbosity=verbosity,
            tunit=data.units.temperature_str
        )
    elif data.ensemble.ensemble == 'NPT':
        pvconvert = 6.022140857e-2
        pvconvert *= (data.units.pressure_conversion *
                      data.units.volume_conversion)
        pvconvert /= data.units.energy_conversion
        result = ensemble.estimate_interval(
            ens_string='NPT',
            ens_temp=data.ensemble.temperature,
            energy=ene,
            kb=data.units.kb,
            ens_press=data.ensemble.pressure,
            volume=data.observables.volume,
            pvconvert=pvconvert,
            verbosity=verbosity,
            tunit=data.units.temperature_str,
            punit=data.units.pressure_str
        )
    else:
        raise NotImplementedError('estimate_interval() not implemented for ensemble ' +
                                  data.ensemble.ensemble)

    return result
