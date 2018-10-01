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
"""
The `integratorconvergence` module is part of the physical_validation
package, and consists of checks of the convergence of the MD integrator.
"""
# py2 compatibility
from __future__ import absolute_import, division, print_function

from .util import integrator as util_integ
from .data import SimulationData
from .util import error as pv_error


def convergence(simulations,
                convergence_test='max_deviation',
                verbose=True, slope=False,
                screen=False, filename=None):
    r"""
    Compares the convergence of the fluctuations of conserved quantities
    with decreasing simulation time step to theoretical expectations.

    Parameters
    ----------
    simulations : list of SimulationData
        The (otherwise identical) simulations performed using different
        timesteps
    convergence_test : str, optional
        A function defining the convergence test. Currently, only one
        test is implemented:
        `max_deviation`
    verbose : bool, optional
    screen : bool
        Plot convergence on screen. Default: False.
    filename : string
        Plot convergence to `filename`.pdf. Default: None.

    Returns
    -------
    result : float

    Notes
    -----
    For a simplectic integration algorithm, the fluctuations
    :math:`\delta E` of a constant of motion :math:`E` (such as, for
    example, the total energy in a NVE simulations) are theoretically
    expected to scale like the squared timestep of the integration.
    When comparing two otherwise identical simulations performed at
    different time step :math:`\Delta t`, the following equality is
    hence expected to hold:

    .. math::
        \frac{\Delta t_1}{\Delta t_2} = \frac{\delta E_1}{\delta E_2}

    This function calculates the ratio of the fluctuation for simulations
    performed at different timesteps and compares it to the analytically
    expected value. If the deviation is larger than `tol`, the test is
    considered failed.

    """
    constant_of_motion = {}

    convergence_tests = {
        'max_deviation': util_integ.max_deviation
    }

    if convergence_test not in convergence_tests:
        raise pv_error.InputError('convergence_test',
                                  'Unknown convergence test.')

    convergence_test = convergence_tests[convergence_test]

    for s in simulations:
        if not isinstance(s, SimulationData):
            raise pv_error.InputError('simulations',
                                      'Expected a list of objects of type SimulationData')
        if s.dt <= 0:
            raise pv_error.InputError('simulations',
                                      'Found SimulationData with timestep dt<=0')
        key = str(s.dt)

        if key in constant_of_motion:
            raise pv_error.InputError('simulations',
                                      'Found two simulations with identical timestep')

        constant_of_motion[key] = s.observables.constant_of_motion

    return util_integ.check_convergence(constant_of_motion,
                                        convergence_test=convergence_test,
                                        verbose=verbose, slope=slope,
                                        screen=screen, filename=filename)
