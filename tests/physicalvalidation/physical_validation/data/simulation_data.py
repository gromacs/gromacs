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
Data structures carrying simulation data.
"""
from ..util import error as pv_error
# py2.7 compatibility
from .unit_data import UnitData
from .ensemble_data import EnsembleData
from .system_data import SystemData
from .observable_data import ObservableData
from .trajectory_data import TrajectoryData
# replace lines above by this when py2.7 support is dropped:
# from . import UnitData, EnsembleData, SystemData, ObservableData, TrajectoryData


class SimulationData(object):
    r"""SimulationData: System information and simulation results

    The SimulationData class holds both the information on the system and
    the results of a simulation run of that system. SimulationData contains
    all information on a simulation run needed by the physical validation
    tests. SimulationData objects can either be created directly by calling
    the class constructor, or by using a parser returning a SimulationData
    object.
    """

    @staticmethod
    def compatible(data_1, data_2):
        r"""Checks whether two simulations are compatible for common validation.

        Parameters
        ----------
        data_1 : SimulationData
        data_2 : SimulationData

        Returns
        -------
        result : bool

        """
        if not isinstance(data_1, SimulationData):
            raise pv_error.InputError('data_1',
                                      'Expected type SimulationData')
        if not isinstance(data_2, SimulationData):
            raise pv_error.InputError('data_2',
                                      'Expected type SimulationData')

        return data_1.units == data_2.units

    def __init__(self, units=None, dt=None,
                 system=None, ensemble=None,
                 observables=None, trajectory=None):
        self.__units = None
        if units is not None:
            self.units = units
        self.__dt = 0
        if dt is not None:
            self.dt = dt
        self.__topology = None
        if system is not None:
            self.system = system
        self.__ensemble = None
        if ensemble is not None:
            self.ensemble = ensemble
        self.__observables = None
        if observables is not None:
            self.observables = observables
        self.__trajectory = None
        if trajectory is not None:
            self.trajectory = trajectory

    @property
    def ensemble(self):
        r"""EnsembleData: Information on the sampled ensemble

        Returns
        -------
        ensemble : EnsembleData
        """
        return self.__ensemble

    @ensemble.setter
    def ensemble(self, ensemble):
        if not isinstance(ensemble, EnsembleData):
            raise TypeError('No known conversion from ' + type(ensemble) +
                            'to EnsembleData')
        self.__ensemble = ensemble

    @property
    def units(self):
        r"""UnitsData: Information on the sampled units

        Returns
        -------
        units : UnitData
        """
        return self.__units

    @units.setter
    def units(self, units):
        if not isinstance(units, UnitData):
            raise TypeError('No known conversion from ' + type(units) +
                            'to UnitData')
        self.__units = units

    @property
    def observables(self):
        r"""ObservableData: Observables collected during the simulation

        Returns
        -------
        observables : ObservableData
        """
        return self.__observables

    @observables.setter
    def observables(self, observables):
        if not isinstance(observables, ObservableData):
            raise TypeError('No known conversion from ' + type(observables) +
                            'to ObservableData')
        self.__observables = observables

    @property
    def trajectory(self):
        r"""TrajectoryData: Trajectories collected during the simulation

        Returns
        -------
        trajectory : TrajectoryData
        """
        return self.__trajectory

    @trajectory.setter
    def trajectory(self, trajectory):
        if not isinstance(trajectory, TrajectoryData):
            raise TypeError('No known conversion from ' + type(trajectory) +
                            'to TrajectoryData')
        self.__trajectory = trajectory

    @property
    def system(self):
        r"""SystemData: Information on the system's system

        Returns
        -------
        system : SystemData
        """
        return self.__topology

    @system.setter
    def system(self, topology):
        if not isinstance(topology, SystemData):
            raise TypeError('No known conversion from ' + type(topology) +
                            'to SystemData')
        self.__topology = topology

    @property
    def dt(self):
        r""" The timestep of the simulation run.

        Returns
        -------
        timestep : float
        """
        return self.__dt

    @dt.setter
    def dt(self, dt):
        dt = float(dt)
        self.__dt = dt

    def set_ensemble(self, ensemble,
                     natoms=None, mu=None,
                     volume=None, pressure=None,
                     energy=None, temperature=None):
        self.__ensemble = EnsembleData(ensemble,
                                       natoms, mu,
                                       volume, pressure,
                                       energy, temperature)
