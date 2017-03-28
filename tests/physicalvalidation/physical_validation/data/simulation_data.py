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
import warnings
import numpy as np

import physical_validation.util.error as pv_error


class SimulationData(object):
    r"""
    TODO!
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

    def __init__(self):
        self.__ensemble = None
        self.__observables = ObservableData()
        self.__trajectory = TrajectoryData()
        self.__topology = TopologyData()
        self.__units = None
        self.__dt = 0

    @property
    def ensemble(self):
        r"""EnsembleData: Information on the sampled ensemble
    
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
    
        """
        return self.__trajectory

    @trajectory.setter
    def trajectory(self, trajectory):
        if not isinstance(trajectory, TrajectoryData):
            raise TypeError('No known conversion from ' + type(trajectory) +
                            'to TrajectoryData')
        self.__trajectory = trajectory

    @property
    def topology(self):
        r"""TopologyData: Information on the system's topology 
    
        """
        return self.__topology

    @topology.setter
    def topology(self, topology):
        if not isinstance(topology, TopologyData):
            raise TypeError('No known conversion from ' + type(topology) +
                            'to TopologyData')
        self.__topology = topology

    @property
    def dt(self):
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


class UnitData(object):
    r"""
    Class holding information about the units used
    """

    def __init__(self, kb,
                 energy_str, energy_conversion,
                 length_str, length_conversion,
                 volume_str, volume_conversion,
                 pressure_str, pressure_conversion,
                 time_str, time_conversion):

        self.__kb = float(kb)
        self.__energy_str = str(energy_str)
        self.__energy_conversion = float(energy_conversion)
        self.__length_str = str(length_str)
        self.__length_conversion = float(length_conversion)
        self.__volume_str = str(volume_str)
        self.__volume_conversion = float(volume_conversion)
        self.__pressure_str = str(pressure_str)
        self.__pressure_conversion = float(pressure_conversion)
        self.__time_str = str(time_str)
        self.__time_conversion = float(time_conversion)

    def __eq__(self, other):
        if not isinstance(other, UnitData):
            return False

        return (self.kb == other.kb and
                self.energy_conversion == other.energy_conversion and
                self.length_conversion == other.length_conversion and
                self.volume_conversion == other.volume_conversion and
                self.pressure_conversion == other.pressure_conversion and
                self.time_conversion == other.time_conversion)

    @property
    def kb(self):
        """float: The value of the Boltzmann constant"""
        return self.__kb

    @property
    def energy_str(self):
        """str: Energy unit"""
        return self.__energy_str

    @property
    def length_str(self):
        """str: Length unit"""
        return self.__length_str

    @property
    def volume_str(self):
        """str: Volume unit"""
        return self.__volume_str

    @property
    def pressure_str(self):
        """str: Pressure unit"""
        return self.__pressure_str

    @property
    def time_str(self):
        """str: Time unit"""
        return self.__time_str

    @property
    def energy_conversion(self):
        """float: Energy conversion factor: 1 ene_unit = energy_conversion * kJ/mol"""
        return self.__energy_conversion

    @property
    def length_conversion(self):
        """float: Length conversion factor: 1 length_unit = length_conversion * nm"""
        return self.__length_conversion

    @property
    def volume_conversion(self):
        """float: Volume conversion factor: 1 volume_unit = volume_conversion * nm^3"""
        return self.__volume_conversion

    @property
    def pressure_conversion(self):
        """float: Pressure conversion factor: 1 pressure_unit = pressure_conversion * bar"""
        return self.__pressure_conversion

    @property
    def time_conversion(self):
        """float: Time conversion factor: 1 time_unit = time_conversion * ps"""
        return self.__time_conversion


class EnsembleData(object):
    r"""
    Class holding ensemble information
    """

    @staticmethod
    def ensembles():
        return ('NVE',
                'NVT',
                'NPT',
                'muVT')

    def __init__(self, ensemble,
                 natoms=None, mu=None,
                 volume=None, pressure=None,
                 energy=None, temperature=None):
        self._ensemble = None
        self._n = None
        self._mu = None
        self._v = None
        self._p = None
        self._e = None
        self._t = None

        if ensemble not in self.ensembles():
            raise pv_error.InputError('ensemble',
                                      'Given ensemble unknown.')
        self._ensemble = ensemble

        if ensemble == 'NVE':
            if natoms is None:
                warnings.warn(ensemble + ' with undefined natoms.')
            if volume is None:
                warnings.warn(ensemble + ' with undefined volume.')
            # if energy is None:
            #     warnings.warn(ensemble + ' with undefined energy.')
            self._n = natoms
            self._v = volume
            self._e = energy
        if ensemble == 'NVT':
            if natoms is None:
                warnings.warn(ensemble + ' with undefined natoms.')
            if volume is None:
                warnings.warn(ensemble + ' with undefined volume.')
            if temperature is None:
                warnings.warn(ensemble + ' with undefined temperature.')
            self._n = natoms
            self._v = volume
            self._t = temperature
        if ensemble == 'NPT':
            if natoms is None:
                warnings.warn(ensemble + ' with undefined natoms.')
            if pressure is None:
                warnings.warn(ensemble + ' with undefined pressure.')
            if temperature is None:
                warnings.warn(ensemble + ' with undefined temperature.')
            self._n = natoms
            self._p = pressure
            self._t = temperature
        if ensemble == 'muVT':
            if mu is None:
                warnings.warn(ensemble + ' with undefined mu.')
            if volume is None:
                warnings.warn(ensemble + ' with undefined volume.')
            if temperature is None:
                warnings.warn(ensemble + ' with undefined temperature.')
            self._mu = mu
            self._v = volume
            self._t = temperature

    @property
    def ensemble(self):
        """Get ensemble"""
        return self._ensemble

    @property
    def natoms(self):
        """Get natoms"""
        return self._n

    @property
    def mu(self):
        """Get mu"""
        return self._mu

    @property
    def volume(self):
        """Get volume"""
        return self._v

    @property
    def pressure(self):
        """Get pressure"""
        return self._p

    @property
    def energy(self):
        """Get energy"""
        return self._e

    @property
    def temperature(self):
        """Get temperature"""
        return self._t


class TrajectoryData(object):
    r"""
    Class holding trajectory data.
    """

    @staticmethod
    def trajectories():
        return ('position',
                'velocity')

    def __init__(self, position=None, velocity=None):
        self._position = None
        self._velocity = None
        self._nframes = 0

        self.__getters = {
            'position': TrajectoryData.position.__get__,
            'velocity': TrajectoryData.velocity.__get__
        }

        if position is not None:
            position = np.array(position)
            if position.ndim == 2:
                # create 3-dimensional array
                position = np.array([position])
            if position.ndim != 3:
                warnings.warn('Expected 2- or 3-dimensional array.')
            self._nframes = position.shape[0]
            self._position = position

        if velocity is not None:
            velocity = np.array(velocity)
            if velocity.ndim == 2:
                # create 3-dimensional array
                velocity = np.array([velocity])
            if velocity.ndim != 3:
                warnings.warn('Expected 2- or 3-dimensional array.')
            if self._nframes != velocity.shape[0] and position is not None:
                raise pv_error.InputError(['position', 'velocity'],
                                          'Expected equal number of frames.')

            self._velocity = velocity

    def get(self, key):
        return self[key]

    def __getitem__(self, key):
        if key not in self.trajectories():
            raise KeyError
        return self.__getters[key](self)

    @property
    def position(self):
        """Get position"""
        return self._position

    @property
    def velocity(self):
        """Get velocity"""
        return self._velocity

    @property
    def nframes(self):
        """Get number of frames"""
        return self._nframes


class ObservableData(object):
    @staticmethod
    def observables():
        return ['kinetic_energy',
                'potential_energy',
                'total_energy',
                'volume',
                'pressure',
                'temperature',
                'constant_of_motion',
                'kin_per_molec']

    def __init__(self):
        self.__kinetic_energy = None
        self.__potential_energy = None
        self.__total_energy = None
        self.__volume = None
        self.__pressure = None
        self.__temperature = None
        self.__constant_of_motion = None
        self.__nframes = -1
        self.__kinetic_energy_per_molec = None

        self.__getters = {
            'kinetic_energy': ObservableData.kinetic_energy.__get__,
            'potential_energy': ObservableData.potential_energy.__get__,
            'total_energy': ObservableData.total_energy.__get__,
            'volume': ObservableData.volume.__get__,
            'pressure': ObservableData.pressure.__get__,
            'temperature': ObservableData.temperature.__get__,
            'constant_of_motion': ObservableData.constant_of_motion.__get__,
            'kin_per_molec': ObservableData.kinetic_energy_per_molecule.__get__
        }

        self.__setters = {
            'kinetic_energy': ObservableData.kinetic_energy.__set__,
            'potential_energy': ObservableData.potential_energy.__set__,
            'total_energy': ObservableData.total_energy.__set__,
            'volume': ObservableData.volume.__set__,
            'pressure': ObservableData.pressure.__set__,
            'temperature': ObservableData.temperature.__set__,
            'constant_of_motion': ObservableData.constant_of_motion.__set__,
            'kin_per_molec': ObservableData.kinetic_energy_per_molecule.__set__
        }

    def get(self, key):
        return self[key]

    def __getitem__(self, key):
        if key not in self.observables():
            raise KeyError
        return self.__getters[key](self)

    def set(self, key, value):
        self[key] = value

    def __setitem__(self, key, value):
        if key not in self.observables():
            raise KeyError
        self.__setters[key](self, value)

    @property
    def kinetic_energy(self):
        """Get kinetic_energy"""
        return self.__kinetic_energy

    @kinetic_energy.setter
    def kinetic_energy(self, kinetic_energy):
        """Set kinetic_energy"""
        if kinetic_energy is None:
            self.__kinetic_energy = None
            return
        kinetic_energy = np.array(kinetic_energy)
        if kinetic_energy.ndim != 1:
            raise pv_error.InputError('kinetic_energy',
                                      'Expected 1-dimensional array.')
        if self.nframes == -1:
            self.__nframes = kinetic_energy.size
        elif self.nframes != kinetic_energy.size:
            warnings.warn('Mismatch in number of frames. '
                          'Setting `nframes = None`.')
            self.__nframes = None
        self.__kinetic_energy = kinetic_energy

    @property
    def potential_energy(self):
        """Get potential_energy"""
        return self.__potential_energy

    @potential_energy.setter
    def potential_energy(self, potential_energy):
        """Set potential_energy"""
        if potential_energy is None:
            self.__potential_energy = None
            return
        potential_energy = np.array(potential_energy)
        if potential_energy.ndim != 1:
            raise pv_error.InputError('potential_energy',
                                      'Expected 1-dimensional array.')
        if self.nframes == -1:
            self.__nframes = potential_energy.size
        elif self.nframes != potential_energy.size:
            warnings.warn('Mismatch in number of frames. '
                          'Setting `nframes = None`.')
            self.__nframes = None
        self.__potential_energy = potential_energy

    @property
    def total_energy(self):
        """Get total_energy"""
        return self.__total_energy

    @total_energy.setter
    def total_energy(self, total_energy):
        """Set total_energy"""
        if total_energy is None:
            self.__total_energy = None
            return
        total_energy = np.array(total_energy)
        if total_energy.ndim != 1:
            raise pv_error.InputError('total_energy',
                                      'Expected 1-dimensional array.')
        if self.nframes == -1:
            self.__nframes = total_energy.size
        elif self.nframes != total_energy.size:
            warnings.warn('Mismatch in number of frames. '
                          'Setting `nframes = None`.')
            self.__nframes = None
        self.__total_energy = total_energy

    @property
    def volume(self):
        """Get volume"""
        return self.__volume

    @volume.setter
    def volume(self, volume):
        """Set volume"""
        if volume is None:
            self.__volume = None
            return
        volume = np.array(volume)
        if volume.ndim != 1:
            raise pv_error.InputError('volume',
                                      'Expected 1-dimensional array.')
        if self.nframes == -1:
            self.__nframes = volume.size
        elif self.nframes != volume.size:
            warnings.warn('Mismatch in number of frames. '
                          'Setting `nframes = None`.')
            self.__nframes = None
        self.__volume = volume

    @property
    def pressure(self):
        """Get pressure"""
        return self.__pressure

    @pressure.setter
    def pressure(self, pressure):
        """Set pressure"""
        if pressure is None:
            self.__pressure = None
            return
        pressure = np.array(pressure)
        if pressure.ndim != 1:
            raise pv_error.InputError('pressure',
                                      'Expected 1-dimensional array.')
        if self.nframes == -1:
            self.__nframes = pressure.size
        elif self.nframes != pressure.size:
            warnings.warn('Mismatch in number of frames. '
                          'Setting `nframes = None`.')
            self.__nframes = None
        self.__pressure = pressure

    @property
    def temperature(self):
        """Get temperature"""
        return self.__temperature

    @temperature.setter
    def temperature(self, temperature):
        """Set temperature"""
        if temperature is None:
            self.__temperature = None
            return
        temperature = np.array(temperature)
        if temperature.ndim != 1:
            raise pv_error.InputError('temperature',
                                      'Expected 1-dimensional array.')
        if self.nframes == -1:
            self.__nframes = temperature.size
        elif self.nframes != temperature.size:
            warnings.warn('Mismatch in number of frames. '
                          'Setting `nframes = None`.')
            self.__nframes = None
        self.__temperature = temperature

    @property
    def constant_of_motion(self):
        """Get constant_of_motion"""
        return self.__constant_of_motion

    @constant_of_motion.setter
    def constant_of_motion(self, constant_of_motion):
        """Set constant_of_motion"""
        if constant_of_motion is None:
            self.__constant_of_motion = None
            return
        constant_of_motion = np.array(constant_of_motion)
        if constant_of_motion.ndim != 1:
            raise pv_error.InputError('constant_of_motion',
                                      'Expected 1-dimensional array.')
        if self.nframes == -1:
            self.__nframes = constant_of_motion.size
        elif self.nframes != constant_of_motion.size:
            warnings.warn('Mismatch in number of frames. '
                          'Setting `nframes = None`.')
            self.__nframes = None
        self.__constant_of_motion = constant_of_motion

    @property
    def nframes(self):
        """Get number of frames"""
        if self.__nframes is None:
            warnings.warn('A mismatch in the number of frames between observables '
                          'was detected. Setting `nframes = None`.')
        return self.__nframes

    @property
    def kinetic_energy_per_molecule(self):
        """Get kinetic_energy"""
        return self.__kinetic_energy_per_molec

    @kinetic_energy_per_molecule.setter
    def kinetic_energy_per_molecule(self, kinetic_energy):
        """Set kinetic_energy"""
        if kinetic_energy is None:
            self.__kinetic_energy_per_molec = None
            return
        # used internally - check for consistency?
        self.__kinetic_energy_per_molec = kinetic_energy


class TopologyData(object):
    r"""Class holding toplogical information.

    """

    def __init__(self):
        self._natoms = None
        self._mass = None
        self._nconstraints = None
        self._ndof_total = None
        self._ndof_reduction_tra = None
        self._ndof_reduction_rot = None
        self._molecule_idx = None
        self._nconstraints_per_molecule = None
        self._ndof_per_molecule = None
        self._bonds = None
        self._constrained_bonds = None

    @property
    def natoms(self):
        """int: Number of atoms in the system

        """
        return self._natoms

    @natoms.setter
    def natoms(self, natoms):
        self._natoms = int(natoms)

    @property
    def mass(self):
        """nd-array: Mass vector for the atoms

        Setter accepts array-like objects.

        """
        return self._mass

    @mass.setter
    def mass(self, mass):
        mass = np.asarray(mass)
        if mass.ndim != 1:
            raise pv_error.InputError('mass',
                                      'Expected 1-dimensional array.')
        if self.natoms is None:
            self.natoms = mass.size
        elif mass.size != self.natoms:
            raise pv_error.InputError('mass',
                                      'Mass vector does not have length == natoms.')
        self._mass = mass

    @property
    def nconstraints(self):
        """float: Total number of constraints in the system 
        
        Does not include the reduction of degrees of freedom in the absence of
        external forces.

        """
        return self._nconstraints

    @nconstraints.setter
    def nconstraints(self, nconstraints):
        self._nconstraints = float(nconstraints)

    @property
    def ndof_total(self):
        """float: Total number of degrees of freedom in the system 
        
        """
        return self._ndof_total

    @ndof_total.setter
    def ndof_total(self, ndof_total):
        self._ndof_total = float(ndof_total)

    @property
    def ndof_reduction_tra(self):
        """float: Number of translational degrees of freedom deducted 
        from 3*[# of molecules]

        """
        return self._ndof_reduction_tra

    @ndof_reduction_tra.setter
    def ndof_reduction_tra(self, ndof_reduction_tra):
        self._ndof_reduction_tra = float(ndof_reduction_tra)

    @property
    def ndof_reduction_rot(self):
        """float: Number of rotational degrees of freedom deducted 
        from 3*[# of molecules]

        """
        return self._ndof_reduction_rot

    @ndof_reduction_rot.setter
    def ndof_reduction_rot(self, ndof_reduction_rot):
        self._ndof_reduction_rot = float(ndof_reduction_rot)

    @property
    def molecule_idx(self):
        """nd-array: List of index of first atom of each molecule

        Setter accepts array-like objects.

        """
        return self._molecule_idx

    @molecule_idx.setter
    def molecule_idx(self, molecule_idx):
        molecule_idx = np.asarray(molecule_idx)
        if molecule_idx.ndim != 1:
            raise pv_error.InputError('molecule_idx',
                                      'Expected 1-dimensional array.')
        if (self.nconstraints_per_molecule is not None and
           self.nconstraints_per_molecule.shape != molecule_idx.shape):
            warnings.warn('New `molecule_idx` does not have the same'
                          'shape as previously set `nconstraints_per_molecule`.'
                          'Setting `nconstraints_per_molecule = None` to avoid'
                          'errors.')
            self._nconstraints_per_molecule = None
        self._molecule_idx = molecule_idx

    @property
    def nconstraints_per_molecule(self):
        """nd-array: List of number of constraints per molecule

        Setter accepts array-like objects.

        """
        return self._nconstraints_per_molecule

    @nconstraints_per_molecule.setter
    def nconstraints_per_molecule(self, nconstraints_per_molecule):
        nconstraints_per_molecule = np.array(nconstraints_per_molecule)
        if nconstraints_per_molecule.ndim != 1:
            raise pv_error.InputError('nconstraints_per_molecule',
                                      'Expected 1-dimensional array.')
        if self.molecule_idx is not None:
            if nconstraints_per_molecule.shape != self.molecule_idx.shape:
                raise pv_error.InputError('nconstraints_per_molecule',
                                          'Expected `nconstraints_per_molecule` to have'
                                          'the same shape as `moldecule_idx`.')

        self._nconstraints_per_molecule = nconstraints_per_molecule

    @property
    def ndof_per_molecule(self):
        """nd-array: List of number of degrees of freedom per molecule

        Setter accepts array-like objects.

        """
        return self._ndof_per_molecule

    @ndof_per_molecule.setter
    def ndof_per_molecule(self, ndof_per_molecule):
        # used internally - check for consistency?
        self._ndof_per_molecule = ndof_per_molecule

    @property
    def bonds(self):
        """List[List[int]]: List of bonds per molecule
        """
        return self._bonds

    @bonds.setter
    def bonds(self, bonds):
        self._bonds = bonds

    @property
    def constrained_bonds(self):
        """List[List[int]]: List of constrained bonds per molecule
        """
        return self._constrained_bonds

    @constrained_bonds.setter
    def constrained_bonds(self, constrained_bonds):
        self._constrained_bonds = constrained_bonds
