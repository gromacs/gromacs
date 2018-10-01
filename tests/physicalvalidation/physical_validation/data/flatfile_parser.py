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
flatfile_parser.py
"""
from . import parser
from . import SimulationData, TrajectoryData, ObservableData


class FlatfileParser(parser.Parser):
    """
    FlatfileParser
    """

    def __init__(self):
        super(FlatfileParser, self).__init__()

    def get_simulation_data(self, units=None, ensemble=None, system=None, dt=None,
                            position_file=None, velocity_file=None,
                            kinetic_ene_file=None, potential_ene_file=None,
                            total_ene_file=None, volume_file=None,
                            pressure_file=None, temperature_file=None,
                            const_of_mot_file=None):
        r"""Read simulation data from flat files

        Returns a SimulationData object created from (optionally) provided UnitData, EnsembleData
        and SystemData, as well as TrajectoryData and ObservableData objects created from flat
        files. The files are expected to be in one of the following formats:

        * xyz-format
          trajectory files (position_file, velocity_file)
          - three numbers per line, separated by white space
          - frames delimited by a completely blank line
          - any character after (and including) a '#' are ignored
        * 1d-format
          all other files
          - one number per line
          - any character after (and including) a '#' are ignored

        Parameters
        ----------
        units: UnitData, optional
            A UnitData object representing the units used in the simulation
        ensemble: EnsembleData, optional
            A EnsembleData object representing the ensemble the simulation has been performed in
        system: SystemData, optional
            A SystemData object representing the atoms and molecules in the system
        dt: float, optional
            The time step used in the simulation
        position_file: str, optional
            Path to a file in xyz-format containing the position trajectory
        velocity_file: str, optional
            Path to a file in xyz-format containing the velocity trajectory
        kinetic_ene_file: str, optional
            Path to a file in 1d-format containing the kinetic energy trajectory
        potential_ene_file: str, optional
            Path to a file in 1d-format containing the potential energy trajectory
        total_ene_file: str, optional
            Path to a file in 1d-format containing the total energy trajectory
        volume_file: str, optional
            Path to a file in 1d-format containing the volume trajectory
        pressure_file: str, optional
            Path to a file in 1d-format containing the pressure trajectory
        temperature_file: str, optional
            Path to a file in 1d-format containing the temperature trajectory
        const_of_mot_file: str, optional
            Path to a file in 1d-format containing the constant of motion trajectory

        Returns
        -------
        result: SimulationData
            A SimulationData filled with the provided ensemble and
            system objects as well as the trajectory data found in the
            edr and trr / gro files.

        """

        trj_dict = {
            'position': position_file,
            'velocity': velocity_file
        }

        if any(trj_dict.values()):
            trajectory = TrajectoryData()
            for key, filename in trj_dict.items():
                if filename is None:
                    continue
                trajectory[key] = self.__read_xyz(filename)
        else:
            trajectory = None

        obs_dict = {
            'kinetic_energy': kinetic_ene_file,
            'potential_energy': potential_ene_file,
            'total_energy': total_ene_file,
            'volume': volume_file,
            'pressure': pressure_file,
            'temperature': temperature_file,
            'constant_of_motion': const_of_mot_file
        }

        if any(obs_dict.values()):
            observables = ObservableData()
            for key, filename in obs_dict.items():
                if filename is None:
                    continue
                observables[key] = self.__read_1d(filename)
        else:
            observables = None

        result = SimulationData(units=units, dt=dt, system=system, ensemble=ensemble,
                                observables=observables, trajectory=trajectory)

        return result

    @staticmethod
    def __read_xyz(filename):
        result = []
        with open(filename) as f:
            frame = []
            for line in f:
                if not line:
                    # blank line
                    if frame:
                        result.append(frame)
                        frame = []
                    continue
                line = line.split('#', maxsplit=1)[0].strip()
                if not line:
                    # only comment on this line
                    continue
                xyz = line.split()[0:3]
                frame.append([float(n) for n in xyz])
            if frame:
                result.append(frame)
        return result

    @staticmethod
    def __read_1d(filename):
        result = []
        with open(filename) as f:
            for line in f:
                line = line.split('#', maxsplit=1)[0].strip()
                if not line:
                    # blank or comment-only line
                    continue
                result.append(float(line.strip()))
        return result
