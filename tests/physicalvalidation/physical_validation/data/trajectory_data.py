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


class TrajectoryData(object):
    r"""TrajectoryData: The position and velocity trajectory along the simulation

    The full trajectory is needed to calculate the equipartition of the kinetic energy.
    As they are used in connection, the position and velocity trajectories are expected
    to have the same shape and number of frames.

    The position and velocity trajectories can be accessed either using the getters
    of an object, as in
        trajectory.position
        trajectory.velocity
    or using the key notation, as in
        trajectory['position']
        trajectory['velocity']
    """

    @staticmethod
    def trajectories():
        return ('position',
                'velocity')

    def __init__(self, position=None, velocity=None):
        self.__position = None
        self.__velocity = None
        self.__nframes = 0

        if position is not None:
            self.position = position
        if velocity is not None:
            self.velocity = velocity

        self.__getters = {
            'position': TrajectoryData.position.__get__,
            'velocity': TrajectoryData.velocity.__get__
        }

        self.__setters = {
            'position': TrajectoryData.position.__set__,
            'velocity': TrajectoryData.velocity.__set__
        }

    def get(self, key):
        return self[key]

    def __getitem__(self, key):
        if key not in self.trajectories():
            raise KeyError
        return self.__getters[key](self)

    def set(self, key, value):
        self[key] = value

    def __setitem__(self, key, value):
        if key not in self.trajectories():
            raise KeyError
        self.__setters[key](self, value)

    @property
    def position(self):
        """Get position"""
        return self.__position

    @position.setter
    def position(self, pos):
        """Set position"""
        pos = np.array(pos)
        if pos.ndim == 2:
            # create 3-dimensional array
            pos = np.array([pos])
        if pos.ndim != 3:
            warnings.warn('Expected 2- or 3-dimensional array.')
        if self.__nframes == 0 and self.__velocity is None:
            self.__nframes = pos.shape[0]
        elif self.__nframes != pos.shape[0]:
            raise pv_error.InputError(['pos'],
                                      'Expected equal number of frames as in velocity trajectory.')
        self.__position = pos

    @property
    def velocity(self):
        """Get velocity"""
        return self.__velocity

    @velocity.setter
    def velocity(self, vel):
        """Set velocity"""
        vel = np.array(vel)
        if vel.ndim == 2:
            # create 3-dimensional array
            vel = np.array([vel])
        if vel.ndim != 3:
            warnings.warn('Expected 2- or 3-dimensional array.')
        if self.__nframes == 0 and self.__position is None:
            self.__nframes = vel.shape[0]
        elif self.__nframes != vel.shape[0]:
            raise pv_error.InputError(['vel'],
                                      'Expected equal number of frames as in position trajectory.')
        self.__velocity = vel

    @property
    def nframes(self):
        """Get number of frames"""
        return self.__nframes
