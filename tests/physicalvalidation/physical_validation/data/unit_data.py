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


class UnitData(object):
    r"""UnitData: Information about the units used

    The information about units consists of different parts:
    * The name of the units (energy_str, length_str, volume_str, pressure_str, time_str),
    * the value of kB in the used energy units, and
    * the conversion factor to GROMACS units (kJ/mol, nm, nm^3, bar, ps).
    The names are only used for output (console printing and plotting), and are optional.
    The conversion factors and kB are, on the other hand, used in computations and need
    to be given.
    """
    def __init__(self, kb, energy_conversion, length_conversion,
                 volume_conversion, pressure_conversion, time_conversion,
                 energy_str='ENE', length_str='LEN',
                 volume_str='VOL', pressure_str='PRESS', time_str='TIME'):

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
