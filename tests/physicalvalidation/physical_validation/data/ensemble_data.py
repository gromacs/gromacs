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

from ..util import error as pv_error


class EnsembleData(object):
    r"""EnsembleData: Holds data defining the ensemble

    The ensemble is a string indicating the thermodynamical ensemble a simulation was
    performed in, and is any of 'NVE', 'NVT', 'NPT', 'muVT'.
    Depending on the ensemble, EnsembleData then holds additional information defining
    the ensemble, such as the number of particles N, the chemical potential mu, the
    volume V, the pressure P, the constant energy E or the temperature T. While any
    of these additional information are optional, most of them are needed by certain
    tests, such that not fully defining the ensemble results in warnings. The notable
    exception to this rule is the constant energy E for the NVE, which is not needed
    by any test and can hence be omitted without raising a warning.
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
        self.__ensemble = None
        self.__n = None
        self.__mu = None
        self.__v = None
        self.__p = None
        self.__e = None
        self.__t = None

        if ensemble not in self.ensembles():
            raise pv_error.InputError('ensemble',
                                      'Given ensemble unknown.')
        self.__ensemble = ensemble

        if ensemble == 'NVE':
            if natoms is None:
                warnings.warn(ensemble + ' with undefined natoms.')
            if volume is None:
                warnings.warn(ensemble + ' with undefined volume.')
            # if energy is None:
            #     warnings.warn(ensemble + ' with undefined energy.')
            self.__n = natoms
            self.__v = volume
            self.__e = energy
        if ensemble == 'NVT':
            if natoms is None:
                warnings.warn(ensemble + ' with undefined natoms.')
            if volume is None:
                warnings.warn(ensemble + ' with undefined volume.')
            if temperature is None:
                warnings.warn(ensemble + ' with undefined temperature.')
            self.__n = natoms
            self.__v = volume
            self.__t = temperature
        if ensemble == 'NPT':
            if natoms is None:
                warnings.warn(ensemble + ' with undefined natoms.')
            if pressure is None:
                warnings.warn(ensemble + ' with undefined pressure.')
            if temperature is None:
                warnings.warn(ensemble + ' with undefined temperature.')
            self.__n = natoms
            self.__p = pressure
            self.__t = temperature
        if ensemble == 'muVT':
            if mu is None:
                warnings.warn(ensemble + ' with undefined mu.')
            if volume is None:
                warnings.warn(ensemble + ' with undefined volume.')
            if temperature is None:
                warnings.warn(ensemble + ' with undefined temperature.')
            self.__mu = mu
            self.__v = volume
            self.__t = temperature

    @property
    def ensemble(self):
        """Get ensemble"""
        return self.__ensemble

    @property
    def natoms(self):
        """Get natoms"""
        return self.__n

    @property
    def mu(self):
        """Get mu"""
        return self.__mu

    @property
    def volume(self):
        """Get volume"""
        return self.__v

    @property
    def pressure(self):
        """Get pressure"""
        return self.__p

    @property
    def energy(self):
        """Get energy"""
        return self.__e

    @property
    def temperature(self):
        """Get temperature"""
        return self.__t
