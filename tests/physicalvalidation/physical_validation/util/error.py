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
Module containing the custom exception classes for the physical_validation
package.
"""


class PhysicalValidationError(Exception):
    r"""Base class for exceptions in the physical_validation module."""
    pass


class InputError(PhysicalValidationError):
    r"""Exception raised for input errors"""

    def __init__(self, argument, message):
        r"""

        Parameters
        ----------
        argument : string or list of strings
        message : string
        """
        self.argument = argument
        self.message = message


class ParserValueNotSetError(PhysicalValidationError):
    r"""
    Exception raised if a requested data value
    was not set by the user previously
    """

    def __init__(self, message):
        r"""

        Parameters
        ----------
        message : string
        """
        self.message = message


class FileFormatError(PhysicalValidationError):
    r"""Exception raised for files not following expected format"""

    def __init__(self, argument, message):
        r"""

        Parameters
        ----------
        argument : string or list of strings
        message : string
        """
        self.argument = argument
        self.message = message
