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
Data structure carrying information on the simulated system.
"""
import warnings
import numpy as np

from ..util import error as pv_error


class SystemData(object):
    r"""SystemData: Information about the atoms and molecules in the system.

    The information stored in SystemData objects describes the atom and molecules
    in the system as far as the physical validation tests need it.

    The system is described in terms of
    natoms: the total number of atoms in the system
    nconstraints: the total number of constraints in the system
    ndof_reduction_tra: global reduction of translational degrees of freedom (e.g.
                        due to constraining the center of mass of the system)
    ndof_reduction_rot: global reduction of rotational degrees of freedom (e.g.
                        due to constraining the center of mass of the system)

    The atoms are described in terms of
    mass: a list of the mass of every atom in the system

    The molecules are described by
    molecule_idx: a list with the indices first atoms of every molecule (this assumes
                  that the atoms are sorted by molecule)
    nconstraints_per_molecule: a list with the number of constraints in every molecule

    Only used internally:
    ndof_per_molecule: a list with the number of degrees of freedom of every molecule

    Reserved for future use:
    bonds
    constrained_bonds

    Notes:
    ------
    kinetic_energy.mb_ensemble() only requires information on the system
        (natoms, nconstraints, ndof_reduction_tra, ndof_reduction_rot)
    kinetic_energy.equipartition() additionally requires information on the atoms and molecules
        (mass, molecule_idx, nconstraints_per_molecule)
    All other tests do not require and information from SystemData.
    """

    def __init__(self,
                 natoms=None, nconstraints=None,
                 ndof_reduction_tra=None, ndof_reduction_rot=None,
                 mass=None, molecule_idx=None, nconstraints_per_molecule=None):
        self.__natoms = None
        self.__nconstraints = None
        self.__ndof_reduction_tra = None
        self.__ndof_reduction_rot = None
        self.__mass = None
        self.__molecule_idx = None
        self.__nconstraints_per_molecule = None
        self.__ndof_per_molecule = None
        self.__bonds = None
        self.__constrained_bonds = None

        if natoms is not None:
            self.natoms = natoms
        if nconstraints is not None:
            self.nconstraints = nconstraints
        if ndof_reduction_tra is not None:
            self.ndof_reduction_tra = ndof_reduction_tra
        if ndof_reduction_rot is not None:
            self.ndof_reduction_rot = ndof_reduction_rot
        if mass is not None:
            self.mass = mass
        if molecule_idx is not None:
            self.molecule_idx = molecule_idx
        if nconstraints_per_molecule is not None:
            self.nconstraints_per_molecule = nconstraints_per_molecule

    @property
    def natoms(self):
        """int: Number of atoms in the system

        """
        return self.__natoms

    @natoms.setter
    def natoms(self, natoms):
        self.__natoms = int(natoms)

    @property
    def nconstraints(self):
        """float: Total number of constraints in the system

        Does not include the reduction of degrees of freedom in the absence of
        external forces.

        """
        return self.__nconstraints

    @nconstraints.setter
    def nconstraints(self, nconstraints):
        self.__nconstraints = float(nconstraints)

    @property
    def ndof_reduction_tra(self):
        """float: Number of translational degrees of freedom deducted
        from 3*[# of molecules]

        """
        return self.__ndof_reduction_tra

    @ndof_reduction_tra.setter
    def ndof_reduction_tra(self, ndof_reduction_tra):
        self.__ndof_reduction_tra = float(ndof_reduction_tra)

    @property
    def ndof_reduction_rot(self):
        """float: Number of rotational degrees of freedom deducted
        from 3*[# of molecules]

        """
        return self.__ndof_reduction_rot

    @ndof_reduction_rot.setter
    def ndof_reduction_rot(self, ndof_reduction_rot):
        self.__ndof_reduction_rot = float(ndof_reduction_rot)

    @property
    def mass(self):
        """nd-array: Mass vector for the atoms

        Setter accepts array-like objects.

        """
        return self.__mass

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
        self.__mass = mass

    @property
    def molecule_idx(self):
        """nd-array: List of index of first atom of each molecule

        Setter accepts array-like objects.

        """
        return self.__molecule_idx

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
            self.__nconstraints_per_molecule = None
        self.__molecule_idx = molecule_idx

    @property
    def nconstraints_per_molecule(self):
        """nd-array: List of number of constraints per molecule

        Setter accepts array-like objects.

        """
        return self.__nconstraints_per_molecule

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

        self.__nconstraints_per_molecule = nconstraints_per_molecule

    @property
    def ndof_per_molecule(self):
        """nd-array: List of number of degrees of freedom per molecule

        Setter accepts array-like objects.

        """
        return self.__ndof_per_molecule

    @ndof_per_molecule.setter
    def ndof_per_molecule(self, ndof_per_molecule):
        # used internally - check for consistency?
        self.__ndof_per_molecule = ndof_per_molecule

    @property
    def bonds(self):
        """List[List[int]]: List of bonds per molecule
        """
        return self.__bonds

    @bonds.setter
    def bonds(self, bonds):
        self.__bonds = bonds

    @property
    def constrained_bonds(self):
        """List[List[int]]: List of constrained bonds per molecule
        """
        return self.__constrained_bonds

    @constrained_bonds.setter
    def constrained_bonds(self, constrained_bonds):
        self.__constrained_bonds = constrained_bonds
