#!/usr/bin/env python
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2019, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.

import unittest

import gmxapi as gmx


class ImmediateResultTestCase(unittest.TestCase):
    """Test data model and data flow for basic operations."""

    def test_scalar(self):
        operation = gmx.operation.make_constant(42)
        assert isinstance(operation.dtype, type)
        assert operation.dtype == int
        assert operation.result() == 42

    def test_list(self):
        list_a = [1, 2, 3]

        # TODO: test input validation
        list_result = gmx.operation.concatenate_lists(sublists=[list_a])
        # TODO: should be NDArray
        assert list_result.dtype == type(list_a)
        # Note: this is specifically for the built-in tuple type.
        # Equality comparison may work differently for different sequence types.
        assert tuple(list_result.result()) == tuple(list_a)
        assert len(list_result.result()) == len(list_a)

        list_result = gmx.operation.concatenate_lists([list_a, list_a])
        assert len(list_result.result()) == len(list_a) * 2
        assert tuple(list_result.result()) == tuple(list_a + list_a)

        list_b = gmx.operation.make_constant([42])

        list_result = gmx.operation.concatenate_lists(sublists=[list_b])
        assert list_result.result()[0] == 42

        list_result = gmx.operation.append_list(list_a, list_b)
        assert len(list_result.result()) == len(list_a) + 1
        assert tuple(list_result.result()) == tuple(list(list_a) + [42])


if __name__ == '__main__':
    unittest.main()
