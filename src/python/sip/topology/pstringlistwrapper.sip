/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

// This class is needed to transparently map 'atomname' and similar fields of struct t_atoms.
// Avoids copying and mapping entire array to python.

%ModuleHeaderCode
class PStringListWrapper {
public:
    PStringListWrapper(char ***stringList, size_t size);
    const char* operator[] (size_t i);
private:
    char ***stringList;
    size_t size;
};
%End

%ModuleCode
PStringListWrapper::PStringListWrapper(char ***stringList, size_t size) :
    stringList(stringList), size(size)
    {}
const char* PStringListWrapper::operator[] (size_t i) {
    if (i >= size)
        return NULL;
    return *stringList[i];
}
%End

class PStringListWrapper /NoDefaultCtors/ {
public:
    const char* operator[] (int i);
    %MethodCode
        sipRes = (*sipCpp)[a0];

        if (!sipRes) {
            SIP_BLOCK_THREADS;
            PyErr_SetString(PyExc_IndexError, "Index out of range");
            SIP_UNBLOCK_THREADS;
            sipIsErr = 1;
        }
    %End
};
