/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

#ifndef GMX_GMXPREPROCESS_SPECBOND_H
#define GMX_GMXPREPROCESS_SPECBOND_H

#include <string>
#include <vector>

#include "gromacs/topology/symtab.h"
#include "gromacs/utility/vectypes.h"

struct t_atoms;
struct SpecialBond;

struct DisulfideBond
{
    int         firstResidue = -1, secondResidue = -1;
    std::string firstAtom, secondAtom;
    /* Vector used to read custom improper dihedral defintion
     * from specbond.dat. Each entry will a string which has
     * the format [residue]-[atom name] (e.g. B-SG). The letter
     * can be A or B, corresponding to residue A or B in
     * the specbond.dat file.
     */
    std::vector<std::string> customImproper;
};

std::vector<DisulfideBond>
makeDisulfideBonds(t_atoms* pdba, t_symtab* symtab, rvec x[], bool bInteractive, bool bVerbose);

std::vector<SpecialBond> generateSpecialBonds();

#endif
