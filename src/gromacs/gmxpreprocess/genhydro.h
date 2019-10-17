/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2018,2019, by the GROMACS development team, led by
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

#ifndef GMX_GMXPREPROCESS_GENHYDRO_H
#define GMX_GMXPREPROCESS_GENHYDRO_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"

struct t_atoms;
struct t_symtab;
struct MoleculePatchDatabase;

/*! \brief
 * Generate hydrogen atoms and N and C terminal patches.
 *
 * \param[inout] initialAtoms The input atoms data structure to be modified.
 * \param[inout] localAtoms The extra atoms for reassigning the new entries.
 * \param[inout] xptr Coordinates to be updated with those for new atoms.
 * \param[in] globalPatches The atom modifications to use.
 * \param[inout] symtab Global symbol table for atom names.
 * \param[in] nterpairs Number of termini pairs in the molecule.
 * \param[in] ntdb Entries for N-terminus in each chain, each entry can be valid or nullptr.
 * \param[in] ctdb Entries for C-terminus in each cahin, each entry can be valid or nullptr.
 * \param[in] rN Residue number of the N-terminus of each chain.
 * \param[in] rC Residue number of the C-terminus of each chain.
 * \param[in] bMissing If routine should continue if atoms are not found.
 * \returns New total number of atoms.
 */
int add_h(t_atoms**                                   initialAtoms,
          t_atoms**                                   localAtoms,
          std::vector<gmx::RVec>*                     xptr,
          gmx::ArrayRef<const MoleculePatchDatabase>  globalPatches,
          t_symtab*                                   symtab,
          int                                         nterpairs,
          gmx::ArrayRef<MoleculePatchDatabase* const> ntdb,
          gmx::ArrayRef<MoleculePatchDatabase* const> ctdb,
          gmx::ArrayRef<const int>                    rN,
          gmx::ArrayRef<const int>                    rC,
          bool                                        bMissing);
#endif
