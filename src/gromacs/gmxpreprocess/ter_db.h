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

#ifndef GMX_GMXPREPROCESS_TER_DB_H
#define GMX_GMXPREPROCESS_TER_DB_H

#include <filesystem>
#include <vector>

class PreprocessingAtomTypes;
struct MoleculePatchDatabase;

namespace gmx
{
template<typename>
class ArrayRef;
}

/*! \brief
 * Read database for N&C terminal modifications.
 *
 * \param[in] ffdir Directory for files.
 * \param[in] ter Which terminal side to read.
 * \param[inout] tbptr Database for terminii entry to populate.
 * \param[in] atype Database for atomtype information.
 * \returns Number of entries entered into database.
 */
int read_ter_db(const std::filesystem::path&        ffdir,
                char                                ter,
                std::vector<MoleculePatchDatabase>* tbptr,
                PreprocessingAtomTypes*             atype);

/*! \brief
 * Return entries for modification blocks that match a residue name.
 *
 * \param[in] tb Complete modification database.
 * \param[in] resname Residue name for terminus.
 * \returns A list of pointers to entries that match, or of nullptr for no matching entry.
 */
std::vector<MoleculePatchDatabase*> filter_ter(gmx::ArrayRef<MoleculePatchDatabase> tb, const char* resname);

/*! \brief
 * Interactively select one terminus.
 *
 * \param[in] tb List of possible entries, with pointer to actual entry or nullptr.
 * \param[in] title Name of entry.
 * \returns The modification block selected.
 */
MoleculePatchDatabase* choose_ter(gmx::ArrayRef<MoleculePatchDatabase*> tb, const char* title);

#endif
