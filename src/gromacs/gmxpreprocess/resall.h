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
/*! \file
 * \internal \brief
 * Methods to get residue information during preprocessing.
 */
#ifndef GMX_GMXPREPROCESS_RESALL_H
#define GMX_GMXPREPROCESS_RESALL_H

#include <cstdio>

#include <filesystem>
#include <string>
#include <vector>

#include "gromacs/utility/arrayref.h"

class PreprocessingAtomTypes;

namespace gmx
{
class MDLogger;
}
struct PreprocessResidue;
struct t_symtab;

/*! \brief
 * Search for an entry in the rtp database.
 *
 * A mismatch of one character is allowed,
 * if there is only one nearly matching entry in the database,
 * a warning will be generated.
 *
 * \param[in] key The atomname to search for.
 * \param[in] rtpDBEntry Database with residue information.
 * \param[in] logger Logging object.
 * \returns The rtp residue name.
 */
std::string searchResidueDatabase(const std::string&                     key,
                                  gmx::ArrayRef<const PreprocessResidue> rtpDBEntry,
                                  const gmx::MDLogger&                   logger);

/*! \brief
 * Returns matching entry in database.
 *
 * \param[in] rtpname Name of the entry looked for.
 * \param[in] rtpDBEntry Database to search.
 * \throws If the name can not be found in the database.
 */
gmx::ArrayRef<const PreprocessResidue>::const_iterator
getDatabaseEntry(const std::string& rtpname, gmx::ArrayRef<const PreprocessResidue> rtpDBEntry);

/*! \brief
 * Read atom types into database.
 *
 * \param[in] ffdir Force field directory.
 * \returns Atom type database.
 */
PreprocessingAtomTypes read_atype(const std::filesystem::path& ffdir);

/*! \brief
 * Read in database, append to exisiting.
 *
 * \param[in] resdb Name of database file.
 * \param[inout] rtpDBEntry Database to populate.
 * \param[inout] atype Atomtype information.
 * \param[inout] tab Symbol table for names.
 * \param[in] logger MDLogger interface.
 * \param[in] bAllowOverrideRTP If entries can be overwritten in the database.
 */
void readResidueDatabase(const std::filesystem::path&    resdb,
                         std::vector<PreprocessResidue>* rtpDBEntry,
                         PreprocessingAtomTypes*         atype,
                         t_symtab*                       tab,
                         const gmx::MDLogger&            logger,
                         bool                            bAllowOverrideRTP);

/*! \brief
 * Print out database.
 *
 * \param[in] out File to write to.
 * \param[in] rtpDBEntry Database to write out.
 * \param[in] atype Atom type information.
 */
void print_resall(FILE*                                  out,
                  gmx::ArrayRef<const PreprocessResidue> rtpDBEntry,
                  const PreprocessingAtomTypes&          atype);

#endif
