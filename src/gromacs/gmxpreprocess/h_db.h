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

#ifndef GMX_GMXPREPROCESS_H_DB_H
#define GMX_GMXPREPROCESS_H_DB_H

#include <cstdio>

#include <vector>

#include "gromacs/utility/arrayref.h"

struct MoleculePatch;
struct MoleculePatchDatabase;

/* functions for the h-database */

void read_ab(char* line, const char* fn, MoleculePatch* ab);
/* Read one add block */

/*! \brief
 * Read the databse from hdb file(s).
 *
 * \param[in] ffdir Directory for files.
 * \param[inout] globalPatches The database for atom modifications to populate.
 * \returns The number of modifications stored.
 */
int read_h_db(const char* ffdir, std::vector<MoleculePatchDatabase>* globalPatches);

void print_ab(FILE* out, const MoleculePatch& ab, const char* nname);
/* print one add block */

/*! \brief
 * Search for an entry.
 *
 * \param[in] globalPatches Database to search.
 * \param[in] key Name to search for.
 */
gmx::ArrayRef<const MoleculePatchDatabase>::iterator
search_h_db(gmx::ArrayRef<const MoleculePatchDatabase> globalPatches, const char* key);
/* Search for an entry in the database */

#endif
