/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2014,2018,2019, by the GROMACS development team, led by
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
#ifndef GMX_TOPOLOGY_RESIDUETYPES_H
#define GMX_TOPOLOGY_RESIDUETYPES_H

#include <string>
#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"

struct ResidueType
{
    //! Name of the stored residue.
    std::string resname;
    //! Type of the stored residue (Protein, DNA, Lipid, Other).
    std::string restype;
};

//! Convenience definition for vector of residue types.
using ResidueTypes = std::vector<ResidueType>;
//! Convenience definition for const reference to residue types.
using ConstResidueTypesRef = gmx::ArrayRef<const ResidueType>;

//! Prepare a ResidueTypes datastructure with the basic information.
ResidueTypes initializeResidueTypes();

/*! \brief
 * Find out if entry already exists in the datbase for \p resname.
 *
 * \param[in] rt Current state of the database.
 * \param[in] resname Name of new residue that should be checked.
 * \returns If the residue is already defined or not.
 */
bool isResidueInResidueTypes (const ResidueTypes *rt, const char *resname);

/*! \brief
 * Find type previously defined residue.
 *
 * \param[in] rt ResidueTypes data.
 * \param[in] resname Name of residue.
 * \returns Name of the entry.
 */
std::string previouslyDefinedType(const ResidueTypes *rt, const char *resname);

/*! \brief
 * Add new entry to ResidueTypes.
 *
 * Performs checking if entry is already defined.
 *
 * \param[in] rt Database to add entry to.
 * \param[in] resname Name of new residue.
 * \param[in] restype Type of enw residue.
 */
void addResidue(ResidueTypes *rt, const char *resname, const char *restype);

/*! \brief
 * Is the current residue of type protein.
 *
 * \param[in] rt ResidueTypes data.
 * \param[in] resnm Residue name to check.
 * \returns Result of check true or false.
 */
bool isResidueTypeProtein(const ResidueTypes *rt, const char *resnm);

/*! \brief
 * Is the current residue of type DNA.
 *
 * \param[in] rt ResidueTypes data.
 * \param[in] resnm Residue name to check.
 * \returns Result of check true or false.
 */
bool isResidueTypeDNA(const ResidueTypes *rt, const char *resnm);

/*! \brief
 * Is the current residue of type RNA.
 *
 * \param[in] rt ResidueTypes data.
 * \param[in] resnm Residue name to check.
 * \returns Result of check true or false.
 */
bool isResidueTypeRNA(const ResidueTypes *rt, const char *resnm);

/*! \brief
 * Find index for \p resnm into \p rt.
 *
 * \param[in] rt ResidueTypes data.
 * \param[in] resnm Name of residue to find index for.
 * \returns Index for the residue or -1 if not found.
 */
int findResidueIndex(ConstResidueTypesRef rt, const char *resnm);

/*! \brief
 * Find name of residue with \p index in \p rt.
 *
 * \param[in] rt ResidueTypes data.
 * \param[in] index Index to check in \p rt.
 * \returns Name of residue with \p index or nullptr if not found.
 */
const char *findResidueName(ConstResidueTypesRef rt, int index);

#endif
