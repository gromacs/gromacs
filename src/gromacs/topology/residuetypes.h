/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2014,2018,2019,2020,2021, by the GROMACS development team, led by
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
#include <unordered_map>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/stringutil.h"

/*! \brief Convenience type aliases
 *
 * These are not as useful as strong types, but they will
 * help clarify usage to humans in some cases. */
//! \{
using ResidueName = std::string;
using ResidueType = std::string;
//! \}

/*! \brief Maps residue names to residue types
 *
 * The contents are typically loaded from share/top/residuetypes.dat
 * or similar file provided in the users's working directory.
 */
using ResidueTypeMap =
        std::unordered_map<ResidueName, ResidueType, std::hash<ResidueName>, gmx::EqualCaseInsensitive>;

/*! \brief
 * Add entry to ResidueTypeMap if unique.
 *
 * \param[in] residueTypeMap Map to which to add new name+type entry
 * \param[in] residueName    Name of new residue.
 * \param[in] residueType    Type of new residue.
 */
void addResidue(ResidueTypeMap* residueTypeMap, const ResidueName& residueName, const ResidueType& residueType);

/*! \brief Returns a ResidueTypeMap filled from a file
 *
 * The value of the parameter is typically "residuetypes.dat" which
 * treats that as a GROMACS library file, ie. loads it from the working
 * directory or from "share/top" corresponding to the sourced GMXRC.
 *
 * \param[in] residueTypesDatFilename Library file to read and from which to fill the returned map
 */
ResidueTypeMap residueTypeMapFromLibraryFile(const std::string& residueTypesDatFilename);

/*! \brief
 * Checks if the indicated \p residueName is of \p residueType.
 *
 * \param[in] residueTypeMap Map to search
 * \param[in] residueName    Residue that should be checked.
 * \param[in] residueType    Which ResidueType the residue should have.
 * \returns If the check was successful.
 */
bool namedResidueHasType(const ResidueTypeMap& residueTypeMap,
                         const ResidueName&    residueName,
                         const ResidueType&    residueType);

/*! \brief
 * Return the residue type if a residue with that name exists, or "Other"
 *
 * \param[in] residueTypeMap Map to search
 * \param[in] residueName    Name of the residue to search for.
 * \returns The residue type of any matching residue, or "Other"
 */
ResidueType typeOfNamedDatabaseResidue(const ResidueTypeMap& residueTypeMap, const ResidueName& residueName);

#endif
