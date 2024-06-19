/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
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
#include "gmxpre.h"

#include "gromacs/topology/residuetypes.h"

#include <cstdio>

#include <filesystem>
#include <string>
#include <unordered_map>
#include <utility>

#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/fileptr.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/strdb.h"
#include "gromacs/utility/stringutil.h"

//! Definition for residue type that is not known.
const ResidueType c_undefinedResidueType = "Other";

void addResidue(ResidueTypeMap* residueTypeMap, const ResidueName& residueName, const ResidueType& residueType)
{
    if (auto [foundIt, insertionTookPlace] = residueTypeMap->insert({ residueName, residueType });
        !insertionTookPlace)
    {
        if (!gmx::equalCaseInsensitive(foundIt->second, residueType))
        {
            fprintf(stderr,
                    "Warning: Residue '%s' already present with type '%s' in database, ignoring "
                    "new type '%s'.\n",
                    residueName.c_str(),
                    foundIt->second.c_str(),
                    residueType.c_str());
        }
    }
}

ResidueTypeMap residueTypeMapFromLibraryFile(const std::string& residueTypesDatFilename)
{
    char line[STRLEN];
    char resname[STRLEN], restype[STRLEN], dum[STRLEN];

    gmx::FilePtr db = gmx::openLibraryFile(residueTypesDatFilename);

    ResidueTypeMap residueTypeMap;
    while (get_a_line(db.get(), line, STRLEN))
    {
        strip_comment(line);
        trim(line);
        if (line[0] != '\0')
        {
            if (sscanf(line, "%1000s %1000s %1000s", resname, restype, dum) != 2)
            {
                gmx_fatal(
                        FARGS,
                        "Incorrect number of columns (2 expected) for line in residuetypes.dat  ");
            }
            addResidue(&residueTypeMap, resname, restype);
        }
    }
    return residueTypeMap;
}

bool namedResidueHasType(const ResidueTypeMap& residueTypeMap,
                         const ResidueName&    residueName,
                         const ResidueType&    residueType)
{
    if (auto foundIt = residueTypeMap.find(residueName); foundIt != residueTypeMap.end())
    {
        return gmx::equalCaseInsensitive(residueType, foundIt->second);
    }
    return false;
}

ResidueType typeOfNamedDatabaseResidue(const ResidueTypeMap& residueTypeMap, const ResidueName& residueName)
{
    if (auto foundIt = residueTypeMap.find(residueName); foundIt != residueTypeMap.end())
    {
        return foundIt->second;
    }
    return c_undefinedResidueType;
}
