/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2013,2014,2015,2016 by the GROMACS development team.
 * Copyright (c) 2017,2018,2019,2020,2021, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "residuetypes.h"

#include <cassert>
#include <cstdio>

#include <algorithm>
#include <optional>
#include <string>
#include <unordered_map>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strdb.h"

//! Definition for residue type that is not known.
const ResidueType c_undefinedResidueType = "Other";

//! Function object for comparisons used in std::unordered_map
class EqualCaseInsensitive
{
public:
    bool operator()(const ResidueName& lhs, const ResidueName& rhs) const
    {
        return gmx::equalCaseInsensitive(lhs, rhs);
    }
};

//! Implementation detail for ResidueTypeMap
class ResidueTypeMap::Impl
{
public:
    //! Storage object for entries.
    std::unordered_map<ResidueName, ResidueType, std::hash<ResidueName>, EqualCaseInsensitive> entries;
};

ResidueTypeMap::ResidueTypeMap() : impl_(new Impl)
{
    char line[STRLEN];
    char resname[STRLEN], restype[STRLEN], dum[STRLEN];

    gmx::FilePtr db = gmx::openLibraryFile("residuetypes.dat");

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
            addResidue(resname, restype);
        }
    }
}

ResidueTypeMap::~ResidueTypeMap() = default;

bool ResidueTypeMap::nameIndexedInResidueTypeMap(const ResidueName& residueName)
{
    return impl_->entries.find(residueName) != impl_->entries.end();
}

void ResidueTypeMap::addResidue(const ResidueName& residueName, const ResidueType& residueType)
{
    if (auto [foundIt, insertionTookPlace] = impl_->entries.insert({ residueName, residueType });
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

bool ResidueTypeMap::namedResidueHasType(const ResidueName& residueName, const ResidueType& residueType)
{
    if (auto foundIt = impl_->entries.find(residueName); foundIt != impl_->entries.end())
    {
        return gmx::equalCaseInsensitive(residueType, foundIt->second);
    }
    return false;
}

ResidueType ResidueTypeMap::typeOfNamedDatabaseResidue(const ResidueName& residueName)
{
    if (auto foundIt = impl_->entries.find(residueName); foundIt != impl_->entries.end())
    {
        return foundIt->second;
    }
    return c_undefinedResidueType;
}

std::optional<ResidueType> ResidueTypeMap::optionalTypeOfNamedDatabaseResidue(const ResidueName& residueName)
{
    if (auto foundIt = impl_->entries.find(residueName); foundIt != impl_->entries.end())
    {
        return std::make_optional(foundIt->second);
    }
    return std::nullopt;
}
