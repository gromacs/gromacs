/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2013,2014,2015,2016 by the GROMACS development team.
 * Copyright (c) 2017,2018,2019,2020, by the GROMACS development team, led by
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
#include <iterator>
#include <optional>
#include <string>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strdb.h"

//! Definition for residue type that is not known.
const std::string c_undefinedResidueType = "Other";

//! Single ResidueType entry object
struct ResidueTypeEntry
{
    //! Default constructor creates complete object.
    ResidueTypeEntry(const std::string& rName, const std::string& rType) :
        residueName(rName),
        residueType(rType)
    {
    }
    //! Name of the residue in the entry.
    std::string residueName;
    //! Type of the residue in the entry.
    std::string residueType;
};

//! Implementation detail for ResidueTypes
class ResidueType::Impl
{
public:
    //! Storage object for entries.
    std::vector<ResidueTypeEntry> entry;
};

ResidueType::ResidueType() : impl_(new Impl)
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

ResidueType::~ResidueType() {}

/*! \brief
 * Return an optional const iterator to a residue entry that matches the given name.
 *
 * \param[in] entries Currently registered residue entries in the database.
 * \param[in] residueName Name of a residue to compare to database.
 * \returns An optional iterator to the residue entry that was found.
 */
static std::optional<gmx::ArrayRef<const ResidueTypeEntry>::const_iterator>
findResidueEntryWithName(gmx::ArrayRef<const ResidueTypeEntry> entries, const std::string& residueName)
{
    auto foundIt =
            std::find_if(entries.begin(), entries.end(), [&residueName](const ResidueTypeEntry& old) {
                return gmx::equalCaseInsensitive(residueName, old.residueName);
            });
    return (foundIt != entries.end()) ? std::make_optional(foundIt) : std::nullopt;
}

bool ResidueType::nameIndexedInResidueTypes(const std::string& residueName)
{
    return findResidueEntryWithName(impl_->entry, residueName).has_value();
}

void ResidueType::addResidue(const std::string& residueName, const std::string& residueType)
{
    if (auto foundIt = findResidueEntryWithName(impl_->entry, residueName))
    {
        if (!gmx::equalCaseInsensitive((*foundIt)->residueType, residueType))
        {
            fprintf(stderr,
                    "Warning: Residue '%s' already present with type '%s' in database, ignoring "
                    "new type '%s'.\n",
                    residueName.c_str(), (*foundIt)->residueType.c_str(), residueType.c_str());
        }
    }
    else
    {
        impl_->entry.emplace_back(residueName, residueType);
    }
}

bool ResidueType::namedResidueHasType(const std::string& residueName, const std::string& residueType)
{
    auto foundIt = findResidueEntryWithName(impl_->entry, residueName);
    return foundIt ? gmx::equalCaseInsensitive(residueType, (*foundIt)->residueType) : false;
}

int ResidueType::numberOfEntries() const
{
    return impl_->entry.size();
}

int ResidueType::indexFromResidueName(const std::string& residueName) const
{
    gmx::ArrayRef<const ResidueTypeEntry> temp(impl_->entry);
    auto                                  foundIt = findResidueEntryWithName(temp, residueName);
    return foundIt ? std::distance(temp.begin(), *foundIt) : -1;
}

std::string ResidueType::nameFromResidueIndex(int index) const
{
    if (index >= 0 && index < gmx::ssize(impl_->entry))
    {
        return impl_->entry[index].residueName;
    }
    else
    {
        return "";
    }
}

std::string ResidueType::typeOfNamedDatabaseResidue(const std::string& residueName)
{
    auto foundIt = findResidueEntryWithName(impl_->entry, residueName);
    return foundIt ? (*foundIt)->residueType : c_undefinedResidueType;
}

std::optional<std::string> ResidueType::optionalTypeOfNamedDatabaseResidue(const std::string& residueName)
{
    auto foundIt = findResidueEntryWithName(impl_->entry, residueName);
    return foundIt ? std::make_optional((*foundIt)->residueType) : std::nullopt;
}
