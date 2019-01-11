/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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

#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strdb.h"

//! Definition for residue type that is not known.
const char undefinedResidueType[] = "Other";

//! Single ResidueType entry object
struct ResidueTypeEntry
{
    //! Name of residue.
    char *residueName;
    //! Type name for residue entry.
    char *residueType;
};

//! Implementation detail for ResidueTypes
class ResidueType::Impl
{
    public:
        //! Number of entries
        int               n     = 0;
        //! Storage object for entries.
        ResidueTypeEntry *entry = nullptr;
};

ResidueType::ResidueType()
    : impl_(new Impl)
{
    char                    line[STRLEN];
    char                    resname[STRLEN], restype[STRLEN], dum[STRLEN];

    gmx::FilePtr            db = gmx::openLibraryFile("residuetypes.dat");

    while (get_a_line(db.get(), line, STRLEN))
    {
        strip_comment(line);
        trim(line);
        if (line[0] != '\0')
        {
            if (sscanf(line, "%1000s %1000s %1000s", resname, restype, dum) != 2)
            {
                gmx_fatal(FARGS, "Incorrect number of columns (2 expected) for line in residuetypes.dat  ");
            }
            addResidue(resname, restype);
        }
    }
}

ResidueType::~ResidueType()
{
    for (int i = 0; i < impl_->n; i++)
    {
        sfree(impl_->entry[i].residueName);
        sfree(impl_->entry[i].residueType);
    }
    sfree(impl_->entry);
}

bool ResidueType::nameIndexedInResidueTypes(const char *residueName, const char **residueTypePointer)
{
    int    i, rc;

    rc = -1;
    for (i = 0; i < impl_->n && rc; i++)
    {
        rc = gmx_strcasecmp(impl_->entry[i].residueName, residueName);
    }

    *residueTypePointer = (rc == 0) ? impl_->entry[i-1].residueType : undefinedResidueType;

    return (rc == 0);
}

void ResidueType::addResidue(const char *residueName, const char *residueType)
{
    const char *  p_oldtype;

    bool          found = nameIndexedInResidueTypes(residueName, &p_oldtype);

    if (found && gmx_strcasecmp(p_oldtype, residueType))
    {
        fprintf(stderr, "Warning: Residue '%s' already present with type '%s' in database, ignoring new type '%s'.",
                residueName, p_oldtype, residueType);
    }

    if (!found)
    {
        srenew(impl_->entry, impl_->n+1);
        impl_->entry[impl_->n].residueName = gmx_strdup(residueName);
        impl_->entry[impl_->n].residueType = gmx_strdup(residueType);
        impl_->n++;
    }
}

bool ResidueType::namedResidueHasType(const char *residueName, const char *residueType)
{
    bool        rc;
    const char *p_type;

    rc = nameIndexedInResidueTypes(residueName, &p_type) &&
        gmx_strcasecmp(p_type, residueType) == 0;

    return rc;
}

int ResidueType::numberOfEntries() const
{
    return impl_->n;
}

int ResidueType::indexFromResidueName(const char *residueName) const
{
    int i, rc;

    rc = -1;
    for (i = 0; i < impl_->n && rc; i++)
    {
        rc = gmx_strcasecmp(impl_->entry[i].residueName, residueName);
    }

    return (0 == rc) ? i-1 : -1;
}

const char *ResidueType::nameFromResidueIndex(int index) const
{
    if (index >= 0 && index < impl_->n)
    {
        return impl_->entry[index].residueName;
    }
    else
    {
        return nullptr;
    }
}
