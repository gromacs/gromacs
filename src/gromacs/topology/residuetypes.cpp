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

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strdb.h"

//! Definition for residue type that is not known.
const std::string undefinedResidueType = "Other";

ResidueTypes
initializeResidueTypes()
{
    char                    line[STRLEN];
    char                    resname[STRLEN], restype[STRLEN], dum[STRLEN];
    ResidueTypes            readTypes;

    gmx::FilePtr            db = gmx::openLibraryFile("residuetypes.dat");

    while (get_a_line(db.get(), line, STRLEN))
    {
        strip_comment(line);
        trim(line);
        if (line[0] != '\0')
        {
            if (sscanf(line, "%1000s %1000s %1000s", resname, restype, dum) != 2)
            {
                gmx_fatal(FARGS, "Incorrect number of columns (2 expected) for line in residuetypes.dat");
            }
            addResidue(&readTypes, resname, restype);
        }
    }

    return readTypes;
}

int isResidueInResidueTypes(const ResidueTypes *rt, const char *resname)
{

    int    rc = -1;
    size_t i;
    for (i = 0; i < rt->size() && rc; i++)
    {
        rc = gmx_strcasecmp(rt->at(i).resname.c_str(), resname);
    }
    if (rc == 0)
    {
        return i - 1;
    }
    else
    {
        return rc;
    }
}

std::string previouslyDefinedType(const ResidueTypes *rt, const char *resname)
{
    int entry = isResidueInResidueTypes(rt, resname);
    if (entry != -1)
    {
        return rt->at(entry).restype;
    }
    else
    {
        return undefinedResidueType;
    }
}


void addResidue(ResidueTypes *rt, const char *resname, const char *restype)
{
    int found = isResidueInResidueTypes(rt, resname);

    if (found != -1)
    {
        std::string oldentry = previouslyDefinedType(rt, restype);
        if (gmx_strcasecmp(oldentry.c_str(), restype))
        {
            fprintf(stderr, "Warning: Residue '%s' already present with type '%s' in database, ignoring new type '%s'.\n",
                    resname, oldentry.c_str(), restype);
        }
    }
    else
    {
        ResidueType newType;
        newType.resname = resname;
        newType.restype = restype;
        rt->push_back(newType);
    }
}

bool isResidueTypeProtein(const ResidueTypes *rt, const char *resnm)
{

    std::string type = previouslyDefinedType(rt, resnm);
    bool        rc   = gmx_strcasecmp(type.c_str(), "Protein") == 0;
    return rc;
}

bool isResidueTypeDNA(const ResidueTypes *rt, const char *resnm)
{
    std::string type = previouslyDefinedType(rt, resnm);
    bool        rc   = gmx_strcasecmp(type.c_str(), "DNA") == 0;

    return rc;
}

bool isResidueTypeRNA(const ResidueTypes *rt, const char *resnm)
{
    std::string type = previouslyDefinedType(rt, resnm);
    bool        rc   = gmx_strcasecmp(type.c_str(), "RNA") == 0;
    return rc;
}

int findResidueIndex(ConstResidueTypesRef rt, const char *resnm)
{

    int rc = -1;
    int i  = 0;
    for (; i < rt.size() && rc; i++)
    {
        rc = gmx_strcasecmp(rt[i].resname.c_str(), resnm);
    }

    return (0 == rc) ? i-1 : -1;
}

const char *
findResidueName(ConstResidueTypesRef rt, int index)
{
    if (index >= 0 && index < rt.size())
    {
        return rt[index].resname.c_str();
    }
    else
    {
        return nullptr;
    }
}
