/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2013,2014, by the GROMACS development team, led by
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

#include "gromacs/fileio/strdb.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

const char gmx_residuetype_undefined[] = "Other";

struct gmx_residuetype_t
{
    int      n;
    char **  resname;
    char **  restype;
};

int
gmx_residuetype_init(gmx_residuetype_t **prt)
{
    FILE                 *  db;
    char                    line[STRLEN];
    char                    resname[STRLEN], restype[STRLEN], dum[STRLEN];
    gmx_residuetype_t      *rt;

    snew(rt, 1);
    *prt = rt;

    rt->n        = 0;
    rt->resname  = NULL;
    rt->restype  = NULL;

    db = libopen("residuetypes.dat");

    while (get_a_line(db, line, STRLEN))
    {
        strip_comment(line);
        trim(line);
        if (line[0] != '\0')
        {
            if (sscanf(line, "%1000s %1000s %1000s", resname, restype, dum) != 2)
            {
                gmx_fatal(FARGS, "Incorrect number of columns (2 expected) for line in residuetypes.dat");
            }
            gmx_residuetype_add(rt, resname, restype);
        }
    }

    fclose(db);

    return 0;
}

int
gmx_residuetype_destroy(gmx_residuetype_t *rt)
{
    int i;

    for (i = 0; i < rt->n; i++)
    {
        sfree(rt->resname[i]);
        sfree(rt->restype[i]);
    }
    sfree(rt->resname);
    sfree(rt->restype);
    sfree(rt);

    return 0;
}

/* Return 0 if the name was found, otherwise -1.
 * p_restype is set to a pointer to the type name, or 'Other' if we did not find it.
 */
int
gmx_residuetype_get_type(gmx_residuetype_t *rt, const char * resname, const char ** p_restype)
{
    int    i, rc;

    rc = -1;
    for (i = 0; i < rt->n && rc; i++)
    {
        rc = gmx_strcasecmp(rt->resname[i], resname);
    }

    *p_restype = (rc == 0) ? rt->restype[i-1] : gmx_residuetype_undefined;

    return rc;
}

int
gmx_residuetype_add(gmx_residuetype_t *rt, const char *newresname, const char *newrestype)
{
    int           found;
    const char *  p_oldtype;

    found = !gmx_residuetype_get_type(rt, newresname, &p_oldtype);

    if (found && gmx_strcasecmp(p_oldtype, newrestype))
    {
        fprintf(stderr, "Warning: Residue '%s' already present with type '%s' in database, ignoring new type '%s'.",
                newresname, p_oldtype, newrestype);
    }

    if (found == 0)
    {
        srenew(rt->resname, rt->n+1);
        srenew(rt->restype, rt->n+1);
        rt->resname[rt->n] = gmx_strdup(newresname);
        rt->restype[rt->n] = gmx_strdup(newrestype);
        rt->n++;
    }

    return 0;
}

int
gmx_residuetype_get_alltypes(gmx_residuetype_t   *rt,
                             const char ***       p_typenames,
                             int *                ntypes)
{
    int            i, n;
    const char **  my_typename;

    n           = 0;
    my_typename = NULL;
    for (i = 0; i < rt->n; i++)
    {
        const char *const p      = rt->restype[i];
        bool              bFound = false;
        for (int j = 0; j < n && !bFound; j++)
        {
            assert(my_typename != NULL);
            bFound = !gmx_strcasecmp(p, my_typename[j]);
        }
        if (!bFound)
        {
            srenew(my_typename, n+1);
            my_typename[n] = p;
            n++;
        }
    }
    *ntypes      = n;
    *p_typenames = my_typename;

    return 0;
}

gmx_bool
gmx_residuetype_is_protein(gmx_residuetype_t *rt, const char *resnm)
{
    gmx_bool    rc;
    const char *p_type;

    if (gmx_residuetype_get_type(rt, resnm, &p_type) == 0 &&
        gmx_strcasecmp(p_type, "Protein") == 0)
    {
        rc = TRUE;
    }
    else
    {
        rc = FALSE;
    }
    return rc;
}

gmx_bool
gmx_residuetype_is_dna(gmx_residuetype_t *rt, const char *resnm)
{
    gmx_bool    rc;
    const char *p_type;

    if (gmx_residuetype_get_type(rt, resnm, &p_type) == 0 &&
        gmx_strcasecmp(p_type, "DNA") == 0)
    {
        rc = TRUE;
    }
    else
    {
        rc = FALSE;
    }
    return rc;
}

gmx_bool
gmx_residuetype_is_rna(gmx_residuetype_t *rt, const char *resnm)
{
    gmx_bool    rc;
    const char *p_type;

    if (gmx_residuetype_get_type(rt, resnm, &p_type) == 0 &&
        gmx_strcasecmp(p_type, "RNA") == 0)
    {
        rc = TRUE;
    }
    else
    {
        rc = FALSE;
    }
    return rc;
}

/* Return the size of the arrays */
int
gmx_residuetype_get_size(gmx_residuetype_t *rt)
{
    return rt->n;
}

/* Search for a residuetype with name resnm within the
 * gmx_residuetype database. Return the index if found,
 * otherwise -1.
 */
int
gmx_residuetype_get_index(gmx_residuetype_t *rt, const char *resnm)
{
    int i, rc;

    rc = -1;
    for (i = 0; i < rt->n && rc; i++)
    {
        rc = gmx_strcasecmp(rt->resname[i], resnm);
    }

    return (0 == rc) ? i-1 : -1;
}

/* Return the name of the residuetype with the given index, or
 * NULL if not found. */
const char *
gmx_residuetype_get_name(gmx_residuetype_t *rt, int index)
{
    if (index >= 0 && index < rt->n)
    {
        return rt->resname[index];
    }
    else
    {
        return NULL;
    }
}
