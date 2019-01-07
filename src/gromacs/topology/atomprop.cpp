/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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

#include "atomprop.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>

#include "gromacs/compat/make_unique.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/topology/residuetypes.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strdb.h"

/* NOTFOUND should be smallest, others larger in increasing priority */
enum {
    NOTFOUND = -4, WILDCARD, WILDPROT, PROTEIN
};

/* return number of matching characters,
   or NOTFOUND if not at least all characters in char *database match */
static int dbcmp_len(const char *search, const char *database)
{
    int i;

    i = 0;
    while (search[i] && database[i] && (search[i] == database[i]) )
    {
        i++;
    }

    if (database[i])
    {
        i = NOTFOUND;
    }
    return i;
}

static int get_prop_index(AtomProperty *ap, gmx::ArrayRef<const ResidueType> restype,
                          char *resnm, char *atomnm,
                          gmx_bool *bExact)
{
    int      j = NOTFOUND;
    long int alen, rlen;
    long int malen, mrlen;
    gmx_bool bProtein, bProtWild;

    bProtein  = isResidueTypeProtein(restype, resnm);
    bProtWild = (strcmp(resnm, "AAA") == 0);
    malen     = NOTFOUND;
    mrlen     = NOTFOUND;
    size_t i = 0;
    for (; (i < ap->entry.size()); i++)
    {
        rlen = dbcmp_len(resnm, ap->entry[i].resnm.c_str());
        if (rlen == NOTFOUND)
        {
            if ( (strcmp(ap->entry[i].resnm.c_str(), "*") == 0) ||
                 (strcmp(ap->entry[i].resnm.c_str(), "???") == 0) )
            {
                rlen = WILDCARD;
            }
            else if (strcmp(ap->entry[i].resnm.c_str(), "AAA") == 0)
            {
                rlen = WILDPROT;
            }
        }
        alen = dbcmp_len(atomnm, ap->entry[i].atomnm.c_str());
        if ( (alen > NOTFOUND) && (rlen > NOTFOUND))
        {
            if ( ( (alen > malen) && (rlen >= mrlen)) ||
                 ( (rlen > mrlen) && (alen >= malen) ) )
            {
                malen = alen;
                mrlen = rlen;
                j     = i;
            }
        }
    }

    *bExact = ((malen == static_cast<long int>(strlen(atomnm))) &&
               ((mrlen == static_cast<long int>(strlen(resnm))) ||
                ((mrlen == WILDPROT) && bProtWild) ||
                ((mrlen == WILDCARD) && !bProtein && !bProtWild)));

    if (debug)
    {
        fprintf(debug, "searching residue: %4s atom: %4s\n", resnm, atomnm);
        if (j == NOTFOUND)
        {
            fprintf(debug, " not successful\n");
        }
        else
        {
            fprintf(debug, " match: %4s %4s\n", ap->entry[i].resnm.c_str(), ap->entry[i].atomnm.c_str());
        }
    }
    return j;
}

static void add_prop(AtomProperty *ap, gmx::ArrayRef<const ResidueType> restype,
                     char *resnm, char *atomnm,
                     real p, int line)
{
    bool bExact;

    int  j = get_prop_index(ap, restype, resnm, atomnm, &bExact);

    if (!bExact)
    {
        BaseEntry entry;
        entry.atomnm = atomnm;
        entry.resnm  = resnm;
        ap->entry.push_back(entry);
        j                     = ap->entry.size() - 1;
    }
    if (ap->entry[j].isAvailable)
    {
        if (ap->entry[j].value == p)
        {
            fprintf(stderr, "Warning double identical entries for %s %s %g on line %d in file %s\n",
                    resnm, atomnm, p, line, ap->db.c_str());
        }
        else
        {
            fprintf(stderr, "Warning double different entries %s %s %g and %g on line %d in file %s\n"
                    "Using last entry (%g)\n",
                    resnm, atomnm, p, ap->entry[j].value, line, ap->db.c_str(), p);
            ap->entry[j].value = p;
        }
    }
    else
    {
        ap->entry[j].isAvailable = TRUE;
        ap->entry[j].value       = p;
    }
}

static void read_prop(AtomProperties *aps, int eprop, double factor)
{
    char          line[STRLEN], resnm[32], atomnm[32];
    double        pp;
    int           line_no;
    AtomProperty *ap = &aps->prop[eprop];

    gmx::FilePtr  fp = gmx::openLibraryFile(ap->db);
    line_no = 0;
    while (get_a_line(fp.get(), line, STRLEN))
    {
        line_no++;
        if (sscanf(line, "%31s %31s %20lf", resnm, atomnm, &pp) == 3)
        {
            pp *= factor;
            add_prop(ap, aps->restype, resnm, atomnm, pp, line_no);
        }
        else
        {
            fprintf(stderr, "WARNING: Error in file %s at line %d ignored\n",
                    ap->db.c_str(), line_no);
        }
    }
    ap->isSet = TRUE;
}

static void set_prop(AtomProperties *aps, int eprop)
{
    const char   *fns[epropNR]  = { "atommass.dat", "vdwradii.dat", "dgsolv.dat", "electroneg.dat", "elements.dat" };
    double        fac[epropNR]  = { 1.0,    1.0,  418.4, 1.0, 1.0 };
    double        def[epropNR]  = { 12.011, 0.14, 0.0, 2.2, -1 };
    AtomProperty *ap            = &aps->prop[eprop];
    if (!ap->isSet)
    {
        ap->db  = fns[eprop];
        ap->def = def[eprop];
        read_prop(aps, eprop, fac[eprop]);

        if (debug)
        {
            fprintf(debug, "Entries in %s: %zu\n", ap->db.c_str(), ap->entry.size());
        }

        if ( ( (!aps->bWarned) && (eprop == epropMass) ) || (eprop == epropVDW))
        {
            printf("\n"
                   "WARNING: Masses and atomic (Van der Waals) radii will be guessed\n"
                   "         based on residue and atom names, since they could not be\n"
                   "         definitively assigned from the information in your input\n"
                   "         files. These guessed numbers might deviate from the mass\n"
                   "         and radius of the atom type. Please check the output\n"
                   "         files if necessary.\n\n");
            aps->bWarned = TRUE;
        }
    }
}

static void vdw_warning(FILE *fp)
{
    if (nullptr != fp)
    {
        fprintf(fp, "NOTE: From version 5.0 %s uses the Van der Waals radii\n",
                gmx::getProgramContext().displayName());
        fprintf(fp, "from the source below. This means the results may be different\n");
        fprintf(fp, "compared to previous GROMACS versions.\n");
        please_cite(fp, "Bondi1964a");
    }
}

bool setAtomProperty(AtomProperties *aps,
                     int eprop, const char *resnm, const char *atomnm,
                     real *value)
{
    int           j;
#define MAXQ 32
    char          atomname[MAXQ], resname[MAXQ];
    gmx_bool      bExact;

    set_prop(aps, eprop);
    if ((strlen(atomnm) > MAXQ-1) || (strlen(resnm) > MAXQ-1))
    {
        if (debug)
        {
            fprintf(debug, "WARNING: will only compare first %d characters\n",
                    MAXQ-1);
        }
    }
    if (isdigit(atomnm[0]))
    {
        int i;
        /* put digit after atomname */
        for (i = 1; i < MAXQ-1 && atomnm[i] != '\0'; i++)
        {
            atomname[i-1] = atomnm[i];
        }
        atomname[i-1] = atomnm[0];
        atomname[i]   = '\0';
    }
    else
    {
        strncpy(atomname, atomnm, MAXQ-1);
    }
    strncpy(resname, resnm, MAXQ-1);

    j = get_prop_index(&(aps->prop[eprop]), aps->restype, resname,
                       atomname, &bExact);

    if (eprop == epropVDW && !aps->bWarnVDW)
    {
        vdw_warning(stdout);
        aps->bWarnVDW = TRUE;
    }
    if (j >= 0)
    {
        *value = aps->prop[eprop].entry[j].value;
        return TRUE;
    }
    else
    {
        *value = aps->prop[eprop].def;
        return FALSE;
    }
}

std::string elementFromAtomnumber(AtomProperties *aps, int atomnumber)
{
    set_prop(aps, epropElement);
    for (size_t i = 0; (i < aps->prop[epropElement].entry.size()); i++)
    {
        if (std::round(aps->prop[epropElement].entry[i].value) == atomnumber)
        {
            return aps->prop[epropElement].entry[i].atomnm;
        }
    }
    return "";
}

int atomnumberFromElement(AtomProperties *aps, const char *elem)
{
    set_prop(aps, epropElement);
    for (size_t i = 0; (i < aps->prop[epropElement].entry.size()); i++)
    {
        if (gmx_strcasecmp(aps->prop[epropElement].entry[i].atomnm.c_str(), elem) == 0)
        {
            return gmx::roundToInt(aps->prop[epropElement].entry[i].value);
        }
    }
    return -1;
}
