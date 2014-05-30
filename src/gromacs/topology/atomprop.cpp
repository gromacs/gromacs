/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "gromacs/fileio/strdb.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/math/utilities.h"
#include "gromacs/topology/residuetypes.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

typedef struct {
    gmx_bool    bSet;
    int         nprop, maxprop;
    char       *db;
    double      def;
    char      **atomnm;
    char      **resnm;
    gmx_bool   *bAvail;
    real       *value;
} aprop_t;

typedef struct gmx_atomprop {
    gmx_bool           bWarned, bWarnVDW;
    aprop_t            prop[epropNR];
    gmx_residuetype_t *restype;
} gmx_atomprop;



/* NOTFOUND should be smallest, others larger in increasing priority */
enum {
    NOTFOUND = -4, WILDCARD, WILDPROT, PROTEIN
};

/* return number of matching characters,
   or NOTFOUND if not at least all characters in char *database match */
static int dbcmp_len(char *search, char *database)
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

static int get_prop_index(aprop_t *ap, gmx_residuetype_t *restype,
                          char *resnm, char *atomnm,
                          gmx_bool *bExact)
{
    int      i, j = NOTFOUND;
    long int alen, rlen;
    long int malen, mrlen;
    gmx_bool bProtein, bProtWild;

    bProtein  = gmx_residuetype_is_protein(restype, resnm);
    bProtWild = (strcmp(resnm, "AAA") == 0);
    malen     = NOTFOUND;
    mrlen     = NOTFOUND;
    for (i = 0; (i < ap->nprop); i++)
    {
        rlen = dbcmp_len(resnm, ap->resnm[i]);
        if (rlen == NOTFOUND)
        {
            if ( (strcmp(ap->resnm[i], "*") == 0) ||
                 (strcmp(ap->resnm[i], "???") == 0) )
            {
                rlen = WILDCARD;
            }
            else if (strcmp(ap->resnm[i], "AAA") == 0)
            {
                rlen = WILDPROT;
            }
        }
        alen = dbcmp_len(atomnm, ap->atomnm[i]);
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

    *bExact = ((malen == (long int)strlen(atomnm)) &&
               ((mrlen == (long int)strlen(resnm)) ||
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
            fprintf(debug, " match: %4s %4s\n", ap->resnm[j], ap->atomnm[j]);
        }
    }
    return j;
}

static void add_prop(aprop_t *ap, gmx_residuetype_t *restype,
                     char *resnm, char *atomnm,
                     real p, int line)
{
    int      i, j;
    gmx_bool bExact;

    j = get_prop_index(ap, restype, resnm, atomnm, &bExact);

    if (!bExact)
    {
        if (ap->nprop >= ap->maxprop)
        {
            ap->maxprop += 10;
            srenew(ap->resnm, ap->maxprop);
            srenew(ap->atomnm, ap->maxprop);
            srenew(ap->value, ap->maxprop);
            srenew(ap->bAvail, ap->maxprop);
            for (i = ap->nprop; (i < ap->maxprop); i++)
            {
                ap->atomnm[i] = NULL;
                ap->resnm[i]  = NULL;
                ap->value[i]  = 0;
                ap->bAvail[i] = FALSE;
            }
        }
        ap->atomnm[ap->nprop] = gmx_strdup(atomnm);
        ap->resnm[ap->nprop]  = gmx_strdup(resnm);
        j                     = ap->nprop;
        ap->nprop++;
    }
    if (ap->bAvail[j])
    {
        if (ap->value[j] == p)
        {
            fprintf(stderr, "Warning double identical entries for %s %s %g on line %d in file %s\n",
                    resnm, atomnm, p, line, ap->db);
        }
        else
        {
            fprintf(stderr, "Warning double different entries %s %s %g and %g on line %d in file %s\n"
                    "Using last entry (%g)\n",
                    resnm, atomnm, p, ap->value[j], line, ap->db, p);
            ap->value[j] = p;
        }
    }
    else
    {
        ap->bAvail[j] = TRUE;
        ap->value[j]  = p;
    }
}

static void read_prop(gmx_atomprop_t aps, int eprop, double factor)
{
    gmx_atomprop *ap2 = (gmx_atomprop*) aps;
    FILE         *fp;
    char          line[STRLEN], resnm[32], atomnm[32];
    double        pp;
    int           line_no;
    aprop_t      *ap;

    ap = &ap2->prop[eprop];

    fp      = libopen(ap->db);
    line_no = 0;
    while (get_a_line(fp, line, STRLEN))
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
                    ap->db, line_no);
        }
    }

    /* for libraries we can use the low-level close routines */
    gmx_ffclose(fp);

    ap->bSet = TRUE;
}

static void set_prop(gmx_atomprop_t aps, int eprop)
{
    gmx_atomprop *ap2           = (gmx_atomprop*) aps;
    const char   *fns[epropNR]  = { "atommass.dat", "vdwradii.dat", "dgsolv.dat", "electroneg.dat", "elements.dat" };
    double        fac[epropNR]  = { 1.0,    1.0,  418.4, 1.0, 1.0 };
    double        def[epropNR]  = { 12.011, 0.14, 0.0, 2.2, -1 };
    aprop_t      *ap;

    ap = &ap2->prop[eprop];
    if (!ap->bSet)
    {
        ap->db  = gmx_strdup(fns[eprop]);
        ap->def = def[eprop];
        read_prop(aps, eprop, fac[eprop]);

        if (debug)
        {
            fprintf(debug, "Entries in %s: %d\n", ap->db, ap->nprop);
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

gmx_atomprop_t gmx_atomprop_init(void)
{
    gmx_atomprop *aps;

    snew(aps, 1);

    gmx_residuetype_init(&aps->restype);
    aps->bWarned  = FALSE;
    aps->bWarnVDW = FALSE;

    return (gmx_atomprop_t)aps;
}

static void destroy_prop(aprop_t *ap)
{
    int i;

    if (ap->bSet)
    {
        sfree(ap->db);

        for (i = 0; i < ap->nprop; i++)
        {
            sfree(ap->atomnm[i]);
            sfree(ap->resnm[i]);
        }
        sfree(ap->atomnm);
        sfree(ap->resnm);
        sfree(ap->bAvail);
        sfree(ap->value);
    }
}

void gmx_atomprop_destroy(gmx_atomprop_t aps)
{
    gmx_atomprop *ap = (gmx_atomprop*) aps;
    int           p;

    if (aps == NULL)
    {
        printf("\nWARNING: gmx_atomprop_destroy called with a NULL pointer\n\n");
        return;
    }

    for (p = 0; p < epropNR; p++)
    {
        destroy_prop(&ap->prop[p]);
    }

    gmx_residuetype_destroy(ap->restype);

    sfree(ap);
}

static void vdw_warning(FILE *fp)
{
    if (NULL != fp)
    {
        fprintf(fp, "NOTE: From version 5.0 %s uses the Van der Waals radii\n",
                ShortProgram());
        fprintf(fp, "from the source below. This means the results may be different\n");
        fprintf(fp, "compared to previous GROMACS versions.\n");
        please_cite(fp, "Bondi1964a");
    }
}

gmx_bool gmx_atomprop_query(gmx_atomprop_t aps,
                            int eprop, const char *resnm, const char *atomnm,
                            real *value)
{
    gmx_atomprop *ap = (gmx_atomprop*) aps;
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

    j = get_prop_index(&(ap->prop[eprop]), ap->restype, resname,
                       atomname, &bExact);

    if (eprop == epropVDW && !ap->bWarnVDW)
    {
        vdw_warning(stdout);
        ap->bWarnVDW = TRUE;
    }
    if (j >= 0)
    {
        *value = ap->prop[eprop].value[j];
        return TRUE;
    }
    else
    {
        *value = ap->prop[eprop].def;
        return FALSE;
    }
}

char *gmx_atomprop_element(gmx_atomprop_t aps, int atomnumber)
{
    gmx_atomprop *ap = (gmx_atomprop*) aps;
    int           i;

    set_prop(aps, epropElement);
    for (i = 0; (i < ap->prop[epropElement].nprop); i++)
    {
        if (gmx_nint(ap->prop[epropElement].value[i]) == atomnumber)
        {
            return ap->prop[epropElement].atomnm[i];
        }
    }
    return NULL;
}

int gmx_atomprop_atomnumber(gmx_atomprop_t aps, const char *elem)
{
    gmx_atomprop *ap = (gmx_atomprop*) aps;
    int           i;

    set_prop(aps, epropElement);
    for (i = 0; (i < ap->prop[epropElement].nprop); i++)
    {
        if (gmx_strcasecmp(ap->prop[epropElement].atomnm[i], elem) == 0)
        {
            return gmx_nint(ap->prop[epropElement].value[i]);
        }
    }
    return -1;
}
