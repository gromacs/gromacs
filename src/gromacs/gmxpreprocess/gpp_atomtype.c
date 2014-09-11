/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2014, by the GROMACS development team, led by
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

#include "gpp_atomtype.h"

#include <math.h>
#include <string.h>

#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

typedef struct gpp_atomtype {
    int              nr;           /* The number of atomtypes		*/
    t_atom          *atom;         /* Array of atoms			*/
    char          ***atomname;     /* Names of the atomtypes		*/
    t_param         *nb;           /* Nonbonded force default params	*/
    int             *bondatomtype; /* The bond_atomtype for each atomtype  */
    real            *radius;       /* Radius for GBSA stuff                */
    real            *vol;          /* Effective volume for GBSA            */
    real            *surftens;     /* Surface tension with water, for GBSA */
    real            *gb_radius;    /* Radius for Still model               */
    real            *S_hct;        /* Overlap factor for HCT model         */
    int             *atomnumber;   /* Atomic number, used for QM/MM        */
} t_gpp_atomtype;

int get_atomtype_type(const char *str, gpp_atomtype_t ga)
{
    int i;

    /* Atom types are always case sensitive */
    for (i = 0; (i < ga->nr); i++)
    {
        if (strcmp(str, *(ga->atomname[i])) == 0)
        {
            return i;
        }
    }

    return NOTSET;
}

int get_atomtype_ntypes(gpp_atomtype_t ga)
{
    return ga->nr;
}

char *get_atomtype_name(int nt, gpp_atomtype_t ga)
{
    if ((nt < 0) || (nt >= ga->nr))
    {
        return NULL;
    }

    return *(ga->atomname[nt]);
}

real get_atomtype_massA(int nt, gpp_atomtype_t ga)
{
    if ((nt < 0) || (nt >= ga->nr))
    {
        return NOTSET;
    }

    return ga->atom[nt].m;
}

real get_atomtype_massB(int nt, gpp_atomtype_t ga)
{
    if ((nt < 0) || (nt >= ga->nr))
    {
        return NOTSET;
    }

    return ga->atom[nt].mB;
}

real get_atomtype_qA(int nt, gpp_atomtype_t ga)
{
    if ((nt < 0) || (nt >= ga->nr))
    {
        return NOTSET;
    }

    return ga->atom[nt].q;
}

real get_atomtype_qB(int nt, gpp_atomtype_t ga)
{
    if ((nt < 0) || (nt >= ga->nr))
    {
        return NOTSET;
    }

    return ga->atom[nt].qB;
}

int get_atomtype_ptype(int nt, gpp_atomtype_t ga)
{
    if ((nt < 0) || (nt >= ga->nr))
    {
        return NOTSET;
    }

    return ga->atom[nt].ptype;
}

int get_atomtype_batype(int nt, gpp_atomtype_t ga)
{
    if ((nt < 0) || (nt >= ga->nr))
    {
        return NOTSET;
    }

    return ga->bondatomtype[nt];
}

int get_atomtype_atomnumber(int nt, gpp_atomtype_t ga)
{
    if ((nt < 0) || (nt >= ga->nr))
    {
        return NOTSET;
    }

    return ga->atomnumber[nt];
}

real get_atomtype_radius(int nt, gpp_atomtype_t ga)
{
    if ((nt < 0) || (nt >= ga->nr))
    {
        return NOTSET;
    }

    return ga->radius[nt];
}

real get_atomtype_vol(int nt, gpp_atomtype_t ga)
{
    if ((nt < 0) || (nt >= ga->nr))
    {
        return NOTSET;
    }

    return ga->vol[nt];
}

real get_atomtype_surftens(int nt, gpp_atomtype_t ga)
{
    if ((nt < 0) || (nt >= ga->nr))
    {
        return NOTSET;
    }

    return ga->surftens[nt];
}

real get_atomtype_gb_radius(int nt, gpp_atomtype_t ga)
{
    if ((nt < 0) || (nt >= ga->nr))
    {
        return NOTSET;
    }

    return ga->gb_radius[nt];
}

real get_atomtype_S_hct(int nt, gpp_atomtype_t ga)
{
    if ((nt < 0) || (nt >= ga->nr))
    {
        return NOTSET;
    }

    return ga->S_hct[nt];
}

real get_atomtype_nbparam(int nt, int param, gpp_atomtype_t ga)
{
    if ((nt < 0) || (nt >= ga->nr))
    {
        return NOTSET;
    }
    if ((param < 0) || (param >= MAXFORCEPARAM))
    {
        return NOTSET;
    }
    return ga->nb[nt].c[param];
}

gpp_atomtype_t init_atomtype(void)
{
    gpp_atomtype_t ga;

    snew(ga, 1);

    ga->nr           = 0;
    ga->atom         = NULL;
    ga->atomname     = NULL;
    ga->nb           = NULL;
    ga->bondatomtype = NULL;
    ga->radius       = NULL;
    ga->vol          = NULL;
    ga->surftens     = NULL;
    ga->atomnumber   = NULL;
    ga->gb_radius    = NULL;
    ga->S_hct        = NULL;

    return ga;
}

int
set_atomtype_gbparam(gpp_atomtype_t ga, int i,
                     real radius, real vol, real surftens,
                     real gb_radius, real S_hct)
{
    if ( (i < 0) || (i >= ga->nr))
    {
        return NOTSET;
    }

    ga->radius[i]    = radius;
    ga->vol[i]       = vol;
    ga->surftens[i]  = surftens;
    ga->gb_radius[i] = gb_radius;
    ga->S_hct[i]     = S_hct;

    return i;
}


int set_atomtype(int nt, gpp_atomtype_t ga, t_symtab *tab,
                 t_atom *a, const char *name, t_param *nb,
                 int bondatomtype,
                 real radius, real vol, real surftens, int atomnumber,
                 real gb_radius, real S_hct)
{
    if ((nt < 0) || (nt >= ga->nr))
    {
        return NOTSET;
    }

    ga->atom[nt]         = *a;
    ga->atomname[nt]     = put_symtab(tab, name);
    ga->nb[nt]           = *nb;
    ga->bondatomtype[nt] = bondatomtype;
    ga->radius[nt]       = radius;
    ga->vol[nt]          = vol;
    ga->surftens[nt]     = surftens;
    ga->atomnumber[nt]   = atomnumber;
    ga->gb_radius[nt]    = gb_radius;
    ga->S_hct[nt]        = S_hct;

    return nt;
}

int add_atomtype(gpp_atomtype_t ga, t_symtab *tab,
                 t_atom *a, const char *name, t_param *nb,
                 int bondatomtype,
                 real radius, real vol, real surftens, real atomnumber,
                 real gb_radius, real S_hct)
{
    int i;

    for (i = 0; (i < ga->nr); i++)
    {
        if (strcmp(*ga->atomname[i], name) == 0)
        {
            if (NULL != debug)
            {
                fprintf(debug, "Trying to add atomtype %s again. Skipping it.\n", name);
            }
            break;
        }
    }
    if (i == ga->nr)
    {
        ga->nr++;
        srenew(ga->atom, ga->nr);
        srenew(ga->atomname, ga->nr);
        srenew(ga->nb, ga->nr);
        srenew(ga->bondatomtype, ga->nr);
        srenew(ga->radius, ga->nr);
        srenew(ga->vol, ga->nr);
        srenew(ga->surftens, ga->nr);
        srenew(ga->atomnumber, ga->nr);
        srenew(ga->gb_radius, ga->nr);
        srenew(ga->S_hct, ga->nr);

        return set_atomtype(ga->nr-1, ga, tab, a, name, nb, bondatomtype, radius,
                            vol, surftens, atomnumber, gb_radius, S_hct);
    }
    else
    {
        return i;
    }
}

void print_at (FILE * out, gpp_atomtype_t ga)
{
    int         i;
    t_atom     *atom = ga->atom;
    t_param    *nb   = ga->nb;

    fprintf (out, "[ %s ]\n", dir2str(d_atomtypes));
    fprintf (out, "; %6s  %8s  %8s  %8s  %12s  %12s\n",
             "type", "mass", "charge", "particle", "c6", "c12");
    for (i = 0; (i < ga->nr); i++)
    {
        fprintf(out, "%8s  %8.3f  %8.3f  %8s  %12e  %12e\n",
                *(ga->atomname[i]), atom[i].m, atom[i].q, "A",
                nb[i].C0, nb[i].C1);
    }

    fprintf (out, "\n");
}

void done_atomtype(gpp_atomtype_t ga)
{
    sfree(ga->atom);
    sfree(ga->atomname);
    sfree(ga->nb);
    sfree(ga->bondatomtype);
    sfree(ga->radius);
    sfree(ga->vol);
    sfree(ga->gb_radius);
    sfree(ga->S_hct);
    sfree(ga->surftens);
    sfree(ga->atomnumber);
    ga->nr = 0;
    sfree(ga);
}

static int search_atomtypes(gpp_atomtype_t ga, int *n, int typelist[],
                            int thistype,
                            t_param param[], int ftype)
{
    int      i, nn, nrfp, j, k, ntype, tli;
    gmx_bool bFound = FALSE;

    nn    = *n;
    nrfp  = NRFP(ftype);
    ntype = get_atomtype_ntypes(ga);

    for (i = 0; (i < nn); i++)
    {
        if (typelist[i] == thistype)
        {
            /* This type number has already been added */
            break;
        }

        /* Otherwise, check if the parameters are identical to any previously added type */

        bFound = TRUE;
        for (j = 0; j < ntype && bFound; j++)
        {
            /* Check nonbonded parameters */
            for (k = 0; k < nrfp && bFound; k++)
            {
                bFound = (param[ntype*typelist[i]+j].c[k] == param[ntype*thistype+j].c[k]);
            }

            /* Check radius, volume, surftens */
            tli    = typelist[i];
            bFound = bFound &&
                (get_atomtype_radius(tli, ga) == get_atomtype_radius(thistype, ga)) &&
                (get_atomtype_vol(tli, ga) == get_atomtype_vol(thistype, ga)) &&
                (get_atomtype_surftens(tli, ga) == get_atomtype_surftens(thistype, ga)) &&
                (get_atomtype_atomnumber(tli, ga) == get_atomtype_atomnumber(thistype, ga)) &&
                (get_atomtype_gb_radius(tli, ga) == get_atomtype_gb_radius(thistype, ga)) &&
                (get_atomtype_S_hct(tli, ga) == get_atomtype_S_hct(thistype, ga));
        }
        if (bFound)
        {
            break;
        }
    }

    if (i == nn)
    {
        if (debug)
        {
            fprintf(debug, "Renumbering atomtype %d to %d\n", thistype, nn);
        }
        if (nn == ntype)
        {
            gmx_fatal(FARGS, "Atomtype horror n = %d, %s, %d", nn, __FILE__, __LINE__);
        }
        typelist[nn] = thistype;
        nn++;
    }
    *n = nn;

    return i;
}

void renum_atype(t_params plist[], gmx_mtop_t *mtop,
                 int *wall_atomtype,
                 gpp_atomtype_t ga, gmx_bool bVerbose)
{
    int         i, j, k, l, molt, mi, mj, nat, nrfp, ftype, ntype;
    t_atoms    *atoms;
    t_param    *nbsnew;
    int        *typelist;
    real       *new_radius;
    real       *new_vol;
    real       *new_surftens;
    real       *new_gb_radius;
    real       *new_S_hct;
    int        *new_atomnumber;
    char     ***new_atomname;

    ntype = get_atomtype_ntypes(ga);
    snew(typelist, ntype);

    if (bVerbose)
    {
        fprintf(stderr, "renumbering atomtypes...\n");
    }

    /* Since the bonded interactions have been assigned now,
     * we want to reduce the number of atom types by merging
     * ones with identical nonbonded interactions, in addition
     * to removing unused ones.
     *
     * With Generalized-Born electrostatics, or implicit solvent
     * we also check that the atomtype radius, effective_volume
     * and surface tension match.
     *
     * With QM/MM we also check that the atom numbers match
     */

    /* Get nonbonded interaction type */
    if (plist[F_LJ].nr > 0)
    {
        ftype = F_LJ;
    }
    else
    {
        ftype = F_BHAM;
    }

    /* Renumber atomtypes by first making a list of which ones are actually used.
     * We provide the list of nonbonded parameters so search_atomtypes
     * can determine if two types should be merged.
     */
    nat = 0;
    for (molt = 0; molt < mtop->nmoltype; molt++)
    {
        atoms = &mtop->moltype[molt].atoms;
        for (i = 0; (i < atoms->nr); i++)
        {
            atoms->atom[i].type =
                search_atomtypes(ga, &nat, typelist, atoms->atom[i].type,
                                 plist[ftype].param, ftype);
            atoms->atom[i].typeB =
                search_atomtypes(ga, &nat, typelist, atoms->atom[i].typeB,
                                 plist[ftype].param, ftype);
        }
    }

    for (i = 0; i < 2; i++)
    {
        if (wall_atomtype[i] >= 0)
        {
            wall_atomtype[i] = search_atomtypes(ga, &nat, typelist, wall_atomtype[i],
                                                plist[ftype].param, ftype);
        }
    }

    snew(new_radius, nat);
    snew(new_vol, nat);
    snew(new_surftens, nat);
    snew(new_atomnumber, nat);
    snew(new_gb_radius, nat);
    snew(new_S_hct, nat);
    snew(new_atomname, nat);

    /* We now have a list of unique atomtypes in typelist */

    if (debug)
    {
        pr_ivec(debug, 0, "typelist", typelist, nat, TRUE);
    }

    /* Renumber nlist */
    nbsnew = NULL;
    snew(nbsnew, plist[ftype].nr);

    nrfp  = NRFP(ftype);

    for (i = k = 0; (i < nat); i++)
    {
        mi = typelist[i];
        for (j = 0; (j < nat); j++, k++)
        {
            mj = typelist[j];
            for (l = 0; (l < nrfp); l++)
            {
                nbsnew[k].c[l] = plist[ftype].param[ntype*mi+mj].c[l];
            }
        }
        new_radius[i]     = get_atomtype_radius(mi, ga);
        new_vol[i]        = get_atomtype_vol(mi, ga);
        new_surftens[i]   = get_atomtype_surftens(mi, ga);
        new_atomnumber[i] = get_atomtype_atomnumber(mi, ga);
        new_gb_radius[i]  = get_atomtype_gb_radius(mi, ga);
        new_S_hct[i]      = get_atomtype_S_hct(mi, ga);
        new_atomname[i]   = ga->atomname[mi];
    }

    for (i = 0; (i < nat*nat); i++)
    {
        for (l = 0; (l < nrfp); l++)
        {
            plist[ftype].param[i].c[l] = nbsnew[i].c[l];
        }
    }
    plist[ftype].nr     = i;
    mtop->ffparams.atnr = nat;

    sfree(ga->radius);
    sfree(ga->vol);
    sfree(ga->surftens);
    sfree(ga->atomnumber);
    sfree(ga->gb_radius);
    sfree(ga->S_hct);
    /* Dangling atomname pointers ? */
    sfree(ga->atomname);

    ga->radius     = new_radius;
    ga->vol        = new_vol;
    ga->surftens   = new_surftens;
    ga->atomnumber = new_atomnumber;
    ga->gb_radius  = new_gb_radius;
    ga->S_hct      = new_S_hct;
    ga->atomname   = new_atomname;

    ga->nr = nat;

    sfree(nbsnew);
    sfree(typelist);
}

void copy_atomtype_atomtypes(gpp_atomtype_t ga, t_atomtypes *atomtypes)
{
    int i, ntype;

    /* Copy the atomtype data to the topology atomtype list */
    ntype         = get_atomtype_ntypes(ga);
    atomtypes->nr = ntype;
    snew(atomtypes->radius, ntype);
    snew(atomtypes->vol, ntype);
    snew(atomtypes->surftens, ntype);
    snew(atomtypes->atomnumber, ntype);
    snew(atomtypes->gb_radius, ntype);
    snew(atomtypes->S_hct, ntype);

    for (i = 0; i < ntype; i++)
    {
        atomtypes->radius[i]     = ga->radius[i];
        atomtypes->vol[i]        = ga->vol[i];
        atomtypes->surftens[i]   = ga->surftens[i];
        atomtypes->atomnumber[i] = ga->atomnumber[i];
        atomtypes->gb_radius[i]  = ga->gb_radius[i];
        atomtypes->S_hct[i]      = ga->S_hct[i];
    }
}
