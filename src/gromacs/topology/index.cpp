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

#include "index.h"

#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/strdb.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/invblock.h"
#include "gromacs/topology/residuetypes.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

static gmx_bool gmx_ask_yesno(gmx_bool bASK)
{
    char c;

    if (bASK)
    {
        do
        {
            c = toupper(fgetc(stdin));
        }
        while ((c != 'Y') && (c != 'N'));

        return (c == 'Y');
    }
    else
    {
        return FALSE;
    }
}

void write_index(const char *outf, t_blocka *b, char **gnames, gmx_bool bDuplicate, int natoms)
{
    FILE *out;
    int   i, j, k;

    out = gmx_fio_fopen(outf, "w");
    /* fprintf(out,"%5d  %5d\n",b->nr,b->nra); */
    for (i = 0; (i < b->nr); i++)
    {
        fprintf(out, "[ %s ]\n", gnames[i]);
        for (k = 0, j = b->index[i]; j < b->index[i+1]; j++, k++)
        {
            fprintf(out, "%4d ", b->a[j]+1);
            if ((k % 15) == 14)
            {
                fprintf(out, "\n");
            }
        }
        fprintf(out, "\n");
    }

    /* Duplicate copy, useful for computational electrophysiology double-layer setups */
    if (bDuplicate)
    {
        fprintf(stderr, "Duplicating the whole system with an atom offset of %d atoms.\n", natoms);
        for (i = 0; (i < b->nr); i++)
        {
            fprintf(out, "[ %s_copy ]\n", gnames[i]);
            for (k = 0, j = b->index[i]; j < b->index[i+1]; j++, k++)
            {
                fprintf(out, "%4d ", b->a[j]+1 + natoms );
                if ((k % 15) == 14)
                {
                    fprintf(out, "\n");
                }
            }
            fprintf(out, "\n");
        }
    }

    gmx_fio_fclose(out);
}

void add_grp(t_blocka *b, char ***gnames, int nra, atom_id a[], const char *name)
{
    int i;

    srenew(b->index, b->nr+2);
    srenew(*gnames, b->nr+1);
    (*gnames)[b->nr] = gmx_strdup(name);

    srenew(b->a, b->nra+nra);
    for (i = 0; (i < nra); i++)
    {
        b->a[b->nra++] = a[i];
    }
    b->nr++;
    b->index[b->nr] = b->nra;
}

/* compare index in `a' with group in `b' at `index',
   when `index'<0 it is relative to end of `b' */
static gmx_bool grp_cmp(t_blocka *b, int nra, atom_id a[], int index)
{
    int i;

    if (index < 0)
    {
        index = b->nr-1+index;
    }
    if (index >= b->nr)
    {
        gmx_fatal(FARGS, "no such index group %d in t_blocka (nr=%d)", index, b->nr);
    }
    /* compare sizes */
    if (nra != b->index[index+1] - b->index[index])
    {
        return FALSE;
    }
    for (i = 0; i < nra; i++)
    {
        if (a[i] != b->a[b->index[index]+i])
        {
            return FALSE;
        }
    }
    return TRUE;
}

static void
p_status(const char *const *restype, int nres,
         const char *const *typenames, int ntypes)
{
    int   i, j;
    int * counter;

    snew(counter, ntypes);
    for (i = 0; i < ntypes; i++)
    {
        counter[i] = 0;
    }
    for (i = 0; i < nres; i++)
    {
        for (j = 0; j < ntypes; j++)
        {
            if (!gmx_strcasecmp(restype[i], typenames[j]))
            {
                counter[j]++;
            }
        }
    }

    for (i = 0; (i < ntypes); i++)
    {
        if (counter[i] > 0)
        {
            printf("There are: %5d %10s residues\n", counter[i], typenames[i]);
        }
    }

    sfree(counter);
}


static atom_id *
mk_aid(t_atoms *atoms, const char ** restype, const char * typestring, int *nra, gmx_bool bMatch)
/* Make an array of atom_ids for all atoms with residuetypes matching typestring, or the opposite if bMatch is false */
{
    atom_id *a;
    int      i;
    int      res;

    snew(a, atoms->nr);
    *nra = 0;
    for (i = 0; (i < atoms->nr); i++)
    {
        res = !gmx_strcasecmp(restype[atoms->atom[i].resind], typestring);
        if (bMatch == FALSE)
        {
            res = !res;
        }
        if (res)
        {
            a[(*nra)++] = i;
        }
    }

    return a;
}

typedef struct {
    char    *rname;
    gmx_bool bNeg;
    char    *gname;
} restp_t;

static void analyse_other(const char ** restype, t_atoms *atoms,
                          t_blocka *gb, char ***gn, gmx_bool bASK, gmx_bool bVerb)
{
    restp_t *restp = NULL;
    char   **attp  = NULL;
    char    *rname, *aname;
    atom_id *aid, *aaid;
    int      i, j, k, l, resind, naid, naaid, natp, nrestp = 0;

    for (i = 0; (i < atoms->nres); i++)
    {
        if (gmx_strcasecmp(restype[i], "Protein") && gmx_strcasecmp(restype[i], "DNA") && gmx_strcasecmp(restype[i], "RNA") && gmx_strcasecmp(restype[i], "Water"))
        {
            break;
        }
    }
    if (i < atoms->nres)
    {
        /* we have others */
        if (bVerb)
        {
            printf("Analysing residues not classified as Protein/DNA/RNA/Water and splitting into groups...\n");
        }
        for (k = 0; (k < atoms->nr); k++)
        {
            resind = atoms->atom[k].resind;
            rname  = *atoms->resinfo[resind].name;
            if (gmx_strcasecmp(restype[resind], "Protein") && gmx_strcasecmp(restype[resind], "DNA") &&
                gmx_strcasecmp(restype[resind], "RNA") && gmx_strcasecmp(restype[resind], "Water"))
            {

                for (l = 0; (l < nrestp); l++)
                {
                    if (strcmp(restp[l].rname, rname) == 0)
                    {
                        break;
                    }
                }
                if (l == nrestp)
                {
                    srenew(restp, nrestp+1);
                    restp[nrestp].rname = gmx_strdup(rname);
                    restp[nrestp].bNeg  = FALSE;
                    restp[nrestp].gname = gmx_strdup(rname);
                    nrestp++;
                }
            }
        }
        for (i = 0; (i < nrestp); i++)
        {
            snew(aid, atoms->nr);
            naid = 0;
            for (j = 0; (j < atoms->nr); j++)
            {
                rname = *atoms->resinfo[atoms->atom[j].resind].name;
                if ((strcmp(restp[i].rname, rname) == 0 && !restp[i].bNeg) ||
                    (strcmp(restp[i].rname, rname) != 0 &&  restp[i].bNeg))
                {
                    aid[naid++] = j;
                }
            }
            add_grp(gb, gn, naid, aid, restp[i].gname);
            if (bASK)
            {
                printf("split %s into atoms (y/n) ? ", restp[i].gname);
                fflush(stdout);
                if (gmx_ask_yesno(bASK))
                {
                    natp = 0;
                    for (k = 0; (k < naid); k++)
                    {
                        aname = *atoms->atomname[aid[k]];
                        for (l = 0; (l < natp); l++)
                        {
                            if (strcmp(aname, attp[l]) == 0)
                            {
                                break;
                            }
                        }
                        if (l == natp)
                        {
                            srenew(attp, ++natp);
                            attp[natp-1] = aname;
                        }
                    }
                    if (natp > 1)
                    {
                        for (l = 0; (l < natp); l++)
                        {
                            snew(aaid, naid);
                            naaid = 0;
                            for (k = 0; (k < naid); k++)
                            {
                                aname = *atoms->atomname[aid[k]];
                                if (strcmp(aname, attp[l]) == 0)
                                {
                                    aaid[naaid++] = aid[k];
                                }
                            }
                            add_grp(gb, gn, naaid, aaid, attp[l]);
                            sfree(aaid);
                        }
                    }
                    sfree(attp);
                    attp = NULL;
                }
            }
            sfree(aid);
            sfree(restp[i].rname);
            sfree(restp[i].gname);
        }
        sfree(restp);
    }
}

/*! \brief
 * Cata necessary to construct a single (protein) index group in
 * analyse_prot().
 */
typedef struct gmx_help_make_index_group
{
    /** The set of atom names that will be used to form this index group */
    const char **defining_atomnames;
    /** Size of the defining_atomnames array */
    int          num_defining_atomnames;
    /** Name of this index group */
    const char  *group_name;
    /** Whether the above atom names name the atoms in the group, or
        those not in the group */
    gmx_bool bTakeComplement;
    /** The index in wholename gives the first item in the arrays of
       atomnames that should be tested with 'gmx_strncasecmp' in stead of
       gmx_strcasecmp, or -1 if all items should be tested with strcasecmp
       This is comparable to using a '*' wildcard at the end of specific
       atom names, but that is more involved to implement...
     */
    int wholename;
    /** Only create this index group if it differs from the one specified in compareto,
       where -1 means to always create this group. */
    int compareto;
} t_gmx_help_make_index_group;

static void analyse_prot(const char ** restype, t_atoms *atoms,
                         t_blocka *gb, char ***gn, gmx_bool bASK, gmx_bool bVerb)
{
    /* lists of atomnames to be used in constructing index groups: */
    static const char *pnoh[]    = { "H", "HN" };
    static const char *pnodum[]  = {
        "MN1",  "MN2",  "MCB1", "MCB2", "MCG1", "MCG2",
        "MCD1", "MCD2", "MCE1", "MCE2", "MNZ1", "MNZ2"
    };
    static const char *calpha[]  = { "CA" };
    static const char *bb[]      = { "N", "CA", "C" };
    static const char *mc[]      = { "N", "CA", "C", "O", "O1", "O2", "OC1", "OC2", "OT", "OXT" };
    static const char *mcb[]     = { "N", "CA", "CB", "C", "O", "O1", "O2", "OC1", "OC2", "OT", "OXT" };
    static const char *mch[]     = {
        "N", "CA", "C", "O", "O1", "O2", "OC1", "OC2", "OT", "OXT",
        "H1", "H2", "H3", "H", "HN"
    };

    static const t_gmx_help_make_index_group constructing_data[] =
    {{ NULL,   0, "Protein",      TRUE,  -1, -1},
     { pnoh,   asize(pnoh),   "Protein-H",    TRUE,  0,  -1},
     { calpha, asize(calpha), "C-alpha",      FALSE, -1, -1},
     { bb,     asize(bb),     "Backbone",     FALSE, -1, -1},
     { mc,     asize(mc),     "MainChain",    FALSE, -1, -1},
     { mcb,    asize(mcb),    "MainChain+Cb", FALSE, -1, -1},
     { mch,    asize(mch),    "MainChain+H",  FALSE, -1, -1},
     { mch,    asize(mch),    "SideChain",    TRUE,  -1, -1},
     { mch,    asize(mch),    "SideChain-H",  TRUE,  11, -1},
     { pnodum, asize(pnodum), "Prot-Masses",  TRUE,  -1, 0}, };
    const int   num_index_groups = asize(constructing_data);

    int         n, j;
    atom_id    *aid;
    int         nra, npres;
    gmx_bool    match;
    char        ndx_name[STRLEN], *atnm;
    int         i;

    if (bVerb)
    {
        printf("Analysing Protein...\n");
    }
    snew(aid, atoms->nr);

    /* calculate the number of protein residues */
    npres = 0;
    for (i = 0; (i < atoms->nres); i++)
    {
        if (0 == gmx_strcasecmp(restype[i], "Protein"))
        {
            npres++;
        }
    }
    /* find matching or complement atoms */
    for (i = 0; (i < (int)num_index_groups); i++)
    {
        nra = 0;
        for (n = 0; (n < atoms->nr); n++)
        {
            if (0 == gmx_strcasecmp(restype[atoms->atom[n].resind], "Protein"))
            {
                match = FALSE;
                for (j = 0; (j < constructing_data[i].num_defining_atomnames); j++)
                {
                    /* skip digits at beginning of atomname, e.g. 1H */
                    atnm = *atoms->atomname[n];
                    while (isdigit(atnm[0]))
                    {
                        atnm++;
                    }
                    if ( (constructing_data[i].wholename == -1) || (j < constructing_data[i].wholename) )
                    {
                        if (0 == gmx_strcasecmp(constructing_data[i].defining_atomnames[j], atnm))
                        {
                            match = TRUE;
                        }
                    }
                    else
                    {
                        if (0 == gmx_strncasecmp(constructing_data[i].defining_atomnames[j], atnm, strlen(constructing_data[i].defining_atomnames[j])))
                        {
                            match = TRUE;
                        }
                    }
                }
                if (constructing_data[i].bTakeComplement != match)
                {
                    aid[nra++] = n;
                }
            }
        }
        /* if we want to add this group always or it differs from previous
           group, add it: */
        if (-1 == constructing_data[i].compareto || !grp_cmp(gb, nra, aid, constructing_data[i].compareto-i) )
        {
            add_grp(gb, gn, nra, aid, constructing_data[i].group_name);
        }
    }

    if (bASK)
    {
        for (i = 0; (i < (int)num_index_groups); i++)
        {
            printf("Split %12s into %5d residues (y/n) ? ", constructing_data[i].group_name, npres);
            if (gmx_ask_yesno(bASK))
            {
                int resind;
                nra = 0;
                for (n = 0; ((atoms->atom[n].resind < npres) && (n < atoms->nr)); )
                {
                    resind = atoms->atom[n].resind;
                    for (; ((atoms->atom[n].resind == resind) && (n < atoms->nr)); n++)
                    {
                        match = FALSE;
                        for (j = 0; (j < constructing_data[i].num_defining_atomnames); j++)
                        {
                            if (0 == gmx_strcasecmp(constructing_data[i].defining_atomnames[j], *atoms->atomname[n]))
                            {
                                match = TRUE;
                            }
                        }
                        if (constructing_data[i].bTakeComplement != match)
                        {
                            aid[nra++] = n;
                        }
                    }
                    /* copy the residuename to the tail of the groupname */
                    if (nra > 0)
                    {
                        t_resinfo *ri;
                        ri = &atoms->resinfo[resind];
                        sprintf(ndx_name, "%s_%s%d%c",
                                constructing_data[i].group_name, *ri->name, ri->nr, ri->ic == ' ' ? '\0' : ri->ic);
                        add_grp(gb, gn, nra, aid, ndx_name);
                        nra = 0;
                    }
                }
            }
        }
        printf("Make group with sidechain and C=O swapped (y/n) ? ");
        if (gmx_ask_yesno(bASK))
        {
            /* Make swap sidechain C=O index */
            int resind, hold;
            nra = 0;
            for (n = 0; ((atoms->atom[n].resind < npres) && (n < atoms->nr)); )
            {
                resind = atoms->atom[n].resind;
                hold   = -1;
                for (; ((atoms->atom[n].resind == resind) && (n < atoms->nr)); n++)
                {
                    if (strcmp("CA", *atoms->atomname[n]) == 0)
                    {
                        aid[nra++] = n;
                        hold       = nra;
                        nra       += 2;
                    }
                    else if (strcmp("C", *atoms->atomname[n]) == 0)
                    {
                        if (hold == -1)
                        {
                            gmx_incons("Atom naming problem");
                        }
                        aid[hold] = n;
                    }
                    else if (strcmp("O", *atoms->atomname[n]) == 0)
                    {
                        if (hold == -1)
                        {
                            gmx_incons("Atom naming problem");
                        }
                        aid[hold+1] = n;
                    }
                    else if (strcmp("O1", *atoms->atomname[n]) == 0)
                    {
                        if (hold == -1)
                        {
                            gmx_incons("Atom naming problem");
                        }
                        aid[hold+1] = n;
                    }
                    else
                    {
                        aid[nra++] = n;
                    }
                }
            }
            /* copy the residuename to the tail of the groupname */
            if (nra > 0)
            {
                add_grp(gb, gn, nra, aid, "SwapSC-CO");
            }
        }
    }
    sfree(aid);
}


void analyse(t_atoms *atoms, t_blocka *gb, char ***gn, gmx_bool bASK, gmx_bool bVerb)
{
    gmx_residuetype_t*rt = NULL;
    char             *resnm;
    atom_id          *aid;
    const char    **  restype;
    int               nra;
    int               i, k;
    int               ntypes;
    char           ** p_typename;
    int               iwater, iion;
    int               nwater, nion;
    int               found;

    if (bVerb)
    {
        printf("Analysing residue names:\n");
    }
    /* Create system group, every single atom */
    snew(aid, atoms->nr);
    for (i = 0; i < atoms->nr; i++)
    {
        aid[i] = i;
    }
    add_grp(gb, gn, atoms->nr, aid, "System");
    sfree(aid);

    /* For every residue, get a pointer to the residue type name */
    gmx_residuetype_init(&rt);
    assert(rt);

    snew(restype, atoms->nres);
    ntypes     = 0;
    p_typename = NULL;
    for (i = 0; i < atoms->nres; i++)
    {
        resnm = *atoms->resinfo[i].name;
        gmx_residuetype_get_type(rt, resnm, &(restype[i]));

        /* Note that this does not lead to a N*N loop, but N*K, where
         * K is the number of residue _types_, which is small and independent of N.
         */
        found = 0;
        for (k = 0; k < ntypes && !found; k++)
        {
            assert(p_typename != NULL);
            found = !strcmp(restype[i], p_typename[k]);
        }
        if (!found)
        {
            srenew(p_typename, ntypes+1);
            p_typename[ntypes++] = gmx_strdup(restype[i]);
        }
    }

    if (bVerb)
    {
        p_status(restype, atoms->nres, p_typename, ntypes);
    }

    for (k = 0; k < ntypes; k++)
    {
        aid = mk_aid(atoms, restype, p_typename[k], &nra, TRUE);

        /* Check for special types to do fancy stuff with */

        if (!gmx_strcasecmp(p_typename[k], "Protein") && nra > 0)
        {
            sfree(aid);
            /* PROTEIN */
            analyse_prot(restype, atoms, gb, gn, bASK, bVerb);

            /* Create a Non-Protein group */
            aid = mk_aid(atoms, restype, "Protein", &nra, FALSE);
            if ((nra > 0) && (nra < atoms->nr))
            {
                add_grp(gb, gn, nra, aid, "non-Protein");
            }
            sfree(aid);
        }
        else if (!gmx_strcasecmp(p_typename[k], "Water") && nra > 0)
        {
            add_grp(gb, gn, nra, aid, p_typename[k]);
            /* Add this group as 'SOL' too, for backward compatibility with older gromacs versions */
            add_grp(gb, gn, nra, aid, "SOL");

            sfree(aid);

            /* Solvent, create a negated group too */
            aid = mk_aid(atoms, restype, "Water", &nra, FALSE);
            if ((nra > 0) && (nra < atoms->nr))
            {
                add_grp(gb, gn, nra, aid, "non-Water");
            }
            sfree(aid);
        }
        else if (nra > 0)
        {
            /* Other groups */
            add_grp(gb, gn, nra, aid, p_typename[k]);
            sfree(aid);
            analyse_other(restype, atoms, gb, gn, bASK, bVerb);
        }
        sfree(p_typename[k]);
    }

    sfree(p_typename);
    sfree(restype);
    gmx_residuetype_destroy(rt);

    /* Create a merged water_and_ions group */
    iwater = -1;
    iion   = -1;
    nwater = 0;
    nion   = 0;

    for (i = 0; i < gb->nr; i++)
    {
        if (!gmx_strcasecmp((*gn)[i], "Water"))
        {
            iwater = i;
            nwater = gb->index[i+1]-gb->index[i];
        }
        else if (!gmx_strcasecmp((*gn)[i], "Ion"))
        {
            iion = i;
            nion = gb->index[i+1]-gb->index[i];
        }
    }

    if (nwater > 0 && nion > 0)
    {
        srenew(gb->index, gb->nr+2);
        srenew(*gn, gb->nr+1);
        (*gn)[gb->nr] = gmx_strdup("Water_and_ions");
        srenew(gb->a, gb->nra+nwater+nion);
        if (nwater > 0)
        {
            for (i = gb->index[iwater]; i < gb->index[iwater+1]; i++)
            {
                gb->a[gb->nra++] = gb->a[i];
            }
        }
        if (nion > 0)
        {
            for (i = gb->index[iion]; i < gb->index[iion+1]; i++)
            {
                gb->a[gb->nra++] = gb->a[i];
            }
        }
        gb->nr++;
        gb->index[gb->nr] = gb->nra;
    }
}


void check_index(char *gname, int n, atom_id index[], char *traj, int natoms)
{
    int i;

    for (i = 0; i < n; i++)
    {
        if (index[i] >= natoms)
        {
            gmx_fatal(FARGS, "%s atom number (index[%d]=%d) is larger than the number of atoms in %s (%d)",
                      gname ? gname : "Index", i+1, index[i]+1,
                      traj ? traj : "the trajectory", natoms);
        }
        else if (index[i] < 0)
        {
            gmx_fatal(FARGS, "%s atom number (index[%d]=%d) is less than zero",
                      gname ? gname : "Index", i+1, index[i]+1);
        }
    }
}

t_blocka *init_index(const char *gfile, char ***grpname)
{
    FILE      *in;
    t_blocka  *b;
    int        maxentries;
    int        i, j;
    char       line[STRLEN], *pt, str[STRLEN];

    in = gmx_fio_fopen(gfile, "r");
    snew(b, 1);
    b->nr      = 0;
    b->index   = NULL;
    b->nra     = 0;
    b->a       = NULL;
    *grpname   = NULL;
    maxentries = 0;
    while (get_a_line(in, line, STRLEN))
    {
        if (get_header(line, str))
        {
            b->nr++;
            srenew(b->index, b->nr+1);
            srenew(*grpname, b->nr);
            if (b->nr == 1)
            {
                b->index[0] = 0;
            }
            b->index[b->nr]     = b->index[b->nr-1];
            (*grpname)[b->nr-1] = gmx_strdup(str);
        }
        else
        {
            if (b->nr == 0)
            {
                gmx_fatal(FARGS, "The first header of your indexfile is invalid");
            }
            pt = line;
            while (sscanf(pt, "%s", str) == 1)
            {
                i = b->index[b->nr];
                if (i >= maxentries)
                {
                    maxentries += 1024;
                    srenew(b->a, maxentries);
                }
                assert(b->a != NULL); // for clang analyzer
                b->a[i] = strtol(str, NULL, 10)-1;
                b->index[b->nr]++;
                (b->nra)++;
                pt = strstr(pt, str)+strlen(str);
            }
        }
    }
    gmx_fio_fclose(in);

    for (i = 0; (i < b->nr); i++)
    {
        assert(b->a != NULL); // for clang analyzer
        for (j = b->index[i]; (j < b->index[i+1]); j++)
        {
            if (b->a[j] < 0)
            {
                fprintf(stderr, "\nWARNING: negative index %d in group %s\n\n",
                        b->a[j], (*grpname)[i]);
            }
        }
    }

    return b;
}

static void minstring(char *str)
{
    int i;

    for (i = 0; (i < (int)strlen(str)); i++)
    {
        if (str[i] == '-')
        {
            str[i] = '_';
        }
    }
}

int find_group(char s[], int ngrps, char **grpname)
{
    int      aa, i, n;
    char     string[STRLEN];
    gmx_bool bMultiple;

    bMultiple = FALSE;
    n         = strlen(s);
    aa        = NOTSET;
    /* first look for whole name match */
    if (aa == NOTSET)
    {
        for (i = 0; i < ngrps; i++)
        {
            if (gmx_strcasecmp_min(s, grpname[i]) == 0)
            {
                if (aa != NOTSET)
                {
                    bMultiple = TRUE;
                }
                aa = i;
            }
        }
    }
    /* second look for first string match */
    if (aa == NOTSET)
    {
        for (i = 0; i < ngrps; i++)
        {
            if (gmx_strncasecmp_min(s, grpname[i], n) == 0)
            {
                if (aa != NOTSET)
                {
                    bMultiple = TRUE;
                }
                aa = i;
            }
        }
    }
    /* last look for arbitrary substring match */
    if (aa == NOTSET)
    {
        upstring(s);
        minstring(s);
        for (i = 0; i < ngrps; i++)
        {
            strcpy(string, grpname[i]);
            upstring(string);
            minstring(string);
            if (strstr(string, s) != NULL)
            {
                if (aa != NOTSET)
                {
                    bMultiple = TRUE;
                }
                aa = i;
            }
        }
    }
    if (bMultiple)
    {
        printf("Error: Multiple groups '%s' selected\n", s);
        aa = NOTSET;
    }
    return aa;
}

static int qgroup(int *a, int ngrps, char **grpname)
{
    char     s[STRLEN];
    int      aa;
    gmx_bool bInRange;
    char    *end;

    do
    {
        fprintf(stderr, "Select a group: ");
        do
        {
            if (scanf("%s", s) != 1)
            {
                gmx_fatal(FARGS, "Cannot read from input");
            }
            trim(s); /* remove spaces */
        }
        while (strlen(s) == 0);
        aa = strtol(s, &end, 10);
        if (aa == 0 && end[0] != '\0') /* string entered */
        {
            aa = find_group(s, ngrps, grpname);
        }
        bInRange = (aa >= 0 && aa < ngrps);
        if (!bInRange)
        {
            printf("Error: No such group '%s'\n", s);
        }
    }
    while (!bInRange);
    printf("Selected %d: '%s'\n", aa, grpname[aa]);
    *a = aa;
    return aa;
}

static void rd_groups(t_blocka *grps, char **grpname, char *gnames[],
                      int ngrps, int isize[], atom_id *index[], int grpnr[])
{
    int i, j, gnr1;

    if (grps->nr == 0)
    {
        gmx_fatal(FARGS, "Error: no groups in indexfile");
    }
    for (i = 0; (i < grps->nr); i++)
    {
        fprintf(stderr, "Group %5d (%15s) has %5d elements\n", i, grpname[i],
                grps->index[i+1]-grps->index[i]);
    }
    for (i = 0; (i < ngrps); i++)
    {
        if (grps->nr > 1)
        {
            do
            {
                gnr1 = qgroup(&grpnr[i], grps->nr, grpname);
                if ((gnr1 < 0) || (gnr1 >= grps->nr))
                {
                    fprintf(stderr, "Select between %d and %d.\n", 0, grps->nr-1);
                }
            }
            while ((gnr1 < 0) || (gnr1 >= grps->nr));
        }
        else
        {
            fprintf(stderr, "There is one group in the index\n");
            gnr1 = 0;
        }
        gnames[i] = gmx_strdup(grpname[gnr1]);
        isize[i]  = grps->index[gnr1+1]-grps->index[gnr1];
        snew(index[i], isize[i]);
        for (j = 0; (j < isize[i]); j++)
        {
            index[i][j] = grps->a[grps->index[gnr1]+j];
        }
    }
}

void rd_index(const char *statfile, int ngrps, int isize[],
              atom_id *index[], char *grpnames[])
{
    char    **gnames;
    t_blocka *grps;
    int      *grpnr;

    snew(grpnr, ngrps);
    if (!statfile)
    {
        gmx_fatal(FARGS, "No index file specified");
    }
    grps = init_index(statfile, &gnames);
    rd_groups(grps, gnames, grpnames, ngrps, isize, index, grpnr);
}

void rd_index_nrs(char *statfile, int ngrps, int isize[],
                  atom_id *index[], char *grpnames[], int grpnr[])
{
    char    **gnames;
    t_blocka *grps;

    if (!statfile)
    {
        gmx_fatal(FARGS, "No index file specified");
    }
    grps = init_index(statfile, &gnames);

    rd_groups(grps, gnames, grpnames, ngrps, isize, index, grpnr);
}

void get_index(t_atoms *atoms, const char *fnm, int ngrps,
               int isize[], atom_id *index[], char *grpnames[])
{
    char    ***gnames;
    t_blocka  *grps = NULL;
    int       *grpnr;

    snew(grpnr, ngrps);
    snew(gnames, 1);
    if (fnm != NULL)
    {
        grps = init_index(fnm, gnames);
    }
    else if (atoms)
    {
        snew(grps, 1);
        snew(grps->index, 1);
        analyse(atoms, grps, gnames, FALSE, FALSE);
    }
    else
    {
        gmx_incons("You need to supply a valid atoms structure or a valid index file name");
    }

    rd_groups(grps, *gnames, grpnames, ngrps, isize, index, grpnr);
}

t_cluster_ndx *cluster_index(FILE *fplog, const char *ndx)
{
    t_cluster_ndx *c;
    int            i;

    snew(c, 1);
    c->clust     = init_index(ndx, &c->grpname);
    c->maxframe  = -1;
    for (i = 0; (i < c->clust->nra); i++)
    {
        c->maxframe = std::max(c->maxframe, c->clust->a[i]);
    }
    fprintf(fplog ? fplog : stdout,
            "There are %d clusters containing %d structures, highest framenr is %d\n",
            c->clust->nr, c->clust->nra, c->maxframe);
    if (debug)
    {
        pr_blocka(debug, 0, "clust", c->clust, TRUE);
        for (i = 0; (i < c->clust->nra); i++)
        {
            if ((c->clust->a[i] < 0) || (c->clust->a[i] > c->maxframe))
            {
                gmx_fatal(FARGS, "Range check error for c->clust->a[%d] = %d\n"
                          "should be within 0 and %d", i, c->clust->a[i], c->maxframe+1);
            }
        }
    }
    c->inv_clust = make_invblocka(c->clust, c->maxframe);

    return c;
}
