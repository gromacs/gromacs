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

#include "resall.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/fileio/strdb.h"
#include "gromacs/gmxpreprocess/fflibutil.h"
#include "gromacs/gmxpreprocess/pgutil.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

gpp_atomtype_t read_atype(const char *ffdir, t_symtab *tab)
{
    int            nfile, f;
    char         **file;
    FILE          *in;
    char           buf[STRLEN], name[STRLEN];
    double         m;
    int            nratt = 0;
    gpp_atomtype_t at;
    t_atom        *a;
    t_param       *nb;

    nfile = fflib_search_file_end(ffdir, ".atp", TRUE, &file);
    at    = init_atomtype();
    snew(a, 1);
    snew(nb, 1);

    for (f = 0; f < nfile; f++)
    {
        in = fflib_open(file[f]);
        while (!feof(in))
        {
            /* Skip blank or comment-only lines */
            do
            {
                if (fgets2(buf, STRLEN, in) != NULL)
                {
                    strip_comment(buf);
                    trim(buf);
                }
            }
            while (!feof(in) && strlen(buf) == 0);

            if (sscanf(buf, "%s%lf", name, &m) == 2)
            {
                a->m = m;
                add_atomtype(at, tab, a, name, nb, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
                fprintf(stderr, "\rAtomtype %d", ++nratt);
            }
            else
            {
                fprintf(stderr, "\nInvalid format: %s\n", buf);
            }
        }
        gmx_ffclose(in);
        sfree(file[f]);
    }
    fprintf(stderr, "\n");
    sfree(file);

    return at;
}

static void print_resatoms(FILE *out, gpp_atomtype_t atype, t_restp *rtp)
{
    int   j, tp;
    char *tpnm;

    /* fprintf(out,"%5s\n",rtp->resname);
       fprintf(out,"%5d\n",rtp->natom); */
    fprintf(out, "[ %s ]\n", rtp->resname);
    fprintf(out, " [ atoms ]\n");

    for (j = 0; (j < rtp->natom); j++)
    {
        tp   = rtp->atom[j].type;
        tpnm = get_atomtype_name(tp, atype);
        if (tpnm == NULL)
        {
            gmx_fatal(FARGS, "Incorrect atomtype (%d)", tp);
        }
        fprintf(out, "%6s  %6s  %8.3f  %6d\n",
                *(rtp->atomname[j]), tpnm, rtp->atom[j].q, rtp->cgnr[j]);
    }
}

static gmx_bool read_atoms(FILE *in, char *line,
                           t_restp *r0, t_symtab *tab, gpp_atomtype_t atype)
{
    int    i, j, cg, maxentries;
    char   buf[256], buf1[256];
    double q;

    /* Read Atoms */
    maxentries   = 0;
    r0->atom     =     NULL;
    r0->atomname = NULL;
    r0->cgnr     =     NULL;
    i            = 0;
    while (get_a_line(in, line, STRLEN) && (strchr(line, '[') == NULL))
    {
        if (sscanf(line, "%s%s%lf%d", buf, buf1, &q, &cg) != 4)
        {
            return FALSE;
        }
        if (i >= maxentries)
        {
            maxentries += 100;
            srenew(r0->atom,     maxentries);
            srenew(r0->atomname, maxentries);
            srenew(r0->cgnr,     maxentries);
        }
        r0->atomname[i] = put_symtab(tab, buf);
        r0->atom[i].q   = q;
        r0->cgnr[i]     = cg;
        j               = get_atomtype_type(buf1, atype);
        if (j == NOTSET)
        {
            gmx_fatal(FARGS, "Atom type %s (residue %s) not found in atomtype "
                      "database", buf1, r0->resname);
        }
        r0->atom[i].type = j;
        r0->atom[i].m    = get_atomtype_massA(j, atype);
        i++;
    }
    r0->natom = i;
    srenew(r0->atom, i);
    srenew(r0->atomname, i);
    srenew(r0->cgnr, i);

    return TRUE;
}

gmx_bool read_bondeds(int bt, FILE *in, char *line, t_restp *rtp)
{
    char str[STRLEN];
    int  j, n, ni, maxrb;

    maxrb = rtp->rb[bt].nb;
    while (get_a_line(in, line, STRLEN) && (strchr(line, '[') == NULL))
    {
        if (rtp->rb[bt].nb >= maxrb)
        {
            maxrb += 100;
            srenew(rtp->rb[bt].b, maxrb);
        }
        n = 0;
        for (j = 0; j < btsNiatoms[bt]; j++)
        {
            if (sscanf(line+n, "%s%n", str, &ni) == 1)
            {
                rtp->rb[bt].b[rtp->rb[bt].nb].a[j] = gmx_strdup(str);
            }
            else
            {
                return FALSE;
            }
            n += ni;
        }
        for (; j < MAXATOMLIST; j++)
        {
            rtp->rb[bt].b[rtp->rb[bt].nb].a[j] = NULL;
        }
        while (isspace(line[n]))
        {
            n++;
        }
        rtrim(line+n);
        rtp->rb[bt].b[rtp->rb[bt].nb].s = gmx_strdup(line+n);
        rtp->rb[bt].nb++;
    }
    /* give back unused memory */
    srenew(rtp->rb[bt].b, rtp->rb[bt].nb);

    return TRUE;
}

static void print_resbondeds(FILE *out, int bt, t_restp *rtp)
{
    int i, j;

    if (rtp->rb[bt].nb)
    {
        fprintf(out, " [ %s ]\n", btsNames[bt]);

        for (i = 0; i < rtp->rb[bt].nb; i++)
        {
            for (j = 0; j < btsNiatoms[bt]; j++)
            {
                fprintf(out, "%6s ", rtp->rb[bt].b[i].a[j]);
            }
            if (rtp->rb[bt].b[i].s[0])
            {
                fprintf(out, "    %s", rtp->rb[bt].b[i].s);
            }
            fprintf(out, "\n");
        }
    }
}

static void check_rtp(int nrtp, t_restp rtp[], char *libfn)
{
    int i;

    /* check for double entries, assuming list is already sorted */
    for (i = 1; (i < nrtp); i++)
    {
        if (gmx_strcasecmp(rtp[i-1].resname, rtp[i].resname) == 0)
        {
            fprintf(stderr, "WARNING double entry %s in file %s\n",
                    rtp[i].resname, libfn);
        }
    }
}

static int comprtp(const void *a, const void *b)
{
    t_restp *ra, *rb;

    ra = (t_restp *)a;
    rb = (t_restp *)b;

    return gmx_strcasecmp(ra->resname, rb->resname);
}

int get_bt(char* header)
{
    int i;

    for (i = 0; i < ebtsNR; i++)
    {
        if (gmx_strcasecmp(btsNames[i], header) == 0)
        {
            return i;
        }
    }
    return NOTSET;
}

void clear_t_restp(t_restp *rrtp)
{
    memset((void *)rrtp, 0, sizeof(t_restp));
}

/* print all the ebtsNR type numbers */
void print_resall_header(FILE *out, t_restp rtp[])
{
    fprintf(out, "[ bondedtypes ]\n");
    fprintf(out, "; bonds  angles  dihedrals  impropers all_dihedrals nr_exclusions  HH14  remove_dih\n");
    fprintf(out, " %5d  %6d  %9d  %9d  %14d  %14d %14d %14d\n\n",
            rtp[0].rb[0].type,
            rtp[0].rb[1].type,
            rtp[0].rb[2].type,
            rtp[0].rb[3].type,
            rtp[0].bKeepAllGeneratedDihedrals,
            rtp[0].nrexcl,
            rtp[0].bGenerateHH14Interactions,
            rtp[0].bRemoveDihedralIfWithImproper);
}

void print_resall(FILE *out, int nrtp, t_restp rtp[],
                  gpp_atomtype_t atype)
{
    int i, bt;

    if (nrtp == 0)
    {
        return;
    }

    print_resall_header(out, rtp);

    for (i = 0; i < nrtp; i++)
    {
        if (rtp[i].natom > 0)
        {
            print_resatoms(out, atype, &rtp[i]);
            for (bt = 0; bt < ebtsNR; bt++)
            {
                print_resbondeds(out, bt, &rtp[i]);
            }
        }
    }
}

void read_resall(char *rrdb, int *nrtpptr, t_restp **rtp,
                 gpp_atomtype_t atype, t_symtab *tab,
                 gmx_bool bAllowOverrideRTP)
{
    FILE         *in;
    char          filebase[STRLEN], *ptr, line[STRLEN], header[STRLEN];
    int           i, nrtp, maxrtp, bt, nparam;
    int           dum1, dum2, dum3;
    t_restp      *rrtp, *header_settings;
    gmx_bool      bNextResidue, bError;
    int           firstrtp;

    fflib_filename_base(rrdb, filebase, STRLEN);

    in = fflib_open(rrdb);

    if (debug)
    {
        fprintf(debug, "%9s %5s", "Residue", "atoms");
        for (i = 0; i < ebtsNR; i++)
        {
            fprintf(debug, " %10s", btsNames[i]);
        }
        fprintf(debug, "\n");
    }
    snew(header_settings, 1);

    /* these bonded parameters will overwritten be when  *
     * there is a [ bondedtypes ] entry in the .rtp file */
    header_settings->rb[ebtsBONDS].type  = 1; /* normal bonds     */
    header_settings->rb[ebtsANGLES].type = 1; /* normal angles    */
    header_settings->rb[ebtsPDIHS].type  = 1; /* normal dihedrals */
    header_settings->rb[ebtsIDIHS].type  = 2; /* normal impropers */
    header_settings->rb[ebtsEXCLS].type  = 1; /* normal exclusions */
    header_settings->rb[ebtsCMAP].type   = 1; /* normal cmap torsions */

    header_settings->bKeepAllGeneratedDihedrals    = FALSE;
    header_settings->nrexcl                        = 3;
    header_settings->bGenerateHH14Interactions     = TRUE;
    header_settings->bRemoveDihedralIfWithImproper = TRUE;

    /* Column 5 & 6 aren't really bonded types, but we include
     * them here to avoid introducing a new section:
     * Column 5 : This controls the generation of dihedrals from the bonding.
     *            All possible dihedrals are generated automatically. A value of
     *            1 here means that all these are retained. A value of
     *            0 here requires generated dihedrals be removed if
     *              * there are any dihedrals on the same central atoms
     *                specified in the residue topology, or
     *              * there are other identical generated dihedrals
     *                sharing the same central atoms, or
     *              * there are other generated dihedrals sharing the
     *                same central bond that have fewer hydrogen atoms
     * Column 6: Number of bonded neighbors to exclude.
     * Column 7: Generate 1,4 interactions between two hydrogen atoms
     * Column 8: Remove proper dihedrals if centered on the same bond
     *           as an improper dihedral
     */
    get_a_line(in, line, STRLEN);
    if (!get_header(line, header))
    {
        gmx_fatal(FARGS, "in .rtp file at line:\n%s\n", line);
    }
    if (gmx_strncasecmp("bondedtypes", header, 5) == 0)
    {
        get_a_line(in, line, STRLEN);
        if ((nparam = sscanf(line, "%d %d %d %d %d %d %d %d",
                             &header_settings->rb[ebtsBONDS].type, &header_settings->rb[ebtsANGLES].type,
                             &header_settings->rb[ebtsPDIHS].type, &header_settings->rb[ebtsIDIHS].type,
                             &dum1, &header_settings->nrexcl, &dum2, &dum3)) < 4)
        {
            gmx_fatal(FARGS, "need 4 to 8 parameters in the header of .rtp file %s at line:\n%s\n", rrdb, line);
        }
        header_settings->bKeepAllGeneratedDihedrals    = (dum1 != 0);
        header_settings->bGenerateHH14Interactions     = (dum2 != 0);
        header_settings->bRemoveDihedralIfWithImproper = (dum3 != 0);
        get_a_line(in, line, STRLEN);
        if (nparam < 5)
        {
            fprintf(stderr, "Using default: not generating all possible dihedrals\n");
            header_settings->bKeepAllGeneratedDihedrals = FALSE;
        }
        if (nparam < 6)
        {
            fprintf(stderr, "Using default: excluding 3 bonded neighbors\n");
            header_settings->nrexcl = 3;
        }
        if (nparam < 7)
        {
            fprintf(stderr, "Using default: generating 1,4 H--H interactions\n");
            header_settings->bGenerateHH14Interactions = TRUE;
        }
        if (nparam < 8)
        {
            fprintf(stderr, "Using default: removing proper dihedrals found on the same bond as a proper dihedral\n");
            header_settings->bRemoveDihedralIfWithImproper = TRUE;
        }
    }
    else
    {
        fprintf(stderr,
                "Reading .rtp file without '[ bondedtypes ]' directive,\n"
                "Will proceed as if the entry was:\n");
        print_resall_header(stderr, header_settings);
    }
    /* We don't know the current size of rrtp, but simply realloc immediately */
    nrtp   = *nrtpptr;
    rrtp   = *rtp;
    maxrtp = nrtp;
    while (!feof(in))
    {
        if (nrtp >= maxrtp)
        {
            maxrtp += 100;
            srenew(rrtp, maxrtp);
        }

        /* Initialise rtp entry structure */
        rrtp[nrtp] = *header_settings;
        if (!get_header(line, header))
        {
            gmx_fatal(FARGS, "in .rtp file at line:\n%s\n", line);
        }
        rrtp[nrtp].resname  = gmx_strdup(header);
        rrtp[nrtp].filebase = gmx_strdup(filebase);

        get_a_line(in, line, STRLEN);
        bError       = FALSE;
        bNextResidue = FALSE;
        do
        {
            if (!get_header(line, header))
            {
                bError = TRUE;
            }
            else
            {
                bt = get_bt(header);
                if (bt != NOTSET)
                {
                    /* header is an bonded directive */
                    bError = !read_bondeds(bt, in, line, &rrtp[nrtp]);
                }
                else if (gmx_strncasecmp("atoms", header, 5) == 0)
                {
                    /* header is the atoms directive */
                    bError = !read_atoms(in, line, &(rrtp[nrtp]), tab, atype);
                }
                else
                {
                    /* else header must be a residue name */
                    bNextResidue = TRUE;
                }
            }
            if (bError)
            {
                gmx_fatal(FARGS, "in .rtp file in residue %s at line:\n%s\n",
                          rrtp[nrtp].resname, line);
            }
        }
        while (!feof(in) && !bNextResidue);

        if (rrtp[nrtp].natom == 0)
        {
            gmx_fatal(FARGS, "No atoms found in .rtp file in residue %s\n",
                      rrtp[nrtp].resname);
        }
        if (debug)
        {
            fprintf(debug, "%3d %5s %5d",
                    nrtp+1, rrtp[nrtp].resname, rrtp[nrtp].natom);
            for (i = 0; i < ebtsNR; i++)
            {
                fprintf(debug, " %10d", rrtp[nrtp].rb[i].nb);
            }
            fprintf(debug, "\n");
        }

        firstrtp = -1;
        for (i = 0; i < nrtp; i++)
        {
            if (gmx_strcasecmp(rrtp[i].resname, rrtp[nrtp].resname) == 0)
            {
                firstrtp = i;
            }
        }

        if (firstrtp == -1)
        {
            nrtp++;
            fprintf(stderr, "\rResidue %d", nrtp);
        }
        else
        {
            if (firstrtp >= *nrtpptr)
            {
                gmx_fatal(FARGS, "Found a second entry for '%s' in '%s'",
                          rrtp[nrtp].resname, rrdb);
            }
            if (bAllowOverrideRTP)
            {
                fprintf(stderr, "\n");
                fprintf(stderr, "Found another rtp entry for '%s' in '%s', ignoring this entry and keeping the one from '%s.rtp'\n",
                        rrtp[nrtp].resname, rrdb, rrtp[firstrtp].filebase);
                /* We should free all the data for this entry.
                 * The current code gives a lot of dangling pointers.
                 */
                clear_t_restp(&rrtp[nrtp]);
            }
            else
            {
                gmx_fatal(FARGS, "Found rtp entries for '%s' in both '%s' and '%s'. If you want the first definition to override the second one, set the -rtpo option of pdb2gmx.", rrtp[nrtp].resname, rrtp[firstrtp].filebase, rrdb);
            }
        }
    }
    gmx_ffclose(in);
    /* give back unused memory */
    srenew(rrtp, nrtp);

    fprintf(stderr, "\nSorting it all out...\n");
    qsort(rrtp, nrtp, (size_t)sizeof(rrtp[0]), comprtp);

    check_rtp(nrtp, rrtp, rrdb);

    *nrtpptr = nrtp;
    *rtp     = rrtp;
}

/************************************************************
 *
 *                  SEARCH   ROUTINES
 *
 ***********************************************************/
static gmx_bool is_sign(char c)
{
    return (c == '+' || c == '-');
}

/* Compares if the strins match except for a sign at the end
 * of (only) one of the two.
 */
static int neq_str_sign(const char *a1, const char *a2)
{
    int l1, l2, lm;

    l1 = (int)strlen(a1);
    l2 = (int)strlen(a2);
    lm = min(l1, l2);

    if (lm >= 1 &&
        ((l1 == l2+1 && is_sign(a1[l1-1])) ||
         (l2 == l1+1 && is_sign(a2[l2-1]))) &&
        gmx_strncasecmp(a1, a2, lm) == 0)
    {
        return lm;
    }
    else
    {
        return 0;
    }
}

char *search_rtp(const char *key, int nrtp, t_restp rtp[])
{
    int  i, n, nbest, best, besti;
    char bestbuf[STRLEN];

    nbest =  0;
    besti = -1;
    /* We want to match at least one character */
    best  =  1;
    for (i = 0; (i < nrtp); i++)
    {
        if (gmx_strcasecmp(key, rtp[i].resname) == 0)
        {
            besti = i;
            nbest = 1;
            break;
        }
        else
        {
            /* Allow a mismatch of at most a sign character (with warning) */
            n = neq_str_sign(key, rtp[i].resname);
            if (n >= best &&
                n+1 >= (int)strlen(key) &&
                n+1 >= (int)strlen(rtp[i].resname))
            {
                if (n == best)
                {
                    if (nbest == 1)
                    {
                        strcpy(bestbuf, rtp[besti].resname);
                    }
                    if (nbest >= 1)
                    {
                        strcat(bestbuf, " or ");
                        strcat(bestbuf, rtp[i].resname);
                    }
                }
                else
                {
                    nbest = 0;
                }
                besti = i;
                best  = n;
                nbest++;
            }
        }
    }
    if (nbest > 1)
    {
        gmx_fatal(FARGS, "Residue '%s' not found in residue topology database, looks a bit like %s", key, bestbuf);
    }
    else if (besti == -1)
    {
        gmx_fatal(FARGS, "Residue '%s' not found in residue topology database", key);
    }
    if (gmx_strcasecmp(rtp[besti].resname, key) != 0)
    {
        fprintf(stderr,
                "\nWARNING: '%s' not found in residue topology database, "
                "trying to use '%s'\n\n", key, rtp[besti].resname);
    }

    return rtp[besti].resname;
}

t_restp *get_restp(const char *rtpname, int nrtp, t_restp rtp[])
{
    int  i;

    i = 0;
    while (i < nrtp && gmx_strcasecmp(rtpname, rtp[i].resname) != 0)
    {
        i++;
    }
    if (i >= nrtp)
    {
        /* This should never happen, since search_rtp should have been called
         * before calling get_restp.
         */
        gmx_fatal(FARGS, "Residue type '%s' not found in residue topology database", rtpname);
    }

    return &rtp[i];
}
