/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017, by the GROMACS development team, led by
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

#include "ter_db.h"

#include <ctype.h>
#include <string.h>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/gmxpreprocess/fflibutil.h"
#include "gromacs/gmxpreprocess/h_db.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/resall.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strdb.h"

/* use bonded types definitions in hackblock.h */
#define ekwRepl ebtsNR+1
#define ekwAdd  ebtsNR+2
#define ekwDel  ebtsNR+3
#define ekwNR   3
const char *kw_names[ekwNR] = {
    "replace", "add", "delete"
};

static int find_kw(char *keyw)
{
    int i;

    for (i = 0; i < ebtsNR; i++)
    {
        if (gmx_strcasecmp(btsNames[i], keyw) == 0)
        {
            return i;
        }
    }
    for (i = 0; i < ekwNR; i++)
    {
        if (gmx_strcasecmp(kw_names[i], keyw) == 0)
        {
            return ebtsNR + 1 + i;
        }
    }

    return NOTSET;
}

#define FATAL() gmx_fatal(FARGS, "Reading Termini Database: not enough items on line\n%s", line)

static void read_atom(char *line, gmx_bool bAdd,
                      char **nname, t_atom *a, gpp_atomtype_t atype, int *cgnr)
{
    int    nr, i;
    char   buf[5][30];
    double m, q;

    /* This code is messy, because of support for different formats:
     * for replace: [new name] <atom type> <m> <q> [cgnr (old format)]
     * for add:                <atom type> <m> <q> [cgnr]
     */
    nr = sscanf(line, "%s %s %s %s %s", buf[0], buf[1], buf[2], buf[3], buf[4]);

    /* Here there an ambiguity due to the old replace format with cgnr,
     * which was read for years, but ignored in the rest of the code.
     * We have to assume that the atom type does not start with a digit
     * to make a line with 4 entries uniquely interpretable.
     */
    if (!bAdd && nr == 4 && isdigit(buf[1][0]))
    {
        nr = 3;
    }

    if (nr < 3 || nr > 4)
    {
        gmx_fatal(FARGS, "Reading Termini Database: expected %d or %d items of atom data in stead of %d on line\n%s", 3, 4, nr, line);
    }
    i = 0;
    if (!bAdd)
    {
        if (nr == 4)
        {
            *nname = gmx_strdup(buf[i++]);
        }
        else
        {
            *nname = nullptr;
        }
    }
    a->type = get_atomtype_type(buf[i++], atype);
    sscanf(buf[i++], "%lf", &m);
    a->m = m;
    sscanf(buf[i++], "%lf", &q);
    a->q = q;
    if (bAdd && nr == 4)
    {
        sscanf(buf[i++], "%d", cgnr);
    }
    else
    {
        *cgnr = NOTSET;
    }
}

static void print_atom(FILE *out, t_atom *a, gpp_atomtype_t atype)
{
    fprintf(out, "\t%s\t%g\t%g\n",
            get_atomtype_name(a->type, atype), a->m, a->q);
}

static void print_ter_db(const char *ff, char C, int nb, t_hackblock tb[],
                         gpp_atomtype_t atype)
{
    FILE *out;
    int   i, j, k, bt, nrepl, nadd, ndel;
    char  buf[STRLEN];

    sprintf(buf, "%s-%c.tdb", ff, C);
    out = gmx_fio_fopen(buf, "w");

    for (i = 0; (i < nb); i++)
    {
        fprintf(out, "[ %s ]\n", tb[i].name);

        /* first count: */
        nrepl = 0;
        nadd  = 0;
        ndel  = 0;
        for (j = 0; j < tb[i].nhack; j++)
        {
            if (tb[i].hack[j].oname != nullptr && tb[i].hack[j].nname != nullptr)
            {
                nrepl++;
            }
            else if (tb[i].hack[j].oname == nullptr && tb[i].hack[j].nname != nullptr)
            {
                nadd++;
            }
            else if (tb[i].hack[j].oname != nullptr && tb[i].hack[j].nname == nullptr)
            {
                ndel++;
            }
            else if (tb[i].hack[j].oname == nullptr && tb[i].hack[j].nname == nullptr)
            {
                gmx_fatal(FARGS, "invalid hack (%s) in termini database", tb[i].name);
            }
        }
        if (nrepl)
        {
            fprintf(out, "[ %s ]\n", kw_names[ekwRepl-ebtsNR-1]);
            for (j = 0; j < tb[i].nhack; j++)
            {
                if (tb[i].hack[j].oname != nullptr && tb[i].hack[j].nname != nullptr)
                {
                    fprintf(out, "%s\t", tb[i].hack[j].oname);
                    print_atom(out, tb[i].hack[j].atom, atype);
                }
            }
        }
        if (nadd)
        {
            fprintf(out, "[ %s ]\n", kw_names[ekwAdd-ebtsNR-1]);
            for (j = 0; j < tb[i].nhack; j++)
            {
                if (tb[i].hack[j].oname == nullptr && tb[i].hack[j].nname != nullptr)
                {
                    print_ab(out, &(tb[i].hack[j]), tb[i].hack[j].nname);
                    print_atom(out, tb[i].hack[j].atom, atype);
                }
            }
        }
        if (ndel)
        {
            fprintf(out, "[ %s ]\n", kw_names[ekwDel-ebtsNR-1]);
            for (j = 0; j < tb[i].nhack; j++)
            {
                if (tb[i].hack[j].oname != nullptr && tb[i].hack[j].nname == nullptr)
                {
                    fprintf(out, "%s\n", tb[i].hack[j].oname);
                }
            }
        }
        for (bt = 0; bt < ebtsNR; bt++)
        {
            if (tb[i].rb[bt].nb)
            {
                fprintf(out, "[ %s ]\n", btsNames[bt]);
                for (j = 0; j < tb[i].rb[bt].nb; j++)
                {
                    for (k = 0; k < btsNiatoms[bt]; k++)
                    {
                        fprintf(out, "%s%s", k ? "\t" : "", tb[i].rb[bt].b[j].a[k]);
                    }
                    if (tb[i].rb[bt].b[j].s)
                    {
                        fprintf(out, "\t%s", tb[i].rb[bt].b[j].s);
                    }
                    fprintf(out, "\n");
                }
            }
        }
        fprintf(out, "\n");
    }
    gmx_fio_fclose(out);
}

static void read_ter_db_file(char *fn,
                             int *ntbptr, t_hackblock **tbptr,
                             gpp_atomtype_t atype)
{
    char         filebase[STRLEN], *ptr;
    FILE        *in;
    char         header[STRLEN], buf[STRLEN], line[STRLEN];
    t_hackblock *tb;
    int          i, j, n, ni, kwnr, nb, maxnb, nh;

    fflib_filename_base(fn, filebase, STRLEN);
    /* Remove the C/N termini extension */
    ptr = strrchr(filebase, '.');
    if (ptr != nullptr)
    {
        ptr[0] = '\0';
    }

    in = fflib_open(fn);
    if (debug)
    {
        fprintf(debug, "Opened %s\n", fn);
    }

    tb    = *tbptr;
    nb    = *ntbptr - 1;
    maxnb = 0;
    kwnr  = NOTSET;
    get_a_line(in, line, STRLEN);
    while (!feof(in))
    {
        if (get_header(line, header))
        {
            /* this is a new block, or a new keyword */
            kwnr = find_kw(header);

            if (kwnr == NOTSET)
            {
                nb++;
                /* here starts a new block */
                if (nb >= maxnb)
                {
                    maxnb = nb + 100;
                    srenew(tb, maxnb);
                }
                clear_t_hackblock(&tb[nb]);
                tb[nb].name     = gmx_strdup(header);
                tb[nb].filebase = gmx_strdup(filebase);
            }
        }
        else
        {
            if (nb < 0)
            {
                gmx_fatal(FARGS, "reading termini database: "
                          "directive expected before line:\n%s\n"
                          "This might be a file in an old format.", line);
            }
            /* this is not a header, so it must be data */
            if (kwnr >= ebtsNR)
            {
                /* this is a hack: add/rename/delete atoms */
                /* make space for hacks */
                if (tb[nb].nhack >= tb[nb].maxhack)
                {
                    tb[nb].maxhack += 10;
                    srenew(tb[nb].hack, tb[nb].maxhack);
                }
                nh = tb[nb].nhack;
                clear_t_hack(&(tb[nb].hack[nh]));
                for (i = 0; i < 4; i++)
                {
                    tb[nb].hack[nh].a[i] = nullptr;
                }
                tb[nb].nhack++;

                /* get data */
                n = 0;
                if (kwnr == ekwRepl || kwnr == ekwDel)
                {
                    if (sscanf(line, "%s%n", buf, &n) != 1)
                    {
                        gmx_fatal(FARGS, "Reading Termini Database '%s': "
                                  "expected atom name on line\n%s", fn, line);
                    }
                    tb[nb].hack[nh].oname = gmx_strdup(buf);
                    /* we only replace or delete one atom at a time */
                    tb[nb].hack[nh].nr = 1;
                }
                else if (kwnr == ekwAdd)
                {
                    read_ab(line, fn, &(tb[nb].hack[nh]));
                    get_a_line(in, line, STRLEN);
                }
                else
                {
                    gmx_fatal(FARGS, "unimplemented keyword number %d (%s:%d)",
                              kwnr, __FILE__, __LINE__);
                }
                if (kwnr == ekwRepl || kwnr == ekwAdd)
                {
                    snew(tb[nb].hack[nh].atom, 1);
                    read_atom(line+n, kwnr == ekwAdd,
                              &tb[nb].hack[nh].nname, tb[nb].hack[nh].atom, atype,
                              &tb[nb].hack[nh].cgnr);
                    if (tb[nb].hack[nh].nname == nullptr)
                    {
                        if (tb[nb].hack[nh].oname != nullptr)
                        {
                            tb[nb].hack[nh].nname = gmx_strdup(tb[nb].hack[nh].oname);
                        }
                        else
                        {
                            gmx_fatal(FARGS, "Reading Termini Database '%s': don't know which name the new atom should have on line\n%s", fn, line);
                        }
                    }
                }
            }
            else if (kwnr >= 0 && kwnr < ebtsNR)
            {
                /* this is bonded data: bonds, angles, dihedrals or impropers */
                srenew(tb[nb].rb[kwnr].b, tb[nb].rb[kwnr].nb+1);
                n = 0;
                for (j = 0; j < btsNiatoms[kwnr]; j++)
                {
                    if (sscanf(line+n, "%s%n", buf, &ni) == 1)
                    {
                        tb[nb].rb[kwnr].b[tb[nb].rb[kwnr].nb].a[j] = gmx_strdup(buf);
                    }
                    else
                    {
                        gmx_fatal(FARGS, "Reading Termini Database '%s': expected %d atom names (found %d) on line\n%s", fn, btsNiatoms[kwnr], j-1, line);
                    }
                    n += ni;
                }
                for (; j < MAXATOMLIST; j++)
                {
                    tb[nb].rb[kwnr].b[tb[nb].rb[kwnr].nb].a[j] = nullptr;
                }
                strcpy(buf, "");
                sscanf(line+n, "%s", buf);
                tb[nb].rb[kwnr].b[tb[nb].rb[kwnr].nb].s = gmx_strdup(buf);
                tb[nb].rb[kwnr].nb++;
            }
            else
            {
                gmx_fatal(FARGS, "Reading Termini Database: Expecting a header at line\n"
                          "%s", line);
            }
        }
        get_a_line(in, line, STRLEN);
    }
    nb++;
    srenew(tb, nb);

    gmx_ffclose(in);

    *ntbptr = nb;
    *tbptr  = tb;
}

int read_ter_db(const char *ffdir, char ter,
                t_hackblock **tbptr, gpp_atomtype_t atype)
{
    char   ext[STRLEN];
    int    ntdbf, f;
    char **tdbf;
    int    ntb;

    sprintf(ext, ".%c.tdb", ter);

    /* Search for termini database files.
     * Do not generate an error when none are found.
     */
    ntdbf  = fflib_search_file_end(ffdir, ext, FALSE, &tdbf);
    ntb    = 0;
    *tbptr = nullptr;
    for (f = 0; f < ntdbf; f++)
    {
        read_ter_db_file(tdbf[f], &ntb, tbptr, atype);
        sfree(tdbf[f]);
    }
    sfree(tdbf);

    if (debug)
    {
        print_ter_db("new", ter, ntb, *tbptr, atype);
    }

    return ntb;
}

t_hackblock **filter_ter(int nrtp, t_restp rtp[],
                         int nb, t_hackblock tb[],
                         const char *resname,
                         const char *rtpname,
                         int *nret)
{
    // TODO Four years later, no force fields have ever used this, so decide status of this feature
    /* Since some force fields (e.g. OPLS) needs different
     * atomtypes for different residues there could be a lot
     * of entries in the databases for specific residues
     * (e.g. GLY-NH3+, SER-NH3+, ALA-NH3+).
     *
     * To reduce the database size, we assume that a terminus specifier liker
     *
     * [ GLY|SER|ALA-NH3+ ]
     *
     * would cover all of the three residue types above.
     * Wildcards (*,?) are not OK. Don't worry about e.g. GLU vs. GLUH since
     * pdb2gmx only uses the first 3 letters when calling this routine.
     *
     * To automate this, this routines scans a list of termini
     * for the residue name "resname" and returns an allocated list of
     * pointers to the termini that could be applied to the
     * residue in question. The variable pointed to by nret will
     * contain the number of valid pointers in the list.
     * Remember to free the list when you are done with it...
     */

    t_restp     *   restp;
    int             i, j, n, none_idx;
    gmx_bool        found;
    char           *rtpname_match, *s;
    t_hackblock   **list;

    rtpname_match = search_rtp(rtpname, nrtp, rtp);
    restp         = get_restp(rtpname_match, nrtp, rtp);

    n    = 0;
    list = nullptr;

    for (i = 0; i < nb; i++)
    {
        s     = tb[i].name;
        found = FALSE;
        do
        {
            /* The residue name should appear in a tdb file with the same base name
             * as the file containing the rtp entry.
             * This makes termini selection for different molecule types
             * much cleaner.
             */
            if (gmx_strcasecmp(restp->filebase, tb[i].filebase) == 0 &&
                gmx_strncasecmp(resname, s, 3) == 0)
            {
                found = TRUE;
                srenew(list, n+1);
                list[n] = &(tb[i]);
                n++;
            }
            else
            {
                /* advance to next |-separated field */
                s = strchr(s, '|');
                if (s != nullptr)
                {
                    s++;
                }
            }
        }
        while (!found && s != nullptr);
    }

    /* All residue-specific termini have been added. We might have to fall
     * back on generic termini, which are characterized by not having
     * '-' in the name prior to the last position (which indicates charge).
     * The [ None ] alternative is special since we don't want that
     * to be the default, so we put it last in the list we return.
     */
    none_idx = -1;
    for (i = 0; i < nb; i++)
    {
        s = tb[i].name;
        /* The residue name should appear in a tdb file with the same base name
         * as the file containing the rtp entry.
         * This makes termini selection for different molecule types
         * much cleaner.
         */
        if (gmx_strcasecmp(restp->filebase, tb[i].filebase) == 0)
        {
            if (!gmx_strcasecmp("None", s))
            {
                none_idx = i;
            }
            else
            {
                /* Time to see if there's a generic terminus that matches.
                   Is there a hyphen? */
                char *c = strchr(s, '-');

                /* A conjunction hyphen normally indicates a residue-specific
                   terminus, which is named like "GLY-COOH". A generic terminus
                   won't have a hyphen. */
                bool bFoundAnyHyphen = (c != nullptr);
                /* '-' as the last character indicates charge, so if that's
                   the only one found e.g. "COO-", then it was not a conjunction
                   hyphen, so this is a generic terminus */
                bool bOnlyFoundChargeHyphen = (bFoundAnyHyphen &&
                                               *(c+1) == '\0');
                /* Thus, "GLY-COO-" is not recognized as a generic terminus. */
                bool bFoundGenericTerminus = !bFoundAnyHyphen || bOnlyFoundChargeHyphen;
                if (bFoundGenericTerminus)
                {
                    /* Check that we haven't already added a residue-specific version
                     * of this terminus.
                     */
                    for (j = 0; j < n && strstr((*list[j]).name, s) == nullptr; j++)
                    {
                        ;
                    }
                    if (j == n)
                    {
                        srenew(list, n+1);
                        list[n] = &(tb[i]);
                        n++;
                    }
                }
            }
        }
    }
    if (none_idx >= 0)
    {
        srenew(list, n+1);
        list[n] = &(tb[none_idx]);
        n++;
    }

    *nret = n;
    return list;
}


t_hackblock *choose_ter(int nb, t_hackblock **tb, const char *title)
{
    int i, sel, ret;

    printf("%s\n", title);
    for (i = 0; (i < nb); i++)
    {
        bool bIsZwitterion = (0 == gmx_wcmatch("*ZWITTERION*", (*tb[i]).name));
        printf("%2d: %s%s\n", i, (*tb[i]).name,
               bIsZwitterion ? " (only use with zwitterions containing exactly one residue)" : "");
    }
    do
    {
        ret = fscanf(stdin, "%d", &sel);
    }
    while ((ret != 1) || (sel < 0) || (sel >= nb));

    return tb[sel];
}
