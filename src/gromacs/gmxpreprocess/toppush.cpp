/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "toppush.h"

#include <ctype.h>
#include <stdlib.h>

#include <cmath>
#include <cstring>

#include <algorithm>

#include "gromacs/fileio/warninp.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/gpp_bond_atomtype.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

void generate_nbparams(int comb, int ftype, t_params *plist, gpp_atomtype_t atype,
                       warninp_t wi)
{
    int   i, j, k = -1, nf;
    int   nr, nrfp;
    real  c, bi, bj, ci, cj, ci0, ci1, ci2, cj0, cj1, cj2;
    char  errbuf[STRLEN];

    /* Lean mean shortcuts */
    nr   = get_atomtype_ntypes(atype);
    nrfp = NRFP(ftype);
    snew(plist->param, nr*nr);
    plist->nr = nr*nr;

    /* Fill the matrix with force parameters */
    switch (ftype)
    {
        case F_LJ:
            switch (comb)
            {
                case eCOMB_GEOMETRIC:
                    /* Gromos rules */
                    for (i = k = 0; (i < nr); i++)
                    {
                        for (j = 0; (j < nr); j++, k++)
                        {
                            for (nf = 0; (nf < nrfp); nf++)
                            {
                                ci = get_atomtype_nbparam(i, nf, atype);
                                cj = get_atomtype_nbparam(j, nf, atype);
                                c  = std::sqrt(ci * cj);
                                plist->param[k].c[nf]      = c;
                            }
                        }
                    }
                    break;

                case eCOMB_ARITHMETIC:
                    /* c0 and c1 are sigma and epsilon */
                    for (i = k = 0; (i < nr); i++)
                    {
                        for (j = 0; (j < nr); j++, k++)
                        {
                            ci0                  = get_atomtype_nbparam(i, 0, atype);
                            cj0                  = get_atomtype_nbparam(j, 0, atype);
                            ci1                  = get_atomtype_nbparam(i, 1, atype);
                            cj1                  = get_atomtype_nbparam(j, 1, atype);
                            plist->param[k].c[0] = (fabs(ci0) + fabs(cj0))*0.5;
                            /* Negative sigma signals that c6 should be set to zero later,
                             * so we need to propagate that through the combination rules.
                             */
                            if (ci0 < 0 || cj0 < 0)
                            {
                                plist->param[k].c[0] *= -1;
                            }
                            plist->param[k].c[1] = std::sqrt(ci1*cj1);
                        }
                    }

                    break;
                case eCOMB_GEOM_SIG_EPS:
                    /* c0 and c1 are sigma and epsilon */
                    for (i = k = 0; (i < nr); i++)
                    {
                        for (j = 0; (j < nr); j++, k++)
                        {
                            ci0                  = get_atomtype_nbparam(i, 0, atype);
                            cj0                  = get_atomtype_nbparam(j, 0, atype);
                            ci1                  = get_atomtype_nbparam(i, 1, atype);
                            cj1                  = get_atomtype_nbparam(j, 1, atype);
                            plist->param[k].c[0] = std::sqrt(fabs(ci0*cj0));
                            /* Negative sigma signals that c6 should be set to zero later,
                             * so we need to propagate that through the combination rules.
                             */
                            if (ci0 < 0 || cj0 < 0)
                            {
                                plist->param[k].c[0] *= -1;
                            }
                            plist->param[k].c[1] = std::sqrt(ci1*cj1);
                        }
                    }

                    break;
                default:
                    sprintf(errbuf, "No such combination rule %d", comb);
                    warning_error_and_exit(wi, errbuf, FARGS);
            }
            if (plist->nr != k)
            {
                gmx_incons("Topology processing, generate nb parameters");
            }
            break;

        case F_BHAM:
            /* Buckingham rules */
            for (i = k = 0; (i < nr); i++)
            {
                for (j = 0; (j < nr); j++, k++)
                {
                    ci0                  = get_atomtype_nbparam(i, 0, atype);
                    cj0                  = get_atomtype_nbparam(j, 0, atype);
                    ci2                  = get_atomtype_nbparam(i, 2, atype);
                    cj2                  = get_atomtype_nbparam(j, 2, atype);
                    bi                   = get_atomtype_nbparam(i, 1, atype);
                    bj                   = get_atomtype_nbparam(j, 1, atype);
                    plist->param[k].c[0] = std::sqrt(ci0 * cj0);
                    if ((bi == 0) || (bj == 0))
                    {
                        plist->param[k].c[1] = 0;
                    }
                    else
                    {
                        plist->param[k].c[1] = 2.0/(1/bi+1/bj);
                    }
                    plist->param[k].c[2] = std::sqrt(ci2 * cj2);
                }
            }

            break;
        default:
            sprintf(errbuf, "Invalid nonbonded type %s",
                    interaction_function[ftype].longname);
            warning_error(wi, errbuf);
    }
}

static void realloc_nb_params(gpp_atomtype_t at,
                              t_nbparam ***nbparam, t_nbparam ***pair)
{
    /* Add space in the non-bonded parameters matrix */
    int atnr = get_atomtype_ntypes(at);
    srenew(*nbparam, atnr);
    snew((*nbparam)[atnr-1], atnr);
    if (pair)
    {
        srenew(*pair, atnr);
        snew((*pair)[atnr-1], atnr);
    }
}

static void copy_B_from_A(int ftype, double *c)
{
    int nrfpA, nrfpB, i;

    nrfpA = NRFPA(ftype);
    nrfpB = NRFPB(ftype);

    /* Copy the B parameters from the first nrfpB A parameters */
    for (i = 0; (i < nrfpB); i++)
    {
        c[nrfpA+i] = c[i];
    }
}

void push_at (t_symtab *symtab, gpp_atomtype_t at, t_bond_atomtype bat,
              char *line, int nb_funct,
              t_nbparam ***nbparam, t_nbparam ***pair,
              warninp_t wi)
{
    typedef struct {
        const char *entry;
        int         ptype;
    } t_xlate;
    t_xlate    xl[eptNR] = {
        { "A",   eptAtom },
        { "N",   eptNucleus },
        { "S",   eptShell },
        { "B",   eptBond },
        { "V",   eptVSite },
    };

    int        nr, i, nfields, j, pt, nfp0 = -1;
    int        batype_nr, nread;
    char       type[STRLEN], btype[STRLEN], ptype[STRLEN];
    double     m, q;
    double     c[MAXFORCEPARAM];
    double     radius, vol, surftens, gb_radius, S_hct;
    char       tmpfield[12][100]; /* Max 12 fields of width 100 */
    char       errbuf[STRLEN];
    t_atom    *atom;
    t_param   *param;
    int        atomnr;
    gmx_bool   have_atomic_number;
    gmx_bool   have_bonded_type;

    snew(atom, 1);
    snew(param, 1);

    /* First assign input line to temporary array */
    nfields = sscanf(line, "%s%s%s%s%s%s%s%s%s%s%s%s",
                     tmpfield[0], tmpfield[1], tmpfield[2], tmpfield[3], tmpfield[4], tmpfield[5],
                     tmpfield[6], tmpfield[7], tmpfield[8], tmpfield[9], tmpfield[10], tmpfield[11]);

    /* Comments on optional fields in the atomtypes section:
     *
     * The force field format is getting a bit old. For OPLS-AA we needed
     * to add a special bonded atomtype, and for Gerrit Groenhofs QM/MM stuff
     * we also needed the atomic numbers.
     * To avoid making all old or user-generated force fields unusable we
     * have introduced both these quantities as optional columns, and do some
     * acrobatics to check whether they are present or not.
     * This will all look much nicer when we switch to XML... sigh.
     *
     * Field 0 (mandatory) is the nonbonded type name. (string)
     * Field 1 (optional)  is the bonded type (string)
     * Field 2 (optional)  is the atomic number (int)
     * Field 3 (mandatory) is the mass (numerical)
     * Field 4 (mandatory) is the charge (numerical)
     * Field 5 (mandatory) is the particle type (single character)
     * This is followed by a number of nonbonded parameters.
     *
     * The safest way to identify the format is the particle type field.
     *
     * So, here is what we do:
     *
     * A. Read in the first six fields as strings
     * B. If field 3 (starting from 0) is a single char, we have neither
     *    bonded_type or atomic numbers.
     * C. If field 5 is a single char we have both.
     * D. If field 4 is a single char we check field 1. If this begins with
     *    an alphabetical character we have bonded types, otherwise atomic numbers.
     */

    if (nfields < 6)
    {
        too_few(wi);
        return;
    }

    if ( (strlen(tmpfield[5]) == 1) && isalpha(tmpfield[5][0]) )
    {
        have_bonded_type   = TRUE;
        have_atomic_number = TRUE;
    }
    else if ( (strlen(tmpfield[3]) == 1) && isalpha(tmpfield[3][0]) )
    {
        have_bonded_type   = FALSE;
        have_atomic_number = FALSE;
    }
    else
    {
        have_bonded_type   = ( isalpha(tmpfield[1][0]) != 0 );
        have_atomic_number = !have_bonded_type;
    }

    /* optional fields */
    surftens  = -1;
    vol       = -1;
    radius    = -1;
    gb_radius = -1;
    atomnr    = -1;
    S_hct     = -1;

    switch (nb_funct)
    {

        case F_LJ:
            nfp0 = 2;

            if (have_atomic_number)
            {
                if (have_bonded_type)
                {
                    nread = sscanf(line, "%s%s%d%lf%lf%s%lf%lf%lf%lf%lf%lf",
                                   type, btype, &atomnr, &m, &q, ptype, &c[0], &c[1],
                                   &radius, &vol, &surftens, &gb_radius);
                    if (nread < 8)
                    {
                        too_few(wi);
                        return;
                    }
                }
                else
                {
                    /* have_atomic_number && !have_bonded_type */
                    nread = sscanf(line, "%s%d%lf%lf%s%lf%lf%lf%lf%lf%lf",
                                   type, &atomnr, &m, &q, ptype, &c[0], &c[1],
                                   &radius, &vol, &surftens, &gb_radius);
                    if (nread < 7)
                    {
                        too_few(wi);
                        return;
                    }
                }
            }
            else
            {
                if (have_bonded_type)
                {
                    /* !have_atomic_number && have_bonded_type */
                    nread = sscanf(line, "%s%s%lf%lf%s%lf%lf%lf%lf%lf%lf",
                                   type, btype, &m, &q, ptype, &c[0], &c[1],
                                   &radius, &vol, &surftens, &gb_radius);
                    if (nread < 7)
                    {
                        too_few(wi);
                        return;
                    }
                }
                else
                {
                    /* !have_atomic_number && !have_bonded_type */
                    nread = sscanf(line, "%s%lf%lf%s%lf%lf%lf%lf%lf%lf",
                                   type, &m, &q, ptype, &c[0], &c[1],
                                   &radius, &vol, &surftens, &gb_radius);
                    if (nread < 6)
                    {
                        too_few(wi);
                        return;
                    }
                }
            }

            if (!have_bonded_type)
            {
                strcpy(btype, type);
            }

            if (!have_atomic_number)
            {
                atomnr = -1;
            }

            break;

        case F_BHAM:
            nfp0 = 3;

            if (have_atomic_number)
            {
                if (have_bonded_type)
                {
                    nread = sscanf(line, "%s%s%d%lf%lf%s%lf%lf%lf%lf%lf%lf%lf",
                                   type, btype, &atomnr, &m, &q, ptype, &c[0], &c[1], &c[2],
                                   &radius, &vol, &surftens, &gb_radius);
                    if (nread < 9)
                    {
                        too_few(wi);
                        return;
                    }
                }
                else
                {
                    /* have_atomic_number && !have_bonded_type */
                    nread = sscanf(line, "%s%d%lf%lf%s%lf%lf%lf%lf%lf%lf%lf",
                                   type, &atomnr, &m, &q, ptype, &c[0], &c[1], &c[2],
                                   &radius, &vol, &surftens, &gb_radius);
                    if (nread < 8)
                    {
                        too_few(wi);
                        return;
                    }
                }
            }
            else
            {
                if (have_bonded_type)
                {
                    /* !have_atomic_number && have_bonded_type */
                    nread = sscanf(line, "%s%s%lf%lf%s%lf%lf%lf%lf%lf%lf%lf",
                                   type, btype, &m, &q, ptype, &c[0], &c[1], &c[2],
                                   &radius, &vol, &surftens, &gb_radius);
                    if (nread < 8)
                    {
                        too_few(wi);
                        return;
                    }
                }
                else
                {
                    /* !have_atomic_number && !have_bonded_type */
                    nread = sscanf(line, "%s%lf%lf%s%lf%lf%lf%lf%lf%lf%lf",
                                   type, &m, &q, ptype, &c[0], &c[1], &c[2],
                                   &radius, &vol, &surftens, &gb_radius);
                    if (nread < 7)
                    {
                        too_few(wi);
                        return;
                    }
                }
            }

            if (!have_bonded_type)
            {
                strcpy(btype, type);
            }

            if (!have_atomic_number)
            {
                atomnr = -1;
            }

            break;

        default:
            sprintf(errbuf, "Invalid function type %d in push_at", nb_funct);
            warning_error_and_exit(wi, errbuf, FARGS);
    }
    for (j = nfp0; (j < MAXFORCEPARAM); j++)
    {
        c[j] = 0.0;
    }

    if (strlen(type) == 1 && isdigit(type[0]))
    {
        warning_error_and_exit(wi, "Atom type names can't be single digits.", FARGS);
    }

    if (strlen(btype) == 1 && isdigit(btype[0]))
    {
        warning_error_and_exit(wi, "Bond atom type names can't be single digits.", FARGS);
    }

    /* Hack to read old topologies */
    if (gmx_strcasecmp(ptype, "D") == 0)
    {
        sprintf(ptype, "V");
    }
    for (j = 0; (j < eptNR); j++)
    {
        if (gmx_strcasecmp(ptype, xl[j].entry) == 0)
        {
            break;
        }
    }
    if (j == eptNR)
    {
        sprintf(errbuf, "Invalid particle type %s", ptype);
        warning_error_and_exit(wi, errbuf, FARGS);
    }
    pt = xl[j].ptype;
    if (debug)
    {
        fprintf(debug, "ptype: %s\n", ptype_str[pt]);
    }

    atom->q     = q;
    atom->m     = m;
    atom->ptype = pt;
    for (i = 0; (i < MAXFORCEPARAM); i++)
    {
        param->c[i] = c[i];
    }

    if ((batype_nr = get_bond_atomtype_type(btype, bat)) == NOTSET)
    {
        add_bond_atomtype(bat, symtab, btype);
    }
    batype_nr = get_bond_atomtype_type(btype, bat);

    if ((nr = get_atomtype_type(type, at)) != NOTSET)
    {
        sprintf(errbuf, "Overriding atomtype %s", type);
        warning(wi, errbuf);
        if ((nr = set_atomtype(nr, at, symtab, atom, type, param, batype_nr,
                               radius, vol, surftens, atomnr, gb_radius, S_hct)) == NOTSET)
        {
            sprintf(errbuf, "Replacing atomtype %s failed", type);
            warning_error_and_exit(wi, errbuf, FARGS);
        }
    }
    else if ((add_atomtype(at, symtab, atom, type, param,
                           batype_nr, radius, vol,
                           surftens, atomnr, gb_radius, S_hct)) == NOTSET)
    {
        sprintf(errbuf, "Adding atomtype %s failed", type);
        warning_error_and_exit(wi, errbuf, FARGS);
    }
    else
    {
        /* Add space in the non-bonded parameters matrix */
        realloc_nb_params(at, nbparam, pair);
    }
    sfree(atom);
    sfree(param);
}

//! Return whether the contents of \c a and \c b are the same, considering also reversed order.
template <typename T>
static bool equalEitherForwardOrBackward(gmx::ArrayRef<const T> a, gmx::ArrayRef<const T> b)
{
    return (std::equal(a.begin(), a.end(), b.begin()) ||
            std::equal(a.begin(), a.end(), b.rbegin()));
}

static void push_bondtype(t_params     *       bt,
                          const t_param *      b,
                          int                  nral,
                          int                  ftype,
                          gmx_bool             bAllowRepeat,
                          const char *         line,
                          warninp_t            wi)
{
    int      nr   = bt->nr;
    int      nrfp = NRFP(ftype);
    char     errbuf[STRLEN];

    /* If bAllowRepeat is TRUE, we allow multiple entries as long as they
       are on directly _adjacent_ lines.
     */

    /* First check if our atomtypes are _identical_ (not reversed) to the previous
       entry. If they are not identical we search for earlier duplicates. If they are
       we can skip it, since we already searched for the first line
       in this group.
     */

    bool isContinuationOfBlock = false;
    if (bAllowRepeat && nr > 1)
    {
        isContinuationOfBlock = true;
        for (int j = 0; j < nral; j++)
        {
            if (b->a[j] != bt->param[nr - 2].a[j])
            {
                isContinuationOfBlock = false;
            }
        }
    }

    /* Search for earlier duplicates if this entry was not a continuation
       from the previous line.
     */
    bool addBondType = true;
    bool haveWarned  = false;
    bool haveErrored = false;
    for (int i = 0; (i < nr); i++)
    {
        gmx::ArrayRef<const int> bParams(b->a, b->a + nral);
        gmx::ArrayRef<const int> testParams(bt->param[i].a, bt->param[i].a + nral);
        if (equalEitherForwardOrBackward(bParams, testParams))
        {
            GMX_ASSERT(nrfp <= MAXFORCEPARAM, "This is ensured in other places, but we need this assert to keep the clang analyzer happy");
            // TODO consider improving the following code by using:
            // bool identicalParameters = std::equal(bt->param[i].c, bt->param[i].c + nrfp, b->c);
            bool identicalParameters = true;
            for (int j = 0; (j < nrfp); j++)
            {
                identicalParameters = identicalParameters && (bt->param[i].c[j] == b->c[j]);
            }

            if (!bAllowRepeat || identicalParameters)
            {
                addBondType = false;
            }

            if (!identicalParameters)
            {
                if (bAllowRepeat)
                {
                    /* With dihedral type 9 we only allow for repeating
                     * of the same parameters with blocks with 1 entry.
                     * Allowing overriding is too complex to check.
                     */
                    if (!isContinuationOfBlock && !haveErrored)
                    {
                        warning_error(wi, "Encountered a second block of parameters for dihedral type 9 for the same atoms, with either different parameters and/or the first block has multiple lines. This is not supported.");
                        haveErrored = true;
                    }
                }
                else if (!haveWarned)
                {
                    sprintf(errbuf, "Overriding %s parameters.%s",
                            interaction_function[ftype].longname,
                            (ftype == F_PDIHS) ?
                            "\nUse dihedraltype 9 to allow several multiplicity terms. Only consecutive lines are combined. Non-consective lines overwrite each other."
                            : "");
                    warning(wi, errbuf);

                    fprintf(stderr, "  old:                                         ");
                    for (int j = 0; (j < nrfp); j++)
                    {
                        fprintf(stderr, " %g", bt->param[i].c[j]);
                    }
                    fprintf(stderr, " \n  new: %s\n\n", line);

                    haveWarned = true;
                }
            }

            if (!identicalParameters && !bAllowRepeat)
            {
                /* Overwrite the parameters with the latest ones */
                // TODO considering improving the following code by replacing with:
                // std::copy(b->c, b->c + nrfp, bt->param[i].c);
                for (int j = 0; (j < nrfp); j++)
                {
                    bt->param[i].c[j] = b->c[j];
                }
            }
        }
    }

    if (addBondType)
    {
        /* alloc */
        pr_alloc (2, bt);

        /* fill the arrays up and down */
        memcpy(bt->param[bt->nr].c,  b->c, sizeof(b->c));
        memcpy(bt->param[bt->nr].a,  b->a, sizeof(b->a));
        memcpy(bt->param[bt->nr+1].c, b->c, sizeof(b->c));

        /* The definitions of linear angles depend on the order of atoms,
         * that means that for atoms i-j-k, with certain parameter a, the
         * corresponding k-j-i angle will have parameter 1-a.
         */
        if (ftype == F_LINEAR_ANGLES)
        {
            bt->param[bt->nr+1].c[0] = 1-bt->param[bt->nr+1].c[0];
            bt->param[bt->nr+1].c[2] = 1-bt->param[bt->nr+1].c[2];
        }

        for (int j = 0; (j < nral); j++)
        {
            bt->param[bt->nr+1].a[j] = b->a[nral-1-j];
        }

        bt->nr += 2;
    }
}

void push_bt(directive d, t_params bt[], int nral,
             gpp_atomtype_t at,
             t_bond_atomtype bat, char *line,
             warninp_t wi)
{
    const char *formal[MAXATOMLIST+1] = {
        "%s",
        "%s%s",
        "%s%s%s",
        "%s%s%s%s",
        "%s%s%s%s%s",
        "%s%s%s%s%s%s",
        "%s%s%s%s%s%s%s"
    };
    const char *formnl[MAXATOMLIST+1] = {
        "%*s",
        "%*s%*s",
        "%*s%*s%*s",
        "%*s%*s%*s%*s",
        "%*s%*s%*s%*s%*s",
        "%*s%*s%*s%*s%*s%*s",
        "%*s%*s%*s%*s%*s%*s%*s"
    };
    const char *formlf = "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf";
    int         i, ft, ftype, nn, nrfp, nrfpA;
    char        f1[STRLEN];
    char        alc[MAXATOMLIST+1][20];
    /* One force parameter more, so we can check if we read too many */
    double      c[MAXFORCEPARAM+1];
    t_param     p;
    char        errbuf[STRLEN];

    if ((bat && at) || (!bat && !at))
    {
        gmx_incons("You should pass either bat or at to push_bt");
    }

    /* Make format string (nral ints+functype) */
    if ((nn = sscanf(line, formal[nral],
                     alc[0], alc[1], alc[2], alc[3], alc[4], alc[5])) != nral+1)
    {
        sprintf(errbuf, "Not enough atomtypes (%d instead of %d)", nn-1, nral);
        warning_error(wi, errbuf);
        return;
    }

    ft    = strtol(alc[nral], nullptr, 10);
    ftype = ifunc_index(d, ft);
    nrfp  = NRFP(ftype);
    nrfpA = interaction_function[ftype].nrfpA;
    strcpy(f1, formnl[nral]);
    strcat(f1, formlf);
    if ((nn = sscanf(line, f1, &c[0], &c[1], &c[2], &c[3], &c[4], &c[5], &c[6], &c[7], &c[8], &c[9], &c[10], &c[11], &c[12]))
        != nrfp)
    {
        if (nn == nrfpA)
        {
            /* Copy the B-state from the A-state */
            copy_B_from_A(ftype, c);
        }
        else
        {
            if (nn < nrfpA)
            {
                warning_error(wi, "Not enough parameters");
            }
            else if (nn > nrfpA && nn < nrfp)
            {
                warning_error(wi, "Too many parameters or not enough parameters for topology B");
            }
            else if (nn > nrfp)
            {
                warning_error(wi, "Too many parameters");
            }
            for (i = nn; (i < nrfp); i++)
            {
                c[i] = 0.0;
            }
        }
    }
    for (i = 0; (i < nral); i++)
    {
        if (at && ((p.a[i] = get_atomtype_type(alc[i], at)) == NOTSET))
        {
            sprintf(errbuf, "Unknown atomtype %s\n", alc[i]);
            warning_error_and_exit(wi, errbuf, FARGS);
        }
        else if (bat && ((p.a[i] = get_bond_atomtype_type(alc[i], bat)) == NOTSET))
        {
            sprintf(errbuf, "Unknown bond_atomtype %s\n", alc[i]);
            warning_error_and_exit(wi, errbuf, FARGS);
        }
    }
    for (i = 0; (i < nrfp); i++)
    {
        p.c[i] = c[i];
    }
    push_bondtype (&(bt[ftype]), &p, nral, ftype, FALSE, line, wi);
}


void push_dihedraltype(directive d, t_params bt[],
                       t_bond_atomtype bat, char *line,
                       warninp_t wi)
{
    const char  *formal[MAXATOMLIST+1] = {
        "%s",
        "%s%s",
        "%s%s%s",
        "%s%s%s%s",
        "%s%s%s%s%s",
        "%s%s%s%s%s%s",
        "%s%s%s%s%s%s%s"
    };
    const char  *formnl[MAXATOMLIST+1] = {
        "%*s",
        "%*s%*s",
        "%*s%*s%*s",
        "%*s%*s%*s%*s",
        "%*s%*s%*s%*s%*s",
        "%*s%*s%*s%*s%*s%*s",
        "%*s%*s%*s%*s%*s%*s%*s"
    };
    const char  *formlf[MAXFORCEPARAM] = {
        "%lf",
        "%lf%lf",
        "%lf%lf%lf",
        "%lf%lf%lf%lf",
        "%lf%lf%lf%lf%lf",
        "%lf%lf%lf%lf%lf%lf",
        "%lf%lf%lf%lf%lf%lf%lf",
        "%lf%lf%lf%lf%lf%lf%lf%lf",
        "%lf%lf%lf%lf%lf%lf%lf%lf%lf",
        "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
        "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
        "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
    };
    int          i, ft, ftype, nn, nrfp, nrfpA, nral;
    char         f1[STRLEN];
    char         alc[MAXATOMLIST+1][20];
    double       c[MAXFORCEPARAM];
    t_param      p;
    gmx_bool     bAllowRepeat;
    char         errbuf[STRLEN];

    /* This routine accepts dihedraltypes defined from either 2 or 4 atoms.
     *
     * We first check for 2 atoms with the 3th column being an integer
     * defining the type. If this isn't the case, we try it with 4 atoms
     * and the 5th column defining the dihedral type.
     */
    nn = sscanf(line, formal[4], alc[0], alc[1], alc[2], alc[3], alc[4]);
    if (nn >= 3 && strlen(alc[2]) == 1 && isdigit(alc[2][0]))
    {
        nral  = 2;
        ft    = strtol(alc[nral], nullptr, 10);
        /* Move atom types around a bit and use 'X' for wildcard atoms
         * to create a 4-atom dihedral definition with arbitrary atoms in
         * position 1 and 4.
         */
        if (alc[2][0] == '2')
        {
            /* improper - the two atomtypes are 1,4. Use wildcards for 2,3 */
            strcpy(alc[3], alc[1]);
            sprintf(alc[2], "X");
            sprintf(alc[1], "X");
            /* alc[0] stays put */
        }
        else
        {
            /* proper - the two atomtypes are 2,3. Use wildcards for 1,4 */
            sprintf(alc[3], "X");
            strcpy(alc[2], alc[1]);
            strcpy(alc[1], alc[0]);
            sprintf(alc[0], "X");
        }
    }
    else if (nn == 5 && strlen(alc[4]) == 1 && isdigit(alc[4][0]))
    {
        nral  = 4;
        ft    = strtol(alc[nral], nullptr, 10);
    }
    else
    {
        sprintf(errbuf, "Incorrect number of atomtypes for dihedral (%d instead of 2 or 4)", nn-1);
        warning_error(wi, errbuf);
        return;
    }

    if (ft == 9)
    {
        /* Previously, we have always overwritten parameters if e.g. a torsion
           with the same atomtypes occurs on multiple lines. However, CHARMM and
           some other force fields specify multiple dihedrals over some bonds,
           including cosines with multiplicity 6 and somethimes even higher.
           Thus, they cannot be represented with Ryckaert-Bellemans terms.
           To add support for these force fields, Dihedral type 9 is identical to
           normal proper dihedrals, but repeated entries are allowed.
         */
        bAllowRepeat = TRUE;
        ft           = 1;
    }
    else
    {
        bAllowRepeat = FALSE;
    }


    ftype = ifunc_index(d, ft);
    nrfp  = NRFP(ftype);
    nrfpA = interaction_function[ftype].nrfpA;

    strcpy(f1, formnl[nral]);
    strcat(f1, formlf[nrfp-1]);

    /* Check number of parameters given */
    if ((nn = sscanf(line, f1, &c[0], &c[1], &c[2], &c[3], &c[4], &c[5], &c[6], &c[7], &c[8], &c[9], &c[10], &c[11]))
        != nrfp)
    {
        if (nn == nrfpA)
        {
            /* Copy the B-state from the A-state */
            copy_B_from_A(ftype, c);
        }
        else
        {
            if (nn < nrfpA)
            {
                warning_error(wi, "Not enough parameters");
            }
            else if (nn > nrfpA && nn < nrfp)
            {
                warning_error(wi, "Too many parameters or not enough parameters for topology B");
            }
            else if (nn > nrfp)
            {
                warning_error(wi, "Too many parameters");
            }
            for (i = nn; (i < nrfp); i++)
            {
                c[i] = 0.0;
            }
        }
    }

    for (i = 0; (i < 4); i++)
    {
        if (!strcmp(alc[i], "X"))
        {
            p.a[i] = -1;
        }
        else
        {
            if ((p.a[i] = get_bond_atomtype_type(alc[i], bat)) == NOTSET)
            {
                sprintf(errbuf, "Unknown bond_atomtype %s", alc[i]);
                warning_error_and_exit(wi, errbuf, FARGS);
            }
        }
    }
    for (i = 0; (i < nrfp); i++)
    {
        p.c[i] = c[i];
    }
    /* Always use 4 atoms here, since we created two wildcard atoms
     * if there wasn't of them 4 already.
     */
    push_bondtype (&(bt[ftype]), &p, 4, ftype, bAllowRepeat, line, wi);
}


void push_nbt(directive d, t_nbparam **nbt, gpp_atomtype_t atype,
              char *pline, int nb_funct,
              warninp_t wi)
{
    /* swap the atoms */
    const char *form3 = "%*s%*s%*s%lf%lf%lf";
    const char *form4 = "%*s%*s%*s%lf%lf%lf%lf";
    const char *form5 = "%*s%*s%*s%lf%lf%lf%lf%lf";
    char        a0[80], a1[80];
    int         i, f, n, ftype, nrfp;
    double      c[4], dum;
    real        cr[4];
    int         ai, aj;
    t_nbparam  *nbp;
    gmx_bool    bId;
    char        errbuf[STRLEN];

    if (sscanf (pline, "%s%s%d", a0, a1, &f) != 3)
    {
        too_few(wi);
        return;
    }

    ftype = ifunc_index(d, f);

    if (ftype != nb_funct)
    {
        sprintf(errbuf, "Trying to add %s while the default nonbond type is %s",
                interaction_function[ftype].longname,
                interaction_function[nb_funct].longname);
        warning_error(wi, errbuf);
        return;
    }

    /* Get the force parameters */
    nrfp = NRFP(ftype);
    if (ftype == F_LJ14)
    {
        n = sscanf(pline, form4, &c[0], &c[1], &c[2], &c[3]);
        if (n < 2)
        {
            too_few(wi);
            return;
        }
        /* When the B topology parameters are not set,
         * copy them from topology A
         */
        GMX_ASSERT(nrfp <= 4, "LJ-14 cannot have more than 4 parameters");
        for (i = n; i < nrfp; i++)
        {
            c[i] = c[i-2];
        }
    }
    else if (ftype == F_LJC14_Q)
    {
        n = sscanf(pline, form5, &c[0], &c[1], &c[2], &c[3], &dum);
        if (n != 4)
        {
            incorrect_n_param(wi);
            return;
        }
    }
    else if (nrfp == 2)
    {
        if (sscanf(pline, form3, &c[0], &c[1], &dum) != 2)
        {
            incorrect_n_param(wi);
            return;
        }
    }
    else if (nrfp == 3)
    {
        if (sscanf(pline, form4, &c[0], &c[1], &c[2], &dum) != 3)
        {
            incorrect_n_param(wi);
            return;
        }
    }
    else
    {
        sprintf(errbuf, "Number of force parameters for nonbonded interactions is %d", nrfp);
        warning_error_and_exit(wi, errbuf, FARGS);
    }
    for (i = 0; (i < nrfp); i++)
    {
        cr[i] = c[i];
    }

    /* Put the parameters in the matrix */
    if ((ai = get_atomtype_type (a0, atype)) == NOTSET)
    {
        sprintf(errbuf, "Atomtype %s not found", a0);
        warning_error_and_exit(wi, errbuf, FARGS);
    }
    if ((aj = get_atomtype_type (a1, atype)) == NOTSET)
    {
        sprintf(errbuf, "Atomtype %s not found", a1);
        warning_error_and_exit(wi, errbuf, FARGS);
    }
    nbp = &(nbt[std::max(ai, aj)][std::min(ai, aj)]);

    if (nbp->bSet)
    {
        bId = TRUE;
        for (i = 0; i < nrfp; i++)
        {
            bId = bId && (nbp->c[i] == cr[i]);
        }
        if (!bId)
        {
            sprintf(errbuf, "Overriding non-bonded parameters,");
            warning(wi, errbuf);
            fprintf(stderr, "  old:");
            for (i = 0; i < nrfp; i++)
            {
                fprintf(stderr, " %g", nbp->c[i]);
            }
            fprintf(stderr, " new\n%s\n", pline);
        }
    }
    nbp->bSet = TRUE;
    for (i = 0; i < nrfp; i++)
    {
        nbp->c[i] = cr[i];
    }
}

void
push_gb_params (gpp_atomtype_t at, char *line,
                warninp_t wi)
{
    int    atype;
    double radius, vol, surftens, gb_radius, S_hct;
    char   atypename[STRLEN];
    char   errbuf[STRLEN];

    if ( (sscanf(line, "%s%lf%lf%lf%lf%lf", atypename, &radius, &vol, &surftens, &gb_radius, &S_hct)) != 6)
    {
        sprintf(errbuf, "Too few gb parameters for type %s\n", atypename);
        warning(wi, errbuf);
    }

    /* Search for atomtype */
    atype = get_atomtype_type(atypename, at);

    if (atype == NOTSET)
    {
        printf("Couldn't find topology match for atomtype %s\n", atypename);
        abort();
    }

    set_atomtype_gbparam(at, atype, radius, vol, surftens, gb_radius, S_hct);
}

void
push_cmaptype(directive d, t_params bt[], int nral, gpp_atomtype_t at,
              t_bond_atomtype bat, char *line,
              warninp_t wi)
{
    const char  *formal = "%s%s%s%s%s%s%s%s%n";

    int          i, ft, ftype, nn, nrfp, nrfpA, nrfpB;
    int          start, nchar_consumed;
    int          nxcmap, nycmap, ncmap, read_cmap, sl, nct;
    char         s[20], alc[MAXATOMLIST+2][20];
    t_param      p;
    char         errbuf[STRLEN];

    /* Keep the compiler happy */
    read_cmap = 0;
    start     = 0;

    GMX_ASSERT(nral == 5, "CMAP requires 5 atoms per interaction");

    /* Here we can only check for < 8 */
    if ((nn = sscanf(line, formal, alc[0], alc[1], alc[2], alc[3], alc[4], alc[5], alc[6], alc[7], &nchar_consumed)) < nral+3)
    {
        sprintf(errbuf, "Incorrect number of atomtypes for cmap (%d instead of 5)", nn-1);
        warning_error(wi, errbuf);
        return;
    }
    start += nchar_consumed;

    ft     = strtol(alc[nral], nullptr, 10);
    nxcmap = strtol(alc[nral+1], nullptr, 10);
    nycmap = strtol(alc[nral+2], nullptr, 10);

    /* Check for equal grid spacing in x and y dims */
    if (nxcmap != nycmap)
    {
        sprintf(errbuf, "Not the same grid spacing in x and y for cmap grid: x=%d, y=%d", nxcmap, nycmap);
        warning_error(wi, errbuf);
    }

    ncmap  = nxcmap*nycmap;
    ftype  = ifunc_index(d, ft);
    nrfpA  = strtol(alc[6], nullptr, 10)*strtol(alc[6], nullptr, 10);
    nrfpB  = strtol(alc[7], nullptr, 10)*strtol(alc[7], nullptr, 10);
    nrfp   = nrfpA+nrfpB;

    /* Allocate memory for the CMAP grid */
    bt[F_CMAP].ncmap += nrfp;
    srenew(bt[F_CMAP].cmap, bt[F_CMAP].ncmap);

    /* Read in CMAP parameters */
    sl = 0;
    for (i = 0; i < ncmap; i++)
    {
        while (isspace(*(line+start+sl)))
        {
            sl++;
        }
        nn  = sscanf(line+start+sl, " %s ", s);
        sl += strlen(s);
        bt[F_CMAP].cmap[i+(bt[F_CMAP].ncmap)-nrfp] = strtod(s, nullptr);

        if (nn == 1)
        {
            read_cmap++;
        }
        else
        {
            sprintf(errbuf, "Error in reading cmap parameter for angle %s %s %s %s %s", alc[0], alc[1], alc[2], alc[3], alc[4]);
            warning_error(wi, errbuf);
        }

    }

    /* Check do that we got the number of parameters we expected */
    if (read_cmap == nrfpA)
    {
        for (i = 0; i < ncmap; i++)
        {
            bt[F_CMAP].cmap[i+ncmap] = bt[F_CMAP].cmap[i];
        }
    }
    else
    {
        if (read_cmap < nrfpA)
        {
            warning_error(wi, "Not enough cmap parameters");
        }
        else if (read_cmap > nrfpA && read_cmap < nrfp)
        {
            warning_error(wi, "Too many cmap parameters or not enough parameters for topology B");
        }
        else if (read_cmap > nrfp)
        {
            warning_error(wi, "Too many cmap parameters");
        }
    }


    /* Set grid spacing and the number of grids (we assume these numbers to be the same for all grids
     * so we can safely assign them each time
     */
    bt[F_CMAP].grid_spacing = nxcmap;            /* Or nycmap, they need to be equal */
    bt[F_CMAP].nc           = bt[F_CMAP].nc + 1; /* Since we are incrementing here, we need to subtract later, see (*****) */
    nct                     = (nral+1) * bt[F_CMAP].nc;

    /* Allocate memory for the cmap_types information */
    srenew(bt[F_CMAP].cmap_types, nct);

    for (i = 0; (i < nral); i++)
    {
        if (at && ((p.a[i] = get_bond_atomtype_type(alc[i], bat)) == NOTSET))
        {
            sprintf(errbuf, "Unknown atomtype %s\n", alc[i]);
            warning_error(wi, errbuf);
        }
        else if (bat && ((p.a[i] = get_bond_atomtype_type(alc[i], bat)) == NOTSET))
        {
            sprintf(errbuf, "Unknown bond_atomtype %s\n", alc[i]);
            warning_error(wi, errbuf);
        }

        /* Assign a grid number to each cmap_type */
        bt[F_CMAP].cmap_types[bt[F_CMAP].nct++] = get_bond_atomtype_type(alc[i], bat);
    }

    /* Assign a type number to this cmap */
    bt[F_CMAP].cmap_types[bt[F_CMAP].nct++] = bt[F_CMAP].nc-1; /* Since we inremented earlier, we need to subtrac here, to get the types right (****) */

    /* Check for the correct number of atoms (again) */
    if (bt[F_CMAP].nct != nct)
    {
        sprintf(errbuf, "Incorrect number of atom types (%d) in cmap type %d\n", nct, bt[F_CMAP].nc);
        warning_error(wi, errbuf);
    }

    /* Is this correct?? */
    for (i = 0; i < MAXFORCEPARAM; i++)
    {
        p.c[i] = NOTSET;
    }

    /* Push the bond to the bondlist */
    push_bondtype (&(bt[ftype]), &p, nral, ftype, FALSE, line, wi);
}


static void push_atom_now(t_symtab *symtab, t_atoms *at, int atomnr,
                          int atomicnumber,
                          int type, char *ctype, int ptype,
                          char *resnumberic,
                          char *resname, char *name, real m0, real q0,
                          int typeB, char *ctypeB, real mB, real qB,
                          warninp_t wi)
{
    int           j, resind = 0, resnr;
    unsigned char ric;
    int           nr = at->nr;
    char          errbuf[STRLEN];

    if (((nr == 0) && (atomnr != 1)) || (nr && (atomnr != at->nr+1)))
    {
        sprintf(errbuf, "Atoms in the .top are not numbered consecutively from 1 (rather, atomnr = %d, while at->nr = %d)", atomnr, at->nr);
        warning_error_and_exit(wi, errbuf, FARGS);
    }

    j = strlen(resnumberic) - 1;
    if (isdigit(resnumberic[j]))
    {
        ric = ' ';
    }
    else
    {
        ric = resnumberic[j];
        if (j == 0 || !isdigit(resnumberic[j-1]))
        {
            sprintf(errbuf, "Invalid residue number '%s' for atom %d",
                    resnumberic, atomnr);
            warning_error_and_exit(wi, errbuf, FARGS);
        }
    }
    resnr = strtol(resnumberic, nullptr, 10);

    if (nr > 0)
    {
        resind = at->atom[nr-1].resind;
    }
    if (nr == 0 || strcmp(resname, *at->resinfo[resind].name) != 0 ||
        resnr != at->resinfo[resind].nr ||
        ric   != at->resinfo[resind].ic)
    {
        if (nr == 0)
        {
            resind = 0;
        }
        else
        {
            resind++;
        }
        at->nres = resind + 1;
        srenew(at->resinfo, at->nres);
        at->resinfo[resind].name = put_symtab(symtab, resname);
        at->resinfo[resind].nr   = resnr;
        at->resinfo[resind].ic   = ric;
    }
    else
    {
        resind = at->atom[at->nr-1].resind;
    }

    /* New atom instance
     * get new space for arrays
     */
    srenew(at->atom, nr+1);
    srenew(at->atomname, nr+1);
    srenew(at->atomtype, nr+1);
    srenew(at->atomtypeB, nr+1);

    /* fill the list */
    at->atom[nr].type  = type;
    at->atom[nr].ptype = ptype;
    at->atom[nr].q     = q0;
    at->atom[nr].m     = m0;
    at->atom[nr].typeB = typeB;
    at->atom[nr].qB    = qB;
    at->atom[nr].mB    = mB;

    at->atom[nr].resind     = resind;
    at->atom[nr].atomnumber = atomicnumber;
    at->atomname[nr]        = put_symtab(symtab, name);
    at->atomtype[nr]        = put_symtab(symtab, ctype);
    at->atomtypeB[nr]       = put_symtab(symtab, ctypeB);
    at->nr++;
}

static void push_cg(t_block *block, int *lastindex, int index, int a)
{
    if (debug)
    {
        fprintf (debug, "Index %d, Atom %d\n", index, a);
    }

    if (((block->nr) && (*lastindex != index)) || (!block->nr))
    {
        /* add a new block */
        block->nr++;
        srenew(block->index, block->nr+1);
    }
    block->index[block->nr] = a + 1;
    *lastindex              = index;
}

void push_atom(t_symtab *symtab, t_block *cgs,
               t_atoms *at, gpp_atomtype_t atype, char *line, int *lastcg,
               warninp_t wi)
{
    int           nr, ptype;
    int           cgnumber, atomnr, type, typeB, nscan;
    char          id[STRLEN], ctype[STRLEN], ctypeB[STRLEN],
                  resnumberic[STRLEN], resname[STRLEN], name[STRLEN], check[STRLEN];
    double        m, q, mb, qb;
    real          m0, q0, mB, qB;
    char          errbuf[STRLEN];

    /* Make a shortcut for writing in this molecule  */
    nr = at->nr;

    /* Fixed parameters */
    if (sscanf(line, "%s%s%s%s%s%d",
               id, ctype, resnumberic, resname, name, &cgnumber) != 6)
    {
        too_few(wi);
        return;
    }
    sscanf(id, "%d", &atomnr);
    if ((type  = get_atomtype_type(ctype, atype)) == NOTSET)
    {
        sprintf(errbuf, "Atomtype %s not found", ctype);
        warning_error_and_exit(wi, errbuf, FARGS);
    }
    ptype = get_atomtype_ptype(type, atype);

    /* Set default from type */
    q0    = get_atomtype_qA(type, atype);
    m0    = get_atomtype_massA(type, atype);
    typeB = type;
    qB    = q0;
    mB    = m0;

    /* Optional parameters */
    nscan = sscanf(line, "%*s%*s%*s%*s%*s%*s%lf%lf%s%lf%lf%s",
                   &q, &m, ctypeB, &qb, &mb, check);

    /* Nasty switch that falls thru all the way down! */
    if (nscan > 0)
    {
        q0 = qB = q;
        if (nscan > 1)
        {
            m0 = mB = m;
            if (nscan > 2)
            {
                if ((typeB = get_atomtype_type(ctypeB, atype)) == NOTSET)
                {
                    sprintf(errbuf, "Atomtype %s not found", ctypeB);
                    warning_error_and_exit(wi, errbuf, FARGS);
                }
                qB = get_atomtype_qA(typeB, atype);
                mB = get_atomtype_massA(typeB, atype);
                if (nscan > 3)
                {
                    qB = qb;
                    if (nscan > 4)
                    {
                        mB = mb;
                        if (nscan > 5)
                        {
                            warning_error(wi, "Too many parameters");
                        }
                    }
                }
            }
        }
    }
    if (debug)
    {
        fprintf(debug, "mB=%g, qB=%g, typeB=%d\n", mB, qB, typeB);
    }

    push_cg(cgs, lastcg, cgnumber, nr);

    push_atom_now(symtab, at, atomnr, get_atomtype_atomnumber(type, atype),
                  type, ctype, ptype, resnumberic,
                  resname, name, m0, q0, typeB,
                  typeB == type ? ctype : ctypeB, mB, qB, wi);
}

void push_molt(t_symtab *symtab, int *nmol, t_molinfo **mol, char *line,
               warninp_t wi)
{
    char       type[STRLEN];
    int        nrexcl, i;
    t_molinfo *newmol;
    char       errbuf[STRLEN];

    if ((sscanf(line, "%s%d", type, &nrexcl)) != 2)
    {
        warning_error(wi, "Expected a molecule type name and nrexcl");
    }

    /* Test if this moleculetype overwrites another */
    i    = 0;
    while (i < *nmol)
    {
        if (strcmp(*((*mol)[i].name), type) == 0)
        {
            sprintf(errbuf, "moleculetype %s is redefined", type);
            warning_error_and_exit(wi, errbuf, FARGS);
        }
        i++;
    }

    (*nmol)++;
    srenew(*mol, *nmol);
    newmol = &((*mol)[*nmol-1]);
    init_molinfo(newmol);

    /* Fill in the values */
    newmol->name     = put_symtab(symtab, type);
    newmol->nrexcl   = nrexcl;
    newmol->excl_set = FALSE;
}

static gmx_bool default_nb_params(int ftype, t_params bt[], t_atoms *at,
                                  t_param *p, int c_start, gmx_bool bB, gmx_bool bGenPairs)
{
    int          i, j, ti, tj, ntype;
    gmx_bool     bFound;
    t_param     *pi    = nullptr;
    int          nr    = bt[ftype].nr;
    int          nral  = NRAL(ftype);
    int          nrfp  = interaction_function[ftype].nrfpA;
    int          nrfpB = interaction_function[ftype].nrfpB;

    if ((!bB && nrfp == 0) || (bB && nrfpB == 0))
    {
        return TRUE;
    }

    bFound = FALSE;
    if (bGenPairs)
    {
        /* First test the generated-pair position to save
         * time when we have 1000*1000 entries for e.g. OPLS...
         */
        ntype = static_cast<int>(std::sqrt(static_cast<double>(nr)));
        GMX_ASSERT(ntype * ntype == nr, "Number of pairs of generated non-bonded parameters should be a perfect square");
        if (bB)
        {
            ti = at->atom[p->a[0]].typeB;
            tj = at->atom[p->a[1]].typeB;
        }
        else
        {
            ti = at->atom[p->a[0]].type;
            tj = at->atom[p->a[1]].type;
        }
        pi     = &(bt[ftype].param[ntype*ti+tj]);
        bFound = ((ti == pi->a[0]) && (tj == pi->a[1]));
    }

    /* Search explicitly if we didnt find it */
    if (!bFound)
    {
        for (i = 0; ((i < nr) && !bFound); i++)
        {
            pi = &(bt[ftype].param[i]);
            if (bB)
            {
                for (j = 0; ((j < nral) &&
                             (at->atom[p->a[j]].typeB == pi->a[j])); j++)
                {
                    ;
                }
            }
            else
            {
                for (j = 0; ((j < nral) &&
                             (at->atom[p->a[j]].type == pi->a[j])); j++)
                {
                    ;
                }
            }
            bFound = (j == nral);
        }
    }

    if (bFound)
    {
        if (bB)
        {
            if (nrfp+nrfpB > MAXFORCEPARAM)
            {
                gmx_incons("Too many force parameters");
            }
            for (j = c_start; (j < nrfpB); j++)
            {
                p->c[nrfp+j] = pi->c[j];
            }
        }
        else
        {
            for (j = c_start; (j < nrfp); j++)
            {
                p->c[j] = pi->c[j];
            }
        }
    }
    else
    {
        for (j = c_start; (j < nrfp); j++)
        {
            p->c[j] = 0.0;
        }
    }
    return bFound;
}

static gmx_bool default_cmap_params(t_params bondtype[],
                                    t_atoms *at, gpp_atomtype_t atype,
                                    t_param *p, gmx_bool bB,
                                    int *cmap_type, int *nparam_def,
                                    warninp_t wi)
{
    int        i, nparam_found;
    int        ct;
    gmx_bool   bFound = FALSE;
    char       errbuf[STRLEN];

    nparam_found = 0;
    ct           = 0;

    /* Match the current cmap angle against the list of cmap_types */
    for (i = 0; i < bondtype[F_CMAP].nct && !bFound; i += 6)
    {
        if (bB)
        {

        }
        else
        {
            if (
                (get_atomtype_batype(at->atom[p->a[0]].type, atype) == bondtype[F_CMAP].cmap_types[i])   &&
                (get_atomtype_batype(at->atom[p->a[1]].type, atype) == bondtype[F_CMAP].cmap_types[i+1]) &&
                (get_atomtype_batype(at->atom[p->a[2]].type, atype) == bondtype[F_CMAP].cmap_types[i+2]) &&
                (get_atomtype_batype(at->atom[p->a[3]].type, atype) == bondtype[F_CMAP].cmap_types[i+3]) &&
                (get_atomtype_batype(at->atom[p->a[4]].type, atype) == bondtype[F_CMAP].cmap_types[i+4]))
            {
                /* Found cmap torsion */
                bFound       = TRUE;
                ct           = bondtype[F_CMAP].cmap_types[i+5];
                nparam_found = 1;
            }
        }
    }

    /* If we did not find a matching type for this cmap torsion */
    if (!bFound)
    {
        sprintf(errbuf, "Unknown cmap torsion between atoms %d %d %d %d %d",
                p->a[0]+1, p->a[1]+1, p->a[2]+1, p->a[3]+1, p->a[4]+1);
        warning_error_and_exit(wi, errbuf, FARGS);
    }

    *nparam_def = nparam_found;
    *cmap_type  = ct;

    return bFound;
}

/* Returns the number of exact atom type matches, i.e. non wild-card matches,
 * returns -1 when there are no matches at all.
 */
static int natom_match(t_param *pi,
                       int type_i, int type_j, int type_k, int type_l,
                       const gpp_atomtype_t atype)
{
    if ((pi->ai() == -1 || get_atomtype_batype(type_i, atype) == pi->ai()) &&
        (pi->aj() == -1 || get_atomtype_batype(type_j, atype) == pi->aj()) &&
        (pi->ak() == -1 || get_atomtype_batype(type_k, atype) == pi->ak()) &&
        (pi->al() == -1 || get_atomtype_batype(type_l, atype) == pi->al()))
    {
        return
            (pi->ai() == -1 ? 0 : 1) +
            (pi->aj() == -1 ? 0 : 1) +
            (pi->ak() == -1 ? 0 : 1) +
            (pi->al() == -1 ? 0 : 1);
    }
    else
    {
        return -1;
    }
}

static gmx_bool default_params(int ftype, t_params bt[],
                               t_atoms *at, gpp_atomtype_t atype,
                               t_param *p, gmx_bool bB,
                               t_param **param_def,
                               int *nparam_def)
{
    int          nparam_found;
    gmx_bool     bFound, bSame;
    t_param     *pi    = nullptr;
    t_param     *pj    = nullptr;
    int          nr    = bt[ftype].nr;
    int          nral  = NRAL(ftype);
    int          nrfpA = interaction_function[ftype].nrfpA;
    int          nrfpB = interaction_function[ftype].nrfpB;

    if ((!bB && nrfpA == 0) || (bB && nrfpB == 0))
    {
        return TRUE;
    }


    bFound       = FALSE;
    nparam_found = 0;
    if (ftype == F_PDIHS || ftype == F_RBDIHS || ftype == F_IDIHS || ftype == F_PIDIHS)
    {
        int nmatch_max = -1;
        int i          = -1;
        int t;

        /* For dihedrals we allow wildcards. We choose the first type
         * that has the most real matches, i.e. non-wildcard matches.
         */
        for (t = 0; ((t < nr) && nmatch_max < 4); t++)
        {
            int      nmatch;
            t_param *pt;

            pt = &(bt[ftype].param[t]);
            if (bB)
            {
                nmatch = natom_match(pt, at->atom[p->ai()].typeB, at->atom[p->aj()].typeB, at->atom[p->ak()].typeB, at->atom[p->al()].typeB, atype);
            }
            else
            {
                nmatch = natom_match(pt, at->atom[p->ai()].type, at->atom[p->aj()].type, at->atom[p->ak()].type, at->atom[p->al()].type, atype);
            }
            if (nmatch > nmatch_max)
            {
                nmatch_max = nmatch;
                i          = t;
                bFound     = TRUE;
            }
        }

        if (bFound == TRUE)
        {
            int j;

            pi    = &(bt[ftype].param[i]);
            nparam_found++;

            /* Find additional matches for this dihedral - necessary
             * for ftype==9.
             * The rule in that case is that additional matches
             * HAVE to be on adjacent lines!
             */
            bSame = TRUE;
            /* Continue from current i value */
            for (j = i + 2; j < nr && bSame; j += 2)
            {
                pj    = &(bt[ftype].param[j]);
                bSame = (pi->ai() == pj->ai() && pi->aj() == pj->aj() && pi->ak() == pj->ak() && pi->al() == pj->al());
                if (bSame)
                {
                    nparam_found++;
                }
                /* nparam_found will be increased as long as the numbers match */
            }
        }
    }
    else   /* Not a dihedral */
    {
        int i, j;

        for (i = 0; ((i < nr) && !bFound); i++)
        {
            pi = &(bt[ftype].param[i]);
            if (bB)
            {
                for (j = 0; ((j < nral) &&
                             (get_atomtype_batype(at->atom[p->a[j]].typeB, atype) == pi->a[j])); j++)
                {
                    ;
                }
            }
            else
            {
                for (j = 0; ((j < nral) &&
                             (get_atomtype_batype(at->atom[p->a[j]].type, atype) == pi->a[j])); j++)
                {
                    ;
                }
            }
            bFound = (j == nral);
        }
        if (bFound)
        {
            nparam_found = 1;
        }
    }

    *param_def  = pi;
    *nparam_def = nparam_found;

    return bFound;
}



void push_bond(directive d, t_params bondtype[], t_params bond[],
               t_atoms *at, gpp_atomtype_t atype, char *line,
               gmx_bool bBonded, gmx_bool bGenPairs, real fudgeQQ,
               gmx_bool bZero, gmx_bool *bWarn_copy_A_B,
               warninp_t wi)
{
    const char  *aaformat[MAXATOMLIST] = {
        "%d%d",
        "%d%d%d",
        "%d%d%d%d",
        "%d%d%d%d%d",
        "%d%d%d%d%d%d",
        "%d%d%d%d%d%d%d"
    };
    const char  *asformat[MAXATOMLIST] = {
        "%*s%*s",
        "%*s%*s%*s",
        "%*s%*s%*s%*s",
        "%*s%*s%*s%*s%*s",
        "%*s%*s%*s%*s%*s%*s",
        "%*s%*s%*s%*s%*s%*s%*s"
    };
    const char  *ccformat = "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf";
    int          nr, i, j, nral, nral_fmt, nread, ftype;
    char         format[STRLEN];
    /* One force parameter more, so we can check if we read too many */
    double       cc[MAXFORCEPARAM+1];
    int          aa[MAXATOMLIST+1];
    t_param      param, *param_defA, *param_defB;
    gmx_bool     bFoundA = FALSE, bFoundB = FALSE, bDef, bPert, bSwapParity = FALSE;
    int          nparam_defA, nparam_defB;
    char         errbuf[STRLEN];

    nparam_defA = nparam_defB = 0;

    ftype = ifunc_index(d, 1);
    nral  = NRAL(ftype);
    for (j = 0; j < MAXATOMLIST; j++)
    {
        aa[j] = NOTSET;
    }
    bDef = (NRFP(ftype) > 0);

    if (ftype == F_SETTLE)
    {
        /* SETTLE acts on 3 atoms, but the topology format only specifies
         * the first atom (for historical reasons).
         */
        nral_fmt = 1;
    }
    else
    {
        nral_fmt = nral;
    }

    nread = sscanf(line, aaformat[nral_fmt-1],
                   &aa[0], &aa[1], &aa[2], &aa[3], &aa[4], &aa[5]);

    if (ftype == F_SETTLE)
    {
        aa[3] = aa[1];
        aa[1] = aa[0] + 1;
        aa[2] = aa[0] + 2;
    }

    if (nread < nral_fmt)
    {
        too_few(wi);
        return;
    }
    else if (nread > nral_fmt)
    {
        /* this is a hack to allow for virtual sites with swapped parity */
        bSwapParity = (aa[nral] < 0);
        if (bSwapParity)
        {
            aa[nral] = -aa[nral];
        }
        ftype = ifunc_index(d, aa[nral]);
        if (bSwapParity)
        {
            switch (ftype)
            {
                case F_VSITE3FAD:
                case F_VSITE3OUT:
                    break;
                default:
                    sprintf(errbuf, "Negative function types only allowed for %s and %s",
                            interaction_function[F_VSITE3FAD].longname,
                            interaction_function[F_VSITE3OUT].longname);
                    warning_error_and_exit(wi, errbuf, FARGS);
            }
        }
    }


    /* Check for double atoms and atoms out of bounds */
    for (i = 0; (i < nral); i++)
    {
        if (aa[i] < 1 || aa[i] > at->nr)
        {
            sprintf(errbuf, "Atom index (%d) in %s out of bounds (1-%d).\n"
                    "This probably means that you have inserted topology section \"%s\"\n"
                    "in a part belonging to a different molecule than you intended to.\n"
                    "In that case move the \"%s\" section to the right molecule.",
                    aa[i], dir2str(d), at->nr, dir2str(d), dir2str(d));
            warning_error_and_exit(wi, errbuf, FARGS);
        }
        for (j = i+1; (j < nral); j++)
        {
            GMX_ASSERT(j < MAXATOMLIST + 1, "Values from nral=NRAL() will satisfy this, we assert to keep gcc 4 happy");
            if (aa[i] == aa[j])
            {
                sprintf(errbuf, "Duplicate atom index (%d) in %s", aa[i], dir2str(d));
                if (ftype == F_ANGRES)
                {
                    /* Since the angle restraints uses 2 pairs of atoms to
                     * defines an angle between vectors, it can be useful
                     * to use one atom twice, so we only issue a note here.
                     */
                    warning_note(wi, errbuf);
                }
                else
                {
                    warning_error(wi, errbuf);
                }
            }
        }
    }

    /* default force parameters  */
    for (j = 0; (j < MAXATOMLIST); j++)
    {
        param.a[j] = aa[j]-1;
    }
    for (j = 0; (j < MAXFORCEPARAM); j++)
    {
        param.c[j] = 0.0;
    }

    /* Get force params for normal and free energy perturbation
     * studies, as determined by types!
     */

    if (bBonded)
    {
        bFoundA = default_params(ftype, bondtype, at, atype, &param, FALSE, &param_defA, &nparam_defA);
        if (bFoundA)
        {
            /* Copy the A-state and B-state default parameters. */
            GMX_ASSERT(NRFPA(ftype)+NRFPB(ftype) <= MAXFORCEPARAM, "Bonded interactions may have at most 12 parameters");
            for (j = 0; (j < NRFPA(ftype)+NRFPB(ftype)); j++)
            {
                param.c[j] = param_defA->c[j];
            }
        }
        bFoundB = default_params(ftype, bondtype, at, atype, &param, TRUE, &param_defB, &nparam_defB);
        if (bFoundB)
        {
            /* Copy only the B-state default parameters */
            for (j = NRFPA(ftype); (j < NRFP(ftype)); j++)
            {
                param.c[j] = param_defB->c[j];
            }
        }
    }
    else if (ftype == F_LJ14)
    {
        bFoundA = default_nb_params(ftype, bondtype, at, &param, 0, FALSE, bGenPairs);
        bFoundB = default_nb_params(ftype, bondtype, at, &param, 0, TRUE, bGenPairs);
    }
    else if (ftype == F_LJC14_Q)
    {
        param.c[0] = fudgeQQ;
        /* Fill in the A-state charges as default parameters */
        param.c[1] = at->atom[param.a[0]].q;
        param.c[2] = at->atom[param.a[1]].q;
        /* The default LJ parameters are the standard 1-4 parameters */
        bFoundA = default_nb_params(F_LJ14, bondtype, at, &param, 3, FALSE, bGenPairs);
        bFoundB = TRUE;
    }
    else if (ftype == F_LJC_PAIRS_NB)
    {
        /* Defaults are not supported here */
        bFoundA = FALSE;
        bFoundB = TRUE;
    }
    else
    {
        gmx_incons("Unknown function type in push_bond");
    }

    if (nread > nral_fmt)
    {
        /* Manually specified parameters - in this case we discard multiple torsion info! */

        strcpy(format, asformat[nral_fmt-1]);
        strcat(format, ccformat);

        nread = sscanf(line, format, &cc[0], &cc[1], &cc[2], &cc[3], &cc[4], &cc[5],
                       &cc[6], &cc[7], &cc[8], &cc[9], &cc[10], &cc[11], &cc[12]);

        if ((nread == NRFPA(ftype)) && (NRFPB(ftype) != 0))
        {
            /* We only have to issue a warning if these atoms are perturbed! */
            bPert = FALSE;
            for (j = 0; (j < nral); j++)
            {
                bPert = bPert || PERTURBED(at->atom[param.a[j]]);
            }

            if (bPert && *bWarn_copy_A_B)
            {
                sprintf(errbuf,
                        "Some parameters for bonded interaction involving perturbed atoms are specified explicitly in state A, but not B - copying A to B");
                warning(wi, errbuf);
                *bWarn_copy_A_B = FALSE;
            }

            /* If only the A parameters were specified, copy them to the B state */
            /* The B-state parameters correspond to the first nrfpB
             * A-state parameters.
             */
            for (j = 0; (j < NRFPB(ftype)); j++)
            {
                cc[nread++] = cc[j];
            }
        }

        /* If nread was 0 or EOF, no parameters were read => use defaults.
         * If nread was nrfpA we copied above so nread=nrfp.
         * If nread was nrfp we are cool.
         * For F_LJC14_Q we allow supplying fudgeQQ only.
         * Anything else is an error!
         */
        if ((nread != 0) && (nread != EOF) && (nread != NRFP(ftype)) &&
            !(ftype == F_LJC14_Q && nread == 1))
        {
            sprintf(errbuf, "Incorrect number of parameters - found %d, expected %d or %d for %s (after the function type).",
                    nread, NRFPA(ftype), NRFP(ftype), interaction_function[ftype].longname);
            warning_error_and_exit(wi, errbuf, FARGS);
        }

        for (j = 0; (j < nread); j++)
        {
            param.c[j] = cc[j];
        }

        /* Check whether we have to use the defaults */
        if (nread == NRFP(ftype))
        {
            bDef = FALSE;
        }
    }
    else
    {
        nread = 0;
    }
    /* nread now holds the number of force parameters read! */

    if (bDef)
    {
        /* Use defaults */
        /* When we have multiple terms it would be very dangerous to allow perturbations to a different atom type! */
        if (ftype == F_PDIHS)
        {
            if ((nparam_defA != nparam_defB) || ((nparam_defA > 1 || nparam_defB > 1) && (param_defA != param_defB)))
            {
                sprintf(errbuf,
                        "Cannot automatically perturb a torsion with multiple terms to different form.\n"
                        "Please specify perturbed parameters manually for this torsion in your topology!");
                warning_error(wi, errbuf);
            }
        }

        if (nread > 0 && nread < NRFPA(ftype))
        {
            /* Issue an error, do not use defaults */
            sprintf(errbuf, "Not enough parameters, there should be at least %d (or 0 for defaults)", NRFPA(ftype));
            warning_error(wi, errbuf);
        }

        if (nread == 0 || nread == EOF)
        {
            if (!bFoundA)
            {
                if (interaction_function[ftype].flags & IF_VSITE)
                {
                    /* set them to NOTSET, will be calculated later */
                    for (j = 0; (j < MAXFORCEPARAM); j++)
                    {
                        param.c[j] = NOTSET;
                    }

                    if (bSwapParity)
                    {
                        param.c1() = -1; /* flag to swap parity of vsite construction */
                    }
                }
                else
                {
                    if (bZero)
                    {
                        fprintf(stderr, "NOTE: No default %s types, using zeroes\n",
                                interaction_function[ftype].longname);
                    }
                    else
                    {
                        sprintf(errbuf, "No default %s types", interaction_function[ftype].longname);
                        warning_error(wi, errbuf);
                    }
                }
            }
            else
            {
                if (bSwapParity)
                {
                    switch (ftype)
                    {
                        case F_VSITE3FAD:
                            param.c0() = 360 - param.c0();
                            break;
                        case F_VSITE3OUT:
                            param.c2() = -param.c2();
                            break;
                    }
                }
            }
            if (!bFoundB)
            {
                /* We only have to issue a warning if these atoms are perturbed! */
                bPert = FALSE;
                for (j = 0; (j < nral); j++)
                {
                    bPert = bPert || PERTURBED(at->atom[param.a[j]]);
                }

                if (bPert)
                {
                    sprintf(errbuf, "No default %s types for perturbed atoms, "
                            "using normal values", interaction_function[ftype].longname);
                    warning(wi, errbuf);
                }
            }
        }
    }

    if ((ftype == F_PDIHS || ftype == F_ANGRES || ftype == F_ANGRESZ)
        && param.c[5] != param.c[2])
    {
        sprintf(errbuf, "%s multiplicity can not be perturbed %f!=%f",
                interaction_function[ftype].longname,
                param.c[2], param.c[5]);
        warning_error_and_exit(wi, errbuf, FARGS);
    }

    if (IS_TABULATED(ftype) && param.c[2] != param.c[0])
    {
        sprintf(errbuf, "%s table number can not be perturbed %d!=%d",
                interaction_function[ftype].longname,
                (int)(param.c[0]+0.5), (int)(param.c[2]+0.5));
        warning_error_and_exit(wi, errbuf, FARGS);
    }

    /* Dont add R-B dihedrals where all parameters are zero (no interaction) */
    if (ftype == F_RBDIHS)
    {
        nr = 0;
        for (i = 0; i < NRFP(ftype); i++)
        {
            if (param.c[i] != 0)
            {
                nr++;
            }
        }
        if (nr == 0)
        {
            return;
        }
    }

    /* Put the values in the appropriate arrays */
    add_param_to_list (&bond[ftype], &param);

    /* Push additional torsions from FF for ftype==9 if we have them.
     * We have already checked that the A/B states do not differ in this case,
     * so we do not have to double-check that again, or the vsite stuff.
     * In addition, those torsions cannot be automatically perturbed.
     */
    if (bDef && ftype == F_PDIHS)
    {
        for (i = 1; i < nparam_defA; i++)
        {
            /* Advance pointer! */
            param_defA += 2;
            for (j = 0; (j < NRFPA(ftype)+NRFPB(ftype)); j++)
            {
                param.c[j] = param_defA->c[j];
            }
            /* And push the next term for this torsion */
            add_param_to_list (&bond[ftype], &param);
        }
    }
}

void push_cmap(directive d, t_params bondtype[], t_params bond[],
               t_atoms *at, gpp_atomtype_t atype, char *line,
               warninp_t wi)
{
    const char *aaformat[MAXATOMLIST+1] =
    {
        "%d",
        "%d%d",
        "%d%d%d",
        "%d%d%d%d",
        "%d%d%d%d%d",
        "%d%d%d%d%d%d",
        "%d%d%d%d%d%d%d"
    };

    int      i, j, ftype, nral, nread, ncmap_params;
    int      cmap_type;
    int      aa[MAXATOMLIST+1];
    char     errbuf[STRLEN];
    gmx_bool bFound;
    t_param  param;

    ftype        = ifunc_index(d, 1);
    nral         = NRAL(ftype);
    ncmap_params = 0;

    nread = sscanf(line, aaformat[nral-1],
                   &aa[0], &aa[1], &aa[2], &aa[3], &aa[4], &aa[5]);

    if (nread < nral)
    {
        too_few(wi);
        return;
    }
    else if (nread == nral)
    {
        ftype = ifunc_index(d, 1);
    }

    /* Check for double atoms and atoms out of bounds */
    for (i = 0; i < nral; i++)
    {
        if (aa[i] < 1 || aa[i] > at->nr)
        {
            sprintf(errbuf, "Atom index (%d) in %s out of bounds (1-%d).\n"
                    "This probably means that you have inserted topology section \"%s\"\n"
                    "in a part belonging to a different molecule than you intended to.\n"
                    "In that case move the \"%s\" section to the right molecule.",
                    aa[i], dir2str(d), at->nr, dir2str(d), dir2str(d));
            warning_error_and_exit(wi, errbuf, FARGS);
        }

        for (j = i+1; (j < nral); j++)
        {
            if (aa[i] == aa[j])
            {
                sprintf(errbuf, "Duplicate atom index (%d) in %s", aa[i], dir2str(d));
                warning_error(wi, errbuf);
            }
        }
    }

    /* default force parameters  */
    for (j = 0; (j < MAXATOMLIST); j++)
    {
        param.a[j] = aa[j]-1;
    }
    for (j = 0; (j < MAXFORCEPARAM); j++)
    {
        param.c[j] = 0.0;
    }

    /* Get the cmap type for this cmap angle */
    bFound = default_cmap_params(bondtype, at, atype, &param, FALSE, &cmap_type, &ncmap_params, wi);

    /* We want exactly one parameter (the cmap type in state A (currently no state B) back */
    if (bFound && ncmap_params == 1)
    {
        /* Put the values in the appropriate arrays */
        param.c[0] = cmap_type;
        add_param_to_list(&bond[ftype], &param);
    }
    else
    {
        /* This is essentially the same check as in default_cmap_params() done one more time */
        sprintf(errbuf, "Unable to assign a cmap type to torsion %d %d %d %d and %d\n",
                param.a[0]+1, param.a[1]+1, param.a[2]+1, param.a[3]+1, param.a[4]+1);
        warning_error_and_exit(wi, errbuf, FARGS);
    }
}



void push_vsitesn(directive d, t_params bond[],
                  t_atoms *at, char *line,
                  warninp_t wi)
{
    char   *ptr;
    int     type, ftype, j, n, ret, nj, a;
    int    *atc    = nullptr;
    double *weight = nullptr, weight_tot;
    t_param param;
    char    errbuf[STRLEN];

    /* default force parameters  */
    for (j = 0; (j < MAXATOMLIST); j++)
    {
        param.a[j] = NOTSET;
    }
    for (j = 0; (j < MAXFORCEPARAM); j++)
    {
        param.c[j] = 0.0;
    }

    ptr  = line;
    ret  = sscanf(ptr, "%d%n", &a, &n);
    ptr += n;
    if (ret == 0)
    {
        sprintf(errbuf, "Expected an atom index in section \"%s\"",
                dir2str(d));
        warning_error_and_exit(wi, errbuf, FARGS);
    }

    param.a[0] = a - 1;

    sscanf(ptr, "%d%n", &type, &n);
    ptr  += n;
    ftype = ifunc_index(d, type);

    weight_tot = 0;
    nj         = 0;
    do
    {
        ret  = sscanf(ptr, "%d%n", &a, &n);
        ptr += n;
        if (ret > 0)
        {
            if (nj % 20 == 0)
            {
                srenew(atc, nj+20);
                srenew(weight, nj+20);
            }
            atc[nj] = a - 1;
            switch (type)
            {
                case 1:
                    weight[nj] = 1;
                    break;
                case 2:
                    /* Here we use the A-state mass as a parameter.
                     * Note that the B-state mass has no influence.
                     */
                    weight[nj] = at->atom[atc[nj]].m;
                    break;
                case 3:
                    weight[nj] = -1;
                    ret        = sscanf(ptr, "%lf%n", &(weight[nj]), &n);
                    ptr       += n;
                    if (weight[nj] < 0)
                    {
                        sprintf(errbuf, "No weight or negative weight found for vsiten constructing atom %d (atom index %d)",
                                nj+1, atc[nj]+1);
                        warning_error_and_exit(wi, errbuf, FARGS);
                    }
                    break;
                default:
                    sprintf(errbuf, "Unknown vsiten type %d", type);
                    warning_error_and_exit(wi, errbuf, FARGS);
            }
            weight_tot += weight[nj];
            nj++;
        }
    }
    while (ret > 0);

    if (nj == 0)
    {
        sprintf(errbuf, "Expected more than one atom index in section \"%s\"",
                dir2str(d));
        warning_error_and_exit(wi, errbuf, FARGS);
    }

    if (weight_tot == 0)
    {
        warning_error_and_exit(wi, "The total mass of the construting atoms is zero", FARGS);
    }

    for (j = 0; j < nj; j++)
    {
        param.a[1] = atc[j];
        param.c[0] = nj;
        param.c[1] = weight[j]/weight_tot;
        /* Put the values in the appropriate arrays */
        add_param_to_list (&bond[ftype], &param);
    }

    sfree(atc);
    sfree(weight);
}

void push_mol(int nrmols, t_molinfo mols[], char *pline, int *whichmol,
              int *nrcopies,
              warninp_t wi)
{
    char type[STRLEN];
    char errbuf[STRLEN];

    if (sscanf(pline, "%s%d", type, nrcopies) != 2)
    {
        too_few(wi);
        return;
    }

    /* Search moleculename.
     * Here we originally only did case insensitive matching. But because
     * some PDB files can have many chains and use case to generate more
     * chain-identifiers, which in turn end up in our moleculetype name,
     * we added support for case-sensitivity.
     */
    int nrcs    = 0;
    int nrci    = 0;
    int matchci = -1;
    int matchcs = -1;
    for (int i = 0; i < nrmols; i++)
    {
        if (strcmp(type, *(mols[i].name)) == 0)
        {
            nrcs++;
            matchcs = i;
        }
        if (gmx_strcasecmp(type, *(mols[i].name)) == 0)
        {
            nrci++;
            matchci = i;
        }
    }

    if (nrcs == 1)
    {
        // select the case sensitive match
        *whichmol = matchcs;
    }
    else
    {
        // avoid matching case-insensitive when we have multiple matches
        if (nrci > 1)
        {
            sprintf(errbuf, "For moleculetype '%s' in [ system ] %d case insensitive matches, but %d case sensitive matches were found. Check the case of the characters in the moleculetypes.", type, nrci, nrcs);
            warning_error_and_exit(wi, errbuf, FARGS);
        }
        if (nrci == 1)
        {
            // select the unique case insensitive match
            *whichmol = matchci;
        }
        else
        {
            sprintf(errbuf, "No such moleculetype %s", type);
            warning_error_and_exit(wi, errbuf, FARGS);
        }
    }
}

void init_block2(t_block2 *b2, int natom)
{
    int i;

    b2->nr = natom;
    snew(b2->nra, b2->nr);
    snew(b2->a, b2->nr);
    for (i = 0; (i < b2->nr); i++)
    {
        b2->a[i] = nullptr;
    }
}

void done_block2(t_block2 *b2)
{
    int i;

    if (b2->nr)
    {
        for (i = 0; (i < b2->nr); i++)
        {
            sfree(b2->a[i]);
        }
        sfree(b2->a);
        sfree(b2->nra);
        b2->nr = 0;
    }
}

void push_excl(char *line, t_block2 *b2, warninp_t wi)
{
    int  i, j;
    int  n;
    char base[STRLEN], format[STRLEN];
    char errbuf[STRLEN];

    if (sscanf(line, "%d", &i) == 0)
    {
        return;
    }

    if ((1 <= i) && (i <= b2->nr))
    {
        i--;
    }
    else
    {
        if (debug)
        {
            fprintf(debug, "Unbound atom %d\n", i-1);
        }
        return;
    }
    strcpy(base, "%*d");
    do
    {
        strcpy(format, base);
        strcat(format, "%d");
        n = sscanf(line, format, &j);
        if (n == 1)
        {
            if ((1 <= j) && (j <= b2->nr))
            {
                j--;
                srenew(b2->a[i], ++(b2->nra[i]));
                b2->a[i][b2->nra[i]-1] = j;
                /* also add the reverse exclusion! */
                srenew(b2->a[j], ++(b2->nra[j]));
                b2->a[j][b2->nra[j]-1] = i;
                strcat(base, "%*d");
            }
            else
            {
                sprintf(errbuf, "Invalid Atomnr j: %d, b2->nr: %d\n", j, b2->nr);
                warning_error_and_exit(wi, errbuf, FARGS);
            }
        }
    }
    while (n == 1);
}

void b_to_b2(t_blocka *b, t_block2 *b2)
{
    int     i;
    int     j, a;

    for (i = 0; (i < b->nr); i++)
    {
        for (j = b->index[i]; (j < b->index[i+1]); j++)
        {
            a = b->a[j];
            srenew(b2->a[i], ++b2->nra[i]);
            b2->a[i][b2->nra[i]-1] = a;
        }
    }
}

void b2_to_b(t_block2 *b2, t_blocka *b)
{
    int     i, nra;
    int     j;

    nra = 0;
    for (i = 0; (i < b2->nr); i++)
    {
        b->index[i] = nra;
        for (j = 0; (j < b2->nra[i]); j++)
        {
            b->a[nra+j] = b2->a[i][j];
        }
        nra += b2->nra[i];
    }
    /* terminate list */
    b->index[i] = nra;
}

static int icomp(const void *v1, const void *v2)
{
    return (*((int *) v1))-(*((int *) v2));
}

void merge_excl(t_blocka *excl, t_block2 *b2, warninp_t wi)
{
    int     i, k;
    int     j;
    int     nra;
    char    errbuf[STRLEN];

    if (!b2->nr)
    {
        return;
    }
    else if (b2->nr != excl->nr)
    {
        sprintf(errbuf, "DEATH HORROR: b2->nr = %d, while excl->nr = %d",
                b2->nr, excl->nr);
        warning_error_and_exit(wi, errbuf, FARGS);
    }
    else if (debug)
    {
        fprintf(debug, "Entering merge_excl\n");
    }

    /* First copy all entries from excl to b2 */
    b_to_b2(excl, b2);

    /* Count and sort the exclusions */
    nra = 0;
    for (i = 0; (i < b2->nr); i++)
    {
        if (b2->nra[i] > 0)
        {
            /* remove double entries */
            qsort(b2->a[i], (size_t)b2->nra[i], (size_t)sizeof(b2->a[i][0]), icomp);
            k = 1;
            for (j = 1; (j < b2->nra[i]); j++)
            {
                if (b2->a[i][j] != b2->a[i][k-1])
                {
                    b2->a[i][k] = b2->a[i][j];
                    k++;
                }
            }
            b2->nra[i] = k;
            nra       += b2->nra[i];
        }
    }
    excl->nra = nra;
    srenew(excl->a, excl->nra);

    b2_to_b(b2, excl);
}

int add_atomtype_decoupled(t_symtab *symtab, gpp_atomtype_t at,
                           t_nbparam ***nbparam, t_nbparam ***pair)
{
    t_atom  atom;
    t_param param;
    int     i, nr;

    /* Define an atom type with all parameters set to zero (no interactions) */
    atom.q = 0.0;
    atom.m = 0.0;
    /* Type for decoupled atoms could be anything,
     * this should be changed automatically later when required.
     */
    atom.ptype = eptAtom;
    for (i = 0; (i < MAXFORCEPARAM); i++)
    {
        param.c[i] = 0.0;
    }

    nr = add_atomtype(at, symtab, &atom, "decoupled", &param, -1, 0.0, 0.0, 0.0, 0, 0, 0);

    /* Add space in the non-bonded parameters matrix */
    realloc_nb_params(at, nbparam, pair);

    return nr;
}

static void convert_pairs_to_pairsQ(t_params *plist,
                                    real fudgeQQ, t_atoms *atoms)
{
    t_param *paramp1, *paramp2, *paramnew;
    int      i, j, p1nr, p2nr, p2newnr;

    /* Add the pair list to the pairQ list */
    p1nr    = plist[F_LJ14].nr;
    p2nr    = plist[F_LJC14_Q].nr;
    p2newnr = p1nr + p2nr;
    snew(paramnew, p2newnr);

    paramp1             = plist[F_LJ14].param;
    paramp2             = plist[F_LJC14_Q].param;

    /* Fill in the new F_LJC14_Q array with the old one. NOTE:
       it may be possible to just ADD the converted F_LJ14 array
       to the old F_LJC14_Q array, but since we have to create
       a new sized memory structure, better just to deep copy it all.
     */

    for (i = 0; i < p2nr; i++)
    {
        /* Copy over parameters */
        for (j = 0; j < 5; j++) /* entries are 0-4 for F_LJC14_Q */
        {
            paramnew[i].c[j] = paramp2[i].c[j];
        }

        /* copy over atoms */
        for (j = 0; j < 2; j++)
        {
            paramnew[i].a[j] = paramp2[i].a[j];
        }
    }

    for (i = p2nr; i < p2newnr; i++)
    {
        j             = i-p2nr;

        /* Copy over parameters */
        paramnew[i].c[0] = fudgeQQ;
        paramnew[i].c[1] = atoms->atom[paramp1[j].a[0]].q;
        paramnew[i].c[2] = atoms->atom[paramp1[j].a[1]].q;
        paramnew[i].c[3] = paramp1[j].c[0];
        paramnew[i].c[4] = paramp1[j].c[1];

        /* copy over atoms */
        paramnew[i].a[0] = paramp1[j].a[0];
        paramnew[i].a[1] = paramp1[j].a[1];
    }

    /* free the old pairlists */
    sfree(plist[F_LJC14_Q].param);
    sfree(plist[F_LJ14].param);

    /* now assign the new data to the F_LJC14_Q structure */
    plist[F_LJC14_Q].param   = paramnew;
    plist[F_LJC14_Q].nr      = p2newnr;

    /* Empty the LJ14 pairlist */
    plist[F_LJ14].nr    = 0;
    plist[F_LJ14].param = nullptr;
}

static void generate_LJCpairsNB(t_molinfo *mol, int nb_funct, t_params *nbp, warninp_t wi)
{
    int       n, ntype, i, j, k;
    t_atom   *atom;
    t_blocka *excl;
    gmx_bool  bExcl;
    t_param   param;
    char      errbuf[STRLEN];

    n    = mol->atoms.nr;
    atom = mol->atoms.atom;

    ntype = static_cast<int>(std::sqrt(static_cast<double>(nbp->nr)));
    GMX_ASSERT(ntype * ntype == nbp->nr, "Number of pairs of generated non-bonded parameters should be a perfect square");

    for (i = 0; i < MAXATOMLIST; i++)
    {
        param.a[i] = NOTSET;
    }
    for (i = 0; i < MAXFORCEPARAM; i++)
    {
        param.c[i] = NOTSET;
    }

    /* Add a pair interaction for all non-excluded atom pairs */
    excl = &mol->excls;
    for (i = 0; i < n; i++)
    {
        for (j = i+1; j < n; j++)
        {
            bExcl = FALSE;
            for (k = excl->index[i]; k < excl->index[i+1]; k++)
            {
                if (excl->a[k] == j)
                {
                    bExcl = TRUE;
                }
            }
            if (!bExcl)
            {
                if (nb_funct != F_LJ)
                {
                    sprintf(errbuf, "Can only generate non-bonded pair interactions for Van der Waals type Lennard-Jones");
                    warning_error_and_exit(wi, errbuf, FARGS);
                }
                param.a[0] = i;
                param.a[1] = j;
                param.c[0] = atom[i].q;
                param.c[1] = atom[j].q;
                param.c[2] = nbp->param[ntype*atom[i].type+atom[j].type].c[0];
                param.c[3] = nbp->param[ntype*atom[i].type+atom[j].type].c[1];
                add_param_to_list(&mol->plist[F_LJC_PAIRS_NB], &param);
            }
        }
    }
}

static void set_excl_all(t_blocka *excl)
{
    int nat, i, j, k;

    /* Get rid of the current exclusions and exclude all atom pairs */
    nat       = excl->nr;
    excl->nra = nat*nat;
    srenew(excl->a, excl->nra);
    k = 0;
    for (i = 0; i < nat; i++)
    {
        excl->index[i] = k;
        for (j = 0; j < nat; j++)
        {
            excl->a[k++] = j;
        }
    }
    excl->index[nat] = k;
}

static void decouple_atoms(t_atoms *atoms, int atomtype_decouple,
                           int couple_lam0, int couple_lam1,
                           const char *mol_name, warninp_t wi)
{
    int  i;
    char errbuf[STRLEN];

    for (i = 0; i < atoms->nr; i++)
    {
        t_atom *atom;

        atom = &atoms->atom[i];

        if (atom->qB != atom->q || atom->typeB != atom->type)
        {
            sprintf(errbuf, "Atom %d in molecule type '%s' has different A and B state charges and/or atom types set in the topology file as well as through the mdp option '%s'. You can not use both these methods simultaneously.",
                    i + 1, mol_name, "couple-moltype");
            warning_error_and_exit(wi, errbuf, FARGS);
        }

        if (couple_lam0 == ecouplamNONE || couple_lam0 == ecouplamVDW)
        {
            atom->q     = 0.0;
        }
        if (couple_lam0 == ecouplamNONE || couple_lam0 == ecouplamQ)
        {
            atom->type  = atomtype_decouple;
        }
        if (couple_lam1 == ecouplamNONE || couple_lam1 == ecouplamVDW)
        {
            atom->qB    = 0.0;
        }
        if (couple_lam1 == ecouplamNONE || couple_lam1 == ecouplamQ)
        {
            atom->typeB = atomtype_decouple;
        }
    }
}

void convert_moltype_couple(t_molinfo *mol, int atomtype_decouple, real fudgeQQ,
                            int couple_lam0, int couple_lam1,
                            gmx_bool bCoupleIntra, int nb_funct, t_params *nbp,
                            warninp_t wi)
{
    convert_pairs_to_pairsQ(mol->plist, fudgeQQ, &mol->atoms);
    if (!bCoupleIntra)
    {
        generate_LJCpairsNB(mol, nb_funct, nbp, wi);
        set_excl_all(&mol->excls);
    }
    decouple_atoms(&mol->atoms, atomtype_decouple, couple_lam0, couple_lam1,
                   *mol->name, wi);
}
