/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014,2015,2017,2018 by the GROMACS development team.
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
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

#include "toputil.h"

#include <climits>
#include <cmath>
#include <cstring>

#include <algorithm>

#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

/* UTILITIES */

void add_param_to_list(InteractionsOfType* list, const InteractionOfType& b)
{
    list->interactionTypes.emplace_back(b);
}

/* PRINTING STRUCTURES */

static void print_bt(FILE*                                   out,
                     Directive                               d,
                     PreprocessingAtomTypes*                 at,
                     int                                     ftype,
                     int                                     fsubtype,
                     gmx::ArrayRef<const InteractionsOfType> plist,
                     bool                                    bFullDih)
{
    /* This dihp is a DIRTY patch because the dih-types do not use
     * all four atoms to determine the type.
     */
    const int dihp[2][2] = { { 1, 2 }, { 0, 3 } };
    int       nral, nrfp;
    bool      bDih = false, bSwapParity;

    const InteractionsOfType* bt = &(plist[ftype]);

    if (bt->size() == 0)
    {
        return;
    }

    int f = 0;
    switch (ftype)
    {
        case F_G96ANGLES: // Intended to fall through
        case F_G96BONDS: f = 1; break;
        case F_MORSE: f = 2; break;
        case F_CUBICBONDS: f = 3; break;
        case F_CONNBONDS: f = 4; break;
        case F_HARMONIC: f = 5; break;
        case F_CROSS_BOND_ANGLES: f = 2; break;
        case F_CROSS_BOND_BONDS: f = 3; break;
        case F_UREY_BRADLEY: f = 4; break;
        case F_PDIHS:  // Intended to fall through
        case F_RBDIHS: // Intended to fall through
        case F_FOURDIHS: bDih = TRUE; break;
        case F_IDIHS:
            f    = 1;
            bDih = TRUE;
            break;
        case F_CONSTRNC: // Intended to fall through
        case F_VSITE3FD: f = 1; break;
        case F_VSITE3FAD: f = 2; break;
        case F_VSITE3OUT: f = 3; break;
        case F_VSITE4FDN: // Intended to fall through
        case F_CMAP: f = 1; break;

        default: bDih = FALSE;
    }
    if (bFullDih)
    {
        bDih = FALSE;
    }
    if (fsubtype)
    {
        f = fsubtype - 1;
    }

    nral = NRAL(ftype);
    nrfp = NRFP(ftype);

    /* header */
    fprintf(out, "[ %s ]\n", dir2str(d));
    fprintf(out, "; ");
    if (!bDih)
    {
        fprintf(out, "%3s  %4s", "ai", "aj");
        for (int j = 2; (j < nral); j++)
        {
            fprintf(out, "  %3c%c", 'a', 'i' + j);
        }
    }
    else
    {
        for (int j = 0; (j < 2); j++)
        {
            fprintf(out, "%3c%c", 'a', 'i' + dihp[f][j]);
        }
    }

    fprintf(out, " funct");
    for (int j = 0; (j < nrfp); j++)
    {
        fprintf(out, " %12c%1d", 'c', j);
    }
    fprintf(out, "\n");

    /* print bondtypes */
    for (const auto& parm : bt->interactionTypes)
    {
        bSwapParity                    = (parm.c0() == NOTSET) && (parm.c1() == -1);
        gmx::ArrayRef<const int> atoms = parm.atoms();
        if (!bDih)
        {
            for (int j = 0; (j < nral); j++)
            {
                fprintf(out, "%5s ", at->atomNameFromAtomType(atoms[j]));
            }
        }
        else
        {
            for (int j = 0; (j < 2); j++)
            {
                fprintf(out, "%5s ", at->atomNameFromAtomType(atoms[dihp[f][j]]));
            }
        }
        fprintf(out, "%5d ", bSwapParity ? -f - 1 : f + 1);

        if (!parm.interactionTypeName().empty())
        {
            fprintf(out, "   %s", parm.interactionTypeName().c_str());
        }
        else
        {
            gmx::ArrayRef<const real> forceParam = parm.forceParam();
            for (int j = 0; (j < nrfp) && (forceParam[j] != NOTSET); j++)
            {
                fprintf(out, "%13.6e ", forceParam[j]);
            }
        }

        fprintf(out, "\n");
    }
    fprintf(out, "\n");
    fflush(out);
}

void print_excl(FILE* out, int natoms, t_excls excls[])
{
    int  i;
    bool have_excl;
    int  j;

    have_excl = FALSE;
    for (i = 0; i < natoms && !have_excl; i++)
    {
        have_excl = (excls[i].nr > 0);
    }

    if (have_excl)
    {
        fprintf(out, "[ %s ]\n", dir2str(Directive::d_exclusions));
        fprintf(out, "; %4s    %s\n", "i", "excluded from i");
        for (i = 0; i < natoms; i++)
        {
            if (excls[i].nr > 0)
            {
                fprintf(out, "%6d ", i + 1);
                for (j = 0; j < excls[i].nr; j++)
                {
                    fprintf(out, " %5d", excls[i].e[j] + 1);
                }
                fprintf(out, "\n");
            }
        }
        fprintf(out, "\n");
        fflush(out);
    }
}

static double get_residue_charge(const t_atoms* atoms, int at)
{
    int    ri;
    double q;

    ri = atoms->atom[at].resind;
    q  = 0;
    while (at < atoms->nr && atoms->atom[at].resind == ri)
    {
        q += atoms->atom[at].q;
        at++;
    }

    return q;
}

void print_atoms(FILE* out, PreprocessingAtomTypes* atype, t_atoms* at, int* cgnr, bool bRTPresname)
{
    int         i, ri;
    int         tpA, tpB;
    const char* as;
    const char *tpnmA, *tpnmB;
    double      qres, qtot;

    as = dir2str(Directive::d_atoms);
    fprintf(out, "[ %s ]\n", as);
    fprintf(out, "; %4s %10s %6s %7s%6s %6s %10s %10s %6s %10s %10s\n", "nr", "type", "resnr",
            "residue", "atom", "cgnr", "charge", "mass", "typeB", "chargeB", "massB");

    qtot = 0;

    if (at->nres)
    {
        /* if the information is present... */
        for (i = 0; (i < at->nr); i++)
        {
            ri = at->atom[i].resind;
            if ((i == 0 || ri != at->atom[i - 1].resind) && at->resinfo[ri].rtp != nullptr)
            {
                qres = get_residue_charge(at, i);
                fprintf(out, "; residue %3d %-3s rtp %-4s q ", at->resinfo[ri].nr,
                        *at->resinfo[ri].name, *at->resinfo[ri].rtp);
                if (fabs(qres) < 0.001)
                {
                    fprintf(out, " %s", "0.0");
                }
                else
                {
                    fprintf(out, "%+3.1f", qres);
                }
                fprintf(out, "\n");
            }
            tpA = at->atom[i].type;
            if ((tpnmA = atype->atomNameFromAtomType(tpA)) == nullptr)
            {
                gmx_fatal(FARGS, "tpA = %d, i= %d in print_atoms", tpA, i);
            }

            /* This is true by construction, but static analysers don't know */
            GMX_ASSERT(!bRTPresname || at->resinfo[at->atom[i].resind].rtp,
                       "-rtpres did not have residue name available");
            fprintf(out, "%6d %10s %6d%c %5s %6s %6d %10g %10g", i + 1, tpnmA, at->resinfo[ri].nr,
                    at->resinfo[ri].ic,
                    bRTPresname ? *(at->resinfo[at->atom[i].resind].rtp)
                                : *(at->resinfo[at->atom[i].resind].name),
                    *(at->atomname[i]), cgnr[i], at->atom[i].q, at->atom[i].m);
            if (PERTURBED(at->atom[i]))
            {
                tpB = at->atom[i].typeB;
                if ((tpnmB = atype->atomNameFromAtomType(tpB)) == nullptr)
                {
                    gmx_fatal(FARGS, "tpB = %d, i= %d in print_atoms", tpB, i);
                }
                fprintf(out, " %6s %10g %10g", tpnmB, at->atom[i].qB, at->atom[i].mB);
            }
            // Accumulate the total charge to help troubleshoot issues.
            qtot += static_cast<double>(at->atom[i].q);
            // Round it to zero if it is close to zero, because
            // printing -9.34e-5 confuses users.
            if (fabs(qtot) < 0.0001)
            {
                qtot = 0;
            }
            // Write the total charge for the last atom of the system
            // and/or residue, because generally that's where it is
            // expected to be an integer.
            if (i == at->nr - 1 || ri != at->atom[i + 1].resind)
            {
                fprintf(out, "   ; qtot %.4g\n", qtot);
            }
            else
            {
                fputs("\n", out);
            }
        }
    }
    fprintf(out, "\n");
    fflush(out);
}

void print_bondeds(FILE* out, int natoms, Directive d, int ftype, int fsubtype, gmx::ArrayRef<const InteractionsOfType> plist)
{
    t_symtab stab;
    t_atom*  a;

    PreprocessingAtomTypes atype;
    snew(a, 1);
    open_symtab(&stab);
    for (int i = 0; (i < natoms); i++)
    {
        char buf[12];
        sprintf(buf, "%4d", (i + 1));
        atype.addType(&stab, *a, buf, InteractionOfType({}, {}), 0, 0);
    }
    print_bt(out, d, &atype, ftype, fsubtype, plist, TRUE);

    done_symtab(&stab);
    sfree(a);
}
