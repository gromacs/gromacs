/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include <cstdio>
#include <cstring>

#include <filesystem>
#include <string>
#include <vector>

#include "gromacs/gmxana/gstat.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"

std::vector<t_dlist> mk_dlist(FILE*          log,
                              const t_atoms* atoms,
                              gmx_bool       bPhi,
                              gmx_bool       bPsi,
                              gmx_bool       bChi,
                              gmx_bool       bHChi,
                              int            maxchi,
                              int            r0,
                              int            rN)
{
    int       i, j, ii;
    t_dihatms atm, prev;
    int       nl = 0, nc[edMax];
    char*     thisres;
    // Initially, size this to all possible residues. Later it might
    // be reduced to only handle those residues identified to contain
    // dihedrals.
    std::vector<t_dlist> dl(atoms->nres + 1);

    /* Enforcing indexing rules for last residue */
    if (rN != -1 && rN < r0)
    {
        gmx_fatal(
                FARGS,
                "Ending residue index (-rN) must be greater than beginning residue index (-r0).\n");
    }
    if (rN == -1)
    {
        /* rN is same default as size of dl */
        rN = atoms->nres + 1;
    }

    fprintf(stderr, "Analyzing from residue %d to residue %d", r0, rN);

    prev.C = prev.Cn[1] = -1; /* Keep the compiler quiet */
    for (i = 0; (i < edMax); i++)
    {
        nc[i] = 0;
    }
    i = 0;
    while (i < atoms->nr)
    {
        int ires = atoms->atom[i].resind;

        /* Initiate all atom numbers to -1 */
        atm.minC = atm.maxN = atm.H = atm.N = atm.C = atm.O = atm.minCalpha = -1;
        for (j = 0; (j < MAXCHI + 3); j++)
        {
            atm.Cn[j] = -1;
        }

        /* Look for atoms in this residue */
        /* maybe should allow for chis to hydrogens? */
        while ((i < atoms->nr) && (atoms->atom[i].resind == ires))
        {
            /* look at residues within the region specified by
               r0 - 1 and rN + 1 (includes previous and subsequent
               residues for phi/psi/omega calculations). This does not affect
               default behavior. Use ires + 1 because atoms->atom[i].resind
               starts at 0. */
            if ((ires + 1 >= r0 - 1) && (ires + 1 <= rN + 1))
            {
                if ((std::strcmp(*(atoms->atomname[i]), "H") == 0)
                    || (std::strcmp(*(atoms->atomname[i]), "H1") == 0)
                    || (std::strcmp(*(atoms->atomname[i]), "HN") == 0))
                {
                    atm.H = i;
                }
                else if (std::strcmp(*(atoms->atomname[i]), "N") == 0)
                {
                    atm.N = i;
                    /* specify N for psi of previous residue */
                    if (nl != 0)
                    {
                        dl[nl - 1].atm.maxN = i;
                    }
                }
                else if (std::strcmp(*(atoms->atomname[i]), "C") == 0)
                {
                    atm.C = i;
                    /* specify prev.C for nondefault r0 (phi & omega) */
                    if (ires + 1 == r0 - 1)
                    {
                        prev.C = i;
                    }
                }
                else if ((std::strcmp(*(atoms->atomname[i]), "O") == 0)
                         || (std::strcmp(*(atoms->atomname[i]), "O1") == 0)
                         || (std::strcmp(*(atoms->atomname[i]), "OC1") == 0)
                         || (std::strcmp(*(atoms->atomname[i]), "OT1") == 0))
                {
                    atm.O = i;
                }
                else if (std::strcmp(*(atoms->atomname[i]), "CA") == 0)
                {
                    atm.Cn[1] = i;
                    /* specify prev.Cn[1] (Calpha) for nondefault r0 */
                    if (ires + 1 == r0 - 1)
                    {
                        prev.Cn[1] = i;
                    }
                }
                else if (std::strcmp(*(atoms->atomname[i]), "CB") == 0)
                {
                    atm.Cn[2] = i;
                }
                else if ((std::strcmp(*(atoms->atomname[i]), "CG") == 0)
                         || (std::strcmp(*(atoms->atomname[i]), "CG1") == 0)
                         || (std::strcmp(*(atoms->atomname[i]), "OG") == 0)
                         || (std::strcmp(*(atoms->atomname[i]), "OG1") == 0)
                         || (std::strcmp(*(atoms->atomname[i]), "SG") == 0))
                {
                    atm.Cn[3] = i;
                }
                else if ((std::strcmp(*(atoms->atomname[i]), "CD") == 0)
                         || (std::strcmp(*(atoms->atomname[i]), "CD1") == 0)
                         || (std::strcmp(*(atoms->atomname[i]), "SD") == 0)
                         || (std::strcmp(*(atoms->atomname[i]), "OD1") == 0)
                         || (std::strcmp(*(atoms->atomname[i]), "ND1") == 0))
                { // NOLINT bugprone-branch-clone
                    atm.Cn[4] = i;
                }
                /* by grs - split the Cn[4] into 2 bits to check allowing dih to H */
                else if (bHChi
                         && ((std::strcmp(*(atoms->atomname[i]), "HG") == 0)
                             || (std::strcmp(*(atoms->atomname[i]), "HG1") == 0)))
                {
                    atm.Cn[4] = i;
                }
                else if ((std::strcmp(*(atoms->atomname[i]), "CE") == 0)
                         || (std::strcmp(*(atoms->atomname[i]), "CE1") == 0)
                         || (std::strcmp(*(atoms->atomname[i]), "OE1") == 0)
                         || (std::strcmp(*(atoms->atomname[i]), "NE") == 0))
                {
                    atm.Cn[5] = i;
                }
                else if ((std::strcmp(*(atoms->atomname[i]), "CZ") == 0)
                         || (std::strcmp(*(atoms->atomname[i]), "NZ") == 0))
                {
                    atm.Cn[6] = i;
                }
                /* HChi flag here too */
                else if (bHChi && (std::strcmp(*(atoms->atomname[i]), "NH1") == 0))
                {
                    atm.Cn[7] = i;
                }
            }
            i++;
        }

        thisres = *(atoms->resinfo[ires].name);

        /* added by grs - special case for aromatics, whose chis above 2 are
           not real and produce rubbish output - so set back to -1 */
        if (std::strcmp(thisres, "PHE") == 0 || std::strcmp(thisres, "TYR") == 0
            || std::strcmp(thisres, "PTR") == 0 || std::strcmp(thisres, "TRP") == 0
            || std::strcmp(thisres, "HIS") == 0 || std::strcmp(thisres, "HISA") == 0
            || std::strcmp(thisres, "HISB") == 0)
        {
            for (ii = 5; ii <= 7; ii++)
            {
                atm.Cn[ii] = -1;
            }
        }
        /* end fixing aromatics */

        /* Special case for Pro, has no H */
        if (std::strcmp(thisres, "PRO") == 0)
        {
            atm.H = atm.Cn[4];
        }
        /* Carbon from previous residue */
        if (prev.C != -1)
        {
            atm.minC = prev.C;
        }
        /* Alpha-carbon from previous residue */
        if (prev.Cn[1] != -1)
        {
            atm.minCalpha = prev.Cn[1];
        }
        prev = atm;

        /* Check how many dihedrals we have. Use ires + 1 because
           atoms->atom[i].resind starts at 0. Only consider residues
           within range defined by r0 and rN in this case. */
        if ((ires + 1 >= r0) && (ires + 1 <= rN))
        {
            if ((atm.N != -1) && (atm.Cn[1] != -1) && (atm.C != -1) && (atm.O != -1)
                && ((atm.H != -1) || (atm.minC != -1)))
            {
                dl[nl].resnr     = ires + 1;
                dl[nl].atm       = atm;
                dl[nl].atm.Cn[0] = atm.N;
                if ((atm.Cn[3] != -1) && (atm.Cn[2] != -1) && (atm.Cn[1] != -1))
                {
                    nc[0]++;
                    if (atm.Cn[4] != -1)
                    {
                        nc[1]++;
                        if (atm.Cn[5] != -1)
                        {
                            nc[2]++;
                            if (atm.Cn[6] != -1)
                            {
                                nc[3]++;
                                if (atm.Cn[7] != -1)
                                {
                                    nc[4]++;
                                    if (atm.Cn[8] != -1)
                                    {
                                        nc[5]++;
                                    }
                                }
                            }
                        }
                    }
                }
                if ((atm.minC != -1) && (atm.minCalpha != -1))
                {
                    nc[6]++;
                }

                dl[nl].residueName = thisres;

                /* ires + 1 to correctly number residue files
                   (atoms->atom[i].resind starts at zero)     */
                sprintf(dl[nl].name, "%s%d", thisres, ires + 1);
                nl++;
            }
        }
        else if (debug)
        {
            fprintf(debug,
                    "Could not find N atom but could find other atoms"
                    " in residue %s%d\n",
                    thisres,
                    ires + 1);
        }
    }
    // Leave only the residues that were recognized to contain dihedrals
    dl.resize(nl);

    fprintf(stderr, "\n");
    fprintf(log, "\n");
    fprintf(log, "There are %d residues with dihedrals\n", nl);
    j = 0;
    if (bPhi)
    {
        j += nl;
    }
    if (bPsi)
    {
        j += nl;
    }
    if (bChi)
    {
        for (i = 0; (i < maxchi); i++)
        {
            j += nc[i];
        }
    }
    fprintf(log, "There are %d dihedrals\n", j);
    fprintf(log, "Dihedral: ");
    if (bPhi)
    {
        fprintf(log, " Phi  ");
    }
    if (bPsi)
    {
        fprintf(log, " Psi  ");
    }
    if (bChi)
    {
        for (i = 0; (i < maxchi); i++)
        {
            fprintf(log, "Chi%d  ", i + 1);
        }
    }
    fprintf(log, "\nNumber:   ");
    if (bPhi)
    {
        fprintf(log, "%4d  ", nl);
    }
    if (bPsi)
    {
        fprintf(log, "%4d  ", nl);
    }
    if (bChi)
    {
        for (i = 0; (i < maxchi); i++)
        {
            fprintf(log, "%4d  ", nc[i]);
        }
    }
    fprintf(log, "\n");

    return dl;
}

gmx_bool has_dihedral(int Dih, const t_dlist& dl)
{
    gmx_bool b = FALSE;
    int      ddd;

    switch (Dih)
    {
        case edPhi:
            b = ((dl.atm.H != -1) && (dl.atm.N != -1) && (dl.atm.Cn[1] != -1) && (dl.atm.C != -1));
            break;
        case edPsi:
            b = ((dl.atm.N != -1) && (dl.atm.Cn[1] != -1) && (dl.atm.C != -1) && (dl.atm.O != -1));
            break;
        case edOmega:
            b = ((dl.atm.minCalpha != -1) && (dl.atm.minC != -1) && (dl.atm.N != -1)
                 && (dl.atm.Cn[1] != -1));
            break;
        case edChi1:
        case edChi2:
        case edChi3:
        case edChi4:
        case edChi5:
        case edChi6:
            ddd = Dih - edChi1;
            b = ((dl.atm.Cn[ddd] != -1) && (dl.atm.Cn[ddd + 1] != -1) && (dl.atm.Cn[ddd + 2] != -1)
                 && (dl.atm.Cn[ddd + 3] != -1));
            break;
        default:
            pr_dlist(stdout, gmx::constArrayRefFromArray(&dl, 1), 1, 0, TRUE, TRUE, TRUE, TRUE, MAXCHI);
            gmx_fatal(FARGS, "Non existent dihedral %d in file %s, line %d", Dih, __FILE__, __LINE__);
    }
    return b;
}

static void pr_one_ro(FILE* fp, const t_dlist& dl, int nDih, real gmx_unused dt)
{
    int k;
    for (k = 0; k < NROT; k++)
    {
        fprintf(fp, "  %6.2f", dl.rot_occ[nDih][k]);
    }
    fprintf(fp, "\n");
}

static void pr_ntr_s2(FILE* fp, const t_dlist& dl, int nDih, real dt)
{
    fprintf(fp, "  %6.2f  %6.2f\n", (dt == 0) ? 0 : dl.ntr[nDih] / dt, dl.S2[nDih]);
}

void pr_dlist(FILE*                        fp,
              gmx::ArrayRef<const t_dlist> dlist,
              real                         dt,
              int                          printtype,
              gmx_bool                     bPhi,
              gmx_bool                     bPsi,
              gmx_bool                     bChi,
              gmx_bool                     bOmega,
              int                          maxchi)
{
    void (*pr_props)(FILE*, const t_dlist&, int, real);

    /* Analysis of dihedral transitions etc */

    if (printtype == edPrintST)
    {
        pr_props = pr_ntr_s2;
        fprintf(stderr, "Now printing out transitions and OPs...\n");
    }
    else
    {
        pr_props = pr_one_ro;
        fprintf(stderr, "Now printing out rotamer occupancies...\n");
        fprintf(fp, "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n");
    }

    /* change atom numbers from 0 based to 1 based */
    for (const auto& dihedral : dlist)
    {
        fprintf(fp, "Residue %s\n", dihedral.name);
        if (printtype == edPrintST)
        {
            fprintf(fp,
                    " Angle [   AI,   AJ,   AK,   AL]  #tr/ns  S^2D  \n"
                    "--------------------------------------------\n");
        }
        else
        {
            fprintf(fp,
                    " Angle [   AI,   AJ,   AK,   AL]  rotamers  0  g(-)  t  g(+)\n"
                    "--------------------------------------------\n");
        }
        if (bPhi)
        {
            fprintf(fp,
                    "   Phi [%5d,%5d,%5d,%5d]",
                    (dihedral.atm.H == -1) ? 1 + dihedral.atm.minC : 1 + dihedral.atm.H,
                    1 + dihedral.atm.N,
                    1 + dihedral.atm.Cn[1],
                    1 + dihedral.atm.C);
            pr_props(fp, dihedral, edPhi, dt);
        }
        if (bPsi)
        {
            fprintf(fp,
                    "   Psi [%5d,%5d,%5d,%5d]",
                    1 + dihedral.atm.N,
                    1 + dihedral.atm.Cn[1],
                    1 + dihedral.atm.C,
                    1 + dihedral.atm.O);
            pr_props(fp, dihedral, edPsi, dt);
        }
        if (bOmega && has_dihedral(edOmega, dihedral))
        {
            fprintf(fp,
                    " Omega [%5d,%5d,%5d,%5d]",
                    1 + dihedral.atm.minCalpha,
                    1 + dihedral.atm.minC,
                    1 + dihedral.atm.N,
                    1 + dihedral.atm.Cn[1]);
            pr_props(fp, dihedral, edOmega, dt);
        }
        for (int Xi = 0; Xi < MAXCHI; Xi++)
        {
            if (bChi && (Xi < maxchi) && (dihedral.atm.Cn[Xi + 3] != -1))
            {
                fprintf(fp,
                        "   Chi%d[%5d,%5d,%5d,%5d]",
                        Xi + 1,
                        1 + dihedral.atm.Cn[Xi],
                        1 + dihedral.atm.Cn[Xi + 1],
                        1 + dihedral.atm.Cn[Xi + 2],
                        1 + dihedral.atm.Cn[Xi + 3]);
                pr_props(fp, dihedral, Xi + edChi1, dt); /* Xi+2 was wrong here */
            }
        }
        fprintf(fp, "\n");
    }
}
