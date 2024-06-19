/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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

#include "gromacs/pbcutil/pbcmethods.h"

#include <cstdio>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <memory>

#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

enum class PbcType : int;

void calc_pbc_cluster(int ecenter, int nrefat, t_topology* top, PbcType pbcType, rvec x[], const int index[], matrix box)
{
    int       m, i, j, j0, j1, jj, ai, aj;
    int       imin, jmin;
    real      fac, min_dist2;
    rvec      dx, xtest, box_center;
    int       nmol, imol_center;
    int*      molind;
    gmx_bool *bMol, *bTmp;
    rvec *    m_com, *m_shift;
    t_pbc     pbc;
    int*      cluster;
    int*      added;
    int       ncluster, nadded;
    real      tmp_r2;

    calc_box_center(ecenter, box, box_center);

    /* Initiate the pbc structure */
    std::memset(&pbc, 0, sizeof(pbc));
    set_pbc(&pbc, pbcType, box);

    /* Convert atom index to molecular */
    nmol   = top->mols.nr;
    molind = top->mols.index;
    snew(bMol, nmol);
    snew(m_com, nmol);
    snew(m_shift, nmol);
    snew(cluster, nmol);
    snew(added, nmol);
    snew(bTmp, top->atoms.nr);

    for (i = 0; (i < nrefat); i++)
    {
        /* Mark all molecules in the index */
        ai       = index[i];
        bTmp[ai] = TRUE;
        /* Binary search assuming the molecules are sorted */
        j0 = 0;
        j1 = nmol - 1;
        while (j0 < j1)
        {
            if (ai < molind[j0 + 1])
            {
                j1 = j0;
            }
            else if (ai >= molind[j1])
            {
                j0 = j1;
            }
            else
            {
                jj = (j0 + j1) / 2;
                if (ai < molind[jj + 1])
                {
                    j1 = jj;
                }
                else
                {
                    j0 = jj;
                }
            }
        }
        bMol[j0] = TRUE;
    }
    /* Double check whether all atoms in all molecules that are marked are part
     * of the cluster. Simultaneously compute the center of geometry.
     */
    min_dist2   = 10 * gmx::square(trace(box));
    imol_center = -1;
    ncluster    = 0;
    for (i = 0; i < nmol; i++)
    {
        for (j = molind[i]; j < molind[i + 1]; j++)
        {
            if (bMol[i] && !bTmp[j])
            {
                gmx_fatal(FARGS,
                          "Molecule %d marked for clustering but not atom %d in it - check your "
                          "index!",
                          i + 1,
                          j + 1);
            }
            else if (!bMol[i] && bTmp[j])
            {
                gmx_fatal(FARGS,
                          "Atom %d marked for clustering but not molecule %d - this is an internal "
                          "error...",
                          j + 1,
                          i + 1);
            }
            else if (bMol[i])
            {
                /* Make molecule whole, move 2nd and higher atom to same periodicity as 1st atom in molecule */
                if (j > molind[i])
                {
                    pbc_dx(&pbc, x[j], x[j - 1], dx);
                    rvec_add(x[j - 1], dx, x[j]);
                }
                /* Compute center of geometry of molecule - m_com[i] was zeroed when we did snew() on it! */
                rvec_inc(m_com[i], x[j]);
            }
        }
        if (bMol[i])
        {
            /* Normalize center of geometry */
            fac = 1.0 / (molind[i + 1] - molind[i]);
            for (m = 0; (m < DIM); m++)
            {
                m_com[i][m] *= fac;
            }
            /* Determine which molecule is closest to the center of the box */
            pbc_dx(&pbc, box_center, m_com[i], dx);
            tmp_r2 = iprod(dx, dx);

            if (tmp_r2 < min_dist2)
            {
                min_dist2   = tmp_r2;
                imol_center = i;
            }
            cluster[ncluster++] = i;
        }
    }
    sfree(bTmp);

    if (ncluster <= 0)
    {
        fprintf(stderr, "No molecules selected in the cluster\n");
        return;
    }
    else if (imol_center == -1)
    {
        fprintf(stderr, "No central molecules could be found\n");
        return;
    }

    nadded            = 0;
    added[nadded++]   = imol_center;
    bMol[imol_center] = FALSE;

    while (nadded < ncluster)
    {
        /* Find min distance between cluster molecules and those remaining to be added */
        min_dist2 = 10 * gmx::square(trace(box));
        imin      = -1;
        jmin      = -1;
        /* Loop over added mols */
        for (i = 0; i < nadded; i++)
        {
            ai = added[i];
            /* Loop over all mols */
            for (j = 0; j < ncluster; j++)
            {
                aj = cluster[j];
                /* check those remaining to be added */
                if (bMol[aj])
                {
                    pbc_dx(&pbc, m_com[aj], m_com[ai], dx);
                    tmp_r2 = iprod(dx, dx);
                    if (tmp_r2 < min_dist2)
                    {
                        min_dist2 = tmp_r2;
                        imin      = ai;
                        jmin      = aj;
                    }
                }
            }
        }

        /* Add the best molecule */
        added[nadded++] = jmin;
        bMol[jmin]      = FALSE;
        /* Calculate the shift from the ai molecule */
        pbc_dx(&pbc, m_com[jmin], m_com[imin], dx);
        rvec_add(m_com[imin], dx, xtest);
        rvec_sub(xtest, m_com[jmin], m_shift[jmin]);
        rvec_inc(m_com[jmin], m_shift[jmin]);

        for (j = molind[jmin]; j < molind[jmin + 1]; j++)
        {
            rvec_inc(x[j], m_shift[jmin]);
        }
        fprintf(stdout, "\rClustering iteration %d of %d...", nadded, ncluster);
        fflush(stdout);
    }

    sfree(added);
    sfree(cluster);
    sfree(bMol);
    sfree(m_com);
    sfree(m_shift);

    fprintf(stdout, "\n");
}

void put_molecule_com_in_box(int      unitcell_enum,
                             int      ecenter,
                             t_block* mols,
                             int      natoms,
                             t_atom   atom[],
                             PbcType  pbcType,
                             matrix   box,
                             rvec     x[])
{
    int    i, j;
    int    d;
    rvec   com, shift, box_center;
    real   m;
    double mtot;
    t_pbc  pbc;

    calc_box_center(ecenter, box, box_center);
    set_pbc(&pbc, pbcType, box);
    if (mols->nr <= 0)
    {
        gmx_fatal(FARGS,
                  "There are no molecule descriptions. I need a .tpr file for this pbc option.");
    }
    for (i = 0; (i < mols->nr); i++)
    {
        /* calc COM */
        clear_rvec(com);
        mtot = 0;
        for (j = mols->index[i]; (j < mols->index[i + 1] && j < natoms); j++)
        {
            m = atom[j].m;
            for (d = 0; d < DIM; d++)
            {
                com[d] += m * x[j][d];
            }
            mtot += m;
        }
        /* calculate final COM */
        svmul(1.0 / mtot, com, com);

        /* check if COM is outside box */
        gmx::RVec newCom;
        copy_rvec(com, newCom);
        auto newComArrayRef = gmx::arrayRefFromArray(&newCom, 1);
        switch (unitcell_enum)
        {
            case euRect: put_atoms_in_box(pbcType, box, newComArrayRef); break;
            case euTric: put_atoms_in_triclinic_unitcell(ecenter, box, newComArrayRef); break;
            case euCompact:
                put_atoms_in_compact_unitcell(pbcType, ecenter, box, newComArrayRef);
                break;
        }
        rvec_sub(newCom, com, shift);
        if (norm2(shift) > 0)
        {
            if (debug)
            {
                fprintf(debug,
                        "\nShifting position of molecule %d "
                        "by %8.3f  %8.3f  %8.3f\n",
                        i + 1,
                        shift[XX],
                        shift[YY],
                        shift[ZZ]);
            }
            for (j = mols->index[i]; (j < mols->index[i + 1] && j < natoms); j++)
            {
                rvec_inc(x[j], shift);
            }
        }
    }
}

void put_residue_com_in_box(int     unitcell_enum,
                            int     ecenter,
                            int     natoms,
                            t_atom  atom[],
                            PbcType pbcType,
                            matrix  box,
                            rvec    x[])
{
    int              i, j, res_start, res_end;
    int              d, presnr;
    real             m;
    double           mtot;
    rvec             box_center, com, shift;
    static const int NOTSET = -12347;
    calc_box_center(ecenter, box, box_center);

    presnr    = NOTSET;
    res_start = 0;
    clear_rvec(com);
    mtot = 0;
    for (i = 0; i < natoms + 1; i++)
    {
        if (i == natoms || (presnr != atom[i].resind && presnr != NOTSET))
        {
            /* calculate final COM */
            res_end = i;
            svmul(1.0 / mtot, com, com);

            /* check if COM is outside box */
            gmx::RVec newCom;
            copy_rvec(com, newCom);
            auto newComArrayRef = gmx::arrayRefFromArray(&newCom, 1);
            switch (unitcell_enum)
            {
                case euRect: put_atoms_in_box(pbcType, box, newComArrayRef); break;
                case euTric: put_atoms_in_triclinic_unitcell(ecenter, box, newComArrayRef); break;
                case euCompact:
                    put_atoms_in_compact_unitcell(pbcType, ecenter, box, newComArrayRef);
                    break;
            }
            rvec_sub(newCom, com, shift);
            if (norm2(shift) != 0.0F)
            {
                if (debug)
                {
                    fprintf(debug,
                            "\nShifting position of residue %d (atoms %d-%d) "
                            "by %g,%g,%g\n",
                            atom[res_start].resind + 1,
                            res_start + 1,
                            res_end + 1,
                            shift[XX],
                            shift[YY],
                            shift[ZZ]);
                }
                for (j = res_start; j < res_end; j++)
                {
                    rvec_inc(x[j], shift);
                }
            }
            clear_rvec(com);
            mtot = 0;

            /* remember start of new residue */
            res_start = i;
        }
        if (i < natoms)
        {
            /* calc COM */
            m = atom[i].m;
            for (d = 0; d < DIM; d++)
            {
                com[d] += m * x[i][d];
            }
            mtot += m;

            presnr = atom[i].resind;
        }
    }
}

void center_x(int ecenter, rvec x[], matrix box, int n, int nc, const int ci[])
{
    int  i, m, ai;
    rvec cmin, cmax, box_center, dx;

    if (nc > 0)
    {
        copy_rvec(x[ci[0]], cmin);
        copy_rvec(x[ci[0]], cmax);
        for (i = 0; i < nc; i++)
        {
            ai = ci[i];
            for (m = 0; m < DIM; m++)
            {
                if (x[ai][m] < cmin[m])
                {
                    cmin[m] = x[ai][m];
                }
                else if (x[ai][m] > cmax[m])
                {
                    cmax[m] = x[ai][m];
                }
            }
        }
        calc_box_center(ecenter, box, box_center);
        for (m = 0; m < DIM; m++)
        {
            dx[m] = box_center[m] - (cmin[m] + cmax[m]) * 0.5;
        }

        for (i = 0; i < n; i++)
        {
            rvec_inc(x[i], dx);
        }
    }
}
