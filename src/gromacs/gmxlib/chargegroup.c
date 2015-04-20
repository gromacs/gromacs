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

#include "gromacs/legacyheaders/chargegroup.h"

#include <math.h>

#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"


void calc_chargegroup_radii(const gmx_mtop_t *mtop, rvec *x,
                            real *rvdw1, real *rvdw2,
                            real *rcoul1, real *rcoul2)
{
    real            r2v1, r2v2, r2c1, r2c2, r2;
    int             ntype, i, j, mb, m, cg, a_mol, a0, a1, a;
    gmx_bool       *bLJ;
    gmx_molblock_t *molb;
    gmx_moltype_t  *molt;
    t_block        *cgs;
    t_atom         *atom;
    rvec            cen;

    r2v1 = 0;
    r2v2 = 0;
    r2c1 = 0;
    r2c2 = 0;

    ntype = mtop->ffparams.atnr;
    snew(bLJ, ntype);
    for (i = 0; i < ntype; i++)
    {
        bLJ[i] = FALSE;
        for (j = 0; j < ntype; j++)
        {
            if (mtop->ffparams.iparams[i*ntype+j].lj.c6  != 0 ||
                mtop->ffparams.iparams[i*ntype+j].lj.c12 != 0)
            {
                bLJ[i] = TRUE;
            }
        }
    }

    a_mol = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molb = &mtop->molblock[mb];
        molt = &mtop->moltype[molb->type];
        cgs  = &molt->cgs;
        atom = molt->atoms.atom;
        for (m = 0; m < molb->nmol; m++)
        {
            for (cg = 0; cg < cgs->nr; cg++)
            {
                a0 = cgs->index[cg];
                a1 = cgs->index[cg+1];
                if (a1 - a0 > 1)
                {
                    clear_rvec(cen);
                    for (a = a0; a < a1; a++)
                    {
                        rvec_inc(cen, x[a_mol+a]);
                    }
                    svmul(1.0/(a1-a0), cen, cen);
                    for (a = a0; a < a1; a++)
                    {
                        r2 = distance2(cen, x[a_mol+a]);
                        if (r2 > r2v2 && (bLJ[atom[a].type ] ||
                                          bLJ[atom[a].typeB]))
                        {
                            if (r2 > r2v1)
                            {
                                r2v2 = r2v1;
                                r2v1 = r2;
                            }
                            else
                            {
                                r2v2 = r2;
                            }
                        }
                        if (r2 > r2c2 &&
                            (atom[a].q != 0 || atom[a].qB != 0))
                        {
                            if (r2 > r2c1)
                            {
                                r2c2 = r2c1;
                                r2c1 = r2;
                            }
                            else
                            {
                                r2c2 = r2;
                            }
                        }
                    }
                }
            }
            a_mol += molb->natoms_mol;
        }
    }

    sfree(bLJ);

    *rvdw1  = sqrt(r2v1);
    *rvdw2  = sqrt(r2v2);
    *rcoul1 = sqrt(r2c1);
    *rcoul2 = sqrt(r2c2);
}

void calc_cgcm(FILE gmx_unused *fplog, int cg0, int cg1, t_block *cgs,
               rvec pos[], rvec cg_cm[])
{
    int      icg, k, k0, k1, d;
    rvec     cg;
    real     nrcg, inv_ncg;
    atom_id *cgindex;

#ifdef DEBUG
    fprintf(fplog, "Calculating centre of geometry for charge groups %d to %d\n",
            cg0, cg1);
#endif
    cgindex = cgs->index;

    /* Compute the center of geometry for all charge groups */
    for (icg = cg0; (icg < cg1); icg++)
    {
        k0      = cgindex[icg];
        k1      = cgindex[icg+1];
        nrcg    = k1-k0;
        if (nrcg == 1)
        {
            copy_rvec(pos[k0], cg_cm[icg]);
        }
        else
        {
            inv_ncg = 1.0/nrcg;

            clear_rvec(cg);
            for (k = k0; (k < k1); k++)
            {
                for (d = 0; (d < DIM); d++)
                {
                    cg[d] += pos[k][d];
                }
            }
            for (d = 0; (d < DIM); d++)
            {
                cg_cm[icg][d] = inv_ncg*cg[d];
            }
        }
    }
}

void put_charge_groups_in_box(FILE gmx_unused *fplog, int cg0, int cg1,
                              int ePBC, matrix box, t_block *cgs,
                              rvec pos[], rvec cg_cm[])

{
    int      npbcdim, icg, k, k0, k1, d, e;
    rvec     cg;
    real     nrcg, inv_ncg;
    atom_id *cgindex;
    gmx_bool bTric;

    if (ePBC == epbcNONE)
    {
        gmx_incons("Calling put_charge_groups_in_box for a system without PBC");
    }

#ifdef DEBUG
    fprintf(fplog, "Putting cgs %d to %d in box\n", cg0, cg1);
#endif
    cgindex = cgs->index;

    if (ePBC == epbcXY)
    {
        npbcdim = 2;
    }
    else
    {
        npbcdim = 3;
    }

    bTric = TRICLINIC(box);

    for (icg = cg0; (icg < cg1); icg++)
    {
        /* First compute the center of geometry for this charge group */
        k0      = cgindex[icg];
        k1      = cgindex[icg+1];
        nrcg    = k1-k0;

        if (nrcg == 1)
        {
            copy_rvec(pos[k0], cg_cm[icg]);
        }
        else
        {
            inv_ncg = 1.0/nrcg;

            clear_rvec(cg);
            for (k = k0; (k < k1); k++)
            {
                for (d = 0; d < DIM; d++)
                {
                    cg[d] += pos[k][d];
                }
            }
            for (d = 0; d < DIM; d++)
            {
                cg_cm[icg][d] = inv_ncg*cg[d];
            }
        }
        /* Now check pbc for this cg */
        if (bTric)
        {
            for (d = npbcdim-1; d >= 0; d--)
            {
                while (cg_cm[icg][d] < 0)
                {
                    for (e = d; e >= 0; e--)
                    {
                        cg_cm[icg][e] += box[d][e];
                        for (k = k0; (k < k1); k++)
                        {
                            pos[k][e] += box[d][e];
                        }
                    }
                }
                while (cg_cm[icg][d] >= box[d][d])
                {
                    for (e = d; e >= 0; e--)
                    {
                        cg_cm[icg][e] -= box[d][e];
                        for (k = k0; (k < k1); k++)
                        {
                            pos[k][e] -= box[d][e];
                        }
                    }
                }
            }
        }
        else
        {
            for (d = 0; d < npbcdim; d++)
            {
                while (cg_cm[icg][d] < 0)
                {
                    cg_cm[icg][d] += box[d][d];
                    for (k = k0; (k < k1); k++)
                    {
                        pos[k][d] += box[d][d];
                    }
                }
                while (cg_cm[icg][d] >= box[d][d])
                {
                    cg_cm[icg][d] -= box[d][d];
                    for (k = k0; (k < k1); k++)
                    {
                        pos[k][d] -= box[d][d];
                    }
                }
            }
        }
#ifdef DEBUG_PBC
        for (d = 0; (d < npbcdim); d++)
        {
            if ((cg_cm[icg][d] < 0) || (cg_cm[icg][d] >= box[d][d]))
            {
                gmx_fatal(FARGS, "cg_cm[%d] = %15f  %15f  %15f\n"
                          "box = %15f  %15f  %15f\n",
                          icg, cg_cm[icg][XX], cg_cm[icg][YY], cg_cm[icg][ZZ],
                          box[XX][XX], box[YY][YY], box[ZZ][ZZ]);
            }
        }
#endif
    }
}
