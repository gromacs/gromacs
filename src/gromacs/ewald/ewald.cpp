/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
/*! \internal \file
 *
 * \brief This file contains function definitions necessary for
 * computing energies and forces for the plain-Ewald long-ranged part,
 * and the correction for overall system charge for all Ewald-family
 * methods.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "ewald.h"

#include <math.h>
#include <stdio.h>

#include <cstdlib>

#include <algorithm>

#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/legacyheaders/types/inputrec.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/smalloc.h"

struct gmx_ewald_tab_t
{
    int        nx, ny, nz, kmax;
    cvec     **eir;
    t_complex *tab_xy, *tab_qxyz;
};

void init_ewald_tab(struct gmx_ewald_tab_t **et, const t_inputrec *ir, FILE *fp)
{
    snew(*et, 1);
    if (fp)
    {
        fprintf(fp, "Will do ordinary reciprocal space Ewald sum.\n");
    }

    (*et)->nx       = ir->nkx+1;
    (*et)->ny       = ir->nky+1;
    (*et)->nz       = ir->nkz+1;
    (*et)->kmax     = std::max((*et)->nx, std::max((*et)->ny, (*et)->nz));
    (*et)->eir      = NULL;
    (*et)->tab_xy   = NULL;
    (*et)->tab_qxyz = NULL;
}

//! Make tables for the structure factor parts
static void tabulateStructureFactors(int natom, rvec x[], int kmax, cvec **eir, rvec lll)
{
    int  i, j, m;

    if (kmax < 1)
    {
        printf("Go away! kmax = %d\n", kmax);
        exit(1);
    }

    for (i = 0; (i < natom); i++)
    {
        for (m = 0; (m < 3); m++)
        {
            eir[0][i][m].re = 1;
            eir[0][i][m].im = 0;
        }

        for (m = 0; (m < 3); m++)
        {
            eir[1][i][m].re = cos(x[i][m]*lll[m]);
            eir[1][i][m].im = sin(x[i][m]*lll[m]);
        }
        for (j = 2; (j < kmax); j++)
        {
            for (m = 0; (m < 3); m++)
            {
                eir[j][i][m] = cmul(eir[j-1][i][m], eir[1][i][m]);
            }
        }
    }
}

real do_ewald(t_inputrec *ir,
              rvec x[],        rvec f[],
              real chargeA[],  real chargeB[],
              rvec box,
              t_commrec *cr,   int natoms,
              matrix lrvir,    real ewaldcoeff,
              real lambda,     real *dvdlambda,
              struct gmx_ewald_tab_t *et)
{
    real     factor     = -1.0/(4*ewaldcoeff*ewaldcoeff);
    real     scaleRecip = 4.0*M_PI/(box[XX]*box[YY]*box[ZZ])*ONE_4PI_EPS0/ir->epsilon_r; /* 1/(Vol*e0) */
    real    *charge, energy_AB[2], energy;
    rvec     lll;
    int      lowiy, lowiz, ix, iy, iz, n, q;
    real     tmp, cs, ss, ak, akv, mx, my, mz, m2, scale;
    gmx_bool bFreeEnergy;

    if (cr != NULL)
    {
        if (PAR(cr))
        {
            gmx_fatal(FARGS, "No parallel Ewald. Use PME instead.\n");
        }
    }


    if (!et->eir) /* allocate if we need to */
    {
        snew(et->eir, et->kmax);
        for (n = 0; n < et->kmax; n++)
        {
            snew(et->eir[n], natoms);
        }
        snew(et->tab_xy, natoms);
        snew(et->tab_qxyz, natoms);
    }

    bFreeEnergy = (ir->efep != efepNO);

    clear_mat(lrvir);

    calc_lll(box, lll);
    tabulateStructureFactors(natoms, x, et->kmax, et->eir, lll);

    for (q = 0; q < (bFreeEnergy ? 2 : 1); q++)
    {
        if (!bFreeEnergy)
        {
            charge = chargeA;
            scale  = 1.0;
        }
        else if (q == 0)
        {
            charge = chargeA;
            scale  = 1.0 - lambda;
        }
        else
        {
            charge = chargeB;
            scale  = lambda;
        }
        lowiy        = 0;
        lowiz        = 1;
        energy_AB[q] = 0;
        for (ix = 0; ix < et->nx; ix++)
        {
            mx = ix*lll[XX];
            for (iy = lowiy; iy < et->ny; iy++)
            {
                my = iy*lll[YY];
                if (iy >= 0)
                {
                    for (n = 0; n < natoms; n++)
                    {
                        et->tab_xy[n] = cmul(et->eir[ix][n][XX], et->eir[iy][n][YY]);
                    }
                }
                else
                {
                    for (n = 0; n < natoms; n++)
                    {
                        et->tab_xy[n] = cmul(et->eir[ix][n][XX], conjugate(et->eir[-iy][n][YY]));
                    }
                }
                for (iz = lowiz; iz < et->nz; iz++)
                {
                    mz  = iz*lll[ZZ];
                    m2  = mx*mx+my*my+mz*mz;
                    ak  = exp(m2*factor)/m2;
                    akv = 2.0*ak*(1.0/m2-factor);
                    if (iz >= 0)
                    {
                        for (n = 0; n < natoms; n++)
                        {
                            et->tab_qxyz[n] = rcmul(charge[n], cmul(et->tab_xy[n],
                                                                    et->eir[iz][n][ZZ]));
                        }
                    }
                    else
                    {
                        for (n = 0; n < natoms; n++)
                        {
                            et->tab_qxyz[n] = rcmul(charge[n], cmul(et->tab_xy[n],
                                                                    conjugate(et->eir[-iz][n][ZZ])));
                        }
                    }

                    cs = ss = 0;
                    for (n = 0; n < natoms; n++)
                    {
                        cs += et->tab_qxyz[n].re;
                        ss += et->tab_qxyz[n].im;
                    }
                    energy_AB[q]  += ak*(cs*cs+ss*ss);
                    tmp            = scale*akv*(cs*cs+ss*ss);
                    lrvir[XX][XX] -= tmp*mx*mx;
                    lrvir[XX][YY] -= tmp*mx*my;
                    lrvir[XX][ZZ] -= tmp*mx*mz;
                    lrvir[YY][YY] -= tmp*my*my;
                    lrvir[YY][ZZ] -= tmp*my*mz;
                    lrvir[ZZ][ZZ] -= tmp*mz*mz;
                    for (n = 0; n < natoms; n++)
                    {
                        /*tmp=scale*ak*(cs*tab_qxyz[n].im-ss*tab_qxyz[n].re);*/
                        tmp       = scale*ak*(cs*et->tab_qxyz[n].im-ss*et->tab_qxyz[n].re);
                        f[n][XX] += tmp*mx*2*scaleRecip;
                        f[n][YY] += tmp*my*2*scaleRecip;
                        f[n][ZZ] += tmp*mz*2*scaleRecip;
#if 0
                        f[n][XX] += tmp*mx;
                        f[n][YY] += tmp*my;
                        f[n][ZZ] += tmp*mz;
#endif
                    }
                    lowiz = 1-et->nz;
                }
                lowiy = 1-et->ny;
            }
        }
    }

    if (!bFreeEnergy)
    {
        energy = energy_AB[0];
    }
    else
    {
        energy      = (1.0 - lambda)*energy_AB[0] + lambda*energy_AB[1];
        *dvdlambda += scaleRecip*(energy_AB[1] - energy_AB[0]);
    }

    lrvir[XX][XX] = -0.5*scaleRecip*(lrvir[XX][XX]+energy);
    lrvir[XX][YY] = -0.5*scaleRecip*(lrvir[XX][YY]);
    lrvir[XX][ZZ] = -0.5*scaleRecip*(lrvir[XX][ZZ]);
    lrvir[YY][YY] = -0.5*scaleRecip*(lrvir[YY][YY]+energy);
    lrvir[YY][ZZ] = -0.5*scaleRecip*(lrvir[YY][ZZ]);
    lrvir[ZZ][ZZ] = -0.5*scaleRecip*(lrvir[ZZ][ZZ]+energy);

    lrvir[YY][XX] = lrvir[XX][YY];
    lrvir[ZZ][XX] = lrvir[XX][ZZ];
    lrvir[ZZ][YY] = lrvir[YY][ZZ];

    energy *= scaleRecip;

    return energy;
}

real ewald_charge_correction(t_commrec *cr, t_forcerec *fr, real lambda,
                             matrix box,
                             real *dvdlambda, tensor vir)

{
    real vol, fac, qs2A, qs2B, vc, enercorr;
    int  d;

    if (MASTER(cr))
    {
        /* Apply charge correction */
        vol = box[XX][XX]*box[YY][YY]*box[ZZ][ZZ];

        fac = M_PI*ONE_4PI_EPS0/(fr->epsilon_r*2.0*vol*vol*sqr(fr->ewaldcoeff_q));

        qs2A = fr->qsum[0]*fr->qsum[0];
        qs2B = fr->qsum[1]*fr->qsum[1];

        vc = (qs2A*(1 - lambda) + qs2B*lambda)*fac;

        enercorr = -vol*vc;

        *dvdlambda += -vol*(qs2B - qs2A)*fac;

        for (d = 0; d < DIM; d++)
        {
            vir[d][d] += vc;
        }

        if (debug)
        {
            fprintf(debug, "Total charge correction: Vcharge=%g\n", enercorr);
        }
    }
    else
    {
        enercorr = 0;
    }

    return enercorr;
}
