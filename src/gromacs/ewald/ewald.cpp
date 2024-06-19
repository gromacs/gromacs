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

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <array>
#include <filesystem>

#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

using cvec = std::array<t_complex, DIM>;

gmx_ewald_tab_t::gmx_ewald_tab_t(const t_inputrec& ir, FILE* fp)
{
    if (fp)
    {
        fprintf(fp, "Will do ordinary reciprocal space Ewald sum.\n");
    }

    nx   = ir.nkx + 1;
    ny   = ir.nky + 1;
    nz   = ir.nkz + 1;
    kmax = std::max(nx, std::max(ny, nz));
}

gmx_ewald_tab_t::~gmx_ewald_tab_t() = default;

//! Calculates wave vectors.
static void calc_lll(const rvec box, rvec lll)
{
    lll[XX] = 2.0 * M_PI / box[XX];
    lll[YY] = 2.0 * M_PI / box[YY];
    lll[ZZ] = 2.0 * M_PI / box[ZZ];
}

//! Make tables for the structure factor parts
static void tabulateStructureFactors(int natom, gmx::ArrayRef<const gmx::RVec> x, int kmax, cvec** eir, const rvec lll)
{
    int i, j, m;

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
            eir[1][i][m].re = std::cos(x[i][m] * lll[m]);
            eir[1][i][m].im = std::sin(x[i][m] * lll[m]);
        }
        for (j = 2; (j < kmax); j++)
        {
            for (m = 0; (m < 3); m++)
            {
                eir[j][i][m] = cmul(eir[j - 1][i][m], eir[1][i][m]);
            }
        }
    }
}

real do_ewald(bool                           havePbcXY2Walls,
              real                           wallEwaldZfac,
              real                           epsilonR,
              FreeEnergyPerturbationType     freeEnergyPerturbationType,
              gmx::ArrayRef<const gmx::RVec> coords,
              gmx::ArrayRef<gmx::RVec>       forces,
              gmx::ArrayRef<const real>      chargeA,
              gmx::ArrayRef<const real>      chargeB,
              const matrix                   box,
              const t_commrec*               commrec,
              int                            natoms,
              matrix                         lrvir,
              real                           ewaldcoeff,
              real                           lambda,
              real*                          dvdlambda,
              gmx_ewald_tab_t*               et)
{
    real   factor = -1.0 / (4 * ewaldcoeff * ewaldcoeff);
    real   energy_AB[2], energy;
    rvec   lll;
    int    lowiy, lowiz, ix, iy, iz, n, q;
    real   tmp, cs, ss, ak, akv, mx, my, mz, m2, scale;
    cvec** eir;
    bool   bFreeEnergy;

    if (commrec != nullptr)
    {
        if (PAR(commrec))
        {
            gmx_fatal(FARGS, "No parallel Ewald. Use PME instead.\n");
        }
    }

    /* Scale box with Ewald wall factor */
    matrix          scaledBox;
    EwaldBoxZScaler boxScaler(havePbcXY2Walls, wallEwaldZfac);
    boxScaler.scaleBox(box, scaledBox);

    rvec boxDiag;
    for (int i = 0; (i < DIM); i++)
    {
        boxDiag[i] = scaledBox[i][i];
    }

    /* 1/(Vol*e0) */
    real scaleRecip = 4.0 * M_PI / (boxDiag[XX] * boxDiag[YY] * boxDiag[ZZ]) * gmx::c_one4PiEps0 / epsilonR;

    snew(eir, et->kmax);
    for (n = 0; n < et->kmax; n++)
    {
        snew(eir[n], natoms);
    }
    et->tab_xy.resize(natoms);
    et->tab_qxyz.resize(natoms);

    bFreeEnergy = (freeEnergyPerturbationType != FreeEnergyPerturbationType::No);

    clear_mat(lrvir);

    calc_lll(boxDiag, lll);
    tabulateStructureFactors(natoms, coords, et->kmax, eir, lll);

    gmx::ArrayRef<const real> charge;
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
            mx = ix * lll[XX];
            for (iy = lowiy; iy < et->ny; iy++)
            {
                my = iy * lll[YY];
                if (iy >= 0)
                {
                    for (n = 0; n < natoms; n++)
                    {
                        et->tab_xy[n] = cmul(eir[ix][n][XX], eir[iy][n][YY]);
                    }
                }
                else
                {
                    for (n = 0; n < natoms; n++)
                    {
                        et->tab_xy[n] = cmul(eir[ix][n][XX], conjugate(eir[-iy][n][YY]));
                    }
                }
                for (iz = lowiz; iz < et->nz; iz++)
                {
                    mz  = iz * lll[ZZ];
                    m2  = mx * mx + my * my + mz * mz;
                    ak  = std::exp(m2 * factor) / m2;
                    akv = 2.0 * ak * (1.0 / m2 - factor);
                    if (iz >= 0)
                    {
                        for (n = 0; n < natoms; n++)
                        {
                            et->tab_qxyz[n] = rcmul(charge[n], cmul(et->tab_xy[n], eir[iz][n][ZZ]));
                        }
                    }
                    else
                    {
                        for (n = 0; n < natoms; n++)
                        {
                            et->tab_qxyz[n] =
                                    rcmul(charge[n], cmul(et->tab_xy[n], conjugate(eir[-iz][n][ZZ])));
                        }
                    }

                    cs = ss = 0;
                    for (n = 0; n < natoms; n++)
                    {
                        cs += et->tab_qxyz[n].re;
                        ss += et->tab_qxyz[n].im;
                    }
                    energy_AB[q] += ak * (cs * cs + ss * ss);
                    tmp = scale * akv * (cs * cs + ss * ss);
                    lrvir[XX][XX] -= tmp * mx * mx;
                    lrvir[XX][YY] -= tmp * mx * my;
                    lrvir[XX][ZZ] -= tmp * mx * mz;
                    lrvir[YY][YY] -= tmp * my * my;
                    lrvir[YY][ZZ] -= tmp * my * mz;
                    lrvir[ZZ][ZZ] -= tmp * mz * mz;
                    for (n = 0; n < natoms; n++)
                    {
                        /*tmp=scale*ak*(cs*tab_qxyz[n].im-ss*tab_qxyz[n].re);*/
                        tmp = scale * ak * (cs * et->tab_qxyz[n].im - ss * et->tab_qxyz[n].re);
                        forces[n][XX] += tmp * mx * 2 * scaleRecip;
                        forces[n][YY] += tmp * my * 2 * scaleRecip;
                        forces[n][ZZ] += tmp * mz * 2 * scaleRecip;
                    }
                    lowiz = 1 - et->nz;
                }
                lowiy = 1 - et->ny;
            }
        }
    }

    if (!bFreeEnergy)
    {
        energy = energy_AB[0];
    }
    else
    {
        energy = (1.0 - lambda) * energy_AB[0] + lambda * energy_AB[1];
        *dvdlambda += scaleRecip * (energy_AB[1] - energy_AB[0]);
    }

    lrvir[XX][XX] = -0.5 * scaleRecip * (lrvir[XX][XX] + energy);
    lrvir[XX][YY] = -0.5 * scaleRecip * (lrvir[XX][YY]);
    lrvir[XX][ZZ] = -0.5 * scaleRecip * (lrvir[XX][ZZ]);
    lrvir[YY][YY] = -0.5 * scaleRecip * (lrvir[YY][YY] + energy);
    lrvir[YY][ZZ] = -0.5 * scaleRecip * (lrvir[YY][ZZ]);
    lrvir[ZZ][ZZ] = -0.5 * scaleRecip * (lrvir[ZZ][ZZ] + energy);

    lrvir[YY][XX] = lrvir[XX][YY];
    lrvir[ZZ][XX] = lrvir[XX][ZZ];
    lrvir[ZZ][YY] = lrvir[YY][ZZ];

    energy *= scaleRecip;

    return energy;
}

real ewald_charge_correction(const t_commrec*            commrec,
                             const real                  epsilonR,
                             const real                  ewaldcoeffQ,
                             gmx::ArrayRef<const double> qsum,
                             const real                  lambda,
                             const matrix                box,
                             real*                       dvdlambda,
                             tensor                      vir)

{
    real enercorr = 0;

    if (MAIN(commrec))
    {
        /* Apply charge correction */
        real vol = box[XX][XX] * box[YY][YY] * box[ZZ][ZZ];

        real fac = M_PI * gmx::c_one4PiEps0 / (epsilonR * 2.0 * vol * vol * gmx::square(ewaldcoeffQ));

        real qs2A = qsum[0] * qsum[0];
        real qs2B = qsum[1] * qsum[1];

        real vc = (qs2A * (1 - lambda) + qs2B * lambda) * fac;

        enercorr = -vol * vc;

        *dvdlambda += -vol * (qs2B - qs2A) * fac;

        for (int d = 0; d < DIM; d++)
        {
            vir[d][d] += vc;
        }
    }

    return enercorr;
}
