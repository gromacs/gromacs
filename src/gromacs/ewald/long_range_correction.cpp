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

#include "long_range_correction.h"

#include <cmath>

#include <array>
#include <filesystem>
#include <string>

#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

/* There's nothing special to do here if just masses are perturbed,
 * but if either charge or type is perturbed then the implementation
 * requires that B states are defined for both charge and type, and
 * does not optimize for the cases where only one changes.
 *
 * The parameter vectors for B states are left undefined in atoms2md()
 * when either FEP is inactive, or when there are no mass/charge/type
 * perturbations. The parameter vectors for LJ-PME are likewise
 * undefined when LJ-PME is not active. This works because
 * bHaveChargeOrTypePerturbed handles the control flow. */
void ewald_LRcorrection(const int                      numAtomsLocal,
                        const t_commrec*               commrec,
                        int                            numThreads,
                        int                            thread,
                        const real                     epsilonR,
                        gmx::ArrayRef<const double>    qsum,
                        EwaldGeometry                  ewaldGeometry,
                        const real                     epsilonSurface,
                        bool                           havePbcXY2Walls,
                        real                           wallEwaldZfac,
                        gmx::ArrayRef<const real>      chargeA,
                        gmx::ArrayRef<const real>      chargeB,
                        bool                           bHaveChargePerturbed,
                        gmx::ArrayRef<const gmx::RVec> coords,
                        const matrix                   box,
                        gmx::ArrayRef<const gmx::RVec> mu_tot,
                        gmx::ArrayRef<gmx::RVec>       forces,
                        real*                          Vcorr_q,
                        real                           lambda_q,
                        real*                          dvdlambda_q)
{
    /* We need to correct only self interactions */
    const int start = (numAtomsLocal * thread) / numThreads;
    const int end   = (numAtomsLocal * (thread + 1)) / numThreads;

    std::array<rvec, 2> mutot;
    gmx::RVec           dipcorrA   = { 0, 0, 0 };
    gmx::RVec           dipcorrB   = { 0, 0, 0 };
    std::array<real, 2> chargecorr = { 0 };
    std::array<real, 2> Vdipole    = { 0 };

    /* Scale the Ewald unit cell when dimension z is not periodic */
    matrix          scaledBox;
    EwaldBoxZScaler boxScaler(havePbcXY2Walls, wallEwaldZfac);
    boxScaler.scaleBox(box, scaledBox);

    const real one_4pi_eps = gmx::c_one4PiEps0 / epsilonR;
    const real L1_q        = 1.0 - lambda_q;
    /* Note that we have to transform back to gromacs units, since
     * mu_tot contains the dipole in debye units (for output).
     */
    for (int i = 0; (i < DIM); i++)
    {
        mutot[0][i] = mu_tot[0][i] * gmx::c_debye2Enm;
        mutot[1][i] = mu_tot[1][i] * gmx::c_debye2Enm;
    }
    real dipole_coeff = 0;

    real boxVolume = scaledBox[XX][XX] * scaledBox[YY][YY] * scaledBox[ZZ][ZZ];
    switch (ewaldGeometry)
    {
        case EwaldGeometry::ThreeD:
            if (epsilonSurface != 0)
            {
                dipole_coeff =
                        2 * M_PI * gmx::c_one4PiEps0 / ((2 * epsilonSurface + epsilonR) * boxVolume);
                for (int i = 0; (i < DIM); i++)
                {
                    dipcorrA[i] = 2 * dipole_coeff * mutot[0][i];
                    dipcorrB[i] = 2 * dipole_coeff * mutot[1][i];
                }
            }
            break;
        case EwaldGeometry::ThreeDC:
            dipole_coeff = 2 * M_PI * one_4pi_eps / boxVolume;
            dipcorrA[ZZ] = 2 * dipole_coeff * mutot[0][ZZ];
            dipcorrB[ZZ] = 2 * dipole_coeff * mutot[1][ZZ];
            for (int q = 0; q < (bHaveChargePerturbed ? 2 : 1); q++)
            {
                /* Avoid charge corrections with near-zero net charge */
                if (std::fabs(qsum[q]) > 1e-4)
                {
                    chargecorr[q] = 2 * dipole_coeff * qsum[q];
                }
            }
            break;
        default: gmx_incons("Unsupported Ewald geometry");
    }
    const bool bNeedLongRangeCorrection = (dipole_coeff != 0);
    if (bNeedLongRangeCorrection && !bHaveChargePerturbed)
    {
        for (int i = start; (i < end); i++)
        {
            for (int j = 0; (j < DIM); j++)
            {
                forces[i][j] -= dipcorrA[j] * chargeA[i];
            }
            if (chargecorr[0] != 0)
            {
                forces[i][ZZ] += chargecorr[0] * chargeA[i] * coords[i][ZZ];
            }
        }
    }
    else if (bNeedLongRangeCorrection)
    {
        for (int i = start; (i < end); i++)
        {
            for (int j = 0; (j < DIM); j++)
            {
                forces[i][j] -= L1_q * dipcorrA[j] * chargeA[i] + lambda_q * dipcorrB[j] * chargeB[i];
            }
            if (chargecorr[0] != 0 || chargecorr[1] != 0)
            {
                forces[i][ZZ] +=
                        (L1_q * chargecorr[0] * chargeA[i] + lambda_q * chargecorr[1]) * coords[i][ZZ];
            }
        }
    }

    /* Global corrections only on main process */
    if (MAIN(commrec) && thread == 0)
    {
        for (int q = 0; q < (bHaveChargePerturbed ? 2 : 1); q++)
        {
            /* Apply surface and charged surface dipole correction:
             * correction = dipole_coeff * ( (dipole)^2
             *              - qsum*sum_i q_i z_i^2 - qsum^2 * box_z^2 / 12 )
             */
            if (dipole_coeff != 0)
            {
                if (ewaldGeometry == EwaldGeometry::ThreeD)
                {
                    Vdipole[q] = dipole_coeff * iprod(mutot[q], mutot[q]);
                }
                else if (ewaldGeometry == EwaldGeometry::ThreeDC)
                {
                    Vdipole[q] = dipole_coeff * mutot[q][ZZ] * mutot[q][ZZ];

                    if (chargecorr[q] != 0)
                    {
                        /* Here we use a non thread-parallelized loop,
                         * because this is the only loop over atoms for
                         * energies and they need reduction (unlike forces).
                         * We could implement a reduction over threads,
                         * but this case is rarely used.
                         */
                        gmx::ArrayRef<const real> charge = (q == 0 ? chargeA : chargeB);
                        real                      sumQZ2 = 0;
                        for (int i = 0; i < numAtomsLocal; i++)
                        {
                            sumQZ2 += charge[i] * coords[i][ZZ] * coords[i][ZZ];
                        }
                        Vdipole[q] -= dipole_coeff * qsum[q]
                                      * (sumQZ2 + qsum[q] * box[ZZ][ZZ] * box[ZZ][ZZ] / 12);
                    }
                }
            }
        }
    }
    if (!bHaveChargePerturbed)
    {
        *Vcorr_q = Vdipole[0];
    }
    else
    {
        *Vcorr_q = L1_q * Vdipole[0] + lambda_q * Vdipole[1];
        *dvdlambda_q += Vdipole[1] - Vdipole[0];
    }
}
