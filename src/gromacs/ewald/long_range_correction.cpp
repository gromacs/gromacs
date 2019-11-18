/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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

#include "long_range_correction.h"

#include <cmath>

#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

#include "pme_internal.h"

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
void ewald_LRcorrection(const int         numAtomsLocal,
                        const t_commrec*  cr,
                        int               numThreads,
                        int               thread,
                        const t_forcerec& fr,
                        const t_inputrec& ir,
                        const real*       chargeA,
                        const real*       chargeB,
                        gmx_bool          bHaveChargePerturbed,
                        const rvec        x[],
                        const matrix      box,
                        const rvec        mu_tot[],
                        rvec*             f,
                        real*             Vcorr_q,
                        real              lambda_q,
                        real*             dvdlambda_q)
{
    /* We need to correct only self interactions */
    const int start = (numAtomsLocal * thread) / numThreads;
    const int end   = (numAtomsLocal * (thread + 1)) / numThreads;

    int    i, j, q;
    double Vexcl_q, dvdl_excl_q; /* Necessary for precision */
    real   one_4pi_eps;
    real   Vself_q[2], Vdipole[2];
    rvec   mutot[2], dipcorrA, dipcorrB;
    real   L1_q, dipole_coeff;
    real   chargecorr[2] = { 0, 0 };

    /* Scale the Ewald unit cell when dimension z is not periodic */
    matrix          scaledBox;
    EwaldBoxZScaler boxScaler(ir);
    boxScaler.scaleBox(box, scaledBox);

    one_4pi_eps = ONE_4PI_EPS0 / fr.ic->epsilon_r;
    Vexcl_q     = 0;
    dvdl_excl_q = 0;
    Vdipole[0]  = 0;
    Vdipole[1]  = 0;
    L1_q        = 1.0 - lambda_q;
    /* Note that we have to transform back to gromacs units, since
     * mu_tot contains the dipole in debye units (for output).
     */
    for (i = 0; (i < DIM); i++)
    {
        mutot[0][i] = mu_tot[0][i] * DEBYE2ENM;
        mutot[1][i] = mu_tot[1][i] * DEBYE2ENM;
        dipcorrA[i] = 0;
        dipcorrB[i] = 0;
    }
    dipole_coeff = 0;

    real boxVolume = scaledBox[XX][XX] * scaledBox[YY][YY] * scaledBox[ZZ][ZZ];
    switch (ir.ewald_geometry)
    {
        case eewg3D:
            if (ir.epsilon_surface != 0)
            {
                dipole_coeff = 2 * M_PI * ONE_4PI_EPS0
                               / ((2 * ir.epsilon_surface + fr.ic->epsilon_r) * boxVolume);
                for (i = 0; (i < DIM); i++)
                {
                    dipcorrA[i] = 2 * dipole_coeff * mutot[0][i];
                    dipcorrB[i] = 2 * dipole_coeff * mutot[1][i];
                }
            }
            break;
        case eewg3DC:
            dipole_coeff = 2 * M_PI * one_4pi_eps / boxVolume;
            dipcorrA[ZZ] = 2 * dipole_coeff * mutot[0][ZZ];
            dipcorrB[ZZ] = 2 * dipole_coeff * mutot[1][ZZ];
            for (int q = 0; q < (bHaveChargePerturbed ? 2 : 1); q++)
            {
                /* Avoid charge corrections with near-zero net charge */
                if (fabs(fr.qsum[q]) > 1e-4)
                {
                    chargecorr[q] = 2 * dipole_coeff * fr.qsum[q];
                }
            }
            break;
        default: gmx_incons("Unsupported Ewald geometry");
    }
    if (debug)
    {
        fprintf(debug, "dipcorr = %8.3f  %8.3f  %8.3f\n", dipcorrA[XX], dipcorrA[YY], dipcorrA[ZZ]);
        fprintf(debug, "mutot   = %8.3f  %8.3f  %8.3f\n", mutot[0][XX], mutot[0][YY], mutot[0][ZZ]);
    }
    const bool bNeedLongRangeCorrection = (dipole_coeff != 0);
    if (bNeedLongRangeCorrection && !bHaveChargePerturbed)
    {
        for (i = start; (i < end); i++)
        {
            for (j = 0; (j < DIM); j++)
            {
                f[i][j] -= dipcorrA[j] * chargeA[i];
            }
            if (chargecorr[0] != 0)
            {
                f[i][ZZ] += chargecorr[0] * chargeA[i] * x[i][ZZ];
            }
        }
    }
    else if (bNeedLongRangeCorrection)
    {
        for (i = start; (i < end); i++)
        {
            for (j = 0; (j < DIM); j++)
            {
                f[i][j] -= L1_q * dipcorrA[j] * chargeA[i] + lambda_q * dipcorrB[j] * chargeB[i];
            }
            if (chargecorr[0] != 0 || chargecorr[1] != 0)
            {
                f[i][ZZ] += (L1_q * chargecorr[0] * chargeA[i] + lambda_q * chargecorr[1]) * x[i][ZZ];
            }
        }
    }

    Vself_q[0] = 0;
    Vself_q[1] = 0;

    /* Global corrections only on master process */
    if (MASTER(cr) && thread == 0)
    {
        for (q = 0; q < (bHaveChargePerturbed ? 2 : 1); q++)
        {
            /* Apply surface and charged surface dipole correction:
             * correction = dipole_coeff * ( (dipole)^2
             *              - qsum*sum_i q_i z_i^2 - qsum^2 * box_z^2 / 12 )
             */
            if (dipole_coeff != 0)
            {
                if (ir.ewald_geometry == eewg3D)
                {
                    Vdipole[q] = dipole_coeff * iprod(mutot[q], mutot[q]);
                }
                else if (ir.ewald_geometry == eewg3DC)
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
                        const real* qPtr   = (q == 0 ? chargeA : chargeB);
                        real        sumQZ2 = 0;
                        for (int i = 0; i < numAtomsLocal; i++)
                        {
                            sumQZ2 += qPtr[i] * x[i][ZZ] * x[i][ZZ];
                        }
                        Vdipole[q] -= dipole_coeff * fr.qsum[q]
                                      * (sumQZ2 + fr.qsum[q] * box[ZZ][ZZ] * box[ZZ][ZZ] / 12);
                    }
                }
            }
        }
    }
    if (!bHaveChargePerturbed)
    {
        *Vcorr_q = Vdipole[0] - Vself_q[0] - Vexcl_q;
    }
    else
    {
        *Vcorr_q = L1_q * (Vdipole[0] - Vself_q[0]) + lambda_q * (Vdipole[1] - Vself_q[1]) - Vexcl_q;
        *dvdlambda_q += Vdipole[1] - Vself_q[1] - (Vdipole[0] - Vself_q[0]) - dvdl_excl_q;
    }

    if (debug)
    {
        fprintf(debug, "Long Range corrections for Ewald interactions:\n");
        fprintf(debug, "q2sum = %g, Vself_q=%g\n", L1_q * fr.q2sum[0] + lambda_q * fr.q2sum[1],
                L1_q * Vself_q[0] + lambda_q * Vself_q[1]);
        fprintf(debug, "Electrostatic Long Range correction: Vexcl=%g\n", Vexcl_q);
        if (MASTER(cr) && thread == 0)
        {
            if (ir.epsilon_surface > 0 || ir.ewald_geometry == eewg3DC)
            {
                fprintf(debug, "Total dipole correction: Vdipole=%g\n",
                        L1_q * Vdipole[0] + lambda_q * Vdipole[1]);
            }
        }
    }
}
