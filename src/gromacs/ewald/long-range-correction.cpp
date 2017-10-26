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

#include "long-range-correction.h"

#include <cmath>

#include "gromacs/ewald/ewald-utils.h"
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

#include "pme-internal.h"

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
void ewald_LRcorrection(int numAtomsLocal,
                        t_commrec *cr,
                        int numThreads, int thread,
                        t_forcerec *fr,
                        const t_inputrec *ir,
                        real *chargeA, real *chargeB,
                        real *C6A, real *C6B,
                        real *sigmaA, real *sigmaB,
                        real *sigma3A, real *sigma3B,
                        gmx_bool bHaveChargeOrTypePerturbed,
                        gmx_bool calc_excl_corr,
                        t_blocka *excl, rvec x[],
                        matrix box, rvec mu_tot[],
                        int ewald_geometry, real epsilon_surface,
                        rvec *f, tensor vir_q, tensor vir_lj,
                        real *Vcorr_q, real *Vcorr_lj,
                        real lambda_q, real lambda_lj,
                        real *dvdlambda_q, real *dvdlambda_lj)
{
    int numAtomsToBeCorrected;
    if (calc_excl_corr)
    {
        /* We need to correct all exclusion pairs (cutoff-scheme = group) */
        numAtomsToBeCorrected = excl->nr;

        GMX_RELEASE_ASSERT(numAtomsToBeCorrected >= numAtomsLocal, "We might need to do self-corrections");
    }
    else
    {
        /* We need to correct only self interactions */
        numAtomsToBeCorrected = numAtomsLocal;
    }
    int         start =  (numAtomsToBeCorrected* thread     )/numThreads;
    int         end   =  (numAtomsToBeCorrected*(thread + 1))/numThreads;

    int         i, i1, i2, j, k, m, iv, jv, q;
    int        *AA;
    double      Vexcl_q, dvdl_excl_q, dvdl_excl_lj; /* Necessary for precision */
    double      Vexcl_lj;
    real        one_4pi_eps;
    real        v, vc, qiA, qiB, dr2, rinv;
    real        Vself_q[2], Vself_lj[2], Vdipole[2], rinv2, ewc_q = fr->ic->ewaldcoeff_q, ewcdr;
    real        ewc_lj = fr->ic->ewaldcoeff_lj, ewc_lj2 = ewc_lj * ewc_lj;
    real        c6Ai   = 0, c6Bi = 0, c6A = 0, c6B = 0, ewcdr2, ewcdr4, c6L = 0, rinv6;
    rvec        df, dx, mutot[2], dipcorrA, dipcorrB;
    tensor      dxdf_q = {{0}}, dxdf_lj = {{0}};
    real        L1_q, L1_lj, dipole_coeff, qqA, qqB, qqL, vr0_q, vr0_lj = 0;
    real        chargecorr[2] = { 0, 0 };
    gmx_bool    bMolPBC       = fr->bMolPBC;
    gmx_bool    bDoingLBRule  = (fr->ljpme_combination_rule == eljpmeLB);
    gmx_bool    bNeedLongRangeCorrection;

    GMX_ASSERT(ir, "Invalid inputrec pointer");
    matrix          scaledBox;
    EwaldBoxZScaler boxScaler(*ir);
    boxScaler.scaleBox(box, scaledBox);

    /* This routine can be made faster by using tables instead of analytical interactions
     * However, that requires a thorough verification that they are correct in all cases.
     */

    bool vdwPme   = EVDW_PME(fr->ic->vdwtype);

    one_4pi_eps   = ONE_4PI_EPS0/fr->ic->epsilon_r;
    vr0_q         = ewc_q*M_2_SQRTPI;
    if (vdwPme)
    {
        vr0_lj    = -gmx::power6(ewc_lj)/6.0;
    }

    AA           = excl->a;
    Vexcl_q      = 0;
    Vexcl_lj     = 0;
    dvdl_excl_q  = 0;
    dvdl_excl_lj = 0;
    Vdipole[0]   = 0;
    Vdipole[1]   = 0;
    L1_q         = 1.0-lambda_q;
    L1_lj        = 1.0-lambda_lj;
    /* Note that we have to transform back to gromacs units, since
     * mu_tot contains the dipole in debye units (for output).
     */
    for (i = 0; (i < DIM); i++)
    {
        mutot[0][i] = mu_tot[0][i]*DEBYE2ENM;
        mutot[1][i] = mu_tot[1][i]*DEBYE2ENM;
        dipcorrA[i] = 0;
        dipcorrB[i] = 0;
    }
    dipole_coeff = 0;

    real boxVolume = scaledBox[XX][XX]*scaledBox[YY][YY]*scaledBox[ZZ][ZZ];
    switch (ewald_geometry)
    {
        case eewg3D:
            if (epsilon_surface != 0)
            {
                dipole_coeff =
                    2*M_PI*ONE_4PI_EPS0/((2*epsilon_surface + fr->ic->epsilon_r)*boxVolume);
                for (i = 0; (i < DIM); i++)
                {
                    dipcorrA[i] = 2*dipole_coeff*mutot[0][i];
                    dipcorrB[i] = 2*dipole_coeff*mutot[1][i];
                }
            }
            break;
        case eewg3DC:
            dipole_coeff  = 2*M_PI*one_4pi_eps/boxVolume;
            dipcorrA[ZZ]  = 2*dipole_coeff*mutot[0][ZZ];
            dipcorrB[ZZ]  = 2*dipole_coeff*mutot[1][ZZ];
            for (int q = 0; q < (bHaveChargeOrTypePerturbed ? 2 : 1); q++)
            {
                /* Avoid charge corrections with near-zero net charge */
                if (fabs(fr->qsum[q]) > 1e-4)
                {
                    chargecorr[q] = 2*dipole_coeff*fr->qsum[q];
                }
            }
            break;
        default:
            gmx_incons("Unsupported Ewald geometry");
            break;
    }
    if (debug)
    {
        fprintf(debug, "dipcorr = %8.3f  %8.3f  %8.3f\n",
                dipcorrA[XX], dipcorrA[YY], dipcorrA[ZZ]);
        fprintf(debug, "mutot   = %8.3f  %8.3f  %8.3f\n",
                mutot[0][XX], mutot[0][YY], mutot[0][ZZ]);
    }
    bNeedLongRangeCorrection = (calc_excl_corr || dipole_coeff != 0);
    if (bNeedLongRangeCorrection && !bHaveChargeOrTypePerturbed)
    {
        for (i = start; (i < end); i++)
        {
            /* Initiate local variables (for this i-particle) to 0 */
            qiA = chargeA[i]*one_4pi_eps;
            if (vdwPme)
            {
                c6Ai = C6A[i];
                if (bDoingLBRule)
                {
                    c6Ai *= sigma3A[i];
                }
            }
            if (calc_excl_corr)
            {
                i1  = excl->index[i];
                i2  = excl->index[i+1];

                /* Loop over excluded neighbours */
                for (j = i1; (j < i2); j++)
                {
                    k = AA[j];
                    /*
                     * First we must test whether k <> i, and then,
                     * because the exclusions are all listed twice i->k
                     * and k->i we must select just one of the two.  As
                     * a minor optimization we only compute forces when
                     * the charges are non-zero.
                     */
                    if (k > i)
                    {
                        qqA = qiA*chargeA[k];
                        if (vdwPme)
                        {
                            c6A  = c6Ai * C6A[k];
                            if (bDoingLBRule)
                            {
                                c6A *= gmx::power6(0.5*(sigmaA[i]+sigmaA[k]))*sigma3A[k];
                            }
                        }
                        if (qqA != 0.0 || c6A != 0.0)
                        {
                            rvec_sub(x[i], x[k], dx);
                            if (bMolPBC)
                            {
                                /* Cheap pbc_dx, assume excluded pairs are at short distance. */
                                for (m = DIM-1; (m >= 0); m--)
                                {
                                    if (dx[m] > 0.5*box[m][m])
                                    {
                                        rvec_dec(dx, box[m]);
                                    }
                                    else if (dx[m] < -0.5*box[m][m])
                                    {
                                        rvec_inc(dx, box[m]);
                                    }
                                }
                            }
                            dr2 = norm2(dx);
                            /* Distance between two excluded particles
                             * may be zero in the case of shells
                             */
                            if (dr2 != 0)
                            {
                                rinv              = gmx::invsqrt(dr2);
                                rinv2             = rinv*rinv;
                                if (qqA != 0.0)
                                {
                                    real dr, fscal;

                                    dr       = 1.0/rinv;
                                    ewcdr    = ewc_q*dr;
                                    vc       = qqA*std::erf(ewcdr)*rinv;
                                    Vexcl_q += vc;
#if GMX_DOUBLE
                                    /* Relative accuracy at R_ERF_R_INACC of 3e-10 */
#define       R_ERF_R_INACC 0.006
#else
                                    /* Relative accuracy at R_ERF_R_INACC of 2e-5 */
#define       R_ERF_R_INACC 0.1
#endif
                                    /* fscal is the scalar force pre-multiplied by rinv,
                                     * to normalise the relative position vector dx */
                                    if (ewcdr > R_ERF_R_INACC)
                                    {
                                        fscal = rinv2*(vc - qqA*ewc_q*M_2_SQRTPI*exp(-ewcdr*ewcdr));
                                    }
                                    else
                                    {
                                        /* Use a fourth order series expansion for small ewcdr */
                                        fscal = ewc_q*ewc_q*qqA*vr0_q*(2.0/3.0 - 0.4*ewcdr*ewcdr);
                                    }

                                    /* The force vector is obtained by multiplication with
                                     * the relative position vector
                                     */
                                    svmul(fscal, dx, df);
                                    rvec_inc(f[k], df);
                                    rvec_dec(f[i], df);
                                    for (iv = 0; (iv < DIM); iv++)
                                    {
                                        for (jv = 0; (jv < DIM); jv++)
                                        {
                                            dxdf_q[iv][jv] += dx[iv]*df[jv];
                                        }
                                    }
                                }

                                if (c6A != 0.0)
                                {
                                    real fscal;

                                    rinv6     = rinv2*rinv2*rinv2;
                                    ewcdr2    = ewc_lj2*dr2;
                                    ewcdr4    = ewcdr2*ewcdr2;
                                    /* We get the excluded long-range contribution from -C6*(1-g(r))
                                     * g(r) is also defined in the manual under LJ-PME
                                     */
                                    vc        = -c6A*rinv6*(1.0 - exp(-ewcdr2)*(1 + ewcdr2 + 0.5*ewcdr4));
                                    Vexcl_lj += vc;
                                    /* The force is the derivative of the potential vc.
                                     * fscal is the scalar force pre-multiplied by rinv,
                                     * to normalise the relative position vector dx */
                                    fscal     = 6.0*vc*rinv2 + c6A*rinv6*exp(-ewcdr2)*ewc_lj2*ewcdr4;

                                    /* The force vector is obtained by multiplication with
                                     * the relative position vector
                                     */
                                    svmul(fscal, dx, df);
                                    rvec_inc(f[k], df);
                                    rvec_dec(f[i], df);
                                    for (iv = 0; (iv < DIM); iv++)
                                    {
                                        for (jv = 0; (jv < DIM); jv++)
                                        {
                                            dxdf_lj[iv][jv] += dx[iv]*df[jv];
                                        }
                                    }
                                }
                            }
                            else
                            {
                                Vexcl_q  += qqA*vr0_q;
                                Vexcl_lj += c6A*vr0_lj;
                            }
                        }
                    }
                }
            }
            /* Dipole correction on force */
            if (dipole_coeff != 0 && i < numAtomsLocal)
            {
                for (j = 0; (j < DIM); j++)
                {
                    f[i][j] -= dipcorrA[j]*chargeA[i];
                }
                if (chargecorr[0] != 0)
                {
                    f[i][ZZ] += chargecorr[0]*chargeA[i]*x[i][ZZ];
                }
            }
        }
    }
    else if (bNeedLongRangeCorrection)
    {
        for (i = start; (i < end); i++)
        {
            /* Initiate local variables (for this i-particle) to 0 */
            qiA = chargeA[i]*one_4pi_eps;
            qiB = chargeB[i]*one_4pi_eps;
            if (vdwPme)
            {
                c6Ai = C6A[i];
                c6Bi = C6B[i];
                if (bDoingLBRule)
                {
                    c6Ai *= sigma3A[i];
                    c6Bi *= sigma3B[i];
                }
            }
            if (calc_excl_corr)
            {
                i1  = excl->index[i];
                i2  = excl->index[i+1];

                /* Loop over excluded neighbours */
                for (j = i1; (j < i2); j++)
                {
                    k = AA[j];
                    if (k > i)
                    {
                        qqA = qiA*chargeA[k];
                        qqB = qiB*chargeB[k];
                        if (vdwPme)
                        {
                            c6A = c6Ai*C6A[k];
                            c6B = c6Bi*C6B[k];
                            if (bDoingLBRule)
                            {
                                c6A *= gmx::power6(0.5*(sigmaA[i]+sigmaA[k]))*sigma3A[k];
                                c6B *= gmx::power6(0.5*(sigmaB[i]+sigmaB[k]))*sigma3B[k];
                            }
                        }
                        if (qqA != 0.0 || qqB != 0.0 || c6A != 0.0 || c6B != 0.0)
                        {
                            real fscal;

                            qqL   = L1_q*qqA + lambda_q*qqB;
                            if (vdwPme)
                            {
                                c6L = L1_lj*c6A + lambda_lj*c6B;
                            }
                            rvec_sub(x[i], x[k], dx);
                            if (bMolPBC)
                            {
                                /* Cheap pbc_dx, assume excluded pairs are at short distance. */
                                for (m = DIM-1; (m >= 0); m--)
                                {
                                    if (dx[m] > 0.5*box[m][m])
                                    {
                                        rvec_dec(dx, box[m]);
                                    }
                                    else if (dx[m] < -0.5*box[m][m])
                                    {
                                        rvec_inc(dx, box[m]);
                                    }
                                }
                            }
                            dr2 = norm2(dx);
                            if (dr2 != 0)
                            {
                                rinv    = gmx::invsqrt(dr2);
                                rinv2   = rinv*rinv;
                                if (qqA != 0.0 || qqB != 0.0)
                                {
                                    real dr;

                                    dr           = 1.0/rinv;
                                    v            = std::erf(ewc_q*dr)*rinv;
                                    vc           = qqL*v;
                                    Vexcl_q     += vc;
                                    /* fscal is the scalar force pre-multiplied by rinv,
                                     * to normalise the relative position vector dx */
                                    fscal        = rinv2*(vc-qqL*ewc_q*M_2_SQRTPI*exp(-ewc_q*ewc_q*dr2));
                                    dvdl_excl_q += (qqB - qqA)*v;

                                    /* The force vector is obtained by multiplication with
                                     * the relative position vector
                                     */
                                    svmul(fscal, dx, df);
                                    rvec_inc(f[k], df);
                                    rvec_dec(f[i], df);
                                    for (iv = 0; (iv < DIM); iv++)
                                    {
                                        for (jv = 0; (jv < DIM); jv++)
                                        {
                                            dxdf_q[iv][jv] += dx[iv]*df[jv];
                                        }
                                    }
                                }

                                if ((c6A != 0.0 || c6B != 0.0) && vdwPme)
                                {
                                    rinv6         = rinv2*rinv2*rinv2;
                                    ewcdr2        = ewc_lj2*dr2;
                                    ewcdr4        = ewcdr2*ewcdr2;
                                    v             = -rinv6*(1.0 - exp(-ewcdr2)*(1 + ewcdr2 + 0.5*ewcdr4));
                                    vc            = c6L*v;
                                    Vexcl_lj     += vc;
                                    /* fscal is the scalar force pre-multiplied by rinv,
                                     * to normalise the relative position vector dx */
                                    fscal         = 6.0*vc*rinv2 + c6L*rinv6*exp(-ewcdr2)*ewc_lj2*ewcdr4;
                                    dvdl_excl_lj += (c6B - c6A)*v;

                                    /* The force vector is obtained by multiplication with
                                     * the relative position vector
                                     */
                                    svmul(fscal, dx, df);
                                    rvec_inc(f[k], df);
                                    rvec_dec(f[i], df);
                                    for (iv = 0; (iv < DIM); iv++)
                                    {
                                        for (jv = 0; (jv < DIM); jv++)
                                        {
                                            dxdf_lj[iv][jv] += dx[iv]*df[jv];
                                        }
                                    }
                                }
                            }
                            else
                            {
                                Vexcl_q      += qqL*vr0_q;
                                dvdl_excl_q  += (qqB - qqA)*vr0_q;
                                Vexcl_lj     += c6L*vr0_lj;
                                dvdl_excl_lj += (c6B - c6A)*vr0_lj;
                            }
                        }
                    }
                }
            }
            /* Dipole correction on force */
            if (dipole_coeff != 0 && i < numAtomsLocal)
            {
                for (j = 0; (j < DIM); j++)
                {
                    f[i][j] -= L1_q*dipcorrA[j]*chargeA[i]
                        + lambda_q*dipcorrB[j]*chargeB[i];
                }
                if (chargecorr[0] != 0 || chargecorr[1] != 0)
                {
                    f[i][ZZ] += (L1_q*chargecorr[0]*chargeA[i]
                                 + lambda_q*chargecorr[1])*x[i][ZZ];
                }
            }
        }
    }
    for (iv = 0; (iv < DIM); iv++)
    {
        for (jv = 0; (jv < DIM); jv++)
        {
            vir_q[iv][jv]  += 0.5*dxdf_q[iv][jv];
            vir_lj[iv][jv] += 0.5*dxdf_lj[iv][jv];
        }
    }

    Vself_q[0]  = 0;
    Vself_q[1]  = 0;
    Vself_lj[0] = 0;
    Vself_lj[1] = 0;

    /* Global corrections only on master process */
    if (MASTER(cr) && thread == 0)
    {
        for (q = 0; q < (bHaveChargeOrTypePerturbed ? 2 : 1); q++)
        {
            if (calc_excl_corr)
            {
                /* Self-energy correction */
                Vself_q[q] = ewc_q*one_4pi_eps*fr->q2sum[q]*M_1_SQRTPI;
                if (vdwPme)
                {
                    Vself_lj[q] =  fr->c6sum[q]*0.5*vr0_lj;
                }
            }

            /* Apply surface and charged surface dipole correction:
             * correction = dipole_coeff * ( (dipole)^2
             *              - qsum*sum_i q_i z_i^2 - qsum^2 * box_z^2 / 12 )
             */
            if (dipole_coeff != 0)
            {
                if (ewald_geometry == eewg3D)
                {
                    Vdipole[q] = dipole_coeff*iprod(mutot[q], mutot[q]);
                }
                else if (ewald_geometry == eewg3DC)
                {
                    Vdipole[q] = dipole_coeff*mutot[q][ZZ]*mutot[q][ZZ];

                    if (chargecorr[q] != 0)
                    {
                        /* Here we use a non thread-parallelized loop,
                         * because this is the only loop over atoms for
                         * energies and they need reduction (unlike forces).
                         * We could implement a reduction over threads,
                         * but this case is rarely used.
                         */
                        const real *qPtr   = (q == 0 ? chargeA : chargeB);
                        real        sumQZ2 = 0;
                        for (int i = 0; i < numAtomsLocal; i++)
                        {
                            sumQZ2 += qPtr[i]*x[i][ZZ]*x[i][ZZ];
                        }
                        Vdipole[q] -= dipole_coeff*fr->qsum[q]*(sumQZ2 + fr->qsum[q]*box[ZZ][ZZ]*box[ZZ][ZZ]/12);
                    }
                }
            }
        }
    }
    if (!bHaveChargeOrTypePerturbed)
    {
        *Vcorr_q = Vdipole[0] - Vself_q[0] - Vexcl_q;
        if (vdwPme)
        {
            *Vcorr_lj = -Vself_lj[0] - Vexcl_lj;
        }
    }
    else
    {
        *Vcorr_q = L1_q*(Vdipole[0] - Vself_q[0])
            + lambda_q*(Vdipole[1] - Vself_q[1])
            - Vexcl_q;
        *dvdlambda_q += Vdipole[1] - Vself_q[1]
            - (Vdipole[0] - Vself_q[0]) - dvdl_excl_q;
        if (vdwPme)
        {
            *Vcorr_lj      = -(L1_lj*Vself_lj[0] + lambda_lj*Vself_lj[1]) - Vexcl_lj;
            *dvdlambda_lj += -Vself_lj[1] + Vself_lj[0] - dvdl_excl_lj;
        }
    }

    if (debug)
    {
        fprintf(debug, "Long Range corrections for Ewald interactions:\n");
        fprintf(debug, "q2sum = %g, Vself_q=%g c6sum = %g, Vself_lj=%g\n",
                L1_q*fr->q2sum[0]+lambda_q*fr->q2sum[1], L1_q*Vself_q[0]+lambda_q*Vself_q[1], L1_lj*fr->c6sum[0]+lambda_lj*fr->c6sum[1], L1_lj*Vself_lj[0]+lambda_lj*Vself_lj[1]);
        fprintf(debug, "Electrostatic Long Range correction: Vexcl=%g\n", Vexcl_q);
        fprintf(debug, "Lennard-Jones Long Range correction: Vexcl=%g\n", Vexcl_lj);
        if (MASTER(cr) && thread == 0)
        {
            if (epsilon_surface > 0 || ewald_geometry == eewg3DC)
            {
                fprintf(debug, "Total dipole correction: Vdipole=%g\n",
                        L1_q*Vdipole[0]+lambda_q*Vdipole[1]);
            }
        }
    }
}
