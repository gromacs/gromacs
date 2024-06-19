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
 * \brief Defines SETTLE code.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "settle.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <filesystem>

#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/pbc_simd.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/*! \brief Initializes a projection matrix.
 *
 * \param[in]  invmO                  Reciprocal oxygen mass
 * \param[in]  invmH                  Reciprocal hydrogen mass
 * \param[in]  dOH                    Target O-H bond length
 * \param[in]  dHH                    Target H-H bond length
 * \param[out] inverseCouplingMatrix  Inverse bond coupling matrix for the projection version of SETTLE
 */
static void initializeProjectionMatrix(const real invmO,
                                       const real invmH,
                                       const real dOH,
                                       const real dHH,
                                       matrix     inverseCouplingMatrix)
{
    // We normalize the inverse masses with invmO for the matrix inversion.
    // so we can keep using masses of almost zero for frozen particles,
    // without running out of the float range in invertMatrix.
    double invmORelative = 1.0;
    double invmHRelative = invmH / static_cast<double>(invmO);
    double distanceRatio = dHH / static_cast<double>(dOH);

    /* Construct the constraint coupling matrix */
    matrix mat;
    mat[0][0] = invmORelative + invmHRelative;
    mat[0][1] = invmORelative * (1.0 - 0.5 * gmx::square(distanceRatio));
    mat[0][2] = invmHRelative * 0.5 * distanceRatio;
    mat[1][1] = mat[0][0];
    mat[1][2] = mat[0][2];
    mat[2][2] = invmHRelative + invmHRelative;
    mat[1][0] = mat[0][1];
    mat[2][0] = mat[0][2];
    mat[2][1] = mat[1][2];

    invertMatrix(mat, inverseCouplingMatrix);

    msmul(inverseCouplingMatrix, 1 / invmO, inverseCouplingMatrix);
}

SettleParameters
settleParameters(const real mO, const real mH, const real invmO, const real invmH, const real dOH, const real dHH)
{
    SettleParameters params;

    // We calculate parameters in double precision to minimize errors.
    // The velocity correction applied during SETTLE coordinate constraining
    // introduces a systematic error of approximately 1 bit per atom,
    // depending on what the compiler does with the code.
    double wohh;

    params.mO   = mO;
    params.mH   = mH;
    wohh        = mO + 2.0 * mH;
    params.wh   = mH / wohh;
    params.dOH  = dOH;
    params.dHH  = dHH;
    double rc   = dHH / 2.0;
    double ra   = 2.0 * mH * std::sqrt(dOH * dOH - rc * rc) / wohh;
    params.rb   = std::sqrt(dOH * dOH - rc * rc) - ra;
    params.rc   = rc;
    params.ra   = ra;
    params.irc2 = 1.0 / dHH;

    // For projection: inverse masses and coupling matrix inversion
    params.imO = invmO;
    params.imH = invmH;

    params.invdOH = 1.0 / dOH;
    params.invdHH = 1.0 / dHH;

    initializeProjectionMatrix(invmO, invmH, dOH, dHH, params.invmat);

    if (debug)
    {
        fprintf(debug, "wh =%g, rc = %g, ra = %g\n", params.wh, params.rc, params.ra);
        fprintf(debug, "rb = %g, irc2 = %g, dHH = %g, dOH = %g\n", params.rb, params.irc2, params.dHH, params.dOH);
    }

    return params;
}

SettleData::SettleData(const gmx_mtop_t& mtop) :
    useSimd_(getenv("GMX_DISABLE_SIMD_KERNELS") == nullptr)
{
    /* Check that we have only one settle type */
    int       settle_type = -1;
    const int nral1       = 1 + NRAL(F_SETTLE);
    for (const auto ilists : IListRange(mtop))
    {
        const InteractionList& ilist = ilists.list()[F_SETTLE];
        for (int i = 0; i < ilist.size(); i += nral1)
        {
            if (settle_type == -1)
            {
                settle_type = ilist.iatoms[i];
            }
            else if (ilist.iatoms[i] != settle_type)
            {
                gmx_fatal(FARGS,
                          "The [molecules] section of your topology specifies more than one block "
                          "of\n"
                          "a [moleculetype] with a [settles] block. Only one such is allowed.\n"
                          "If you are trying to partition your solvent into different *groups*\n"
                          "(e.g. for freezing, T-coupling, etc.), you are using the wrong "
                          "approach. Index\n"
                          "files specify groups. Otherwise, you may wish to change the least-used\n"
                          "block of molecules with SETTLE constraints into 3 normal constraints.");
            }
        }
    }
    GMX_RELEASE_ASSERT(settle_type >= 0, "settle_init called without settles");

    /* We will not initialize the normal SETTLE parameters here yet,
     * since the atom (inv)masses can depend on the integrator and
     * free-energy perturbation. We set mO=-1 to trigger later initialization.
     */
    parametersMassWeighted_.mO = -1;

    real dOH              = mtop.ffparams.iparams[settle_type].settle.doh;
    real dHH              = mtop.ffparams.iparams[settle_type].settle.dhh;
    parametersAllMasses1_ = settleParameters(1.0, 1.0, 1.0, 1.0, dOH, dHH);
}

void SettleData::setConstraints(const InteractionList&    il_settle,
                                const int                 numHomeAtoms,
                                gmx::ArrayRef<const real> masses,
                                gmx::ArrayRef<const real> inverseMasses)
{
#if GMX_SIMD_HAVE_REAL
    const int pack_size = GMX_SIMD_REAL_WIDTH;
#else
    const int pack_size = 1;
#endif

    const int nral1   = 1 + NRAL(F_SETTLE);
    int       nsettle = il_settle.size() / nral1;
    numSettles_       = nsettle;

    if (nsettle > 0)
    {
        ArrayRef<const int> iatoms = il_settle.iatoms;

        /* Here we initialize the normal SETTLE parameters */
        if (parametersMassWeighted_.mO < 0)
        {
            int firstO              = iatoms[1];
            int firstH              = iatoms[2];
            parametersMassWeighted_ = settleParameters(masses[firstO],
                                                       masses[firstH],
                                                       inverseMasses[firstO],
                                                       inverseMasses[firstH],
                                                       parametersAllMasses1_.dOH,
                                                       parametersAllMasses1_.dHH);
        }

        const int paddedSize = divideRoundUp(nsettle, pack_size) * pack_size;
        ow1_.resize(paddedSize);
        hw2_.resize(paddedSize);
        hw3_.resize(paddedSize);
        virfac_.resize(paddedSize);

        for (int i = 0; i < nsettle; i++)
        {
            ow1_[i] = iatoms[i * nral1 + 1];
            hw2_[i] = iatoms[i * nral1 + 2];
            hw3_[i] = iatoms[i * nral1 + 3];
            /* We should avoid double counting of virial contributions for
             * SETTLEs that appear in multiple DD domains, so we only count
             * the contribution on the home range of the oxygen atom.
             */
            virfac_[i] = (iatoms[i * nral1 + 1] < numHomeAtoms ? 1 : 0);
        }

        /* Pad the index array to the full SIMD width with copies from
         * the last normal entry, but with no virial contribution.
         */
        for (int i = nsettle; i < paddedSize; i++)
        {
            ow1_[i]    = ow1_[nsettle - 1];
            hw2_[i]    = hw2_[nsettle - 1];
            hw3_[i]    = hw3_[nsettle - 1];
            virfac_[i] = 0;
        }
    }
}

void settle_proj(const SettleData&    settled,
                 ConstraintVariable   econq,
                 int                  nsettle,
                 const t_iatom        iatoms[],
                 const t_pbc*         pbc,
                 ArrayRef<const RVec> x,
                 ArrayRef<RVec>       der,
                 ArrayRef<RVec>       derp,
                 int                  calcvir_atom_end,
                 tensor               vir_r_m_dder)
{
    /* Settle for projection out constraint components
     * of derivatives of the coordinates.
     * Berk Hess 2008-1-10
     */

    const SettleParameters* p;
    real                    imO, imH, dOH, dHH, invdOH, invdHH;
    matrix                  invmat;
    int                     i, m, m2, ow1, hw2, hw3;
    rvec                    roh2, roh3, rhh, dc, fc;

    calcvir_atom_end *= DIM;

    if (econq == ConstraintVariable::Force)
    {
        p = &settled.parametersAllMasses1();
    }
    else
    {
        p = &settled.parametersMassWeighted();
    }
    imO = p->imO;
    imH = p->imH;
    copy_mat(p->invmat, invmat);
    dOH    = p->dOH;
    dHH    = p->dHH;
    invdOH = p->invdOH;
    invdHH = p->invdHH;

    const int nral1 = 1 + NRAL(F_SETTLE);

    for (i = 0; i < nsettle; i++)
    {
        ow1 = iatoms[i * nral1 + 1];
        hw2 = iatoms[i * nral1 + 2];
        hw3 = iatoms[i * nral1 + 3];

        if (pbc == nullptr)
        {
            rvec_sub(x[ow1], x[hw2], roh2);
            rvec_sub(x[ow1], x[hw3], roh3);
            rvec_sub(x[hw2], x[hw3], rhh);
        }
        else
        {
            pbc_dx_aiuc(pbc, x[ow1], x[hw2], roh2);
            pbc_dx_aiuc(pbc, x[ow1], x[hw3], roh3);
            pbc_dx_aiuc(pbc, x[hw2], x[hw3], rhh);
        }
        svmul(invdOH, roh2, roh2);
        svmul(invdOH, roh3, roh3);
        svmul(invdHH, rhh, rhh);
        /* 18 flops */

        /* Determine the projections of der on the bonds */
        clear_rvec(dc);
        for (m = 0; m < DIM; m++)
        {
            dc[0] += (der[ow1][m] - der[hw2][m]) * roh2[m];
            dc[1] += (der[ow1][m] - der[hw3][m]) * roh3[m];
            dc[2] += (der[hw2][m] - der[hw3][m]) * rhh[m];
        }
        /* 27 flops */

        /* Determine the correction for the three bonds */
        mvmul(invmat, dc, fc);
        /* 15 flops */

        /* Subtract the corrections from derp */
        for (m = 0; m < DIM; m++)
        {
            derp[ow1][m] -= imO * (fc[0] * roh2[m] + fc[1] * roh3[m]);
            derp[hw2][m] -= imH * (-fc[0] * roh2[m] + fc[2] * rhh[m]);
            derp[hw3][m] -= imH * (-fc[1] * roh3[m] - fc[2] * rhh[m]);
        }

        /* 45 flops */

        if (ow1 < calcvir_atom_end)
        {
            /* Determining r \dot m der is easy,
             * since fc contains the mass weighted corrections for der.
             */

            for (m = 0; m < DIM; m++)
            {
                for (m2 = 0; m2 < DIM; m2++)
                {
                    vir_r_m_dder[m][m2] += dOH * roh2[m] * roh2[m2] * fc[0]
                                           + dOH * roh3[m] * roh3[m2] * fc[1]
                                           + dHH * rhh[m] * rhh[m2] * fc[2];
                }
            }
        }
    }
}


/*! \brief The actual settle code, templated for real/SimdReal and for optimization */
template<typename T, typename TypeBool, int packSize, typename TypePbc, bool bCorrectVelocity, bool bCalcVirial>
static void settleTemplate(const SettleData&  settled,
                           int                settleStart,
                           int                settleEnd,
                           const TypePbc      pbc,
                           const real*        x,
                           real*              xprime,
                           real               invdt,
                           real* gmx_restrict v,
                           tensor             vir_r_m_dr,
                           bool*              bErrorHasOccurred)
{
    /* ******************************************************************* */
    /*                                                                  ** */
    /*    Original code by Shuichi Miyamoto, last update Oct. 1, 1992   ** */
    /*                                                                  ** */
    /*    Algorithm changes by Berk Hess:                               ** */
    /*    2004-07-15 Convert COM to double precision to avoid drift     ** */
    /*    2006-10-16 Changed velocity update to use differences         ** */
    /*    2012-09-24 Use oxygen as reference instead of COM             ** */
    /*    2016-02    Complete rewrite of the code for SIMD              ** */
    /*    2020-06    Completely remove use of COM to minimize drift     ** */
    /*                                                                  ** */
    /*    Reference for the SETTLE algorithm                            ** */
    /*           S. Miyamoto et al., J. Comp. Chem., 13, 952 (1992).    ** */
    /*                                                                  ** */
    /* ******************************************************************* */

    assert(settleStart % packSize == 0);
    assert(settleEnd % packSize == 0);

    TypeBool bError = TypeBool(false);

    const SettleParameters* p    = &settled.parametersMassWeighted();
    T                       wh   = T(p->wh);
    T                       rc   = T(p->rc);
    T                       ra   = T(p->ra);
    T                       rb   = T(p->rb);
    T                       irc2 = T(p->irc2);
    T                       mO   = T(p->mO);
    T                       mH   = T(p->mH);

    T almost_zero = T(1e-12);

    T sum_r_m_dr[DIM][DIM];

    if (bCalcVirial)
    {
        for (int d2 = 0; d2 < DIM; d2++)
        {
            for (int d = 0; d < DIM; d++)
            {
                sum_r_m_dr[d2][d] = T(0);
            }
        }
    }

    for (int i = settleStart; i < settleEnd; i += packSize)
    {
        /* Here we pad up to packSize with copies from the last valid entry.
         * This gives correct results, since we store (not increment) all
         * output, so we store the same output multiple times.
         */
        const int* ow1 = settled.ow1() + i;
        const int* hw2 = settled.hw2() + i;
        const int* hw3 = settled.hw3() + i;

        T x_ow1[DIM], x_hw2[DIM], x_hw3[DIM];

        gatherLoadUTranspose<3>(x, ow1, &x_ow1[XX], &x_ow1[YY], &x_ow1[ZZ]);
        gatherLoadUTranspose<3>(x, hw2, &x_hw2[XX], &x_hw2[YY], &x_hw2[ZZ]);
        gatherLoadUTranspose<3>(x, hw3, &x_hw3[XX], &x_hw3[YY], &x_hw3[ZZ]);

        T xprime_ow1[DIM], xprime_hw2[DIM], xprime_hw3[DIM];

        gatherLoadUTranspose<3>(xprime, ow1, &xprime_ow1[XX], &xprime_ow1[YY], &xprime_ow1[ZZ]);
        gatherLoadUTranspose<3>(xprime, hw2, &xprime_hw2[XX], &xprime_hw2[YY], &xprime_hw2[ZZ]);
        gatherLoadUTranspose<3>(xprime, hw3, &xprime_hw3[XX], &xprime_hw3[YY], &xprime_hw3[ZZ]);

        T dist21[DIM], dist31[DIM];
        T doh2[DIM], doh3[DIM];

        pbc_dx_aiuc(pbc, x_hw2, x_ow1, dist21);

        pbc_dx_aiuc(pbc, x_hw3, x_ow1, dist31);

        pbc_dx_aiuc(pbc, xprime_hw2, xprime_ow1, doh2);

        pbc_dx_aiuc(pbc, xprime_hw3, xprime_ow1, doh3);
        /* 4 * 18 flops (would be 4 * 3 without PBC) */

        /* Note that we completely avoid computing the center of mass and
         * only use distances. This minimizes energy drift and also makes
         * the computation slightly cheaper.
         * Straightforward computation of the COM, as in the original algorithm,
         * makes SETTLE the largest source of energy drift for simulations of water,
         * as then the oxygen coordinate is multiplied by 0.89 at every step,
         * which can then transfer a systematic rounding to the oxygen velocity.
         * For some time we computed the COM using offsets from the oxygen, this
         * significantly reduces the energy drift, but not using the COM at all,
         * as we do now, is optimal.
         */
        T a1[DIM];
        for (int d = 0; d < DIM; d++)
        {
            a1[d] = -(doh2[d] + doh3[d]) * wh;
        }
        T b1[DIM];
        for (int d = 0; d < DIM; d++)
        {
            b1[d] = doh2[d] + a1[d];
        }
        T c1[DIM];
        for (int d = 0; d < DIM; d++)
        {
            c1[d] = doh3[d] + a1[d];
        }
        /* 12 flops */

        T xakszd = dist21[YY] * dist31[ZZ] - dist21[ZZ] * dist31[YY];
        T yakszd = dist21[ZZ] * dist31[XX] - dist21[XX] * dist31[ZZ];
        T zakszd = dist21[XX] * dist31[YY] - dist21[YY] * dist31[XX];
        T xaksxd = a1[YY] * zakszd - a1[ZZ] * yakszd;
        T yaksxd = a1[ZZ] * xakszd - a1[XX] * zakszd;
        T zaksxd = a1[XX] * yakszd - a1[YY] * xakszd;
        T xaksyd = yakszd * zaksxd - zakszd * yaksxd;
        T yaksyd = zakszd * xaksxd - xakszd * zaksxd;
        T zaksyd = xakszd * yaksxd - yakszd * xaksxd;
        /* 27 flops */

        T axlng = gmx::invsqrt(xaksxd * xaksxd + yaksxd * yaksxd + zaksxd * zaksxd);
        T aylng = gmx::invsqrt(xaksyd * xaksyd + yaksyd * yaksyd + zaksyd * zaksyd);
        T azlng = gmx::invsqrt(xakszd * xakszd + yakszd * yakszd + zakszd * zakszd);

        T trns1[DIM], trns2[DIM], trns3[DIM];

        trns1[XX] = xaksxd * axlng;
        trns2[XX] = yaksxd * axlng;
        trns3[XX] = zaksxd * axlng;
        trns1[YY] = xaksyd * aylng;
        trns2[YY] = yaksyd * aylng;
        trns3[YY] = zaksyd * aylng;
        trns1[ZZ] = xakszd * azlng;
        trns2[ZZ] = yakszd * azlng;
        trns3[ZZ] = zakszd * azlng;
        /* 24 flops */

        T b0d[2], c0d[2];

        for (int d = 0; d < 2; d++)
        {
            b0d[d] = trns1[d] * dist21[XX] + trns2[d] * dist21[YY] + trns3[d] * dist21[ZZ];
            c0d[d] = trns1[d] * dist31[XX] + trns2[d] * dist31[YY] + trns3[d] * dist31[ZZ];
        }

        T a1d_z, b1d[DIM], c1d[DIM];

        a1d_z = trns1[ZZ] * a1[XX] + trns2[ZZ] * a1[YY] + trns3[ZZ] * a1[ZZ];
        for (int d = 0; d < DIM; d++)
        {
            b1d[d] = trns1[d] * b1[XX] + trns2[d] * b1[YY] + trns3[d] * b1[ZZ];
            c1d[d] = trns1[d] * c1[XX] + trns2[d] * c1[YY] + trns3[d] * c1[ZZ];
        }
        /* 65 flops */

        T tmp, tmp2;

        T sinphi = a1d_z * gmx::invsqrt(ra * ra);
        tmp2     = 1.0 - sinphi * sinphi;

        /* If tmp2 gets close to or beyond zero we have severly distorted
         * water molecules and we should terminate the simulation.
         * Below we take the max with almost_zero to continue the loop.
         */
        bError = bError || (tmp2 <= almost_zero);

        tmp2     = max(tmp2, almost_zero);
        tmp      = gmx::invsqrt(tmp2);
        T cosphi = tmp2 * tmp;
        T sinpsi = (b1d[ZZ] - c1d[ZZ]) * irc2 * tmp;
        tmp2     = 1.0 - sinpsi * sinpsi;

        T cospsi = tmp2 * gmx::invsqrt(tmp2);
        /* 46 flops */

        T a2d_y = ra * cosphi;
        T b2d_x = -rc * cospsi;
        T t1    = -rb * cosphi;
        T t2    = rc * sinpsi * sinphi;
        T b2d_y = t1 - t2;
        T c2d_y = t1 + t2;
        /* 7 flops */

        /*     --- Step3  al,be,ga            --- */
        T alpha  = b2d_x * (b0d[XX] - c0d[XX]) + b0d[YY] * b2d_y + c0d[YY] * c2d_y;
        T beta   = b2d_x * (c0d[YY] - b0d[YY]) + b0d[XX] * b2d_y + c0d[XX] * c2d_y;
        T gamma  = b0d[XX] * b1d[YY] - b1d[XX] * b0d[YY] + c0d[XX] * c1d[YY] - c1d[XX] * c0d[YY];
        T al2be2 = alpha * alpha + beta * beta;
        tmp2     = (al2be2 - gamma * gamma);
        T sinthe = (alpha * gamma - beta * tmp2 * gmx::invsqrt(tmp2)) * gmx::invsqrt(al2be2 * al2be2);
        /* 47 flops */

        /*  --- Step4  A3' --- */
        tmp2     = 1.0 - sinthe * sinthe;
        T costhe = tmp2 * gmx::invsqrt(tmp2);

        T a3d[DIM], b3d[DIM], c3d[DIM];

        a3d[XX] = -a2d_y * sinthe;
        a3d[YY] = a2d_y * costhe;
        a3d[ZZ] = a1d_z;
        b3d[XX] = b2d_x * costhe - b2d_y * sinthe;
        b3d[YY] = b2d_x * sinthe + b2d_y * costhe;
        b3d[ZZ] = b1d[ZZ];
        c3d[XX] = -b2d_x * costhe - c2d_y * sinthe;
        c3d[YY] = -b2d_x * sinthe + c2d_y * costhe;
        c3d[ZZ] = c1d[ZZ];
        /* 26 flops */

        /*    --- Step5  A3 --- */
        T a3[DIM], b3[DIM], c3[DIM];

        a3[XX] = trns1[XX] * a3d[XX] + trns1[YY] * a3d[YY] + trns1[ZZ] * a3d[ZZ];
        a3[YY] = trns2[XX] * a3d[XX] + trns2[YY] * a3d[YY] + trns2[ZZ] * a3d[ZZ];
        a3[ZZ] = trns3[XX] * a3d[XX] + trns3[YY] * a3d[YY] + trns3[ZZ] * a3d[ZZ];
        b3[XX] = trns1[XX] * b3d[XX] + trns1[YY] * b3d[YY] + trns1[ZZ] * b3d[ZZ];
        b3[YY] = trns2[XX] * b3d[XX] + trns2[YY] * b3d[YY] + trns2[ZZ] * b3d[ZZ];
        b3[ZZ] = trns3[XX] * b3d[XX] + trns3[YY] * b3d[YY] + trns3[ZZ] * b3d[ZZ];
        c3[XX] = trns1[XX] * c3d[XX] + trns1[YY] * c3d[YY] + trns1[ZZ] * c3d[ZZ];
        c3[YY] = trns2[XX] * c3d[XX] + trns2[YY] * c3d[YY] + trns2[ZZ] * c3d[ZZ];
        c3[ZZ] = trns3[XX] * c3d[XX] + trns3[YY] * c3d[YY] + trns3[ZZ] * c3d[ZZ];
        /* 45 flops */

        /* Compute and store the corrected new coordinate */
        T dxOw1[DIM];
        for (int d = 0; d < DIM; d++)
        {
            dxOw1[d]      = a3[d] - a1[d];
            xprime_ow1[d] = xprime_ow1[d] + dxOw1[d];
        }
        T dxHw2[DIM];
        for (int d = 0; d < DIM; d++)
        {
            dxHw2[d]      = b3[d] - b1[d];
            xprime_hw2[d] = xprime_hw2[d] + dxHw2[d];
        }
        T dxHw3[DIM];
        for (int d = 0; d < DIM; d++)
        {
            dxHw3[d]      = c3[d] - c1[d];
            xprime_hw3[d] = xprime_hw3[d] + dxHw3[d];
        }
        /* 9 + 9 flops */

        transposeScatterStoreU<3>(xprime, ow1, xprime_ow1[XX], xprime_ow1[YY], xprime_ow1[ZZ]);
        transposeScatterStoreU<3>(xprime, hw2, xprime_hw2[XX], xprime_hw2[YY], xprime_hw2[ZZ]);
        transposeScatterStoreU<3>(xprime, hw3, xprime_hw3[XX], xprime_hw3[YY], xprime_hw3[ZZ]);

        if (bCorrectVelocity)
        {
            T v_ow1[DIM], v_hw2[DIM], v_hw3[DIM];

            gatherLoadUTranspose<3>(v, ow1, &v_ow1[XX], &v_ow1[YY], &v_ow1[ZZ]);
            gatherLoadUTranspose<3>(v, hw2, &v_hw2[XX], &v_hw2[YY], &v_hw2[ZZ]);
            gatherLoadUTranspose<3>(v, hw3, &v_hw3[XX], &v_hw3[YY], &v_hw3[ZZ]);

            /* Add the position correction divided by dt to the velocity */
            for (int d = 0; d < DIM; d++)
            {
                v_ow1[d] = gmx::fma(dxOw1[d], invdt, v_ow1[d]);
            }
            for (int d = 0; d < DIM; d++)
            {
                v_hw2[d] = gmx::fma(dxHw2[d], invdt, v_hw2[d]);
            }
            for (int d = 0; d < DIM; d++)
            {
                v_hw3[d] = gmx::fma(dxHw3[d], invdt, v_hw3[d]);
            }
            /* 3*6 flops */

            transposeScatterStoreU<3>(v, ow1, v_ow1[XX], v_ow1[YY], v_ow1[ZZ]);
            transposeScatterStoreU<3>(v, hw2, v_hw2[XX], v_hw2[YY], v_hw2[ZZ]);
            transposeScatterStoreU<3>(v, hw3, v_hw3[XX], v_hw3[YY], v_hw3[ZZ]);
        }

        if (bCalcVirial)
        {
            /* Filter out the non-local settles */
            T filter = load<T>(settled.virfac() + i);
            T mOf    = filter * mO;
            T mHf    = filter * mH;

            T mdo[DIM], mdb[DIM], mdc[DIM];

            for (int d = 0; d < DIM; d++)
            {
                mdb[d] = mHf * dxHw2[d];
                mdc[d] = mHf * dxHw3[d];
                mdo[d] = mOf * dxOw1[d] + mdb[d] + mdc[d];
            }

            for (int d2 = 0; d2 < DIM; d2++)
            {
                for (int d = 0; d < DIM; d++)
                {
                    sum_r_m_dr[d2][d] =
                            sum_r_m_dr[d2][d]
                            - (x_ow1[d2] * mdo[d] + dist21[d2] * mdb[d] + dist31[d2] * mdc[d]);
                }
            }
            /* 71 flops */
        }
    }

    if (bCalcVirial)
    {
        for (int d2 = 0; d2 < DIM; d2++)
        {
            for (int d = 0; d < DIM; d++)
            {
                vir_r_m_dr[d2][d] += reduce(sum_r_m_dr[d2][d]);
            }
        }
    }

    *bErrorHasOccurred = anyTrue(bError);
}

/*! \brief Wrapper template function that divides the settles over threads
 * and instantiates the core template with instantiated booleans.
 */
template<typename T, typename TypeBool, int packSize, typename TypePbc>
static void settleTemplateWrapper(const SettleData& settled,
                                  int               nthread,
                                  int               thread,
                                  TypePbc           pbc,
                                  const real        x[],
                                  real              xprime[],
                                  real              invdt,
                                  real*             v,
                                  bool              bCalcVirial,
                                  tensor            vir_r_m_dr,
                                  bool*             bErrorHasOccurred)
{
    /* We need to assign settles to threads in groups of pack_size */
    int numSettlePacks = divideRoundUp(settled.numSettles(), packSize);
    /* Round the end value up to give thread 0 more work */
    int settleStart = divideRoundUp(numSettlePacks * thread, nthread) * packSize;
    int settleEnd   = divideRoundUp(numSettlePacks * (thread + 1), nthread) * packSize;

    if (v != nullptr)
    {
        if (!bCalcVirial)
        {
            settleTemplate<T, TypeBool, packSize, TypePbc, true, false>(
                    settled, settleStart, settleEnd, pbc, x, xprime, invdt, v, nullptr, bErrorHasOccurred);
        }
        else
        {
            settleTemplate<T, TypeBool, packSize, TypePbc, true, true>(
                    settled, settleStart, settleEnd, pbc, x, xprime, invdt, v, vir_r_m_dr, bErrorHasOccurred);
        }
    }
    else
    {
        if (!bCalcVirial)
        {
            settleTemplate<T, TypeBool, packSize, TypePbc, false, false>(
                    settled, settleStart, settleEnd, pbc, x, xprime, invdt, v, nullptr, bErrorHasOccurred);
        }
        else
        {
            settleTemplate<T, TypeBool, packSize, TypePbc, false, true>(
                    settled, settleStart, settleEnd, pbc, x, xprime, invdt, v, vir_r_m_dr, bErrorHasOccurred);
        }
    }
}

void csettle(const SettleData&               settled,
             int                             nthread,
             int                             thread,
             const t_pbc*                    pbc,
             ArrayRefWithPadding<const RVec> x,
             ArrayRefWithPadding<RVec>       xprime,
             real                            invdt,
             ArrayRefWithPadding<RVec>       v,
             bool                            bCalcVirial,
             tensor                          vir_r_m_dr,
             bool*                           bErrorHasOccurred)
{
    const real* xPtr      = as_rvec_array(x.paddedArrayRef().data())[0];
    real*       xprimePtr = as_rvec_array(xprime.paddedArrayRef().data())[0];
    real*       vPtr      = as_rvec_array(v.paddedArrayRef().data())[0];

#if GMX_SIMD_HAVE_REAL
    if (settled.useSimd())
    {
        /* Convert the pbc struct for SIMD */
        alignas(GMX_SIMD_ALIGNMENT) real pbcSimd[9 * GMX_SIMD_REAL_WIDTH];
        set_pbc_simd(pbc, pbcSimd);

        settleTemplateWrapper<SimdReal, SimdBool, GMX_SIMD_REAL_WIDTH, const real*>(
                settled, nthread, thread, pbcSimd, xPtr, xprimePtr, invdt, vPtr, bCalcVirial, vir_r_m_dr, bErrorHasOccurred);
    }
    else
#endif
    {
        /* This construct is needed because pbc_dx_aiuc doesn't accept pbc=NULL */
        t_pbc        pbcNo;
        const t_pbc* pbcNonNull;

        if (pbc != nullptr)
        {
            pbcNonNull = pbc;
        }
        else
        {
            set_pbc(&pbcNo, PbcType::No, nullptr);
            pbcNonNull = &pbcNo;
        }

        settleTemplateWrapper<real, bool, 1, const t_pbc*>(
                settled, nthread, thread, pbcNonNull, &xPtr[0], &xprimePtr[0], invdt, &vPtr[0], bCalcVirial, vir_r_m_dr, bErrorHasOccurred);
    }
}

} // namespace gmx
