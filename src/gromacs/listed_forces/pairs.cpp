/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
 * \brief This file defines functions for "pair" interactions
 * (i.e. listed non-bonded interactions, e.g. 1-4 interactions)
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_listed_forces
 */
#include "gmxpre.h"

#include "pairs.h"

#include <cmath>
#include <cstdint>
#include <cstdio>

#include <filesystem>
#include <memory>

#include "gromacs/listed_forces/bonded.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/pbc_simd.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"
#include "gromacs/simd/vector_operations.h"
#include "gromacs/tables/forcetable.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

#include "listed_internal.h"

using namespace gmx; // TODO: Remove when this file is moved into gmx namespace

/*! \brief Issue a warning if a listed interaction is beyond a table limit */
static void warning_rlimit(const rvec* x, int ai, int aj, int* global_atom_index, real r, real rlimit)
{
    gmx_warning(
            "Listed nonbonded interaction between particles %d and %d\n"
            "at distance %.3f which is larger than the table limit %.3f nm.\n\n"
            "This is likely either a 1,4 interaction, or a listed interaction inside\n"
            "a smaller molecule you are decoupling during a free energy calculation.\n"
            "Since interactions at distances beyond the table cannot be computed,\n"
            "they are skipped until they are inside the table limit again. You will\n"
            "only see this message once, even if it occurs for several interactions.\n\n"
            "IMPORTANT: This should not happen in a stable simulation, so there is\n"
            "probably something wrong with your system. Only change the table-extension\n"
            "distance in the mdp file if you are really sure that is the reason.\n",
            glatnr(global_atom_index, ai),
            glatnr(global_atom_index, aj),
            r,
            rlimit);

    if (debug)
    {
        fprintf(debug,
                "%8f %8f %8f\n%8f %8f %8f\n1-4 (%d,%d) interaction not within cut-off! r=%g. "
                "Ignored\n",
                x[ai][XX],
                x[ai][YY],
                x[ai][ZZ],
                x[aj][XX],
                x[aj][YY],
                x[aj][ZZ],
                glatnr(global_atom_index, ai),
                glatnr(global_atom_index, aj),
                r);
    }
}

/*! \brief Compute the energy and force for a single pair interaction */
static real evaluate_single(real        r2,
                            real        tabscale,
                            const real* vftab,
                            real        tableStride,
                            real        qq,
                            real        c6,
                            real        c12,
                            real*       velec,
                            real*       vvdw)
{
    real rinv, r, rtab, eps, eps2, Y, F, Geps, Heps2, Fp, VVe, FFe, VVd, FFd, VVr, FFr, fscal;
    int  ntab;

    /* Do the tabulated interactions - first table lookup */
    rinv = gmx::invsqrt(r2);
    r    = r2 * rinv;
    rtab = r * tabscale;
    ntab = static_cast<int>(rtab);
    eps  = rtab - ntab;
    eps2 = eps * eps;
    ntab = static_cast<int>(tableStride * ntab);
    /* Electrostatics */
    Y     = vftab[ntab];
    F     = vftab[ntab + 1];
    Geps  = eps * vftab[ntab + 2];
    Heps2 = eps2 * vftab[ntab + 3];
    Fp    = F + Geps + Heps2;
    VVe   = Y + eps * Fp;
    FFe   = Fp + Geps + 2.0 * Heps2;
    /* Dispersion */
    Y     = vftab[ntab + 4];
    F     = vftab[ntab + 5];
    Geps  = eps * vftab[ntab + 6];
    Heps2 = eps2 * vftab[ntab + 7];
    Fp    = F + Geps + Heps2;
    VVd   = Y + eps * Fp;
    FFd   = Fp + Geps + 2.0 * Heps2;
    /* Repulsion */
    Y     = vftab[ntab + 8];
    F     = vftab[ntab + 9];
    Geps  = eps * vftab[ntab + 10];
    Heps2 = eps2 * vftab[ntab + 11];
    Fp    = F + Geps + Heps2;
    VVr   = Y + eps * Fp;
    FFr   = Fp + Geps + 2.0 * Heps2;

    *velec = qq * VVe;
    *vvdw  = c6 * VVd + c12 * VVr;

    fscal = -(qq * FFe + c6 * FFd + c12 * FFr) * tabscale * rinv;

    return fscal;
}

static inline real sixthRoot(const real r)
{
    return gmx::invsqrt(std::cbrt(r));
}

/*! \brief Compute the energy and force for a single pair interaction under FEP */
template<KernelSoftcoreType softcoreType>
static real free_energy_evaluate_single(real                                           r2,
                                        real                                           rCoulCutoff,
                                        const interaction_const_t::SoftCoreParameters& scParams,
                                        real                                           tabscale,
                                        const real*                                    vftab,
                                        real                                           tableStride,
                                        real                                           qqA,
                                        real                                           c6A,
                                        real                                           c12A,
                                        real                                           qqB,
                                        real                                           c6B,
                                        real                                           c12B,
                                        real                                           facel,
                                        const real                                     LFC[2],
                                        const real                                     LFV[2],
                                        const real                                     DLF[2],
                                        const real                                     lfac_coul[2],
                                        const real                                     lfac_vdw[2],
                                        const real dlfac_coul[2],
                                        const real dlfac_vdw[2],
                                        real*      velectot,
                                        real*      vvdwtot,
                                        real*      dvdl)
{
    real       rtab, eps, eps2, Y, F, Geps, Heps2, Fp, VV, FF, fscal;
    real       qq[2], c6[2], c12[2], sigma6[2], sigma_pow[2];
    real       alpha_coul_eff, alpha_vdw_eff, dvdl_coul_sum, dvdl_vdw_sum;
    real       rpinv, r_coul, r_vdw, velecsum, vvdwsum;
    real       fscal_vdw[2], fscal_elec[2];
    real       velec[2], vvdw[2];
    real       dvdl_elec[2], dvdl_vdw[2];
    real       gapsysScaleLinpointCoul, gapsysScaleLinpointVdW, gapsysSigma6VdW[2];
    real       rQ, rLJ;
    real       scaleDvdlRCoul;
    int        i, ntab;
    const real half = 0.5_real;
    const real one  = 1.0_real;
    const real two  = 2.0_real;

    qq[0]  = qqA;
    qq[1]  = qqB;
    c6[0]  = c6A;
    c6[1]  = c6B;
    c12[0] = c12A;
    c12[1] = c12B;

    const real rpm2 = r2 * r2;       /* r4 */
    const real rp   = rpm2 * r2;     /* r6 */
    const real r    = std::sqrt(r2); /* r1 */

    /* Loop over state A(0) and B(1) */
    for (i = 0; i < 2; i++)
    {
        if constexpr (softcoreType == KernelSoftcoreType::Beutler)
        {
            if ((c6[i] > 0) && (c12[i] > 0))
            {
                /* The c6 & c12 coefficients now contain the constants 6.0 and 12.0, respectively.
                 * Correct for this by multiplying with (1/12.0)/(1/6.0)=6.0/12.0=0.5.
                 */
                sigma6[i] = half * c12[i] / c6[i];
                if (sigma6[i] < scParams.sigma6Minimum) /* for disappearing coul and vdw with soft core at the same time */
                {
                    sigma6[i] = scParams.sigma6Minimum;
                }
            }
            else
            {
                sigma6[i] = scParams.sigma6WithInvalidSigma;
            }
            sigma_pow[i] = sigma6[i];
        }
        // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
        if constexpr (softcoreType == KernelSoftcoreType::Gapsys)
        {
            if ((c6[i] > 0) && (c12[i] > 0))
            {
                /* The c6 & c12 coefficients now contain the constants 6.0 and 12.0, respectively.
                 * Correct for this by multiplying with (1/12.0)/(1/6.0)=6.0/12.0=0.5.
                 */
                gapsysSigma6VdW[i] = half * c12[i] / c6[i];
            }
            else
            {
                gapsysSigma6VdW[i] = scParams.gapsysSigma6VdW;
            }
        }
    }

    // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
    if constexpr (softcoreType == KernelSoftcoreType::Beutler)
    {
        /* only use softcore if one of the states has a zero endstate - softcore is for avoiding infinities!*/
        if ((c12[0] > 0) && (c12[1] > 0))
        {
            alpha_vdw_eff  = 0;
            alpha_coul_eff = 0;
        }
        else
        {
            alpha_vdw_eff  = scParams.alphaVdw;
            alpha_coul_eff = scParams.alphaCoulomb;
        }
    }
    // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
    if constexpr (softcoreType == KernelSoftcoreType::Gapsys)
    {
        /* only use softcore if one of the states has a zero endstate - softcore is for avoiding infinities!*/
        if ((c12[0] > 0) && (c12[1] > 0))
        {
            gapsysScaleLinpointVdW  = 0;
            gapsysScaleLinpointCoul = 0;
        }
        else
        {
            gapsysScaleLinpointVdW  = scParams.gapsysScaleLinpointVdW;
            gapsysScaleLinpointCoul = scParams.gapsysScaleLinpointCoul;
        }
    }

    /* Loop over A and B states again */
    // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
    for (i = 0; i < 2; i++)
    {
        fscal_elec[i] = 0;
        fscal_vdw[i]  = 0;
        velec[i]      = 0;
        vvdw[i]       = 0;
        dvdl_elec[i]  = 0;
        dvdl_vdw[i]   = 0;
        rQ            = 0;
        rLJ           = 0;

        /* Only spend time on A or B state if it is non-zero */
        if ((qq[i] != 0) || (c6[i] != 0) || (c12[i] != 0))
        {
            /* Coulomb */
            if constexpr (softcoreType == KernelSoftcoreType::Beutler)
            {
                rpinv  = one / (alpha_coul_eff * lfac_coul[i] * sigma_pow[i] + rp);
                r_coul = sixthRoot(rpinv);
            }
            else
            {
                rpinv  = one / rp;
                r_coul = r;
            }

            // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
            if constexpr (softcoreType == KernelSoftcoreType::Gapsys)
            {
                if ((facel != 0) && (LFC[i] < 1))
                {
                    rQ = gmx::sixthroot(one - LFC[i]) * (one + std::fabs(qq[i] / facel));
                    rQ *= gapsysScaleLinpointCoul;
                }
                else
                {
                    rQ = 0;
                }
                scaleDvdlRCoul = 1;
                if (rQ > rCoulCutoff)
                {
                    rQ             = rCoulCutoff;
                    scaleDvdlRCoul = 0;
                }
            }

            // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
            if ((softcoreType == KernelSoftcoreType::Gapsys) && (r < rQ))
            {
                real rInvQ    = one / rQ;
                real constFac = qq[i] * rInvQ;
                real linFac   = constFac * r * rInvQ;
                real quadrFac = linFac * r * rInvQ;

                /* Computing Coulomb force and potential energy */
                fscal_elec[i] = -2 * quadrFac + 3 * linFac;
                fscal_elec[i] *= rpinv;

                velec[i] = quadrFac - 3 * (linFac - constFac);

                dvdl_elec[i] += scaleDvdlRCoul * DLF[i] * half * (LFC[i] / (1 - LFC[i]))
                                * (quadrFac - 2 * linFac + constFac);
            }
            else // Beutler, resp. hardcore
            {
                /* Electrostatics table lookup data */
                rtab = r_coul * tabscale;
                ntab = static_cast<int>(rtab);
                eps  = rtab - ntab;
                eps2 = eps * eps;
                ntab = static_cast<int>(tableStride * ntab);
                /* Electrostatics */
                Y             = vftab[ntab];
                F             = vftab[ntab + 1];
                Geps          = eps * vftab[ntab + 2];
                Heps2         = eps2 * vftab[ntab + 3];
                Fp            = F + Geps + Heps2;
                VV            = Y + eps * Fp;
                FF            = Fp + Geps + two * Heps2;
                velec[i]      = qq[i] * VV;
                fscal_elec[i] = -qq[i] * FF * r_coul * rpinv * tabscale;
            }

            /* Vdw */
            if constexpr (softcoreType == KernelSoftcoreType::Beutler)
            {
                rpinv = one / (alpha_vdw_eff * lfac_vdw[i] * sigma_pow[i] + rp);
                r_vdw = sixthRoot(rpinv);
            }
            else
            {
                rpinv = one / rp;
                r_vdw = r;
            }

            if constexpr (softcoreType == KernelSoftcoreType::Gapsys)
            {
                constexpr real c_twentySixSeventh = 26.0_real / 7.0_real;

                if (LFV[i] < 1)
                {

                    rLJ = gmx::sixthroot(c_twentySixSeventh * gapsysSigma6VdW[i] * (one - LFV[i]));
                    rLJ *= gapsysScaleLinpointVdW;
                }
                else
                {
                    rLJ = 0;
                }
            }

            // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
            if ((softcoreType == KernelSoftcoreType::Gapsys) && (r < rLJ))
            {
                // scaled values for c6 and c12
                real c6s, c12s;
                c6s  = c6[i] / 6.0_real;
                c12s = c12[i] / 12.0_real;

                /* Temporary variables for inverted values */
                real rInvLJ = one / rLJ;
                real rInv14, rInv13, rInv12;
                real rInv8, rInv7, rInv6;
                rInv6 = rInvLJ * rInvLJ * rInvLJ;
                rInv6 *= rInv6;
                rInv7  = rInv6 * rInvLJ;
                rInv8  = rInv7 * rInvLJ;
                rInv14 = c12s * rInv7 * rInv7 * r2;
                rInv13 = c12s * rInv7 * rInv6 * r;
                rInv12 = c12s * rInv6 * rInv6;
                rInv8 *= c6s * r2;
                rInv7 *= c6s * r;
                rInv6 *= c6s;

                /* Temporary variables for A and B */
                real quadrFac, linearFac, constFac;
                quadrFac  = 156 * rInv14 - 42 * rInv8;
                linearFac = 168 * rInv13 - 48 * rInv7;
                constFac  = 91 * rInv12 - 28 * rInv6;

                /* Computing LJ force and potential energy*/
                fscal_vdw[i] = -quadrFac + linearFac;
                fscal_vdw[i] *= rpinv;

                vvdw[i] = 0.5_real * quadrFac - linearFac + constFac;

                dvdl_vdw[i] += DLF[i] * 28 * (LFV[i] / (one - LFV[i]))
                               * ((6.5_real * rInv14 - rInv8) - (13 * rInv13 - 2 * rInv7)
                                  + (6.5_real * rInv12 - rInv6));
            }
            else // Beutler, resp. hardcore
            {
                /* Vdw table lookup data */
                rtab = r_vdw * tabscale;
                ntab = static_cast<int>(rtab);
                eps  = rtab - ntab;
                eps2 = eps * eps;
                ntab = 12 * ntab;
                /* Dispersion */
                Y            = vftab[ntab + 4];
                F            = vftab[ntab + 5];
                Geps         = eps * vftab[ntab + 6];
                Heps2        = eps2 * vftab[ntab + 7];
                Fp           = F + Geps + Heps2;
                VV           = Y + eps * Fp;
                FF           = Fp + Geps + two * Heps2;
                vvdw[i]      = c6[i] * VV;
                fscal_vdw[i] = -c6[i] * FF;

                /* Repulsion */
                Y     = vftab[ntab + 8];
                F     = vftab[ntab + 9];
                Geps  = eps * vftab[ntab + 10];
                Heps2 = eps2 * vftab[ntab + 11];
                Fp    = F + Geps + Heps2;
                VV    = Y + eps * Fp;
                FF    = Fp + Geps + two * Heps2;
                vvdw[i] += c12[i] * VV;
                fscal_vdw[i] -= c12[i] * FF;
                fscal_vdw[i] *= r_vdw * rpinv * tabscale;
            }
        }
    }
    /* Now we have velec[i], vvdw[i], and fscal[i] for both states */
    /* Assemble A and B states */
    velecsum      = 0;
    vvdwsum       = 0;
    dvdl_coul_sum = 0;
    dvdl_vdw_sum  = 0;
    fscal         = 0;
    for (i = 0; i < 2; i++)
    {
        velecsum += LFC[i] * velec[i];
        vvdwsum += LFV[i] * vvdw[i];

        fscal += (LFC[i] * fscal_elec[i] + LFV[i] * fscal_vdw[i]) * rpm2;

        if constexpr (softcoreType == KernelSoftcoreType::Gapsys)
        {
            dvdl_coul_sum += dvdl_elec[i];
            dvdl_vdw_sum += dvdl_vdw[i];
        }
        // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
        dvdl_coul_sum += velec[i] * DLF[i];
        dvdl_vdw_sum += vvdw[i] * DLF[i];
        if constexpr (softcoreType == KernelSoftcoreType::Beutler)
        {
            dvdl_coul_sum += LFC[i] * alpha_coul_eff * dlfac_coul[i] * fscal_elec[i] * sigma_pow[i];
            dvdl_vdw_sum += LFV[i] * alpha_vdw_eff * dlfac_vdw[i] * fscal_vdw[i] * sigma_pow[i];
        }
    }

    dvdl[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)] += dvdl_coul_sum;
    dvdl[static_cast<int>(FreeEnergyPerturbationCouplingType::Vdw)] += dvdl_vdw_sum;

    *velectot = velecsum;
    *vvdwtot  = vvdwsum;

    return fscal;
}

/*! \brief Calculate pair interactions, supports all types and conditions. */
template<BondedKernelFlavor flavor>
static real do_pairs_general(int                                 ftype,
                             int                                 nbonds,
                             const t_iatom                       iatoms[],
                             const t_iparams                     iparams[],
                             const rvec                          x[],
                             rvec4                               f[],
                             rvec                                fshift[],
                             const struct t_pbc*                 pbc,
                             const real*                         lambda,
                             real*                               dvdl,
                             gmx::ArrayRef<const real>           chargeA,
                             gmx::ArrayRef<const real>           chargeB,
                             gmx::ArrayRef<const bool>           atomIsPerturbed,
                             gmx::ArrayRef<const unsigned short> cENER,
                             int                                 numEnergyGroups,
                             const t_forcerec*                   fr,
                             gmx_grppairener_t*                  grppener,
                             int*                                global_atom_index)
{
    real            qq, c6, c12;
    rvec            dx;
    int             i, itype, ai, aj, gid;
    int             fshift_index;
    real            r2;
    real            fscal, velec, vvdw;
    real*           energygrp_elec;
    real*           energygrp_vdw;
    static gmx_bool warned_rlimit = FALSE;
    /* Free energy stuff */
    gmx_bool   bFreeEnergy;
    real       LFC[2], LFV[2], DLF[2], lfac_coul[2], lfac_vdw[2], dlfac_coul[2], dlfac_vdw[2];
    real       qqB, c6B, c12B;
    const real oneSixth = 1.0_real / 6.0_real;

    switch (ftype)
    {
        case F_LJ14:
        case F_LJC14_Q:
            energygrp_elec = grppener->energyGroupPairTerms[NonBondedEnergyTerms::Coulomb14].data();
            energygrp_vdw  = grppener->energyGroupPairTerms[NonBondedEnergyTerms::LJ14].data();
            break;
        case F_LJC_PAIRS_NB:
            energygrp_elec = grppener->energyGroupPairTerms[NonBondedEnergyTerms::CoulombSR].data();
            energygrp_vdw  = grppener->energyGroupPairTerms[NonBondedEnergyTerms::LJSR].data();
            break;
        default:
            energygrp_elec = nullptr; /* Keep compiler happy */
            energygrp_vdw  = nullptr; /* Keep compiler happy */
            gmx_fatal(FARGS, "Unknown function type %d in do_nonbonded14", ftype);
    }

    if (fr->efep != FreeEnergyPerturbationType::No)
    {
        if (atomIsPerturbed.empty())
        {
            // There are no perturbed atoms, chargeB is the same as chargeA
            chargeB = chargeA;
        }
        else
        {
            GMX_ASSERT(nbonds == 0 || !chargeB.empty(), "With perturbed atoms we need chargeB");
        }

        /* Lambda factor for state A=1-lambda and B=lambda */
        LFC[0] = 1.0 - lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)];
        LFV[0] = 1.0 - lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Vdw)];
        LFC[1] = lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)];
        LFV[1] = lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Vdw)];

        /*derivative of the lambda factor for state A and B */
        DLF[0] = -1;
        DLF[1] = 1;

        GMX_ASSERT(fr->ic->softCoreParameters, "We need soft-core parameters");
        const auto& scParams = *fr->ic->softCoreParameters;

        for (i = 0; i < 2; i++)
        {
            lfac_coul[i] = (scParams.lambdaPower == 2 ? (1 - LFC[i]) * (1 - LFC[i]) : (1 - LFC[i]));
            dlfac_coul[i] = DLF[i] * scParams.lambdaPower * oneSixth
                            * (scParams.lambdaPower == 2 ? (1 - LFC[i]) : 1);
            lfac_vdw[i]  = (scParams.lambdaPower == 2 ? (1 - LFV[i]) * (1 - LFV[i]) : (1 - LFV[i]));
            dlfac_vdw[i] = DLF[i] * scParams.lambdaPower * oneSixth
                           * (scParams.lambdaPower == 2 ? (1 - LFV[i]) : 1);
        }
    }

    /* TODO This code depends on the logic in tables.c that constructs
       the table layout, which should be made explicit in future
       cleanup. */
    GMX_ASSERT(
            etiNR == 3,
            "Pair-interaction code that uses GROMACS interaction tables supports exactly 3 tables");
    GMX_ASSERT(
            fr->pairsTable->interaction_ == TableInteraction::ElectrostaticVdwRepulsionVdwDispersion,
            "Pair interaction kernels need a table with Coulomb, repulsion and dispersion entries");

    const real epsfac = fr->ic->epsfac;

    bFreeEnergy = FALSE;
    for (i = 0; (i < nbonds);)
    {
        itype = iatoms[i++];
        ai    = iatoms[i++];
        aj    = iatoms[i++];
        gid   = GID(cENER[ai], cENER[aj], numEnergyGroups);

        /* Get parameters */
        switch (ftype)
        {
            case F_LJ14:
                bFreeEnergy =
                        (fr->efep != FreeEnergyPerturbationType::No
                         && ((!atomIsPerturbed.empty() && (atomIsPerturbed[ai] || atomIsPerturbed[aj]))
                             || iparams[itype].lj14.c6A != iparams[itype].lj14.c6B
                             || iparams[itype].lj14.c12A != iparams[itype].lj14.c12B));
                qq  = chargeA[ai] * chargeA[aj] * epsfac * fr->fudgeQQ;
                c6  = iparams[itype].lj14.c6A;
                c12 = iparams[itype].lj14.c12A;
                break;
            case F_LJC14_Q:
                qq = iparams[itype].ljc14.qi * iparams[itype].ljc14.qj * epsfac
                     * iparams[itype].ljc14.fqq;
                c6  = iparams[itype].ljc14.c6;
                c12 = iparams[itype].ljc14.c12;
                break;
            case F_LJC_PAIRS_NB:
                qq  = iparams[itype].ljcnb.qi * iparams[itype].ljcnb.qj * epsfac;
                c6  = iparams[itype].ljcnb.c6;
                c12 = iparams[itype].ljcnb.c12;
                break;
            default:
                /* Cannot happen since we called gmx_fatal() above in this case */
                qq = c6 = c12 = 0; /* Keep compiler happy */
                break;
        }

        /* To save flops in the optimized kernels, c6/c12 have 6.0/12.0 derivative prefactors
         * included in the general nfbp array now. This means the tables are scaled down by the
         * same factor, so when we use the original c6/c12 parameters from iparams[] they must
         * be scaled up.
         */
        c6 *= 6.0;
        c12 *= 12.0;

        /* Do we need to apply full periodic boundary conditions? */
        if (fr->bMolPBC)
        {
            fshift_index = pbc_dx_aiuc(pbc, x[ai], x[aj], dx);
        }
        else
        {
            fshift_index = c_centralShiftIndex;
            rvec_sub(x[ai], x[aj], dx);
        }
        r2 = norm2(dx);

        if (r2 >= fr->pairsTable->interactionRange * fr->pairsTable->interactionRange)
        {
            /* This check isn't race free. But it doesn't matter because if a race occurs the only
             * disadvantage is that the warning is printed twice */
            if (!warned_rlimit)
            {
                warning_rlimit(x, ai, aj, global_atom_index, std::sqrt(r2), fr->pairsTable->interactionRange);
                warned_rlimit = TRUE;
            }
            continue;
        }

        if (bFreeEnergy)
        {
            /* Currently free energy is only supported for F_LJ14, so no need to check for that if we got here */
            qqB  = chargeB[ai] * chargeB[aj] * epsfac * fr->fudgeQQ;
            c6B  = iparams[itype].lj14.c6B * 6.0;
            c12B = iparams[itype].lj14.c12B * 12.0;

            const auto& scParams = *fr->ic->softCoreParameters;
            if (scParams.softcoreType == SoftcoreType::Beutler)
            {
                if (scParams.alphaCoulomb == 0 && scParams.alphaVdw == 0)
                {
                    fscal = free_energy_evaluate_single<KernelSoftcoreType::None>(
                            r2,
                            fr->ic->rcoulomb,
                            *fr->ic->softCoreParameters,
                            fr->pairsTable->scale,
                            fr->pairsTable->data.data(),
                            fr->pairsTable->stride,
                            qq,
                            c6,
                            c12,
                            qqB,
                            c6B,
                            c12B,
                            epsfac,
                            LFC,
                            LFV,
                            DLF,
                            lfac_coul,
                            lfac_vdw,
                            dlfac_coul,
                            dlfac_vdw,
                            &velec,
                            &vvdw,
                            dvdl);
                }
                else
                {
                    fscal = free_energy_evaluate_single<KernelSoftcoreType::Beutler>(
                            r2,
                            fr->ic->rcoulomb,
                            *fr->ic->softCoreParameters,
                            fr->pairsTable->scale,
                            fr->pairsTable->data.data(),
                            fr->pairsTable->stride,
                            qq,
                            c6,
                            c12,
                            qqB,
                            c6B,
                            c12B,
                            epsfac,
                            LFC,
                            LFV,
                            DLF,
                            lfac_coul,
                            lfac_vdw,
                            dlfac_coul,
                            dlfac_vdw,
                            &velec,
                            &vvdw,
                            dvdl);
                }
            }
            else // Gapsys
            {
                if (scParams.gapsysScaleLinpointCoul == 0 && scParams.gapsysScaleLinpointVdW == 0)
                {
                    fscal = free_energy_evaluate_single<KernelSoftcoreType::None>(
                            r2,
                            fr->ic->rcoulomb,
                            *fr->ic->softCoreParameters,
                            fr->pairsTable->scale,
                            fr->pairsTable->data.data(),
                            fr->pairsTable->stride,
                            qq,
                            c6,
                            c12,
                            qqB,
                            c6B,
                            c12B,
                            epsfac,
                            LFC,
                            LFV,
                            DLF,
                            lfac_coul,
                            lfac_vdw,
                            dlfac_coul,
                            dlfac_vdw,
                            &velec,
                            &vvdw,
                            dvdl);
                }
                else
                {
                    fscal = free_energy_evaluate_single<KernelSoftcoreType::Gapsys>(
                            r2,
                            fr->ic->rcoulomb,
                            *fr->ic->softCoreParameters,
                            fr->pairsTable->scale,
                            fr->pairsTable->data.data(),
                            fr->pairsTable->stride,
                            qq,
                            c6,
                            c12,
                            qqB,
                            c6B,
                            c12B,
                            epsfac,
                            LFC,
                            LFV,
                            DLF,
                            lfac_coul,
                            lfac_vdw,
                            dlfac_coul,
                            dlfac_vdw,
                            &velec,
                            &vvdw,
                            dvdl);
                }
            }
        }
        else
        {
            /* Evaluate tabulated interaction without free energy */
            fscal = evaluate_single(r2,
                                    fr->pairsTable->scale,
                                    fr->pairsTable->data.data(),
                                    fr->pairsTable->stride,
                                    qq,
                                    c6,
                                    c12,
                                    &velec,
                                    &vvdw);
        }

        energygrp_elec[gid] += velec;
        energygrp_vdw[gid] += vvdw;
        svmul(fscal, dx, dx);

        /* Add the forces */
        rvec_inc(f[ai], dx);
        rvec_dec(f[aj], dx);

        if (computeVirial(flavor))
        {
            if (fshift_index != c_centralShiftIndex)
            {
                rvec_inc(fshift[fshift_index], dx);
                rvec_dec(fshift[c_centralShiftIndex], dx);
            }
        }
    }
    return 0.0;
}

/*! \brief Calculate pairs, only for plain-LJ + plain Coulomb normal type.
 *
 * This function is templated for real/SimdReal and for optimization.
 */
template<typename T, int pack_size, typename pbc_type>
static void do_pairs_simple(int                       nbonds,
                            const t_iatom             iatoms[],
                            const t_iparams           iparams[],
                            const rvec                x[],
                            rvec4                     f[],
                            const pbc_type            pbc,
                            gmx::ArrayRef<const real> charge,
                            const real                scale_factor)
{
    const int nfa1 = 1 + 2;

    T six(6);
    T twelve(12);
    T ef(scale_factor);

#if GMX_SIMD_HAVE_REAL
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t ai[pack_size];
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t aj[pack_size];
    alignas(GMX_SIMD_ALIGNMENT) real         coeff[3 * pack_size];
#else
    std::int32_t ai[pack_size];
    std::int32_t aj[pack_size];
    real         coeff[3 * pack_size];
#endif

    /* nbonds is #pairs*nfa1, here we step pack_size pairs */
    for (int i = 0; i < nbonds; i += pack_size * nfa1)
    {
        /* Collect atoms for pack_size pairs.
         * iu indexes into iatoms, we should not let iu go beyond nbonds.
         */
        int iu = i;
        for (int s = 0; s < pack_size; s++)
        {
            int itype = iatoms[iu];
            ai[s]     = iatoms[iu + 1];
            aj[s]     = iatoms[iu + 2];

            if (i + s * nfa1 < nbonds)
            {
                coeff[0 * pack_size + s] = iparams[itype].lj14.c6A;
                coeff[1 * pack_size + s] = iparams[itype].lj14.c12A;
                coeff[2 * pack_size + s] = charge[ai[s]] * charge[aj[s]];

                /* Avoid indexing the iatoms array out of bounds.
                 * We pad the coordinate indices with the last atom pair.
                 */
                if (iu + nfa1 < nbonds)
                {
                    iu += nfa1;
                }
            }
            else
            {
                /* Pad the coefficient arrays with zeros to get zero forces */
                coeff[0 * pack_size + s] = 0;
                coeff[1 * pack_size + s] = 0;
                coeff[2 * pack_size + s] = 0;
            }
        }

        /* Load the coordinates */
        T xi[DIM], xj[DIM];
        gatherLoadUTranspose<3>(reinterpret_cast<const real*>(x), ai, &xi[XX], &xi[YY], &xi[ZZ]);
        gatherLoadUTranspose<3>(reinterpret_cast<const real*>(x), aj, &xj[XX], &xj[YY], &xj[ZZ]);

        T c6  = load<T>(coeff + 0 * pack_size);
        T c12 = load<T>(coeff + 1 * pack_size);
        T qq  = load<T>(coeff + 2 * pack_size);

        /* We could save these operations by storing 6*C6,12*C12 */
        c6  = six * c6;
        c12 = twelve * c12;

        T dr[DIM];
        pbc_dx_aiuc(pbc, xi, xj, dr);

        T rsq   = dr[XX] * dr[XX] + dr[YY] * dr[YY] + dr[ZZ] * dr[ZZ];
        T rinv  = gmx::invsqrt(rsq);
        T rinv2 = rinv * rinv;
        T rinv6 = rinv2 * rinv2 * rinv2;

        /* Calculate the Coulomb force * r */
        T cfr = ef * qq * rinv;

        /* Calculate the LJ force * r and add it to the Coulomb part */
        T fr = gmx::fma(fms(c12, rinv6, c6), rinv6, cfr);

        T finvr = fr * rinv2;
        T fx    = finvr * dr[XX];
        T fy    = finvr * dr[YY];
        T fz    = finvr * dr[ZZ];

        /* Add the pair forces to the force array.
         * Note that here we might add multiple force components for some atoms
         * due to the SIMD padding. But the extra force components are zero.
         */
        transposeScatterIncrU<4>(reinterpret_cast<real*>(f), ai, fx, fy, fz);
        transposeScatterDecrU<4>(reinterpret_cast<real*>(f), aj, fx, fy, fz);
    }
}

/*! \brief Calculate all listed pair interactions */
void do_pairs(int                                 ftype,
              int                                 nbonds,
              const t_iatom                       iatoms[],
              const t_iparams                     iparams[],
              const rvec                          x[],
              rvec4                               f[],
              rvec                                fshift[],
              const struct t_pbc*                 pbc,
              const real*                         lambda,
              real*                               dvdl,
              gmx::ArrayRef<const real>           chargeA,
              gmx::ArrayRef<const real>           chargeB,
              gmx::ArrayRef<const bool>           atomIsPerturbed,
              gmx::ArrayRef<const unsigned short> cENER,
              const int                           numEnergyGroups,
              const t_forcerec*                   fr,
              const bool                          havePerturbedInteractions,
              const gmx::StepWorkload&            stepWork,
              gmx_grppairener_t*                  grppener,
              int*                                global_atom_index)
{
    if (ftype == F_LJ14 && fr->ic->vdwtype != VanDerWaalsType::User
        && !usingUserTableElectrostatics(fr->ic->eeltype) && !havePerturbedInteractions
        && (!stepWork.computeVirial && !stepWork.computeEnergy))
    {
        /* We use a fast code-path for plain LJ 1-4 without FEP.
         *
         * TODO: Add support for energies (straightforward) and virial
         * in the SIMD template. For the virial it's inconvenient to store
         * the force sums for the shifts and we should directly calculate
         * and sum the virial for the shifts. But we should do this
         * at once for the angles and dihedrals as well.
         */
#if GMX_SIMD_HAVE_REAL
        if (fr->use_simd_kernels)
        {
            alignas(GMX_SIMD_ALIGNMENT) real pbc_simd[9 * GMX_SIMD_REAL_WIDTH];
            set_pbc_simd(pbc, pbc_simd);

            do_pairs_simple<SimdReal, GMX_SIMD_REAL_WIDTH, const real*>(
                    nbonds, iatoms, iparams, x, f, pbc_simd, chargeA, fr->ic->epsfac * fr->fudgeQQ);
        }
        else
#endif
        {
            /* This construct is needed because pbc_dx_aiuc doesn't accept pbc=NULL */
            t_pbc        pbc_no;
            const t_pbc* pbc_nonnull;

            if (pbc != nullptr)
            {
                pbc_nonnull = pbc;
            }
            else
            {
                set_pbc(&pbc_no, PbcType::No, nullptr);
                pbc_nonnull = &pbc_no;
            }

            do_pairs_simple<real, 1, const t_pbc*>(
                    nbonds, iatoms, iparams, x, f, pbc_nonnull, chargeA, fr->ic->epsfac * fr->fudgeQQ);
        }
    }
    else if (stepWork.computeVirial)
    {
        do_pairs_general<BondedKernelFlavor::ForcesAndVirialAndEnergy>(ftype,
                                                                       nbonds,
                                                                       iatoms,
                                                                       iparams,
                                                                       x,
                                                                       f,
                                                                       fshift,
                                                                       pbc,
                                                                       lambda,
                                                                       dvdl,
                                                                       chargeA,
                                                                       chargeB,
                                                                       makeArrayRef(atomIsPerturbed),
                                                                       cENER,
                                                                       numEnergyGroups,
                                                                       fr,
                                                                       grppener,
                                                                       global_atom_index);
    }
    else
    {
        do_pairs_general<BondedKernelFlavor::ForcesAndEnergy>(ftype,
                                                              nbonds,
                                                              iatoms,
                                                              iparams,
                                                              x,
                                                              f,
                                                              fshift,
                                                              pbc,
                                                              lambda,
                                                              dvdl,
                                                              chargeA,
                                                              chargeB,
                                                              atomIsPerturbed,
                                                              cENER,
                                                              numEnergyGroups,
                                                              fr,
                                                              grppener,
                                                              global_atom_index);
    }
}
