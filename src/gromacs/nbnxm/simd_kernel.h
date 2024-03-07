/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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

#include "config.h"

#include "gromacs/pbcutil/ishift.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"
#include "gromacs/simd/vector_operations.h"

#include "kernel_common.h"
#include "nbnxm_simd.h"
#include "simd_coulomb_functions.h"
#include "simd_diagonal_masker.h"
#include "simd_energy_accumulator.h"
#include "simd_lennardjones_functions.h"
#include "simd_load_store_functions.h"

/*! \internal \file
 *
 * \brief
 * Declares and defines the NBNxM SIMD pair interaction kernel function.
 *
 * This definition is in an include file and not in a source file to enable
 * parallel compilation of different template instantiations in different
 * compilation units. The compilation units and a kernel lookup table are
 * currently manually generated offline by a python script (not at build time) any time
 * something changes in the interface.
 *
 * This kernel makes uses of many, templated, helper classes, so this file
 * can contain only the skeleton of the outer part of the kernel.
 * The inner loop of the SIMD kernel is defined in simd_kernel_inner.h,
 * which is currently included six times in the \p simdKernelInner()
 * function to optimize performance for different properties on pair list
 * entries.
 *
 * Many, but not all, combinations of the Cartesian product of the different
 * possible values of all template parameters are instantiated for this kernel.
 * In particular, \p vdwCutoffCheck is only enabled with Ewald type electrostatics
 * and only certain combinations of \p ljCombinationRule and \p vdwModifier are used.
 *
 * Which kernel layouts are instantiated depends on the extent of SIMD support of
 * the architecture and on which layouts are assumed to produce the best peformance.
 * Currently we have two kernel layouts:
 * - 4xM: this stores one j-cluster of width M in a SIMD register of width M and
 *        uses 4 registers for each i-atom variable, one for each i-atom.
 * - 2xMM: this stores two j-clusters of width M in a SIMD register of width 2*M and
 *         uses 2 registers for each i-atom variable, one for two i-atoms.
 * Most of the logic concerning the kernel layouts is not directly present in this
 * function, but is hidden in functions that are called through template specialization.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

namespace gmx
{

//! Whether we calculate shift forces, always true, because it's cheap anyhow
static constexpr bool sc_calculateShiftForces = true;

namespace
{

//! Energy accumulator without data members can be static
static EnergyAccumulator<false, false> s_energyAccumulator;

//! General template for EnergyAccumulator getter, only specializations are used
template<bool haveEnergyGroups, bool computeEnergies>
class EnergyAccumulatorGetter;

//! Specialized EnergyAccumulator getter for no energy output
template<>
class EnergyAccumulatorGetter<false, false>
{
public:
    inline gmx_unused EnergyAccumulatorGetter(nbnxn_atomdata_output_t gmx_unused* out) {}

    inline gmx_unused EnergyAccumulator<false, false>& get() { return s_energyAccumulator; }
};

//! Specialized EnergyAccumulator getter for single energy group output
template<>
class EnergyAccumulatorGetter<false, true>
{
public:
    inline gmx_unused EnergyAccumulatorGetter(nbnxn_atomdata_output_t* out) : out_(*out) {}

    inline gmx_unused EnergyAccumulator<false, true>& get()
    {
        return *out_.accumulatorSingleEnergies;
    }

private:
    nbnxn_atomdata_output_t& out_;
};

//! Specialized EnergyAccumulator getter for multiple energy groups output
template<>
class EnergyAccumulatorGetter<true, true>
{
public:
    inline gmx_unused EnergyAccumulatorGetter(nbnxn_atomdata_output_t* out) : out_(*out) {}

    inline gmx_unused EnergyAccumulator<true, true>& get()
    {
        return *out_.accumulatorGroupEnergies;
    }

private:
    nbnxn_atomdata_output_t& out_;
};

} // namespace

/*! \brief The actual NBNxM SIMD kernel
 *
 * \tparam kernelLayout    The kernel layout: either 2xMM or 4xM
 * \tparam coulombType     The type of Coulomb interaction
 * \tparam vdwCutoffCheck  Whether a separate check for the VdW cutoff is needed
 * \tparam ljCombinationRule  The combination rule for the LJ parameters
 * \tparam vdwModifier     The modifier for the LJ interactions
 * \tparam ljEwald         The type of LJ Ewald treatment, can be none
 * \tparam energyOutput    Which types of output are requested
 * \param[in] nbl          The cluster pair list
 * \param[in] nbat         Input data for atoms, including charges and LJ parameters
 * \param[in] ic           The interaction constants
 * \param[in] shift_vec    A list of PBC shift vectors
 * \param[in,out] out      Struct for accumulating forces, energies and virial
 */
template<KernelLayout         kernelLayout,
         KernelCoulombType    coulombType,
         VdwCutoffCheck       vdwCutoffCheck,
         LJCombinationRule    ljCombinationRule,
         InteractionModifiers vdwModifier,
         LJEwald              ljEwald,
         EnergyOutput         energyOutput>
void nbnxmKernelSimd(const NbnxnPairlistCpu*    nbl,
                     const nbnxn_atomdata_t*    nbat,
                     const interaction_const_t* ic,
                     const rvec*                shift_vec,
                     nbnxn_atomdata_output_t*   out)
{
    constexpr int c_numJClustersPerSimdRegister = (kernelLayout == KernelLayout::r2xMM ? 2 : 1);

    constexpr int c_iClusterSize = sc_iClusterSize(kernelLayout);
    constexpr int c_jClusterSize = sc_jClusterSize(kernelLayout);

    static_assert(c_numJClustersPerSimdRegister * c_jClusterSize == GMX_SIMD_REAL_WIDTH);

    // The ratio of cluster sizes
    constexpr KernelLayoutClusterRatio clusterRatio = kernelLayoutClusterRatio<kernelLayout>();

    // The stride of all atom data arrays
    constexpr int c_stride = std::max(c_iClusterSize, c_jClusterSize);

    /* The number of 'i' SIMD registers */
    static_assert(c_iClusterSize % c_numJClustersPerSimdRegister == 0);
    constexpr int nR = c_iClusterSize / c_numJClustersPerSimdRegister;

    constexpr bool haveVdwCutoffCheck = (vdwCutoffCheck != VdwCutoffCheck::No);

    constexpr bool haveLJEwaldGeometric = (ljEwald == LJEwald::CombGeometric);

    constexpr bool calculateEnergies = (energyOutput != EnergyOutput::None);
    constexpr bool useEnergyGroups   = (energyOutput == EnergyOutput::GroupPairs);

    /* Unpack pointers for output */
    real* f                 = out->f.data();
    real gmx_unused* fshift = out->fshift.data();

    const SimdReal zero_S(0.0);

#ifdef COUNT_PAIRS
    int npair = 0;
#endif

    const nbnxn_atomdata_t::Params& nbatParams = nbat->params();

    static_assert(!(haveLJEwaldGeometric && ljCombinationRule == LJCombinationRule::LorentzBerthelot),
                  "Can not have LJ-PME with LB combination rule");

    const real gmx_unused* gmx_restrict ljc;
    if constexpr (ljCombinationRule != LJCombinationRule::None || haveLJEwaldGeometric)
    {
        ljc = nbatParams.lj_comb.data();
    }
    const real gmx_unused* gmx_restrict nbfp_ptr;
    const int gmx_unused* gmx_restrict type;
    if constexpr (ljCombinationRule == LJCombinationRule::None)
    {
        /* No combination rule used */
        nbfp_ptr = nbatParams.nbfp_aligned.data();
        type     = nbatParams.type.data();
    }

    /* Set up the diagonal exclusion masks */
    const DiagonalMasker<nR, kernelLayout, clusterRatio> diagonalMasker(nbat->simdMasks());

#if GMX_DOUBLE && !GMX_SIMD_HAVE_INT32_LOGICAL
    const std::uint64_t* gmx_restrict exclusion_filter = nbat->simdMasks().exclusion_filter64.data();
#else
    const std::uint32_t* gmx_restrict exclusion_filter = nbat->simdMasks().exclusion_filter.data();
#endif

    /* Here we cast the exclusion filters from unsigned * to int * or real *.
     * Since we only check bits, the actual value they represent does not
     * matter, as long as both filter and mask data are treated the same way.
     */
    SimdBitMask exclusionFilterV[nR];
    for (int i = 0; i < nR; i++)
    {
#if GMX_SIMD_HAVE_INT32_LOGICAL
        exclusionFilterV[i] = load<SimdBitMask>(reinterpret_cast<const int*>(
                exclusion_filter + i * c_numJClustersPerSimdRegister * c_jClusterSize));
#else
        exclusionFilterV[i] = load<SimdBitMask>(reinterpret_cast<const real*>(
                exclusion_filter + i * c_numJClustersPerSimdRegister * c_jClusterSize));
#endif
    }

    CoulombCalculator<coulombType> coulombCalculator(*ic);

    gmx_unused SimdReal ewaldShift;
    if constexpr (coulombType != KernelCoulombType::RF && calculateEnergies)
    {
        ewaldShift = SimdReal(ic->sh_ewald);
    }

    /* LJ function constants, only actually needed with energies or potential switching */
    SimdReal sixth_S(1.0_real / 6.0_real);
    SimdReal twelveth_S(1.0_real / 12.0_real);

    static_assert(!(haveLJEwaldGeometric && vdwModifier != InteractionModifiers::PotShift),
                  "LJ-PME only supports potential-shift");

    LennardJonesCalculator<calculateEnergies, vdwModifier> ljCalculator(*ic);

    std::array<SimdReal, haveLJEwaldGeometric ? 5 : 0> gmx_unused ljEwaldParams;
    real                                                          lj_ewaldcoeff6_6;
    if constexpr (haveLJEwaldGeometric)
    {
        ljEwaldParams[0]          = SimdReal(1.0_real);
        ljEwaldParams[1]          = SimdReal(0.5_real);
        const real lj_ewaldcoeff2 = ic->ewaldcoeff_lj * ic->ewaldcoeff_lj;
        lj_ewaldcoeff6_6          = lj_ewaldcoeff2 * lj_ewaldcoeff2 * lj_ewaldcoeff2 / 6;
        ljEwaldParams[2]          = SimdReal(lj_ewaldcoeff2);
        ljEwaldParams[3]          = SimdReal(lj_ewaldcoeff6_6);
        /* Determine the grid potential at the cut-off */
        ljEwaldParams[4] = ic->sh_lj_ewald;
    }

    /* The kernel either supports rcoulomb = rvdw or rcoulomb >= rvdw */
    const SimdReal cutoffSquared(ic->rcoulomb * ic->rcoulomb);
    SimdReal       vdwCutoffSquared;
    if constexpr (haveVdwCutoffCheck)
    {
        vdwCutoffSquared = SimdReal(ic->rvdw * ic->rvdw);
    }

    const SimdReal minDistanceSquared(c_nbnxnMinDistanceSquared);

    const real* gmx_restrict q        = nbatParams.q.data();
    const real               facel    = ic->epsfac;
    const real* gmx_restrict shiftvec = shift_vec[0];
    const real* gmx_restrict x        = nbat->x().data();

    EnergyAccumulator<useEnergyGroups, calculateEnergies>& energyAccumulator =
            EnergyAccumulatorGetter<useEnergyGroups, calculateEnergies>(out).get();

    const nbnxn_cj_t* l_cj = nbl->cj.list_.data();

    for (const nbnxn_ci_t& ciEntry : nbl->ci)
    {
        const int ish    = (ciEntry.shift & NBNXN_CI_SHIFT);
        const int ish3   = ish * 3;
        const int cjind0 = ciEntry.cj_ind_start;
        const int cjind1 = ciEntry.cj_ind_end;
        const int ci     = ciEntry.ci;
        const int ci_sh  = (ish == c_centralShiftIndex ? ci : -1);

        // Load the periodic shift vector for the i-atoms
        const SimdReal iShiftX(shiftvec[ish3]);
        const SimdReal iShiftY(shiftvec[ish3 + 1]);
        const SimdReal iShiftZ(shiftvec[ish3 + 2]);

        int sci;
        int scix;
        int sci2;
        if constexpr (c_jClusterSize <= 4)
        {
            sci  = ci * c_stride;
            scix = sci * DIM;
            sci2 = sci * 2;
        }
        else
        {
            sci  = (ci >> 1) * c_stride;
            scix = sci * DIM + (ci & 1) * (c_stride >> 1);
            sci2 = sci * 2 + (ci & 1) * (c_stride >> 1);
            sci += (ci & 1) * (c_stride >> 1);
        }

        /* We have 5 LJ/C combinations, but use only three inner loops,
         * as the other combinations are unlikely and/or not much faster:
         * inner half-LJ + C for half-LJ + C / no-LJ + C
         * inner LJ + C      for full-LJ + C
         * inner LJ          for full-LJ + no-C / half-LJ + no-C
         */
        const bool do_LJ   = ((ciEntry.shift & NBNXN_CI_DO_LJ(0)) != 0);
        const bool do_coul = ((ciEntry.shift & NBNXN_CI_DO_COUL(0)) != 0);
        const bool half_LJ = (((ciEntry.shift & NBNXN_CI_HALF_LJ(0)) != 0) || !do_LJ) && do_coul;

        energyAccumulator.template initICluster<c_iClusterSize>(ci);

        if constexpr (calculateEnergies)
        {
            // Compute self interaction energies, when present
            const bool do_self = haveLJEwaldGeometric || do_coul;

            if (do_self && l_cj[ciEntry.cj_ind_start].cj == cjFromCi<clusterRatio>(ci_sh))
            {
                if (do_coul)
                {
                    const real Vc_sub_self = coulombCalculator.selfEnergy();

                    for (int ia = 0; ia < c_iClusterSize; ia++)
                    {
                        const real qi = q[sci + ia];

                        energyAccumulator.addCoulombEnergy(ia, -facel * qi * qi * Vc_sub_self);
                    }
                }

                if constexpr (haveLJEwaldGeometric)
                {
                    for (int ia = 0; ia < c_iClusterSize; ia++)
                    {
                        real c6_i =
                                nbatParams.nbfp[nbatParams.type[sci + ia] * (nbatParams.numTypes + 1) * 2]
                                / 6;
                        energyAccumulator.addVdwEnergy(ia, 0.5 * c6_i * lj_ewaldcoeff6_6);
                    }
                }
            }

        } // calulateEnergies

        /* Load i atom data */
        const int sciy = scix + c_stride;
        const int sciz = sciy + c_stride;

        const auto ixV =
                genArr<nR>([&](int i) { return loadIAtomData<kernelLayout>(x, scix, i) + iShiftX; });
        const auto iyV =
                genArr<nR>([&](int i) { return loadIAtomData<kernelLayout>(x, sciy, i) + iShiftY; });
        const auto izV =
                genArr<nR>([&](int i) { return loadIAtomData<kernelLayout>(x, sciz, i) + iShiftZ; });

        std::array<SimdReal, nR> chargeIV;
        if (do_coul)
        {
            chargeIV =
                    genArr<nR>([&](int i) { return facel * loadIAtomData<kernelLayout>(q, sci, i); });
        }

        constexpr bool c_ljCombLB = (ljCombinationRule == LJCombinationRule::LorentzBerthelot);
        // Note that when half_lj==true we actually only need nR/2 LJ parameters
        std::array<SimdReal, c_ljCombLB ? nR : 0> gmx_unused halfSigmaIV;
        std::array<SimdReal, c_ljCombLB ? nR : 0> gmx_unused sqrtEpsilonIV;
        std::array<SimdReal, (ljCombinationRule == LJCombinationRule::Geometric || haveLJEwaldGeometric) ? nR : 0> gmx_unused c6GeomV;
        std::array<SimdReal, (ljCombinationRule == LJCombinationRule::Geometric) ? nR : 0> gmx_unused c12GeomV;
        std::array<const real*, (ljCombinationRule == LJCombinationRule::None) ? c_iClusterSize : 0> gmx_unused nbfpI;
        if constexpr (c_ljCombLB)
        {
            for (int i = 0; i < nR; i++)
            {
                halfSigmaIV[i]   = loadIAtomData<kernelLayout>(ljc, sci2, i);
                sqrtEpsilonIV[i] = loadIAtomData<kernelLayout>(ljc, sci2 + c_stride, i);
            }
        }
        else if constexpr (ljCombinationRule == LJCombinationRule::Geometric)
        {
            for (int i = 0; i < nR; i++)
            {
                c6GeomV[i]  = loadIAtomData<kernelLayout>(ljc, sci2, i);
                c12GeomV[i] = loadIAtomData<kernelLayout>(ljc, sci2 + c_stride, i);
            }
        }
        else
        {
            const int numTypes = nbatParams.numTypes;
            for (int i = 0; i < c_iClusterSize; i++)
            {
                nbfpI[i] = nbfp_ptr + type[sci + i] * numTypes * c_simdBestPairAlignment;
            }
        }
        if constexpr (haveLJEwaldGeometric && ljCombinationRule != LJCombinationRule::Geometric)
        {
            /* We need the geometrically combined C6 for the PME grid correction */
            for (int i = 0; i < nR; i++)
            {
                c6GeomV[i] = loadIAtomData<kernelLayout>(ljc, sci2, i);
            }
        }

        /* Declare and clear i atom forces */
        auto forceIXV = genArr<nR>([&](int gmx_unused i) {
            SimdReal tmp = setZero();
            return tmp;
        });
        auto forceIYV = genArr<nR>([&](int gmx_unused i) {
            SimdReal tmp = setZero();
            return tmp;
        });
        auto forceIZV = genArr<nR>([&](int gmx_unused i) {
            SimdReal tmp = setZero();
            return tmp;
        });


        int cjind = cjind0;

        /* Currently all kernels use (at least half) LJ */
        if (half_LJ)
        {
            /* Coulomb: all i-atoms, LJ: first half i-atoms */
            constexpr bool            c_calculateCoulombInteractions = true;
            constexpr ILJInteractions c_iLJInteractions              = ILJInteractions::Half;
            {
                constexpr bool c_needToCheckExclusions = true;
                while (cjind < cjind1 && nbl->cj.excl(cjind) != NBNXN_INTERACTION_MASK_ALL)
                {
#include "simd_kernel_inner.h"
                    cjind++;
                }
            }
            {
                constexpr bool c_needToCheckExclusions = false;
                for (; (cjind < cjind1); cjind++)
                {
#include "simd_kernel_inner.h"
                }
            }
        }
        else if (do_coul)
        {
            /* Coulomb: all i-atoms, LJ: all i-atoms */
            constexpr bool            c_calculateCoulombInteractions = true;
            constexpr ILJInteractions c_iLJInteractions              = ILJInteractions::All;
            {
                constexpr bool c_needToCheckExclusions = true;
                while (cjind < cjind1 && nbl->cj.excl(cjind) != NBNXN_INTERACTION_MASK_ALL)
                {
#include "simd_kernel_inner.h"
                    cjind++;
                }
            }
            {
                constexpr bool c_needToCheckExclusions = false;
                for (; (cjind < cjind1); cjind++)
                {
#include "simd_kernel_inner.h"
                }
            }
        }
        else
        {
            /* Coulomb: none, LJ: all i-atoms */
            constexpr bool            c_calculateCoulombInteractions = false;
            constexpr ILJInteractions c_iLJInteractions              = ILJInteractions::All;
            {
                constexpr bool c_needToCheckExclusions = true;
                while (cjind < cjind1 && nbl->cj.excl(cjind) != NBNXN_INTERACTION_MASK_ALL)
                {
#include "simd_kernel_inner.h"
                    cjind++;
                }
            }
            {
                constexpr bool c_needToCheckExclusions = false;
                for (; (cjind < cjind1); cjind++)
                {
#include "simd_kernel_inner.h"
                }
            }
        }
        /* Add accumulated i-forces to the force array */
        real fShiftX;
        real fShiftY;
        real fShiftZ;
        if constexpr (c_numJClustersPerSimdRegister == 1)
        {
            fShiftX = reduceIncr4ReturnSum(f + scix, forceIXV[0], forceIXV[1], forceIXV[2], forceIXV[3]);
            fShiftY = reduceIncr4ReturnSum(f + sciy, forceIYV[0], forceIYV[1], forceIYV[2], forceIYV[3]);
            fShiftZ = reduceIncr4ReturnSum(f + sciz, forceIZV[0], forceIZV[1], forceIZV[2], forceIZV[3]);
        }
        else
        {
            fShiftX = reduceIncr4ReturnSumHsimd(f + scix, forceIXV[0], forceIXV[1]);
            fShiftY = reduceIncr4ReturnSumHsimd(f + sciy, forceIYV[0], forceIYV[1]);
            fShiftZ = reduceIncr4ReturnSumHsimd(f + sciz, forceIZV[0], forceIZV[1]);
        }

        if constexpr (sc_calculateShiftForces)
        {
            fshift[ish3 + 0] += fShiftX;
            fshift[ish3 + 1] += fShiftY;
            fshift[ish3 + 2] += fShiftZ;
        }

        energyAccumulator.reduceIEnergies(do_coul);

        /* Outer loop uses 6 flops/iteration */
    }

#ifdef COUNT_PAIRS
    printf("atom pairs %d\n", npair);
#endif
}

} // namespace gmx
