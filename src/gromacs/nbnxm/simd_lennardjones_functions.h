/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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

/*! \internal
 * \ingroup __module_nbnxm
 *
 * \brief Defines function for computing Lennard-Jones interaction using SIMD
 *
 * The philosophy is as follows. The different flavors of LJ interactions types
 * and calculating energies or not are set using constexpr variables. These are
 * then used to select the appropriate templated functions in this files through
 * template specialization. Functions in this file take C-style arrays of SIMD
 * of size \p nR registers as arguments, indicated by a suffix with the letter 'V',
 * which are a list of registers with one for each i-particle (4xM kernels) or
 * pair of i-particles (2xMM kernels) that have LJ interactions. This can be half
 * the number of total registers when only part of the i-atoms have LJ.
 * Note that we do not use functions for single SIMD registers because this limits
 * the instruction parallelism that compilers can extract.
 *
 * A specialized templated LennardJonesCalculator class is used.
 * For each cluster pair, \p forceC612() or \p forceSigmaEpsilon() is called,
 * depending on the combination rule. Additionaly \p addLennardJonesEwaldCorrections()
 * needs to be called when using Ewald for Lennard-Jones.
 *
 * Note that when atoms can get very close to each other, which we assume only
 * happens when atoms are excluded, the masking template parameter should be set
 * to true to avoid overflows when calculating r^-6.
 *
 * Note that only plain or potential-shifted LJ interactions are supported with
 * Lorentz-Berthelot combination rules. For switched LJ interations choose no
 * combination rule.
 *
 * \author Berk Hess <hess@kth.se>
 */

#ifndef GMX_NBNXM_SIMD_LENNARDJONES_FUNCTIONS_H
#define GMX_NBNXM_SIMD_LENNARDJONES_FUNCTIONS_H

#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/simd/simd.h"

#include "atomdata.h"

namespace gmx
{

//! The fraction of i-particles for which LJ interactions need to be computed
enum class ILJInteractions
{
    All,  //!< all i-particles
    Half, //!< the first half of the i-particles
    None  //!< none of i-particles
};

//! Base LJ calculator class, only specializations are used
template<bool calculateEnergies, InteractionModifiers vdwModifier>
class LennardJonesCalculator;

//! Computes r^-6 and r^-12, masked when requested
template<int nR, bool maskInteractions, std::size_t inputSize, std::size_t interactSize>
inline void rInvSixAndRInvTwelve(const std::array<SimdReal, inputSize>&    rInvSquaredV,
                                 const std::array<SimdBool, interactSize>& interactV,
                                 std::array<SimdReal, nR>&                 rInvSixV,
                                 std::array<SimdReal, nR>&                 rInvTwelveV)
{
    rInvSixV = genArr<nR>([&](int i) { return rInvSquaredV[i] * rInvSquaredV[i] * rInvSquaredV[i]; });

    // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
    if constexpr (maskInteractions)
    {
        rInvSixV = genArr<nR>([&](int i) { return selectByMask(rInvSixV[i], interactV[i]); });
    }

    // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
    rInvTwelveV = genArr<nR>([&](int i) { return rInvSixV[i] * rInvSixV[i]; });
}

//! Returns F*r and optionally the potential for LJ with (un)shifted potential with sigma/epsilon
template<int nR, bool maskInteractions, bool haveCutoffCheck, bool calculateEnergies, std::size_t inputSize, std::size_t interactSize, std::size_t vljvSize>
inline void lennardJonesInteractionsSigmaEpsilon(const std::array<SimdReal, inputSize>&    rInvV,
                                                 const std::array<SimdBool, interactSize>& interactV,
                                                 const SimdBool* const           withinCutoffV,
                                                 const std::array<SimdReal, nR>& sigmaV,
                                                 const std::array<SimdReal, nR>& epsilonV,
                                                 const SimdReal                  dispersionShift,
                                                 const SimdReal                  repulsionShift,
                                                 const SimdReal                  sixth,
                                                 const SimdReal                  twelfth,
                                                 std::array<SimdReal, nR>&       frLJV,
                                                 std::array<SimdReal, vljvSize>& vLJV)
{
    static_assert(inputSize >= nR);
    static_assert(!calculateEnergies || vljvSize == nR);

    const auto sigmaInvRV  = genArr<nR>([&](int i) { return sigmaV[i] * rInvV[i]; });
    const auto sigmaInvR2V = genArr<nR>([&](int i) { return sigmaInvRV[i] * sigmaInvRV[i]; });
    auto       sigmaInvR6V =
            genArr<nR>([&](int i) { return sigmaInvR2V[i] * sigmaInvR2V[i] * sigmaInvR2V[i]; });
    if constexpr (maskInteractions)
    {
        sigmaInvR6V = genArr<nR>([&](int i) { return selectByMask(sigmaInvR6V[i], interactV[i]); });
    }

    // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
    if constexpr (haveCutoffCheck)
    {
        sigmaInvR6V =
                genArr<nR>([&](int i) { return selectByMask(sigmaInvR6V[i], withinCutoffV[i]); });
    }
    else
    {
        GMX_UNUSED_VALUE(withinCutoffV);
    }

    const auto frLJ6V  = genArr<nR>([&](int i) { return epsilonV[i] * sigmaInvR6V[i]; });
    const auto frLJ12V = genArr<nR>([&](int i) { return frLJ6V[i] * sigmaInvR6V[i]; });
    frLJV              = genArr<nR>([&](int i) { return frLJ12V[i] - frLJ6V[i]; });

    if constexpr (calculateEnergies)
    {
        /* We need C6 and C12 to calculate the LJ potential shift */
        const auto sigma2V = genArr<nR>([&](int i) { return sigmaV[i] * sigmaV[i]; });
        const auto sigma6V = genArr<nR>([&](int i) { return sigma2V[i] * sigma2V[i] * sigma2V[i]; });
        const auto c6V     = genArr<nR>([&](int i) { return epsilonV[i] * sigma6V[i]; });
        const auto c12V    = genArr<nR>([&](int i) { return c6V[i] * sigma6V[i]; });

        /* Calculate the LJ energies, with constant potential shift */
        vLJV = genArr<nR>([&](int i) { return sixth * fma(c6V[i], dispersionShift, frLJ6V[i]); });
        vLJV = genArr<nR>([&](int i) {
            return fms(twelfth, fma(c12V[i], repulsionShift, frLJ12V[i]), vLJV[i]);
        });
    }
    else
    {
        GMX_UNUSED_VALUE(dispersionShift);
        GMX_UNUSED_VALUE(repulsionShift);
        GMX_UNUSED_VALUE(sixth);
        GMX_UNUSED_VALUE(twelfth);
    }
}

//! Specialized calculator for LJ with potential shift and no energy calculation
template<>
class LennardJonesCalculator<false, InteractionModifiers::PotShift>
{
public:
    inline LennardJonesCalculator(const interaction_const_t gmx_unused& ic) {}

    //! Computes F*r for LJ with (un)shifted potential with C6/C12 parameters
    template<int nR, bool maskInteractions, std::size_t inputSize, std::size_t interactSize, std::size_t vljvSize>
    inline void forceC6C12(const std::array<SimdReal, inputSize>&    rSquaredV,
                           const std::array<SimdReal, inputSize>&    rInvV,
                           const std::array<SimdReal, inputSize>&    rInvSquaredV,
                           const std::array<SimdBool, interactSize>& interactV,
                           const std::array<SimdReal, nR>&           c6V,
                           const std::array<SimdReal, nR>&           c12V,
                           SimdReal                                  sixth,
                           SimdReal                                  twelfth,
                           std::array<SimdReal, nR>&                 frLJV,
                           std::array<SimdReal, vljvSize>&           vLJV)
    {
        std::array<SimdReal, nR> rInvSixV;
        std::array<SimdReal, nR> rInvTwelveV;
        rInvSixAndRInvTwelve<nR, maskInteractions>(rInvSquaredV, interactV, rInvSixV, rInvTwelveV);

        frLJV = genArr<nR>([&](int i) { return fms(c12V[i], rInvTwelveV[i], c6V[i] * rInvSixV[i]); });

        GMX_UNUSED_VALUE(rSquaredV);
        GMX_UNUSED_VALUE(rInvV);
        GMX_UNUSED_VALUE(rInvSquaredV);
        GMX_UNUSED_VALUE(sixth);
        GMX_UNUSED_VALUE(twelfth);
        GMX_UNUSED_VALUE(vLJV);
    }

    template<int nR, bool maskInteractions, bool haveCutoffCheck, std::size_t inputSize, std::size_t interactSize, std::size_t vljvSize>
    inline void forceSigmaEpsilon(const std::array<SimdReal, inputSize>&    rInvV,
                                  const std::array<SimdBool, interactSize>& interactV,
                                  SimdBool*                                 withinCutoffV,
                                  const std::array<SimdReal, nR>&           sigmaV,
                                  const std::array<SimdReal, nR>&           epsilonV,
                                  SimdReal                                  sixth,
                                  SimdReal                                  twelfth,
                                  std::array<SimdReal, nR>&                 frLJV,
                                  std::array<SimdReal, vljvSize>&           vLJV)
    {
        const SimdReal dummy = setZero();

        lennardJonesInteractionsSigmaEpsilon<nR, maskInteractions, haveCutoffCheck, false>(
                rInvV, interactV, withinCutoffV, sigmaV, epsilonV, dummy, dummy, sixth, twelfth, frLJV, vLJV);
    }
};

//! Specialized calculator for LJ with potential shift and energy calculation
template<>
class LennardJonesCalculator<true, InteractionModifiers::PotShift>
{
public:
    inline LennardJonesCalculator(const interaction_const_t& ic) :
        dispersionShift_(ic.dispersion_shift.cpot), repulsionShift_(ic.repulsion_shift.cpot)
    {
    }

    //! Computes F*r and the potential for LJ with (un)shifted potential with C6/C12 parameters
    template<int nR, bool maskInteractions, std::size_t inputSize, std::size_t interactSize, std::size_t vljvSize>
    inline void forceC6C12(const std::array<SimdReal, inputSize>&    rSquaredV,
                           const std::array<SimdReal, inputSize>&    rInvV,
                           const std::array<SimdReal, inputSize>&    rInvSquaredV,
                           const std::array<SimdBool, interactSize>& interactV,
                           const std::array<SimdReal, nR>&           c6V,
                           const std::array<SimdReal, nR>&           c12V,
                           SimdReal                                  sixth,
                           SimdReal                                  twelfth,
                           std::array<SimdReal, nR>&                 frLJV,
                           std::array<SimdReal, vljvSize>&           vLJV)
    {
        static_assert(inputSize >= nR);
        static_assert(vljvSize == nR);

        std::array<SimdReal, nR> frLJ6V;
        std::array<SimdReal, nR> frLJ12V;
        rInvSixAndRInvTwelve<nR, maskInteractions>(rInvSquaredV, interactV, frLJ6V, frLJ12V);

        frLJ6V  = genArr<nR>([&](int i) { return c6V[i] * frLJ6V[i]; });
        frLJ12V = genArr<nR>([&](int i) { return c12V[i] * frLJ12V[i]; });
        frLJV   = genArr<nR>([&](int i) { return frLJ12V[i] - frLJ6V[i]; });

        vLJV = genArr<nR>([&](int i) { return sixth * fma(c6V[i], dispersionShift_, frLJ6V[i]); });
        vLJV = genArr<nR>([&](int i) {
            return fms(twelfth, fma(c12V[i], repulsionShift_, frLJ12V[i]), vLJV[i]);
        });

        GMX_UNUSED_VALUE(rSquaredV);
        GMX_UNUSED_VALUE(rInvV);
    }

    template<int nR, bool maskInteractions, bool haveCutoffCheck, std::size_t inputSize, std::size_t interactSize, std::size_t vljvSize>
    inline void forceSigmaEpsilon(const std::array<SimdReal, inputSize>&    rInvV,
                                  const std::array<SimdBool, interactSize>& interactV,
                                  SimdBool*                                 withinCutoffV,
                                  const std::array<SimdReal, nR>&           sigmaV,
                                  const std::array<SimdReal, nR>&           epsilonV,
                                  SimdReal                                  sixth,
                                  SimdReal                                  twelfth,
                                  std::array<SimdReal, nR>&                 frLJV,
                                  std::array<SimdReal, vljvSize>&           vLJV)
    {
        lennardJonesInteractionsSigmaEpsilon<nR, maskInteractions, haveCutoffCheck, true>(
                rInvV, interactV, withinCutoffV, sigmaV, epsilonV, dispersionShift_, repulsionShift_, sixth, twelfth, frLJV, vLJV);
    }

private:
    const SimdReal dispersionShift_;
    const SimdReal repulsionShift_;
};

//! Computes (r - r_switch), (r - r_switch)^2 and (r - r_switch)^2 * r
template<int nR, std::size_t inputSize>
inline void computeForceSwitchVariables(const std::array<SimdReal, inputSize>& rSquaredV,
                                        const std::array<SimdReal, inputSize>& rInvV,
                                        SimdReal                               rSwitch,
                                        std::array<SimdReal, nR>&              rSwitchedV,
                                        std::array<SimdReal, nR>&              rSwitchedSquaredV,
                                        std::array<SimdReal, nR>& rSwitchedSquaredTimesRV)
{
    static_assert(inputSize >= nR);

    auto rV           = genArr<nR>([&](int i) { return rSquaredV[i] * rInvV[i]; });
    rSwitchedV        = genArr<nR>([&](int i) { return max(rV[i] - rSwitch, setZero()); });
    rSwitchedSquaredV = genArr<nR>([&](int i) { return rSwitchedV[i] * rSwitchedV[i]; });

    rSwitchedSquaredTimesRV = genArr<nR>([&](int i) { return rSwitchedSquaredV[i] * rV[i]; });
}

//! Adds the force switch term to \p force
inline SimdReal addLJForceSwitch(SimdReal force,
                                 SimdReal rSwitched,
                                 SimdReal rSwitchedSquaredTimesR,
                                 SimdReal c2,
                                 SimdReal c3)
{
    return fma(fma(c3, rSwitched, c2), rSwitchedSquaredTimesR, force);
}

//! Returns the LJ force switch function for the potential
inline SimdReal ljForceSwitchPotential(SimdReal rSwitched,
                                       SimdReal rSwitchedSquaredTimesR,
                                       SimdReal c0,
                                       SimdReal c3,
                                       SimdReal c4)
{
    return fma(fma(c4, rSwitched, c3), rSwitchedSquaredTimesR * rSwitched, c0);
}

//! Specialized calculator for LJ with force switch
template<bool calculateEnergies>
class LennardJonesCalculator<calculateEnergies, InteractionModifiers::ForceSwitch>
{
public:
    inline LennardJonesCalculator(const interaction_const_t& ic) :
        rSwitch_(ic.rvdw_switch),
        dispersionShiftC2_(ic.dispersion_shift.c2),
        dispersionShiftC3_(ic.dispersion_shift.c3),
        repulsionShiftC2_(ic.repulsion_shift.c2),
        repulsionShiftC3_(ic.repulsion_shift.c3)
    {
        if constexpr (calculateEnergies)
        {
            SimdReal mthird_S(-1.0_real / 3.0_real);
            SimdReal mfourth_S(-1.0_real / 4.0_real);

            potentialParams_[0] = mthird_S * ic.dispersion_shift.c2;
            potentialParams_[1] = mfourth_S * ic.dispersion_shift.c3;
            potentialParams_[2] = SimdReal(ic.dispersion_shift.cpot / 6.0_real);
            potentialParams_[3] = mthird_S * ic.repulsion_shift.c2;
            potentialParams_[4] = mfourth_S * ic.repulsion_shift.c3;
            potentialParams_[5] = SimdReal(ic.repulsion_shift.cpot / 12.0_real);
        }
    }

    //! Computes F*r and optionally the potential for LJ with force switch and C6/C12 parameters
    template<int nR, bool maskInteractions, std::size_t inputSize, std::size_t interactSize, std::size_t vljvSize>
    inline void forceC6C12(const std::array<SimdReal, inputSize>&    rSquaredV,
                           const std::array<SimdReal, inputSize>&    rInvV,
                           const std::array<SimdReal, inputSize>&    rInvSquaredV,
                           const std::array<SimdBool, interactSize>& interactV,
                           const std::array<SimdReal, nR>&           c6V,
                           const std::array<SimdReal, nR>&           c12V,
                           SimdReal                                  sixth,
                           SimdReal                                  twelfth,
                           std::array<SimdReal, nR>&                 frLJV,
                           std::array<SimdReal, vljvSize>&           vLJV)
    {
        static_assert(inputSize >= nR);
        static_assert(!calculateEnergies || vljvSize == nR);

        std::array<SimdReal, nR> rInvSixV;
        std::array<SimdReal, nR> rInvTwelveV;
        rInvSixAndRInvTwelve<nR, maskInteractions>(rInvSquaredV, interactV, rInvSixV, rInvTwelveV);

        std::array<SimdReal, nR> rSwitchedV;
        std::array<SimdReal, nR> rSwitchedSquaredV;
        std::array<SimdReal, nR> rSwitchedSquaredTimesRV;
        computeForceSwitchVariables<nR>(
                rSquaredV, rInvV, rSwitch_, rSwitchedV, rSwitchedSquaredV, rSwitchedSquaredTimesRV);

        frLJV = genArr<nR>([&](int i) {
            return c6V[i]
                   * addLJForceSwitch(rInvSixV[i],
                                      rSwitchedV[i],
                                      rSwitchedSquaredTimesRV[i],
                                      dispersionShiftC2_,
                                      dispersionShiftC3_);
        });

        frLJV = genArr<nR>([&](int i) {
            return c12V[i]
                           * addLJForceSwitch(rInvTwelveV[i],
                                              rSwitchedV[i],
                                              rSwitchedSquaredTimesRV[i],
                                              repulsionShiftC2_,
                                              repulsionShiftC3_)
                   - frLJV[i];
        });

        if constexpr (calculateEnergies)
        {
            vLJV = genArr<nR>([&](int i) {
                return c6V[i]
                       * fma(sixth,
                             rInvSixV[i],
                             ljForceSwitchPotential(rSwitchedV[i],
                                                    rSwitchedSquaredV[i],
                                                    potentialParams_[2],
                                                    potentialParams_[0],
                                                    potentialParams_[1]));
            });

            vLJV = genArr<nR>([&](int i) {
                return c12V[i]
                               * fma(twelfth,
                                     rInvSixV[i] * rInvSixV[i],
                                     ljForceSwitchPotential(rSwitchedV[i],
                                                            rSwitchedSquaredV[i],
                                                            potentialParams_[5],
                                                            potentialParams_[3],
                                                            potentialParams_[4]))
                       - vLJV[i];
            });
        }
        else
        {
            GMX_UNUSED_VALUE(sixth);
            GMX_UNUSED_VALUE(twelfth);
        }
    }

    //! Never used, only declared to enable compilation
    template<int nR, bool maskInteractions, bool haveCutoffCheck, std::size_t inputSize, std::size_t interactSize, std::size_t vljvSize>
    void forceSigmaEpsilon(const std::array<SimdReal, inputSize>&    rInvV,
                           const std::array<SimdBool, interactSize>& interactV,
                           SimdBool*                                 withinCutoffV,
                           const std::array<SimdReal, nR>&           sigmaV,
                           const std::array<SimdReal, nR>&           epsilonV,
                           SimdReal                                  sixth,
                           SimdReal                                  twelfth,
                           std::array<SimdReal, nR>&                 frLJV,
                           std::array<SimdReal, vljvSize>&           vLJV);

private:
    const SimdReal                                  rSwitch_;
    const SimdReal                                  dispersionShiftC2_;
    const SimdReal                                  dispersionShiftC3_;
    const SimdReal                                  repulsionShiftC2_;
    const SimdReal                                  repulsionShiftC3_;
    std::array<SimdReal, calculateEnergies ? 6 : 0> potentialParams_;
};

//! Computes (r - r_switch) and (r - r_switch)^2
template<int nR, std::size_t inputSize>
inline void computePotentialSwitchVariables(const std::array<SimdReal, inputSize>& rSquaredV,
                                            const std::array<SimdReal, inputSize>& rInvV,
                                            SimdReal                               rSwitch,
                                            std::array<SimdReal, nR>&              rSwitchedV,
                                            std::array<SimdReal, nR>& rSwitchedSquaredV)
{
    static_assert(inputSize >= nR);

    auto rV           = genArr<nR>([&](int i) { return rSquaredV[i] * rInvV[i]; });
    rSwitchedV        = genArr<nR>([&](int i) { return max(rV[i] - rSwitch, setZero()); });
    rSwitchedSquaredV = genArr<nR>([&](int i) { return rSwitchedV[i] * rSwitchedV[i]; });
}

//! Returns the potential switch function
inline SimdReal potentialSwitchFunction(SimdReal rsw, SimdReal rsw2, SimdReal c3, SimdReal c4, SimdReal c5)
{
    return fma(fma(fma(c5, rsw, c4), rsw, c3), rsw2 * rsw, SimdReal(1.0_real));
}

//! Returns the derivative of the potential switch function
inline SimdReal potentialSwitchFunctionDerivative(SimdReal rsw, SimdReal rsw2, SimdReal c2, SimdReal c3, SimdReal c4)
{
    return fma(fma(c4, rsw, c3), rsw, c2) * rsw2;
}

//! Specialized calculator for LJ with potential switch
template<bool calculateEnergies>
class LennardJonesCalculator<calculateEnergies, InteractionModifiers::PotSwitch>
{
public:
    inline LennardJonesCalculator(const interaction_const_t& ic) :
        rSwitch_(ic.rvdw_switch),
        c3_(ic.vdw_switch.c3),
        c4_(ic.vdw_switch.c4),
        c5_(ic.vdw_switch.c5),
        c3times3_(3.0_real * ic.vdw_switch.c3),
        c4times4_(4.0_real * ic.vdw_switch.c4),
        c5times5_(5.0_real * ic.vdw_switch.c5)
    {
    }

    //! Computes F*r and optionally the potential for LJ with potential switch and C6/C12 parameters
    template<int nR, bool maskInteractions, std::size_t inputSize, std::size_t interactSize, std::size_t vljvSize>
    inline void forceC6C12(const std::array<SimdReal, inputSize>&    rSquaredV,
                           const std::array<SimdReal, inputSize>&    rInvV,
                           const std::array<SimdReal, inputSize>&    rInvSquaredV,
                           const std::array<SimdBool, interactSize>& interactV,
                           const std::array<SimdReal, nR>&           c6V,
                           const std::array<SimdReal, nR>&           c12V,
                           SimdReal                                  sixth,
                           SimdReal                                  twelfth,
                           std::array<SimdReal, nR>&                 frLJV,
                           std::array<SimdReal, vljvSize>&           vLJV)
    {
        static_assert(inputSize >= nR);
        static_assert(!calculateEnergies || vljvSize == nR);

        // Compute the plain LJ force
        std::array<SimdReal, nR> frLJ6V;
        std::array<SimdReal, nR> frLJ12V;
        rInvSixAndRInvTwelve<nR, maskInteractions>(rInvSquaredV, interactV, frLJ6V, frLJ12V);

        frLJ6V  = genArr<nR>([&](int i) { return c6V[i] * frLJ6V[i]; });
        frLJ12V = genArr<nR>([&](int i) { return c12V[i] * frLJ12V[i]; });
        frLJV   = genArr<nR>([&](int i) { return frLJ12V[i] - frLJ6V[i]; });

        // We always need the potential, since it is needed for the switched force
        auto vLJTmpV = genArr<nR>([&](int i) { return sixth * frLJ6V[i]; });
        vLJTmpV      = genArr<nR>([&](int i) { return fms(twelfth, frLJ12V[i], vLJTmpV[i]); });

        std::array<SimdReal, nR> rSwitchedV;
        std::array<SimdReal, nR> rSwitchedSquaredV;
        computePotentialSwitchVariables<nR>(rSquaredV, rInvV, rSwitch_, rSwitchedV, rSwitchedSquaredV);

        const auto switchV  = genArr<nR>([&](int i) {
            return potentialSwitchFunction(rSwitchedV[i], rSwitchedSquaredV[i], c3_, c4_, c5_);
        });
        const auto dSwitchV = genArr<nR>([&](int i) {
            return potentialSwitchFunctionDerivative(
                    rSwitchedV[i], rSwitchedSquaredV[i], c3times3_, c4times4_, c5times5_);
        });

        frLJV = genArr<nR>([&](int i) {
            SimdReal r = rSquaredV[i] * rInvV[i];
            return fnma(dSwitchV[i] * vLJTmpV[i], r, switchV[i] * frLJV[i]);
        });

        if constexpr (calculateEnergies)
        {
            vLJV = genArr<nR>([&](int i) { return vLJV[i] = switchV[i] * vLJTmpV[i]; });
        }
    }

    //! Never used, only declared to enable compilation
    template<int nR, bool maskInteractions, bool haveCutoffCheck, std::size_t inputSize, std::size_t interactSize, std::size_t vljvSize>
    void forceSigmaEpsilon(const std::array<SimdReal, inputSize>&    rInvV,
                           const std::array<SimdBool, interactSize>& interactV,
                           SimdBool*                                 withinCutoffV,
                           const std::array<SimdReal, nR>&           sigmaV,
                           const std::array<SimdReal, nR>&           epsilonV,
                           SimdReal                                  sixth,
                           SimdReal                                  twelfth,
                           std::array<SimdReal, nR>&                 frLJV,
                           std::array<SimdReal, vljvSize>&           vLJV);

private:
    const SimdReal rSwitch_;
    const SimdReal c3_;
    const SimdReal c4_;
    const SimdReal c5_;
    const SimdReal c3times3_;
    const SimdReal c4times4_;
    const SimdReal c5times5_;
};

//! Adds the Ewald long-range correction for r^-6
template<int nR, bool maskInteractions, bool calculateEnergies, std::size_t inputSize, std::size_t interactSize, std::size_t ljepSize, std::size_t vljvSize>
inline void addLennardJonesEwaldCorrections(const std::array<SimdReal, inputSize>&    rSquaredV,
                                            const std::array<SimdReal, inputSize>&    rInvSquaredV,
                                            const std::array<SimdBool, interactSize>& interactV,
                                            const SimdBool*                           withinCutoffV,
                                            const std::array<SimdReal, nR>&           c6GridV,
                                            const std::array<SimdReal, ljepSize>&     ljEwaldParams,
                                            SimdReal                                  sixth,
                                            std::array<SimdReal, nR>&                 frLJV,
                                            std::array<SimdReal, vljvSize>&           vLJV)
{
    static_assert(inputSize >= nR);
    static_assert(ljepSize == 5);
    static_assert(!calculateEnergies || vljvSize == nR);

    /* Recalculate r^-6 not masked for exclusions.
     * Note that we could reuse the previously calculated r^-6 which is unmasked
     * for exclusions when not calculating energies.
     */
    const auto rInvSixV =
            genArr<nR>([&](int i) { return rInvSquaredV[i] * rInvSquaredV[i] * rInvSquaredV[i]; });

    /* Mask for the cut-off to avoid overflow of cr2^2 */
    const auto rSquaredMaskedV = genArr<nR>(
            [&](int i) { return ljEwaldParams[2] * selectByMask(rSquaredV[i], withinCutoffV[i]); });

    // Unsafe version of our exp() should be fine, since these arguments should never
    // be smaller than -127 for any reasonable choice of cutoff or ewald coefficients.
    const auto expRSquaredMaskedV =
            genArr<nR>([&](int i) { return exp<MathOptimization::Unsafe>(-rSquaredMaskedV[i]); });

    /* 1 + cr2 + 1/2*cr2^2 */
    const auto polyV = genArr<nR>([&](int i) {
        return fma(fma(ljEwaldParams[1], rSquaredMaskedV[i], ljEwaldParams[0]),
                   rSquaredMaskedV[i],
                   ljEwaldParams[0]);
    });

    /* We calculate LJ F*r = (6*C6)*(r^-6 - F_mesh/6), we use:
     * r^-6*cexp*(1 + cr2 + cr2^2/2 + cr2^3/6) = cexp*(r^-6*poly + c^6/6)
     */
    frLJV = genArr<nR>([&](int i) {
        return fma(c6GridV[i],
                   fnma(expRSquaredMaskedV[i], fma(rInvSixV[i], polyV[i], ljEwaldParams[3]), rInvSixV[i]),
                   frLJV[i]);
    });

    if constexpr (calculateEnergies)
    {
        std::array<SimdReal, nR> shiftMaskedV;
        if constexpr (maskInteractions)
        {
            shiftMaskedV =
                    genArr<nR>([&](int i) { return selectByMask(ljEwaldParams[4], interactV[i]); });
        }
        else
        {
            shiftMaskedV = genArr<nR>([&](int gmx_unused i) { return ljEwaldParams[4]; });
        }

        vLJV = genArr<nR>([&](int i) {
            return fma(sixth * c6GridV[i],
                       fma(rInvSixV[i],
                           fnma(expRSquaredMaskedV[i], polyV[i], ljEwaldParams[0]),
                           shiftMaskedV[i]),
                       vLJV[i]);
        });
    }
    else
    {
        GMX_UNUSED_VALUE(sixth);
    }
}

} // namespace gmx

#endif // GMX_NBNXM_SIMD_LENNARDJONES_FUNCTIONS_H
