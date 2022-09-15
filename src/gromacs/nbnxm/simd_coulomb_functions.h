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
 * \brief Defines function for computing Coulomb interaction using SIMD
 *
 * The philosophy is as follows. The flavor of Coulomb interactions types
 * is set using constexpr variables. This is then used to select the appropriate
 * templated class with functions in this files through template specialization.
 * Functions in this file take C-style arrays of SIMD of size \p nR registers
 * as arguments, indicated by a suffix with the letter 'V', which are a list
 * of registers with one for each i-particle (4xM kernels) or
 * pair of i-particles (2xMM kernels) that have LJ interactions.
 * Note that we do not use functions for single SIMD registers because this limits
 * the instruction parallelism that compilers can extract.
 *
 * \author Berk Hess <hess@kth.se>
 */

#ifndef GMX_NBNXM_SIMD_COULOMB_FUNCTIONS_H
#define GMX_NBNXM_SIMD_COULOMB_FUNCTIONS_H

#include "gromacs/math/units.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/simd/simd.h"

#include "atomdata.h"

namespace gmx
{

//! List of type of Nbnxm kernel coulomb type implementations
enum class KernelCoulombType
{
    RF,              //!< Reaction-field, also used for plain cut-off
    EwaldAnalytical, //!< Ewald with analytical reciprocal contribution correction
    EwaldTabulated   //!< Ewald with tabulated reciprocal contribution correction
};

//! Base Coulomb calculator class, only specializations are used
template<KernelCoulombType coulombType>
class CoulombCalculator;

//! Specialized calculator for RF
template<>
class CoulombCalculator<KernelCoulombType::RF>
{
public:
    inline CoulombCalculator(const interaction_const_t& ic) :
        minusTwoTimesRFCoeff_(-2.0_real * ic.reactionFieldCoefficient),
        rfOffset_(ic.reactionFieldShift),
        selfEnergy_(0.5_real * ic.reactionFieldShift)
    {
    }

    //! Returns the self energy
    inline real selfEnergy() const { return selfEnergy_; }

    //! Returns the force
    template<int nR>
    inline std::array<SimdReal, nR> force(const std::array<SimdReal, nR>& rSquaredV,
                                          const std::array<SimdReal, nR> gmx_unused& dummyRInvV,
                                          const std::array<SimdReal, nR>&            rInvExclV,
                                          const std::array<SimdBool, nR> gmx_unused& withinCutoffV)
    {
        return genArr<nR>(
                [&](int i) { return fma(rSquaredV[i], minusTwoTimesRFCoeff_, rInvExclV[i]); });
    }

    //! Computes forces and energies without 1/r term for reaction-field
    template<int nR, std::size_t energySize>
    inline void forceAndCorrectionEnergy(const std::array<SimdReal, nR>& rSquaredV,
                                         const std::array<SimdReal, nR> gmx_unused& dummyRInvV,
                                         const std::array<SimdReal, nR>&            rInvExclV,
                                         const std::array<SimdBool, nR> gmx_unused& withinCutoffV,
                                         std::array<SimdReal, nR>&                  forceV,
                                         std::array<SimdReal, energySize>& correctionEnergyV)
    {
        forceV = genArr<nR>(
                [&](int i) { return fma(rSquaredV[i], minusTwoTimesRFCoeff_, rInvExclV[i]); });

        const SimdReal factor = SimdReal(0.5_real) * minusTwoTimesRFCoeff_;
        correctionEnergyV = genArr<nR>([&](int i) { return fma(rSquaredV[i], factor, rfOffset_); });
    }

private:
    //! -2 * RF-coeff
    const SimdReal minusTwoTimesRFCoeff_;
    //! The RF potential shift
    const SimdReal rfOffset_;
    //! The reaction-field self-energy
    const real selfEnergy_;
};

//! Specialized calculator for Ewald using an analytic approximation
template<>
class CoulombCalculator<KernelCoulombType::EwaldAnalytical>
{
public:
    inline CoulombCalculator(const interaction_const_t& ic) :
        beta_(ic.ewaldcoeff_q),
        betaSquared_(ic.ewaldcoeff_q * ic.ewaldcoeff_q),
        selfEnergy_(0.5_real * ic.ewaldcoeff_q * M_2_SQRTPI) // beta/sqrt(pi)
    {
    }

    //! Returns the self energy
    inline real selfEnergy() const { return selfEnergy_; }

    template<int nR>
    inline std::array<SimdReal, nR> force(const std::array<SimdReal, nR>& rSquaredV,
                                          const std::array<SimdReal, nR> gmx_unused& dummyRInvV,
                                          const std::array<SimdReal, nR>&            rInvExclV,
                                          const std::array<SimdBool, nR>&            withinCutoffV)
    {
        const auto brsqV = genArr<nR>(
                [&](int i) { return betaSquared_ * selectByMask(rSquaredV[i], withinCutoffV[i]); });

        const auto ewcorrV = genArr<nR>([&](int i) { return beta_ * pmeForceCorrection(brsqV[i]); });

        return genArr<nR>([&](int i) { return fma(ewcorrV[i], brsqV[i], rInvExclV[i]); });
    }

    //! Computes the Coulomb force and the correction energy for the Ewald reciprocal part
    template<int nR, std::size_t energySize>
    inline void forceAndCorrectionEnergy(const std::array<SimdReal, nR>&   rSquaredV,
                                         const std::array<SimdReal, nR>&   rInvV,
                                         const std::array<SimdReal, nR>&   rInvExclV,
                                         const std::array<SimdBool, nR>&   withinCutoffV,
                                         std::array<SimdReal, nR>&         forceV,
                                         std::array<SimdReal, energySize>& correctionEnergyV)
    {
        const auto brsqV = genArr<nR>(
                [&](int i) { return betaSquared_ * selectByMask(rSquaredV[i], withinCutoffV[i]); });

        const auto ewcorrV = genArr<nR>([&](int i) { return beta_ * pmeForceCorrection(brsqV[i]); });

        forceV = genArr<nR>([&](int i) { return fma(ewcorrV[i], brsqV[i], rInvExclV[i]); });

        correctionEnergyV =
                genArr<nR>([&](int i) { return beta_ * pmePotentialCorrection(brsqV[i]); });

        GMX_UNUSED_VALUE(rInvV);
    }

private:
    //! Ewald beta
    const SimdReal beta_;
    //! Ewald beta^2
    const SimdReal betaSquared_;
    //! The self energy of the reciprocal part
    const real selfEnergy_;
};

//! Specialized calculator for Ewald using tabulated functions
template<>
class CoulombCalculator<KernelCoulombType::EwaldTabulated>
{
public:
    inline CoulombCalculator(const interaction_const_t& ic) :
        invTableSpacing_(ic.coulombEwaldTables->scale),
        minusHalfTableSpacing_(-0.5_real / ic.coulombEwaldTables->scale),
        tableForce_(c_useTableFormatFDV0 ? ic.coulombEwaldTables->tableFDV0.data()
                                         : ic.coulombEwaldTables->tableF.data()),
        tablePotential_(ic.coulombEwaldTables->tableV.data()),
        selfEnergy_(0.5_real
                    * (c_useTableFormatFDV0 ? ic.coulombEwaldTables->tableFDV0[2]
                                            : ic.coulombEwaldTables->tableV[0]))
    {
    }

    //! Returns the self energy
    inline real selfEnergy() const { return selfEnergy_; }

    //! Returns the force
    template<int nR>
    inline std::array<SimdReal, nR> force(const std::array<SimdReal, nR>& rSquaredV,
                                          const std::array<SimdReal, nR>& rInvV,
                                          const std::array<SimdReal, nR>& rInvExclV,
                                          const std::array<SimdBool, nR> gmx_unused& withinCutoffV)
    {
        /* We use separate registers for r for tabulated Ewald and LJ to keep the code simpler */
        const auto rV = genArr<nR>([&](int i) { return rSquaredV[i] * rInvV[i]; });

        /* Convert r to scaled table units */
        const auto rScaledV = genArr<nR>([&](int i) { return rV[i] * invTableSpacing_; });

        /* Truncate scaled r to an int */
        std::array<SimdInt32, nR> tableIndexV;
        for (int i = 0; i < nR; i++)
        {
            tableIndexV[i] = cvttR2I(rScaledV[i]);
        }

        /* Convert r to scaled table units */
        const auto rScaledTruncatedV = genArr<nR>([&](int i) { return trunc(rScaledV[i]); });
        const auto rScaledFractionV =
                genArr<nR>([&](int i) { return rScaledV[i] - rScaledTruncatedV[i]; });

        /* Load and interpolate table forces and possibly energies.
         * Force and energy can be combined in one table, stride 4: FDV0
         * or in two separate tables with stride 1: F and V
         * Currently single precision uses FDV0, double F and V.
         */
        std::array<SimdReal, nR> coulombTable0V;
        std::array<SimdReal, nR> coulombTable1V;
        if constexpr (c_useTableFormatFDV0)
        {
            for (int i = 0; i < nR; i++)
            {
                gatherLoadBySimdIntTranspose<4>(
                        tableForce_, tableIndexV[i], &coulombTable0V[i], &coulombTable1V[i]);
            }
        }
        else
        {
            for (int i = 0; i < nR; i++)
            {
                gatherLoadUBySimdIntTranspose<1>(
                        tableForce_, tableIndexV[i], &coulombTable0V[i], &coulombTable1V[i]);
            }
            coulombTable1V = genArr<nR>([&](int i) { return coulombTable1V[i] - coulombTable0V[i]; });
        }
        const auto forceCorrectionV = genArr<nR>([&](int i) {
            return fma(rScaledFractionV[i], coulombTable1V[i], coulombTable0V[i]);
        });

        const auto forceV =
                genArr<nR>([&](int i) { return fnma(forceCorrectionV[i], rV[i], rInvExclV[i]); });

        return forceV;
    }

    //! Computes the Coulomb force and the Ewald reciprocal pot correction energy
    template<int nR, std::size_t energySize>
    inline void forceAndCorrectionEnergy(const std::array<SimdReal, nR>& rSquaredV,
                                         const std::array<SimdReal, nR>& rInvV,
                                         const std::array<SimdReal, nR>& rInvExclV,
                                         const std::array<SimdBool, nR> gmx_unused& withinCutoffV,
                                         std::array<SimdReal, nR>&                  forceV,
                                         std::array<SimdReal, energySize>& correctionEnergyV)
    {
        /* We use separate registers for r for tabulated Ewald and LJ to keep the code simpler */
        const auto rV = genArr<nR>([&](int i) { return rSquaredV[i] * rInvV[i]; });

        /* Convert r to scaled table units */
        const auto rScaledV = genArr<nR>([&](int i) { return rV[i] * invTableSpacing_; });

        /* Truncate scaled r to an int */
        std::array<SimdInt32, nR> tableIndexV;
        for (int i = 0; i < nR; i++)
        {
            tableIndexV[i] = cvttR2I(rScaledV[i]);
        }

        /* Convert r to scaled table units */
        const auto rScaledTruncatedV = genArr<nR>([&](int i) { return trunc(rScaledV[i]); });
        const auto rScaledFractionV =
                genArr<nR>([&](int i) { return rScaledV[i] - rScaledTruncatedV[i]; });

        /* Load and interpolate table forces and possibly energies.
         * Force and energy can be combined in one table, stride 4: FDV0
         * or in two separate tables with stride 1: F and V
         * Currently single precision uses FDV0, double F and V.
         */
        std::array<SimdReal, nR> coulombTable0V;
        std::array<SimdReal, nR> coulombTable1V;
        std::array<SimdReal, nR> coulombTablePotV;
        std::array<SimdReal, nR> dumV;
        if constexpr (c_useTableFormatFDV0)
        {
            for (int i = 0; i < nR; i++)
            {
                gatherLoadBySimdIntTranspose<4>(tableForce_,
                                                tableIndexV[i],
                                                &coulombTable0V[i],
                                                &coulombTable1V[i],
                                                &coulombTablePotV[i],
                                                &dumV[i]);
            }
        }
        else
        {
            for (int i = 0; i < nR; i++)
            {
                gatherLoadUBySimdIntTranspose<1>(
                        tableForce_, tableIndexV[i], &coulombTable0V[i], &coulombTable1V[i]);
                gatherLoadUBySimdIntTranspose<1>(
                        tablePotential_, tableIndexV[i], &coulombTablePotV[i], &dumV[i]);
            }
            coulombTable1V = genArr<nR>([&](int i) { return coulombTable1V[i] - coulombTable0V[i]; });
        }
        const auto forceCorrectionV = genArr<nR>([&](int i) {
            return fma(rScaledFractionV[i], coulombTable1V[i], coulombTable0V[i]);
        });

        forceV = genArr<nR>([&](int i) { return fnma(forceCorrectionV[i], rV[i], rInvExclV[i]); });

        correctionEnergyV = genArr<nR>([&](int i) {
            return fma((minusHalfTableSpacing_ * rScaledFractionV[i]),
                       (coulombTable0V[i] + forceCorrectionV[i]),
                       coulombTablePotV[i]);
        });
    }

private:
    //! 1 / table-spacing
    const SimdReal invTableSpacing_;
    //! -0.5 * table-spacing
    const SimdReal minusHalfTableSpacing_;
    //! The force table, either a list of forces or a list of F / DeltaF / V / 0
    const real* const tableForce_;
    //! The potential table
    const real* const tablePotential_;
    //! The self energy of the reciprocal part
    const real selfEnergy_;
};

} // namespace gmx

#endif // GMX_NBNXM_SIMD_COULOMB_FUNCTIONS_H
