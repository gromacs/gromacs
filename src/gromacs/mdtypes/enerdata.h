/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
#ifndef GMX_MDTYPES_TYPES_ENERDATA_H
#define GMX_MDTYPES_TYPES_ENERDATA_H

#include <array>
#include <utility>
#include <vector>

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

struct t_commrec;
struct t_lambda;

// The non-bonded energy terms accumulated for energy group pairs
enum class NonBondedEnergyTerms : int
{
    CoulombSR,
    LJSR,
    BuckinghamSR,
    Coulomb14,
    LJ14,
    Count
};

// Struct for accumulating non-bonded energies between energy group pairs
struct gmx_grppairener_t
{
    gmx_grppairener_t(int numEnergyGroups) : nener(numEnergyGroups * numEnergyGroups)
    {
        for (auto& term : energyGroupPairTerms)
        {
            term.resize(nener);
        }
    }

    void clear();

    int nener; /* The number of energy group pairs */
    gmx::EnumerationArray<NonBondedEnergyTerms, std::vector<real>> energyGroupPairTerms; /* Energy terms for each pair of groups */
};

//! Accumulates free-energy foreign lambda energies and dH/dlamba
class ForeignLambdaTerms
{
public:
    /*! \brief Constructor
     *
     * \param[in] allLambdas  The list of lambda values for all lambda components, can be nullptr
     */
    ForeignLambdaTerms(const gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, std::vector<double>>* allLambdas);

    //! Returns the number of foreign lambda values
    int numLambdas() const { return numLambdas_; }

    //! Returns the H(lambdaIndex) - H(lambda_current)
    double deltaH(int lambdaIndex) const { return energies_[1 + lambdaIndex] - energies_[0]; }

    /*! \brief Returns a list of partial energies, the part which depends on lambda),
     * current lambda in entry 0, foreign lambda i in entry 1+i
     *
     * Note: the potential terms needs to be finalized before calling this method.
     */
    gmx::ArrayRef<double> energies()
    {
        GMX_ASSERT(finalizedPotentialContributions_, "Should be finalized");
        return energies_;
    }

    /*! \brief Returns a list of partial energies, the part which depends on lambda),
     * current lambda in entry 0, foreign lambda i in entry 1+i
     *
     * Note: the potential terms needs to be finalized before calling this method.
     */
    gmx::ArrayRef<const double> energies() const
    {
        GMX_ASSERT(finalizedPotentialContributions_, "Should be finalized");
        return energies_;
    }

    /*! \brief Adds an energy and dV/dl constribution to lambda list index \p listIndex
     *
     * This should only be used for terms with non-linear dependence on lambda
     * The value passed as listIndex should be 0 for the current lambda
     * and 1+i for foreign lambda index i.
     */
    void accumulate(int listIndex, FreeEnergyPerturbationCouplingType couplingType, double energy, real dvdl)
    {
        GMX_ASSERT(!finalizedPotentialContributions_,
                   "Can only accumulate with an unfinalized object");

        energies_[listIndex] += energy;
        dhdl_[listIndex][couplingType] += dvdl;
    }

    /*! \brief Adds an energy and dV/dl constribution to lambda list index \p listIndex
     *
     * This should only be used for terms with non-linear dependence on lambda
     * The value passed as listIndex should be 0 for the current lambda
     * and 1+i for foreign lambda index i.
     */
    void accumulate(int    listIndex,
                    double energy,
                    const gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, real>& dvdl)
    {
        GMX_ASSERT(!finalizedPotentialContributions_,
                   "Can only accumulate with an unfinalized object");

        energies_[listIndex] += energy;
        for (auto fepct : gmx::EnumerationWrapper<FreeEnergyPerturbationCouplingType>{})
        {
            dhdl_[listIndex][fepct] += dvdl[fepct];
        }
    }

    /*! \brief Finalizes the potential (non-kinetic) terms
     *
     * Note: This can be called multiple times during the same force calculations
     * without affecting the results.
     *
     * \param[in] dvdlLinear  List of dV/dlambda contributions of size efptNR with depend linearly on lambda
     * \param[in] lambda      Lambda values for the efptNR contribution types
     * \param[in] fepvals     Free-energy parameters
     */
    void finalizePotentialContributions(
            const gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, double>& dvdlLinear,
            gmx::ArrayRef<const real>                                                lambda,
            const t_lambda&                                                          fepvals);

    /*! \brief Accumulates the kinetic and constraint free-energy contributions
     *
     * \param[in] energyTerms  List of energy terms, pass \p term in \p gmx_enerdata_t
     * \param[in] dhdlMass     The mass dependent contribution to dH/dlambda
     * \param[in] lambda       Lambda values for the efptNR contribution types
     * \param[in] fepvals      Free-energy parameters
     */
    void finalizeKineticContributions(gmx::ArrayRef<const real> energyTerms,
                                      double                    dhdlMass,
                                      gmx::ArrayRef<const real> lambda,
                                      const t_lambda&           fepvals);

    /*! \brief Returns a pair of lists of deltaH and dH/dlambda
     *
     * Both lists are of size numLambdas() and are indexed with the lambda index.
     *
     * Note: should only be called after the object has been finalized by a call to
     * accumulateLinearPotentialComponents() (is asserted).
     *
     * \param[in] cr  Communication record, used to reduce the terms when !=nullptr
     */
    std::pair<std::vector<double>, std::vector<double>> getTerms(const t_commrec* cr) const;

    //! Sets all terms to 0
    void zeroAllTerms();

private:
    //! As accumulate(), but for kinetic contributions
    void accumulateKinetic(int listIndex, double energy, double dhdl);

    //! Add a dH/dl contribution that does not depend on lambda to all foreign dH/dl terms
    void addConstantDhdl(FreeEnergyPerturbationCouplingType couplingType, double dhdl);

    //! The number of foreign lambdas
    int numLambdas_;
    //! The lambda vectors for all components
    const gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, std::vector<double>>* allLambdas_;
    //! Storage for foreign lambda energies
    std::vector<double> energies_;
    //! Storage for foreign lambda dH/dlambda
    std::vector<gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, double>> dhdl_;
    //! Tells whether all potential energy contributions have been accumulated
    bool finalizedPotentialContributions_ = false;
};

//! Struct for accumulating all potential energy terms and some kinetic energy terms
struct gmx_enerdata_t
{
    /*! \brief
     * Constructor with specific number of energy groups and lambdas.
     *
     * \param[in] numEnergyGroups Number of energy groups used.
     * \param[in] allLambdas      The lambda vectors for every component, can be nullptr
     */
    gmx_enerdata_t(int numEnergyGroups,
                   const gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, std::vector<double>>* allLambdas);

    //! The energies for all different interaction types
    std::array<real, F_NRE> term = { 0 };
    //! Energy group pair non-bonded energies
    struct gmx_grppairener_t grpp;
    //! Contributions to dV/dlambda with linear dependence on lambda
    gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, double> dvdl_lin = { 0 };
    //! Contributions to dV/dlambda with non-linear dependence on lambda
    gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, double> dvdl_nonlin = { 0 };
    /* The idea is that dvdl terms with linear lambda dependence will be added
     * automatically to enerpart_lambda. Terms with non-linear lambda dependence
     * should explicitly determine the energies at foreign lambda points
     * when n_lambda > 0. */

    //! Foreign lambda energies and dH/dl
    ForeignLambdaTerms foreignLambdaTerms;
};

#endif
