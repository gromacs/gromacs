/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * \brief This file defines functions to reduce and check
 * exclusion counts.
 *
 * \ingroup module_nbnxm
 */

#include "gmxpre.h"

#include "gromacs/nbnxm/exclusionchecker.h"

#include <functional>
#include <utility>
#include <vector>

#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/observablesreducer.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/listoflists.h"

/*! \brief Data to help check local topology construction
 *
 * Partitioning could incorrectly miss a bonded interaction.
 * However, checking for that requires a global communication
 * stage, which does not otherwise happen during partitioning. So,
 * for performance, we do that alongside the first global energy
 * reduction after a new DD is made. These variables handle
 * whether the check happens, its input for this domain, output
 * across all domains, and the expected value it should match. */
class ExclusionChecker::Impl
{
public:
    //! Constructor
    Impl(const t_commrec* cr, const gmx_mtop_t& mtop);

    //! Checks the count passed with the expected number and exits with a fatal error at mismatch
    void check(int numTotalPerturbedExclusionsFound);

    //! Object used when reporting that exclusions are missing
    //! {
    //! Communication record
    const t_commrec* cr_;
    //! }

    /*! \brief View used for computing the global number of bonded interactions.
     *
     * Can be written any time, but that is only useful when followed
     * by a call of the callbackToRequireReduction. Useful to read
     * only from the callback that the ObservablesReducer will later
     * make after reduction. */
    gmx::ArrayRef<double> reductionBuffer_;
    /*! \brief Callback used after repartitioning to require reduction
     * of numBondedInteractionsToReduce so that the total number of
     * bonded interactions can be checked. */
    gmx::ObservablesReducerBuilder::CallbackToRequireReduction callbackToRequireReduction_;
    /*! \brief The expected number of global non-bonded perturbed pair exclusions */
    int expectedNumGlobalPerturbedExclusions_;
};


/*! \brief Compute the number of exclusions involving perturbed atoms
 *
 * \param[in] mtop              The global system topology
 *
 * \returns the total number of exclusions exclusions (one-way) that
 *          involve at least one perturbed atom, including self-exclusions.
 */
static int computeNumGlobalPerturbedExclusions(const gmx_mtop_t& mtop)
{
    int numPerturbedExclusions = 0;

    for (const auto& molblock : mtop.molblock)
    {
        const gmx_moltype_t&         molt  = mtop.moltype[molblock.type];
        const gmx::ListOfLists<int>& excls = molt.excls;
        const t_atoms&               atoms = molt.atoms;

        int numPerturbedExclusionsInMol = 0;
        for (int atomI = 0; atomI < ssize(excls); atomI++)
        {
            const bool atomIIsPerturbed = PERTURBED(atoms.atom[atomI]);
            for (int atomJ : excls[atomI])
            {
                if (atomJ >= atomI && (atomIIsPerturbed || PERTURBED(atoms.atom[atomJ])))
                {
                    numPerturbedExclusionsInMol++;
                }
            }
        }

        numPerturbedExclusions += molblock.nmol * numPerturbedExclusionsInMol;
    }

    return numPerturbedExclusions;
}

ExclusionChecker::Impl::Impl(const t_commrec* cr, const gmx_mtop_t& mtop) :
    cr_(cr), expectedNumGlobalPerturbedExclusions_(computeNumGlobalPerturbedExclusions(mtop))
{
}

ExclusionChecker::ExclusionChecker(const t_commrec*                cr,
                                   const gmx_mtop_t&               mtop,
                                   gmx::ObservablesReducerBuilder* observablesReducerBuilder) :
    impl_(std::make_unique<Impl>(cr, mtop))
{
    if (cr == nullptr || !havePPDomainDecomposition(cr))
    {
        // No reduction required
        return;
    }

    GMX_RELEASE_ASSERT(observablesReducerBuilder,
                       "With DD an ObservablesReducerBuilder is required");

    Impl*                                               impl = impl_.get();
    gmx::ObservablesReducerBuilder::CallbackFromBuilder callbackFromBuilder =
            [impl](gmx::ObservablesReducerBuilder::CallbackToRequireReduction c, gmx::ArrayRef<double> v) {
                impl->callbackToRequireReduction_ = std::move(c);
                impl->reductionBuffer_            = v;
            };

    // Make the callback that runs afer reduction.
    gmx::ObservablesReducerBuilder::CallbackAfterReduction callbackAfterReduction = [impl](gmx::Step /*step*/) {
        // Pass the total after reduction to the check
        impl->check(impl->reductionBuffer_[0]);
    };

    observablesReducerBuilder->addSubscriber(
            1, std::move(callbackFromBuilder), std::move(callbackAfterReduction));
}

ExclusionChecker::~ExclusionChecker() = default;

ExclusionChecker::ExclusionChecker(ExclusionChecker&&) noexcept = default;

ExclusionChecker& ExclusionChecker::operator=(ExclusionChecker&& other) noexcept
{
    impl_ = std::move(other.impl_);
    return *this;
}

void ExclusionChecker::Impl::check(const int numTotalPerturbedExclusionsFound)
{
    if (numTotalPerturbedExclusionsFound != expectedNumGlobalPerturbedExclusions_)
    {
        // Give error and exit
        gmx_fatal_collective(
                FARGS,
                cr_->mpi_comm_mygroup,
                MAIN(cr_),
                "There are %d perturbed, excluded non-bonded pair interactions beyond the "
                "pair-list "
                "cut-off, which is not supported. This can happen because the system is "
                "unstable or because intra-molecular interactions at long distances are "
                "excluded. If the "
                "latter is the case, you can try to increase nstlist or rlist to avoid this."
                "The error is likely triggered by the use of couple-intramol=no "
                "and the maximal distance in the decoupled molecule exceeding rlist.",
                expectedNumGlobalPerturbedExclusions_ - numTotalPerturbedExclusionsFound);
    }
}

void ExclusionChecker::scheduleCheckOfExclusions(const int numPerturbedExclusionsToReduce)
{
    // When we have a single domain, we don't need to reduce and we algorithmically can not miss
    // any interactions, so we can assert here.
    if (!havePPDomainDecomposition(impl_->cr_))
    {
        impl_->check(numPerturbedExclusionsToReduce);
    }
    else
    {
        // Fill the reduction buffer with the value from this domain to reduce
        impl_->reductionBuffer_[0] = double(numPerturbedExclusionsToReduce);

        // Pass the post-reduction callback to the ObservablesReducer via
        // the callback it gave us for the purpose.
        //
        // Note that it's possible that the callbackAfterReduction is already
        // outstanding, e.g. if repartitioning was triggered before
        // observables were reduced. This could happen for example when
        // replica exchange took place soon after a partition. If so, the
        // callback will be called again. So long as there is no race
        // between the calls to this function and the calls to
        // ObservablesReducer for reduction, this will work correctly. It
        // could be made safer e.g. with checks against duplicate
        // callbacks, but there is no problem to solve.
        //
        // There is no need to check the return value from this callback,
        // as it is not an error to request reduction at a future step.
        impl_->callbackToRequireReduction_(gmx::ReductionRequirement::Eventually);
    }
}
