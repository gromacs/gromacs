/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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

/*! \internal \file
 *
 * \brief
 * Declares the PairlistSets class
 *
 * This class holds the local and non-local pairlist sets.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_PAIRLISTSETS_H
#define GMX_NBNXM_PAIRLISTSETS_H

#include <memory>

#include "gromacs/mdtypes/locality.h"

#include "pairlistparams.h"

struct nbnxn_atomdata_t;
class PairlistSet;
enum class PairlistType;
class PairSearch;
struct t_blocka;
struct t_nrnb;


class PairlistSets
{
public:
    PairlistSets(const PairlistParams& pairlistParams,
                 bool                  haveMultipleDomains,
                 int                   minimumIlistCountForGpuBalancing);

    //! Construct the pairlist set for the given locality
    void construct(gmx::InteractionLocality iLocality,
                   PairSearch*              pairSearch,
                   nbnxn_atomdata_t*        nbat,
                   const t_blocka*          excl,
                   int64_t                  step,
                   t_nrnb*                  nrnb);

    //! Dispatches the dynamic pruning kernel for the given locality
    void dispatchPruneKernel(gmx::InteractionLocality iLocality,
                             const nbnxn_atomdata_t*  nbat,
                             const rvec*              shift_vec);

    //! Returns the pair list parameters
    const PairlistParams& params() const { return params_; }

    //! Returns the number of steps performed with the current pair list
    int numStepsWithPairlist(int64_t step) const
    {
        return static_cast<int>(step - outerListCreationStep_);
    }

    //! Returns whether step is a dynamic list pruning step, for CPU lists
    bool isDynamicPruningStepCpu(int64_t step) const
    {
        return (params_.useDynamicPruning && numStepsWithPairlist(step) % params_.nstlistPrune == 0);
    }

    //! Returns whether step is a dynamic list pruning step, for GPU lists
    bool isDynamicPruningStepGpu(int64_t step) const
    {
        const int age = numStepsWithPairlist(step);

        return (params_.useDynamicPruning && age > 0 && age < params_.lifetime
                && (params_.haveMultipleDomains || age % 2 == 0));
    }

    //! Changes the pair-list outer and inner radius
    void changePairlistRadii(real rlistOuter, real rlistInner)
    {
        params_.rlistOuter = rlistOuter;
        params_.rlistInner = rlistInner;
    }

    //! Returns the pair-list set for the given locality
    const PairlistSet& pairlistSet(gmx::InteractionLocality iLocality) const
    {
        if (iLocality == gmx::InteractionLocality::Local)
        {
            return *localSet_;
        }
        else
        {
            GMX_ASSERT(nonlocalSet_, "Need a non-local set when requesting access");
            return *nonlocalSet_;
        }
    }

private:
    //! Returns the pair-list set for the given locality
    PairlistSet& pairlistSet(gmx::InteractionLocality iLocality)
    {
        if (iLocality == gmx::InteractionLocality::Local)
        {
            return *localSet_;
        }
        else
        {
            GMX_ASSERT(nonlocalSet_, "Need a non-local set when requesting access");
            return *nonlocalSet_;
        }
    }

    //! Parameters for the search and list pruning setup
    PairlistParams params_;
    //! Pair list balancing parameter for use with GPU
    int minimumIlistCountForGpuBalancing_;
    //! Local pairlist set
    std::unique_ptr<PairlistSet> localSet_;
    //! Non-local pairlist set
    std::unique_ptr<PairlistSet> nonlocalSet_;
    //! MD step at with the outer lists in pairlistSets_ were created
    int64_t outerListCreationStep_;
};

#endif
