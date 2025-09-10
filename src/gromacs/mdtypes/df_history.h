/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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

#ifndef GMX_MDTYPES_DF_HISTORY_H
#define GMX_MDTYPES_DF_HISTORY_H

#include <vector>

#include "gromacs/math/multidimarray.h"
#include "gromacs/mdtypes/checkpointdata.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/real.h"

namespace gmx
{
enum class CheckpointDataOperation;
template<CheckpointDataOperation operation>
class CheckpointData;
} // namespace gmx

using TwoDLambdaArray =
        gmx::MultiDimArray<std::vector<real>, gmx::extents<gmx::dynamic_extent, gmx::dynamic_extent>>;

/*! \brief Free-energy sampling history struct
 *
 * Note that an intended invariant of this structure is that all
 * vectors (and the dimensions of the matrices) have the same size,
 * ie. nlambda.
 *
 * \todo Split out into microstate and observables history.
 */
struct df_history_t
{
    explicit df_history_t(int numLambdaValues);
    //! Default copy constructor.
    df_history_t(const df_history_t& src) = default;
    //! Default copy assignment operator.
    df_history_t& operator=(const df_history_t& v) = default;
    //! Default move constructor.
    df_history_t(df_history_t&& src) noexcept = default;
    //! Default move assignment operator.
    df_history_t& operator=(df_history_t&& v) noexcept = default;

    int nlambda = 0; //!< total number of lambda states - useful to history as the number of lambdas determines the size of arrays.

    bool bEquil = false; //!< Have we reached equilibration yet, where the weights stop updating?
    std::vector<int> numSamplesAtLambdaForStatistics; //!< The number of points observed at each lambda up to the current time, in this simulation, for calculating statistics
    std::vector<int> numSamplesAtLambdaForEquilibration; //!< The number of points observed at each lambda up to the current time, over a set of simulations, for determining equilibration
    std::vector<real> wl_histo; //!< The histogram for WL flatness determination.  Can be preserved between simulations with input options.
    real wl_delta; //!< The current wang-landau delta, used to increment each state when visited.

    std::vector<real> sum_weights; //!< Sum of weights of each state over all states.
    std::vector<real> sum_dg; //!< Sum of the free energies of the states -- not actually used for weighting, but informational
    std::vector<real> sum_minvar;   //!< corrections to weights for minimum variance
    std::vector<real> sum_variance; //!< variances of the states

    TwoDLambdaArray accum_p;  //!< accumulated bennett weights for n+1
    TwoDLambdaArray accum_m;  //!< accumulated bennett weights for n-1
    TwoDLambdaArray accum_p2; //!< accumulated squared bennett weights for n+1
    TwoDLambdaArray accum_m2; //!< accumulated squared bennett weights for n-1

    TwoDLambdaArray Tij; //!< Transition matrix, estimated from probabilities of transitions.
    TwoDLambdaArray Tij_empirical; //!< Empirical transition matrix, estimated from only counts of transitions.

    /*! \brief Allows to read and write checkpoint within modular simulator
     *
     * \tparam operation  Whether we're reading or writing
     * \param checkpointData  The CheckpointData object
     * \param elamstats  How the lambda weights are calculated
     */
    template<gmx::CheckpointDataOperation operation>
    void doCheckpoint(gmx::CheckpointData<operation> checkpointData, LambdaWeightCalculation elamstats);
};

// extern template declarations, defined in source file
extern template void df_history_t::doCheckpoint(gmx::CheckpointData<gmx::CheckpointDataOperation::Read> checkpointData,
                                                LambdaWeightCalculation elamstats);
extern template void df_history_t::doCheckpoint(gmx::CheckpointData<gmx::CheckpointDataOperation::Write> checkpointData,
                                                LambdaWeightCalculation elamstats);

#endif
