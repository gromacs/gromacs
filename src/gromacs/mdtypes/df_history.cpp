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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "df_history.h"

#include "gromacs/utility/arrayref.h"

df_history_t::df_history_t(const int numLambdaValues)
{
    nlambda = numLambdaValues;

    sum_weights.resize(numLambdaValues);
    sum_dg.resize(numLambdaValues);
    sum_minvar.resize(numLambdaValues);
    sum_variance.resize(numLambdaValues);
    numSamplesAtLambdaForStatistics.resize(numLambdaValues);
    numSamplesAtLambdaForEquilibration.resize(numLambdaValues);
    wl_histo.resize(numLambdaValues);

    // Allocate transition matrices here
    Tij.resize(numLambdaValues, numLambdaValues);
    Tij_empirical.resize(numLambdaValues, numLambdaValues);

    // Allocate accumulators for various transition-matrix
    // free energy methods here.
    accum_p.resize(numLambdaValues, numLambdaValues);
    accum_m.resize(numLambdaValues, numLambdaValues);
    accum_p2.resize(numLambdaValues, numLambdaValues);
    accum_m2.resize(numLambdaValues, numLambdaValues);
}

namespace
{

/*!
 * \brief Enum describing the contents df_history_t writes to modular checkpoint
 *
 * When changing the checkpoint content, add a new element just above Count, and adjust the
 * checkpoint functionality.
 */
enum class DFHistoryCheckpointVersion
{
    Base, //!< First version of modular checkpointing
    Count //!< Number of entries. Add new versions right above this!
};
constexpr auto c_dfHistoryCurrentVersion =
        DFHistoryCheckpointVersion(int(DFHistoryCheckpointVersion::Count) - 1);
} // namespace

template<gmx::CheckpointDataOperation operation>
void df_history_t::doCheckpoint(gmx::CheckpointData<operation> checkpointData, LambdaWeightCalculation elamstats)
{
    gmx::checkpointVersion(&checkpointData, "df_history_t version", c_dfHistoryCurrentVersion);

    auto numLambdas = nlambda;
    checkpointData.scalar("nlambda", &numLambdas);
    if (operation == gmx::CheckpointDataOperation::Read)
    {
        // If this isn't matching, we haven't allocated the right amount of data
        GMX_RELEASE_ASSERT(numLambdas == nlambda,
                           "df_history_t checkpoint reading: Lambda vectors size mismatch.");
    }

    checkpointData.scalar("bEquil", &bEquil);
    checkpointData.arrayRef("numSamplesAtLambdaForStatistics",
                            gmx::makeCheckpointArrayRef<operation>(numSamplesAtLambdaForStatistics));
    checkpointData.arrayRef("numSamplesAtLambdaForEquilibration",
                            gmx::makeCheckpointArrayRef<operation>(numSamplesAtLambdaForEquilibration));
    checkpointData.arrayRef("sum_weights", gmx::makeCheckpointArrayRef<operation>(sum_weights));
    checkpointData.arrayRef("sum_dg", gmx::makeCheckpointArrayRef<operation>(sum_dg));
    checkpointData.arrayRef("Tij", gmx::makeCheckpointArrayRef<operation>(Tij.toArrayRef()));
    checkpointData.arrayRef("Tij_empirical",
                            gmx::makeCheckpointArrayRef<operation>(Tij_empirical.toArrayRef()));

    if (EWL(elamstats))
    {
        checkpointData.arrayRef("wl_histo", gmx::makeCheckpointArrayRef<operation>(wl_histo));
        checkpointData.scalar("wl_delta", &wl_delta);
    }
    if ((elamstats == LambdaWeightCalculation::Minvar) || (elamstats == LambdaWeightCalculation::Barker)
        || (elamstats == LambdaWeightCalculation::Metropolis))
    {
        checkpointData.arrayRef("sum_minvar", gmx::makeCheckpointArrayRef<operation>(sum_minvar));
        checkpointData.arrayRef("sum_variance", gmx::makeCheckpointArrayRef<operation>(sum_variance));
        checkpointData.arrayRef("accum_p", gmx::makeCheckpointArrayRef<operation>(accum_p.toArrayRef()));
        checkpointData.arrayRef("accum_m", gmx::makeCheckpointArrayRef<operation>(accum_m.toArrayRef()));
        checkpointData.arrayRef("accum_p2", gmx::makeCheckpointArrayRef<operation>(accum_p2.toArrayRef()));
        checkpointData.arrayRef("accum_m2", gmx::makeCheckpointArrayRef<operation>(accum_m2.toArrayRef()));
    }
}

// extern template specializations
template void df_history_t::doCheckpoint(gmx::CheckpointData<gmx::CheckpointDataOperation::Read> checkpointData,
                                         LambdaWeightCalculation elamstats);
template void df_history_t::doCheckpoint(gmx::CheckpointData<gmx::CheckpointDataOperation::Write> checkpointData,
                                         LambdaWeightCalculation elamstats);
