/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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

#include "energyhistory.h"

#include <cstddef>

#include <string>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

#include "checkpointdata.h"

//! \cond INTERNAL
// mirroring the \cond from energyhistory.h to avoid Doxygen errors

namespace
{
/*!
 * \brief Enum describing the contents delta_h_history_t writes to modular checkpoint
 *
 * When changing the checkpoint content, add a new element just above Count, and adjust the
 * checkpoint functionality.
 */
enum class DeltaHHistoryCheckpointVersion
{
    Base, //!< First version of modular checkpointing
    Count //!< Number of entries. Add new versions right above this!
};
constexpr auto c_currentVersionDeltaHH =
        DeltaHHistoryCheckpointVersion(int(DeltaHHistoryCheckpointVersion::Count) - 1);
} // namespace

//! Read / write vector size from / to checkpoint and resize vector if reading
template<gmx::CheckpointDataOperation operation, typename T>
static void checkpointVectorSize(gmx::CheckpointData<operation>* checkpointData,
                                 const std::string&              name,
                                 std::vector<T>*                 vector)
{
    auto size = static_cast<int64_t>(vector->size());
    checkpointData->scalar(name, &size);
    if (operation == gmx::CheckpointDataOperation::Read)
    {
        vector->resize(size);
    }
};

template<gmx::CheckpointDataOperation operation>
void delta_h_history_t::doCheckpoint(gmx::CheckpointData<operation> checkpointData)
{
    gmx::checkpointVersion(&checkpointData, "delta_h_history_t version", c_currentVersionDeltaHH);

    checkpointVectorSize(&checkpointData, "numDeltaH", &dh);
    checkpointData.scalar("start_time", &start_time);
    checkpointData.scalar("start_lambda", &start_lambda);
    checkpointData.scalar("start_lambda_set", &start_lambda_set);
    for (std::size_t idx = 0; idx < dh.size(); ++idx)
    {
        checkpointVectorSize(&checkpointData, gmx::formatString("vecSize %zu", idx), &dh[idx]);
        checkpointData.arrayRef(gmx::formatString("vec %zu", idx),
                                gmx::makeCheckpointArrayRef<operation>(dh[idx]));
    }
}

namespace
{
/*!
 * \brief Enum describing the contents energyhistory_t writes to modular checkpoint
 *
 * When changing the checkpoint content, add a new element just above Count, and adjust the
 * checkpoint functionality.
 */
enum class EnergyHistoryCheckpointVersion
{
    Base, //!< First version of modular checkpointing
    Count //!< Number of entries. Add new versions right above this!
};
constexpr auto c_currentVersionEnergyHistory =
        EnergyHistoryCheckpointVersion(int(EnergyHistoryCheckpointVersion::Count) - 1);
} // namespace

template<gmx::CheckpointDataOperation operation>
void energyhistory_t::doCheckpoint(gmx::CheckpointData<operation> checkpointData)
{
    gmx::checkpointVersion(&checkpointData, "energyhistory_t version", c_currentVersionEnergyHistory);

    bool useCheckpoint = (nsum <= 0 && nsum_sim <= 0);
    checkpointData.scalar("useCheckpoint", &useCheckpoint);

    if (!useCheckpoint)
    {
        return;
    }

    checkpointVectorSize(&checkpointData, "enerAveSize", &ener_ave);
    checkpointVectorSize(&checkpointData, "enerSumSize", &ener_sum);
    checkpointVectorSize(&checkpointData, "enerSumSimSize", &ener_sum_sim);

    checkpointData.scalar("nsteps", &nsteps);
    checkpointData.scalar("nsteps_sim", &nsteps_sim);

    checkpointData.scalar("nsum", &nsum);
    checkpointData.scalar("nsum_sim", &nsum_sim);

    auto hasForeignLambdas = (deltaHForeignLambdas != nullptr);
    checkpointData.scalar("has foreign lambdas", &hasForeignLambdas);
    if (hasForeignLambdas && deltaHForeignLambdas == nullptr)
    {
        deltaHForeignLambdas = std::make_unique<delta_h_history_t>();
    }

    if (nsum > 0)
    {
        checkpointData.arrayRef("ener_ave", gmx::makeCheckpointArrayRef<operation>(ener_ave));
        checkpointData.arrayRef("ener_sum", gmx::makeCheckpointArrayRef<operation>(ener_sum));
    }
    if (nsum_sim > 0)
    {
        checkpointData.arrayRef("ener_sum_sim", gmx::makeCheckpointArrayRef<operation>(ener_sum_sim));
    }
    if (hasForeignLambdas)
    {
        deltaHForeignLambdas->doCheckpoint<operation>(
                checkpointData.subCheckpointData("deltaHForeignLambdas"));
    }
}

// explicit template instatiation
template void energyhistory_t::doCheckpoint(gmx::CheckpointData<gmx::CheckpointDataOperation::Read> checkpointData);
template void energyhistory_t::doCheckpoint(gmx::CheckpointData<gmx::CheckpointDataOperation::Write> checkpointData);


//! \endcond
