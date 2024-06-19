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
/*! \internal \file
 *
 * \brief Implements the multi-simulation support routines.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrunutility
 */
#include "gmxpre.h"

#include "multisim.h"

#include "config.h"

#include <cinttypes>

#include <filesystem>

#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/smalloc.h"

std::unique_ptr<gmx_multisim_t> buildMultiSimulation(MPI_Comm                         worldComm,
                                                     gmx::ArrayRef<const std::string> multidirs)
{
    if (multidirs.empty())
    {
        return nullptr;
    }

    if (!GMX_LIB_MPI && !multidirs.empty())
    {
        GMX_THROW(gmx::NotImplementedError(
                "Multi-simulations are only supported when GROMACS has been "
                "configured with a proper external MPI library."));
    }

    if (multidirs.size() == 1)
    {
        /* NOTE: It would be nice if this special case worked, but this requires checks/tests. */
        GMX_THROW(gmx::NotImplementedError(
                "To run mdrun in multi-simulation mode, more then one "
                "actual simulation is required. The single simulation case is not supported."));
    }

#if GMX_LIB_MPI
    int numRanks;
    MPI_Comm_size(worldComm, &numRanks);
    if (numRanks % multidirs.size() != 0)
    {
        auto message = gmx::formatString(
                "The number of ranks (%d) is not a multiple of the number of simulations (%td)",
                numRanks,
                multidirs.ssize());
        GMX_THROW(gmx::InconsistentInputError(message));
    }

    int numRanksPerSimulation = numRanks / multidirs.size();
    int rankWithinWorldComm;
    MPI_Comm_rank(worldComm, &rankWithinWorldComm);

    if (debug)
    {
        fprintf(debug,
                "We have %td simulations, %d ranks per simulation, local simulation is %d\n",
                multidirs.ssize(),
                numRanksPerSimulation,
                rankWithinWorldComm / numRanksPerSimulation);
    }

    int numSimulations = multidirs.size();
    // Create a communicator for the main ranks of each simulation
    std::vector<int> ranksOfMains(numSimulations);
    for (int i = 0; i < numSimulations; i++)
    {
        ranksOfMains[i] = i * numRanksPerSimulation;
    }
    MPI_Group worldGroup;
    // No need to free worldGroup later, we didn't create it.
    MPI_Comm_group(worldComm, &worldGroup);

    MPI_Group mainsGroup = MPI_GROUP_NULL;
    MPI_Group_incl(worldGroup, numSimulations, ranksOfMains.data(), &mainsGroup);
    MPI_Comm mainRanksComm = MPI_COMM_NULL;
    MPI_Comm_create(worldComm, mainsGroup, &mainRanksComm);
    if (mainsGroup != MPI_GROUP_NULL)
    {
        MPI_Group_free(&mainsGroup);
    }

    int      simulationIndex = rankWithinWorldComm / numRanksPerSimulation;
    MPI_Comm simulationComm  = MPI_COMM_NULL;
    MPI_Comm_split(worldComm, simulationIndex, rankWithinWorldComm, &simulationComm);

    try
    {
        gmx_chdir(multidirs[simulationIndex].c_str());
    }
    catch (gmx::GromacsException& e)
    {
        e.prependContext("While changing directory for multi-simulation to " + multidirs[simulationIndex]);
        throw;
    }
    return std::make_unique<gmx_multisim_t>(numSimulations, simulationIndex, mainRanksComm, simulationComm);
#else
    GMX_UNUSED_VALUE(worldComm);
    return nullptr;
#endif
}

gmx_multisim_t::gmx_multisim_t(int numSimulations, int simulationIndex, MPI_Comm mainRanksComm, MPI_Comm simulationComm) :
    numSimulations_(numSimulations),
    simulationIndex_(simulationIndex),
    mainRanksComm_(mainRanksComm),
    simulationComm_(simulationComm)
{
}

gmx_multisim_t::~gmx_multisim_t()
{
#if GMX_LIB_MPI
    // TODO This would work better if the result of MPI_Comm_split was
    // put into an RAII-style guard, such as gmx::unique_cptr.
    if (mainRanksComm_ != MPI_COMM_NULL && mainRanksComm_ != MPI_COMM_WORLD)
    {
        MPI_Comm_free(&mainRanksComm_);
    }
    if (simulationComm_ != MPI_COMM_NULL && simulationComm_ != MPI_COMM_WORLD)
    {
        MPI_Comm_free(&simulationComm_);
    }
#endif
}

#if GMX_MPI
static void gmx_sumd_comm(int nr, double r[], MPI_Comm mpi_comm)
{
    MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_DOUBLE, MPI_SUM, mpi_comm);
}
#endif

#if GMX_MPI
static void gmx_sumf_comm(int nr, float r[], MPI_Comm mpi_comm)
{
    MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_FLOAT, MPI_SUM, mpi_comm);
}
#endif

void gmx_sumd_sim(int gmx_unused nr, double gmx_unused r[], const gmx_multisim_t gmx_unused* ms)
{
#if !GMX_MPI
    GMX_RELEASE_ASSERT(false, "Invalid call to gmx_sumd_sim");
#else
    gmx_sumd_comm(nr, r, ms->mainRanksComm_);
#endif
}

void gmx_sumf_sim(int gmx_unused nr, float gmx_unused r[], const gmx_multisim_t gmx_unused* ms)
{
#if !GMX_MPI
    GMX_RELEASE_ASSERT(false, "Invalid call to gmx_sumf_sim");
#else
    gmx_sumf_comm(nr, r, ms->mainRanksComm_);
#endif
}

void gmx_sumi_sim(int gmx_unused nr, int gmx_unused r[], const gmx_multisim_t gmx_unused* ms)
{
#if !GMX_MPI
    GMX_RELEASE_ASSERT(false, "Invalid call to gmx_sumi_sim");
#else
    MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_INT, MPI_SUM, ms->mainRanksComm_);
#endif
}

void gmx_sumli_sim(int gmx_unused nr, int64_t gmx_unused r[], const gmx_multisim_t gmx_unused* ms)
{
#if !GMX_MPI
    GMX_RELEASE_ASSERT(false, "Invalid call to gmx_sumli_sim");
#else
    MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_INT64_T, MPI_SUM, ms->mainRanksComm_);
#endif
}

std::vector<int> gatherIntFromMultiSimulation(const gmx_multisim_t* ms, const int localValue)
{
    std::vector<int> valuesFromAllRanks;
#if GMX_MPI
    if (ms != nullptr)
    {
        valuesFromAllRanks.resize(ms->numSimulations_);
        valuesFromAllRanks[ms->simulationIndex_] = localValue;
        gmx_sumi_sim(ms->numSimulations_, valuesFromAllRanks.data(), ms);
    }
    else
    {
        valuesFromAllRanks.emplace_back(localValue);
    }
#else
    GMX_UNUSED_VALUE(ms);
    valuesFromAllRanks.emplace_back(localValue);
#endif
    return valuesFromAllRanks;
}

void check_multi_int(FILE* log, const gmx_multisim_t* ms, int val, const char* name, gmx_bool bQuiet)
{
    int *    ibuf, p;
    gmx_bool bCompatible;

    if (nullptr != log && !bQuiet)
    {
        fprintf(log, "Multi-checking %s ... ", name);
    }

    if (ms == nullptr)
    {
        gmx_fatal(FARGS, "check_multi_int called with a NULL communication pointer");
    }

    snew(ibuf, ms->numSimulations_);
    ibuf[ms->simulationIndex_] = val;
    gmx_sumi_sim(ms->numSimulations_, ibuf, ms);

    bCompatible = TRUE;
    for (p = 1; p < ms->numSimulations_; p++)
    {
        bCompatible = bCompatible && (ibuf[p - 1] == ibuf[p]);
    }

    if (bCompatible)
    {
        if (nullptr != log && !bQuiet)
        {
            fprintf(log, "OK\n");
        }
    }
    else
    {
        if (nullptr != log)
        {
            fprintf(log, "\n%s is not equal for all subsystems\n", name);
            for (p = 0; p < ms->numSimulations_; p++)
            {
                fprintf(log, "  subsystem %d: %d\n", p, ibuf[p]);
            }
        }
        gmx_fatal(FARGS, "The %d subsystems are not compatible\n", ms->numSimulations_);
    }

    sfree(ibuf);
}

void check_multi_int64(FILE* log, const gmx_multisim_t* ms, int64_t val, const char* name, gmx_bool bQuiet)
{
    int64_t* ibuf;
    int      p;
    gmx_bool bCompatible;

    if (nullptr != log && !bQuiet)
    {
        fprintf(log, "Multi-checking %s ... ", name);
    }

    if (ms == nullptr)
    {
        gmx_fatal(FARGS, "check_multi_int called with a NULL communication pointer");
    }

    snew(ibuf, ms->numSimulations_);
    ibuf[ms->simulationIndex_] = val;
    gmx_sumli_sim(ms->numSimulations_, ibuf, ms);

    bCompatible = TRUE;
    for (p = 1; p < ms->numSimulations_; p++)
    {
        bCompatible = bCompatible && (ibuf[p - 1] == ibuf[p]);
    }

    if (bCompatible)
    {
        if (nullptr != log && !bQuiet)
        {
            fprintf(log, "OK\n");
        }
    }
    else
    {
        // TODO Part of this error message would also be good to go to
        // stderr (from one rank of one sim only)
        if (nullptr != log)
        {
            fprintf(log, "\n%s is not equal for all subsystems\n", name);
            for (p = 0; p < ms->numSimulations_; p++)
            {
                char strbuf[255];
                /* first make the format string */
                snprintf(strbuf, 255, "  subsystem %%d: %s\n", "%" PRId64);
                fprintf(log, strbuf, p, ibuf[p]);
            }
        }
        gmx_fatal(FARGS, "The %d subsystems are not compatible\n", ms->numSimulations_);
    }

    sfree(ibuf);
}

bool findIsSimulationMainRank(const gmx_multisim_t* ms, MPI_Comm communicator)
{
    if (GMX_LIB_MPI)
    {
        // Ranks of multi-simulations know whether they are a main
        // rank. Ranks of non-multi simulation do not know until a
        // t_commrec is available.
        if ((ms != nullptr) && (ms->numSimulations_ > 1))
        {
            return ms->mainRanksComm_ != MPI_COMM_NULL;
        }
        else
        {
            int rank = 0;
#if GMX_LIB_MPI
            MPI_Comm_rank(communicator, &rank);
#endif
            return (rank == 0);
        }
    }
    else if (GMX_THREAD_MPI)
    {
        GMX_RELEASE_ASSERT(communicator == MPI_COMM_NULL || communicator == MPI_COMM_WORLD,
                           "Invalid communicator");
        // Spawned threads have MPI_COMM_WORLD upon creation, so if
        // the communicator is MPI_COMM_NULL this is not a spawned thread,
        // ie is the main thread
        return (communicator == MPI_COMM_NULL);
    }
    else
    {
        // No MPI means it must be the main (and only) rank.
        return true;
    }
}

bool isMainSim(const gmx_multisim_t* ms)
{
    return !isMultiSim(ms) || ms->simulationIndex_ == 0;
}

bool isMainSimMainRank(const gmx_multisim_t* ms, const bool isMain)
{
    return (isMain && isMainSim(ms));
}

static bool multisim_int_all_are_equal(const gmx_multisim_t* ms, int64_t value)
{
    bool     allValuesAreEqual = true;
    int64_t* buf;

    GMX_RELEASE_ASSERT(ms, "Invalid use of multi-simulation pointer");

    snew(buf, ms->numSimulations_);
    /* send our value to all other main ranks, receive all of theirs */
    buf[ms->simulationIndex_] = value;
    gmx_sumli_sim(ms->numSimulations_, buf, ms);

    for (int s = 0; s < ms->numSimulations_; s++)
    {
        if (buf[s] != value)
        {
            allValuesAreEqual = false;
            break;
        }
    }

    sfree(buf);

    return allValuesAreEqual;
}

void logInitialMultisimStatus(const gmx_multisim_t* ms,
                              const t_commrec*      cr,
                              const gmx::MDLogger&  mdlog,
                              const bool            simulationsShareState,
                              const int             numSteps,
                              const int             initialStep)
{
    if (!multisim_int_all_are_equal(ms, numSteps))
    {
        GMX_LOG(mdlog.warning)
                .appendText(
                        "Note: The number of steps is not consistent across multi "
                        "simulations,\n"
                        "but we are proceeding anyway!");
    }
    if (!multisim_int_all_are_equal(ms, initialStep))
    {
        if (simulationsShareState)
        {
            if (MAIN(cr))
            {
                gmx_fatal(FARGS,
                          "The initial step is not consistent across multi simulations which "
                          "share the state");
            }
            gmx_barrier(cr->mpi_comm_mygroup);
        }
        else
        {
            GMX_LOG(mdlog.warning)
                    .appendText(
                            "Note: The initial step is not consistent across multi "
                            "simulations,\n"
                            "but we are proceeding anyway!");
        }
    }
}
