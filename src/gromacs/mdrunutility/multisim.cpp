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
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_mdrunutility
 */
#include "gmxpre.h"

#include "multisim.h"

#include "config.h"

#include <cinttypes>

#include <filesystem>
#include <limits>

#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/mpitypes.h"
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

//! Sums array \p r of size \p nr over \p mpi_comm, the result is available on all ranks
template<typename T>
static void mpiAllSumReduce(const int nr, T r[], MPI_Comm mpi_comm)
{
#if GMX_MPI
    MPI_Allreduce(MPI_IN_PLACE, r, nr, gmx::mpiType<T>(), MPI_SUM, mpi_comm);
#else
    GMX_RELEASE_ASSERT(false, "Invalid call to mpiAllReduce()");

    GMX_UNUSED_VALUE(nr);
    GMX_UNUSED_VALUE(r);
    GMX_UNUSED_VALUE(mpi_comm);
#endif
}

//! Sums array \p r over \p mpi_comm, the result is available on all ranks, size of \p r should not exceed MAX_INT
template<typename T>
static void mpiAllSumReduce(gmx::ArrayRef<T> r, MPI_Comm mpi_comm)
{
    GMX_RELEASE_ASSERT(r.size() <= std::numeric_limits<int>::max(), "Size should fit in an int");

    mpiAllSumReduce(r.size(), r.data(), mpi_comm);
}

void gmx_sumd_sim(int nr, double r[], const gmx_multisim_t* ms)
{
    mpiAllSumReduce(nr, r, ms->mainRanksComm_);
}

void gmx_sumf_sim(int nr, float r[], const gmx_multisim_t* ms)
{
    mpiAllSumReduce(nr, r, ms->mainRanksComm_);
}

void gmx_sumi_sim(int nr, int r[], const gmx_multisim_t* ms)
{
    mpiAllSumReduce(nr, r, ms->mainRanksComm_);
}

void gmx_sumli_sim(gmx::ArrayRef<int64_t> r, const gmx_multisim_t& ms)
{
    mpiAllSumReduce(r, ms.mainRanksComm_);
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

/*! \brief Returns whether \p value is equal on all main ranks in \p ms
 *
 * Optionally an external buffer can be passed in \p externaBuffer. This will then hold
 * values of value of each main rank in \p ms.
 */
template<typename T>
static bool multiSimValuesAllAreEqual(const gmx_multisim_t& ms,
                                      const T               value,
                                      std::vector<T>*       externalBuffer = nullptr)
{
    std::vector<T>  localBuffer;
    std::vector<T>& buffer = externalBuffer ? *externalBuffer : localBuffer;

    buffer.resize(ms.numSimulations_, 0);
    /* send our value to all other main ranks, receive all of theirs */
    buffer[ms.simulationIndex_] = value;
    mpiAllSumReduce<T>(buffer, ms.mainRanksComm_);

    bool allValuesAreEqual = true;
    for (int s = 0; s < ms.numSimulations_; s++)
    {
        if (buffer[s] != value)
        {
            allValuesAreEqual = false;
            break;
        }
    }

    return allValuesAreEqual;
}

/*! \brief Checks whether values \p val are equal on all main ranks of \p ms
 *
 * \param[in] log   Pointer to the loga file, can be nullptr
 * \param[in] ms    Multi-simulation communication setup
 * \param[in] val   The value to check
 * \param[in] name  The name to print to log, can be nullptr when \p log==nullptr
 */
template<typename T>
static void checkMultiInt(FILE* log, const gmx_multisim_t& ms, const T val, const char* name)
{
    if (nullptr != log)
    {
        fprintf(log, "Multi-checking %s ... ", name);
    }

    std::vector<T> buf;
    if (multiSimValuesAllAreEqual(ms, val, &buf))
    {
        if (nullptr != log)
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
            for (int p = 0; p < ms.numSimulations_; p++)
            {
                char strbuf[255];
                /* first make the format string */
                snprintf(strbuf, 255, "  subsystem %%d: %s\n", std::is_same_v<T, int> ? "%d" : "%" PRId64);
                fprintf(log, strbuf, p, buf[p]);
            }
        }
        gmx_fatal(FARGS, "The %d subsystems are not compatible\n", ms.numSimulations_);
    }
}

void check_multi_int(FILE* log, const gmx_multisim_t& ms, int val, const char* name, bool beQuiet)
{
    checkMultiInt(beQuiet ? nullptr : log, ms, val, name);
}

void check_multi_int64(FILE* log, const gmx_multisim_t& ms, int64_t val, const char* name, bool beQuiet)
{
    checkMultiInt(beQuiet ? nullptr : log, ms, val, name);
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

void logInitialMultisimStatus(const gmx_multisim_t& ms,
                              const t_commrec*      cr,
                              const gmx::MDLogger&  mdlog,
                              const bool            simulationsShareState,
                              const int             numSteps,
                              const int             initialStep)
{
    if (!multiSimValuesAllAreEqual(ms, numSteps))
    {
        GMX_LOG(mdlog.warning)
                .appendText(
                        "Note: The number of steps is not consistent across multi "
                        "simulations,\n"
                        "but we are proceeding anyway!");
    }
    if (!multiSimValuesAllAreEqual(ms, initialStep))
    {
        if (simulationsShareState)
        {
            if (cr->commMySim.isMainRank())
            {
                gmx_fatal(FARGS,
                          "The initial step is not consistent across multi simulations which "
                          "share the state");
            }
            gmx_barrier(cr->commMyGroup.comm());
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
