/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
 * \brief Implements the multi-simulation support routines.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrunutility
 */
#include "gmxpre.h"

#include "multisim.h"

#include "config.h"

#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/mpiinplacebuffers.h"
#include "gromacs/utility/smalloc.h"

gmx_multisim_t::gmx_multisim_t() = default;

gmx_multisim_t::gmx_multisim_t(MPI_Comm comm, gmx::ArrayRef<const std::string> multidirs)
{
    if (multidirs.empty())
    {
        return;
    }

    if (!GMX_LIB_MPI && !multidirs.empty())
    {
        gmx_fatal(FARGS,
                  "mdrun -multidir is only supported when GROMACS has been "
                  "configured with a proper external MPI library.");
    }

    if (multidirs.size() == 1)
    {
        /* NOTE: It would be nice if this special case worked, but this requires checks/tests. */
        gmx_fatal(FARGS,
                  "To run mdrun in multiple simulation mode, more then one "
                  "actual simulation is required. The single simulation case is not supported.");
    }

#if GMX_MPI
    int numRanks;
    MPI_Comm_size(comm, &numRanks);
    if (numRanks % multidirs.size() != 0)
    {
        gmx_fatal(FARGS,
                  "The number of ranks (%d) is not a multiple of the number of simulations (%td)",
                  numRanks, multidirs.ssize());
    }

    int numRanksPerSim = numRanks / multidirs.size();
    int rankWithinComm;
    MPI_Comm_rank(comm, &rankWithinComm);

    if (debug)
    {
        fprintf(debug, "We have %td simulations, %d ranks per simulation, local simulation is %d\n",
                multidirs.ssize(), numRanksPerSim, rankWithinComm / numRanksPerSim);
    }

    nsim = multidirs.size();
    sim  = rankWithinComm / numRanksPerSim;
    /* Create a communicator for the master nodes */
    std::vector<int> rank(nsim);
    for (int i = 0; i < nsim; i++)
    {
        rank[i] = i * numRanksPerSim;
    }
    MPI_Group mpi_group_world;
    MPI_Comm_group(comm, &mpi_group_world);
    MPI_Group_incl(mpi_group_world, nsim, rank.data(), &mpi_group_masters);
    MPI_Comm_create(comm, mpi_group_masters, &mpi_comm_masters);

#    if !MPI_IN_PLACE_EXISTS
    /* initialize the MPI_IN_PLACE replacement buffers */
    snew(mpb, 1);
    mpb->ibuf        = nullptr;
    mpb->libuf       = nullptr;
    mpb->fbuf        = nullptr;
    mpb->dbuf        = nullptr;
    mpb->ibuf_alloc  = 0;
    mpb->libuf_alloc = 0;
    mpb->fbuf_alloc  = 0;
    mpb->dbuf_alloc  = 0;
#    endif

    // TODO This should throw upon error
    gmx_chdir(multidirs[sim].c_str());
#else
    GMX_UNUSED_VALUE(comm);
#endif
}

gmx_multisim_t::~gmx_multisim_t()
{
    done_mpi_in_place_buf(mpb);

#if GMX_MPI
    // TODO This would work better if the result of MPI_Comm_split was
    // put into an RAII-style guard, such as gmx::unique_cptr.
    if (mpi_comm_masters != MPI_COMM_NULL && mpi_comm_masters != MPI_COMM_WORLD)
    {
        MPI_Comm_free(&mpi_comm_masters);
    }
    if (mpi_group_masters != MPI_GROUP_NULL)
    {
        MPI_Group_free(&mpi_group_masters);
    }
#endif
}

#if GMX_MPI
static void gmx_sumd_comm(int nr, double r[], MPI_Comm mpi_comm)
{
#    if MPI_IN_PLACE_EXISTS
    MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_DOUBLE, MPI_SUM, mpi_comm);
#    else
    /* this function is only used in code that is not performance critical,
       (during setup, when comm_rec is not the appropriate communication
       structure), so this isn't as bad as it looks. */
    double* buf;
    int     i;

    snew(buf, nr);
    MPI_Allreduce(r, buf, nr, MPI_DOUBLE, MPI_SUM, mpi_comm);
    for (i = 0; i < nr; i++)
    {
        r[i] = buf[i];
    }
    sfree(buf);
#    endif
}
#endif

#if GMX_MPI
static void gmx_sumf_comm(int nr, float r[], MPI_Comm mpi_comm)
{
#    if MPI_IN_PLACE_EXISTS
    MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_FLOAT, MPI_SUM, mpi_comm);
#    else
    /* this function is only used in code that is not performance critical,
       (during setup, when comm_rec is not the appropriate communication
       structure), so this isn't as bad as it looks. */
    float* buf;
    int    i;

    snew(buf, nr);
    MPI_Allreduce(r, buf, nr, MPI_FLOAT, MPI_SUM, mpi_comm);
    for (i = 0; i < nr; i++)
    {
        r[i] = buf[i];
    }
    sfree(buf);
#    endif
}
#endif

void gmx_sumd_sim(int gmx_unused nr, double gmx_unused r[], const gmx_multisim_t gmx_unused* ms)
{
#if !GMX_MPI
    gmx_call("gmx_sumd_sim");
#else
    gmx_sumd_comm(nr, r, ms->mpi_comm_masters);
#endif
}

void gmx_sumf_sim(int gmx_unused nr, float gmx_unused r[], const gmx_multisim_t gmx_unused* ms)
{
#if !GMX_MPI
    gmx_call("gmx_sumf_sim");
#else
    gmx_sumf_comm(nr, r, ms->mpi_comm_masters);
#endif
}

void gmx_sumi_sim(int gmx_unused nr, int gmx_unused r[], const gmx_multisim_t gmx_unused* ms)
{
#if !GMX_MPI
    gmx_call("gmx_sumi_sim");
#else
#    if MPI_IN_PLACE_EXISTS
    MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_INT, MPI_SUM, ms->mpi_comm_masters);
#    else
    /* this is thread-unsafe, but it will do for now: */
    int i;

    if (nr > ms->mpb->ibuf_alloc)
    {
        ms->mpb->ibuf_alloc = nr;
        srenew(ms->mpb->ibuf, ms->mpb->ibuf_alloc);
    }
    MPI_Allreduce(r, ms->mpb->ibuf, nr, MPI_INT, MPI_SUM, ms->mpi_comm_masters);
    for (i = 0; i < nr; i++)
    {
        r[i] = ms->mpb->ibuf[i];
    }
#    endif
#endif
}

void gmx_sumli_sim(int gmx_unused nr, int64_t gmx_unused r[], const gmx_multisim_t gmx_unused* ms)
{
#if !GMX_MPI
    gmx_call("gmx_sumli_sim");
#else
#    if MPI_IN_PLACE_EXISTS
    MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_INT64_T, MPI_SUM, ms->mpi_comm_masters);
#    else
    /* this is thread-unsafe, but it will do for now: */
    int i;

    if (nr > ms->mpb->libuf_alloc)
    {
        ms->mpb->libuf_alloc = nr;
        srenew(ms->mpb->libuf, ms->mpb->libuf_alloc);
    }
    MPI_Allreduce(r, ms->mpb->libuf, nr, MPI_INT64_T, MPI_SUM, ms->mpi_comm_masters);
    for (i = 0; i < nr; i++)
    {
        r[i] = ms->mpb->libuf[i];
    }
#    endif
#endif
}

std::vector<int> gatherIntFromMultiSimulation(const gmx_multisim_t* ms, const int localValue)
{
    std::vector<int> valuesFromAllRanks;
    if (GMX_MPI && ms != nullptr)
    {
        valuesFromAllRanks.resize(ms->nsim);
        valuesFromAllRanks[ms->sim] = localValue;
        gmx_sumi_sim(ms->nsim, valuesFromAllRanks.data(), ms);
    }
    else
    {
        valuesFromAllRanks.emplace_back(localValue);
    }
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

    snew(ibuf, ms->nsim);
    ibuf[ms->sim] = val;
    gmx_sumi_sim(ms->nsim, ibuf, ms);

    bCompatible = TRUE;
    for (p = 1; p < ms->nsim; p++)
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
            for (p = 0; p < ms->nsim; p++)
            {
                fprintf(log, "  subsystem %d: %d\n", p, ibuf[p]);
            }
        }
        gmx_fatal(FARGS, "The %d subsystems are not compatible\n", ms->nsim);
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

    snew(ibuf, ms->nsim);
    ibuf[ms->sim] = val;
    gmx_sumli_sim(ms->nsim, ibuf, ms);

    bCompatible = TRUE;
    for (p = 1; p < ms->nsim; p++)
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
            for (p = 0; p < ms->nsim; p++)
            {
                char strbuf[255];
                /* first make the format string */
                snprintf(strbuf, 255, "  subsystem %%d: %s\n", "%" PRId64);
                fprintf(log, strbuf, p, ibuf[p]);
            }
        }
        gmx_fatal(FARGS, "The %d subsystems are not compatible\n", ms->nsim);
    }

    sfree(ibuf);
}

bool findIsSimulationMasterRank(const gmx_multisim_t* ms, MPI_Comm communicator)
{
    if (GMX_LIB_MPI)
    {
        // Ranks of multi-simulations know whether they are a master
        // rank. Ranks of non-multi simulation do not know until a
        // t_commrec is available.
        if ((ms != nullptr) && (ms->nsim > 1))
        {
            return ms->mpi_comm_masters != MPI_COMM_NULL;
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
        // ie is the master thread
        return (communicator == MPI_COMM_NULL);
    }
    else
    {
        // No MPI means it must be the master (and only) rank.
        return true;
    }
}

bool isMasterSim(const gmx_multisim_t* ms)
{
    return !isMultiSim(ms) || ms->sim == 0;
}

bool isMasterSimMasterRank(const gmx_multisim_t* ms, const bool isMaster)
{
    return (isMaster && isMasterSim(ms));
}
