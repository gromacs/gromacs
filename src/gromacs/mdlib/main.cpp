/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "main.h"

#include "config.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/binaryinformation.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/sysinfo.h"

/* The source code in this file should be thread-safe.
         Please keep it that way. */

// TODO move this to multi-sim module
void check_multi_int(FILE *log, const gmx_multisim_t *ms, int val,
                     const char *name,
                     gmx_bool bQuiet)
{
    int     *ibuf, p;
    gmx_bool bCompatible;

    if (nullptr != log && !bQuiet)
    {
        fprintf(log, "Multi-checking %s ... ", name);
    }

    if (ms == nullptr)
    {
        gmx_fatal(FARGS,
                  "check_multi_int called with a NULL communication pointer");
    }

    snew(ibuf, ms->nsim);
    ibuf[ms->sim] = val;
    gmx_sumi_sim(ms->nsim, ibuf, ms);

    bCompatible = TRUE;
    for (p = 1; p < ms->nsim; p++)
    {
        bCompatible = bCompatible && (ibuf[p-1] == ibuf[p]);
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

// TODO move this to multi-sim module
void check_multi_int64(FILE *log, const gmx_multisim_t *ms,
                       gmx_int64_t val, const char *name,
                       gmx_bool bQuiet)
{
    gmx_int64_t      *ibuf;
    int               p;
    gmx_bool          bCompatible;

    if (nullptr != log && !bQuiet)
    {
        fprintf(log, "Multi-checking %s ... ", name);
    }

    if (ms == nullptr)
    {
        gmx_fatal(FARGS,
                  "check_multi_int called with a NULL communication pointer");
    }

    snew(ibuf, ms->nsim);
    ibuf[ms->sim] = val;
    gmx_sumli_sim(ms->nsim, ibuf, ms);

    bCompatible = TRUE;
    for (p = 1; p < ms->nsim; p++)
    {
        bCompatible = bCompatible && (ibuf[p-1] == ibuf[p]);
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
                snprintf(strbuf, 255, "  subsystem %%d: %s\n",
                         "%" GMX_PRId64);
                fprintf(log, strbuf, p, ibuf[p]);
            }
        }
        gmx_fatal(FARGS, "The %d subsystems are not compatible\n", ms->nsim);
    }

    sfree(ibuf);
}


void gmx_log_open(const char *lognm, const t_commrec *cr,
                  gmx_bool bAppendFiles, FILE** fplog)
{
    int    pid;
    char   host[256];
    char   timebuf[STRLEN];
    FILE  *fp = *fplog;

    if (!bAppendFiles)
    {
        fp = gmx_fio_fopen(lognm, bAppendFiles ? "a+" : "w+" );
    }

    gmx_fatal_set_log_file(fp);

    /* Get some machine parameters */
    gmx_gethostname(host, 256);
    pid = gmx_getpid();
    gmx_format_current_time(timebuf, STRLEN);

    if (bAppendFiles)
    {
        fprintf(fp,
                "\n"
                "\n"
                "-----------------------------------------------------------\n"
                "Restarting from checkpoint, appending to previous log file.\n"
                "\n"
                );
    }

    fprintf(fp,
            "Log file opened on %s"
            "Host: %s  pid: %d  rank ID: %d  number of ranks:  %d\n",
            timebuf, host, pid, cr->nodeid, cr->nnodes);
    try
    {
        gmx::BinaryInformationSettings settings;
        settings.extendedInfo(true);
        settings.copyright(!bAppendFiles);
        gmx::printBinaryInformation(fp, gmx::getProgramContext(), settings);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    fprintf(fp, "\n");

    fflush(fp);

    *fplog = fp;
}

void gmx_log_close(FILE *fp)
{
    if (fp)
    {
        gmx_fatal_set_log_file(nullptr);
        gmx_fio_fclose(fp);
    }
}

// TODO move this to multi-sim module
gmx_multisim_t *init_multisystem(MPI_Comm                         comm,
                                 gmx::ArrayRef<const std::string> multidirs)
{
    gmx_multisim_t *ms;
#if GMX_MPI
    MPI_Group       mpi_group_world;
    int            *rank;
#endif

    if (multidirs.empty())
    {
        return nullptr;
    }

    if (!GMX_LIB_MPI && multidirs.size() >= 1)
    {
        gmx_fatal(FARGS, "mdrun -multidir is only supported when GROMACS has been "
                  "configured with a proper external MPI library.");
    }

    if (multidirs.size() == 1)
    {
        /* NOTE: It would be nice if this special case worked, but this requires checks/tests. */
        gmx_fatal(FARGS, "To run mdrun in multiple simulation mode, more then one "
                  "actual simulation is required. The single simulation case is not supported.");
    }

#if GMX_MPI
    int numRanks;
    MPI_Comm_size(comm, &numRanks);
    if (numRanks % multidirs.size() != 0)
    {
        gmx_fatal(FARGS, "The number of ranks (%d) is not a multiple of the number of simulations (%zu)", numRanks, multidirs.size());
    }

    int numRanksPerSim = numRanks/multidirs.size();
    int rankWithinComm;
    MPI_Comm_rank(comm, &rankWithinComm);

    if (debug)
    {
        fprintf(debug, "We have %zu simulations, %d ranks per simulation, local simulation is %d\n", multidirs.size(), numRanksPerSim, rankWithinComm/numRanksPerSim);
    }

    ms       = new gmx_multisim_t;
    ms->nsim = multidirs.size();
    ms->sim  = rankWithinComm/numRanksPerSim;
    /* Create a communicator for the master nodes */
    snew(rank, ms->nsim);
    for (int i = 0; i < ms->nsim; i++)
    {
        rank[i] = i*numRanksPerSim;
    }
    MPI_Comm_group(comm, &mpi_group_world);
    MPI_Group_incl(mpi_group_world, ms->nsim, rank, &ms->mpi_group_masters);
    sfree(rank);
    MPI_Comm_create(MPI_COMM_WORLD, ms->mpi_group_masters,
                    &ms->mpi_comm_masters);

#if !MPI_IN_PLACE_EXISTS
    /* initialize the MPI_IN_PLACE replacement buffers */
    snew(ms->mpb, 1);
    ms->mpb->ibuf        = NULL;
    ms->mpb->libuf       = NULL;
    ms->mpb->fbuf        = NULL;
    ms->mpb->dbuf        = NULL;
    ms->mpb->ibuf_alloc  = 0;
    ms->mpb->libuf_alloc = 0;
    ms->mpb->fbuf_alloc  = 0;
    ms->mpb->dbuf_alloc  = 0;
#endif

    // TODO This should throw upon error
    gmx_chdir(multidirs[ms->sim].c_str());
#else
    GMX_UNUSED_VALUE(comm);
    ms = nullptr;
#endif

    return ms;
}

void done_multisim(gmx_multisim_t *ms)
{
    if (nullptr != ms)
    {
        done_mpi_in_place_buf(ms->mpb);
        delete ms;
    }
}
