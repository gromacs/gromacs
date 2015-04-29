/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#include "gromacs/legacyheaders/main.h"

#include "config.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/fileio/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"
#include "gromacs/utility/sysinfo.h"

/* The source code in this file should be thread-safe.
         Please keep it that way. */

#define BUFSIZE 1024

static void par_fn(char *base, int ftp, const t_commrec *cr,
                   gmx_bool bAppendSimId, gmx_bool bAppendNodeId,
                   char buf[], int bufsize)
{
    if ((size_t)bufsize < (strlen(base)+10))
    {
        gmx_mem("Character buffer too small!");
    }

    /* Copy to buf, and strip extension */
    strcpy(buf, base);
    buf[strlen(base) - strlen(ftp2ext(fn2ftp(base))) - 1] = '\0';

    if (bAppendSimId)
    {
        sprintf(buf+strlen(buf), "%d", cr->ms->sim);
    }
    if (bAppendNodeId)
    {
        strcat(buf, "_rank");
        sprintf(buf+strlen(buf), "%d", cr->nodeid);
    }
    strcat(buf, ".");

    /* Add extension again */
    strcat(buf, (ftp == efTPR) ? "tpr" : (ftp == efEDR) ? "edr" : ftp2ext(ftp));
    if (debug)
    {
        fprintf(debug, "rank %d par_fn '%s'\n", cr->nodeid, buf);
        if (fn2ftp(buf) == efLOG)
        {
            fprintf(debug, "log\n");
        }
    }
}

void check_multi_int(FILE *log, const gmx_multisim_t *ms, int val,
                     const char *name,
                     gmx_bool bQuiet)
{
    int     *ibuf, p;
    gmx_bool bCompatible;

    if (NULL != log && !bQuiet)
    {
        fprintf(log, "Multi-checking %s ... ", name);
    }

    if (ms == NULL)
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
        if (NULL != log && !bQuiet)
        {
            fprintf(log, "OK\n");
        }
    }
    else
    {
        if (NULL != log)
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

void check_multi_int64(FILE *log, const gmx_multisim_t *ms,
                       gmx_int64_t val, const char *name,
                       gmx_bool bQuiet)
{
    gmx_int64_t      *ibuf;
    int               p;
    gmx_bool          bCompatible;

    if (NULL != log && !bQuiet)
    {
        fprintf(log, "Multi-checking %s ... ", name);
    }

    if (ms == NULL)
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
        if (NULL != log && !bQuiet)
        {
            fprintf(log, "OK\n");
        }
    }
    else
    {
        if (NULL != log)
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

    debug_gmx();

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
    debug_gmx();

    *fplog = fp;
}

void gmx_log_close(FILE *fp)
{
    if (fp)
    {
        gmx_fatal_set_log_file(NULL);
        gmx_fio_fclose(fp);
    }
}

void init_multisystem(t_commrec *cr, int nsim, char **multidirs,
                      int nfile, const t_filenm fnm[], gmx_bool bParFn)
{
    gmx_multisim_t *ms;
    int             nnodes, nnodpersim, sim, i, ftp;
    char            buf[256];
#ifdef GMX_MPI
    MPI_Group       mpi_group_world;
    int            *rank;
#endif

#ifndef GMX_MPI
    if (nsim > 1)
    {
        gmx_fatal(FARGS, "This binary is compiled without MPI support, can not do multiple simulations.");
    }
#endif

    nnodes  = cr->nnodes;
    if (nnodes % nsim != 0)
    {
        gmx_fatal(FARGS, "The number of ranks (%d) is not a multiple of the number of simulations (%d)", nnodes, nsim);
    }

    nnodpersim = nnodes/nsim;
    sim        = cr->nodeid/nnodpersim;

    if (debug)
    {
        fprintf(debug, "We have %d simulations, %d ranks per simulation, local simulation is %d\n", nsim, nnodpersim, sim);
    }

    snew(ms, 1);
    cr->ms   = ms;
    ms->nsim = nsim;
    ms->sim  = sim;
#ifdef GMX_MPI
    /* Create a communicator for the master nodes */
    snew(rank, ms->nsim);
    for (i = 0; i < ms->nsim; i++)
    {
        rank[i] = i*nnodpersim;
    }
    MPI_Comm_group(MPI_COMM_WORLD, &mpi_group_world);
    MPI_Group_incl(mpi_group_world, nsim, rank, &ms->mpi_group_masters);
    sfree(rank);
    MPI_Comm_create(MPI_COMM_WORLD, ms->mpi_group_masters,
                    &ms->mpi_comm_masters);

#if !defined(MPI_IN_PLACE_EXISTS)
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

#endif

    /* Reduce the intra-simulation communication */
    cr->sim_nodeid = cr->nodeid % nnodpersim;
    cr->nnodes     = nnodpersim;
#ifdef GMX_MPI
    MPI_Comm_split(MPI_COMM_WORLD, sim, cr->sim_nodeid, &cr->mpi_comm_mysim);
    cr->mpi_comm_mygroup = cr->mpi_comm_mysim;
    cr->nodeid           = cr->sim_nodeid;
#endif

    if (debug)
    {
        fprintf(debug, "This is simulation %d", cr->ms->sim);
        if (PAR(cr))
        {
            fprintf(debug, ", local number of ranks %d, local rank ID %d",
                    cr->nnodes, cr->sim_nodeid);
        }
        fprintf(debug, "\n\n");
    }

    if (multidirs)
    {
        if (debug)
        {
            fprintf(debug, "Changing to directory %s\n", multidirs[cr->ms->sim]);
        }
        gmx_chdir(multidirs[cr->ms->sim]);
    }
    else if (bParFn)
    {
        /* Patch output and tpx, cpt and rerun input file names */
        for (i = 0; (i < nfile); i++)
        {
            /* Because of possible multiple extensions per type we must look
             * at the actual file name
             */
            if (is_output(&fnm[i]) ||
                fnm[i].ftp == efTPR || fnm[i].ftp == efCPT ||
                strcmp(fnm[i].opt, "-rerun") == 0)
            {
                ftp = fn2ftp(fnm[i].fns[0]);
                par_fn(fnm[i].fns[0], ftp, cr, TRUE, FALSE, buf, 255);
                sfree(fnm[i].fns[0]);
                fnm[i].fns[0] = gmx_strdup(buf);
            }
        }
    }
}
