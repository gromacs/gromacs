/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
#include "gmx_fatal.h"
#include "gmx_fatal_collective.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <string.h>

#include "thread_mpi/threads.h"

#include "main.h"
#include "types/commrec.h"
#include "network.h"
#include "copyrite.h"
#include "macros.h"

#include "gromacs/fileio/futil.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/smalloc.h"

static gmx_bool            bDebug         = FALSE;
static FILE               *log_file       = NULL;

static tMPI_Thread_mutex_t debug_mutex     = TMPI_THREAD_MUTEX_INITIALIZER;
static tMPI_Thread_mutex_t where_mutex     = TMPI_THREAD_MUTEX_INITIALIZER;

gmx_bool bDebugMode(void)
{
    gmx_bool ret;
#if 0
    tMPI_Thread_mutex_lock(&debug_mutex);
#endif
    ret = bDebug;
#if 0
    tMPI_Thread_mutex_unlock(&debug_mutex);
#endif
    return bDebug;
}

void gmx_fatal_set_log_file(FILE *fp)
{
    log_file = fp;
}

void _where(const char *file, int line)
{
    static gmx_bool bFirst = TRUE;
    static int      nskip  = -1;
    static int      nwhere =  0;
    FILE           *fp;
    char           *temp;

    if (bFirst)
    {
        tMPI_Thread_mutex_lock(&where_mutex);
        if (bFirst) /* we repeat the check in the locked section because things
                       might have changed */
        {
            if ((temp = getenv("GMX_PRINT_DEBUG_LINES")) != NULL)
            {
                nskip = strtol(temp, NULL, 10);
            }
            bFirst = FALSE;
        }
        tMPI_Thread_mutex_unlock(&where_mutex);
    }

    if (nskip >= 0)
    {
        /* Skip the first n occasions, this allows to see where it goes wrong */
        if (nwhere >= nskip)
        {
            if (log_file)
            {
                fp = log_file;
            }
            else
            {
                fp = stderr;
            }
            fprintf(fp, "WHERE %d, file %s - line %d\n", nwhere, file, line);
        }
        nwhere++;
    }
}

static int fatal_errno = 0;

static void quit_gmx(const char *msg)
{
    tMPI_Thread_mutex_lock(&debug_mutex);
    if (fatal_errno == 0)
    {
        if (log_file)
        {
            fprintf(log_file, "%s\n", msg);
        }
        fprintf(stderr, "%s\n", msg);
        /* we set it to no-zero because if this function is called, something
           has gone wrong */
        fatal_errno = 255;
    }
    else
    {
        if (fatal_errno != -1)
        {
            errno = fatal_errno;
        }
        perror(msg);
    }

#ifdef GMX_LIB_MPI
    if (gmx_mpi_initialized())
    {
        int  nnodes;
        int  noderank;

        nnodes   = gmx_node_num();
        noderank = gmx_node_rank();

        if (nnodes > 1)
        {
            fprintf(stderr, "Error on rank %d, will try to stop all ranks\n",
                    noderank);
        }
        gmx_abort(noderank, nnodes, -1);
    }
#endif

    if (debug)
    {
        fflush(debug);
    }
    if (bDebugMode())
    {
        fprintf(stderr, "dump core (y/n):");
        fflush(stderr);
        if (toupper(getc(stdin)) != 'N')
        {
            (void) abort();
        }
    }

    exit(fatal_errno);
    tMPI_Thread_mutex_unlock(&debug_mutex);
}

/* The function below should be identical to quit_gmx,
 * except that is does not actually quit and call gmx_abort.
 */
static void quit_gmx_noquit(const char *msg)
{
    tMPI_Thread_mutex_lock(&debug_mutex);
    if (!fatal_errno)
    {
        if (log_file)
        {
            fprintf(log_file, "%s\n", msg);
        }
        fprintf(stderr, "%s\n", msg);
        /* we set it to no-zero because if this function is called, something
           has gone wrong */
        fatal_errno = 255;
    }
    else
    {
        if (fatal_errno != -1)
        {
            errno = fatal_errno;
        }
        perror(msg);
    }

#ifndef GMX_LIB_MPI
    if (debug)
    {
        fflush(debug);
    }
    if (bDebugMode())
    {
        fprintf(stderr, "dump core (y/n):");
        fflush(stderr);
        if (toupper(getc(stdin)) != 'N')
        {
            (void) abort();
        }
    }
#endif

    tMPI_Thread_mutex_unlock(&debug_mutex);
}

void gmx_fatal(int f_errno, const char *file, int line, const char *fmt, ...)
{
    va_list ap;
    char    msg[STRLEN];

    va_start(ap, fmt);
    vsprintf(msg, fmt, ap);
    va_end(ap);

    tMPI_Thread_mutex_lock(&debug_mutex);
    fatal_errno = f_errno;
    tMPI_Thread_mutex_unlock(&debug_mutex);

    _gmx_error("fatal", msg, file, line);
}

void gmx_fatal_collective(int f_errno, const char *file, int line,
                          const t_commrec *cr, gmx_domdec_t *dd,
                          const char *fmt, ...)
{
    gmx_bool    bFinalize;
    va_list     ap;
    char        msg[STRLEN];
#ifdef GMX_MPI
    int         result;
#endif

    bFinalize = TRUE;

#ifdef GMX_MPI
    /* Check if we are calling on all processes in MPI_COMM_WORLD */
    if (cr != NULL)
    {
        MPI_Comm_compare(cr->mpi_comm_mysim, MPI_COMM_WORLD, &result);
    }
    else
    {
        MPI_Comm_compare(dd->mpi_comm_all, MPI_COMM_WORLD, &result);
    }
    /* Any result except MPI_UNEQUAL allows us to call MPI_Finalize */
    bFinalize = (result != MPI_UNEQUAL);
#endif

    if ((cr != NULL && MASTER(cr)  ) ||
        (dd != NULL && DDMASTER(dd)))
    {
        va_start(ap, fmt);
        vsprintf(msg, fmt, ap);
        va_end(ap);

        tMPI_Thread_mutex_lock(&debug_mutex);
        fatal_errno = f_errno;
        tMPI_Thread_mutex_unlock(&debug_mutex);

        if (bFinalize)
        {
            /* Use an error handler that does not quit */
            set_gmx_error_handler(quit_gmx_noquit);
        }

        _gmx_error("fatal", msg, file, line);
    }

#ifdef GMX_MPI
    if (bFinalize)
    {
        /* Broadcast the fatal error number possibly modified
         * on the master process, in case the user would like
         * to use the return status on a non-master process.
         * The master process in cr and dd always has global rank 0.
         */
        MPI_Bcast(&fatal_errno, sizeof(fatal_errno), MPI_BYTE,
                  0, MPI_COMM_WORLD);

        /* Finalize nicely instead of aborting */
        MPI_Finalize();
    }
    else
    {
        /* Let all other processes wait till the master has printed
         * the error message and issued MPI_Abort.
         */
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    exit(fatal_errno);
}

/*
 * These files are global variables in the gromacs preprocessor
 * Every routine in a file that includes gmx_fatal.h can write to these
 * debug channels. Depending on the debuglevel used
 * 0 to 3 of these filed are redirected to /dev/null
 *
 */
FILE    *debug           = NULL;
gmx_bool gmx_debug_at    = FALSE;

void init_debug(const int dbglevel, const char *dbgfile)
{
    tMPI_Thread_mutex_lock(&debug_mutex);
    if (!bDebug) /* another thread hasn't already run this*/
    {
        no_buffers();
        debug  = gmx_fio_fopen(dbgfile, "w+");
        bDebug = TRUE;
        if (dbglevel >= 2)
        {
            gmx_debug_at = TRUE;
        }
    }
    tMPI_Thread_mutex_unlock(&debug_mutex);
}

static const char *gmxuser = "Please report this to the mailing list (gmx-users@gromacs.org)";

static void        (*gmx_error_handler)(const char *msg) = quit_gmx;

void set_gmx_error_handler(void (*func)(const char *msg))
{
    tMPI_Thread_mutex_lock(&debug_mutex);
    gmx_error_handler = func;
    tMPI_Thread_mutex_unlock(&debug_mutex);
}

char *gmx_strerror(const char *key)
{
    typedef struct {
        const char *key, *msg;
    } error_msg_t;
    error_msg_t msg[] = {
        { "bug",    "Possible bug" },
        { "call",   "Routine should not have been called" },
        { "comm",   "Communication (parallel processing) problem" },
        { "fatal",  "Fatal error" },
        { "cmd",    "Invalid command line argument" },
        { "file",   "File input/output error" },
        { "impl",   "Implementation restriction" },
        { "incons", "Software inconsistency error" },
        { "input",  "Input error or input inconsistency" },
        { "mem",    "Memory allocation/freeing error" },
        { "open",   "Can not open file" },
        { "range",  "Range checking error" }
    };
#define NMSG asize(msg)
    char        buf[1024];
    size_t      i;

    if (key == NULL)
    {
        return strdup("Empty message");
    }
    else
    {
        for (i = 0; (i < NMSG); i++)
        {
            if (strcmp(key, msg[i].key) == 0)
            {
                break;
            }
        }
        if (i == NMSG)
        {
            sprintf(buf, "No error message associated with key %s\n%s", key, gmxuser);
            return strdup(buf);
        }
        else
        {
            return strdup(msg[i].msg);
        }
    }
}


void _gmx_error(const char *key, const char *msg, const char *file, int line)
{
    char        buf[10240], errerrbuf[1024];
    const char *llines = "-------------------------------------------------------";
    char       *strerr;

    /* protect the audience from suggestive discussions */

    if (msg == NULL)
    {
        sprintf(errerrbuf, "Empty fatal_error message. %s", gmxuser);
    }

    strerr = gmx_strerror(key);
    sprintf(buf, "\n%s\nProgram %s, %s\n"
            "Source code file: %s, line: %d\n\n"
            "%s:\n%s\nFor more information and tips for troubleshooting, please check the GROMACS\n"
            "website at http://www.gromacs.org/Documentation/Errors\n%s\n",
            llines, ShortProgram(), GromacsVersion(), file, line,
            strerr, msg ? msg : errerrbuf, llines);
    free(strerr);

    gmx_error_handler(buf);
}

void _range_check(int n, int n_min, int n_max, const char *warn_str,
                  const char *var, const char *file, int line)
{
    char buf[1024];

    if ((n < n_min) || (n >= n_max))
    {
        if (warn_str != NULL)
        {
            strcpy(buf, warn_str);
            strcat(buf, "\n");
        }
        else
        {
            buf[0] = '\0';
        }

        sprintf(buf+strlen(buf), "Variable %s has value %d. It should have been "
                "within [ %d .. %d ]\n", var, n, n_min, n_max);

        _gmx_error("range", buf, file, line);
    }
}

void gmx_warning(const char *fmt, ...)
{
    va_list ap;
    char    msg[STRLEN];

    va_start(ap, fmt);
    vsprintf(msg, fmt, ap);
    va_end(ap);

    fprintf(stderr, "\nWARNING: %s\n\n", msg);
}
