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
#include "gmxpre.h"

#include "fatalerror.h"

#include "config.h"

#include <cerrno>
#include <cstddef>
#include <cstdlib>
#include <cstring>

#include <exception>

#include "thread_mpi/threads.h"

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/programcontext.h"

#ifdef GMX_MPI
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/gmxmpi.h"
#endif

static bool                bDebug         = false;
static tMPI_Thread_mutex_t where_mutex    = TMPI_THREAD_MUTEX_INITIALIZER;

FILE                      *debug          = NULL;
gmx_bool                   gmx_debug_at   = FALSE;

static FILE               *log_file       = NULL;
static tMPI_Thread_mutex_t error_mutex    = TMPI_THREAD_MUTEX_INITIALIZER;
static const char *const   gmxuser
    = "Please report this to the mailing list (gmx-users@gromacs.org)";

void gmx_init_debug(const int dbglevel, const char *dbgfile)
{
    if (!bDebug)
    {
        gmx_disable_file_buffering();
        debug  = gmx_ffopen(dbgfile, "w+");
        bDebug = true;
        if (dbglevel >= 2)
        {
            gmx_debug_at = TRUE;
        }
    }
}

gmx_bool bDebugMode(void)
{
    return bDebug;
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

void gmx_fatal_set_log_file(FILE *fp)
{
    log_file = fp;
}

static int fatal_errno = 0;

static void default_error_handler(const char *msg)
{
    tMPI_Thread_mutex_lock(&error_mutex);
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
    tMPI_Thread_mutex_unlock(&error_mutex);
}

static void (*gmx_error_handler)(const char *msg) = default_error_handler;

void set_gmx_error_handler(void (*func)(const char *msg))
{
    // TODO: Either this is unnecessary, or also reads to the handler should be
    // protected by a mutex.
    tMPI_Thread_mutex_lock(&error_mutex);
    gmx_error_handler = func;
    tMPI_Thread_mutex_unlock(&error_mutex);
}

static void call_error_handler(const char *key, const char *file, int line, const char *msg)
{
    char        buf[10240], errerrbuf[1024];
    const char *llines = "-------------------------------------------------------";
    char       *strerr;

    if (msg == NULL)
    {
        sprintf(errerrbuf, "Empty fatal_error message. %s", gmxuser);
    }
    // In case ProgramInfo is not initialized and there is an issue with the
    // initialization, fall back to "GROMACS".
    const char *programName = "GROMACS";
    try
    {
        programName = gmx::getProgramContext().displayName();
    }
    catch (const std::exception &)
    {
    }

    strerr = gmx_strerror(key);
    sprintf(buf, "\n%s\nProgram %s, %s\n"
            "Source code file: %s, line: %d\n\n"
            "%s:\n%s\nFor more information and tips for troubleshooting, please check the GROMACS\n"
            "website at http://www.gromacs.org/Documentation/Errors\n%s\n",
            llines, programName, gmx_version(), file, line,
            strerr, msg ? msg : errerrbuf, llines);
    free(strerr);

    gmx_error_handler(buf);
}

gmx_noreturn static void do_exit(bool bMaster, bool bFinalize)
{
    if (debug)
    {
        fflush(debug);
    }

#ifdef GMX_MPI
    if (gmx_mpi_initialized())
    {
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
        else if (bMaster)
        {
#ifdef GMX_LIB_MPI
            gmx_abort(1);
#endif
        }
        else
        {
            /* Let all other processes wait till the master has printed
             * the error message and issued MPI_Abort.
             */
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
#else
    GMX_UNUSED_VALUE(bMaster);
    GMX_UNUSED_VALUE(bFinalize);
#endif

    if (bDebugMode())
    {
        std::abort();
    }
    std::exit(1);
}

void gmx_fatal_mpi_va(int f_errno, const char *file, int line,
                      gmx_bool bMaster, gmx_bool bFinalize,
                      const char *fmt, va_list ap)
{
    if (bMaster)
    {
        char msg[STRLEN];
        vsprintf(msg, fmt, ap);

        tMPI_Thread_mutex_lock(&error_mutex);
        fatal_errno = f_errno;
        tMPI_Thread_mutex_unlock(&error_mutex);

        call_error_handler("fatal", file, line, msg);
    }

    do_exit(bMaster, bFinalize);
}

void gmx_fatal(int f_errno, const char *file, int line, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    gmx_fatal_mpi_va(f_errno, file, line, TRUE, FALSE, fmt, ap);
    va_end(ap);
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
        { "range",  "Range checking error" },
        { NULL,     NULL}
    };

    if (key == NULL)
    {
        return strdup("Empty message");
    }
    else
    {
        for (size_t i = 0; msg[i].key != NULL; ++i)
        {
            if (strcmp(key, msg[i].key) == 0)
            {
                return strdup(msg[i].msg);
            }
        }
        char buf[1024];
        sprintf(buf, "No error message associated with key %s\n%s", key, gmxuser);
        return strdup(buf);
    }
}


void _gmx_error(const char *key, const char *msg, const char *file, int line)
{
    call_error_handler(key, file, line, msg);
    do_exit(true, false);
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
