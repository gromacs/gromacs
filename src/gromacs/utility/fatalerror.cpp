/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/errorcodes.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/mutex.h"
#include "gromacs/utility/programcontext.h"

#if GMX_MPI
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/gmxmpi.h"
#endif

#include "errorformat.h"

static bool       bDebug         = false;
static gmx::Mutex where_mutex;

FILE             *debug          = nullptr;
gmx_bool          gmx_debug_at   = FALSE;

static FILE      *log_file       = nullptr;
static gmx::Mutex error_mutex;

using Lock = gmx::lock_guard<gmx::Mutex>;

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
        Lock lock(where_mutex);
        if (bFirst) /* we repeat the check in the locked section because things
                       might have changed */
        {
            if ((temp = getenv("GMX_PRINT_DEBUG_LINES")) != nullptr)
            {
                nskip = strtol(temp, nullptr, 10);
            }
            bFirst = FALSE;
        }
    }

    // TODO None of this is thread safe, and presumably it was only
    // meant to run when debugging. But it runs many times every MD
    // step. Nice. See Redmine #2122.
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

static void default_error_handler(const char *title, const char *msg,
                                  const char *file, int line)
{
    if (log_file)
    {
        gmx::internal::printFatalErrorHeader(log_file, title, nullptr, file, line);
        gmx::internal::printFatalErrorMessageLine(log_file, msg, 0);
        gmx::internal::printFatalErrorFooter(log_file);
    }
    gmx::internal::printFatalErrorHeader(stderr, title, nullptr, file, line);
    gmx::internal::printFatalErrorMessageLine(stderr, msg, 0);
    gmx::internal::printFatalErrorFooter(stderr);
}

static gmx_error_handler_t gmx_error_handler = default_error_handler;

void gmx_set_error_handler(gmx_error_handler_t func)
{
    Lock lock(error_mutex);
    gmx_error_handler = func;
}

static const char *gmx_strerror(const char *key)
{
    struct ErrorKeyEntry {
        const char *key;
        const char *msg;
    };
    ErrorKeyEntry map[] = {
        { "call",   "Routine should not have been called" },
        { "comm",   "Communication (parallel processing) problem" },
        { "fatal",  "Fatal error" },
        { "file",   "File input/output error" },
        { "impl",   "Implementation restriction" },
        { "incons", "Software inconsistency error" },
        { "input",  "Input error or input inconsistency" },
        { "mem",    "Memory allocation/freeing error" },
        { "open",   "Cannot open file" },
        { "range",  "Range checking error" }
    };

    if (key == nullptr)
    {
        return "NULL error type (should not occur)";
    }
    for (const ErrorKeyEntry &entry : map)
    {
        if (std::strcmp(key, entry.key) == 0)
        {
            return entry.msg;
        }
    }
    return gmx::getErrorCodeString(gmx::eeUnknownError);
}

static void call_error_handler(const char *key, const char *file, int line, const char *msg)
{
    if (msg == nullptr)
    {
        msg = "Empty gmx_fatal message (bug).";
    }
    Lock lock(error_mutex);
    gmx_error_handler(gmx_strerror(key), msg, file, line);
}

void gmx_exit_on_fatal_error(ExitType exitType, int returnValue)
{
    if (log_file)
    {
        std::fflush(log_file);
    }
    if (debug)
    {
        std::fflush(debug);
    }
    std::fflush(stdout);
    std::fflush(stderr);

#if GMX_MPI
    if (gmx_mpi_initialized())
    {
        switch (exitType)
        {
            case ExitType_CleanExit:
                MPI_Finalize();
                break;
            case ExitType_Abort:
#if GMX_LIB_MPI
                gmx_abort(returnValue);
#endif
                break;
            case ExitType_NonMasterAbort:
                // Let all other processes wait till the master has printed
                // the error message and issued MPI_Abort.
                MPI_Barrier(MPI_COMM_WORLD);
                break;
        }
    }
#endif

    if (exitType == ExitType_CleanExit)
    {
        std::exit(returnValue);
    }
    // We cannot use std::exit() if other threads may still be executing, since that would cause destructors to be
    // called for global objects that may still be in use elsewhere.
    std::_Exit(returnValue);
}

void gmx_fatal_mpi_va(int /*f_errno*/, const char *file, int line,
                      gmx_bool bMaster, gmx_bool bFinalize,
                      const char *fmt, va_list ap)
{
    if (bMaster)
    {
        char msg[STRLEN];
        vsprintf(msg, fmt, ap);
        call_error_handler("fatal", file, line, msg);
    }

    ExitType exitType = ExitType_CleanExit;
    if (!bFinalize)
    {
        exitType = bMaster ? ExitType_Abort : ExitType_NonMasterAbort;
    }
    gmx_exit_on_fatal_error(exitType, 1);
}

void gmx_fatal(int f_errno, const char *file, int line, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    gmx_fatal_mpi_va(f_errno, file, line, TRUE, FALSE, fmt, ap);
    va_end(ap);
}

void _gmx_error(const char *key, const char *msg, const char *file, int line)
{
    call_error_handler(key, file, line, msg);
    gmx_exit_on_fatal_error(ExitType_Abort, 1);
}

void _range_check(int n, int n_min, int n_max, const char *warn_str,
                  const char *var, const char *file, int line)
{
    char buf[1024];

    if ((n < n_min) || (n >= n_max))
    {
        if (warn_str != nullptr)
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
