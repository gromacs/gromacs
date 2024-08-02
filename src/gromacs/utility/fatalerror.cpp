/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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

#include "gromacs/utility/fatalerror.h"

#include "config.h"

#include <cerrno>
#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <exception>
#include <filesystem>
#include <mutex>
#include <string>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/stringutil.h"

#include "errorcodes.h"

#if GMX_MPI
#    include "gromacs/utility/basenetwork.h"
#    include "gromacs/utility/gmxmpi.h"
#endif

#include "errorformat.h"

static bool bDebug = false; // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)

FILE* debug        = nullptr; // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)
bool  gmx_debug_at = false;   // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)

static FILE*      log_file = nullptr; // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)
static std::mutex error_mutex;        // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)

using Lock = std::lock_guard<std::mutex>;

void gmx_init_debug(const int dbglevel, const std::filesystem::path& dbgfile)
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

gmx_bool bDebugMode()
{
    return bDebug;
}

void gmx_fatal_set_log_file(FILE* fp)
{
    log_file = fp;
}

static void default_error_handler(const char*                  title,
                                  const std::string&           msg,
                                  const std::filesystem::path& file,
                                  int                          line)
{
    if (log_file)
    {
        gmx::internal::printFatalErrorHeader(log_file, title, nullptr, file.string().c_str(), line);
        gmx::internal::printFatalErrorMessageLine(log_file, msg.c_str(), 0);
        gmx::internal::printFatalErrorFooter(log_file);
    }
    gmx::internal::printFatalErrorHeader(stderr, title, nullptr, file.string().c_str(), line);
    gmx::internal::printFatalErrorMessageLine(stderr, msg.c_str(), 0);
    gmx::internal::printFatalErrorFooter(stderr);
}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static gmx_error_handler_t gmx_error_handler = default_error_handler;

void gmx_set_error_handler(gmx_error_handler_t func)
{
    Lock lock(error_mutex);
    gmx_error_handler = func;
}

static const char* gmx_strerror(const char* key)
{
    struct ErrorKeyEntry
    {
        const char* key;
        const char* msg;
    };
    ErrorKeyEntry map[] = { { "call", "Routine should not have been called" },
                            { "comm", "Communication (parallel processing) problem" },
                            { "fatal", "Fatal error" },
                            { "file", "File input/output error" },
                            { "impl", "Implementation restriction" },
                            { "incons", "Software inconsistency error" },
                            { "input", "Input error or input inconsistency" },
                            { "mem", "Memory allocation/freeing error" },
                            { "open", "Cannot open file" },
                            { "range", "Range checking error" } };

    if (key == nullptr)
    {
        return "NULL error type (should not occur)";
    }
    for (const ErrorKeyEntry& entry : map)
    {
        if (std::strcmp(key, entry.key) == 0)
        {
            return entry.msg;
        }
    }
    return gmx::getErrorCodeString(gmx::eeUnknownError);
}

static void call_error_handler(const char* key, const std::filesystem::path& file, int line, const std::string& msg)
{
    Lock lock(error_mutex);
    gmx_error_handler(gmx_strerror(key), msg.empty() ? "Empty gmx_fatal message (bug)." : msg, file, line);
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
            case ExitType_CleanExit: MPI_Finalize(); break;
            case ExitType_Abort:
#    if GMX_LIB_MPI
                gmx_abort(returnValue);
#    else
                break;
#    endif
            case ExitType_NonMainAbort:
                // Let all other processes wait till the main has printed
                // the error message and issued MPI_Abort.
                MPI_Barrier(MPI_COMM_WORLD);
                break;
        }
    }
#endif

    if (!GMX_FAHCORE)
    {
        if (exitType == ExitType_CleanExit)
        {
            std::exit(returnValue);
        }
        // We cannot use std::exit() if other threads may still be executing, since that would cause
        // destructors to be called for global objects that may still be in use elsewhere.
        std::_Exit(returnValue);
    }
}

void gmx_fatal_mpi_va(int /*f_errno*/,
                      const std::filesystem::path& file,
                      int                          line,
                      gmx_bool                     bMain,
                      gmx_bool                     bFinalize,
                      const char*                  fmt,
                      std::va_list                 ap)
{
    if (bMain)
    {
        std::string msg = gmx::formatStringV(fmt, ap);
        call_error_handler("fatal", file, line, msg);
    }

    ExitType exitType = ExitType_CleanExit;
    if (!bFinalize)
    {
        exitType = bMain ? ExitType_Abort : ExitType_NonMainAbort;
    }
    gmx_exit_on_fatal_error(exitType, 1);
}

void gmx_fatal(int f_errno, const std::filesystem::path& file, int line, gmx_fmtstr const char* fmt, ...)
{
    std::va_list ap;
    va_start(ap, fmt);
    gmx_fatal_mpi_va(f_errno, file, line, TRUE, FALSE, fmt, ap);
    va_end(ap);
}

void gmx_error_function(const char* key, const std::string& msg, const std::filesystem::path& file, int line)
{
    call_error_handler(key, file, line, msg);
    gmx_exit_on_fatal_error(ExitType_Abort, 1);
}

void range_check_function(int                          n,
                          int                          n_min,
                          int                          n_max,
                          const char*                  warn_str,
                          const char*                  var,
                          const std::filesystem::path& file,
                          int                          line)
{
    if ((n < n_min) || (n >= n_max))
    {
        std::string buf;
        if (warn_str != nullptr)
        {
            buf = warn_str;
            buf += "\n";
        }

        buf += gmx::formatString(
                "Variable %s has value %d. It should have been "
                "within [ %d .. %d ]\n",
                var,
                n,
                n_min,
                n_max);

        gmx_error_function("range", buf, file, line);
    }
}

void gmx_warning(gmx_fmtstr const char* fmt, ...)
{
    std::va_list ap;
    char         msg[STRLEN];

    va_start(ap, fmt);
    vsprintf(msg, fmt, ap);
    va_end(ap);

    fprintf(stderr, "\nWARNING: %s\n\n", msg);
}
