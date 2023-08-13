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
/*! \file
 * \brief
 * Declares fatal error handling and debugging routines for C code.
 *
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_FATALERROR_H
#define GMX_UTILITY_FATALERROR_H

#include <cstdarg>
#include <cstdio>

#include <filesystem>

#include "gromacs/libgromacs_export.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/stringutil.h"

/*! \brief
 * Debug log file.
 *
 * Functions can write to this file for debug info.
 * Before writing to it, it should be checked whether the file is not NULL:
 * \code
   if (debug)
   {
        std::fprintf(debug, "%s", "Debug text");
   }
   \endcode
 */
LIBGROMACS_EXPORT extern FILE* debug; //NOLINT(cppcoreguidelines-avoid-non-const-global-variables,-warnings-as-errors)
/** Whether extra debugging is enabled. */
extern gmx_bool gmx_debug_at; //NOLINT(cppcoreguidelines-avoid-non-const-global-variables,-warnings-as-errors)

/*! \brief
 * Initializes debugging variables.
 *
 * This function is not threadsafe.  It should be called as part of
 * initializing \Gromacs, before any other thread accesses the library.
 * For command line programs, gmx::CommandLineModuleManager takes care
 * of this if the user requests debugging.
 */
void gmx_init_debug(int dbglevel, const std::filesystem::path& dbgfile);

/** Returns TRUE when the program was started in debug mode */
gmx_bool bDebugMode();

/** Sets the log file for printing error messages. */
void gmx_fatal_set_log_file(FILE* fp);

/** Function pointer type for fatal error handler callback. */
typedef void (*gmx_error_handler_t)(const char*                  title,
                                    const std::string&           msg,
                                    const std::filesystem::path& file,
                                    int                          line);

/*! \brief
 * Sets an error handler for gmx_fatal() and other fatal error routines.
 *
 * The default handler prints the message.
 * \Gromacs will terminate the program after the error handler returns.
 * To make gmx_fatal_collective() work, the error handler should not terminate
 * the program, as it cannot know what is the desired way of termination.
 * The message passed to the handler may be a multi-line string.
 *
 * \see gmx_fatal()
 */
void gmx_set_error_handler(gmx_error_handler_t func);

/** Identifies the state of the program on a fatal error. */
enum ExitType
{
    /*! \brief
     * Clean exit is possible.
     *
     * There should be no concurrently executing code that might be accessing
     * global objects, and all MPI ranks should reach the same fatal error.
     */
    ExitType_CleanExit,
    /*! \brief
     * Program needs to be aborted.
     *
     * There are no preconditions for this state.
     */
    ExitType_Abort,
    /*! \brief
     * Program needs to be aborted, but some other rank is responsible of it.
     *
     * There should be some other MPI rank that reaches the same fatal error,
     * but uses ExitType_Abort.  The other ranks can then use
     * ExitType_NonMainAbort to wait for that one rank to issue the abort.
     */
    ExitType_NonMainAbort
};

/*! \brief
 * Helper function to terminate the program on a fatal error.
 *
 * \param[in] exitType  Identifies the state of the program at the time of the
 *    call, determining how the program can be terminated.
 * \param[in] returnValue  Exit code for the program, for cases where it can be
 *    used.
 */
[[noreturn]] void gmx_exit_on_fatal_error(enum ExitType exitType, int returnValue);

/*! \brief
 * Low-level fatal error reporting routine for collective MPI errors.
 *
 * This function works as gmx_fatal(), but provides additional control for
 * cases where it is known that the same error occurs on multiple MPI ranks.
 * The error handler is called only if \p bMain is `TRUE`, and MPI_Finalize()
 * is called instead of MPI_Abort() in MPI-enabled \Gromacs if \p bFinalize is
 * `TRUE`.
 *
 * This is used to implement gmx_fatal_collective() (which cannot be declared
 * here, since it would bring with it mdrun-specific dependencies).
 *
 * This function is deprecated and no new calls should be made to it.
 */
[[noreturn]] void gmx_fatal_mpi_va(int                          fatal_errno,
                                   const std::filesystem::path& file,
                                   int                          line,
                                   gmx_bool                     bMain,
                                   gmx_bool                     bFinalize,
                                   const char*                  fmt,
                                   std::va_list                 ap);

/*! \brief
 * Fatal error reporting routine for \Gromacs.
 *
 * This function prints a fatal error message with a header that contains the
 * source file and line number of the call, followed by the string specified by
 * \p fmt and supplied parameters.
 * If \p fatal_errno is 0, only the message and arguments are printed.
 * If \p fatal_errno is a legal system errno or -1, a perror()-like message is
 * printed after the first message; if fatal_errno is -1, the last system errno
 * will be used.
 * The format of \p fmt uses printf()-like formatting.
 *
 * In case all MPI processes want to stop with the same fatal error,
 * use gmx_fatal_collective(), declared in network.h,
 * to avoid having as many error messages as processes.
 *
 * This function is deprecated and no new calls should be made to it.
 *
 * The first three parameters can be provided through ::FARGS:
 * \code
   gmx_fatal(FARGS, fmt, ...);
   \endcode
 */
[[noreturn]] void
gmx_fatal(int fatal_errno, const std::filesystem::path& file, int line, gmx_fmtstr const char* fmt, ...)
        gmx_format(printf, 4, 5);
/** Helper macro to pass first three parameters to gmx_fatal(). */
#define FARGS 0, __FILE__, __LINE__

/*! \brief Implementation for gmx_error().
 *
 * This function is deprecated and no new calls should be made to it. */
[[noreturn]] void gmx_error_function(const char*                  key,
                                     const std::string&           msg,
                                     const std::filesystem::path& file,
                                     int                          line);
/*! \brief
 * Alternative fatal error routine with canned messages.
 *
 * This works as gmx_fatal(), except that a generic error message is added
 * based on a string key, and printf-style formatting is not supported.
 * Should not typically be called directly, but through gmx_call() etc.
 *
 * This macro is deprecated and no new calls should be made to it.
 */
#define gmx_error(key, msg) gmx_error_function(key, msg, __FILE__, __LINE__)

/*! \name Fatal error routines for certain types of errors
 *
 * These wrap gmx_error() and provide the \p key parameter as one of the
 * recognized strings.
 *
 * These macros are deprecated and no new calls should be made to them.
 */
/*! \{ */
#define gmx_call(msg) gmx_error("call", msg)
#define gmx_comm(msg) gmx_error("comm", msg)
#define gmx_file(msg) gmx_error("file", msg)
#define gmx_impl(msg) gmx_error("impl", msg)
#define gmx_incons(msg) gmx_error("incons", msg)
#define gmx_input(msg) gmx_error("input", msg)
#define gmx_mem(msg) gmx_error("mem", msg)
#define gmx_open(fn) gmx_error("open", fn)
/*! \} */

/*! \brief
 * Implementation for range_check() and range_check_mesg().
 *
 * \p warn_str can be NULL.
 */
void range_check_function(int                          n,
                          int                          n_min,
                          int                          n_max,
                          const char*                  warn_str,
                          const char*                  var,
                          const std::filesystem::path& file,
                          int                          line);

/*! \brief
 * Checks that a variable is within a range.
 *
 * If \p n is not in range [n_min, n_max), a fatal error is raised.
 * \p n_min is inclusive, but \p n_max is not.
 */
#define range_check_mesg(n, n_min, n_max, str) \
    range_check_function(n, n_min, n_max, str, #n, __FILE__, __LINE__)

/*! \brief
 * Checks that a variable is within a range.
 *
 * This works as range_check_mesg(), but with a default error message.
 */
#define range_check(n, n_min, n_max) \
    range_check_function(n, n_min, n_max, NULL, #n, __FILE__, __LINE__)

/*! \brief
 * Prints a warning message to stderr.
 *
 * The format of \p fmt uses printf()-like formatting.
 * The message string should NOT start with "WARNING"
 * and should NOT end with a newline.
 */
void gmx_warning(gmx_fmtstr const char* fmt, ...) gmx_format(printf, 1, 2);

#endif
