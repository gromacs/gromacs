/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Declares fatal error handling and debugging routines for C code.
 *
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_FATALERROR_H
#define GMX_UTILITY_FATALERROR_H

#include <stdarg.h>
#include <stdio.h>

#include "gromacs/utility/basedefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief
 * Debug log file.
 *
 * Functions can write to this file for debug info.
 * Before writing to it, it should be checked whether the file is not NULL:
 * \code
   if (debug)
   {
       fprintf(debug, "%s", "Debug text");
   }
   \endcode
 */
extern FILE    *debug;
/** Whether extra debugging is enabled. */
extern gmx_bool gmx_debug_at;

/*! \brief
 * Initializes debugging variables.
 *
 * This function is not threadsafe.  It should be called as part of
 * initializing \Gromacs, before any other thread accesses the library.
 * For command line programs, gmx::CommandLineModuleManager takes care
 * of this if the user requests debugging.
 */
void gmx_init_debug(const int dbglevel, const char *dbgfile);

/** Returns TRUE when the program was started in debug mode */
gmx_bool bDebugMode(void);

/** Implementation for where(). */
void
_where(const char *file, int line);
/** Prints filename and line to stdlog. */
#define where() _where(__FILE__, __LINE__)

/** Sets the log file for printing error messages. */
void
gmx_fatal_set_log_file(FILE *fp);

/*! \brief
 * Sets an error handler for gmx_fatal() and other fatal error routines.
 *
 * The default handler prints the message.
 * \Gromacs will terminate the program after the error handler returns.
 * To make gmx_fatal_collective() work, the error handler should not terminate
 * the program, as it cannot know what is the desired way of termination.
 * The string passed to the handler may be a multi-line string.
 *
 * \see gmx_fatal()
 */
void
set_gmx_error_handler(void (*func)(const char *msg));

/*! \brief
 * Low-level fatal error reporting routine for collective MPI errors.
 *
 * This function works as gmx_fatal(), but provides additional control for
 * cases where it is known that the same error occurs on multiple MPI ranks.
 * The error handler is called only if \p bMaster is `TRUE`, and MPI_Finalize()
 * is called instead of MPI_Abort() in MPI-enabled \Gromacs if \p bFinalize is
 * `TRUE`.
 *
 * This is used to implement gmx_fatal_collective() (which cannot be declared
 * here, since it would bring with it mdrun-specific dependencies).
 */
gmx_noreturn void
gmx_fatal_mpi_va(int fatal_errno, const char *file, int line, gmx_bool bMaster,
                 gmx_bool bFinalize, const char *fmt, va_list ap);

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
 * The first three parameters can be provided through ::FARGS:
 * \code
   gmx_fatal(FARGS, fmt, ...);
   \endcode
 */
gmx_noreturn void
gmx_fatal(int fatal_errno, const char *file, int line, const char *fmt, ...);
/** Helper macro to pass first three parameters to gmx_fatal(). */
#define FARGS 0, __FILE__, __LINE__

/*! \brief
 * Returns error message corresponding to a string key.
 *
 * This maps the strings used by gmx_error() to actual error messages.
 * Caller is responsible of freeing the returned string.
 */
char *gmx_strerror(const char *key);

/** Implementation for gmx_error(). */
gmx_noreturn void _gmx_error(const char *key, const char *msg, const char *file, int line);
/*! \brief
 * Alternative fatal error routine with canned messages.
 *
 * This works as gmx_fatal(), except that a generic error message is added
 * based on a string key, and printf-style formatting is not supported.
 * Should not typically be called directly, but through gmx_bug(), gmx_call()
 * etc.
 */
#define gmx_error(key, msg) _gmx_error(key, msg, __FILE__, __LINE__)

/*! \name Fatal error routines for certain types of errors
 *
 * These wrap gmx_error() and provide the \p key parameter as one of the
 * recognized strings.
 */
/*! \{ */
#define gmx_bug(msg)    gmx_error("bug", msg)
#define gmx_call(msg)   gmx_error("call", msg)
#define gmx_comm(msg)   gmx_error("comm", msg)
#define gmx_file(msg)   gmx_error("file", msg)
#define gmx_cmd(msg)    gmx_error("cmd", msg)
#define gmx_impl(msg)   gmx_error("impl", msg)
#define gmx_incons(msg) gmx_error("incons", msg)
#define gmx_input(msg)  gmx_error("input", msg)
#define gmx_mem(msg)    gmx_error("mem", msg)
#define gmx_open(fn)    gmx_error("open", fn)
/*! \} */

/*! \brief
 * Implementation for range_check() and range_check_mesg().
 *
 * \p warn_str can be NULL.
 */
void _range_check(int n, int n_min, int n_max, const char *warn_str,
                  const char *var,
                  const char *file, int line);

/*! \brief
 * Checks that a variable is within a range.
 *
 * If \p n is not in range [n_min, n_max), a fatal error is raised.
 * \p n_min is inclusive, but \p n_max is not.
 */
#define range_check_mesg(n, n_min, n_max, str) _range_check(n, n_min, n_max, str,#n, __FILE__, __LINE__)

/*! \brief
 * Checks that a variable is within a range.
 *
 * This works as range_check_mesg(), but with a default error message.
 */
#define range_check(n, n_min, n_max) _range_check(n, n_min, n_max, NULL,#n, __FILE__, __LINE__)

/*! \brief
 * Prints a warning message to stderr.
 *
 * The format of \p fmt uses printf()-like formatting.
 * The message string should NOT start with "WARNING"
 * and should NOT end with a newline.
 */
void gmx_warning(const char *fmt, ...);

#ifdef __cplusplus
}
#endif

#endif
