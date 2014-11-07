/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014, by the GROMACS development team, led by
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
/*! \defgroup module_utility Low-Level Utilities (utility)
 * \ingroup group_utilitymodules
 * \brief
 * Provides low-level utilities for error handling and other tasks.
 *
 * This module provides various low-level utilities for tasks such as error
 * handling and string formatting, as well as helper classes and common custom
 * containers to simplify implementation of other code.  Contents of the module
 * are discussed in more details under the different headings below.
 * Some of the code in installed headers in the module is intended for use
 * directly from code outside the Gromacs library, but a significant portion is
 * exposed only because other public headers depend on it.
 *
 * Since this module implements error handling, it should be at the lowest
 * level: it should not depend on other modules.  Any functionality needed by
 * the error handling code should also be kept in this module.
 *
 * <H3>Error Handling</H3>
 *
 * Exception classes used in the library are declared in the exceptions.h header
 * file.  Most Gromacs-specific exceptions derive from gmx::GromacsException.
 *
 * This header also declares a ::GMX_THROW macro that should be used for
 * throwing exceptions.  ::GMX_THROW_WITH_ERRNO is also provided for reporting
 * syscall errors, but its use should be mostly limited to within the library.
 * This header also declares helper functions printFatalErrorMessage(),
 * formatExceptionMessageToString(), formatExceptionMessageToFile(), and
 * translateException() for creating standard error messages and translating
 * exceptions to error return codes.  processExceptionAtExit() provides
 * clean-up code before exiting the program after an exception.
 * To help in cases where bottom-up conversion to C++ is appropriate, macro
 * ::GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR is also provided to catch all
 * exceptions at C++ to C boundary.
 *
 * Use of error return codes should be avoided in new code except in C wrappers
 * and similar, but to ease migration (and for these cases where they are
 * necessary), facilities for handling them are provided by the errorcodes.h
 * header file.
 * It provides a set of error codes (the enum \ref gmx::ErrorCode) that should
 * be used for return codes in functions.
 * It also provides macros ::GMX_ERROR and ::GMX_ERROR_NORET that should be
 * used for returning an error code.  setFatalErrorHandler() is provided to
 * alter the behavior of ::GMX_ERROR and ::GMX_ERROR_NORET.  The default
 * handler prints the reason of the error to \c stderr and aborts the
 * execution.
 *
 * Header file gmxassert.h is also provided for assertions.  It declares macros
 * ::GMX_ASSERT and ::GMX_RELEASE_ASSERT that should be used for assertions.
 *
 * \if internal
 * Internally, functions from errorformat.h are used for all the above cases to
 * format error messages to \c stderr.
 * \endif
 *
 *
 * <H3>Basic %File Handling</H3>
 *
 * The header file.h declares a gmx::File class for basic I/O support.
 *
 * The header path.h declares helpers for manipulating paths and for managing
 * directories.
 *
 * The fate of these headers depends on what is decided in Redmine issue #950.
 *
 *
 * <H3>Implementation Helpers</H3>
 *
 * The header basedefinitions.h contains common definitions and macros used
 * throughout \Gromacs.  It includes fixed-width integer types (`gmx_int64_t`
 * and friends), `gmx_bool` for C code, some macros for compiler-specific
 * attributes, and ::GMX_UNUSED_VALUE and ::GMX_IGNORE_RETURN_VALUE for
 * handling warnings about unused values.
 *
 * The header classhelpers.h implements a gmx::PrivateImplPointer template for easily
 * writing classes that use the private implementation idiom.  This header also
 * declares ::GMX_DISALLOW_COPY_AND_ASSIGN and ::GMX_DISALLOW_ASSIGN macros for
 * class declarations.
 *
 * The header flags.h implements a gmx::FlagsTemplate template for better type
 * safety when using bit flag fields.
 *
 *
 * <H3>Other Functionality</H3>
 *
 * The header init.h declares gmx::init() and gmx::finalize() for initializing
 * and deinitializing the \Gromacs library.
 *
 * The header arrayref.h implements a gmx::ConstArrayRef class for exposing a
 * C array or part of a std::vector (basically, any continuous stretch of
 * memory) throuh a std::vector-like interface.
 *
 * The header stringutil.h declares various string utility routines.
 *
 * \if libapi
 *
 * The header gmxmpi.h abstracts away the mechanism of including either MPI or
 * thread-MPI headers depending on compilation options.
 * Similarly, gmxomp.h removes the need to use conditional compilation for code
 * that needs to include omp.h for OpenMP functions.
 *
 * The header gmxregex.h declares gmx::Regex and regexMatch() for basic regular
 * expression matching using an interface similar to C++11.
 *
 * The header messagestringcollector.h declares a gmx::MessageStringCollector
 * class for composing messages with context information.
 *
 * The header sysinfo.h declares gmx_getpid() for getting the current process
 * id.
 *
 * The header programcontext.h declares a gmx::ProgramContextInterface that is
 * used to
 * initialize and access information about the running program, such as the
 * name and path of the executable.  This information is used, e.g., by the
 * error handling code in formatting standard error messages.
 *
 * The header qsort_threadsafe.h provides a guaranteed threadsafe
 * implementation for qsort().
 *
 * The header uniqueptr.h declares gmx::gmx_unique_ptr, which is intended for
 * declaring smart pointer types with unique ownership.
 *
 * \endif
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 */
/*! \file
 * \brief
 * Public API convenience header for low-level utilities.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_H
#define GMX_UTILITY_H

#include "gromacs/utility/errorcodes.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/init.h"
#include "gromacs/utility/programcontext.h"

#endif
