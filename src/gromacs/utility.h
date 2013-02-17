/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
/*! \defgroup module_utility Low-level Utility Module
 * \ingroup group_utilitymodules
 * \brief
 * Provides various low-level utilities.
 *
 * <H3>Handling fatal errors</H3>
 *
 * Exception classes used in the library are defined in the exceptions.h header
 * file.  This header also declares a ::GMX_THROW macro that should be used for
 * throwing exceptions.  It also declares helper functions formatErrorMessage()
 * and translateException() for creating standard error messages and
 * translating exceptions to error return codes.
 *
 * Use of error return codes should be avoided in new code except in C wrappers
 * and similar, but for compatibility, facilities for handling them are
 * provided by the errorcodes.h header file.  It provides a set of error codes
 * (the enum \ref gmx::ErrorCode) that should be used for return codes in functions.
 * It also provides macros ::GMX_ERROR and ::GMX_ERROR_NORET that should be
 * used for returning an error code.  setFatalErrorHandler() is provided to
 * alter the behavior of ::GMX_ERROR and ::GMX_ERROR_NORET.  The default
 * handler prints the reason of the error to standard error and aborts the
 * execution.
 *
 * Header file gmxassert.h is also provided for assertions.  It declares macros
 * ::GMX_ASSERT and ::GMX_RELEASE_ASSERT that should be used for assertions.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 */
/*! \file
 * \brief
 * Public API convenience header for low-level utilities.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_H
#define GMX_UTILITY_H

#include "utility/errorcodes.h"
#include "utility/exceptions.h"
#include "utility/gmxassert.h"

#endif
