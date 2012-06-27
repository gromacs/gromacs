/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief
 * Implements assertion handlers.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_utility
 */
#include "gromacs/utility/gmxassert.h"

#include <cstdio>
#include <cstdlib>

#include <string>

#include "gromacs/utility/stringutil.h"

#include "errorformat.h"

namespace gmx
{

/*! \cond internal */
namespace internal
{

void assertHandler(const char *condition, const char *msg,
                   const char *func, const char *file, int line)
{
    try
    {
        std::string title = formatString("Condition: %s\n%s", condition, msg);
        printFatalError(stderr, "Assertion failed", title.c_str(),
                        func, file, line);
    }
    catch (const std::bad_alloc &)
    {
        printFatalError(stderr, "Assertion failed",
                "(memory allocation failed while formatting the error message)",
                func, file, line);
    }
    std::abort();
}

} // namespace internal
//! \endcond

} // namespace gmx
