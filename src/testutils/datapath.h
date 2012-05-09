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
/*! \libinternal \file
 * \brief
 * Functions for accessing test input files.
 *
 * Functions in this header provide methods to access data files that are
 * located in the test source directory.  This is typically used to provide
 * input files for the tests.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_testutils
 */
/*! \libinternal \defgroup module_testutils Unit Testing Utilities
 * \brief
 * Common helper classes and functions for writing unit tests.
 *
 * To build unit tests using these utilities, libxml2 is required.
 *
 * \ingroup group_utilitymodules
 */
#ifndef GMX_TESTUTILS_DATAPATH_H
#define GMX_TESTUTILS_DATAPATH_H

#include <string>

namespace gmx
{
/*! \libinternal \brief
 * Namespace for unit testing utilities.
 *
 * This namespace contains utilities that are shared between unit tests.
 * Most members are declared in the \ref module_testutils module, but some
 * are also declared within individual tests (these are typically candidates
 * for using in other tests as well).
 *
 * \ingroup module_testutils
 */
namespace test
{

/*! \libinternal \brief
 * Returns the path to a test input file.
 *
 * \param[in] filename  Relative path/filename to a test input file.
 * \returns Path to \p filename under the test input data directory.
 *
 * \inlibraryapi
 */
std::string getTestFilePath(const char *filename);
/*! \libinternal \brief
 * Returns the path to the test input directory.
 *
 * \returns Path to input data directory for the test executable.
 *
 * \inlibraryapi
 */
const char *getTestDataPath();
/*! \libinternal \brief
 * Sets the test input directory.
 *
 * \param[in] path  Path from which test input data is looked up from.
 *
 * \p path must name an existing directory.
 *
 * This function is automatically called by test_main_gtest.cpp and
 * test_main_gmock.cpp.
 *
 * \inlibraryapi
 */
void setTestDataPath(const char *path);

} // namespace test
} // namespace gmx

#endif
