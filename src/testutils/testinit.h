/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Functions for initialing \Gromacs unit test executables.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_TESTINIT_H
#define GMX_TESTUTILS_TESTINIT_H

namespace gmx
{

namespace test
{

//! \cond internal
/*! \internal
 * \brief
 * Initializes the test utilities library.
 *
 * Does not throw.  Terminates the program with a non-zero error code if an
 * error occurs.
 *
 * This function is automatically called by unittest_main.cpp.
 *
 * \ingroup module_testutils
 */
void initTestUtils(const char *dataPath, const char *tempPath, bool usesMpi,
                   bool usesHardwareDetection, int *argc, char ***argv);

/*! \internal
 * \brief
 * Finalizes the test utilities library.
 *
 * Does not throw.  Terminates the program with a non-zero error code if an
 * error occurs.
 *
 * This function is automatically called by unittest_main.cpp.
 *
 * \ingroup module_testutils
 */
void finalizeTestUtils();
//! \endcond

/*! \brief Declare a function that all unit test implementations can use
 * to set up any environment that they need.
 *
 * When registering the unit test in CMake, the HARDWARE_DETECTION
 * flag requires that the code for that unit test implements this
 * function.  Otherwise, a default stub implementation is provided.
 *
 * This approach conforms to the recommendation by GoogleTest to
 * arrange for the code that sets up the global test environment to be
 * called from main, rather than potentially rely on brittle static
 * initialization order. */
void callAddGlobalTestEnvironment();

} // namespace test
} // namespace gmx

#endif
