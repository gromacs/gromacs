/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * Dummy header for \ref module_testutils documentation.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
/*! \libinternal \defgroup module_testutils Testing Utilities (testutils)
 * \brief
 * Common helper classes and functions for writing tests using Google Test.
 *
 * General documentation for the testing in \Gromacs can be found in the
 * \linktodevmanual{index,developer guide}.  This page provides an overview of
 * the actual API provided by the `testutils` module.
 *
 * There are several distinct functionalities provided:
 *  - gmx::test::TestFileManager (in testfilemanager.h) provides functionality
 *    for locating test input files from the source directory and managing
 *    temporary files that need to be created during the test.
 *  - gmx::test::TestFileInputRedirector (in testfileredirector.h) provides
 *    functionality for capturing file existence checks in code that uses
 *    gmx::FileInputRedirectorInterface.
 *  - #GMX_TEST_OPTIONS macro provides facilities for adding custom command
 *    line options for the test binary.
 *  - testasserts.h provides several custom test assertions for better
 *    exception and floating-point handling than built-in Google Test
 *    assertions.
 *  - gmx::test::TestReferenceData and related classes (in refdata.h) provide
 *    utilities to write regression-style tests that check that the test
 *    produces the same results as an earlier run of the same test.  The
 *    reference data is stored as XML files.  For certain types of tests, the
 *    XML can also be easier to inspect manually for correctness than writing
 *    the checks in C++, providing an alternative method to write assertions.
 *  - gmx::test::CommandLine and related classes (in cmdlinetest.h) provide
 *    utilities for constructing command line argument arrays for use in tests
 *    that invoke actual commands.
 *  - gmx::test::StringTestBase provides a test fixture for tests that need to
 *    test long strings for correctness.
 *  - gmx::test::IntegrationTestFixture provides a test fixture for tests that
 *    execute legacy code where `stdin` reading etc. cannot be easily mocked.
 *
 * Additionally, testinit.h and mpi-printer.h, and their corresponding source
 * files, provide functionality that is not visible on the API level: they
 * provide initialization routines for the above functionality, which are
 * automatically called by the %main() function provided in unittest_main.cpp.
 *
 * mpi-printer.h provides a Google Test listener that is installed when the
 * tests are compiled with MPI.  This listener allows the test binary to be run
 * on multiple MPI ranks, and synchronizes the execution and output from the
 * test cases, as well as makes the test fail on even if an assertion fails
 * only on one rank.
 *
 * \ingroup group_utilitymodules
 */
