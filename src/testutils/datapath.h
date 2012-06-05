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
 * Functions for constructing file names for test files.
 *
 * Functions getTestFilePath() and getTestDataPath() provide means to access
 * data files that are located in the test source directory.  This is typically
 * used to provide input files for the tests.
 *
 * For temporary files used within a single test (typically in testing code
 * that writes into files), gmx::test::TestTemporaryFileManager is provided.
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

#include "gromacs/utility/common.h"

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
 * Helper for tests that need temporary output files.
 *
 * To be used as a member in a test fixture class, this class provides
 * getTemporaryFilePath() method that returns a path for creating file names
 * for temporary files.  The returned path contains the name of the running
 * test, making it unique across tests.  Additionally, this class takes care of
 * removing any temporary files (i.e., all paths returned by
 * getTemporaryFilePath()) at test teardown (i.e., when the
 * TestTemporaryFileManager is desctructed).
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class TestFileManager
{
    public:
        TestFileManager();
        /*! \brief
         * Frees internal storage and deletes any accessed file paths.
         *
         * Any errors (e.g., missing files) encountered while deleting the
         * files are ignored.
         */
        ~TestFileManager();

        /*! \brief
         * Creates a name for a temporary file within a single unit test.
         *
         * \param[in] suffix  Suffix to add to the file name (should contain an
         *      extension if one is desired).
         * \returns   Temporary file name that includes the test name and
         *      \p suffix.
         *
         * This method should only be called from within a Google Test test.
         * Two calls with the same \p suffix return the same string within the
         * same test.
         */
        std::string getTemporaryFilePath(const char *suffix);

        /*! \brief
         * Creates a name for reference data or temporary file within a single unit test.
         *
         * \param[in] suffix  Suffix to add to the file name (should contain an
         *      extension if one is desired).
         * \returns   File name that includes the test name and
         *      \p suffix.
         *
         * This method should only be called from within a Google Test test.
         * Two calls with the same \p suffix return the same string within the
         * same test.
         */
        static std::string getTestSpecificFileName(const char *suffix);

        /*! \libinternal \brief
         * Returns the path to a test input file.
         *
         * \param[in] filename  Relative path/filename to a test input file.
         * \returns Path to \p filename under the test input data directory.
         *
         * \inlibraryapi
         */
        static std::string getTestFilePath(const char *filename);

        /*! \libinternal \brief
         * Returns the path to the test input directory.
         *
         * \returns Path to input data directory for the test executable.
         *
         * \inlibraryapi
         */
        static const char *getTestDataPath();

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
        static void setTestDataPath(const char *path);

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace test
} // namespace gmx

#endif
