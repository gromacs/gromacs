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
 * Declares gmx::test::TestFileManager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
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
#ifndef GMX_TESTUTILS_TESTFILEMANAGER_H
#define GMX_TESTUTILS_TESTFILEMANAGER_H

#include <string>

#include "gromacs/utility/common.h"

namespace gmx
{
/*! \libinternal \brief
 * Testing utilities namespace.
 *
 * This namespace contains utilities for writing unit tests, mostly from the
 * \ref module_testutils module.
 */
namespace test
{

/*! \libinternal \brief
 * Helper for tests that need output files.
 *
 * To be used as a member in a test fixture class, this class provides
 * getTemporaryFilePath() method that returns a path for creating file names
 * for temporary files.  The returned path contains the name of the running
 * test, making it unique across tests.  Additionally, this class takes care of
 * removing any temporary files (i.e., all paths returned by
 * getTemporaryFilePath()) at test teardown (i.e., when the
 * TestFileManager is destructed).
 *
 * Functions getInputFilePath() and getInputDataDirectory() provide means to
 * access data files that are located in the test source directory.
 * This is used to provide input files for the tests, and also to store test
 * reference data persistently (see TestReferenceData).
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
         * Creates a file name for use within a single unit test.
         *
         * \param[in] suffix  Suffix to add to the file name (should contain an
         *      extension if one is desired).
         * \returns   File name that includes the test name and
         *      \p suffix.
         *
         * This method should only be called from within a Google Test test.
         * Two calls with the same \p suffix return the same string within the
         * same test.
         * Intended to produce distinct names for files that may be stored in
         * the same directory for multiple tests.
         */
        static std::string getTestSpecificFileName(const char *suffix);

        /*! \brief
         * Returns the path to a test input file.
         *
         * \param[in] filename  Relative path/filename to a test input file.
         * \returns Path to \p filename under the test input data directory.
         */
        static std::string getInputFilePath(const char *filename);

        /*! \brief
         * Returns the path to the test input directory.
         *
         * \returns Path to input data directory for the test executable.
         */
        static const char *getInputDataDirectory();

        /*! \brief
         * Sets the test input directory.
         *
         * \param[in] path  Path from which test input data is looked up from.
         *
         * \p path must name an existing directory.
         *
         * This function is automatically called by unittest_main.cpp through
         * initTestUtils().
         */
        static void setInputDataDirectory(const char *path);

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace test
} // namespace gmx

#endif
