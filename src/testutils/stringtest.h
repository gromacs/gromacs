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
 * Declares gmx::test::StringTestBase.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_STRINGTEST_H
#define GMX_TESTUTILS_STRINGTEST_H

#include <string>

#include <boost/scoped_ptr.hpp>
#include <gtest/gtest.h>

#include "refdata.h"

namespace gmx
{
namespace test
{

/*! \libinternal \brief
 * Test fixture for tests that check string formatting.
 *
 * For development, tests that use this fixture as their base can be run with a
 * '-stdout' command-line option to print out the tested strings to stdout.
 * If this flag is not given, they check the strings using the XML reference
 * framework (see TestReferenceData).
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class StringTestBase : public ::testing::Test
{
    public:
        //! Static fixture setup to parse command-line options.
        static void SetUpTestCase();

        StringTestBase();
        ~StringTestBase();

        /*! \brief
         * Returns the root checker for this test's reference data.
         *
         * Can be used to perform custom checks against reference data (e.g.,
         * if the test needs to check some other values than plain strings.
         */
        TestReferenceChecker &checker();

        /*! \brief
         * Check a string.
         *
         * \param[in] text  String to check.
         * \param[in] id    Unique (within a single test) id for the string.
         */
        void checkText(const std::string &text, const char *id);
        /*! \brief
         * Check contents of a file as a single string.
         *
         * \param[in] filename  Name of the file to check.
         * \param[in] id        Unique (within a single test) id for the string.
         *
         * Provided for convenience.  Reads the contents of \p filename into a
         * single string and calls checkText().
         */
        void checkFileContents(const std::string &filename, const char *id);

    private:
        static bool                             s_bWriteToStdOut;

        TestReferenceData                       data_;
        boost::scoped_ptr<TestReferenceChecker> checker_;
};

} // namespace test
} // namespace gmx

#endif
