/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2017, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Provides helper classes for testing interactive prompts.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_INTERACTIVETEST_H
#define GMX_TESTUTILS_INTERACTIVETEST_H

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class TextInputStream;
class TextOutputStream;

namespace test
{

class TestReferenceChecker;

/*! \libinternal \brief
 * Helper class for testing interactive sessions.
 *
 * The calling test can set the user input using setInputLines() (and possibly
 * setLastNewline()), pass the streams from inputStream() and outputStream() to
 * the code that executes the interactive session, and then call checkSession()
 * after the session is finished.
 * The input is provided from the array set with setInputLines(), and all
 * output is checked using the reference data framework.
 * The reference XML data can be viewed with the XSLT stylesheet to show
 * exactly how the session went.
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class InteractiveTestHelper
{
    public:
        /*! \brief
         * Initializes the helper.
         *
         * \param[in] checker  Parent reference checker to use.
         *
         * The helper creates a compound item under \p checker for the
         * interactive session it tests.
         */
        explicit InteractiveTestHelper(gmx::test::TestReferenceChecker checker);
        ~InteractiveTestHelper();

        //! Sets whether the last input line contains a newline (by default, it does).
        void setLastNewline(bool bInclude);
        /*! \brief
         * Sets the input lines for the interactive session.
         *
         * Calls to TextInputStream::readLine() will return strings from this
         * array in sequence.
         * Newlines are added at the end automatically (except for the last
         * line if `setLastNewLine(false)` has been called).
         * If there are more `readLine()` calls than there are input lines,
         * the remaining calls return end-of-input.
         */
        void setInputLines(const ArrayRef<const char *const> &inputLines);

        //! Returns the input stream for the session.
        TextInputStream  &inputStream();
        //! Returns the output stream for the session.
        TextOutputStream &outputStream();

        /*! \brief
         * Finalizes the checking for the session.
         *
         * This must be called after all input and output from a session has
         * occurred, as the helper will not otherwise know when output after
         * the last input has finished.  This method also checks that the
         * required number of input lines were read in the session.
         */
        void checkSession();

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};
} // namespace test
} // namespace gmx

#endif
