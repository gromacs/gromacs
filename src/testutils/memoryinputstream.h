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
/*! \libinternal \file
 * \brief
 * Provides helper class for testing reading of text files.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_MEMORYINPUTSTREAM_H
#define GMX_TESTUTILS_MEMORYINPUTSTREAM_H

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class TextInputStream;

namespace test
{

/*! \libinternal \brief
 * Helper class for testing file reading
 *
 * The calling test can set the user input using setInputLines(), and
 * pass the stream from inputStream() to the code that needs one.
 *
 * \todo This class duplicates part of `InteractiveTestHelper`, but it
 * is unclear if anything can or should be done about it.
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class MemoryInputStream
{
    public:
        /*! \brief Initializes the memory stream. */
        MemoryInputStream();
        ~MemoryInputStream();

        /*! \brief
         * Sets the input lines for the interactive session.
         *
         * Calls to TextInputStream::readLine() will return strings from this
         * array in sequence.
         * Newlines are added at the end automatically.
         * If there are more `readLine()` calls than there are input lines,
         * the remaining calls return end-of-input.
         */
        void setInputLines(const ConstArrayRef<const char *> &inputLines);

        //! Returns the input stream for the session.
        TextInputStream  &inputStream();

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};
} // namespace test
} // namespace gmx

#endif
