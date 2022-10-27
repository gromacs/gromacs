/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2013- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \libinternal \file
 * \brief
 * Declares gmx::test::StdioTestHelper.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_STDIOHELPER_H
#define GMX_TESTUTILS_STDIOHELPER_H

#include <memory>

#include "gromacs/utility/classhelpers.h"

namespace gmx
{
namespace test
{

class TestFileManager;

/*! \libinternal \brief
 * Helper class for tests where code reads directly from `stdin`.
 *
 * Any method in this class may throw std::bad_alloc if out of memory.
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class StdioTestHelper
{
public:
    //! Creates a helper using the given file manager.
    explicit StdioTestHelper(TestFileManager* fileManager) : fileManager_(*fileManager) {}
    //! Destructor
    ~StdioTestHelper();

    /*! \brief Accepts a string as input, writes it to a temporary
     * file and then reopens stdin to read the contents of that
     * string.
     *
     * \throws FileIOError  when the freopen() fails
     */
    void redirectStringToStdin(const char* theString);

private:
    TestFileManager& fileManager_;
    bool             redirected = false;

    GMX_DISALLOW_COPY_AND_ASSIGN(StdioTestHelper);
};

} // namespace test
} // namespace gmx

#endif
