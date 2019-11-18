/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2019, by the GROMACS development team, led by
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
 * Declares gmx::test::LoggerTestHelper.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_LOGGERTEST_H
#define GMX_TESTUTILS_LOGGERTEST_H

#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/logger.h"

namespace gmx
{

namespace test
{

/*! \libinternal \brief
 * Helper class for tests to check output written to a logger.
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class LoggerTestHelper
{
public:
    LoggerTestHelper();
    ~LoggerTestHelper();

    //! Returns the logger to pass to code under test.
    const MDLogger& logger();

    /*! \brief
     * Expects a log entry at a given level matching a given regex.
     *
     * Currently, the order of the entries is not checked, and if this
     * method is called once for a log level, then it needs to be called
     * for all entries produced by the test.
     *
     * If not called for a log level, all entries for that level are
     * accepted.
     */
    void expectEntryMatchingRegex(gmx::MDLogger::LogLevel level, const char* re);

private:
    class Impl;

    PrivateImplPointer<Impl> impl_;
};

} // namespace test
} // namespace gmx

#endif
