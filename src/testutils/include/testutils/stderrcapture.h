/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * Declares gmx::test::StderrCapture.
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_STDERRCAPTURE_H
#define GMX_TESTUTILS_STDERRCAPTURE_H

#include <string>

#include "gromacs/utility/classhelpers.h"

namespace gmx
{
namespace test
{

/*! \libinternal \brief
 * Captures everything written to stderr while an object of this class
 * is alive, or until stop() is called.
 *
 * This is useful for testing code, such as grompp, that reports
 * notes and warnings only by writing them to stderr. Only one object
 * of this class may be alive at any time. If the test fails before
 * stop() is called, the destructor stops the capture and discards
 * what was captured.
 *
 * The implementation depends on GoogleTest internals, namely
 * internal::testing::CaptureStderr because this works now, may be
 * portable, and can be reconsidered if/when GoogleTest changes that
 * internal functionality. Maybe we can drop tests that were using
 * this functionality, or implement this class with whatever
 * implementation GoogleTest had been using.
 *
 * A better long-term solution in cases of this type is to have the
 * test configure e.g. grompp with a WarningHandler whose contents can
 * later be checked, rather than try to capture the terminal output
 * produced by the default WarningHandler. However that solution is a
 * bit complex and should be deferred until we have enough tests that
 * would use it, or a test that is valuable enough that we want the
 * assurance that it will keep working.
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class StderrCapture
{
public:
    //! Starts capturing stderr
    StderrCapture();
    //! Stops capturing stderr, if stop() was not called
    ~StderrCapture();

    /*! \brief Stops capturing stderr and returns what was captured
     *
     * May only be called once.
     */
    std::string stop();

private:
    //! Whether stderr is still being captured
    bool isCapturing_ = true;

    GMX_DISALLOW_COPY_AND_ASSIGN(StderrCapture);
};

} // namespace test
} // namespace gmx

#endif
