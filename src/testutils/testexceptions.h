/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2015,2018,2019, by the GROMACS development team, led by
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
 * Exception classes for errors in tests.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_TESTEXCEPTIONS_H
#define GMX_TESTUTILS_TESTEXCEPTIONS_H

#include <string>

#include "gromacs/utility/exceptions.h"

namespace gmx
{
namespace test
{

/*! \libinternal \brief
 * Exception class for reporting errors in tests.
 *
 * This exception should be used for error conditions that are internal to the
 * test, i.e., do not indicate errors in the tested code.
 *
 * \ingroup module_testutils
 */
class TestException : public GromacsException
{
public:
    /*! \brief
     * Creates a test exception object with the provided detailed reason.
     *
     * \param[in] reason Detailed reason for the exception.
     */
    explicit TestException(const std::string& reason) : GromacsException(reason) {}
    /*! \brief
     * Creates a test exception based on another GromacsException object.
     *
     * \param[in] base  Exception to wrap.
     *
     * \see GMX_THROW_WRAPPER_TESTEXCEPTION
     */
    explicit TestException(const GromacsException& base) : GromacsException(base) {}

    int errorCode() const override { return -1; }
};

/*! \brief
 * Macro for throwing a TestException that wraps another exception.
 *
 * \param[in] e    Exception object to wrap.
 *
 * This macro is intended for wrapping exceptions thrown by Gromacs methods
 * that are called from a test for the test's internal purposes.  It wraps the
 * exception in a TestException to make it possible to tell from the type of
 * the exception whether the exception was thrown by the code under test, or by
 * the test code itself.
 *
 * \p e should evaluate to an instance of an object derived from
 * GromacsException.
 *
 * Typical usage in test code:
 * \code
   try
   {
       // some code that may throw a GromacsException
   }
   catch (const GromacsException &ex)
   {
       GMX_THROW_WRAPPER_TESTEXCEPTION(ex);
   }
 * \endcode
 */
#define GMX_THROW_WRAPPER_TESTEXCEPTION(e) throw ::gmx::test::TestException(e)

} // namespace test
} // namespace gmx

#endif
