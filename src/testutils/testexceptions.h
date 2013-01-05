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
 * Exception classes for errors in tests.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
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
        explicit TestException(const std::string &reason)
            : GromacsException(reason) {}
        /*! \brief
         * Creates a test exception based on another GromacsException object.
         *
         * \param[in] base  Exception to wrap.
         *
         * \see GMX_THROW_WRAPPER_TESTEXCEPTION
         */
        explicit TestException(const GromacsException &base)
            : GromacsException(base) {}

        virtual int errorCode() const { return -1; }
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
#define GMX_THROW_WRAPPER_TESTEXCEPTION(e) \
    throw ::boost::enable_current_exception(::gmx::test::TestException(e))

} // namespace test
} // namespace gmx

#endif
