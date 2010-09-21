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
/*! \file
 * \brief
 * Defines ::gmx::ErrorContext.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_errorreporting
 */
#ifndef GMX_ERRORREPORTING_ERRORCONTEXT_H
#define GMX_ERRORREPORTING_ERRORCONTEXT_H

#include <string>

#include "abstracterrorreporter.h"

namespace gmx
{

/*! \brief
 * Convenience class for creating an error context.
 *
 * This class provides a RAII-style interface to the
 * AbstractErrorReporter::startContext() and
 * AbstractErrorReporter::finishContext() methods: finishContext() is called
 * upon destruction of the object.  This avoids the need to call
 * AbstractErrorReporter::finishContext() on every possible exit point.
 *
 * Example usage:
 * \code
int function(::gmx::AbstractErrorReporter *errors)
{
    ::gmx::ErrorContext errcontext(errors, "In function()");
    int rc;
    rc = function2(errors);
    if (rc != 0)
    {
        return rc;
    }
    rc = function3(errors);
    if (rc != 0)
    {
        return rc;
    }
    // <more processing>
    return 0;
}
 * \endcode
 *
 * \see AbstractErrorReporter
 * \inpublicapi
 * \ingroup module_errorreporting
 */
class ErrorContext
{
    public:
        /*! \brief
         * Adds a context for the given reporter.
         */
        ErrorContext(AbstractErrorReporter *reporter, const char *name)
            : _reporter(*reporter)
        {
            _reporter.startContext(name);
        }
        /*! \brief
         * Adds a context for the given reporter.
         */
        ErrorContext(AbstractErrorReporter *reporter, const std::string &name)
            : _reporter(*reporter)
        {
            _reporter.startContext(name);
        }
        /*! \brief
         * Calls AbstractReporter::finishContext() on the wrapped reporter.
         */
        ~ErrorContext()
        {
            _reporter.finishContext();
        }

    private:
        //! The wrapped reporter object.
        AbstractErrorReporter  &_reporter;

        // Disallow copy and assign.
        ErrorContext(const ErrorContext &);
        void operator =(const ErrorContext &);
};

} // namespace gmx

#endif
