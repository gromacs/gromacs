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
 * Defines ::gmx::AbstractErrorReporter.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_errorreporting
 */
#ifndef GMX_ERRORREPORTING_ABSTRACTERRORREPORTER_H
#define GMX_ERRORREPORTING_ABSTRACTERRORREPORTER_H

#include <string>

namespace gmx
{

/*! \brief
 * Abstract base class for error reporters.
 *
 * This class provides an interface for reporting non-fatal errors from a
 * complex function.  Such a function should take a pointer to an
 * AbstractReporter object, and use the provided methods to report any errors
 * it encounters.  If the function calls other functions that can also detect
 * errors, it can pass the reporter object to these functions as well, possibly
 * after adding context information using startContext()/finishContext() or the
 * ErrorContext class.  The caller of such a function can then create an error
 * reporter of their choice and pass it to the function to alter how errors are
 * reported to the user.
 *
 * This abstract implementation provides basic facilities to implement counting
 * of errors of different type.  Derived classes should call the
 * incrementCount() method from their implementation of the add() method for
 * this to work properly.
 *
 * \see ErrorContext
 * \inpublicapi
 * \ingroup module_errorreporting
 */
class AbstractErrorReporter
{
    public:
        /*! \brief
         * Type of error.
         */
        enum ErrorType
        {
            etNote,
            etWarning,
            etError,
            etNR
        };

        virtual ~AbstractErrorReporter() {}

        /*! \brief
         * Starts a context for errors.
         *
         * \param[in] name  Short description of the context.
         *
         * Derived classes should duplicate the string if they need to store
         * it.
         *
         * \see ErrorContext
         */
        virtual void startContext(const char *name) = 0;
        /*! \brief
         * Reports a new error.
         *
         * \param[in] type    Type of the error.
         * \param[in] reason  Short description of the error.
         *
         * Derived classes should duplicate the string if they need to store
         * it.
         */
        virtual void add(ErrorType type, const char *reason) = 0;
        /*! \brief
         * Ends a context started with startContext().
         *
         * \see ErrorContext
         */
        virtual void finishContext() = 0;

        /*! \brief
         * Convenience wrapper for startContext(const char *).
         */
        void startContext(const std::string &name) { startContext(name.c_str()); }
        /*! \brief
         * Convenience wrapper for add(ErrorType, const char *).
         */
        void add(ErrorType type, const std::string &reason) { add(type, reason.c_str()); }
        /*! \brief
         * Convenience wrapper for add() for adding notes.
         */
        void note(const char *reason) { add(etNote, reason); }
        //! \copydoc note(const char *)
        void note(const std::string &reason) { add(etNote, reason); }
        /*! \brief
         * Convenience wrapper for add() for adding warnings.
         */
        void warning(const char *reason) { add(etWarning, reason); }
        //! \copydoc warning(const char *)
        void warning(const std::string &reason) { add(etWarning, reason); }
        /*! \brief
         * Convenience wrapper for add() for adding errors.
         */
        void error(const char *reason) { add(etError, reason); }
        //! \copydoc error(const char *)
        void error(const std::string &reason) { add(etError, reason); }

        /*! \brief
         * Returns the number of errors of particular type that have occurred.
         */
        int errorCount(ErrorType type) const { return _count[type]; }

    protected:
        /*! \brief
         * Initializes a reporter with zero error counts.
         */
        AbstractErrorReporter()
        {
            for (int i = 0; i < etNR; ++i)
            {
                _count[i] = 0;
            }
        }

        /*! \brief
         * Increment the count of a particular error type.
         *
         * Derived classes should call this function from their implementation
         * of the add() method.
         */
        void incrementCount(ErrorType type) { _count[type]++; }

    private:
        //! Number of errors occurred for each type.
        int                      _count[etNR];
};

} // namespace gmx

#endif
