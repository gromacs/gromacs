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
 * Declares ::gmx::MessageStringCollector.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_MESSAGESTRINGCOLLECTOR_H
#define GMX_UTILITY_MESSAGESTRINGCOLLECTOR_H

#include <string>

#include "../utility/common.h"

namespace gmx
{

/*! \libinternal \brief
 * Helper class for collecting message strings, optionally with context.
 *
 * After strings have been collected, they can be formatted into one long
 * string for, e.g., printing out or for including in an exception.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class MessageStringCollector
{
    public:
        MessageStringCollector();
        ~MessageStringCollector();

        /*! \brief
         * Starts a context for messages.
         *
         * \param[in] name  Short description of the context.
         *
         * \see finishContext()
         * \see MessageStringContext
         */
        void startContext(const char *name);
        /*! \brief
         * Convenience wrapper for startContext(const char *).
         */
        void startContext(const std::string &name)
        { startContext(name.c_str()); }
        /*! \brief
         * Adds a new message.
         */
        void append(const char *message)
        { append(std::string(message)); }
        /*! \brief
         * Adds a new message.
         */
        void append(const std::string &message);
        /*! \brief
         * Ends a context started with startContext().
         *
         * \see MessageStringContext
         */
        void finishContext();
        /*! \brief
         * Clears all collected messages.
         */
        void clear();

        /*! \brief
         * Returns true if any messages have been added.
         *
         * \returns true if append() has been called at least once.
         *
         * The return value is identical to \c toString().empty().
         * Calls to startContext()/finishContext() only do not cause this
         * function to return true.
         */
        bool isEmpty() const;
        /*! \brief
         * Returns all collected messages as one string.
         */
        std::string toString() const;

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

/*! \libinternal \brief
 * Convenience class for creating a message context.
 *
 * This class provides a RAII-style interface to the
 * MessageStringCollector::startContext() and
 * MessageStringCollector::finishContext() methods: finishContext() is called
 * upon destruction of the object.  This avoids the need to call
 * MessageStringCollector::finishContext() on every possible exit point.
 *
 * Example usage:
 * \code
   bool function(::gmx::MessageStringCollector *errors)
   {
       ::gmx::MessageStringContext errcontext(errors, "In function()");
       bool bOk = function2(errors);
       bOk = function3(errors) && bOk;
       // <more processing>
       return bOk;
   }
 * \endcode
 *
 * \see MessageStringCollector
 * \inlibraryapi
 * \ingroup module_utility
 */
class MessageStringContext
{
    public:
        /*! \brief
         * Adds a context for the given object.
         */
        MessageStringContext(MessageStringCollector *collector, const char *name)
            : collector_(*collector)
        {
            collector_.startContext(name);
        }
        /*! \brief
         * Adds a context for the given object.
         */
        MessageStringContext(MessageStringCollector *collector,
                             const std::string      &name)
            : collector_(*collector)
        {
            collector_.startContext(name);
        }
        /*! \brief
         * Calls MessageStringCollector::finishContext() on the wrapped object.
         */
        ~MessageStringContext()
        {
            collector_.finishContext();
        }

    private:
        //! The wrapped object.
        MessageStringCollector &collector_;

        GMX_DISALLOW_COPY_AND_ASSIGN(MessageStringContext);
};

} // namespace gmx

#endif
