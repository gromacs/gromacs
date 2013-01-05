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
 * Declares gmx::HelpManager.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_onlinehelp
 */
#ifndef GMX_ONLINEHELP_HELPMANAGER_H
#define GMX_ONLINEHELP_HELPMANAGER_H

#include <string>

#include "../utility/common.h"

namespace gmx
{

class HelpTopicInterface;
class HelpWriterContext;

/*! \libinternal \brief
 * Helper for providing interactive online help.
 *
 * \inlibraryapi
 * \ingroup module_onlinehelp
 */
class HelpManager
{
    public:
        /*! \brief
         * Creates a manager that uses a given root topic.
         *
         * \param[in] rootTopic  Help topic that can be accessed through this
         *      manager.
         * \param[in] context    Context object for writing the help.
         * \throws    std::bad_alloc if out of memory.
         *
         * The provided topic and context objects must remain valid for the
         * lifetime of this manager object.
         */
        HelpManager(const HelpTopicInterface &rootTopic,
                    const HelpWriterContext  &context);
        ~HelpManager();

        /*! \brief
         * Enters a subtopic with the given name under the active topic.
         *
         * \param[in] name  Subtopic name to enter.
         * \throws    std::bad_allod if out of memory.
         * \throws    InvalidInputError if topic with \p name is not found.
         */
        void enterTopic(const char *name);
        //! \copydoc enterTopic(const char *)
        void enterTopic(const std::string &name);

        /*! \brief
         * Writes out the help for the currently active topic.
         *
         * \throws  std::bad_alloc if out of memory.
         * \throws  FileIOError on any I/O error.
         */
        void writeCurrentTopic() const;

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
