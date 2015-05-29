/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014,2015, by the GROMACS development team, led by
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
 * Declares gmx::HelpTopicInterface.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_onlinehelp
 */
#ifndef GMX_ONLINEHELP_HELPTOPICINTERFACE_H
#define GMX_ONLINEHELP_HELPTOPICINTERFACE_H

#include "gromacs/utility/uniqueptr.h"

namespace gmx
{

class HelpWriterContext;

/*! \libinternal \brief
 * Provides a single online help topic.
 *
 * Implementations of these methods should not throw, except that writeHelp()
 * is allowed to throw on out-of-memory or I/O errors since those it cannot
 * avoid.
 *
 * Header helptopic.h contains classes that implement this interface and make
 * it simple to write concrete help topic classes.
 *
 * \inlibraryapi
 * \ingroup module_onlinehelp
 */
class HelpTopicInterface
{
    public:
        virtual ~HelpTopicInterface() {}

        /*! \brief
         * Returns the name of the topic.
         *
         * This should be a single lowercase word, used to identify the topic.
         * It is not used for the root of the help topic tree.
         */
        virtual const char *name() const = 0;
        /*! \brief
         * Returns a title for the topic.
         *
         * May return NULL, in which case the topic is omitted from normal
         * subtopic lists and no title is printed by the methods provided in
         * helptopic.h.
         */
        virtual const char *title() const = 0;

        //! Returns whether the topic has any subtopics.
        virtual bool hasSubTopics() const = 0;
        /*! \brief
         * Finds a subtopic by name.
         *
         * \param[in] name  Name of subtopic to find.
         * \returns   Pointer to the found subtopic, or NULL if matching topic
         *      is not found.
         */
        virtual const HelpTopicInterface *findSubTopic(const char *name) const = 0;

        /*! \brief
         * Prints the help text for this topic.
         *
         * \param[in] context  Context object for writing the help.
         * \throws    std::bad_alloc if out of memory.
         * \throws    FileIOError on any I/O error.
         */
        virtual void writeHelp(const HelpWriterContext &context) const = 0;
};

//! Smart pointer type to manage a HelpTopicInterface object.
typedef gmx_unique_ptr<HelpTopicInterface>::type HelpTopicPointer;

} // namespace gmx

#endif
