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
/*! \internal \file
 * \brief
 * Implements gmx::HelpManager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_onlinehelp
 */
#include "gmxpre.h"

#include "helpmanager.h"

#include <string>
#include <vector>

#include "gromacs/onlinehelp/helptopicinterface.h"
#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

/********************************************************************
 * HelpManager::Impl
 */

/*! \internal \brief
 * Private implementation class for HelpManager.
 *
 * \ingroup module_onlinehelp
 */
class HelpManager::Impl
{
    public:
        //! Container type for keeping the stack of active topics.
        typedef std::vector<const HelpTopicInterface *> TopicStack;

        //! Initializes a new manager with the given context.
        explicit Impl(const HelpWriterContext &context)
            : rootContext_(context)
        {
        }

        //! Whether the active topic is the root topic.
        bool isAtRootTopic() const { return topicStack_.size() == 1; }
        //! Returns the active topic.
        const HelpTopicInterface &currentTopic() const
        {
            return *topicStack_.back();
        }
        //! Formats the active topic as a string, including its parent topics.
        std::string currentTopicAsString() const;

        //! Context with which the manager was initialized.
        const HelpWriterContext &rootContext_;
        /*! \brief
         * Stack of active topics.
         *
         * The first item is always the root topic, and each item is a subtopic
         * of the preceding item.  The last item is the currently active topic.
         */
        TopicStack               topicStack_;
};

std::string HelpManager::Impl::currentTopicAsString() const
{
    std::string                result;
    TopicStack::const_iterator topic;
    for (topic = topicStack_.begin() + 1; topic != topicStack_.end(); ++topic)
    {
        if (!result.empty())
        {
            result.append(" ");
        }
        result.append((*topic)->name());
    }
    return result;
}

/********************************************************************
 * HelpManager
 */

HelpManager::HelpManager(const HelpTopicInterface &rootTopic,
                         const HelpWriterContext  &context)
    : impl_(new Impl(context))
{
    impl_->topicStack_.push_back(&rootTopic);
}

HelpManager::~HelpManager()
{
}

void HelpManager::enterTopic(const char *name)
{
    const HelpTopicInterface &topic = impl_->currentTopic();
    if (!topic.hasSubTopics())
    {
        GMX_THROW(InvalidInputError(
                          formatString("Help topic '%s' has no subtopics",
                                       impl_->currentTopicAsString().c_str())));
    }
    const HelpTopicInterface *newTopic = topic.findSubTopic(name);
    if (newTopic == NULL)
    {
        if (impl_->isAtRootTopic())
        {
            GMX_THROW(InvalidInputError(
                              formatString("No help available for '%s'", name)));
        }
        else
        {
            GMX_THROW(InvalidInputError(
                              formatString("Help topic '%s' has no subtopic '%s'",
                                           impl_->currentTopicAsString().c_str(), name)));
        }
    }
    impl_->topicStack_.push_back(newTopic);
}

void HelpManager::enterTopic(const std::string &name)
{
    enterTopic(name.c_str());
}

void HelpManager::writeCurrentTopic() const
{
    const HelpTopicInterface &topic = impl_->currentTopic();
    const char               *title = topic.title();
    HelpWriterContext         context(impl_->rootContext_);
    context.enterSubSection(title != NULL ? title : "");
    topic.writeHelp(context);
}

} // namespace gmx
