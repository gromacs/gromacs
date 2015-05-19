/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
 * Implements classes and functions from helptopic.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_onlinehelp
 */
#include "gmxpre.h"

#include "helptopic.h"

#include <map>
#include <utility>

#include "gromacs/onlinehelp/helpformat.h"
#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

/********************************************************************
 * AbstractSimpleHelpTopic
 */

bool AbstractSimpleHelpTopic::hasSubTopics() const
{
    return false;
}

const HelpTopicInterface *
AbstractSimpleHelpTopic::findSubTopic(const char * /* name */) const
{
    return NULL;
}

void AbstractSimpleHelpTopic::writeHelp(const HelpWriterContext &context) const
{
    context.writeTextBlock(helpText());
}

/********************************************************************
 * AbstractCompositeHelpTopic::Impl
 */

/*! \internal \brief
 * Private implementation class for AbstractCompositeHelpTopic.
 *
 * \ingroup module_onlinehelp
 */
class AbstractCompositeHelpTopic::Impl
{
    public:
        //! Container for subtopics.
        typedef std::vector<HelpTopicPointer> SubTopicList;
        //! Container for mapping subtopic names to help topic objects.
        typedef std::map<std::string, const HelpTopicInterface *> SubTopicMap;

        /*! \brief
         * Subtopics in the order they were added.
         *
         * Owns the contained subtopics.
         */
        SubTopicList            subTopics_;
        /*! \brief
         * Maps subtopic names to help topic objects.
         *
         * Points to objects in the \a subTopics_ map.
         */
        SubTopicMap             subTopicMap_;
};

/********************************************************************
 * AbstractCompositeHelpTopic
 */

AbstractCompositeHelpTopic::AbstractCompositeHelpTopic()
    : impl_(new Impl)
{
}

AbstractCompositeHelpTopic::~AbstractCompositeHelpTopic()
{
}

bool AbstractCompositeHelpTopic::hasSubTopics() const
{
    return !impl_->subTopics_.empty();
}

const HelpTopicInterface *
AbstractCompositeHelpTopic::findSubTopic(const char *name) const
{
    Impl::SubTopicMap::const_iterator topic = impl_->subTopicMap_.find(name);
    if (topic == impl_->subTopicMap_.end())
    {
        return NULL;
    }
    return topic->second;
}

void AbstractCompositeHelpTopic::writeHelp(const HelpWriterContext &context) const
{
    context.writeTextBlock(helpText());
    writeSubTopicList(context, "\nAvailable subtopics:");
}

bool
AbstractCompositeHelpTopic::writeSubTopicList(const HelpWriterContext &context,
                                              const std::string       &title) const
{
    if (context.outputFormat() != eHelpOutputFormat_Console)
    {
        Impl::SubTopicList::const_iterator topic;
        for (topic = impl_->subTopics_.begin(); topic != impl_->subTopics_.end(); ++topic)
        {
            const char *const title = (*topic)->title();
            if (!isNullOrEmpty(title))
            {
                context.outputFile().writeLine();
                HelpWriterContext subContext(context);
                subContext.enterSubSection(title);
                (*topic)->writeHelp(subContext);
            }
        }
        return true;
    }
    int maxNameLength = 0;
    Impl::SubTopicMap::const_iterator topic;
    for (topic = impl_->subTopicMap_.begin(); topic != impl_->subTopicMap_.end(); ++topic)
    {
        const char *const title = topic->second->title();
        if (!isNullOrEmpty(title))
        {
            int nameLength = static_cast<int>(topic->first.length());
            if (nameLength > maxNameLength)
            {
                maxNameLength = nameLength;
            }
        }
    }
    if (maxNameLength == 0)
    {
        return false;
    }
    File              &file = context.outputFile();
    TextTableFormatter formatter;
    formatter.addColumn(NULL, maxNameLength + 1, false);
    formatter.addColumn(NULL, 72 - maxNameLength, true);
    formatter.setFirstColumnIndent(4);
    file.writeLine(title);
    for (topic = impl_->subTopicMap_.begin(); topic != impl_->subTopicMap_.end(); ++topic)
    {
        const char *const name  = topic->first.c_str();
        const char *const title = topic->second->title();
        if (!isNullOrEmpty(title))
        {
            formatter.clear();
            formatter.addColumnLine(0, name);
            formatter.addColumnLine(1, title);
            file.writeString(formatter.formatRow());
        }
    }
    return true;
}

void AbstractCompositeHelpTopic::addSubTopic(HelpTopicPointer topic)
{
    GMX_ASSERT(impl_->subTopicMap_.find(topic->name()) == impl_->subTopicMap_.end(),
               "Attempted to register a duplicate help topic name");
    const HelpTopicInterface *topicPtr = topic.get();
    impl_->subTopics_.reserve(impl_->subTopics_.size() + 1);
    impl_->subTopicMap_.insert(std::make_pair(std::string(topicPtr->name()), topicPtr));
    impl_->subTopics_.push_back(move(topic));
}

} // namespace gmx
