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
/*! \internal \file
 * \brief
 * Implements classes and functions from helptopic.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_onlinehelp
 */
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

/*! \cond libapi */
void writeBasicHelpTopic(const HelpWriterContext &context,
                         const HelpTopicInterface &topic,
                         const std::string &text)
{
    const char *title = topic.title();
    if (title != NULL && title[0] != '\0')
    {
        context.writeTitle(title);
    }
    context.writeTextBlock(text);
}
//! \endcond

/********************************************************************
 * AbstractSimpleHelpTopic
 */

bool AbstractSimpleHelpTopic::hasSubTopics() const
{
    return false;
}

const HelpTopicInterface *
AbstractSimpleHelpTopic::findSubTopic(const char *name) const
{
    return NULL;
}

void AbstractSimpleHelpTopic::writeHelp(const HelpWriterContext &context) const
{
    writeBasicHelpTopic(context, *this, helpText());
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
        //! Container for mapping subtopic names to help topic objects.
        typedef std::map<std::string, HelpTopicPointer> SubTopicMap;

        /*! \brief
         * Maps subtopic names to help topic objects.
         *
         * Owns the contained subtopics.
         */
        SubTopicMap             subtopics_;
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
    return !impl_->subtopics_.empty();
}

const HelpTopicInterface *
AbstractCompositeHelpTopic::findSubTopic(const char *name) const
{
    Impl::SubTopicMap::const_iterator topic = impl_->subtopics_.find(name);
    if (topic == impl_->subtopics_.end())
    {
        return NULL;
    }
    return topic->second.get();
}

void AbstractCompositeHelpTopic::writeHelp(const HelpWriterContext &context) const
{
    writeBasicHelpTopic(context, *this, helpText());
    writeSubTopicList(context, "\nAvailable subtopics:");
}

bool
AbstractCompositeHelpTopic::writeSubTopicList(const HelpWriterContext &context,
                                              const std::string &title) const
{
    if (context.outputFormat() != eHelpOutputFormat_Console)
    {
        // TODO: Implement once the situation with Redmine issue #969 is more
        // clear.
        GMX_THROW(NotImplementedError(
                      "Subtopic listing is not implemented for this output format"));
    }
    int maxNameLength = 0;
    Impl::SubTopicMap::const_iterator topic;
    for (topic = impl_->subtopics_.begin(); topic != impl_->subtopics_.end(); ++topic)
    {
        const char *title = topic->second->title();
        if (title == NULL || title[0] == '\0')
        {
            continue;
        }
        int nameLength = static_cast<int>(topic->first.length());
        if (nameLength > maxNameLength)
        {
            maxNameLength = nameLength;
        }
    }
    if (maxNameLength == 0)
    {
        return false;
    }
    File &file = context.outputFile();
    TextTableFormatter formatter;
    formatter.addColumn(NULL, maxNameLength + 1, false);
    formatter.addColumn(NULL, 72 - maxNameLength, true);
    formatter.setFirstColumnIndent(4);
    file.writeLine(title);
    for (topic = impl_->subtopics_.begin(); topic != impl_->subtopics_.end(); ++topic)
    {
        const char *name = topic->first.c_str();
        const char *title = topic->second->title();
        if (title != NULL && title[0] != '\0')
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
    GMX_ASSERT(impl_->subtopics_.find(topic->name()) == impl_->subtopics_.end(),
               "Attempted to register a duplicate help topic name");
    impl_->subtopics_.insert(std::make_pair(std::string(topic->name()),
                                            move(topic)));
}

} // namespace gmx
