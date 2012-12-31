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
 * Declares helper classes and functions for implementing
 * gmx::HelpTopicInterface.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_onlinehelp
 */
#ifndef GMX_ONLINEHELP_HELPTOPIC_H
#define GMX_ONLINEHELP_HELPTOPIC_H

#include "../utility/common.h"
#include "../utility/stringutil.h"
#include "../utility/uniqueptr.h"

#include "helptopicinterface.h"

namespace gmx
{

/*! \cond libapi */
/*! \libinternal \brief
 * Helper for writing simple help text.
 *
 * \param[in] context Context for writing the help.
 * \param[in] topic   Topic to write the help for (used for title).
 * \param[in] text    Text to write for the topic.
 * \throws    std::bad_alloc if out of memory.
 * \throws    FileIOError on any I/O error.
 *
 * Formats basic help by writing a title (obtained from \p topic), followed by
 * \p text with markup substituted and lines properly wrapped.
 *
 * \inlibraryapi
 */
void writeBasicHelpTopic(const HelpWriterContext &context,
                         const HelpTopicInterface &topic,
                         const std::string &text);
//! \endcond

/*! \libinternal \brief
 * Abstract base class for help topics that have simple text and no subtopics.
 *
 * This class implements subtopic-related methods from HelpTopicInterface such
 * that there are no subtopics.  writeHelp() is also implemented such that it
 * uses writeBasicHelpTopic() to write out the text returned by a new virtual
 * method helpText().
 *
 * \see SimpleHelpTopic
 *
 * \inlibraryapi
 * \ingroup module_onlinehelp
 */
class AbstractSimpleHelpTopic : public HelpTopicInterface
{
    public:
        virtual const char *name() const = 0;
        virtual const char *title() const = 0;

        virtual bool hasSubTopics() const;
        virtual const HelpTopicInterface *findSubTopic(const char *name) const;

        virtual void writeHelp(const HelpWriterContext &context) const;

    protected:
        /*! \brief
         * Returns the help text for this topic.
         *
         * writeHelp() calls this method to obtain the actual text to format
         * for the topic.  Markup substitution etc. is done automatically by
         * writeHelp().
         */
        virtual std::string helpText() const = 0;
};

/*! \libinternal \brief
 * Abstract base class for help topics that have simple text and subtopics.
 *
 * This class implements an internal container for subtopics and provides
 * public methods for adding subtopics (as HelpTopicInterface objects).
 * Subtopic-related methods from HelpTopicInterface are implemented to access
 * the internal container.  writeHelp() is also implemented such that it
 * uses writeBasicHelpTopic() to write out the text returned by a new virtual
 * method helpText(), and a list of subtopics is written after the actual text.
 *
 * \see CompositeHelpTopic
 *
 * \inlibraryapi
 * \ingroup module_onlinehelp
 */
class AbstractCompositeHelpTopic : public HelpTopicInterface
{
    public:
        AbstractCompositeHelpTopic();
        virtual ~AbstractCompositeHelpTopic();

        virtual const char *name() const = 0;
        virtual const char *title() const = 0;

        virtual bool hasSubTopics() const;
        virtual const HelpTopicInterface *findSubTopic(const char *name) const;

        virtual void writeHelp(const HelpWriterContext &context) const;

        /*! \brief
         * Adds a given topic as a subtopic of this topic.
         *
         * \param   topic  Topis to add.
         * \throws  std::bad_alloc if out of memory.
         *
         * This topic takes ownership of the object.
         *
         * \see registerSubTopic()
         */
        void addSubTopic(HelpTopicPointer topic);
        /*! \brief
         * Registers a subtopic of a certain type to this topic.
         *
         * \tparam  Topic  Type of topic to register.
         * \throws  std::bad_alloc if out of memory.
         *
         * \p Topic must be default-constructible and implement
         * HelpTopicInterface.
         *
         * This method is provided as a convenient alternative to addSubTopic()
         * for cases where each topic is implemented by a different type
         * (which is a common case outside unit tests).
         */
        template <class Topic>
        void registerSubTopic()
        {
            addSubTopic(HelpTopicPointer(new Topic));
        }

    protected:
        //! \copydoc gmx::AbstractSimpleHelpTopic::helpText()
        virtual std::string helpText() const = 0;

        /*! \brief
         * Writes the list of subtopics.
         *
         * \param[in] context Context for writing the help.
         * \param[in] title  Title for the written list.
         * \returns   true if anything was printed.
         * \throws    std::bad_alloc if out of memory.
         * \throws    FileIOError on any I/O error.
         *
         * Subtopics with empty titles are skipped from the list.
         * If there would be no subtopics in the list, \p title is not printed
         * either.
         *
         * This method is provided for cases where helpText() does not provide
         * the needed flexibility and the derived class needs to override
         * writeHelp().  This method can then be called to print the same
         * subtopic list that is printed by the default writeHelp()
         * implementation.
         */
        bool writeSubTopicList(const HelpWriterContext &context,
                               const std::string &title) const;

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

/*! \cond libapi */
/*! \libinternal \brief
 * Smart pointer type to manage a AbstractCompositeHelpTopic object.
 *
 * \inlibraryapi
 */
typedef gmx_unique_ptr<AbstractCompositeHelpTopic>::type
CompositeHelpTopicPointer;
//! \endcond

/*! \libinternal \brief
 * Template for simple implementation of AbstractSimpleHelpTopic.
 *
 * \tparam HelpText Struct that defines the data for the topic.
 *
 * \p HelpText should have public static members \c "const char name[]",
 * \c "const char title[]" and \c "const char *const text[]".
 *
 * Typical use:
 * \code
struct ExampleHelpText
{
    static const char name[];
    static const char title[];
    static const char *const text[];
};

const char ExampleHelpText::name[] = "example";
const char ExampleHelpText::title[] =
    "Example title";
const char *const ExampleHelpText::text[] = {
    "Text for the topic.",
    "More text for the topic."
};

typedef SimpleHelpTopic<ExampleHelpText> ExampleHelpTopic;
 * \endcode
 *
 * \inlibraryapi
 * \ingroup module_onlinehelp
 */
template <class HelpText>
class SimpleHelpTopic : public AbstractSimpleHelpTopic
{
    public:
        virtual const char *name() const
        {
            return HelpText::name;
        }
        virtual const char *title() const
        {
            return HelpText::title;
        }

    protected:
        virtual std::string helpText() const
        {
            return concatenateStrings(HelpText::text);
        }
};

/*! \libinternal \brief
 * Template for simple implementation of AbstractCompositeHelpTopic.
 *
 * \tparam HelpText Struct that defines the data for the topic.
 *
 * Used similarly to SimpleHelpTopic.
 * \p HelpText should satisfy the same criteria as for SimpleHelpTopic.
 *
 * \see SimpleHelpTopic
 *
 * \inlibraryapi
 * \ingroup module_onlinehelp
 */
template <class HelpText>
class CompositeHelpTopic : public AbstractCompositeHelpTopic
{
    public:
        virtual const char *name() const
        {
            return HelpText::name;
        }
        virtual const char *title() const
        {
            return HelpText::title;
        }

    protected:
        virtual std::string helpText() const
        {
            return concatenateStrings(HelpText::text);
        }
};

} // namespace gmx

#endif
