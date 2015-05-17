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
 * Declares helper classes for implementing gmx::HelpTopicInterface.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_onlinehelp
 */
#ifndef GMX_ONLINEHELP_HELPTOPIC_H
#define GMX_ONLINEHELP_HELPTOPIC_H

#include "gromacs/onlinehelp/helptopicinterface.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/uniqueptr.h"

namespace gmx
{

/*! \libinternal \brief
 * Abstract base class for help topics that have simple text and no subtopics.
 *
 * This class implements subtopic-related methods from HelpTopicInterface such
 * that there are no subtopics.  writeHelp() is also implemented such that it
 * uses HelpTopicContext::writeTextBlock() to write out the text returned by a
 * new virtual method helpText().
 *
 * \see SimpleHelpTopic
 *
 * \inlibraryapi
 * \ingroup module_onlinehelp
 */
class AbstractSimpleHelpTopic : public HelpTopicInterface
{
    public:
        virtual const char *name() const  = 0;
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
 * uses HelpTopicContext::writeTextBlock() to write out the text returned by a
 * new virtual method helpText(), and a list of subtopics is written after the
 * actual text.
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

        virtual const char *name() const  = 0;
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
                               const std::string       &title) const;

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
            return joinStrings(HelpText::text, "\n");
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
            return joinStrings(HelpText::text, "\n");
        }
};

} // namespace gmx

#endif
