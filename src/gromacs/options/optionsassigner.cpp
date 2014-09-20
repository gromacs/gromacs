/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014, by the GROMACS development team, led by
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
 * Implements gmx::OptionsAssigner.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#include "gmxpre.h"

#include "optionsassigner.h"

#include <deque>

#include "gromacs/options/abstractoptionstorage.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

#include "options-impl.h"

namespace gmx
{

/********************************************************************
 * OptionsAssigner::Impl
 */

/*! \internal \brief
 * Private implementation class for OptionsAssigner.
 *
 * \ingroup module_options
 */
class OptionsAssigner::Impl
{
    public:
        //! Sets the option object to assign to.
        explicit Impl(Options *options);

        //! Returns true if a subsection has been set.
        bool inSubSection() const { return sectionStack_.size() > 1; }
        //! Returns the Options object for the current section.
        Options &currentSection() const { return *sectionStack_.back(); }
        /*! \brief
         * Finds an option by the given name.
         *
         * \param[in] name  Name of the option to look for.
         * \returns Pointer to the found option, or NULL if none found.
         *
         * This function takes into account the flags specified, and may change
         * the internal state of the assigner to match the option found.
         * If no option is found, the internal state is not modified.
         */
        AbstractOptionStorage *findOption(const char *name);

        //! Options object to assign to.
        Options                &options_;
        //! Recognize boolean option "name" also as "noname".
        bool                    bAcceptBooleanNoPrefix_;
        //! Look for options in all sections, not just the current one.
        bool                    bNoStrictSectioning_;
        /*! \brief
         * List of (sub)sections being assigned to.
         *
         * The first element always points to \a options_.
         */
        std::vector<Options *>  sectionStack_;
        //! Current option being assigned to, or NULL if none.
        AbstractOptionStorage  *currentOption_;
        /*! \brief
         * Number of values assigned so far to the current option.
         *
         * Counts the number of attempted assignments, whether they have been
         * successful or not.
         */
        int                     currentValueCount_;
        //! If true, a "no" prefix was given for the current boolean option.
        bool                    reverseBoolean_;
};

OptionsAssigner::Impl::Impl(Options *options)
    : options_(*options), bAcceptBooleanNoPrefix_(false),
      bNoStrictSectioning_(false), currentOption_(NULL),
      currentValueCount_(0), reverseBoolean_(false)
{
    sectionStack_.push_back(&options_);
}

AbstractOptionStorage *
OptionsAssigner::Impl::findOption(const char *name)
{
    GMX_RELEASE_ASSERT(currentOption_ == NULL,
                       "Cannot search for another option while processing one");
    AbstractOptionStorage *option  = NULL;
    Options               *section = NULL;
    Options               *root    = &currentSection();
    Options               *oldRoot = NULL;
    int                    upcount = 0;
    std::deque<Options *>  searchList;
    searchList.push_back(root);
    while (option == NULL && !searchList.empty())
    {
        section = searchList.front();
        option  = section->impl_->findOption(name);
        if (option == NULL && bAcceptBooleanNoPrefix_)
        {
            if (name[0] == 'n' && name[1] == 'o')
            {
                option = section->impl_->findOption(name + 2);
                if (option != NULL && option->isBoolean())
                {
                    reverseBoolean_ = true;
                }
                else
                {
                    option = NULL;
                }
            }
        }
        searchList.pop_front();
        if (bNoStrictSectioning_)
        {
            Options::Impl::SubSectionList::const_iterator i;
            for (i = section->impl_->subSections_.begin();
                 i != section->impl_->subSections_.end(); ++i)
            {
                if (*i != oldRoot)
                {
                    searchList.push_back(*i);
                }
            }
            if (searchList.empty() && root != &options_)
            {
                root = root->impl_->parent_;
                ++upcount;
                searchList.push_back(root);
            }
        }
    }
    if (bNoStrictSectioning_ && option != NULL)
    {
        while (upcount > 0)
        {
            sectionStack_.pop_back();
            --upcount;
        }
        std::vector<Options *> sections;
        while (section != &currentSection())
        {
            sections.push_back(section);
            section = section->impl_->parent_;
        }
        while (!sections.empty())
        {
            sectionStack_.push_back(sections.back());
            sections.pop_back();
        }
    }
    return option;
}

/********************************************************************
 * OptionsAssigner
 */

OptionsAssigner::OptionsAssigner(Options *options)
    : impl_(new Impl(options))
{
}

OptionsAssigner::~OptionsAssigner()
{
}

void OptionsAssigner::setAcceptBooleanNoPrefix(bool bEnabled)
{
    impl_->bAcceptBooleanNoPrefix_ = bEnabled;
}

void OptionsAssigner::setNoStrictSectioning(bool bEnabled)
{
    impl_->bNoStrictSectioning_ = bEnabled;
}

void OptionsAssigner::start()
{
    impl_->options_.impl_->startSource();
}

void OptionsAssigner::startSubSection(const char *name)
{
    Options *section = impl_->currentSection().impl_->findSubSection(name);
    if (section == NULL)
    {
        GMX_THROW(InvalidInputError("Unknown subsection"));
    }
    impl_->sectionStack_.push_back(section);
}

void OptionsAssigner::startOption(const char *name)
{
    if (!tryStartOption(name))
    {
        GMX_THROW(InvalidInputError("Unknown option " + std::string(name)));
    }
}

bool OptionsAssigner::tryStartOption(const char *name)
{
    GMX_RELEASE_ASSERT(impl_->currentOption_ == NULL, "finishOption() not called");
    AbstractOptionStorage *option = impl_->findOption(name);
    if (option == NULL)
    {
        return false;
    }
    option->startSet();
    impl_->currentOption_     = option;
    impl_->currentValueCount_ = 0;
    return true;
}

void OptionsAssigner::appendValue(const std::string &value)
{
    AbstractOptionStorage *option = impl_->currentOption_;
    GMX_RELEASE_ASSERT(option != NULL, "startOption() not called");
    ++impl_->currentValueCount_;
    option->appendValue(value);
}

void OptionsAssigner::finishOption()
{
    AbstractOptionStorage *option = impl_->currentOption_;
    GMX_RELEASE_ASSERT(option != NULL, "startOption() not called");
    bool                   bBoolReverseValue = false;
    if (option->isBoolean())
    {
        if (impl_->currentValueCount_ == 0)
        {
            // Should not throw, otherwise something is wrong.
            // TODO: Get rid of the hard-coded values.
            option->appendValue(impl_->reverseBoolean_ ? "0" : "1");
        }
        else if (impl_->reverseBoolean_)
        {
            bBoolReverseValue = true;
        }
    }
    impl_->currentOption_  = NULL;
    impl_->reverseBoolean_ = false;
    option->finishSet();
    if (bBoolReverseValue)
    {
        GMX_THROW(InvalidInputError("Cannot specify a value together with 'no' prefix"));
    }
}

void OptionsAssigner::finishSubSection()
{
    // Should only be called if we are in a subsection.
    GMX_RELEASE_ASSERT(impl_->inSubSection(), "startSubSection() not called");
    impl_->sectionStack_.pop_back();
}

void OptionsAssigner::finish()
{
    GMX_RELEASE_ASSERT(impl_->currentOption_ == NULL, "finishOption() not called");
    if (impl_->bNoStrictSectioning_)
    {
        while (impl_->inSubSection())
        {
            finishSubSection();
        }
    }
    GMX_RELEASE_ASSERT(!impl_->inSubSection(), "finishSubSection() not called");
}

} // namespace gmx
