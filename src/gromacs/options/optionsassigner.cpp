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
 * Implements gmx::OptionsAssigner.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_options
 */
#include "gromacs/options/optionsassigner.h"

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
        //! Possible flags for controlling assignment behavior.
        enum Flag
        {
            //! Recognize boolean option "name" also as "noname".
            efAcceptBooleanNoPrefix     = 1<<0,
            //! Look for options in all sections, not just the current one.
            efNoStrictSectioning        = 1<<1,
        };
        //! Sets the option object to assign to.
        Impl(Options *options);

        //! Sets or clears the given flag.
        void setFlag(Flag flag, bool bSet);

        //! Returns true if the given flag is set.
        bool hasFlag(Flag flag) const { return _flags & flag; }
        //! Returns true if a subsection has been set.
        bool inSubSection() const { return _sectionStack.size() > 1; }
        //! Returns the Options object for the current section.
        Options &currentSection() const { return *_sectionStack.back(); }
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
        Options                &_options;
        //! Flags that control assignment behavior.
        unsigned long           _flags;
        /*! \brief
         * List of (sub)sections being assigned to.
         *
         * The first element always points to \a _options.
         */
        std::vector<Options *>  _sectionStack;
        //! Current option being assigned to, or NULL if none.
        AbstractOptionStorage  *_currentOption;
        //! Number of values assigned so far to the current option.
        int                     _currentValueCount;
        //! If true, a "no" prefix was given for the current boolean option.
        bool                    _reverseBoolean;
};

OptionsAssigner::Impl::Impl(Options *options)
    : _options(*options), _flags(0), _currentOption(NULL),
      _currentValueCount(0), _reverseBoolean(false)
{
    _sectionStack.push_back(&_options);
}

void OptionsAssigner::Impl::setFlag(OptionsAssigner::Impl::Flag flag, bool bSet)
{
    if (bSet)
    {
        _flags |= flag;
    }
    else
    {
        _flags &= ~flag;
    }
}

AbstractOptionStorage *
OptionsAssigner::Impl::findOption(const char *name)
{
    GMX_RELEASE_ASSERT(_currentOption == NULL,
                       "Cannot search for another option while processing one");
    AbstractOptionStorage *option = NULL;
    Options *section = NULL;
    Options *root = &currentSection();
    Options *oldRoot = NULL;
    int      upcount = 0;
    std::deque<Options *> searchList;
    searchList.push_back(root);
    while (option == NULL && !searchList.empty())
    {
        section = searchList.front();
        option = section->_impl->findOption(name);
        if (option == NULL && hasFlag(efAcceptBooleanNoPrefix))
        {
            if (name[0] == 'n' && name[1] == 'o')
            {
                option = section->_impl->findOption(name + 2);
                if (option != NULL && option->isBoolean())
                {
                    _reverseBoolean = true;
                }
                else
                {
                    option = NULL;
                }
            }
        }
        searchList.pop_front();
        if (hasFlag(efNoStrictSectioning))
        {
            Options::Impl::SubSectionList::const_iterator i;
            for (i = section->_impl->_subSections.begin();
                 i != section->_impl->_subSections.end(); ++i)
            {
                if (*i != oldRoot)
                {
                    searchList.push_back(*i);
                }
            }
            if (searchList.empty() && root != &_options)
            {
                root = root->_impl->_parent;
                ++upcount;
                searchList.push_back(root);
            }
        }
    }
    if (hasFlag(efNoStrictSectioning) && option != NULL)
    {
        while (upcount > 0)
        {
            _sectionStack.pop_back();
            --upcount;
        }
        std::vector<Options *> sections;
        while (section != &currentSection())
        {
            sections.push_back(section);
            section = section->_impl->_parent;
        }
        while (!sections.empty())
        {
            _sectionStack.push_back(sections.back());
            sections.pop_back();
        }
    }
    return option;
}

/********************************************************************
 * OptionsAssigner
 */

OptionsAssigner::OptionsAssigner(Options *options)
    : _impl(new Impl(options))
{
}

OptionsAssigner::~OptionsAssigner()
{
}

void OptionsAssigner::setAcceptBooleanNoPrefix(bool enabled)
{
    _impl->setFlag(Impl::efAcceptBooleanNoPrefix, enabled);
}

void OptionsAssigner::setNoStrictSectioning(bool enabled)
{
    _impl->setFlag(Impl::efNoStrictSectioning, enabled);
}

void OptionsAssigner::start()
{
    _impl->_options._impl->startSource();
}

void OptionsAssigner::startSubSection(const char *name)
{
    Options *section = _impl->currentSection()._impl->findSubSection(name);
    if (section == NULL)
    {
        GMX_THROW(InvalidInputError("Unknown subsection"));
    }
    _impl->_sectionStack.push_back(section);
}

void OptionsAssigner::startOption(const char *name)
{
    GMX_RELEASE_ASSERT(_impl->_currentOption == NULL, "finishOption() not called");
    AbstractOptionStorage *option = _impl->findOption(name);
    if (option == NULL)
    {
        GMX_THROW(InvalidInputError("Unknown option"));
    }
    option->startSet();
    _impl->_currentOption = option;
    _impl->_currentValueCount = 0;
}

void OptionsAssigner::appendValue(const std::string &value)
{
    AbstractOptionStorage *option = _impl->_currentOption;
    GMX_RELEASE_ASSERT(option != NULL, "startOption() not called");
    // Does not count correctly, but the actual count is not really used.
    // TODO: Rename the variable to better reflect the usage.
    ++_impl->_currentValueCount;
    option->appendValue(value);
}

void OptionsAssigner::finishOption()
{
    AbstractOptionStorage *option = _impl->_currentOption;
    GMX_RELEASE_ASSERT(option != NULL, "startOption() not called");
    bool bBoolReverseValue = false;
    if (option->isBoolean())
    {
        if (_impl->_currentValueCount == 0)
        {
            // Should not throw, otherwise something is wrong.
            // TODO: Get rid of the hard-coded values.
            option->appendValue(_impl->_reverseBoolean ? "0" : "1");
        }
        else if (_impl->_reverseBoolean)
        {
            bBoolReverseValue = true;
        }
    }
    _impl->_currentOption = NULL;
    _impl->_reverseBoolean = false;
    option->finishSet();
    if (bBoolReverseValue)
    {
        GMX_THROW(InvalidInputError("Cannot specify a value together with 'no' prefix"));
    }
}

void OptionsAssigner::finishSubSection()
{
    // Should only be called if we are in a subsection.
    GMX_RELEASE_ASSERT(_impl->inSubSection(), "startSubSection() not called");
    _impl->_sectionStack.pop_back();
}

void OptionsAssigner::finish()
{
    GMX_RELEASE_ASSERT(_impl->_currentOption == NULL, "finishOption() not called");
    if (_impl->hasFlag(Impl::efNoStrictSectioning))
    {
        while (_impl->inSubSection())
        {
            finishSubSection();
        }
    }
    GMX_RELEASE_ASSERT(!_impl->inSubSection(), "finishSubSection() not called");
}

} // namespace gmx
