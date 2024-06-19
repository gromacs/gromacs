/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
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
#include <vector>

#include "gromacs/options/abstractoptionstorage.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/any.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

#include "options_impl.h"

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
    //! Shorthand for the internal type used to represent a section.
    typedef internal::OptionSectionImpl Section;

    //! Sets the option object to assign to.
    explicit Impl(Options* options);

    //! Returns true if a subsection has been set.
    bool inSection() const { return sectionStack_.size() > 1; }
    //! Returns the Options object for the current section.
    Section& currentSection() const { return *sectionStack_.back(); }
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
    AbstractOptionStorage* findOption(const char* name);

    //! Options object to assign to.
    Options& options_;
    //! Recognize boolean option "name" also as "noname".
    bool bAcceptBooleanNoPrefix_;
    /*! \brief
     * List of (sub)sections being assigned to.
     *
     * The first element always points to \a options_.
     */
    std::vector<Section*> sectionStack_;
    //! Current option being assigned to, or NULL if none.
    AbstractOptionStorage* currentOption_;
    /*! \brief
     * Number of values assigned so far to the current option.
     *
     * Counts the number of attempted assignments, whether they have been
     * successful or not.
     */
    int currentValueCount_;
    //! If true, a "no" prefix was given for the current boolean option.
    bool reverseBoolean_;
};

OptionsAssigner::Impl::Impl(Options* options) :
    options_(*options),
    bAcceptBooleanNoPrefix_(false),
    currentOption_(nullptr),
    currentValueCount_(0),
    reverseBoolean_(false)
{
    sectionStack_.push_back(&options_.impl_->rootSection_);
}

AbstractOptionStorage* OptionsAssigner::Impl::findOption(const char* name)
{
    GMX_RELEASE_ASSERT(currentOption_ == nullptr,
                       "Cannot search for another option while processing one");
    const Section&         section = currentSection();
    AbstractOptionStorage* option  = section.findOption(name);
    if (option == nullptr && bAcceptBooleanNoPrefix_)
    {
        if (name[0] == 'n' && name[1] == 'o')
        {
            option = section.findOption(name + 2);
            if (option != nullptr && option->isBoolean())
            {
                reverseBoolean_ = true;
            }
            else
            {
                option = nullptr;
            }
        }
    }
    return option;
}

/********************************************************************
 * OptionsAssigner
 */

OptionsAssigner::OptionsAssigner(Options* options) : impl_(new Impl(options)) {}

OptionsAssigner::~OptionsAssigner() {}

void OptionsAssigner::setAcceptBooleanNoPrefix(bool bEnabled)
{
    impl_->bAcceptBooleanNoPrefix_ = bEnabled;
}

void OptionsAssigner::start()
{
    impl_->options_.impl_->rootSection_.start();
}

void OptionsAssigner::startSection(const char* name)
{
    Impl::Section* section = impl_->currentSection().findSection(name);
    if (section == nullptr)
    {
        GMX_THROW(InvalidInputError("Unknown subsection"));
    }
    impl_->sectionStack_.push_back(section);
    section->start();
}

void OptionsAssigner::startOption(const char* name)
{
    if (!tryStartOption(name))
    {
        GMX_THROW(InvalidInputError("Unknown option " + std::string(name)));
    }
}

bool OptionsAssigner::tryStartOption(const char* name)
{
    GMX_RELEASE_ASSERT(impl_->currentOption_ == nullptr, "finishOption() not called");
    AbstractOptionStorage* option = impl_->findOption(name);
    if (option == nullptr)
    {
        return false;
    }
    option->startSet();
    impl_->currentOption_     = option;
    impl_->currentValueCount_ = 0;
    return true;
}

void OptionsAssigner::appendValue(const std::string& value)
{
    appendValue(Any(value));
}

void OptionsAssigner::appendValue(const Any& value)
{
    AbstractOptionStorage* option = impl_->currentOption_;
    GMX_RELEASE_ASSERT(option != nullptr, "startOption() not called");
    ++impl_->currentValueCount_;
    option->appendValue(value);
}

void OptionsAssigner::finishOption()
{
    AbstractOptionStorage* option = impl_->currentOption_;
    GMX_RELEASE_ASSERT(option != nullptr, "startOption() not called");
    bool bBoolReverseValue = false;
    if (option->isBoolean())
    {
        if (impl_->currentValueCount_ == 0)
        {
            // Should not throw, otherwise something is wrong.
            option->appendValue(Any::create<bool>(!impl_->reverseBoolean_));
        }
        else if (impl_->reverseBoolean_)
        {
            bBoolReverseValue = true;
        }
    }
    impl_->currentOption_  = nullptr;
    impl_->reverseBoolean_ = false;
    option->finishSet();
    if (bBoolReverseValue)
    {
        GMX_THROW(InvalidInputError("Cannot specify a value together with 'no' prefix"));
    }
}

void OptionsAssigner::finishSection()
{
    // Should only be called if we are in a subsection.
    GMX_RELEASE_ASSERT(impl_->inSection(), "startSection() not called");
    Impl::Section* section = impl_->sectionStack_.back();
    section->finish();
    impl_->sectionStack_.pop_back();
}

void OptionsAssigner::finish()
{
    GMX_RELEASE_ASSERT(impl_->currentOption_ == nullptr, "finishOption() not called");
    GMX_RELEASE_ASSERT(!impl_->inSection(), "finishSection() not called");
}

} // namespace gmx
