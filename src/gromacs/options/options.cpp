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
 * Implements gmx::Options.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#include "gmxpre.h"

#include "gromacs/options/options.h"

#include <map>
#include <memory>
#include <string>
#include <utility>

#include "gromacs/options/abstractoption.h"
#include "gromacs/options/abstractoptionstorage.h"
#include "gromacs/options/abstractsection.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/ioptionscontainerwithsections.h"
#include "gromacs/options/isectionstorage.h"
#include "gromacs/options/optionmanagercontainer.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "options_impl.h"

namespace gmx
{

/********************************************************************
 * IOptionManager
 */

IOptionManager::~IOptionManager() {}

/********************************************************************
 * IOptionsContainer
 */

IOptionsContainer::~IOptionsContainer() {}

/********************************************************************
 * IOptionsContainerWithSections
 */

IOptionsContainerWithSections::~IOptionsContainerWithSections() {}

/********************************************************************
 * IOptionSectionStorage
 */

IOptionSectionStorage::~IOptionSectionStorage() {}

/********************************************************************
 * OptionsImpl
 */

namespace internal
{

OptionsImpl::OptionsImpl() : rootSection_(managers_, nullptr, "") {}

/********************************************************************
 * OptionSectionImpl
 */

OptionSectionImpl* OptionSectionImpl::addSectionImpl(const AbstractOptionSection& section)
{
    const char* name = section.name_;
    // Make sure that there are no duplicate sections.
    GMX_RELEASE_ASSERT(findSection(name) == nullptr, "Duplicate subsection name");
    std::unique_ptr<IOptionSectionStorage> storage(section.createStorage());
    subsections_.push_back(std::make_unique<OptionSectionImpl>(managers_, std::move(storage), name));
    return subsections_.back().get();
}

IOptionsContainer& OptionSectionImpl::addGroup()
{
    return rootGroup_.addGroup();
}

OptionInfo* OptionSectionImpl::addOptionImpl(const AbstractOption& settings)
{
    return rootGroup_.addOptionImpl(settings);
}

OptionSectionImpl* OptionSectionImpl::findSection(const char* name) const
{
    for (const auto& section : subsections_)
    {
        if (section->name_ == name)
        {
            return section.get();
        }
    }
    return nullptr;
}

AbstractOptionStorage* OptionSectionImpl::findOption(const char* name) const
{
    OptionMap::const_iterator i = optionMap_.find(name);
    if (i == optionMap_.end())
    {
        return nullptr;
    }
    return i->second.get();
}

void OptionSectionImpl::start()
{
    for (const auto& entry : optionMap_)
    {
        entry.second->startSource();
    }
    if (storage_ != nullptr)
    {
        if (!storageInitialized_)
        {
            storage_->initStorage();
            storageInitialized_ = true;
        }
        storage_->startSection();
    }
}

void OptionSectionImpl::finish()
{
    // TODO: Consider how to customize these error messages based on context.
    ExceptionInitializer errors("Invalid input values");
    for (const auto& entry : optionMap_)
    {
        AbstractOptionStorage& option = *entry.second;
        try
        {
            option.finish();
        }
        catch (UserInputError& ex)
        {
            ex.prependContext("In option " + option.name());
            errors.addCurrentExceptionAsNested();
        }
    }
    if (errors.hasNestedExceptions())
    {
        // TODO: This exception type may not always be appropriate.
        GMX_THROW(InvalidInputError(errors));
    }
    if (storage_ != nullptr)
    {
        storage_->finishSection();
    }
}

/********************************************************************
 * OptionSectionImpl::Group
 */

IOptionsContainer& OptionSectionImpl::Group::addGroup()
{
    subgroups_.emplace_back(parent_);
    return subgroups_.back();
}

OptionInfo* OptionSectionImpl::Group::addOptionImpl(const AbstractOption& settings)
{
    OptionSectionImpl::AbstractOptionStoragePointer option(settings.createStorage(parent_->managers_));
    options_.reserve(options_.size() + 1);
    auto insertionResult = parent_->optionMap_.insert(std::make_pair(option->name(), std::move(option)));
    if (!insertionResult.second)
    {
        const std::string& name = insertionResult.first->second->name();
        GMX_THROW(APIError("Duplicate option: " + name));
    }
    AbstractOptionStorage& insertedOption = *insertionResult.first->second;
    options_.push_back(&insertedOption);
    return &insertedOption.optionInfo();
}

} // namespace internal

using internal::OptionsImpl;

/********************************************************************
 * Options
 */

Options::Options() : impl_(new OptionsImpl) {}

Options::~Options() {}


void Options::addManager(IOptionManager* manager)
{
    // This ensures that all options see the same set of managers.
    GMX_RELEASE_ASSERT(impl_->rootSection_.optionMap_.empty(),
                       "Can only add a manager before options");
    // This check could be relaxed if we instead checked that the subsections
    // do not have options.
    GMX_RELEASE_ASSERT(impl_->rootSection_.subsections_.empty(),
                       "Can only add a manager before subsections");
    impl_->managers_.add(manager);
}

internal::OptionSectionImpl* Options::addSectionImpl(const AbstractOptionSection& section)
{
    return impl_->rootSection_.addSectionImpl(section);
}

IOptionsContainer& Options::addGroup()
{
    return impl_->rootSection_.addGroup();
}

OptionInfo* Options::addOptionImpl(const AbstractOption& settings)
{
    return impl_->rootSection_.addOptionImpl(settings);
}

OptionSectionInfo& Options::rootSection()
{
    return impl_->rootSection_.info();
}

const OptionSectionInfo& Options::rootSection() const
{
    return impl_->rootSection_.info();
}

void Options::finish()
{
    impl_->rootSection_.finish();
}

} // namespace gmx
