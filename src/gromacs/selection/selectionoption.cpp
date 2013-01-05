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
 * Implements classes in selectionoption.h and selectionoptionstorage.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_selection
 */
#include "selectionoption.h"
#include "selectionfileoption.h"
#include "selectionoptionstorage.h"
#include "selectionfileoptionstorage.h"

#include <string>

#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoptionmanager.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/messagestringcollector.h"

namespace gmx
{

/********************************************************************
 * SelectionOptionStorage
 */

SelectionOptionStorage::SelectionOptionStorage(const SelectionOption &settings)
    : MyBase(settings, OptionFlags() | efOption_NoDefaultValue
             | efOption_DontCheckMinimumCount),
      info_(this), manager_(NULL), selectionFlags_(settings.selectionFlags_)
{
    GMX_RELEASE_ASSERT(!hasFlag(efOption_MultipleTimes),
                       "allowMultiple() is not supported for selection options");
}


void SelectionOptionStorage::setManager(SelectionOptionManager *manager)
{
    GMX_RELEASE_ASSERT(manager_ == NULL || manager_ == manager,
                       "Manager cannot be changed once set");
    if (manager_ == NULL)
    {
        manager->registerOption(this);
        manager_ = manager;
    }
}


std::string SelectionOptionStorage::formatSingleValue(const Selection &value) const
{
    return value.selectionText();
}


void SelectionOptionStorage::addSelections(
    const SelectionList &selections,
    bool                 bFullValue)
{
    if (bFullValue && selections.size() < static_cast<size_t>(minValueCount()))
    {
        GMX_THROW(InvalidInputError("Too few selections provided"));
    }
    if (bFullValue)
    {
        clearSet();
    }
    SelectionList::const_iterator i;
    for (i = selections.begin(); i != selections.end(); ++i)
    {
        // TODO: Having this check in the parser would make interactive input
        // behave better.
        if (selectionFlags_.test(efSelection_OnlyStatic) && i->isDynamic())
        {
            GMX_THROW(InvalidInputError("Dynamic selections not supported"));
        }
        addValue(*i);
    }
    if (bFullValue)
    {
        commitValues();
        markAsSet();
    }
}


void SelectionOptionStorage::convertValue(const std::string &value)
{
    GMX_RELEASE_ASSERT(manager_ != NULL, "Manager is not set");

    manager_->convertOptionValue(this, value);
}

void SelectionOptionStorage::processSetValues(ValueList *values)
{
    GMX_RELEASE_ASSERT(manager_ != NULL, "Manager is not set");

    if (values->size() == 0)
    {
        manager_->requestOptionDelayedParsing(this);
    }
    else if (values->size() < static_cast<size_t>(minValueCount()))
    {
        GMX_THROW(InvalidInputError("Too few (valid) values provided"));
    }
    ValueList::iterator i;
    for (i = values->begin(); i != values->end(); ++i)
    {
        i->data().setFlags(selectionFlags_);
    }
}

void SelectionOptionStorage::processAll()
{
    if (isRequired() && !isSet())
    {
        GMX_RELEASE_ASSERT(manager_ != NULL, "Manager is not set");

        manager_->requestOptionDelayedParsing(this);
        markAsSet();
    }
}

void SelectionOptionStorage::setAllowedValueCount(int count)
{
    // TODO: It should be possible to have strong exception safety here.
    MessageStringCollector errors;
    errors.startContext("In option '" + name() + "'");
    if (count >= 0)
    {
        // Should not throw because efOption_DontCheckMinimumCount is set.
        setMinValueCount(count);
        if (valueCount() > 0 && valueCount() < count)
        {
            errors.append("Too few (valid) values provided");
        }
    }
    try
    {
        setMaxValueCount(count);
    }
    catch (const UserInputError &ex)
    {
        errors.append(ex.what());
    }
    errors.finishContext();
    if (!errors.isEmpty())
    {
        GMX_THROW(InvalidInputError(errors.toString()));
    }
}

void SelectionOptionStorage::setSelectionFlag(SelectionFlag flag, bool bSet)
{
    ValueList::iterator i;
    for (i = values().begin(); i != values().end(); ++i)
    {
        if (flag == efSelection_OnlyStatic && bSet && i->isDynamic())
        {
            MessageStringCollector errors;
            errors.startContext("In option '" + name() + "'");
            errors.append("Dynamic selections not supported");
            errors.finishContext();
            GMX_THROW(InvalidInputError(errors.toString()));
        }
    }
    selectionFlags_.set(flag, bSet);
    for (i = values().begin(); i != values().end(); ++i)
    {
        i->data().setFlags(selectionFlags_);
    }
}


/********************************************************************
 * SelectionOptionInfo
 */

SelectionOptionInfo::SelectionOptionInfo(SelectionOptionStorage *option)
    : OptionInfo(option)
{
}

SelectionOptionStorage &SelectionOptionInfo::option()
{
    return static_cast<SelectionOptionStorage &>(OptionInfo::option());
}

const SelectionOptionStorage &SelectionOptionInfo::option() const
{
    return static_cast<const SelectionOptionStorage &>(OptionInfo::option());
}

void SelectionOptionInfo::setManager(SelectionOptionManager *manager)
{
    option().setManager(manager);
}

void SelectionOptionInfo::setValueCount(int count)
{
    option().setAllowedValueCount(count);
}

void SelectionOptionInfo::setEvaluateVelocities(bool bEnabled)
{
    option().setSelectionFlag(efSelection_EvaluateVelocities, bEnabled);
}

void SelectionOptionInfo::setEvaluateForces(bool bEnabled)
{
    option().setSelectionFlag(efSelection_EvaluateForces, bEnabled);
}

void SelectionOptionInfo::setOnlyAtoms(bool bEnabled)
{
    option().setSelectionFlag(efSelection_OnlyAtoms, bEnabled);
}

void SelectionOptionInfo::setOnlyStatic(bool bEnabled)
{
    option().setSelectionFlag(efSelection_OnlyStatic, bEnabled);
}

void SelectionOptionInfo::setDynamicMask(bool bEnabled)
{
    option().setSelectionFlag(efSelection_DynamicMask, bEnabled);
}

void SelectionOptionInfo::setDynamicOnlyWhole(bool bEnabled)
{
    option().setSelectionFlag(efSelection_DynamicOnlyWhole, bEnabled);
}


/********************************************************************
 * SelectionOption
 */

AbstractOptionStoragePointer SelectionOption::createStorage() const
{
    return AbstractOptionStoragePointer(new SelectionOptionStorage(*this));
}


/********************************************************************
 * SelectionFileOptionStorage
 */

SelectionFileOptionStorage::SelectionFileOptionStorage(const SelectionFileOption &settings)
    : AbstractOptionStorage(settings, OptionFlags() | efOption_MultipleTimes
                            | efOption_DontCheckMinimumCount),
      info_(this), manager_(NULL), bValueParsed_(false)
{
}

void SelectionFileOptionStorage::clearSet()
{
    bValueParsed_ = false;
}

void SelectionFileOptionStorage::convertValue(const std::string &value)
{
    GMX_RELEASE_ASSERT(manager_ != NULL, "Manager is not set");

    if (bValueParsed_)
    {
        GMX_THROW(InvalidInputError("More than one file name provided"));
    }
    bValueParsed_ = true;
    // TODO: Should we throw an InvalidInputError if the file does not exist?
    manager_->parseRequestedFromFile(value);
}

void SelectionFileOptionStorage::processSet()
{
    if (!bValueParsed_)
    {
        GMX_THROW(InvalidInputError("No file name provided"));
    }
}


/********************************************************************
 * SelectionFileOptionInfo
 */

SelectionFileOptionInfo::SelectionFileOptionInfo(SelectionFileOptionStorage *option)
    : OptionInfo(option)
{
}

SelectionFileOptionStorage &SelectionFileOptionInfo::option()
{
    return static_cast<SelectionFileOptionStorage &>(OptionInfo::option());
}

const SelectionFileOptionStorage &SelectionFileOptionInfo::option() const
{
    return static_cast<const SelectionFileOptionStorage &>(OptionInfo::option());
}

void SelectionFileOptionInfo::setManager(SelectionOptionManager *manager)
{
    option().setManager(manager);
}


/********************************************************************
 * SelectionFileOption
 */

SelectionFileOption::SelectionFileOption(const char *name)
    : AbstractOption(name)
{
    setDescription("Provide selections from files");
}

AbstractOptionStoragePointer SelectionFileOption::createStorage() const
{
    return AbstractOptionStoragePointer(new SelectionFileOptionStorage(*this));
}

} // namespace gmx
