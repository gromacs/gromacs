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
 * Implements classes in selectionoption.h and selectionoptionstorage.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "gromacs/selection/selectionoption.h"

#include <cstddef>

#include <string>
#include <vector>

#include "gromacs/options/abstractoption.h"
#include "gromacs/options/abstractoptionstorage.h"
#include "gromacs/options/optionflags.h"
#include "gromacs/options/optionmanagercontainer.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionenums.h"
#include "gromacs/selection/selectionfileoption.h"
#include "gromacs/selection/selectionoptionmanager.h"
#include "gromacs/utility/any.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/flags.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/message_string_collector.h"

#include "selectionfileoptionstorage.h"
#include "selectionoptionstorage.h"

namespace gmx
{

/********************************************************************
 * SelectionOptionStorage
 */

SelectionOptionStorage::SelectionOptionStorage(const SelectionOption&  settings,
                                               SelectionOptionManager* manager) :
    MyBase(settings, OptionFlags() | efOption_NoDefaultValue | efOption_DontCheckMinimumCount),
    info_(this),
    manager_(*manager),
    defaultText_(settings.defaultText_),
    selectionFlags_(settings.selectionFlags_)
{
    GMX_RELEASE_ASSERT(manager != nullptr,
                       "SelectionOptionManager must be added before SelectionOption");
    GMX_RELEASE_ASSERT(!hasFlag(efOption_MultipleTimes),
                       "allowMultiple() is not supported for selection options");
    manager_.registerOption(this);
}


std::string SelectionOptionStorage::formatSingleValue(const Selection& value) const
{
    return value.selectionText();
}


std::vector<Any> SelectionOptionStorage::normalizeValues(const std::vector<Any>& /*values*/) const
{
    GMX_THROW(NotImplementedError("Selection options not supported in this context"));
}


void SelectionOptionStorage::addSelections(const SelectionList& selections, bool bFullValue)
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
        // Create a copy to allow modifications.
        Selection sel(*i);
        sel.data().setFlags(selectionFlags_);
        addValue(sel);
    }
    if (bFullValue)
    {
        commitValues();
        markAsSet();
    }
}


void SelectionOptionStorage::convertValue(const Any& value)
{
    manager_.convertOptionValue(this, value.cast<std::string>(), false);
}

void SelectionOptionStorage::processSetValues(ValueList* values)
{
    if (values->empty())
    {
        manager_.requestOptionDelayedParsing(this);
    }
    else if (values->size() < static_cast<size_t>(minValueCount()))
    {
        GMX_THROW(InvalidInputError("Too few (valid) values provided"));
    }
}

void SelectionOptionStorage::processAll()
{
    if (!isSet() && !defaultText_.empty())
    {
        manager_.convertOptionValue(this, defaultText_, true);
    }
    if (isRequired() && !isSet())
    {
        manager_.requestOptionDelayedParsing(this);
        markAsSet();
    }
}

void SelectionOptionStorage::setAllowedValueCount(int count)
{
    // TODO: It should be possible to have strong exception safety here.
    // TODO: Use ExceptionInitializer here.
    MessageStringCollector errors;
    errors.startContext("In option '" + name() + "'");
    if (count >= 0)
    {
        // Should not throw because efOption_DontCheckMinimumCount is set.
        setMinValueCount(count);
        errors.appendIf((valueCount() > 0 && valueCount() < count),
                        "Too few (valid) values provided");
    }
    try
    {
        setMaxValueCount(count);
    }
    catch (const UserInputError& ex)
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
    for (const Selection& value : values())
    {
        if (flag == efSelection_OnlyStatic && bSet && value.isDynamic())
        {
            MessageStringCollector errors;
            errors.startContext("In option '" + name() + "'");
            errors.append("Dynamic selections not supported");
            errors.finishContext();
            GMX_THROW(InvalidInputError(errors.toString()));
        }
    }
    selectionFlags_.set(flag, bSet);
    for (Selection& value : values())
    {
        value.data().setFlags(selectionFlags_);
    }
}


/********************************************************************
 * SelectionOptionInfo
 */

SelectionOptionInfo::SelectionOptionInfo(SelectionOptionStorage* option) : OptionInfo(option) {}

SelectionOptionStorage& SelectionOptionInfo::option()
{
    return static_cast<SelectionOptionStorage&>(OptionInfo::option());
}

const SelectionOptionStorage& SelectionOptionInfo::option() const
{
    return static_cast<const SelectionOptionStorage&>(OptionInfo::option());
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


/********************************************************************
 * SelectionOption
 */

AbstractOptionStorage* SelectionOption::createStorage(const OptionManagerContainer& managers) const
{
    return new SelectionOptionStorage(*this, managers.get<SelectionOptionManager>());
}


/********************************************************************
 * SelectionFileOptionStorage
 */

SelectionFileOptionStorage::SelectionFileOptionStorage(const SelectionFileOption& settings,
                                                       SelectionOptionManager*    manager) :
    AbstractOptionStorage(settings, OptionFlags() | efOption_MultipleTimes | efOption_DontCheckMinimumCount),
    info_(this),
    manager_(*manager),
    bValueParsed_(false)
{
    GMX_RELEASE_ASSERT(manager != nullptr,
                       "SelectionOptionManager must be added before SelectionFileOption");
}

void SelectionFileOptionStorage::clearSet()
{
    bValueParsed_ = false;
}

void SelectionFileOptionStorage::convertValue(const Any& value)
{
    if (bValueParsed_)
    {
        GMX_THROW(InvalidInputError("More than one file name provided"));
    }
    bValueParsed_ = true;
    // TODO: Should we throw an InvalidInputError if the file does not exist?
    manager_.parseRequestedFromFile(value.cast<std::string>());
}

void SelectionFileOptionStorage::processSet()
{
    if (!bValueParsed_)
    {
        GMX_THROW(InvalidInputError("No file name provided"));
    }
}

std::vector<Any> SelectionFileOptionStorage::defaultValues() const
{
    return {};
}

std::vector<Any> SelectionFileOptionStorage::normalizeValues(const std::vector<Any>& values) const
{
    return values;
}


/********************************************************************
 * SelectionFileOptionInfo
 */

SelectionFileOptionInfo::SelectionFileOptionInfo(SelectionFileOptionStorage* option) :
    OptionInfo(option)
{
}


/********************************************************************
 * SelectionFileOption
 */

SelectionFileOption::SelectionFileOption(const char* name) : AbstractOption(name)
{
    setDescription("Provide selections from files");
}

AbstractOptionStorage* SelectionFileOption::createStorage(const OptionManagerContainer& managers) const
{
    return new SelectionFileOptionStorage(*this, managers.get<SelectionOptionManager>());
}

} // namespace gmx
