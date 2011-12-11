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
#include "gromacs/selection/selectionoption.h"

#include <string>
#include <vector>

#include "gromacs/fatalerror/exceptions.h"
#include "gromacs/fatalerror/gmxassert.h"
#include "gromacs/fatalerror/messagestringcollector.h"
#include "gromacs/options/globalproperties.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectioncollection.h"

#include "selectioncollection-impl.h"
#include "selectionoptionstorage.h"

namespace gmx
{

/********************************************************************
 * SelectionOptionStorage
 */

SelectionOptionStorage::SelectionOptionStorage(const SelectionOption &settings,
                                               Options *options)
    : MyBase(settings, options,
             OptionFlags() | efNoDefaultValue | efDontCheckMinimumCount),
      _selectionFlags(settings._selectionFlags), _adjuster(NULL)
{
    options->globalProperties().request(eogpSelectionCollection);
    if (settings._adjuster != NULL)
    {
        _adjuster = new SelectionOptionAdjuster(this);
        *settings._adjuster = _adjuster;
    }
}


SelectionOptionStorage::~SelectionOptionStorage()
{
    delete _adjuster;
}


std::string SelectionOptionStorage::formatValue(int i) const
{
    Selection *sel = values().at(i);
    return (sel != NULL ? sel->selectionText() : "");
}


void SelectionOptionStorage::addSelections(
        const std::vector<Selection *> &selections,
        bool bFullValue)
{
    if (bFullValue && selections.size() < static_cast<size_t>(minValueCount()))
    {
        GMX_THROW(InvalidInputError("Too few selections provided"));
    }
    if (bFullValue)
    {
        clearSet();
    }
    std::vector<Selection *>::const_iterator i;
    for (i = selections.begin(); i != selections.end(); ++i)
    {
        // TODO: Having this check in the parser would make interactive input
        // behave better.
        if (_selectionFlags.test(efOnlyStatic) && (*i)->isDynamic())
        {
            GMX_THROW(InvalidInputError("Dynamic selections not supported"));
        }
        (*i)->setFlags(_selectionFlags);
        addValue(*i);
    }
    if (bFullValue)
    {
        commitValues();
    }
}


void SelectionOptionStorage::convertValue(const std::string &value)
{
    SelectionCollection *sc =
        hostOptions().globalProperties().selectionCollection();
    GMX_RELEASE_ASSERT(sc != NULL, "Selection collection is not set");

    std::vector<Selection *> selections;
    // TODO: Implement reading from a file.
    sc->parseFromString(value, &selections);
    addSelections(selections, false);
}

void SelectionOptionStorage::processSetValues(ValueList *values)
{
    if (values->size() > 0 && values->size() < static_cast<size_t>(minValueCount()))
    {
        GMX_THROW(InvalidInputError("Too few (valid) values provided"));
    }
}

void SelectionOptionStorage::processAll()
{
    if ((hasFlag(efRequired) || hasFlag(efSet)) && valueCount() == 0)
    {
        SelectionCollection *sc =
            hostOptions().globalProperties().selectionCollection();
        GMX_RELEASE_ASSERT(sc != NULL, "Selection collection is not set");

        sc->_impl->requestSelections(name(), description(), this);
        setFlag(efSet);
    }
}

void SelectionOptionStorage::setAllowedValueCount(int count)
{
    MessageStringCollector errors;
    errors.startContext("In option '" + name() + "'");
    if (count >= 0)
    {
        // Should not throw because efDontCheckMinimumCount is set
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
    catch (UserInputError &ex)
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
    _selectionFlags.set(flag, bSet);
    ValueList::const_iterator i;
    for (i = values().begin(); i != values().end(); ++i)
    {
        if (_selectionFlags.test(efOnlyStatic) && (*i)->isDynamic())
        {
            MessageStringCollector errors;
            errors.startContext("In option '" + name() + "'");
            errors.append("Dynamic selections not supported");
            errors.finishContext();
            GMX_THROW(InvalidInputError(errors.toString()));
        }
        (*i)->setFlags(_selectionFlags);
    }
}


/********************************************************************
 * SelectionOptionAdjuster
 */

SelectionOptionAdjuster::SelectionOptionAdjuster(SelectionOptionStorage *storage)
    : _storage(*storage)
{
}

void SelectionOptionAdjuster::setValueCount(int count)
{
    return storage().setAllowedValueCount(count);
}

void SelectionOptionAdjuster::setEvaluateVelocities(bool bEnabled)
{
    return storage().setSelectionFlag(efEvaluateVelocities, bEnabled);
}

void SelectionOptionAdjuster::setEvaluateForces(bool bEnabled)
{
    return storage().setSelectionFlag(efEvaluateForces, bEnabled);
}

void SelectionOptionAdjuster::setOnlyAtoms(bool bEnabled)
{
    return storage().setSelectionFlag(efOnlyAtoms, bEnabled);
}

void SelectionOptionAdjuster::setOnlyStatic(bool bEnabled)
{
    return storage().setSelectionFlag(efOnlyStatic, bEnabled);
}

void SelectionOptionAdjuster::setDynamicMask(bool bEnabled)
{
    return storage().setSelectionFlag(efDynamicMask, bEnabled);
}

void SelectionOptionAdjuster::setDynamicOnlyWhole(bool bEnabled)
{
    return storage().setSelectionFlag(efDynamicOnlyWhole, bEnabled);
}


/********************************************************************
 * SelectionOption
 */

AbstractOptionStorage *SelectionOption::createDefaultStorage(Options *options) const
{
    return new SelectionOptionStorage(*this, options);
}

} // namespace gmx
