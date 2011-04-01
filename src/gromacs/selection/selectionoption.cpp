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

#include <cassert>

#include <string>
#include <vector>

#include "gromacs/errorreporting/abstracterrorreporter.h"
#include "gromacs/errorreporting/errorcontext.h"
#include "gromacs/fatalerror/fatalerror.h"
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

SelectionOptionStorage::SelectionOptionStorage()
    : _adjuster(NULL)
{
    MyBase::setFlag(efNoDefaultValue);
    MyBase::setFlag(efConversionMayNotAddValues);
    MyBase::setFlag(efDontCheckMinimumCount);
}


SelectionOptionStorage::~SelectionOptionStorage()
{
    delete _adjuster;
}


int SelectionOptionStorage::init(const SelectionOption &settings,
                                 Options *options)
{
    _selectionFlags = settings._selectionFlags;
    int rc = MyBase::init(settings, options);
    if (rc == 0)
    {
        options->globalProperties().request(eogpSelectionCollection);
        if (settings._adjuster)
        {
            _adjuster = new SelectionOptionAdjuster(this);
            *settings._adjuster = _adjuster;
        }
    }
    return rc;
}


std::string SelectionOptionStorage::formatValue(int i) const
{
    Selection *sel = values().at(i);
    return (sel != NULL ? sel->selectionText() : "");
}


int SelectionOptionStorage::addSelections(
        const std::vector<Selection *> &selections,
        bool bFullValue, AbstractErrorReporter *errors)
{
    if (bFullValue && selections.size() < static_cast<size_t>(minValueCount()))
    {
        errors->error("Too few selections provided");
        return eeInvalidInput;
    }
    std::vector<Selection *>::const_iterator i;
    for (i = selections.begin(); i != selections.end(); ++i)
    {
        // TODO: Having this check in the parser would make interactive input
        // behave better.
        if (_selectionFlags.test(efOnlyStatic) && (*i)->isDynamic())
        {
            errors->error("Dynamic selections not supported");
            return eeInvalidInput;
        }
        (*i)->setFlags(_selectionFlags);
        int rc = addValue(*i);
        if (rc != 0)
        {
            return rc;
        }
    }
    if (bFullValue)
    {
        processValues(selections.size(), true);
    }
    return 0;
}


int SelectionOptionStorage::convertValue(const std::string &value,
                                         AbstractErrorReporter *errors)
{
    SelectionCollection *sc =
        hostOptions().globalProperties().selectionCollection();
    assert(sc != NULL);

    std::vector<Selection *> selections;
    // TODO: Implement reading from a file.
    int rc = sc->parseFromString(value, errors, &selections);
    if (rc == 0)
    {
        rc = addSelections(selections, false, errors);
    }
    return rc;
}

int SelectionOptionStorage::processSet(int nvalues,
                                       AbstractErrorReporter *errors)
{
    if (nvalues > 0 && nvalues < minValueCount())
    {
        // TODO: Remove the invalid values
        errors->error("Too few (valid) values provided");
        return eeInvalidInput;
    }
    return MyBase::processSet(nvalues, errors);
}

int SelectionOptionStorage::processAll(AbstractErrorReporter *errors)
{
    if ((hasFlag(efRequired) || hasFlag(efSet)) && valueCount() == 0)
    {
        SelectionCollection *sc =
            hostOptions().globalProperties().selectionCollection();
        assert(sc != NULL);

        sc->_impl->requestSelections(name(), description(), this);
        setFlag(efSet);
    }
    return MyBase::processAll(errors);
}

int SelectionOptionStorage::setAllowedValueCount(int count,
                                                 AbstractErrorReporter *errors)
{
    ErrorContext context(errors, "In option '" + name() + "'");
    int rc = 0;
    if (count > 0)
    {
        rc = setMinValueCount(count, errors);
        if (rc == 0 && valueCount() > 0 && valueCount() < count)
        {
            errors->error("Too few (valid) values provided");
            rc = eeInvalidInput;
        }
    }
    int rc1 = setMaxValueCount(count, errors);
    return rc != 0 ? rc : rc1;
}

int SelectionOptionStorage::setSelectionFlag(SelectionFlag flag, bool bSet,
                                             AbstractErrorReporter *errors)
{
    ErrorContext context(errors, "In option '" + name() + "'");
    _selectionFlags.set(flag, bSet);
    ValueList::const_iterator i;
    for (i = values().begin(); i != values().end(); ++i)
    {
        if (_selectionFlags.test(efOnlyStatic) && (*i)->isDynamic())
        {
            errors->error("Dynamic selections not supported");
            return eeInvalidInput;
        }
        (*i)->setFlags(_selectionFlags);
    }
    return 0;
}


/********************************************************************
 * SelectionOptionAdjuster
 */

SelectionOptionAdjuster::SelectionOptionAdjuster(SelectionOptionStorage *storage)
    : _storage(*storage), _errors(NULL)
{
}

AbstractErrorReporter *
SelectionOptionAdjuster::setErrorReporter(AbstractErrorReporter *errors)
{
    AbstractErrorReporter *old = _errors;
    _errors = errors;
    return old;
}

int SelectionOptionAdjuster::setValueCount(int count)
{
    return storage().setAllowedValueCount(count, errors());
}

int SelectionOptionAdjuster::setEvaluateVelocities(bool bEnabled)
{
    return storage().setSelectionFlag(efEvaluateVelocities, bEnabled, errors());
}

int SelectionOptionAdjuster::setEvaluateForces(bool bEnabled)
{
    return storage().setSelectionFlag(efEvaluateForces, bEnabled, errors());
}

int SelectionOptionAdjuster::setOnlyAtoms(bool bEnabled)
{
    return storage().setSelectionFlag(efOnlyAtoms, bEnabled, errors());
}

int SelectionOptionAdjuster::setOnlyStatic(bool bEnabled)
{
    return storage().setSelectionFlag(efOnlyStatic, bEnabled, errors());
}

int SelectionOptionAdjuster::setDynamicMask(bool bEnabled)
{
    return storage().setSelectionFlag(efDynamicMask, bEnabled, errors());
}

int SelectionOptionAdjuster::setDynamicOnlyWhole(bool bEnabled)
{
    return storage().setSelectionFlag(efDynamicOnlyWhole, bEnabled, errors());
}


/********************************************************************
 * SelectionOption
 */

int SelectionOption::createDefaultStorage(Options *options,
                                          AbstractOptionStorage **storage) const
{
    return createOptionStorage<SelectionOption, SelectionOptionStorage>(this, options, storage);
}

} // namespace gmx
