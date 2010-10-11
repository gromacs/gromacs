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
 * Implements gmx::SelectionOptionStorage.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_selection
 */
#include "selectionoptionstorage.h"

#include <cassert>

#include <string>
#include <vector>

#include "gromacs/errorreporting/abstracterrorreporter.h"
#include "gromacs/fatalerror/fatalerror.h"
#include "gromacs/options/globalproperties.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionoption.h"

#include "selectioncollection-impl.h"

namespace gmx
{

/********************************************************************
 * SelectionOptionStorage
 */

SelectionOptionStorage::SelectionOptionStorage()
{
    MyBase::setFlag(efNoDefaultValue);
    MyBase::setFlag(efConversionMayNotAddValues);
    MyBase::setFlag(efDontCheckMinimumCount);
}


int SelectionOptionStorage::init(const SelectionOption &settings,
                                 Options *options)
{
    _selectionFlags = settings._selectionFlags;
    options->globalProperties().request(eogpSelectionCollection);
    return MyBase::init(settings, options);
}


std::string SelectionOptionStorage::formatValue(int i) const
{
    Selection *sel = values().at(i);
    return (sel != NULL ? sel->selectionText() : "");
}


int SelectionOptionStorage::addSelections(
        const std::vector<Selection *> &selections,
        bool bFullValue)
{
    if (bFullValue && selections.size() < static_cast<size_t>(minValueCount()))
    {
        return eeInvalidInput;
    }
    std::vector<Selection *>::const_iterator i;
    for (i = selections.begin(); i != selections.end(); ++i)
    {
        // TODO: Having this check in the parser would make interactive input
        // behave better.
        if (_selectionFlags.test(efOnlyStatic) && (*i)->isDynamic())
        {
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
                                         AbstractErrorReporter * /*errors*/)
{
    SelectionCollection *sc =
        hostOptions().globalProperties().selectionCollection();
    assert(sc != NULL);

    std::vector<Selection *> selections;
    // TODO: Implement reading from a file.
    int rc = sc->parseFromString(value, &selections);
    if (rc == 0)
    {
        rc = addSelections(selections, false);
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

        sc->_impl->requestSelections(name(), description(),
                                     maxValueCount(), this);
        setFlag(efSet);
    }
    return MyBase::processAll(errors);
}

} // namespace gmx
