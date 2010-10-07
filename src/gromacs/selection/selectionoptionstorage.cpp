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

#include "gromacs/options/globalproperties.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionoption.h"

namespace gmx
{

/********************************************************************
 * SelectionOptionStorage
 */

SelectionOptionStorage::SelectionOptionStorage()
{
    MyBase::setFlag(efConversionMayNotAddValues);
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
    return values().at(i)->selectionText();
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
        std::vector<Selection *>::const_iterator i;
        for (i = selections.begin(); i != selections.end(); ++i)
        {
            (*i)->setFlags(_selectionFlags);
            addValue(*i);
            // FIXME: Check the return value
        }
    }

    return rc;
}

int SelectionOptionStorage::processSet(int nvalues,
                                       AbstractErrorReporter *errors)
{
    return MyBase::processSet(nvalues, errors);
}

int SelectionOptionStorage::processAll(AbstractErrorReporter *errors)
{
    if ((hasFlag(efRequired) || hasFlag(efSet)) && valueCount() == 0)
    {
        SelectionCollection *sc =
            hostOptions().globalProperties().selectionCollection();
        assert(sc != NULL);

        // TODO: The actual interactive input should really happen later,
        // because otherwise some completely unrelated errors in user input may
        // be detected only after the user has been prompted for selections.
        // TODO: Check whether the input is a pipe.
        bool bInteractive = true;
        if (bInteractive)
        {
            std::fprintf(stderr, "\nSpecify ");
            if (maxValueCount() < 0)
            {
                std::fprintf(stderr, "any number of selections");
            }
            else if (maxValueCount() == 1)
            {
                std::fprintf(stderr, "a selection");
            }
            else
            {
                std::fprintf(stderr, "%d selections", maxValueCount());
            }
            std::fprintf(stderr, " for option '%s' (%s):\n",
                         name().c_str(), description().c_str());
            std::fprintf(stderr, "(one selection per line, 'help' for help%s)\n",
                         maxValueCount() < 0 ? ", Ctrl-D to end" : "");
        }
        std::vector<Selection *> selections;
        int rc = sc->parseFromStdin(maxValueCount(), bInteractive, &selections);
        if (rc != 0)
        {
            return rc;
        }
        std::vector<Selection *>::const_iterator i;
        for (i = selections.begin(); i != selections.end(); ++i)
        {
            addValue(*i);
            // FIXME: Check return value.
        }
        setFlag(efSet);
    }
    return MyBase::processAll(errors);
}

} // namespace gmx
