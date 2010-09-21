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
 * Implements gmx::OptionsGlobalProperties.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_options
 */
#include "gromacs/options/globalproperties.h"

#include <cassert>
#include <cstddef>

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"

namespace gmx
{

static const char *const timeUnits[] = {
    "fs", "ps", "ns", "us", "ms",  "s", NULL
};
static const double timeScaleFactors[] = {
    1e-3,    1,  1e3,  1e6,  1e9, 1e12
};

OptionsGlobalProperties::OptionsGlobalProperties()
    : _timeUnit(1), _usedProperties(0)
{
}

double OptionsGlobalProperties::timeScaleFactor() const
{
    assert(_timeUnit >= 0
           && (size_t)_timeUnit < sizeof(timeScaleFactors)/sizeof(timeScaleFactors[0]));
    return timeScaleFactors[_timeUnit];
}

void OptionsGlobalProperties::addDefaultOptions(Options *options)
{
    if (isPropertyUsed(eogpTimeScaleFactor))
    {
        options->addOption(StringOption("tu").enumValue(timeUnits)
                               .defaultValue("ps")
                               .storeEnumIndex(&_timeUnit)
                               .description("Unit for time values"));
    }
}

} // namespace gmx
