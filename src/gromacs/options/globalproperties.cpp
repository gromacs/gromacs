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

#include <cstddef>

#include <smalloc.h>
#include <statutil.h>

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"

namespace gmx
{

static const char *const plotFormats[] = {
    "none", "xmgrace", "xmgr", NULL
};


OptionsGlobalProperties::OptionsGlobalProperties()
    : _usedProperties(0), _plotFormat(1),
      _oenv(NULL)
{
    // TODO: If/when this is refactored, exception safety should be considered
    snew(_oenv, 1);
    output_env_init_default(_oenv);
}


OptionsGlobalProperties::~OptionsGlobalProperties()
{
    if (_oenv != NULL)
    {
        output_env_done(_oenv);
    }
}


void OptionsGlobalProperties::addDefaultOptions(Options *options)
{
    if (isPropertyUsed(eogpPlotFormat))
    {
        options->addOption(StringOption("xvg").enumValue(plotFormats)
                               .defaultValue("xmgrace")
                               .storeEnumIndex(&_plotFormat)
                               .description("Plot formatting"));
    }
}


void OptionsGlobalProperties::finish()
{
    if (isPropertyUsed(eogpPlotFormat))
    {
        if (_plotFormat == 0)
        {
            _oenv->xvg_format = exvgNONE;
        }
        else
        {
            _oenv->xvg_format = static_cast<xvg_format_t>(_plotFormat);
        }
    }
}

} // namespace gmx
