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
 * Implements createTrajectoryAnalysisModule().
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_trajectoryanalysis
 */
#include "gromacs/trajectoryanalysis/modules.h"

#include "string2.h"

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/format.h"

#include "modules/angle.h"
#include "modules/distance.h"
#include "modules/select.h"

namespace
{

using namespace gmx::analysismodules;

struct module_map_t
{
    const char                            *name;
    gmx::TrajectoryAnalysisModulePointer (*creator)(void);
};

const module_map_t modules[] =
{
    {gmx::analysismodules::angle,    Angle::create},
    {gmx::analysismodules::distance, Distance::create},
    {gmx::analysismodules::select,   Select::create},
    {NULL,                           NULL},
};

} // namespace

namespace gmx
{

TrajectoryAnalysisModulePointer
createTrajectoryAnalysisModule(const char *name)
{
    size_t len = strlen(name);
    int match_i = -1;

    for (int i = 0; modules[i].name != NULL; ++i)
    {
        if (gmx_strncasecmp(name, modules[i].name, len) == 0)
        {
            if (strlen(modules[i].name) == len)
            {
                match_i = i;
                break;
            }
            else if (match_i == -1)
            {
                match_i = i;
            }
            else
            {
                GMX_THROW(InvalidInputError(
                            gmx::formatString("Requested analysis module '%s' is ambiguous", name)));
            }
        }
    }
    if (match_i != -1)
    {
        return modules[match_i].creator();
    }
    GMX_THROW(InvalidInputError(
                gmx::formatString("Unknown analysis module: %s", name)));
}

} // namespace gmx
