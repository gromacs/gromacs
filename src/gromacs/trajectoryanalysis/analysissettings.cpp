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
 * Implements classes in analysissettings.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_trajectoryanalysis
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "smalloc.h"
#include "statutil.h"
#include "vec.h"

#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/exceptions.h"

#include "analysissettings-impl.h"

namespace gmx
{


/********************************************************************
 * TrajectoryAnalysisSettings
 */

TrajectoryAnalysisSettings::TrajectoryAnalysisSettings()
    : _impl(new Impl)
{
    _impl->frflags |= TRX_NEED_X;
}


TrajectoryAnalysisSettings::~TrajectoryAnalysisSettings()
{
}


const TimeUnitManager &
TrajectoryAnalysisSettings::timeUnitManager() const
{
    return _impl->timeUnitManager;
}


const AnalysisDataPlotSettings &
TrajectoryAnalysisSettings::plotSettings() const
{
    return _impl->plotSettings;
}


unsigned long
TrajectoryAnalysisSettings::flags() const
{
    return _impl->flags;
}


bool
TrajectoryAnalysisSettings::hasFlag(unsigned long flag) const
{
    return _impl->flags & flag;
}


bool
TrajectoryAnalysisSettings::hasPBC() const
{
    return _impl->bPBC;
}


bool
TrajectoryAnalysisSettings::hasRmPBC() const
{
    return _impl->bRmPBC;
}


int
TrajectoryAnalysisSettings::frflags() const
{
    return _impl->frflags;
}


void
TrajectoryAnalysisSettings::setFlags(unsigned long flags)
{
    _impl->flags = flags;
}


void
TrajectoryAnalysisSettings::setFlag(unsigned long flag, bool bSet)
{
    if (bSet)
    {
        _impl->flags |= flag;
    }
    else
    {
        _impl->flags &= ~flag;
    }
}


void
TrajectoryAnalysisSettings::setPBC(bool bPBC)
{
    _impl->bPBC = bPBC;
}


void
TrajectoryAnalysisSettings::setRmPBC(bool bRmPBC)
{
    _impl->bRmPBC = bRmPBC;
}


void
TrajectoryAnalysisSettings::setFrameFlags(int frflags)
{
    _impl->frflags = frflags;
}


/********************************************************************
 * TopologyInformation
 */

TopologyInformation::TopologyInformation()
    : _top(NULL), _bTop(false), _xtop(NULL), _ePBC(-1)
{
    clear_mat(_boxtop);
}


TopologyInformation::~TopologyInformation()
{
    if (_top)
    {
        done_top(_top);
        sfree(_top);
    }
    sfree(_xtop);
}


void
TopologyInformation::getTopologyConf(rvec **x, matrix box) const
{
    if (box)
    {
        copy_mat(const_cast<rvec *>(_boxtop), box);
    }
    if (x)
    {
        if (!_xtop)
        {
            *x = NULL;
            GMX_THROW(APIError("Topology coordinates requested without setting efUseTopX"));
        }
        *x = _xtop;
    }
}

} // namespace gmx
