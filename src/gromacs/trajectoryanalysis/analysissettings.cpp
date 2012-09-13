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
    : impl_(new Impl)
{
    impl_->frflags |= TRX_NEED_X;
}


TrajectoryAnalysisSettings::~TrajectoryAnalysisSettings()
{
}


const TimeUnitManager &
TrajectoryAnalysisSettings::timeUnitManager() const
{
    return impl_->timeUnitManager;
}


const AnalysisDataPlotSettings &
TrajectoryAnalysisSettings::plotSettings() const
{
    return impl_->plotSettings;
}


unsigned long
TrajectoryAnalysisSettings::flags() const
{
    return impl_->flags;
}


bool
TrajectoryAnalysisSettings::hasFlag(unsigned long flag) const
{
    return impl_->flags & flag;
}


bool
TrajectoryAnalysisSettings::hasPBC() const
{
    return impl_->bPBC;
}


bool
TrajectoryAnalysisSettings::hasRmPBC() const
{
    return impl_->bRmPBC;
}


int
TrajectoryAnalysisSettings::frflags() const
{
    return impl_->frflags;
}


void
TrajectoryAnalysisSettings::setFlags(unsigned long flags)
{
    impl_->flags = flags;
}


void
TrajectoryAnalysisSettings::setFlag(unsigned long flag, bool bSet)
{
    if (bSet)
    {
        impl_->flags |= flag;
    }
    else
    {
        impl_->flags &= ~flag;
    }
}


void
TrajectoryAnalysisSettings::setPBC(bool bPBC)
{
    impl_->bPBC = bPBC;
}


void
TrajectoryAnalysisSettings::setRmPBC(bool bRmPBC)
{
    impl_->bRmPBC = bRmPBC;
}


void
TrajectoryAnalysisSettings::setFrameFlags(int frflags)
{
    impl_->frflags = frflags;
}


/********************************************************************
 * TopologyInformation
 */

TopologyInformation::TopologyInformation()
    : top_(NULL), bTop_(false), xtop_(NULL), ePBC_(-1)
{
    clear_mat(boxtop_);
}


TopologyInformation::~TopologyInformation()
{
    if (top_)
    {
        free_t_atoms(&top_->atoms, TRUE);
        done_top(top_);
        sfree(top_);
    }
    sfree(xtop_);
}


void
TopologyInformation::getTopologyConf(rvec **x, matrix box) const
{
    if (box)
    {
        copy_mat(const_cast<rvec *>(boxtop_), box);
    }
    if (x)
    {
        if (!xtop_)
        {
            *x = NULL;
            GMX_THROW(APIError("Topology coordinates requested without setting efUseTopX"));
        }
        *x = xtop_;
    }
}

} // namespace gmx
