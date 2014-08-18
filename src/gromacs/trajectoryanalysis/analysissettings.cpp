/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements classes in analysissettings.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "analysissettings.h"

#include "gromacs/fileio/trxio.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"

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
