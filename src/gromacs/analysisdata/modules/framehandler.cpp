/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
 * Implements classes in framehandler.h.
 *
 * \ingroup module_trajectorydata
 * \author 
 */
#include "gmxpre.h"

#include "framehandler.h"
#include "settings.h"

#include <cstdio>
#include <cstring>

#include <string>
#include <vector>

#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/timeunitmanager.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"

namespace gmx
{

Framehandler::Framehandler(){}
Framehandler::Framehandler(TrajectoryDataWriteSettings *settings)
: settings_(settings)
{}

Framehandler::~Framehandler()
{}


void
Framehandler::setSettings(TrajectoryDataWriteSettings *settings)
{
    settings_ = settings;
}

void
Framehandler::modifyFrame(t_trxframe *newFrame, const t_trxframe *oldFrame)
{
    const Selection *sel = settings_->getInputSel();
    int natoms = sel->atomCount();

    *newFrame = *oldFrame;

    newFrame->time   = settings_->getTime();
    newFrame->bV     = (oldFrame->bV && settings_->getbVel());
    newFrame->bF     = (oldFrame->bF && settings_->getbForce());
    newFrame->natoms = natoms;
    newFrame->bPrec  = (oldFrame->bPrec && settings_->getbPrec());
    newFrame->prec   = settings_->getPrecision();
    newFrame->atoms  = settings_->getAtoms();
    newFrame->bAtoms = settings_->getbAtoms();

    rvec *xmem = nullptr;
    rvec *vmem = nullptr;
    rvec *fmem = nullptr;
    snew(xmem,natoms);
    if (newFrame->bV)
    {
        snew(vmem,natoms);
    }
    if (newFrame->bF)
    {
        snew(fmem,natoms);
    }
    newFrame->x = xmem;
    newFrame->v = vmem;
    newFrame->f = fmem;


    for (int i = 0; i < natoms; i++)
    {
        int pos = sel->position(i).refId();
        copy_rvec(oldFrame->x[pos], newFrame->x[i]);
        if (newFrame->bV)
        {
            copy_rvec(oldFrame->v[pos], newFrame->v[i]);
        }
        if (newFrame->bF)
        {
            copy_rvec(oldFrame->f[pos], newFrame->f[i]);
        }
    }
}
    



} // namespace gmx
