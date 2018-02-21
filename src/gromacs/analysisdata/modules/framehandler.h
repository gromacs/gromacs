/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Declares gmx::TrajectoryDataWriteModule for writing trajectory data (into a file).
 *
 * \inpublicapi
 * \ingroup module_trajectorydata
 * \author 
 */
#ifndef GMX_TRAJECTORYDATA_FRAMEHANDLER_H
#define GMX_TRAJECTORYDATA_FRAMEHANDLER_H

#include <memory>
#include <string>

#include "settings.h"

#include "gromacs/options/timeunitmanager.h"
#include "gromacs/utility/classhelpers.h"

struct t_atoms;
struct t_topology;
struct gmx_mtop_t;
struct t_trxstatus;
struct t_trxframe;

namespace gmx
{

class IOptionsContainer;
class SelectionCollection;
class Selection;

/*! \brief
 * Abstract data module for writing data into a file.
 *
 * Implements features common to all plotting modules.  Subclasses implement
 * features specific to certain applications (TrajectoryDataWriteModule implements
 * straightforward plotting).
 *
 * By default, the data is written into an xvgr file, according to the
 * options read from the TrajectoryDataWriteSettings object given to the
 * constructor.
 * For non-xvgr data, it's possible to skip all headers by calling
 * setPlainOutput().
 *
 * A single output line corresponds to a single frame.  In most cases with
 * multipoint data, setPlainOutput() should be called since the output does not
 * make sense as an xvgr file, but this is not enforced.
 *
 * Multipoint data and multiple data sets are both supported, in which case all
 * the points are written to the output, in the order in which they are added
 * to the data.
 *
 * \ingroup module_analysisdata
 */

class Framehandler 
{
    public:
        virtual ~Framehandler();

        explicit Framehandler(TrajectoryDataWriteSettings *settings);
        Framehandler();


        /*! \brief
         * Set common settings for the plotting.
         */
        void setSettings(TrajectoryDataWriteSettings *settings);

        void modifyFrame(t_trxframe *newFrame, const t_trxframe *oldFrame);


    private:
        TrajectoryDataWriteSettings *settings_;
};


} // namespace gmx

#endif
