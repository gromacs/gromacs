/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
 * \ingroup module_analysisdata
 * \author
 */
#ifndef GMX_ANALYSISDATA_MODULES_WRITE_H
#define GMX_ANALYSISDATA_MODULES_WRITE_H

#include <memory>
#include <string>

#include "gromacs/analysisdata/abstractdata.h"
#include "gromacs/analysisdata/arraydata.h"
#include "gromacs/analysisdata/datamodule.h"
#include "gromacs/analysisdata/modules/settings.h"
#include "gromacs/options/timeunitmanager.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class AnalysisDataValue;
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
class AbstractWriteModule : public AbstractAnalysisArrayData,
                            public AnalysisDataModuleParallel
{
    public:
        virtual ~AbstractWriteModule();

        /*! \brief
         * Set common settings for the plotting.
         */
        void setSettings(const TrajectoryDataWriteSettings &settings);

        virtual int flags() const;

        virtual bool parallelDataStarted(AbstractAnalysisData              *data,
                                         const AnalysisDataParallelOptions &options);
        virtual void frameStarted(const AnalysisDataFrameHeader & /*header*/);
        virtual void pointsAdded(const AnalysisDataPointSetRef &points);
        virtual void frameFinished(const AnalysisDataFrameHeader &header);
        virtual void frameFinishedSerial(int frameIndex);
        virtual void dataFinished();

        void setExternal(const Selection *sel, std::string name, const gmx_mtop_t *mtop, const t_topology *top);

    protected:
        /*! \cond libapi */
        AbstractWriteModule();
        //! Creates AbstractWriteModule and assign common settings.
        explicit AbstractWriteModule(const TrajectoryDataWriteSettings &settings);

        //! \endcond

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};


/*! \brief
 * TrjWriteting module for straightforward plotting of data.
 *
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class TrajectoryDataWriteModule : public AbstractWriteModule
{
    public:
        TrajectoryDataWriteModule();
        //! Creates AnalysisDataWriteModule and assign common settings.
        explicit TrajectoryDataWriteModule(const TrajectoryDataWriteSettings &settings);

        //       virtual void dataStarted(AbstractAnalysisData *data);
        //       virtual void frameStarted(const AnalysisDataFrameHeader &header);
        //       virtual void pointsAdded(const AnalysisDataPointSetRef &points);
        //       virtual void frameFinished(const AnalysisDataFrameHeader &header);
        //       virtual void dataFinished();

        // Copy and assign disallowed by base.
};


//! Smart pointer to manage an AnalysisDataWriteModule object.
typedef std::shared_ptr<TrajectoryDataWriteModule>
    TrajectoryDataWriteModulePointer;

} // namespace gmx

#endif
