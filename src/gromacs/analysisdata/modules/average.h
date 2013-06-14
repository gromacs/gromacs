/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * Declares gmx::AnalysisDataAverageModule.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_MODULES_AVERAGE_H
#define GMX_ANALYSISDATA_MODULES_AVERAGE_H

#include <vector>

#include "../abstractdata.h"
#include "../arraydata.h"
#include "../datamodule.h"
#include "../../utility/common.h"

namespace gmx
{

/*! \brief
 * Data module for independently averaging each column in input data.
 *
 * Computes the average and standard deviation independently for each column in
 * the input data.  Multipoint data and missing data points are both supported.
 * The average is always calculated over all frames and data points for a
 * column.
 * Multiple input data sets are currently not supported.
 *
 * Output data contains a frame for each column in the input data.
 * The first column of each output frame is the average of the corresponding
 * input column.  The second output column is the standard deviation of the
 * corresponding input column.
 * average() and stddev() methods are also provided for convenient access to
 * these properties.
 *
 * The output data becomes available only after the input data has been
 * finished.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisDataAverageModule : public AbstractAnalysisArrayData,
                                  public AnalysisDataModuleInterface
{
    public:
        AnalysisDataAverageModule();
        virtual ~AnalysisDataAverageModule();

        using AbstractAnalysisArrayData::setXAxis;

        virtual int flags() const;

        virtual void dataStarted(AbstractAnalysisData *data);
        virtual void frameStarted(const AnalysisDataFrameHeader &header);
        virtual void pointsAdded(const AnalysisDataPointSetRef &points);
        virtual void frameFinished(const AnalysisDataFrameHeader &header);
        virtual void dataFinished();

        //! Convenience access to the average of a data column.
        real average(int index) const;
        //! Convenience access to the standard deviation of a data column.
        real stddev(int index) const;

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

//! Smart pointer to manage an AnalysisDataAverageModule object.
typedef boost::shared_ptr<AnalysisDataAverageModule>
    AnalysisDataAverageModulePointer;

/*! \brief
 * Data module for averaging of columns for each frame.
 *
 * Output data has the same number of frames as the input data.
 * The number of columns in the output data is the same as the number of data
 * sets in the input data.
 * Each frame in the output contains the average of the column values for each
 * data set in the corresponding frame of the input data.
 *
 * Multipoint data and missing data points are both supported.  The average
 * is always calculated over all data points present in a column for a data
 * set.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisDataFrameAverageModule : public AbstractAnalysisData,
                                       public AnalysisDataModuleInterface
{
    public:
        AnalysisDataFrameAverageModule();
        virtual ~AnalysisDataFrameAverageModule();

        virtual int flags() const;

        virtual void dataStarted(AbstractAnalysisData *data);
        virtual void frameStarted(const AnalysisDataFrameHeader &header);
        virtual void pointsAdded(const AnalysisDataPointSetRef &points);
        virtual void frameFinished(const AnalysisDataFrameHeader &header);
        virtual void dataFinished();

    private:
        virtual AnalysisDataFrameRef tryGetDataFrameInternal(int index) const;
        virtual bool requestStorageInternal(int nframes);

        class Impl;

        PrivateImplPointer<Impl> impl_;
};

//! Smart pointer to manage an AnalysisDataFrameAverageModule object.
typedef boost::shared_ptr<AnalysisDataFrameAverageModule>
    AnalysisDataFrameAverageModulePointer;

} // namespace gmx

#endif
