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

#include "gromacs/analysisdata/abstractdata.h"
#include "gromacs/analysisdata/arraydata.h"
#include "gromacs/analysisdata/datamodule.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

/*! \brief
 * Data module for independently averaging each column in input data.
 *
 * Computes the average and standard deviation independently for each column in
 * the input data.  Multipoint data, multiple data sets, and missing data
 * points are all supported.
 * The average is always calculated over all frames and data points for a
 * column.
 *
 * Output data contains a column for each data set in the input data, and a
 * frame for each column in the input data.  If different data sets have
 * different number of columns, the frame count accomodates the largest data
 * set.  Other columns are padded with zero values that are additionally marked
 * as missing.
 * Each value in the output data is the average of the corresponding
 * input column in the corresponding input data set.  The error value for each
 * value provides the standard deviation of the corresponding input column.
 * average(), standardDeviation(), and sampleCount() methods are also
 * provided for convenient access to these properties.
 *
 * The output data becomes available only after the input data has been
 * finished.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisDataAverageModule : public AbstractAnalysisArrayData,
                                  public AnalysisDataModuleSerial
{
    public:
        AnalysisDataAverageModule();
        virtual ~AnalysisDataAverageModule();

        using AbstractAnalysisArrayData::setXAxis;
        using AbstractAnalysisArrayData::setXAxisValue;

        /*! \brief
         * Sets the averaging to happen over entire data sets.
         *
         * If \p bDataSets is false (the default), the module averages each
         * column separately.  The output will have a column for each data set,
         * and a row for each column.
         *
         * If \p bDataSets is true, the module averages all values within
         * a single data set into a single average/standard deviation.
         * The output will have only one column, with one row for each data
         * set.
         */
        void setAverageDataSets(bool bDataSets);

        virtual int flags() const;

        virtual void dataStarted(AbstractAnalysisData *data);
        virtual void frameStarted(const AnalysisDataFrameHeader &header);
        virtual void pointsAdded(const AnalysisDataPointSetRef &points);
        virtual void frameFinished(const AnalysisDataFrameHeader &header);
        virtual void dataFinished();

        /*! \brief
         * Convenience access to the average of a data column.
         *
         * Note that the interpretation of the parameters follows their naming:
         * with \c setAverageDataSets(false), \p dataSet corresponds to a
         * column in the output, but with \c setAverageDataSets(false) it
         * corresponds to an output row.  In both cases, it selects the data
         * set; with \c setAverageDataSets(false), \p column should always be
         * zero as there is only one value per data set.
         */
        real average(int dataSet, int column) const;
        /*! \brief
         * Convenience access to the standard deviation of a data column.
         *
         * See average() for the interpretation of the parameters.
         */
        real standardDeviation(int dataSet, int column) const;
        /*! \brief
         * Access the number of samples for a data column.
         *
         * See average() for the interpretation of the parameters.
         */
        int sampleCount(int dataSet, int column) const;

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
                                       public AnalysisDataModuleSerial
{
    public:
        AnalysisDataFrameAverageModule();
        virtual ~AnalysisDataFrameAverageModule();

        virtual int frameCount() const;

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
