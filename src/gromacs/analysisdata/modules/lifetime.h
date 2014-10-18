/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
 * Declares gmx::AnalysisDataLifetimeModule.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_MODULES_LIFETIME_H
#define GMX_ANALYSISDATA_MODULES_LIFETIME_H

#include "gromacs/analysisdata/arraydata.h"
#include "gromacs/analysisdata/datamodule.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

/*! \brief
 * Data module for computing lifetime histograms for columns in input data.
 *
 * The input data set is treated as a boolean array: each value that is present
 * (AnalysisDataValue::isPresent() returns true) and is >0 is treated as
 * present, other values are treated as absent.
 * For each input data set, analyzes the columns to identify the intervals
 * where a column is continuously present.
 * Produces a histogram from the lengths of these intervals.
 * Input data should have frames with evenly spaced x values.
 *
 * Output data contains one column for each data set in the input data.
 * This column gives the lifetime histogram for the corresponding data set.
 * x axis in the output is spaced the same as in the input data, and extends
 * as long as required to cover all the histograms.
 * Histograms are padded with zeros as required to be of the same length.
 * setCumulative() can be used to alter the handling of subintervals in the
 * output histogram.
 *
 * The output data becomes available only after the input data has been
 * finished.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisDataLifetimeModule : public AbstractAnalysisArrayData,
                                   public AnalysisDataModuleSerial
{
    public:
        AnalysisDataLifetimeModule();
        virtual ~AnalysisDataLifetimeModule();

        /*! \brief
         * Sets a cumulative histogram mode.
         *
         * \param[in] bCumulative If true, all subintervals of a long
         *   interval are also explicitly added into the histogram.
         *
         * Does not throw.
         */
        void setCumulative(bool bCumulative);

        virtual int flags() const;

        virtual void dataStarted(AbstractAnalysisData *data);
        virtual void frameStarted(const AnalysisDataFrameHeader &header);
        virtual void pointsAdded(const AnalysisDataPointSetRef &points);
        virtual void frameFinished(const AnalysisDataFrameHeader &header);
        virtual void dataFinished();

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

//! Smart pointer to manage an AnalysisDataLifetimeModule object.
typedef boost::shared_ptr<AnalysisDataLifetimeModule>
    AnalysisDataLifetimeModulePointer;

} // namespace gmx

#endif
