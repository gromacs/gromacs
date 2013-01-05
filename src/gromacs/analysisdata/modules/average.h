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
/*! \file
 * \brief
 * Declares gmx::AnalysisDataAverageModule.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
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
 * Data module for simple averaging of columns.
 *
 * Output data contains a frame for each column of input data.
 * There are two columns: the average and standard deviation of
 * that column.
 * The data becomes available only after the original data has been
 * finished.
 *
 * Multipoint data and missing data points are both supported. The average
 * is always calculated over all data points present in a column.
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
        std::vector<int>        nsamples_;

        // Copy and assign disallowed by base.
};

//! Smart pointer to manage an AnalysisDataAverageModule object.
typedef boost::shared_ptr<AnalysisDataAverageModule>
    AnalysisDataAverageModulePointer;

/*! \brief
 * Data module for averaging of columns for each frame.
 *
 * Output data has the same number of frames as the input data, but only one
 * column.
 * Each frame in the output contains the average of the column values in the
 * corresponding frame of the input data.
 *
 * Multipoint data and missing data points are both supported. The average
 * is always calculated over all data points present in a column.
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
