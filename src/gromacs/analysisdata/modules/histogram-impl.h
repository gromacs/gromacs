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
 * Declares private implementation classes for classes in histogram.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_MODULES_HISTOGRAM_IMPL_H
#define GMX_ANALYSISDATA_MODULES_HISTOGRAM_IMPL_H

#include <vector>

#include "histogram.h"

#include "../datastorage.h"

namespace gmx
{

namespace internal
{

/*! \internal \brief
 * Represents copies of average histograms.
 *
 * Methods in AbstractAverageHistogram that return new histogram instances
 * return objects of this class.
 * Initialization of values is handled in those methods.
 *
 * \ingroup module_analysisdata
 */
class StaticAverageHistogram : public AbstractAverageHistogram
{
    public:
        StaticAverageHistogram();
        //! Creates an average histogram module with defined bin parameters.
        explicit StaticAverageHistogram(const AnalysisHistogramSettings &settings);

        // Copy and assign disallowed by base.
};

/*! \internal \brief
 * Implements average histogram module that averages per-frame histograms.
 *
 * This class is used for accumulating average histograms in per-frame
 * histogram modules (those that use BasicHistogramImpl as their implementation
 * class).
 * There are two columns, first for the average and second for standard
 * deviation.
 *
 * \ingroup module_analysisdata
 */
class BasicAverageHistogramModule : public AbstractAverageHistogram,
                                    public AnalysisDataModuleInterface
{
    public:
        BasicAverageHistogramModule();
        //! Creates an average histogram module with defined bin parameters.
        explicit BasicAverageHistogramModule(const AnalysisHistogramSettings &settings);

        using AbstractAverageHistogram::init;

        virtual int flags() const;

        virtual void dataStarted(AbstractAnalysisData *data);
        virtual void frameStarted(const AnalysisDataFrameHeader &header);
        virtual void pointsAdded(const AnalysisDataPointSetRef &points);
        virtual void frameFinished(const AnalysisDataFrameHeader &header);
        virtual void dataFinished();

    private:
        //! Number of frames accumulated so far.
        int                     frameCount_;

        // Copy and assign disallowed by base.
};

//! Smart pointer to manage an BasicAverageHistogramModule object.
typedef boost::shared_ptr<BasicAverageHistogramModule>
    BasicAverageHistogramModulePointer;

/*! \internal \brief
 * Private implementation class for AnalysisDataSimpleHistogramModule and
 * AnalysisDataWeightedHistogramModule.
 *
 * \ingroup module_analysisdata
 */
class BasicHistogramImpl
{
    public:
        BasicHistogramImpl();
        //! Creates an histogram impl with defined bin parameters.
        explicit BasicHistogramImpl(const AnalysisHistogramSettings &settings);
        ~BasicHistogramImpl();

        /*! \brief
         * (Re)initializes the histogram from settings.
         */
        void init(const AnalysisHistogramSettings &settings);
        /*! \brief
         * Initializes data storage frame when a new frame starts.
         */
        void initFrame(AnalysisDataStorageFrame *frame);

        //! Storage implementation object.
        AnalysisDataStorage                  storage_;
        //! Settings for the histogram object.
        AnalysisHistogramSettings            settings_;
        //! Averager module.
        BasicAverageHistogramModulePointer   averager_;
};

} // namespace internal

} // namespace gmx

#endif
