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
 * Declares gmx::AnalysisDataProxy.
 *
 * This header is only meant for internal use of the gmx::AbstractAnalysisData
 * class to implement modules that handle only a subset of columns.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_DATAPROXY_H
#define GMX_ANALYSISDATA_DATAPROXY_H

#include "abstractdata.h"
#include "datamodule.h"

namespace gmx
{

/*! \internal \brief
 * Internal implementation class used to implement column modules.
 *
 * \ingroup module_analysisdata
 */
class AnalysisDataProxy : public AbstractAnalysisData,
                          public AnalysisDataModuleInterface
{
    public:
        /*! \brief
         * Creates a proxy object that only presents certain columns.
         *
         * \param[in] col   First column to present.
         * \param[in] span  Number of columns to present.
         * \param[in] data  Data object that should be wrapped.
         */
        AnalysisDataProxy(int col, int span, AbstractAnalysisData *data);

        virtual int frameCount() const;
        virtual bool getDataWErr(int index, real *x, real *dx,
                                 const real **y, const real **dy,
                                 const bool **missing = 0) const;
        virtual bool requestStorage(int nframes = -1);

        virtual int flags() const;

        virtual void dataStarted(AbstractAnalysisData *data);
        virtual void frameStarted(const AnalysisDataFrameHeader &frame);
        virtual void pointsAdded(const AnalysisDataPointSetRef &points);
        virtual void frameFinished(const AnalysisDataFrameHeader &header);
        virtual void dataFinished();

    private:
        AbstractAnalysisData   &_source;
        int                     _col;
        int                     _span;

        // Copy and assign disallowed by base.
};

} // namespace gmx

#endif
