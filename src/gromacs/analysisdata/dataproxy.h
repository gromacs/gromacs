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
 * This header is only meant for internal use to implement
 * gmx::AbstractAnalysisData::setColumnModule().
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
 * This class serves as a proxy between AbstractAnalysisData and the attached
 * AnalysisDataModuleInterface object.  For each notification that
 * AbstractAnalysisData sends, it maps it such that only the relevant columns
 * are visible to the AnalysisDataModuleInterface.  Similarly, it implements
 * the frame access methods of AbstractAnalysisData such that only the relevant
 * columns are returned.
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
         *
         * Does not throw.
         */
        AnalysisDataProxy(int col, int span, AbstractAnalysisData *data);

        virtual int flags() const;

        virtual void dataStarted(AbstractAnalysisData *data);
        virtual void frameStarted(const AnalysisDataFrameHeader &frame);
        virtual void pointsAdded(const AnalysisDataPointSetRef &points);
        virtual void frameFinished(const AnalysisDataFrameHeader &header);
        virtual void dataFinished();

    private:
        virtual AnalysisDataFrameRef tryGetDataFrameInternal(int index) const;
        virtual bool requestStorageInternal(int nframes);

        AbstractAnalysisData   &_source;
        int                     _col;
        int                     _span;

        // Copy and assign disallowed by base.
};

} // namespace gmx

#endif
