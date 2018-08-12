/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2015,2018, by the GROMACS development team, led by
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
 * Declares gmx::AnalysisDataProxy.
 *
 * This header is only meant for internal use to implement
 * gmx::AbstractAnalysisData::setColumnModule().
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_DATAPROXY_H
#define GMX_ANALYSISDATA_DATAPROXY_H

#include "gromacs/analysisdata/abstractdata.h"
#include "gromacs/analysisdata/datamodule.h"

namespace gmx
{

/*! \internal
 * \brief
 * Internal implementation class used to implement column modules.
 *
 * This class serves as a proxy between AbstractAnalysisData and the attached
 * IAnalysisDataModule object.  For each notification that
 * AbstractAnalysisData sends, it maps it such that only the relevant columns
 * are visible to the IAnalysisDataModule.  Similarly, it implements
 * the frame access methods of AbstractAnalysisData such that only the relevant
 * columns are returned.
 *
 * \ingroup module_analysisdata
 */
class AnalysisDataProxy : public AbstractAnalysisData,
                          public IAnalysisDataModule
{
    public:
        /*! \brief
         * Creates a proxy object that only presents certain columns.
         *
         * \param[in] firstColumn  First column to present.
         * \param[in] columnSpan   Number of columns to present.
         * \param[in] data         Data object that should be wrapped.
         *
         * Does not throw.
         */
        AnalysisDataProxy(int firstColumn, int columnSpan,
                          AbstractAnalysisData *data);

        int frameCount() const override;

        int flags() const override;

        void dataStarted(AbstractAnalysisData *data) override;
        bool parallelDataStarted(
            AbstractAnalysisData              *data,
            const AnalysisDataParallelOptions &options) override;
        void frameStarted(const AnalysisDataFrameHeader &frame) override;
        void pointsAdded(const AnalysisDataPointSetRef &points) override;
        void frameFinished(const AnalysisDataFrameHeader &header) override;
        void frameFinishedSerial(int frameIndex) override;
        void dataFinished() override;

    private:
        AnalysisDataFrameRef tryGetDataFrameInternal(int index) const override;
        bool requestStorageInternal(int nframes) override;

        AbstractAnalysisData   &source_;
        int                     firstColumn_;
        int                     columnSpan_;
        bool                    bParallel_;

        // Copy and assign disallowed by base.
};

} // namespace gmx

#endif
