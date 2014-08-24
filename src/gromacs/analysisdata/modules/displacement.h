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
 * Declares gmx::AnalysisDataDisplacementModule.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_MODULES_DISPLACEMENT_H
#define GMX_ANALYSISDATA_MODULES_DISPLACEMENT_H

#include "gromacs/analysisdata/abstractdata.h"
#include "gromacs/analysisdata/datamodule.h"
#include "gromacs/utility/real.h"

namespace gmx
{

class AnalysisDataBinAverageModule;

/*! \brief
 * Data module for calculating displacements.
 *
 * Output data contains a frame for each frame in the input data except the
 * first one.  For each frame, there can be multiple points, each of which
 * describes displacement for a certain time difference ending that that frame.
 * The first column contains the time difference (backwards from the current
 * frame), and the remaining columns the sizes of the displacements.
 *
 * Current implementation is not very generic, but should be easy to extend.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisDataDisplacementModule : public AbstractAnalysisData,
                                       public AnalysisDataModuleSerial
{
    public:
        AnalysisDataDisplacementModule();
        virtual ~AnalysisDataDisplacementModule();

        /*! \brief
         * Sets the largest displacement time to be calculated.
         */
        void setMaxTime(real tmax);
        /*! \brief
         * Sets an histogram module that will receive a MSD histogram.
         *
         * If this function is not called, no histogram is calculated.
         */
        void setMSDHistogram(boost::shared_ptr<AnalysisDataBinAverageModule> histm);

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

        PrivateImplPointer<Impl> _impl;
};

//! Smart pointer to manage an AnalysisDataDisplacementModule object.
typedef boost::shared_ptr<AnalysisDataDisplacementModule>
    AnalysisDataDisplacementModulePointer;

} // namespace gmx

#endif
