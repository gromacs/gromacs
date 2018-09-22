/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2016,2017,2018, by the GROMACS development team, led by
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
 * Implements gmx::AnalysisDataDisplacementModule.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#include "gmxpre.h"

#include "displacement.h"

#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/datamodulemanager.h"
#include "gromacs/analysisdata/modules/histogram.h"
#include "gromacs/math/utilities.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

/********************************************************************
 * AnalysisDataDisplacementModule::Impl
 */

/*! \internal \brief
 * Private implementation class for AnalysisDataDisplacementModule.
 *
 * \ingroup module_analysisdata
 */
class AnalysisDataDisplacementModule::Impl
{
    public:
        Impl();
        ~Impl();

        //! Maximum number of particles for which the displacements are calculated.
        int                     nmax;
        //! Maximum time for which the displacements are needed.
        real                    tmax;
        //! Number of dimensions per data point.
        int                     ndim;

        //! true if no frames have been read.
        bool                    bFirst;
        //! Stores the time of the first frame.
        real                    t0;
        //! Stores the time interval between frames.
        real                    dt;
        //! Stores the time of the current frame.
        real                    t;
        //! Stores the index in the store for the current positions.
        int                     ci;

        //! Maximum number of positions to store for a particle.
        int                            max_store;
        //! The total number of positions ever stored (can be larger than \p max_store).
        int                            nstored;
        //! Old values.
        real                          *oldval;
        //! The most recently calculated displacements.
        std::vector<AnalysisDataValue> currValues_;

        //! Histogram module for calculating MSD histograms, or NULL if not set.
        AnalysisDataBinAverageModule *histm;
};

AnalysisDataDisplacementModule::Impl::Impl()
    : nmax(0), tmax(0.0), ndim(3),
      bFirst(true), t0(0.0), dt(0.0), t(0.0), ci(0),
      max_store(-1), nstored(0), oldval(nullptr),
      histm(nullptr)
{
}

AnalysisDataDisplacementModule::Impl::~Impl()
{
    sfree(oldval);
}

/********************************************************************
 * AnalysisDataDisplacementModule
 */

AnalysisDataDisplacementModule::AnalysisDataDisplacementModule()
    : impl_(new Impl())
{
    setMultipoint(true);
}


AnalysisDataDisplacementModule::~AnalysisDataDisplacementModule()
{
}


void
AnalysisDataDisplacementModule::setMaxTime(real tmax)
{
    impl_->tmax = tmax;
}


void
AnalysisDataDisplacementModule::setMSDHistogram(
        const AnalysisDataBinAverageModulePointer &histm)
{
    GMX_RELEASE_ASSERT(impl_->histm == nullptr, "Can only set MSD histogram once");
    impl_->histm = histm.get();
    addModule(histm);
}


AnalysisDataFrameRef
AnalysisDataDisplacementModule::tryGetDataFrameInternal(int /*index*/) const
{
    return AnalysisDataFrameRef();
}


bool
AnalysisDataDisplacementModule::requestStorageInternal(int /*nframes*/)
{
    return false;
}


int
AnalysisDataDisplacementModule::flags() const
{
    return efAllowMulticolumn;
}


void
AnalysisDataDisplacementModule::dataStarted(AbstractAnalysisData *data)
{
    if (data->columnCount() % impl_->ndim != 0)
    {
        GMX_THROW(APIError("Data has incorrect number of columns"));
    }
    impl_->nmax = data->columnCount();
    snew(impl_->oldval, impl_->nmax);
    impl_->ci = -impl_->nmax;

    int ncol = impl_->nmax / impl_->ndim + 1;
    impl_->currValues_.reserve(ncol);
    setColumnCount(0, ncol);
}


void
AnalysisDataDisplacementModule::frameStarted(const AnalysisDataFrameHeader &header)
{
    // Initialize times.
    if (impl_->bFirst)
    {
        impl_->t0 = header.x();
    }
    else if (impl_->dt <= 0)
    {
        impl_->dt = header.x() - impl_->t0;
        if (impl_->dt < 0 || gmx_within_tol(impl_->dt, 0.0, GMX_REAL_EPS))
        {
            GMX_THROW(APIError("Identical or decreasing frame times"));
        }
    }
    else
    {
        if (!gmx_within_tol(header.x() - impl_->t, impl_->dt, GMX_REAL_EPS))
        {
            GMX_THROW(APIError("Frames not evenly spaced"));
        }
    }
    impl_->t = header.x();

    // Allocate memory for all the positions once it is possible.
    if (impl_->max_store == -1 && !impl_->bFirst)
    {
        impl_->max_store = impl_->nmax * static_cast<int>(impl_->tmax/impl_->dt + 1);
        srenew(impl_->oldval, impl_->max_store);
    }

    // Increment the index where current positions are stored.
    impl_->ci += impl_->nmax;
    if (impl_->ci >= impl_->max_store)
    {
        impl_->ci = 0;
    }

/*
    for (int i = 0; i < impl_->nmax; ++i)
    {
        impl_->p[impl_->ci + i].bPres = false;
    }
 */
    impl_->nstored++;
    impl_->bFirst = false;
}


void
AnalysisDataDisplacementModule::pointsAdded(const AnalysisDataPointSetRef &points)
{
    if (points.firstColumn() % impl_->ndim != 0
        || points.columnCount() % impl_->ndim != 0)
    {
        GMX_THROW(APIError("Partial data points"));
    }
    for (int i = 0; i < points.columnCount(); ++i)
    {
        impl_->oldval[impl_->ci + points.firstColumn() + i] = points.y(i);
    }
}


void
AnalysisDataDisplacementModule::frameFinished(const AnalysisDataFrameHeader & /*header*/)
{
    if (impl_->nstored <= 1)
    {
        return;
    }

    int step, i;

    if (impl_->nstored == 2)
    {
        if (impl_->histm)
        {
            impl_->histm->init(histogramFromBins(0, impl_->max_store / impl_->nmax,
                                                 impl_->dt).integerBins());
        }
        moduleManager().notifyDataStart(this);
    }
    AnalysisDataFrameHeader header(impl_->nstored - 2, impl_->t, 0);
    moduleManager().notifyFrameStart(header);

    for (i = impl_->ci - impl_->nmax, step = 1;
         step < impl_->nstored && i != impl_->ci;
         i -= impl_->nmax, ++step)
    {
        if (i < 0)
        {
            i += impl_->max_store;
        }
        impl_->currValues_.clear();
        impl_->currValues_.emplace_back(step * impl_->dt);
        int k = 1;
        for (int j = 0; j < impl_->nmax; j += impl_->ndim, ++k)
        {
            real dist2 = 0.0;

            for (int d = 0; d < impl_->ndim; ++d)
            {
                real displ = impl_->oldval[impl_->ci + j + d]
                    - impl_->oldval[i + j + d];
                dist2 += displ * displ;
            }
            impl_->currValues_.emplace_back(dist2);
        }
        moduleManager().notifyPointsAdd(AnalysisDataPointSetRef(header, impl_->currValues_));
    }

    moduleManager().notifyFrameFinish(header);
}


void
AnalysisDataDisplacementModule::dataFinished()
{
    if (impl_->nstored >= 2)
    {
        moduleManager().notifyDataFinish();
    }
}

} // namespace gmx
