/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements gmx::AnalysisDataDisplacementModule.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#include "gmxpre.h"

#include "gromacs/analysisdata/modules/displacement.h"

#include <memory>
#include <vector>

#include "gromacs/analysisdata/abstractdata.h"
#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/datamodulemanager.h"
#include "gromacs/analysisdata/modules/histogram.h"
#include "gromacs/math/utilities.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
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
    int nmax;
    //! Maximum time for which the displacements are needed.
    real tmax;
    //! Number of dimensions per data point.
    int ndim;

    //! true if no frames have been read.
    bool bFirst;
    //! Stores the time of the first frame.
    real t0;
    //! Stores the time interval between frames.
    real dt;
    //! Stores the time of the current frame.
    real t;
    //! Stores the index in the store for the current positions.
    int ci;

    //! Maximum number of positions to store for a particle.
    int max_store;
    //! The total number of positions ever stored (can be larger than \p max_store).
    int nstored;
    //! Old values.
    real* oldval;
    //! The most recently calculated displacements.
    std::vector<AnalysisDataValue> currValues_;

    //! Histogram module for calculating MSD histograms, or NULL if not set.
    AnalysisDataBinAverageModule* histm;
};

AnalysisDataDisplacementModule::Impl::Impl() :
    nmax(0),
    tmax(0.0),
    ndim(3),
    bFirst(true),
    t0(0.0),
    dt(0.0),
    t(0.0),
    ci(0),
    max_store(-1),
    nstored(0),
    oldval(nullptr),
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

AnalysisDataDisplacementModule::AnalysisDataDisplacementModule() : _impl(new Impl())
{
    setMultipoint(true);
}


AnalysisDataDisplacementModule::~AnalysisDataDisplacementModule() {}


void AnalysisDataDisplacementModule::setMaxTime(real tmax)
{
    _impl->tmax = tmax;
}


void AnalysisDataDisplacementModule::setMSDHistogram(const AnalysisDataBinAverageModulePointer& histm)
{
    GMX_RELEASE_ASSERT(_impl->histm == nullptr, "Can only set MSD histogram once");
    _impl->histm = histm.get();
    addModule(histm);
}


AnalysisDataFrameRef AnalysisDataDisplacementModule::tryGetDataFrameInternal(int /*index*/) const
{
    return AnalysisDataFrameRef();
}


bool AnalysisDataDisplacementModule::requestStorageInternal(int /*nframes*/)
{
    return false;
}


int AnalysisDataDisplacementModule::flags() const
{
    return efAllowMulticolumn;
}


void AnalysisDataDisplacementModule::dataStarted(AbstractAnalysisData* data)
{
    if (data->columnCount() % _impl->ndim != 0)
    {
        GMX_THROW(APIError("Data has incorrect number of columns"));
    }
    _impl->nmax = data->columnCount();
    snew(_impl->oldval, _impl->nmax);
    _impl->ci = -_impl->nmax;

    int ncol = _impl->nmax / _impl->ndim + 1;
    _impl->currValues_.reserve(ncol);
    setColumnCount(0, ncol);
}


void AnalysisDataDisplacementModule::frameStarted(const AnalysisDataFrameHeader& header)
{
    // Initialize times.
    if (_impl->bFirst)
    {
        _impl->t0 = header.x();
    }
    else if (_impl->dt <= 0)
    {
        _impl->dt = header.x() - _impl->t0;
        if (_impl->dt < 0 || gmx_within_tol(_impl->dt, 0.0, GMX_REAL_EPS))
        {
            GMX_THROW(APIError("Identical or decreasing frame times"));
        }
    }
    else
    {
        if (!gmx_within_tol(header.x() - _impl->t, _impl->dt, GMX_REAL_EPS))
        {
            GMX_THROW(APIError("Frames not evenly spaced"));
        }
    }
    _impl->t = header.x();

    // Allocate memory for all the positions once it is possible.
    if (_impl->max_store == -1 && !_impl->bFirst)
    {
        _impl->max_store = _impl->nmax * static_cast<int>(_impl->tmax / _impl->dt + 1);
        srenew(_impl->oldval, _impl->max_store);
    }

    // Increment the index where current positions are stored.
    _impl->ci += _impl->nmax;
    if (_impl->ci >= _impl->max_store)
    {
        _impl->ci = 0;
    }

    /*
        for (int i = 0; i < _impl->nmax; ++i)
        {
            _impl->p[_impl->ci + i].bPres = false;
        }
     */
    _impl->nstored++;
    _impl->bFirst = false;
}


void AnalysisDataDisplacementModule::pointsAdded(const AnalysisDataPointSetRef& points)
{
    if (points.firstColumn() % _impl->ndim != 0 || points.columnCount() % _impl->ndim != 0)
    {
        GMX_THROW(APIError("Partial data points"));
    }
    for (int i = 0; i < points.columnCount(); ++i)
    {
        _impl->oldval[_impl->ci + points.firstColumn() + i] = points.y(i);
    }
}


void AnalysisDataDisplacementModule::frameFinished(const AnalysisDataFrameHeader& /*header*/)
{
    if (_impl->nstored <= 1)
    {
        return;
    }

    if (_impl->nstored == 2)
    {
        if (_impl->histm)
        {
            _impl->histm->init(
                    histogramFromBins(0, _impl->max_store / _impl->nmax, _impl->dt).integerBins());
        }
        moduleManager().notifyDataStart(this);
    }
    AnalysisDataFrameHeader header(_impl->nstored - 2, _impl->t, 0);
    moduleManager().notifyFrameStart(header);

    for (int i = _impl->ci - _impl->nmax, step = 1; step < _impl->nstored && i != _impl->ci;
         i -= _impl->nmax, ++step)
    {
        if (i < 0)
        {
            i += _impl->max_store;
        }
        _impl->currValues_.clear();
        _impl->currValues_.emplace_back(step * _impl->dt);
        int k = 1;
        for (int j = 0; j < _impl->nmax; j += _impl->ndim, ++k)
        {
            real dist2 = 0.0;

            for (int d = 0; d < _impl->ndim; ++d)
            {
                real displ = _impl->oldval[_impl->ci + j + d] - _impl->oldval[i + j + d];
                dist2 += displ * displ;
            }
            _impl->currValues_.emplace_back(dist2);
        }
        moduleManager().notifyPointsAdd(AnalysisDataPointSetRef(header, _impl->currValues_));
    }

    moduleManager().notifyFrameFinish(header);
}


void AnalysisDataDisplacementModule::dataFinished()
{
    if (_impl->nstored >= 2)
    {
        moduleManager().notifyDataFinish();
    }
}

} // namespace gmx
