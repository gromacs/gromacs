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
 * Implements gmx::AnalysisDataDisplacementModule.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#include "gromacs/analysisdata/modules/displacement.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

// Legacy include.
#include "smalloc.h"

#include "gromacs/basicmath.h"
#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/modules/histogram.h"
#include "gromacs/fatalerror/exceptions.h"
#include "gromacs/fatalerror/gmxassert.h"

#include "displacement-impl.h"

namespace gmx
{

/********************************************************************
 * AnalysisDataDisplacementModule::Impl
 */

AnalysisDataDisplacementModule::Impl::Impl()
    : nmax(0), tmax(0.0), ndim(3),
      bFirst(true), t0(0.0), dt(0.0), t(0.0),
      max_store(-1), nstored(0), oldval(NULL), currd(NULL),
      histm(NULL)
{
}

AnalysisDataDisplacementModule::Impl::~Impl()
{
    sfree(oldval);
    sfree(currd);
}

/********************************************************************
 * AnalysisDataDisplacementModule
 */

AnalysisDataDisplacementModule::AnalysisDataDisplacementModule()
    : _impl(new Impl())
{
    setMultipoint(true);
}


AnalysisDataDisplacementModule::~AnalysisDataDisplacementModule()
{
    delete _impl;
}


void
AnalysisDataDisplacementModule::setMaxTime(real tmax)
{
    _impl->tmax = tmax;
}


void
AnalysisDataDisplacementModule::setMSDHistogram(AnalysisDataBinAverageModule *histm)
{
    GMX_RELEASE_ASSERT(!_impl->histm, "Can only set MSD histogram once");
    _impl->histm = histm;
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
    if (data->columnCount() % _impl->ndim != 0)
    {
        GMX_THROW(APIError("Data has incorrect number of columns"));
    }
    _impl->nmax = data->columnCount();
    snew(_impl->oldval, _impl->nmax);
    _impl->ci = -_impl->nmax;

    int ncol = _impl->nmax / _impl->ndim + 1;
    snew(_impl->currd, ncol);
    setColumnCount(ncol);
}


void
AnalysisDataDisplacementModule::frameStarted(const AnalysisDataFrameHeader &header)
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
        _impl->max_store = _impl->nmax * (int)(_impl->tmax/_impl->dt + 1);
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


void
AnalysisDataDisplacementModule::pointsAdded(const AnalysisDataPointSetRef &points)
{
    if (points.firstColumn() % _impl->ndim != 0
        || points.columnCount() % _impl->ndim != 0)
    {
        GMX_THROW(APIError("Partial data points"));
    }
    for (int i = 0; i < points.columnCount(); ++i)
    {
        _impl->oldval[_impl->ci + points.firstColumn() + i] = points.y(i);
    }
}


void
AnalysisDataDisplacementModule::frameFinished(const AnalysisDataFrameHeader & /*header*/)
{
    if (_impl->nstored <= 1)
    {
        return;
    }

    int step, i;
    int rc;

    if (_impl->nstored == 2)
    {
        if (_impl->histm)
        {
            _impl->histm->init(histogramFromBins(0, _impl->max_store / _impl->nmax,
                                                 _impl->dt).integerBins());
        }
        notifyDataStart();
    }
    AnalysisDataFrameHeader header(_impl->nstored - 2, _impl->t, 0);
    notifyFrameStart(header);

    for (i = _impl->ci - _impl->nmax, step = 1;
         step < _impl->nstored && i != _impl->ci;
         i -= _impl->nmax, ++step)
    {
        if (i < 0)
        {
            i += _impl->max_store;
        }
        _impl->currd[0] = step * _impl->dt;
        int k = 1;
        for (int j = 0; j < _impl->nmax; j += _impl->ndim, ++k)
        {
            real dist2 = 0.0;

            for (int d = 0; d < _impl->ndim; ++d)
            {
                dist2 += sqr(_impl->oldval[_impl->ci + j + d]
                             - _impl->oldval[i + j + d]);
            }
            _impl->currd[k] = dist2;
        }
        notifyPointsAdd(AnalysisDataPointSetRef(
                header, 0, k, _impl->currd, NULL, NULL));
    }

    notifyFrameFinish(header);
}


void
AnalysisDataDisplacementModule::dataFinished()
{
    if (_impl->nstored >= 2)
    {
        notifyDataFinish();
    }
}

} // namespace gmx
