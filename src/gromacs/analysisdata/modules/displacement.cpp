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

#include "gromacs/analysisdata/modules/histogram.h"
#include "gromacs/basicmath.h"
#include "gromacs/fatalerror/fatalerror.h"

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


int
AnalysisDataDisplacementModule::setMSDHistogram(AnalysisDataBinAverageModule *histm)
{
    assert(!_impl->histm);
    _impl->histm = histm;
    histm->setIgnoreMissing(true);
    return addModule(histm);
}


int
AnalysisDataDisplacementModule::frameCount() const
{
    return _impl->nstored > 1 ? _impl->nstored - 1 : 0;
}


int
AnalysisDataDisplacementModule::getDataWErr(int index, real *x, real *dx,
                                            const real **y, const real **dy,
                                            const bool **present) const
{
    return eedataDataNotAvailable;
}


int
AnalysisDataDisplacementModule::requestStorage(int nframes)
{
    GMX_ERROR(eeNotImplemented, "Displacement storage not supported");
}


int
AnalysisDataDisplacementModule::flags() const
{
    return efAllowMulticolumn;
}


int
AnalysisDataDisplacementModule::dataStarted(AbstractAnalysisData *data)
{
    if (data->columnCount() % _impl->ndim != 0)
    {
        GMX_ERROR(eeInvalidValue, "Data has incorrect number of columns");
    }
    _impl->nmax = data->columnCount();
    snew(_impl->oldval, _impl->nmax);
    _impl->ci = -_impl->nmax;

    int ncol = _impl->nmax / _impl->ndim + 1;
    snew(_impl->currd, ncol);
    setColumnCount(ncol);

    return 0;
}


int
AnalysisDataDisplacementModule::frameStarted(real x, real dx)
{
    // Initialize times.
    if (_impl->bFirst)
    {
        _impl->t0 = x;
    }
    else if (_impl->dt <= 0)
    {
        _impl->dt = x - _impl->t0;
        if (_impl->dt < 0 || gmx_within_tol(_impl->dt, 0.0, GMX_REAL_EPS))
        {
            GMX_ERROR(eeInvalidInput, "Identical or decreasing frame times");
        }
    }
    else
    {
        if (!gmx_within_tol(x - _impl->t, _impl->dt, GMX_REAL_EPS))
        {
            GMX_ERROR(eeInvalidInput, "Frames not evenly spaced");
        }
    }
    _impl->t = x;

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
    return 0;
}


int
AnalysisDataDisplacementModule::pointsAdded(real x, real dx, int firstcol, int n,
                                            const real *y, const real *dy,
                                            const bool *present)
{
    if (firstcol % _impl->ndim != 0 || n % _impl->ndim != 0)
    {
        GMX_ERROR(eeInvalidValue, "Partial data points");
    }
    for (int i = firstcol; i < firstcol + n; ++i)
    {
        _impl->oldval[_impl->ci + i] = y[i];
    }
    return 0;
}


int
AnalysisDataDisplacementModule::frameFinished()
{
    if (_impl->nstored <= 1)
    {
        return 0;
    }

    int step, i;
    int rc;

    if (_impl->nstored == 2)
    {
        if (_impl->histm)
        {
            _impl->histm->initNBins(0, _impl->dt,
                                    _impl->max_store / _impl->nmax, true);
        }
        rc = notifyDataStart();
        if (rc != 0)
        {
            _impl->nstored = 1;
            return rc;
        }
    }
    rc = notifyFrameStart(_impl->t, 0);
    if (rc != 0)
    {
        return rc;
    }

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
        rc = notifyPointsAdd(0, k, _impl->currd, NULL, NULL);
        if (rc != 0)
        {
            return rc;
        }
    }

    return notifyFrameFinish();
}


int
AnalysisDataDisplacementModule::dataFinished()
{
    if (_impl->nstored >= 2)
    {
        return notifyDataFinish();
    }
    return 0;
}

} // namespace gmx
