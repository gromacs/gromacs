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
 * Implements classes in histogram.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#include "gromacs/analysisdata/modules/histogram.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cassert>
#include <cmath>

// Legacy include.
#include "smalloc.h"

#include "gromacs/basicmath.h"
#include "gromacs/fatalerror/fatalerror.h"

namespace gmx
{

/********************************************************************
 * AbstractHistogramModule
 */

AbstractHistogramModule::AbstractHistogramModule()
    : _hist(NULL), _averager(NULL), _nbins(0)
{
}

AbstractHistogramModule::~AbstractHistogramModule()
{
    sfree(_hist);
}


int
AbstractHistogramModule::initNBins(real miny, real binw, int nbins,
                                   bool bIntegerBins)
{
    assert(nbins > 0 && binw > 0);
    if (bIntegerBins)
    {
        miny -= 0.5*binw;
    }
    _nbins    = nbins;
    _miny     = miny;
    _binwidth = binw;
    _maxy     = miny + nbins * binw;
    _invbw    = 1.0/binw;
    setColumnCount(_nbins);
    return 0;
}


int
AbstractHistogramModule::initRange(real miny, real maxy, real binw,
                                   bool bIntegerBins)
{
    assert(miny < maxy && binw > 0);
    if (bIntegerBins)
    {
        _nbins = (int)ceil((maxy - miny) / binw) + 1;
        miny -= 0.5 * binw;
        maxy  = miny + _nbins * binw;
    }
    else
    {
        miny = binw * floor(miny / binw);
        maxy = binw * ceil(maxy / binw);
        if (miny != 0)
        {
            miny -= binw;
        }
        maxy += binw;
        _nbins = (int)((maxy - miny) / binw + 0.5);
    }
    _miny     = miny;
    _maxy     = maxy;
    _binwidth = binw;
    _invbw    = 1.0/binw;
    setColumnCount(_nbins);
    return 0;
}


void
AbstractHistogramModule::setAll(bool bAll)
{
    _bAll = bAll;
}


int
AbstractHistogramModule::findBin(real y) const
{
    if (y < _miny)
    {
        return _bAll ? 0 : -1;
    }
    else if (y >= _maxy)
    {
        return _bAll ? _nbins-1 : -1;
    }
    return (int)((y - _miny) * _invbw);
}


HistogramAverageModule *
AbstractHistogramModule::averager()
{
    if (!_averager)
    {
        createAverager();
    }
    return _averager;
}


int
AbstractHistogramModule::flags() const
{
    return efAllowMultipoint;
}


int
AbstractHistogramModule::dataStarted(AbstractAnalysisData *data)
{
    if (!_averager)
    {
        createAverager();
    }
    _averager->setXAxis(_miny + 0.5 * _binwidth, _binwidth);
    snew(_hist, nbins());
    return startDataStore();
}


int
AbstractHistogramModule::frameStarted(real x, real dx)
{
    for (int i = 0; i < nbins(); ++i)
    {
        _hist[i] = 0.0;
    }
    return startNextFrame(x, dx);
}


int
AbstractHistogramModule::frameFinished()
{
    return storeThisFrame(_hist, NULL, NULL);
}


int
AbstractHistogramModule::dataFinished()
{
    return notifyDataFinish();
}


void
AbstractHistogramModule::createAverager()
{
    _averager = new HistogramAverageModule();
    addModule(_averager);
    _averager->setXAxis(_miny + 0.5 * _binwidth, _binwidth);
}


/********************************************************************
 * HistogramAverageModule
 */

HistogramAverageModule::HistogramAverageModule()
    : _nframes(0), _bIgnoreMissing(false)
{
    setColumnCount(2);
}


void
HistogramAverageModule::setIgnoreMissing(bool bIgnoreMissing)
{
    _bIgnoreMissing = bIgnoreMissing;
    setColumnCount(bIgnoreMissing ? 3 : 2);
}


int
HistogramAverageModule::flags() const
{
    return efAllowMulticolumn | efAllowMissing;
}


int
HistogramAverageModule::dataStarted(AbstractAnalysisData *data)
{
    setRowCount(data->columnCount());
    return 0;
}


int
HistogramAverageModule::frameStarted(real x, real dx)
{
    return 0;
}


int
HistogramAverageModule::pointsAdded(real x, real dx, int firstcol, int n,
                                    const real *y, const real *dy,
                                    const bool *present)
{
    for (int i = 0; i < n; ++i)
    {
        value(firstcol + i, 0) += y[i];
        value(firstcol + i, 1) += y[i] * y[i];
    }
    if (_bIgnoreMissing)
    {
        assert(present != NULL);
        for (int i = 0; i < n; ++i)
        {
            if (present[i])
            {
                value(firstcol + i, 2) += 1;
            }
        }
    }
    return 0;
}


int
HistogramAverageModule::frameFinished()
{
    ++_nframes;
    return 0;
}


int
HistogramAverageModule::dataFinished()
{
    for (int i = 0; i < rowCount(); ++i)
    {
        real ave = 0.0;
        real std = 0.0;
        if (_bIgnoreMissing)
        {
            if (value(i, 2) > 0)
            {
                ave = value(i, 0) / value(i, 2);
                std = sqrt(value(i, 1) / value(i, 2) - ave * ave);
            }
        }
        else
        {
            ave = value(i, 0) / _nframes;
            std = sqrt(value(i, 1) / _nframes - ave * ave);
        }
        setValue(i, 0, ave);
        setValue(i, 1, std);
    }
    return 0;
}


HistogramAverageModule *
HistogramAverageModule::resampleDoubleBinWidth(bool bIntegerBins) const
{
    HistogramAverageModule *dest = new HistogramAverageModule();
    int nbins = rowCount() / 2;
    dest->setRowCount(nbins);
    real minx = xstart() + xstep() / 2;
    if (bIntegerBins)
    {
        minx -= xstep();
    }
    dest->setXAxis(minx, xstep() * 2);

    int  i, j;
    for (i = j = 0; i < nbins; ++i)
    {
        real  v, ve;
        if (bIntegerBins && i == 0)
        {
            v  = value(0, 0);
            ve = sqr(value(0, 1));
            ++j;
        }
        else
        {
            v  =     value(j, 0)  +     value(j+1, 0);
            ve = sqr(value(j, 1)) + sqr(value(j+1, 1));
            j += 2;
        }
        ve = sqrt(ve);
        dest->setValue(i, 0, v);
        dest->setValue(i, 1, ve);
    }
    return dest;
}


HistogramAverageModule *
HistogramAverageModule::clone() const
{
    HistogramAverageModule *dest = new HistogramAverageModule();
    copyContents(this, dest);
    return dest;
}


void
HistogramAverageModule::normalizeProbability()
{
    real sum = 0;
    for (int i = 0; i < rowCount(); ++i)
    {
        sum += value(i, 0);
    }
    scale(1.0 / (sum * xstep()));
}


void
HistogramAverageModule::scale(real norm)
{
    for (int i = 0; i < rowCount(); ++i)
    {
        value(i, 0) *= norm;
        value(i, 1) *= norm;
    }
}


void
HistogramAverageModule::scaleVector(real norm[])
{
    for (int i = 0; i < rowCount(); ++i)
    {
        value(i, 0) *= norm[i];
        value(i, 1) *= norm[i];
    }
}


/********************************************************************
 * AnalysisDataSimpleHistogramModule
 */

int
AnalysisDataSimpleHistogramModule::pointsAdded(real /*x*/, real /*dx*/,
                                               int /*firstcol*/, int n,
                                               const real *y, const real * /*dy*/,
                                               const bool * /*present*/)
{
    for (int i = 0; i < n; ++i)
    {
        int bin = findBin(y[i]);
        _hist[bin] += 1;
    }
    return 0;
}


/********************************************************************
 * AnalysisDataWeightedHistogramModule
 */

int
AnalysisDataWeightedHistogramModule::flags() const
{
    return AbstractHistogramModule::flags() | efAllowMulticolumn;
}


int
AnalysisDataWeightedHistogramModule::pointsAdded(real x, real dx, int firstcol, int n,
                                                 const real *y, const real *dy,
                                                 const bool *present)
{
    if (firstcol != 0 || n < 2)
    {
        GMX_ERROR(eeInvalidValue, "Invalid data layout");
    }
    int bin = findBin(y[0]);
    for (int i = 1; i < n; ++i)
    {
        _hist[bin] += y[i];
    }
    return 0;
}


/********************************************************************
 * AnalysisDataBinAverageModule
 */

AnalysisDataBinAverageModule::AnalysisDataBinAverageModule()
    : _n(NULL), _present(NULL), _bIgnoreMissing(false)
{
}

AnalysisDataBinAverageModule::~AnalysisDataBinAverageModule()
{
    sfree(_n);
    sfree(_present);
}


void
AnalysisDataBinAverageModule::setIgnoreMissing(bool bIgnoreMissing)
{
    // Changes can only be made before there is data.
    assert(!_n);
    _bIgnoreMissing = bIgnoreMissing;
    averager()->setIgnoreMissing(bIgnoreMissing);
}


int
AnalysisDataBinAverageModule::flags() const
{
    return AbstractHistogramModule::flags() | efAllowMulticolumn;
}


int
AnalysisDataBinAverageModule::dataStarted(AbstractAnalysisData *data)
{
    snew(_n, nbins());
    snew(_present, nbins());
    return AbstractHistogramModule::dataStarted(data);
}


int
AnalysisDataBinAverageModule::frameStarted(real x, real dx)
{
    for (int i = 0; i < nbins(); ++i)
    {
        _n[i] = 0;
    }
    return AbstractHistogramModule::frameStarted(x, dx);
}


int
AnalysisDataBinAverageModule::pointsAdded(real x, real dx, int firstcol, int n,
                                          const real *y, const real *dy,
                                          const bool *present)
{
    if (firstcol != 0 || n < 2)
    {
        GMX_ERROR(eeInvalidValue, "Invalid data layout");
    }
    int bin = findBin(y[0]);
    for (int i = 1; i < n; ++i)
    {
        _hist[bin] += y[i];
    }
    _n[bin] += n - 1;
    return 0;
}


int
AnalysisDataBinAverageModule::frameFinished()
{
    for (int i = 0; i < nbins(); ++i)
    {
        _present[i] = (_n[i] > 0);
        if (_n[i] > 0)
        {
            _hist[i] /= _n[i];
        }
    }
    if (!_bIgnoreMissing)
    {
        return AbstractHistogramModule::frameFinished();
    }
    return storeThisFrame(_hist, NULL, _present);
}

} // namespace gmx
