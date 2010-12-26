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
 * Implements gmx::AnalysisDataProxy.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#include "dataproxy.h"

#include <cassert>

namespace gmx
{

AnalysisDataProxy::AnalysisDataProxy(int col, int span,
                                     AbstractAnalysisData *data)
    : _source(*data), _col(col), _span(span)
{
    assert(data);
    assert(col >= 0 && span > 0);
    setColumnCount(span);
    setMultipoint(_source.isMultipoint());
}


int
AnalysisDataProxy::frameCount() const
{
    return _source.frameCount();
}


int
AnalysisDataProxy::getDataWErr(int index, real *x, real *dx,
                               const real **y, const real **dy,
                               const bool **missing) const
{
    int rc = _source.getDataWErr(index, x, dx, y, dy, missing);
    if (rc == 0)
    {
        if (y && *y)
        {
            *y += _col;
        }
        if (dy && *dy)
        {
            *dy += _col;
        }
        if (missing && *missing)
        {
            *missing += _col;
        }
    }
    return rc;
}


int
AnalysisDataProxy::requestStorage(int nframes)
{
    return _source.requestStorage(nframes);
}


int
AnalysisDataProxy::flags() const
{
    return efAllowMultipoint | efAllowMulticolumn | efAllowMissing;
}


int
AnalysisDataProxy::dataStarted(AbstractAnalysisData *data)
{
    assert(data == &_source);
    assert(_col + _span <= _source.columnCount());
    return notifyDataStart();
}


int
AnalysisDataProxy::frameStarted(real x, real dx)
{
    return notifyFrameStart(x, dx);
}


int
AnalysisDataProxy::pointsAdded(real x, real dx, int firstcol, int n,
                               const real *y, const real *dy,
                               const bool *missing)
{
    if (firstcol + n <= _col || firstcol >= _col + _span)
    {
        return 0;
    }
    firstcol -= _col;
    if (firstcol < 0)
    {
        if (y)
        {
            y +=  -firstcol;
        }
        if (dy)
        {
            dy += -firstcol;
        }
        if (missing)
        {
            missing += -firstcol;
        }
        n -= -firstcol;
        firstcol = 0;
    }
    if (firstcol + n > _span)
    {
        n = _span - firstcol;
    }
    return notifyPointsAdd(firstcol, n, y, dy, missing);
}


int
AnalysisDataProxy::frameFinished()
{
    return notifyFrameFinish();
}


int
AnalysisDataProxy::dataFinished()
{
    return notifyDataFinish();
}

} // namespace gmx
