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
 * Implements classes in arraydata.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#include "gromacs/analysisdata/arraydata.h"

#include <cassert>

// Legacy header.
#include "smalloc.h"

#include "gromacs/fatalerror/fatalerror.h"

namespace gmx
{

AbstractAnalysisArrayData::AbstractAnalysisArrayData()
    : _nrows(0), _value(NULL), _xstart(0.0), _xstep(1.0), _bReady(false)
{
}

AbstractAnalysisArrayData::~AbstractAnalysisArrayData()
{
    sfree(_value);
}


int
AbstractAnalysisArrayData::frameCount() const
{
    return _bReady ? _nrows : 0;
}


int
AbstractAnalysisArrayData::getDataWErr(int index, real *x, real *dx,
                                       const real **y, const real **dy,
                                       const bool **present) const
{
    if (index < 0)
    {
        index += _nrows;
        if (index < 0)
        {
            GMX_ERROR(eeInvalidValue, "Frame index out of range");
        }
    }
    if (index >= frameCount())
    {
        GMX_ERROR(eeInvalidValue, "Frame index out of range");
    }
    if (x)
    {
        *x = _xstart + index * _xstep;
    }
    if (dx)
    {
        *dx = 0.0;
    }
    if (y)
    {
        *y = _value + (index * columnCount());
    }
    if (dy)
    {
        // TODO: Implement
        *dy = NULL;
    }
    if (present)
    {
        // TODO: Implement
        *present = NULL;
    }
    return 0;
}


int
AbstractAnalysisArrayData::requestStorage(int nframes)
{
    return 0;
}


void
AbstractAnalysisArrayData::setColumnCount(int ncols)
{
    assert(!_value);
    AbstractAnalysisData::setColumnCount(ncols);
}


void
AbstractAnalysisArrayData::setRowCount(int nrows)
{
    assert(nrows > 0);
    assert(!_value && columnCount() > 0);
    _nrows = nrows;
    snew(_value, _nrows * columnCount());
}


void
AbstractAnalysisArrayData::setXAxis(real start, real step)
{
    assert(!_bReady);
    _xstart = start;
    _xstep = step;
}


int
AbstractAnalysisArrayData::valuesReady()
{
    assert(columnCount() > 0 && _nrows > 0);
    if (_bReady)
    {
        return 0;
    }
    _bReady = true;

    int rc = notifyDataStart();
    if (rc != 0)
    {
        return rc;
    }
    for (int i = 0; i < _nrows; ++i)
    {
        rc = notifyFrameStart(_xstart + i * _xstep, 0);
        if (rc != 0)
        {
            return rc;
        }
        rc = notifyPointsAdd(0, columnCount(), _value + (i * columnCount()),
                             NULL, NULL);
        if (rc != 0)
        {
            return rc;
        }
        rc = notifyFrameFinish();
        if (rc != 0)
        {
            return rc;
        }
    }
    return notifyDataFinish();
}


void
AbstractAnalysisArrayData::copyContents(const AbstractAnalysisArrayData *src,
                                        AbstractAnalysisArrayData *dest)
{
    assert(src->columnCount() > 0 && src->_nrows > 0 && src->_value);
    assert(!dest->_value);
    dest->setColumnCount(src->columnCount());
    dest->setRowCount(src->_nrows);
    dest->setXAxis(src->_xstart, src->_xstep);
    for (int i = 0; i < src->_nrows * src->columnCount(); ++i)
    {
        dest->_value[i] = src->_value[i];
    }
}

} // namespace gmx
