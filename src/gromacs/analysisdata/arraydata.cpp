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

// Legacy header.
#include "smalloc.h"

#include "gromacs/fatalerror/exceptions.h"
#include "gromacs/fatalerror/gmxassert.h"

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


bool
AbstractAnalysisArrayData::getDataWErr(int index, real *x, real *dx,
                                       const real **y, const real **dy,
                                       const bool **present) const
{
    if (index < 0)
    {
        index += _nrows;
        if (index < 0)
        {
            return false;
        }
    }
    if (index >= frameCount())
    {
        return false;
    }
    if (x != NULL)
    {
        *x = _xstart + index * _xstep;
    }
    if (dx != NULL)
    {
        *dx = 0.0;
    }
    if (y != NULL)
    {
        *y = _value + (index * columnCount());
    }
    if (dy != NULL)
    {
        // TODO: Implement
        *dy = NULL;
    }
    if (present != NULL)
    {
        // TODO: Implement
        *present = NULL;
    }
    return true;
}


bool
AbstractAnalysisArrayData::requestStorage(int /*nframes*/)
{
    return true;
}


void
AbstractAnalysisArrayData::setColumnCount(int ncols)
{
    GMX_RELEASE_ASSERT(!_value,
                       "Cannot change column count after data has been allocated");
    AbstractAnalysisData::setColumnCount(ncols);
}


void
AbstractAnalysisArrayData::setRowCount(int nrows)
{
    GMX_RELEASE_ASSERT(nrows > 0, "Invalid number of rows");
    GMX_RELEASE_ASSERT(!_value,
                       "Cannot change row count after data has been allocated");
    GMX_RELEASE_ASSERT(columnCount() > 0, "Column count must be set before row count");
    _nrows = nrows;
    snew(_value, _nrows * columnCount());
}


void
AbstractAnalysisArrayData::setXAxis(real start, real step)
{
    GMX_RELEASE_ASSERT(!_bReady, "X axis cannot be set after data is finished");
    _xstart = start;
    _xstep = step;
}


void
AbstractAnalysisArrayData::valuesReady()
{
    GMX_RELEASE_ASSERT(columnCount() > 0 && _nrows > 0 && _value,
                       "There must be some data");
    if (_bReady)
    {
        return;
    }
    _bReady = true;

    notifyDataStart();
    for (int i = 0; i < _nrows; ++i)
    {
        notifyFrameStart(_xstart + i * _xstep, 0);
        notifyPointsAdd(0, columnCount(), _value + (i * columnCount()),
                        NULL, NULL);
        notifyFrameFinish();
    }
    notifyDataFinish();
}


void
AbstractAnalysisArrayData::copyContents(const AbstractAnalysisArrayData *src,
                                        AbstractAnalysisArrayData *dest)
{
    GMX_RELEASE_ASSERT(src->columnCount() > 0 && src->_nrows > 0 && src->_value,
                       "Source data must not be empty");
    GMX_RELEASE_ASSERT(!dest->_value, "Destination data must not be allocated");
    dest->setColumnCount(src->columnCount());
    dest->setRowCount(src->_nrows);
    dest->setXAxis(src->_xstart, src->_xstep);
    for (int i = 0; i < src->_nrows * src->columnCount(); ++i)
    {
        dest->_value[i] = src->_value[i];
    }
}

} // namespace gmx
