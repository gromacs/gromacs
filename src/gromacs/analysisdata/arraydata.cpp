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

#include <algorithm>

#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

AbstractAnalysisArrayData::AbstractAnalysisArrayData()
    : _nrows(0), _xstart(0.0), _xstep(1.0), _bReady(false)
{
}

AbstractAnalysisArrayData::~AbstractAnalysisArrayData()
{
}


AnalysisDataFrameRef
AbstractAnalysisArrayData::tryGetDataFrameInternal(int index) const
{
    if (!isAllocated())
    {
        return AnalysisDataFrameRef();
    }
    std::vector<AnalysisDataValue>::const_iterator begin
        = _value.begin() + index * columnCount();
    return AnalysisDataFrameRef(
                AnalysisDataFrameHeader(index, xvalue(index), 0.0),
                AnalysisDataValuesRef(begin, begin + columnCount()));
}


bool
AbstractAnalysisArrayData::requestStorageInternal(int /*nframes*/)
{
    return true;
}


void
AbstractAnalysisArrayData::setColumnCount(int ncols)
{
    GMX_RELEASE_ASSERT(!isAllocated(),
                       "Cannot change column count after data has been allocated");
    AbstractAnalysisData::setColumnCount(ncols);
}


void
AbstractAnalysisArrayData::setRowCount(int nrows)
{
    GMX_RELEASE_ASSERT(nrows > 0, "Invalid number of rows");
    GMX_RELEASE_ASSERT(!isAllocated(),
                       "Cannot change row count after data has been allocated");
    _nrows = nrows;
}


void
AbstractAnalysisArrayData::allocateValues()
{
    GMX_RELEASE_ASSERT(!isAllocated(), "Can only allocate values once");
    GMX_RELEASE_ASSERT(rowCount() > 0 && columnCount() > 0,
                       "Row and column counts must be set before allocating values");
    _value.resize(rowCount() * columnCount());
    std::vector<AnalysisDataValue>::iterator i;
    for (i = _value.begin(); i != _value.end(); ++i)
    {
        i->setValue(0.0);
    }
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
    GMX_RELEASE_ASSERT(isAllocated(), "There must be some data");
    if (_bReady)
    {
        return;
    }
    _bReady = true;

    std::vector<AnalysisDataValue>::const_iterator valueIter = _value.begin();
    notifyDataStart();
    for (int i = 0; i < rowCount(); ++i, valueIter += columnCount())
    {
        AnalysisDataFrameHeader header(i, xvalue(i), 0);
        notifyFrameStart(header);
        notifyPointsAdd(AnalysisDataPointSetRef(header, 0,
                            AnalysisDataValuesRef(valueIter,
                                                  valueIter + columnCount())));
        notifyFrameFinish(header);
    }
    notifyDataFinish();
}


void
AbstractAnalysisArrayData::copyContents(const AbstractAnalysisArrayData *src,
                                        AbstractAnalysisArrayData *dest)
{
    GMX_RELEASE_ASSERT(src->isAllocated(), "Source data must not be empty");
    GMX_RELEASE_ASSERT(!dest->isAllocated(),
                       "Destination data must not be allocated");
    dest->setColumnCount(src->columnCount());
    dest->setRowCount(src->rowCount());
    dest->allocateValues();
    dest->setXAxis(src->xstart(), src->xstep());
    std::copy(src->_value.begin(), src->_value.end(), dest->_value.begin());
}

} // namespace gmx
