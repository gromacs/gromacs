/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2015,2017, by the GROMACS development team, led by
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
 * Implements classes in arraydata.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#include "gmxpre.h"

#include "arraydata.h"

#include <algorithm>

#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/datamodulemanager.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

AbstractAnalysisArrayData::AbstractAnalysisArrayData()
    : rowCount_(0), pointSetInfo_(0, 0, 0, 0), xstep_(1.0),
      bUniformX_(true), bReady_(false)
{
    xvalue_.push_back(0);
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
    return AnalysisDataFrameRef(
            AnalysisDataFrameHeader(index, xvalue(index), 0.0),
            makeConstArrayRef(value_).
                subArray(index * columnCount(), columnCount()),
            constArrayRefFromArray(&pointSetInfo_, 1));
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
    AbstractAnalysisData::setColumnCount(0, ncols);
    pointSetInfo_ = AnalysisDataPointSetInfo(0, ncols, 0, 0);
}


void
AbstractAnalysisArrayData::setRowCount(int rowCount)
{
    GMX_RELEASE_ASSERT(rowCount > 0, "Invalid number of rows");
    GMX_RELEASE_ASSERT(!isAllocated(),
                       "Cannot change row count after data has been allocated");
    GMX_RELEASE_ASSERT(bUniformX_ || xvalue_.empty()
                       || rowCount == static_cast<int>(xvalue_.size()),
                       "X axis set with setXAxisValue() does not match the row count");
    xvalue_.resize(rowCount);
    if (bUniformX_ && rowCount > rowCount_)
    {
        for (int i = rowCount_; i < rowCount; ++i)
        {
            xvalue_[i] = xvalue_[0] + i * xstep_;
        }
    }
    rowCount_ = rowCount;
}


void
AbstractAnalysisArrayData::allocateValues()
{
    GMX_RELEASE_ASSERT(!isAllocated(), "Can only allocate values once");
    GMX_RELEASE_ASSERT(rowCount() > 0 && columnCount() > 0,
                       "Row and column counts must be set before allocating values");
    value_.resize(rowCount() * columnCount());
    std::vector<AnalysisDataValue>::iterator i;
    for (i = value_.begin(); i != value_.end(); ++i)
    {
        i->setValue(0.0);
    }
}


void
AbstractAnalysisArrayData::setXAxis(real start, real step)
{
    GMX_RELEASE_ASSERT(!bReady_, "X axis cannot be set after data is finished");
    xvalue_[0] = start;
    xstep_     = step;
    bUniformX_ = true;
    for (int i = 0; i < rowCount_; ++i)
    {
        xvalue_[i] = start + i * xstep_;
    }
}


void
AbstractAnalysisArrayData::setXAxisValue(int row, real value)
{
    GMX_RELEASE_ASSERT(!bReady_, "X axis cannot be set after data is finished");
    if (rowCount_ > 0)
    {
        GMX_RELEASE_ASSERT(row >= 0 && row < rowCount(), "Row index out of range");
    }
    else if (row >= static_cast<int>(xvalue_.size()))
    {
        xvalue_.resize(row + 1);
    }
    bUniformX_   = false;
    xstep_       = 0.0;
    xvalue_[row] = value;
}


void
AbstractAnalysisArrayData::valuesReady()
{
    GMX_RELEASE_ASSERT(isAllocated(), "There must be some data");
    if (bReady_)
    {
        return;
    }
    bReady_ = true;

    AnalysisDataModuleManager                     &modules   = moduleManager();
    modules.notifyDataStart(this);
    for (int i = 0; i < rowCount(); ++i)
    {
        AnalysisDataFrameHeader header(i, xvalue(i), 0);
        modules.notifyFrameStart(header);
        modules.notifyPointsAdd(
                AnalysisDataPointSetRef(
                        header, pointSetInfo_,
                        makeConstArrayRef(value_).
                            subArray(i*columnCount(), columnCount())));
        modules.notifyFrameFinish(header);
    }
    modules.notifyDataFinish();
}


void
AbstractAnalysisArrayData::copyContents(const AbstractAnalysisArrayData *src,
                                        AbstractAnalysisArrayData       *dest)
{
    GMX_RELEASE_ASSERT(src->isAllocated(), "Source data must not be empty");
    GMX_RELEASE_ASSERT(!dest->isAllocated(),
                       "Destination data must not be allocated");
    dest->setColumnCount(src->columnCount());
    dest->setRowCount(src->rowCount());
    dest->allocateValues();
    dest->xstep_     = src->xstep_;
    dest->bUniformX_ = src->bUniformX_;
    std::copy(src->xvalue_.begin(), src->xvalue_.end(), dest->xvalue_.begin());
    std::copy(src->value_.begin(), src->value_.end(), dest->value_.begin());
}

} // namespace gmx
