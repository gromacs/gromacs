/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * Implements classes in dataframe.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#include "dataframe.h"

#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/********************************************************************
 * AnalysisDataFrameHeader
 */

AnalysisDataFrameHeader::AnalysisDataFrameHeader()
    : index_(-1), x_(0.0), dx_(0.0)
{
}


AnalysisDataFrameHeader::AnalysisDataFrameHeader(int index, real x, real dx)
    : index_(index), x_(x), dx_(dx)
{
    GMX_ASSERT(index >= 0, "Invalid frame index");
}


/********************************************************************
 * AnalysisDataPointSetRef
 */

AnalysisDataPointSetRef::AnalysisDataPointSetRef(
        const AnalysisDataFrameHeader &header, int firstColumn,
        const AnalysisDataValuesRef &values)
    : header_(header), firstColumn_(firstColumn), values_(values)
{
    GMX_ASSERT(header_.isValid(),
               "Invalid point set reference should not be constructed");
    GMX_ASSERT(firstColumn >= 0, "Invalid first column");
}


AnalysisDataPointSetRef::AnalysisDataPointSetRef(
        const AnalysisDataFrameHeader        &header,
        const std::vector<AnalysisDataValue> &values)
    : header_(header), firstColumn_(0), values_(values.begin(), values.end())
{
    GMX_ASSERT(header_.isValid(),
               "Invalid point set reference should not be constructed");
}


AnalysisDataPointSetRef::AnalysisDataPointSetRef(
        const AnalysisDataPointSetRef &points, int firstColumn, int columnCount)
    : header_(points.header()), firstColumn_(0)
{
    GMX_ASSERT(firstColumn >= 0, "Invalid first column");
    GMX_ASSERT(columnCount >= 0, "Invalid column count");
    if (points.lastColumn() < firstColumn
        || points.firstColumn() >= firstColumn + columnCount
        || columnCount == 0)
    {
        return;
    }
    AnalysisDataValuesRef::const_iterator begin = points.values().begin();
    int newFirstColumn = firstColumn - points.firstColumn();
    if (newFirstColumn > 0)
    {
        // Offset pointer if the first column is not the first in points.
        begin         += newFirstColumn;
        newFirstColumn = 0;
    }
    else
    {
        // Take into account if first column is before the first in points.
        columnCount -= -newFirstColumn;
    }
    // Decrease column count if there are not enough columns in points.
    AnalysisDataValuesRef::const_iterator end = begin + columnCount;
    if (newFirstColumn + columnCount > points.columnCount())
    {
        end = points.values().end();
    }
    values_ = AnalysisDataValuesRef(begin, end);
}


bool AnalysisDataPointSetRef::allPresent() const
{
    AnalysisDataValuesRef::const_iterator i;
    for (i = values_.begin(); i != values_.end(); ++i)
    {
        if (!i->isPresent())
        {
            return false;
        }
    }
    return true;
}


/********************************************************************
 * AnalysisDataFrameRef
 */

AnalysisDataFrameRef::AnalysisDataFrameRef()
{
}


AnalysisDataFrameRef::AnalysisDataFrameRef(
        const AnalysisDataFrameHeader &header,
        const AnalysisDataValuesRef   &values)
    : header_(header), values_(values)
{
}


AnalysisDataFrameRef::AnalysisDataFrameRef(
        const AnalysisDataFrameHeader        &header,
        const std::vector<AnalysisDataValue> &values)
    : header_(header), values_(values.begin(), values.end())
{
}


AnalysisDataFrameRef::AnalysisDataFrameRef(
        const AnalysisDataFrameRef &frame, int firstColumn, int columnCount)
    : header_(frame.header()), values_(&frame.values_[firstColumn], columnCount)
{
    GMX_ASSERT(firstColumn >= 0, "Invalid first column");
    GMX_ASSERT(columnCount >= 0, "Invalid column count");
    GMX_ASSERT(firstColumn + columnCount <= frame.columnCount(),
               "Invalid last column");
}


bool AnalysisDataFrameRef::allPresent() const
{
    GMX_ASSERT(isValid(), "Invalid data frame accessed");
    AnalysisDataValuesRef::const_iterator i;
    for (i = values_.begin(); i != values_.end(); ++i)
    {
        if (!i->isPresent())
        {
            return false;
        }
    }
    return true;
}

} // namespace gmx
