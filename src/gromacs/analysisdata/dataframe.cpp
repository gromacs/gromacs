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
 * Implements classes in dataframe.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#include "dataframe.h"

#include "gromacs/fatalerror/gmxassert.h"

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
        int index, real x, real dx, int firstColumn, int columnCount,
        const real *y, const real *dy, const bool *present)
    : header_(index, x, dx), firstColumn_(firstColumn), columnCount_(columnCount),
      y_(y), dy_(dy), present_(present)
{
    GMX_ASSERT(firstColumn >= 0, "Invalid first column");
    GMX_ASSERT(columnCount >= 0, "Invalid column count");
    GMX_ASSERT(columnCount == 0 || y_ != NULL,
               "Values must be provided if there are columns");
}


AnalysisDataPointSetRef::AnalysisDataPointSetRef(
        const AnalysisDataFrameHeader &header, int firstColumn, int columnCount,
        const real *y, const real *dy, const bool *present)
    : header_(header), firstColumn_(firstColumn), columnCount_(columnCount),
      y_(y), dy_(dy), present_(present)
{
    GMX_ASSERT(firstColumn >= 0, "Invalid first column");
    GMX_ASSERT(columnCount >= 0, "Invalid column count");
    GMX_ASSERT(columnCount == 0 || y_ != NULL,
               "Values must be provided if there are columns");
}


AnalysisDataPointSetRef::AnalysisDataPointSetRef(
        const AnalysisDataPointSetRef &points, int firstColumn, int columnCount)
    : header_(points.header()), firstColumn_(0), columnCount_(columnCount),
      y_(points.y_), dy_(points.dy_), present_(points.present_)
{
    GMX_ASSERT(firstColumn >= 0, "Invalid first column");
    GMX_ASSERT(columnCount >= 0, "Invalid column count");
    if (points.lastColumn() < firstColumn
        || points.firstColumn() >= firstColumn + columnCount
        || columnCount == 0)
    {
        columnCount_ = 0;
        return;
    }
    int newFirstColumn = firstColumn - points.firstColumn();
    if (newFirstColumn > 0)
    {
        // Offset pointers if the first column is not the first in points.
        y_ += newFirstColumn;
        if (dy_ != NULL)
        {
            dy_ += newFirstColumn;
        }
        if (present_ != NULL)
        {
            present_ += newFirstColumn;
        }
        newFirstColumn = 0;
    }
    else
    {
        // Take into account if first column is before the first in points.
        columnCount_ -= -newFirstColumn;
    }
    // Decrease column count if there are not enough columns in points.
    if (newFirstColumn + columnCount_ > points.columnCount())
    {
        columnCount_ = points.columnCount() - newFirstColumn;
    }
}


bool AnalysisDataPointSetRef::allPresent() const
{
    if (present_ != NULL)
    {
        for (int i = 0; i < columnCount(); ++i)
        {
            if (!present_[i])
            {
                return false;
            }
        }
    }
    return true;
}


/********************************************************************
 * AnalysisDataFrameRef
 */

AnalysisDataFrameRef::AnalysisDataFrameRef()
    : points_(AnalysisDataFrameHeader(), 0, 0, NULL, NULL, NULL)
{
}


AnalysisDataFrameRef::AnalysisDataFrameRef(
        int index, real x, real dx, int columnCount,
        const real *y, const real *dy, const bool *present)
    : points_(index, x, dx, 0, columnCount, y, dy, present)
{
}


AnalysisDataFrameRef::AnalysisDataFrameRef(
        const AnalysisDataFrameHeader &header, int columnCount,
        const real *y, const real *dy, const bool *present)
    : points_(header, 0, columnCount, y, dy, present)
{
}


AnalysisDataFrameRef::AnalysisDataFrameRef(
        const AnalysisDataFrameRef &frame, int firstColumn, int columnCount)
    : points_(frame.points(), firstColumn, columnCount)
{
}

} // namespace gmx
