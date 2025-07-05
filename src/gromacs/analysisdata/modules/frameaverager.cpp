/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2013- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements gmx::AnalysisDataFrameAverager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#include "gmxpre.h"

#include "frameaverager.h"

#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

void AnalysisDataFrameAverager::setColumnCount(size_t columnCount)
{
    GMX_RELEASE_ASSERT(values_.empty(), "Cannot initialize multiple times");
    values_.resize(columnCount);
}

void AnalysisDataFrameAverager::addValue(size_t index, real value)
{
    AverageItem& item  = values_[index];
    const double delta = value - item.average;
    item.samples += 1;
    item.average += delta / item.samples;
    item.squaredSum += delta * (value - item.average);
}

void AnalysisDataFrameAverager::addPoints(const AnalysisDataPointSetRef& points)
{
    const int firstColumn = points.firstColumn();
    GMX_ASSERT(static_cast<size_t>(firstColumn + points.columnCount()) <= values_.size(),
               "Initialized with too few columns");
    for (int i = 0; i < points.columnCount(); ++i)
    {
        if (points.present(i))
        {
            addValue(firstColumn + i, points.y(i));
        }
    }
}

void AnalysisDataFrameAverager::finish()
{
    bFinished_ = true;
}

} // namespace gmx
