/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
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
 * Implements gmx::AnalysisDataFrameAverager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#include "frameaverager.h"

#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

void AnalysisDataFrameAverager::addPoints(const AnalysisDataPointSetRef &points)
{
    GMX_ASSERT(!bFinished_, "Cannot add data after finished");
    const int firstColumn = points.firstColumn();
    GMX_ASSERT(static_cast<size_t>(firstColumn + points.columnCount()) <= values_.size(),
               "Initialized with too few columns");
    for (int i = 0; i < points.columnCount(); ++i)
    {
        if (points.present(i))
        {
            const real y                       = points.y(i);
            values_[firstColumn + i].average  += y;
            values_[firstColumn + i].variance += y * y;
            values_[firstColumn + i].samples  += 1;
        }
    }

}

void AnalysisDataFrameAverager::finish()
{
    if (bFinished_)
    {
        return;
    }
    bFinished_ = true;
    for (size_t i = 0; i < values_.size(); ++i)
    {
        if (values_[i].samples > 0)
        {
            const real average  = values_[i].average / values_[i].samples;
            const real variance =
                values_[i].variance / values_[i].samples - average * average;
            values_[i].average  = average;
            values_[i].variance = (variance > 0.0 ? variance : 0.0);
        }
    }

}

} // namespace gmx
