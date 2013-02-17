/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012, by the GROMACS development team, led by
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
 * Implements gmx::AnalysisDataProxy.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#include "dataproxy.h"

#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

AnalysisDataProxy::AnalysisDataProxy(int firstColumn, int columnSpan,
                                     AbstractAnalysisData *data)
    : source_(*data), firstColumn_(firstColumn), columnSpan_(columnSpan)
{
    GMX_RELEASE_ASSERT(data, "Source data must not be NULL");
    GMX_RELEASE_ASSERT(firstColumn >= 0 && columnSpan > 0, "Invalid proxy column");
    setColumnCount(columnSpan);
    setMultipoint(source_.isMultipoint());
}


AnalysisDataFrameRef
AnalysisDataProxy::tryGetDataFrameInternal(int index) const
{
    AnalysisDataFrameRef frame = source_.tryGetDataFrame(index);
    if (!frame.isValid())
    {
        return AnalysisDataFrameRef();
    }
    return AnalysisDataFrameRef(frame, firstColumn_, columnSpan_);
}


bool
AnalysisDataProxy::requestStorageInternal(int nframes)
{
    return source_.requestStorage(nframes);
}


int
AnalysisDataProxy::flags() const
{
    return efAllowMultipoint | efAllowMulticolumn | efAllowMissing;
}


void
AnalysisDataProxy::dataStarted(AbstractAnalysisData *data)
{
    GMX_RELEASE_ASSERT(data == &source_, "Source data mismatch");
    GMX_RELEASE_ASSERT(firstColumn_ + columnSpan_ <= source_.columnCount(),
                       "Invalid column(s) specified");
    notifyDataStart();
}


void
AnalysisDataProxy::frameStarted(const AnalysisDataFrameHeader &frame)
{
    notifyFrameStart(frame);
}


void
AnalysisDataProxy::pointsAdded(const AnalysisDataPointSetRef &points)
{
    AnalysisDataPointSetRef columns(points, firstColumn_, columnSpan_);
    if (columns.columnCount() > 0)
    {
        notifyPointsAdd(columns);
    }
}


void
AnalysisDataProxy::frameFinished(const AnalysisDataFrameHeader &header)
{
    notifyFrameFinish(header);
}


void
AnalysisDataProxy::dataFinished()
{
    notifyDataFinish();
}

} // namespace gmx
