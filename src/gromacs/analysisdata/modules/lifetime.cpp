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
 * Implements gmx::AnalysisDataLifetimeModule.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#include "lifetime.h"

#include <cmath>

#include <deque>
#include <vector>

#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/datastorage.h"

namespace gmx
{

/********************************************************************
 * AnalysisDataLifetimeModule
 */

class AnalysisDataLifetimeModule::Impl
{
    public:
        typedef std::deque<int> LifetimeHistogram;

        Impl() : firstx_(0.0), lastx_(0.0), frameCount_(0), bCumulative_(false)
        {
        }

        void addLifetime(int dataSet, int lifetime)
        {
            if (lifetime > 0)
            {
                LifetimeHistogram &histogram = lifetimeHistograms_[dataSet];
                if (histogram.size() < static_cast<unsigned>(lifetime))
                {
                    histogram.resize(lifetime, 0);
                }
                ++histogram[lifetime - 1];
            }
        }

        real                            firstx_;
        real                            lastx_;
        int                             frameCount_;
        bool                            bCumulative_;
        std::vector<std::vector<int> >  currentLifetimes_;
        std::vector<LifetimeHistogram>  lifetimeHistograms_;
};

AnalysisDataLifetimeModule::AnalysisDataLifetimeModule()
    : impl_(new Impl())
{
}

AnalysisDataLifetimeModule::~AnalysisDataLifetimeModule()
{
}

void AnalysisDataLifetimeModule::setCumulative(bool bCumulative)
{
    impl_->bCumulative_ = bCumulative;
}

int AnalysisDataLifetimeModule::flags() const
{
    return efAllowMulticolumn | efAllowMissing | efAllowMultipleDataSets;
}

void
AnalysisDataLifetimeModule::dataStarted(AbstractAnalysisData *data)
{
    impl_->currentLifetimes_.reserve(data->dataSetCount());
    impl_->lifetimeHistograms_.reserve(data->dataSetCount());
    for (int i = 0; i < data->dataSetCount(); ++i)
    {
        impl_->currentLifetimes_.push_back(std::vector<int>(data->columnCount(i), 0));
        impl_->lifetimeHistograms_.push_back(std::deque<int>());
    }
}

void
AnalysisDataLifetimeModule::frameStarted(const AnalysisDataFrameHeader &header)
{
    if (header.index() == 0)
    {
        impl_->firstx_ = header.x();
    }
    impl_->lastx_ = header.x();
    ++impl_->frameCount_;
    // TODO: Check the input for even spacing.
}

void
AnalysisDataLifetimeModule::pointsAdded(const AnalysisDataPointSetRef &points)
{
    const int dataSet = points.dataSetIndex();
    // This assumption is strictly not necessary, but this is how the
    // framework works currently, and makes the code below simpler.
    GMX_ASSERT(points.firstColumn() == 0
               && points.lastColumn() == static_cast<int>(impl_->currentLifetimes_[dataSet].size()) - 1,
               "Point set should cover all columns");
    for (int i = 0; i < points.columnCount(); ++i)
    {
        // TODO: Perhaps add control over how this is determined?
        const bool bPresent = points.present(i) && points.y(i) > 0.0;
        if (bPresent)
        {
            ++impl_->currentLifetimes_[dataSet][i];
        }
        else if (impl_->currentLifetimes_[dataSet][i] > 0)
        {
            impl_->addLifetime(dataSet, impl_->currentLifetimes_[dataSet][i]);
            impl_->currentLifetimes_[dataSet][i] = 0;
        }
    }
}

void
AnalysisDataLifetimeModule::frameFinished(const AnalysisDataFrameHeader & /*header*/)
{
}

void
AnalysisDataLifetimeModule::dataFinished()
{
    // Count the elements present in the last frame.
    for (size_t i = 0; i < impl_->currentLifetimes_.size(); ++i)
    {
        for (size_t j = 0; j < impl_->currentLifetimes_[i].size(); ++j)
        {
            impl_->addLifetime(i, impl_->currentLifetimes_[i][j]);
        }
    }
    impl_->currentLifetimes_.clear();

    if (impl_->bCumulative_)
    {
        // Sum up all subintervals for longer intervals.
        std::vector<Impl::LifetimeHistogram>::iterator histogram;
        for (histogram  = impl_->lifetimeHistograms_.begin();
             histogram != impl_->lifetimeHistograms_.end();
             ++histogram)
        {
            Impl::LifetimeHistogram::iterator i, j;
            for (i = histogram->begin(); i != histogram->end(); ++i)
            {
                int n = 2; // Number of subintervals.
                for (j = i + 1; j != histogram->end(); ++j, ++n)
                {
                    *i += n * (*j);
                }
            }
        }
    }

    const real spacing =
        (impl_->frameCount_ > 1)
        ? (impl_->lastx_ - impl_->firstx_) / (impl_->frameCount_ - 1)
        : 0.0;
    setXAxis(0.0, spacing);

    setColumnCount(impl_->lifetimeHistograms_.size());
    std::vector<Impl::LifetimeHistogram>::const_iterator histogram;
    size_t maxLifetime = 1;
    for (histogram  = impl_->lifetimeHistograms_.begin();
         histogram != impl_->lifetimeHistograms_.end();
         ++histogram)
    {
        maxLifetime = std::max(maxLifetime, histogram->size());
    }
    setRowCount(maxLifetime);

    allocateValues();
    int column = 0;
    for (histogram  = impl_->lifetimeHistograms_.begin();
         histogram != impl_->lifetimeHistograms_.end();
         ++histogram, ++column)
    {
        int row = 0;
        Impl::LifetimeHistogram::const_iterator i;
        for (i = histogram->begin(); i != histogram->end(); ++i, ++row)
        {
            value(row, column).setValue(*i / static_cast<real>(impl_->frameCount_ - row));
        }
        for (; row < rowCount(); ++row)
        {
            value(row, column).setValue(0.0);
        }
    }
    impl_->lifetimeHistograms_.clear();
    valuesReady();
}

} // namespace gmx
