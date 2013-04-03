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
 * Declares gmx::AnalysisDataFrameAverager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_MODULES_FRAMEAVERAGER_H
#define GMX_ANALYSISDATA_MODULES_FRAMEAVERAGER_H

#include <vector>

#include "../../legacyheaders/types/simple.h"

#include "../../utility/gmxassert.h"

namespace gmx
{

class AnalysisDataPointSetRef;

/*! \internal \brief
 * Helper class for modules that average values over frames.
 *
 * This class implements common functionality for analysis data modules that
 * need to average a set of values over frames.  Currently, it is designed for
 * computing averages for each input column independently, but should be
 * relatively easy to make more general if required.
 *
 * This class takes care of accumulating the values and computing their
 * variance.  It allows different number of samples for each input column.
 * Accumulation is always in double precision.  For now, does nothing fancy,
 * but provides ground for other implementation (e.g., related to
 * parallelization) that would benefit all such modules.
 *
 * Methods in this class do not throw unless otherwise indicated.
 *
 * \ingroup module_analysisdata
 */
class AnalysisDataFrameAverager
{
    public:
        AnalysisDataFrameAverager() : bFinished_(false) {}

        /*! \brief
         * Sets the number of columns in the input data.
         *
         * \throws std::bad_alloc if out of memory.
         *
         * Typically called from AnalysisDataModuleInterface::dataStarted().
         *
         * Must be called exactly once, before setting calling any other method
         * in the class.
         */
        void setColumnCount(int columnCount)
        {
            GMX_RELEASE_ASSERT(columnCount >= 0, "Invalid column count");
            GMX_RELEASE_ASSERT(values_.empty(),
                               "Cannot initialize multiple times");
            values_.resize(columnCount);
        }
        /*! \brief
         * Accumulates data from a given point set into the average.
         *
         * Typically called from AnalysisDataModuleInterface::pointsAdded().
         *
         * Each call accumulates the values for those columns that are present
         * in the point set.  Can be called multiple times for a frame, and
         * does not need to be called for every frame.
         */
        void addPoints(const AnalysisDataPointSetRef &points);
        /*! \brief
         * Finalizes the calculation of the averages and variances.
         *
         * Does any computation that is not done during the accumulation in
         * addPoints(), e.g., normalizing by the number of samples.
         *
         * Typically called from AnalysisDataModuleInterface::dataFinished().
         */
        void finish();

        /*! \brief
         * Returns the computed average for a given column.
         *
         * If called before finish(), the results are undefined.
         */
        real average(int index) const
        {
            GMX_ASSERT(index >= 0 && index <= static_cast<int>(values_.size()),
                       "Invalid column index");
            GMX_ASSERT(bFinished_,
                       "Values available only after finished() has been called");
            return values_[index].average;
        }
        /*! \brief
         * Returns the computed variance for a given column.
         *
         * If called before finish(), the results are undefined.
         */
        real variance(int index) const
        {
            GMX_ASSERT(index >= 0 && index <= static_cast<int>(values_.size()),
                       "Invalid column index");
            GMX_ASSERT(bFinished_,
                       "Values available only after finished() has been called");
            return values_[index].variance;
        }

    private:
        struct AverageItem
        {
            AverageItem() : average(0.0), variance(0.0), samples(0) {}

            double               average;
            double               variance;
            int                  samples;
        };

        std::vector<AverageItem> values_;
        bool                     bFinished_;
};

} // namespace gmx

#endif
