/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
#ifndef GMX_MATH_REALFIELDMEASURE_H
#define GMX_MATH_REALFIELDMEASURE_H

#include <vector>
#include "gromacs/math/griddata/griddata.h"
#include "gromacs/math/vec.h"

namespace gmx
{

class RealGridDataMeasure
{
    public:
        explicit RealGridDataMeasure(const GridDataReal3D &realfield);
        ~RealGridDataMeasure() = default;
        RVec center_of_mass() const;
        float entropy() const;
    private:
        const GridDataReal3D &realfield_;

};

class DataVectorMeasure
{
    public:
        explicit DataVectorMeasure(const std::vector<float> &vecData);
        /*! \brief The sum of all grid values. */
        float sum() const;
        /*! \brief The minimum grid data value. */
        float min() const;
        /*! \brief The maximum grid data value. */
        float max() const;
        /*! \brief The mean grid data value. */
        float mean() const;
        /*! \brief
         * The size of the griddata in bytes.
         */
        size_t data_size() const;

        float sumOfSquareDeviation() const;

        float normSquared() const;

        float norm() const;
        /*! \brief The variance of data values.
         *
         * Two pass algorithm, where mean is calculated fist, due to large expected
         * amount of data. (Typical data size=1,000,000)
         */
        float var() const;
        /*! \brief The root mean square deviation of grid data values. */
        float rms() const;

    private:
        const std::vector<float> &vecData_;
};

}      // namespace gmx
#endif
