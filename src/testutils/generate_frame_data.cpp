/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
/*! \brief Test helper utilities for generating frame data.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "gmxpre.h"

#include "testutils/generate_frame_data.h"

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/matrix.h"
#include "gromacs/utility/vectypes.h"

namespace gmx
{
namespace test
{
namespace
{

template<typename T>
static void fillVectorFrame(ArrayRef<T[DIM]>                      values,
                            TrajectoryFrameDataGenerator::Offset* offset,
                            const TrajectoryFrameMode             mode)
{
    T sign = mode == TrajectoryFrameMode::Positive ? 1.0 : -1.0;
    for (int i = 0; i < gmx::ssize(values); ++i)
    {
        for (int d = 0; d < DIM; ++d)
        {
            values[i][d] = sign * (offset->current_ + (i * offset->value_) + ((d + 1) * offset->dim_));
        }

        if (mode == TrajectoryFrameMode::Alternating)
        {
            sign *= -1.0;
        }
    }
    offset->current_ += offset->frame_;
}

template<typename T>
static void fillMatrix(ArrayRef<T> values, MatrixFrameDataGenerator::Offset* offset)
{
    T value = offset->value_;
    for (int i = 0; i < gmx::ssize(values); ++i)
    {
        values[i] = offset->current_ + value;
        value += offset->value_;
    }
    offset->current_ += offset->frame_;
}

}; // namespace

TrajectoryFrameDataGenerator::TrajectoryFrameDataGenerator(const TrajectoryFrameMode mode,
                                                           const Offset&             offset) :
    mode_{ mode }, offset_{ offset }
{
}

void TrajectoryFrameDataGenerator::operator()(ArrayRef<BasicVector<float>> values)
{
    fillVectorFrame(arrayRefFromArray(as_vec_array(values.data()), values.size()), &offset_, mode_);
}

void TrajectoryFrameDataGenerator::operator()(ArrayRef<BasicVector<double>> values)
{
    fillVectorFrame(arrayRefFromArray(as_vec_array(values.data()), values.size()), &offset_, mode_);
}

void TrajectoryFrameDataGenerator::operator()(ArrayRef<float[3]> values)
{
    fillVectorFrame(values, &offset_, mode_);
}

void TrajectoryFrameDataGenerator::operator()(ArrayRef<double[3]> values)
{
    fillVectorFrame(values, &offset_, mode_);
}

TrajectoryFrameDataGenerator::Offset::Offset(const double start,
                                             const double frame,
                                             const double value,
                                             const double dim) :
    current_{ start }, frame_{ frame }, value_{ value }, dim_{ dim }
{
}

MatrixFrameDataGenerator::MatrixFrameDataGenerator(const Offset& offset) : offset_{ offset } {}

void MatrixFrameDataGenerator::operator()(float box[3][3])
{
    fillMatrix(arrayRefFromArray(box[0], DIM * DIM), &offset_);
}

void MatrixFrameDataGenerator::operator()(double box[3][3])
{
    fillMatrix(arrayRefFromArray(box[0], DIM * DIM), &offset_);
}

void MatrixFrameDataGenerator::operator()(BasicMatrix3x3<float>& box)
{
    fillMatrix(box.toArrayRef(), &offset_);
}

void MatrixFrameDataGenerator::operator()(BasicMatrix3x3<double>& box)
{
    fillMatrix(box.toArrayRef(), &offset_);
}

MatrixFrameDataGenerator::Offset::Offset(const double start, const double frame, const double value) :
    current_{ start }, frame_{ frame }, value_{ value }
{
}

} // namespace test
} // namespace gmx
