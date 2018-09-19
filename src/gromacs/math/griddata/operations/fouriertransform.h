/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
/*!  \file
 * \brief
 * Defines volume data containers.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 */
#ifndef GMX_MATH_FOURIERTRANSFORM_H
#define GMX_MATH_FOURIERTRANSFORM_H

#include <functional>
#include <memory>

#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/griddata/griddata.h"
#include "gromacs/utility/real.h"

struct gmx_parallel_3dfft;

namespace gmx
{
bounds<DIM> realGridExtendFromFourierTransfrom(const bounds<DIM> &extend);
bounds<DIM> fourierTransformGridExtendfromRealExtend(const bounds<DIM> &extend);

/*! \brief
 * Convert Lattice to its corresponding lattice in reciprocal space.
 */
Grid <DIM> convertGridToReciprocalSpace(const Grid < DIM, GridWithTranslation < DIM>> &grid );

class FourierTransform3D
{
    public:
        FourierTransform3D(const FourierTransform3D &other) = default;
        ~FourierTransform3D();
        std::array<int, DIM> columnMajorExtendToRowMajorExtend(const bounds<DIM> &extend) const;

    protected:
        FourierTransform3D()  = default;
        gmx_parallel_3dfft * plan_      = nullptr;
        ivec                 rsize;
        ivec                 csize;
        real               * realData_;
        ivec                 ndata;
        t_complex          * complexData_;
        bool                 bNormalize = false;
};

class FourierTransformRealToComplex3D : public FourierTransform3D
{
    public:
        FourierTransformRealToComplex3D();
        std::unique_ptr < GridDataComplex3D> result(const GridDataFloat3D &realInputGridData);
        void  result(const GridDataFloat3D &realInputGridData, GridDataComplex3D &complexTransformedGridData);
        FourierTransformRealToComplex3D &normalize();
};

class FourierTransformComplexToReal3D : public FourierTransform3D
{
    public:
        FourierTransformComplexToReal3D();
        std::unique_ptr < GridDataFloat3D> result(const GridDataComplex3D &complexInputGridData);
        void result(const GridDataComplex3D &complexInputGridData, GridDataFloat3D &realTransformedGridData);
        FourierTransformComplexToReal3D &normalize();
};

class ApplyToUnshiftedFourierTransform
{
    public:
        typedef const std::function<void(t_complex &, RVec)> &FunctionOnComplexGridData;
        ApplyToUnshiftedFourierTransform(GridDataComplex3D &field);
        const ApplyToUnshiftedFourierTransform &
        apply(FunctionOnComplexGridData appliedFunction);

    private:
        void applyToAllColumnsWithinSection_(
            RVec coordinateSectionBegin, int nRowsPerColumn, int nColumns,
            RVec deltakRow, RVec deltakColumn,
            GridDataComplex3D::iterator itColumnStart,
            GridDataComplex3D::iterator itColumnEnd,
            FunctionOnComplexGridData appliedFunction);

        void applyToAllRowsWithinColumn_(RVec coordinateColumnBegin, RVec deltakRow,
                                         GridDataComplex3D::iterator itRowStart,
                                         GridDataComplex3D::iterator itRowEnd,
                                         FunctionOnComplexGridData appliedFunction);
        GridDataComplex3D &field_;
};
}
#endif /* end of include guard: GMX_MATH_FOURIERTRANSFORM_H */
