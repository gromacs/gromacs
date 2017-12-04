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
/*!  \file
 * \brief
 * Defines volume data containers.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 */
#ifndef GMX_MATH_FOURIERTRANSFORM_H
#define GMX_MATH_FOURIERTRANSFORM_H

#include "gromacs/math/griddata/griddata.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/real.h"
#include <fftw3.h>
#include <memory>
#include <vector>
#include <functional>

namespace gmx
{

std::array<int, 3> realGridExtendFromFourierTransfrom(const std::array<int, 3> &extend);
std::array<int, 3> fourierTransformGridExtendfromRealExtend(const std::array<int, 3> &extend);

/*! \brief
 * Convert Lattice to its corresponding lattice in reciprocal space.
 */
std::unique_ptr < IGrid < DIM>> convertGridToReciprocalSpace(const IGrid<DIM> &grid );

class FourierTransform3D
{
    public:
        FourierTransform3D(const FourierTransform3D &other) = default;
        ~FourierTransform3D();
        std::array<int, 3> columnMajorExtendToRowMajorExtend(const std::array<int, 3> &extend) const;
        void execute();

    protected:
        FourierTransform3D()  = default;
        fftwf_plan_s * plan_      = nullptr;
        bool           bNormalize = false;
};

class FourierTransformRealToComplex3D : public FourierTransform3D
{
    public:
        FourierTransformRealToComplex3D(const GridDataReal3D &realInputGridData);
        std::unique_ptr < GridDataComplex3D> result();
        void  result(GridDataComplex3D &complexTransformedGridData);
        FourierTransformRealToComplex3D &normalize();

    private:
        std::unique_ptr < GridDataComplex3D> createComplexTransformedGridDataFromInput_() const;
        const GridDataReal3D &realInputGridData_;
};

class FourierTransformComplexToReal3D : public FourierTransform3D
{
    public:
        FourierTransformComplexToReal3D(const GridDataComplex3D &complexInputGridData);
        std::unique_ptr < GridDataReal3D> result();
        void result(GridDataReal3D &realTransformedGridData);
        FourierTransformComplexToReal3D &normalize();

    private:
        std::unique_ptr < GridDataReal3D>  createRealTransformedGridDataFromInput_() const;
        const GridDataComplex3D &complexInputGridData_;
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
            std::vector<t_complex>::iterator itColumnStart,
            std::vector<t_complex>::iterator itColumnEnd,
            FunctionOnComplexGridData appliedFunction);

        void applyToAllRowsWithinColumn_(RVec coordinateColumnBegin, RVec deltakRow,
                                         std::vector<t_complex>::iterator itRowStart,
                                         std::vector<t_complex>::iterator itRowEnd,
                                         FunctionOnComplexGridData appliedFunction);
        GridDataComplex3D &field_;
};
}
#endif /* end of include guard: GMX_MATH_FOURIERTRANSFORM_H */
