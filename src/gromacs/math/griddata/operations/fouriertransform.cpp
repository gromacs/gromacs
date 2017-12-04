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
 * Fourier transform wrapper implementation.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 */

#include "fouriertransform.h"
#include "densitypadding.h"
#include "gromacs/math/griddata/griddata.h"
#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/fileio/griddataio.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/real.h"
#include <algorithm>
#include <cassert>

namespace gmx
{

std::array<int, 3> FourierTransform3D::columnMajorExtendToRowMajorExtend(const std::array<int, 3> &extend) const
{
    return { {
                 extend[ZZ], extend[YY], extend[XX]
             } };
}
std::unique_ptr < IGrid < DIM>> convertGridToReciprocalSpace(const IGrid<DIM> &grid )
{
    const auto         &extend = grid.lattice().extend();
    Grid<DIM>::NdVector scale;
    std::copy(std::begin(extend), std::end(extend), std::begin(scale));
    return std::unique_ptr < IGrid < DIM>>(new Grid<DIM>(grid.cell().inverse().scaledCopy(scale), grid.lattice()));
}


FourierTransformRealToComplex3D &
FourierTransformRealToComplex3D::normalize()
{
    bNormalize = true;
    return *this;
};

FourierTransformComplexToReal3D &
FourierTransformComplexToReal3D::normalize()
{
    bNormalize = true;
    return *this;
};

FourierTransform3D::~FourierTransform3D() { gmx_parallel_3dfft_destroy(plan_); }

std::array<int, 3> fourierTransformGridExtendfromRealExtend(const std::array<int, 3> &extend)
{
    return {
               {
                   extend[XX] / 2 + 1, extend[YY], extend[ZZ]
               }
    };
};

std::array<int, 3> realGridExtendFromFourierTransfrom(const std::array<int, 3> &extend)
{
    return {
               {
                   (extend[XX] - 1) * 2, extend[YY], extend[ZZ]
               }
    };
};

FourierTransformRealToComplex3D::FourierTransformRealToComplex3D( )
{ };
std::unique_ptr < GridDataComplex3D> FourierTransformRealToComplex3D::result(const GridDataReal3D &realInputGridData)
{
    /*
     * set to the size of the real grid first for correct reciprocal basis vectors,
     * then convert to abriged FFT size
     */
    auto complexTransformedGrid = convertGridToReciprocalSpace(realInputGridData.getGrid());
    complexTransformedGrid->setLatticeAndRescaleCell(
            fourierTransformGridExtendfromRealExtend(realInputGridData.getGrid().lattice().extend()));
    auto complexTransformedGridData = std::unique_ptr < GridDataComplex3D>(new GridDataComplex3D {std::move(complexTransformedGrid)});
    result(realInputGridData, *complexTransformedGridData);
    return complexTransformedGridData;
}

void FourierTransformRealToComplex3D::result(const GridDataReal3D &realInputGridData,
                                             GridDataComplex3D    &complexTransformedGridData)
{

    if (plan_ == nullptr)
    {
        // auto x = columnMajorExtendToRowMajorExtend(realInputGridData.getGrid().lattice().extend()).data();

        // (real **)realInputGridData.data()
        //(t_complex **)complexTransformedGridData.data()
        MPI_Comm comm[2] = {MPI_COMM_NULL, MPI_COMM_NULL};
        gmx_parallel_3dfft_init(&plan_,
                                columnMajorExtendToRowMajorExtend(realInputGridData.getGrid().lattice().extend()).data(),
                                &realData_, &complexData_, comm, false, 1);
        // gmx_parallel_3dfft_real_limits(plan_, );
    }
    gmx_parallel_3dfft_execute(plan_, GMX_FFT_REAL_TO_COMPLEX, 0, nullptr);
    if (bNormalize)
    {
        auto orthoscale = 1.0 / sqrt(complexTransformedGridData.size());

        auto normalizeComplexTransform = [orthoscale](t_complex &value) {
                value.re *= orthoscale;
                value.im *= orthoscale;
            };

        std::for_each(complexTransformedGridData.begin(),
                      complexTransformedGridData.end(),
                      normalizeComplexTransform);
    }
}

FourierTransformComplexToReal3D::FourierTransformComplexToReal3D()
{};

std::unique_ptr < GridDataReal3D>
FourierTransformComplexToReal3D::result(const GridDataComplex3D &complexInputGridData)
{
    Grid<DIM> realGrid {
        complexInputGridData.getGrid().cell(), complexInputGridData.getGrid().lattice()
    };
    realGrid.setLatticeAndRescaleCell(realGridExtendFromFourierTransfrom(complexInputGridData.getGrid().lattice().extend()));
    auto realGridPointer         = convertGridToReciprocalSpace(realGrid);
    auto realTransformedGridData = std::unique_ptr < GridDataReal3D>(new GridDataReal3D(std::move(realGridPointer)));

    result(complexInputGridData, *realTransformedGridData);
    return realTransformedGridData;
}

void
FourierTransformComplexToReal3D::result(const GridDataComplex3D  & /*complexInputGridData*/,
                                        GridDataReal3D          &realTransformedGridData)
{
    if (plan_ == nullptr)
    {
        // plan_ = fftwf_plan_dft_c2r(
        //             3, columnMajorExtendToRowMajorExtend(realTransformedGridData.getGrid().lattice().extend()).data(),
        //             (fftwf_complex *)complexInputGridData.data(),
        //             (float *)realTransformedGridData.data(), FFTW_ESTIMATE);
    }
    // fftwf_execute_dft_c2r(
    //         plan_, (fftwf_complex *)complexInputGridData.data(),
    //         (float *)realTransformedGridData.data());

    if (bNormalize)
    {
        auto orthoscale         = 1 / sqrt(realTransformedGridData.size());
        auto normalizeTransform = [orthoscale](real &value) {
                value *= orthoscale;
            };
        std::for_each(realTransformedGridData.begin(),
                      realTransformedGridData.end(), normalizeTransform);
    }
};

ApplyToUnshiftedFourierTransform::ApplyToUnshiftedFourierTransform(
        GridDataComplex3D &field)
    : field_ (field)
{};

const ApplyToUnshiftedFourierTransform &ApplyToUnshiftedFourierTransform::apply(
        FunctionOnComplexGridData appliedFunction)
{

    /*
     * Assume that second halves of indices denote negative frequencies
     * Use iterators to step though the grid.
     * This relies on the assumption that the grid is stored with z the slowest
     * and x the fasted changing dimension with no padding
     * (x,y,z not being linked to any coordiante system, but short-hand for
     * first, second, third dimension)
     * Loosing generality through this approach, we save substantial time when
     * we don't have to calculate the grid index.
     */

    const auto &fieldGrid    = field_.getGrid();
    auto        extend       = fieldGrid.lattice().extend();

    auto        k             = RVec {
        0., 0., 0.
    };
    RVec deltakRow     = {fieldGrid.unitCell().basisVectorLength(XX), 0, 0};
    RVec deltakColumn  = {0, fieldGrid.unitCell().basisVectorLength(YY), 0};
    RVec deltakSection = {0, 0, fieldGrid.unitCell().basisVectorLength(ZZ)};
    auto itMidSection  = field_.iteratorAtMultiIndex({{0, 0, extend[ZZ] / 2 + 1}});
    for (auto itValue = field_.iteratorAtMultiIndex({{0, 0, 0}}); itValue != itMidSection;
         itValue += extend[XX]*extend[YY])
    {
        applyToAllColumnsWithinSection_(
                k, extend[XX], extend[YY], deltakRow, deltakColumn,
                itValue, itValue + extend[XX]*extend[YY], appliedFunction);
        rvec_inc(k, deltakSection);
    }

    k = {0., 0., 0.};
    rvec_dec(k, deltakSection);

    for (auto itValue = field_.end(); itValue != itMidSection;
         itValue -= extend[XX]*extend[YY])
    {
        applyToAllColumnsWithinSection_(
                k, extend[XX], extend[YY], deltakRow, deltakColumn,
                itValue-extend[XX]*extend[YY], itValue, appliedFunction);
        rvec_dec(k, deltakSection);
    }
    return *this;
};

void ApplyToUnshiftedFourierTransform::applyToAllColumnsWithinSection_(
        RVec coordinateSectionBegin, int nRowsPerColumn, int nColumns,
        RVec deltakRow, RVec deltakColumn,
        GridDataComplex3D::iterator itColumnStart,
        GridDataComplex3D::iterator itColumnEnd,
        FunctionOnComplexGridData appliedFunction)
{

    auto itMiddleColumn = itColumnStart + ((nColumns + 1) / 2) * nRowsPerColumn;
    auto k              = coordinateSectionBegin;
    for (auto itValue = itColumnStart; itValue != itMiddleColumn;
         itValue += nRowsPerColumn)
    {
        applyToAllRowsWithinColumn_(k, deltakRow, itValue, itValue + nRowsPerColumn,
                                    appliedFunction);
        rvec_inc(k, deltakColumn);
    }

    k = coordinateSectionBegin;
    rvec_dec(k, deltakColumn);
    for (auto itValue = itColumnEnd; itValue != itMiddleColumn;
         itValue -= nRowsPerColumn)
    {
        applyToAllRowsWithinColumn_(k, deltakRow, itValue - nRowsPerColumn, itValue,
                                    appliedFunction);
        rvec_dec(k, deltakColumn);
    }
};

void ApplyToUnshiftedFourierTransform::applyToAllRowsWithinColumn_(
        RVec coordinateColumnBegin, RVec deltakRow,
        GridDataComplex3D::iterator itRowStart,
        GridDataComplex3D::iterator itRowEnd,
        FunctionOnComplexGridData appliedFunction)
{

    auto k = coordinateColumnBegin;
    for (auto itValue = itRowStart; itValue != itRowEnd; ++itValue)
    {
        appliedFunction(*itValue, k);
        rvec_inc(k, deltakRow);
    }
}
}
