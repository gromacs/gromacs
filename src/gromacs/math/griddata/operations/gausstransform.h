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
#ifndef _GAUSSTRANSFORM_H
#define _GAUSSTRANSFORM_H

#include <array>
#include <memory>
#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/math/griddata/griddata.h"
#include "gromacs/math/griddata/mdarrayindexing.h"
 #include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

template <size_t N> class CanonicalVectorBasis;

class ISpreader
{
    public:
        virtual ArrayRef<real> spread(real weight, real x) = 0;
};


class GaussTransform
{
    public:
        GaussTransform() {};
        void setSigma(real sigma);
        void setGrid(GridDataFloat3D * grid, bool gaussianIsIntegrated = false);
        void setSpreadWidthInNSigma(real nSigma);
        /*! \brief Perform gaussian spreading of one source with a weight.
         *
         * Feed one source at a time.
         */
        void addTransform(const rvec x, real weight);
        offset<DIM> minimumUsedGridIndex();
        offset<DIM> maximumUsedGridIndex();
    private:
        void prepare_2d_grid(const rvec x, const real weight);
        offset<DIM> atomGridIndex_;
        int                                          m_spread_ = 0;
        std::array<std::unique_ptr<ISpreader>, 3>    OneDSpreader_;
        std::array<ArrayRef<real>, 3>                spread_1d_;
        std::vector < std::vector < real, gmx::AlignedAllocator < real>>> spread_2d_;
        void tensor_product_2d_();
        void tensor_product_();
        bool allVectorsSameLength(const CanonicalVectorBasis<3> &basis, real ftol, real abstol) const;

        real nu_ = 0;         // spacing/sigma
        std::vector < std::vector < std::vector < real, gmx::AlignedAllocator < real>> >> spread_block_;
        std::vector < std::vector < int>> ceilSqrtLUT_;


        GridDataFloat3D * grid_;
        real              sigma_  = 0;
        real              nSigma_ = 0;
        offset<DIM>       minimumUsedGridIndex_;
        offset<DIM>       maximumUsedGridIndex_;
};


class GaussTransformOneD : public ISpreader
{
    public:
        GaussTransformOneD(int m_spread, real nu, real sigma);
        ArrayRef<real> spread(real weight, real dx) override;
    private:
        int               m_spread_;
        real              nu_;
        real              sigma_;

        std::vector < real, gmx::AlignedAllocator < real>> spread_;
        std::vector < real, gmx::AlignedAllocator < real>> E3_;  //< exp(-l^2*nu^2/2) , following the naming convention of Greengard et al., ;

};

class GaussIntergralOneD : public ISpreader
{
    public:
        GaussIntergralOneD(int m_spread, real nu);
        ArrayRef<real> spread(real weight, real dx) override;
    private:
        int  m_spread_;

        std::vector < real, gmx::AlignedAllocator < real>> erfBuffer_;
        std::vector < real, gmx::AlignedAllocator < real>> voxelBoundaryBuffer_;

};

/*! \brief Efficient spreading of sources on a grid with a Gaussian kernel.
 * For reference see:
 * Leslie Greengard and June-Yub Lee, Accelerating the Nonuniform Fast Fourier Transform
 * SIAM REV 2004 Vol. 46, No. 3, pp. 443–454
 * */

}      // namespace gmx

#endif /* end of include guard: _GAUSSTRANSFORM_H */
