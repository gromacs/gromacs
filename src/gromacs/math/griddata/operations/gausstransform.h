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
#ifndef _GAUSSTRANSFORM_H
#define _GAUSSTRANSFORM_H

#include "gromacs/math/vec.h"

#include "gromacs/math/griddata/columnmajorlattice.h"
#include "gromacs/math/griddata/griddata.h"
#include <array>
#include <vector>
#include <memory>

namespace gmx
{

typedef BasicVector<int> IVec;
template <int N> class CanonicalVectorBasis;

class GaussTransform
{
    public:
        GaussTransform()  = default;
        ~GaussTransform() = default;
        virtual void set_grid(std::unique_ptr < GridDataReal3D> grid) = 0;
        void set_sigma(real sigma);
        void set_n_sigma(real n_sigma);
        virtual void transform(const rvec x, real weight) = 0;
        virtual std::unique_ptr < GridDataReal3D> finish_and_return_grid() = 0;
        IVec getMinimumUsedGridIndex();
        IVec getMaximumUsedGridIndex();
    protected:
        // no other object should have access to the grid while Gauss transform is in progress
        std::unique_ptr < GridDataReal3D> grid_;
        real                              sigma_;
        real                              n_sigma_;
        IVec minimumUsedGridIndex_;
        IVec maximumUsedGridIndex_;
};

/*! \brief Efficient spreading of sources on a grid with a Gaussian kernel.
 * For reference see:
 * Leslie Greengard and June-Yub Lee, Accelerating the Nonuniform Fast Fourier Transform
 * SIAM REV 2004 Vol. 46, No. 3, pp. 443–454
 * */
class FastGaussianGridding : public GaussTransform
{

    public:
        FastGaussianGridding()  = default;
        ~FastGaussianGridding() = default;
        /*! \brief Checks if grid is rectangular and equispaced.
         */
        void set_grid(std::unique_ptr < GridDataReal3D> grid);
        /*! \brief Perform gaussian spreading of one source with a weight.
         *
         * Feed one source at a time.
         */
        void transform(const rvec x, real weight);
        /*! \brief Perform any outstanding caluclations, then hand back ownership of the grid.
         */
        std::unique_ptr < GridDataReal3D> finish_and_return_grid();
    protected:
        void prepare_2d_grid(const rvec x, const real weight);
        ColumnMajorLattice<DIM>::MultiIndex grid_index_of_spread_atom_;
        int                                 m_spread_;
        std::array<std::vector<real>, 3>    spread_1d_;
        std::vector < std::vector < real>> spread_2d_;
    private:
        void spread_1d(real weight, int m_spread, rvec dx, real nu, const std::vector<real> &E3, int dimension);
        void tensor_product_2d_();
        void tensor_product_();
        bool allVectorsSameLength(const CanonicalVectorBasis<3> &basis, real ftol, real abstol) const;

        real nu_; // spacing/sigma
        std::vector < std::vector < std::vector<real>>> spread_block_;

        std::vector<real>                E3_;      //< exp(-l^2*nu^2/2) , following the naming convention of Greengard et al., ;
};

}      // namespace gmx

#endif /* end of include guard: _GAUSSTRANSFORM_H */
