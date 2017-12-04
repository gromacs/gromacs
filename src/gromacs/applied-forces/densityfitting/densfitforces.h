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
#ifndef GMX_APPLIEDFORCES_DENSITYFITTING_DENSFITFORCES_H_
#define GMX_APPLIEDFORCES_DENSITYFITTING_DENSFITFORCES_H_
#include <memory>
#include "gromacs/math/vectypes.h"
#include "gromacs/math/griddata/griddata.h"
#include "gromacs/math/paddedvector.h"

namespace gmx
{
template <int N> class IGrid;

class DensfitForces
{
    public:
        DensfitForces(const IGrid<DIM> &grid, real sigma, real nSigma);
        RVec force(const RVec &x, const GridDataReal3D &densityDensityDerivative);
        void setSigma(real sigma, real nSigma);

    private:
        std::unique_ptr < IGrid < DIM>> grid_;
        int          voxrange_; /**< Max. number of voxels to be computed for a single atom in a single dimension x, y, or z   */
        const double dRhoDxPrefactor_ = sqrt(2. / M_PI) * 1 / (2. * 2. * 2.);
        real         nu_;
        /* The following two temporary vectors (one for each OpenMP thread) store
         * erf values around a single atoms, thus we can compute them all in one go,
         * and with SIMD acceleration */
        using AlignedRealVector = std::vector < real, gmx::AlignedAllocator < real>>;
        std::vector<AlignedRealVector>
        erfVector;                                /**< vector of vectors of erf values */
        std::vector<AlignedRealVector> expVector; /**< same for exp values */

        real sigma_;
};
}
#endif /* end of include guard: \
          GMX_APPLIEDFORCES_DENSITYFITTING_DENSFITFORCES_H_ */
