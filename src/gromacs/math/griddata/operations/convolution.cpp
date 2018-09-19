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
#include "convolution.h"

#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/griddata/operations/densitypadding.h"
#include "gromacs/math/griddata/operations/fouriertransform.h"
#include "gromacs/utility/real.h"

namespace gmx
{

GaussConvolution::GaussConvolution(const GridDataFloat3D &input)
    : input_(input), padded_input_
          (nullptr)
{}
GaussConvolution &GaussConvolution::pad(const MdFloatVector<DIM> &paddingFactor)
{
    padded_input_ = padDensity(input_, paddingFactor);
    return *this;
}

std::unique_ptr < GridDataFloat3D> GaussConvolution::convolute(real sigma)
{
    if (padded_input_ != nullptr)
    {
        fourierTransform_ =
            FourierTransformRealToComplex3D().normalize().result(*padded_input_);
    }
    else
    {
        fourierTransform_ =
            FourierTransformRealToComplex3D().normalize().result(input_);
    }
    auto sigmaSquared          = gmx::square(sigma);
    auto convoluteWithGaussian = [sigmaSquared](t_complex &value, RVec k) {
            auto prefactor = 1/(sqrt(2)) * exp(-2.0 * M_PI * M_PI * sigmaSquared * norm2(k));
            value.re *= prefactor;
            value.im *= prefactor;
        };

    ApplyToUnshiftedFourierTransform(*fourierTransform_).apply(convoluteWithGaussian);
    auto result =
        FourierTransformComplexToReal3D().normalize().result(*fourierTransform_);

    if (padded_input_ != nullptr)
    {
        result = unpadDensity(*result, input_.getGrid().lattice());
    }
    return result;
};
}
