/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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

/*! \libinternal \file
 *
 * \brief
 * Common interface that the implementations of the Fast Multipole Methods (FMM) can use.
 *
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 *
 * \inlibraryapi
 * \ingroup module_fmm
 */
#ifndef GMX_FMM_FMM_IMPL_H
#define GMX_FMM_FMM_IMPL_H

#include <memory>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

struct fmmInputParameters;
struct gmx_mtop_t;
struct t_commrec;
struct t_inputrec;

namespace gmx
{

/*! \brief
 *  Identification string that can be prepended to FMM-related output
 */
constexpr char fmmStr[] = "FMM:";


class FmmImpl
{
    public:
        virtual void calculateForcesAndEnergies(const t_commrec *cr,
                                                const matrix     box,
                                                const real      *q,
                                                const rvec      *x,
                                                ArrayRef<RVec>   force,
                                                double          *coulombEnergy) = 0;

        virtual ~FmmImpl() { };
};

/*! \brief
 * Creates an instance of a specific FMM implementation (e.g., FmSolvr, GPU-FMM, or ExaFMM)
 *
 * Which FMM is instantiated depends on how the FMM implementation class is defined in the source
 * file that includes this header. The underlying FMM implementation class should be derived as
 * "class SpecificFmmImpl final : public FmmImpl" from the FmmImpl interface defined here.
 */
std::unique_ptr<FmmImpl> createFmmImpl(const t_inputrec *ir, const fmmInputParameters *fmmParameters, const gmx_mtop_t *mtop);

}      // namespace gmx

#endif // GMX_FMM_FMM_IMPL_H
