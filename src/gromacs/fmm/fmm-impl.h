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

#include "gromacs/fmm/fmm.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/utility/arrayref.h"

struct t_commrec;

namespace gmx
{

struct ForceProviderInitOptions;
struct fmmInputParameters;
class FmmImpl;

/*! \brief
 * Creates an instance of a specific FMM implementation (e.g., FmSolvr, GPU-FMM, or ExaFMM)
 *
 * Any underlying FMM needs to implement the methods in fmm-impl.h. The FMM implementation
 * class should be derived from IForceProvider as "class FmmImpl final : public IForceProvider"
 */
std::unique_ptr<FmmImpl> createFmmImpl(const fmmInputParameters *fmmParameters, const ForceProviderInitOptions *options);


class FmmImpl : public IForceProvider
{
    public:
        FmmImpl(gmx_unused const fmmInputParameters       *fmmParameters,
                gmx_unused const ForceProviderInitOptions *options) { };

        void calculateForces(gmx_unused const t_commrec *cr,
                             gmx_unused const t_mdatoms *mdatoms,
                             gmx_unused const matrix     box,
                             gmx_unused double           t,
                             gmx_unused const rvec      *x,
                             gmx_unused ForceWithVirial *forceWithVirial,
                             gmx_unused gmx_enerdata_t  *enerd) override { };
        ~FmmImpl() { };
};

}      // namespace gmx

#endif // GMX_FMM_FMM_IMPL_H
