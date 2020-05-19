/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief Defines the box deformation code.
 *
 * \todo The .mdp specification should have a boolean for this module.
 *
 * \todo grompp should set up fields in the tpr file that carry the
 * information about the original box, then the deform module would
 * can be build alongside the update module, rather than need to be
 * set up before the checkpoint is read.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "boxdeformation.h"

#include <memory>

#include "gromacs/gmxlib/network.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{

std::unique_ptr<BoxDeformation> prepareBoxDeformation(const matrix&     initialBox,
                                                      DDRole            ddRole,
                                                      NumRanks          numRanks,
                                                      MPI_Comm          communicator,
                                                      const t_inputrec& inputrec)
{
    if (!inputrecDeform(&inputrec))
    {
        return nullptr;
    }
    if (!EI_DYNAMICS(inputrec.eI))
    {
        GMX_THROW(NotImplementedError(
                "Box deformation is only supported with dynamical integrators"));
    }

    matrix box;
    // Only the rank that read the tpr has the global state, and thus
    // the initial box, so we pass that around.
    // (numRanks != NumRanks::Multiple helps clang static analyzer to
    // understand that box is defined in all cases)
    if (ddRole == DDRole::Master || numRanks != NumRanks::Multiple)
    {
        copy_mat(initialBox, box);
    }
    if (numRanks == NumRanks::Multiple)
    {
        gmx_bcast(sizeof(box), box, communicator);
    }

    return std::make_unique<BoxDeformation>(inputrec.delta_t, inputrec.init_step, inputrec.deform, box);
}

BoxDeformation::BoxDeformation(double        timeStep,
                               int64_t       initialStep,
                               const tensor& deformationTensor,
                               const matrix& referenceBox) :
    timeStep_(timeStep),
    initialStep_(initialStep)
{
    copy_mat(deformationTensor, deformationTensor_);
    copy_mat(referenceBox, referenceBox_);
}

void BoxDeformation::apply(ArrayRef<RVec> x, matrix box, int64_t step)
{
    matrix updatedBox, invbox, mu;

    double elapsedTime = (step + 1 - initialStep_) * timeStep_;
    copy_mat(box, updatedBox);
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            if (deformationTensor_[i][j] != 0)
            {
                updatedBox[i][j] = referenceBox_[i][j] + elapsedTime * deformationTensor_[i][j];
            }
        }
    }
    /* We correct the off-diagonal elements,
     * which can grow indefinitely during shearing,
     * so the shifts do not get messed up.
     */
    for (int i = 1; i < DIM; i++)
    {
        for (int j = i - 1; j >= 0; j--)
        {
            while (updatedBox[i][j] - box[i][j] > 0.5 * updatedBox[j][j])
            {
                rvec_dec(updatedBox[i], updatedBox[j]);
            }
            while (updatedBox[i][j] - box[i][j] < -0.5 * updatedBox[j][j])
            {
                rvec_inc(updatedBox[i], updatedBox[j]);
            }
        }
    }
    invertBoxMatrix(box, invbox);
    // Return the updated box
    copy_mat(updatedBox, box);
    mmul_ur0(box, invbox, mu);

    for (auto& thisX : x)
    {
        thisX[XX] = mu[XX][XX] * thisX[XX] + mu[YY][XX] * thisX[YY] + mu[ZZ][XX] * thisX[ZZ];
        thisX[YY] = mu[YY][YY] * thisX[YY] + mu[ZZ][YY] * thisX[ZZ];
        thisX[ZZ] = mu[ZZ][ZZ] * thisX[ZZ];
    }
}

} // namespace gmx
