/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
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

#include <array>
#include <memory>

#include "gromacs/gmxlib/network.h"
#include "gromacs/math/boxmatrix.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"

namespace gmx
{

std::unique_ptr<BoxDeformation> buildBoxDeformation(const Matrix3x3&  initialBox,
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

    // Do not allow reading old tpr files from versions where deform scaled the flow field,
    // as the velocities will not obey the flow field.
    // There tpr file was not updated with the deform change, but version < 130 separates releases.
    if (inputrec.tpxFileVersion < 130)
    {
        GMX_THROW(InvalidInputError(
                "The implementation of the deform functionality has changed between the GROMACS "
                "versions used to generate the tpr file and this binary. Please read the "
                "documentation and regenerate the tpr file with a newer version of GROMACS to set "
                "up the initial flow field."));
    }

    Matrix3x3 box;
    // Only the rank that read the tpr has the global state, and thus
    // the initial box, so we pass that around.
    // (numRanks != NumRanks::Multiple helps clang static analyzer to
    // understand that box is defined in all cases)
    if (ddRole == DDRole::Main || numRanks != NumRanks::Multiple)
    {
        box = initialBox;
    }
    if (numRanks == NumRanks::Multiple)
    {
        auto boxArrayRef = box.toArrayRef();
        gmx_bcast(boxArrayRef.size() * sizeof(boxArrayRef[0]), boxArrayRef.data(), communicator);
    }

    return std::make_unique<BoxDeformation>(
            inputrec.delta_t, inputrec.init_step, createMatrix3x3FromLegacyMatrix(inputrec.deform), box);
}

BoxDeformation::BoxDeformation(const double     timeStep,
                               const int64_t    initialStep,
                               const Matrix3x3& deformationTensor,
                               const Matrix3x3& referenceBox) :
    timeStep_(timeStep),
    initialStep_(initialStep),
    deformationTensor_(deformationTensor),
    referenceBox_(referenceBox)
{
}

void BoxDeformation::apply(Matrix3x3* box, const int64_t step)
{
    const real elapsedTime = (step + 1 - initialStep_) * timeStep_;
    Matrix3x3  updatedBox  = *box;
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            if (deformationTensor_(i, j) != 0)
            {
                updatedBox(i, j) = referenceBox_(i, j) + elapsedTime * deformationTensor_(i, j);
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
            while (updatedBox(i, j) - (*box)(i, j) > 0.5_real * updatedBox(j, j))
            {
                updatedBox(i, XX) -= updatedBox(j, XX);
                updatedBox(i, YY) -= updatedBox(j, YY);
                updatedBox(i, ZZ) -= updatedBox(j, ZZ);
            }
            while (updatedBox(i, j) - (*box)(i, j) < -0.5_real * updatedBox(j, j))
            {
                updatedBox(i, XX) += updatedBox(j, XX);
                updatedBox(i, YY) += updatedBox(j, YY);
                updatedBox(i, ZZ) += updatedBox(j, ZZ);
            }
        }
    }

    // Return the updated box
    *box = updatedBox;
}

void setBoxDeformationFlowMatrix(const matrix boxDeformationVelocity, const matrix box, matrix flowMatrix)
{
    for (int d1 = 0; d1 < DIM; d1++)
    {
        for (int d2 = 0; d2 < DIM; d2++)
        {
            // The flow matrix is transposed with respect to the deform matrix
            flowMatrix[d1][d2] = boxDeformationVelocity[d2][d1] / box[d2][d2];
        }
    }
}

} // namespace gmx
