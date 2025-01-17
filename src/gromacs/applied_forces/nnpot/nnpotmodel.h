/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * \brief
 * Declares an abstract wrapper class for a neural network model to predict energies/forces.
 *
 * \author Lukas Müllender <lukas.muellender@gmail.com>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_NNPOTMODEL_H
#define GMX_APPLIED_FORCES_NNPOTMODEL_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"

struct t_commrec;
struct gmx_enerdata_t;
enum class PbcType;

namespace gmx
{

/*! \brief NNPot Module
 *
 * Abstract class for a neural network potential model.
 * Inherit from this class to implement a specific neural network backend/framework.
 * Refer to torchmodel.h for an example implementation.
 */
class INNPotModel
{
public:
    //! Initialize the neural network model
    virtual void initModel() = 0;

    /*! \brief Prepare inputs for NN model.
     *
     * Currently supported inputs:
     *   - atom positions (std::vector<RVec>): atomic positions
     *   - atom numbers (std::vector<int>): atomic numbers
     *   - box (matrix): simulation box vectors
     *   - pbc type (PbcType): boolean flags for periodic boundary conditions in x, y, z
     */
    //! \{
    virtual void prepareAtomPositions(std::vector<RVec>&) = 0;
    virtual void prepareAtomNumbers(std::vector<int>&)    = 0;
    virtual void prepareBox(matrix&)                      = 0;
    virtual void preparePbcType(PbcType&)                 = 0;
    //! \}

    //! call inference on NN model
    virtual void evaluateModel() = 0;

    //! retrieve NN model outputs
    virtual void getOutputs(std::vector<int>&, gmx_enerdata_t&, const ArrayRef<RVec>&) = 0;

    //! set communication record for possible communication of input/output data between ranks
    virtual void setCommRec(const t_commrec*) = 0;

    //! helper function to check if model outputs forces
    virtual bool outputsForces() const = 0;
};

} // namespace gmx

#endif // GMX_APPLIED_FORCES_NNPOTMODEL_H
