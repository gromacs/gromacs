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
 * Declares the NNPot Force Provider class
 *
 * \author Lukas Müllender <lukas.muellender@gmail.com>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_NNPOTFORCEPROVIDER_H
#define GMX_APPLIED_FORCES_NNPOTFORCEPROVIDER_H

#include <string.h>

#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/topology/atoms.h"

namespace gmx
{

struct NNPotParameters;
class INNPotModel;

class MDLogger;

/*! \brief \internal
 * NNPotForceProvider class
 *
 * Implements the IForceProvider interface for the NNPot force provider.
 */
class NNPotForceProvider final : public IForceProvider
{
public:
    NNPotForceProvider(const NNPotParameters&, const MDLogger* logger);

    //! Destroy force provider for NNPot
    ~NNPotForceProvider();

    /*! \brief Calculate forces of NNPot.
     *
     * Prepares the input for the neural network model triggers model inference.
     * \param[in] fInput input for force provider
     * \param[out] fOutput output for force provider
     */
    void calculateForces(const ForceProviderInput& fInput, ForceProviderOutput* fOutput) override;

    //! Gather atom numbers and indices. Triggered on AtomsRedistributed signal.
    void gatherAtomNumbersIndices();

    //! Gather atom positions for NN input.
    void gatherAtomPositions(ArrayRef<const RVec> pos);

private:
    //! reference to NNPot parameters
    const NNPotParameters& params_;

    //! neural network model
    std::shared_ptr<INNPotModel> model_;

    //! vector storing all atom positions
    std::vector<RVec> positions_;

    //! vector storing all atomic numbers
    std::vector<int> atomNumbers_;

    //! global index lookup table to map indices from model input to global atom indices
    std::vector<int> idxLookup_;

    //! local copy of simulation box
    matrix box_;

    //! logger
    const MDLogger* logger_;

    //! stores communication record
    const t_commrec* cr_;
};

} // namespace gmx

#endif
