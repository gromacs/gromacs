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
 * Declares a wrapper class for a TorchScript-compiled neural network model.
 *
 * \author Lukas Müllender <lukas.muellender@gmail.com>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_TORCHMODEL_H
#define GMX_APPLIED_FORCES_TORCHMODEL_H

// undef needed because of macro clashes with libtorch/ATen
#ifdef DIM
#    undef DIM
#endif
#include <torch/cuda.h>
#include <torch/script.h>
#include <torch/torch.h>
#ifdef DIM
#    undef DIM
#endif

#include <string>

#include "gromacs/applied_forces/nnpot/nnpotmodel.h"

namespace gmx
{

enum class NNPotEmbedding;

/*! \brief Define the torch datatype according to GMX_DOUBLE.
 *
 * Important for converting data types, as model inference is always done in float32.
 */
static constexpr auto torchRealType = GMX_DOUBLE ? torch::kFloat64 : torch::kFloat32;

class MpiComm;
class MDLogger;
class LinkFrontierAtom;

/*! \brief
 * Class responsible for loading and evaluating a TorchScript-compiled neural network model.
 * Inherits from NNPotModel.
 */
class TorchModel : public INNPotModel
{
public:
    /*! \brief Constructor for TorchModel.
     * \param[in] filename path to the TorchScript model file
     * \param[in] embedding embedding scheme used for NNP/MM interaction
     * \param[in] logger handle to MDLogger
     * \param[in] mpiComm handle to MpiComm
     */
    TorchModel(const std::string& filename, NNPotEmbedding embedding, const MDLogger& logger, const MpiComm& mpiComm);

    /*! Call inference on NN model and retrieve outputs
     * \param[out] enerd energy data struct
     * \param[out] forces forces on atoms
     * \param[in] indexLookup lookup table for atom indices
     * \param[in] mmIndices indices of MM atoms
     * \param[in] inputs list of strings specifying input data
     * \param[in] positions atom positions
     * \param[in] atomNumbers atom numbers
     * \param[in] atomPairs list of all input atom pairs within cutoff
     * \param[in] pairShifts list of periodic shift vectors corresponding to atom pairs
     * \param[in] positionsMM MM atom positions
     * \param[in] chargesMM MM atom charges
     * \param[in] nnpCharge total charge of NNP region
     * \param[in] linkFrontier link frontier atoms
     * \param[in] box simulation box
     * \param[in] pbcType periodic boundary conditions
     */
    void evaluateModel(gmx_enerdata_t*                  enerd,
                       ArrayRef<RVec>                   forces,
                       ArrayRef<const int>              indexLookup,
                       ArrayRef<const int>              mmIndices,
                       ArrayRef<const std::string>      inputs,
                       ArrayRef<RVec>                   positions,
                       ArrayRef<int>                    atomNumbers,
                       ArrayRef<int>                    atomPairs,
                       ArrayRef<RVec>                   pairShifts,
                       ArrayRef<RVec>                   positionsMM,
                       ArrayRef<real>                   chargesMM,
                       real                             nnpCharge,
                       ArrayRef<const LinkFrontierAtom> linkFrontier,
                       matrix&                          box,
                       PbcType&                         pbcType) override;

    //! helper function to check if model outputs forces
    bool outputsForces() const override;

private:
    //! Functions to prepare inputs for NN model. Create input torch::Tensors for the model.
    //! \{
    void prepareAtomPositions(ArrayRef<RVec> positions);
    void prepareAtomNumbers(ArrayRef<int> atomTypes);
    void prepareBox(matrix& box);
    void preparePbcType(PbcType& pbcType);
    void prepareAtomPairs(ArrayRef<int> atomPairs);
    void preparePairShifts(ArrayRef<RVec> pairShifts);
    void prepareMMPositions(ArrayRef<RVec> pos);
    void prepareMMCharges(ArrayRef<real> charges);
    void prepareNNPCharge(real charge);
    //! \}

    //! pointer to the communication object
    const MpiComm& mpiComm_;
    //! MDLogger during mdrun
    const MDLogger& logger_;

    const NNPotEmbedding embeddingScheme_;

    //! device to run the model on
    torch::Device device_;
    //! TorchScript model
    torch::jit::script::Module model_;
    //! input tensors for the model
    std::vector<torch::jit::IValue> inputs_;
    //! output tensors for the model
    c10::intrusive_ptr<c10::ivalue::Tuple> outputs_;
};

} // namespace gmx

#endif
