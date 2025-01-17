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

/*! \brief Define the torch datatype according to GMX_DOUBLE.
 *
 * Important for converting data types, as model inference is always done in float32.
 */
static constexpr auto torchRealType = GMX_DOUBLE ? torch::kFloat64 : torch::kFloat32;

class MDLogger;

/*! \brief
 * Class responsible for loading and evaluating a TorchScript-compiled neural network model.
 * Inherits from NNPotModel.
 */
class TorchModel : public INNPotModel
{
public:
    /*! \brief Constructor for TorchModel.
     * \param[in] filename path to the TorchScript model file
     * \param[in] logger pointer to the MDLogger
     */
    TorchModel(const std::string& filename, const MDLogger* logger);

    ~TorchModel();

    //! Initialize the neural network model
    void initModel() override;

    //! Functions to prepare inputs for NN model. Create input torch::Tensors for the model.
    //! \{
    void prepareAtomPositions(std::vector<RVec>& positions) override;
    void prepareAtomNumbers(std::vector<int>& atomTypes) override;
    void prepareBox(matrix& box) override;
    void preparePbcType(PbcType& pbcType) override;
    //! \}

    //! Call inference on NN model
    void evaluateModel() override;

    //! Retrieve NN model outputs
    void getOutputs(std::vector<int>& indices, gmx_enerdata_t& enerd, const ArrayRef<RVec>& forces) override;

    //! Set communication record for possible communication of input/output data between ranks
    void setCommRec(const t_commrec* cr) override;

    //! helper function to check if model outputs forces
    bool outputsForces() const override;

    //! determine which device to use depending on GMX_NN_DEVICE environment variable
    void setDevice();

    //! load custom extensions used to compile the TorchScript model from extra_files
    void loadModelExtensions(std::string& extension_libs);

private:
    //! flag to check if model is initialized
    bool isInit_ = false;
    //! flag to check if model output is ready
    bool outputReady_ = false;

    //! pointer to the communication record
    const t_commrec* cr_ = nullptr;
    //! pointer to the MDLogger
    const MDLogger* logger_ = nullptr;

    //! TorchScript model
    std::shared_ptr<torch::Device> device_;
    //! device to run the model on
    torch::jit::script::Module model_;
    //! input tensors for the model
    std::vector<torch::jit::IValue> inputs_;
    //! output tensors for the model
    c10::intrusive_ptr<c10::ivalue::Tuple> outputs_;
};

} // namespace gmx

#endif
