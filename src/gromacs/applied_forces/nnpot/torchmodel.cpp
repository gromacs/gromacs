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
 * Implements a wrapper class for a TorchScript-compiled neural network model.
 *
 * \author Lukas Müllender <lukas.muellender@gmail.com>
 * \ingroup module_applied_forces
 */

#include "torchmodel.h"

#include <dlfcn.h>

#include <iostream>

#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/mpicomm.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

//! \brief \internal Helper function to convert a torch::DataPtr to a string.
static std::string recordToString(std::tuple<at::DataPtr, size_t> data)
{
    return std::string(static_cast<char*>(std::get<0>(data).get()), std::get<1>(data));
}

//! \brief load custom extensions used to compile the TorchScript model from extra_files
static void loadModelExtensions(const std::string& ext_libs, const MDLogger& logger)
{
    if (!ext_libs.empty())
    {
        // split string by colon
        std::vector<std::string> ext_libs_vec = gmx::splitDelimitedString(ext_libs, ':');

        // load extension libraries
        for (const auto& lib : ext_libs_vec)
        {
            auto* ext_lib = dlopen(lib.c_str(), RTLD_LAZY);
            if (!ext_lib)
            {
                GMX_LOG(logger.warning)
                        .appendTextFormatted(
                                "Could not load pytorch extension at %s. Error was: '%s'",
                                lib.c_str(),
                                dlerror());
            }
            else
            {
                GMX_LOG(logger.info).appendText("Loaded PyTorch extension library at " + lib);
            }
        }
    }
}

/*! \brief Determine which device to use depending on GMX_NN_DEVICE environment variable
 *
 * Defaults to CPU if no environment variable is set. Throws an error if CUDA is requested but no CUDA device is available.
 * Also throws an error if the environment variable is set to an invalid value.
 */
static torch::DeviceType determineDeviceType(const MDLogger& logger)
{
    torch::DeviceType device;

    // check if environment variable GMX_NN_DEVICE is set (should be something like cuda:0 or cpu)
    if (const char* env = std::getenv("GMX_NN_DEVICE"))
    {
        std::string dev = std::string(env);
        if (dev == "cpu")
        {
            device = torch::kCPU;
        }
        else if (dev.find("cuda") != std::string::npos)
        {
            if (!torch::cuda::is_available())
            {
                GMX_THROW(
                        InternalError("Environment variable GMX_NN_DEVICE was set to cuda, but no "
                                      "CUDA device is available."));
            }
            device = torch::kCUDA;
        }
        else
        {
            GMX_THROW(InvalidInputError(
                    "Environment variable GMX_NN_DEVICE was set to an invalid value."));
        }
        GMX_LOG(logger.info)
                .asParagraph()
                .appendText(
                        "GMX_NN_DEVICE environment variable found."
                        " Using device: "
                        + dev + ".");
    }
    else
    {
        // todo: default to gpu if available
        device = torch::kCPU;
        GMX_LOG(logger.info)
                .asParagraph()
                .appendText(
                        "No GMX_NN_DEVICE environment variable found."
                        " Using CPU.");
    }
    return device;
}

TorchModel::TorchModel(const std::string& fileName, const MDLogger& logger) :
    logger_(logger), device_(determineDeviceType(logger))
{
    // load pytorch extensions from extra_files
    // we can't use extra_files map from the torch::jit::load function
    // because the extensions are required before loading the model
    auto reader = caffe2::serialize::PyTorchStreamReader(fileName);
    if (reader.hasRecord("extra/extension_libs"))
    {
        std::string ext_libs = recordToString(reader.getRecord("extra/extension_libs"));
        loadModelExtensions(ext_libs, logger_);
    }

    // try loading the model
    try
    {
        model_ = torch::jit::load(fileName, device_);
        model_.eval();
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
}

TorchModel::~TorchModel() {}

void TorchModel::prepareAtomPositions(ArrayRef<RVec>& positions)
{
    int N = positions.size();

    auto options = torch::TensorOptions().dtype(torchRealType).requires_grad(true);
    // preferrable to use from_blob here because it doesn't allocate new memory
    torch::Tensor posTensor = torch::from_blob(positions.data()->as_vec(), { N, DIM }, options);
    posTensor               = posTensor.to(torch::kFloat32);

    // important to set requires_grad to true again after moving to GPU
    // otherwise, the tensor is a leaf node of the computation graph anymore
    posTensor = posTensor.to(device_).requires_grad_(true);

    inputs_.push_back(posTensor);
}

void TorchModel::prepareAtomNumbers(ArrayRef<int>& atomTypes)
{
    const int     N = atomTypes.size();
    torch::Tensor typesTensor =
            torch::from_blob(atomTypes.data(), { N }, torch::TensorOptions().dtype(torch::kInt32));
    typesTensor = typesTensor.to(device_);

    inputs_.push_back(typesTensor);
}

void TorchModel::prepareBox(matrix* box)
{
    torch::Tensor boxTensor =
            torch::from_blob(*box, { DIM, DIM }, torch::TensorOptions().dtype(torchRealType));
    boxTensor = boxTensor.to(torch::kFloat32).to(device_);
    inputs_.push_back(boxTensor);
}

void TorchModel::preparePbcType(PbcType* pbcType)
{
    torch::Tensor pbcTensor =
            torch::tensor({ true, true, true }, torch::TensorOptions().dtype(torch::kBool));
    if (*pbcType == PbcType::XY)
    {
        pbcTensor[2] = false;
    }
    else if (*pbcType != PbcType::Xyz)
    {
        GMX_THROW(InconsistentInputError(
                "Option use_pbc was set to true, but PBC type is not supported."));
    }
    pbcTensor = pbcTensor.to(device_);
    inputs_.push_back(pbcTensor);
}

void TorchModel::evaluateModel(gmx_enerdata_t*              enerd,
                               const ArrayRef<RVec>&        forces,
                               ArrayRef<const int>&         indexLookup,
                               ArrayRef<const std::string>& inputs,
                               ArrayRef<RVec>&              positions,
                               ArrayRef<int>&               atomNumbers,
                               matrix*                      box /* = nullptr*/,
                               PbcType*                     pbcType /* = nullptr*/)
{
    // prepare inputs for NN model
    // order in input vector is the same as in mdp file
    for (const std::string& input : inputs)
    {
        if (input.empty())
        {
            continue;
        }
        else if (input == "atom-positions")
        {
            prepareAtomPositions(positions);
        }
        else if (input == "atom-numbers")
        {
            prepareAtomNumbers(atomNumbers);
        }
        else if (input == "box")
        {
            prepareBox(box);
        }
        else if (input == "pbc")
        {
            preparePbcType(pbcType);
        }
        else
        {
            GMX_THROW(InconsistentInputError("Unknown input to NN model: " + input));
        }
    }

    // for now: only do inference on main rank
    if (mpiComm_->isMainRank())
    {
        // run model
        c10::IValue out = model_.forward(inputs_);
        outputs_ = out.isTuple() ? out.toTuple() : c10::ivalue::Tuple::create({ out.toTensor() });
    }
    else
    {
        // empty tuple for non-main ranks
        outputs_ = c10::ivalue::Tuple::create({});
    }

    // check if model outputs forces after the forward pass
    const bool modelOutputsForces = outputsForces();

    // for now: only evaluate gradient on main rank
    torch::Tensor energyTensor;
    const int     N           = indexLookup.size();
    torch::Tensor forceTensor = torch::zeros({ N, DIM });

    if (mpiComm_->isMainRank())
    {
        // get energy
        energyTensor = outputs_->elements()[0].toTensor();

        // if the model outputs forces, retrieve them
        // otherwise, we calculate them from the energy
        if (modelOutputsForces)
        {
            forceTensor = outputs_->elements()[1].toTensor();
        }
        else
        {
            // compute gradient wrt 1st input (positions) to get forces
            // .retain_grad() is necessary to keep the gradient in the computation graph!
            // not sure why, since inputs_[0].toTensor().requires_grad() is still true
            inputs_[0].toTensor().retain_grad();
            energyTensor.backward();
            forceTensor = -1. * inputs_[0].toTensor().grad().to(torch::kCPU);
        }
        // set energy
        energyTensor          = energyTensor.to(torchRealType).to(torch::kCPU);
        enerd->term[F_ENNPOT] = energyTensor.item<real>();

        forceTensor = forceTensor.to(torchRealType).to(torch::kCPU);
    }

    // distribute forces
    if (mpiComm_->isParallel())
    {
        mpiComm_->sumReduce(3 * N, static_cast<real*>(forceTensor.data_ptr()));
    }

    // accumulate forces only on local atoms
    auto forceAccessor = forceTensor.accessor<real, 2>();
    for (int m = 0; m < DIM; ++m)
    {
        for (int i = 0; i < N; ++i)
        {
            // if value in lookup table is -1, the atom is not local
            if (indexLookup[i] == -1)
            {
                continue;
            }
            forces[indexLookup[i]][m] += forceAccessor[i][m];
        }
    }

    // after output is retrieved: reset input vector
    inputs_.resize(0);
}

void TorchModel::setComm(const MpiComm& mpiComm)
{
    mpiComm_ = &mpiComm;
}

bool TorchModel::outputsForces() const
{
    return outputs_->elements().size() > 1;
}

} // namespace gmx
