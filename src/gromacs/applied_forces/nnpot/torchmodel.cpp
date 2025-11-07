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
#include "gromacs/hardware/device_information.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/embedded_system_preprocessing.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/mpicomm.h"
#include "gromacs/utility/stringutil.h"

#include "nnpotoptions.h"

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

/*! \brief \internal Helper function to get active device used by mdrun
 *
 * This should be called before loading the model, or calling any torch functions that might
 * overwrite the active device.
 * \returns Tuple of device type string and optional device index.
 */
static std::tuple<std::string, std::optional<int>> getMdrunActiveDevice()
{
#if GMX_GPU_CUDA || (GMX_SYCL_ACPP && GMX_ACPP_HAVE_CUDA_TARGET) // AdaptiveCpp uses CUDA Runtime API
    // Annoyingly, USE_CUDA doesn't seem to be defined correctly by libtorch when using CUDA,
    // while for HIP, torch::hasHIP() doesn't work, so we have different checks for both cases.
    GMX_RELEASE_ASSERT(
            torch::hasCUDA(),
            "Libtorch was not compiled with CUDA support. Ensure that GROMACS and Libtorch "
            "are compiled with matching GPU support.");
    int  activeDevice;
    auto success = cudaGetDevice(&activeDevice);
    if (success != cudaSuccess)
    {
        GMX_THROW(InternalError("Could not get active device: cudaGetDevice failed."));
    }
    return { "cuda", activeDevice };
#elif GMX_GPU_HIP || (GMX_SYCL_ACPP && GMX_ACPP_HAVE_HIP_TARGET) // AdaptiveCpp uses HIP Runtime API
#    ifndef USE_ROCM                                             // Defined by torch
    GMX_THROW(InternalError(
            "Libtorch was not compiled with HIP support. Ensure that GROMACS and Libtorch are "
            "compiled with matching GPU support."));
#    endif
    int  activeDevice;
    auto success = hipGetDevice(&activeDevice);
    if (success != hipSuccess)
    {
        GMX_THROW(InternalError("Could not get active device: hipGetDevice failed."));
    }
    return { "hip", activeDevice };
#else
    return { "cpu", std::nullopt };
#endif
}

/*! \brief Determine which device to use depending on GMX_NN_DEVICE environment variable
 *
 * Should be "gpu" or "cpu". "cuda" option is deprecated. Defaults to GPU if no environment variable is set.
 * Throws an error if CUDA/HIP is requested but no compatible device is available.
 * Also throws an error if the environment variable is set to an invalid value.
 */
static torch::Device determineDevice(const MDLogger& logger, const MpiComm& mpiComm)
{
    torch::Device device(torch::kCPU);

    // if not on the main rank, just return CPU device
    if (!mpiComm.isMainRank())
    {
        return device;
    }

    // check active device and type
    auto [torchDeviceType, activeDevice] = getMdrunActiveDevice();

    // check if environment variable GMX_NN_DEVICE is set (should be gpu or cpu)
    if (const char* env = std::getenv("GMX_NN_DEVICE"))
    {
        const std::string devLC = toLowerCase(env);
        if (devLC == "gpu" || devLC == "cuda")
        {
            if (!torch::cuda::is_available())
            {
                GMX_THROW(InternalError(formatString(
                        "GMX_NN_DEVICE was set to '%s', but no matching device is available.", env)));
            }
            GMX_RELEASE_ASSERT(
                    activeDevice.has_value(),
                    "Active device couldn't be determined. Ensure that Libtorch and GROMACS "
                    "are compiled with matching GPU support.");
            device = torch::Device(torch::kCUDA, activeDevice.value());
        }
        else if (devLC != "cpu")
        {
            GMX_THROW(InvalidInputError(formatString(
                    "Environment variable GMX_NN_DEVICE was set to an invalid value: '%s'.", env)));
        }
        GMX_LOG(logger.info)
                .asParagraph()
                .appendTextFormatted(
                        "GMX_NN_DEVICE environment variable found. Using device: '%s'.", env);
    }
    else // not set: automatically determine device
    {
        if (torch::cuda::is_available() && activeDevice.has_value())
        {
            GMX_LOG(logger.info)
                    .asParagraph()
                    .appendText(
                            "GMX_NN_DEVICE environment variable not found."
                            " Using "
                            + toUpperCase(torchDeviceType) + " device.");
            device = torch::Device(torch::kCUDA, activeDevice.value());
        }
        else
        {
            GMX_LOG(logger.info)
                    .asParagraph()
                    .appendText(
                            "GMX_NN_DEVICE environment variable not found, and no GPU available."
                            " Using CPU.");
        }
    }
    return device;
}

TorchModel::TorchModel(const std::string& fileName,
                       NNPotEmbedding     embedding,
                       const MDLogger&    logger,
                       const MpiComm&     mpiComm) :
    mpiComm_(mpiComm), logger_(logger), embeddingScheme_(embedding), device_(determineDevice(logger, mpiComm))
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
    if (mpiComm_.isMainRank())
    {
        try
        {
            model_ = torch::jit::load(fileName, device_);
            model_.eval();
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }
}

void TorchModel::prepareAtomPositions(ArrayRef<RVec> positions)
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

void TorchModel::prepareAtomNumbers(ArrayRef<int> atomTypes)
{
    const int     N = atomTypes.size();
    torch::Tensor typesTensor =
            torch::from_blob(atomTypes.data(), { N }, torch::TensorOptions().dtype(torch::kInt32));
    typesTensor = typesTensor.to(device_);

    inputs_.push_back(typesTensor);
}

void TorchModel::prepareMMPositions(ArrayRef<RVec> positions)
{
    int N = positions.size();

    auto options = torch::TensorOptions().dtype(torchRealType).requires_grad(true);
    // preferable to use from_blob here because it doesn't allocate new memory
    torch::Tensor posTensor = torch::from_blob(positions.data()->as_vec(), { N, DIM }, options);
    posTensor               = posTensor.to(torch::kFloat32);

    // important to set requires_grad to true again after moving to GPU
    // otherwise, the tensor is not a leaf node of the computation graph anymore
    posTensor = posTensor.to(device_).requires_grad_(true);

    inputs_.push_back(posTensor);
}

void TorchModel::prepareMMCharges(ArrayRef<real> charges)
{
    const int     N = charges.size();
    torch::Tensor chargesTensor =
            torch::from_blob(charges.data(), { N }, torch::TensorOptions().dtype(torchRealType));
    chargesTensor = chargesTensor.to(torch::kFloat32).to(device_).requires_grad_(true);

    inputs_.push_back(chargesTensor);
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

void TorchModel::prepareAtomPairs(ArrayRef<int> atomPairs)
{
    const int     numPairs    = atomPairs.size() / 2;
    torch::Tensor pairsTensor = torch::from_blob(
            atomPairs.data(), { numPairs, 2 }, torch::TensorOptions().dtype(torch::kInt32));
    pairsTensor = pairsTensor.to(device_);

    inputs_.push_back(pairsTensor);
}

void TorchModel::preparePairShifts(ArrayRef<RVec> pairShifts)
{
    const int     N            = pairShifts.size();
    torch::Tensor shiftsTensor = torch::from_blob(
            pairShifts.data()->as_vec(), { N, DIM }, torch::TensorOptions().dtype(torchRealType));
    shiftsTensor = shiftsTensor.to(device_);
    inputs_.push_back(shiftsTensor);
}

void TorchModel::prepareNNPCharge(real charge)
{
    torch::Tensor chargeTensor = torch::tensor({ charge }, torch::TensorOptions().dtype(torchRealType));
    chargeTensor = chargeTensor.to(torch::kFloat32).to(device_);
    inputs_.push_back(chargeTensor);
}

void TorchModel::evaluateModel(gmx_enerdata_t*                  enerd,
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
                               matrix*                          box /* = nullptr*/,
                               PbcType*                         pbcType /* = nullptr*/)
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
        else if (input == "atom-pairs")
        {
            prepareAtomPairs(atomPairs);
        }
        else if (input == "pair-shifts")
        {
            preparePairShifts(pairShifts);
        }
        else if (input == "atom-positions-mm")
        {
            prepareMMPositions(positionsMM);
        }
        else if (input == "atom-charges-mm")
        {
            prepareMMCharges(chargesMM);
        }
        else if (input == "nnp-charge")
        {
            prepareNNPCharge(nnpCharge);
        }
        else
        {
            GMX_THROW(InconsistentInputError("Unknown input to NN model: " + input));
        }
    }

    // for now: only do inference on main rank
    if (mpiComm_.isMainRank())
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
    const int     numMM       = mmIndices.size();
    torch::Tensor forceTensor = torch::zeros({ N, DIM });
    // electrostatic embedding produces additional force corrections on MM atoms
    torch::Tensor mmForceTensor = torch::zeros({ numMM, DIM });

    if (mpiComm_.isMainRank())
    {
        // get energy
        energyTensor = outputs_->elements()[0].toTensor();

        // if the model outputs forces, retrieve them
        // otherwise, we calculate them from the energy
        if (modelOutputsForces)
        {
            forceTensor = outputs_->elements()[1].toTensor().to(torch::kCPU);
            if (embeddingScheme_ == NNPotEmbedding::ElectrostaticModel)
            {
                mmForceTensor = outputs_->elements()[2].toTensor().to(torch::kCPU);
            }
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
        energyTensor = energyTensor.to(torchRealType).to(torch::kCPU);
        enerd->term[InteractionFunction::NeuralNetworkPotentialEnergy] = energyTensor.item<real>();

        forceTensor = forceTensor.to(torchRealType).to(torch::kCPU);
    }

    // distribute forces
    if (mpiComm_.isParallel())
    {
        mpiComm_.sumReduce(3 * N, static_cast<real*>(forceTensor.data_ptr()));
    }

    // accumulate forces only on local atoms
    auto forceAccessor = forceTensor.accessor<real, 2>();
    // redistribute forces on link atoms: force on link atom becomes force on MM atom
    t_pbc pbc;
    if (!linkFrontier.empty())
    {
        GMX_RELEASE_ASSERT(box && pbcType,
                           "PBC information required when using link atoms with NNP");
        // set pbc struct
        set_pbc(&pbc, *pbcType, *box);
    }
    for (const auto& link : linkFrontier)
    {
        // distribute force on link atom to the two adjacent atoms
        RVec fL, fNNP, fMM;
        for (int d = 0; d < DIM; ++d)
        {
            fL[d] = forceAccessor[link.getInputIndexMM()][d];
        }
        std::tie(fNNP, fMM) = link.spreadForce(fL, pbc);
        for (int d = 0; d < DIM; ++d)
        {
            forceAccessor[link.getInputIndexEmb()][d] += fNNP[d];
            forceAccessor[link.getInputIndexMM()][d] = fMM[d];
        }
    }
    for (int i = 0; i < N; ++i)
    {
        // if value in lookup table is -1, the atom is not local
        if (indexLookup[i] == -1)
        {
            continue;
        }
        for (int d = 0; d < DIM; ++d)
        {
            forces[indexLookup[i]][d] += forceAccessor[i][d];
        }
    }

    // accumulate forces on MM atoms
    if (embeddingScheme_ == NNPotEmbedding::ElectrostaticModel)
    {
        auto mmForceAccessor = mmForceTensor.accessor<real, 2>();
        for (int m = 0; m < DIM; ++m)
        {
            for (int i = 0; i < numMM; ++i)
            {
                forces[mmIndices[i]][m] += mmForceAccessor[i][m];
            }
        }
    }

    // after output is retrieved: reset input vector
    inputs_.resize(0);
}

bool TorchModel::outputsForces() const
{
    return outputs_->elements().size() > 1;
}

} // namespace gmx
