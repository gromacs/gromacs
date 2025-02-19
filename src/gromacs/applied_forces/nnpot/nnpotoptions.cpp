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
 * Implements the options for NNPot MDModule class.
 *
 * \author Lukas Müllender <lukas.muellender@gmail.com>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "nnpotoptions.h"

#include "gromacs/domdec/localatomset.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"

#include "nnpottopologypreprocessor.h"

#if GMX_TORCH
#    include "gromacs/applied_forces/nnpot/torchmodel.h"
#endif

namespace gmx
{

/*! \brief \internal Helper to declare mdp transform rules.
 *
 * Enforces uniform mdp options that are always prepended with the correct
 * string for the NNPot mdp options.
 *
 * \tparam ToType type to be transformed to
 * \tparam TransformWithFunctionType type of transformation function to be used
 *
 * \param[in] rules KVT transformation rules
 * \param[in] transformationFunction the function to transform the flat kvt tree
 * \param[in] optionTag string tag that describes the mdp option, appended to the
 *                      default string for the QMMM simulation
 */
template<class ToType, class TransformWithFunctionType>
void NNPotMdpTransformFromString(IKeyValueTreeTransformRules* rules,
                                 TransformWithFunctionType    transformationFunction,
                                 const std::string&           optionTag)
{
    rules->addRule()
            .from<std::string>("/" + c_nnpotModuleName + "-" + optionTag)
            .to<ToType>("/" + c_nnpotModuleName + "/" + optionTag)
            .transformWith(transformationFunction);
}

void NNPotOptions::initMdpTransform(IKeyValueTreeTransformRules* rules)
{
    const auto& stringIdentityTransform = [](std::string s) { return s; };
    NNPotMdpTransformFromString<bool>(rules, &fromStdString<bool>, c_activeTag_);
    NNPotMdpTransformFromString<std::string>(rules, stringIdentityTransform, c_modelFileNameTag_);
    NNPotMdpTransformFromString<std::string>(rules, stringIdentityTransform, c_inputGroupTag_);
    NNPotMdpTransformFromString<std::string>(rules, stringIdentityTransform, c_modelInput1Tag_);
    NNPotMdpTransformFromString<std::string>(rules, stringIdentityTransform, c_modelInput2Tag_);
    NNPotMdpTransformFromString<std::string>(rules, stringIdentityTransform, c_modelInput3Tag_);
    NNPotMdpTransformFromString<std::string>(rules, stringIdentityTransform, c_modelInput4Tag_);
}

void NNPotOptions::initMdpOptions(IOptionsContainerWithSections* options)
{
    auto section = options->addSection(OptionSection(c_nnpotModuleName.c_str()));
    section.addOption(BooleanOption(c_activeTag_.c_str()).store(&params_.active_));
    section.addOption(StringOption(c_modelFileNameTag_.c_str()).store(&params_.modelFileName_));
    section.addOption(StringOption(c_inputGroupTag_.c_str()).store(&params_.inputGroup_));
    section.addOption(StringOption(c_modelInput1Tag_.c_str()).store(&params_.modelInput_[0]));
    section.addOption(StringOption(c_modelInput2Tag_.c_str()).store(&params_.modelInput_[1]));
    section.addOption(StringOption(c_modelInput3Tag_.c_str()).store(&params_.modelInput_[2]));
    section.addOption(StringOption(c_modelInput4Tag_.c_str()).store(&params_.modelInput_[3]));
}

void NNPotOptions::buildMdpOutput(KeyValueTreeObjectBuilder* builder) const
{
    // new empty line before writing nnpot mdp values
    builder->addValue<std::string>("comment-" + c_nnpotModuleName + "empty-line", "");

    builder->addValue<std::string>("comment-" + c_nnpotModuleName + "-module",
                                   "; Neural Network potential");
    builder->addValue<bool>(c_nnpotModuleName + "-" + c_activeTag_, params_.active_);

    if (params_.active_)
    {
        builder->addValue<std::string>(c_nnpotModuleName + "-" + c_modelFileNameTag_, params_.modelFileName_);
        builder->addValue<std::string>(c_nnpotModuleName + "-" + c_inputGroupTag_, params_.inputGroup_);
        builder->addValue<std::string>(c_nnpotModuleName + "-" + c_modelInput1Tag_,
                                       params_.modelInput_[0]);
        builder->addValue<std::string>(c_nnpotModuleName + "-" + c_modelInput2Tag_,
                                       params_.modelInput_[1]);
        builder->addValue<std::string>(c_nnpotModuleName + "-" + c_modelInput3Tag_,
                                       params_.modelInput_[2]);
        builder->addValue<std::string>(c_nnpotModuleName + "-" + c_modelInput4Tag_,
                                       params_.modelInput_[3]);
    }
}

bool NNPotOptions::isActive() const
{
    return params_.active_;
}

std::string NNPotOptions::getModelFileName() const
{
    return params_.modelFileName_;
}

void NNPotOptions::setInputGroupIndices(const IndexGroupsAndNames& indexGroupsAndNames)
{
    // Exit if NNPot module is not active
    if (!params_.active_)
    {
        return;
    }

    // Create input index
    params_.inpIndices_ = indexGroupsAndNames.indices(params_.inputGroup_);

    // Check that group is not empty
    if (params_.inpIndices_.empty())
    {
        GMX_THROW(InconsistentInputError(
                formatString("Group %s defining NN potential input atoms should not be empty.",
                             params_.inputGroup_.c_str())));
    }

    // Create temporary index for the whole System
    auto systemIndices = indexGroupsAndNames.indices("System");

    // Sort inpIndices_ and sysIndices_
    std::sort(params_.inpIndices_.begin(), params_.inpIndices_.end());
    std::sort(systemIndices.begin(), systemIndices.end());

    // Create MM index
    params_.mmIndices_.reserve(systemIndices.size());

    // Position in inpIndices_
    size_t j = 0;
    // Write to mmIndices_ only the atoms which do not belong to NNP input region
    for (size_t i = 0; i < systemIndices.size(); i++)
    {
        if (systemIndices[i] != params_.inpIndices_[j])
        {
            params_.mmIndices_.push_back(systemIndices[i]);
        }
        else
        {
            if (j < params_.inpIndices_.size() - 1)
            {
                j++;
            }
        }
    }
}

void NNPotOptions::setLocalInputAtomSet(const LocalAtomSet& localInputAtomSet)
{
    params_.inpAtoms_ = std::make_unique<LocalAtomSet>(localInputAtomSet);
}

void NNPotOptions::setLocalMMAtomSet(const LocalAtomSet& localMMAtomSet)
{
    params_.mmAtoms_ = std::make_unique<LocalAtomSet>(localMMAtomSet);
}

void NNPotOptions::modifyTopology(gmx_mtop_t* top)
{
    // Exit if module is not active
    if (!params_.active_)
    {
        return;
    }

    // subclassing the qmmm topology preprocessor as it has virtually the exact functionality we need
    //! \todo separate this from QMMM module as its reused in multiple places now
    NNPotTopologyPreprocessor topPrep(params_.inpIndices_);
    topPrep.preprocess(top);

    // Get info about modifications
    QMMMTopologyInfo topInfo = topPrep.topInfo();

    // Check that logger and warning handler are valid
    GMX_ASSERT(logger_, "Logger not set.");
    GMX_ASSERT(wi_, "WarningHandler not set.");

    // Inform the user about performed modifications, issue warning if necessary
    GMX_LOG(logger_->info)
            .appendText("Neural network potential Interface is active, topology was modified!");

    GMX_LOG(logger_->info)
            .appendTextFormatted("Number of NN input atoms: %d\nNumber of regular atoms: %d",
                                 topInfo.numQMAtoms,
                                 topInfo.numMMAtoms);

    if (topInfo.numBondsRemoved > 0)
    {
        GMX_LOG(logger_->info).appendTextFormatted("Bonds removed: %d", topInfo.numBondsRemoved);
    }

    if (topInfo.numAnglesRemoved > 0)
    {
        GMX_LOG(logger_->info).appendTextFormatted("Angles removed: %d", topInfo.numAnglesRemoved);
    }

    if (topInfo.numDihedralsRemoved > 0)
    {
        GMX_LOG(logger_->info).appendTextFormatted("Dihedrals removed: %d", topInfo.numDihedralsRemoved);
    }

    if (topInfo.numSettleRemoved > 0)
    {
        GMX_LOG(logger_->info).appendTextFormatted("Settles removed: %d", topInfo.numSettleRemoved);
    }

    if (topInfo.numConnBondsAdded > 0)
    {
        GMX_LOG(logger_->info).appendTextFormatted("Connection-only (type 5) bonds added: %d", topInfo.numConnBondsAdded);
    }

    // Warn in case of broken covalent bonds between NNP input atoms and MM atoms
    if (topInfo.numLinkBonds > 0)
    {
        wi_->addWarning(formatString(
                "%d broken bonds found between NN input atoms and MM atoms."
                " This is can lead to unexpected behavior and is discouraged as of now."
                " Set the -maxwarn option if you want to proceed anyway.",
                topInfo.numLinkBonds));
    }

    // If there are many constrained bonds in NNP input region then we should also warn the user
    if (topInfo.numConstrainedBondsInQMSubsystem > 2)
    {
        wi_->addWarning(
                "Your neural network potential subsystem has a lot of constrained bonds. "
                "They probably have been generated automatically. "
                "That could produce artifacts in the simulation. "
                "Consider constraints = none in the mdp file.");
    }
}

void NNPotOptions::writeParamsToKvt(KeyValueTreeObjectBuilder treeBuilder)
{
    // Check if active
    if (!params_.active_)
    {
        return;
    }

    // Write input atom indices
    auto GroupIndexAdder =
            treeBuilder.addUniformArray<std::int64_t>(c_nnpotModuleName + "-" + c_inputGroupTag_);
    for (const auto& indexValue : params_.inpIndices_)
    {
        GroupIndexAdder.addValue(indexValue);
    }

    // Write MM atoms index
    GroupIndexAdder = treeBuilder.addUniformArray<std::int64_t>(c_nnpotModuleName + "-" + c_mmGroupTag_);
    for (const auto& indexValue : params_.mmIndices_)
    {
        GroupIndexAdder.addValue(indexValue);
    }

    // check that the model and model inputs are valid
    checkNNPotModel();
}

void NNPotOptions::readParamsFromKvt(const KeyValueTreeObject& tree)
{
    // Check if active
    if (!params_.active_)
    {
        return;
    }

    // Try to read input atoms index
    std::string key = c_nnpotModuleName + "-" + c_inputGroupTag_;
    if (!tree.keyExists(key))
    {
        GMX_THROW(InconsistentInputError(
                "Cannot find input atoms index vector required for neural network potential.\n"
                "This could be caused by incompatible or corrupted tpr input file."));
    }
    auto kvtIndexArray = tree[key].asArray().values();
    params_.inpIndices_.resize(kvtIndexArray.size());
    std::transform(std::begin(kvtIndexArray),
                   std::end(kvtIndexArray),
                   std::begin(params_.inpIndices_),
                   [](const KeyValueTreeValue& val) { return val.cast<std::int64_t>(); });

    // Try to read MM atoms index
    key = c_nnpotModuleName + "-" + c_mmGroupTag_;
    if (!tree.keyExists(key))
    {
        GMX_THROW(InconsistentInputError(
                "Cannot find MM atoms index vector required for neural network potential.\n"
                "This could be caused by incompatible or corrupted tpr input file."));
    }
    kvtIndexArray = tree[key].asArray().values();
    params_.mmIndices_.resize(kvtIndexArray.size());
    std::transform(std::begin(kvtIndexArray),
                   std::end(kvtIndexArray),
                   std::begin(params_.mmIndices_),
                   [](const KeyValueTreeValue& val) { return val.cast<std::int64_t>(); });
}

const NNPotParameters& NNPotOptions::parameters()
{
    return params_;
}

void NNPotOptions::setTopology(const gmx_mtop_t& top)
{
    params_.atoms_    = gmx_mtop_global_atoms(top);
    params_.numAtoms_ = params_.atoms_.nr;
}

void NNPotOptions::setPbcType(const PbcType& pbcType)
{
    params_.pbcType_ = std::make_unique<PbcType>(pbcType);
}

void NNPotOptions::setLogger(const MDLogger& logger)
{
    logger_ = &logger;
}

void NNPotOptions::setCommRec(const t_commrec& cr)
{
    params_.cr_ = &cr;
}

const MDLogger* NNPotOptions::logger()
{
    return logger_;
}

void NNPotOptions::setWarninp(WarningHandler* wi)
{
    wi_ = wi;
}

#if !GMX_TORCH
[[noreturn]]
#endif
void NNPotOptions::checkNNPotModel()
{
#if GMX_TORCH

    // try initializing the neural network model
    std::filesystem::path        modelPath(params_.modelFileName_);
    std::unique_ptr<INNPotModel> model;
    if (!std::filesystem::exists(modelPath))
    {
        GMX_THROW(FileIOError("Model file does not exist: " + params_.modelFileName_));
    }
    else if (modelPath.extension() == ".pt")
    {
        model = std::make_unique<TorchModel>(params_.modelFileName_, logger_);
    }
    else
    {
        GMX_THROW(FileIOError("Unrecognized extension for model file: " + params_.modelFileName_));
    }
    model->initModel();

    // check if model accepts inputs
    // prepare dummy inputs for NN model
    int               inputsFound = 0;
    std::vector<RVec> positions(1, RVec({ 0.0, 0.0, 0.0 }));
    std::vector<int>  atomNumbers(1, 1);
    matrix            box = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };
    PbcType           pbc = PbcType();
    t_commrec         cr;
    model->setCommRec(&cr);
    for (const std::string& input : params_.modelInput_)
    {
        if (input.empty())
        {
            continue;
        }
        else if (input == "atom-positions")
        {
            inputsFound++;
            model->prepareAtomPositions(positions);
        }
        else if (input == "atom-numbers")
        {
            inputsFound++;
            model->prepareAtomNumbers(atomNumbers);
        }
        else if (input == "box")
        {
            inputsFound++;
            model->prepareBox(box);
        }
        else if (input == "pbc")
        {
            inputsFound++;
            model->preparePbcType(pbc);
        }
        else
        {
            GMX_THROW(InconsistentInputError("Unknown input to NN model: " + input));
        }
    }
    // check that inputs are not empty
    if (inputsFound == 0)
    {
        GMX_THROW(InconsistentInputError("No inputs to NN model provided."));
    }

    // Check that logger and warning handler are valid
    GMX_ASSERT(logger_, "Logger not set.");
    GMX_ASSERT(wi_, "WarningHandler not set.");
    bool modelOutputsForces = false;
    try
    {
        // might throw a runtime error if the model is not compatible with the dummy input
        model->evaluateModel();

        // check if model outputs forces after forward pass
        modelOutputsForces = model->outputsForces();

        // prepare dummy output
        std::vector<int>  indices(1, 0);
        gmx_enerdata_t    enerd(1, nullptr);
        std::vector<RVec> forcesVec(1, RVec({ 0.0, 0.0, 0.0 }));
        ArrayRef<RVec>    forces(forcesVec);
        model->getOutputs(indices, enerd, forces);

        // log wheter model outputs forces or we need to compute them
        if (modelOutputsForces)
        {
            GMX_LOG(logger_->info).appendText("Will use forces from NNP model.");
        }
        else
        {
            GMX_LOG(logger_->info)
                    .appendText(
                            "NNP model does not output forces. They will be computed as gradients "
                            "of the energy w.r.t. the first input tensor (atom positions).");
        }
    }
    catch (const GromacsException& e)
    {
        // rethrow exceptions issued by our code
        throw;
    }
    catch (const std::exception& e)
    {
        // we only issue a warning here instead of throwing an error, as a torch runtime error might
        // simply be due to mismatched dummy input shapes
        wi_->addWarning("There was an error while checking NN model with a dummy input: " + std::string(e.what()) + "\n"
                        "I can't verify that the model works correctly. This might lead to errors during mdrun.");
    }

    // check if first input is atom-positions
    if ((params_.modelInput_[0] != "atom-positions") && !modelOutputsForces)
    {
        wi_->addWarning("Gradients will be computed with respect to first input to NN model "
                        + params_.modelInput_[0] + " instead of atom positions. Is this intended?");
    }

#else
    GMX_THROW(InternalError(
            "Libtorch/NN backend is not linked into GROMACS, NNPot simulation is not possible."
            " Please, reconfigure GROMACS with -DGMX_NNPOT=TORCH\n"));
#endif // GMX_TORCH
}

} // namespace gmx
