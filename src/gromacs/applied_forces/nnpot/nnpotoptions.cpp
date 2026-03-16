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
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/imdpoptionprovider_helpers.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/embedded_system_preprocessing.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/mpicomm.h"
#include "gromacs/utility/stringutil.h"

#include "nnpot.h"

#if GMX_TORCH
#    include "gromacs/applied_forces/nnpot/torchmodel.h"
#endif

namespace gmx
{

namespace
{

//! Helper function to make a std::string containing the module name
std::string moduleName()
{
    return std::string(NNPotModuleInfo::sc_name);
}

/*! \brief Following Tags denote names of parameters from .mdp file
 * \note Changing these strings will break .tpr backwards compability
 */
//! \{
const std::string c_activeTag_        = "active";
const std::string c_modelFileNameTag_ = "modelfile";
const std::string c_inputGroupTag_    = "input-group";
const std::string c_linkTypeTag_      = "link-type";
const std::string c_linkDistanceTag_  = "link-distance";
const std::string c_pairCutoffTag_    = "pair-cutoff";
const std::string c_nnpChargeTag_     = "nnp-charge";
const std::string c_embeddingTag_     = "embedding";

//! The names of the supported embedding schemes
const EnumerationArray<NNPotEmbedding, const char*> c_embeddingSchemeNames = {
    { "mechanical", "electrostatic-model" }
};

/*! \brief User defined input to NN model.
 *
 *  Possible values:
 * - "atom-positions" vector of atom positions
 * - "atom-numbers" vector of atom numbers
 * - "atom-pairs" list of atom pairs within cutoff set by user
 * - "pair-shifts" list of periodic shifts for atom pairs
 * - "box" unit vectors of simulation box
 * - "pbc" boolean vector indicating periodic boundary conditions
 * - "atom-positions-mm" vector of atom positions
 * - "atom-charges-mm" vector of atomic charges
 * - "nnp-charge" charge of the NNP region
 */
const std::string c_modelInput1Tag_ = "model-input1";
const std::string c_modelInput2Tag_ = "model-input2";
const std::string c_modelInput3Tag_ = "model-input3";
const std::string c_modelInput4Tag_ = "model-input4";
const std::string c_modelInput5Tag_ = "model-input5";
const std::string c_modelInput6Tag_ = "model-input6";
const std::string c_modelInput7Tag_ = "model-input7";
const std::string c_modelInput8Tag_ = "model-input8";
const std::string c_modelInput9Tag_ = "model-input9";
//! \}

/*! \brief Following tags are needed to write parameters generated
 * during preprocessing (grompp) to the .tpr file via KVT
 */
//! \{
const std::string c_mmGroupTag_ = "mm-group";
const std::string c_nnpLinkTag_ = "nnp-link";
const std::string c_mmLinkTag_  = "mm-link";
//! \}

/*! \brief Helper function to preprocess topology for NNP
 *
 * This function performs the following modifications:
 * - Excludes non-bonded interactions between NNP atoms (LJ and Coulomb)
 * - In case of electrostatic embedding: removes classical charges on NNP atoms
 * - Removes bonds containing 1 or more NNP atoms
 * - Removes angles and settles containing 2 or more NNP atoms
 * - Removes dihedrals containing 3 or more NNP atoms
 */
std::vector<LinkFrontierAtom> preprocessNNPotTopology(gmx_mtop_t*           mtop,
                                                      ArrayRef<const Index> nnpIndices,
                                                      const NNPotEmbedding& embedding,
                                                      const real&           nnpCharge,
                                                      const MDLogger&       logger,
                                                      WarningHandler*       wi)
{
    // convert nnpIndices to set for faster lookup
    std::set<int> nnpIndicesSet(nnpIndices.begin(), nnpIndices.end());
    const int     numNNPAtoms     = nnpIndices.size();
    const int     numRegularAtoms = mtop->natoms - numNNPAtoms;

    GMX_LOG(logger.info)
            .appendText("Neural network potential interface is active, topology was modified!");
    GMX_LOG(logger.info)
            .appendTextFormatted("Number of embedded NNP atoms: %d\nNumber of regular atoms: %d\n",
                                 numNNPAtoms,
                                 numRegularAtoms);

    // 1) Split QM-containing molecules from other molecules in blocks
    std::vector<bool> isNNPBlock = splitEmbeddedBlocks(mtop, nnpIndicesSet);

    // 2) Exclude non-bonded interactions between QM atoms
    addEmbeddedNBExclusions(mtop, nnpIndicesSet, logger);

    // 3) Remove classical charges from embedded atoms if electrostatic embedding is used
    if (embedding == NNPotEmbedding::ElectrostaticModel)
    {
        GMX_LOG(logger.info).appendText("Electrostatic embedding scheme is used.\n");
        removeEmbeddedClassicalCharges(mtop, nnpIndicesSet, isNNPBlock, nnpCharge, logger, wi);
    }

    // 4) Make F_CONNBOND between atoms within QM region
    modifyEmbeddedTwoCenterInteractions(mtop, nnpIndicesSet, isNNPBlock, logger);

    // 5) Remove angles and settles containing 2 or more QM atoms
    modifyEmbeddedThreeCenterInteractions(mtop, nnpIndicesSet, isNNPBlock, logger);

    // 6) Remove dihedrals containing 3 or more QM atoms
    modifyEmbeddedFourCenterInteractions(mtop, nnpIndicesSet, isNNPBlock, logger);

    // 7) Check for constrained bonds in subsystem
    checkConstrainedBonds(mtop, nnpIndicesSet, isNNPBlock, wi);

    // 8) Build link frontier information
    std::vector<LinkFrontierAtom> linkFrontier =
            buildLinkFrontier(mtop, nnpIndicesSet, isNNPBlock, logger);

    // finalize topology
    mtop->finalize();

    return linkFrontier;
}

} // namespace

void NNPotOptions::initMdpTransform(IKeyValueTreeTransformRules* rules)
{
    const auto& stringIdentityTransform = [](std::string s) { return s; };
    addMdpTransformFromString<bool>(rules, &fromStdString<bool>, NNPotModuleInfo::sc_name, c_activeTag_);
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, NNPotModuleInfo::sc_name, c_modelFileNameTag_);
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, NNPotModuleInfo::sc_name, c_inputGroupTag_);
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, NNPotModuleInfo::sc_name, c_linkTypeTag_);
    addMdpTransformFromString<real>(
            rules, &fromStdString<real>, NNPotModuleInfo::sc_name, c_linkDistanceTag_);
    addMdpTransformFromString<real>(rules, &fromStdString<real>, NNPotModuleInfo::sc_name, c_pairCutoffTag_);
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, NNPotModuleInfo::sc_name, c_embeddingTag_);
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, NNPotModuleInfo::sc_name, c_modelInput1Tag_);
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, NNPotModuleInfo::sc_name, c_modelInput2Tag_);
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, NNPotModuleInfo::sc_name, c_modelInput3Tag_);
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, NNPotModuleInfo::sc_name, c_modelInput4Tag_);
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, NNPotModuleInfo::sc_name, c_modelInput5Tag_);
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, NNPotModuleInfo::sc_name, c_modelInput6Tag_);
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, NNPotModuleInfo::sc_name, c_modelInput7Tag_);
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, NNPotModuleInfo::sc_name, c_modelInput8Tag_);
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, NNPotModuleInfo::sc_name, c_modelInput9Tag_);
}

void NNPotOptions::initMdpOptions(IOptionsContainerWithSections* options)
{
    auto section = options->addSection(OptionSection(moduleName().c_str()));
    section.addOption(BooleanOption(c_activeTag_.c_str()).store(&params_.active_));
    section.addOption(StringOption(c_modelFileNameTag_.c_str()).store(&params_.modelFileName_));
    section.addOption(StringOption(c_inputGroupTag_.c_str()).store(&params_.inputGroup_));
    section.addOption(StringOption(c_linkTypeTag_.c_str()).store(&params_.linkType_));
    section.addOption(RealOption(c_linkDistanceTag_.c_str()).store(&params_.linkDistance_));
    section.addOption(RealOption(c_pairCutoffTag_.c_str()).store(&params_.pairCutoff_));
    section.addOption(EnumOption<NNPotEmbedding>(c_embeddingTag_.c_str())
                              .enumValue(c_embeddingSchemeNames)
                              .store(&params_.embeddingScheme_));
    section.addOption(RealOption(c_nnpChargeTag_.c_str()).store(&params_.nnpCharge_));
    section.addOption(StringOption(c_modelInput1Tag_.c_str()).store(&params_.modelInput_[0]));
    section.addOption(StringOption(c_modelInput2Tag_.c_str()).store(&params_.modelInput_[1]));
    section.addOption(StringOption(c_modelInput3Tag_.c_str()).store(&params_.modelInput_[2]));
    section.addOption(StringOption(c_modelInput4Tag_.c_str()).store(&params_.modelInput_[3]));
    section.addOption(StringOption(c_modelInput5Tag_.c_str()).store(&params_.modelInput_[4]));
    section.addOption(StringOption(c_modelInput6Tag_.c_str()).store(&params_.modelInput_[5]));
    section.addOption(StringOption(c_modelInput7Tag_.c_str()).store(&params_.modelInput_[6]));
    section.addOption(StringOption(c_modelInput8Tag_.c_str()).store(&params_.modelInput_[7]));
    section.addOption(StringOption(c_modelInput9Tag_.c_str()).store(&params_.modelInput_[8]));
}

void NNPotOptions::buildMdpOutput(KeyValueTreeObjectBuilder* builder) const
{
    // new empty line before writing nnpot mdp values
    addMdpOutputComment(builder, NNPotModuleInfo::sc_name, "empty-line", "");
    addMdpOutputComment(builder, NNPotModuleInfo::sc_name, "module", "; Neural Network potential");
    addMdpOutputValue(builder, NNPotModuleInfo::sc_name, c_activeTag_, params_.active_);

    if (params_.active_)
    {
        addMdpOutputValue<std::string>(
                builder, NNPotModuleInfo::sc_name, c_modelFileNameTag_, params_.modelFileName_);
        addMdpOutputValue<std::string>(
                builder, NNPotModuleInfo::sc_name, c_inputGroupTag_, params_.inputGroup_);
        addMdpOutputValue<std::string>(builder, NNPotModuleInfo::sc_name, c_linkTypeTag_, params_.linkType_);
        addMdpOutputValue<real>(builder, NNPotModuleInfo::sc_name, c_linkDistanceTag_, params_.linkDistance_);
        addMdpOutputValue<real>(builder, NNPotModuleInfo::sc_name, c_pairCutoffTag_, params_.pairCutoff_);
        addMdpOutputValue<std::string>(builder,
                                       NNPotModuleInfo::sc_name,
                                       c_embeddingTag_,
                                       c_embeddingSchemeNames[params_.embeddingScheme_]);
        addMdpOutputValue<std::string>(
                builder, NNPotModuleInfo::sc_name, c_modelInput1Tag_, params_.modelInput_[0]);
        addMdpOutputValue<std::string>(
                builder, NNPotModuleInfo::sc_name, c_modelInput2Tag_, params_.modelInput_[1]);
        addMdpOutputValue<std::string>(
                builder, NNPotModuleInfo::sc_name, c_modelInput3Tag_, params_.modelInput_[2]);
        addMdpOutputValue<std::string>(
                builder, NNPotModuleInfo::sc_name, c_modelInput4Tag_, params_.modelInput_[3]);
        addMdpOutputValue<std::string>(
                builder, NNPotModuleInfo::sc_name, c_modelInput5Tag_, params_.modelInput_[4]);
        addMdpOutputValue<std::string>(
                builder, NNPotModuleInfo::sc_name, c_modelInput6Tag_, params_.modelInput_[5]);
        addMdpOutputValue<std::string>(
                builder, NNPotModuleInfo::sc_name, c_modelInput7Tag_, params_.modelInput_[6]);
        addMdpOutputValue<std::string>(
                builder, NNPotModuleInfo::sc_name, c_modelInput8Tag_, params_.modelInput_[7]);
        addMdpOutputValue<std::string>(
                builder, NNPotModuleInfo::sc_name, c_modelInput9Tag_, params_.modelInput_[8]);
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
    params_.nnpIndices_ = indexGroupsAndNames.indices(params_.inputGroup_);

    // Check that group is not empty
    if (params_.nnpIndices_.empty())
    {
        GMX_THROW(InconsistentInputError(
                formatString("Group %s defining NN potential input atoms should not be empty.",
                             params_.inputGroup_.c_str())));
    }

    // Create temporary index for the whole System
    auto systemIndices = indexGroupsAndNames.indices("System");

    // Sort nnpIndices_ and sysIndices_
    std::sort(params_.nnpIndices_.begin(), params_.nnpIndices_.end());
    std::sort(systemIndices.begin(), systemIndices.end());

    // Create MM index
    params_.mmIndices_.reserve(systemIndices.size());

    // Position in nnpIndices_
    size_t j = 0;
    // Write to mmIndices_ only the atoms which do not belong to NNP input region
    for (size_t i = 0; i < systemIndices.size(); i++)
    {
        if (systemIndices[i] != params_.nnpIndices_[j])
        {
            params_.mmIndices_.push_back(systemIndices[i]);
        }
        else
        {
            if (j < params_.nnpIndices_.size() - 1)
            {
                j++;
            }
        }
    }
}

void NNPotOptions::setLocalInputAtomSet(const LocalAtomSet& localInputAtomSet)
{
    params_.nnpAtoms_ = std::make_unique<LocalAtomSet>(localInputAtomSet);
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

    params_.linkFrontier_ = preprocessNNPotTopology(
            top, params_.nnpIndices_, params_.embeddingScheme_, params_.nnpCharge_, *logger_, wi_);
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
            treeBuilder.addUniformArray<std::int64_t>(moduleName() + "-" + c_inputGroupTag_);
    for (const auto& indexValue : params_.nnpIndices_)
    {
        GroupIndexAdder.addValue(indexValue);
    }

    // Write MM atoms index
    GroupIndexAdder = treeBuilder.addUniformArray<std::int64_t>(moduleName() + "-" + c_mmGroupTag_);
    for (const auto& indexValue : params_.mmIndices_)
    {
        GroupIndexAdder.addValue(indexValue);
    }

    // Write link
    GroupIndexAdder = treeBuilder.addUniformArray<std::int64_t>(moduleName() + "-" + c_nnpLinkTag_);
    for (const auto& link : params_.linkFrontier_)
    {
        GroupIndexAdder.addValue(link.getEmbeddedIndex());
    }
    GroupIndexAdder = treeBuilder.addUniformArray<std::int64_t>(moduleName() + "-" + c_mmLinkTag_);
    for (const auto& link : params_.linkFrontier_)
    {
        GroupIndexAdder.addValue(link.getMMIndex());
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
    std::string key = moduleName() + "-" + c_inputGroupTag_;
    if (!tree.keyExists(key))
    {
        GMX_THROW(InconsistentInputError(
                "Cannot find input atoms index vector required for neural network potential.\n"
                "This could be caused by incompatible or corrupted tpr input file."));
    }
    auto kvtIndexArray = tree[key].asArray().values();
    params_.nnpIndices_.resize(kvtIndexArray.size());
    std::transform(std::begin(kvtIndexArray),
                   std::end(kvtIndexArray),
                   std::begin(params_.nnpIndices_),
                   [](const KeyValueTreeValue& val) { return val.cast<std::int64_t>(); });

    // Try to read MM atoms index
    key = moduleName() + "-" + c_mmGroupTag_;
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

    // Try to read Link Frontier (two separate vectors and then combine)
    std::vector<Index> nnpLink;
    std::vector<Index> mmLink;

    if (!tree.keyExists(moduleName() + "-" + c_nnpLinkTag_))
    {
        GMX_THROW(
                InconsistentInputError("Cannot find NNP Link Frontier vector required for QM/MM "
                                       "simulation.\nThis could be "
                                       "caused by incompatible or corrupted tpr input file."));
    }
    kvtIndexArray = tree[moduleName() + "-" + c_nnpLinkTag_].asArray().values();
    nnpLink.resize(kvtIndexArray.size());
    std::transform(std::begin(kvtIndexArray),
                   std::end(kvtIndexArray),
                   std::begin(nnpLink),
                   [](const KeyValueTreeValue& val) { return val.cast<std::int64_t>(); });

    if (!tree.keyExists(moduleName() + "-" + c_mmLinkTag_))
    {
        GMX_THROW(InconsistentInputError(
                "Cannot find MM Link Frontier vector required for QM/MM simulation.\nThis could be "
                "caused by incompatible or corrupted tpr input file."));
    }
    kvtIndexArray = tree[moduleName() + "-" + c_mmLinkTag_].asArray().values();
    mmLink.resize(kvtIndexArray.size());
    std::transform(std::begin(kvtIndexArray),
                   std::end(kvtIndexArray),
                   std::begin(mmLink),
                   [](const KeyValueTreeValue& val) { return val.cast<std::int64_t>(); });

    params_.linkFrontier_.reserve(nnpLink.size());
    for (size_t i = 0; i < nnpLink.size(); i++)
    {
        params_.linkFrontier_.emplace_back(nnpLink[i], mmLink[i]);
        params_.linkFrontier_.back().setLinkDistance(params_.linkDistance_);
        AtomProperties atomProp;
        const int      linkAtomNumber = atomProp.atomNumberFromElement(params_.linkType_.c_str());
        GMX_RELEASE_ASSERT(
                linkAtomNumber != -1,
                formatString("Unrecognized link atom type symbol: %s", params_.linkType_.c_str()).c_str());
        params_.linkFrontier_.back().setLinkAtomNumber(linkAtomNumber);
    }

    // add MM link atoms to nnp index
    for (const auto& link : params_.linkFrontier_)
    {
        params_.nnpIndices_.push_back(link.getMMIndex());
    }
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

void NNPotOptions::setComm(const MpiComm& mpiComm)
{
    mpiComm_ = &mpiComm;
}

const MDLogger& NNPotOptions::logger() const
{
    GMX_RELEASE_ASSERT(logger_, "Logger not set for NNPotOptions.");
    return *logger_;
}

const MpiComm& NNPotOptions::mpiComm() const
{
    GMX_RELEASE_ASSERT(mpiComm_, "MPI communicator not set for NNPotOptions.");
    return *mpiComm_;
}

void NNPotOptions::setWarninp(WarningHandler* wi)
{
    wi_ = wi;
}

void NNPotOptions::checkNNPotModel()
{
    if (std::getenv("GMX_NNPOT_SKIP_MODEL_CHECK"))
    {
        GMX_LOG(logger().info)
                .appendText("Skipping NNP model check because GMX_NNPOT_SKIP_MODEL_CHECK is set.");
        return;
    }
#if GMX_TORCH

    // try initializing the neural network model
    std::filesystem::path        modelPath(params_.modelFileName_);
    std::unique_ptr<INNPotModel> model;
    MpiComm                      comm(MpiComm::SingleRank{});
    if (!std::filesystem::exists(modelPath))
    {
        GMX_THROW(FileIOError("Model file does not exist: " + params_.modelFileName_));
    }
    else if (modelPath.extension() == ".pt")
    {
        model = std::make_unique<TorchModel>(
                params_.modelFileName_, params_.embeddingScheme_, logger(), comm);
    }
    else
    {
        GMX_THROW(FileIOError("Unrecognized extension for model file: " + params_.modelFileName_));
    }

    // check that link atom type is valid
    AtomProperties atomProp;
    const int      linkAtomNumber = atomProp.atomNumberFromElement(params_.linkType_.c_str());
    if (linkAtomNumber == -1)
    {
        GMX_THROW(InconsistentInputError(
                formatString("Unrecognized link atom type symbol: %s", params_.linkType_.c_str())));
    }

    // check if model accepts inputs
    // Prepare somewhat realistic dummy input
    const size_t     numNNPAtoms = params_.nnpIndices_.size();
    std::vector<int> indices(numNNPAtoms);
    std::iota(indices.begin(), indices.end(), 0);
    std::vector<RVec> positions(numNNPAtoms);
    // set positions to arbitrary values to avoid issues with certain models
    for (size_t i = 0; i < numNNPAtoms; i++)
    {
        // TODO: replace 3 here with DIM (or replacement) when macro conflicts are resolved
        for (size_t j = 0; j < 3; j++)
        {
            positions[i][j] = 0.01 * static_cast<real>(i + 1) * static_cast<real>(j + 1);
        }
    }
    std::vector<int> atomNumbers(
            numNNPAtoms, linkAtomNumber); // to check that the model can accept the link atom type
    std::vector<int> atomPairs;
    for (size_t i = 0; i < numNNPAtoms - 1; i++)
    {
        for (size_t j = i + 1; j < numNNPAtoms; j++)
        {
            atomPairs.insert(atomPairs.end(), { indices[i], indices[j] });
        }
    }
    std::vector<RVec>             pairShifts(atomPairs.size() / 2, RVec({ 0.0, 0.0, 0.0 }));
    std::vector<real>             mmCharges(1, 1.0);
    std::vector<int>              mmIndices(1, 1);
    matrix                        box = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };
    PbcType                       pbc = PbcType();
    std::vector<LinkFrontierAtom> link;

    // check that inputs are not empty
    if (params_.modelInput_.size() == 0)
    {
        GMX_THROW(InconsistentInputError("No inputs to NN model provided."));
    }
    // check that sensible cutoff was provided
    const bool pairlistNeeded =
            params_.modelNeedsInput("atom-pairs") || params_.modelNeedsInput("pair-shifts");
    if (pairlistNeeded && (params_.pairCutoff_ <= 0.0))
    {
        GMX_THROW(InconsistentInputError(formatString(
                "List of atom pairs was requested as input to the NNP model, but pair-cutoff is "
                "%.1f. Please specify a valid cutoff radius or disable pair input.",
                params_.pairCutoff_)));
    }

    const auto& theLogger = logger();
    // Check that warning handler is valid
    GMX_ASSERT(wi_, "WarningHandler not set.");
    bool modelOutputsForces = false;
    try
    {
        // prepare dummy output (NNP Atoms and 1 MM atom)
        gmx_enerdata_t    enerd(numNNPAtoms + 1, nullptr);
        std::vector<RVec> forcesVec(numNNPAtoms + 1, RVec({ 0.0, 0.0, 0.0 }));
        ArrayRef<RVec>    forces(forcesVec);

        // might throw a runtime error if the model is not compatible with the dummy input
        model->evaluateModel(&enerd,
                             forces,
                             indices,
                             mmIndices,
                             params_.modelInput_,
                             positions,
                             atomNumbers,
                             atomPairs,
                             pairShifts,
                             positions,
                             mmCharges,
                             params_.nnpCharge_,
                             link,
                             box,
                             pbc);

        // check if model outputs forces after forward pass
        modelOutputsForces = model->outputsForces();

        // log wheter model outputs forces or we need to compute them
        if (modelOutputsForces)
        {
            GMX_LOG(theLogger.info).appendText("Will use forces from NNP model.");
        }
        else
        {
            if (params_.embeddingScheme_ == NNPotEmbedding::ElectrostaticModel)
            {
                GMX_THROW(
                        InconsistentInputError("NNP model does not output forces, but they are "
                                               "expected to be computed by the model for "
                                               "electrostatic embedding scheme."));
            }
            GMX_LOG(theLogger.info)
                    .appendText(
                            "NNP model does not output forces. They will be computed as gradients "
                            "of the energy w.r.t. the first input tensor (atom positions).");
        }
    }
    catch (const std::exception& e)
    {
        GMX_THROW(
                InconsistentInputError(
                        "There was an error while checking NN model with a dummy input: "
                        + std::string(e.what())
                        + "\n"
                          "Can't verify that the model works correctly. Make sure that the NNP "
                          "model and inputs are compatible."));
    }

    // check if first input is atom-positions
    if ((params_.modelInput_[0] != "atom-positions") && !modelOutputsForces)
    {
        wi_->addWarning("Gradients will be computed with respect to first input to NN model "
                        + params_.modelInput_[0] + " instead of atom positions. Is this intended?");
    }

    // check if model might need MM atom positions and/or charges
    if (params_.embeddingScheme_ == NNPotEmbedding::ElectrostaticModel)
    {
        if (!params_.modelNeedsInput("atom-positions-mm") && !params_.modelNeedsInput("atom-charges-mm"))
        {
            wi_->addWarning(
                    "Embedding scheme is set to Electrostatic model, but MM positions and/or "
                    "charges are not requested as model input. Is this intended?");
        }
    }

#else
    GMX_THROW(InternalError(
            "Libtorch/NN backend is not linked into GROMACS, NNPot simulation is not possible."
            " Please, reconfigure GROMACS with -DGMX_NNPOT=TORCH\n"));
#endif // GMX_TORCH
}

} // namespace gmx
