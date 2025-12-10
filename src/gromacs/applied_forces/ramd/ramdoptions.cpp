/*! \internal \file
 * \brief
 * Implements force provider for RAMD
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "ramdoptions.h"

#include "gromacs/applied_forces/ramd/ramd.h"
#include "gromacs/math/densityfit.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/imdpoptionprovider_helpers.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/textreader.h"

namespace gmx
{

namespace
{

//! Helper function to make a std::string containing the module name
std::string moduleName()
{
    return std::string(RAMDModuleInfo::sc_name);
}

const std::string c_activeTag = "active";
const std::string c_seedTag = "seed";
const std::string c_evalFreqTag = "eval-freq";
const std::string c_outFreqTag = "out-freq";
const std::string c_groupsFileTag = "groups-file";
const std::string c_groupsStringTag = "groups-string";
const std::string c_groupReceptorTag = "receptor";
const std::string c_groupReceptorIndicesTag = "receptor-indices";
const std::string c_groupLigandTag = "ligand";
const std::string c_groupLigandIndicesTag = "ligand-indices";
const std::string c_groupForceTag = "force";
const std::string c_groupMaxDistTag = "max-dist";
const std::string c_groupRMinDistTag = "r-min-dist";
const std::string c_pbcRefPrevStepComTag = "pbc-ref-prev-step-com";
const std::string c_oldAngleDistTag = "old-angle-dist";
const std::string c_connectedLigandsTag = "connected-ligands";

} // namespace

void RAMDOptions::initMdpTransform(IKeyValueTreeTransformRules* rules)
{
    const auto& stringIdentityTransform = [](std::string s) { return s; };
    addMdpTransformFromString<bool>(rules, &fromStdString<bool>, RAMDModuleInfo::sc_name, c_activeTag);
    addMdpTransformFromString<std::int64_t>(rules, &fromStdString<std::int64_t>, RAMDModuleInfo::sc_name, c_seedTag);
    addMdpTransformFromString<int>(rules, &fromStdString<int>, RAMDModuleInfo::sc_name, c_evalFreqTag);
    addMdpTransformFromString<int>(rules, &fromStdString<int>, RAMDModuleInfo::sc_name, c_outFreqTag);
    addMdpTransformFromString<std::string>(rules, stringIdentityTransform, RAMDModuleInfo::sc_name, c_groupsFileTag);
    addMdpTransformFromString<bool>(rules, &fromStdString<bool>, RAMDModuleInfo::sc_name, c_pbcRefPrevStepComTag);
    addMdpTransformFromString<bool>(rules, &fromStdString<bool>, RAMDModuleInfo::sc_name, c_oldAngleDistTag);
    addMdpTransformFromString<bool>(rules, &fromStdString<bool>, RAMDModuleInfo::sc_name, c_connectedLigandsTag);
}

void RAMDOptions::buildMdpOutput(KeyValueTreeObjectBuilder* builder) const
{
    addMdpOutputComment(builder, RAMDModuleInfo::sc_name, "empty-line", "");
    addMdpOutputComment(builder, RAMDModuleInfo::sc_name, "module", "; RAMD");
    addMdpOutputValue(builder, RAMDModuleInfo::sc_name, c_activeTag, parameters_.active_);

    if (parameters_.active_)
    {
        addMdpOutputValue(builder, RAMDModuleInfo::sc_name, c_seedTag, parameters_.seed_);
        addMdpOutputValue(builder, RAMDModuleInfo::sc_name, c_evalFreqTag, parameters_.eval_freq_);
        addMdpOutputValue(builder, RAMDModuleInfo::sc_name, c_outFreqTag, parameters_.out_freq_);
        addMdpOutputValue(builder, RAMDModuleInfo::sc_name, c_groupsFileTag, groupsFile_);
        addMdpOutputValue(
                builder, RAMDModuleInfo::sc_name, c_pbcRefPrevStepComTag, parameters_.pbc_ref_prev_step_com_);
        addMdpOutputValue(
                builder, RAMDModuleInfo::sc_name, c_connectedLigandsTag, parameters_.connected_ligands_);
    }
}

void RAMDOptions::initMdpOptions(IOptionsContainerWithSections* options)
{
    auto section = options->addSection(OptionSection(moduleName().c_str()));
    section.addOption(BooleanOption(c_activeTag.c_str()).store(&parameters_.active_));
    section.addOption(Int64Option(c_seedTag.c_str()).store(&parameters_.seed_));
    section.addOption(IntegerOption(c_evalFreqTag.c_str()).store(&parameters_.eval_freq_));
    section.addOption(IntegerOption(c_outFreqTag.c_str()).store(&parameters_.out_freq_));
    section.addOption(StringOption(c_groupsFileTag.c_str()).store(&groupsFile_));
    section.addOption(BooleanOption(c_pbcRefPrevStepComTag.c_str()).store(&parameters_.pbc_ref_prev_step_com_));
    section.addOption(BooleanOption(c_oldAngleDistTag.c_str()).store(&parameters_.old_angle_dist_));
    section.addOption(BooleanOption(c_connectedLigandsTag.c_str()).store(&parameters_.connected_ligands_));
}

bool RAMDOptions::active() const
{
    return parameters_.active_;
}

const RAMDParameters& RAMDOptions::parameters()
{
    return parameters_;
}

void RAMDOptions::setLogger(const MDLogger& logger)
{
    // Exit if RAMD module is not active
    if (!parameters_.active_)
    {
        return;
    }

    logger_ = &logger;
}

const MDLogger& RAMDOptions::logger() const
{
    GMX_RELEASE_ASSERT(logger_, "Logger not set for RAMDOptions.");
    return *logger_;
}

void RAMDOptions::setInputGroupIndices(const IndexGroupsAndNames& indexGroupsAndNames)
{
    // Exit if RAMD module is not active
    if (!parameters_.active_)
    {
        return;
    }

    // Copy the content of the RAMD input file into a string for latter save in KVT
    if (!groupsFile_.empty())
    {
        groupsString_ = TextReader::readFileToString(groupsFile_);
        readConfigString();
    }

    // Create input index
    for (int g = 0; g < parameters_.ngroups_; ++g)
    {
        parameters_.groups_[g].ligand_indices_ = indexGroupsAndNames.indices(parameters_.groups_[g].ligand_);
        parameters_.groups_[g].receptor_indices_ = indexGroupsAndNames.indices(parameters_.groups_[g].receptor_);

        // Check that group is not empty
        if (parameters_.groups_[g].ligand_indices_.empty())
        {
            GMX_THROW(InconsistentInputError(
                formatString("Group %s defining RAMD ligand atoms should not be empty.",
                             parameters_.groups_[g].ligand_.c_str())));
        }
        if (parameters_.groups_[g].receptor_indices_.empty())
        {
            GMX_THROW(InconsistentInputError(
                formatString("Group %s defining RAMD receptor atoms should not be empty.",
                             parameters_.groups_[g].receptor_.c_str())));
        }
    }
}

void RAMDOptions::writeInternalParametersToKvt(KeyValueTreeObjectBuilder treeBuilder)
{
    // Write groups input file as a string
    treeBuilder.addValue<std::string>(moduleName() + "-" + c_groupsStringTag, groupsString_);
    for (size_t g = 0; g < parameters_.groups_.size(); ++g)
    {
        std::string groupName = moduleName() + "-group-" + std::to_string(g);
        auto arrayBuilder = treeBuilder.addUniformArray<Index>(groupName + "-" + c_groupReceptorIndicesTag);
        for (const auto& val : parameters_.groups_[g].receptor_indices_)
        {
            arrayBuilder.addValue(val);
        }
        auto arrayBuilder2 = treeBuilder.addUniformArray<Index>(groupName + "-" + c_groupLigandIndicesTag);
        for (const auto& val : parameters_.groups_[g].ligand_indices_)
        {
            arrayBuilder2.addValue(val);
        }
    }
}

void RAMDOptions::readInternalParametersFromKvt(const KeyValueTreeObject& tree)
{
    // Check if active
    if (!parameters_.active_) return;

    if (!tree.keyExists(moduleName() + "-" + c_groupsStringTag))
    {
        GMX_THROW(InconsistentInputError(
                "Cannot find ramd-groups-string required for colvars simulation."));
    }

    groupsString_ = tree[moduleName() + "-" + c_groupsStringTag].cast<std::string>();
    readConfigString();

    for (size_t g = 0; g < parameters_.groups_.size(); ++g)
    {
        std::string groupName = moduleName() + "-group-" + std::to_string(g);        
        parameters_.groups_[g].receptor_indices_.resize(
            tree[groupName + "-" + c_groupReceptorIndicesTag].asArray().values().size());
        for (size_t i = 0; i < parameters_.groups_[g].receptor_indices_.size(); ++i)
        {
            parameters_.groups_[g].receptor_indices_[i] =
                tree[groupName + "-" + c_groupReceptorIndicesTag].asArray().values()[i].cast<Index>();
        }
        parameters_.groups_[g].ligand_indices_.resize(
            tree[groupName + "-" + c_groupLigandIndicesTag].asArray().values().size());
        for (size_t i = 0; i < parameters_.groups_[g].ligand_indices_.size(); ++i)
        {
            parameters_.groups_[g].ligand_indices_[i] =
                tree[groupName + "-" + c_groupLigandIndicesTag].asArray().values()[i].cast<Index>();
        }
    }
}

void RAMDOptions::readConfigString()
{
    std::string line, key, value;
    std::istringstream ss(groupsString_);

    while (std::getline(ss, line)) {
        if (line.find_first_not_of(" \t") == std::string::npos) continue;
        if (line.find_first_of("#;") != std::string::npos) continue;
        std::istringstream lineStream(line);
        lineStream >> key;
        if (key == "ramd-group")
        {
            gmx::RAMDGroup newGroup;
            while (std::getline(ss, line)) {
                if (line.find_first_not_of(" \t") == std::string::npos) continue;
                if (line.find_first_of("#;") != std::string::npos) continue;
                std::istringstream lineStream(line);
                lineStream >> key;
                if (key == "}") {
                    parameters_.groups_.push_back(newGroup);
                    break;
                }
                lineStream >> value;
                if (key == "receptor")
                {
                    newGroup.receptor_ = value;
                }
                if (key == "receptor-pbcatom")
                {
                    newGroup.receptor_pbcatom_ = std::stoi(value);
                }
                if (key == "ligand")
                {
                    newGroup.ligand_ = value;
                }
                if (key == "ligand-pbcatom")
                {
                    newGroup.ligand_pbcatom_ = std::stoi(value);
                }
                if (key == "force")
                {
                    newGroup.force_ = std::stod(value);
                }
                if (key == "max-dist")
                {
                    newGroup.max_dist_ = std::stod(value);
                }
                if (key == "r-min-dist")
                {
                    newGroup.r_min_dist_ = std::stod(value);
                }
            }
        }
    }

    // Write number of groups for reading KeyValueTreeObject
    parameters_.ngroups_ = static_cast<int>(parameters_.groups_.size());
}

} // namespace gmx
