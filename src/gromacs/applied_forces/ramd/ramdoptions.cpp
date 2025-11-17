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
const std::string c_groupsFileTag = "groups-file";
const std::string c_groupsStringTag = "groups-string";


// const std::string c_ngroupsTag = "ngroups";
// const std::string c_groupReceptorTag = "group1-receptor";
// const std::string c_groupLigandTag = "group1-ligand";
// const std::string c_groupForceTag = "group1-force";
// const std::string c_groupMaxDistTag = "group1-max-dist";
// const std::string c_groupRMinDistTag = "group1-r-min-dist";

} // namespace

void RAMDOptions::initMdpTransform(IKeyValueTreeTransformRules* rules)
{
    const auto& stringIdentityTransform = [](std::string s) { return s; };
    addMdpTransformFromString<bool>(rules, &fromStdString<bool>, RAMDModuleInfo::sc_name, c_activeTag);
    addMdpTransformFromString<std::int64_t>(rules, &fromStdString<std::int64_t>, RAMDModuleInfo::sc_name, c_seedTag);
    addMdpTransformFromString<int>(rules, &fromStdString<int>, RAMDModuleInfo::sc_name, c_evalFreqTag);
    addMdpTransformFromString<std::string>(rules, stringIdentityTransform, RAMDModuleInfo::sc_name, c_groupsFileTag);
}

void RAMDOptions::buildMdpOutput(KeyValueTreeObjectBuilder* builder) const
{
    addMdpOutputComment(builder, RAMDModuleInfo::sc_name, "empty-line", "");
    addMdpOutputComment(builder, RAMDModuleInfo::sc_name, "module", "; RAMD options");
    addMdpOutputValue(builder, RAMDModuleInfo::sc_name, c_activeTag, parameters_.active_);

    if (parameters_.active_)
    {
        addMdpOutputValue(builder, RAMDModuleInfo::sc_name, c_seedTag, parameters_.seed_);
        addMdpOutputValue(builder, RAMDModuleInfo::sc_name, c_evalFreqTag, parameters_.eval_freq_);
        addMdpOutputValue(builder, RAMDModuleInfo::sc_name, c_groupsFileTag, groupsFile_);
    }
}

void RAMDOptions::initMdpOptions(IOptionsContainerWithSections* options)
{
    auto section = options->addSection(OptionSection(moduleName().c_str()));
    section.addOption(BooleanOption(c_activeTag.c_str()).store(&parameters_.active_));
    section.addOption(Int64Option(c_seedTag.c_str()).store(&parameters_.seed_));
    section.addOption(IntegerOption(c_evalFreqTag.c_str()).store(&parameters_.eval_freq_));
    section.addOption(StringOption(c_groupsFileTag.c_str()).store(&groupsFile_));
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

void RAMDOptions::writeInternalParametersToKvt(KeyValueTreeObjectBuilder treeBuilder)
{
    // Copy the content of the colvars input file into a string for latter save in KVT
    if (!groupsFile_.empty())
    {
        groupsString_ = TextReader::readFileToString(groupsFile_);
    }

    // Write groups input file as a string
    treeBuilder.addValue<std::string>(moduleName() + "-" + c_groupsStringTag, groupsString_);

    //TODO: read groups from file
    // For now, just create one dummy group
    parameters_.ngroups_ = 1;
    parameters_.groups_.resize(parameters_.ngroups_);
    for (int g = 0; g < parameters_.ngroups_; ++g)
    {
        parameters_.groups_[g].ligand_group_ = "1SOL";
        parameters_.groups_[g].receptor_group_ = "2SOL";
    }
}

void RAMDOptions::readInternalParametersFromKvt(const KeyValueTreeObject& tree)
{
    // Check if active
    if (!parameters_.active_)
    {
        return;
    }

    if (!tree.keyExists(moduleName() + "-" + c_groupsStringTag))
    {
        GMX_THROW(InconsistentInputError(
                "Cannot find ramd-groups-string required for colvars simulation."));
    }
    groupsString_ = tree[moduleName() + "-" + c_groupsStringTag].cast<std::string>();

    parameters_.ngroups_ = 1;
    parameters_.groups_.resize(parameters_.ngroups_);
    for (int g = 0; g < parameters_.ngroups_; ++g)
    {
        parameters_.groups_[g].ligand_group_ = "1SOL";
        parameters_.groups_[g].ligand_indices_ = {1, 2, 3};
        parameters_.groups_[g].receptor_group_ = "2SOL";
        parameters_.groups_[g].receptor_indices_ = {4, 5, 6};
        parameters_.groups_[g].force_ = 100.0;
        parameters_.groups_[g].max_dist_ = 1.0;
        parameters_.groups_[g].r_min_dist_ = 0.0025;
    }
}

void RAMDOptions::setInputGroupIndices(const IndexGroupsAndNames& indexGroupsAndNames)
{
    // Exit if RAMD module is not active
    if (!parameters_.active_)
    {
        return;
    }

    // Create input index
    for (int g = 0; g < parameters_.ngroups_; ++g)
    {
        parameters_.groups_[g].ligand_indices_ = indexGroupsAndNames.indices(parameters_.groups_[g].ligand_group_);
        parameters_.groups_[g].receptor_indices_ = indexGroupsAndNames.indices(parameters_.groups_[g].receptor_group_);

        // Check that group is not empty
        if (parameters_.groups_[g].ligand_indices_.empty())
        {
            GMX_THROW(InconsistentInputError(
                formatString("Group %s defining RAMD ligand atoms should not be empty.",
                             parameters_.groups_[g].ligand_group_.c_str())));
        }
        if (parameters_.groups_[g].receptor_indices_.empty())
        {
            GMX_THROW(InconsistentInputError(
                formatString("Group %s defining RAMD receptor atoms should not be empty.",
                             parameters_.groups_[g].receptor_group_.c_str())));
        }
    }

    // // Check that group is not empty
    // if (params_.nnpIndices_.empty())
    // {
    //     GMX_THROW(InconsistentInputError(
    //             formatString("Group %s defining NN potential input atoms should not be empty.",
    //                          params_.inputGroup_.c_str())));
    // }

    // // Create temporary index for the whole System
    // auto systemIndices = indexGroupsAndNames.indices("System");

    // // Sort nnpIndices_ and sysIndices_
    // std::sort(params_.nnpIndices_.begin(), params_.nnpIndices_.end());
    // std::sort(systemIndices.begin(), systemIndices.end());

    // // Create MM index
    // params_.mmIndices_.reserve(systemIndices.size());

    // // Position in nnpIndices_
    // size_t j = 0;
    // // Write to mmIndices_ only the atoms which do not belong to NNP input region
    // for (size_t i = 0; i < systemIndices.size(); i++)
    // {
    //     if (systemIndices[i] != params_.nnpIndices_[j])
    //     {
    //         params_.mmIndices_.push_back(systemIndices[i]);
    //     }
    //     else
    //     {
    //         if (j < params_.nnpIndices_.size() - 1)
    //         {
    //             j++;
    //         }
    //     }
    // }
}

} // namespace gmx
