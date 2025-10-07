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

} // namespace

void RAMDOptions::initMdpTransform(IKeyValueTreeTransformRules* rules)
{
    addMdpTransformFromString<bool>(
        rules, &fromStdString<bool>, RAMDModuleInfo::sc_name, c_activeTag);
    addMdpTransformFromString<std::int64_t>(
        rules, &fromStdString<std::int64_t>, RAMDModuleInfo::sc_name, c_seedTag);
    addMdpTransformFromString<int>(
        rules, &fromStdString<int>, RAMDModuleInfo::sc_name, c_evalFreqTag);
}

void RAMDOptions::buildMdpOutput(KeyValueTreeObjectBuilder* builder) const
{
    addMdpOutputComment(builder, RAMDModuleInfo::sc_name, "empty-line", "");
    addMdpOutputComment(
        builder, RAMDModuleInfo::sc_name, "module", "; Density guided simulation");
    addMdpOutputValue(builder, RAMDModuleInfo::sc_name, c_activeTag, parameters_.active_);

    if (parameters_.active_)
    {
        // Seed
        addMdpOutputComment(
                builder, RAMDModuleInfo::sc_name, c_seedTag, "; Seed");
        addMdpOutputValue(builder, RAMDModuleInfo::sc_name, c_seedTag, parameters_.seed_);
    }
}

void RAMDOptions::initMdpOptions(IOptionsContainerWithSections* options)
{
    auto section = options->addSection(OptionSection(moduleName().c_str()));
    section.addOption(BooleanOption(c_activeTag.c_str()).store(&parameters_.active_));
    section.addOption(Int64Option(c_seedTag.c_str()).store(&parameters_.seed_));
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

    // Create input index
    parameters_.groups_[0].receptor_indices_ = indexGroupsAndNames.indices(parameters_.groups_[0].receptor_group_);

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
