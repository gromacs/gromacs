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

} // namespace gmx
