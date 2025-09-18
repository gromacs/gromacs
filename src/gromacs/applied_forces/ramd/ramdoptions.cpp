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

} // namespace

void RAMDOptions::initMdpTransform(IKeyValueTreeTransformRules* rules)
{

}

void RAMDOptions::buildMdpOutput(KeyValueTreeObjectBuilder* builder) const
{
    addMdpOutputComment(builder, RAMDModuleInfo::sc_name, "empty-line", "");
    addMdpOutputComment(
            builder, RAMDModuleInfo::sc_name, "module", "; Density guided simulation");
    addMdpOutputValue(builder, RAMDModuleInfo::sc_name, c_activeTag, parameters_.active_);
}

void RAMDOptions::initMdpOptions(IOptionsContainerWithSections* options)
{
    auto section = options->addSection(OptionSection(moduleName().c_str()));
    section.addOption(BooleanOption(c_activeTag.c_str()).store(&parameters_.active_));
}

bool RAMDOptions::active() const
{
    return parameters_.active_;
}

} // namespace gmx
