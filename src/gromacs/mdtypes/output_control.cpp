/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * Defines data structure and utilities for output control parameters
 *
 * Eventually the declaration of OutputControl will be
 * found here, once t_inputrec is fully reformed.
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "output_control.h"

#include <memory>

#include "gromacs/gmxpreprocess/inputrecstrings.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdoutputprovider.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainerwithsections.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/strconvert.h"

struct gmx_output_env_t;

namespace gmx
{
struct MDModulesNotifiers;

namespace
{

/*! \internal
 * \brief Manages output control parameters
 *
 * Handles reading and writing of output frequency parameters through
 * the key-value tree system. This module does not own the OutputControl
 * data; instead it writes to an external OutputControl object (typically
 * embedded in t_inputrec) set via setOutputControl().
 */
class OutputControlModule final : public IMDModule, public IMdpOptionProvider
{
public:
    OutputControlModule() = default;

    // From IMDModule
    IMdpOptionProvider* mdpOptionProvider() override { return this; }
    IMDOutputProvider*  outputProvider() override { return nullptr; }
    void                initForceProviders(ForceProviders* /* forceProviders */) override {}

    // From IMdpOptionProvider
    void initMdpTransform(IKeyValueTreeTransformRules* transform) override;
    void initMdpOptions(IOptionsContainerWithSections* options) override;
    void buildMdpOutput(KeyValueTreeObjectBuilder* builder) const override;

    void subscribeToSimulationSetupNotifications(MDModulesNotifiers* /* notifiers */) override {}
    void subscribeToSimulationRunNotifications(MDModulesNotifiers* /* notifiers */) override {}
    void subscribeToPreProcessingNotifications(MDModulesNotifiers* /* notifiers */) override {}

    /*! \brief Set the target OutputControl to write to
     *
     * This must be called before initMdpOptions is used, typically after
     * creating t_inputrec and before parsing MDP options.
     *
     * \param[in] outputControl Pointer to the OutputControl to write to
     */
    void setOutputControl(OutputControl* outputControl) { outputControl_ = outputControl; }

    /*! \brief Enable MDP output writing (for testing)
     *
     * By default, buildMdpOutput() does nothing to avoid duplication with
     * the inp vector mechanism. This method enables output writing so tests
     * can verify the module produces correct output.
     *
     * \param[in] enable Whether to enable MDP output writing
     */
    void setWriteToMdpOutput(bool enable) { writeToMdpOutput_ = enable; }

    /*! \brief Set storage for preprocessing-only string parameters
     *
     * Group specifications like compressed-x-grps and energygrps are only
     * used during grompp preprocessing. This method sets where to store
     * those strings (separate from the runtime OutputControl struct).
     *
     * REQUIRED: Must be called during grompp MDP parsing before assignOptionsToModules().
     * NOT NEEDED: Not required when reading TPR files in mdrun/tools.
     *
     * \param[in] preprocessingStrings Pointer to storage for preprocessing-only strings (must not be null)
     */
    void setPreprocessingStrings(gmx_inputrec_strings* preprocessingStrings)
    {
        GMX_RELEASE_ASSERT(preprocessingStrings,
                           "setPreprocessingStrings() called with null pointer");
        preprocessingGroupNames_ = preprocessingStrings;
    }

private:
    //! Pointer to external OutputControl (typically embedded in t_inputrec)
    OutputControl* outputControl_ = nullptr;
    //! Pointer to preprocessing-only string storage (local to grompp)
    gmx_inputrec_strings* preprocessingGroupNames_ = nullptr;
    //! Whether to write output in buildMdpOutput (default false, for testing only)
    bool writeToMdpOutput_ = false;
};

void OutputControlModule::initMdpTransform(IKeyValueTreeTransformRules* rules)
{
    // Transform flat MDP keys to hierarchical structure
    // MDP files are text, so values are initially strings that need conversion
    rules->addRule().from<std::string>("/nstlog").to<int>("/output-control/nstlog").transformWith(&fromStdString<int>);
    rules->addRule().from<std::string>("/nstxout").to<int>("/output-control/nstxout").transformWith(&fromStdString<int>);
    rules->addRule().from<std::string>("/nstvout").to<int>("/output-control/nstvout").transformWith(&fromStdString<int>);
    rules->addRule().from<std::string>("/nstfout").to<int>("/output-control/nstfout").transformWith(&fromStdString<int>);
    rules->addRule().from<std::string>("/nstenergy").to<int>("/output-control/nstenergy").transformWith(&fromStdString<int>);
    rules->addRule()
            .from<std::string>("/nstxout-compressed")
            .to<int>("/output-control/nstxout-compressed")
            .transformWith(&fromStdString<int>);
    rules->addRule()
            .from<std::string>("/compressed-x-precision")
            .to<real>("/output-control/x-compression-precision")
            .transformWith(&fromStdString<real>);
    rules->addRule()
            .from<std::string>("/nstcalcenergy")
            .to<int>("/output-control/nstcalcenergy")
            .transformWith(&fromStdString<int>);

    // Preprocessing-only string parameters (identity transform for strings)
    rules->addRule()
            .from<std::string>("/compressed-x-grps")
            .to<std::string>("/output-control/compressed-x-groups")
            .transformWith([](const std::string& value) { return value; });
    rules->addRule()
            .from<std::string>("/energygrps")
            .to<std::string>("/output-control/energy-groups")
            .transformWith([](const std::string& value) { return value; });
}

void OutputControlModule::initMdpOptions(IOptionsContainerWithSections* options)
{
    GMX_RELEASE_ASSERT(outputControl_, "OutputControl must be set before initMdpOptions");
    auto section = options->addSection(OptionSection("output-control"));
    section.addOption(IntegerOption("nstlog").store(&outputControl_->nstlog));
    section.addOption(IntegerOption("nstxout").store(&outputControl_->nstxout));
    section.addOption(IntegerOption("nstvout").store(&outputControl_->nstvout));
    section.addOption(IntegerOption("nstfout").store(&outputControl_->nstfout));
    section.addOption(IntegerOption("nstenergy").store(&outputControl_->nstenergy));
    section.addOption(IntegerOption("nstxout-compressed").store(&outputControl_->nstxout_compressed));
    section.addOption(RealOption("x-compression-precision").store(&outputControl_->x_compression_precision));
    section.addOption(IntegerOption("nstcalcenergy").store(&outputControl_->nstcalcenergy));

    // Preprocessing-only string options (store in preprocessingGroupNames_)
    // These are only registered when preprocessingGroupNames_ is set (grompp path).
    // For mdrun/tools reading TPR, these options are not registered since the
    // preprocessing has already happened and the group indices are stored.
    if (preprocessingGroupNames_)
    {
        section.addOption(StringOption("compressed-x-groups")
                                  .store(&preprocessingGroupNames_->compressedXGroups)
                                  .description("Atom groups for compressed trajectory output"));
        section.addOption(StringOption("energy-groups")
                                  .store(&preprocessingGroupNames_->energyGroups)
                                  .description("Atom groups for energy calculation"));
    }
}

void OutputControlModule::buildMdpOutput(KeyValueTreeObjectBuilder* builder) const
{
    // By default, do nothing. Output-control fields are written via the inp vector mechanism
    // to maintain their historical position between TEST PARTICLE INSERTION and NEIGHBORSEARCHING.
    if (!writeToMdpOutput_)
    {
        return;
    }

    // When writeToMdpOutput_ is enabled (for testing), we write the fields here to verify
    // the module can correctly represent its state in MDP output format. Then when all
    // mdp options are handled with KVT processing we can easily remove the legacy processing
    // while keeping the format of the output mdp file stable.

    // Test-only path: preprocessingStrings must be set
    GMX_RELEASE_ASSERT(
            preprocessingGroupNames_,
            "buildMdpOutput enabled for testing but preprocessingStrings not set. "
            "Tests must call setOutputControlPreprocessingStrings() before enabling MDP output.");

    GMX_RELEASE_ASSERT(outputControl_, "OutputControl not set");

    // Section header comment
    builder->addValue<std::string>("comment-output-control-section", "\n; OUTPUT CONTROL OPTIONS");

    // Coordinate/velocity/force output comment and fields
    builder->addValue<std::string>(
            "comment-output-control-coords",
            "; Output frequency for coords (x), velocities (v) and forces (f)");
    builder->addValue<int>("nstxout", outputControl_->nstxout);
    builder->addValue<int>("nstvout", outputControl_->nstvout);
    builder->addValue<int>("nstfout", outputControl_->nstfout);

    // Energy output comment and fields
    builder->addValue<std::string>("comment-output-control-energy",
                                   "; Output frequency for energies to log file and energy file");
    builder->addValue<int>("nstlog", outputControl_->nstlog);
    builder->addValue<int>("nstcalcenergy", outputControl_->nstcalcenergy);
    builder->addValue<int>("nstenergy", outputControl_->nstenergy);

    // Compressed trajectory comment and fields
    builder->addValue<std::string>("comment-output-control-compressed",
                                   "; Output frequency and precision for .xtc file");
    builder->addValue<int>("nstxout-compressed", outputControl_->nstxout_compressed);
    builder->addValue<real>("compressed-x-precision", outputControl_->x_compression_precision);

    // Group specification comments and fields (accessed directly from preprocessingGroupNames_)
    // In production, these are also written via the legacy inp vector mechanism for proper ordering.
    // When writeToMdpOutput_ is enabled (testing only), they are written here as well.
    // Assertion above guarantees preprocessingGroupNames_ is non-null in this test-only path.
    builder->addValue<std::string>("comment-output-control-compressed-groups",
                                   "; This selects the subset of atoms for the compressed\n"
                                   "; trajectory file. You can select multiple groups. By\n"
                                   "; default, all atoms will be written.");
    builder->addValue<std::string>("compressed-x-grps", preprocessingGroupNames_->compressedXGroups);

    builder->addValue<std::string>("comment-output-control-energy-groups",
                                   "; Selection of energy groups");
    builder->addValue<std::string>("energygrps", preprocessingGroupNames_->energyGroups);
}

} // namespace

std::unique_ptr<IMDModule> OutputControlModuleInfo::create()
{
    return std::make_unique<OutputControlModule>();
}

void setOutputControlTarget(IMDModule* module, OutputControl* outputControl)
{
    auto* outputControlModule = dynamic_cast<OutputControlModule*>(module);
    GMX_RELEASE_ASSERT(outputControlModule, "Module must be an OutputControlModule");
    outputControlModule->setOutputControl(outputControl);
}

void setOutputControlWriteToMdpOutput(IMDModule* module, bool enable)
{
    auto* outputControlModule = dynamic_cast<OutputControlModule*>(module);
    GMX_RELEASE_ASSERT(outputControlModule, "Module must be an OutputControlModule");
    outputControlModule->setWriteToMdpOutput(enable);
}

void setOutputControlPreprocessingStrings(IMDModule* module, gmx_inputrec_strings* preprocessingStrings)
{
    auto* outputControlModule = dynamic_cast<OutputControlModule*>(module);
    GMX_RELEASE_ASSERT(outputControlModule, "Module must be an OutputControlModule");
    outputControlModule->setPreprocessingStrings(preprocessingStrings);
}

} // namespace gmx
