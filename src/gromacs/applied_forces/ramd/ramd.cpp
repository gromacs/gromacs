/*! \internal \file
 * \brief
 * Implements RAMD class that implements IMDModule interface
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 * \ingroup module_applied_forces
 */

#include "ramd.h"
#include "ramdoptions.h"
#include "ramdoutputprovider.h"
#include "ramdforceprovider.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"

namespace gmx
{

namespace
{

/*! \internal
 * \brief Helper class that holds simulation data and
 * callback functions for simulation setup time notifications
 */
class RAMDSimulationParameterSetup
{
public:
    RAMDSimulationParameterSetup() = default;

    /*! \brief Set the periodic boundary condition via MdModuleNotifier.
     *
     * The pbc type is wrapped in PeriodicBoundaryConditionType to
     * allow the MdModuleNotifier to statically distinguish the callback
     * function type from other 'int' function callbacks.
     *
     * \param[in] pbc MdModuleNotification class that contains a variable
     *                that enumerates the periodic boundary condition.
     */
    void setPeriodicBoundaryConditionType(const PbcType& pbc)
    {
        pbcType_ = std::make_unique<PbcType>(pbc);
    }

    //! Get the periodic boundary conditions
    PbcType periodicBoundaryConditionType()
    {
        if (pbcType_ == nullptr)
        {
            GMX_THROW(InternalError("Periodic boundary condition enum not set for RAMD."));
        }
        return *pbcType_;
    }

    /*! \brief Set the logger for RAMD during mdrun
     * \param[in] logger Logger instance to be used for output
     */
    void setLogger(const MDLogger& logger) { logger_ = &logger; }

    //! Get the logger instance
    const MDLogger& logger() const
    {
        GMX_RELEASE_ASSERT(logger_, "Logger not set for RAMD.");
        return *logger_;
    }

private:
    //! The type of periodic boundary conditions in the simulation
    std::unique_ptr<PbcType> pbcType_;
    /*! \brief MDLogger during mdrun
     *
     * This is a pointer only because we need an "optional reference"
     * to a const MDLogger before the notification always provides the
     * actual reference. */
    const MDLogger* logger_ = nullptr;

    GMX_DISALLOW_COPY_AND_ASSIGN(RAMDSimulationParameterSetup);
};


/*! \internal
 * \brief RAMD module
 */
class RAMD final : public IMDModule, public RAMDOutputProvider
{
public:
    //! \brief Construct the RAMD module.
    explicit RAMD() = default;

    //! Returns an interface for handling mdp input (and tpr I/O).
    IMdpOptionProvider* mdpOptionProvider() override { return &ramdOptions_; }

    //! Returns an interface for handling output files during simulation.
    IMDOutputProvider* outputProvider() override { return this; }

    //! Initializes force providers from this module.
    void initForceProviders(ForceProviders* forceProviders) override
    {
        if (ramdOptions_.active())
        {
            const auto& parameters = ramdOptions_.parameters();
            forceProvider_ = std::make_unique<RAMDForceProvider>(
                parameters,
                ramdSimulationParameters_.periodicBoundaryConditionType(),
                ramdSimulationParameters_.logger()
            );
            forceProviders->addForceProvider(forceProvider_.get(), "RAMD");
        }
    }

    // From IMDOutputProvider
    void initOutput(FILE* fplog, int nfile, const t_filenm fnm[], bool bAppendFiles, const gmx_output_env_t* oenv) override;
    void finishOutput() override;

    //! Subscribe to pre processing notifications
    void subscribeToPreProcessingNotifications(MDModulesNotifiers* notifiers) override
    {
        if (!ramdOptions_.active())
        {
            return;
        }

        // Set Logger during pre-processing
        const auto setLoggerFunction = [this](const MDLogger& logger)
        { ramdOptions_.setLogger(logger); };
        notifiers->preProcessingNotifier_.subscribe(setLoggerFunction);
    }

    //! Subscribe to simulation setup notifications
    void subscribeToSimulationSetupNotifications(MDModulesNotifiers* notifiers) override
    {
        if (!ramdOptions_.active())
        {
            return;
        }
        
        // Reading PBC parameters during simulation setup
        const auto setPeriodicBoundaryContionsFunction = [this](const PbcType& pbc)
        { this->ramdSimulationParameters_.setPeriodicBoundaryConditionType(pbc); };
        notifiers->simulationSetupNotifier_.subscribe(setPeriodicBoundaryContionsFunction);

        // Saving MDLogger during simulation setup
        const auto setLoggerFunction = [this](const MDLogger& logger)
        { this->ramdSimulationParameters_.setLogger(logger); };
        notifiers->simulationSetupNotifier_.subscribe(setLoggerFunction);
    }

    //! Subscribe to simulation run notifications
    void subscribeToSimulationRunNotifications(MDModulesNotifiers* /* notifiers */) override
    {}

private:
    //! The output provider
    // RAMDOutputProvider ramdOutputProvider_;

    //! The options provided for RAMD
    RAMDOptions ramdOptions_;

    //! Object that evaluates the forces
    std::unique_ptr<RAMDForceProvider> forceProvider_;

    /*! \brief Parameters for RAMD that become available at
     * simulation setup time.
     */
    RAMDSimulationParameterSetup ramdSimulationParameters_;

    GMX_DISALLOW_COPY_AND_ASSIGN(RAMD);
};

void RAMD::initOutput(FILE* /*fplog*/,
                      int nfile,
                      const t_filenm fnm[],
                      bool bAppendFiles,
                      const gmx_output_env_t* oenv)
{
    if (opt2bSet("-ramd", nfile, fnm))
    {
        std::string filename = opt2fn("-ramd", nfile, fnm);
        if (bAppendFiles)
        {
            fpRAMD_ = gmx_fio_fopen(filename.c_str(), "a+");
        }
        else
        {
            fpRAMD_ = xvgropen(filename.c_str(),
                               "Test RAMD Ligand-receptor COM distance",
                               "Time (ps)",
                               "Distance (nm)",
                               oenv);
            std::vector<std::string> setnames;
            for (int g = 1; g <= this->ramdOptions_.parameters().ngroups_; ++g)
            {
                setnames.push_back(std::to_string(g));
            }
            xvgrLegend(fpRAMD_, setnames, oenv);
        }
    }
}

void RAMD::finishOutput() {}

} // namespace 

std::unique_ptr<IMDModule> RAMDModuleInfo::create()
{
    return std::make_unique<RAMD>();
}

} // namespace gmx
