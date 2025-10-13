/*! \internal \file
 * \brief
 * Declares force provider for RAMD
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_RAMDFORCEPROVIDER_H
#define GMX_APPLIED_FORCES_RAMDFORCEPROVIDER_H

#include <string>

#include "gromacs/domdec/localatomset.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/logger.h"

#include "ramdoutputprovider.h"
#include "ramdparameters.h"
#include "randomsphericaldirectiongenerator.h"

namespace gmx
{

/*! \internal \brief
 * Declares IForceProvider for RAMD.
 */
class RAMDForceProvider final : public IForceProvider
{
public:
    RAMDForceProvider(const RAMDParameters& parameters,
                      const std::vector<std::unique_ptr<LocalAtomSet>>& localAtoms,
                      PbcType               pbcType,
                      const MDLogger&       logger,
                      RAMDOutputProvider&   ramdOutputProvider);

    //! Destruct force provider for RAMD
    ~RAMDForceProvider();

    /*!\brief Calculate forces of RAMD.
     * \param[in] fInput input for force provider
     * \param[out] fOutput output for force provider
     */
    void calculateForces(const ForceProviderInput& fInput, ForceProviderOutput* fOutput) override;

private:

    const RAMDParameters& parameters_;

    //! Reference to local atom sets
    const std::vector<std::unique_ptr<LocalAtomSet>>& localAtoms_;

    const PbcType         pbcType_;
    const MDLogger&       logger_;
    RAMDOutputProvider&   ramdOutputProvider_;

    //! Random pull direction
    RandomSphericalDirectionGenerator random_spherical_direction_generator;

    //! Current pull direction
    std::vector<DVec> direction_;

    //! COM of receptor of last RAMD evaluation step
    std::vector<DVec> com_rec_prev_;

    //! COM of ligand of last RAMD evaluation step
    std::vector<DVec> com_lig_prev_;

    /// Has the ligand left his binding site?
    std::vector<int> ligand_exited_;

    /// Control trajectory output
    gmx_bool write_trajectory_;

};

} // namespace gmx

#endif // GMX_APPLIED_FORCES_RAMDFORCEPROVIDER_H
