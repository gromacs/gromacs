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

#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/logger.h"

namespace gmx
{

/*! \internal \brief
 * Implements IForceProvider for RAMD.
 */
class RAMDForceProvider final : public IForceProvider
{
public:
    RAMDForceProvider(const RAMDParameters& parameters,
                      PbcType               pbcType,
                      const MDLogger&       logger);

    //! Destruct force provider for RAMD
    ~RAMDForceProvider();

    /*!\brief Calculate forces of RAMD.
     * \param[in] fInput input for force provider
     * \param[out] fOutput output for force provider
     */
    void calculateForces(const ForceProviderInput& fInput, ForceProviderOutput* fOutput) override;

};

} // namespace gmx

#endif // GMX_APPLIED_FORCES_RAMDFORCEPROVIDER_H
