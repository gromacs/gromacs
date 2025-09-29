/*! \internal \file
 * \brief
 * Declares output to file for RAMD
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_RAMDOUTPUTPROVIDER_H
#define GMX_APPLIED_FORCES_RAMDOUTPUTPROVIDER_H

#include <cstdio>

#include "gromacs/mdtypes/imdoutputprovider.h"

namespace gmx
{

/*! \internal
 * \brief Handle file output for density guided simulations.
 */
class RAMDOutputProvider : public IMDOutputProvider
{
protected:
    FILE* fpRAMD_;
};

} // namespace gmx

#endif
