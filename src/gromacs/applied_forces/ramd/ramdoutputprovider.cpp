/*! \internal \file
 * \brief
 * Declares data structure for RAMD output
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 * \ingroup module_applied_forces
 */

#include "gmxpre.h"

#include "ramdoutputprovider.h"

struct gmx_output_env_t;

namespace gmx
{

void RAMDOutputProvider::initOutput(FILE* /*fplog*/,
                                    int /*nfile*/,
                                    const t_filenm /*fnm*/[],
                                    bool /*bAppendFiles*/,
                                    const gmx_output_env_t* /*oenv*/)
{}

void RAMDOutputProvider::finishOutput() {}

} // namespace gmx
