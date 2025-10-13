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
#include <string>

#include "gromacs/mdtypes/imdoutputprovider.h"

struct gmx_output_env_t;
struct t_filenm;

namespace gmx
{

/*! \internal
 * \brief Handle file output for density guided simulations.
 */
class RAMDOutputProvider final : public IMDOutputProvider
{
public:
    //! Initialize output
    void initOutput(FILE* fplog,
                    int nfile,
                    const t_filenm fnm[],
                    bool bAppendFiles,
                    const gmx_output_env_t* oenv) override;

    //! Add a string to the RAMD output file
    void addLine(const std::string& str);

    //! Finalizes output from a simulation run.
    void finishOutput() override;

private:

    FILE* fpRAMD_;
    
};

} // namespace gmx

#endif
