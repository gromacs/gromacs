/*! \internal \file
 * \brief
 * Declares data structure for RAMD output
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 * \ingroup module_applied_forces
 */

#include "gmxpre.h"

#include "ramdoutputprovider.h"
#include "gromacs/commandline/filenm.h"

struct gmx_output_env_t;

namespace gmx
{

void RAMDOutputProvider::initOutput(FILE* /*fplog*/,
                                    int nfile,
                                    const t_filenm fnm[],
                                    bool /*bAppendFiles*/,
                                    const gmx_output_env_t* /*oenv*/)
{
    // if (mpiComm.isMainRank() and opt2bSet("-ramd", nfile, fnm))
    // {
    //     auto filename = std::string(opt2fn("-ramd", nfile, fnm));

    //     if (startingBehavior == gmx::StartingBehavior::RestartWithAppending)
    //     {
    //         out = gmx_fio_fopen(filename.c_str(), "a+");
    //     }
    //     else
    //     {
    //         out = gmx_fio_fopen(filename.c_str(), "w+");
    //         xvgr_header(out, "RAMD Ligand-receptor COM distance", "Time (ps)", "Distance (nm)", exvggtXNY, oenv);

    //         std::vector<std::string> setnames;
    //         for (int g = 1; g <= params.ngroup; ++g)
    //         {
    //             setnames.push_back(std::to_string(g));
    //         }
    //         xvgrLegend(out, setnames, oenv);
    //     }
    //     fflush(out);
    // }
}

void RAMDOutputProvider::finishOutput() {}

} // namespace gmx
