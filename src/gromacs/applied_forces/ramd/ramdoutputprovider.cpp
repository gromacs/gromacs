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
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"

struct gmx_output_env_t;

namespace gmx
{

void RAMDOutputProvider::initOutput(FILE* fplog,
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
                               "RAMD Ligand-receptor COM distance",
                               "Time (ps)",
                               "Distance (nm)",
                               oenv);
            // std::vector<std::string> setnames;
            // for (int g = 1; g <= params.ngroup; ++g)
            // {
            //     setnames.push_back(std::to_string(g));
            // }
            // xvgrLegend(fpRAMD_, setnames, oenv);
        }
    }
}

void RAMDOutputProvider::addTimePoint(double time)
{
    if (fpRAMD_)
    {
        fprintf(fpRAMD_, "%.4f", time);
    }
}

void RAMDOutputProvider::finishOutput() {}

} // namespace gmx
