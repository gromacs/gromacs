/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <stdlib.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/oenv.h"

#include "poldata.h"
#include "poldata_xml.h"

int alex_poldata_test(int argc, char*argv[])
{
    static const char               *desc[] = {
        "poldata_test reads a poldata (force field) file and writes a new one.",
    };
    gmx_output_env_t                *oenv;
    t_filenm                         fnm[] = {
        { efDAT, "-f", "pdin", ffREAD },
        { efDAT, "-o", "pdout", ffWRITE }
    };
#define NFILE sizeof(fnm)/sizeof(fnm[0])

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, 0, NULL,
                           1, desc, 0, NULL, &oenv))
    {
        return 0;
    }

    gmx_atomprop_t aps  = gmx_atomprop_init();

    alexandria::Poldata pd;
    try 
    {
        alexandria::readPoldata(opt2fn("-f", NFILE, fnm), pd, aps);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    alexandria::writePoldata(opt2fn("-o", NFILE, fnm), pd, 0);

    return 0;
}
