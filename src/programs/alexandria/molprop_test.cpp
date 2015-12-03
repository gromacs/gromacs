/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include <stdio.h>
#include <stdlib.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/topology.h"

#include "molprop.h"
#include "molprop_util.h"
#include "molprop_xml.h"
#include "poldata_xml.h"

int alex_molprop_test(int argc, char*argv[])
{
    static const char               *desc[] = {
        "molprop_test reads a molprop file and writes a new one.",
    };
    gmx_output_env_t                *oenv;
    std::vector<alexandria::MolProp> mpt;
    t_filenm                         fnm[] = {
        { efDAT, "-f", "molin", ffREAD },
        { efDAT, "-o", "molout", ffWRITE }
    };
#define NFILE sizeof(fnm)/sizeof(fnm[0])

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, 0, NULL,
                           1, desc, 0, NULL, &oenv))
    {
        return 0;
    }

    MolPropRead(opt2fn("-f", NFILE, fnm), mpt);
    printf("Read %d molecules from %s\n", (int)mpt.size(), argv[1]);
    MolPropWrite(opt2fn("-o", NFILE, fnm), mpt, 1);

    return 0;
}
