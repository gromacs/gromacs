/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/copyrite.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#include "molprop.h"
#include "molprop_util.h"
#include "molprop_xml.h"
#include "poldata_xml.h"

int alex_molprop_check(int argc, char*argv[])
{
    static const char               *desc[] = {
        "molprop_check checks calculations for missing hydrogens" 
    };
    t_filenm                         fnm[] =
    {
        { efDAT, "-f",  "allmols",  ffREAD }
    };
    int                              NFILE   = (sizeof(fnm)/sizeof(fnm[0]));
    std::vector<alexandria::MolProp> mp;
    gmx_output_env_t                *oenv;

    if (!parse_common_args(&argc, argv, PCA_NOEXIT_ON_ARGS, NFILE, fnm,
                           0, NULL,
                           sizeof(desc)/sizeof(desc[0]), desc,
                           0, NULL, &oenv))
    {
        return 0;
    }
    MolPropRead(opt2fn("-f", NFILE, fnm), mp);
    
    for(auto &m : mp) 
    {
        for(CalculationIterator ci = m.BeginCalculation(); ci < m.EndCalculation(); ++ci)
        {
            int nH = 0, nC = 0;
            for(CalcAtomIterator cai = ci->BeginAtom(); cai < ci->EndAtom(); ++cai)
            {
                std::string name = cai->getName();
                if (name.compare("H") == 0)
                {
                    nH++;
                }
                else if (name.compare("C") == 0)
                {
                    nC++;
                }
            }
            if (nC > 0 && nH == 0)
            {
                printf("%s #C %d #H %d\n", 
                       ci->getDatafile().c_str(),
                       nC, nH);
            }
        }
    }
    
    return 0;
}
