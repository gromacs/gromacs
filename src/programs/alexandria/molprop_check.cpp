/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2018 
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Paul J. van Maaren, 
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */
 
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/commandline/pargs.h"

#include "molprop.h"
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
                           0, nullptr,
                           sizeof(desc)/sizeof(desc[0]), desc,
                           0, nullptr, &oenv))
    {
        return 0;
    }
    MolPropRead(opt2fn("-f", NFILE, fnm), mp);

    for (auto &m : mp)
    {
        for (alexandria::ExperimentIterator ci = m.BeginExperiment(); ci < m.EndExperiment(); ++ci)
        {
            int nH = 0, nC = 0;
            for (alexandria::CalcAtomIterator cai = ci->BeginAtom(); cai < ci->EndAtom(); ++cai)
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
