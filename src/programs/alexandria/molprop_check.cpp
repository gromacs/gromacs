/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/commandline/pargs.h"
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

    for (auto &m : mp)
    {
        for (ExperimentIterator ci = m.BeginExperiment(); ci < m.EndExperiment(); ++ci)
        {
            int nH = 0, nC = 0;
            for (CalcAtomIterator cai = ci->BeginAtom(); cai < ci->EndAtom(); ++cai)
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
