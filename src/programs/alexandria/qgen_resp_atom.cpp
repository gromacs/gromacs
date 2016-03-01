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
#include "qgen_resp_atom.h"

#include "gromacs/topology/atoms.h"

#include "poldata.h"
#include "qgen_resp.h"
#include "stringutil.h"

namespace alexandria
{

RespAtomType::RespAtomType(int                             atype,
                           int                             particleType,
                           bool                            hasShell,
                           const char                     *atomtype,
                           const Poldata                  &pd,
                           ChargeDistributionModel         iDistributionModel,
                           const std::vector<std::string> &dzatoms)
{
    bool bRestr = false;
    if (!dzatoms.empty())
    {
        int k = 0;
        while (dzatoms[k].size() > 0 && !bRestr)
        {
            bRestr = (strcasecmp(atomtype, dzatoms[k].c_str()) == 0);
            k++;
        }
    }
    atype_       = atype;
    atomtype_    = atomtype;
    bRestrained_ = bRestr;
    bHasShell_   = hasShell;
    int nZeta    = std::max(1, pd.getNzeta(iDistributionModel, atomtype_));
    if (particleType == eptShell)
    {
        int         shell        = nZeta-1;
        size_t      shell_name   = atomtype_.find("_s");
        std::string atomtype_new = atomtype_;
        if (shell_name != std::string::npos)
        {
            shell        = 1;
            atomtype_new = atomtype_.substr(0, shell_name);
        }
        rz_.push_back(RowZetaQ(pd.getRow(iDistributionModel, atomtype_new, shell),
                               pd.getZeta(iDistributionModel, atomtype_new, shell),
                               pd.getQ(iDistributionModel, atomtype_new, shell)));
    }
    else
    {
        if (hasShell)
        {
            nZeta--;
        }
        for (int i = 0; i < nZeta; i++)
        {
            rz_.push_back(RowZetaQ(pd.getRow(iDistributionModel, atomtype, i),
                                   pd.getZeta(iDistributionModel, atomtype, i),
                                   pd.getQ(iDistributionModel, atomtype, i)));
        }
    }
}

} // namespace
