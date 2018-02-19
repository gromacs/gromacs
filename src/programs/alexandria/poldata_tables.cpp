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
 * \author  Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/fileio/xvgr.h"
#include "gromacs/linearalgebra/matrix.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/utility/cstringutil.h"

#include "categories.h"
#include "composition.h"
#include "latex_util.h"
#include "poldata_low.h" 
 
namespace alexandria
{

static void eemprops_zeta_header(LongTable &lt)
{
    char             longbuf[STRLEN];
    CompositionSpecs cs;

    lt.setColumns("ccc");

    snprintf(longbuf, STRLEN, "The optimized orbital exponent for the polarizable Gaussian and the Slater charges represented by $\\beta$ and $\\zeta$ in nm$^{-1}$, rescpectively. The atom types are according to the General Amber Force Field~\\cite{Wang2004a}.");
    lt.setCaption(longbuf);
    lt.setLabel("screeningfactor");
    //snprintf(longbuf, STRLEN, "Atom type  & \\multicolumn{2}{c}{Non-polarizable} & \\multicolumn{2}{c}{Polarizable}");
    //lt.addHeadLine(longbuf);
    snprintf(longbuf, STRLEN, " Atom Type & $\\beta$($\\sigma$) & $\\zeta$($\\sigma$)");
    lt.addHeadLine(longbuf);
    lt.printHeader();
}


void alexandria_poldata_eemprops_zeta_table(FILE           *fp,
                                            const Poldata  &pd)
{
    char       longbuf[STRLEN];
    LongTable  lt(fp, false, nullptr);
    
    eemprops_zeta_header(lt); 
    auto ztypes = pd.ztype_names();      
    for (auto ztp = ztypes.begin(); ztp < ztypes.end(); ztp++)
    {
        auto AXpg = pd.ztype2Eem(eqdAXpg, ztp->c_str());
        auto AXps = pd.ztype2Eem(eqdAXps, ztp->c_str());
        
        if (AXpg != pd.EndEemprops() && AXps != pd.EndEemprops())
        {
            size_t      pos   = ztp->find("z_");
            std::string ztype = ztp->c_str();
            if (pos != std::string::npos)
            {
                ztype = ztp->substr(pos+2);
            }
            
            snprintf(longbuf, STRLEN, "%s & %0.3f (%0.3f) & %0.3f (%0.3f)",
                     ztype.c_str(),
                     AXpg->getZeta(0),
                     atof(gmx::splitString(AXpg->getZeta_sigma()).begin()->c_str()),
                     AXps->getZeta(0),
                     atof(gmx::splitString(AXps->getZeta_sigma()).begin()->c_str()));
            lt.printLine(longbuf);
        }
    }
    lt.printFooter();
    fflush(fp);
}

void alexandria_poldata_eemprops_table(FILE               *fp, 
                                       bool                bzeta,
                                       bool                bchiJ00,
                                       const Poldata       &pd)
{
    if (bzeta)
    {
        alexandria_poldata_eemprops_zeta_table(fp, pd);
    }
    else if(bchiJ00)
    {
        //alexandria_poldata_eemprops_chiJ00_table(fp, pd);
    }
}
                                       
} //namespace
