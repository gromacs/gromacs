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

    lt.setColumns("lcc");

    snprintf(longbuf, STRLEN, "The optimized exponent for the polarizable Gaussian and Slater $s$-type orbitals represented by $\\beta$ and $\\zeta$ in nm$^{-1}$, rescpectively. Slater 3s orbital has been optimized rather than the valence Slater $s$-orbital for Bromine and Iodine ({\\it See} THEORY).");
    lt.setCaption(longbuf);
    lt.setLabel("orbitalexpoenent");
    snprintf(longbuf, STRLEN, "Polarizability Type & $\\beta$($\\sigma$) & $\\zeta$($\\sigma$)");
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
            
            snprintf(longbuf, STRLEN, "%s & %0.2f (%0.2f) & %0.2f (%0.2f)",
                     ztype.c_str(),
                     AXpg->getZeta(1),
                     atof(gmx::splitString(AXpg->getZeta_sigma()).back().c_str()) + 0.005,
                     AXps->getZeta(1),
                     atof(gmx::splitString(AXps->getZeta_sigma()).back().c_str()) + 0.005);
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
