/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2019
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

#include "poldata_tables.h"

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

    lt.setColumns("lccc");

    snprintf(longbuf, STRLEN, "The optimized parameters for the Alexandria charge model. The exponent of the $1s$-Gaussian density function is represented by $\\beta$ in nm$^{-1}$. The atomic electronegativity and absolute hardness are represented by $\\chi$ and $\\eta$, respectively, in eV.");
    lt.setCaption(longbuf);
    lt.setLabel("eemprop");
    snprintf(longbuf, STRLEN, "Alexandria Type & $\\chi$($\\sigma$) & $\\eta$($\\sigma$) & $\\beta$($\\sigma$)");
    lt.addHeadLine(longbuf);
    lt.printHeader();
}

void alexandria_poldata_eemprops_table(FILE                   *fp,
                                       const Poldata          &pd,
                                       ChargeDistributionModel qdist)
{
    char       longbuf[STRLEN];
    LongTable  lt(fp, false, nullptr);

    eemprops_zeta_header(lt);
    auto ztypes = pd.ztype_names();
    for (auto ztp = ztypes.begin(); ztp < ztypes.end(); ztp++)
    {
        auto qDist  = pd.ztype2Eem(qdist,  ztp->c_str());
        auto nzeta  = qDist->getNzeta();
        if (qDist != pd.EndEemprops())
        {
            size_t      pos   = ztp->find("z_");
            std::string ztype = ztp->c_str();
            if (pos != std::string::npos)
            {
                ztype = ztp->substr(pos+2);
            }
            snprintf(longbuf, STRLEN, "%s & %0.2f (%0.2f) & %0.2f (%0.2f) & %0.2f (%0.2f)",
                     ztype.c_str(),
                     qDist->getChi0(),
                     qDist->getChi0_sigma() + 0.005,
                     qDist->getJ0(),
                     qDist->getJ0_sigma() + 0.005,
                     qDist->getZeta(nzeta-1),
                     atof(gmx::splitString(qDist->getZeta_sigma()).back().c_str()) + 0.005);
            lt.printLine(longbuf);
        }
    }
    lt.printFooter();
    fflush(fp);
}

} //namespace
