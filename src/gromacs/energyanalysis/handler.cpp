/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/trxio.h"
#include "handler.h"

namespace gmx
{

bool EnergyHandler::readFiles()
{
    gmx_bool     bCont;
    ener_file_t  fp = open_enx(fileName_.begin()->c_str(), "r");
    t_enxframe   frame;
    int          timecheck = 0;
    bool         bOK;

    printf("There are %d analysis tools registered in the energy handler.\n",
           (int) eap_.size());
    {
        int          nre;
        gmx_enxnm_t *enm = NULL;
        do_enxnms(fp, &nre, &enm);
        bOK = true;
        for (EnergyAnalysisPtrIterator eapi = eap_.begin();
             bOK && (eapi < eap_.end()); ++eapi)
        {
            EnergyAnalysisPtr eap = *eapi;
            bOK = bOK && eap->initAnalysis(nre, enm);
        }
        free_enxnms(nre, enm);
    }
    if (!bOK)
    {
        return false;
    }

    init_enxframe(&frame);
    do
    {
        /* This loop searches for the first frame (when -b option is given),
         * or when this has been found it reads just one energy frame
         */
        do
        {
            bCont = do_enx(fp, &frame);
            if (bCont)
            {
                timecheck = check_times(frame.t);
            }
        }
        while (bCont && (timecheck < 0));

        if ((timecheck == 0) && bCont)
        {
            for (EnergyAnalysisPtrIterator eapi = eap_.begin();
                 bCont && (eapi < eap_.end()); ++eapi)
            {
                EnergyAnalysisPtr eap = *eapi;
                bCont = bCont && eap->addAnalysisFrame(&frame);
            }
        }
    }
    while (bCont);
    close_enx(fp);
    /* Printing a new line, just because the gromacs library prints step info
       while reading. */
    fprintf(stderr, "\n");

    for (EnergyAnalysisPtrIterator eapi = eap_.begin();
         bOK && (eapi < eap_.end()); ++eapi)
    {
        EnergyAnalysisPtr eap = *eapi;
        bOK = bOK && eap->finalizeAnalysis();
    }
    return true;
}

}
