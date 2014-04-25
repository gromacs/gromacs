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
#include "handler.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/options.h"
#include "gromacs/commandline/cmdlineparser.h"

namespace gmx
{

EnergyHandler::EnergyHandler() : options_("handler", "Energy file hander options")
{
    // Options for input files
    fnEnergy_.push_back("ener.edr");
    options_.addOption(FileNameOption("f")
                           .filetype(eftEnergy)
                           .inputFile()
                           .storeVector(&fnEnergy_)
                           .defaultBasename("ener")
                           .multiValue(true)
                           .description("Energy file(s)"));
    t0_ = 0;
    // Add options for energy file time control.
    options_.addOption(DoubleOption("b").store(&t0_).timeValue()
                           .description("First frame (%t) to read from energy file"));
    t1_ = 0;
    options_.addOption(DoubleOption("e").store(&t1_).timeValue()
                           .description("Last frame (%t) to read from energy file"));
    tDelta_ = 0;
    options_.addOption(DoubleOption("dt").store(&tDelta_).timeValue()
                           .description("Only use frame if t MOD dt == first time (%t)"));
    oenv_ = NULL;
}

void EnergyHandler::prepare(int *argc, char *argv[])
{
    for (EnergyAnalysisModulePointerIterator eapi = eap_.begin();
         (eapi < eap_.end()); ++eapi)
    {
        (*eapi)->initOptions(&options_);
    }

    // Stuff below copied from trajectoryanalysis/runnercommon.cpp
    CommandLineParser  parser(&options_);
    // TODO: Print the help if user provides an invalid option?
    // Or just add a message advising the user to invoke the help?
    parser.parse(argc, argv);
    // TODO: Implement time scaling.
    //common->scaleTimeOptions(&options_);
    options_.finish();

    // TODO: use the full blown initialization call for this structure
    output_env_init_default(&oenv_);
    for (EnergyAnalysisModulePointerIterator eapi = eap_.begin();
         (eapi < eap_.end()); ++eapi)
    {
        (*eapi)->setOutputEnvironment(oenv_);
    }
}

int EnergyHandler::readFiles()
{
    ener_file_t fp;
    t_enxframe  frame;
    int         timecheck = 0;
    bool        bCont, bOK = true;

    printf("There are %d analysis tools and %d energy files registered in the energy handler.\n",
           (int)eap_.size(), (int)fnEnergy_.size());
    if ((fnEnergy_.size() == 0) || (eap_.size() == 0))
    {
        printf("Nothing to do!\n");
        return -1;
    }
    bool bFirstFile = true;
    for (std::vector<std::string>::iterator fn = fnEnergy_.begin();
         (fn < fnEnergy_.end()); ++fn)
    {
        int          nre;
        gmx_enxnm_t *enm = NULL;

        bOK = true;

        // Notify tools that we're reading a new data set
        for (EnergyAnalysisModulePointerIterator eapi = eap_.begin();
             bOK && (eapi < eap_.end()); ++eapi)
        {
            bOK = bOK && (*eapi)->addDataSet(*fn);
        }

        // Set the energy terms
        fp = open_enx(fn->c_str(), "r");
        do_enxnms(fp, &nre, &enm);
        if (bFirstFile)
        {
            std::vector<std::string> eName, eUnit;
            for (int i = 0; (i < nre); i++)
            {
                eName.push_back(enm[i].name);
                eUnit.push_back(enm[i].unit);
            }
            for (EnergyAnalysisModulePointerIterator eapi = eap_.begin();
                 bOK && (eapi < eap_.end()); ++eapi)
            {
                EnergyAnalysisModulePointer eap = *eapi;
                bOK = bOK && eap->initAnalysis(eName, eUnit);
            }
            bFirstFile = false;
        }
        free_enxnms(nre, enm);
        init_enxframe(&frame);

        if (!bOK)
        {
            return -1;
        }

        do
        {
            /* This loop searches for the first frame (when -b option is given),
             * or when this has been found it reads just one energy frame
             */
            if ((bCont = do_enx(fp, &frame)))
            {
                timecheck = check_times(frame.t);
            }
            while (bCont && (timecheck < 0))
            {
                if ((bCont = do_enx(fp, &frame)))
                {
                    timecheck = check_times(frame.t);
                }
            }

            if ((timecheck == 0) && bCont)
            {
                for (EnergyAnalysisModulePointerIterator eapi = eap_.begin();
                     bCont && (eapi < eap_.end()); ++eapi)
                {
                    EnergyAnalysisModulePointer eap = *eapi;
                    bCont = bCont && eap->addAnalysisFrame(&frame);
                }
            }
        }
        while (bCont);
        close_enx(fp);
        /* Printing a new line, just because the gromacs library prints step info
         * while reading.
         */
        fprintf(stderr, "\n");
    }
    // Finally finish the analysis!
    for (EnergyAnalysisModulePointerIterator eapi = eap_.begin();
         bOK && (eapi < eap_.end()); ++eapi)
    {
        EnergyAnalysisModulePointer eap = *eapi;
        bOK = bOK && eap->finalizeAnalysis();
    }
    return 0;
}

}
