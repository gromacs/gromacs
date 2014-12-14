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
#include "energyhandler.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/options.h"
#include "gromacs/commandline/cmdlineparser.h"

namespace gmx
{

EnergyHandler::EnergyHandler(EnergyAnalysisModulePointer module)
{
    // Options for input files
    t0_     = 0;
    t1_     = 0;
    tDelta_ = 0;
    oenv_   = NULL;
    module_ = module;
}

void EnergyHandler::initOptions(Options *options)
{
    options->addOption(FileNameOption("f")
                           .filetype(eftEnergy)
                           .inputFile()
                           .storeVector(&fnEnergy_)
                           .defaultBasename("ener")
                           .multiValue(true)
                           .description("Energy file(s)"));
    // Add options for energy file time control.
    options->addOption(DoubleOption("b").store(&t0_).timeValue()
                           .description("First frame (%t) to read from energy file"));
    options->addOption(DoubleOption("e").store(&t1_).timeValue()
                           .description("Last frame (%t) to read from energy file"));
    options->addOption(DoubleOption("dt").store(&tDelta_).timeValue()
                           .description("Only use frame if t MOD dt == first time (%t)"));

    // Call the module to do it's bit
    module_->initOptions(options);

    // TODO: use the full blown initialization call for this structure
    output_env_init_default(&oenv_);

    module_->setOutputEnvironment(oenv_);
}

int EnergyHandler::readFiles()
{
    ener_file_t fp;
    t_enxframe  frame;
    int         timecheck = 0;
    bool        bCont, bOK = true;

    printf("There are %d energy files registered in the energy handler.\n",
           (int)fnEnergy_.size());
    if ((fnEnergy_.size() == 0) || (NULL == module_))
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

        // Notify tools that we're reading a new data set
        bOK = module_->addDataSet(*fn);

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
            bOK = bOK && module_->initAnalysis(eName, eUnit);

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
                bCont = bCont && module_->addAnalysisFrame(&frame);
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
    bOK = bOK && module_->finalizeAnalysis();

    return 0;
}

}
