/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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
#include "gromacs/options.h"
#include "gromacs/commandline/cmdlineparser.h"
#include "gromacs/commandline/cmdlineprogramcontext.h"
#include "gromacs/options/timeunitmanager.h"

#include "energyhandler.h"

namespace gmx
{

namespace energyanalysis
{

EnergyHandler::EnergyHandler(EnergyAnalysisModulePointer module)
{
    // Options for input files
    t0_       = -1;
    t1_       = -1;
    tDelta_   = 0;
    oenv_     = NULL;
    bVerbose_ = true;
    bView_    = false;
    module_   = module;
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
    options->addOption(BooleanOption("w").store(&bView_)
                           .description("View output [TT].xvg[tt], [TT].xpm[tt], "
                                        "[TT].eps[tt] and [TT].pdb[tt] files"));
    options->addOption(BooleanOption("v").store(&bVerbose_)
                           .description("Verbose output"));
    timeUnitManager_.addTimeUnitOption(options, "tu");

    // Call the module to do it's bit
    module_->initOptions(options);
}

void EnergyHandler::optionsFinished(gmx_unused Options *options)
{
    output_env_init(&oenv_, gmx::getProgramContext(),
                    (time_unit_t)(timeUnitManager_.timeUnit() + 1),
                    bView_, exvgXMGRACE, bVerbose_ ? 1 : 0);

    module_->setOutputEnvironment(oenv_);
}

int EnergyHandler::checkTime(double t)
{
    if ((t0_ >= 0) && (t < t0_))
    {
        return -1;
    }
    else if ((t1_ >= 0) && (t > t1_))
    {
        return 1;
    }
    return 0;
}

int EnergyHandler::run()
{
    ener_file_t fp;
    t_enxframe  frame;
    int         timecheck = 0;
    bool        bCont, bOK = true;

    printf("There are %d energy files registered in the energy handler.\n",
           (int)fnEnergy_.size());
    if (t0_ >= 0)
    {
        printf("Will start reading at %g ps\n", t0_);
    }
    if (t1_ >= 0)
    {
        printf("Will end reading at %g ps\n", t1_);
    }
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
        gmx_int64_t nframes = 0;
        do
        {
            /* This loop searches for the first frame (when -b option is given),
             * or when this has been found it reads just one energy frame
             */
            bCont = do_enx(fp, &frame);
            if (bCont)
            {
                timecheck = checkTime(frame.t);
            }

            if ((timecheck == 0) && bCont)
            {
                bCont = bCont && module_->addAnalysisFrame(&frame);
                nframes++;
            }
        }
        while (bCont && (timecheck <= 0));
        close_enx(fp);
        /* Printing a new line, just because the gromacs library prints step info
         * while reading.
         */
        char buf[256];
        fprintf(stderr, "\nRead %s frames from %s\n",
                gmx_step_str(nframes, buf), fn->c_str() );
    }

    // Finally finish the analysis!
    bOK = bOK && module_->finalizeAnalysis();

    return bOK;
}

}

}
