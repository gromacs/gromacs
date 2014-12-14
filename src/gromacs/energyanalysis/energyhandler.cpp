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
#include "gmxpre.h"

#include "energyhandler.h"

#include "gromacs/options.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/commandline/cmdlineparser.h"
#include "gromacs/commandline/cmdlineprogramcontext.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/options/timeunitmanager.h"
#include "gromacs/utility/cstringutil.h"

namespace gmx
{

namespace energyanalysis
{

EnergyAnalysisRunner::EnergyAnalysisRunner(EnergyAnalysisModulePointer module) : module_(std::move(module)), timeUnitBehavior_(new TimeUnitBehavior())

{
    // Options for input files
    t0_       = -1;
    t1_       = -1;
    tDelta_   = 0;
    oenv_     = NULL;
    bVerbose_ = true;
    bView_    = false;
}

void EnergyAnalysisRunner::initOptions(IOptionsContainer                 *options,
                                       ICommandLineOptionsModuleSettings *settings)
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

    timeUnitBehavior_->addTimeUnitOption(options, "tu");
    settings->addOptionsBehavior(timeUnitBehavior_);
    plotSettings_.initOptions(options);

    // Call the module to do it's bit
    module_->initOptions(options, settings);
}

int EnergyAnalysisRunner::checkTime(double t)
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

int EnergyAnalysisRunner::run()
{
    ener_file_t fp;
    t_enxframe  frame;

    output_env_init(&oenv_, gmx::getProgramContext(),
                    (time_unit_t)(timeUnitBehavior_->timeUnit() + 1),
                    bView_,
                    static_cast<xvg_format_t>(plotSettings_.plotFormat()),
                    bVerbose_ ? 1 : 0);

    module_->setOutputEnvironment(oenv_);

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

        // Set the energy terms
        fp = open_enx(fn->c_str(), "r");
        do_enxnms(fp, &nre, &enm);
        if (bFirstFile)
        {
            for (int i = 0; (i < nre); i++)
            {
                energyNameUnit enu;
                enu.energyName = enm[i].name;
                enu.energyUnit = enm[i].unit;
                eNU_.push_back(enu);
            }
            module_->initAnalysis(eNU_);

            bFirstFile = false;
        }
        else
        {
            GMX_RELEASE_ASSERT((static_cast<unsigned int>(nre) == eNU_.size()),
                               "Energy files should have the same number of terms");
            for (int i = 0; (i < nre); i++)
            {
                GMX_RELEASE_ASSERT((0 == eNU_[i].energyName.compare(enm[i].name)),
                                   "Energy names mismatch between files");
                GMX_RELEASE_ASSERT((0 == eNU_[i].energyUnit.compare(enm[i].unit)),
                                   "Energy units mismatch between files");
            }
        }
        free_enxnms(nre, enm);
        init_enxframe(&frame);

        gmx_int64_t nframes = 0;
        while (do_enx(fp, &frame))
        {
            /* This loop searches for the first frame (when -b option is given),
             * or when this has been found it reads just one energy frame
             */
            if (0 == checkTime(frame.t))
            {
                module_->analyzeFrame(&frame);
                nframes++;
            }
        }

        close_enx(fp);
        /* Printing a new line, just because the gromacs library prints step info
         * while reading.
         */
        char buf[256];
        fprintf(stderr, "\nRead %s frames from %s\n",
                gmx_step_str(nframes, buf), fn->c_str() );
    }

    // Finally finish the analysis!
    module_->finalizeAnalysis();

    return 0;
}

} // namespace energyanalysis

} // namespace gmx
