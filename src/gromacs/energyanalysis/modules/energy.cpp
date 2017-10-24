/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements classes in energy.h.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_energyanalysis
 */
#include "gmxpre.h"

#include "energy.h"

#include <cstdio>
#include <cstring>

#include <memory>

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/energyanalysis/energytermcontainer.h"
#include "gromacs/energyanalysis/select_energy.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/unique_cptr.h"

namespace gmx
{

namespace energyanalysis
{

/*! \libinternal
 * \brief
 * Extract energy terms from file and print them.
 */
class EnergyEnergyModule : public EnergyAnalysisModule
{
    public:
        //! Constructor
        EnergyEnergyModule();

        /*! \brief
         * Initialize the command line settings. Does nothing.
         */
        virtual void init(CommandLineModuleSettings * /*settings*/) {}

        virtual void initOptions(IOptionsContainer                 *options,
                                 ICommandLineOptionsModuleSettings *settings);

        virtual void initAnalysis(ArrayRef<const EnergyNameUnit>  eNU,
                                  const gmx_output_env_t         *oenv);

        virtual void analyzeFrame(t_enxframe *fr, const gmx_output_env_t *oenv);

        virtual void finalizeAnalysis(const gmx_output_env_t *oenv);

        virtual void viewOutput(const gmx_output_env_t *oenv);

    private:
        //! Output file
        std::string               fnEnergy_;

        //! File pointer for storing energies
        gmx::unique_cptr<FILE, xvgrclose> fp_;

        //! Fluctuation convergence output (typically a xvg file)
        std::string               fnFluctConv_;

        //! File pointer for storing fluctuations
        gmx::unique_cptr<FILE, xvgrclose> fc_;

        //! Boolean instructing us whether to sum the energy terms
        bool                      bSum_;

        //! Boolean telling whether or not to print high precision
        bool                      bDouble_;

        //! The list of energy terms to extract
        std::vector<std::string>  termList_;

        //! Energy helper class for low level stuff
        EnergyTermContainer       ehelper_;
};

EnergyEnergyModule::EnergyEnergyModule() : fp_(nullptr), fc_(nullptr), bSum_(false), bDouble_(false)
{
    ehelper_.setStoreData(true);
}

void EnergyEnergyModule::initOptions(IOptionsContainer                 *options,
                                     ICommandLineOptionsModuleSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] extracts energy components",
        "data from an energy file. The user is prompted to interactively",
        "select the desired energy terms.",
        "",
        "Average, RMSD, and drift are calculated with full precision from the",
        "simulation (see printed manual). Drift is calculated by performing",
        "a least-squares fit of the data to a straight line. The reported total drift",
        "is the difference of the fit at the first and last point.",
        "An error estimate of the average is given based on a block averages",
        "over 5 blocks using the full-precision averages. The error estimate",
        "can be performed over multiple block lengths with the options",
        "[TT]-nbmin[tt] and [TT]-nbmax[tt].",
        "",
        "[BB]Note[bb] that in most cases the energy files contains averages over all",
        "MD steps, or over many more points than the number of frames in",
        "energy file. This makes the [THISMODULE] statistics output more accurate",
        "than the [TT].xvg[tt] output. When exact averages are not present in the energy",
        "file, the statistics mentioned above are simply over the single, per-frame",
        "energy values."
    };
    settings->setHelpText(desc);

    // Add option for output files
    fnEnergy_.assign("energy.xvg");
    options->addOption(FileNameOption("o").filetype(eftPlot).outputFile()
                           .store(&fnEnergy_).defaultBasename("energy")
                           .description("Energy terms as a function of time"));
    options->addOption(FileNameOption("convergence")
                           .filetype(eftPlot)
                           .outputFile()
                           .store(&fnFluctConv_)
                           .defaultBasename("fluct_conv")
                           .description("Convergence of fluctuation properties shown as standard deviation"));
    options->addOption(BooleanOption("sum")
                           .store(&bSum_)
                           .description("Sum the energy terms selected rather than display them all"));
    options->addOption(BooleanOption("double")
                           .store(&bDouble_)
                           .description("Print results with 15 digits instead of 10"));
    options->addOption(StringOption("term")
                           .storeVector(&termList_).multiValue()
                           .description("Provide list of energy terms on the command line rather than interactively"));
    ehelper_.initOptions(options);
}

void EnergyEnergyModule::initAnalysis(ArrayRef<const EnergyNameUnit>  eNU,
                                      const gmx_output_env_t         *oenv)
{
    std::vector<int> set;

    // Select which energies to use
    bool bVerbose = (output_env_get_verbosity(oenv) > 0);
    if (termList_.size() == 0)
    {
        try
        {
            select_energies(eNU, bVerbose, &StandardInputStream::instance(), set);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }
    else
    {
        try
        {
            StringInputStream input(termList_);
            select_energies(eNU, false, &input, set);
            input.close();
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }

    for (unsigned int i = 0; i < set.size(); i++)
    {
        EnergyTerm et(set[i], ehelper_.storeData(),
                      eNU[set[i]].energyName, eNU[set[i]].energyUnit);
        ehelper_.addEnergyTerm(et);
    }

    if (fnFluctConv_.size() > 0)
    {
        fc_.reset(xvgropen(fnFluctConv_.c_str(),
                           "Standard Deviation", "Time (ps)",
                           "", oenv));
    }
    if (fnEnergy_.size() > 0)
    {
        std::string yaxis;
        yAxis(ehelper_.begin(), ehelper_.end(), &yaxis);
        fp_.reset(xvgropen(fnEnergy_.c_str(),
                           "GROMACS Energies", "Time (ps)",
                           yaxis.c_str(), oenv));
        printXvgLegend(fp_.get(), ehelper_.begin(), ehelper_.end(), oenv);
    }
}

void EnergyEnergyModule::analyzeFrame(t_enxframe *fr, gmx_unused const gmx_output_env_t *oenv)
{
    /* We read a valid frame, so we can use it */
    if (fr->nre > 0)
    {
        if (nullptr != fc_)
        {
            fprintf(fc_.get(), "%10g", fr->t);
        }
        if (nullptr != fp_)
        {
            fprintf(fp_.get(), "%10g", fr->t);
        }
        ehelper_.addFrame(fr);
        for (auto &eti : ehelper_)
        {
            unsigned int findex = eti.fileIndex();
            if (findex < (unsigned int)fr->nre)
            {
                if (nullptr != fc_)
                {
                    fprintf(fc_.get(), "  %10g", eti.standardDeviation());
                }
                if (nullptr != fp_)
                {
                    double ee = fr->ener[findex].e;
                    if (eti.isEner())
                    {
                        ee /= ehelper_.nMol();
                    }
                    if (bDouble_)
                    {
                        fprintf(fp_.get(), "  %15e", ee);
                    }
                    else
                    {
                        fprintf(fp_.get(), "  %10g", ee);
                    }
                }
            }
        }
        if (nullptr != fc_)
        {
            fprintf(fc_.get(), "\n");
        }
        if (nullptr != fp_)
        {
            fprintf(fp_.get(), "\n");
        }
    }
}

static void sumEnergies(EnergyTermContainer *ehelper)
{
    if (ehelper->storeData())
    {
        EnergyTerm etsum(0, true, "Total", "kJ/mol");
        etsum.setStoreData(ehelper->storeData());

        gmx_int64_t nEner = ehelper->begin()->numFrames();

        for (gmx_int64_t i = 0; (i < nEner); i++)
        {
            double sum = 0;
            for (auto &eti : *ehelper)
            {
                sum += eti.findFrame(i)->energy();
            }
            // TODO:
            // It may be possible to sum the variance and sum of energy terms
            // Please check! Does not look like it, but leaving the comment in
            // anyway. DvdS 2014-05-08.
            etsum.addData(ehelper->begin()->findFrame(i)->time(),
                          i, 1, 0.0, 0.0, sum);
        }
        ehelper->addEnergyTerm(etsum);
    }
    else
    {
        fprintf(stderr, "WARNING: Energies not stored, so I can not sum them\n");
    }
}

void EnergyEnergyModule::finalizeAnalysis(const gmx_output_env_t *oenv)
{
    if (bSum_)
    {
        sumEnergies(&ehelper_);
        printEnergies(fnEnergy_, ehelper_.begin(), ehelper_.end(), bDouble_, oenv);
    }
    printStatistics(stdout, ehelper_.begin(), ehelper_.end(), ehelper_.nBlocks());
}

void EnergyEnergyModule::viewOutput(const gmx_output_env_t *oenv)
{
    if (fnEnergy_.size() > 0)
    {
        do_view(oenv, fnEnergy_.c_str(), "-nxy");
    }
    if (fnFluctConv_.size() > 0)
    {
        do_view(oenv, fnFluctConv_.c_str(), "-nxy");
    }
}

const char EnergyInfo::name[] = "energy";

const char EnergyInfo::shortDescription[] = "Extract energy from data files and print graphs";

EnergyAnalysisModulePointer EnergyInfo::create()
{
    return EnergyAnalysisModulePointer(new EnergyEnergyModule());
}

} // namespace energyanalysis

} // namespace gmx
