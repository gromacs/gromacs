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
/*! \internal \file
 * \brief
 * Implements classes in simple.h.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_energyanalysis
 */
#include "gmxpre.h"

#include "simple.h"

#include <cstdio>
#include <cstring>

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/viewit.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"

#include "energyhandler.h"
#include "energyhelper.h"
#include "energyinfo.h"
#include "select_energy.h"

namespace gmx
{

namespace energyanalysis
{

/*! \libinternal
 * \brief
 * Extract energy terms from file and print them.
 */
class SimpleEnergyModule : public EnergyAnalysisModule
{
    public:
        //! Constructor
        SimpleEnergyModule();

        //! Get the output file name
        std::string outputFile() { return fnEnergy_; }

        //! Set the output file name
        void setOutputFile(std::string fnEnergy) { fnEnergy_ = fnEnergy; }

        //! Get the fluctuation convergence file name
        std::string fluctConvFile() { return fnFluctConv_; }

        //! Set the summing
        void setSumming(bool bSum) { bSum_ = bSum; }

        //! Get the summing status
        bool sumStatus() { return bSum_; }

        //! Pass the output environment on for printing etc.
        virtual void setOutputEnvironment(output_env_t oenv)
        {
            ehelper_->setOutputEnvironment(oenv);
        }

        //! Initiate the class based on the settings.
        virtual void init(CommandLineModuleSettings * /*settings*/)
        {
        }

        //! Initiate the command line options
        virtual void initOptions(Options                           *options,
                                 ICommandLineOptionsModuleSettings *settings);

        //! Initiate the command line options
        virtual void optionsFinished(Options * /*options*/)
        {
        }

        /*! \brief
         * Does the initiation of the analysis of the file
         * \param[in] eName Names of the energy terms etc.
         * \param[in] eUnit Units of the energy terms etc.
         */
        virtual void initAnalysis(const std::vector<std::string> &eName,
                                  const std::vector<std::string> &eUnit);

        /*! \brief
         * Analyse one frame and stores the results in memory
         * \param[in] fr The energy data frame
         */
        virtual void addAnalysisFrame(t_enxframe *fr);

        //! Finalize reading
        virtual void finalizeAnalysis();

        //! Sum all the energy terms and delete the original data sets
        void sumEnergies();

        //! View the output file(s)
        void viewOutput();

    private:
        //! Output file
        std::string fnEnergy_;

        //! File pointer for storing energies
        FILE        *fp_;

        //! Fluctuation convergence output (typically a xvg file)
        std::string  fnFluctConv_;

        //! File pointer for storing flucutuations
        FILE        *fc_;

        //! Boolean instructing us whether to sum the energy terms
        bool         bSum_;

        //! String to store terms to select
        std::string    term_;
        //! Energy helper class for low level stuff
        EnergyHelper  *ehelper_;

};

SimpleEnergyModule::SimpleEnergyModule()
{
    fc_        = NULL;
    fp_        = NULL;
    bSum_      = false;
    ehelper_   = new(EnergyHelper);
}

void SimpleEnergyModule::initOptions(Options                           *options,
                                     ICommandLineOptionsModuleSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] extracts energy components",
        "data from an energy file. The user is prompted to interactively",
        "select the desired energy terms.[PAR]",

        "Average, RMSD, and drift are calculated with full precision from the",
        "simulation (see printed manual). Drift is calculated by performing",
        "a least-squares fit of the data to a straight line. The reported total drift",
        "is the difference of the fit at the first and last point.",
        "An error estimate of the average is given based on a block averages",
        "over 5 blocks using the full-precision averages. The error estimate",
        "can be performed over multiple block lengths with the options",
        "[TT]-nbmin[tt] and [TT]-nbmax[tt].",
        "[BB]Note[bb] that in most cases the energy files contains averages over all",
        "MD steps, or over many more points than the number of frames in",
        "energy file. This makes the [THISMODULE] statistics output more accurate",
        "than the [TT].xvg[tt] output. When exact averages are not present in the energy",
        "file, the statistics mentioned above are simply over the single, per-frame",
        "energy values.[PAR]"
    };
    settings->setHelpText(desc);

    // Add option for output files
    setOutputFile("energy.xvg");
    options->addOption(FileNameOption("o").filetype(eftPlot).outputFile()
                           .store(&fnEnergy_).defaultBasename("energy")
                           .description("Energy terms as a function of time"));
    options->addOption(FileNameOption("convergence")
                           .filetype(eftPlot)
                           .outputFile()
                           .store(&fnFluctConv_)
                           .defaultBasename("fluct_conv")
                           .description("Convergence of fluctuation properties"));
    options->addOption(BooleanOption("sum")
                           .store(&bSum_)
                           .description("Sum the energy terms selected rather than display them all"));
    options->addOption(StringOption("term")
                           .store(&term_)
                           .description("List the energy terms within quotes on the command line rather than interactively, e.g. '3 4' or 'Coulomb Angle'"));
    ehelper_->initOptions(options);
}

void SimpleEnergyModule::initAnalysis(const std::vector<std::string> &eName,
                                      const std::vector<std::string> &eUnit)
{
    std::vector<int> set;

    // Select which energies to use
    bool bVerbose = (output_env_get_verbosity(ehelper_->outputEnvironment()) > 0);
    set = select_energies(eName, bVerbose, term_);

    for (EnergyTermIterator eti = ehelper_->etBegin(); (eti < ehelper_->etEnd()); ++eti)
    {
        eti->reset();
    }
    for (unsigned int i = 0; (i < set.size()); i++)
    {
        EnergyTerm et(set[i], ehelper_->storeData(),
                      eName[set[i]], eUnit[set[i]]);
        ehelper_->addEnergyTerm(et);
    }

    if (fluctConvFile().size() > 0)
    {
        fc_ = xvgropen(fluctConvFile().c_str(),
                       "Standard Deviation", "Time (ps)",
                       "", ehelper_->outputEnvironment());
    }
    if (outputFile().size() > 0)
    {
        fp_ = xvgropen(outputFile().c_str(),
                       "Energy", "Time (ps)",
                       "", ehelper_->outputEnvironment());
        ehelper_->printXvgLegend(fp_);
    }
}

void SimpleEnergyModule::addAnalysisFrame(t_enxframe *fr)
{
    /* We read a valid frame, so we can use it */
    if (fr->nre > 0)
    {
        if (NULL != fc_)
        {
            fprintf(fc_, "%10g", fr->t);
        }
        if (NULL != fp_)
        {
            fprintf(fp_, "%10g", fr->t);
        }
        for (EnergyTermIterator eti = ehelper_->etBegin(); (eti < ehelper_->etEnd()); ++eti)
        {
            unsigned int findex = eti->fileIndex();
            if (findex < (unsigned int)fr->nre)
            {
                eti->addData(fr->t,
                             fr->step,
                             fr->nsum,
                             fr->ener[findex].esum,
                             fr->ener[findex].eav,
                             fr->ener[findex].e);
                if (NULL != fc_)
                {
                    fprintf(fc_, "  %10g", eti->standardDeviation());
                }
                if (NULL != fp_)
                {
                    double ee = fr->ener[findex].e;
                    if (eti->isEner())
                    {
                        ee /= ehelper_->nMol();
                    }
                    fprintf(fp_, "  %10g", ee);
                }
            }
        }
        if (NULL != fc_)
        {
            fprintf(fc_, "\n");
        }
        if (NULL != fp_)
        {
            fprintf(fp_, "\n");
        }
    }
}

void SimpleEnergyModule::sumEnergies()
{
    if (ehelper_->storeData())
    {
        EnergyTerm etsum(0, true, "Total", "kJ/mol");
        etsum.setStoreData(ehelper_->storeData());

        gmx_int64_t nEner = ehelper_->etBegin()->nEnergy();

        for (gmx_int64_t i = 0; (i < nEner); i++)
        {
            double sum = 0;
            for (EnergyTermIterator eti = ehelper_->etBegin(); (eti < ehelper_->etEnd()); ++eti)
            {
                sum += eti->searchEF(i)->energy();
            }
            // It may be possible to sum the variance and sum of energy terms
            // Please check! Does not look like it, but leaving the comment in
            // anyway. DvdS 2014-05-08.
            etsum.addData(ehelper_->etBegin()->searchEF(i)->time(),
                          i, 1, 0.0, 0.0, sum);
        }
        ehelper_->addEnergyTerm(etsum);
    }
    else
    {
        fprintf(stderr, "WARNING: Energies not stored, so I can not sum them\n");
    }
}

void SimpleEnergyModule::finalizeAnalysis()
{
    if (NULL != fc_)
    {
        xvgrclose(fc_);
    }
    if (NULL != fp_)
    {
        xvgrclose(fp_);
    }
    if (sumStatus())
    {
        sumEnergies();
        ehelper_->printEnergies(outputFile());
    }
    ehelper_->printStatistics(stdout);
    viewOutput();
}

void SimpleEnergyModule::viewOutput()
{
    if (outputFile().size() > 0)
    {
        do_view(ehelper_->outputEnvironment(), outputFile().c_str(), "-nxy");
    }
    if (fluctConvFile().size() > 0)
    {
        do_view(ehelper_->outputEnvironment(), fluctConvFile().c_str(), "-nxy");
    }
}

}   // namespace energyanalysis

const char SimpleInfo::name[] = "energy";

const char SimpleInfo::shortDescription[] = "Extract energy from data files and print graphs";

EnergyInfoPointer SimpleInfo::create()
{
    energyanalysis::EnergyAnalysisModulePointer eamp(new energyanalysis::SimpleEnergyModule);
    return new EnergyInfo(eamp);
}

} // namespace gmx
