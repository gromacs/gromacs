/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
 * Implements classes in viscosity.h.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_energyanalysis
 */
#include "gmxpre.h"

#include "viscosity.h"

#include <cstdio>
#include <cstring>

#include <algorithm>
#include <array>

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/correlationfunctions/autocorr.h"
#include "gromacs/correlationfunctions/expfit.h"
#include "gromacs/energyanalysis/analysismodule.h"
#include "gromacs/energyanalysis/energytermcontainer.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/unique_cptr.h"

namespace gmx
{

namespace energyanalysis
{

//! Index into energy array.
enum vIndex {
    vPresXX      = 0, vPresXY = 1, vPresXZ = 2,
    vPresYX      = 3, vPresYY = 4, vPresYZ = 5,
    vPresZX      = 6, vPresZY = 7, vPresZZ = 8,
    vTemperature = 9, vVolume = 10, vPressure = 11,
    vNR          = 12
};

class ViscosityModule : public EnergyAnalysisModule
{
    public:
        //! Constructor
        ViscosityModule();

        virtual void initOptions(IOptionsContainer                 *options,
                                 ICommandLineOptionsModuleSettings *settings);

        virtual void initAnalysis(ArrayRef<const EnergyNameUnit>  eNU,
                                  const gmx_output_env_t         *oenv);

        virtual void analyzeFrame(t_enxframe *fr, const gmx_output_env_t *oenv);

        virtual void finalizeAnalysis(const gmx_output_env_t *oenv);

        virtual void viewOutput(gmx_unused const gmx_output_env_t *oenv) {}

    private:
        //! Filename for Einstein viscosity
        std::string                        fnEinstein_;
        //! Filename for Einstein viscosity integral
        std::string                        fnEinsteinIntegral_;
        //! Filename for Pressure Autocorrelation functions
        std::string                        fnPressureAcf_;
        //! User given volume of the box
        double                             volume_;
        //! Use differences between P trace components (Pxx-Pyy, Pyy-Pzz);
        bool                               usePressureTraceDiff_;
        //! Names of energy terms corresponding to what is needed in this analysis
        const std::array<std::string, vNR> setnm_;
        /*! \brief
         * Compute viscosity using Einstein relation, see GROMACS manual.
         * The output is written to two files.
         * \param[in] sum The sum data sets
         * \param[in] V The average volume
         * \param[in] T The average temperature
         * \param[in] oenv GROMACS output environment
         */
        void doEinstein(std::vector < std::vector < real> >  sum,
                        real                                 V,
                        real                                 T,
                        const gmx_output_env_t              *oenv);
        //! Energy helper class for low level stuff
        EnergyTermContainer   ehelper_;

};

//! The strings corresponding to the energy terms needed for viscosity
static const std::array<std::string, vNR> setnames = {
    { "Pres-XX",  "Pres-XY",  "Pres-XZ",
      "Pres-YX",  "Pres-YY",  "Pres-YZ",
      "Pres-ZX",  "Pres-ZY",  "Pres-ZZ",
      "Temperature",  "Volume",  "Pressure" }
};

ViscosityModule::ViscosityModule() : volume_(0.0),
                                     usePressureTraceDiff_(true),
                                     setnm_(setnames)
{
}

//! Set the options and setting
void ViscosityModule::initOptions(IOptionsContainer                 *options,
                                  ICommandLineOptionsModuleSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] computes shear- and bulk viscosity using a",
        "number of different algorithms. The best value for the",
        "viscosity is obtained from the integral of the pressure",
        "autocorrelation function. Convergence plots can be obtained",
        "from the options [TT]evisco[tt] and [TT]eviscoi[tt] or",
        "the option [TT]-pcorr[tt]."
    };
    // Add the descriptive text (program help text) to the options
    settings->setHelpText(desc);

    // Add option for optional output files
    options->addOption(FileNameOption("evisco").filetype(eftPlot).outputFile()
                           .store(&fnEinstein_).defaultBasename("evisco")
                           .description("Shear viscosity from Einstein relation"));
    options->addOption(FileNameOption("eviscoi").filetype(eftPlot).outputFile()
                           .store(&fnEinsteinIntegral_).defaultBasename("eviscointegral")
                           .description("Integral of shear viscosity from Einstein relation"));
    options->addOption(BooleanOption("usePressureTraceDiff").store(&usePressureTraceDiff_)
                           .description("Use differences between P trace components (Pxx-Pyy, Pyy-Pzz) to compute shear viscosity"));
    options->addOption(DoubleOption("volume")
                           .store(&volume_)
                           .description("Volume of the computational box (nm^3) in order to be able to use the tool with constant volume simulations."));
    ehelper_.initOptions(options);
}

void ViscosityModule::initAnalysis(ArrayRef<const EnergyNameUnit>    eNU,
                                   gmx_unused const gmx_output_env_t *)
{
    bool bFound[vNR];
    int  nfound = 0;
    for (int j = 0; (j < vNR); j++)
    {
        bFound[j] = false;
    }
    unsigned int i = 0;
    for (auto enu : eNU)
    {
        for (int j = 0; (j < vNR); j++)
        {
            if (enu.energyName == setnm_[j])
            {
                EnergyTerm et(i, ehelper_.storeData(), enu.energyName,
                              enu.energyUnit);
                ehelper_.addEnergyTerm(et);
                bFound[j] = true;
                nfound++;
            }
        }
        i++;
    }
    if (!bFound[vVolume] && (volume_ > 0))
    {
        // User set the volume on the command line.
        bFound[vVolume] = true;
        nfound++;
    }
    for (int j = 0; (j < vNR); j++)
    {
        if (!bFound[j])
        {
            fprintf(stderr, "Energy file does not contain %s. Can not compute viscosity.\n",
                    setnm_[j].c_str());
        }
    }
    // Tell the helper to store all it's data
    ehelper_.setStoreData(true);

    if (nfound != vNR)
    {
        GMX_THROW(InvalidInputError("Not enough energy terms in the input file"));
    }
}

void ViscosityModule::analyzeFrame(t_enxframe *fr, gmx_unused const gmx_output_env_t *oenv)
{
    ehelper_.addFrame(fr);
}

//! Compute the average of two energy components
static double calcEintSum(int64_t i, const EnergyTermIterator a, const EnergyTermIterator b)
{
    return 0.5*(a->findFrame(i)->energySum() + b->findFrame(i)->energySum())/a->findFrame(i)->nSum();
}

//! Compute half the difference of two energy coponents
static double calcEintDiff(int64_t i, const EnergyTermIterator a, const EnergyTermIterator b)
{
    return 0.5*(a->findFrame(i)->energySum() - b->findFrame(i)->energySum())/a->findFrame(i)->nSum();
}

void ViscosityModule::finalizeAnalysis(const gmx_output_env_t *oenv)
{
    double             eAver[vNR], eStddev[vNR];
    EnergyTermIterator etxyz[vNR];
    double             Dt;
    int64_t            nsteps, nframes;

    for (int i = 0; (i < vNR); i++)
    {
        if (!ehelper_.energyTerm(setnm_[i], &eAver[i], &eStddev[i]))
        {
            if (vVolume == i)
            {
                if (volume_ > 0)
                {
                    eAver[vVolume]   = volume_;
                    eStddev[vVolume] = 0;
                }
                else
                {
                    fprintf(stderr, "Could not extract %s from energy file.\n",
                            setnm_[i].c_str());
                }
            }
        }
        else
        {
            etxyz[i] = ehelper_.etSearch(setnm_[i]);
            if (nullptr != debug)
            {
                fprintf(debug, "Found %s for %s\n",
                        etxyz[i]->name().c_str(), setnm_[i].c_str());
            }
        }
    }

    nsteps  = ehelper_.begin()->numSteps();
    nframes = ehelper_.begin()->numFrames();
    if ((nsteps < 2) || (nframes < 2))
    {
        char buf1[64], buf2[64];
        fprintf(stderr, "Not enough data: %s steps and %s frames.\n",
                gmx_step_str(nsteps, buf1),
                gmx_step_str(nframes, buf2));
        if (!ehelper_.storeData())
        {
            GMX_THROW(APIError("You forgot to turn on storing data"));
        }
    }
    Dt      = ehelper_.begin()->timeSpan() / (nframes-1);

    int maxNumberComponents = 5;
    int numberComponents    = 3;
    if  (usePressureTraceDiff_)
    {
        numberComponents = maxNumberComponents;
    }
    // Do viscosity calculation using the Einstein method
    if (!fnEinstein_.empty() || !fnEinsteinIntegral_.empty())
    {
        std::vector < std::vector < real> > eneint;
        std::vector<real>                   zero;
        zero.push_back(0.0);
        for (int i = 0; (i < numberComponents); i++)
        {
            eneint.push_back(zero);
        }
        // Assume all components are stored equally often
        int    nSum = etxyz[vPresXY]->findFrame(0)->nSum();
        GMX_ASSERT (0 < nSum, "Number of data points should be > 0");
        for (int64_t i = 0; (i < nframes); i++)
        {
            eneint[0].push_back(eneint[0].back() +
                                Dt*calcEintSum(i, etxyz[vPresXY], etxyz[vPresYX]));
            eneint[1].push_back(eneint[1].back() +
                                Dt*calcEintSum(i, etxyz[vPresXZ], etxyz[vPresZX]));
            eneint[2].push_back(eneint[2].back() +
                                Dt*calcEintSum(i, etxyz[vPresYZ], etxyz[vPresZY]));
            if (usePressureTraceDiff_)
            {
                eneint[3].push_back(eneint[3].back() +
                                    Dt*calcEintDiff(i, etxyz[vPresXX], etxyz[vPresYY]));
                eneint[4].push_back(eneint[4].back() +
                                    Dt*calcEintDiff(i, etxyz[vPresYY], etxyz[vPresZZ]));
            }
        }
        doEinstein(eneint, eAver[vVolume], eAver[vTemperature], oenv);
    }

#undef maxNumberComponents
}

void ViscosityModule::doEinstein(std::vector < std::vector < real> > sum,
                                 real V, real T,
                                 const gmx_output_env_t *oenv)
{
    std::vector<double> av, avold, gamma;
    double              avoldtot, gammatot;
    int                 i, nf4;

    int64_t             nframes = ehelper_.begin()->numFrames()+1;
    double              t0      = ehelper_.begin()->timeBegin();
    if (nframes < 4)
    {
        fprintf(stderr, "Number of frames %d too few for computing Einstein viscosity.\n", (int)nframes);
        return;
    }

    double              Dt      = ehelper_.begin()->timeSpan() / (nframes-2);
    if (Dt <= 0)
    {
        fprintf(stderr, "Delta T incorrect (%g) should be > 0\n", Dt);
        return;
    }
    nf4 = nframes/4+1;

    const char *legend[6];
    int         numberGraphs = 0;
    legend[numberGraphs++] = "(P\\sxy\\N+P\\syx\\N)/2";
    legend[numberGraphs++] = "(P\\sxz\\N+P\\szx\\N)/2";
    legend[numberGraphs++] = "(P\\syz\\N+P\\szy\\N)/2";
    if (usePressureTraceDiff_)
    {
        legend[numberGraphs++] = "P\\sxx\\N-P\\syy\\N";
        legend[numberGraphs++] = "P\\syy\\N-P\\szz\\N";
    }
    legend[numberGraphs++] = "Average";

    gmx::unique_cptr<FILE, xvgrclose> fpi(nullptr);
    if (!fnEinsteinIntegral_.empty())
    {
        fpi.reset(xvgropen(fnEinsteinIntegral_.c_str(),
                           "Shear viscosity integral",
                           "Time (ps)", "(kg m\\S-1\\N s\\S-1\\N ps)",
                           oenv));
        xvgr_legend(fpi.get(), numberGraphs, legend, oenv);
    }
    gmx::unique_cptr<FILE, xvgrclose> fp(nullptr);
    if (!fnEinstein_.empty())
    {
        fp.reset(xvgropen(fnEinstein_.c_str(),
                          "Shear viscosity using Einstein relation",
                          "Time (ps)", "(kg m\\S-1\\N s\\S-1\\N)",
                          oenv));
        xvgr_legend(fp.get(), numberGraphs, legend, oenv);
    }
    for (unsigned int m = 0; m < sum.size(); m++)
    {
        av.push_back(0);
        avold.push_back(0);
        gamma.push_back(0);
    }
    avoldtot = gammatot = 0;
    for (i = 0; i < nf4; i++)
    {
        double avtot = 0;
        for (size_t m = 0; m < sum.size(); m++)
        {
            av[m] = 0;
            for (int j = 0; j < nframes-i; j++)
            {
                av[m] += square(sum[m][j+i]-sum[m][j]);
            }
            avtot += av[m];
        }
        avtot /= sum.size();
        /* Convert to SI for the viscosity */
        double fac = (V*NANO*NANO*NANO*PICO*1e10)/(2*BOLTZMANN*T)/(nframes-i);
        if (nullptr != fpi)
        {
            fprintf(fpi.get(), "%10g", i*Dt);
        }
        for (size_t m = 0; (m < sum.size()); m++)
        {
            av[m] = fac*av[m];
            if (nullptr != fpi)
            {
                fprintf(fpi.get(), "  %10g", av[m]);
            }
        }
        avtot *= fac;
        if (nullptr != fpi)
        {
            fprintf(fpi.get(), "  %10g\n", avtot);
        }
        // Normal viscosity graph
        if (nullptr != fp)
        {
            fprintf(fp.get(), "%10g", (i+0.5)*Dt-t0);
        }
        for (unsigned int m = 0; (m < sum.size()); m++)
        {
            gamma[m] = (av[m]-avold[m])/Dt;
            if (nullptr != fp)
            {
                fprintf(fp.get(), "  %10g", gamma[m]);
            }
            avold[m] = av[m];
        }
        gammatot = (avtot-avoldtot)/Dt;
        avoldtot = avtot;
        if (nullptr != fp)
        {
            fprintf(fp.get(), "  %10g\n", gammatot);
        }
    }
    printf("\n");
    printf("Computed viscosity using Einstein relation based on %u pressure terms.\n",
           static_cast<unsigned int>(sum.size()));
    printf("Shear viscosity: %g cP\n", 1000*gammatot);
    please_cite(stdout, "Allen1987a");
}

const char ViscosityInfo::name[] = "viscosity";

const char ViscosityInfo::shortDescription[] = "Compute viscosity in many ways";

EnergyAnalysisModulePointer ViscosityInfo::create()
{
    return EnergyAnalysisModulePointer(new ViscosityModule());
}

} // namespace energyanalysis

} // namespace gmx
