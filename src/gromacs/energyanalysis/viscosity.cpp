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
#include <stdio.h>
#include <string.h>

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/correlationfunctions/autocorr.h"
#include "gromacs/correlationfunctions/expfit.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/viewit.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/units.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/smalloc.h"

#include "energyhandler.h"
#include "energyhelper.h"
#include "energyinfo.h"
#include "select.h"
#include "simple.h"
#include "viscosity.h"

namespace gmx
{

class ViscosityModule : public SimpleEnergyModule
{
    private:
        //! Filename for Shear and Bulk viscosity
        std::string fnVisco_;
        //! Filename for Einstein viscosity
        std::string fnEinstein_;
        //! Filename for Einstein viscosity integral
        std::string fnEinsteinIntegral_;
        //! Filename for Pressure Autocorrelation functions
        std::string fnPressureAcf_;
        /*! \brief
         * Compute viscosity using Einstein relation, see GROMACS manual.
         * The output is written to two files.
         * \param[in] nsets Number of sum data sets (should be three)
         * \param[in] sum The sum data sets
         * \param[in] V The average volume
         * \param[in] T The average temperature
         */
        void doEinstein(int nsets, real **sum, real V, real T);
        //! Energy helper class for low level stuff
        EnergyHelper  *ehelper_;

    public:
        //! Constructor
        ViscosityModule();

        //! Destructor must be virtual since there are other virtual functions
        virtual ~ViscosityModule() {};

        //! Set output file name for Einstein viscosity
        void setFnEinstein(std::string fn) { fnEinstein_ = fn; }

        //! Return the output file name for Einstein viscosity
        std::string getFnEinstein() { return fnEinstein_; }

        //! Set output file name for Einstein integral
        void setFnEinsteinIntegral(std::string fn) { fnEinsteinIntegral_ = fn; }

        //! Return the output file name for Einstein viscosity integral
        std::string getFnEinsteinIntegral() { return fnEinsteinIntegral_; }

        //! Set output file name for the Pressure autocorrelation
        void setFnPressureAcf(std::string fn) { fnPressureAcf_ = fn; }

        //! Return the output file name for the Pressure autocorrelation
        std::string getFnPressureAcf() { return fnPressureAcf_; }

        //! Pass the output environment on for printing etc.
        virtual void setOutputEnvironment(output_env_t oenv)
        {
            ehelper_->setOutputEnvironment(oenv);
        }

        //! Initiate the command line options
        virtual void initOptions(Options *options);

        /*! \brief
         * Does the initiation of the analysis of the file
         * \param[in] eName Names of the energy terms etc.
         * \param[in] eUnit Units of the energy terms etc.
         * \return true if OK, false otherwise
         */
        virtual bool initAnalysis(std::vector<std::string> eName,
                                  std::vector<std::string> eUnit);

        //! Analyse one frame and stores the results in memory
        virtual bool addAnalysisFrame(t_enxframe *fr);

        //! Finalize reading
        virtual bool finalizeAnalysis();
};

enum {
    vPresXX      = 0, vPresXY = 1, vPresXZ = 2,
    vPresYX      = 3, vPresYY = 4, vPresYZ = 5,
    vPresZX      = 6, vPresZY = 7, vPresZZ = 8,
    vTemperature = 9, vVolume = 10, vPressure = 11,
    vNR          = 12
};

/*! \brief
 * Strings that should appear in the energy file in order to
 * do viscosity analysis
 */
static const char *setnm[] = {
    "Pres-XX",     "Pres-XY", "Pres-XZ",
    "Pres-YX",     "Pres-YY", "Pres-YZ",
    "Pres-ZX",     "Pres-ZY", "Pres-ZZ",
    "Temperature", "Volume",  "Pressure"
};

ViscosityModule::ViscosityModule()
{
    ehelper_   = new(EnergyHelper);
}

//! Set the options and setting
void ViscosityModule::initOptions(Options *options)
{
    static const char *const desc[] = {
        "[THISMODULE] computes shear- and bulk viscosity using a",
        "number of different algorithms. The best value for the",
        "viscosity is obtained from the integral of the pressure",
        "autocorrelation function, which is computed by setting",
        "the option [TT]-pcorr[tt]."
    };
    // Add the descriptive text (program help text) to the options
    options->setDescription(desc);

    // Add option for optional output files
    options->addOption(FileNameOption("vis").filetype(eftPlot).outputFile()
                           .store(&fnVisco_).defaultBasename("visco")
                           .description("Computed shear and bulk viscosity"));
    options->addOption(FileNameOption("evis").filetype(eftPlot).outputFile()
                           .store(&fnEinstein_).defaultBasename("evisco")
                           .description("Shear viscosity from Einstein relation"));
    options->addOption(FileNameOption("eivis").filetype(eftPlot).outputFile()
                           .store(&fnEinsteinIntegral_).defaultBasename("eviscointegral")
                           .description("Integral of shear viscosity from Einstein relation"));
    options->addOption(FileNameOption("pcorr").filetype(eftPlot).outputFile()
                           .store(&fnPressureAcf_).defaultBasename("presscorr")
                           .description("Pressure autocorrelation"));
}

bool ViscosityModule::initAnalysis(std::vector<std::string> eName,
                                   std::vector<std::string> eUnit)
{
    bool bFound[vNR];
    int  nfound = 0;
    for (int j = 0; (j < vNR); j++)
    {
        bFound[j] = false;
    }
    for (unsigned int i = 0; (i < eName.size()); i++)
    {
        for (int j = 0; (j < vNR); j++)
        {
            if (gmx_strcasecmp(eName[i].c_str(), setnm[j]) == 0)
            {
                EnergyTerm et(i, eName[i], eUnit[i]);
                ehelper_->addEnergyTerm(et);
                bFound[j] = true;
                nfound++;
            }
        }
    }
    ehelper_->setStoreData(true);
    for (int j = 0; (j < vNR); j++)
    {
        if (!bFound[j])
        {
            fprintf(stderr, "Energy file does not contain %s. Can not compute viscosity.\n",
                    setnm[j]);
        }
    }
    return (nfound == vNR);
}

bool ViscosityModule::addAnalysisFrame(t_enxframe *fr)
{
    return SimpleEnergyModule::addAnalysisFrame(fr);
}

bool ViscosityModule::finalizeAnalysis()
{
    double      eAver[vNR], eStddev[vNR];
    double      Dt;
    gmx_int64_t nsteps, nframes;
    int         nout;
    gmx_bool    bVerbose;

    bVerbose = output_env_get_verbosity(ehelper_->getOutputEnvironment()) > 0;
    for (int i = 0; (i < vNR); i++)
    {
        if (!ehelper_->getEnergyTerm(setnm[i], &eAver[i], &eStddev[i]))
        {
            fprintf(stderr, "Could not extract %s from energy file.\n",
                    setnm[i]);
            return false;
        }
    }
    nsteps  = ehelper_->etBegin()->nSteps();
    if (nsteps < 2)
    {
        char buf[256];
        fprintf(stderr, "Only %s steps, that's not enough.\n",
                gmx_step_str(nsteps, buf));
    }
    nframes = ehelper_->etBegin()->nEnergy();
    Dt      = ehelper_->etBegin()->timeSpan() / (nframes-1);

    // Symmetrise tensor! (and store in first three elements)
#define NPCOMP 6
    real      *eneset[NPCOMP];
    for (int i = 0; i < NPCOMP; i++)
    {
        snew(eneset[i], nframes);
    }
    real      **eneint;
    snew(eneint, 3);
    for (int i = 0; i < 3; i++)
    {
        snew(eneint[i], nframes);
    }
    EnergyTermIterator etxyz[vNR];
    for (int i = 0; (i < vNR); i++)
    {
        etxyz[i] = ehelper_->etSearch(setnm[i]);
        printf("Found %s for %s\n", etxyz[i]->getEterm().c_str(), setnm[i]);
    }
    for (gmx_int64_t i = 0; (i < nframes); i++)
    {
        /* See Fanourgakis J. Phys. Chem. A 2012, 116, 2564-2570 */
        eneset[0][i] = 0.5*(etxyz[vPresXX]->searchEF(i)->getE() -
                            etxyz[vPresYY]->searchEF(i)->getE());
        eneset[1][i] = 0.5*(etxyz[vPresYY]->searchEF(i)->getE() -
                            etxyz[vPresZZ]->searchEF(i)->getE());
        eneset[2][i] = etxyz[vPresXY]->searchEF(i)->getE() - eAver[vPresXY];
        eneset[3][i] = etxyz[vPresYZ]->searchEF(i)->getE() - eAver[vPresYZ];
        eneset[4][i] = etxyz[vPresZX]->searchEF(i)->getE() - eAver[vPresZX];
        // Average pressure fluctuations
        eneset[5][i] = etxyz[vPressure]->searchEF(i)->getE() - eAver[vPressure];
    }
    for (gmx_int64_t i = 0; (i < nframes-1); i++)
    {
        // Assume all components are store equally often
        int nSum = etxyz[vPresXY]->searchEF(i)->getNsum();
        if (0 >= nSum)
        {
            gmx_fatal(FARGS, "nSum = %d", nSum);
        }
        double fac = 0.5*Dt/nSum;

        eneint[0][i+1] = eneint[0][i] +
            fac*(etxyz[vPresXY]->searchEF(i)->getEsum()+
                 etxyz[vPresYX]->searchEF(i)->getEsum());
        eneint[1][i+1] = eneint[1][i] +
            fac*(etxyz[vPresXZ]->searchEF(i)->getEsum()+
                 etxyz[vPresZX]->searchEF(i)->getEsum());
        eneint[2][i+1] = eneint[2][i] +
            fac*(etxyz[vPresYZ]->searchEF(i)->getEsum()+
                 etxyz[vPresZY]->searchEF(i)->getEsum());
    }
    doEinstein(3, eneint,
               etxyz[vVolume]->average(),
               etxyz[vTemperature]->average());

    /* Do it for shear viscosity */
#define NLEG 2
    const char* leg[NLEG] = { "Shear", "Bulk" };
    nout = nframes/2+1;
    //Todo t_pargs not used but need to be initilized
    int       n        = 0;
    t_pargs * tempArgs = add_acf_pargs(&n, NULL); // init args for autocorr

    low_do_autocorr(NULL, ehelper_->getOutputEnvironment(), leg[0], nframes, NPCOMP-1,
                    nout, eneset, Dt,
                    eacNormal, 1, TRUE, FALSE, bVerbose, 0.0, 0.0, 0);

    /* Now for bulk viscosity */
    low_do_autocorr(NULL, ehelper_->getOutputEnvironment(), leg[1], nframes, 1,
                    nout, &(eneset[NPCOMP-1]), Dt,
                    eacNormal, 1, TRUE, FALSE, bVerbose, 0.0, 0.0, 0);

    free(tempArgs);
    double factor = (eAver[vVolume]*1e-26/(BOLTZMANN*eAver[vTemperature]));
    if (getFnPressureAcf().size() > 0)
    {
        real   Pintegral;
        int    fix = 0;
        double bulk0, shear0;
        double bulkfit[6]  = {
            0.5, /* C */ 10, /* (1/ps) omega */ 4,                    /* (ps) tau_Kf */
            1, /* beta_f */ 0.5, /* (ps) tau_Ks */ 1                  /* beta_s */
        };
        double shearfit[6] = {
            0.5, /* C */ 10, /* (1/ps) omega */ 4,                    /* (ps) tau_Kf */
            1, /* beta_f */ 0.5, /* (ps) tau_Ks */ 1                  /* beta_s */
        };

        /* Run the analysis based on the bulk pressure flutuations */
        bulk0 = eneset[NPCOMP-1][0];
        for (int i = 0; (i < nout); i++)
        {
            eneset[NPCOMP-1][i] /= bulk0;
        }
        Pintegral = bulk0*do_lmfit(nout, eneset[NPCOMP-1], NULL, Dt, NULL, 0, nout*Dt,
                                   ehelper_->getOutputEnvironment(),
                                   FALSE, effnPRES, bulkfit, fix);
        printf("Integral of fit to pressure fluctuation ACF: %g\n", Pintegral);

        /* Run the shear analysis */
        for (int i = 0; (i < nout); i++)
        {
            real sum = 0.2*(eneset[0][i]+eneset[1][i]+eneset[2][i]+eneset[3][i]+eneset[4][i]);
            if (i == 0)
            {
                shear0 = sum;
            }
            eneset[0][i] = sum/shear0;
        }
        Pintegral = shear0*do_lmfit(nout, eneset[0], NULL, Dt, NULL, 0, nout*Dt,
                                    ehelper_->getOutputEnvironment(),
                                    FALSE, effnPRES, shearfit, fix);
        printf("Integral of fit to off-diagonal pressure ACF: %g\n", Pintegral);

        FILE *fp = xvgropen(getFnPressureAcf().c_str(),
                            "Pressure ACF", "Time (ps)",
                            "(units)", ehelper_->getOutputEnvironment());
        if (NULL != fp)
        {
            double bInt, sInt;
#define NPLEG 4
            const char *pleg[NPLEG] = { "Bulk Pacf", "Shear Pacf", "Fit to Bulk", "Fit to Shear" };
            xvgr_legend(fp, NPLEG, pleg, ehelper_->getOutputEnvironment());
#undef NPLEG
            bInt = sInt = 0;
            for (int i = 0; (i < nout); i++)
            {
                double x         = i*Dt;
                double bfit      = bulk0*fit_function(effnPRES, bulkfit, x);
                double sfit      = shear0*fit_function(effnPRES, shearfit, x);
                double trapezoid = (i == 0 || i == nout-1) ? 0.5 : 1.0;

                fprintf(fp, "%10g  %10g  %10g  %10g  %10g\n",
                        x,
                        bulk0*eneset[NPCOMP-1][i], bfit,
                        shear0*eneset[0][i], sfit);
                bInt += x*bfit*trapezoid;
                sInt += x*sfit*trapezoid;
            }
            xvgrclose(fp);
            printf("Bulk viscosity:  %g Pa s\n",  bInt*factor);
            printf("Shear viscosity: %g Pa s\n", sInt*factor);
        }
    }
    FILE * fp = NULL;
    if (getOutputFile().size() > 0)
    {
        fp = xvgropen(getOutputFile().c_str(), "Viscosity", "Time (ps)",
                      "\\8h\\4 (cp)", ehelper_->getOutputEnvironment());
        if (NULL != fp)
        {
            xvgr_legend(fp, NLEG, leg, ehelper_->getOutputEnvironment());
        }
    }
    /* Use trapezium rule for integration */
    double intShear = 0;
    double intBulk  = 0;
    double avShear  = 0, avS2 = 0;
    double avBulk   = 0, avB2 = 0;
    int    nAv      = 0;
    for (int i = 1; (i < nout); i++)
    {
        double dShear = 0.5*(eneset[0][i-1]  + eneset[0][i])*factor*Dt;
        double dBulk  = 0.5*(eneset[3][i-1] + eneset[3][i])*factor*Dt;
        intShear += dShear;
        intBulk  += dBulk;
        if (i >= nout/2)
        {
            avShear += dShear;
            avS2    += dShear*dShear;
            avBulk  += dBulk;
            avB2    += dBulk*dBulk;
            nAv++;
        }
        if (NULL != fp)
        {
            fprintf(fp, "%10g  %10g  %10g\n", (i*Dt), intShear, intBulk);
        }
    }
    for (int i = 0; (i < NPCOMP); i++)
    {
        sfree(eneset[i]);
    }
    avShear /= nAv;
    avBulk  /= nAv;
    printf("Shear viscosity %g +/- %g\n", avShear,
           sqrt(avS2/nAv - avShear*avShear));
    printf("Bulk viscosity  %g +/- %g\n", avBulk,
           sqrt(avB2/nAv - avBulk*avBulk));
    for (int i = 0; (i < 3); i++)
    {
        sfree(eneint[i]);
    }
    sfree(eneint);
    if (NULL != fp)
    {
        xvgrclose(fp);
    }

    return true;
}


void ViscosityModule::doEinstein(int nsets, real **sum,
                                 real V, real T)
{
#define NSETS 3
    double      av[1+NSETS], avold[1+NSETS], gamma[1+NSETS];
    double      fac, di;
    int         i, j, m, nf4;

    gmx_int64_t nframes = ehelper_->etBegin()->nEnergy();
    double      t0      = ehelper_->etBegin()->timeBegin();
    double      Dt      = ehelper_->etBegin()->timeSpan() / (nframes-1);

    if (nframes < 4)
    {
        fprintf(stderr, "Number of frames %d too few for computing Einstein viscosity.\n", (int)nframes);
        return;
    }
    if (nsets != NSETS)
    {
        fprintf(stderr, "Number of sets passed to doEinstein should be exactly %d.\n", NSETS);
        return;
    }
    nf4 = nframes/4+1;

    std::string fni = getFnEinsteinIntegral();
    FILE       *fpi = NULL;
    if (fni.size() > 0)
    {
        fpi = xvgropen(fni.c_str(), "Shear viscosity integral",
                       "Time (ps)", "(kg m\\S-1\\N s\\S-1\\N ps)",
                       ehelper_->getOutputEnvironment());
    }
    std::string fn = getFnEinstein();
    FILE       *fp = NULL;
    if (fn.size() > 0)
    {
        fp = xvgropen(fn.c_str(), "Shear viscosity using Einstein relation",
                      "Time (ps)", "(kg m\\S-1\\N s\\S-1\\N)",
                      ehelper_->getOutputEnvironment());
    }
    for (m = 0; m <= NSETS; m++)
    {
        av[m] = avold[m] = gamma[m] = 0;
    }
    for (i = 1; i < nf4; i++)
    {
        for (m = 0; m <= NSETS; m++)
        {
            av[m] = 0;
        }
        for (j = 0; j < nframes-i; j++)
        {
            for (m = 0; m < NSETS; m++)
            {
                di   = sqr(sum[m][j+i]-sum[m][j]);

                av[m]     += di;
                av[NSETS] += di/NSETS;
            }
        }
        /* Convert to SI for the viscosity */
        fac = (V*NANO*NANO*NANO*PICO*1e10)/(2*BOLTZMANN*T)/(nframes-i);
        if (NULL != fpi)
        {
            fprintf(fpi, "%10g", i*Dt);
        }
        for (m = 0; (m <= NSETS); m++)
        {
            av[m] = fac*av[m];
            if (NULL != fpi)
            {
                fprintf(fpi, "  %10g", av[m]);
            }
        }
        if (NULL != fpi)
        {
            fprintf(fpi, "\n");
        }
        if (NULL != fp)
        {
            fprintf(fp, "%10g", (i+0.5)*Dt-t0);
        }
        for (m = 0; (m <= NSETS); m++)
        {
            gamma[m] = (av[m]-avold[m])/Dt;
            if (NULL != fp)
            {
                fprintf(fp, "  %10g", gamma[m]);
            }
            avold[m] = av[m];
        }
        if (NULL != fp)
        {
            fprintf(fp, "\n");
        }
    }
    printf("Viscosity determined from Einstein relation %.3f cP\n",
           1000*gamma[NSETS]);
    if (NULL != fpi)
    {
        xvgrclose(fpi);
    }
    if (NULL != fp)
    {
        xvgrclose(fp);
    }
#undef NSETS
}

const char ViscosityInfo::name[] = "viscosity";

const char ViscosityInfo::shortDescription[] = "Compute viscosity in many ways";

CommandLineOptionsModuleInterface *ViscosityInfo::create()
{
    EnergyAnalysisModulePointer eamp(new ViscosityModule);
    return (CommandLineOptionsModuleInterface *)(new EnergyInfo(eamp));
}

}
