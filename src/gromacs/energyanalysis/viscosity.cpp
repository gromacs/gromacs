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
 * Implements classes in viscosity.h.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_energyanalysis
 */
#include "gmxpre.h"

#include "viscosity.h"

#include <stdio.h>
#include <string.h>

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/correlationfunctions/autocorr.h"
#include "gromacs/correlationfunctions/expfit.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/viewit.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"

#include "energyhandler.h"
#include "energyhelper.h"
#include "energyinfo.h"
#include "simple.h"

namespace gmx
{

namespace energyanalysis
{

class ViscosityModule : public SimpleEnergyModule
{
    private:
        //! Filename for Shear and Bulk viscosity
        std::string              fnVisco_;
        //! Filename for Einstein viscosity
        std::string              fnEinstein_;
        //! Filename for Einstein viscosity integral
        std::string              fnEinsteinIntegral_;
        //! Filename for Pressure Autocorrelation functions
        std::string              fnPressureAcf_;
        //! Names of energy terms corresponding to what is needed in this analysis
        std::vector<std::string> setnm_;
        //! User given volume of the box
        double                   volume_;
        /*! \brief
         * Compute viscosity using Einstein relation, see GROMACS manual.
         * The output is written to two files.
         * \param[in] sum The sum data sets
         * \param[in] V The average volume
         * \param[in] T The average temperature
         */
        void doEinstein(std::vector < std::vector < real> > sum, real V, real T);

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
        virtual void initOptions(Options                           *options,
                                 ICommandLineOptionsModuleSettings *settings);

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

ViscosityModule::ViscosityModule()
{
    volume_ = 0;

    setnm_.push_back("Pres-XX");
    setnm_.push_back("Pres-XY");
    setnm_.push_back("Pres-XZ");
    setnm_.push_back("Pres-YX");
    setnm_.push_back("Pres-YY");
    setnm_.push_back("Pres-YZ");
    setnm_.push_back("Pres-ZX");
    setnm_.push_back("Pres-ZY");
    setnm_.push_back("Pres-ZZ");
    setnm_.push_back("Temperature");
    setnm_.push_back("Volume");
    setnm_.push_back("Pressure");
};

//! Set the options and setting
void ViscosityModule::initOptions(Options                           *options,
                                  ICommandLineOptionsModuleSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] computes shear- and bulk viscosity using a",
        "number of different algorithms. The best value for the",
        "viscosity is obtained from the integral of the pressure",
        "autocorrelation function, which is computed by setting",
        "the option [TT]-pcorr[tt]."
    };
    // Add the descriptive text (program help text) to the options
    settings->setHelpText(desc);

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
    options->addOption(DoubleOption("volume")
                           .store(&volume_)
                           .description("Volume of the computational box (nm^3) in order to be able to use the tool with constant volume simulations."));
    ehelper_->initOptions(options);
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
            if (eName[i] == setnm_[j])
            {
                EnergyTerm et(i, ehelper_->getStoreData(), eName[i], eUnit[i]);
                ehelper_->addEnergyTerm(et);
                bFound[j] = true;
                nfound++;
            }
        }
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
    ehelper_->setStoreData(true);

    return (nfound == vNR);
}

bool ViscosityModule::addAnalysisFrame(t_enxframe *fr)
{
    return SimpleEnergyModule::addAnalysisFrame(fr);
}

//! Estimate omega, the period of rapid fluctuations in the pressure acf.
static void estimateOmega(double Dt, std::vector<real> data,
                          double fitparam[])
{
    std::vector<real> maxima;

    for (unsigned int i = 1; (i < data.size()) && (maxima.size() < 10); i++)
    {
        if ((data[i] > data[i-1]) && (data[i] > data[i+1]))
        {
            // Found a maximum
            maxima.push_back(Dt*i);
        }
    }
    // Estimate period from difference between first and last top
    if (maxima.size() > 2)
    {
        fitparam[1] = (maxima.size()-1)*2*M_PI/(maxima.back() - maxima.front());
        printf("Estimating omega to be %g (1/ps)\n", fitparam[1]);
    }
    else
    {
        // Default return 1/ps
        fitparam[1] = 1;
    }
    fitparam[0] = 0.5; /* C */
    fitparam[2] = 0.1; /* (ps) tau_Kf */
    fitparam[3] = 0.8; /* beta_f */
    fitparam[4] = 0.5; /* (ps) tau_Ks */
    fitparam[5] = 0.6; /* beta_s */
}

//! Utility function to dump parameters to a file.
static void print_params(FILE *fp,
                         const char *header, double params[])
{
    fprintf(fp, "%s:\n", header);
    fprintf(fp, "C      = %12g    Omega  = %12g 1/ps\n", params[0], params[1]);
    fprintf(fp, "tau_Kf = %12g ps beta_f = %12g\n", params[2], params[3]);
    fprintf(fp, "tau_Ks = %12g ps beta_s = %12g\n", params[4], params[5]);
}

//! Utility function fit the pressure acf to the function from Fanourgakis.
static void fit_pressure_acf(int nout, double Dt, std::vector<real> P,
                             double fit[6], output_env_t oenv)
{
    real *dy;

    snew(dy, 1+nout);
    estimateOmega(Dt, P, fit);
    print_params(stdout, "Parameters before", fit);
    double period = 1/(2*M_PI*fit[1]);
    double t_end  = std::min(std::max(100*Dt, 50*period), nout*Dt);
    printf("Using t_end = %g for fitting to pressure ACF.\n", t_end);
    (void) do_lmfit(nout, &P[0], dy, Dt, NULL,
                    0, t_end, oenv,
                    FALSE, effnPRES, fit, 0, NULL);
    print_params(stdout, "Parameters after", fit);
    sfree(dy);
}

bool ViscosityModule::finalizeAnalysis()
{
    double             eAver[vNR], eStddev[vNR];
    EnergyTermIterator etxyz[vNR];
    double             Dt;
    gmx_int64_t        nsteps, nframes;
    int                nout;
    const char        *leg[2] = { "Shear", "Bulk" };
    int                nleg   = sizeof(leg)/sizeof(leg[0]);

    for (int i = 0; (i < vNR); i++)
    {
        if (!ehelper_->getEnergyTerm(setnm_[i], &eAver[i], &eStddev[i]))
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
                    return false;
                }
            }
        }
        else
        {
            etxyz[i] = ehelper_->etSearch(setnm_[i]);
            if (NULL != debug)
            {
                fprintf(debug, "Found %s for %s\n",
                        etxyz[i]->getEterm().c_str(), setnm_[i].c_str());
            }
        }
    }

    nsteps  = ehelper_->etBegin()->nSteps();
    nframes = ehelper_->etBegin()->nEnergy();
    if ((nsteps < 2) || (nframes < 2))
    {
        char buf1[64], buf2[64];
        fprintf(stderr, "Not enough data: %s steps and %s frames.\n",
                gmx_step_str(nsteps, buf1),
                gmx_step_str(nframes, buf2));
        if (!ehelper_->getStoreData())
        {
            fprintf(stderr, "You forgot to turn on storing data\n");
        }
        return false;
    }
    Dt      = ehelper_->etBegin()->timeSpan() / (nframes-1);

    // Do viscosity calculation using the Einstein method
    if ((fnEinstein_.size() > 0) || (fnEinsteinIntegral_.size() > 0))
    {
        std::vector < std::vector < real> > eneint;
        std::vector<real>                   zero;
        zero.push_back(0.0);
        for (int i = 0; (i < 5); i++)
        {
            eneint.push_back(zero);
        }
        // Assume all components are stored equally often
        int nSum = etxyz[vPresXY]->searchEF(0)->getNsum();
        if (0 >= nSum)
        {
            gmx_fatal(FARGS, "nSum = %d", nSum);
        }
        double fac = Dt/nSum;
        for (gmx_int64_t i = 0; (i < nframes-1); i++)
        {
            eneint[0].push_back(eneint[0].back() +
                                fac*(etxyz[vPresXY]->searchEF(i)->getEsum() -
                                     eAver[vPresXY]));
            eneint[1].push_back(eneint[1].back() +
                                fac*(etxyz[vPresXZ]->searchEF(i)->getEsum() -
                                     eAver[vPresXZ]));
            eneint[2].push_back(eneint[2].back() +
                                fac*(etxyz[vPresYZ]->searchEF(i)->getEsum() -
                                     eAver[vPresYZ]));
            eneint[3].push_back(eneint[2].back() +
                                0.5*fac*(etxyz[vPresXX]->searchEF(i)->getEsum() -
                                         etxyz[vPresYY]->searchEF(i)->getEsum()));
            eneint[4].push_back(eneint[2].back() +
                                0.5*fac*(etxyz[vPresYY]->searchEF(i)->getEsum() -
                                         etxyz[vPresZZ]->searchEF(i)->getEsum()));
        }
        doEinstein(eneint, eAver[vVolume], eAver[vTemperature]);
    }

    // Do viscosity calculation using a fit to the autocorrelation function
    std::vector<real> shearP[5], bulkP;
    for (gmx_int64_t i = 0; (i < nframes); i++)
    {
        /* See Fanourgakis J. Phys. Chem. A 2012, 116, 2564-2570 */
        shearP[0].push_back(0.5*(etxyz[vPresXX]->searchEF(i)->getE() -
                                 etxyz[vPresYY]->searchEF(i)->getE()));
        shearP[1].push_back(0.5*(etxyz[vPresYY]->searchEF(i)->getE() -
                                 etxyz[vPresZZ]->searchEF(i)->getE()));
        shearP[2].push_back(etxyz[vPresXY]->searchEF(i)->getE() - eAver[vPresXY]);
        shearP[3].push_back(etxyz[vPresYZ]->searchEF(i)->getE() - eAver[vPresYZ]);
        shearP[4].push_back(etxyz[vPresZX]->searchEF(i)->getE() - eAver[vPresZX]);
        // Average pressure fluctuations for Bulk viscosity
        bulkP.push_back(etxyz[vPressure]->searchEF(i)->getE() - eAver[vPressure]);
    }

    nout = nframes/2+1;
    // Todo t_pargs not used but need to be initilized
    int       n        = 0;
    t_pargs  *tempArgs = add_acf_pargs(&n, NULL); // init args for autocorr
    real     *sp[5], *bp;
    for (int i = 0; (i < 5); i++)
    {
        sp[i] = &(shearP[i][0]);
    }
    bp = &(bulkP[0]);
    low_do_autocorr(NULL, ehelper_->getOutputEnvironment(), leg[0], nframes,
                    5, nout, sp, Dt,
                    eacNormal, 1, TRUE, FALSE, FALSE, 0.0, 0.0, 0);

    /* Now for bulk viscosity */
    low_do_autocorr(NULL, ehelper_->getOutputEnvironment(), leg[1], nframes,
                    1, nout, &bp, Dt,
                    eacNormal, 1, TRUE, FALSE, FALSE, 0.0, 0.0, 0);
    sfree(tempArgs);

    double factor = (eAver[vVolume]*1e-26/(BOLTZMANN*eAver[vTemperature]));
    if ((nout > 1) && (getFnPressureAcf().size() > 0))
    {
        double bulk0, bulkfit[6];
        double shear0, shearfit[6];

        /* Run the analysis based on the bulk pressure flutuations */
        bulk0 = bulkP[0];
        for (int i = 0; (i < nout); i++)
        {
            bulkP[i] /= bulk0;
        }
        fit_pressure_acf(nout, Dt, bulkP, bulkfit,
                         ehelper_->getOutputEnvironment());

        /* Run the shear analysis */
        shear0 = 1; // In case nout == 0
        for (int i = 0; (i < nout); i++)
        {
            real sum = 0.2*(shearP[0][i]+shearP[1][i]+shearP[2][i]+shearP[3][i]+shearP[4][i]);
            if (i == 0)
            {
                shear0 = sum;
            }
            shearP[0][i] = sum/shear0;
        }
        fit_pressure_acf(nout, Dt, shearP[0], shearfit,
                         ehelper_->getOutputEnvironment());

        FILE *fp = xvgropen(getFnPressureAcf().c_str(),
                            "Pressure ACF", "Time (ps)",
                            "(units)", ehelper_->getOutputEnvironment());
        if (NULL != fp)
        {
            double      bInt, sInt;
            const char *pleg[] = { "Bulk Pacf", "Shear Pacf", "Fit to Bulk", "Fit to Shear" };
            int         npleg  = sizeof(pleg)/sizeof(pleg[0]);
            xvgr_legend(fp, npleg, pleg, ehelper_->getOutputEnvironment());
            bInt = sInt = 0;
            for (int i = 0; (i < nout); i++)
            {
                double x         = i*Dt;
                double bfit      = bulk0*fit_function(effnPRES, bulkfit, x);
                double sfit      = shear0*fit_function(effnPRES, shearfit, x);
                double trapezoid = (i == 0 || i == nout-1) ? 0.5 : 1.0;

                fprintf(fp, "%10g  %10g  %10g  %10g  %10g\n",
                        x,
                        bulk0*bulkP[i], shear0*shearP[0][i],
                        bfit, sfit);
                bInt += x*bfit*trapezoid;
                sInt += x*sfit*trapezoid;
            }
            xvgrclose(fp);
            printf("Computed viscosity based on integral of fit to pressure autocorrelation.\n");
            printf("Bulk viscosity:  %g cP\n", 1000*bInt*factor);
            printf("Shear viscosity: %g cP\n", 1000*sInt*factor);
            please_cite(stdout, "Fanourgakis2012a");
        }
    }

    if (0)
    {
        FILE * fp = NULL;
        if (getOutputFile().size() > 0)
        {
            fp = xvgropen(getOutputFile().c_str(), "Viscosity", "Time (ps)",
                          "\\8h\\4 (cp)", ehelper_->getOutputEnvironment());
            if (NULL != fp)
            {
                xvgr_legend(fp, nleg, leg, ehelper_->getOutputEnvironment());
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
            double dShear = 0.5*(shearP[0][i-1]  + shearP[0][i])*factor*Dt;
            double dBulk  = 0.5*(bulkP[i-1] + bulkP[i])*factor*Dt;
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
        if (nAv > 1)
        {
            avShear /= nAv;
            avBulk  /= nAv;
            printf("Shear viscosity %g +/- %g\n", avShear,
                   sqrt(avS2/nAv - avShear*avShear));
            printf("Bulk viscosity  %g +/- %g\n", avBulk,
                   sqrt(avB2/nAv - avBulk*avBulk));
        }
        else
        {
            fprintf(stderr, "Too few data points to compute viscosities\n");
            return false;
        }
        if (NULL != fp)
        {
            xvgrclose(fp);
        }
    }
    return true;
}

void ViscosityModule::doEinstein(std::vector < std::vector < real> > sum,
                                 real V, real T)
{
    std::vector<double> av, avold, gamma;
    double              avoldtot, gammatot;
    double              fac, di;
    int                 i, nf4;

    gmx_int64_t         nframes = ehelper_->etBegin()->nEnergy();
    double              t0      = ehelper_->etBegin()->timeBegin();
    double              Dt      = ehelper_->etBegin()->timeSpan() / (nframes-1);

    if (nframes < 4)
    {
        fprintf(stderr, "Number of frames %d too few for computing Einstein viscosity.\n", (int)nframes);
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
    for (unsigned int m = 0; m < sum.size(); m++)
    {
        av.push_back(0);
        avold.push_back(0);
        gamma.push_back(0);
    }
    avoldtot = gammatot = 0;
    for (i = 1; i < nf4; i++)
    {
        double avtot = 0;
        for (unsigned int m = 0; m < sum.size(); m++)
        {
            av[m] = 0;
            for (int j = 0; j < nframes-i; j++)
            {
                di   = sqr(sum[m][j+i]-sum[m][j]);

                av[m] += di;
                avtot += di/sum.size();
            }
        }
        /* Convert to SI for the viscosity */
        fac = (V*NANO*NANO*NANO*PICO*1e10)/(2*BOLTZMANN*T)/(nframes-i);
        if (NULL != fpi)
        {
            fprintf(fpi, "%10g", i*Dt);
        }
        for (unsigned int m = 0; (m < sum.size()); m++)
        {
            av[m] = fac*av[m];
            if (NULL != fpi)
            {
                fprintf(fpi, "  %10g", av[m]);
            }
        }
        avtot *= fac;
        if (NULL != fpi)
        {
            fprintf(fpi, "  %10g\n", avtot);
        }
        if (NULL != fp)
        {
            fprintf(fp, "%10g", (i+0.5)*Dt-t0);
        }
        for (unsigned int m = 0; (m < sum.size()); m++)
        {
            gamma[m] = (av[m]-avold[m])/Dt;
            if (NULL != fp)
            {
                fprintf(fp, "  %10g", gamma[m]);
            }
            avold[m] = av[m];
        }
        gammatot = (avtot-avoldtot)/Dt;
        avoldtot = avtot;
        if (NULL != fp)
        {
            fprintf(fp, "  %10g\n", gammatot);
        }
    }
    printf("Computed viscosity using Einstein relation based on %u pressure terms.\n",
           static_cast<unsigned int>(sum.size()));
    printf("Shear viscosity: %g cP\n", 1000*gammatot);
    if (NULL != fpi)
    {
        xvgrclose(fpi);
    }
    if (NULL != fp)
    {
        xvgrclose(fp);
    }
}

}   // namespace energyanalysis

const char ViscosityInfo::name[] = "viscosity";

const char ViscosityInfo::shortDescription[] = "Compute viscosity in many ways";

ICommandLineOptionsModule *ViscosityInfo::create()
{
    energyanalysis::EnergyAnalysisModulePointer eamp(new energyanalysis::ViscosityModule);
    return new EnergyInfo(eamp);
}

} // namespace gmx
