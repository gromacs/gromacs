/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016, by the GROMACS development team, led by
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
 * Implements classes in fluctprops.h.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_energyanalysis
 */
#include "gmxpre.h"

#include "fluctprops.h"

#include <cstdio>
#include <cstring>

#include <memory>

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/energyanalysis/analysismodule.h"
#include "gromacs/energyanalysis/energytermcontainer.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/unique_cptr.h"

namespace gmx
{

namespace energyanalysis
{

class FluctPropsModule : public EnergyAnalysisModule
{
    public:
        //! Constructor
        FluctPropsModule();

        /*! \brief
         * Add a term to the data structure and return the energy file index
         * \param[in] eNU   The names and units of energy terms
         * \param[in] term  The name of the term to search for
         * \return The energy file index for this term or INT_MAX if not succesfull
         */
        unsigned int addTerm(const std::vector<EnergyNameUnit> &eNU,
                             const char                        *term);

        virtual void setOutputEnvironment(OutputEnvPointer &&oenv)
        {
            ehelper_.setOutputEnvironment(std::move(oenv));
        }

        void initOptions(IOptionsContainer                 *options,
                         ICommandLineOptionsModuleSettings *settings);

        virtual void initAnalysis(const std::vector<EnergyNameUnit> &eNU);

        virtual void analyzeFrame(t_enxframe *fr);

        virtual void finalizeAnalysis();

        //! View the output file(s)
        void viewOutput();

    private:
        //! Fluctuation convergence output (typically a xvg file)
        std::string                 fnFluctConv_;

        //! File pointer for storing flucutuations
        gmx::unique_cptr<FILE, xvgrclose> fc_;

        //! Index to Enthalpy (not found if < 0)
        unsigned int                iEnthalpy_;

        //! Index to Volume (not found if < 0)
        unsigned int                iVolume_;

        //! Index to Temperature (not found if < 0)
        unsigned int                iTemp_;

        //! Index to Total Energy (not found if < 0)
        unsigned int                iTotal_;

        //! Index to Enthalpy*Volume (not found if < 0)
        unsigned int                iHV_;

        //! Should we correct for the drift in the data?
        bool                        driftCorr_;

        //! Energy helper class for low level stuff
        EnergyTermContainer         ehelper_;

        //! Special EnergyTerm for HV
        std::unique_ptr<EnergyTerm> energyTermHV_;
};

//! String describing the product of Enthalpy and Volume
static const char *HVstring = "HV";

FluctPropsModule::FluctPropsModule() : fc_(nullptr), iEnthalpy_(INT_MAX), iVolume_(INT_MAX),
                                       iTemp_(INT_MAX), iTotal_(INT_MAX), iHV_(INT_MAX),
                                       driftCorr_(false), energyTermHV_(nullptr)
{}

unsigned int FluctPropsModule::addTerm(const std::vector<EnergyNameUnit> &eNU,
                                       const char                        *term)
{
    unsigned int i = 0;
    for (std::vector<EnergyNameUnit>::const_iterator enu = eNU.begin(); (enu < eNU.end()); ++enu, ++i)
    {
        if (enu->energyName == term)
        {
            EnergyTerm et(i, ehelper_.storeData(), term, enu->energyUnit);
            ehelper_.addEnergyTerm(et);
            return i;
        }
    }
    return INT_MAX;
}

void FluctPropsModule::initOptions(IOptionsContainer                 *options,
                                   ICommandLineOptionsModuleSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] will compute the following properties:",
        "",
        "===============================  ===================",
        "Property                         Energy terms needed",
        "===============================  ===================",
        "Heat capacity C[SUB]p[sub] (NPT sims):    Enthalpy, Temp",
        "Heat capacity C[SUB]v[sub] (NVT sims):    Etot, Temp",
        "Thermal expansion coeff. (NPT):  Enthalpy, Vol, Temp",
        "Isothermal compressibility:      Vol, Temp",
        "Adiabatic bulk modulus:          Vol, Temp",
        "===============================  ===================",
        "",
        "You always need to set the number of molecules [TT]-nmol[tt].",
        "The C[SUB]p[sub]/C[SUB]v[sub] computations do [BB]not[bb] include any corrections",
        "for quantum effects. Use the [gmx-dos] program if you need that (and you do).",
    };
    settings->setHelpText(desc);

    options->addOption(BooleanOption("driftcorrections")
                           .store(&driftCorr_)
                           .description("The drift in the observables will be subtracted before computing the fluctuation properties"));
    // Add option for optional output files
    options->addOption(FileNameOption("convergence")
                           .filetype(eftPlot)
                           .outputFile()
                           .store(&fnFluctConv_)
                           .defaultBasename("fluct_conv")
                           .description("Convergence of fluctuation properties by reporting the standard deviation"));
    ehelper_.initOptions(options);
}

void FluctPropsModule::initAnalysis(const std::vector<EnergyNameUnit> &eNU)
{
    // Select which energies to use
    iTemp_     = addTerm(eNU, "Temperature");
    iVolume_   = addTerm(eNU, "Volume");
    iTotal_    = addTerm(eNU, "Total Energy");
    iEnthalpy_ = addTerm(eNU, "Enthalpy");
    iHV_       = addTerm(eNU, HVstring);

    if (((iVolume_ < INT_MAX) && (iEnthalpy_ < INT_MAX)) &&
        (iHV_ == INT_MAX))
    {
        iHV_          = eNU.size();
        energyTermHV_ = std::unique_ptr<EnergyTerm>(new EnergyTerm(iHV_, true,
                                                                   HVstring,   "kJ/mol nm^3"));
    }
    if (fnFluctConv_.size() > 0)
    {
        fc_.reset(xvgropen(fnFluctConv_.c_str(),
                           "Standard Deviation", "Time (ps)",
                           "", ehelper_.outputEnvironment()));
        ehelper_.printXvgLegend(fc_.get());
    }
    // We need to store all the energy data for later processing.
    ehelper_.setStoreData(true);
}

void FluctPropsModule::analyzeFrame(t_enxframe *fr)
{
    /* We read a valid frame, so we can use it */
    if (fr->nre > 0)
    {
        if (nullptr != fc_)
        {
            fprintf(fc_.get(), "%10g", fr->t);
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
            }
        }
        if (iHV_ != INT_MAX)
        {
            /* Special case for old energy files where HV is not stored. */
            double hv = fr->ener[iEnthalpy_].e * fr->ener[iVolume_].e;
            energyTermHV_->addData(fr->t, fr->step, 1, 0, 0, hv);
            if (nullptr != fc_)
            {
                fprintf(fc_.get(), "  %10g", energyTermHV_->standardDeviation());
            }
        }
        if (nullptr != fc_)
        {
            fprintf(fc_.get(), "\n");
        }
    }
}

void FluctPropsModule::finalizeAnalysis()
{
    bool   bTemp, bTotal, bEnthalpy, bVolume, bHV;
    double Temp, TempStddev;
    double Total, TotalStddev;
    double Enthalpy, EnthalpyStddev;
    double Volume, VolumeStddev;
    double HV, HVStddev;
    double kappa;
    int    nMol = ehelper_.nMol();
    /* Compute it all! */
    kappa = 0;

    if (nMol < 2)
    {
        printf("\nWARNING: #molecules = %d, this may not be what you want.\n",
               nMol);
    }
    printf("\nTemperature dependent fluctuation properties.\n");
    printf("\nHeat capacities obtained from fluctuations do *not* include\n");
    printf("quantum corrections. If you want to get a more accurate estimate\n");
    printf("please use the gmx dos tool.\n\n");
    printf("WARNING: Please verify that your simulations are converged and perform\n"
           "a block-averaging error analysis (not implemented in gmx fluctprops yet)\n\n");

    /* Temperature */
    if (true == (bTemp = ehelper_.energyTerm(F_TEMP, &Temp, &TempStddev)))
    {
        printf("Temperature                              = %10g +/- %10g K\n",
               Temp, TempStddev);
    }

    if (true == (bEnthalpy = ehelper_.energyTerm("Enthalpy", &Enthalpy, &EnthalpyStddev)))
    {
        printf("Enthalpy                                 = %10g +/- %10g kJ/mol\n",
               Enthalpy/nMol, EnthalpyStddev/nMol);
        bTotal = false;
    }
    else
    {
        if (true == (bTotal = ehelper_.energyTerm(F_ETOT, &Total, &TotalStddev)))
        {
            printf("Total energy                         = %10g +/- %10g kJ/mol\n",
                   Total/nMol, TotalStddev/nMol);
        }
        else
        {
            fprintf(stderr, "WARNING: Neither Enthalpy nor Total-Energy in your file.\n");
        }
    }

    if (true == (bVolume = ehelper_.energyTerm("Volume", &Volume, &VolumeStddev)))
    {
        printf("Volume                                   = %10g +/- %10g nm^3/molecule\n",
               Volume/nMol, VolumeStddev/nMol);
    }

    bHV = ehelper_.energyTerm(HVstring, &HV, &HVStddev);

    if (bTemp)
    {
        if (bVolume)
        {
            /* Compute kappa */
            double varv    = square(VolumeStddev);
            kappa   = pow(NANO, 3.0)*(varv/Volume)/(BOLTZMANN*Temp);
            printf("Isothermal Compressibility Kappa         = %10g (1/GPa)\n",
                   kappa/1.0e-9);
            printf("Adiabatic bulk modulus                   = %10g (GPa)\n",
                   1.0e-9/kappa);
        }
        if (bEnthalpy)
        {
            /* Compute cP */
            double varh = square(EnthalpyStddev)/nMol;
            double cp   = KILO*(varh/(BOLTZ*Temp*Temp));
            printf("Heat capacity at constant pressure Cp    = %10g J/mol K\n", cp);
        }
        else if (bTotal)
        {
            /* Total energy.
             * Only compute cv in constant volume runs, which we can test
             * by checking whether the enthalpy was computed.
             */
            double varet = square(TotalStddev);
            double cv    = KILO*((varet/nMol)/(BOLTZ*Temp*Temp));
            printf("Heat capacity at constant volume Cv      = %10g J/mol K\n", cv);
        }
        /* Alpha, dcp */
        if (bHV && kappa > 0)
        {
            double alpha = ((KILO/AVOGADRO)*(HV-Volume*Enthalpy)/
                            (Volume*BOLTZMANN*Temp*Temp));
            double dcp   = ((pow(NANO, 3.0)*Volume*AVOGADRO/nMol)*Temp*
                            square(alpha)/(kappa));
            printf("Coefficient of Thermal Expansion Alpha_P = %10g (1/K)\n",
                   alpha);
            printf("Cp-Cv                                    =  %10g J/mol K\n",
                   dcp);
        }

        please_cite(stdout, "Allen1987a");
    }
    else
    {
        fprintf(stderr, "WARNING: No temperature found in your energy file.\n");
    }
    //ehelper_.printStatistics(stdout);
    viewOutput();
}

void FluctPropsModule::viewOutput()
{
    if (fnFluctConv_.size() > 0)
    {
        do_view(ehelper_.outputEnvironment(), fnFluctConv_.c_str(), "-nxy");
    }
}

const char FluctPropsInfo::name[] = "fluctprops";

const char FluctPropsInfo::shortDescription[] = "Compute heat capacities, and other fluctuation properties";

EnergyAnalysisModulePointer FluctPropsInfo::create()
{
    return EnergyAnalysisModulePointer(new FluctPropsModule());
}

} // namespace

} // namespace gmx
