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
#include <stdio.h>
#include <string.h>

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
#include "gromacs/utility/smalloc.h"

#include "analysismodule.h"
#include "energyhandler.h"
#include "energyinfo.h"
#include "energyhelper.h"
#include "fluctprops.h"

namespace gmx
{

namespace energyanalysis
{

class FluctPropsModule : public EnergyAnalysisModule
{
    private:
        //! Fluctuation convergence output (typically a xvg file)
        std::string  fnFluctConv_;

        //! File pointer for storing flucutuations
        FILE                 *fc_;

        //! Index to Enthalpy (not found if < 0)
        unsigned int          iEnthalpy_;

        //! Index to Volume (not found if < 0)
        unsigned int          iVolume_;

        //! Index to Temperature (not found if < 0)
        unsigned int          iTemp_;

        //! Index to Total Energy (not found if < 0)
        unsigned int          iTotal_;

        //! Index to Enthalpy*Volume (not found if < 0)
        unsigned int          iHV_;

        //! Number of blocks for error estimate
        int                   nBlock_;

        //! Should we correct for the drift in the data?
        bool                  bDriftCorr_;

        //! Energy helper class for low level stuff
        EnergyHelper  *ehelper_;
    public:
        //! Constructor
        FluctPropsModule();

        //! Destructor must be virtual since there are other virtual functions
        virtual ~FluctPropsModule() {};

        /*! \brief
         * Add a term to the data structure and return the energy file index
         * \param[in] eName The names of energy terms
         * \param[in] eUnit The units of energy terms
         * \param[in] term  The name of the term to search for
         * \return The energy file index for this term or INT_MAX if not succesfull
         */
        unsigned int addTerm(std::vector<std::string> eName,
                             std::vector<std::string> eUnit,
                             std::string              term);

        //! Set the fluctuation convergence file name
        void setFluctConvFile(std::string fluctConvFile)
        {
            fnFluctConv_ = fluctConvFile;
        }

        //! Get the fluctuation convergence file name
        std::string getFluctConvFile() { return fnFluctConv_; }

        //! Pass the output environment on for printing etc.
        virtual void setOutputEnvironment(output_env_t oenv)
        {
            ehelper_->setOutputEnvironment(oenv);
        }

        //! Initiate the command line options
        void initOptions(Options *options);

        /*! \brief
         * Does the initiation of the analysis of the file
         * \param[in] eName Names of the energy terms etc.
         * \param[in] eUnit Units of the energy terms etc.
         * \return true if OK, false otherwise
         */
        virtual bool initAnalysis(std::vector<std::string> eName,
                                  std::vector<std::string> eUnit);

        /*! \brief
         * Start a new data set (either from a file or something else).
         * \param name String describing the new data set
         * \return true if OK, false otherwise
         */
        virtual bool addDataSet(std::string name);

        /*! \brief
         * Analyse one frame and stores the results in memory
         * \param[in] fr The energy data frame
         * \return true if OK, false otherwise
         */
        virtual bool addAnalysisFrame(t_enxframe *fr);

        //! Finalize reading
        virtual bool finalizeAnalysis();

        //! View the output file(s)
        void viewOutput();
};

//! String describing the product of Enthalpy and Volume
static const char *HVstring = "HV";

FluctPropsModule::FluctPropsModule()
{
    fc_         = NULL;
    iEnthalpy_  = INT_MAX;
    iVolume_    = INT_MAX;
    iTemp_      = INT_MAX;
    iTotal_     = INT_MAX;
    iHV_        = INT_MAX;
    nBlock_     = 5;
    bDriftCorr_ = false;
    ehelper_    = new(EnergyHelper);
};

unsigned int FluctPropsModule::addTerm(std::vector<std::string> eName,
                                       std::vector<std::string> eUnit,
                                       std::string              term)
{
    for (unsigned int i = 0; (i < eName.size()); i++)
    {
        if (eName[i].compare(term) == 0)
        {
            EnergyTerm et(i, ehelper_->getStoreData(), term, eUnit[i]);
            ehelper_->addEnergyTerm(et);
            return i;
        }
    }
    return INT_MAX;
}

bool FluctPropsModule::addDataSet(std::string name)
{
    ehelper_->addDataSet(name);
    return true;
}

void FluctPropsModule::initOptions(Options *options)
{
    static const char *const desc[] = {
        "[THISMODUE] will compute the following properties:[BR]",
        "Property                        Energy terms needed[BR]",
        "---------------------------------------------------[BR]",
        "Heat capacity C[SUB]p[sub] (NPT sims):    Enthalpy, Temp     [BR]",
        "Heat capacity C[SUB]v[sub] (NVT sims):    Etot, Temp         [BR]",
        "Thermal expansion coeff. (NPT): Enthalpy, Vol, Temp[BR]",
        "Isothermal compressibility:     Vol, Temp          [BR]",
        "Adiabatic bulk modulus:         Vol, Temp          [BR]",
        "---------------------------------------------------[BR]",
        "You always need to set the number of molecules [TT]-nmol[tt].",
        "The C[SUB]p[sub]/C[SUB]v[sub] computations do [BB]not[bb] include any corrections",
        "for quantum effects. Use the [gmx-dos] program if you need that (and you do).[PAR]"
    };
    options->setDescription(desc);

    options->addOption(BooleanOption("driftcorr")
                           .store(&bDriftCorr_)
                           .description("The drift in the observables will be subtracted before computing the fluctuation properties"));
    options->addOption(IntegerOption("nblocks")
                           .store(&nBlock_)
                           .description("Number of blocks for error estimate"));

    // Add option for optional output files
    options->addOption(FileNameOption("convergence")
                           .filetype(eftPlot)
                           .outputFile()
                           .store(&fnFluctConv_)
                           .defaultBasename("fluct_conv")
                           .description("Convergence of fluctuation properties"));
    ehelper_->initOptions(options);
}

bool FluctPropsModule::initAnalysis(std::vector<std::string> eName,
                                    std::vector<std::string> eUnit)
{
    // Select which energies to use
    for (EnergyTermIterator eti = ehelper_->etBegin(); (eti < ehelper_->etEnd()); ++eti)
    {
        eti->reset();
    }
    iTemp_     = addTerm(eName, eUnit, "Temperature");
    iVolume_   = addTerm(eName, eUnit, "Volume");
    iTotal_    = addTerm(eName, eUnit, "Total Energy");
    iEnthalpy_ = addTerm(eName, eUnit, "Enthalpy");
    iHV_       = addTerm(eName, eUnit, HVstring);

    if (((iVolume_ < INT_MAX) && (iEnthalpy_ < INT_MAX)) &&
        (iHV_ == INT_MAX))
    {
        iHV_ = eName.size();
        EnergyTerm et(iHV_, true, HVstring, "kJ/mol");
        ehelper_->addEnergyTerm(et);
    }
    if (getFluctConvFile().size() > 0)
    {
        fc_ = xvgropen(getFluctConvFile().c_str(),
                       "Standard Deviation", "Time (ps)",
                       "", ehelper_->getOutputEnvironment());
        ehelper_->printXvgLegend(fc_);
    }
    // We need to store all the energy data for later processing.
    ehelper_->setStoreData(true);

    return true;
}

bool FluctPropsModule::addAnalysisFrame(t_enxframe *fr)
{
    /* We read a valid frame, so we can use it */
    if (fr->nre > 0)
    {
        if (NULL != fc_)
        {
            fprintf(fc_, "%10g", fr->t);
        }
        for (EnergyTermIterator eti = ehelper_->etBegin(); (eti < ehelper_->etEnd()); ++eti)
        {
            unsigned int findex = eti->getIndex();
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
            }
            else if ((iHV_ == findex) && (iEnthalpy_ < INT_MAX) && (iVolume_ < INT_MAX))
            {
                /* Special case for old energy files where HV is not stored. */
                double hv = fr->ener[iEnthalpy_].e * fr->ener[iVolume_].e;
                eti->addData(fr->t, fr->step, 1, 0, 0, hv);
                if (NULL != fc_)
                {
                    fprintf(fc_, "  %10g", eti->standardDeviation());
                }
            }
        }
        if (NULL != fc_)
        {
            fprintf(fc_, "\n");
        }
    }
    return true;
}

bool FluctPropsModule::finalizeAnalysis()
{
    bool   bTemp, bTotal, bEnthalpy, bVolume, bHV;
    double Temp, TempStddev;
    double Total, TotalStddev;
    double Enthalpy, EnthalpyStddev;
    double Volume, VolumeStddev;
    double HV, HVStddev;
    double kappa, cp, dcp, cv;
    int    nMol = ehelper_->getNmol();
    /* Compute it all! */
    kappa = cp = dcp = cv = NOTSET;

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
    if (true == (bTemp = ehelper_->getEnergyTerm(F_TEMP, &Temp, &TempStddev)))
    {
        printf("Temperature                              = %10g +/- %10g K\n",
               Temp, TempStddev);
    }

    if (true == (bEnthalpy = ehelper_->getEnergyTerm("Enthalpy", &Enthalpy, &EnthalpyStddev)))
    {
        printf("Enthalpy                                 = %10g +/- %10g kJ/mol\n",
               Enthalpy/nMol, EnthalpyStddev/nMol);
        bTotal = false;
    }
    else
    {
        if (true == (bTotal = ehelper_->getEnergyTerm(F_ETOT, &Total, &TotalStddev)))
        {
            printf("Total energy                         = %10g +/- %10g kJ/mol\n",
                   Total/nMol, TotalStddev/nMol);
        }
        else
        {
            fprintf(stderr, "WARNING: Neither Enthalpy nor Total-Energy in your file.\n");
        }
    }

    if (true == (bVolume = ehelper_->getEnergyTerm("Volume", &Volume, &VolumeStddev)))
    {
        printf("Volume                                   = %10g +/- %10g nm^3/molecule\n",
               Volume/nMol, VolumeStddev/nMol);
    }

    bHV = ehelper_->getEnergyTerm(HVstring, &HV, &HVStddev);

    if (bTemp)
    {
        if (bVolume)
        {
            /* Compute kappa */
            double varv    = dsqr(VolumeStddev);
            kappa   = pow(NANO, 3.0)*(varv/Volume)/(BOLTZMANN*Temp);
            printf("Isothermal Compressibility Kappa         = %10g (1/GPa)\n",
                   kappa/1.0e-9);
            printf("Adiabatic bulk modulus                   = %10g (GPa)\n",
                   1.0e-9/kappa);
        }
        if (bEnthalpy)
        {
            /* Compute cP */
            double varh = dsqr(EnthalpyStddev)/nMol;
            cp          = KILO*(varh/(BOLTZ*Temp*Temp));
            printf("Heat capacity at constant pressure Cp    = %10g J/mol K\n", cp);
        }
        else if (bTotal)
        {
            /* Total energy.
             * Only compute cv in constant volume runs, which we can test
             * by checking whether the enthalpy was computed.
             */
            double varet = dsqr(TotalStddev);
            cv    = KILO*((varet/nMol)/(BOLTZ*Temp*Temp));
            printf("Heat capacity at constant volume Cv      = %10g J/mol K\n", cv);
        }
        /* Alpha, dcp */
        if (bHV)
        {
            double alpha   = ((KILO/AVOGADRO)*(HV-Volume*Enthalpy)/
                              (Volume*BOLTZMANN*Temp*Temp));
            printf("Coefficient of Thermal Expansion Alpha_P = %10g (1/K)\n",
                   alpha);
            dcp     = ((pow(NANO, 3.0)*Volume*AVOGADRO/nMol)*Temp*
                       dsqr(alpha)/(kappa));
            printf("Cp-Cv                                    =  %10g J/mol K\n",
                   dcp);
        }

        please_cite(stdout, "Allen1987a");
    }
    else
    {
        fprintf(stderr, "WARNING: No temperature found in your energy file.\n");
    }
    if (NULL != fc_)
    {
        xvgrclose(fc_);
    }
    //ehelper_->printStatistics(stdout);
    viewOutput();
    return true;
}

void FluctPropsModule::viewOutput()
{
    if (getFluctConvFile().size() > 0)
    {
        do_view(ehelper_->getOutputEnvironment(), getFluctConvFile().c_str(), "-nxy");
    }
}

const char FluctPropsInfo::name[] = "fluctprops";

const char FluctPropsInfo::shortDescription[] = "Compute heat capacities, and other fluctuation properties";

gmx::CommandLineOptionsModuleInterface *FluctPropsInfo::create()
{
    EnergyAnalysisModulePointer eamp(new FluctPropsModule);
    return (gmx::CommandLineOptionsModuleInterface *)(new EnergyInfo(eamp));
}

}

}
