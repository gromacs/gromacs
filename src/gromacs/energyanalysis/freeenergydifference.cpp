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
#include "gromacs/utility/smalloc.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/units.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/utility/arrayref.h"

#include "energyhandler.h"
#include "energyhelper.h"
#include "energyinfo.h"
#include "freeenergydifference.h"
#include "select.h"
#include "simple.h"

namespace gmx
{

class FreeEnergyDifferenceModule : public SimpleEnergyModule
{
    private:
        //! Running average for anlysis
        std::string runAverFile_;

        //! Reference temperature
        double      refTemp_;

        //! Energy helper class for low level stuff. Need two of them.
        EnergyHelper eh_[2];

        //! Data set index
        int fileIndex_;
    public:
        //! Constructor
        FreeEnergyDifferenceModule();

        //! Destructor must be virtual since there are other virtual functions
        virtual ~FreeEnergyDifferenceModule() {};

        //! Return pointer to the EnergyHelper for cleaner code downstream
        EnergyHelper *helper() { return &(eh_[fileIndex_]); }

        /*! \brief
         * Some option
         * \param[in] reftemp Reference Temperature
         */
        void setParameters(double reftemp) { refTemp_ = reftemp; }

        //! Set the fluctuation convergence file name
        void setRunAverFile(std::string runAverFile)
        {
            runAverFile_ = runAverFile;
        }

        //! Get the fluctuation convergence file name
        std::string getRunAverFile() { return runAverFile_; }

        //! Pass the output environment on for printing etc.
        virtual void setOutputEnvironment(output_env_t oenv)
        {
            helper()->setOutputEnvironment(oenv);
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
};

FreeEnergyDifferenceModule::FreeEnergyDifferenceModule()
{
    refTemp_   = 0;
    fileIndex_ = -1;
    ehelper_   = new(EnergyHelper);
}

void FreeEnergyDifferenceModule::initOptions(Options *options)
{
    static const char *const desc[] = {
        "With [TT]-fee[tt] an estimate is calculated for the free-energy",
        "difference with an ideal gas state: [BR]",
        "  [GRK]Delta[grk] A = A(N,V,T) - A[SUB]idealgas[sub](N,V,T) = kT [LN][CHEVRON][EXP]U[SUB]pot[sub]/kT[exp][chevron][ln][BR]",
        "  [GRK]Delta[grk] G = G(N,p,T) - G[SUB]idealgas[sub](N,p,T) = kT [LN][CHEVRON][EXP]U[SUB]pot[sub]/kT[exp][chevron][ln][BR]",
        "where k is Boltzmann's constant, T is set by [TT]-fetemp[tt] and",
        "the average is over the ensemble (or time in a trajectory).",
        "Note that this is in principle",
        "only correct when averaging over the whole (Boltzmann) ensemble",
        "and using the potential energy. This also allows for an entropy",
        "estimate using:[BR]",
        "  [GRK]Delta[grk] S(N,V,T) = S(N,V,T) - S[SUB]idealgas[sub](N,V,T) = ([CHEVRON]U[SUB]pot[sub][chevron] - [GRK]Delta[grk] A)/T[BR]",
        "  [GRK]Delta[grk] S(N,p,T) = S(N,p,T) - S[SUB]idealgas[sub](N,p,T) = ([CHEVRON]U[SUB]pot[sub][chevron] + pV - [GRK]Delta[grk] G)/T",
        "[PAR]",

        "When two energy files are provided and the ([TT]-fediff[tt]) flag is set, a free energy",
        "difference is calculated [BR] dF = -kT [LN][CHEVRON][EXP]-(E[SUB]B[sub]-E[SUB]A[sub])/kT[exp][chevron][SUB]A[sub][ln] ,",
        "where E[SUB]A[sub] and E[SUB]B[sub] are the energies from the first and second energy",
        "files, and the average is over the ensemble A. The running average",
        "of the free energy difference is printed to a file specified by [TT]-ravg[tt].",
        "[BB]Note[bb] that the energies must both be calculated from the same trajectory."
    };
    options->setDescription(desc);


}

bool FreeEnergyDifferenceModule::initAnalysis(std::vector<std::string> eName,
                                              std::vector<std::string> eUnit)
{
    if ((fileIndex_ < 2) && SimpleEnergyModule::initAnalysis(eName, eUnit))
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool FreeEnergyDifferenceModule::addDataSet(std::string name)
{
    if (fileIndex_ < 1)
    {
        fileIndex_++;
        helper()->addDataSet(name);

        return true;
    }
    return false;
}

bool FreeEnergyDifferenceModule::addAnalysisFrame(t_enxframe *fr)
{
    if ((fileIndex_ < 0) || (fileIndex_ > 1))
    {
        return false;
    }
    else
    {
        return SimpleEnergyModule::addAnalysisFrame(fr);
    }
}

bool FreeEnergyDifferenceModule::finalizeAnalysis()
{
    const char * ravgleg[] = {
        "\\8D\\4E = E\\sB\\N-E\\sA\\N",
        "<e\\S-\\8D\\4E/kT\\N>\\s0..t\\N"
    };
    FILE        *fp;
    double       beta;
    bool         bOK;

    /* check */
    if (eh_[0].nEnergyTerm() != eh_[1].nEnergyTerm())
    {
        fprintf(stderr, "Number of energy terms extracted from files differs\n");
        return false;
    }

    // 1/T Boltzmann constant
    beta = 1.0/(BOLTZ*refTemp_);

    /* calculate fe difference dF = -kT ln < exp(-(E_B-E_A)/kT) >_A */
    fp = NULL;
    if (getRunAverFile().size() > 0)
    {
        fp = xvgropen(getRunAverFile().c_str(),
                      "Running average free energy difference",
                      "Time (" unit_time ")", "\\8D\\4E (" unit_energy ")",
                      helper()->getOutputEnvironment());
        xvgr_legend(fp, asize(ravgleg), ravgleg, helper()->getOutputEnvironment());
    }
    fprintf(stdout, "\n%-24s %10s\n",
            "Energy", "dF = -kT ln < exp(-(EB-EA)/kT) >A");

    bOK = true;
    for (EnergyTermIterator et0 = eh_[0].etBegin();
         bOK && (et0 < eh_[0].etEnd()); ++et0)
    {
        EnergyTermIterator et1 = eh_[1].etSearch(et0->getEterm());
        if (eh_[1].etEnd() == et1)
        {
            fprintf(stderr, "Can not find %s in second file.\n",
                    et0->getEterm().c_str());
            bOK = false;
        }
        if (bOK && (et0->nEnergy() != et1->nEnergy()))
        {
            fprintf(stderr, "Not same length of energy files\n");
            bOK = false;
        }
        double sum = 0;
        for (int i = 0; (i < et0->nEnergy()); i++)
        {
            EnergyFrameIterator ef0  = et0->searchEF(i);
            EnergyFrameIterator ef1  = et1->searchEF(i);
            double              dE   = ef0->getE() - ef1->getE();
            sum        += exp(-dE*beta);
            if (NULL != fp)
            {
                fprintf(fp, "%10g %10g %10g\n",
                        ef0->getT(), dE, -BOLTZ*refTemp_*log(sum/(i+1)) );
            }
        }
        double aver = -BOLTZ*refTemp_*log(sum/et0->nEnergy());
        fprintf(stdout, "%-24s %10g\n", et0->getEterm().c_str(), aver);
    }
    if (fp)
    {
        gmx_ffclose(fp);
    }

    return bOK;
}

const char FreeEnergyDifferenceInfo::name[] = "fediff";

const char FreeEnergyDifferenceInfo::shortDescription[] = "Compute free energy estimate from two energy files";

CommandLineOptionsModuleInterface *FreeEnergyDifferenceInfo::create()
{
    EnergyAnalysisModulePointer eamp(new FreeEnergyDifferenceModule);
    return (CommandLineOptionsModuleInterface *)(new EnergyInfo(eamp));
}
}
