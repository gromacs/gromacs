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
/*! \libinternal \file
 * \brief
 * Declares gmx::energyanalysis::SimpleEnergy and gmx::energyanalysis::SimpleInfo
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_energyanalysis
 */
#ifndef GMX_ENERGYANALYSIS_SIMPLE_H
#define GMX_ENERGYANALYSIS_SIMPLE_H

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "analysismodule.h"
#include "energyhelper.h"

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
    protected:
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

    public:
        //! Constructor
        SimpleEnergyModule();

        //! Destructor must be virtual since there are other virtual functions
        virtual ~SimpleEnergyModule() {};

        //! Get the output file name
        std::string getOutputFile() { return fnEnergy_; }

        //! Set the output file name
        void setOutputFile(std::string fnEnergy) { fnEnergy_ = fnEnergy; }

        //! Get the fluctuation convergence file name
        std::string getFluctConvFile() { return fnFluctConv_; }

        //! Set the summing
        void setSumming(bool bSum) { bSum_ = bSum; }

        //! Get the summing status
        bool getSumming() { return bSum_; }

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
        virtual void initOptions(Options *options);

        //! Initiate the command line options
        virtual void optionsFinished(Options */*options*/)
        {
        }

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

        //! Sum all the energy terms and delete the original data sets
        void sumEnergies();

        //! View the output file(s)
        void viewOutput();
};

class SimpleInfo
{
    public:
        static const char name[];
        static const char shortDescription[];
        static gmx::CommandLineOptionsModuleInterface *create();
};

}

}

#endif
