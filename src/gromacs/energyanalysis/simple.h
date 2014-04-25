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
/*! \internal \file
 * \brief
 * Declares gmx::SimpleEnergy
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_energyanalysis
 */
#ifndef GMX_ENERGYANALYSIS_SIMPLE_H
#define GMX_ENERGYANALYSIS_SIMPLE_H

#include <stdio.h>
#include "helper.h"
#include "energyanalysis.h"

namespace gmx
{
class SimpleEnergy : public EnergyAnalysis
{
    private:
        //! File pointer for storing energies
        FILE        *fp_;

        //! Fluctuation convergence output (typically a xvg file)
        std::string  fluctConvFile_;

        //! File pointer for storing flucutuations
        FILE        *fc_;

        //! Boolean instructing us whether to sum the energy terms
        bool         bSum_;

        //! Energy helper class for low level stuff
        EnergyHelper eh_;

    public:
        //! Constructor
        SimpleEnergy();

        //! Destructor must be virtual since there are other virtual functions
        virtual ~SimpleEnergy() {};

        //! Return pointer to the EnergyHelper for cleaner code downstream
        EnergyHelper *helper() { return &eh_; }

        //! Set the fluctuation convergence file name
        void setFluctConvFile(std::string fluctConvFile)
        {
            fluctConvFile_ = fluctConvFile;
        }

        //! Get the fluctuation convergence file name
        std::string getFluctConvFile() { return fluctConvFile_; }

        //! Set the summing
        void setSumming(bool bSum) { bSum_ = bSum; }

        //! Get the summing status
        bool getSumming() { return bSum_; }

        /*! \brief
         * Does the initiation of the analysis of the file
         * \param[in] nre Number of energy terms in the file
         * \param[in] enm Names of the energy terms etc.
         * \return true if OK, false otherwise
         */
        virtual bool initAnalysis(int nre, gmx_enxnm_t enm[]);

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

};

#endif
