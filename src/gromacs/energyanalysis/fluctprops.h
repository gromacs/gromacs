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
 * Declares gmx::FluctProps
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_energyanalysis
 */
#ifndef GMX_ENERGYANALYSIS_FLUCTPROPS_H
#define GMX_ENERGYANALYSIS_FLUCTPROPS_H

#include <stdio.h>
#include "energyanalysis.h"
#include "helper.h"

namespace gmx
{
class FluctProps : public EnergyAnalysis
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

        //! Energy helper class for low level stuff
        EnergyHelper          eh_;

        //! Number of blocks for error estimate
        int                   nBlock_;

        //! Number of molecules
        int                   nMol_;

        //! Should we correct for the drift in the data?
        bool                  bDriftCorr_;
    public:
        //! Constructor
        FluctProps();

        //! Destructor must be virtual since there are other virtual functions
        virtual ~FluctProps() {};

        //! Return pointer to the EnergyHelper for cleaner code downstream
        EnergyHelper *helper() { return &eh_; }

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
            helper()->setOutputEnvironment(oenv);
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

};

#endif
