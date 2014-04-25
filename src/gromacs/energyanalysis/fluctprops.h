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
#ifndef GMX_FLUCTPROPS_H
#define GMX_FLUCTPROPS_H

#include <stdio.h>
#include "energyanalysis.h"

namespace gmx
{
class FluctProps : public EnergyAnalysis
{
    private:
        //! Fluctuation convergence output (typically a xvg file)
        std::string  fluctConvFile_;

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

    public:
        //! Constructor
        FluctProps();

        //! Destructor
        ~FluctProps() {};

        /*! \brief
         * Add a term to the data structure and return the energy file index
         * \param[in] nre Number of energy terms in the input file
         * \param[in] enm The names of energy terms
         * \param[in] term The name of the term to search for
         * \return The energy file index for this term or INT_MAX if not succesfull
         */
        unsigned int searchTerm(int nre, gmx_enxnm_t enm[], const char *term);

        //! Set the fluctuation convergence file name
        void setFluctConvFile(std::string fluctConvFile)
        {
            fluctConvFile_ = fluctConvFile;
        }

        //! Get the fluctuation convergence file name
        std::string getFluctConvFile() { return fluctConvFile_; }

        //! Does the initiation of the reading of the file
        virtual bool initAnalysis(int nre, gmx_enxnm_t enm[]);

        //! Analyse one frame and stores the results in memory
        virtual bool addAnalysisFrame(t_enxframe *fr);

        //! Finalize reading
        virtual bool finalizeAnalysis();
};

};

#endif
