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
 * Declares gmx::EnergyHandler
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_energyyanalysis
 */
#ifndef GMX_ENERGYANALYSIS_HANDLER_H
#define GMX_ENERGYANALYSIS_HANDLER_H

#include <string>
#include <vector>

#include "gromacs/legacyheaders/typedefs.h"
#include "energyanalysis.h"

namespace gmx
{
/*! \brief
 * Class doing the actual reading of an energy file, and passing the
 * information to different energy analysis tools.
 */
class EnergyHandler
{
    private:
        //! The energy files
        std::vector<std::string>       fileName_;
        //! The energy analysis tools
        std::vector<EnergyAnalysisPtr> eap_;
    public:
        //! Constructor
        EnergyHandler() {};
        //! Destructor
        ~EnergyHandler() {}
        //! Add a handler
        void addAnalysisTool(EnergyAnalysisPtr eap) { eap_.push_back(eap); }
        //! \return the variance of energies
        void addEnergyFile(std::string efile) { fileName_.push_back(efile); }
        /*! \brief
         * Read the files and call the tools to analyze them
         * \return true if everything went smoothly
         */
        bool readFiles();
};

}
#endif
