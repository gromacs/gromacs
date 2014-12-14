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
 * Declares gmx::energyanalysis::EnergyInfo
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_energyanalysis
 */
#ifndef GMX_ENERGYANALYSIS_ENERGYINFO_H
#define GMX_ENERGYANALYSIS_ENERGYINFO_H
#include "gromacs/commandline/cmdlineoptionsmodule.h"

#include "analysismodule.h"
#include "energyhandler.h"

namespace gmx
{

/*! \libinternal
 * \brief
 * Low level helper class.
 */
class EnergyInfo : public ICommandLineOptionsModule
{
    public:
        /*! \brief Constructor
         *
         * \param[in] eamp Pointer to EnergyAnalysisModule
         */
        EnergyInfo(energyanalysis::EnergyAnalysisModulePointer eamp) : ehandler_(eamp)
        {
            eamp_     = eamp;
        }

        virtual void init(CommandLineModuleSettings * /*settings*/) {};

        virtual void initOptions(IOptionsContainer                 *options,
                                 ICommandLineOptionsModuleSettings *settings)
        {
            ehandler_.initOptions(options, settings);
        }

        virtual void optionsFinished() {}

        //! Does the actual work.
        virtual int run() { return ehandler_.run(); }

    private:
        //! The module doing the real work
        energyanalysis::EnergyAnalysisModulePointer eamp_;
        //! The object handling file I/O etc.
        energyanalysis::EnergyHandler               ehandler_;
};

//! Type to make sure these pointers can not be copied
//typedef gmx_unique_ptr<EnergyInfo>::type EnergyInfoPointer;
typedef ICommandLineOptionsModule *EnergyInfoPointer;
//typedef boost::shared_ptr<EnergyInfo> EnergyInfoPointer;

}

#endif
