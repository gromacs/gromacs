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
#include "gromacs/commandline/cmdlineoptionsmodule.h"

#include "analysismodule.h"
#include "energyhandler.h"

namespace gmx
{

namespace energyanalysis
{

/*! \libinternal
 * \brief
 * Low level helper class.
 */
class EnergyInfo : public gmx::CommandLineOptionsModuleInterface
{
    public:
        /*! \brief Constructor
         *
         * \param[in] eamp Pointer to EnergyAnalysisModule
         */
        EnergyInfo(EnergyAnalysisModulePointer eamp) : ehandler_(eamp)
        {
            eamp_     = eamp;
        }

        //! Initiate the command line module settings
        virtual void init(CommandLineModuleSettings * /*settings*/) {};

        //! Initiate the command line options
        virtual void initOptions(Options *options)
        {
            ehandler_.initOptions(options);
        }

        //! Called when parsing options has been finished
        virtual void optionsFinished(Options *options)
        {
            ehandler_.optionsFinished(options);
        }

        //! Does the actual work.
        virtual int run() { return ehandler_.run(); }

    private:
        //! The module doing the real work
        EnergyAnalysisModulePointer eamp_;
        //! The object handling file I/O etc.
        EnergyHandler               ehandler_;
};

}

}
