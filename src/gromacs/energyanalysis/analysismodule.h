/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
 * Declares gmx::energyanalysis::EnergyAnalysisModule
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_energyanalysis
 */
#ifndef GMX_ENERGYANALYSIS_ANALYSISMODULE_H
#define GMX_ENERGYANALYSIS_ANALYSISMODULE_H

#include <string>
#include <vector>

#include "gromacs/options.h"
#include "gromacs/commandline/cmdlineoptionsmodule.h"

struct gmx_output_env_t;
struct t_enxframe;

namespace gmx
{

/*! \libinternal \brief
 * Convenience structure for keeping energy name and unit together
 */
struct EnergyNameUnit
{
    //! Name of this energy term
    std::string energyName;
    //! Unit of this energy term
    std::string energyUnit;
};

/*! \libinternal
 * \brief
 * Pure class overloaded by the separate energy modules.
 */
class EnergyAnalysisModule
{
    public:
        //! Obligatory virtual destructor
        virtual ~EnergyAnalysisModule() {};

        //! Initiate the command line options
        virtual void initOptions(IOptionsContainer                 *options,
                                 ICommandLineOptionsModuleSettings *settings) = 0;

        /*! \brief
         * Does the initiation of the analysis of the file
         * \param[in] eNU   Names and units of the energy terms.
         * \param[in] oenv  GROMACS output environment
         */
        virtual void initAnalysis(ArrayRef<const EnergyNameUnit> eNU,
                                  const gmx_output_env_t        *oenv) = 0;

        /*! \brief
         * Analyse one frame and stores the results in memory
         * \param[in] fr The energy data frame
         * \param[in] oenv  GROMACS output environment. This is needed in some cases
         *                  where the output to be written depends on the content of
         *                  the energy file.
         */
        virtual void analyzeFrame(t_enxframe *fr, const gmx_output_env_t *oenv) = 0;

        /*! \brief
         * Finalize reading and write output files.
         * \param[in] oenv  GROMACS output environment.
         */
        virtual void finalizeAnalysis(const gmx_output_env_t *oenv) = 0;

        /*! \brief
         * View the output file(s)
         * \param[in] oenv  GROMACS output environment.
         */
        virtual void viewOutput(const gmx_output_env_t *oenv) = 0;
};

//! Pointer to the EnergyAnalysisModule classes.
typedef std::unique_ptr<EnergyAnalysisModule> EnergyAnalysisModulePointer;

} // namespace gmx

#endif
