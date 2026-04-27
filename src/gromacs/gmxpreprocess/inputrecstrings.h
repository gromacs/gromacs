/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \libinternal \file
 * \brief
 * Declares gmx_inputrec_strings for preprocessing-only MDP parameters
 *
 * This structure holds string representations of MDP parameters that are
 * only used during grompp preprocessing. These are not serialized to TPR
 * files because they are converted to numeric indices or other runtime
 * representations during preprocessing.
 *
 * \ingroup module_gmxpreprocess
 * \inlibraryapi
 */
#ifndef GMX_GMXPREPROCESS_INPUTRECSTRINGS_H
#define GMX_GMXPREPROCESS_INPUTRECSTRINGS_H

#include <string>
#include <vector>

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"

/*! \libinternal
 * \brief Storage for preprocessing-only string parameters from MDP files
 *
 * These strings are parsed by grompp to set up groups, tables, and other
 * preprocessing configuration. They are not included in the TPR file because
 * the processed results (e.g., group names and indices) are stored instead.
 *
 * In the long term, storage for these variables should be provided by the
 * relevant modules, who e.g. provide them to the index-group machinery via
 * callback and can also write them back to mdp output when required.
 */
struct gmx_inputrec_strings
{
    char tcgrps[STRLEN], tau_t[STRLEN], ref_t[STRLEN], accelerationGroups[STRLEN],
            acceleration[STRLEN], freeze[STRLEN], frdim[STRLEN], user1[STRLEN], user2[STRLEN],
            vcm[STRLEN], couple_moltype[STRLEN], orirefitgrp[STRLEN], egptable[STRLEN],
            egpexcl[STRLEN], wall_atomtype[STRLEN], wall_density[STRLEN], deform[STRLEN],
            QMMM[STRLEN], imd_grp[STRLEN];
    gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, std::string> fep_lambda;
    char                                                                   lambdaWeights[STRLEN];
    char                                                                   lambdaCounts[STRLEN];
    char                     wlHistogramCounts[STRLEN];
    std::vector<std::string> pullGroupNames;
    std::vector<std::string> rotateGroupNames;
    char anneal[STRLEN], anneal_npoints[STRLEN], anneal_time[STRLEN], anneal_temp[STRLEN];

    std::string compressedXGroups; //!< Groups for compressed trajectory output
    std::string energyGroups;      //!< Groups for energy calculation
};

#endif // GMX_GMXPREPROCESS_INPUTRECSTRINGS_H
