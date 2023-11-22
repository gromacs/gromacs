/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 *
 * \brief Functionality for testing whether calls to mdrun produce the
 * same energy and force quantities when they should do so.
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_SIMULATIONDATABASE_H
#define GMX_TESTUTILS_SIMULATIONDATABASE_H

#include <map>
#include <string>

namespace gmx
{

namespace test
{

/*! \brief Return whether the number of ranks is supported by
 * the simulation \c simulationName in the database.
 *
 * This method lets test runners understand when end-to-end tests
 * should be expected to work. */
bool isNumberOfPpRanksSupported(const std::string& simulationName, int possibleNumberOfPpRanks);

/*! \brief Return a string describing the numbers of ranks supported
 * for the simulation \c simulationName in the database. */
std::string reportNumbersOfPpRanksSupported(const std::string& simulationName);

//! Helper typedef
using MdpFieldValues = std::map<std::string, std::string>;

//! Strong type denoting mdp parameters needed to test specific algorithms
enum class MdpParameterDatabase
{
    Default,
    Pull,
    Awh,
    ExpandedEnsemble,
    Count
};

/*! \brief Set up values for an .mdp file that permits a highly
 * reproducible simulation.
 *
 * An internal database of several kinds of simulation useful for such
 * comparisons is available, whose \c simulationName keys are
 *     - argon12
 *     - argon5832
 *     - tip3p5
 *     - spc2
 *     - spc216
 *     - alanine_vsite_vacuo
 *     - alanine_vsite_solvated
 *     - glycine_vacuo
 *     - glycine_no_constraints_vacuo
 *     - nonanol_vacuo
 *
 * Some of these systems are pretty minimal, because having
 * few atoms means few interactions, highly reproducible
 * forces, and allows tests to focus on the correctness of the
 * implementation of high-level mdrun features. The boxes are
 * of a reasonable size so that domain decomposition is
 * possible. The pressure-coupling parameters are isotropic,
 * and set up so that there will not be dramatic collapse of
 * volume over the handful of MD steps that will be run. A
 * single temperature-coupling group is used.
 *
 * \param[in]  simulationName           The name of the simulation, which indexes the database
 * \param[in]  integrator               The integrator to use
 * \param[in]  tcoupl                   The temperature-coupling algorithm to use
 * \param[in]  pcoupl                   The pressure-coupling algorithm to use
 * \param[in]  additionalMdpParameters  Additional mdp parameters to be added
 * \return                              Mdp file values
 *
 * \throws  std::bad_alloc     if out of memory
 * \throws  std::out_of_range  if \c simulationName is not in the database */
MdpFieldValues prepareMdpFieldValues(const std::string& simulationName,
                                     const std::string& integrator,
                                     const std::string& tcoupl,
                                     const std::string& pcoupl,
                                     MdpParameterDatabase additionalMdpParameters = MdpParameterDatabase::Default);

//! \copydoc prepareMdpFieldValues()
MdpFieldValues prepareMdpFieldValues(const char* simulationName,
                                     const char* integrator,
                                     const char* tcoupl,
                                     const char* pcoupl,
                                     MdpParameterDatabase additionalMdpParameters = MdpParameterDatabase::Default);

/*! \brief Make a string containing an .mdp file from the \c mdpFieldValues.
 *
 * \throws  std::bad_alloc     if out of memory */
std::string prepareMdpFileContents(const MdpFieldValues& mdpFieldValues);

} // namespace test
} // namespace gmx

#endif
