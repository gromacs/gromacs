/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
#ifndef GMX_MDLIB_FORCEREC_H
#define GMX_MDLIB_FORCEREC_H

#include <cstdio>

#include <string>
#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

struct gmx_hw_info_t;
struct t_commrec;
struct t_forcerec;
struct t_filenm;
struct t_inputrec;
struct gmx_localtop_t;
struct gmx_mtop_t;
struct gmx_wallcycle;
struct interaction_const_t;
union t_iparams;
enum class LongRangeVdW : int;

namespace gmx
{
class MDLogger;
class PhysicalNodeCommunicator;
class SimulationWorkload;
} // namespace gmx

/*! \brief Create nonbonded parameter lists
 *
 * \param[in] numAtomTypes           The number of atom types
 * \param[in] iparams                The LJ parameters
 * \param[in] useBuckinghamPotential Use Buckingham potential
 */
std::vector<real> makeNonBondedParameterLists(int                            numAtomTypes,
                                              gmx::ArrayRef<const t_iparams> iparams,
                                              bool useBuckinghamPotential);

/*! \brief Calculate c6 parameters for grid correction
 *
 * \param[in] numAtomTypes           The number of atom types
 * \param[in] iparams                The LJ parameters
 * \param[in] ljpme_combination_rule How long range LJ is treated
 */
std::vector<real> makeLJPmeC6GridCorrectionParameters(int                            numAtomTypes,
                                                      gmx::ArrayRef<const t_iparams> iparams,
                                                      LongRangeVdW ljpme_combination_rule);

/*! \brief Set the number of charge groups and atoms.
 *
 * The force calculation needs information on which atoms it
 * should do work.
 * \param[inout] fr                  The forcerec
 * \param[in]    natoms_force        Number of atoms to compute force on
 * \param[in]    natoms_force_constr Number of atoms involved in constraints
 * \param[in]    natoms_f_novirsum   Number of atoms for which
 *                                   force is to be compute but no virial
 */
void forcerec_set_ranges(t_forcerec* fr, int natoms_force, int natoms_force_constr, int natoms_f_novirsum);

/*! \brief Initiate table constants
 *
 * Initializes the tables in the interaction constant data structure.
 * \param[in] fp                     File for debugging output
 * \param[in] ic                     Structure holding the table constant
 * \param[in] rlist                  Length of the neighbour list
 * \param[in] tableExtensionLength   Length by which to extend the tables. Taken from the input record.
 */
void init_interaction_const_tables(FILE* fp, interaction_const_t* ic, real rlist, real tableExtensionLength);

/*! \brief Initialize forcerec structure.
 *
 * \param[in]  fplog              File for printing
 * \param[in]  mdlog              File for printing
 * \param[out] forcerec                 The forcerec
 * \param[in]  simulationWork           Simulation workload flags
 * \param[in]  inputrec                 Inputrec structure
 * \param[in]  mtop               Molecular topology
 * \param[in]  commrec                 Communication structures
 * \param[in]  box                Simulation box
 * \param[in]  tabfn              Table potential file for non-bonded interactions
 * \param[in]  tabpfn             Table potential file for pair interactions
 * \param[in]  tabbfnm            Table potential files for bonded interactions
 * \param[in]  print_force        Print forces for atoms with force >= print_force
 */
void init_forcerec(FILE*                            fplog,
                   const gmx::MDLogger&             mdlog,
                   const gmx::SimulationWorkload&   simulationWork,
                   t_forcerec*                      forcerec,
                   const t_inputrec&                inputrec,
                   const gmx_mtop_t&                mtop,
                   const t_commrec*                 commrec,
                   matrix                           box,
                   const char*                      tabfn,
                   const char*                      tabpfn,
                   gmx::ArrayRef<const std::string> tabbfnm,
                   real                             print_force);

#endif
