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

#ifndef GMX_GMXPREPROCESS_TOPIO_H
#define GMX_GMXPREPROCESS_TOPIO_H

#include <filesystem>
#include <memory>
#include <optional>
#include <vector>

#include "gromacs/utility/real.h"

struct gmx_molblock_t;
struct gmx_mtop_t;
class PreprocessingAtomTypes;
struct t_gromppopts;
struct t_inputrec;
struct MoleculeInformation;
struct InteractionsOfType;
struct t_symtab;
class WarningHandler;
enum class CombinationRule : int;

namespace gmx
{
template<typename>
class ArrayRef;
class MDLogger;
} // namespace gmx

double check_mol(const gmx_mtop_t* mtop, WarningHandler* wi);
/* Check mass and charge */

/*! \brief Check that RB/Fourier style dihedral sums to zero in free-energy computation
 *
 * Searches for, and emits a warning when dihedral of function type 3 (Ryckaert-Bellemans or
 * Fourier) have a sum of coefficient different from zero, when performing a free-energy
 * computation, in stateA and stateB, as this is likely a parameter mistake that can
 * produce erroneous dHdl values.
 *
 */
void checkRBDihedralSum(const gmx_mtop_t& mtop, const t_inputrec& ir, WarningHandler* wi);


char** do_top(bool                                        bVerbose,
              const char*                                 topfile,
              const std::optional<std::filesystem::path>& topppfile,
              t_gromppopts*                               opts,
              bool                                        bZero,
              t_symtab*                                   symtab,
              gmx::ArrayRef<InteractionsOfType>           plist,
              CombinationRule*                            combination_rule,
              double*                                     repulsion_power,
              real*                                       fudgeQQ,
              PreprocessingAtomTypes*                     atype,
              std::vector<MoleculeInformation>*           molinfo,
              std::unique_ptr<MoleculeInformation>*       intermolecular_interactions,
              const t_inputrec*                           ir,
              std::vector<gmx_molblock_t>*                molblock,
              bool*                                       ffParametrizedWithHBondConstraints,
              WarningHandler*                             wi,
              const gmx::MDLogger&                        logger);

/* This routine expects sys->molt[m].ilist to be of size F_NRE and ordered. */
void generate_qmexcl(gmx_mtop_t* sys, t_inputrec* ir, const gmx::MDLogger& logger);

#endif
