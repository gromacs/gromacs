/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014,2015,2016,2018,2019, by the GROMACS development team, led by
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

#ifndef GMX_GMXPREPROCESS_TOPIO_H
#define GMX_GMXPREPROCESS_TOPIO_H

#include <memory>
#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

struct gmx_molblock_t;
struct gmx_mtop_t;
class PreprocessingAtomTypes;
struct t_gromppopts;
struct t_inputrec;
struct MoleculeInformation;
struct InteractionsOfType;
struct t_symtab;
struct warninp;
enum struct GmxQmmmMode;
typedef warninp* warninp_t;

double check_mol(const gmx_mtop_t* mtop, warninp_t wi);
/* Check mass and charge */

char** do_top(bool                                  bVerbose,
              const char*                           topfile,
              const char*                           topppfile,
              t_gromppopts*                         opts,
              bool                                  bZero,
              t_symtab*                             symtab,
              gmx::ArrayRef<InteractionsOfType>     plist,
              int*                                  combination_rule,
              double*                               repulsion_power,
              real*                                 fudgeQQ,
              PreprocessingAtomTypes*               atype,
              std::vector<MoleculeInformation>*     molinfo,
              std::unique_ptr<MoleculeInformation>* intermolecular_interactions,
              const t_inputrec*                     ir,
              std::vector<gmx_molblock_t>*          molblock,
              bool*                                 ffParametrizedWithHBondConstraints,
              warninp_t                             wi);

/* This routine expects sys->molt[m].ilist to be of size F_NRE and ordered. */
void generate_qmexcl(gmx_mtop_t* sys, t_inputrec* ir, warninp_t wi, GmxQmmmMode qmmmMode);

#endif
