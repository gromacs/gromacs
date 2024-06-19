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
/*! \libinternal \file
 * \brief
 * Declares functions for handling orientation restraints.
 *
 * \inlibraryapi
 * \ingroup module_listed_forces
 */
#ifndef GMX_LISTED_FORCES_ORIRES_H
#define GMX_LISTED_FORCES_ORIRES_H

#include <cstdio>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/real.h"

struct gmx_mtop_t;
struct gmx_multisim_t;
class history_t;
struct t_inputrec;
struct t_pbc;
struct t_commrec;
struct t_oriresdata;
struct t_disresdata;
struct t_fcdata;
class t_state;
struct t_mdatoms;
union t_iparams;

namespace gmx
{
template<typename>
class ArrayRef;
} // namespace gmx

/*! \brief Extends \p globalState with orientation restraint history
 * when there are restraints and time averaging is used.
 */
void extendStateWithOriresHistory(const gmx_mtop_t& mtop, const t_inputrec& ir, t_state* globalState);

/*! \brief
 * Calculates the time averaged D matrices, the S matrix for each experiment.
 *
 * Returns the weighted RMS deviation of the orientation restraints.
 */
real calc_orires_dev(const gmx_multisim_t*          ms,
                     int                            nfa,
                     const t_iatom                  fa[],
                     const t_iparams                ip[],
                     gmx::ArrayRef<const gmx::RVec> xWholeMolecules,
                     const rvec                     x[],
                     const t_pbc*                   pbc,
                     t_oriresdata*                  oriresdata);

/*! \brief
 * Diagonalizes the order tensor(s) of the orienation restraints.
 *
 * For each experiment eig containts first 3 eigenvalues and then
 * the 3 eigenvectors. The eigenvalues are ordered on magnitude.
 */
void diagonalize_orires_tensors(t_oriresdata* od);

//! Prints order parameter, eigenvalues and eigenvectors to the log file.
void print_orires_log(FILE* log, t_oriresdata* od);

//! Calculates the orientation restraint forces.
real orires(int                       nfa,
            const t_iatom             forceatoms[],
            const t_iparams           ip[],
            const rvec                x[],
            rvec4                     f[],
            rvec                      fshift[],
            const t_pbc*              pbc,
            real                      lambda,
            real*                     dvdlambda,
            gmx::ArrayRef<const real> charge,
            t_fcdata*                 fcd,
            t_disresdata*             disresdata,
            t_oriresdata*             oriresdata,
            int*                      global_atom_index);

#endif
