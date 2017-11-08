/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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
#ifndef GMX_MDLIB_MDSETUP_H
#define GMX_MDLIB_MDSETUP_H

#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/shellfc.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/topology/topology.h"

/*! \brief Sets atom data for several MD algorithms
 *
 * Most MD algorithms require two different setup calls:
 * one for initialization and parameter setting and one for atom data setup.
 * This routine sets the atom data for the (locally available) atoms.
 * This is called at the start of serial runs and during domain decomposition.
 *
 * \param[in]     cr         Communication record
 * \param[in]     ir         Input parameter record
 * \param[in]     top_global The global topology
 * \param[in,out] top        The local topology
 * \param[in,out] fr         The force calculation parameter/data record
 * \param[out]    graph      The molecular graph, can be NULL
 * \param[out]    mdAtoms    The MD atom data
 * \param[in,out] vsite      The virtual site data, can be NULL
 * \param[in,out] shellfc    The shell/flexible-constraint data, can be NULL
 */
void mdAlgorithmsSetupAtomData(t_commrec         *cr,
                               const t_inputrec  *ir,
                               const gmx_mtop_t  *top_global,
                               gmx_localtop_t    *top,
                               t_forcerec        *fr,
                               t_graph          **graph,
                               gmx::MDAtoms      *mdAtoms,
                               gmx_vsite_t       *vsite,
                               gmx_shellfc_t     *shellfc);

#endif
