/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
 *
 * \brief
 * Declares functions to calculate both essential dynamics constraints
 * as well as flooding potentials and forces.
 *
 * \authors Bert de Groot <bgroot@gwdg.de>, Oliver Lange <oliver.lange@tum.de>,
 * Carsten Kutzner <ckutzne@gwdg.de>
 *
 * \inlibraryapi
 */
#ifndef GMX_ESSENTIALDYNAMICS_EDSAM_H
#define GMX_ESSENTIALDYNAMICS_EDSAM_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"


/*! \brief Abstract type for essential dynamics
 *
 * The main type is defined only in edsam.cpp
 */
typedef struct gmx_edsam *gmx_edsam_t;

struct gmx_domdec_t;
struct gmx_mtop_t;
struct gmx_output_env_t;
struct ObservablesHistory;
struct t_commrec;
struct t_filenm;
struct t_inputrec;
class t_state;

/*! \brief Applies essential dynamics constrains as defined in the .edi input file.
 *
 * \param ir                MD input parameter record.
 * \param step              Number of the time step.
 * \param cr                Data needed for MPI communication.
 * \param xs                The local positions on this processor.
 * \param v                 The local velocities.
 * \param box               The simulation box.
 * \param ed                The essential dynamics data.
 */
void do_edsam(const t_inputrec *ir, gmx_int64_t step,
              t_commrec *cr, rvec xs[], rvec v[], matrix box, gmx_edsam_t ed);


/*! \brief Initializes the essential dynamics and flooding module.
 *
 * \param ediFileName       Essential dynamics input file.
 * \param edoFileName       Output file for essential dynamics data.
 * \param mtop              Molecular topology.
 * \param ir                MD input parameter record.
 * \param cr                Data needed for MPI communication.
 * \param constr            Data structure keeping the constraint information.
 * \param globalState       The global state, only used on the master rank.
 * \param oh                The observables history container.
 * \param oenv              The output environment information.
 * \param bAppend           Append to existing output files?
 *
 * \returns                 A pointer to the ED data structure.
 */
gmx_edsam_t init_edsam(
        const char             *ediFileName,
        const char             *edoFileName,
        const gmx_mtop_t       *mtop,
        const t_inputrec       *ir,
        t_commrec              *cr,
        struct gmx_constr      *constr,
        const t_state          *globalState,
        ObservablesHistory     *oh,
        const gmx_output_env_t *oenv,
        gmx_bool                bAppend);

/*! \brief Make a selection of the home atoms for the ED groups.
 *
 * Should be called at every domain decomposition.
 *
 * \param dd                Domain decomposition data.
 * \param ed                Essential dynamics and flooding data.
 */
void dd_make_local_ed_indices(gmx_domdec_t *dd, gmx_edsam_t ed);


/*! \brief Evaluate the flooding potential(s) and forces as requested in the .edi input file.
 *
 * \param cr                Data needed for MPI communication.
 * \param ir                MD input parameter record.
 * \param x                 Positions on the local processor.
 * \param force             Forcefield forces to which the flooding forces are added.
 * \param ed                The essential dynamics data.
 * \param box               The simulation box.
 * \param step              Number of the time step.
 * \param bNS               Are we in a neighbor searching step?
 */
void do_flood(t_commrec *cr, const t_inputrec *ir, rvec x[], rvec force[], gmx_edsam_t ed,
              matrix box, gmx_int64_t step, gmx_bool bNS);

/*! \brief Clean up
 *
 * \param ed                The essential dynamics data
 */
void done_ed(gmx_edsam_t *ed);

#endif
