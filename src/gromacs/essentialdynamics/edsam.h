/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#include "gromacs/fileio/filenm.h"
#include "gromacs/legacyheaders/typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

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
void do_edsam(t_inputrec *ir, gmx_int64_t step,
              t_commrec *cr, rvec xs[], rvec v[], matrix box, gmx_edsam_t ed);


/*! \brief Reads in the .edi file containing the essential dynamics and flooding data.
 *
 * This function opens the ED input and output files, reads in all datasets it finds
 * in the input file, and cross-checks whether the .edi file information is consistent
 * with the essential dynamics data found in the checkpoint file (if present).
 * gmx make_edi can be used to create an .edi input file.
 *
 * \param natoms            Number of atoms of the whole MD system.
 * \param EDstate           Essential dynamics and flooding data stored in the checkpoint file.
 * \param nfile             Number of entries (files) in the fnm structure.
 * \param fnm               The filenames struct; it contains also the names of the
 *                          essential dynamics and flooding in + output files.
 * \param Flags             Flags passed over from main, used to determine
 *                          whether we are appending.
 * \param oenv              Needed to open the output xvgr file.
 * \param cr                Data needed for MPI communication.
 * \returns                 Pointer to the initialized essential dynamics / flooding data.
 */
gmx_edsam_t ed_open(int natoms, edsamstate_t *EDstate, int nfile, const t_filenm fnm[],
                    unsigned long Flags, const output_env_t oenv, t_commrec *cr);


/*! \brief Initializes the essential dynamics and flooding module.
 *
 * \param mtop              Molecular topology.
 * \param ir                MD input parameter record.
 * \param cr                Data needed for MPI communication.
 * \param ed                The essential dynamics data.
 * \param x                 Positions of the whole MD system.
 * \param box               The simulation box.
 * \param EDstate           ED data stored in the checkpoint file.
 */
void init_edsam(gmx_mtop_t *mtop, t_inputrec *ir, t_commrec *cr,
                gmx_edsam_t ed, rvec x[], matrix box, edsamstate_t *EDstate);


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
void do_flood(t_commrec *cr, t_inputrec *ir, rvec x[], rvec force[], gmx_edsam_t ed,
              matrix box, gmx_int64_t step, gmx_bool bNS);

/*! \brief Clean up
 *
 * \param ed                The essential dynamics data
 */
void done_ed(gmx_edsam_t *ed);

#ifdef __cplusplus
}
#endif

#endif
