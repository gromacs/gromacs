/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2005,2006,2007,2008,2009,2010,2011,2012,2013,2014, by the GROMACS development team, led by
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
 *  \brief
 * Halo communication for the force calculation.
 *
 * This file contains functions to set up the halo communication
 * for the eighth shell domain decomposition and to execute
 * the halo communication.
 *
 * \inlibraryapi
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_DOMDEC_HALO_H
#define GMX_DOMDEC_DOMDEC_HALO_H

#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/legacyheaders/types/forcerec.h"

/*! \brief Initiate a non-blocking halo coordinate receive.
 *
 * Iniatiate the non-blocking receive of the halo coordinates \p x.
 * Should be called as early as possible for fast communication.
 *
 * \param[in]  dd Pointer to the domain decomposition setup.
 * \param[out] x  Pointer to the start of the coordinate array.
 */
void dd_halo_initiate_recv_x(gmx_domdec_t *dd, rvec *x);

/*! \brief Initiate a non-blocking halo coordinate send.
 *
 * Initiate the non-blocking send of the halo coordinates \p x.
 * Should be called as early as possible when \p x is ready.
 *
 * \param[in]  dd  Pointer to the domain decomposition setup.
 * \param[in]  box Unit-cell matrix.
 * \param[in]  x   Pointer to the start of the coordinate array.
 */
void dd_halo_initiate_send_x(gmx_domdec_t *dd, matrix box, rvec x[]);

/*! \brief Ensures that the coordinates have been received.
 *
 * Needs to be called after dd_halo_initiate_recv_x.
 * Call as late as possible.
 *
 * \param[in]  dd Pointer to the domain decomposition setup.
 */
void dd_halo_complete_recv_x(gmx_domdec_t *dd);

/*! \brief Ensures that the coordinates have been received.
 *
 * Needs to be called after dd_halo_initiate_send_x.
 * Call as late as possible.
 *
 * \param[in]  dd Pointer to the domain decomposition setup.
 */
void dd_halo_complete_send_x(gmx_domdec_t *dd);

/*! \brief Send and receive the halo coordinates.
 *
 * Internally calls the 4 functions above.
 *
 * \param[in]     dd  Pointer to the domain decomposition setup.
 * \param[in]     box Unit-cell matrix.
 * \param[in,out] x   Pointer to the start of the coordinate array.
 */
void dd_halo_move_x(gmx_domdec_t *dd, matrix box, rvec *x);

/*! \brief Initiate a non-blocking halo force receive.
 *
 * Iniatiate the non-blocking receive of the halo forces.
 * Should be called as early as possible for fast communication.
 *
 * \param[in]  dd Pointer to the domain decomposition setup.
 */
void dd_halo_initiate_recv_f(gmx_domdec_t *dd);

/*! \brief Initiate a non-blocking halo force send.
 *
 * Initiate the non-blocking send of the halo forces \p f
 * Should be called as early as possible when \p f is ready.
 *
 * \param[in]  dd  Pointer to the domain decomposition setup.
 * \param[in]  f   Pointer to the start of the force array.
 */
void dd_halo_initiate_send_f(gmx_domdec_t *dd, rvec *f);

/*! \brief Ensures that the forces have been received and reduces the forces.
 *
 * If \p fshift!=NULL, also updates the shift forces \p fshift.
 * Needs to be called after dd_halo_initiate_recv_f.
 * Call as late as possible.
 *
 * \param[in]     dd     Pointer to the domain decomposition setup.
 * \param[in,out] f      Pointer to the start of the force array.
 * \param[in,out] fshift Pointer to the start of the shift force array.
 */
void dd_halo_complete_recv_f(gmx_domdec_t *dd, rvec *f, rvec *fshift);

/*! \brief Completes the force send.
 *
 * Needs to be called after dd_halo_initiate_send_f.
 * Call as late as possible.
 *
 * \param[in]  dd Pointer to the domain decomposition setup.
 */
void dd_halo_complete_send_f(gmx_domdec_t *dd);

/*! \brief Communicates and reduces the halo forces.
 *
 * If \p fshift!=NULL, also updates the shift forces.
 * Internally calls the 4 functions above.
 *
 * \param[in]     dd     Pointer to the domain decomposition setup.
 * \param[in,out] f      Pointer to the start of the force array.
 * \param[in,out] fshift Pointer to the start of the shift force array.
 */
void dd_halo_move_f(gmx_domdec_t *dd, rvec *f, rvec *fshift);

/*! \brief Set up the eighth-shell halo coordinate/force communcation
 *
 * Set up the eighth-shell halo communcation for non-bonded + bonded
 * interactions
 * \p bCellsChanged indicates if the domain decomposition cells changed,
 * either due to dynamic load balancing or pressure scaling.
 *
 * \param[in,out] dd     Pointer to the domain decomposition setup.
 * \param[in]     box    Unit-cell matrix.
 * \param[in]     ddbox  Domain decomposition unit-cell and PBC data.
 * \param[in]     fr     Pointr to the force record struct.
 * \param[in]     bCellsChanged Tells if the domain decomposition cells changed.
 */
void setup_halo_communication(gmx_domdec_t *dd,
                              const matrix box, const gmx_ddbox_t *ddbox,
                              t_forcerec *fr,
                              gmx_bool bCellsChanged);

#endif
