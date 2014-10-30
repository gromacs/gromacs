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

#ifdef __cplusplus
extern "C" {
#endif

/* Iniatiate the non-blocking receive of the halo coordinates x.
 * Should be called as early as possible for fast communication.
 */
void dd_halo_initiate_recv_x(gmx_domdec_t *dd, rvec x[]);

/* Initiate the non-blocking send of the halo coordinates x.
 * Should be called as early as possible when x is ready.
 */
void dd_halo_initiate_send_x(gmx_domdec_t *dd, matrix box, rvec x[]);

/* Ensures that x has been received.
 * Needs to be called after dd_halo_initiate_recv_x. Call as late as possible.
 */
void dd_halo_complete_recv_x(gmx_domdec_t *dd);

/* Ensures that x has been sent.
 * Needs to be called after dd_halo_initiate_send_x. Call as late as possible.
 */
void dd_halo_complete_send_x(gmx_domdec_t *dd);

/* Communicates the halo coordinates x.
 * Internally calls the 4 functions above.
 */
void dd_halo_move_x(gmx_domdec_t *dd, matrix box, rvec x[]);

/* Iniatiate the non-blocking receive of the halo forces f.
 * Should be called as early as possible for fast communication.
 */
void dd_halo_initiate_recv_f(gmx_domdec_t *dd);

/* Initiate the non-blocking send of the halo forces f.
 * Should be called as early as possible when f is ready.
 */
void dd_halo_initiate_send_f(gmx_domdec_t *dd, rvec *f);

/* Ensures that x has been received.
 * Needs to be called after dd_halo_initiate_recv_x. Call as late as possible.
 */
void dd_halo_complete_recv_x(gmx_domdec_t *dd);

/* Ensures that x has been sent.
 * Needs to be called after dd_halo_initiate_send_x. Call as late as possible.
 */
void dd_halo_complete_send_x(gmx_domdec_t *dd);

/* Communicates the halo coordinates x.
 * Internally calls the 4 functions above.
 */
void dd_halo_move_x(gmx_domdec_t *dd, matrix box, rvec *x);

/* Iniatiate the non-blocking receive of the halo forces f.
 * Should be called as early as possible for fast communication.
 */
void dd_halo_initiate_recv_f(gmx_domdec_t *dd);

/* Initiate the non-blocking send of the halo forces f.
 * Should be called as early as possible when f is ready.
 */
void dd_halo_initiate_send_f(gmx_domdec_t *dd, rvec *f);

/* Ensures that f has been received and reduced the received forces.
 * If fshift!=NULL, also updates the shift forces.
 * Needs to be called after dd_halo_initiate_recv_f. Call as late as possible.
 */
void dd_halo_complete_recv_f(gmx_domdec_t *dd, rvec *f, rvec *fshift);

/* Ensures that f has been sent.
 * Needs to be called after dd_halo_initiate_send_f. Call as late as possible.
 */
void dd_halo_complete_send_f(gmx_domdec_t *dd);

/* Communicates the halo coordinates f.
 * Internally calls the 4 functions above.
 */
void dd_halo_move_f(gmx_domdec_t *dd, rvec *f, rvec *fshift);

/* Set up the eighth-shell halo communcation for non-bonded + bonded int.
 * bCellsChanged indicates if the domain deomposition cells changed,
 * either due to dynamic load balancing or pressure scaling.
 */
void setup_halo_communication(gmx_domdec_t *dd,
                              const matrix box, const gmx_ddbox_t *ddbox,
                              t_forcerec *fr,
                              gmx_bool bCellsChanged);

#ifdef __cplusplus
}
#endif
