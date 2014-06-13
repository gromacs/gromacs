/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
/*! \file
 * \brief
 * Declares code for calling GROMACS routines from an external program
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \inpublicapi
 * \ingroup module_qmmm
 */
#ifndef GMX_MMSLAVE_H
#define GMX_MMSLAVE_H

#include "gromacs/legacyheaders/types/simple.h"
#include "gromacs/legacyheaders/network.h"

#ifdef __cplusplus
extern "C" {
#endif

//! Abstract type for the mmslave code
typedef struct gmx_mmslave *gmx_mmslave_t;

/*! \brief
 * Create the data structure and returns it
 * \param[in] gms  the data structure to be freed
 */
extern gmx_mmslave_t mmslave_init(const t_commrec *cr);

/*! \brief
 * Function to clean up
 * \param[in] gms  the data structure to be freed
 */
extern void mmslave_done(gmx_mmslave_t gms);

/*! \brief
 * Function to read a tpr file and store the information in the
 * data structure
 * \param[in] tpr  the file name
 * \param[in] gms  the data structure to be read
 * \return 1 if successful, 0 otherwise
 */
extern int mmslave_read_tpr(const char   *tpr,
                            gmx_mmslave_t gms);

/*! \brief
 * The number of groups are numbered 0 to N-1. The group
 * index of each atom can be extracted using mmslave_get_group_id
 * index in subsequent calls to related routines.
 * \param[in] gms  the data structure
 * \return the number N of atom-groups in the system
 */
extern int mmslave_ngroups(gmx_mmslave_t gms);

/*! \brief
 * Return number of atoms in specified group
 * \param[in] gms    the data structure
 * \return the number of atoms in group or 0 if group is out of range
 */
extern int mmslave_natoms(gmx_mmslave_t gms);

/*! \brief
 * Copy internal coordinate array
 * \param[in]  gms    the data structure
 * \param[in]  natoms length of array
 * \param[out] x      array of rvec (must be allocated by the caller)
 * \return 1 on success, 0 otherwise (typically if natoms is too small or x is NULL)
 */
extern int mmslave_copyX(gmx_mmslave_t gms,
                         int           natoms,
                         rvec         *x);

/*! \brief
 * Copy internal velocity array
 * \param[in]  gms    the data structure containing the data
 * \param[in]  natoms length of array
 * \param[out] v      array of rvecs (must be allocated by the caller)
 * \return 1 on success, 0 otherwise (typically if natoms is too small or v is NULL)
 */
extern int mmslave_copyV(gmx_mmslave_t gms,
                         int           natoms,
                         rvec         *v);

/*! \brief
 * Copy internal force array
 * \param[in]  gms    the data structure containing the data
 * \param[in]  natoms length of array
 * \param[out] f      array of rvecs (must be allocated by the caller)
 * \return 1 on success, 0 otherwise (typically if natoms is too small or f is NULL)
 */
extern int mmslave_copyF(gmx_mmslave_t gms,
                         int           natoms,
                         rvec         *f);

/*! \brief
 * Function to cleanup the innards of the data structure
 */
extern void mmslave_clean(gmx_mmslave_t gms);

/*! \brief
 * Function to modify the charge of an atom
 * data structure
 * \param[in]  gms    the data structure to be modified
 * \param[in]  id     the atom id
 * \param[in]  q      the new charge
 * \return 1 if successful, 0 otherwise
 */
extern int mmslave_set_q(gmx_mmslave_t gms,
                         atom_id       id,
                         double        q);

/*! \brief
 * Function to retrieve the charge of an atom
 * data structure
 * \param[in] gms    the data structure to be returned
 * \param[in] id     the atom id
 * \param[out] q     the charge of the atom
 * \return 1 on success 0 otherwise
 */
extern int mmslave_get_q(gmx_mmslave_t  gms,
                         atom_id        id,
                         double        *q);

/*! \brief
 * Function to retrieve the atom number of an atom
 * data structure
 * \param[in] gms   the data structure to be returned
 * \param[in] id    the atom id
 * \return the atom number or NOTSET
 */
extern int mmslave_get_atomnumber(gmx_mmslave_t gms,
                                  atom_id       id);

/*! \brief
 * Function to retrieve the group ID of an atom
 * data structure
 * \param[in] gms   the data structure to be returned
 * \param[in] id    the atom id
 * \return the group_id or NOTSET
 */
extern int mmslave_get_group_id(gmx_mmslave_t gms,
                                atom_id       id);

/*! \brief
 * Function to compute the energy and forces
 * \param[in]  gms    the data structure to be modified
 * \param[in]  fplog  a File pointer for (debug) output. May be NULL, or stdout.
 * \param[in]  x      the atomic coordinates for the whole system (MM+QM)
 * \param[out] f      the forces on all atoms
 * \param[out] A      the electric field on all atoms
 * \param[out] phi    the electrostatic potential on all atoms
 * \param[out] energy the total MM energy
 * \return 1 if successful, 0 otherwise
 */
extern int mmslave_calc_energy(gmx_mmslave_t gms,
                               FILE         *fplog,
                               const rvec   *x,
                               rvec         *f,
                               rvec         *A,
                               real         *phi,
                               double       *energy);


#ifdef __cplusplus
}
#endif

#endif
