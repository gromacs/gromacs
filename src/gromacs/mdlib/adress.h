/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009 Christoph Junghans, Brad Lambeth.
 * Copyright (c) 2011,2014, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief Implementation of the AdResS method.
 */
#ifndef _adress_h_
#define _adress_h_

#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/simple.h"

#ifdef __cplusplus
extern "C"
{
#endif

struct t_pbc;

/** \brief calculates the AdResS weight of a particle
 *
 * \param[in] x position of the particle
 * \param[in] adresstype type of address weight function
 *                       eAdressOff - explicit simulation
 *                       eAdressConst - constant weight all over the box
 *                       eAdressXSplit - split in x direction with ref as center
 *                       eAdressSphere - spherical splitting with ref as center
 *                       else - weight = 1 - explicit simulation
 * \param[in] adressr radius/size of the explicit zone
 * \param[in] adressw size of the hybrid zone
 * \param[in] ref center of the explicit zone
 *                for adresstype 1 - unused
 *                for adresstype 2 - only ref[0] is used
 * \param[in] pbc pbc struct for calculating shortest distance
 * \param[in] fr the forcerec containing all the parameters
 *
 * \return weight of the particle
 *
 */
real
adress_weight(rvec                 x,
              int                  adresstype,
              real                 adressr,
              real                 adressw,
              rvec     *           ref,
              struct t_pbc     *   pbc,
              t_forcerec *         fr);

/** \brief update the weight of all coarse-grained particles in several charge groups for com vsites
 *
 * \param[in,out] fplog log file in case of debug
 * \param[in] cg0 first charge group to update
 * \param[in] cg1 last+1 charge group to update
 * \param[in] cgs block containing the cg index
 * \param[in] x array with all the particle positions
 * \param[in] fr the forcerec containing all the parameters
 * \param[in,out] mdatoms the struct containing all the atoms properties
 * \param[in] pbc for shortest distance in adress_weight
 */
void
update_adress_weights_com(FILE *               fplog,
                          int                  cg0,
                          int                  cg1,
                          t_block *            cgs,
                          rvec                 x[],
                          t_forcerec *         fr,
                          t_mdatoms *          mdatoms,
                          struct t_pbc *       pbc);

/** \brief update the weight of all coarse-grained particles for cog vsites
 *
 * \param[in] ip contains interaction parameters, in this case the number of constructing atoms n for vsitesn
 * \param[in] ilist list of interaction types, in this case the virtual site types are what's important
 * \param[in] x array with all the particle positions
 * \param[in] fr the forcerec containing all the parameters
 * \param[in,out] mdatoms the struct containing all the atoms properties
 * \param[in] pbc for shortest distance in adress_weight
 */
void
update_adress_weights_cog(t_iparams            ip[],
                          t_ilist              ilist[],
                          rvec                 x[],
                          t_forcerec *         fr,
                          t_mdatoms *          mdatoms,
                          struct t_pbc *       pbc);

/** \brief update the weight of all coarse-grained particles in several charge groups for atom vsites
 *
 * \param[in] cg0 first charge group to update
 * \param[in] cg1 last+1 charge group to update
 * \param[in] cgs block containing the cg index
 * \param[in] x array with all the particle positions
 * \param[in] fr the forcerec containing all the parameters
 * \param[in,out] mdatoms the struct containing all the atoms properties
 * \param[in] pbc for shortest distance in adress_weight
 */
void
update_adress_weights_atom(int                  cg0,
                           int                  cg1,
                           t_block *            cgs,
                           rvec                 x[],
                           t_forcerec *         fr,
                           t_mdatoms *          mdatoms,
                           struct t_pbc *       pbc);

/** \brief update the weight on per atom basis of all coarse-grained particles in several charge groups for atom vsites
 *
 * \param[in] cg0 first charge group to update
 * \param[in] cg1 last+1 charge group to update
 * \param[in] cgs block containing the cg index
 * \param[in] x array with all the particle positions
 * \param[in] fr the forcerec containing all the parameters
 * \param[in,out] mdatoms the struct containing all the atoms properties
 * \param[in] pbc for shortest distance in adress_weight
 */
void
update_adress_weights_atom_per_atom(int                  cg0,
                                    int                  cg1,
                                    t_block *            cgs,
                                    rvec                 x[],
                                    t_forcerec *         fr,
                                    t_mdatoms *          mdatoms,
                                    struct t_pbc *       pbc);

/** \brief add AdResS IC thermodynamic force to f_novirsum
 *
 * \param[in] cg0 first charge group to update
 * \param[in] cg1 last+1 charge group to update
 * \param[in] cgs block containing the cg index
 * \param[in] x array with all the particle positions
 * \param[in,out] f the force array pointing at f_novirsum from sim_util.c
 * \param[in] fr the forcerec containing all the parameters
 * \param[in] mdatoms the struct containing all the atoms properties
 * \param[in] pbc for shortest distance to fr->adress_refs
 */
void
adress_thermo_force(int                  cg0,
                    int                  cg1,
                    t_block *            cgs,
                    rvec                 x[],
                    rvec                 f[],
                    t_forcerec *         fr,
                    t_mdatoms *          mdatoms,
                    struct t_pbc *       pbc);


/** \brief checks weather a cpu calculates only coarse-grained or explicit interactions
 *
 * \param[in] n_ex number of explicit particles
 * \param[in] n_hyb number of hybrid particles
 * \param[in] n_cg number of coarse-grained particles
 * \param[in,out] mdatoms the struct containing all the atoms properties
 */
void adress_set_kernel_flags(int n_ex, int n_hyb, int n_cg, t_mdatoms * mdatoms);

/** \brief looks up  if a energy group is explicit
 * \param[in] fr the forcerec containing all the parameters
 * \param[in] egp_nr energy group number
 * \return boolean if explicit or not
 */
gmx_bool egp_explicit(t_forcerec *   fr, int egp_nr);

/** \brief looks up  if a energy group is coarse-grained
 * \param[in] fr the forcerec containing all the parameters
 * \param[in] egp_nr energy group number
 * \return boolean if coarse-grained or not
 */
gmx_bool egp_coarsegrained(t_forcerec *   fr, int egp_nr);

#ifdef __cplusplus
}
#endif

#endif
