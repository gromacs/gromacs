/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.5
 * Written by Christoph Junghans, Brad Lambeth, and possibly others.
 * Copyright (c) 2009 Christoph Junghans, Brad Lambeth.
 * All rights reserved.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * GROningen Mixture of Alchemy and Childrens' Stories
 */

#ifndef _adress_h_
#define _adress_h_

/** \file adress.h
 *
 * \brief Implementation of the AdResS method
 *
 */

#include "types/simple.h"
#include "typedefs.h"

/** \brief calculates the AdResS weight of a particle
 *
 * \param[in] x position of the particle
 * \param[in] adresstype type of address weight function
 *                       1 - constant weight all over the box
 *                       2 - split in x direction with ref as center
 *                       3 - spherical splitting with ref as center
 *                       4 - sphere moving along the last charge group
 *                       else - weight = 1 - explicit simulation
 * \param[in] adressr radius/size of the explicit zone
 * \param[in] adressw size of the hybrid zone 
 *                    for adresstype 1, the value of the weight
 * \param[in] bnew_wf new weighting function or not
 * \param[in] ref center of the explicit zone
 *                for adresstype 1 - unused
 *                for adresstype 2 - only ref[0] is used
 * \param[in] box2 array[3] containing half lengths of the rectangular box
 *                 only used by adresstype 4
 * \param[in] box matrix containing lengths of the box
 *                 only used by adresstype 4
 * \return weight of the particle
 *
 * \todo calc AdResS weight for non-rectangular boxes 
 */
real 
adress_weight(rvec             x,
              int              adresstype,
              real             adressr,
              real             adressw,
              bool             bnew_wf,
              rvec *           ref,
              rvec             box2,
              matrix           box);

/** \brief update the weight of all coarse-grained particles in several charge groups for com vsites
 *
 * \param[in,out] fplog log file in case of debug
 * \param[in] cg0 first charge group to update
 * \param[in] cg1 last+1 charge group to update
 * \param[in] cgs block containing the cg index 
 * \param[in] x array with all the particle positions  
 * \param[in] fr the frocerec containing all the parameters
 * \param[in,out] mdatoms the struct containing all the atoms properties
 * \param[in] box matrix containing lengths of the box
 */
void
update_adress_weights_com(FILE *               fplog,
                          int                  cg0,
                          int                  cg1,
                          t_block *            cgs,
                          rvec                 x[],
                          t_forcerec *         fr,
                          t_mdatoms *          mdatoms,
                          matrix               box);

/** \brief update the weight of all coarse-grained particles for cog vsites
 *
 * \param[in] ip contains interaction parameters, in this case the number of constructing atoms n for vsitesn and vsitesrefmol 
 * \param[in] ilist list of interaction types, in this case the virtual site types are what's important
 * \param[in] x array with all the particle positions  
 * \param[in] fr the frocerec containing all the parameters
 * \param[in,out] mdatoms the struct containing all the atoms properties
 * \param[in] box matrix containing lengths of the box
 */
void
update_adress_weights_cog(t_iparams            ip[],
                          t_ilist              ilist[],
                          rvec                 x[],
                          t_forcerec *         fr,
                          t_mdatoms *          mdatoms,
                          matrix               box);

#endif
