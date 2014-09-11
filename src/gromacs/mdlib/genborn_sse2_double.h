/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
#ifndef _genborn_sse2_double_h
#define _genborn_sse2_double_h

#include "gromacs/legacyheaders/typedefs.h"

int
calc_gb_rad_still_sse2_double(t_commrec *cr, t_forcerec *fr, int natoms, gmx_localtop_t *top,
                              double *x, t_nblist *nl, gmx_genborn_t *born);

int
calc_gb_chainrule_sse2_double(int natoms, t_nblist *nl, double *dadx, double *dvda, double *xd, double *f,
                              double *fshift, double *shift_vec, int gb_algorithm,
                              gmx_genborn_t *born, t_mdatoms *md);

int
calc_gb_rad_hct_obc_sse2_double(t_commrec *cr, t_forcerec *fr, int natoms, gmx_localtop_t *top,
                                double *x, t_nblist *nl, gmx_genborn_t *born, t_mdatoms *md, int gb_algorithm);

#endif /* _genborn_sse2_double_h */
