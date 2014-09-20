/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team.
 * Copyright (c) 2010,2014, by the GROMACS development team, led by
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
#ifndef _GENBORN_ALLVSALL_SSE2_DOUBLE_H
#define _GENBORN_ALLVSALL_SSE2_DOUBLE_H

#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/simple.h"

int
genborn_allvsall_calc_still_radii_sse2_double(t_forcerec *           fr,
                                              t_mdatoms *            mdatoms,
                                              gmx_genborn_t *        born,
                                              gmx_localtop_t *       top,
                                              double *               x,
                                              t_commrec *            cr,
                                              void *                 work);

int
genborn_allvsall_calc_hct_obc_radii_sse2_double(t_forcerec *           fr,
                                                t_mdatoms *            mdatoms,
                                                gmx_genborn_t *        born,
                                                int                    gb_algorithm,
                                                gmx_localtop_t *       top,
                                                double *               x,
                                                t_commrec *            cr,
                                                void *                 work);

int
genborn_allvsall_calc_chainrule_sse2_double(t_forcerec *           fr,
                                            t_mdatoms *            mdatoms,
                                            gmx_genborn_t *        born,
                                            double *               x,
                                            double *               f,
                                            int                    gb_algorithm,
                                            void *                 work);

#endif
