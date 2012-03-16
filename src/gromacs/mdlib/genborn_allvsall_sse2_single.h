/*
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team
 *
 * Gromacs is a library for molecular simulation and trajectory analysis,
 * written by Erik Lindahl, David van der Spoel, Berk Hess, and others - for
 * a full list of developers and information, check out http://www.gromacs.org
 *
 * This program is free software; you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the Free 
 * Software Foundation; either version 2 of the License, or (at your option) any 
 * later version.
 * As a special exception, you may use this file as part of a free software
 * library without restriction.  Specifically, if other files instantiate
 * templates or use macros or inline functions from this file, or you compile
 * this file and link it with other files to produce an executable, this
 * file does not by itself cause the resulting executable to be covered by
 * the GNU Lesser General Public License.  
 *
 * In plain-speak: do not worry about classes/macros/templates either - only
 * changes to the library have to be LGPL, not an application linking with it.
 *
 * To help fund GROMACS development, we humbly ask that you cite
 * the papers people have written on it - you can find them on the website!
 */
#ifndef _GENBORN_ALLVSALL_SSE2_SINGLE_H
#define _GENBORN_ALLVSALL_SSE2_SINGLE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "types/simple.h"
#include "typedefs.h"

int
genborn_allvsall_calc_still_radii_sse2_single(t_forcerec *           fr,
                                              t_mdatoms *            mdatoms,
                                              gmx_genborn_t *        born,
                                              gmx_localtop_t *       top,
                                              real *                 x,
                                              t_commrec *            cr,
                                              void *                 work);

int
genborn_allvsall_calc_hct_obc_radii_sse2_single(t_forcerec *           fr,
                                                t_mdatoms *            mdatoms,
                                                gmx_genborn_t *        born,
                                                int                    gb_algorithm,
                                                gmx_localtop_t *       top,
                                                real *                 x,
                                                t_commrec *            cr,
                                                void *                 work);

int
genborn_allvsall_calc_chainrule_sse2_single(t_forcerec *           fr,
                                            t_mdatoms *            mdatoms,
                                            gmx_genborn_t *        born,
                                            real *                 x,
                                            real *                 f,
                                            int                    gb_algorithm,
                                            void *                 work);

#endif

