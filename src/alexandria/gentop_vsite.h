/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: poldata.h,v 1.13 2009/02/03 16:08:37 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.99
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

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
 * Groningen Machine for Chemical Simulation
 */

#ifndef _gentop_vsite_h
#define _gentop_vsite_h

#include "grompp.h"
#include "poldata.h"

enum { egvtNO, egvtLINEAR, egvtPLANAR, egvtALL, egvtNR };

typedef struct gentop_vsite *gentop_vsite_t;

extern gentop_vsite_t gentop_vsite_init(int egvt);

extern void gentop_vsite_done(gentop_vsite_t *gvt);

extern void gentop_vsite_add_linear(gentop_vsite_t gvt,int a1,int a2,int a3);

extern void gentop_vsite_add_planar(gentop_vsite_t gvt,int a1,int a2,int a3,int a4,int nbonds[]);

extern void gentop_vsite_generate_special(gentop_vsite_t gvt,gmx_bool bGenVsites,
                                          t_atoms *atoms,rvec **x,
                                          t_params plist[],
                                          t_symtab *symtab,gpp_atomtype_t atype,
                                          t_excls **excls,gmx_poldata_t pd);

#endif
