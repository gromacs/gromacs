/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * Grunge ROck MAChoS
 */

#ifndef _axp_asm_h
#define _axp_asm_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* Function definitions for alpha AXP assembly routines */

/* Double precision inverse square root */
void sqrtiv_(double *indata, double *outdata, int *n);

/* Double precision square root */
void sqrtv_(double *indata, double *outdata, int *n);

/* Single precision inverse square root */
void ssqrtiv_(float *indata, float *outdata, int *n);

/* Single precision square root */
void ssqrtv_(float *indata, float *outdata, int *n);

#endif

