/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Good ROcking Metal Altar for Chronical Sinners
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define D_BOX_Z 1
#ifdef ALLOW_OFFDIAG_LT_HALFDIAG
#define D_BOX_Y 2
#define D_BOX_X 2
#else
#define D_BOX_Y 1
#define D_BOX_X 1
#endif 
#define N_BOX_Z (2*D_BOX_Z+1)
#define N_BOX_Y (2*D_BOX_Y+1)
#define N_BOX_X (2*D_BOX_X+1)
#define N_IVEC  (N_BOX_Z*N_BOX_Y*N_BOX_X)
#define CENTRAL (N_IVEC/2)
#define SHIFTS  N_IVEC

#define XYZ2IS(x,y,z) (N_BOX_X*(N_BOX_Y*((z)+D_BOX_Z)+(y)+D_BOX_Y)+(x)+D_BOX_X)
#define IVEC2IS(iv)   (XYZ2IS((iv)[XX],(iv)[YY],(iv)[ZZ]))
