/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
 * GRoups of Organic Molecules in ACtion for Science
 */

#ifdef __cplusplus
extern "C" {
#endif

#include "simple.h"

typedef struct {
  real r,g,b;
} t_rgb;

typedef struct {
  char c1; /* should all be non-zero (and printable and not '"') */
  char c2; /* 
	    * should all be zero (single char color names: smaller xpm's) 
	    * or should all be non-zero (double char color names: more colors)
	    */
} t_xpmelmt;

typedef short t_matelmt;

typedef struct {
  t_xpmelmt code; /* see comment for t_xpmelmt */
  const char *desc;
  t_rgb rgb;
} t_mapping;

#define MAT_SPATIAL_X (1<<0)			    
#define MAT_SPATIAL_Y (1<<1)
/* Defines if x and y are spatial dimensions,
 * when not, there are n axis ticks at the middle of the elements,
 * when set, there are n+1 axis ticks at the edges of the elements.
 */

typedef struct {
  unsigned int flags; /* The possible flags are defined above */
  int  nx,ny;
  int  y0;
  char title[256];
  char legend[256];
  char label_x[256];
  char label_y[256];
  gmx_bool bDiscrete;
  real *axis_x;
  real *axis_y;
  t_matelmt **matrix;
  int nmap;
  t_mapping *map;
} t_matrix;
  /* title      matrix title
   * legend     label for the continuous legend
   * label_x    label for the x-axis
   * label_y    label for the y-axis
   * nx, ny     size of the matrix
   * axis_x[]   the x-ticklabels
   * axis_y[]   the y-ticklables
   * *matrix[]  element x,y is matrix[x][y]
   * nmap       number of color levels for the output(?)
   */

#ifdef __cplusplus
}
#endif

