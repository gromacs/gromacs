/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

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
  char *desc;
  t_rgb rgb;
} t_mapping;

typedef struct {
  int  nx,ny;
  int  y0;
  char title[256];
  char legend[256];
  char label_x[256];
  char label_y[256];
  bool bDiscrete;
  real *axis_x;
  real *axis_y;
  t_matelmt **matrix;
  int nmap;
  t_mapping *map;
} t_matrix;
