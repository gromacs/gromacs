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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

typedef struct {
  int  nr;			/* Row- and Column-length of vectors	*/
  int  nrfp;			/* Number of force parameters		*/
  real ***c;		 	/* Matrix with arrays with force	*/
				/* constants (dim: nr x nr x nrfp)	*/
} t_nbs;

typedef struct
{
  int		maxbt;		/* The number of nbs types	*/
  t_nbs		*nbs;		/* The array of nbs structs	*/
} t_nbslist;

