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
 * Getting the Right Output Means no Artefacts in Calculating Stuff
 */

#ifndef _sheader_h
#define _sheader_h

static char *SRCID_sheader_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_IDENT
#ident	"@(#) sheader.h 1.4 11/23/92"
#endif /* HAVE_IDENT */

#include <stdio.h>
#include "typedefs.h"
#include "tpxio.h"

typedef struct		/* This struct describes the order and the	*/
  /* sizes of the structs in a statusfile, sizes are given in bytes.	*/
{
  int	ir_size;	/* Non zero if input_rec is present		*/
  int	e_size;		/* Non zero if energies are present		*/
  int	box_size;	/* Non zero if a box is present			*/
  int   vir_size;       /* Non zero if a virial is present              */
  int   pres_size;      /* Non zero if a pressure is present            */
  int	top_size;	/* Non zero if a topology is present		*/
  int	sym_size;	/* Non zero if a symbol table is present	*/
  int	x_size;		/* Non zero if coordinates are present		*/
  int	v_size;		/* Non zero if velocities are present		*/
  int	f_size;		/* Non zero if forces are present		*/

  int	natoms;		/* The total number of atoms			*/
  int	step;		/* Current step number				*/
  int	nre;		/* Nr of energies 				*/
  real	t;		/* Current time					*/
  real	lambda;		/* Current value of lambda			*/
} t_statheader;

extern void pr_header(FILE *fp,int indent,char *title,t_tpxheader *sh);
     /*
      * This routine prints out a (human) readable representation of
      * a header to the file fp. Ident specifies the number of spaces
      * the text should be indented. Title is used to print a header text.
      */

#endif	/* _sheader_h */
