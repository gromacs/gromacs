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

#ifndef _sheader_h
#define _sheader_h

static char *SRCID_sheader_h = "$Id$";

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
