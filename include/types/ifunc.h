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

#ifndef _ifunc_h
#define _ifunc_h


typedef real t_ifunc(int nbonds,t_iatom iatoms[],t_iparams *iparams,
                     rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
		     matrix box,real lambd,real *dvdlambda,
		     t_mdatoms *md,int ngrp,real egnb[],real egcoul[]);
/*
 * The function type t_ifunc() calculates one interaction, using iatoms[] 
 * and iparams. Within the function the number of atoms to be used is 
 * known. Within the function only the atomid part of the iatoms[] array 
 * is supplied, not the type field (see also t_ilist). The function 
 * returns the potential energy. The coordinates in x are such that
 * no calculation of PBC is necessary.
 */

#define IF_NULL       0
#define IF_BOND       1
#define IF_DUMMY      1<<1
#define IF_CONSTRAINT 1<<2
#define IF_CONNECT    1<<3
#define IF_BTYPE      1<<5
#define IF_ATYPE      1<<6
/* These flags tell to some of the routines what can be done with this
 * item in the list. If flags & IF_BOND, then bonded interactions will
 * be calculated. If flags & IF_CONNECT this link specifies a connection 
 * (chemical bond) between two particles. By specifying this here, we can 
 * keep all the information in one place.
 */
typedef struct
{
  char    *name;	/* the name of this function			*/
  char    *longname;    /* The name for printing etc.                   */
  int     nratoms;	/* nr of atoms needed for this function		*/
  int     nrfpA,nrfpB;  /* number of parameters for this function.      */
                        /* this corresponds to the number of params in  */
                        /* iparams struct! (see idef.h)                 */
  /* A and B are for normal and free energy components respectively.    */
  unsigned long   flags;        /* Flags (see above)                            */
  int     nrnb_ind;     /* index for nrnb (-1 if unknown)               */
  t_ifunc *ifunc;	/* the function it self				*/
} t_interaction_function;

#define NRFP(ftype) (interaction_function[(ftype)].nrfpA+interaction_function[(ftype)].nrfpB)
#define NRAL(ftype) (interaction_function[(ftype)].nratoms)

extern t_interaction_function interaction_function[F_NRE];
/* initialised interaction functions descriptor				*/

#endif

