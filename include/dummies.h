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
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _dummies_h
#define _dummies_h

static char *SRCID_dummies_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "typedefs.h"

typedef struct {
  int nprevdum; /* how many dummy particles are nonlocal */     
  int nnextdum;
  int *idxprevdum; /* index of nonlocal dummy particles */
  int *idxnextdum;
  int nprevconstr; /* how many constr. atoms are nonlocal */
  int nnextconstr;
  int *idxprevconstr; /* indices of nonlocal constructing atoms */
  int *idxnextconstr;
} t_comm_dummies;

/* Communication routines for dummies. The coordinates and
 * forces are only move on a need-to-know basis, usually only
 * 2-3 atoms per processor. To achieve this small amount of
 * communication, and to limit it to nearest neighbour messages,
 * we demand that dummies are not spread over nonadjacent nodes.
 * Thus, keep your dummies close to your constructing atoms.
 * (mdrun & grompp will report an error otherwise)
 */ 

extern void move_construct_x(t_comm_dummies *dummycomm,rvec x[],t_commrec *cr); 
/* Move coords of nonlocal constructing atoms */

extern void move_dummy_xv(t_comm_dummies *dummycomm,rvec x[],rvec v[],t_commrec *cr);
/* Send the coordinates and velocity of a constructed dummy to the home node */

extern void move_dummy_f(t_comm_dummies *dummycomm,rvec f[],t_commrec *cr);
/* Get dummy forces from the home node */

extern void move_construct_f(t_comm_dummies *dummycomm,rvec f[],t_commrec *cr);
/* Send spreaded forces to nonlocal constructing atoms */


extern void construct_dummies(FILE *log,rvec x[],t_nrnb *nrnb,
			      real dt,rvec v[],t_idef *idef);
/* Create positions of dummy atoms based on surrounding atoms.
 */
 
extern void spread_dummy_f(FILE *log,rvec x[],rvec f[],rvec buf[],
			   t_nrnb *nrnb,t_idef *idef);
/* Spread the force operating on the dummy atoms on the surrounding atoms.
 */

#endif

