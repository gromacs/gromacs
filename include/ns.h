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

#ifndef _ns_h
#define _ns_h

static char *SRCID_ns_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_IDENT
#ident	"@(#) ns.h 1.50 2/2/97"
#endif /* HAVE_IDENT */

#include <stdio.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "pbc.h"
#include "tgroup.h"
#include "nsb.h"
#include "network.h"

/****************************************************
 *
 *    U T I L I T I E S May be found in ns.c
 *
 ****************************************************/
extern void init_neighbor_list(FILE *log,t_forcerec *fr,int homenr);
/* 
 * nn is the number of energy terms in the energy matrix
 * (ngener*(ngener-1))/2
 * start is the first atom on this processor
 * homenr is the number of atoms on this processor
 */
 
extern int calc_naaj(int icg,int cgtot);
/* Calculate the number of charge groups to interact with for icg */

/****************************************************
 *
 *    N E I G H B O R  S E A R C H I N G
 *
 *    Calls either ns5_core (when grid selected in .mdp file)
 *    or ns_simple_core (when simple selected in .mdp file)
 *
 *    Return total number of pairs searched 
 *
 ****************************************************/
extern int search_neighbours(FILE *log,t_forcerec *fr,
			     rvec x[],matrix box,
			     t_topology *top,t_groups *grps,
			     t_commrec *cr,t_nsborder *nsb,t_nrnb *nrnb,
			     t_mdatoms *md,real lambda,real *dvdlambda);
 

#endif	/* _ns_h */
