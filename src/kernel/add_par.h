/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * GRoups of Organic Molecules in ACtion for Science
 */
#ifndef _add_par_h
#define _add_par_h

static char *SRCID_add_par_h = "$Id$";

#include "typedefs.h"

extern void add_param(t_params *ps,int ai,int aj, real *c);

extern void add_imp_param(t_params *ps,int ai,int aj,int ak,int al,
			  real c0, real c1);

extern void add_dum2_2_param(t_params *ps,int ai,int aj,int ak,real a);

extern void add_dum3_2_param(t_params *ps,int ai,int aj,int ak,int al,
			     real a, real b);

extern void add_dum3_3_param(t_params *ps,int ai,int aj,int ak,int al,
			     real a, real b, real c);

extern int search_jtype(t_restp *rp,char *name,bool bFirstRes);

#endif	/* _add_par_h */
