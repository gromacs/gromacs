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
 * GRowing Old MAkes el Chrono Sweat
 */

#ifndef _add_par_h
#define _add_par_h

static char *SRCID_add_par_h = "$Id$";

#include "typedefs.h"
#include "pdb2top.h"

extern void add_param(t_params *ps, int ai, int aj, real *c, char *s);

extern void add_imp_param(t_params *ps, int ai, int aj, int ak, int al,
			  real c0, real c1, char *s);
			  
extern void add_dih_param(t_params *ps,int ai,int aj,int ak,int al,
			  real c0, real c1, real c2, char *s);

extern void add_dum2_atoms(t_params *ps, int ai, int aj, int ak);

extern void add_dum3_atoms(t_params *ps, int ai, int aj, int ak, int al, 
			   bool bSwapParity);

extern void add_dum3_param(t_params *ps, int ai, int aj, int ak, int al, 
			   real c0, real c1);

extern void add_dum4_atoms(t_params *ps, int ai, int aj, int ak, int al, 
			   int am);

extern int search_jtype(t_restp *rp,char *name,bool bFirstRes);

extern void cp_param(t_param *dest,t_param *src);

#endif	/* _add_par_h */
