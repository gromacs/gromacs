/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * S  C  A  M  O  R  G
 */

#ifndef	_tconf_h
#define	_tconf_h

#ifdef HAVE_IDENT
#ident	"@(#) tconf.h 1.6 2/2/97"
#endif /* HAVE_IDENT */

#include <stdio.h>
#include "typedefs.h"
#include "statusio.h"

typedef struct
{
  t_statheader header;
  t_inputrec *ir;
  matrix *box;
  tensor *vir,*pres;
  rvec *x;
  rvec *v;
  rvec *f;
  t_energy *e;
  t_topology *top;

  int step;
  real t;
  real lambda;
  int natoms;
  int nre;
} t_status_parm;

typedef struct
{
  t_status_parm data;	/* Passed to read_status for reading in data      */
  t_status_parm config; /* Complete allocated configuration (first step?) */
} t_config;

extern void set_configuration(t_config *config);

extern char *rd_configuration(FILE *fp,t_config *cf);

extern void wr_configuration(FILE *fp,t_config *cf);

extern char *init_configuration(FILE *fp,t_config *cf);

extern char *read_configuration(char *fn,t_config *cf);

extern void pr_configuration(FILE *fp,int indent,char *title,
                             t_status_parm *sp);
                             
#endif	/* _tconf_h */
