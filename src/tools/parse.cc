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
 * GROtesk MACabre and Sinister
 */

#ifndef _parse_cc
#define _parse_cc

static char *SRCID_parse_cc = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) parse.cc 1.20 9/30/97"
#endif /* HAVE_IDENT */

#include "parse.h"
#include "typedefs.h"
#include "statutil.h"
#include "physics.h"
#include "sysstuff.h"
#include "pbc.h"

/*
 *for parse_common_args the following variables must be filled
 * -f status.trj     ; the trajectory
 * -p status.tpb     ; the topology
 * -n stat.ndx       ; the selectionable groups
 * -b begin time     ;
 * -e end time       ;
 * -r cut off radius ;
 * -a cut off angle  ;
 * -o output file    ;
 * -v verbosive      ; print a new hydrogen bond to stdout when it is found
 */
 
void init_topology(char *topology)
{
  t_statheader header;
  int          step,nre;
  real         t,lambda;
  t_inputrec   ir;
  int          natoms;

  read_status_header(topology,&header);
  x = (rvec *)malloc(header.natoms*sizeof(rvec));

  read_status(topology,
              &step,&t,&lambda,&ir,
              box,NULL,NULL,
              &natoms,
              x,NULL,NULL,&nre,NULL,
              top);

  init_pbc(box,FALSE);
}

void user_input(int *nr_groups)
{
  int a;
  fprintf(stderr,"1. Analyse Hydrogen Bonds of 1 group\n");
  fprintf(stderr,"2. Analyse Hydrogen Bond interaction of 2 (or more) groups\n");
  fprintf(stderr,"3. Analyse Selected Hydrogen bonds (Formatted input)\n");
  fprintf(stderr,"4. Analyse Water insertion (Formatted input) \n");
  fprintf(stderr,"5. Count Hydrogen bonds \n");
  fprintf(stderr,"\nSelect a number: ");
  do {
    scanf("%d",&a);
  } while ((a<1)||(a>5));
  switch (a) {
  case 1:
    mode=INTRA_GRP;
    break;
  case 2:
    mode=INTER_GRP;
    fprintf(stderr,"\nHow many groups: ");
    scanf("%d",nr_groups);
    if ((*nr_groups<2)||(*nr_groups>100)) {
      fprintf(stderr,"Illegal choice\n");
      exit(0);
    }
    break;
  case 3:
    mode=SELECTED;
    break;
  case 4:
    mode=INSERT;
    *nr_groups=2;
    break;
  case 5:
    mode=COUNT;
    *nr_groups=2;
    break;
  default:
    fprintf(stderr,"Ilegal choice, program terminated");
    exit(0);
    break;
  }
}

#endif	/* _parse_cc */





