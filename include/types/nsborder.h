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

typedef struct {
  int  nodeid;			/* Node id	         		*/
  int  nnodes;			/* The number of nodes    		*/
  int  cgtotal; 		/* Total number of charge groups	*/
  int  natoms;			/* Total number of atoms		*/
  int  nstDlb;                  /* Every how many steps must we do load */
                                /* balancing                            */
  int  shift,bshift;		/* Coordinates are shifted left for     */
                                /* 'shift' systolic pulses, and right   */
				/* for 'bshift' pulses. Forces are      */
				/* shifted right for 'shift' pulses     */
				/* and left for 'bshift' pulses         */
				/* This way is not necessary to shift   */
				/* the coordinates over the entire ring */
  int  homenr[MAXNODES];       	/* The number of home particles		*/
  int  index[MAXNODES];		/* The starting of the home atoms	*/
  int  cgload[MAXNODES];        /* Division of charge groups over CPUS  */
                                /* This is static, i.e. it does not     */
				/* change during the simulation         */
  int  workload[MAXNODES];      /* This is the load for neighbor-       */
                                /* searching, this is initially the same*/
				/* as cgload, but may change due to     */
				/* dynamic load balancing               */
} t_nsborder;

#define START(nsb)  ((nsb)->index[(nsb)->nodeid])
#define HOMENR(nsb) ((nsb)->homenr[(nsb)->nodeid])
#define CG0(nsb)    (((nsb)->nodeid == 0) ? 0 : (nsb)->cgload[(nsb)->nodeid-1])
#define CG1(nsb)    ((nsb)->cgload[(nsb)->nodeid])
