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
 * Grunge ROck MAChoS
 */
enum { egCOUL, egLJ, egBHAM, egLR, egCOUL14, egLJ14, egNR };
	
typedef struct {
  real    T;		/* Temperature	    */
  real    lambda;       /* T Scaling factor */
  tensor  ekin;		/* Kinetic energy   */
} t_grp_tcstat;

typedef struct {
  int           nn;             /* Number of terms 			*/
  real 		*ee[egNR];	/* Arrays of energy terms for THIS 	*/
  				/* group with ALL other groups		*/
} t_grp_ener;

typedef struct {
  int     nat;		/* Number of atoms in this group		*/
  atom_id *aid;		/* Atom ids of the atoms in this group		*/
  real    M;		/* Total mass of group				*/
  rvec	  u;           	/* Mean velocities of home particles    	*/
  rvec	  uold;        	/* Old velocities of home particles    	        */
  rvec	  ut;          	/* Mean velocities of home particles    	*/
} t_grp_acc;

typedef struct {
  t_grp_ener   estat;		/* Energy logging stuff			*/
  t_grp_tcstat *tcstat;         /* T Coupling Data 			*/
  t_grp_acc    *grpstat;	/* Acceleration data			*/
} t_groups;

#define GID(igid,jgid,gnr) ((igid < jgid) ? (igid*gnr+jgid) : (jgid*gnr+igid))

