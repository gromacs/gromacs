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
 * Gyas ROwers Mature At Cryogenic Speed
 */
#ifdef SHORT_NL
typedef unsigned short t_nl_j;
#else
typedef int t_nl_j;
#endif

typedef struct {
  t_ishift shift;		/* The shift index			*/
  atom_nr  nj;			/* The number of j-particles 		*/
  int      j_index;             /* Starting index in j_array            */
  atom_id  i_atom;		/* The i-atom 				*/
  bool     bWater;              /* TRUE if this is a water molecule     */
} t_nl_i;

typedef struct {
  atom_nr nri;			/* The number of i particles		*/
  atom_nr nrj;			/* The number of j particles		*/
  atom_nr maxnri;		/* Max number of i particles		*/
  atom_nr maxnrj;               /* Max number of j particles            */
  t_nl_i  *nl_i;		/* The i-elements			*/
  t_nl_j  *nl_j;		/* The j-atom list (in ints or shorts)	*/
} t_nblist;

/* Structures for buffering in neighbour searching */
#define MAXNB_LR 1024

typedef struct {
  int     nj;
  int     nlj[2*MAXNB_LR];
} t_nblist_lr;












