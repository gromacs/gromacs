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

#define MAXSHIFT 6
 
typedef struct {
  int     il_code;              /* Code that determines the innerloop   */
                                /* corresponding to codes in nrnb.h     */
				/* Currently there are 19 different ones*/
  int     nri,maxnri;           /* Current/max number of i particles	*/
  int     nrj,maxnrj;		/* Current/max number of j particles	*/
  int     *iinr;		/* The i-elements			*/
  int     *gid;                 /* Index in energy arrays               */
  int     *shift;               /* Shift vector index (maxnri*MAXSHIFT) */
  int     *jindex;              /* Index in jjnr (maxnri*MAXSHIFT)+1    */
  int     *jjnr;		/* The j-atom list      	        */
} t_nblist;

/* For atom I =  nblist->iinr[N] (0 <= N < nblist->nri) we define
 * nblist->sindex[N+1] - nblist->sindex[N] different shift vector indices 
 * SI, i.e. (nblist->sindex[N] <= SI < nblist->sindex[N+1])
 * corresponding to atom I.
 * For shift vector S = nblist->shift[SI], the corresponding J atoms are
 * nblist->jjnr[JI]
 * with nblist->jindex[SI] <= JI < nblist->jindex[SI+1]
 *
 * Clear?
 */










