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
  int     il_code;             /* Code that determines the innerloop    */
                               /* corresponding to codes in nrnb.h      */
  int     nri,maxnri;          /* Current/max number of i particles	*/
  int     nrj,maxnrj;	       /* Current/max number of j particles	*/
  int     maxlen;              /* maxnr of j atoms for a single i atom 	*/
  int     solvent;             /* type of solvent optimization          */
  int     *iinr;	       /* The i-elements			*/
  int     *gid;                /* Index in energy arrays                */
  int     *shift;              /* Shift vector index                    */
  int     *jindex;             /* Index in jjnr                         */
  int     *jjnr;	       /* The j-atom list                       */
  int     *nsatoms;            /* list with number of atoms for general */
                               /* solvents. There are two entries for 	*/
                               /* each molecule - first is total natoms */
                               /* and second how many at the beginning 	*/
                               /* have LJ interactions.                 */
                               /* This is NOT used for water!           */
#ifdef USE_THREADS
  int      count;              /* counter to multithread the innerloops */
  pthread_mutex_t *mtx;        /* mutex to lock the counter             */
#else
  int      pad1,*pad2;         /* padding to make size constant         */
#endif
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










