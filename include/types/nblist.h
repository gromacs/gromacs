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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

typedef struct {
  int     il_code;             /* Code that determines the innerloop    */
                               /* corresponding to codes in nrnb.h      */
  int     nri,maxnri;          /* Current/max number of i particles	*/
  int     nrj,maxnrj;	       /* Current/max number of j particles	*/
  int     maxlen;              /* maxnr of j atoms for a single i atom 	*/
  int     solvent;             /* type of solvent optimization          */
  atom_id *iinr;	       /* The i-elements			*/
  int     *gid;                /* Index in energy arrays                */
  int     *shift;              /* Shift vector index                    */
  int     *jindex;             /* Index in jjnr                         */
  atom_id *jjnr;	       /* The j-atom list                       */
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

/* For atom I =  nblist->iinr[N] (0 <= N < nblist->nri) there can be
 * several neighborlists (N's), for different energy groups (gid) and
 * different shifts (shift).
 * For corresponding J atoms for each list are are:
 * nblist->jjnr[JI]
 * with nblist->jindex[N] <= JI < nblist->jindex[N+1]
 *
 * Clear?
 */










