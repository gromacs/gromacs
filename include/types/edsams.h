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
  bool 		bEdsam;		/* Do ED sampling?			*/
  char          *edinam; 	/* name of ED sampling input file       */
  char          *edonam;        /*                     output           */
} t_edsamyn;

typedef struct {
  int 		neig;		/* nr of eigenvectors            	*/
  int           *ieig;          /* index nrs of eigenvectors            */
  real          *stpsz;         /* stepsizes (per eigenvector)          */
  rvec          **vec;          /* eigenvector components               */
  real          *xproj;         /* instantaneous x projections          */
  real          *vproj;         /* instantaneous v projections          */
  real          *fproj;         /* instantaneous f projections          */
  real          *refproj;       /* starting or target projecions        */
  real          radius;         /* instantaneous radius                 */
} t_eigvec;

typedef struct {
  t_eigvec 	mon;		/* only monitored, no constraints       */
  t_eigvec 	linfix;		/* fixed linear constraints             */
  t_eigvec 	linacc;		/* acceptance linear constraints        */
  t_eigvec 	radfix;		/* fixed radial constraints (exp)       */
  t_eigvec 	radacc;		/* acceptance radial constraints (exp)  */  
  t_eigvec 	radcon;		/* acceptance rad. contraction constr.  */
} t_edvecs;

typedef struct {
  int 		nr;		/* Nr. of atoms        			*/
  int           *anrs;          /* Index numbers                        */
  rvec          *x;             /* Positions                            */
  matrix        box;            /* Box lenghts                          */
} t_edx;


typedef struct {
  int 		nini;		/* Total Nr of atoms    		*/
  int           npro;           /* Nr of protein atoms                  */
  int           ned;            /* Nr of atoms in essdyn                */
  bool          selmas;         /* true if trans fit with cm            */
  int           outfrq;         /* freq (in steps) of writing output    */
  int           logfrq;         /* freq (in steps) of writing to log    */
  int           maxedsteps;     /* max nr of steps per cycle            */
  t_edx         sref;           /* reference positions                  */
  t_edx         sav;            /* average positions                    */
  t_edvecs      vecs;           /* eigenvectors                         */
  real          slope;          /* minimal slope in acceptance radexp   */
  t_edx         star;           /* target positions                     */
  t_edx         sori;           /* origin positions                     */
  int           nmass;          /* Nr of masses                         */
  int           *masnrs;        /* index nrs for atoms with masses      */
  real          *mass;          /* atomic masses                        */
  real          tmass;          /* total mass                           */
  int           nfit;           /* Number of atoms to use for rot fit   */
  int           *fitnrs;        /* index nrs of atoms to use for rot fit*/
  FILE          *edo;           /* output file                          */
} t_edpar;

