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
static char *SRCID_topshake_c = "$Id$";

#include <ctype.h>

#include "sysstuff.h"
#include "physics.h"
#include "macros.h"
#include "typedefs.h"
#include "topshake.h"
#include "toppush.h"
#include "toputil.h"
#include "topdirs.h"

static void copy_bond (t_params *pr, int to, int from)
/* copies an entry in a bond list to another position.
 * does no allocing or freeing of memory
 */
{
  memcpy((char*) &(pr->param[to]),(char*) &(pr->param[from]),
	 (size_t)sizeof(pr->param[from]));
}

static int count_hydrogens (char ***atomname, int nra, atom_id a[])
{
  int  i,nh;
  
  if (!atomname) 
    fatal_error(0,"Error in any_hydrogen");
    
  nh=0;
  for (i=0; (i<nra); i++) 
    if (toupper(**(atomname[a[i]]))=='H')
      nh++;
 
  return nh;
}

void make_shake (t_params plist[],t_atoms *atoms,t_atomtype *at,int nshake)
{
  char         ***info=atoms->atomname;
  t_params     *pr,*shake;
  t_params     *bonds;
  t_param      p,*bond,*ang;
  real         b_ij,b_jk;
  int          nb,b,i,j;
  bool         bFound;

  switch (nshake) {
  case eshNONE:
    return;
  case eshHBONDS:
    printf("turning H bonds into constraints...\n");
    break;
  case eshALLBONDS:
    printf("turning all bonds into constraints...\n");
    break;
  case eshHANGLES:
    printf("turning all bonds and H angles into constraints...\n");
    break;
  case eshALLANGLES:
    printf("turning all bonds and angles into constraints...\n");
    break;
  default:
    fatal_error(0,"Invalid option for make_shake (%d)",nshake);
  }
  
  /* Add all the angles with hydrogens to the shake list
   * and remove them from the bond list
   */
  bonds=&(plist[F_BONDS]);
  shake=&(plist[F_SHAKE]);
  
  if ((nshake == eshHANGLES) || (nshake == eshALLANGLES)) {
    /* horrible shortcut */
    pr = &(plist[F_ANGLES]);
    for (i=0; (i < pr->nr); ) {
      ang=&(pr->param[i]);
#ifdef DEBUG
      printf("Angle: %d-%d-%d\n",ang->AI,ang->AJ,ang->AK); fflush(stdout);
#endif
      if ((nshake == eshALLANGLES) || 
	  (count_hydrogens(info,3,ang->a) > 0)) {
	/* Can only add hydrogen angle shake, if the two bonds
	 * are constrained.
	 * append this angle to the shake list 
	 */
	p.AI = ang->AI;
	p.AJ = ang->AK;
	
	/* Calculate length of constraint */
	bFound=FALSE;
	b_ij=b_jk=0.0;
	for (j=0; !bFound && (j<bonds->nr); j++) {
	  bond=&(bonds->param[j]);
	  if (((bond->AI==ang->AI) && 
	       (bond->AJ==ang->AJ)) ||
	      ((bond->AI==ang->AJ) && 
	       (bond->AJ==ang->AI)))
	    b_ij=bond->C0;
	  if (((bond->AI==ang->AK) && 
	       (bond->AJ==ang->AJ)) ||
	      ((bond->AI==ang->AJ) && 
	       (bond->AJ==ang->AK)))
	    b_jk=bond->C0;
	  bFound = (b_ij!=0.0) && (b_jk!=0.0);
	}
	/* apply law of cosines */
	if (!bFound) {
	  fprintf(stderr,"No bond information for bond %s-%s or %s-%s\n",
		  *info[ang->AI],*info[ang->AJ],
		  *info[ang->AJ],*info[ang->AK]);
	  exit(1);
	}
	p.C0 = sqrt(b_ij*b_ij+b_jk*b_jk-2.0*b_ij*b_jk*cos(DEG2RAD*ang->C0));
	p.C1 = p.C0;
#ifdef DEBUG
	printf("p: %d, q: %d, dist: %12.5e\n",p.AI,p.AJ,p.C0);
#endif
	push_bondnow (shake,&p);
	/* move the last bond to this position */
	copy_bond (pr,i,pr->nr-1);
	/* should free memory here!! */
	pr->nr--;
      }
      else
	i++;
    }
  }
  if ((nshake == eshALLANGLES) || (nshake == eshHANGLES))
    nshake = eshALLBONDS;
  
  /* Add all the bonds with hydrogens to the shake list
   * and remove them from the bond list
   */
  if ((nshake == eshHBONDS) || (nshake == eshALLBONDS)) {
    /* horrible shortcut */
    pr = &(plist[F_BONDS]);
    for (i=0; (i < pr->nr); ) {
      if ((nshake == eshALLBONDS) || 
	  (count_hydrogens (info,2,pr->param[i].a) > 0)) {
	/* append this bond to the shake list */
	p.AI = pr->param[i].AI;
	p.AJ = pr->param[i].AJ;
	p.C0 = pr->param[i].C0;
	p.C1 = pr->param[i].C2;
	push_bondnow (shake,&p);
	
	/* move the last bond to this position */
	copy_bond (pr,i,pr->nr-1);
	
	/* should free memory here!! */
	pr->nr--;
      }
      else
	i++;
    }
  }
}
