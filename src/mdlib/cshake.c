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
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_cshake_c = "$Id$";

#include <math.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "pbc.h"
#include "txtdump.h"
#include "constr.h"
#include "vec.h"
#include "nrnb.h"

void cshake(atom_id iatom[],int ncon,int *nnit,int maxnit,
	    real dist2[],real xp[],real rij[],real m2[],
	    real invmass[],real tt[],real lagr[],int *nerror)
{
  const   real mytol=1e-6;
  
  int     ll,i,j,i3,j3,l3;
  int     ix,iy,iz,jx,jy,jz;
  real    toler,rpij2,rrpr,tx,ty,tz,diff,acor,im,jm;
  real    xh,yh,zh,rijx,rijy,rijz;
  real    tix,tiy,tiz;
  real    tjx,tjy,tjz;
  int     nit,error,iconv,nconv;
  
  error=0;
  nconv=1;
  for (nit=0; (nit<maxnit) && (nconv != 0) && (error == 0); nit++) {
    nconv=0;
    for(ll=0; (ll<ncon) && (error == 0); ll++) {
      l3    = 3*ll;
      rijx  = rij[l3+XX];
      rijy  = rij[l3+YY];
      rijz  = rij[l3+ZZ];
      i     = iatom[l3+1];
      j     = iatom[l3+2];
      i3    = 3*i;
      j3    = 3*j;
      ix    = i3+XX;
      iy    = i3+YY;
      iz    = i3+ZZ;
      jx    = j3+XX;
      jy    = j3+YY;
      jz    = j3+ZZ;
      
      tx      = xp[ix]-xp[jx];
      ty      = xp[iy]-xp[jy];
      tz      = xp[iz]-xp[jz];
      rpij2   = tx*tx+ty*ty+tz*tz;
      toler   = dist2[ll];
      diff    = toler-rpij2;
      
      /* iconv is zero when the error is smaller than a bound */
      iconv   = fabs(diff)*tt[ll];
      
      if (iconv != 0) {
	nconv   = nconv + iconv;
	rrpr    = rijx*tx+rijy*ty+rijz*tz;
	
	if (rrpr < toler*mytol) 
	  error=ll;
	else {
	  acor      = diff*m2[ll]/rrpr;
	  lagr[ll] += acor;
	  xh        = rijx*acor;
	  yh        = rijy*acor;
	  zh        = rijz*acor;
	  im        = invmass[i];
	  jm        = invmass[j];
	  if ((im != 0) && (jm != 0)) {
	    xp[ix] += xh*im;
	    xp[iy] += yh*im;
	    xp[iz] += zh*im;
	    xp[jx] -= xh*jm;
	    xp[jy] -= yh*jm;
	    xp[jz] -= zh*jm;
	  }
	  else if (im == 0) {
	    xp[ix] += xh*jm;
	    xp[iy] += yh*jm;
	    xp[iz] += zh*jm;
	  }
	  else if (jm == 0) {
	    xp[jx] -= xh*im;
	    xp[jy] -= yh*im;
	    xp[jz] -= zh*im;
	  }
	  else 
	    fatal_error(0,"Constraint between two massless particles %d and %",
			im,jm);
	}
      }
    }
  }
  *nnit=nit;
  *nerror=error;
}

