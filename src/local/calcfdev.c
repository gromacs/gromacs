#include <math.h>
#include "typedefs.h"
#include "main.h"
#include "vec.h"
#include "txtdump.h"

void calc_force(int natom,rvec f[],rvec fff[])
{
  int  i,j,m;
  int  jindex[] = { 0, 5, 10};
  rvec dx,df;
  real msf1,msf2;
  
  for(j=0; (j<2); j++) {
    clear_rvec(fff[j]);
    for(i=jindex[j]; (i<jindex[j+1]); i++) {
      for(m=0; (m<DIM); m++) {
	fff[j][m] += f[i][m];
	xcm[j][m] += 
      }
    }
  }
  
  msf1 = iprod(fff[0],fff[0]);
  msf2 = iprod(fff[1],fff[1]);
  if (debug) {
    pr_rvecs(debug,0,"force",f,natom);
    
    fprintf(debug,"FMOL:  %10.3f  %10.3f  %10.3f  %10.3f  %10.3f  %10.3f\n",
	    fff[0][XX],fff[0][YY],fff[0][ZZ],fff[1][XX],fff[1][YY],fff[1][ZZ]);
    fprintf(debug,"RMSF:  %10.3e  %10.3e\n",msf1,msf2);
  }
}


void calc_f_dev(int natoms,real charge[],rvec x[],rvec f[],
		t_idef *idef,real *xiH,real *xiS)
{
  enum { wwwO, wwwH1, wwwH2, wwwS, wwwNR };
  int  lj_index[wwwNR] = { 0,    4,     4,     8 };
  real rmsf,dFH,dFS,dFO,dr_14;
  real q[wwwNR],c6[wwwNR],c12[wwwNR],c12ratio;
  rvec fff[2],dx;
  int  i,j,aj;
  
  for(i=0; (i<wwwNR); i++) {
    q[i]   = charge[i];
    c12[i] = idef->iparams[lj_index[i]].lj.c12;
    c6[i]  = idef->iparams[lj_index[i]].lj.c6;
  }
  
  calc_force(natoms,f,fff);
  rmsf = norm(fff[0]);
  dFS = 0;
  dFH = 0;
  dFO = 0;
  for(i=0; (i<4); i++) {
    for(j=4; (j<8); j++) {
      if (c12[i] != 0) {
	rvec_sub(x[i],x[j],dx);
	aj       = j % 4;
	dr_14    = pow(iprod(dx,dx),-7);
	c12ratio = -12*sqrt(c12[aj]/c12[i])*dr_14*iprod(fff[0],dx);
	
	switch (i) {
	case wwwH1:
	case wwwH2:
	  dFH += c12ratio;
	  break;
	case wwwS:
	  dFS += c12ratio;
	  break;
	case wwwO:
	  dFS += c12ratio;
	  break;
	}
      }
    }
  }

  if (debug) {    
    fprintf(debug,"FFF: dFS=%10.3e,  dFH=%10.3e,  dFO=%10.3e, rmsf=%10.3e\n",
	    dFS,dFH,dFO,rmsf);
  }
  if (dFH == 0)
    *xiH = 1;
  else
    *xiH=rmsf/(10*c12[wwwH1]*dFH);
    
  if (dFS == 0)
    *xiS = 1;
  else
    *xiS=rmsf/(10*c12[wwwS]*dFS);
}
