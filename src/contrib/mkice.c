/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.99_development_20071104
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2006, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <math.h>
#include "typedefs.h"
#include "gromacs/commandline/pargs.h"
#include "copyrite.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/fileio/pdbio.h"
#include "macros.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/math/units.h"
#include "names.h"
#include "txtdump.h"
#include "gromacs/fileio/trnio.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/fileio/strdb.h"
#include "gromacs/fileio/confio.h"

#define TET   109.47
#define DCONS 0.117265878

typedef struct {
  int n,aa[4];
} t_bbb;

static char *watname[] = { "OW ", "HW1", "HW2", "DW", "SW" };
static char *diamname[] = { "C1", "H2", "H3", "H4", "H5", "H2", "H3", "H4", "H5" };
static real qspc[]     = { -0.82, 0.41, 0.41 };
static real qyaw[]     = { 1.24588, 0.62134, 0.62134, 0.0, -2.48856 };
static real spc_lj[3][6] = {
  { 2.6171e-3, 2.6331e-6, 0, 0, 0, 0 },
  { 0,      0,      0, 0, 0, 0 },
  { 0,      0,      0, 0, 0, 0 }
};
#define CHH 2e-9
#define CHS 2e-9
#define COS 2e-6
static real yaw_lj[5][10] = {
  { 0, 0,   0, 0,   0,   0, 0, 0, 0,      COS },
  { 0, 0,   0, CHH, 0, CHH, 0, 0, 0,      CHS },
  { 0, 0,   0, CHH, 0, CHH, 0, 0, 0,      CHS },
  { 0, 0,   0, 0,   0,   0, 0, 0, 0,      0   },
  { 0, COS, 0, CHS, 0, CHS, 0, 0, 2.6e-3, 0   }
};

void unitcell(rvec x[],rvec box,gmx_bool bYaw,real odist,real hdist)
{
#define cx  0.81649658
#define cy  0.47140452
#define cy2 0.94280904
#define cz  0.33333333
  
  rvec xx[24] = {
    { 0,   0,         0 }, /* O1 */
    { 0,   0,         1 }, /* H relative to Oxygen */
    { cx, cy,       -cz },
    { cx, cy,       -cz }, /* O2 */
    { 0, 0,       -1    }, /* H relative to Oxygen */
    { cx,-cy,       +cz },
    { cx, cy+cy2,     0 }, /* O3 */
    { -cx, cy,    -cz   }, /* H relative to Oxygen */
    { 0,   -cy2,    -cz },
    { 0,  2*cy+cy2, -cz }, /* O4 */
    {-cx,-cy,       +cz }, /* H relative to Oxygen */
    { 0 , cy2,      +cz },
    { 0,   0,         1 }, /* O5 */
    {-cx, cy,       +cz }, /* H relative to Oxygen */
    { 0 , -cy2,     +cz },
    { cx, cy,      1+cz }, /* O6 */
    { -cx, -cy,     -cz }, /* H relative to Oxygen */
    { 0,   cy2,     -cz },
    { cx, cy+cy2,     1 }, /* O7 */
    { 0,  0,       -1   }, /* H relative to Oxygen */
    { cx, cy,       +cz },
    { 0,  2*cy+cy2,1+cz }, /* O8 */
    { 0,   0,         1 }, /* H relative to Oxygen */
    { cx,   -cy,    -cz }
  };
  int  i,iin,iout,j,m;
  rvec tmp,t2,dip;
  
  clear_rvec(dip);
  for(i=0; (i<8); i++) {
    iin = 3*i;
    if (bYaw)
      iout = 5*i;
    else
      iout = iin;
    svmul(odist,xx[iin],x[iout]);
    svmul(-0.82,x[iout],t2);
    rvec_inc(dip,t2);
    for(j=1; (j<=2); j++) {
      svmul(hdist,xx[iin+j],tmp);
      rvec_add(x[iout],tmp,x[iout+j]);
      svmul(0.41,x[iout+j],t2);
      rvec_inc(dip,t2);
    }
    if (bYaw)
      for(m=0; (m<DIM); m++) 
	x[iout+3][m] = x[iout+4][m] = 
	  (1-2*DCONS)*x[iout][m]+DCONS*(x[iout+1][m]+x[iout+2][m]);
  }
  
  box[XX] = 2*cx;
  box[YY] = 2*(cy2+cy);
  box[ZZ] = 2*(1+cz);
  for(i=0; (i<DIM); i++)
    box[i] *= odist;
    
  printf("Unitcell:  %10.5f  %10.5f  %10.5f\n",box[XX],box[YY],box[ZZ]);
  printf("Dipole:    %10.5f  %10.5f  %10.5f (e nm)\n",dip[XX],dip[YY],dip[ZZ]);
}

void random_h_coords(int natmol,int nmol,rvec x[],rvec box,
		     gmx_bool bYaw,real odist,real hdist)
{
#define cx  0.81649658
#define cy  0.47140452
#define cy2 0.94280904
#define cz  0.33333333
  
  rvec xx[24] = {
    { 0,   0,         0 }, /* O1 */
    { 0,   0,         1 }, /* H relative to Oxygen */
    { cx, cy,       -cz },
    { cx, cy,       -cz }, /* O2 */
    { 0, 0,       -1    }, /* H relative to Oxygen */
    { cx,-cy,       +cz },
    { cx, cy+cy2,     0 }, /* O3 */
    { -cx, cy,    -cz   }, /* H relative to Oxygen */
    { 0,   -cy2,    -cz },
    { 0,  2*cy+cy2, -cz }, /* O4 */
    {-cx,-cy,       +cz }, /* H relative to Oxygen */
    { 0 , cy2,      +cz },
    { 0,   0,         1 }, /* O5 */
    {-cx, cy,       +cz }, /* H relative to Oxygen */
    { 0 , -cy2,     +cz },
    { cx, cy,      1+cz }, /* O6 */
    { -cx, -cy,     -cz }, /* H relative to Oxygen */
    { 0,   cy2,     -cz },
    { cx, cy+cy2,     1 }, /* O7 */
    { 0,  0,       -1   }, /* H relative to Oxygen */
    { cx, cy,       +cz },
    { 0,  2*cy+cy2,1+cz }, /* O8 */
    { 0,   0,         1 }, /* H relative to Oxygen */
    { cx,   -cy,    -cz }
  };
  int  i,iin,iout,j,m;
  rvec tmp,t2,dip;
  
  clear_rvec(dip);
  for(i=0; (i<nmol); i++) {
    iin = natmol*i;
    iout = iin;
    svmul(odist,x[iin],x[iout]);
    svmul(-0.82,x[iout],t2);
    rvec_inc(dip,t2);
    for(j=1; (j<=2); j++) {
      svmul(hdist,xx[3*(i % 8)+j],tmp);
      rvec_add(x[iout],tmp,x[iout+j]);
      svmul(0.41,x[iout+j],t2);
      rvec_inc(dip,t2);
    }
  }
  
  box[XX] = 2*cx;
  box[YY] = 2*(cy2+cy);
  box[ZZ] = 2*(1+cz);
  for(i=0; (i<DIM); i++)
    box[i] *= odist;
    
  printf("Unitcell:  %10.5f  %10.5f  %10.5f\n",box[XX],box[YY],box[ZZ]);
  printf("Dipole:    %10.5f  %10.5f  %10.5f (e nm)\n",dip[XX],dip[YY],dip[ZZ]);
}

void unitcell_d(rvec x[],rvec box,real odist)
{
  rvec cc[8] = {
    { 0,   0,         0 }, /* C1 */
    { cx, cy,       -cz }, /* C2 */
    { cx, cy+cy2,     0 }, /* C3 */
    { 0,  2*cy+cy2, -cz }, /* C4 */
    { 0,   0,         1 }, /* C5 */
    { cx, cy,      1+cz }, /* C6 */
    { cx, cy+cy2,     1 }, /* C7 */
    { 0,  2*cy+cy2,1+cz }, /* C8 */
  };
  rvec hh[4] = {
    { 0,   0,         1  }, /* H relative to C */
    { cx,  cy,       -cz },
    { cx, -cy,       -cz }, 
    {-cy2,  0,       -cz }
  };
  int  i,iin,iout,j,m;
  rvec tmp,t2,dip;
  
  clear_rvec(dip);
  for(i=0; (i<8); i++) {
    iin  = i;
    iout = i;
    svmul(odist,cc[iin],x[iout]);
  }  
  box[XX] = 2*cx;
  box[YY] = 2*(cy2+cy);
  box[ZZ] = 2*(1+cz);
  for(i=0; (i<DIM); i++)
    box[i] *= odist;
    
  printf("Unitcell:  %10.5f  %10.5f  %10.5f\n",box[XX],box[YY],box[ZZ]);
}

static t_bbb *mk_bonds(int natoms,rvec x[],real odist,
		       gmx_bool bPBC,matrix box)
{
  real  od2 = odist*odist+1e-5;
  t_pbc pbc;
  t_bbb *bbb;
  int   i,j;
  rvec  dx;
  
  if (bPBC)
    set_pbc(&pbc,box);
  snew(bbb,natoms);
  for(i=0; (i<natoms); i++) {
    for(j=i+1; (j<natoms); j++) {
      if (bPBC)
	pbc_dx(&pbc,x[i],x[j],dx);
      else
	rvec_sub(x[i],x[j],dx);
      if (iprod(dx,dx) <= od2) {
	bbb[i].aa[bbb[i].n++] = j;
	bbb[j].aa[bbb[j].n++] = i;
      }
    }
  }
  if (debug) 
#define PRB(nn) (bbb[(i)].n >= (nn)) ? bbb[i].aa[nn-1] : -1
    for(i=0; (i<natoms); i++)
      fprintf(debug,"bonds from %3d:  %d %d %d %d\n",
	      i,PRB(1),PRB(2),PRB(3),PRB(4));
#undef PRB
  return bbb;
}

static void mk_diamond(t_atoms *a,rvec x[],real odist,t_symtab *symtab,
		       gmx_bool bPBC,matrix box)
{
  
  int   i,ib,j,k,l,m,nrm=0;
  t_bbb *bbb;
  gmx_bool  *bRemove;
  rvec  dx;
  
  do {
    nrm = 0;
    bbb = mk_bonds(a->nr,x,odist,bPBC,box);
    
    for(i=0; (i<a->nr); i++) {
      if (bbb[i].n < 2) {
	for(k=0; (k<bbb[i].n); k++) {
	  ib = bbb[i].aa[k];
	  for(j=0; (j<bbb[ib].n); j++)
	    if (bbb[ib].aa[j] == i)
	      break;
	  if (j == bbb[ib].n)
	    gmx_fatal(FARGS,"Bond inconsistency (%d not in list of %d)!\n",i,ib);
	  for( ; (j<bbb[ib].n-1); j++)
	    bbb[ib].aa[j] = bbb[ib].aa[j+1];
	  bbb[ib].n--;
	  nrm++;
	}
	bbb[i].n = 0;
      }
    }
  
    for(i=j=0; (i<a->nr); i++) {
      if (bbb[i].n >= 2) {
	copy_rvec(x[i],x[j]);
	j++;
      }
    }
    fprintf(stderr,"Kicking out %d carbon atoms (out of %d)\n",
	    a->nr-j,a->nr);
    a->nr = j;
    sfree(bbb);
  } while (nrm > 0);
  /* Rename atoms */
  bbb = mk_bonds(a->nr,x,odist,bPBC,box);
  for(i=0; (i<a->nr); i++) {
    switch (bbb[i].n) {
    case 4:
      a->atomname[i] = put_symtab(symtab,"C");
      break;
    case 3:
      a->atomname[i] = put_symtab(symtab,"CH1");
      break;
    case 2:
      a->atomname[i] = put_symtab(symtab,"CH2");
      break;
    default:
      gmx_fatal(FARGS,"This atom (%d) has %d bonds only",i,bbb[i].n);
    }
  }
  sfree(bbb);
}

static real calc_ener(real c6,real c12,rvec dx,tensor vir)
{
  real r,e,f;
  int  m,n;
  
  r  = norm(dx);
  e  = c12*pow(r,-12.0) - c6*pow(r,-6.0);
  f  = 12*c12*pow(r,-14.0) - 6*c6*pow(r,-8.0);
  for(m=0; (m<DIM); m++) 
    for(n=0; (n<DIM); n++)
      vir[m][n] -= 0.5*f*dx[m]*dx[n];
      
  return e;
}

static int read_rel_coords(char *fn,rvec **xx,int natmol)
{
  int    i,nline;
  double x,y,z;
  char   **strings=NULL;
  
  nline = get_file(fn,&strings);
  printf("Read %d lines from %s\n",nline,fn);
  snew(*xx,nline*natmol);
  for(i=0; (i<nline); i++) {
    if (sscanf(strings[i],"%lf%lf%lf",&x,&y,&z) != 3)
      gmx_fatal(FARGS,"Not enough arguments on line %d of file %s (should be 3)",i,fn);
    (*xx)[natmol*i][XX] = x;
    (*xx)[natmol*i][YY] = y;
    (*xx)[natmol*i][ZZ] = z;
  }
  return natmol*nline;
}

void virial(FILE *fp,gmx_bool bFull,int nmol,rvec x[],matrix box,real rcut,
	    gmx_bool bYaw,real q[],gmx_bool bLJ)
{
  int  i,j,im,jm,natmol,ik,jk,m,ninter;
  rvec dx,f,ftot,dvir,vir,pres,xcmi,xcmj,*force;
  real dx6,dx2,dx1,fscal,c6,c12,vcoul,v12,v6,vctot,v12tot,v6tot;
  t_pbc pbc;
  
  set_pbc(&pbc,box);
  fprintf(fp,"%3s   -  %3s: %6s %6s %6s  %6s %8s %8s %8s\n",
	  "ai","aj","dx","dy","dz","|d|","virx","viry","virz");
  clear_rvec(ftot);
  clear_rvec(vir);
  ninter = 0;
  vctot  = 0;
  v12tot = 0;
  v6tot  = 0;
  natmol = bYaw ? 5 : 3;
  snew(force,nmol*natmol);
  
  for(i=0; (i<nmol); i++) {
    im = natmol*i;
    /* Center of geometry */
    clear_rvec(xcmi);
    for(ik=0; (ik<natmol); ik++)
      rvec_inc(xcmi,x[im+ik]);
    for(m=0; (m<DIM); m++)
      xcmi[m] /= natmol;

    for(j=i+1; (j<nmol); j++) {
      jm = natmol*j;
      /* Center of geometry */
      clear_rvec(xcmj);
      for(jk=0; (jk<natmol); jk++)
	rvec_inc(xcmj,x[jm+jk]);
      for(m=0; (m<DIM); m++)
	xcmj[m] /= natmol;

      /* First check COM-COM distance */
      pbc_dx(&pbc,xcmi,xcmj,dx);
      if (norm(dx) < rcut) {
	ninter++;
	/* Neirest neighbour molecules! */
	clear_rvec(dvir);
	for(ik=0; (ik<natmol); ik++) {
	  for(jk=0; (jk<natmol); jk++) {
	    pbc_dx(&pbc,x[im+ik],x[jm+jk],dx);
	    dx2    = iprod(dx,dx);
	    dx1    = sqrt(dx2);
	    vcoul  = q[ik]*q[jk]*ONE_4PI_EPS0/dx1;
	    vctot += vcoul;
	    
	    if (bLJ) {
	      if (bYaw) {
		c6  = yaw_lj[ik][2*jk];
		c12 = yaw_lj[ik][2*jk+1];
	      }
	      else {
		c6  = spc_lj[ik][2*jk];
		c12 = spc_lj[ik][2*jk+1];
	      }
	      dx6    = dx2*dx2*dx2;
	      v6     = c6/dx6;
	      v12    = c12/(dx6*dx6);
	      v6tot -= v6;
	      v12tot+= v12;
	      fscal  = (vcoul+12*v12-6*v6)/dx2;
	    }
	    else
	      fscal  = vcoul/dx2;
	    for(m=0; (m<DIM); m++) {
	      f[m]     = dx[m]*fscal;
	      dvir[m] -= 0.5*dx[m]*f[m];
	    }
	    rvec_inc(force[ik+im],f);
	    rvec_dec(force[jk+jm],f);
	    /*if (bFull)
	      fprintf(fp,"%3s%4d-%3s%4d: %6.3f %6.3f %6.3f %6.3f"
		      " %8.3f %8.3f %8.3f\n",
		      watname[ik],im+ik,watname[jk],jm+jk,
		      dx[XX],dx[YY],dx[ZZ],norm(dx),
		      dvir[XX],dvir[YY],dvir[ZZ]);*/
	  }
	}
	if (bFull)
	  fprintf(fp,"%3s%4d-%3s%4d: "
		  " %8.3f %8.3f %8.3f\n",
		  "SOL",i,"SOL",j,dvir[XX],dvir[YY],dvir[ZZ]);
	rvec_inc(vir,dvir);
      }
    }
  }
  fprintf(fp,"There were %d interactions between the %d molecules (%.2f %%)\n",
	  ninter,nmol,(real)ninter/(0.5*nmol*(nmol-1)));
  fprintf(fp,"Vcoul: %10.4e  V12: %10.4e  V6: %10.4e  Vtot: %10.4e (kJ/mol)\n",
	  vctot/nmol,v12tot/nmol,v6tot/nmol,(vctot+v12tot+v6tot)/nmol);
  pr_rvec(fp,0,"vir ",vir,DIM,TRUE);
  
  for(m=0; (m<DIM); m++) 
    pres[m] = -2*PRESFAC/(det(box))*vir[m];
  pr_rvec(fp,0,"pres",pres,DIM,TRUE);
  pr_rvecs(fp,0,"force",force,natmol*nmol);
  sfree(force);
}



int main(int argc,char *argv[])
{
  static char *desc[] = {
    "[TT]mkice[tt] generates an ice crystal in the Ih crystal form which is the",
    "most stable form. The rectangular unitcell contains eight molecules",
    "and all oxygens are tetrahedrally coordinated.[PAR]",
    "If an input file is given it is interpreted as a series of oxygen",
    "coordinates the distance between which can be scaled by the odist flag.",
    "The program then adds hydrogens to the oxygens in random orientation",
    "but with proper OH distances and HOH angle. This feature allows to",
    "build water clusters based on oxygen coordinates only."
  };
  static int nx=1,ny=1,nz=1;
  static gmx_bool bYaw=FALSE,bLJ=TRUE,bFull=TRUE,bSeries=FALSE;
  static gmx_bool bOrdered=TRUE,bDiamond=FALSE,bPBC=TRUE;
  static real rcut=0.3,odist=0.274,hdist=0.09572;
  t_pargs pa[] = {
    { "-nx",    FALSE, etINT,  {&nx}, "nx" },
    { "-ny",    FALSE, etINT,  {&ny}, "ny" },
    { "-nz",    FALSE, etINT,  {&nz}, "nz" },
    { "-yaw",   FALSE, etBOOL, {&bYaw},
      "Generate virtual sites and shell positions" },
    { "-lj",    FALSE, etBOOL, {&bLJ},
      "Use LJ as well as coulomb for virial calculation" },
    { "-rcut",  FALSE,etREAL,  {&rcut},"Cut-off for virial calculations" },
    { "-full",  FALSE,etBOOL,  {&bFull},"Full virial output" },
    { "-odist", FALSE, etREAL, {&odist}, "Distance (nm) between oxygens" },
    { "-hdist", FALSE, etREAL, {&hdist}, "Bondlength (nm) for OH bond" },
    { "-diamond",FALSE,etBOOL, {&bDiamond}, "Make a diamond instead" },
    { "-pbc",   FALSE, etBOOL, {&bPBC},  "Make a periodic diamond" },
    { "-order", FALSE,etBOOL,  {&bOrdered}, "Make a proton-ordered ice lattice" },
    { "-series",FALSE, etBOOL, {&bSeries}, 
      "Do a series of virial calculations with different cut-off (from 0.3 up till the specified one)" }
  };
  t_filenm fnm[] = {
    { efSTO, "-p", "ice", ffWRITE },
    { efSTO, "-c", NULL,  ffOPTRD },
    { efDAT, "-f", NULL,  ffOPTRD },
    { efTRN, "-o", "ice", ffOPTWR }
  };
#define NFILE asize(fnm)

  FILE      *fp;
  char      *fn,quote[256];
  int       i,j,k,n,nmax,m,natom,natmol;
  t_atoms   *pdba;
  t_atoms   atoms;
  t_symtab  symtab;
  rvec      box,tmp,*xx;
  matrix    boxje;
  
  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,asize(pa),pa,asize(desc),
		    desc,0,NULL);
  if (debug) {
    fprintf(debug,"nx  = %3d, ny  = %3d,  nz   = %3d\n",nx,ny,nz);
    fprintf(debug,"YAW = %3s, LJ  = %3s,  rcut = %g\n",yesno_names[bYaw],
	    yesno_names[bLJ],rcut);
  }

  if (bYaw)
    natmol = 5;
  else if (bDiamond)
    natmol = 1;
  else
    natmol = 3;
    
  if (opt2bSet("-f",NFILE,fnm)) {
    natom = read_rel_coords(opt2fn("-f",NFILE,fnm),&xx,natmol);
    nmax  = natom;
  }
  else {
    natom = natmol*8;
    nmax = natom*nx*ny*nz;
    snew(xx,nmax);
  }
  snew(pdba,1);
  init_t_atoms(pdba,nmax,TRUE);
  pdba->nr = nmax;
  open_symtab(&symtab);
  for(n=0; (n<nmax); n++) {
    pdba->pdbinfo[n].type   = epdbATOM;
    pdba->pdbinfo[n].atomnr = 1+n;
    pdba->atom[n].resnr     = 1+(n/natmol);
    pdba->atomname[n] = put_symtab(&symtab,
				   bDiamond ? diamname[(n % natmol)] : watname[(n % natmol)]);
    if (bDiamond)
      pdba->resname[n] = put_symtab(&symtab,"DIA");
    else
      pdba->resname[n] = put_symtab(&symtab,"SOL");
    sprintf(pdba->pdbinfo[n].pdbresnr,"%d",n);
    pdba->atom[n].chain = ' ';
  }
  
  /* Generate the unit cell */
  if (bDiamond)
    unitcell_d(xx,box,odist); 
  else if (opt2bSet("-f",NFILE,fnm)) {
    random_h_coords(natmol,natom/natmol,xx,box,bYaw,odist,hdist);
  }
  else
    unitcell(xx,box,bYaw,odist,hdist);
  if (debug) {
    clear_mat(boxje);
    boxje[XX][XX] = box[XX];
    boxje[YY][YY] = box[YY];
    boxje[ZZ][ZZ] = box[ZZ];
  }
  n=0;
  for(i=0; (i<nx); i++) {
    tmp[XX] = i*box[XX];
    for(j=0; (j<ny); j++) {
      tmp[YY] = j*box[YY];
      for(k=0; (k<nz); k++) {
	tmp[ZZ] = k*box[ZZ];
	for(m=0; (m<natom); m++,n++) {
	  if ((!bOrdered && ((m % natmol) == 0)) || bOrdered)
	    rvec_add(xx[n % natom],tmp,xx[n]);
	  else
	    ;
	}
      }
    }
  }
    
  clear_mat(boxje);
  boxje[XX][XX] = box[XX]*nx;
  boxje[YY][YY] = box[YY]*ny;
  boxje[ZZ][ZZ] = box[ZZ]*nz;
  
  printf("Crystal:   %10.5f  %10.5f  %10.5f\n",
	 nx*box[XX],ny*box[YY],nz*box[ZZ]);
  
  if (debug && !bDiamond) {
    if (bSeries)
      for(i=3; (i<=10*rcut); i++) {
	fprintf(debug,"This is with rcut = %g\n",i*0.1);
	virial(debug,bFull,nmax/natmol,xx,boxje,
	       0.1*i,bYaw,bYaw ? qyaw : qspc,bLJ);
      }    
    else
      virial(debug,bFull,nmax/natmol,xx,boxje,
	     rcut,bYaw,bYaw ? qyaw : qspc,bLJ);
  }
  
  if (bDiamond) 
    mk_diamond(pdba,xx,odist,&symtab,bPBC,boxje);

  fn = ftp2fn(efSTO,NFILE,fnm);
  if (fn2ftp(fn) == efPDB) {
    fp = gmx_ffopen(fn,"w");
    if (bDiamond)
      fprintf(fp,"HEADER    This is a *diamond*\n");
    else
      fprintf(fp,"HEADER    A beautiful Ice Ih crystal\n");
    fprintf(fp,"REMARK    Generated by mkice with the following options:\n"
	    "REMARK    nx = %d, ny = %d, nz = %d, odist = %g, hdist = %g\n",
	    nx,ny,nz,odist,hdist);
	bromacs(quote,255);
    write_pdbfile(fp,quote,pdba,xx,boxje,' ',-1);
    gmx_ffclose(fp);
  }
  else {
    bromacs(quote,255);
    write_sto_conf(fn,quote,pdba,xx,NULL,boxje);
  }
  
  if (ftp2bSet(efTRN,NFILE,fnm))
    write_trn(ftp2fn(efTRN,NFILE,fnm),0,0,0,boxje,nmax,xx,NULL,NULL);
	       
  return 0;
}

