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
 * Good ROcking Metal Altar for Chronical Sinners
 */
static char *SRCID_g_rdens_c = "$Id$";

#include <math.h>
#include "sysstuff.h"
#include "string.h"
#include "string2.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "gstat.h"
#include "vec.h"
#include "xvgr.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "rdgroup.h"
#include "pbc.h"

real sphere_vol(real r)
{
  return (4.0*M_PI/3.0)*(r*r*r);
}

/* rdf_add just counts nr. of atoms in each slice and total mass 
   in each slice */
void rdf_add(int bin[], real mbin[],real width,real hb2,rvec x[],rvec xcm,
	     int nx2,atom_id index2[], t_atom atom[])
{
  int     j,ind,nout;
  int     jx;
  rvec    dx;
  real    r2;

  for(j=0; (j < nx2); j++) {
    jx=index2[j];
    rvec_sub(x[jx],xcm,dx);
    r2=iprod(dx,dx);
    if (r2 < hb2) { 
      ind=sqrt(r2)/width;
      mbin[ind]+=atom[jx].m;
      bin[ind]++;
    }
    else
      nout++;
  }
}

void rdf_calc(char *fn,char *pdens, char *rdens, char *ndens, 
	      real width, int n1,char *grp1,atom_id index1[],
	      int ngrps,int n2[],char **grpname,atom_id *index[],
	      t_atom atom[])
{
  FILE         *pout,*rout,*nout;
  int          **bin;
  real         **mbin; /* keeps mass per bin */
  char         buf[256];
  real         t,hb2,nf_1;
  int          i,j,natoms,status,maxbin;
  rvec         *x0,xcm;
  matrix       box;
  
  if ((natoms=read_first_x(&status,fn,&t,&x0,box))==0)
    fatal_error(0,"Could not read coordinates from statusfile\n");
  
  hb2=min(box[XX][XX],min(box[YY][YY],box[ZZ][ZZ]))/2;
  maxbin=(hb2/width)+2;
  snew(bin,ngrps);
  snew(mbin,ngrps);
  for(i=0; (i<ngrps); i++)
  {
    snew(bin[i],maxbin);
    snew(mbin[i],maxbin);
  }

  fprintf(stderr,"maxbin = %d, hb2 = %g\n",maxbin,hb2);
      
  j=0;
  do {
    if ((j % 10) == 0)
      fprintf(stderr,"\rframe: %5d",j);
      
    /* Must init pbc every step because of pressure coupling */
    init_pbc(box,FALSE);
    
    calc_xcm(x0,n1,index1,atom,xcm,FALSE);
    
    for(i=0; (i<ngrps); i++)
      rdf_add(bin[i],mbin[i],width,sqr(hb2),x0,xcm,n2[i],index[i],
	      atom);
    j++;
  } while (read_next_x(status,&t,natoms,x0,box));
  fprintf(stderr,"\n");
  close_trj(status);
  sfree(x0);
  nf_1=1.0/(real)j;
  
  /* now bin is a mega array with the total nr. of atoms of each type 
     in each bin. bin*nf_1 is the nr. of atoms of each type in each bin
     per frame.
     */
  for(i=0; (i<ngrps); i++) {
    sprintf(buf,"prob:%s-%s.xvg",grp1,grpname[i]);
    pout=xvgropen(buf,"Radial Probability Plot","R (nm)","nm\\S-1\\N");
    sprintf(buf,"dens:%s-%s.xvg",grp1,grpname[i]);
    rout=xvgropen(buf,"Radial Density Plot","R (nm)","g cm\\S-3\\N");
    sprintf(buf,"ndens:%s-%s.xvg",grp1,grpname[i]);
    nout=xvgropen(buf,"Radial Number Density Plot","R (nm)","Atoms nm\\S-3\\N");

    fprintf(pout,"@ subtitle \"%s-%s\"\n",grpname[0],grpname[1]);
    fprintf(nout,"@ subtitle \"%s-%s\"\n",grpname[0],grpname[1]);
    fprintf(rout,"@ subtitle \"%s-%s\"\n",grpname[0],grpname[1]);

    for(j=0; (j<maxbin); j++) {
      fprintf(nout,"%10g  %10g\n",width*j,bin[i][j]*nf_1/
		(4*M_PI*width*(j+0.5)*width*(j+0.5)));
      fprintf(pout,"%10g  %10g\n",width*j,bin[i][j]*nf_1);
      fprintf(rout,"%10g  %10g\n",width*j,mbin[i][j]*nf_1/
	      (60.22*4*M_PI*width*(j+0.5)*width*(j+0.5)));
    }
    fclose(nout); fclose(pout); fclose(rout);
  }
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "Compute radial densities across the box, in three flavors:"
    "probability density, number density, real density"
  };
  static real      width=0.12;
  t_pargs pa[] = {
    { "-width", FALSE, etREAL, &width, "bin width for radial axis" }
  };
  
  char      *grp1,**grpname;                /* groupnames                 */
  int       gnx1,*gnx;                      /* sizes of groups            */
  atom_id   *ind1,**index;             	    /* indices for all groups     */
  int       ngrps;
  t_topology *top;                	    /* topology 		  */ 
  t_filenm  fnm[] = {             	    /* files for g_order 	  */
    { efTRX, "-f", NULL,  ffREAD },    	    /* trajectory file 	          */
    { efNDX, NULL, NULL,  ffREAD },    	    /* index file 		  */
    { efTPX, NULL, NULL,  ffREAD },    	    /* topology file           	  */
    { efXVG,"-op","p_rdens", ffWRITE },     /* xvgr output: prob. dens.   */
    { efXVG,"-on","n_rdens", ffWRITE },     /* xvgr output: number. dens. */
    { efXVG,"-or","r_rdens", ffWRITE },     /* xvgr output: real dens. 	  */
  };

#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME,TRUE,NFILE,fnm,
		    asize(pa),pa,asize(desc),desc,0,NULL);
		      
  top = read_top(ftp2fn(efTPX,NFILE,fnm));     /* read topology file */
  fprintf(stderr,"Choose first group for Center of Mass computation!\n");
  
  fprintf(stderr,"Select group for Center of Mass calculation:\n");
  rd_index(ftp2fn(efNDX,NFILE,fnm),1,&gnx1,&ind1,&grp1);
  
  fprintf(stderr,"How many groups do you want to calc the density plot of?\n");
  do {
    scanf("%d",&ngrps);
  } while (ngrps < 0);
  
  snew(gnx,ngrps);
  snew(index,ngrps);
  snew(grpname,ngrps);
  fprintf(stderr,"Select groups for Density Plot:\n");
  rd_index(ftp2fn(efNDX,NFILE,fnm),ngrps,gnx,index,grpname);
  
  rdf_calc(ftp2fn(efTRX,NFILE,fnm),opt2fn("-op",NFILE,fnm), 
	   opt2fn("-or",NFILE,fnm),opt2fn("-on",NFILE,fnm), 
	   width,gnx1,grp1,ind1,ngrps,gnx,grpname,index,top->atoms.atom);
  
  thanx(stdout);
  
  return 0;
}



