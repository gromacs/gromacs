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
 * S  C  A  M  O  R  G
 */
#include <math.h>
#include "sysstuff.h"
#include "string.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "xvgr.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "grid.h"
#include "statutil.h"
#include "rdgroup.h"
#include "ns.h"

static void gc1d(int icx,int nr,int delta,int *xmin,int *xmax)
{
  int dx;
  
  if (nr <= 2*delta+1) {
    *xmin=0;
    *xmax=nr-1;
  }
  else {
    dx=icx-delta;
    *xmin=(dx < 0) ? dx+nr : dx;
    *xmax=*xmin+min(2*delta,nr-1);
  }
}

void get_cells(FILE *log,t_grid *grid,int icg,
               int *lcxmin,int *lcxmax,
               int *lcymin,int *lcymax,
               int *lczmin,int *lczmax)
{
  int delta,icx,icy,icz;
  
  ci2xyz(grid,icg,&icx,&icy,&icz);

  delta = grid->delta;
  gc1d(icx,grid->nrx,delta,lcxmin,lcxmax);
  gc1d(icy,grid->nry,delta,lcymin,lcymax);
  gc1d(icz,grid->nrz,delta,lczmin,lczmax);

#ifdef DEBUG  
  fprintf(log,"icg=%d",icg);
  fprintf(log,"  lcx=(%2d-%2d)",*lcxmin,*lcxmax);
  fprintf(log,"  lcy=(%2d-%2d)",*lcymin,*lcymax);
  fprintf(log,"  lcz=(%2d-%2d)\n",*lczmin,*lczmax);
#endif
}


#define CELLSZ 0.6 /* nm */

void do_wi(FILE *log,FILE *out[],real t,
	   int nhb,atom_id hbind[],t_grid *grid,rvec x[],matrix box,
	   int nw,atom_id wind[])
{
  int    ci,cj,ci0,lx,ly,lz,cx,cy,cz,nrj,cgj0,jjcg;
  int    lcxmin,lcxmax,lcymin,lcymax,lczmin,lczmax,nrij;
  int    nrx,nry,nrz;
  ivec   nbox,cell,nmin;
  rvec   box_size,hbox,dmin,n_b,dx;
  real   d2,dHA,dDA;
  int    h,i,j,k,m,n,ai,ak;
  int    *grida,*gridnra,*gridind;
  
  nbox[XX]=grid->nrx,nbox[YY]=grid->nry,nbox[ZZ]=grid->nrz;
  nrx=grid->nrx,nry=grid->nry,nrz=grid->nrz;
  grida=grid->a;
  gridnra=grid->nra;
  gridind=grid->index;
  for(m=0; (m<DIM); m++) {
    box_size[m]=box[m][m];
    hbox[m]=0.5*box_size[m];
    n_b[m]=nbox[m]/box_size[m];
  }
  
  /* Calc index for hbond particles */
  for(n=0; (n<nhb); n++) {
    dmin[XX]=dmin[YY]=dmin[ZZ]=sqr(2*CELLSZ);
    nmin[XX]=nmin[YY]=nmin[ZZ]=-1;
    for(h=0; (h<3); h++) {
      ai=hbind[n*3+h];
      
      for(m=0; (m<DIM); m++) 
	cell[m]=((int)(n_b[m]*x[ai][m])+nbox[m]) % nbox[m];
  
      ci0=xyz2ci(nry,nrz,cell[XX],cell[YY],cell[ZZ]);
  
      get_cells(log,grid,ci0,&lcxmin,&lcxmax,&lcymin,&lcymax,&lczmin,&lczmax);
      for(lx=lcxmin; (lx <= lcxmax); lx++) {
	cx=(lx+nrx) % nrx;
	for(ly=lcymin; (ly <= lcymax); ly++) {
	  cy=(ly+nry) % nry;
	  for(lz=lczmin; (lz <= lczmax); lz++) {
	    cz=(lz+nrz) % nrz;
	    
	    cj=xyz2ci(nry,nrz,cx,cy,cz);
	    nrj=gridnra[cj];
	    cgj0=gridind[cj];
	    for (j=0; (j<nrj); j++) {
	      jjcg=grida[cgj0+j];
	      for(k=0; (k<3); k++) {
		ak=wind[3*jjcg+k];
		pbc_dx(box,x[ai],x[ak],dx);
		d2=iprod(dx,dx);
		if (d2 < dmin[h]) {
		  dmin[h]=d2;
		  nmin[h]=ak;
		}
	      }
	    }
	  }
	}
      }
    }
    /* Now we have the closest atoms for each of D,H,A */
    /* Find the closest molecule */
    for(m=0; (m<DIM); m++)
      dmin[m]=sqrt(dmin[m]);
    pbc_dx(box,x[hbind[3*n]],x[hbind[3*n+2]],dx);
    dDA=norm(dx);
    pbc_dx(box,x[hbind[3*n+1]],x[hbind[3*n+2]],dx);
    dHA=norm(dx);
    fprintf(out[n],"%10g  %10g  %10g  %10g  %10g  %10g  %5d  %5d  %5d\n",
	    t,dDA,dHA,dmin[0],dmin[1],dmin[2],nmin[0],nmin[1],nmin[2]);
  }
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "waterins assumes first group in index contains h-bonds in order",
    "D-H-A, and second group in index contains waters in order O-H-H.",
    "If this is not the case, you will die."
  };
  static int res0=1;
  t_pargs pa[] = {
    { "-r0", FALSE, etINT, &res0,
      "Sequence number of the first residue in your topology." }
  };
  static char *sets[] = { "D-A", "H-A", "D-S1", "H-S2", "A-S3", 
			    "S1 #","S2 #", "S3 #" };
#define MAXF 40
  FILE       *log,*out[MAXF];
  char       buf[256],bef[256],tit[256];
  t_topology *top;
  char       *grpname[2];
  int        status,gnx[2];
  rvec       *x,*cg_cm;
  matrix     box;
  int        natoms,nhb,nw,teller,i,wsize,rN,rO;
  real       t;
  t_grid     *grid;
  atom_id    *index[2],aN,aH,aO;
  t_filenm   fnm[] = {
    { efTRJ, "-f", NULL,  ffREAD },
    { efTPB, NULL, NULL,  ffREAD },
    { efNDX, NULL, NULL,  ffREAD },
    { efXVG, NULL, "dist",ffWRITE }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
		        
  top=read_top(ftp2fn(efTPB,NFILE,fnm));
  rd_index(ftp2fn(efNDX,NFILE,fnm),2,gnx,index,grpname);
  wsize=3;
  nhb=gnx[0]/3;
  nw=gnx[1]/wsize;
  fprintf(stderr,"There are %d hbonds\n",nhb);
  fprintf(stderr,"There are %d waters\n",nw);
  if (nhb > MAXF) {
    fprintf(stderr,"Can only handle %d hbonds at the same time\n",MAXF);
    fprintf(stderr,"Will truncate hbond list...\n");
    nhb=MAXF;
  }
  sprintf(bef,"%s",ftp2fn(efXVG,NFILE,fnm));
  bef[strlen(bef)-4]='\0';
  for(i=0; (i<nhb); i++) {
    sprintf(buf,"%s%d.xvg",bef,i);
   
    aN=index[0][3*i];
    aH=index[0][3*i+1];
    aO=index[0][3*i+2];
    rN=top->atoms.atom[aN].resnr;
    rO=top->atoms.atom[aO].resnr;
    sprintf(tit,"Water Insertion %s%d:%s%s-%s%d:%s",
	    *(top->atoms.resname[rN]),rN+res0,
	    *(top->atoms.atomname[aN]),
	    *(top->atoms.atomname[aH]),
	    *(top->atoms.resname[rO]),rO+res0,
	    *(top->atoms.atomname[aO])
	    );
    out[i]=xvgropen(buf,tit,"Time (ps)","Distance (nm)");
    xvgr_legend(out[i],asize(sets),sets);
  }
    
  /* Read first x */
  natoms=read_first_x(&status,ftp2fn(efTRJ,NFILE,fnm),&t,&x,box);
  
  /* Initiate grid stuff */
  log=ffopen("waterins.log","w");
  snew(cg_cm,nw);
  snew(grid,1);
  init_grid(log,grid,1,box,CELLSZ,nw);
  
  teller=0;
  do {
    if ((teller++ % 10) == 0)
      fprintf(stderr,"\rt=%.2f",t);
    /*put_atoms_in_box(log,0,nw,FALSE,box,box_size,cgs,x,shift_vec,cg_cm);*/
    for(i=0; (i<nw); i++)
      copy_rvec(x[index[1][3*i]],cg_cm[i]);
        
    /* Store everything on a grid */
    grid_first(log,grid,box,CELLSZ);
    fill_grid(log,grid,box,nw,0,nw,cg_cm);
    calc_elemnr(log,grid,0,nw,nw);
    calc_ptrs(grid);
    grid_last(log,grid,0,nw,nw);
    
    do_wi(log,out,t,nhb,index[0],grid,x,box,nw,index[1]);
  } while (read_next_x(status,&t,natoms,x,box));
  fprintf(stderr,"\n");
  
  close_trj(status);
  
  for(i=0; (i<nhb); i++) 
    fclose(out[i]);
  
  thanx(stdout);
  
  return 0;
}

