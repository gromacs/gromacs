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
 * Giant Rising Ordinary Mutants for A Clerical Setup
 */
#include "copyrite.h"
#include "string2.h"
#include "smalloc.h"
#include "sysstuff.h"
#include "confio.h"
#include "paramio.h"
#include "statutil.h"
#include "physics.h"
#include "vec.h"
#include "random.h"
#include "pbc.h"
#include "txtdump.h"
#include "gbutil.h"
#include "3dview.h"

#define GCIN \
  RTYPE("dx",                   dx,      0.15) \
  RTYPE("dy",                   dy,      0.15) \
  RTYPE("dz",                   dz,      0.15) \
  ITYPE("seed",                 seed,    1993) \
  NULL

static void rand_rot(int natoms,rvec x[],rvec v[],vec4 xrot[],vec4 vrot[],
                     int *seed)
{
  mat4 mt1,mt2,mr[DIM],mtemp1,mtemp2,mtemp3,mxtot,mvtot;
  rvec xcm;
  real phi;
  int  i,m;
  
  clear_rvec(xcm);
  for(i=0; (i<natoms); i++)
    for(m=0; (m<DIM); m++) {
      xcm[m]+=x[i][m]/natoms;
    }
  
  translate(-xcm[XX],-xcm[YY],-xcm[ZZ],mt1);
  for(m=0; (m<DIM); m++) {
    phi=2*M_PI*rando(seed);
    rotate(m,phi,mr[m]);
  }
  translate(xcm[XX],xcm[YY],xcm[ZZ],mt2);

  mult_matrix(mtemp1,mt2,mr[XX]);
  mult_matrix(mtemp2,mr[YY],mr[ZZ]);
  mult_matrix(mtemp3,mtemp1,mtemp2);
  mult_matrix(mxtot,mtemp3,mt1);
  mult_matrix(mvtot,mr[XX],mtemp2);
  
  for(i=0; (i<natoms); i++) {
    m4_op(mxtot,x[i],xrot[i]);
    m4_op(mvtot,v[i],vrot[i]);
  }
}

static void move_x(int natoms,rvec x[],matrix box)
{
  int  i,m;
  rvec xcm;

  clear_rvec(xcm);
  for(i=0; (i<natoms); i++)
    for(m=0; (m<DIM); m++)
      xcm[m]+=x[i][m];
  for(m=0; (m<DIM); m++)
    xcm[m] = 0.5*box[m][m]-xcm[m]/natoms;
  for(i=0; (i<natoms); i++)
    for(m=0; (m<DIM); m++)
      x[i][m] += xcm[m];
}

void readconf(char *fn,t_atoms *atoms,rvec **x,rvec **v,matrix box)
{
  int natoms;
  char title[256];
  
  get_coordnum(fn,&natoms);
  snew(atoms->resname,natoms);
  snew(atoms->atomname,natoms);
  snew(atoms->atom,natoms);
  snew(*x,natoms);
  snew(*v,natoms);
  read_whole_conf(fn,title,atoms,*x,*v,box);
  fprintf(stderr,"Read \"%s\"\n",title);
}

static void y_rot(int natoms,rvec x[],rvec v[])
{
  mat4 mt1,mt2,mv1,mv2,mr,mtemp1,mxtot,mvtot;
  vec4 xrot,vrot;
  rvec xcm,vcm;
  int  i,m;
  
  clear_rvec(xcm);
  for(i=0; (i<natoms); i++)
    for(m=0; (m<DIM); m++) {
      xcm[m]+=x[i][m]/natoms;
      vcm[m]+=v[i][m]/natoms;
    }
  
  translate(-xcm[XX],-xcm[YY],-xcm[ZZ],mt1);
  translate(-vcm[XX],-vcm[YY],-vcm[ZZ],mv1);
  rotate(XX,DEG2RAD*90,mr);
  translate(xcm[XX],xcm[YY],xcm[ZZ],mt2);
  translate(vcm[XX],vcm[YY],vcm[ZZ],mv2);
  mult_matrix(mtemp1,mt2,mr);
  mult_matrix(mxtot,mtemp1,mt1);
  
  mult_matrix(mtemp1,mv2,mr);
  mult_matrix(mvtot,mtemp1,mv1);
  
  for(i=0; (i<natoms); i++) {
    m4_op(mxtot,x[i],xrot);
    m4_op(mvtot,v[i],vrot);
    for(m=0; (m<DIM); m++) {
      x[i][m]=xrot[m];
      v[i][m]=vrot[m];
    }
  }
}

static void orient_pep(int natoms,rvec x[],rvec v[],matrix box,
		       bool bTrans)
{
  rvec angle;
  rvec xmin,xmax;
  int  i,m;
  
  init_pbc(box,FALSE);
  
  clear_rvec(angle);
  orient(natoms,x,v,angle,box);
  if (bTrans) {
    /* Orient has put the peptide along the Z-axis, now we'll rotate it
     * over 90 deg to be along the Y-axis 
     */
     fprintf(stderr,"Now orienting along Y-axis...\n");
     y_rot(natoms,x,v);
  }
  fprintf(stderr,"Recomputing box...\n");
  copy_rvec(x[0],xmin);
  copy_rvec(x[0],xmax);
  for(i=1; (i<natoms); i++) {
    for(m=0; (m<DIM); m++) {
      xmin[m]=min(xmin[m],x[i][m]);
      xmax[m]=max(xmax[m],x[i][m]);
    }
  }
  fprintf(stderr,"Putting peptide in first octant of coord system...\n");
  for(i=0; (i<natoms); i++)
    rvec_dec(x[i],xmin);
  rvec_dec(xmax,xmin);
  for(m=0; (m<DIM); m++)
    box[m][m]=xmax[m];
  fprintf(stderr,"Box is: %g,%g,%g\n",xmax[XX],xmax[YY],xmax[ZZ]);
}

void main(int argc, char *argv[])
{
  static char *desc[] = {
    "mk_membr puts a peptide [BB]in[bb] or [BB]on[bb] a membrane.",
    "It is assumed that the lipid file contains a lipid or oil layer",
    "with the interface in the X-Y plane. This means that a trans membrane",
    "peptide will be placed with the primary axis (eg. the helical axis)",
    "in the Z-direction,",
    "a parallel peptide will be in the X-Y plane."
  };
  static char *opts[] = {
    "-trans"
  };
  static char *odesc[] = {
    "Make a trans-membrane peptide, instead of parallel to the membrane."
  };
  t_manual man = {asize(desc),desc,asize(opts),opts,odesc,0,NULL};

  int     nx,ny,nz,vol;
  
  t_atoms atoms_p,atoms_l;
  rvec    *x_p,*v_p,*x_l,*v_l;
  real    *r_l;
  matrix  box_p,box_l;
  
  real    dx,dy,dz;
  bool    bTrans=FALSE;
  rvec    shift;
  int     seed;
  ivec    n_box;
  int     natoms,nres;
  int     i,j,k,l,m,ndx,nrdx;
  t_filenm fnm[] = {
    { efGCP, NULL,  NULL,       ffREAD  },
    { efGCP, "-po", "gcout",    ffWRITE },
    { efGRO, "-cp", "conf_pep", ffREAD  },
    { efGRO, "-cl", "conf_lip", ffREAD  },
    { efGRO, "-co", "confout",  ffWRITE }
  };
#define NFILE asize(fnm)
    
  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,FALSE,&man);

  for(i=1; (i < argc); i++) {
    if (strcmp(argv[i],"-trans") == 0)
      bTrans=TRUE;
  }
  if (bTrans)
    fprintf(stderr,"Will try to make a trans-membrane peptide\n");
  
  read_params (ftp2fn(efGCP,NFILE,fnm),GCIN);
  write_params(opt2fn("-po",NFILE,fnm),GCIN);

  readconf(opt2fn("-cp",NFILE,fnm),&atoms_p,&x_p,&v_p,box_p);
  orient_pep(atoms_p.nr,x_p,v_p,box_p,bTrans);
  
  readconf(opt2fn("-cl",NFILE,fnm),&atoms_l,&x_l,&v_l,box_l);
  if (!bTrans) {
    n_box[XX]=1+(box_l[XX][XX]/(box_p[XX][XX]+2*dx));
    n_box[YY]=1+(box_l[YY][YY]/(box_p[YY][YY]+2*dy));
    n_box[ZZ]=1+(box_l[ZZ][ZZ]/(box_p[ZZ][ZZ]+2*dz));
    vol=n_box[XX]*n_box[YY]*n_box[ZZ];
    natoms=vol*atoms_l.nr+atoms_p.nr;
    
    fprintf(stderr,"Lipid box will be multiplied %d times\n",vol);
    
    srenew(atoms_p.resname,natoms);
    srenew(atoms_p.atomname,natoms);
    srenew(atoms_p.atom,natoms);
    srenew(x_p,natoms);
    srenew(v_p,natoms);
    snew(r_p,natoms);
    genconf(atoms_p,x_p,v_p,r_p,box_p,n_box);
    
  }
  
  write_conf(opt2fn("-co",NFILE,fnm),"Peppie",&atoms_p,x_p,v_p,box_p);
  
  exit(1);
  
  
}

