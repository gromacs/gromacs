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
 * Great Red Owns Many ACres of Sand 
 */
static char *SRCID_genion_c = "$Id$";

#include "copyrite.h"
#include "string2.h"
#include "smalloc.h"
#include "sysstuff.h"
#include "confio.h"
#include "assert.h"
#include "statutil.h"
#include "pbc.h"
#include "force.h"
#include "fatal.h"
#include "futil.h"
#include "maths.h"
#include "macros.h"
#include "physics.h"
#include "vec.h"

void get_params(char *giin,char *giout,
		int *p_num,char p_name[],real *p_q,
		int *n_num,char n_name[],real *n_q,
		int *w1,int *nw,real *rcut)
{
  t_inpfile *inp;
  char      *tmp;
  int       ninp;

  inp=read_inpfile(giin,&ninp);
  ITYPE("n+",           *p_num,      1);
  STYPE("name+",    	p_name,	    "");
  RTYPE("q+",		*p_q,	    1.0);
  ITYPE("n-",           *n_num,       1);
  STYPE("name-",	n_name,	    "");
  RTYPE("q-",		*n_q,	   -1.0);
  ITYPE("water-1",	*w1,	    0);
  ITYPE("n-water",	*nw,	    0);
  RTYPE("rexcl",	*rcut,	    0.5);
  
  write_inpfile(giout,ninp,inp);
}

typedef struct {
  int     nr;		/* The number of neighbours	*/
  atom_id *j;		/* The array of j-part		*/
} t_nl;

static real dr2(rvec xi,rvec xj,matrix box)
{
  rvec r_ij;
  real rn;
  
  pbc_dx(box,xi,xj,r_ij);
  rn=iprod(r_ij,r_ij);

  return rn;
}

static void insert(t_nl *nl,int j)
{
#define HOPS 64 

  if ((nl->nr % HOPS)==0) 
    srenew(nl->j,nl->nr+HOPS);
  nl->j[nl->nr++]=j;
}

static void mk_nl(int nw,int index[],matrix box,rvec x[],real rcut,t_nl nl[])
{

  int  i,j,ii,jj;
  real rij2,rc2=rcut*rcut;
  real avnb;

  fprintf(stderr,"Making neighbourlist with cut-off %g for exclusions\n",rcut);
  for(i=0; (i<nw); i++) {
    ii=index[i];
    for(j=i+1; (j<nw); j++) {
      jj=index[j];
      rij2=dr2(x[ii],x[jj],box);
      if (rij2 < rc2) {
	insert(&(nl[i]),j);
	insert(&(nl[j]),i);
      }
    }
  }
  avnb=0;
  for(i=0; (i<nw); i++)
    avnb+=nl[i].nr;
  avnb/=nw;
  fprintf(stderr,"Average number of neighbours: %g\n",avnb);
}

static int calc_pot(char *infile,
		    t_topology *top,rvec **x0,rvec **v0,
		    real **coulomb,matrix box,real *rl2)
{
  FILE           *status;
  t_inputrec     ir;
  t_statheader   sh;
  
  rvec        *x,*x_s;
  real        *q,qi,rij,rij2,rlong2;
  real        *coul;
  int         i,j;

  /* Some dummies */
  int         step,natoms;
  real        t;

  status=ffopen(infile,"r");
  rd_header(status,&sh);
  snew(x,sh.natoms);
  snew(*v0,sh.natoms);
  snew(x_s,sh.natoms);
  fprintf(stderr,"Read statusfile version %s\n",
	  rd_hstatus(status,&sh,&step,&t,&t,
		     &ir,box,NULL,NULL,
		     &natoms,x,*v0,NULL,&step,NULL,top));
  fclose(status);

  /* Calc the force */
  fprintf(stderr,"Doing single force calculation...\n");

  snew(coul,sh.natoms);
  snew(q,natoms);
  init_pbc(box,FALSE);
  for(i=0; (i<natoms); i++)
    q[i]=top->atoms.atom[i].q;

  rlong2=ir.rlong*ir.rlong;
  *rl2=rlong2;
  for(i=0; (i<natoms); i++) {
    qi=q[i];
    for(j=i+1; (j<natoms); j++) {
      rij2=dr2(x[i],x[j],box);
      if (rij2 < rlong2) {
	rij=invsqrt(rij2);
	coul[i]+=q[j]*rij;
	coul[j]+=qi*rij;
      }
    }
  }
  sfree(q);

  *coulomb=coul;
  *x0=x;

  return sh.natoms;
}

static int *mk_windex(int w1,int nw)
{
  int *index;
  int i;

  snew(index,nw);
  for(i=0; (i<nw); i++)
    index[i]=w1+3*i;
  
  return index;
}

static void update_coul(real q,int nw,int index[],int ei,
			real coulomb[],rvec x[],matrix box,real rlong2)
{
  int  i,ii,j;
  real rij2;
  real *xow;
  real qO=-0.82*ONE_4PI_EPS0;
  real qH;

  qH=-qO/2;
  xow=x[ei];

  q*=ONE_4PI_EPS0;
  for(i=0; (i<nw); i++) {
    ii=index[i];
    
    if (ii != ei) {
      rij2=dr2(xow,x[ii],box);
      if (rij2 < rlong2) {
	coulomb[ii]+=(q-qO)*invsqrt(rij2);
	for(j=1; (j<3); j++) {
	  rij2=dr2(x[ei+j],x[ii],box);
	  coulomb[ii]-=qH*invsqrt(rij2);
	}
      }
    }
  }
}

static void insert_ion(real q,int nw,t_nl nl[],bool bSet[],
		       int index[],int nSubs[],real coulomb[],
		       rvec x[],matrix box,real rlong2,int natoms)
{
  int  i,ii,ei;
  real extr_e,qii;
  bool bSub=FALSE;

  ei=-1;

  for(i=0; (i<nw); i++) {
    if (!bSet[i]) {
      ii=index[i];
      if (nSubs[ii] == 0) {
	qii=coulomb[ii];
	if (q > 0) {
	  if ((qii <= extr_e) || !bSub) {
	    extr_e=qii;
	    ei=i;
	    bSub=TRUE;
	  }
	}
	else {
	  if ((qii >= extr_e) || !bSub) {
	    extr_e=qii;
	    ei=i;
	    bSub=TRUE;
	  }
	}
      }
    }
  }
  if (ei == -1)
    fatal_error(0,"No More replacable waters!");

  /* Now update the coulomb energy... */
  update_coul(q,nw,index,index[ei],coulomb,x,box,rlong2);
  nSubs[index[ei]]=(q < 0) ? -1 : 1;
  for(i=0; (i<nl[ei].nr); i++) {
    int jj=nl[ei].j[i];
    assert((jj >= 0) && (jj < nw));
    bSet[jj]=TRUE;
  }
}

static void p_xvn(FILE *out,int i,rvec x[],rvec v[])
{
  int m;

  for(m=0; (m<DIM); m++)
    fprintf(out,"%8.3f",x[i][m]);
  for(m=0; (m<DIM); m++)
    fprintf(out,"%8.4f",v[i][m]);
  fprintf(out,"\n");
  fflush(out);
}

static void p_atom(FILE *out,int i,t_atoms *atoms,rvec x[],rvec v[],
		   int anr,int rnr)
{
  int resnr;

  resnr=atoms->atom[i].resnr;
  fprintf(out,"%5d%-5.5s%5.5s%5d",
	  rnr+1,*(atoms->resname[resnr]),
	  *(atoms->atomname[i]),anr+1);
  p_xvn(out,i,x,v);
}

static void print_nsub(char *outfile,int nion,int w1,int nw,
		       int index[],int nSubs[],
		       char *p_name,char *n_name,
		       t_topology *top,rvec x[],rvec v[],
		       matrix box)
{
  FILE        *out;
  t_atoms     *atoms;
  int         i,j,m,ii;
  int         aind,rind;

  out=ffopen(outfile,"w");
  atoms=&(top->atoms);
  fprintf(out,"%s\n",*top->name);
  fprintf(out,"%5d\n",atoms->nr-2*nion);
  for(i=0; (i<w1); i++) 
    p_atom(out,i,atoms,x,v,i,atoms->atom[i].resnr);
  
  aind=w1;
  rind=atoms->atom[aind].resnr;
  for(i=0; (i<nw); i++) {
    ii=index[i];
    if (nSubs[ii] == 0) {
      for(m=0; (m<3); m++) 
	p_atom(out,ii+m,atoms,x,v,aind++,rind);
      rind++;
    }
  }
  for(i=0; (i<nw); i++) {
    ii=index[i];
    if (nSubs[ii] > 0) {
      fprintf(stdout,"Water Molec. %5d has been replaced by %s\n",i,p_name);
      fprintf(out,"%5d%-5.5s%5.5s%5d",
	      rind+1,p_name,p_name,aind+1);
      p_xvn(out,ii,x,v);
      rind++;
      aind++;
    }
  }
  for(i=0; (i<nw); i++) {
    ii=index[i];
    if (nSubs[ii] < 0) {
      fprintf(stdout,"Water Molec. %5d has been replaced by %s\n",i,n_name);
      fprintf(out,"%5d%-5.5s%5.5s%5d",
	      rind+1,n_name,n_name,aind+1);
      p_xvn(out,ii,x,v);
      rind++;
      aind++;
    }
  }
  for(i=aind+2*nion; (i<atoms->nr); i++) {
    rind=atoms->atom[i].resnr;
    p_atom(out,i,atoms,x,v,aind++,rind);
  }
  for(i=0; (i<DIM); i++) 
    fprintf(out,"%10.5f",box[i][i]);
  for(i=0; (i<DIM); i++)
    for(j=0; (j<DIM); j++)
      if (i != j)
	fprintf(out,"%10.5f",box[i][j]);
  fprintf(out,"\n");
  fclose(out);
}

int main(int argc, char *argv[])
{
  static char *desc[] = {
    "Generate ions in the positions of water molecules. To avoid",
    "clustering of ions it is advisable to set rexcl (the radius around",
    "an ion in which no other ion will be placed) to a value high enough",
    "to allow solvent around the ion (> 0.6 nm)."
  };
  static char *bugs[] = {
    "Only monatomic ions can be used. For larger ions, e.g. sulfate we recommended to use genbox."
  };
  
  int         p_num,n_num;
  char        p_name[STRLEN],n_name[STRLEN];
  real        p_q,n_q,rcut;
  t_topology  top;
  rvec        *x,*v;
  real        *coulomb;
  real        rlong2;
  matrix      box;
  int         w1,nw;
  int         *index;
  int         *nSubs;
  bool        *bSet;
  t_nl        *nl;
  int         i,nion,natoms;
  t_filenm fnm[] = {
    { efGIP, "-f",  NULL,     ffREAD },
    { efGIP, "-po", "gi-out", ffWRITE },
    { efTPB, NULL,  NULL,     ffREAD },
    { efGRO, "-o",  NULL,     ffWRITE }
  };
#define NFILE asize(fnm)
  
  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,TRUE,NFILE,fnm,0,NULL,asize(desc),desc,
		    asize(bugs),bugs);

  get_params(opt2fn("-f",NFILE,fnm),opt2fn("-po",NFILE,fnm),
	     &p_num,p_name,&p_q,&n_num,n_name,&n_q,&w1,&nw,&rcut);
  
  nion=p_num+n_num;

  natoms=calc_pot(fnm[2].fn,&top,&x,&v,&coulomb,box,&rlong2);
  index=mk_windex(w1,nw);
  snew(nSubs,natoms);
  snew(bSet,nw);
  snew(nl,nw);
  mk_nl(nw,index,box,x,rcut,nl);

  do {
    for(i=0; (i<nw); i++)
      printf("coul[%5d]=%g\n",i,coulomb[index[i]]);
    if (p_num >= n_num) {
      insert_ion(p_q,nw,nl,bSet,index,nSubs,coulomb,x,box,rlong2,natoms);
      fprintf(stderr,"+");
      p_num--;
    }
    else {
      insert_ion(n_q,nw,nl,bSet,index,nSubs,coulomb,x,box,rlong2,natoms);
      fprintf(stderr,"-");
      n_num--;
    }
  } while (p_num+n_num > 0);
  fprintf(stderr,"\n");

  print_nsub(fnm[3].fn,nion,w1,nw,index,nSubs,p_name,n_name,&top,x,v,box);

  thanx(stdout);
  
  return 0;
}
