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
 * S  C  A  M  O  R  G
 */
static char *SRCID_tpbconv_c = "$Id$";

#include <math.h>
#include "rdgroup.h"
#include "fatal.h"
#include "string2.h"
#include "sysstuff.h"
#include "smalloc.h"
#include "macros.h"
#include "names.h"
#include "typedefs.h"
#include "tpxio.h"
#include "trnio.h"
#include "readir.h"
#include "statutil.h"
#include "copyrite.h"
#include "futil.h"
#include "assert.h"
#include "vec.h"
#include "rdgroup.h"


static bool *bKeepIt(int gnx,int natoms,atom_id index[])
{
  bool *b;
  int  i;
  
  snew(b,natoms);
  for(i=0; (i<gnx); i++) {
    assert(index[i] < natoms);
    b[index[i]] = TRUE;
  }
  
  return b;
}

static atom_id *invind(int gnx,int natoms,atom_id index[])
{
  atom_id *inv;
  int     i;
  
  snew(inv,natoms);
  for(i=0; (i<gnx); i++) {
    assert(index[i] < natoms);
    inv[index[i]] = i;
  }
  
  return inv;
}

static void reduce_block(atom_id invindex[],bool bKeep[],t_block *block,
			 char *name,bool bExcl)
{
  atom_id *index,*a;
  int i,j,k,newi,newj;
  
  snew(index,block->nr);
  snew(a,block->nra);
  
  newi = newj = 0;
  for(i=0; (i<block->nr); i++) {
    if (!bExcl || bKeep[i]) {
      for(j=block->index[i]; (j<block->index[i+1]); j++) {
	k = block->a[j];
	if (bKeep[k]) {
	  a[newj] = invindex[k];
	  newj++;
	}
      }
      if (newj > index[newi]) {
	newi++;
	index[newi] = newj;
      }
    }
  }
  
  fprintf(stderr,"Reduced block %8s from %6d to %6d index-, %6d to %6d a-entries\n",
	  name,block->nr,newi,block->nra,newj);
  block->index = index;
  block->a     = a;
  block->nr    = newi;
  block->nra   = newj;
}

static void reduce_rvec(int gnx,atom_id index[],rvec vv[])
{
  rvec *ptr;
  int  i;
  
  snew(ptr,gnx);
  for(i=0; (i<gnx); i++)
    copy_rvec(vv[index[i]],ptr[i]);
  for(i=0; (i<gnx); i++)
    copy_rvec(ptr[i],vv[i]);
  sfree(ptr);
}

static void reduce_atom(int gnx,atom_id index[],t_atom atom[],char ***atomname,
			int *nres, char ***resname)
{
  t_atom *ptr;
  char   ***aname,***rname;
  int    i,nr;
  
  snew(ptr,gnx);
  snew(aname,gnx);
  snew(rname,atom[index[gnx-1]].resnr+1);
  for(i=0; (i<gnx); i++) {
    ptr[i]   = atom[index[i]];
    aname[i] = atomname[index[i]];
  }
  nr=-1;   
  for(i=0; (i<gnx); i++) {
    atom[i]     = ptr[i];
    atomname[i] = aname[i];
    if ((i==0) || (atom[i].resnr != atom[i-1].resnr)) {
      nr++;
      rname[nr]=resname[atom[i].resnr];
    }
    atom[i].resnr=nr;
  }
  nr++;
  for(i=0; (i<nr); i++)
    resname[i]=rname[i];
  *nres=nr;

  sfree(aname);
  sfree(ptr);
  sfree(rname);
}

static void reduce_ilist(atom_id invindex[],bool bKeep[],
			 t_ilist *il,int nratoms,char *name)
{
  t_iatom *ia;
  int i,j,newnr;
  bool bB;

  if (il->nr) {  
    snew(ia,il->nr);
    newnr=0;
    for(i=0; (i<il->nr); i+=nratoms+1) {
      bB = TRUE;
      for(j=1; (j<=nratoms); j++) {
	bB = bB && bKeep[il->iatoms[i+j]];
      }
      if (bB) {
	ia[newnr++] = il->iatoms[i];
	for(j=1; (j<=nratoms); j++)
	  ia[newnr++] = invindex[il->iatoms[i+j]];
      }
    }
    fprintf(stderr,"Reduced ilist %8s from %6d to %6d entries\n",
	    name,il->nr/(nratoms+1),
	  newnr/(nratoms+1));
    
    il->nr = newnr;
    for(i=0; (i<newnr); i++)
      il->iatoms[i] = ia[i];
    
    for(i=0; (i<MAXPROC); i++)
      il->multinr[i] = newnr;
  
    sfree(ia);
  }
}

static void reduce_topology_x(int gnx,atom_id index[],
			      t_topology *top,rvec x[],rvec v[])
{
  bool    *bKeep;
  atom_id *invindex;
  int     i;
  
  bKeep    = bKeepIt(gnx,top->atoms.nr,index);
  invindex = invind(gnx,top->atoms.nr,index);
  
  for(i=0; (i<ebNR); i++)
    reduce_block(invindex,bKeep,&(top->blocks[i]),eblock_names[i],FALSE);
  reduce_block(invindex,bKeep,&(top->atoms.excl),"EXCL",TRUE);
  reduce_rvec(gnx,index,x);
  reduce_rvec(gnx,index,v);
  reduce_atom(gnx,index,top->atoms.atom,top->atoms.atomname,
	      &(top->atoms.nres),top->atoms.resname);

  for(i=0; (i<F_NRE); i++) {
    reduce_ilist(invindex,bKeep,&(top->idef.il[i]),
		 interaction_function[i].nratoms,
		 interaction_function[i].name);
  }
    
  top->atoms.nr = gnx;
}

int main (int argc, char *argv[])
{
  static char *desc[] = {
    "tpbconv can edit tpx files in multiple ways.[PAR]"
    "[BB]1st.[bb] by creating a run input file",
    "for a continuation run when your simulation has crashed due to e.g.",
    "a full disk, or by making a continuation run input file.",
    "Note that both velocities and coordinates are needed,",
    "which means that when you never write velocities, you can not use",
    "tpbconv and you have to start the run again from the beginning.[PAR]",
    "[BB]2nd.[bb] by creating a tpx file for a subset of your original",
    "tpx file, which is useful when you want to remove the solvent from",
    "your tpx file, or when you want to make e.g. a pure Ca tpx file.",
    "[BB]WARNING: this tpx file is not fully functional[bb].[PAR]",
    "[BB]3rd.[bb] If the optional mdp file is [BB]explicitly[bb] given, ",
    "the inputrec",
    "will be read from there, and some options may be overruled this way.",
    "Note that this is a fast way to make a continuation tpx",
    "as the topology does not have to be processed again."
  };

  char         *fn;
  int          fp;
  t_tpxheader  tpx;
  int          i,nstep,natoms;
  real         t,newt,lambda;
  bool         bCont;
  t_topology   top;
  t_inputrec   *ir,*irnew;
  t_gromppopts *gopts;
  rvec         *x=NULL,*v=NULL,*newx,*newv,*tmpx,*tmpv;
  matrix       box,newbox,tmpbox;
  int          gnx;
  char         *grpname;
  atom_id      *index=NULL;
  t_filenm fnm[] = {
    { efMDP, "-pi",  NULL,   ffOPTRD },
    { efMDP, "-po", "mdout", ffOPTWR },
    { efTRN, "-f",  NULL,    ffOPTRD },
    { efTPX, NULL,  NULL,    ffREAD  },
    { efTPX, "-o",  "tpxout",ffWRITE },
    { efNDX, NULL,  NULL,    ffOPTRD },
  };
#define NFILE asize(fnm)

  /* Command line options */
  static real max_t = -1.0;
  static bool bSel=FALSE;
  static t_pargs pa[] = {
    { "-time", FALSE, etREAL, &max_t, 
      "time in the trajectory from which you want to continue." },
  };
  
  CopyRight(stdout,argv[0]);
  
  /* Parse the command line */
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);

  bSel = (bSel || ftp2bSet(efNDX,NFILE,fnm));
  
  fn = ftp2fn(efTPX,NFILE,fnm);
  printf("Reading toplogy and shit from %s\n",fn);
  
  read_tpxheader(fn,&tpx);
  snew(x,tpx.natoms);
  snew(v,tpx.natoms);
  snew(ir,1);
  read_tpx(fn,&nstep,&t,&lambda,ir,box,&natoms,x,v,NULL,&top);
  
  bCont=(ftp2bSet(efMDP,NFILE,fnm));
  if (bCont) {
    printf("Reading new inputrec from %s\n",ftp2fn(efMDP,NFILE,fnm));
    
    /* Initiate some variables */
    snew(irnew,1);
    snew(gopts,1);
    init_ir(irnew,gopts);
    
    /* PARAMETER file processing */
    get_ir(ftp2fn(efMDP,NFILE,fnm),opt2fn("-po",NFILE,fnm),irnew,gopts);
  }
  
  if (ftp2bSet(efTRN,NFILE,fnm)) {
    fn = ftp2fn(efTRN,NFILE,fnm);
    printf("\nREADING COORDS, VELS AND BOX FROM TRAJECTORY %s...\n\n",fn);
    
    sfree(x);
    sfree(v);
    natoms = read_first_x_v(&fp,ftp2fn(efTRN,NFILE,fnm),&t,&x,&v,box);
    snew(newx,tpx.natoms);
    snew(newv,tpx.natoms);

    if (top.atoms.nr != natoms) 
      fatal_error(0,"Number of atoms in Topology (%d) "
		  "is not the same as in Trajectory (%d)\n",
		  top.atoms.nr,natoms);
    
    /* Now scan until the last set of x and v (step == 0)
     * or the ones at step step.
     */
    while (read_next_x_v(fp,&newt,natoms,newx,newv,newbox)) {
      tmpx=newx;
      newx=x;
      x=tmpx;
      tmpv=newv;
      newv=v;
      v=tmpv;
      t=newt;
      copy_mat(newbox,tmpbox);
      copy_mat(box,newbox);
      copy_mat(tmpbox,box);

      if ((max_t != -1.0) && (t >= max_t))
	break;
    }
    close_trj(fp);
    fprintf(stderr,"\n");
  } else 
    printf("\nUSING COORDS, VELS AND BOX FROM TPX FILE %s...\n\n",
	   ftp2fn(efTPX,NFILE,fnm));
  
  /* change the input record to the actual data */
  if (bCont) {
    ir=irnew;
    if (ir->init_t != t)
      printf("WARNING: I'm using t_init in mdp-file (%g)"
	     " while in trajectory it is (%g)\n",ir->init_t,t);
    if (ir->init_lambda != lambda)
      printf("WARNING: I'm using lambda_init in mdp-file (%g)"
	     " while in trajectory it is (%g)\n",ir->init_lambda,lambda);
  }
  else {
    ir->nsteps     -= nstep;
    ir->init_t      = t;
    ir->init_lambda = lambda;
  }
  
  if (!ftp2bSet(efTRN,NFILE,fnm)) {
    get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,
	      &gnx,&index,&grpname);
    bSel=(gnx!=natoms);
    for (i=0; ((i<gnx) && (!bSel)); i++)
      bSel = (i!=index[i]);
    if (bSel) {
      fprintf(stderr,"Will write subset %s of original tpx containg %d "
	      "atoms\n",grpname,gnx);
      reduce_topology_x(gnx,index,&top,x,v);
      natoms = gnx;
    } else
      fprintf(stderr,"Will write full tpx file (no selection)\n");
  }    
  
  printf("writing statusfile with starting time %g...\n",ir->init_t);
  write_tpx(opt2fn("-o",NFILE,fnm),
	    0,ir->init_t,ir->init_lambda,ir,box,
	    natoms,x,v,NULL,&top);
  thanx(stdout);
  
  return 0;
}
