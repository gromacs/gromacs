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
    
    for(i=0; (i<MAXNODES); i++)
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
    "tpbconv can edit run input files in two ways.[PAR]"
    "[BB]1st.[bb] by creating a run input file",
    "for a continuation run when your simulation has crashed due to e.g.",
    "a full disk, or by making a continuation run input file.",
    "Note that a frame with coordinates and velocities is needed,",
    "which means that when you never write velocities, you can not use",
    "tpbconv and you have to start the run again from the beginning.[PAR]",
    "[BB]2nd.[bb] by creating a tpx file for a subset of your original",
    "tpx file, which is useful when you want to remove the solvent from",
    "your tpx file, or when you want to make e.g. a pure Ca tpx file.",
    "[BB]WARNING: this tpx file is not fully functional[bb]."
  };

  char         *top_fn,*frame_fn;
  int          fp;
  t_tpxheader  tpx;
  t_trnheader head;
  int          i,natoms,frame,step,run_step,nsteps_org;
  real         run_t,run_lambda;
  bool         bOK,bFrame,bTime,bSel;
  t_topology   top;
  t_inputrec   *ir,*irnew=NULL;
  t_gromppopts *gopts;
  rvec         *x=NULL,*v=NULL,*newx,*newv,*tmpx,*tmpv;
  matrix       box,newbox;
  int          gnx;
  char         *grpname;
  atom_id      *index=NULL;
  t_filenm fnm[] = {
    { efTPX, NULL,  NULL,    ffREAD  },
    { efTRN, "-f",  NULL,    ffOPTRD },
    { efNDX, NULL,  NULL,    ffOPTRD },
    { efTPX, "-o",  "tpxout",ffWRITE }
  };
#define NFILE asize(fnm)

  /* Command line options */
  static real start_t = -1.0, extend_t = 0.0, until_t = 0.0;
  static t_pargs pa[] = {
    { "-time", FALSE, etREAL, {&start_t}, 
      "Continue from frame at this time (ps) instead of the last frame" },
    { "-extend", FALSE, etREAL, {&extend_t}, 
      "Extend runtime by this amount (ps)" },
    { "-until", FALSE, etREAL, {&until_t}, 
      "Extend runtime until this ending time (ps)" },
  };
  int nerror = 0;
  
  CopyRight(stdout,argv[0]);
  
  /* Parse the command line */
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);

  bTime = opt2parg_bSet("-time",asize(pa),pa);
  
  top_fn = ftp2fn(efTPX,NFILE,fnm);
  fprintf(stderr,"Reading toplogy and shit from %s\n",top_fn);
  
  read_tpxheader(top_fn,&tpx);
  snew(x,tpx.natoms);
  snew(v,tpx.natoms);
  snew(ir,1);
  read_tpx(top_fn,&step,&run_t,&run_lambda,ir,box,&natoms,x,v,NULL,&top);
  run_step   = 0;

  if (ftp2bSet(efTRN,NFILE,fnm)) {
    frame_fn = ftp2fn(efTRN,NFILE,fnm);
    fprintf(stderr,
	    "\nREADING COORDS, VELS AND BOX FROM TRAJECTORY %s...\n\n",
	    frame_fn);
    
    fp=open_trn(frame_fn,"r");
    fread_trnheader(fp,&head,&bOK);
    if (top.atoms.nr != head.natoms) 
      fatal_error(0,"Number of atoms in Topology (%d) "
		  "is not the same as in Trajectory (%d)\n",
		  top.atoms.nr,head.natoms);
    snew(newx,head.natoms);
    snew(newv,head.natoms);
    run_t      = head.t;
    run_step   = head.step;
    run_lambda = head.lambda;
    bOK        = fread_htrn(fp,&head,box,x,v,NULL);

    /* Now scan until the last set of x and v (step == 0)
     * or the ones at step step.
     */
    bFrame = bOK;
    frame  = 0;
    while (bFrame) {
      fprintf(stderr,"\rRead frame %6d: step %6d time %8.3f",
	      frame,run_step,run_t);
      bFrame=fread_trnheader(fp,&head,&bOK);
      if (bFrame || !bOK)
	frame++;
      bFrame=bFrame && bOK;
      if (bFrame)
	bOK=fread_htrn(fp,&head,newbox,newx,newv,NULL);
      bFrame=bFrame && bOK;
      if (bFrame && (head.x_size) && (head.v_size)) {
	tmpx = newx;
	newx = x;
	x    = tmpx;
	tmpv = newv;
	newv = v;
	v    = tmpv;
	run_t      = head.t;
	run_step   = head.step;
	run_lambda = head.lambda;
	copy_mat(newbox,box);
      }
      if (bTime && (head.t >= start_t))
	bFrame=FALSE;
    }
    close_trn(fp);
    fprintf(stderr,"\n");
    if (!bOK)
      fprintf(stderr,"Frame %d (step %d, time %g) is incomplete\n",
	      frame,head.step,head.t);
    fprintf(stderr,"\nUsing frame of step %d time %g\n",run_step,run_t);
  } 
  else {
    frame_fn = ftp2fn(efTPX,NFILE,fnm);
    fprintf(stderr,"\nUSING COORDS, VELS AND BOX FROM TPX FILE %s...\n\n",
	    ftp2fn(efTPX,NFILE,fnm));
  }

  /* Determine total number of steps remaining */
  if (extend_t) {
    ir->nsteps = ir->nsteps - run_step + (int)(extend_t/ir->delta_t + 0.5);
    printf("Extending remaining runtime of by %g ps (now %d steps)\n",
	   extend_t,ir->nsteps);
  }
  else if (until_t) {
    printf("nsteps = %d, run_step = %d, current_t = %g, until = %g\n",
	   ir->nsteps,run_step,run_t,until_t);
    ir->nsteps = (int)((until_t-run_t)/ir->delta_t + 0.5);
    printf("Extending remaining runtime until %g ps (now %d steps)\n",
	   until_t,ir->nsteps);
  }
  else {
    ir->nsteps -= run_step; 
    /* Print message */
    printf("%d steps (%g ps) remaining from first run.\n",
	   ir->nsteps,ir->nsteps*ir->delta_t);
	  
  }
  if (ir->nsteps > 0) {

    ir->init_t      = run_t;
    ir->init_lambda = run_lambda;
    
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
      } 
      else
	fprintf(stderr,"Will write full tpx file (no selection)\n");
    }    
    
    fprintf(stderr,"Writing statusfile with starting time %g and %d steps...\n",
	    ir->init_t,ir->nsteps);
    write_tpx(opt2fn("-o",NFILE,fnm),
	      0,ir->init_t,ir->init_lambda,ir,box,
	      natoms,x,v,NULL,&top);
  }
  else
    printf("You've simulated long enough. Not writing tpr file\n");
	      
  thanx(stderr);
  
  return 0;
}
