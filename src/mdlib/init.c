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
static char *SRCID_init_c = "$Id$";

#include <stdio.h>
#include "typedefs.h"
#include "tpxio.h"
#include "smalloc.h"
#include "vec.h"
#include "vveclib.h"
#include "main.h"
#include "mvdata.h"
#include "fatal.h"
#include "symtab.h"
#include "txtdump.h"
#include "splittop.h"
#include "mdatoms.h"
#include "mdrun.h"
#include "statutil.h"

#define BUFSIZE	256

#define NOT_FINISHED(l1,l2) \
  printf("not finished yet: lines %d .. %d in %s\n",l1,l2,__FILE__)

void check_nnodes_top(char *fn,t_topology *top,int nnodes)
{
  int i,np=0;

  for(i=MAXNODES-1; (i>0) && (top->blocks[ebCGS].multinr[i] == 0); i--)
    ;
  np = i+1;
  
  if (np != nnodes)
    fatal_error(0,"run input file %s was made for %d nodes,\n"
		"             while %s expected it to be for %d nodes.",
		fn,np,ShortProgram(),nnodes);
}

static char *int_title(char *title,int nodeid)
{
  static char buf[BUFSIZE];

  sprintf(buf,"%s (%d)",title,nodeid);
  return buf;
}

static void rm_idef(t_idef *idef)
{
  int i;
  
  sfree(idef->functype);
  sfree(idef->iparams);
  for(i=0; (i<F_NRE); i++)
    sfree(idef->il[i].iatoms);
}

static void rm_block(t_block *block)
{
  sfree(block->index);
  sfree(block->a);
}

static void rm_atoms(t_atoms *atoms)
{
  sfree(atoms->atom);
  sfree(atoms->atomname);
  sfree(atoms->resname);
  sfree(atoms->grpname);
  rm_block(&atoms->excl);
}

void init_single(FILE *log,t_parm *parm,
		 char *tpxfile,t_topology *top, 
                 rvec **x,rvec **v,t_mdatoms **mdatoms,
		 t_nsborder *nsb)
{
  int         step,natoms;
  real        t,lambda;
  t_tpxheader *tpx;

  snew(tpx,1);
  
  read_tpxheader(tpxfile,tpx);
  snew(*x,tpx->natoms);
  snew(*v,tpx->natoms);
  
  sfree(tpx);
  
  read_tpx(tpxfile,&step,&t,&lambda,&parm->ir,
	   parm->box,&natoms,*x,*v,NULL,top);
  check_nnodes_top(tpxfile,top,1);

  *mdatoms=atoms2md(log,&top->atoms,parm->ir.opts.nFreeze,
		    parm->ir.eI==eiBD && parm->ir.bd_fric!=0,
		    parm->ir.efep!=efepNO,FALSE);
  
  pr_inputrec(log,0,"Input Parameters",&parm->ir);
  calc_nsb(log,&(top->blocks[ebCGS]),1,nsb,0);
  print_nsb(log,"Neighbor Search Blocks",nsb);
}

void distribute_parts(int left,int right,int nodeid,int nnodes,t_parm *parm,
		      char *tpxfile,int nstDlb)
{
  int         natoms,step;
  real        t,lambda;
  t_tpxheader tpx;
  t_topology  top;
  t_nsborder  nsb;
  rvec        *x,*v;
  
  read_tpxheader(tpxfile,&tpx);
  snew(x,tpx.natoms);
  snew(v,tpx.natoms);
  read_tpx(tpxfile,&step,&t,&lambda,&parm->ir,parm->box,
	   &natoms,x,v,NULL,&top);
  check_nnodes_top(tpxfile,&top,nnodes);
  
  calc_nsb(stdlog,&(top.blocks[ebCGS]),nnodes,&nsb,nstDlb);
  mv_data(left,right,parm,&nsb,&top,x,v);
  done_top(&top);
  sfree(x);
  sfree(v);
}

static void count_ones(FILE *log,t_commrec *cr)
{
#define NMAX 6
  int i;
  int n[NMAX];
  
  where();  
  fprintf(log,"Topology loaded succesfully\n"); fflush(log);
  for(i=0; (i<NMAX); i++)
    n[i]=i;
  gmx_sumi(6,n,cr);
  for(i=0; (i<NMAX); i++) {
    if (n[i] != (i*cr->nnodes)) {
      fprintf(log,"Error in count_ones: n[%d] = %d, should be %d\n",
	      i,n[i],i*cr->nnodes);
      fflush(log);
    }
  }
  fprintf(log,"Completed confidence test\n"); fflush(log);
  where();
#undef NMAX
}

void init_parts(FILE *log,t_commrec *cr,
		t_parm *parm,t_topology *top,
		rvec **x,rvec **v,t_mdatoms **mdatoms,
		t_nsborder *nsb,int list)
{
  ld_data(cr->left,cr->right,parm,nsb,top,x,v);
  if (cr->nodeid != 0)
    mv_data(cr->left,cr->right,parm,nsb,top,*x,*v);
#ifdef _amb_
  count_ones(log,cr);
#endif

  mdsplit_top(log,top,cr);

  if (list) {
    if (list&LIST_SCALARS) 
      print_nsb(log,"Listing Scalars",nsb);
    if (list&LIST_PARM)
      write_parm(log,"parameters of the run",cr->nodeid,parm);
    if (list&LIST_X)
      pr_rvecs(log,0,int_title("x",0),*x,nsb->natoms);
    if (list&LIST_V)
      pr_rvecs(log,0,int_title("v",0),*v,nsb->natoms);
    if (list&LIST_TOP)
      pr_top(log,0,int_title("topology",cr->nodeid),top);
    fflush(log);
  }
  *mdatoms=atoms2md(log,&(top->atoms),parm->ir.opts.nFreeze,
		    parm->ir.eI==eiBD && parm->ir.bd_fric!=0,
		    parm->ir.efep!=efepNO,FALSE);
}

void write_parm(FILE *log,char *title,int nodeid,t_parm *parm)
{
  fprintf(log,"%s (nodeid=%d):\n",title,nodeid);
  pr_inputrec(log,0,"input record",&parm->ir);
  pr_rvecs(log,0,"box",parm->box,DIM);
  pr_rvecs(log,0,"ekin",parm->ekin,DIM);
  pr_rvecs(log,0,"pres",parm->pres,DIM);
  pr_rvecs(log,0,"vir",parm->vir,DIM);
}

