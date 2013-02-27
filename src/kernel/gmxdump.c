/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "main.h"
#include "macros.h"
#include "futil.h"
#include "statutil.h"
#include "copyrite.h"
#include "sysstuff.h"
#include "txtdump.h"
#include "gmx_fatal.h"
#include "xtcio.h"
#include "enxio.h"
#include "smalloc.h"
#include "names.h"
#include "gmxfio.h"
#include "tpxio.h"
#include "trnio.h"
#include "txtdump.h"
#include "gmxcpp.h"
#include "checkpoint.h"
#include "mtop_util.h"
#include "sparsematrix.h"
#include "mtxio.h"

static void dump_top(FILE *fp,t_topology *top,char *tpr)
{
  int i,j,k,*types;
  
  fprintf(fp,"; Topology generated from %s by program %s\n",tpr,Program());
  fprintf(fp,"[ defaults ]\n 1 1 no 1.0 1.0\n\n");
  fprintf(fp,"[ atomtypes ]\n");
  fprintf(fp,";name  at.num    mass      charge ptype        c6        c12\n");
  snew(types,top->atomtypes.nr);
  for(i=0; (i<top->atomtypes.nr); i++) {
    for(j=0; (j<top->atoms.nr) && (top->atoms.atom[j].type != i); j++)
      ;
    if (j<top->atoms.nr) {
      types[i] = j;
      fprintf(fp,"%5s  %4d   %8.4f   %8.4f  %2s  %8.3f  %8.3f\n",
	      *(top->atoms.atomtype[j]),top->atomtypes.atomnumber[i],
	      0.0,0.0,"A",0.0,0.0);
    }
  }
  fprintf(fp,"[ nonbonded_params ]\n");
  for(i=k=0; (i<top->idef.ntypes); i++) {
    for(j=0; (j<top->idef.ntypes); j++,k++) {
      fprintf(fp,"%12s  %12s  1  %12.5e  %12.5e\n",
	      *(top->atoms.atomtype[types[i]]),
	      *(top->atoms.atomtype[types[j]]),
	      top->idef.iparams[k].lj.c12,top->idef.iparams[k].lj.c6);
    }
  }
  sfree(types);
}

static void list_tpx(const char *fn, gmx_bool bShowNumbers,const char *mdpfn,
                     gmx_bool bSysTop)
{
  FILE *gp;
  int         fp,indent,i,j,**gcount,atot;
  t_state     state;
  rvec        *f=NULL;
  t_inputrec  ir;
  t_tpxheader tpx;
  gmx_mtop_t  mtop;
  gmx_groups_t *groups;
  t_topology  top;

  read_tpxheader(fn,&tpx,TRUE,NULL,NULL);

  read_tpx_state(fn,
		 tpx.bIr  ? &ir : NULL,
		 &state,tpx.bF ? f : NULL,
		 tpx.bTop ? &mtop: NULL);
  
  if (mdpfn && tpx.bIr) {
    gp = gmx_fio_fopen(mdpfn,"w");
    pr_inputrec(gp,0,NULL,&(ir),TRUE);
    gmx_fio_fclose(gp);
  }

  if (!mdpfn) {  
    if (bSysTop)
      top = gmx_mtop_t_to_t_topology(&mtop);

    if (available(stdout,&tpx,0,fn)) {
      indent=0;
      indent=pr_title(stdout,indent,fn);
      pr_inputrec(stdout,0,"inputrec",tpx.bIr ? &(ir) : NULL,FALSE);
      
      indent = 0;
      pr_header(stdout,indent,"header",&(tpx));
      
      if (!bSysTop)
	pr_mtop(stdout,indent,"topology",&(mtop),bShowNumbers);
      else
	pr_top(stdout,indent,"topology",&(top),bShowNumbers);

      pr_rvecs(stdout,indent,"box",tpx.bBox ? state.box : NULL,DIM);
      pr_rvecs(stdout,indent,"box_rel",tpx.bBox ? state.box_rel : NULL,DIM);
      pr_rvecs(stdout,indent,"boxv",tpx.bBox ? state.boxv : NULL,DIM);
      pr_rvecs(stdout,indent,"pres_prev",tpx.bBox ? state.pres_prev : NULL,DIM);
      pr_rvecs(stdout,indent,"svir_prev",tpx.bBox ? state.svir_prev : NULL,DIM);
      pr_rvecs(stdout,indent,"fvir_prev",tpx.bBox ? state.fvir_prev : NULL,DIM);
      /* leave nosehoover_xi in for now to match the tpr version */
      pr_doubles(stdout,indent,"nosehoover_xi",state.nosehoover_xi,state.ngtc);
      /*pr_doubles(stdout,indent,"nosehoover_vxi",state.nosehoover_vxi,state.ngtc);*/
      /*pr_doubles(stdout,indent,"therm_integral",state.therm_integral,state.ngtc);*/
      pr_rvecs(stdout,indent,"x",tpx.bX ? state.x : NULL,state.natoms);
      pr_rvecs(stdout,indent,"v",tpx.bV ? state.v : NULL,state.natoms);
      if (state,tpx.bF) {
	pr_rvecs(stdout,indent,"f",f,state.natoms);
      }
    }
    
    groups = &mtop.groups;

    snew(gcount,egcNR);
    for(i=0; (i<egcNR); i++) 
      snew(gcount[i],groups->grps[i].nr);
    
    for(i=0; (i<mtop.natoms); i++) {
      for(j=0; (j<egcNR); j++) 
	gcount[j][ggrpnr(groups,j,i)]++;
    }
    printf("Group statistics\n");
    for(i=0; (i<egcNR); i++) {
      atot=0;
      printf("%-12s: ",gtypes[i]);
      for(j=0; (j<groups->grps[i].nr); j++) {
	printf("  %5d",gcount[i][j]);
	atot+=gcount[i][j];
      }
      printf("  (total %d atoms)\n",atot);
      sfree(gcount[i]);
    }
    sfree(gcount);
  }
  done_state(&state);
  sfree(f);
}

static void list_top(const char *fn)
{
  int status,done;
#define BUFLEN 256
  char buf[BUFLEN];
  gmx_cpp_t handle;
  char *cppopts[] = { NULL };

  status = cpp_open_file(fn,&handle,cppopts);
  if (status != 0) 
    gmx_fatal(FARGS,cpp_error(&handle,status));
  do {
    status = cpp_read_line(&handle,BUFLEN,buf);
    done = (status == eCPP_EOF);
    if (!done) {
      if (status != eCPP_OK)
	gmx_fatal(FARGS,cpp_error(&handle,status));
      else 
	printf("%s\n",buf);
    }
  } while (!done);
  status = cpp_close_file(&handle);
  if (status != eCPP_OK) 
    gmx_fatal(FARGS,cpp_error(&handle,status));
}

static void list_trn(const char *fn)
{
  t_fileio    *fpread, *fpwrite;
  int         nframe,indent;
  char        buf[256];
  rvec        *x,*v,*f;
  matrix      box;
  t_trnheader trn;
  gmx_bool        bOK;

  fpread  = open_trn(fn,"r"); 
  fpwrite = open_tpx(NULL,"w");
  gmx_fio_setdebug(fpwrite,TRUE);
  
  nframe = 0;
  while (fread_trnheader(fpread,&trn,&bOK)) {
    snew(x,trn.natoms);
    snew(v,trn.natoms);
    snew(f,trn.natoms);
    if (fread_htrn(fpread,&trn,
		   trn.box_size ? box : NULL,
		   trn.x_size   ? x : NULL,
		   trn.v_size   ? v : NULL,
		   trn.f_size   ? f : NULL)) {
      sprintf(buf,"%s frame %d",fn,nframe);
      indent=0;
      indent=pr_title(stdout,indent,buf);
      pr_indent(stdout,indent);
      fprintf(stdout,"natoms=%10d  step=%10d  time=%12.7e  lambda=%10g\n",
	      trn.natoms,trn.step,trn.t,trn.lambda);
      if (trn.box_size)
	pr_rvecs(stdout,indent,"box",box,DIM);
      if (trn.x_size)
	pr_rvecs(stdout,indent,"x",x,trn.natoms);
      if (trn.v_size)
	pr_rvecs(stdout,indent,"v",v,trn.natoms);
      if (trn.f_size)
	pr_rvecs(stdout,indent,"f",f,trn.natoms);
    } 
    else
      fprintf(stderr,"\nWARNING: Incomplete frame: nr %d, t=%g\n",
	      nframe,trn.t);
    
    sfree(x);
    sfree(v);
    sfree(f);
    nframe++;
  }
  if (!bOK)
    fprintf(stderr,"\nWARNING: Incomplete frame header: nr %d, t=%g\n",
	    nframe,trn.t);
  close_tpx(fpwrite);
  close_trn(fpread);
}

void list_xtc(const char *fn, gmx_bool bXVG)
{
  t_fileio *xd;
  int    indent;
  char   buf[256];
  rvec   *x;
  matrix box;
  int    nframe,natoms,step;
  real   prec,time;
  gmx_bool   bOK;
  
  xd = open_xtc(fn,"r");
  read_first_xtc(xd,&natoms,&step,&time,box,&x,&prec,&bOK);
		
  nframe=0;
  do {
    if (bXVG) {
      int i,d;
      
      fprintf(stdout,"%g",time);
      for(i=0; i<natoms; i++)
	for(d=0; d<DIM; d++)
	  fprintf(stdout," %g",x[i][d]);
      fprintf(stdout,"\n");
    } else {
      sprintf(buf,"%s frame %d",fn,nframe);
      indent=0;
      indent=pr_title(stdout,indent,buf);
      pr_indent(stdout,indent);
      fprintf(stdout,"natoms=%10d  step=%10d  time=%12.7e  prec=%10g\n",
	    natoms,step,time,prec);
      pr_rvecs(stdout,indent,"box",box,DIM);
      pr_rvecs(stdout,indent,"x",x,natoms);
    }
    nframe++;
  } while (read_next_xtc(xd,natoms,&step,&time,box,x,&prec,&bOK));
  if (!bOK)
    fprintf(stderr,"\nWARNING: Incomplete frame at time %g\n",time);
  close_xtc(xd);
}

void list_trx(const char *fn,gmx_bool bXVG)
{
  int ftp;
  
  ftp = fn2ftp(fn);
  if (ftp == efXTC)
    list_xtc(fn,bXVG);
  else if ((ftp == efTRR) || (ftp == efTRJ))
    list_trn(fn);
  else
    fprintf(stderr,"File %s is of an unsupported type. Try using the command\n 'less %s'\n",
	    fn,fn);
}

void list_ene(const char *fn)
{
    int        ndr;
    ener_file_t in;
    gmx_bool       bCont;
    gmx_enxnm_t *enm=NULL;
    t_enxframe *fr;
    int        i,j,nre,b;
    real       rav,minthird;
    char       buf[22];

    printf("gmxdump: %s\n",fn);
    in = open_enx(fn,"r");
    do_enxnms(in,&nre,&enm);

    printf("energy components:\n");
    for(i=0; (i<nre); i++) 
        printf("%5d  %-24s (%s)\n",i,enm[i].name,enm[i].unit);

    minthird=-1.0/3.0;
    snew(fr,1);
    do {
        bCont=do_enx(in,fr);

        if (bCont) 
        {
            printf("\n%24s  %12.5e  %12s  %12s\n","time:",
                   fr->t,"step:",gmx_step_str(fr->step,buf));
            printf("%24s  %12s  %12s  %12s\n",
                   "","","nsteps:",gmx_step_str(fr->nsteps,buf));
            printf("%24s  %12.5e  %12s  %12s\n",
                   "delta_t:",fr->dt,"sum steps:",gmx_step_str(fr->nsum,buf));
            if (fr->nre == nre) {
                printf("%24s  %12s  %12s  %12s\n",
                       "Component","Energy","Av. Energy","Sum Energy");
                if (fr->nsum > 0) {
                    for(i=0; (i<nre); i++) 
                        printf("%24s  %12.5e  %12.5e  %12.5e\n",
                               enm[i].name,fr->ener[i].e,fr->ener[i].eav,
                               fr->ener[i].esum);
                } else {
                    for(i=0; (i<nre); i++) 
                        printf("%24s  %12.5e\n",
                               enm[i].name,fr->ener[i].e);
                }
            }
            for(b=0; b<fr->nblock; b++)
            {
                const char *typestr="";

                t_enxblock *eb=&(fr->block[b]);
                printf("Block data %2d (%3d subblocks, id=%d)\n",
		       b, eb->nsub, eb->id);

                if (eb->id < enxNR)
                    typestr=enx_block_id_name[eb->id];
                printf("  id='%s'\n", typestr);
                for(i=0;i<eb->nsub;i++)
                {
                    t_enxsubblock *sb=&(eb->sub[i]);
                    printf("  Sub block %3d (%5d elems, type=%s) values:\n", 
                           i, sb->nr, xdr_datatype_names[sb->type]);

                    switch(sb->type)
                    {
                        case xdr_datatype_float:
                            for(j=0;j<sb->nr;j++)
                                printf("%14d   %8.4f\n",j, sb->fval[j]);
                            break;
                        case xdr_datatype_double:
                            for(j=0;j<sb->nr;j++)
                                printf("%14d   %10.6f\n",j, sb->dval[j]);
                            break;
                        case xdr_datatype_int:
                            for(j=0;j<sb->nr;j++)
                                printf("%14d %10d\n",j, sb->ival[j]);
                            break;
                        case xdr_datatype_large_int:
                            for(j=0;j<sb->nr;j++)
                                printf("%14d %s\n",
				       j, gmx_step_str(sb->lval[j],buf));
                            break;
                        case xdr_datatype_char:
                            for(j=0;j<sb->nr;j++)
                                printf("%14d %1c\n",j, sb->cval[j]);
                            break;
                        case xdr_datatype_string:
                            for(j=0;j<sb->nr;j++)
                                printf("%14d %80s\n",j, sb->sval[j]);
                            break;
                        default:
                            gmx_incons("Unknown subblock type");
                    }
                }
            }
        }
    } while (bCont);

    close_enx(in);

    free_enxframe(fr);
    sfree(fr);
    sfree(enm);
}

static void list_mtx(const char *fn)
{
  int  nrow,ncol,i,j,k;
  real *full=NULL,value;
  gmx_sparsematrix_t * sparse=NULL;

  gmx_mtxio_read(fn,&nrow,&ncol,&full,&sparse);

  if (full == NULL) {
    snew(full,nrow*ncol);
    for(i=0;i<nrow*ncol;i++) {
      full[i] = 0;
    }
    
    for(i=0;i<sparse->nrow;i++) {
	for(j=0;j<sparse->ndata[i];j++) {
	  k     = sparse->data[i][j].col;
	  value = sparse->data[i][j].value;
	  full[i*ncol+k] = value;
	  full[k*ncol+i] = value;
	}
    }
    gmx_sparsematrix_destroy(sparse);
  }

  printf("%d %d\n",nrow,ncol);
  for(i=0; i<nrow; i++) {
    for(j=0; j<ncol; j++) {
      printf(" %g",full[i*ncol+j]);
    }
    printf("\n");
  }

  sfree(full);
}

int main(int argc,char *argv[])
{
  const char *desc[] = {
    "[TT]gmxdump[tt] reads a run input file ([TT].tpa[tt]/[TT].tpr[tt]/[TT].tpb[tt]),",
    "a trajectory ([TT].trj[tt]/[TT].trr[tt]/[TT].xtc[tt]), an energy",
    "file ([TT].ene[tt]/[TT].edr[tt]), or a checkpoint file ([TT].cpt[tt])",
    "and prints that to standard output in a readable format.",
    "This program is essential for checking your run input file in case of",
    "problems.[PAR]",
    "The program can also preprocess a topology to help finding problems.",
    "Note that currently setting [TT]GMXLIB[tt] is the only way to customize",
    "directories used for searching include files.",
  };
  const char *bugs[] = {
    "Position restraint output from -sys -s is broken"
  };
  t_filenm fnm[] = {
    { efTPX, "-s", NULL, ffOPTRD },
    { efTRX, "-f", NULL, ffOPTRD },
    { efEDR, "-e", NULL, ffOPTRD },
    { efCPT, NULL, NULL, ffOPTRD },
    { efTOP, "-p", NULL, ffOPTRD },
    { efMTX, "-mtx", "hessian", ffOPTRD }, 
    { efMDP, "-om", NULL, ffOPTWR }
  };
#define NFILE asize(fnm)

  output_env_t oenv;
  /* Command line options */
  static gmx_bool bXVG=FALSE;
  static gmx_bool bShowNumbers=TRUE;
  static gmx_bool bSysTop=FALSE;
  t_pargs pa[] = {
    { "-xvg", FALSE, etBOOL, {&bXVG}, "HIDDENXVG layout for xtc" },
    { "-nr",FALSE, etBOOL, {&bShowNumbers},"Show index numbers in output (leaving them out makes comparison easier, but creates a useless topology)" },
    { "-sys", FALSE, etBOOL, {&bSysTop}, "List the atoms and bonded interactions for the whole system instead of for each molecule type" }
  };
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,asize(pa),pa,
                    asize(desc),desc,asize(bugs),bugs,&oenv);


  if (ftp2bSet(efTPX,NFILE,fnm))
    list_tpx(ftp2fn(efTPX,NFILE,fnm),bShowNumbers,
	     ftp2fn_null(efMDP,NFILE,fnm),bSysTop);
  else if (ftp2bSet(efTRX,NFILE,fnm)) 
    list_trx(ftp2fn(efTRX,NFILE,fnm),bXVG);
  else if (ftp2bSet(efEDR,NFILE,fnm))
    list_ene(ftp2fn(efEDR,NFILE,fnm));
  else if (ftp2bSet(efCPT,NFILE,fnm))
    list_checkpoint(ftp2fn(efCPT,NFILE,fnm),stdout);
  else if (ftp2bSet(efTOP,NFILE,fnm))
    list_top(ftp2fn(efTOP,NFILE,fnm));
  else if (ftp2bSet(efMTX,NFILE,fnm))
    list_mtx(ftp2fn(efMTX,NFILE,fnm));
    
  thanx(stderr);

  return 0;
}
