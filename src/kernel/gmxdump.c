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
#include "invblock.h"
#include "txtdump.h"
#include "gmxcpp.h"
#include "checkpoint.h"
#include "mtop_util.h"
#include "sparsematrix.h"
#include "mtxio.h"

char *ft2category(int ft)
{
  switch(ft) {
  case F_BONDS: 
  case F_G96BONDS:
  case F_MORSE:
  case F_CUBICBONDS:
  case F_CONNBONDS:
  case F_HARMONIC:
  case F_FENEBONDS:
  case F_TABBONDS:
  case F_TABBONDSNC:
  case F_RESTRBONDS:
    return "bonds";
  case F_ANGLES: 
  case F_G96ANGLES:
  case F_CROSS_BOND_BONDS:
  case F_CROSS_BOND_ANGLES:
  case F_UREY_BRADLEY:
  case F_QUARTIC_ANGLES:
  case F_TABANGLES:
    return "angles";
  case F_PDIHS:
  case F_RBDIHS: 
  case F_FOURDIHS:
  case F_IDIHS: 
  case F_PIDIHS: 
  case F_TABDIHS:
    return "dihedrals";
  case F_CMAP:
  case F_GB12:
  case F_GB13:
  case F_GB14:
  case F_GBPOL:
  case F_NPSOLVATION:
    return "unknown";
  case F_LJ14:
  case F_COUL14:
  case F_LJC14_Q:
  case F_LJC_PAIRS_NB:
    return "pairs";
  case F_LJ:
  case F_BHAM:
  case F_LJ_LR:
  case F_BHAM_LR:
  case F_DISPCORR:
  case F_COUL_SR:
  case F_COUL_LR:
  case F_RF_EXCL:
  case F_COUL_RECIP:
  case F_DPD:
    return "unknown";
  case F_POLARIZATION:
  case F_WATER_POL:
  case F_THOLE_POL:
    return "polarization";
  case F_POSRES:
    return "position_restraints";
  case F_DISRES:
    return "distance_restraints";
  case F_DISRESVIOL:
  case F_ORIRES:
    return "orientation_restraints";
  case F_ORIRESDEV:
    return "unknown";
  case F_ANGRES:
  case F_ANGRESZ:
    return "angle_restraints";
  case F_DIHRES:
    return "dihedral_restraints";
  case F_DIHRESVIOL:
    return "unknown";
  case F_CONSTR:
  case F_CONSTRNC:
    return "constraints";
  case F_SETTLE:
    return "settle";
  case F_VSITE2:
    return "virtual_sites2";
  case F_VSITE3:
  case F_VSITE3FD:
  case F_VSITE3FAD:
  case F_VSITE3OUT:
    return "virtual_sites3";
  case F_VSITE4FD:
  case F_VSITE4FDN:
    return "virtual_sites4";
  case F_VSITEN:
  case F_COM_PULL:
  case F_EQM:
  case F_EPOT:
  case F_EKIN:
  case F_ETOT:
  case F_ECONSERVED:
  case F_TEMP:
  case F_VTEMP:
  case F_PDISPCORR:
  case F_PRES:
  case F_DVDL:
  case F_DKDL:
  case F_DHDL_CON:
  default:
    return "poep";
  }
}

static void dump_ilist(FILE *fp,int ft,t_ilist *il)
{
  if (il->nr > 0) {
    fprintf(fp,"[ %s ]\n",ft2category(ft));
  }
  fprintf(fp,"\n");
 
}

static void dump_moltype(FILE *fp,gmx_moltype_t mt)
{
  atom_id *invcg;
  int i,j,resnr;
  
  fprintf(fp,"[ molecule_type ]\n");
  fprintf(fp,"%s 1\n\n",*mt.name);
  
  invcg = make_invblock(&mt.cgs,mt.atoms.nr);
  fprintf(fp,"[ atoms ]\n");
  fprintf(fp,";   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB\n");
  for(i=0; (i<mt.atoms.nr); i++) {
    resnr = mt.atoms.atom[i].resind;
    fprintf(fp,"%5d  %8s  %5d  %5s  %5s  %5d  %10g  %10g  %8s  %10g  %10g\n",
	    i+1,*mt.atoms.atomtype[i],1+resnr,*mt.atoms.resinfo[resnr].name,
	    *mt.atoms.atomname[i],invcg[i]+1,mt.atoms.atom[i].q,mt.atoms.atom[i].m,
	    *mt.atoms.atomtypeB[i],mt.atoms.atom[i].qB,mt.atoms.atom[i].mB);
  }
  fprintf(fp,"\n");
  for(j=0; (j<F_NRE); j++) 
    dump_ilist(fp,j,&mt.ilist[j]);
  if (mt.excls.nr > 0) {
    fprintf(fp,"[ exclusions ]\n");
    for(i=0; (i<mt.excls.nr); i++) {
      if (mt.excls.index[i] < mt.excls.index[i+1]) {
	fprintf(fp,"%d",i+1);
	for(j=mt.excls.index[i]; (j<mt.excls.index[i+1]); j++) {
	  fprintf(fp,"  %d",mt.excls.a[j]+1);
	}
	fprintf(fp,"\n");
      }
    }
  }
}

static void dump_atomtypes(FILE *fp,t_atomtypes at)
{
  int i;
  
  fprintf(fp,"[ atomtypes ]\n");
  fprintf(fp,";name  at.num    mass      charge ptype        c6        c12\n");
  
  
  
}

static void dump_top(const char *tpr,const char *top)
{
  FILE        *fp;
  int         indent,i,j,**gcount,atot;
  t_state     state;
  t_tpxheader tpx;
  gmx_mtop_t  mtop;

  read_tpxheader(tpr,&tpx,TRUE,NULL,NULL);
  
  if (!tpx.bTop) {
    printf("Need a tpr file containing a topology in order to write a top file\n");
    return;
  }
  read_tpx_state(tpr,NULL,&state,NULL,&mtop);
  fp = ffopen(top,"w");
  fprintf(fp,"; Topology generated from %s by program %s\n",tpr,Program());
  fprintf(fp,"[ defaults ]\n 1 1 no 1.0 1.0\n\n");
  dump_atomtypes(fp,mtop.atomtypes);
     
  for(j=0; (j<mtop.nmoltype); j++) {
    dump_moltype(fp,mtop.moltype[j]);
  }

  fprintf(fp,"[ system ]\n");
  fprintf(fp,"%s\n\n",*mtop.name);
  for(j=0; (j<mtop.nmolblock); j++) {
    fprintf(fp,"%s  %d\n",*(mtop.moltype[mtop.molblock[j].type].name),mtop.molblock[j].nmol);
  }
  fclose(fp);
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
    "Note that currently setting GMXLIB is the only way to customize",
    "directories used for searching include files.",
  };
  t_filenm fnm[] = {
    { efTPX, "-s", NULL, ffOPTRD },
    { efTRX, "-f", NULL, ffOPTRD },
    { efEDR, "-e", NULL, ffOPTRD },
    { efCPT, NULL, NULL, ffOPTRD },
    { efTOP, "-p", NULL, ffOPTRD },
    { efMTX, "-mtx", "hessian", ffOPTRD }, 
    { efMDP, "-om", NULL, ffOPTWR },
    { efTOP, "-op", NULL, ffOPTWR }
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
		    asize(desc),desc,0,NULL,&oenv);


  if (ftp2bSet(efTPX,NFILE,fnm)) {
    if (opt2bSet("-op",NFILE,fnm))
      dump_top(ftp2fn(efTPX,NFILE,fnm),opt2fn("-op",NFILE,fnm));
    else
      list_tpx(ftp2fn(efTPX,NFILE,fnm),bShowNumbers,
	       ftp2fn_null(efMDP,NFILE,fnm),bSysTop);
  }
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
