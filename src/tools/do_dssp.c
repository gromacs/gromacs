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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "sysstuff.h"
#include "typedefs.h"
#include "string2.h"
#include "strdb.h"
#include "macros.h"
#include "smalloc.h"
#include "mshift.h"
#include "statutil.h"
#include "copyrite.h"
#include "pdbio.h"
#include "gmx_fatal.h"
#include "xvgr.h"
#include "matio.h"
#include "index.h"
#include "gstat.h"
#include "tpxio.h"
#include "viewit.h"

static int strip_dssp(char *dsspfile,int nres,
		       gmx_bool bPhobres[],real t,
		       real *acc,FILE *fTArea,
		       t_matrix *mat,int average_area[],
                       const output_env_t oenv)
{
  static gmx_bool bFirst=TRUE;
  static char *ssbuf;
  FILE *tapeout;
  static int xsize,frame;
  char buf[STRLEN+1];
  char SSTP;
  int  i,nr,iacc,nresidues;
  int  naccf,naccb; /* Count hydrophobic and hydrophilic residues */
  real iaccf,iaccb;
  t_xpmelmt c;
  
  tapeout=ffopen(dsspfile,"r");
  
  /* Skip header */
  do {
    fgets2(buf,STRLEN,tapeout);
  } while (strstr(buf,"KAPPA") == NULL);
  if (bFirst) {
    /* Since we also have empty lines in the dssp output (temp) file,
     * and some line content is saved to the ssbuf variable,
     * we need more memory than just nres elements. To be shure,
     * we allocate 2*nres-1, since for each chain there is a 
     * separating line in the temp file. (At most each residue 
     * could have been defined as a separate chain.) */
    snew(ssbuf,2*nres-1);
  }
  
  iaccb=iaccf=0;
  nresidues = 0;
  naccf = 0;
  naccb = 0;
  for(nr=0; (fgets2(buf,STRLEN,tapeout) != NULL); nr++) {
      if (buf[13] == '!')   /* Chain separator line has '!' at pos. 13 */
          SSTP='=';         /* Chain separator sign '=' */
      else
          SSTP=buf[16]==' ' ? '~' : buf[16];
    ssbuf[nr]=SSTP;
    
    buf[39]='\0';
    
    /* Only calculate solvent accessible area if needed */
    if ((NULL != acc) && (buf[13] != '!'))
    {
      sscanf(&(buf[34]),"%d",&iacc);
      acc[nr]=iacc;
      /* average_area and bPhobres are counted from 0...nres-1 */
      average_area[nresidues]+=iacc;
      if (bPhobres[nresidues])
      {
        naccb++;
        iaccb+=iacc;
      }
      else
      {
        naccf++;
        iaccf+=iacc;
      }
      /* Keep track of the residue number (do not count chain separator lines '!') */
      nresidues++;
    }
  }
  ssbuf[nr]='\0';
  
  if (bFirst) {
    if (0 != acc)
        fprintf(stderr, "%d residues were classified as hydrophobic and %d as hydrophilic.\n", naccb, naccf);
    
    sprintf(mat->title,"Secondary structure");
    mat->legend[0]=0;
    sprintf(mat->label_x,"%s",output_env_get_time_label(oenv));
    sprintf(mat->label_y,"Residue");
    mat->bDiscrete=TRUE;
    mat->ny=nr;
    snew(mat->axis_y,nr);
    for(i=0; i<nr; i++)
      mat->axis_y[i]=i+1;
    mat->axis_x=NULL;
    mat->matrix=NULL;
    xsize=0;
    frame=0;
    bFirst=FALSE;
  }
  if (frame>=xsize) {
    xsize+=10;
    srenew(mat->axis_x,xsize);
    srenew(mat->matrix,xsize);
  }
  mat->axis_x[frame]=t;
  snew(mat->matrix[frame],nr);
  c.c2=0;
  for(i=0; i<nr; i++) {
    c.c1=ssbuf[i];
    mat->matrix[frame][i] = max(0,searchcmap(mat->nmap,mat->map,c));
  }
  frame++;
  mat->nx=frame;
  
  if (fTArea)
    fprintf(fTArea,"%10g  %10g  %10g\n",t,0.01*iaccb,0.01*iaccf);
  ffclose(tapeout);
  
  /* Return the number of lines found in the dssp file (i.e. number
   * of redidues plus chain separator lines).
   * This is the number of y elements needed for the area xpm file */
  return nr;
}

gmx_bool *bPhobics(t_atoms *atoms)
{
  int  i,nb;
  char **cb;
  gmx_bool *bb;
  
  nb=get_strings("phbres.dat",&cb);
  snew(bb,atoms->nres);
  
  for(i=0; (i<atoms->nres); i++) {
    if (search_str(nb,cb,*atoms->resinfo[i].name) != -1)
      bb[i]=TRUE;
  }
  return bb;
}
 
static void check_oo(t_atoms *atoms)
{
  char *OOO;

  int i;

  OOO=strdup("O");

  for(i=0; (i<atoms->nr); i++) {
    if (strcmp(*(atoms->atomname[i]),"OXT")==0)
       *atoms->atomname[i]=OOO;
    else if (strcmp(*(atoms->atomname[i]),"O1")==0)
      *atoms->atomname[i]=OOO;
    else if (strcmp(*(atoms->atomname[i]),"OC1")==0)
      *atoms->atomname[i]=OOO;
  }
}

static void norm_acc(t_atoms *atoms, int nres, 
		     real av_area[], real norm_av_area[])
{
  int i,n,n_surf;
  
  char surffn[]="surface.dat";
  char **surf_res, **surf_lines;
  double *surf;
  
  n_surf = get_lines(surffn, &surf_lines);
  snew(surf, n_surf);
  snew(surf_res, n_surf);
  for(i=0; (i<n_surf); i++) {
    snew(surf_res[i], 5);
    sscanf(surf_lines[i],"%s %lf",surf_res[i],&surf[i]);
  }
  
  for(i=0; (i<nres); i++) {
    n = search_str(n_surf,surf_res,*atoms->resinfo[i].name);
    if ( n != -1)
      norm_av_area[i] = av_area[i] / surf[n];
    else 
      fprintf(stderr,"Residue %s not found in surface database (%s)\n",
	      *atoms->resinfo[i].name,surffn);
  }
}

void prune_ss_legend(t_matrix *mat)
{
  gmx_bool *present;
  int  *newnum;
  int i,r,f,newnmap;
  t_mapping *newmap;
  
  snew(present,mat->nmap);
  snew(newnum,mat->nmap);

  for(f=0; f<mat->nx; f++)
    for(r=0; r<mat->ny; r++)
      present[mat->matrix[f][r]]=TRUE;
      
  newnmap=0;
  newmap=NULL;
  for (i=0; i<mat->nmap; i++) {
    newnum[i]=-1;
    if (present[i]) {
      newnum[i]=newnmap;
      newnmap++;
      srenew(newmap,newnmap);
      newmap[newnmap-1]=mat->map[i];
    }
  }
  if (newnmap!=mat->nmap) {
    mat->nmap=newnmap;
    mat->map=newmap;
    for(f=0; f<mat->nx; f++)
      for(r=0; r<mat->ny; r++)
	mat->matrix[f][r]=newnum[mat->matrix[f][r]];
  }
}

void write_sas_mat(const char *fn,real **accr,int nframe,int nres,t_matrix *mat)
{
  real lo,hi;
  int i,j,nlev;
  t_rgb rlo={1,1,1}, rhi={0,0,0};
  FILE *fp;
  
  if(fn) {
    hi=lo=accr[0][0];
    for(i=0; i<nframe; i++)
      for(j=0; j<nres; j++) {
	lo=min(lo,accr[i][j]);
	hi=max(hi,accr[i][j]);
      }
    fp=ffopen(fn,"w");
    nlev=hi-lo+1;
    write_xpm(fp,0,"Solvent Accessible Surface","Surface (A^2)",
	      "Time","Residue Index",nframe,nres,
	      mat->axis_x,mat->axis_y,accr,lo,hi,rlo,rhi,&nlev);
    ffclose(fp);
  }
}

void analyse_ss(const char *outfile, t_matrix *mat, const char *ss_string,
                const output_env_t oenv)
{
  FILE *fp;
  t_mapping *map;
  int f,r,*count,ss_count;
  size_t s;
  const char** leg;
  
  map=mat->map;
  snew(count,mat->nmap);
  snew(leg,mat->nmap+1);
  leg[0]="Structure";
  for(s=0; s<mat->nmap; s++)
    leg[s+1]=strdup(map[s].desc);
  
  fp=xvgropen(outfile,"Secondary Structure",
	      output_env_get_xvgr_tlabel(oenv),"Number of Residues",oenv);
  if (output_env_get_print_xvgr_codes(oenv))
    fprintf(fp,"@ subtitle \"Structure = ");
  for(s=0; s<strlen(ss_string); s++) {
    if (s>0)
      fprintf(fp," + ");
    for(f=0; f<mat->nmap; f++)
      if (ss_string[s]==map[f].code.c1)
	fprintf(fp,"%s",map[f].desc);
  }
  fprintf(fp,"\"\n");
  xvgr_legend(fp,mat->nmap+1,leg,oenv);
  
  for(f=0; f<mat->nx; f++) {
    ss_count=0;
    for(s=0; s<mat->nmap; s++)
      count[s]=0;
    for(r=0; r<mat->ny; r++)
      count[mat->matrix[f][r]]++;
    for(s=0; s<mat->nmap; s++) {
      if (strchr(ss_string,map[s].code.c1))
	ss_count+=count[s];
    }
    fprintf(fp,"%8g %5d",mat->axis_x[f],ss_count);
    for(s=0; s<mat->nmap; s++)
      fprintf(fp," %5d",count[s]);
    fprintf(fp,"\n");
  }
  
  ffclose(fp);
  sfree(leg);
  sfree(count);
}

int main(int argc,char *argv[])
{
  const char *desc[] = {
    "do_dssp ", 
    "reads a trajectory file and computes the secondary structure for",
    "each time frame ",
    "calling the dssp program. If you do not have the dssp program,",
    "get it. do_dssp assumes that the dssp executable is",
    "/usr/local/bin/dssp. If this is not the case, then you should",
    "set an environment variable [BB]DSSP[bb] pointing to the dssp",
    "executable, e.g.: [PAR]",
    "[TT]setenv DSSP /opt/dssp/bin/dssp[tt][PAR]",
    "The structure assignment for each residue and time is written to an",
    "[TT].xpm[tt] matrix file. This file can be visualized with for instance",
    "[TT]xv[tt] and can be converted to postscript with [TT]xpm2ps[tt].",
    "Individual chains are separated by light grey lines in the xpm and",
    "postscript files.",
    "The number of residues with each secondary structure type and the",
    "total secondary structure ([TT]-sss[tt]) count as a function of",
    "time are also written to file ([TT]-sc[tt]).[PAR]",
    "Solvent accessible surface (SAS) per residue can be calculated, both in",
    "absolute values (A^2) and in fractions of the maximal accessible",
    "surface of a residue. The maximal accessible surface is defined as",
    "the accessible surface of a residue in a chain of glycines.",
    "[BB]Note[bb] that the program [TT]g_sas[tt] can also compute SAS",
    "and that is more efficient.[PAR]",
    "Finally, this program can dump the secondary structure in a special file",
    "[TT]ssdump.dat[tt] for usage in the program [TT]g_chi[tt]. Together",
    "these two programs can be used to analyze dihedral properties as a",
    "function of secondary structure type."
  };
  static gmx_bool bVerbose;
  static const char *ss_string="HEBT"; 
  t_pargs pa[] = {
    { "-v",  FALSE, etBOOL, {&bVerbose},
      "HIDDENGenerate miles of useless information" },
    { "-sss", FALSE, etSTR, {&ss_string},
      "Secondary structures for structure count"}
  };
  
  t_trxstatus *status;
  FILE       *tapein;
  FILE       *ss,*acc,*fTArea,*tmpf;
  const char *fnSCount,*fnArea,*fnTArea,*fnAArea;
  const char *leg[] = { "Phobic", "Phylic" };
  t_topology top;
  int        ePBC;
  t_atoms    *atoms;
  t_matrix   mat;
  int        nres,nr0,naccr,nres_plus_separators;
  gmx_bool       *bPhbres,bDoAccSurf;
  real       t;
  int        i,natoms,nframe=0;
  matrix     box;
  int        gnx;
  char       *grpnm,*ss_str;
  atom_id    *index;
  rvec       *xp,*x;
  int        *average_area;
  real       **accr,*accr_ptr=NULL,*av_area, *norm_av_area;
  char       pdbfile[32],tmpfile[32],title[256];
  char       dssp[256];
  const char *dptr;
  output_env_t oenv;
  gmx_rmpbc_t  gpbc=NULL;
  
  t_filenm   fnm[] = {
    { efTRX, "-f",   NULL,      ffREAD },
    { efTPS, NULL,   NULL,      ffREAD },
    { efNDX, NULL,   NULL,      ffOPTRD },
    { efDAT, "-ssdump", "ssdump", ffOPTWR },
    { efMAP, "-map", "ss",      ffLIBRD },
    { efXPM, "-o",   "ss",      ffWRITE },
    { efXVG, "-sc",  "scount",  ffWRITE },
    { efXPM, "-a",   "area",    ffOPTWR },
    { efXVG, "-ta",  "totarea", ffOPTWR },
    { efXVG, "-aa",  "averarea",ffOPTWR }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,
                    PCA_CAN_TIME | PCA_CAN_VIEW | PCA_TIME_UNIT | PCA_BE_NICE ,
		    NFILE,fnm, asize(pa),pa, asize(desc),desc,0,NULL,&oenv);
  fnSCount= opt2fn("-sc",NFILE,fnm);
  fnArea  = opt2fn_null("-a", NFILE,fnm);
  fnTArea = opt2fn_null("-ta",NFILE,fnm);
  fnAArea = opt2fn_null("-aa",NFILE,fnm);
  bDoAccSurf=(fnArea || fnTArea || fnAArea);
  
  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xp,NULL,box,FALSE);
  atoms=&(top.atoms);
  check_oo(atoms);
  bPhbres=bPhobics(atoms);
  
  get_index(atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&gnx,&index,&grpnm);
  nres=0;
  nr0=-1;
  for(i=0; (i<gnx); i++) {
    if (atoms->atom[index[i]].resind != nr0) {
      nr0=atoms->atom[index[i]].resind;
      nres++;
    }
  }
  fprintf(stderr,"There are %d residues in your selected group\n",nres);

  strcpy(pdbfile,"ddXXXXXX");
  gmx_tmpnam(pdbfile);
  if ((tmpf = fopen(pdbfile,"w")) == NULL) {
    sprintf(pdbfile,"%ctmp%cfilterXXXXXX",DIR_SEPARATOR,DIR_SEPARATOR);
    gmx_tmpnam(pdbfile);
    if ((tmpf = fopen(pdbfile,"w")) == NULL) 
      gmx_fatal(FARGS,"Can not open tmp file %s",pdbfile);
  }
  else
    fclose(tmpf);
    
  strcpy(tmpfile,"ddXXXXXX");
  gmx_tmpnam(tmpfile);
  if ((tmpf = fopen(tmpfile,"w")) == NULL) {
    sprintf(tmpfile,"%ctmp%cfilterXXXXXX",DIR_SEPARATOR,DIR_SEPARATOR);
    gmx_tmpnam(tmpfile);
    if ((tmpf = fopen(tmpfile,"w")) == NULL) 
      gmx_fatal(FARGS,"Can not open tmp file %s",tmpfile);
  }
  else
    fclose(tmpf);
  
  if ((dptr=getenv("DSSP")) == NULL)
    dptr="/usr/local/bin/dssp";
  if (!gmx_fexist(dptr))
    gmx_fatal(FARGS,"DSSP executable (%s) does not exist (use setenv DSSP)",
		dptr);
  sprintf(dssp,"%s %s %s %s > /dev/null %s",
	  dptr,bDoAccSurf?"":"-na",pdbfile,tmpfile,bVerbose?"":"2> /dev/null");
  if (bVerbose)
    fprintf(stderr,"dssp cmd='%s'\n",dssp);
  
  if (fnTArea) {
    fTArea=xvgropen(fnTArea,"Solvent Accessible Surface Area",
		    output_env_get_xvgr_tlabel(oenv),"Area (nm\\S2\\N)",oenv);
    xvgr_legend(fTArea,2,leg,oenv);
  } else
    fTArea=NULL;
  
  mat.map=NULL;
  mat.nmap=getcmap(libopen(opt2fn("-map",NFILE,fnm)),
		   opt2fn("-map",NFILE,fnm),&(mat.map));
  
  natoms=read_first_x(oenv,&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  if (natoms > atoms->nr) 
    gmx_fatal(FARGS,"\nTrajectory does not match topology!");
  if (gnx > natoms)
    gmx_fatal(FARGS,"\nTrajectory does not match selected group!");
  
  snew(average_area, atoms->nres);
  snew(av_area     , atoms->nres);
  snew(norm_av_area, atoms->nres);
  accr=NULL;
  naccr=0;
  
  gpbc = gmx_rmpbc_init(&top.idef,ePBC,natoms,box);
  do {
    t = output_env_conv_time(oenv,t);
    if (bDoAccSurf && nframe>=naccr) {
      naccr+=10;
      srenew(accr,naccr);
      for(i=naccr-10; i<naccr; i++)
        snew(accr[i],2*atoms->nres-1);
    }
    gmx_rmpbc(gpbc,natoms,box,x);
    tapein=ffopen(pdbfile,"w");
    write_pdbfile_indexed(tapein,NULL,atoms,x,ePBC,box,' ',-1,gnx,index,NULL,TRUE);
    ffclose(tapein);

#ifdef GMX_NO_SYSTEM
    printf("Warning-- No calls to system(3) supported on this platform.");
    printf("Warning-- Skipping execution of 'system(\"%s\")'.", dssp);
    exit(1);
#else
    if(0 != system(dssp))
    {
	gmx_fatal(FARGS,"Failed to execute command: %s",dssp);
    }
#endif

    /* strip_dssp returns the number of lines found in the dssp file, i.e.
     * the number of residues plus the separator lines */
    
    if (bDoAccSurf)
        accr_ptr = accr[nframe];

    nres_plus_separators = strip_dssp(tmpfile,nres,bPhbres,t,
                                      accr_ptr,fTArea,&mat,average_area,oenv);
    remove(tmpfile);
    remove(pdbfile);
    nframe++;
  } while(read_next_x(oenv,status,&t,natoms,x,box));
  fprintf(stderr,"\n");
  close_trj(status);
  if (fTArea)
    ffclose(fTArea);
  gmx_rmpbc_done(gpbc);

  prune_ss_legend(&mat);
  
  ss=opt2FILE("-o",NFILE,fnm,"w");
  mat.flags = 0;  
  write_xpm_m(ss,mat);
  ffclose(ss);
  
  if (opt2bSet("-ssdump",NFILE,fnm)) {
    snew(ss_str,nres+1);
    for(i=0; (i<nres); i++)
      ss_str[i] = mat.map[mat.matrix[0][i]].code.c1;
    ss_str[i] = '\0';
    ss = opt2FILE("-ssdump",NFILE,fnm,"w");
    fprintf(ss,"%d\n%s\n",nres,ss_str);
    ffclose(ss);
    sfree(ss_str);
  }
  analyse_ss(fnSCount,&mat,ss_string,oenv);

  if (bDoAccSurf) {
    write_sas_mat(fnArea,accr,nframe,nres_plus_separators,&mat);
  
    for(i=0; i<atoms->nres; i++)
      av_area[i] = (average_area[i] / (real)nframe);
    
    norm_acc(atoms, nres, av_area, norm_av_area);
    
    if (fnAArea) {
      acc=xvgropen(fnAArea,"Average Accessible Area",
		   "Residue","A\\S2",oenv);
      for(i=0; (i<nres); i++)
	fprintf(acc,"%5d  %10g %10g\n",i+1,av_area[i], norm_av_area[i]);
      ffclose(acc);
    }
  }

  view_all(oenv, NFILE, fnm);

  thanx(stderr);
  
  return 0;
}
