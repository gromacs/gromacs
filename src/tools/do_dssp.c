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
static char *SRCID_do_dssp_c = "$Id$";

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
#include "fatal.h"
#include "xvgr.h"
#include "matio.h"
#include "rdgroup.h"
#include "gstat.h"
#include "tpxio.h"

#ifdef MY_DSSP
extern void dssp_main(bool bDoAcc, bool bVerbose);
extern FILE *tapein, *tapeout;
#endif
		       
static void strip_dssp(char *dsspfile,int nres,
		       bool bPhobres[],real t,
		       real *acc,FILE *fTArea,
		       t_matrix *mat,int average_area[])
{
  static bool bFirst=TRUE;
  static char *ssbuf;
#ifndef MY_DSSP
  FILE *tapeout;
#endif
  static int xsize,frame;
  char buf[STRLEN+1];
  char SSTP;
  int  i,nr,iacc;
  real iaccf,iaccb;
  t_xpmelmt c;
  
#ifndef MY_DSSP
  tapeout=ffopen(dsspfile,"r");
#endif
  
  /* Skip header */
  do {
    fgets2(buf,STRLEN,tapeout);
  } while (strstr(buf,"KAPPA") == NULL);
  if (bFirst) {
    snew(ssbuf,nres+10);
  }
  
  iaccb=iaccf=0;
  for(nr=0; (fgets2(buf,STRLEN,tapeout) != NULL); nr++) {
    SSTP=buf[16]==' ' ? '~' : buf[16];
    ssbuf[nr]=SSTP;
    
    buf[39]='\0';
    sscanf(&(buf[34]),"%d",&iacc);
    acc[nr]=iacc;
    average_area[nr]+=iacc;
    if (bPhobres[nr])
      iaccb+=iacc;
    else
      iaccf+=iacc;
  }
  ssbuf[nr]='\0';
  
  if (bFirst) {
    sprintf(mat->title,"Secondary structure");
    mat->legend[0]=0;
    sprintf(mat->label_x,"Time (ps)");
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
    mat->matrix[frame][i]=searchcmap(mat->nmap,mat->map,c);
  }
  frame++;
  mat->nx=frame;
  
  if (fTArea)
    fprintf(fTArea,"%10g  %10g  %10g\n",t,0.01*iaccb,0.01*iaccf);
#ifndef MY_DSSP
  fclose(tapeout);
#endif
}

bool *bPhobics(t_atoms *atoms)
{
  int  i,nb;
  char **cb;
  bool *bb;
  
  nb=get_strings("phbres.dat",&cb);
  snew(bb,atoms->nres);
  
  for(i=0; (i<atoms->nres); i++) {
    if (search_str(nb,cb,*atoms->resname[i]) != -1)
      bb[i]=TRUE;
  }
  return bb;
}
 
static void check_oo(t_atoms *atoms)
{
  static char *OOO="O";
  int i;
  
  for(i=0; (i<atoms->nr); i++) {
    if (strcmp(*(atoms->atomname[i]),"OXT")==0)
      atoms->atomname[i]=&OOO;
    else if (strcmp(*(atoms->atomname[i]),"O1")==0)
      atoms->atomname[i]=&OOO;
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
    n = search_str(n_surf,surf_res,*atoms->resname[i]);
    if ( n != -1)
      norm_av_area[i] = av_area[i] / surf[n];
    else 
      fprintf(stderr,"Residue %s not found in surface database (%s)\n",
	      *atoms->resname[i],surffn);
  }
}

void prune_ss_legend(t_matrix *mat)
{
  bool *present;
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

void write_sas_mat(char *fn,real **accr,int nframe,int nres,t_matrix *mat)
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
    write_xpm(fp,"Solvent Accessible Surface","Surface (A^2)",
	      "Time","Residue Index",nframe,nres,
	      mat->axis_x,mat->axis_y,accr,lo,hi,rlo,rhi,&nlev);
    ffclose(fp);
  }
}

void analyse_ss(char *outfile, t_matrix *mat, char *ss_string)
{
  FILE *fp;
  t_mapping *map;
  int s,f,r,*count,ss_count;
  char **leg;
  
  map=mat->map;
  snew(count,mat->nmap);
  snew(leg,mat->nmap+1);
  leg[0]="Structure";
  for(s=0; s<mat->nmap; s++)
    leg[s+1]=map[s].desc;
  
  fp=xvgropen(outfile,"Secondary Structure",
	      mat->label_x,"Number of Residues");
  fprintf(fp,"@ subtitle \"Structure = ");
  for(s=0; s<strlen(ss_string); s++) {
    if (s>0)
      fprintf(fp," + ");
    for(f=0; f<mat->nmap; f++)
      if (ss_string[s]==map[f].code.c1)
	fprintf(fp,map[f].desc);
  }
  fprintf(fp,"\"\n");
  xvgr_legend(fp,mat->nmap+1,leg);
  
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
  
  fclose(fp);
  sfree(leg);
  sfree(count);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
#ifdef MY_DSSP
    "my_dssp ", 
#else
    "do_dssp ", 
#endif
    "reads a trajectory file and computes the secondary structure for",
    "each time frame ",
#ifdef MY_DSSP
    "using the dssp program.[PAR]",
#else
    "calling the dssp program. If you do not have the dssp program,",
    "get it. do_dssp assumes that the dssp executable is in",
    "/home/mdgroup/dssp/dssp. If that is not the case, then you should",
    "set an environment variable [BB]DSSP[bb] pointing to the dssp",
    "executable as in: [PAR]",
    "[TT]setenv DSSP /usr/local/bin/dssp[tt][PAR]",
#endif
    "The structure assignment for each residue and time is written to an",
    "[TT].xpm[tt] matrix file. This file can be visualized with for instance",
    "[TT]xv[tt] and can be converted to postscript with [TT]xpm2ps[tt].",
    "The number of residues with each secondary structure type and the",
    "total secondary structure ([TT]-sss[tt]) count as a function of",
    "time are also written to file ([TT]-sc[tt]).[PAR]",
    "Solvent accessible surface (SAS) per residue can be calculated, both in",
    "absolute values (A^2) and in fractions of the maximal accessible",
    "surface of a residue. The maximal accessible surface is defined as",
    "the accessible surface of a residue in a chain of glycines.",
    "[B]Note[b] that the program [TT]g_sas[tt] can also compute SAS",
    "and that is more efficient.[PAR]",
    "Finally, this program can dump the secondary structure in a special file",
    "[TT]ssdump.dat[tt] for usage in the program [TT]g_chi[tt]. Together",
    "these two programs can be used to analyze dihedral properties as a",
    "function of secondary structure type."
  };
  static bool bVerbose;
  static char *ss_string="HEBT"; 
  t_pargs pa[] = {
    { "-v",  FALSE, etBOOL, {&bVerbose},
      "HIDDENGenerate miles of useless information" },
    { "-sss", FALSE, etSTR, {&ss_string},
      "Secondary structures for structure count"}
  };
  
#ifndef MY_DSSP
  static char *bugs[] = { 
    "The program is very slow"
  };
#endif
  
  int        status;
#ifndef MY_DSSP
  FILE       *tapein;
#endif
  FILE       *ss,*acc,*fTArea;
  char       *fnSCount,*fnArea,*fnTArea,*fnAArea;
  char       *leg[] = { "Phobic", "Phylic" };
  t_topology top;
  t_atoms    *atoms;
  t_matrix   mat;
  int        nres,nr0,naccr;
  bool       *bPhbres,bDoAccSurf;
  real       t;
  int        i,j,natoms,nframe=0;
  matrix     box;
  int        gnx;
  char       *grpnm,*ss_str;
  atom_id    *index;
  rvec       *xp,*x;
  int        *average_area;
  real       **accr,*av_area, *norm_av_area;
  char       pdbfile[L_tmpnam],tmpfile[L_tmpnam],title[256];
#ifdef MY_DSSP
#define MAXBUF 1000000
  char       inbuf[MAXBUF],outbuf[MAXBUF];
#else
  char       dssp[256],*dptr;
#endif
  
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

  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,TRUE,
	 	    NFILE,fnm,asize(pa),pa,asize(desc),desc,
#ifdef MY_DSSP
		    0,NULL
#else
		    asize(bugs),bugs
#endif
		    );
  fnSCount= opt2fn("-sc",NFILE,fnm);
  fnArea  = opt2fn_null("-a", NFILE,fnm);
  fnTArea = opt2fn_null("-ta",NFILE,fnm);
  fnAArea = opt2fn_null("-aa",NFILE,fnm);
  bDoAccSurf=(fnArea || fnTArea || fnAArea);
  
  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&xp,NULL,box,FALSE);
  atoms=&(top.atoms);
  check_oo(atoms);
  bPhbres=bPhobics(atoms);
  
  get_index(atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&gnx,&index,&grpnm);
  nres=0;
  nr0=-1;
  for(i=0; (i<gnx); i++) {
    if (atoms->atom[index[i]].resnr != nr0) {
      nr0=atoms->atom[index[i]].resnr;
      nres++;
    }
  }
  fprintf(stderr,"There are %d residues in your selected group\n",nres);
  
  (void) tmpnam(pdbfile);
  (void) tmpnam(tmpfile);
#ifdef MY_DSSP
  /* Open all files read-write */
  tapein=ffopen(pdbfile,"w+");
  setvbuf(tapein,inbuf,_IOFBF,MAXBUF);
  tapeout=ffopen(tmpfile,"w+");
  setvbuf(tapeout,outbuf,_IOFBF,MAXBUF);
#else
  if ((dptr=getenv("DSSP")) == NULL)
    dptr="/home/mdgroup/dssp/dssp";
  if (!fexist(dptr))
    fatal_error(0,"DSSP executable (%s) does not exist (use setenv DSSP)",
		dptr);
  sprintf(dssp,"%s %s %s %s > /dev/null %s",
	  dptr,bDoAccSurf?"":"-na",pdbfile,tmpfile,bVerbose?"":"2> /dev/null");
  if (bVerbose)
    fprintf(stderr,"dssp cmd='%s'\n",dssp);
#endif
  
  if (fnTArea) {
    fTArea=xvgropen(fnTArea,"Solvent Accessible Surface Area",
		  "Time (ps)","Area (nm\\S2\\N)");
    xvgr_legend(fTArea,2,leg);
  } else
    fTArea=NULL;
  
  mat.map=NULL;
  mat.nmap=getcmap(libopen(opt2fn("-map",NFILE,fnm)),
		   opt2fn("-map",NFILE,fnm),&(mat.map));
  
  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  if (natoms > atoms->nr) 
    fatal_error(0,"\nTrajectory does not match topology!");
  if (gnx > natoms)
    fatal_error(0,"\nTrajectory does not match selected group!");
  
  snew(average_area,atoms->nres+10);
  snew(av_area,atoms->nres+10);
  snew(norm_av_area,atoms->nres+10);
  accr=NULL;
  naccr=0;
  do {
    if (nframe>=naccr) {
      naccr+=10;
      srenew(accr,naccr);
      for(i=naccr-10; i<naccr; i++)
	snew(accr[i],atoms->nres);
    }
    rm_pbc(&(top.idef),natoms,box,x,x);
#ifndef MY_DSSP
    tapein=ffopen(pdbfile,"w");
#endif
    write_pdbfile_indexed(tapein,NULL,atoms,x,box,0,-1,gnx,index);
#ifdef MY_DSSP
    rewind(tapein);
    dssp_main(bDoAccSurf,bVerbose);
    rewind(tapein);
    rewind(tapeout);
#else
    fclose(tapein);
    system(dssp);
#endif
    strip_dssp(tmpfile,nres,bPhbres,t,
	       accr[nframe],fTArea,&mat,average_area);
#ifdef MY_DSSP
    rewind(tapeout);
#else
    remove(tmpfile);
    remove(pdbfile);
#endif
    nframe++;
  } while(read_next_x(status,&t,natoms,x,box));
  fprintf(stderr,"\n");
  close_trj(status);
  if (fTArea)
    ffclose(fTArea);
  
  prune_ss_legend(&mat);
  
  ss=opt2FILE("-o",NFILE,fnm,"w");
  write_xpm_m(ss,mat);
  ffclose(ss);
  
  if (opt2bSet("-ssdump",NFILE,fnm)) {
    snew(ss_str,nres+1);
    for(i=0; (i<nres); i++)
      ss_str[i] = mat.map[mat.matrix[0][i]].code.c1;
    ss_str[i] = '\0';
    ss = opt2FILE("-ssdump",NFILE,fnm,"w");
    fprintf(ss,"%d\n%s\n",nres,ss_str);
    fclose(ss);
    sfree(ss_str);
  }
  analyse_ss(fnSCount,&mat,ss_string);

  if (bDoAccSurf) {
    write_sas_mat(fnArea,accr,nframe,nres,&mat);
  
    for(i=0; i<atoms->nres; i++)
      av_area[i] = (average_area[i] / (real)nframe);
    
    norm_acc(atoms, nres, av_area, norm_av_area);
    
    if (fnAArea) {
      acc=xvgropen(fnAArea,"Average Accessible Area",
		   "Residue","A\\S2");
      for(i=0; (i<nres); i++)
	fprintf(acc,"%5d  %10g %10g\n",i+1,av_area[i], norm_av_area[i]);
      ffclose(acc);
    }
  }
      
  if (fnTArea)  xvgr_file(fnTArea ,NULL);
  if (fnAArea)  xvgr_file(fnAArea ,NULL);
  if (fnSCount) xvgr_file(fnSCount,NULL);

  thanx(stderr);
  
  return 0;
}
