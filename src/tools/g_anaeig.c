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
 * Good ROcking Metal Altar for Chronical Sinners
 */
static char *SRCID_g_nmeig_c = "$Id$";

#include <math.h>
#include <string.h>
#include "statutil.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "fatal.h"
#include "vec.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "rdgroup.h"
#include "pdbio.h"
#include "tpxio.h"
#include "trnio.h"
#include "matio.h"
#include "mshift.h"
#include "xvgr.h"
#include "gstat.h"
#include "txtdump.h"

real tick_spacing(real range,int minticks)
{
  real sp;

  sp = 0.2*exp(log(10)*ceil(log(range)/log(10)));
  while (range/sp < minticks-1)
    sp = sp/2;

  return sp;
}

void write_xvgr_graphs(char *file,int ngraphs,
		       char *title,char *xlabel,char **ylabel,
		       int n,real *x, real **y,bool bZero)
{
  FILE *out;
  int g,i;
  real min,max,xsp,ysp;
  
  out=ffopen(file,"w"); 
  for(g=0; g<ngraphs; g++) {
    min=y[g][0];
    max=y[g][0];
    for(i=0; i<n; i++) {
      if (y[g][i]<min) min=y[g][i];
      if (y[g][i]>max) max=y[g][i];
    }
    if (bZero)
      min=0;
    else
      min=min-0.1*(max-min);
    max=max+0.1*(max-min);
    xsp=tick_spacing(x[n-1]-x[0],4);
    ysp=tick_spacing(max-min,3);
    fprintf(out,"@ with g%d\n@ g%d on\n",ngraphs-1-g,ngraphs-1-g);
    fprintf(out,"@ g%d autoscale type AUTO\n",ngraphs-1-g);
    if (g==0)
      fprintf(out,"@ title \"%s\"\n",title);
    if (g==ngraphs-1)
      fprintf(out,"@ xaxis  label \"%s\"\n",xlabel);
    else 
      fprintf(out,"@ xaxis  ticklabel off\n");
    fprintf(out,"@ world xmin %g\n",x[0]);
    fprintf(out,"@ world xmax %g\n",x[n-1]);
    fprintf(out,"@ world ymin %g\n",min);
    fprintf(out,"@ world ymax %g\n",max);
    fprintf(out,"@ view xmin 0.15\n");
    fprintf(out,"@ view xmax 0.85\n");
    fprintf(out,"@ view ymin %g\n",0.15+(ngraphs-1-g)*0.7/ngraphs);
    fprintf(out,"@ view ymax %g\n",0.15+(ngraphs-g)*0.7/ngraphs);
    fprintf(out,"@ yaxis  label \"%s\"\n",ylabel[g]);
    fprintf(out,"@ xaxis tick major %g\n",xsp);
    fprintf(out,"@ xaxis tick minor %g\n",xsp/2);
    fprintf(out,"@ xaxis ticklabel start type spec\n");
    fprintf(out,"@ xaxis ticklabel start %g\n",ceil(min/xsp)*xsp);
    fprintf(out,"@ yaxis tick major %g\n",ysp);
    fprintf(out,"@ yaxis tick minor %g\n",ysp/2);
    fprintf(out,"@ yaxis ticklabel start type spec\n");
    fprintf(out,"@ yaxis ticklabel start %g\n",ceil(min/ysp)*ysp);
    if ((min<0) && (max>0)) {
      fprintf(out,"@ zeroxaxis bar on\n");
      fprintf(out,"@ zeroxaxis bar linestyle 3\n");
    }
    for(i=0; i<n; i++) 
      fprintf(out,"%10.4f %10.5f\n",x[i],y[g][i]);
    fprintf(out,"&\n");
  }
  fclose(out);
}

void inprod_matrix(char *matfile,int natoms,
		   int nvec1,int *eignr1,rvec **eigvec1,
		   int nvec2,int *eignr2,rvec **eigvec2)
{
  FILE *out;
  real **mat;
  int i,x,y,nlevels;
  real inp,*t_x,*t_y,max;
  t_rgb rlo,rhi;

  fprintf(stderr,"Calculating inner-product matrix of %dx%d eigenvectors...\n",
	  nvec1,nvec2);
  
  snew(mat,nvec1);
  for(x=0; x<nvec1; x++)
    snew(mat[x],nvec2);

  max=0;
  for(x=0; x<nvec1; x++)
    for(y=0; y<nvec2; y++) {
      inp=0;
      for(i=0; i<natoms; i++)
	inp+=iprod(eigvec1[x][i],eigvec2[y][i]);
      mat[x][y]=fabs(inp);
      if (mat[x][y]>max)
	max=mat[x][y];
    }
  snew(t_x,nvec1);
  for(i=0; i<nvec1; i++)
    t_x[i]=eignr1[i]+1;
  snew(t_y,nvec2);
  for(i=0; i<nvec2; i++)
    t_y[i]=eignr2[i]+1;
  rlo.r = 1; rlo.g = 1; rlo.b = 1;
  rhi.r = 0; rhi.g = 0; rhi.b = 0;
  nlevels=41;
  out=ffopen(matfile,"w");
  write_xpm(out,"Eigenvector inner-products","in.prod.","run 1","run 2",
	    nvec1,nvec2,t_x,t_y,mat,0.0,max,rlo,rhi,&nlevels);
  fclose(out);
}

void overlap(char *outfile,int natoms,
	     int nvec1,int *eignr1,rvec **eigvec1,
	     int nvec2,int *eignr2,rvec **eigvec2,
	     int noutvec,int *outvec)
{
  FILE *out;
  int i,v,vec,x;
  real overlap,inp;

  fprintf(stderr,"Calculating overlap between eigenvectors of set 2 with eigenvectors\n");
  for(i=0; i<noutvec; i++)
    fprintf(stderr,"%d ",outvec[i]+1);
  fprintf(stderr,"\n");

  out=xvgropen(outfile,"Subspace overlap",
	       "Eigenvectors of trajectory 2","Overlap");
  fprintf(out,"@ subtitle \"using %d eigenvectors of trajectory 1\"\n",
	  noutvec);
  overlap=0;
  for(x=0; x<nvec2; x++) {
    for(v=0; v<noutvec; v++) {
      vec=outvec[v];
      inp=0;
      for(i=0; i<natoms; i++)
	inp+=iprod(eigvec1[vec][i],eigvec2[x][i]);
      overlap+=sqr(inp);
    }
    fprintf(out,"%5d  %5.3f\n",eignr2[x]+1,overlap/noutvec);
  }

  fclose(out);
}

void project(char *trajfile,t_topology *top,matrix topbox,rvec *xtop,
	     char *projfile,char *twodplotfile,char *filterfile,int skip,
	     char *extremefile,bool bExtrAll,real extreme,int nextr,
	     t_atoms *atoms,int natoms,atom_id *index,
	     rvec *xref,int nfit,atom_id *ifit,real *w_rls,
	     real *sqrtm,rvec *xav,
	     int nvec,int *eignr,rvec **eigvec,
	     int noutvec,int *outvec)
{
  FILE    *xvgrout;
  int     status,out,nat,i,j,d,v,vec,nfr,nframes,snew_size,frame;
  int     *imin,*imax;
  atom_id *all_at;
  matrix  box;
  rvec    *xread,*x;
  real    t,inp,**inprod,min,max;
  char    str[STRLEN],str2[STRLEN],**ylabel,*c;
  
  snew(x,natoms);
  
  if (!bExtrAll)
    noutvec=1;
  if (trajfile) {
    snew(inprod,noutvec+1);
    
    if (filterfile) {
      fprintf(stderr,"Writing a filtered trajectory to %s using eigenvectors\n",
	      filterfile);
      for(i=0; i<noutvec; i++)
	fprintf(stderr,"%d ",outvec[i]+1);
      fprintf(stderr,"\n");
      out=open_trx(filterfile,"w");
    }
    snew_size=0;
    nfr=0;
    nframes=0;
    nat=read_first_x(&status,trajfile,&t,&xread,box);
    snew(all_at,nat);
    for(i=0; (i<nat); i++)
      all_at[i]=i;
    do {
      if (nfr % skip == 0) {
	if (top)
	  rm_pbc(&(top->idef),nat,box,xread,xread);
	if (nframes>=snew_size) {
	  snew_size+=100;
	  for(i=0; i<noutvec+1; i++)
	    srenew(inprod[i],snew_size);
	}
	inprod[noutvec][nframes]=t;
	/* calculate x: a fitted struture of the selected atoms */
	if (xref==NULL) {
	  reset_x(nfit,ifit,nat,all_at,xread,w_rls);
	  do_fit(atoms->nr,w_rls,xtop,xread);
	}
	for (i=0; i<natoms; i++)
	  copy_rvec(xread[index[i]],x[i]);
	if (xref) {
	  reset_x(natoms,all_at,natoms,all_at,x,w_rls);
	  do_fit(natoms,w_rls,xref,x);
	}
	
	for(v=0; v<noutvec; v++) {
	  vec=outvec[v];
	  /* calculate (mass-weighted) projection */
	  inp=0;
	  for (i=0; i<natoms; i++) {
	    inp+=(eigvec[vec][i][0]*(x[i][0]-xav[i][0])+
	    eigvec[vec][i][1]*(x[i][1]-xav[i][1])+
	    eigvec[vec][i][2]*(x[i][2]-xav[i][2]))*sqrtm[i];
	  }
	  inprod[v][nframes]=inp;
	  if (filterfile) 
	    for(i=0; i<natoms; i++)
	      for(d=0; d<DIM; d++)
		xread[index[i]][d] = xav[i][d]+
		  inprod[v][nframes]*eigvec[outvec[v]][i][d]/sqrtm[i];
	}
	if (filterfile && (nframes % skip == 0)) 
	  write_trx(out,natoms,index,atoms,0,t,box,xread,NULL);
	nframes++;
      }
      nfr++;
    } while (read_next_x(status,&t,nat,xread,box));
    close_trj(status);
     sfree(x);
     if (filterfile)
       close_trx(out);
  }
  else
    snew(xread,atoms->nr);
  

  if (projfile) {
    snew(ylabel,noutvec);
    for(v=0; v<noutvec; v++) {
      sprintf(str,"vec %d",eignr[outvec[v]]+1);
      ylabel[v]=strdup(str);
    }
    write_xvgr_graphs(projfile,noutvec,
		      "Projection on eigenvectors (nm)","Time (ps)",
		      ylabel,nframes,inprod[noutvec],inprod,FALSE);
  }
  if (twodplotfile) {
    sprintf(str,"projection on eigenvector %d (nm)",
	    eignr[outvec[0]]+1);
    sprintf(str2,"projection on eigenvector %d (nm)",
	    eignr[outvec[noutvec-1]]+1); 
    xvgrout=xvgropen(twodplotfile,"2D projection of trajectory",str,str2);
    for(i=0; i<nframes; i++)
      fprintf(xvgrout,"%10.5f %10.5f\n",inprod[0][i],inprod[noutvec-1][i]);
    fclose(xvgrout);
  }
  if (extremefile) {
    if (extreme==0) {
      fprintf(stderr,"%11s %17s %17s\n","eigenvector","Minimum","Maximum");
      fprintf(stderr,
	      "%11s %10s %10s %10s %10s\n","","value","time","value","time");
      snew(imin,noutvec);
      snew(imax,noutvec);
      for(v=0; v<noutvec; v++) {
	for(i=0; i<nframes; i++) {
	  if (inprod[v][i]<inprod[v][imin[v]])
	    imin[v]=i;
	  if (inprod[v][i]>inprod[v][imax[v]])
	    imax[v]=i;
	}
	min=inprod[v][imin[v]];
	max=inprod[v][imax[v]];
	fprintf(stderr,"%7d     %10.6f %10.1f %10.6f %10.1f\n",
		eignr[outvec[v]]+1,
		min,inprod[noutvec][imin[v]],max,inprod[noutvec][imax[v]]); 
      }
    }
    else {
      min=-extreme;
      max=+extreme;
    }
    /* build format string for filename: */
    strcpy(str,extremefile);/* copy filename */
    c=strrchr(str,'.'); /* find where extention begins */
    strcpy(str2,c); /* get extention */
    sprintf(c,"%%d%s",str2); /* append '%s' and extention to filename */
    for(v=0; v<noutvec; v++) {
      /* make filename using format string */
      if (noutvec==1)
	strcpy(str2,extremefile);
      else
	sprintf(str2,str,eignr[outvec[v]]+1);
      fprintf(stderr,"Writing %d frames between the extremes to %s\n",
	      nextr,str2);
      out=open_trx(str2,"w");
      for(frame=0; frame<nextr; frame++) {
	if ((extreme==0) && (nextr<=3))
	  for(i=0; i<natoms; i++)
	    atoms->atom[index[i]].chain='A'+frame;
	for(i=0; i<natoms; i++)
	  for(d=0; d<DIM; d++) 
	    xread[index[i]][d] = 
	      (xav[i][d] + (min*(nextr-frame-1)+max*frame)/(nextr-1)
	      *eigvec[outvec[v]][i][d]/sqrtm[i]);
	write_trx(out,natoms,index,atoms,0,frame,topbox,xread,NULL);
      }
      close_trx(out);
    }
  }
  fprintf(stderr,"\n");
}

void components(char *outfile,int natoms,real *sqrtm,
		int nvec,int *eignr,rvec **eigvec,
		int noutvec,int *outvec)
{
  int g,v,i;
  real *x,**y;
  char str[STRLEN],**ylabel;

  fprintf(stderr,"Writing eigenvector components to %s\n",outfile);

  snew(ylabel,noutvec);
  snew(y,noutvec);
  snew(x,natoms);
  for(i=0; i<natoms; i++)
    x[i]=i+1;
  for(g=0; g<noutvec; g++) {
    v=outvec[g];
    sprintf(str,"vec %d",eignr[v]+1);
    ylabel[g]=strdup(str);
    snew(y[g],natoms);
    for(i=0; i<natoms; i++)
      y[g][i] = norm(eigvec[v][i])/sqrtm[i];
  }
  write_xvgr_graphs(outfile,noutvec,"Eigenvector components","Atom number",
		    ylabel,natoms,x,y,TRUE);
  fprintf(stderr,"\n");
}

void read_eigenvectors(char *file,int *natoms,rvec **xref,bool *bDMR,
		       rvec **xav,bool *bDMA,
		       int *nvec, int **eignr, rvec ***eigvec)
{
  t_trnheader head;
  int status,i,snew_size;
  rvec *x;
  matrix box;
  bool bOK;

  *bDMR=FALSE;

  /* read (reference (t=-1) and) average (t=0) structure */
  status=open_trn(file,"r");
  fread_trnheader(status,&head,&bOK);
  *natoms=head.natoms;
  snew(*xav,*natoms);
  fread_htrn(status,&head,box,*xav,NULL,NULL);
  if ((head.t>=-1.1) && (head.t<=-0.9)) {
    snew(*xref,*natoms);
    for(i=0; i<*natoms; i++)
      copy_rvec((*xav)[i],(*xref)[i]);
    *bDMR = (head.lambda > 0.5);
    fprintf(stderr,"Read %smass weighted reference structure with %d atoms from %s\n",
	    *bDMR ? "" : "non ",*natoms,file);
    fread_trnheader(status,&head,&bOK);
    fread_htrn(status,&head,box,*xav,NULL,NULL);
  }
  else
    *xref=NULL;
  *bDMA = (head.lambda > 0.5);
  if ((head.t<=-0.01) || (head.t>=0.01))
    fatal_error(0,"%s does not start with t=0, which should be the average "
		"structure. This might not be a eigenvector file.",file);
  fprintf(stderr,"Read %smass weighted average structure with %d atoms from %s\n",
	  *bDMA ? "" : "non ",*natoms,file);
  
  snew(x,*natoms);
  snew_size=0;
  *nvec=0;
  while (fread_trnheader(status,&head,&bOK)) {
    fread_htrn(status,&head,box,x,NULL,NULL);
    if (*nvec >= snew_size) {
      snew_size+=10;
      srenew(*eignr,snew_size);
      srenew(*eigvec,snew_size);
    }
    i=(int)(head.t+0.01);
    if ((head.t-i<=-0.01) || (head.t-i>=0.01))
      fatal_error(0,"%s contains a frame with non-integer time (%f), this "
		  "time should be an eigenvector index. "
		  "This might not be a eigenvector file.",file,head.t);
    (*eignr)[*nvec]=i-1;
    snew((*eigvec)[*nvec],*natoms);
    for(i=0; i<*natoms; i++)
      copy_rvec(x[i],(*eigvec)[*nvec][i]);
    (*nvec)++;
  }
  sfree(x);
  fprintf(stderr,"Read %d eigenvectors (dim=%d)\n\n",*nvec,*natoms*DIM);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "[TT]g_anaeig[tt] analyzes eigenvectors of a covariance matrix which",
    "can be calculated with [TT]g_covar[tt].",
    "All structures are fitted to the structure in the eigenvector file,",
    "if present, otherwise to the structure in the structure file.",
    "When no run input file is supplied, periodicity will not be taken into",
    "account.",  
    "Most analyses are done on eigenvectors [TT]-first[tt] to [TT]-last[tt],",
    "but when [TT]-first[tt] is set to -1 you will be prompted for a",
    "selection.[PAR]",
    "[TT]-comp[tt]: plot eigenvector components per atom of eigenvectors",
    "[TT]-first[tt] to [TT]-last[tt].[PAR]",
    "[TT]-proj[tt]: calculate projections of a trajectory on eigenvectors",
    "[TT]-first[tt] to [TT]-last[tt].[PAR]",
    "[TT]-2d[tt]: calculate a 2d projection of a trajectory on eigenvectors",
    "[TT]-first[tt] and [TT]-last[tt].[PAR]",
    "[TT]-filt[tt]: filter the trajectory to show only the motion along",
    "eigenvectors [TT]-first[tt] to [TT]-last[tt].[PAR]",
    "[TT]-extr[tt]: calculate the two extreme projections along a trajectory",
    "on the average structure and interpolate [TT]-nframes[tt] frames between",
    "them, or set your own extremes with [TT]-max[tt]. The eigenvector",
    "[TT]-first[tt] will be written unless [TT]-first[tt] and [TT]-last[tt]",
    "are explicitly set, in which case all eigenvectors will be written",
    "to separate files.",
    "Chain identifiers will be added when",
    "writing a [TT].pdb[tt] file with two or three strcutures",
    "(you can use [TT]rasmol -nmrpdb[tt] to view such a pdb file).[PAR]",
    "[TT]-over[tt]: calculate the subspace overlap of the eigenvectors in",
    "file [TT]-v2[tt] with eigenvectors [TT]-first[tt] to [TT]-last[tt]",
    "in file [TT]-v[tt].[PAR]",
    "[TT]-inpr[tt]: calculate a matrix of inner-products between eigenvectors",
    "in files [TT]-v[tt] and [TT]-v2[tt]."
  };
  static int  first=1,last=8,skip=1,nextr=2;
  static real max=0.0;
  t_pargs pa[] = {
    { "-first", FALSE, etINT, &first,     
      "First eigenvector for analysis (-1 is select)" },
    { "-last",  FALSE, etINT, &last, 
      "Last eigenvector for analysis (-1 is till the last)" },
     { "-skip",  FALSE, etINT, &skip,
      "Only analyse every nr-th frame" },
    { "-max",  FALSE, etREAL, &max, 
      "Maximum for projection of the eigenvector on the average structure, max=0 gives the extremes" },
    { "-nframes",  FALSE, etINT, &nextr, 
      "Number of frames for the extremes output" }
  };
#define NPA asize(pa)
  
  FILE       *out;
  int        status,trjout;
  t_topology top;
  t_atoms    *atoms;
  rvec       *xtop,*xref1,*xref2;
  bool       bDMR1,bDMA1,bDMR2,bDMA2;
  int        nvec1,nvec2,*eignr1=NULL,*eignr2=NULL;
  rvec       *x,*xread,*xav1,*xav2,**eigvec1=NULL,**eigvec2=NULL;
  matrix     topbox;
  real       xid,totmass,*sqrtm,avsqrtm,*w_rls,t,lambda;
  int        natoms,step;
  char       *grpname,*indexfile,title[STRLEN];
  int        i,j,d;
  int        nout,*iout,noutvec,*outvec,nfit;
  atom_id    *index,*ifit;
  char       *Vec2File,*topfile,*CompFile,*ProjOnVecFile,*TwoDPlotFile;
  char       *FilterFile,*ExtremeFile;
  char       *OverlapFile,*InpMatFile;
  bool       bM,bIndex,bTPS,bTop,bVec2,bProj,bFirstToLast,bExtremeAll,bTraj;
  t_filenm fnm[] = { 
    { efTRN, "-v",    "eigenvec",    ffREAD  },
    { efTRN, "-v2",   "eigenvec2",   ffOPTRD },
    { efTRX, "-f",    NULL,          ffOPTRD }, 
    { efTPS, NULL,    NULL,          ffOPTRD },
    { efNDX, NULL,    NULL,          ffOPTRD },
    { efXVG, "-comp", "eigcomp",     ffOPTWR },
    { efXVG, "-proj", "proj",        ffOPTWR },
    { efXVG, "-2d",   "2dproj",      ffOPTWR },
    { efTRX, "-filt", "filtered",    ffOPTWR },
    { efTRX, "-extr", "extreme.pdb", ffOPTWR },
    { efXVG, "-over", "overlap",     ffOPTWR },
    { efXPM, "-inpr", "inprod",      ffOPTWR }
  }; 
#define NFILE asize(fnm) 

  CopyRight(stderr,argv[0]); 
  parse_common_args(&argc,argv,PCA_CAN_TIME,TRUE,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL); 

  indexfile=ftp2fn_null(efNDX,NFILE,fnm);

  Vec2File        = opt2fn_null("-v2",NFILE,fnm);
  topfile         = ftp2fn(efTPS,NFILE,fnm); 
  CompFile        = opt2fn_null("-comp",NFILE,fnm);
  ProjOnVecFile   = opt2fn_null("-proj",NFILE,fnm);
  TwoDPlotFile    = opt2fn_null("-2d",NFILE,fnm);
  FilterFile      = opt2fn_null("-filt",NFILE,fnm);
  ExtremeFile     = opt2fn_null("-extr",NFILE,fnm);
  OverlapFile     = opt2fn_null("-over",NFILE,fnm);
  InpMatFile      = ftp2fn_null(efXPM,NFILE,fnm);
  bTop   = fn2bTPX(topfile);
  bProj  = ProjOnVecFile || FilterFile || ExtremeFile || TwoDPlotFile;
  bExtremeAll  = 
    opt2parg_bSet("-first",NPA,pa) && opt2parg_bSet("-last",NPA,pa);
  bFirstToLast = CompFile || ProjOnVecFile || FilterFile || OverlapFile || 
    ( ExtremeFile && bExtremeAll );
  bVec2  = Vec2File || OverlapFile || InpMatFile;
  bM     = CompFile || bProj;
  bTraj  = ProjOnVecFile || FilterFile || (ExtremeFile && (max==0))
    || TwoDPlotFile;
  bIndex = bM || bProj;
  bTPS   = ftp2bSet(efTPS,NFILE,fnm) || bM || bTraj ||
    FilterFile  || (bIndex && indexfile);

  read_eigenvectors(opt2fn("-v",NFILE,fnm),&natoms,&xref1,&bDMR1,&xav1,&bDMA1,
		    &nvec1,&eignr1,&eigvec1);
  if (bVec2) {
    read_eigenvectors(Vec2File,&i,&xref2,&bDMR2,&xav2,&bDMA2,
		      &nvec2,&eignr2,&eigvec2);
    if (i!=natoms)
      fatal_error(0,"Dimensions in the eigenvector files don't match");
  }

  if (xref1 && !bDMR1 && !bDMA1) 
    bM=FALSE;
  if ((xref1==NULL) && (bM || bTraj))
    bTPS=TRUE;
    
  xtop=NULL;
  nfit=0;
  ifit=NULL;
  w_rls=NULL;
  if (bTPS) {
    bTop=read_tps_conf(ftp2fn(efTPS,NFILE,fnm),
		       title,&top,&xtop,NULL,topbox,TRUE);
    atoms=&top.atoms;
    if (bTop)
      rm_pbc(&(top.idef),atoms->nr,topbox,xtop,xtop);
    if (xref1==NULL || (bM && bDMR1)) {
      printf("\nNote: the structure in %s should be the same\n"
	       "      as the one used for the fit in g_covar\n",topfile);
      printf("\nSelect the index group that was used for the least squares fit in g_covar\n",natoms);
      get_index(atoms,indexfile,1,&nfit,&ifit,&grpname);
      if (xref1==NULL) {
	snew(w_rls,atoms->nr);
	for(i=0; (i<nfit); i++)
	  if (bM)
	    w_rls[ifit[i]]=atoms->atom[ifit[i]].m;
	  else
	    w_rls[ifit[i]]=1.0;
      }
      else {
	/* make the fit index in xref instead of xtop */
	snew(w_rls,nfit);
	for(i=0; (i<nfit); i++) {
	  w_rls[i]=atoms->atom[ifit[i]].m;
	  ifit[i]=i;
	}
      }
    }
    else {
      /* make the fit non mass weighted on xref */
      nfit=natoms;
      snew(ifit,nfit);
      snew(w_rls,nfit);
      for(i=0; i<nfit; i++) {
	ifit[i]=i;
	w_rls[i]=1.0;
      }
    }
  }
  else
    bTop=FALSE;

  if (bIndex) {
    printf("\nSelect an index group of %d elements that corresponds to the eigenvectors\n",natoms);
    get_index(atoms,indexfile,1,&i,&index,&grpname);
    if (i!=natoms)
      fatal_error(0,"you selected a group with %d elements instead of %d",
		  i,natoms);
    printf("\n");
  }

  snew(sqrtm,natoms);
  if (bM && bDMA1) {
    avsqrtm=0;
    for(i=0; (i<natoms); i++) {
      sqrtm[i]=sqrt(atoms->atom[index[i]].m);
      avsqrtm+=sqrtm[i];
    }
    avsqrtm=natoms/avsqrtm;
    for(i=0; (i<natoms); i++) 
      sqrtm[i]=sqrtm[i]*avsqrtm;
  } else
    for(i=0; (i<natoms); i++)
      sqrtm[i]=1.0;
  
  if (bVec2) {
    t=0;
    totmass=0;
    for(i=0; (i<natoms); i++)
      for(d=0;(d<DIM);d++) {
	t+=sqr((xav1[i][d]-xav2[i][d])*sqrtm[i]);
	totmass+=sqr(sqrtm[i]);
      }
    fprintf(stderr,"RMSD (without fit) between the two average structures:"
	    " %.3f (nm)\n\n",sqrt(t/totmass));
  }
  
  if (last==-1)
    last=natoms*DIM;
  if (first>-1) {
    if (bFirstToLast) {
      /* make an index from first to last */
      nout=last-first+1;
      snew(iout,nout);
      for(i=0; i<nout; i++)
	iout[i]=first-1+i;
    } else {
      /* make an index of first and last */
      nout=2;
      snew(iout,nout);
      iout[0]=first-1;
      iout[1]=last-1;
    }
  }
  else {
    printf("Select eigenvectors for output, end your selection with 0\n");
    nout=-1;
    iout=NULL;
    do {
      nout++;
      srenew(iout,nout+1);
      scanf("%d",&iout[nout]);
      iout[nout]--;
    } while (iout[nout]>=0);
    printf("\n");
  }
  /* make an index of the eigenvectors which are present */
  snew(outvec,nout);
  noutvec=0;
  for(i=0; i<nout; i++) {
    j=0;
    while ((j<nvec1) && (eignr1[j]!=iout[i]))
      j++;
    if ((j<nvec1) && (eignr1[j]==iout[i])) {
      outvec[noutvec]=j;
      noutvec++;
    }
  }
  fprintf(stderr,"%d eigenvectors selected for output:",noutvec);
  for(j=0; j<noutvec; j++)
    fprintf(stderr," %d",eignr1[outvec[j]]+1);
  fprintf(stderr,"\n");
  
  if (CompFile)
    components(CompFile,natoms,sqrtm,nvec1,eignr1,eigvec1,noutvec,outvec);
  
  if (bProj)
    project(bTraj ? opt2fn("-f",NFILE,fnm) : NULL,
	    bTop ? &top : NULL,topbox,xtop,
	    ProjOnVecFile,TwoDPlotFile,FilterFile,skip,
	    ExtremeFile,bExtremeAll,max,nextr,atoms,natoms,index,
	    xref1,nfit,ifit,w_rls,
	    sqrtm,xav1,nvec1,eignr1,eigvec1,noutvec,outvec);
  
  if (OverlapFile)
    overlap(OverlapFile,natoms,
	    nvec1,eignr1,eigvec1,nvec2,eignr2,eigvec2,noutvec,outvec);
  
  if (InpMatFile)
    inprod_matrix(InpMatFile,natoms,
		  nvec1,eignr1,eigvec1,nvec2,eignr2,eigvec2);
  
  return 0;
}
  
