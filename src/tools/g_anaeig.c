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
#include "trnio.h"
#include "matio.h"
#include "mshift.h"
#include "xvgr.h"
#include "gstat.h"
#include "txtdump.h"

void write_xvgr_graphs(char *file,int ngraphs,
		       char *title,char *xlabel,char **ylabel,
		       int n,real *x, real **y,bool bZero)
{
  FILE *out;
  int g,i;
  real min,max;
  
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
    fprintf(out,"@ with g%d\n@ g%d on\n",ngraphs-1-g,ngraphs-1-g);
    fprintf(out,"@g%d autoscale type AUTO\n",ngraphs-1-g);
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
    fprintf(out,"@ view xmax 0.90\n");
    fprintf(out,"@ view ymin %g\n",0.15+(ngraphs-1-g)*0.7/ngraphs);
    fprintf(out,"@ view ymax %g\n",0.15+(ngraphs-g)*0.7/ngraphs);
    fprintf(out,"@ yaxis  label \"%s\"\n",ylabel[g]);
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

  out=xvgropen(outfile,"Cumulative inner-products",
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

void project(char *trajfile,char *projfile,
	     char *twodplotfile,char *filterfile,int skip,
	     char *extremefile,real extreme,int nextr,
	     t_atoms *atoms,int natoms,atom_id *index,
	     int nref,rvec *xref,int nfit,atom_id *ifit,real *w_rls,
	     real *mass,real *sqrtm,rvec *xav,
	     int nvec,int *eignr,rvec **eigvec,
	     int noutvec,int *outvec)
{
  FILE    *xvgrout;
  int     status,out,nat,i,j,d,v,vec,nframes,snew_size,min,max,frame;
  atom_id *all_at;
  matrix  box,zerobox;
  rvec    *xread,*x,tmp;
  real    t,inp,**inprod;
  char    str[STRLEN],str2[STRLEN],**ylabel;

  clear_mat(zerobox);

  snew(x,natoms);
  snew(all_at,natoms);
  for(i=0; (i<natoms); i++)
    all_at[i]=i;
  
  snew(inprod,noutvec+1);

  for(i=0; (i<natoms); i++)
    all_at[i]=i; 
  if (filterfile) {
    fprintf(stderr,"Writing a filtered trajectory to %s using eigenvectors\n",
	    filterfile);
    for(i=0; i<noutvec; i++)
      fprintf(stderr,"%d ",outvec[i]+1);
    fprintf(stderr,"\n");
    out=open_trx(filterfile,"w");
  }
  snew_size=0;
  nframes=0;
  nat=read_first_x(&status,trajfile,&t,&xread,box);
  snew(all_at,nat);
  do {
    if (nframes>=snew_size) {
      snew_size+=100;
      for(i=0; i<noutvec+1; i++)
	srenew(inprod[i],snew_size);
    }
    inprod[noutvec][nframes]=t;
    /* calculate x: a fitted struture of the selected atoms */
    if (nref==-1) {
      for (i=0; i<natoms; i++)
	copy_rvec(xread[index[i]],x[i]);
      reset_x(natoms,all_at,natoms,all_at,x,mass);
      do_fit(natoms,mass,xav,x);
    } 
    else {
      reset_x(nfit,ifit,nat,all_at,xread,w_rls);
      do_fit(nat,w_rls,xref,xread);
      for (i=0; i<natoms; i++)
	copy_rvec(xread[index[i]],x[i]);
    }

    clear_rvec(tmp);
    for (i=0; i<natoms; i++)
      rvec_inc(tmp,x[i]);

    for(v=0; v<noutvec; v++) {
      vec=outvec[v];
      /* calculate (mass-weighted) projection */
      inp=0;
      for (i=0; i<natoms; i++)
	inp+=(eigvec[vec][i][0]*(x[i][0]-xav[i][0])+
              eigvec[vec][i][1]*(x[i][1]-xav[i][1])+
	      eigvec[vec][i][2]*(x[i][2]-xav[i][2]))*sqrtm[i];
      inprod[v][nframes]=inp;
      if (filterfile && (nframes % skip == 0)) 
	for(i=0; i<natoms; i++)
	  for(d=0; d<DIM; d++)
	    xread[index[i]][d] = xav[i][d]+
	      inprod[v][nframes]*eigvec[outvec[v]][i][d]/sqrtm[i];
    }
    if (filterfile && (nframes % skip == 0)) 
      write_trx(out,natoms,index,atoms,0,t,zerobox,xread,NULL);
    
    nframes++;
  } while (read_next_x(status,&t,nat,xread,box));
  close_trj(status);
  if (filterfile)
    close_trx(out);
  sfree(x);
  
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
    sprintf(str,"projection on eigenvector %d",eignr[outvec[0]]+1);
    sprintf(str2,"projection on eigenvector %d",eignr[outvec[noutvec-1]]+1); 
    xvgrout=xvgropen(twodplotfile,"2D projection of trajectory",str,str2);
    for(i=0; i<nframes; i++)
      fprintf(xvgrout,"%10.5f %10.5f\n",inprod[0][i],inprod[noutvec-1][i]);
    fclose(xvgrout);
  }
  if (extremefile) {
    if (extreme==0) {
      min=inprod[0][0];
      max=inprod[0][0];
      for(i=0; i<nframes; i++) {
	if (inprod[0][i]<inprod[0][min])
	  min=i;
	if (inprod[0][i]>inprod[0][max])
	  max=i;
      }
      fprintf(stderr,"\nMinimum along eigenvector %d is %g at t=%g\n",
	      eignr[outvec[0]]+1,inprod[0][min],inprod[noutvec][min]); 
      fprintf(stderr,"Maximum along eigenvector %d is %g at t=%g\n",
	      eignr[outvec[0]]+1,inprod[0][max],inprod[noutvec][max]); 
    }
    else {
      min=-extreme;
      max=+extreme;
    }
    out=open_trx(extremefile,"w");
    for(frame=0; frame<nextr; frame++) {
      if ((extreme==0) && (nextr<=3))
	for(i=0; i<natoms; i++)
	  atoms->atom[index[i]].chain='A'+frame;
      for(i=0; i<natoms; i++)
	for(d=0; d<DIM; d++)
	  xread[index[i]][d] = xav[i][d]+
	    (inprod[0][min]*(nextr-frame-1)+inprod[0][max]*frame)/(nextr-1)
	    *eigvec[outvec[0]][i][d]/sqrtm[i];
      write_trx(out,natoms,index,atoms,0,inprod[noutvec][j],
		zerobox,xread,NULL);
    }
    close_trx(out);
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
}

void read_eigenvectors(char *file,int *natoms,rvec **xref,rvec **xav,
		       int *nvec, int **eignr, rvec ***eigvec)
{
  int status,i,snew_size;
  rvec *x;
  matrix box;
  real t;

  /* read (reference (t=-1) and) average (t=0) structure */
  *natoms=read_first_x(&status,file,&t,xav,box);
  if ((t>=-1.1) && (t<=-0.9)) {
    snew(*xref,*natoms);
    for(i=0; i<*natoms; i++)
      copy_rvec((*xav)[i],(*xref)[i]);
    fprintf(stderr,"\nRead reference structure with %d atoms from %s\n",
	    *natoms,file);
    read_next_x(status,&t,*natoms,*xav,box);
  }
  else
    *xref=NULL;
  if ((t<=-0.01) || (t>=0.01))
    fatal_error(0,"\n%s does not start with t=0, which should be the average "
		"structure. This might not be a eigenvector file.",file);
  fprintf(stderr,"\nRead average structure with %d atoms from %s\n",
	  *natoms,file);
  
  snew(x,*natoms);
  snew_size=0;
  *nvec=0;
  while (read_next_x(status,&t,*natoms,x,box)) {
    if (*nvec >= snew_size) {
      snew_size+=10;
      srenew(*eignr,snew_size);
      srenew(*eigvec,snew_size);
    }
    i=(int)(t+0.01);
    if ((t-i<=-0.01) || (t-i>=0.01))
      fatal_error(0,"\n%s contains a frame with non-integer time (%f), this "
		  "time should be an eigenvector index. "
		  "This might not be a eigenvector file.",file,t);
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
    "if present, otherwise to the structure in the run input file.[PAR]",
    "[TT]-comp[tt]: plot eigenvector components per atom of eigenvectors",
    "[TT]-first[tt] to [TT]-last[tt].[PAR]",
    "[TT]-proj[tt]: calculate projections of a trajectory on eigenvectors",
    "[TT]-first[tt] to [TT]-last[tt].[PAR]",
    "[TT]-2d[tt]: calculate a 2d projection of a trajectory on eigenvectors",
    "[TT]-first[tt] and [TT]-last[tt].[PAR]",
    "[TT]-filt[tt]: filter the trajectory to show only the motion along",
    "eigenvectors [TT]-first[tt] to [TT]-last[tt].[PAR]",
    "[TT]-extr[tt]: calculate the two extreme projections along a trajectory",
    "on the average structure and interpolate between them. The eigenvector",
    "is specified with [TT]-first[tt].[PAR]",
    "[TT]-over[tt]: calculate the cumulative overlap (inner-products)",
    "of the eigenvectors in file [TT]-v2[tt] with eigenvectors",
    "[TT]-first[tt] to [TT]-last[tt].[PAR]",
    "[TT]-inpr[tt]: calculate a matrix of inner-products between two sets",
    "of eigenvectors."
  };
  static bool bM=TRUE;
  static int  first=1,last=10,skip=1,nextr=2;
  static real max=0.0;
  t_pargs pa[] = {
    { "-m",  FALSE, etBOOL, &bM,
      "mass weighted eigenvectors"},
    { "-first", FALSE, etINT, &first,     
      "first eigenvector for analysis" },
    { "-last",  FALSE, etINT, &last, 
      "last eigenvector for analysis (-1 is till the last)" },
     { "-skip",  FALSE, etINT, &skip,
      "only write a filtered frame every nr-th frame" },
    { "-max",  FALSE, etREAL, &max, 
      "maximum for projection of the eigenvector on the average structure, max=0 gives the extremes" },
    { "-nframes",  FALSE, etINT, &nextr, 
      "number of frames for the extremes output" }
  };
  FILE       *out;
  int        status,trjout;
  t_tpxheader  header;
  t_inputrec   ir;
  t_topology   top;
  t_idef       *idef;
  int        nref1,nref2;
  rvec       *xref1,*xref2;
  int        nvec1,nvec2,*eignr1=NULL,*eignr2=NULL;
  rvec       *x,*xread,*xav1,*xav2,**eigvec1=NULL,**eigvec2=NULL;
  matrix     box;
  real       xid,*mass,*sqrtm,avsqrtm,*w_rls,t,lambda;
  int        natoms,ntopatoms,step;
  char       *grpname,*indexfile;
  int        i,j,d;
  int        nout,*iout,noutvec,*outvec,nfit;
  atom_id    *index,*ifit;
  char       *Vec2File,*CompFile,*ProjOnVecFile,*TwoDPlotFile;
  char       *FilterFile,*ExtremeFile;
  char       *OverlapFile,*InpMatFile;
  bool       bIndex,bTop,bVec2,bProj,bFit,bFirstToLast;
  t_filenm fnm[] = { 
    { efTRN, "-v", "eigvec", ffREAD },
    { efTRN, "-v2", "eigvec2", ffOPTRD },
    { efTRX, "-f", NULL, ffOPTRD }, 
    { efTPX, NULL, NULL, ffOPTRD },
    { efNDX, NULL, NULL, ffOPTRD },
    { efXVG, "-comp", "eigcomp", ffOPTWR },
    { efXVG, "-proj", "proj", ffOPTWR },
    { efXVG, "-2d", "2dplot", ffOPTWR },
    { efTRX, "-filt", "filtered", ffOPTWR },
    { efTRX, "-extr", "extreme.pdb", ffOPTWR },
    { efXVG, "-over", "overlap", ffOPTWR },
    { efXPM, "-inpr", "inprod", ffOPTWR }
  }; 
#define NFILE asize(fnm) 

  CopyRight(stderr,argv[0]); 
  parse_common_args(&argc,argv,PCA_CAN_TIME,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL); 

  indexfile=ftp2fn_null(efNDX,NFILE,fnm);

  Vec2File        = opt2fn_null("-v2",NFILE,fnm);
  CompFile        = opt2fn_null("-comp",NFILE,fnm);
  ProjOnVecFile   = opt2fn_null("-proj",NFILE,fnm);
  TwoDPlotFile    = opt2fn_null("-2d",NFILE,fnm);
  FilterFile      = opt2fn_null("-filt",NFILE,fnm);
  ExtremeFile     = opt2fn_null("-extr",NFILE,fnm);
  OverlapFile     = opt2fn_null("-over",NFILE,fnm);
  InpMatFile      = ftp2fn_null(efXPM,NFILE,fnm);
  bProj  = ProjOnVecFile || FilterFile || ExtremeFile || TwoDPlotFile;
  bFirstToLast    = CompFile || ProjOnVecFile || FilterFile || OverlapFile;
  bVec2  = Vec2File || OverlapFile || InpMatFile;
  bM     = bM && (CompFile || bProj);
  bFit   = bProj;
  bIndex = bM || bProj;
  bTop   = ftp2bSet(efTPX,NFILE,fnm) || bM || bFit ||
    FilterFile  || (bIndex && indexfile);

  read_eigenvectors(opt2fn("-v",NFILE,fnm),&natoms,&xref1,&xav1,
		    &nvec1,&eignr1,&eigvec1);
  if (bVec2) {
    read_eigenvectors(Vec2File,&i,&xref2,&xav2,&nvec2,&eignr2,&eigvec2);
    if (i!=natoms)
      fatal_error(0,"Dimensions in the eigenvector files don't match");
  }
  if ((xref1==NULL) && bFit) { 
    bTop=TRUE;
    nref1=0;
  }
  else
    nref1=-1;
    
  nfit=0;
  ifit=NULL;
  w_rls=NULL;
  if (bTop) {
    read_tpxheader(ftp2fn(efTPX,NFILE,fnm),&header);
    ntopatoms=header.natoms;
    snew(xref1,ntopatoms);
    read_tpx(ftp2fn(efTPX,NFILE,fnm),&step,&t,&lambda,&ir,box,
	   &ntopatoms,xref1,NULL,NULL,&top);
    if (nref1==0) {
      nref1=ntopatoms;
      printf("\nNote: the structure in %s should be the same\n"
               "      as the one used for the fit in g_covar\n",
	     ftp2fn(efTPX,NFILE,fnm));
      printf("\nSelect the index group that was used for the least squares fit in g_covar\n",natoms);
      get_index(&top.atoms,indexfile,1,&nfit,&ifit,&grpname);
      snew(w_rls,ntopatoms);
      for(i=0; (i<nfit); i++)
	w_rls[ifit[i]]=top.atoms.atom[ifit[i]].m;
    }
  }

  if (bIndex) {
    printf("\nSelect an index group of %d elements that corresponds to the eigenvectors\n",natoms);
    get_index(&top.atoms,indexfile,1,&i,&index,&grpname);
    if (i!=natoms)
      fatal_error(0,"you selected a group with %d elements instead of %d",
		  i,natoms);
    printf("\n");
  }

  snew(mass,natoms);
  snew(sqrtm,natoms);
  if (bM) {
    avsqrtm=0;
    for(i=0; (i<natoms); i++) {
      mass[i]=top.atoms.atom[index[i]].m;
      sqrtm[i]=sqrt(mass[i]);
      avsqrtm+=sqrtm[i];
    }
    avsqrtm=natoms/avsqrtm;
    for(i=0; (i<natoms); i++) 
      sqrtm[i]=sqrtm[i]*avsqrtm;
  } else
    for(i=0; (i<natoms); i++) {
      mass[i]=1.0;
      sqrtm[i]=1.0;
    }

  if (bVec2) {
    t=0;
    for(i=0; (i<natoms); i++)
      for(d=0;(d<DIM);d++)
	t+=sqr(xav1[i][d]-xav2[i][d])*mass[i];
    fprintf(stderr,"RMSD (without fit) between the two average structures:"
	    " %.3f (nm)\n\n",sqrt(t/natoms));
  }
  
  if (last==-1)
    last=natoms*DIM;
  if (bFirstToLast) {
    /* make an index from first to last */
    nout=last-first+1;
    snew(iout,nout);
    for(i=0; i<nout; i++)
      iout[i]=first-1+i;
  } 
  else {
    nout=2;
    snew(iout,nout);
    iout[0]=first-1;
    iout[1]=last-1;
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

  if (CompFile)
    components(CompFile,natoms,sqrtm,nvec1,eignr1,eigvec1,noutvec,outvec);

  if (bProj)
    project(opt2fn("-f",NFILE,fnm),ProjOnVecFile,TwoDPlotFile,FilterFile,skip,
	    ExtremeFile,max,nextr,&(top.atoms),natoms,index,
	    nref1,xref1,nfit,ifit,w_rls,
	    mass,sqrtm,xav1,nvec1,eignr1,eigvec1,noutvec,outvec);

  if (OverlapFile)
    overlap(OverlapFile,natoms,
	    nvec1,eignr1,eigvec1,nvec2,eignr2,eigvec2,noutvec,outvec);
  
  if (InpMatFile)
    inprod_matrix(InpMatFile,natoms,
		  nvec1,eignr1,eigvec1,nvec2,eignr2,eigvec2);
    
  return 0;
}
  
