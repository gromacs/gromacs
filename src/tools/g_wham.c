/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Gyas ROwers Mature At Cryogenic Speed
 */

/* This line is only for CVS version info */
static char *SRCID_template_c = "$Id$";
static float Temperature=298.0;
#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"


static real Tolerance;

typedef struct {
  int nSkip;
  char Reference[256];
  int nPull;
  int nDim;
  rvec Dims;
  char PullName[4][256];
  double UmbPos[4][3];
  double UmbCons[4][3];
} t_UmbrellaHeader;

typedef struct {
  int nPull;
  int nBin;
  double **Histo;
  double *k;
  double *pos;
  double *z;
  double * N;
} t_UmbrellaWindow;

void read_umbrella_header(FILE * file,t_UmbrellaHeader * header)
{
  char Buffer1[256];
  char Buffer2[256];
  char Buffer3[256];
  int i,j;
	
  fscanf(file,"%s%s",Buffer1,Buffer2);

  if(strcmp(Buffer1,"UMBRELLA")) 
    fatal_error(0,"This does not appear to be a valid pdo file");
  if(strcmp(Buffer2,"3.0"))
    fatal_error(0,"This does not appear to be a version 3.0 pdo file");

  fscanf(file,"%f%f%f",&(header->Dims[0]),&(header->Dims[1]),&(header->Dims[2]));
  fscanf(file,"%d",&(header->nSkip));
  fscanf(file,"%s",header->Reference);
  fscanf(file,"%d%d",&(header->nPull),&(header->nDim));
  if(header->nDim!=1)
    fatal_error(0,"Currently only supports one dimension");
	
  for(i=0;i<header->nPull;++i) {
    fscanf(file,"%s",header->PullName[i]);
    for(j=0;j<header->nDim;++j) {
      fscanf(file,"%lf%lf",&(header->UmbPos[i][j]),&(header->UmbCons[i][j]));
    }
  }
	
  fscanf(file,"%s",Buffer3);
}

void read_umbrella_data(FILE * file, t_UmbrellaHeader * header,
			int bins, double min, double max,
			int fileno, t_UmbrellaWindow * win)
{
  int i,j;
  double temp;
  t_UmbrellaWindow * window;
  window=win+(fileno-1);
  /* Need to alocate memory and set up structure*/
  window->nPull=header->nPull;
  window->nBin=bins;
	
  snew(window->Histo,window->nPull);
  snew(window->z,window->nPull);
  snew(window->k,window->nPull);
  snew(window->pos,window->nPull);
  snew(window->N, window->nPull);

  for(i=0;i<window->nPull;++i) {
    window->z[i]=1;
    snew(window->Histo[i],bins);
    window->k[i]=header->UmbCons[i][0];
    window->pos[i]=header->UmbPos[i][0];
    window->N[i]=0;
  }

  /* Done with setup */
					
  while(!feof(file)) {
    for(i=0;i<header->nPull;++i) {
      if(fscanf(file,"%lf",&temp)) {
	temp+=window->pos[i];
	temp-=min;
	temp/=(max-min);
	temp*=bins;
	temp=floor(temp);
	if((int)temp>=0 && (int)temp<bins) {
	  window->Histo[i][(int)temp]+=1;
	  window->N[i]++;
	}
      }
    }
  }
  for(i=0;i<header->nPull;++i) {
    for(j=0;j<bins;++j) {
      window->Histo[i][j]/=window->N[i];
    }
  }
  /*printf("N %f\b",window->N);*/
}

void calc_profile(double *profile,t_UmbrellaWindow * window, int nWindows, double min, double max)
{
  int i,j,k,l,m,n;
  double num;
  double denom;
  double U=0,temp=0;
  double TOTAL=0;
	
  for(i=0;i<window[0].nBin;++i) {	
    num=denom=0;
    for(j=0;j<nWindows;++j) {
      for(k=0;k<window[j].nPull;++k) {
	temp=(double)(i+0.5)*(max-min)/window[j].nBin+min;
	U=0.5*window[j].k[k]*(window[j].pos[k]-temp)*(window[j].pos[k]-temp);
	num+=window[j].N[k]*window[j].Histo[k][i];
	denom+=window[j].N[k]*exp(- U/(8.314e-3*Temperature))/window[j].z[k];
      }
    }
    profile[i]=num/denom;
    TOTAL+=profile[i];
  }
}

double calc_z(double * profile,t_UmbrellaWindow * window, int nWindows, double min, double max)
{
  int i,j,k;
  double U=0,dist=0;
  double MAX=-1e20;
  double total=0;
  double log_total;
  double temp;
  double crap;
	
  for(i=0;i<nWindows;++i) {
    for(j=0;j<window[i].nPull;++j) {
      total=0;
      for(k=0;k<window[i].nBin;++k) {
	dist=(double)(k+0.5)*(max-min)/window[i].nBin+min;
	U=window[i].k[j]*0.5*(window[i].pos[j]-dist)*(window[i].pos[j]-dist);
	total+=profile[k]*exp(-U/(8.314e-3*Temperature));
      }
      log_total=-log(total);
      temp=fabs(log_total+log(window[i].z[j]));
      if(temp>MAX) MAX=temp;
      window[i].z[j]=total;
    }
  }
  return MAX;
}


int main(int argc,char *argv[])
{
  static char *desc[] = {
    "This is an analysis program that implements the Weighted",
    "Histogram Analysis Method (WHAM).  It is intended to analyze",
    ".pdo files generated by mdrun using umbrella sampling. The",
    "options are as follows:\n",
    "   -o name of the potential of mean force output file\n",
    "   -hist name of the histograms output file\n",
    "   -min the minimum coordinate to construct PMF over.\n",   		
    "   -max like min\n",
    "Note: the program will throw out any data that is outside",
    "of min - max.\n",
    "\n",
    "   -noprof only calculate the min and max\n",
    "   -bins   the number of bins to use in calculation of histograms\n",
    "           and free energy profile\n",
    "   -tol    the calculation will stop once it has converged to\n",
    "           within this value\n",
    "\n",
    "***NOTE: THIS PROGRAM IS EXPERIMENTAL AND MAY CONTAIN BUGS OR",
    "INCOMPLETE FEATURES***"};
  
  static float min=0;
  static float max=0;
  static int bins=100;
  static bool noprof=TRUE;
	
  static int n=1;
  /* Extra arguments - but note how you always get the begin/end
   * options when running the program, without mentioning them here!
   */
  
  t_pargs pa[] = {
    { "-min", FALSE, etREAL, {&min},
      "Minimum coordinate in profile"},
    { "-max", FALSE, etREAL, {&max},
      "Maximum coordinate in profile"},
    { "-bins",FALSE, etINT, {&bins},
      "Number of bins in profile"},
    { "-prof", FALSE, etBOOL, {&noprof},
      "Only calculate min and max"},
    { "-temp", FALSE, etREAL, {&Temperature},
      "Temperature"},
    { "-tol", FALSE, etREAL, {&Tolerance},
      "Tolerance"}
  };
  
  t_filenm fnm[] = {
    { efXVG, "-o", "profile", ffWRITE },     	/* output file for profile */
    { efXVG, "-hist","histo", ffWRITE}		/* output file for histograms */
  };
  
  int i,j,k,l; 
  t_UmbrellaHeader header;
  t_UmbrellaWindow * window=NULL;
  double *profile;
  bool flag=FALSE;
  
#define NFILE asize(fnm)

  Tolerance=0.01;
  CopyRight(stderr,argv[0]);

  /* This is the routine responsible for adding default options,
   * calling the X/motif interface, etc. */
  parse_common_args(&argc,argv,PCA_CAN_VIEW|PCA_NOEXIT_ON_ARGS,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
		    
  /* Check to see what we get for command line arguments */
  if(argc<2) fatal_error(0,"You need to specify a series of pdo files as input");
  for(i=1;i<argc;++i) {
    if(!fexist(argv[i]))
      fatal_error(0,"Unable to open file %s.",argv[i]);
  }

  /* Now we need to process the files */
			
  if(!noprof) {
    min=1e20;
    max=-1e20;
  }
  else {
    snew(window,argc-1);
  }
	
  /* Loop over all files */
  for(i=1;i<argc;++i) {
    /* read in the headers */
    FILE  * file;
    char Buffer[256];
    char * Path=NULL;
		
    if(!(Path=getenv("GMX_PATH_GZIP")))
      sprintf(Buffer,"/bin/gunzip -c < %s",argv[i]);
    else
      sprintf(Buffer,"%s/gunzip -c < %s",Path,argv[i]);

    if((file=popen(Buffer,"r"))==NULL)
      fatal_error(0,"Unable to open file %s",argv[i]);
    printf("Opening file %s.\n",argv[i]);
    read_umbrella_header(file,&header);
		
    if(!noprof) {				/* We only want to calculate the min and max */
      while(!feof(file)) {
	for(j=0;j<header.nPull;++j)	{	/* Loop over all pulled groups */
	  double temp=0;
	  int result=0;
	  result=fscanf(file,"%lf",&temp);
	  if(result==1) {
	    temp+=header.UmbPos[j][0];
	    if(temp<min)
	      min=temp;
	    if(temp>max)
	      max=temp;
	  }
	}
      }
    }
		
    else {	
      read_umbrella_data(file,&header,bins,min,max,i,window);
    }
		
    pclose(file); 
  }
  /* If we're only calculated min and max output that now */
  if(!noprof) {
    printf("Final min: %e\n",min);
    printf("Final max: %e\n",max);
  }
  else {
    /* output histograms */
    FILE * histout;
    FILE * profout;
    int counter=0;
    char buffer[255];
    int iterations=0;
    double conv;
    histout=ffopen(opt2fn("-hist",NFILE,fnm),"w");
 		
    /* Do output here */
    for(l=0;l<bins;++l) {
      fprintf(histout,"%e\t",(double)(l+0.5)/bins*(max-min)+min);
      for(i=1;i<argc;++i) {
	for(j=0;j<window[i-1].nPull;++j) {
	  fprintf(histout,"%e\t",window[i-1].Histo[j][l]);
	}
      }
      fprintf(histout,"\n");
    }
	
  
    /* Calculate profile */
    snew(profile,bins);
	
    do {
      calc_profile(profile,window,argc-1,min,max);
      conv=calc_z(profile, window, argc-1, min, max);
      if(! (iterations%100)) {
	printf("Iteration %d. Converged to %f.\n", iterations, conv);
      }
    } while(conv > Tolerance);
		
    profout=ffopen(opt2fn("-o",NFILE,fnm),"w");
    for(i=0;i<bins;++i)
      fprintf(profout,"%e\t%e\n",(double)(i+0.5)/bins*(max-min)+min,profile[i]);
    ffclose(profout);	
  }
  thanx(stderr);
  
  return 0;
}
