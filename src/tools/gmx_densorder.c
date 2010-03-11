/*
 * $Id: densorder.c,v 0.9
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "tpxio.h"
#include "sysstuff.h"
#include "string2.h"
#include "macros.h"
#include "gstat.h"
#include "xvgr.h"
#include "pbc.h"
#include "futil.h"
#include "index.h"
#include "pdbio.h"
#include "physics.h"
#include "gstat.h"
#include "matio.h"
#include "binsearch.h"
#include "powerspect.h"

enum {methSEL, methBISECT, methFUNCFIT, methNR};

void center_coords(t_atoms *atoms,matrix box,rvec x0[],int axis)
{
  int  i,m;
  real tmass,mm;
  rvec com,shift,box_center;
  
  tmass = 0;
  clear_rvec(com);
  for(i=0; (i<atoms->nr); i++) {
    mm     = atoms->atom[i].m;
    tmass += mm;
    for(m=0; (m<DIM); m++) 
      com[m] += mm*x0[i][m];
  }
  for(m=0; (m<DIM); m++) 
    com[m] /= tmass;
  calc_box_center(ecenterDEF,box,box_center);
  rvec_sub(box_center,com,shift);
  shift[axis] -= box_center[axis];
  
  for(i=0; (i<atoms->nr); i++) 
    rvec_dec(x0[i],shift);
}



/*Code needs bounds checking and array-pointer debugging . 29/4*/

static void density_in_time (const char *fn, const char *field,
			     atom_id **index ,int nndx, int grpn, real binwidth, 
			     int nsttblock, real *****Densdevel, int *xslices, int *yslices, int *zslices,
			     int *tblock, t_topology *top, int ePBC, int axis, bool bCenter,output_env_t oenv)

{
/*  
 * *****Densdevel pointer to array of density values in slices and frame-blocks Densdevel[*nsttblock][*xslices][*yslices][*zslices]
 * Densslice[x][y][z]
 * nsttblock - nr of frames in each time-block 
 * binwidth  widths of normal slices 
 *
 * axis	 - axis direction (normal to slices)
 * nndx - number ot atoms in **index
 * grpn	 - group number in index
 */
	matrix box; /* Box - 3x3 -each step*/
	rvec *x0; /* List of Coord without PBC*/ 
	int natoms, status,i,j,k,n, /* loop indices, checks etc*/
	ax1=0,ax2=0, /* tangent directions */
	framenr=0, /* frame number in trajectory*/
	  slicex, slicey, slicez; /*slice # of x y z position */
	real ***Densslice; /* Density-slice in one frame*/
	real dscale; /*physical scaling factor*/
	real t,x,y,z; /* time and coordinates*/
	*tblock=0;/* blocknr in block average - initialise to 0*/

	/* Axis: X=0, Y=1,Z=2 */
	switch(axis)	
	{
	case 0:
		ax1=YY; ax2=ZZ; /*Surface: YZ*/
		break;
	case 1:
		ax1=ZZ; ax2=XX; /* Surface : XZ*/
		break;
 	case 2:
    		ax1 = XX; ax2 = YY; /* Surface XY*/
    		break;
  	default:
    		gmx_fatal(FARGS,"Invalid axes. Terminating\n");
	}
	
	if( (natoms= read_first_x(oenv,&status,fn,&t,&x0,box)==0))
	  gmx_fatal(FARGS, "Could not read init. coordinates!"); /* Open trajectory for read*/
	

	*zslices=(int)(box[axis][axis]/binwidth+0.5);
	*yslices=(int)(box[ax2][ax2]/binwidth+0.5);
	*xslices=(int)(box[ax1][ax1]/binwidth+0.5);  
	fprintf(stderr,
	"\nDividing the box in %5d x %5d x %5d slices with binw %f along axis %d\n",*xslices,*yslices,*zslices,binwidth,axis );
	
	
	
	/****Start trajectory processing***/
	
	/*Initialize Densdevel*/
	
	*Densdevel=NULL;		
	
	do 	{

	  /*Reset Densslice every nsttblock steps*/
		if   ( framenr % nsttblock==0  ){ 
			snew(Densslice,*xslices);
				for (i=0;i<*xslices;i++) {
					snew(Densslice[i],*yslices);
					for (j=0;j<*yslices;j++){
						snew(Densslice[i][j],*zslices);
					}
				}
	
		/*Allocate Memory to  extra frame in Densdevel -  rather stupid approach:				  	  			*A single frame each time, although only every nsttblock steps.*/
		srenew(*Densdevel,*tblock+1);
				
		}
	
		

		rm_pbc(&(top->idef),ePBC,top->atoms.nr,box,x0,x0);
	
		dscale=(*xslices)*(*yslices)*(*zslices)*AMU/ (box[ax1][ax1]*box[ax2][ax2]*box[axis][axis]*nsttblock*(NANO*NANO*NANO));
		
		if (bCenter)
			center_coords(&top->atoms,box,x0,axis);


		for (j=0;j<nndx;j++) { /*Loop over all atoms in selected index*/
			x=x0[index[grpn][j]][ax1];
			y=x0[index[grpn][j]][ax2];
			z=x0[index[grpn][j]][axis];
			while (x<0)
				x+=box[ax1][ax1];
			while(x>box[ax1][ax1])
				x-=box[ax1][ax1];
			
			while (y<0)
				y+=box[ax2][ax2];
			while(y>box[ax2][ax2])
				y-=box[ax2][ax2];
			
			while (z<0)
				z+=box[axis][axis];
			while(z>box[axis][axis])
				z-=box[axis][axis];
			
			slicex=((int) (x/box[ax1][ax1]* *xslices )) % *xslices;
			slicey=((int) (y/box[ax2][ax2]* *yslices)) % *yslices;
			slicez=((int) (z/box[axis][axis]* *zslices)) % *zslices;
			Densslice[slicex][slicey][slicez]+=(top->atoms.atom[index[grpn][j]].m*dscale);
		
			
		}
			
		framenr++;
	
		if(framenr % nsttblock == 0){
		  /*Implicit incrementation of Densdevel via renewal of Densslice*/
		  /*only every nsttblock steps*/
			(*Densdevel)[*tblock]=Densslice; 			     
			(*tblock)++;
		}	
			
	} while(read_next_x(oenv,status,&t,natoms,x0,box));
		
				
	

	/*Need to sum and divide Densdevel array with constant factor blocknr for each element to average in time.*/
	if(NULL != field){
	  /*Debug-filename and filehandle*/
	  FILE *fldH;
	  char fff[STRLEN];
	  int nnn=1;
	  real ddd;

	  sprintf(fff,"%s%%6.2f%%6.2f\n",pdbformat);
	  fldH=ffopen(field,"w");
	  for(i=0;i<*xslices;i++){
	    for (j=0;j<*yslices;j++){
	      for (k=0;k<*zslices;k++){
		ddd = 0;
		for(n=0;n<*tblock;n++){
		  ddd += (*Densdevel)[n][i][j][k];
		}
		ddd = 0.01*ddd/(*tblock);
		fprintf(fldH,fff,"ATOM",nnn++,"HE","HE",' ',1,'A',4.0*i,4.0*j,4.0*k,1.0,ddd);
	      }
	    }
	  }
	  
	  ffclose(fldH);
	}
	/*Free memory we no longer need.*/
	sfree(x0);
		
}



static void printdensavg( real ****Densmap, real binwidth, int xslices, int yslices, int zslices, int tdim)
{
	int n,i,j,k;
	real totdens=0;

	/*Allocate memory to array*/
	
		for(n=0;n<tdim;n++){
			for(i=0;i<xslices;i++){
				for (j=0;j<yslices;j++){
					for (k=0;k<zslices;k++){
						totdens+=((Densmap[n][i][j][k])/(xslices*yslices*zslices*tdim));
					}
				}
			}
		}
	
	fprintf(stderr,"Total density [kg/m^3]  %8f",totdens);

}


static void interfaces_txy (real ****Densmap, int xslices, int yslices, int zslices,
			    int tblocks, real binwidth,int method, real dens1, real dens2, t_interf ****intf1, t_interf ****intf2,output_env_t oenv)
{

  /*Returns two pointers to 3D arrays of t_interf structs containing (position,thickness) of the interface(s)*/
	FILE *xvg;
	real *zDensavg; /* zDensavg[z]*/
	int i,j,k,n;
	int xysize;
	int  ndx1, ndx2, deltandx, *zperm;
	real densmid, densl, densr, alpha, pos, spread;
	real splitpoint, startpoint, endpoint;
	real *sigma1, *sigma2;
	const real onehalf= 1.00/2.00;
	t_interf ***int1=NULL,***int2=NULL; /*Interface matrices [t][x,y] - last index in row-major order*/
		/*Create int1(t,xy) and int2(t,xy) arrays with correct number of interf_t elements*/
	xysize=xslices*yslices;
	snew(int1,tblocks);
	snew(int2,tblocks);
	for (i=0;i<tblocks;i++){
		snew(int1[i],xysize);
		snew(int2[i],xysize);
		for(j=0;j<xysize;j++){
			snew(int1[i][j],1);
			snew(int2[i][j],1);
			init_interf(int1[i][j]);
			init_interf(int2[i][j]);
		}		
	}	

if(method==methBISECT){
   densmid= onehalf*(dens1+dens2);	
   snew(zperm,zslices);
   for(n=0;n<tblocks;n++){
	for(i=0;i<xslices;i++){
		for (j=0;j<yslices;j++){
		  rangeArray(zperm,zslices); /*reset permutation array to identity*/
		/*Binsearch returns slice-nr where the order param is  <= setpoint sgmid*/
			ndx1=start_binsearch(Densmap[n][i][j],zperm,0,zslices/2-1,densmid,1);
			ndx2=start_binsearch(Densmap[n][i][j],zperm,zslices/2,zslices-1,densmid,-1);
		/* rho_11= Densmap[n][i][j][zperm[ndx1]]
 		 * rho_12 =Densmap[n][i][j][zperm[ndx1+1]] - in worst case might be far off
 		 * deltandx=zperm[ndx1+1]-zperm[ndx+1] 
 * 		 * alpha=(densmid-rho_11)/(rho_12-rho_11) */
			densl= Densmap[n][i][j][zperm[ndx1]];
			densr= Densmap[n][i][j][zperm[ndx1+1]];
			alpha=(densmid-densl)/(densr-densl);
			deltandx=zperm[ndx1+1]-zperm[ndx1];
			printf("Alpha, Deltandx  %f %i\n", alpha,deltandx);
			if(abs(alpha)>1.0 || abs(deltandx)>3){
					pos=zperm[ndx1];
					spread=-1;
			}
			else {
				pos=zperm[ndx1]+alpha*deltandx;
				spread=binwidth*deltandx;
				}
			int1[n][j+(i*yslices)]->Z=(pos+onehalf)*binwidth;
			int1[n][j+(i*yslices)]->t=spread;
		/*Leave the 2nd interface as test at the moment -needs to be made similar to 1st later.	
		 *We can use the same formulation, since alpha should become negative ie: 
			  alpha=(densmid-Densmap[n][i][j][zperm[ndx2]])/
                                (Densmap[n][i][j][zperm[nxd2+1]]-Densmap[n][i][j][zperm[ndx2]]);
                        deltandx=zperm[ndx2+1]-zperm[ndx2];
                        pos=zperm[ndx2]+alpha*deltandx;   */
			int2[n][j+(i*yslices)]->Z=(zperm[ndx2]+onehalf)*binwidth;
			int2[n][j+(i*yslices)]->t=binwidth;
		}
	}
  }				
}	

if(method==methFUNCFIT){
  /*Assume a box divided in 2 along midpoint of z for starters*/
	startpoint=0.0;
	endpoint = binwidth*zslices;
	splitpoint = (startpoint+endpoint)/2.0;
	/*Initial fit proposals*/
	real beginfit1[4] = {dens1, dens2, (splitpoint/2), 0.5};
	real beginfit2[4] = {dens2, dens1, (3*splitpoint/2), 0.5};
	real *fit1=NULL,*fit2=NULL;

	snew(zDensavg,zslices);
	snew(sigma1,zslices);
	snew(sigma2,zslices);

	for(k=0;k<zslices;k++) sigma1[k]=sigma2[k]=1;	
	/*Calculate average density along z - avoid smoothing by using coarse-grained-mesh*/
	for(k=0;k<zslices;k++){
		for(n=0;n<tblocks;n++){
			for (i=0;i<xslices;i++){
				for(j=0;j<yslices;j++){
			zDensavg[k]+=(Densmap[n][i][j][k]/(xslices*yslices*tblocks));
				}
			}
		}
	}

	if(debug){
	  xvg=xvgropen("DensprofileonZ.xvg", "Averaged Densityprofile on Z","z[nm]","Density[kg/m^3]",oenv);
		for(k=0;k<zslices;k++) fprintf(xvg, "%4f.3   %8f.4\n", k*binwidth,zDensavg[k]);
		fclose(xvg);
	}
	
	/*Fit average density in z over whole trajectory to obtain tentative fit-parameters in fit1 and fit2*/
	
	/*Fit 1st half of box*/
	do_lmfit(zslices, zDensavg,sigma1,binwidth,NULL,startpoint,splitpoint,FALSE,effnERF,beginfit1,3,oenv);
	/*Fit 2nd half of box*/
	do_lmfit(zslices ,zDensavg,sigma2,binwidth,NULL,splitpoint,endpoint,FALSE,effnERF,beginfit2,3,oenv);
	
	/*Initialise the const arrays for storing the average fit parameters*/
		const real * const avgfit1=beginfit1;
		const real * const avgfit2=beginfit2;

		
		
		/*Now do fit over each x  y and t slice to get Zint(x,y,t) - loop is very large, we potentially should average over time directly*/
	for(n=0;n<tblocks;n++){
		for(i=0;i<xslices;i++){
			for (j=0;j<yslices;j++){
			  /*Reinitialise fit for each mesh-point*/
				srenew(fit1,4);
				srenew(fit2,4);
				for (k=0;k<4;k++){
					fit1[k]=avgfit1[k];
					fit2[k]=avgfit2[k];
				}
				/*Now fit and store in structures in row-major order int[n][i][j]*/
				do_lmfit(zslices,Densmap[n][i][j],sigma1,binwidth,NULL,startpoint,splitpoint,FALSE,effnERF,fit1,1,oenv);
				int1[n][j+(yslices*i)]->Z=fit1[2];
				int1[n][j+(yslices*i)]->t=fit1[3];
				do_lmfit(zslices,Densmap[n][i][j],sigma2,binwidth,NULL,splitpoint,endpoint, FALSE,effnERF,fit2,2,oenv);
				int2[n][j+(yslices*i)]->Z=fit2[2];
				int2[n][j+(yslices*i)]->t=fit2[3];
			}
		}
	}
}


*intf1=int1;
*intf2=int2;

}

void writesurftoxpms(t_interf ***surf1,t_interf ***surf2, int tblocks,int xbins, int ybins, int zbins, 
		     real bw,char **outfiles,int maplevels ) 
{
  char numbuf[16];
  int n, i, j;
  real **profile1, **profile2;
  real max1, max2, min1, min2, *xticks, *yticks;
  t_rgb lo={0,0,0};
  t_rgb hi={1,1,1};
  FILE *xpmfile1, *xpmfile2;

/*Prepare xpm structures for output*/

/*Allocate memory to tick's and matrices*/
snew (xticks,xbins+1);
snew (yticks,ybins+1);

profile1=mk_matrix(xbins,ybins,FALSE);
profile2=mk_matrix(xbins,ybins,FALSE);
 
for (i=0;i<xbins+1;i++) xticks[i]+=bw;			
for (j=0;j<ybins+1;j++) yticks[j]+=bw;

xpmfile1 = ffopen(outfiles[0],"w");
xpmfile2 = ffopen(outfiles[1],"w");

max1=max2=0.0;
min1=min2=zbins*bw;

for(n=0;n<tblocks;n++){
sprintf(numbuf,"tblock: %4i",n);
/*Filling matrices for inclusion in xpm-files*/
	for(i=0;i<xbins;i++){ 
		for(j=0;j<ybins;j++){	
				profile1[i][j]=(surf1[n][j+ybins*i])->Z;								 
				profile2[i][j]=(surf2[n][j+ybins*i])->Z;
				/*Finding max and min values*/
				if(profile1[i][j]>max1) max1=profile1[i][j];
				if(profile1[i][j]<min1) min1=profile1[i][j];
				if(profile2[i][j]>max2) max2=profile2[i][j];
				if(profile2[i][j]<min2) min2=profile2[i][j];
		}	
	}

	write_xpm(xpmfile1,3,numbuf,"Height","x[nm]","y[nm]",xbins,ybins,xticks,yticks,profile1,min1,max1,lo,hi,&maplevels);
	write_xpm(xpmfile2,3,numbuf,"Height","x[nm]","y[nm]",xbins,ybins,xticks,yticks,profile2,min2,max2,lo,hi,&maplevels);
}	

ffclose(xpmfile1);
ffclose(xpmfile2); 


sfree(profile1);
sfree(profile2);
sfree(xticks);
sfree(yticks);
}

void writeraw(t_interf ***int1, t_interf ***int2, int tblocks,int xbins, int ybins,char **fnms){
	FILE *raw1, *raw2;
    	int i,j,n;
	
	raw1=ffopen(fnms[0],"w");
	raw2=ffopen(fnms[1],"w");
	fprintf(raw1,"#Produced by: %s #\n", command_line());
	fprintf(raw2,"#Produced by: %s #\n", command_line());
	fprintf(raw1,"#Legend\n#TBlock\n#Xbin Ybin Z t\n");
	fprintf(raw2,"#Legend\n#TBlock\n#Xbin Ybin Z t\n");
	for (n=0;n<tblocks;n++){
		fprintf(raw1,"%5d\n",n);
		fprintf(raw2,"%5d\n",n);
		for(i=0;i<xbins;i++){
			for(j=0;j<ybins;j++){
				fprintf(raw1,"%i  %i  %8.5f  %6.4f\n",i,j,(int1[n][j+ybins*i])->Z,(int1[n][j+ybins*i])->t);
				fprintf(raw2,"%i  %i  %8.5f  %6.4f\n",i,j,(int2[n][j+ybins*i])->Z,(int2[n][j+ybins*i])->t);
			}
		}
	}

	ffclose(raw1);
	ffclose(raw2);
}
		

	
int gmx_densorder(int argc,char *argv[])
{
  static const char *desc[] = {
    "A small program to reduce a two-phase density distribution", 
    "along an axis, computed over a MD trajectory",
    "to 2D surfaces fluctuating in time, by a fit to",
    "a functional profile for interfacial densities",
    "A time-averaged spatial representation of the" ,
    "interfaces can be output with the option -tavg" 
  };
  
  /* Extra arguments - but note how you always get the begin/end
   * options when running the program, without mentioning them here!
   */
  
   
  t_topology *top;
  atom_id **index; /* Index list for single group*/
  char       title[STRLEN],**grpname;
  int        ePBC, *ngx;
  static real binwidth=0.2;
  static real dens1=0.00;
  static real dens2=1000.00;
  static int nsttblock=100;
  static int axis= 2;
  static char *axtitle="Z";
  static int ngrps=1;
  int xslices, yslices, zslices, tblock;
  static bool bCenter = FALSE;
  static bool bFourier=FALSE;
  static bool bAvgt  = FALSE;
  static bool bRawOut=FALSE;
  static int nlevels=100;
  /*Densitymap - Densmap[t][x][y][z]*/
  real ****Densmap=NULL;
  /* Surfaces surf[t][surf_x,surf_y]*/
  t_interf ***surf1, ***surf2;

 static const char *meth[methNR+1]={NULL,"bisect","functional",NULL};
 int eMeth;	

  char **outfiles, **rawfiles, **spectra; /* Filenames for xpm-surface maps, rawdata and powerspectra */
  int nfxpm,nfraw, nfspect; /* # files for interface maps and spectra = # interfaces */
 
t_pargs pa[] = {
    { "-tavg", FALSE, etBOOL, {&bAvgt},
      "Plot time averaged interface profile"},
    { "-bw",FALSE,etREAL,{&binwidth},
	"Binwidth of density distribution along axis"},
    {"-axis", FALSE, etSTR,{&axtitle},
	"Axis Direction - X, Y or Z"},
    {"-method",FALSE ,etENUM,{meth},
	"Interface location method"},
    {"-d1", FALSE, etREAL,{&dens1},
      "Bulk density phase 1 (at small z)"},
    {"-d2", FALSE, etREAL,{&dens2},
      "Bulk density phase 2 (at large z)"},
    { "-tblock",FALSE,etINT,{&nsttblock},
	"Number of frames in one time-block average"},
    { "-nlevel", FALSE,etINT, {&nlevels},
	"Number of Height levels in 2D - XPixMaps"}
  };


  t_filenm fnm[] = {
    { efTPX, "-s",  NULL, ffREAD },   /* this is for the topology */
    { efTRX, "-f", NULL, ffREAD },      /* and this for the trajectory */
    { efNDX, "-n", NULL, ffREAD}, /* this is to select groups */
    { efOUT, "-or", NULL,ffOPTWRMULT}, /* This is for writing out the entire information in the t_interf arrays */
    { efXPM, "-o" ,"interface",ffWRMULT},/* This is for writing out the interface meshes - one xpm-file per tblock*/ 
    { efOUT, "-Spect","intfspect",ffOPTWRMULT}, /* This is for the trajectory averaged Fourier-spectra*/
    { efPDB, "-field","densfield",ffOPTWR}		
  };
  
#define NFILE asize(fnm)
  output_env_t oenv;
  
  CopyRight(stderr,argv[0]);

  /* This is the routine responsible for adding default options,
   * calling the X/motif interface, etc. */
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL,&oenv);


eMeth=nenum(meth);	
bFourier=opt2bSet("-Spect",NFILE,fnm);
bRawOut=opt2bSet("-or",NFILE,fnm);
  
top=read_top(ftp2fn(efTPX,NFILE,fnm),&ePBC);
snew(grpname,ngrps);
snew(index,ngrps);
snew(ngx,ngrps);

/* Calculate axis */
  axis = toupper(axtitle[0]) - 'X';

get_index(&top->atoms,ftp2fn_null(efNDX,NFILE,fnm),ngrps,ngx,index,grpname);

 density_in_time(ftp2fn(efTRX,NFILE,fnm),ftp2fn_null(efPDB,NFILE,fnm),
		 index,ngx[0],0,binwidth,nsttblock,&Densmap,
		 &xslices,&yslices,&zslices,&tblock, top,ePBC,axis,bCenter,oenv);

 interfaces_txy(Densmap,xslices,yslices,zslices,tblock,binwidth,eMeth,dens1,dens2,&surf1,&surf2,oenv);

 /*Output surface-xpms*/
nfxpm=opt2fns(&outfiles,"-o",NFILE,fnm);
if(nfxpm!=2){
	gmx_fatal(FARGS,"No or not correct number (2) of output-files: %d",nfxpm);
}
writesurftoxpms(surf1, surf2, tblock,xslices,yslices,zslices,binwidth , outfiles,zslices);


if (debug){
printdensavg(Densmap,binwidth,xslices,yslices,zslices,tblock);
}

	


/*Output raw-data*/
if (bRawOut){
	nfraw=opt2fns(&rawfiles,"-or",NFILE,fnm);
	if(nfraw!=2){
		gmx_fatal(FARGS,"No or not correct number (2) of output-files: %d",nfxpm);
	}
	writeraw(surf1,surf2,tblock,xslices,yslices,rawfiles);
}



if(bFourier){
	nfspect=opt2fns(&spectra,"-Spect",NFILE,fnm);
	if(nfspect!=2){
		gmx_fatal(FARGS,"No or not correct number (2) of output-file-series: %d"
		,nfspect);
	}
	powerspectavg_intf(surf1, surf2, tblock,xslices,yslices,spectra);
}

sfree(Densmap);
sfree(surf1);
sfree(surf2);


  thanx(stderr);
  return 0;
}
