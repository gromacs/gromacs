/*
 * $Id$
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

#include <stdio.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"
#include "names.h"

#ifndef HAVE_STRDUP
#define HAVE_STRDUP
#endif
#include "string2.h"
#include "xvgr.h"


/* enum for energy units */
enum { enSel, en_kJ, en_kCal, en_kT, enNr };
/* enum for type of input files (pdos, tpr, or pullf) */
enum { whamin_unknown, whamin_tpr, whamin_pullxf, whamin_pdo };
/* enum for methods to make profile cyclic/periodic */
enum { enCycl, enCycl_no, enCycl_yes, enCycl_weighted, enCycl_nr};

typedef struct
{
    /* umbrella with pull code of gromacs 4 */
    int npullgrps;      /* nr of pull groups in tpr file         */
    int pull_geometry;  /* such as distance, position            */
    ivec pull_dim;      /* pull dimension with geometry distance */
    int  pull_ndim;     /* nr of pull_dim != 0                   */
    real *k;            /* force constants in tpr file           */
    rvec *init_dist;    /* reference displacements               */
    real *umbInitDist;  /* referebce displacement in umbrella direction */

    /* From here, old pdo stuff */
    int nSkip;
    char Reference[256];
    int nPull;
    int nDim;
    ivec Dims;
    char PullName[4][256];
    double UmbPos[4][3];
    double UmbCons[4][3];
    bool Flipped[4];
} t_UmbrellaHeader;

typedef struct
{
    int nPull;
    int nBin;
    double **Histo,**cum;
    double *k;
    double *pos;
    double *z;
    double * N, *Ntot;
    bool * Flipped;
    double dt;
    bool **bContrib;
} t_UmbrellaWindow;

typedef struct
{
    char *fnTpr,*fnPullf,*fnPdo,*fnPullx;
    bool bTpr,bPullf,bPdo,bPullx;
    int bins,cycl;
    bool verbose,bShift,bAuto,bBoundsOnly;
    bool bFlipProf;
    real tmin, tmax, dt;
    real Temperature,Tolerance;
    int nBootStrap,histBootStrapBlockLength;
    real dtBootStrap,zProfZero,alpha;
    int bsSeed,stepchange;
    bool bHistBootStrap,bWeightedCycl,bHistOutOnly;
    bool bAutobounds,bNoprof;
    real min,max,dz;
    bool bLog;
    int unit;
    real zProf0;
    bool bProf0Set,bs_verbose;
    bool bHistEq, bTab;
    double *tabX,*tabY,tabMin,tabMax,tabDz;
    int tabNbins;
} t_UmbrellaOptions;


/* Return j such that xx[j] <= x < xx[j+1] */
void searchOrderedTable(double xx[], int n, double x, int *j)
{
    int ju,jm,jl;
    int ascending;


    jl=-1;
    ju=n;
    ascending=(xx[n-1] > xx[0]);
    while (ju-jl > 1)
    {
        jm=(ju+jl) >> 1;
        if ((x >= xx[jm]) == ascending)
            jl=jm;
        else
            ju=jm;
    }
    if (x==xx[0]) *j=0;
    else if (x==xx[n-1]) *j=n-2;
    else *j=jl;
}


/* Read and setup tabulated umbrella potential */
void setup_tab(char *fn,t_UmbrellaOptions *opt)
{
    int i,ny,nl;
    double **y;


    printf("Setting up tabulated potential from file %s\n",fn);
    nl=read_xvg(fn,&y,&ny);
    opt->tabNbins=nl;
    if (ny!=2)
        gmx_fatal(FARGS,"Found %d columns in %s. Expected 2.\n",ny,fn);
    opt->tabMin=y[0][0];
    opt->tabMax=y[0][nl-1];
    opt->tabDz=(opt->tabMax-opt->tabMin)/(nl-1);
    if (opt->tabDz<=0)
        gmx_fatal(FARGS,"The tabulated potential in %s must be provided in \n"
                "ascending z-direction",fn);
    for (i=0;i<nl-1;i++)
        if  (fabs(y[0][i+1]-y[0][i]-opt->tabDz) > opt->tabDz/1e6)
            gmx_fatal(FARGS,"z-values in %s are not equally spaced.\n",ny,fn);
    snew(opt->tabY,nl);
    snew(opt->tabX,nl);
    for (i=0;i<nl;i++){
        opt->tabX[i]=y[0][i];
        opt->tabY[i]=y[1][i];
    }
    printf("Found equally spaced tabulated potential from %g to %g, spacing %g\n",
            opt->tabMin,opt->tabMax,opt->tabDz);
}


void read_pdo_header(FILE * file,t_UmbrellaHeader * header, t_UmbrellaOptions *opt)
{
    char Buffer0[256];
    char Buffer1[256];
    char Buffer2[256];
    char Buffer3[256];
    char Buffer4[256];
    int i,j;


    /*  line 1 */
    if(3 != fscanf(file,"%s%s%s",Buffer0,Buffer1,Buffer2))
    {
	gmx_fatal(FARGS,"Error reading header from pdo file");
    }
    if(strcmp(Buffer1,"UMBRELLA"))
        gmx_fatal(FARGS,"This does not appear to be a valid pdo file. Found %s, expected %s",
                Buffer1, "UMBRELLA");
    if(strcmp(Buffer2,"3.0"))
        gmx_fatal(FARGS,"This does not appear to be a version 3.0 pdo file");

    /*  line 2 */
    if(6 != fscanf(file,"%s%s%s%d%d%d",Buffer0,Buffer1,Buffer2,
		   &(header->Dims[0]),&(header->Dims[1]),&(header->Dims[2])))
    { 
	gmx_fatal(FARGS,"Error reading dimensions in header from pdo file");
    }

    /* printf("%d %d %d\n", header->Dims[0],header->Dims[1],header->Dims[2]); */

    header->nDim = header->Dims[0] + header->Dims[1] + header->Dims[2];
    if(header->nDim!=1)
        gmx_fatal(FARGS,"Currently only supports one dimension");

    /* line3 */
    if(3 != fscanf(file,"%s%s%d",Buffer0,Buffer1,&(header->nSkip)))
    { 
	gmx_fatal(FARGS,"Error reading header from pdo file");
    }

    /* line 4 */
    if(4 != fscanf(file,"%s%s%s%s",Buffer0,Buffer1,Buffer2,header->Reference))
    { 
	gmx_fatal(FARGS,"Error reading header from pdo file");
    }

    /* line 5 */
    if(6 != fscanf(file,"%s%s%s%s%s%d",Buffer0,Buffer1,Buffer2,Buffer3,Buffer4,&(header->nPull)))
    { 
	gmx_fatal(FARGS,"Error reading header from pdo file");
    }

    if (opt->verbose)
        printf("Found nPull=%d , nSkip=%d, ref=%s\n",header->nPull,header->nSkip,
                header->Reference);

    for(i=0;i<header->nPull;++i)
    {
      if(4 != fscanf(file,"%s%s%s%s",Buffer0,Buffer1,Buffer2,header->PullName[i]))
      { 
	  gmx_fatal(FARGS,"Error reading header from pdo file");
      }
        if (opt->verbose)
            printf("pullgroup %d, pullname = %s\n",i,header->PullName[i]);
        for(j=0;j<header->nDim;++j)
        {
	  if(6 != fscanf(file,"%s%s%lf%s%s%lf",Buffer0,Buffer1,&(header->UmbPos[i][j]),
			 Buffer2,Buffer3,&(header->UmbCons[i][j])))
          { 
	      gmx_fatal(FARGS,"Error reading header from pdo file");
	  }

            if (opt->bFlipProf)
            {
                /* We want to combine both halves of a profile into one */
                if(header->UmbPos[i][j]<0)
                {
                    header->UmbPos[i][j]= -header->UmbPos[i][j];
                    header->Flipped[i]=TRUE;
                }
            }
            else header->Flipped[i]=FALSE;
            /*printf("%f\t%f\n",header->UmbPos[i][j],header->UmbCons[i][j]);*/
        }
    }

    if(1 != fscanf(file,"%s",Buffer3))
    { 
	gmx_fatal(FARGS,"Error reading header from pdo file");
    }

    if (strcmp(Buffer3,"#####") != 0)
        gmx_fatal(FARGS,"Expected '#####', found %s. Hick.\n",Buffer3);
}


static char *fgets3(FILE *fp,char ptr[],int *len)
{
    char *p;
    int  slen;


    if (fgets(ptr,*len-1,fp) == NULL)
        return NULL;
    p = ptr;
    while ((strchr(ptr,'\n') == NULL) && (!feof(fp)))
    {
        /* This line is longer than len characters, let's increase len! */
        *len += STRLEN;
        p    += STRLEN;
        srenew(ptr,*len);
        if (fgets(p-1,STRLEN,fp) == NULL)
            break;
    }
    slen = strlen(ptr);
    if (ptr[slen-1] == '\n')
        ptr[slen-1] = '\0';

    return ptr;
}


void read_pdo_data(FILE * file, t_UmbrellaHeader * header,
        int fileno, t_UmbrellaWindow * win,
        t_UmbrellaOptions *opt,
        bool bGetMinMax,real *mintmp,real *maxtmp)
{
    int i,inttemp,bins,count;
    real min,max,minfound,maxfound;
    double temp,time,time0=0,dt;
    char *ptr;
    t_UmbrellaWindow * window=0;
    bool timeok,dt_ok=1;
    char  *tmpbuf,fmt[256],fmtign[256];
    int    len=STRLEN,dstep=1;

	minfound=1e20;
	maxfound=-1e20;
	
    if (!bGetMinMax)
    {
        bins=opt->bins;
        min=opt->min;
        max=opt->max;

        window=win+fileno;
        /* Need to alocate memory and set up structure */
        window->nPull=header->nPull;
        window->nBin=bins;

        snew(window->Histo,window->nPull);
        snew(window->z,window->nPull);
        snew(window->k,window->nPull);
        snew(window->pos,window->nPull);
        snew(window->Flipped,window->nPull);
        snew(window->N, window->nPull);
        snew(window->Ntot, window->nPull);

        for(i=0;i<window->nPull;++i)
        {
            window->z[i]=1;
            snew(window->Histo[i],bins);
            window->k[i]=header->UmbCons[i][0];
            window->pos[i]=header->UmbPos[i][0];
            window->Flipped[i]=header->Flipped[i];
            window->N[i]=0;
            window->Ntot[i]=0;
        }
        /* Done with setup */
    }
    else
    {
        minfound=1e20;
        maxfound=-1e20;
        min=max=bins=0; /* Get rid of warnings */
    }

    count=0;
    snew(tmpbuf,len);
    while ( (ptr=fgets3(file,tmpbuf,&len)) != NULL)
    {
        trim(ptr);

        if (ptr[0] == '#' || strlen(ptr)<2)
            continue;

        /* Initiate format string */
        fmtign[0] = '\0';
        strcat(fmtign,"%*s");

        sscanf(ptr,"%lf",&time); /* printf("Time %f\n",time); */
        /* Round time to fs */
        time=1.0/1000*( (int) (time*1000+0.5) );

        /* get time step of pdo file */
        if (count==0)
            time0=time;
        else if (count==1)
        {
            dt=time-time0;
            if (opt->dt>0.0)
            {
	      dstep=(int)(opt->dt/dt+0.5);
                if (dstep==0)
                    dstep=1;
            }
            if (!bGetMinMax)
                window->dt=dt*dstep;
        }
        count++;

        dt_ok=((count-1)%dstep == 0);
        timeok=(dt_ok && time >= opt->tmin && time <= opt->tmax);
        /* if (opt->verbose)
      printf(" time = %f, (tmin,tmax)=(%e,%e), dt_ok=%d timeok=%d\n",
      time,opt->tmin, opt->tmax, dt_ok,timeok); */

        if (timeok)
        {
            for(i=0;i<header->nPull;++i)
            {
                strcpy(fmt,fmtign);
                strcat(fmt,"%lf");      /* Creating a format stings such as "%*s...%*s%lf" */
                strcat(fmtign,"%*s");   /* ignoring one more entry in the next loop */
                if(sscanf(ptr,fmt,&temp))
                {
                    if(opt->bFlipProf)
                    {
                        if(header->Flipped[i]) temp=-temp;
                    }

                    temp+=header->UmbPos[i][0];
                    if (bGetMinMax){
                        if (temp<minfound)
                            minfound=temp;
                        if (temp>maxfound)
                            maxfound=temp;
                    }
                    else
                    {
                        temp-=min;
                        temp/=(max-min);
                        temp*=bins;
                        temp=floor(temp);

                        inttemp = (int)temp;
                        if (opt->cycl==enCycl_yes)
                        {
                            if (inttemp < 0)
                                inttemp+=bins;
                            else if (inttemp >= bins)
                                inttemp-=bins;
                        }

                        if(inttemp >= 0 && inttemp < bins)
                        {
                            window->Histo[i][inttemp]+=1;
                            window->N[i]++;
                        }
                        window->Ntot[i]++;
                    }
                }
            }
        }
        if (time>opt->tmax)
        {
            if (opt->verbose)
                printf("time %f larger than tmax %f, stop reading pdo file\n",time,opt->tmax);
            break;
        }
    }

    if (bGetMinMax)
    {
        *mintmp=minfound;
        *maxtmp=maxfound;
    }
}


void enforceEqualWeights(t_UmbrellaWindow * window,int nWindows)
{
    int i,k,j,NEnforced;
    double ratio;


    NEnforced=window[0].Ntot[0];
    printf("\nFound -hist-eq. Enforcing equal weights for all histograms, \ni.e. doing a "
            "non-weighted histogram analysis method. Ndata = %d\n",NEnforced);
    /* enforce all histograms to have the same weight as the very first histogram */

    for(j=0;j<nWindows;++j)
        for(k=0;k<window[j].nPull;++k)
        {
            ratio=1.0*NEnforced/window[j].Ntot[k];
            for(i=0;i<window[0].nBin;++i)
                window[j].Histo[k][i]*=ratio;
            window[j].N[k]=(int)(ratio*window[j].N[k]+0.5);
        }
}


/* Simple linear interpolation between two given tabulated points */
double tabulated_pot(double dist, t_UmbrellaOptions *opt)
{
    int jl,ju;
    double pl,pu,dz,dp;


    jl=floor((dist-opt->tabMin)/opt->tabDz);
    ju=jl+1;
    if (jl<0 || ju>=opt->tabNbins)
        gmx_fatal(FARGS,"Distance %f out of bounds of tabulated potential (jl=%d, ju=%d).\n"
                "Provide an extended table.",dist,jl,ju);
    pl=opt->tabY[jl];
    pu=opt->tabY[ju];
    dz=dist-opt->tabX[jl];
    dp=(pu-pl)*dz/opt->tabDz;
    return pl+dp;
}


/* Don't worry, that routine does not mean we compute the PMF in limited precision.
   After rapid convergence (using only substiantal contributions), we always switch to
   full precision. */
#define WHAM_CONTRIB_LIM 1e-10
void setup_acc_wham(t_UmbrellaWindow * window,int nWindows, t_UmbrellaOptions *opt)
{
    int i,j,k;
    double U,min=opt->min,dz=opt->dz,temp,ztot_half,distance,ztot,contrib;
    bool bAnyContrib;


    ztot=opt->max-opt->min;
    ztot_half=ztot/2;

    for(i=0;i<nWindows;++i)
    {
        snew(window[i].bContrib,window[i].nPull);
        for(j=0;j<window[i].nPull;++j)
        {
            snew(window[i].bContrib[j],opt->bins);
            bAnyContrib=FALSE;
            for(k=0;k<opt->bins;++k)
            {
                temp=(1.0*k+0.5)*dz+min;
                distance = temp - window[i].pos[j];   /* distance to umbrella center */
                if (opt->cycl==enCycl_yes)
                {                                     /* in cyclic wham:             */
                    if (distance > ztot_half)           /*    |distance| < ztot_half   */
                        distance-=ztot;
                    else if (distance < -ztot_half)
                        distance+=ztot;
                }
                if (!opt->bTab)
                    U=0.5*window[i].k[j]*sqr(distance);       /* harmonic potential assumed. */
                else
                    U=tabulated_pot(distance,opt);            /* Use tabulated potential     */

                contrib=exp(- U/(8.314e-3*opt->Temperature));
                window[i].bContrib[j][k] = (contrib > WHAM_CONTRIB_LIM);
                bAnyContrib = (bAnyContrib | window[i].bContrib[j][k]);
            }
            /* If this histo is far outside min and max all bContrib may be FALSE,
             * causing a floating point exception later on. To avoid that, switch
             * them all to true.*/
            if (!bAnyContrib)
                for(k=0;k<opt->bins;++k)
                    window[i].bContrib[j][k]=TRUE;
        }
    }
    printf("Initialized rapid wham stuff.\n");
}


void calc_profile(double *profile,t_UmbrellaWindow * window, int nWindows, t_UmbrellaOptions *opt,
        bool bExact)
{
    int i,k,j;
    double num,ztot_half,ztot,distance,min=opt->min,dz=opt->dz;
    double denom,U=0,temp=0;


    ztot=opt->max-opt->min;
    ztot_half=ztot/2;


    for(i=0;i<opt->bins;++i)
    {
        num=denom=0;
        for(j=0;j<nWindows;++j)
        {
            for(k=0;k<window[j].nPull;++k)
            {
                temp=(1.0*i+0.5)*dz+min;
                num+=window[j].Histo[k][i];

                if (! (bExact || window[j].bContrib[k][i]))
                    continue;
                distance = temp - window[j].pos[k];   /* distance to umbrella center */
                if (opt->cycl==enCycl_yes){           /* in cyclic wham:             */
                    if (distance > ztot_half)           /*    |distance| < ztot_half   */
                        distance-=ztot;
                    else if (distance < -ztot_half)
                        distance+=ztot;
                }

                if (!opt->bTab)
                    U=0.5*window[j].k[k]*sqr(distance);       /* harmonic potential assumed. */
                else
                    U=tabulated_pot(distance,opt);            /* Use tabulated potential     */
                denom+=window[j].N[k]*exp(- U/(8.314e-3*opt->Temperature) + window[j].z[k]);
            }
        }
        profile[i]=num/denom;
    }
}


double calc_z(double * profile,t_UmbrellaWindow * window, int nWindows, t_UmbrellaOptions *opt,
        bool bExact)
{
    int i,j,k;
    double U=0,min=opt->min,dz=opt->dz,temp,ztot_half,distance,ztot;
    double MAX=-1e20, total=0;


    ztot=opt->max-opt->min;
    ztot_half=ztot/2;

    for(i=0;i<nWindows;++i)
    {
        for(j=0;j<window[i].nPull;++j)
        {
            total=0;
            for(k=0;k<window[i].nBin;++k)
            {
                if (! (bExact || window[i].bContrib[j][k]))
                    continue;
                temp=(1.0*k+0.5)*dz+min;
                distance = temp - window[i].pos[j];   /* distance to umbrella center */
                if (opt->cycl==enCycl_yes)
                {                                     /* in cyclic wham:             */
                    if (distance > ztot_half)           /*    |distance| < ztot_half   */
                        distance-=ztot;
                    else if (distance < -ztot_half)
                        distance+=ztot;
                }

                if (!opt->bTab)
                    U=0.5*window[i].k[j]*sqr(distance);       /* harmonic potential assumed. */
                else
                    U=tabulated_pot(distance,opt);            /* Use tabulated potential     */

                total+=profile[k]*exp(-U/(8.314e-3*opt->Temperature));
            }
            if (total > 0.0)
                total = -log(total);
            else
                total = 1000.0;
            temp = fabs(total - window[i].z[j]);
            if(temp > MAX) MAX=temp;
            window[i].z[j] = total;
        }
    }
    return MAX;
}


void cyclicProfByWeightedCorr(double *profile,t_UmbrellaWindow *window,
        int nWindows, t_UmbrellaOptions * opt,
        bool bAppendCorr2File, char *fn)
{
    int i,j,k,bins=opt->bins;
    static int first=1;
    double *weight,sum=0.,diff,*histsum,*corr,sumCorr=0.,dCorr;
    FILE *fp;
    char buf[256];

    if (first)
    {
        printf("\nEnforcing a periodic profile by a sampling wheighted correction.");
        please_cite(stdout,"Hub2008");
    }

    snew(weight,bins-1);
    snew(histsum,bins);
    snew(corr,bins-1);

    /* generate weights proportional to 1/(n(i)*n(i+1))^alpha
     where n(i) is the total nur of data points in bin i from all histograms */
    for(i=0;i<nWindows;++i)
        for(j=0;j<window[i].nPull;++j)
            for(k=0;k<bins;++k)
                histsum[k]+=window[i].Histo[j][k];

    for(k=0,sum=0.;k<bins-1;++k)
    {
        weight[k]=1./pow(histsum[k]*histsum[k+1],opt->alpha);
        sum+=weight[k];
    }
    for(k=0;k<bins-1;++k)
        weight[k]/=sum;

    /* difference between last and first bin */
    diff=profile[bins-1]-profile[0];
    printf("Distributing %f between adjacent bins to enforce a cyclic profile\n",diff);

    for (i=0;i<bins-1;i++)
    {
        dCorr=weight[i]*diff;
        sumCorr+=dCorr;
        corr[i]=sumCorr;
    }

    for (i=0;i<bins-1;i++)
        profile[i+1]-=corr[i];

    if (bAppendCorr2File)
    {
        fp=xvgropen(fn,"Corrections to enforce periodicity","z","\\f{12}D\\f{}G(z)");
        sprintf(buf,"corrections propotional to 1/(n\\si\\Nn\\si+1\\N)\\S%.2f",opt->alpha);
        xvgr_subtitle(fp,buf);
        for (i=0;i<bins-1;i++)
            fprintf(fp,"%g %g\n",opt->min+opt->dz*(i+1),-corr[i]);
        fclose(fp);
    }

    sfree(histsum);
    sfree(corr);
    sfree(weight);
    first=0;
}


void prof_normalization_and_unit(double * profile, t_UmbrellaOptions *opt)
{
    int i,bins,imin;
    double unit_factor=1., R_MolarGasConst, diff;


    R_MolarGasConst=8.314472e-3; /* in kJ/(mol*K) */
    bins=opt->bins;

    /* No log? Nothing to do! */
    if (!opt->bLog)
        return;

    /* Get profile in units of kT, kJ/mol, or kCal/mol */
    if (opt->unit == en_kT)
        unit_factor=1.0;
    else if (opt->unit == en_kJ)
        unit_factor=R_MolarGasConst*opt->Temperature;
    else if (opt->unit == en_kCal)
        unit_factor=R_MolarGasConst*opt->Temperature/4.1868;
    else
        gmx_fatal(FARGS,"Sorry I don't know this energy unit.");

    for (i=0;i<bins;i++)
        if (profile[i]>0.0)
            profile[i]=-log(profile[i])*unit_factor;

    /* shift to zero at z=opt->zProf0 */
    if (!opt->bProf0Set)
        diff=profile[0];
    else{
        /* Get bin with shortest distance to opt->zProf0 */
      imin=(int)((opt->zProf0-opt->min)/opt->dz);
        if (imin<0)
            imin=0;
        else if (imin>=bins)
            imin=bins-1;
        diff=profile[imin];
    }

    /* Shift to zero */
    for (i=0;i<bins;i++)
        profile[i]-=diff;
}


void getRandomIntArray(int nPull,int blockLength,int* randomArray)
{
    int ipull,blockBase,nr,ipullRandom;


    if (blockLength==0)
        blockLength=nPull;

    for (ipull=0; ipull<nPull; ipull++)
    {
        blockBase=(ipull/blockLength)*blockLength;
        do{      /* make sure nothing bad happens in the last block */
            nr=(int)((1.0*rand()/RAND_MAX)*blockLength);
            ipullRandom = blockBase + nr;
        } while (ipullRandom >= nPull);
        if (ipullRandom<0 || ipullRandom>=nPull)
            gmx_fatal(FARGS,"Ups, random iWin = %d, nPull = %d, nr = %d, "
                    "blockLength = %d, blockBase = %d\n",
                    ipullRandom,nPull,nr,blockLength,blockBase);
        randomArray[ipull]=ipullRandom;
    }
    /*for (ipull=0; ipull<nPull; ipull++)
    printf("%d ",randomArray[ipull]); printf("\n"); */
}


void copy_pullgrp_to_synthwindow(t_UmbrellaWindow *synthWindow,
        t_UmbrellaWindow *thisWindow,int pullid)
{
    synthWindow->N       [0]=thisWindow->N        [pullid];
    synthWindow->Histo   [0]=thisWindow->Histo    [pullid];
    synthWindow->pos     [0]=thisWindow->pos      [pullid];
    synthWindow->z       [0]=thisWindow->z        [pullid];
    synthWindow->k       [0]=thisWindow->k        [pullid];
    synthWindow->bContrib[0]=thisWindow->bContrib [pullid];
}


/* Calculate cummulative of all histograms. They allow to create random numbers
   which are distributed according to the histograms. Required to generate
   the "synthetic" histograms for the Bootstrap method */
void calc_cummulants(t_UmbrellaWindow *window,int nWindows,
        t_UmbrellaOptions *opt,char *fnhist)
{
    int i,j,k,nbin;
    double last;
    char *fn=0,*buf=0;
    FILE *fp=0;


    if (opt->bs_verbose)
    {
        snew(fn,strlen(fnhist)+10);
        snew(buf,strlen(fnhist)+10);
        sprintf(fn,"%s_cumul.xvg",strncpy(buf,fnhist,strlen(fnhist)-4));
        fp=xvgropen(fn,"Cummulants of umbrella histograms","z","cummulant");
    }

    nbin=opt->bins;
    for (i=0; i<nWindows; i++)
    {
        snew(window[i].cum,window[i].nPull);
        for (j=0; j<window[i].nPull; j++)
        {
            snew(window[i].cum[j],nbin+1);
            window[i].cum[j][0]=0.;
            for (k=1; k<=nbin; k++)
                window[i].cum[j][k] = window[i].cum[j][k-1]+window[i].Histo[j][k-1];

            /* normalize cummulant. Ensure cum[nbin]==1 */
            last = window[i].cum[j][nbin];
            for (k=0; k<=nbin; k++)
                window[i].cum[j][k] /= last;
        }
    }

    printf("Cumulants of all histograms created.\n");
    if (opt->bs_verbose)
    {
        for (k=0; k<=nbin; k++)
        {
            fprintf(fp,"%g\t",opt->min+k*opt->dz);
            for (i=0; i<nWindows; i++)
                for (j=0; j<window[i].nPull; j++)
                    fprintf(fp,"%g\t",window[i].cum[j][k]);
            fprintf(fp,"\n");
        }
        printf("Wrote cumulants to %s\n",fn);
        ffclose(fp);
        sfree(fn);
        sfree(buf);
    }
}


/* Return j such that xx[j] <= x < xx[j+1] */
void searchCummulant(double xx[], int n, double x, int *j)
{
    int ju,jm,jl;


    jl=-1;
    ju=n;
    while (ju-jl > 1)
    {
        jm=(ju+jl) >> 1;
        if (x >= xx[jm])
            jl=jm;
        else
            ju=jm;
    }
    if (x==xx[0])
        *j=0;
    else if (x==xx[n-1])
        *j=n-2;
    else
        *j=jl;
}


void create_synthetic_histo(t_UmbrellaWindow *synthWindow, t_UmbrellaWindow *thisWindow,
        int pullid,t_UmbrellaOptions *opt)
{
    int nsynth,N,i,nbins,r_index;
    double r;
    static bool bWarnout=0;


    N=thisWindow->N[pullid];
    nbins=thisWindow->nBin;

    /* nsynth = nr of data points in synthetic histo */
    if (opt->dtBootStrap==0.0)
        nsynth=N;
    else
    {
      nsynth=(int)(thisWindow->N[pullid]*thisWindow->dt/opt->dtBootStrap+0.5);
        if (nsynth>N)
            nsynth=N;
    }

    if (!bWarnout && nsynth<10)
    {
        printf("\n++++ WARNING ++++\n\tOnly %d data points per synthetic histogram!\n"
                "\tYou may want to consider option -bs-dt.\n\n",nsynth);
        bWarnout=1;
    }

    synthWindow->N       [0]=nsynth;
    synthWindow->pos     [0]=thisWindow->pos[pullid];
    synthWindow->z       [0]=thisWindow->z[pullid];
    synthWindow->k       [0]=thisWindow->k[pullid];
    synthWindow->bContrib[0]=thisWindow->bContrib[pullid];

    for (i=0;i<nbins;i++)
        synthWindow->Histo[0][i]=0.;

    for (i=0;i<nsynth;i++)
    {
        r=1.0*rand()/RAND_MAX;
        searchCummulant(thisWindow->cum[pullid],nbins+1 ,r,&r_index);
        synthWindow->Histo[0][r_index]+=1.;
    }
}


void print_histograms(char *fnhist, t_UmbrellaWindow * window, int nWindows,
        int bs_index,t_UmbrellaOptions *opt)
{
    char *fn,*buf=0,title[256];
    FILE *fp;
    int bins,l,i,j;


    if (bs_index<0)
    {
        fn=fnhist;
        strcpy(title,"Umbrella histograms");
    }
    else
    {
        snew(fn,strlen(fnhist)+6);
        snew(buf,strlen(fnhist)+1);
        sprintf(fn,"%s_bs%d.xvg",strncpy(buf,fnhist,strlen(fnhist)-4),bs_index);
        sprintf(title,"Umbrella histograms. Bootstrap #%d",bs_index);
    }

    fp=xvgropen(fn,title,"z","count");
    bins=opt->bins;

    /* Write histograms */
    for(l=0;l<bins;++l)
    {
        fprintf(fp,"%e\t",(double)(l+0.5)*opt->dz+opt->min);
        for(i=0;i<nWindows;++i)
        {
            for(j=0;j<window[i].nPull;++j)
            {
                fprintf(fp,"%e\t",window[i].Histo[j][l]);
            }
        }
        fprintf(fp,"\n");
    }

    ffclose(fp);
    printf("Wrote %s\n",fn);
    if (buf)
    {
        sfree(buf);
        sfree(fn);
    }
}


void do_bootstrapping(char *fnres, char* fnprof, char *fnhist,
        char* ylabel, double *profile,
        t_UmbrellaWindow * window, int nWindows, t_UmbrellaOptions *opt)
{
    t_UmbrellaWindow * synthWindow;
    double *bsProfile,*bsProfiles_av, *bsProfiles_av2,maxchange=1e20,tmp,stddev;
    int i,j,*randomArray=0,winid,pullid,ib;
    int iAllPull,nAllPull,*allPull_winId,*allPull_pullId;
    FILE *fp;
    bool bExact=FALSE;


    /* init random */
    if (opt->bsSeed==-1)
        srand(time(NULL));
    else
        srand(opt->bsSeed);

    snew(bsProfile,     opt->bins);
    snew(bsProfiles_av, opt->bins);
    snew(bsProfiles_av2,opt->bins);

    /* Create array of all pull groups. Note that different windows
     may have different nr of pull groups
     First: Get total nr of pull groups */
    nAllPull=0;
    for (i=0;i<nWindows;i++)
        nAllPull+=window[i].nPull;
    snew(allPull_winId,nAllPull);
    snew(allPull_pullId,nAllPull);
    iAllPull=0;
    /* Setup one array of all pull groups */
    for (i=0;i<nWindows;i++)
        for (j=0;j<window[i].nPull;j++)
        {
            allPull_winId[iAllPull]=i;
            allPull_pullId[iAllPull]=j;
            iAllPull++;
        }

    /* setup stuff for synthetic windows */
    snew(synthWindow,nAllPull);
    for (i=0;i<nAllPull;i++)
    {
        synthWindow[i].nPull=1;
        synthWindow[i].nBin=opt->bins;
        snew(synthWindow[i].Histo,1);
        if (!opt->bHistBootStrap)
            snew(synthWindow[i].Histo[0],opt->bins);
        snew(synthWindow[i].N,1);
        snew(synthWindow[i].pos,1);
        snew(synthWindow[i].z,1);
        snew(synthWindow[i].k,1);
        snew(synthWindow[i].bContrib,1);
    }

    if (opt->bHistBootStrap)
    {
        snew(randomArray,nAllPull);
        printf("\n\nWhen computing statistical errors by bootstrapping entire histograms:\n");
        please_cite(stdout,"Hub2006");
    }
    else
    {
        calc_cummulants(window,nWindows,opt,fnhist);
    }

    /* do bootstrapping */
    fp=xvgropen(fnprof,"Boot strap profiles","z",ylabel);
    for (ib=0;ib<opt->nBootStrap;ib++)
    {
        printf("  *******************************************\n"
                "  ******** Start bootstrap nr %d ************\n"
                "  *******************************************\n",ib+1);

        if (opt->bHistBootStrap)
        {
            /* only mix given histos */
            getRandomIntArray(nAllPull,opt->histBootStrapBlockLength,randomArray);
            for (i=0;i<nAllPull;i++)
            {
                winid =allPull_winId [randomArray[i]];
                pullid=allPull_pullId[randomArray[i]];
                copy_pullgrp_to_synthwindow(synthWindow+i,window+winid,pullid);
            }
        }
        else
        {
            /* create new histos from given histos */
            for (i=0;i<nAllPull;i++)
            {
                winid=allPull_winId[i];
                pullid=allPull_pullId[i];
                create_synthetic_histo(synthWindow+i,window+winid,pullid,opt);
            }
        }

        /* print histos in case of verbose output */
        if (opt->bs_verbose)
            print_histograms(fnhist,synthWindow,nAllPull,ib,opt);

        /* do wham */
        i=0;
        bExact=FALSE;
        maxchange=1e20;
        memcpy(bsProfile,profile,opt->bins*sizeof(double)); /* use profile as guess */
        do {
            if (maxchange<opt->Tolerance)
                bExact=TRUE;
            if (((i%opt->stepchange) == 0 || i==1) && !i==0)
                printf("\t%4d) Maximum change %e\n",i,maxchange);
            calc_profile(bsProfile,synthWindow,nAllPull,opt,bExact);
            i++;
        } while( (maxchange=calc_z(bsProfile, synthWindow, nAllPull, opt,bExact)) > opt->Tolerance || !bExact);
        printf("\tConverged in %d iterations. Final maximum change %g\n",i,maxchange);

        if (opt->bLog)
            prof_normalization_and_unit(bsProfile,opt);
        /* Force cyclic profile by wheighted correction */
        if (opt->cycl==enCycl_weighted)
            cyclicProfByWeightedCorr(bsProfile,synthWindow,nAllPull,opt, FALSE, 0);

        /* save stuff to get average and stddev */
        for (i=0;i<opt->bins;i++)
        {
            tmp=bsProfile[i];
            bsProfiles_av[i]+=tmp;
            bsProfiles_av2[i]+=tmp*tmp;
            fprintf(fp,"%e\t%e\n",(i+0.5)*opt->dz+opt->min,tmp);
        }
        fprintf(fp,"&\n");
    }
    ffclose(fp);

    /* write average and stddev */
    fp=ffopen(fnres,"w");
    fprintf(fp,"@    title \"%s\"\n","Average and stddev from bootstrapping");
    fprintf(fp,"@    xaxis  label \"%s\"\n","z");
    fprintf(fp,"@    yaxis  label \"%s\"\n",ylabel);
    fprintf(fp,"@TYPE xydy\n");
    for (i=0;i<opt->bins;i++)
    {
        bsProfiles_av [i]/=opt->nBootStrap;
        bsProfiles_av2[i]/=opt->nBootStrap;
        tmp=bsProfiles_av2[i]-sqr(bsProfiles_av[i]);
        stddev=(tmp>=0.) ? sqrt(tmp) : 0.; /* Catch rouding errors */
        fprintf(fp,"%e\t%e\t%e\n",(i+0.5)*opt->dz+opt->min,bsProfiles_av [i],stddev);
    }
    ffclose(fp);
    printf("Wrote boot strap result to %s\n",fnres);
}


int whaminFileType(char *fn)
{
    int len;
    len=strlen(fn);
    if (strcmp(fn+len-3,"tpr")==0)
        return whamin_tpr;
    else if (strcmp(fn+len-3,"xvg")==0 || strcmp(fn+len-6,"xvg.gz")==0)
        return whamin_pullxf;
    else if (strcmp(fn+len-3,"pdo")==0 || strcmp(fn+len-6,"pdo.gz")==0)
        return whamin_pdo;
    else
        gmx_fatal(FARGS,"Unknown file type of %s. Should be tpr, xvg, or pdo.\n",fn);
    return whamin_unknown;
}


void read_wham_in(char *fn,char ***filenamesRet, int *nfilesRet,
        t_UmbrellaOptions *opt)
{
    char **filename,tmp[STRLEN];
    int nread,sizenow,i,block=10;
    FILE *fp;
#define MAXFILELEN 512


    fp=ffopen(fn,"r");
    sizenow=block;
    snew(filename,sizenow);
    for (i=0;i<sizenow;i++)
        snew(filename[i],MAXFILELEN);
    nread=0;
    while (fscanf(fp,"%s",tmp) != EOF)
    {
        if (strlen(tmp)>=MAXFILELEN)
            gmx_fatal(FARGS,"Filename too long. Only %d characters allowed\n",MAXFILELEN);
        strcpy(filename[nread],tmp);
        if (opt->verbose)
            printf("Found file %s in %s\n",filename[nread],fn);
        nread++;
        if (nread>=sizenow)
        {
            sizenow+=block;
            srenew(filename,sizenow);
            for (i=sizenow-block;i<sizenow;i++)
                snew(filename[i],MAXFILELEN);
        }
    }
    *filenamesRet=filename;
    *nfilesRet=nread;
}


FILE *pdo_open_file(char *fn)
{
    char Buffer[1024],gunzip[1024],*Path=0;
    FILE *fp;

    if (!fexist(fn))
	{
        gmx_fatal(FARGS,"File %s does not exist.\n",fn);
	}
	
    /* gzipped pdo file? */
    if (strcmp(fn+strlen(fn)-3,".gz")==0)
    {
#ifdef HAVE_PIPES
        if(!(Path=getenv("GMX_PATH_GZIP")))
            sprintf(gunzip,"%s","/bin/gunzip");
        else
            sprintf(gunzip,"%s/gunzip",Path);
        if (!fexist(gunzip))
            gmx_fatal(FARGS,"Cannot find executable %s. You may want to define the path to gunzip "
                    "with the environment variable GMX_PATH_GZIP.",gunzip);
        sprintf(Buffer,"%s -c < %s",gunzip,fn);
		if((fp=popen(Buffer,"r"))==NULL)
		{
			gmx_fatal(FARGS,"Unable to open pipe to `%s'\n",Buffer);
		}
#else
		gmx_fatal(FARGS,"Cannot open a compressed file on platform without pipe support");
#endif
    }
    else
	{
		if((fp=fopen(fn,"r"))==NULL)
		{
			gmx_fatal(FARGS,"Unable to open file %s\n",fn);
		}		
	}
	return fp;
}

void
pdo_close_file(FILE *fp)
{
#ifdef HAVE_PIPES
	pclose(fp);
#else
	fclose(fp);
#endif
}

/* Reading pdo files */
void read_pdo_files(char **fn, int nfiles, t_UmbrellaHeader* header,
        t_UmbrellaWindow **window, t_UmbrellaOptions *opt)
{
    FILE  * file;
    real mintmp,maxtmp;
    int i;


    if(nfiles<1)
        gmx_fatal(FARGS,"No files found. Hick.");

    /* if min and max are not given, get min and max from the input files */
    if (opt->bAuto)
    {
        printf("Automatic determination of boundaries from %d pdo files...\n",nfiles);
        opt->min=1e20;
        opt->max=-1e20;
        for(i=0;i<nfiles;++i)
        {
            file=pdo_open_file(fn[i]);
            printf("\rOpening %s ...",fn[i]); fflush(stdout);
            if (opt->verbose)
                printf("\n");
            read_pdo_header(file,header,opt);
            /* here only determine min and max of this window */
            read_pdo_data(file,header,i,NULL,opt,TRUE,&mintmp,&maxtmp);
            if (maxtmp>opt->max)
                opt->max=maxtmp;
            if (mintmp<opt->min)
                opt->min=mintmp;
            pdo_close_file(file);
        }
        printf("\n");
        printf("\nDetermined boundaries to %f and %f\n\n",opt->min,opt->max);
        if (opt->bBoundsOnly)
        {
            printf("Found option -boundsonly, now exiting.\n");
            exit (0);
        }
    }
    /* store stepsize in profile */
    opt->dz=(opt->max-opt->min)/opt->bins;

    snew(*window,nfiles);

    /* Having min and max, we read in all files */
    /* Loop over all files */
    for(i=0;i<nfiles;++i)
    {
        printf("\rOpening %s ...",fn[i]); fflush(stdout);
        if (opt->verbose)
            printf("\n");
        file=pdo_open_file(fn[i]);
        /* read in the headers */
        read_pdo_header(file,header,opt);
        /* load data into window */
        read_pdo_data(file,header,i,*window,opt,FALSE,NULL,NULL);
        pdo_close_file(file);
    }
    printf("\n");
}


#define int2YN(a) (((a)==0)?("N"):("Y"))

void read_tpr_header(char *fn,t_UmbrellaHeader* header,t_UmbrellaOptions *opt)
{
    t_inputrec  ir;
    int         i,ngrp,d;
    t_state     state;
    static int first=1;


    /* printf("Reading %s \n",fn); */
    read_tpx_state(fn,&ir,&state,NULL,NULL);

    if (ir.ePull != epullUMBRELLA)
        gmx_fatal(FARGS,"This is not a tpr of an umbrella simulation. Found ir.ePull = %s\n",
                epullg_names[ir.ePull]);

    /* nr of pull groups */
    ngrp=ir.pull->ngrp;
    if (ngrp < 1)
        gmx_fatal(FARGS,"This is not a tpr of umbrella simulation. Found only %d pull groups\n",ngrp);

    header->npullgrps=ir.pull->ngrp;
    header->pull_geometry=ir.pull->eGeom;
    copy_ivec(ir.pull->dim,header->pull_dim);
    header->pull_ndim=header->pull_dim[0]+header->pull_dim[1]+header->pull_dim[2];
    if (header->pull_geometry==epullgPOS && header->pull_ndim>1)
        gmx_fatal(FARGS,"Found pull geometry 'position' and more than 1 pull dimension (%d).\n"
                "Hence, the pull potential does not correspond to a one-dimensional umbrella potential.\n"
                "If you have some special umbrella setup you may want to write your own pdo files\n"
                "and feed them into g_wham. Check g_wham -h !\n",header->pull_ndim);
    snew(header->k,ngrp);
    snew(header->init_dist,ngrp);
    snew(header->umbInitDist,ngrp);

    for (i=0;i<ngrp;i++)
    {
        header->k[i]=ir.pull->grp[i+1].k;
        if (header->k[i]==0.0)
            gmx_fatal(FARGS,"Pull group %d has force constant of of 0.0 in %s.\n"
                    "That doesn't seem to be an Umbrella tpr.\n",
                    i,fn);
        copy_rvec(ir.pull->grp[i+1].init,header->init_dist[i]);
        header->Flipped[i]=opt->bFlipProf;

        /* initial distance to reference */
        switch(header->pull_geometry)
        {
        case epullgPOS:
            for (d=0;d<DIM;d++)
                if (header->pull_dim[d])
                    header->umbInitDist[i]=header->init_dist[i][d];
            break;
        case epullgDIST:
            header->umbInitDist[i]=header->init_dist[i][0];
            break;
        default:
            gmx_fatal(FARGS,"Pull geometry %s not supported\n",epullg_names[header->pull_geometry]);
        }
    }

    if (opt->verbose || first)
    {
        printf("File %s, %d groups, geometry \"%s\", dimensions [%s %s %s], (%d dimensions)\n",
                fn,header->npullgrps,epullg_names[header->pull_geometry],
                int2YN(header->pull_dim[0]),int2YN(header->pull_dim[1]),int2YN(header->pull_dim[2]),
                header->pull_ndim);
        for (i=0;i<ngrp;i++)
            printf("\tgrp %d) k = %.3f  inittial distance = %g\n",i,header->k[i],header->umbInitDist[i]);
    }
    first=0;
}


double dist_ndim(double **dx,int ndim,int line)
{
    int i;
    double r2=0.;


    for (i=0;i<ndim;i++)
        r2+=sqr(dx[i][line]);

    return sqrt(r2);
}


void read_pull_xf(char *fn, char *fntpr, t_UmbrellaHeader * header,
        t_UmbrellaWindow * window,
        t_UmbrellaOptions *opt,
        bool bGetMinMax,real *mintmp,real *maxtmp)
{
    double **y,pos=0.,t,force,time0=0.,dt;
    int ny,nt,bins,ibin,i,g,dstep=1,nColPerGrp,nColRef,nColExpect;
    real min,max,minfound,maxfound;
    bool dt_ok,timeok,bHaveForce;
    char *quantity;

	minfound=1e20;
	maxfound=-1e20;

    /*
     in force    output pullf.xvg: No   reference, one  column  per pull group
     in position output pullx.xvg: ndim reference, ndim columns per pull group
     */
    nColPerGrp = opt->bPullx ? header->pull_ndim : 1;
    nColRef    = opt->bPullx ? header->pull_ndim : 0;
    quantity   = opt->bPullx ? "position" : "force";
    nColExpect = 1 + nColRef + header->npullgrps*nColPerGrp;
    bHaveForce = opt->bPullf;

    nt=read_xvg(fn,&y,&ny);

    /* Check consistency */
    if (nt<1)
        gmx_fatal(FARGS,"Empty pull %s file %s\n",quantity,fn);
    if (ny != nColExpect)
        gmx_fatal(FARGS,"Found %d pull groups in %s,\n but %d data columns in %s (expected %d)\n",
                header->npullgrps,fntpr,ny-1,fn,nColExpect-1);

    if (opt->verbose)
        printf("Found %d times and %d %s sets %s\n",nt,(ny-1)/nColPerGrp,quantity,fn);

    if (!bGetMinMax)
    {
        bins=opt->bins;
        min=opt->min;
        max=opt->max;
        if (nt>1)
            window->dt=y[0][1]-y[0][0];
        else if (opt->nBootStrap && opt->dtBootStrap!=0.0)
            fprintf(stderr,"\n *** WARNING, Could not determine time step in %s\n",fn);

        /* Need to alocate memory and set up structure */
        window->nPull=header->npullgrps;
        window->nBin=bins;

        snew(window->Histo,window->nPull);
        snew(window->z,window->nPull);
        snew(window->k,window->nPull);
        snew(window->pos,window->nPull);
        snew(window->Flipped,window->nPull);
        snew(window->N, window->nPull);
        snew(window->Ntot, window->nPull);

        for(g=0;g<window->nPull;++g)
        {
            window->z[g]=1;
            snew(window->Histo[g],bins);
            window->k[g]=header->k[g];
            window->Flipped[g]=header->Flipped[g];
            window->N[g]=0;
            window->Ntot[g]=0;
            window->pos[g]=header->umbInitDist[g];
        }
    }
    else
    {
        /* only determine min and max */
        minfound=1e20;
        maxfound=-1e20;
        min=max=bins=0; /* Get rid of warnings */
    }

    if(header->Flipped[0])
        gmx_fatal(FARGS,"Sorry, flipping not supported for gmx4 output\n");

    for (i=0;i<nt;i++)
    {
        /* Do you want that time frame? */
      t=1.0/1000*((int)(0.5+y[0][i]*1000)); /* round time to fs */

        /* get time step of pdo file and get dstep from opt->dt */
        if (i==0)
            time0=t;
        else if (i==1)
        {
            dt=t-time0;
            if (opt->dt>0.0)
            {
                dstep=(int)(opt->dt/dt+0.5);
                if (dstep==0)
                    dstep=1;
            }
            if (!bGetMinMax)
                window->dt=dt*dstep;
        }

        dt_ok=(i%dstep == 0);
        timeok=(dt_ok && t >= opt->tmin && t <= opt->tmax);
        /*if (opt->verbose)
      printf(" time = %f, (tmin,tmax)=(%e,%e), dt_ok=%d timeok=%d\n",
      t,opt->tmin, opt->tmax, dt_ok,timeok); */

        if (timeok)
        {
            for(g=0;g<header->npullgrps;++g)
            {
                if (bHaveForce)
                {
                    /* y has 1 time column y[0] and one column per force y[1],...,y[nGrps] */
                    force=y[g+1][i];
                    pos= -force/header->k[g] + header->umbInitDist[g];
                }
                else
                {
                    switch (header->pull_geometry)
                    {
                    case epullgDIST:
                        /* y has 1 time column y[0] and nColPerGrps columns per pull group;
	                       Distance to reference: */
                        pos=dist_ndim(y+1+nColRef+g*nColPerGrp,header->pull_ndim,i);
                        break;
                    case epullgPOS:
                        /* with geometry==position, we have always one column per group;
	                       Distance to reference: */
                        pos=y[1+nColRef+g][i];
                        break;
                    default:
                        gmx_fatal(FARGS,"Bad error, this error should have been catched before. Ups.\n");
                    }
                }

                /* printf("grp %d dpos %f poseq %f pos %f \n",g,dpos,poseq,pos); */
                if (bGetMinMax)
                {
                    if (pos<minfound)
                        minfound=pos;
                    if (pos>maxfound)
                        maxfound=pos;
                }
                else
                {
                    ibin=(int) floor((pos-min)/(max-min)*bins);
                    if (opt->cycl==enCycl_yes)
                    {
                        if (ibin<0)
                            while ( (ibin+=bins) < 0);
                        else if (ibin>=bins)
                            while ( (ibin-=bins) >= bins);
                    }
                    if(ibin >= 0 && ibin < bins)
                    {
                        window->Histo[g][ibin]+=1.;
                        window->N[g]++;
                    }
                    window->Ntot[g]++;
                }
            }
        }
        else if (t>opt->tmax)
        {
            if (opt->verbose)
                printf("time %f larger than tmax %f, stop reading this pullx/pullf file\n",t,opt->tmax);
            break;
        }
    }

    if (bGetMinMax)
    {
        *mintmp=minfound;
        *maxtmp=maxfound;
    }
}


void read_tpr_pullxf_files(char **fnTprs,char **fnPull,int nfiles,
        t_UmbrellaHeader* header,
        t_UmbrellaWindow **window, t_UmbrellaOptions *opt)
{
    int i;
    real mintmp,maxtmp;


    printf("Reading %d tpr and pullf files\n",nfiles/2);

    /* min and max not given? */
    if (opt->bAuto)
    {
        printf("Automatic determination of boundaries...\n");
        opt->min=1e20;
        opt->max=-1e20;
        for (i=0;i<nfiles; i++)
        {
            if (whaminFileType(fnTprs[i]) != whamin_tpr)
                gmx_fatal(FARGS,"Expected the %d'th file in input file to be a tpr file\n",i);
            read_tpr_header(fnTprs[i],header,opt);
            if (whaminFileType(fnPull[i]) != whamin_pullxf)
                gmx_fatal(FARGS,"Expected the %d'th file in input file to be a xvg (pullx/pullf) file\n",i);
            read_pull_xf(fnPull[i],fnTprs[i],header,NULL,opt,TRUE,&mintmp,&maxtmp);
            if (maxtmp>opt->max)
                opt->max=maxtmp;
            if (mintmp<opt->min)
                opt->min=mintmp;
        }
        printf("\nDetermined boundaries to %f and %f\n\n",opt->min,opt->max);
        if (opt->bBoundsOnly){
            printf("Found option -boundsonly, now exiting.\n");
            exit (0);
        }
    }
    /* store stepsize in profile */
    opt->dz=(opt->max-opt->min)/opt->bins;

    snew(*window,nfiles);
    for (i=0;i<nfiles; i++)
    {
        if (whaminFileType(fnTprs[i]) != whamin_tpr)
            gmx_fatal(FARGS,"Expected the %d'th file in input file to be a tpr file\n",i);
        read_tpr_header(fnTprs[i],header,opt);
        if (whaminFileType(fnPull[i]) != whamin_pullxf)
            gmx_fatal(FARGS,"Expected the %d'th file in input file to be a xvg (pullx/pullf) file\n",i);
        read_pull_xf(fnPull[i],fnTprs[i],header,*window+i,opt,FALSE,NULL,NULL);
    }
}


int gmx_wham(int argc,char *argv[])
{
    static char *desc[] = {
            "This is an analysis program that implements the Weighted",
            "Histogram Analysis Method (WHAM). It is intended to analyze",
            "output files generated by umbrella sampling simulations to ",
            "compute a potential of mean force (PMF). [PAR] ",
            "At present, three input modes are supported.[BR]",
            "[TT]* With option -it, the user provides a file which contains the[BR]",
            "[TT]  filenames of the umbrella simulation run-input files (tpr files),[BR]",
            "[TT]  AND, with option -ix, a file which contains filenames of [BR]",
            "[TT]  the pullx mdrun output files. The tpr and pullx files must [BR]",
            "[TT]  be in corresponding order, i.e. the first tpr created the [BR]",
            "[TT]  first pullx, etc.[BR]",
            "[TT]* Same as the previous input mode, except that the the user [BR]",
            "[TT]  provides the pull force ouput file names (pullf.xvg) with option -if.[BR]",
            "[TT]  From the pull force the position in the ubrella potential is [BR]",
            "[TT]  computed. This does not work with tabulated umbrella potentials.[BR]"
            "[TT]* With option -ip, the user provides filenames of (gzipped) pdo files, i.e.[BR]",
            "[TT]  the gromacs 3.3 umbrella output files. If you have some unusual[BR]"
            "[TT]  reaction coordinate you may also generate your own pdo files and [BR]",
            "[TT]  feed them with the -ip option into to g_wham. The pdo file header [BR]",
            "[TT]  must be similar to the folowing:[BR]",
            "# UMBRELLA      3.0[BR]",
            "# Component selection: 0 0 1[BR]",
            "# nSkip 1[BR]",
            "# Ref. Group 'TestAtom'[BR]",
            "# Nr. of pull groups 2[BR]",
            "# Group 1 'GR1'  Umb. Pos. 5.0 Umb. Cons. 1000.0[BR]",
            "# Group 2 'GR2'  Umb. Pos. 2.0 Umb. Cons. 500.0[BR]",
            "#####[BR]",
            "[TT]  Nr of pull groups, umbrella positions, force constants, and names[BR]",
            "[TT]  may (of course) differ. Following the header, a time column and[BR]",
            "[TT]  a data columns for each pull group follow (i.e. the displacement [BR]",
            "[TT]  with respect to the umbrella center). Up to four pull groups are possible [BR]",
            "[TT]  at present.[PAR]",
            "By default, the output files are[BR]",
            "  [TT]-o[tt]      PMF output file[BR]",
            "  [TT]-hist[tt]   histograms output file[PAR]",
            "The umbrella potential is assumed to be harmonic and the force constants are ",
            "read from the tpr or pdo files. If a non-harmonic umbrella force was applied ",
            "a tabulated potential can be provied with -tab.[PAR]",
            "WHAM OPTIONS[PAR]",
            "  [TT]-bins[tt]   Nr of bins used in analysis[BR]",
            "  [TT]-temp[tt]   Temperature in the simulations[BR]",
            "  [TT]-tol[tt]    Stop iteration if profile (probability) changed less than tolerance[BR]",
            "  [TT]-auto[tt]   Automatic determination of boudndaries[BR]",
            "  [TT]-min,-max[tt]   Boundaries of the profile [BR]",
            "The data points which are used ",
            "to compute the profile can be restricted with options -b, -e, and -dt. ",
            "Play particularly with -b to ensure sufficient equilibration in each ",
            "umbrella window![PAR]",
            "With -log (default) the profile is written in energy units, otherwise (-nolog) as ",
            "probability. The unit can be specified with -unit. With energy output, ",
            "the energy in the first bin is defined to be zero. If you want the free energy at a different ",
            "position to be zero, choose with -zprof0 (useful with bootstrapping, see below).[PAR]",
            "For cyclic (or periodic) reaction coordinates (dihedral angle, channel PMF",
            "without osmotic gradient), -cycl is useful.[BR]",
            " [TT]-cycl yes[tt]        min and max are assumed to "
            "be neighboring points and histogram points outside min and max are mapped into ",
            "the interval [min,max] (compare histogram output). [BR]",
            " [TT]-cycl weighted[tt]   First, a non-cyclic profile is computed. Subsequently, ",
            "periodicity is enforced by adding corrections dG(i) between neighboring bins",
            "i and i+1. The correction is chosen proportional to 1/[n(i)*n(i+1)]^alpha, where",
            "n(i) denotes the total nr of data points in bin i as collected from all histograms.",
            "alpha is defined with -alpha. The corrections are written to the file defined by -wcorr.",
            " (Compare Hub and de Groot, PNAS 105:1198 (2008))[PAR]",
            "ERROR ANALYSIS[BR]",
            "Statistical errors may be estimated with bootstrap analysis. Use it with care, ",
            "otherwise the statistical error may be substantially undererstimated !![BR]",
            "-nBootstrap defines the nr of bootstraps. Two bootstrapping modes are supported.[BR]",
            "  [TT]-histbs[tt]    Complete histograms are considered as independent data points (default). For each",
            "bootstrap, N histograms are randomly chosen from the N given histograms (allowing duplication).",
            "To avoid gaps without data along the reaction coordinate blocks of histograms (-histbs-block)",
            "may be defined. In that case, the given histograms are divided into blocks and ",
            "only histograms within each block are mixed. Note that the histograms",
            "within each block must be representative for all possible histograms, otherwise the",
            "statistical error is undererstimated![BR]",
            "  [TT]-nohistbs[tt]  The given histograms are used to generate new random histograms,",
            "such that the generated data points are distributed according the given histograms. The number",
            "of points generated for each bootstrap histogram can be controlled with -bs-dt.",
            "Note that one data point should be generated for each *independent* point in the given",
            "histograms. With the long autocorrelations in MD simulations, this procedure may ",
            "easily understimate the error![BR]",
            "Bootstrapping output:[BR]",
            "  [TT]-bsres[tt]   Average profile and standard deviations[BR]",
            "  [TT]-bsprof[tt]  All bootstrapping profiles[BR]",
            "With -vbs (verbose bootstrapping), the histograms of each bootstrap are written, and, ",
            "with -nohistBS, the cummulants of the histogram.",
    };

    static t_UmbrellaOptions opt;
    static bool bHistOnly=FALSE;

    static char *en_unit[]={NULL,"kJ","kCal","kT",NULL};
    static char *en_unit_label[]={"","E (kJ mol\\S-1\\N)","E (kcal mol\\S-1\\N)","E (kT)",};
    static char *en_cycl[]={NULL,"no","yes","weighted",NULL};

    t_pargs pa[] = {
            { "-min", FALSE, etREAL, {&opt.min},
              "Minimum coordinate in profile"},
            { "-max", FALSE, etREAL, {&opt.max},
              "Maximum coordinate in profile"},
            { "-auto", FALSE, etBOOL, {&opt.bAuto},
              "determine min and max automatically"},
            { "-bins",FALSE, etINT, {&opt.bins},
              "Number of bins in profile"},
            { "-temp", FALSE, etREAL, {&opt.Temperature},
              "Temperature"},
            { "-tol", FALSE, etREAL, {&opt.Tolerance},
              "Tolerance"},
            { "-v", FALSE, etBOOL, {&opt.verbose},
              "verbose mode"},
            { "-b", FALSE, etREAL, {&opt.tmin},
              "first time to analyse (ps)"},
            { "-e", FALSE, etREAL, {&opt.tmax},
              "last time to analyse (ps)"},
            { "-dt", FALSE, etREAL, {&opt.dt},
              "Analyse only every dt ps"},
            { "-histonly", FALSE, etBOOL, {&bHistOnly},
              "Write histograms and exit"},
            { "-boundsonly", FALSE, etBOOL, {&opt.bBoundsOnly},
              "Determine min and max and exit (with -auto)"},
            { "-log", FALSE, etBOOL, {&opt.bLog},
              "Calculate the log of the profile before printing"},
            { "-unit", FALSE,  etENUM, {en_unit},
              "energy unit in case of log output" },
            { "-zprof0", FALSE, etREAL, {&opt.zProf0},
              "Define profile to 0.0 at this position (with -log)"},
            { "-cycl", FALSE, etENUM, {en_cycl},
              "Create cyclic/periodic profile. Assumes min and max are the same point."},
            { "-alpha", FALSE, etREAL, {&opt.alpha},
              "for '-cycl weighted', set parameter alpha"},
            { "-flip", FALSE, etBOOL, {&opt.bFlipProf},
              "Combine halves of profile (not supported)"},
            { "-hist-eq", FALSE, etBOOL, {&opt.bHistEq},
              "Enforce equal weight for all histograms. (Non-Weighed-HAM)"},
            { "-nBootstrap", FALSE,  etINT, {&opt.nBootStrap},
              "nr of bootstraps to estimate statistical uncertainty" },
            { "-bs-dt", FALSE, etREAL, {&opt.dtBootStrap},
              "timestep for synthetic bootstrap histograms (ps). Ensure independent data points!"},
            { "-bs-seed", FALSE, etINT, {&opt.bsSeed},
              "seed for bootstrapping. (-1 = use time)"},
            { "-histbs", FALSE, etBOOL, {&opt.bHistBootStrap},
              "In bootstrapping, consider complete histograms as one data point."
              "Accounts better for long autocorrelations."},
            { "-histbs-block", FALSE, etINT, {&opt.histBootStrapBlockLength},
              "when mixin histograms only mix within blocks of -histBS_block."},
            { "-vbs", FALSE, etBOOL, {&opt.bs_verbose},                                                                                                                                                                              "verbose bootstrapping. Print the cummulants and a histogram file for each bootstrap."},
    };

    t_filenm fnm[] = {
            { efDAT, "-ix","pullx-files",ffOPTRD},     /* wham input: pullf.xvg's and tprs    */
            { efDAT, "-if","pullf-files",ffOPTRD},     /* wham input: pullf.xvg's and tprs    */
            { efDAT, "-it","tpr-files",ffOPTRD},       /* wham input: tprs                    */
            { efDAT, "-ip","pdo-files",ffOPTRD},       /* wham input: pdo files (gmx3 style)  */
            { efXVG, "-o", "profile", ffWRITE },       /* output file for profile */
            { efXVG, "-hist","histo", ffWRITE},	       /* output file for histograms */
            { efXVG, "-bsres","bsResult", ffOPTWR},    /* average and errors of bootstrap analysis */
            { efXVG, "-bsprof","bsProfs", ffOPTWR},    /* output file for bootstrap profiles       */
            { efDAT, "-tab","umb-pot",ffOPTRD},        /* Tabulated umbrella potential (if not harmonic) */
            { efXVG, "-wcorr","cycl-corr",ffOPTRD},    /* Corrections to profile in case -cycl weighted */
    };

    int i,j,l,nfiles,nwins,nfiles2;
    t_UmbrellaHeader header;
    t_UmbrellaWindow * window=NULL;
    double *profile,maxchange=1e20;
    bool bMinSet,bMaxSet,bAutoSet,bExact=FALSE;
    char **fninTpr,**fninPull,**fninPdo,*fnPull;
    FILE *histout,*profout;
    char ylabel[256],title[256];

    opt.bins=200;
    opt.verbose=FALSE;
    opt.cycl=enCycl_no;
    opt.tmin=50;
    opt.tmax=1e20;
    opt.dt=0.0;
    opt.bShift=TRUE;
    opt.nBootStrap=0;
    opt.dtBootStrap=0.0;
    opt.bsSeed=-1;
    opt.bHistBootStrap=TRUE;
    opt.histBootStrapBlockLength=12;
    opt.zProfZero=0.0;
    opt.bWeightedCycl=FALSE;
    opt.alpha=2;
    opt.bHistOutOnly=FALSE;
    opt.min=0;
    opt.max=0;
    opt.bLog=TRUE;
    opt.unit=en_kJ;
    opt.zProf0=0.0;
    opt.nBootStrap=0;
    opt.bsSeed=-1;
    opt.bHistBootStrap=TRUE;
    opt.histBootStrapBlockLength=8;
    opt.bs_verbose=FALSE;
    opt.Temperature=298;
    opt.bFlipProf=FALSE;
    opt.Tolerance=1e-6;
    opt.bAuto=TRUE;
    opt.bBoundsOnly=FALSE;


#define NFILE asize(fnm)

    CopyRight(stderr,argv[0]);
    parse_common_args(&argc,argv,PCA_BE_NICE,
            NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

    opt.unit=nenum(en_unit);
    opt.cycl=nenum(en_cycl);

    opt.bProf0Set=opt2parg_bSet("-zprof0",  asize(pa), pa);

    opt.bTab=opt2bSet("-tab",NFILE,fnm);
    opt.bPdo=opt2bSet("-ip",NFILE,fnm);
    opt.bTpr=opt2bSet("-it",NFILE,fnm);
    opt.bPullx=opt2bSet("-ix",NFILE,fnm);
    opt.bPullf=opt2bSet("-if",NFILE,fnm);
    if  (opt.bTab && opt.bPullf)
        gmx_fatal(FARGS,"Force input does not work with tabulated potentials. "
                "Provide pullx.xvg or pdo files!");

#define BOOLXOR(a,b) ( ((!(a))&&(b)) || ((a)&&(!(b))))
    if (!opt.bPdo && !BOOLXOR(opt.bPullx,opt.bPullf))
        gmx_fatal(FARGS,"Give either pullx (-ix) OR pullf (-if) data. Not both.");
    if ( !opt.bPdo && !(opt.bTpr || opt.bPullf || opt.bPullx))
        gmx_fatal(FARGS,"g_wham supports three input modes, pullx, pullf, or pdo file input."
                "\n\n Check g_wham -h !");

    opt.fnPdo=opt2fn("-ip",NFILE,fnm);
    opt.fnTpr=opt2fn("-it",NFILE,fnm);
    opt.fnPullf=opt2fn("-if",NFILE,fnm);
    opt.fnPullx=opt2fn("-ix",NFILE,fnm);

    bMinSet = opt2parg_bSet("-min",  asize(pa), pa);
    bMaxSet = opt2parg_bSet("-max",  asize(pa), pa);
    bAutoSet = opt2parg_bSet("-auto",  asize(pa), pa);
    if ( (bMinSet || bMaxSet) && bAutoSet)
        gmx_fatal(FARGS,"With -auto, do not give -min or -max\n");

    if ( (bMinSet && !bMaxSet) || (!bMinSet && bMaxSet))
        gmx_fatal(FARGS,"When giving -min, you must give -max (and vice versa), too\n");

    if (bMinSet && opt.bAuto)
    {
        printf("Note: min and max given, switching off -auto.\n");
        opt.bAuto=FALSE;
    }

    /* Reading gmx4 pull output and tpr files */
    if (opt.bTpr || opt.bPullf || opt.bPullx)
    {
        read_wham_in(opt.fnTpr,&fninTpr,&nfiles,&opt);

        fnPull=opt.bPullf ? opt.fnPullf : opt.fnPullx;
        read_wham_in(fnPull,&fninPull,&nfiles2,&opt);
        printf("Found %d tpr and %d pull %s files in %s and %s, respectively\n",
                nfiles,nfiles2,opt.bPullf ? "force" : "position",opt.fnTpr,fnPull);
        if (nfiles!=nfiles2)
            gmx_fatal(FARGS,"Found %d file names in %s, but %d in %s\n",nfiles,
                    opt.fnTpr,nfiles2,fnPull);
        read_tpr_pullxf_files(fninTpr,fninPull,nfiles, &header, &window, &opt);
    }
    else
    {
        /* reading pdo files */
        read_wham_in(opt.fnPdo,&fninPdo,&nfiles,&opt);
        printf("Found %d pdo files in %s\n",nfiles,opt.fnPdo);
        read_pdo_files(fninPdo,nfiles, &header, &window, &opt);
    }
    nwins=nfiles;

    /* enforce equal weight for all histograms? */
    if (opt.bHistEq)
        enforceEqualWeights(window,nwins);

    /* write histograms */
    histout=xvgropen(opt2fn("-hist",NFILE,fnm),"Umbrella histograms",
            "z","count");
    for(l=0;l<opt.bins;++l)
    {
        fprintf(histout,"%e\t",(double)(l+0.5)/opt.bins*(opt.max-opt.min)+opt.min);
        for(i=0;i<nwins;++i)
        {
            for(j=0;j<window[i].nPull;++j)
            {
                fprintf(histout,"%e\t",window[i].Histo[j][l]);
            }
        }
        fprintf(histout,"\n");
    }
    ffclose(histout);
    if (bHistOnly)
    {
        printf("Wrote histograms to %s, now exiting.\n",opt2fn("-hist",NFILE,fnm));
        return 0;
    }

    /* Using tabulated umbrella potential */
    if (opt.bTab)
        setup_tab(opt2fn("-tab",NFILE,fnm),&opt);

    setup_acc_wham(window,nwins,&opt);

    /* Calculate profile */
    snew(profile,opt.bins);
    opt.stepchange=(opt.verbose)? 1 : 100;
    i=0;
    do {
        if (maxchange<opt.Tolerance)
        {
            bExact=TRUE;
            /* if (opt.verbose) */
            printf("Switched to exact iteration in iteration %d\n",i);
        }
        calc_profile(profile,window,nwins,&opt,bExact);
        if (((i%opt.stepchange) == 0 || i==1) && !i==0)
            printf("\t%4d) Maximum change %e\n",i,maxchange);
        i++;
    } while( (maxchange=calc_z(profile, window, nwins, &opt,bExact)) > opt.Tolerance || !bExact);
    printf("Converged in %d iterations. Final maximum change %g\n",i,maxchange);

    /* Write profile in energy units? */
    if (opt.bLog)
    {
        prof_normalization_and_unit(profile,&opt);
        strcpy(ylabel,en_unit_label[opt.unit]);
        strcpy(title,"Umbrella potential");
    }
    else
    {
        strcpy(ylabel,"Density of states");
        strcpy(title,"Density of states");
    }
    /* Force cyclic profile by wheighted correction? */
    if (opt.cycl==enCycl_weighted)
        cyclicProfByWeightedCorr(profile,window,nwins,&opt,TRUE,
                opt2fn("-wcorr",NFILE,fnm));

    profout=xvgropen(opt2fn("-o",NFILE,fnm),title,"z",ylabel);
    for(i=0;i<opt.bins;++i)
        fprintf(profout,"%e\t%e\n",(double)(i+0.5)/opt.bins*(opt.max-opt.min)+opt.min,profile[i]);
    ffclose(profout);
    printf("Wrote %s\n",opt2fn("-o",NFILE,fnm));

    /* Bootstrap Method */
    if (opt.nBootStrap)
        do_bootstrapping(opt2fn("-bsres",NFILE,fnm),opt2fn("-bsprof",NFILE,fnm),
                opt2fn("-hist",NFILE,fnm),
                ylabel, profile, window, nwins, &opt);

    thanx(stderr);
    return 0;
}
