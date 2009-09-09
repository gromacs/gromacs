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

#include <string.h>
#include <math.h>
#include <ctype.h>

#include "typedefs.h"
#include "gmx_fatal.h"
#include "vec.h"
#include "string2.h"
#include "smalloc.h"
#include "enxio.h"
#include "statutil.h"
#include "names.h"
#include "copyrite.h"
#include "macros.h"
#include "xvgr.h"
#include "gstat.h"
#include "physics.h"
#include "tpxio.h"
#include "viewit.h"
#include "mtop_util.h"
#include "gmx_statistics.h"

static real       minthird=-1.0/3.0,minsixth=-1.0/6.0;

typedef struct {
  int    nframes;
  double *sum;
  double *sum2;
} enersum_t;

static double mypow(double x,double y)
{
  if (x > 0)
    return pow(x,y);
  else 
    return 0.0;
}

static int *select_it(int nre,char *nm[],int *nset)
{
  bool *bE;
  int  n,k,j,i;
  int  *set;
  bool bVerbose = TRUE;
  
  if ((getenv("VERBOSE")) != NULL)
    bVerbose = FALSE;
  
  fprintf(stderr,"Select the terms you want from the following list\n");
  fprintf(stderr,"End your selection with 0\n");

  if ( bVerbose ) {
    for(k=0; (k<nre); ) {
      for(j=0; (j<4) && (k<nre); j++,k++) 
	fprintf(stderr," %3d=%14s",k+1,nm[k]);
      fprintf(stderr,"\n");
    }
  }

  snew(bE,nre);
  do {
    if(1 != scanf("%d",&n))
    {
      gmx_fatal(FARGS,"Error reading user input");
    }
    if ((n>0) && (n<=nre))
      bE[n-1]=TRUE;
  } while (n != 0);

  snew(set,nre);
  for(i=(*nset)=0; (i<nre); i++)
    if (bE[i])
      set[(*nset)++]=i;
 
  sfree(bE);
  
  return set;
}

static int strcount(const char *s1,const char *s2)
{
  int n=0;
  while (s1 && s2 && (toupper(s1[n]) == toupper(s2[n])))
    n++;
  return n;
}

static void chomp(char *buf)
{
  int len = strlen(buf);
  
  while ((len > 0) && (buf[len-1] == '\n')) {
    buf[len-1] = '\0';
    len--;
  }
}

static int *select_by_name(int nre,gmx_enxnm_t *nm,int *nset)
{
  bool *bE;
  int  n,k,kk,j,i,nmatch,nind,nss;
  int  *set;
  bool bEOF,bVerbose = TRUE,bLong=FALSE;
  char *ptr,buf[STRLEN];
  const char *fm4="%3d  %-14s";
  const char *fm2="%3d  %-34s";
  char **newnm=NULL;
  
  if ((getenv("VERBOSE")) != NULL)
    bVerbose = FALSE;
  
  fprintf(stderr,"\n");
  fprintf(stderr,"Select the terms you want from the following list by\n");
  fprintf(stderr,"selecting either (part of) the name or the number or a combination.\n");
  fprintf(stderr,"End your selection with an empty line or a zero.\n");
  fprintf(stderr,"-------------------------------------------------------------------\n");
  
  snew(newnm,nre);
  j = 0;
  for(k=0; k<nre; k++) {
    newnm[k] = strdup(nm[k].name);
    /* Insert dashes in all the names */
    while ((ptr = strchr(newnm[k],' ')) != NULL) {
      *ptr = '-';
    }
    if ( bVerbose ) {
      if (j == 0) {
	if (k > 0) {
	  fprintf(stderr,"\n");
	}
	bLong = FALSE;
	for(kk=k; kk<k+4; kk++) {
	  if (kk < nre && strlen(nm[kk].name) > 14) {
	    bLong = TRUE;
	  }
	}
      } else {
	fprintf(stderr," ");
      }
      if (!bLong) {
	fprintf(stderr,fm4,k+1,newnm[k]);
	j++;
	if (j == 4) {
	  j = 0;
	}
      } else {
	fprintf(stderr,fm2,k+1,newnm[k]);
	j++;
	if (j == 2) {
	  j = 0;
	}
      }
    }
  }
  if ( bVerbose ) {
    fprintf(stderr,"\n\n");
  }
  
  snew(bE,nre);
  
  bEOF = FALSE;
  while (!bEOF && (fgets2(buf,STRLEN-1,stdin))) {
    /* Remove newlines */
    chomp(buf);
    
    /* Remove spaces */
    trim(buf);
    
    /* Empty line means end of input */
    bEOF = (strlen(buf) == 0);
    if (!bEOF) {
      ptr = buf;
      do {
	if (!bEOF) {
	  /* First try to read an integer */
	  nss   = sscanf(ptr,"%d",&nind);
	  if (nss == 1) {
	    /* Zero means end of input */
	    if (nind == 0) {
	      bEOF = TRUE;
	    } else if ((1<=nind) && (nind<=nre)) {
	      bE[nind-1] = TRUE;
	    } else {
	      fprintf(stderr,"number %d is out of range\n",nind);
	    }
	  }
	  else {
	    /* Now try to read a string */
	    i = strlen(ptr);
	    nmatch = 0;
	    for(nind=0; nind<nre; nind++) {
	      if (strcasecmp(newnm[nind],ptr) == 0) {
		bE[nind] = TRUE;
		nmatch++;
	      }
	    }
	    if (nmatch == 0) {
	      i = strlen(ptr);
	      nmatch = 0;
	      for(nind=0; nind<nre; nind++) {
		if (strncasecmp(newnm[nind],ptr,i) == 0) {
		  bE[nind] = TRUE;
		  nmatch++;
		}
	      }
	      if (nmatch == 0) {
		fprintf(stderr,"String '%s' does not match anything\n",ptr);
	      }
	    }
	  }
	}
	/* Look for the first space, and remove spaces from there */
	if ((ptr = strchr(ptr,' ')) != NULL) {
	  trim(ptr);
	}
      } while (!bEOF && (ptr && (strlen(ptr) > 0)));
    }
  }
  
  snew(set,nre);
  for(i=(*nset)=0; (i<nre); i++)
    if (bE[i])
      set[(*nset)++]=i;
 
  sfree(bE);
  
  if (*nset == 0)
    gmx_fatal(FARGS,"No energy terms selected");

  for(i=0; (i<nre); i++) 
    sfree(newnm[i]);
  sfree(newnm);
  
  return set;
}

static void get_orires_parms(const char *topnm,
			     int *nor,int *nex,int **label,real **obs)
{
  gmx_mtop_t mtop;
  gmx_localtop_t *top;
  t_inputrec ir;
  t_iparams  *ip;
  int        natoms,i;
  t_iatom    *iatom;
  int        nb;
  matrix     box;

  read_tpx(topnm,&ir,box,&natoms,NULL,NULL,NULL,&mtop);
  top = gmx_mtop_generate_local_top(&mtop,&ir);

  ip       = top->idef.iparams;
  iatom    = top->idef.il[F_ORIRES].iatoms;
  
  /* Count how many distance restraint there are... */
  nb = top->idef.il[F_ORIRES].nr;
  if (nb == 0)
    gmx_fatal(FARGS,"No orientation restraints in topology!\n");
  
  *nor = nb/3;
  *nex = 0;
  snew(*label,*nor);
  snew(*obs,*nor);
  for(i=0; i<nb; i+=3) {
    (*label)[i/3] = ip[iatom[i]].orires.label;
    (*obs)[i/3]   = ip[iatom[i]].orires.obs;
    if (ip[iatom[i]].orires.ex >= *nex)
      *nex = ip[iatom[i]].orires.ex+1;
  }
  fprintf(stderr,"Found %d orientation restraints with %d experiments",
	  *nor,*nex);
}

static int get_bounds(const char *topnm,
		      real **bounds,int **index,int **dr_pair,int *npairs,
		      gmx_mtop_t *mtop,gmx_localtop_t **ltop,t_inputrec *ir)
{
  gmx_localtop_t *top;
  t_functype *functype;
  t_iparams  *ip;
  int        natoms,i,j,k,type,ftype,natom;
  t_ilist    *disres;
  t_iatom    *iatom;
  real       *b;
  int        *ind,*pair;
  int        nb,label1;
  matrix     box;

  read_tpx(topnm,ir,box,&natoms,NULL,NULL,NULL,mtop);
  snew(*ltop,1);
  top = gmx_mtop_generate_local_top(mtop,ir);
  *ltop = top;

  functype = top->idef.functype;
  ip       = top->idef.iparams;
  
  /* Count how many distance restraint there are... */
  nb=top->idef.il[F_DISRES].nr;
  if (nb == 0)
    gmx_fatal(FARGS,"No distance restraints in topology!\n");
  
  /* Allocate memory */
  snew(b,nb);
  snew(ind,nb);
  snew(pair,nb+1);
  
  /* Fill the bound array */
  nb=0;
  for(i=0; (i<top->idef.ntypes); i++) {
    ftype = functype[i];
    if (ftype == F_DISRES) {

      label1 = ip[i].disres.label;
      b[nb]   = ip[i].disres.up1;
      ind[nb] = label1;
      nb++;
    }
  }
  *bounds = b;
  
  /* Fill the index array */
  label1  = -1;
  disres  = &(top->idef.il[F_DISRES]);
  iatom   = disres->iatoms;
  for(i=j=k=0; (i<disres->nr); ) {
    type  = iatom[i];
    ftype = top->idef.functype[type];
    natom = interaction_function[ftype].nratoms+1;
    if (label1 != top->idef.iparams[type].disres.label) {
      pair[j] = k;
      label1  = top->idef.iparams[type].disres.label; 
      j ++;
    }
    k++;
    i += natom;
  }
  pair[j]  = k;
  *npairs = k;
  if (j != nb)
    gmx_incons("get_bounds for distance restraints");

  *index   = ind;
  *dr_pair = pair;
  
  return nb;
}

static void calc_violations(real rt[],real rav3[],int nb,int index[],
			    real bounds[],real *viol,double *st,double *sa)
{
  const   real sixth=1.0/6.0;
  int     i,j;
  double  rsum,rav,sumaver,sumt;
  
  sumaver = 0;
  sumt    = 0;
  for(i=0; (i<nb); i++) {
    rsum = 0.0;
    rav  = 0.0;
    for(j=index[i]; (j<index[i+1]); j++) {
      if (viol)
	viol[j] += mypow(rt[j],-3.0);
      rav     += sqr(rav3[j]);
      rsum    += mypow(rt[j],-6);
    }
    rsum    = max(0.0,mypow(rsum,-sixth)-bounds[i]);
    rav     = max(0.0,mypow(rav, -sixth)-bounds[i]);
    
    sumt    += rsum;
    sumaver += rav;
  }
  *st = sumt;
  *sa = sumaver;
}

static void analyse_disre(const char *voutfn,    int nframes,
			  real violaver[], real bounds[], int index[],
			  int pair[],      int nbounds,
                          const output_env_t oenv)
{
    FILE   *vout;
    double sum,sumt,sumaver;
    int    i,j;

    /* Subtract bounds from distances, to calculate violations */
    calc_violations(violaver,violaver,
            nbounds,pair,bounds,NULL,&sumt,&sumaver);

#ifdef DEBUG
    fprintf(stdout,"\nSum of violations averaged over simulation: %g nm\n",
            sumaver);
    fprintf(stdout,"Largest violation averaged over simulation: %g nm\n\n",
            sumt);
#endif		    
    vout=xvgropen(voutfn,"r\\S-3\\N average violations","DR Index","nm",
            oenv);
    sum  = 0.0;
    sumt = 0.0;
    for(i=0; (i<nbounds); i++) {
        /* Do ensemble averaging */
        sumaver = 0;
        for(j=pair[i]; (j<pair[i+1]); j++) 
            sumaver += sqr(violaver[j]/nframes); 
        sumaver = max(0.0,mypow(sumaver,minsixth)-bounds[i]);

        sumt   += sumaver;
        sum     = max(sum,sumaver);
        fprintf(vout,"%10d  %10.5e\n",index[i],sumaver);
    }
#ifdef DEBUG
    for(j=0; (j<dr.ndr); j++)
        fprintf(vout,"%10d  %10.5e\n",j,mypow(violaver[j]/nframes,minthird));
#endif
    ffclose(vout);

    fprintf(stdout,"\nSum of violations averaged over simulation: %g nm\n",
            sumt);
    fprintf(stdout,"Largest violation averaged over simulation: %g nm\n\n",sum);

    do_view(oenv,voutfn,"-graphtype bar");
}

static void einstein_visco(const char *fn,const char *fni,int nsets,
                           int nframes,real **sum,
                           real V,real T,int nsteps,double time[],
                           const output_env_t oenv)
{
    FILE *fp0,*fp1;
    double av[4],avold[4];
    double fac,dt,di;
    int  i,j,m,nf4;

    if (nframes < 1)
        return;

    dt  = (time[1]-time[0]);
    nf4 = nframes/4+1;

    for(i=0; i<=nsets; i++)
        avold[i] = 0;
    fp0=xvgropen(fni,"Shear viscosity integral",
            "Time (ps)","(kg m\\S-1\\N s\\S-1\\N ps)",oenv);
    fp1=xvgropen(fn,"Shear viscosity using Einstein relation",
            "Time (ps)","(kg m\\S-1\\N s\\S-1\\N)",oenv);
    for(i=1; i<nf4; i++) {
        fac = dt*nframes/nsteps;
        for(m=0; m<=nsets; m++)
            av[m] = 0;
        for(j=0; j<nframes-i; j++) {
            for(m=0; m<nsets; m++) {
                di   = sqr(fac*(sum[m][j+i]-sum[m][j]));

                av[m]     += di;
                av[nsets] += di/nsets;
            }
        }
        /* Convert to SI for the viscosity */
        fac = (V*NANO*NANO*NANO*PICO*1e10)/(2*BOLTZMANN*T)/(nframes-i);
        fprintf(fp0,"%10g",time[i]-time[0]);
        for(m=0; (m<=nsets); m++) {
            av[m] = fac*av[m];
            fprintf(fp0,"  %10g",av[m]);
        }
        fprintf(fp0,"\n");
        fprintf(fp1,"%10g",0.5*(time[i]+time[i-1])-time[0]);
        for(m=0; (m<=nsets); m++) {
            fprintf(fp1,"  %10g",(av[m]-avold[m])/dt);
            avold[m] = av[m];
        }
        fprintf(fp1,"\n");
    }
    fclose(fp0);
    fclose(fp1);
}

static void analyse_ener(bool bCorr,const char *corrfn,
			 bool bFee,bool bSum,bool bFluct,
			 bool bVisco,const char *visfn,int  nmol,int ndf,
			 gmx_large_int_t start_step,double start_t,
			 gmx_large_int_t step,double t,
			 double time[], real reftemp,
			 gmx_large_int_t ee_nsum,
			 t_energy ee_sum[],enersum_t *enersum,
			 int nset,int set[],int nenergy,real **eneset,
			 real **enesum,
			 char **leg,gmx_enxnm_t *enm,
			 real Vaver,real ezero, const output_env_t oenv)
{
  FILE *fp;
  /* Check out the printed manual for equations! */
  double Dt,aver,avertot,stddev,delta_t,totaldrift;
  real a,b,r;
  real xxx,integral,intBulk;
  real sfrac,oldfrac,diffsum,diffav,fstep,pr_aver,pr_stddev,fluct2;
  double beta=0,expE,expEtot,*fee=NULL;
  gmx_large_int_t nsteps;
  int  nexact,nnotexact,iset;
  double x1m,x1mk;
  real Temp=-1,Pres=-1,VarV=-1,VarT=-1;
  int  i,j;
  real chi2;
  bool bIsEner;
  char buf[256];

  nsteps  = step - start_step + 1;
  if (nsteps < 1) {
    fprintf(stdout,"Not enough steps (%s) for statistics\n",
	    gmx_step_str(nsteps,buf));
  }
  else {
    /* Calculate the time difference */
    delta_t = t - start_t;
    
    fprintf(stdout,"\nStatistics over %s steps [ %.4f thru %.4f ps ], %d data sets\n",
	    gmx_step_str(nsteps,buf),start_t,t,nset);

    if (ee_nsum == 0) {
      nexact    = 0;
      nnotexact = nset;
    } else {
      nexact    = 0;
      nnotexact = 0;
      for(i=0; (i<nset); i++) {
	if (ee_sum[set[i]].esum != 0) {
	  nexact++;
	} else {
	  nnotexact++;
	}
      }
    }
    
    if (nnotexact == 0) {
      fprintf(stdout,"All averages are over %s frames\n",
	      gmx_step_str(ee_nsum,buf));
    } else if (nexact == 0 || ee_nsum == enersum->nframes) {
      fprintf(stdout,"All averages are over %d frames\n",enersum->nframes);
    } else {
      fprintf(stdout,"The term%s",nnotexact==1 ? "" : "s");
      for(i=0; (i<nset); i++) {
	if (ee_sum[set[i]].esum == 0) {
	  fprintf(stdout," '%s'",leg[i]);
	}
      }
      fprintf(stdout," %s averaged over %d frames\n",
	      nnotexact==1 ? "is" : "are",enersum->nframes);
      fprintf(stdout,"All other averages are over %s frames\n",
	      gmx_step_str(ee_nsum,buf));
    }
    fprintf(stdout,"\n");

    fprintf(stdout,"%-24s %10s %10s %10s %10s",
	    "Energy","Average","RMSD","Fluct.","Tot-Drift");
    if (bFee)
      fprintf(stdout,"  %10s\n","-kT ln<e^(E/kT)>");
    else
      fprintf(stdout,"\n");
    fprintf(stdout,"-------------------------------------------------------------------------------\n");
    
    /* Initiate locals, only used with -sum */
    avertot=0;
    expEtot=0;
    if (bFee) {
      beta = 1.0/(BOLTZ*reftemp);
      snew(fee,nset);
    }
    for(i=0; (i<nset); i++) {
      iset = set[i];
      if (ee_nsum > 0 && ee_sum[iset].esum != 0) {
	/* We have exact sums over all the MD steps */
	aver   = ee_sum[iset].esum/ee_nsum;
	stddev = sqrt(ee_sum[iset].eav/ee_nsum);
      } else {
	/* We only have data at energy file frame steps */
	aver   = enersum->sum[i]/enersum->nframes;
	if (enersum->nframes > 1) {
	  stddev = sqrt(enersum->sum2[i]/enersum->nframes - aver*aver);
	} else {
	  stddev = 0;
	}
      }

      if (bSum) 
	avertot+=aver;
      if (bFee) {
	expE = 0;
	for(j=0; (j<nenergy); j++) {
	  expE += exp(beta*(eneset[i][j]-aver)/nmol);
	}
	if (bSum) 
	  expEtot+=expE/nenergy;
	
	fee[i] = log(expE/nenergy)/beta + aver/nmol;
      }
      if (strstr(leg[i],"empera") != NULL) {
	VarT = sqr(stddev);
	Temp = aver;
      } else if (strstr(leg[i],"olum") != NULL) {
	VarV = sqr(stddev);
	Vaver= aver;
      } else if (strstr(leg[i],"essure") != NULL) {
	Pres = aver;
      }
      bIsEner = FALSE;
      for (j=0; (j <= F_ETOT); j++)
	bIsEner = bIsEner || 
	  (strcasecmp(interaction_function[j].longname,leg[i]) == 0);
      if (bIsEner) {
	pr_aver   = aver/nmol-ezero;
	pr_stddev = stddev/nmol;
      }
      else {
	pr_aver   = aver;
	pr_stddev = stddev;
      }
      if (nenergy > 1) {
	lsq_y_ax_b_xdouble(nenergy,time,eneset[i],&a,&b,&r,&chi2);
	totaldrift = a * delta_t;
      } else {
	totaldrift = 0;
      }
      fluct2 = sqr(pr_stddev) - sqr(totaldrift)/12;
      if (fluct2 < 0)
	fluct2 = 0;
      fprintf(stdout,"%-24s %10g %10g %10g %10g",
	      leg[i],pr_aver,pr_stddev,sqrt(fluct2),totaldrift);
      if (bFee) 
	fprintf(stdout,"  %10g",fee[i]);
      
      fprintf(stdout,"  (%s)\n",enm[set[i]].unit);

      if (bFluct) {
	for(j=0; (j<nenergy); j++)
	  eneset[i][j] -= aver;
      }
    }
    if (bSum) {
      fprintf(stdout,"%-24s %10g %10s %10s %10s %10s",
	      "Total",avertot/nmol,"--","--","--","--");
      /* pr_aver,pr_stddev,a,totaldrift */
      if (bFee) 
	fprintf(stdout,"  %10g  %10g\n",
		log(expEtot)/beta + avertot/nmol,log(expEtot)/beta);
      else
	fprintf(stdout,"\n");
    }
    if (Temp != -1) {
      real factor;
      
      factor = nmol*ndf*VarT/(3.0*sqr(Temp));
      fprintf(stdout,"Heat Capacity Cv:   %10g J/mol K (factor = %g)\n",
	      1000*BOLTZ/(2.0/3.0 - factor),factor);
    }
    if ((VarV != -1) && (Temp != -1)) {
      real tmp = VarV/(Vaver*BOLTZ*Temp*PRESFAC);
      
      fprintf(stdout,"Isothermal Compressibility: %10g /%s\n",
	      tmp,unit_pres_bar);
      fprintf(stdout,"Adiabatic bulk modulus:     %10g  %s\n",
	      1.0/tmp,unit_pres_bar);
    }
    /* Do correlation function */
    Dt = delta_t/nenergy;
    if (bVisco) {
      char *leg[] = { "Shear", "Bulk" };
      real factor;
    
      /* Assume pressure tensor is in Pxx Pxy Pxz Pyx Pyy Pyz Pzx Pzy Pzz */
      
      /* Symmetrise tensor! (and store in first three elements) 
       * And subtract average pressure!
       */
      for(i=0; (i<nenergy); i++) {
	eneset[0][i] = 0.5*(eneset[1][i]+eneset[3][i]);
	eneset[1][i] = 0.5*(eneset[2][i]+eneset[6][i]);
	eneset[2][i] = 0.5*(eneset[5][i]+eneset[7][i]);
	eneset[11][i] -= Pres;
	enesum[0][i] = 0.5*(enesum[1][i]+enesum[3][i]);
	enesum[1][i] = 0.5*(enesum[2][i]+enesum[6][i]);
	enesum[2][i] = 0.5*(enesum[5][i]+enesum[7][i]);
      }
      
      einstein_visco("evisco.xvg","eviscoi.xvg",
		     3,nenergy,enesum,Vaver,Temp,nsteps,time,oenv);
      
      /*do_autocorr(corrfn,buf,nenergy,3,eneset,Dt,eacNormal,TRUE);*/
      /* Do it for shear viscosity */
      strcpy(buf,"Shear Viscosity");
      low_do_autocorr(corrfn,oenv,buf,nenergy,3,(nenergy+1)/2,eneset,Dt,
		      eacNormal,1,TRUE,FALSE,FALSE,0.0,0.0,0,1);
	
      /* Now for bulk viscosity */
      strcpy(buf,"Bulk Viscosity");
      low_do_autocorr(corrfn,oenv,buf,nenergy,1,(nenergy+1)/2,&(eneset[11]),Dt,
		      eacNormal,1,TRUE,FALSE,FALSE,0.0,0.0,0,1);
      
      factor = (Vaver*1e-26/(BOLTZMANN*Temp))*Dt;
      fp=xvgropen(visfn,buf,"Time (ps)","\\8h\\4 (cp)",oenv);
      xvgr_legend(fp,asize(leg),leg,oenv);
      
      /* Use trapezium rule for integration */
      integral = 0;
      intBulk  = 0;
      for(i=1; (i<nenergy/2); i++) {
	integral += 0.5*(eneset[0][i-1]  + eneset[0][i])*factor;
	intBulk  += 0.5*(eneset[11][i-1] + eneset[11][i])*factor;
	fprintf(fp,"%10g  %10g  %10g\n",(i*Dt),integral,intBulk);
      }
      fclose(fp);
    }
    else if (bCorr) {
      if (bFluct)
	strcpy(buf,"Autocorrelation of Energy Fluctuations");
      else
	strcpy(buf,"Energy Autocorrelation");
      do_autocorr(corrfn,oenv,buf,nenergy,
		  bSum ? 1                 : nset,
		  bSum ? &(eneset[nset-1]) : eneset,
		  (delta_t/nenergy),eacNormal,FALSE);
    }
  }
}

static void print_time(FILE *fp,double t)
{
  fprintf(fp,"%12.6f",t);
}

static void print1(FILE *fp,bool bDp,real e)
{
  if (bDp)
    fprintf(fp,"  %16.12f",e);
  else
    fprintf(fp,"  %10.6f",e);

}

static void fec(const char *ene2fn, const char *runavgfn, 
		real reftemp, int nset, int set[], char *leg[], 
		int nenergy, real **eneset, double time[],
                const output_env_t oenv)
{
  char *ravgleg[] = { "\\8D\\4E = E\\sB\\N-E\\sA\\N", 
	   	      "<e\\S-\\8D\\4E/kT\\N>\\s0..t\\N" };
  FILE *fp;
  ener_file_t enx;
  int  nre,timecheck,step,nenergy2,maxenergy;
  int  i,j;
  bool bCont;
  real aver, beta;
  real **eneset2;
  double dE, sum;
  gmx_enxnm_t *enm=NULL;
  t_enxframe *fr;
  char buf[22];
  
  /* read second energy file */
  snew(fr,1);
  enm = NULL;
  enx = open_enx(ene2fn,"r");
  do_enxnms(enx,&(fr->nre),&enm);
  
  snew(eneset2,nset+1);
  nenergy2=0;
  maxenergy=0;
  timecheck=0;
  do {
    /* This loop searches for the first frame (when -b option is given), 
     * or when this has been found it reads just one energy frame
     */
    do {
      bCont = do_enx(enx,fr);
      
      if (bCont)
	timecheck = check_times(fr->t);
      
    } while (bCont && (timecheck < 0));
    
    /* Store energies for analysis afterwards... */
    if ((timecheck == 0) && bCont) {
      if (fr->nre > 0) {
	if ( nenergy2 >= maxenergy ) {
	  maxenergy += 1000;
	  for(i=0; i<=nset; i++)
	    srenew(eneset2[i],maxenergy);
	}
	if (fr->t != time[nenergy2])
	  fprintf(stderr,"\nWARNING time mismatch %g!=%g at frame %s\n",
		  fr->t, time[nenergy2], gmx_step_str(fr->step,buf));
	for(i=0; i<nset; i++)
	  eneset2[i][nenergy2] = fr->ener[set[i]].e;
	nenergy2++;
      }
    }
  } while (bCont && (timecheck == 0));
  
  /* check */
  if(nenergy!=nenergy2)
    fprintf(stderr,"\nWARNING file length mismatch %d!=%d\n",nenergy,nenergy2);
  nenergy=min(nenergy,nenergy2);
  
  /* calculate fe difference dF = -kT ln < exp(-(E_B-E_A)/kT) >_A */
  fp=NULL;
  if (runavgfn) {
    fp=xvgropen(runavgfn,"Running average free energy difference",
		"Time (" unit_time ")","\\8D\\4E (" unit_energy ")",oenv);
    xvgr_legend(fp,asize(ravgleg),ravgleg,oenv);
  }
  fprintf(stdout,"\n%-24s %10s\n",
	  "Energy","dF = -kT ln < exp(-(EB-EA)/kT) >A");
  sum=0;
  beta = 1.0/(BOLTZ*reftemp);
  for(i=0; i<nset; i++) {
    if (strcasecmp(leg[i],enm[set[i]].name)!=0)
      fprintf(stderr,"\nWARNING energy set name mismatch %s!=%s\n",
	      leg[i],enm[set[i]].name);
    for(j=0; j<nenergy; j++) {
      dE = eneset2[i][j]-eneset[i][j];
      sum += exp(-dE*beta);
      if (fp)
	fprintf(fp,"%10g %10g %10g\n", 
		time[j], dE, -BOLTZ*reftemp*log(sum/(j+1)) );
    }
    aver = -BOLTZ*reftemp*log(sum/nenergy);
    fprintf(stdout,"%-24s %10g\n",leg[i],aver);
  }
  if(fp) ffclose(fp);
  sfree(fr);
}

int gmx_energy(int argc,char *argv[])
{
  const char *desc[] = {
    
    "g_energy extracts energy components or distance restraint",
    "data from an energy file. The user is prompted to interactively",
    "select the energy terms she wants.[PAR]",
    
    "Average and RMSD are calculated with full precision from the",
    "simulation (see printed manual). Drift is calculated by performing",
    "a LSQ fit of the data to a straight line. The reported total drift",
    "is the difference of the fit at the first and last point.",
    "The term fluctuation gives the RMSD around the LSQ fit.[PAR]",
    
    "When the [TT]-viol[tt] option is set, the time averaged",
    "violations are plotted and the running time-averaged and",
    "instantaneous sum of violations are recalculated. Additionally",
    "running time-averaged and instantaneous distances between",
    "selected pairs can be plotted with the [TT]-pairs[tt] option.[PAR]",

    "Options [TT]-ora[tt], [TT]-ort[tt], [TT]-oda[tt], [TT]-odr[tt] and",
    "[TT]-odt[tt] are used for analyzing orientation restraint data.",
    "The first two options plot the orientation, the last three the",
    "deviations of the orientations from the experimental values.",
    "The options that end on an 'a' plot the average over time",
    "as a function of restraint. The options that end on a 't'",
    "prompt the user for restraint label numbers and plot the data",
    "as a function of time. Option [TT]-odr[tt] plots the RMS",
    "deviation as a function of restraint.",
    "When the run used time or ensemble averaged orientation restraints,",
    "option [TT]-orinst[tt] can be used to analyse the instantaneous,",
    "not ensemble-averaged orientations and deviations instead of",
    "the time and ensemble averages.[PAR]",

    "Option [TT]-oten[tt] plots the eigenvalues of the molecular order",
    "tensor for each orientation restraint experiment. With option",
    "[TT]-ovec[tt] also the eigenvectors are plotted.[PAR]",

    "With [TT]-fee[tt] an estimate is calculated for the free-energy",
    "difference with an ideal gas state: [BR]",
    "  Delta A = A(N,V,T) - A_idgas(N,V,T) = kT ln < e^(Upot/kT) >[BR]",
    "  Delta G = G(N,p,T) - G_idgas(N,p,T) = kT ln < e^(Upot/kT) >[BR]",
    "where k is Boltzmann's constant, T is set by [TT]-fetemp[tt] and"
    "the average is over the ensemble (or time in a trajectory).",
    "Note that this is in principle",
    "only correct when averaging over the whole (Boltzmann) ensemble",
    "and using the potential energy. This also allows for an entropy",
    "estimate using:[BR]",
    "  Delta S(N,V,T) = S(N,V,T) - S_idgas(N,V,T) = (<Upot> - Delta A)/T[BR]",
    "  Delta S(N,p,T) = S(N,p,T) - S_idgas(N,p,T) = (<Upot> + pV - Delta G)/T",
    "[PAR]",
    
    "When a second energy file is specified ([TT]-f2[tt]), a free energy",
    "difference is calculated dF = -kT ln < e ^ -(EB-EA)/kT >A ,",
    "where EA and EB are the energies from the first and second energy",
    "files, and the average is over the ensemble A. [BB]NOTE[bb] that",
    "the energies must both be calculated from the same trajectory."
    
  };
  static bool bSum=FALSE,bFee=FALSE,bAll=FALSE,bFluct=FALSE;
  static bool bDp=FALSE,bMutot=FALSE,bOrinst=FALSE,bOvec=FALSE;
  static int  skip=0,nmol=1,ndf=3;
  static real reftemp=300.0,ezero=0;
  t_pargs pa[] = {
    { "-fee",   FALSE, etBOOL,  {&bFee},
      "Do a free energy estimate" },
    { "-fetemp", FALSE, etREAL,{&reftemp},
      "Reference temperature for free energy calculation" },
    { "-zero", FALSE, etREAL, {&ezero},
      "Subtract a zero-point energy" },
    { "-sum",  FALSE, etBOOL, {&bSum},
      "Sum the energy terms selected rather than display them all" },
    { "-dp",   FALSE, etBOOL, {&bDp},
      "Print energies in high precision" },
    { "-mutot",FALSE, etBOOL, {&bMutot},
      "Compute the total dipole moment from the components" },
    { "-skip", FALSE, etINT,  {&skip},
      "Skip number of frames between data points" },
    { "-aver", FALSE, etBOOL, {&bAll},
      "Print also the X1,t and sigma1,t, only if only 1 energy is requested" },
    { "-nmol", FALSE, etINT,  {&nmol},
      "Number of molecules in your sample: the energies are divided by this number" },
    { "-ndf",  FALSE, etINT,  {&ndf},
      "Number of degrees of freedom per molecule. Necessary for calculating the heat capacity" },
    { "-fluc", FALSE, etBOOL, {&bFluct},
      "Calculate autocorrelation of energy fluctuations rather than energy itself" },
    { "-orinst", FALSE, etBOOL, {&bOrinst},
      "Analyse instantaneous orientation data" },
    { "-ovec", FALSE, etBOOL, {&bOvec},
      "Also plot the eigenvectors with -oten" }
  };
  char *drleg[] = {
    "Running average",
    "Instantaneous"
  };
  static const char *setnm[] = {
    "Pres-XX", "Pres-XY", "Pres-XZ", "Pres-YX", "Pres-YY",
    "Pres-YZ", "Pres-ZX", "Pres-ZY", "Pres-ZZ", "Temperature",
    "Volume",  "Pressure"
  };
  
  FILE       *out,*fp_pairs=NULL,*fort=NULL,*fodt=NULL,*foten=NULL;
  FILE       **drout;
  ener_file_t fp;
  int        timecheck=0;
  gmx_mtop_t mtop;
  gmx_localtop_t *top=NULL;
  t_inputrec ir;
  t_energy   *ee_sum,**ee;
  enersum_t  enersum;
  gmx_enxnm_t *enm=NULL;
  t_enxframe *frame,*fr=NULL;
  int        cur=0;
#define NEXT (1-cur)
  int        nre,teller,teller_disre;
  gmx_large_int_t start_step,ee_nsteps,ee_nsum;
  int        nor=0,nex=0,norfr=0,enx_i=0;
  real       start_t;
  real       *bounds=NULL,*violaver=NULL,*oobs=NULL,*orient=NULL,*odrms=NULL;
  int        *index=NULL,*pair=NULL,norsel=0,*orsel=NULL,*or_label=NULL;
  int        nbounds=0,npairs;
  bool       bDisRe,bDRAll,bORA,bORT,bODA,bODR,bODT,bORIRE,bOTEN;
  bool       bFoundStart,bCont,bEDR,bVisco;
  double     sum,sumaver,sumt,ener,dbl;
  double     *time=NULL;
  real       **eneset=NULL, **enesum=NULL,Vaver;
  int        *set=NULL,i,j,k,nset,sss,nenergy;
  char       **pairleg,**odtleg,**otenleg;
  char       **leg=NULL;
  char       **nms;
  char       *anm_j,*anm_k,*resnm_j,*resnm_k;
  int        resnr_j,resnr_k;
  const char *orinst_sub = "@ subtitle \"instantaneous\"\n";
  char       buf[256];
  output_env_t oenv;
  t_filenm   fnm[] = {
    { efEDR, "-f",    NULL,      ffREAD  },
    { efEDR, "-f2",   NULL,      ffOPTRD },
    { efTPX, "-s",    NULL,      ffOPTRD },
    { efXVG, "-o",    "energy",  ffWRITE },
    { efXVG, "-viol", "violaver",ffOPTWR },
    { efXVG, "-pairs","pairs",   ffOPTWR },
    { efXVG, "-ora",  "orienta", ffOPTWR },
    { efXVG, "-ort",  "orientt", ffOPTWR },
    { efXVG, "-oda",  "orideva", ffOPTWR },
    { efXVG, "-odr",  "oridevr", ffOPTWR },
    { efXVG, "-odt",  "oridevt", ffOPTWR },
    { efXVG, "-oten", "oriten",  ffOPTWR },
    { efXVG, "-corr", "enecorr", ffOPTWR },
    { efXVG, "-vis",  "visco",   ffOPTWR },
    { efXVG, "-ravg", "runavgdf",ffOPTWR }
  };
#define NFILE asize(fnm)
  int     npargs;
  t_pargs *ppa;
  
  CopyRight(stderr,argv[0]);
  npargs = asize(pa);
  ppa    = add_acf_pargs(&npargs,pa);
  parse_common_args(&argc,argv,
                    PCA_CAN_VIEW | PCA_CAN_BEGIN | PCA_CAN_END | PCA_BE_NICE,
		    NFILE,fnm,npargs,ppa,asize(desc),desc,0,NULL,&oenv);
  
  bDRAll = opt2bSet("-pairs",NFILE,fnm);
  bDisRe = opt2bSet("-viol",NFILE,fnm) || bDRAll;
  bORA   = opt2bSet("-ora",NFILE,fnm);
  bORT   = opt2bSet("-ort",NFILE,fnm);
  bODA   = opt2bSet("-oda",NFILE,fnm);
  bODR   = opt2bSet("-odr",NFILE,fnm);
  bODT   = opt2bSet("-odt",NFILE,fnm);
  bORIRE = bORA || bORT || bODA || bODR || bODT;
  bOTEN  = opt2bSet("-oten",NFILE,fnm);

  nset = 0;

  snew(frame,2);
  fp = open_enx(ftp2fn(efEDR,NFILE,fnm),"r");
  do_enxnms(fp,&nre,&enm);

  /* Initiate energies and set them to zero */
  snew(ee_sum,nre);
  enersum.nframes = 0;
  snew(enersum.sum,nre);
  snew(enersum.sum2,nre);
  nenergy = 0;
  Vaver = -1;
  
  bVisco = opt2bSet("-vis",NFILE,fnm);
  
  if (!bDisRe) {
    if (bVisco) {
      nset=asize(setnm);
      snew(set,nset);
      /* This is nasty code... To extract Pres tensor, Volume and Temperature */
      for(j=0; j<nset; j++) {
	for(i=0; i<nre; i++) {
	  if (strstr(enm[i].name,setnm[j])) {
	    set[j]=i;
	    break;
	  }
	}
        if (i == nre) {
	  if (strcasecmp(setnm[j],"Volume")==0) {
	    printf("Enter the box volume (" unit_volume "): ");
	    if(1 != scanf("%lf",&dbl))
	    {
	      gmx_fatal(FARGS,"Error reading user input");
            }
	    Vaver = dbl;
	  } else
	    gmx_fatal(FARGS,"Could not find term %s for viscosity calculation",
			setnm[j]);
        }
      }
    }
    else {
      set=select_by_name(nre,enm,&nset);
    }
    /* Print all the different units once */
    sprintf(buf,"(%s)",enm[set[0]].unit);
    for(i=1; i<nset; i++) {
      for(j=0; j<i; j++) {
	if (strcmp(enm[set[i]].unit,enm[set[j]].unit) == 0) {
	  break;
	}
      }
      if (j == i) {
	strcat(buf,", (");
	strcat(buf,enm[set[i]].unit);
	strcat(buf,")");
      }
    }
    out=xvgropen(opt2fn("-o",NFILE,fnm),"Gromacs Energies","Time (ps)",buf,
                 oenv);
    
    snew(leg,nset+1);
    for(i=0; (i<nset); i++)
      leg[i] = enm[set[i]].name;
    if (bSum) {
      leg[nset]=strdup("Sum");
      xvgr_legend(out,nset+1,leg,oenv);
    }
    else
      xvgr_legend(out,nset,leg,oenv);
    
    snew(eneset,nset+1);
    if (bVisco)
      snew(enesum,nset+1);
    time = NULL;

    if (bORIRE || bOTEN)
      get_orires_parms(ftp2fn(efTPX,NFILE,fnm),&nor,&nex,&or_label,&oobs);
    
    if (bORIRE) {
      if (bOrinst)
	enx_i = enxORI;
      else
	enx_i = enxOR;

      if (bORA || bODA)
	snew(orient,nor);
      if (bODR)
	snew(odrms,nor);
      if (bORT || bODT) {
	fprintf(stderr,"Select the orientation restraint labels you want (-1 is all)\n");
	fprintf(stderr,"End your selection with 0\n");
	j = -1;
	orsel = NULL;
	do {
	  j++;
	  srenew(orsel,j+1);
	  if(1 != scanf("%d",&(orsel[j])))
	  {
	    gmx_fatal(FARGS,"Error reading user input");
	  }
	} while (orsel[j] > 0);
	if (orsel[0] == -1) {
	  fprintf(stderr,"Selecting all %d orientation restraints\n",nor);
	  norsel = nor;
	  srenew(orsel,nor);
	  for(i=0; i<nor; i++)
	    orsel[i] = i;
	} else {
	  /* Build the selection */
	  norsel=0;
	  for(i=0; i<j; i++) {
	    for(k=0; k<nor; k++)
	      if (or_label[k] == orsel[i]) {
		orsel[norsel] = k;
		norsel++;
		break;
	      }
	    if (k == nor)
	      fprintf(stderr,"Orientation restraint label %d not found\n",
		      orsel[i]);
	  }
	}
	snew(odtleg,norsel);
	for(i=0; i<norsel; i++) {
	  snew(odtleg[i],256);
	  sprintf(odtleg[i],"%d",or_label[orsel[i]]);
	}
	if (bORT) {
	  fort=xvgropen(opt2fn("-ort",NFILE,fnm), "Calculated orientations",
			"Time (ps)","",oenv);
	  if (bOrinst)
	    fprintf(fort,"%s",orinst_sub);
	  xvgr_legend(fort,norsel,odtleg,oenv);
	}
	if (bODT) {
	  fodt=xvgropen(opt2fn("-odt",NFILE,fnm),
			"Orientation restraint deviation",
			"Time (ps)","",oenv);
	  if (bOrinst)
	    fprintf(fodt,"%s",orinst_sub);
	  xvgr_legend(fodt,norsel,odtleg,oenv);
	}
      }
    }
    if (bOTEN) {
      foten=xvgropen(opt2fn("-oten",NFILE,fnm),
		     "Order tensor","Time (ps)","",oenv);
      snew(otenleg,bOvec ? nex*12 : nex*3);
      for(i=0; i<nex; i++) {
	for(j=0; j<3; j++) {
	  sprintf(buf,"eig%d",j+1);
	  otenleg[(bOvec ? 12 : 3)*i+j] = strdup(buf);
	}
	if (bOvec) {
	  for(j=0; j<9; j++) {
	    sprintf(buf,"vec%d%s",j/3+1,j%3==0 ? "x" : (j%3==1 ? "y" : "z"));
	    otenleg[12*i+3+j] = strdup(buf);
	  }
	}
      }
      xvgr_legend(foten,bOvec ? nex*12 : nex*3,otenleg,oenv);
    }
  }
  else {
    nbounds=get_bounds(ftp2fn(efTPX,NFILE,fnm),&bounds,&index,&pair,&npairs,
		       &mtop,&top,&ir);
    snew(violaver,npairs);
    out=xvgropen(opt2fn("-o",NFILE,fnm),"Sum of Violations",
		 "Time (ps)","nm",oenv);
    xvgr_legend(out,2,drleg,oenv);  
    if (bDRAll) { 
      fp_pairs=xvgropen(opt2fn("-pairs",NFILE,fnm),"Pair Distances",
			"Time (ps)","Distance (nm)",oenv);
      if (get_print_xvgr_codes(oenv))
	fprintf(fp_pairs,"@ subtitle \"averaged (tau=%g) and instantaneous\"\n",
		ir.dr_tau);
    }
  }
  
  /* Initiate counters */
  teller       = 0;
  teller_disre = 0;
  bFoundStart  = FALSE;
  start_step   = 0;
  start_t      = 0;
  ee_nsteps    = 0;
  ee_nsum      = 0;
  do {
    /* This loop searches for the first frame (when -b option is given), 
     * or when this has been found it reads just one energy frame
     */
    do {
      bCont = do_enx(fp,&(frame[NEXT]));
      
      if (bCont) {
	timecheck = check_times(frame[NEXT].t);
      }      
    } while (bCont && (timecheck < 0));
    
    if ((timecheck == 0) && bCont) {
      /* We read a valid frame, so we can use it */
      fr = &(frame[NEXT]);
      
      if (fr->nre > 0) {
	/* The frame contains energies, so update cur */
	cur  = NEXT;
	
	if (!bFoundStart) {
	  bFoundStart = TRUE;
	  /* Initiate the previous step data */
	  start_step = fr->step;
	  start_t    = fr->t;
	  /* Initiate the energy sums */
	  for(i=0; i<fr->nre; i++) {
	    ee_sum[i].esum = fr->ener[i].e;
	    ee_sum[i].eav  = 0;
	  }
	  ee_nsteps = 1;
	  ee_nsum   = 1;
	} else {
	  if (ee_nsum > 0) {
	    if (fr->step - start_step + 1 == ee_nsteps + fr->nsteps) {
	      if (fr->nsum <= 1) {
		/* We have a sum of a single frame:
		 * add the energy to the sums.
		 */
		for(i=0; i<fr->nre; i++) {
		  ee_sum[i].eav  +=
		    dsqr(ee_sum[i].esum/ee_nsum - (ee_sum[i].esum + fr->ener[i].e)/(ee_nsum + 1))*
		    ee_nsum*(ee_nsum + 1);
		  ee_sum[i].esum += fr->ener[i].e;
		}
		ee_nsum += 1;
	      } else {
		/* Add the sums to the total */
		for(i=0; i<fr->nre; i++) {
		  ee_sum[i].eav  +=
		    fr->ener[i].eav +
		    dsqr(ee_sum[i].esum/ee_nsum - (ee_sum[i].esum + fr->ener[i].esum)/(ee_nsum + fr->nsum))*
		    ee_nsum*(ee_nsum + fr->nsum)/(double)fr->nsum;
		ee_sum[i].esum += fr->ener[i].esum;
		}
		ee_nsum += fr->nsum;
	      }
	    } else {
	      /* The interval does not match fr->nsteps:
	       * can not do exact averages.
	       */
	      ee_nsum = 0;
	    }
	    ee_nsteps = fr->step - start_step + 1;
	  }
	}
      }
      /*
       * Define distance restraint legends. Can only be done after
       * the first frame has been read... (Then we know how many there are)
       */
      if (bDisRe && bDRAll && !leg && (fr->ndisre > 0)) {
	t_iatom   *fa;
	t_iparams *ip;
	
	fa = top->idef.il[F_DISRES].iatoms; 
	ip = top->idef.iparams;

	if (fr->ndisre != top->idef.il[F_DISRES].nr/3)
	  gmx_fatal(FARGS,"Number of disre pairs in the energy file (%d) does not match the number in the run input file (%d)\n",
		      fr->ndisre,top->idef.il[F_DISRES].nr/3);
	
	snew(pairleg,fr->ndisre);
	for(i=0; i<fr->ndisre; i++) {
	  snew(pairleg[i],30);
	  j=fa[3*i+1];
	  k=fa[3*i+2];
	  gmx_mtop_atominfo_global(&mtop,j,&anm_j,&resnr_j,&resnm_j);
	  gmx_mtop_atominfo_global(&mtop,k,&anm_k,&resnr_k,&resnm_k);
	  sprintf(pairleg[i],"%d %s %d %s (%d)",
		  resnr_j+1,anm_j,resnr_k+1,anm_k,
		  ip[fa[3*i]].disres.label);
	}
	set=select_it(fr->ndisre,pairleg,&nset);
	snew(leg,2*nset);
	for(i=0; (i<nset); i++) {
	  snew(leg[2*i],32);
	  sprintf(leg[2*i],  "a %s",pairleg[set[i]]);
	  snew(leg[2*i+1],32);
	  sprintf(leg[2*i+1],"i %s",pairleg[set[i]]);
	}
	xvgr_legend(fp_pairs,2*nset,leg,oenv);    
      }
      
      /* 
       * Store energies for analysis afterwards... 
       */
      if (!bDisRe && (fr->nre > 0)) {
	if ((nenergy % 1000) == 0) {
	  srenew(time,nenergy+1000);
	  for(i=0; (i<=nset); i++) {
	    srenew(eneset[i],nenergy+1000);
	    if (bVisco)
	      srenew(enesum[i],nenergy+1000);
	  }
	}
	time[nenergy] = fr->t;
	sum=0;
	for(i=0; (i<nset); i++) {
	  ener = fr->ener[set[i]].e;
	  eneset[i][nenergy] = ener;
	  sum += ener;
	  if (bVisco) {
	    enesum[i][nenergy] = fr->ener[set[i]].esum;
	  }
	  /* Sum the actual frame energies,
	   * for in case we do not have exact sums in the energy file.
	   */
	  enersum.sum[i]  += ener;
	  enersum.sum2[i] += ener*ener;
	}
	if (bSum) {
	  eneset[nset][nenergy] = sum;
	}
	nenergy++;
	enersum.nframes++;
      }
      /* 
       * Printing time, only when we do not want to skip frames
       */
      if (!skip || teller % skip == 0) {
	if (bDisRe) {
	  /*******************************************
	   * D I S T A N C E   R E S T R A I N T S  
	   *******************************************/
	  if (fr->ndisre > 0) {
	    print_time(out,fr->t);
	    if (violaver == NULL)
	      snew(violaver,fr->ndisre);
	    
	    /* Subtract bounds from distances, to calculate violations */
	    calc_violations(fr->disre_rt,fr->disre_rm3tav,
			    nbounds,pair,bounds,violaver,&sumt,&sumaver);

	    fprintf(out,"  %8.4f  %8.4f\n",sumaver,sumt);
	    if (bDRAll) {
	      print_time(fp_pairs,fr->t);
	      for(i=0; (i<nset); i++) {
		sss=set[i];
		fprintf(fp_pairs,"  %8.4f",
			mypow(fr->disre_rm3tav[sss],minthird));
		fprintf(fp_pairs,"  %8.4f",
			fr->disre_rt[sss]);
	      }
	      fprintf(fp_pairs,"\n");
	    }
	    teller_disre++;
	  }
	}
	/*******************************************
	 * E N E R G I E S
	 *******************************************/
	else {
	  if (fr->nre > 0) {
	    print_time(out,fr->t);
	    if (bSum) 
	      print1(out,bDp,(eneset[nset][nenergy-1])/nmol-ezero);
	    else if ((nset == 1) && bAll) {
	      print1(out,bDp,fr->ener[set[0]].e);
	      print1(out,bDp,fr->ener[set[0]].esum);
	      print1(out,bDp,fr->ener[set[0]].eav);
	    }
	    else for(i=0; (i<nset); i++)
	      print1(out,bDp,(fr->ener[set[i]].e)/nmol-ezero);

	    fprintf(out,"\n");
	  }
	  if (bORIRE && fr->nblock>enx_i && fr->nr[enx_i]>0) {
	    if (fr->nr[enx_i] != nor)
	      gmx_fatal(FARGS,"Number of orientation restraints in energy file (%d) does not match with the topology (%d)",fr->nr[enx_i],nor);
	    if (bORA || bODA)
	      for(i=0; i<nor; i++)
		orient[i] += fr->block[enx_i][i];
	    if (bODR)
	      for(i=0; i<nor; i++)
		odrms[i] += sqr(fr->block[enx_i][i]-oobs[i]);
	    if (bORT) {
	      fprintf(fort,"  %10f",fr->t);
	      for(i=0; i<norsel; i++)
		fprintf(fort," %g",fr->block[enx_i][orsel[i]]); 
	      fprintf(fort,"\n");
	    }
	    if (bODT) {
	      fprintf(fodt,"  %10f",fr->t);
	      for(i=0; i<norsel; i++)
		fprintf(fodt," %g",fr->block[enx_i][orsel[i]]-oobs[orsel[i]]); 
	      fprintf(fodt,"\n");
	    }
	    norfr++;
	  }
	  if (bOTEN && fr->nblock>enxORT) {
	    if (fr->nr[enxORT] != nex*12)
	      gmx_fatal(FARGS,"Number of orientation experiments in energy file (%g) does not match with the topology (%d)",fr->nr[enxORT]/12,nex);
	    fprintf(foten,"  %10f",fr->t);
	    for(i=0; i<nex; i++)
	      for(j=0; j<(bOvec?12:3); j++)
	      	fprintf(foten," %g",fr->block[enxORT][i*12+j]);
	    fprintf(foten,"\n");
	  }
	}
      }
    }
  } while (bCont && (timecheck == 0));
  
  fprintf(stderr,"\n");
  close_enx(fp);
  
  ffclose(out);

  if (bDRAll)
    ffclose(fp_pairs);

  if (bORT)
    ffclose(fort);
  if (bODT)
    ffclose(fodt);
  if (bORA) {
    out = xvgropen(opt2fn("-ora",NFILE,fnm),
		   "Average calculated orientations",
		   "Restraint label","",oenv);
    if (bOrinst)
      fprintf(out,"%s",orinst_sub);
    for(i=0; i<nor; i++)
      fprintf(out,"%5d  %g\n",or_label[i],orient[i]/norfr);
    ffclose(out);
  }
  if (bODA) {
    out = xvgropen(opt2fn("-oda",NFILE,fnm),
		   "Average restraint deviation",
		   "Restraint label","",oenv);
    if (bOrinst)
      fprintf(out,"%s",orinst_sub);
    for(i=0; i<nor; i++)
      fprintf(out,"%5d  %g\n",or_label[i],orient[i]/norfr-oobs[i]);
    ffclose(out);
  }
  if (bODR) {
    out = xvgropen(opt2fn("-odr",NFILE,fnm),
		   "RMS orientation restraint deviations",
		   "Restraint label","",oenv);
    if (bOrinst)
      fprintf(out,"%s",orinst_sub);
    for(i=0; i<nor; i++)
      fprintf(out,"%5d  %g\n",or_label[i],sqrt(odrms[i]/norfr));
    ffclose(out);
  }
  if (bOTEN)
    ffclose(foten);

  if (bDisRe) {
    analyse_disre(opt2fn("-viol",NFILE,fnm),
		  teller_disre,violaver,bounds,index,pair,nbounds,oenv);
  } else {
    analyse_ener(opt2bSet("-corr",NFILE,fnm),opt2fn("-corr",NFILE,fnm),
		 bFee,bSum,bFluct,bVisco,opt2fn("-vis",NFILE,fnm),
		 nmol,ndf,start_step,start_t,frame[cur].step,frame[cur].t,
		 time,reftemp,ee_nsum,ee_sum,&enersum,
		 nset,set,nenergy,eneset,enesum,leg,enm,Vaver,ezero,oenv);
  }
  if (opt2bSet("-f2",NFILE,fnm)) {
    fec(opt2fn("-f2",NFILE,fnm), opt2fn("-ravg",NFILE,fnm), 
	reftemp, nset, set, leg, nenergy, eneset, time ,oenv);
  }
  
  {
    const char *nxy = "-nxy";
    
    do_view(oenv,opt2fn("-o",NFILE,fnm),nxy);
    do_view(oenv,opt2fn_null("-ravg",NFILE,fnm),nxy);
    do_view(oenv,opt2fn_null("-ora",NFILE,fnm),nxy);
    do_view(oenv,opt2fn_null("-ort",NFILE,fnm),nxy);
    do_view(oenv,opt2fn_null("-oda",NFILE,fnm),nxy);
    do_view(oenv,opt2fn_null("-odr",NFILE,fnm),nxy);
    do_view(oenv,opt2fn_null("-odt",NFILE,fnm),nxy);
    do_view(oenv,opt2fn_null("-oten",NFILE,fnm),nxy);
  }
  thanx(stderr);
  
  return 0;
}
