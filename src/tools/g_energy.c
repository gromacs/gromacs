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
 * GROningen MAchine for Chemical Simulation
 */
static char *SRCID_g_energy_c = "$Id$";

#include <string.h>
#include <math.h>
#include "typedefs.h"
#include "fatal.h"
#include "vec.h"
#include "smalloc.h"
#include "enerio.h"
#include "statutil.h"
#include "assert.h"
#include "names.h"
#include "copyrite.h"
#include "macros.h"
#include "xvgr.h"
#include "gstat.h"
#include "physics.h"

static real       minthird=-1.0/3.0,minsixth=-1.0/6.0;

double mypow(double x,double y)
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
  
  printf("Select the terms you want from the following list\n");
  printf("End your selection with 0\n");

  if ( bVerbose ) {
    for(k=0; (k<nre); ) {
      for(j=0; (j<4) && (k<nre); j++,k++) 
	printf(" %3d=%14s",k+1,nm[k]);
      printf("\n");
    }
  }

  snew(bE,nre);
  do {
    scanf("%d",&n);
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

int get_bounds(char *topnm,real **bounds,int **dr_index)
{
  t_topology *top;
  t_functype *functype;
  t_iparams  *ip;
  int        i,j,k,type,ftype,natom;
  t_ilist    *disres;
  t_iatom    *iatom;
  real       *b;
  int        *ind;
  int        nb,index;
  
  top      = read_top(topnm);
  functype = top->idef.functype;
  ip       = top->idef.iparams;
  
  /* Count how many distance restraint there are... */
  nb=top->idef.il[F_DISRES].nr;
  if (nb == 0) {
    fprintf(stderr,"No distance restraints in topology!\n");
    exit(1);
  }
  
  /* Allocate memory */
  snew(b,nb);
  snew(ind,nb+1);
  
  /* Fill the bound array */
  nb=0;
  for(i=0; (i<top->idef.ntypes); i++) {
    ftype = functype[i];
    if (ftype == F_DISRES) {
      index = ip[i].disres.index;
      if (b[index] != 0.0)
	fprintf(stderr,"warning index %d occurs multiple times in topology."
		" i=%d, nb=%d\n",index,i,nb);
      else {
	b[nb] = ip[i].disres.rx0;
	nb++;
      }
    }
  }
  *bounds = b;
  
  /* Fill the index array */
  index   = -1;
  disres  = &(top->idef.il[F_DISRES]);
  iatom   = disres->iatoms;
  for(i=j=k=0; (i<disres->nr); ) {
    type  = iatom[i];
    ftype = top->idef.functype[type];
    natom = interaction_function[ftype].nratoms+1;
    if (index != top->idef.iparams[type].disres.index) {
      ind[j] = k;
      index  = top->idef.iparams[type].disres.index; 
      j ++;
    }
    k++;
    i += natom;
  }
  ind[j]  = k;
  assert(j == nb);
  
  *dr_index = ind;
  
  return nb;
}

void calc_violations(real rt[],real rav3[],int nb,int index[],
		     real bounds[],real viol[],double *st,double *sa)
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

static void analyse_disre(char *voutfn,    int nframes,
			  real violaver[], real bounds[],
			  int index[],     int nbounds)
{
  FILE   *vout;
  double sum,sumt,sumaver;
  int    i,j;
  
  /* Subtract bounds from distances, to calculate violations */
  calc_violations(violaver,violaver,
		  nbounds,index,bounds,NULL,&sumt,&sumaver);
  
#ifdef DEBUG
  fprintf(stderr,"\nSum of violations averaged over simulation: %g nm\n",
	  sumaver);
  fprintf(stderr,"Largest violation averaged over simulation: %g nm\n\n",
	  sumt);
#endif		    
  vout=xvgropen(voutfn,"r\\S-3\\N average violations","DR Index","nm");
  sum  = 0.0;
  sumt = 0.0;
  for(i=0; (i<nbounds); i++) {
    /* Do ensemble averaging */
    sumaver = 0;
    for(j=index[i]; (j<index[i+1]); j++) 
      sumaver += sqr(violaver[j]/nframes); 
    sumaver = max(0.0,mypow(sumaver,minsixth)-bounds[i]);
    
    sumt   += sumaver;
    sum     = max(sum,sumaver);
    fprintf(vout,"%10d  %10.5e\n",i,sumaver);
  }
#ifdef DEBUG
  for(j=0; (j<dr.ndr); j++)
    fprintf(vout,"%10d  %10.5e\n",j,mypow(violaver[j]/nframes,minthird));
#endif
  ffclose(vout);
  
  fprintf(stderr,"\nSum of violations averaged over simulation: %g nm\n",sumt);
  fprintf(stderr,"Largest violation averaged over simulation: %g nm\n\n",sum);
  
  xvgr_file(voutfn,"-graphtype bar");
}

static void einstein_visco(char *fn,int nsets,int nframes,real **set,
			   real V,real T,real time[])
{
  FILE *fp;
  real fac,**integral,dt,di,r_nf2;
  int  i,j,m,nf2;
  
  if (nframes < 1)
    return;
    
  dt  = (time[1]-time[0]);
  
  /* Store the average in integral[nsets] */
  nf2 = nframes/2;
  snew(integral,nsets+1);
  for(m=0; (m<nsets+1); m++)
    snew(integral[m],nf2);
    
  /* Use trapezium rule for integration and average over time origins */  
  for(j=0; (j<nf2); j++) {
    for(i=1; (i<nf2); i++) {
      for(m=0; (m<nsets); m++) {
	di   = 0.5*dt*(set[m][i+j]+set[m][i-1+j]);
	
	integral[m][i]     += di;
	integral[nsets][i] += di/nsets;
      }
    }
  }

  fp=xvgropen(fn,"Shear viscosity using Einstein relation","Time (ps)","cp");
  fac   = (V*1e-26)/(2*BOLTZMANN*T);
  r_nf2 = 1.0/nf2;
  for(i=0; (i<nf2); i++) {
    fprintf(fp,"%10g",time[i]);
    for(m=0; (m<=nsets); m++) 
      fprintf(fp,"  %10g",fac*sqr(r_nf2*integral[m][i]));
    fprintf(fp,"\n");
  }
  fclose(fp);
  for(m=0; (m<nsets+1); m++)
    sfree(integral[m]);
  sfree(integral);
}

static void analyse_ener(bool bCorr,char *corrfn,
			 bool bDG,bool bSum,bool bFluct,
			 bool bVisco,char *visfn,int  nmol,int ndf,
			 int  oldstep,int step,real oldt,real t,
			 real time[], real reftemp,
			 t_energy oldee[],t_energy ee[],
			 int nset,int set[],int nenergy,real **eneset,
			 char *leg[])
{
  FILE *fp;
  /* Check out the printed manual for equations! */
  real Dt,a,b,aver,avertot,stddev,delta_t,sigma,totaldrift;
  real xxx,integral,intBulk;
  real sfrac,oldfrac,diffsum,diffav,fstep,pr_aver,pr_stddev;
  double beta,expE,expEtot,*deltag;
  int  nsteps,iset;
  real x1m,x1mk,VarE=-1,Temp=-1,Pres=-1,VarV=-1,Vaver=-1,VarT=-1;
  int  i,j,m,k;
  char buf[256];
  
  nsteps  = step - oldstep;
  if (nsteps < 2) {
    fprintf(stderr,"Not enough steps (%d) for statistics\n",nsteps);
  }
  else {
    /* Calculate the time difference */
    delta_t = ((t - oldt)*nsteps)/(nsteps-1);
    
    fprintf(stderr,"\nStatistics over %d steps [ %g ps ], %d sets, %d energy frames\n\n",
	    nsteps,delta_t,nset,nenergy);
    
#ifdef DEBUG
    fprintf(stderr,"oldstep = %d, step = %d, oldt = %g, t = %g\n",
	    oldstep,step,oldt,t);
#endif
    fprintf(stderr,"%-35s  %10s  %10s  %10s  %10s",
	    "Energy","Average","RMS","Drift","Tot-Drift");
    if (bDG)
      fprintf(stderr,"  %10s  %10s\n","-Ln<e^(b E)>/b","Entropy");
    else
      fprintf(stderr,"\n");
    fprintf(stderr,"-------------------------------------------------------------------\n");
    
    if (bDG) {
      beta = 1.0/(BOLTZ*reftemp);
      snew(deltag,nset);
      if (bSum) {
	expEtot=0;
	avertot=0;
      }
    }
    for(i=0; (i<nset); i++) {
      iset = set[i];
      if (oldstep == 0) {
	aver   = ee[iset].esum/nsteps;
	stddev = sqrt(ee[iset].eav/nsteps);
      }
      else {
	m      = oldstep;
	k      = nsteps-1;
	aver   = (ee[iset].esum  - oldee[iset].esum)/k;
	fstep  = ((real) m) * ((real) (m+k))/((real) k);
	x1m    = oldee[iset].esum/m;
	x1mk   = ee[iset].esum/(m+k); 
	xxx    = sqr(x1m - x1mk);
	sigma  = ee[iset].eav - oldee[iset].eav - xxx * fstep;
	stddev = sqrt(sigma/k);
      }
      if (bDG) {
	expE = 0;
	for(j=0; (j<nenergy); j++) {
	  expE += exp(beta*(eneset[i][j]-aver)/nmol);
	}
	if (bSum) {
	  avertot+=aver;
	  expEtot+=expE/nenergy;
	}
	deltag[i] = log(expE/nenergy)/beta + aver/nmol;
      }
      if (strstr(leg[i],"Kin") != NULL) 
	VarE = sqr(stddev);
      else if (strstr(leg[i],"empera") != NULL) {
	VarT = sqr(stddev);
	Temp = aver;
      }
      else if (strstr(leg[i],"olum") != NULL) {
	VarV = sqr(stddev);
	Vaver= aver;
      }
      else if (strstr(leg[i],"essure") != NULL) {
	Pres = aver;
      }
      if (iset < F_TEMP) {
	pr_aver   = aver/nmol;
	pr_stddev = stddev/nmol;
      }
      else {
	pr_aver   = aver;
	pr_stddev = stddev;
      }
      lsq_y_ax_b(nenergy,time,eneset[i],&a,&b);
      totaldrift = a * delta_t;
      fprintf(stderr,"%-35s  %10g  %10g  %10g  %10g",
	      leg[i],pr_aver,pr_stddev,a,totaldrift);
      if (bDG) 
	fprintf(stderr,"  %10g  %10s\n",deltag[i],deltag[i]-pr_aver);
      else
	fprintf(stderr,"\n");
      if (bFluct) {
	for(j=0; (j<nenergy); j++)
	  eneset[i][j] -= aver;
      }
    }
    if (bSum) {
      fprintf(stderr,"%-35s  %10g  %10s  %10s  %10s",
	      "Total",avertot/nmol,"--","--","--");
      /* pr_aver,pr_stddev,a,totaldrift */
      if (bDG) 
	fprintf(stderr,"  %10g  %10g\n",
		log(expEtot)/beta + avertot/nmol,log(expEtot)/beta);
      else
	fprintf(stderr,"\n");
    }
    if (Temp != -1) {
      real factor;
      
      factor = nmol*ndf*VarT/(3.0*sqr(Temp));
      fprintf(stderr,"Heat Capacity Cv:   %10g J/mol K (factor = %g)\n",
	      1000*BOLTZ/(2.0/3.0 - factor),factor);
    }
    if ((VarV != -1) && (Temp != -1)) {
      real tmp = VarV/(Vaver*BOLTZ*Temp*PRESFAC);
      
      fprintf(stderr,"Isothermal Compressibility: %10g /bar\n",tmp);
      fprintf(stderr,"Adiabatic bulk modulus:     %10g  bar\n",1.0/tmp);
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
      }
      
      einstein_visco("evisje.xvg",3,nenergy,eneset,Vaver,Temp,time);
      
      /*do_autocorr(corrfn,buf,nenergy,3,eneset,Dt,eacNormal,TRUE,NULL,NULL);*/
      /* Do it for shear viscosity */
      strcpy(buf,"Shear Viscosity");
      low_do_autocorr(corrfn,buf,nenergy,3,eneset,Dt,
		      eacNormal,1,FALSE,TRUE,TRUE,FALSE,
		      "shear-fit.xvg","fitting",FALSE,0.0,0.0,0);
	
      /* Now for bulk viscosity */
      strcpy(buf,"Bulk Viscosity");
      low_do_autocorr(corrfn,buf,nenergy,1,&(eneset[11]),Dt,
		      eacNormal,1,FALSE,TRUE,TRUE,FALSE,
		      "bulk-fit.xvg","fitting",FALSE,0.0,0.0,0);
      
      factor = (Vaver*1e-26/(BOLTZMANN*Temp))*Dt;
      fp=xvgropen(visfn,buf,"Time (ps)","\\8h\\4 (cp)");
      xvgr_legend(fp,asize(leg),leg);
      
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
      do_autocorr(corrfn,buf,nenergy,nset,eneset,(delta_t/nenergy),
		  eacNormal,FALSE,NULL,NULL);
    }
  }
}

void print_one(FILE *fp,bool bDp,real e)
{
  if (bDp)
    fprintf(fp,"  %16.10e",e);
  else
    fprintf(fp,"  %10.4e",e);

}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_energy extracts energy components or distance restraint data",
    "from an energy file.",
    "The user is prompted to interactively select the energy terms",
    "she wants.[PAR]",
    "In the case of distance rstraints it is also possible to calculate",
    "the sum of violations as a time average and instantaneous.[PAR]",
    "Average and RMS are calculated with full precision from",
    "the simulation (see printed manual). Drift is calculated",
    "by performing a LSQ fit of the data to a straight line.",
    "Total drift is drift multiplied by total time."
  };
  static bool bDisRe=FALSE,bDRAll=FALSE;
  static bool bSum=FALSE,bDp=FALSE,bAll=FALSE,bFluct=FALSE;
  static bool bVisco=FALSE,bDG=FALSE;
  static int  skip=0,nmol=1,ndf=3;
  static real reftemp=300.0;
  t_pargs pa[] = {
    { "-rall", FALSE, etBOOL, &bDRAll,
      "Extract individual distance restraint data rather than energy terms" },
    { "-rsum", FALSE, etBOOL, &bDisRe,
      "Extract sum of violations (instantaneous and time-averaged)" },
    { "-dg",   FALSE, etBOOL, &bDG,
      "Do a free energy estimate" },
    { "-temp", FALSE, etREAL, &reftemp,
      "Reference temperature for free energy calculation" },
    { "-sum",  FALSE, etBOOL, &bSum,
      "Sum the energy terms selected rather than display them all" },
    { "-dp",   FALSE, etBOOL, &bDp,
      "Print energies in high precision (%16.10e)" },
    { "-skip", FALSE, etINT,  &skip,
      "Skip number of frames between data points" },
    { "-aver", FALSE, etBOOL, &bAll,
      "Print also the X1,t and sigma1,t, only if only 1 energy is requested" },
    { "-nmol", FALSE, etINT,  &nmol,
      "Number of molecules in your sample: the energies are divided by this number" },
    { "-ndf",  FALSE, etINT,  &ndf,
      "Number of degrees of freedom per molecule. Necessary for calculating the heat capacity" },
    { "-fluc", FALSE, etBOOL, &bFluct,
      "Calculate autocorrelation of energy fluctuations rather than energy itself" }
  };

  static char *drleg[] = {
    "Running average",
    "Instantaneous"
  };
  
  FILE       *in,*out;
  FILE       **drout;
  int        ndrout;
  int        timecheck;
  XDR        xdr;
  t_energy   *ee,*oldee;
  t_drblock  dr;
  int        teller=0,nre,step,oldstep;
  real       t,oldt;
  real       *bounds,*violaver=NULL;
  int        *index;
  int        nbounds;
  bool       bStarted,bCont,bEDR;
  double     sum,sumaver,sumt;
  real       **eneset,*time;
  int        *set,i,j,k,nset,sss,nenergy;
  char       **enm,**leg=NULL;
  char       **nms;
  t_filenm   fnm[] = {
    { efENE, "-f", NULL, ffOPTRD },
    { efEDR, "-d", NULL, ffOPTRD },
    { efTPB, "-s", NULL, ffOPTRD },
    { efXVG, "-o", NULL, ffWRITE },
    { efXVG, "-v",   "violaver", ffOPTWR },
    { efXVG, "-corr", "enecorr", ffOPTWR },
    { efXVG, "-vis",  "visco.xvg", ffOPTWR }
  };
#define NFILE asize(fnm)
  int     npargs;
  t_pargs *ppa;

  CopyRight(stderr,argv[0]);
  npargs = asize(pa);
  ppa    = add_acf_pargs(&npargs,pa);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,npargs,ppa,asize(desc),desc,0,NULL);
  if (bDRAll)
    bDisRe=TRUE;

  bEDR=ftp2bSet(efEDR,NFILE,fnm);
  if (bEDR) {
    xdropen(&xdr,ftp2fn(efEDR,NFILE,fnm),"r");
    enm=NULL;
    edr_nms(&xdr,&nre,&enm);
  }
  else {
    in=ftp2FILE(efENE,NFILE,fnm,"r");
    rd_ener_nms(in,&nre,&enm);
  }
  
  if (nre == 0) {
    fprintf(stderr,"No energies!\n");
    exit(1);
  }
  
  dr.ndr=0;  
  snew(ee,nre);
  snew(oldee,nre);
  nenergy = 0;

  bVisco = opt2bSet("-vis",NFILE,fnm);

  if (!bDisRe) {
    if (bVisco) {
      static char *setnm[] = {
	"Pres-XX", "Pres-XY", "Pres-XZ", "Pres-YX", "Pres-YY",
	"Pres-YZ", "Pres-ZX", "Pres-ZY", "Pres-ZZ", "Temperature",
	"Volume",  "Pressure"
      };
      nset=asize(setnm);
      snew(set,nset);
      /* This is nasty code... To extract Pres tensor, Volume and Temperature */
      for(j=0; (j<nset); j++) {
	for(i=0; (i<nre); i++) {
	  if (strstr(enm[i],setnm[j])) {
	    set[j]=i;
	    break;
	  }
	}
	if (i == nre)
	  fatal_error(0,"Could not find term %s for viscosity calculation",
		      setnm[j]);
      }
    }
    else {
      set=select_it(nre,enm,&nset);
    }
    out=xvgropen(opt2fn("-o",NFILE,fnm),"Gromacs Energies","Time (ps)",
		 "E (kJ mol\\S-1\\N)");

    snew(leg,nset+1);
    for(i=0; (i<nset); i++)
      leg[i]=enm[set[i]];
    if (bSum) {
      leg[nset]="Sum";
      xvgr_legend(out,nset+1,leg);
    }
    else
      xvgr_legend(out,nset,leg);
      
    snew(eneset,nset+1);
    time = NULL;
  }
  else {
    nbounds=get_bounds(ftp2fn(efTPB,NFILE,fnm),&bounds,&index);
    snew(violaver,nbounds);
    out=xvgropen(opt2fn("-o",NFILE,fnm),"Sum of Violations","Time (ps)","nm");
    if (!bDRAll)
      xvgr_legend(out,2,drleg);    
  }
  step     = 0;
  oldstep  = -1;
  t        = 0;
  oldt     = 0;
  do {
    do {
      if (bEDR)
	bCont=edr_io(&xdr,&t,&step,&nre,ee,&dr);
      else
	bCont=rd_ener(in,&t,&step,ee,&dr);
	
      if (bCont) {
	timecheck=check_times(t);
      
	/* It is necessary for statistics to start counting from 1 */
	step += 1; 
      }
      
    } while (bCont && (timecheck < 0));
    
    if (timecheck == 0) {
      if (oldstep == -1) {
	/* Initiate the previous step data 
	 * Since step ia incremented after reading, it is always >= 1,
	 * therefore this will be executed only once.
	 */
	oldstep = step - 1;
	oldt    = t;
	if (oldstep > 0)
	  for(i=0; (i<nre); i++) {
	    oldee[i].esum = ee[i].esum - ee[i].e;
	    oldee[i].eav  = ee[i].eav  - (sqr(oldee[i].esum - oldstep*ee[i].e)/
					  (oldstep*step));
	  }
      }
      
      if (bDisRe && bDRAll && !leg) {
	snew(leg,dr.ndr*2);
	for(i=0,k=dr.ndr; (i<dr.ndr); i++,k++) {
	  snew(leg[i],12);
	  sprintf(leg[i],  "<r>  %2d",i);
	  snew(leg[k],12);
	  sprintf(leg[k],"r(t) %2d",i);
	}
	set=select_it(dr.ndr*2,leg,&nset);
	for(i=0; (i<nset); i++)
	  leg[i]=leg[set[i]];
	xvgr_legend(out,nset,leg);    
      }
      
#define DONTSKIP(cnt) (skip) ? ((cnt % skip) == 0) : TRUE
      
      if (bCont) {
	/*******************************************
	 * D I S T A N C E   R E S T R A I N T S  
	 *******************************************/
	if (!bDisRe) {
	  if ((nenergy % 1000) == 0) {
	    srenew(time,nenergy+1000);
	    for(i=0; (i<=nset); i++)
	      srenew(eneset[i],nenergy+1000);
	  }
	  time[nenergy] = t;
	  sum=0;
	  for(i=0; (i<nset); i++) {
	    eneset[i][nenergy] = ee[set[i]].e;
	    sum+=ee[set[i]].e;
	  }
	  if (bSum) 
	    eneset[nset][nenergy] = sum;
	  nenergy++;
	}
	if (DONTSKIP(teller)) {
	  print_one(out,bDp,t);

	  if (bDisRe) {
	    if (violaver == NULL)
	      snew(violaver,dr.ndr);
	    
	    /* Subtract bounds from distances, to calculate violations */
	    calc_violations(dr.rt,dr.rav,
			    nbounds,index,bounds,violaver,&sumt,&sumaver);
	    
	    if (bDRAll) {
	      for(i=0; (i<nset); i++) {
		sss=set[i];
		if (sss < dr.ndr)
		  fprintf(out,"  %8.4f",dr.rav[sss]);
		else
		  fprintf(out,"  %8.4f",dr.rt[sss-dr.ndr]);
	      }
	    }
	    else {
	      fprintf(out,"  %8.4f  %8.4f",sumaver,sumt);
	    }
	  }
	  /*******************************************
	   * E N E R G I E S
	   *******************************************/
	  else {
	    if (bSum) 
	      print_one(out,bDp,eneset[nset][nenergy-1]/nmol);
	    else if ((nset == 1) && bAll) {
	      print_one(out,bDp,ee[set[0]].e);
	      print_one(out,bDp,ee[set[0]].esum);
	      print_one(out,bDp,ee[set[0]].eav);
	    }
	    else for(i=0; (i<nset); i++)
	      print_one(out,bDp,ee[set[i]].e);
	  }
	  fprintf(out,"\n");
	}
      }
      teller++;
    }
  } while (bCont && (timecheck == 0));
  
  fprintf(stderr,"\n");
  if (bEDR)
    xdrclose(&xdr);
  else
    ffclose(in);
    
  ffclose(out);

  if (bDisRe) 
    analyse_disre(opt2fn("-v",NFILE,fnm),teller,violaver,bounds,index,nbounds);
  else 
    analyse_ener(opt2bSet("-corr",NFILE,fnm),opt2fn("-corr",NFILE,fnm),
		 bDG,bSum,bFluct,bVisco,opt2fn("-vis",NFILE,fnm),
		 nmol,ndf,oldstep,step,oldt,t,time,reftemp,
		 oldee,ee,nset,set,nenergy,eneset,leg);
  
  xvgr_file(opt2fn("-o",NFILE,fnm),"-nxy");
    
  thanx(stdout);
  
  return 0;
}
