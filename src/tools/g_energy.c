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
static char *SRCID_g_energy_c = "$Id$";

#include <string.h>
#include <math.h>
#include "typedefs.h"
#include "fatal.h"
#include "vec.h"
#include "smalloc.h"
#include "enxio.h"
#include "statutil.h"
#include "assert.h"
#include "names.h"
#include "copyrite.h"
#include "macros.h"
#include "xvgr.h"
#include "gstat.h"
#include "physics.h"
#include "tpxio.h"

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

int get_bounds(char *topnm,real **bounds,int **index,int **dr_pair,int *npairs,
	       t_topology *top,t_inputrec *ir)
{
  t_functype *functype;
  t_iparams  *ip;
  int        natoms,i,j,k,type,ftype,natom;
  t_ilist    *disres;
  t_iatom    *iatom;
  real       *b,t;
  int        *ind,*pair;
  int        nb,index1;
  matrix     box;

  read_tpx(topnm,&i,&t,&t,ir,box,&natoms,NULL,NULL,NULL,top);

  functype = top->idef.functype;
  ip       = top->idef.iparams;
  
  /* Count how many distance restraint there are... */
  nb=top->idef.il[F_DISRES].nr;
  if (nb == 0)
    fatal_error(0,"No distance restraints in topology!\n");
  
  /* Allocate memory */
  snew(b,nb);
  snew(ind,nb);
  snew(pair,nb+1);
  
  /* Fill the bound array */
  nb=0;
  for(i=0; (i<top->idef.ntypes); i++) {
    ftype = functype[i];
    if (ftype == F_DISRES) {
      index1 = ip[i].disres.index;
      b[nb]   = ip[i].disres.up1;
      ind[nb] = index1;
      nb++;
    }
  }
  *bounds = b;
  
  /* Fill the index array */
  index1  = -1;
  disres  = &(top->idef.il[F_DISRES]);
  iatom   = disres->iatoms;
  for(i=j=k=0; (i<disres->nr); ) {
    type  = iatom[i];
    ftype = top->idef.functype[type];
    natom = interaction_function[ftype].nratoms+1;
    if (index1 != top->idef.iparams[type].disres.index) {
      pair[j] = k;
      index1  = top->idef.iparams[type].disres.index; 
      j ++;
    }
    k++;
    i += natom;
  }
  pair[j]  = k;
  *npairs = k;
  assert(j == nb);

  *index   = ind;
  *dr_pair = pair;
  
  return nb;
}

void calc_violations(real rt[],real rav3[],int nb,int index[],
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

static void analyse_disre(char *voutfn,    int nframes,
			  real violaver[], real bounds[], int index[],
			  int pair[],      int nbounds)
{
  FILE   *vout;
  double sum,sumt,sumaver;
  int    i,j;
  
  /* Subtract bounds from distances, to calculate violations */
  calc_violations(violaver,violaver,
		  nbounds,pair,bounds,NULL,&sumt,&sumaver);
  
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
  
  fprintf(stderr,"\nSum of violations averaged over simulation: %g nm\n",sumt);
  fprintf(stderr,"Largest violation averaged over simulation: %g nm\n\n",sum);
  
  xvgr_file(voutfn,"-graphtype bar");
}

static void einstein_visco(char *fn,char *fni,int nsets,int nframes,real **sum,
			   real V,real T,int nsteps,real time[])
{
  FILE *fp0,*fp1;
  real av[4],avold[4];
  real fac,dt,di;
  int  i,j,m,nf4;
  
  if (nframes < 1)
    return;
    
  dt  = (time[1]-time[0]);
  nf4 = nframes/4+1;
  
  for(i=0; i<=nsets; i++)
    avold[i] = 0;
  fp0=xvgropen(fni,"Shear viscosity integral",
	       "Time (ps)","(kg m\\S-1\\N s\\S-1\\N ps)");
  fp1=xvgropen(fn,"Shear viscosity using Einstein relation",
	       "Time (ps)","(kg m\\S-1\\N s\\S-1\\N)");
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

static void analyse_ener(bool bCorr,char *corrfn,
			 bool bGibbs,bool bSum,bool bFluct,
			 bool bVisco,char *visfn,int  nmol,int ndf,
			 int  oldstep,real oldt,int step,real t,
			 real time[], real reftemp,
			 t_energy oldee[],t_energy ee[],
			 int nset,int set[],int nenergy,real **eneset,
			 real **enesum,
			 char *leg[],real Vaver)
{
  FILE *fp;
  /* Check out the printed manual for equations! */
  real Dt,a,b,aver,avertot,stddev,delta_t,sigma,totaldrift;
  real xxx,integral,intBulk;
  real sfrac,oldfrac,diffsum,diffav,fstep,pr_aver,pr_stddev,fluct2;
  double beta=0,expE,expEtot,*gibbs=NULL;
  int  nsteps,iset;
  real x1m,x1mk,Temp=-1,Pres=-1,VarV=-1,VarT=-1;
  int  i,j,m,k;
  char buf[256];

  nsteps  = step - oldstep;
#ifdef DEBUG
  fprintf(stderr,"oldstep: %d, oldt: %g, step: %d, t: %g, nenergy: %d\n",
	  oldstep,oldt,step,t,nenergy);
#endif
  if (nsteps < 2) {
    fprintf(stderr,"Not enough steps (%d) for statistics\n",nsteps);
  }
  else {
    /* Calculate the time difference */
    delta_t = t - oldt;
    
    fprintf(stderr,"\nStatistics over %d steps [ %.4f thru %.4f ps ], %d data sets\n\n",
	    nsteps,oldt,t,nset);
    
    fprintf(stderr,"%-24s %10s %10s %10s %10s %10s",
	    "Energy","Average","RMSD","Fluct.","Drift","Tot-Drift");
    if (bGibbs)
      fprintf(stderr,"  %10s\n","-ln<e^(E/kT)>*kT");
    else
      fprintf(stderr,"\n");
    fprintf(stderr,"-------------------------------------------------------------------------------\n");
    
    /* Initiate locals, only used with -sum */
    avertot=0;
    expEtot=0;
    if (bGibbs) {
      beta = 1.0/(BOLTZ*reftemp);
      snew(gibbs,nset);
    }
    for(i=0; (i<nset); i++) {
      iset = set[i];
      m      = oldstep;
      k      = nsteps;
#ifdef DEBUG
      fprintf(stderr,"sum: %g, oldsum: %g, k: %d\n",
	      ee[iset].esum,oldee[iset].esum,k);
#endif
      aver   = (ee[iset].esum  - oldee[iset].esum)/k;
      fstep  = ((real) m) * ((real) (m+k))/((real) k);
      x1m    = (m > 0) ? oldee[iset].esum/m : 0.0;
      x1mk   = ee[iset].esum/(m+k); 
      xxx    = sqr(x1m - x1mk);
      sigma  = ee[iset].eav - oldee[iset].eav - xxx * fstep;
      if((sigma/k)<GMX_REAL_EPS)
	sigma=0;
      stddev = sqrt(sigma/k);
      
      if (bSum) 
	avertot+=aver;
      if (bGibbs) {
	expE = 0;
	for(j=0; (j<nenergy); j++) {
	  expE += exp(beta*(eneset[i][j]-aver)/nmol);
	}
	if (bSum) 
	  expEtot+=expE/nenergy;
	
	gibbs[i] = log(expE/nenergy)/beta + aver/nmol;
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
      if (iset < F_TEMP) {
	pr_aver   = aver/nmol;
	pr_stddev = stddev/nmol;
      }
      else {
	pr_aver   = aver;
	pr_stddev = stddev;
      }
      lsq_y_ax_b(nenergy,time,eneset[i],&a,&b);
      if(fabs(a)<GMX_REAL_EPS)
	a=0;
      totaldrift = a * delta_t * (nsteps+1)/nsteps;
      fluct2 = sqr(pr_stddev) - sqr(totaldrift)/12;
      if (fluct2 < 0)
	fluct2 = 0;
      fprintf(stderr,"%-24s %10g %10g %10g %10g %10g",
	      leg[i],pr_aver,pr_stddev,sqrt(fluct2),a,totaldrift);
      if (bGibbs) 
	fprintf(stderr,"  %10g\n",gibbs[i]);
      else
	fprintf(stderr,"\n");
      if (bFluct) {
	for(j=0; (j<nenergy); j++)
	  eneset[i][j] -= aver;
      }
    }
    if (bSum) {
      fprintf(stderr,"%-24s %10g %10s %10s %10s %10s",
	      "Total",avertot/nmol,"--","--","--","--");
      /* pr_aver,pr_stddev,a,totaldrift */
      if (bGibbs) 
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
	enesum[0][i] = 0.5*(enesum[1][i]+enesum[3][i]);
	enesum[1][i] = 0.5*(enesum[2][i]+enesum[6][i]);
	enesum[2][i] = 0.5*(enesum[5][i]+enesum[7][i]);
      }
      
      einstein_visco("evisco.xvg","eviscoi.xvg",
		     3,nenergy,enesum,Vaver,Temp,nsteps,time);
      
      /*do_autocorr(corrfn,buf,nenergy,3,eneset,Dt,eacNormal,TRUE);*/
      /* Do it for shear viscosity */
      strcpy(buf,"Shear Viscosity");
      low_do_autocorr(corrfn,buf,nenergy,3,(nenergy+1)/2,eneset,Dt,
		      eacNormal,1,TRUE,TRUE,FALSE,FALSE,0.0,0.0,0);
	
      /* Now for bulk viscosity */
      strcpy(buf,"Bulk Viscosity");
      low_do_autocorr(corrfn,buf,nenergy,1,(nenergy+1)/2,&(eneset[11]),Dt,
		      eacNormal,1,TRUE,TRUE,FALSE,FALSE,0.0,0.0,0);
      
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
		  eacNormal,FALSE);
    }
  }
}

void print_one(FILE *fp,bool bDp,real e)
{
  if (bDp)
    fprintf(fp,"  %16g",e);
  else
    fprintf(fp,"  %10g",e);

}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    
    "g_energy extracts energy components or distance restraint",
    "data from an energy file. The user is prompted to interactively",
    "select the energy terms she wants.[PAR]",
    
    "When the [TT]-viol[tt] option is set, the time averaged",
    "violations are plotted and the running time-averaged and",
    "instantaneous sum of violations are recalculated. Additionally",
    "running time-averaged and instantaneous distances between",
    "selected pairs can be plotted with the [TT]-pairs[tt] option.[PAR]",
    
    "Average and RMSD are calculated with full precision from the",
    "simulation (see printed manual). Drift is calculated by performing",
    "a LSQ fit of the data to a straight line. Total drift is drift",
    "multiplied by total time.[PAR]",
    
    "With [TT]-G[tt] a Gibbs free energy estimate is calculated using",
    "the formula: G = -ln < e ^ (E/kT) > * kT, where k is Boltzmann's",
    "constant, T is set by [TT]-Gtemp[tt] and the average is over the",
    "ensemble (or time in a trajectory). Note that this is in principle",
    "only correct when averaging over the whole (Boltzmann) ensemble",
    "and using the potential energy. This also allows for an entropy",
    "estimate using G = H - T S, where H is the enthalpy (H = U + p V)",
    "and S entropy."
  };
  static bool bSum=FALSE,bGibbs=FALSE,bAll=FALSE,bFluct=FALSE;
  static bool bDp=FALSE,bMutot=FALSE;
  static int  skip=0,nmol=1,ndf=3;
  static real reftemp=300.0,ezero=0;
  t_pargs pa[] = {
    { "-G",   FALSE, etBOOL,  {&bGibbs},
      "Do a free energy estimate" },
    { "-Gtemp", FALSE, etREAL,{&reftemp},
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
      "Calculate autocorrelation of energy fluctuations rather than energy itself" }
  };
  static char *drleg[] = {
    "Running average",
    "Instantaneous"
  };
  static char *setnm[] = {
    "Pres-XX", "Pres-XY", "Pres-XZ", "Pres-YX", "Pres-YY",
    "Pres-YZ", "Pres-ZX", "Pres-ZY", "Pres-ZZ", "Temperature",
    "Volume",  "Pressure"
  };
  
  FILE       *out,*fp_pairs=NULL;
  FILE       **drout;
  int        fp;
  int        ndr;
  int        timecheck=0;
  t_topology top;
  t_inputrec ir;
  t_energy   *oldee,**ee;
  t_drblock  dr;
  int        teller,teller_disre,nre,step[2],oldstep;
  real       t[2],oldt;
  int        cur=0;
#define NEXT (1-cur)
  real       *bounds,*violaver=NULL;
  int        *index,*pair;
  int        nbounds=0,npairs;
  bool       bDisRe,bDRAll,bStarted,bCont,bEDR,bVisco;
  double     sum,sumaver,sumt,dbl;
  real       **eneset=NULL, **enesum=NULL,*time=NULL,Vaver;
  int        *set=NULL,i,j,k,nset,sss,nenergy;
  char       **enm=NULL,**leg=NULL,**pairleg;
  char       **nms;
  t_filenm   fnm[] = {
    { efENX, "-f",    NULL,      ffOPTRD },
    { efTPX, "-s",    NULL,      ffOPTRD },
    { efXVG, "-o",    "energy",  ffWRITE },
    { efXVG, "-viol", "violaver",ffOPTWR },
    { efXVG, "-pairs","pairs",   ffOPTWR },
    { efXVG, "-corr", "enecorr", ffOPTWR },
    { efXVG, "-vis",  "visco",   ffOPTWR }
  };
#define NFILE asize(fnm)
  int     npargs;
  t_pargs *ppa;
  
  CopyRight(stderr,argv[0]);
  npargs = asize(pa);
  ppa    = add_acf_pargs(&npargs,pa);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_BEGIN | PCA_CAN_END,TRUE,
		    NFILE,fnm,npargs,ppa,asize(desc),desc,0,NULL);
  
  bDRAll = opt2bSet("-pairs",NFILE,fnm);
  bDisRe = opt2bSet("-viol",NFILE,fnm) || bDRAll;
  
  fp = open_enx(ftp2fn(efENX,NFILE,fnm),"r");
  do_enxnms(fp,&nre,&enm);
  
  /* Initiate energies and set them to zero */
  dr.ndr=0;  
  snew(ee,2);
  snew(ee[cur],nre);
  snew(ee[NEXT],nre);
  snew(oldee,nre);
  nenergy = 0;
  Vaver = -1;
  
  bVisco = opt2bSet("-vis",NFILE,fnm);

  if (!bDisRe) {
    if (bVisco) {
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
	  if (strcmp(setnm[j],"Volume")==0) {
	    printf("Enter the box volume (nm^3): ");
	    scanf("%lf",&dbl);
	    Vaver = dbl;
	  } else
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
    if (bVisco)
      snew(enesum,nset+1);
    time = NULL;
  }
  else {
    nbounds=get_bounds(ftp2fn(efTPX,NFILE,fnm),&bounds,&index,&pair,&npairs,
		       &top,&ir);
    snew(violaver,npairs);
    out=xvgropen(opt2fn("-o",NFILE,fnm),"Sum of Violations",
		 "Time (ps)","nm");
    xvgr_legend(out,2,drleg);  
    if (bDRAll) { 
      fp_pairs=xvgropen(opt2fn("-pairs",NFILE,fnm),"Pair Distances",
			"Time (ps)","Distance (nm)");
      fprintf(fp_pairs,"@ subtitle \"averaged (tau=%g) and instantaneous\"\n",
	      ir.dr_tau);
    }
  }
  
  /* Initiate counters */
  teller       = 0;
  teller_disre = 0;
  step[cur]    = 0;
  t[cur]       = 0;
  oldstep      = -1;
  oldt         = 0;
  do {
    /* This loop searches for the first frame (when -b option is given), 
     * or when this has been found it reads just one energy frame
     */
    do {
      bCont = do_enx(fp,&(t[NEXT]),&(step[NEXT]),&nre,ee[NEXT],&ndr,&dr);
      
      if (bCont)
	timecheck = check_times(t[NEXT],t[NEXT]);
      
    } while (bCont && (timecheck < 0));
    
    if ((timecheck == 0) && bCont) {
      if (nre > 0) {
	/* 
	 * Only copy the values we just read, if we are within the time bounds
	 * It is necessary for statistics to start counting from 1 
	 */
	cur  = NEXT;
	step[cur]++;
	
	if (oldstep == -1) {
	  /* 
	   * Initiate the previous step data 
	   * Since step is incremented after reading, it is always >= 1,
	   * therefore this will be executed only once.
	   */
	  oldstep = step[cur] - 1;
	  oldt    = t[cur];
	  /* 
	   * If we did not start at the first step, oldstep will be > 0
	   * and we must initiate the data in the array. Otherwise it remains 0
	   */
	  if (oldstep > 0)
	    for(i=0; (i<nre); i++) {
	      oldee[i].esum = ee[cur][i].esum - ee[cur][i].e;
	      oldee[i].eav  = ee[cur][i].eav  - 
		(sqr(oldee[i].esum - oldstep*ee[cur][i].e)/(oldstep*step[cur]));
	    }
	}
      }
      /*
       * Define distance restraint legends. Can only be done after
       * the first frame has been read... (Then we know how many there are)
       */
      if (bDisRe && bDRAll && !leg && (ndr > 0)) {
	t_iatom   *fa;
	t_iparams *ip;
	
	fa = top.idef.il[F_DISRES].iatoms; 
	ip = top.idef.iparams;

	if (ndr != top.idef.il[F_DISRES].nr/3)
	  fatal_error(0,"Number of disre pairs in the energy file (%d) does not match the number in the run input file (%d)\n",
		      ndr,top.idef.il[F_DISRES].nr/3);
	
	snew(pairleg,dr.ndr);
	for(i=0; (i<dr.ndr); i++) {
	  snew(pairleg[i],30);
	  j=fa[3*i+1];
	  k=fa[3*i+2];
	  sprintf(pairleg[i],"%d %s %d %s (%d)",
		  top.atoms.atom[j].resnr+1,*top.atoms.atomname[j],
		  top.atoms.atom[k].resnr+1,*top.atoms.atomname[k],
		  ip[fa[3*i]].disres.index);
	}
	set=select_it(dr.ndr,pairleg,&nset);
	snew(leg,2*nset);
	for(i=0; (i<nset); i++) {
	  snew(leg[2*i],32);
	  sprintf(leg[2*i],  "a %s",pairleg[set[i]]);
	  snew(leg[2*i+1],32);
	  sprintf(leg[2*i+1],"i %s",pairleg[set[i]]);
	}
	xvgr_legend(fp_pairs,2*nset,leg);    
      }
      
      /* 
       * Store energies for analysis afterwards... 
       */
      if (!bDisRe && (nre > 0)) {
	if ((nenergy % 1000) == 0) {
	  srenew(time,nenergy+1000);
	  for(i=0; (i<=nset); i++) {
	    srenew(eneset[i],nenergy+1000);
	    if (bVisco)
	      srenew(enesum[i],nenergy+1000);
	  }
	}
	time[nenergy] = t[cur];
	sum=0;
	  for(i=0; (i<nset); i++) {
	    eneset[i][nenergy] = ee[cur][set[i]].e;
	    sum += ee[cur][set[i]].e;
	    if (bVisco)
	      enesum[i][nenergy] = ee[cur][set[i]].esum;
	  }
	if (bSum) 
	  eneset[nset][nenergy] = sum;
	nenergy++;
      }
      /* 
       * Printing time, only when we do not want to skip frames
       */
      if ((!skip) || ((teller % skip) == 0)) {
	if (bDisRe) {
	  /*******************************************
	   * D I S T A N C E   R E S T R A I N T S  
	   *******************************************/
	  if (ndr > 0) {
	    print_one(out,bDp,t[cur]);
	    if (violaver == NULL)
	      snew(violaver,dr.ndr);
	    
	    /* Subtract bounds from distances, to calculate violations */
	    calc_violations(dr.rt,dr.rav,
			    nbounds,pair,bounds,violaver,&sumt,&sumaver);

	    fprintf(out,"  %8.4f  %8.4f\n",sumaver,sumt);
	    if (bDRAll) {
	      print_one(fp_pairs,bDp,t[cur]);
	      for(i=0; (i<nset); i++) {
		sss=set[i];
		fprintf(fp_pairs,"  %8.4f",mypow(dr.rav[sss],minthird));
		fprintf(fp_pairs,"  %8.4f",dr.rt[sss]);
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
	  if (nre > 0) {
	    print_one(out,bDp,t[cur]);
	    if (bSum) 
	      print_one(out,bDp,(eneset[nset][nenergy-1]-ezero)/nmol);
	    else if ((nset == 1) && bAll) {
	      print_one(out,bDp,ee[cur][set[0]].e);
	      print_one(out,bDp,ee[cur][set[0]].esum);
	      print_one(out,bDp,ee[cur][set[0]].eav);
	    }
	    else for(i=0; (i<nset); i++)
	      print_one(out,bDp,(ee[cur][set[i]].e-ezero)/nmol);

	    fprintf(out,"\n");
	  }
	}
      }
      teller++;
    }
  } while (bCont && (timecheck == 0));
  
  fprintf(stderr,"\n");
  close_enx(fp);
  
  ffclose(out);
  if (bDRAll)
    ffclose(fp_pairs);

  if (bDisRe) 
    analyse_disre(opt2fn("-viol",NFILE,fnm),
		  teller_disre,violaver,bounds,index,pair,nbounds);
  else 
    analyse_ener(opt2bSet("-corr",NFILE,fnm),opt2fn("-corr",NFILE,fnm),
		 bGibbs,bSum,bFluct,bVisco,opt2fn("-vis",NFILE,fnm),
		 nmol,ndf,oldstep,oldt,step[cur],t[cur],time,reftemp,
		 oldee,ee[cur],nset,set,nenergy,eneset,enesum,leg,Vaver);
  
  xvgr_file(opt2fn("-o",NFILE,fnm),"-nxy");
    
  thanx(stdout);
  
  return 0;
}
