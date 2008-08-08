/*
 * $Id$
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *                        VERSION 3.0
 *
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 *
 * This program is free software; you can redistribute it and/or
 *
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
 *
 * finished FD 09/07/08
 *
 */

#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"
#include "xvgr.h"
#include "rmpbc.h"
#include "pbc.h"
#include "physics.h"
#include "index.h"


#define SQR(x) (pow(x,2.0))

static void precalc(t_topology top, t_trxframe fr,real mass2[]){

  real mtot;
  real qtot;
  int i,j,k,l;
  for(i=0;i<top.mols.nr;i++){
    k=top.mols.index[i];
    l=top.mols.index[i+1];
    mtot=0.0;
    qtot=0.0;
    for(j=k;j<l;j++){
      top.atoms.atom[j].q*=ENM2DEBYE;
      mtot+=top.atoms.atom[j].m;
      qtot+=top.atoms.atom[j].q;
    }
    for(j=k;j<l;j++){
      top.atoms.atom[j].q-=top.atoms.atom[j].m*qtot/mtot;
      mass2[j]=qtot*top.atoms.atom[j].m/mtot;
    }
  }
}


static void calc_mj(int ePBC, matrix avbox,int isize,int index0[],rvec fr[], rvec mj,real mass2[]){

  int   i,j,k,l;
  rvec  tmp;
  rvec  ran;
  rvec  rcat;
  t_pbc pbc;

  clear_rvec(ran);
  clear_rvec(rcat);
  clear_rvec(tmp);

  for(j=0;j<isize;j++){
    svmul(mass2[index0[j]],fr[index0[j]],tmp);
    (mass2[j]>0.0) ? rvec_inc(rcat,tmp):rvec_inc(ran,tmp);
  }

  svmul(-1.0,ran,ran);
  set_pbc(&pbc,ePBC,avbox);
  pbc_dx(&pbc,rcat,ran,mj);

}

static void remove_jump(matrix box,int natoms,rvec xp[],rvec x[]){

  rvec hbox;
  int d,i,m;

  for(d=0; d<DIM; d++)
    hbox[d] = 0.5*box[d][d];
  for(i=0; i<natoms; i++)
    for(m=DIM-1; m>=0; m--) {
      while (x[i][m]-xp[i][m] <= -hbox[m])
	for(d=0; d<=m; d++)
	  x[i][d] += box[m][d];
      while (x[i][m]-xp[i][m] > hbox[m])
	for(d=0; d<=m; d++)
	  x[i][d] -= box[m][d];
    }
}

static real calc_integral(bool bInt,FILE *outi,real *mc,real *time,int nsteps,real trust,real prefactor,real m2av,real mj2){

  real 	tmperr;
  real	averr;
  real	eps=0.0;
  real 	avfr;
  real	jca;
  real 	corint;
  real 	deltat;
  real 	int_trust=0.0;
  int  	itrust;
  int  	i;
  
  averr=0.0;
  jca=0.0;
  deltat=(time[nsteps-1]-time[0])/(real)(nsteps-1);
  corint=0.5*deltat*mc[0];
  itrust=(int)(trust*(real)(nsteps));
  eps=m2av-2.0*corint+mj2;
  eps*=prefactor;
  eps+=1.0;
  jca+=eps;
  if (bInt)
    fprintf(outi,"%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\n",time[0]-time[0],corint,eps,jca,averr);
  for (i=1;i<nsteps-1;i++){
    avfr=(real)(i+1);
    deltat=time[i]-time[i-1];
    
    corint+=deltat*mc[i];
    eps=m2av-2.0*corint+mj2;
    eps*=prefactor;
    eps+=1.0;
    
    tmperr=jca-(real)i*eps;
    tmperr=SQR(tmperr);
    tmperr/=avfr*avfr-avfr;
    averr+=tmperr;
    tmperr=averr/avfr;
    tmperr=sqrt(tmperr);
    jca+=eps;
    if (bInt)
      fprintf(outi,"%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\n",time[i]-time[0],corint,eps,jca/avfr,tmperr);
    if (itrust==i)
      int_trust=corint;
    
  }
  avfr=(real)(nsteps);
  deltat=time[nsteps-1]-time[nsteps-2];
  corint+=0.5*deltat*mc[nsteps-1];
  eps=m2av-2.0*corint+mj2;
  eps*=prefactor;
  eps+=1.0;
  tmperr=jca-(real)(nsteps-1)*eps;
  tmperr=SQR(tmperr);
  tmperr/=avfr*avfr-avfr;
  averr+=tmperr;
  jca+=eps;
  if (bInt)
    fprintf(outi,"%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\n",time[nsteps-1]-time[0],corint,eps,jca/avfr,sqrt(averr/avfr));
  
  if (itrust==nsteps)
    int_trust=corint;
  printf("\nTimestep: %f ps\n",deltat);
  printf("CorrelationIntegral (trust at %d) (complete) = %f %f\n",itrust,int_trust,corint);
  
  return int_trust;
}


static void do_corr(int n,int nsteps,int vfr[],real *jc,rvec *v0,rvec *v1,real *time,bool bCor,FILE *outf){

  int	j,k,i;
  int 	*avc;
  real 	tmperr;
  real	averr;
  real	jca;
  real 	avfr;

  i=0;
  snew(avc,nsteps);

  do{
    for(j=i;j<nsteps;j++){
      k=vfr[j];
      if(k>=0){
	jc[j-i]+=iprod(v0[i],v1[k]);
	avc[j-i]+=1;
      }
    }
    
    if ((i%1000)==0)
      printf("\r step %d / %d",i/n,nsteps/n);
    i+=n;
    
  }while (i<nsteps);

  printf("\r step %d / %d\n",nsteps/n,nsteps/n);
  
  i=0;
  averr=0.0;
  jca=0.0;
  
  if (bCor){
    while(i<nsteps){
      avfr=(real)(i+1);
      jc[i]/=(real)(avc[i]);
      tmperr=jca-(real)i*jc[i];
      tmperr=SQR(tmperr);
      if(i!=0) tmperr/=avfr*avfr-avfr;
      averr+=tmperr;
      tmperr=averr/avfr;
      tmperr=sqrt(tmperr);
      jca+=jc[i];
      fprintf(outf,"%8.3f\t%8.5f\t%8.5f\t%8.5f\n",time[i]-time[0],jc[i],jca/avfr,tmperr);
      i++;
    }
  }
}

static void dielectric(FILE *fmj,FILE *fmd,FILE *outi,FILE *outf,FILE *mcor,bool bInt,bool bCor,bool bTRR,int ePBC,t_topology top,t_trxframe fr,real temp,real trust,int n,int status,int isize, atom_id *index0,real mass2[],real eps_rf)
{
  int   i,j,k,l,f;
  int	nalloc,nfr,nvfr=0,m;
  int	*vfr=NULL;
  real	refr;
  real	*jc=NULL;
  real	*time=NULL;
  real	*mc=NULL;
  real	corint=0.0;
  real	prefactorav=0.0;
  real	prefactor=0.0;
  real	volume;
  real	volume_av=0.0;
  real	dk_s,dk_d;
  real	dm_s,dm_d;
  real  mj=0.0;
  real  mj2=0.0;
  real  mjd=0.0;
  real  mjdav=0.0;
  real	md2=0.0;
  real	mdav2=0.0;
  real	mderr=0.0;
  real	mjerr=0.0;
  real	tmperr;
  rvec 	**sx=NULL;
  rvec  mja_tmp;
  rvec  mj_tmp;
  rvec  mjd_tmp;
  rvec  mdvec;
  rvec	*mu=NULL;
  rvec	*xp=NULL;
  rvec	*v0=NULL;
  matrix  *sbox=NULL;
  matrix  avbox;
  matrix  corr;

  clear_mat(avbox);
  
  /* This is the main loop over frames */
  
  nfr = 0;
  nalloc = 0;
  
  do{
    
    if(nfr >= nalloc){
      nalloc+=100;
      srenew(sbox,nalloc);
      srenew(sx,nalloc);
      srenew(time,nalloc);
      srenew(mu,nalloc);
      srenew(vfr,nalloc);
    }
    
    copy_mat(fr.box,sbox[nfr]);
    m_add(avbox,fr.box,avbox);
    snew(sx[nfr],isize);
    for(i=0; i<isize;i++)
      copy_rvec(fr.x[index0[i]],sx[nfr][i]);
    
    
    
    if (fr.v!=NULL){
      srenew(v0,(nvfr+1));
      vfr[nfr]=nvfr;
      clear_rvec(v0[nvfr]);
      for(i=0;i<isize;i++){
	j=index0[i];
	svmul(mass2[j],fr.v[j],fr.v[j]);
	rvec_inc(v0[nvfr],fr.v[j]);
      }
      nvfr++;
    }
    else{
      vfr[nfr]=-1;
    }
    
    time[nfr]=fr.time;
    clear_rvec(mu[nfr]);
    if (xp)
      remove_jump(fr.box,fr.natoms,xp,fr.x);
    else{
      snew(xp,fr.natoms);
    }
    
    for(i=0; i<fr.natoms;i++){
      copy_rvec(fr.x[i],xp[i]);
    }
    
    for(i=0;i<isize;i++){
      j=index0[i];
      svmul(top.atoms.atom[j].q,fr.x[j],fr.x[j]);
      rvec_inc(mu[nfr],fr.x[j]);
    }

    volume = det(fr.box);
    volume_av += volume;
    
    nfr++;
    
  }while(read_next_frame(status,&fr));
  
  
  /* treatment of pressure scaling  and error estimation*/
  
  msmul(avbox,1.0/nfr,avbox);
  clear_rvec(mj_tmp);
  clear_rvec(mja_tmp);
  clear_rvec(mjd_tmp);
  clear_rvec(mdvec);
  for(f=0;f<nfr;f++){
    refr=(real)(f+1);
    if (f % 100 == 0)
      fprintf(stderr,"\rProcessing frame %d",f);
    
    m_inv(sbox[f],corr);
    mmul(avbox,corr,corr);
    for(i=0;i<isize;i++)
      mvmul(corr,sx[f][index0[i]],sx[f][index0[i]]);
    mvmul(corr,mu[f],mu[f]);
    calc_mj(ePBC,corr,isize,index0,sx[f],mj_tmp,mass2);
    rvec_inc(mja_tmp,mj_tmp);
    mjd+=iprod(mu[f],mj_tmp);
    rvec_inc(mdvec,mu[f]);
    tmperr=iprod(mu[f],mu[f]);
    tmperr=md2-(real)f*tmperr;
    tmperr*=tmperr;
    if(f!=0) tmperr/=refr*refr-refr;
    mderr+=tmperr;
    tmperr=iprod(mj_tmp,mj_tmp);
    tmperr=mj2-(real)f*tmperr;
    tmperr*=tmperr;
    if(f!=0) tmperr/=refr*refr-refr;
    mjerr+=tmperr;
    mj2+=iprod(mj_tmp,mj_tmp);
    md2+=iprod(mu[f],mu[f]);
    
    fprintf(fmj,"%8.3f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\n",time[f]-time[0],mj_tmp[XX],mj_tmp[YY],mj_tmp[ZZ],mj2/refr,sqrt(mjerr/refr));
    fprintf(fmd,"%8.3f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\n",time[f]-time[0],mu[f][XX],mu[f][YY],mu[f][ZZ],md2/refr,sqrt(mderr/refr));
  }
  
  mjd/=(real)nfr;
  md2/=(real)nfr;
  svmul(1.0/(real)nfr,mdvec,mdvec);
  mdav2=iprod(mdvec,mdvec);
  mjdav=iprod(mdvec,mja_tmp);
  mjdav/=(real)(nfr*nfr);
  mj2/=(real)nfr;
  mj=iprod(mja_tmp,mja_tmp);
  mj/=(real)(nfr*nfr);
  
  printf("\n\nAverage translational dipole moment M_J [D] after %d frames (|M|^2): %f %f %f (%f)\n",nfr,mja_tmp[XX],mja_tmp[YY],mja_tmp[ZZ],mj2);
  printf("\n\nAverage molecular dipole moment M_D [D] after %d frames (|M|^2): %f %f %f (%f)\n",nfr,mdvec[XX],mdvec[YY],mdvec[ZZ],md2);
  
  volume_av/=(real)nfr;
  prefactor=4.0*M_PI*1.0e-22/(9.0*volume_av*BOLTZMANN*temp);
  
  printf("\n\nAverage volume: %f nm^3\n\n",volume_av);
  
  if (v0!=NULL){
    snew(mc,nfr);
    printf("\nCalculating M_D - current correlation ... \n");
    do_corr(n,nfr,vfr,mc,mu,v0,time,bCor,mcor);
    corint=calc_integral(bInt,outi,mc,time,nfr,trust,prefactor,md2,mj2);
    printf("\nCalculating current autocorrelation ... \n");
    snew(jc,nfr);
    for(f=0;f<nfr;f++)
      svmul(1.0/ENM2DEBYE,v0[f],v0[f]);
    do_corr(n,nfr,vfr, jc, v0,v0,time,TRUE,outf);
    sfree(mc);
  }
  
  /* Calculation of the dielectric constant */
  
  dm_s=(md2-mdav2+mj2-mj+2.0*(mjd-mjdav));
  dm_d=md2-2.0*corint+mj2;
  if (eps_rf==0.0){
    dk_s=1.0+3.0*prefactor*dm_s;
    dk_d=1.0+3.0*prefactor*dm_d;
  }
  else{
    dk_s=2.0*eps_rf+1.0+6.0*eps_rf*prefactor*dm_s;
    dk_s/=2.0*eps_rf+1.0-3.0*prefactor*dm_s;
    dk_d=2.0*eps_rf+1.0+6.0*eps_rf*prefactor*dm_d;
    dk_d/=2.0*eps_rf+1.0-3.0*prefactor*dm_d;
  }
  
  printf("\nStatic dielectric constant using fluctuations (deltaM_D , deltaM_J, deltaM_JD) : %f (%f,%f,%f)\n",	 dk_s,md2-mdav2,mj2-mj,mjd-mjdav);
  if (v0)	  printf("\nStatic dielectric constant using integral     : %f\n",dk_d);
  
  sfree(mu);
  sfree(sx);
  sfree(sbox);
  sfree(time);
  if (v0 || jc){
    sfree(v0);
    sfree(jc);
  }
}


int gmx_current(int argc,char *argv[])
{

  static int  sh=10;
  static real temp=300.0;
  static real trust=0.25;
  static real eps_rf=0.0;
  t_pargs pa[] = {
    { "-sh", FALSE, etINT, {&sh},
      "Shift of the frames for averaging"},
    { "-eps", FALSE, etREAL, {&eps_rf},
     "Dielectric constant of the surrounding medium used for reaction field or Ewald summation, eps=0 corresponds to eps=infinity"},
    { "-tr", FALSE, etREAL, {&trust},
      "Fraction of the trajectory taken into account for the integral"},
    { "-temp", FALSE, etREAL, {&temp},
      "Temperature for calculating epsilon"
    }
  };

  t_topology top;
  char       title[STRLEN];
  char       **grpname=NULL;
  char	     *indexfn;
  t_trxframe fr;
  real       *mass2=NULL;
  rvec       *xtop,*vtop;
  matrix     box;
  atom_id    *index0=NULL;
  int        isize;
  int        status;
  int        flags = 0;
  bool	     bTop;
  bool	     bTRR;
  bool       bInt;
  bool       bCor;
  int	     ePBC=-1;
  int	     natoms;
  int 	     i,j,k=0,l;
  int        step;
  real	     t;
  real       lambda;
  FILE	     *outf=NULL;
  FILE       *outi=NULL;
  FILE       *tfile=NULL;
  FILE       *mcor=NULL;
  FILE       *fmj=NULL;
  FILE       *fmd=NULL;
  t_filenm fnm[] = {
    { efTPS,  NULL,  NULL, ffREAD },   /* this is for the topology */
    { efNDX, NULL, NULL, ffOPTRD },
    { efTRX, "-f", NULL, ffREAD },     /* and this for the trajectory */
    { efXVG, "-o", "caf.xvg", ffWRITE },
    { efXVG, "-mj", "mj.xvg", ffWRITE },
    { efXVG, "-md", "md.xvg", ffWRITE },
    { efXVG, "-ci", "cint.xvg", ffOPTWR }, /* for the output */
    { efXVG, "-mc", "mc.xvg", ffOPTWR }
  };

#define NFILE asize(fnm)


    static char *desc[] = {
    "This is a small tool for calculating the current autocorrelation function, the correlation",
    "of the rotational and translational dipole moment of the system, and the resulting static",
    "dielectric constant. To obtain a reasonable result the index group has to be neutral.",
    "[PAR]",
    "Options [TT]-rc[tt] and [TT]-tr[tt] are responsible for the averaging and integration of the",
    "autocorrelation functions. Since averaging proceeds by shifting the starting point",
    "through the trajectory, the shift can be modified to enable the choice of uncorrelated",
    "starting points. Towards the end, statistical inaccuracy grows and integrating the",
    "correlation function only yields reliable values until a certain point, depending on",
    "the number of frames. The option [TT]-rc[t] controls the region of the integral taken into account",
    "for calculating the static dielectric constant.",
    "[PAR]",
    "Option [TT]-temp[tt] sets the temperature required for the computation of the static dielectric constant.",
    "[PAR]",
    "Option [TT]-eps[tt] controls the dielectric constant of the surrounding medium for simulations using",
    "a reaction field or dipole corrections of the Ewald summation (eps=0 corresponds to",
    "tin-foil boundary conditions)."
  };





  /* At first the arguments will be parsed and the system information processed */


  CopyRight(stderr,argv[0]);

  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  bInt  = opt2bSet("-ci",NFILE,fnm);
  bCor  = opt2bSet("-mc",NFILE,fnm);

  bTop=read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xtop,&vtop,box,TRUE);

  sfree(xtop);
  sfree(vtop);
  indexfn = ftp2fn_null(efNDX,NFILE,fnm);
  snew(grpname,1);

  get_index(&(top.atoms),indexfn,1,&isize,&index0,grpname);

  flags = flags | TRX_READ_X | TRX_READ_V;

  read_first_frame(&status,ftp2fn(efTRX,NFILE,fnm),&fr,flags);

  snew(mass2,top.atoms.nr);

  precalc(top,fr,mass2);

  if (fr.v!=NULL){
    if (bInt){
      outi = xvgropen(opt2fn("-ci",NFILE,fnm),
		      "M\\sJ\\N-M\\sD\\N correlation Integral",xvgr_tlabel(),"CI (D\\S2\\N)");
      fprintf(outi,"# time\t integral\t average \t std.dev\n");
    }
    outf = xvgropen(opt2fn("-o",NFILE,fnm),
		    "Current autocorrelation function",xvgr_tlabel(),"ACF (e nm/ps)\\S2");
    fprintf(outf,"# time\t acf\t average \t std.dev\n");
    if (bCor){
      mcor = xvgropen(opt2fn("-mc",NFILE,fnm),
		      "M\\sD\\N - current  autocorrelation function",xvgr_tlabel(),
		      "< M\\sD\\N (0)\\c7\\CJ(t) >  (e nm/ps)\\S2");
      fprintf(mcor,"# time\t x\t y \t z \t average \t std.dev\n");
    }
    
  }
  
  fmj = xvgropen(opt2fn("-mj",NFILE,fnm),
		 "Averaged translational part of M",xvgr_tlabel(),"< M\\sJ\\N > (D)");
  fprintf(fmj,"# time\t x\t y \t z \t average of M_J^2 \t std.dev\n");
  fmd = xvgropen(opt2fn("-md",NFILE,fnm),
		 "Averaged rotational part of M",xvgr_tlabel(),"< M\\sD\\N > (D)");
  fprintf(fmd,"# time\t x\t y \t z \t average of M_D^2 \t std.dev\n");

  /* System information is read and prepared, dielectric() processes the frames and calculates the dielectric constant */

  dielectric(fmj,fmd,outi,outf,mcor,bInt,bCor,TRUE,ePBC,top,fr,temp,trust,sh,status,isize,index0,mass2,eps_rf);

  fclose(fmj);
  fclose(fmd);
  if (outf || mcor || outi){
    fclose(outf);
    if (bCor)
      fclose(mcor);
    if (bInt)
      fclose(outi);
  }

  thanx(stderr);

  return 0;
}
