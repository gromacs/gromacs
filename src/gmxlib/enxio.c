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
 * Great Red Owns Many ACres of Sand 
 */
static char *SRCID_enxio_c = "$Id$";
#include "futil.h"
#include "string2.h"
#include "fatal.h"
#include "smalloc.h"
#include "gmxfio.h"
#include "enxio.h"
	
typedef struct {
  real t;		/* Time frame				*/
  int  step;		/* MD step				*/
  int  nre;		/* Number of energies			*/
  int  ndisre;		/* Number of disre blocks		*/
  int  nuser;		/* User definable number		*/
  int  e_size;		/* Size (in bytes) of energies		*/
  int  d_size;		/* Size (in bytes) of disre blocks	*/
  int  u_size;		/* Size (in bytes) of user blocks	*/
} t_eheader;

static void wr_ener_nms(FILE *out,int nre,char *nms[])
{
  int i;
  
  fprintf(out,"%d\n",nre);
  for(i=0; (i<nre); i++)
    fprintf(out,"%s\n",nms[i]);
}

static void rd_ener_nms(FILE *in,int *nre,char ***nm)
{
  char line[256];
  int  i;
  
  fgets2(line,255,in);
  if (sscanf(line,"%d",nre) == 0) {
    *nre=0;
    return;
  }
  snew((*nm),*nre);
  for(i=0; (i< (*nre)); i++) {
    fgets2(line,255,in);
    trim(line);
    (*nm)[i]=strdup(line);
  }
}

static void edr_nms(int fp,int *nre,char ***nms)
{
  XDR  *xdr;
  bool bRead = fio_getread(fp);
  int  i;
  char **NM;

  xdr = fio_getxdr(fp);

  NM=*nms;
  
  if (!xdr_int(xdr,nre)) {
    *nre=0;
    return;
  }
  if (NM == NULL) {
    snew(NM,*nre);
  }
  for(i=0; i<*nre; i++) {
    if (bRead && NM[i]) {
      sfree(NM[i]);
      NM[i] = NULL;
    }
    xdr_string(xdr,&(NM[i]),STRLEN);
  }
  *nms=NM;
}

static bool do_eheader(int fp,t_eheader *eh,bool *bOK)
{
  bool bRead = fio_getread(fp);
  
  *bOK=TRUE;
  if (!do_real(eh->t))      return FALSE;
  if (!do_int (eh->step))   *bOK = FALSE;
  if (!do_int (eh->nre))    *bOK = FALSE;
  if (!do_int (eh->ndisre)) *bOK = FALSE;
  if (!do_int (eh->nuser))  *bOK = FALSE;
  if (!do_int (eh->e_size)) *bOK = FALSE;
  if (!do_int (eh->d_size)) *bOK = FALSE;
  if (!do_int (eh->u_size)) *bOK = FALSE;
  
  return *bOK;
}

void do_enxnms(int fp,int *nre,char ***nms)
{
  bool bRead;
  
  bRead = fio_getread(fp);
  if (fio_getftp(fp) == efEDR) {
    fio_select(fp);
    edr_nms(fp,nre,nms);
  }
  else if (bRead)
    rd_ener_nms(fio_getfp(fp),nre,nms);
  else
    wr_ener_nms(fio_getfp(fp),*nre,*nms);
}

void close_enx(int fp)
{
  fio_close(fp);
}

static bool empty_file(char *fn)
{
  FILE *fp;
  char dum;
  bool bEmpty;
  
  fp = ffopen(fn,"r");
  fread(&dum,sizeof(dum),1,fp);
  bEmpty = feof(fp);
  fclose(fp);

  return bEmpty;
}

static int framenr;

int open_enx(char *fn,char *mode)
{
  int       i,fp;
  char      **nm=NULL;
  int       nre;
  t_eheader *eh;
  bool bDum=TRUE;

  if (mode[0]=='r') {
    fp=fio_open(fn,mode);
    fio_select(fp);
    fio_setprecision(fp,FALSE);
    do_enxnms(fp,&nre,&nm);
    snew(eh,1);
    do_eheader(fp,eh,&bDum);
    
    /* Now check whether this file is in single precision */
    if (((eh->e_size && (eh->nre == nre) && 
	  (nre*4*sizeof(float) == eh->e_size)) ||
	 (eh->d_size && 
	  (eh->ndisre*sizeof(float)*2+sizeof(int) == eh->d_size)))){
      fprintf(stderr,"Opened %s as single precision energy file\n",fn);
      for(i=0; (i<nre); i++)
	sfree(nm[i]);
      sfree(nm);
    }
    else {
      fio_rewind(fp);
      fio_select(fp);
      fio_setprecision(fp,TRUE);
      do_enxnms(fp,&nre,&nm);
      do_eheader(fp,eh,&bDum);
      if (((eh->e_size && (eh->nre == nre) && 
	    (nre*4*sizeof(double) == eh->e_size)) ||
	   (eh->d_size && 
	    (eh->ndisre*sizeof(double)*2+sizeof(int) == eh->d_size))))
	fprintf(stderr,"Opened %s as double precision energy file\n",fn);
      else {
	if (empty_file(fn))
	  fatal_error(0,"File %s is empty",fn);
	else
	  fatal_error(0,"Energy file %s not recognized, maybe different CPU?",
		      fn);
      }
      for(i=0; (i<nre); i++)
	  sfree(nm[i]);
      sfree(nm);
    }
    sfree(eh);
    fio_rewind(fp);
  }
  else 
    fp = fio_open(fn,mode);
    
  framenr=0;
    
  return fp;
}

bool do_enx(int fp,real *t,int *step,int *nre,t_energy ener[],
	    int *ndr,t_drblock *drblock)
{
  int       i,ndisre=0,nuser=0;
  t_eheader eh;
  bool      bRead,bOK,bOK1;
  real      tmp1,tmp2;

  bOK = TRUE;
  bRead = fio_getread(fp);
  if (!bRead) {  
    eh.t      = *t;
    eh.step   = *step;
    eh.nre    = *nre;
    eh.ndisre = drblock ? (drblock->ndr) : 0;
    eh.nuser  = 0;
    eh.e_size = (*nre)*sizeof(ener[0].e)*4;
    eh.d_size = eh.ndisre ? 
      (drblock->ndr*(sizeof(drblock->rav[0])+sizeof(drblock->rt[0]))) : 0;
    eh.u_size = 0;
  }
  fio_select(fp);

  if (!do_eheader(fp,&eh,&bOK)) {
    if (bRead) {
	fprintf(stderr,"\rLast frame read %d                          ",framenr-1);
	if (!bOK)
	  fprintf(stderr,"\nWARNING: Incomplete frame: nr %6d time %8.3f\n",
		  framenr,eh.t);
    }
    return FALSE;
  }
  if (bRead) {
    if ( ( framenr<10 ) || ( framenr%10 == 0) ) 
      fprintf(stderr,"\rReading frame %6d time %8.3f           ",framenr,eh.t);
    framenr++;
  }
  /* Check sanity of this header */
  if ((eh.step >= 0) && (eh.nre > 0 || eh.ndisre > 0) && (eh.nuser == 0)) {
    *t     = eh.t;
    *step  = eh.step;
    *nre   = eh.nre;
    ndisre = eh.ndisre;
    nuser  = 0;
  }
  else {
    fprintf(stderr,"\nWARNING: there may be something wrong with energy file %s\n",
	    fio_getname(fp));
    fprintf(stderr,"Found: step=%d, nre=%d, ndisre=%d, nuser=%d time=%g.\n"
	    "Trying to skip frame expect a crash though\n",
	    eh.step,eh.nre,eh.ndisre,eh.nuser,eh.t);
    ndisre = 0;
    nuser  = 0;
  }
  for(i=0; (i<*nre); i++) {
    bOK = bOK && do_real(ener[i].e);
    
    tmp1=ener[i].eav;
    if((tmp1/(*step+1))<GMX_REAL_EPS)
      tmp1=0;
    bOK = bOK && do_real(tmp1);
    ener[i].eav=tmp1;
    
    /* This is to save only in single precision (unless compiled in DP) */
    tmp2 = ener[i].esum;
    bOK = bOK && do_real(tmp2);
    ener[i].esum = tmp2;
    
    bOK = bOK && do_real(ener[i].e2sum);
  }
  if (ndisre) {
    if (drblock->ndr == 0) {
      snew(drblock->rav,eh.ndisre);
      snew(drblock->rt,eh.ndisre);
      drblock->ndr = eh.ndisre;
    }
    ndo_real(drblock->rav,drblock->ndr,bOK1);
    bOK = bOK && bOK1;
    ndo_real(drblock->rt,drblock->ndr,bOK1);
    bOK = bOK && bOK1;
    if (bRead)
      *ndr = eh.ndisre;
  }
  else
    if (bRead)
      *ndr = 0;
  if (nuser) {
    fatal_error(0,"Can't handle user blocks");
  }
  if (!bOK) {
    if (bRead) {
      fprintf(stderr,"\nLast frame read %d                               ",framenr-1);
      fprintf(stderr,"\nWARNING: Incomplete frame: nr %6d time %8.3f     \n",
	      framenr,*t);
    } else 
      fatal_error(-1,"could not write energies");
    return FALSE; 
  }
  
  return TRUE;
}



