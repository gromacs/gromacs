/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * Great Red Owns Many ACres of Sand 
 */
static char *SRCID_enxio_c = "$Id$";
#include "futil.h"
#include "string2.h"
#include "fatal.h"
#include "smalloc.h"
#include "gmxfio.h"
#include "enxio.h"

void free_enxframe(t_enxframe *fr)
{
  int b;

  if (fr->e_alloc)
    sfree(fr->ener);
  if (fr->d_alloc) {
    sfree(fr->rav);
    sfree(fr->rt);
  }
  for(b=0; b<fr->nblock; b++)
    sfree(fr->block[b]);
  sfree(fr->block);
  sfree(fr->b_alloc);
  sfree(fr->nr);
}

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

static bool do_eheader(int fp,t_enxframe *fr,bool *bOK)
{
  int  block,i,dum=0;
  bool bRead = fio_getread(fp);
  int  tempfix_nr=0;

  *bOK=TRUE;
  if (!do_real(fr->t))       return FALSE;
  if (!do_int (fr->step))    *bOK = FALSE;
  if (!do_int (fr->nre))     *bOK = FALSE;
  if (!do_int (fr->ndisre))  *bOK = FALSE;
  if (!do_int (fr->nblock))  *bOK = FALSE;
  if (bRead && fr->nblock>70) {
    /* Temporary fix for intermediate file version, only used by B. Hess */
    tempfix_nr = fr->nblock;
    fr->nblock = 1;
  }
  if (bRead && fr->nblock>fr->nr_alloc) {
    srenew(fr->nr,fr->nblock);
    srenew(fr->b_alloc,fr->nblock);
    srenew(fr->block,fr->nblock);
    for(i=fr->nr_alloc; i<fr->nblock; i++) {
      fr->block[i]   = NULL;
      fr->b_alloc[i] = 0;
    }
    fr->nr_alloc = fr->nblock;
  }
  if (tempfix_nr)
    fr->nr[0] = tempfix_nr;
  else {
    for(block=0; block<fr->nblock; block++)
      if (!do_int (fr->nr[block])) 
	*bOK = FALSE;
  }
  if (!do_int (fr->e_size))  *bOK = FALSE;
  if (!do_int (fr->d_size))  *bOK = FALSE;
  /* Do a dummy int to keep the format compatible with the old code */
  if (!do_int (dum))         *bOK = FALSE;
  
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
  int        fp,nre,i;
  char       **nm=NULL;
  t_enxframe *fr;
  bool       bDum=TRUE;
  
  /* Energy files should always be opened as binary files,
   * but that is checked in fio_open.
   */

  if (mode[0]=='r') {
    fp=fio_open(fn,mode);
    fio_select(fp);
    fio_setprecision(fp,FALSE);
    do_enxnms(fp,&nre,&nm);
    snew(fr,1);
    do_eheader(fp,fr,&bDum);
    
    /* Now check whether this file is in single precision */
    if (((fr->e_size && (fr->nre == nre) && 
	  (nre*4*sizeof(float) == fr->e_size)) ||
	 (fr->d_size && 
	  (fr->ndisre*sizeof(float)*2+sizeof(int) == fr->d_size)))){
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
      do_eheader(fp,fr,&bDum);
      if (((fr->e_size && (fr->nre == nre) && 
	    (nre*4*sizeof(double) == fr->e_size)) ||
	   (fr->d_size && 
	    (fr->ndisre*sizeof(double)*2+sizeof(int) == fr->d_size))))
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
    free_enxframe(fr);
    sfree(fr);
    fio_rewind(fp);
  }
  else 
    fp = fio_open(fn,mode);
    
  framenr=0;
    
  return fp;
}

bool do_enx(int fp,t_enxframe *fr)
{
  int       i,block;
  bool      bRead,bOK,bOK1,bSane;
  real      tmp1,tmp2;

  bOK = TRUE;
  bRead = fio_getread(fp);
  if (!bRead) {  
    fr->e_size = fr->nre*sizeof(fr->ener[0].e)*4;
    fr->d_size = fr->ndisre*(sizeof(fr->rav[0]) + 
			   sizeof(fr->rt[0]));
  }
  fio_select(fp);

  if (!do_eheader(fp,fr,&bOK)) {
    if (bRead) {
	fprintf(stderr,"\rLast frame read %d                          ",framenr-1);
	if (!bOK)
	  fprintf(stderr,"\nWARNING: Incomplete frame: nr %6d time %8.3f\n",
		  framenr,fr->t);
    }
    return FALSE;
  }
  if (bRead) {
    if ( ( framenr<10 ) || ( framenr%10 == 0) ) 
      fprintf(stderr,"\rReading frame %6d time %8.3f           ",framenr,fr->t);
    framenr++;
  }
  /* Check sanity of this header */
  bSane = (fr->nre > 0 || fr->ndisre > 0);
  for(block=0; block<fr->nblock; block++)
    bSane = bSane || (fr->nr[block] > 0);
  if (!((fr->step >= 0) && bSane)) {
    fprintf(stderr,"\nWARNING: there may be something wrong with energy file %s\n",
	    fio_getname(fp));
    fprintf(stderr,"Found: step=%d, nre=%d, ndisre=%d, nblock=%d, time=%g.\n"
	    "Trying to skip frame expect a crash though\n",
	    fr->step,fr->nre,fr->ndisre,fr->nblock,fr->t);
  }
  if (bRead && fr->nre>fr->e_alloc) {
    srenew(fr->ener,fr->nre);
    fr->e_alloc = fr->nre;
  }
  for(i=0; i<fr->nre; i++) {
    bOK = bOK && do_real(fr->ener[i].e);
    
    tmp1 = fr->ener[i].eav;
    if((tmp1/(fr->step+1))<GMX_REAL_EPS)
      tmp1=0;
    bOK = bOK && do_real(tmp1);
    fr->ener[i].eav = tmp1;
    
    /* This is to save only in single precision (unless compiled in DP) */
    tmp2 = fr->ener[i].esum;
    bOK = bOK && do_real(tmp2);
    fr->ener[i].esum = tmp2;
    
    bOK = bOK && do_real(fr->ener[i].e2sum);
  }
  if (fr->ndisre) {
    if (bRead && fr->ndisre>fr->d_alloc) {
      srenew(fr->rav,fr->ndisre);
      srenew(fr->rt,fr->ndisre);
      fr->d_alloc = fr->ndisre;
    }
    ndo_real(fr->rav,fr->ndisre,bOK1);
    bOK = bOK && bOK1;
    ndo_real(fr->rt,fr->ndisre,bOK1);
    bOK = bOK && bOK1;
  }
  for(block=0; block<fr->nblock; block++) {
    if (bRead && fr->nr[block]>fr->b_alloc[block]) {
      srenew(fr->block[block],fr->nr[block]);
      fr->b_alloc[block] = fr->nr[block];
    }
    ndo_real(fr->block[block],fr->nr[block],bOK1);
    bOK = bOK && bOK1;
  }
  if (!bOK) {
    if (bRead) {
      fprintf(stderr,"\nLast frame read %d                               ",
	      framenr-1);
      fprintf(stderr,"\nWARNING: Incomplete frame: nr %6d time %8.3f     \n",
	      framenr,fr->t);
    } else 
      fatal_error(-1,"could not write energies");
    return FALSE; 
  }
  
  return TRUE;
}



