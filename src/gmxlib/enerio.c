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
static char *SRCID_enerio_c = "$Id$";

#include "errno.h"
#include "sysstuff.h"
#include "smalloc.h"
#include "string2.h"
#include "futil.h"
#include "typedefs.h"
#include "names.h"
#include "macros.h"
#include "enerio.h"
#include "fatal.h"
#include "xdrf.h"
#include "xtcio.h"
<<<<<<< enerio.c
#include "tpaio.h"
=======

void wr_ener_nms(FILE *out,int nre,char *nms[])
{
  int i;
  
  fprintf(out,"%d\n",nre);
  for(i=0; (i<nre); i++)
    fprintf(out,"%s\n",nms[i]);
}
>>>>>>> 1.2

static void low_wr_ener(FILE *out,t_eheader *eh,
			t_energy ener[],t_drblock *drblock,void *ublock)
{
  int size;
  
  fwrite(eh,sizeof(*eh),1,out);
  if (eh->e_size > 0) 
    fwrite(ener,1,eh->e_size,out);
  if (eh->d_size > 0) {
    size=drblock->ndr*sizeof(drblock->rav[0]);
    fwrite(drblock->rav,1,size,out);
    size=drblock->ndr*sizeof(drblock->rt[0]);
    fwrite(drblock->rt, 1,size,out);
  }
  if (eh->u_size > 0) 
    fwrite(ublock,1,eh->u_size,out);
}

void wr_ener(FILE *out,real t,int step,int nre,t_energy ener[],
	     t_drblock *drblock)
{
  t_eheader eh;
  
  eh.t=t;
  eh.step=step;
  eh.nre=nre;
  eh.ndisre=drblock ? (drblock->ndr) : 0;
  eh.nuser=0;
  eh.e_size=nre*sizeof(ener[0]);
  eh.d_size=eh.ndisre ? (drblock->ndr*
			 (sizeof(drblock->rav[0])+sizeof(drblock->rt[0]))) : 0;
  eh.u_size=0;
  low_wr_ener(out,&eh,ener,drblock,NULL);
}

void wr_ener_nms(FILE *out,int nre,char *nms[])
{
  int i;
  
  fprintf(out,"%d\n",nre);
  for(i=0; (i<nre); i++)
    fprintf(out,"%s\n",nms[i]);
}

void rd_ener_nms(FILE *in,int *nre,char ***nm)
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

bool gmx_fread(void *ptr,char *s,size_t size,size_t n,
	       FILE *fp,char *fn,int line)
{
  size_t  nfr;
  
  nfr=fread(ptr,size,n,fp);
  if (nfr != n) {
    fprintf(stderr,
	    "Error reading %s in file %s, line %d (only %d items iso %d)\n",
	    s,fn,line,nfr,n);
    perror(s);
    return FALSE;
  }
  return TRUE;
}

#define FREAD(p,s,n,fp) gmx_fread(p,#p,s,n,fp,__FILE__,__LINE__)

static ulong lowest_read(FILE *in,long fpos,int size,void *ptr,char *nm,
			 bool bCanSeek)
{
  if (ptr) {
    if (size > 0) {
      if (!FREAD(ptr,size,1,in))
	return 0;
    }
    else
      fatal_error(0,"%s not found in efile",nm);
  }
  else if (size > 0) {
    if (!bCanSeek)
      fatal_error(0,"Can not seek on this file (compressed ?)");
    fseek(in,fpos+size,SEEK_CUR);
  }
  return size;
}

static bool low_rd_ener(FILE *in,t_eheader *eh,
			t_energy ener[],t_drblock *drblock,void *ublock)
{
  bool  bCont,bCanSeek;
  int   size;
  long  fpos;

  if (!FREAD(eh,sizeof(*eh),1,in))
    return FALSE;

  bCanSeek = !is_pipe(in);
  
  if (bCanSeek)
    fpos=ftell(in);
  else
    fpos=0;
    
#define LREAD(size,ptr) lowest_read(in,fpos,size,ptr,#ptr,bCanSeek)
  if (eh->e_size > 0)
    fpos+=LREAD(eh->e_size,ener ? (void *)ener : NULL);
    
  if (eh->d_size > 0) {
    /* Only read distance restraint data when requested (drblock != NULL)
     * and when present (eh->d_size > 0).
     */
    if (drblock) 
      /* Allocate memory if necessary */
      if (drblock->ndr==0) {
	snew(drblock->rav,eh->ndisre);
	snew(drblock->rt,eh->ndisre);
	drblock->ndr=eh->ndisre;
      }
    size=eh->ndisre*sizeof(drblock->rav[0]);
    fpos+=LREAD(size,drblock ? (void *)drblock->rav : NULL);
    
    size=eh->ndisre*sizeof(drblock->rt[0]);
    fpos+=LREAD(size,drblock ? (void *)drblock->rt : NULL);
    
  }
    
  if (eh->u_size > 0)
    fpos+=LREAD(eh->u_size,ublock ? (void *)ublock : NULL);
#undef LREAD
  
  
  return TRUE;
}

bool rd_ener(FILE *in,real *t,int *step,t_energy ener[],t_drblock *drblock)
{
  bool      bCont;
  t_eheader eh;
  
  bCont = low_rd_ener(in,&eh,ener,drblock,NULL);
  if (bCont) {
    *t    = eh.t;
    *step = eh.step;
  }
  
  return bCont;
}

/************************************************
 *
 *   X D R format energies
 *
 ************************************************/
 
static bool low_edr_io(XDR *xdr,t_eheader *eh,
		       t_energy ener[],t_drblock *drblock,void *ublock)
{
  int i;
  
  if (xdr_real(xdr,&eh->t) == 0)
    return FALSE;
    
  XTC_CHECK_FAT_ERR("eh->step",xdr_int(xdr,&eh->step));
  XTC_CHECK_FAT_ERR("eh->nre",xdr_int(xdr,&eh->nre));
  XTC_CHECK_FAT_ERR("eh->ndisre",xdr_int(xdr,&eh->ndisre));
  XTC_CHECK_FAT_ERR("eh->nuser",xdr_int(xdr,&eh->nuser));
  XTC_CHECK_FAT_ERR("eh->e_size",xdr_int(xdr,&eh->e_size));
  XTC_CHECK_FAT_ERR("eh->d_size",xdr_int(xdr,&eh->d_size));
  XTC_CHECK_FAT_ERR("eh->u_size",xdr_int(xdr,&eh->u_size));
  
  if (eh->e_size > 0) {
    for(i=0; (i<eh->nre); i++) {
      XTC_CHECK_FAT_ERR("ener.e",    xdr_real(xdr,&ener[i].e));
      XTC_CHECK_FAT_ERR("ener.eav",  xdr_real(xdr,&ener[i].eav));
      XTC_CHECK_FAT_ERR("ener.esum", xdr_real(xdr,&ener[i].esum));
      XTC_CHECK_FAT_ERR("ener.e2sum",xdr_real(xdr,&ener[i].e2sum));
    }
  }
  if (eh->d_size > 0) {
    if (drblock->ndr == 0) {
      snew(drblock->rav,eh->ndisre);
      snew(drblock->rt,eh->ndisre);
      drblock->ndr=eh->ndisre;
    }
    for(i=0; (i<eh->ndisre); i++) 
      XTC_CHECK_FAT_ERR("drblock.rav",xdr_real(xdr,&(drblock->rav[i])));
    for(i=0; (i<eh->ndisre); i++) 
      XTC_CHECK_FAT_ERR("drblock.rt",xdr_real(xdr,&(drblock->rt[i])));
  }
  if (eh->u_size > 0) {
    fatal_error(0,"Can't handle user-blocks");
  }
  
  return TRUE;
}

bool edr_io(XDR *xdr,real *t,int *step,int *nre,t_energy ener[],
	    t_drblock *drblock)
{
  t_eheader eh;
  bool bSucces;
  
  eh.t=*t;
  eh.step=*step;
  eh.nre=*nre;
  eh.nuser=0;
  eh.e_size=*nre*sizeof(ener[0]);
  if (drblock) {
    eh.ndisre=drblock->ndr;
    if (eh.ndisre > 0)
      eh.d_size=eh.ndisre*(sizeof(drblock->rav[0]+drblock->rt[0]));
    else
      eh.d_size=0;
  }
  eh.u_size=0;

  bSucces=low_edr_io(xdr,&eh,ener,drblock,NULL);
  
  *t=eh.t;
  *step=eh.step;
  *nre=eh.nre;
  
  return bSucces;
}

void edr_nms(XDR *xdr,int *nre,char ***nms)
{
  int i;
  char **NM;
  
  NM=*nms;
  
  if (!xdr_int(xdr,nre)) {
    *nre=0;
    return;
  }
  if (NM == NULL) {
    snew(NM,*nre);
  }
  for(i=0; (i<*nre); i++)
    xdr_string(xdr,&(NM[i]),STRLEN);
  *nms=NM;
}

void do_enernms(int fp,bool bRead,int *nre,char ***nms)
{
  int i;
  char buf[64];

  if (fio_ftp(fp) == efEDR) {
    select_tpa(fp);
    edr_nms(xdr_tpa(fp),nre,nms);
  }
  else if (bRead)
    rd_ener_nms(fp_tpa(fp),nre,nms);
  else
    wr_ener_nms(fp_tpa(fp),*nre,*nms);
}

static bool do_eheader(int fp,bool bRead,t_eheader *eh)
{
  if (!do_real(eh->t))      return FALSE;
  if (!do_int (eh->step))   return FALSE;
  if (!do_int (eh->nre))    return FALSE;
  if (!do_int (eh->ndisre)) return FALSE;
  if (!do_int (eh->nuser))  return FALSE;
  if (!do_int (eh->e_size)) return FALSE;
  if (!do_int (eh->d_size)) return FALSE;
  if (!do_int (eh->u_size)) return FALSE;
  
  return TRUE;
}

void bswap4(void *n)
{
  ushort lo,hi;
  int    *nn;
  
  nn  = (int *)n;
  hi  = (*nn) & 0xff;
  lo  = (*nn-hi) >> 16;
  
  *nn = (hi << 16) + lo;
}

int open_enx(char *fn,char *mode)
{
  int       i,fp;
  char      **nm=NULL;
  int       nre;
  t_eheader *eh;
  
  if (strcmp(mode,"r") == 0) {
    fp=open_tpa(fn,mode);
    select_tpa(fp);
    set_precision(fp,FALSE);
    do_enernms(fp,TRUE,&nre,&nm);
    snew(eh,1);
    do_eheader(fp,TRUE,eh);
    
    /* Now check whether this file is in single precision */
    if (((eh->e_size && (eh->nre == nre) && (nre*4*sizeof(float) == eh->e_size)) ||
	 (eh->d_size && (eh->ndisre*sizeof(float)*2+sizeof(int) == eh->d_size)))){
      fprintf(stderr,"Opened %s as single precision energy file\n",fn);
      for(i=0; (i<nre); i++)
	sfree(nm[i]);
      sfree(nm);
    }
    else {
      rewind_tpa(fp);
      select_tpa(fp);
      set_precision(fp,TRUE);
      do_enernms(fp,TRUE,&nre,&nm);
      do_eheader(fp,TRUE,eh);
      if (((eh->e_size && (eh->nre == nre) && (nre*4*sizeof(float) == eh->e_size)) ||
	   (eh->d_size && (eh->ndisre*sizeof(float)*2+sizeof(int) == eh->d_size))))
	fprintf(stderr,"Opened %s as double precision energy file\n",fn);
      else
	fatal_error(0,"Energy file %s not recognized, maybe different CPU?",fn);
      for(i=0; (i<nre); i++)
	sfree(nm[i]);
      sfree(nm);
    }
    sfree(eh);
    rewind_tpa(fp);
  }
  else 
    fp = open_tpa(fn,mode);
    
  return fp;
}

bool do_energy(int fp,bool bRead,
	       real *t,int *step,int *nre,t_energy ener[],
	       t_drblock *drblock)
{
  int       i,nre_dum;
  t_eheader eh;

  if (!bRead) {  
    eh.t      = *t;
    eh.step   = *step;
    eh.nre    = *nre;
    eh.ndisre = drblock ? (drblock->ndr) : 0;
    eh.nuser  = 0;
    eh.e_size = (*nre)*sizeof(ener[0]);
    eh.d_size = eh.ndisre ? (drblock->ndr*
		(sizeof(drblock->rav[0])+sizeof(drblock->rt[0]))) : 0;
    eh.u_size = 0;
  }
  select_tpa(fp);

  if (!do_eheader(fp,bRead,&eh))
    return FALSE;
  *t    = eh.t;
  *step = eh.step;
  *nre  = eh.nre;
  for(i=0; (i<*nre); i++) {
    do_real(ener[i].e);
    do_real(ener[i].eav);
    do_real(ener[i].esum);
    do_real(ener[i].e2sum);
  }
  if (eh.ndisre) {
    if (drblock->ndr == 0) {
      snew(drblock->rav,eh.ndisre);
      snew(drblock->rt,eh.ndisre);
      drblock->ndr = eh.ndisre;
    }
    ndo_real(drblock->rav,drblock->ndr);
    ndo_real(drblock->rt,drblock->ndr);
  }
  if (eh.u_size) {
    fatal_error(0,"Can't handle user blocks");
  }
  return TRUE;
}



