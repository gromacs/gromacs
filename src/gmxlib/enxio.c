#include <typedefs.h>
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

static void edr_nms(XDR *xdr,int *nre,char ***nms)
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

static bool do_eheader(int fp,t_eheader *eh)
{
  bool bRead = fio_getread(fp);
  
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

void do_enxnms(int fp,int *nre,char ***nms)
{
  int  i;
  char buf[64];
  bool bRead;
  
  bRead = fio_getread(fp);
  if (fio_getftp(fp) == efEDR) {
    fio_select(fp);
    edr_nms(fio_getxdr(fp),nre,nms);
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

int open_enx(char *fn,char *mode)
{
  int       i,fp;
  char      **nm=NULL;
  int       nre;
  t_eheader *eh;
  
  if (strcmp(mode,"r") == 0) {
    fp=fio_open(fn,mode);
    fio_select(fp);
    fio_setprecision(fp,FALSE);
    do_enxnms(fp,&nre,&nm);
    snew(eh,1);
    do_eheader(fp,eh);
    
    /* Now check whether this file is in single precision */
    if (((eh->e_size && (eh->nre == nre) && (nre*4*sizeof(float) == eh->e_size)) ||
	 (eh->d_size && (eh->ndisre*sizeof(float)*2+sizeof(int) == eh->d_size)))){
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
      do_eheader(fp,eh);
      if (((eh->e_size && (eh->nre == nre) && (nre*4*sizeof(double) == eh->e_size)) ||
	   (eh->d_size && (eh->ndisre*sizeof(double)*2+sizeof(int) == eh->d_size))))
	fprintf(stderr,"Opened %s as double precision energy file\n",fn);
      else
	fatal_error(0,"Energy file %s not recognized, maybe different CPU?",fn);
      for(i=0; (i<nre); i++)
	sfree(nm[i]);
      sfree(nm);
    }
    sfree(eh);
    fio_rewind(fp);
  }
  else 
    fp = fio_open(fn,mode);
    
  return fp;
}

bool do_enx(int fp,real *t,int *step,int *nre,t_energy ener[],
	    t_drblock *drblock)
{
  int       i,nre_dum;
  t_eheader eh;
  bool      bRead;

  bRead = fio_getread(fp);
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
  fio_select(fp);

  if (!do_eheader(fp,&eh))
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



