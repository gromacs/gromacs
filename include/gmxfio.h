#ifndef _gmxfio_h
#define _gmxfio_h

#include "typedefs.h"
#include "xdrf.h"

/* Enumerated for different items in files */
enum { eitemHEADER, eitemIR, eitemBOX, 
       eitemTOP, eitemX, eitemV, eitemF, eitemNR };
       
/* Enumerated for data types in files */
enum { eioREAL, eioINT,   eioNUCHAR, eioUSHORT, 
       eioRVEC, eioNRVEC, eioIVEC,  eioSTRING, eioNR };

/* Functions for reading and writing data */
typedef bool do_func(void *item,int nitem,int eio,
		     char *desc,char *srcfile,int line);
		     
/* Global variables defined in gmxfio.h */
extern do_func *do_read;
extern do_func *do_write;
extern char *itemstr[eitemNR];
extern char *comment_str[eitemNR];

/********************************************************
 * Open and Close 
 ********************************************************/

extern int fio_open(char *fn,char *mode);
/* Open a new file for reading or writing.
 * The file type will be deduced from the file name.
 * If fn is NULL, stdin / stdout will be used for Ascii I/O (TPA type)
 * mode may be "r", "w", "a"
 */
 
extern void fio_close(int fp);
/* Close the file corresponding to fp (if not stdio)
 * The routine will exit when an invalid fio is handled.
 */

extern void fio_select(int fp);
/* This routine sets the global variables do_read and do_write
 * to point to the correct routines for fp.
 */

/********************************************************
 * Change properties of the open file
 ********************************************************/

extern void fio_setprecision(int fio,bool bDouble);
/* Select the floating point precision for reading and writing files */

extern char *fio_getname(int fio);
/* Return the filename corresponding to the fio index */

extern int fio_getftp(int fio);
/* Return the filetype corresponding to the fio index */

extern void fio_setftp_fio(int fio,int ftp);
/* And set it */

extern void fio_setdebug(int fio,bool bDebug);
/* Set the debug mode */

extern bool fio_getdebug(int fio);
/* Return  whether debug mode is on in fio  */

extern bool fio_getread(int fio);
/* Return  whether read mode is on in fio  */

/***************************************************
 * FILE Operations
 ***************************************************/

extern void fio_rewind(int fio);
/* Rewind the tpa file in fio */

extern void fio_flush(int fio);
/* Flush the fio */

extern long fio_ftell(int fio);
/* Return file position if possible */

extern void fio_seek(int fio,long fpos);
/* Set file position if possible, quit otherwise */

extern FILE *fio_getfp(int fio);
/* Return the file pointer itself */

extern XDR *fio_getxdr(int fio);
/* Return the file pointer itself */

/********************************************************
 * Dirty C macros... Try this in FORTRAN 
 * (Oh, and you can do structured programming in C too) 
 *********************************************************/
#define do_real(item)         (bRead ?\
  do_read ((void *)&(item),1,eioREAL,(#item),__FILE__,__LINE__) : \
  do_write((void *)&(item),1,eioREAL,(#item),__FILE__,__LINE__))
  
#define do_int(item)          (bRead ?\
  do_read ((void *)&(item),1,eioINT,(#item),__FILE__,__LINE__) :\
  do_write((void *)&(item),1,eioINT,(#item),__FILE__,__LINE__))
  
#define do_nuchar(item,n)     (bRead ?\
  do_read ((void *)(item),n,eioNUCHAR,(#item),__FILE__,__LINE__) :\
  do_write((void *)(item),n,eioNUCHAR,(#item),__FILE__,__LINE__))
  
#define do_ushort(item)          (bRead ?\
  do_read ((void *)&(item),1,eioUSHORT,(#item),__FILE__,__LINE__) :\
  do_write((void *)&(item),1,eioUSHORT,(#item),__FILE__,__LINE__))
  
#define do_rvec(item)         (bRead ?\
  do_read ((void *)(item),1,eioRVEC,(#item),__FILE__,__LINE__) :\
  do_write((void *)(item),1,eioRVEC,(#item),__FILE__,__LINE__))
  
#define do_ivec(item)         (bRead ?\
  do_read ((void *)(item),1,eioIVEC,(#item),__FILE__,__LINE__) :\
  do_write((void *)(item),1,eioIVEC,(#item),__FILE__,__LINE__))
  
#define do_string(item)       (bRead ?\
  do_read ((void *)(item),1,eioSTRING,(#item),__FILE__,__LINE__) :\
  do_write((void *)(item),1,eioSTRING,(#item),__FILE__,__LINE__))
  
#define ndo_real(item,n) \
  for(i=0; (i<n); i++) {\
    char buf[128];\
    sprintf(buf,"%s[%d]",#item,i);\
    (bRead ?\
     do_read ((void *)&((item)[i]),1,eioREAL,buf,__FILE__,__LINE__):\
     do_write((void *)&(item[i]),1,eioREAL,buf,__FILE__,__LINE__));\
  }
     
#define ndo_int(item,n)  \
  for(i=0; (i<n); i++) {\
    char buf[128];\
    sprintf(buf,"%s[%d]",#item,i);\
    (bRead ?\
     do_read ((void *)&(item[i]),1,eioINT,buf,__FILE__,__LINE__):\
     do_write((void *)&(item[i]),1,eioINT,buf,__FILE__,__LINE__));\
  }
  
#define ndo_rvec(item,n)      (bRead ?\
  do_read ((void *)(item),n,eioNRVEC,(#item),__FILE__,__LINE__) :\
  do_write((void *)(item),n,eioNRVEC,(#item),__FILE__,__LINE__))
  
#define ndo_ivec(item,n) \
  for(i=0; (i<n); i++) {\
    char buf[128];\
    sprintf(buf,"%s[%d]",#item,i);\
    (bRead ?\
     do_read ((void *)(item)[i],1,eioIVEC,buf,__FILE__,__LINE__):\
     do_write((void *)(item)[i],1,eioIVEC,buf,__FILE__,__LINE__));\
  }
  
#define ndo_string(item,n) \
  for(i=0; (i<n); i++) {\
    char buf[128];\
    sprintf(buf,"%s[%d]",#item,i);\
    (bRead ?\
     do_read ((void *)(item)[i],1,eioSTRING,buf,__FILE__,__LINE__):\
     do_write((void *)(item)[i],1,eioSTRING,buf,__FILE__,__LINE__));\
  }     

#endif
