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
 * Gyas ROwers Mature At Cryogenic Speed
 */
static char *SRCID_statusio_c = "$Id$";

 
#include <string.h>
#include "sysstuff.h"
#include "statusio.h"
#include "fatal.h"
#include "names.h"
#include "binio.h"
#include "recio.h"
#include "rwtop.h"
#include "futil.h"

#define offset(base,field) (((char *)&(base).field)-((char *)&(base)))
#define can_read(fp,ptr,len) _can_read((fp),(ptr),(len),#ptr,__FILE__,__LINE__)

static int _can_read(FILE *fp,void *buf,int len,char *what,char *file,int line)
{
  int ok;
  
  ok=0;
  if (buf==NULL)
    {
      if (len!=0)
        if (fseek(fp,len,SEEK_CUR)!=0)
          fatal_error(errno,"skip %s in %s, line %d",what,file,line);
    }
  else
    {
      if (len==0)
        fatal_error(0,"%s not available in %s, line %d",what,file,line);
      else
        {
          ok=1;
          /* (void) fprintf(stderr,"can read %s\n",what); */
        }
    }
  return ok;
}

size_t wr_status(FILE *fp,int step,real t,real lambda,
		 t_inputrec *ir,rvec *box,rvec *vir,rvec *pres,int natoms,
		 rvec *x,rvec *v,rvec *f,
		 int nre,t_energy *e,t_topology *top)
{
  size_t fpos,ir_pos;
  int patched;
  t_statheader *sh;

  snew(sh,1);
  sh->ir_size=0;
  sh->e_size=(e)?(nre*sizeof(e[0])):0;
  sh->box_size=(box)?sizeof(matrix):0;
  sh->vir_size=(vir)?sizeof(tensor):0;
  sh->pres_size=(pres)?sizeof(tensor):0;
  sh->top_size=0;
  sh->sym_size=0;
  sh->x_size=((x)?(natoms*sizeof(x[0])):0);
  sh->v_size=((v)?(natoms*sizeof(v[0])):0);
  sh->f_size=((f)?(natoms*sizeof(f[0])):0);
  sh->natoms=natoms;
  sh->step=step;
  sh->nre=nre;
  sh->t=t;
  sh->lambda=lambda;
  fpos=ftell(fp);
  wr_header(fp,sh);
  patched=0;
  if (ir) 
    {
      ir_pos=ftell(fp); 
      inputrec_blockio(fp,write,*ir);
      sh->ir_size=ftell(fp)-ir_pos; 
      patched=1;
    }
  if (sh->e_size!=0) cblockwrite(fp,e,sh->e_size);
  if (sh->box_size!=0) cblockwrite(fp,box,sh->box_size);
  if (sh->vir_size!=0) cblockwrite(fp,vir,sh->vir_size);
  if (sh->pres_size!=0) cblockwrite(fp,pres,sh->pres_size);
  
  if (top!=NULL) { sh->top_size=wr_top(fp,top); patched=1; }
  if (sh->x_size!=0) nblockwrite(fp,natoms,x);
  if (sh->v_size!=0) nblockwrite(fp,natoms,v);
  if (sh->f_size!=0) nblockwrite(fp,natoms,f);
    
  if (patched) patch(fp,fpos,wr_header(fp,sh));
  
  sfree(sh);
  return (ftell(fp)-fpos);
}

char *rd_hstatus(FILE *fp,t_statheader *sh,int *step,real *t,real *lambda,
                 t_inputrec *ir,rvec *box,rvec *vir,rvec *pres,
		 int *natoms,rvec *x,rvec *v,
                 rvec *f,int *nre,t_energy *e,t_topology *top)
{
  static char *version="@(#) statusio.c 1.57 9/30/97";

  if (can_read(fp,ir,sh->ir_size)) 
    inputrec_blockio(fp,read,*ir);
  if (can_read(fp,e,sh->e_size)) 
    cblockread(fp,e,sh->e_size);
  if (can_read(fp,box,sh->box_size)) 
    cblockread(fp,box,sh->box_size);
  if (can_read(fp,vir,sh->vir_size)) 
    cblockread(fp,vir,sh->vir_size);
  if (can_read(fp,pres,sh->pres_size)) 
    cblockread(fp,pres,sh->pres_size);
  if (can_read(fp,top,sh->top_size)) 
    (void) rd_top(fp,top);
  if (can_read(fp,x,sh->x_size)) 
    cblockread(fp,x,sh->x_size);
  if (can_read(fp,v,sh->v_size)) 
    cblockread(fp,v,sh->v_size);
  if (can_read(fp,f,sh->f_size)) 
    cblockread(fp,f,sh->f_size);
  *step=sh->step;
  *t=sh->t;
  *lambda=sh->lambda;
  *natoms=sh->natoms;
  *nre=sh->nre;
  return version;
}

char *rd_status(FILE *fp,int *step,real *t,real *lambda,
                t_inputrec *ir,rvec *box,rvec *vir,rvec *pres,
		int *natoms,rvec *x,
                rvec *v,rvec *f,int *nre,t_energy *e,
                t_topology *top)
{
  char *version;
  t_statheader sh;

  (void) rd_header(fp,&sh);
  version=rd_hstatus(fp,&sh,step,t,lambda,ir,box,vir,pres,natoms,x,v,f,nre,e,
                     top);
  return version;
}

void write_status(char *fn,int step,real t,real lambda,t_inputrec *ir,
                  rvec *box,rvec *vir,rvec *pres,
		  int natoms,rvec *x,rvec *v,rvec *f,
                  int nre,t_energy *e,t_topology *top)
{
  FILE *fp;
  
  fp=ffopen(fn,"wb");
  (void) wr_status(fp,step,t,lambda,ir,box,vir,pres,natoms,x,v,f,nre,e,
                   top);
  ffclose(fp);
}

char *read_status(char *fn,int *step,real *t,real *lambda,t_inputrec *ir,
                  rvec *box,rvec *vir,rvec *pres,
		  int *natoms,rvec *x,rvec *v,rvec *f,int *nre,
                  t_energy *e,t_topology *top)
{
  char *version;
  FILE *fp;
  
  fp=ffopen(fn,"rb");
  version=rd_status(fp,step,t,lambda,ir,box,vir,pres,natoms,x,v,f,nre,e,
                    top);
  ffclose(fp);
  return version;
}

void read_status_header(char *fn,t_statheader *header)
{
  FILE *fp;
  
  fp=ffopen(fn,"rb");
  (void) rd_header(fp,header);
  ffclose(fp);
}

