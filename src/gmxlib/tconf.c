/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * Grunge ROck MAChoS
 */

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <sysstuff.h>
#include <string2.h>
#include <txtdump.h>
#include <rwtop.h>
#include <smalloc.h>
#include <typedefs.h>
#include <fatal.h>
#include <tconf.h>

void set_configuration(t_config *cf)
{
  cf->data.ir=(cf->data.header.ir_size)?cf->config.ir:NULL;
  cf->data.e=(cf->data.header.e_size)?cf->config.e:NULL;
  cf->data.box=(cf->data.header.box_size)?cf->config.box:NULL;
  cf->data.vir=(cf->data.header.vir_size)?cf->config.vir:NULL;
  cf->data.pres=(cf->data.header.pres_size)?cf->config.pres:NULL;
  cf->data.top=(cf->data.header.top_size)?cf->config.top:NULL;
  cf->data.x=(cf->data.header.x_size)?cf->config.x:NULL;
  cf->data.v=(cf->data.header.v_size)?cf->config.v:NULL;
  cf->data.f=(cf->data.header.f_size)?cf->config.f:NULL;
}

char *rd_configuration(FILE *fp,t_config *cf)
{
  char *version;

  version=rd_header(fp,&cf->data.header);
  set_configuration(cf);
  version=rd_hstatus(fp,&cf->data.header,&cf->data.step,&cf->data.t,
                     &cf->data.lambda,cf->data.ir,*cf->data.box,
		     *cf->data.vir,*cf->data.pres,
                     &cf->data.natoms,cf->data.x,cf->data.v,cf->data.f,
                     &cf->data.nre,cf->data.e,cf->data.top);
  return version;
}

void wr_configuration(FILE *fp,t_config *cf)
{
  wr_status(fp,cf->data.step,cf->data.t,cf->data.lambda,cf->data.ir,
            *cf->data.box,*cf->data.vir,*cf->data.pres,
	    cf->data.natoms,cf->data.x,cf->data.v,cf->data.f,
            cf->data.nre,cf->data.e,cf->data.top);
}

char *init_configuration(FILE *fp,t_config *cf)
{
  char *version;
  
  memset((char *)&cf->config,0,sizeof(cf->config));
  version=rd_header(fp,&cf->config.header);
  rewind(fp);
  snew(cf->config.ir,1);
  snew(cf->config.box,1);
  snew(cf->config.vir,1);
  snew(cf->config.pres,1);
  snew(cf->config.x,cf->config.header.natoms);
  snew(cf->config.v,cf->config.header.natoms);
  snew(cf->config.f,cf->config.header.natoms);
  snew(cf->config.e,F_NRE);
  snew(cf->config.top,1); 
  return version;
}

char *read_configuration(char *fn,t_config *cf)
{
  char *version;
  FILE *fp;
  
  if ((fp=fopen(fn,"r"))==NULL) fatal_error(errno,"reading %s",fn);
  version=init_configuration(fp,cf);
  version=rd_configuration(fp,cf);
  fclose(fp);
  return version;
}

void pr_configuration(FILE *fp,int indent,char *title,t_status_parm *sp)
{
  if (available(fp,sp,title))
    {
      indent=pr_title(fp,indent,title);
      pr_header(fp,indent,"header",&(sp->header));
      pr_inputrec(fp,indent,"ir",sp->ir);
      pr_rvecs(fp,indent,"box",*(sp->box),DIM);
      pr_rvecs(fp,indent,"vir",*(sp->vir),DIM);
      pr_rvecs(fp,indent,"pres",*(sp->pres),DIM);
      pr_rvecs(fp,indent,"x",sp->x,sp->natoms);
      pr_rvecs(fp,indent,"v",sp->v,sp->natoms);
      pr_rvecs(fp,indent,"f",sp->f,sp->natoms);
      pr_energies(fp,indent,"e",sp->e,sp->nre);
      pr_top(fp,indent,"topology",sp->top);
    }
}

