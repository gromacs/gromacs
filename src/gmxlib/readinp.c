/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * Gromacs Runs On Most of All Computer Systems
 */
static char *SRCID_readinp_c = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include "typedefs.h"
#include "string2.h"
#include "futil.h"
#include "smalloc.h"
#include "readinp.h"
#include "macros.h"

static int inp_count = 1;

t_inpfile *read_inpfile(char *fn,int *ninp)
{
  FILE      *in;
  char      buf[STRLEN],lbuf[STRLEN],rbuf[STRLEN];
  char      *ptr,*isptr,*cptr;
  t_inpfile *inp=NULL;
  int       nin,lc,i,j,k;

  inp_count = 1;
  in        = ffopen(fn,"r");
  nin = lc  = 0;
  do {
    ptr=fgets2(buf,STRLEN-1,in);
    lc++;
    if (ptr) {
      /* Strip comment */
      if ((cptr=strchr(buf,COMMENTSIGN)) != NULL)
	*cptr='\0';
      /* Strip spaces */
      trim(buf);
      
      for(j=0; (buf[j] != '=') && (buf[j] != '\0'); j++)
	;
      if (buf[j] == '\0') {
	if (j > 0)
	  fprintf(stderr,"No = on line %d in file %s, ignored\n",lc,fn);
      }
      else {
	for(i=0; (i<j); i++)
	  lbuf[i]=buf[i];
	lbuf[i]='\0';
	trim(lbuf);
	if (lbuf[0] == '\0')
	  fprintf(stderr,"Empty left hand side on line %d in file %s, ignored\n",lc,fn);
	else {
	  for(i=j+1,k=0; (buf[i] != '\0'); i++,k++)
	    rbuf[k]=buf[i];
	  rbuf[k]='\0';
	  trim(rbuf);
	  if (rbuf[0] == '\0')
	    fprintf(stderr,"Empty right hand side on line %d in file %s, ignored\n",lc,fn);
	  else {
	    /* Now finally something sensible */
	    srenew(inp,++nin);
	    inp[nin-1].count = 0;
	    inp[nin-1].bSet  = FALSE;
	    inp[nin-1].name  = strdup(lbuf);
	    inp[nin-1].value = strdup(rbuf);
	  }
	}
      }
    }
  } while (ptr);
  fclose(in);
  
  *ninp=nin;
  return inp;
}

static int inp_comp(const void *a,const void *b)
{
  return ((t_inpfile *)a)->count - ((t_inpfile *)b)->count;
}

static void sort_inp(int ninp,t_inpfile inp[])
{
  int i,mm;
  
  mm=-1;
  for(i=0; (i<ninp); i++)
    mm=max(mm,inp[i].count);
  for(i=0; (i<ninp); i++) {
    if (inp[i].count == 0)
      inp[i].count = mm++;
  }
  qsort(inp,ninp,sizeof(inp[0]),inp_comp);
}

void write_inpfile(char *fn,int ninp,t_inpfile inp[])
{
  FILE *out;
  int  i;

  sort_inp(ninp,inp);  
  out=ffopen(fn,"w");
  for(i=0; (i<ninp); i++) {
    if (inp[i].bSet)
      fprintf(out,"%-24s = %s\n",inp[i].name,inp[i].value ? inp[i].value : "");
    else
      fprintf(stderr,"Warning: unknown left-hand %s in parameter file\n",
	      inp[i].name);
  }
  fclose(out);
}

static int get_einp(int *ninp,t_inpfile **inp,char *name)
{
  static int  ival;
  static real rval;
  int    i;
  
  for(i=0; (i<(*ninp)); i++)
    if (gmx_strcasecmp(name,(*inp)[i].name) == 0)
      break;
  if (i == (*ninp)) {
    (*ninp)++;
    srenew(*inp,(*ninp));
    (*inp)[i].name=strdup(name);
    (*inp)[i].bSet=TRUE;
  }
  (*inp)[i].count = inp_count++;
    (*inp)[i].bSet  = TRUE;
  if (debug) 
    fprintf(debug,"Inp %d = %s\n",(*inp)[i].count,(*inp)[i].name);

    if (i == (*ninp)-1)
      return -1;
    else
      return i;
}

int get_eint(int *ninp,t_inpfile **inp,char *name,int def)
{
  char buf[32];
  int  ii;
  
  ii=get_einp(ninp,inp,name);
  
  if (ii == -1) {
    sprintf(buf,"%d",def);
    (*inp)[(*ninp)-1].value=strdup(buf);
    
    return def;
  }
  else 
    return atoi((*inp)[ii].value);
}

real get_ereal(int *ninp,t_inpfile **inp,char *name,real def)
{
  char buf[32];
  int  ii;
  
  ii=get_einp(ninp,inp,name);
  
  if (ii == -1) {
    sprintf(buf,"%g",def);
    (*inp)[(*ninp)-1].value=strdup(buf);
    
    return def;
  }
  else 
    return atof((*inp)[ii].value);
}

char *get_estr(int *ninp,t_inpfile **inp,char *name,char *def)
{
  char buf[32];
  int  ii;
  
  ii=get_einp(ninp,inp,name);
  
  if (ii == -1) {
    if (def) {
      sprintf(buf,"%s",def);
      (*inp)[(*ninp)-1].value=strdup(buf);
    }
    else
      (*inp)[(*ninp)-1].value=NULL;
      
    return def;
  }
  else 
    return (*inp)[ii].value;
}

int get_eenum(int *ninp,t_inpfile **inp,char *name,char **defs)
{
  int  ii,i,j;
  
  ii=get_einp(ninp,inp,name);
  
  if (ii == -1) {
    (*inp)[(*ninp)-1].value=strdup(defs[0]);
  
    return 0;
  }

  for(i=0; (defs[i] != NULL); i++)
    if (gmx_strcasecmp(defs[i],(*inp)[ii].value) == 0)
      break;
      
  if (defs[i] == NULL) {
    (*inp)[ii].value=strdup(defs[0]);
    fprintf(stderr,"Invalid enum %s for variable %s, using %s\n",
	    (*inp)[ii].value,name,defs[0]);
    fprintf(stderr,"Next time use one of:");
    j=0;
    while (defs[j]) {
      fprintf(stderr," '%s'",defs[j]);
      j++;
    }
    fprintf(stderr,"\n");
    
    return 0;
  }

  return i;
}
