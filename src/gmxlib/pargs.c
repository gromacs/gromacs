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
 * Good ROcking Metal Altar for Chronical Sinners
 */
static char *SRCID_pargs_c = "$Id$";

#include <stdio.h>
#include <string.h>
#include "typedefs.h"
#include "fatal.h"
#include "statutil.h"
#include "readinp.h"
#include "smalloc.h"
#include "names.h"
#include "string2.h"

void get_pargs(int *argc,char *argv[],int nparg,t_pargs pa[],bool bKeepArgs)
{
  int  i,j;
  bool *bKeep;
  char buf[32];
  
  snew(bKeep,*argc+1);
  bKeep[0]     = TRUE;
  bKeep[*argc] = TRUE;
  
  for(i=1; (i<*argc); i++) {
    bKeep[i] = TRUE;
    for(j=0; (j<nparg); j++) {
      if (pa[j].type == etBOOL) {
	sprintf(buf,"-no%s",pa[j].option+1);
	if (strcmp(pa[j].option,argv[i])== 0) {
	  *pa[j].u.b = TRUE;
	  pa[j].bSet = TRUE;
	  bKeep[i] = FALSE;
	}
	else if (strcmp(buf,argv[i])== 0) {
	  *pa[j].u.b = FALSE;
	  pa[j].bSet = TRUE;
	  bKeep[i] = FALSE;
	}
      }
      else if (strcmp(pa[j].option,argv[i])== 0) {
	if (pa[j].bSet)
	  fprintf(stderr,"Setting option %s more than once!\n",pa[j].option);
	pa[j].bSet = TRUE;
	bKeep[i] = FALSE;
	switch(pa[j].type) {
	case etINT:
	  *pa[j].u.i = iscan(*argc,argv,&i);
	  break;
	case etREAL:
	  *pa[j].u.r = dscan(*argc,argv,&i);
	  break;
	case etSTR:
	  *(pa[j].u.c) = sscan(*argc,argv,&i);
	  break;
	default:
	  fatal_error(0,"Invalid type %d in pargs",pa[j].type);
	}
	/* i may be incremented, so set it to not keep */
	bKeep[i] = FALSE;
      }
    }
  }
  if (!bKeepArgs) {
    /* Remove used entries */
    for(i=j=0; (i<=*argc); i++) {
      if (bKeep[i])
	argv[j++]=argv[i];
    }
    (*argc)=j-1;
  }
  sfree(bKeep);
}

int opt2parg_int(char *option,int nparg,t_pargs pa[])
{
  int i;
  
  for(i=0; (i<nparg); i++)
    if (strcmp(pa[i].option,option) == 0)
      return *pa[i].u.i;
  
  fatal_error(0,"No integer option %s in pargs",option);
  
  return 0;
}

bool opt2parg_bool(char *option,int nparg,t_pargs pa[])
{
  int i;
  
  for(i=0; (i<nparg); i++)
    if (strcmp(pa[i].option,option) == 0)
      return *pa[i].u.b;
  
  fatal_error(0,"No boolean option %s in pargs",option);
  
  return FALSE;
}

real opt2parg_real(char *option,int nparg,t_pargs pa[])
{
  int i;
  
  for(i=0; (i<nparg); i++)
    if (strcmp(pa[i].option,option) == 0)
      return *pa[i].u.r;
  
  fatal_error(0,"No real option %s in pargs",option);
  
  return 0.0;
}

char *opt2parg_str(char *option,int nparg,t_pargs pa[])
{
  int i;
  
  for(i=0; (i<nparg); i++)
    if (strcmp(pa[i].option,option) == 0)
      return *(pa[i].u.c);
  
  fatal_error(0,"No string option %s in pargs",option);
  
  return NULL;
}

bool opt2parg_bSet(char *option,int nparg,t_pargs pa[])
{
  int i;
  
  for(i=0; (i<nparg); i++)
    if (strcmp(pa[i].option,option) == 0)
      return pa[i].bSet;
  
  fatal_error(0,"No such option %s in pargs",option);
  
  return FALSE; /* Too make some compilers happy */
}

char *pa_val(t_pargs *pa)
{
  static char buf[256];
  
  buf[0]='\0';
  switch(pa->type) {
  case etINT:
    sprintf(buf,"%d",*(pa->u.i));
    break;
  case etREAL:
    sprintf(buf,"%6g",*(pa->u.r));
    break;
  case etBOOL:
    sprintf(buf,"%6s",bool_names[*(pa->u.b)]);
    break;
  case etSTR:
    if (*(pa->u.c))
      if (strlen(*(pa->u.c)) >= 256)
	fatal_error(0,"Argument to long: \"%d\"\n",*(pa->u.c));
      else
	strcpy(buf,*(pa->u.c));
    break;
  }
  return buf;
}

void print_pargs(FILE *fp, int npargs,t_pargs pa[])
{
  static char *argtp[etNR] = { "int", "real", "string", "bool" };
  bool bShowHidden;
  char buf[32],buf2[256];
  char *ptr,*desc,*wdesc;
  int  i,j;
  
  /* Cannot call opt2parg_bSet here, because it crashes when the option
   * is not in the list (mdrun)
   */
  bShowHidden = FALSE;
  for(i=0; (i<npargs); i++) 
    if ((strcmp(pa[i].option,"-hidden")==0) && (pa[i].bSet))
      bShowHidden = TRUE;
  
  if (npargs > 0) {
    fprintf(fp,"%12s %6s %6s  %s\n","Opt","Type","Value","Description");
    fprintf(fp,"------------------------------------------------------\n");
    for(i=0; (i<npargs); i++) {
      ptr   = strstr(pa[i].desc,"HIDDEN");
      if (bShowHidden && ptr) {
	snew(desc,strlen(ptr)+4);
	sprintf(desc,"[hidden] %s",ptr+6);
      }
      else
	desc=strdup(pa[i].desc);
	
      if (bShowHidden || (ptr == NULL)) {
	if (pa[i].type == etBOOL)
	  sprintf(buf,"-[no]%s",pa[i].option+1);
	else
	  strcpy(buf,pa[i].option);
	sprintf(buf2,"%12s %6s %6s  %s\n",buf,argtp[pa[i].type],
		pa_val(&(pa[i])),desc);
	wdesc=wrap_lines(buf2,80,28);
	fprintf(fp,wdesc);
      }
      sfree(desc);
    }
    fprintf(fp,"\n");
  }
}
