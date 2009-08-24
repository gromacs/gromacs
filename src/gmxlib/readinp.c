/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include "typedefs.h"
#include "string2.h"
#include "futil.h"
#include "smalloc.h"
#include "readinp.h"
#include "macros.h"
#include "statutil.h"
#include "gmxfio.h"
#include "gmxcpp.h"
#include "names.h"

static int inp_count = 1;

/* find an entry; return index, or -1 if not found */
static int search_einp(int ninp, const t_inpfile *inp, const char *name);


t_inpfile *read_inpfile(const char *fn,int *ninp, 
			char **cppopts)
{
  /*FILE      *in;*/
  gmx_cpp_t in;
  char      buf[STRLEN],lbuf[STRLEN],rbuf[STRLEN];
  char      *ptr,*cptr;
  t_inpfile *inp=NULL;
  int       nin,lc,i,j,k,status;
  /* setting cppopts from command-line options would be cooler */
  char *cppopts_null[] = { NULL }; 
  char **cppopts_given;
  bool allow_override=FALSE;

    
  if (debug)
    fprintf(debug,"Reading MDP file %s\n",fn);
  inp_count = 1;
  if (!cppopts)
    cppopts_given=cppopts;
  else
    cppopts_given=cppopts_null;

  status = cpp_open_file(fn, &in, cppopts_given);
    
  if (status != eCPP_OK)
  {
    sprintf(warn_buf, "%s", cpp_error(&in, status));
    warning_error(NULL);
    return NULL;
  }
  nin = lc  = 0;
  do {
    status=cpp_read_line(&in, STRLEN-1, buf);
    
    if (status !=eCPP_OK)
    {
      if (status != eCPP_EOF)
      {
	sprintf(warn_buf, "%s", cpp_error(&in, status));
	warning_error(NULL);
	return NULL;
      }
      else
	ptr=0;
    }
    else
    {
      ptr=buf;
    }
    
    
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
	if (j > 0) {
	  if (debug)
	    fprintf(debug,"No = on line %d in file %s, ignored\n",lc,fn);
	}
      }
      else {
	for(i=0; (i<j); i++)
	  lbuf[i]=buf[i];
	lbuf[i]='\0';
	trim(lbuf);
	if (lbuf[0] == '\0') {
	  if (debug) 
	    fprintf(debug,"Empty left hand side on line %d in file %s, ignored\n",lc,fn);
	}
	else {
	  for(i=j+1,k=0; (buf[i] != '\0'); i++,k++)
	    rbuf[k]=buf[i];
	  rbuf[k]='\0';
	  trim(rbuf);
	  if (rbuf[0] == '\0') {
	    if (debug)
	      fprintf(debug,"Empty right hand side on line %d in file %s, ignored\n",lc,fn);
	  }
	  else {
	    /* Now finally something sensible */
	    int found_index;
	    
	    /* first check whether we hit the 'multiple_entries' option */
	    if (strcasecmp_min(eMultentOpt_names[eMultentOptName], lbuf)==0) 
	    {
	      /* we now check whether to allow overrides from here or not */
	      if (strcasecmp_min(eMultentOpt_names[eMultentOptNo], rbuf)==0)
	      {
		allow_override=FALSE;
	      }
	      else if (strcasecmp_min(eMultentOpt_names[eMultentOptLast], rbuf)==0)
	      {
		allow_override=TRUE;
	      }
	      else
	      {
		sprintf(warn_buf,
			"Parameter \"%s\" should either be %s or %s\n",
			lbuf,
			eMultentOpt_names[eMultentOptNo],
			eMultentOpt_names[eMultentOptLast]);
		warning_error(NULL);
	      }
	    }
	    else
	    {
	      /* it is a regular option; check for duplicates */
	      found_index=search_einp(nin, inp, lbuf);
	      
	      if (found_index == -1)
	      {
		/* add a new item */
		srenew(inp,++nin);
		inp[nin-1].count      = 0;
		inp[nin-1].bObsolete  = FALSE;
		inp[nin-1].bSet       = FALSE;
		inp[nin-1].name       = strdup(lbuf);
		inp[nin-1].value      = strdup(rbuf);
	      }
	      else
	      {
		if (!allow_override)
		{
		  sprintf(warn_buf,
			  "Parameter \"%s\" doubly defined (and multiple assignments not allowed)\n",
			  lbuf);
		  warning_error(NULL);
		}
		else
		{
		  /* override */
		  sfree(inp[found_index].value);
		  inp[found_index].value = strdup(rbuf);
		  sprintf(warn_buf,
			  "Overriding existing parameter \"%s\" with value \"%s\"\n",
			  lbuf,rbuf);
		  warning_note(warn_buf);
		}
	      }
	    }
	  }
	}
      }
    }
  } while (ptr);
  cpp_close_file(&in);

  if (debug)
    fprintf(debug,"Done reading MDP file, there were %d entries in there\n",
	    nin);
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
  qsort(inp,ninp,(size_t)sizeof(inp[0]),inp_comp);
}

void write_inpfile(const char *fn,int ninp,t_inpfile inp[],bool bHaltOnUnknown)
{
  FILE *out;
  int  i;

  sort_inp(ninp,inp);  
  out=gmx_fio_fopen(fn,"w");
  nice_header(out,fn);
  for(i=0; (i<ninp); i++) {
    if (inp[i].bSet) {
      if(inp[i].name[0]==';' || (strlen(inp[i].name)>2 && inp[i].name[1]==';'))
	fprintf(out,"%-24s\n",inp[i].name);
      else
	fprintf(out,"%-24s = %s\n",inp[i].name,inp[i].value ? inp[i].value : "");
    } else if (!inp[i].bObsolete) {
      sprintf(warn_buf,"Unknown left-hand '%s' in parameter file\n",
	      inp[i].name);
      if (bHaltOnUnknown) {
	warning_error(NULL);
      } else {
	warning(NULL);
      }
    }
  }
  gmx_fio_fclose(out);

  check_warning_error(FARGS);
}

void replace_inp_entry(int ninp,t_inpfile *inp,const char *old_entry,const char *new_entry)
{
  int  i;
  
  for(i=0; (i<ninp); i++) {
    if (strcasecmp_min(old_entry,inp[i].name) == 0) {
      if (new_entry) {
	fprintf(stderr,"Replacing old mdp entry '%s' by '%s'\n",
		inp[i].name,new_entry);
	sfree(inp[i].name);
	inp[i].name = strdup(new_entry);
      } else {
	fprintf(stderr,"Ignoring obsolete mdp entry '%s'\n",
		inp[i].name);
	inp[i].bObsolete = TRUE;
      }
    }
  }
}

static int search_einp(int ninp, const t_inpfile *inp, const char *name)
{
  int i;
  
  if (inp==NULL)
    return -1;
  for(i=0; i<ninp; i++)
    if (strcasecmp_min(name,inp[i].name) == 0)
      return i;
  return -1;
}

static int get_einp(int *ninp,t_inpfile **inp,const char *name)
{
  int    i;
  int	 notfound=FALSE;

/*  if (inp==NULL)
    return -1;
  for(i=0; (i<(*ninp)); i++)
    if (strcasecmp_min(name,(*inp)[i].name) == 0)
      break;
  if (i == (*ninp)) {*/
  i=search_einp(*ninp, *inp, name);
  if (i == -1)
  {
    notfound=TRUE;
    i=(*ninp)++;
    srenew(*inp,(*ninp));
    (*inp)[i].name=strdup(name);
    (*inp)[i].bSet=TRUE;
  }
  (*inp)[i].count = inp_count++;
  (*inp)[i].bSet  = TRUE;
  if (debug) 
    fprintf(debug,"Inp %d = %s\n",(*inp)[i].count,(*inp)[i].name);
  
  /*if (i == (*ninp)-1)*/
  if (notfound)
    return -1;
  else
    return i;
}

int get_eint(int *ninp,t_inpfile **inp,const char *name,int def)
{
  char buf[32],*ptr;
  int  ii;
  int  ret;
  
  ii=get_einp(ninp,inp,name);
  
  if (ii == -1) {
    sprintf(buf,"%d",def);
    (*inp)[(*ninp)-1].value=strdup(buf);
    
    return def;
  }
  else {
    ret = strtol((*inp)[ii].value,&ptr,10);
    if (ptr == (*inp)[ii].value) {
      sprintf(warn_buf,"Right hand side '%s' for parameter '%s' in parameter file is not an integer value\n",(*inp)[ii].value,(*inp)[ii].name);
      warning_error(NULL);
    }

    return ret;
  }
}

gmx_step_t get_egmx_step_t(int *ninp,t_inpfile **inp,
			   const char *name,gmx_step_t def)
{
  char buf[32],*ptr;
  int  ii;
  gmx_step_t ret;
  
  ii=get_einp(ninp,inp,name);
  
  if (ii == -1) {
    sprintf(buf,gmx_step_pfmt,def);
    (*inp)[(*ninp)-1].value=strdup(buf);
    
    return def;
  }
  else {
    ret = str_to_gmx_step_t((*inp)[ii].value,&ptr);
    if (ptr == (*inp)[ii].value) {
      sprintf(warn_buf,"Right hand side '%s' for parameter '%s' in parameter file is not an integer value\n",(*inp)[ii].value,(*inp)[ii].name);
      warning_error(NULL);
    }

    return ret;
  }
}

double get_ereal(int *ninp,t_inpfile **inp,const char *name,double def)
{
  char buf[32],*ptr;
  int  ii;
  double ret;
  
  ii=get_einp(ninp,inp,name);
  
  if (ii == -1) {
    sprintf(buf,"%g",def);
    (*inp)[(*ninp)-1].value=strdup(buf);
    
    return def;
  }
  else {
    ret = strtod((*inp)[ii].value,&ptr);
    if (ptr == (*inp)[ii].value) {
      sprintf(warn_buf,"Right hand side '%s' for parameter '%s' in parameter file is not a real value\n",(*inp)[ii].value,(*inp)[ii].name);
      warning_error(NULL);
    }

    return ret;
  }
}

const char *get_estr(int *ninp,t_inpfile **inp,const char *name,const char *def)
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

int get_eeenum(int *ninp,t_inpfile **inp,const char *name,const char **defs,
	       int *nerror,bool bPrintError)
{
  int  ii,i,j;
  
  ii=get_einp(ninp,inp,name);
  
  if (ii == -1) {
    (*inp)[(*ninp)-1].value=strdup(defs[0]);
    
    return 0;
  }
  
  for(i=0; (defs[i] != NULL); i++)
    if (strcasecmp_min(defs[i],(*inp)[ii].value) == 0)
      break;
  
  if (defs[i] == NULL) {
    fprintf(stderr,"%snvalid enum '%s' for variable %s, using '%s'\n",
	    bPrintError ? "ERROR: i" : "I",(*inp)[ii].value,name,defs[0]);
    fprintf(stderr,"Next time use one of:");
    (*nerror)++;
    j=0;
    while (defs[j]) {
      fprintf(stderr," '%s'",defs[j]);
      j++;
    }
    fprintf(stderr,"\n");
    (*inp)[ii].value=strdup(defs[0]);
    
    return 0;
  }
  
  return i;
}

int get_eenum(int *ninp,t_inpfile **inp,const char *name,const char **defs)
{
  int dum=0;

  return get_eeenum(ninp,inp,name,defs,&dum,FALSE);
}

