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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _readinp_h
#define _readinp_h

#include "typedefs.h"
#include "warninp.h"


#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  int  count; /* sort order for output  */
  gmx_bool bObsolete; /* whether it is an obsolete param value */
  gmx_bool bSet; /* whether it it has been read out */
  char *name; /* name of the parameter */
  char *value; /* parameter value string */
  int inp_count; /* number of einps read. Only valid for the first item
                                          in the inpfile list. */
} t_inpfile;
/* entry in input files (like .mdp files). 
 Initally read in with read_inpfile, then filled in with missing values
 through get_eint, get_ereal, etc. */



t_inpfile *read_inpfile(const char *fn,int *ninp,
			       char **cppopts,
			       warninp_t wi);
/* Create & populate a t_inpfile struct from values in file fn. 
   fn = the file name
   ninp = the number of read parameters
   cppopts = the cpp-style options for #include paths and #defines */

void write_inpfile(const char *fn,int ninp,t_inpfile inp[],
			  gmx_bool bHaltOnUnknown,
			  warninp_t wi);

void replace_inp_entry(int ninp,t_inpfile *inp,
			      const char *old_entry,const char *new_entry);

int get_eint(int *ninp,t_inpfile **inp,const char *name,int def,
		      warninp_t wi);

gmx_large_int_t get_egmx_large_int(int *ninp,t_inpfile **inp,
					  const char *name,gmx_large_int_t def,
					  warninp_t);
  
double get_ereal(int *ninp,t_inpfile **inp,const char *name,double def,
			warninp_t wi);

const char *get_estr(int *ninp,t_inpfile **inp,const char *name,const char *def);

int get_eeenum(int *ninp,t_inpfile **inp,const char *name,const char **defs,
		      warninp_t wi);
/* defs must be NULL terminated */

int get_eenum(int *ninp,t_inpfile **inp,const char *name,const char **defs);
/* defs must be NULL terminated */

/* Here are some dirty macros to extract data from the inp structures.
 * Most macros assume the variables ninp, inp and wi are present.
 * Elements are removed from the list after reading.
 */
#define REM_TYPE(name)       replace_inp_entry(ninp,inp,name,NULL)
#define REPL_TYPE(old,new)   replace_inp_entry(ninp,inp,old,new)
#define STYPE(name,var,def)  if ((tmp=get_estr(&ninp,&inp,name,def)) != NULL) strcpy(var,tmp)
#define STYPENC(name,def) get_estr(&ninp,&inp,name,def)
#define ITYPE(name,var,def)  var=get_eint(&ninp,&inp,name,def,wi)
#define STEPTYPE(name,var,def)  var=get_egmx_large_int(&ninp,&inp,name,def,wi)
#define RTYPE(name,var,def)  var=get_ereal(&ninp,&inp,name,def,wi)
#define ETYPE(name,var,defs) var=get_eenum(&ninp,&inp,name,defs)
#define EETYPE(name,var,defs) var=get_eeenum(&ninp,&inp,name,defs,wi)
#define CCTYPE(s) STYPENC("\n; "s,NULL)
#define CTYPE(s)  STYPENC("; "s,NULL)
/* This last one prints a comment line where you can add some explanation */

/* This structure is used for parsing arguments off the comand line */
enum { 
  etINT, etGMX_LARGE_INT, etREAL, etTIME, etSTR,    etBOOL, etRVEC,   etENUM, etNR
};

/* names to print in help info */
static const char *argtp[etNR] = {
  "int", "step", "real", "time", "string", "bool", "vector", "enum" 
};

typedef struct {
  const char *option;
  gmx_bool bSet;
  int  type;
  union {
    void *v;   /* This is a nasty workaround, to be able to use initialized */
    int  *i;   /* arrays */
    gmx_large_int_t *is;
    real *r;
    const char **c; /* Must be pointer to string (when type == etSTR)         */
               /* or null terminated list of enums (when type == etENUM) */
    gmx_bool *b;
    rvec *rv;
  } u;
  const char *desc;
} t_pargs;

void get_pargs(int *argc,char *argv[],int nparg,t_pargs pa[],
		      gmx_bool bKeepArgs);
/* Read a number of arguments from the command line. 
 * For etINT, etREAL and etCHAR an extra argument is read (when present)
 * for etBOOL the gmx_boolean option is changed to the negate value
 * If !bKeepArgs, the command line arguments are removed from the command line
 */

gmx_bool is_hidden(t_pargs *pa);
/* Return TRUE when the option is a secret one */

char *pa_val(t_pargs *pa,char *buf, int sz);
/* Return the value of pa in the provided buffer buf, of size sz.
 * The return value is also a pointer to buf.
 */

int opt2parg_int(const char *option,int nparg,t_pargs pa[]);

gmx_bool opt2parg_gmx_bool(const char *option,int nparg,t_pargs pa[]);

real opt2parg_real(const char *option,int nparg,t_pargs pa[]);

const char *opt2parg_str(const char *option,int nparg,t_pargs pa[]);

const char *opt2parg_enum(const char *option,int nparg,t_pargs pa[]);

gmx_bool opt2parg_bSet(const char *option,int nparg,t_pargs pa[]);

void print_pargs(FILE *fp, int npargs,t_pargs pa[],gmx_bool bLeadingSpace);

char *pargs_print_line(t_pargs *pa,gmx_bool bLeadingSpace);

void pr_enums(FILE *fp, int npargs,t_pargs pa[],int shell);

#ifdef __cplusplus
}
#endif

#endif
