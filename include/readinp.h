/*
 * $Id$
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"

typedef struct {
  int  count;
  bool bSet;
  char *name;
  char *value;
} t_inpfile;

extern t_inpfile *read_inpfile(char *fn,int *ninp);

extern void write_inpfile(char *fn,int ninp,t_inpfile inp[],
			  bool bHaltOnUnknown);

extern int get_eint(int *ninp,t_inpfile **inp,const char *name,int def);

extern real get_ereal(int *ninp,t_inpfile **inp,const char *name,real def);

extern char *get_estr(int *ninp,t_inpfile **inp,const char *name,char *def);

extern int get_eeenum(int *ninp,t_inpfile **inp,const char *name,const char **defs,
		      int *nerror,bool bPrintError);
/* defs must be NULL terminated, 
 * Add errors to nerror 
 * When bPrintError=TRUE and invalid enum: print "ERROR: ..."
 */

extern int get_eenum(int *ninp,t_inpfile **inp,const char *name,const char **defs);
/* defs must be NULL terminated */

/* Here are some macros to extract data from the inp structures.
 * Elements that are  removed  from the list after reading
 */
#define STYPE(name,var,def)  if ((tmp=get_estr(&ninp,&inp,name,def)) != NULL) strcpy(var,tmp)
#define ITYPE(name,var,def)  var=get_eint(&ninp,&inp,name,def)
#define RTYPE(name,var,def)  var=get_ereal(&ninp,&inp,name,def)
#define ETYPE(name,var,defs) var=get_eenum(&ninp,&inp,name,defs)
#define EETYPE(name,var,defs,nerr,bErr) var=get_eeenum(&ninp,&inp,name,defs,nerr,bErr)
#define CCTYPE(s) STYPE("\n; "s,dummy,NULL)
#define CTYPE(s)  STYPE("; "s,dummy,NULL)
/* This last one prints a comment line where you can add some explanation */

/* This structure is used for parsing arguments off the comand line */
enum { 
  etINT, etREAL, etTIME, etSTR,    etBOOL, etRVEC,   etENUM, etNR
};
/* names to print in help info */
static char *argtp[etNR] = {
  "int", "real", "time", "string", "bool", "vector", "enum" 
};

typedef struct {
  char *option;
  bool bSet;
  int  type;
  union {
    void *v;   /* This is a nasty workaround, to be able to use initialized */
    int  *i;   /* arrays */
    real *r;
    char **c;  /* Must be pointer to string (when type == etSTR)         */
               /* or null terminated list of enums (when type == etENUM) */
    bool *b;
    rvec *rv;
  } u;
  char *desc;
} t_pargs;

extern void get_pargs(int *argc,char *argv[],int nparg,t_pargs pa[],
		      bool bKeepArgs);
/* Read a number of arguments from the command line. 
 * For etINT, etREAL and etCHAR an extra argument is read (when present)
 * for etBOOL the boolean option is changed to the negate value
 * If !bKeepArgs, the command line arguments are removed from the command line
 */

extern bool is_hidden(t_pargs *pa);
/* Return TRUE when the option is a secret one */

extern char *pa_val(t_pargs *pa,char *buf, int sz);
/* Return the value of pa in the provided buffer buf, of size sz.
 * The return value is also a pointer to buf.
 */

extern int opt2parg_int(char *option,int nparg,t_pargs pa[]);

extern bool opt2parg_bool(char *option,int nparg,t_pargs pa[]);

extern real opt2parg_real(char *option,int nparg,t_pargs pa[]);

extern char *opt2parg_str(char *option,int nparg,t_pargs pa[]);

extern char *opt2parg_enum(char *option,int nparg,t_pargs pa[]);

extern bool opt2parg_bSet(char *option,int nparg,t_pargs pa[]);

extern void print_pargs(FILE *fp, int npargs,t_pargs pa[]);

extern void pr_enums(FILE *fp, int npargs,t_pargs pa[],int shell);

#endif
