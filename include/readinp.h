/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifndef _readinp_h
#define _readinp_h

static char *SRCID_readinp_h = "$Id$";

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

extern void write_inpfile(char *fn,int ninp,t_inpfile inp[]);

extern int get_eint(int *ninp,t_inpfile **inp,char *name,int def);

extern real get_ereal(int *ninp,t_inpfile **inp,char *name,real def);

extern char *get_estr(int *ninp,t_inpfile **inp,char *name,char *def);

extern int get_eeenum(int *ninp,t_inpfile **inp,char *name,char **defs,
		      int *nerror,bool bPrintError);
/* defs must be NULL terminated, 
 * Add errors to nerror 
 * When bPrintError=TRUE and invalid enum: print "ERROR: ..."
 */

extern int get_eenum(int *ninp,t_inpfile **inp,char *name,char **defs);
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
enum { etINT, etREAL, etSTR, etBOOL, etRVEC, etENUM, etNR };
/* names to print in help info */
static char *argtp[etNR] = { "int", "real", "string", "bool", "vector", "enum" };

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

extern char *pa_val(t_pargs *pa);
/* Return a pointer to a static buffer containing the value of pa */

extern int opt2parg_int(char *option,int nparg,t_pargs pa[]);

extern bool opt2parg_bool(char *option,int nparg,t_pargs pa[]);

extern real opt2parg_real(char *option,int nparg,t_pargs pa[]);

extern char *opt2parg_str(char *option,int nparg,t_pargs pa[]);

extern char *opt2parg_enum(char *option,int nparg,t_pargs pa[]);

extern bool opt2parg_bSet(char *option,int nparg,t_pargs pa[]);

extern void print_pargs(FILE *fp, int npargs,t_pargs pa[]);

extern void pr_enums(FILE *fp, int npargs,t_pargs pa[]);

#endif
