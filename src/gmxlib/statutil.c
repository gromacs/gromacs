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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <ctype.h>
#include "sysstuff.h"
#include "macros.h"
#include "string2.h"
#include "smalloc.h"
#include "pbc.h"
#include "statutil.h"
#include "names.h"
#include "vec.h"
#include "futil.h"
#include "wman.h"
#include "tpxio.h"
#include "fatal.h"
#include "network.h"
#include "vec.h"

/* used for npri */
#ifdef __sgi
#include <sys/schedctl.h>
#include <sys/sysmp.h>
#endif

/******************************************************************
 *
 *             T R A J E C T O R Y   S T U F F
 *
 ******************************************************************/

/* Globals for trajectory input */
real         tbegin = -1;
real         tend   = -1;
real         tdelta = -1;
real         timefactor = NOTSET;
char         *timelabel = NULL;
static char *timestr[] = { NULL, "fs", "ps", "ns", "us", "ms", "s",   
			   "m",              "h",                NULL };
real timefactors[]     = { 0,    1e3,  1,    1e-3, 1e-6, 1e-9, 1e-12, 
			   (1.0/60.0)*1e-12, (1.0/3600.0)*1e-12, 0 };
static char *xvgrtimestr[] = { NULL, "fs", "ps", "ns", "\\8m\\4s", "ms", "s",
			       "m", "h", NULL };
static bool  bView=FALSE;
static unsigned long uFlags=0;
static char  *program=NULL;

#define FF(arg) ((uFlags & arg)==arg)

char *ShortProgram(void)
{
  char *pr;
  
  if (program) {
    if ((pr=strrchr(program,'/')) != NULL)
      return pr+1;
    else
      return program;
  }
  else
    return "GROMACS";
}

char *Program(void)
{
  if (program)
    return program;
  else
    return "GROMACS";
}

void set_program_name(char *argvzero)
{
  /* When you run a dynamically linked program before installing
   * it, libtool uses wrapper scripts and prefixes the name with "lt-".
   * Until libtool is fixed to set argv[0] right, rip away the prefix:
   */
  if(program==NULL) {
    if(strlen(argvzero)>3 && !strncmp(argvzero,"lt-",3))
      program = strdup(argvzero+3);
    else
      program = strdup(argvzero);
  }
}

/****************************************************************
 *
 *            E X P O R T E D   F U N C T I O N S
 *
 ****************************************************************/

bool bRmod(double a, double b, double c)
{
  int iq;
  double tol;

  tol = 5*GMX_REAL_EPS;

  iq = (a - b + tol*a)/c;

  if (fabs(a - b - c*iq) <= tol*fabs(a))
    return TRUE;
  else
    return FALSE;
}

int check_times2(real t,real t0,real tp, real tpp)
{
  int  r;
  real margin;
  
  if (t-tp>0 && tp-tpp>0)
    margin = 0.1*min(t-tp,tp-tpp);
  else
    margin = 0;

  r=-1;
  if ((((tbegin >= 0.0) && (t >= tbegin)) || (tbegin == -1.0)) &&
      (((tend   >= 0.0) && (t <= tend+margin))   || (tend   == -1.0))) {
    if (tdelta > 0 && !bRmod(t,t0,tdelta))
      r = -1;
    else
      r = 0;
  }
  else if ((tend != -1.0) && (t>=tend))
    r = 1;
  if (debug) fprintf(debug,"t=%g, t0=%g, b=%g, e=%g, dt=%g: r=%d\n",
		     t,t0,tbegin,tend,tdelta,r);
  return r;
}

int check_times(real t)
{
  return check_times2(t,t,t,t);
}

char *time_unit(void)
{
  return timestr[0];
}

char *time_label(void)
{
  static char label[20];

  sprintf(label,"Time (%s)",timestr[0] ? timestr[0] : "ps");

  return label;
}

char *xvgr_tlabel(void)
{
  static char label[20];

  sprintf(label,"Time (%s)",
	  nenum(timestr) ? xvgrtimestr[nenum(timestr)] : "ps");

  return label;
}

static void init_time_factor()
{
  if (timefactor == NOTSET) 
    timefactor = timefactors[nenum(timestr)];
}

real time_factor(void)
{
  init_time_factor();
  
  return timefactor;
}

real convert_time(real time)
{
  init_time_factor();
  
  return (time*timefactor);
}


void convert_times(int n, real *time)
{
  int i;
  
  init_time_factor();
 
  if (timefactor!=1)
    for(i=0; i<n; i++)
      time[i] *= timefactor;
}

void default_time(void)
{
  timestr[0] = timestr[1];
  timefactor = timefactors[1];
  xvgrtimestr[0] = xvgrtimestr[1];
}

static void set_default_time_unit(char *select)
{
  int i,j;
  
  i=1;
  while(timestr[i] && strcmp(timestr[i], select)!=0)
    i++;
  if (strcmp(timestr[i], select)==0) {
    timestr[0] = timestr[i];
    timefactors[0] = timefactors[i];
    xvgrtimestr[0] = xvgrtimestr[i];
    for(j=i; j>1; j--) {
      timestr[j]=timestr[j-1];
      timefactors[j]=timefactors[j-1];
      xvgrtimestr[j]=xvgrtimestr[j-1];
    }
    timestr[1]=timestr[0];
    timefactors[1]=timefactors[0];
    xvgrtimestr[1]=xvgrtimestr[0];
  }
}

/***** T O P O L O G Y   S T U F F ******/

t_topology *read_top(char *fn)
{
  int        step,natoms;
  real       t,lambda;
  t_topology *top;

  snew(top,1);
  read_tpx(fn,&step,&t,&lambda,NULL,NULL,
	   &natoms,NULL,NULL,NULL,top);
  
  return top;
}

void mk_single_top(t_topology *top)
{
  int i;

  for(i=0; (i<ebNR); i++)
    top->blocks[i].multinr[0]=top->blocks[i].multinr[MAXNODES-1];
  for(i=0; (i<F_NRE); i++)
    top->idef.il[i].multinr[0]=top->idef.il[i].multinr[MAXNODES-1];
}

/*************************************************************
 *
 *           P A R S I N G   S T U F F
 *
 *************************************************************/

static void usage(char *type,char *arg)
{
  if (arg != NULL)
    gmx_fatal(FARGS,"Expected %s argument for option %s\n",type,arg);
}
 
int iscan(int argc,char *argv[],int *i)
{
  int var;

  if (argc > (*i)+1) {
    if (!sscanf(argv[++(*i)],"%d",&var))
      usage("an integer",argv[(*i)-1]);
  } else
    usage("an integer",argv[*i]);

  return var;
}

double dscan(int argc,char *argv[],int *i)
{
  double var;

  if (argc > (*i)+1) {
    if (!sscanf(argv[++(*i)],"%lf",&var))
      usage("a real",argv[(*i)-1]);
  } else
    usage("a real",argv[*i]);

  return var;
}

char *sscan(int argc,char *argv[],int *i)
{
  if (argc > (*i)+1) {
    if ( (argv[(*i)+1][0]=='-') && (argc > (*i)+2) && (argv[(*i)+2][0]!='-') )
      fprintf(stderr,"Possible missing string argument for option %s\n\n",
	      argv[*i]);
  } else
    usage("a string",argv[*i]);
  
  return argv[++(*i)];
}

int nenum(char *enumc[])
{
  int i;
  
  i=1;
  /* we *can* compare pointers directly here! */
  while(enumc[i] && enumc[0]!=enumc[i])
    i++;
  
  return i;
}

static void pdesc(char *desc)
{
  char *ptr,*nptr;
  
  ptr=desc;
  if ((int)strlen(ptr) < 70)
    fprintf(stderr,"\t%s\n",ptr);
  else {
    for(nptr=ptr+70; (nptr != ptr) && (!isspace(*nptr)); nptr--)
      ;
    if (nptr == ptr)
      fprintf(stderr,"\t%s\n",ptr);
    else {
      *nptr='\0';
      nptr++;
      fprintf(stderr,"\t%s\n",ptr);
      pdesc(nptr);
    }
  }
}

bool bDoView(void)
{
  return bView;
}

static FILE *man_file(char *program,char *mantp)
{
  FILE   *fp;
  char   buf[256],*pr;

  if ((pr=strrchr(program,'/')) == NULL)
    pr=program;
  else 
    pr+=1;
    
  if (strcmp(mantp,"ascii") != 0)
    sprintf(buf,"%s.%s",pr,mantp);
  else
    sprintf(buf,"%s.txt",pr);
  fp = ffopen(buf,"w");
  
  return fp;
}

static int add_parg(int npargs,t_pargs *pa,t_pargs *pa_add)
{
  memcpy(&(pa[npargs]),pa_add,sizeof(*pa_add));
  
  return npargs+1;
}

static char *mk_desc(t_pargs *pa, char *time_unit_str)
{
  char *newdesc=NULL,*ndesc=NULL,*ptr=NULL;
  int  len,k;
  
  /* First compute length for description */
  len = strlen(pa->desc)+1;
  if ((ptr = strstr(pa->desc,"HIDDEN")) != NULL)
    len += 4;
  if (pa->type == etENUM) {
    len += 10;
    for(k=1; (pa->u.c[k] != NULL); k++) {
      len += strlen(pa->u.c[k])+12;
    }
  }
  snew(newdesc,len);
  
  /* add label for hidden options */
  if (is_hidden(pa)) 
    sprintf(newdesc,"[hidden] %s",ptr+6);
  else
    strcpy(newdesc,pa->desc);
  
  /* change '%t' into time_unit */
#define TUNITLABEL "%t"
#define NTUNIT strlen(TUNITLABEL)
  if (pa->type == etTIME)
    while( (ptr=strstr(newdesc,TUNITLABEL)) != NULL ) {
      ptr[0]='\0';
      ptr+=NTUNIT;
      len+=strlen(time_unit_str)-NTUNIT;
      snew(ndesc,len);
      strcpy(ndesc,newdesc);
      strcat(ndesc,time_unit_str);
      strcat(ndesc,ptr);
      sfree(newdesc);
      newdesc=ndesc;
      ndesc=NULL;
    }
#undef TUNITLABEL
#undef NTUNIT
  
  /* Add extra comment for enumerateds */
  if (pa->type == etENUM) {
    strcat(newdesc,": ");
    for(k=1; (pa->u.c[k] != NULL); k++) {
      strcat(newdesc,"[TT]");
      strcat(newdesc,pa->u.c[k]);
      strcat(newdesc,"[tt]");
      /* Print a comma everywhere but at the last one */
      if (pa->u.c[k+1] != NULL) {
	if (pa->u.c[k+2] == NULL)
	  strcat(newdesc," or ");
	else
	  strcat(newdesc,", ");
      }
    }
  }
  return newdesc;
}

void parse_common_args(int *argc,char *argv[],unsigned long Flags,
		       int nfile,t_filenm fnm[],int npargs,t_pargs *pa,
		       int ndesc,char **desc,int nbugs,char **bugs)
{
  static bool bHelp=FALSE,bHidden=FALSE,bQuiet=FALSE;
  static char *manstr[]      = { NULL, "no", "html", "tex", "nroff", "ascii", "completion", "py", "xml", NULL };
  static char *not_nicestr[] = { NULL, "0", "4", "10", "19", NULL };
  static char *nicestr[]     = { NULL, "19", "10", "4", "0", NULL };
  static char *not_npristr[] = { NULL, "0", "128", "100", "200", "250", NULL };
  static char *npristr[]     = { NULL, "128", "250", "200", "100", "0", NULL };
  static int  nicelevel=0,mantp=0,npri=0;
  static bool bGUI=FALSE,bDebug=FALSE;
  static char *deffnm=NULL;
     
  t_pargs *all_pa=NULL;
  
  t_pargs motif_pa  = { "-X",    FALSE, etBOOL,  {&bGUI},
		       "Use dialog box GUI to edit command line options" };
  t_pargs npri_paX  = { "-npri", FALSE, etENUM,  {not_npristr},
		       "Set non blocking priority" };
  t_pargs npri_pa   = { "-npri", FALSE, etINT,   {&npri},
		       "HIDDEN Set non blocking priority (try 128)" };
  t_pargs nice_paX  = { "-nice", FALSE, etENUM,  {not_nicestr}, 
		       "Set the nicelevel" };
  t_pargs nice_pa   = { "-nice", FALSE, etINT,   {&nicelevel}, 
		       "Set the nicelevel" };
  t_pargs deffnm_pa = { "-deffnm", FALSE, etSTR, {&deffnm}, 
		       "Set the default filename for all file options" };
  t_pargs begin_pa  = { "-b",    FALSE, etTIME,  {&tbegin},        
		       "First frame (%t) to read from trajectory" };
  t_pargs end_pa    = { "-e",    FALSE, etTIME,  {&tend},        
		       "Last frame (%t) to read from trajectory" };
  t_pargs dt_pa     = { "-dt",   FALSE, etTIME,  {&tdelta},        
		       "Only use frame when t MOD dt = first time (%t)" };
  t_pargs view_pa   = { "-w",    FALSE, etBOOL,  {&bView},     
		       "View output xvg, xpm, eps and pdb files" };
  t_pargs time_pa   = { "-tu",   FALSE, etENUM,  {timestr},
			"Time unit" };
  /* Maximum number of extra arguments */
#define EXTRA_PA 16

  t_pargs pca_pa[] = {
    { "-h",    FALSE, etBOOL, {&bHelp},     
      "Print help info and quit" }, 
    { "-hidden", FALSE, etBOOL, {&bHidden},
      "HIDDENPrint hidden options" },
    { "-quiet",FALSE, etBOOL, {&bQuiet},
      "HIDDENDo not print help info" },
    { "-man",  FALSE, etENUM,  {manstr},
      "HIDDENWrite manual and quit" },
    { "-debug",FALSE, etBOOL, {&bDebug},
      "HIDDENWrite file with debug information" },
  };
#define NPCA_PA asize(pca_pa)
  FILE *fp;  
  bool bPrint,bExit;
  int  i,j,k,npall,max_pa;
  char *ptr,*newdesc;
  char *envstr;

  /* Check for double arguments */
  for (i=1; (i<*argc); i++) {
    if (argv[i] && (strlen(argv[i]) > 1) && (!isdigit(argv[i][1]))) {
      for (j=i+1; (j<*argc); j++) {
	if ( (argv[i][0]=='-') && (argv[j][0]=='-') && 
	     (strcmp(argv[i],argv[j])==0) ) {
	  if (FF(PCA_NOEXIT_ON_ARGS))
	    fprintf(stderr,"Double command line argument %s\n",argv[i]);
	  else
	    gmx_fatal(FARGS,"Double command line argument %s\n",argv[i]);
	}
      }
    }
  }
  debug_gmx();

  /* Handle the flags argument, which is a bit field 
   * The FF macro returns whether or not the bit is set
   */
  uFlags        = Flags;
  bPrint        = !FF(PCA_SILENT);
  
  /* Check whether we should have GUI or not */
#ifdef HAVE_MOTIF
  bGUI = (getenv("GMXMOTIF") != NULL);
  for(i=1; (i<*argc); i++) {
    if (strcmp(argv[i],"-X") == 0)
      bGUI = TRUE;
    else if (strcmp(argv[i],"-noX") == 0)
      bGUI = FALSE;
  }
  if (bGUI)
    bQuiet = TRUE;
#else
  bGUI = FALSE;
#endif

  set_program_name(argv[0]);

  /* Check ALL the flags ... */
  max_pa = NPCA_PA + EXTRA_PA + npargs;
  snew(all_pa,max_pa);
  
  for(i=npall=0; (i<NPCA_PA); i++)
    npall = add_parg(npall,all_pa,&(pca_pa[i]));

#ifdef HAVE_MOTIF
  /* Motif options */
  if (!bGUI)
    npall = add_parg(npall,all_pa,&motif_pa);
#endif

#ifdef __sgi
  envstr = getenv("GMXNPRIALL");
  if (envstr)
    npri=atoi(envstr);
  if (FF(PCA_BE_NICE)) {
    envstr = getenv("GMXNPRI");
    if (envstr)
      npri=atoi(envstr);
  }
  if (bGUI) {
    if (npri)
      npri_paX.u.c = npristr;
    npall = add_parg(npall,all_pa,&npri_paX);
  }
  else
    npall = add_parg(npall,all_pa,&npri_pa);
#endif

  if (bGUI) {
    /* Automatic nice or scheduling options */
    if (FF(PCA_BE_NICE)) 
      nice_paX.u.c = nicestr;
    npall = add_parg(npall,all_pa,&nice_paX);
  }
  else {
    if (FF(PCA_BE_NICE)) 
      nicelevel=19;
    npall = add_parg(npall,all_pa,&nice_pa);
  }

  if (FF(PCA_CAN_SET_DEFFNM)) 
    npall = add_parg(npall,all_pa,&deffnm_pa);   
  if (FF(PCA_CAN_BEGIN)) 
    npall = add_parg(npall,all_pa,&begin_pa);
  if (FF(PCA_CAN_END))
    npall = add_parg(npall,all_pa,&end_pa);
  if (FF(PCA_CAN_DT))
    npall = add_parg(npall,all_pa,&dt_pa);
  if (FF(PCA_TIME_UNIT)) {
    envstr = getenv("GMXTIMEUNIT");
    if ( envstr == NULL )
      envstr="ps";
    set_default_time_unit(envstr);
    npall = add_parg(npall,all_pa,&time_pa);
  } else
    set_default_time_unit("ps");
  if (FF(PCA_CAN_VIEW))
    npall = add_parg(npall,all_pa,&view_pa);

  /* Now append the program specific arguments */
  for(i=0; (i<npargs); i++)
    npall = add_parg(npall,all_pa,&(pa[i]));

  /* set etENUM options to default */
  for(i=0; (i<npall); i++)
    if (all_pa[i].type==etENUM)
      all_pa[i].u.c[0]=all_pa[i].u.c[1];
  
  /* Now parse all the command-line options */
  get_pargs(argc,argv,npall,all_pa,FF(PCA_KEEP_ARGS));

  if (FF(PCA_CAN_SET_DEFFNM) && (deffnm!=NULL))
    set_default_file_name(deffnm);

  /* Parse the file args */
  parse_file_args(argc,argv,nfile,fnm,FF(PCA_KEEP_ARGS));

  /* Open the debug file */
  if (bDebug) {
    char buf[256];

    if (gmx_node_num() > 1)
      sprintf(buf,"%s%d.log",ShortProgram(),gmx_node_id());
    else
      sprintf(buf,"%s.log",ShortProgram());
      
    init_debug(buf);
    fprintf(stderr,"Opening debug file %s (src code file %s, line %d)\n",
	    buf,__FILE__,__LINE__);
  }

  /* Now we have parsed the command line arguments. If the user wants it
   * we can now plop up a GUI dialog box to edit options.
   */
  if (bGUI) {
#ifdef HAVE_MOTIF
    gmx_gui(argc,argv,nfile,fnm,npall,all_pa,ndesc,desc,nbugs,bugs);
#else
    gmx_fatal(FARGS,"GROMACS compiled without MOTIF support - can't use X interface");
#endif
  }

  /* Now copy the results back... */
  for(i=0,k=npall-npargs; (i<npargs); i++,k++) 
    memcpy(&(pa[i]),&(all_pa[k]),(size_t)sizeof(pa[i]));
  
  for(i=0; (i<npall); i++)
    all_pa[i].desc = mk_desc(&(all_pa[i]), time_unit() );

  bExit = bHelp || (strcmp(manstr[0],"no") != 0);

#if (defined __sgi && USE_SGI_FPE)
  doexceptions();
#endif

  /* Set the nice level */
#ifdef __sgi
  if (bGUI)
    if (npri)
      sscanf(npristr[0],"%d",&npri);
    else
      sscanf(not_npristr[0],"%d",&npri);
  if (npri != 0 && !bExit) {
    (void) schedctl(MPTS_RTPRI,0,npri);
  }
  else
#endif 

#ifdef HAVE_UNISTD_H
    if (bGUI) {
      if (FF(PCA_BE_NICE))
	sscanf(nicestr[0],"%d",&nicelevel);
      else
	sscanf(not_nicestr[0],"%d",&nicelevel);
    }
  if (nicelevel != 0 && !bExit)
    nice(nicelevel);
#endif
  
  if (!(FF(PCA_QUIET) || bQuiet )) {
    if (bHelp)
      write_man(stderr,"help",program,ndesc,desc,nfile,fnm,npall,all_pa,
		nbugs,bugs,bHidden);
    else if (bPrint) {
      pr_fns(stderr,nfile,fnm);
      print_pargs(stderr,npall,all_pa);
    }
  }

  if (strcmp(manstr[0],"no") != 0) {
    if(!strcmp(manstr[0],"completion")) {
      /* one file each for csh, bash and zsh if we do completions */
      fp=man_file(program,"completion-zsh");
      write_man(fp,"completion-zsh",program,ndesc,desc,nfile,fnm,npall,all_pa,nbugs,bugs,bHidden);
      fclose(fp);
      fp=man_file(program,"completion-bash");
      write_man(fp,"completion-bash",program,ndesc,desc,nfile,fnm,npall,all_pa,nbugs,bugs,bHidden);
      fclose(fp);
      fp=man_file(program,"completion-csh");
      write_man(fp,"completion-csh",program,ndesc,desc,nfile,fnm,npall,all_pa,nbugs,bugs,bHidden);
      fclose(fp);
    } else {
      fp=man_file(program,manstr[0]);
      write_man(fp,manstr[0],program,ndesc,desc,nfile,fnm,npall,all_pa,nbugs,bugs,bHidden);
      fclose(fp);
    }
  }
  
  /* convert time options, must be done after printing! */
  init_time_factor();
  for(i=0; i<npall; i++) {
    if ((all_pa[i].type == etTIME) && (*all_pa[i].u.r >= 0)) {
      *all_pa[i].u.r /= timefactor;
    }
  }
  
  /* clear memory */
  for(i=0; i<npall; i++)
    sfree(all_pa[i].desc);
  sfree(all_pa);
  
  if (!FF(PCA_NOEXIT_ON_ARGS)) {
    if (*argc > 1) {
      gmx_cmd(argv[1]);
    }
  } 
  if (bExit) {
    if (gmx_parallel_env)
      gmx_abort(gmx_node_id(),gmx_node_num(),0);
    exit(0);
  }
}

