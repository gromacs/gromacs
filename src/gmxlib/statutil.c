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
static char *SRCID_statutil_c = "$Id$";

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
#include "assert.h"
#include "fatal.h"

/******************************************************************
 *
 *             T R A J E C T O R Y   S T U F F
 *
 ******************************************************************/

/* Globals for trajectory input */
real         tbegin=-1.0,tend=-1.0;
static bool  bView=FALSE;
static ulong uFlags=0;
static char  *program=NULL;

#define FF(arg) ((uFlags & arg)==arg)

static char *ShortProgram(void)
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

/****************************************************************
 *
 *            E X P O R T E D   F U N C T I O N S
 *
 ****************************************************************/

int check_times(real t) 
{
  if ((((tbegin >= 0.0) && (t >= tbegin)) || (tbegin == -1.0)) &&
      (((tend   >= 0.0) && (t <= tend))   || (tend   == -1.0))) {
    return 0;
  }
  else if ((tend != -1.0) && (t>=tend))
    return 1;
  return -1;
}

#ifdef NOTOBSOLETE

/**** C O O R D I N A T E S   A N D   V  E L O C  I T I E S ****/
int read_first_x_v(FILE *status,real *t,rvec **x,rvec **v,matrix box)
{
  t_statheader sh;
  real lambda;
  int  step,nre;

  fprintf(stderr,"\nReading statusfile, version: %s\n",rd_header(status,&sh));
  snew(*v,sh.natoms);
  snew(*x,sh.natoms);
  rd_hstatus(status,&sh,&step,t,&lambda,NULL,NULL,NULL,NULL,
	     &sh.natoms,*x,*v,NULL,
	     &nre,NULL,NULL);

  return sh.natoms;
}

bool read_next_x_v(FILE *status,real *t,int natoms,rvec x[],rvec v[],matrix box)
{
  t_statheader sh;
  real lambda;
  int  step,nre;
  bool bV;
  bool bX;

  while (!eof(status)) {
    rd_header(status,&sh);
    bX=sh.x_size;
    bV=sh.v_size;
    rd_hstatus(status,&sh,&step,t,&lambda,NULL,NULL,NULL,NULL,
	       &sh.natoms,bX ? x : NULL,bV ? v : NULL,NULL,
	       &nre,NULL,NULL);
    if ((check_times(*t)==0) && (bV) && (bX))
      return TRUE;
    if (check_times(*t) > 0)
      return FALSE;
  }
  return FALSE;
}

#endif

/***** T O P O L O G Y   S T U F F ******/

t_topology *read_top(char *fn)
{
  int        step,nre,natoms;
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
    top->blocks[i].multinr[0]=top->blocks[i].multinr[MAXPROC-1];
  for(i=0; (i<F_NRE); i++)
    top->idef.il[i].multinr[0]=top->idef.il[i].multinr[MAXPROC-1];
}

/*************************************************************
 *
 *           P A R S I N G   S T U F F
 *
 *************************************************************/
 
int iscan(int argc,char *argv[],int *i)
{
  int var;

  if (argc > (*i)+1) {
    if (!sscanf(argv[++(*i)],"%d",&var))
      usage(argv[0],argv[*i]);
  }
  else
    usage(argv[0],argv[*i]);

  return var;
}

double dscan(int argc,char *argv[],int *i)
{
  double var;

  if (argc > (*i)+1) {
    if (!sscanf(argv[++(*i)],"%lf",&var))
      usage(argv[0],argv[*i]);
  }
  else
    usage(argv[0],argv[*i]);

  return var;
}

char *sscan(int argc,char *argv[],int *i)
{
  /* 
     if (strlen(argv[*i]) > 2) 
     return &(argv[*i][2]);
     else
     */
  if (argc > (*i)+1) 
    return argv[++(*i)];
  else
    usage(argv[0],argv[*i]);

  return NULL; /* For the compiler! */
}

char *common_args(void)
{
  static char ca[1024];

  ca[0]='\0';  
#define CAT(s) { strcat(ca,s); strcat(ca," "); }

  if (FF(PCA_CAN_VIEW))  CAT("[-w (view it) ]");
  if (FF(PCA_CAN_BEGIN)) CAT("[-b begin time ]");
  if (FF(PCA_CAN_END))   CAT("[-e end time ]");
  CAT("[ -h (this message) ]");
#undef CAT    

  return ca;
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

void usage(char *prog,char *arg)
{
  if (arg != NULL)
    fprintf(stderr,"Conflicting argument for program %s: %s\n",prog,arg);
  exit(1);
}

bool bDoView(void)
{
  return bView;
}

static FILE *man_file(char *program,int mantp)
{
  FILE   *fp;
  char   buf[256],*pr;
  static char *manext[eotNR] = 
  { NULL, ".tex", ".html", ".1", ".txt", ".java" };
  
  if ((pr=strrchr(program,'/')) == NULL)
    pr=program;
  else 
    pr+=1;
    
  sprintf(buf,"%s%s",pr,manext[mantp]);
  fp = ffopen(buf,"w");
  
  return fp;
}

void parse_common_args(int *argc,char *argv[],ulong Flags,bool bNice,
		       int nfile,t_filenm fnm[],int npargs,t_pargs pa[],
		       int ndesc,char *desc[],int nbugs,char *bugs[])
{
  static bool bHelp=FALSE,bHidden=FALSE;
  static int  nicelevel=0;
  static int  mantp=0;
  static int  npri=0;
  static bool bDebug=FALSE;
  
  FILE *fp;  
  bool *bKeep,bPrint,bNo_get_pargs;
  int  i,k,b,npall;
  t_pargs *all_pa;
  int      nall_pa;
  t_pargs pca_pa[] = {
    { "-b",    FALSE, etREAL, &tbegin,        
      "first frame (ps) to read from trajectory" },
    { "-e",    FALSE, etREAL, &tend,        
      "last frame (ps) to read from trajectory" },
    { "-w",    FALSE, etBOOL, &bView,     
      "View output using xvgr or ghostview" },
    { "-nice", FALSE, etINT,  &nicelevel, 
      "Set the nicelevel" },
    { "-debug",FALSE, etBOOL, &bDebug,
      "HIDDENwrite file with debug information" },
    { "-h",    FALSE, etBOOL, &bHelp,     
      "Print help info and quit" },
    { "-hidden", FALSE, etBOOL, &bHidden,
      "HIDDENPrint hidden options" },
    { "-man",  FALSE, etINT,  &mantp,
      "HIDDENManual type: 0=none, 1=tex, 2=html, 3=nroff, 4=ascii, 5=java" }
#ifdef _SGI_
#ifdef USE_SGI_FPE
    ,
    { "-exception", FALSE, etBOOL, &bExcept,
      "HIDDENTurn on exception handling" }
#endif
    ,
    { "-npri", FALSE, etINT,  &npri,
      "Set non blocking priority (SGI only) try 250" }
#endif
  };
#define NPCA_PA asize(pca_pa)
  bool bFlags[NPCA_PA] = 
#ifdef _SGI_
#ifdef USE_SGI_FPE
  { 0, 0, 0, 1, 1, 1, 1, 1, 1, 1 };
#else
  { 0, 0, 0, 1, 1, 1, 1, 1, 1 };
#endif
#else
  { 0, 0, 0, 1, 1, 1, 1, 1 };
#endif
  
  /* First do file stuff */
  if (!program)
    program  = strdup(argv[0]);
  
  /* Now other options */
  uFlags  = Flags;
  bPrint  = (Flags & PCA_SILENT)   != PCA_SILENT;
  bNo_get_pargs = (Flags & PCA_NOGET_PARGS) != PCA_NOGET_PARGS;
  
  /* Check ALL the flags ... */
  if (bNice)
    nicelevel = 19;
    
  /* Check whether we have to add -b -e or -w options */
  bFlags[0] = (bool) (Flags & PCA_CAN_BEGIN);
  bFlags[1] = (bool) (Flags & PCA_CAN_END);
  bFlags[2] = (bool) (Flags & PCA_CAN_VIEW);
    
  snew(all_pa,NPCA_PA+npargs);
  for(i=npall=0; (i<NPCA_PA); i++)
    if (bFlags[i]) {
      memcpy(&(all_pa[npall]),&(pca_pa[i]),sizeof(pca_pa[i]));
      npall++;
    }
  for(i=0; (i<npargs); i++,npall++)
    memcpy(&(all_pa[npall]),&(pa[i]),sizeof(pa[i]));
  
  /* Parse the file args */
  parse_file_args(argc,argv,nfile,fnm,FF(PCA_KEEP_ARGS));

  /* Now parse the other options */
  get_pargs(argc,argv,npall,all_pa,FF(PCA_KEEP_ARGS));

  if (bNo_get_pargs) {  
    /* Now copy the results back... */
    for(i=0,k=npall-npargs; (i<npargs); i++,k++) 
      memcpy(&(pa[i]),&(all_pa[k]),sizeof(pa[i]));
  }
#ifdef _SGI_
#ifdef USE_SGI_FPE
  /* Install exception handler if necessary */
  if (bExcept)
    doexceptions();
#endif
#endif
  
#ifndef NO_NICE
  /* Set the nice level */
#ifdef _SGI_
  if (npri != 0) {
#include <sys/schedctl.h>
#include <sys/sysmp.h>
    (void) schedctl(MPTS_RTPRI,0,npri);
  }
  else
#endif
    if (nicelevel != 0)
      nice(nicelevel);
#endif

  if (bDebug) {
    char buf[256];
    
    sprintf(buf,"%s.log",ShortProgram());
    init_debug(buf);
    fprintf(debug,"%s (this file) opened in file %s, line %d\n",
	    buf,__FILE__,__LINE__);
  }
  
  if (!FF(PCA_QUIET)) {
    if (bHelp)
      write_man(stdout,eotHelp,program,ndesc,desc,nfile,fnm,npall,all_pa,nbugs,bugs,bHidden);
    else if (bPrint) {
      pr_fns(stdout,nfile,fnm);
      print_pargs(stdout,npall,all_pa);
    }
  }

  if (mantp != 0) {
    fp=man_file(program,mantp);
    write_man(fp,mantp,program,ndesc,desc,nfile,fnm,npall,all_pa,nbugs,bugs,bHidden);
    fclose(fp);
    exit(0);
  }
  sfree(all_pa);

  if (!FF(PCA_NOEXIT_ON_ARGS)) {
    if (*argc > 1) {
      for(i=1; (i<*argc); i++) 
	fprintf(stderr,"Unknown argument: %s\n",argv[i]);
      fprintf(stderr,"\nProgram %s halted\n",Program());
      exit(1);
    }
  } 
  if (bHelp)
    exit(0);
}
