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
 * GROningen MAchine for Chemical Simulation
 */
static char *SRCID_filenm_c = "$Id$";

#include <string.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "string2.h"
#include "fatal.h"
#include "filenm.h"
#include "futil.h"
#include "wman.h"
#include "macros.h"

/* Use bitflag ... */
#define is_set(fn) ((fn.flag & ffSET) != 0)
#define is_opt(fn) ((fn.flag & ffOPT) != 0)
#define un_set(fn) (fn.flag = (fn.flag & ~ffSET))
#define do_set(fn) (fn.flag = (fn.flag |  ffSET))

typedef struct {
  bool bBinary;
  char *ext;
  char *defnm;
  char *defopt;
  char *descr;
} t_deffile;

static t_deffile deffile[efNR] = {
  { FALSE, ".mdp", "grompp", "-f", "grompp input file with MD parameters"    },
  { FALSE, ".gbp", "genbox", "-f", "genbox input file with parameters etc."  },
  { FALSE, ".gcp", "genconf","-f", "genconf input file with parameters etc." },
  { FALSE, ".gip", "genion", "-f", "genion input file with parameters etc."  },
  { FALSE, ".gdp", "gendock","-f", "gendock input file with parameters etc." },
  { FALSE, ".wdp", "widom",  "-f", "g_widom input parameters."               },
  { FALSE, ".gct", "gct",    "-f", "general coupling stuff"                  },
  { FALSE, ".gpp", "group",  "-g", "group related input for grompp"          },
  { TRUE,  ".???", "xtraj",  NULL, "Generic trajectory file format"          },
  { TRUE,  ".trr", "traj",   NULL, "Trajectory in portable xdr format" },
  { TRUE,  ".trj", "btraj",  NULL, "Trajectory file (cpu specific)"    },
  { TRUE , ".xtc", "ctraj",  NULL, "Compressed trajectory (portable xdr format)"},
  { FALSE, ".g87", "gtraj",  NULL, "Gromos-87 ASCII trajectory format" },
  { FALSE, ".asc", "atraj",  NULL, "Trajectory file in ascii converted by stat2asc" },
  { TRUE,  ".???", "ener",   NULL, "Generic energy format" },
  { TRUE,  ".edr", "ener",   NULL, "Energy file in portable XDR format"      },
  { TRUE,  ".ene", "ener",   NULL, "Energy file"                             },
  { FALSE, ".???", "xconf",  NULL, "Generic structure format"                },
  { FALSE, ".gro", "conf",   "-c", "Coordinate file in Gromos-87 format"     },
  { FALSE, ".pdb", "eiwit",  "-f", "Protein data bank file"                  },
  { FALSE, ".brk", "eiwit",  "-f", "Brookhaven data bank file"               },
  { FALSE, ".ent", "eiwit",  "-f", "Entry in the protein date bank"          },
  { TRUE , ".xdr", "cconf",  NULL, "Compressed coordinate file (portable xdr format)"},
  { FALSE, ".log", "run",    "-l", "Log file from MD/LD/EM/NM run"           },
  { FALSE, ".xvg", "graph",  "-o", "xvgr input file as produced by analysis tools"  },
  { FALSE, ".out", "hello",  "-o", "generic output file"                     },
  { FALSE, ".ndx", "index",  "-n", "GROMACS indexfile",                      },
  { FALSE, ".top", "topol",  "-p", "GROMACS topology file"                   },
  { FALSE, ".itp", "topinc", NULL, "Include file for GROMACS topology"       },
  { TRUE,  ".???", "topol",  "-s", "Generic GROMACS processed topology"      },
  { TRUE,  ".tpr", "topol",  "-s", "XDR GROMACS processed topology"          },
  { FALSE, ".tpa", "topol",  "-s", "Ascii GROMACS processed topology"        },
  { TRUE,  ".tpb", "topol",  "-s", "Binary GROMACS processed topology"       },
  { FALSE, ".tex", "doc",    "-o", "LaTeX file, suitable for inclusion in article" },
  { FALSE, ".rtp", "residue",NULL, "Residue Type file used by pdb2gmx"       },
  { FALSE, ".atp", "atomtp", NULL, "Atomtype file used by pdb2gmx"           },
  { FALSE, ".hdb", "polar",  NULL, "Hydrogen data base"                      },
  { FALSE, ".dat", "nnnice", NULL, "Generic data file"                       },
  { FALSE, ".dlg", "user",   NULL, "Dialog Box data for ngmx"                },
  { FALSE, ".gld", "atoms",  "-l", "Grompp Load Description"		     },
  { FALSE, ".map", "ss",     NULL, "File that maps matrix data to colors"    },
  { FALSE, ".eps", "plot",   NULL, "Encapsulated PostScript (tm) file"       },
  { FALSE, ".mat", "ss",     NULL, "Matrix Data file"			     },
  { FALSE, ".m2p", "ps",     NULL, "Input file for mat2ps"                   },
  { FALSE, ".vdw", "radii",  NULL, "Database containing Van der Waals radii" },
  { TRUE , ".mtx", "hessian","-m", "Hessian matrix"                          },
  { FALSE, ".vec", "eigvec", NULL, "(Eigen)Vectors"                          },
  { FALSE, ".edi", "sam",    NULL, "ED sampling input"                       },
  { FALSE, ".edo", "sam",    NULL, "ED sampling output"                      },
  { FALSE, ".hat", "gk",     NULL, "Fourier transform of spread function"    },
  { FALSE, ".xpm", "root",   NULL, "X PixMap compatible matrix file"         }
};

char *ftp2ext(int ftp)
{
  if ((0 <= ftp) && (ftp < efNR))
    return deffile[ftp].ext+1;
  else
    return NULL;
}

char *ftp2desc(int ftp)
{
  if ((0 <= ftp) && (ftp < efNR))
    return deffile[ftp].descr;
  else
    return NULL;
}

char *ftp2ftype(int ftp)
{
  if ((0 <= ftp) && (ftp < efNR)) {
    if (deffile[ftp].bBinary)
      return "Binary";
    else
      return "ASCII";
  }
  else
    return NULL;
}

char *ftp2defnm(int ftp)
{
  static char buf[256];
  
  if ((0 <= ftp) && (ftp < efNR)) {
    sprintf(buf,"%s%s",deffile[ftp].defnm,deffile[ftp].ext);
    return buf;
  }
  else
    return NULL;
}

void pr_def(FILE *fp,int ftp)
{
  t_deffile *df;
  
  df=&(deffile[ftp]);
  fprintf(fp,"%8s & %5s & %3s & %5s & ",
	  df->defnm,df->ext,df->bBinary ? "B" : "A",
	  df->defopt ? df->defopt : "");
  fprintf(fp,"%s \\\\\n",check_tex(df->descr));
}

void pr_fn(FILE *fp,t_filenm *tfn)
{
}

void pr_fns(FILE *fp,int nf,t_filenm tfn[])
{
  int  i,j;
  char buf[256],*wbuf;
#define OPTLEN 4
#define NAMELEN 14
  fprintf(fp,"%4s %14s  %-12s  %s\n",
	  "Opt","Filename","Type","Description");
  fprintf(fp,"------------------------------------------------------------\n");
  for(i=0; (i<nf); i++) {
    sprintf(buf, "%4s %14s  %-12s  %s\n", tfn[i].opt,tfn[i].fn,
	    fileopt(tfn[i].flag),deffile[tfn[i].ftp].descr);
    wbuf=wrap_lines(buf,80,35);
    if ( (strlen(tfn[i].opt)>OPTLEN) && 
	 (strlen(tfn[i].opt)<=((OPTLEN+NAMELEN)-strlen(tfn[i].fn))) ) {
      for(j=strlen(tfn[i].opt); 
	  j<strlen(wbuf)-(strlen(tfn[i].opt)-OPTLEN)+1; j++)
	wbuf[j]=wbuf[j+strlen(tfn[i].opt)-OPTLEN];
    }
    fprintf(fp,wbuf);
    sfree(wbuf);
  }
  fprintf(fp,"\n");
  fflush(fp);
}

static void check_opts(int nf,t_filenm fnm[])
{
  int       i;
  t_deffile *df;
  
  for(i=0; (i<nf); i++) {
    df=&(deffile[fnm[i].ftp]);
    if (fnm[i].opt == NULL) {
      if (df->defopt == NULL)
	fatal_error(0,"No default cmd-line option for %s",fnm[i].fn);
      else
	fnm[i].opt=df->defopt;
    }
  }
}

int fn2ftp(char *fn)
{
  int  i,len;
  char *feptr,*eptr;
  
  if (!fn)
    return efNR;

  len=strlen(fn);
  if ((len >= 4) && (fn[len-4] == '.'))
    feptr=&(fn[len-4]);
  else
    return efNR;
  
  for(i=0; (i<efNR); i++)
    if ((eptr=deffile[i].ext) != NULL)
      if (strcasecmp(feptr,eptr)==0)
	break;
      
  return i;
}

static void set_extension(char *buf,int ftp)
{
  int len;
  t_deffile *df;
  
  df=&(deffile[ftp]);
  len=strlen(buf);
  
  /* check if extension is already at end of filename */
  if (len >= 4) 
    if (strcasecmp(&(buf[len-4]),df->ext) == 0) 
      return;
  
  strcat(buf,df->ext);
}

void set_grpfnm(t_filenm *fnm,char *name,int nopts,int ftps[])
{
  char      buf[256];
  int       i,ti;

  /* First check whether we have a valid filename already */
  /* ti is the extension type */
  ti=fn2ftp(name); 
  
  for(i=0; (i<nopts); i++) {
    /* check if it is a trajectory */
    if (ti == ftps[i]) 
      break;
  }
  if (i != nopts) {
    /* If it is a trajectory file, return */
    fnm->fn=strdup(name);
    
    return;
  }
  for(i=0; (i<nopts); i++) {
    ti=ftps[i];
    if (name) {
      strcpy(buf,name);
      set_extension(buf,ti);
    }
    else
      strcpy(buf,ftp2defnm(ti));
    
    if (fexist(buf))
      break;
  }
  if (i == nopts) {
    ti=ftps[0];
    if (name) {
      strcpy(buf,name);
      set_extension(buf,ti);
    }
    else
      strcpy(buf,ftp2defnm(ti));
  }
  fnm->fn=strdup(buf);
}

static void set_trxnm(t_filenm *fnm,char *name)
{
  static    int trjs[]={ efTRR, efTRJ, efXTC, efG87, efPDB, efGRO };
#define NTRJS asize(trjs)

  set_grpfnm(fnm,name,NTRJS,trjs);
}

static void set_stxnm(t_filenm *fnm,char *name)
{
  static    int stxs[]={ efGRO, efPDB, efBRK, efENT, efTPB, efTPA, efTPX };
#define NSTXS asize(stxs)
  
  set_grpfnm(fnm,name,NSTXS,stxs);
}

static void set_enxnm(t_filenm *fnm,char *name)
{
  static    int enxs[]={ efEDR, efENE };
#define ENTXS asize(enxs)
  
  set_grpfnm(fnm,name,ENTXS,enxs);
}

static void set_tpxnm(t_filenm *fnm,char *name)
{
  static    int tpxs[]={ efTPR, efTPB, efTPA };
#define NTPXS asize(tpxs)
  
  set_grpfnm(fnm,name,NTPXS,tpxs);
}


static void set_filenm(t_filenm *fnm,char *name)
{
  /* Set the default filename, extension and option for those fields that 
   * are not already set. An extension is added if not present, if fn = NULL
   * or empty, the default filename is given.
   */
  char      buf[256];

  if ((fnm->ftp < 0) || (fnm->ftp >= efNR))
    fatal_error(0,"file type out of range (%d)",fnm->ftp);

  if (fnm->ftp == efTRX) {
    set_trxnm(fnm,name);
  }
  else if (fnm->ftp == efSTX) {
    set_stxnm(fnm,name);
  }
  else if (fnm->ftp == efTPX) {
    set_tpxnm(fnm,name);
  }
  else if (fnm->ftp == efENX) {
    set_enxnm(fnm,name);
  }
  else {
    if (name != NULL) {
      strcpy(buf,name);
    }
    else {
      strcpy(buf,deffile[fnm->ftp].defnm);
    }
    set_extension(buf,fnm->ftp);
    
    fnm->fn=strdup(buf);
  }
}

static void set_filenms(int nf,t_filenm fnm[])
{
  int i;

  for(i=0; (i<nf); i++)
    if (!is_set(fnm[i]))
      set_filenm(&(fnm[i]),fnm[i].fn);
}

void parse_file_args(int *argc,char *argv[],int nf,t_filenm fnm[],
		     bool bKeep)
{
  int  i,j;
  bool *bRemove;

  check_opts(nf,fnm);
  
  for(i=0; (i<nf); i++)
    un_set(fnm[i]);
    
  if (*argc > 1) {
    snew(bRemove,(*argc)+1);
    i=1;
    do {
      for(j=0; (j<nf); j++) {
	if (strcmp(argv[i],fnm[j].opt) == 0) {
	  do_set(fnm[j]);
	  bRemove[i]=TRUE;
	  i++;
	  if ((i < *argc) && (argv[i][0] != '-')) {
	    set_filenm(&fnm[j],argv[i]);
	    bRemove[i]=TRUE;
	    i++;
	  }
	  else
	    set_filenm(&fnm[j],fnm[j].fn);

	  break;
	}
      }
      /* No file found corresponding to option argv[i] */
      if (j == nf)
	i++;
    } while (i < *argc);
    
    if (!bKeep) {
      /* Remove used entries */
      for(i=j=0; (i<=*argc); i++) {
	if (!bRemove[i])
	  argv[j++]=argv[i];
      }
      (*argc)=j-1;
    }
    sfree(bRemove);
  }
  
  set_filenms(nf,fnm);
}

char *opt2fn(char *opt,int nfile,t_filenm fnm[])
{
  int i;
  
  for(i=0; (i<nfile); i++)
    if (strcmp(opt,fnm[i].opt)==0)
      return fnm[i].fn;

  fprintf(stderr,"No option %s\n",opt);
  return NULL;
}

char *ftp2fn(int ftp,int nfile,t_filenm fnm[])
{
  int i;
  
  for(i=0; (i<nfile); i++)
    if (ftp == fnm[i].ftp)
      return fnm[i].fn;
  
  fprintf(stderr,"ftp2fn: No filetype %s\n",deffile[ftp].ext);
  return NULL;
}

bool ftp2bSet(int ftp,int nfile,t_filenm fnm[])
{
  int i;
  
  for(i=0; (i<nfile); i++)
    if (ftp == fnm[i].ftp)
      return (bool) is_set(fnm[i]);
      
  fprintf(stderr,"ftp2bSet: No filetype %s\n",deffile[ftp].ext);
  
  return FALSE;
}

bool opt2bSet(char *opt,int nfile,t_filenm fnm[])
{
  int i;
  
  for(i=0; (i<nfile); i++)
    if (strcmp(opt,fnm[i].opt)==0)
      return (bool) is_set(fnm[i]);

  fprintf(stderr,"No option %s\n",opt);
  
  return FALSE;
}

char *opt2fn_null(char *opt,int nfile,t_filenm fnm[])
{
  int i;
  
  for(i=0; (i<nfile); i++)
    if (strcmp(opt,fnm[i].opt)==0)
      if (is_opt(fnm[i]) && !is_set(fnm[i]))
	return NULL;
      else
	return fnm[i].fn;

  fprintf(stderr,"No option %s\n",opt);
  return NULL;
}

char *ftp2fn_null(int ftp,int nfile,t_filenm fnm[])
{
  int i;
  
  for(i=0; (i<nfile); i++)
    if (ftp == fnm[i].ftp)
      if (is_opt(fnm[i]) && !is_set(fnm[i]))
	return NULL;
      else
	return fnm[i].fn;
  
  fprintf(stderr,"ftp2fn: No filetype %s\n",deffile[ftp].ext);
  return NULL;
}






