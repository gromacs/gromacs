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
 * Gromacs Runs On Most of All Computer Systems
 */
static char *SRCID_wman_c = "$Id$";

#include "string2.h"
#include "smalloc.h"
#include "sysstuff.h"
#include "filenm.h"
#include "macros.h"
#include "replace.h"
#include "wman.h"
#include "statutil.h"
#include "copyrite.h"
#include "strdb.h"
#include "time.h"
#include "readinp.h"
#include "javaio.h"

static char *argtp[etNR] = { "int", "real", "string", "bool" };

typedef struct {
  char *search,*replace;
} t_sandr;

t_sandr sandrTeX[] = {
  { "[TT]", "{\\tt " },
  { "[tt]", "}"      },
  { "[BB]", "{\\bf " },
  { "[bb]", "}"      },
  { "[IT]", "{\\em " },
  { "[it]", "}"      },
  { "[PAR]","\n\n"   },
  { "_",    "\\_"    },
  { "$",    "\\$"    },
  { "<",    "$<$"    },
  { ">",    "$>$"    },
  { "^",    "\\^"    },
  { "\\^2",   "$^2$" },
  { "#",    "\\#"    },
  { "[BR]", "\\\\"   },
  { "%",    "\\%"    }
};
#define NSRTEX asize(sandrTeX)

t_sandr sandrTty[] = {
  { "[TT]", "" },
  { "[tt]", "" },
  { "[BB]", "" },
  { "[bb]", "" },
  { "[IT]", "" },
  { "[it]", "" },
  { "[PAR]","\n\n" },
/*   { "_",    "_" }, */
/*   { "<",    "" }, */
/*   { ">",    ">" }, */
/*   { "^",    "^" }, */
/*   { "#",    "#" }, */
/*   { "$",    "$"  }, */
  { "[BR]", "\n"}
};
#define NSRTTY asize(sandrTty)

t_sandr sandrNROFF[] = {
  { "[TT]", "\n.B " },
  { "[tt]", "\n" },
  { "[BB]", "\n.B " },
  { "[bb]", "\n" },
  { "[IT]", "\n.I " },
  { "[it]", "\n" },
  { "[PAR]","\n\n" },
  { "\n ",    "\n" },
  { "<",    "" },
  { ">",    "" },
  { "^",    "" },
  { "#",    "" },
  { "[BR]", "\n"}
};
#define NSRNROFF asize(sandrNROFF)

t_sandr sandrHTML[] = {
  { "[TT]", "<tt>" },
  { "[tt]", "</tt>" },
  { "[BB]", "<b>" },
  { "[bb]", "</b>" },
  { "[IT]", "<it>" },
  { "[it]", "</it>" },
  { "[PAR]", "<p>" },
  { "[BR]",  "<br>" }
};
#define NSRHTML asize(sandrHTML)

static char *mydate(void)
{
  time_t now;
  char   *cnow;
  
  now=time(NULL);
  
  cnow=strdup(ctime(&now));
  
  return cnow;
}

static char *repall(char *s,int nsr,t_sandr sa[])
{
  int  i;
  char *buf,*buf2;
  
  buf=s;
  for(i=0; (i<nsr); i++) {
    buf2=replace(buf,sa[i].search,sa[i].replace);
    buf=buf2;
  }
  
  return buf;
}

static char *repallww(char *s,int nsr,t_sandr sa[])
{
  int  i;
  char *buf,*buf2;
  
  buf=s;
  for(i=0; (i<nsr); i++) {
    buf2=replaceww(buf,sa[i].search,sa[i].replace);
    buf=buf2;
  }
  
  return buf;
}

char *check_tex(char *s)
{
  return repall(s,NSRTEX,sandrTeX);
}

char *check_nroff(char *s)
{
  return repall(s,NSRNROFF,sandrNROFF);
}

static char *html_xref(char *s,char *program)
{
  static int     nstr;
  static char    **str;
  static t_sandr *sr=NULL;
  char   buf[256];
  int    i,j;
  
  if (sr == NULL) {
    nstr=get_file("links.dat",&str);
    snew(sr,nstr);
    for(i=j=0; (i<nstr); i++) {
      if (!program || (strcasecmp(program,str[i]) != 0)) {
	sr[j].search=str[i];
	sprintf(buf,"<a href=\"%s.html\">%s</a>",str[i],str[i]);
	sr[j].replace=strdup(buf);
	j++;
      }
    }
    nstr=j;
  }
  return repallww(s,nstr,sr);
}

char *fileopt(ulong flag)
{
  static char buf[32];
  
  if ((flag & ffRW) == ffRW)
    strcpy(buf,"In/Out");
  else if ((flag & ffREAD) == ffREAD)
    strcpy(buf,"Input");
  else if ((flag & ffWRITE) == ffWRITE)
    strcpy(buf,"Output");
  else
    strcpy(buf,"Dunno");
  if ((flag & ffOPT) == ffOPT)
    strcat(buf,", Optional");
  if ((flag & ffLIB) == ffLIB)
    strcat(buf,", Library");
    
  return buf;
}

static void write_texman(FILE *out,char *program,
			 int nldesc,char *desc[],
			 int nfile,t_filenm fnm[],
			 int npargs,t_pargs pa[],
			 int nbug,char *bugs[])
{
  int i;

  fprintf(out,"\\newpage\n");
  fprintf(out,"\\section{%s.}\n",check_tex(program));
  fprintf(out,"{\\bf %s\\\\%s}\n\n",GromacsVersion(),mydate());
  
  if (nldesc > 0) {
    fprintf(out,"\\subsubsection*{Description.}\n");
    for(i=0; (i<nldesc); i++) 
      fprintf(out,"%s\n",check_tex(desc[i]));
  }
  if (nfile > 0) {
    fprintf(out,"\\subsubsection*{Files.}\n");
    fprintf(out,"\\begin{table}[htp]\n");
    fprintf(out,"\\begin{tabularx}{\\linewidth}{lllX}\n");
    for(i=0; (i<nfile); i++)
      fprintf(out,"%s & %s & %s & %s \\\\\n",
	      fnm[i].opt,fnm[i].fn,
	      fileopt(fnm[i].flag),check_tex(ftp2desc(fnm[i].ftp)));
    fprintf(out,"\\end{tabularx}\n");
    fprintf(out,"\\end{table}\n");
    fprintf(out,"Remember that filenames are not fixed, but \n");
    fprintf(out,"file extensions are.\n");
  }
  if (npargs > 0) {
    fprintf(out,"\\subsubsection*{Other options.}\n");
    fprintf(out,"\\begin{table}[htp]\n");
    fprintf(out,"\\begin{tabularx}{\\linewidth}{lllX}\n");
    for(i=0; (i<npargs); i++) {
      fprintf(out,"%s & %s & %s &%s\\\\\n",
	      check_tex(pa[i].option),argtp[pa[i].type],
	      check_tex(pa_val(&(pa[i]))),
	      check_tex(pa[i].desc));
    }
    fprintf(out,"\\end{tabularx}\n");
    fprintf(out,"\\end{table}\n");
  }
  if (nbug > 0) {
    fprintf(out,"\\subsubsection*{Diagnostics.}\n");
    fprintf(out,"\\begin{itemize}\n");
    for(i=0; (i<nbug); i++)
      fprintf(out,"\\item\t%s\n",check_tex(bugs[i]));
    fprintf(out,"\\end{itemize}\n");
  }
}

static void write_nroffman(FILE *out,
			   char *program,
			   int nldesc,char *desc[],
			   int nfile,t_filenm fnm[],
			   int npargs,t_pargs pa[],
			   int nbug,char *bugs[])

{
  int i; /* counter */

  fprintf(out,".TH %s 1 \"15 apr 2012\"\n",program);
  fprintf(out,".SH NAME\n");
  fprintf(out,"%s\n",program);
  fprintf(out,".B %s\n.B%s\n",GromacsVersion(),mydate());
  
  fprintf(out,".SH SYNOPSIS\n");
  fprintf(out,"\\f3%s\\fP\n",program);

  /* command line arguments */
  if (nfile > 0) {
    for(i=0; (i<nfile); i++)
      fprintf(out,".BI \"%s\" \" %s \"\n",fnm[i].opt,
	      fnm[i].fn);
  }
  
  /* description */
  if (nldesc > 0) {
    fprintf(out,".SH DESCRIPTION\n");
    for(i=0; (i<nldesc); i++) 
      fprintf(out,"%s\n",check_nroff(desc[i]));
  }

  /* FILES */
  if (nfile > 0) {
    fprintf(out,".SH FILES\n");
    for(i=0; (i<nfile); i++)
      fprintf(out,".BI \"%s\" \" %s\" \n.B %s\n %s \n\n",
	      fnm[i].opt,fnm[i].fn,fileopt(fnm[i].flag),
	      check_nroff(ftp2desc(fnm[i].ftp)));
    fprintf(out,"Remember that filenames are not fixed, but \n");
    fprintf(out,"file extensions are.\n");
  }
  
  /* other options */
  fprintf(out,".SH OPTIONS\n");
  if ( npargs > 0 ) {
    for(i=0; (i<npargs); i++) {
      fprintf(out,".BI %s & %s & %s &%s\\\\\n",
	      check_nroff(pa[i].option),argtp[pa[i].type],
	      pa_val(&(pa[i])),check_nroff(pa[i].desc));
    }
  }

  if (nbug > 0) {
    fprintf(out,".SH DIAGNOSTICS\n");
    for(i=0; (i<nbug); i++)
      fprintf(out,"\\- %s\n\n",check_nroff(bugs[i]));
  }

}

char *check_tty(char *s)
{
  return repall(s,NSRTTY,sandrTty);
}

void print_tty_formatted(FILE *out, int nldesc, char *desc[])
{
  char *buf,*temp;
  int i,j;
#define LINE_WIDTH 79

  /* Just to be sure */
  j=0;
  for(i=0; (i<nldesc); i++) 
    j+=strlen(desc[i])+10;
  snew(buf,j);
  for(i=0; (i<nldesc); i++) {
    if ((strlen(buf)>0) && 
	(buf[strlen(buf)-1] !=' ') && (buf[strlen(buf)-1] !='\n'))
      strcat(buf," ");
    temp=check_tty(desc[i]);
    strcat(buf,temp);
    sfree(temp);
  }
  i=0;
  while (i+LINE_WIDTH<strlen(buf)) {
    j=strchr(&(buf[i]),'\n')-buf;
    if ((j<i) || (j>i+LINE_WIDTH))
      j=i+LINE_WIDTH;
    while ((j>i) && (buf[j]!=' ') && (buf[j]!='\n'))
      j--;
    if (j!=i)
      buf[j]='\n';
    i=j+1;
  }
  fprintf(out,"%s\n",buf);
  sfree(buf);
}

static void write_ttyman(FILE *out,
			 char *program,
			 int nldesc,char *desc[],
			 int nfile,t_filenm fnm[],
			 int npargs,t_pargs pa[],
			 int nbug,char *bugs[],bool bHeader)
{
  if (bHeader) {
    fprintf(out,"%s\n\n",check_tty(program));
    fprintf(out,"%s\n%s\n",GromacsVersion(),mydate());
  }
  if (nldesc > 0) {
    fprintf(out,"DESCRIPTION:\n\n");
    print_tty_formatted(out,nldesc,desc);
  }
  if (nbug > 0) {
    fprintf(out,"\nDIAGNOSTICS\n");
    print_tty_formatted(out,nbug,bugs);
  }
  if (nfile > 0) {
    fprintf(out,"\n");
    pr_fns(out,nfile,fnm);
    fprintf(out,"Remember that filenames are not fixed, but file extensions are.\n");
  }
  if (npargs > 0) {
    fprintf(out,"\n");
    print_pargs(out,npargs,pa);
  }
}

char *check_html(char *s,char *program)
{
  char *buf;
  
  buf=html_xref(s,program);
  buf=repall(buf,NSRHTML,sandrHTML);
  
  return buf;
}

static void write_htmlman(FILE *out,
			  char *program,
			  int nldesc,char *desc[],
			  int nfile,t_filenm fnm[],
			  int npargs,t_pargs pa[],
			  int nbug,char *bugs[])
{
  int i;
  
#define NSR(s) check_html(s,program)
  
  fprintf(out,"<title>%s</title>\n",program);
  fprintf(out,"<h2>%s</h2>\n",program);
  fprintf(out,"<b>%s</b><br>\n<b>%s</b><p>",GromacsVersion(),mydate());
  
  if (nldesc > 0) {
    fprintf(out,"<h3>Description</h3>\n");
    for(i=0; (i<nldesc); i++) 
      fprintf(out,"%s\n",NSR(desc[i]));
  }
  if (nfile > 0) {
    fprintf(out,"<p>\n");
    fprintf(out,"<h3>Files</h3>\n");
    fprintf(out,"<dl>\n");
    for(i=0; (i<nfile); i++)
      fprintf(out,"<dt>%s <a href=\"%s.html\">%12s</a> <b>%s</b><dd>%s\n",
	      fnm[i].opt,
	      ftp2ext(fnm[i].ftp),
	      fnm[i].fn,fileopt(fnm[i].flag),
	      NSR(ftp2desc(fnm[i].ftp)));
    fprintf(out,"</dl>\n");
    fprintf(out,"Remember that filenames are not fixed, but \n");
    fprintf(out,"file extensions are.\n");
  }
  if (npargs > 0) {
    fprintf(out,"<p>\n");
    fprintf(out,"<h3>Other options</h3>\n");
    fprintf(out,"<dl>\n");
    for(i=0; (i<npargs); i++)
      fprintf(out,"<dt><b>%s</b> %s <i>%s</i><dd>%s\n",pa[i].option,
	      argtp[pa[i].type],pa_val(&(pa[i])),NSR(pa[i].desc));
    fprintf(out,"</dl>\n");
  }
  if (nbug > 0) {
    fprintf(out,"<p>\n");
    fprintf(out,"<h3>Diagnostics</h3>\n");
    fprintf(out,"<ul>\n");
    for(i=0; (i<nbug); i++)
      fprintf(out,"<li>%s\n",NSR(bugs[i]));
    fprintf(out,"</ul>\n");
  }
  fprintf(out,"<p>\n");
}

void write_man(FILE *out,int otype,
	       char *program,
	       int nldesc,char *desc[],
	       int nfile,t_filenm fnm[],
	       int npargs,t_pargs pa[],
	       int nbug,char *bugs[],bool bHidden)
{
  char *pr;
  int     i,npar;
  t_pargs *par;

  if (bHidden) {
    npar=npargs;
    par=pa;
  }    
  else {
    snew(par,npargs);
    npar=0;
    for(i=0;i<npargs;i++)
      if (strstr(pa[i].desc,"HIDDEN") == NULL) {
	par[npar]=pa[i];
	npar++;
      }
  }
  
  if ((pr=strrchr(program,'/')) == NULL)
    pr=program;
  else
    pr+=1;
  switch (otype) {
  case eotLaTeX:
    write_texman(out,pr,nldesc,desc,nfile,fnm,npar,par,nbug,bugs);
    break;
  case eotNroff:
    write_nroffman(out,pr,nldesc,desc,nfile,fnm,npar,par,nbug,bugs);
    break;
  case eotAscii:
    write_ttyman(out,pr,nldesc,desc,nfile,fnm,npar,par,nbug,bugs,TRUE);
    break;
  case eotHelp:
    write_ttyman(out,pr,nldesc,desc,nfile,fnm,npar,par,nbug,bugs,FALSE);
    break;
  case eotHTML:
    write_htmlman(out,pr,nldesc,desc,nfile,fnm,npar,par,nbug,bugs);
    break;
  case eotJava:
    write_java(out,pr,nldesc,desc,nfile,fnm,npar,par,nbug,bugs);
    break;
  default:
    break;
  }
  if (!bHidden)
    sfree(par);
}


