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
  { "<",    "&lt;" },
  { ">",    "&gt;" },
  { "[TT]", "<tt>" },
  { "[tt]", "</tt>" },
  { "[BB]", "<b>" },
  { "[bb]", "</b>" },
  { "[IT]", "<it>" },
  { "[it]", "</it>" },
  { "[PAR]","<p>" },
  { "[BR]", "<br>" }
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
    if (i && buf)
      sfree(buf);
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
    if (i && buf)
      sfree(buf);
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
    strcat(buf,", Opt.");
  if ((flag & ffLIB) == ffLIB)
    strcat(buf,", Lib.");
    
  return buf;
}

static void write_texman(FILE *out,char *program,
			 int nldesc,char **desc,
			 int nfile,t_filenm *fnm,
			 int npargs,t_pargs *pa,
			 int nbug,char **bugs)
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
			   int nldesc,char **desc,
			   int nfile,t_filenm *fnm,
			   int npargs,t_pargs *pa,
			   int nbug,char **bugs)

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

void print_tty_formatted(FILE *out, int nldesc, char **desc)
{
  char *buf,*temp;
  int i,j;

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
  /* Make lines of at most 79 characters */
  temp = wrap_lines(buf,80,0);
  fprintf(out,"%s\n",temp);
  sfree(temp);
  sfree(buf);
}

static void write_ttyman(FILE *out,
			 char *program,
			 int nldesc,char **desc,
			 int nfile,t_filenm *fnm,
			 int npargs,t_pargs *pa,
			 int nbug,char **bugs,bool bHeader)
{
  int i;
  char *tmp;
  
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
    for(i=0; i<nbug; i++) {
      snew(tmp,strlen(bugs[i])+3);
      strcpy(tmp,"* ");
      strcpy(tmp+2,bugs[i]);
      fprintf(out,"%s\n",wrap_lines(tmp,80,2));
      sfree(tmp);
    }
  }
  if (nfile > 0) {
    fprintf(out,"\n");
    pr_fns(out,nfile,fnm);
  }
  if (npargs > 0) {
    print_pargs(out,npargs,pa);
  }
}

char *check_html(char *s,char *program)
{
  char *buf;
  
  buf=repall(s,NSRHTML,sandrHTML);
  buf=html_xref(buf,program);
  
  return buf;
}

static void write_htmlman(FILE *out,
			  char *program,
			  int nldesc,char **desc,
			  int nfile,t_filenm *fnm,
			  int npargs,t_pargs *pa,
			  int nbug,char **bugs)
{
  int i;
  char link[10];
  
#define NSR(s) check_html(s,program)
  
  fprintf(out,"<TITLE>%s</TITLE>\n",program);
  fprintf(out,"<BODY>\n");
  fprintf(out,"<H2>%s</H2>\n",program);
  fprintf(out,"<B>%s</B><BR>\n<B>%s</B><P>",GromacsVersion(),mydate());
  
  if (nldesc > 0) {
    fprintf(out,"<H3>Description</H3>\n");
    for(i=0; (i<nldesc); i++) 
      fprintf(out,"%s\n",NSR(desc[i]));
  }
  if (nfile > 0) {
    fprintf(out,"<P>\n");
    fprintf(out,"<H3>Files</H3>\n");
    fprintf(out,
	    "<TABLE BORDER=1 CELLSPACING=0 CELLPADDING=2>\n"
	    "<TR>"
	    "<TH>option</TH>"
	    "<TH>filename</TH>"
	    "<TH>type</TH>"
	    "<TH>description</TH>"
	    "</TR>\n");
    for(i=0; (i<nfile); i++) {
      strcpy(link,ftp2ext(fnm[i].ftp));
      if (strcmp(link,"???")==0)
	strcpy(link,"files");
      fprintf(out,
	      "<TR>"
	      "<TD ALIGN=RIGHT> <b>%s</b> </TD>"
	      "<TD ALIGN=RIGHT> <a href=\"%s.html\">%12s</a> </TD>"
	      "<TD> %s </TD>"
	      "<TD> %s </TD>"
	      "</TR>\n",
	      fnm[i].opt,link,fnm[i].fn,fileopt(fnm[i].flag),
	      NSR(ftp2desc(fnm[i].ftp)));
    }
    fprintf(out,"</TABLE>\n");
  }
  if (npargs > 0) {
    fprintf(out,"<P>\n");
    fprintf(out,"<H3>Other options</H3>\n");
    fprintf(out,
	    "<TABLE BORDER=1 CELLSPACING=0 CELLPADDING=2>\n"
	    "<TR>"
	    "<TH>option</TH>"
	    "<TH>type</TH>"
	    "<TH>default</TH>"
	    "<TH>description</TH>"
	    "</TR>\n");
    for(i=0; (i<npargs); i++)
      fprintf(out,
	      "<TR>"
	      "<TD ALIGN=RIGHT> <b>%s%s</b> </TD>"
	      "<TD ALIGN=RIGHT> %s </TD>"
	      "<TD ALIGN=RIGHT> <i>%s</i> </TD>"
	      "<TD> %s </TD>"
	      "</TD>\n",
	      (pa[i].type == etBOOL)?"-</b>[no]<b>":"-",pa[i].option+1,
	      argtp[pa[i].type],pa_val(&(pa[i])),NSR(pa[i].desc));
    fprintf(out,"</TABLE>\n");
  }
  if (nbug > 0) {
    fprintf(out,"<P>\n");
    fprintf(out,"<H3>Diagnostics</H3>\n");
    fprintf(out,"<UL>\n");
    for(i=0; (i<nbug); i++)
      fprintf(out,"<LI>%s\n",NSR(bugs[i]));
    fprintf(out,"</UL>\n");
  }
  fprintf(out,"<P>\n");
  fprintf(out,"</BODY>\n");
}

static void pr_opts(FILE *fp, 
		    int nfile,  t_filenm *fnm, 
		    int npargs, t_pargs pa[])
{
  int i;
  
  fprintf(fp," \"c/-/(");
  for (i=0; i<nfile; i++)
    fprintf(fp," %s",fnm[i].opt+1);
  for (i=0; i<npargs; i++)
    if ( (pa[i].type==etBOOL) && *(pa[i].u.b) )
      fprintf(fp," no%s",pa[i].option+1);
    else
      fprintf(fp," %s",pa[i].option+1);
  fprintf(fp,")/\"");
    
}

static void write_compl(FILE *out,
			char *program,
			int nldesc, char **desc,
			int nfile,  t_filenm *fnm,
			int npargs, t_pargs *pa,
			int nbug,   char **bugs)
{
  fprintf(out,"complete %s",ShortProgram());
  pr_enums(out,npargs,pa);
  pr_fopts(out,nfile,fnm);
  pr_opts(out,nfile,fnm,npargs,pa);
  fprintf(out,"\n");
}

void write_man(FILE *out,char *mantp,
	       char *program,
	       int nldesc,char **desc,
	       int nfile,t_filenm *fnm,
	       int npargs,t_pargs *pa,
	       int nbug,char **bugs,bool bHidden)
{
  char *pr;
  int     i,npar;
  t_pargs *par;
  
  if (bHidden || (strcmp(mantp,"completion")==0) ) {
    npar=npargs;
    par=pa;
  }    
  else {
    snew(par,npargs);
    npar=0;
    for(i=0;i<npargs;i++)
      if (!is_hidden(&pa[i])) {
	par[npar]=pa[i];
	npar++;
      }
  }
  
  if ((pr=strrchr(program,'/')) == NULL)
    pr=program;
  else
    pr+=1;
  if (strcmp(mantp,"tex")==0)
    write_texman(out,pr,nldesc,desc,nfile,fnm,npar,par,nbug,bugs);
  if (strcmp(mantp,"nroff")==0)
    write_nroffman(out,pr,nldesc,desc,nfile,fnm,npar,par,nbug,bugs);
  if (strcmp(mantp,"ascii")==0)
    write_ttyman(out,pr,nldesc,desc,nfile,fnm,npar,par,nbug,bugs,TRUE);
  if (strcmp(mantp,"help")==0)
    write_ttyman(out,pr,nldesc,desc,nfile,fnm,npar,par,nbug,bugs,FALSE);
  if (strcmp(mantp,"html")==0)
    write_htmlman(out,pr,nldesc,desc,nfile,fnm,npar,par,nbug,bugs);
  if (strcmp(mantp,"java")==0)
    write_java(out,pr,nldesc,desc,nfile,fnm,npar,par,nbug,bugs);
  if (strcmp(mantp,"completion")==0)
    write_compl(out,pr,nldesc,desc,nfile,fnm,npar,par,nbug,bugs);

  if (!bHidden)
    sfree(par);
}
