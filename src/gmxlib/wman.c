/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Great Red Owns Many ACres of Sand 
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
  { "\\^3",   "$^3$" },
  { "\\^6",   "$^6$" },
  { "#",    "\\#"    },
  { "[BR]", "\\\\"   },
  { "%",    "\\%"    },
  { "&",    "\\&"    },
  { "||",    "or"    },
  { "|",     "or"    }
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
  static char *mon[] = { "Jan", "Feb", "Mar", "Apr", "May", "Jun", 
			 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };
  static char *day[] = { "Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat" };
  time_t now;
  static char tbuf[128];
  struct tm *tm;
  
  (void) time(&now);
  tm = localtime(&now);
  sprintf(tbuf,"%s %d %s %d",day[tm->tm_wday],tm->tm_mday,
	  mon[tm->tm_mon],tm->tm_year+1900);
  
  return tbuf;
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

char *fileopt(unsigned long flag)
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
  if ((flag & ffOPT) == ffOPT) {
    strcat(buf,", Opt");
    if ((flag & ffSET) == ffSET) 
      strcat(buf,"!");
    else
      strcat(buf,".");
  }
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
  
  fprintf(out,"\\section{\\normindex{%s}}\n\n",check_tex(program));
  
  if (nldesc > 0)
    for(i=0; (i<nldesc); i++) 
      fprintf(out,"%s\n",check_tex(desc[i]));

  if (nfile > 0) {
    fprintf(out,"\\vspace{-2ex}\\begin{tabbing}\n");
    fprintf(out,"\n{\\normalsize \\bf Files}\\nopagebreak\\\\\n");
    fprintf(out,"{\\tt ~~~~~~~} \\= {\\tt ~~~~~~~~~~~~~~} \\= "
	    "~~~~~~~~~~~~~~~~~~~~~~ \\= \\nopagebreak\\kill\n");
    for(i=0; (i<nfile); i++)
      fprintf(out,"\\>{\\tt %s} \\'\\> {\\tt %s} \\' %s \\> "
	      "\\parbox[t]{0.55\\linewidth}{%s} \\\\\n",
	      check_tex(fnm[i].opt),check_tex(fnm[i].fn),
	      check_tex(fileopt(fnm[i].flag)),
	      check_tex(ftp2desc(fnm[i].ftp)));
    fprintf(out,"\\end{tabbing}\\vspace{-4ex}\n");
  }
  if (npargs > 0) {
    fprintf(out,"\\vspace{-2ex}\\begin{tabbing}\n");
    fprintf(out,"\n{\\normalsize \\bf Other options}\\nopagebreak\\\\\n");
    fprintf(out,"{\\tt ~~~~~~~~~~} \\= vector \\= "
	    "{\\tt ~~~~~~~} \\= \\nopagebreak\\kill\n");
    for(i=0; (i<npargs); i++) {
      if (strlen(check_tex(pa_val(&(pa[i])))) <= 8)
	fprintf(out,"\\> {\\tt %s} \\'\\> %s \\'\\> {\\tt %s} \\' "
		"\\parbox[t]{0.68\\linewidth}{%s}\\\\\n",
		check_tex(pa[i].option),argtp[pa[i].type],
		check_tex(pa_val(&(pa[i]))),
		check_tex(pa[i].desc));
      else
      	fprintf(out,"\\> {\\tt %s} \\'\\> %s \\'\\>\\\\\n"
		"\\> \\'\\> \\'\\> {\\tt %s} \\' "
		"\\parbox[t]{0.7\\linewidth}{%s}\\\\\n",
		check_tex(pa[i].option),argtp[pa[i].type],
		check_tex(pa_val(&(pa[i]))),
		check_tex(pa[i].desc));
    }
    fprintf(out,"\\end{tabbing}\\vspace{-4ex}\n");
  }
  if (nbug > 0) {
    fprintf(out,"\n");
    fprintf(out,"\\begin{itemize}\n");
    for(i=0; (i<nbug); i++)
      fprintf(out,"\\item %s\n",check_tex(bugs[i]));
    fprintf(out,"\\end{itemize}\n");
  }
/*   fprintf(out,"\n\\newpage\n"); */
}

static void write_nroffman(FILE *out,
			   char *program,
			   int nldesc,char **desc,
			   int nfile,t_filenm *fnm,
			   int npargs,t_pargs *pa,
			   int nbug,char **bugs)

{
  int i;
  
  fprintf(out,".TH %s 1 \"%s\"\n",program,mydate());
  fprintf(out,".SH NAME\n");
  fprintf(out,"%s\n",program);
  fprintf(out,".B %s\n",GromacsVersion());
  
  fprintf(out,".SH SYNOPSIS\n");
  fprintf(out,"\\f3%s\\fP\n",program);

  /* command line arguments */
  if (nfile > 0) {
    for(i=0; (i<nfile); i++)
      fprintf(out,".BI \"%s\" \" %s \"\n",fnm[i].opt,
	      fnm[i].fn);
  }
  if (npargs > 0) {
    for(i=0; (i<npargs); i++)
      if (pa[i].type == etBOOL)
	fprintf(out,".BI \"-[no]%s\" \"\"\n",pa[i].option+1);
      else
	fprintf(out,".BI \"%s\" \" %s \"\n",pa[i].option,
		argtp[pa[i].type]);
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
  fprintf(out,".SH OTHER OPTIONS\n");
  if ( npargs > 0 ) {
    for(i=0; (i<npargs); i++) {
      if (pa[i].type == etBOOL)
	fprintf(out,".BI \"-[no]%s\"  \"%s\"\n %s\n\n",
		check_nroff(pa[i].option+1),
		pa_val(&(pa[i])),check_nroff(pa[i].desc));
      else
	fprintf(out,".BI \"%s\"  \" %s\" \" %s\" \n %s\n\n",
		check_nroff(pa[i].option),argtp[pa[i].type],
		pa_val(&(pa[i])),check_nroff(pa[i].desc));
    }
  }

  if (nbug > 0) {
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
    fprintf(out,"\n");
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
  
  fprintf(out,"<HTML>\n<HEAD>\n<TITLE>%s</TITLE>\n",program);
  fprintf(out,"<LINK rel=stylesheet href=\"style.css\" type=\"text/css\">\n");
  fprintf(out,"<BODY text=\"#000000\" bgcolor=\"#FFFFFF\" link=\"#0000FF\" vlink=\"#990000\" alink=\"#FF0000\">\n");
  fprintf(out,"<table WIDTH=\"800\" NOBORDER >\n<TR>\n");
  fprintf(out,"<td WIDTH=\"120\" HEIGHT=\"133\">\n"
	  "<a href=\"http://www.gromacs.org/\">"
	  "<img SRC=\"../images/gmxlogo_small.jpg\""
	  "BORDER=0 height=133 width=116></a></td>");
  fprintf(out,"<td ALIGN=LEFT VALIGN=TOP WIDTH=480>"
	  "<br><br><h2>GROMACS Online Reference:<br>%s</h2>",program);
  fprintf(out,"<font size=-1><A HREF=\"../online.html\">Main Table of Contents</A></font><br>");
  fprintf(out,"<br></td>\n<TD ALIGN=RIGHT VALIGN=BOTTOM><B>%s<br>\n",GromacsVersion());
  fprintf(out,"%s</B></td></tr></TABLE>\n<HR>\n",mydate());
  
  if (nldesc > 0) {
    fprintf(out,"<H3>Description</H3>\n<p>\n");
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
	      "<TD ALIGN=RIGHT> <b><tt>%s</tt></b> </TD>"
	      "<TD ALIGN=RIGHT> <tt><a href=\"%s.html\">%12s</a></tt> </TD>"
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
	      "<TD ALIGN=RIGHT> <b><tt>%s%s</tt></b> </TD>"
	      "<TD ALIGN=RIGHT> %s </TD>"
	      "<TD ALIGN=RIGHT> <tt>%s</tt> </TD>"
	      "<TD> %s </TD>"
	      "</TD>\n",
	      (pa[i].type == etBOOL)?"-[no]":"-",pa[i].option+1,
	      argtp[pa[i].type],pa_val(&(pa[i])),NSR(pa[i].desc));
    fprintf(out,"</TABLE>\n");
  }
  if (nbug > 0) {
    fprintf(out,"<P>\n");
    fprintf(out,"<UL>\n");
    for(i=0; (i<nbug); i++)
      fprintf(out,"<LI>%s\n",NSR(bugs[i]));
    fprintf(out,"</UL>\n");
  }
  fprintf(out,"<P>\n");
  fprintf(out,"<hr>\n<div ALIGN=RIGHT>\n");
  fprintf(out,"<font size=\"-1\"><a href=\"http://www.gromacs.org\">"
	  "http://www.gromacs.org</a></font><br>\n");
  fprintf(out,"<font size=\"-1\"><a href=\"mailto:gromacs@gromacs.org\">"
	  "gromacs@gromacs.org</a></font><br>\n");
  fprintf(out,"</div>\n");
  fprintf(out,"</BODY>\n");
}

static void pr_opts(FILE *fp, 
		    int nfile,  t_filenm *fnm, 
		    int npargs, t_pargs pa[], int shell)
{
  int i;
  
  switch (shell) {
  case eshellCSH:
    fprintf(fp," \"c/-/(");
    for (i=0; i<nfile; i++)
      fprintf(fp," %s",fnm[i].opt+1);
    for (i=0; i<npargs; i++)
      if ( (pa[i].type==etBOOL) && *(pa[i].u.b) )
	fprintf(fp," no%s",pa[i].option+1);
      else
	fprintf(fp," %s",pa[i].option+1);
    fprintf(fp,")/\"");
    break;
  case eshellBASH:
    fprintf(fp,"if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W '");
    for (i=0; i<nfile; i++)
      fprintf(fp," -%s",fnm[i].opt+1);
    for (i=0; i<npargs; i++)
      if ( (pa[i].type==etBOOL) && *(pa[i].u.b) )
	fprintf(fp," -no%s",pa[i].option+1);
      else
	fprintf(fp," -%s",pa[i].option+1);
    fprintf(fp,"' -- $c)); return 0; fi\n");
    break;
  case eshellZSH:
    fprintf(fp," -x 's[-]' -s \"");
    for (i=0; i<nfile; i++)
      fprintf(fp," %s",fnm[i].opt+1);
    for (i=0; i<npargs; i++)
      if ( (pa[i].type==etBOOL) && *(pa[i].u.b) )
	fprintf(fp," no%s",pa[i].option+1);
      else
	fprintf(fp," %s",pa[i].option+1);
    fprintf(fp,"\" ");
    break;
  }
}

static void write_cshcompl(FILE *out,
			   int nfile,  t_filenm *fnm,
			   int npargs, t_pargs *pa)
{
  fprintf(out,"complete %s",ShortProgram());
  pr_enums(out,npargs,pa,eshellCSH);
  pr_fopts(out,nfile,fnm,eshellCSH);
  pr_opts(out,nfile,fnm,npargs,pa,eshellCSH);
  fprintf(out,"\n");
}

static void write_zshcompl(FILE *out,
			   int nfile,  t_filenm *fnm,
			   int npargs, t_pargs *pa)
{
  fprintf(out,"compctl ");

  /* start with options, since they are always present */
  pr_opts(out,nfile,fnm,npargs,pa,eshellZSH);
  pr_enums(out,npargs,pa,eshellZSH);
  pr_fopts(out,nfile,fnm,eshellZSH);
  fprintf(out,"-- %s\n",ShortProgram());
}

static void write_bashcompl(FILE *out,
			    int nfile,  t_filenm *fnm,
			    int npargs, t_pargs *pa)
{
  /* Advanced bash completions are handled by shell functions.
   * p and c hold the previous and current word on the command line.
   * We need to use extended globbing, so write it in each completion file */
  fprintf(out,"shopt -s extglob\n");
  fprintf(out,"_%s_compl() {\nlocal p c\n",ShortProgram());
  fprintf(out,"COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}\n");
  pr_opts(out,nfile,fnm,npargs,pa,eshellBASH);
  fprintf(out,"case \"$p\" in\n");
  pr_enums(out,npargs,pa,eshellBASH);
  pr_fopts(out,nfile,fnm,eshellBASH);
  fprintf(out,"esac }\ncomplete -F _%s_compl %s\n",ShortProgram(),ShortProgram());
}

void write_man(FILE *out,char *mantp,
	       char *program,
	       int nldesc,char **desc,
	       int nfile,t_filenm *fnm,
	       int npargs,t_pargs *pa,
	       int nbug,char **bugs,bool bHidden)
{
  char    *pr;
  int     i,npar;
  t_pargs *par;
 
  /* Don't write hidden options to completions, it just
   * makes the options more complicated for normal users
   */

  if (bHidden) {
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
  if (strcmp(mantp,"completion-zsh")==0)
    write_zshcompl(out,nfile,fnm,npar,par);
  if (strcmp(mantp,"completion-bash")==0)
    write_bashcompl(out,nfile,fnm,npar,par);
  if (strcmp(mantp,"completion-csh")==0)
    write_cshcompl(out,nfile,fnm,npar,par);

  if (!bHidden)
    sfree(par);
}
