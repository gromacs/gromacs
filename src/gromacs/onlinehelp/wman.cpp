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

#include <time.h>

#include <string>

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/onlinehelp/wman.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/stringutil.h"

#include "gmx_fatal.h"
#include "string2.h"
#include "smalloc.h"
#include "sysstuff.h"
#include "filenm.h"
#include "macros.h"
#include "statutil.h"
#include "copyrite.h"
#include "strdb.h"
#include "readinp.h"

/* The source code in this file should be thread-safe.
         Please keep it that way. */


typedef struct {
    const char *search, *replace;
} t_sandr_const;

typedef struct {
    char *search, *replace;
} t_sandr;

/* The order of these arrays is significant. Text search and replace
 * for each element occurs in order, so earlier changes can induce
 * subsequent changes even though the original text might not appear
 * to invoke the latter changes. */

const t_sandr_const sandrTeX[] = {
    { "[TT]", "{\\tt " },
    { "[tt]", "}"      },
    { "[BB]", "{\\bf " },
    { "[bb]", "}"      },
    { "[IT]", "{\\em " },
    { "[it]", "}"      },
    { "[PAR]", "\n\n"   },
    /* Escaping underscore for LaTeX is no longer necessary, and it breaks
     * text searching and the index if you do. */
    /*
       { "_",    "\\_"    },
     */
    { "$",    "\\$"    },
    { "<=",   "\\ensuremath{\\leq{}}"},
    { ">=",   "\\ensuremath{\\geq{}}"},
    { "<",    "\\textless{}" },
    { ">",    "\\textgreater{}" },
    { "^",    "\\^{}"    },
    { "\\^{}t", "\\ensuremath{^t}" },
    { "\\^{}a", "\\ensuremath{^a}" },
    { "\\^{}b", "\\ensuremath{^b}" },
    { "\\^{}2", "\\ensuremath{^2}" },
    { "\\^{}3", "\\ensuremath{^3}" },
    { "\\^{}6", "\\ensuremath{^6}" },
    { "#",    "\\#"    },
    { "[BR]", "\\\\"   },
    { "%",    "\\%"    },
    { "&",    "\\&"    },
    /* The next couple of lines allow true Greek symbols to be written to the
       manual, which makes it look pretty */
    { "[GRK]", "\\ensuremath{\\" },
    { "[grk]", "}" },
    { "[MATH]", "\\ensuremath{" },
    { "[math]", "}" },
    { "[CHEVRON]", "\\ensuremath{<}" },
    { "[chevron]", "\\ensuremath{>}" },
    { "[MAG]", "\\ensuremath{|}" },
    { "[mag]", "\\ensuremath{|}" },
    { "[INT]", "\\ensuremath{\\int" },
    { "[FROM]", "_" },
    { "[from]", "" },
    { "[TO]", "^" },
    { "[to]", "" },
    { "[int]", "}" },
    { "[SUM]", "\\ensuremath{\\sum" },
    { "[sum]", "}" },
    { "[SUB]", "\\ensuremath{_{" },
    { "[sub]", "}}" },
    { "[SQRT]", "\\ensuremath{\\sqrt{" },
    { "[sqrt]", "}}" },
    { "[EXP]", "\\ensuremath{\\exp{(" },
    { "[exp]", ")}}" },
    { "[LN]", "\\ensuremath{\\ln{(" },
    { "[ln]", ")}}" },
    { "[LOG]", "\\ensuremath{\\log{(" },
    { "[log]", ")}}" },
    { "[COS]", "\\ensuremath{\\cos{(" },
    { "[cos]", ")}}" },
    { "[SIN]", "\\ensuremath{\\sin{(" },
    { "[sin]", ")}}" },
    { "[TAN]", "\\ensuremath{\\tan{(" },
    { "[tan]", ")}}" },
    { "[COSH]", "\\ensuremath{\\cosh{(" },
    { "[cosh]", ")}}" },
    { "[SINH]", "\\ensuremath{\\sinh{(" },
    { "[sinh]", ")}}" },
    { "[TANH]", "\\ensuremath{\\tanh{(" },
    { "[tanh]", ")}}" }
};
#define NSRTEX asize(sandrTeX)

const t_sandr_const sandrTty[] = {
    { "[TT]", "" },
    { "[tt]", "" },
    { "[BB]", "" },
    { "[bb]", "" },
    { "[IT]", "" },
    { "[it]", "" },
    { "[MATH]", "" },
    { "[math]", "" },
    { "[CHEVRON]", "<" },
    { "[chevron]", ">" },
    { "[MAG]", "|" },
    { "[mag]", "|" },
    { "[INT]", "integral" },
    { "[FROM]", " from " },
    { "[from]", "" },
    { "[TO]", " to " },
    { "[to]", " of" },
    { "[int]", "" },
    { "[SUM]", "sum" },
    { "[sum]", "" },
    { "[SUB]", "_" },
    { "[sub]", "" },
    { "[SQRT]", "sqrt(" },
    { "[sqrt]", ")" },
    { "[EXP]", "exp(" },
    { "[exp]", ")" },
    { "[LN]", "ln(" },
    { "[ln]", ")" },
    { "[LOG]", "log(" },
    { "[log]", ")" },
    { "[COS]", "cos(" },
    { "[cos]", ")" },
    { "[SIN]", "sin(" },
    { "[sin]", ")" },
    { "[TAN]", "tan(" },
    { "[tan]", ")" },
    { "[COSH]", "cosh(" },
    { "[cosh]", ")" },
    { "[SINH]", "sinh(" },
    { "[sinh]", ")" },
    { "[TANH]", "tanh(" },
    { "[tanh]", ")" },
    { "[PAR]", "\n\n" },
    { "[BR]", "\n"},
    { "[GRK]", "" },
    { "[grk]", "" }
};
#define NSRTTY asize(sandrTty)

const t_sandr_const sandrNROFF[] = {
    { "[TT]", "\\fB " },
    { "[tt]", "\\fR" },
    { "[BB]", "\\fB " },
    { "[bb]", "\\fR" },
    { "[IT]", "\\fI " },
    { "[it]", "\\fR" },
    { "[MATH]", "" },
    { "[math]", "" },
    { "[CHEVRON]", "<" },
    { "[chevron]", ">" },
    { "[MAG]", "|" },
    { "[mag]", "|" },
    { "[INT]", "integral" },
    { "[FROM]", " from " },
    { "[from]", "" },
    { "[TO]", " to " },
    { "[to]", " of" },
    { "[int]", "" },
    { "[SUM]", "sum" },
    { "[sum]", "" },
    { "[SUB]", "_" },
    { "[sub]", "" },
    { "[SQRT]", "sqrt(" },
    { "[sqrt]", ")", },
    { "[EXP]", "exp(" },
    { "[exp]", ")" },
    { "[LN]", "ln(" },
    { "[ln]", ")" },
    { "[LOG]", "log(" },
    { "[log]", ")" },
    { "[COS]", "cos(" },
    { "[cos]", ")" },
    { "[SIN]", "sin(" },
    { "[sin]", ")" },
    { "[TAN]", "tan(" },
    { "[tan]", ")" },
    { "[COSH]", "cosh(" },
    { "[cosh]", ")" },
    { "[SINH]", "sinh(" },
    { "[sinh]", ")" },
    { "[TANH]", "tanh(" },
    { "[tanh]", ")" },
    { "[PAR]", "\n\n" },
    { "\n ",    "\n" },
    { "<",    "" },
    { ">",    "" },
    { "^",    "" },
    { "#",    "" },
    { "[BR]", "\n"},
    { "-",    "\\-"},
    { "[GRK]", "" },
    { "[grk]", "" }
};
#define NSRNROFF asize(sandrNROFF)

const t_sandr_const sandrHTML[] = {
    { "<",    "&lt;" },
    { ">",    "&gt;" },
    { "[TT]", "<tt>" },
    { "[tt]", "</tt>" },
    { "[BB]", "<b>" },
    { "[bb]", "</b>" },
    { "[IT]", "<it>" },
    { "[it]", "</it>" },
    { "[MATH]", "" },
    { "[math]", "" },
    { "[CHEVRON]", "<" },
    { "[chevron]", ">" },
    { "[MAG]", "|" },
    { "[mag]", "|" },
    { "[INT]", "integral" },
    { "[FROM]", " from " },
    { "[from]", "" },
    { "[TO]", " to " },
    { "[to]", " of" },
    { "[int]", "" },
    { "[SUM]", "sum" },
    { "[sum]", "" },
    { "[SUB]", "_" },
    { "[sub]", "" },
    { "[SQRT]", "sqrt(" },
    { "[sqrt]", ")", },
    { "[EXP]", "exp(" },
    { "[exp]", ")" },
    { "[LN]", "ln(" },
    { "[ln]", ")" },
    { "[LOG]", "log(" },
    { "[log]", ")" },
    { "[COS]", "cos(" },
    { "[cos]", ")" },
    { "[SIN]", "sin(" },
    { "[sin]", ")" },
    { "[TAN]", "tan(" },
    { "[tan]", ")" },
    { "[COSH]", "cosh(" },
    { "[cosh]", ")" },
    { "[SINH]", "sinh(" },
    { "[sinh]", ")" },
    { "[TANH]", "tanh(" },
    { "[tanh]", ")" },
    { "[PAR]", "<p>" },
    { "[BR]", "<br>" },
    { "[GRK]", "&"  },
    { "[grk]", ";"  }
};
#define NSRHTML asize(sandrHTML)


static char *mydate(char buf[], int maxsize)
{
    const char *mon[] = {
        "Jan", "Feb", "Mar", "Apr", "May", "Jun",
        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
    };
    const char *day[] = { "Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat" };
    time_t      now;
    struct tm   tm;

    time(&now);
#ifdef GMX_NATIVE_WINDOWS
    /* Native windows */
    localtime_s(&tm, &now);
#else
    localtime_r(&now, &tm);
#endif

    sprintf(buf, "%s %d %s %d", day[tm.tm_wday], tm.tm_mday,
            mon[tm.tm_mon], tm.tm_year+1900);

    return buf;
}

/* Data structure for saved HTML links */
typedef struct t_linkdata {
    int      nsr;
    t_sandr *sr;
} t_linkdata;

static t_linkdata *init_linkdata()
{
    t_linkdata *p;
    snew(p, 1);
    p->sr  = NULL;
    p->nsr = 0;

    return p;
}

static void finish_linkdata(t_linkdata *p)
{
    int i;

    for (i = 0; i < p->nsr; i++)
    {
        sfree(p->sr[i].search);
        sfree(p->sr[i].replace);
    }
    sfree(p->sr);
    sfree(p);
}

static char *repall(const char *s, int nsr, const t_sandr_const sa[])
{
    try
    {
        std::string result(s);
        for (int i = 0; i < nsr; ++i)
        {
            result = gmx::replaceAll(result, sa[i].search, sa[i].replace);
        }
        return gmx_strdup(result.c_str());
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
}

static char *repallww(const char *s, int nsr, const t_sandr sa[])
{
    try
    {
        std::string result(s);
        for (int i = 0; i < nsr; ++i)
        {
            result = gmx::replaceAllWords(result, sa[i].search, sa[i].replace);
        }
        return gmx_strdup(result.c_str());
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
}

static char *html_xref(char *s, const char *program, t_linkdata *links)
{
    char   buf[256], **filestr;
    int    i, j, n;

    if (links->sr == NULL)
    {
        n          = get_file("links.dat", &(filestr));
        links->nsr = n;
        snew(links->sr, n);
        for (i = 0, j = 0; (i < n); i++)
        {
            if (!program || (gmx_strcasecmp(program, filestr[i])  != 0))
            {
                links->sr[j].search = gmx_strdup(filestr[i]);
                sprintf(buf, "<a href=\"%s.html\">%s</a>", filestr[i], filestr[i]);
                links->sr[j].replace = gmx_strdup(buf);
                j++;
            }
        }
        links->nsr = j;
        for (i = 0; i < n; i++)
        {
            sfree(filestr[i]);
        }
        sfree(filestr);
    }
    return repallww(s, links->nsr, links->sr);
}

char *check_tex(const char *s)
{
    return repall(s, NSRTEX, sandrTeX);
}

static char *check_nroff(const char *s)
{
    return repall(s, NSRNROFF, sandrNROFF);
}

static char *check_html(const char *s, const char *program, t_linkdata *links)
{
    char *buf;

    buf = repall(s, NSRHTML, sandrHTML);
    buf = html_xref(buf, program, links);

    return buf;
}

#define NSR(s) check_html(s, program, links)

#define FLAG_SET(flag, mask) ((flag &mask) == mask)
char *fileopt(unsigned long flag, char buf[], int maxsize)
{
    char tmp[256];

    if (FLAG_SET(flag, ffRW))
    {
        sprintf(tmp, "In/Out");
    }
    else if (FLAG_SET(flag, ffREAD))
    {
        sprintf(tmp, "Input");
    }
    else if (FLAG_SET(flag, ffWRITE))
    {
        sprintf(tmp, "Output");
    }
    else
    {
        sprintf(tmp, "Dunno");
    }

    if (FLAG_SET(flag, ffOPT))
    {
        strcat(tmp, ", Opt");
        if (FLAG_SET(flag, ffSET))
        {
            strcat(tmp, "!");
        }
        else
        {
            strcat(tmp, ".");
        }
    }
    if (FLAG_SET(flag, ffLIB))
    {
        strcat(tmp, ", Lib.");
    }
    if (FLAG_SET(flag, ffMULT))
    {
        strcat(tmp, ", Mult.");
    }

    sprintf(buf, "%s", tmp);

    return buf;
}

static void write_texman(FILE *out, const char *program,
                         int nldesc, const char **desc,
                         int nfile, t_filenm *fnm,
                         int npargs, t_pargs *pa,
                         int nbug, const char **bugs,
                         t_linkdata *links)
{
    int  i;
    char tmp[256];

    fprintf(out, "\\section{\\normindex{%s}}\\label{%s}\n\n", check_tex(program), check_tex(program));

    if (nldesc > 0)
    {
        for (i = 0; (i < nldesc); i++)
        {
            fprintf(out, "%s\n", check_tex(desc[i]));
        }
    }

    if (nfile > 0)
    {
        fprintf(out, "\\vspace{-2ex}\\begin{tabbing}\n");
        fprintf(out, "\n{\\normalsize \\bf Files}\\nopagebreak\\\\\n");
        fprintf(out, "{\\tt ~~~~~~~} \\= {\\tt ~~~~~~~~~~~~~~} \\= "
                "~~~~~~~~~~~~~~~~~~~~~~ \\= \\nopagebreak\\kill\n");
        for (i = 0; (i < nfile); i++)
        {
            fprintf(out, "\\>{\\tt %s} \\'\\> {\\tt %s} \\' %s \\> "
                    "\\parbox[t]{0.55\\linewidth}{%s} \\\\\n",
                    check_tex(fnm[i].opt), check_tex(fnm[i].fns[0]),
                    check_tex(fileopt(fnm[i].flag, tmp, 255)),
                    check_tex(ftp2desc(fnm[i].ftp)));
        }
        fprintf(out, "\\end{tabbing}\\vspace{-4ex}\n");
    }
    if (npargs > 0)
    {
        fprintf(out, "\\vspace{-2ex}\\begin{tabbing}\n");
        fprintf(out, "\n{\\normalsize \\bf Other options}\\nopagebreak\\\\\n");
        fprintf(out, "{\\tt ~~~~~~~~~~} \\= vector \\= "
                "{\\tt ~~~~~~~} \\= \\nopagebreak\\kill\n");
        for (i = 0; (i < npargs); i++)
        {
            if (strlen(check_tex(pa_val(&(pa[i]), tmp, 255))) <= 8)
            {
                fprintf(out, "\\> {\\tt %s} \\'\\> %s \\'\\> {\\tt %s} \\' "
                        "\\parbox[t]{0.68\\linewidth}{%s}\\\\\n",
                        check_tex(pa[i].option), get_arg_desc(pa[i].type),
                        check_tex(pa_val(&(pa[i]), tmp, 255)),
                        check_tex(pa[i].desc));
            }
            else
            {
                fprintf(out, "\\> {\\tt %s} \\'\\> %s \\'\\>\\\\\n"
                        "\\> \\'\\> \\'\\> {\\tt %s} \\' "
                        "\\parbox[t]{0.7\\linewidth}{%s}\\\\\n",
                        check_tex(pa[i].option), get_arg_desc(pa[i].type),
                        check_tex(pa_val(&(pa[i]), tmp, 255)),
                        check_tex(pa[i].desc));
            }
        }
        fprintf(out, "\\end{tabbing}\\vspace{-4ex}\n");
    }
    if (nbug > 0)
    {
        fprintf(out, "\n");
        fprintf(out, "\\begin{itemize}\n");
        for (i = 0; (i < nbug); i++)
        {
            fprintf(out, "\\item %s\n", check_tex(bugs[i]));
        }
        fprintf(out, "\\end{itemize}\n");
    }
/*   fprintf(out,"\n\\newpage\n"); */
}

static void write_nroffman(FILE *out,
                           const char *program,
                           int nldesc, const char **desc,
                           int nfile, t_filenm *fnm,
                           int npargs, t_pargs *pa,
                           int nbug, const char **bugs,
                           t_linkdata *links)

{
    int  i;
    char tmp[256];


    fprintf(out, ".TH %s 1 \"%s\" \"\" \"GROMACS suite, %s\"\n", program, mydate(tmp, 255), GromacsVersion());
    fprintf(out, ".SH NAME\n");
    fprintf(out, "%s@DESC@\n\n", program);
    fprintf(out, ".B %s\n", GromacsVersion());

    fprintf(out, ".SH SYNOPSIS\n");
    fprintf(out, "\\f3%s\\fP\n", program);

    /* command line arguments */
    if (nfile > 0)
    {
        for (i = 0; (i < nfile); i++)
        {
            fprintf(out, ".BI \"%s\" \" %s \"\n", check_nroff(fnm[i].opt),
                    check_nroff(fnm[i].fns[0]));
        }
    }
    if (npargs > 0)
    {
        for (i = 0; (i < npargs); i++)
        {
            if (pa[i].type == etBOOL)
            {
                fprintf(out, ".BI \"\\-[no]%s\" \"\"\n", check_nroff(pa[i].option+1));
            }
            else
            {
                fprintf(out, ".BI \"%s\" \" %s \"\n", check_nroff(pa[i].option),
                        check_nroff(get_arg_desc(pa[i].type)));
            }
        }
    }

    /* description */
    if (nldesc > 0)
    {
        fprintf(out, ".SH DESCRIPTION\n");
        for (i = 0; (i < nldesc); i++)
        {
            fprintf(out, "\\&%s\n", check_nroff(desc[i]));
        }
    }

    /* FILES */
    if (nfile > 0)
    {
        fprintf(out, ".SH FILES\n");
        for (i = 0; (i < nfile); i++)
        {
            fprintf(out, ".BI \"%s\" \" %s\" \n.B %s\n %s \n\n",
                    check_nroff(fnm[i].opt),
                    check_nroff(fnm[i].fns[0]),
                    check_nroff(fileopt(fnm[i].flag, tmp, 255)),
                    check_nroff(ftp2desc(fnm[i].ftp)));
        }
    }

    /* other options */
    fprintf(out, ".SH OTHER OPTIONS\n");
    if (npargs > 0)
    {
        for (i = 0; (i < npargs); i++)
        {
            if (pa[i].type == etBOOL)
            {
                fprintf(out, ".BI \"\\-[no]%s\"  \"%s\"\n %s\n\n",
                        check_nroff(pa[i].option+1),
                        check_nroff(pa_val(&(pa[i]), tmp, 255)),
                        check_nroff(pa[i].desc));
            }
            else
            {
                fprintf(out, ".BI \"%s\"  \" %s\" \" %s\" \n %s\n\n",
                        check_nroff(pa[i].option),
                        check_nroff(get_arg_desc(pa[i].type)),
                        check_nroff(pa_val(&(pa[i]), tmp, 255)),
                        check_nroff(pa[i].desc));
            }
        }
    }

    if (nbug > 0)
    {
        fprintf(out, ".SH KNOWN PROBLEMS\n");
        for (i = 0; (i < nbug); i++)
        {
            fprintf(out, "\\- %s\n\n", check_nroff(bugs[i]));
        }
    }

    fprintf(out, ".SH SEE ALSO\n.BR gromacs(7)\n\n");
    fprintf(out, "More information about \\fBGROMACS\\fR is available at <\\fIhttp://www.gromacs.org/\\fR>.\n");

}

char *check_tty(const char *s)
{
    return repall(s, NSRTTY, sandrTty);
}

static void
print_tty_formatted(FILE *out, int nldesc, const char **desc, int indent,
                    t_linkdata *links, const char *program)
{
    char *buf;
    char *temp;
    int   buflen, i;

    buflen = 80*nldesc;
    snew(buf, buflen);
    for (i = 0; (i < nldesc); i++)
    {
        if ((strlen(buf) > 0) &&
            (buf[strlen(buf)-1] != ' ') && (buf[strlen(buf)-1] != '\n'))
        {
            strcat(buf, " ");
        }
        temp = check_tty(desc[i]);
        if (strlen(buf) + strlen(temp) >= (size_t)(buflen-2))
        {
            buflen += strlen(temp);
            srenew(buf, buflen);
        }
        strcat(buf, temp);
        sfree(temp);
    }
    /* Make lines of at most 79 characters */
    temp = wrap_lines(buf, 78, indent, FALSE);
    fprintf(out, "%s\n", temp);
    sfree(temp);
    sfree(buf);
}

static void write_ttyman(FILE *out,
                         const char *program,
                         int nldesc, const char **desc,
                         int nfile, t_filenm *fnm,
                         int npargs, t_pargs *pa,
                         int nbug, const char **bugs,
                         t_linkdata *links)
{
    int   i;
    char *tmp;

    if (nldesc > 0)
    {
        fprintf(out, "DESCRIPTION\n-----------\n");
        print_tty_formatted(out, nldesc, desc, 0, links, program);
    }
    if (nbug > 0)
    {
        fprintf(out, "\n");
        fprintf(out, "KNOWN PROBLEMS\n----------\n");
        for (i = 0; i < nbug; i++)
        {
            snew(tmp, strlen(bugs[i])+3);
            strcpy(tmp, "* ");
            strcpy(tmp+2, check_tty(bugs[i]));
            fprintf(out, "%s\n", wrap_lines(tmp, 78, 2, FALSE));
            sfree(tmp);
        }
    }
    if (nfile > 0)
    {
        fprintf(out, "\n");
        pr_fns(out, nfile, fnm);
    }
    if (npargs > 0)
    {
        print_pargs(out, npargs, pa, FALSE);
    }
}

static void pr_html_files(FILE *out, int nfile, t_filenm fnm[],
                          const char *program, t_linkdata *links)
{
    int  i;
    char link[10], tmp[255];

    fprintf(out,
            "<TABLE BORDER=1 CELLSPACING=0 CELLPADDING=2>\n"
            "<TR>"
            "<TH>option</TH>"
            "<TH>filename</TH>"
            "<TH>type</TH>"
            "<TH>description</TH>"
            "</TR>\n");

    for (i = 0; (i < nfile); i++)
    {
        strcpy(link, ftp2ext(fnm[i].ftp));
        if (strcmp(link, "???") == 0)
        {
            strcpy(link, "files");
        }
        fprintf(out,
                "<TR>"
                "<TD ALIGN=RIGHT> <b><tt>%s</tt></b> </TD>"
                "<TD ALIGN=RIGHT> <tt><a href=\"%s.html\">%12s</a></tt> </TD>"
                "<TD> %s </TD>"
                "<TD> %s </TD>"
                "</TR>\n",
                fnm[i].opt, link, fnm[i].fns[0], fileopt(fnm[i].flag, tmp, 255),
                NSR(ftp2desc(fnm[i].ftp)));
    }

    fprintf(out, "</TABLE>\n");
}

static void write_htmlman(FILE *out,
                          const char *program,
                          int nldesc, const char **desc,
                          int nfile, t_filenm *fnm,
                          int npargs, t_pargs *pa,
                          int nbug, const char **bugs,
                          t_linkdata *links)
{
    int  i;
    char tmp[255];

    fprintf(out, "<HTML>\n<HEAD>\n<TITLE>%s</TITLE>\n", program);
    fprintf(out, "<LINK rel=stylesheet href=\"style.css\" type=\"text/css\">\n");
    fprintf(out, "<BODY text=\"#000000\" bgcolor=\"#FFFFFF\" link=\"#0000FF\" vlink=\"#990000\" alink=\"#FF0000\">\n");
    fprintf(out, "<TABLE WIDTH=\"98%%\" NOBORDER >\n<TR><TD WIDTH=400>\n");
    fprintf(out, "<TABLE WIDTH=400 NOBORDER>\n<TD WIDTH=116>\n");
    fprintf(out, "<a href=\"http://www.gromacs.org/\">"
            "<img SRC=\"../images/gmxlogo_small.png\""
            "BORDER=0 </a></td>\n");
    fprintf(out, "<td ALIGN=LEFT VALIGN=TOP WIDTH=280>"
            "<br><h2>%s</h2>", program);
    fprintf(out, "<font size=-1><A HREF=\"../online.html\">Main Table of Contents</A></font><br>");
    fprintf(out, "<br></td>\n</TABLE></TD><TD WIDTH=\"*\" ALIGN=RIGHT VALIGN=BOTTOM><p><B>%s<br>\n", GromacsVersion());
    fprintf(out, "%s</B></td></tr></TABLE>\n<HR>\n", mydate(tmp, 255));

    if (nldesc > 0)
    {
        fprintf(out, "<H3>Description</H3>\n<p>\n");
        for (i = 0; (i < nldesc); i++)
        {
            fprintf(out, "%s\n", NSR(desc[i]));
        }
    }
    if (nfile > 0)
    {
        fprintf(out, "<P>\n");
        fprintf(out, "<H3>Files</H3>\n");
        pr_html_files(out, nfile, fnm, program, links);
    }
    if (npargs > 0)
    {
        fprintf(out, "<P>\n");
        fprintf(out, "<H3>Other options</H3>\n");
        fprintf(out,
                "<TABLE BORDER=1 CELLSPACING=0 CELLPADDING=2>\n"
                "<TR>"
                "<TH>option</TH>"
                "<TH>type</TH>"
                "<TH>default</TH>"
                "<TH>description</TH>"
                "</TR>\n");
        for (i = 0; (i < npargs); i++)
        {
            fprintf(out,
                    "<TR>"
                    "<TD ALIGN=RIGHT> <b><tt>%s%s</tt></b> </TD>"
                    "<TD ALIGN=RIGHT> %s </TD>"
                    "<TD ALIGN=RIGHT> <tt>%s</tt> </TD>"
                    "<TD> %s </TD>"
                    "</TD>\n",
                    (pa[i].type == etBOOL) ? "-[no]" : "-", pa[i].option+1,
                    get_arg_desc(pa[i].type), pa_val(&(pa[i]), tmp, 255), NSR(pa[i].desc));
        }
        fprintf(out, "</TABLE>\n");
    }
    if (nbug > 0)
    {
        fprintf(out, "<P>\n");
        fprintf(out, "<H3>Known problems</H3>\n");
        fprintf(out, "<UL>\n");
        for (i = 0; (i < nbug); i++)
        {
            fprintf(out, "<LI>%s\n", NSR(bugs[i]));
        }
        fprintf(out, "</UL>\n");
    }
    fprintf(out, "<P>\n");
    fprintf(out, "<hr>\n<div ALIGN=RIGHT>\n");
    fprintf(out, "<font size=\"-1\"><a href=\"http://www.gromacs.org\">"
            "http://www.gromacs.org</a></font><br>\n");
    fprintf(out, "<font size=\"-1\"><a href=\"mailto:gromacs@gromacs.org\">"
            "gromacs@gromacs.org</a></font><br>\n");
    fprintf(out, "</div>\n");
    fprintf(out, "</BODY>\n");
}

static void pr_opts(FILE *fp,
                    int nfile,  t_filenm *fnm,
                    int npargs, t_pargs pa[], int shell)
{
    int i;

    switch (shell)
    {
        case eshellCSH:
            fprintf(fp, " \"c/-/(");
            for (i = 0; i < nfile; i++)
            {
                fprintf(fp, " %s", fnm[i].opt+1);
            }
            for (i = 0; i < npargs; i++)
            {
                if ( (pa[i].type == etBOOL) && *(pa[i].u.b) )
                {
                    fprintf(fp, " no%s", pa[i].option+1);
                }
                else
                {
                    fprintf(fp, " %s", pa[i].option+1);
                }
            }
            fprintf(fp, ")/\"");
            break;
        case eshellBASH:
            fprintf(fp, "if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W '");
            for (i = 0; i < nfile; i++)
            {
                fprintf(fp, " -%s", fnm[i].opt+1);
            }
            for (i = 0; i < npargs; i++)
            {
                if ( (pa[i].type == etBOOL) && *(pa[i].u.b) )
                {
                    fprintf(fp, " -no%s", pa[i].option+1);
                }
                else
                {
                    fprintf(fp, " -%s", pa[i].option+1);
                }
            }
            fprintf(fp, "' -- $c)); return 0; fi\n");
            break;
        case eshellZSH:
            fprintf(fp, " -x 's[-]' -s \"");
            for (i = 0; i < nfile; i++)
            {
                fprintf(fp, " %s", fnm[i].opt+1);
            }
            for (i = 0; i < npargs; i++)
            {
                if ( (pa[i].type == etBOOL) && *(pa[i].u.b) )
                {
                    fprintf(fp, " no%s", pa[i].option+1);
                }
                else
                {
                    fprintf(fp, " %s", pa[i].option+1);
                }
            }
            fprintf(fp, "\" ");
            break;
    }
}

static void write_cshcompl(FILE *out,
                           int nfile,  t_filenm *fnm,
                           int npargs, t_pargs *pa)
{
    fprintf(out, "complete %s", ShortProgram());
    pr_enums(out, npargs, pa, eshellCSH);
    pr_fopts(out, nfile, fnm, eshellCSH);
    pr_opts(out, nfile, fnm, npargs, pa, eshellCSH);
    fprintf(out, "\n");
}

static void write_zshcompl(FILE *out,
                           int nfile,  t_filenm *fnm,
                           int npargs, t_pargs *pa)
{
    fprintf(out, "compctl ");

    /* start with options, since they are always present */
    pr_opts(out, nfile, fnm, npargs, pa, eshellZSH);
    pr_enums(out, npargs, pa, eshellZSH);
    pr_fopts(out, nfile, fnm, eshellZSH);
    fprintf(out, "-- %s\n", ShortProgram());
}

static void write_bashcompl(FILE *out,
                            int nfile,  t_filenm *fnm,
                            int npargs, t_pargs *pa)
{
    /* Advanced bash completions are handled by shell functions.
     * p and c hold the previous and current word on the command line.
     * We need to use extended globbing, so write it in each completion file */
    fprintf(out, "shopt -s extglob\n");
    fprintf(out, "_%s_compl() {\nlocal p c\n", ShortProgram());
    fprintf(out, "COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}\n");
    pr_opts(out, nfile, fnm, npargs, pa, eshellBASH);
    fprintf(out, "case \"$p\" in\n");

    pr_enums(out, npargs, pa, eshellBASH);
    pr_fopts(out, nfile, fnm, eshellBASH);
    fprintf(out, "esac }\ncomplete -F _%s_compl %s\n", ShortProgram(), ShortProgram());
}

void write_man(const char *mantp,
               const char *program,
               int nldesc, const char **desc,
               int nfile, t_filenm *fnm,
               int npargs, t_pargs *pa,
               int nbug, const char **bugs,
               gmx_bool bHidden)
{
    int         i, npar;
    t_pargs    *par;

    t_linkdata *links;

    links = init_linkdata();

    const gmx::CommandLineHelpContext *context
        = gmx::GlobalCommandLineHelpContext::get();
    bool  bFileOpened = false;
    FILE *out;
    if (context != NULL)
    {
        out = context->writerContext().outputFile().handle();
    }
    else if (strcmp(mantp, "help") == 0)
    {
        out = stderr;
    }
    else
    {
        char buf[256];
        sprintf(buf, "%s.%s", program, mantp);
        out         = gmx_fio_fopen(buf, "w");
        bFileOpened = true;
    }

    if (bHidden)
    {
        npar = npargs;
        par  = pa;
    }
    else
    {
        snew(par, npargs);
        npar = 0;
        for (i = 0; i < npargs; i++)
        {
            if (!is_hidden(&pa[i]))
            {
                par[npar] = pa[i];
                npar++;
            }
        }
    }

    if (strcmp(mantp, "tex") == 0)
    {
        write_texman(out, program, nldesc, desc, nfile, fnm, npar, par, nbug, bugs, links);
    }
    if (strcmp(mantp, "nroff") == 0)
    {
        write_nroffman(out, program, nldesc, desc, nfile, fnm, npar, par, nbug, bugs, links);
    }
    if (strcmp(mantp, "help") == 0)
    {
        write_ttyman(out, program, nldesc, desc, nfile, fnm, npar, par, nbug, bugs, links);
    }
    if (strcmp(mantp, "html") == 0)
    {
        write_htmlman(out, program, nldesc, desc, nfile, fnm, npar, par, nbug, bugs, links);
    }
    if (strcmp(mantp, "completion-zsh") == 0)
    {
        write_zshcompl(out, nfile, fnm, npar, par);
    }
    if (strcmp(mantp, "completion-bash") == 0)
    {
        write_bashcompl(out, nfile, fnm, npar, par);
    }
    if (strcmp(mantp, "completion-csh") == 0)
    {
        write_cshcompl(out, nfile, fnm, npar, par);
    }

    if (bFileOpened)
    {
        gmx_fio_fclose(out);
    }

    if (!bHidden)
    {
        sfree(par);
    }

    finish_linkdata(links);
}

const char *get_arg_desc(int type)
{
    static const char *argtp[etNR] = {
        "int", "step", "real", "time", "string", "bool", "vector", "enum"
    };
    return argtp[type];
}
