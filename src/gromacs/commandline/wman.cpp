/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "gromacs/commandline/wman.h"

#include <cstdio>
#include <cstring>

#include <string>

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/fileio/filenm.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "gmx_fatal.h"
#include "string2.h"
#include "smalloc.h"

/* The source code in this file should be thread-safe.
 * Please keep it that way. */

static std::string check(const char *s, const gmx::HelpWriterContext &context)
{
    return context.substituteMarkupAndWrapToString(gmx::TextLineWrapperSettings(), s);
}

static std::string check(const char *s, const gmx::CommandLineHelpContext &context)
{
    return check(s, context.writerContext());
}

#define FLAG_SET(flag, mask) ((flag &mask) == mask)
/* Return a string describing the file type in flag.
 * flag should the flag field of a filenm struct.
 * You have to provide a buffer and buffer length in which
 * the result will be written. The returned pointer is just
 * a pointer to this buffer.
 */
static char *fileopt(unsigned long flag, char buf[])
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

#define OPTLEN 4
#define NAMELEN 14
static void pr_fns(FILE *fp, int nf, const t_filenm tfn[])
{
    int    i, f;
    size_t j;
    char   buf[256], *wbuf, opt_buf[32];

    fprintf(fp, "%6s %12s  %-12s %s\n", "Option", "Filename", "Type",
            "Description");
    fprintf(fp,
            "------------------------------------------------------------\n");
    for (i = 0; (i < nf); i++)
    {
        for (f = 0; (f < tfn[i].nfiles); f++)
        {
            sprintf(buf, "%4s %14s  %-12s ", (f == 0) ? tfn[i].opt : "",
                    tfn[i].fns[f], (f == 0) ? fileopt(tfn[i].flag, opt_buf)
                    : "");
            if (f < tfn[i].nfiles - 1)
            {
                fprintf(fp, "%s\n", buf);
            }
        }
        if (tfn[i].nfiles > 0)
        {
            strcat(buf, ftp2desc(tfn[i].ftp));
            if ((strlen(tfn[i].opt) > OPTLEN)
                && (strlen(tfn[i].opt) <= ((OPTLEN + NAMELEN)
                                           - strlen(tfn[i].fns[tfn[i].nfiles - 1]))))
            {
                for (j = strlen(tfn[i].opt); j < strlen(buf)
                     - (strlen(tfn[i].opt) - OPTLEN) + 1; j++)
                {
                    buf[j] = buf[j + strlen(tfn[i].opt) - OPTLEN];
                }
            }
            wbuf = wrap_lines(buf, 78, 35, FALSE);
            fprintf(fp, "%s\n", wbuf);
            sfree(wbuf);
        }
    }
    fprintf(fp, "\n");
    fflush(fp);
}
#undef OPTLEN
#undef NAMELEN

/* name to print in help info for command line arguments
 * (defined in enum in readinp.h) */
static const char *get_arg_desc(int type)
{
    const char *const argtp[etNR] = {
        "int", "step", "real", "time", "string", "bool", "vector", "enum"
    };
    return argtp[type];
}

/* Return the value of pa in the provided buffer buf, of size sz.
 * The return value is also a pointer to buf.
 */
static char *pa_val(t_pargs *pa, char buf[], int sz)
{
    real r;
    char buf_str[1256]; buf_str[0] = '\0';

    buf[0] = '\0';

    GMX_RELEASE_ASSERT(sz >= 255, "Buffer must be at least 255 chars");

    switch (pa->type)
    {
        case etINT:
            sprintf(buf, "%-d", *(pa->u.i));
            break;
        case etINT64:
            sprintf(buf, "%" GMX_PRId64, *(pa->u.is));
            break;
        case etTIME:
        case etREAL:
            r = *(pa->u.r);
            sprintf(buf_str, "%-6g", r);
            strcpy(buf, buf_str);
            break;
        case etBOOL:
            sprintf(buf, "%-6s", *(pa->u.b) ? "yes" : "no");
            break;
        case etSTR:
            if (*(pa->u.c))
            {
                if (strlen(*(pa->u.c)) >= (size_t)sz)
                {
                    gmx_fatal(FARGS, "Argument too long: \"%d\"\n", *(pa->u.c));
                }
                else
                {
                    strcpy(buf, *(pa->u.c));
                }
            }
            break;
        case etENUM:
            strcpy(buf, *(pa->u.c));
            break;
        case etRVEC:
            sprintf(buf, "%g %g %g", (*pa->u.rv)[0],
                    (*pa->u.rv)[1],
                    (*pa->u.rv)[2]);
            break;
    }
    return buf;
}

#define OPTLEN 12
#define TYPELEN 6
#define LONGSTR 1024
static char *pargs_print_line(t_pargs *pa, const gmx::HelpWriterContext &context)
{
    char buf[LONGSTR], *buf2, *tmp;

    snew(buf2, LONGSTR+strlen(pa->desc));
    snew(tmp, LONGSTR+strlen(pa->desc));

    if (pa->type == etBOOL)
    {
        sprintf(buf, "-[no]%s", pa->option+1);
    }
    else
    {
        strcpy(buf, pa->option);
    }
    std::string desc = check(pa->desc, context);
    if ((int)strlen(buf) > ((OPTLEN+TYPELEN)-std::max((int)strlen(get_arg_desc(pa->type)), 4)))
    {
        sprintf(buf2, "%s %-6s %-6s  %-s\n",
                buf, get_arg_desc(pa->type), pa_val(pa, tmp, LONGSTR-1),
                desc.c_str());
    }
    else if (strlen(buf) > OPTLEN)
    {
        /* so type can be 3 or 4 char's, this fits in the %4s */
        sprintf(buf2, "%-14s %-4s %-6s  %-s\n",
                buf, get_arg_desc(pa->type), pa_val(pa, tmp, LONGSTR-1),
                desc.c_str());
    }
    else
    {
        sprintf(buf2, "%-12s %-6s %-6s  %-s\n",
                buf, get_arg_desc(pa->type), pa_val(pa, tmp, LONGSTR-1),
                desc.c_str());
    }
    sfree(tmp);

    tmp = wrap_lines(buf2, 78, 28, FALSE);

    sfree(buf2);

    return tmp;
}
#undef OPTLEN
#undef TYPELEN
#undef LONGSTR

static void print_pargs(FILE *fp, int npargs, t_pargs pa[],
                        const gmx::HelpWriterContext &context)
{
    if (npargs > 0)
    {
        fprintf(fp, "%-12s %-6s %-6s  %-s\n",
                "Option", "Type", "Value", "Description");
        fprintf(fp, "------------------------------------------------------\n");
        for (int i = 0; i < npargs; i++)
        {
            char *wdesc = pargs_print_line(&pa[i], context);
            fprintf(fp, "%s", wdesc);
            sfree(wdesc);
        }
        fprintf(fp, "\n");
    }
}

static void write_nroffman(FILE *out,
                           int nldesc, const char **desc,
                           int nfile, t_filenm *fnm,
                           int npargs, t_pargs *pa,
                           int nbug, const char **bugs,
                           const gmx::CommandLineHelpContext &context)
{
    int  i;
    char tmp[256];

    fprintf(out, ".SH SYNOPSIS\n");
    fprintf(out, "\\f3%s\\fP\n", context.moduleDisplayName());

    /* command line arguments */
    if (nfile > 0)
    {
        for (i = 0; (i < nfile); i++)
        {
            fprintf(out, ".BI \"%s\" \" %s \"\n",
                    check(fnm[i].opt, context).c_str(),
                    check(fnm[i].fns[0], context).c_str());
        }
    }
    if (npargs > 0)
    {
        for (i = 0; (i < npargs); i++)
        {
            if (pa[i].type == etBOOL)
            {
                fprintf(out, ".BI \"\\-[no]%s\" \"\"\n",
                        check(pa[i].option+1, context).c_str());
            }
            else
            {
                fprintf(out, ".BI \"%s\" \" %s \"\n",
                        check(pa[i].option, context).c_str(),
                        check(get_arg_desc(pa[i].type), context).c_str());
            }
        }
    }

    /* description */
    if (nldesc > 0)
    {
        fprintf(out, ".SH DESCRIPTION\n");
        for (i = 0; (i < nldesc); i++)
        {
            fprintf(out, "\\&%s\n", check(desc[i], context).c_str());
        }
    }

    /* FILES */
    if (nfile > 0)
    {
        fprintf(out, ".SH FILES\n");
        for (i = 0; (i < nfile); i++)
        {
            fprintf(out, ".BI \"%s\" \" %s\" \n.B %s\n %s \n\n",
                    check(fnm[i].opt, context).c_str(),
                    check(fnm[i].fns[0], context).c_str(),
                    check(fileopt(fnm[i].flag, tmp), context).c_str(),
                    check(ftp2desc(fnm[i].ftp), context).c_str());
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
                        check(pa[i].option+1, context).c_str(),
                        check(pa_val(&(pa[i]), tmp, 255), context).c_str(),
                        check(pa[i].desc, context).c_str());
            }
            else
            {
                fprintf(out, ".BI \"%s\"  \" %s\" \" %s\" \n %s\n\n",
                        check(pa[i].option, context).c_str(),
                        check(get_arg_desc(pa[i].type), context).c_str(),
                        check(pa_val(&(pa[i]), tmp, 255), context).c_str(),
                        check(pa[i].desc, context).c_str());
            }
        }
    }

    if (nbug > 0)
    {
        fprintf(out, ".SH KNOWN PROBLEMS\n");
        for (i = 0; (i < nbug); i++)
        {
            fprintf(out, "\\- %s\n\n", check(bugs[i], context).c_str());
        }
    }
}

static void
print_tty_formatted(FILE *out, int nldesc, const char **desc,
                    const gmx::HelpWriterContext &context)
{
    char *buf;
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
        std::string temp = check(desc[i], context);
        if (strlen(buf) + temp.length() >= (size_t)(buflen-2))
        {
            buflen += temp.length();
            srenew(buf, buflen);
        }
        strcat(buf, temp.c_str());
    }
    /* Make lines of at most 79 characters */
    char *temp = wrap_lines(buf, 78, 0, FALSE);
    fprintf(out, "%s\n", temp);
    sfree(temp);
    sfree(buf);
}

static void write_ttyman(FILE *out,
                         int nldesc, const char **desc,
                         int nfile, t_filenm *fnm,
                         int npargs, t_pargs *pa,
                         int nbug, const char **bugs,
                         const gmx::HelpWriterContext &context)
{
    int   i;
    char *tmp;

    if (nldesc > 0)
    {
        fprintf(out, "DESCRIPTION\n-----------\n");
        print_tty_formatted(out, nldesc, desc, context);
    }
    if (nbug > 0)
    {
        fprintf(out, "\n");
        fprintf(out, "KNOWN PROBLEMS\n----------\n");
        for (i = 0; i < nbug; i++)
        {
            snew(tmp, strlen(bugs[i])+3);
            strcpy(tmp, "* ");
            strcpy(tmp+2, check(bugs[i], context).c_str());
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
        print_pargs(out, npargs, pa, context);
    }
}

static void pr_html_files(FILE *out, int nfile, t_filenm fnm[],
                          const gmx::HelpWriterContext &context)
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
                fnm[i].opt, link, fnm[i].fns[0], fileopt(fnm[i].flag, tmp),
                check(ftp2desc(fnm[i].ftp), context).c_str());
    }
    fprintf(out, "</TABLE>\n");
}

static void write_htmlman(FILE *out,
                          int nldesc, const char **desc,
                          int nfile, t_filenm *fnm,
                          int npargs, t_pargs *pa,
                          int nbug, const char **bugs,
                          const gmx::HelpWriterContext &context)
{
    int  i;
    char tmp[255];

    if (nldesc > 0)
    {
        fprintf(out, "<H3>Description</H3>\n<p>\n");
        for (i = 0; (i < nldesc); i++)
        {
            fprintf(out, "%s\n", check(desc[i], context).c_str());
        }
    }
    if (nfile > 0)
    {
        fprintf(out, "<P>\n");
        fprintf(out, "<H3>Files</H3>\n");
        pr_html_files(out, nfile, fnm, context);
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
                    get_arg_desc(pa[i].type), pa_val(&(pa[i]), tmp, 255),
                    check(pa[i].desc, context).c_str());
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
            fprintf(out, "<LI>%s\n", check(bugs[i], context).c_str());
        }
        fprintf(out, "</UL>\n");
    }
}

void write_man(const gmx::CommandLineHelpContext &context,
               int nldesc, const char **desc,
               int nfile, t_filenm *fnm,
               int npargs, t_pargs *pa,
               int nbug, const char **bugs)
{
    FILE       *out     = context.writerContext().outputFile().handle();
    const bool  bHidden = context.showHidden();

    int         npar;
    t_pargs    *par;

    if (bHidden)
    {
        npar = npargs;
        par  = pa;
    }
    else
    {
        snew(par, npargs);
        npar = 0;
        for (int i = 0; i < npargs; i++)
        {
            if (!is_hidden(&pa[i]))
            {
                par[npar] = pa[i];
                npar++;
            }
        }
    }

    switch (context.writerContext().outputFormat())
    {
        case gmx::eHelpOutputFormat_Man:
            write_nroffman(out, nldesc, desc, nfile, fnm, npar, par, nbug, bugs,
                           context);
            break;
        case gmx::eHelpOutputFormat_Console:
            write_ttyman(out, nldesc, desc, nfile, fnm, npar, par, nbug, bugs,
                         context.writerContext());
            break;
        case gmx::eHelpOutputFormat_Html:
            write_htmlman(out, nldesc, desc, nfile, fnm, npar, par, nbug, bugs,
                          context.writerContext());
            break;
        default:
            GMX_THROW(gmx::NotImplementedError("Help format not implemented"));
    }

    if (!bHidden)
    {
        sfree(par);
    }
}
