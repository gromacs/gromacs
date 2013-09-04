/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
#include "gromacs/commandline/shellcompletions.h"

#include <cstdio>

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/filenm.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"

#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/smalloc.h"

// Shell types for completion.
enum {
    eshellCSH, eshellBASH, eshellZSH
};

// TODO: Don't duplicate this from filenm.c/futil.c.
#define NZEXT 2
static const char *z_ext[NZEXT] = { ".gz", ".Z" };

static void pr_fopts(FILE *fp, int nf, const t_filenm tfn[], int shell)
{
    switch (shell)
    {
        case eshellCSH:
            for (int i = 0; i < nf; i++)
            {
                fprintf(fp, " \"n/%s/f:*.", tfn[i].opt);
                const int ftp          = tfn[i].ftp;
                const int genericCount = ftp2generic_count(ftp);
                if (genericCount > 0)
                {
                    fprintf(fp, "{");
                    const int *const genericTypes = ftp2generic_list(ftp);
                    for (int j = 0; j < genericCount; j++)
                    {
                        if (j > 0)
                        {
                            fprintf(fp, ",");
                        }
                        fprintf(fp, "%s", ftp2ext(genericTypes[j]));
                    }
                    fprintf(fp, "}");
                }
                else
                {
                    fprintf(fp, "%s", ftp2ext(ftp));
                }
                fprintf(fp, "{");
                for (int j = 0; j < NZEXT; j++)
                {
                    fprintf(fp, ",%s", z_ext[j]);
                }
                fprintf(fp, "}/\"");
            }
            break;
        case eshellBASH:
            for (int i = 0; i < nf; i++)
            {
                fprintf(fp, "%s) COMPREPLY=( $(compgen -X '!*.", tfn[i].opt);
                const int ftp          = tfn[i].ftp;
                const int genericCount = ftp2generic_count(ftp);
                if (genericCount > 0)
                {
                    fprintf(fp, "+(");
                    const int *const genericTypes = ftp2generic_list(ftp);
                    for (int j = 0; j < genericCount; j++)
                    {
                        if (j > 0)
                        {
                            fprintf(fp, "|");
                        }
                        fprintf(fp, "%s", ftp2ext(genericTypes[j]));
                    }
                    fprintf(fp, ")");
                }
                else
                {
                    fprintf(fp, "%s", ftp2ext(ftp));
                }
                fprintf(fp, "*(");
                for (int j = 0; j < NZEXT; j++)
                {
                    if (j > 0)
                    {
                        fprintf(fp, "|");
                    }
                    fprintf(fp, "%s", z_ext[j]);
                }
                fprintf(fp, ")' -f $c ; compgen -S '/' -X '.*' -d $c ));;\n");
            }
            break;
        case eshellZSH:
            for (int i = 0; i < nf; i++)
            {
                fprintf(fp, "- 'c[-1,%s]' -g '*.", tfn[i].opt);
                const int ftp          = tfn[i].ftp;
                const int genericCount = ftp2generic_count(ftp);
                if (genericCount > 0)
                {
                    fprintf(fp, "(");
                    const int *const genericTypes = ftp2generic_list(ftp);
                    for (int j = 0; j < genericCount; j++)
                    {
                        if (j > 0)
                        {
                            fprintf(fp, "|");
                        }
                        fprintf(fp, "%s", ftp2ext(genericTypes[j]));
                    }
                    fprintf(fp, ")");
                }
                else
                {
                    fprintf(fp, "%s", ftp2ext(ftp));
                }
                fprintf(fp, "(");
                for (int j = 0; j < NZEXT; j++)
                {
                    fprintf(fp, "|%s", z_ext[j]);
                }
                fprintf(fp, ") *(/)' ");
            }
            break;
    }
}

static void pr_opts(FILE *fp,
                    int nfile,  t_filenm *fnm,
                    int npargs, t_pargs pa[], int shell)
{
    switch (shell)
    {
        case eshellCSH:
            fprintf(fp, " \"c/-/(");
            for (int i = 0; i < nfile; i++)
            {
                fprintf(fp, " %s", fnm[i].opt+1);
            }
            for (int i = 0; i < npargs; i++)
            {
                if (pa[i].type == etBOOL && *(pa[i].u.b))
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
            for (int i = 0; i < nfile; i++)
            {
                fprintf(fp, " -%s", fnm[i].opt+1);
            }
            for (int i = 0; i < npargs; i++)
            {
                if (pa[i].type == etBOOL && *(pa[i].u.b))
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
            for (int i = 0; i < nfile; i++)
            {
                fprintf(fp, " %s", fnm[i].opt+1);
            }
            for (int i = 0; i < npargs; i++)
            {
                if (pa[i].type == etBOOL && *(pa[i].u.b))
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

static void pr_enums(FILE *fp, int npargs, t_pargs pa[], int shell)
{
    int i, j;

    switch (shell)
    {
        case eshellCSH:
            for (i = 0; i < npargs; i++)
            {
                if (pa[i].type == etENUM)
                {
                    fprintf(fp, " \"n/%s/(", pa[i].option);
                    for (j = 1; pa[i].u.c[j]; j++)
                    {
                        fprintf(fp, " %s", pa[i].u.c[j]);
                    }
                    fprintf(fp, ")/\"");
                }
            }
            break;
        case eshellBASH:
            for (i = 0; i < npargs; i++)
            {
                if (pa[i].type == etENUM)
                {
                    fprintf(fp, "%s) COMPREPLY=( $(compgen -W '", pa[i].option);
                    for (j = 1; pa[i].u.c[j]; j++)
                    {
                        fprintf(fp, " %s", pa[i].u.c[j]);
                    }
                    fprintf(fp, " ' -- $c ));;\n");
                }
            }
            break;
        case eshellZSH:
            for (i = 0; i < npargs; i++)
            {
                if (pa[i].type == etENUM)
                {
                    fprintf(fp, "- 'c[-1,%s]' -s \"", pa[i].option);
                    for (j = 1; pa[i].u.c[j]; j++)
                    {
                        fprintf(fp, " %s", pa[i].u.c[j]);
                    }
                    fprintf(fp, "\" ");
                }
            }
            break;
    }
}

static void write_bashcompl(FILE *out,
                            const char *moduleName,
                            int nfile,  t_filenm *fnm,
                            int npargs, t_pargs *pa)
{
    /* Advanced bash completions are handled by shell functions.
     * p and c hold the previous and current word on the command line.
     * We need to use extended globbing, so write it in each completion file */
    fprintf(out, "_%s_compl() {\n", moduleName);
    fprintf(out, "local p c\n");
    fprintf(out, "COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}\n");
    pr_opts(out, nfile, fnm, npargs, pa, eshellBASH);
    fprintf(out, "case \"$p\" in\n");

    pr_enums(out, npargs, pa, eshellBASH);
    pr_fopts(out, nfile, fnm, eshellBASH);
    fprintf(out, "esac }\n");
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

void write_completions(const gmx::CommandLineHelpContext &context,
                       int nfile,  t_filenm *fnm,
                       int npargs, t_pargs *pa)
{
    FILE    *out = context.writerContext().outputFile().handle();

    int      npar;
    t_pargs *par;

    // Remove hidden arguments.
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

    switch (context.writerContext().outputFormat())
    {
        case gmx::eHelpOutputFormat_CompletionBash:
            write_bashcompl(out, context.moduleDisplayName(), nfile, fnm, npar, par);
            break;
        case gmx::eHelpOutputFormat_CompletionCsh:
            write_cshcompl(out, nfile, fnm, npar, par);
            break;
        case gmx::eHelpOutputFormat_CompletionZsh:
            write_zshcompl(out, nfile, fnm, npar, par);
            break;
        default:
            GMX_THROW(gmx::NotImplementedError("Help format not implemented"));
    }

    sfree(par);
}
