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

#include <string>

#include <boost/scoped_ptr.hpp>

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/filenm.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "gromacs/legacyheaders/smalloc.h"

// TODO: Don't duplicate this from filenm.c/futil.c.
#define NZEXT 2
static const char *z_ext[NZEXT] = { ".gz", ".Z" };

static void pr_fopts(FILE *fp, int nf, const t_filenm tfn[])
{
    for (int i = 0; i < nf; i++)
    {
        const int   ftp          = tfn[i].ftp;
        const char *multiplicity = "(( $n == 1 )) && ";
        if (tfn[i].flag & ffMULT)
        {
            multiplicity = "";
        }
        if (ftp == efRND)
        {
            fprintf(fp, "%s) %sCOMPREPLY=( $(compgen -S ' ' -d $c) );;\n",
                    tfn[i].opt, multiplicity);
            continue;
        }
        fprintf(fp, "%s) %sCOMPREPLY=( $(compgen -S ' ' -X '!*.",
                tfn[i].opt, multiplicity);
        const int genericCount = ftp2generic_count(ftp);
        if (genericCount > 0)
        {
            fprintf(fp, "@(");
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
        fprintf(fp, "?(");
        for (int j = 0; j < NZEXT; j++)
        {
            if (j > 0)
            {
                fprintf(fp, "|");
            }
            fprintf(fp, "%s", z_ext[j]);
        }
        fprintf(fp, ")' -f $c ; compgen -S '/' -d $c ));;\n");
    }
}

static void pr_opts(FILE *fp,
                    int nfile,  t_filenm *fnm,
                    int npargs, t_pargs pa[])
{
    fprintf(fp, "if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen -S ' '  -W $'");
    const char *sep = "";
    for (int i = 0; i < nfile; i++)
    {
        fprintf(fp, "%s-%s", sep, fnm[i].opt+1);
        sep = "\\n";
    }
    for (int i = 0; i < npargs; i++)
    {
        if (pa[i].type == etBOOL && *(pa[i].u.b))
        {
            fprintf(fp, "%s-no%s", sep, pa[i].option + 1);
        }
        else
        {
            fprintf(fp, "%s-%s", sep, pa[i].option + 1);
        }
        sep = "\\n";
    }
    fprintf(fp, "' -- $c)); return 0; fi\n");
}

static void pr_enums(FILE *fp, int npargs, t_pargs pa[])
{
    for (int i = 0; i < npargs; i++)
    {
        if (pa[i].type == etENUM)
        {
            fprintf(fp, "%s) (( $n == 1 )) && COMPREPLY=( $(compgen -S ' ' -W $'", pa[i].option);
            for (int j = 1; pa[i].u.c[j]; j++)
            {
                fprintf(fp, "%s%s", (j == 1 ? "" : "\\n"), pa[i].u.c[j]);
            }
            fprintf(fp, "' -- $c ));;\n");
        }
    }
}

static void write_bashcompl(FILE *out,
                            const char *funcName,
                            int nfile,  t_filenm *fnm,
                            int npargs, t_pargs *pa)
{
    /* Advanced bash completions are handled by shell functions.
     * p and c hold the previous and current word on the command line.
     */
    fprintf(out, "%s() {\n", funcName);
    fprintf(out, "local IFS=$'\\n'\n");
    fprintf(out, "local c=${COMP_WORDS[COMP_CWORD]}\n");
    fprintf(out, "local n\n");
    fprintf(out, "for ((n=1;n<COMP_CWORD;++n)) ; do [[ \"${COMP_WORDS[COMP_CWORD-n]}\" == -* ]] && break ; done\n");
    fprintf(out, "local p=${COMP_WORDS[COMP_CWORD-n]}\n");
    fprintf(out, "COMPREPLY=()\n");

    pr_opts(out, nfile, fnm, npargs, pa);
    fprintf(out, "case \"$p\" in\n");

    pr_enums(out, npargs, pa);
    pr_fopts(out, nfile, fnm);
    fprintf(out, "esac }\n");
}

namespace gmx
{

class ShellCompletionWriter::Impl
{
    public:
        Impl(const std::string &binaryName, ShellCompletionFormat /*format*/)
            : binaryName_(binaryName)
        {
        }

        std::string completionFunctionName(const char *moduleName) const
        {
            // TODO: Consider if some characters need to be escaped.
            return formatString("_%s_%s_compl", binaryName_.c_str(), moduleName);
        }

        std::string             binaryName_;
        boost::scoped_ptr<File> file_;
};

ShellCompletionWriter::ShellCompletionWriter(const std::string     &binaryName,
                                             ShellCompletionFormat  format)
    : impl_(new Impl(binaryName, format))
{
}

ShellCompletionWriter::~ShellCompletionWriter()
{
}

File *ShellCompletionWriter::outputFile()
{
    return impl_->file_.get();
}

void ShellCompletionWriter::startCompletions()
{
    impl_->file_.reset(new File(impl_->binaryName_ + "-completion.bash", "w"));
    impl_->file_->writeLine("shopt -s extglob");
}

void ShellCompletionWriter::writeLegacyModuleCompletions(
        const char *moduleName,
        int nfile,  t_filenm *fnm,
        int npargs, t_pargs *pa)
{
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

    write_bashcompl(impl_->file_->handle(),
                    impl_->completionFunctionName(moduleName).c_str(),
                    nfile, fnm, npar, par);

    sfree(par);
}

void ShellCompletionWriter::writeModuleCompletions(
        const char    *moduleName,
        const Options  & /*options*/)
{
    // TODO: Implement.
    impl_->file_->writeLine(
            impl_->completionFunctionName(moduleName) + "() {\nCOMPREPLY=()\n}\n");
}

void ShellCompletionWriter::writeWrapperCompletions(
        const ModuleNameList &modules)
{
    impl_->file_->writeLine("_" + impl_->binaryName_ + "_compl() {");
    impl_->file_->writeLine("local i c m");
    impl_->file_->writeLine("local IFS=$'\\n'\n");
    impl_->file_->writeLine("COMPREPLY=()");
    impl_->file_->writeLine("unset COMP_WORDS[0]");
    impl_->file_->writeLine("for ((i=1;i<COMP_CWORD;++i)) ; do");
    impl_->file_->writeLine("[[ \"${COMP_WORDS[i]}\" != -* ]] && break");
    impl_->file_->writeLine("unset COMP_WORDS[i]");
    impl_->file_->writeLine("done");
    impl_->file_->writeLine("if (( i == COMP_CWORD )); then");
    impl_->file_->writeLine("c=${COMP_WORDS[COMP_CWORD]}");
    // TODO: Get rid of these hard-coded options.
    std::string completions("-h\\n-quiet\\n-version\\n-nocopyright");
    for (ModuleNameList::const_iterator i = modules.begin();
         i != modules.end(); ++i)
    {
        completions.append("\\n");
        completions.append(*i);
    }
    impl_->file_->writeLine("COMPREPLY=( $(compgen -S ' ' -W $'" + completions + "' -- $c) )");
    impl_->file_->writeLine("return 0");
    impl_->file_->writeLine("fi");
    impl_->file_->writeLine("m=${COMP_WORDS[i]}");
    impl_->file_->writeLine("COMP_WORDS=( \"${COMP_WORDS[@]}\" )");
    impl_->file_->writeLine("COMP_CWORD=$((COMP_CWORD-i))");
    impl_->file_->writeLine("case \"$m\" in");
    for (ModuleNameList::const_iterator i = modules.begin();
         i != modules.end(); ++i)
    {
        const char *const name = i->c_str();
        impl_->file_->writeLine(formatString("%s) %s ;;", name,
                                             impl_->completionFunctionName(name).c_str()));
    }
    impl_->file_->writeLine("esac }");
}

void ShellCompletionWriter::finishCompletions()
{
    impl_->file_->close();
    impl_->file_.reset();
}

} // namespace gmx
