/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include "filenm.h"

#include <cctype>
#include <cstdio>
#include <cstring>

#include <filesystem>
#include <string_view>

#include "gromacs/fileio/filetypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

/* Use bitflag ... */
static bool IS_SET(const t_filenm& fileOption)
{
    return (fileOption.flag & ffSET) != 0;
}

static bool IS_OPT(const t_filenm& fileOption)
{
    return (fileOption.flag & ffOPT) != 0;
}

static const t_filenm* getFileOption(const char* opt, int nfile, const t_filenm fnm[])
{
    GMX_RELEASE_ASSERT(nfile == 0 || fnm, "need a valid list of filenames");

    for (int i = 0; i < nfile; i++)
    {
        if ((fnm[i].opt != nullptr && strcmp(opt, fnm[i].opt) == 0)
            || (fnm[i].opt == nullptr && strcmp(opt, ftp2defopt(fnm[i].ftp)) == 0))
        {
            return &fnm[i];
        }
    }

    return nullptr;
}

const char* opt2fn(const char* opt, int nfile, const t_filenm fnm[])
{
    const t_filenm* fileOption = getFileOption(opt, nfile, fnm);

    if (fileOption)
    {
        return fileOption->filenames[0].c_str();
    }

    GMX_RELEASE_ASSERT(false, "opt2fn should be called with a valid option");

    return nullptr;
}

gmx::ArrayRef<const std::string> opt2fns(const char* opt, int nfile, const t_filenm fnm[])
{
    const t_filenm* fileOption = getFileOption(opt, nfile, fnm);

    if (fileOption)
    {
        return fileOption->filenames;
    }

    GMX_RELEASE_ASSERT(false, "opt2fns should be called with a valid option");

    return {};
}

gmx::ArrayRef<const std::string> opt2fnsIfOptionSet(const char* opt, int nfile, const t_filenm fnm[])
{
    if (opt2bSet(opt, nfile, fnm))
    {
        return opt2fns(opt, nfile, fnm);
    }
    else
    {
        return {};
    }
}

const char* ftp2fn(int ftp, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
    {
        if (ftp == fnm[i].ftp)
        {
            return fnm[i].filenames[0].c_str();
        }
    }

    GMX_RELEASE_ASSERT(false, "ftp2fn should be called with a valid option");

    return nullptr;
}

gmx::ArrayRef<const std::string> ftp2fns(int ftp, int nfile, const t_filenm fnm[])
{
    for (int i = 0; (i < nfile); i++)
    {
        if (ftp == fnm[i].ftp)
        {
            return fnm[i].filenames;
        }
    }

    GMX_RELEASE_ASSERT(false, "ftp2fns should be called with a valid option");

    return {};
}

gmx_bool ftp2bSet(int ftp, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
    {
        if (ftp == fnm[i].ftp)
        {
            return static_cast<gmx_bool>(IS_SET(fnm[i]));
        }
    }

    GMX_RELEASE_ASSERT(false, "ftp2bSet should be called with a valid option");

    return FALSE;
}

gmx_bool opt2bSet(const char* opt, int nfile, const t_filenm fnm[])
{
    const t_filenm* fileOption = getFileOption(opt, nfile, fnm);

    if (fileOption)
    {
        return static_cast<gmx_bool>(IS_SET(*fileOption));
    }

    GMX_RELEASE_ASSERT(false, "opt2bSet should be called with a valid option");

    return FALSE;
}

std::optional<std::filesystem::path> opt2path_optional(const char* opt, int nfile, const t_filenm fnm[])
{
    const char* r = opt2fn_null(opt, nfile, fnm);
    if (r != nullptr)
    {
        return std::make_optional<std::filesystem::path>(r);
    }
    return std::nullopt;
}

const char* opt2fn_null(const char* opt, int nfile, const t_filenm fnm[])
{
    const t_filenm* fileOption = getFileOption(opt, nfile, fnm);

    if (fileOption)
    {
        if (IS_OPT(*fileOption) && !IS_SET(*fileOption))
        {
            return nullptr;
        }
        else
        {
            return fileOption->filenames[0].c_str();
        }
    }

    GMX_RELEASE_ASSERT(false, "opt2fn_null should be called with a valid option");

    return nullptr;
}

std::optional<std::filesystem::path> ftp2path_optional(int ftp, int nfile, const t_filenm fnm[])
{
    const char* r = ftp2fn_null(ftp, nfile, fnm);
    if (r != nullptr)
    {
        return std::make_optional<std::filesystem::path>(r);
    }
    return std::nullopt;
}

const char* ftp2fn_null(int ftp, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
    {
        if (ftp == fnm[i].ftp)
        {
            if (IS_OPT(fnm[i]) && !IS_SET(fnm[i]))
            {
                return nullptr;
            }
            else
            {
                return fnm[i].filenames[0].c_str();
            }
        }
    }

    GMX_RELEASE_ASSERT(false, "ftp2fn_null should be called with a valid option");

    return nullptr;
}

gmx_bool is_optional(const t_filenm* fnm)
{
    return ((fnm->flag & ffOPT) == ffOPT);
}

gmx_bool is_output(const t_filenm* fnm)
{
    return ((fnm->flag & ffWRITE) == ffWRITE);
}

gmx_bool is_set(const t_filenm* fnm)
{
    return ((fnm->flag & ffSET) == ffSET);
}

namespace
{

/*! \brief Return the first position within \c filename of the ".partNNNN"
 * interior sequence produced by mdrun -noappend, or npos if not found. */
size_t findSuffixFromNoAppendPosition(const std::string_view filename)
{
    size_t partPosition = filename.find(".part");
    if ((partPosition != decltype(filename)::npos) && (filename.length() - partPosition >= 10)
        && (std::isdigit(filename[partPosition + 5])) && (std::isdigit(filename[partPosition + 6]))
        && (std::isdigit(filename[partPosition + 7])) && (std::isdigit(filename[partPosition + 8]))
        && filename[partPosition + 9] == '.')
    {
        return partPosition;
    }
    return decltype(filename)::npos;
}

} // namespace

bool hasSuffixFromNoAppend(const std::string_view filename)
{
    return (findSuffixFromNoAppendPosition(filename) != decltype(filename)::npos);
}

int add_suffix_to_output_names(t_filenm* fnm, int nfile, const char* suffix)
{
    for (int i = 0; i < nfile; i++)
    {
        if (is_output(&fnm[i]) && fnm[i].ftp != efCPT)
        {
            /* We never use multiple _outputs_, but we might as well check
               for it, just in case... */
            for (std::string& filename : fnm[i].filenames)
            {
                // mdrun should not generate files like
                // md.part0002.part0003.log. mdrun should permit users
                // to name files like md.equil.part0002.log. So,
                // before we use concatenateBeforeExtension to
                // add the requested suffix, we need to check for
                // files matching mdrun's pattern for adding part
                // numbers. Then we can remove that if needed.
                for (size_t partPosition;
                     (partPosition = findSuffixFromNoAppendPosition(filename)) != std::string::npos;)
                {
                    // Remove the ".partNNNN" that we have found,
                    // and then run the loop again to make sure
                    // there isn't another one to remove, somehow.
                    std::string temporary = filename.substr(0, partPosition);
                    temporary += filename.substr(partPosition + 9);
                    filename.swap(temporary);
                }
                filename = gmx::concatenateBeforeExtension(filename, suffix).string();
            }
        }
    }
    return 0;
}
