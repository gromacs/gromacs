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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "viewit.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <array>
#include <filesystem>
#include <string>
#include <type_traits>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"

static constexpr std::array<int, 5> canViewFileType = { 0, efEPS, efXPM, efXVG, efPDB };

static constexpr int numberOfPossibleFiles = canViewFileType.size();

static int can_view(int ftp)
{
    for (int i = 1; i < numberOfPossibleFiles; i++)
    {
        if (ftp == canViewFileType[i])
        {
            return i;
        }
    }

    return 0;
}

void do_view(const gmx_output_env_t* oenv, const char* fn, const char* opts)
{
    std::array<const char*, 5> viewProgram = {
        nullptr, "ghostview", "display", nullptr, "xterm -e rasmol"
    };
    char        buf[STRLEN], env[STRLEN];
    const char* cmd;
    int         ftp, n;

    if (output_env_get_view(oenv) && fn)
    {
        if (getenv("DISPLAY") == nullptr)
        {
            fprintf(stderr, "Can not view %s, no DISPLAY environment variable.\n", fn);
        }
        else
        {
            ftp = fn2ftp(fn);
            sprintf(env, "GMX_VIEW_%s", ftp2ext(ftp));
            upstring(env);
            switch (ftp)
            {
                case efXVG:
                    if (!(cmd = getenv(env)))
                    {
                        if (getenv("GMX_USE_XMGR"))
                        {
                            cmd = "xmgr";
                        }
                        else
                        {
                            cmd = "xmgrace";
                        }
                    }
                    break;
                default:
                    if ((n = can_view(ftp)))
                    {
                        if (!(cmd = getenv(env)))
                        {
                            cmd = viewProgram[n];
                        }
                    }
                    else
                    {
                        fprintf(stderr, "Don't know how to view file %s", fn);
                        return;
                    }
            }
            if (strlen(cmd))
            {
                sprintf(buf, "%s %s %s &", cmd, opts ? opts : "", fn);
                fprintf(stderr, "Executing '%s'\n", buf);
                if (0 != system(buf))
                {
                    gmx_fatal(FARGS, "Failed executing command: %s", buf);
                }
            }
        }
    }
}

void view_all(const gmx_output_env_t* oenv, int nf, t_filenm fnm[])
{
    int i;

    for (i = 0; i < nf; i++)
    {
        if (can_view(fnm[i].ftp) && is_output(&(fnm[i])) && (!is_optional(&(fnm[i])) || is_set(&(fnm[i]))))
        {
            do_view(oenv, fnm[i].filenames[0].c_str(), nullptr);
        }
    }
}
