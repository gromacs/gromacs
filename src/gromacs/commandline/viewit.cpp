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
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/path.h"

//! A program that can display a file, and how to tell whether it is installed.
struct Viewer
{
    //! Command line that views the file.
    const char* command;
    /*! \brief Executable whose presence decides whether \c command can run.
     *
     * Defaults to the first word of \c command, which is the right test unless
     * the command wraps the real viewer in another program, as
     * "xterm -e rasmol" does.
     */
    const char* program = nullptr;
};

/*! \name Viewers for each file type, most preferred first.
 *
 * The viewer GROMACS has always used for a type heads its list, so that anyone
 * who has it keeps getting it. The alternatives only take effect on machines
 * where the established viewer is absent and \c -w therefore used to fail.
 * \{
 */
static constexpr std::array<Viewer, 5> c_epsViewers = {
    { { "ghostview" }, { "gv" }, { "evince" }, { "okular" }, { "xdg-open" } }
};
static constexpr std::array<Viewer, 4> c_xpmViewers = {
    { { "display" }, { "eog" }, { "feh" }, { "xdg-open" } },
};
//! GMX_USE_XMGR still forces xmgr ahead of all of these; see do_view() below.
static constexpr std::array<Viewer, 4> c_xvgViewers = {
    { { "xmgrace" }, { "grace" }, { "qtgrace" }, { "xmgr" } }
};
//! RasMol is driven from a terminal, so test for RasMol, not the xterm wrapping it.
static constexpr std::array<Viewer, 5> c_pdbViewers = {
    { { "xterm -e rasmol", "rasmol" }, { "pymol" }, { "chimerax" }, { "jmol" }, { "vmd" } }
};
static constexpr std::array<Viewer, 3> c_cifViewers = {
    { { "pymol" }, { "chimerax" }, { "jmol" } },
};
static constexpr std::array<Viewer, 4> c_groViewers = {
    { { "pymol" }, { "chimerax" }, { "jmol" }, { "vmd" } },
};
//! \}

/*! \brief A file type do_view() can display, and the viewers that display it.
 *
 * Pairing the file type with its viewers means the two cannot fall out of step,
 * as separate parallel arrays indexed by the same number did.
 */
struct ViewerForFileType
{
    //! The file type, as enumerated in filetypes.h.
    int fileType;
    //! The viewers that can display it, most preferred first.
    gmx::ArrayRef<const Viewer> viewers;
};

//! The file types that do_view() knows how to display.
static const std::array<ViewerForFileType, 6> c_viewerForFileType = {
    { { efEPS, c_epsViewers },
      { efXPM, c_xpmViewers },
      { efXVG, c_xvgViewers },
      { efPDB, c_pdbViewers },
      { efCIF, c_cifViewers },
      { efGRO, c_groViewers } },
};

//! Whether the program \p viewer needs is installed.
static bool viewerIsAvailable(const Viewer& viewer)
{
    if (viewer.command == nullptr)
    {
        return false;
    }

    std::string program;
    if (viewer.program != nullptr)
    {
        program = viewer.program;
    }
    else
    {
        program                 = viewer.command;
        const size_t firstSpace = program.find(' ');
        if (firstSpace != std::string::npos)
        {
            program.resize(firstSpace);
        }
    }

    std::error_code errorCode;
    if (program.find('/') != std::string::npos)
    {
        return std::filesystem::is_regular_file(program, errorCode);
    }
    for (const std::filesystem::path& directory : gmx::getSystemExecutablePaths())
    {
        if (std::filesystem::is_regular_file(directory / program, errorCode))
        {
            return true;
        }
    }

    return false;
}

/*! \brief Return how to view \p ftp, or nullptr if there is no way to view it.
 *
 * Prefers a viewer that is actually installed. If none of them is, returns the
 * most preferred one anyway, so that the user is told about a viewer we do
 * suggest rather than nothing at all.
 */
static const Viewer* findViewer(int ftp)
{
    for (const ViewerForFileType& entry : c_viewerForFileType)
    {
        if (entry.fileType != ftp)
        {
            continue;
        }
        for (const Viewer& viewer : entry.viewers)
        {
            if (viewerIsAvailable(viewer))
            {
                return &viewer;
            }
        }
        return entry.viewers.empty() ? nullptr : &entry.viewers.front();
    }

    return nullptr;
}

static bool can_view(int ftp)
{
    return findViewer(ftp) != nullptr;
}

void do_view(const gmx_output_env_t* oenv, const char* fn, const char* opts)
{
    char        buf[STRLEN], env[STRLEN];
    const char* cmd;
    int         ftp;

    if (output_env_get_view(oenv) && fn)
    {
        if (std::getenv("DISPLAY") == nullptr)
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
                    if (!(cmd = std::getenv(env)))
                    {
                        if (std::getenv("GMX_USE_XMGR"))
                        {
                            cmd = "xmgr";
                        }
                        else
                        {
                            // xmgrace heads the efXVG entries, so this still
                            // selects it whenever it is installed.
                            const Viewer* viewer = findViewer(ftp);
                            cmd = (viewer != nullptr) ? viewer->command : "xmgrace";
                        }
                    }
                    break;
                default:
                {
                    const Viewer* viewer = findViewer(ftp);
                    if (viewer == nullptr)
                    {
                        fprintf(stderr, "Don't know how to view file %s", fn);
                        return;
                    }
                    if (!(cmd = std::getenv(env)))
                    {
                        cmd = viewer->command;
                    }
                    break;
                }
            }
            if (cmd != nullptr && std::strlen(cmd))
            {
                sprintf(buf, "%s %s %s &", cmd, opts ? opts : "", fn);
                fprintf(stderr, "Executing '%s'\n", buf);
                if (0 != std::system(buf))
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
