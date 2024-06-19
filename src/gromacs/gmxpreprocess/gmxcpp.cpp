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

#include "gmxcpp.h"

#include <cctype>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <memory>
#include <unordered_set>
#include <vector>

#include <sys/types.h>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"

struct t_define
{
    std::string name;
    std::string def;
};

/* enum used for handling ifdefs */
enum
{
    eifTRUE,
    eifFALSE,
    eifIGNORE,
    eifNR
};

struct gmx_cpp
{
    std::shared_ptr<std::vector<t_define>>              defines;
    std::shared_ptr<std::vector<std::filesystem::path>> includes;
    std::unordered_set<std::string>                     unmatched_defines;
    FILE*                                               fp = nullptr;
    std::filesystem::path                               path;
    std::filesystem::path                               cwd;
    std::filesystem::path                               fn;
    std::string                                         line;
    int                                                 line_nr;
    std::vector<int>                                    ifdefs;
    struct gmx_cpp*                                     child  = nullptr;
    struct gmx_cpp*                                     parent = nullptr;
};

static bool is_word_end(char c)
{
    return !((isalnum(c) != 0) || c == '_');
}

static const char* strstrw(const char* buf, const char* word)
{
    const char* ptr;

    while ((ptr = strstr(buf, word)) != nullptr)
    {
        /* Check if we did not find part of a longer word */
        if (ptr && is_word_end(ptr[strlen(word)])
            && (((ptr > buf) && is_word_end(ptr[-1])) || (ptr == buf)))
        {
            return ptr;
        }

        buf = ptr + strlen(word);
    }
    return nullptr;
}

/* Finds a preprocessor directive, whose name (after the '#') is
 * returned in *name, and the remainder of the line after leading
 * whitespace, without trailing whitespace, is returned in *val
 */
static bool find_directive(const char* buf, std::string* name, std::string* val)
{
    /* Skip initial whitespace */
    while (isspace(*buf))
    {
        ++buf;
    }
    /* Check if this is a directive */
    if (*buf != '#')
    {
        return FALSE;
    }
    /* Skip the hash and any space after it */
    ++buf;
    while (isspace(*buf))
    {
        ++buf;
    }
    /* Set the name pointer and find the next space */
    name->clear();
    while (*buf != '\0' && !isspace(*buf))
    {
        *name += *buf;
        ++buf;
    }
    /* Set the end of the name here, and skip any space */
    if (*buf != '\0')
    {
        ++buf;
        while (isspace(*buf))
        {
            ++buf;
        }
    }
    /* Check if anything is remaining */
    if (*buf != '\0')
    {
        *val = buf;
        // Remove trailing whitespace
        while (!val->empty() && isspace(val->back()))
        {
            val->resize(val->size() - 1);
        }
    }
    else
    {
        val->clear();
    }

    return TRUE;
}

static bool is_ifdeffed_out(gmx::ArrayRef<const int> ifdefs)
{
    return (!ifdefs.empty() && ifdefs.back() != eifTRUE);
}

static void add_include(std::vector<std::filesystem::path>* includes, const std::filesystem::path& includePath)
{
    GMX_RELEASE_ASSERT(includes, "Need valid includes");
    GMX_RELEASE_ASSERT(!includePath.empty(), "Need a valid include path");

    for (const auto& include : *includes)
    {
        if (include.string() == includePath.string())
        {
            return;
        }
    }

    includes->emplace_back(includePath);
}

static void add_define(std::vector<t_define>* defines, const std::string& name, const char* value)
{
    GMX_RELEASE_ASSERT(defines, "Need defines");
    GMX_RELEASE_ASSERT(value, "Need a value");

    for (t_define& define : *defines)
    {
        if (define.name == name)
        {
            define.def = value;
            return;
        }
    }

    defines->push_back({ name, value });
}

/* Open the file to be processed. The handle variable holds internal
   info for the cpp emulator. Return integer status */
static int cpp_open_file(const std::filesystem::path&                         filenm,
                         gmx_cpp_t*                                           handle,
                         char**                                               cppopts,
                         std::shared_ptr<std::vector<t_define>>*              definesFromParent,
                         std::shared_ptr<std::vector<std::filesystem::path>>* includesFromParent)
{
    // TODO: We should avoid new/delete, we should use Pimpl instead
    gmx_cpp* cpp = new gmx_cpp;
    *handle      = cpp;

    if (definesFromParent)
    {
        cpp->defines = *definesFromParent;
    }
    else
    {
        cpp->defines = std::make_shared<std::vector<t_define>>();
    }

    if (includesFromParent)
    {
        cpp->includes = *includesFromParent;
    }
    else
    {
        cpp->includes = std::make_shared<std::vector<std::filesystem::path>>();
    }

    /* First process options, they might be necessary for opening files
       (especially include statements). */
    int i = 0;
    if (cppopts)
    {
        while (cppopts[i])
        {
            if (strstr(cppopts[i], "-I") == cppopts[i])
            {
                add_include(cpp->includes.get(), cppopts[i] + 2);
            }
            if (strstr(cppopts[i], "-D") == cppopts[i])
            {
                /* If the option contains a =, split it into name and value. */
                char* ptr = strchr(cppopts[i], '=');
                if (ptr)
                {
                    std::string buf = cppopts[i] + 2;
                    buf.resize(ptr - cppopts[i] - 2);
                    add_define(cpp->defines.get(), buf, ptr + 1);
                    cpp->unmatched_defines.insert(buf);
                }
                else
                {
                    add_define(cpp->defines.get(), cppopts[i] + 2, "");
                    cpp->unmatched_defines.insert(cppopts[i] + 2);
                }
            }
            i++;
        }
    }

    /* Find the file. First check whether it is in the current directory. */
    if (gmx_fexist(filenm))
    {
        cpp->fn = filenm;
    }
    else
    {
        /* If not, check all the paths given with -I. */
        for (const auto& include : *cpp->includes)
        {
            auto buf = include;
            buf.append(filenm.string());
            if (gmx_fexist(buf))
            {
                cpp->fn = buf;
                break;
            }
        }
        /* If still not found, check the Gromacs library search path. */
        if (cpp->fn.empty())
        {
            cpp->fn = gmx::findLibraryFile(filenm, false, false);
        }
    }
    if (cpp->fn.empty())
    {
        gmx_fatal(FARGS, "Topology include file \"%s\" not found", filenm.string().c_str());
    }
    /* If the file name has a path component, we need to change to that
     * directory.
     */
    if (cpp->fn.has_parent_path())
    {
        cpp->path = cpp->fn.parent_path();
        cpp->fn   = cpp->fn.filename();
        cpp->cwd  = gmx_getcwd();
        gmx_chdir(cpp->path);
    }
    cpp->line.clear();
    cpp->line_nr = 0;
    cpp->ifdefs.clear();
    cpp->child  = nullptr;
    cpp->parent = nullptr;
    if (cpp->fp == nullptr)
    {
        cpp->fp = fopen(cpp->fn.string().c_str(), "r");
    }
    if (cpp->fp == nullptr)
    {
        switch (errno)
        {
            case EINVAL:
            default: return eCPP_UNKNOWN;
        }
    }
    return eCPP_OK;
}

/* Open the file to be processed. The handle variable holds internal
   info for the cpp emulator. Return integer status */
int cpp_open_file(const std::filesystem::path& filenm, gmx_cpp_t* handle, char** cppopts)
{
    return cpp_open_file(filenm, handle, cppopts, nullptr, nullptr);
}

/* Note that dval might be null, e.g. when handling a line like '#define */
static int process_directive(gmx_cpp_t* handlep, const std::string& dname, const std::string& dval)
{
    gmx_cpp_t handle = *handlep;

    std::vector<int>& ifdefs = handle->ifdefs;

    /* #ifdef or ifndef statement */
    bool bIfdef  = (dname == "ifdef");
    bool bIfndef = (dname == "ifndef");
    if (bIfdef || bIfndef)
    {
        if (is_ifdeffed_out(ifdefs))
        {
            handle->ifdefs.push_back(eifIGNORE);
        }
        else
        {
            // A bare '#ifdef' or '#ifndef' is invalid
            if (dval.empty())
            {
                return eCPP_SYNTAX;
            }
            bool found = false;
            for (const t_define& define : *handle->defines)
            {
                if (define.name == dval)
                {
                    // erase from unmatched_defines in original handle
                    gmx_cpp_t root = handle;
                    while (root->parent != nullptr)
                    {
                        root = root->parent;
                    }
                    root->unmatched_defines.erase(dval);

                    found = true;
                    break;
                }
            }
            if ((bIfdef && found) || (bIfndef && !found))
            {
                ifdefs.push_back(eifTRUE);
            }
            else
            {
                ifdefs.push_back(eifFALSE);
            }
        }
        return eCPP_OK;
    }

    /* #else statement */
    if (dname == "else")
    {
        if (ifdefs.empty())
        {
            return eCPP_SYNTAX;
        }
        if (ifdefs.back() == eifTRUE)
        {
            ifdefs.back() = eifFALSE;
        }
        else if (ifdefs.back() == eifFALSE)
        {
            ifdefs.back() = eifTRUE;
        }
        return eCPP_OK;
    }

    /* #endif statement */
    if (dname == "endif")
    {
        if (ifdefs.empty())
        {
            return eCPP_SYNTAX;
        }
        ifdefs.erase(ifdefs.end() - 1);
        return eCPP_OK;
    }

    /* Check whether we're not ifdeffed out. The order of this statement
       is important. It has to come after #ifdef, #else and #endif, but
       anything else should be ignored. */
    if (is_ifdeffed_out(ifdefs))
    {
        return eCPP_OK;
    }

    /* Check for include statements */
    if (dname == "include")
    {
        int len = -1;
        int i0  = 0;
        // A bare '#include' is an invalid line
        if (dval.empty())
        {
            return eCPP_SYNTAX;
        }
        // An include needs to be followed by either a '"' or a '<' as a first character.
        if ((dval[0] != '"') && (dval[0] != '<'))
        {
            return eCPP_INVALID_INCLUDE_DELIMITER;
        }
        for (size_t i1 = 0; i1 < dval.size(); i1++)
        {
            if ((dval[i1] == '"') || (dval[i1] == '<') || (dval[i1] == '>'))
            {
                if (len == -1)
                {
                    i0  = i1 + 1;
                    len = 0;
                }
                else
                {
                    break;
                }
            }
            else if (len >= 0)
            {
                len++;
            }
        }
        if (len == -1)
        {
            return eCPP_SYNTAX;
        }
        std::string inc_fn = dval.substr(i0, len);

        /* Open include file and store it as a child in the handle structure */
        int status = cpp_open_file(
                inc_fn.c_str(), &(handle->child), nullptr, &handle->defines, &handle->includes);
        if (status != eCPP_OK)
        {
            handle->child = nullptr;
            return status;
        }
        /* Make a linked list of open files and move on to the include file */
        handle->child->parent = handle;
        *handlep              = handle->child;
        return eCPP_OK;
    }

    /* #define statement */
    if (dname == "define")
    {
        // A bare '#define' is an invalid line
        if (dval.empty())
        {
            return eCPP_SYNTAX;
        }
        /* Split it into name and value. */
        const char* ptr = dval.c_str();
        while ((*ptr != '\0') && !isspace(*ptr))
        {
            ptr++;
        }
        std::string name = dval.substr(0, ptr - dval.c_str());

        while ((*ptr != '\0') && isspace(*ptr))
        {
            ptr++;
        }

        add_define(handle->defines.get(), name, ptr);
        return eCPP_OK;
    }

    /* #undef statement */
    if (dname == "undef")
    {
        // A bare '#undef' is an invalid line
        if (dval.empty())
        {
            return eCPP_SYNTAX;
        }
        std::vector<t_define>& defines = *handle->defines;
        for (size_t i = 0; i < defines.size(); i++)
        {
            if (defines[i].name == dval)
            {
                defines.erase(defines.begin() + i);
                break;
            }
        }

        return eCPP_OK;
    }

    /* If we haven't matched anything, this is an unknown directive */
    return eCPP_SYNTAX;
}

/* Return one whole line from the file into buf which holds at most n
   characters, for subsequent processing. Returns integer status. This
   routine also does all the "intelligent" work like processing cpp
   directives and so on. Note that often the routine is called
   recursively and no cpp directives are printed. */
int cpp_read_line(gmx_cpp_t* handlep, int n, char buf[])
{
    gmx_cpp_t handle = *handlep;
    int       status;
    bool      bEOF;

    if (!handle)
    {
        return eCPP_INVALID_HANDLE;
    }
    if (!handle->fp)
    {
        return eCPP_FILE_NOT_OPEN;
    }

    bEOF = (feof(handle->fp) != 0);
    if (!bEOF)
    {
        /* Read the actual line now. */
        if (fgets2(buf, n - 1, handle->fp) == nullptr)
        {
            /* Recheck EOF, since we could have been at the end before
             * the fgets2 call, but we need to read past the end to know.
             */
            bEOF = (feof(handle->fp) != 0);
            if (!bEOF)
            {
                /* Something strange happened, fgets returned NULL,
                 * but we are not at EOF. Maybe wrong line endings?
                 */
                return eCPP_UNKNOWN;
            }
        }
    }

    if (bEOF)
    {
        if (handle->parent == nullptr)
        {
            return eCPP_EOF;
        }
        cpp_close_file(handlep);
        *handlep = handle->parent;
        delete handle;
        return cpp_read_line(handlep, n, buf);
    }
    else
    {
        handle->line = buf;
        handle->line_nr++;
    } /* Now we've read a line! */

    /* Process directives if this line contains one */
    std::string dname;
    std::string dval;
    if (find_directive(buf, &dname, &dval))
    {
        status = process_directive(handlep, dname, dval);
        if (status != eCPP_OK)
        {
            return status;
        }
        /* Don't print lines with directives, go on to the next */
        return cpp_read_line(handlep, n, buf);
    }

    /* Check whether we're not ifdeffed out. The order of this statement
       is important. It has to come after #ifdef, #else and #endif, but
       anything else should be ignored. */
    if (is_ifdeffed_out(handle->ifdefs))
    {
        return cpp_read_line(handlep, n, buf);
    }

    /* Check whether we have any defines that need to be replaced. Note
       that we have to use a best fit algorithm, rather than first come
       first go. We do this by sorting the defines on length first, and
       then on alphabetical order. */
    for (t_define& define : *handle->defines)
    {
        if (!define.def.empty())
        {
            int         nn  = 0;
            const char* ptr = buf;
            while ((ptr = strstrw(ptr, define.name.c_str())) != nullptr)
            {
                nn++;
                ptr += strlen(define.name.c_str());
            }
            if (nn > 0)
            {
                // Need to erase  unmatched define in original handle
                gmx_cpp_t root = handle;
                while (root->parent != nullptr)
                {
                    root = root->parent;
                }
                root->unmatched_defines.erase(define.name);

                std::string name;
                const char* ptr = buf;
                const char* ptr2;
                while ((ptr2 = strstrw(ptr, define.name.c_str())) != nullptr)
                {
                    name.append(ptr, ptr2 - ptr);
                    name += define.def;
                    ptr = ptr2 + define.name.size();
                }
                name += ptr;
                GMX_RELEASE_ASSERT(name.size() < static_cast<size_t>(n),
                                   "The line should fit in buf");
                strcpy(buf, name.c_str());
            }
        }
    }

    return eCPP_OK;
}

std::filesystem::path cpp_cur_file(const gmx_cpp_t* handlep)
{
    return (*handlep)->fn;
}

int cpp_cur_linenr(const gmx_cpp_t* handlep)
{
    return (*handlep)->line_nr;
}

/* Close the file! Return integer status. */
int cpp_close_file(gmx_cpp_t* handlep)
{
    gmx_cpp_t handle = *handlep;

    if (!handle)
    {
        return eCPP_INVALID_HANDLE;
    }
    if (!handle->fp)
    {
        return eCPP_FILE_NOT_OPEN;
    }
    fclose(handle->fp);

    if (!handle->cwd.empty())
    {
        gmx_chdir(handle->cwd);
    }

    handle->fp      = nullptr;
    handle->line_nr = 0;
    handle->line.clear();

    return eCPP_OK;
}

const std::string* cpp_find_define(const gmx_cpp_t* handlep, const std::string& defineName)
{
    for (const t_define& define : *(*handlep)->defines)
    {
        if (define.name == defineName)
        {
            return &define.def;
        }
    }

    return nullptr;
}

void cpp_done(gmx_cpp_t handle)
{
    int status = cpp_close_file(&handle);
    if (status != eCPP_OK)
    {
        gmx_fatal(FARGS, "%s", cpp_error(&handle, status));
    }
    delete handle;
}

/* Return a string containing the error message coresponding to status
   variable */
char* cpp_error(gmx_cpp_t* handlep, int status)
{
    char        buf[256];
    const char* ecpp[] = { "OK",
                           "File not found",
                           "End of file",
                           "Syntax error",
                           "Interrupted",
                           "Invalid file handle",
                           "Invalid delimiter for filename in #include statement",
                           "File not open",
                           "Unknown error, perhaps your text file uses wrong line endings?",
                           "Error status out of range" };
    gmx_cpp_t   handle = *handlep;

    if (!handle)
    {
        return const_cast<char*>(ecpp[eCPP_INVALID_HANDLE]);
    }

    if ((status < 0) || (status >= eCPP_NR))
    {
        status = eCPP_NR;
    }

    sprintf(buf,
            "%s - File %s, line %d\nLast line read:\n'%s'",
            ecpp[status],
            (handle && !handle->fn.empty()) ? handle->fn.string().c_str() : "unknown",
            (handle) ? handle->line_nr : -1,
            !handle->line.empty() ? handle->line.c_str() : "");

    return gmx_strdup(buf);
}

std::string checkAndWarnForUnusedDefines(const gmx_cpp& handle)
{
    std::string warning;
    if (!handle.unmatched_defines.empty())
    {
        warning =
                "The following macros were defined in the 'define' mdp field with the -D prefix, "
                "but "
                "were not used in the topology:\n";
        for (const auto& str : handle.unmatched_defines)
        {
            warning += ("    " + str + "\n");
        }
        warning +=
                "If you haven't made a spelling error, either use the macro you defined, "
                "or don't define the macro";
    }
    return warning;
}
