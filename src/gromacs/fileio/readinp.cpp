/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "readinp.h"

#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "gromacs/fileio/warninp.h"
#include "gromacs/utility/binaryinformation.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/niceheader.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"
#include "gromacs/utility/textwriter.h"

std::vector<t_inpfile>
read_inpfile(gmx::TextInputStream *stream, const char *fn,
             warninp_t wi)
{
    std::vector<t_inpfile> inp;

    if (debug)
    {
        fprintf(debug, "Reading MDP file %s\n", fn);
    }

    int             indexOfLineReadFromFile = 0;
    std::string     line;
    gmx::TextReader reader(stream);
    reader.setTrimTrailingWhiteSpace(true);
    reader.setTrimTrailingComment(true, ';');
    while (reader.readLine(&line))
    {
        indexOfLineReadFromFile++;
        set_warning_line(wi, fn, indexOfLineReadFromFile);

        if (line.empty())
        {
            continue;
        }

        auto tokens = gmx::splitAndTrimDelimitedString(line, '=');
        if (tokens.size() < 2)
        {
            // TODO this seems like it silently ignores the user accidentally deleting an equals sign...
            if (debug)
            {
                fprintf(debug, "No = on line %d in file %s, ignored\n", indexOfLineReadFromFile, fn);
            }
            continue;
        }
        if (tokens.size() > 2)
        {
            // More than one equals symbol in the original line is
            // valid if the RHS is a free string, and needed for
            // "define = -DBOOLVAR -DVAR=VALUE".
            //
            // First, drop all the fields on the RHS of the first equals symbol.
            tokens.resize(1);
            // This find cannot return std::string::npos.
            auto firstEqualsPos = line.find('=');
            tokens.emplace_back(gmx::stripString(line.substr(firstEqualsPos + 1)));
        }
        if (tokens[0].empty())
        {
            // TODO ignoring such lines does not seem like good behaviour
            if (debug)
            {
                fprintf(debug, "Empty left hand side on line %d in file %s, ignored\n", indexOfLineReadFromFile, fn);
            }
            continue;
        }
        if (tokens[1].empty())
        {
            // TODO ignoring such lines does not seem like good behaviour
            if (debug)
            {
                fprintf(debug, "Empty right hand side on line %d in file %s, ignored\n", indexOfLineReadFromFile, fn);
            }
            continue;
        }

        /* Now finally something sensible; check for duplicates */
        int found_index = search_einp(inp, tokens[0].c_str());

        if (found_index == -1)
        {
            /* add a new item */
            inp.emplace_back(0, 1, false, false, false,
                             tokens[0], tokens[1]);
        }
        else
        {
            auto message = gmx::formatString("Parameter \"%s\" doubly defined\n",
                                             tokens[0].c_str());
            warning_error(wi, message.c_str());
        }
    }
    /* This preserves the behaviour of the old code, which issues some
       warnings after completing parsing. Regenerating regressiontest
       warning files is not worth the effort. */
    indexOfLineReadFromFile++;
    set_warning_line(wi, fn, indexOfLineReadFromFile);

    if (debug)
    {
        fprintf(debug, "Done reading MDP file, there were %zu entries in there\n",
                inp.size());
    }

    return inp;
}

gmx::KeyValueTreeObject flatKeyValueTreeFromInpFile(gmx::ArrayRef<const t_inpfile> inp)
{
    gmx::KeyValueTreeBuilder  builder;
    auto                      root = builder.rootObject();
    for (auto &local : inp)
    {
        root.addValue<std::string>(local.name_, !local.value_.empty() ? local.value_ : "");
    }
    return builder.build();
}


struct inp_comp
{
    bool operator()(t_inpfile const &a, t_inpfile const &b)
    {
        return a.count_ < b.count_;
    }
};

static void sort_inp(std::vector<t_inpfile> *inp)
{
    std::vector<t_inpfile> &inpRef = *inp;
    int                     mm;

    mm = -1;
    for (const auto &local : inpRef)
    {
        mm = std::max(mm, local.count_);
    }
    for (auto &local : inpRef)
    {
        if (local.count_ == 0)
        {
            local.count_ = mm++;
        }
    }
    std::sort(inpRef.begin(), inpRef.end(), inp_comp());
}

void write_inpfile(gmx::TextOutputStream *stream, const char *fn, std::vector<t_inpfile> *inp,
                   gmx_bool bHaltOnUnknown,
                   WriteMdpHeader writeHeader,
                   warninp_t wi)
{
    using gmx::formatString;

    sort_inp(inp);

    gmx::TextWriter writer(stream);
    if (writeHeader == WriteMdpHeader::yes)
    {
        gmx::niceHeader(&writer, fn, ';');

        gmx::BinaryInformationSettings settings;
        settings.generatedByHeader(true);
        settings.linePrefix(";\t");
        gmx::printBinaryInformation(&writer, gmx::getProgramContext(), settings);
    }
    for (const auto &local : *inp)
    {
        if (local.bHandledAsKeyValueTree_)
        {
        }
        else if (local.bSet_)
        {
            if (local.name_[0] == ';' || (local.name_.length() > 2 && local.name_[1] == ';'))
            {
                writer.writeLine(formatString("%-24s", local.name_.c_str()));
            }
            else
            {
                writer.writeLine(formatString("%-24s = %s", local.name_.c_str(), !local.value_.empty() ? local.value_.c_str() : ""));
            }
        }
        else if (!local.bObsolete_)
        {
            auto message = formatString("Unknown left-hand '%s' in parameter file\n",
                                        local.name_.c_str());
            if (bHaltOnUnknown)
            {
                warning_error(wi, message.c_str());
            }
            else
            {
                warning(wi, message.c_str());
            }
        }
    }

    check_warning_error(wi, FARGS);
}

void replace_inp_entry(gmx::ArrayRef<t_inpfile> inp, const char *old_entry, const char *new_entry)
{
    for (auto &local : inp)
    {
        if (gmx_strcasecmp_min(old_entry, local.name_.c_str()) == 0)
        {
            if (new_entry)
            {
                fprintf(stderr, "Replacing old mdp entry '%s' by '%s'\n",
                        local.name_.c_str(), new_entry);

                int foundIndex = search_einp(inp, new_entry);
                if (foundIndex >= 0)
                {
                    gmx_fatal(FARGS, "A parameter is present with both the old name '%s' and the new name '%s'.", local.name_.c_str(), inp[foundIndex].name_.c_str());
                }

                local.name_.assign(new_entry);
            }
            else
            {
                fprintf(stderr, "Ignoring obsolete mdp entry '%s'\n",
                        local.name_.c_str());
                local.bObsolete_ = TRUE;
            }
        }
    }
}

int search_einp(gmx::ArrayRef<const t_inpfile> inp, const char *name)
{
    if (inp.empty())
    {
        return -1;
    }
    for (gmx::index i = 0; i < inp.size(); i++)
    {
        if (gmx_strcasecmp_min(name, inp[i].name_.c_str()) == 0)
        {
            return i;
        }
    }
    return -1;
}

void mark_einp_set(gmx::ArrayRef<t_inpfile> inp, const char *name)
{
    int                     i      = search_einp(inp, name);
    if (i != -1)
    {
        inp[i].count_ = inp.front().inp_count_++;
        inp[i].bSet_  = TRUE;
        /* Prevent mdp lines being written twice for
           options that are handled via key-value trees. */
        inp[i].bHandledAsKeyValueTree_ = TRUE;
    }
}

static int get_einp(std::vector<t_inpfile> *inp, const char *name)
{
    std::vector<t_inpfile> &inpRef   = *inp;
    bool                    notfound = false;

    int                     i = search_einp(inpRef, name);

    if (i == -1)
    {
        notfound = true;
        inpRef.emplace_back(0, 0, false, true, false,
                            name, "");
        i = inpRef.size() - 1;
        if (inpRef.size()  == 1)
        {
            inpRef.front().inp_count_ = 1;
        }
    }
    inpRef[i].count_ = inpRef.front().inp_count_++;
    inpRef[i].bSet_  = TRUE;
    if (debug)
    {
        fprintf(debug, "Inp %d = %s\n", inpRef[i].count_, inpRef[i].name_.c_str());
    }

    if (notfound)
    {
        return -1;
    }
    else
    {
        return i;
    }
}

/* Note that sanitizing the trailing part of inp[ii].value was the responsibility of read_inpfile() */
int get_eint(std::vector<t_inpfile> *inp, const char *name, int def,
             warninp_t wi)
{
    std::vector<t_inpfile> &inpRef = *inp;
    char                    buf[32], *ptr, warn_buf[STRLEN];

    int                     ii = get_einp(inp, name);

    if (ii == -1)
    {
        sprintf(buf, "%d", def);
        inpRef.back().value_.assign(buf);

        return def;
    }
    else
    {
        int ret = std::strtol(inpRef[ii].value_.c_str(), &ptr, 10);
        if (*ptr != '\0')
        {
            sprintf(warn_buf, "Right hand side '%s' for parameter '%s' in parameter file is not an integer value\n", inpRef[ii].value_.c_str(), inpRef[ii].name_.c_str());
            warning_error(wi, warn_buf);
        }

        return ret;
    }
}

int get_eint(std::vector<t_inpfile> *inp, const std::string &name, int def,
             warninp_t wi)
{
    return get_eint(inp, name.c_str(), def, wi);
}

/* Note that sanitizing the trailing part of inp[ii].value was the responsibility of read_inpfile() */
int64_t get_eint64(std::vector<t_inpfile> *inp,
                   const char *name, int64_t def,
                   warninp_t wi)
{
    std::vector<t_inpfile> &inpRef = *inp;
    char                    buf[32], *ptr, warn_buf[STRLEN];

    int                     ii = get_einp(inp, name);

    if (ii == -1)
    {
        sprintf(buf, "%" PRId64, def);
        inpRef.back().value_.assign(buf);

        return def;
    }
    else
    {
        int64_t ret = str_to_int64_t(inpRef[ii].value_.c_str(), &ptr);
        if (*ptr != '\0')
        {
            sprintf(warn_buf, "Right hand side '%s' for parameter '%s' in parameter file is not an integer value\n", inpRef[ii].value_.c_str(), inpRef[ii].name_.c_str());
            warning_error(wi, warn_buf);
        }

        return ret;
    }
}

int64_t get_eint64(std::vector<t_inpfile> *inp,
                   const std::string &name, int64_t def,
                   warninp_t wi)
{
    return get_eint64(inp, name.c_str(), def, wi);
}

/* Note that sanitizing the trailing part of inp[ii].value was the responsibility of read_inpfile() */
double get_ereal(std::vector<t_inpfile> *inp, const char *name, double def,
                 warninp_t wi)
{
    std::vector<t_inpfile> &inpRef = *inp;
    char                    buf[32], *ptr, warn_buf[STRLEN];

    int                     ii = get_einp(inp, name);

    if (ii == -1)
    {
        sprintf(buf, "%g", def);
        inpRef.back().value_.assign(buf);

        return def;
    }
    else
    {
        double ret = strtod(inpRef[ii].value_.c_str(), &ptr);
        if (*ptr != '\0')
        {
            sprintf(warn_buf, "Right hand side '%s' for parameter '%s' in parameter file is not a real value\n", inpRef[ii].value_.c_str(), inpRef[ii].name_.c_str());
            warning_error(wi, warn_buf);
        }

        return ret;
    }
}

double get_ereal(std::vector<t_inpfile> *inp, const std::string &name, double def,
                 warninp_t wi)
{
    return get_ereal(inp, name.c_str(), def, wi);
}

/* Note that sanitizing the trailing part of inp[ii].value was the responsibility of read_inpfile() */
const char *get_estr(std::vector<t_inpfile> *inp, const char *name, const char *def)
{
    std::vector<t_inpfile> &inpRef = *inp;
    char                    buf[32];

    int                     ii = get_einp(inp, name);

    if (ii == -1)
    {
        if (def)
        {
            sprintf(buf, "%s", def);
            inpRef.back().value_.assign(buf);
        }
        else
        {
            inpRef.back().value_.clear();
        }

        return def;
    }
    else
    {
        return inpRef[ii].value_.c_str();
    }
}

const char *get_estr(std::vector<t_inpfile> *inp, const std::string &name, const char *def)
{
    return get_estr(inp, name.c_str(), def);
}

/* Note that sanitizing the trailing part of inp[ii].value was the responsibility of read_inpfile() */
int get_eeenum(std::vector<t_inpfile> *inp, const char *name, const char **defs,
               warninp_t wi)
{
    std::vector<t_inpfile> &inpRef = *inp;
    int                     n      = 0;
    char                    buf[STRLEN];

    int                     ii = get_einp(inp, name);

    if (ii == -1)
    {
        inpRef.back().value_.assign(defs[0]);

        return 0;
    }
    int i = 0;
    for (i = 0; (defs[i] != nullptr); i++)
    {
        if (gmx_strcasecmp_min(defs[i], inpRef[ii].value_.c_str()) == 0)
        {
            break;
        }
    }

    if (defs[i] == nullptr)
    {
        n += sprintf(buf, "Invalid enum '%s' for variable %s, using '%s'\n",
                     inpRef[ii].value_.c_str(), name, defs[0]);
        n += sprintf(buf+n, "Next time use one of:");
        int j  = 0;
        while (defs[j])
        {
            n += sprintf(buf+n, " '%s'", defs[j]);
            j++;
        }
        if (wi != nullptr)
        {
            warning_error(wi, buf);
        }
        else
        {
            fprintf(stderr, "%s\n", buf);
        }

        inpRef[ii].value_ = gmx_strdup(defs[0]);

        return 0;
    }

    return i;
}

int get_eeenum(std::vector<t_inpfile> *inp, const std::string &name, const char **defs,
               warninp_t wi)
{
    return get_eeenum(inp, name.c_str(), defs, wi);
}

int get_eenum(std::vector<t_inpfile> *inp, const char *name, const char **defs)
{
    return get_eeenum(inp, name, defs, nullptr);
}

void
printStringNewline(std::vector<t_inpfile> *inp, const char *line)
{
    std::string tmp("\n; ");
    tmp.append(line);
    get_estr(inp, tmp.c_str(), nullptr);
}

void
printStringNoNewline(std::vector<t_inpfile> *inp, const char *line)
{
    std::string tmp("; ");
    tmp.append(line);
    get_estr(inp, tmp.c_str(), nullptr);
}
void
setStringEntry(std::vector<t_inpfile> *inp, const char *name, char *newName, const char *def)
{
    const char *found = nullptr;
    found = get_estr(inp, name, def);
    if (found != nullptr)
    {
        std::strcpy(newName, found);
    }
}
