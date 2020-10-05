/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
#ifndef GMX_FILEIO_READINP_H
#define GMX_FILEIO_READINP_H

#include <cstring>

#include <string>
#include <utility>
#include <vector>

#include "gromacs/utility/basedefinitions.h"

struct warninp;
typedef warninp* warninp_t;

namespace gmx
{
template<typename>
class ArrayRef;
class KeyValueTreeObject;
class TextInputStream;
class TextOutputStream;
} // namespace gmx

/* !\brief Input file structure that is populated with entries read from a file.
 *
 * This structure contains the information read from the mdp file that is later used
 * to build the ir from it. It is first constructed with a set of names and values,
 * and later populated when checking against the available options in readir.
 * Uses the functions below to both check entries and set the information.
 */
struct t_inpfile
{
    /*!\brief Minimum allowed constructor sets all elements */
    t_inpfile(int         count,
              int         inp_count,
              bool        bObsolete,
              bool        bSet,
              bool        bHandledAsKeyValueTree,
              std::string name,
              std::string value) :
        count_(count),
        bObsolete_(bObsolete),
        bSet_(bSet),
        bHandledAsKeyValueTree_(bHandledAsKeyValueTree),
        name_(std::move(name)),
        value_(std::move(value)),
        inp_count_(inp_count)
    {
    }
    int  count_;                  /* sort order for output  */
    bool bObsolete_;              /* whether it is an obsolete param value */
    bool bSet_;                   /* whether it it has been read out */
    bool bHandledAsKeyValueTree_; /* whether it it has been handled with key-value machinery */
    std::string name_;            /* name of the parameter */
    std::string value_;           /* parameter value string */
    int         inp_count_;       /* number of einps read. Only valid for the first item
                                                          in the inpfile list. */
};

/*! \brief Create and return a vector of t_inpfile structs
 * from "key = value" lines in \c stream corresponding to file \c fn.
 *
 * \param[in]  stream          Text stream to read.
 * \param[in]  fn              Filename corresponding to \c reader.
 * \param[out] wi              Handler for context-sensitive warnings.
 * \throws     std::bad_alloc  If out of memory.
 * \throws     Anything the stream underlying \c reader can throw. */
std::vector<t_inpfile> read_inpfile(gmx::TextInputStream* stream, const char* fn, warninp_t wi);

gmx::KeyValueTreeObject flatKeyValueTreeFromInpFile(gmx::ArrayRef<const t_inpfile> inp);

enum class WriteMdpHeader
{
    no,
    yes
};

/*! \brief Write "key = value" lines from \c inp to \c stream.
 *
 * \param[in]  stream          Text stream to write.
 * \param[in]  fn              Filename corresponding to \c stream.
 * \param[in]  inp             vector of key-value pairs.
 * \param[in]  bHaltOnUnknown  Whether to issue a fatal error if an unknown key is found.
 * \param[in]  writeHeader     Whether to write a header recording some context a user might like.
 * \param[out] wi              Handler for context-sensitive warnings.
 * \throws     std::bad_alloc  If out of memory.
 * \throws     Anything the stream underlying \c writer can throw. */
void write_inpfile(gmx::TextOutputStream*  stream,
                   const char*             fn,
                   std::vector<t_inpfile>* inp,
                   gmx_bool                bHaltOnUnknown,
                   WriteMdpHeader          writeHeader,
                   warninp_t               wi);
/* Write inp to fn, warning (and perhaps halting) if any fields are
 * unknown. The helpful header contains irreproducible content, so
 * its writing can be suppressed to make testing more useful. */

void replace_inp_entry(gmx::ArrayRef<t_inpfile> inp, const char* old_entry, const char* new_entry);

int search_einp(gmx::ArrayRef<const t_inpfile> inp, const char* name);
/* Return the index of an .mdp field with the given name within the
 * inp vector, if it exists. Return -1 if it does not exist. */

void mark_einp_set(gmx::ArrayRef<t_inpfile> inp, const char* name);

int get_eint(std::vector<t_inpfile>* inp, const char* name, int def, warninp_t wi);
int get_eint(std::vector<t_inpfile>* inp, const std::string& name, int def, warninp_t wi);

int64_t get_eint64(std::vector<t_inpfile>* inp, const char* name, int64_t def, warninp_t wi);
int64_t get_eint64(std::vector<t_inpfile>* inp, const std::string& name, int64_t def, warninp_t wi);

double get_ereal(std::vector<t_inpfile>* inp, const char* name, double def, warninp_t wi);
double get_ereal(std::vector<t_inpfile>* inp, const std::string& name, double def, warninp_t wi);

const char* get_estr(std::vector<t_inpfile>* inp, const char* name, const char* def);
const char* get_estr(std::vector<t_inpfile>* inp, const std::string& name, const char* def);

int get_eeenum(std::vector<t_inpfile>* inp, const char* name, const char** defs, warninp_t wi);
/* defs must be NULL terminated */
int get_eeenum(std::vector<t_inpfile>* inp, const std::string& name, const char** defs, warninp_t wi);
/* defs must be NULL terminated */

int get_eenum(std::vector<t_inpfile>* inp, const char* name, const char** defs);
/* defs must be NULL terminated */

//! Replace for macro CCTYPE, prints comment string after newline
void printStringNewline(std::vector<t_inpfile>* inp, const char* line);
//! Replace for macro CTYPE, prints comment string
void printStringNoNewline(std::vector<t_inpfile>* inp, const char* line);
//! Replace for macro STYPE, checks for existing string entry and if possible replaces it
void setStringEntry(std::vector<t_inpfile>* inp, const char* name, char* newName, const char* def);

/*! \brief
 * Returns a string value and sets the value in \p inp
 *
 * The value is either from \p inp when \p name is found or \p def otherwise.
 *
 * \note this is a wrapper function for g_estr()
 */
std::string setStringEntry(std::vector<t_inpfile>* inp, const std::string& name, const std::string& def);

#endif
