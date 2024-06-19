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
#ifndef GMX_FILEIO_READINP_H
#define GMX_FILEIO_READINP_H

#include <cstdint>
#include <cstdio>
#include <cstring>

#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include "gromacs/fileio/warninp.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/stringutil.h"

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
std::vector<t_inpfile> read_inpfile(gmx::TextInputStream*        stream,
                                    const std::filesystem::path& fn,
                                    WarningHandler*              wi);

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
void write_inpfile(gmx::TextOutputStream*       stream,
                   const std::filesystem::path& fn,
                   std::vector<t_inpfile>*      inp,
                   gmx_bool                     bHaltOnUnknown,
                   WriteMdpHeader               writeHeader,
                   WarningHandler*              wi);
/* Write inp to fn, warning (and perhaps halting) if any fields are
 * unknown. The helpful header contains irreproducible content, so
 * its writing can be suppressed to make testing more useful. */

void replace_inp_entry(gmx::ArrayRef<t_inpfile> inp, const char* old_entry, const char* new_entry);

int search_einp(gmx::ArrayRef<const t_inpfile> inp, const char* name);
/* Return the index of an .mdp field with the given name within the
 * inp vector, if it exists. Return -1 if it does not exist. */

void mark_einp_set(gmx::ArrayRef<t_inpfile> inp, const char* name);

int get_eint(std::vector<t_inpfile>* inp, const char* name, int def, WarningHandler* wi);
int get_eint(std::vector<t_inpfile>* inp, const std::string& name, int def, WarningHandler* wi);

int64_t get_eint64(std::vector<t_inpfile>* inp, const char* name, int64_t def, WarningHandler* wi);
int64_t get_eint64(std::vector<t_inpfile>* inp, const std::string& name, int64_t def, WarningHandler* wi);

double get_ereal(std::vector<t_inpfile>* inp, const char* name, double def, WarningHandler* wi);
double get_ereal(std::vector<t_inpfile>* inp, const std::string& name, double def, WarningHandler* wi);

const char* get_estr(std::vector<t_inpfile>* inp, const char* name, const char* def);
const char* get_estr(std::vector<t_inpfile>* inp, const std::string& name, const char* def);

int get_eeenum(std::vector<t_inpfile>* inp, const char* name, const char* const* defs, WarningHandler* wi);
/* defs must be NULL terminated */
int get_eeenum(std::vector<t_inpfile>* inp, const std::string& name, const char* const* defs, WarningHandler* wi);
/* defs must be NULL terminated */

int get_eenum(std::vector<t_inpfile>* inp, const char* name, const char* const* defs);
/* defs must be NULL terminated */

//! Get index of option `name`. Exposed here so that `getEnum` can access it.
int get_einp(std::vector<t_inpfile>* inp, const char* name);

/*! \brief Read option from input and return corresponding enum value
 *
 * If the option is not set, return the first value of the enum as default.
 * Defined here so we don't need to instantiate the templates in the source file.
 *
 * \tparam EnumType  The type of enum to be returned
 * \param[in]  inp   The input file vector
 * \param[in]  name  The name of the option to be read
 * \param[out] wi    Handler for context-sensitive warnings.
 * \return  Enum value corresponding to read input
 */
template<typename EnumType>
EnumType getEnum(std::vector<t_inpfile>* inp, const char* name, WarningHandler* wi)
{
    // If there's no valid option, we'll use the EnumType::Default.
    // Note, this assumes the enum is zero based, which is also assumed by
    // EnumerationWrapper and EnumerationArray.
    const auto  defaultEnumValue = EnumType::Default;
    const auto& defaultName      = enumValueToString(defaultEnumValue);
    // Get index of option in input
    const auto ii = get_einp(inp, name);
    if (ii == -1)
    {
        // If the option wasn't set, we return EnumType::Default
        inp->back().value_.assign(defaultName);
        return defaultEnumValue;
    }

    // Check if option string can be mapped to a valid enum value.
    //
    // Note that this cannot be replaced with
    // StringToEnumValueConverter until all instantiations of this
    // function have a matching enumValueToString, and all of the
    // latter are in the same namespace. Currently some of those
    // function declarations are in gmx namespace and some are not.
    const auto* optionString = (*inp)[ii].value_.c_str();
    for (auto enumValue : gmx::EnumerationWrapper<EnumType>{})
    {
        if (gmx_strcasecmp_min(enumValueToString(enumValue), optionString) == 0)
        {
            return enumValue;
        }
    }

    // If we get here, the option set is invalid. Print error.
    std::string errorMessage = gmx::formatString(
            "Invalid enum '%s' for variable %s, using '%s'\n", optionString, name, defaultName);
    errorMessage += gmx::formatString("Next time, use one of:");
    for (auto enumValue : gmx::EnumerationWrapper<EnumType>{})
    {
        errorMessage += gmx::formatString(" '%s'", enumValueToString(enumValue));
    }
    if (wi != nullptr)
    {
        wi->addError(errorMessage);
    }
    else
    {
        fprintf(stderr, "%s\n", errorMessage.c_str());
    }
    (*inp)[ii].value_.assign(defaultName);
    return defaultEnumValue;
}


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
