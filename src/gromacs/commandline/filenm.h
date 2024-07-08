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
/*! \file
 * \brief
 * Declares t_filenm for old-style command-line parsing of file name options.
 *
 * \inpublicapi
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_FILENM_H
#define GMX_COMMANDLINE_FILENM_H

#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "gromacs/fileio/filetypes.h"


//! \addtogroup module_commandline
//! \{

namespace gmx
{
template<typename>
class ArrayRef;
} // namespace gmx

/*! \brief
 * File name option definition for C code.
 *
 * \inpublicapi
 */
struct t_filenm
{
    int                      ftp; //!< File type, see enum in filetypes.h
    const char*              opt; //!< Command line option, can be nullptr in which case the commandline module, including all opt2??? functions below, will use the default option for the file type
    const char*              fn; //!< File name (as set in source code), can be nullptr in which case the commandline module will use the default file name for the file type
    unsigned long            flag;      //!< Flag for all kinds of info (see defs)
    std::vector<std::string> filenames; //!< File names
};

//! Whether a file name option is set.
#define ffSET 1 << 0
//! Whether a file name option specifies an input file.
#define ffREAD 1 << 1
//! Whether a file name option specifies an output file.
#define ffWRITE 1 << 2
//! Whether a file name option specifies an optional file.
#define ffOPT 1 << 3
//! Whether a file name option specifies a library file.
#define ffLIB 1 << 4
//! Whether a file name option accepts multiple file names.
#define ffMULT 1 << 5
//! Whether an input file name option accepts non-existent files.
#define ffALLOW_MISSING 1 << 6
//! Convenience flag for an input/output file.
#define ffRW (ffREAD | ffWRITE)
//! Convenience flag for an optional input file.
#define ffOPTRD (ffREAD | ffOPT)
//! Convenience flag for an optional output file.
#define ffOPTWR (ffWRITE | ffOPT)
//! Convenience flag for an optional input/output file.
#define ffOPTRW (ffRW | ffOPT)
//! Convenience flag for a library input file.
#define ffLIBRD (ffREAD | ffLIB)
//! Convenience flag for an optional library input file.
#define ffLIBOPTRD (ffOPTRD | ffLIB)
//! Convenience flag for an input file that accepts multiple files.
#define ffRDMULT (ffREAD | ffMULT)
//! Convenience flag for an optional input file that accepts multiple files.
#define ffOPTRDMULT (ffRDMULT | ffOPT)
//! Convenience flag for an output file that accepts multiple files.
#define ffWRMULT (ffWRITE | ffMULT)
//! Convenience flag for an optional output file that accepts multiple files.
#define ffOPTWRMULT (ffWRMULT | ffOPT)

/*! \brief
 * Returns the filename belonging to cmd-line option opt, or NULL when
 * no such option.
 */
const char* opt2fn(const char* opt, int nfile, const t_filenm fnm[]);

/*! \brief
 * Returns the filenames belonging to cmd-line option opt.
 *
 * An assertion will fail when the option does not exist.
 */
gmx::ArrayRef<const std::string> opt2fns(const char* opt, int nfile, const t_filenm fnm[]);

/*! \brief
 * Returns the filenames belonging to cmd-line option opt when set,
 * returns an empty vector when the option is not set.
 *
 * An assertion will fail when the option does not exist.
 */
gmx::ArrayRef<const std::string> opt2fnsIfOptionSet(const char* opt, int nfile, const t_filenm fnm[]);

//! Returns a file pointer from the filename.
#define opt2FILE(opt, nfile, fnm, mode) gmx_ffopen(opt2fn(opt, nfile, fnm), mode)

//! Returns the first file name with type ftp, or NULL when none found.
const char* ftp2fn(int ftp, int nfile, const t_filenm fnm[]);

/*! \brief
 * Returns the filenames for the first option with type ftp.
 *
 * An assertion will fail when when none found.
 */
gmx::ArrayRef<const std::string> ftp2fns(int ftp, int nfile, const t_filenm fnm[]);

//! Returns a file pointer from the file type.
#define ftp2FILE(ftp, nfile, fnm, mode) gmx_ffopen(ftp2fn(ftp, nfile, fnm), mode)

//! Returns TRUE when this file type has been found on the cmd-line.
bool ftp2bSet(int ftp, int nfile, const t_filenm fnm[]);

//! Returns TRUE when this option has been found on the cmd-line.
bool opt2bSet(const char* opt, int nfile, const t_filenm fnm[]);

/*! \brief
 * DEPRECATED Returns the file name belonging top cmd-line option opt,
 * or NULL when no such option.
 *
 * Also return NULL when opt is optional and option is not set.
 */
const char* opt2fn_null(const char* opt, int nfile, const t_filenm fnm[]);

/*! \brief
 * Returns the file name belonging top cmd-line option opt, or
 * std::nullopt when no such option.
 *
 * Also return std::nullopt when opt is optional and option is not set.
 */
std::optional<std::filesystem::path> opt2path_optional(const char* opt, int nfile, const t_filenm fnm[]);

/*! \brief
 * DEPRECATED Returns the first file name with type ftp, or NULL when none found.
 *
 * Also return NULL when ftp is optional and option is not set.
 */
const char* ftp2fn_null(int ftp, int nfile, const t_filenm fnm[]);

/*! \brief
 * Returns the first file name with type ftp, or std::nullopt when none found.
 *
 * Also return std::nullopt when ftp is optional and option is not set.
 */
std::optional<std::filesystem::path> ftp2path_optional(int ftp, int nfile, const t_filenm fnm[]);

//! Returns whether or not this filenm is optional.
bool is_optional(const t_filenm* fnm);

//! Returns whether or not this filenm is output.
bool is_output(const t_filenm* fnm);

//! Returns whether or not this filenm is set.
bool is_set(const t_filenm* fnm);

/*! \brief Return whether \c filename might have been produced by mdrun -noappend.
 *
 * If so, it must match "prefix.partNNNN.extension", for four decimal
 * digits N and non-empty prefix and extension. */
bool hasSuffixFromNoAppend(std::string_view filename);

/*! \brief
 * When we do checkpointing, this routine is called to check for previous
 * output files and append a '.partNNNN' suffix before the (output) file extensions.
 * If there was already a '.partNNNN' suffix before the file extension, that
 * is removed before the new suffix is added.
 */
int add_suffix_to_output_names(t_filenm* fnm, int nfile, const char* suffix);

//! \}

#endif
