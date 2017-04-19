/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Declares t_filenm for old-style command-line parsing of file name options.
 *
 * \inpublicapi
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_FILENM_H
#define GMX_COMMANDLINE_FILENM_H

#include "gromacs/fileio/filetypes.h"
#include "gromacs/utility/basedefinitions.h"


//! \addtogroup module_commandline
//! \{

/*! \brief
 * File name option definition for C code.
 *
 * \inpublicapi
 */
struct t_filenm {
    int           ftp;    //!< File type (see enum in filetypes.h)
    const char   *opt;    //!< Command line option
    const char   *fn;     //!< File name (as set in source code)
    unsigned long flag;   //!< Flag for all kinds of info (see defs)
    int           nfiles; //!< number of files
    char        **fns;    //!< File names
};

//! Whether a file name option is set.
#define ffSET   1<<0
//! Whether a file name option specifies an input file.
#define ffREAD  1<<1
//! Whether a file name option specifies an output file.
#define ffWRITE 1<<2
//! Whether a file name option specifies an optional file.
#define ffOPT   1<<3
//! Whether a file name option specifies a library file.
#define ffLIB   1<<4
//! Whether a file name option accepts multiple file names.
#define ffMULT  1<<5
//! Whether an input file name option accepts non-existent files.
#define ffALLOW_MISSING 1<<6
//! Convenience flag for an input/output file.
#define ffRW    (ffREAD | ffWRITE)
//! Convenience flag for an optional input file.
#define ffOPTRD (ffREAD | ffOPT)
//! Convenience flag for an optional output file.
#define ffOPTWR (ffWRITE| ffOPT)
//! Convenience flag for an optional input/output file.
#define ffOPTRW (ffRW   | ffOPT)
//! Convenience flag for a library input file.
#define ffLIBRD (ffREAD | ffLIB)
//! Convenience flag for an optional library input file.
#define ffLIBOPTRD (ffOPTRD | ffLIB)
//! Convenience flag for an input file that accepts multiple files.
#define ffRDMULT   (ffREAD  | ffMULT)
//! Convenience flag for an optional input file that accepts multiple files.
#define ffOPTRDMULT   (ffRDMULT | ffOPT)
//! Convenience flag for an output file that accepts multiple files.
#define ffWRMULT   (ffWRITE  | ffMULT)
//! Convenience flag for an optional output file that accepts multiple files.
#define ffOPTWRMULT   (ffWRMULT | ffOPT)

/*! \brief
 * Returns the filename belonging to cmd-line option opt, or NULL when
 * no such option.
 */
const char *opt2fn(const char *opt, int nfile, const t_filenm fnm[]);

/*! \brief
 * Returns the filenames belonging to cmd-line option opt, or NULL when
 * no such option.
 */
int opt2fns(char **fns[], const char *opt, int nfile,
            const t_filenm fnm[]);

/*! \brief
 * Return a pointer to the t_filenm data structure of filenames belonging to
 * command-line option opt, or NULL when no such option was used.
 */
const t_filenm *getFilenm(const char *opt, int nfile, const t_filenm fnm[]);

//! Returns a file pointer from the filename.
#define opt2FILE(opt, nfile, fnm, mode) gmx_ffopen(opt2fn(opt, nfile, fnm), mode)

//! Returns the first file name with type ftp, or NULL when none found.
const char *ftp2fn(int ftp, int nfile, const t_filenm fnm[]);

/*! \brief
 * Returns the number of files for the first option with type ftp
 * and the files in **fns[] (will be allocated), or NULL when none found.
 */
int ftp2fns(char **fns[], int ftp, int nfile, const t_filenm fnm[]);

//! Returns a file pointer from the file type.
#define ftp2FILE(ftp, nfile, fnm, mode) gmx_ffopen(ftp2fn(ftp, nfile, fnm), mode)

//! Returns TRUE when this file type has been found on the cmd-line.
gmx_bool ftp2bSet(int ftp, int nfile, const t_filenm fnm[]);

//! Returns TRUE when this option has been found on the cmd-line.
gmx_bool opt2bSet(const char *opt, int nfile, const t_filenm fnm[]);

/*! \brief
 * Returns the file name belonging top cmd-line option opt, or NULL when
 * no such option.
 *
 * Also return NULL when opt is optional and option is not set.
 */
const char *opt2fn_null(const char *opt, int nfile, const t_filenm fnm[]);

/*! \brief
 * Returns the first file name with type ftp, or NULL when none found.
 *
 * Also return NULL when ftp is optional and option is not set.
 */
const char *ftp2fn_null(int ftp, int nfile, const t_filenm fnm[]);

//! Returns whether or not this filenm is optional.
gmx_bool is_optional(const t_filenm *fnm);

//! Returns whether or not this filenm is output.
gmx_bool is_output(const t_filenm *fnm);

//! Returns whether or not this filenm is set.
gmx_bool is_set(const t_filenm *fnm);

/*! \brief
 * When we do checkpointing, this routine is called to check for previous
 * output files and append a '.partNNNN' suffix before the (output) file extensions.
 */
int add_suffix_to_output_names(t_filenm *fnm, int nfile, const char *suffix);

/*! \brief
 * Duplicates the filename list (to make a private copy for each thread,
 * for example).
 */
t_filenm *dup_tfn(int nf, const t_filenm tfn[]);

//! Frees memory allocated for file names by parse_common_args().
void done_filenms(int nf, t_filenm fnm[]);

//! \}

#endif
