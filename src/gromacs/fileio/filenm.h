/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#ifndef GMX_FILEIO_FILENM_H
#define GMX_FILEIO_FILENM_H

#include "gromacs/legacyheaders/types/commrec_fwd.h"
#include "gromacs/utility/basedefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

/* this enum should correspond to the array deffile in gmxlib/filenm.c */
enum {
    efMDP,
    efTRX, efTRO, efTRN, efTRR, efCOMPRESSED, efXTC, efTNG,
    efEDR,
    efSTX, efSTO, efGRO, efG96, efPDB, efBRK, efENT, efESP, efPQR,
    efCPT,
    efLOG, efXVG, efOUT,
    efNDX,
    efTOP, efITP,
    efTPS, efTPR,
    efTEX, efRTP, efATP, efHDB,
    efDAT, efDLG,
    efMAP, efEPS, efMAT, efM2P,
    efMTX,
    efEDI,
    efCUB,
    efXPM,
    efRND,
    efNR
};

typedef struct {
    int           ftp;    /* File type (see enum above)		*/
    const char   *opt;    /* Command line option			*/
    const char   *fn;     /* File name (as set in source code)	*/
    unsigned long flag;   /* Flag for all kinds of info (see defs)*/
    int           nfiles; /* number of files			*/
    char        **fns;    /* File names				*/
} t_filenm;

#define ffSET   1<<0
#define ffREAD  1<<1
#define ffWRITE 1<<2
#define ffOPT   1<<3
#define ffLIB   1<<4
#define ffMULT  1<<5
#define ffALLOW_MISSING 1<<6
#define ffRW    (ffREAD | ffWRITE)
#define ffOPTRD (ffREAD | ffOPT)
#define ffOPTWR (ffWRITE| ffOPT)
#define ffOPTRW (ffRW   | ffOPT)
#define ffLIBRD (ffREAD | ffLIB)
#define ffLIBOPTRD (ffOPTRD | ffLIB)
#define ffRDMULT   (ffREAD  | ffMULT)
#define ffOPTRDMULT   (ffRDMULT | ffOPT)
#define ffWRMULT   (ffWRITE  | ffMULT)
#define ffOPTWRMULT   (ffWRMULT | ffOPT)

const char *ftp2ext(int ftp);
/* Return extension for filetype */

const char *ftp2ext_generic(int ftp);
/* Return extension for filetype, and a generic name for generic types
   (e.g. trx)*/

const char *ftp2ext_with_dot(int ftp);
/* Return extension for filetype with a leading dot */

int ftp2generic_count(int ftp);
/* Return the number of filetypes for a generic filetype */

const int *ftp2generic_list(int ftp);
/* Return the list of filetypes for a generic filetype */

const char *ftp2desc(int ftp);
/* Return description for file type */

const char *ftp2defnm(int ftp);
/* Return default file name for file type */

const char *ftp2defopt(int ftp);
/* Return default option name for file type */

const char *ftp2ftype(int ftp);
/* Return Binary or ASCII depending on file type */

const char *opt2fn(const char *opt, int nfile, const t_filenm fnm[]);
/* Return the filename belonging to cmd-line option opt, or NULL when
 * no such option. */

const char *opt2fn_master(const char *opt, int nfile,
                          const t_filenm fnm[], t_commrec *cr);
/* Return the filename belonging to cmd-line option opt, or NULL when
 * no such option or not running on master */


int opt2fns(char **fns[], const char *opt, int nfile,
            const t_filenm fnm[]);
/* Return the filenames belonging to cmd-line option opt, or NULL when
 * no such option. */

#define opt2FILE(opt, nfile, fnm, mode) gmx_ffopen(opt2fn(opt, nfile, fnm), mode)
/* Return a file pointer from the filename (see above) */

int fn2ftp(const char *fn);
/* Return the filetype corrsponding to filename */

const char *ftp2fn(int ftp, int nfile, const t_filenm fnm[]);
/* Return the first file name with type ftp, or NULL when none found. */

int ftp2fns(char **fns[], int ftp, int nfile, const t_filenm fnm[]);
/* Return the number of files for the first option with type ftp
   and the files in **fns[] (will be allocated), or NULL when none found. */

#define ftp2FILE(ftp, nfile, fnm, mode) gmx_ffopen(ftp2fn(ftp, nfile, fnm), mode)
/* Return a file pointer from the filename (see above) */

gmx_bool ftp2bSet(int ftp, int nfile, const t_filenm fnm[]);
/* Return TRUE when this file type has been found on the cmd-line */

gmx_bool opt2bSet(const char *opt, int nfile, const t_filenm fnm[]);
/* Return TRUE when this option has been found on the cmd-line */

const char *opt2fn_null(const char *opt, int nfile, const t_filenm fnm[]);
/* Return the filenm belonging top cmd-line option opt, or NULL when
 * no such option.
 * Also return NULL when opt is optional and option is not set.
 */

const char *ftp2fn_null(int ftp, int nfile, const t_filenm fnm[]);
/* Return the first file name with type ftp, or NULL when none found.
 * Also return NULL when ftp is optional and option is not set.
 */

gmx_bool is_optional(const t_filenm *fnm);
/* Return whether or not this filenm is optional */

gmx_bool is_output(const t_filenm *fnm);
/* Return whether or not this filenm is output */

gmx_bool is_set(const t_filenm *fnm);
/* Return whether or not this filenm is set */

/* When we do checkpointing, this routine is called to check for previous
 * output files and append a '.partNNNN' suffix before the (output) file extensions.
 */
int add_suffix_to_output_names(t_filenm *fnm, int nfile, const char *suffix);

/* duplicate the filename list (to make a private copy for each thread,
   for example) */
t_filenm *dup_tfn(int nf, const t_filenm tfn[]);

/* Free memory allocated for file names by parse_file_args(). */
void done_filenms(int nf, t_filenm fnm[]);

#ifdef __cplusplus
}
#endif

#endif  /* GMX_FILEIO_FILENM_H */
