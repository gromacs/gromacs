/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017, by the GROMACS development team, led by
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
 * Declares \c t_pargs, parse_common_args() and related methods.
 *
 * \inpublicapi
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_PARGS_H
#define GMX_COMMANDLINE_PARGS_H

#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_output_env_t;

#ifdef __cplusplus
extern "C"
{
#endif

/*! \addtogroup module_commandline
 * \{
 */

/** Command line argument type. */
enum
{
    etINT, etINT64, etREAL, etTIME, etSTR, etBOOL, etRVEC, etENUM, etNR
};

/*! \brief
 * Command-line argument definition for C code.
 *
 * \inpublicapi
 */
typedef struct
{
    /** Name of the argument (with leading dash included). */
    const char *option;
    /** Whether the argument is set (should be initialized to `FALSE`). */
    gmx_bool    bSet;
    /** Type of the argument (one of the enums in pargs.h). */
    int         type;
    /*! \brief
     * Pointer to variable that is to receive the value.
     *
     * The expected type depends on the value of \a type.
     * If an argument is not set by the user, then the pointed value is not
     * changed.  In other words, the initial value for the variable defines the
     * default value.
     */
    union
    {
        /*! \brief
         * Generic pointer for operations that do not need type information.
         *
         * Needs to be the first member to use initialized arrays.
         */
        void            *v;
        /** Integer value for etINT. */
        int             *i;
        /** Integer value for etINT64. */
        gmx_int64_t     *is;
        /** Real value for etREAL and etTIME. */
        real            *r;
        /*! \brief
         * String value for etSTR and etENUM.
         *
         * For etSTR, this should point to a `const char *` variable, which
         * will receive a pointer to the string value (caller should not free
         * the string).
         *
         * For etENUM, this should be an array of `const char *` values, where
         * the first and last values are `NULL`, and the other values define
         * the allowed enum values.  The first non-NULL value is the default
         * value.  After the arguments are parsed, the first element in the array
         * points to the selected enum value (pointers will be equal).
         */
        const char     **c;
        /** Boolean value for etBOOL. */
        gmx_bool        *b;
        /** Vector value for etRVEC. */
        rvec            *rv;
    }           u;
    /*! \brief
     * Description for the argument.
     *
     * If the string starts with `HIDDEN`, then the argument is hidden from
     * normal help listing and shell completions.
     */
    const char *desc;
} t_pargs;

/*! \brief
 * Returns ordinal number for an etENUM argument.
 *
 * \param[in] enumc  Array passed to `t_pargs` for an etENUM argument.
 * \returns   Index of selected enum value in the array.
 *
 * See t_pargs::u::c for the expected format of the array, including how the
 * first element should be initialized.
 * Note that the return value starts at one instead of zero: if the first enum
 * value is selected, this returns 1.
 */
int nenum(const char *const enumc[]);

/*! \brief
 * Returns value of an etINT option.
 *
 * \param[in] option Name of etINT argument to query.
 * \param[in] nparg  Number of elements in \p pa.
 * \param[in] pa     Array of arguments.
 * \returns   Value of \p option.
 *
 * \p option must specify a valid argument in \p pa of the correct type.
 */
int opt2parg_int(const char *option, int nparg, t_pargs pa[]);

/*! \brief
 * Returns value of an etBOOL option.
 *
 * \param[in] option Name of etBOOL argument to query.
 * \param[in] nparg  Number of elements in \p pa.
 * \param[in] pa     Array of arguments.
 * \returns   Value of \p option.
 *
 * \p option must specify a valid argument in \p pa of the correct type.
 */
gmx_bool opt2parg_bool(const char *option, int nparg, t_pargs pa[]);

/*! \brief
 * Returns value of an etREAL/etTIME option.
 *
 * \param[in] option Name of etREAL/etTIME argument to query.
 * \param[in] nparg  Number of elements in \p pa.
 * \param[in] pa     Array of arguments.
 * \returns   Value of \p option.
 *
 * \p option must specify a valid argument in \p pa of the correct type.
 */
real opt2parg_real(const char *option, int nparg, t_pargs pa[]);

/*! \brief
 * Returns value of an etSTR option.
 *
 * \param[in] option Name of etSTR argument to query.
 * \param[in] nparg  Number of elements in \p pa.
 * \param[in] pa     Array of arguments.
 * \returns   Value of \p option.
 *
 * \p option must specify a valid argument in \p pa of the correct type.
 */
const char *opt2parg_str(const char *option, int nparg, t_pargs pa[]);

/*! \brief
 * Returns value of an etENUM option.
 *
 * \param[in] option Name of etENUM argument to query.
 * \param[in] nparg  Number of elements in \p pa.
 * \param[in] pa     Array of arguments.
 * \returns   Value of \p option.
 *
 * \p option must specify a valid argument in \p pa of the correct type.
 */
const char *opt2parg_enum(const char *option, int nparg, t_pargs pa[]);

/*! \brief
 * Returns whether an argument has been set.
 *
 * \param[in] option Name of argument to check.
 * \param[in] nparg  Number of elements in \p pa.
 * \param[in] pa     Array of arguments.
 * \returns   `true` if \p option has been set.
 *
 * \p option must specify a valid argument in \p pa.
 */
gmx_bool opt2parg_bSet(const char *option, int nparg, t_pargs pa[]);


/** Add option -w to view output files (must be implemented in program). */
#define PCA_CAN_VIEW       (1<<5)
/** Add option to set begin time for trajectory reading. */
#define PCA_CAN_BEGIN      (1<<6)
/** Add option to set end time for trajectory reading. */
#define PCA_CAN_END        (1<<7)
/** Add option to set time step for trajectory reading. */
#define PCA_CAN_DT         (1<<14)
/** Add all options for trajectory time control. */
#define PCA_CAN_TIME       (PCA_CAN_BEGIN | PCA_CAN_END | PCA_CAN_DT)
/** Add option -tu to set time unit for output. */
#define PCA_TIME_UNIT      (1<<15)
/** Add option -deffnm to set default for all file options. */
#define PCA_CAN_SET_DEFFNM (1<<10)
/** Do not raise a fatal error when invalid options are encountered. */
#define PCA_NOEXIT_ON_ARGS (1<<11)
/** Is this node not reading: for parallel all nodes but the master */
#define PCA_NOT_READ_NODE  (1<<16)
/** Don't do any special processing for ffREAD files */
#define PCA_DISABLE_INPUT_FILE_CHECKING (1<<17)

/*! \brief
 * Parse command-line arguments.
 *
 * Some common default arguments are also recognized in addition to those
 * provided through \p pa.  The set of recognized default arguments is affected
 * by \p Flags.
 *
 * Recognized arguments are removed from the list.
 *
 * For full functionality, this function needs to be used within a function
 * that is passed to gmx_run_cmain().  It should be called as the first thing in
 * that function.  Initialization code can be executed before it, but you need
 * to be aware that if the program is executed with -h and MPI, the code before
 * parse_common_args() only executes on the master node.
 *
 * If the return value is `FALSE`, the program should return immediately (this
 * is necessary for -h and a few other cases).
 *
 * \see gmx_run_cmain().
 */
gmx_bool parse_common_args(int *argc, char *argv[], unsigned long Flags,
                           int nfile, t_filenm fnm[], int npargs, t_pargs *pa,
                           int ndesc, const char **desc,
                           int nbugs, const char **bugs,
                           gmx_output_env_t **oenv);

/*! \} */

#ifdef __cplusplus
}
#endif

#endif
