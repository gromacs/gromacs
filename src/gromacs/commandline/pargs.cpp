/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
/* This file is completely threadsafe - keep it that way! */
#include "gromacs/commandline/pargs.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include "thread_mpi/threads.h"

#include "gromacs/legacyheaders/gmx_fatal.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/network.h"

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/commandline/cmdlinehelpwriter.h"
#include "gromacs/commandline/shellcompletions.h"
#include "gromacs/fileio/timecontrol.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

/* The source code in this file should be thread-safe.
      Please keep it that way. */

static void usage(const char *type, const char *arg)
{
    GMX_ASSERT(arg != NULL, "NULL command-line argument should not occur");
    gmx_fatal(FARGS, "Expected %s argument for option %s\n", type, arg);
}

/* Scan an int for argument argv[*i] from argument at argv[*i + 1].
 * eg: -p 32.  argv[*i] is only used for error reporting.
 * If there is no value, or the conversion is not successful, the
 * routine exits with an error, otherwise it returns the value found.
 * *i is incremented once.
 */
static int iscan(int argc, char *argv[], int *i)
{
    const char *const arg = argv[*i];
    if (argc <= (*i)+1)
    {
        usage("an integer", arg);
    }
    const char *const value = argv[++(*i)];
    char             *endptr;
    int               var = std::strtol(value, &endptr, 10);
    if (*value == '\0' || *endptr != '\0')
    {
        usage("an integer", arg);
    }
    return var;
}

/* Same as above, but for large integer values */
static gmx_int64_t istepscan(int argc, char *argv[], int *i)
{
    const char *const arg = argv[*i];
    if (argc <= (*i)+1)
    {
        usage("an integer", arg);
    }
    const char *const value = argv[++(*i)];
    char             *endptr;
    gmx_int64_t       var = str_to_int64_t(value, &endptr);
    if (*value == '\0' || *endptr != '\0')
    {
        usage("an integer", arg);
    }
    return var;
}

/* Routine similar to the above, but working on doubles. */
static double dscan(int argc, char *argv[], int *i)
{
    const char *const arg = argv[*i];
    if (argc <= (*i)+1)
    {
        usage("a real", arg);
    }
    const char *const value = argv[++(*i)];
    char             *endptr;
    double            var = std::strtod(value, &endptr);
    if (*value == '\0' || *endptr != '\0')
    {
        usage("a real", arg);
    }
    return var;
}

/* Routine similar to the above, but working on strings. The pointer
 * returned is a pointer to the argv field.
 */
static char *sscan(int argc, char *argv[], int *i)
{
    if (argc > (*i)+1)
    {
        if ( (argv[(*i)+1][0] == '-') && (argc > (*i)+2) &&
             (argv[(*i)+2][0] != '-') )
        {
            fprintf(stderr, "Possible missing string argument for option %s\n\n",
                    argv[*i]);
        }
    }
    else
    {
        usage("a string", argv[*i]);
    }

    return argv[++(*i)];
}

static gmx_bool is_hidden(t_pargs *pa)
{
    return (strstr(pa->desc, "HIDDEN") != NULL);
}

int nenum(const char *const enumc[])
{
    int i;

    i = 1;
    /* we *can* compare pointers directly here! */
    while (enumc[i] && enumc[0] != enumc[i])
    {
        i++;
    }

    return i;
}

/* Read a number of arguments from the command line.
 * For etINT, etREAL and etCHAR an extra argument is read (when present)
 * for etBOOL the gmx_boolean option is changed to the negate value
 */
static void get_pargs(int *argc, char *argv[], int nparg, t_pargs pa[])
{
    int       i, j, k, match;
    gmx_bool *bKeep;
    char      buf[32];
    char     *ptr;

    snew(bKeep, *argc+1);
    bKeep[0]     = TRUE;
    bKeep[*argc] = TRUE;

    for (i = 1; (i < *argc); i++)
    {
        bKeep[i] = TRUE;
        for (j = 0; (j < nparg); j++)
        {
            if (pa[j].type == etBOOL)
            {
                sprintf(buf, "-no%s", pa[j].option+1);
                if (strcmp(pa[j].option, argv[i]) == 0)
                {
                    *pa[j].u.b = TRUE;
                    pa[j].bSet = TRUE;
                    bKeep[i]   = FALSE;
                }
                else if (strcmp(buf, argv[i]) == 0)
                {
                    *pa[j].u.b = FALSE;
                    pa[j].bSet = TRUE;
                    bKeep[i]   = FALSE;
                }
            }
            else if (strcmp(pa[j].option, argv[i]) == 0)
            {
                if (pa[j].bSet)
                {
                    fprintf(stderr, "Setting option %s more than once!\n",
                            pa[j].option);
                }
                pa[j].bSet = TRUE;
                bKeep[i]   = FALSE;
                switch (pa[j].type)
                {
                    case etINT:
                        *pa[j].u.i = iscan(*argc, argv, &i);
                        break;
                    case etINT64:
                        *pa[j].u.is = istepscan(*argc, argv, &i);
                        break;
                    case etTIME:
                    case etREAL:
                        *pa[j].u.r = dscan(*argc, argv, &i);
                        break;
                    case etSTR:
                        *(pa[j].u.c) = sscan(*argc, argv, &i);
                        break;
                    case etENUM:
                        match = -1;
                        ptr   = sscan(*argc, argv, &i);
                        for (k = 1; (pa[j].u.c[k] != NULL); k++)
                        {
                            /* only check ptr against beginning of
                               pa[j].u.c[k] */
                            if (gmx_strncasecmp(ptr, pa[j].u.c[k], strlen(ptr)) == 0)
                            {
                                if ( ( match == -1 ) ||
                                     ( strlen(pa[j].u.c[k]) <
                                       strlen(pa[j].u.c[match]) ) )
                                {
                                    match = k;
                                }
                            }
                        }
                        if (match != -1)
                        {
                            pa[j].u.c[0] = pa[j].u.c[match];
                        }
                        else
                        {
                            gmx_fatal(FARGS, "Invalid argument %s for option %s",
                                      ptr, pa[j].option);
                        }
                        break;
                    case etRVEC:
                        (*pa[j].u.rv)[0] = dscan(*argc, argv, &i);
                        if ( (i+1 == *argc) ||
                             ( (argv[i+1][0] == '-') &&
                               !isdigit(argv[i+1][1]) ) )
                        {
                            (*pa[j].u.rv)[1]     =
                                (*pa[j].u.rv)[2] =
                                    (*pa[j].u.rv)[0];
                        }
                        else
                        {
                            bKeep[i]         = FALSE;
                            (*pa[j].u.rv)[1] = dscan(*argc, argv, &i);
                            if ( (i+1 == *argc) ||
                                 ( (argv[i+1][0] == '-') &&
                                   !isdigit(argv[i+1][1]) ) )
                            {
                                gmx_fatal(FARGS,
                                          "%s: vector must have 1 or 3 real parameters",
                                          pa[j].option);
                            }
                            bKeep[i]         = FALSE;
                            (*pa[j].u.rv)[2] = dscan(*argc, argv, &i);
                        }
                        break;
                    default:
                        gmx_fatal(FARGS, "Invalid type %d in pargs", pa[j].type);
                }
                /* i may be incremented, so set it to not keep */
                bKeep[i] = FALSE;
            }
        }
    }

    /* Remove used entries */
    for (i = j = 0; (i <= *argc); i++)
    {
        if (bKeep[i])
        {
            argv[j++] = argv[i];
        }
    }
    (*argc) = j-1;
    sfree(bKeep);
}

int opt2parg_int(const char *option, int nparg, t_pargs pa[])
{
    int i;

    for (i = 0; (i < nparg); i++)
    {
        if (strcmp(pa[i].option, option) == 0)
        {
            return *pa[i].u.i;
        }
    }

    gmx_fatal(FARGS, "No integer option %s in pargs", option);

    return 0;
}

gmx_bool opt2parg_bool(const char *option, int nparg, t_pargs pa[])
{
    int i;

    for (i = 0; (i < nparg); i++)
    {
        if (strcmp(pa[i].option, option) == 0)
        {
            return *pa[i].u.b;
        }
    }

    gmx_fatal(FARGS, "No boolean option %s in pargs", option);

    return FALSE;
}

real opt2parg_real(const char *option, int nparg, t_pargs pa[])
{
    int i;

    for (i = 0; (i < nparg); i++)
    {
        if (strcmp(pa[i].option, option) == 0)
        {
            return *pa[i].u.r;
        }
    }

    gmx_fatal(FARGS, "No real option %s in pargs", option);

    return 0.0;
}

const char *opt2parg_str(const char *option, int nparg, t_pargs pa[])
{
    int i;

    for (i = 0; (i < nparg); i++)
    {
        if (strcmp(pa[i].option, option) == 0)
        {
            return *(pa[i].u.c);
        }
    }

    gmx_fatal(FARGS, "No string option %s in pargs", option);

    return NULL;
}

gmx_bool opt2parg_bSet(const char *option, int nparg, t_pargs pa[])
{
    int i;

    for (i = 0; (i < nparg); i++)
    {
        if (strcmp(pa[i].option, option) == 0)
        {
            return pa[i].bSet;
        }
    }

    gmx_fatal(FARGS, "No such option %s in pargs", option);

    return FALSE; /* Too make some compilers happy */
}

const char *opt2parg_enum(const char *option, int nparg, t_pargs pa[])
{
    int i;

    for (i = 0; (i < nparg); i++)
    {
        if (strcmp(pa[i].option, option) == 0)
        {
            return pa[i].u.c[0];
        }
    }

    gmx_fatal(FARGS, "No such option %s in pargs", option);

    return NULL;
}

/********************************************************************
 * parse_common_args()
 */

static void set_default_time_unit(const char *time_list[], gmx_bool bCanTime)
{
    int         i      = 0;
    const char *select = NULL;

    if (bCanTime)
    {
        select = getenv("GMXTIMEUNIT");
        if (select != NULL)
        {
            i = 1;
            while (time_list[i] && strcmp(time_list[i], select) != 0)
            {
                i++;
            }
        }
    }
    if (!bCanTime || select == NULL ||
        time_list[i] == NULL || strcmp(time_list[i], select) != 0)
    {
        /* Set it to the default: ps */
        i = 1;
        while (time_list[i] && strcmp(time_list[i], "ps") != 0)
        {
            i++;
        }

    }
    time_list[0] = time_list[i];
}

static void set_default_xvg_format(const char *xvg_list[])
{
    int         i;
    const char *select;

    select = getenv("GMX_VIEW_XVG");
    if (select == NULL)
    {
        /* The default is the first option */
        xvg_list[0] = xvg_list[1];
    }
    else
    {
        i = 1;
        while (xvg_list[i] && strcmp(xvg_list[i], select) != 0)
        {
            i++;
        }
        if (xvg_list[i] != NULL)
        {
            xvg_list[0] = xvg_list[i];
        }
        else
        {
            xvg_list[0] = xvg_list[exvgNONE];
        }
    }
}

static int add_parg(int npargs, t_pargs *pa, t_pargs *pa_add)
{
    memcpy(&(pa[npargs]), pa_add, sizeof(*pa_add));

    return npargs+1;
}

namespace gmx
{

namespace
{

/*! \brief
 * Converts a t_filenm option into an Options option.
 *
 * \param     options Options object to add the new option to.
 * \param[in] fnm     t_filenm option to convert.
 *
 * \ingroup module_commandline
 */
void filenmToOptions(Options *options, const t_filenm *fnm)
{
    const bool        bRead     = ((fnm->flag & ffREAD)  != 0);
    const bool        bWrite    = ((fnm->flag & ffWRITE) != 0);
    const bool        bOptional = ((fnm->flag & ffOPT)   != 0);
    const bool        bLibrary  = ((fnm->flag & ffLIB)   != 0);
    const bool        bMultiple = ((fnm->flag & ffMULT)  != 0);
    const char *const name      = &fnm->opt[1];
    const char *      defName   = fnm->fn;
    if (defName == NULL)
    {
        defName = ftp2defnm(fnm->ftp);
    }
    // Since we are not (yet) using this for actual parsing, we don't need to
    // set any storage.
    options->addOption(
            FileNameOption(name).defaultBasename(defName).legacyType(fnm->ftp)
                .readWriteFlags(bRead, bWrite).required(!bOptional)
                .libraryFile(bLibrary).multiValue(bMultiple)
                .description(ftp2desc(fnm->ftp)));
}

/*! \brief
 * Converts a t_pargs option into an Options option.
 *
 * \param     options Options object to add the new option to.
 * \param[in] pa      t_pargs option to convert.
 *
 * \ingroup module_commandline
 */
void pargsToOptions(Options *options, t_pargs *pa)
{
    const bool        bHidden = is_hidden(pa);
    const char *const name    = &pa->option[1];
    const char *const desc    = (bHidden ? &pa->desc[6] : pa->desc);
    // Since we are not (yet) using this for actual parsing, we can take some
    // shortcuts and not set any storage where there is no direct
    // correspondence in the types.
    switch (pa->type)
    {
        case etINT:
            options->addOption(
                IntegerOption(name).store(pa->u.i)
                    .description(desc).hidden(bHidden));
            return;
        case etINT64:
            options->addOption(
                Int64Option(name).store(pa->u.is)
                    .description(desc).hidden(bHidden));
            return;
        case etREAL:
            options->addOption(
                RealOption(name).store(pa->u.r)
                    .description(desc).hidden(bHidden));
            return;
        case etTIME:
            options->addOption(
                RealOption(name).store(pa->u.r).timeValue()
                    .description(desc).hidden(bHidden));
            return;
        case etSTR:
        {
            const char *const defValue = (*pa->u.c != NULL ? *pa->u.c : "");
            options->addOption(
                    StringOption(name).defaultValue(defValue)
                        .description(desc).hidden(bHidden));
            return;
        }
        case etBOOL:
            options->addOption(
                BooleanOption(name).defaultValue(*pa->u.b)
                    .description(desc).hidden(bHidden));
            return;
        case etRVEC:
            options->addOption(
                RealOption(name).store(*pa->u.rv).vector()
                    .description(desc).hidden(bHidden));
            return;
        case etENUM:
            options->addOption(
                StringOption(name).defaultEnumIndex(nenum(pa->u.c) - 1)
                    .enumValueFromNullTerminatedArray(pa->u.c + 1)
                    .description(desc).hidden(bHidden));
            return;
    }
    GMX_THROW(NotImplementedError("Argument type not implemented"));
}

} // namespace

} // namespace gmx

gmx_bool parse_common_args(int *argc, char *argv[], unsigned long Flags,
                           int nfile, t_filenm fnm[], int npargs, t_pargs *pa,
                           int ndesc, const char **desc,
                           int nbugs, const char **bugs,
                           output_env_t *oenv)
{
    /* This array should match the order of the enum in oenv.h */
    const char *xvg_format[] = { NULL, "xmgrace", "xmgr", "none", NULL };
    /* This array should match the order of the enum in oenv.h */
    const char *time_units[] = {
        NULL, "fs", "ps", "ns", "us", "ms", "s",
        NULL
    };
    int         nicelevel = 0, debug_level = 0;
    char       *deffnm    = NULL;
    real        tbegin    = 0, tend = 0, tdelta = 0;
    gmx_bool    bView     = FALSE;

    t_pargs    *all_pa = NULL;

    t_pargs     nice_pa   = {
        "-nice", FALSE, etINT,   {&nicelevel},
        "Set the nicelevel"
    };
    t_pargs     deffnm_pa = {
        "-deffnm", FALSE, etSTR, {&deffnm},
        "Set the default filename for all file options"
    };
    t_pargs     begin_pa  = {
        "-b",    FALSE, etTIME,  {&tbegin},
        "First frame (%t) to read from trajectory"
    };
    t_pargs     end_pa    = {
        "-e",    FALSE, etTIME,  {&tend},
        "Last frame (%t) to read from trajectory"
    };
    t_pargs     dt_pa     = {
        "-dt",   FALSE, etTIME,  {&tdelta},
        "Only use frame when t MOD dt = first time (%t)"
    };
    t_pargs     view_pa   = {
        "-w",    FALSE, etBOOL,  {&bView},
        "View output [TT].xvg[tt], [TT].xpm[tt], [TT].eps[tt] and [TT].pdb[tt] files"
    };
    t_pargs     xvg_pa    = {
        "-xvg",  FALSE, etENUM,  {xvg_format},
        "xvg plot formatting"
    };
    t_pargs     time_pa   = {
        "-tu",   FALSE, etENUM,  {time_units},
        "Time unit"
    };
    /* Maximum number of extra arguments */
#define EXTRA_PA 16

    t_pargs  pca_pa[] = {
        { "-debug", FALSE, etINT, {&debug_level},
          "HIDDENWrite file with debug information, 1: short, 2: also x and f" },
    };
#define NPCA_PA asize(pca_pa)
    gmx_bool bXvgr;
    int      i, j, k, npall, max_pa;

    // Handle the flags argument, which is a bit field
    // The FF macro returns whether or not the bit is set
#define FF(arg) ((Flags & arg) == arg)

    /* Check for double arguments */
    for (i = 1; (i < *argc); i++)
    {
        if (argv[i] && (strlen(argv[i]) > 1) && (!std::isdigit(argv[i][1])))
        {
            for (j = i+1; (j < *argc); j++)
            {
                if ( (argv[i][0] == '-') && (argv[j][0] == '-') &&
                     (strcmp(argv[i], argv[j]) == 0) )
                {
                    if (FF(PCA_NOEXIT_ON_ARGS))
                    {
                        fprintf(stderr, "Double command line argument %s\n",
                                argv[i]);
                    }
                    else
                    {
                        gmx_fatal(FARGS, "Double command line argument %s\n",
                                  argv[i]);
                    }
                }
            }
        }
    }
    debug_gmx();

    /* Check ALL the flags ... */
    max_pa = NPCA_PA + EXTRA_PA + npargs+1;
    snew(all_pa, max_pa);

    for (i = npall = 0; (i < static_cast<int>(NPCA_PA)); i++)
    {
        npall = add_parg(npall, all_pa, &(pca_pa[i]));
    }

    if (FF(PCA_BE_NICE))
    {
        nicelevel = 19;
    }
    npall = add_parg(npall, all_pa, &nice_pa);

    if (FF(PCA_CAN_SET_DEFFNM))
    {
        npall = add_parg(npall, all_pa, &deffnm_pa);
    }
    if (FF(PCA_CAN_BEGIN))
    {
        npall = add_parg(npall, all_pa, &begin_pa);
    }
    if (FF(PCA_CAN_END))
    {
        npall = add_parg(npall, all_pa, &end_pa);
    }
    if (FF(PCA_CAN_DT))
    {
        npall = add_parg(npall, all_pa, &dt_pa);
    }
    if (FF(PCA_TIME_UNIT))
    {
        npall = add_parg(npall, all_pa, &time_pa);
    }
    if (FF(PCA_CAN_VIEW))
    {
        npall = add_parg(npall, all_pa, &view_pa);
    }

    bXvgr = FALSE;
    for (i = 0; (i < nfile); i++)
    {
        bXvgr = bXvgr ||  (fnm[i].ftp == efXVG);
    }
    if (bXvgr)
    {
        npall = add_parg(npall, all_pa, &xvg_pa);
    }

    /* Now append the program specific arguments */
    for (i = 0; (i < npargs); i++)
    {
        npall = add_parg(npall, all_pa, &(pa[i]));
    }

    /* set etENUM options to default */
    for (i = 0; (i < npall); i++)
    {
        if (all_pa[i].type == etENUM)
        {
            all_pa[i].u.c[0] = all_pa[i].u.c[1];
        }
    }
    set_default_time_unit(time_units, FF(PCA_TIME_UNIT));
    set_default_xvg_format(xvg_format);

    /* Now parse all the command-line options */
    get_pargs(argc, argv, npall, all_pa);

    /* set program name, command line, and default values for output options */
    output_env_init(oenv, gmx::getProgramContext(), (time_unit_t)nenum(time_units), bView,
                    (xvg_format_t)nenum(xvg_format), 0, debug_level);

    /* Parse the file args */
    parse_file_args(argc, argv, nfile, fnm, deffnm, !FF(PCA_NOT_READ_NODE));

    /* Open the debug file */
    if (debug_level > 0)
    {
        char buf[256];

        if (gmx_mpi_initialized())
        {
            sprintf(buf, "%s%d.debug", output_env_get_short_program_name(*oenv),
                    gmx_node_rank());
        }
        else
        {
            sprintf(buf, "%s.debug", output_env_get_short_program_name(*oenv));
        }

        init_debug(debug_level, buf);
        fprintf(stderr, "Opening debug file %s (src code file %s, line %d)\n",
                buf, __FILE__, __LINE__);
    }

    /* Now copy the results back... */
    for (i = 0, k = npall-npargs; (i < npargs); i++, k++)
    {
        memcpy(&(pa[i]), &(all_pa[k]), (size_t)sizeof(pa[i]));
    }

    bool bExit = false;
    try
    {
        const gmx::CommandLineHelpContext *context =
            gmx::GlobalCommandLineHelpContext::get();
        bExit = (context != NULL);
        if (context != NULL && !(FF(PCA_QUIET)))
        {
            gmx::Options options(NULL, NULL);
            options.setDescription(gmx::constArrayRefFromArray(desc, ndesc));
            for (i = 0; i < nfile; i++)
            {
                gmx::filenmToOptions(&options, &fnm[i]);
            }
            for (i = 0; i < npall; i++)
            {
                gmx::pargsToOptions(&options, &all_pa[i]);
            }
            gmx::CommandLineHelpWriter(options)
                .setShowDescriptions(true)
                .setTimeUnitString(output_env_get_time_unit(*oenv))
                .setKnownIssues(gmx::constArrayRefFromArray(bugs, nbugs))
                .writeHelp(*context);
        }
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    /* Set the nice level */
#if defined(HAVE_UNISTD_H) && !defined(__MINGW32__)
#ifndef GMX_NO_NICE
    /* The some system, e.g. the catamount kernel on cray xt3 do not have nice(2). */
    if (nicelevel != 0 && !bExit)
    {
        static gmx_bool            nice_set   = FALSE; /* only set it once */
        static tMPI_Thread_mutex_t init_mutex = TMPI_THREAD_MUTEX_INITIALIZER;
        tMPI_Thread_mutex_lock(&init_mutex);
        if (!nice_set)
        {
            if (nice(nicelevel) == -1)
            {
                /* Do nothing, but use the return value to avoid warnings. */
            }
            nice_set = TRUE;
        }
        tMPI_Thread_mutex_unlock(&init_mutex);
    }
#endif
#endif

    /* convert time options, must be done after printing! */

    for (i = 0; i < npall; i++)
    {
        if (all_pa[i].type == etTIME && all_pa[i].bSet)
        {
            *all_pa[i].u.r *= output_env_get_time_invfactor(*oenv);
        }
    }

    /* Extract Time info from arguments */
    if (FF(PCA_CAN_BEGIN) && opt2parg_bSet("-b", npall, all_pa))
    {
        setTimeValue(TBEGIN, opt2parg_real("-b", npall, all_pa));
    }

    if (FF(PCA_CAN_END) && opt2parg_bSet("-e", npall, all_pa))
    {
        setTimeValue(TEND, opt2parg_real("-e", npall, all_pa));
    }

    if (FF(PCA_CAN_DT) && opt2parg_bSet("-dt", npall, all_pa))
    {
        setTimeValue(TDELTA, opt2parg_real("-dt", npall, all_pa));
    }

    /* clear memory */
    sfree(all_pa);

    if (!FF(PCA_NOEXIT_ON_ARGS))
    {
        if (*argc > 1)
        {
            gmx_cmd(argv[1]);
        }
    }
    return !bExit;
#undef FF
}
