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
#include "gmxpre.h"

#include "pargs.h"

#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <list>

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/commandline/cmdlinehelpwriter.h"
#include "gromacs/commandline/cmdlineparser.h"
#include "gromacs/fileio/timecontrol.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/filenameoptionmanager.h"
#include "gromacs/options/options.h"
#include "gromacs/options/timeunitmanager.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

/* The source code in this file should be thread-safe.
      Please keep it that way. */

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

namespace gmx
{

namespace
{

/*! \brief
 * Returns the index of the default xvg format.
 *
 * \ingroup module_commandline
 */
int getDefaultXvgFormat(gmx::ConstArrayRef<const char *> xvgFormats)
{
    const char *const select = getenv("GMX_VIEW_XVG");
    if (select != NULL)
    {
        ConstArrayRef<const char *>::const_iterator i =
            std::find(xvgFormats.begin(), xvgFormats.end(), std::string(select));
        if (i != xvgFormats.end())
        {
            return i - xvgFormats.begin();
        }
        else
        {
            return exvgNONE - 1;
        }
    }
    /* The default is the first option */
    return 0;
}

/*! \brief
 * Conversion helper between t_pargs/t_filenm and Options.
 *
 * This class holds the necessary mapping between the old C structures and
 * the new C++ options to allow copying values back after parsing for cases
 * where the C++ options do not directly provide the type of value required for
 * the C structures.
 *
 * \ingroup module_commandline
 */
class OptionsAdapter
{
    public:
        /*! \brief
         * Initializes the adapter to convert from a specified command line.
         *
         * The command line is required, because t_pargs wants to return
         * strings by reference to the original command line.
         * OptionsAdapter creates a copy of the `argv` array (but not the
         * strings) to make this possible, even if the parser removes
         * options it has recognized.
         */
        OptionsAdapter(int argc, const char *const argv[])
            : argv_(argv, argv + argc)
        {
        }

        /*! \brief
         * Converts a t_filenm option into an Options option.
         *
         * \param options Options object to add the new option to.
         * \param fnm     t_filenm option to convert.
         */
        void filenmToOptions(Options *options, t_filenm *fnm);
        /*! \brief
         * Converts a t_pargs option into an Options option.
         *
         * \param     options Options object to add the new option to.
         * \param     pa      t_pargs option to convert.
         */
        void pargsToOptions(Options *options, t_pargs *pa);

        /*! \brief
         * Copies values back from options to t_pargs/t_filenm.
         */
        void copyValues(bool bReadNode);

    private:
        struct FileNameData
        {
            //! Creates a conversion helper for a given `t_filenm` struct.
            explicit FileNameData(t_filenm *fnm) : fnm(fnm), optionInfo(NULL)
            {
            }

            //! t_filenm structure to receive the final values.
            t_filenm                 *fnm;
            //! Option info object for the created FileNameOption.
            FileNameOptionInfo       *optionInfo;
            //! Value storage for the created FileNameOption.
            std::vector<std::string>  values;
        };
        struct ProgramArgData
        {
            //! Creates a conversion helper for a given `t_pargs` struct.
            explicit ProgramArgData(t_pargs *pa)
                : pa(pa), optionInfo(NULL), enumIndex(0), boolValue(false)
            {
            }

            //! t_pargs structure to receive the final values.
            t_pargs                 *pa;
            //! Option info object for the created option.
            OptionInfo              *optionInfo;
            //! Value storage for a non-enum StringOption (unused for other types).
            std::string              stringValue;
            //! Value storage for an enum option (unused for other types).
            int                      enumIndex;
            //! Value storage for a BooleanOption (unused for other types).
            bool                     boolValue;
        };

        std::vector<const char *>    argv_;
        // These are lists instead of vectors to avoid relocating existing
        // objects in case the container is reallocated (the Options object
        // contains pointes to members of the objects, which would get
        // invalidated).
        std::list<FileNameData>      fileNameOptions_;
        std::list<ProgramArgData>    programArgs_;

        GMX_DISALLOW_COPY_AND_ASSIGN(OptionsAdapter);
};

void OptionsAdapter::filenmToOptions(Options *options, t_filenm *fnm)
{
    if (fnm->opt == NULL)
    {
        // Existing code may use opt2fn() instead of ftp2fn() for
        // options that use the default option name, so we need to
        // keep the old behavior instead of fixing opt2fn().
        // TODO: Check that this is not the case, remove this, and make
        // opt2*() work even if fnm->opt is NULL for some options.
        fnm->opt = ftp2defopt(fnm->ftp);
    }
    const bool        bRead     = ((fnm->flag & ffREAD)  != 0);
    const bool        bWrite    = ((fnm->flag & ffWRITE) != 0);
    const bool        bOptional = ((fnm->flag & ffOPT)   != 0);
    const bool        bLibrary  = ((fnm->flag & ffLIB)   != 0);
    const bool        bMultiple = ((fnm->flag & ffMULT)  != 0);
    const bool        bMissing  = ((fnm->flag & ffALLOW_MISSING) != 0);
    const char *const name      = &fnm->opt[1];
    const char *      defName   = fnm->fn;
    int               defType   = -1;
    if (defName == NULL)
    {
        defName = ftp2defnm(fnm->ftp);
    }
    else if (Path::hasExtension(defName))
    {
        defType = fn2ftp(defName);
        GMX_RELEASE_ASSERT(defType != efNR,
                           "File name option specifies an invalid extension");
    }
    fileNameOptions_.push_back(FileNameData(fnm));
    FileNameData &data = fileNameOptions_.back();
    data.optionInfo = options->addOption(
                FileNameOption(name).storeVector(&data.values)
                    .defaultBasename(defName).defaultType(defType)
                    .legacyType(fnm->ftp).legacyOptionalBehavior()
                    .readWriteFlags(bRead, bWrite).required(!bOptional)
                    .libraryFile(bLibrary).multiValue(bMultiple)
                    .allowMissing(bMissing)
                    .description(ftp2desc(fnm->ftp)));
}

void OptionsAdapter::pargsToOptions(Options *options, t_pargs *pa)
{
    const bool        bHidden = startsWith(pa->desc, "HIDDEN");
    const char *const name    = &pa->option[1];
    const char *const desc    = (bHidden ? &pa->desc[6] : pa->desc);
    programArgs_.push_back(ProgramArgData(pa));
    ProgramArgData   &data = programArgs_.back();
    switch (pa->type)
    {
        case etINT:
            data.optionInfo = options->addOption(
                        IntegerOption(name).store(pa->u.i)
                            .description(desc).hidden(bHidden));
            return;
        case etINT64:
            data.optionInfo = options->addOption(
                        Int64Option(name).store(pa->u.is)
                            .description(desc).hidden(bHidden));
            return;
        case etREAL:
            data.optionInfo = options->addOption(
                        RealOption(name).store(pa->u.r)
                            .description(desc).hidden(bHidden));
            return;
        case etTIME:
            data.optionInfo = options->addOption(
                        RealOption(name).store(pa->u.r).timeValue()
                            .description(desc).hidden(bHidden));
            return;
        case etSTR:
        {
            const char *const defValue = (*pa->u.c != NULL ? *pa->u.c : "");
            data.optionInfo = options->addOption(
                        StringOption(name).store(&data.stringValue)
                            .defaultValue(defValue)
                            .description(desc).hidden(bHidden));
            return;
        }
        case etBOOL:
            data.optionInfo = options->addOption(
                        BooleanOption(name).store(&data.boolValue)
                            .defaultValue(*pa->u.b)
                            .description(desc).hidden(bHidden));
            return;
        case etRVEC:
            data.optionInfo = options->addOption(
                        RealOption(name).store(*pa->u.rv).vector()
                            .description(desc).hidden(bHidden));
            return;
        case etENUM:
        {
            const int defaultIndex = (pa->u.c[0] != NULL ? nenum(pa->u.c) - 1 : 0);
            data.optionInfo = options->addOption(
                        StringOption(name).storeEnumIndex(&data.enumIndex)
                            .defaultEnumIndex(defaultIndex)
                            .enumValueFromNullTerminatedArray(pa->u.c + 1)
                            .description(desc).hidden(bHidden));
            return;
        }
    }
    GMX_THROW(NotImplementedError("Argument type not implemented"));
}

void OptionsAdapter::copyValues(bool bReadNode)
{
    std::list<FileNameData>::const_iterator file;
    for (file = fileNameOptions_.begin(); file != fileNameOptions_.end(); ++file)
    {
        if (!bReadNode && (file->fnm->flag & ffREAD))
        {
            continue;
        }
        if (file->optionInfo->isSet())
        {
            file->fnm->flag |= ffSET;
        }
        file->fnm->nfiles = file->values.size();
        snew(file->fnm->fns, file->fnm->nfiles);
        for (int i = 0; i < file->fnm->nfiles; ++i)
        {
            file->fnm->fns[i] = gmx_strdup(file->values[i].c_str());
        }
    }
    std::list<ProgramArgData>::const_iterator arg;
    for (arg = programArgs_.begin(); arg != programArgs_.end(); ++arg)
    {
        arg->pa->bSet = arg->optionInfo->isSet();
        switch (arg->pa->type)
        {
            case etSTR:
            {
                if (arg->pa->bSet)
                {
                    std::vector<const char *>::const_iterator pos =
                        std::find(argv_.begin(), argv_.end(), arg->stringValue);
                    GMX_RELEASE_ASSERT(pos != argv_.end(),
                                       "String argument got a value not in argv");
                    *arg->pa->u.c = *pos;
                }
                break;
            }
            case etBOOL:
                *arg->pa->u.b = arg->boolValue;
                break;
            case etENUM:
                *arg->pa->u.c = arg->pa->u.c[arg->enumIndex + 1];
                break;
            default:
                // For other types, there is nothing type-specific to do.
                break;
        }
    }
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
    const char *const xvg_formats[] = { "xmgrace", "xmgr", "none" };

    // Handle the flags argument, which is a bit field
    // The FF macro returns whether or not the bit is set
#define FF(arg) ((Flags & arg) == arg)

    try
    {
        double                     tbegin    = 0.0, tend = 0.0, tdelta = 0.0;
        bool                       bView     = false;
        int                        xvgFormat = 0;
        gmx::TimeUnitManager       timeUnitManager;
        gmx::OptionsAdapter        adapter(*argc, argv);
        gmx::Options               options(NULL, NULL);
        gmx::FileNameOptionManager fileOptManager;

        fileOptManager.disableInputOptionChecking(
                FF(PCA_NOT_READ_NODE) || FF(PCA_DISABLE_INPUT_FILE_CHECKING));
        options.addManager(&fileOptManager);
        options.setDescription(gmx::constArrayRefFromArray<const char *>(desc, ndesc));

        if (FF(PCA_CAN_SET_DEFFNM))
        {
            fileOptManager.addDefaultFileNameOption(&options, "deffnm");
        }
        if (FF(PCA_CAN_BEGIN))
        {
            options.addOption(
                    gmx::DoubleOption("b").store(&tbegin).timeValue()
                        .description("First frame (%t) to read from trajectory"));
        }
        if (FF(PCA_CAN_END))
        {
            options.addOption(
                    gmx::DoubleOption("e").store(&tend).timeValue()
                        .description("Last frame (%t) to read from trajectory"));
        }
        if (FF(PCA_CAN_DT))
        {
            options.addOption(
                    gmx::DoubleOption("dt").store(&tdelta).timeValue()
                        .description("Only use frame when t MOD dt = first time (%t)"));
        }
        if (FF(PCA_TIME_UNIT))
        {
            timeUnitManager.setTimeUnitFromEnvironment();
            timeUnitManager.addTimeUnitOption(&options, "tu");
        }
        if (FF(PCA_CAN_VIEW))
        {
            options.addOption(
                    gmx::BooleanOption("w").store(&bView)
                        .description("View output [REF].xvg[ref], [REF].xpm[ref], "
                                     "[REF].eps[ref] and [REF].pdb[ref] files"));
        }

        bool bXvgr = false;
        for (int i = 0; i < nfile; i++)
        {
            bXvgr = bXvgr || (fnm[i].ftp == efXVG);
        }
        xvgFormat = gmx::getDefaultXvgFormat(xvg_formats);
        if (bXvgr)
        {
            options.addOption(
                    gmx::StringOption("xvg").enumValue(xvg_formats)
                        .storeEnumIndex(&xvgFormat)
                        .description("xvg plot formatting"));
        }

        /* Now append the program specific arguments */
        for (int i = 0; i < nfile; i++)
        {
            adapter.filenmToOptions(&options, &fnm[i]);
        }
        for (int i = 0; i < npargs; i++)
        {
            adapter.pargsToOptions(&options, &pa[i]);
        }

        const gmx::CommandLineHelpContext *context =
            gmx::GlobalCommandLineHelpContext::get();
        if (context != NULL)
        {
            GMX_RELEASE_ASSERT(gmx_node_rank() == 0,
                               "Help output should be handled higher up and "
                               "only get called only on the master rank");
            gmx::CommandLineHelpWriter(options)
                .setShowDescriptions(true)
                .setTimeUnitString(timeUnitManager.timeUnitAsString())
                .setKnownIssues(gmx::constArrayRefFromArray(bugs, nbugs))
                .writeHelp(*context);
            return FALSE;
        }

        /* Now parse all the command-line options */
        gmx::CommandLineParser(&options).skipUnknown(FF(PCA_NOEXIT_ON_ARGS))
            .parse(argc, argv);
        options.finish();

        /* set program name, command line, and default values for output options */
        output_env_init(oenv, gmx::getProgramContext(),
                        (time_unit_t)(timeUnitManager.timeUnit() + 1), bView,
                        (xvg_format_t)(xvgFormat + 1), 0);

        timeUnitManager.scaleTimeOptions(&options);

        /* Extract Time info from arguments */
        // TODO: Use OptionInfo objects instead of string constants
        if (FF(PCA_CAN_BEGIN) && options.isSet("b"))
        {
            setTimeValue(TBEGIN, tbegin);
        }
        if (FF(PCA_CAN_END) && options.isSet("e"))
        {
            setTimeValue(TEND, tend);
        }
        if (FF(PCA_CAN_DT) && options.isSet("dt"))
        {
            setTimeValue(TDELTA, tdelta);
        }

        adapter.copyValues(!FF(PCA_NOT_READ_NODE));

        return TRUE;
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
#undef FF
}
