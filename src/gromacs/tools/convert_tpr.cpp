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
#include "gmxpre.h"

#include "convert_tpr.h"

#include <cmath>

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trrio.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/random/seed.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

static void rangeCheck(int numberInIndexFile, int maxAtomNumber)
{
    if ((numberInIndexFile) >= (maxAtomNumber))
    {
        gmx_fatal(FARGS,
                  "Your index file contains atomnumbers (e.g. %d)\nthat are larger than the number "
                  "of atoms in the tpr file (%d)",
                  (numberInIndexFile),
                  (maxAtomNumber));
    }
}

static std::vector<bool> bKeepIt(int gnx, int natoms, int index[])
{
    std::vector<bool> b(natoms);

    for (int i = 0; (i < gnx); i++)
    {
        rangeCheck(index[i], natoms);
        b[index[i]] = TRUE;
    }

    return b;
}

static std::vector<int> invind(int gnx, int natoms, int index[])
{
    std::vector<int> inv(natoms);

    for (int i = 0; (i < gnx); i++)
    {
        rangeCheck(index[i], natoms);
        inv[index[i]] = i;
    }

    return inv;
}

static gmx::ListOfLists<int> reduce_listoflists(gmx::ArrayRef<const int>     invindex,
                                                const std::vector<bool>&     bKeep,
                                                const gmx::ListOfLists<int>& src,
                                                const char*                  name)
{
    gmx::ListOfLists<int> lists;

    std::vector<int> exclusionsForAtom;
    for (gmx::index i = 0; i < src.ssize(); i++)
    {
        if (bKeep[i])
        {
            exclusionsForAtom.clear();
            for (const int j : src[i])
            {
                if (bKeep[j])
                {
                    exclusionsForAtom.push_back(invindex[j]);
                }
            }
            lists.pushBack(exclusionsForAtom);
        }
    }

    fprintf(stderr,
            "Reduced block %8s from %6zu to %6zu index-, %6d to %6d a-entries\n",
            name,
            src.size(),
            lists.size(),
            src.numElements(),
            lists.numElements());

    return lists;
}

static void reduce_rvec(int gnx, const int index[], rvec vv[])
{
    rvec* ptr;
    int   i;

    snew(ptr, gnx);
    for (i = 0; (i < gnx); i++)
    {
        copy_rvec(vv[index[i]], ptr[i]);
    }
    for (i = 0; (i < gnx); i++)
    {
        copy_rvec(ptr[i], vv[i]);
    }
    sfree(ptr);
}

static void reduce_atom(int gnx, const int index[], t_atom atom[], char*** atomname, int* nres, t_resinfo* resinfo)
{
    t_atom*    ptr;
    char***    aname;
    t_resinfo* rinfo;
    int        i, nr;

    snew(ptr, gnx);
    snew(aname, gnx);
    snew(rinfo, atom[index[gnx - 1]].resind + 1);
    for (i = 0; (i < gnx); i++)
    {
        ptr[i]   = atom[index[i]];
        aname[i] = atomname[index[i]];
    }
    nr = -1;
    for (i = 0; (i < gnx); i++)
    {
        atom[i]     = ptr[i];
        atomname[i] = aname[i];
        if ((i == 0) || (atom[i].resind != atom[i - 1].resind))
        {
            nr++;
            rinfo[nr] = resinfo[atom[i].resind];
        }
        atom[i].resind = nr;
    }
    nr++;
    for (i = 0; (i < nr); i++)
    {
        resinfo[i] = rinfo[i];
    }
    *nres = nr;

    sfree(aname);
    sfree(ptr);
    sfree(rinfo);
}

static void reduce_ilist(gmx::ArrayRef<const int> invindex,
                         const std::vector<bool>& bKeep,
                         InteractionList*         il,
                         int                      nratoms,
                         const char*              name)
{
    if (!il->empty())
    {
        std::vector<int> newAtoms(nratoms);
        InteractionList  ilReduced;
        for (int i = 0; i < il->size(); i += nratoms + 1)
        {
            bool bB = true;
            for (int j = 0; j < nratoms; j++)
            {
                bB = bB && bKeep[il->iatoms[i + 1 + j]];
            }
            if (bB)
            {
                for (int j = 0; j < nratoms; j++)
                {
                    newAtoms[j] = invindex[il->iatoms[i + 1 + j]];
                }
                ilReduced.push_back(il->iatoms[i], nratoms, newAtoms.data());
            }
        }
        fprintf(stderr,
                "Reduced ilist %8s from %6d to %6d entries\n",
                name,
                il->size() / (nratoms + 1),
                ilReduced.size() / (nratoms + 1));

        *il = std::move(ilReduced);
    }
}

static void reduce_topology_x(int gnx, int index[], gmx_mtop_t* mtop, rvec x[], rvec v[])
{
    gmx_localtop_t top(mtop->ffparams);
    gmx_mtop_generate_local_top(*mtop, &top, false);
    t_atoms atoms = gmx_mtop_global_atoms(mtop);

    const std::vector<bool> bKeep    = bKeepIt(gnx, atoms.nr, index);
    const std::vector<int>  invindex = invind(gnx, atoms.nr, index);

    reduce_rvec(gnx, index, x);
    reduce_rvec(gnx, index, v);
    reduce_atom(gnx, index, atoms.atom, atoms.atomname, &(atoms.nres), atoms.resinfo);

    for (int i = 0; (i < F_NRE); i++)
    {
        reduce_ilist(invindex,
                     bKeep,
                     &(top.idef.il[i]),
                     interaction_function[i].nratoms,
                     interaction_function[i].name);
    }

    atoms.nr = gnx;

    mtop->moltype.resize(1);
    mtop->moltype[0].name  = mtop->name;
    mtop->moltype[0].atoms = atoms;
    mtop->moltype[0].excls = reduce_listoflists(invindex, bKeep, top.excls, "excls");
    for (int i = 0; i < F_NRE; i++)
    {
        mtop->moltype[0].ilist[i] = std::move(top.idef.il[i]);
    }

    mtop->molblock.resize(1);
    mtop->molblock[0].type = 0;
    mtop->molblock[0].nmol = 1;

    mtop->natoms = atoms.nr;
}

static void zeroq(const int index[], gmx_mtop_t* mtop)
{
    for (gmx_moltype_t& moltype : mtop->moltype)
    {
        for (int i = 0; i < moltype.atoms.nr; i++)
        {
            moltype.atoms.atom[index[i]].q  = 0;
            moltype.atoms.atom[index[i]].qB = 0;
        }
    }
}

namespace gmx
{

namespace
{

class ConvertTpr : public ICommandLineOptionsModule
{
public:
    ConvertTpr() {}

    // From ICommandLineOptionsModule
    void init(CommandLineModuleSettings* /*settings*/) override {}
    void initOptions(IOptionsContainer* options, ICommandLineOptionsModuleSettings* settings) override;
    void optionsFinished() override;
    int  run() override;

private:
    //! Name of input tpr file.
    std::string inputTprFileName_;
    //! Name of input index file.
    std::string inputIndexFileName_;
    //! Name of output tpr file.
    std::string outputTprFileName_;
    //! If we have read in an index file.
    bool haveReadIndexFile_ = false;
    //! Time to extend simulation by.
    real extendTime_ = 0;
    //! If the option to extend simulation time is set.
    bool extendTimeIsSet_ = false;
    //! Final run time value.
    real runToMaxTime_ = 0;
    //! If the option to run simulation until specified time is set.
    bool runToMaxTimeIsSet_ = false;
    //! Maximum number of steps to run.
    int64_t maxSteps_ = 0;
    //! If the option to use maximumstep number is set.
    bool maxStepsIsSet_ = false;
    //! If the option to zero charge is set.
    bool zeroQIsSet_ = false;
};

void ConvertTpr::initOptions(IOptionsContainer* options, ICommandLineOptionsModuleSettings* settings)
{
    std::vector<const char*> desc = {
        "[THISMODULE] can edit run input files in three ways.[PAR]",
        "[BB]1.[bb] by modifying the number of steps in a run input file",
        "with options [TT]-extend[tt], [TT]-until[tt] or [TT]-nsteps[tt]",
        "(nsteps=-1 means unlimited number of steps)[PAR]",
        "[BB]2.[bb] by creating a [REF].tpx[ref] file for a subset of your original",
        "tpx file, which is useful when you want to remove the solvent from",
        "your [REF].tpx[ref] file, or when you want to make e.g. a pure C[GRK]alpha[grk] ",
        "[REF].tpx[ref] file.",
        "Note that you may need to use [TT]-nsteps -1[tt] (or similar) to get",
        "this to work.",
        "[BB]WARNING: this [REF].tpx[ref] file is not fully functional[bb].[PAR]",
        "[BB]3.[bb] by setting the charges of a specified group",
        "to zero. This is useful when doing free energy estimates",
        "using the LIE (Linear Interaction Energy) method."
    };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("s")
                               .filetype(eftTopology)
                               .inputFile()
                               .required()
                               .store(&inputTprFileName_)
                               .defaultBasename("topol")
                               .description("Run input file to modify"));
    options->addOption(FileNameOption("n")
                               .filetype(eftIndex)
                               .inputFile()
                               .store(&inputIndexFileName_)
                               .storeIsSet(&haveReadIndexFile_)
                               .defaultBasename("index")
                               .description("File containing additional index groups"));
    options->addOption(FileNameOption("o")
                               .filetype(eftTopology)
                               .outputFile()
                               .store(&outputTprFileName_)
                               .defaultBasename("tprout")
                               .description("Generated modified run input file"));
    options->addOption(RealOption("extend")
                               .store(&extendTime_)
                               .storeIsSet(&extendTimeIsSet_)
                               .timeValue()
                               .description("Extend runtime by this amount (ps)"));
    options->addOption(RealOption("until")
                               .store(&runToMaxTime_)
                               .storeIsSet(&runToMaxTimeIsSet_)
                               .timeValue()
                               .description("Extend runtime until this ending time (ps)"));
    options->addOption(
            Int64Option("nsteps").store(&maxSteps_).storeIsSet(&maxStepsIsSet_).description("Change the number of steps"));
    options->addOption(
            BooleanOption("zeroq").store(&zeroQIsSet_).description("Set the charges of a group (from the index) to zero"));
}

void ConvertTpr::optionsFinished() {}

int ConvertTpr::run()
{
    gmx_mtop_t mtop;
    t_atoms    atoms;
    t_state    state;
    char       buf[200], buf2[200];

    fprintf(stderr, "Reading toplogy and stuff from %s\n", inputTprFileName_.c_str());

    t_inputrec  irInstance;
    t_inputrec* ir = &irInstance;
    read_tpx_state(inputTprFileName_.c_str(), ir, &state, &mtop);
    int64_t currentMaxStep    = ir->init_step;
    double  currentRunTime    = ir->init_step * ir->delta_t + ir->init_t;
    real    currentMaxRunTime = 0.0;

    if (maxStepsIsSet_)
    {
        fprintf(stderr, "Setting nsteps to %s\n", gmx_step_str(maxSteps_, buf));
        ir->nsteps = maxSteps_;
    }
    else
    {
        /* Determine total number of steps remaining */
        if (extendTimeIsSet_)
        {
            ir->nsteps = ir->nsteps - (currentMaxStep - ir->init_step)
                         + gmx::roundToInt64(extendTime_ / ir->delta_t);
            printf("Extending remaining runtime of by %g ps (now %s steps)\n",
                   extendTime_,
                   gmx_step_str(ir->nsteps, buf));
        }
        else if (runToMaxTimeIsSet_)
        {
            printf("nsteps = %s, run_step = %s, current_t = %g, until = %g\n",
                   gmx_step_str(ir->nsteps, buf),
                   gmx_step_str(currentMaxStep, buf2),
                   currentRunTime,
                   runToMaxTime_);
            ir->nsteps = gmx::roundToInt64((currentMaxRunTime - currentRunTime) / ir->delta_t);
            printf("Extending remaining runtime until %g ps (now %s steps)\n",
                   currentMaxRunTime,
                   gmx_step_str(ir->nsteps, buf));
        }
        else
        {
            ir->nsteps -= currentMaxStep - ir->init_step;
            /* Print message */
            printf("%s steps (%g ps) remaining from first run.\n",
                   gmx_step_str(ir->nsteps, buf),
                   ir->nsteps * ir->delta_t);
        }
    }

    if (maxStepsIsSet_ || zeroQIsSet_ || (ir->nsteps > 0))
    {
        ir->init_step = currentMaxStep;

        if (haveReadIndexFile_ || !(maxStepsIsSet_ || extendTimeIsSet_ || runToMaxTimeIsSet_))
        {
            atoms         = gmx_mtop_global_atoms(&mtop);
            int   gnx     = 0;
            int*  index   = nullptr;
            char* grpname = nullptr;
            get_index(&atoms, inputIndexFileName_.c_str(), 1, &gnx, &index, &grpname);
            bool bSel = false;
            if (!zeroQIsSet_)
            {
                bSel = (gnx != state.natoms);
                for (int i = 0; ((i < gnx) && (!bSel)); i++)
                {
                    bSel = (i != index[i]);
                }
            }
            else
            {
                bSel = false;
            }
            if (bSel)
            {
                fprintf(stderr,
                        "Will write subset %s of original tpx containing %d "
                        "atoms\n",
                        grpname,
                        gnx);
                reduce_topology_x(gnx, index, &mtop, state.x.rvec_array(), state.v.rvec_array());
                state.natoms = gnx;
            }
            else if (zeroQIsSet_)
            {
                zeroq(index, &mtop);
                fprintf(stderr, "Zero-ing charges for group %s\n", grpname);
            }
            else
            {
                fprintf(stderr, "Will write full tpx file (no selection)\n");
            }
        }

        double stateTime = ir->init_t + ir->init_step * ir->delta_t;
        sprintf(buf,
                "Writing statusfile with starting step %s%s and length %s%s steps...\n",
                "%10",
                PRId64,
                "%10",
                PRId64);
        fprintf(stderr, buf, ir->init_step, ir->nsteps);
        fprintf(stderr,
                "                                 time %10.3f and length %10.3f ps\n",
                stateTime,
                ir->nsteps * ir->delta_t);
        write_tpx_state(outputTprFileName_.c_str(), ir, &state, &mtop);
    }
    else
    {
        printf("You've simulated long enough. Not writing tpr file\n");
    }

    return 0;
}

} // namespace

const char ConvertTprInfo::name[]             = "convert-tpr";
const char ConvertTprInfo::shortDescription[] = "Make a modifed run-input file";
ICommandLineOptionsModulePointer ConvertTprInfo::create()
{
    return ICommandLineOptionsModulePointer(std::make_unique<ConvertTpr>());
}

} // namespace gmx
