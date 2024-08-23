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
/*! \internal \file
 * \brief Implements gmx dump utility.
 *
 * \ingroup module_tools
 */
#include "gmxpre.h"

#include "dump.h"

#include "config.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/mtxio.h"
#include "gromacs/fileio/tngio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trrio.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/gmxpreprocess/gmxcpp.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/mdrun/mdmodules.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/txtdump.h"

namespace gmx
{

namespace
{

//! Dump a TPR file
void list_tpr(const char* fn,
              gmx_bool    bShowNumbers,
              gmx_bool    bShowParameters,
              const char* mdpfn,
              gmx_bool    bSysTop,
              gmx_bool    bOriginalInputrec)
{
    FILE*      gp;
    int        indent, atot;
    t_state    state;
    gmx_mtop_t mtop;
    t_topology top;

    TpxFileHeader tpx = readTpxHeader(fn, true);
    t_inputrec    ir;

    read_tpx_state(fn, tpx.bIr ? &ir : nullptr, &state, tpx.bTop ? &mtop : nullptr);
    if (tpx.bIr && !bOriginalInputrec)
    {
        MDModules().adjustInputrecBasedOnModules(&ir);
    }

    if (mdpfn && tpx.bIr)
    {
        gp = gmx_fio_fopen(mdpfn, "w");
        pr_inputrec(gp, 0, nullptr, &ir, TRUE);
        gmx_fio_fclose(gp);
    }

    if (!mdpfn)
    {
        if (bSysTop)
        {
            top = gmx_mtop_t_to_t_topology(&mtop, false);
        }

        if (available(stdout, &tpx, 0, fn))
        {
            indent = 0;
            pr_title(stdout, indent, fn);
            pr_inputrec(stdout, 0, "inputrec", tpx.bIr ? &ir : nullptr, FALSE);

            pr_tpxheader(stdout, indent, "header", &(tpx));

            if (!bSysTop)
            {
                pr_mtop(stdout, indent, "topology", &(mtop), bShowNumbers, bShowParameters);
            }
            else
            {
                pr_top(stdout, indent, "topology", &(top), bShowNumbers, bShowParameters);
            }

            pr_rvecs(stdout, indent, "box", tpx.bBox ? state.box : nullptr, DIM);
            pr_rvecs(stdout, indent, "box_rel", tpx.bBox ? state.box_rel : nullptr, DIM);
            pr_rvecs(stdout, indent, "boxv", tpx.bBox ? state.boxv : nullptr, DIM);
            pr_rvecs(stdout, indent, "pres_prev", tpx.bBox ? state.pres_prev : nullptr, DIM);
            pr_rvecs(stdout, indent, "svir_prev", tpx.bBox ? state.svir_prev : nullptr, DIM);
            pr_rvecs(stdout, indent, "fvir_prev", tpx.bBox ? state.fvir_prev : nullptr, DIM);
            /* leave nosehoover_xi in for now to match the tpr version */
            pr_doubles(stdout, indent, "nosehoover_xi", state.nosehoover_xi.data(), state.ngtc);
            /*pr_doubles(stdout,indent,"nosehoover_vxi",state.nosehoover_vxi,state.ngtc);*/
            /*pr_doubles(stdout,indent,"therm_integral",state.therm_integral,state.ngtc);*/
            pr_rvecs(stdout, indent, "x", tpx.bX ? state.x.rvec_array() : nullptr, state.numAtoms());
            pr_rvecs(stdout, indent, "v", tpx.bV ? state.v.rvec_array() : nullptr, state.numAtoms());
        }

        const SimulationGroups& groups = mtop.groups;

        gmx::EnumerationArray<SimulationAtomGroupType, std::vector<int>> gcount;
        for (auto group : keysOf(gcount))
        {
            gcount[group].resize(groups.groups[group].size());
        }

        for (int i = 0; (i < mtop.natoms); i++)
        {
            for (auto group : keysOf(gcount))
            {
                gcount[group][getGroupType(groups, group, i)]++;
            }
        }
        printf("Group statistics\n");
        for (auto group : keysOf(gcount))
        {
            atot = 0;
            printf("%-12s: ", shortName(group));
            for (const auto& entry : gcount[group])
            {
                printf("  %5d", entry);
                atot += entry;
            }
            printf("  (total %d atoms)\n", atot);
        }
    }
}

//! Dump a topology file
void list_top(const char* fn)
{
    int status, done;
    // Legacy string length macro
    char      buf[STRLEN];
    gmx_cpp_t handle;
    char*     cppopts[] = { nullptr };

    status = cpp_open_file(fn, &handle, cppopts);
    if (status != 0)
    {
        gmx_fatal(FARGS, "%s", cpp_error(&handle, status));
    }
    do
    {
        status = cpp_read_line(&handle, STRLEN, buf);
        done   = static_cast<int>(status == eCPP_EOF);
        if (!done)
        {
            if (status != eCPP_OK)
            {
                gmx_fatal(FARGS, "%s", cpp_error(&handle, status));
            }
            else
            {
                printf("%s\n", buf);
            }
        }
    } while (done == 0);
    status = cpp_close_file(&handle);
    if (status != eCPP_OK)
    {
        gmx_fatal(FARGS, "%s", cpp_error(&handle, status));
    }
}

//! Dump a TRR file
void list_trr(const char* fn)
{
    t_fileio*        fpread;
    int              nframe, indent;
    char             buf[256];
    rvec *           x, *v, *f;
    matrix           box;
    gmx_trr_header_t trrheader;
    gmx_bool         bOK;

    fpread = gmx_trr_open(fn, "r");

    nframe = 0;
    while (gmx_trr_read_frame_header(fpread, &trrheader, &bOK))
    {
        snew(x, trrheader.natoms);
        snew(v, trrheader.natoms);
        snew(f, trrheader.natoms);
        if (gmx_trr_read_frame_data(fpread,
                                    &trrheader,
                                    trrheader.box_size ? box : nullptr,
                                    trrheader.x_size ? x : nullptr,
                                    trrheader.v_size ? v : nullptr,
                                    trrheader.f_size ? f : nullptr))
        {
            sprintf(buf, "%s frame %d", fn, nframe);
            indent = 0;
            indent = pr_title(stdout, indent, buf);
            pr_indent(stdout, indent);
            fprintf(stdout,
                    "natoms=%10d  step=%10" PRId64 "  time=%12.7e  lambda=%10g\n",
                    trrheader.natoms,
                    trrheader.step,
                    trrheader.t,
                    trrheader.lambda);
            if (trrheader.box_size)
            {
                pr_rvecs(stdout, indent, "box", box, DIM);
            }
            if (trrheader.x_size)
            {
                pr_rvecs(stdout, indent, "x", x, trrheader.natoms);
            }
            if (trrheader.v_size)
            {
                pr_rvecs(stdout, indent, "v", v, trrheader.natoms);
            }
            if (trrheader.f_size)
            {
                pr_rvecs(stdout, indent, "f", f, trrheader.natoms);
            }
        }
        else
        {
            fprintf(stderr, "\nWARNING: Incomplete frame: nr %d, t=%g\n", nframe, trrheader.t);
        }

        sfree(x);
        sfree(v);
        sfree(f);
        nframe++;
    }
    if (!bOK)
    {
        fprintf(stderr, "\nWARNING: Incomplete frame header: nr %d, t=%g\n", nframe, trrheader.t);
    }
    gmx_trr_close(fpread);
}

//! Dump an xtc file
void list_xtc(const char* fn)
{
    t_fileio* xd;
    int       indent;
    char      buf[256];
    rvec*     x;
    matrix    box;
    int       nframe, natoms;
    int64_t   step;
    real      prec, time;
    gmx_bool  bOK;

    xd = open_xtc(fn, "r");
    read_first_xtc(xd, &natoms, &step, &time, box, &x, &prec, &bOK);

    nframe = 0;
    do
    {
        sprintf(buf, "%s frame %d", fn, nframe);
        indent = 0;
        indent = pr_title(stdout, indent, buf);
        pr_indent(stdout, indent);
        fprintf(stdout, "natoms=%10d  step=%10" PRId64 "  time=%12.7e  prec=%10g\n", natoms, step, time, prec);
        pr_rvecs(stdout, indent, "box", box, DIM);
        pr_rvecs(stdout, indent, "x", x, natoms);
        nframe++;
    } while (read_next_xtc(xd, natoms, &step, &time, box, x, &prec, &bOK) != 0);
    if (!bOK)
    {
        fprintf(stderr, "\nWARNING: Incomplete frame at time %g\n", time);
    }
    sfree(x);
    close_xtc(xd);
}

#if GMX_USE_TNG

/*! \brief Callback used by list_tng_for_gmx_dump. */
void list_tng_inner(const char* fn,
                    gmx_bool    bFirstFrame,
                    real*       values,
                    int64_t     step,
                    double      frame_time,
                    int64_t     n_values_per_frame,
                    int64_t     n_atoms,
                    real        prec,
                    int64_t     nframe,
                    char*       block_name)
{
    char buf[256];
    int  indent = 0;

    if (bFirstFrame)
    {
        sprintf(buf, "%s frame %" PRId64, fn, nframe);
        indent = 0;
        indent = pr_title(stdout, indent, buf);
        pr_indent(stdout, indent);
        fprintf(stdout, "natoms=%10" PRId64 "  step=%10" PRId64 "  time=%12.7e", n_atoms, step, frame_time);
        if (prec > 0)
        {
            fprintf(stdout, "  prec=%10g", prec);
        }
        fprintf(stdout, "\n");
    }
    pr_reals_of_dim(stdout, indent, block_name, values, n_atoms, n_values_per_frame);
}

#endif

//! Dump a TNG file
void list_tng(const char* fn)
{
#if GMX_USE_TNG
    gmx_tng_trajectory_t tng;
    int64_t              nframe = 0;
    int64_t              i, *block_ids = nullptr, step, ndatablocks;
    gmx_bool             bOK;
    real*                values = nullptr;

    gmx_tng_open(fn, 'r', &tng);
    gmx_print_tng_molecule_system(tng, stdout);

    bOK = gmx_get_tng_data_block_types_of_next_frame(tng, -1, 0, nullptr, &step, &ndatablocks, &block_ids);
    do
    {
        for (i = 0; i < ndatablocks; i++)
        {
            double  frame_time;
            real    prec;
            int64_t n_values_per_frame, n_atoms;
            char    block_name[STRLEN];

            gmx_get_tng_data_next_frame_of_block_type(tng,
                                                      block_ids[i],
                                                      &values,
                                                      &step,
                                                      &frame_time,
                                                      &n_values_per_frame,
                                                      &n_atoms,
                                                      &prec,
                                                      block_name,
                                                      STRLEN,
                                                      &bOK);
            if (!bOK)
            {
                /* Can't write any output because we don't know what
                   arrays are valid. */
                fprintf(stderr, "\nWARNING: Incomplete frame at time %g, will not write output\n", frame_time);
            }
            else
            {
                list_tng_inner(
                        fn, (0 == i), values, step, frame_time, n_values_per_frame, n_atoms, prec, nframe, block_name);
            }
        }
        nframe++;
    } while (gmx_get_tng_data_block_types_of_next_frame(
            tng, step, 0, nullptr, &step, &ndatablocks, &block_ids));

    if (block_ids)
    {
        sfree(block_ids);
    }
    sfree(values);
    gmx_tng_close(&tng);
#else
    GMX_UNUSED_VALUE(fn);
#endif
}

//! Dump a trajectory file
void list_trx(const char* fn)
{
    switch (fn2ftp(fn))
    {
        case efXTC: list_xtc(fn); break;
        case efTRR: list_trr(fn); break;
        case efTNG: list_tng(fn); break;
        default:
            fprintf(stderr, "File %s is of an unsupported type. Try using the command\n 'less %s'\n", fn, fn);
    }
}

//! Dump an energy file
void list_ene(const char* fn)
{
    ener_file_t  in;
    gmx_bool     bCont;
    gmx_enxnm_t* enm = nullptr;
    t_enxframe*  fr;
    int          i, j, nre, b;
    char         buf[22];

    printf("gmx dump: %s\n", fn);
    in = open_enx(fn, "r");
    do_enxnms(in, &nre, &enm);
    assert(enm);

    printf("energy components:\n");
    for (i = 0; (i < nre); i++)
    {
        printf("%5d  %-24s (%s)\n", i, enm[i].name, enm[i].unit);
    }

    snew(fr, 1);
    do
    {
        bCont = do_enx(in, fr);

        if (bCont)
        {
            printf("\n%24s  %12.5e  %12s  %12s\n", "time:", fr->t, "step:", gmx_step_str(fr->step, buf));
            printf("%24s  %12s  %12s  %12s\n", "", "", "nsteps:", gmx_step_str(fr->nsteps, buf));
            printf("%24s  %12.5e  %12s  %12s\n", "delta_t:", fr->dt, "sum steps:", gmx_step_str(fr->nsum, buf));
            if (fr->nre == nre)
            {
                printf("%24s  %12s  %12s  %12s\n",
                       "Component",
                       "Energy",
                       "Av. Energy",
                       "Sum Energy");
                if (fr->nsum > 0)
                {
                    for (i = 0; (i < nre); i++)
                    {
                        printf("%24s  %12.5e  %12.5e  %12.5e\n",
                               enm[i].name,
                               fr->ener[i].e,
                               fr->ener[i].eav,
                               fr->ener[i].esum);
                    }
                }
                else
                {
                    for (i = 0; (i < nre); i++)
                    {
                        printf("%24s  %12.5e\n", enm[i].name, fr->ener[i].e);
                    }
                }
            }
            for (b = 0; b < fr->nblock; b++)
            {
                const char* typestr = "";

                t_enxblock* eb = &(fr->block[b]);
                printf("Block data %2d (%3d subblocks, id=%d)\n", b, eb->nsub, eb->id);

                if (eb->id < enxNR)
                {
                    typestr = enx_block_id_name[eb->id];
                }
                printf("  id='%s'\n", typestr);
                for (i = 0; i < eb->nsub; i++)
                {
                    t_enxsubblock* sb = &(eb->sub[i]);
                    printf("  Sub block %3d (%5d elems, type=%s) values:\n",
                           i,
                           sb->nr,
                           enumValueToString(sb->type));

                    switch (sb->type)
                    {
                        case XdrDataType::Float:
                            for (j = 0; j < sb->nr; j++)
                            {
                                printf("%14d   %8.4f\n", j, sb->fval[j]);
                            }
                            break;
                        case XdrDataType::Double:
                            for (j = 0; j < sb->nr; j++)
                            {
                                printf("%14d   %10.6f\n", j, sb->dval[j]);
                            }
                            break;
                        case XdrDataType::Int:
                            for (j = 0; j < sb->nr; j++)
                            {
                                printf("%14d %10d\n", j, sb->ival[j]);
                            }
                            break;
                        case XdrDataType::Int64:
                            for (j = 0; j < sb->nr; j++)
                            {
                                printf("%14d %s\n", j, gmx_step_str(sb->lval[j], buf));
                            }
                            break;
                        case XdrDataType::Char:
                            for (j = 0; j < sb->nr; j++)
                            {
                                printf("%14d %1c\n", j, sb->cval[j]);
                            }
                            break;
                        case XdrDataType::String:
                            for (j = 0; j < sb->nr; j++)
                            {
                                printf("%14d %80s\n", j, sb->sval[j]);
                            }
                            break;
                        default: gmx_incons("Unknown subblock type");
                    }
                }
            }
        }
    } while (bCont);

    close_enx(in);

    free_enxframe(fr);
    sfree(fr);
    sfree(enm);
}

//! Dump a (Hessian) matrix file
void list_mtx(const char* fn)
{
    int                 nrow, ncol, i, j, k;
    real *              full   = nullptr, value;
    gmx_sparsematrix_t* sparse = nullptr;

    gmx_mtxio_read(fn, &nrow, &ncol, &full, &sparse);

    if (full == nullptr)
    {
        snew(full, nrow * ncol);
        for (i = 0; i < nrow * ncol; i++)
        {
            full[i] = 0;
        }

        for (i = 0; i < sparse->nrow; i++)
        {
            for (j = 0; j < sparse->ndata[i]; j++)
            {
                k                  = sparse->data[i][j].col;
                value              = sparse->data[i][j].value;
                full[i * ncol + k] = value;
                full[k * ncol + i] = value;
            }
        }
        gmx_sparsematrix_destroy(sparse);
    }

    printf("%d %d\n", nrow, ncol);
    for (i = 0; i < nrow; i++)
    {
        for (j = 0; j < ncol; j++)
        {
            printf(" %g", full[i * ncol + j]);
        }
        printf("\n");
    }

    sfree(full);
}

class Dump : public ICommandLineOptionsModule
{
public:
    Dump() {}

    // From ICommandLineOptionsModule
    void init(CommandLineModuleSettings* /*settings*/) override {}

    void initOptions(IOptionsContainer* options, ICommandLineOptionsModuleSettings* settings) override;

    void optionsFinished() override;

    int run() override;

private:
    //! Commandline options
    //! \{
    bool bShowNumbers_      = true;
    bool bShowParams_       = false;
    bool bSysTop_           = false;
    bool bOriginalInputrec_ = false;
    //! \}
    //! Commandline file options
    //! \{
    std::string inputTprFilename_;
    std::string inputTrajectoryFilename_;
    std::string inputEnergyFilename_;
    std::string inputCheckpointFilename_;
    std::string inputTopologyFilename_;
    std::string inputMatrixFilename_;
    std::string outputMdpFilename_;
    //! \}
};

void Dump::initOptions(IOptionsContainer* options, ICommandLineOptionsModuleSettings* settings)
{
    const char* desc[] = { "[THISMODULE] reads a run input file ([REF].tpr[ref]),",
                           "a trajectory ([REF].trr[ref]/[REF].xtc[ref]/[TT]tng[tt]), an energy",
                           "file ([REF].edr[ref]), a checkpoint file ([REF].cpt[ref])",
                           "or topology file ([REF].top[ref])",
                           "and prints that to standard output in a readable format.",
                           "This program is essential for checking your run input file in case of",
                           "problems." };
    settings->setHelpText(desc);

    const char* bugs[] = {
        "The [REF].mdp[ref] file produced by [TT]-om[tt] can not be read by grompp."
    };
    settings->setBugText(bugs);
    // TODO If this ancient note acknowledging a bug is still true,
    // fix it or block that run path:
    //   Position restraint output from -sys -s is broken

    options->addOption(FileNameOption("s")
                               .filetype(OptionFileType::RunInput)
                               .inputFile()
                               .store(&inputTprFilename_)
                               .description("Run input file to dump"));
    options->addOption(FileNameOption("f")
                               .filetype(OptionFileType::Trajectory)
                               .inputFile()
                               .store(&inputTrajectoryFilename_)
                               .description("Trajectory file to dump"));
    options->addOption(FileNameOption("e")
                               .filetype(OptionFileType::Energy)
                               .inputFile()
                               .store(&inputEnergyFilename_)
                               .description("Energy file to dump"));
    options->addOption(
            FileNameOption("cp").legacyType(efCPT).inputFile().store(&inputCheckpointFilename_).description("Checkpoint file to dump"));
    options->addOption(
            FileNameOption("p").legacyType(efTOP).inputFile().store(&inputTopologyFilename_).description("Topology file to dump"));
    options->addOption(
            FileNameOption("mtx").legacyType(efMTX).inputFile().store(&inputMatrixFilename_).description("Hessian matrix to dump"));
    options->addOption(FileNameOption("om")
                               .legacyType(efMDP)
                               .outputFile()
                               .store(&outputMdpFilename_)
                               .description("grompp input file from run input file"));

    options->addOption(
            BooleanOption("nr").store(&bShowNumbers_).defaultValue(true).description("Show index numbers in output (leaving them out makes comparison easier, but creates a useless topology)"));
    options->addOption(
            BooleanOption("param").store(&bShowParams_).defaultValue(false).description("Show parameters for each bonded interaction (for comparing dumps, it is useful to combine this with -nonr)"));
    options->addOption(BooleanOption("sys").store(&bSysTop_).defaultValue(false).description(
            "List the atoms and bonded interactions for the whole system instead of for each "
            "molecule type"));
    options->addOption(
            BooleanOption("orgir").store(&bOriginalInputrec_).defaultValue(false).description("Show input parameters from tpr as they were written by the version that produced the file, instead of how the current version reads them"));
}

void Dump::optionsFinished()
{
    // TODO Currently gmx dump ignores user input that seeks to dump
    // multiple files. Here, we could enforce that the user only asks
    // to dump one file.
}

int Dump::run()
{
    if (!inputTprFilename_.empty())
    {
        list_tpr(inputTprFilename_.c_str(),
                 bShowNumbers_,
                 bShowParams_,
                 outputMdpFilename_.empty() ? nullptr : outputMdpFilename_.c_str(),
                 bSysTop_,
                 bOriginalInputrec_);
    }
    else if (!inputTrajectoryFilename_.empty())
    {
        list_trx(inputTrajectoryFilename_.c_str());
    }
    else if (!inputEnergyFilename_.empty())
    {
        list_ene(inputEnergyFilename_.c_str());
    }
    else if (!inputCheckpointFilename_.empty())
    {
        list_checkpoint(inputCheckpointFilename_.c_str(), stdout);
    }
    else if (!inputTopologyFilename_.empty())
    {
        list_top(inputTopologyFilename_.c_str());
    }
    else if (!inputMatrixFilename_.empty())
    {
        list_mtx(inputMatrixFilename_.c_str());
    }

    return 0;
}

} // namespace

LIBGROMACS_EXPORT const char     DumpInfo::name[]             = "dump";
LIBGROMACS_EXPORT const char     DumpInfo::shortDescription[] = "Make binary files human readable";
ICommandLineOptionsModulePointer DumpInfo::create()
{
    return std::make_unique<Dump>();
}

} // namespace gmx
