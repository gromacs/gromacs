/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands.
 * Copyright (c) 2018, by the GROMACS development team, led by
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

#include "g_fca.h"

#include "config.h"

#include <ctype.h>
#include <math.h>
#include <signal.h>
#include <string.h>

#include <cassert>

#include <utility>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/index.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/smalloc.h"

#include "eigenvec.hpp"
#include "entropy_basic_histo.hpp"
#include "fca_minimizing.h"
#include "frame_set.hpp"
#include "input_functions.hpp"
#include "utils.hpp"


namespace FCA
{
static void write_fca_valfile_fca(FILE* out, const int dim, const real* mat, const int natoms,
                                  const real* amplitude, const real* anharmonicity)
{

    std::unique_ptr<real[]> collectivity(new real[dim]);

    if (mat)
    {
        /** computes the collectivity of 1..dim eigvecs ( mat[k*ndim+j*DIM+d] ) k-mode j-atom d-[xyz]
           this computes the squared sum of triples and then the collectivity
           ASSUMPTION: one atom is three consecutive entries in vecs
           RETURN: collectivity[0..dim-1] **/
        const int  ndim    = DIM * natoms;
        const real invlogn = 1.0 / log(natoms);
        for (int k = 0; k < dim; k++)
        {
            collectivity[k] = 0;
            for (int j = 0; j < natoms; j++)
            {
                real sqratom = 0;
                for (int d = 0; d < DIM; d++)
                {
                    sqratom += utils::squareof(mat[k * ndim + j * DIM + d]);
                }
                collectivity[k] = -invlogn * sqratom * log(sqratom);
            }
        }
    }
    for (int k = 0; k < dim; k++)
    {
        if (mat)
        {
            fprintf(out, "%10.5f %10.5f %10.5f\n", amplitude[k],
                    anharmonicity[k], collectivity[k]);
        }
        else
        {
            fprintf(out, "%10.5f %10.5f\n", amplitude[k], anharmonicity[k]);
        }
    }
}

static std::vector<int> select_eigen_vectors(const int nbAtoms, const int nvec, const int eignr[], const int first, int last)
{
    if (last == -1)
    {
        last = nbAtoms * DIM;
    }

    std::vector<int> iout;

    if (first > -1)
    {
        /* make an index from first to last */
        iout.resize(last - first + 1);
        for (int i = 0; i < int(iout.size()); i++)
        {
            iout[i] = first - 1 + i;
        }
    }
    else
    {
        printf("Select eigenvectors for output, end your selection with 0\n");
        while (true)
        {
            int inputnb;
            scanf("%d", &inputnb);
            if (inputnb <= 0)
            {
                break;
            }
            iout.emplace_back(inputnb);
        }
        printf("\n");
    }
    /* make an index of the eigenvectors which are present */
    std::vector<int> outvec;

    fprintf(stderr,
            "%d eigenvectors in file -- looking for  %d eigenvectors \n", nvec,
            int(iout.size()));
    for (int i = 0; i < int(iout.size()); i++)
    {
        int j = 0;
        while ((j < nvec) && eignr && (eignr[j] != iout[i]))
        {
            j++;
        }
        if ((j < nvec) && eignr && (eignr[j] == iout[i]))
        {
            outvec.emplace_back(j);
        }
    }
    fprintf(stderr, "%zu eigenvectors selected for output", outvec.size());
    if (outvec.size() <= 100 && eignr)
    {
        fprintf(stderr, ":");
        for (int j = 0; j < int(outvec.size()); j++)
        {
            fprintf(stderr, " %d", eignr[outvec[j]] + 1);
        }
    }
    fprintf(stderr, "\n");
    return outvec;
}

}

int gmx_fca(int argc, char* argv[])
{
    const char* desc[] =
    {
        "TO DO: flag for fit in fca_vec.trr is not set if fit-group is different from PCA group",
        "g_fca writes out collective coordinates describing the ensemble motion",
        "in fca_vec.trr (c.f., g_covar,eigvec.trr)."
        "[PAR]FCA (Full Correlation Analysis) minimizes the mutual information between",
        "the collective coordinates, and thus generalizes PCA (g_covar), which minimizes just the covariations",
        "between the coordinates. The improvement over PCA is that a huge part",
        "of the correlation between coordinates remains undetected by PCA, such",
        "that the principal components often do not resolve collective motions as",
        "optimally as possible with FCA. See also",
        "[PAR]Full Correlation Analysis, Oliver F. Lange and Helmut Grubmueller, submitted to Proteins[PAR]",
        "[PAR]Note that, FCA coordinates are gained by an orthogonal transformation of the atomic coordintates",
        "in configurational space. Thus, they do not change the length-scale or other geometric properties",
        "[PAR]g_fca can operate on two differnt kinds of inputs",
        "the most often used is probably that you define an input trajectory with -f and the eigvec.trr gained from",
        "g_covar. The PCA coordinates of eigvec.trr are used as a first guess.",
        "with -last  you select the 1..last  principal components on which you will perform minimization",
        "the remaining principal components will be copied untouched to the output <fca_vec.trr>",
        "From our experience the small-amplitude (e.g., all but the first 100) PCA modes are sufficiently uncorrelated",
        "such that the FCA analysis does not improve considerably (but slows down), if small amplitude modes are included",
        "[PAR]The second input alternative uses the option -pos to read any kind of high dimensional input data",
        "g_fca will interpret the rows as high-dimensional samples, i.e., a column contains samples of a single coordinate",
        "[PAR]The ic_matrix.xvg contains the matrix for the coordinate transformation from the input data (i.e., PCA-modes or -pos)",
        "to the optimal set of FCA modes.",
        "The file fca_vec.trr is only written if g_fca was started with -v. It contains the coordinate transformation from cartesian",
        "atomic coordinates to the FCA modes and can be used, for example, with g_anaeig to analyze the trajectory",
        "The file fca_val.xvg contains characteristics of the FCA modes. The first column contains the variance of the mode <x^2>, this",
        "is the same information written in eigval.xvg of g_covar",
        "The second column of fca_val.xvg contains the neg-entropy (anharmonicity) of a FCA-mode (see paper). A high value",
        "indicates that this mode resolves a transition between conformational states",
        "the third column of fca_val.xvg contains the collectivity of a FCA-mode (only with -v option).",
        "This number quantifies how many atoms are contributing to the motion described by the respective FCA-mode",
        "Note, that g_fca sorts the modes by neg-entropy and variance. You can set the weights for the sorting criterium",
        "with the options -negw and -ampw for neg-entropy and variance, respectively",
        "[PAR]Note, that you can use g_fca with the option -nominimize to resort the modes given in the input file.",
        "With Ctrl-C you can stop the minimization and with -r ic_matrix.xvg you can restart the minimization with the last",
        "transformation matrix taken from the -ic output.",
        "[PAR]The restart file is not written in full precision (trr) this can lead to differences",
        "in the MI-sum after a restart. TODO: write restart file as .trr trajectory, if multiple frames",
        "ask which to be used for restart",

    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    FCA::utils::mpi mpi(MPI_COMM_WORLD);
    FCA::utils::Log logger;

    if (mpi.isMaster())
    {
        logger.MI = gmx_ffopen("MI_dump.log", "w");
        setbuf(logger.MI, nullptr);
        logger.ic_matrix = gmx_ffopen("ic_matrix.log", "w");
        setbuf(logger.ic_matrix, nullptr);
        //        logger.move_log    = fopen("moves.log","w");
        //        logger.move_matrix = fopen("move_matrix.log","w");
        //        logger.MI_matrix   = fopen("MI_matrix.log","w");
        //        logger.icx         = fopen("ic_lastproj.log","w");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    t_filenm fnm[] =
    {
        { efTRN, "-v", "eigenvec", ffOPTRD, 0, nullptr },
        { efTRX, "-f", nullptr, ffOPTRD, 0, nullptr },
        { efTPS, nullptr, nullptr, ffOPTRD, 0, nullptr },
        { efNDX, nullptr, nullptr, ffOPTRD, 0, nullptr },
        { efXVG, "-pos", "pos", ffOPTRD, 0, nullptr },
        { efXVG, "-r", "restart_ic_matrix", ffOPTRD, 0, nullptr },
        { efXVG, "-ic", "ic_matrix", ffWRITE, 0, nullptr },
        { efTRN, "-o", "fca_vec", ffWRITE, 0, nullptr },
        { efLOG, "-l", "fca.log", ffWRITE, 0, nullptr },
        { efXVG, "-val", "fca_val", ffWRITE, 0, nullptr },
        { efXVG, "-mi", "mi_matrix", ffOPTWR, 0, nullptr },
    };

    std::unique_ptr< FCA::EigenVec >         eigVecManager;
    std::vector< std::unique_ptr< real[] > > projx;
    FCA::utils::TrxInfo trxinfo;

    int                 number_of_moves  = 0;
    int                 dim              = -1;
    int                 nframes          = 0;
    double              log_kappa1       = -2.0;
    double              log_kappa2       = -0.2;
    gmx_bool            bMinimize        = TRUE;
    gmx_bool            bSort            = TRUE;
    real                amplitude_weight = 1;
    gmx_bool            bBasicHisto      = FALSE;
    const char        * PosFile          = nullptr;

    if (mpi.isMaster())
    {
        int     last            = 8;
        int     load_frame_step = 1;

        t_pargs pa[] =
        { { "-last", FALSE, etINT, { &last }, "Last eigenvector for analysis (-1 is till the last)" },
          { "-load_frame_step", FALSE, etINT, { &load_frame_step }, "Only analyse every nr-th frame" },
          { "-minimize", FALSE, etBOOL, { &bMinimize }, "minimize mutual information (do FCA analysis)" },
          { "-ampw", FALSE, etREAL, { &amplitude_weight }, "weighting of mode-amplitude for sorting of FCA-modes" },
          { "-sort", FALSE, etBOOL, { &bSort }, "sort the modes by mode-amplitude/neg-entropy (c.f., -negw -ampw)" },
          { "-basic", FALSE, etBOOL, { &bBasicHisto }, "use a basic histogram approach" },
          { "-kappa1", FALSE, etREAL, { &log_kappa1 }, "log of kappa1 for basic histogram" },
          { "-kappa2", FALSE, etREAL, { &log_kappa2 }, "log of kappa2 for basic histogram" } };

        gmx_output_env_t* oenv;
        output_env_init_default(&oenv);

        parse_common_args(&argc, argv,
                          PCA_CAN_TIME | PCA_TIME_UNIT | PCA_CAN_VIEW /* | PCA_BE_NICE*/,
                          asize(fnm), fnm, asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv);

        PosFile          = opt2fn_null("-pos", asize(fnm), fnm);
        trxinfo.setTrj(opt2fn_null("-f", asize(fnm), fnm));
        trxinfo.setPca(opt2fn_null("-v", asize(fnm), fnm));

        // -----------------------------  do the projection -----------------------------------------------

        fprintf(stderr, "start projection \n");
        if (opt2bSet("-f", asize(fnm), fnm) && opt2bSet("-v", asize(fnm), fnm))
        {
            eigVecManager.reset(new FCA::EigenVec(opt2fn_null("-v", asize(fnm), fnm),
                                                  ftp2fn(efTPS, asize(fnm), fnm),
                                                  ftp2fn_null(efNDX, asize(fnm), fnm),
                                                  nullptr));

            std::vector<int> proj_info = FCA::select_eigen_vectors(eigVecManager->getNbAtoms(), eigVecManager->getNbVec(),
                                                                   eigVecManager->geteignr(), 1, last);


            dim            = int(proj_info.size());
            trxinfo.setDimPca(eigVecManager->getNbVec());

            FCA::FrameSet frameSet(opt2fn("-f", asize(fnm), fnm), oenv, load_frame_step);

            auto          frames = frameSet.releaseFrames();
            FCA::utils::apply_rm_pbc_plain_to_all(frameSet.getNbAtoms(), &frames, frameSet.getBox());

            projx = eigVecManager->project_vectors(frames, nullptr, nullptr, nullptr, nullptr, frameSet.getTimeFrame(),
                                                   proj_info.size(), proj_info.data());
        }
        else if (PosFile)
        {
            projx   = FCA::input_functions::read_fca_proj(PosFile, &dim);
            nframes = projx.size();
        }
        else
        {
            fprintf(stderr, "you need to provide either a trajectory and eigenvec.trr or with -pos a file of any set of vectors\n");
            return 1;
        }
        fprintf(stderr, "start processing...\n");
        fprintf(stderr, "init FCA...\n");

        if (mpi.getNumNodes() > 1)
        {
            number_of_moves = lrint(ceil(10.0 / (mpi.getNumNodes() - 1))) * (mpi.getNumNodes() - 1);
        }
        else
        {
            number_of_moves = 5;
        }
    } //if MASTER

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

#if GMX_MPI == 1
    FCA::utils::AssertMpi(MPI_Bcast(&dim, 1, MPI_INT, FCA_MASTER_RANK, MPI_COMM_WORLD), __FILE__, __LINE__);
    FCA::utils::AssertMpi(MPI_Bcast(&nframes, 1, MPI_INT, FCA_MASTER_RANK, MPI_COMM_WORLD), __FILE__, __LINE__);
#endif
    FCA::FcaMaster fca(dim, nframes, std::move(projx), mpi);

    if (mpi.isMaster() && opt2bSet("-r", asize(fnm), fnm))
    {
        fprintf(stderr, "accessing restart file...\n");
        fca.read_ic_matrix(opt2fn("-r", asize(fnm), fnm));
        fca.dump_ic_matrix(logger.ic_matrix);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
#if GMX_MPI == 1
    FCA::utils::AssertMpi(MPI_Bcast(&number_of_moves, 1, MPI_INT, FCA_MASTER_RANK, MPI_COMM_WORLD), __FILE__, __LINE__);
    FCA::utils::AssertMpi(MPI_Bcast(&bMinimize, 1, MPI_INT, FCA_MASTER_RANK, MPI_COMM_WORLD), __FILE__, __LINE__);
#endif
    /// ---- Minimize
    if (bMinimize)
    {
        fprintf(stderr, "minimize correlation for %d frames, %d trial moves on %d slave-nodes per step\n", nframes, number_of_moves, mpi.getNumNodes() - 1);

        if (mpi.isMaster())
        {
            fprintf(stderr, "distribute data on nodes...\n");
        }

        fca.broadcast_ic_matrix(); // to sync restart ic_matrix

        if (mpi.isMaster())
        {
            fprintf(stderr, "start minimizing...\n");
        }

        static bool stop = false;
        signal(SIGTERM, [](int signalp){
                   if (signalp == SIGTERM)
                   {
                       fprintf(stderr, "Received sigterm...\n");
                       stop = true;
                   }
               });

        // Loop
        while (fca.minimization(number_of_moves, logger) && !stop)
        {
            if (mpi.isMaster())
            {
                fprintf(stderr, "MI-sum: %10.5f\r", fca.getMIsum());
            }
        } // while (minimize)

        signal(SIGTERM, SIG_IGN);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
#if GMX_MPI == 1
    FCA::utils::AssertMpi(MPI_Bcast(&bSort, 1, MPI_INT, FCA_MASTER_RANK, MPI_COMM_WORLD), __FILE__, __LINE__);
    FCA::utils::AssertMpi(MPI_Bcast(&amplitude_weight, 1, MPI_INT, FCA_MASTER_RANK, MPI_COMM_WORLD), __FILE__, __LINE__);
#endif
    /// --- POST PROCESSING ----- single processor
    if (mpi.isMaster())
    {
        std::unique_ptr<real[]> amplitude(new real[dim]());
        std::unique_ptr<real[]> anharmonicity(new real[dim]());

        fca.project_ic_signal(); //* to be sure to have the actual data
        fca.compute_S1();
        fca.compute_mode_amplitude(amplitude.get());
        fca.compute_mode_anharmonicity(amplitude.get(), anharmonicity.get());
        if (bSort)
        {
            fca.sort_modes(amplitude.get(), anharmonicity.get(), amplitude_weight,
                           1.0 - amplitude_weight);
        }

        // output ic-matrix
        {
            FILE* icout = gmx_ffopen(opt2fn("-ic", asize(fnm), fnm), "w");
            fca.dump_ic_matrix(icout);
            fclose(icout);
        }
        //* output log-file
        trxinfo.write_fcalog(opt2fn("-l", asize(fnm), fnm), fca, PosFile,
                             opt2fn_null("-v", asize(fnm), fnm));

        // output icx for testing
        int z, y;
        if (logger.icx)
        {
            for (y = 0; y < fca.getNFrames(); y++)
            {
                for (z = 0; z < fca.getDim(); z++)
                {
                    fprintf(logger.icx, "%10.5f", fca.getIcx()[z][y]);
                }
                fprintf(logger.icx, "\n");
            }
        }
        std::unique_ptr< real[] > fca_eigvec;

        // output fca_vec.trr
        if (opt2bSet("-v", asize(fnm), fnm))
        {
            assert(eigVecManager != nullptr);

            fca_eigvec = eigVecManager->produce_new_eigvecs(fca.getIcMatrix(), fca.getDim());
            eigVecManager->setVal(amplitude.get(), dim);
            eigVecManager->writeEigvecsWithMat(fca_eigvec.get(), opt2fn("-o", asize(fnm), fnm));

            // if (opt2bSet("-ascii",asize(fnm),fnm))
            // write_eigvecs_ascii_fca(mat,eigvec->getNbAtoms(), eigVecManager->getNVec(),opt2fn("-ascii",asize(fnm),fnm));
        } // if opt2bSet("-v")

        // output fca_val.xvg
        {
            FILE* val_out = gmx_ffopen(opt2fn("-val", asize(fnm), fnm), "w");
            FCA::write_fca_valfile_fca(val_out, fca.getDim(), fca_eigvec.get(),
                                       eigVecManager ? eigVecManager->getNbAtoms() : 0,
                                       amplitude.get(), anharmonicity.get());
            fclose(val_out);
        }
    } // endif MASTER

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    /// ----  More POST-PROCESSING --- in parallel this time...

    gmx_bool bWriteMI = FALSE;

    if (mpi.isMaster())
    {
        bWriteMI = opt2bSet("-mi", asize(fnm), fnm);
    }
#if GMX_MPI == 1
    // distribute flag
    FCA::utils::AssertMpi(MPI_Bcast(&bWriteMI, 1, MPI_INT, FCA_MASTER_RANK, MPI_COMM_WORLD), __FILE__, __LINE__);
    FCA::utils::AssertMpi(MPI_Bcast(&log_kappa1, 1, MPI_DOUBLE, FCA_MASTER_RANK, MPI_COMM_WORLD), __FILE__, __LINE__);
    FCA::utils::AssertMpi(MPI_Bcast(&log_kappa2, 1, MPI_DOUBLE, FCA_MASTER_RANK, MPI_COMM_WORLD), __FILE__, __LINE__);
    FCA::utils::AssertMpi(MPI_Bcast(&bBasicHisto, 1, MPI_INT, FCA_MASTER_RANK, MPI_COMM_WORLD), __FILE__, __LINE__);
#endif
    // compute MI matrix
    if (bWriteMI)
    {
        if (mpi.isMaster())
        {
            fprintf(stderr, "going to compute MI Matrix for output ....\n");
        }

        fca.broadcast_ic_matrix();

        if (bBasicHisto)
        {
            FCA::EntropyHisto::Basic_compute_MI_matrix(&fca, log_kappa1, log_kappa2);
        }
        else
        {
            gmx::ThreeFry2x64<64> rng(123456, gmx::RandomDomain::Other);
            fca.compute_MI_matrix(nullptr, &rng);
        }

        if (mpi.isMaster())
        {
            FILE* out = gmx_ffopen(opt2fn("-mi", asize(fnm), fnm), "w");
            fca.write_MI_matrix(out);
            gmx_ffclose(out);
        }
    }

    return 0;
}
