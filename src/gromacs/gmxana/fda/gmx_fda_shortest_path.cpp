/*
 * gmx_fda_shortest_path.cpp
 *
 *  Created on: Feb 13, 2015
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/gmxana/fda/BoostGraph.h"
#include "gromacs/gmxana/fda/EnumParser.h"
#include "gromacs/gmxana/fda/FrameType.h"
#include "gromacs/gmxana/fda/Helpers.h"
#include "gromacs/gmxana/fda/PDB.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

using namespace fda_analysis;

#define PRINT_DEBUG

int gmx_fda_shortest_path(int argc, char *argv[])
{
    const char *desc[] = {
        "[THISMODULE] calculates the k-shortest paths between a source node "
    	"and a destination node of a FDA force network. "
        "If the optional file [TT]-ipf-diff[tt] is used the differences of the pairwise forces will be taken. "
        "The option [TT]-pymol[tt] can be used to generate a Pymol script, which can be directly called by Pymol. "
        "The graph will printed in the PDB-format, which allow an easy visualization with "
    	"an program of your choice. "
    	"Each path will be determined and segment names will be assign to each "
        "of them, thus coloring them by segment id will help the analysis "
        "(32 different colors). The Bfactor column will be used for the value of "
        "the force and helps the coloring as a function of the force magnitude. "
        "The CONNECT header will be used to create bonds between nodes. "
    };

    gmx_output_env_t *oenv;
    static const char* frameString = "average 1";
    static int source = 0;
    static int dest = 0;
    static int numberOfShortestPaths = 1;
    static bool convert = false;

    t_pargs pa[] = {
        { "-frame", FALSE, etSTR, {&frameString}, "Specify a single frame number or \"average n\" to take the mean over every n-th frame"
                " or \"skip n\" to take every n-th frame or \"all\" to take all frames (e.g. for movies)" },
        { "-source", FALSE, etINT, {&source}, "Source node of the path" },
        { "-dest", FALSE, etINT, {&dest}, "Destination point of the path" },
        { "-nk", FALSE, etINT, {&numberOfShortestPaths}, "Number of shortest paths" },
        { "-convert", FALSE, etBOOL, {&convert}, "Convert force unit from kJ/mol/nm into pN" }
    };

    t_filenm fnm[] = {
        { efPFX, "-ipf", NULL, ffREAD },
        { efPFX, "-ipf-diff", NULL, ffOPTRD },
        { efTPS, NULL, NULL, ffREAD },
        { efTRX, "-traj", NULL, ffOPTRD },
        { efNDX, NULL, NULL, ffOPTRD },
        { efPDB, "-o", "result", ffWRITE },
        { efPML, "-pymol", "result", ffOPTWR }
    };

#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, PCA_CAN_TIME,
        NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv)) return 0;

    if (opt2bSet("-ipf-diff", NFILE, fnm) and (fn2ftp(opt2fn("-ipf-diff", NFILE, fnm)) != fn2ftp(opt2fn("-ipf", NFILE, fnm))))
        gmx_fatal(FARGS, "Type of the file (-ipf-diff) does not match the type of the file (-ipf).");

    if (fn2ftp(opt2fn("-ipf", NFILE, fnm)) == efPFR and !opt2bSet("-n", NFILE, fnm))
        gmx_fatal(FARGS, "Index file is needed for residuebased pairwise forces.");

    // Get number of particles
    int nbParticles = getMaxIndexSecondColumnFirstFrame(opt2fn("-ipf", NFILE, fnm)) + 1;
    int nbParticles2 = nbParticles * nbParticles;

    // Interactive input of group name for residue model points
    int isize = 0;
    int *index;
    char *grpname;
    fprintf(stderr, "\nSelect group for residue model points:\n");
    rd_index(ftp2fn(efNDX, NFILE, fnm), 1, &isize, &index, &grpname);
    if (isize != nbParticles) gmx_fatal(FARGS, "Number of atoms in group does not match number of FDA points.");

    int frameValue;
    FrameType frameType = getFrameTypeAndSkipValue(frameString, frameValue);

	#ifdef PRINT_DEBUG
		std::cerr << "frameType = " << EnumParser<FrameType>()(frameType) << std::endl;
		std::cerr << "frameValue = " << frameValue << std::endl;
		std::cerr << "Number of particles (np) = " << nbParticles << std::endl;
		std::cerr << "source = " << source << std::endl;
		std::cerr << "dest = " << dest << std::endl;
		std::cerr << "Number of shortest paths (nk) = " << numberOfShortestPaths << std::endl;
		std::cerr << "convert = " << convert << std::endl;
		std::cerr << "pfx filename = " << opt2fn("-ipf", NFILE, fnm) << std::endl;
		if (opt2bSet("-ipf-diff", NFILE, fnm)) std::cerr << "pfx-diff filename = " << opt2fn("-ipf-diff", NFILE, fnm) << std::endl;
	    std::cerr << "structure filename = " << opt2fn("-s", NFILE, fnm) << std::endl;
		std::cerr << "result filename = " << opt2fn("-o", NFILE, fnm) << std::endl;
		std::cerr << "trajectory filename = " << opt2fn("-traj", NFILE, fnm) << std::endl;
    #endif

    // Read input structure coordinates
    rvec *coord;
    t_topology top;
    int ePBC;
    matrix box;
    read_tps_conf(ftp2fn(efTPS, NFILE, fnm), &top, &ePBC, &coord, NULL, box, TRUE);

    PDB pdb(opt2fn("-s", NFILE, fnm), std::vector<int>(index, index + isize));

    std::vector<double> forceMatrix, forceMatrix2;

    // Pymol pml-file
    std::string molecularTrajectoryFilename = "traj.pdb";
    FILE *molecularTrajectoryFile = NULL;
    if (opt2bSet("-pymol", NFILE, fnm)) {
        std::ofstream pmlFile(opt2fn("-pymol", NFILE, fnm));
        pmlFile << "load " << molecularTrajectoryFilename << ", object=trajectory" << std::endl;
        pmlFile << "set connect_mode, 1" << std::endl;
        pmlFile << "load " << opt2fn("-o", NFILE, fnm) << ", object=network, discrete=1, multiplex=1" << std::endl;
        pmlFile << "spectrum segi, blue_white_red, network" << std::endl;
        molecularTrajectoryFile = gmx_ffopen(molecularTrajectoryFilename.c_str(), "w");
    }

    if (frameType == SINGLE) {

        int frame = atoi(frameString);
        forceMatrix = parseScalarFileFormat(opt2fn("-ipf", NFILE, fnm), nbParticles, frame);
        if (opt2bSet("-ipf-diff", NFILE, fnm)) {
            forceMatrix2 = parseScalarFileFormat(opt2fn("-ipf-diff", NFILE, fnm), nbParticles, frame);
            for (int i = 0; i < nbParticles2; ++i) forceMatrix[i] -= forceMatrix2[i];
        }

        for (auto & f : forceMatrix) f = std::abs(f);

        // Convert from kJ/mol/nm into pN
        if (convert) for (auto & f : forceMatrix) f *= 1.66;

        BoostGraph graph(forceMatrix);
        BoostGraph::PathList shortestPaths = graph.findKShortestPaths(source, dest, numberOfShortestPaths);

        pdb.writePaths(opt2fn("-o", NFILE, fnm), shortestPaths, forceMatrix, false);

    } else {

        // Read trajectory coordinates
        t_trxstatus *status;
        real time;
        rvec *coord_traj;
        matrix box;

        int nbFrames = getNumberOfFrames(opt2fn("-ipf", NFILE, fnm));
        for (int frame = 0; frame < nbFrames; ++frame)
        {
            if (frame == 0) read_first_x(oenv, &status, opt2fn("-traj", NFILE, fnm), &time, &coord_traj, box);
            else read_next_x(oenv, status, &time, coord_traj, box);

            if (frameType == SKIP and frame%frameValue) continue;

			forceMatrix = parseScalarFileFormat(opt2fn("-ipf", NFILE, fnm), nbParticles, frame);
			if (opt2bSet("-ipf-diff", NFILE, fnm)) {
			    forceMatrix2 = parseScalarFileFormat(opt2fn("-ipf-diff", NFILE, fnm), nbParticles, frame);
			    for (int i = 0; i < nbParticles2; ++i) forceMatrix[i] -= forceMatrix2[i];
			}
    		for (auto & f : forceMatrix) f = std::abs(f);

    		// Convert from kJ/mol/nm into pN
    		if (convert) for (auto & f : forceMatrix) f *= 1.66;

            if (frameType == AVERAGE) {
                for (int frameAvg = 0; frameAvg < frameValue - 1; ++frameAvg)
                {
                    std::vector<double> forceMatrixAvg, forceMatrixAvg2;
                    forceMatrixAvg = parseScalarFileFormat(opt2fn("-ipf", NFILE, fnm), nbParticles, frame);
                    if (opt2bSet("-ipf-diff", NFILE, fnm)) {
                        forceMatrixAvg2 = parseScalarFileFormat(opt2fn("-ipf-diff", NFILE, fnm), nbParticles, frame);
                        for (int i = 0; i < nbParticles2; ++i) forceMatrixAvg[i] -= forceMatrixAvg2[i];
                    }

                    for (auto & f : forceMatrixAvg) f = std::abs(f);

                    // Convert from kJ/mol/nm into pN
                    if (convert) for (auto & f : forceMatrixAvg) f *= 1.66;

                    for (int i = 0; i < nbParticles2; ++i) forceMatrix[i] += forceMatrixAvg[i];
                }
                for (int i = 0; i < nbParticles2; ++i) forceMatrix[i] /= frameValue;
            }

    		BoostGraph graph(forceMatrix);
    		BoostGraph::PathList shortestPaths = graph.findKShortestPaths(source, dest, numberOfShortestPaths);

    		pdb.updateCoordinates(coord_traj);
    		pdb.writePaths(opt2fn("-o", NFILE, fnm), shortestPaths, forceMatrix, frame);

            // Write moleculare trajectory for pymol script
            if (opt2bSet("-pymol", NFILE, fnm))
                write_pdbfile(molecularTrajectoryFile, "FDA trajectory for Pymol visualization", &top.atoms, coord_traj, ePBC, box, ' ', 0, NULL, TRUE);

            if (frameType == AVERAGE) {
                frame += frameValue - 1;
                // It would be better to skip the next frameValue - 1 frames instead of reading them.
                for (int frameAvg = 0; frameAvg < frameValue - 1; ++frameAvg) {
                    read_next_x(oenv, status, &time, coord_traj, box);
                }
            }
        }
        close_trx(status);
    }

    if (opt2bSet("-pymol", NFILE, fnm)) gmx_ffclose(molecularTrajectoryFile);

    std::cout << "All done." << std::endl;
    return 0;

}
