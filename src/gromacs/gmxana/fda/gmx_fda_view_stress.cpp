/*
 * gmx_fda_view_stress.cpp
 *
 *  Created on: Feb 13, 2015
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/gmxana/fda/EnumParser.h"
#include "gromacs/gmxana/fda/FrameType.h"
#include "gromacs/gmxana/fda/Helpers.h"
#include "gromacs/gmxana/fda/ParticleType.h"
#include "gromacs/gmxana/fda/StressType.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

using namespace fda_analysis;

#define PRINT_DEBUG

int gmx_fda_view_stress(int argc, char *argv[])
{
    const char *desc[] = {
        "[THISMODULE] illustrates the punctual or von Mises virial stress of FDA "
        "as xpm or pdb-file. "
        "For the xpm-file the number of different colors can be set with the "
        "[TT]-nbColors[tt] option. "
        "The x-axis of the xpm file represent the particle number and the y-axis the "
        "frame number. "
        "For the pdb-file the Bfactor column will be used for the value of "
        "the stress and helps the coloring as a function of the stress magnitude. "
        "The pdb-file can be visualized with an program of your choice. "
    };

    gmx_output_env_t *oenv;
    static const char* frameString = "average 1";
    static bool convert = false;
    static int nbColors = 10;

    t_pargs pa[] = {
        { "-frame", FALSE, etSTR, {&frameString}, "Specify a single frame number or \"average n\" to take the mean over every n-th frame"
              " or \"skip n\" to take every n-th frame or \"all\" to take all frames (e.g. for movies)" },
        { "-convert", FALSE, etBOOL, {&convert}, "Convert force unit from kJ/mol/nm into pN (only for punctual stress)" },
        { "-nbColors", FALSE, etINT, {&nbColors}, "Number of colors for xpm-files" }
    };

    t_filenm fnm[] = {
		{ efSTR, NULL, NULL, ffREAD },
        { efTPS, NULL, NULL, ffOPTRD },
        { efTRX, "-traj", NULL, ffOPTRD },
        { efNDX, NULL, NULL, ffOPTRD },
        { efVST, "-o", "result", ffWRITE }
    };

#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, PCA_CAN_TIME,
        NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

	StressType stressType = PUNCTUAL;
	ParticleType particleType = ATOM;
    if (fn2ftp(opt2fn("-f", NFILE, fnm)) == efPSA) {
        stressType = PUNCTUAL;
        particleType = ATOM;
    } else if (fn2ftp(opt2fn("-f", NFILE, fnm)) == efPSR) {
        stressType = PUNCTUAL;
        particleType = RESIDUE;
    } else if (fn2ftp(opt2fn("-f", NFILE, fnm)) == efVMA) {
        stressType = VIRIAL_VON_MISES;
        particleType = ATOM;
    } else {
    	gmx_fatal(FARGS, "Unkown stress or particle type of stress file.");
    }

    int frameValue = 0;
    FrameType frameType = getFrameTypeAndSkipValue(frameString, frameValue);

    int nbFrames, nbParticles;
    std::vector<real> stressMatrix = readStress(opt2fn("-f", NFILE, fnm), nbFrames, nbParticles);

    // Interactive input of group name for residue model points
    int isize = 0;
    int *index;
    char *grpname;
    if (ftp2fn_null(efNDX, NFILE, fnm)) {
		fprintf(stderr, "\nSelect group for residue model points:\n");
		rd_index(ftp2fn(efNDX, NFILE, fnm), 1, &isize, &index, &grpname);
        if (isize != nbParticles) gmx_fatal(FARGS, "Number of atoms in group does not match number of FDA points.");
    }

	#ifdef PRINT_DEBUG
	    std::cerr << "stress filename = " << opt2fn("-f", NFILE, fnm) << std::endl;
        std::cerr << "result filename = " << opt2fn("-o", NFILE, fnm) << std::endl;
	    std::cerr << "frameType = " << EnumParser<FrameType>()(frameType) << std::endl;
		std::cerr << "frameValue = " << frameValue << std::endl;
	    std::cerr << "stressType = " << EnumParser<StressType>()(stressType) << std::endl;
	    std::cerr << "particleType = " << EnumParser<ParticleType>()(particleType) << std::endl;
	    std::cerr << "convert = " << convert << std::endl;
    #endif

    if (fn2ftp(opt2fn("-o", NFILE, fnm)) == efPDB and !ftp2fn_null(efTPS, NFILE, fnm)) gmx_fatal(FARGS, "Input structure is missing.");
    if (stressType != PUNCTUAL and convert) gmx_fatal(FARGS, "Option -convert makes only sense for punctual stress.");

	for (auto & elem : stressMatrix) elem = std::abs(elem);

	// Convert from kJ/mol/nm into pN
	if (convert) for (auto & elem : stressMatrix) elem *= 1.66;

	std::string title;
	if (stressType == PUNCTUAL) title += "Punctual ";
	else if (stressType == VIRIAL) title += "Virial ";
	else if (stressType == VIRIAL_VON_MISES) title += "Von Mises virial ";
	title += "stress";
	if (particleType == ATOM) title += " over atoms";
	else if (particleType == RESIDUE) title += " over residues";

	bool valueToLargeForPDB = false;

	if (fn2ftp(opt2fn("-o", NFILE, fnm)) == efPDB) {

        rvec *xp;
        t_topology top;
        int ePBC;
        matrix box;

        read_tps_conf(ftp2fn(efTPS, NFILE, fnm), &top, &ePBC, &xp, NULL, box, TRUE);

        real currentStress;

		if (frameType == SINGLE) {

			FILE *fp = gmx_ffopen(opt2fn("-o", NFILE, fnm), "w");

			for (int i = 0; i < nbParticles; ++i) {
			    currentStress = stressMatrix[frameValue*nbParticles + i];
			    if (currentStress > 999.99) {
			        top.atoms.pdbinfo[i].bfac = 999.99;
			        valueToLargeForPDB = true;
			    }
			    else top.atoms.pdbinfo[i].bfac = currentStress;
			}

			write_pdbfile(fp, title.c_str(), &top.atoms, xp, ePBC, box, ' ', 0, NULL, TRUE);
			gmx_ffclose(fp);

		} else {

		    // Read trajectory coordinates
	        t_trxstatus *status;
	        real time;
	        rvec *coord_traj;
	        matrix box;
	        read_first_x(oenv, &status, opt2fn("-traj", NFILE, fnm), &time, &coord_traj, box);

	        for (int frame = 0; frame < nbFrames; ++frame)
	        {
                read_next_x(oenv, status, &time, coord_traj, box);
	            if (frame%frameValue) continue;

	            FILE *fp = NULL;
	            if (frame) fp = gmx_ffopen(opt2fn("-o", NFILE, fnm), "a");
	            else fp = gmx_ffopen(opt2fn("-o", NFILE, fnm), "w");

	            for (int i = 0; i < nbParticles; ++i) {
	                currentStress = stressMatrix[frameValue*nbParticles + i];
	                top.atoms.pdbinfo[i].bfac = currentStress;
	                if (currentStress > 999.99) {
	                    valueToLargeForPDB = true;
	                }
	            }

	            write_pdbfile(fp, title.c_str(), &top.atoms, coord_traj, ePBC, box, ' ', 0, NULL, TRUE);
	            gmx_ffclose(fp);
	        }
	        close_trx(status);
		}
	} else if (fn2ftp(opt2fn("-o", NFILE, fnm)) == efXPM) {

	    // Reorder stressMatrix for writing xpm
	    real **stressMatrix2 = NULL;
	    snew(stressMatrix2, nbParticles);
	    real minValue = std::numeric_limits<real>::max();
	    real maxValue = 0;
	    int nbFramesForOutput = 1;
	    if (frameType == SINGLE) {
	        for (int i = 0; i < nbParticles; ++i) {
	            snew(stressMatrix2[i], nbFramesForOutput);
	            real value = stressMatrix[frameValue*nbParticles + i];
	            stressMatrix2[i][0] = value;
	            if (value < minValue) minValue = value;
	            if (value > maxValue) maxValue = value;
	        }
	    } else {
	        nbFramesForOutput = ceil(static_cast<real>(nbFrames) / frameValue);
	        for (int i = 0; i < nbParticles; ++i) {
	            snew(stressMatrix2[i], nbFramesForOutput);
	            if (frameType == AVERAGE) {
	                for (int j = 0, js = 0; j < nbFramesForOutput; ++j) {
	                    real value = 0.0;
	                    for (int k = 0; k < frameValue and js < nbFrames; ++k, ++js) value += stressMatrix[js*nbParticles + i];
	                    value /= frameValue;
	                    stressMatrix2[i][j] = value;
	                    if (value < minValue) minValue = value;
	                    if (value > maxValue) maxValue = value;
	                }
	            } else {
	                for (int j = 0, js = 0; j < nbFramesForOutput; ++j, js += frameValue) {
	                    real value = stressMatrix[js*nbParticles + i];
	                    stressMatrix2[i][j] = value;
	                    if (value < minValue) minValue = value;
	                    if (value > maxValue) maxValue = value;
	                }
	            }
	        }
	    }

        FILE *out = gmx_ffopen(opt2fn("-o", NFILE, fnm), "w");
        t_rgb rlo = {1, 1, 1}, rhi = {0, 0, 0};
        write_xpm(out, 0, title.c_str(), "", "Particle", "Frame", nbParticles, nbFramesForOutput,
            NULL, NULL, stressMatrix2, minValue, maxValue, rlo, rhi, &nbColors);
        gmx_ffclose(out);

	} else gmx_fatal(FARGS, "Missing output filename -opdb or -oxpm.");

    if (valueToLargeForPDB)
        gmx_warning("Stress values larger than 999.99 are detected. Therefore, the general PDB format of the b-factor column of Real(6.2) is broken. "
                    "It is tested that it works for Pymol and VMD, but it is not guaranteed that it will work for other visualization programs.");

    std::cout << "All done." << std::endl;
    return 0;

}
