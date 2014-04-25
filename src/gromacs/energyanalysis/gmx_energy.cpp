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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/gmxana/gstat.h"
#include "gmx_energy.h"
#include "select.h"
#include "dhdl.h"
#include "simple.h"
#include "viscosity.h"
#include "freeenergyestimate.h"
#include "nmr.h"
#include "fluctprops.h"
#include "handler.h"

int gmx_energy(int argc, char *argv[])
{
    const char        *desc[] = {
        "[THISMODULE] extracts energy components or distance restraint",
        "data from an energy file. The user is prompted to interactively",
        "select the desired energy terms.[PAR]",

        "Average, RMSD, and drift are calculated with full precision from the",
        "simulation (see printed manual). Drift is calculated by performing",
        "a least-squares fit of the data to a straight line. The reported total drift",
        "is the difference of the fit at the first and last point.",
        "An error estimate of the average is given based on a block averages",
        "over 5 blocks using the full-precision averages. The error estimate",
        "can be performed over multiple block lengths with the options",
        "[TT]-nbmin[tt] and [TT]-nbmax[tt].",
        "[BB]Note[bb] that in most cases the energy files contains averages over all",
        "MD steps, or over many more points than the number of frames in",
        "energy file. This makes the [THISMODULE] statistics output more accurate",
        "than the [TT].xvg[tt] output. When exact averages are not present in the energy",
        "file, the statistics mentioned above are simply over the single, per-frame",
        "energy values.[PAR]",

        "The term fluctuation gives the RMSD around the least-squares fit.[PAR]",

        "Some fluctuation-dependent properties can be calculated provided",
        "the correct energy terms are selected, and that the command line option",
        "[TT]-fluct_props[tt] is given. The following properties",
        "will be computed:[BR]",
        "Property                        Energy terms needed[BR]",
        "---------------------------------------------------[BR]",
        "Heat capacity C[SUB]p[sub] (NPT sims):    Enthalpy, Temp     [BR]",
        "Heat capacity C[SUB]v[sub] (NVT sims):    Etot, Temp         [BR]",
        "Thermal expansion coeff. (NPT): Enthalpy, Vol, Temp[BR]",
        "Isothermal compressibility:     Vol, Temp          [BR]",
        "Adiabatic bulk modulus:         Vol, Temp          [BR]",
        "---------------------------------------------------[BR]",
        "You always need to set the number of molecules [TT]-nmol[tt].",
        "The C[SUB]p[sub]/C[SUB]v[sub] computations do [BB]not[bb] include any corrections",
        "for quantum effects. Use the [gmx-dos] program if you need that (and you do).[PAR]"
        "When the [TT]-viol[tt] option is set, the time averaged",
        "violations are plotted and the running time-averaged and",
        "instantaneous sum of violations are recalculated. Additionally",
        "running time-averaged and instantaneous distances between",
        "selected pairs can be plotted with the [TT]-pairs[tt] option.[PAR]",

        "Options [TT]-ora[tt], [TT]-ort[tt], [TT]-oda[tt], [TT]-odr[tt] and",
        "[TT]-odt[tt] are used for analyzing orientation restraint data.",
        "The first two options plot the orientation, the last three the",
        "deviations of the orientations from the experimental values.",
        "The options that end on an 'a' plot the average over time",
        "as a function of restraint. The options that end on a 't'",
        "prompt the user for restraint label numbers and plot the data",
        "as a function of time. Option [TT]-odr[tt] plots the RMS",
        "deviation as a function of restraint.",
        "When the run used time or ensemble averaged orientation restraints,",
        "option [TT]-orinst[tt] can be used to analyse the instantaneous,",
        "not ensemble-averaged orientations and deviations instead of",
        "the time and ensemble averages.[PAR]",

        "Option [TT]-oten[tt] plots the eigenvalues of the molecular order",
        "tensor for each orientation restraint experiment. With option",
        "[TT]-ovec[tt] also the eigenvectors are plotted.[PAR]",

        "Option [TT]-odh[tt] extracts and plots the free energy data",
        "(Hamiltoian differences and/or the Hamiltonian derivative dhdl)",
        "from the [TT]ener.edr[tt] file.[PAR]",

        "With [TT]-fee[tt] an estimate is calculated for the free-energy",
        "difference with an ideal gas state: [BR]",
        "  [GRK]Delta[grk] A = A(N,V,T) - A[SUB]idealgas[sub](N,V,T) = kT [LN][CHEVRON][EXP]U[SUB]pot[sub]/kT[exp][chevron][ln][BR]",
        "  [GRK]Delta[grk] G = G(N,p,T) - G[SUB]idealgas[sub](N,p,T) = kT [LN][CHEVRON][EXP]U[SUB]pot[sub]/kT[exp][chevron][ln][BR]",
        "where k is Boltzmann's constant, T is set by [TT]-fetemp[tt] and",
        "the average is over the ensemble (or time in a trajectory).",
        "Note that this is in principle",
        "only correct when averaging over the whole (Boltzmann) ensemble",
        "and using the potential energy. This also allows for an entropy",
        "estimate using:[BR]",
        "  [GRK]Delta[grk] S(N,V,T) = S(N,V,T) - S[SUB]idealgas[sub](N,V,T) = ([CHEVRON]U[SUB]pot[sub][chevron] - [GRK]Delta[grk] A)/T[BR]",
        "  [GRK]Delta[grk] S(N,p,T) = S(N,p,T) - S[SUB]idealgas[sub](N,p,T) = ([CHEVRON]U[SUB]pot[sub][chevron] + pV - [GRK]Delta[grk] G)/T",
        "[PAR]",

        "When a second energy file is specified ([TT]-f2[tt]), a free energy",
        "difference is calculated [BR] dF = -kT [LN][CHEVRON][EXP]-(E[SUB]B[sub]-E[SUB]A[sub])/kT[exp][chevron][SUB]A[sub][ln] ,",
        "where E[SUB]A[sub] and E[SUB]B[sub] are the energies from the first and second energy",
        "files, and the average is over the ensemble A. The running average",
        "of the free energy difference is printed to a file specified by [TT]-ravg[tt].",
        "[BB]Note[bb] that the energies must both be calculated from the same trajectory."

    };
    static gmx_bool    bSum    = FALSE, bFee = FALSE, bPrAll = FALSE, bFluct = FALSE, bDriftCorr = FALSE;
    static gmx_bool    bDp     = FALSE, bMutot = FALSE, bOrinst = FALSE, bOvec = FALSE, bFluctProps = FALSE;
    static int         skip    = 0, nmol = 1, nbmin = 5, nbmax = 5;
    static real        reftemp = 300.0, ezero = 0;
    t_pargs            pa[]    = {
        { "-fee",   FALSE, etBOOL,  {&bFee},
          "Do a free energy estimate" },
        { "-fetemp", FALSE, etREAL, {&reftemp},
          "Reference temperature for free energy calculation" },
        { "-zero", FALSE, etREAL, {&ezero},
          "Subtract a zero-point energy" },
        { "-sum",  FALSE, etBOOL, {&bSum},
          "Sum the energy terms selected rather than display them all" },
        { "-dp",   FALSE, etBOOL, {&bDp},
          "Print energies in high precision" },
        { "-nbmin", FALSE, etINT, {&nbmin},
          "Minimum number of blocks for error estimate" },
        { "-nbmax", FALSE, etINT, {&nbmax},
          "Maximum number of blocks for error estimate" },
        { "-mutot", FALSE, etBOOL, {&bMutot},
          "Compute the total dipole moment from the components" },
        { "-skip", FALSE, etINT,  {&skip},
          "Skip number of frames between data points" },
        { "-aver", FALSE, etBOOL, {&bPrAll},
          "Also print the exact average and rmsd stored in the energy frames (only when 1 term is requested)" },
        { "-nmol", FALSE, etINT,  {&nmol},
          "Number of molecules in your sample: the energies are divided by this number" },
        { "-fluct_props", FALSE, etBOOL, {&bFluctProps},
          "Compute properties based on energy fluctuations, like heat capacity" },
        { "-driftcorr", FALSE, etBOOL, {&bDriftCorr},
          "Useful only for calculations of fluctuation properties. The drift in the observables will be subtracted before computing the fluctuation properties."},
        { "-fluc", FALSE, etBOOL, {&bFluct},
          "Calculate autocorrelation of energy fluctuations rather than energy itself" },
        { "-orinst", FALSE, etBOOL, {&bOrinst},
          "Analyse instantaneous orientation data" },
        { "-ovec", FALSE, etBOOL, {&bOvec},
          "Also plot the eigenvectors with [TT]-oten[tt]" }
    };

    t_filenm           fnm[] = {
        { efEDR, "-f",    NULL,      ffREAD  },
        { efEDR, "-f2",   NULL,      ffOPTRD },
        { efTPX, "-s",    NULL,      ffOPTRD },
        { efXVG, "-o",    "energy",  ffWRITE },
        { efXVG, "-fluct_conv", "fluct", ffOPTWR },
        { efXVG, "-viol", "violaver", ffOPTWR },
        { efXVG, "-pairs", "pairs",   ffOPTWR },
        { efXVG, "-ora",  "orienta", ffOPTWR },
        { efXVG, "-ort",  "orientt", ffOPTWR },
        { efXVG, "-oda",  "orideva", ffOPTWR },
        { efXVG, "-odr",  "oridevr", ffOPTWR },
        { efXVG, "-odt",  "oridevt", ffOPTWR },
        { efXVG, "-oten", "oriten",  ffOPTWR },
        { efXVG, "-corr", "enecorr", ffOPTWR },
        { efXVG, "-vis",  "visco",   ffOPTWR },
        { efXVG, "-evis", "evisco",  ffOPTWR },
        { efXVG, "-eivis", "eivisco", ffOPTWR },
        { efXVG, "-ravg", "runavgdf", ffOPTWR },
        { efXVG, "-odh",  "dhdl", ffOPTWR }
    };
#define NFILE asize(fnm)
    int                npargs;
    t_pargs           *ppa;
    output_env_t       oenv;

    npargs = asize(pa);
    ppa    = add_acf_pargs(&npargs, pa);
    if (!parse_common_args(&argc, argv,
                           PCA_CAN_VIEW | PCA_CAN_BEGIN | PCA_CAN_END | PCA_BE_NICE,
                           NFILE, fnm, npargs, ppa, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    gmx::EnergyHandler eh;
    eh.addEnergyFile(opt2fn("-f", NFILE, fnm));
    if (opt2bSet("-odh", NFILE, fnm))
    {
        // DHDL analysis
        gmx::DhdlEnergy *dh = new(gmx::DhdlEnergy);

        dh->helper()->setOutputFile(opt2fn("-odh", NFILE, fnm));
        dh->setParameters(opt2fn("-s", NFILE, fnm));
        dh->helper()->setDoublePrecision(bDp);
        dh->helper()->setOutputEnvironment(oenv);
        eh.addAnalysisTool(dh);
    }
    else if (opt2bSet("-f2", NFILE, fnm))
    {
        // Free energy estimate
        gmx::FreeEnergyEstimate *fee = new(gmx::FreeEnergyEstimate);

        eh.addEnergyFile(opt2fn("-f2", NFILE, fnm));
        fee->helper()->setOutputFile(opt2fn("-o", NFILE, fnm));
        fee->setParameters(reftemp);
        eh.addAnalysisTool(fee);
    }
    else if (opt2parg_bSet("-fluct_props", npargs, ppa))
    {
        // Fluctuation property analysis
        gmx::FluctProps *fp = new(gmx::FluctProps);

        fp->helper()->setOutputEnvironment(oenv);
        fp->helper()->setNmol(nmol);
        fp->helper()->setNblocks(nbmin);
        const char *fc = opt2fn_null("-fluct_conv", NFILE, fnm);
        if (NULL != fc)
        {
            fp->setFluctConvFile(fc);
        }
        eh.addAnalysisTool(fp);
    }
    else if (opt2bSet("-vis", NFILE, fnm) ||
             opt2bSet("-evis", NFILE, fnm) ||
             opt2bSet("-eivis", NFILE, fnm))
    {
        const char     *fn;
        // Viscosity analysis
        gmx::Viscosity *vis = new(gmx::Viscosity);

        vis->helper()->setOutputEnvironment(oenv);
        fn = opt2fn_null("-vis", NFILE, fnm);
        if (NULL != fn)
        {
            vis->helper()->setOutputFile(fn);
        }
        fn = opt2fn_null("-evis", NFILE, fnm);
        if (NULL != fn)
        {
            vis->setFnEinstein(fn);
        }
        fn = opt2fn_null("-eivis", NFILE, fnm);
        if (NULL != fn)
        {
            vis->setFnEinsteinIntegral(fn);
        }
        vis->helper()->setNblocks(nbmin);
        eh.addAnalysisTool(vis);
    }
    else
    {
        // Simple energy analysis
        gmx::SimpleEnergy *se = new(gmx::SimpleEnergy);

        se->helper()->setOutputEnvironment(oenv);
        se->helper()->setNmol(nmol);
        se->helper()->setOutputFile(opt2fn("-o", NFILE, fnm));
        se->helper()->setStoreData(bSum);
        se->helper()->setDoublePrecision(bDp);
        se->helper()->setNblocks(nbmin);
        se->setSumming(bSum);
        const char *fc = opt2fn_null("-fluct_conv", NFILE, fnm);
        if (NULL != fc)
        {
            se->setFluctConvFile(fc);
        }
        eh.addAnalysisTool(se);
    }

    if (!eh.readFiles())
    {
        fprintf(stderr, "Too bad. Something went wrong reading %s.\n",
                opt2fn("-f", NFILE, fnm));
    }

    return 0;
}
