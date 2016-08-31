/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "gauss_io.h"
#include "molprop.h"
#include "molprop_util.h"
#include "molprop_xml.h"
#include "poldata.h"
#include "poldata_xml.h"

int alex_gauss2molprop(int argc, char *argv[])
{
    static const char               *desc[] = {
        "gauss2molprop reads a series of Gaussian output files, and collects",
        "useful information, and saves it to molprop file."
    };

    t_filenm                         fnm[] = {
        { efLOG, "-g03",  "gauss",  ffRDMULT },
        { efDAT, "-o",    "molprop", ffWRITE }
    };
#define NFILE sizeof(fnm)/sizeof(fnm[0])
    static gmx_bool                  bVerbose   = FALSE;
    static char                     *molnm      = NULL, *iupac = NULL, *jobtype = (char *)"Opt";
    static char                     *conf       = (char *)"minimum", *basis = NULL;
    static const char               *forcefield = "GAFF";
    static int                       maxpot     = 0;
    static int                       nsymm      = 0;
    static gmx_bool                  compress   = FALSE;
    t_pargs                          pa[]       = {
        { "-v",      FALSE, etBOOL, {&bVerbose},
          "Generate verbose terminal output." },
        { "-compress", FALSE, etBOOL, {&compress},
          "Compress output XML files" },
        { "-molnm", FALSE, etSTR, {&molnm},
          "Name of the molecule in *all* input files. Do not use if you have different molecules in the input files." },
        { "-nsymm", FALSE, etINT, {&nsymm},
          "Symmetry number of the molecule can be supplied here if you know there is an error in the input file" },
        { "-iupac", FALSE, etSTR, {&iupac},
          "IUPAC name of the molecule in *all* input files. Do not use if you have different molecules in the input files." },
        { "-conf",  FALSE, etSTR, {&conf},
          "Conformation of the molecule" },
        { "-basis",  FALSE, etSTR, {&basis},
          "Basis-set used in this calculation for those case where it is difficult to extract from a Gaussian file" },
        { "-jobtype",  FALSE, etSTR, {&jobtype},
          "The job type used in the Gaussian calculation: Opt, Polar, SP, and etc." },
        { "-ff", FALSE, etSTR, {&forcefield},
          "Force field for basic atom typing available in OpenBabel" },
        { "-maxpot", FALSE, etINT, {&maxpot},
          "Max number of potential points to add to the molprop file. If 0 all points are registered, else a selection of points evenly spread over the range of values is taken" }
    };
    gmx_output_env_t                *oenv;
    gmx_atomprop_t                   aps;
    alexandria::Poldata              pd;
    std::vector<alexandria::MolProp> mp;
    alexandria::GaussAtomProp        gap;
    char **fns = NULL;
    int    i, nfn;

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm,
                           sizeof(pa)/sizeof(pa[0]), pa,
                           sizeof(desc)/sizeof(desc[0]), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    /* Read standard atom properties */
    aps = gmx_atomprop_init();

    /* Read force field stuff */
    try
    {
        readPoldata("", pd, aps);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    nfn = ftp2fns(&fns, efLOG, NFILE, fnm);
    for (i = 0; (i < nfn); i++)
    {
        alexandria::MolProp mmm;

        ReadGauss(fns[i], mmm, molnm, iupac, conf, basis,
                  maxpot, nsymm, pd.getForceField().c_str(), jobtype);
        mp.push_back(mmm);
    }

    printf("Succesfully read %d molprops from %d Gaussian files.\n",
           (int)mp.size(), nfn);
    alexandria::MolSelect gms;
    MolPropSort(mp, MPSA_MOLNAME, NULL, gms);
    merge_doubles(mp, NULL, TRUE);
    if (mp.size() > 0)
    {
        MolPropWrite(opt2fn("-o", NFILE, fnm), mp, (int)compress);
    }

    return 0;
}
