/*
 * This source file is part of the Alexandria project.
 *
 * Copyright (C) 2014 David van der Spoel and Paul J. van Maaren
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "gromacs/utility/real.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/bonded/bonded.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/fileio/filenm.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/readinp.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atomprop.h"
#include "poldata.h"
#include "poldata_xml.h"
#include "molprop.h"
#include "molprop_util.h"
#include "molprop_xml.h"
#include "gauss_io.h"

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
#define NFILE asize(fnm)
    static gmx_bool                  bVerbose   = FALSE;
    static char                     *molnm      = NULL, *iupac = NULL, *conf = (char *)"minimum", *basis = NULL;
    static const char               *forcefield = "GAFF";
    static int                       maxpot     = 0;
    static gmx_bool                  compress   = FALSE;
    static gmx_bool                  bBabel     = TRUE;
    t_pargs                          pa[]       = {
        { "-v",      FALSE, etBOOL, {&bVerbose},
          "Generate verbose terminal output." },
        { "-compress", FALSE, etBOOL, {&compress},
          "Compress output XML files" },
#ifdef HAVE_LIBOPENBABEL2
        { "-babel", FALSE, etBOOL, {&bBabel},
          "Use the OpenBabel engine to process gaussian input files" },
#endif
        { "-molnm", FALSE, etSTR, {&molnm},
          "Name of the molecule in *all* input files. Do not use if you have different molecules in the input files." },
        { "-iupac", FALSE, etSTR, {&iupac},
          "IUPAC name of the molecule in *all* input files. Do not use if you have different molecules in the input files." },
        { "-conf",  FALSE, etSTR, {&conf},
          "Conformation of the molecule" },
        { "-basis",  FALSE, etSTR, {&basis},
          "Basis-set used in this calculation for those case where it is difficult to extract from a Gaussian file" },
        { "-ff", FALSE, etSTR, {&forcefield},
          "Force field for basic atom typing available in OpenBabel" },
        { "-maxpot", FALSE, etINT, {&maxpot},
          "Max number of potential points to add to the molprop file. If 0 all points are registered, else a selection of points evenly spread over the range of values is taken" }
    };
    output_env_t                     oenv;
    gmx_atomprop_t                   aps;
    gmx_poldata_t                    pd;
    std::vector<alexandria::MolProp> mp;
    alexandria::GaussAtomProp        gap;
    char **fns = NULL;
    int    i, nfn;

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, asize(pa), pa,
                           asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    /* Read standard atom properties */
    aps = gmx_atomprop_init();

    /* Read polarization stuff */
    if ((pd = gmx_poldata_read(NULL, aps)) == NULL)
    {
        gmx_fatal(FARGS, "Can not read the force field information. File missing or incorrect.");
    }

    nfn = ftp2fns(&fns, efLOG, NFILE, fnm);
    for (i = 0; (i < nfn); i++)
    {
        alexandria::MolProp mmm;

        ReadGauss(fns[i], mmm, gap, bBabel, aps, pd, molnm, iupac, conf, basis,
                  maxpot, bVerbose, gmx_poldata_get_force_field(pd));
        mp.push_back(mmm);
    }

    printf("Succesfully read %d molprops from %d Gaussian files.\n",
           (int)mp.size(), nfn);
    MolPropSort(mp, MPSA_MOLNAME, NULL, NULL);
    merge_doubles(mp, NULL, TRUE);
    if (mp.size() > 0)
    {
        MolPropWrite(opt2fn("-o", NFILE, fnm), mp, (int)compress);
    }

    return 0;
}
