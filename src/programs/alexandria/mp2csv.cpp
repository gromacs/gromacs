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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#include "molprop.h"
#include "molprop_util.h"
#include "molprop_xml.h"
#include "poldata_xml.h"

static void gmx_molprop_csv(const char *fn,
                            std::vector<alexandria::MolProp> mp,
                            const char *dip_str, const char *pol_str, const char *ener_str)
{
    alexandria::MolPropIterator mpi;
    FILE                       *fp;
    int                         k, ll;
    double                      T, d, err, vec[3];
    tensor                      quadrupole;
#define NEMP 3
    MolPropObservable           mpo[NEMP]   = { MPO_DIPOLE, MPO_POLARIZABILITY, MPO_ENERGY  };
    const char                 *ename[NEMP] = { "Dipole", "Polarizability", "Heat of formation" };
    alexandria::QmCount         qmc[NEMP];

    find_calculations(mp, mpo[0], dip_str, &qmc[0]);
    find_calculations(mp, mpo[1], pol_str, &qmc[1]);
    find_calculations(mp, mpo[2], ener_str, &qmc[2]);
    for (k = 0; (k < NEMP); k++)
    {
        printf("--------------------------------------------------\n");
        printf("      Some statistics for %s\n", mpo_name[mpo[k]]);
        for (auto q = qmc[k].beginCalc(); q < qmc[k].endCalc(); ++q)
        {
            printf("There are %d calculation results using %s type %s\n",
                   q->count(), q->lot().c_str(), q->type().c_str());
        }
    }
    printf("--------------------------------------------------\n");

    fp = gmx_ffopen(fn, "w");
    fprintf(fp, "\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\"",
            "Molecule", "Formula", "InChi", "Charge", "Multiplicity", "Mass");
    for (k = 0; (k < NEMP); k++)
    {
        for (size_t j = 0; (j < qmc[k].nCalc()+2); j++)
        {
            fprintf(fp, ",\"%s\"", ename[k]);
        }
    }
    fprintf(fp, "\n");
    for (ll = 0; (ll < NEMP); ll++)
    {
        fprintf(fp, "\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\"", "", "", "", "", "", "");
        for (k = 0; (k < 3); k++)
        {
            if (ll == 0)
            {
                fprintf(fp, ",\"Experiment\",\"Reference\"");
            }
            else
            {
                fprintf(fp, ",\"\",\"\"");
            }
            for (auto j = qmc[k].beginCalc(); j < qmc[k].endCalc(); ++j)
            {
                switch (ll)
                {
                    case 0:
                        fprintf(fp, ",\"%s\"", j->method().c_str());
                        break;
                    case 1:
                        fprintf(fp, ",\"%s\"", j->basis().c_str());
                        break;
                    case 2:
                        fprintf(fp, ",\"%s\"", j->type().c_str());
                        break;
                    default:
                        fprintf(stderr, "BOE\n");
                        exit(1);
                }
            }
        }
        fprintf(fp, "\n");
    }
    for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        fprintf(fp, "\"%s\",\"%s\",\"%s\",\"%d\",\"%d\",\"%g\"",
                mpi->getIupac().c_str(),
                mpi->formula().c_str(),
                mpi->getInchi().c_str(),
                mpi->getCharge(),
                mpi->getMultiplicity(),
                mpi->getMass());
        for (k = 0; (k < NEMP); k++)
        {
            std::string ref, mylot;
            if (mpi->getPropRef(mpo[k], iqmExp,
                                NULL, NULL, NULL, &d, &err, &T, ref, mylot, vec,
                                quadrupole))
            {
                fprintf(fp, ",\"%.4f\",\"%s\"", d, ref.c_str());
            }
            else
            {
                fprintf(fp, ",\"\",\"\"");
            }
            for (auto j = qmc[k].beginCalc(); j < qmc[k].endCalc(); j++)
            {
                if (mpi->getProp(mpo[k], iqmQM, j->lot(), NULL, j->type(),
                                 &T, &d, NULL))
                {
                    fprintf(fp, ",\"%.4f\"", d);
                }
                else
                {
                    fprintf(fp, ",\"\"");
                }
            }
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

int alex_mp2csv(int argc, char*argv[])
{
    static const char               *desc[] = {
        "mp2csv converts a molprop database into a spreadsheet"
    };
    t_filenm                         fnm[] =
    {
        { efDAT, "-f",  "allmols",  ffREAD },
        { efDAT, "-o",  "csvout",   ffWRITE }
    };
    int                              NFILE   = (sizeof(fnm)/sizeof(fnm[0]));
    static const char               *sort[]  = { NULL, "molname", "formula", "composition", NULL };
    static const char               *dip_str = "", *pol_str = "", *ener_str = "";
    static gmx_bool                  bMerge  = FALSE;
    t_pargs                          pa[]    =
    {
        { "-sort",   FALSE, etENUM, {sort},
          "Key to sort the final data file on." },
        { "-merge",  FALSE, etBOOL, {&bMerge},
          "Merge molecule records in the input file" },
        { "-dip_str", FALSE, etSTR, {&dip_str},
          "Selection of the dipole stuff you want in the tables, given as a single string with spaces like: method1/basis1/type1:method2/basis2/type2 (you may have to put quotes around the whole thing in order to prevent the shell from interpreting it)." },
        { "-pol_str", FALSE, etSTR, {&pol_str},
          "Same but for polarizabilities" },
        { "-ener_str", FALSE, etSTR, {&ener_str},
          "Same but for energies" }
    };
    std::vector<alexandria::MolProp> mp;
    gmx_atomprop_t                   ap;
    gmx_output_env_t                *oenv;

    if (!parse_common_args(&argc, argv, PCA_NOEXIT_ON_ARGS, NFILE, fnm,
                           sizeof(pa)/sizeof(pa[0]), pa,
                           sizeof(desc)/sizeof(desc[0]), desc,
                           0, NULL, &oenv))
    {
        return 0;
    }
    MolPropRead(opt2fn("-f", NFILE, fnm), mp);
    ap = gmx_atomprop_init();

    alexandria::MolSelect gms;
    MolPropSort(mp, MPSA_COMPOSITION, ap, gms);

    gmx_molprop_csv(opt2fn("-o", NFILE, fnm), mp,
                    strlen(dip_str) > 0  ? dip_str : NULL,
                    strlen(pol_str) > 0  ? pol_str : NULL,
                    strlen(ener_str) > 0 ? ener_str : NULL);

    return 0;
}
