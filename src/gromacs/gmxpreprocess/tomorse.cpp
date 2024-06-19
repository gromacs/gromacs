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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "tomorse.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <array>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/fileptr.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

typedef struct
{
    char *ai, *aj;
    real  e_diss;
} t_2morse;

static t_2morse* read_dissociation_energies(int* n2morse)
{
    char        ai[32], aj[32];
    double      e_diss;
    const char* fn     = "edissoc.dat";
    t_2morse*   t2m    = nullptr;
    int         maxn2m = 0, n2m = 0;
    int         nread;

    /* Open the file with dissociation energies */
    gmx::FilePtr fp = gmx::openLibraryFile(fn);
    do
    {
        /* Try and read two atom names and an energy term from it */
        nread = fscanf(fp.get(), "%s%s%lf", ai, aj, &e_diss);
        if (nread == 3)
        {
            /* If we got three terms, it probably was OK, no further checking */
            if (n2m >= maxn2m)
            {
                /* Increase memory for 16 at once, some mallocs are stupid */
                maxn2m += 16;
                srenew(t2m, maxn2m);
            }
            /* Copy the values */
            t2m[n2m].ai     = gmx_strdup(ai);
            t2m[n2m].aj     = gmx_strdup(aj);
            t2m[n2m].e_diss = e_diss;
            /* Increment counter */
            n2m++;
        }
        /* If we did not read three items, quit reading */
    } while (nread == 3);

    /* Set the return values */
    *n2morse = n2m;

    return t2m;
}

static int nequal(const char* a1, const char* a2)
{
    int i;

    /* Count the number of (case insensitive) characters that are equal in
     * two strings. If they are equally long their respective null characters are
     * counted also.
     */
    for (i = 0; (a1[i] != '\0') && (a2[i] != '\0'); i++)
    {
        if (toupper(a1[i]) != toupper(a2[i]))
        {
            break;
        }
    }
    if ((a1[i] == '\0') && (a2[i] == '\0'))
    {
        i++;
    }

    return i;
}

static real search_e_diss(int n2m, t_2morse t2m[], const char* ai, const char* aj)
{
    int  i;
    int  ibest = -1;
    int  nii, njj, nbstii = 0, nbstjj = 0;
    real ediss = 400;

    /* Do a best match search for dissociation energies */
    for (i = 0; (i < n2m); i++)
    {
        /* Check for a perfect match */
        if (((gmx_strcasecmp(t2m[i].ai, ai) == 0) && (gmx_strcasecmp(t2m[i].aj, aj) == 0))
            || ((gmx_strcasecmp(t2m[i].aj, ai) == 0) && (gmx_strcasecmp(t2m[i].ai, aj) == 0)))
        {
            ibest = i;
            break;
        }
        else
        {
            /* Otherwise count the number of equal characters in the strings ai and aj
             * and the ones from the file
             */
            nii = nequal(t2m[i].ai, ai);
            njj = nequal(t2m[i].aj, aj);
            if (((nii > nbstii) && (njj >= nbstjj)) || ((nii >= nbstii) && (njj > nbstjj)))
            {
                if ((nii > 0) && (njj > 0))
                {
                    ibest  = i;
                    nbstii = nii;
                    nbstjj = njj;
                }
            }
            else
            {
                /* Swap ai and aj (at least in counting the number of equal chars) */
                nii = nequal(t2m[i].ai, aj);
                njj = nequal(t2m[i].aj, ai);
                if (((nii > nbstii) && (njj >= nbstjj)) || ((nii >= nbstii) && (njj > nbstjj)))
                {
                    if ((nii > 0) && (njj > 0))
                    {
                        ibest  = i;
                        nbstii = nii;
                        nbstjj = njj;
                    }
                }
            }
        }
    }
    /* Return the dissocation energy corresponding to the best match, if we have
     * found one.
     */
    if (ibest == -1)
    {
        return ediss;
    }
    else
    {
        return t2m[ibest].e_diss;
    }
}

void convert_harmonics(gmx::ArrayRef<MoleculeInformation> mols, PreprocessingAtomTypes* atype)
{
    int       n2m;
    t_2morse* t2m;

    /* First get the data */
    t2m = read_dissociation_energies(&n2m);
    if (n2m <= 0)
    {
        fprintf(stderr, "No dissocation energies read\n");
        return;
    }

    /* For all the molecule types */
    int i = 0;
    for (auto& mol : mols)
    {
        /* Check how many morse and harmonic BONDSs there are, increase size of
         * morse with the number of harmonics
         */
        for (int bb = 0; (bb < F_NRE); bb++)
        {
            if ((interaction_function[bb].flags & IF_BTYPE) && (bb != F_MORSE))
            {
                int nrharm = mol.interactions[bb].size();

                /* Now loop over the harmonics, trying to convert them */
                for (auto harmonic = mol.interactions[bb].interactionTypes.begin();
                     harmonic != mol.interactions[bb].interactionTypes.end();)
                {
                    int  ni   = harmonic->ai();
                    int  nj   = harmonic->aj();
                    real edis = search_e_diss(
                            n2m,
                            t2m,
                            atype->atomNameFromAtomType(mol.atoms.atom[ni].type)->c_str(),
                            atype->atomNameFromAtomType(mol.atoms.atom[nj].type)->c_str());
                    if (edis != 0)
                    {
                        real              b0         = harmonic->c0();
                        real              kb         = harmonic->c1();
                        real              beta       = std::sqrt(kb / (2 * edis));
                        std::vector<int>  atoms      = { ni, nj };
                        std::vector<real> forceParam = { b0, edis, beta };
                        mol.interactions[F_MORSE].interactionTypes.emplace_back(
                                InteractionOfType(atoms, forceParam));
                        harmonic = mol.interactions[bb].interactionTypes.erase(harmonic);
                    }
                    else
                    {
                        ++harmonic;
                    }
                }

                int newHarmonics = mol.interactions[bb].size();
                fprintf(stderr,
                        "Converted %d out of %d %s to morse bonds for mol %d\n",
                        nrharm - newHarmonics,
                        nrharm,
                        interaction_function[bb].name,
                        i);
            }
        }
        i++;
    }
    sfree(t2m);
}
