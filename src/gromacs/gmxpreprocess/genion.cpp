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
#include "gmxpre.h"

#include "genion.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <iterator>
#include <numeric>
#include <string>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/random/seed.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformintdistribution.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

enum class PbcType : int;
struct gmx_output_env_t;


/*! \brief Return whether any atoms of two groups are below minimum distance.
 *
 * \param[in] pbc the periodic boundary conditions
 * \param[in] x the coordinates
 * \param[in] smallerGroupIndices the atom indices of the first group to check
 * \param[in] largerGroupIndices the atom indices of the second group to check
 * \param[in] minimumDistance the minimum required distance betwenn any atom in
 *                            the first and second group
 * \returns true if any distance between an atom from group A and group B is
 *               smaller than a minimum distance.
 */
static bool groupsCloserThanCutoffWithPbc(t_pbc*                   pbc,
                                          rvec                     x[],
                                          gmx::ArrayRef<const int> smallerGroupIndices,
                                          gmx::ArrayRef<const int> largerGroupIndices,
                                          real                     minimumDistance)
{
    const real minimumDistance2 = minimumDistance * minimumDistance;
    for (int aIndex : largerGroupIndices)
    {
        for (int bIndex : smallerGroupIndices)
        {
            rvec dx;
            pbc_dx(pbc, x[aIndex], x[bIndex], dx);
            if (norm2(dx) < minimumDistance2)
            {
                return true;
            }
        }
    }
    return false;
}

/*! \brief Calculate the solvent molecule atom indices from molecule number.
 *
 * \note the solvent group index has to be continuous
 *
 * \param[in] solventMoleculeNumber the number of the solvent molecule
 * \param[in] numberAtomsPerSolventMolecule how many atoms each solvent molecule contains
 * \param[in] solventGroupIndex continuous index of solvent atoms
 *
 * \returns atom indices of the specified solvent molecule
 */
static std::vector<int> solventMoleculeIndices(int solventMoleculeNumber,
                                               int numberAtomsPerSolventMolecule,
                                               gmx::ArrayRef<const int> solventGroupIndex)
{
    std::vector<int> indices(numberAtomsPerSolventMolecule);
    for (int solventAtomNumber = 0; solventAtomNumber < numberAtomsPerSolventMolecule; ++solventAtomNumber)
    {
        indices[solventAtomNumber] =
                solventGroupIndex[numberAtomsPerSolventMolecule * solventMoleculeNumber + solventAtomNumber];
    }
    return indices;
}

static void insert_ion(int                      nsa,
                       std::vector<int>*        solventMoleculesForReplacement,
                       int                      repl[],
                       gmx::ArrayRef<const int> index,
                       rvec                     x[],
                       t_pbc*                   pbc,
                       int                      sign,
                       int                      q,
                       const char*              ionname,
                       t_atoms*                 atoms,
                       real                     rmin,
                       std::vector<int>*        notSolventGroup)
{
    std::vector<int> solventMoleculeAtomsToBeReplaced =
            solventMoleculeIndices(solventMoleculesForReplacement->back(), nsa, index);

    if (rmin > 0.0)
    {
        // check for proximity to non-solvent
        while (groupsCloserThanCutoffWithPbc(pbc, x, solventMoleculeAtomsToBeReplaced, *notSolventGroup, rmin)
               && !solventMoleculesForReplacement->empty())
        {
            solventMoleculesForReplacement->pop_back();
            solventMoleculeAtomsToBeReplaced =
                    solventMoleculeIndices(solventMoleculesForReplacement->back(), nsa, index);
        }
    }

    if (solventMoleculesForReplacement->empty())
    {
        gmx_fatal(FARGS, "No more replaceable solvent!");
    }

    fprintf(stderr,
            "Replacing solvent molecule %d (atom %d) with %s\n",
            solventMoleculesForReplacement->back(),
            solventMoleculeAtomsToBeReplaced[0],
            ionname);

    /* Replace solvent molecule charges with ion charge */
    notSolventGroup->push_back(solventMoleculeAtomsToBeReplaced[0]);
    repl[solventMoleculesForReplacement->back()] = sign;

    // The first solvent molecule atom is replaced with an ion and the respective
    // charge while the rest of the solvent molecule atoms is set to 0 charge.
    atoms->atom[solventMoleculeAtomsToBeReplaced.front()].q = q;
    for (auto replacedMoleculeAtom = solventMoleculeAtomsToBeReplaced.begin() + 1;
         replacedMoleculeAtom != solventMoleculeAtomsToBeReplaced.end();
         ++replacedMoleculeAtom)
    {
        atoms->atom[*replacedMoleculeAtom].q = 0;
    }
    solventMoleculesForReplacement->pop_back();
}


static char* aname(const char* mname)
{
    char* str;
    int   i;

    str = gmx_strdup(mname);
    i   = std::strlen(str) - 1;
    while (i > 1 && (std::isdigit(str[i]) || (str[i] == '+') || (str[i] == '-')))
    {
        str[i] = '\0';
        i--;
    }

    return str;
}

static void sort_ions(int                      nsa,
                      int                      nw,
                      const int                repl[],
                      gmx::ArrayRef<const int> index,
                      t_atoms*                 atoms,
                      rvec                     x[],
                      char**                   pptr,
                      char**                   nptr,
                      char**                   paptr,
                      char**                   naptr)
{
    int   i, j, k, r, np, nn, starta, startr, npi, nni;
    rvec* xt;

    snew(xt, atoms->nr);

    /* Put all the solvent in front and count the added ions */
    np = 0;
    nn = 0;
    j  = index[0];
    for (i = 0; i < nw; i++)
    {
        r = repl[i];
        if (r == 0)
        {
            for (k = 0; k < nsa; k++)
            {
                copy_rvec(x[index[nsa * i + k]], xt[j++]);
            }
        }
        else if (r > 0)
        {
            np++;
        }
        else if (r < 0)
        {
            nn++;
        }
    }

    if (np + nn > 0)
    {
        /* Put the positive and negative ions at the end */
        starta = index[nsa * (nw - np - nn)];
        startr = atoms->atom[starta].resind;

        npi = 0;
        nni = 0;
        for (i = 0; i < nw; i++)
        {
            r = repl[i];
            if (r > 0)
            {
                j = starta + npi;
                k = startr + npi;
                copy_rvec(x[index[nsa * i]], xt[j]);
                atoms->atomname[j]     = paptr;
                atoms->atom[j].resind  = k;
                atoms->resinfo[k].name = pptr;
                npi++;
            }
            else if (r < 0)
            {
                j = starta + np + nni;
                k = startr + np + nni;
                copy_rvec(x[index[nsa * i]], xt[j]);
                atoms->atomname[j]     = naptr;
                atoms->atom[j].resind  = k;
                atoms->resinfo[k].name = nptr;
                nni++;
            }
        }
        for (i = index[nsa * nw - 1] + 1; i < atoms->nr; i++)
        {
            j                  = i - (nsa - 1) * (np + nn);
            atoms->atomname[j] = atoms->atomname[i];
            atoms->atom[j]     = atoms->atom[i];
            copy_rvec(x[i], xt[j]);
        }
        atoms->nr -= (nsa - 1) * (np + nn);

        /* Copy the new positions back */
        for (i = index[0]; i < atoms->nr; i++)
        {
            copy_rvec(xt[i], x[i]);
        }
        sfree(xt);
    }
}

static void update_topol(const char* topinout, int p_num, int n_num, const char* p_name, const char* n_name, char* grpname)
{
    FILE *   fpin, *fpout;
    char     buf[STRLEN], buf2[STRLEN], *temp, **mol_line = nullptr;
    int      i, nmol_line, sol_line, nsol_last;
    gmx_bool bMolecules;
    char     temporary_filename[STRLEN];

    printf("\nProcessing topology\n");
    fpin = gmx_ffopen(topinout, "r");
    std::strncpy(temporary_filename, "temp.topXXXXXX", STRLEN);
    fpout = gmx_fopen_temporary(temporary_filename);

    bMolecules = FALSE;
    nmol_line  = 0;
    sol_line   = -1;
    nsol_last  = -1;
    while (fgets(buf, STRLEN, fpin))
    {
        std::strcpy(buf2, buf);
        if ((temp = std::strchr(buf2, '\n')) != nullptr)
        {
            temp[0] = '\0';
        }
        ltrim(buf2);
        if (buf2[0] == '[')
        {
            buf2[0] = ' ';
            if ((temp = std::strchr(buf2, '\n')) != nullptr)
            {
                temp[0] = '\0';
            }
            rtrim(buf2);
            if (buf2[std::strlen(buf2) - 1] == ']')
            {
                buf2[std::strlen(buf2) - 1] = '\0';
                ltrim(buf2);
                rtrim(buf2);
                bMolecules = (gmx_strcasecmp(buf2, "molecules") == 0);
            }
            fprintf(fpout, "%s", buf);
        }
        else if (!bMolecules)
        {
            fprintf(fpout, "%s", buf);
        }
        else
        {
            /* Check if this is a line with solvent molecules */
            sscanf(buf, "%s", buf2);
            if (gmx_strcasecmp(buf2, grpname) == 0)
            {
                sol_line = nmol_line;
                sscanf(buf, "%*s %d", &nsol_last);
            }
            /* Store this molecules section line */
            srenew(mol_line, nmol_line + 1);
            mol_line[nmol_line] = gmx_strdup(buf);
            nmol_line++;
        }
    }
    gmx_ffclose(fpin);

    if (sol_line == -1)
    {
        gmx_ffclose(fpout);
        gmx_fatal(FARGS,
                  "No line with moleculetype '%s' found the [ molecules ] section of file '%s'",
                  grpname,
                  topinout);
    }
    if (nsol_last < p_num + n_num)
    {
        gmx_ffclose(fpout);
        gmx_fatal(FARGS,
                  "The last entry for moleculetype '%s' in the [ molecules ] section of file '%s' "
                  "has less solvent molecules (%d) than were replaced (%d)",
                  grpname,
                  topinout,
                  nsol_last,
                  p_num + n_num);
    }

    /* Print all the molecule entries */
    for (i = 0; i < nmol_line; i++)
    {
        if (i != sol_line)
        {
            fprintf(fpout, "%s", mol_line[i]);
        }
        else
        {
            printf("Replacing %d solute molecules in topology file (%s) "
                   " by %d %s and %d %s ions.\n",
                   p_num + n_num,
                   topinout,
                   p_num,
                   p_name,
                   n_num,
                   n_name);
            nsol_last -= p_num + n_num;
            if (nsol_last > 0)
            {
                fprintf(fpout, "%-10s  %d\n", grpname, nsol_last);
            }
            if (p_num > 0)
            {
                fprintf(fpout, "%-15s  %d\n", p_name, p_num);
            }
            if (n_num > 0)
            {
                fprintf(fpout, "%-15s  %d\n", n_name, n_num);
            }
        }
    }
    gmx_ffclose(fpout);
    make_backup(topinout);
    gmx_file_rename(temporary_filename, topinout);
}

/*! \brief Return all atom indices that do not belong to an index group.
 * \param[in] nrAtoms the total number of atoms
 * \param[in] indexGroup the index group to be inverted. Note that a copy is
 *                       made in order to sort the group indices
 * \returns the inverted group indices
 */
static std::vector<int> invertIndexGroup(int nrAtoms, std::vector<int> indexGroup)
{
    // Add the indices -1 and nrAtoms, so that all intervals 0 .. firstofIndexGroup
    // as well as lastOfIndexGroup .. nrAtoms are covered for the inverted indexgroup
    indexGroup.push_back(-1);
    indexGroup.push_back(nrAtoms);
    std::sort(indexGroup.begin(), indexGroup.end());

    // construct the inverted index group by adding all indicies between two
    // indices of indexGroup
    std::vector<int> invertedGroup;
    for (auto indexGroupIt = std::begin(indexGroup); indexGroupIt != std::end(indexGroup) - 1; ++indexGroupIt)
    {
        const int firstToAddToInvertedGroup = *indexGroupIt + 1;
        const int numIndicesToAdd           = *(indexGroupIt + 1) - firstToAddToInvertedGroup;
        if (numIndicesToAdd > 0)
        {
            invertedGroup.resize(invertedGroup.size() + numIndicesToAdd);
            std::iota(std::end(invertedGroup) - numIndicesToAdd, std::end(invertedGroup), firstToAddToInvertedGroup);
        }
    }

    return invertedGroup;
}

int gmx_genion(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] randomly replaces solvent molecules with monoatomic ions.",
        "The group of solvent molecules should be continuous and all molecules",
        "should have the same number of atoms.",
        "The user should add the ion molecules to the topology file or use",
        "the [TT]-p[tt] option to automatically modify the topology.[PAR]",
        "The ion molecule type, residue and atom names in all force fields",
        "are the capitalized element names without sign. This molecule name",
        "should be given with [TT]-pname[tt] or [TT]-nname[tt], and the",
        "[TT][molecules][tt] section of your topology updated accordingly,",
        "either by hand or with [TT]-p[tt]. Do not use an atom name instead!",
        "[PAR]Ions which can have multiple charge states get the multiplicity",
        "added, without sign, for the uncommon states only.[PAR]",
        "For larger ions, e.g. sulfate we recommended using [gmx-insert-molecules]."
    };
    const char* bugs[] = {
        "If you specify a salt concentration existing ions are not taken into "
        "account. In effect you therefore specify the amount of salt to be added.",
    };
    int         p_num = 0, n_num = 0, p_q = 1, n_q = -1;
    const char *p_name = "NA", *n_name = "CL";
    real        rmin = 0.6, conc = 0;
    int         seed     = 0;
    gmx_bool    bNeutral = FALSE;
    t_pargs     pa[]     = {
        { "-np", FALSE, etINT, { &p_num }, "Number of positive ions" },
        { "-pname", FALSE, etSTR, { &p_name }, "Name of the positive ion" },
        { "-pq", FALSE, etINT, { &p_q }, "Charge of the positive ion" },
        { "-nn", FALSE, etINT, { &n_num }, "Number of negative ions" },
        { "-nname", FALSE, etSTR, { &n_name }, "Name of the negative ion" },
        { "-nq", FALSE, etINT, { &n_q }, "Charge of the negative ion" },
        { "-rmin", FALSE, etREAL, { &rmin }, "Minimum distance between ions and non-solvent" },
        { "-seed", FALSE, etINT, { &seed }, "Seed for random number generator (0 means generate)" },
        { "-conc",
          FALSE,
          etREAL,
          { &conc },
          "Specify salt concentration (mol/liter). This will add sufficient ions to reach up to "
          "the specified concentration as computed from the volume of the cell in the input "
          "[REF].tpr[ref] file. Overrides the [TT]-np[tt] and [TT]-nn[tt] options." },
        { "-neutral",
          FALSE,
          etBOOL,
          { &bNeutral },
          "This option will add enough ions to neutralize the system. These ions are added on top "
          "of those specified with [TT]-np[tt]/[TT]-nn[tt] or [TT]-conc[tt]. " }
    };
    t_topology        top;
    rvec*             x;
    real              vol;
    matrix            box;
    t_atoms           atoms;
    t_pbc             pbc;
    int*              repl;
    PbcType           pbcType;
    int               nw, nsa, nsalt, iqtot;
    gmx_output_env_t* oenv  = nullptr;
    t_filenm          fnm[] = { { efTPR, nullptr, nullptr, ffREAD },
                       { efNDX, nullptr, nullptr, ffOPTRD },
                       { efSTO, "-o", nullptr, ffWRITE },
                       { efTOP, "-p", "topol", ffOPTRW } };
#define NFILE asize(fnm)

    if (!parse_common_args(
                &argc, argv, 0, NFILE, fnm, asize(pa), pa, asize(desc), desc, asize(bugs), bugs, &oenv))
    {
        if (oenv != nullptr)
        {
            output_env_done(oenv);
        }
        return 0;
    }

    /* Check input for something sensible */
    if ((p_num < 0) || (n_num < 0))
    {
        gmx_fatal(FARGS, "Negative number of ions to add?");
    }

    if (conc > 0 && (p_num > 0 || n_num > 0))
    {
        fprintf(stderr, "WARNING: -conc specified, overriding -nn and -np.\n");
    }

    /* Read atom positions and charges */
    read_tps_conf(ftp2fn(efTPR, NFILE, fnm), &top, &pbcType, &x, nullptr, box, FALSE);
    atoms = top.atoms;

    /* Compute total charge */
    double qtot = 0;
    for (int i = 0; (i < atoms.nr); i++)
    {
        qtot += atoms.atom[i].q;
    }
    iqtot = gmx::roundToInt(qtot);


    if (conc > 0)
    {
        /* Compute number of ions to be added */
        vol   = det(box);
        nsalt = gmx::roundToInt(conc * vol * gmx::c_avogadro / 1e24);
        p_num = std::abs(nsalt * n_q);
        n_num = std::abs(nsalt * p_q);
    }
    if (bNeutral)
    {
        int qdelta = p_num * p_q + n_num * n_q + iqtot;

        /* Check if the system is neutralizable
         * is (qdelta == p_q*p_num + n_q*n_num) solvable for p_num and n_num? */
        int gcd = std::gcd(n_q, p_q);
        if ((qdelta % gcd) != 0)
        {
            gmx_fatal(FARGS,
                      "Can't neutralize this system using -nq %d and"
                      " -pq %d.\n",
                      n_q,
                      p_q);
        }

        while (qdelta != 0)
        {
            while (qdelta < 0)
            {
                p_num++;
                qdelta += p_q;
            }
            while (qdelta > 0)
            {
                n_num++;
                qdelta += n_q;
            }
        }
    }

    char* pptr  = gmx_strdup(p_name);
    char* paptr = aname(p_name);
    char* nptr  = gmx_strdup(n_name);
    char* naptr = aname(n_name);

    if ((p_num == 0) && (n_num == 0))
    {
        fprintf(stderr, "No ions to add, will just copy input configuration.\n");
    }
    else
    {
        char* grpname = nullptr;

        printf("Will try to add %d %s ions and %d %s ions.\n", p_num, p_name, n_num, n_name);
        printf("Select a continuous group of solvent molecules\n");

        std::vector<int> solventGroup;
        {
            int* index = nullptr;
            int  nwa;
            get_index(&atoms, ftp2path_optional(efNDX, NFILE, fnm), 1, &nwa, &index, &grpname);
            solventGroup.assign(index, index + nwa);
            sfree(index);
        }

        for (gmx::Index i = 1; i < gmx::ssize(solventGroup); i++)
        {
            if (solventGroup[i] != solventGroup[i - 1] + 1)
            {
                gmx_fatal(FARGS,
                          "The solvent group %s is not continuous: "
                          "index[%d]=%d, index[%d]=%d",
                          grpname,
                          int(i),
                          solventGroup[i - 1] + 1,
                          int(i + 1),
                          solventGroup[i] + 1);
            }
        }
        nsa = 1;
        while ((nsa < gmx::ssize(solventGroup))
               && (atoms.atom[solventGroup[nsa]].resind == atoms.atom[solventGroup[nsa - 1]].resind))
        {
            nsa++;
        }
        if (solventGroup.size() % nsa != 0)
        {
            gmx_fatal(FARGS,
                      "Your solvent group size (%td) is not a multiple of %d",
                      gmx::ssize(solventGroup),
                      nsa);
        }
        nw = solventGroup.size() / nsa;
        fprintf(stderr, "Number of (%d-atomic) solvent molecules: %d\n", nsa, nw);
        if (p_num + n_num > nw)
        {
            gmx_fatal(FARGS, "Not enough solvent for adding ions");
        }

        if (opt2bSet("-p", NFILE, fnm))
        {
            update_topol(opt2fn("-p", NFILE, fnm), p_num, n_num, p_name, n_name, grpname);
        }

        snew(repl, nw);
        set_pbc(&pbc, pbcType, box);


        if (seed == 0)
        {
            // For now we make do with 32 bits to avoid changing the user input to 64 bit hex
            seed = static_cast<int>(gmx::makeRandomSeed());
        }
        fprintf(stderr, "Using random seed %d.\n", seed);


        std::vector<int> notSolventGroup = invertIndexGroup(atoms.nr, solventGroup);

        std::vector<int> solventMoleculesForReplacement(nw);
        std::iota(std::begin(solventMoleculesForReplacement), std::end(solventMoleculesForReplacement), 0);

        // Randomly shuffle the solvent molecules that shall be replaced by ions
        // then pick molecules from the back of the list as replacement candidates
        gmx::DefaultRandomEngine rng(seed);
        std::shuffle(
                std::begin(solventMoleculesForReplacement), std::end(solventMoleculesForReplacement), rng);

        /* Now loop over the ions that have to be placed */
        while (p_num-- > 0)
        {
            insert_ion(
                    nsa, &solventMoleculesForReplacement, repl, solventGroup, x, &pbc, 1, p_q, p_name, &atoms, rmin, &notSolventGroup);
        }
        while (n_num-- > 0)
        {
            insert_ion(
                    nsa, &solventMoleculesForReplacement, repl, solventGroup, x, &pbc, -1, n_q, n_name, &atoms, rmin, &notSolventGroup);
        }
        fprintf(stderr, "\n");

        if (nw)
        {
            sort_ions(nsa, nw, repl, solventGroup, &atoms, x, &pptr, &nptr, &paptr, &naptr);
        }

        sfree(repl);
        sfree(grpname);
    }

    sfree(atoms.pdbinfo);
    atoms.pdbinfo = nullptr;
    write_sto_conf(ftp2fn(efSTO, NFILE, fnm), *top.name, &atoms, x, nullptr, pbcType, box);

    sfree(pptr);
    sfree(paptr);
    sfree(nptr);
    sfree(naptr);

    sfree(x);
    output_env_done(oenv);
    done_top(&top);

    return 0;
}
