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

#include "solvate.h"

#include <cstdio>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <random>
#include <string>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/gmxlib/conformation_utilities.h"
#include "gromacs/gmxpreprocess/makeexclusiondistances.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/boxutilities.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/random/seed.h"
#include "gromacs/selection/nbsearch.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/atomsbuilder.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

struct gmx_output_env_t;
struct t_symtab;

using gmx::RVec;

/*! \brief Describes a molecule type, and keeps track of the number of these molecules
 *
 *  Used for sorting coordinate file data after solvation
 */
struct MoleculeType
{
    //! molecule name
    std::string name;
    //! number of atoms in the molecule
    int numAtoms = 0;
    //! number of occurences of molecule
    int numMolecules = 0;
};

static void sort_molecule(t_atoms** atoms_solvt, t_atoms** newatoms, std::vector<RVec>* x, std::vector<RVec>* v)
{

    fprintf(stderr, "Sorting configuration\n");
    t_atoms* atoms = *atoms_solvt;

    /* copy each residue from *atoms to a molecule in *molecule */
    std::vector<MoleculeType> molTypes;
    for (int i = 0; i < atoms->nr; i++)
    {
        if ((i == 0) || (atoms->atom[i].resind != atoms->atom[i - 1].resind))
        {
            /* see if this was a molecule type we haven't had yet: */
            auto matchingMolType = std::find_if(
                    molTypes.begin(), molTypes.end(), [atoms, i](const MoleculeType& molecule) {
                        return molecule.name == *atoms->resinfo[atoms->atom[i].resind].name;
                    });
            if (matchingMolType == molTypes.end())
            {
                int numAtomsInMolType = 0;
                while ((i + numAtomsInMolType < atoms->nr)
                       && (atoms->atom[i].resind == atoms->atom[i + numAtomsInMolType].resind))
                {
                    numAtomsInMolType++;
                }
                molTypes.emplace_back(MoleculeType{
                        *atoms->resinfo[atoms->atom[i].resind].name, numAtomsInMolType, 1 });
            }
            else
            {
                matchingMolType->numMolecules++;
            }
        }
    }

    fprintf(stderr,
            "Found %zu%s molecule type%s:\n",
            molTypes.size(),
            molTypes.size() == 1 ? "" : " different",
            molTypes.size() == 1 ? "" : "s");
    for (const auto& molType : molTypes)
    {
        fprintf(stderr, "%7s (%4d atoms): %5d residues\n", molType.name.c_str(), molType.numAtoms, molType.numMolecules);
    }

    /* if we have only 1 moleculetype, we don't have to sort */
    if (molTypes.size() > 1)
    {
        /* now put them there: */
        snew(*newatoms, 1);
        init_t_atoms(*newatoms, atoms->nr, FALSE);
        (*newatoms)->nres = atoms->nres;
        srenew((*newatoms)->resinfo, atoms->nres);
        std::vector<RVec> newx(x->size());
        std::vector<RVec> newv(v->size());

        int residIndex = 0;
        int atomIndex  = 0;
        for (const auto& moleculeType : molTypes)
        {
            int i = 0;
            while (i < atoms->nr)
            {
                int residOfCurrAtom = atoms->atom[i].resind;
                if (moleculeType.name == *atoms->resinfo[residOfCurrAtom].name)
                {
                    /* Copy the residue info */
                    (*newatoms)->resinfo[residIndex] = atoms->resinfo[residOfCurrAtom];
                    // Residue numbering starts from 1, so +1 from the index
                    (*newatoms)->resinfo[residIndex].nr = residIndex + 1;
                    /* Copy the atom info */
                    do
                    {
                        (*newatoms)->atom[atomIndex]        = atoms->atom[i];
                        (*newatoms)->atomname[atomIndex]    = atoms->atomname[i];
                        (*newatoms)->atom[atomIndex].resind = residIndex;
                        copy_rvec((*x)[i], newx[atomIndex]);
                        if (!v->empty())
                        {
                            copy_rvec((*v)[i], newv[atomIndex]);
                        }
                        i++;
                        atomIndex++;
                    } while (i < atoms->nr && atoms->atom[i].resind == residOfCurrAtom);
                    /* Increase the new residue counters */
                    residIndex++;
                }
                else
                {
                    /* Skip this residue */
                    do
                    {
                        i++;
                    } while (i < atoms->nr && atoms->atom[i].resind == residOfCurrAtom);
                }
            }
        }

        /* put them back into the original arrays and throw away temporary arrays */
        done_atom(atoms);
        *atoms_solvt = (*newatoms);
        std::swap(*x, newx);
        std::swap(*v, newv);
    }
}

static void rm_res_pbc(const t_atoms* atoms, std::vector<RVec>* x, matrix box)
{
    int  start, nat;
    rvec xcg;

    start = 0;
    nat   = 0;
    clear_rvec(xcg);
    for (int n = 0; n < atoms->nr; n++)
    {
        if (!is_hydrogen(*(atoms->atomname[n])))
        {
            nat++;
            rvec_inc(xcg, (*x)[n]);
        }
        if ((n + 1 == atoms->nr) || (atoms->atom[n + 1].resind != atoms->atom[n].resind))
        {
            /* if nat==0 we have only hydrogens in the solvent,
               we take last coordinate as cg */
            if (nat == 0)
            {
                nat = 1;
                copy_rvec((*x)[n], xcg);
            }
            svmul(1.0 / nat, xcg, xcg);
            for (int d = 0; d < DIM; d++)
            {
                while (xcg[d] < 0)
                {
                    for (int i = start; i <= n; i++)
                    {
                        (*x)[i][d] += box[d][d];
                    }
                    xcg[d] += box[d][d];
                }
                while (xcg[d] >= box[d][d])
                {
                    for (int i = start; i <= n; i++)
                    {
                        (*x)[i][d] -= box[d][d];
                    }
                    xcg[d] -= box[d][d];
                }
            }
            start = n + 1;
            nat   = 0;
            clear_rvec(xcg);
        }
    }
}

/*! \brief
 * Generates a solvent configuration of desired size by stacking solvent boxes.
 *
 * \param[in,out] atoms      Solvent atoms.
 * \param[in,out] x          Solvent positions.
 * \param[in,out] v          Solvent velocities (`*v` can be NULL).
 * \param[in,out] r          Solvent exclusion radii.
 * \param[in]     box        Initial solvent box.
 * \param[in]     boxTarget  Target box size.
 *
 * The solvent box of desired size is created by stacking the initial box in
 * the smallest k*l*m array that covers the box, and then removing any residue
 * where all atoms are outside the target box (with a small margin).
 * This function does not remove overlap between solvent atoms across the
 * edges.
 *
 * Note that the input configuration should be in the rectangular unit cell and
 * have whole residues.
 */
static void replicateSolventBox(t_atoms*           atoms,
                                std::vector<RVec>* x,
                                std::vector<RVec>* v,
                                std::vector<real>* r,
                                const matrix       box,
                                const matrix       boxTarget)
{
    // Calculate the box multiplication factors.
    ivec n_box;
    int  nmol = 1;
    for (int i = 0; i < DIM; ++i)
    {
        n_box[i] = 1;
        while (n_box[i] * box[i][i] < boxTarget[i][i])
        {
            n_box[i]++;
        }
        nmol *= n_box[i];
    }
    fprintf(stderr, "Will generate new solvent configuration of %dx%dx%d boxes\n", n_box[XX], n_box[YY], n_box[ZZ]);

    // Create arrays for storing the generated system (cannot be done in-place
    // in case the target box is smaller than the original in one dimension,
    // but not in all).
    t_atoms newAtoms;
    init_t_atoms(&newAtoms, 0, FALSE);
    gmx::AtomsBuilder builder(&newAtoms, nullptr);
    builder.reserve(atoms->nr * nmol, atoms->nres * nmol);
    std::vector<RVec> newX(atoms->nr * nmol);
    std::vector<RVec> newV(!v->empty() ? atoms->nr * nmol : 0);
    std::vector<real> newR(atoms->nr * nmol);

    const real maxRadius = *std::max_element(r->begin(), r->end());
    rvec       boxWithMargin;
    for (int i = 0; i < DIM; ++i)
    {
        // The code below is only interested about the box diagonal.
        boxWithMargin[i] = boxTarget[i][i] + 3 * maxRadius;
    }

    for (int ix = 0; ix < n_box[XX]; ++ix)
    {
        rvec delta;
        delta[XX] = ix * box[XX][XX];
        for (int iy = 0; iy < n_box[YY]; ++iy)
        {
            delta[YY] = iy * box[YY][YY];
            for (int iz = 0; iz < n_box[ZZ]; ++iz)
            {
                delta[ZZ]         = iz * box[ZZ][ZZ];
                bool bKeepResidue = false;
                for (int i = 0; i < atoms->nr; ++i)
                {
                    const int newIndex  = builder.currentAtomCount();
                    bool      bKeepAtom = true;
                    for (int m = 0; m < DIM; ++m)
                    {
                        const real newCoord = delta[m] + (*x)[i][m];
                        bKeepAtom           = bKeepAtom && (newCoord < boxWithMargin[m]);
                        newX[newIndex][m]   = newCoord;
                    }
                    bKeepResidue = bKeepResidue || bKeepAtom;
                    if (!v->empty())
                    {
                        copy_rvec((*v)[i], newV[newIndex]);
                    }
                    newR[newIndex] = (*r)[i];
                    builder.addAtom(*atoms, i);
                    if (i == atoms->nr - 1 || atoms->atom[i + 1].resind != atoms->atom[i].resind)
                    {
                        if (bKeepResidue)
                        {
                            builder.finishResidue(atoms->resinfo[atoms->atom[i].resind]);
                        }
                        else
                        {
                            builder.discardCurrentResidue();
                        }
                        // Reset state for the next residue.
                        bKeepResidue = false;
                    }
                }
            }
        }
    }
    sfree(atoms->atom);
    sfree(atoms->atomname);
    sfree(atoms->resinfo);
    atoms->nr       = newAtoms.nr;
    atoms->nres     = newAtoms.nres;
    atoms->atom     = newAtoms.atom;
    atoms->atomname = newAtoms.atomname;
    atoms->resinfo  = newAtoms.resinfo;
    if (atoms->havePdbInfo)
    {
        sfree(atoms->pdbinfo);
        // gmx::AtomsBuilder does not fill pdbinfo, but let's copy nulls just in case.
        atoms->pdbinfo     = newAtoms.pdbinfo;
        atoms->havePdbInfo = newAtoms.havePdbInfo;
    }

    newX.resize(atoms->nr);
    std::swap(*x, newX);
    if (!v->empty())
    {
        newV.resize(atoms->nr);
        std::swap(*v, newV);
    }
    newR.resize(atoms->nr);
    std::swap(*r, newR);

    fprintf(stderr, "Solvent box contains %d atoms in %d residues\n", atoms->nr, atoms->nres);
}

/*! \brief
 * Removes overlap of solvent atoms across the edges.
 *
 * \param[in,out] atoms      Solvent atoms.
 * \param[in,out] x          Solvent positions.
 * \param[in,out] v          Solvent velocities (can be empty).
 * \param[in,out] r          Solvent exclusion radii.
 * \param[in]     pbc        PBC information.
 *
 * Solvent residues that lay on the edges that do not touch the origin are
 * removed if they overlap with other solvent atoms across the PBC.
 * This is done in this way as the assumption is that the input solvent
 * configuration is already equilibrated, and so does not contain any
 * undesirable overlap.  The only overlap that should be removed is caused by
 * cutting the box in half in replicateSolventBox() and leaving a margin of
 * solvent outside those box edges; these atoms can then overlap with those on
 * the opposite box edge in a way that is not part of the pre-equilibrated
 * configuration.
 */
static void removeSolventBoxOverlap(t_atoms*           atoms,
                                    std::vector<RVec>* x,
                                    std::vector<RVec>* v,
                                    std::vector<real>* r,
                                    const t_pbc&       pbc)
{
    gmx::AtomsRemover remover(*atoms);

    // TODO: We could limit the amount of pairs searched significantly,
    // since we are only interested in pairs where the positions are on
    // opposite edges.
    const real                maxRadius = *std::max_element(r->begin(), r->end());
    gmx::AnalysisNeighborhood nb;
    nb.setCutoff(2 * maxRadius);
    gmx::AnalysisNeighborhoodPositions  pos(*x);
    gmx::AnalysisNeighborhoodSearch     search     = nb.initSearch(&pbc, pos);
    gmx::AnalysisNeighborhoodPairSearch pairSearch = search.startPairSearch(pos);
    gmx::AnalysisNeighborhoodPair       pair;
    while (pairSearch.findNextPair(&pair))
    {
        const int i1 = pair.refIndex();
        const int i2 = pair.testIndex();
        if (remover.isMarked(i2))
        {
            pairSearch.skipRemainingPairsForTestPosition();
            continue;
        }
        if (remover.isMarked(i1) || atoms->atom[i1].resind == atoms->atom[i2].resind)
        {
            continue;
        }
        if (pair.distance2() < gmx::square((*r)[i1] + (*r)[i2]))
        {
            rvec dx;
            rvec_sub((*x)[i2], (*x)[i1], dx);
            bool bCandidate1 = false, bCandidate2 = false;
            // To satisfy Clang static analyzer.
            GMX_ASSERT(pbc.ndim_ePBC <= DIM, "Too many periodic dimensions");
            for (int d = 0; d < pbc.ndim_ePBC; ++d)
            {
                // If the distance in some dimension is larger than the
                // cutoff, then it means that the distance has been computed
                // over the PBC.  Mark the position with a larger coordinate
                // for potential removal.
                if (dx[d] > maxRadius)
                {
                    bCandidate2 = true;
                }
                else if (dx[d] < -maxRadius)
                {
                    bCandidate1 = true;
                }
            }
            // Only mark one of the positions for removal if both were
            // candidates.
            if (bCandidate2 && (!bCandidate1 || i2 > i1))
            {
                remover.markResidue(*atoms, i2, true);
                pairSearch.skipRemainingPairsForTestPosition();
            }
            else if (bCandidate1)
            {
                remover.markResidue(*atoms, i1, true);
            }
        }
    }

    remover.removeMarkedElements(x);
    if (!v->empty())
    {
        remover.removeMarkedElements(v);
    }
    remover.removeMarkedElements(r);
    const int originalAtomCount = atoms->nr;
    remover.removeMarkedAtoms(atoms);
    fprintf(stderr, "Removed %d solvent atoms due to solvent-solvent overlap\n", originalAtomCount - atoms->nr);
}

/*! \brief
 * Remove all solvent molecules outside a give radius from the solute.
 *
 * \param[in,out] atoms     Solvent atoms.
 * \param[in,out] x_solvent Solvent positions.
 * \param[in,out] v_solvent Solvent velocities.
 * \param[in,out] r         Atomic exclusion radii.
 * \param[in]     pbc       PBC information.
 * \param[in]     x_solute  Solute positions.
 * \param[in]     rshell    The radius outside the solute molecule.
 */
static void removeSolventOutsideShell(t_atoms*                 atoms,
                                      std::vector<RVec>*       x_solvent,
                                      std::vector<RVec>*       v_solvent,
                                      std::vector<real>*       r,
                                      const t_pbc&             pbc,
                                      const std::vector<RVec>& x_solute,
                                      real                     rshell)
{
    gmx::AtomsRemover         remover(*atoms);
    gmx::AnalysisNeighborhood nb;
    nb.setCutoff(rshell);
    gmx::AnalysisNeighborhoodPositions  posSolute(x_solute);
    gmx::AnalysisNeighborhoodSearch     search = nb.initSearch(&pbc, posSolute);
    gmx::AnalysisNeighborhoodPositions  pos(*x_solvent);
    gmx::AnalysisNeighborhoodPairSearch pairSearch = search.startPairSearch(pos);
    gmx::AnalysisNeighborhoodPair       pair;

    // Remove everything
    remover.markAll();
    // Now put back those within the shell without checking for overlap
    while (pairSearch.findNextPair(&pair))
    {
        remover.markResidue(*atoms, pair.testIndex(), false);
        pairSearch.skipRemainingPairsForTestPosition();
    }
    remover.removeMarkedElements(x_solvent);
    if (!v_solvent->empty())
    {
        remover.removeMarkedElements(v_solvent);
    }
    remover.removeMarkedElements(r);
    const int originalAtomCount = atoms->nr;
    remover.removeMarkedAtoms(atoms);
    fprintf(stderr,
            "Removed %d solvent atoms more than %f nm from solute.\n",
            originalAtomCount - atoms->nr,
            rshell);
}

/*! \brief
 * Removes solvent molecules that overlap with the solute.
 *
 * \param[in,out] atoms    Solvent atoms.
 * \param[in,out] x        Solvent positions.
 * \param[in,out] v        Solvent velocities (can be empty).
 * \param[in,out] r        Solvent exclusion radii.
 * \param[in]     pbc      PBC information.
 * \param[in]     x_solute Solute positions.
 * \param[in]     r_solute Solute exclusion radii.
 */
static void removeSolventOverlappingWithSolute(t_atoms*                 atoms,
                                               std::vector<RVec>*       x,
                                               std::vector<RVec>*       v,
                                               std::vector<real>*       r,
                                               const t_pbc&             pbc,
                                               const std::vector<RVec>& x_solute,
                                               const std::vector<real>& r_solute)
{
    gmx::AtomsRemover remover(*atoms);
    const real        maxRadius1 = *std::max_element(r->begin(), r->end());
    const real        maxRadius2 = *std::max_element(r_solute.begin(), r_solute.end());

    // Now check for overlap.
    gmx::AnalysisNeighborhood     nb;
    gmx::AnalysisNeighborhoodPair pair;
    nb.setCutoff(maxRadius1 + maxRadius2);
    gmx::AnalysisNeighborhoodPositions  posSolute(x_solute);
    gmx::AnalysisNeighborhoodSearch     search = nb.initSearch(&pbc, posSolute);
    gmx::AnalysisNeighborhoodPositions  pos(*x);
    gmx::AnalysisNeighborhoodPairSearch pairSearch = search.startPairSearch(pos);
    while (pairSearch.findNextPair(&pair))
    {
        if (remover.isMarked(pair.testIndex()))
        {
            pairSearch.skipRemainingPairsForTestPosition();
            continue;
        }
        const real r1      = r_solute[pair.refIndex()];
        const real r2      = (*r)[pair.testIndex()];
        const bool bRemove = (pair.distance2() < gmx::square(r1 + r2));
        remover.markResidue(*atoms, pair.testIndex(), bRemove);
    }

    remover.removeMarkedElements(x);
    if (!v->empty())
    {
        remover.removeMarkedElements(v);
    }
    remover.removeMarkedElements(r);
    const int originalAtomCount = atoms->nr;
    remover.removeMarkedAtoms(atoms);
    fprintf(stderr, "Removed %d solvent atoms due to solute-solvent overlap\n", originalAtomCount - atoms->nr);
}

/*! \brief
 * Removes a given number of solvent residues.
 *
 * \param[in,out] atoms           Solvent atoms.
 * \param[in,out] x               Solvent positions.
 * \param[in,out] v               Solvent velocities (can be empty).
 * \param[in]     numberToRemove  Number of residues to remove.
 *
 * This function is called last in the process of creating the solvent box,
 * so it does not operate on the exclusion radii, as no code after this needs
 * them.
 */
static void removeExtraSolventMolecules(t_atoms* atoms, std::vector<RVec>* x, std::vector<RVec>* v, int numberToRemove)
{
    gmx::AtomsRemover               remover(*atoms);
    std::mt19937                    randomNumberGenerator(gmx::makeRandomSeed());
    std::uniform_int_distribution<> randomDistribution(0, atoms->nr - 1);
    while (numberToRemove > 0)
    {
        int atomIndex = randomDistribution(randomNumberGenerator);
        if (!remover.isMarked(atomIndex))
        {
            remover.markResidue(*atoms, atomIndex, true);
            numberToRemove--;
        }
    }
    remover.removeMarkedElements(x);
    if (!v->empty())
    {
        remover.removeMarkedElements(v);
    }
    remover.removeMarkedAtoms(atoms);
}

static void add_solv(const char*        filename,
                     t_atoms*           atoms,
                     t_symtab*          symtab,
                     std::vector<RVec>* x,
                     std::vector<RVec>* v,
                     PbcType            pbcType,
                     matrix             box,
                     AtomProperties*    aps,
                     real               defaultDistance,
                     real               scaleFactor,
                     real               rshell,
                     int                max_sol)
{
    gmx_mtop_t        topSolvent;
    std::vector<RVec> xSolvent, vSolvent;
    matrix            boxSolvent = { { 0 } };
    PbcType           pbcTypeSolvent;

    fprintf(stderr, "Reading solvent configuration\n");
    bool  bTprFileWasRead;
    rvec *temporaryX = nullptr, *temporaryV = nullptr;
    readConfAndTopology(gmx::findLibraryFile(filename).c_str(),
                        &bTprFileWasRead,
                        &topSolvent,
                        &pbcTypeSolvent,
                        &temporaryX,
                        &temporaryV,
                        boxSolvent);
    t_atoms* atomsSolvent;
    snew(atomsSolvent, 1);
    *atomsSolvent = gmx_mtop_global_atoms(topSolvent);
    xSolvent.assign(temporaryX, temporaryX + topSolvent.natoms);
    sfree(temporaryX);
    vSolvent.assign(temporaryV, temporaryV + topSolvent.natoms);
    sfree(temporaryV);
    if (gmx::boxIsZero(boxSolvent))
    {
        gmx_fatal(FARGS,
                  "No box information for solvent in %s, please use a properly formatted file\n",
                  filename);
    }
    if (0 == atomsSolvent->nr)
    {
        gmx_fatal(FARGS, "No solvent in %s, please check your input\n", filename);
    }
    fprintf(stderr, "\n");

    /* initialise distance arrays for solvent configuration */
    fprintf(stderr, "Initialising inter-atomic distances...\n");
    const std::vector<real> exclusionDistances(
            makeExclusionDistances(atoms, aps, defaultDistance, scaleFactor));
    std::vector<real> exclusionDistances_solvt(
            makeExclusionDistances(atomsSolvent, aps, defaultDistance, scaleFactor));

    /* generate a new solvent configuration */
    fprintf(stderr, "Generating solvent configuration\n");
    t_pbc pbc;
    set_pbc(&pbc, pbcType, box);
    if (!gmx::boxesAreEqual(boxSolvent, box))
    {
        if (TRICLINIC(boxSolvent))
        {
            gmx_fatal(FARGS,
                      "Generating from non-rectangular solvent boxes is currently not supported.\n"
                      "You can try to pass the same box for -cp and -cs.");
        }
        /* apply pbc for solvent configuration for whole molecules */
        rm_res_pbc(atomsSolvent, &xSolvent, boxSolvent);
        replicateSolventBox(atomsSolvent, &xSolvent, &vSolvent, &exclusionDistances_solvt, boxSolvent, box);
        if (pbcType != PbcType::No)
        {
            removeSolventBoxOverlap(atomsSolvent, &xSolvent, &vSolvent, &exclusionDistances_solvt, pbc);
        }
    }
    if (atoms->nr > 0)
    {
        if (rshell > 0.0)
        {
            removeSolventOutsideShell(
                    atomsSolvent, &xSolvent, &vSolvent, &exclusionDistances_solvt, pbc, *x, rshell);
        }
        removeSolventOverlappingWithSolute(
                atomsSolvent, &xSolvent, &vSolvent, &exclusionDistances_solvt, pbc, *x, exclusionDistances);
    }

    if (max_sol > 0 && atomsSolvent->nres > max_sol)
    {
        const int numberToRemove = atomsSolvent->nres - max_sol;
        removeExtraSolventMolecules(atomsSolvent, &xSolvent, &vSolvent, numberToRemove);
    }

    /* Sort the solvent mixture, not the protein... */
    t_atoms* newatoms = nullptr;
    // The sort_molecule function does something creative with the
    // t_atoms pointers. We need to make sure we neither leak, nor
    // double-free, so make a shallow pointer that is fine for it to
    // change.
    t_atoms* sortedAtomsSolvent = atomsSolvent;
    sort_molecule(&sortedAtomsSolvent, &newatoms, &xSolvent, &vSolvent);

    // Merge the two configurations.
    x->insert(x->end(), xSolvent.begin(), xSolvent.end());
    if (!v->empty())
    {
        v->insert(v->end(), vSolvent.begin(), vSolvent.end());
    }
    {
        gmx::AtomsBuilder builder(atoms, symtab);
        builder.mergeAtoms(*sortedAtomsSolvent);
    }
    fprintf(stderr,
            "Generated solvent containing %d atoms in %d residues\n",
            atomsSolvent->nr,
            atomsSolvent->nres);

    if (newatoms)
    {
        done_atom(newatoms);
        sfree(newatoms);
    }
    if (atomsSolvent)
    {
        done_atom(atomsSolvent);
        sfree(atomsSolvent);
    }
}

static void update_top(t_atoms*        atoms,
                       int             firstSolventResidueIndex,
                       matrix          box,
                       int             NFILE,
                       t_filenm        fnm[],
                       AtomProperties* aps)
{
    FILE *      fpin, *fpout;
    char        buf[STRLEN * 2], buf2[STRLEN], *temp;
    const char* topinout;
    bool        bSystem;
    int         i;
    double      mtot;
    real        vol, mm;

    int nsol = atoms->nres - firstSolventResidueIndex;

    mtot = 0;
    for (i = 0; (i < atoms->nr); i++)
    {
        aps->setAtomProperty(epropMass,
                             std::string(*atoms->resinfo[atoms->atom[i].resind].name),
                             std::string(*atoms->atomname[i]),
                             &mm);
        mtot += mm;
    }

    vol = det(box);

    fprintf(stderr, "Volume                 :  %10g (nm^3)\n", vol);
    fprintf(stderr, "Density                :  %10g (g/l)\n", (mtot * 1e24) / (gmx::c_avogadro * vol));
    fprintf(stderr, "Number of solvent molecules:  %5d   \n\n", nsol);

    /* open topology file and append sol molecules */
    topinout = ftp2fn(efTOP, NFILE, fnm);
    if (ftp2bSet(efTOP, NFILE, fnm))
    {
        char temporary_filename[STRLEN];
        strncpy(temporary_filename, "temp.topXXXXXX", STRLEN);

        fprintf(stderr, "Processing topology\n");
        fpin    = gmx_ffopen(topinout, "r");
        fpout   = gmx_fopen_temporary(temporary_filename);
        bSystem = false;
        while (fgets(buf, STRLEN, fpin))
        {
            strcpy(buf2, buf);
            if ((temp = strchr(buf2, '\n')) != nullptr)
            {
                temp[0] = '\0';
            }
            ltrim(buf2);
            if (buf2[0] == '[')
            {
                buf2[0] = ' ';
                if ((temp = strchr(buf2, '\n')) != nullptr)
                {
                    temp[0] = '\0';
                }
                rtrim(buf2);
                if (buf2[strlen(buf2) - 1] == ']')
                {
                    buf2[strlen(buf2) - 1] = '\0';
                    ltrim(buf2);
                    rtrim(buf2);
                    bSystem = (gmx_strcasecmp(buf2, "system") == 0);
                }
            }
            else if (bSystem && nsol && (buf[0] != ';'))
            {
                /* if sol present, append "in water" to system name */
                rtrim(buf2);
                if (buf2[0] && (!strstr(buf2, " water")))
                {
                    sprintf(buf, "%s in water\n", buf2);
                    bSystem = false;
                }
            }
            fprintf(fpout, "%s", buf);
        }
        gmx_ffclose(fpin);

        // Add new solvent molecules to the topology
        if (nsol > 0)
        {
            std::string currRes  = *atoms->resinfo[firstSolventResidueIndex].name;
            int         resCount = 0;

            // Iterate through solvent molecules and increment a count until new resname found
            for (int i = firstSolventResidueIndex; i < atoms->nres; i++)
            {
                if ((currRes == *atoms->resinfo[i].name))
                {
                    resCount += 1;
                }
                else
                {
                    // Change topology and restart count
                    fprintf(stdout,
                            "Adding line for %d solvent molecules with resname (%s) to "
                            "topology file (%s)\n",
                            resCount,
                            currRes.c_str(),
                            topinout);
                    fprintf(fpout, "%-15s %5d\n", currRes.c_str(), resCount);
                    currRes  = *atoms->resinfo[i].name;
                    resCount = 1;
                }
            }
            // One more print needed for last residue type
            fprintf(stdout,
                    "Adding line for %d solvent molecules with resname (%s) to "
                    "topology file (%s)\n",
                    resCount,
                    currRes.c_str(),
                    topinout);
            fprintf(fpout, "%-15s %5d\n", currRes.c_str(), resCount);
        }
        gmx_ffclose(fpout);
        make_backup(topinout);
        gmx_file_rename(temporary_filename, topinout);
    }
}

int gmx_solvate(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] can do one of 2 things:[PAR]",

        "1) Generate a box of solvent. Specify [TT]-cs[tt] and [TT]-box[tt].",
        "Or specify [TT]-cs[tt] and [TT]-cp[tt] with a structure file with",
        "a box, but without atoms.[PAR]",

        "2) Solvate a solute configuration, e.g. a protein, in a bath of solvent ",
        "molecules. Specify [TT]-cp[tt] (solute) and [TT]-cs[tt] (solvent). ",
        "The box specified in the solute coordinate file ([TT]-cp[tt]) is used,",
        "unless [TT]-box[tt] is set.",
        "If you want the solute to be centered in the box,",
        "the program [gmx-editconf] has sophisticated options",
        "to change the box dimensions and center the solute.",
        "Solvent molecules are removed from the box where the ",
        "distance between any atom of the solute molecule(s) and any atom of ",
        "the solvent molecule is less than the sum of the scaled van der Waals",
        "radii of both atoms. A database ([TT]vdwradii.dat[tt]) of van der",
        "Waals radii is read by the program, and the resulting radii scaled",
        "by [TT]-scale[tt]. If radii are not found in the database, those",
        "atoms are assigned the (pre-scaled) distance [TT]-radius[tt].",
        "Note that the usefulness of those radii depends on the atom names,",
        "and thus varies widely with force field.",
        "",
        "The default solvent is Simple Point Charge water (SPC), with coordinates ",
        "from [TT]$GMXLIB/spc216.gro[tt]. These coordinates can also be used",
        "for other 3-site water models, since a short equibilibration will remove",
        "the small differences between the models.",
        "Other solvents are also supported, as well as mixed solvents. The",
        "only restriction to solvent types is that a solvent molecule consists",
        "of exactly one residue. The residue information in the coordinate",
        "files is used, and should therefore be more or less consistent.",
        "In practice this means that two subsequent solvent molecules in the ",
        "solvent coordinate file should have different residue number.",
        "The box of solute is built by stacking the coordinates read from",
        "the coordinate file. This means that these coordinates should be ",
        "equlibrated in periodic boundary conditions to ensure a good",
        "alignment of molecules on the stacking interfaces.",
        "The [TT]-maxsol[tt] option simply adds only the first [TT]-maxsol[tt]",
        "solvent molecules and leaves out the rest that would have fitted",
        "into the box. This can create a void that can cause problems later.",
        "Choose your volume wisely.[PAR]",

        "Setting [TT]-shell[tt] larger than zero will place a layer of water of",
        "the specified thickness (nm) around the solute. Hint: it is a good",
        "idea to put the protein in the center of a box first (using [gmx-editconf]).",
        "[PAR]",

        "Finally, [THISMODULE] will optionally remove lines from your topology file in ",
        "which a number of solvent molecules is already added, and adds a ",
        "line with the total number of solvent molecules in your coordinate file."
    };

    const char* bugs[] = {
        "Molecules must be whole in the initial configurations.",
    };

    /* parameter data */
    gmx_bool    bProt, bBox;
    const char *conf_prot, *confout;

    t_filenm fnm[] = {
        { efSTX, "-cp", "protein", ffOPTRD },
        { efSTX, "-cs", "spc216", ffLIBRD },
        { efSTO, nullptr, nullptr, ffWRITE },
        { efTOP, nullptr, nullptr, ffOPTRW },
    };
#define NFILE asize(fnm)

    real              defaultDistance = 0.105, r_shell = 0, scaleFactor = 0.57;
    rvec              new_box                  = { 0.0, 0.0, 0.0 };
    gmx_bool          bReadV                   = FALSE;
    int               max_sol                  = 0;
    int               firstSolventResidueIndex = 0;
    gmx_output_env_t* oenv;
    t_pargs           pa[] = {
        { "-box", FALSE, etRVEC, { new_box }, "Box size (in nm)" },
        { "-radius", FALSE, etREAL, { &defaultDistance }, "Default van der Waals distance" },
        { "-scale",
          FALSE,
          etREAL,
          { &scaleFactor },
          "Scale factor to multiply Van der Waals radii from the database in "
          "share/gromacs/top/vdwradii.dat. The default value of 0.57 yields density close to 1000 "
          "g/l for proteins in water." },
        { "-shell", FALSE, etREAL, { &r_shell }, "Thickness of optional water layer around solute" },
        { "-maxsol",
          FALSE,
          etINT,
          { &max_sol },
          "Maximum number of solvent molecules to add if they fit in the box. If zero (default) "
          "this is ignored" },
        { "-vel", FALSE, etBOOL, { &bReadV }, "Keep velocities from input solute and solvent" },
    };

    if (!parse_common_args(
                &argc, argv, 0, NFILE, fnm, asize(pa), pa, asize(desc), desc, asize(bugs), bugs, &oenv))
    {
        return 0;
    }

    const char* solventFileName = opt2fn("-cs", NFILE, fnm);
    bProt                       = opt2bSet("-cp", NFILE, fnm);
    bBox                        = opt2parg_bSet("-box", asize(pa), pa);

    /* check input */
    if (!bProt && !bBox)
    {
        gmx_fatal(FARGS,
                  "When no solute (-cp) is specified, "
                  "a box size (-box) must be specified");
    }

    AtomProperties aps;

    /* solute configuration data */
    gmx_mtop_t        top;
    std::vector<RVec> x, v;
    matrix            box     = { { 0 } };
    PbcType           pbcType = PbcType::Unset;
    t_atoms*          atoms;
    snew(atoms, 1);
    if (bProt)
    {
        /* Generate a solute configuration */
        conf_prot = opt2fn("-cp", NFILE, fnm);
        fprintf(stderr, "Reading solute configuration%s\n", bReadV ? " and velocities" : "");
        bool  bTprFileWasRead;
        rvec *temporaryX = nullptr, *temporaryV = nullptr;
        readConfAndTopology(
                conf_prot, &bTprFileWasRead, &top, &pbcType, &temporaryX, bReadV ? &temporaryV : nullptr, box);
        *atoms = gmx_mtop_global_atoms(top);
        x.assign(temporaryX, temporaryX + top.natoms);
        sfree(temporaryX);
        if (temporaryV)
        {
            v.assign(temporaryV, temporaryV + top.natoms);
            sfree(temporaryV);
        }
        else if (bReadV)
        {
            fprintf(stderr, "Note: no velocities found\n");
        }
        if (atoms->nr == 0)
        {
            fprintf(stderr, "Note: no atoms in %s\n", conf_prot);
            bProt = FALSE;
        }
        else
        {
            firstSolventResidueIndex = atoms->nres;
        }
    }
    PbcType pbcTypeForOutput = pbcType;
    if (bBox)
    {
        pbcTypeForOutput = PbcType::Xyz;
        clear_mat(box);
        box[XX][XX] = new_box[XX];
        box[YY][YY] = new_box[YY];
        box[ZZ][ZZ] = new_box[ZZ];
    }
    if (det(box) == 0)
    {
        gmx_fatal(FARGS,
                  "Undefined solute box.\nCreate one with gmx editconf "
                  "or give explicit -box command line option");
    }

    add_solv(solventFileName, atoms, &top.symtab, &x, &v, pbcTypeForOutput, box, &aps, defaultDistance, scaleFactor, r_shell, max_sol);

    /* write new configuration 1 to file confout */
    confout = ftp2fn(efSTO, NFILE, fnm);
    fprintf(stderr, "Writing generated configuration to %s\n", confout);
    const char* outputTitle = (bProt ? *top.name : "Generated by gmx solvate");
    write_sto_conf(confout,
                   outputTitle,
                   atoms,
                   as_rvec_array(x.data()),
                   !v.empty() ? as_rvec_array(v.data()) : nullptr,
                   pbcTypeForOutput,
                   box);

    /* print size of generated configuration */
    fprintf(stderr, "\nOutput configuration contains %d atoms in %d residues\n", atoms->nr, atoms->nres);
    update_top(atoms, firstSolventResidueIndex, box, NFILE, fnm, &aps);

    done_atom(atoms);
    sfree(atoms);
    output_env_done(oenv);

    return 0;
}
