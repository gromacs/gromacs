/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#include "solvate.h"

#include <string.h>

#include <algorithm>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/gmxlib/conformation-utilities.h"
#include "gromacs/gmxpreprocess/read-conformation.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/nbsearch.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/atomsbuilder.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

typedef struct {
    char *name;
    int   natoms;
    int   nmol;
    int   i, i0;
    int   res0;
} t_moltypes;

static void sort_molecule(t_atoms **atoms_solvt, rvec *x, rvec *v, real *r)
{
    int         atnr, i, j, moltp = 0, nrmoltypes, resi_o, resi_n, resnr;
    t_moltypes *moltypes;
    t_atoms    *atoms, *newatoms;
    rvec       *newx, *newv = NULL;
    real       *newr;

    fprintf(stderr, "Sorting configuration\n");

    atoms = *atoms_solvt;

    /* copy each residue from *atoms to a molecule in *molecule */
    moltypes   = NULL;
    nrmoltypes = 0;
    for (i = 0; i < atoms->nr; i++)
    {
        if ( (i == 0) || (atoms->atom[i].resind != atoms->atom[i-1].resind) )
        {
            /* see if this was a molecule type we haven't had yet: */
            moltp = -1;
            for (j = 0; (j < nrmoltypes) && (moltp == -1); j++)
            {
                /* cppcheck-suppress nullPointer
                 * moltypes is guaranteed to be allocated because otherwise
                 * nrmoltypes is 0. */
                if (strcmp(*(atoms->resinfo[atoms->atom[i].resind].name), moltypes[j].name) == 0)
                {
                    moltp = j;
                }
            }
            if (moltp == -1)
            {
                moltp = nrmoltypes;
                nrmoltypes++;
                srenew(moltypes, nrmoltypes);
                moltypes[moltp].name = *(atoms->resinfo[atoms->atom[i].resind].name);
                atnr                 = 0;
                while ((i+atnr < atoms->nr) &&
                       (atoms->atom[i].resind == atoms->atom[i+atnr].resind))
                {
                    atnr++;
                }
                moltypes[moltp].natoms = atnr;
                moltypes[moltp].nmol   = 0;
            }
            moltypes[moltp].nmol++;
        }
    }

    fprintf(stderr, "Found %d%s molecule type%s:\n",
            nrmoltypes, nrmoltypes == 1 ? "" : " different", nrmoltypes == 1 ? "" : "s");
    for (j = 0; j < nrmoltypes; j++)
    {
        if (j == 0)
        {
            moltypes[j].res0 = 0;
        }
        else
        {
            moltypes[j].res0 = moltypes[j-1].res0+moltypes[j-1].nmol;
        }
        fprintf(stderr, "%7s (%4d atoms): %5d residues\n",
                moltypes[j].name, moltypes[j].natoms, moltypes[j].nmol);
    }

    /* if we have only 1 moleculetype, we don't have to sort */
    if (nrmoltypes > 1)
    {
        /* find out which molecules should go where: */
        moltypes[0].i = moltypes[0].i0 = 0;
        for (j = 1; j < nrmoltypes; j++)
        {
            moltypes[j].i      =
                moltypes[j].i0 =
                    moltypes[j-1].i0+moltypes[j-1].natoms*moltypes[j-1].nmol;
        }

        /* now put them there: */
        snew(newatoms, 1);
        init_t_atoms(newatoms, atoms->nr, FALSE);
        newatoms->nres = atoms->nres;
        snew(newatoms->resinfo, atoms->nres);
        snew(newx, atoms->nr);
        if (v)
        {
            snew(newv, atoms->nr);
        }
        snew(newr, atoms->nr);

        resi_n = 0;
        resnr  = 1;
        j      = 0;
        for (moltp = 0; moltp < nrmoltypes; moltp++)
        {
            i = 0;
            while (i < atoms->nr)
            {
                resi_o = atoms->atom[i].resind;
                if (strcmp(*atoms->resinfo[resi_o].name, moltypes[moltp].name) == 0)
                {
                    /* Copy the residue info */
                    newatoms->resinfo[resi_n]    = atoms->resinfo[resi_o];
                    newatoms->resinfo[resi_n].nr = resnr;
                    /* Copy the atom info */
                    do
                    {
                        newatoms->atom[j]        = atoms->atom[i];
                        newatoms->atomname[j]    = atoms->atomname[i];
                        newatoms->atom[j].resind = resi_n;
                        copy_rvec(x[i], newx[j]);
                        if (v != NULL)
                        {
                            copy_rvec(v[i], newv[j]);
                        }
                        newr[j] = r[i];
                        i++;
                        j++;
                    }
                    while (i < atoms->nr && atoms->atom[i].resind == resi_o);
                    /* Increase the new residue counters */
                    resi_n++;
                    resnr++;
                }
                else
                {
                    /* Skip this residue */
                    do
                    {
                        i++;
                    }
                    while (i < atoms->nr && atoms->atom[i].resind == resi_o);
                }
            }
        }

        /* put them back into the original arrays and throw away temporary arrays */
        sfree(atoms->atomname);
        sfree(atoms->resinfo);
        sfree(atoms->atom);
        sfree(atoms);
        *atoms_solvt = newatoms;
        for (i = 0; i < (*atoms_solvt)->nr; i++)
        {
            copy_rvec(newx[i], x[i]);
            if (v)
            {
                copy_rvec(newv[i], v[i]);
            }
            r[i] = newr[i];
        }
        sfree(newx);
        if (v)
        {
            sfree(newv);
        }
        sfree(newr);
    }
    sfree(moltypes);
}

static void rm_res_pbc(t_atoms *atoms, rvec *x, matrix box)
{
    int  i, start, n, d, nat;
    rvec xcg;

    start = 0;
    nat   = 0;
    clear_rvec(xcg);
    for (n = 0; n < atoms->nr; n++)
    {
        if (!is_hydrogen(*(atoms->atomname[n])))
        {
            nat++;
            rvec_inc(xcg, x[n]);
        }
        if ( (n+1 == atoms->nr) ||
             (atoms->atom[n+1].resind != atoms->atom[n].resind) )
        {
            /* if nat==0 we have only hydrogens in the solvent,
               we take last coordinate as cg */
            if (nat == 0)
            {
                nat = 1;
                copy_rvec(x[n], xcg);
            }
            svmul(1.0/nat, xcg, xcg);
            for (d = 0; d < DIM; d++)
            {
                while (xcg[d] < 0)
                {
                    for (i = start; i <= n; i++)
                    {
                        x[i][d] += box[d][d];
                    }
                    xcg[d] += box[d][d];
                }
                while (xcg[d] >= box[d][d])
                {
                    for (i = start; i <= n; i++)
                    {
                        x[i][d] -= box[d][d];
                    }
                    xcg[d] -= box[d][d];
                }
            }
            start = n+1;
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
static void replicateSolventBox(t_atoms *atoms, rvec **x, rvec **v, real **r,
                                const matrix box, const matrix boxTarget)
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
    fprintf(stderr, "Will generate new solvent configuration of %dx%dx%d boxes\n",
            n_box[XX], n_box[YY], n_box[ZZ]);

    // Create arrays for storing the generated system (cannot be done in-place
    // in case the target box is smaller than the original in one dimension,
    // but not in all).
    t_atoms           newAtoms;
    init_t_atoms(&newAtoms, 0, FALSE);
    gmx::AtomsBuilder builder(&newAtoms);
    builder.reserve(atoms->nr * nmol, atoms->nres * nmol);
    rvec             *newX;
    rvec             *newV = NULL;
    real             *newR;
    snew(newX,              atoms->nr * nmol);
    snew(newR,              atoms->nr * nmol);
    if (*v)
    {
        snew(newV,          atoms->nr * nmol);
    }

    const real maxRadius = *std::max_element(*r, *r + atoms->nr);
    rvec       boxWithMargin;
    for (int i = 0; i < DIM; ++i)
    {
        // The code below is only interested about the box diagonal.
        boxWithMargin[i] = boxTarget[i][i] + 3*maxRadius;
    }

    for (int ix = 0; ix < n_box[XX]; ++ix)
    {
        rvec delta;
        delta[XX] = ix*box[XX][XX];
        for (int iy = 0; iy < n_box[YY]; ++iy)
        {
            delta[YY] = iy*box[YY][YY];
            for (int iz = 0; iz < n_box[ZZ]; ++iz)
            {
                delta[ZZ] = iz*box[ZZ][ZZ];
                bool bKeepResidue     = false;
                for (int i = 0; i < atoms->nr; ++i)
                {
                    const int newIndex  = builder.currentAtomCount();
                    bool      bKeepAtom = true;
                    for (int m = 0; m < DIM; ++m)
                    {
                        const real newCoord = delta[m] + (*x)[i][m];
                        bKeepAtom         = bKeepAtom && (newCoord < boxWithMargin[m]);
                        newX[newIndex][m] = newCoord;
                    }
                    bKeepResidue = bKeepResidue || bKeepAtom;
                    if (newV)
                    {
                        copy_rvec((*v)[i], newV[newIndex]);
                    }
                    newR[newIndex] = (*r)[i];
                    builder.addAtom(*atoms, i);
                    if (i == atoms->nr - 1
                        || atoms->atom[i+1].resind != atoms->atom[i].resind)
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
                        bKeepResidue     = false;
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
    sfree(*x);
    sfree(*v);
    sfree(*r);
    *x = newX;
    *v = newV;
    *r = newR;
    fprintf(stderr, "Solvent box contains %d atoms in %d residues\n",
            atoms->nr, atoms->nres);
}

/*! \brief
 * Removes overlap of solvent atoms across the edges.
 *
 * \param[in,out] atoms      Solvent atoms.
 * \param[in,out] x          Solvent positions.
 * \param[in,out] v          Solvent velocities (can be NULL).
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
static void removeSolventBoxOverlap(t_atoms *atoms, rvec *x, rvec *v, real *r,
                                    const t_pbc &pbc)
{
    gmx::AtomsRemover remover(*atoms);

    // TODO: We could limit the amount of pairs searched significantly,
    // since we are only interested in pairs where the positions are on
    // opposite edges.
    const real maxRadius = *std::max_element(r, r + atoms->nr);
    gmx::AnalysisNeighborhood           nb;
    nb.setCutoff(2*maxRadius);
    gmx::AnalysisNeighborhoodPositions  pos(x, atoms->nr);
    gmx::AnalysisNeighborhoodSearch     search     = nb.initSearch(&pbc, pos);
    gmx::AnalysisNeighborhoodPairSearch pairSearch = search.startPairSearch(pos);
    gmx::AnalysisNeighborhoodPair       pair;
    while (pairSearch.findNextPair(&pair))
    {
        const int  i1 = pair.refIndex();
        const int  i2 = pair.testIndex();
        if (remover.isMarked(i2))
        {
            pairSearch.skipRemainingPairsForTestPosition();
            continue;
        }
        if (remover.isMarked(i1) || atoms->atom[i1].resind == atoms->atom[i2].resind)
        {
            continue;
        }
        if (pair.distance2() < sqr(r[i1] + r[i2]))
        {
            rvec dx;
            rvec_sub(x[i2], x[i1], dx);
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

    remover.removeMarkedVectors(x);
    if (v != NULL)
    {
        remover.removeMarkedVectors(v);
    }
    remover.removeMarkedValues(r);
    const int originalAtomCount = atoms->nr;
    remover.removeMarkedAtoms(atoms);
    fprintf(stderr, "Removed %d solvent atoms due to solvent-solvent overlap\n",
            originalAtomCount - atoms->nr);
}

/*! \brief
 * Removes solvent molecules that overlap with the solute, and optionally also
 * those that are outside a given shell radius from the solute.
 *
 * \param[in,out] atoms      Solvent atoms.
 * \param[in,out] x          Solvent positions.
 * \param[in,out] v          Solvent velocities (can be NULL).
 * \param[in,out] r          Solvent exclusion radii.
 * \param[in]     pbc        PBC information.
 * \param[in]     soluteAtomCount Number of solute atoms.
 * \param[in]     x_solute   Solute positions.
 * \param[in]     r_solute   Solute exclusion radii.
 * \param[in]     rshell     If >0, only keep solvent atoms within a shell of
 *     this size from the solute.
 */
static void removeSoluteOverlap(t_atoms *atoms, rvec *x, rvec *v, real *r,
                                const t_pbc &pbc, int soluteAtomCount,
                                const rvec *x_solute, const real *r_solute,
                                real rshell)
{
    const real                          maxRadius1
        = *std::max_element(r, r + atoms->nr);
    const real                          maxRadius2
        = *std::max_element(r_solute, r_solute + soluteAtomCount);

    gmx::AtomsRemover                   remover(*atoms);
    // If rshell is >0, the neighborhood search looks at all pairs
    // within rshell, and unmarks those that are within the cutoff.
    // This line marks everything, so that solvent outside rshell remains
    // marked after the loop.
    // Without rshell, the neighborhood search only marks the overlapping
    // solvent atoms, and all others are left alone.
    if (rshell > 0.0)
    {
        remover.markAll();
    }

    gmx::AnalysisNeighborhood           nb;
    nb.setCutoff(std::max(maxRadius1 + maxRadius2, rshell));
    gmx::AnalysisNeighborhoodPositions  posSolute(x_solute, soluteAtomCount);
    gmx::AnalysisNeighborhoodSearch     search     = nb.initSearch(&pbc, posSolute);
    gmx::AnalysisNeighborhoodPositions  pos(x, atoms->nr);
    gmx::AnalysisNeighborhoodPairSearch pairSearch = search.startPairSearch(pos);
    gmx::AnalysisNeighborhoodPair       pair;
    while (pairSearch.findNextPair(&pair))
    {
        if (remover.isMarked(pair.testIndex()))
        {
            pairSearch.skipRemainingPairsForTestPosition();
            continue;
        }
        const real r1      = r_solute[pair.refIndex()];
        const real r2      = r[pair.testIndex()];
        const bool bRemove = (pair.distance2() < sqr(r1 + r2));
        remover.markResidue(*atoms, pair.testIndex(), bRemove);
    }

    remover.removeMarkedVectors(x);
    if (v != NULL)
    {
        remover.removeMarkedVectors(v);
    }
    remover.removeMarkedValues(r);
    const int originalAtomCount = atoms->nr;
    remover.removeMarkedAtoms(atoms);
    fprintf(stderr, "Removed %d solvent atoms due to solute-solvent overlap\n",
            originalAtomCount - atoms->nr);
}

/*! \brief
 * Removes a given number of solvent residues.
 *
 * \param[in,out] atoms           Solvent atoms.
 * \param[in,out] x               Solvent positions.
 * \param[in,out] v               Solvent velocities (can be NULL).
 * \param[in]     numberToRemove  Number of residues to remove.
 *
 * This function is called last in the process of creating the solvent box,
 * so it does not operate on the exclusion radii, as no code after this needs
 * them.
 */
static void removeExtraSolventMolecules(t_atoms *atoms, rvec *x, rvec *v,
                                        int numberToRemove)
{
    gmx::AtomsRemover remover(*atoms);
    // TODO: It might be nicer to remove a random set of residues, but
    // in practice this should give a roughly uniform spatial distribution.
    const int stride = atoms->nr / numberToRemove;
    for (int i = 0; i < numberToRemove; ++i)
    {
        int atomIndex = (i+1)*stride - 1;
        while (remover.isMarked(atomIndex))
        {
            ++atomIndex;
            if (atomIndex == atoms->nr)
            {
                atomIndex = 0;
            }
        }
        remover.markResidue(*atoms, atomIndex, true);
    }
    remover.removeMarkedVectors(x);
    if (v != NULL)
    {
        remover.removeMarkedVectors(v);
    }
    remover.removeMarkedAtoms(atoms);
}

static void add_solv(const char *fn, t_atoms *atoms, rvec **x, rvec **v,
                     int ePBC, matrix box, gmx_atomprop_t aps,
                     real defaultDistance, real scaleFactor,
                     real rshell, int max_sol)
{
    t_atoms *atoms_solvt;
    rvec    *x_solvt, *v_solvt = NULL;
    int      ePBC_solvt;
    matrix   box_solvt;

    char    *filename = gmxlibfn(fn);
    snew(atoms_solvt, 1);
    char    *title_solvt
        = readConformation(filename, atoms_solvt, &x_solvt, &v_solvt,
                           &ePBC_solvt, box_solvt, "solvent");
    sfree(title_solvt);
    if (0 == atoms_solvt->nr)
    {
        gmx_fatal(FARGS, "No solvent in %s, please check your input\n", filename);
    }
    sfree(filename);
    fprintf(stderr, "\n");

    /* apply pbc for solvent configuration for whole molecules */
    rm_res_pbc(atoms_solvt, x_solvt, box_solvt);

    /* initialise distance arrays for solvent configuration */
    fprintf(stderr, "Initialising inter-atomic distances...\n");
    real *exclusionDistances
        = makeExclusionDistances(atoms, aps, defaultDistance, scaleFactor);
    real *exclusionDistances_solvt
        = makeExclusionDistances(atoms_solvt, aps, defaultDistance, scaleFactor);

    /* generate a new solvent configuration */
    fprintf(stderr, "Generating solvent configuration\n");
    t_pbc pbc;
    set_pbc(&pbc, ePBC, box);
    replicateSolventBox(atoms_solvt, &x_solvt, &v_solvt, &exclusionDistances_solvt,
                        box_solvt, box);
    if (ePBC != epbcNONE)
    {
        removeSolventBoxOverlap(atoms_solvt, x_solvt, v_solvt, exclusionDistances_solvt, pbc);
    }
    if (atoms->nr > 0)
    {
        removeSoluteOverlap(atoms_solvt, x_solvt, v_solvt, exclusionDistances_solvt, pbc,
                            atoms->nr, *x, exclusionDistances, rshell);
    }

    if (max_sol > 0 && atoms_solvt->nres > max_sol)
    {
        const int numberToRemove = atoms_solvt->nres - max_sol;
        removeExtraSolventMolecules(atoms_solvt, x_solvt, v_solvt, numberToRemove);
    }

    /* Sort the solvent mixture, not the protein... */
    sort_molecule(&atoms_solvt, x_solvt, v_solvt, exclusionDistances_solvt);

    // Merge the two configurations.
    srenew(*x, atoms->nr + atoms_solvt->nr);
    if (v != NULL)
    {
        srenew(*v, atoms->nr + atoms_solvt->nr);
    }
    for (int i = 0; i < atoms_solvt->nr; ++i)
    {
        const int index = atoms->nr + i;
        copy_rvec(x_solvt[i], (*x)[index]);
        if (v != NULL)
        {
            copy_rvec(v_solvt[i], (*v)[index]);
        }
    }
    {
        gmx::AtomsBuilder builder(atoms);
        builder.mergeAtoms(*atoms_solvt);
    }
    fprintf(stderr, "Generated solvent containing %d atoms in %d residues\n",
            atoms_solvt->nr, atoms_solvt->nres);

    sfree(x_solvt);
    sfree(v_solvt);
    sfree(exclusionDistances);
    sfree(exclusionDistances_solvt);
    done_atom(atoms_solvt);
    sfree(atoms_solvt);
}

static void update_top(t_atoms *atoms, matrix box, int NFILE, t_filenm fnm[],
                       gmx_atomprop_t aps)
{
    FILE       *fpin, *fpout;
    char        buf[STRLEN], buf2[STRLEN], *temp;
    const char *topinout;
    int         line;
    gmx_bool    bSystem, bMolecules, bSkip;
    int         i, nsol = 0;
    double      mtot;
    real        vol, mm;

    for (i = 0; (i < atoms->nres); i++)
    {
        /* calculate number of SOLvent molecules */
        if ( (strcmp(*atoms->resinfo[i].name, "SOL") == 0) ||
             (strcmp(*atoms->resinfo[i].name, "WAT") == 0) ||
             (strcmp(*atoms->resinfo[i].name, "HOH") == 0) )
        {
            nsol++;
        }
    }
    mtot = 0;
    for (i = 0; (i < atoms->nr); i++)
    {
        gmx_atomprop_query(aps, epropMass,
                           *atoms->resinfo[atoms->atom[i].resind].name,
                           *atoms->atomname[i], &mm);
        mtot += mm;
    }

    vol = det(box);

    fprintf(stderr, "Volume                 :  %10g (nm^3)\n", vol);
    fprintf(stderr, "Density                :  %10g (g/l)\n",
            (mtot*1e24)/(AVOGADRO*vol));
    fprintf(stderr, "Number of SOL molecules:  %5d   \n\n", nsol);

    /* open topology file and append sol molecules */
    topinout  = ftp2fn(efTOP, NFILE, fnm);
    if (ftp2bSet(efTOP, NFILE, fnm) )
    {
        char temporary_filename[STRLEN];
        strncpy(temporary_filename, "temp.topXXXXXX", STRLEN);

        fprintf(stderr, "Processing topology\n");
        fpin    = gmx_ffopen(topinout, "r");
        gmx_tmpnam(temporary_filename);
        fpout   = gmx_ffopen(temporary_filename, "w");
        line    = 0;
        bSystem = bMolecules = FALSE;
        while (fgets(buf, STRLEN, fpin))
        {
            bSkip = FALSE;
            line++;
            strcpy(buf2, buf);
            if ((temp = strchr(buf2, '\n')) != NULL)
            {
                temp[0] = '\0';
            }
            ltrim(buf2);
            if (buf2[0] == '[')
            {
                buf2[0] = ' ';
                if ((temp = strchr(buf2, '\n')) != NULL)
                {
                    temp[0] = '\0';
                }
                rtrim(buf2);
                if (buf2[strlen(buf2)-1] == ']')
                {
                    buf2[strlen(buf2)-1] = '\0';
                    ltrim(buf2);
                    rtrim(buf2);
                    bSystem    = (gmx_strcasecmp(buf2, "system") == 0);
                    bMolecules = (gmx_strcasecmp(buf2, "molecules") == 0);
                }
            }
            else if (bSystem && nsol && (buf[0] != ';') )
            {
                /* if sol present, append "in water" to system name */
                rtrim(buf2);
                if (buf2[0] && (!strstr(buf2, " water")) )
                {
                    sprintf(buf, "%s in water\n", buf2);
                    bSystem = FALSE;
                }
            }
            else if (bMolecules)
            {
                /* check if this is a line with solvent molecules */
                sscanf(buf, "%4095s", buf2);
                if (strcmp(buf2, "SOL") == 0)
                {
                    sscanf(buf, "%*4095s %20d", &i);
                    nsol -= i;
                    if (nsol < 0)
                    {
                        bSkip = TRUE;
                        nsol += i;
                    }
                }
            }
            if (bSkip)
            {
                if ((temp = strchr(buf, '\n')) != NULL)
                {
                    temp[0] = '\0';
                }
                fprintf(stdout, "Removing line #%d '%s' from topology file (%s)\n",
                        line, buf, topinout);
            }
            else
            {
                fprintf(fpout, "%s", buf);
            }
        }
        gmx_ffclose(fpin);
        if (nsol)
        {
            fprintf(stdout, "Adding line for %d solvent molecules to "
                    "topology file (%s)\n", nsol, topinout);
            fprintf(fpout, "%-15s %5d\n", "SOL", nsol);
        }
        gmx_ffclose(fpout);
        /* use gmx_ffopen to generate backup of topinout */
        fpout = gmx_ffopen(topinout, "w");
        gmx_ffclose(fpout);
        rename(temporary_filename, topinout);
    }
}

int gmx_solvate(int argc, char *argv[])
{
    const char *desc[] = {
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
        "by [TT]-scale[tt]. If radii are not found in the database, those"
        "atoms are assigned the (pre-scaled) distance [TT]-radius[tt].[PAR]",

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

    const char *bugs[] = {
        "Molecules must be whole in the initial configurations.",
    };

    /* parameter data */
    gmx_bool       bProt, bBox;
    const char    *conf_prot, *confout;
    gmx_atomprop_t aps;

    /* protein configuration data */
    char    *title = NULL;
    t_atoms *atoms;
    rvec    *x    = NULL, *v = NULL;
    int      ePBC = -1;
    matrix   box;

    t_filenm fnm[] = {
        { efSTX, "-cp", "protein", ffOPTRD },
        { efSTX, "-cs", "spc216",  ffLIBRD},
        { efSTO, NULL,  NULL,      ffWRITE},
        { efTOP, NULL,  NULL,      ffOPTRW},
    };
#define NFILE asize(fnm)

    static real     defaultDistance = 0.105, r_shell = 0, scaleFactor = 0.57;
    static rvec     new_box         = {0.0, 0.0, 0.0};
    static gmx_bool bReadV          = FALSE;
    static int      max_sol         = 0;
    output_env_t    oenv;
    t_pargs         pa[]              = {
        { "-box",    FALSE, etRVEC, {new_box},
          "Box size (in nm)" },
        { "-radius",   FALSE, etREAL, {&defaultDistance},
          "Default van der Waals distance"},
        { "-scale", FALSE, etREAL, {&scaleFactor},
          "Scale factor to multiply Van der Waals radii from the database in share/gromacs/top/vdwradii.dat. The default value of 0.57 yields density close to 1000 g/l for proteins in water." },
        { "-shell",  FALSE, etREAL, {&r_shell},
          "Thickness of optional water layer around solute" },
        { "-maxsol", FALSE, etINT,  {&max_sol},
          "Maximum number of solvent molecules to add if they fit in the box. If zero (default) this is ignored" },
        { "-vel",    FALSE, etBOOL, {&bReadV},
          "Keep velocities from input solute and solvent" },
    };

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, asize(pa), pa,
                           asize(desc), desc, asize(bugs), bugs, &oenv))
    {
        return 0;
    }

    const char *solventFileName = opt2fn("-cs", NFILE, fnm);
    bProt     = opt2bSet("-cp", NFILE, fnm);
    bBox      = opt2parg_bSet("-box", asize(pa), pa);

    /* check input */
    if (!bProt && !bBox)
    {
        gmx_fatal(FARGS, "When no solute (-cp) is specified, "
                  "a box size (-box) must be specified");
    }

    aps = gmx_atomprop_init();

    snew(atoms, 1);
    init_t_atoms(atoms, 0, FALSE);
    if (bProt)
    {
        /* Generate a solute configuration */
        conf_prot = opt2fn("-cp", NFILE, fnm);
        title     = readConformation(conf_prot, atoms, &x,
                                     bReadV ? &v : NULL, &ePBC, box, "solute");
        if (bReadV && !v)
        {
            fprintf(stderr, "Note: no velocities found\n");
        }
        if (atoms->nr == 0)
        {
            fprintf(stderr, "Note: no atoms in %s\n", conf_prot);
            bProt = FALSE;
        }
    }
    if (bBox)
    {
        ePBC = epbcXYZ;
        clear_mat(box);
        box[XX][XX] = new_box[XX];
        box[YY][YY] = new_box[YY];
        box[ZZ][ZZ] = new_box[ZZ];
    }
    if (det(box) == 0)
    {
        gmx_fatal(FARGS, "Undefined solute box.\nCreate one with gmx editconf "
                  "or give explicit -box command line option");
    }

    add_solv(solventFileName, atoms, &x, v ? &v : NULL, ePBC, box,
             aps, defaultDistance, scaleFactor, r_shell, max_sol);

    /* write new configuration 1 to file confout */
    confout = ftp2fn(efSTO, NFILE, fnm);
    fprintf(stderr, "Writing generated configuration to %s\n", confout);
    if (bProt)
    {
        write_sto_conf(confout, title, atoms, x, v, ePBC, box);
    }
    else
    {
        write_sto_conf(confout, "Generated by gmx solvate", atoms, x, v, ePBC, box);
    }

    /* print size of generated configuration */
    fprintf(stderr, "\nOutput configuration contains %d atoms in %d residues\n",
            atoms->nr, atoms->nres);
    update_top(atoms, box, NFILE, fnm, aps);

    gmx_atomprop_destroy(aps);
    sfree(x);
    sfree(v);
    done_atom(atoms);
    sfree(atoms);
    sfree(title);
    output_env_done(oenv);
    done_filenms(NFILE, fnm);

    return 0;
}
