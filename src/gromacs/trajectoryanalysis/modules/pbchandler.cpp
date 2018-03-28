/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Helper classes to handle PBC changes
 *
 * \author
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "pbchandler.h"

#include <algorithm>

#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"

namespace gmx
{

void
Pbchandler::initPBCOptions(IOptionsContainer *options)
{
    options->addOption(SelectionOption("selcluster").store(&selClust_).dynamicMask()
                           .storeIsSet(&bHaveSelClust_)
                           .description("Reference selection for clustering"));

    options->addOption(SelectionOption("selcenter").store(&selCent_).dynamicMask()
                           .storeIsSet(&bHaveSelCent_)
                           .description("Selection for centering"));

    options->addOption(SelectionOption("selcom").store(&selCOM_).dynamicMask()
                           .storeIsSet(&bHaveSelCom_)
                           .description("Selection for COM fitting, default is ouput selection"));

    options->addOption(EnumOption<PBCType>("setpbc").store(&selPBCType_)
                           .enumValue(c_pbcTypes).defaultValue(ePBCTypeWhole)
                           .description("Type of output PBC option"));

    options->addOption(EnumOption<UnitCellType>("ur").store(&selUnitCellType_)
                           .enumValue(c_unitCellTypes).defaultValueIfSet(eUnitCellTypeCompact)
                           .description("Type of Unit Cell representation for the output"));

    options->addOption(EnumOption<CenteringType>("boxcenter").store(&selCenteringType_)
                           .enumValue(c_centeringTypes).storeIsSet(&bCenter_).defaultValueIfSet(eCenteringTypeTriclinic)
                           .description("How should the system be centered"));
}


void
Pbchandler::checkOptions(const Selection *sel, TrajectoryAnalysisSettings *settings)
{
    if ((selCenteringType_ == eCenteringNotSet || selCenteringType_ ==  eCenteringTypeZero) &&
        !bHaveSelCent_)
    {
        if (selCenteringType_ != eCenteringTypeZero)
        {
            GMX_THROW(InconsistentInputError("Box centering requested but no selection provided for it with -selcenter"));
        }
        selCent_ = *sel;
    }
    if (selPBCType_ == ePBCTypeCluster &&
        !bHaveSelClust_)
    {
        GMX_THROW(InconsistentInputError("PBC clustering requested but no selection provided for it with -selcluster"));
    }
    if ((selPBCType_ == ePBCTypeCluster || selPBCType_ == ePBCTypeMolecule || selPBCType_ == ePBCTypeResidue) &&
        selCenteringType_ == eCenteringNotSet)
    {
        selCenteringType_ = eCenteringTypeZero;
    }
    if (selPBCType_ == ePBCTypeWhole || selPBCType_ == ePBCTypeMolecule || selPBCType_ == ePBCTypeResidue)
    {
        settings->setPBC(true);
    }
    else
    {
        settings->setPBC(false);
    }
    if ((selPBCType_ == ePBCTypeMolecule || selPBCType_ == ePBCTypeResidue || selPBCType_ == ePBCTypeAtom) && !bHaveSelCom_)
    {
        selCOM_ = *sel;
    }
    sel_ = *sel;
}

void
Pbchandler::doPBC(const t_trxframe *input, t_trxframe *output, const gmx_mtop_t *mtop)
{
    switch (selPBCType_)
    {
        case (ePBCTypeNoJump):
            removeJump(output);
            centerSystem(output, selCenteringType_);
            break;
        case (ePBCTypeCluster):
            cluster(output, mtop);
            centerSystem(output, selCenteringType_);
            break;
        case (ePBCTypeAtom):
            centerSystem(output, selCenteringType_);
            comAtom(output);
            break;
        case (ePBCTypeResidue):
            centerSystem(output, selCenteringType_);
            comResidue(input, output, mtop);
            break;
        case (ePBCTypeMolecule):
            centerSystem(output, selCenteringType_);
            comMolecule(input, output, mtop);
            break;
        case (ePBCTypeWhole):
        case (ePBCNotSet):
            // do nothing
            break;
        default:
            GMX_THROW(InconsistentInputError("No such PBC type"));
    }
}


void
Pbchandler::shiftCoord(const int atomStart, const int atomEnd, rvec shift, rvec *coord, const bool checkValid)
{
    ArrayRef<const int> indices;
    // if we need to check the indicies, means we are operating on global (input) data,
    // not already sanitized data from only the output frame
    if (checkValid)
    {
        indices = sel_.atomIndices();
        for (int i = atomStart; i < atomEnd; i++)
        {
            int pos = std::find(indices.begin(), indices.end(), i) - indices.begin();
            if (pos < indices.size())
            {
                rvec_inc(coord[pos], shift);
            }
        }
    }
    // means we are only operating on the already corrected number of output coordinates
    else
    {
        for (int i = atomStart; i < atomEnd; i++)
        {
            rvec_inc(coord[i], shift);
        }

    }
}



void
Pbchandler::setReferenceCoordinates(const TopologyInformation &top)
{
    t_trxframe                 fr;
    clear_trxframe(&fr, TRUE);
    const TopologyInformation *localTop = &top;
    localTop->getTopologyConf(&fr.x, fr.box);
    const int                  natoms = localTop->mtop()->natoms;
    refCoord_.resize(natoms);

    for (int i = 0; (i < natoms); i++)
    {
        refCoord_[i] = fr.x[i];
    }
}

void
Pbchandler::centerSystem(t_trxframe *output, int centeringType)
{
    if (centeringType == eCenteringNotSet || centeringType == eCenteringTypeZero)
    {
        return;
    }

    Selection sel    = selCent_;
    int       natoms = sel.atomCount();

    RVec      cmin = refCoord_[0];
    RVec      cmax = cmin;
    matrix    box;
    copy_mat(output->box, box);
    // we are iterating over all atoms the selCom selection, but need to know
    // positions from refCoord for iteration
    for (int i = 0; i < natoms; i++)
    {
        for (int m = 0; m < DIM; m++)
        {
            if (output->x[i][m] < cmin[m])
            {
                cmin[m] = output->x[i][m];
            }
            else if (output->x[i][m] > cmax[m])
            {
                cmax[m] = output->x[i][m];
            }
        }
    }
    RVec boxCenter(0, 0, 0);
    calc_box_center(centeringType, box, boxCenter);
    RVec dx(0, 0, 0);
    for (int m = 0; m < DIM; m++)
    {
        dx[m] = boxCenter[m]-(cmin[m]+cmax[m])*0.5;
    }
    // now only work on atoms from final output
    natoms = sel_.atomCount();
    shiftCoord(0, natoms, dx, output->x, false);
}


void
Pbchandler::removeJump(t_trxframe *output)
{
    RVec hbox(0, 0, 0);
    // load topology coordinates

    for (int d = 0; d < DIM; d++)
    {
        hbox[d] = 0.5*output->box[d][d];
    }
    int natoms = sel_.atomCount();
    for (int i = 0; i < natoms; i++)
    {
        int pos = sel_.position(i).refId();
        for (int m = DIM-1; m >= 0; m--)
        {
            if (hbox[m] > 0)
            {
                while (output->x[i][m]-refCoord_[pos][m] <= -hbox[m])
                {
                    for (int d = 0; d <= m; d++)
                    {
                        output->x[i][d] += output->box[m][d];
                    }
                }
                while (output->x[i][m]-refCoord_[pos][m] > hbox[m])
                {
                    for (int d = 0; d <= m; d++)
                    {
                        output->x[i][d] -= output->box[m][d];
                    }
                }
            }
        }
    }
}

void
Pbchandler::cluster(t_trxframe *output, const gmx_mtop_t *mtop)
{
    t_pbc     pbc;

    RVec      boxCenter(0, 0, 0);

    calc_box_center(selCenteringType_, output->box, boxCenter);

    /* Initiate the pbc structure */
    std::memset(&pbc, 0, sizeof(pbc));
    set_pbc(&pbc, output->ePBC, output->box);

    /* Convert atom index to molecular */
    int              *molind = mtop->mols.index;

    int               allAtoms = mtop->natoms;
    int               allMols  = mtop->mols.nr;

    std::vector<bool> bTmp;
    bTmp.resize(allAtoms);
    std::vector<bool> bMol;
    bMol.resize(allMols);
    std::vector<int>  cluster;

    for (int i = 0; (i < allAtoms); i++)
    {
        /* Mark all molecules in the index */
        int pos       = selClust_.position(i).refId();
        bTmp[pos] = TRUE;
        /* Binary search assuming the molecules are sorted */
        int j0 = 0;
        int j1 = allMols-1;
        while (j0 < j1)
        {
            if (pos < molind[j0+1])
            {
                j1 = j0;
            }
            else if (pos >= molind[j1])
            {
                j0 = j1;
            }
            else
            {
                int jj = (j0+j1)/2;
                if (pos < molind[jj+1])
                {
                    j1 = jj;
                }
                else
                {
                    j0 = jj;
                }
            }
        }
        bMol[j0] = TRUE;
    }
    /* Double check whether all atoms in all molecules that are marked are part
     * of the cluster. Simultaneously compute the center of geometry.
     */
    real              min_dist2   = 10*gmx::square(trace(output->box));
    int               imol_center = -1;
    int               ncluster    = 0;
    // local buffer for coordinates of complete system.
    std::vector<RVec> local;
    local.resize(allAtoms);
    std::vector<RVec> mCom;
    std::vector<RVec> mShift;
    for (int i = 0; i < allMols; i++)
    {
        mCom.push_back(RVec(0, 0, 0));
        mShift.push_back(RVec(0, 0, 0));
        for (int j = molind[i]; j < molind[i+1]; j++)
        {
            local[j] = output->x[j];
            if (bMol[i] && !bTmp[j])
            {
                gmx_fatal(FARGS, "Molecule %d marked for clustering but not atom %d in it - check your index!", i+1, j+1);
            }
            else if (!bMol[i] && bTmp[j])
            {
                gmx_fatal(FARGS, "Atom %d marked for clustering but not molecule %d - this is an internal error...", j+1, i+1);
            }
            else if (bMol[i])
            {
                /* Make molecule whole, move 2nd and higher atom to same periodicity as 1st atom in molecule */
                if (j > molind[i])
                {
                    RVec dx(0, 0, 0);
                    pbc_dx(&pbc, local[j], local[j-1], dx);
                    rvec_add(local[j-1], dx, local[j]);
                }
                /* Compute center of geometry of molecule - mCom[i] was zeroed when we did snew() on it! */
                rvec_inc(mCom[i], local[j]);
            }
        }
        if (bMol[i])
        {
            /* Normalize center of geometry */
            real fac = 1.0/(molind[i+1]-molind[i]);
            for (int m = 0; (m < DIM); m++)
            {
                mCom[i][m] *= fac;
            }
            /* Determine which molecule is closest to the center of the box */
            RVec dx(0, 0, 0);
            pbc_dx(&pbc, boxCenter, mCom[i], dx);
            real tmp_r2 = iprod(dx, dx);

            if (tmp_r2 < min_dist2)
            {
                min_dist2   = tmp_r2;
                imol_center = i;
            }
            cluster.push_back(i);
            ncluster++;
        }
    }

    if (ncluster <= 0)
    {
        fprintf(stderr, "No molecules selected in the cluster\n");
        return;
    }
    else if (imol_center == -1)
    {
        fprintf(stderr, "No central molecules could be found\n");
        return;
    }

    int              nadded            = 0;
    std::vector<int> added;
    added.push_back(imol_center);
    nadded++;
    bMol[imol_center] = FALSE;

    while (nadded < ncluster)
    {
        /* Find min distance between cluster molecules and those remaining to be added */
        min_dist2   = 10*gmx::square(trace(output->box));
        int imin        = -1;
        int jmin        = -1;
        /* Loop over added mols */
        for (int i = 0; i < nadded; i++)
        {
            int ai = added[i];
            /* Loop over all mols */
            for (int j = 0; j < ncluster; j++)
            {
                int aj = cluster[j];
                /* check those remaining to be added */
                if (bMol[aj])
                {
                    RVec dx(0, 0, 0);
                    pbc_dx(&pbc, mCom[aj], mCom[ai], dx);
                    real tmp_r2 = iprod(dx, dx);
                    if (tmp_r2 < min_dist2)
                    {
                        min_dist2   = tmp_r2;
                        imin        = ai;
                        jmin        = aj;
                    }
                }
            }
        }

        /* Add the best molecule */
        added.push_back(jmin);
        nadded++;
        bMol[jmin]        = FALSE;
        /* Calculate the shift from the ai molecule */
        RVec dx(0, 0, 0);
        pbc_dx(&pbc, mCom[jmin], mCom[imin], dx);
        RVec xtest(0, 0, 0);
        rvec_add(mCom[imin], dx, xtest);
        rvec_sub(xtest, mCom[jmin], mShift[jmin]);
        rvec_inc(mCom[jmin], mShift[jmin]);

        shiftCoord(molind[jmin], molind[jmin+1], mShift[jmin], gmx::as_rvec_array(local.data()), false);
        fprintf(stdout, "\rClustering iteration %d of %d...", nadded, ncluster);
        fflush(stdout);
    }

    fprintf(stdout, "\n");
}

void
Pbchandler::comAtom(t_trxframe *output)
{
    auto positionsArrayRef = arrayRefFromArray(reinterpret_cast<gmx::RVec *>(output->x), output->natoms);
    switch (selUnitCellType_)
    {
        case eUnitCellTypeRectangular:
            put_atoms_in_box(output->ePBC, output->box, positionsArrayRef);
            break;
        case eUnitCellTypeTriclinic:
            put_atoms_in_triclinic_unitcell(selCenteringType_, output->box, positionsArrayRef);
            break;
        case eUnitCellTypeCompact:
            put_atoms_in_compact_unitcell(output->ePBC, selCenteringType_, output->box, positionsArrayRef);
            break;
        case eUnitCellNotSet:
            // handle this somehow
            break;
    }
}

void
Pbchandler::comMolecule(const t_trxframe *input, t_trxframe *output, const gmx_mtop_t *mtop)
{
    t_pbc   pbc;
    RVec    boxCenter(0, 0, 0);

    calc_box_center(selCenteringType_, output->box, boxCenter);
    set_pbc(&pbc, output->ePBC, output->box);
    int       allAtoms = mtop->natoms;
    int       nMols    = mtop->mols.nr;
    int      *molind   = mtop->mols.index;
    // calculate com only on atoms selected for it
    Selection sel = selCOM_;

    for (int i = 0; (i < nMols); i++)
    {
        /* calc COM
         *
         * Need to make sure we stop iterating over molecules after we have included all atoms that the COM
         * should be calculated for, or this becomes insanely slow
         *
         * */
        RVec                com(0, 0, 0);
        double              mtot       = 0;
        int                 mblock     = 0;
        ArrayRef<const int> comIndices = sel.atomIndices();
        for (int j = molind[i]; (j < molind[i+1] && j < allAtoms); j++)
        {
            if (std::binary_search(comIndices.begin(), comIndices.end(), j))
            {
                real m = mtopGetAtomMass(mtop, j, &mblock);


                for (int d = 0; d < DIM; d++)
                {
                    com[d] += m*input->x[j][d];
                }
                mtot += m;
            }
        }
        /* calculate final COM */
        if (mtot != 0)
        {
            svmul(1.0/mtot, com, com);
        }

        /* check if COM is outside box */
        RVec      newCom(com);
        auto      newComArrayRef = gmx::arrayRefFromArray(&newCom, 1);
        switch (selUnitCellType_)
        {
            case eUnitCellTypeRectangular:
                put_atoms_in_box(output->ePBC, output->box, newComArrayRef);
                break;
            case eUnitCellTypeTriclinic:
                put_atoms_in_triclinic_unitcell(selCenteringType_, output->box, newComArrayRef);
                break;
            case eUnitCellTypeCompact:
                put_atoms_in_compact_unitcell(output->ePBC, selCenteringType_, output->box, newComArrayRef);
                break;
        }
        RVec shift(0, 0, 0);
        rvec_sub(newCom, com, shift);
        if (norm2(shift) > 0)
        {
            if (debug)
            {
                fprintf(debug, "\nShifting position of molecule %d "
                        "by %8.3f  %8.3f  %8.3f\n", i+1,
                        shift[XX], shift[YY], shift[ZZ]);
            }
            shiftCoord(molind[i], molind[i+1], shift, output->x, true);
        }
    }
}

void
Pbchandler::comResidue(const t_trxframe *input, t_trxframe *output, const gmx_mtop_t *mtop)
{
    RVec             boxCenter(0, 0, 0);
    static const int NOTSET = -12347;
    calc_box_center(selCenteringType_, output->box, boxCenter);

    int    res_start = 0;
    int    res_end   = 0;
    RVec   com(0, 0, 0);
    double mtot = 0;

    int    natoms = selCOM_.atomCount();
    int    presnr = NOTSET;
    for (int i = 0; i < natoms+1; i++)
    {
        int molb = 0;
        int pos  = 0;
        if (i < natoms)
        {
            pos = selCOM_.position(i).refId();
        }
        t_atom atom = mtopGetAtomParameters(mtop, pos, &molb);
        if (i == natoms || (presnr != atom.resind && presnr != NOTSET))
        {
            /* calculate final COM */
            res_end = i;
            if (mtot != 0)
            {
                svmul(1.0/mtot, com, com);
            }

            /* check if COM is outside box */
            RVec      newCom(com);
            auto      newComArrayRef = gmx::arrayRefFromArray(&newCom, 1);
            switch (selUnitCellType_)
            {
                case eUnitCellTypeRectangular:
                    put_atoms_in_box(output->ePBC, output->box, newComArrayRef);
                    break;
                case eUnitCellTypeTriclinic:
                    put_atoms_in_triclinic_unitcell(selCenteringType_, output->box, newComArrayRef);
                    break;
                case eUnitCellTypeCompact:
                    put_atoms_in_compact_unitcell(output->ePBC, selCenteringType_, output->box, newComArrayRef);
                    break;
            }
            RVec shift(0, 0, 0);
            rvec_sub(newCom, com, shift);
            if (norm2(shift))
            {
                if (debug)
                {
                    fprintf(debug, "\nShifting position of residue %d (atoms %d-%d) "
                            "by %g,%g,%g\n", atom.resind,
                            res_start+1, res_end+1, shift[XX], shift[YY], shift[ZZ]);
                }
                shiftCoord(res_start, res_end, shift, output->x, true);
            }
            clear_rvec(com);
            mtot = 0;

            /* remember start of new residue */
            res_start = i;
        }
        if (i < natoms)
        {
            /* calc COM */
            real m = atom.m;
            for (int d = 0; d < DIM; d++)
            {
                com[d] += m*input->x[pos][d];
            }
            mtot += m;

            presnr = atom.resind;
        }
    }
}

} // namespace gmx
