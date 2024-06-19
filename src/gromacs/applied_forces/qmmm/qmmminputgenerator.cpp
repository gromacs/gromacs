/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Implements input generator class for CP2K QMMM
 *
 * \author Dmitry Morozov <dmitry.morozov@jyu.fi>
 * \ingroup module_applied_forces
 */

#include "gmxpre.h"

#include "qmmminputgenerator.h"

#include <cmath>
#include <cstddef>

#include <vector>

#include "gromacs/applied_forces/qmmm/qmmmtypes.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/stringutil.h"

enum class PbcType : int;

namespace gmx
{

QMMMInputGenerator::QMMMInputGenerator(const QMMMParameters& parameters,
                                       PbcType               pbcType,
                                       const matrix          box,
                                       ArrayRef<const real>  q,
                                       ArrayRef<const RVec>  x) :
    parameters_(parameters),
    pbc_(pbcType),
    qmBox_{ { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } },
    qmCenter_{ 0.0, 0.0, 0.0 },
    qmTrans_{ 0.0, 0.0, 0.0 },
    q_(q),
    x_(x)
{
    // Copy box
    copy_mat(box, box_);

    // Fill qmAtoms_ set
    for (const auto& val : parameters_.qmIndices_)
    {
        qmAtoms_.emplace(val);
    }

    // Compute QM box
    computeQMBox(sc_qmBoxScale, sc_qmBoxMinLength);
}

bool QMMMInputGenerator::isQMAtom(Index globalAtomIndex) const
{
    return (qmAtoms_.find(globalAtomIndex) != qmAtoms_.end());
}

void QMMMInputGenerator::computeQMBox(real scale, real minNorm)
{
    // Init atom numbers
    size_t nQm = parameters_.qmIndices_.size();

    // If there is only one QM atom, then just copy the box_
    if (nQm < 2)
    {
        copy_mat(box_, qmBox_);
        qmCenter_ = x_[parameters_.qmIndices_[0]];
        qmTrans_  = RVec(qmBox_[0]) / 2 + RVec(qmBox_[1]) / 2 + RVec(qmBox_[2]) / 2 - qmCenter_;
        return;
    }

    // Initialize pbc
    t_pbc pbc;
    set_pbc(&pbc, pbc_, box_);

    /* To compute qmBox_:
     * 1) Compute maximum dimension - maxDist betweeen QM atoms within PBC
     * 2) Make projection of the each box_ vector onto the cross prod of other two vectors
     * 3) Calculate scales so that norm will be maxDist for each box_ vector
     * 4) Apply scales to the box_ to get qmBox_
     */
    RVec dx(0.0, 0.0, 0.0);
    real maxDist = 0.0;

    // Search for the maxDist - maximum distance within QM system
    for (size_t i = 0; i < nQm - 1; i++)
    {
        for (size_t j = i + 1; j < nQm; j++)
        {
            pbc_dx(&pbc, x_[parameters_.qmIndices_[i]], x_[parameters_.qmIndices_[j]], dx);
            maxDist = dx.norm() > maxDist ? dx.norm() : maxDist;
        }
    }

    // Apply scale factor: qmBox_ should be *scale times bigger than maxDist
    maxDist *= scale;

    // Compute QM Box vectors
    RVec vec0 = computeQMBoxVec(box_[0], box_[1], box_[2], maxDist, minNorm, norm(box_[0]));
    RVec vec1 = computeQMBoxVec(box_[1], box_[0], box_[2], maxDist, minNorm, norm(box_[1]));
    RVec vec2 = computeQMBoxVec(box_[2], box_[0], box_[1], maxDist, minNorm, norm(box_[2]));

    // Form qmBox_
    copy_rvec(vec0, qmBox_[0]);
    copy_rvec(vec1, qmBox_[1]);
    copy_rvec(vec2, qmBox_[2]);

    /* Now we need to also compute translation vector.
     * In order to center QM atoms in the computed qmBox_
     *
     * First compute center of QM system by averaging PBC-aware distance vectors
     * with respect to the first QM atom.
     */
    for (size_t i = 1; i < nQm; i++)
    {
        // compute pbc-aware distance vector between QM atom 0 and QM atom i
        pbc_dx(&pbc, x_[parameters_.qmIndices_[i]], x_[parameters_.qmIndices_[0]], dx);

        // add that to the qmCenter_
        qmCenter_ += dx;
    }

    // Average over all nQm atoms and add first atom coordiantes
    qmCenter_ = qmCenter_ / nQm + x_[parameters_.qmIndices_[0]];

    // Translation vector will be center of qmBox_ - qmCenter_
    qmTrans_ = vec0 / 2 + vec1 / 2 + vec2 / 2 - qmCenter_;
}

std::string QMMMInputGenerator::generateGlobalSection()
{
    std::string res;

    res += "&GLOBAL\n";
    res += "  PRINT_LEVEL LOW\n";
    res += "  PROJECT GROMACS\n";
    res += "  RUN_TYPE ENERGY_FORCE\n";
    res += "&END GLOBAL\n";

    return res;
}

std::string QMMMInputGenerator::generateDFTSection() const
{
    std::string res;

    res += "  &DFT\n";

    // write charge and multiplicity
    res += formatString("    CHARGE %d\n", parameters_.qmCharge_);
    res += formatString("    MULTIPLICITY %d\n", parameters_.qmMultiplicity_);

    // If multiplicity is not 1 then we should use unrestricted Kohn-Sham
    if (parameters_.qmMultiplicity_ > 1)
    {
        res += "    UKS\n";
    }

    // Basis files, Grid setup and SCF parameters
    res += "    BASIS_SET_FILE_NAME  BASIS_MOLOPT\n";
    res += "    POTENTIAL_FILE_NAME  POTENTIAL\n";
    res += "    &MGRID\n";
    res += "      NGRIDS 5\n";
    res += "      CUTOFF 450\n";
    res += "      REL_CUTOFF 50\n";
    res += "      COMMENSURATE\n";
    res += "    &END MGRID\n";
    res += "    &SCF\n";
    res += "      SCF_GUESS RESTART\n";
    res += "      EPS_SCF 5.0E-8\n";
    res += "      MAX_SCF 20\n";
    res += "      &OT  T\n";
    res += "        MINIMIZER  DIIS\n";
    res += "        STEPSIZE   0.15\n";
    res += "        PRECONDITIONER FULL_ALL\n";
    res += "      &END OT\n";
    res += "      &OUTER_SCF  T\n";
    res += "        MAX_SCF 20\n";
    res += "        EPS_SCF 5.0E-8\n";
    res += "      &END OUTER_SCF\n";
    res += "    &END SCF\n";

    // DFT functional parameters
    res += "    &XC\n";
    res += "      DENSITY_CUTOFF     1.0E-12\n";
    res += "      GRADIENT_CUTOFF    1.0E-12\n";
    res += "      TAU_CUTOFF         1.0E-12\n";
    res += formatString("      &XC_FUNCTIONAL %s\n", c_qmmmQMMethodNames[parameters_.qmMethod_]);
    res += "      &END XC_FUNCTIONAL\n";
    res += "    &END XC\n";
    res += "    &QS\n";
    res += "     METHOD GPW\n";
    res += "     EPS_DEFAULT 1.0E-10\n";
    res += "     EXTRAPOLATION ASPC\n";
    res += "     EXTRAPOLATION_ORDER  4\n";
    res += "    &END QS\n";
    res += "  &END DFT\n";

    return res;
}

std::string QMMMInputGenerator::generateQMMMSection() const
{
    std::string res;

    // Init some numbers
    const size_t nQm = parameters_.qmIndices_.size();

    // Count the numbers of individual QM atoms per type
    std::vector<int> num_atoms(periodic_system.size(), 0);
    for (size_t i = 0; i < nQm; i++)
    {
        num_atoms[parameters_.atomNumbers_[parameters_.qmIndices_[i]]]++;
    }

    res += "  &QMMM\n";

    // QM cell
    res += "    &CELL\n";
    res += formatString(
            "      A %.3lf %.3lf %.3lf\n", qmBox_[0][0] * 10, qmBox_[0][1] * 10, qmBox_[0][2] * 10);
    res += formatString(
            "      B %.3lf %.3lf %.3lf\n", qmBox_[1][0] * 10, qmBox_[1][1] * 10, qmBox_[1][2] * 10);
    res += formatString(
            "      C %.3lf %.3lf %.3lf\n", qmBox_[2][0] * 10, qmBox_[2][1] * 10, qmBox_[2][2] * 10);
    res += "      PERIODIC XYZ\n";
    res += "    &END CELL\n";

    res += "    CENTER EVERY_STEP\n";
    res += "    CENTER_GRID TRUE\n";
    res += "    &WALLS\n";
    res += "      TYPE REFLECTIVE\n";
    res += "    &END WALLS\n";

    res += "    ECOUPL GAUSS\n";
    res += "    USE_GEEP_LIB 12\n";
    res += "    &PERIODIC\n";
    res += "      GMAX     1.0E+00\n";
    res += "      &MULTIPOLE ON\n";
    res += "         RCUT     1.0E+01\n";
    res += "         EWALD_PRECISION     1.0E-06\n";
    res += "      &END\n";
    res += "    &END PERIODIC\n";

    // Print indicies of QM atoms
    // Loop over counter of QM atom types
    for (size_t i = 0; i < num_atoms.size(); i++)
    {
        if (num_atoms[i] > 0)
        {
            res += formatString("    &QM_KIND %3s\n", periodic_system[i].c_str());
            res += "      MM_INDEX";
            // Loop over all QM atoms indexes
            for (size_t j = 0; j < nQm; j++)
            {
                if (parameters_.atomNumbers_[parameters_.qmIndices_[j]] == static_cast<Index>(i))
                {
                    res += formatString(" %d", static_cast<int>(parameters_.qmIndices_[j] + 1));
                }
            }

            res += "\n";
            res += "    &END QM_KIND\n";
        }
    }

    // Print &LINK groups
    // Loop over parameters_.link_
    for (size_t i = 0; i < parameters_.link_.size(); i++)
    {
        res += "    &LINK\n";
        res += formatString("      QM_INDEX %d\n", static_cast<int>(parameters_.link_[i].qm) + 1);
        res += formatString("      MM_INDEX %d\n", static_cast<int>(parameters_.link_[i].mm) + 1);
        res += "    &END LINK\n";
    }

    res += "  &END QMMM\n";

    return res;
}

std::string QMMMInputGenerator::generateMMSection()
{
    std::string res;

    res += "  &MM\n";
    res += "    &FORCEFIELD\n";
    res += "      DO_NONBONDED FALSE\n";
    res += "    &END FORCEFIELD\n";
    res += "    &POISSON\n";
    res += "      &EWALD\n";
    res += "        EWALD_TYPE NONE\n";
    res += "      &END EWALD\n";
    res += "    &END POISSON\n";
    res += "  &END MM\n";

    return res;
}

std::string QMMMInputGenerator::generateSubsysSection() const
{
    std::string res;

    // Init some numbers
    size_t nQm = parameters_.qmIndices_.size();
    size_t nMm = parameters_.mmIndices_.size();
    size_t nAt = nQm + nMm;

    // Count the numbers of individual QM atoms per type
    std::vector<int> num_atoms(periodic_system.size(), 0);
    for (size_t i = 0; i < nQm; i++)
    {
        num_atoms[parameters_.atomNumbers_[parameters_.qmIndices_[i]]]++;
    }

    res += "  &SUBSYS\n";

    // Print cell parameters
    res += "    &CELL\n";
    res += formatString("      A %.3lf %.3lf %.3lf\n", box_[0][0] * 10, box_[0][1] * 10, box_[0][2] * 10);
    res += formatString("      B %.3lf %.3lf %.3lf\n", box_[1][0] * 10, box_[1][1] * 10, box_[1][2] * 10);
    res += formatString("      C %.3lf %.3lf %.3lf\n", box_[2][0] * 10, box_[2][1] * 10, box_[2][2] * 10);
    res += "      PERIODIC XYZ\n";
    res += "    &END CELL\n";

    // Print topology section
    res += "    &TOPOLOGY\n";

    // Placeholder for PDB file name that will be replaced with actuall name during mdrun
    res += "      COORD_FILE_NAME %s\n";

    res += "      COORD_FILE_FORMAT PDB\n";
    res += "      CHARGE_EXTENDED TRUE\n";
    res += "      CONNECTIVITY OFF\n";
    res += "      &GENERATE\n";
    res += "         &ISOLATED_ATOMS\n";
    res += formatString("            LIST %d..%d\n", 1, static_cast<int>(nAt));
    res += "         &END\n";
    res += "      &END GENERATE\n";
    res += "    &END TOPOLOGY\n";

    // Now we will print basises for all types of QM atoms
    // Loop over counter of QM atom types
    for (size_t i = 0; i < num_atoms.size(); i++)
    {
        if (num_atoms[i] > 0)
        {
            res += "    &KIND " + periodic_system[i] + "\n";
            res += "      ELEMENT " + periodic_system[i] + "\n";
            res += "      BASIS_SET DZVP-MOLOPT-GTH\n";
            res += "      POTENTIAL GTH-" + std::string(c_qmmmQMMethodNames[parameters_.qmMethod_]);
            res += "\n";
            res += "    &END KIND\n";
        }
    }

    // Add element kind X - they are represents virtual sites
    res += "    &KIND X\n";
    res += "      ELEMENT H\n";
    res += "    &END KIND\n";
    res += "  &END SUBSYS\n";

    return res;
}


std::string QMMMInputGenerator::generateCP2KInput() const
{

    std::string inp;

    // Begin CP2K input generation
    inp += generateGlobalSection();

    inp += "&FORCE_EVAL\n";
    inp += "  METHOD QMMM\n";

    // Add DFT parameters
    inp += generateDFTSection();

    // Add QMMM parameters
    inp += generateQMMMSection();

    // Add MM parameters
    inp += generateMMSection();

    // Add SUBSYS parameters
    inp += generateSubsysSection();

    inp += "&END FORCE_EVAL\n";

    return inp;
}

std::string QMMMInputGenerator::generateCP2KPdb() const
{

    std::string pdb;

    /* Generate *.pdb formatted lines
     * and append to std::string pdb
     */
    for (size_t i = 0; i < x_.size(); i++)
    {
        pdb += "ATOM  ";

        // Here we need to print i % 100000 because atom counter in *.pdb has only 5 digits
        pdb += formatString("%5d ", static_cast<int>(i % 100000));

        // Atom name
        pdb += formatString(" %3s ", periodic_system[parameters_.atomNumbers_[i]].c_str());

        // Lable atom as QM or MM residue
        if (isQMAtom(i))
        {
            pdb += " QM     1     ";
        }
        else
        {
            pdb += " MM     2     ";
        }

        // Coordinates
        pdb += formatString("%7.3lf %7.3lf %7.3lf  1.00  0.00         ",
                            (x_[i][XX] + qmTrans_[XX]) * 10,
                            (x_[i][YY] + qmTrans_[YY]) * 10,
                            (x_[i][ZZ] + qmTrans_[ZZ]) * 10);

        // Atom symbol
        pdb += formatString(" %3s ", periodic_system[parameters_.atomNumbers_[i]].c_str());

        // Point charge for MM atoms or 0 for QM atoms
        pdb += formatString("%lf\n", q_[i]);
    }

    return pdb;
}

const RVec& QMMMInputGenerator::qmTrans() const
{
    return qmTrans_;
}

const matrix& QMMMInputGenerator::qmBox() const
{
    return qmBox_;
}


RVec computeQMBoxVec(const RVec& a, const RVec& b, const RVec& c, real h, real minNorm, real maxNorm)
{
    RVec vec0(a);
    RVec vec1(b);
    RVec vec2(c);
    RVec dx(0.0, 0.0, 0.0);
    RVec res(0.0, 0.0, 0.0);

    // Compute scales sc that will transform a
    dx = vec1.cross(vec2);
    dx /= dx.norm();

    // Transform a
    res = h / std::fabs(vec0.dot(dx)) * a;

    // If vector is smaller than minL then scale it up
    if (res.norm() < minNorm)
    {
        res *= minNorm / res.norm();
    }

    // If vector is longer than maxL then scale it down
    if (res.norm() > maxNorm)
    {
        res *= maxNorm / res.norm();
    }

    return res;
}

} // namespace gmx
