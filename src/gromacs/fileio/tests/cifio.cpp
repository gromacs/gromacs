/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * Tests for mmCIF/PDBx structure reading and writing (cifio).
 *
 * These exercise the round trip through the confio structure-I/O dispatch, and
 * specifically that the extended PDB entry id and the expanded (>3 char)
 * component/ligand names survive a write/read cycle.
 *
 * \author Amr Alhossary
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "gromacs/fileio/cifio.h"

#include <cstring>

#include <filesystem>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/vec.h"

#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

//! The (extended, 12-char) PDB entry id used as the structure title.
const char* const c_extendedEntryId = "pdb_00001abc";
//! A 5-character ligand/component name (expanded CCD format).
const char* const c_longLigandName = "A1LU6";

/*! \brief Fixture building a small topology: ALA, GLY and a 5-char ligand.
 *
 * Two protein residues on chain A and one HETATM ligand on chain B, in an
 * orthorhombic 3x3x3 nm box, titled with an extended (12-char) entry id.
 */
class CifStructureIoTest : public ::testing::Test
{
protected:
    CifStructureIoTest()
    {
        snew(refTop_, 1);
        open_symtab(&refTop_->symtab);
        refTop_->name = put_symtab(&refTop_->symtab, c_extendedEntryId);

        struct AtomSpec
        {
            const char* atomName;
            const char* resName;
            int         resNr;
            char        chain;
            bool        hetatm;
        };
        const AtomSpec spec[] = {
            { "N", "ALA", 1, 'A', false },  { "CA", "ALA", 1, 'A', false },
            { "C", "ALA", 1, 'A', false },  { "N", "GLY", 2, 'A', false },
            { "CA", "GLY", 2, 'A', false }, { "C1", c_longLigandName, 1, 'B', true },
        };
        numAtoms_ = static_cast<int>(sizeof(spec) / sizeof(spec[0]));

        init_t_atoms(&refTop_->atoms, numAtoms_, TRUE);
        int  resind    = -1;
        int  prevResNr = -1;
        char prevChain = 0;
        for (int i = 0; i < numAtoms_; ++i)
        {
            refTop_->atoms.atomname[i] = put_symtab(&refTop_->symtab, spec[i].atomName);
            const bool newResidue      = (spec[i].resNr != prevResNr || spec[i].chain != prevChain);
            if (newResidue)
            {
                ++resind;
            }
            // resind must be set before t_atoms_set_resinfo(), which indexes
            // resinfo through atoms->atom[i].resind.
            refTop_->atoms.atom[i].resind = resind;
            if (newResidue)
            {
                t_atoms_set_resinfo(
                        &refTop_->atoms, i, &refTop_->symtab, spec[i].resName, spec[i].resNr, ' ', 0, spec[i].chain);
                prevResNr = spec[i].resNr;
                prevChain = spec[i].chain;
            }
            refTop_->atoms.pdbinfo[i].type = spec[i].hetatm ? PdbRecordType::Hetatm : PdbRecordType::Atom;
            refTop_->atoms.pdbinfo[i].occup = 1.0;
            refTop_->atoms.pdbinfo[i].bfac  = 0.0;
        }
        refTop_->atoms.nres = resind + 1;
        close_symtab(&refTop_->symtab);

        clear_mat(refBox_);
        refBox_[XX][XX] = 3.0;
        refBox_[YY][YY] = 3.0;
        refBox_[ZZ][ZZ] = 3.0;
        snew(refX_, numAtoms_);
        for (int i = 0; i < numAtoms_; ++i)
        {
            refX_[i][XX] = 0.1 * (i + 1);
            refX_[i][YY] = 0.2 * (i + 1);
            refX_[i][ZZ] = 0.15 * (i + 1);
        }
    }
    ~CifStructureIoTest() override
    {
        done_top(refTop_);
        sfree(refTop_);
        sfree(refX_);
    }

    //! Write the reference structure to a temporary .cif file and return its path.
    std::string writeReference(const char* suffix)
    {
        const std::string cifFile = fileManager_.getTemporaryFilePath(suffix).string();
        write_sto_conf(cifFile, *refTop_->name, &refTop_->atoms, refX_, nullptr, PbcType::Xyz, refBox_);
        return cifFile;
    }

    TestFileManager fileManager_;
    t_topology*     refTop_;
    rvec*           refX_;
    matrix          refBox_;
    int             numAtoms_;
};

TEST_F(CifStructureIoTest, RoundTripPreservesAtomsCoordinatesAndBox)
{
    const std::string cifFile = writeReference("conf.cif");

    t_topology* top;
    snew(top, 1);
    PbcType pbcType = PbcType::Unset;
    rvec*   x       = nullptr;
    matrix  box;
    read_tps_conf(cifFile, top, &pbcType, &x, nullptr, box, FALSE);

    EXPECT_EQ(top->atoms.nr, numAtoms_);
    for (int i = 0; i < numAtoms_; ++i)
    {
        EXPECT_NEAR(x[i][XX], refX_[i][XX], 1e-3);
        EXPECT_NEAR(x[i][YY], refX_[i][YY], 1e-3);
        EXPECT_NEAR(x[i][ZZ], refX_[i][ZZ], 1e-3);
    }
    for (int d = 0; d < DIM; ++d)
    {
        EXPECT_NEAR(box[d][d], refBox_[d][d], 1e-3);
    }

    sfree(x);
    done_top(top);
    sfree(top);
}

TEST_F(CifStructureIoTest, RoundTripPreservesExtendedEntryId)
{
    const std::string cifFile = writeReference("id.cif");

    t_topology* top;
    snew(top, 1);
    PbcType pbcType = PbcType::Unset;
    rvec*   x       = nullptr;
    matrix  box;
    read_tps_conf(cifFile, top, &pbcType, &x, nullptr, box, FALSE);

    ASSERT_NE(top->name, nullptr);
    // The full 12-character extended id must survive, not be truncated to 4.
    EXPECT_STREQ(*top->name, c_extendedEntryId);

    sfree(x);
    done_top(top);
    sfree(top);
}

TEST_F(CifStructureIoTest, RoundTripPreservesLongComponentName)
{
    const std::string cifFile = writeReference("comp.cif");

    t_topology* top;
    snew(top, 1);
    PbcType pbcType = PbcType::Unset;
    rvec*   x       = nullptr;
    matrix  box;
    read_tps_conf(cifFile, top, &pbcType, &x, nullptr, box, FALSE);

    // The last residue is the 5-character ligand; it must not be truncated.
    bool foundLongLigand = false;
    for (int r = 0; r < top->atoms.nres; ++r)
    {
        if (std::strcmp(*top->atoms.resinfo[r].name, c_longLigandName) == 0)
        {
            foundLongLigand = true;
        }
    }
    EXPECT_TRUE(foundLongLigand) << "5-character component name '" << c_longLigandName
                                 << "' was lost on read/write";

    sfree(x);
    done_top(top);
    sfree(top);
}

} // namespace
} // namespace test
} // namespace gmx
