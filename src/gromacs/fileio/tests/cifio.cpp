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
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/pdbio.h"
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

public:
    ~CifStructureIoTest() override
    {
        done_top(refTop_);
        sfree(refTop_);
        sfree(refX_);
    }

protected:
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

TEST(CifResidueTest, ChainChangeStartsNewResidue)
{
    // Two ligands with the same name and number on adjacent rows, as when each
    // chain numbers its own ligands; only the chain tells them apart.
    TestFileManager   fileManager;
    const std::string cifFile = fileManager.getTemporaryFilePath("chains.cif").string();
    {
        std::ofstream out(cifFile);
        out << "data_chains\n"
               "loop_\n"
               "_atom_site.group_PDB\n"
               "_atom_site.id\n"
               "_atom_site.type_symbol\n"
               "_atom_site.label_atom_id\n"
               "_atom_site.label_comp_id\n"
               "_atom_site.auth_asym_id\n"
               "_atom_site.auth_seq_id\n"
               "_atom_site.Cartn_x\n"
               "_atom_site.Cartn_y\n"
               "_atom_site.Cartn_z\n"
               "HETATM 1 C C1 NAG A 401 1.0 1.0 1.0\n"
               "HETATM 2 C C1 NAG H 401 2.0 2.0 2.0\n";
    }

    t_topology* top;
    snew(top, 1);
    PbcType pbcType = PbcType::Unset;
    rvec*   x       = nullptr;
    matrix  box;
    read_tps_conf(cifFile, top, &pbcType, &x, nullptr, box, FALSE);

    ASSERT_EQ(top->atoms.nr, 2);
    ASSERT_EQ(top->atoms.nres, 2);
    EXPECT_EQ(top->atoms.atom[1].resind, 1);
    EXPECT_EQ(top->atoms.resinfo[0].chainid, 'A');
    EXPECT_EQ(top->atoms.resinfo[1].chainid, 'H');
    EXPECT_NE(top->atoms.resinfo[0].chainnum, top->atoms.resinfo[1].chainnum);

    sfree(x);
    done_top(top);
    sfree(top);
}

/*! \brief One atom of the small structures used by the connectivity tests.
 *
 * These are written directly with CifIO rather than through write_sto_conf(),
 * because only CifIO takes connectivity.
 */
struct ConnectivityAtomSpec
{
    //! Atom name, as written to label_atom_id.
    const char* atomName;
    //! Residue name.
    const char* resName;
    //! Residue number.
    int resNr;
    //! Chain identifier.
    char chain;
    //! Element symbol, or "" when unknown.
    const char* element;
    //! Coordinates in nm.
    rvec x;
};

//! Structure built from a ConnectivityAtomSpec list, with the memory it owns.
class ConnectivityStructure
{
public:
    //! Build the atoms and coordinates described by \p spec.
    explicit ConnectivityStructure(const std::vector<ConnectivityAtomSpec>& spec)
    {
        const int numAtoms = static_cast<int>(spec.size());
        open_symtab(&symtab_);
        init_t_atoms(&atoms_, numAtoms, TRUE);
        snew(x_, numAtoms);

        int  resind    = -1;
        int  prevResNr = -1;
        char prevChain = 0;
        for (int i = 0; i < numAtoms; ++i)
        {
            atoms_.atomname[i]    = put_symtab(&symtab_, spec[i].atomName);
            const bool newResidue = (spec[i].resNr != prevResNr || spec[i].chain != prevChain);
            if (newResidue)
            {
                ++resind;
            }
            atoms_.atom[i].resind = resind;
            if (newResidue)
            {
                t_atoms_set_resinfo(
                        &atoms_, i, &symtab_, spec[i].resName, spec[i].resNr, ' ', 0, spec[i].chain);
                prevResNr = spec[i].resNr;
                prevChain = spec[i].chain;
            }
            std::strncpy(atoms_.atom[i].elem, spec[i].element, sizeof(atoms_.atom[i].elem) - 1);
            atoms_.atom[i].elem[sizeof(atoms_.atom[i].elem) - 1] = '\0';
            atoms_.pdbinfo[i].type                               = PdbRecordType::Atom;
            atoms_.pdbinfo[i].occup                              = 1.0;
            atoms_.pdbinfo[i].bfac                               = 0.0;
            copy_rvec(spec[i].x, x_[i]);
        }
        atoms_.nres = resind + 1;
        close_symtab(&symtab_);

        clear_mat(box_);
        box_[XX][XX] = 3.0;
        box_[YY][YY] = 3.0;
        box_[ZZ][ZZ] = 3.0;
    }

    ~ConnectivityStructure()
    {
        done_atom(&atoms_);
        done_symtab(&symtab_);
        sfree(x_);
    }

    //! The atoms.
    t_atoms* atoms() { return &atoms_; }
    //! The coordinates.
    const rvec* x() const { return x_; }
    //! The number of atoms.
    int numAtoms() const { return atoms_.nr; }

    //! The (orthorhombic) box.
    matrix box_;

private:
    t_symtab symtab_;
    t_atoms  atoms_;
    rvec*    x_;
};

//! The whole content of \p fileName.
std::string readFile(const std::string& fileName)
{
    std::ifstream     stream(fileName);
    std::stringstream buffer;
    buffer << stream.rdbuf();
    return buffer.str();
}

//! A two-residue peptide, with its omega dihedral in the trans conformation.
std::vector<ConnectivityAtomSpec> peptideSpec()
{
    // CA(1)-C(1)-N(2)-CA(2) with the C-N bond along x, the two CA on opposite
    // sides of that bond.
    const real caNextY = -0.125;
    return { { "CA", "ALA", 1, 'A', "C", { -0.075, 0.125, 0.0 } },
             { "C", "ALA", 1, 'A', "C", { 0.0, 0.0, 0.0 } },
             { "N", "PRO", 2, 'A', "N", { 0.133, 0.0, 0.0 } },
             { "CA", "PRO", 2, 'A', "C", { 0.208, caNextY, 0.0 } } };
}

//! Two SG sulfurs that form a disulfide bridge, and a plain carbon bonded to one.
std::vector<ConnectivityAtomSpec> disulfideSpec()
{
    return { { "SG", "CYS", 1, 'A', "S", { 0.0, 0.0, 0.0 } },
             { "SG", "CYS", 2, 'A', "S", { 0.204, 0.0, 0.0 } },
             { "CB", "CYS", 1, 'A', "C", { -0.18, 0.0, 0.0 } } };
}

TEST(CifConnectivityTest, RoundTripPreservesBonds)
{
    TestFileManager       fileManager;
    ConnectivityStructure structure(peptideSpec());
    const std::string     cifFile = fileManager.getTemporaryFilePath("conn.cif").string();

    gmx_conect written = gmx_conect_init();
    gmx_conect_add(written, 0, 1);
    gmx_conect_add(written, 1, 2);
    gmx_conect_add(written, 2, 3);
    CifIO::write(cifFile, "conn", structure.atoms(), structure.x(), PbcType::Xyz, structure.box_, written);

    int natoms = 0;
    CifIO::countAtoms(cifFile, &natoms);
    ASSERT_EQ(natoms, structure.numAtoms());

    t_symtab symtab;
    open_symtab(&symtab);
    t_atoms atoms;
    init_t_atoms(&atoms, natoms, TRUE);
    rvec* x = nullptr;
    snew(x, natoms);
    matrix     box;
    PbcType    pbcType  = PbcType::Unset;
    gmx_conect readBack = gmx_conect_init();
    CifIO::read(cifFile, &symtab, nullptr, &atoms, x, &pbcType, box, readBack);

    EXPECT_EQ(gmx_conect_nconect(readBack), gmx_conect_nconect(written));
    EXPECT_TRUE(gmx_conect_exist(readBack, 0, 1));
    EXPECT_TRUE(gmx_conect_exist(readBack, 1, 2));
    EXPECT_TRUE(gmx_conect_exist(readBack, 2, 3));
    // A bond that was never written must not appear.
    EXPECT_FALSE(gmx_conect_exist(readBack, 0, 3));

    gmx_conect_done(readBack);
    gmx_conect_done(written);
    sfree(x);
    done_atom(&atoms);
    done_symtab(&symtab);
}

TEST(CifConnectivityTest, DisulfideBondIsTypedAsSuch)
{
    TestFileManager       fileManager;
    ConnectivityStructure structure(disulfideSpec());
    const std::string     cifFile = fileManager.getTemporaryFilePath("ss.cif").string();

    gmx_conect conect = gmx_conect_init();
    gmx_conect_add(conect, 0, 1); // SG-SG
    gmx_conect_add(conect, 0, 2); // SG-CB
    CifIO::write(cifFile, "ss", structure.atoms(), structure.x(), PbcType::Xyz, structure.box_, conect);

    const std::string content = readFile(cifFile);
    EXPECT_NE(content.find("disulf1"), std::string::npos) << content;
    EXPECT_NE(content.find("covale1"), std::string::npos) << content;
    // Both connection types must be declared.
    EXPECT_NE(content.find("_struct_conn_type.id"), std::string::npos) << content;

    gmx_conect_done(conect);
}

TEST(CifConnectivityTest, BondsToUnselectedAtomsAreDropped)
{
    TestFileManager       fileManager;
    ConnectivityStructure structure(peptideSpec());
    const std::string     cifFile = fileManager.getTemporaryFilePath("subset.cif").string();

    gmx_conect conect = gmx_conect_init();
    gmx_conect_add(conect, 0, 1);
    gmx_conect_add(conect, 1, 2);
    gmx_conect_add(conect, 2, 3); // atom 3 is left out below

    const int index[] = { 0, 1, 2 };
    CifIO::writeIndexed(
            cifFile, "subset", structure.atoms(), structure.x(), PbcType::Xyz, structure.box_, 3, index, conect);

    const std::string content = readFile(cifFile);
    EXPECT_NE(content.find("covale1"), std::string::npos) << content;
    EXPECT_NE(content.find("covale2"), std::string::npos) << content;
    EXPECT_EQ(content.find("covale3"), std::string::npos)
            << "a bond to an atom outside the index group was written\n"
            << content;

    gmx_conect_done(conect);
}

TEST(CifConnectivityTest, WithoutConnectivityNothingIsWritten)
{
    TestFileManager       fileManager;
    ConnectivityStructure structure(peptideSpec());
    const std::string     cifFile = fileManager.getTemporaryFilePath("plain.cif").string();

    CifIO::write(cifFile, "plain", structure.atoms(), structure.x(), PbcType::Xyz, structure.box_);

    const std::string content = readFile(cifFile);
    EXPECT_EQ(content.find("_struct_conn."), std::string::npos) << content;
}

TEST(CifAtomSiteTest, AlternateLocationAndEntityColumnsAreWritten)
{
    TestFileManager       fileManager;
    ConnectivityStructure structure(peptideSpec());
    const std::string     cifFile = fileManager.getTemporaryFilePath("altloc.cif").string();

    structure.atoms()->pdbinfo[1].altloc = 'B';
    CifIO::write(cifFile, "altloc", structure.atoms(), structure.x(), PbcType::Xyz, structure.box_);

    // Readers are entitled to both columns; gemmi reads no atoms without them.
    const std::string content = readFile(cifFile);
    EXPECT_NE(content.find("_atom_site.label_alt_id"), std::string::npos) << content;
    EXPECT_NE(content.find("_atom_site.label_entity_id"), std::string::npos) << content;

    // The alternate location survives the round trip; the reader has always
    // looked for label_alt_id, only the writer did not produce it.
    int natoms = 0;
    CifIO::countAtoms(cifFile, &natoms);
    ASSERT_EQ(natoms, structure.numAtoms());
    t_symtab symtab;
    open_symtab(&symtab);
    t_atoms atoms;
    init_t_atoms(&atoms, natoms, TRUE);
    rvec* x = nullptr;
    snew(x, natoms);
    matrix  box;
    PbcType pbcType = PbcType::Unset;
    CifIO::read(cifFile, &symtab, nullptr, &atoms, x, &pbcType, box);

    EXPECT_EQ(atoms.pdbinfo[1].altloc, 'B');
    EXPECT_EQ(atoms.pdbinfo[0].altloc, ' ');

    sfree(x);
    done_atom(&atoms);
    done_symtab(&symtab);
}

} // namespace
} // namespace test
} // namespace gmx
