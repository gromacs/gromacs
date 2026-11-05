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
 * Tests for multi-model mmCIF trajectory reading and writing (cifio).
 *
 * Each mmCIF model (pdbx_PDB_model_num) is treated as a trajectory frame; these
 * check that models written frame-by-frame are recovered as separate frames
 * with the correct per-model coordinates.
 *
 * \author Amr Alhossary
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include <cstdio>

#include <array>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fileio/cifio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/textreader.h"
#include "gromacs/utility/vec.h"

#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

class CifTrajectoryIoTest : public ::testing::Test
{
protected:
    static constexpr int c_numAtoms = 3;

    CifTrajectoryIoTest()
    {
        // Minimal single-residue t_atoms for the writer.
        open_symtab(&symtab_);
        init_t_atoms(&atoms_, c_numAtoms, FALSE);
        const char* const names[c_numAtoms] = { "N", "CA", "C" };
        for (int i = 0; i < c_numAtoms; ++i)
        {
            atoms_.atomname[i]    = put_symtab(&symtab_, names[i]);
            atoms_.atom[i].resind = 0;
        }
        t_atoms_set_resinfo(&atoms_, 0, &symtab_, "ALA", 1, ' ', 0, 'A');
        atoms_.nres = 1;
        close_symtab(&symtab_);

        clear_mat(box_);
        box_[XX][XX] = 3.0;
        box_[YY][YY] = 3.0;
        box_[ZZ][ZZ] = 3.0;
    }

public:
    ~CifTrajectoryIoTest() override
    {
        done_atom(&atoms_);
        done_symtab(&symtab_);
    }

protected:
    //! Coordinates for model \p m (models differ by a per-model x shift).
    static std::array<RVec, c_numAtoms> modelCoords(int m)
    {
        std::array<RVec, c_numAtoms> x;
        for (int i = 0; i < c_numAtoms; ++i)
        {
            x[i] = { 1.0F + (0.1F * i) + static_cast<real>(m), 1.0F, 1.0F };
        }
        return x;
    }

    TestFileManager fileManager_;
    t_symtab        symtab_;
    t_atoms         atoms_;
    matrix          box_;
};

TEST_F(CifTrajectoryIoTest, WriteReadTwoModelsRoundTrips)
{
    const std::string cifFile = fileManager_.getTemporaryFilePath("traj.cif").string();
    const std::array<int, c_numAtoms> index     = { 0, 1, 2 };
    const int                         numModels = 2;

    // Stream two models to the file; the header is written with the first only.
    FILE* out = std::fopen(cifFile.c_str(), "w");
    ASSERT_NE(out, nullptr);
    for (int m = 0; m < numModels; ++m)
    {
        const auto x = modelCoords(m);
        CifIO::writeTrajectoryFrame(out,
                                    "traj",
                                    &atoms_,
                                    as_rvec_array(x.data()),
                                    box_,
                                    c_numAtoms,
                                    index.data(),
                                    /*modelNr*/ m + 1,
                                    /*writeHeader*/ m == 0);
    }
    std::fclose(out);

    // Read the models back as frames.
    const std::unique_ptr<CifTrajectory> trajectory = CifIO::openTrajectory(cifFile);
    ASSERT_NE(trajectory, nullptr);
    EXPECT_EQ(trajectory->numAtoms(), c_numAtoms);

    t_trxframe fr;
    clear_trxframe(&fr, TRUE);
    snew(fr.x, c_numAtoms);

    for (int m = 0; m < numModels; ++m)
    {
        SCOPED_TRACE("model " + std::to_string(m));
        ASSERT_TRUE(trajectory->readNextFrame(&fr));
        const auto expected = modelCoords(m);
        for (int i = 0; i < c_numAtoms; ++i)
        {
            EXPECT_NEAR(fr.x[i][XX], expected[i][XX], 1e-3);
            EXPECT_NEAR(fr.x[i][YY], expected[i][YY], 1e-3);
            EXPECT_NEAR(fr.x[i][ZZ], expected[i][ZZ], 1e-3);
        }
        EXPECT_TRUE(fr.bBox);
        EXPECT_NEAR(fr.box[XX][XX], box_[XX][XX], 1e-3);
    }
    // No more models.
    EXPECT_FALSE(trajectory->readNextFrame(&fr));

    sfree(fr.x);
}

TEST_F(CifTrajectoryIoTest, ConnectivityIsWrittenOnceForAllModels)
{
    const std::string cifFile = fileManager_.getTemporaryFilePath("conn_traj.cif").string();
    const std::array<int, c_numAtoms> index     = { 0, 1, 2 };
    const int                         numModels = 3;

    gmx_conect conect = gmx_conect_init();
    gmx_conect_add(conect, 0, 1);
    gmx_conect_add(conect, 1, 2);

    FILE* out = std::fopen(cifFile.c_str(), "w");
    ASSERT_NE(out, nullptr);
    for (int m = 0; m < numModels; ++m)
    {
        const auto x = modelCoords(m);
        CifIO::writeTrajectoryFrame(out,
                                    "traj",
                                    &atoms_,
                                    as_rvec_array(x.data()),
                                    box_,
                                    c_numAtoms,
                                    index.data(),
                                    /*modelNr*/ m + 1,
                                    /*writeHeader*/ m == 0,
                                    conect);
    }
    std::fclose(out);
    gmx_conect_done(conect);

    // The connectivity does not change between models, so it is written with the
    // header only; the atom_site loop that follows stays open for all models.
    const std::string content = TextReader::readFileToString(cifFile);

    int numLoopHeaders = 0;
    for (size_t pos = content.find("_struct_conn.id"); pos != std::string::npos;
         pos        = content.find("_struct_conn.id", pos + 1))
    {
        numLoopHeaders++;
    }
    EXPECT_EQ(numLoopHeaders, 1) << content;
    EXPECT_NE(content.find("covale1"), std::string::npos) << content;
    EXPECT_NE(content.find("covale2"), std::string::npos) << content;

    // The models are still readable afterwards.
    const std::unique_ptr<CifTrajectory> trajectory = CifIO::openTrajectory(cifFile);
    ASSERT_NE(trajectory, nullptr);
    EXPECT_EQ(trajectory->numAtoms(), c_numAtoms);
    t_trxframe fr;
    clear_trxframe(&fr, TRUE);
    snew(fr.x, c_numAtoms);
    for (int m = 0; m < numModels; ++m)
    {
        SCOPED_TRACE("model " + std::to_string(m));
        ASSERT_TRUE(trajectory->readNextFrame(&fr));
        const auto expected = modelCoords(m);
        EXPECT_NEAR(fr.x[0][XX], expected[0][XX], 1e-3);
    }
    EXPECT_FALSE(trajectory->readNextFrame(&fr));
    sfree(fr.x);
}

TEST_F(CifTrajectoryIoTest, TitleWithSpacesIsReadBack)
{
    const std::string cifFile = fileManager_.getTemporaryFilePath("title.cif").string();
    const std::array<int, c_numAtoms> index        = { 0, 1, 2 };
    const auto                        x            = modelCoords(0);
    const char* const                 writtenTitle = "a title with spaces";

    FILE* out = std::fopen(cifFile.c_str(), "w");
    ASSERT_NE(out, nullptr);
    CifIO::writeTrajectoryFrame(out,
                                writtenTitle,
                                &atoms_,
                                as_rvec_array(x.data()),
                                box_,
                                c_numAtoms,
                                index.data(),
                                /*modelNr*/ 1,
                                /*writeHeader*/ true);
    std::fclose(out);

    t_symtab symtab;
    open_symtab(&symtab);
    t_atoms atoms;
    init_t_atoms(&atoms, c_numAtoms, FALSE);
    std::array<RVec, c_numAtoms> readX;
    char*                        title = nullptr;
    matrix                       box;
    PbcType                      pbcType = PbcType::Unset;
    CifIO::read(cifFile, &symtab, &title, &atoms, as_rvec_array(readX.data()), &pbcType, box);
    EXPECT_STREQ(title, writtenTitle);
    sfree(title);
    done_atom(&atoms);
    done_symtab(&symtab);
}

TEST_F(CifTrajectoryIoTest, EmptyFileYieldsNoTrajectory)
{
    // A data block with no atom_site table.
    const std::unique_ptr<CifTrajectory> trajectory =
            CifIO::openTrajectory(TestFileManager::getInputFilePath("empty.cif"));
    EXPECT_EQ(trajectory, nullptr);
}

} // namespace
} // namespace test
} // namespace gmx
