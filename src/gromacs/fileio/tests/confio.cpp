/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
 * Tests for reading/writing different structure file formats.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "gromacs/fileio/confio.h"

#include <iostream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fileio/filetypes.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/stringtest.h"
#include "testutils/testfilemanager.h"

// TODO: These should really appear somewhere centralized.
/*! \brief
 * Google Test formatter for GromacsFileType values.
 */
static void PrintTo(const GromacsFileType& ftp, std::ostream* os)
{
    *os << "'" << ftp2ext(ftp) << "'";
}

namespace
{

class StructureIORoundtripTest :
    public gmx::test::StringTestBase,
    public ::testing::WithParamInterface<GromacsFileType>
{
public:
    StructureIORoundtripTest()
    {
        generateReferenceTopology();
        generateReferenceCoordinates();
        testTop_ = nullptr;
        testX_   = nullptr;
        clear_mat(testBox_);
        referenceFilename_ = fileManager_.getTemporaryFilePath(getFileSuffix("ref"));
        testFilename_      = fileManager_.getTemporaryFilePath(getFileSuffix("test"));
    }
    ~StructureIORoundtripTest() override
    {
        if (testTop_ != nullptr)
        {
            done_top(testTop_);
            sfree(testTop_);
        }
        sfree(testX_);
        done_top(refTop_);
        sfree(refTop_);
    }

    void writeReferenceFile()
    {
        write_sto_conf(referenceFilename_.c_str(), *refTop_->name, &refTop_->atoms,
                       as_rvec_array(refX_.data()), nullptr, -1, refBox_);
    }

    void readReferenceFileTps()
    {
        snew(testTop_, 1);
        int ePBC = -2;
        read_tps_conf(referenceFilename_.c_str(), testTop_, &ePBC, &testX_, nullptr, testBox_, FALSE);
    }

    void testTopologies()
    {
        // TODO: Compare the topologies.
    }

    void writeTestFileAndTest()
    {
        write_sto_conf(testFilename_.c_str(), *testTop_->name, &testTop_->atoms, testX_, nullptr,
                       -1, testBox_);
        testFilesEqual(referenceFilename_, testFilename_);
    }

private:
    std::string getFileSuffix(const char* type)
    {
        return std::string(type) + "." + ftp2ext(GetParam());
    }

    void generateReferenceTopology()
    {
        snew(refTop_, 1);
        open_symtab(&refTop_->symtab);
        if (GetParam() == efESP)
        {
            // Titles cannot be read from an .esp file...
            refTop_->name = put_symtab(&refTop_->symtab, "");
        }
        else
        {
            refTop_->name = put_symtab(&refTop_->symtab, "Test title");
        }
        const int atomCount = 10;
        init_t_atoms(&refTop_->atoms, atomCount, FALSE);
        for (int i = 0; i < atomCount; ++i)
        {
            char name[3];
            name[0]                       = 'A';
            name[1]                       = 'A' + i % 3;
            name[2]                       = '\0';
            refTop_->atoms.atomname[i]    = put_symtab(&refTop_->symtab, name);
            refTop_->atoms.atom[i].resind = i / 3;
            if (i % 3 == 0)
            {
                char resname[3];
                resname[0] = 'R';
                resname[1] = 'A' + i / 3;
                resname[2] = '\0';
                t_atoms_set_resinfo(&refTop_->atoms, i, &refTop_->symtab, resname, i / 3 + 1, ' ', 0, ' ');
            }
        }
        refTop_->atoms.nres = 4;
        close_symtab(&refTop_->symtab);
    }

    void generateReferenceCoordinates()
    {
        clear_mat(refBox_);
        refBox_[XX][XX]     = 1;
        refBox_[YY][YY]     = 2;
        refBox_[ZZ][ZZ]     = 3;
        const int atomCount = refTop_->atoms.nr;
        refX_.reserve(atomCount);
        for (int i = 0; i < atomCount; ++i)
        {
            refX_.emplace_back(i % 4, i / 4, (i / 2) % 3);
        }
    }

    gmx::test::TestFileManager fileManager_;
    std::string                referenceFilename_;
    std::string                testFilename_;
    t_topology*                refTop_;
    std::vector<gmx::RVec>     refX_;
    matrix                     refBox_;
    t_topology*                testTop_;
    rvec*                      testX_;
    matrix                     testBox_;
};

TEST_P(StructureIORoundtripTest, ReadWriteTpsConf)
{
    writeReferenceFile();
    readReferenceFileTps();
    testTopologies();
    writeTestFileAndTest();
}

INSTANTIATE_TEST_CASE_P(WithDifferentFormats,
                        StructureIORoundtripTest,
                        ::testing::Values(efGRO, efG96, efPDB, efESP));

} // namespace
