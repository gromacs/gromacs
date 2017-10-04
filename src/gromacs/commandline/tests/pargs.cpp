/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
 * Tests parse_common_args().
 *
 * Currently, negative tests are not possible, because those trigger
 * gmx_fatal() and terminate the whole test binary.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include "gmxpre.h"

#include "gromacs/commandline/pargs.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace
{

using gmx::test::CommandLine;

class ParseCommonArgsTest : public ::testing::Test
{
    public:
        enum FileArgumentType
        {
            efFull,
            efNoExtension,
            efEmptyValue
        };

        ParseCommonArgsTest()
            : oenv_(nullptr), fileCount_(0)
        {
        }
        virtual ~ParseCommonArgsTest()
        {
            output_env_done(oenv_);
        }

        int nfile() const { return fileCount_; }

        void parseFromArgs(unsigned long           flags,
                           gmx::ArrayRef<t_filenm> fnm,
                           gmx::ArrayRef<t_pargs>  pa)
        {
            fileCount_ = fnm.size();
            bool bOk = parse_common_args(&args_.argc(), args_.argv(), flags,
                                         fnm.size(), fnm.data(),
                                         pa.size(), pa.data(),
                                         0, nullptr, 0, nullptr, &oenv_);
            EXPECT_TRUE(bOk);
        }
        void parseFromArray(gmx::ArrayRef<const char *const> cmdline,
                            unsigned long                    flags,
                            gmx::ArrayRef<t_filenm>          fnm,
                            gmx::ArrayRef<t_pargs>           pa)
        {
            args_.initFromArray(cmdline);
            parseFromArgs(flags, fnm, pa);
        }
        std::string addFileArg(const char *name, const char *extension,
                               FileArgumentType type)
        {
            std::string filename(tempFiles_.getTemporaryFilePath(extension));
            gmx::TextWriter::writeFileFromString(filename, "Dummy file");
            if (name != nullptr)
            {
                args_.append(name);
                switch (type)
                {
                    case efFull:
                        args_.append(filename);
                        break;
                    case efNoExtension:
                        args_.append(gmx::Path::stripExtension(filename));
                        break;
                    case efEmptyValue:
                        break;
                }
            }
            return filename;
        }

        // This must be a member that persists until the end of the test,
        // because string arguments are not duplicated in the output.
        CommandLine                 args_;

    private:
        gmx_output_env_t           *oenv_;
        size_t                      fileCount_;
        gmx::test::TestFileManager  tempFiles_;
};

/********************************************************************
 * Tests for different types of options
 */

TEST_F(ParseCommonArgsTest, ParsesIntegerArgs)
{
    int               value1 = 0, value2 = 0, value3 = 3;
    t_pargs           pa[]   = {
        { "-i1", FALSE, etINT, {&value1}, "Description" },
        { "-i2", FALSE, etINT, {&value2}, "Description" },
        { "-i3", FALSE, etINT, {&value3}, "Description" }
    };
    const char *const cmdline[] = {
        "test", "-i1", "2", "-i2", "-3"
    };
    parseFromArray(cmdline, 0, gmx::EmptyArrayRef(), pa);
    EXPECT_EQ( 2, value1);
    EXPECT_EQ(-3, value2);
    EXPECT_EQ( 3, value3);
}

TEST_F(ParseCommonArgsTest, ParsesInt64Args)
{
    gmx_int64_t       value1 = 0, value2 = 0, value3 = 3;
    t_pargs           pa[]   = {
        { "-i1", FALSE, etINT64, {&value1}, "Description" },
        { "-i2", FALSE, etINT64, {&value2}, "Description" },
        { "-i3", FALSE, etINT64, {&value3}, "Description" }
    };
    const char *const cmdline[] = {
        "test", "-i1", "2", "-i2", "-3"
    };
    parseFromArray(cmdline, 0, gmx::EmptyArrayRef(), pa);
    EXPECT_EQ( 2, value1);
    EXPECT_EQ(-3, value2);
    EXPECT_EQ( 3, value3);
}

TEST_F(ParseCommonArgsTest, ParsesRealArgs)
{
    real              value1 = 0.0, value2 = 0.0, value3 = 2.5;
    t_pargs           pa[]   = {
        { "-r1", FALSE, etREAL, {&value1}, "Description" },
        { "-r2", FALSE, etREAL, {&value2}, "Description" },
        { "-r3", FALSE, etREAL, {&value3}, "Description" }
    };
    const char *const cmdline[] = {
        "test", "-r1", "2", "-r2", "-.5"
    };
    parseFromArray(cmdline, 0, gmx::EmptyArrayRef(), pa);
    EXPECT_EQ( 2.0, value1);
    EXPECT_EQ(-0.5, value2);
    EXPECT_EQ( 2.5, value3);
}

TEST_F(ParseCommonArgsTest, ParsesStringArgs)
{
    const char       *value1 = "def", *value2 = "", *value3 = "default";
    t_pargs           pa[]   = {
        { "-s1", FALSE, etSTR, {&value1}, "Description" },
        { "-s2", FALSE, etSTR, {&value2}, "Description" },
        { "-s3", FALSE, etSTR, {&value3}, "Description" }
    };
    const char *const cmdline[] = {
        "test", "-s1", "", "-s2", "test"
    };
    parseFromArray(cmdline, 0, gmx::EmptyArrayRef(), pa);
    EXPECT_STREQ("", value1);
    EXPECT_STREQ("test", value2);
    EXPECT_STREQ("default", value3);
}

TEST_F(ParseCommonArgsTest, ParsesBooleanArgs)
{
    gmx_bool          value1 = TRUE, value2 = FALSE, value3 = TRUE;
    t_pargs           pa[]   = {
        { "-b1", FALSE, etBOOL, {&value1}, "Description" },
        { "-b2", FALSE, etBOOL, {&value2}, "Description" },
        { "-b3", FALSE, etBOOL, {&value3}, "Description" }
    };
    const char *const cmdline[] = {
        "test", "-nob1", "-b2"
    };
    parseFromArray(cmdline, 0, gmx::EmptyArrayRef(), pa);
    EXPECT_FALSE(value1);
    EXPECT_TRUE(value2);
    EXPECT_TRUE(value3);
}

TEST_F(ParseCommonArgsTest, ParsesVectorArgs)
{
    rvec              value1 = {0, 0, 0}, value2 = {0, 0, 0}, value3 = {1, 2, 3};
    t_pargs           pa[]   = {
        { "-v1", FALSE, etRVEC, {&value1}, "Description" },
        { "-v2", FALSE, etRVEC, {&value2}, "Description" },
        { "-v3", FALSE, etRVEC, {&value3}, "Description" }
    };
    const char *const cmdline[] = {
        "test", "-v1", "2", "1", "3", "-v2", "1"
    };
    parseFromArray(cmdline, 0, gmx::EmptyArrayRef(), pa);
    EXPECT_EQ(2.0, value1[XX]);
    EXPECT_EQ(1.0, value1[YY]);
    EXPECT_EQ(3.0, value1[ZZ]);
    EXPECT_EQ(1.0, value2[XX]);
    EXPECT_EQ(1.0, value2[YY]);
    EXPECT_EQ(1.0, value2[ZZ]);
    EXPECT_EQ(1.0, value3[XX]);
    EXPECT_EQ(2.0, value3[YY]);
    EXPECT_EQ(3.0, value3[ZZ]);
}

TEST_F(ParseCommonArgsTest, ParsesTimeArgs)
{
    real              value1 = 1.0, value2 = 2.0, value3 = 2.5;
    t_pargs           pa[]   = {
        { "-t1", FALSE, etTIME, {&value1}, "Description" },
        { "-t2", FALSE, etTIME, {&value2}, "Description" },
        { "-t3", FALSE, etTIME, {&value3}, "Description" }
    };
    const char *const cmdline[] = {
        "test", "-t1", "2", "-t2", "-.5"
    };
    parseFromArray(cmdline, 0, gmx::EmptyArrayRef(), pa);
    EXPECT_EQ( 2.0, value1);
    EXPECT_EQ(-0.5, value2);
    EXPECT_EQ( 2.5, value3);
}

TEST_F(ParseCommonArgsTest, ParsesTimeArgsWithTimeUnit)
{
    real              value1 = 1.0, value2 = 2.0, value3 = 2.5;
    t_pargs           pa[]   = {
        { "-t1", FALSE, etTIME, {&value1}, "Description" },
        { "-t2", FALSE, etTIME, {&value2}, "Description" },
        { "-t3", FALSE, etTIME, {&value3}, "Description" }
    };
    const char *const cmdline[] = {
        "test", "-t1", "2", "-t2", "-.5", "-tu", "ns"
    };
    parseFromArray(cmdline, PCA_TIME_UNIT, gmx::EmptyArrayRef(), pa);
    EXPECT_EQ( 2000.0, value1);
    EXPECT_EQ(-500.0, value2);
    EXPECT_EQ( 2.5, value3);
}

TEST_F(ParseCommonArgsTest, ParsesEnumArgs)
{
    const char       *value1[] = {nullptr, "none", "on", "off", nullptr };
    const char       *value2[] = {nullptr, "def", "value", "value_other", nullptr };
    const char       *value3[] = {nullptr, "default", "value", nullptr };
    t_pargs           pa[]     = {
        { "-s1", FALSE, etENUM, {value1}, "Description" },
        { "-s2", FALSE, etENUM, {value2}, "Description" },
        { "-s3", FALSE, etENUM, {value3}, "Description" }
    };
    const char *const cmdline[] = {
        "test", "-s1", "off", "-s2", "val"
    };
    parseFromArray(cmdline, 0, gmx::EmptyArrayRef(), pa);
    EXPECT_STREQ("off", value1[0]);
    EXPECT_STREQ("value", value2[0]);
    EXPECT_STREQ("default", value3[0]);
    EXPECT_EQ(value1[nenum(value1)], value1[0]);
    EXPECT_EQ(value2[nenum(value2)], value2[0]);
    EXPECT_EQ(value3[nenum(value3)], value3[0]);
}

/********************************************************************
 * Tests for file name options (output files, not dependent on file system)
 */

TEST_F(ParseCommonArgsTest, ParsesFileArgs)
{
    t_filenm          fnm[] = {
        { efXVG, "-o1", "out1", ffOPTWR },
        { efXVG, "-o2", "out2", ffOPTWR },
        { efXVG, "-om", "outm", ffWRMULT },
        { efXVG, "-om2", "outm2", ffWRMULT }
    };
    const char *const cmdline[] = {
        "test", "-o1", "-o2", "test", "-om", "test1", "test2.xvg", "-om2"
    };
    parseFromArray(cmdline, 0, fnm, gmx::EmptyArrayRef());
    EXPECT_STREQ("out1.xvg", opt2fn_null("-o1", nfile(), fnm));
    EXPECT_STREQ("test.xvg", opt2fn_null("-o2", nfile(), fnm));
    char **files;
    EXPECT_EQ(2, opt2fns(&files, "-om", nfile(), fnm));
    EXPECT_TRUE(files != nullptr);
    if (files != nullptr)
    {
        EXPECT_STREQ("test1.xvg", files[0]);
        EXPECT_STREQ("test2.xvg", files[1]);
    }
    EXPECT_STREQ("outm2.xvg", opt2fn("-om2", nfile(), fnm));
    done_filenms(nfile(), fnm);
}

TEST_F(ParseCommonArgsTest, ParsesFileArgsWithDefaults)
{
    t_filenm          fnm[] = {
        { efTPS, nullptr,  nullptr,   ffWRITE },
        { efTRX, "-f2", nullptr,   ffOPTWR },
        { efTRX, "-f3", "trj3", ffWRITE },
        { efXVG, "-o",  "out",  ffWRITE },
        { efXVG, "-om", "outm", ffWRMULT }
    };
    const char *const cmdline[] = {
        "test"
    };
    parseFromArray(cmdline, 0, fnm, gmx::EmptyArrayRef());
    EXPECT_STREQ("topol.tpr", ftp2fn(efTPS, nfile(), fnm));
    EXPECT_STREQ("traj.xtc", opt2fn("-f2", nfile(), fnm));
    EXPECT_EQ(nullptr, opt2fn_null("-f2", nfile(), fnm));
    EXPECT_STREQ("trj3.xtc", opt2fn("-f3", nfile(), fnm));
    EXPECT_STREQ("out.xvg", opt2fn("-o", nfile(), fnm));
    EXPECT_STREQ("outm.xvg", opt2fn("-om", nfile(), fnm));
    done_filenms(nfile(), fnm);
}

TEST_F(ParseCommonArgsTest, ParsesFileArgsWithDefaultFileName)
{
    t_filenm          fnm[] = {
        { efTPS, "-s",  nullptr,   ffWRITE },
        { efTRX, "-f2", nullptr,   ffWRITE },
        { efTRX, "-f3", "trj3", ffWRITE },
        { efXVG, "-o",  "out",  ffWRITE },
        { efXVG, "-om", "outm", ffWRMULT }
    };
    const char *const cmdline[] = {
        "test", "-deffnm", "def", "-f2", "other", "-o"
    };
    parseFromArray(cmdline, PCA_CAN_SET_DEFFNM, fnm, gmx::EmptyArrayRef());
    EXPECT_STREQ("def.tpr", ftp2fn(efTPS, nfile(), fnm));
    EXPECT_STREQ("other.xtc", opt2fn("-f2", nfile(), fnm));
    EXPECT_STREQ("def.xtc", opt2fn("-f3", nfile(), fnm));
    EXPECT_STREQ("def.xvg", opt2fn("-o", nfile(), fnm));
    EXPECT_STREQ("def.xvg", opt2fn("-om", nfile(), fnm));
    done_filenms(nfile(), fnm);
}

TEST_F(ParseCommonArgsTest, ParseFileArgsWithCustomDefaultExtension)
{
    t_filenm          fnm[] = {
        { efTRX, "-o1", "conf1.gro", ffWRITE },
        { efTRX, "-o2", "conf2.pdb", ffWRITE },
        { efTRX, "-o3", "conf3.gro", ffWRITE }
    };
    const char *const cmdline[] = {
        "test", "-o2", "-o3", "test"
    };
    parseFromArray(cmdline, PCA_CAN_SET_DEFFNM, fnm, gmx::EmptyArrayRef());
    EXPECT_STREQ("conf1.gro", opt2fn("-o1", nfile(), fnm));
    EXPECT_STREQ("conf2.pdb", opt2fn("-o2", nfile(), fnm));
    EXPECT_STREQ("test.gro", opt2fn("-o3", nfile(), fnm));
    done_filenms(nfile(), fnm);
}

/********************************************************************
 * Tests for file name options (input files, dependent on file system contents)
 */

TEST_F(ParseCommonArgsTest, HandlesNonExistentInputFiles)
{
    t_filenm          fnm[] = {
        { efTPS, "-s",  nullptr,   ffREAD },
        { efTRX, "-f",  "trj",  ffREAD },
        { efTRX, "-f2", "trj2", ffREAD },
        { efTRX, "-f3", "trj3", ffREAD },
        { efTRX, "-f4", "trj4", ffREAD },
        { efGRO, "-g",  "cnf",  ffREAD },
        { efGRO, "-g2", "cnf2", ffREAD }
    };
    const char *const cmdline[] = {
        "test", "-f2", "-f3", "other", "-f4", "trj.gro", "-g2", "foo"
    };
    parseFromArray(cmdline, PCA_DISABLE_INPUT_FILE_CHECKING, fnm, gmx::EmptyArrayRef());
    EXPECT_STREQ("topol.tpr", ftp2fn(efTPS, nfile(), fnm));
    EXPECT_STREQ("trj.xtc", opt2fn("-f", nfile(), fnm));
    EXPECT_STREQ("trj2.xtc", opt2fn("-f2", nfile(), fnm));
    EXPECT_STREQ("other.xtc", opt2fn("-f3", nfile(), fnm));
    EXPECT_STREQ("trj.gro", opt2fn("-f4", nfile(), fnm));
    EXPECT_STREQ("cnf.gro", opt2fn("-g", nfile(), fnm));
    EXPECT_STREQ("foo.gro", opt2fn("-g2", nfile(), fnm));
    done_filenms(nfile(), fnm);
}

TEST_F(ParseCommonArgsTest, HandlesNonExistentOptionalInputFiles)
{
    t_filenm          fnm[] = {
        { efTPS, "-s",  nullptr,   ffOPTRD },
        { efTRX, "-f",  "trj",  ffOPTRD }
    };
    const char *const cmdline[] = {
        "test"
    };
    parseFromArray(cmdline, 0, fnm, gmx::EmptyArrayRef());
    EXPECT_STREQ("topol.tpr", ftp2fn(efTPS, nfile(), fnm));
    EXPECT_STREQ("trj.xtc", opt2fn("-f", nfile(), fnm));
    done_filenms(nfile(), fnm);
}

TEST_F(ParseCommonArgsTest, AcceptsNonExistentInputFilesIfSpecified)
{
    t_filenm          fnm[] = {
        { efCPT, "-c",  "file1", ffOPTRD | ffALLOW_MISSING },
        { efCPT, "-c2", "file2", ffOPTRD | ffALLOW_MISSING },
        { efCPT, "-c3", "file3", ffOPTRD | ffALLOW_MISSING },
        { efCPT, "-c4", "file4", ffOPTRD | ffALLOW_MISSING },
        { efTRX, "-f",  "trj",   ffOPTRD | ffALLOW_MISSING }
    };
    const char *const cmdline[] = {
        "test", "-c2", "-c3", "nonexistent", "-c4", "nonexistent.cpt", "-f", "nonexistent"
    };
    parseFromArray(cmdline, 0, fnm, gmx::EmptyArrayRef());
    EXPECT_STREQ("file1.cpt", opt2fn("-c", nfile(), fnm));
    EXPECT_STREQ("file2.cpt", opt2fn("-c2", nfile(), fnm));
    EXPECT_STREQ("nonexistent.cpt", opt2fn("-c3", nfile(), fnm));
    EXPECT_STREQ("nonexistent.cpt", opt2fn("-c4", nfile(), fnm));
    EXPECT_STREQ("nonexistent.xtc", opt2fn("-f", nfile(), fnm));
    done_filenms(nfile(), fnm);
}

TEST_F(ParseCommonArgsTest, HandlesCompressedFiles)
{
    t_filenm          fnm[] = {
        { efTRX, "-f", nullptr, ffREAD },
        { efGRO, "-g", nullptr, ffREAD }
    };
    args_.append("test");
    std::string       expectedF = addFileArg("-f", ".pdb.gz", efFull);
    std::string       expectedG = addFileArg("-g", ".gro.Z", efFull);
    expectedF = gmx::Path::stripExtension(expectedF);
    expectedG = gmx::Path::stripExtension(expectedG);
    parseFromArgs(0, fnm, gmx::EmptyArrayRef());
    EXPECT_EQ(expectedF, opt2fn("-f", nfile(), fnm));
    EXPECT_EQ(expectedG, opt2fn("-g", nfile(), fnm));
    done_filenms(nfile(), fnm);
}

TEST_F(ParseCommonArgsTest, AcceptsUnknownTrajectoryExtension)
{
    t_filenm          fnm[] = {
        { efTRX, "-f", nullptr, ffREAD }
    };
    args_.append("test");
    std::string       expected = addFileArg("-f", ".foo", efFull);
    parseFromArgs(0, fnm, gmx::EmptyArrayRef());
    EXPECT_EQ(expected, opt2fn("-f", nfile(), fnm));
    done_filenms(nfile(), fnm);
}

TEST_F(ParseCommonArgsTest, CompletesExtensionFromExistingFile)
{
    t_filenm          fnm[] = {
        { efTRX, "-f1", nullptr, ffREAD },
        { efTRX, "-f2", nullptr, ffREAD },
        { efTRX, "-f3", nullptr, ffREAD },
        { efTRX, "-f4", nullptr, ffREAD }
    };
    args_.append("test");
    std::string       expected1 = addFileArg("-f1", "1.xtc", efNoExtension);
    std::string       expected2 = addFileArg("-f2", "2.gro", efNoExtension);
    std::string       expected3 = addFileArg("-f3", "3.tng", efNoExtension);
    std::string       expected4 = addFileArg("-f4", ".gro", efEmptyValue);
    std::string       def4      = gmx::Path::stripExtension(expected4);
    fnm[3].fn = def4.c_str();
    parseFromArgs(0, fnm, gmx::EmptyArrayRef());
    EXPECT_EQ(expected1, opt2fn("-f1", nfile(), fnm));
    EXPECT_EQ(expected2, opt2fn("-f2", nfile(), fnm));
    EXPECT_EQ(expected3, opt2fn("-f3", nfile(), fnm));
    EXPECT_EQ(expected4, opt2fn("-f4", nfile(), fnm));
    done_filenms(nfile(), fnm);
}

TEST_F(ParseCommonArgsTest, CompletesExtensionFromExistingFileWithDefaultFileName)
{
    t_filenm          fnm[] = {
        { efTRX, "-f1", nullptr,  ffREAD },
        { efSTO, "-f2", "foo", ffREAD },
        { efTRX, "-f3", nullptr,  ffREAD },
        { efSTX, "-f4", nullptr,  ffREAD }
    };
    args_.append("test");
    std::string       expected1 = addFileArg("-f1", "1.trr", efNoExtension);
    std::string       expected2 = addFileArg("-f2", ".pdb", efEmptyValue);
    std::string       expected3 = addFileArg("-f3", ".trr", efEmptyValue);
    std::string       expected4 = addFileArg(nullptr, ".pdb", efEmptyValue);
    std::string       deffnm    = gmx::Path::stripExtension(expected3);
    args_.append("-deffnm");
    args_.append(deffnm);
    parseFromArgs(PCA_CAN_SET_DEFFNM, fnm, gmx::EmptyArrayRef());
    EXPECT_EQ(expected1, opt2fn("-f1", nfile(), fnm));
    EXPECT_EQ(expected2, opt2fn("-f2", nfile(), fnm));
    EXPECT_EQ(expected3, opt2fn("-f3", nfile(), fnm));
    EXPECT_EQ(expected4, opt2fn("-f4", nfile(), fnm));
    done_filenms(nfile(), fnm);
}

/********************************************************************
 * Tests for general behavior
 */

TEST_F(ParseCommonArgsTest, HandlesNonReadNode)
{
    t_filenm          fnm[] = {
        { efTPS, "-s",  nullptr,  ffREAD },
        { efTRX, "-f",  nullptr,  ffREAD },
        { efTRX, "-f2", nullptr,  ffREAD }
    };
    const char *const cmdline[] = {
        "test", "-f", "-f2", "other"
    };
    parseFromArray(cmdline, PCA_NOT_READ_NODE, fnm, gmx::EmptyArrayRef());
    EXPECT_EQ(nullptr, fnm[0].fns);
    EXPECT_EQ(nullptr, fnm[1].fns);
    EXPECT_EQ(nullptr, fnm[2].fns);
    done_filenms(nfile(), fnm);
}

TEST_F(ParseCommonArgsTest, HandlesNonReadNodeWithDefaultFileName)
{
    t_filenm          fnm[] = {
        { efTPS, "-s",  nullptr,  ffREAD },
        { efTRX, "-f",  nullptr,  ffREAD },
        { efTRX, "-f2", nullptr,  ffREAD }
    };
    const char *const cmdline[] = {
        "test", "-deffnm", "def", "-f", "-f2", "other"
    };
    parseFromArray(cmdline, PCA_CAN_SET_DEFFNM | PCA_NOT_READ_NODE, fnm, gmx::EmptyArrayRef());
    EXPECT_EQ(nullptr, fnm[0].fns);
    EXPECT_EQ(nullptr, fnm[1].fns);
    EXPECT_EQ(nullptr, fnm[2].fns);
    done_filenms(nfile(), fnm);
}

TEST_F(ParseCommonArgsTest, CanKeepUnknownArgs)
{
    int               ivalue = 0;
    gmx_bool          bvalue = FALSE;
    t_pargs           pa[]   = {
        { "-i", FALSE, etINT, {&ivalue}, "Description" },
        { "-b", FALSE, etBOOL, {&bvalue}, "Description" },
    };
    t_filenm          fnm[] = {
        { efXVG, "-o1", "out1", ffOPTWR },
        { efXVG, "-o2", "out2", ffOPTWR }
    };
    const char *const cmdline[] = {
        "test", "foo", "-unk", "-o1", "-unk2", "-o2", "test",
        "-i", "2", "-unk3", "-b", "-unk4"
    };
    parseFromArray(cmdline, PCA_NOEXIT_ON_ARGS, fnm, pa);
    EXPECT_EQ(2, ivalue);
    EXPECT_TRUE(bvalue);
    EXPECT_STREQ("out1.xvg", opt2fn_null("-o1", nfile(), fnm));
    EXPECT_STREQ("test.xvg", opt2fn_null("-o2", nfile(), fnm));
    EXPECT_EQ(6, args_.argc());
    EXPECT_STREQ(cmdline[0],  args_.arg(0));
    EXPECT_STREQ(cmdline[1],  args_.arg(1));
    EXPECT_STREQ(cmdline[2],  args_.arg(2));
    EXPECT_STREQ(cmdline[4],  args_.arg(3));
    EXPECT_STREQ(cmdline[9],  args_.arg(4));
    EXPECT_STREQ(cmdline[11], args_.arg(5));
    done_filenms(nfile(), fnm);
}

} // namespace
