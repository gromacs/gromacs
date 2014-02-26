/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
 * Implements classes in moduletest.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "moduletest.h"

#include <map>
#include <string>
#include <vector>

#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/cmdlinerunner.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/stringutil.h"

#include "gromacs/analysisdata/tests/datatest.h"
#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{

/********************************************************************
 * AbstractTrajectoryAnalysisModuleTestFixture::Impl
 */

class AbstractTrajectoryAnalysisModuleTestFixture::Impl
{
    public:
        struct DatasetInfo
        {
            DatasetInfo()
                : bCheck(true), tolerance(defaultRealTolerance())
            {
            }

            bool                   bCheck;
            FloatingPointTolerance tolerance;
        };

        struct OutputFileInfo
        {
            OutputFileInfo(const char *option, const std::string &path)
                : option(option), path(path)
            {
            }

            std::string         option;
            std::string         path;
        };

        typedef std::map<std::string, DatasetInfo> DatasetList;
        typedef std::vector<OutputFileInfo>        OutputFileList;

        explicit Impl(AbstractTrajectoryAnalysisModuleTestFixture *parent);

        TrajectoryAnalysisModule &module();
        void ensureModuleCreated();
        bool hasCheckedDatasets() const;

        AbstractTrajectoryAnalysisModuleTestFixture    &parent_;
        TrajectoryAnalysisModulePointer                 module_;
        TestReferenceData                               data_;
        CommandLine                                     cmdline_;
        TestFileManager                                 tempFiles_;
        DatasetList                                     datasets_;
        OutputFileList                                  outputFiles_;
        bool                                            bDatasetsIncluded_;
};

AbstractTrajectoryAnalysisModuleTestFixture::Impl::Impl(
        AbstractTrajectoryAnalysisModuleTestFixture *parent)
    : parent_(*parent), bDatasetsIncluded_(false)
{
    cmdline_.append("module");
}

TrajectoryAnalysisModule &
AbstractTrajectoryAnalysisModuleTestFixture::Impl::module()
{
    ensureModuleCreated();
    return *module_;
}

void
AbstractTrajectoryAnalysisModuleTestFixture::Impl::ensureModuleCreated()
{
    if (module_.get() == NULL)
    {
        module_ = parent_.createModule();
        const std::vector<std::string>          &datasetNames(module_->datasetNames());
        datasets_.clear();
        std::vector<std::string>::const_iterator i;
        for (i = datasetNames.begin(); i != datasetNames.end(); ++i)
        {
            datasets_[*i] = DatasetInfo();
        }
    }
}

bool
AbstractTrajectoryAnalysisModuleTestFixture::Impl::hasCheckedDatasets() const
{
    DatasetList::const_iterator dataset;
    for (dataset = datasets_.begin(); dataset != datasets_.end(); ++dataset)
    {
        if (dataset->second.bCheck)
        {
            return true;
        }
    }
    return false;
}

/********************************************************************
 * AbstractTrajectoryAnalysisModuleTestFixture
 */

AbstractTrajectoryAnalysisModuleTestFixture::AbstractTrajectoryAnalysisModuleTestFixture()
    : impl_(new Impl(this))
{
}

AbstractTrajectoryAnalysisModuleTestFixture::~AbstractTrajectoryAnalysisModuleTestFixture()
{
}

void
AbstractTrajectoryAnalysisModuleTestFixture::setTopology(const char *filename)
{
    impl_->cmdline_.append("-s");
    impl_->cmdline_.append(TestFileManager::getInputFilePath(filename));
}

void
AbstractTrajectoryAnalysisModuleTestFixture::setTrajectory(const char *filename)
{
    impl_->cmdline_.append("-f");
    impl_->cmdline_.append(TestFileManager::getInputFilePath(filename));
}

void
AbstractTrajectoryAnalysisModuleTestFixture::setOutputFile(const char *option,
                                                           const char *filename)
{
    std::string fullFilename = impl_->tempFiles_.getTemporaryFilePath(filename);
    impl_->cmdline_.append(option);
    impl_->cmdline_.append(fullFilename);
    impl_->outputFiles_.push_back(Impl::OutputFileInfo(option, fullFilename));
}

void
AbstractTrajectoryAnalysisModuleTestFixture::setOutputFileNoTest(
        const char *option, const char *extension)
{
    std::string fullFilename = impl_->tempFiles_.getTemporaryFilePath(
                formatString("%d.%s", impl_->cmdline_.argc(), extension));
    impl_->cmdline_.append(option);
    impl_->cmdline_.append(fullFilename);
}

void
AbstractTrajectoryAnalysisModuleTestFixture::includeDataset(const char *name)
{
    impl_->ensureModuleCreated();
    if (!impl_->bDatasetsIncluded_)
    {
        Impl::DatasetList::iterator i;
        for (i = impl_->datasets_.begin(); i != impl_->datasets_.end(); ++i)
        {
            i->second.bCheck = false;
        }
    }
    Impl::DatasetList::iterator dataset = impl_->datasets_.find(name);
    const bool                  bFound  = (dataset != impl_->datasets_.end());
    GMX_RELEASE_ASSERT(bFound, "Attempted to include a non-existent dataset");
    dataset->second.bCheck = true;
}

void
AbstractTrajectoryAnalysisModuleTestFixture::excludeDataset(const char *name)
{
    impl_->ensureModuleCreated();
    Impl::DatasetList::iterator dataset = impl_->datasets_.find(name);
    const bool                  bFound  = (dataset != impl_->datasets_.end());
    GMX_RELEASE_ASSERT(bFound, "Attempted to exclude a non-existent dataset");
    dataset->second.bCheck = false;
}

void
AbstractTrajectoryAnalysisModuleTestFixture::setDatasetTolerance(
        const char *name, const FloatingPointTolerance &tolerance)
{
    impl_->ensureModuleCreated();
    Impl::DatasetList::iterator dataset = impl_->datasets_.find(name);
    const bool                  bFound  = (dataset != impl_->datasets_.end());
    GMX_RELEASE_ASSERT(bFound, "Attempted to set a tolerance for a non-existent dataset");
    dataset->second.tolerance = tolerance;
}

void
AbstractTrajectoryAnalysisModuleTestFixture::runTest(const CommandLine &args)
{
    TrajectoryAnalysisModule &module = impl_->module();
    // Skip first argument if it is the module name.
    int firstArg = (args.arg(0)[0] == '-' ? 0 : 1);
    for (int i = firstArg; i < args.argc(); ++i)
    {
        impl_->cmdline_.append(args.arg(i));
    }

    TestReferenceChecker rootChecker(impl_->data_.rootChecker());

    rootChecker.checkString(args.toString(), "CommandLine");

    if (impl_->hasCheckedDatasets())
    {
        TestReferenceChecker               dataChecker(
                rootChecker.checkCompound("OutputData", "Data"));
        Impl::DatasetList::const_iterator  dataset;
        for (dataset = impl_->datasets_.begin();
             dataset != impl_->datasets_.end();
             ++dataset)
        {
            if (dataset->second.bCheck)
            {
                const char *const     name = dataset->first.c_str();
                AbstractAnalysisData &data = module.datasetFromName(name);
                AnalysisDataTestFixture::addReferenceCheckerModule(
                        dataChecker, name, &data, dataset->second.tolerance);
            }
        }
    }

    TrajectoryAnalysisCommandLineRunner runner(&module);
    runner.setUseDefaultGroups(false);
    int rc = 0;
    EXPECT_NO_THROW_GMX(rc = runner.run(impl_->cmdline_.argc(), impl_->cmdline_.argv()));
    EXPECT_EQ(0, rc);

    if (!impl_->outputFiles_.empty())
    {
        TestReferenceChecker                 outputChecker(
                rootChecker.checkCompound("OutputFiles", "Files"));
        Impl::OutputFileList::const_iterator outfile;
        for (outfile = impl_->outputFiles_.begin();
             outfile != impl_->outputFiles_.end();
             ++outfile)
        {
            std::string output = File::readToString(outfile->path);
            outputChecker.checkStringBlock(output, outfile->option.c_str());
        }
    }
}

} // namespace test
} // namespace gmx
