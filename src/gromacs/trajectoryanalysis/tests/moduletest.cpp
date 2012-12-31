/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief
 * Implements classes in moduletest.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_trajectoryanalysis
 */
#include "moduletest.h"

#include <set>
#include <string>
#include <vector>

#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/cmdlinerunner.h"
#include "gromacs/utility/file.h"

#include "testutils/cmdlinetest.h"
#include "testutils/datatest.h"
#include "testutils/refdata.h"
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
        struct OutputFileInfo
        {
            OutputFileInfo(const char *option, const std::string &path)
                : option(option), path(path)
            {
            }

            std::string         option;
            std::string         path;
        };

        typedef std::set<std::string>       DatasetNames;
        typedef std::vector<OutputFileInfo> OutputFileList;

        explicit Impl(AbstractTrajectoryAnalysisModuleTestFixture *parent);

        TrajectoryAnalysisModule &module();
        void ensureModuleCreated();

        AbstractTrajectoryAnalysisModuleTestFixture    &parent_;
        TrajectoryAnalysisModulePointer                 module_;
        TestReferenceData               data_;
        CommandLine                     cmdline_;
        TestFileManager                 tempFiles_;
        DatasetNames                    moduleDatasets_;
        DatasetNames                    outputDatasets_;
        OutputFileList                  outputFiles_;
        bool                            bDatasetsIncluded_;
};

AbstractTrajectoryAnalysisModuleTestFixture::Impl::Impl(
    AbstractTrajectoryAnalysisModuleTestFixture *parent)
    : parent_(*parent), bDatasetsIncluded_(false)
{
    cmdline_.append("module");
    cmdline_.append("-quiet");
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
        const std::vector<std::string> &datasetNames(module_->datasetNames());
        moduleDatasets_.clear();
        moduleDatasets_.insert(datasetNames.begin(), datasetNames.end());
        outputDatasets_ = moduleDatasets_;
    }
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
AbstractTrajectoryAnalysisModuleTestFixture::includeDataset(const char *name)
{
    impl_->ensureModuleCreated();
    if (!impl_->bDatasetsIncluded_)
    {
        impl_->outputDatasets_.clear();
    }
    bool bFound = (impl_->moduleDatasets_.find(name) != impl_->moduleDatasets_.end());
    GMX_RELEASE_ASSERT(bFound, "Attempted to include a non-existent dataset");
    impl_->outputDatasets_.insert(name);
}

void
AbstractTrajectoryAnalysisModuleTestFixture::excludeDataset(const char *name)
{
    impl_->ensureModuleCreated();
    bool bFound = (impl_->outputDatasets_.erase(name) > 0);
    GMX_RELEASE_ASSERT(bFound, "Attempted to exclude a non-existent dataset");
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

    if (!impl_->outputDatasets_.empty())
    {
        TestReferenceChecker dataChecker(
            rootChecker.checkCompound("OutputData", "Data"));
        Impl::DatasetNames::const_iterator dataset;
        for (dataset = impl_->outputDatasets_.begin();
             dataset != impl_->outputDatasets_.end();
             ++dataset)
        {
            const char *name = dataset->c_str();
            AbstractAnalysisData &dataset = module.datasetFromName(name);
            AnalysisDataTestFixture::addReferenceCheckerModule(
                dataChecker, name, &dataset);
        }
    }

    TrajectoryAnalysisCommandLineRunner runner(&module);
    runner.setPrintCopyright(false);
    int rc = 0;
    EXPECT_NO_THROW(rc = runner.run(impl_->cmdline_.argc(), impl_->cmdline_.argv()));
    EXPECT_EQ(0, rc);

    if (!impl_->outputFiles_.empty())
    {
        TestReferenceChecker outputChecker(
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
