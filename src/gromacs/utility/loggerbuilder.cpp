/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
#include "gmxpre.h"

#include "gromacs/utility/loggerbuilder.h"

#include <cstdio>

#include <memory>
#include <utility>
#include <vector>

#include "gromacs/utility/filestream.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textstream.h"
#include "gromacs/utility/textwriter.h"

namespace gmx
{

class LogTargetCollection : public ILogTarget
{
public:
    void addTarget(ILogTarget* target) { targets_.push_back(target); }

    void writeEntry(const LogEntry& entry) override
    {
        for (ILogTarget* target : targets_)
        {
            target->writeEntry(entry);
        }
    }

private:
    std::vector<ILogTarget*> targets_;
};

class LogTargetFormatter : public ILogTarget
{
public:
    explicit LogTargetFormatter(TextOutputStream* stream) : writer_(stream) {}

    void writeEntry(const LogEntry& entry) override;

private:
    TextWriter writer_;
};


void LogTargetFormatter::writeEntry(const LogEntry& entry)
{
    if (entry.asParagraph)
    {
        writer_.ensureEmptyLine();
    }
    writer_.writeLine(entry.text);
    if (entry.asParagraph)
    {
        writer_.ensureEmptyLine();
    }
}

/********************************************************************
 * LoggerOwner::Impl
 */

class LoggerOwner::Impl
{
public:
    explicit Impl(ILogTarget* loggerTargets[MDLogger::LogLevelCount]) : logger_(loggerTargets) {}

    MDLogger                                       logger_;
    std::vector<std::unique_ptr<TextOutputStream>> streams_;
    std::vector<std::unique_ptr<ILogTarget>>       targets_;
};

/********************************************************************
 * LoggerOwner
 */

LoggerOwner::LoggerOwner(std::unique_ptr<Impl> impl) :
    impl_(impl.release()), logger_(&impl_->logger_)
{
}

LoggerOwner::LoggerOwner(LoggerOwner&& other) noexcept :
    impl_(std::move(other.impl_)), logger_(&impl_->logger_)
{
}

LoggerOwner& LoggerOwner::operator=(LoggerOwner&& other) noexcept
{
    impl_   = std::move(other.impl_);
    logger_ = &impl_->logger_;
    return *this;
}

LoggerOwner::~LoggerOwner() {}

/********************************************************************
 * LoggerBuilder::Impl
 */

class LoggerBuilder::Impl
{
public:
    std::vector<std::unique_ptr<TextOutputStream>> streams_;
    std::vector<std::unique_ptr<ILogTarget>>       targets_;
    std::vector<ILogTarget*>                       loggerTargets_[MDLogger::LogLevelCount];
};

/********************************************************************
 * LoggerBuilder
 */

LoggerBuilder::LoggerBuilder() : impl_(new Impl) {}

LoggerBuilder::~LoggerBuilder() {}

void LoggerBuilder::addTargetStream(MDLogger::LogLevel level, TextOutputStream* stream)
{
    impl_->targets_.push_back(std::unique_ptr<ILogTarget>(new LogTargetFormatter(stream)));
    ILogTarget* target = impl_->targets_.back().get();
    for (int i = 0; i <= static_cast<int>(level); ++i)
    {
        impl_->loggerTargets_[i].push_back(target);
    }
}

void LoggerBuilder::addTargetFile(MDLogger::LogLevel level, FILE* fp)
{
    std::unique_ptr<TextOutputStream> stream(new TextOutputFile(fp));
    addTargetStream(level, stream.get());
    impl_->streams_.push_back(std::move(stream));
}

LoggerOwner LoggerBuilder::build()
{
    ILogTarget* loggerTargets[MDLogger::LogLevelCount];
    for (int i = 0; i < MDLogger::LogLevelCount; ++i)
    {
        auto& levelTargets = impl_->loggerTargets_[i];
        loggerTargets[i]   = nullptr;
        if (!levelTargets.empty())
        {
            if (levelTargets.size() == 1)
            {
                loggerTargets[i] = levelTargets[0];
            }
            else
            {
                std::unique_ptr<LogTargetCollection> collection(new LogTargetCollection);
                for (auto& target : levelTargets)
                {
                    collection->addTarget(target);
                }
                loggerTargets[i] = collection.get();
                impl_->targets_.push_back(std::move(collection));
            }
        }
        levelTargets.clear();
    }
    std::unique_ptr<LoggerOwner::Impl> data(new LoggerOwner::Impl(loggerTargets));
    data->targets_ = std::move(impl_->targets_);
    data->streams_ = std::move(impl_->streams_);
    return LoggerOwner(std::move(data));
}

} // namespace gmx
