#include "debugtrace.h"

#include <ctime>

#include "file.h"
#include "stringutil.h"

namespace gmx
{

class DebugTracer::Impl
{
    public:
        explicit Impl(const std::string &filename)
            : file_(filename, "w"), eventNumber_(0), startTime_(clock())
        {
        }

        File                    file_;
        int                     eventNumber_;
        clock_t                 startTime_;
};

DebugTracer::DebugTracer()
    : impl_(NULL)
{
}

DebugTracer::~DebugTracer()
{
}

void DebugTracer::startTrace(const std::string &filename)
{
    if (filename.empty())
    {
        return;
    }
    impl_.reset(new Impl(filename));
}

void DebugTracer::write(const char *message)
{
    ++impl_->eventNumber_;
    double elapsedTime
        = static_cast<double>(clock() - impl_->startTime_) / CLOCKS_PER_SEC;
    impl_->file_.writeLine(formatString("%d\t#%d\t%s",
                                        static_cast<int>(elapsedTime * 1000),
                                        impl_->eventNumber_,
                                        message));
}

} // namespace gmx
