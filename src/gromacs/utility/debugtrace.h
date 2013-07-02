#ifndef GMX_UTILITY_TRACE_H
#define GMX_UTILITY_TRACE_H

#include "common.h"
#include "stringutil.h"

namespace gmx
{

class DebugTracer
{
    public:
        DebugTracer();
        ~DebugTracer();

        bool enabled() const { return impl_.get() != NULL; }
        void startTrace(const std::string &filename);

        void write(const char *message);
        void write(const std::string &message)
        {
            write(message.c_str());
        }

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

#define GMX_TRACE(tracer, message) \
    if (!tracer.enabled()) \
        ; \
    else \
        tracer.write(message)

/*
class DebugTraceGroup
{
    public:
        DebugTraceGroup(const char *name, DebugTrace *trace);
        ~DebugTraceGroup();

        bool enabled() const { return trace_.enabled(); }

        void write(const char *message)
        {
            trace_.write(message);
        }
        void write(const std::string &message)
        {
            trace_.write(message.c_str());
        }
        void write(const char *fmt, ...)
        {
            if (enabled())
            {
                va_list ap;
                va_start(ap, fmt);
                try
                {
                    trace_.write(formatStringV(fmt, ap).c_str());
                }
                catch (...)
                {
                    va_end(ap);
                    throw;
                }
                va_end(ap);
            }
        }

    private:
        class Impl;

        DebugTrace               &trace_;
        PrivateImplPointer<Impl>  impl_;
};
*/

} // namespace gmx

#endif
