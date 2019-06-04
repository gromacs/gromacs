//
// Created by Eric Irrgang on 6/11/18.
//

#include "sessionresources.h"

#include <cassert>

#include <memory>

#include "gmxapi/exceptions.h"
#include "gmxapi/md/mdsignals.h"

namespace plugin
{

// Explicit instantiation.
template
class ::plugin::Matrix<double>;

void EnsembleResourceHandle::reduce(const Matrix<double>& send,
                                    Matrix<double>* receive) const
{
    assert(reduce_);
    if (*reduce_)
    {
        (*reduce_)(send,
               receive);
    }
    else
    {
        throw gmxapi::ProtocolError("'reduce' functor was not initialized before use.");
    }
}

void EnsembleResourceHandle::stop()
{
    assert(session_);
    auto signaller = gmxapi::getMdrunnerSignal(session_,
                                               gmxapi::md::signals::STOP);

    // Should probably check that the function object has been initialized...
    signaller();
}
//
//gmxapi::session::OutputStream* EnsembleResourceHandle::ostream()
//{
//    return ostream_.get();
//}

EnsembleResourceHandle EnsembleResources::getHandle() const
{
    auto handle = EnsembleResourceHandle();

    if (!bool(reduce_))
    {
        throw gmxapi::ProtocolError("reduce operation functor is not set, which should not happen...");
    }
    handle.reduce_ = &reduce_;

    if (!session_)
    {
        throw gmxapi::ProtocolError("EnsembleResources::getHandle() must not be called before setSession() has been called.");
    }
    handle.session_ = session_;

    return handle;
}

void EnsembleResources::setSession(gmxapi::SessionResources* session)
{
    if (!session)
    {
        throw gmxapi::ProtocolError("EnsembleResources::setSession received a null SessionResources pointer.");
    }
    session_ = session;
}
//
//void EnsembleResources::setOutputStream(std::unique_ptr<gmxapi::session::OutputStream> ostream)
//{
//    ostream_ = std::move(ostream);
//}

} // end namespace myplugin

