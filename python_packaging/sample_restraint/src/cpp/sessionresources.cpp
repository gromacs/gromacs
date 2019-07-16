/*! \file
 * \brief Definitions for some useful types and templates for GROMACS restraints.
 *
 * \todo This should be part of a template library installed with GROMACS.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 */

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

void ResourcesHandle::reduce(const Matrix<double>& send,
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

void ResourcesHandle::stop()
{
    assert(session_);
    auto signaller = gmxapi::getMdrunnerSignal(session_,
                                               gmxapi::md::signals::STOP);

    // Should probably check that the function object has been initialized...
    signaller();
}

ResourcesHandle Resources::getHandle() const
{
    auto handle = ResourcesHandle();

    if (!bool(reduce_))
    {
        throw gmxapi::ProtocolError("reduce operation functor is not set, which should not happen...");
    }
    handle.reduce_ = &reduce_;

    if (!session_)
    {
        throw gmxapi::ProtocolError("Resources::getHandle() must not be called before setSession() has been called.");
    }
    handle.session_ = session_;

    return handle;
}

void Resources::setSession(gmxapi::SessionResources* session)
{
    if (!session)
    {
        throw gmxapi::ProtocolError("Resources::setSession received a null SessionResources pointer.");
    }
    session_ = session;
}

} // end namespace myplugin

