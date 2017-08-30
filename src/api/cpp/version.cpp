#include "apiversion.h"
#include "gmxapi/version.h"

namespace gmxapi
{

unsigned int Version::major()
{
    return GMXAPI_MAJOR;
}

unsigned int Version::minor()
{
    return GMXAPI_MINOR;
}

unsigned int Version::patch()
{
    return GMXAPI_PATCH;
}

std::string Version::release()
{
    return GMXAPI_RELEASE;
}

bool Version::has_feature(std::string featurename)
{
    // For features introduced after version 1.0, we can consult a map somewhere.
    (void)featurename;
    return false;
}

bool Version::is_at_least(unsigned int major, unsigned int minor, unsigned int patch)
{
    if (Version::major() >= major && Version::minor() >= minor && Version::patch() >= patch)
    {
        return true;
    }
    else
    {
        return false;
    };
}

} // end namespace gmxapi
