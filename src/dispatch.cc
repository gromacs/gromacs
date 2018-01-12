#include <dlfcn.h>
#include <stdlib.h>
#include "gromacs/hardware/cpuinfo.h"

// the __attribute__ ensures that this function is called when the library is loaded
__attribute__((constructor)) void init()
{
    // manually load the appropriate shared library based upon what the CPU supports
    // at runtime
    gmx::CpuInfo  cpuInfo(gmx::CpuInfo::detect());

    void *handle;
    if(cpuInfo.feature(gmx::CpuInfo::Feature::X86_Avx512F)) handle=dlopen("libgromacs_AVX_512.so", RTLD_NOW | RTLD_GLOBAL);
    else if(cpuInfo.feature(gmx::CpuInfo::Feature::X86_Avx2)) handle=dlopen("libgromacs_AVX2_256.so", RTLD_NOW | RTLD_GLOBAL);
    else if(cpuInfo.feature(gmx::CpuInfo::Feature::X86_Avx)) handle=dlopen("libgromacs_AVX_256.so", RTLD_NOW | RTLD_GLOBAL);
    else handle=dlopen("libgromacs_SSE2.so", RTLD_NOW | RTLD_GLOBAL);
    if (!handle) {
        fprintf(stderr, "%s\n", dlerror());
    }
    // NOTE: this is just an example; you should check the return values from
    // dlopen() above and handle errors accordingly
}
