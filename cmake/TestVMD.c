#include "molfile_plugin.h"
#include "vmddlopen.c"
#include "stdio.h"

static int register_cb(void *v, vmdplugin_t *p) {
    *(molfile_plugin_t**)v = (molfile_plugin_t *)p;
    return VMDPLUGIN_SUCCESS;
}

typedef int (*initfunc)(void);
typedef int (*regfunc)(void *, vmdplugin_register_cb);

/*run: gcc TestVMD.c -DGMX_USE_PLUGINS -Wall -ldl src/gmxlib/vmddlopen.c -I src/gmxlib && ./a.out .../xyzplugin.so ; echo $?*/
int main(int argc, char** argv)
{
    void *handle, *ifunc, *registerfunc;
    molfile_plugin_t* api;
    if (argc!=2) return -1;
    handle = vmddlopen(argv[1]);
    if (!handle)
    {
        fprintf(stderr,"%s\n",vmddlerror());
        return 1;
    }
    ifunc = vmddlsym(handle, "vmdplugin_init");
    if (!ifunc || ((initfunc)(ifunc))()) return 2;
    registerfunc = vmddlsym(handle, "vmdplugin_register");
    if (!registerfunc) return 3;
    ((regfunc)registerfunc)(&api, register_cb);
    if (!api) return 4;
    if (api->abiversion<10) return 5;
    return 0;
}
