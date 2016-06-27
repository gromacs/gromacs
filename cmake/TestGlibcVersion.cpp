#include <features.h>

int main()
{
    // This compiles if the C++ compiler is using glibc with version 2.23+
#if (__GLIBC__ > 2) || (__GLIBC__ == 2 && __GLIBC_MINOR__ >= 23)
    return 0;
#else
#error A glibc version prior to 2.23 was found
#endif
}
