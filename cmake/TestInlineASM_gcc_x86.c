#define __STDC_LIMIT_MACROS
#include <stdint.h>
int
main()
{
  unsigned int _eax,_ebx,_ecx,_edx;
  unsigned int level = 0;

  /* Test gcc inline asm for x86 */
#if defined(UINTPTR_MAX) && defined(UINT32_MAX) && (UINTPTR_MAX==UINT32_MAX)
    __asm__("push %%ebx       \n\t"
            "cpuid            \n\t"
            "movl %%ebx, %1   \n\t"
            "pop %%ebx        \n\t"
            : "=a"(_eax), "=r"(_ebx), "=c"(_ecx), "=d"(_edx) : "0"(level));
#elif defined(UINTPTR_MAX) && defined(UINT32_MAX) && (UINTPTR_MAX>UINT32_MAX)
    __asm__("push %%rbx       \n\t"
            "cpuid            \n\t"
            "movl %%ebx, %1   \n\t"
            "pop %%rbx        \n\t"
            : "=a"(_eax), "=r"(_ebx), "=c"(_ecx), "=d"(_edx) : "0"(level));
#else
#    error Cannot detect whether this is a 32-bit or 64-bit x86 build.
#endif

  return 0;
}    
