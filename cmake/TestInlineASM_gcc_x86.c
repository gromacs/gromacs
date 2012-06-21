int
main()
{
  unsigned int _eax,_ebx,_ecx,_edx;
  unsigned int level = 0;

  /* Test gcc inline asm for x86 */
#if defined (__x86_64__) || defined (_M_X64)
    __asm__("push %%rbx       \n\t"
            "cpuid            \n\t"
            "movl %%ebx, %1   \n\t"
            "pop %%rbx        \n\t"
            : "=a"(_eax), "=r"(_ebx), "=c"(_ecx), "=d"(_edx) : "0"(level));
#else
    __asm__("push %%ebx       \n\t"
            "cpuid            \n\t"
            "movl %%ebx, %1   \n\t"
            "pop %%ebx        \n\t"
            : "=a"(_eax), "=r"(_ebx), "=c"(_ecx), "=d"(_edx) : "0"(level));
#endif

  return 0;
}    
