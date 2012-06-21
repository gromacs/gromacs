int
main()
{
  unsigned int _eax,_ebx,_ecx,_edx;
  unsigned int level = 0;

  /* Test gcc inline asm for x86 */
  __asm__("pushl %%ebx      \n\t"
          "cpuid            \n\t"
          "movl %%ebx, %1   \n\t"
          "popl %%ebx       \n\t"
	  : "=a"(_eax), "=r"(_ebx), "=c"(_ecx), "=d"(_edx) : "0"(level));

  return 0;
}    
