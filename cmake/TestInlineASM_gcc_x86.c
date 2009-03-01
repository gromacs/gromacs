int
main()
{
  float f;
  int i; 
  /* Test gcc inline asm for x86 */
  asm("fld %1\nfistpl %0\n" : "=m" (*&i) : "f" (f));
  return 0;
}    
