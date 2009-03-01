int
main()
{
  float f;
  int i; 
  /* Test microsoft visual studio inline asm for x86 */
  _asm { fld f } ; 
  _asm { fistpl i };
  return 0;
}    
