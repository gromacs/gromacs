int
main()
{
  /* Check that a double is 8 bytes - compilation dies if it isnt */
  extern char xyz [sizeof(double) == 8 ? 1 : -1];
  int i;
  double d;

  double abc [10] = {
    /* Zero-terminated strings encoded as floating-point numbers */
    /* "GROMACSX" in ascii    */
    (double)  3.80279098314984902657e+35 , (double) 0.0,
    /* "GROMACSX" in ebcdic   */
    (double) -1.37384666579378297437e+38 , (double) 0.0,
    /* "D__float" (vax)       */
    (double)  3.53802595280598432000e+18 , (double) 0.0,
    /* "IBMHEXFP" s390/ascii  */
    (double)  1.77977764695171661377e+10 , (double) 0.0,
    /* "IBMHEXFP" s390/ebcdic */
    (double) -5.22995989424860458374e+10 , (double) 0.0,
  };

  /* Make sure some compilers do not optimize away the entire structure
   * with floating-point data by using it to produce a return value.
   */
  for(i=0,d=0;i<10;i++)
  {
    d+=abc[i];
  }

  return (d==12345.0);
}
