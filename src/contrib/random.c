/*
 * Print a lookup table for Gaussian numbers with 4 entries on each
 * line, formatted for inclusion in this file. Size is 2^bits.
 */

void
print_gaussian_table(int bits)
{
  int n,nh,i,j;
  double invn,fac,x,invgauss,det,dx;
  real  *table;
  
  n = 1 << bits;
  table = (real *)malloc(n*sizeof(real));
  
  /* Fill a table of size n such that random draws from it
    * produce a Gaussian distribution.
    * We integrate the Gaussian distribution G approximating:
    *   integral(x->x+dx) G(y) dy
    * with:
    *   G(x) dx + G'(x) dx^2/2 = G(x) dx - G(x) x dx^2/2
    * Then we need to find dx such that the integral is 1/n.
    * The last step uses dx = 1/x as the approximation is not accurate enough.
    */
  invn = 1.0/n;
  fac = sqrt(2*M_PI);
  x = 0.5*fac*invn;
  nh = n/2;
  for(i=0; i<nh; i++) {
    if (i > 0) {
      if (i < nh-1) {
	invgauss = fac*exp(0.5*x*x);
	/* det is larger than 0 for all x, except for the last */
	det = 1 - 2*invn*x*invgauss;
	dx = (1 - sqrt(det))/x;
      } else {
	dx = 1/x;
      }
      x = x + dx;
    }
    table[nh-1-i] = -x;
    table[nh+i]   =  x;
  }
  printf("static const real *\ngaussian_table[%d] = {\n",n);
  for(i=0;i<n;i+=4) {
    printf("  ");
    for(j=0;j<4;j++) {
      printf("%14.7e",table[i+j]);
      if(i+j<(n-1))
	printf(",");
    }
    printf("\n");
  }
  printf("};\n");
  free(table);
}

