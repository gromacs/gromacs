/* Lambert W function. 
   Was ~/C/LambertW.c written K M Briggs Keith dot Briggs at bt dot com 97 May 21.  
   Revised KMB 97 Nov 20; 98 Feb 11, Nov 24, Dec 28; 99 Jan 13; 00 Feb 23; 01 Apr 09

   Computes Lambert W function, principal branch.
   See LambertW1.c for -1 branch.

   Returned value W(z) satisfies W(z)*exp(W(z))=z
   test data...
      W(1)= 0.5671432904097838730
      W(2)= 0.8526055020137254914
      W(20)=2.2050032780240599705
   To solve (a+b*R)*exp(-c*R)-d=0 for R, use
   R=-(b*W(-exp(-a*c/b)/b*d*c)+a*c)/b/c
*/

/* Added and adapted to gromacs by Erik Marklund, 20010 Nov 26.
 * Some cleaning up was also done in the process. */

#ifndef _lambertw_h

#ifdef __cplusplus
extern "C" {
#endif

extern real LambertW(const real z);
  
#ifdef __cplusplus
           }
#endif

#endif /* _lambertw_h */
