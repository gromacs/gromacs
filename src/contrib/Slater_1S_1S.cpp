#include "slater_low.h"

cl_R Slater_1S_1S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = (5LL*xi)/8LL

    ; } else {  S = (1LL/r)*((-24LL + 24LL*exp(2LL*rxi) - 33LL*rxi - 18LL*Power(rxi,2LL) - 4LL*Power(rxi,3LL))/

      (24LL*exp(2LL*rxi))

    ); }
 
  }
  else {
      if (r == 0LL) { S = (xi*xj*(Power(xi,2LL) + 3LL*xi*xj + Power(xj,2LL)))/Power(xi + xj,3LL)

    ; } else { S = (1LL/r)*((exp(2LL*(rxi + rxj))*Power(Power(rxi,2LL) - Power(rxj,2LL),3LL) + 

        exp(2LL*rxj)*Power(rxj,4LL)*

         (-3LL*Power(rxi,2LL) - Power(rxi,3LL) + Power(rxj,2LL) + rxi*Power(rxj,2LL)) - 

        exp(2LL*rxi)*Power(rxi,4LL)*

         (Power(rxi,2LL)*(1LL + rxj) - Power(rxj,2LL)*(3LL + rxj)))/

      (exp(2LL*(rxi + rxj))*Power(rxi - rxj,3LL)*Power(rxi + rxj,3LL))

    ); }
   
  }
  return S;
}

