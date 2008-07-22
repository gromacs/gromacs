#include "slater_low.h"

cl_R Slater_1S_2S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = (7LL*xi)/16LL

    ; } else {  S = (1LL/r)*((-240LL + 240LL*exp(2LL*rxi) - 375LL*rxi - 270LL*Power(rxi,2LL) - 115LL*Power(rxi,3LL) - 

        30LL*Power(rxi,4LL) - 4LL*Power(rxi,5LL))/(240LL*exp(2LL*rxi))

    ); }
 
  }
  else {
      if (r == 0LL) { S = (xi*xj*(Power(xi,4LL) + 5LL*Power(xi,3LL)*xj + 10LL*Power(xi,2LL)*Power(xj,2LL) + 

          10LL*xi*Power(xj,3LL) + 2LL*Power(xj,4LL)))/(2LL*Power(xi + xj,5LL))

    ; } else { S = (1LL/r)*((6LL*exp(2LL*(rxi + rxj))*Power(Power(rxi,2LL) - Power(rxj,2LL),5LL) + 

        6LL*exp(2LL*rxj)*Power(rxj,6LL)*

         (-4LL*Power(rxi,4LL) - Power(rxi,5LL) - 5LL*Power(rxi,2LL)*Power(rxj,2LL) + 

           Power(rxj,4LL) + rxi*Power(rxj,4LL)) - 

        exp(2LL*rxi)*Power(rxi,4LL)*

         (Power(rxi,6LL)*(6LL + 9LL*rxj + 6LL*Power(rxj,2LL) + 2LL*Power(rxj,3LL)) - 

           3LL*Power(rxi,4LL)*Power(rxj,2LL)*

            (10LL + 15LL*rxj + 10LL*Power(rxj,2LL) + 2LL*Power(rxj,3LL)) + 

           3LL*Power(rxi,2LL)*Power(rxj,4LL)*

            (20LL + 33LL*rxj + 14LL*Power(rxj,2LL) + 2LL*Power(rxj,3LL)) - 

           Power(rxj,6LL)*(84LL + 63LL*rxj + 18LL*Power(rxj,2LL) + 2LL*Power(rxj,3LL))))/

      (6LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj,5LL)*Power(rxi + rxj,5LL))

    ); }
   
  }
  return S;
}


cl_R Slater_2S_1S(cl_R r,cl_R xi,cl_R xj)
{
  return Slater_1S_2S(r,xj,xi);
}

