// slater_integrals.cpp (c) 2008 Paul J. van Maaren and David van der Spoel
#include <iostream>
using namespace std;
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef HAVE_CLN_CLN_H
#include <cln/cln.h>
#include "slater_integrals.h"
using namespace cln;

#define PRECISION 80

static cl_R     ZERO = "0.0_80";
static cl_R      ONE = "1.0_80";
static cl_R      TWO = "2.0_80";
static cl_R    THREE = "3.0_80";
static cl_R     FOUR = "4.0_80";
static cl_R     FIVE = "5.0_80";
static cl_R      SIX = "6.0_80";
static cl_R    SEVEN = "7.0_80";
static cl_R    EIGHT = "8.0_80";
static cl_R     NINE = "9.0_80";
static float_format_t precision = float_format(80);

cl_R Power(cl_R a,int b)
{
  if (b < 0) { cerr << "negative exponent in Power" << endl; exit(1); }
  if (b == 0) return ONE;
  if (a == ZERO) return ZERO;
  if ((b % 2) == 0) return Power(a*a,b/2);
  else if ((b % 2) == 1) return a*Power(a*a,b/2);
  return ZERO;
}

cl_R Slater_1S_1S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = (5LL*xi)/8LL

    ; } else {  S = (1LL/r)*(-24LL + 24LL*exp(2LL*rxi) - 33LL*rxi - 18LL*Power(rxi,2LL) - 4LL*Power(rxi,3LL))/

      (24LL*exp(2LL*rxi))

    ; }
 
  }
  else {
      if (r == 0LL) { S = (xi*xj*(Power(xi,2LL) + 3LL*xi*xj + Power(xj,2LL)))/Power(xi + xj,3LL)

    ; } else { S = (1LL/r)*(exp(2LL*(rxi + rxj))*Power(Power(rxi,2LL) - Power(rxj,2LL),3LL) + 

        exp(2LL*rxj)*Power(rxj,4LL)*

         (-3LL*Power(rxi,2LL) - Power(rxi,3LL) + Power(rxj,2LL) + rxi*Power(rxj,2LL)) - 

        exp(2LL*rxi)*Power(rxi,4LL)*

         (Power(rxi,2LL)*(1LL + rxj) - Power(rxj,2LL)*(3LL + rxj)))/

      (exp(2LL*(rxi + rxj))*Power(rxi - rxj,3LL)*Power(rxi + rxj,3LL))

     ; }
   
  }
  return S;
}

cl_R Slater_1S_2S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = (7LL*xi)/16LL

    ; } else {  S = (1LL/r)*(-240LL + 240LL*exp(2LL*rxi) - 375LL*rxi - 270LL*Power(rxi,2LL) - 115LL*Power(rxi,3LL) - 

        30LL*Power(rxi,4LL) - 4LL*Power(rxi,5LL))/(240LL*exp(2LL*rxi))

    ; }
 
  }
  else {
      if (r == 0LL) { S = (xi*xj*(Power(xi,4LL) + 5LL*Power(xi,3LL)*xj + 10LL*Power(xi,2LL)*Power(xj,2LL) + 

          10LL*xi*Power(xj,3LL) + 2LL*Power(xj,4LL)))/(2LL*Power(xi + xj,5LL))

    ; } else { S = (1LL/r)*(6LL*exp(2LL*(rxi + rxj))*Power(Power(rxi,2LL) - Power(rxj,2LL),5LL) + 

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

     ; }
   
  }
  return S;
}

cl_R Slater_1S_3S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = (41LL*xi)/128LL

    ; } else {  S = (1LL/r)*(-120960LL + 120960LL*exp(2LL*rxi) - 203175LL*rxi - 164430LL*Power(rxi,2LL) - 

        84420LL*Power(rxi,3LL) - 30240LL*Power(rxi,4LL) - 7728LL*Power(rxi,5LL) - 

        1344LL*Power(rxi,6LL) - 128LL*Power(rxi,7LL))/(120960LL*exp(2LL*rxi))

    ; }
 
  }
  else {
      if (r == 0LL) { S = (xi*xj*(Power(xi,6LL) + 7LL*Power(xi,5LL)*xj + 21LL*Power(xi,4LL)*Power(xj,2LL) + 

          35LL*Power(xi,3LL)*Power(xj,3LL) + 35LL*Power(xi,2LL)*Power(xj,4LL) + 

          21LL*xi*Power(xj,5LL) + 3LL*Power(xj,6LL)))/(3LL*Power(xi + xj,7LL))

    ; } else { S = (1LL/r)*(45LL*exp(2LL*(rxi + rxj))*Power(Power(rxi,2LL) - Power(rxj,2LL),7LL) + 

        15LL*exp(2LL*rxj)*Power(rxj,8LL)*

         (-15LL*Power(rxi,6LL) - 3LL*Power(rxi,7LL) - 63LL*Power(rxi,4LL)*Power(rxj,2LL) - 

           7LL*Power(rxi,5LL)*Power(rxj,2LL) - 21LL*Power(rxi,2LL)*Power(rxj,4LL) + 

           7LL*Power(rxi,3LL)*Power(rxj,4LL) + 3LL*Power(rxj,6LL) + 3LL*rxi*Power(rxj,6LL)) + 

        exp(2LL*rxi)*Power(rxi,4LL)*

         (-10LL*Power(rxi,2LL)*Power(rxj,8LL)*

            (135LL + 333LL*rxj + 228LL*Power(rxj,2LL) + 75LL*Power(rxj,3LL) + 

              13LL*Power(rxj,4LL) + Power(rxj,5LL)) + 

           2LL*Power(rxj,10LL)*(945LL + 945LL*rxj + 420LL*Power(rxj,2LL) + 

              105LL*Power(rxj,3LL) + 15LL*Power(rxj,4LL) + Power(rxj,5LL)) - 

           Power(rxi,10LL)*(45LL + 75LL*rxj + 60LL*Power(rxj,2LL) + 30LL*Power(rxj,3LL) + 

              10LL*Power(rxj,4LL) + 2LL*Power(rxj,5LL)) + 

           5LL*Power(rxi,8LL)*Power(rxj,2LL)*

            (63LL + 105LL*rxj + 84LL*Power(rxj,2LL) + 42LL*Power(rxj,3LL) + 

              14LL*Power(rxj,4LL) + 2LL*Power(rxj,5LL)) - 

           5LL*Power(rxi,6LL)*Power(rxj,4LL)*

            (189LL + 315LL*rxj + 252LL*Power(rxj,2LL) + 132LL*Power(rxj,3LL) + 

              36LL*Power(rxj,4LL) + 4LL*Power(rxj,5LL)) + 

           5LL*Power(rxi,4LL)*Power(rxj,6LL)*

            (315LL + 513LL*rxj + 468LL*Power(rxj,2LL) + 204LL*Power(rxj,3LL) + 

              44LL*Power(rxj,4LL) + 4LL*Power(rxj,5LL))))/

      (45LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj,7LL)*Power(rxi + rxj,7LL))

     ; }
   
  }
  return S;
}

cl_R Slater_1S_4S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = (253LL*xi)/1024LL

    ; } else {  S = (1LL/r)*(-2903040LL + 2903040LL*exp(2LL*rxi) - 5088825LL*rxi - 4371570LL*Power(rxi,2LL) - 

        2439990LL*Power(rxi,3LL) - 986580LL*Power(rxi,4LL) - 303912LL*Power(rxi,5LL) - 

        72576LL*Power(rxi,6LL) - 13248LL*Power(rxi,7LL) - 1728LL*Power(rxi,8LL) - 

        128LL*Power(rxi,9LL))/(2.90304e6*exp(2LL*rxi))

    ; }
 
  }
  else {
      if (r == 0LL) { S = (xi*xj*(Power(xi,8LL) + 9LL*Power(xi,7LL)*xj + 36LL*Power(xi,6LL)*Power(xj,2LL) + 

          84LL*Power(xi,5LL)*Power(xj,3LL) + 126LL*Power(xi,4LL)*Power(xj,4LL) + 

          126LL*Power(xi,3LL)*Power(xj,5LL) + 84LL*Power(xi,2LL)*Power(xj,6LL) + 

          36LL*xi*Power(xj,7LL) + 4LL*Power(xj,8LL)))/(4LL*Power(xi + xj,9LL))

    ; } else { S = (1LL/r)*(1260LL*exp(2LL*(rxi + rxj))*Power(Power(rxi,2LL) - Power(rxj,2LL),9LL) + 

        1260LL*exp(2LL*rxj)*Power(rxj,10LL)*

         (-6LL*Power(rxi,8LL) - Power(rxi,9LL) - 51LL*Power(rxi,6LL)*Power(rxj,2LL) - 

           6LL*Power(rxi,7LL)*Power(rxj,2LL) - 63LL*Power(rxi,4LL)*Power(rxj,4LL) - 

           9LL*Power(rxi,2LL)*Power(rxj,6LL) + 6LL*Power(rxi,3LL)*Power(rxj,6LL) + 

           Power(rxj,8LL) + rxi*Power(rxj,8LL)) - 

        exp(2LL*rxi)*Power(rxi,4LL)*

         (42LL*Power(rxi,10LL)*Power(rxj,4LL)*

            (1080LL + 1890LL*rxj + 1620LL*Power(rxj,2LL) + 900LL*Power(rxj,3LL) + 

              360LL*Power(rxj,4LL) + 111LL*Power(rxj,5LL) + 22LL*Power(rxj,6LL) + 

              2LL*Power(rxj,7LL)) - 70LL*Power(rxi,8LL)*Power(rxj,6LL)*

            (1512LL + 2646LL*rxj + 2268LL*Power(rxj,2LL) + 1248LL*Power(rxj,3LL) + 

              528LL*Power(rxj,4LL) + 153LL*Power(rxj,5LL) + 26LL*Power(rxj,6LL) + 

              2LL*Power(rxj,7LL)) + 14LL*Power(rxi,2LL)*Power(rxj,12LL)*

            (2970LL + 16335LL*rxj + 15390LL*Power(rxj,2LL) + 7110LL*Power(rxj,3LL) + 

              1980LL*Power(rxj,4LL) + 351LL*Power(rxj,5LL) + 38LL*Power(rxj,6LL) + 

              2LL*Power(rxj,7LL)) - 2LL*Power(rxj,14LL)*

            (62370LL + 72765LL*rxj + 39690LL*Power(rxj,2LL) + 13230LL*Power(rxj,3LL) + 

              2940LL*Power(rxj,4LL) + 441LL*Power(rxj,5LL) + 42LL*Power(rxj,6LL) + 

              2LL*Power(rxj,7LL)) + Power(rxi,14LL)*

            (1260LL + 2205LL*rxj + 1890LL*Power(rxj,2LL) + 1050LL*Power(rxj,3LL) + 

              420LL*Power(rxj,4LL) + 126LL*Power(rxj,5LL) + 28LL*Power(rxj,6LL) + 

              4LL*Power(rxj,7LL)) - 7LL*Power(rxi,12LL)*Power(rxj,2LL)*

            (1620LL + 2835LL*rxj + 2430LL*Power(rxj,2LL) + 1350LL*Power(rxj,3LL) + 

              540LL*Power(rxj,4LL) + 162LL*Power(rxj,5LL) + 36LL*Power(rxj,6LL) + 

              4LL*Power(rxj,7LL)) + 35LL*Power(rxi,6LL)*Power(rxj,8LL)*

            (4536LL + 7983LL*rxj + 6534LL*Power(rxj,2LL) + 4014LL*Power(rxj,3LL) + 

              1644LL*Power(rxj,4LL) + 414LL*Power(rxj,5LL) + 60LL*Power(rxj,6LL) + 

              4LL*Power(rxj,7LL)) - 21LL*Power(rxi,4LL)*Power(rxj,10LL)*

            (7920LL + 11385LL*rxj + 12330LL*Power(rxj,2LL) + 7410LL*Power(rxj,3LL) + 

              2580LL*Power(rxj,4LL) + 546LL*Power(rxj,5LL) + 68LL*Power(rxj,6LL) + 

              4LL*Power(rxj,7LL))))/

      (1260LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj,9LL)*Power(rxi + rxj,9LL))

     ; }
   
  }
  return S;
}

cl_R Slater_1S_5S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = (2041LL*xi)/10240LL

    ; } else {  S = (1LL/r)*(-1596672000LL + 1596672000LL*exp(2LL*rxi) - 2875101075LL*rxi - 

        2556858150LL*Power(rxi,2LL) - 1492929900LL*Power(rxi,3LL) - 

        641163600LL*Power(rxi,4LL) - 214719120LL*Power(rxi,5LL) - 

        57879360LL*Power(rxi,6LL) - 12735360LL*Power(rxi,7LL) - 2280960LL*Power(rxi,8LL) - 

        323840LL*Power(rxi,9LL) - 33792LL*Power(rxi,10LL) - 2048LL*Power(rxi,11LL))/

      (1.596672e9*exp(2LL*rxi))

    ; }
 
  }
  else {
      if (r == 0LL) { S = (xi*xj*(Power(xi,10LL) + 11LL*Power(xi,9LL)*xj + 55LL*Power(xi,8LL)*Power(xj,2LL) + 

          165LL*Power(xi,7LL)*Power(xj,3LL) + 330LL*Power(xi,6LL)*Power(xj,4LL) + 

          462LL*Power(xi,5LL)*Power(xj,5LL) + 462LL*Power(xi,4LL)*Power(xj,6LL) + 

          330LL*Power(xi,3LL)*Power(xj,7LL) + 165LL*Power(xi,2LL)*Power(xj,8LL) + 

          55LL*xi*Power(xj,9LL) + 5LL*Power(xj,10LL)))/(5LL*Power(xi + xj,11LL))

    ; } else { S = (1LL/r)*(14175LL*exp(2LL*(rxi + rxj))*Power(Power(rxi,2LL) - Power(rxj,2LL),11LL) + 

        2835LL*exp(2LL*rxj)*Power(rxj,12LL)*

         (-35LL*Power(rxi,10LL) - 5LL*Power(rxi,11LL) - 495LL*Power(rxi,8LL)*Power(rxj,2LL) - 

           55LL*Power(rxi,9LL)*Power(rxj,2LL) - 1254LL*Power(rxi,6LL)*Power(rxj,4LL) - 

           66LL*Power(rxi,7LL)*Power(rxj,4LL) - 726LL*Power(rxi,4LL)*Power(rxj,6LL) + 

           66LL*Power(rxi,5LL)*Power(rxj,6LL) - 55LL*Power(rxi,2LL)*Power(rxj,8LL) + 

           55LL*Power(rxi,3LL)*Power(rxj,8LL) + 5LL*Power(rxj,10LL) + 5LL*rxi*Power(rxj,10LL)) 

    - exp(2LL*rxi)*Power(rxi,4LL)*(Power(rxi,18LL)*

            (14175LL + 25515LL*rxj + 22680LL*Power(rxj,2LL) + 13230LL*Power(rxj,3LL) + 

              5670LL*Power(rxj,4LL) + 1890LL*Power(rxj,5LL) + 504LL*Power(rxj,6LL) + 

              108LL*Power(rxj,7LL) + 18LL*Power(rxj,8LL) + 2LL*Power(rxj,9LL)) - 

           9LL*Power(rxi,16LL)*Power(rxj,2LL)*

            (17325LL + 31185LL*rxj + 27720LL*Power(rxj,2LL) + 16170LL*Power(rxj,3LL) + 

              6930LL*Power(rxj,4LL) + 2310LL*Power(rxj,5LL) + 616LL*Power(rxj,6LL) + 

              132LL*Power(rxj,7LL) + 22LL*Power(rxj,8LL) + 2LL*Power(rxj,9LL)) + 

           126LL*Power(rxi,10LL)*Power(rxj,8LL)*

            (37125LL + 66825LL*rxj + 59400LL*Power(rxj,2LL) + 34725LL*Power(rxj,3LL) + 

              14625LL*Power(rxj,4LL) + 5043LL*Power(rxj,5LL) + 1396LL*Power(rxj,6LL) + 

              276LL*Power(rxj,7LL) + 34LL*Power(rxj,8LL) + 2LL*Power(rxj,9LL)) - 

           126LL*Power(rxi,8LL)*Power(rxj,10LL)*

            (51975LL + 93420LL*rxj + 84240LL*Power(rxj,2LL) + 46815LL*Power(rxj,3LL) + 

              20835LL*Power(rxj,4LL) + 7485LL*Power(rxj,5LL) + 1964LL*Power(rxj,6LL) + 

              348LL*Power(rxj,7LL) + 38LL*Power(rxj,8LL) + 2LL*Power(rxj,9LL)) + 

           9LL*Power(rxi,2LL)*Power(rxj,16LL)*

            (-135135LL + 405405LL*rxj + 582120LL*Power(rxj,2LL) + 346500LL*Power(rxj,3LL) + 

              124740LL*Power(rxj,4LL) + 30492LL*Power(rxj,5LL) + 5264LL*Power(rxj,6LL) + 

              636LL*Power(rxj,7LL) + 50LL*Power(rxj,8LL) + 2LL*Power(rxj,9LL)) - 

           Power(rxj,18LL)*(2837835LL + 3648645LL*rxj + 2245320LL*Power(rxj,2LL) + 

              873180LL*Power(rxj,3LL) + 238140LL*Power(rxj,4LL) + 47628LL*Power(rxj,5LL) + 

              7056LL*Power(rxj,6LL) + 756LL*Power(rxj,7LL) + 54LL*Power(rxj,8LL) + 

              2LL*Power(rxj,9LL)) + 9LL*Power(rxi,14LL)*Power(rxj,4LL)*

            (86625LL + 155925LL*rxj + 138600LL*Power(rxj,2LL) + 80850LL*Power(rxj,3LL) + 

              34650LL*Power(rxj,4LL) + 11550LL*Power(rxj,5LL) + 3080LL*Power(rxj,6LL) + 

              672LL*Power(rxj,7LL) + 104LL*Power(rxj,8LL) + 8LL*Power(rxj,9LL)) - 

           21LL*Power(rxi,12LL)*Power(rxj,6LL)*

            (111375LL + 200475LL*rxj + 178200LL*Power(rxj,2LL) + 103950LL*Power(rxj,3LL) + 

              44550LL*Power(rxj,4LL) + 14778LL*Power(rxj,5LL) + 4056LL*Power(rxj,6LL) + 

              864LL*Power(rxj,7LL) + 120LL*Power(rxj,8LL) + 8LL*Power(rxj,9LL)) + 

           21LL*Power(rxi,6LL)*Power(rxj,12LL)*

            (307125LL + 594945LL*rxj + 456840LL*Power(rxj,2LL) + 281790LL*Power(rxj,3LL) + 

              137430LL*Power(rxj,4LL) + 47250LL*Power(rxj,5LL) + 11064LL*Power(rxj,6LL) + 

              1728LL*Power(rxj,7LL) + 168LL*Power(rxj,8LL) + 8LL*Power(rxj,9LL)) - 

           9LL*Power(rxi,4LL)*Power(rxj,14LL)*

            (675675LL + 675675LL*rxj + 748440LL*Power(rxj,2LL) + 561330LL*Power(rxj,3LL) + 

              256410LL*Power(rxj,4LL) + 76230LL*Power(rxj,5LL) + 15400LL*Power(rxj,6LL) + 

              2112LL*Power(rxj,7LL) + 184LL*Power(rxj,8LL) + 8LL*Power(rxj,9LL))))/

      (14175LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj,11LL)*Power(rxi + rxj,11LL))

     ; }
   
  }
  return S;
}

cl_R Slater_2S_2S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = (93LL*xi)/256LL

    ; } else {  S = (1LL/r)*(-80640LL + 80640LL*exp(2LL*rxi) - 131985LL*rxi - 102690LL*Power(rxi,2LL) - 

        49980LL*Power(rxi,3LL) - 16800LL*Power(rxi,4LL) - 4032LL*Power(rxi,5LL) - 

        672LL*Power(rxi,6LL) - 64LL*Power(rxi,7LL))/(80640LL*exp(2LL*rxi))

    ; }
 
  }
  else {
      if (r == 0LL) { S = (xi*xj*(Power(xi,6LL) + 7LL*Power(xi,5LL)*xj + 21LL*Power(xi,4LL)*Power(xj,2LL) + 

          35LL*Power(xi,3LL)*Power(xj,3LL) + 21LL*Power(xi,2LL)*Power(xj,4LL) + 

          7LL*xi*Power(xj,5LL) + Power(xj,6LL)))/(2LL*Power(xi + xj,7LL))

    ; } else { S = (1LL/r)*(6LL*exp(2LL*(rxi + rxj))*Power(Power(rxi,2LL) - Power(rxj,2LL),7LL) - 

        exp(2LL*rxi)*Power(rxi,6LL)*

         (21LL*Power(rxi,4LL)*Power(rxj,4LL)*(6LL + 11LL*rxj + 2LL*Power(rxj,2LL)) - 

           2LL*Power(rxj,8LL)*(90LL + 54LL*rxj + 12LL*Power(rxj,2LL) + Power(rxj,3LL)) + 

           Power(rxi,8LL)*(6LL + 9LL*rxj + 6LL*Power(rxj,2LL) + 2LL*Power(rxj,3LL)) + 

           Power(rxi,2LL)*Power(rxj,6LL)*

            (-390LL - 69LL*rxj + 18LL*Power(rxj,2LL) + 4LL*Power(rxj,3LL)) - 

           Power(rxi,6LL)*Power(rxj,2LL)*

            (42LL + 63LL*rxj + 42LL*Power(rxj,2LL) + 4LL*Power(rxj,3LL))) + 

        exp(2LL*rxj)*Power(rxj,6LL)*

         (-24LL*Power(rxi,10LL) - 2LL*Power(rxi,11LL) - 69LL*Power(rxi,7LL)*Power(rxj,2LL) + 

           6LL*Power(rxj,8LL) + 9LL*rxi*Power(rxj,8LL) + 

           4LL*Power(rxi,9LL)*(-27LL + Power(rxj,2LL)) + 

           18LL*Power(rxi,8LL)*(-10LL + Power(rxj,2LL)) + 

           6LL*Power(rxi,2LL)*Power(rxj,6LL)*(-7LL + Power(rxj,2LL)) - 

           42LL*Power(rxi,4LL)*Power(rxj,4LL)*(-3LL + Power(rxj,2LL)) + 

           Power(rxi,3LL)*Power(rxj,6LL)*(-63LL + 2LL*Power(rxj,2LL)) + 

           6LL*Power(rxi,6LL)*Power(rxj,2LL)*(-65LL + 7LL*Power(rxj,2LL)) + 

           Power(rxi,5LL)*(231LL*Power(rxj,4LL) - 4LL*Power(rxj,6LL))))/

      (6LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj,7LL)*Power(rxi + rxj,7LL))

     ; }
   
  }
  return S;
}

cl_R Slater_2S_3S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = (451LL*xi)/1536LL

    ; } else {  S = (1LL/r)*(-4354560LL + 4354560LL*exp(2LL*rxi) - 7430535LL*rxi - 6151950LL*Power(rxi,2LL) - 

        3275370LL*Power(rxi,3LL) - 1251180LL*Power(rxi,4LL) - 361368LL*Power(rxi,5LL) - 

        80640LL*Power(rxi,6LL) - 13824LL*Power(rxi,7LL) - 1728LL*Power(rxi,8LL) - 

        128LL*Power(rxi,9LL))/(4.35456e6*exp(2LL*rxi))

    ; }
 
  }
  else {
      if (r == 0LL) { S = (xi*xj*(2LL*Power(xi,8LL) + 18LL*Power(xi,7LL)*xj + 72LL*Power(xi,6LL)*Power(xj,2LL) + 

          168LL*Power(xi,5LL)*Power(xj,3LL) + 252LL*Power(xi,4LL)*Power(xj,4LL) + 

          252LL*Power(xi,3LL)*Power(xj,5LL) + 108LL*Power(xi,2LL)*Power(xj,6LL) + 

          27LL*xi*Power(xj,7LL) + 3LL*Power(xj,8LL)))/(6LL*Power(xi + xj,9LL))

    ; } else { S = (1LL/r)*(90LL*exp(2LL*(rxi + rxj))*Power(Power(rxi,2LL) - Power(rxj,2LL),9LL) + 

        5LL*exp(2LL*rxj)*Power(rxj,8LL)*

         (-90LL*Power(rxi,12LL) - 6LL*Power(rxi,13LL) + 18LL*Power(rxj,10LL) + 

           27LL*rxi*Power(rxj,10LL) + 18LL*Power(rxi,2LL)*Power(rxj,8LL)*

            (-9LL + Power(rxj,2LL)) - 162LL*Power(rxi,4LL)*Power(rxj,6LL)*

            (-4LL + Power(rxj,2LL)) - 198LL*Power(rxi,10LL)*(5LL + Power(rxj,2LL)) - 

           108LL*Power(rxi,6LL)*Power(rxj,4LL)*(36LL + Power(rxj,2LL)) + 

           2LL*Power(rxi,5LL)*Power(rxj,6LL)*(675LL + Power(rxj,2LL)) - 

           18LL*Power(rxi,7LL)*Power(rxj,4LL)*(-81LL + 2LL*Power(rxj,2LL)) + 

           3LL*Power(rxi,3LL)*Power(rxj,8LL)*(-81LL + 2LL*Power(rxj,2LL)) - 

           Power(rxi,11LL)*(495LL + 2LL*Power(rxj,2LL)) + 

           9LL*Power(rxi,9LL)*Power(rxj,2LL)*(-233LL + 4LL*Power(rxj,2LL)) + 

           6LL*Power(rxi,8LL)*Power(rxj,2LL)*(-1063LL + 90LL*Power(rxj,2LL))) - 

        2LL*exp(2LL*rxi)*Power(rxi,6LL)*

         (-90LL*Power(rxi,6LL)*Power(rxj,6LL)*

            (42LL + 65LL*rxj + 76LL*Power(rxj,2LL) + 22LL*Power(rxj,3LL) + 2LL*Power(rxj,4LL)) - 

           2LL*Power(rxj,12LL)*(2970LL + 2475LL*rxj + 900LL*Power(rxj,2LL) + 

              180LL*Power(rxj,3LL) + 20LL*Power(rxj,4LL) + Power(rxj,5LL)) + 

           10LL*Power(rxi,8LL)*Power(rxj,4LL)*

            (162LL + 270LL*rxj + 216LL*Power(rxj,2LL) + 122LL*Power(rxj,3LL) + 

              22LL*Power(rxj,4LL) + Power(rxj,5LL)) - 

           5LL*Power(rxi,4LL)*Power(rxj,8LL)*

            (-639LL - 3555LL*rxj - 1452LL*Power(rxj,2LL) - 174LL*Power(rxj,3LL) + 

              6LL*Power(rxj,4LL) + 2LL*Power(rxj,5LL)) + 

           Power(rxi,12LL)*(45LL + 75LL*rxj + 60LL*Power(rxj,2LL) + 30LL*Power(rxj,3LL) + 

              10LL*Power(rxj,4LL) + 2LL*Power(rxj,5LL)) - 

           Power(rxi,10LL)*Power(rxj,2LL)*

            (405LL + 675LL*rxj + 540LL*Power(rxj,2LL) + 270LL*Power(rxj,3LL) + 

              90LL*Power(rxj,4LL) + 8LL*Power(rxj,5LL)) + 

           Power(rxi,2LL)*Power(rxj,10LL)*

            (-21615LL - 9075LL*rxj - 300LL*Power(rxj,2LL) + 490LL*Power(rxj,3LL) + 

              110LL*Power(rxj,4LL) + 8LL*Power(rxj,5LL))))/

      (90LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj,9LL)*Power(rxi + rxj,9LL))

     ; }
   
  }
  return S;
}

cl_R Slater_2S_4S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = (975LL*xi)/4096LL

    ; } else {  S = (1LL/r)*(-638668800LL + 638668800LL*exp(2LL*rxi) - 1125310725LL*rxi - 

        973283850LL*Power(rxi,2LL) - 549063900LL*Power(rxi,3LL) - 

        226195200LL*Power(rxi,4LL) - 72099720LL*Power(rxi,5LL) - 18350640LL*Power(rxi,6LL) - 

        3785760LL*Power(rxi,7LL) - 633600LL*Power(rxi,8LL) - 84480LL*Power(rxi,9LL) - 

        8448LL*Power(rxi,10LL) - 512LL*Power(rxi,11LL))/(6.386688e8*exp(2LL*rxi))

    ; }
 
  }
  else {
      if (r == 0LL) { S = (xi*xj*(Power(xi,10LL) + 11LL*Power(xi,9LL)*xj + 55LL*Power(xi,8LL)*Power(xj,2LL) + 

          165LL*Power(xi,7LL)*Power(xj,3LL) + 330LL*Power(xi,6LL)*Power(xj,4LL) + 

          462LL*Power(xi,5LL)*Power(xj,5LL) + 462LL*Power(xi,4LL)*Power(xj,6LL) + 

          330LL*Power(xi,3LL)*Power(xj,7LL) + 110LL*Power(xi,2LL)*Power(xj,8LL) + 

          22LL*xi*Power(xj,9LL) + 2LL*Power(xj,10LL)))/(4LL*Power(xi + xj,11LL))

    ; } else { S = (1LL/r)*(1260LL*exp(2LL*(rxi + rxj))*Power(Power(rxi,2LL) - Power(rxj,2LL),11LL) + 

        210LL*exp(2LL*rxj)*Power(rxj,10LL)*

         (-36LL*Power(rxi,14LL) - 2LL*Power(rxi,15LL) - 

           1287LL*Power(rxi,9LL)*Power(rxj,4LL) + 6LL*Power(rxj,12LL) + 

           9LL*rxi*Power(rxj,12LL) - 22LL*Power(rxi,7LL)*Power(rxj,6LL)*

            (-135LL + Power(rxj,2LL)) + 

           6LL*Power(rxi,2LL)*Power(rxj,10LL)*(-11LL + Power(rxj,2LL)) - 

           66LL*Power(rxi,4LL)*Power(rxj,8LL)*(-5LL + Power(rxj,2LL)) + 

           8LL*Power(rxi,5LL)*Power(rxj,8LL)*(99LL + Power(rxj,2LL)) + 

           Power(rxi,3LL)*Power(rxj,10LL)*(-99LL + 2LL*Power(rxj,2LL)) - 

           132LL*Power(rxi,6LL)*Power(rxj,6LL)*(27LL + 2LL*Power(rxj,2LL)) - 

           78LL*Power(rxi,12LL)*(7LL + 3LL*Power(rxj,2LL)) - 

           2LL*Power(rxi,13LL)*(117LL + 4LL*Power(rxj,2LL)) + 

           66LL*Power(rxi,8LL)*Power(rxj,4LL)*(-191LL + 6LL*Power(rxj,2LL)) + 

           Power(rxi,11LL)*Power(rxj,2LL)*(-2151LL + 22LL*Power(rxj,2LL)) + 

           6LL*Power(rxi,10LL)*Power(rxj,2LL)*(-1099LL + 33LL*Power(rxj,2LL))) - 

        exp(2LL*rxi)*Power(rxi,6LL)*

         (385LL*Power(rxi,8LL)*Power(rxj,8LL)*

            (1080LL + 1935LL*rxj + 1350LL*Power(rxj,2LL) + 1170LL*Power(rxj,3LL) + 

              420LL*Power(rxj,4LL) + 66LL*Power(rxj,5LL) + 4LL*Power(rxj,6LL)) - 

           4LL*Power(rxj,16LL)*(135135LL + 135135LL*rxj + 62370LL*Power(rxj,2LL) + 

              17325LL*Power(rxj,3LL) + 3150LL*Power(rxj,4LL) + 378LL*Power(rxj,5LL) + 

              28LL*Power(rxj,6LL) + Power(rxj,7LL)) + 

           Power(rxi,16LL)*(1260LL + 2205LL*rxj + 1890LL*Power(rxj,2LL) + 

              1050LL*Power(rxj,3LL) + 420LL*Power(rxj,4LL) + 126LL*Power(rxj,5LL) + 

              28LL*Power(rxj,6LL) + 4LL*Power(rxj,7LL)) + 

           7LL*Power(rxi,6LL)*Power(rxj,10LL)*

            (-99540LL - 58095LL*rxj - 190710LL*Power(rxj,2LL) - 100950LL*Power(rxj,3LL) - 

              21660LL*Power(rxj,4LL) - 1938LL*Power(rxj,5LL) - 4LL*Power(rxj,6LL) + 

              8LL*Power(rxj,7LL)) - 7LL*Power(rxi,4LL)*Power(rxj,12LL)*

            (114660LL - 343395LL*rxj - 242910LL*Power(rxj,2LL) - 61950LL*Power(rxj,3LL) - 

              6060LL*Power(rxj,4LL) + 282LL*Power(rxj,5LL) + 116LL*Power(rxj,6LL) + 

              8LL*Power(rxj,7LL)) + 7LL*Power(rxi,12LL)*Power(rxj,4LL)*

            (9900LL + 17325LL*rxj + 14850LL*Power(rxj,2LL) + 8250LL*Power(rxj,3LL) + 

              3300LL*Power(rxj,4LL) + 1074LL*Power(rxj,5LL) + 164LL*Power(rxj,6LL) + 

              8LL*Power(rxj,7LL)) - 7LL*Power(rxi,10LL)*Power(rxj,6LL)*

            (29700LL + 51975LL*rxj + 44550LL*Power(rxj,2LL) + 23850LL*Power(rxj,3LL) + 

              11700LL*Power(rxj,4LL) + 2814LL*Power(rxj,5LL) + 284LL*Power(rxj,6LL) + 

              8LL*Power(rxj,7LL)) - Power(rxi,14LL)*Power(rxj,2LL)*

            (13860LL + 24255LL*rxj + 20790LL*Power(rxj,2LL) + 11550LL*Power(rxj,3LL) + 

              4620LL*Power(rxj,4LL) + 1386LL*Power(rxj,5LL) + 308LL*Power(rxj,6LL) + 

              24LL*Power(rxj,7LL)) + Power(rxi,2LL)*Power(rxj,14LL)*

            (-3063060LL - 1936935LL*rxj - 408870LL*Power(rxj,2LL) + 11550LL*Power(rxj,3LL) + 

              23100LL*Power(rxj,4LL) + 5082LL*Power(rxj,5LL) + 532LL*Power(rxj,6LL) + 

              24LL*Power(rxj,7LL))))/

      (1260LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj,11LL)*Power(rxi + rxj,11LL))

     ; }
   
  }
  return S;
}

cl_R Slater_2S_5S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = (2011LL*xi)/10240LL

    ; } else {  S = (1LL/r)*(-124540416000LL + 124540416000LL*exp(2LL*rxi) - 224622748350LL*rxi - 

        200164664700LL*Power(rxi,2LL) - 117249207075LL*Power(rxi,3LL) - 

        50639138550LL*Power(rxi,4LL) - 17132415300LL*Power(rxi,5LL) - 

        4704860160LL*Power(rxi,6LL) - 1071195840LL*Power(rxi,7LL) - 

        204478560LL*Power(rxi,8LL) - 32809920LL*Power(rxi,9LL) - 4392960LL*Power(rxi,10LL) - 

        479232LL*Power(rxi,11LL) - 39936LL*Power(rxi,12LL) - 2048LL*Power(rxi,13LL))/

      (1.24540416e11*exp(2LL*rxi))

    ; }
 
  }
  else {
      if (r == 0LL) { S = (xi*xj*(2LL*Power(xi,12LL) + 26LL*Power(xi,11LL)*xj + 156LL*Power(xi,10LL)*Power(xj,2LL) + 

          572LL*Power(xi,9LL)*Power(xj,3LL) + 1430LL*Power(xi,8LL)*Power(xj,4LL) + 

          2574LL*Power(xi,7LL)*Power(xj,5LL) + 3432LL*Power(xi,6LL)*Power(xj,6LL) + 

          3432LL*Power(xi,5LL)*Power(xj,7LL) + 2574LL*Power(xi,4LL)*Power(xj,8LL) + 

          1430LL*Power(xi,3LL)*Power(xj,9LL) + 390LL*Power(xi,2LL)*Power(xj,10LL) + 

          65LL*xi*Power(xj,11LL) + 5LL*Power(xj,12LL)))/(10LL*Power(xi + xj,13LL))

    ; } else { S = (1LL/r)*(28350LL*exp(2LL*(rxi + rxj))*Power(Power(rxi,2LL) - Power(rxj,2LL),13LL) + 

        945LL*exp(2LL*rxj)*Power(rxj,12LL)*

         (-210LL*Power(rxi,16LL) - 10LL*Power(rxi,17LL) + 30LL*Power(rxj,14LL) + 

           45LL*rxi*Power(rxj,14LL) + 39LL*Power(rxi,7LL)*Power(rxj,8LL)*

            (1309LL - 2LL*Power(rxj,2LL)) + 

           858LL*Power(rxi,8LL)*Power(rxj,6LL)*(-305LL + Power(rxj,2LL)) + 

           30LL*Power(rxi,2LL)*Power(rxj,12LL)*(-13LL + Power(rxj,2LL)) - 

           390LL*Power(rxi,4LL)*Power(rxj,10LL)*(-6LL + Power(rxj,2LL)) - 

           143LL*Power(rxi,9LL)*Power(rxj,6LL)*(-153LL + 2LL*Power(rxj,2LL)) + 

           5LL*Power(rxi,3LL)*Power(rxj,12LL)*(-117LL + 2LL*Power(rxj,2LL)) - 

           45LL*Power(rxi,15LL)*(35LL + 2LL*Power(rxj,2LL)) - 

           138LL*Power(rxi,12LL)*Power(rxj,2LL)*(580LL + 13LL*Power(rxj,2LL)) - 

           150LL*Power(rxi,14LL)*(28LL + 17LL*Power(rxj,2LL)) + 

           13LL*Power(rxi,11LL)*Power(rxj,4LL)*(-4071LL + 22LL*Power(rxj,2LL)) + 

           3LL*Power(rxi,13LL)*Power(rxj,2LL)*(-8135LL + 26LL*Power(rxj,2LL)) + 

           3LL*Power(rxi,5LL)*Power(rxj,10LL)*(2171LL + 30LL*Power(rxj,2LL)) + 

           234LL*Power(rxi,10LL)*Power(rxj,4LL)*(-1235LL + 33LL*Power(rxj,2LL)) - 

           78LL*Power(rxi,6LL)*Power(rxj,8LL)*(550LL + 47LL*Power(rxj,2LL))) - 

        2LL*exp(2LL*rxi)*Power(rxi,6LL)*

         (-819LL*Power(rxi,10LL)*Power(rxj,10LL)*

            (22275LL + 39780LL*rxj + 38160LL*Power(rxj,2LL) + 16560LL*Power(rxj,3LL) + 

              9840LL*Power(rxj,4LL) + 3900LL*Power(rxj,5LL) + 816LL*Power(rxj,6LL) + 

              88LL*Power(rxj,7LL) + 4LL*Power(rxj,8LL)) + 

           Power(rxi,20LL)*(14175LL + 25515LL*rxj + 22680LL*Power(rxj,2LL) + 

              13230LL*Power(rxj,3LL) + 5670LL*Power(rxj,4LL) + 1890LL*Power(rxj,5LL) + 

              504LL*Power(rxj,6LL) + 108LL*Power(rxj,7LL) + 18LL*Power(rxj,8LL) + 

              2LL*Power(rxj,9LL)) - Power(rxj,20LL)*

            (16216200LL + 18243225LL*rxj + 9729720LL*Power(rxj,2LL) + 

              3243240LL*Power(rxj,3LL) + 748440LL*Power(rxj,4LL) + 

              124740LL*Power(rxj,5LL) + 15120LL*Power(rxj,6LL) + 1296LL*Power(rxj,7LL) + 

              72LL*Power(rxj,8LL) + 2LL*Power(rxj,9LL)) + 

           18LL*Power(rxi,16LL)*Power(rxj,4LL)*

            (61425LL + 110565LL*rxj + 98280LL*Power(rxj,2LL) + 57330LL*Power(rxj,3LL) + 

              24570LL*Power(rxj,4LL) + 8190LL*Power(rxj,5LL) + 2184LL*Power(rxj,6LL) + 

              496LL*Power(rxj,7LL) + 64LL*Power(rxj,8LL) + 3LL*Power(rxj,9LL)) - 

           18LL*Power(rxi,4LL)*Power(rxj,16LL)*

            (6572475LL - 3161340LL*rxj - 4782960LL*Power(rxj,2LL) - 

              1912365LL*Power(rxj,3LL) - 378105LL*Power(rxj,4LL) - 34125LL*Power(rxj,5LL) + 

              1092LL*Power(rxj,6LL) + 650LL*Power(rxj,7LL) + 71LL*Power(rxj,8LL) + 

              3LL*Power(rxj,9LL)) - 21LL*Power(rxi,8LL)*Power(rxj,12LL)*

            (-1063800LL - 2775735LL*rxj - 862920LL*Power(rxj,2LL) - 

              1132020LL*Power(rxj,3LL) - 698580LL*Power(rxj,4LL) - 

              196920LL*Power(rxj,5LL) - 28992LL*Power(rxj,6LL) - 2064LL*Power(rxj,7LL) - 

              24LL*Power(rxj,8LL) + 4LL*Power(rxj,9LL)) + 

           21LL*Power(rxi,12LL)*Power(rxj,8LL)*

            (482625LL + 868725LL*rxj + 772200LL*Power(rxj,2LL) + 455400LL*Power(rxj,3LL) + 

              178200LL*Power(rxj,4LL) + 72180LL*Power(rxj,5LL) + 19920LL*Power(rxj,6LL) + 

              2952LL*Power(rxj,7LL) + 204LL*Power(rxj,8LL) + 4LL*Power(rxj,9LL)) + 

           6LL*Power(rxi,6LL)*Power(rxj,14LL)*

            (-10357200LL + 5071815LL*rxj - 6463800LL*Power(rxj,2LL) - 

              7151130LL*Power(rxj,3LL) - 2572290LL*Power(rxj,4LL) - 

              468720LL*Power(rxj,5LL) - 42672LL*Power(rxj,6LL) - 648LL*Power(rxj,7LL) + 

              228LL*Power(rxj,8LL) + 16LL*Power(rxj,9LL)) - 

           Power(rxi,18LL)*Power(rxj,2LL)*

            (184275LL + 331695LL*rxj + 294840LL*Power(rxj,2LL) + 171990LL*Power(rxj,3LL) + 

              73710LL*Power(rxj,4LL) + 24570LL*Power(rxj,5LL) + 6552LL*Power(rxj,6LL) + 

              1404LL*Power(rxj,7LL) + 234LL*Power(rxj,8LL) + 16LL*Power(rxj,9LL)) + 

           Power(rxi,2LL)*Power(rxj,18LL)*

            (-133783650LL - 107432325LL*rxj - 35675640LL*Power(rxj,2LL) - 

              5135130LL*Power(rxj,3LL) + 270270LL*Power(rxj,4LL) + 

              270270LL*Power(rxj,5LL) + 57960LL*Power(rxj,6LL) + 6948LL*Power(rxj,7LL) + 

              486LL*Power(rxj,8LL) + 16LL*Power(rxj,9LL)) - 

           6LL*Power(rxi,14LL)*Power(rxj,6LL)*

            (675675LL + 1216215LL*rxj + 1081080LL*Power(rxj,2LL) + 630630LL*Power(rxj,3LL) + 

              270270LL*Power(rxj,4LL) + 88200LL*Power(rxj,5LL) + 26544LL*Power(rxj,6LL) + 

              5160LL*Power(rxj,7LL) + 492LL*Power(rxj,8LL) + 16LL*Power(rxj,9LL))))/

      (28350LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj,13LL)*Power(rxi + rxj,13LL))

     ; }
   
  }
  return S;
}

cl_R Slater_2S_1S(cl_R r,cl_R xi,cl_R xj)
{
  return Slater_1S_2S(r,xj,xi);
}

cl_R Slater_3S_3S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = (793LL*xi)/3072LL

    ; } else {  S = (1LL/r)*(-1437004800LL + 1437004800LL*exp(2LL*rxi) - 2503064025LL*rxi - 

        2132118450LL*Power(rxi,2LL) - 1180664100LL*Power(rxi,3LL) - 

        476506800LL*Power(rxi,4LL) - 148856400LL*Power(rxi,5LL) - 

        37255680LL*Power(rxi,6LL) - 7603200LL*Power(rxi,7LL) - 1267200LL*Power(rxi,8LL) - 

        168960LL*Power(rxi,9LL) - 16896LL*Power(rxi,10LL) - 1024LL*Power(rxi,11LL))/

      (1.4370048e9*exp(2LL*rxi))

    ; }
 
  }
  else {
      if (r == 0LL) { S = (xi*xj*(Power(xi,10LL) + 11LL*Power(xi,9LL)*xj + 55LL*Power(xi,8LL)*Power(xj,2LL) + 

          165LL*Power(xi,7LL)*Power(xj,3LL) + 330LL*Power(xi,6LL)*Power(xj,4LL) + 

          462LL*Power(xi,5LL)*Power(xj,5LL) + 330LL*Power(xi,4LL)*Power(xj,6LL) + 

          165LL*Power(xi,3LL)*Power(xj,7LL) + 55LL*Power(xi,2LL)*Power(xj,8LL) + 

          11LL*xi*Power(xj,9LL) + Power(xj,10LL)))/(3LL*Power(xi + xj,11LL))

    ; } else { S = (1LL/r)*(135LL*exp(2LL*(rxi + rxj))*Power(Power(rxi,2LL) - Power(rxj,2LL),11LL) + 

        exp(2LL*rxj)*Power(rxj,8LL)*

         (-150LL*Power(rxi,18LL) - 6LL*Power(rxi,19LL) + 135LL*Power(rxj,14LL) + 

           225LL*rxi*Power(rxj,14LL) + 10LL*Power(rxi,17LL)*(-165LL + Power(rxj,2LL)) - 

           30LL*Power(rxi,16LL)*(330LL + Power(rxj,2LL)) + 

           45LL*Power(rxi,3LL)*Power(rxj,12LL)*(-55LL + 2LL*Power(rxj,2LL)) + 

           45LL*Power(rxi,2LL)*Power(rxj,12LL)*(-33LL + 4LL*Power(rxj,2LL)) + 

           Power(rxi,9LL)*Power(rxj,6LL)*

            (234135LL - 4950LL*Power(rxj,2LL) - 34LL*Power(rxj,4LL)) - 

           5LL*Power(rxi,7LL)*Power(rxj,8LL)*

            (6237LL - 1242LL*Power(rxj,2LL) + 2LL*Power(rxj,4LL)) + 

           3LL*Power(rxi,5LL)*Power(rxj,10LL)*

            (4125LL - 330LL*Power(rxj,2LL) + 2LL*Power(rxj,4LL)) + 

           15LL*Power(rxi,4LL)*Power(rxj,10LL)*

            (495LL - 132LL*Power(rxj,2LL) + 2LL*Power(rxj,4LL)) - 

           165LL*Power(rxi,6LL)*Power(rxj,8LL)*

            (135LL - 60LL*Power(rxj,2LL) + 2LL*Power(rxj,4LL)) - 

           5LL*Power(rxi,13LL)*Power(rxj,2LL)*

            (43875LL - 3438LL*Power(rxj,2LL) + 22LL*Power(rxj,4LL)) + 

           5LL*Power(rxi,11LL)*Power(rxj,4LL)*

            (7695LL - 2442LL*Power(rxj,2LL) + 22LL*Power(rxj,4LL)) + 

           15LL*Power(rxi,8LL)*Power(rxj,6LL)*

            (-33LL - 3564LL*Power(rxj,2LL) + 26LL*Power(rxj,4LL)) + 

           Power(rxi,15LL)*(-32175LL - 3690LL*Power(rxj,2LL) + 34LL*Power(rxj,4LL)) + 

           15LL*Power(rxi,10LL)*Power(rxj,4LL)*

            (-32277LL + 1364LL*Power(rxj,2LL) + 66LL*Power(rxj,4LL)) + 

           15LL*Power(rxi,14LL)*(-3003LL - 2932LL*Power(rxj,2LL) + 94LL*Power(rxj,4LL)) - 

           15LL*Power(rxi,12LL)*Power(rxj,2LL)*

            (28119LL - 5252LL*Power(rxj,2LL) + 154LL*Power(rxj,4LL))) + 

        exp(2LL*rxi)*Power(rxi,8LL)*

         (-5LL*Power(rxi,2LL)*Power(rxj,12LL)*

            (-84357LL - 43875LL*rxj - 8796LL*Power(rxj,2LL) - 738LL*Power(rxj,3LL) - 

              6LL*Power(rxj,4LL) + 2LL*Power(rxj,5LL)) - 

           3LL*Power(rxi,14LL)*(45LL + 75LL*rxj + 60LL*Power(rxj,2LL) + 30LL*Power(rxj,3LL) + 

              10LL*Power(rxj,4LL) + 2LL*Power(rxj,5LL)) - 

           55LL*Power(rxi,8LL)*Power(rxj,6LL)*

            (-405LL - 567LL*rxj - 972LL*Power(rxj,2LL) - 90LL*Power(rxj,3LL) + 

              18LL*Power(rxj,4LL) + 2LL*Power(rxj,5LL)) + 

           55LL*Power(rxi,6LL)*Power(rxj,8LL)*

            (9LL - 4257LL*rxj - 372LL*Power(rxj,2LL) + 222LL*Power(rxj,3LL) + 

              42LL*Power(rxj,4LL) + 2LL*Power(rxj,5LL)) + 

           3LL*Power(rxj,14LL)*(15015LL + 10725LL*rxj + 3300LL*Power(rxj,2LL) + 

              550LL*Power(rxj,3LL) + 50LL*Power(rxj,4LL) + 2LL*Power(rxj,5LL)) + 

           5LL*Power(rxi,12LL)*Power(rxj,2LL)*

            (297LL + 495LL*rxj + 396LL*Power(rxj,2LL) + 198LL*Power(rxj,3LL) + 

              66LL*Power(rxj,4LL) + 2LL*Power(rxj,5LL)) + 

           Power(rxi,10LL)*Power(rxj,4LL)*

            (-7425LL - 12375LL*rxj - 9900LL*Power(rxj,2LL) - 6210LL*Power(rxj,3LL) - 

              390LL*Power(rxj,4LL) + 34LL*Power(rxj,5LL)) - 

           Power(rxi,4LL)*Power(rxj,10LL)*

            (-484155LL + 38475LL*rxj + 78780LL*Power(rxj,2LL) + 17190LL*Power(rxj,3LL) + 

              1410LL*Power(rxj,4LL) + 34LL*Power(rxj,5LL))))/

      (135LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj,11LL)*Power(rxi + rxj,11LL))

     ; }
   
  }
  return S;
}

cl_R Slater_3S_4S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = (1363LL*xi)/6144LL

    ; } else {  S = (1LL/r)*(-74724249600LL + 74724249600LL*exp(2LL*rxi) - 132871488750LL*rxi - 

        116294478300LL*Power(rxi,2LL) - 66678987375LL*Power(rxi,3LL) - 

        28114836750LL*Power(rxi,4LL) - 9274044780LL*Power(rxi,5LL) - 

        2484321840LL*Power(rxi,6LL) - 553204080LL*Power(rxi,7LL) - 

        103783680LL*Power(rxi,8LL) - 16473600LL*Power(rxi,9LL) - 2196480LL*Power(rxi,10LL) - 

        239616LL*Power(rxi,11LL) - 19968LL*Power(rxi,12LL) - 1024LL*Power(rxi,13LL))/

      (7.47242496e10*exp(2LL*rxi))

    ; }
 
  }
  else {
      if (r == 0LL) { S = (xi*xj*(3LL*Power(xi,12LL) + 39LL*Power(xi,11LL)*xj + 234LL*Power(xi,10LL)*Power(xj,2LL) + 

          858LL*Power(xi,9LL)*Power(xj,3LL) + 2145LL*Power(xi,8LL)*Power(xj,4LL) + 

          3861LL*Power(xi,7LL)*Power(xj,5LL) + 5148LL*Power(xi,6LL)*Power(xj,6LL) + 

          5148LL*Power(xi,5LL)*Power(xj,7LL) + 2860LL*Power(xi,4LL)*Power(xj,8LL) + 

          1144LL*Power(xi,3LL)*Power(xj,9LL) + 312LL*Power(xi,2LL)*Power(xj,10LL) + 

          52LL*xi*Power(xj,11LL) + 4LL*Power(xj,12LL)))/(12LL*Power(xi + xj,13LL))

    ; } else { S = (1LL/r)*(3780LL*exp(2LL*(rxi + rxj))*Power(Power(rxi,2LL) - Power(rxj,2LL),13LL) + 

        84LL*exp(2LL*rxj)*Power(rxj,10LL)*

         (-60LL*Power(rxi,20LL) - 2LL*Power(rxi,21LL) + 45LL*Power(rxj,16LL) + 

           75LL*rxi*Power(rxj,16LL) - 4LL*Power(rxi,19LL)*(195LL + Power(rxj,2LL)) + 

           15LL*Power(rxi,3LL)*Power(rxj,14LL)*(-65LL + 2LL*Power(rxj,2LL)) + 

           15LL*Power(rxi,2LL)*Power(rxj,14LL)*(-39LL + 4LL*Power(rxj,2LL)) - 

           30LL*Power(rxi,18LL)*(182LL + 9LL*Power(rxj,2LL)) + 

           30LL*Power(rxi,13LL)*Power(rxj,4LL)*(-13047LL + 377LL*Power(rxj,2LL)) + 

           2LL*Power(rxi,5LL)*Power(rxj,12LL)*

            (2925LL - 195LL*Power(rxj,2LL) + Power(rxj,4LL)) + 

           10LL*Power(rxi,4LL)*Power(rxj,12LL)*

            (351LL - 78LL*Power(rxj,2LL) + Power(rxj,4LL)) - 

           130LL*Power(rxi,6LL)*Power(rxj,10LL)*

            (99LL - 36LL*Power(rxj,2LL) + Power(rxj,4LL)) + 

           13LL*Power(rxi,11LL)*Power(rxj,6LL)*

            (30735LL - 1650LL*Power(rxj,2LL) + 4LL*Power(rxj,4LL)) + 

           Power(rxi,7LL)*Power(rxj,10LL)*

            (-15015LL + 3330LL*Power(rxj,2LL) + 4LL*Power(rxj,4LL)) + 

           210LL*Power(rxi,16LL)*(-156LL - 262LL*Power(rxj,2LL) + 5LL*Power(rxj,4LL)) - 

           6LL*Power(rxi,9LL)*Power(rxj,8LL)*

            (-48620LL - 715LL*Power(rxj,2LL) + 6LL*Power(rxj,4LL)) + 

           3LL*Power(rxi,17LL)*(-6825LL - 1870LL*Power(rxj,2LL) + 12LL*Power(rxj,4LL)) - 

           30LL*Power(rxi,14LL)*Power(rxj,2LL)*

            (17934LL - 12LL*Power(rxj,2LL) + 13LL*Power(rxj,4LL)) - 

           15LL*Power(rxi,8LL)*Power(rxj,8LL)*

            (2145LL + 2860LL*Power(rxj,2LL) + 14LL*Power(rxj,4LL)) + 

           65LL*Power(rxi,10LL)*Power(rxj,6LL)*

            (-13725LL - 792LL*Power(rxj,2LL) + 22LL*Power(rxj,4LL)) - 

           10LL*Power(rxi,12LL)*Power(rxj,4LL)*

            (153630LL - 15054LL*Power(rxj,2LL) + 143LL*Power(rxj,4LL)) + 

           Power(rxi,15LL)*(-269325LL*Power(rxj,2LL) + 9270LL*Power(rxj,4LL) - 

              52LL*Power(rxj,6LL))) + exp(2LL*rxi)*Power(rxi,8LL)*

         (Power(rxi,2LL)*Power(rxj,16LL)*

            (70073640LL + 47669895LL*rxj + 13931190LL*Power(rxj,2LL) + 

              2170350LL*Power(rxj,3LL) + 169260LL*Power(rxj,4LL) + 1638LL*Power(rxj,5LL) - 

              756LL*Power(rxj,6LL) - 44LL*Power(rxj,7LL)) + 

           364LL*Power(rxi,10LL)*Power(rxj,8LL)*

            (-7425LL - 13860LL*rxj - 5940LL*Power(rxj,2LL) - 11880LL*Power(rxj,3LL) - 

              2640LL*Power(rxj,4LL) - 45LL*Power(rxj,5LL) + 30LL*Power(rxj,6LL) + 

              2LL*Power(rxj,7LL)) - 364LL*Power(rxi,8LL)*Power(rxj,10LL)*

            (-20925LL + 18270LL*rxj - 58320LL*Power(rxj,2LL) - 17730LL*Power(rxj,3LL) - 

              300LL*Power(rxj,4LL) + 423LL*Power(rxj,5LL) + 54LL*Power(rxj,6LL) + 

              2LL*Power(rxj,7LL)) - 3LL*Power(rxi,18LL)*

            (1260LL + 2205LL*rxj + 1890LL*Power(rxj,2LL) + 1050LL*Power(rxj,3LL) + 

              420LL*Power(rxj,4LL) + 126LL*Power(rxj,5LL) + 28LL*Power(rxj,6LL) + 

              4LL*Power(rxj,7LL)) + 3LL*Power(rxj,18LL)*

            (1801800LL + 1576575LL*rxj + 630630LL*Power(rxj,2LL) + 

              150150LL*Power(rxj,3LL) + 23100LL*Power(rxj,4LL) + 2310LL*Power(rxj,5LL) + 

              140LL*Power(rxj,6LL) + 4LL*Power(rxj,7LL)) + 

           2LL*Power(rxi,14LL)*Power(rxj,4LL)*

            (-147420LL - 257985LL*rxj - 221130LL*Power(rxj,2LL) - 122850LL*Power(rxj,3LL) - 

              49140LL*Power(rxj,4LL) - 17388LL*Power(rxj,5LL) - 1512LL*Power(rxj,6LL) + 

              8LL*Power(rxj,7LL)) - 42LL*Power(rxi,12LL)*Power(rxj,6LL)*

            (-25740LL - 45045LL*rxj - 38610LL*Power(rxj,2LL) - 19470LL*Power(rxj,3LL) - 

              12540LL*Power(rxj,4LL) - 1836LL*Power(rxj,5LL) - 8LL*Power(rxj,6LL) + 

              8LL*Power(rxj,7LL)) + 42LL*Power(rxi,6LL)*Power(rxj,12LL)*

            (921600LL - 1640835LL*rxj - 546030LL*Power(rxj,2LL) + 20730LL*Power(rxj,3LL) + 

              30180LL*Power(rxj,4LL) + 5028LL*Power(rxj,5LL) + 344LL*Power(rxj,6LL) + 

              8LL*Power(rxj,7LL)) - 2LL*Power(rxi,4LL)*Power(rxj,14LL)*

            (-67767840LL - 13377735LL*rxj + 6601770LL*Power(rxj,2LL) + 

              3115350LL*Power(rxj,3LL) + 548940LL*Power(rxj,4LL) + 48132LL*Power(rxj,5LL) + 

              1848LL*Power(rxj,6LL) + 8LL*Power(rxj,7LL)) + 

           Power(rxi,16LL)*Power(rxj,2LL)*

            (49140LL + 85995LL*rxj + 73710LL*Power(rxj,2LL) + 40950LL*Power(rxj,3LL) + 

              16380LL*Power(rxj,4LL) + 4914LL*Power(rxj,5LL) + 1092LL*Power(rxj,6LL) + 

              44LL*Power(rxj,7LL))))/

      (3780LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj,13LL)*Power(rxi + rxj,13LL))

     ; }
   
  }
  return S;
}

cl_R Slater_3S_5S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = (31059LL*xi)/163840LL

    ; } else {  S = (1LL/r)*(-313841848320000LL + 313841848320000LL*exp(2LL*rxi) - 568188982486125LL*rxi - 

        508694268332250LL*Power(rxi,2LL) - 299892470377500LL*Power(rxi,3LL) - 

        130753815192000LL*Power(rxi,4LL) - 44881155118800LL*Power(rxi,5LL) - 

        12601803614400LL*Power(rxi,6LL) - 2967953788800LL*Power(rxi,7LL) - 

        596237241600LL*Power(rxi,8LL) - 103264761600LL*Power(rxi,9LL) - 

        15498362880LL*Power(rxi,10LL) - 2012774400LL*Power(rxi,11LL) - 

        223641600LL*Power(rxi,12LL) - 20643840LL*Power(rxi,13LL) - 

        1474560LL*Power(rxi,14LL) - 65536LL*Power(rxi,15LL))/

      (3.1384184832e14*exp(2LL*rxi))

    ; }
 
  }
  else {
      if (r == 0LL) { S = (xi*xj*(3LL*Power(xi,14LL) + 45LL*Power(xi,13LL)*xj + 315LL*Power(xi,12LL)*Power(xj,2LL) + 

          1365LL*Power(xi,11LL)*Power(xj,3LL) + 4095LL*Power(xi,10LL)*Power(xj,4LL) + 

          9009LL*Power(xi,9LL)*Power(xj,5LL) + 15015LL*Power(xi,8LL)*Power(xj,6LL) + 

          19305LL*Power(xi,7LL)*Power(xj,7LL) + 19305LL*Power(xi,6LL)*Power(xj,8LL) + 

          15015LL*Power(xi,5LL)*Power(xj,9LL) + 6825LL*Power(xi,4LL)*Power(xj,10LL) + 

          2275LL*Power(xi,3LL)*Power(xj,11LL) + 525LL*Power(xi,2LL)*Power(xj,12LL) + 

          75LL*xi*Power(xj,13LL) + 5LL*Power(xj,14LL)))/(15LL*Power(xi + xj,15LL))

    ; } else { S = (1LL/r)*(42525LL*exp(2LL*(rxi + rxj))*Power(Power(rxi,2LL) - Power(rxj,2LL),15LL) + 

        189LL*exp(2LL*rxj)*Power(rxj,12LL)*

         (-350LL*Power(rxi,22LL) - 10LL*Power(rxi,23LL) + 225LL*Power(rxj,18LL) + 

           375LL*rxi*Power(rxj,18LL) - 70LL*Power(rxi,21LL)*(75LL + Power(rxj,2LL)) + 

           75LL*Power(rxi,3LL)*Power(rxj,16LL)*(-75LL + 2LL*Power(rxj,2LL)) + 

           75LL*Power(rxi,2LL)*Power(rxj,16LL)*(-45LL + 4LL*Power(rxj,2LL)) - 

           50LL*Power(rxi,20LL)*(840LL + 71LL*Power(rxj,2LL)) + 

           Power(rxi,9LL)*Power(rxj,10LL)*

            (4694625LL + 124800LL*Power(rxj,2LL) - 248LL*Power(rxj,4LL)) + 

           20LL*Power(rxi,17LL)*Power(rxj,2LL)*

            (-185895LL - 948LL*Power(rxj,2LL) + 2LL*Power(rxj,4LL)) + 

           5LL*Power(rxi,5LL)*Power(rxj,14LL)*

            (7875LL - 450LL*Power(rxj,2LL) + 2LL*Power(rxj,4LL)) + 

           25LL*Power(rxi,4LL)*Power(rxj,14LL)*

            (945LL - 180LL*Power(rxj,2LL) + 2LL*Power(rxj,4LL)) - 

           375LL*Power(rxi,6LL)*Power(rxj,12LL)*

            (273LL - 84LL*Power(rxj,2LL) + 2LL*Power(rxj,4LL)) - 

           5LL*Power(rxi,11LL)*Power(rxj,8LL)*

            (-2803125LL + 49140LL*Power(rxj,2LL) + 8LL*Power(rxj,4LL)) + 

           5LL*Power(rxi,7LL)*Power(rxj,12LL)*

            (-16965LL + 5152LL*Power(rxj,2LL) + 14LL*Power(rxj,4LL)) + 

           325LL*Power(rxi,10LL)*Power(rxj,8LL)*

            (-60117LL - 5340LL*Power(rxj,2LL) + 40LL*Power(rxj,4LL)) - 

           15LL*Power(rxi,15LL)*Power(rxj,4LL)*

            (845085LL - 22960LL*Power(rxj,2LL) + 52LL*Power(rxj,4LL)) + 

           15LL*Power(rxi,13LL)*Power(rxj,6LL)*

            (-139125LL - 10140LL*Power(rxj,2LL) + 52LL*Power(rxj,4LL)) + 

           75LL*Power(rxi,12LL)*Power(rxj,6LL)*

            (-729687LL + 25532LL*Power(rxj,2LL) + 52LL*Power(rxj,4LL)) + 

           60LL*Power(rxi,18LL)*(-5355LL - 11940LL*Power(rxj,2LL) + 86LL*Power(rxj,4LL)) + 

           2LL*Power(rxi,19LL)*(-89250LL - 35425LL*Power(rxj,2LL) + 124LL*Power(rxj,4LL)) + 

           100LL*Power(rxi,16LL)*Power(rxj,2LL)*

            (-79713LL - 13311LL*Power(rxj,2LL) + 146LL*Power(rxj,4LL)) - 

           5LL*Power(rxi,8LL)*Power(rxj,10LL)*

            (157365LL + 95940LL*Power(rxj,2LL) + 952LL*Power(rxj,4LL)) - 

           15LL*Power(rxi,14LL)*Power(rxj,4LL)*

            (2638467LL - 157500LL*Power(rxj,2LL) + 1820LL*Power(rxj,4LL))) - 

        exp(2LL*rxi)*Power(rxi,8LL)*

         (3LL*Power(rxi,14LL)*Power(rxj,8LL)*

            (19348875LL + 34827975LL*rxj + 30958200LL*Power(rxj,2LL) + 

              18689580LL*Power(rxj,3LL) + 5847660LL*Power(rxj,4LL) + 

              3723300LL*Power(rxj,5LL) + 845040LL*Power(rxj,6LL) + 58680LL*Power(rxj,7LL) - 

              1548LL*Power(rxj,8LL) - 236LL*Power(rxj,9LL)) - 

           42LL*Power(rxi,4LL)*Power(rxj,18LL)*

            (251336925LL + 104824125LL*rxj + 340200LL*Power(rxj,2LL) - 

              9122085LL*Power(rxj,3LL) - 2798145LL*Power(rxj,4LL) - 

              433755LL*Power(rxj,5LL) - 39060LL*Power(rxj,6LL) - 1890LL*Power(rxj,7LL) - 

              27LL*Power(rxj,8LL) + Power(rxj,9LL)) - 

           6LL*Power(rxj,22LL)*(34459425LL + 34459425LL*rxj + 16216200LL*Power(rxj,2LL) + 

              4729725LL*Power(rxj,3LL) + 945945LL*Power(rxj,4LL) + 

              135135LL*Power(rxj,5LL) + 13860LL*Power(rxj,6LL) + 990LL*Power(rxj,7LL) + 

              45LL*Power(rxj,8LL) + Power(rxj,9LL)) + 

           3LL*Power(rxi,22LL)*(14175LL + 25515LL*rxj + 22680LL*Power(rxj,2LL) + 

              13230LL*Power(rxj,3LL) + 5670LL*Power(rxj,4LL) + 1890LL*Power(rxj,5LL) + 

              504LL*Power(rxj,6LL) + 108LL*Power(rxj,7LL) + 18LL*Power(rxj,8LL) + 

              2LL*Power(rxj,9LL)) + 21LL*Power(rxi,18LL)*Power(rxj,4LL)*

            (212625LL + 382725LL*rxj + 340200LL*Power(rxj,2LL) + 198450LL*Power(rxj,3LL) + 

              85050LL*Power(rxj,4LL) + 28350LL*Power(rxj,5LL) + 7560LL*Power(rxj,6LL) + 

              1836LL*Power(rxj,7LL) + 162LL*Power(rxj,8LL) + 2LL*Power(rxj,9LL)) - 

           54LL*Power(rxi,6LL)*Power(rxj,16LL)*

            (133451955LL - 73700865LL*rxj - 54096840LL*Power(rxj,2LL) - 

              8306235LL*Power(rxj,3LL) + 966945LL*Power(rxj,4LL) + 

              516747LL*Power(rxj,5LL) + 80724LL*Power(rxj,6LL) + 6434LL*Power(rxj,7LL) + 

              251LL*Power(rxj,8LL) + 3LL*Power(rxj,9LL)) + 

           315LL*Power(rxi,12LL)*Power(rxj,10LL)*

            (-405405LL - 710073LL*rxj - 805896LL*Power(rxj,2LL) - 101556LL*Power(rxj,3LL) - 

              258804LL*Power(rxj,4LL) - 90972LL*Power(rxj,5LL) - 9744LL*Power(rxj,6LL) + 

              120LL*Power(rxj,7LL) + 84LL*Power(rxj,8LL) + 4LL*Power(rxj,9LL)) - 

           315LL*Power(rxi,10LL)*Power(rxj,12LL)*

            (-482895LL - 2656395LL*rxj + 1186920LL*Power(rxj,2LL) - 

              1155420LL*Power(rxj,3LL) - 643356LL*Power(rxj,4LL) - 93492LL*Power(rxj,5LL) + 

              336LL*Power(rxj,6LL) + 1368LL*Power(rxj,7LL) + 132LL*Power(rxj,8LL) + 

              4LL*Power(rxj,9LL)) + 27LL*Power(rxi,16LL)*Power(rxj,6LL)*

            (-716625LL - 1289925LL*rxj - 1146600LL*Power(rxj,2LL) - 

              668850LL*Power(rxj,3LL) - 286650LL*Power(rxj,4LL) - 90006LL*Power(rxj,5LL) - 

              32872LL*Power(rxj,6LL) - 4812LL*Power(rxj,7LL) - 178LL*Power(rxj,8LL) + 

              6LL*Power(rxj,9LL)) + 2LL*Power(rxi,2LL)*Power(rxj,20LL)*

            (-1782492075LL - 1449175455LL*rxj - 533365560LL*Power(rxj,2LL) - 

              114631335LL*Power(rxj,3LL) - 15221115LL*Power(rxj,4LL) - 

              1142505LL*Power(rxj,5LL) - 18396LL*Power(rxj,6LL) + 5238LL*Power(rxj,7LL) + 

              513LL*Power(rxj,8LL) + 17LL*Power(rxj,9LL)) - 

           Power(rxi,20LL)*Power(rxj,2LL)*

            (637875LL + 1148175LL*rxj + 1020600LL*Power(rxj,2LL) + 

              595350LL*Power(rxj,3LL) + 255150LL*Power(rxj,4LL) + 85050LL*Power(rxj,5LL) + 

              22680LL*Power(rxj,6LL) + 4860LL*Power(rxj,7LL) + 810LL*Power(rxj,8LL) + 

              34LL*Power(rxj,9LL)) + 3LL*Power(rxi,8LL)*Power(rxj,14LL)*

            (-593408025LL + 946053675LL*rxj - 394427880LL*Power(rxj,2LL) - 

              315870660LL*Power(rxj,3LL) - 53891460LL*Power(rxj,4LL) + 

              910980LL*Power(rxj,5LL) + 1409520LL*Power(rxj,6LL) + 192168LL*Power(rxj,7LL) + 

              11196LL*Power(rxj,8LL) + 236LL*Power(rxj,9LL))))/

      (42525LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj,15LL)*Power(rxi + rxj,15LL))

     ; }
   
  }
  return S;
}

cl_R Slater_3S_1S(cl_R r,cl_R xi,cl_R xj)
{
  return Slater_1S_3S(r,xj,xi);
}

cl_R Slater_3S_2S(cl_R r,cl_R xi,cl_R xj)
{
  return Slater_2S_3S(r,xj,xi);
}

cl_R Slater_4S_4S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = (26333LL*xi)/131072LL

    ; } else {  S = (1LL/r)*(-83691159552000LL + 83691159552000LL*exp(2LL*rxi) - 150568359566625LL*rxi - 

        133754400029250LL*Power(rxi,2LL) - 78142908343500LL*Power(rxi,3LL) - 

        33740723016000LL*Power(rxi,4LL) - 11470756096800LL*Power(rxi,5LL) - 

        3193358968800LL*Power(rxi,6LL) - 747112766400LL*Power(rxi,7LL) - 

        149448499200LL*Power(rxi,8LL) - 25830604800LL*Power(rxi,9LL) - 

        3874590720LL*Power(rxi,10LL) - 503193600LL*Power(rxi,11LL) - 

        55910400LL*Power(rxi,12LL) - 5160960LL*Power(rxi,13LL) - 368640LL*Power(rxi,14LL) - 

        16384LL*Power(rxi,15LL))/(8.3691159552e13*exp(2LL*rxi))

    ; }
 
  }
  else {
      if (r == 0LL) { S = (xi*xj*(Power(xi,14LL) + 15LL*Power(xi,13LL)*xj + 105LL*Power(xi,12LL)*Power(xj,2LL) + 

          455LL*Power(xi,11LL)*Power(xj,3LL) + 1365LL*Power(xi,10LL)*Power(xj,4LL) + 

          3003LL*Power(xi,9LL)*Power(xj,5LL) + 5005LL*Power(xi,8LL)*Power(xj,6LL) + 

          6435LL*Power(xi,7LL)*Power(xj,7LL) + 5005LL*Power(xi,6LL)*Power(xj,8LL) + 

          3003LL*Power(xi,5LL)*Power(xj,9LL) + 1365LL*Power(xi,4LL)*Power(xj,10LL) + 

          455LL*Power(xi,3LL)*Power(xj,11LL) + 105LL*Power(xi,2LL)*Power(xj,12LL) + 

          15LL*xi*Power(xj,13LL) + Power(xj,14LL)))/(4LL*Power(xi + xj,15LL))

    ; } else { S = (1LL/r)*(1260LL*exp(2LL*(rxi + rxj))*Power(Power(rxi,2LL) - Power(rxj,2LL),15LL) + 

        exp(2LL*rxj)*Power(rxj,10LL)*

         (-3276LL*Power(rxi,25LL) - 168LL*Power(rxi,26LL) - 4LL*Power(rxi,27LL) + 

           1260LL*Power(rxj,20LL) + 2205LL*rxi*Power(rxj,20LL) + 

           1890LL*Power(rxi,2LL)*Power(rxj,18LL)*(-10LL + Power(rxj,2LL)) - 

           420LL*Power(rxi,24LL)*(91LL + Power(rxj,2LL)) + 

           525LL*Power(rxi,3LL)*Power(rxj,18LL)*(-63LL + 2LL*Power(rxj,2LL)) + 

           42LL*Power(rxi,23LL)*(-6825LL - 405LL*Power(rxj,2LL) + 2LL*Power(rxj,4LL)) + 

           63LL*Power(rxi,5LL)*Power(rxj,16LL)*

            (3675LL - 250LL*Power(rxj,2LL) + 2LL*Power(rxj,4LL)) + 

           210LL*Power(rxi,4LL)*Power(rxj,16LL)*

            (630LL - 135LL*Power(rxj,2LL) + 2LL*Power(rxj,4LL)) + 

           252LL*Power(rxi,22LL)*(-5460LL - 1225LL*Power(rxj,2LL) + 17LL*Power(rxj,4LL)) - 

           1260LL*Power(rxi,17LL)*Power(rxj,4LL)*

            (141729LL - 10145LL*Power(rxj,2LL) + 116LL*Power(rxj,4LL)) + 

           21LL*Power(rxi,9LL)*Power(rxj,12LL)*

            (164775LL - 18460LL*Power(rxj,2LL) + 828LL*Power(rxj,4LL)) + 

           14LL*Power(rxi,6LL)*Power(rxj,14LL)*

            (-40950LL + 14175LL*Power(rxj,2LL) - 450LL*Power(rxj,4LL) + 2LL*Power(rxj,6LL)) - 

           210LL*Power(rxi,8LL)*Power(rxj,12LL)*

            (-8190LL + 4095LL*Power(rxj,2LL) - 210LL*Power(rxj,4LL) + 2LL*Power(rxj,6LL)) + 

           42LL*Power(rxi,10LL)*Power(rxj,10LL)*

            (-209430LL - 2925LL*Power(rxj,2LL) - 8840LL*Power(rxj,4LL) + 4LL*Power(rxj,6LL)) 

    + Power(rxi,7LL)*Power(rxj,14LL)*(-1003275LL + 110250LL*Power(rxj,2LL) - 

              1890LL*Power(rxj,4LL) + 4LL*Power(rxj,6LL)) - 

           21LL*Power(rxi,11LL)*Power(rxj,10LL)*

            (-1033695LL - 218400LL*Power(rxj,2LL) + 552LL*Power(rxj,4LL) + 

              4LL*Power(rxj,6LL)) + 280LL*Power(rxi,18LL)*Power(rxj,2LL)*

            (-385560LL - 73953LL*Power(rxj,2LL) + 2370LL*Power(rxj,4LL) + 4LL*Power(rxj,6LL)) 

    - 35LL*Power(rxi,15LL)*Power(rxj,6LL)*

            (-1565613LL + 359520LL*Power(rxj,2LL) - 7020LL*Power(rxj,4LL) + 

              8LL*Power(rxj,6LL)) + 14LL*Power(rxi,19LL)*Power(rxj,2LL)*

            (-4980150LL + 126765LL*Power(rxj,2LL) - 3852LL*Power(rxj,4LL) + 

              20LL*Power(rxj,6LL)) - 630LL*Power(rxi,14LL)*Power(rxj,6LL)*

            (708714LL - 14385LL*Power(rxj,2LL) - 2340LL*Power(rxj,4LL) + 20LL*Power(rxj,6LL)) 

    + 210LL*Power(rxi,16LL)*Power(rxj,4LL)*

            (-2087532LL + 328491LL*Power(rxj,2LL) - 11740LL*Power(rxj,4LL) + 

              52LL*Power(rxj,6LL)) - 84LL*Power(rxi,20LL)*

            (59670LL + 236250LL*Power(rxj,2LL) - 8745LL*Power(rxj,4LL) + 92LL*Power(rxj,6LL)) 

    - 2LL*Power(rxi,21LL)*(1949220LL + 1598625LL*Power(rxj,2LL) - 41391LL*Power(rxj,4LL) + 

              128LL*Power(rxj,6LL)) + Power(rxi,13LL)*Power(rxj,8LL)*

            (173037375LL - 2784600LL*Power(rxj,2LL) - 112140LL*Power(rxj,4LL) + 

              256LL*Power(rxj,6LL)) + 14LL*Power(rxi,12LL)*Power(rxj,8LL)*

            (-7260750LL - 2521935LL*Power(rxj,2LL) + 19500LL*Power(rxj,4LL) + 

              344LL*Power(rxj,6LL))) + 

        exp(2LL*rxi)*Power(rxi,10LL)*

         (210LL*Power(rxi,2LL)*Power(rxj,18LL)*

            (514080LL + 332010LL*rxj + 94500LL*Power(rxj,2LL) + 15225LL*Power(rxj,3LL) + 

              1470LL*Power(rxj,4LL) + 81LL*Power(rxj,5LL) + 2LL*Power(rxj,6LL)) + 

           105LL*Power(rxi,18LL)*Power(rxj,2LL)*

            (180LL + 315LL*rxj + 270LL*Power(rxj,2LL) + 150LL*Power(rxj,3LL) + 

              60LL*Power(rxj,4LL) + 18LL*Power(rxj,5LL) + 4LL*Power(rxj,6LL)) - 

           1365LL*Power(rxi,10LL)*Power(rxj,10LL)*

            (-6444LL + 15903LL*rxj - 25866LL*Power(rxj,2LL) - 2040LL*Power(rxj,3LL) + 

              1080LL*Power(rxj,4LL) + 180LL*Power(rxj,5LL) + 8LL*Power(rxj,6LL)) + 

           Power(rxi,14LL)*Power(rxj,6LL)*

            (573300LL + 1003275LL*rxj + 859950LL*Power(rxj,2LL) + 387660LL*Power(rxj,3LL) + 

              371280LL*Power(rxj,4LL) + 11592LL*Power(rxj,5LL) - 4816LL*Power(rxj,6LL) - 

              256LL*Power(rxj,7LL)) + 2LL*Power(rxj,20LL)*

            (2506140LL + 1949220LL*rxj + 687960LL*Power(rxj,2LL) + 

              143325LL*Power(rxj,3LL) + 19110LL*Power(rxj,4LL) + 1638LL*Power(rxj,5LL) + 

              84LL*Power(rxj,6LL) + 2LL*Power(rxj,7LL)) - 

           42LL*Power(rxi,4LL)*Power(rxj,16LL)*

            (-10437660LL - 4251870LL*rxj - 493020LL*Power(rxj,2LL) + 

              42255LL*Power(rxj,3LL) + 17490LL*Power(rxj,4LL) + 1971LL*Power(rxj,5LL) + 

              102LL*Power(rxj,6LL) + 2LL*Power(rxj,7LL)) + 

           21LL*Power(rxi,16LL)*Power(rxj,4LL)*

            (-6300LL - 11025LL*rxj - 9450LL*Power(rxj,2LL) - 5250LL*Power(rxj,3LL) - 

              2100LL*Power(rxj,4LL) - 828LL*Power(rxj,5LL) - 8LL*Power(rxj,6LL) + 

              4LL*Power(rxj,7LL)) - Power(rxi,20LL)*

            (1260LL + 2205LL*rxj + 1890LL*Power(rxj,2LL) + 1050LL*Power(rxj,3LL) + 

              420LL*Power(rxj,4LL) + 126LL*Power(rxj,5LL) + 28LL*Power(rxj,6LL) + 

              4LL*Power(rxj,7LL)) - 35LL*Power(rxi,8LL)*Power(rxj,12LL)*

            (-2904300LL + 4943925LL*rxj + 258930LL*Power(rxj,2LL) - 

              359520LL*Power(rxj,3LL) - 70440LL*Power(rxj,4LL) - 4176LL*Power(rxj,5LL) + 

              32LL*Power(rxj,6LL) + 8LL*Power(rxj,7LL)) + 

           35LL*Power(rxi,12LL)*Power(rxj,8LL)*

            (-49140LL - 98865LL*rxj + 3510LL*Power(rxj,2LL) - 131040LL*Power(rxj,3LL) - 

              7800LL*Power(rxj,4LL) + 3204LL*Power(rxj,5LL) + 360LL*Power(rxj,6LL) + 

              8LL*Power(rxj,7LL)) + Power(rxi,6LL)*Power(rxj,14LL)*

            (446489820LL - 54796455LL*rxj - 68983110LL*Power(rxj,2LL) - 

              12782700LL*Power(rxj,3LL) - 663600LL*Power(rxj,4LL) + 53928LL*Power(rxj,5LL) + 

              7728LL*Power(rxj,6LL) + 256LL*Power(rxj,7LL))))/

      (1260LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj,15LL)*Power(rxi + rxj,15LL))

     ; }
   
  }
  return S;
}

cl_R Slater_4S_5S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = (234137LL*xi)/1.31072e6

    ; } else {  S = (1LL/r)*(-14227497123840000LL + 14227497123840000LL*exp(2LL*rxi) - 

        25913502934444125LL*rxi - 23372011621208250LL*Power(rxi,2LL) - 

        13907709869303250LL*Power(rxi,3LL) - 6137735659555500LL*Power(rxi,4LL) - 

        2140857388870200LL*Power(rxi,5LL) - 614116575072000LL*Power(rxi,6LL) - 

        148809580920000LL*Power(rxi,7LL) - 31036639233600LL*Power(rxi,8LL) - 

        5645342102400LL*Power(rxi,9LL) - 903333150720LL*Power(rxi,10LL) - 

        127744081920LL*Power(rxi,11LL) - 15968010240LL*Power(rxi,12LL) - 

        1754726400LL*Power(rxi,13LL) - 167116800LL*Power(rxi,14LL) - 

        13369344LL*Power(rxi,15LL) - 835584LL*Power(rxi,16LL) - 32768LL*Power(rxi,17LL))/

      (1.422749712384e16*exp(2LL*rxi))

    ; }
 
  }
  else {
      if (r == 0LL) { S = (xi*xj*(4LL*Power(xi,16LL) + 68LL*Power(xi,15LL)*xj + 544LL*Power(xi,14LL)*Power(xj,2LL) + 

          2720LL*Power(xi,13LL)*Power(xj,3LL) + 9520LL*Power(xi,12LL)*Power(xj,4LL) + 

          24752LL*Power(xi,11LL)*Power(xj,5LL) + 49504LL*Power(xi,10LL)*Power(xj,6LL) + 

          77792LL*Power(xi,9LL)*Power(xj,7LL) + 97240LL*Power(xi,8LL)*Power(xj,8LL) + 

          97240LL*Power(xi,7LL)*Power(xj,9LL) + 61880LL*Power(xi,6LL)*Power(xj,10LL) + 

          30940LL*Power(xi,5LL)*Power(xj,11LL) + 11900LL*Power(xi,4LL)*Power(xj,12LL) + 

          3400LL*Power(xi,3LL)*Power(xj,13LL) + 680LL*Power(xi,2LL)*Power(xj,14LL) + 

          85LL*xi*Power(xj,15LL) + 5LL*Power(xj,16LL)))/(20LL*Power(xi + xj,17LL))

    ; } else { S = (1LL/r)*(56700LL*exp(2LL*(rxi + rxj))*Power(Power(rxi,2LL) - Power(rxj,2LL),17LL) + 

        9LL*exp(2LL*rxj)*Power(rxj,12LL)*

         (-980LL*Power(rxi,28LL) - 20LL*Power(rxi,29LL) + 6300LL*Power(rxj,22LL) + 

           11025LL*rxi*Power(rxj,22LL) - 50LL*Power(rxi,27LL)*(441LL + 2LL*Power(rxj,2LL)) + 

           3150LL*Power(rxi,2LL)*Power(rxj,20LL)*(-34LL + 3LL*Power(rxj,2LL)) + 

           525LL*Power(rxi,3LL)*Power(rxj,20LL)*(-357LL + 10LL*Power(rxj,2LL)) - 

           420LL*Power(rxi,26LL)*(700LL + 19LL*Power(rxj,2LL)) + 

           1050LL*Power(rxi,4LL)*Power(rxj,18LL)*

            (816LL - 153LL*Power(rxj,2LL) + 2LL*Power(rxj,4LL)) + 

           210LL*Power(rxi,5LL)*Power(rxj,18LL)*

            (7140LL - 425LL*Power(rxj,2LL) + 3LL*Power(rxj,4LL)) + 

           42LL*Power(rxi,25LL)*(-59500LL - 6035LL*Power(rxj,2LL) + 18LL*Power(rxj,4LL)) + 

           84LL*Power(rxi,24LL)*(-160650LL - 52700LL*Power(rxj,2LL) + 397LL*Power(rxj,4LL)) - 

           28LL*Power(rxi,12LL)*Power(rxj,10LL)*

            (100849950LL + 27100125LL*Power(rxj,2LL) + 186150LL*Power(rxj,4LL) - 

              2177LL*Power(rxj,6LL)) + 

           140LL*Power(rxi,6LL)*Power(rxj,16LL)*

            (-30600LL + 9180LL*Power(rxj,2LL) - 255LL*Power(rxj,4LL) + Power(rxj,6LL)) - 

           2380LL*Power(rxi,8LL)*Power(rxj,14LL)*

            (-6300LL + 2700LL*Power(rxj,2LL) - 120LL*Power(rxj,4LL) + Power(rxj,6LL)) + 

           10LL*Power(rxi,7LL)*Power(rxj,16LL)*

            (-749700LL + 71400LL*Power(rxj,2LL) - 1071LL*Power(rxj,4LL) + 2LL*Power(rxj,6LL)) 

    + 204LL*Power(rxi,15LL)*Power(rxj,8LL)*

            (28962255LL - 1744750LL*Power(rxj,2LL) + 9555LL*Power(rxj,4LL) + 

              6LL*Power(rxj,6LL)) - 42LL*Power(rxi,11LL)*Power(rxj,12LL)*

            (-12911925LL - 1634550LL*Power(rxj,2LL) - 7103LL*Power(rxj,4LL) + 

              18LL*Power(rxj,6LL)) + 2LL*Power(rxi,9LL)*Power(rxj,14LL)*

            (16948575LL - 1184400LL*Power(rxj,2LL) + 63861LL*Power(rxj,4LL) + 

              50LL*Power(rxj,6LL)) + 28LL*Power(rxi,22LL)*

            (-2180250LL - 10993050LL*Power(rxj,2LL) + 14925LL*Power(rxj,4LL) + 

              73LL*Power(rxj,6LL)) - 952LL*Power(rxi,14LL)*Power(rxj,8LL)*

            (16966215LL + 725175LL*Power(rxj,2LL) - 36075LL*Power(rxj,4LL) + 

              79LL*Power(rxj,6LL)) - 84LL*Power(rxi,10LL)*Power(rxj,12LL)*

            (1723800LL + 279225LL*Power(rxj,2LL) + 45600LL*Power(rxj,4LL) + 

              107LL*Power(rxj,6LL)) - 35LL*Power(rxi,17LL)*Power(rxj,6LL)*

            (132637869LL - 2205240LL*Power(rxj,2LL) - 48348LL*Power(rxj,4LL) + 

              136LL*Power(rxj,6LL)) - 6LL*Power(rxi,21LL)*Power(rxj,2LL)*

            (192298050LL + 12644275LL*Power(rxj,2LL) - 218029LL*Power(rxj,4LL) + 

              204LL*Power(rxj,6LL)) + 4LL*Power(rxi,13LL)*Power(rxj,10LL)*

            (1259522775LL + 15895425LL*Power(rxj,2LL) - 493017LL*Power(rxj,4LL) + 

              263LL*Power(rxj,6LL)) - 140LL*Power(rxi,16LL)*Power(rxj,6LL)*

            (180826281LL - 15101406LL*Power(rxj,2LL) + 160140LL*Power(rxj,4LL) + 

              442LL*Power(rxj,6LL)) - 2LL*Power(rxi,23LL)*

            (21366450LL + 23526300LL*Power(rxj,2LL) - 246729LL*Power(rxj,4LL) + 

              526LL*Power(rxj,6LL)) + 7LL*Power(rxi,19LL)*Power(rxj,4LL)*

            (-811081215LL + 39095550LL*Power(rxj,2LL) - 515916LL*Power(rxj,4LL) + 

              680LL*Power(rxj,6LL)) + 70LL*Power(rxi,18LL)*Power(rxj,4LL)*

            (-180554454LL + 9873711LL*Power(rxj,2LL) - 414120LL*Power(rxj,4LL) + 

              2924LL*Power(rxj,6LL)) - 

           14LL*Power(rxi,20LL)*Power(rxj,2LL)*

            (136919700LL + 71867115LL*Power(rxj,2LL) - 2154150LL*Power(rxj,4LL) + 

              10268LL*Power(rxj,6LL))) - 

        4LL*exp(2LL*rxi)*Power(rxi,10LL)*

         (-10710LL*Power(rxi,12LL)*Power(rxj,12LL)*

            (-3555LL - 127008LL*rxj + 138384LL*Power(rxj,2LL) - 74556LL*Power(rxj,3LL) - 

              22284LL*Power(rxj,4LL) + 408LL*Power(rxj,5LL) + 576LL*Power(rxj,6LL) + 

              60LL*Power(rxj,7LL) + 2LL*Power(rxj,8LL)) + 

           2LL*Power(rxi,20LL)*Power(rxj,4LL)*

            (963900LL + 1735020LL*rxj + 1542240LL*Power(rxj,2LL) + 

              899640LL*Power(rxj,3LL) + 385560LL*Power(rxj,4LL) + 128520LL*Power(rxj,5LL) + 

              34272LL*Power(rxj,6LL) + 9126LL*Power(rxj,7LL) + 333LL*Power(rxj,8LL) - 

              20LL*Power(rxj,9LL)) - 2LL*Power(rxj,24LL)*

            (119041650LL + 107137485LL*rxj + 45110520LL*Power(rxj,2LL) + 

              11695320LL*Power(rxj,3LL) + 2063880LL*Power(rxj,4LL) + 

              257985LL*Power(rxj,5LL) + 22932LL*Power(rxj,6LL) + 1404LL*Power(rxj,7LL) + 

              54LL*Power(rxj,8LL) + Power(rxj,9LL)) + 

           2LL*Power(rxi,2LL)*Power(rxj,22LL)*

            (-3264488325LL - 2505368880LL*rxj - 881390160LL*Power(rxj,2LL) - 

              185775660LL*Power(rxj,3LL) - 25639740LL*Power(rxj,4LL) - 

              2361555LL*Power(rxj,5LL) - 139356LL*Power(rxj,6LL) - 4482LL*Power(rxj,7LL) - 

              27LL*Power(rxj,8LL) + 2LL*Power(rxj,9LL)) + 

           Power(rxi,24LL)*(14175LL + 25515LL*rxj + 22680LL*Power(rxj,2LL) + 

              13230LL*Power(rxj,3LL) + 5670LL*Power(rxj,4LL) + 1890LL*Power(rxj,5LL) + 

              504LL*Power(rxj,6LL) + 108LL*Power(rxj,7LL) + 18LL*Power(rxj,8LL) + 

              2LL*Power(rxj,9LL)) - 102LL*Power(rxi,10LL)*Power(rxj,14LL)*

            (44986725LL - 97433280LL*rxj + 44467920LL*Power(rxj,2LL) + 

              15857100LL*Power(rxj,3LL) - 457380LL*Power(rxj,4LL) - 

              620550LL*Power(rxj,5LL) - 83160LL*Power(rxj,6LL) - 4068LL*Power(rxj,7LL) - 

              6LL*Power(rxj,8LL) + 4LL*Power(rxj,9LL)) + 

           102LL*Power(rxi,14LL)*Power(rxj,10LL)*

            (-859950LL - 1437345LL*rxj - 2260440LL*Power(rxj,2LL) + 

              810810LL*Power(rxj,3LL) - 1056510LL*Power(rxj,4LL) - 

              217854LL*Power(rxj,5LL) + 6552LL*Power(rxj,6LL) + 3852LL*Power(rxj,7LL) + 

              258LL*Power(rxj,8LL) + 4LL*Power(rxj,9LL)) - 

           Power(rxi,22LL)*Power(rxj,2LL)*

            (240975LL + 433755LL*rxj + 385560LL*Power(rxj,2LL) + 224910LL*Power(rxj,3LL) + 

              96390LL*Power(rxj,4LL) + 32130LL*Power(rxj,5LL) + 8568LL*Power(rxj,6LL) + 

              1836LL*Power(rxj,7LL) + 306LL*Power(rxj,8LL) + 4LL*Power(rxj,9LL)) + 

           2LL*Power(rxi,4LL)*Power(rxj,20LL)*

            (-18032978565LL - 9823683240LL*rxj - 2047323600LL*Power(rxj,2LL) - 

              129098340LL*Power(rxj,3LL) + 26410860LL*Power(rxj,4LL) + 

              7094304LL*Power(rxj,5LL) + 788256LL*Power(rxj,6LL) + 48654LL*Power(rxj,7LL) + 

              1593LL*Power(rxj,8LL) + 20LL*Power(rxj,9LL)) - 

           6LL*Power(rxi,16LL)*Power(rxj,8LL)*

            (-5622750LL - 10120950LL*rxj - 8996400LL*Power(rxj,2LL) - 

              5698350LL*Power(rxj,3LL) - 897750LL*Power(rxj,4LL) - 

              1641591LL*Power(rxj,5LL) - 211932LL*Power(rxj,6LL) + 10224LL*Power(rxj,7LL) + 

              2364LL*Power(rxj,8LL) + 73LL*Power(rxj,9LL)) + 

           2LL*Power(rxi,18LL)*Power(rxj,6LL)*

            (-4819500LL - 8675100LL*rxj - 7711200LL*Power(rxj,2LL) - 

              4498200LL*Power(rxj,3LL) - 1927800LL*Power(rxj,4LL) - 

              561519LL*Power(rxj,5LL) - 279468LL*Power(rxj,6LL) - 20682LL*Power(rxj,7LL) + 

              1305LL*Power(rxj,8LL) + 106LL*Power(rxj,9LL)) + 

           3LL*Power(rxi,8LL)*Power(rxj,16LL)*

            (-9364244085LL + 6940428705LL*rxj + 2117684520LL*Power(rxj,2LL) - 

              230268150LL*Power(rxj,3LL) - 149610510LL*Power(rxj,4LL) - 

              21824334LL*Power(rxj,5LL) - 1223208LL*Power(rxj,6LL) + 

              12708LL*Power(rxj,7LL) + 4470LL*Power(rxj,8LL) + 146LL*Power(rxj,9LL)) - 

           Power(rxi,6LL)*Power(rxj,18LL)*

            (57304872765LL + 7147185255LL*rxj - 5801702760LL*Power(rxj,2LL) - 

              2053388610LL*Power(rxj,3LL) - 271655370LL*Power(rxj,4LL) - 

              10864854LL*Power(rxj,5LL) + 1337112LL*Power(rxj,6LL) + 

              202716LL*Power(rxj,7LL) + 10746LL*Power(rxj,8LL) + 212LL*Power(rxj,9LL))))/

      (56700LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj,17LL)*Power(rxi + rxj,17LL))

     ; }
   
  }
  return S;
}

cl_R Slater_4S_1S(cl_R r,cl_R xi,cl_R xj)
{
  return Slater_1S_4S(r,xj,xi);
}

cl_R Slater_4S_2S(cl_R r,cl_R xi,cl_R xj)
{
  return Slater_2S_4S(r,xj,xi);
}

cl_R Slater_4S_3S(cl_R r,cl_R xi,cl_R xj)
{
  return Slater_3S_4S(r,xj,xi);
}

cl_R Slater_5S_5S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = (43191LL*xi)/262144LL

    ; } else {  S = (1LL/r)*1LL + (-1LL - (481097LL*rxi)/262144LL - (218953LL*Power(rxi,2LL))/131072LL - 

         (988003LL*Power(rxi,3LL))/983040LL - (110459LL*Power(rxi,4LL))/245760LL - 

         (65243LL*Power(rxi,5LL))/409600LL - (4769LL*Power(rxi,6LL))/102400LL - 

         (186229LL*Power(rxi,7LL))/1.6128e7 - (713LL*Power(rxi,8LL))/288000LL - 

         (33791LL*Power(rxi,9LL))/7.2576e7 - (11LL*Power(rxi,10LL))/141750LL - 

         Power(rxi,11LL)/86625LL - (2LL*Power(rxi,12LL))/1.299375e6 - 

         (4LL*Power(rxi,13LL))/2.1718125e7 - Power(rxi,14LL)/5.0675625e7 - 

         (2LL*Power(rxi,15LL))/1.064188125e9 - Power(rxi,16LL)/6.38512875e9 - 

         Power(rxi,17LL)/9.0455990625e10 - Power(rxi,18LL)/1.62820783125e12 - 

         Power(rxi,19LL)/4.6403923190625e13)/exp(2LL*rxi)

    ; }
 
  }
  else {
      if (r == 0LL) { S = (xi*xj*(Power(xi,18LL) + 19LL*Power(xi,17LL)*xj + 171LL*Power(xi,16LL)*Power(xj,2LL) + 

          969LL*Power(xi,15LL)*Power(xj,3LL) + 3876LL*Power(xi,14LL)*Power(xj,4LL) + 

          11628LL*Power(xi,13LL)*Power(xj,5LL) + 27132LL*Power(xi,12LL)*Power(xj,6LL) + 

          50388LL*Power(xi,11LL)*Power(xj,7LL) + 75582LL*Power(xi,10LL)*Power(xj,8LL) + 

          92378LL*Power(xi,9LL)*Power(xj,9LL) + 75582LL*Power(xi,8LL)*Power(xj,10LL) + 

          50388LL*Power(xi,7LL)*Power(xj,11LL) + 27132LL*Power(xi,6LL)*Power(xj,12LL) + 

          11628LL*Power(xi,5LL)*Power(xj,13LL) + 3876LL*Power(xi,4LL)*Power(xj,14LL) + 

          969LL*Power(xi,3LL)*Power(xj,15LL) + 171LL*Power(xi,2LL)*Power(xj,16LL) + 

          19LL*xi*Power(xj,17LL) + Power(xj,18LL)))/(5LL*Power(xi + xj,19LL))

    ; } else { S = (1LL/r)*(70875LL*exp(2LL*(rxi + rxj))*Power(Power(rxi,2LL) - Power(rxj,2LL),19LL) + 

        exp(2LL*rxj)*Power(rxj,12LL)*

         (-630LL*Power(rxi,34LL) - 10LL*Power(rxi,35LL) + 70875LL*Power(rxj,26LL) + 

           127575LL*rxi*Power(rxj,26LL) - 30LL*Power(rxi,33LL)*(630LL + Power(rxj,2LL)) + 

           14175LL*Power(rxi,2LL)*Power(rxj,24LL)*(-95LL + 8LL*Power(rxj,2LL)) + 

           4725LL*Power(rxi,3LL)*Power(rxj,24LL)*(-513LL + 14LL*Power(rxj,2LL)) - 

           90LL*Power(rxi,32LL)*(3920LL + 43LL*Power(rxj,2LL)) + 

           4725LL*Power(rxi,5LL)*Power(rxj,22LL)*

            (4617LL - 266LL*Power(rxj,2LL) + 2LL*Power(rxj,4LL)) + 

           14175LL*Power(rxi,4LL)*Power(rxj,22LL)*

            (855LL - 152LL*Power(rxj,2LL) + 2LL*Power(rxj,4LL)) + 

           36LL*Power(rxi,31LL)*(-124950LL - 4985LL*Power(rxj,2LL) + 13LL*Power(rxj,4LL)) + 

           36LL*Power(rxi,30LL)*(-1124550LL - 127960LL*Power(rxj,2LL) + 

              863LL*Power(rxj,4LL)) + 135LL*Power(rxi,7LL)*Power(rxj,20LL)*

            (-915705LL + 83790LL*Power(rxj,2LL) - 1330LL*Power(rxj,4LL) + 4LL*Power(rxj,6LL)) 

    + 315LL*Power(rxi,6LL)*Power(rxj,20LL)*

            (-218025LL + 61560LL*Power(rxj,2LL) - 1710LL*Power(rxj,4LL) + 8LL*Power(rxj,6LL)) 

    - 36LL*Power(rxi,29LL)*(7122150LL + 2102730LL*Power(rxj,2LL) - 23294LL*Power(rxj,4LL) + 

              37LL*Power(rxj,6LL)) - 36LL*Power(rxi,28LL)*

            (30523500LL + 23401350LL*Power(rxj,2LL) - 299250LL*Power(rxj,4LL) + 

              1297LL*Power(rxj,6LL)) + 

           Power(rxi,17LL)*Power(rxj,10LL)*

            (1073961177975LL - 21753487980LL*Power(rxj,2LL) - 

              745994340LL*Power(rxj,4LL) + 5307156LL*Power(rxj,6LL) - 818LL*Power(rxj,8LL)) 

    + 10LL*Power(rxi,9LL)*Power(rxj,18LL)*

            (49448070LL - 6409935LL*Power(rxj,2LL) + 161595LL*Power(rxj,4LL) - 

              1026LL*Power(rxj,6LL) + Power(rxj,8LL)) + 

           90LL*Power(rxi,8LL)*Power(rxj,18LL)*

            (3052350LL - 1220940LL*Power(rxj,2LL) + 53865LL*Power(rxj,4LL) - 

              532LL*Power(rxj,6LL) + Power(rxj,8LL)) - 

           1710LL*Power(rxi,10LL)*Power(rxj,16LL)*

            (481950LL - 257040LL*Power(rxj,2LL) + 16065LL*Power(rxj,4LL) - 

              252LL*Power(rxj,6LL) + Power(rxj,8LL)) + 

           6LL*Power(rxi,11LL)*Power(rxj,16LL)*

            (-207559800LL + 50390550LL*Power(rxj,2LL) - 1165815LL*Power(rxj,4LL) + 

              21396LL*Power(rxj,6LL) + 5LL*Power(rxj,8LL)) - 

           18LL*Power(rxi,13LL)*Power(rxj,14LL)*

            (-1703720025LL - 155669850LL*Power(rxj,2LL) - 7410270LL*Power(rxj,4LL) - 

              1532LL*Power(rxj,6LL) + 26LL*Power(rxj,8LL)) + 

           18LL*Power(rxi,15LL)*Power(rxj,12LL)*

            (19380896325LL + 1329128850LL*Power(rxj,2LL) - 7608930LL*Power(rxj,4LL) - 

              116238LL*Power(rxj,6LL) + 74LL*Power(rxj,8LL)) - 

           18LL*Power(rxi,12LL)*Power(rxj,14LL)*

            (89026875LL + 179071200LL*Power(rxj,2LL) + 1552950LL*Power(rxj,4LL) + 

              295820LL*Power(rxj,6LL) + 146LL*Power(rxj,8LL)) + 

           18LL*Power(rxi,25LL)*Power(rxj,2LL)*

            (-5449970925LL - 1137574935LL*Power(rxj,2LL) + 37834755LL*Power(rxj,4LL) - 

              273062LL*Power(rxj,6LL) + 171LL*Power(rxj,8LL)) - 

           9LL*Power(rxi,19LL)*Power(rxj,8LL)*

            (-37914907275LL + 7613889570LL*Power(rxj,2LL) - 170524620LL*Power(rxj,4LL) + 

              397936LL*Power(rxj,6LL) + 342LL*Power(rxj,8LL)) + 

           Power(rxi,27LL)*(-2884470750LL - 6409935000LL*Power(rxj,2LL) + 

              28332990LL*Power(rxj,4LL) + 58104LL*Power(rxj,6LL) + 818LL*Power(rxj,8LL)) - 

           3LL*Power(rxi,23LL)*Power(rxj,4LL)*

            (219130630425LL - 11118046590LL*Power(rxj,2LL) + 

              327611970LL*Power(rxj,4LL) - 2920908LL*Power(rxj,6LL) + 2584LL*Power(rxj,8LL)

    ) + 3LL*Power(rxi,21LL)*Power(rxj,6LL)*

            (-345162539925LL + 19030764690LL*Power(rxj,2LL) - 

              141976170LL*Power(rxj,4LL) - 1441872LL*Power(rxj,6LL) + 2584LL*Power(rxj,8LL)

    ) + 63LL*Power(rxi,20LL)*Power(rxj,6LL)*

            (-50980542525LL + 6240202920LL*Power(rxj,2LL) - 201314310LL*Power(rxj,4LL) + 

              956080LL*Power(rxj,6LL) + 2584LL*Power(rxj,8LL)) + 

           18LL*Power(rxi,14LL)*Power(rxj,12LL)*

            (-7803332775LL - 2519206200LL*Power(rxj,2LL) - 119719950LL*Power(rxj,4LL) + 

              182280LL*Power(rxj,6LL) + 2734LL*Power(rxj,8LL)) - 

           18LL*Power(rxi,26LL)*(195859125LL + 1794781800LL*Power(rxj,2LL) + 

              67337235LL*Power(rxj,4LL) - 1659700LL*Power(rxj,6LL) + 4089LL*Power(rxj,8LL)) 

    + 9LL*Power(rxi,18LL)*Power(rxj,8LL)*

            (-357591274425LL + 8328390840LL*Power(rxj,2LL) + 

              912042180LL*Power(rxj,4LL) - 12842480LL*Power(rxj,6LL) + 

              10678LL*Power(rxj,8LL)) - 

           9LL*Power(rxi,16LL)*Power(rxj,10LL)*

            (128599724925LL + 21298077360LL*Power(rxj,2LL) - 

              267928500LL*Power(rxj,4LL) - 5458320LL*Power(rxj,6LL) + 

              14722LL*Power(rxj,8LL)) + 

           18LL*Power(rxi,24LL)*Power(rxj,2LL)*

            (-7604930025LL - 8866107180LL*Power(rxj,2LL) + 399272265LL*Power(rxj,4LL) - 

              5925780LL*Power(rxj,6LL) + 17651LL*Power(rxj,8LL)) - 

           9LL*Power(rxi,22LL)*Power(rxj,4LL)*

            (129194933175LL + 3909863160LL*Power(rxj,2LL) + 91420770LL*Power(rxj,4LL) - 

              8762040LL*Power(rxj,6LL) + 43928LL*Power(rxj,8LL))) + 

        exp(2LL*rxi)*Power(rxi,12LL)*

         (Power(rxi,8LL)*Power(rxj,18LL)*

            (3218321469825LL - 341234165475LL*rxj - 393132783960LL*Power(rxj,2LL) - 

              57092294070LL*Power(rxj,3LL) + 822786930LL*Power(rxj,4LL) + 

              982835910LL*Power(rxj,5LL) + 106664040LL*Power(rxj,6LL) + 

              4915116LL*Power(rxj,7LL) + 73602LL*Power(rxj,8LL) - 818LL*Power(rxj,9LL)) + 

           10LL*Power(rxj,26LL)*(352546425LL + 288447075LL*rxj + 

              109884600LL*Power(rxj,2LL) + 25639740LL*Power(rxj,3LL) + 

              4048380LL*Power(rxj,4LL) + 449820LL*Power(rxj,5LL) + 35280LL*Power(rxj,6LL) + 

              1890LL*Power(rxj,7LL) + 63LL*Power(rxj,8LL) + Power(rxj,9LL)) + 

           30LL*Power(rxi,2LL)*Power(rxj,24LL)*

            (4562958015LL + 3269982555LL*rxj + 1076869080LL*Power(rxj,2LL) + 

              213664500LL*Power(rxj,3LL) + 28081620LL*Power(rxj,4LL) + 

              2523276LL*Power(rxj,5LL) + 153552LL*Power(rxj,6LL) + 5982LL*Power(rxj,7LL) + 

              129LL*Power(rxj,8LL) + Power(rxj,9LL)) - 

           15LL*Power(rxi,24LL)*Power(rxj,2LL)*

            (-89775LL - 161595LL*rxj - 143640LL*Power(rxj,2LL) - 83790LL*Power(rxj,3LL) - 

              35910LL*Power(rxj,4LL) - 11970LL*Power(rxj,5LL) - 3192LL*Power(rxj,6LL) - 

              684LL*Power(rxj,7LL) - 114LL*Power(rxj,8LL) + 2LL*Power(rxj,9LL)) - 

           5LL*Power(rxi,26LL)*(14175LL + 25515LL*rxj + 22680LL*Power(rxj,2LL) + 

              13230LL*Power(rxj,3LL) + 5670LL*Power(rxj,4LL) + 1890LL*Power(rxj,5LL) + 

              504LL*Power(rxj,6LL) + 108LL*Power(rxj,7LL) + 18LL*Power(rxj,8LL) + 

              2LL*Power(rxj,9LL)) - 1938LL*Power(rxi,14LL)*Power(rxj,12LL)*

            (-826875LL + 15824025LL*rxj - 23398200LL*Power(rxj,2LL) + 

              12344850LL*Power(rxj,3LL) + 1244250LL*Power(rxj,4LL) - 

              384930LL*Power(rxj,5LL) - 59640LL*Power(rxj,6LL) - 1848LL*Power(rxj,7LL) + 

              84LL*Power(rxj,8LL) + 4LL*Power(rxj,9LL)) + 

           1938LL*Power(rxi,12LL)*Power(rxj,14LL)*

            (72476775LL - 180008325LL*rxj + 98907480LL*Power(rxj,2LL) + 

              11224710LL*Power(rxj,3LL) - 4235490LL*Power(rxj,4LL) - 

              791910LL*Power(rxj,5LL) - 31080LL*Power(rxj,6LL) + 2232LL*Power(rxj,7LL) + 

              204LL*Power(rxj,8LL) + 4LL*Power(rxj,9LL)) + 

           342LL*Power(rxi,16LL)*Power(rxj,10LL)*

            (2409750LL + 3641400LL*rxj + 9424800LL*Power(rxj,2LL) - 

              8193150LL*Power(rxj,3LL) + 6301050LL*Power(rxj,4LL) + 

              400470LL*Power(rxj,5LL) - 143640LL*Power(rxj,6LL) - 15518LL*Power(rxj,7LL) - 

              281LL*Power(rxj,8LL) + 9LL*Power(rxj,9LL)) - 

           171LL*Power(rxi,10LL)*Power(rxj,16LL)*

            (-6768406575LL + 6280474725LL*rxj + 438336360LL*Power(rxj,2LL) - 

              400731030LL*Power(rxj,3LL) - 74168430LL*Power(rxj,4LL) - 

              2490810LL*Power(rxj,5LL) + 461160LL*Power(rxj,6LL) + 51244LL*Power(rxj,7LL) + 

              1858LL*Power(rxj,8LL) + 18LL*Power(rxj,9LL)) + 

           9LL*Power(rxi,22LL)*Power(rxj,4LL)*

            (-1346625LL - 2423925LL*rxj - 2154600LL*Power(rxj,2LL) - 

              1256850LL*Power(rxj,3LL) - 538650LL*Power(rxj,4LL) - 

              179550LL*Power(rxj,5LL) - 47880LL*Power(rxj,6LL) - 14264LL*Power(rxj,7LL) + 

              292LL*Power(rxj,8LL) + 52LL*Power(rxj,9LL)) - 

           9LL*Power(rxi,4LL)*Power(rxj,22LL)*

            (-129194933175LL - 73043543475LL*rxj - 17732214360LL*Power(rxj,2LL) - 

              2275149870LL*Power(rxj,3LL) - 134674470LL*Power(rxj,4LL) + 

              3148110LL*Power(rxj,5LL) + 1197000LL*Power(rxj,6LL) + 

              93176LL*Power(rxj,7LL) + 3452LL*Power(rxj,8LL) + 52LL*Power(rxj,9LL)) + 

           9LL*Power(rxi,6LL)*Power(rxj,20LL)*

            (356863797675LL + 115054179975LL*rxj + 3909863160LL*Power(rxj,2LL) - 

              3706015530LL*Power(rxj,3LL) - 798544530LL*Power(rxj,4LL) - 

              75669510LL*Power(rxj,5LL) - 3319400LL*Power(rxj,6LL) - 

              6456LL*Power(rxj,7LL) + 5188LL*Power(rxj,8LL) + 148LL*Power(rxj,9LL)) - 

           9LL*Power(rxi,20LL)*Power(rxj,6LL)*

            (-7630875LL - 13735575LL*rxj - 12209400LL*Power(rxj,2LL) - 

              7122150LL*Power(rxj,3LL) - 3052350LL*Power(rxj,4LL) - 

              777210LL*Power(rxj,5LL) - 591640LL*Power(rxj,6LL) + 3064LL*Power(rxj,7LL) + 

              5468LL*Power(rxj,8LL) + 148LL*Power(rxj,9LL)) + 

           2LL*Power(rxi,18LL)*Power(rxj,8LL)*

            (-137355750LL - 247240350LL*rxj - 219769200LL*Power(rxj,2LL) - 

              151171650LL*Power(rxj,3LL) + 13976550LL*Power(rxj,4LL) - 

              66692430LL*Power(rxj,5LL) - 1640520LL*Power(rxj,6LL) + 

              1046142LL*Power(rxj,7LL) + 66249LL*Power(rxj,8LL) + 409LL*Power(rxj,9LL))))/

      (70875LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj,19LL)*Power(rxi + rxj,19LL))

     ; }
   
  }
  return S;
}

cl_R Slater_5S_1S(cl_R r,cl_R xi,cl_R xj)
{
  return Slater_1S_5S(r,xj,xi);
}

cl_R Slater_5S_2S(cl_R r,cl_R xi,cl_R xj)
{
  return Slater_2S_5S(r,xj,xi);
}

cl_R Slater_5S_3S(cl_R r,cl_R xi,cl_R xj)
{
  return Slater_3S_5S(r,xj,xi);
}

cl_R Slater_5S_4S(cl_R r,cl_R xi,cl_R xj)
{
  return Slater_4S_5S(r,xj,xi);
}

cl_R DSlater_1S_1S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = 0LL

    ; } else {  S = -(-33LL*xi + 48LL*exp(2LL*r*xi)*xi - 36LL*r*Power(xi,2LL) - 

          12LL*Power(r,2LL)*Power(xi,3LL))/(24LL*exp(2LL*r*xi)*r) + 

      (-24LL + 24LL*exp(2LL*r*xi) - 33LL*r*xi - 18LL*Power(r,2LL)*Power(xi,2LL) - 

         4LL*Power(r,3LL)*Power(xi,3LL))/(24LL*exp(2LL*r*xi)*Power(r,2LL)) + 

      (xi*(-24LL + 24LL*exp(2LL*r*xi) - 33LL*r*xi - 18LL*Power(r,2LL)*Power(xi,2LL) - 

           4LL*Power(r,3LL)*Power(xi,3LL)))/(12LL*exp(2LL*r*xi)*r)

    ; }
 
  }
  else {
      if (r == 0LL) {  S = 0LL

    ; } else {  S = (exp(2LL*r*(xi + xj))*Power(Power(xi,2LL) - Power(xj,2LL),3LL) + 

         exp(2LL*r*xj)*Power(xj,4LL)*

          (-3LL*Power(xi,2LL) - r*Power(xi,3LL) + Power(xj,2LL) + r*xi*Power(xj,2LL)) - 

         exp(2LL*r*xi)*Power(xi,4LL)*

          (Power(xi,2LL)*(1LL + r*xj) - Power(xj,2LL)*(3LL + r*xj)))/

       (exp(2LL*r*(xi + xj))*Power(r,2LL)*Power(xi - xj,3LL)*Power(xi + xj,3LL)) + 

      (2LL*(exp(2LL*r*(xi + xj))*Power(Power(xi,2LL) - Power(xj,2LL),3LL) + 

           exp(2LL*r*xj)*Power(xj,4LL)*

            (-3LL*Power(xi,2LL) - r*Power(xi,3LL) + Power(xj,2LL) + r*xi*Power(xj,2LL)) - 

           exp(2LL*r*xi)*Power(xi,4LL)*

            (Power(xi,2LL)*(1LL + r*xj) - Power(xj,2LL)*(3LL + r*xj))))/

       (exp(2LL*r*(xi + xj))*r*Power(xi - xj,3LL)*Power(xi + xj,2LL)) - 

      (2LL*exp(2LL*r*(xi + xj))*(xi + xj)*Power(Power(xi,2LL) - Power(xj,2LL),3LL) + 

         exp(2LL*r*xj)*Power(xj,4LL)*(-Power(xi,3LL) + xi*Power(xj,2LL)) + 

         2LL*exp(2LL*r*xj)*Power(xj,5LL)*

          (-3LL*Power(xi,2LL) - r*Power(xi,3LL) + Power(xj,2LL) + r*xi*Power(xj,2LL)) - 

         exp(2LL*r*xi)*Power(xi,4LL)*(Power(xi,2LL)*xj - Power(xj,3LL)) - 

         2LL*exp(2LL*r*xi)*Power(xi,5LL)*

          (Power(xi,2LL)*(1LL + r*xj) - Power(xj,2LL)*(3LL + r*xj)))/

       (exp(2LL*r*(xi + xj))*r*Power(xi - xj,3LL)*Power(xi + xj,3LL))

    ; }
   
  }
  return S;
}

cl_R DSlater_1S_2S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = 0LL

    ; } else {  S = -(-375LL*xi + 480LL*exp(2LL*r*xi)*xi - 540LL*r*Power(xi,2LL) - 

          345LL*Power(r,2LL)*Power(xi,3LL) - 120LL*Power(r,3LL)*Power(xi,4LL) - 

          20LL*Power(r,4LL)*Power(xi,5LL))/(240LL*exp(2LL*r*xi)*r) + 

      (-240LL + 240LL*exp(2LL*r*xi) - 375LL*r*xi - 270LL*Power(r,2LL)*Power(xi,2LL) - 

         115LL*Power(r,3LL)*Power(xi,3LL) - 30LL*Power(r,4LL)*Power(xi,4LL) - 

         4LL*Power(r,5LL)*Power(xi,5LL))/(240LL*exp(2LL*r*xi)*Power(r,2LL)) + 

      (xi*(-240LL + 240LL*exp(2LL*r*xi) - 375LL*r*xi - 270LL*Power(r,2LL)*Power(xi,2LL) - 

           115LL*Power(r,3LL)*Power(xi,3LL) - 30LL*Power(r,4LL)*Power(xi,4LL) - 

           4LL*Power(r,5LL)*Power(xi,5LL)))/(120LL*exp(2LL*r*xi)*r)

    ; }
 
  }
  else {
      if (r == 0LL) {  S = 0LL

    ; } else {  S = (6LL*exp(2LL*r*(xi + xj))*Power(Power(xi,2LL) - Power(xj,2LL),5LL) + 

         6LL*exp(2LL*r*xj)*Power(xj,6LL)*

          (-4LL*Power(xi,4LL) - r*Power(xi,5LL) - 5LL*Power(xi,2LL)*Power(xj,2LL) + 

            Power(xj,4LL) + r*xi*Power(xj,4LL)) - 

         exp(2LL*r*xi)*Power(xi,4LL)*

          (Power(xi,6LL)*(6LL + 9LL*r*xj + 6LL*Power(r,2LL)*Power(xj,2LL) + 

               2LL*Power(r,3LL)*Power(xj,3LL)) - 

            3LL*Power(xi,4LL)*Power(xj,2LL)*

             (10LL + 15LL*r*xj + 10LL*Power(r,2LL)*Power(xj,2LL) + 

               2LL*Power(r,3LL)*Power(xj,3LL)) + 

            3LL*Power(xi,2LL)*Power(xj,4LL)*

             (20LL + 33LL*r*xj + 14LL*Power(r,2LL)*Power(xj,2LL) + 

               2LL*Power(r,3LL)*Power(xj,3LL)) - 

            Power(xj,6LL)*(84LL + 63LL*r*xj + 18LL*Power(r,2LL)*Power(xj,2LL) + 

               2LL*Power(r,3LL)*Power(xj,3LL))))/

       (6LL*exp(2LL*r*(xi + xj))*Power(r,2LL)*Power(xi - xj,5LL)*Power(xi + xj,5LL)) + 

      (6LL*exp(2LL*r*(xi + xj))*Power(Power(xi,2LL) - Power(xj,2LL),5LL) + 

         6LL*exp(2LL*r*xj)*Power(xj,6LL)*

          (-4LL*Power(xi,4LL) - r*Power(xi,5LL) - 5LL*Power(xi,2LL)*Power(xj,2LL) + 

            Power(xj,4LL) + r*xi*Power(xj,4LL)) - 

         exp(2LL*r*xi)*Power(xi,4LL)*

          (Power(xi,6LL)*(6LL + 9LL*r*xj + 6LL*Power(r,2LL)*Power(xj,2LL) + 

               2LL*Power(r,3LL)*Power(xj,3LL)) - 

            3LL*Power(xi,4LL)*Power(xj,2LL)*

             (10LL + 15LL*r*xj + 10LL*Power(r,2LL)*Power(xj,2LL) + 

               2LL*Power(r,3LL)*Power(xj,3LL)) + 

            3LL*Power(xi,2LL)*Power(xj,4LL)*

             (20LL + 33LL*r*xj + 14LL*Power(r,2LL)*Power(xj,2LL) + 

               2LL*Power(r,3LL)*Power(xj,3LL)) - 

            Power(xj,6LL)*(84LL + 63LL*r*xj + 18LL*Power(r,2LL)*Power(xj,2LL) + 

               2LL*Power(r,3LL)*Power(xj,3LL))))/

       (3LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,5LL)*Power(xi + xj,4LL)) - 

      (12LL*exp(2LL*r*(xi + xj))*(xi + xj)*Power(Power(xi,2LL) - Power(xj,2LL),5LL) + 

         6LL*exp(2LL*r*xj)*Power(xj,6LL)*(-Power(xi,5LL) + xi*Power(xj,4LL)) + 

         12LL*exp(2LL*r*xj)*Power(xj,7LL)*

          (-4LL*Power(xi,4LL) - r*Power(xi,5LL) - 5LL*Power(xi,2LL)*Power(xj,2LL) + 

            Power(xj,4LL) + r*xi*Power(xj,4LL)) - 

         exp(2LL*r*xi)*Power(xi,4LL)*

          (Power(xi,6LL)*(9LL*xj + 12LL*r*Power(xj,2LL) + 6LL*Power(r,2LL)*Power(xj,3LL)) - 

            3LL*Power(xi,4LL)*Power(xj,2LL)*

             (15LL*xj + 20LL*r*Power(xj,2LL) + 6LL*Power(r,2LL)*Power(xj,3LL)) + 

            3LL*Power(xi,2LL)*Power(xj,4LL)*

             (33LL*xj + 28LL*r*Power(xj,2LL) + 6LL*Power(r,2LL)*Power(xj,3LL)) - 

            Power(xj,6LL)*(63LL*xj + 36LL*r*Power(xj,2LL) + 6LL*Power(r,2LL)*Power(xj,3LL))) - 

         2LL*exp(2LL*r*xi)*Power(xi,5LL)*

          (Power(xi,6LL)*(6LL + 9LL*r*xj + 6LL*Power(r,2LL)*Power(xj,2LL) + 

               2LL*Power(r,3LL)*Power(xj,3LL)) - 

            3LL*Power(xi,4LL)*Power(xj,2LL)*

             (10LL + 15LL*r*xj + 10LL*Power(r,2LL)*Power(xj,2LL) + 

               2LL*Power(r,3LL)*Power(xj,3LL)) + 

            3LL*Power(xi,2LL)*Power(xj,4LL)*

             (20LL + 33LL*r*xj + 14LL*Power(r,2LL)*Power(xj,2LL) + 

               2LL*Power(r,3LL)*Power(xj,3LL)) - 

            Power(xj,6LL)*(84LL + 63LL*r*xj + 18LL*Power(r,2LL)*Power(xj,2LL) + 

               2LL*Power(r,3LL)*Power(xj,3LL))))/

       (6LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,5LL)*Power(xi + xj,5LL))

    ; }
   
  }
  return S;
}

cl_R DSlater_1S_3S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = 0LL

    ; } else {  S = -(-203175LL*xi + 241920LL*exp(2LL*r*xi)*xi - 328860LL*r*Power(xi,2LL) - 

          253260LL*Power(r,2LL)*Power(xi,3LL) - 120960LL*Power(r,3LL)*Power(xi,4LL) - 

          38640LL*Power(r,4LL)*Power(xi,5LL) - 8064LL*Power(r,5LL)*Power(xi,6LL) - 

          896LL*Power(r,6LL)*Power(xi,7LL))/(120960LL*exp(2LL*r*xi)*r) + 

      (-120960LL + 120960LL*exp(2LL*r*xi) - 203175LL*r*xi - 

         164430LL*Power(r,2LL)*Power(xi,2LL) - 84420LL*Power(r,3LL)*Power(xi,3LL) - 

         30240LL*Power(r,4LL)*Power(xi,4LL) - 7728LL*Power(r,5LL)*Power(xi,5LL) - 

         1344LL*Power(r,6LL)*Power(xi,6LL) - 128LL*Power(r,7LL)*Power(xi,7LL))/

       (120960LL*exp(2LL*r*xi)*Power(r,2LL)) + 

      (xi*(-120960LL + 120960LL*exp(2LL*r*xi) - 203175LL*r*xi - 

           164430LL*Power(r,2LL)*Power(xi,2LL) - 84420LL*Power(r,3LL)*Power(xi,3LL) - 

           30240LL*Power(r,4LL)*Power(xi,4LL) - 7728LL*Power(r,5LL)*Power(xi,5LL) - 

           1344LL*Power(r,6LL)*Power(xi,6LL) - 128LL*Power(r,7LL)*Power(xi,7LL)))/

       (60480LL*exp(2LL*r*xi)*r)

    ; }
 
  }
  else {
      if (r == 0LL) {  S = 0LL

    ; } else {  S = (45LL*exp(2LL*r*(xi + xj))*Power(Power(xi,2LL) - Power(xj,2LL),7LL) + 

         15LL*exp(2LL*r*xj)*Power(xj,8LL)*

          (-15LL*Power(xi,6LL) - 3LL*r*Power(xi,7LL) - 63LL*Power(xi,4LL)*Power(xj,2LL) - 

            7LL*r*Power(xi,5LL)*Power(xj,2LL) - 21LL*Power(xi,2LL)*Power(xj,4LL) + 

            7LL*r*Power(xi,3LL)*Power(xj,4LL) + 3LL*Power(xj,6LL) + 3LL*r*xi*Power(xj,6LL)) + 

         exp(2LL*r*xi)*Power(xi,4LL)*

          (-10LL*Power(xi,2LL)*Power(xj,8LL)*

             (135LL + 333LL*r*xj + 228LL*Power(r,2LL)*Power(xj,2LL) + 

               75LL*Power(r,3LL)*Power(xj,3LL) + 13LL*Power(r,4LL)*Power(xj,4LL) + 

               Power(r,5LL)*Power(xj,5LL)) + 

            2LL*Power(xj,10LL)*(945LL + 945LL*r*xj + 420LL*Power(r,2LL)*Power(xj,2LL) + 

               105LL*Power(r,3LL)*Power(xj,3LL) + 15LL*Power(r,4LL)*Power(xj,4LL) + 

               Power(r,5LL)*Power(xj,5LL)) - 

            Power(xi,10LL)*(45LL + 75LL*r*xj + 60LL*Power(r,2LL)*Power(xj,2LL) + 

               30LL*Power(r,3LL)*Power(xj,3LL) + 10LL*Power(r,4LL)*Power(xj,4LL) + 

               2LL*Power(r,5LL)*Power(xj,5LL)) + 

            5LL*Power(xi,8LL)*Power(xj,2LL)*

             (63LL + 105LL*r*xj + 84LL*Power(r,2LL)*Power(xj,2LL) + 

               42LL*Power(r,3LL)*Power(xj,3LL) + 14LL*Power(r,4LL)*Power(xj,4LL) + 

               2LL*Power(r,5LL)*Power(xj,5LL)) - 

            5LL*Power(xi,6LL)*Power(xj,4LL)*

             (189LL + 315LL*r*xj + 252LL*Power(r,2LL)*Power(xj,2LL) + 

               132LL*Power(r,3LL)*Power(xj,3LL) + 36LL*Power(r,4LL)*Power(xj,4LL) + 

               4LL*Power(r,5LL)*Power(xj,5LL)) + 

            5LL*Power(xi,4LL)*Power(xj,6LL)*

             (315LL + 513LL*r*xj + 468LL*Power(r,2LL)*Power(xj,2LL) + 

               204LL*Power(r,3LL)*Power(xj,3LL) + 44LL*Power(r,4LL)*Power(xj,4LL) + 

               4LL*Power(r,5LL)*Power(xj,5LL))))/

       (45LL*exp(2LL*r*(xi + xj))*Power(r,2LL)*Power(xi - xj,7LL)*Power(xi + xj,7LL)) 

    + (2LL*(45LL*exp(2LL*r*(xi + xj))*Power(Power(xi,2LL) - Power(xj,2LL),7LL) + 

           15LL*exp(2LL*r*xj)*Power(xj,8LL)*

            (-15LL*Power(xi,6LL) - 3LL*r*Power(xi,7LL) - 63LL*Power(xi,4LL)*Power(xj,2LL) - 

              7LL*r*Power(xi,5LL)*Power(xj,2LL) - 21LL*Power(xi,2LL)*Power(xj,4LL) + 

              7LL*r*Power(xi,3LL)*Power(xj,4LL) + 3LL*Power(xj,6LL) + 3LL*r*xi*Power(xj,6LL)) 

    + exp(2LL*r*xi)*Power(xi,4LL)*(-10LL*Power(xi,2LL)*Power(xj,8LL)*

               (135LL + 333LL*r*xj + 228LL*Power(r,2LL)*Power(xj,2LL) + 

                 75LL*Power(r,3LL)*Power(xj,3LL) + 13LL*Power(r,4LL)*Power(xj,4LL) + 

                 Power(r,5LL)*Power(xj,5LL)) + 

              2LL*Power(xj,10LL)*(945LL + 945LL*r*xj + 420LL*Power(r,2LL)*Power(xj,2LL) + 

                 105LL*Power(r,3LL)*Power(xj,3LL) + 15LL*Power(r,4LL)*Power(xj,4LL) + 

                 Power(r,5LL)*Power(xj,5LL)) - 

              Power(xi,10LL)*(45LL + 75LL*r*xj + 60LL*Power(r,2LL)*Power(xj,2LL) + 

                 30LL*Power(r,3LL)*Power(xj,3LL) + 10LL*Power(r,4LL)*Power(xj,4LL) + 

                 2LL*Power(r,5LL)*Power(xj,5LL)) + 

              5LL*Power(xi,8LL)*Power(xj,2LL)*

               (63LL + 105LL*r*xj + 84LL*Power(r,2LL)*Power(xj,2LL) + 

                 42LL*Power(r,3LL)*Power(xj,3LL) + 14LL*Power(r,4LL)*Power(xj,4LL) + 

                 2LL*Power(r,5LL)*Power(xj,5LL)) - 

              5LL*Power(xi,6LL)*Power(xj,4LL)*

               (189LL + 315LL*r*xj + 252LL*Power(r,2LL)*Power(xj,2LL) + 

                 132LL*Power(r,3LL)*Power(xj,3LL) + 36LL*Power(r,4LL)*Power(xj,4LL) + 

                 4LL*Power(r,5LL)*Power(xj,5LL)) + 

              5LL*Power(xi,4LL)*Power(xj,6LL)*

               (315LL + 513LL*r*xj + 468LL*Power(r,2LL)*Power(xj,2LL) + 

                 204LL*Power(r,3LL)*Power(xj,3LL) + 44LL*Power(r,4LL)*Power(xj,4LL) + 

                 4LL*Power(r,5LL)*Power(xj,5LL)))))/

       (45LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,7LL)*Power(xi + xj,6LL)) - 

      (90LL*exp(2LL*r*(xi + xj))*(xi + xj)*Power(Power(xi,2LL) - Power(xj,2LL),7LL) + 

         15LL*exp(2LL*r*xj)*Power(xj,8LL)*

          (-3LL*Power(xi,7LL) - 7LL*Power(xi,5LL)*Power(xj,2LL) + 

            7LL*Power(xi,3LL)*Power(xj,4LL) + 3LL*xi*Power(xj,6LL)) + 

         30LL*exp(2LL*r*xj)*Power(xj,9LL)*

          (-15LL*Power(xi,6LL) - 3LL*r*Power(xi,7LL) - 63LL*Power(xi,4LL)*Power(xj,2LL) - 

            7LL*r*Power(xi,5LL)*Power(xj,2LL) - 21LL*Power(xi,2LL)*Power(xj,4LL) + 

            7LL*r*Power(xi,3LL)*Power(xj,4LL) + 3LL*Power(xj,6LL) + 3LL*r*xi*Power(xj,6LL)) + 

         exp(2LL*r*xi)*Power(xi,4LL)*

          (-10LL*Power(xi,2LL)*Power(xj,8LL)*

             (333LL*xj + 456LL*r*Power(xj,2LL) + 225LL*Power(r,2LL)*Power(xj,3LL) + 

               52LL*Power(r,3LL)*Power(xj,4LL) + 5LL*Power(r,4LL)*Power(xj,5LL)) + 

            2LL*Power(xj,10LL)*(945LL*xj + 840LL*r*Power(xj,2LL) + 

               315LL*Power(r,2LL)*Power(xj,3LL) + 60LL*Power(r,3LL)*Power(xj,4LL) + 

               5LL*Power(r,4LL)*Power(xj,5LL)) - 

            Power(xi,10LL)*(75LL*xj + 120LL*r*Power(xj,2LL) + 

               90LL*Power(r,2LL)*Power(xj,3LL) + 40LL*Power(r,3LL)*Power(xj,4LL) + 

               10LL*Power(r,4LL)*Power(xj,5LL)) + 

            5LL*Power(xi,8LL)*Power(xj,2LL)*

             (105LL*xj + 168LL*r*Power(xj,2LL) + 126LL*Power(r,2LL)*Power(xj,3LL) + 

               56LL*Power(r,3LL)*Power(xj,4LL) + 10LL*Power(r,4LL)*Power(xj,5LL)) - 

            5LL*Power(xi,6LL)*Power(xj,4LL)*

             (315LL*xj + 504LL*r*Power(xj,2LL) + 396LL*Power(r,2LL)*Power(xj,3LL) + 

               144LL*Power(r,3LL)*Power(xj,4LL) + 20LL*Power(r,4LL)*Power(xj,5LL)) + 

            5LL*Power(xi,4LL)*Power(xj,6LL)*

             (513LL*xj + 936LL*r*Power(xj,2LL) + 612LL*Power(r,2LL)*Power(xj,3LL) + 

               176LL*Power(r,3LL)*Power(xj,4LL) + 20LL*Power(r,4LL)*Power(xj,5LL))) + 

         2LL*exp(2LL*r*xi)*Power(xi,5LL)*

          (-10LL*Power(xi,2LL)*Power(xj,8LL)*

             (135LL + 333LL*r*xj + 228LL*Power(r,2LL)*Power(xj,2LL) + 

               75LL*Power(r,3LL)*Power(xj,3LL) + 13LL*Power(r,4LL)*Power(xj,4LL) + 

               Power(r,5LL)*Power(xj,5LL)) + 

            2LL*Power(xj,10LL)*(945LL + 945LL*r*xj + 420LL*Power(r,2LL)*Power(xj,2LL) + 

               105LL*Power(r,3LL)*Power(xj,3LL) + 15LL*Power(r,4LL)*Power(xj,4LL) + 

               Power(r,5LL)*Power(xj,5LL)) - 

            Power(xi,10LL)*(45LL + 75LL*r*xj + 60LL*Power(r,2LL)*Power(xj,2LL) + 

               30LL*Power(r,3LL)*Power(xj,3LL) + 10LL*Power(r,4LL)*Power(xj,4LL) + 

               2LL*Power(r,5LL)*Power(xj,5LL)) + 

            5LL*Power(xi,8LL)*Power(xj,2LL)*

             (63LL + 105LL*r*xj + 84LL*Power(r,2LL)*Power(xj,2LL) + 

               42LL*Power(r,3LL)*Power(xj,3LL) + 14LL*Power(r,4LL)*Power(xj,4LL) + 

               2LL*Power(r,5LL)*Power(xj,5LL)) - 

            5LL*Power(xi,6LL)*Power(xj,4LL)*

             (189LL + 315LL*r*xj + 252LL*Power(r,2LL)*Power(xj,2LL) + 

               132LL*Power(r,3LL)*Power(xj,3LL) + 36LL*Power(r,4LL)*Power(xj,4LL) + 

               4LL*Power(r,5LL)*Power(xj,5LL)) + 

            5LL*Power(xi,4LL)*Power(xj,6LL)*

             (315LL + 513LL*r*xj + 468LL*Power(r,2LL)*Power(xj,2LL) + 

               204LL*Power(r,3LL)*Power(xj,3LL) + 44LL*Power(r,4LL)*Power(xj,4LL) + 

               4LL*Power(r,5LL)*Power(xj,5LL))))/

       (45LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,7LL)*Power(xi + xj,7LL))

    ; }
   
  }
  return S;
}

cl_R DSlater_1S_4S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = 0LL

    ; } else {  S = -(-5088825LL*xi + 5806080LL*exp(2LL*r*xi)*xi - 8743140LL*r*Power(xi,2LL) - 

          7319970LL*Power(r,2LL)*Power(xi,3LL) - 3946320LL*Power(r,3LL)*Power(xi,4LL) - 

          1519560LL*Power(r,4LL)*Power(xi,5LL) - 435456LL*Power(r,5LL)*Power(xi,6LL) - 

          92736LL*Power(r,6LL)*Power(xi,7LL) - 13824LL*Power(r,7LL)*Power(xi,8LL) - 

          1152LL*Power(r,8LL)*Power(xi,9LL))/(2.90304e6*exp(2LL*r*xi)*r) + 

      (-2903040LL + 2903040LL*exp(2LL*r*xi) - 5088825LL*r*xi - 

         4371570LL*Power(r,2LL)*Power(xi,2LL) - 2439990LL*Power(r,3LL)*Power(xi,3LL) - 

         986580LL*Power(r,4LL)*Power(xi,4LL) - 303912LL*Power(r,5LL)*Power(xi,5LL) - 

         72576LL*Power(r,6LL)*Power(xi,6LL) - 13248LL*Power(r,7LL)*Power(xi,7LL) - 

         1728LL*Power(r,8LL)*Power(xi,8LL) - 128LL*Power(r,9LL)*Power(xi,9LL))/

       (2.90304e6*exp(2LL*r*xi)*Power(r,2LL)) + 

      (xi*(-2903040LL + 2903040LL*exp(2LL*r*xi) - 5088825LL*r*xi - 

           4371570LL*Power(r,2LL)*Power(xi,2LL) - 2439990LL*Power(r,3LL)*Power(xi,3LL) - 

           986580LL*Power(r,4LL)*Power(xi,4LL) - 303912LL*Power(r,5LL)*Power(xi,5LL) - 

           72576LL*Power(r,6LL)*Power(xi,6LL) - 13248LL*Power(r,7LL)*Power(xi,7LL) - 

           1728LL*Power(r,8LL)*Power(xi,8LL) - 128LL*Power(r,9LL)*Power(xi,9LL)))/

       (1.45152e6*exp(2LL*r*xi)*r)

    ; }
 
  }
  else {
      if (r == 0LL) {  S = 0LL

    ; } else {  S = (1260LL*exp(2LL*r*(xi + xj))*Power(Power(xi,2LL) - Power(xj,2LL),9LL) + 

         1260LL*exp(2LL*r*xj)*Power(xj,10LL)*

          (-6LL*Power(xi,8LL) - r*Power(xi,9LL) - 51LL*Power(xi,6LL)*Power(xj,2LL) - 

            6LL*r*Power(xi,7LL)*Power(xj,2LL) - 63LL*Power(xi,4LL)*Power(xj,4LL) - 

            9LL*Power(xi,2LL)*Power(xj,6LL) + 6LL*r*Power(xi,3LL)*Power(xj,6LL) + 

            Power(xj,8LL) + r*xi*Power(xj,8LL)) + 

         exp(2LL*r*xi)*Power(xi,4LL)*

          (-42LL*Power(xi,10LL)*Power(xj,4LL)*

             (1080LL + 1890LL*r*xj + 1620LL*Power(r,2LL)*Power(xj,2LL) + 

               900LL*Power(r,3LL)*Power(xj,3LL) + 360LL*Power(r,4LL)*Power(xj,4LL) + 

               111LL*Power(r,5LL)*Power(xj,5LL) + 22LL*Power(r,6LL)*Power(xj,6LL) + 

               2LL*Power(r,7LL)*Power(xj,7LL)) + 

            70LL*Power(xi,8LL)*Power(xj,6LL)*

             (1512LL + 2646LL*r*xj + 2268LL*Power(r,2LL)*Power(xj,2LL) + 

               1248LL*Power(r,3LL)*Power(xj,3LL) + 528LL*Power(r,4LL)*Power(xj,4LL) + 

               153LL*Power(r,5LL)*Power(xj,5LL) + 26LL*Power(r,6LL)*Power(xj,6LL) + 

               2LL*Power(r,7LL)*Power(xj,7LL)) - 

            14LL*Power(xi,2LL)*Power(xj,12LL)*

             (2970LL + 16335LL*r*xj + 15390LL*Power(r,2LL)*Power(xj,2LL) + 

               7110LL*Power(r,3LL)*Power(xj,3LL) + 1980LL*Power(r,4LL)*Power(xj,4LL) + 

               351LL*Power(r,5LL)*Power(xj,5LL) + 38LL*Power(r,6LL)*Power(xj,6LL) + 

               2LL*Power(r,7LL)*Power(xj,7LL)) + 

            2LL*Power(xj,14LL)*(62370LL + 72765LL*r*xj + 

               39690LL*Power(r,2LL)*Power(xj,2LL) + 13230LL*Power(r,3LL)*Power(xj,3LL) + 

               2940LL*Power(r,4LL)*Power(xj,4LL) + 441LL*Power(r,5LL)*Power(xj,5LL) + 

               42LL*Power(r,6LL)*Power(xj,6LL) + 2LL*Power(r,7LL)*Power(xj,7LL)) - 

            Power(xi,14LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r,2LL)*Power(xj,2LL) + 

               1050LL*Power(r,3LL)*Power(xj,3LL) + 420LL*Power(r,4LL)*Power(xj,4LL) + 

               126LL*Power(r,5LL)*Power(xj,5LL) + 28LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL)) + 

            7LL*Power(xi,12LL)*Power(xj,2LL)*

             (1620LL + 2835LL*r*xj + 2430LL*Power(r,2LL)*Power(xj,2LL) + 

               1350LL*Power(r,3LL)*Power(xj,3LL) + 540LL*Power(r,4LL)*Power(xj,4LL) + 

               162LL*Power(r,5LL)*Power(xj,5LL) + 36LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL)) - 

            35LL*Power(xi,6LL)*Power(xj,8LL)*

             (4536LL + 7983LL*r*xj + 6534LL*Power(r,2LL)*Power(xj,2LL) + 

               4014LL*Power(r,3LL)*Power(xj,3LL) + 1644LL*Power(r,4LL)*Power(xj,4LL) + 

               414LL*Power(r,5LL)*Power(xj,5LL) + 60LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL)) + 

            21LL*Power(xi,4LL)*Power(xj,10LL)*

             (7920LL + 11385LL*r*xj + 12330LL*Power(r,2LL)*Power(xj,2LL) + 

               7410LL*Power(r,3LL)*Power(xj,3LL) + 2580LL*Power(r,4LL)*Power(xj,4LL) + 

               546LL*Power(r,5LL)*Power(xj,5LL) + 68LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL))))/

       (1260LL*exp(2LL*r*(xi + xj))*Power(r,2LL)*Power(xi - xj,9LL)*

         Power(xi + xj,9LL)) + (1260LL*exp(2LL*r*(xi + xj))*

          Power(Power(xi,2LL) - Power(xj,2LL),9LL) + 

         1260LL*exp(2LL*r*xj)*Power(xj,10LL)*

          (-6LL*Power(xi,8LL) - r*Power(xi,9LL) - 51LL*Power(xi,6LL)*Power(xj,2LL) - 

            6LL*r*Power(xi,7LL)*Power(xj,2LL) - 63LL*Power(xi,4LL)*Power(xj,4LL) - 

            9LL*Power(xi,2LL)*Power(xj,6LL) + 6LL*r*Power(xi,3LL)*Power(xj,6LL) + 

            Power(xj,8LL) + r*xi*Power(xj,8LL)) + 

         exp(2LL*r*xi)*Power(xi,4LL)*

          (-42LL*Power(xi,10LL)*Power(xj,4LL)*

             (1080LL + 1890LL*r*xj + 1620LL*Power(r,2LL)*Power(xj,2LL) + 

               900LL*Power(r,3LL)*Power(xj,3LL) + 360LL*Power(r,4LL)*Power(xj,4LL) + 

               111LL*Power(r,5LL)*Power(xj,5LL) + 22LL*Power(r,6LL)*Power(xj,6LL) + 

               2LL*Power(r,7LL)*Power(xj,7LL)) + 

            70LL*Power(xi,8LL)*Power(xj,6LL)*

             (1512LL + 2646LL*r*xj + 2268LL*Power(r,2LL)*Power(xj,2LL) + 

               1248LL*Power(r,3LL)*Power(xj,3LL) + 528LL*Power(r,4LL)*Power(xj,4LL) + 

               153LL*Power(r,5LL)*Power(xj,5LL) + 26LL*Power(r,6LL)*Power(xj,6LL) + 

               2LL*Power(r,7LL)*Power(xj,7LL)) - 

            14LL*Power(xi,2LL)*Power(xj,12LL)*

             (2970LL + 16335LL*r*xj + 15390LL*Power(r,2LL)*Power(xj,2LL) + 

               7110LL*Power(r,3LL)*Power(xj,3LL) + 1980LL*Power(r,4LL)*Power(xj,4LL) + 

               351LL*Power(r,5LL)*Power(xj,5LL) + 38LL*Power(r,6LL)*Power(xj,6LL) + 

               2LL*Power(r,7LL)*Power(xj,7LL)) + 

            2LL*Power(xj,14LL)*(62370LL + 72765LL*r*xj + 

               39690LL*Power(r,2LL)*Power(xj,2LL) + 13230LL*Power(r,3LL)*Power(xj,3LL) + 

               2940LL*Power(r,4LL)*Power(xj,4LL) + 441LL*Power(r,5LL)*Power(xj,5LL) + 

               42LL*Power(r,6LL)*Power(xj,6LL) + 2LL*Power(r,7LL)*Power(xj,7LL)) - 

            Power(xi,14LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r,2LL)*Power(xj,2LL) + 

               1050LL*Power(r,3LL)*Power(xj,3LL) + 420LL*Power(r,4LL)*Power(xj,4LL) + 

               126LL*Power(r,5LL)*Power(xj,5LL) + 28LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL)) + 

            7LL*Power(xi,12LL)*Power(xj,2LL)*

             (1620LL + 2835LL*r*xj + 2430LL*Power(r,2LL)*Power(xj,2LL) + 

               1350LL*Power(r,3LL)*Power(xj,3LL) + 540LL*Power(r,4LL)*Power(xj,4LL) + 

               162LL*Power(r,5LL)*Power(xj,5LL) + 36LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL)) - 

            35LL*Power(xi,6LL)*Power(xj,8LL)*

             (4536LL + 7983LL*r*xj + 6534LL*Power(r,2LL)*Power(xj,2LL) + 

               4014LL*Power(r,3LL)*Power(xj,3LL) + 1644LL*Power(r,4LL)*Power(xj,4LL) + 

               414LL*Power(r,5LL)*Power(xj,5LL) + 60LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL)) + 

            21LL*Power(xi,4LL)*Power(xj,10LL)*

             (7920LL + 11385LL*r*xj + 12330LL*Power(r,2LL)*Power(xj,2LL) + 

               7410LL*Power(r,3LL)*Power(xj,3LL) + 2580LL*Power(r,4LL)*Power(xj,4LL) + 

               546LL*Power(r,5LL)*Power(xj,5LL) + 68LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL))))/

       (630LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,9LL)*Power(xi + xj,8LL)) - 

      (2520LL*exp(2LL*r*(xi + xj))*(xi + xj)*

          Power(Power(xi,2LL) - Power(xj,2LL),9LL) + 

         1260LL*exp(2LL*r*xj)*Power(xj,10LL)*

          (-Power(xi,9LL) - 6LL*Power(xi,7LL)*Power(xj,2LL) + 

            6LL*Power(xi,3LL)*Power(xj,6LL) + xi*Power(xj,8LL)) + 

         2520LL*exp(2LL*r*xj)*Power(xj,11LL)*

          (-6LL*Power(xi,8LL) - r*Power(xi,9LL) - 51LL*Power(xi,6LL)*Power(xj,2LL) - 

            6LL*r*Power(xi,7LL)*Power(xj,2LL) - 63LL*Power(xi,4LL)*Power(xj,4LL) - 

            9LL*Power(xi,2LL)*Power(xj,6LL) + 6LL*r*Power(xi,3LL)*Power(xj,6LL) + 

            Power(xj,8LL) + r*xi*Power(xj,8LL)) + 

         exp(2LL*r*xi)*Power(xi,4LL)*

          (-42LL*Power(xi,10LL)*Power(xj,4LL)*

             (1890LL*xj + 3240LL*r*Power(xj,2LL) + 2700LL*Power(r,2LL)*Power(xj,3LL) + 

               1440LL*Power(r,3LL)*Power(xj,4LL) + 555LL*Power(r,4LL)*Power(xj,5LL) + 

               132LL*Power(r,5LL)*Power(xj,6LL) + 14LL*Power(r,6LL)*Power(xj,7LL)) + 

            70LL*Power(xi,8LL)*Power(xj,6LL)*

             (2646LL*xj + 4536LL*r*Power(xj,2LL) + 3744LL*Power(r,2LL)*Power(xj,3LL) + 

               2112LL*Power(r,3LL)*Power(xj,4LL) + 765LL*Power(r,4LL)*Power(xj,5LL) + 

               156LL*Power(r,5LL)*Power(xj,6LL) + 14LL*Power(r,6LL)*Power(xj,7LL)) - 

            14LL*Power(xi,2LL)*Power(xj,12LL)*

             (16335LL*xj + 30780LL*r*Power(xj,2LL) + 21330LL*Power(r,2LL)*Power(xj,3LL) + 

               7920LL*Power(r,3LL)*Power(xj,4LL) + 1755LL*Power(r,4LL)*Power(xj,5LL) + 

               228LL*Power(r,5LL)*Power(xj,6LL) + 14LL*Power(r,6LL)*Power(xj,7LL)) + 

            2LL*Power(xj,14LL)*(72765LL*xj + 79380LL*r*Power(xj,2LL) + 

               39690LL*Power(r,2LL)*Power(xj,3LL) + 11760LL*Power(r,3LL)*Power(xj,4LL) + 

               2205LL*Power(r,4LL)*Power(xj,5LL) + 252LL*Power(r,5LL)*Power(xj,6LL) + 

               14LL*Power(r,6LL)*Power(xj,7LL)) - 

            Power(xi,14LL)*(2205LL*xj + 3780LL*r*Power(xj,2LL) + 

               3150LL*Power(r,2LL)*Power(xj,3LL) + 1680LL*Power(r,3LL)*Power(xj,4LL) + 

               630LL*Power(r,4LL)*Power(xj,5LL) + 168LL*Power(r,5LL)*Power(xj,6LL) + 

               28LL*Power(r,6LL)*Power(xj,7LL)) + 

            7LL*Power(xi,12LL)*Power(xj,2LL)*

             (2835LL*xj + 4860LL*r*Power(xj,2LL) + 4050LL*Power(r,2LL)*Power(xj,3LL) + 

               2160LL*Power(r,3LL)*Power(xj,4LL) + 810LL*Power(r,4LL)*Power(xj,5LL) + 

               216LL*Power(r,5LL)*Power(xj,6LL) + 28LL*Power(r,6LL)*Power(xj,7LL)) - 

            35LL*Power(xi,6LL)*Power(xj,8LL)*

             (7983LL*xj + 13068LL*r*Power(xj,2LL) + 12042LL*Power(r,2LL)*Power(xj,3LL) + 

               6576LL*Power(r,3LL)*Power(xj,4LL) + 2070LL*Power(r,4LL)*Power(xj,5LL) + 

               360LL*Power(r,5LL)*Power(xj,6LL) + 28LL*Power(r,6LL)*Power(xj,7LL)) + 

            21LL*Power(xi,4LL)*Power(xj,10LL)*

             (11385LL*xj + 24660LL*r*Power(xj,2LL) + 22230LL*Power(r,2LL)*Power(xj,3LL) + 

               10320LL*Power(r,3LL)*Power(xj,4LL) + 2730LL*Power(r,4LL)*Power(xj,5LL) + 

               408LL*Power(r,5LL)*Power(xj,6LL) + 28LL*Power(r,6LL)*Power(xj,7LL))) + 

         2LL*exp(2LL*r*xi)*Power(xi,5LL)*

          (-42LL*Power(xi,10LL)*Power(xj,4LL)*

             (1080LL + 1890LL*r*xj + 1620LL*Power(r,2LL)*Power(xj,2LL) + 

               900LL*Power(r,3LL)*Power(xj,3LL) + 360LL*Power(r,4LL)*Power(xj,4LL) + 

               111LL*Power(r,5LL)*Power(xj,5LL) + 22LL*Power(r,6LL)*Power(xj,6LL) + 

               2LL*Power(r,7LL)*Power(xj,7LL)) + 

            70LL*Power(xi,8LL)*Power(xj,6LL)*

             (1512LL + 2646LL*r*xj + 2268LL*Power(r,2LL)*Power(xj,2LL) + 

               1248LL*Power(r,3LL)*Power(xj,3LL) + 528LL*Power(r,4LL)*Power(xj,4LL) + 

               153LL*Power(r,5LL)*Power(xj,5LL) + 26LL*Power(r,6LL)*Power(xj,6LL) + 

               2LL*Power(r,7LL)*Power(xj,7LL)) - 

            14LL*Power(xi,2LL)*Power(xj,12LL)*

             (2970LL + 16335LL*r*xj + 15390LL*Power(r,2LL)*Power(xj,2LL) + 

               7110LL*Power(r,3LL)*Power(xj,3LL) + 1980LL*Power(r,4LL)*Power(xj,4LL) + 

               351LL*Power(r,5LL)*Power(xj,5LL) + 38LL*Power(r,6LL)*Power(xj,6LL) + 

               2LL*Power(r,7LL)*Power(xj,7LL)) + 

            2LL*Power(xj,14LL)*(62370LL + 72765LL*r*xj + 39690LL*Power(r,2LL)*Power(xj,2LL) + 

               13230LL*Power(r,3LL)*Power(xj,3LL) + 2940LL*Power(r,4LL)*Power(xj,4LL) + 

               441LL*Power(r,5LL)*Power(xj,5LL) + 42LL*Power(r,6LL)*Power(xj,6LL) + 

               2LL*Power(r,7LL)*Power(xj,7LL)) - 

            Power(xi,14LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r,2LL)*Power(xj,2LL) + 

               1050LL*Power(r,3LL)*Power(xj,3LL) + 420LL*Power(r,4LL)*Power(xj,4LL) + 

               126LL*Power(r,5LL)*Power(xj,5LL) + 28LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL)) + 

            7LL*Power(xi,12LL)*Power(xj,2LL)*

             (1620LL + 2835LL*r*xj + 2430LL*Power(r,2LL)*Power(xj,2LL) + 

               1350LL*Power(r,3LL)*Power(xj,3LL) + 540LL*Power(r,4LL)*Power(xj,4LL) + 

               162LL*Power(r,5LL)*Power(xj,5LL) + 36LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL)) - 

            35LL*Power(xi,6LL)*Power(xj,8LL)*

             (4536LL + 7983LL*r*xj + 6534LL*Power(r,2LL)*Power(xj,2LL) + 

               4014LL*Power(r,3LL)*Power(xj,3LL) + 1644LL*Power(r,4LL)*Power(xj,4LL) + 

               414LL*Power(r,5LL)*Power(xj,5LL) + 60LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL)) + 

            21LL*Power(xi,4LL)*Power(xj,10LL)*

             (7920LL + 11385LL*r*xj + 12330LL*Power(r,2LL)*Power(xj,2LL) + 

               7410LL*Power(r,3LL)*Power(xj,3LL) + 2580LL*Power(r,4LL)*Power(xj,4LL) + 

               546LL*Power(r,5LL)*Power(xj,5LL) + 68LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL))))/

       (1260LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,9LL)*Power(xi + xj,9LL))

    ; }
   
  }
  return S;
}

cl_R DSlater_1S_5S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = 0LL

    ; } else {  S = -(-2875101075LL*xi + 3193344000LL*exp(2LL*r*xi)*xi - 

          5113716300LL*r*Power(xi,2LL) - 4478789700LL*Power(r,2LL)*Power(xi,3LL) - 

          2564654400LL*Power(r,3LL)*Power(xi,4LL) - 

          1073595600LL*Power(r,4LL)*Power(xi,5LL) - 

          347276160LL*Power(r,5LL)*Power(xi,6LL) - 89147520LL*Power(r,6LL)*Power(xi,7LL) - 

          18247680LL*Power(r,7LL)*Power(xi,8LL) - 2914560LL*Power(r,8LL)*Power(xi,9LL) - 

          337920LL*Power(r,9LL)*Power(xi,10LL) - 22528LL*Power(r,10LL)*Power(xi,11LL))/

       (1.596672e9*exp(2LL*r*xi)*r) + 

      (-1596672000LL + 1596672000LL*exp(2LL*r*xi) - 2875101075LL*r*xi - 

         2556858150LL*Power(r,2LL)*Power(xi,2LL) - 

         1492929900LL*Power(r,3LL)*Power(xi,3LL) - 641163600LL*Power(r,4LL)*Power(xi,4LL) - 

         214719120LL*Power(r,5LL)*Power(xi,5LL) - 57879360LL*Power(r,6LL)*Power(xi,6LL) - 

         12735360LL*Power(r,7LL)*Power(xi,7LL) - 2280960LL*Power(r,8LL)*Power(xi,8LL) - 

         323840LL*Power(r,9LL)*Power(xi,9LL) - 33792LL*Power(r,10LL)*Power(xi,10LL) - 

         2048LL*Power(r,11LL)*Power(xi,11LL))/(1.596672e9*exp(2LL*r*xi)*Power(r,2LL)) + 

      (xi*(-1596672000LL + 1596672000LL*exp(2LL*r*xi) - 2875101075LL*r*xi - 

           2556858150LL*Power(r,2LL)*Power(xi,2LL) - 

           1492929900LL*Power(r,3LL)*Power(xi,3LL) - 

           641163600LL*Power(r,4LL)*Power(xi,4LL) - 214719120LL*Power(r,5LL)*Power(xi,5LL) - 

           57879360LL*Power(r,6LL)*Power(xi,6LL) - 12735360LL*Power(r,7LL)*Power(xi,7LL) - 

           2280960LL*Power(r,8LL)*Power(xi,8LL) - 323840LL*Power(r,9LL)*Power(xi,9LL) - 

           33792LL*Power(r,10LL)*Power(xi,10LL) - 2048LL*Power(r,11LL)*Power(xi,11LL)))/

       (7.98336e8*exp(2LL*r*xi)*r)

    ; }
 
  }
  else {
      if (r == 0LL) {  S = 0LL

    ; } else {  S = (14175LL*exp(2LL*r*(xi + xj))*Power(Power(xi,2LL) - Power(xj,2LL),11LL) + 

         2835LL*exp(2LL*r*xj)*Power(xj,12LL)*

          (-35LL*Power(xi,10LL) - 5LL*r*Power(xi,11LL) - 495LL*Power(xi,8LL)*Power(xj,2LL) - 

            55LL*r*Power(xi,9LL)*Power(xj,2LL) - 1254LL*Power(xi,6LL)*Power(xj,4LL) - 

            66LL*r*Power(xi,7LL)*Power(xj,4LL) - 726LL*Power(xi,4LL)*Power(xj,6LL) + 

            66LL*r*Power(xi,5LL)*Power(xj,6LL) - 55LL*Power(xi,2LL)*Power(xj,8LL) + 

            55LL*r*Power(xi,3LL)*Power(xj,8LL) + 5LL*Power(xj,10LL) + 5LL*r*xi*Power(xj,10LL)

    ) + exp(2LL*r*xi)*Power(xi,4LL)*

          (-(Power(xi,18LL)*(14175LL + 25515LL*r*xj + 22680LL*Power(r,2LL)*Power(xj,2LL) + 

                 13230LL*Power(r,3LL)*Power(xj,3LL) + 5670LL*Power(r,4LL)*Power(xj,4LL) + 

                 1890LL*Power(r,5LL)*Power(xj,5LL) + 504LL*Power(r,6LL)*Power(xj,6LL) + 

                 108LL*Power(r,7LL)*Power(xj,7LL) + 18LL*Power(r,8LL)*Power(xj,8LL) + 

                 2LL*Power(r,9LL)*Power(xj,9LL))) + 

            9LL*Power(xi,16LL)*Power(xj,2LL)*

             (17325LL + 31185LL*r*xj + 27720LL*Power(r,2LL)*Power(xj,2LL) + 

               16170LL*Power(r,3LL)*Power(xj,3LL) + 6930LL*Power(r,4LL)*Power(xj,4LL) + 

               2310LL*Power(r,5LL)*Power(xj,5LL) + 616LL*Power(r,6LL)*Power(xj,6LL) + 

               132LL*Power(r,7LL)*Power(xj,7LL) + 22LL*Power(r,8LL)*Power(xj,8LL) + 

               2LL*Power(r,9LL)*Power(xj,9LL)) - 

            126LL*Power(xi,10LL)*Power(xj,8LL)*

             (37125LL + 66825LL*r*xj + 59400LL*Power(r,2LL)*Power(xj,2LL) + 

               34725LL*Power(r,3LL)*Power(xj,3LL) + 14625LL*Power(r,4LL)*Power(xj,4LL) + 

               5043LL*Power(r,5LL)*Power(xj,5LL) + 1396LL*Power(r,6LL)*Power(xj,6LL) + 

               276LL*Power(r,7LL)*Power(xj,7LL) + 34LL*Power(r,8LL)*Power(xj,8LL) + 

               2LL*Power(r,9LL)*Power(xj,9LL)) + 

            126LL*Power(xi,8LL)*Power(xj,10LL)*

             (51975LL + 93420LL*r*xj + 84240LL*Power(r,2LL)*Power(xj,2LL) + 

               46815LL*Power(r,3LL)*Power(xj,3LL) + 20835LL*Power(r,4LL)*Power(xj,4LL) + 

               7485LL*Power(r,5LL)*Power(xj,5LL) + 1964LL*Power(r,6LL)*Power(xj,6LL) + 

               348LL*Power(r,7LL)*Power(xj,7LL) + 38LL*Power(r,8LL)*Power(xj,8LL) + 

               2LL*Power(r,9LL)*Power(xj,9LL)) - 

            9LL*Power(xi,2LL)*Power(xj,16LL)*

             (-135135LL + 405405LL*r*xj + 582120LL*Power(r,2LL)*Power(xj,2LL) + 

               346500LL*Power(r,3LL)*Power(xj,3LL) + 124740LL*Power(r,4LL)*Power(xj,4LL) + 

               30492LL*Power(r,5LL)*Power(xj,5LL) + 5264LL*Power(r,6LL)*Power(xj,6LL) + 

               636LL*Power(r,7LL)*Power(xj,7LL) + 50LL*Power(r,8LL)*Power(xj,8LL) + 

               2LL*Power(r,9LL)*Power(xj,9LL)) + 

            Power(xj,18LL)*(2837835LL + 3648645LL*r*xj + 

               2245320LL*Power(r,2LL)*Power(xj,2LL) + 

               873180LL*Power(r,3LL)*Power(xj,3LL) + 238140LL*Power(r,4LL)*Power(xj,4LL) + 

               47628LL*Power(r,5LL)*Power(xj,5LL) + 7056LL*Power(r,6LL)*Power(xj,6LL) + 

               756LL*Power(r,7LL)*Power(xj,7LL) + 54LL*Power(r,8LL)*Power(xj,8LL) + 

               2LL*Power(r,9LL)*Power(xj,9LL)) - 

            9LL*Power(xi,14LL)*Power(xj,4LL)*

             (86625LL + 155925LL*r*xj + 138600LL*Power(r,2LL)*Power(xj,2LL) + 

               80850LL*Power(r,3LL)*Power(xj,3LL) + 34650LL*Power(r,4LL)*Power(xj,4LL) + 

               11550LL*Power(r,5LL)*Power(xj,5LL) + 3080LL*Power(r,6LL)*Power(xj,6LL) + 

               672LL*Power(r,7LL)*Power(xj,7LL) + 104LL*Power(r,8LL)*Power(xj,8LL) + 

               8LL*Power(r,9LL)*Power(xj,9LL)) + 

            21LL*Power(xi,12LL)*Power(xj,6LL)*

             (111375LL + 200475LL*r*xj + 178200LL*Power(r,2LL)*Power(xj,2LL) + 

               103950LL*Power(r,3LL)*Power(xj,3LL) + 44550LL*Power(r,4LL)*Power(xj,4LL) + 

               14778LL*Power(r,5LL)*Power(xj,5LL) + 4056LL*Power(r,6LL)*Power(xj,6LL) + 

               864LL*Power(r,7LL)*Power(xj,7LL) + 120LL*Power(r,8LL)*Power(xj,8LL) + 

               8LL*Power(r,9LL)*Power(xj,9LL)) - 

            21LL*Power(xi,6LL)*Power(xj,12LL)*

             (307125LL + 594945LL*r*xj + 456840LL*Power(r,2LL)*Power(xj,2LL) + 

               281790LL*Power(r,3LL)*Power(xj,3LL) + 137430LL*Power(r,4LL)*Power(xj,4LL) + 

               47250LL*Power(r,5LL)*Power(xj,5LL) + 11064LL*Power(r,6LL)*Power(xj,6LL) + 

               1728LL*Power(r,7LL)*Power(xj,7LL) + 168LL*Power(r,8LL)*Power(xj,8LL) + 

               8LL*Power(r,9LL)*Power(xj,9LL)) + 

            9LL*Power(xi,4LL)*Power(xj,14LL)*

             (675675LL + 675675LL*r*xj + 748440LL*Power(r,2LL)*Power(xj,2LL) + 

               561330LL*Power(r,3LL)*Power(xj,3LL) + 256410LL*Power(r,4LL)*Power(xj,4LL) + 

               76230LL*Power(r,5LL)*Power(xj,5LL) + 15400LL*Power(r,6LL)*Power(xj,6LL) + 

               2112LL*Power(r,7LL)*Power(xj,7LL) + 184LL*Power(r,8LL)*Power(xj,8LL) + 

               8LL*Power(r,9LL)*Power(xj,9LL))))/

       (14175LL*exp(2LL*r*(xi + xj))*Power(r,2LL)*Power(xi - xj,11LL)*

         Power(xi + xj,11LL)) + (2LL*(14175LL*exp(2LL*r*(xi + xj))*

            Power(Power(xi,2LL) - Power(xj,2LL),11LL) + 

           2835LL*exp(2LL*r*xj)*Power(xj,12LL)*

            (-35LL*Power(xi,10LL) - 5LL*r*Power(xi,11LL) - 

              495LL*Power(xi,8LL)*Power(xj,2LL) - 55LL*r*Power(xi,9LL)*Power(xj,2LL) - 

              1254LL*Power(xi,6LL)*Power(xj,4LL) - 66LL*r*Power(xi,7LL)*Power(xj,4LL) - 

              726LL*Power(xi,4LL)*Power(xj,6LL) + 66LL*r*Power(xi,5LL)*Power(xj,6LL) - 

              55LL*Power(xi,2LL)*Power(xj,8LL) + 55LL*r*Power(xi,3LL)*Power(xj,8LL) + 

              5LL*Power(xj,10LL) + 5LL*r*xi*Power(xj,10LL)) + 

           exp(2LL*r*xi)*Power(xi,4LL)*

            (-(Power(xi,18LL)*(14175LL + 25515LL*r*xj + 

                   22680LL*Power(r,2LL)*Power(xj,2LL) + 

                   13230LL*Power(r,3LL)*Power(xj,3LL) + 

                   5670LL*Power(r,4LL)*Power(xj,4LL) + 1890LL*Power(r,5LL)*Power(xj,5LL) + 

                   504LL*Power(r,6LL)*Power(xj,6LL) + 108LL*Power(r,7LL)*Power(xj,7LL) + 

                   18LL*Power(r,8LL)*Power(xj,8LL) + 2LL*Power(r,9LL)*Power(xj,9LL))) + 

              9LL*Power(xi,16LL)*Power(xj,2LL)*

               (17325LL + 31185LL*r*xj + 27720LL*Power(r,2LL)*Power(xj,2LL) + 

                 16170LL*Power(r,3LL)*Power(xj,3LL) + 6930LL*Power(r,4LL)*Power(xj,4LL) + 

                 2310LL*Power(r,5LL)*Power(xj,5LL) + 616LL*Power(r,6LL)*Power(xj,6LL) + 

                 132LL*Power(r,7LL)*Power(xj,7LL) + 22LL*Power(r,8LL)*Power(xj,8LL) + 

                 2LL*Power(r,9LL)*Power(xj,9LL)) - 

              126LL*Power(xi,10LL)*Power(xj,8LL)*

               (37125LL + 66825LL*r*xj + 59400LL*Power(r,2LL)*Power(xj,2LL) + 

                 34725LL*Power(r,3LL)*Power(xj,3LL) + 14625LL*Power(r,4LL)*Power(xj,4LL) + 

                 5043LL*Power(r,5LL)*Power(xj,5LL) + 1396LL*Power(r,6LL)*Power(xj,6LL) + 

                 276LL*Power(r,7LL)*Power(xj,7LL) + 34LL*Power(r,8LL)*Power(xj,8LL) + 

                 2LL*Power(r,9LL)*Power(xj,9LL)) + 

              126LL*Power(xi,8LL)*Power(xj,10LL)*

               (51975LL + 93420LL*r*xj + 84240LL*Power(r,2LL)*Power(xj,2LL) + 

                 46815LL*Power(r,3LL)*Power(xj,3LL) + 20835LL*Power(r,4LL)*Power(xj,4LL) + 

                 7485LL*Power(r,5LL)*Power(xj,5LL) + 1964LL*Power(r,6LL)*Power(xj,6LL) + 

                 348LL*Power(r,7LL)*Power(xj,7LL) + 38LL*Power(r,8LL)*Power(xj,8LL) + 

                 2LL*Power(r,9LL)*Power(xj,9LL)) - 

              9LL*Power(xi,2LL)*Power(xj,16LL)*

               (-135135LL + 405405LL*r*xj + 582120LL*Power(r,2LL)*Power(xj,2LL) + 

                 346500LL*Power(r,3LL)*Power(xj,3LL) + 

                 124740LL*Power(r,4LL)*Power(xj,4LL) + 

                 30492LL*Power(r,5LL)*Power(xj,5LL) + 5264LL*Power(r,6LL)*Power(xj,6LL) + 

                 636LL*Power(r,7LL)*Power(xj,7LL) + 50LL*Power(r,8LL)*Power(xj,8LL) + 

                 2LL*Power(r,9LL)*Power(xj,9LL)) + 

              Power(xj,18LL)*(2837835LL + 3648645LL*r*xj + 

                 2245320LL*Power(r,2LL)*Power(xj,2LL) + 

                 873180LL*Power(r,3LL)*Power(xj,3LL) + 

                 238140LL*Power(r,4LL)*Power(xj,4LL) + 

                 47628LL*Power(r,5LL)*Power(xj,5LL) + 7056LL*Power(r,6LL)*Power(xj,6LL) + 

                 756LL*Power(r,7LL)*Power(xj,7LL) + 54LL*Power(r,8LL)*Power(xj,8LL) + 

                 2LL*Power(r,9LL)*Power(xj,9LL)) - 

              9LL*Power(xi,14LL)*Power(xj,4LL)*

               (86625LL + 155925LL*r*xj + 138600LL*Power(r,2LL)*Power(xj,2LL) + 

                 80850LL*Power(r,3LL)*Power(xj,3LL) + 34650LL*Power(r,4LL)*Power(xj,4LL) + 

                 11550LL*Power(r,5LL)*Power(xj,5LL) + 3080LL*Power(r,6LL)*Power(xj,6LL) + 

                 672LL*Power(r,7LL)*Power(xj,7LL) + 104LL*Power(r,8LL)*Power(xj,8LL) + 

                 8LL*Power(r,9LL)*Power(xj,9LL)) + 

              21LL*Power(xi,12LL)*Power(xj,6LL)*

               (111375LL + 200475LL*r*xj + 178200LL*Power(r,2LL)*Power(xj,2LL) + 

                 103950LL*Power(r,3LL)*Power(xj,3LL) + 

                 44550LL*Power(r,4LL)*Power(xj,4LL) + 14778LL*Power(r,5LL)*Power(xj,5LL) + 

                 4056LL*Power(r,6LL)*Power(xj,6LL) + 864LL*Power(r,7LL)*Power(xj,7LL) + 

                 120LL*Power(r,8LL)*Power(xj,8LL) + 8LL*Power(r,9LL)*Power(xj,9LL)) - 

              21LL*Power(xi,6LL)*Power(xj,12LL)*

               (307125LL + 594945LL*r*xj + 456840LL*Power(r,2LL)*Power(xj,2LL) + 

                 281790LL*Power(r,3LL)*Power(xj,3LL) + 

                 137430LL*Power(r,4LL)*Power(xj,4LL) + 

                 47250LL*Power(r,5LL)*Power(xj,5LL) + 11064LL*Power(r,6LL)*Power(xj,6LL) + 

                 1728LL*Power(r,7LL)*Power(xj,7LL) + 168LL*Power(r,8LL)*Power(xj,8LL) + 

                 8LL*Power(r,9LL)*Power(xj,9LL)) + 

              9LL*Power(xi,4LL)*Power(xj,14LL)*

               (675675LL + 675675LL*r*xj + 748440LL*Power(r,2LL)*Power(xj,2LL) + 

                 561330LL*Power(r,3LL)*Power(xj,3LL) + 

                 256410LL*Power(r,4LL)*Power(xj,4LL) + 76230LL*Power(r,5LL)*Power(xj,5LL) + 

                 15400LL*Power(r,6LL)*Power(xj,6LL) + 2112LL*Power(r,7LL)*Power(xj,7LL) + 

                 184LL*Power(r,8LL)*Power(xj,8LL) + 8LL*Power(r,9LL)*Power(xj,9LL)))))/

       (14175LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,11LL)*Power(xi + xj,10LL)) - 

      (28350LL*exp(2LL*r*(xi + xj))*(xi + xj)*

          Power(Power(xi,2LL) - Power(xj,2LL),11LL) + 

         2835LL*exp(2LL*r*xj)*Power(xj,12LL)*

          (-5LL*Power(xi,11LL) - 55LL*Power(xi,9LL)*Power(xj,2LL) - 

            66LL*Power(xi,7LL)*Power(xj,4LL) + 66LL*Power(xi,5LL)*Power(xj,6LL) + 

            55LL*Power(xi,3LL)*Power(xj,8LL) + 5LL*xi*Power(xj,10LL)) + 

         5670LL*exp(2LL*r*xj)*Power(xj,13LL)*

          (-35LL*Power(xi,10LL) - 5LL*r*Power(xi,11LL) - 495LL*Power(xi,8LL)*Power(xj,2LL) - 

            55LL*r*Power(xi,9LL)*Power(xj,2LL) - 1254LL*Power(xi,6LL)*Power(xj,4LL) - 

            66LL*r*Power(xi,7LL)*Power(xj,4LL) - 726LL*Power(xi,4LL)*Power(xj,6LL) + 

            66LL*r*Power(xi,5LL)*Power(xj,6LL) - 55LL*Power(xi,2LL)*Power(xj,8LL) + 

            55LL*r*Power(xi,3LL)*Power(xj,8LL) + 5LL*Power(xj,10LL) + 5LL*r*xi*Power(xj,10LL)) 

    + exp(2LL*r*xi)*Power(xi,4LL)*(-(Power(xi,18LL)*

               (25515LL*xj + 45360LL*r*Power(xj,2LL) + 

                 39690LL*Power(r,2LL)*Power(xj,3LL) + 22680LL*Power(r,3LL)*Power(xj,4LL) + 

                 9450LL*Power(r,4LL)*Power(xj,5LL) + 3024LL*Power(r,5LL)*Power(xj,6LL) + 

                 756LL*Power(r,6LL)*Power(xj,7LL) + 144LL*Power(r,7LL)*Power(xj,8LL) + 

                 18LL*Power(r,8LL)*Power(xj,9LL))) + 

            9LL*Power(xi,16LL)*Power(xj,2LL)*

             (31185LL*xj + 55440LL*r*Power(xj,2LL) + 48510LL*Power(r,2LL)*Power(xj,3LL) + 

               27720LL*Power(r,3LL)*Power(xj,4LL) + 11550LL*Power(r,4LL)*Power(xj,5LL) + 

               3696LL*Power(r,5LL)*Power(xj,6LL) + 924LL*Power(r,6LL)*Power(xj,7LL) + 

               176LL*Power(r,7LL)*Power(xj,8LL) + 18LL*Power(r,8LL)*Power(xj,9LL)) - 

            126LL*Power(xi,10LL)*Power(xj,8LL)*

             (66825LL*xj + 118800LL*r*Power(xj,2LL) + 

               104175LL*Power(r,2LL)*Power(xj,3LL) + 58500LL*Power(r,3LL)*Power(xj,4LL) + 

               25215LL*Power(r,4LL)*Power(xj,5LL) + 8376LL*Power(r,5LL)*Power(xj,6LL) + 

               1932LL*Power(r,6LL)*Power(xj,7LL) + 272LL*Power(r,7LL)*Power(xj,8LL) + 

               18LL*Power(r,8LL)*Power(xj,9LL)) + 

            126LL*Power(xi,8LL)*Power(xj,10LL)*

             (93420LL*xj + 168480LL*r*Power(xj,2LL) + 

               140445LL*Power(r,2LL)*Power(xj,3LL) + 83340LL*Power(r,3LL)*Power(xj,4LL) + 

               37425LL*Power(r,4LL)*Power(xj,5LL) + 11784LL*Power(r,5LL)*Power(xj,6LL) + 

               2436LL*Power(r,6LL)*Power(xj,7LL) + 304LL*Power(r,7LL)*Power(xj,8LL) + 

               18LL*Power(r,8LL)*Power(xj,9LL)) - 

            9LL*Power(xi,2LL)*Power(xj,16LL)*

             (405405LL*xj + 1164240LL*r*Power(xj,2LL) + 

               1039500LL*Power(r,2LL)*Power(xj,3LL) + 

               498960LL*Power(r,3LL)*Power(xj,4LL) + 152460LL*Power(r,4LL)*Power(xj,5LL) + 

               31584LL*Power(r,5LL)*Power(xj,6LL) + 4452LL*Power(r,6LL)*Power(xj,7LL) + 

               400LL*Power(r,7LL)*Power(xj,8LL) + 18LL*Power(r,8LL)*Power(xj,9LL)) + 

            Power(xj,18LL)*(3648645LL*xj + 4490640LL*r*Power(xj,2LL) + 

               2619540LL*Power(r,2LL)*Power(xj,3LL) + 

               952560LL*Power(r,3LL)*Power(xj,4LL) + 238140LL*Power(r,4LL)*Power(xj,5LL) + 

               42336LL*Power(r,5LL)*Power(xj,6LL) + 5292LL*Power(r,6LL)*Power(xj,7LL) + 

               432LL*Power(r,7LL)*Power(xj,8LL) + 18LL*Power(r,8LL)*Power(xj,9LL)) - 

            9LL*Power(xi,14LL)*Power(xj,4LL)*

             (155925LL*xj + 277200LL*r*Power(xj,2LL) + 

               242550LL*Power(r,2LL)*Power(xj,3LL) + 138600LL*Power(r,3LL)*Power(xj,4LL) + 

               57750LL*Power(r,4LL)*Power(xj,5LL) + 18480LL*Power(r,5LL)*Power(xj,6LL) + 

               4704LL*Power(r,6LL)*Power(xj,7LL) + 832LL*Power(r,7LL)*Power(xj,8LL) + 

               72LL*Power(r,8LL)*Power(xj,9LL)) + 

            21LL*Power(xi,12LL)*Power(xj,6LL)*

             (200475LL*xj + 356400LL*r*Power(xj,2LL) + 

               311850LL*Power(r,2LL)*Power(xj,3LL) + 178200LL*Power(r,3LL)*Power(xj,4LL) + 

               73890LL*Power(r,4LL)*Power(xj,5LL) + 24336LL*Power(r,5LL)*Power(xj,6LL) + 

               6048LL*Power(r,6LL)*Power(xj,7LL) + 960LL*Power(r,7LL)*Power(xj,8LL) + 

               72LL*Power(r,8LL)*Power(xj,9LL)) - 

            21LL*Power(xi,6LL)*Power(xj,12LL)*

             (594945LL*xj + 913680LL*r*Power(xj,2LL) + 

               845370LL*Power(r,2LL)*Power(xj,3LL) + 549720LL*Power(r,3LL)*Power(xj,4LL) + 

               236250LL*Power(r,4LL)*Power(xj,5LL) + 66384LL*Power(r,5LL)*Power(xj,6LL) + 

               12096LL*Power(r,6LL)*Power(xj,7LL) + 1344LL*Power(r,7LL)*Power(xj,8LL) + 

               72LL*Power(r,8LL)*Power(xj,9LL)) + 

            9LL*Power(xi,4LL)*Power(xj,14LL)*

             (675675LL*xj + 1496880LL*r*Power(xj,2LL) + 

               1683990LL*Power(r,2LL)*Power(xj,3LL) + 

               1025640LL*Power(r,3LL)*Power(xj,4LL) + 381150LL*Power(r,4LL)*Power(xj,5LL) + 

               92400LL*Power(r,5LL)*Power(xj,6LL) + 14784LL*Power(r,6LL)*Power(xj,7LL) + 

               1472LL*Power(r,7LL)*Power(xj,8LL) + 72LL*Power(r,8LL)*Power(xj,9LL))) + 

         2LL*exp(2LL*r*xi)*Power(xi,5LL)*

          (-(Power(xi,18LL)*(14175LL + 25515LL*r*xj + 22680LL*Power(r,2LL)*Power(xj,2LL) + 

                 13230LL*Power(r,3LL)*Power(xj,3LL) + 5670LL*Power(r,4LL)*Power(xj,4LL) + 

                 1890LL*Power(r,5LL)*Power(xj,5LL) + 504LL*Power(r,6LL)*Power(xj,6LL) + 

                 108LL*Power(r,7LL)*Power(xj,7LL) + 18LL*Power(r,8LL)*Power(xj,8LL) + 

                 2LL*Power(r,9LL)*Power(xj,9LL))) + 

            9LL*Power(xi,16LL)*Power(xj,2LL)*

             (17325LL + 31185LL*r*xj + 27720LL*Power(r,2LL)*Power(xj,2LL) + 

               16170LL*Power(r,3LL)*Power(xj,3LL) + 6930LL*Power(r,4LL)*Power(xj,4LL) + 

               2310LL*Power(r,5LL)*Power(xj,5LL) + 616LL*Power(r,6LL)*Power(xj,6LL) + 

               132LL*Power(r,7LL)*Power(xj,7LL) + 22LL*Power(r,8LL)*Power(xj,8LL) + 

               2LL*Power(r,9LL)*Power(xj,9LL)) - 

            126LL*Power(xi,10LL)*Power(xj,8LL)*

             (37125LL + 66825LL*r*xj + 59400LL*Power(r,2LL)*Power(xj,2LL) + 

               34725LL*Power(r,3LL)*Power(xj,3LL) + 14625LL*Power(r,4LL)*Power(xj,4LL) + 

               5043LL*Power(r,5LL)*Power(xj,5LL) + 1396LL*Power(r,6LL)*Power(xj,6LL) + 

               276LL*Power(r,7LL)*Power(xj,7LL) + 34LL*Power(r,8LL)*Power(xj,8LL) + 

               2LL*Power(r,9LL)*Power(xj,9LL)) + 

            126LL*Power(xi,8LL)*Power(xj,10LL)*

             (51975LL + 93420LL*r*xj + 84240LL*Power(r,2LL)*Power(xj,2LL) + 

               46815LL*Power(r,3LL)*Power(xj,3LL) + 20835LL*Power(r,4LL)*Power(xj,4LL) + 

               7485LL*Power(r,5LL)*Power(xj,5LL) + 1964LL*Power(r,6LL)*Power(xj,6LL) + 

               348LL*Power(r,7LL)*Power(xj,7LL) + 38LL*Power(r,8LL)*Power(xj,8LL) + 

               2LL*Power(r,9LL)*Power(xj,9LL)) - 

            9LL*Power(xi,2LL)*Power(xj,16LL)*

             (-135135LL + 405405LL*r*xj + 582120LL*Power(r,2LL)*Power(xj,2LL) + 

               346500LL*Power(r,3LL)*Power(xj,3LL) + 124740LL*Power(r,4LL)*Power(xj,4LL) + 

               30492LL*Power(r,5LL)*Power(xj,5LL) + 5264LL*Power(r,6LL)*Power(xj,6LL) + 

               636LL*Power(r,7LL)*Power(xj,7LL) + 50LL*Power(r,8LL)*Power(xj,8LL) + 

               2LL*Power(r,9LL)*Power(xj,9LL)) + 

            Power(xj,18LL)*(2837835LL + 3648645LL*r*xj + 

               2245320LL*Power(r,2LL)*Power(xj,2LL) + 873180LL*Power(r,3LL)*Power(xj,3LL) + 

               238140LL*Power(r,4LL)*Power(xj,4LL) + 47628LL*Power(r,5LL)*Power(xj,5LL) + 

               7056LL*Power(r,6LL)*Power(xj,6LL) + 756LL*Power(r,7LL)*Power(xj,7LL) + 

               54LL*Power(r,8LL)*Power(xj,8LL) + 2LL*Power(r,9LL)*Power(xj,9LL)) - 

            9LL*Power(xi,14LL)*Power(xj,4LL)*

             (86625LL + 155925LL*r*xj + 138600LL*Power(r,2LL)*Power(xj,2LL) + 

               80850LL*Power(r,3LL)*Power(xj,3LL) + 34650LL*Power(r,4LL)*Power(xj,4LL) + 

               11550LL*Power(r,5LL)*Power(xj,5LL) + 3080LL*Power(r,6LL)*Power(xj,6LL) + 

               672LL*Power(r,7LL)*Power(xj,7LL) + 104LL*Power(r,8LL)*Power(xj,8LL) + 

               8LL*Power(r,9LL)*Power(xj,9LL)) + 

            21LL*Power(xi,12LL)*Power(xj,6LL)*

             (111375LL + 200475LL*r*xj + 178200LL*Power(r,2LL)*Power(xj,2LL) + 

               103950LL*Power(r,3LL)*Power(xj,3LL) + 44550LL*Power(r,4LL)*Power(xj,4LL) + 

               14778LL*Power(r,5LL)*Power(xj,5LL) + 4056LL*Power(r,6LL)*Power(xj,6LL) + 

               864LL*Power(r,7LL)*Power(xj,7LL) + 120LL*Power(r,8LL)*Power(xj,8LL) + 

               8LL*Power(r,9LL)*Power(xj,9LL)) - 

            21LL*Power(xi,6LL)*Power(xj,12LL)*

             (307125LL + 594945LL*r*xj + 456840LL*Power(r,2LL)*Power(xj,2LL) + 

               281790LL*Power(r,3LL)*Power(xj,3LL) + 137430LL*Power(r,4LL)*Power(xj,4LL) + 

               47250LL*Power(r,5LL)*Power(xj,5LL) + 11064LL*Power(r,6LL)*Power(xj,6LL) + 

               1728LL*Power(r,7LL)*Power(xj,7LL) + 168LL*Power(r,8LL)*Power(xj,8LL) + 

               8LL*Power(r,9LL)*Power(xj,9LL)) + 

            9LL*Power(xi,4LL)*Power(xj,14LL)*

             (675675LL + 675675LL*r*xj + 748440LL*Power(r,2LL)*Power(xj,2LL) + 

               561330LL*Power(r,3LL)*Power(xj,3LL) + 256410LL*Power(r,4LL)*Power(xj,4LL) + 

               76230LL*Power(r,5LL)*Power(xj,5LL) + 15400LL*Power(r,6LL)*Power(xj,6LL) + 

               2112LL*Power(r,7LL)*Power(xj,7LL) + 184LL*Power(r,8LL)*Power(xj,8LL) + 

               8LL*Power(r,9LL)*Power(xj,9LL))))/

       (14175LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,11LL)*Power(xi + xj,11LL))

    ; }
   
  }
  return S;
}

cl_R DSlater_2S_2S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = 0LL

    ; } else {  S = -(-131985LL*xi + 161280LL*exp(2LL*r*xi)*xi - 205380LL*r*Power(xi,2LL) - 

          149940LL*Power(r,2LL)*Power(xi,3LL) - 67200LL*Power(r,3LL)*Power(xi,4LL) - 

          20160LL*Power(r,4LL)*Power(xi,5LL) - 4032LL*Power(r,5LL)*Power(xi,6LL) - 

          448LL*Power(r,6LL)*Power(xi,7LL))/(80640LL*exp(2LL*r*xi)*r) + 

      (-80640LL + 80640LL*exp(2LL*r*xi) - 131985LL*r*xi - 

         102690LL*Power(r,2LL)*Power(xi,2LL) - 49980LL*Power(r,3LL)*Power(xi,3LL) - 

         16800LL*Power(r,4LL)*Power(xi,4LL) - 4032LL*Power(r,5LL)*Power(xi,5LL) - 

         672LL*Power(r,6LL)*Power(xi,6LL) - 64LL*Power(r,7LL)*Power(xi,7LL))/

       (80640LL*exp(2LL*r*xi)*Power(r,2LL)) + 

      (xi*(-80640LL + 80640LL*exp(2LL*r*xi) - 131985LL*r*xi - 

           102690LL*Power(r,2LL)*Power(xi,2LL) - 49980LL*Power(r,3LL)*Power(xi,3LL) - 

           16800LL*Power(r,4LL)*Power(xi,4LL) - 4032LL*Power(r,5LL)*Power(xi,5LL) - 

           672LL*Power(r,6LL)*Power(xi,6LL) - 64LL*Power(r,7LL)*Power(xi,7LL)))/

       (40320LL*exp(2LL*r*xi)*r)

    ; }
 
  }
  else {
      if (r == 0LL) {  S = 0LL

    ; } else {  S = (6LL*exp(2LL*r*(xi + xj))*Power(Power(xi,2LL) - Power(xj,2LL),7LL) - 

         exp(2LL*r*xi)*Power(xi,6LL)*

          (21LL*Power(xi,4LL)*Power(xj,4LL)*

             (6LL + 11LL*r*xj + 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            2LL*Power(xj,8LL)*(90LL + 54LL*r*xj + 12LL*Power(r,2LL)*Power(xj,2LL) + 

               Power(r,3LL)*Power(xj,3LL)) + 

            Power(xi,8LL)*(6LL + 9LL*r*xj + 6LL*Power(r,2LL)*Power(xj,2LL) + 

               2LL*Power(r,3LL)*Power(xj,3LL)) + 

            Power(xi,2LL)*Power(xj,6LL)*

             (-390LL - 69LL*r*xj + 18LL*Power(r,2LL)*Power(xj,2LL) + 

               4LL*Power(r,3LL)*Power(xj,3LL)) - 

            Power(xi,6LL)*Power(xj,2LL)*

             (42LL + 63LL*r*xj + 42LL*Power(r,2LL)*Power(xj,2LL) + 

               4LL*Power(r,3LL)*Power(xj,3LL))) + 

         exp(2LL*r*xj)*Power(xj,6LL)*

          (-24LL*Power(r,2LL)*Power(xi,10LL) - 2LL*Power(r,3LL)*Power(xi,11LL) - 

            69LL*r*Power(xi,7LL)*Power(xj,2LL) + 6LL*Power(xj,8LL) + 9LL*r*xi*Power(xj,8LL) + 

            4LL*r*Power(xi,9LL)*(-27LL + Power(r,2LL)*Power(xj,2LL)) + 

            18LL*Power(xi,8LL)*(-10LL + Power(r,2LL)*Power(xj,2LL)) + 

            6LL*Power(xi,2LL)*Power(xj,6LL)*(-7LL + Power(r,2LL)*Power(xj,2LL)) - 

            42LL*Power(xi,4LL)*Power(xj,4LL)*(-3LL + Power(r,2LL)*Power(xj,2LL)) + 

            r*Power(xi,3LL)*Power(xj,6LL)*(-63LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            6LL*Power(xi,6LL)*Power(xj,2LL)*(-65LL + 7LL*Power(r,2LL)*Power(xj,2LL)) + 

            Power(xi,5LL)*(231LL*r*Power(xj,4LL) - 4LL*Power(r,3LL)*Power(xj,6LL))))/

       (6LL*exp(2LL*r*(xi + xj))*Power(r,2LL)*Power(xi - xj,7LL)*Power(xi + xj,7LL)) + 

      (6LL*exp(2LL*r*(xi + xj))*Power(Power(xi,2LL) - Power(xj,2LL),7LL) - 

         exp(2LL*r*xi)*Power(xi,6LL)*

          (21LL*Power(xi,4LL)*Power(xj,4LL)*

             (6LL + 11LL*r*xj + 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            2LL*Power(xj,8LL)*(90LL + 54LL*r*xj + 12LL*Power(r,2LL)*Power(xj,2LL) + 

               Power(r,3LL)*Power(xj,3LL)) + 

            Power(xi,8LL)*(6LL + 9LL*r*xj + 6LL*Power(r,2LL)*Power(xj,2LL) + 

               2LL*Power(r,3LL)*Power(xj,3LL)) + 

            Power(xi,2LL)*Power(xj,6LL)*

             (-390LL - 69LL*r*xj + 18LL*Power(r,2LL)*Power(xj,2LL) + 

               4LL*Power(r,3LL)*Power(xj,3LL)) - 

            Power(xi,6LL)*Power(xj,2LL)*

             (42LL + 63LL*r*xj + 42LL*Power(r,2LL)*Power(xj,2LL) + 

               4LL*Power(r,3LL)*Power(xj,3LL))) + 

         exp(2LL*r*xj)*Power(xj,6LL)*

          (-24LL*Power(r,2LL)*Power(xi,10LL) - 2LL*Power(r,3LL)*Power(xi,11LL) - 

            69LL*r*Power(xi,7LL)*Power(xj,2LL) + 6LL*Power(xj,8LL) + 9LL*r*xi*Power(xj,8LL) + 

            4LL*r*Power(xi,9LL)*(-27LL + Power(r,2LL)*Power(xj,2LL)) + 

            18LL*Power(xi,8LL)*(-10LL + Power(r,2LL)*Power(xj,2LL)) + 

            6LL*Power(xi,2LL)*Power(xj,6LL)*(-7LL + Power(r,2LL)*Power(xj,2LL)) - 

            42LL*Power(xi,4LL)*Power(xj,4LL)*(-3LL + Power(r,2LL)*Power(xj,2LL)) + 

            r*Power(xi,3LL)*Power(xj,6LL)*(-63LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            6LL*Power(xi,6LL)*Power(xj,2LL)*(-65LL + 7LL*Power(r,2LL)*Power(xj,2LL)) + 

            Power(xi,5LL)*(231LL*r*Power(xj,4LL) - 4LL*Power(r,3LL)*Power(xj,6LL))))/

       (3LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,7LL)*Power(xi + xj,6LL)) - 

      (12LL*exp(2LL*r*(xi + xj))*(xi + xj)*Power(Power(xi,2LL) - Power(xj,2LL),7LL) - 

         exp(2LL*r*xi)*Power(xi,6LL)*

          (21LL*Power(xi,4LL)*Power(xj,4LL)*(11LL*xj + 4LL*r*Power(xj,2LL)) - 

            2LL*Power(xj,8LL)*(54LL*xj + 24LL*r*Power(xj,2LL) + 

               3LL*Power(r,2LL)*Power(xj,3LL)) + 

            Power(xi,8LL)*(9LL*xj + 12LL*r*Power(xj,2LL) + 6LL*Power(r,2LL)*Power(xj,3LL)) + 

            Power(xi,2LL)*Power(xj,6LL)*

             (-69LL*xj + 36LL*r*Power(xj,2LL) + 12LL*Power(r,2LL)*Power(xj,3LL)) - 

            Power(xi,6LL)*Power(xj,2LL)*

             (63LL*xj + 84LL*r*Power(xj,2LL) + 12LL*Power(r,2LL)*Power(xj,3LL))) - 

         2LL*exp(2LL*r*xi)*Power(xi,7LL)*

          (21LL*Power(xi,4LL)*Power(xj,4LL)*

             (6LL + 11LL*r*xj + 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            2LL*Power(xj,8LL)*(90LL + 54LL*r*xj + 12LL*Power(r,2LL)*Power(xj,2LL) + 

               Power(r,3LL)*Power(xj,3LL)) + 

            Power(xi,8LL)*(6LL + 9LL*r*xj + 6LL*Power(r,2LL)*Power(xj,2LL) + 

               2LL*Power(r,3LL)*Power(xj,3LL)) + 

            Power(xi,2LL)*Power(xj,6LL)*

             (-390LL - 69LL*r*xj + 18LL*Power(r,2LL)*Power(xj,2LL) + 

               4LL*Power(r,3LL)*Power(xj,3LL)) - 

            Power(xi,6LL)*Power(xj,2LL)*

             (42LL + 63LL*r*xj + 42LL*Power(r,2LL)*Power(xj,2LL) + 

               4LL*Power(r,3LL)*Power(xj,3LL))) + 

         exp(2LL*r*xj)*Power(xj,6LL)*

          (-48LL*r*Power(xi,10LL) - 6LL*Power(r,2LL)*Power(xi,11LL) - 

            69LL*Power(xi,7LL)*Power(xj,2LL) + 36LL*r*Power(xi,8LL)*Power(xj,2LL) + 

            8LL*Power(r,2LL)*Power(xi,9LL)*Power(xj,2LL) + 

            84LL*r*Power(xi,6LL)*Power(xj,4LL) - 84LL*r*Power(xi,4LL)*Power(xj,6LL) + 

            9LL*xi*Power(xj,8LL) + 12LL*r*Power(xi,2LL)*Power(xj,8LL) + 

            4LL*Power(r,2LL)*Power(xi,3LL)*Power(xj,8LL) + 

            4LL*Power(xi,9LL)*(-27LL + Power(r,2LL)*Power(xj,2LL)) + 

            Power(xi,3LL)*Power(xj,6LL)*(-63LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            Power(xi,5LL)*(231LL*Power(xj,4LL) - 12LL*Power(r,2LL)*Power(xj,6LL))) + 

         2LL*exp(2LL*r*xj)*Power(xj,7LL)*

          (-24LL*Power(r,2LL)*Power(xi,10LL) - 2LL*Power(r,3LL)*Power(xi,11LL) - 

            69LL*r*Power(xi,7LL)*Power(xj,2LL) + 6LL*Power(xj,8LL) + 9LL*r*xi*Power(xj,8LL) + 

            4LL*r*Power(xi,9LL)*(-27LL + Power(r,2LL)*Power(xj,2LL)) + 

            18LL*Power(xi,8LL)*(-10LL + Power(r,2LL)*Power(xj,2LL)) + 

            6LL*Power(xi,2LL)*Power(xj,6LL)*(-7LL + Power(r,2LL)*Power(xj,2LL)) - 

            42LL*Power(xi,4LL)*Power(xj,4LL)*(-3LL + Power(r,2LL)*Power(xj,2LL)) + 

            r*Power(xi,3LL)*Power(xj,6LL)*(-63LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            6LL*Power(xi,6LL)*Power(xj,2LL)*(-65LL + 7LL*Power(r,2LL)*Power(xj,2LL)) + 

            Power(xi,5LL)*(231LL*r*Power(xj,4LL) - 4LL*Power(r,3LL)*Power(xj,6LL))))/

       (6LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,7LL)*Power(xi + xj,7LL))

    ; }
   
  }
  return S;
}

cl_R DSlater_2S_3S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = 0LL

    ; } else {  S = -(-7430535LL*xi + 8709120LL*exp(2LL*r*xi)*xi - 12303900LL*r*Power(xi,2LL) - 

          9826110LL*Power(r,2LL)*Power(xi,3LL) - 5004720LL*Power(r,3LL)*Power(xi,4LL) - 

          1806840LL*Power(r,4LL)*Power(xi,5LL) - 483840LL*Power(r,5LL)*Power(xi,6LL) - 

          96768LL*Power(r,6LL)*Power(xi,7LL) - 13824LL*Power(r,7LL)*Power(xi,8LL) - 

          1152LL*Power(r,8LL)*Power(xi,9LL))/(4.35456e6*exp(2LL*r*xi)*r) + 

      (-4354560LL + 4354560LL*exp(2LL*r*xi) - 7430535LL*r*xi - 

         6151950LL*Power(r,2LL)*Power(xi,2LL) - 3275370LL*Power(r,3LL)*Power(xi,3LL) - 

         1251180LL*Power(r,4LL)*Power(xi,4LL) - 361368LL*Power(r,5LL)*Power(xi,5LL) - 

         80640LL*Power(r,6LL)*Power(xi,6LL) - 13824LL*Power(r,7LL)*Power(xi,7LL) - 

         1728LL*Power(r,8LL)*Power(xi,8LL) - 128LL*Power(r,9LL)*Power(xi,9LL))/

       (4.35456e6*exp(2LL*r*xi)*Power(r,2LL)) + 

      (xi*(-4354560LL + 4354560LL*exp(2LL*r*xi) - 7430535LL*r*xi - 

           6151950LL*Power(r,2LL)*Power(xi,2LL) - 3275370LL*Power(r,3LL)*Power(xi,3LL) - 

           1251180LL*Power(r,4LL)*Power(xi,4LL) - 361368LL*Power(r,5LL)*Power(xi,5LL) - 

           80640LL*Power(r,6LL)*Power(xi,6LL) - 13824LL*Power(r,7LL)*Power(xi,7LL) - 

           1728LL*Power(r,8LL)*Power(xi,8LL) - 128LL*Power(r,9LL)*Power(xi,9LL)))/

       (2.17728e6*exp(2LL*r*xi)*r)

    ; }
 
  }
  else {
      if (r == 0LL) {  S = 0LL

    ; } else {  S = (90LL*exp(2LL*r*(xi + xj))*Power(Power(xi,2LL) - Power(xj,2LL),9LL) + 

         5LL*exp(2LL*r*xj)*Power(xj,8LL)*

          (-90LL*Power(r,2LL)*Power(xi,12LL) - 6LL*Power(r,3LL)*Power(xi,13LL) + 

            18LL*Power(xj,10LL) + 27LL*r*xi*Power(xj,10LL) + 

            18LL*Power(xi,2LL)*Power(xj,8LL)*(-9LL + Power(r,2LL)*Power(xj,2LL)) - 

            162LL*Power(xi,4LL)*Power(xj,6LL)*(-4LL + Power(r,2LL)*Power(xj,2LL)) - 

            198LL*Power(xi,10LL)*(5LL + Power(r,2LL)*Power(xj,2LL)) - 

            108LL*Power(xi,6LL)*Power(xj,4LL)*(36LL + Power(r,2LL)*Power(xj,2LL)) + 

            2LL*r*Power(xi,5LL)*Power(xj,6LL)*(675LL + Power(r,2LL)*Power(xj,2LL)) - 

            18LL*r*Power(xi,7LL)*Power(xj,4LL)*(-81LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            3LL*r*Power(xi,3LL)*Power(xj,8LL)*(-81LL + 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            r*Power(xi,11LL)*(495LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            9LL*r*Power(xi,9LL)*Power(xj,2LL)*(-233LL + 4LL*Power(r,2LL)*Power(xj,2LL)) + 

            6LL*Power(xi,8LL)*Power(xj,2LL)*(-1063LL + 90LL*Power(r,2LL)*Power(xj,2LL))) - 

         2LL*exp(2LL*r*xi)*Power(xi,6LL)*

          (-90LL*Power(xi,6LL)*Power(xj,6LL)*

             (42LL + 65LL*r*xj + 76LL*Power(r,2LL)*Power(xj,2LL) + 

               22LL*Power(r,3LL)*Power(xj,3LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) - 

            2LL*Power(xj,12LL)*(2970LL + 2475LL*r*xj + 900LL*Power(r,2LL)*Power(xj,2LL) + 

               180LL*Power(r,3LL)*Power(xj,3LL) + 20LL*Power(r,4LL)*Power(xj,4LL) + 

               Power(r,5LL)*Power(xj,5LL)) + 

            10LL*Power(xi,8LL)*Power(xj,4LL)*

             (162LL + 270LL*r*xj + 216LL*Power(r,2LL)*Power(xj,2LL) + 

               122LL*Power(r,3LL)*Power(xj,3LL) + 22LL*Power(r,4LL)*Power(xj,4LL) + 

               Power(r,5LL)*Power(xj,5LL)) - 

            5LL*Power(xi,4LL)*Power(xj,8LL)*

             (-639LL - 3555LL*r*xj - 1452LL*Power(r,2LL)*Power(xj,2LL) - 

               174LL*Power(r,3LL)*Power(xj,3LL) + 6LL*Power(r,4LL)*Power(xj,4LL) + 

               2LL*Power(r,5LL)*Power(xj,5LL)) + 

            Power(xi,12LL)*(45LL + 75LL*r*xj + 60LL*Power(r,2LL)*Power(xj,2LL) + 

               30LL*Power(r,3LL)*Power(xj,3LL) + 10LL*Power(r,4LL)*Power(xj,4LL) + 

               2LL*Power(r,5LL)*Power(xj,5LL)) - 

            Power(xi,10LL)*Power(xj,2LL)*

             (405LL + 675LL*r*xj + 540LL*Power(r,2LL)*Power(xj,2LL) + 

               270LL*Power(r,3LL)*Power(xj,3LL) + 90LL*Power(r,4LL)*Power(xj,4LL) + 

               8LL*Power(r,5LL)*Power(xj,5LL)) + 

            Power(xi,2LL)*Power(xj,10LL)*

             (-21615LL - 9075LL*r*xj - 300LL*Power(r,2LL)*Power(xj,2LL) + 

               490LL*Power(r,3LL)*Power(xj,3LL) + 110LL*Power(r,4LL)*Power(xj,4LL) + 

               8LL*Power(r,5LL)*Power(xj,5LL))))/

       (90LL*exp(2LL*r*(xi + xj))*Power(r,2LL)*Power(xi - xj,9LL)*Power(xi + xj,9LL)) 

    + (90LL*exp(2LL*r*(xi + xj))*Power(Power(xi,2LL) - Power(xj,2LL),9LL) + 

         5LL*exp(2LL*r*xj)*Power(xj,8LL)*

          (-90LL*Power(r,2LL)*Power(xi,12LL) - 6LL*Power(r,3LL)*Power(xi,13LL) + 

            18LL*Power(xj,10LL) + 27LL*r*xi*Power(xj,10LL) + 

            18LL*Power(xi,2LL)*Power(xj,8LL)*(-9LL + Power(r,2LL)*Power(xj,2LL)) - 

            162LL*Power(xi,4LL)*Power(xj,6LL)*(-4LL + Power(r,2LL)*Power(xj,2LL)) - 

            198LL*Power(xi,10LL)*(5LL + Power(r,2LL)*Power(xj,2LL)) - 

            108LL*Power(xi,6LL)*Power(xj,4LL)*(36LL + Power(r,2LL)*Power(xj,2LL)) + 

            2LL*r*Power(xi,5LL)*Power(xj,6LL)*(675LL + Power(r,2LL)*Power(xj,2LL)) - 

            18LL*r*Power(xi,7LL)*Power(xj,4LL)*(-81LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            3LL*r*Power(xi,3LL)*Power(xj,8LL)*(-81LL + 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            r*Power(xi,11LL)*(495LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            9LL*r*Power(xi,9LL)*Power(xj,2LL)*(-233LL + 4LL*Power(r,2LL)*Power(xj,2LL)) + 

            6LL*Power(xi,8LL)*Power(xj,2LL)*(-1063LL + 90LL*Power(r,2LL)*Power(xj,2LL))) - 

         2LL*exp(2LL*r*xi)*Power(xi,6LL)*

          (-90LL*Power(xi,6LL)*Power(xj,6LL)*

             (42LL + 65LL*r*xj + 76LL*Power(r,2LL)*Power(xj,2LL) + 

               22LL*Power(r,3LL)*Power(xj,3LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) - 

            2LL*Power(xj,12LL)*(2970LL + 2475LL*r*xj + 900LL*Power(r,2LL)*Power(xj,2LL) + 

               180LL*Power(r,3LL)*Power(xj,3LL) + 20LL*Power(r,4LL)*Power(xj,4LL) + 

               Power(r,5LL)*Power(xj,5LL)) + 

            10LL*Power(xi,8LL)*Power(xj,4LL)*

             (162LL + 270LL*r*xj + 216LL*Power(r,2LL)*Power(xj,2LL) + 

               122LL*Power(r,3LL)*Power(xj,3LL) + 22LL*Power(r,4LL)*Power(xj,4LL) + 

               Power(r,5LL)*Power(xj,5LL)) - 

            5LL*Power(xi,4LL)*Power(xj,8LL)*

             (-639LL - 3555LL*r*xj - 1452LL*Power(r,2LL)*Power(xj,2LL) - 

               174LL*Power(r,3LL)*Power(xj,3LL) + 6LL*Power(r,4LL)*Power(xj,4LL) + 

               2LL*Power(r,5LL)*Power(xj,5LL)) + 

            Power(xi,12LL)*(45LL + 75LL*r*xj + 60LL*Power(r,2LL)*Power(xj,2LL) + 

               30LL*Power(r,3LL)*Power(xj,3LL) + 10LL*Power(r,4LL)*Power(xj,4LL) + 

               2LL*Power(r,5LL)*Power(xj,5LL)) - 

            Power(xi,10LL)*Power(xj,2LL)*

             (405LL + 675LL*r*xj + 540LL*Power(r,2LL)*Power(xj,2LL) + 

               270LL*Power(r,3LL)*Power(xj,3LL) + 90LL*Power(r,4LL)*Power(xj,4LL) + 

               8LL*Power(r,5LL)*Power(xj,5LL)) + 

            Power(xi,2LL)*Power(xj,10LL)*

             (-21615LL - 9075LL*r*xj - 300LL*Power(r,2LL)*Power(xj,2LL) + 

               490LL*Power(r,3LL)*Power(xj,3LL) + 110LL*Power(r,4LL)*Power(xj,4LL) + 

               8LL*Power(r,5LL)*Power(xj,5LL))))/

       (45LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,9LL)*Power(xi + xj,8LL)) - 

      (180LL*exp(2LL*r*(xi + xj))*(xi + xj)*Power(Power(xi,2LL) - Power(xj,2LL),9LL) + 

         5LL*exp(2LL*r*xj)*Power(xj,8LL)*

          (-180LL*r*Power(xi,12LL) - 18LL*Power(r,2LL)*Power(xi,13LL) - 

            396LL*r*Power(xi,10LL)*Power(xj,2LL) - 

            4LL*Power(r,2LL)*Power(xi,11LL)*Power(xj,2LL) + 

            1080LL*r*Power(xi,8LL)*Power(xj,4LL) + 

            72LL*Power(r,2LL)*Power(xi,9LL)*Power(xj,4LL) - 

            216LL*r*Power(xi,6LL)*Power(xj,6LL) - 

            72LL*Power(r,2LL)*Power(xi,7LL)*Power(xj,6LL) - 

            324LL*r*Power(xi,4LL)*Power(xj,8LL) + 

            4LL*Power(r,2LL)*Power(xi,5LL)*Power(xj,8LL) + 27LL*xi*Power(xj,10LL) + 

            36LL*r*Power(xi,2LL)*Power(xj,10LL) + 

            12LL*Power(r,2LL)*Power(xi,3LL)*Power(xj,10LL) + 

            2LL*Power(xi,5LL)*Power(xj,6LL)*(675LL + Power(r,2LL)*Power(xj,2LL)) - 

            18LL*Power(xi,7LL)*Power(xj,4LL)*(-81LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            3LL*Power(xi,3LL)*Power(xj,8LL)*(-81LL + 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            Power(xi,11LL)*(495LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            9LL*Power(xi,9LL)*Power(xj,2LL)*(-233LL + 4LL*Power(r,2LL)*Power(xj,2LL))) + 

         10LL*exp(2LL*r*xj)*Power(xj,9LL)*

          (-90LL*Power(r,2LL)*Power(xi,12LL) - 6LL*Power(r,3LL)*Power(xi,13LL) + 

            18LL*Power(xj,10LL) + 27LL*r*xi*Power(xj,10LL) + 

            18LL*Power(xi,2LL)*Power(xj,8LL)*(-9LL + Power(r,2LL)*Power(xj,2LL)) - 

            162LL*Power(xi,4LL)*Power(xj,6LL)*(-4LL + Power(r,2LL)*Power(xj,2LL)) - 

            198LL*Power(xi,10LL)*(5LL + Power(r,2LL)*Power(xj,2LL)) - 

            108LL*Power(xi,6LL)*Power(xj,4LL)*(36LL + Power(r,2LL)*Power(xj,2LL)) + 

            2LL*r*Power(xi,5LL)*Power(xj,6LL)*(675LL + Power(r,2LL)*Power(xj,2LL)) - 

            18LL*r*Power(xi,7LL)*Power(xj,4LL)*(-81LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            3LL*r*Power(xi,3LL)*Power(xj,8LL)*(-81LL + 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            r*Power(xi,11LL)*(495LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            9LL*r*Power(xi,9LL)*Power(xj,2LL)*(-233LL + 4LL*Power(r,2LL)*Power(xj,2LL)) + 

            6LL*Power(xi,8LL)*Power(xj,2LL)*(-1063LL + 90LL*Power(r,2LL)*Power(xj,2LL))) - 

         2LL*exp(2LL*r*xi)*Power(xi,6LL)*

          (-90LL*Power(xi,6LL)*Power(xj,6LL)*

             (65LL*xj + 152LL*r*Power(xj,2LL) + 66LL*Power(r,2LL)*Power(xj,3LL) + 

               8LL*Power(r,3LL)*Power(xj,4LL)) - 

            2LL*Power(xj,12LL)*(2475LL*xj + 1800LL*r*Power(xj,2LL) + 

               540LL*Power(r,2LL)*Power(xj,3LL) + 80LL*Power(r,3LL)*Power(xj,4LL) + 

               5LL*Power(r,4LL)*Power(xj,5LL)) + 

            10LL*Power(xi,8LL)*Power(xj,4LL)*

             (270LL*xj + 432LL*r*Power(xj,2LL) + 366LL*Power(r,2LL)*Power(xj,3LL) + 

               88LL*Power(r,3LL)*Power(xj,4LL) + 5LL*Power(r,4LL)*Power(xj,5LL)) - 

            5LL*Power(xi,4LL)*Power(xj,8LL)*

             (-3555LL*xj - 2904LL*r*Power(xj,2LL) - 522LL*Power(r,2LL)*Power(xj,3LL) + 

               24LL*Power(r,3LL)*Power(xj,4LL) + 10LL*Power(r,4LL)*Power(xj,5LL)) + 

            Power(xi,12LL)*(75LL*xj + 120LL*r*Power(xj,2LL) + 

               90LL*Power(r,2LL)*Power(xj,3LL) + 40LL*Power(r,3LL)*Power(xj,4LL) + 

               10LL*Power(r,4LL)*Power(xj,5LL)) - 

            Power(xi,10LL)*Power(xj,2LL)*

             (675LL*xj + 1080LL*r*Power(xj,2LL) + 810LL*Power(r,2LL)*Power(xj,3LL) + 

               360LL*Power(r,3LL)*Power(xj,4LL) + 40LL*Power(r,4LL)*Power(xj,5LL)) + 

            Power(xi,2LL)*Power(xj,10LL)*

             (-9075LL*xj - 600LL*r*Power(xj,2LL) + 1470LL*Power(r,2LL)*Power(xj,3LL) + 

               440LL*Power(r,3LL)*Power(xj,4LL) + 40LL*Power(r,4LL)*Power(xj,5LL))) - 

         4LL*exp(2LL*r*xi)*Power(xi,7LL)*

          (-90LL*Power(xi,6LL)*Power(xj,6LL)*

             (42LL + 65LL*r*xj + 76LL*Power(r,2LL)*Power(xj,2LL) + 

               22LL*Power(r,3LL)*Power(xj,3LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) - 

            2LL*Power(xj,12LL)*(2970LL + 2475LL*r*xj + 900LL*Power(r,2LL)*Power(xj,2LL) + 

               180LL*Power(r,3LL)*Power(xj,3LL) + 20LL*Power(r,4LL)*Power(xj,4LL) + 

               Power(r,5LL)*Power(xj,5LL)) + 

            10LL*Power(xi,8LL)*Power(xj,4LL)*

             (162LL + 270LL*r*xj + 216LL*Power(r,2LL)*Power(xj,2LL) + 

               122LL*Power(r,3LL)*Power(xj,3LL) + 22LL*Power(r,4LL)*Power(xj,4LL) + 

               Power(r,5LL)*Power(xj,5LL)) - 

            5LL*Power(xi,4LL)*Power(xj,8LL)*

             (-639LL - 3555LL*r*xj - 1452LL*Power(r,2LL)*Power(xj,2LL) - 

               174LL*Power(r,3LL)*Power(xj,3LL) + 6LL*Power(r,4LL)*Power(xj,4LL) + 

               2LL*Power(r,5LL)*Power(xj,5LL)) + 

            Power(xi,12LL)*(45LL + 75LL*r*xj + 60LL*Power(r,2LL)*Power(xj,2LL) + 

               30LL*Power(r,3LL)*Power(xj,3LL) + 10LL*Power(r,4LL)*Power(xj,4LL) + 

               2LL*Power(r,5LL)*Power(xj,5LL)) - 

            Power(xi,10LL)*Power(xj,2LL)*

             (405LL + 675LL*r*xj + 540LL*Power(r,2LL)*Power(xj,2LL) + 

               270LL*Power(r,3LL)*Power(xj,3LL) + 90LL*Power(r,4LL)*Power(xj,4LL) + 

               8LL*Power(r,5LL)*Power(xj,5LL)) + 

            Power(xi,2LL)*Power(xj,10LL)*

             (-21615LL - 9075LL*r*xj - 300LL*Power(r,2LL)*Power(xj,2LL) + 

               490LL*Power(r,3LL)*Power(xj,3LL) + 110LL*Power(r,4LL)*Power(xj,4LL) + 

               8LL*Power(r,5LL)*Power(xj,5LL))))/

       (90LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,9LL)*Power(xi + xj,9LL))

    ; }
   
  }
  return S;
}

cl_R DSlater_2S_4S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = 0LL

    ; } else {  S = -(-1125310725LL*xi + 1277337600LL*exp(2LL*r*xi)*xi - 

          1946567700LL*r*Power(xi,2LL) - 1647191700LL*Power(r,2LL)*Power(xi,3LL) - 

          904780800LL*Power(r,3LL)*Power(xi,4LL) - 360498600LL*Power(r,4LL)*Power(xi,5LL) - 

          110103840LL*Power(r,5LL)*Power(xi,6LL) - 26500320LL*Power(r,6LL)*Power(xi,7LL) - 

          5068800LL*Power(r,7LL)*Power(xi,8LL) - 760320LL*Power(r,8LL)*Power(xi,9LL) - 

          84480LL*Power(r,9LL)*Power(xi,10LL) - 5632LL*Power(r,10LL)*Power(xi,11LL))/

       (6.386688e8*exp(2LL*r*xi)*r) + 

      (-638668800LL + 638668800LL*exp(2LL*r*xi) - 1125310725LL*r*xi - 

         973283850LL*Power(r,2LL)*Power(xi,2LL) - 549063900LL*Power(r,3LL)*Power(xi,3LL) - 

         226195200LL*Power(r,4LL)*Power(xi,4LL) - 72099720LL*Power(r,5LL)*Power(xi,5LL) - 

         18350640LL*Power(r,6LL)*Power(xi,6LL) - 3785760LL*Power(r,7LL)*Power(xi,7LL) - 

         633600LL*Power(r,8LL)*Power(xi,8LL) - 84480LL*Power(r,9LL)*Power(xi,9LL) - 

         8448LL*Power(r,10LL)*Power(xi,10LL) - 512LL*Power(r,11LL)*Power(xi,11LL))/

       (6.386688e8*exp(2LL*r*xi)*Power(r,2LL)) + 

      (xi*(-638668800LL + 638668800LL*exp(2LL*r*xi) - 1125310725LL*r*xi - 

           973283850LL*Power(r,2LL)*Power(xi,2LL) - 549063900LL*Power(r,3LL)*Power(xi,3LL) - 

           226195200LL*Power(r,4LL)*Power(xi,4LL) - 72099720LL*Power(r,5LL)*Power(xi,5LL) - 

           18350640LL*Power(r,6LL)*Power(xi,6LL) - 3785760LL*Power(r,7LL)*Power(xi,7LL) - 

           633600LL*Power(r,8LL)*Power(xi,8LL) - 84480LL*Power(r,9LL)*Power(xi,9LL) - 

           8448LL*Power(r,10LL)*Power(xi,10LL) - 512LL*Power(r,11LL)*Power(xi,11LL)))/

       (3.193344e8*exp(2LL*r*xi)*r)

    ; }
 
  }
  else {
      if (r == 0LL) {  S = 0LL

    ; } else {  S = (1260LL*exp(2LL*r*(xi + xj))*Power(Power(xi,2LL) - Power(xj,2LL),11LL) + 

         210LL*exp(2LL*r*xj)*Power(xj,10LL)*

          (-36LL*Power(r,2LL)*Power(xi,14LL) - 2LL*Power(r,3LL)*Power(xi,15LL) - 

            1287LL*r*Power(xi,9LL)*Power(xj,4LL) + 6LL*Power(xj,12LL) + 

            9LL*r*xi*Power(xj,12LL) - 22LL*r*Power(xi,7LL)*Power(xj,6LL)*

             (-135LL + Power(r,2LL)*Power(xj,2LL)) + 

            6LL*Power(xi,2LL)*Power(xj,10LL)*(-11LL + Power(r,2LL)*Power(xj,2LL)) - 

            66LL*Power(xi,4LL)*Power(xj,8LL)*(-5LL + Power(r,2LL)*Power(xj,2LL)) + 

            8LL*r*Power(xi,5LL)*Power(xj,8LL)*(99LL + Power(r,2LL)*Power(xj,2LL)) + 

            r*Power(xi,3LL)*Power(xj,10LL)*(-99LL + 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            132LL*Power(xi,6LL)*Power(xj,6LL)*(27LL + 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            78LL*Power(xi,12LL)*(7LL + 3LL*Power(r,2LL)*Power(xj,2LL)) - 

            2LL*r*Power(xi,13LL)*(117LL + 4LL*Power(r,2LL)*Power(xj,2LL)) + 

            66LL*Power(xi,8LL)*Power(xj,4LL)*(-191LL + 6LL*Power(r,2LL)*Power(xj,2LL)) + 

            r*Power(xi,11LL)*Power(xj,2LL)*(-2151LL + 22LL*Power(r,2LL)*Power(xj,2LL)) + 

            6LL*Power(xi,10LL)*Power(xj,2LL)*(-1099LL + 33LL*Power(r,2LL)*Power(xj,2LL))) + 

         exp(2LL*r*xi)*Power(xi,6LL)*

          (-385LL*Power(xi,8LL)*Power(xj,8LL)*

             (1080LL + 1935LL*r*xj + 1350LL*Power(r,2LL)*Power(xj,2LL) + 

               1170LL*Power(r,3LL)*Power(xj,3LL) + 420LL*Power(r,4LL)*Power(xj,4LL) + 

               66LL*Power(r,5LL)*Power(xj,5LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) + 

            7LL*Power(xi,6LL)*Power(xj,10LL)*

             (99540LL + 58095LL*r*xj + 190710LL*Power(r,2LL)*Power(xj,2LL) + 

               100950LL*Power(r,3LL)*Power(xj,3LL) + 21660LL*Power(r,4LL)*Power(xj,4LL) + 

               1938LL*Power(r,5LL)*Power(xj,5LL) + 4LL*Power(r,6LL)*Power(xj,6LL) - 

               8LL*Power(r,7LL)*Power(xj,7LL)) + 

            4LL*Power(xj,16LL)*(135135LL + 135135LL*r*xj + 

               62370LL*Power(r,2LL)*Power(xj,2LL) + 17325LL*Power(r,3LL)*Power(xj,3LL) + 

               3150LL*Power(r,4LL)*Power(xj,4LL) + 378LL*Power(r,5LL)*Power(xj,5LL) + 

               28LL*Power(r,6LL)*Power(xj,6LL) + Power(r,7LL)*Power(xj,7LL)) - 

            Power(xi,16LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r,2LL)*Power(xj,2LL) + 

               1050LL*Power(r,3LL)*Power(xj,3LL) + 420LL*Power(r,4LL)*Power(xj,4LL) + 

               126LL*Power(r,5LL)*Power(xj,5LL) + 28LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL)) + 

            7LL*Power(xi,4LL)*Power(xj,12LL)*

             (114660LL - 343395LL*r*xj - 242910LL*Power(r,2LL)*Power(xj,2LL) - 

               61950LL*Power(r,3LL)*Power(xj,3LL) - 6060LL*Power(r,4LL)*Power(xj,4LL) + 

               282LL*Power(r,5LL)*Power(xj,5LL) + 116LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) - 

            7LL*Power(xi,12LL)*Power(xj,4LL)*

             (9900LL + 17325LL*r*xj + 14850LL*Power(r,2LL)*Power(xj,2LL) + 

               8250LL*Power(r,3LL)*Power(xj,3LL) + 3300LL*Power(r,4LL)*Power(xj,4LL) + 

               1074LL*Power(r,5LL)*Power(xj,5LL) + 164LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) + 

            7LL*Power(xi,10LL)*Power(xj,6LL)*

             (29700LL + 51975LL*r*xj + 44550LL*Power(r,2LL)*Power(xj,2LL) + 

               23850LL*Power(r,3LL)*Power(xj,3LL) + 11700LL*Power(r,4LL)*Power(xj,4LL) + 

               2814LL*Power(r,5LL)*Power(xj,5LL) + 284LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) + 

            Power(xi,14LL)*Power(xj,2LL)*

             (13860LL + 24255LL*r*xj + 20790LL*Power(r,2LL)*Power(xj,2LL) + 

               11550LL*Power(r,3LL)*Power(xj,3LL) + 4620LL*Power(r,4LL)*Power(xj,4LL) + 

               1386LL*Power(r,5LL)*Power(xj,5LL) + 308LL*Power(r,6LL)*Power(xj,6LL) + 

               24LL*Power(r,7LL)*Power(xj,7LL)) - 

            Power(xi,2LL)*Power(xj,14LL)*

             (-3063060LL - 1936935LL*r*xj - 408870LL*Power(r,2LL)*Power(xj,2LL) + 

               11550LL*Power(r,3LL)*Power(xj,3LL) + 23100LL*Power(r,4LL)*Power(xj,4LL) + 

               5082LL*Power(r,5LL)*Power(xj,5LL) + 532LL*Power(r,6LL)*Power(xj,6LL) + 

               24LL*Power(r,7LL)*Power(xj,7LL))))/

       (1260LL*exp(2LL*r*(xi + xj))*Power(r,2LL)*Power(xi - xj,11LL)*

         Power(xi + xj,11LL)) + (1260LL*exp(2LL*r*(xi + xj))*

          Power(Power(xi,2LL) - Power(xj,2LL),11LL) + 

         210LL*exp(2LL*r*xj)*Power(xj,10LL)*

          (-36LL*Power(r,2LL)*Power(xi,14LL) - 2LL*Power(r,3LL)*Power(xi,15LL) - 

            1287LL*r*Power(xi,9LL)*Power(xj,4LL) + 6LL*Power(xj,12LL) + 

            9LL*r*xi*Power(xj,12LL) - 22LL*r*Power(xi,7LL)*Power(xj,6LL)*

             (-135LL + Power(r,2LL)*Power(xj,2LL)) + 

            6LL*Power(xi,2LL)*Power(xj,10LL)*(-11LL + Power(r,2LL)*Power(xj,2LL)) - 

            66LL*Power(xi,4LL)*Power(xj,8LL)*(-5LL + Power(r,2LL)*Power(xj,2LL)) + 

            8LL*r*Power(xi,5LL)*Power(xj,8LL)*(99LL + Power(r,2LL)*Power(xj,2LL)) + 

            r*Power(xi,3LL)*Power(xj,10LL)*(-99LL + 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            132LL*Power(xi,6LL)*Power(xj,6LL)*(27LL + 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            78LL*Power(xi,12LL)*(7LL + 3LL*Power(r,2LL)*Power(xj,2LL)) - 

            2LL*r*Power(xi,13LL)*(117LL + 4LL*Power(r,2LL)*Power(xj,2LL)) + 

            66LL*Power(xi,8LL)*Power(xj,4LL)*(-191LL + 6LL*Power(r,2LL)*Power(xj,2LL)) + 

            r*Power(xi,11LL)*Power(xj,2LL)*(-2151LL + 22LL*Power(r,2LL)*Power(xj,2LL)) + 

            6LL*Power(xi,10LL)*Power(xj,2LL)*(-1099LL + 33LL*Power(r,2LL)*Power(xj,2LL))) + 

         exp(2LL*r*xi)*Power(xi,6LL)*

          (-385LL*Power(xi,8LL)*Power(xj,8LL)*

             (1080LL + 1935LL*r*xj + 1350LL*Power(r,2LL)*Power(xj,2LL) + 

               1170LL*Power(r,3LL)*Power(xj,3LL) + 420LL*Power(r,4LL)*Power(xj,4LL) + 

               66LL*Power(r,5LL)*Power(xj,5LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) + 

            7LL*Power(xi,6LL)*Power(xj,10LL)*

             (99540LL + 58095LL*r*xj + 190710LL*Power(r,2LL)*Power(xj,2LL) + 

               100950LL*Power(r,3LL)*Power(xj,3LL) + 21660LL*Power(r,4LL)*Power(xj,4LL) + 

               1938LL*Power(r,5LL)*Power(xj,5LL) + 4LL*Power(r,6LL)*Power(xj,6LL) - 

               8LL*Power(r,7LL)*Power(xj,7LL)) + 

            4LL*Power(xj,16LL)*(135135LL + 135135LL*r*xj + 

               62370LL*Power(r,2LL)*Power(xj,2LL) + 17325LL*Power(r,3LL)*Power(xj,3LL) + 

               3150LL*Power(r,4LL)*Power(xj,4LL) + 378LL*Power(r,5LL)*Power(xj,5LL) + 

               28LL*Power(r,6LL)*Power(xj,6LL) + Power(r,7LL)*Power(xj,7LL)) - 

            Power(xi,16LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r,2LL)*Power(xj,2LL) + 

               1050LL*Power(r,3LL)*Power(xj,3LL) + 420LL*Power(r,4LL)*Power(xj,4LL) + 

               126LL*Power(r,5LL)*Power(xj,5LL) + 28LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL)) + 

            7LL*Power(xi,4LL)*Power(xj,12LL)*

             (114660LL - 343395LL*r*xj - 242910LL*Power(r,2LL)*Power(xj,2LL) - 

               61950LL*Power(r,3LL)*Power(xj,3LL) - 6060LL*Power(r,4LL)*Power(xj,4LL) + 

               282LL*Power(r,5LL)*Power(xj,5LL) + 116LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) - 

            7LL*Power(xi,12LL)*Power(xj,4LL)*

             (9900LL + 17325LL*r*xj + 14850LL*Power(r,2LL)*Power(xj,2LL) + 

               8250LL*Power(r,3LL)*Power(xj,3LL) + 3300LL*Power(r,4LL)*Power(xj,4LL) + 

               1074LL*Power(r,5LL)*Power(xj,5LL) + 164LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) + 

            7LL*Power(xi,10LL)*Power(xj,6LL)*

             (29700LL + 51975LL*r*xj + 44550LL*Power(r,2LL)*Power(xj,2LL) + 

               23850LL*Power(r,3LL)*Power(xj,3LL) + 11700LL*Power(r,4LL)*Power(xj,4LL) + 

               2814LL*Power(r,5LL)*Power(xj,5LL) + 284LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) + 

            Power(xi,14LL)*Power(xj,2LL)*

             (13860LL + 24255LL*r*xj + 20790LL*Power(r,2LL)*Power(xj,2LL) + 

               11550LL*Power(r,3LL)*Power(xj,3LL) + 4620LL*Power(r,4LL)*Power(xj,4LL) + 

               1386LL*Power(r,5LL)*Power(xj,5LL) + 308LL*Power(r,6LL)*Power(xj,6LL) + 

               24LL*Power(r,7LL)*Power(xj,7LL)) - 

            Power(xi,2LL)*Power(xj,14LL)*

             (-3063060LL - 1936935LL*r*xj - 408870LL*Power(r,2LL)*Power(xj,2LL) + 

               11550LL*Power(r,3LL)*Power(xj,3LL) + 23100LL*Power(r,4LL)*Power(xj,4LL) + 

               5082LL*Power(r,5LL)*Power(xj,5LL) + 532LL*Power(r,6LL)*Power(xj,6LL) + 

               24LL*Power(r,7LL)*Power(xj,7LL))))/

       (630LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,11LL)*Power(xi + xj,10LL)) - 

      (2520LL*exp(2LL*r*(xi + xj))*(xi + xj)*

          Power(Power(xi,2LL) - Power(xj,2LL),11LL) + 

         210LL*exp(2LL*r*xj)*Power(xj,10LL)*

          (-72LL*r*Power(xi,14LL) - 6LL*Power(r,2LL)*Power(xi,15LL) - 

            468LL*r*Power(xi,12LL)*Power(xj,2LL) - 

            16LL*Power(r,2LL)*Power(xi,13LL)*Power(xj,2LL) - 

            1287LL*Power(xi,9LL)*Power(xj,4LL) + 396LL*r*Power(xi,10LL)*Power(xj,4LL) + 

            44LL*Power(r,2LL)*Power(xi,11LL)*Power(xj,4LL) + 

            792LL*r*Power(xi,8LL)*Power(xj,6LL) - 528LL*r*Power(xi,6LL)*Power(xj,8LL) - 

            44LL*Power(r,2LL)*Power(xi,7LL)*Power(xj,8LL) - 

            132LL*r*Power(xi,4LL)*Power(xj,10LL) + 

            16LL*Power(r,2LL)*Power(xi,5LL)*Power(xj,10LL) + 9LL*xi*Power(xj,12LL) + 

            12LL*r*Power(xi,2LL)*Power(xj,12LL) + 

            4LL*Power(r,2LL)*Power(xi,3LL)*Power(xj,12LL) - 

            22LL*Power(xi,7LL)*Power(xj,6LL)*(-135LL + Power(r,2LL)*Power(xj,2LL)) + 

            8LL*Power(xi,5LL)*Power(xj,8LL)*(99LL + Power(r,2LL)*Power(xj,2LL)) + 

            Power(xi,3LL)*Power(xj,10LL)*(-99LL + 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            2LL*Power(xi,13LL)*(117LL + 4LL*Power(r,2LL)*Power(xj,2LL)) + 

            Power(xi,11LL)*Power(xj,2LL)*(-2151LL + 22LL*Power(r,2LL)*Power(xj,2LL))) + 

         420LL*exp(2LL*r*xj)*Power(xj,11LL)*

          (-36LL*Power(r,2LL)*Power(xi,14LL) - 2LL*Power(r,3LL)*Power(xi,15LL) - 

            1287LL*r*Power(xi,9LL)*Power(xj,4LL) + 6LL*Power(xj,12LL) + 

            9LL*r*xi*Power(xj,12LL) - 22LL*r*Power(xi,7LL)*Power(xj,6LL)*

             (-135LL + Power(r,2LL)*Power(xj,2LL)) + 

            6LL*Power(xi,2LL)*Power(xj,10LL)*(-11LL + Power(r,2LL)*Power(xj,2LL)) - 

            66LL*Power(xi,4LL)*Power(xj,8LL)*(-5LL + Power(r,2LL)*Power(xj,2LL)) + 

            8LL*r*Power(xi,5LL)*Power(xj,8LL)*(99LL + Power(r,2LL)*Power(xj,2LL)) + 

            r*Power(xi,3LL)*Power(xj,10LL)*(-99LL + 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            132LL*Power(xi,6LL)*Power(xj,6LL)*(27LL + 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            78LL*Power(xi,12LL)*(7LL + 3LL*Power(r,2LL)*Power(xj,2LL)) - 

            2LL*r*Power(xi,13LL)*(117LL + 4LL*Power(r,2LL)*Power(xj,2LL)) + 

            66LL*Power(xi,8LL)*Power(xj,4LL)*(-191LL + 6LL*Power(r,2LL)*Power(xj,2LL)) + 

            r*Power(xi,11LL)*Power(xj,2LL)*(-2151LL + 22LL*Power(r,2LL)*Power(xj,2LL)) + 

            6LL*Power(xi,10LL)*Power(xj,2LL)*(-1099LL + 33LL*Power(r,2LL)*Power(xj,2LL))) + 

         exp(2LL*r*xi)*Power(xi,6LL)*

          (-385LL*Power(xi,8LL)*Power(xj,8LL)*

             (1935LL*xj + 2700LL*r*Power(xj,2LL) + 3510LL*Power(r,2LL)*Power(xj,3LL) + 

               1680LL*Power(r,3LL)*Power(xj,4LL) + 330LL*Power(r,4LL)*Power(xj,5LL) + 

               24LL*Power(r,5LL)*Power(xj,6LL)) + 

            7LL*Power(xi,6LL)*Power(xj,10LL)*

             (58095LL*xj + 381420LL*r*Power(xj,2LL) + 

               302850LL*Power(r,2LL)*Power(xj,3LL) + 86640LL*Power(r,3LL)*Power(xj,4LL) + 

               9690LL*Power(r,4LL)*Power(xj,5LL) + 24LL*Power(r,5LL)*Power(xj,6LL) - 

               56LL*Power(r,6LL)*Power(xj,7LL)) + 

            4LL*Power(xj,16LL)*(135135LL*xj + 124740LL*r*Power(xj,2LL) + 

               51975LL*Power(r,2LL)*Power(xj,3LL) + 12600LL*Power(r,3LL)*Power(xj,4LL) + 

               1890LL*Power(r,4LL)*Power(xj,5LL) + 168LL*Power(r,5LL)*Power(xj,6LL) + 

               7LL*Power(r,6LL)*Power(xj,7LL)) - 

            Power(xi,16LL)*(2205LL*xj + 3780LL*r*Power(xj,2LL) + 

               3150LL*Power(r,2LL)*Power(xj,3LL) + 1680LL*Power(r,3LL)*Power(xj,4LL) + 

               630LL*Power(r,4LL)*Power(xj,5LL) + 168LL*Power(r,5LL)*Power(xj,6LL) + 

               28LL*Power(r,6LL)*Power(xj,7LL)) + 

            7LL*Power(xi,4LL)*Power(xj,12LL)*

             (-343395LL*xj - 485820LL*r*Power(xj,2LL) - 

               185850LL*Power(r,2LL)*Power(xj,3LL) - 24240LL*Power(r,3LL)*Power(xj,4LL) + 

               1410LL*Power(r,4LL)*Power(xj,5LL) + 696LL*Power(r,5LL)*Power(xj,6LL) + 

               56LL*Power(r,6LL)*Power(xj,7LL)) - 

            7LL*Power(xi,12LL)*Power(xj,4LL)*

             (17325LL*xj + 29700LL*r*Power(xj,2LL) + 24750LL*Power(r,2LL)*Power(xj,3LL) + 

               13200LL*Power(r,3LL)*Power(xj,4LL) + 5370LL*Power(r,4LL)*Power(xj,5LL) + 

               984LL*Power(r,5LL)*Power(xj,6LL) + 56LL*Power(r,6LL)*Power(xj,7LL)) + 

            7LL*Power(xi,10LL)*Power(xj,6LL)*

             (51975LL*xj + 89100LL*r*Power(xj,2LL) + 71550LL*Power(r,2LL)*Power(xj,3LL) + 

               46800LL*Power(r,3LL)*Power(xj,4LL) + 14070LL*Power(r,4LL)*Power(xj,5LL) + 

               1704LL*Power(r,5LL)*Power(xj,6LL) + 56LL*Power(r,6LL)*Power(xj,7LL)) + 

            Power(xi,14LL)*Power(xj,2LL)*

             (24255LL*xj + 41580LL*r*Power(xj,2LL) + 34650LL*Power(r,2LL)*Power(xj,3LL) + 

               18480LL*Power(r,3LL)*Power(xj,4LL) + 6930LL*Power(r,4LL)*Power(xj,5LL) + 

               1848LL*Power(r,5LL)*Power(xj,6LL) + 168LL*Power(r,6LL)*Power(xj,7LL)) - 

            Power(xi,2LL)*Power(xj,14LL)*

             (-1936935LL*xj - 817740LL*r*Power(xj,2LL) + 

               34650LL*Power(r,2LL)*Power(xj,3LL) + 92400LL*Power(r,3LL)*Power(xj,4LL) + 

               25410LL*Power(r,4LL)*Power(xj,5LL) + 3192LL*Power(r,5LL)*Power(xj,6LL) + 

               168LL*Power(r,6LL)*Power(xj,7LL))) + 

         2LL*exp(2LL*r*xi)*Power(xi,7LL)*

          (-385LL*Power(xi,8LL)*Power(xj,8LL)*

             (1080LL + 1935LL*r*xj + 1350LL*Power(r,2LL)*Power(xj,2LL) + 

               1170LL*Power(r,3LL)*Power(xj,3LL) + 420LL*Power(r,4LL)*Power(xj,4LL) + 

               66LL*Power(r,5LL)*Power(xj,5LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) + 

            7LL*Power(xi,6LL)*Power(xj,10LL)*

             (99540LL + 58095LL*r*xj + 190710LL*Power(r,2LL)*Power(xj,2LL) + 

               100950LL*Power(r,3LL)*Power(xj,3LL) + 21660LL*Power(r,4LL)*Power(xj,4LL) + 

               1938LL*Power(r,5LL)*Power(xj,5LL) + 4LL*Power(r,6LL)*Power(xj,6LL) - 

               8LL*Power(r,7LL)*Power(xj,7LL)) + 

            4LL*Power(xj,16LL)*(135135LL + 135135LL*r*xj + 

               62370LL*Power(r,2LL)*Power(xj,2LL) + 17325LL*Power(r,3LL)*Power(xj,3LL) + 

               3150LL*Power(r,4LL)*Power(xj,4LL) + 378LL*Power(r,5LL)*Power(xj,5LL) + 

               28LL*Power(r,6LL)*Power(xj,6LL) + Power(r,7LL)*Power(xj,7LL)) - 

            Power(xi,16LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r,2LL)*Power(xj,2LL) + 

               1050LL*Power(r,3LL)*Power(xj,3LL) + 420LL*Power(r,4LL)*Power(xj,4LL) + 

               126LL*Power(r,5LL)*Power(xj,5LL) + 28LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL)) + 

            7LL*Power(xi,4LL)*Power(xj,12LL)*

             (114660LL - 343395LL*r*xj - 242910LL*Power(r,2LL)*Power(xj,2LL) - 

               61950LL*Power(r,3LL)*Power(xj,3LL) - 6060LL*Power(r,4LL)*Power(xj,4LL) + 

               282LL*Power(r,5LL)*Power(xj,5LL) + 116LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) - 

            7LL*Power(xi,12LL)*Power(xj,4LL)*

             (9900LL + 17325LL*r*xj + 14850LL*Power(r,2LL)*Power(xj,2LL) + 

               8250LL*Power(r,3LL)*Power(xj,3LL) + 3300LL*Power(r,4LL)*Power(xj,4LL) + 

               1074LL*Power(r,5LL)*Power(xj,5LL) + 164LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) + 

            7LL*Power(xi,10LL)*Power(xj,6LL)*

             (29700LL + 51975LL*r*xj + 44550LL*Power(r,2LL)*Power(xj,2LL) + 

               23850LL*Power(r,3LL)*Power(xj,3LL) + 11700LL*Power(r,4LL)*Power(xj,4LL) + 

               2814LL*Power(r,5LL)*Power(xj,5LL) + 284LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) + 

            Power(xi,14LL)*Power(xj,2LL)*

             (13860LL + 24255LL*r*xj + 20790LL*Power(r,2LL)*Power(xj,2LL) + 

               11550LL*Power(r,3LL)*Power(xj,3LL) + 4620LL*Power(r,4LL)*Power(xj,4LL) + 

               1386LL*Power(r,5LL)*Power(xj,5LL) + 308LL*Power(r,6LL)*Power(xj,6LL) + 

               24LL*Power(r,7LL)*Power(xj,7LL)) - 

            Power(xi,2LL)*Power(xj,14LL)*

             (-3063060LL - 1936935LL*r*xj - 408870LL*Power(r,2LL)*Power(xj,2LL) + 

               11550LL*Power(r,3LL)*Power(xj,3LL) + 23100LL*Power(r,4LL)*Power(xj,4LL) + 

               5082LL*Power(r,5LL)*Power(xj,5LL) + 532LL*Power(r,6LL)*Power(xj,6LL) + 

               24LL*Power(r,7LL)*Power(xj,7LL))))/

       (1260LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,11LL)*Power(xi + xj,11LL))

    ; }
   
  }
  return S;
}

cl_R DSlater_2S_5S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = 0LL

    ; } else {  S = -(-224622748350LL*xi + 249080832000LL*exp(2LL*r*xi)*xi - 

          400329329400LL*r*Power(xi,2LL) - 351747621225LL*Power(r,2LL)*Power(xi,3LL) - 

          202556554200LL*Power(r,3LL)*Power(xi,4LL) - 

          85662076500LL*Power(r,4LL)*Power(xi,5LL) - 

          28229160960LL*Power(r,5LL)*Power(xi,6LL) - 

          7498370880LL*Power(r,6LL)*Power(xi,7LL) - 

          1635828480LL*Power(r,7LL)*Power(xi,8LL) - 

          295289280LL*Power(r,8LL)*Power(xi,9LL) - 43929600LL*Power(r,9LL)*Power(xi,10LL) - 

          5271552LL*Power(r,10LL)*Power(xi,11LL) - 479232LL*Power(r,11LL)*Power(xi,12LL) - 

          26624LL*Power(r,12LL)*Power(xi,13LL))/(1.24540416e11*exp(2LL*r*xi)*r) + 

      (-124540416000LL + 124540416000LL*exp(2LL*r*xi) - 224622748350LL*r*xi - 

         200164664700LL*Power(r,2LL)*Power(xi,2LL) - 

         117249207075LL*Power(r,3LL)*Power(xi,3LL) - 

         50639138550LL*Power(r,4LL)*Power(xi,4LL) - 

         17132415300LL*Power(r,5LL)*Power(xi,5LL) - 

         4704860160LL*Power(r,6LL)*Power(xi,6LL) - 

         1071195840LL*Power(r,7LL)*Power(xi,7LL) - 204478560LL*Power(r,8LL)*Power(xi,8LL) - 

         32809920LL*Power(r,9LL)*Power(xi,9LL) - 4392960LL*Power(r,10LL)*Power(xi,10LL) - 

         479232LL*Power(r,11LL)*Power(xi,11LL) - 39936LL*Power(r,12LL)*Power(xi,12LL) - 

         2048LL*Power(r,13LL)*Power(xi,13LL))/

       (1.24540416e11*exp(2LL*r*xi)*Power(r,2LL)) + 

      (xi*(-124540416000LL + 124540416000LL*exp(2LL*r*xi) - 224622748350LL*r*xi - 

           200164664700LL*Power(r,2LL)*Power(xi,2LL) - 

           117249207075LL*Power(r,3LL)*Power(xi,3LL) - 

           50639138550LL*Power(r,4LL)*Power(xi,4LL) - 

           17132415300LL*Power(r,5LL)*Power(xi,5LL) - 

           4704860160LL*Power(r,6LL)*Power(xi,6LL) - 

           1071195840LL*Power(r,7LL)*Power(xi,7LL) - 

           204478560LL*Power(r,8LL)*Power(xi,8LL) - 32809920LL*Power(r,9LL)*Power(xi,9LL) - 

           4392960LL*Power(r,10LL)*Power(xi,10LL) - 479232LL*Power(r,11LL)*Power(xi,11LL) - 

           39936LL*Power(r,12LL)*Power(xi,12LL) - 2048LL*Power(r,13LL)*Power(xi,13LL)))/

       (6.2270208e10*exp(2LL*r*xi)*r)

    ; }
 
  }
  else {
      if (r == 0LL) {  S = 0LL

    ; } else {  S = (28350LL*exp(2LL*r*(xi + xj))*Power(Power(xi,2LL) - Power(xj,2LL),13LL) + 

         945LL*exp(2LL*r*xj)*Power(xj,12LL)*

          (-210LL*Power(r,2LL)*Power(xi,16LL) - 10LL*Power(r,3LL)*Power(xi,17LL) + 

            30LL*Power(xj,14LL) + 45LL*r*xi*Power(xj,14LL) + 

            39LL*r*Power(xi,7LL)*Power(xj,8LL)*(1309LL - 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            858LL*Power(xi,8LL)*Power(xj,6LL)*(-305LL + Power(r,2LL)*Power(xj,2LL)) + 

            30LL*Power(xi,2LL)*Power(xj,12LL)*(-13LL + Power(r,2LL)*Power(xj,2LL)) - 

            390LL*Power(xi,4LL)*Power(xj,10LL)*(-6LL + Power(r,2LL)*Power(xj,2LL)) - 

            143LL*r*Power(xi,9LL)*Power(xj,6LL)*(-153LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            5LL*r*Power(xi,3LL)*Power(xj,12LL)*(-117LL + 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            45LL*r*Power(xi,15LL)*(35LL + 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            138LL*Power(xi,12LL)*Power(xj,2LL)*(580LL + 13LL*Power(r,2LL)*Power(xj,2LL)) - 

            150LL*Power(xi,14LL)*(28LL + 17LL*Power(r,2LL)*Power(xj,2LL)) + 

            13LL*r*Power(xi,11LL)*Power(xj,4LL)*

             (-4071LL + 22LL*Power(r,2LL)*Power(xj,2LL)) + 

            3LL*r*Power(xi,13LL)*Power(xj,2LL)*(-8135LL + 26LL*Power(r,2LL)*Power(xj,2LL)) + 

            3LL*r*Power(xi,5LL)*Power(xj,10LL)*(2171LL + 30LL*Power(r,2LL)*Power(xj,2LL)) + 

            234LL*Power(xi,10LL)*Power(xj,4LL)*(-1235LL + 33LL*Power(r,2LL)*Power(xj,2LL)) - 

            78LL*Power(xi,6LL)*Power(xj,8LL)*(550LL + 47LL*Power(r,2LL)*Power(xj,2LL))) - 

         2LL*exp(2LL*r*xi)*Power(xi,6LL)*

          (-819LL*Power(xi,10LL)*Power(xj,10LL)*

             (22275LL + 39780LL*r*xj + 38160LL*Power(r,2LL)*Power(xj,2LL) + 

               16560LL*Power(r,3LL)*Power(xj,3LL) + 9840LL*Power(r,4LL)*Power(xj,4LL) + 

               3900LL*Power(r,5LL)*Power(xj,5LL) + 816LL*Power(r,6LL)*Power(xj,6LL) + 

               88LL*Power(r,7LL)*Power(xj,7LL) + 4LL*Power(r,8LL)*Power(xj,8LL)) + 

            Power(xi,20LL)*(14175LL + 25515LL*r*xj + 22680LL*Power(r,2LL)*Power(xj,2LL) + 

               13230LL*Power(r,3LL)*Power(xj,3LL) + 5670LL*Power(r,4LL)*Power(xj,4LL) + 

               1890LL*Power(r,5LL)*Power(xj,5LL) + 504LL*Power(r,6LL)*Power(xj,6LL) + 

               108LL*Power(r,7LL)*Power(xj,7LL) + 18LL*Power(r,8LL)*Power(xj,8LL) + 

               2LL*Power(r,9LL)*Power(xj,9LL)) - 

            Power(xj,20LL)*(16216200LL + 18243225LL*r*xj + 

               9729720LL*Power(r,2LL)*Power(xj,2LL) + 

               3243240LL*Power(r,3LL)*Power(xj,3LL) + 

               748440LL*Power(r,4LL)*Power(xj,4LL) + 124740LL*Power(r,5LL)*Power(xj,5LL) + 

               15120LL*Power(r,6LL)*Power(xj,6LL) + 1296LL*Power(r,7LL)*Power(xj,7LL) + 

               72LL*Power(r,8LL)*Power(xj,8LL) + 2LL*Power(r,9LL)*Power(xj,9LL)) + 

            18LL*Power(xi,16LL)*Power(xj,4LL)*

             (61425LL + 110565LL*r*xj + 98280LL*Power(r,2LL)*Power(xj,2LL) + 

               57330LL*Power(r,3LL)*Power(xj,3LL) + 24570LL*Power(r,4LL)*Power(xj,4LL) + 

               8190LL*Power(r,5LL)*Power(xj,5LL) + 2184LL*Power(r,6LL)*Power(xj,6LL) + 

               496LL*Power(r,7LL)*Power(xj,7LL) + 64LL*Power(r,8LL)*Power(xj,8LL) + 

               3LL*Power(r,9LL)*Power(xj,9LL)) - 

            18LL*Power(xi,4LL)*Power(xj,16LL)*

             (6572475LL - 3161340LL*r*xj - 4782960LL*Power(r,2LL)*Power(xj,2LL) - 

               1912365LL*Power(r,3LL)*Power(xj,3LL) - 

               378105LL*Power(r,4LL)*Power(xj,4LL) - 34125LL*Power(r,5LL)*Power(xj,5LL) + 

               1092LL*Power(r,6LL)*Power(xj,6LL) + 650LL*Power(r,7LL)*Power(xj,7LL) + 

               71LL*Power(r,8LL)*Power(xj,8LL) + 3LL*Power(r,9LL)*Power(xj,9LL)) - 

            21LL*Power(xi,8LL)*Power(xj,12LL)*

             (-1063800LL - 2775735LL*r*xj - 862920LL*Power(r,2LL)*Power(xj,2LL) - 

               1132020LL*Power(r,3LL)*Power(xj,3LL) - 

               698580LL*Power(r,4LL)*Power(xj,4LL) - 196920LL*Power(r,5LL)*Power(xj,5LL) - 

               28992LL*Power(r,6LL)*Power(xj,6LL) - 2064LL*Power(r,7LL)*Power(xj,7LL) - 

               24LL*Power(r,8LL)*Power(xj,8LL) + 4LL*Power(r,9LL)*Power(xj,9LL)) + 

            21LL*Power(xi,12LL)*Power(xj,8LL)*

             (482625LL + 868725LL*r*xj + 772200LL*Power(r,2LL)*Power(xj,2LL) + 

               455400LL*Power(r,3LL)*Power(xj,3LL) + 178200LL*Power(r,4LL)*Power(xj,4LL) + 

               72180LL*Power(r,5LL)*Power(xj,5LL) + 19920LL*Power(r,6LL)*Power(xj,6LL) + 

               2952LL*Power(r,7LL)*Power(xj,7LL) + 204LL*Power(r,8LL)*Power(xj,8LL) + 

               4LL*Power(r,9LL)*Power(xj,9LL)) + 

            6LL*Power(xi,6LL)*Power(xj,14LL)*

             (-10357200LL + 5071815LL*r*xj - 6463800LL*Power(r,2LL)*Power(xj,2LL) - 

               7151130LL*Power(r,3LL)*Power(xj,3LL) - 

               2572290LL*Power(r,4LL)*Power(xj,4LL) - 

               468720LL*Power(r,5LL)*Power(xj,5LL) - 42672LL*Power(r,6LL)*Power(xj,6LL) - 

               648LL*Power(r,7LL)*Power(xj,7LL) + 228LL*Power(r,8LL)*Power(xj,8LL) + 

               16LL*Power(r,9LL)*Power(xj,9LL)) - 

            Power(xi,18LL)*Power(xj,2LL)*

             (184275LL + 331695LL*r*xj + 294840LL*Power(r,2LL)*Power(xj,2LL) + 

               171990LL*Power(r,3LL)*Power(xj,3LL) + 73710LL*Power(r,4LL)*Power(xj,4LL) + 

               24570LL*Power(r,5LL)*Power(xj,5LL) + 6552LL*Power(r,6LL)*Power(xj,6LL) + 

               1404LL*Power(r,7LL)*Power(xj,7LL) + 234LL*Power(r,8LL)*Power(xj,8LL) + 

               16LL*Power(r,9LL)*Power(xj,9LL)) + 

            Power(xi,2LL)*Power(xj,18LL)*

             (-133783650LL - 107432325LL*r*xj - 35675640LL*Power(r,2LL)*Power(xj,2LL) - 

               5135130LL*Power(r,3LL)*Power(xj,3LL) + 

               270270LL*Power(r,4LL)*Power(xj,4LL) + 270270LL*Power(r,5LL)*Power(xj,5LL) + 

               57960LL*Power(r,6LL)*Power(xj,6LL) + 6948LL*Power(r,7LL)*Power(xj,7LL) + 

               486LL*Power(r,8LL)*Power(xj,8LL) + 16LL*Power(r,9LL)*Power(xj,9LL)) - 

            6LL*Power(xi,14LL)*Power(xj,6LL)*

             (675675LL + 1216215LL*r*xj + 1081080LL*Power(r,2LL)*Power(xj,2LL) + 

               630630LL*Power(r,3LL)*Power(xj,3LL) + 270270LL*Power(r,4LL)*Power(xj,4LL) + 

               88200LL*Power(r,5LL)*Power(xj,5LL) + 26544LL*Power(r,6LL)*Power(xj,6LL) + 

               5160LL*Power(r,7LL)*Power(xj,7LL) + 492LL*Power(r,8LL)*Power(xj,8LL) + 

               16LL*Power(r,9LL)*Power(xj,9LL))))/

       (28350LL*exp(2LL*r*(xi + xj))*Power(r,2LL)*Power(xi - xj,13LL)*

         Power(xi + xj,13LL)) + (28350LL*exp(2LL*r*(xi + xj))*

          Power(Power(xi,2LL) - Power(xj,2LL),13LL) + 

         945LL*exp(2LL*r*xj)*Power(xj,12LL)*

          (-210LL*Power(r,2LL)*Power(xi,16LL) - 10LL*Power(r,3LL)*Power(xi,17LL) + 

            30LL*Power(xj,14LL) + 45LL*r*xi*Power(xj,14LL) + 

            39LL*r*Power(xi,7LL)*Power(xj,8LL)*(1309LL - 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            858LL*Power(xi,8LL)*Power(xj,6LL)*(-305LL + Power(r,2LL)*Power(xj,2LL)) + 

            30LL*Power(xi,2LL)*Power(xj,12LL)*(-13LL + Power(r,2LL)*Power(xj,2LL)) - 

            390LL*Power(xi,4LL)*Power(xj,10LL)*(-6LL + Power(r,2LL)*Power(xj,2LL)) - 

            143LL*r*Power(xi,9LL)*Power(xj,6LL)*(-153LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            5LL*r*Power(xi,3LL)*Power(xj,12LL)*(-117LL + 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            45LL*r*Power(xi,15LL)*(35LL + 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            138LL*Power(xi,12LL)*Power(xj,2LL)*(580LL + 13LL*Power(r,2LL)*Power(xj,2LL)) - 

            150LL*Power(xi,14LL)*(28LL + 17LL*Power(r,2LL)*Power(xj,2LL)) + 

            13LL*r*Power(xi,11LL)*Power(xj,4LL)*

             (-4071LL + 22LL*Power(r,2LL)*Power(xj,2LL)) + 

            3LL*r*Power(xi,13LL)*Power(xj,2LL)*(-8135LL + 26LL*Power(r,2LL)*Power(xj,2LL)) + 

            3LL*r*Power(xi,5LL)*Power(xj,10LL)*(2171LL + 30LL*Power(r,2LL)*Power(xj,2LL)) + 

            234LL*Power(xi,10LL)*Power(xj,4LL)*(-1235LL + 33LL*Power(r,2LL)*Power(xj,2LL)) - 

            78LL*Power(xi,6LL)*Power(xj,8LL)*(550LL + 47LL*Power(r,2LL)*Power(xj,2LL))) - 

         2LL*exp(2LL*r*xi)*Power(xi,6LL)*

          (-819LL*Power(xi,10LL)*Power(xj,10LL)*

             (22275LL + 39780LL*r*xj + 38160LL*Power(r,2LL)*Power(xj,2LL) + 

               16560LL*Power(r,3LL)*Power(xj,3LL) + 9840LL*Power(r,4LL)*Power(xj,4LL) + 

               3900LL*Power(r,5LL)*Power(xj,5LL) + 816LL*Power(r,6LL)*Power(xj,6LL) + 

               88LL*Power(r,7LL)*Power(xj,7LL) + 4LL*Power(r,8LL)*Power(xj,8LL)) + 

            Power(xi,20LL)*(14175LL + 25515LL*r*xj + 22680LL*Power(r,2LL)*Power(xj,2LL) + 

               13230LL*Power(r,3LL)*Power(xj,3LL) + 5670LL*Power(r,4LL)*Power(xj,4LL) + 

               1890LL*Power(r,5LL)*Power(xj,5LL) + 504LL*Power(r,6LL)*Power(xj,6LL) + 

               108LL*Power(r,7LL)*Power(xj,7LL) + 18LL*Power(r,8LL)*Power(xj,8LL) + 

               2LL*Power(r,9LL)*Power(xj,9LL)) - 

            Power(xj,20LL)*(16216200LL + 18243225LL*r*xj + 

               9729720LL*Power(r,2LL)*Power(xj,2LL) + 

               3243240LL*Power(r,3LL)*Power(xj,3LL) + 

               748440LL*Power(r,4LL)*Power(xj,4LL) + 124740LL*Power(r,5LL)*Power(xj,5LL) + 

               15120LL*Power(r,6LL)*Power(xj,6LL) + 1296LL*Power(r,7LL)*Power(xj,7LL) + 

               72LL*Power(r,8LL)*Power(xj,8LL) + 2LL*Power(r,9LL)*Power(xj,9LL)) + 

            18LL*Power(xi,16LL)*Power(xj,4LL)*

             (61425LL + 110565LL*r*xj + 98280LL*Power(r,2LL)*Power(xj,2LL) + 

               57330LL*Power(r,3LL)*Power(xj,3LL) + 24570LL*Power(r,4LL)*Power(xj,4LL) + 

               8190LL*Power(r,5LL)*Power(xj,5LL) + 2184LL*Power(r,6LL)*Power(xj,6LL) + 

               496LL*Power(r,7LL)*Power(xj,7LL) + 64LL*Power(r,8LL)*Power(xj,8LL) + 

               3LL*Power(r,9LL)*Power(xj,9LL)) - 

            18LL*Power(xi,4LL)*Power(xj,16LL)*

             (6572475LL - 3161340LL*r*xj - 4782960LL*Power(r,2LL)*Power(xj,2LL) - 

               1912365LL*Power(r,3LL)*Power(xj,3LL) - 

               378105LL*Power(r,4LL)*Power(xj,4LL) - 34125LL*Power(r,5LL)*Power(xj,5LL) + 

               1092LL*Power(r,6LL)*Power(xj,6LL) + 650LL*Power(r,7LL)*Power(xj,7LL) + 

               71LL*Power(r,8LL)*Power(xj,8LL) + 3LL*Power(r,9LL)*Power(xj,9LL)) - 

            21LL*Power(xi,8LL)*Power(xj,12LL)*

             (-1063800LL - 2775735LL*r*xj - 862920LL*Power(r,2LL)*Power(xj,2LL) - 

               1132020LL*Power(r,3LL)*Power(xj,3LL) - 

               698580LL*Power(r,4LL)*Power(xj,4LL) - 196920LL*Power(r,5LL)*Power(xj,5LL) - 

               28992LL*Power(r,6LL)*Power(xj,6LL) - 2064LL*Power(r,7LL)*Power(xj,7LL) - 

               24LL*Power(r,8LL)*Power(xj,8LL) + 4LL*Power(r,9LL)*Power(xj,9LL)) + 

            21LL*Power(xi,12LL)*Power(xj,8LL)*

             (482625LL + 868725LL*r*xj + 772200LL*Power(r,2LL)*Power(xj,2LL) + 

               455400LL*Power(r,3LL)*Power(xj,3LL) + 178200LL*Power(r,4LL)*Power(xj,4LL) + 

               72180LL*Power(r,5LL)*Power(xj,5LL) + 19920LL*Power(r,6LL)*Power(xj,6LL) + 

               2952LL*Power(r,7LL)*Power(xj,7LL) + 204LL*Power(r,8LL)*Power(xj,8LL) + 

               4LL*Power(r,9LL)*Power(xj,9LL)) + 

            6LL*Power(xi,6LL)*Power(xj,14LL)*

             (-10357200LL + 5071815LL*r*xj - 6463800LL*Power(r,2LL)*Power(xj,2LL) - 

               7151130LL*Power(r,3LL)*Power(xj,3LL) - 

               2572290LL*Power(r,4LL)*Power(xj,4LL) - 

               468720LL*Power(r,5LL)*Power(xj,5LL) - 42672LL*Power(r,6LL)*Power(xj,6LL) - 

               648LL*Power(r,7LL)*Power(xj,7LL) + 228LL*Power(r,8LL)*Power(xj,8LL) + 

               16LL*Power(r,9LL)*Power(xj,9LL)) - 

            Power(xi,18LL)*Power(xj,2LL)*

             (184275LL + 331695LL*r*xj + 294840LL*Power(r,2LL)*Power(xj,2LL) + 

               171990LL*Power(r,3LL)*Power(xj,3LL) + 73710LL*Power(r,4LL)*Power(xj,4LL) + 

               24570LL*Power(r,5LL)*Power(xj,5LL) + 6552LL*Power(r,6LL)*Power(xj,6LL) + 

               1404LL*Power(r,7LL)*Power(xj,7LL) + 234LL*Power(r,8LL)*Power(xj,8LL) + 

               16LL*Power(r,9LL)*Power(xj,9LL)) + 

            Power(xi,2LL)*Power(xj,18LL)*

             (-133783650LL - 107432325LL*r*xj - 35675640LL*Power(r,2LL)*Power(xj,2LL) - 

               5135130LL*Power(r,3LL)*Power(xj,3LL) + 

               270270LL*Power(r,4LL)*Power(xj,4LL) + 270270LL*Power(r,5LL)*Power(xj,5LL) + 

               57960LL*Power(r,6LL)*Power(xj,6LL) + 6948LL*Power(r,7LL)*Power(xj,7LL) + 

               486LL*Power(r,8LL)*Power(xj,8LL) + 16LL*Power(r,9LL)*Power(xj,9LL)) - 

            6LL*Power(xi,14LL)*Power(xj,6LL)*

             (675675LL + 1216215LL*r*xj + 1081080LL*Power(r,2LL)*Power(xj,2LL) + 

               630630LL*Power(r,3LL)*Power(xj,3LL) + 270270LL*Power(r,4LL)*Power(xj,4LL) + 

               88200LL*Power(r,5LL)*Power(xj,5LL) + 26544LL*Power(r,6LL)*Power(xj,6LL) + 

               5160LL*Power(r,7LL)*Power(xj,7LL) + 492LL*Power(r,8LL)*Power(xj,8LL) + 

               16LL*Power(r,9LL)*Power(xj,9LL))))/

       (14175LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,13LL)*Power(xi + xj,12LL)) - 

      (56700LL*exp(2LL*r*(xi + xj))*(xi + xj)*

          Power(Power(xi,2LL) - Power(xj,2LL),13LL) + 

         945LL*exp(2LL*r*xj)*Power(xj,12LL)*

          (-420LL*r*Power(xi,16LL) - 30LL*Power(r,2LL)*Power(xi,17LL) - 

            5100LL*r*Power(xi,14LL)*Power(xj,2LL) - 

            180LL*Power(r,2LL)*Power(xi,15LL)*Power(xj,2LL) - 

            3588LL*r*Power(xi,12LL)*Power(xj,4LL) + 

            156LL*Power(r,2LL)*Power(xi,13LL)*Power(xj,4LL) + 

            15444LL*r*Power(xi,10LL)*Power(xj,6LL) + 

            572LL*Power(r,2LL)*Power(xi,11LL)*Power(xj,6LL) + 

            1716LL*r*Power(xi,8LL)*Power(xj,8LL) - 

            572LL*Power(r,2LL)*Power(xi,9LL)*Power(xj,8LL) - 

            7332LL*r*Power(xi,6LL)*Power(xj,10LL) - 

            156LL*Power(r,2LL)*Power(xi,7LL)*Power(xj,10LL) - 

            780LL*r*Power(xi,4LL)*Power(xj,12LL) + 

            180LL*Power(r,2LL)*Power(xi,5LL)*Power(xj,12LL) + 45LL*xi*Power(xj,14LL) + 

            60LL*r*Power(xi,2LL)*Power(xj,14LL) + 

            20LL*Power(r,2LL)*Power(xi,3LL)*Power(xj,14LL) + 

            39LL*Power(xi,7LL)*Power(xj,8LL)*(1309LL - 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            143LL*Power(xi,9LL)*Power(xj,6LL)*(-153LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            5LL*Power(xi,3LL)*Power(xj,12LL)*(-117LL + 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            45LL*Power(xi,15LL)*(35LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            13LL*Power(xi,11LL)*Power(xj,4LL)*(-4071LL + 22LL*Power(r,2LL)*Power(xj,2LL)) + 

            3LL*Power(xi,13LL)*Power(xj,2LL)*(-8135LL + 26LL*Power(r,2LL)*Power(xj,2LL)) + 

            3LL*Power(xi,5LL)*Power(xj,10LL)*(2171LL + 30LL*Power(r,2LL)*Power(xj,2LL))) + 

         1890LL*exp(2LL*r*xj)*Power(xj,13LL)*

          (-210LL*Power(r,2LL)*Power(xi,16LL) - 10LL*Power(r,3LL)*Power(xi,17LL) + 

            30LL*Power(xj,14LL) + 45LL*r*xi*Power(xj,14LL) + 

            39LL*r*Power(xi,7LL)*Power(xj,8LL)*(1309LL - 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            858LL*Power(xi,8LL)*Power(xj,6LL)*(-305LL + Power(r,2LL)*Power(xj,2LL)) + 

            30LL*Power(xi,2LL)*Power(xj,12LL)*(-13LL + Power(r,2LL)*Power(xj,2LL)) - 

            390LL*Power(xi,4LL)*Power(xj,10LL)*(-6LL + Power(r,2LL)*Power(xj,2LL)) - 

            143LL*r*Power(xi,9LL)*Power(xj,6LL)*(-153LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            5LL*r*Power(xi,3LL)*Power(xj,12LL)*(-117LL + 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            45LL*r*Power(xi,15LL)*(35LL + 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            138LL*Power(xi,12LL)*Power(xj,2LL)*(580LL + 13LL*Power(r,2LL)*Power(xj,2LL)) - 

            150LL*Power(xi,14LL)*(28LL + 17LL*Power(r,2LL)*Power(xj,2LL)) + 

            13LL*r*Power(xi,11LL)*Power(xj,4LL)*(-4071LL + 22LL*Power(r,2LL)*Power(xj,2LL)) + 

            3LL*r*Power(xi,13LL)*Power(xj,2LL)*(-8135LL + 26LL*Power(r,2LL)*Power(xj,2LL)) + 

            3LL*r*Power(xi,5LL)*Power(xj,10LL)*(2171LL + 30LL*Power(r,2LL)*Power(xj,2LL)) + 

            234LL*Power(xi,10LL)*Power(xj,4LL)*(-1235LL + 33LL*Power(r,2LL)*Power(xj,2LL)) - 

            78LL*Power(xi,6LL)*Power(xj,8LL)*(550LL + 47LL*Power(r,2LL)*Power(xj,2LL))) - 

         2LL*exp(2LL*r*xi)*Power(xi,6LL)*

          (-819LL*Power(xi,10LL)*Power(xj,10LL)*

             (39780LL*xj + 76320LL*r*Power(xj,2LL) + 49680LL*Power(r,2LL)*Power(xj,3LL) + 

               39360LL*Power(r,3LL)*Power(xj,4LL) + 19500LL*Power(r,4LL)*Power(xj,5LL) + 

               4896LL*Power(r,5LL)*Power(xj,6LL) + 616LL*Power(r,6LL)*Power(xj,7LL) + 

               32LL*Power(r,7LL)*Power(xj,8LL)) + 

            Power(xi,20LL)*(25515LL*xj + 45360LL*r*Power(xj,2LL) + 

               39690LL*Power(r,2LL)*Power(xj,3LL) + 22680LL*Power(r,3LL)*Power(xj,4LL) + 

               9450LL*Power(r,4LL)*Power(xj,5LL) + 3024LL*Power(r,5LL)*Power(xj,6LL) + 

               756LL*Power(r,6LL)*Power(xj,7LL) + 144LL*Power(r,7LL)*Power(xj,8LL) + 

               18LL*Power(r,8LL)*Power(xj,9LL)) - 

            Power(xj,20LL)*(18243225LL*xj + 19459440LL*r*Power(xj,2LL) + 

               9729720LL*Power(r,2LL)*Power(xj,3LL) + 

               2993760LL*Power(r,3LL)*Power(xj,4LL) + 

               623700LL*Power(r,4LL)*Power(xj,5LL) + 90720LL*Power(r,5LL)*Power(xj,6LL) + 

               9072LL*Power(r,6LL)*Power(xj,7LL) + 576LL*Power(r,7LL)*Power(xj,8LL) + 

               18LL*Power(r,8LL)*Power(xj,9LL)) + 

            18LL*Power(xi,16LL)*Power(xj,4LL)*

             (110565LL*xj + 196560LL*r*Power(xj,2LL) + 

               171990LL*Power(r,2LL)*Power(xj,3LL) + 98280LL*Power(r,3LL)*Power(xj,4LL) + 

               40950LL*Power(r,4LL)*Power(xj,5LL) + 13104LL*Power(r,5LL)*Power(xj,6LL) + 

               3472LL*Power(r,6LL)*Power(xj,7LL) + 512LL*Power(r,7LL)*Power(xj,8LL) + 

               27LL*Power(r,8LL)*Power(xj,9LL)) - 

            18LL*Power(xi,4LL)*Power(xj,16LL)*

             (-3161340LL*xj - 9565920LL*r*Power(xj,2LL) - 

               5737095LL*Power(r,2LL)*Power(xj,3LL) - 

               1512420LL*Power(r,3LL)*Power(xj,4LL) - 

               170625LL*Power(r,4LL)*Power(xj,5LL) + 6552LL*Power(r,5LL)*Power(xj,6LL) + 

               4550LL*Power(r,6LL)*Power(xj,7LL) + 568LL*Power(r,7LL)*Power(xj,8LL) + 

               27LL*Power(r,8LL)*Power(xj,9LL)) - 

            21LL*Power(xi,8LL)*Power(xj,12LL)*

             (-2775735LL*xj - 1725840LL*r*Power(xj,2LL) - 

               3396060LL*Power(r,2LL)*Power(xj,3LL) - 

               2794320LL*Power(r,3LL)*Power(xj,4LL) - 

               984600LL*Power(r,4LL)*Power(xj,5LL) - 173952LL*Power(r,5LL)*Power(xj,6LL) - 

               14448LL*Power(r,6LL)*Power(xj,7LL) - 192LL*Power(r,7LL)*Power(xj,8LL) + 

               36LL*Power(r,8LL)*Power(xj,9LL)) + 

            21LL*Power(xi,12LL)*Power(xj,8LL)*

             (868725LL*xj + 1544400LL*r*Power(xj,2LL) + 

               1366200LL*Power(r,2LL)*Power(xj,3LL) + 

               712800LL*Power(r,3LL)*Power(xj,4LL) + 360900LL*Power(r,4LL)*Power(xj,5LL) + 

               119520LL*Power(r,5LL)*Power(xj,6LL) + 20664LL*Power(r,6LL)*Power(xj,7LL) + 

               1632LL*Power(r,7LL)*Power(xj,8LL) + 36LL*Power(r,8LL)*Power(xj,9LL)) + 

            6LL*Power(xi,6LL)*Power(xj,14LL)*

             (5071815LL*xj - 12927600LL*r*Power(xj,2LL) - 

               21453390LL*Power(r,2LL)*Power(xj,3LL) - 

               10289160LL*Power(r,3LL)*Power(xj,4LL) - 

               2343600LL*Power(r,4LL)*Power(xj,5LL) - 

               256032LL*Power(r,5LL)*Power(xj,6LL) - 4536LL*Power(r,6LL)*Power(xj,7LL) + 

               1824LL*Power(r,7LL)*Power(xj,8LL) + 144LL*Power(r,8LL)*Power(xj,9LL)) - 

            Power(xi,18LL)*Power(xj,2LL)*

             (331695LL*xj + 589680LL*r*Power(xj,2LL) + 

               515970LL*Power(r,2LL)*Power(xj,3LL) + 294840LL*Power(r,3LL)*Power(xj,4LL) + 

               122850LL*Power(r,4LL)*Power(xj,5LL) + 39312LL*Power(r,5LL)*Power(xj,6LL) + 

               9828LL*Power(r,6LL)*Power(xj,7LL) + 1872LL*Power(r,7LL)*Power(xj,8LL) + 

               144LL*Power(r,8LL)*Power(xj,9LL)) + 

            Power(xi,2LL)*Power(xj,18LL)*

             (-107432325LL*xj - 71351280LL*r*Power(xj,2LL) - 

               15405390LL*Power(r,2LL)*Power(xj,3LL) + 

               1081080LL*Power(r,3LL)*Power(xj,4LL) + 

               1351350LL*Power(r,4LL)*Power(xj,5LL) + 

               347760LL*Power(r,5LL)*Power(xj,6LL) + 48636LL*Power(r,6LL)*Power(xj,7LL) + 

               3888LL*Power(r,7LL)*Power(xj,8LL) + 144LL*Power(r,8LL)*Power(xj,9LL)) - 

            6LL*Power(xi,14LL)*Power(xj,6LL)*

             (1216215LL*xj + 2162160LL*r*Power(xj,2LL) + 

               1891890LL*Power(r,2LL)*Power(xj,3LL) + 

               1081080LL*Power(r,3LL)*Power(xj,4LL) + 441000LL*Power(r,4LL)*Power(xj,5LL) + 

               159264LL*Power(r,5LL)*Power(xj,6LL) + 36120LL*Power(r,6LL)*Power(xj,7LL) + 

               3936LL*Power(r,7LL)*Power(xj,8LL) + 144LL*Power(r,8LL)*Power(xj,9LL))) - 

         4LL*exp(2LL*r*xi)*Power(xi,7LL)*

          (-819LL*Power(xi,10LL)*Power(xj,10LL)*

             (22275LL + 39780LL*r*xj + 38160LL*Power(r,2LL)*Power(xj,2LL) + 

               16560LL*Power(r,3LL)*Power(xj,3LL) + 9840LL*Power(r,4LL)*Power(xj,4LL) + 

               3900LL*Power(r,5LL)*Power(xj,5LL) + 816LL*Power(r,6LL)*Power(xj,6LL) + 

               88LL*Power(r,7LL)*Power(xj,7LL) + 4LL*Power(r,8LL)*Power(xj,8LL)) + 

            Power(xi,20LL)*(14175LL + 25515LL*r*xj + 22680LL*Power(r,2LL)*Power(xj,2LL) + 

               13230LL*Power(r,3LL)*Power(xj,3LL) + 5670LL*Power(r,4LL)*Power(xj,4LL) + 

               1890LL*Power(r,5LL)*Power(xj,5LL) + 504LL*Power(r,6LL)*Power(xj,6LL) + 

               108LL*Power(r,7LL)*Power(xj,7LL) + 18LL*Power(r,8LL)*Power(xj,8LL) + 

               2LL*Power(r,9LL)*Power(xj,9LL)) - 

            Power(xj,20LL)*(16216200LL + 18243225LL*r*xj + 

               9729720LL*Power(r,2LL)*Power(xj,2LL) + 

               3243240LL*Power(r,3LL)*Power(xj,3LL) + 748440LL*Power(r,4LL)*Power(xj,4LL) + 

               124740LL*Power(r,5LL)*Power(xj,5LL) + 15120LL*Power(r,6LL)*Power(xj,6LL) + 

               1296LL*Power(r,7LL)*Power(xj,7LL) + 72LL*Power(r,8LL)*Power(xj,8LL) + 

               2LL*Power(r,9LL)*Power(xj,9LL)) + 

            18LL*Power(xi,16LL)*Power(xj,4LL)*

             (61425LL + 110565LL*r*xj + 98280LL*Power(r,2LL)*Power(xj,2LL) + 

               57330LL*Power(r,3LL)*Power(xj,3LL) + 24570LL*Power(r,4LL)*Power(xj,4LL) + 

               8190LL*Power(r,5LL)*Power(xj,5LL) + 2184LL*Power(r,6LL)*Power(xj,6LL) + 

               496LL*Power(r,7LL)*Power(xj,7LL) + 64LL*Power(r,8LL)*Power(xj,8LL) + 

               3LL*Power(r,9LL)*Power(xj,9LL)) - 

            18LL*Power(xi,4LL)*Power(xj,16LL)*

             (6572475LL - 3161340LL*r*xj - 4782960LL*Power(r,2LL)*Power(xj,2LL) - 

               1912365LL*Power(r,3LL)*Power(xj,3LL) - 378105LL*Power(r,4LL)*Power(xj,4LL) - 

               34125LL*Power(r,5LL)*Power(xj,5LL) + 1092LL*Power(r,6LL)*Power(xj,6LL) + 

               650LL*Power(r,7LL)*Power(xj,7LL) + 71LL*Power(r,8LL)*Power(xj,8LL) + 

               3LL*Power(r,9LL)*Power(xj,9LL)) - 

            21LL*Power(xi,8LL)*Power(xj,12LL)*

             (-1063800LL - 2775735LL*r*xj - 862920LL*Power(r,2LL)*Power(xj,2LL) - 

               1132020LL*Power(r,3LL)*Power(xj,3LL) - 698580LL*Power(r,4LL)*Power(xj,4LL) - 

               196920LL*Power(r,5LL)*Power(xj,5LL) - 28992LL*Power(r,6LL)*Power(xj,6LL) - 

               2064LL*Power(r,7LL)*Power(xj,7LL) - 24LL*Power(r,8LL)*Power(xj,8LL) + 

               4LL*Power(r,9LL)*Power(xj,9LL)) + 

            21LL*Power(xi,12LL)*Power(xj,8LL)*

             (482625LL + 868725LL*r*xj + 772200LL*Power(r,2LL)*Power(xj,2LL) + 

               455400LL*Power(r,3LL)*Power(xj,3LL) + 178200LL*Power(r,4LL)*Power(xj,4LL) + 

               72180LL*Power(r,5LL)*Power(xj,5LL) + 19920LL*Power(r,6LL)*Power(xj,6LL) + 

               2952LL*Power(r,7LL)*Power(xj,7LL) + 204LL*Power(r,8LL)*Power(xj,8LL) + 

               4LL*Power(r,9LL)*Power(xj,9LL)) + 

            6LL*Power(xi,6LL)*Power(xj,14LL)*

             (-10357200LL + 5071815LL*r*xj - 6463800LL*Power(r,2LL)*Power(xj,2LL) - 

               7151130LL*Power(r,3LL)*Power(xj,3LL) - 

               2572290LL*Power(r,4LL)*Power(xj,4LL) - 468720LL*Power(r,5LL)*Power(xj,5LL) - 

               42672LL*Power(r,6LL)*Power(xj,6LL) - 648LL*Power(r,7LL)*Power(xj,7LL) + 

               228LL*Power(r,8LL)*Power(xj,8LL) + 16LL*Power(r,9LL)*Power(xj,9LL)) - 

            Power(xi,18LL)*Power(xj,2LL)*

             (184275LL + 331695LL*r*xj + 294840LL*Power(r,2LL)*Power(xj,2LL) + 

               171990LL*Power(r,3LL)*Power(xj,3LL) + 73710LL*Power(r,4LL)*Power(xj,4LL) + 

               24570LL*Power(r,5LL)*Power(xj,5LL) + 6552LL*Power(r,6LL)*Power(xj,6LL) + 

               1404LL*Power(r,7LL)*Power(xj,7LL) + 234LL*Power(r,8LL)*Power(xj,8LL) + 

               16LL*Power(r,9LL)*Power(xj,9LL)) + 

            Power(xi,2LL)*Power(xj,18LL)*

             (-133783650LL - 107432325LL*r*xj - 35675640LL*Power(r,2LL)*Power(xj,2LL) - 

               5135130LL*Power(r,3LL)*Power(xj,3LL) + 270270LL*Power(r,4LL)*Power(xj,4LL) + 

               270270LL*Power(r,5LL)*Power(xj,5LL) + 57960LL*Power(r,6LL)*Power(xj,6LL) + 

               6948LL*Power(r,7LL)*Power(xj,7LL) + 486LL*Power(r,8LL)*Power(xj,8LL) + 

               16LL*Power(r,9LL)*Power(xj,9LL)) - 

            6LL*Power(xi,14LL)*Power(xj,6LL)*

             (675675LL + 1216215LL*r*xj + 1081080LL*Power(r,2LL)*Power(xj,2LL) + 

               630630LL*Power(r,3LL)*Power(xj,3LL) + 270270LL*Power(r,4LL)*Power(xj,4LL) + 

               88200LL*Power(r,5LL)*Power(xj,5LL) + 26544LL*Power(r,6LL)*Power(xj,6LL) + 

               5160LL*Power(r,7LL)*Power(xj,7LL) + 492LL*Power(r,8LL)*Power(xj,8LL) + 

               16LL*Power(r,9LL)*Power(xj,9LL))))/

       (28350LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,13LL)*Power(xi + xj,13LL))

    ; }
   
  }
  return S;
}

cl_R DSlater_2S_1S(cl_R r,cl_R xi,cl_R xj)
{
  return DSlater_1S_2S(r,xj,xi);
}

cl_R DSlater_3S_3S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = 0LL

    ; } else {  S = -(-2503064025LL*xi + 2874009600LL*exp(2LL*r*xi)*xi - 

          4264236900LL*r*Power(xi,2LL) - 3541992300LL*Power(r,2LL)*Power(xi,3LL) - 

          1906027200LL*Power(r,3LL)*Power(xi,4LL) - 

          744282000LL*Power(r,4LL)*Power(xi,5LL) - 223534080LL*Power(r,5LL)*Power(xi,6LL) - 

          53222400LL*Power(r,6LL)*Power(xi,7LL) - 10137600LL*Power(r,7LL)*Power(xi,8LL) - 

          1520640LL*Power(r,8LL)*Power(xi,9LL) - 168960LL*Power(r,9LL)*Power(xi,10LL) - 

          11264LL*Power(r,10LL)*Power(xi,11LL))/(1.4370048e9*exp(2LL*r*xi)*r) + 

      (-1437004800LL + 1437004800LL*exp(2LL*r*xi) - 2503064025LL*r*xi - 

         2132118450LL*Power(r,2LL)*Power(xi,2LL) - 

         1180664100LL*Power(r,3LL)*Power(xi,3LL) - 476506800LL*Power(r,4LL)*Power(xi,4LL) - 

         148856400LL*Power(r,5LL)*Power(xi,5LL) - 37255680LL*Power(r,6LL)*Power(xi,6LL) - 

         7603200LL*Power(r,7LL)*Power(xi,7LL) - 1267200LL*Power(r,8LL)*Power(xi,8LL) - 

         168960LL*Power(r,9LL)*Power(xi,9LL) - 16896LL*Power(r,10LL)*Power(xi,10LL) - 

         1024LL*Power(r,11LL)*Power(xi,11LL))/(1.4370048e9*exp(2LL*r*xi)*Power(r,2LL)) 

    + (xi*(-1437004800LL + 1437004800LL*exp(2LL*r*xi) - 2503064025LL*r*xi - 

           2132118450LL*Power(r,2LL)*Power(xi,2LL) - 

           1180664100LL*Power(r,3LL)*Power(xi,3LL) - 

           476506800LL*Power(r,4LL)*Power(xi,4LL) - 148856400LL*Power(r,5LL)*Power(xi,5LL) - 

           37255680LL*Power(r,6LL)*Power(xi,6LL) - 7603200LL*Power(r,7LL)*Power(xi,7LL) - 

           1267200LL*Power(r,8LL)*Power(xi,8LL) - 168960LL*Power(r,9LL)*Power(xi,9LL) - 

           16896LL*Power(r,10LL)*Power(xi,10LL) - 1024LL*Power(r,11LL)*Power(xi,11LL)))/

       (7.185024e8*exp(2LL*r*xi)*r)

    ; }
 
  }
  else {
      if (r == 0LL) {  S = 0LL

    ; } else {  S = (135LL*exp(2LL*r*(xi + xj))*Power(Power(xi,2LL) - Power(xj,2LL),11LL) + 

         exp(2LL*r*xj)*Power(xj,8LL)*

          (-150LL*Power(r,4LL)*Power(xi,18LL) - 6LL*Power(r,5LL)*Power(xi,19LL) + 

            135LL*Power(xj,14LL) + 225LL*r*xi*Power(xj,14LL) + 

            10LL*Power(r,3LL)*Power(xi,17LL)*(-165LL + Power(r,2LL)*Power(xj,2LL)) - 

            30LL*Power(r,2LL)*Power(xi,16LL)*(330LL + Power(r,2LL)*Power(xj,2LL)) + 

            45LL*r*Power(xi,3LL)*Power(xj,12LL)*(-55LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            45LL*Power(xi,2LL)*Power(xj,12LL)*(-33LL + 4LL*Power(r,2LL)*Power(xj,2LL)) + 

            r*Power(xi,9LL)*Power(xj,6LL)*

             (234135LL - 4950LL*Power(r,2LL)*Power(xj,2LL) - 

               34LL*Power(r,4LL)*Power(xj,4LL)) - 

            5LL*r*Power(xi,7LL)*Power(xj,8LL)*

             (6237LL - 1242LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            3LL*r*Power(xi,5LL)*Power(xj,10LL)*

             (4125LL - 330LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            15LL*Power(xi,4LL)*Power(xj,10LL)*

             (495LL - 132LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) - 

            165LL*Power(xi,6LL)*Power(xj,8LL)*

             (135LL - 60LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) - 

            5LL*r*Power(xi,13LL)*Power(xj,2LL)*

             (43875LL - 3438LL*Power(r,2LL)*Power(xj,2LL) + 22LL*Power(r,4LL)*Power(xj,4LL)) 

    + 5LL*r*Power(xi,11LL)*Power(xj,4LL)*

             (7695LL - 2442LL*Power(r,2LL)*Power(xj,2LL) + 22LL*Power(r,4LL)*Power(xj,4LL)) 

    + 15LL*Power(xi,8LL)*Power(xj,6LL)*(-33LL - 3564LL*Power(r,2LL)*Power(xj,2LL) + 

               26LL*Power(r,4LL)*Power(xj,4LL)) + 

            r*Power(xi,15LL)*(-32175LL - 3690LL*Power(r,2LL)*Power(xj,2LL) + 

               34LL*Power(r,4LL)*Power(xj,4LL)) + 

            15LL*Power(xi,10LL)*Power(xj,4LL)*

             (-32277LL + 1364LL*Power(r,2LL)*Power(xj,2LL) + 

               66LL*Power(r,4LL)*Power(xj,4LL)) + 

            15LL*Power(xi,14LL)*(-3003LL - 2932LL*Power(r,2LL)*Power(xj,2LL) + 

               94LL*Power(r,4LL)*Power(xj,4LL)) - 

            15LL*Power(xi,12LL)*Power(xj,2LL)*

             (28119LL - 5252LL*Power(r,2LL)*Power(xj,2LL) + 154LL*Power(r,4LL)*Power(xj,4LL))

    ) + exp(2LL*r*xi)*Power(xi,8LL)*

          (-5LL*Power(xi,2LL)*Power(xj,12LL)*

             (-84357LL - 43875LL*r*xj - 8796LL*Power(r,2LL)*Power(xj,2LL) - 

               738LL*Power(r,3LL)*Power(xj,3LL) - 6LL*Power(r,4LL)*Power(xj,4LL) + 

               2LL*Power(r,5LL)*Power(xj,5LL)) - 

            3LL*Power(xi,14LL)*(45LL + 75LL*r*xj + 60LL*Power(r,2LL)*Power(xj,2LL) + 

               30LL*Power(r,3LL)*Power(xj,3LL) + 10LL*Power(r,4LL)*Power(xj,4LL) + 

               2LL*Power(r,5LL)*Power(xj,5LL)) - 

            55LL*Power(xi,8LL)*Power(xj,6LL)*

             (-405LL - 567LL*r*xj - 972LL*Power(r,2LL)*Power(xj,2LL) - 

               90LL*Power(r,3LL)*Power(xj,3LL) + 18LL*Power(r,4LL)*Power(xj,4LL) + 

               2LL*Power(r,5LL)*Power(xj,5LL)) + 

            55LL*Power(xi,6LL)*Power(xj,8LL)*

             (9LL - 4257LL*r*xj - 372LL*Power(r,2LL)*Power(xj,2LL) + 

               222LL*Power(r,3LL)*Power(xj,3LL) + 42LL*Power(r,4LL)*Power(xj,4LL) + 

               2LL*Power(r,5LL)*Power(xj,5LL)) + 

            3LL*Power(xj,14LL)*(15015LL + 10725LL*r*xj + 3300LL*Power(r,2LL)*Power(xj,2LL) + 

               550LL*Power(r,3LL)*Power(xj,3LL) + 50LL*Power(r,4LL)*Power(xj,4LL) + 

               2LL*Power(r,5LL)*Power(xj,5LL)) + 

            5LL*Power(xi,12LL)*Power(xj,2LL)*

             (297LL + 495LL*r*xj + 396LL*Power(r,2LL)*Power(xj,2LL) + 

               198LL*Power(r,3LL)*Power(xj,3LL) + 66LL*Power(r,4LL)*Power(xj,4LL) + 

               2LL*Power(r,5LL)*Power(xj,5LL)) + 

            Power(xi,10LL)*Power(xj,4LL)*

             (-7425LL - 12375LL*r*xj - 9900LL*Power(r,2LL)*Power(xj,2LL) - 

               6210LL*Power(r,3LL)*Power(xj,3LL) - 390LL*Power(r,4LL)*Power(xj,4LL) + 

               34LL*Power(r,5LL)*Power(xj,5LL)) - 

            Power(xi,4LL)*Power(xj,10LL)*

             (-484155LL + 38475LL*r*xj + 78780LL*Power(r,2LL)*Power(xj,2LL) + 

               17190LL*Power(r,3LL)*Power(xj,3LL) + 1410LL*Power(r,4LL)*Power(xj,4LL) + 

               34LL*Power(r,5LL)*Power(xj,5LL))))/

       (135LL*exp(2LL*r*(xi + xj))*Power(r,2LL)*Power(xi - xj,11LL)*

         Power(xi + xj,11LL)) + (2LL*(135LL*exp(2LL*r*(xi + xj))*

            Power(Power(xi,2LL) - Power(xj,2LL),11LL) + 

           exp(2LL*r*xj)*Power(xj,8LL)*

            (-150LL*Power(r,4LL)*Power(xi,18LL) - 6LL*Power(r,5LL)*Power(xi,19LL) + 

              135LL*Power(xj,14LL) + 225LL*r*xi*Power(xj,14LL) + 

              10LL*Power(r,3LL)*Power(xi,17LL)*(-165LL + Power(r,2LL)*Power(xj,2LL)) - 

              30LL*Power(r,2LL)*Power(xi,16LL)*(330LL + Power(r,2LL)*Power(xj,2LL)) + 

              45LL*r*Power(xi,3LL)*Power(xj,12LL)*(-55LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

              45LL*Power(xi,2LL)*Power(xj,12LL)*(-33LL + 4LL*Power(r,2LL)*Power(xj,2LL)) + 

              r*Power(xi,9LL)*Power(xj,6LL)*

               (234135LL - 4950LL*Power(r,2LL)*Power(xj,2LL) - 

                 34LL*Power(r,4LL)*Power(xj,4LL)) - 

              5LL*r*Power(xi,7LL)*Power(xj,8LL)*

               (6237LL - 1242LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) 

    + 3LL*r*Power(xi,5LL)*Power(xj,10LL)*

               (4125LL - 330LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) 

    + 15LL*Power(xi,4LL)*Power(xj,10LL)*

               (495LL - 132LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) - 

              165LL*Power(xi,6LL)*Power(xj,8LL)*

               (135LL - 60LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) - 

              5LL*r*Power(xi,13LL)*Power(xj,2LL)*

               (43875LL - 3438LL*Power(r,2LL)*Power(xj,2LL) + 

                 22LL*Power(r,4LL)*Power(xj,4LL)) + 

              5LL*r*Power(xi,11LL)*Power(xj,4LL)*

               (7695LL - 2442LL*Power(r,2LL)*Power(xj,2LL) + 

                 22LL*Power(r,4LL)*Power(xj,4LL)) + 

              15LL*Power(xi,8LL)*Power(xj,6LL)*

               (-33LL - 3564LL*Power(r,2LL)*Power(xj,2LL) + 26LL*Power(r,4LL)*Power(xj,4LL)) 

    + r*Power(xi,15LL)*(-32175LL - 3690LL*Power(r,2LL)*Power(xj,2LL) + 

                 34LL*Power(r,4LL)*Power(xj,4LL)) + 

              15LL*Power(xi,10LL)*Power(xj,4LL)*

               (-32277LL + 1364LL*Power(r,2LL)*Power(xj,2LL) + 

                 66LL*Power(r,4LL)*Power(xj,4LL)) + 

              15LL*Power(xi,14LL)*(-3003LL - 2932LL*Power(r,2LL)*Power(xj,2LL) + 

                 94LL*Power(r,4LL)*Power(xj,4LL)) - 

              15LL*Power(xi,12LL)*Power(xj,2LL)*

               (28119LL - 5252LL*Power(r,2LL)*Power(xj,2LL) + 

                 154LL*Power(r,4LL)*Power(xj,4LL))) + 

           exp(2LL*r*xi)*Power(xi,8LL)*

            (-5LL*Power(xi,2LL)*Power(xj,12LL)*

               (-84357LL - 43875LL*r*xj - 8796LL*Power(r,2LL)*Power(xj,2LL) - 

                 738LL*Power(r,3LL)*Power(xj,3LL) - 6LL*Power(r,4LL)*Power(xj,4LL) + 

                 2LL*Power(r,5LL)*Power(xj,5LL)) - 

              3LL*Power(xi,14LL)*(45LL + 75LL*r*xj + 60LL*Power(r,2LL)*Power(xj,2LL) + 

                 30LL*Power(r,3LL)*Power(xj,3LL) + 10LL*Power(r,4LL)*Power(xj,4LL) + 

                 2LL*Power(r,5LL)*Power(xj,5LL)) - 

              55LL*Power(xi,8LL)*Power(xj,6LL)*

               (-405LL - 567LL*r*xj - 972LL*Power(r,2LL)*Power(xj,2LL) - 

                 90LL*Power(r,3LL)*Power(xj,3LL) + 18LL*Power(r,4LL)*Power(xj,4LL) + 

                 2LL*Power(r,5LL)*Power(xj,5LL)) + 

              55LL*Power(xi,6LL)*Power(xj,8LL)*

               (9LL - 4257LL*r*xj - 372LL*Power(r,2LL)*Power(xj,2LL) + 

                 222LL*Power(r,3LL)*Power(xj,3LL) + 42LL*Power(r,4LL)*Power(xj,4LL) + 

                 2LL*Power(r,5LL)*Power(xj,5LL)) + 

              3LL*Power(xj,14LL)*(15015LL + 10725LL*r*xj + 

                 3300LL*Power(r,2LL)*Power(xj,2LL) + 550LL*Power(r,3LL)*Power(xj,3LL) + 

                 50LL*Power(r,4LL)*Power(xj,4LL) + 2LL*Power(r,5LL)*Power(xj,5LL)) + 

              5LL*Power(xi,12LL)*Power(xj,2LL)*

               (297LL + 495LL*r*xj + 396LL*Power(r,2LL)*Power(xj,2LL) + 

                 198LL*Power(r,3LL)*Power(xj,3LL) + 66LL*Power(r,4LL)*Power(xj,4LL) + 

                 2LL*Power(r,5LL)*Power(xj,5LL)) + 

              Power(xi,10LL)*Power(xj,4LL)*

               (-7425LL - 12375LL*r*xj - 9900LL*Power(r,2LL)*Power(xj,2LL) - 

                 6210LL*Power(r,3LL)*Power(xj,3LL) - 390LL*Power(r,4LL)*Power(xj,4LL) + 

                 34LL*Power(r,5LL)*Power(xj,5LL)) - 

              Power(xi,4LL)*Power(xj,10LL)*

               (-484155LL + 38475LL*r*xj + 78780LL*Power(r,2LL)*Power(xj,2LL) + 

                 17190LL*Power(r,3LL)*Power(xj,3LL) + 1410LL*Power(r,4LL)*Power(xj,4LL) + 

                 34LL*Power(r,5LL)*Power(xj,5LL)))))/

       (135LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,11LL)*Power(xi + xj,10LL)) - 

      (270LL*exp(2LL*r*(xi + xj))*(xi + xj)*

          Power(Power(xi,2LL) - Power(xj,2LL),11LL) + 

         exp(2LL*r*xj)*Power(xj,8LL)*

          (-600LL*Power(r,3LL)*Power(xi,18LL) - 30LL*Power(r,4LL)*Power(xi,19LL) - 

            60LL*Power(r,3LL)*Power(xi,16LL)*Power(xj,2LL) + 

            20LL*Power(r,4LL)*Power(xi,17LL)*Power(xj,2LL) + 225LL*xi*Power(xj,14LL) + 

            360LL*r*Power(xi,2LL)*Power(xj,14LL) + 

            180LL*Power(r,2LL)*Power(xi,3LL)*Power(xj,14LL) + 

            30LL*Power(r,2LL)*Power(xi,17LL)*(-165LL + Power(r,2LL)*Power(xj,2LL)) - 

            60LL*r*Power(xi,16LL)*(330LL + Power(r,2LL)*Power(xj,2LL)) + 

            45LL*Power(xi,3LL)*Power(xj,12LL)*(-55LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            r*Power(xi,9LL)*Power(xj,6LL)*

             (-9900LL*r*Power(xj,2LL) - 136LL*Power(r,3LL)*Power(xj,4LL)) - 

            5LL*r*Power(xi,7LL)*Power(xj,8LL)*

             (-2484LL*r*Power(xj,2LL) + 8LL*Power(r,3LL)*Power(xj,4LL)) + 

            3LL*r*Power(xi,5LL)*Power(xj,10LL)*

             (-660LL*r*Power(xj,2LL) + 8LL*Power(r,3LL)*Power(xj,4LL)) + 

            15LL*Power(xi,4LL)*Power(xj,10LL)*

             (-264LL*r*Power(xj,2LL) + 8LL*Power(r,3LL)*Power(xj,4LL)) - 

            165LL*Power(xi,6LL)*Power(xj,8LL)*

             (-120LL*r*Power(xj,2LL) + 8LL*Power(r,3LL)*Power(xj,4LL)) - 

            5LL*r*Power(xi,13LL)*Power(xj,2LL)*

             (-6876LL*r*Power(xj,2LL) + 88LL*Power(r,3LL)*Power(xj,4LL)) + 

            5LL*r*Power(xi,11LL)*Power(xj,4LL)*

             (-4884LL*r*Power(xj,2LL) + 88LL*Power(r,3LL)*Power(xj,4LL)) + 

            15LL*Power(xi,8LL)*Power(xj,6LL)*

             (-7128LL*r*Power(xj,2LL) + 104LL*Power(r,3LL)*Power(xj,4LL)) + 

            r*Power(xi,15LL)*(-7380LL*r*Power(xj,2LL) + 136LL*Power(r,3LL)*Power(xj,4LL)) + 

            15LL*Power(xi,10LL)*Power(xj,4LL)*

             (2728LL*r*Power(xj,2LL) + 264LL*Power(r,3LL)*Power(xj,4LL)) + 

            15LL*Power(xi,14LL)*(-5864LL*r*Power(xj,2LL) + 

               376LL*Power(r,3LL)*Power(xj,4LL)) - 

            15LL*Power(xi,12LL)*Power(xj,2LL)*

             (-10504LL*r*Power(xj,2LL) + 616LL*Power(r,3LL)*Power(xj,4LL)) + 

            Power(xi,9LL)*Power(xj,6LL)*

             (234135LL - 4950LL*Power(r,2LL)*Power(xj,2LL) - 34LL*Power(r,4LL)*Power(xj,4LL)) 

    - 5LL*Power(xi,7LL)*Power(xj,8LL)*(6237LL - 1242LL*Power(r,2LL)*Power(xj,2LL) + 

               2LL*Power(r,4LL)*Power(xj,4LL)) + 

            3LL*Power(xi,5LL)*Power(xj,10LL)*

             (4125LL - 330LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) - 

            5LL*Power(xi,13LL)*Power(xj,2LL)*

             (43875LL - 3438LL*Power(r,2LL)*Power(xj,2LL) + 22LL*Power(r,4LL)*Power(xj,4LL)) 

    + 5LL*Power(xi,11LL)*Power(xj,4LL)*(7695LL - 2442LL*Power(r,2LL)*Power(xj,2LL) + 

               22LL*Power(r,4LL)*Power(xj,4LL)) + 

            Power(xi,15LL)*(-32175LL - 3690LL*Power(r,2LL)*Power(xj,2LL) + 

               34LL*Power(r,4LL)*Power(xj,4LL))) + 

         2LL*exp(2LL*r*xj)*Power(xj,9LL)*

          (-150LL*Power(r,4LL)*Power(xi,18LL) - 6LL*Power(r,5LL)*Power(xi,19LL) + 

            135LL*Power(xj,14LL) + 225LL*r*xi*Power(xj,14LL) + 

            10LL*Power(r,3LL)*Power(xi,17LL)*(-165LL + Power(r,2LL)*Power(xj,2LL)) - 

            30LL*Power(r,2LL)*Power(xi,16LL)*(330LL + Power(r,2LL)*Power(xj,2LL)) + 

            45LL*r*Power(xi,3LL)*Power(xj,12LL)*(-55LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            45LL*Power(xi,2LL)*Power(xj,12LL)*(-33LL + 4LL*Power(r,2LL)*Power(xj,2LL)) + 

            r*Power(xi,9LL)*Power(xj,6LL)*

             (234135LL - 4950LL*Power(r,2LL)*Power(xj,2LL) - 34LL*Power(r,4LL)*Power(xj,4LL)) 

    - 5LL*r*Power(xi,7LL)*Power(xj,8LL)*(6237LL - 1242LL*Power(r,2LL)*Power(xj,2LL) + 

               2LL*Power(r,4LL)*Power(xj,4LL)) + 

            3LL*r*Power(xi,5LL)*Power(xj,10LL)*

             (4125LL - 330LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            15LL*Power(xi,4LL)*Power(xj,10LL)*

             (495LL - 132LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) - 

            165LL*Power(xi,6LL)*Power(xj,8LL)*

             (135LL - 60LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) - 

            5LL*r*Power(xi,13LL)*Power(xj,2LL)*

             (43875LL - 3438LL*Power(r,2LL)*Power(xj,2LL) + 22LL*Power(r,4LL)*Power(xj,4LL)) 

    + 5LL*r*Power(xi,11LL)*Power(xj,4LL)*

             (7695LL - 2442LL*Power(r,2LL)*Power(xj,2LL) + 22LL*Power(r,4LL)*Power(xj,4LL)) + 

            15LL*Power(xi,8LL)*Power(xj,6LL)*

             (-33LL - 3564LL*Power(r,2LL)*Power(xj,2LL) + 26LL*Power(r,4LL)*Power(xj,4LL)) + 

            r*Power(xi,15LL)*(-32175LL - 3690LL*Power(r,2LL)*Power(xj,2LL) + 

               34LL*Power(r,4LL)*Power(xj,4LL)) + 

            15LL*Power(xi,10LL)*Power(xj,4LL)*

             (-32277LL + 1364LL*Power(r,2LL)*Power(xj,2LL) + 66LL*Power(r,4LL)*Power(xj,4LL)) 

    + 15LL*Power(xi,14LL)*(-3003LL - 2932LL*Power(r,2LL)*Power(xj,2LL) + 

               94LL*Power(r,4LL)*Power(xj,4LL)) - 

            15LL*Power(xi,12LL)*Power(xj,2LL)*

             (28119LL - 5252LL*Power(r,2LL)*Power(xj,2LL) + 154LL*Power(r,4LL)*Power(xj,4LL))) 

    + exp(2LL*r*xi)*Power(xi,8LL)*(-5LL*Power(xi,2LL)*Power(xj,12LL)*

             (-43875LL*xj - 17592LL*r*Power(xj,2LL) - 2214LL*Power(r,2LL)*Power(xj,3LL) - 

               24LL*Power(r,3LL)*Power(xj,4LL) + 10LL*Power(r,4LL)*Power(xj,5LL)) - 

            3LL*Power(xi,14LL)*(75LL*xj + 120LL*r*Power(xj,2LL) + 

               90LL*Power(r,2LL)*Power(xj,3LL) + 40LL*Power(r,3LL)*Power(xj,4LL) + 

               10LL*Power(r,4LL)*Power(xj,5LL)) - 

            55LL*Power(xi,8LL)*Power(xj,6LL)*

             (-567LL*xj - 1944LL*r*Power(xj,2LL) - 270LL*Power(r,2LL)*Power(xj,3LL) + 

               72LL*Power(r,3LL)*Power(xj,4LL) + 10LL*Power(r,4LL)*Power(xj,5LL)) + 

            55LL*Power(xi,6LL)*Power(xj,8LL)*

             (-4257LL*xj - 744LL*r*Power(xj,2LL) + 666LL*Power(r,2LL)*Power(xj,3LL) + 

               168LL*Power(r,3LL)*Power(xj,4LL) + 10LL*Power(r,4LL)*Power(xj,5LL)) + 

            3LL*Power(xj,14LL)*(10725LL*xj + 6600LL*r*Power(xj,2LL) + 

               1650LL*Power(r,2LL)*Power(xj,3LL) + 200LL*Power(r,3LL)*Power(xj,4LL) + 

               10LL*Power(r,4LL)*Power(xj,5LL)) + 

            5LL*Power(xi,12LL)*Power(xj,2LL)*

             (495LL*xj + 792LL*r*Power(xj,2LL) + 594LL*Power(r,2LL)*Power(xj,3LL) + 

               264LL*Power(r,3LL)*Power(xj,4LL) + 10LL*Power(r,4LL)*Power(xj,5LL)) + 

            Power(xi,10LL)*Power(xj,4LL)*

             (-12375LL*xj - 19800LL*r*Power(xj,2LL) - 18630LL*Power(r,2LL)*Power(xj,3LL) - 

               1560LL*Power(r,3LL)*Power(xj,4LL) + 170LL*Power(r,4LL)*Power(xj,5LL)) - 

            Power(xi,4LL)*Power(xj,10LL)*

             (38475LL*xj + 157560LL*r*Power(xj,2LL) + 51570LL*Power(r,2LL)*Power(xj,3LL) + 

               5640LL*Power(r,3LL)*Power(xj,4LL) + 170LL*Power(r,4LL)*Power(xj,5LL))) + 

         2LL*exp(2LL*r*xi)*Power(xi,9LL)*

          (-5LL*Power(xi,2LL)*Power(xj,12LL)*

             (-84357LL - 43875LL*r*xj - 8796LL*Power(r,2LL)*Power(xj,2LL) - 

               738LL*Power(r,3LL)*Power(xj,3LL) - 6LL*Power(r,4LL)*Power(xj,4LL) + 

               2LL*Power(r,5LL)*Power(xj,5LL)) - 

            3LL*Power(xi,14LL)*(45LL + 75LL*r*xj + 60LL*Power(r,2LL)*Power(xj,2LL) + 

               30LL*Power(r,3LL)*Power(xj,3LL) + 10LL*Power(r,4LL)*Power(xj,4LL) + 

               2LL*Power(r,5LL)*Power(xj,5LL)) - 

            55LL*Power(xi,8LL)*Power(xj,6LL)*

             (-405LL - 567LL*r*xj - 972LL*Power(r,2LL)*Power(xj,2LL) - 

               90LL*Power(r,3LL)*Power(xj,3LL) + 18LL*Power(r,4LL)*Power(xj,4LL) + 

               2LL*Power(r,5LL)*Power(xj,5LL)) + 

            55LL*Power(xi,6LL)*Power(xj,8LL)*

             (9LL - 4257LL*r*xj - 372LL*Power(r,2LL)*Power(xj,2LL) + 

               222LL*Power(r,3LL)*Power(xj,3LL) + 42LL*Power(r,4LL)*Power(xj,4LL) + 

               2LL*Power(r,5LL)*Power(xj,5LL)) + 

            3LL*Power(xj,14LL)*(15015LL + 10725LL*r*xj + 3300LL*Power(r,2LL)*Power(xj,2LL) + 

               550LL*Power(r,3LL)*Power(xj,3LL) + 50LL*Power(r,4LL)*Power(xj,4LL) + 

               2LL*Power(r,5LL)*Power(xj,5LL)) + 

            5LL*Power(xi,12LL)*Power(xj,2LL)*

             (297LL + 495LL*r*xj + 396LL*Power(r,2LL)*Power(xj,2LL) + 

               198LL*Power(r,3LL)*Power(xj,3LL) + 66LL*Power(r,4LL)*Power(xj,4LL) + 

               2LL*Power(r,5LL)*Power(xj,5LL)) + 

            Power(xi,10LL)*Power(xj,4LL)*

             (-7425LL - 12375LL*r*xj - 9900LL*Power(r,2LL)*Power(xj,2LL) - 

               6210LL*Power(r,3LL)*Power(xj,3LL) - 390LL*Power(r,4LL)*Power(xj,4LL) + 

               34LL*Power(r,5LL)*Power(xj,5LL)) - 

            Power(xi,4LL)*Power(xj,10LL)*

             (-484155LL + 38475LL*r*xj + 78780LL*Power(r,2LL)*Power(xj,2LL) + 

               17190LL*Power(r,3LL)*Power(xj,3LL) + 1410LL*Power(r,4LL)*Power(xj,4LL) + 

               34LL*Power(r,5LL)*Power(xj,5LL))))/

       (135LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,11LL)*Power(xi + xj,11LL))

    ; }
   
  }
  return S;
}

cl_R DSlater_3S_4S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = 0LL

    ; } else {  S = -(-132871488750LL*xi + 149448499200LL*exp(2LL*r*xi)*xi - 

          232588956600LL*r*Power(xi,2LL) - 200036962125LL*Power(r,2LL)*Power(xi,3LL) - 

          112459347000LL*Power(r,3LL)*Power(xi,4LL) - 

          46370223900LL*Power(r,4LL)*Power(xi,5LL) - 

          14905931040LL*Power(r,5LL)*Power(xi,6LL) - 

          3872428560LL*Power(r,6LL)*Power(xi,7LL) - 

          830269440LL*Power(r,7LL)*Power(xi,8LL) - 148262400LL*Power(r,8LL)*Power(xi,9LL) - 

          21964800LL*Power(r,9LL)*Power(xi,10LL) - 2635776LL*Power(r,10LL)*Power(xi,11LL) - 

          239616LL*Power(r,11LL)*Power(xi,12LL) - 13312LL*Power(r,12LL)*Power(xi,13LL))/

       (7.47242496e10*exp(2LL*r*xi)*r) + 

      (-74724249600LL + 74724249600LL*exp(2LL*r*xi) - 132871488750LL*r*xi - 

         116294478300LL*Power(r,2LL)*Power(xi,2LL) - 

         66678987375LL*Power(r,3LL)*Power(xi,3LL) - 

         28114836750LL*Power(r,4LL)*Power(xi,4LL) - 

         9274044780LL*Power(r,5LL)*Power(xi,5LL) - 

         2484321840LL*Power(r,6LL)*Power(xi,6LL) - 553204080LL*Power(r,7LL)*Power(xi,7LL) - 

         103783680LL*Power(r,8LL)*Power(xi,8LL) - 16473600LL*Power(r,9LL)*Power(xi,9LL) - 

         2196480LL*Power(r,10LL)*Power(xi,10LL) - 239616LL*Power(r,11LL)*Power(xi,11LL) - 

         19968LL*Power(r,12LL)*Power(xi,12LL) - 1024LL*Power(r,13LL)*Power(xi,13LL))/

       (7.47242496e10*exp(2LL*r*xi)*Power(r,2LL)) + 

      (xi*(-74724249600LL + 74724249600LL*exp(2LL*r*xi) - 132871488750LL*r*xi - 

           116294478300LL*Power(r,2LL)*Power(xi,2LL) - 

           66678987375LL*Power(r,3LL)*Power(xi,3LL) - 

           28114836750LL*Power(r,4LL)*Power(xi,4LL) - 

           9274044780LL*Power(r,5LL)*Power(xi,5LL) - 

           2484321840LL*Power(r,6LL)*Power(xi,6LL) - 

           553204080LL*Power(r,7LL)*Power(xi,7LL) - 103783680LL*Power(r,8LL)*Power(xi,8LL) - 

           16473600LL*Power(r,9LL)*Power(xi,9LL) - 2196480LL*Power(r,10LL)*Power(xi,10LL) - 

           239616LL*Power(r,11LL)*Power(xi,11LL) - 19968LL*Power(r,12LL)*Power(xi,12LL) - 

           1024LL*Power(r,13LL)*Power(xi,13LL)))/(3.73621248e10*exp(2LL*r*xi)*r)

    ; }
 
  }
  else {
      if (r == 0LL) {  S = 0LL

    ; } else {  S = (3780LL*exp(2LL*r*(xi + xj))*Power(Power(xi,2LL) - Power(xj,2LL),13LL) + 

         84LL*exp(2LL*r*xj)*Power(xj,10LL)*

          (-60LL*Power(r,4LL)*Power(xi,20LL) - 2LL*Power(r,5LL)*Power(xi,21LL) + 

            45LL*Power(xj,16LL) + 75LL*r*xi*Power(xj,16LL) - 

            4LL*Power(r,3LL)*Power(xi,19LL)*(195LL + Power(r,2LL)*Power(xj,2LL)) + 

            15LL*r*Power(xi,3LL)*Power(xj,14LL)*(-65LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            15LL*Power(xi,2LL)*Power(xj,14LL)*(-39LL + 4LL*Power(r,2LL)*Power(xj,2LL)) - 

            30LL*Power(r,2LL)*Power(xi,18LL)*(182LL + 9LL*Power(r,2LL)*Power(xj,2LL)) + 

            30LL*r*Power(xi,13LL)*Power(xj,4LL)*

             (-13047LL + 377LL*Power(r,2LL)*Power(xj,2LL)) + 

            2LL*r*Power(xi,5LL)*Power(xj,12LL)*

             (2925LL - 195LL*Power(r,2LL)*Power(xj,2LL) + Power(r,4LL)*Power(xj,4LL)) + 

            10LL*Power(xi,4LL)*Power(xj,12LL)*

             (351LL - 78LL*Power(r,2LL)*Power(xj,2LL) + Power(r,4LL)*Power(xj,4LL)) - 

            130LL*Power(xi,6LL)*Power(xj,10LL)*

             (99LL - 36LL*Power(r,2LL)*Power(xj,2LL) + Power(r,4LL)*Power(xj,4LL)) + 

            13LL*r*Power(xi,11LL)*Power(xj,6LL)*

             (30735LL - 1650LL*Power(r,2LL)*Power(xj,2LL) + 4LL*Power(r,4LL)*Power(xj,4LL)) 

    + r*Power(xi,7LL)*Power(xj,10LL)*(-15015LL + 3330LL*Power(r,2LL)*Power(xj,2LL) + 

               4LL*Power(r,4LL)*Power(xj,4LL)) + 

            210LL*Power(xi,16LL)*(-156LL - 262LL*Power(r,2LL)*Power(xj,2LL) + 

               5LL*Power(r,4LL)*Power(xj,4LL)) - 

            6LL*r*Power(xi,9LL)*Power(xj,8LL)*

             (-48620LL - 715LL*Power(r,2LL)*Power(xj,2LL) + 6LL*Power(r,4LL)*Power(xj,4LL)) 

    + 3LL*r*Power(xi,17LL)*(-6825LL - 1870LL*Power(r,2LL)*Power(xj,2LL) + 

               12LL*Power(r,4LL)*Power(xj,4LL)) - 

            30LL*Power(xi,14LL)*Power(xj,2LL)*

             (17934LL - 12LL*Power(r,2LL)*Power(xj,2LL) + 13LL*Power(r,4LL)*Power(xj,4LL)) - 

            15LL*Power(xi,8LL)*Power(xj,8LL)*

             (2145LL + 2860LL*Power(r,2LL)*Power(xj,2LL) + 14LL*Power(r,4LL)*Power(xj,4LL)) 

    + 65LL*Power(xi,10LL)*Power(xj,6LL)*

             (-13725LL - 792LL*Power(r,2LL)*Power(xj,2LL) + 22LL*Power(r,4LL)*Power(xj,4LL)) 

    - 10LL*Power(xi,12LL)*Power(xj,4LL)*

             (153630LL - 15054LL*Power(r,2LL)*Power(xj,2LL) + 

               143LL*Power(r,4LL)*Power(xj,4LL)) + 

            Power(xi,15LL)*(-269325LL*r*Power(xj,2LL) + 

               9270LL*Power(r,3LL)*Power(xj,4LL) - 52LL*Power(r,5LL)*Power(xj,6LL))) + 

         exp(2LL*r*xi)*Power(xi,8LL)*

          (Power(xi,2LL)*Power(xj,16LL)*

             (70073640LL + 47669895LL*r*xj + 13931190LL*Power(r,2LL)*Power(xj,2LL) + 

               2170350LL*Power(r,3LL)*Power(xj,3LL) + 

               169260LL*Power(r,4LL)*Power(xj,4LL) + 1638LL*Power(r,5LL)*Power(xj,5LL) - 

               756LL*Power(r,6LL)*Power(xj,6LL) - 44LL*Power(r,7LL)*Power(xj,7LL)) + 

            364LL*Power(xi,10LL)*Power(xj,8LL)*

             (-7425LL - 13860LL*r*xj - 5940LL*Power(r,2LL)*Power(xj,2LL) - 

               11880LL*Power(r,3LL)*Power(xj,3LL) - 2640LL*Power(r,4LL)*Power(xj,4LL) - 

               45LL*Power(r,5LL)*Power(xj,5LL) + 30LL*Power(r,6LL)*Power(xj,6LL) + 

               2LL*Power(r,7LL)*Power(xj,7LL)) - 

            364LL*Power(xi,8LL)*Power(xj,10LL)*

             (-20925LL + 18270LL*r*xj - 58320LL*Power(r,2LL)*Power(xj,2LL) - 

               17730LL*Power(r,3LL)*Power(xj,3LL) - 300LL*Power(r,4LL)*Power(xj,4LL) + 

               423LL*Power(r,5LL)*Power(xj,5LL) + 54LL*Power(r,6LL)*Power(xj,6LL) + 

               2LL*Power(r,7LL)*Power(xj,7LL)) - 

            3LL*Power(xi,18LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r,2LL)*Power(xj,2LL) + 

               1050LL*Power(r,3LL)*Power(xj,3LL) + 420LL*Power(r,4LL)*Power(xj,4LL) + 

               126LL*Power(r,5LL)*Power(xj,5LL) + 28LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL)) + 

            3LL*Power(xj,18LL)*(1801800LL + 1576575LL*r*xj + 

               630630LL*Power(r,2LL)*Power(xj,2LL) + 150150LL*Power(r,3LL)*Power(xj,3LL) + 

               23100LL*Power(r,4LL)*Power(xj,4LL) + 2310LL*Power(r,5LL)*Power(xj,5LL) + 

               140LL*Power(r,6LL)*Power(xj,6LL) + 4LL*Power(r,7LL)*Power(xj,7LL)) + 

            2LL*Power(xi,14LL)*Power(xj,4LL)*

             (-147420LL - 257985LL*r*xj - 221130LL*Power(r,2LL)*Power(xj,2LL) - 

               122850LL*Power(r,3LL)*Power(xj,3LL) - 49140LL*Power(r,4LL)*Power(xj,4LL) - 

               17388LL*Power(r,5LL)*Power(xj,5LL) - 1512LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) - 

            42LL*Power(xi,12LL)*Power(xj,6LL)*

             (-25740LL - 45045LL*r*xj - 38610LL*Power(r,2LL)*Power(xj,2LL) - 

               19470LL*Power(r,3LL)*Power(xj,3LL) - 12540LL*Power(r,4LL)*Power(xj,4LL) - 

               1836LL*Power(r,5LL)*Power(xj,5LL) - 8LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) + 

            42LL*Power(xi,6LL)*Power(xj,12LL)*

             (921600LL - 1640835LL*r*xj - 546030LL*Power(r,2LL)*Power(xj,2LL) + 

               20730LL*Power(r,3LL)*Power(xj,3LL) + 30180LL*Power(r,4LL)*Power(xj,4LL) + 

               5028LL*Power(r,5LL)*Power(xj,5LL) + 344LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) - 

            2LL*Power(xi,4LL)*Power(xj,14LL)*

             (-67767840LL - 13377735LL*r*xj + 6601770LL*Power(r,2LL)*Power(xj,2LL) + 

               3115350LL*Power(r,3LL)*Power(xj,3LL) + 

               548940LL*Power(r,4LL)*Power(xj,4LL) + 48132LL*Power(r,5LL)*Power(xj,5LL) + 

               1848LL*Power(r,6LL)*Power(xj,6LL) + 8LL*Power(r,7LL)*Power(xj,7LL)) + 

            Power(xi,16LL)*Power(xj,2LL)*

             (49140LL + 85995LL*r*xj + 73710LL*Power(r,2LL)*Power(xj,2LL) + 

               40950LL*Power(r,3LL)*Power(xj,3LL) + 16380LL*Power(r,4LL)*Power(xj,4LL) + 

               4914LL*Power(r,5LL)*Power(xj,5LL) + 1092LL*Power(r,6LL)*Power(xj,6LL) + 

               44LL*Power(r,7LL)*Power(xj,7LL))))/

       (3780LL*exp(2LL*r*(xi + xj))*Power(r,2LL)*Power(xi - xj,13LL)*

         Power(xi + xj,13LL)) + (3780LL*exp(2LL*r*(xi + xj))*

          Power(Power(xi,2LL) - Power(xj,2LL),13LL) + 

         84LL*exp(2LL*r*xj)*Power(xj,10LL)*

          (-60LL*Power(r,4LL)*Power(xi,20LL) - 2LL*Power(r,5LL)*Power(xi,21LL) + 

            45LL*Power(xj,16LL) + 75LL*r*xi*Power(xj,16LL) - 

            4LL*Power(r,3LL)*Power(xi,19LL)*(195LL + Power(r,2LL)*Power(xj,2LL)) + 

            15LL*r*Power(xi,3LL)*Power(xj,14LL)*(-65LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            15LL*Power(xi,2LL)*Power(xj,14LL)*(-39LL + 4LL*Power(r,2LL)*Power(xj,2LL)) - 

            30LL*Power(r,2LL)*Power(xi,18LL)*(182LL + 9LL*Power(r,2LL)*Power(xj,2LL)) + 

            30LL*r*Power(xi,13LL)*Power(xj,4LL)*

             (-13047LL + 377LL*Power(r,2LL)*Power(xj,2LL)) + 

            2LL*r*Power(xi,5LL)*Power(xj,12LL)*

             (2925LL - 195LL*Power(r,2LL)*Power(xj,2LL) + Power(r,4LL)*Power(xj,4LL)) + 

            10LL*Power(xi,4LL)*Power(xj,12LL)*

             (351LL - 78LL*Power(r,2LL)*Power(xj,2LL) + Power(r,4LL)*Power(xj,4LL)) - 

            130LL*Power(xi,6LL)*Power(xj,10LL)*

             (99LL - 36LL*Power(r,2LL)*Power(xj,2LL) + Power(r,4LL)*Power(xj,4LL)) + 

            13LL*r*Power(xi,11LL)*Power(xj,6LL)*

             (30735LL - 1650LL*Power(r,2LL)*Power(xj,2LL) + 4LL*Power(r,4LL)*Power(xj,4LL)) 

    + r*Power(xi,7LL)*Power(xj,10LL)*(-15015LL + 3330LL*Power(r,2LL)*Power(xj,2LL) + 

               4LL*Power(r,4LL)*Power(xj,4LL)) + 

            210LL*Power(xi,16LL)*(-156LL - 262LL*Power(r,2LL)*Power(xj,2LL) + 

               5LL*Power(r,4LL)*Power(xj,4LL)) - 

            6LL*r*Power(xi,9LL)*Power(xj,8LL)*

             (-48620LL - 715LL*Power(r,2LL)*Power(xj,2LL) + 6LL*Power(r,4LL)*Power(xj,4LL)) 

    + 3LL*r*Power(xi,17LL)*(-6825LL - 1870LL*Power(r,2LL)*Power(xj,2LL) + 

               12LL*Power(r,4LL)*Power(xj,4LL)) - 

            30LL*Power(xi,14LL)*Power(xj,2LL)*

             (17934LL - 12LL*Power(r,2LL)*Power(xj,2LL) + 13LL*Power(r,4LL)*Power(xj,4LL)) - 

            15LL*Power(xi,8LL)*Power(xj,8LL)*

             (2145LL + 2860LL*Power(r,2LL)*Power(xj,2LL) + 14LL*Power(r,4LL)*Power(xj,4LL)) 

    + 65LL*Power(xi,10LL)*Power(xj,6LL)*

             (-13725LL - 792LL*Power(r,2LL)*Power(xj,2LL) + 22LL*Power(r,4LL)*Power(xj,4LL)) 

    - 10LL*Power(xi,12LL)*Power(xj,4LL)*

             (153630LL - 15054LL*Power(r,2LL)*Power(xj,2LL) + 

               143LL*Power(r,4LL)*Power(xj,4LL)) + 

            Power(xi,15LL)*(-269325LL*r*Power(xj,2LL) + 

               9270LL*Power(r,3LL)*Power(xj,4LL) - 52LL*Power(r,5LL)*Power(xj,6LL))) + 

         exp(2LL*r*xi)*Power(xi,8LL)*

          (Power(xi,2LL)*Power(xj,16LL)*

             (70073640LL + 47669895LL*r*xj + 13931190LL*Power(r,2LL)*Power(xj,2LL) + 

               2170350LL*Power(r,3LL)*Power(xj,3LL) + 

               169260LL*Power(r,4LL)*Power(xj,4LL) + 1638LL*Power(r,5LL)*Power(xj,5LL) - 

               756LL*Power(r,6LL)*Power(xj,6LL) - 44LL*Power(r,7LL)*Power(xj,7LL)) + 

            364LL*Power(xi,10LL)*Power(xj,8LL)*

             (-7425LL - 13860LL*r*xj - 5940LL*Power(r,2LL)*Power(xj,2LL) - 

               11880LL*Power(r,3LL)*Power(xj,3LL) - 2640LL*Power(r,4LL)*Power(xj,4LL) - 

               45LL*Power(r,5LL)*Power(xj,5LL) + 30LL*Power(r,6LL)*Power(xj,6LL) + 

               2LL*Power(r,7LL)*Power(xj,7LL)) - 

            364LL*Power(xi,8LL)*Power(xj,10LL)*

             (-20925LL + 18270LL*r*xj - 58320LL*Power(r,2LL)*Power(xj,2LL) - 

               17730LL*Power(r,3LL)*Power(xj,3LL) - 300LL*Power(r,4LL)*Power(xj,4LL) + 

               423LL*Power(r,5LL)*Power(xj,5LL) + 54LL*Power(r,6LL)*Power(xj,6LL) + 

               2LL*Power(r,7LL)*Power(xj,7LL)) - 

            3LL*Power(xi,18LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r,2LL)*Power(xj,2LL) + 

               1050LL*Power(r,3LL)*Power(xj,3LL) + 420LL*Power(r,4LL)*Power(xj,4LL) + 

               126LL*Power(r,5LL)*Power(xj,5LL) + 28LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL)) + 

            3LL*Power(xj,18LL)*(1801800LL + 1576575LL*r*xj + 

               630630LL*Power(r,2LL)*Power(xj,2LL) + 150150LL*Power(r,3LL)*Power(xj,3LL) + 

               23100LL*Power(r,4LL)*Power(xj,4LL) + 2310LL*Power(r,5LL)*Power(xj,5LL) + 

               140LL*Power(r,6LL)*Power(xj,6LL) + 4LL*Power(r,7LL)*Power(xj,7LL)) + 

            2LL*Power(xi,14LL)*Power(xj,4LL)*

             (-147420LL - 257985LL*r*xj - 221130LL*Power(r,2LL)*Power(xj,2LL) - 

               122850LL*Power(r,3LL)*Power(xj,3LL) - 49140LL*Power(r,4LL)*Power(xj,4LL) - 

               17388LL*Power(r,5LL)*Power(xj,5LL) - 1512LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) - 

            42LL*Power(xi,12LL)*Power(xj,6LL)*

             (-25740LL - 45045LL*r*xj - 38610LL*Power(r,2LL)*Power(xj,2LL) - 

               19470LL*Power(r,3LL)*Power(xj,3LL) - 12540LL*Power(r,4LL)*Power(xj,4LL) - 

               1836LL*Power(r,5LL)*Power(xj,5LL) - 8LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) + 

            42LL*Power(xi,6LL)*Power(xj,12LL)*

             (921600LL - 1640835LL*r*xj - 546030LL*Power(r,2LL)*Power(xj,2LL) + 

               20730LL*Power(r,3LL)*Power(xj,3LL) + 30180LL*Power(r,4LL)*Power(xj,4LL) + 

               5028LL*Power(r,5LL)*Power(xj,5LL) + 344LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) - 

            2LL*Power(xi,4LL)*Power(xj,14LL)*

             (-67767840LL - 13377735LL*r*xj + 6601770LL*Power(r,2LL)*Power(xj,2LL) + 

               3115350LL*Power(r,3LL)*Power(xj,3LL) + 

               548940LL*Power(r,4LL)*Power(xj,4LL) + 48132LL*Power(r,5LL)*Power(xj,5LL) + 

               1848LL*Power(r,6LL)*Power(xj,6LL) + 8LL*Power(r,7LL)*Power(xj,7LL)) + 

            Power(xi,16LL)*Power(xj,2LL)*

             (49140LL + 85995LL*r*xj + 73710LL*Power(r,2LL)*Power(xj,2LL) + 

               40950LL*Power(r,3LL)*Power(xj,3LL) + 16380LL*Power(r,4LL)*Power(xj,4LL) + 

               4914LL*Power(r,5LL)*Power(xj,5LL) + 1092LL*Power(r,6LL)*Power(xj,6LL) + 

               44LL*Power(r,7LL)*Power(xj,7LL))))/

       (1890LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,13LL)*Power(xi + xj,12LL)) - 

      (7560LL*exp(2LL*r*(xi + xj))*(xi + xj)*

          Power(Power(xi,2LL) - Power(xj,2LL),13LL) + 

         84LL*exp(2LL*r*xj)*Power(xj,10LL)*

          (-240LL*Power(r,3LL)*Power(xi,20LL) - 10LL*Power(r,4LL)*Power(xi,21LL) - 

            540LL*Power(r,3LL)*Power(xi,18LL)*Power(xj,2LL) - 

            8LL*Power(r,4LL)*Power(xi,19LL)*Power(xj,2LL) + 

            22620LL*Power(r,2LL)*Power(xi,13LL)*Power(xj,6LL) + 75LL*xi*Power(xj,16LL) + 

            120LL*r*Power(xi,2LL)*Power(xj,16LL) + 

            60LL*Power(r,2LL)*Power(xi,3LL)*Power(xj,16LL) - 

            12LL*Power(r,2LL)*Power(xi,19LL)*(195LL + Power(r,2LL)*Power(xj,2LL)) + 

            15LL*Power(xi,3LL)*Power(xj,14LL)*(-65LL + 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            60LL*r*Power(xi,18LL)*(182LL + 9LL*Power(r,2LL)*Power(xj,2LL)) + 

            30LL*Power(xi,13LL)*Power(xj,4LL)*(-13047LL + 377LL*Power(r,2LL)*Power(xj,2LL)) + 

            2LL*r*Power(xi,5LL)*Power(xj,12LL)*

             (-390LL*r*Power(xj,2LL) + 4LL*Power(r,3LL)*Power(xj,4LL)) + 

            10LL*Power(xi,4LL)*Power(xj,12LL)*

             (-156LL*r*Power(xj,2LL) + 4LL*Power(r,3LL)*Power(xj,4LL)) - 

            130LL*Power(xi,6LL)*Power(xj,10LL)*

             (-72LL*r*Power(xj,2LL) + 4LL*Power(r,3LL)*Power(xj,4LL)) + 

            13LL*r*Power(xi,11LL)*Power(xj,6LL)*

             (-3300LL*r*Power(xj,2LL) + 16LL*Power(r,3LL)*Power(xj,4LL)) + 

            r*Power(xi,7LL)*Power(xj,10LL)*

             (6660LL*r*Power(xj,2LL) + 16LL*Power(r,3LL)*Power(xj,4LL)) + 

            210LL*Power(xi,16LL)*(-524LL*r*Power(xj,2LL) + 20LL*Power(r,3LL)*Power(xj,4LL)) - 

            6LL*r*Power(xi,9LL)*Power(xj,8LL)*

             (-1430LL*r*Power(xj,2LL) + 24LL*Power(r,3LL)*Power(xj,4LL)) + 

            3LL*r*Power(xi,17LL)*(-3740LL*r*Power(xj,2LL) + 

               48LL*Power(r,3LL)*Power(xj,4LL)) - 

            30LL*Power(xi,14LL)*Power(xj,2LL)*

             (-24LL*r*Power(xj,2LL) + 52LL*Power(r,3LL)*Power(xj,4LL)) - 

            15LL*Power(xi,8LL)*Power(xj,8LL)*

             (5720LL*r*Power(xj,2LL) + 56LL*Power(r,3LL)*Power(xj,4LL)) + 

            65LL*Power(xi,10LL)*Power(xj,6LL)*

             (-1584LL*r*Power(xj,2LL) + 88LL*Power(r,3LL)*Power(xj,4LL)) - 

            10LL*Power(xi,12LL)*Power(xj,4LL)*

             (-30108LL*r*Power(xj,2LL) + 572LL*Power(r,3LL)*Power(xj,4LL)) + 

            2LL*Power(xi,5LL)*Power(xj,12LL)*

             (2925LL - 195LL*Power(r,2LL)*Power(xj,2LL) + Power(r,4LL)*Power(xj,4LL)) + 

            13LL*Power(xi,11LL)*Power(xj,6LL)*

             (30735LL - 1650LL*Power(r,2LL)*Power(xj,2LL) + 4LL*Power(r,4LL)*Power(xj,4LL)) + 

            Power(xi,7LL)*Power(xj,10LL)*

             (-15015LL + 3330LL*Power(r,2LL)*Power(xj,2LL) + 4LL*Power(r,4LL)*Power(xj,4LL)) 

    - 6LL*Power(xi,9LL)*Power(xj,8LL)*(-48620LL - 715LL*Power(r,2LL)*Power(xj,2LL) + 

               6LL*Power(r,4LL)*Power(xj,4LL)) + 

            3LL*Power(xi,17LL)*(-6825LL - 1870LL*Power(r,2LL)*Power(xj,2LL) + 

               12LL*Power(r,4LL)*Power(xj,4LL)) + 

            Power(xi,15LL)*(-269325LL*Power(xj,2LL) + 27810LL*Power(r,2LL)*Power(xj,4LL) - 

               260LL*Power(r,4LL)*Power(xj,6LL))) + 

         168LL*exp(2LL*r*xj)*Power(xj,11LL)*

          (-60LL*Power(r,4LL)*Power(xi,20LL) - 2LL*Power(r,5LL)*Power(xi,21LL) + 

            45LL*Power(xj,16LL) + 75LL*r*xi*Power(xj,16LL) - 

            4LL*Power(r,3LL)*Power(xi,19LL)*(195LL + Power(r,2LL)*Power(xj,2LL)) + 

            15LL*r*Power(xi,3LL)*Power(xj,14LL)*(-65LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            15LL*Power(xi,2LL)*Power(xj,14LL)*(-39LL + 4LL*Power(r,2LL)*Power(xj,2LL)) - 

            30LL*Power(r,2LL)*Power(xi,18LL)*(182LL + 9LL*Power(r,2LL)*Power(xj,2LL)) + 

            30LL*r*Power(xi,13LL)*Power(xj,4LL)*

             (-13047LL + 377LL*Power(r,2LL)*Power(xj,2LL)) + 

            2LL*r*Power(xi,5LL)*Power(xj,12LL)*

             (2925LL - 195LL*Power(r,2LL)*Power(xj,2LL) + Power(r,4LL)*Power(xj,4LL)) + 

            10LL*Power(xi,4LL)*Power(xj,12LL)*

             (351LL - 78LL*Power(r,2LL)*Power(xj,2LL) + Power(r,4LL)*Power(xj,4LL)) - 

            130LL*Power(xi,6LL)*Power(xj,10LL)*

             (99LL - 36LL*Power(r,2LL)*Power(xj,2LL) + Power(r,4LL)*Power(xj,4LL)) + 

            13LL*r*Power(xi,11LL)*Power(xj,6LL)*

             (30735LL - 1650LL*Power(r,2LL)*Power(xj,2LL) + 4LL*Power(r,4LL)*Power(xj,4LL)) + 

            r*Power(xi,7LL)*Power(xj,10LL)*

             (-15015LL + 3330LL*Power(r,2LL)*Power(xj,2LL) + 4LL*Power(r,4LL)*Power(xj,4LL)) 

    + 210LL*Power(xi,16LL)*(-156LL - 262LL*Power(r,2LL)*Power(xj,2LL) + 

               5LL*Power(r,4LL)*Power(xj,4LL)) - 

            6LL*r*Power(xi,9LL)*Power(xj,8LL)*

             (-48620LL - 715LL*Power(r,2LL)*Power(xj,2LL) + 6LL*Power(r,4LL)*Power(xj,4LL)) + 

            3LL*r*Power(xi,17LL)*(-6825LL - 1870LL*Power(r,2LL)*Power(xj,2LL) + 

               12LL*Power(r,4LL)*Power(xj,4LL)) - 

            30LL*Power(xi,14LL)*Power(xj,2LL)*

             (17934LL - 12LL*Power(r,2LL)*Power(xj,2LL) + 13LL*Power(r,4LL)*Power(xj,4LL)) - 

            15LL*Power(xi,8LL)*Power(xj,8LL)*

             (2145LL + 2860LL*Power(r,2LL)*Power(xj,2LL) + 14LL*Power(r,4LL)*Power(xj,4LL)) + 

            65LL*Power(xi,10LL)*Power(xj,6LL)*

             (-13725LL - 792LL*Power(r,2LL)*Power(xj,2LL) + 22LL*Power(r,4LL)*Power(xj,4LL)) 

    - 10LL*Power(xi,12LL)*Power(xj,4LL)*(153630LL - 15054LL*Power(r,2LL)*Power(xj,2LL) + 

               143LL*Power(r,4LL)*Power(xj,4LL)) + 

            Power(xi,15LL)*(-269325LL*r*Power(xj,2LL) + 9270LL*Power(r,3LL)*Power(xj,4LL) - 

               52LL*Power(r,5LL)*Power(xj,6LL))) + 

         exp(2LL*r*xi)*Power(xi,8LL)*

          (Power(xi,2LL)*Power(xj,16LL)*

             (47669895LL*xj + 27862380LL*r*Power(xj,2LL) + 

               6511050LL*Power(r,2LL)*Power(xj,3LL) + 

               677040LL*Power(r,3LL)*Power(xj,4LL) + 8190LL*Power(r,4LL)*Power(xj,5LL) - 

               4536LL*Power(r,5LL)*Power(xj,6LL) - 308LL*Power(r,6LL)*Power(xj,7LL)) + 

            364LL*Power(xi,10LL)*Power(xj,8LL)*

             (-13860LL*xj - 11880LL*r*Power(xj,2LL) - 35640LL*Power(r,2LL)*Power(xj,3LL) - 

               10560LL*Power(r,3LL)*Power(xj,4LL) - 225LL*Power(r,4LL)*Power(xj,5LL) + 

               180LL*Power(r,5LL)*Power(xj,6LL) + 14LL*Power(r,6LL)*Power(xj,7LL)) - 

            364LL*Power(xi,8LL)*Power(xj,10LL)*

             (18270LL*xj - 116640LL*r*Power(xj,2LL) - 53190LL*Power(r,2LL)*Power(xj,3LL) - 

               1200LL*Power(r,3LL)*Power(xj,4LL) + 2115LL*Power(r,4LL)*Power(xj,5LL) + 

               324LL*Power(r,5LL)*Power(xj,6LL) + 14LL*Power(r,6LL)*Power(xj,7LL)) - 

            3LL*Power(xi,18LL)*(2205LL*xj + 3780LL*r*Power(xj,2LL) + 

               3150LL*Power(r,2LL)*Power(xj,3LL) + 1680LL*Power(r,3LL)*Power(xj,4LL) + 

               630LL*Power(r,4LL)*Power(xj,5LL) + 168LL*Power(r,5LL)*Power(xj,6LL) + 

               28LL*Power(r,6LL)*Power(xj,7LL)) + 

            3LL*Power(xj,18LL)*(1576575LL*xj + 1261260LL*r*Power(xj,2LL) + 

               450450LL*Power(r,2LL)*Power(xj,3LL) + 92400LL*Power(r,3LL)*Power(xj,4LL) + 

               11550LL*Power(r,4LL)*Power(xj,5LL) + 840LL*Power(r,5LL)*Power(xj,6LL) + 

               28LL*Power(r,6LL)*Power(xj,7LL)) + 

            2LL*Power(xi,14LL)*Power(xj,4LL)*

             (-257985LL*xj - 442260LL*r*Power(xj,2LL) - 

               368550LL*Power(r,2LL)*Power(xj,3LL) - 196560LL*Power(r,3LL)*Power(xj,4LL) - 

               86940LL*Power(r,4LL)*Power(xj,5LL) - 9072LL*Power(r,5LL)*Power(xj,6LL) + 

               56LL*Power(r,6LL)*Power(xj,7LL)) - 

            42LL*Power(xi,12LL)*Power(xj,6LL)*

             (-45045LL*xj - 77220LL*r*Power(xj,2LL) - 58410LL*Power(r,2LL)*Power(xj,3LL) - 

               50160LL*Power(r,3LL)*Power(xj,4LL) - 9180LL*Power(r,4LL)*Power(xj,5LL) - 

               48LL*Power(r,5LL)*Power(xj,6LL) + 56LL*Power(r,6LL)*Power(xj,7LL)) + 

            42LL*Power(xi,6LL)*Power(xj,12LL)*

             (-1640835LL*xj - 1092060LL*r*Power(xj,2LL) + 

               62190LL*Power(r,2LL)*Power(xj,3LL) + 120720LL*Power(r,3LL)*Power(xj,4LL) + 

               25140LL*Power(r,4LL)*Power(xj,5LL) + 2064LL*Power(r,5LL)*Power(xj,6LL) + 

               56LL*Power(r,6LL)*Power(xj,7LL)) - 

            2LL*Power(xi,4LL)*Power(xj,14LL)*

             (-13377735LL*xj + 13203540LL*r*Power(xj,2LL) + 

               9346050LL*Power(r,2LL)*Power(xj,3LL) + 

               2195760LL*Power(r,3LL)*Power(xj,4LL) + 

               240660LL*Power(r,4LL)*Power(xj,5LL) + 11088LL*Power(r,5LL)*Power(xj,6LL) + 

               56LL*Power(r,6LL)*Power(xj,7LL)) + 

            Power(xi,16LL)*Power(xj,2LL)*

             (85995LL*xj + 147420LL*r*Power(xj,2LL) + 122850LL*Power(r,2LL)*Power(xj,3LL) + 

               65520LL*Power(r,3LL)*Power(xj,4LL) + 24570LL*Power(r,4LL)*Power(xj,5LL) + 

               6552LL*Power(r,5LL)*Power(xj,6LL) + 308LL*Power(r,6LL)*Power(xj,7LL))) + 

         2LL*exp(2LL*r*xi)*Power(xi,9LL)*

          (Power(xi,2LL)*Power(xj,16LL)*

             (70073640LL + 47669895LL*r*xj + 13931190LL*Power(r,2LL)*Power(xj,2LL) + 

               2170350LL*Power(r,3LL)*Power(xj,3LL) + 169260LL*Power(r,4LL)*Power(xj,4LL) + 

               1638LL*Power(r,5LL)*Power(xj,5LL) - 756LL*Power(r,6LL)*Power(xj,6LL) - 

               44LL*Power(r,7LL)*Power(xj,7LL)) + 

            364LL*Power(xi,10LL)*Power(xj,8LL)*

             (-7425LL - 13860LL*r*xj - 5940LL*Power(r,2LL)*Power(xj,2LL) - 

               11880LL*Power(r,3LL)*Power(xj,3LL) - 2640LL*Power(r,4LL)*Power(xj,4LL) - 

               45LL*Power(r,5LL)*Power(xj,5LL) + 30LL*Power(r,6LL)*Power(xj,6LL) + 

               2LL*Power(r,7LL)*Power(xj,7LL)) - 

            364LL*Power(xi,8LL)*Power(xj,10LL)*

             (-20925LL + 18270LL*r*xj - 58320LL*Power(r,2LL)*Power(xj,2LL) - 

               17730LL*Power(r,3LL)*Power(xj,3LL) - 300LL*Power(r,4LL)*Power(xj,4LL) + 

               423LL*Power(r,5LL)*Power(xj,5LL) + 54LL*Power(r,6LL)*Power(xj,6LL) + 

               2LL*Power(r,7LL)*Power(xj,7LL)) - 

            3LL*Power(xi,18LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r,2LL)*Power(xj,2LL) + 

               1050LL*Power(r,3LL)*Power(xj,3LL) + 420LL*Power(r,4LL)*Power(xj,4LL) + 

               126LL*Power(r,5LL)*Power(xj,5LL) + 28LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL)) + 

            3LL*Power(xj,18LL)*(1801800LL + 1576575LL*r*xj + 

               630630LL*Power(r,2LL)*Power(xj,2LL) + 150150LL*Power(r,3LL)*Power(xj,3LL) + 

               23100LL*Power(r,4LL)*Power(xj,4LL) + 2310LL*Power(r,5LL)*Power(xj,5LL) + 

               140LL*Power(r,6LL)*Power(xj,6LL) + 4LL*Power(r,7LL)*Power(xj,7LL)) + 

            2LL*Power(xi,14LL)*Power(xj,4LL)*

             (-147420LL - 257985LL*r*xj - 221130LL*Power(r,2LL)*Power(xj,2LL) - 

               122850LL*Power(r,3LL)*Power(xj,3LL) - 49140LL*Power(r,4LL)*Power(xj,4LL) - 

               17388LL*Power(r,5LL)*Power(xj,5LL) - 1512LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) - 

            42LL*Power(xi,12LL)*Power(xj,6LL)*

             (-25740LL - 45045LL*r*xj - 38610LL*Power(r,2LL)*Power(xj,2LL) - 

               19470LL*Power(r,3LL)*Power(xj,3LL) - 12540LL*Power(r,4LL)*Power(xj,4LL) - 

               1836LL*Power(r,5LL)*Power(xj,5LL) - 8LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) + 

            42LL*Power(xi,6LL)*Power(xj,12LL)*

             (921600LL - 1640835LL*r*xj - 546030LL*Power(r,2LL)*Power(xj,2LL) + 

               20730LL*Power(r,3LL)*Power(xj,3LL) + 30180LL*Power(r,4LL)*Power(xj,4LL) + 

               5028LL*Power(r,5LL)*Power(xj,5LL) + 344LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) - 

            2LL*Power(xi,4LL)*Power(xj,14LL)*

             (-67767840LL - 13377735LL*r*xj + 6601770LL*Power(r,2LL)*Power(xj,2LL) + 

               3115350LL*Power(r,3LL)*Power(xj,3LL) + 548940LL*Power(r,4LL)*Power(xj,4LL) + 

               48132LL*Power(r,5LL)*Power(xj,5LL) + 1848LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) + 

            Power(xi,16LL)*Power(xj,2LL)*

             (49140LL + 85995LL*r*xj + 73710LL*Power(r,2LL)*Power(xj,2LL) + 

               40950LL*Power(r,3LL)*Power(xj,3LL) + 16380LL*Power(r,4LL)*Power(xj,4LL) + 

               4914LL*Power(r,5LL)*Power(xj,5LL) + 1092LL*Power(r,6LL)*Power(xj,6LL) + 

               44LL*Power(r,7LL)*Power(xj,7LL))))/

       (3780LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,13LL)*Power(xi + xj,13LL))

    ; }
   
  }
  return S;
}

cl_R DSlater_3S_5S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = 0LL

    ; } else {  S = -(-568188982486125LL*xi + 627683696640000LL*exp(2LL*r*xi)*xi - 

          1017388536664500LL*r*Power(xi,2LL) - 

          899677411132500LL*Power(r,2LL)*Power(xi,3LL) - 

          523015260768000LL*Power(r,3LL)*Power(xi,4LL) - 

          224405775594000LL*Power(r,4LL)*Power(xi,5LL) - 

          75610821686400LL*Power(r,5LL)*Power(xi,6LL) - 

          20775676521600LL*Power(r,6LL)*Power(xi,7LL) - 

          4769897932800LL*Power(r,7LL)*Power(xi,8LL) - 

          929382854400LL*Power(r,8LL)*Power(xi,9LL) - 

          154983628800LL*Power(r,9LL)*Power(xi,10LL) - 

          22140518400LL*Power(r,10LL)*Power(xi,11LL) - 

          2683699200LL*Power(r,11LL)*Power(xi,12LL) - 

          268369920LL*Power(r,12LL)*Power(xi,13LL) - 

          20643840LL*Power(r,13LL)*Power(xi,14LL) - 983040LL*Power(r,14LL)*Power(xi,15LL))/

       (3.1384184832e14*exp(2LL*r*xi)*r) + 

      (-313841848320000LL + 313841848320000LL*exp(2LL*r*xi) - 

         568188982486125LL*r*xi - 508694268332250LL*Power(r,2LL)*Power(xi,2LL) - 

         299892470377500LL*Power(r,3LL)*Power(xi,3LL) - 

         130753815192000LL*Power(r,4LL)*Power(xi,4LL) - 

         44881155118800LL*Power(r,5LL)*Power(xi,5LL) - 

         12601803614400LL*Power(r,6LL)*Power(xi,6LL) - 

         2967953788800LL*Power(r,7LL)*Power(xi,7LL) - 

         596237241600LL*Power(r,8LL)*Power(xi,8LL) - 

         103264761600LL*Power(r,9LL)*Power(xi,9LL) - 

         15498362880LL*Power(r,10LL)*Power(xi,10LL) - 

         2012774400LL*Power(r,11LL)*Power(xi,11LL) - 

         223641600LL*Power(r,12LL)*Power(xi,12LL) - 

         20643840LL*Power(r,13LL)*Power(xi,13LL) - 1474560LL*Power(r,14LL)*Power(xi,14LL) - 

         65536LL*Power(r,15LL)*Power(xi,15LL))/

       (3.1384184832e14*exp(2LL*r*xi)*Power(r,2LL)) + 

      (xi*(-313841848320000LL + 313841848320000LL*exp(2LL*r*xi) - 

           568188982486125LL*r*xi - 508694268332250LL*Power(r,2LL)*Power(xi,2LL) - 

           299892470377500LL*Power(r,3LL)*Power(xi,3LL) - 

           130753815192000LL*Power(r,4LL)*Power(xi,4LL) - 

           44881155118800LL*Power(r,5LL)*Power(xi,5LL) - 

           12601803614400LL*Power(r,6LL)*Power(xi,6LL) - 

           2967953788800LL*Power(r,7LL)*Power(xi,7LL) - 

           596237241600LL*Power(r,8LL)*Power(xi,8LL) - 

           103264761600LL*Power(r,9LL)*Power(xi,9LL) - 

           15498362880LL*Power(r,10LL)*Power(xi,10LL) - 

           2012774400LL*Power(r,11LL)*Power(xi,11LL) - 

           223641600LL*Power(r,12LL)*Power(xi,12LL) - 

           20643840LL*Power(r,13LL)*Power(xi,13LL) - 

           1474560LL*Power(r,14LL)*Power(xi,14LL) - 65536LL*Power(r,15LL)*Power(xi,15LL)))/

       (1.5692092416e14*exp(2LL*r*xi)*r)

    ; }
 
  }
  else {
      if (r == 0LL) {  S = 0LL

    ; } else {  S = (42525LL*exp(2LL*r*(xi + xj))*Power(Power(xi,2LL) - Power(xj,2LL),15LL) + 

         189LL*exp(2LL*r*xj)*Power(xj,12LL)*

          (-350LL*Power(r,4LL)*Power(xi,22LL) - 10LL*Power(r,5LL)*Power(xi,23LL) + 

            225LL*Power(xj,18LL) + 375LL*r*xi*Power(xj,18LL) - 

            70LL*Power(r,3LL)*Power(xi,21LL)*(75LL + Power(r,2LL)*Power(xj,2LL)) + 

            75LL*r*Power(xi,3LL)*Power(xj,16LL)*(-75LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            75LL*Power(xi,2LL)*Power(xj,16LL)*(-45LL + 4LL*Power(r,2LL)*Power(xj,2LL)) - 

            50LL*Power(r,2LL)*Power(xi,20LL)*(840LL + 71LL*Power(r,2LL)*Power(xj,2LL)) + 

            r*Power(xi,9LL)*Power(xj,10LL)*

             (4694625LL + 124800LL*Power(r,2LL)*Power(xj,2LL) - 

               248LL*Power(r,4LL)*Power(xj,4LL)) + 

            20LL*r*Power(xi,17LL)*Power(xj,2LL)*

             (-185895LL - 948LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) 

    + 5LL*r*Power(xi,5LL)*Power(xj,14LL)*

             (7875LL - 450LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            25LL*Power(xi,4LL)*Power(xj,14LL)*

             (945LL - 180LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) - 

            375LL*Power(xi,6LL)*Power(xj,12LL)*

             (273LL - 84LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) - 

            5LL*r*Power(xi,11LL)*Power(xj,8LL)*

             (-2803125LL + 49140LL*Power(r,2LL)*Power(xj,2LL) + 

               8LL*Power(r,4LL)*Power(xj,4LL)) + 

            5LL*r*Power(xi,7LL)*Power(xj,12LL)*

             (-16965LL + 5152LL*Power(r,2LL)*Power(xj,2LL) + 

               14LL*Power(r,4LL)*Power(xj,4LL)) + 

            325LL*Power(xi,10LL)*Power(xj,8LL)*

             (-60117LL - 5340LL*Power(r,2LL)*Power(xj,2LL) + 

               40LL*Power(r,4LL)*Power(xj,4LL)) - 

            15LL*r*Power(xi,15LL)*Power(xj,4LL)*

             (845085LL - 22960LL*Power(r,2LL)*Power(xj,2LL) + 

               52LL*Power(r,4LL)*Power(xj,4LL)) + 

            15LL*r*Power(xi,13LL)*Power(xj,6LL)*

             (-139125LL - 10140LL*Power(r,2LL)*Power(xj,2LL) + 

               52LL*Power(r,4LL)*Power(xj,4LL)) + 

            75LL*Power(xi,12LL)*Power(xj,6LL)*

             (-729687LL + 25532LL*Power(r,2LL)*Power(xj,2LL) + 

               52LL*Power(r,4LL)*Power(xj,4LL)) + 

            60LL*Power(xi,18LL)*(-5355LL - 11940LL*Power(r,2LL)*Power(xj,2LL) + 

               86LL*Power(r,4LL)*Power(xj,4LL)) + 

            2LL*r*Power(xi,19LL)*(-89250LL - 35425LL*Power(r,2LL)*Power(xj,2LL) + 

               124LL*Power(r,4LL)*Power(xj,4LL)) + 

            100LL*Power(xi,16LL)*Power(xj,2LL)*

             (-79713LL - 13311LL*Power(r,2LL)*Power(xj,2LL) + 

               146LL*Power(r,4LL)*Power(xj,4LL)) - 

            5LL*Power(xi,8LL)*Power(xj,10LL)*

             (157365LL + 95940LL*Power(r,2LL)*Power(xj,2LL) + 

               952LL*Power(r,4LL)*Power(xj,4LL)) - 

            15LL*Power(xi,14LL)*Power(xj,4LL)*

             (2638467LL - 157500LL*Power(r,2LL)*Power(xj,2LL) + 

               1820LL*Power(r,4LL)*Power(xj,4LL))) + 

         exp(2LL*r*xi)*Power(xi,8LL)*

          (2LL*Power(xi,2LL)*Power(xj,20LL)*

             (1782492075LL + 1449175455LL*r*xj + 

               533365560LL*Power(r,2LL)*Power(xj,2LL) + 

               114631335LL*Power(r,3LL)*Power(xj,3LL) + 

               15221115LL*Power(r,4LL)*Power(xj,4LL) + 

               1142505LL*Power(r,5LL)*Power(xj,5LL) + 18396LL*Power(r,6LL)*Power(xj,6LL) - 

               5238LL*Power(r,7LL)*Power(xj,7LL) - 513LL*Power(r,8LL)*Power(xj,8LL) - 

               17LL*Power(r,9LL)*Power(xj,9LL)) + 

            42LL*Power(xi,4LL)*Power(xj,18LL)*

             (251336925LL + 104824125LL*r*xj + 340200LL*Power(r,2LL)*Power(xj,2LL) - 

               9122085LL*Power(r,3LL)*Power(xj,3LL) - 

               2798145LL*Power(r,4LL)*Power(xj,4LL) - 

               433755LL*Power(r,5LL)*Power(xj,5LL) - 39060LL*Power(r,6LL)*Power(xj,6LL) - 

               1890LL*Power(r,7LL)*Power(xj,7LL) - 27LL*Power(r,8LL)*Power(xj,8LL) + 

               Power(r,9LL)*Power(xj,9LL)) + 

            6LL*Power(xj,22LL)*(34459425LL + 34459425LL*r*xj + 

               16216200LL*Power(r,2LL)*Power(xj,2LL) + 

               4729725LL*Power(r,3LL)*Power(xj,3LL) + 

               945945LL*Power(r,4LL)*Power(xj,4LL) + 135135LL*Power(r,5LL)*Power(xj,5LL) + 

               13860LL*Power(r,6LL)*Power(xj,6LL) + 990LL*Power(r,7LL)*Power(xj,7LL) + 

               45LL*Power(r,8LL)*Power(xj,8LL) + Power(r,9LL)*Power(xj,9LL)) - 

            3LL*Power(xi,22LL)*(14175LL + 25515LL*r*xj + 

               22680LL*Power(r,2LL)*Power(xj,2LL) + 13230LL*Power(r,3LL)*Power(xj,3LL) + 

               5670LL*Power(r,4LL)*Power(xj,4LL) + 1890LL*Power(r,5LL)*Power(xj,5LL) + 

               504LL*Power(r,6LL)*Power(xj,6LL) + 108LL*Power(r,7LL)*Power(xj,7LL) + 

               18LL*Power(r,8LL)*Power(xj,8LL) + 2LL*Power(r,9LL)*Power(xj,9LL)) - 

            21LL*Power(xi,18LL)*Power(xj,4LL)*

             (212625LL + 382725LL*r*xj + 340200LL*Power(r,2LL)*Power(xj,2LL) + 

               198450LL*Power(r,3LL)*Power(xj,3LL) + 85050LL*Power(r,4LL)*Power(xj,4LL) + 

               28350LL*Power(r,5LL)*Power(xj,5LL) + 7560LL*Power(r,6LL)*Power(xj,6LL) + 

               1836LL*Power(r,7LL)*Power(xj,7LL) + 162LL*Power(r,8LL)*Power(xj,8LL) + 

               2LL*Power(r,9LL)*Power(xj,9LL)) + 

            54LL*Power(xi,6LL)*Power(xj,16LL)*

             (133451955LL - 73700865LL*r*xj - 54096840LL*Power(r,2LL)*Power(xj,2LL) - 

               8306235LL*Power(r,3LL)*Power(xj,3LL) + 

               966945LL*Power(r,4LL)*Power(xj,4LL) + 516747LL*Power(r,5LL)*Power(xj,5LL) + 

               80724LL*Power(r,6LL)*Power(xj,6LL) + 6434LL*Power(r,7LL)*Power(xj,7LL) + 

               251LL*Power(r,8LL)*Power(xj,8LL) + 3LL*Power(r,9LL)*Power(xj,9LL)) - 

            315LL*Power(xi,12LL)*Power(xj,10LL)*

             (-405405LL - 710073LL*r*xj - 805896LL*Power(r,2LL)*Power(xj,2LL) - 

               101556LL*Power(r,3LL)*Power(xj,3LL) - 258804LL*Power(r,4LL)*Power(xj,4LL) - 

               90972LL*Power(r,5LL)*Power(xj,5LL) - 9744LL*Power(r,6LL)*Power(xj,6LL) + 

               120LL*Power(r,7LL)*Power(xj,7LL) + 84LL*Power(r,8LL)*Power(xj,8LL) + 

               4LL*Power(r,9LL)*Power(xj,9LL)) + 

            315LL*Power(xi,10LL)*Power(xj,12LL)*

             (-482895LL - 2656395LL*r*xj + 1186920LL*Power(r,2LL)*Power(xj,2LL) - 

               1155420LL*Power(r,3LL)*Power(xj,3LL) - 

               643356LL*Power(r,4LL)*Power(xj,4LL) - 93492LL*Power(r,5LL)*Power(xj,5LL) + 

               336LL*Power(r,6LL)*Power(xj,6LL) + 1368LL*Power(r,7LL)*Power(xj,7LL) + 

               132LL*Power(r,8LL)*Power(xj,8LL) + 4LL*Power(r,9LL)*Power(xj,9LL)) - 

            27LL*Power(xi,16LL)*Power(xj,6LL)*

             (-716625LL - 1289925LL*r*xj - 1146600LL*Power(r,2LL)*Power(xj,2LL) - 

               668850LL*Power(r,3LL)*Power(xj,3LL) - 286650LL*Power(r,4LL)*Power(xj,4LL) - 

               90006LL*Power(r,5LL)*Power(xj,5LL) - 32872LL*Power(r,6LL)*Power(xj,6LL) - 

               4812LL*Power(r,7LL)*Power(xj,7LL) - 178LL*Power(r,8LL)*Power(xj,8LL) + 

               6LL*Power(r,9LL)*Power(xj,9LL)) + 

            Power(xi,20LL)*Power(xj,2LL)*

             (637875LL + 1148175LL*r*xj + 1020600LL*Power(r,2LL)*Power(xj,2LL) + 

               595350LL*Power(r,3LL)*Power(xj,3LL) + 255150LL*Power(r,4LL)*Power(xj,4LL) + 

               85050LL*Power(r,5LL)*Power(xj,5LL) + 22680LL*Power(r,6LL)*Power(xj,6LL) + 

               4860LL*Power(r,7LL)*Power(xj,7LL) + 810LL*Power(r,8LL)*Power(xj,8LL) + 

               34LL*Power(r,9LL)*Power(xj,9LL)) + 

            3LL*Power(xi,14LL)*Power(xj,8LL)*

             (-19348875LL - 34827975LL*r*xj - 30958200LL*Power(r,2LL)*Power(xj,2LL) - 

               18689580LL*Power(r,3LL)*Power(xj,3LL) - 

               5847660LL*Power(r,4LL)*Power(xj,4LL) - 

               3723300LL*Power(r,5LL)*Power(xj,5LL) - 

               845040LL*Power(r,6LL)*Power(xj,6LL) - 58680LL*Power(r,7LL)*Power(xj,7LL) + 

               1548LL*Power(r,8LL)*Power(xj,8LL) + 236LL*Power(r,9LL)*Power(xj,9LL)) - 

            3LL*Power(xi,8LL)*Power(xj,14LL)*

             (-593408025LL + 946053675LL*r*xj - 394427880LL*Power(r,2LL)*Power(xj,2LL) - 

               315870660LL*Power(r,3LL)*Power(xj,3LL) - 

               53891460LL*Power(r,4LL)*Power(xj,4LL) + 

               910980LL*Power(r,5LL)*Power(xj,5LL) + 1409520LL*Power(r,6LL)*Power(xj,6LL) + 

               192168LL*Power(r,7LL)*Power(xj,7LL) + 11196LL*Power(r,8LL)*Power(xj,8LL) + 

               236LL*Power(r,9LL)*Power(xj,9LL))))/

       (42525LL*exp(2LL*r*(xi + xj))*Power(r,2LL)*Power(xi - xj,15LL)*

         Power(xi + xj,15LL)) + (2LL*(42525LL*exp(2LL*r*(xi + xj))*

            Power(Power(xi,2LL) - Power(xj,2LL),15LL) + 

           189LL*exp(2LL*r*xj)*Power(xj,12LL)*

            (-350LL*Power(r,4LL)*Power(xi,22LL) - 10LL*Power(r,5LL)*Power(xi,23LL) + 

              225LL*Power(xj,18LL) + 375LL*r*xi*Power(xj,18LL) - 

              70LL*Power(r,3LL)*Power(xi,21LL)*(75LL + Power(r,2LL)*Power(xj,2LL)) + 

              75LL*r*Power(xi,3LL)*Power(xj,16LL)*(-75LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

              75LL*Power(xi,2LL)*Power(xj,16LL)*(-45LL + 4LL*Power(r,2LL)*Power(xj,2LL)) - 

              50LL*Power(r,2LL)*Power(xi,20LL)*(840LL + 71LL*Power(r,2LL)*Power(xj,2LL)) + 

              r*Power(xi,9LL)*Power(xj,10LL)*

               (4694625LL + 124800LL*Power(r,2LL)*Power(xj,2LL) - 

                 248LL*Power(r,4LL)*Power(xj,4LL)) + 

              20LL*r*Power(xi,17LL)*Power(xj,2LL)*

               (-185895LL - 948LL*Power(r,2LL)*Power(xj,2LL) + 

                 2LL*Power(r,4LL)*Power(xj,4LL)) + 

              5LL*r*Power(xi,5LL)*Power(xj,14LL)*

               (7875LL - 450LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) 

    + 25LL*Power(xi,4LL)*Power(xj,14LL)*

               (945LL - 180LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) - 

              375LL*Power(xi,6LL)*Power(xj,12LL)*

               (273LL - 84LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) - 

              5LL*r*Power(xi,11LL)*Power(xj,8LL)*

               (-2803125LL + 49140LL*Power(r,2LL)*Power(xj,2LL) + 

                 8LL*Power(r,4LL)*Power(xj,4LL)) + 

              5LL*r*Power(xi,7LL)*Power(xj,12LL)*

               (-16965LL + 5152LL*Power(r,2LL)*Power(xj,2LL) + 

                 14LL*Power(r,4LL)*Power(xj,4LL)) + 

              325LL*Power(xi,10LL)*Power(xj,8LL)*

               (-60117LL - 5340LL*Power(r,2LL)*Power(xj,2LL) + 

                 40LL*Power(r,4LL)*Power(xj,4LL)) - 

              15LL*r*Power(xi,15LL)*Power(xj,4LL)*

               (845085LL - 22960LL*Power(r,2LL)*Power(xj,2LL) + 

                 52LL*Power(r,4LL)*Power(xj,4LL)) + 

              15LL*r*Power(xi,13LL)*Power(xj,6LL)*

               (-139125LL - 10140LL*Power(r,2LL)*Power(xj,2LL) + 

                 52LL*Power(r,4LL)*Power(xj,4LL)) + 

              75LL*Power(xi,12LL)*Power(xj,6LL)*

               (-729687LL + 25532LL*Power(r,2LL)*Power(xj,2LL) + 

                 52LL*Power(r,4LL)*Power(xj,4LL)) + 

              60LL*Power(xi,18LL)*(-5355LL - 11940LL*Power(r,2LL)*Power(xj,2LL) + 

                 86LL*Power(r,4LL)*Power(xj,4LL)) + 

              2LL*r*Power(xi,19LL)*(-89250LL - 35425LL*Power(r,2LL)*Power(xj,2LL) + 

                 124LL*Power(r,4LL)*Power(xj,4LL)) + 

              100LL*Power(xi,16LL)*Power(xj,2LL)*

               (-79713LL - 13311LL*Power(r,2LL)*Power(xj,2LL) + 

                 146LL*Power(r,4LL)*Power(xj,4LL)) - 

              5LL*Power(xi,8LL)*Power(xj,10LL)*

               (157365LL + 95940LL*Power(r,2LL)*Power(xj,2LL) + 

                 952LL*Power(r,4LL)*Power(xj,4LL)) - 

              15LL*Power(xi,14LL)*Power(xj,4LL)*

               (2638467LL - 157500LL*Power(r,2LL)*Power(xj,2LL) + 

                 1820LL*Power(r,4LL)*Power(xj,4LL))) + 

           exp(2LL*r*xi)*Power(xi,8LL)*

            (2LL*Power(xi,2LL)*Power(xj,20LL)*

               (1782492075LL + 1449175455LL*r*xj + 

                 533365560LL*Power(r,2LL)*Power(xj,2LL) + 

                 114631335LL*Power(r,3LL)*Power(xj,3LL) + 

                 15221115LL*Power(r,4LL)*Power(xj,4LL) + 

                 1142505LL*Power(r,5LL)*Power(xj,5LL) + 

                 18396LL*Power(r,6LL)*Power(xj,6LL) - 5238LL*Power(r,7LL)*Power(xj,7LL) - 

                 513LL*Power(r,8LL)*Power(xj,8LL) - 17LL*Power(r,9LL)*Power(xj,9LL)) + 

              42LL*Power(xi,4LL)*Power(xj,18LL)*

               (251336925LL + 104824125LL*r*xj + 340200LL*Power(r,2LL)*Power(xj,2LL) - 

                 9122085LL*Power(r,3LL)*Power(xj,3LL) - 

                 2798145LL*Power(r,4LL)*Power(xj,4LL) - 

                 433755LL*Power(r,5LL)*Power(xj,5LL) - 

                 39060LL*Power(r,6LL)*Power(xj,6LL) - 1890LL*Power(r,7LL)*Power(xj,7LL) - 

                 27LL*Power(r,8LL)*Power(xj,8LL) + Power(r,9LL)*Power(xj,9LL)) + 

              6LL*Power(xj,22LL)*(34459425LL + 34459425LL*r*xj + 

                 16216200LL*Power(r,2LL)*Power(xj,2LL) + 

                 4729725LL*Power(r,3LL)*Power(xj,3LL) + 

                 945945LL*Power(r,4LL)*Power(xj,4LL) + 

                 135135LL*Power(r,5LL)*Power(xj,5LL) + 

                 13860LL*Power(r,6LL)*Power(xj,6LL) + 990LL*Power(r,7LL)*Power(xj,7LL) + 

                 45LL*Power(r,8LL)*Power(xj,8LL) + Power(r,9LL)*Power(xj,9LL)) - 

              3LL*Power(xi,22LL)*(14175LL + 25515LL*r*xj + 

                 22680LL*Power(r,2LL)*Power(xj,2LL) + 13230LL*Power(r,3LL)*Power(xj,3LL) + 

                 5670LL*Power(r,4LL)*Power(xj,4LL) + 1890LL*Power(r,5LL)*Power(xj,5LL) + 

                 504LL*Power(r,6LL)*Power(xj,6LL) + 108LL*Power(r,7LL)*Power(xj,7LL) + 

                 18LL*Power(r,8LL)*Power(xj,8LL) + 2LL*Power(r,9LL)*Power(xj,9LL)) - 

              21LL*Power(xi,18LL)*Power(xj,4LL)*

               (212625LL + 382725LL*r*xj + 340200LL*Power(r,2LL)*Power(xj,2LL) + 

                 198450LL*Power(r,3LL)*Power(xj,3LL) + 

                 85050LL*Power(r,4LL)*Power(xj,4LL) + 28350LL*Power(r,5LL)*Power(xj,5LL) + 

                 7560LL*Power(r,6LL)*Power(xj,6LL) + 1836LL*Power(r,7LL)*Power(xj,7LL) + 

                 162LL*Power(r,8LL)*Power(xj,8LL) + 2LL*Power(r,9LL)*Power(xj,9LL)) + 

              54LL*Power(xi,6LL)*Power(xj,16LL)*

               (133451955LL - 73700865LL*r*xj - 54096840LL*Power(r,2LL)*Power(xj,2LL) - 

                 8306235LL*Power(r,3LL)*Power(xj,3LL) + 

                 966945LL*Power(r,4LL)*Power(xj,4LL) + 

                 516747LL*Power(r,5LL)*Power(xj,5LL) + 

                 80724LL*Power(r,6LL)*Power(xj,6LL) + 6434LL*Power(r,7LL)*Power(xj,7LL) + 

                 251LL*Power(r,8LL)*Power(xj,8LL) + 3LL*Power(r,9LL)*Power(xj,9LL)) - 

              315LL*Power(xi,12LL)*Power(xj,10LL)*

               (-405405LL - 710073LL*r*xj - 805896LL*Power(r,2LL)*Power(xj,2LL) - 

                 101556LL*Power(r,3LL)*Power(xj,3LL) - 

                 258804LL*Power(r,4LL)*Power(xj,4LL) - 

                 90972LL*Power(r,5LL)*Power(xj,5LL) - 9744LL*Power(r,6LL)*Power(xj,6LL) + 

                 120LL*Power(r,7LL)*Power(xj,7LL) + 84LL*Power(r,8LL)*Power(xj,8LL) + 

                 4LL*Power(r,9LL)*Power(xj,9LL)) + 

              315LL*Power(xi,10LL)*Power(xj,12LL)*

               (-482895LL - 2656395LL*r*xj + 1186920LL*Power(r,2LL)*Power(xj,2LL) - 

                 1155420LL*Power(r,3LL)*Power(xj,3LL) - 

                 643356LL*Power(r,4LL)*Power(xj,4LL) - 

                 93492LL*Power(r,5LL)*Power(xj,5LL) + 336LL*Power(r,6LL)*Power(xj,6LL) + 

                 1368LL*Power(r,7LL)*Power(xj,7LL) + 132LL*Power(r,8LL)*Power(xj,8LL) + 

                 4LL*Power(r,9LL)*Power(xj,9LL)) - 

              27LL*Power(xi,16LL)*Power(xj,6LL)*

               (-716625LL - 1289925LL*r*xj - 1146600LL*Power(r,2LL)*Power(xj,2LL) - 

                 668850LL*Power(r,3LL)*Power(xj,3LL) - 

                 286650LL*Power(r,4LL)*Power(xj,4LL) - 

                 90006LL*Power(r,5LL)*Power(xj,5LL) - 32872LL*Power(r,6LL)*Power(xj,6LL) - 

                 4812LL*Power(r,7LL)*Power(xj,7LL) - 178LL*Power(r,8LL)*Power(xj,8LL) + 

                 6LL*Power(r,9LL)*Power(xj,9LL)) + 

              Power(xi,20LL)*Power(xj,2LL)*

               (637875LL + 1148175LL*r*xj + 1020600LL*Power(r,2LL)*Power(xj,2LL) + 

                 595350LL*Power(r,3LL)*Power(xj,3LL) + 

                 255150LL*Power(r,4LL)*Power(xj,4LL) + 

                 85050LL*Power(r,5LL)*Power(xj,5LL) + 22680LL*Power(r,6LL)*Power(xj,6LL) + 

                 4860LL*Power(r,7LL)*Power(xj,7LL) + 810LL*Power(r,8LL)*Power(xj,8LL) + 

                 34LL*Power(r,9LL)*Power(xj,9LL)) + 

              3LL*Power(xi,14LL)*Power(xj,8LL)*

               (-19348875LL - 34827975LL*r*xj - 30958200LL*Power(r,2LL)*Power(xj,2LL) - 

                 18689580LL*Power(r,3LL)*Power(xj,3LL) - 

                 5847660LL*Power(r,4LL)*Power(xj,4LL) - 

                 3723300LL*Power(r,5LL)*Power(xj,5LL) - 

                 845040LL*Power(r,6LL)*Power(xj,6LL) - 

                 58680LL*Power(r,7LL)*Power(xj,7LL) + 1548LL*Power(r,8LL)*Power(xj,8LL) + 

                 236LL*Power(r,9LL)*Power(xj,9LL)) - 

              3LL*Power(xi,8LL)*Power(xj,14LL)*

               (-593408025LL + 946053675LL*r*xj - 

                 394427880LL*Power(r,2LL)*Power(xj,2LL) - 

                 315870660LL*Power(r,3LL)*Power(xj,3LL) - 

                 53891460LL*Power(r,4LL)*Power(xj,4LL) + 

                 910980LL*Power(r,5LL)*Power(xj,5LL) + 

                 1409520LL*Power(r,6LL)*Power(xj,6LL) + 

                 192168LL*Power(r,7LL)*Power(xj,7LL) + 11196LL*Power(r,8LL)*Power(xj,8LL) + 

                 236LL*Power(r,9LL)*Power(xj,9LL)))))/

       (42525LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,15LL)*Power(xi + xj,14LL)) - 

      (85050LL*exp(2LL*r*(xi + xj))*(xi + xj)*

          Power(Power(xi,2LL) - Power(xj,2LL),15LL) + 

         189LL*exp(2LL*r*xj)*Power(xj,12LL)*

          (-1400LL*Power(r,3LL)*Power(xi,22LL) - 50LL*Power(r,4LL)*Power(xi,23LL) - 

            7100LL*Power(r,3LL)*Power(xi,20LL)*Power(xj,2LL) - 

            140LL*Power(r,4LL)*Power(xi,21LL)*Power(xj,2LL) + 375LL*xi*Power(xj,18LL) + 

            600LL*r*Power(xi,2LL)*Power(xj,18LL) + 

            300LL*Power(r,2LL)*Power(xi,3LL)*Power(xj,18LL) - 

            210LL*Power(r,2LL)*Power(xi,21LL)*(75LL + Power(r,2LL)*Power(xj,2LL)) + 

            75LL*Power(xi,3LL)*Power(xj,16LL)*(-75LL + 2LL*Power(r,2LL)*Power(xj,2LL)) - 

            100LL*r*Power(xi,20LL)*(840LL + 71LL*Power(r,2LL)*Power(xj,2LL)) + 

            r*Power(xi,9LL)*Power(xj,10LL)*

             (249600LL*r*Power(xj,2LL) - 992LL*Power(r,3LL)*Power(xj,4LL)) + 

            20LL*r*Power(xi,17LL)*Power(xj,2LL)*

             (-1896LL*r*Power(xj,2LL) + 8LL*Power(r,3LL)*Power(xj,4LL)) + 

            5LL*r*Power(xi,5LL)*Power(xj,14LL)*

             (-900LL*r*Power(xj,2LL) + 8LL*Power(r,3LL)*Power(xj,4LL)) + 

            25LL*Power(xi,4LL)*Power(xj,14LL)*

             (-360LL*r*Power(xj,2LL) + 8LL*Power(r,3LL)*Power(xj,4LL)) - 

            375LL*Power(xi,6LL)*Power(xj,12LL)*

             (-168LL*r*Power(xj,2LL) + 8LL*Power(r,3LL)*Power(xj,4LL)) - 

            5LL*r*Power(xi,11LL)*Power(xj,8LL)*

             (98280LL*r*Power(xj,2LL) + 32LL*Power(r,3LL)*Power(xj,4LL)) + 

            5LL*r*Power(xi,7LL)*Power(xj,12LL)*

             (10304LL*r*Power(xj,2LL) + 56LL*Power(r,3LL)*Power(xj,4LL)) + 

            325LL*Power(xi,10LL)*Power(xj,8LL)*

             (-10680LL*r*Power(xj,2LL) + 160LL*Power(r,3LL)*Power(xj,4LL)) - 

            15LL*r*Power(xi,15LL)*Power(xj,4LL)*

             (-45920LL*r*Power(xj,2LL) + 208LL*Power(r,3LL)*Power(xj,4LL)) + 

            15LL*r*Power(xi,13LL)*Power(xj,6LL)*

             (-20280LL*r*Power(xj,2LL) + 208LL*Power(r,3LL)*Power(xj,4LL)) + 

            75LL*Power(xi,12LL)*Power(xj,6LL)*

             (51064LL*r*Power(xj,2LL) + 208LL*Power(r,3LL)*Power(xj,4LL)) + 

            60LL*Power(xi,18LL)*(-23880LL*r*Power(xj,2LL) + 

               344LL*Power(r,3LL)*Power(xj,4LL)) + 

            2LL*r*Power(xi,19LL)*(-70850LL*r*Power(xj,2LL) + 

               496LL*Power(r,3LL)*Power(xj,4LL)) + 

            100LL*Power(xi,16LL)*Power(xj,2LL)*

             (-26622LL*r*Power(xj,2LL) + 584LL*Power(r,3LL)*Power(xj,4LL)) - 

            5LL*Power(xi,8LL)*Power(xj,10LL)*

             (191880LL*r*Power(xj,2LL) + 3808LL*Power(r,3LL)*Power(xj,4LL)) - 

            15LL*Power(xi,14LL)*Power(xj,4LL)*

             (-315000LL*r*Power(xj,2LL) + 7280LL*Power(r,3LL)*Power(xj,4LL)) + 

            Power(xi,9LL)*Power(xj,10LL)*

             (4694625LL + 124800LL*Power(r,2LL)*Power(xj,2LL) - 

               248LL*Power(r,4LL)*Power(xj,4LL)) + 

            20LL*Power(xi,17LL)*Power(xj,2LL)*

             (-185895LL - 948LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) 

    + 5LL*Power(xi,5LL)*Power(xj,14LL)*(7875LL - 450LL*Power(r,2LL)*Power(xj,2LL) + 

               2LL*Power(r,4LL)*Power(xj,4LL)) - 

            5LL*Power(xi,11LL)*Power(xj,8LL)*

             (-2803125LL + 49140LL*Power(r,2LL)*Power(xj,2LL) + 

               8LL*Power(r,4LL)*Power(xj,4LL)) + 

            5LL*Power(xi,7LL)*Power(xj,12LL)*

             (-16965LL + 5152LL*Power(r,2LL)*Power(xj,2LL) + 14LL*Power(r,4LL)*Power(xj,4LL)) 

    - 15LL*Power(xi,15LL)*Power(xj,4LL)*(845085LL - 22960LL*Power(r,2LL)*Power(xj,2LL) + 

               52LL*Power(r,4LL)*Power(xj,4LL)) + 

            15LL*Power(xi,13LL)*Power(xj,6LL)*

             (-139125LL - 10140LL*Power(r,2LL)*Power(xj,2LL) + 

               52LL*Power(r,4LL)*Power(xj,4LL)) + 

            2LL*Power(xi,19LL)*(-89250LL - 35425LL*Power(r,2LL)*Power(xj,2LL) + 

               124LL*Power(r,4LL)*Power(xj,4LL))) + 

         378LL*exp(2LL*r*xj)*Power(xj,13LL)*

          (-350LL*Power(r,4LL)*Power(xi,22LL) - 10LL*Power(r,5LL)*Power(xi,23LL) + 

            225LL*Power(xj,18LL) + 375LL*r*xi*Power(xj,18LL) - 

            70LL*Power(r,3LL)*Power(xi,21LL)*(75LL + Power(r,2LL)*Power(xj,2LL)) + 

            75LL*r*Power(xi,3LL)*Power(xj,16LL)*(-75LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            75LL*Power(xi,2LL)*Power(xj,16LL)*(-45LL + 4LL*Power(r,2LL)*Power(xj,2LL)) - 

            50LL*Power(r,2LL)*Power(xi,20LL)*(840LL + 71LL*Power(r,2LL)*Power(xj,2LL)) + 

            r*Power(xi,9LL)*Power(xj,10LL)*

             (4694625LL + 124800LL*Power(r,2LL)*Power(xj,2LL) - 

               248LL*Power(r,4LL)*Power(xj,4LL)) + 

            20LL*r*Power(xi,17LL)*Power(xj,2LL)*

             (-185895LL - 948LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) 

    + 5LL*r*Power(xi,5LL)*Power(xj,14LL)*

             (7875LL - 450LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            25LL*Power(xi,4LL)*Power(xj,14LL)*

             (945LL - 180LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) - 

            375LL*Power(xi,6LL)*Power(xj,12LL)*

             (273LL - 84LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) - 

            5LL*r*Power(xi,11LL)*Power(xj,8LL)*

             (-2803125LL + 49140LL*Power(r,2LL)*Power(xj,2LL) + 

               8LL*Power(r,4LL)*Power(xj,4LL)) + 

            5LL*r*Power(xi,7LL)*Power(xj,12LL)*

             (-16965LL + 5152LL*Power(r,2LL)*Power(xj,2LL) + 14LL*Power(r,4LL)*Power(xj,4LL)) 

    + 325LL*Power(xi,10LL)*Power(xj,8LL)*

             (-60117LL - 5340LL*Power(r,2LL)*Power(xj,2LL) + 40LL*Power(r,4LL)*Power(xj,4LL)) 

    - 15LL*r*Power(xi,15LL)*Power(xj,4LL)*

             (845085LL - 22960LL*Power(r,2LL)*Power(xj,2LL) + 

               52LL*Power(r,4LL)*Power(xj,4LL)) + 

            15LL*r*Power(xi,13LL)*Power(xj,6LL)*

             (-139125LL - 10140LL*Power(r,2LL)*Power(xj,2LL) + 

               52LL*Power(r,4LL)*Power(xj,4LL)) + 

            75LL*Power(xi,12LL)*Power(xj,6LL)*

             (-729687LL + 25532LL*Power(r,2LL)*Power(xj,2LL) + 

               52LL*Power(r,4LL)*Power(xj,4LL)) + 

            60LL*Power(xi,18LL)*(-5355LL - 11940LL*Power(r,2LL)*Power(xj,2LL) + 

               86LL*Power(r,4LL)*Power(xj,4LL)) + 

            2LL*r*Power(xi,19LL)*(-89250LL - 35425LL*Power(r,2LL)*Power(xj,2LL) + 

               124LL*Power(r,4LL)*Power(xj,4LL)) + 

            100LL*Power(xi,16LL)*Power(xj,2LL)*

             (-79713LL - 13311LL*Power(r,2LL)*Power(xj,2LL) + 

               146LL*Power(r,4LL)*Power(xj,4LL)) - 

            5LL*Power(xi,8LL)*Power(xj,10LL)*

             (157365LL + 95940LL*Power(r,2LL)*Power(xj,2LL) + 

               952LL*Power(r,4LL)*Power(xj,4LL)) - 

            15LL*Power(xi,14LL)*Power(xj,4LL)*

             (2638467LL - 157500LL*Power(r,2LL)*Power(xj,2LL) + 

               1820LL*Power(r,4LL)*Power(xj,4LL))) + 

         exp(2LL*r*xi)*Power(xi,8LL)*

          (2LL*Power(xi,2LL)*Power(xj,20LL)*

             (1449175455LL*xj + 1066731120LL*r*Power(xj,2LL) + 

               343894005LL*Power(r,2LL)*Power(xj,3LL) + 

               60884460LL*Power(r,3LL)*Power(xj,4LL) + 

               5712525LL*Power(r,4LL)*Power(xj,5LL) + 

               110376LL*Power(r,5LL)*Power(xj,6LL) - 36666LL*Power(r,6LL)*Power(xj,7LL) - 

               4104LL*Power(r,7LL)*Power(xj,8LL) - 153LL*Power(r,8LL)*Power(xj,9LL)) + 

            42LL*Power(xi,4LL)*Power(xj,18LL)*

             (104824125LL*xj + 680400LL*r*Power(xj,2LL) - 

               27366255LL*Power(r,2LL)*Power(xj,3LL) - 

               11192580LL*Power(r,3LL)*Power(xj,4LL) - 

               2168775LL*Power(r,4LL)*Power(xj,5LL) - 

               234360LL*Power(r,5LL)*Power(xj,6LL) - 13230LL*Power(r,6LL)*Power(xj,7LL) - 

               216LL*Power(r,7LL)*Power(xj,8LL) + 9LL*Power(r,8LL)*Power(xj,9LL)) + 

            6LL*Power(xj,22LL)*(34459425LL*xj + 32432400LL*r*Power(xj,2LL) + 

               14189175LL*Power(r,2LL)*Power(xj,3LL) + 

               3783780LL*Power(r,3LL)*Power(xj,4LL) + 

               675675LL*Power(r,4LL)*Power(xj,5LL) + 83160LL*Power(r,5LL)*Power(xj,6LL) + 

               6930LL*Power(r,6LL)*Power(xj,7LL) + 360LL*Power(r,7LL)*Power(xj,8LL) + 

               9LL*Power(r,8LL)*Power(xj,9LL)) - 

            3LL*Power(xi,22LL)*(25515LL*xj + 45360LL*r*Power(xj,2LL) + 

               39690LL*Power(r,2LL)*Power(xj,3LL) + 22680LL*Power(r,3LL)*Power(xj,4LL) + 

               9450LL*Power(r,4LL)*Power(xj,5LL) + 3024LL*Power(r,5LL)*Power(xj,6LL) + 

               756LL*Power(r,6LL)*Power(xj,7LL) + 144LL*Power(r,7LL)*Power(xj,8LL) + 

               18LL*Power(r,8LL)*Power(xj,9LL)) - 

            21LL*Power(xi,18LL)*Power(xj,4LL)*

             (382725LL*xj + 680400LL*r*Power(xj,2LL) + 

               595350LL*Power(r,2LL)*Power(xj,3LL) + 340200LL*Power(r,3LL)*Power(xj,4LL) + 

               141750LL*Power(r,4LL)*Power(xj,5LL) + 45360LL*Power(r,5LL)*Power(xj,6LL) + 

               12852LL*Power(r,6LL)*Power(xj,7LL) + 1296LL*Power(r,7LL)*Power(xj,8LL) + 

               18LL*Power(r,8LL)*Power(xj,9LL)) + 

            54LL*Power(xi,6LL)*Power(xj,16LL)*

             (-73700865LL*xj - 108193680LL*r*Power(xj,2LL) - 

               24918705LL*Power(r,2LL)*Power(xj,3LL) + 

               3867780LL*Power(r,3LL)*Power(xj,4LL) + 

               2583735LL*Power(r,4LL)*Power(xj,5LL) + 

               484344LL*Power(r,5LL)*Power(xj,6LL) + 45038LL*Power(r,6LL)*Power(xj,7LL) + 

               2008LL*Power(r,7LL)*Power(xj,8LL) + 27LL*Power(r,8LL)*Power(xj,9LL)) - 

            315LL*Power(xi,12LL)*Power(xj,10LL)*

             (-710073LL*xj - 1611792LL*r*Power(xj,2LL) - 

               304668LL*Power(r,2LL)*Power(xj,3LL) - 

               1035216LL*Power(r,3LL)*Power(xj,4LL) - 

               454860LL*Power(r,4LL)*Power(xj,5LL) - 58464LL*Power(r,5LL)*Power(xj,6LL) + 

               840LL*Power(r,6LL)*Power(xj,7LL) + 672LL*Power(r,7LL)*Power(xj,8LL) + 

               36LL*Power(r,8LL)*Power(xj,9LL)) + 

            315LL*Power(xi,10LL)*Power(xj,12LL)*

             (-2656395LL*xj + 2373840LL*r*Power(xj,2LL) - 

               3466260LL*Power(r,2LL)*Power(xj,3LL) - 

               2573424LL*Power(r,3LL)*Power(xj,4LL) - 

               467460LL*Power(r,4LL)*Power(xj,5LL) + 2016LL*Power(r,5LL)*Power(xj,6LL) + 

               9576LL*Power(r,6LL)*Power(xj,7LL) + 1056LL*Power(r,7LL)*Power(xj,8LL) + 

               36LL*Power(r,8LL)*Power(xj,9LL)) - 

            27LL*Power(xi,16LL)*Power(xj,6LL)*

             (-1289925LL*xj - 2293200LL*r*Power(xj,2LL) - 

               2006550LL*Power(r,2LL)*Power(xj,3LL) - 

               1146600LL*Power(r,3LL)*Power(xj,4LL) - 

               450030LL*Power(r,4LL)*Power(xj,5LL) - 197232LL*Power(r,5LL)*Power(xj,6LL) - 

               33684LL*Power(r,6LL)*Power(xj,7LL) - 1424LL*Power(r,7LL)*Power(xj,8LL) + 

               54LL*Power(r,8LL)*Power(xj,9LL)) + 

            Power(xi,20LL)*Power(xj,2LL)*

             (1148175LL*xj + 2041200LL*r*Power(xj,2LL) + 

               1786050LL*Power(r,2LL)*Power(xj,3LL) + 

               1020600LL*Power(r,3LL)*Power(xj,4LL) + 

               425250LL*Power(r,4LL)*Power(xj,5LL) + 136080LL*Power(r,5LL)*Power(xj,6LL) + 

               34020LL*Power(r,6LL)*Power(xj,7LL) + 6480LL*Power(r,7LL)*Power(xj,8LL) + 

               306LL*Power(r,8LL)*Power(xj,9LL)) + 

            3LL*Power(xi,14LL)*Power(xj,8LL)*

             (-34827975LL*xj - 61916400LL*r*Power(xj,2LL) - 

               56068740LL*Power(r,2LL)*Power(xj,3LL) - 

               23390640LL*Power(r,3LL)*Power(xj,4LL) - 

               18616500LL*Power(r,4LL)*Power(xj,5LL) - 

               5070240LL*Power(r,5LL)*Power(xj,6LL) - 

               410760LL*Power(r,6LL)*Power(xj,7LL) + 12384LL*Power(r,7LL)*Power(xj,8LL) + 

               2124LL*Power(r,8LL)*Power(xj,9LL)) - 

            3LL*Power(xi,8LL)*Power(xj,14LL)*

             (946053675LL*xj - 788855760LL*r*Power(xj,2LL) - 

               947611980LL*Power(r,2LL)*Power(xj,3LL) - 

               215565840LL*Power(r,3LL)*Power(xj,4LL) + 

               4554900LL*Power(r,4LL)*Power(xj,5LL) + 

               8457120LL*Power(r,5LL)*Power(xj,6LL) + 

               1345176LL*Power(r,6LL)*Power(xj,7LL) + 89568LL*Power(r,7LL)*Power(xj,8LL) + 

               2124LL*Power(r,8LL)*Power(xj,9LL))) + 

         2LL*exp(2LL*r*xi)*Power(xi,9LL)*

          (2LL*Power(xi,2LL)*Power(xj,20LL)*

             (1782492075LL + 1449175455LL*r*xj + 533365560LL*Power(r,2LL)*Power(xj,2LL) + 

               114631335LL*Power(r,3LL)*Power(xj,3LL) + 

               15221115LL*Power(r,4LL)*Power(xj,4LL) + 

               1142505LL*Power(r,5LL)*Power(xj,5LL) + 18396LL*Power(r,6LL)*Power(xj,6LL) - 

               5238LL*Power(r,7LL)*Power(xj,7LL) - 513LL*Power(r,8LL)*Power(xj,8LL) - 

               17LL*Power(r,9LL)*Power(xj,9LL)) + 

            42LL*Power(xi,4LL)*Power(xj,18LL)*

             (251336925LL + 104824125LL*r*xj + 340200LL*Power(r,2LL)*Power(xj,2LL) - 

               9122085LL*Power(r,3LL)*Power(xj,3LL) - 

               2798145LL*Power(r,4LL)*Power(xj,4LL) - 433755LL*Power(r,5LL)*Power(xj,5LL) - 

               39060LL*Power(r,6LL)*Power(xj,6LL) - 1890LL*Power(r,7LL)*Power(xj,7LL) - 

               27LL*Power(r,8LL)*Power(xj,8LL) + Power(r,9LL)*Power(xj,9LL)) + 

            6LL*Power(xj,22LL)*(34459425LL + 34459425LL*r*xj + 

               16216200LL*Power(r,2LL)*Power(xj,2LL) + 

               4729725LL*Power(r,3LL)*Power(xj,3LL) + 945945LL*Power(r,4LL)*Power(xj,4LL) + 

               135135LL*Power(r,5LL)*Power(xj,5LL) + 13860LL*Power(r,6LL)*Power(xj,6LL) + 

               990LL*Power(r,7LL)*Power(xj,7LL) + 45LL*Power(r,8LL)*Power(xj,8LL) + 

               Power(r,9LL)*Power(xj,9LL)) - 

            3LL*Power(xi,22LL)*(14175LL + 25515LL*r*xj + 22680LL*Power(r,2LL)*Power(xj,2LL) + 

               13230LL*Power(r,3LL)*Power(xj,3LL) + 5670LL*Power(r,4LL)*Power(xj,4LL) + 

               1890LL*Power(r,5LL)*Power(xj,5LL) + 504LL*Power(r,6LL)*Power(xj,6LL) + 

               108LL*Power(r,7LL)*Power(xj,7LL) + 18LL*Power(r,8LL)*Power(xj,8LL) + 

               2LL*Power(r,9LL)*Power(xj,9LL)) - 

            21LL*Power(xi,18LL)*Power(xj,4LL)*

             (212625LL + 382725LL*r*xj + 340200LL*Power(r,2LL)*Power(xj,2LL) + 

               198450LL*Power(r,3LL)*Power(xj,3LL) + 85050LL*Power(r,4LL)*Power(xj,4LL) + 

               28350LL*Power(r,5LL)*Power(xj,5LL) + 7560LL*Power(r,6LL)*Power(xj,6LL) + 

               1836LL*Power(r,7LL)*Power(xj,7LL) + 162LL*Power(r,8LL)*Power(xj,8LL) + 

               2LL*Power(r,9LL)*Power(xj,9LL)) + 

            54LL*Power(xi,6LL)*Power(xj,16LL)*

             (133451955LL - 73700865LL*r*xj - 54096840LL*Power(r,2LL)*Power(xj,2LL) - 

               8306235LL*Power(r,3LL)*Power(xj,3LL) + 966945LL*Power(r,4LL)*Power(xj,4LL) + 

               516747LL*Power(r,5LL)*Power(xj,5LL) + 80724LL*Power(r,6LL)*Power(xj,6LL) + 

               6434LL*Power(r,7LL)*Power(xj,7LL) + 251LL*Power(r,8LL)*Power(xj,8LL) + 

               3LL*Power(r,9LL)*Power(xj,9LL)) - 

            315LL*Power(xi,12LL)*Power(xj,10LL)*

             (-405405LL - 710073LL*r*xj - 805896LL*Power(r,2LL)*Power(xj,2LL) - 

               101556LL*Power(r,3LL)*Power(xj,3LL) - 258804LL*Power(r,4LL)*Power(xj,4LL) - 

               90972LL*Power(r,5LL)*Power(xj,5LL) - 9744LL*Power(r,6LL)*Power(xj,6LL) + 

               120LL*Power(r,7LL)*Power(xj,7LL) + 84LL*Power(r,8LL)*Power(xj,8LL) + 

               4LL*Power(r,9LL)*Power(xj,9LL)) + 

            315LL*Power(xi,10LL)*Power(xj,12LL)*

             (-482895LL - 2656395LL*r*xj + 1186920LL*Power(r,2LL)*Power(xj,2LL) - 

               1155420LL*Power(r,3LL)*Power(xj,3LL) - 643356LL*Power(r,4LL)*Power(xj,4LL) - 

               93492LL*Power(r,5LL)*Power(xj,5LL) + 336LL*Power(r,6LL)*Power(xj,6LL) + 

               1368LL*Power(r,7LL)*Power(xj,7LL) + 132LL*Power(r,8LL)*Power(xj,8LL) + 

               4LL*Power(r,9LL)*Power(xj,9LL)) - 

            27LL*Power(xi,16LL)*Power(xj,6LL)*

             (-716625LL - 1289925LL*r*xj - 1146600LL*Power(r,2LL)*Power(xj,2LL) - 

               668850LL*Power(r,3LL)*Power(xj,3LL) - 286650LL*Power(r,4LL)*Power(xj,4LL) - 

               90006LL*Power(r,5LL)*Power(xj,5LL) - 32872LL*Power(r,6LL)*Power(xj,6LL) - 

               4812LL*Power(r,7LL)*Power(xj,7LL) - 178LL*Power(r,8LL)*Power(xj,8LL) + 

               6LL*Power(r,9LL)*Power(xj,9LL)) + 

            Power(xi,20LL)*Power(xj,2LL)*

             (637875LL + 1148175LL*r*xj + 1020600LL*Power(r,2LL)*Power(xj,2LL) + 

               595350LL*Power(r,3LL)*Power(xj,3LL) + 255150LL*Power(r,4LL)*Power(xj,4LL) + 

               85050LL*Power(r,5LL)*Power(xj,5LL) + 22680LL*Power(r,6LL)*Power(xj,6LL) + 

               4860LL*Power(r,7LL)*Power(xj,7LL) + 810LL*Power(r,8LL)*Power(xj,8LL) + 

               34LL*Power(r,9LL)*Power(xj,9LL)) + 

            3LL*Power(xi,14LL)*Power(xj,8LL)*

             (-19348875LL - 34827975LL*r*xj - 30958200LL*Power(r,2LL)*Power(xj,2LL) - 

               18689580LL*Power(r,3LL)*Power(xj,3LL) - 

               5847660LL*Power(r,4LL)*Power(xj,4LL) - 

               3723300LL*Power(r,5LL)*Power(xj,5LL) - 845040LL*Power(r,6LL)*Power(xj,6LL) - 

               58680LL*Power(r,7LL)*Power(xj,7LL) + 1548LL*Power(r,8LL)*Power(xj,8LL) + 

               236LL*Power(r,9LL)*Power(xj,9LL)) - 

            3LL*Power(xi,8LL)*Power(xj,14LL)*

             (-593408025LL + 946053675LL*r*xj - 394427880LL*Power(r,2LL)*Power(xj,2LL) - 

               315870660LL*Power(r,3LL)*Power(xj,3LL) - 

               53891460LL*Power(r,4LL)*Power(xj,4LL) + 910980LL*Power(r,5LL)*Power(xj,5LL) + 

               1409520LL*Power(r,6LL)*Power(xj,6LL) + 192168LL*Power(r,7LL)*Power(xj,7LL) + 

               11196LL*Power(r,8LL)*Power(xj,8LL) + 236LL*Power(r,9LL)*Power(xj,9LL))))/

       (42525LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,15LL)*Power(xi + xj,15LL))

    ; }
   
  }
  return S;
}

cl_R DSlater_3S_1S(cl_R r,cl_R xi,cl_R xj)
{
  return DSlater_1S_3S(r,xj,xi);
}

cl_R DSlater_3S_2S(cl_R r,cl_R xi,cl_R xj)
{
  return DSlater_2S_3S(r,xj,xi);
}

cl_R DSlater_4S_4S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = 0LL

    ; } else {  S = -(-150568359566625LL*xi + 167382319104000LL*exp(2LL*r*xi)*xi - 

          267508800058500LL*r*Power(xi,2LL) - 

          234428725030500LL*Power(r,2LL)*Power(xi,3LL) - 

          134962892064000LL*Power(r,3LL)*Power(xi,4LL) - 

          57353780484000LL*Power(r,4LL)*Power(xi,5LL) - 

          19160153812800LL*Power(r,5LL)*Power(xi,6LL) - 

          5229789364800LL*Power(r,6LL)*Power(xi,7LL) - 

          1195587993600LL*Power(r,7LL)*Power(xi,8LL) - 

          232475443200LL*Power(r,8LL)*Power(xi,9LL) - 

          38745907200LL*Power(r,9LL)*Power(xi,10LL) - 

          5535129600LL*Power(r,10LL)*Power(xi,11LL) - 

          670924800LL*Power(r,11LL)*Power(xi,12LL) - 

          67092480LL*Power(r,12LL)*Power(xi,13LL) - 

          5160960LL*Power(r,13LL)*Power(xi,14LL) - 245760LL*Power(r,14LL)*Power(xi,15LL))/

       (8.3691159552e13*exp(2LL*r*xi)*r) + 

      (-83691159552000LL + 83691159552000LL*exp(2LL*r*xi) - 

         150568359566625LL*r*xi - 133754400029250LL*Power(r,2LL)*Power(xi,2LL) - 

         78142908343500LL*Power(r,3LL)*Power(xi,3LL) - 

         33740723016000LL*Power(r,4LL)*Power(xi,4LL) - 

         11470756096800LL*Power(r,5LL)*Power(xi,5LL) - 

         3193358968800LL*Power(r,6LL)*Power(xi,6LL) - 

         747112766400LL*Power(r,7LL)*Power(xi,7LL) - 

         149448499200LL*Power(r,8LL)*Power(xi,8LL) - 

         25830604800LL*Power(r,9LL)*Power(xi,9LL) - 

         3874590720LL*Power(r,10LL)*Power(xi,10LL) - 

         503193600LL*Power(r,11LL)*Power(xi,11LL) - 

         55910400LL*Power(r,12LL)*Power(xi,12LL) - 5160960LL*Power(r,13LL)*Power(xi,13LL) - 

         368640LL*Power(r,14LL)*Power(xi,14LL) - 16384LL*Power(r,15LL)*Power(xi,15LL))/

       (8.3691159552e13*exp(2LL*r*xi)*Power(r,2LL)) + 

      (xi*(-83691159552000LL + 83691159552000LL*exp(2LL*r*xi) - 

           150568359566625LL*r*xi - 133754400029250LL*Power(r,2LL)*Power(xi,2LL) - 

           78142908343500LL*Power(r,3LL)*Power(xi,3LL) - 

           33740723016000LL*Power(r,4LL)*Power(xi,4LL) - 

           11470756096800LL*Power(r,5LL)*Power(xi,5LL) - 

           3193358968800LL*Power(r,6LL)*Power(xi,6LL) - 

           747112766400LL*Power(r,7LL)*Power(xi,7LL) - 

           149448499200LL*Power(r,8LL)*Power(xi,8LL) - 

           25830604800LL*Power(r,9LL)*Power(xi,9LL) - 

           3874590720LL*Power(r,10LL)*Power(xi,10LL) - 

           503193600LL*Power(r,11LL)*Power(xi,11LL) - 

           55910400LL*Power(r,12LL)*Power(xi,12LL) - 

           5160960LL*Power(r,13LL)*Power(xi,13LL) - 368640LL*Power(r,14LL)*Power(xi,14LL) - 

           16384LL*Power(r,15LL)*Power(xi,15LL)))/(4.1845579776e13*exp(2LL*r*xi)*r)

    ; }
 
  }
  else {
      if (r == 0LL) {  S = 0LL

    ; } else {  S = (1260LL*exp(2LL*r*(xi + xj))*Power(Power(xi,2LL) - Power(xj,2LL),15LL) + 

         exp(2LL*r*xj)*Power(xj,10LL)*

          (-3276LL*Power(r,5LL)*Power(xi,25LL) - 168LL*Power(r,6LL)*Power(xi,26LL) - 

            4LL*Power(r,7LL)*Power(xi,27LL) + 1260LL*Power(xj,20LL) + 

            2205LL*r*xi*Power(xj,20LL) + 

            1890LL*Power(xi,2LL)*Power(xj,18LL)*(-10LL + Power(r,2LL)*Power(xj,2LL)) - 

            420LL*Power(r,4LL)*Power(xi,24LL)*(91LL + Power(r,2LL)*Power(xj,2LL)) + 

            525LL*r*Power(xi,3LL)*Power(xj,18LL)*(-63LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            42LL*Power(r,3LL)*Power(xi,23LL)*

             (-6825LL - 405LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            63LL*r*Power(xi,5LL)*Power(xj,16LL)*

             (3675LL - 250LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            210LL*Power(xi,4LL)*Power(xj,16LL)*

             (630LL - 135LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            252LL*Power(r,2LL)*Power(xi,22LL)*

             (-5460LL - 1225LL*Power(r,2LL)*Power(xj,2LL) + 17LL*Power(r,4LL)*Power(xj,4LL)) 

    - 1260LL*r*Power(xi,17LL)*Power(xj,4LL)*

             (141729LL - 10145LL*Power(r,2LL)*Power(xj,2LL) + 

               116LL*Power(r,4LL)*Power(xj,4LL)) + 

            21LL*r*Power(xi,9LL)*Power(xj,12LL)*

             (164775LL - 18460LL*Power(r,2LL)*Power(xj,2LL) + 

               828LL*Power(r,4LL)*Power(xj,4LL)) + 

            14LL*Power(xi,6LL)*Power(xj,14LL)*

             (-40950LL + 14175LL*Power(r,2LL)*Power(xj,2LL) - 

               450LL*Power(r,4LL)*Power(xj,4LL) + 2LL*Power(r,6LL)*Power(xj,6LL)) - 

            210LL*Power(xi,8LL)*Power(xj,12LL)*

             (-8190LL + 4095LL*Power(r,2LL)*Power(xj,2LL) - 

               210LL*Power(r,4LL)*Power(xj,4LL) + 2LL*Power(r,6LL)*Power(xj,6LL)) + 

            42LL*Power(xi,10LL)*Power(xj,10LL)*

             (-209430LL - 2925LL*Power(r,2LL)*Power(xj,2LL) - 

               8840LL*Power(r,4LL)*Power(xj,4LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) + 

            r*Power(xi,7LL)*Power(xj,14LL)*

             (-1003275LL + 110250LL*Power(r,2LL)*Power(xj,2LL) - 

               1890LL*Power(r,4LL)*Power(xj,4LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) - 

            21LL*r*Power(xi,11LL)*Power(xj,10LL)*

             (-1033695LL - 218400LL*Power(r,2LL)*Power(xj,2LL) + 

               552LL*Power(r,4LL)*Power(xj,4LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) + 

            280LL*Power(xi,18LL)*Power(xj,2LL)*

             (-385560LL - 73953LL*Power(r,2LL)*Power(xj,2LL) + 

               2370LL*Power(r,4LL)*Power(xj,4LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) - 

            35LL*r*Power(xi,15LL)*Power(xj,6LL)*

             (-1565613LL + 359520LL*Power(r,2LL)*Power(xj,2LL) - 

               7020LL*Power(r,4LL)*Power(xj,4LL) + 8LL*Power(r,6LL)*Power(xj,6LL)) + 

            14LL*r*Power(xi,19LL)*Power(xj,2LL)*

             (-4980150LL + 126765LL*Power(r,2LL)*Power(xj,2LL) - 

               3852LL*Power(r,4LL)*Power(xj,4LL) + 20LL*Power(r,6LL)*Power(xj,6LL)) - 

            630LL*Power(xi,14LL)*Power(xj,6LL)*

             (708714LL - 14385LL*Power(r,2LL)*Power(xj,2LL) - 

               2340LL*Power(r,4LL)*Power(xj,4LL) + 20LL*Power(r,6LL)*Power(xj,6LL)) + 

            210LL*Power(xi,16LL)*Power(xj,4LL)*

             (-2087532LL + 328491LL*Power(r,2LL)*Power(xj,2LL) - 

               11740LL*Power(r,4LL)*Power(xj,4LL) + 52LL*Power(r,6LL)*Power(xj,6LL)) - 

            84LL*Power(xi,20LL)*(59670LL + 236250LL*Power(r,2LL)*Power(xj,2LL) - 

               8745LL*Power(r,4LL)*Power(xj,4LL) + 92LL*Power(r,6LL)*Power(xj,6LL)) - 

            2LL*r*Power(xi,21LL)*(1949220LL + 1598625LL*Power(r,2LL)*Power(xj,2LL) - 

               41391LL*Power(r,4LL)*Power(xj,4LL) + 128LL*Power(r,6LL)*Power(xj,6LL)) + 

            r*Power(xi,13LL)*Power(xj,8LL)*

             (173037375LL - 2784600LL*Power(r,2LL)*Power(xj,2LL) - 

               112140LL*Power(r,4LL)*Power(xj,4LL) + 256LL*Power(r,6LL)*Power(xj,6LL)) + 

            14LL*Power(xi,12LL)*Power(xj,8LL)*

             (-7260750LL - 2521935LL*Power(r,2LL)*Power(xj,2LL) + 

               19500LL*Power(r,4LL)*Power(xj,4LL) + 344LL*Power(r,6LL)*Power(xj,6LL))) + 

         exp(2LL*r*xi)*Power(xi,10LL)*

          (210LL*Power(xi,2LL)*Power(xj,18LL)*

             (514080LL + 332010LL*r*xj + 94500LL*Power(r,2LL)*Power(xj,2LL) + 

               15225LL*Power(r,3LL)*Power(xj,3LL) + 1470LL*Power(r,4LL)*Power(xj,4LL) + 

               81LL*Power(r,5LL)*Power(xj,5LL) + 2LL*Power(r,6LL)*Power(xj,6LL)) + 

            105LL*Power(xi,18LL)*Power(xj,2LL)*

             (180LL + 315LL*r*xj + 270LL*Power(r,2LL)*Power(xj,2LL) + 

               150LL*Power(r,3LL)*Power(xj,3LL) + 60LL*Power(r,4LL)*Power(xj,4LL) + 

               18LL*Power(r,5LL)*Power(xj,5LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) - 

            1365LL*Power(xi,10LL)*Power(xj,10LL)*

             (-6444LL + 15903LL*r*xj - 25866LL*Power(r,2LL)*Power(xj,2LL) - 

               2040LL*Power(r,3LL)*Power(xj,3LL) + 1080LL*Power(r,4LL)*Power(xj,4LL) + 

               180LL*Power(r,5LL)*Power(xj,5LL) + 8LL*Power(r,6LL)*Power(xj,6LL)) + 

            Power(xi,14LL)*Power(xj,6LL)*

             (573300LL + 1003275LL*r*xj + 859950LL*Power(r,2LL)*Power(xj,2LL) + 

               387660LL*Power(r,3LL)*Power(xj,3LL) + 371280LL*Power(r,4LL)*Power(xj,4LL) + 

               11592LL*Power(r,5LL)*Power(xj,5LL) - 4816LL*Power(r,6LL)*Power(xj,6LL) - 

               256LL*Power(r,7LL)*Power(xj,7LL)) + 

            2LL*Power(xj,20LL)*(2506140LL + 1949220LL*r*xj + 

               687960LL*Power(r,2LL)*Power(xj,2LL) + 143325LL*Power(r,3LL)*Power(xj,3LL) + 

               19110LL*Power(r,4LL)*Power(xj,4LL) + 1638LL*Power(r,5LL)*Power(xj,5LL) + 

               84LL*Power(r,6LL)*Power(xj,6LL) + 2LL*Power(r,7LL)*Power(xj,7LL)) - 

            42LL*Power(xi,4LL)*Power(xj,16LL)*

             (-10437660LL - 4251870LL*r*xj - 493020LL*Power(r,2LL)*Power(xj,2LL) + 

               42255LL*Power(r,3LL)*Power(xj,3LL) + 17490LL*Power(r,4LL)*Power(xj,4LL) + 

               1971LL*Power(r,5LL)*Power(xj,5LL) + 102LL*Power(r,6LL)*Power(xj,6LL) + 

               2LL*Power(r,7LL)*Power(xj,7LL)) + 

            21LL*Power(xi,16LL)*Power(xj,4LL)*

             (-6300LL - 11025LL*r*xj - 9450LL*Power(r,2LL)*Power(xj,2LL) - 

               5250LL*Power(r,3LL)*Power(xj,3LL) - 2100LL*Power(r,4LL)*Power(xj,4LL) - 

               828LL*Power(r,5LL)*Power(xj,5LL) - 8LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL)) - 

            Power(xi,20LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r,2LL)*Power(xj,2LL) + 

               1050LL*Power(r,3LL)*Power(xj,3LL) + 420LL*Power(r,4LL)*Power(xj,4LL) + 

               126LL*Power(r,5LL)*Power(xj,5LL) + 28LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL)) - 

            35LL*Power(xi,8LL)*Power(xj,12LL)*

             (-2904300LL + 4943925LL*r*xj + 258930LL*Power(r,2LL)*Power(xj,2LL) - 

               359520LL*Power(r,3LL)*Power(xj,3LL) - 70440LL*Power(r,4LL)*Power(xj,4LL) - 

               4176LL*Power(r,5LL)*Power(xj,5LL) + 32LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) + 

            35LL*Power(xi,12LL)*Power(xj,8LL)*

             (-49140LL - 98865LL*r*xj + 3510LL*Power(r,2LL)*Power(xj,2LL) - 

               131040LL*Power(r,3LL)*Power(xj,3LL) - 7800LL*Power(r,4LL)*Power(xj,4LL) + 

               3204LL*Power(r,5LL)*Power(xj,5LL) + 360LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) + 

            Power(xi,6LL)*Power(xj,14LL)*

             (446489820LL - 54796455LL*r*xj - 68983110LL*Power(r,2LL)*Power(xj,2LL) - 

               12782700LL*Power(r,3LL)*Power(xj,3LL) - 

               663600LL*Power(r,4LL)*Power(xj,4LL) + 53928LL*Power(r,5LL)*Power(xj,5LL) + 

               7728LL*Power(r,6LL)*Power(xj,6LL) + 256LL*Power(r,7LL)*Power(xj,7LL))))/

       (1260LL*exp(2LL*r*(xi + xj))*Power(r,2LL)*Power(xi - xj,15LL)*

         Power(xi + xj,15LL)) + (1260LL*exp(2LL*r*(xi + xj))*

          Power(Power(xi,2LL) - Power(xj,2LL),15LL) + 

         exp(2LL*r*xj)*Power(xj,10LL)*

          (-3276LL*Power(r,5LL)*Power(xi,25LL) - 168LL*Power(r,6LL)*Power(xi,26LL) - 

            4LL*Power(r,7LL)*Power(xi,27LL) + 1260LL*Power(xj,20LL) + 

            2205LL*r*xi*Power(xj,20LL) + 

            1890LL*Power(xi,2LL)*Power(xj,18LL)*(-10LL + Power(r,2LL)*Power(xj,2LL)) - 

            420LL*Power(r,4LL)*Power(xi,24LL)*(91LL + Power(r,2LL)*Power(xj,2LL)) + 

            525LL*r*Power(xi,3LL)*Power(xj,18LL)*(-63LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            42LL*Power(r,3LL)*Power(xi,23LL)*

             (-6825LL - 405LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            63LL*r*Power(xi,5LL)*Power(xj,16LL)*

             (3675LL - 250LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            210LL*Power(xi,4LL)*Power(xj,16LL)*

             (630LL - 135LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            252LL*Power(r,2LL)*Power(xi,22LL)*

             (-5460LL - 1225LL*Power(r,2LL)*Power(xj,2LL) + 17LL*Power(r,4LL)*Power(xj,4LL)) 

    - 1260LL*r*Power(xi,17LL)*Power(xj,4LL)*

             (141729LL - 10145LL*Power(r,2LL)*Power(xj,2LL) + 

               116LL*Power(r,4LL)*Power(xj,4LL)) + 

            21LL*r*Power(xi,9LL)*Power(xj,12LL)*

             (164775LL - 18460LL*Power(r,2LL)*Power(xj,2LL) + 

               828LL*Power(r,4LL)*Power(xj,4LL)) + 

            14LL*Power(xi,6LL)*Power(xj,14LL)*

             (-40950LL + 14175LL*Power(r,2LL)*Power(xj,2LL) - 

               450LL*Power(r,4LL)*Power(xj,4LL) + 2LL*Power(r,6LL)*Power(xj,6LL)) - 

            210LL*Power(xi,8LL)*Power(xj,12LL)*

             (-8190LL + 4095LL*Power(r,2LL)*Power(xj,2LL) - 

               210LL*Power(r,4LL)*Power(xj,4LL) + 2LL*Power(r,6LL)*Power(xj,6LL)) + 

            42LL*Power(xi,10LL)*Power(xj,10LL)*

             (-209430LL - 2925LL*Power(r,2LL)*Power(xj,2LL) - 

               8840LL*Power(r,4LL)*Power(xj,4LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) + 

            r*Power(xi,7LL)*Power(xj,14LL)*

             (-1003275LL + 110250LL*Power(r,2LL)*Power(xj,2LL) - 

               1890LL*Power(r,4LL)*Power(xj,4LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) - 

            21LL*r*Power(xi,11LL)*Power(xj,10LL)*

             (-1033695LL - 218400LL*Power(r,2LL)*Power(xj,2LL) + 

               552LL*Power(r,4LL)*Power(xj,4LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) + 

            280LL*Power(xi,18LL)*Power(xj,2LL)*

             (-385560LL - 73953LL*Power(r,2LL)*Power(xj,2LL) + 

               2370LL*Power(r,4LL)*Power(xj,4LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) - 

            35LL*r*Power(xi,15LL)*Power(xj,6LL)*

             (-1565613LL + 359520LL*Power(r,2LL)*Power(xj,2LL) - 

               7020LL*Power(r,4LL)*Power(xj,4LL) + 8LL*Power(r,6LL)*Power(xj,6LL)) + 

            14LL*r*Power(xi,19LL)*Power(xj,2LL)*

             (-4980150LL + 126765LL*Power(r,2LL)*Power(xj,2LL) - 

               3852LL*Power(r,4LL)*Power(xj,4LL) + 20LL*Power(r,6LL)*Power(xj,6LL)) - 

            630LL*Power(xi,14LL)*Power(xj,6LL)*

             (708714LL - 14385LL*Power(r,2LL)*Power(xj,2LL) - 

               2340LL*Power(r,4LL)*Power(xj,4LL) + 20LL*Power(r,6LL)*Power(xj,6LL)) + 

            210LL*Power(xi,16LL)*Power(xj,4LL)*

             (-2087532LL + 328491LL*Power(r,2LL)*Power(xj,2LL) - 

               11740LL*Power(r,4LL)*Power(xj,4LL) + 52LL*Power(r,6LL)*Power(xj,6LL)) - 

            84LL*Power(xi,20LL)*(59670LL + 236250LL*Power(r,2LL)*Power(xj,2LL) - 

               8745LL*Power(r,4LL)*Power(xj,4LL) + 92LL*Power(r,6LL)*Power(xj,6LL)) - 

            2LL*r*Power(xi,21LL)*(1949220LL + 1598625LL*Power(r,2LL)*Power(xj,2LL) - 

               41391LL*Power(r,4LL)*Power(xj,4LL) + 128LL*Power(r,6LL)*Power(xj,6LL)) + 

            r*Power(xi,13LL)*Power(xj,8LL)*

             (173037375LL - 2784600LL*Power(r,2LL)*Power(xj,2LL) - 

               112140LL*Power(r,4LL)*Power(xj,4LL) + 256LL*Power(r,6LL)*Power(xj,6LL)) + 

            14LL*Power(xi,12LL)*Power(xj,8LL)*

             (-7260750LL - 2521935LL*Power(r,2LL)*Power(xj,2LL) + 

               19500LL*Power(r,4LL)*Power(xj,4LL) + 344LL*Power(r,6LL)*Power(xj,6LL))) + 

         exp(2LL*r*xi)*Power(xi,10LL)*

          (210LL*Power(xi,2LL)*Power(xj,18LL)*

             (514080LL + 332010LL*r*xj + 94500LL*Power(r,2LL)*Power(xj,2LL) + 

               15225LL*Power(r,3LL)*Power(xj,3LL) + 1470LL*Power(r,4LL)*Power(xj,4LL) + 

               81LL*Power(r,5LL)*Power(xj,5LL) + 2LL*Power(r,6LL)*Power(xj,6LL)) + 

            105LL*Power(xi,18LL)*Power(xj,2LL)*

             (180LL + 315LL*r*xj + 270LL*Power(r,2LL)*Power(xj,2LL) + 

               150LL*Power(r,3LL)*Power(xj,3LL) + 60LL*Power(r,4LL)*Power(xj,4LL) + 

               18LL*Power(r,5LL)*Power(xj,5LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) - 

            1365LL*Power(xi,10LL)*Power(xj,10LL)*

             (-6444LL + 15903LL*r*xj - 25866LL*Power(r,2LL)*Power(xj,2LL) - 

               2040LL*Power(r,3LL)*Power(xj,3LL) + 1080LL*Power(r,4LL)*Power(xj,4LL) + 

               180LL*Power(r,5LL)*Power(xj,5LL) + 8LL*Power(r,6LL)*Power(xj,6LL)) + 

            Power(xi,14LL)*Power(xj,6LL)*

             (573300LL + 1003275LL*r*xj + 859950LL*Power(r,2LL)*Power(xj,2LL) + 

               387660LL*Power(r,3LL)*Power(xj,3LL) + 371280LL*Power(r,4LL)*Power(xj,4LL) + 

               11592LL*Power(r,5LL)*Power(xj,5LL) - 4816LL*Power(r,6LL)*Power(xj,6LL) - 

               256LL*Power(r,7LL)*Power(xj,7LL)) + 

            2LL*Power(xj,20LL)*(2506140LL + 1949220LL*r*xj + 

               687960LL*Power(r,2LL)*Power(xj,2LL) + 143325LL*Power(r,3LL)*Power(xj,3LL) + 

               19110LL*Power(r,4LL)*Power(xj,4LL) + 1638LL*Power(r,5LL)*Power(xj,5LL) + 

               84LL*Power(r,6LL)*Power(xj,6LL) + 2LL*Power(r,7LL)*Power(xj,7LL)) - 

            42LL*Power(xi,4LL)*Power(xj,16LL)*

             (-10437660LL - 4251870LL*r*xj - 493020LL*Power(r,2LL)*Power(xj,2LL) + 

               42255LL*Power(r,3LL)*Power(xj,3LL) + 17490LL*Power(r,4LL)*Power(xj,4LL) + 

               1971LL*Power(r,5LL)*Power(xj,5LL) + 102LL*Power(r,6LL)*Power(xj,6LL) + 

               2LL*Power(r,7LL)*Power(xj,7LL)) + 

            21LL*Power(xi,16LL)*Power(xj,4LL)*

             (-6300LL - 11025LL*r*xj - 9450LL*Power(r,2LL)*Power(xj,2LL) - 

               5250LL*Power(r,3LL)*Power(xj,3LL) - 2100LL*Power(r,4LL)*Power(xj,4LL) - 

               828LL*Power(r,5LL)*Power(xj,5LL) - 8LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL)) - 

            Power(xi,20LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r,2LL)*Power(xj,2LL) + 

               1050LL*Power(r,3LL)*Power(xj,3LL) + 420LL*Power(r,4LL)*Power(xj,4LL) + 

               126LL*Power(r,5LL)*Power(xj,5LL) + 28LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL)) - 

            35LL*Power(xi,8LL)*Power(xj,12LL)*

             (-2904300LL + 4943925LL*r*xj + 258930LL*Power(r,2LL)*Power(xj,2LL) - 

               359520LL*Power(r,3LL)*Power(xj,3LL) - 70440LL*Power(r,4LL)*Power(xj,4LL) - 

               4176LL*Power(r,5LL)*Power(xj,5LL) + 32LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) + 

            35LL*Power(xi,12LL)*Power(xj,8LL)*

             (-49140LL - 98865LL*r*xj + 3510LL*Power(r,2LL)*Power(xj,2LL) - 

               131040LL*Power(r,3LL)*Power(xj,3LL) - 7800LL*Power(r,4LL)*Power(xj,4LL) + 

               3204LL*Power(r,5LL)*Power(xj,5LL) + 360LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) + 

            Power(xi,6LL)*Power(xj,14LL)*

             (446489820LL - 54796455LL*r*xj - 68983110LL*Power(r,2LL)*Power(xj,2LL) - 

               12782700LL*Power(r,3LL)*Power(xj,3LL) - 

               663600LL*Power(r,4LL)*Power(xj,4LL) + 53928LL*Power(r,5LL)*Power(xj,5LL) + 

               7728LL*Power(r,6LL)*Power(xj,6LL) + 256LL*Power(r,7LL)*Power(xj,7LL))))/

       (630LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,15LL)*Power(xi + xj,14LL)) - 

      (2520LL*exp(2LL*r*(xi + xj))*(xi + xj)*

          Power(Power(xi,2LL) - Power(xj,2LL),15LL) + 

         exp(2LL*r*xj)*Power(xj,10LL)*

          (-16380LL*Power(r,4LL)*Power(xi,25LL) - 1008LL*Power(r,5LL)*Power(xi,26LL) - 

            28LL*Power(r,6LL)*Power(xi,27LL) - 

            840LL*Power(r,5LL)*Power(xi,24LL)*Power(xj,2LL) + 2205LL*xi*Power(xj,20LL) + 

            3780LL*r*Power(xi,2LL)*Power(xj,20LL) + 

            2100LL*Power(r,2LL)*Power(xi,3LL)*Power(xj,20LL) - 

            1680LL*Power(r,3LL)*Power(xi,24LL)*(91LL + Power(r,2LL)*Power(xj,2LL)) + 

            525LL*Power(xi,3LL)*Power(xj,18LL)*(-63LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            42LL*Power(r,3LL)*Power(xi,23LL)*

             (-810LL*r*Power(xj,2LL) + 8LL*Power(r,3LL)*Power(xj,4LL)) + 

            63LL*r*Power(xi,5LL)*Power(xj,16LL)*

             (-500LL*r*Power(xj,2LL) + 8LL*Power(r,3LL)*Power(xj,4LL)) + 

            210LL*Power(xi,4LL)*Power(xj,16LL)*

             (-270LL*r*Power(xj,2LL) + 8LL*Power(r,3LL)*Power(xj,4LL)) + 

            252LL*Power(r,2LL)*Power(xi,22LL)*

             (-2450LL*r*Power(xj,2LL) + 68LL*Power(r,3LL)*Power(xj,4LL)) - 

            1260LL*r*Power(xi,17LL)*Power(xj,4LL)*

             (-20290LL*r*Power(xj,2LL) + 464LL*Power(r,3LL)*Power(xj,4LL)) + 

            21LL*r*Power(xi,9LL)*Power(xj,12LL)*

             (-36920LL*r*Power(xj,2LL) + 3312LL*Power(r,3LL)*Power(xj,4LL)) + 

            126LL*Power(r,2LL)*Power(xi,23LL)*

             (-6825LL - 405LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            63LL*Power(xi,5LL)*Power(xj,16LL)*

             (3675LL - 250LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            504LL*r*Power(xi,22LL)*(-5460LL - 1225LL*Power(r,2LL)*Power(xj,2LL) + 

               17LL*Power(r,4LL)*Power(xj,4LL)) - 

            1260LL*Power(xi,17LL)*Power(xj,4LL)*

             (141729LL - 10145LL*Power(r,2LL)*Power(xj,2LL) + 

               116LL*Power(r,4LL)*Power(xj,4LL)) + 

            21LL*Power(xi,9LL)*Power(xj,12LL)*

             (164775LL - 18460LL*Power(r,2LL)*Power(xj,2LL) + 

               828LL*Power(r,4LL)*Power(xj,4LL)) + 

            14LL*Power(xi,6LL)*Power(xj,14LL)*

             (28350LL*r*Power(xj,2LL) - 1800LL*Power(r,3LL)*Power(xj,4LL) + 

               12LL*Power(r,5LL)*Power(xj,6LL)) - 

            210LL*Power(xi,8LL)*Power(xj,12LL)*

             (8190LL*r*Power(xj,2LL) - 840LL*Power(r,3LL)*Power(xj,4LL) + 

               12LL*Power(r,5LL)*Power(xj,6LL)) + 

            42LL*Power(xi,10LL)*Power(xj,10LL)*

             (-5850LL*r*Power(xj,2LL) - 35360LL*Power(r,3LL)*Power(xj,4LL) + 

               24LL*Power(r,5LL)*Power(xj,6LL)) + 

            r*Power(xi,7LL)*Power(xj,14LL)*

             (220500LL*r*Power(xj,2LL) - 7560LL*Power(r,3LL)*Power(xj,4LL) + 

               24LL*Power(r,5LL)*Power(xj,6LL)) - 

            21LL*r*Power(xi,11LL)*Power(xj,10LL)*

             (-436800LL*r*Power(xj,2LL) + 2208LL*Power(r,3LL)*Power(xj,4LL) + 

               24LL*Power(r,5LL)*Power(xj,6LL)) + 

            280LL*Power(xi,18LL)*Power(xj,2LL)*

             (-147906LL*r*Power(xj,2LL) + 9480LL*Power(r,3LL)*Power(xj,4LL) + 

               24LL*Power(r,5LL)*Power(xj,6LL)) - 

            35LL*r*Power(xi,15LL)*Power(xj,6LL)*

             (719040LL*r*Power(xj,2LL) - 28080LL*Power(r,3LL)*Power(xj,4LL) + 

               48LL*Power(r,5LL)*Power(xj,6LL)) + 

            14LL*r*Power(xi,19LL)*Power(xj,2LL)*

             (253530LL*r*Power(xj,2LL) - 15408LL*Power(r,3LL)*Power(xj,4LL) + 

               120LL*Power(r,5LL)*Power(xj,6LL)) - 

            630LL*Power(xi,14LL)*Power(xj,6LL)*

             (-28770LL*r*Power(xj,2LL) - 9360LL*Power(r,3LL)*Power(xj,4LL) + 

               120LL*Power(r,5LL)*Power(xj,6LL)) + 

            210LL*Power(xi,16LL)*Power(xj,4LL)*

             (656982LL*r*Power(xj,2LL) - 46960LL*Power(r,3LL)*Power(xj,4LL) + 

               312LL*Power(r,5LL)*Power(xj,6LL)) - 

            84LL*Power(xi,20LL)*(472500LL*r*Power(xj,2LL) - 

               34980LL*Power(r,3LL)*Power(xj,4LL) + 552LL*Power(r,5LL)*Power(xj,6LL)) - 

            2LL*r*Power(xi,21LL)*(3197250LL*r*Power(xj,2LL) - 

               165564LL*Power(r,3LL)*Power(xj,4LL) + 768LL*Power(r,5LL)*Power(xj,6LL)) + 

            r*Power(xi,13LL)*Power(xj,8LL)*

             (-5569200LL*r*Power(xj,2LL) - 448560LL*Power(r,3LL)*Power(xj,4LL) + 

               1536LL*Power(r,5LL)*Power(xj,6LL)) + 

            14LL*Power(xi,12LL)*Power(xj,8LL)*

             (-5043870LL*r*Power(xj,2LL) + 78000LL*Power(r,3LL)*Power(xj,4LL) + 

               2064LL*Power(r,5LL)*Power(xj,6LL)) + 

            Power(xi,7LL)*Power(xj,14LL)*

             (-1003275LL + 110250LL*Power(r,2LL)*Power(xj,2LL) - 

               1890LL*Power(r,4LL)*Power(xj,4LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) - 

            21LL*Power(xi,11LL)*Power(xj,10LL)*

             (-1033695LL - 218400LL*Power(r,2LL)*Power(xj,2LL) + 

               552LL*Power(r,4LL)*Power(xj,4LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) - 

            35LL*Power(xi,15LL)*Power(xj,6LL)*

             (-1565613LL + 359520LL*Power(r,2LL)*Power(xj,2LL) - 

               7020LL*Power(r,4LL)*Power(xj,4LL) + 8LL*Power(r,6LL)*Power(xj,6LL)) + 

            14LL*Power(xi,19LL)*Power(xj,2LL)*

             (-4980150LL + 126765LL*Power(r,2LL)*Power(xj,2LL) - 

               3852LL*Power(r,4LL)*Power(xj,4LL) + 20LL*Power(r,6LL)*Power(xj,6LL)) - 

            2LL*Power(xi,21LL)*(1949220LL + 1598625LL*Power(r,2LL)*Power(xj,2LL) - 

               41391LL*Power(r,4LL)*Power(xj,4LL) + 128LL*Power(r,6LL)*Power(xj,6LL)) + 

            Power(xi,13LL)*Power(xj,8LL)*

             (173037375LL - 2784600LL*Power(r,2LL)*Power(xj,2LL) - 

               112140LL*Power(r,4LL)*Power(xj,4LL) + 256LL*Power(r,6LL)*Power(xj,6LL))) + 

         2LL*exp(2LL*r*xj)*Power(xj,11LL)*

          (-3276LL*Power(r,5LL)*Power(xi,25LL) - 168LL*Power(r,6LL)*Power(xi,26LL) - 

            4LL*Power(r,7LL)*Power(xi,27LL) + 1260LL*Power(xj,20LL) + 

            2205LL*r*xi*Power(xj,20LL) + 

            1890LL*Power(xi,2LL)*Power(xj,18LL)*(-10LL + Power(r,2LL)*Power(xj,2LL)) - 

            420LL*Power(r,4LL)*Power(xi,24LL)*(91LL + Power(r,2LL)*Power(xj,2LL)) + 

            525LL*r*Power(xi,3LL)*Power(xj,18LL)*(-63LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            42LL*Power(r,3LL)*Power(xi,23LL)*

             (-6825LL - 405LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            63LL*r*Power(xi,5LL)*Power(xj,16LL)*

             (3675LL - 250LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            210LL*Power(xi,4LL)*Power(xj,16LL)*

             (630LL - 135LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            252LL*Power(r,2LL)*Power(xi,22LL)*

             (-5460LL - 1225LL*Power(r,2LL)*Power(xj,2LL) + 17LL*Power(r,4LL)*Power(xj,4LL)) 

    - 1260LL*r*Power(xi,17LL)*Power(xj,4LL)*

             (141729LL - 10145LL*Power(r,2LL)*Power(xj,2LL) + 

               116LL*Power(r,4LL)*Power(xj,4LL)) + 

            21LL*r*Power(xi,9LL)*Power(xj,12LL)*

             (164775LL - 18460LL*Power(r,2LL)*Power(xj,2LL) + 

               828LL*Power(r,4LL)*Power(xj,4LL)) + 

            14LL*Power(xi,6LL)*Power(xj,14LL)*

             (-40950LL + 14175LL*Power(r,2LL)*Power(xj,2LL) - 

               450LL*Power(r,4LL)*Power(xj,4LL) + 2LL*Power(r,6LL)*Power(xj,6LL)) - 

            210LL*Power(xi,8LL)*Power(xj,12LL)*

             (-8190LL + 4095LL*Power(r,2LL)*Power(xj,2LL) - 

               210LL*Power(r,4LL)*Power(xj,4LL) + 2LL*Power(r,6LL)*Power(xj,6LL)) + 

            42LL*Power(xi,10LL)*Power(xj,10LL)*

             (-209430LL - 2925LL*Power(r,2LL)*Power(xj,2LL) - 

               8840LL*Power(r,4LL)*Power(xj,4LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) + 

            r*Power(xi,7LL)*Power(xj,14LL)*

             (-1003275LL + 110250LL*Power(r,2LL)*Power(xj,2LL) - 

               1890LL*Power(r,4LL)*Power(xj,4LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) - 

            21LL*r*Power(xi,11LL)*Power(xj,10LL)*

             (-1033695LL - 218400LL*Power(r,2LL)*Power(xj,2LL) + 

               552LL*Power(r,4LL)*Power(xj,4LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) + 

            280LL*Power(xi,18LL)*Power(xj,2LL)*

             (-385560LL - 73953LL*Power(r,2LL)*Power(xj,2LL) + 

               2370LL*Power(r,4LL)*Power(xj,4LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) - 

            35LL*r*Power(xi,15LL)*Power(xj,6LL)*

             (-1565613LL + 359520LL*Power(r,2LL)*Power(xj,2LL) - 

               7020LL*Power(r,4LL)*Power(xj,4LL) + 8LL*Power(r,6LL)*Power(xj,6LL)) + 

            14LL*r*Power(xi,19LL)*Power(xj,2LL)*

             (-4980150LL + 126765LL*Power(r,2LL)*Power(xj,2LL) - 

               3852LL*Power(r,4LL)*Power(xj,4LL) + 20LL*Power(r,6LL)*Power(xj,6LL)) - 

            630LL*Power(xi,14LL)*Power(xj,6LL)*

             (708714LL - 14385LL*Power(r,2LL)*Power(xj,2LL) - 

               2340LL*Power(r,4LL)*Power(xj,4LL) + 20LL*Power(r,6LL)*Power(xj,6LL)) + 

            210LL*Power(xi,16LL)*Power(xj,4LL)*

             (-2087532LL + 328491LL*Power(r,2LL)*Power(xj,2LL) - 

               11740LL*Power(r,4LL)*Power(xj,4LL) + 52LL*Power(r,6LL)*Power(xj,6LL)) - 

            84LL*Power(xi,20LL)*(59670LL + 236250LL*Power(r,2LL)*Power(xj,2LL) - 

               8745LL*Power(r,4LL)*Power(xj,4LL) + 92LL*Power(r,6LL)*Power(xj,6LL)) - 

            2LL*r*Power(xi,21LL)*(1949220LL + 1598625LL*Power(r,2LL)*Power(xj,2LL) - 

               41391LL*Power(r,4LL)*Power(xj,4LL) + 128LL*Power(r,6LL)*Power(xj,6LL)) + 

            r*Power(xi,13LL)*Power(xj,8LL)*

             (173037375LL - 2784600LL*Power(r,2LL)*Power(xj,2LL) - 

               112140LL*Power(r,4LL)*Power(xj,4LL) + 256LL*Power(r,6LL)*Power(xj,6LL)) + 

            14LL*Power(xi,12LL)*Power(xj,8LL)*

             (-7260750LL - 2521935LL*Power(r,2LL)*Power(xj,2LL) + 

               19500LL*Power(r,4LL)*Power(xj,4LL) + 344LL*Power(r,6LL)*Power(xj,6LL))) + 

         exp(2LL*r*xi)*Power(xi,10LL)*

          (210LL*Power(xi,2LL)*Power(xj,18LL)*

             (332010LL*xj + 189000LL*r*Power(xj,2LL) + 

               45675LL*Power(r,2LL)*Power(xj,3LL) + 5880LL*Power(r,3LL)*Power(xj,4LL) + 

               405LL*Power(r,4LL)*Power(xj,5LL) + 12LL*Power(r,5LL)*Power(xj,6LL)) + 

            105LL*Power(xi,18LL)*Power(xj,2LL)*

             (315LL*xj + 540LL*r*Power(xj,2LL) + 450LL*Power(r,2LL)*Power(xj,3LL) + 

               240LL*Power(r,3LL)*Power(xj,4LL) + 90LL*Power(r,4LL)*Power(xj,5LL) + 

               24LL*Power(r,5LL)*Power(xj,6LL)) - 

            1365LL*Power(xi,10LL)*Power(xj,10LL)*

             (15903LL*xj - 51732LL*r*Power(xj,2LL) - 6120LL*Power(r,2LL)*Power(xj,3LL) + 

               4320LL*Power(r,3LL)*Power(xj,4LL) + 900LL*Power(r,4LL)*Power(xj,5LL) + 

               48LL*Power(r,5LL)*Power(xj,6LL)) + 

            Power(xi,14LL)*Power(xj,6LL)*

             (1003275LL*xj + 1719900LL*r*Power(xj,2LL) + 

               1162980LL*Power(r,2LL)*Power(xj,3LL) + 

               1485120LL*Power(r,3LL)*Power(xj,4LL) + 57960LL*Power(r,4LL)*Power(xj,5LL) - 

               28896LL*Power(r,5LL)*Power(xj,6LL) - 1792LL*Power(r,6LL)*Power(xj,7LL)) + 

            2LL*Power(xj,20LL)*(1949220LL*xj + 1375920LL*r*Power(xj,2LL) + 

               429975LL*Power(r,2LL)*Power(xj,3LL) + 76440LL*Power(r,3LL)*Power(xj,4LL) + 

               8190LL*Power(r,4LL)*Power(xj,5LL) + 504LL*Power(r,5LL)*Power(xj,6LL) + 

               14LL*Power(r,6LL)*Power(xj,7LL)) - 

            42LL*Power(xi,4LL)*Power(xj,16LL)*

             (-4251870LL*xj - 986040LL*r*Power(xj,2LL) + 

               126765LL*Power(r,2LL)*Power(xj,3LL) + 69960LL*Power(r,3LL)*Power(xj,4LL) + 

               9855LL*Power(r,4LL)*Power(xj,5LL) + 612LL*Power(r,5LL)*Power(xj,6LL) + 

               14LL*Power(r,6LL)*Power(xj,7LL)) + 

            21LL*Power(xi,16LL)*Power(xj,4LL)*

             (-11025LL*xj - 18900LL*r*Power(xj,2LL) - 15750LL*Power(r,2LL)*Power(xj,3LL) - 

               8400LL*Power(r,3LL)*Power(xj,4LL) - 4140LL*Power(r,4LL)*Power(xj,5LL) - 

               48LL*Power(r,5LL)*Power(xj,6LL) + 28LL*Power(r,6LL)*Power(xj,7LL)) - 

            Power(xi,20LL)*(2205LL*xj + 3780LL*r*Power(xj,2LL) + 

               3150LL*Power(r,2LL)*Power(xj,3LL) + 1680LL*Power(r,3LL)*Power(xj,4LL) + 

               630LL*Power(r,4LL)*Power(xj,5LL) + 168LL*Power(r,5LL)*Power(xj,6LL) + 

               28LL*Power(r,6LL)*Power(xj,7LL)) - 

            35LL*Power(xi,8LL)*Power(xj,12LL)*

             (4943925LL*xj + 517860LL*r*Power(xj,2LL) - 

               1078560LL*Power(r,2LL)*Power(xj,3LL) - 

               281760LL*Power(r,3LL)*Power(xj,4LL) - 20880LL*Power(r,4LL)*Power(xj,5LL) + 

               192LL*Power(r,5LL)*Power(xj,6LL) + 56LL*Power(r,6LL)*Power(xj,7LL)) + 

            35LL*Power(xi,12LL)*Power(xj,8LL)*

             (-98865LL*xj + 7020LL*r*Power(xj,2LL) - 393120LL*Power(r,2LL)*Power(xj,3LL) - 

               31200LL*Power(r,3LL)*Power(xj,4LL) + 16020LL*Power(r,4LL)*Power(xj,5LL) + 

               2160LL*Power(r,5LL)*Power(xj,6LL) + 56LL*Power(r,6LL)*Power(xj,7LL)) + 

            Power(xi,6LL)*Power(xj,14LL)*

             (-54796455LL*xj - 137966220LL*r*Power(xj,2LL) - 

               38348100LL*Power(r,2LL)*Power(xj,3LL) - 

               2654400LL*Power(r,3LL)*Power(xj,4LL) + 269640LL*Power(r,4LL)*Power(xj,5LL) + 

               46368LL*Power(r,5LL)*Power(xj,6LL) + 1792LL*Power(r,6LL)*Power(xj,7LL))) + 

         2LL*exp(2LL*r*xi)*Power(xi,11LL)*

          (210LL*Power(xi,2LL)*Power(xj,18LL)*

             (514080LL + 332010LL*r*xj + 94500LL*Power(r,2LL)*Power(xj,2LL) + 

               15225LL*Power(r,3LL)*Power(xj,3LL) + 1470LL*Power(r,4LL)*Power(xj,4LL) + 

               81LL*Power(r,5LL)*Power(xj,5LL) + 2LL*Power(r,6LL)*Power(xj,6LL)) + 

            105LL*Power(xi,18LL)*Power(xj,2LL)*

             (180LL + 315LL*r*xj + 270LL*Power(r,2LL)*Power(xj,2LL) + 

               150LL*Power(r,3LL)*Power(xj,3LL) + 60LL*Power(r,4LL)*Power(xj,4LL) + 

               18LL*Power(r,5LL)*Power(xj,5LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) - 

            1365LL*Power(xi,10LL)*Power(xj,10LL)*

             (-6444LL + 15903LL*r*xj - 25866LL*Power(r,2LL)*Power(xj,2LL) - 

               2040LL*Power(r,3LL)*Power(xj,3LL) + 1080LL*Power(r,4LL)*Power(xj,4LL) + 

               180LL*Power(r,5LL)*Power(xj,5LL) + 8LL*Power(r,6LL)*Power(xj,6LL)) + 

            Power(xi,14LL)*Power(xj,6LL)*

             (573300LL + 1003275LL*r*xj + 859950LL*Power(r,2LL)*Power(xj,2LL) + 

               387660LL*Power(r,3LL)*Power(xj,3LL) + 371280LL*Power(r,4LL)*Power(xj,4LL) + 

               11592LL*Power(r,5LL)*Power(xj,5LL) - 4816LL*Power(r,6LL)*Power(xj,6LL) - 

               256LL*Power(r,7LL)*Power(xj,7LL)) + 

            2LL*Power(xj,20LL)*(2506140LL + 1949220LL*r*xj + 

               687960LL*Power(r,2LL)*Power(xj,2LL) + 143325LL*Power(r,3LL)*Power(xj,3LL) + 

               19110LL*Power(r,4LL)*Power(xj,4LL) + 1638LL*Power(r,5LL)*Power(xj,5LL) + 

               84LL*Power(r,6LL)*Power(xj,6LL) + 2LL*Power(r,7LL)*Power(xj,7LL)) - 

            42LL*Power(xi,4LL)*Power(xj,16LL)*

             (-10437660LL - 4251870LL*r*xj - 493020LL*Power(r,2LL)*Power(xj,2LL) + 

               42255LL*Power(r,3LL)*Power(xj,3LL) + 17490LL*Power(r,4LL)*Power(xj,4LL) + 

               1971LL*Power(r,5LL)*Power(xj,5LL) + 102LL*Power(r,6LL)*Power(xj,6LL) + 

               2LL*Power(r,7LL)*Power(xj,7LL)) + 

            21LL*Power(xi,16LL)*Power(xj,4LL)*

             (-6300LL - 11025LL*r*xj - 9450LL*Power(r,2LL)*Power(xj,2LL) - 

               5250LL*Power(r,3LL)*Power(xj,3LL) - 2100LL*Power(r,4LL)*Power(xj,4LL) - 

               828LL*Power(r,5LL)*Power(xj,5LL) - 8LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL)) - 

            Power(xi,20LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r,2LL)*Power(xj,2LL) + 

               1050LL*Power(r,3LL)*Power(xj,3LL) + 420LL*Power(r,4LL)*Power(xj,4LL) + 

               126LL*Power(r,5LL)*Power(xj,5LL) + 28LL*Power(r,6LL)*Power(xj,6LL) + 

               4LL*Power(r,7LL)*Power(xj,7LL)) - 

            35LL*Power(xi,8LL)*Power(xj,12LL)*

             (-2904300LL + 4943925LL*r*xj + 258930LL*Power(r,2LL)*Power(xj,2LL) - 

               359520LL*Power(r,3LL)*Power(xj,3LL) - 70440LL*Power(r,4LL)*Power(xj,4LL) - 

               4176LL*Power(r,5LL)*Power(xj,5LL) + 32LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) + 

            35LL*Power(xi,12LL)*Power(xj,8LL)*

             (-49140LL - 98865LL*r*xj + 3510LL*Power(r,2LL)*Power(xj,2LL) - 

               131040LL*Power(r,3LL)*Power(xj,3LL) - 7800LL*Power(r,4LL)*Power(xj,4LL) + 

               3204LL*Power(r,5LL)*Power(xj,5LL) + 360LL*Power(r,6LL)*Power(xj,6LL) + 

               8LL*Power(r,7LL)*Power(xj,7LL)) + 

            Power(xi,6LL)*Power(xj,14LL)*

             (446489820LL - 54796455LL*r*xj - 68983110LL*Power(r,2LL)*Power(xj,2LL) - 

               12782700LL*Power(r,3LL)*Power(xj,3LL) - 663600LL*Power(r,4LL)*Power(xj,4LL) + 

               53928LL*Power(r,5LL)*Power(xj,5LL) + 7728LL*Power(r,6LL)*Power(xj,6LL) + 

               256LL*Power(r,7LL)*Power(xj,7LL))))/

       (1260LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,15LL)*Power(xi + xj,15LL))

    ; }
   
  }
  return S;
}

cl_R DSlater_4S_5S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = 0LL

    ; } else {  S = -(-25913502934444125LL*xi + 28454994247680000LL*exp(2LL*r*xi)*xi - 

          46744023242416500LL*r*Power(xi,2LL) - 

          41723129607909750LL*Power(r,2LL)*Power(xi,3LL) - 

          24550942638222000LL*Power(r,3LL)*Power(xi,4LL) - 

          10704286944351000LL*Power(r,4LL)*Power(xi,5LL) - 

          3684699450432000LL*Power(r,5LL)*Power(xi,6LL) - 

          1041667066440000LL*Power(r,6LL)*Power(xi,7LL) - 

          248293113868800LL*Power(r,7LL)*Power(xi,8LL) - 

          50808078921600LL*Power(r,8LL)*Power(xi,9LL) - 

          9033331507200LL*Power(r,9LL)*Power(xi,10LL) - 

          1405184901120LL*Power(r,10LL)*Power(xi,11LL) - 

          191616122880LL*Power(r,11LL)*Power(xi,12LL) - 

          22811443200LL*Power(r,12LL)*Power(xi,13LL) - 

          2339635200LL*Power(r,13LL)*Power(xi,14LL) - 

          200540160LL*Power(r,14LL)*Power(xi,15LL) - 

          13369344LL*Power(r,15LL)*Power(xi,16LL) - 557056LL*Power(r,16LL)*Power(xi,17LL))/

       (1.422749712384e16*exp(2LL*r*xi)*r) + 

      (-14227497123840000LL + 14227497123840000LL*exp(2LL*r*xi) - 

         25913502934444125LL*r*xi - 23372011621208250LL*Power(r,2LL)*Power(xi,2LL) - 

         13907709869303250LL*Power(r,3LL)*Power(xi,3LL) - 

         6137735659555500LL*Power(r,4LL)*Power(xi,4LL) - 

         2140857388870200LL*Power(r,5LL)*Power(xi,5LL) - 

         614116575072000LL*Power(r,6LL)*Power(xi,6LL) - 

         148809580920000LL*Power(r,7LL)*Power(xi,7LL) - 

         31036639233600LL*Power(r,8LL)*Power(xi,8LL) - 

         5645342102400LL*Power(r,9LL)*Power(xi,9LL) - 

         903333150720LL*Power(r,10LL)*Power(xi,10LL) - 

         127744081920LL*Power(r,11LL)*Power(xi,11LL) - 

         15968010240LL*Power(r,12LL)*Power(xi,12LL) - 

         1754726400LL*Power(r,13LL)*Power(xi,13LL) - 

         167116800LL*Power(r,14LL)*Power(xi,14LL) - 

         13369344LL*Power(r,15LL)*Power(xi,15LL) - 835584LL*Power(r,16LL)*Power(xi,16LL) - 

         32768LL*Power(r,17LL)*Power(xi,17LL))/

       (1.422749712384e16*exp(2LL*r*xi)*Power(r,2LL)) + 

      (xi*(-14227497123840000LL + 14227497123840000LL*exp(2LL*r*xi) - 

           25913502934444125LL*r*xi - 23372011621208250LL*Power(r,2LL)*Power(xi,2LL) - 

           13907709869303250LL*Power(r,3LL)*Power(xi,3LL) - 

           6137735659555500LL*Power(r,4LL)*Power(xi,4LL) - 

           2140857388870200LL*Power(r,5LL)*Power(xi,5LL) - 

           614116575072000LL*Power(r,6LL)*Power(xi,6LL) - 

           148809580920000LL*Power(r,7LL)*Power(xi,7LL) - 

           31036639233600LL*Power(r,8LL)*Power(xi,8LL) - 

           5645342102400LL*Power(r,9LL)*Power(xi,9LL) - 

           903333150720LL*Power(r,10LL)*Power(xi,10LL) - 

           127744081920LL*Power(r,11LL)*Power(xi,11LL) - 

           15968010240LL*Power(r,12LL)*Power(xi,12LL) - 

           1754726400LL*Power(r,13LL)*Power(xi,13LL) - 

           167116800LL*Power(r,14LL)*Power(xi,14LL) - 

           13369344LL*Power(r,15LL)*Power(xi,15LL) - 835584LL*Power(r,16LL)*Power(xi,16LL) - 

           32768LL*Power(r,17LL)*Power(xi,17LL)))/(7.11374856192e15*exp(2LL*r*xi)*r)

    ; }
 
  }
  else {
      if (r == 0LL) {  S = 0LL

    ; } else {  S = (56700LL*exp(2LL*r*(xi + xj))*Power(Power(xi,2LL) - Power(xj,2LL),17LL) + 

         9LL*exp(2LL*r*xj)*Power(xj,12LL)*

          (-980LL*Power(r,6LL)*Power(xi,28LL) - 20LL*Power(r,7LL)*Power(xi,29LL) + 

            6300LL*Power(xj,22LL) + 11025LL*r*xi*Power(xj,22LL) - 

            50LL*Power(r,5LL)*Power(xi,27LL)*(441LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            3150LL*Power(xi,2LL)*Power(xj,20LL)*(-34LL + 3LL*Power(r,2LL)*Power(xj,2LL)) + 

            525LL*r*Power(xi,3LL)*Power(xj,20LL)*

             (-357LL + 10LL*Power(r,2LL)*Power(xj,2LL)) - 

            420LL*Power(r,4LL)*Power(xi,26LL)*(700LL + 19LL*Power(r,2LL)*Power(xj,2LL)) + 

            1050LL*Power(xi,4LL)*Power(xj,18LL)*

             (816LL - 153LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            210LL*r*Power(xi,5LL)*Power(xj,18LL)*

             (7140LL - 425LL*Power(r,2LL)*Power(xj,2LL) + 3LL*Power(r,4LL)*Power(xj,4LL)) + 

            42LL*Power(r,3LL)*Power(xi,25LL)*

             (-59500LL - 6035LL*Power(r,2LL)*Power(xj,2LL) + 

               18LL*Power(r,4LL)*Power(xj,4LL)) + 

            84LL*Power(r,2LL)*Power(xi,24LL)*

             (-160650LL - 52700LL*Power(r,2LL)*Power(xj,2LL) + 

               397LL*Power(r,4LL)*Power(xj,4LL)) - 

            28LL*Power(xi,12LL)*Power(xj,10LL)*

             (100849950LL + 27100125LL*Power(r,2LL)*Power(xj,2LL) + 

               186150LL*Power(r,4LL)*Power(xj,4LL) - 2177LL*Power(r,6LL)*Power(xj,6LL)) + 

            140LL*Power(xi,6LL)*Power(xj,16LL)*

             (-30600LL + 9180LL*Power(r,2LL)*Power(xj,2LL) - 

               255LL*Power(r,4LL)*Power(xj,4LL) + Power(r,6LL)*Power(xj,6LL)) - 

            2380LL*Power(xi,8LL)*Power(xj,14LL)*

             (-6300LL + 2700LL*Power(r,2LL)*Power(xj,2LL) - 

               120LL*Power(r,4LL)*Power(xj,4LL) + Power(r,6LL)*Power(xj,6LL)) + 

            10LL*r*Power(xi,7LL)*Power(xj,16LL)*

             (-749700LL + 71400LL*Power(r,2LL)*Power(xj,2LL) - 

               1071LL*Power(r,4LL)*Power(xj,4LL) + 2LL*Power(r,6LL)*Power(xj,6LL)) + 

            204LL*r*Power(xi,15LL)*Power(xj,8LL)*

             (28962255LL - 1744750LL*Power(r,2LL)*Power(xj,2LL) + 

               9555LL*Power(r,4LL)*Power(xj,4LL) + 6LL*Power(r,6LL)*Power(xj,6LL)) - 

            42LL*r*Power(xi,11LL)*Power(xj,12LL)*

             (-12911925LL - 1634550LL*Power(r,2LL)*Power(xj,2LL) - 

               7103LL*Power(r,4LL)*Power(xj,4LL) + 18LL*Power(r,6LL)*Power(xj,6LL)) + 

            2LL*r*Power(xi,9LL)*Power(xj,14LL)*

             (16948575LL - 1184400LL*Power(r,2LL)*Power(xj,2LL) + 

               63861LL*Power(r,4LL)*Power(xj,4LL) + 50LL*Power(r,6LL)*Power(xj,6LL)) + 

            28LL*Power(xi,22LL)*(-2180250LL - 10993050LL*Power(r,2LL)*Power(xj,2LL) + 

               14925LL*Power(r,4LL)*Power(xj,4LL) + 73LL*Power(r,6LL)*Power(xj,6LL)) - 

            952LL*Power(xi,14LL)*Power(xj,8LL)*

             (16966215LL + 725175LL*Power(r,2LL)*Power(xj,2LL) - 

               36075LL*Power(r,4LL)*Power(xj,4LL) + 79LL*Power(r,6LL)*Power(xj,6LL)) - 

            84LL*Power(xi,10LL)*Power(xj,12LL)*

             (1723800LL + 279225LL*Power(r,2LL)*Power(xj,2LL) + 

               45600LL*Power(r,4LL)*Power(xj,4LL) + 107LL*Power(r,6LL)*Power(xj,6LL)) - 

            35LL*r*Power(xi,17LL)*Power(xj,6LL)*

             (132637869LL - 2205240LL*Power(r,2LL)*Power(xj,2LL) - 

               48348LL*Power(r,4LL)*Power(xj,4LL) + 136LL*Power(r,6LL)*Power(xj,6LL)) - 

            6LL*r*Power(xi,21LL)*Power(xj,2LL)*

             (192298050LL + 12644275LL*Power(r,2LL)*Power(xj,2LL) - 

               218029LL*Power(r,4LL)*Power(xj,4LL) + 204LL*Power(r,6LL)*Power(xj,6LL)) + 

            4LL*r*Power(xi,13LL)*Power(xj,10LL)*

             (1259522775LL + 15895425LL*Power(r,2LL)*Power(xj,2LL) - 

               493017LL*Power(r,4LL)*Power(xj,4LL) + 263LL*Power(r,6LL)*Power(xj,6LL)) - 

            140LL*Power(xi,16LL)*Power(xj,6LL)*

             (180826281LL - 15101406LL*Power(r,2LL)*Power(xj,2LL) + 

               160140LL*Power(r,4LL)*Power(xj,4LL) + 442LL*Power(r,6LL)*Power(xj,6LL)) - 

            2LL*r*Power(xi,23LL)*(21366450LL + 23526300LL*Power(r,2LL)*Power(xj,2LL) - 

               246729LL*Power(r,4LL)*Power(xj,4LL) + 526LL*Power(r,6LL)*Power(xj,6LL)) + 

            7LL*r*Power(xi,19LL)*Power(xj,4LL)*

             (-811081215LL + 39095550LL*Power(r,2LL)*Power(xj,2LL) - 

               515916LL*Power(r,4LL)*Power(xj,4LL) + 680LL*Power(r,6LL)*Power(xj,6LL)) + 

            70LL*Power(xi,18LL)*Power(xj,4LL)*

             (-180554454LL + 9873711LL*Power(r,2LL)*Power(xj,2LL) - 

               414120LL*Power(r,4LL)*Power(xj,4LL) + 2924LL*Power(r,6LL)*Power(xj,6LL)) - 

            14LL*Power(xi,20LL)*Power(xj,2LL)*

             (136919700LL + 71867115LL*Power(r,2LL)*Power(xj,2LL) - 

               2154150LL*Power(r,4LL)*Power(xj,4LL) + 10268LL*Power(r,6LL)*Power(xj,6LL))) 

    - 4LL*exp(2LL*r*xi)*Power(xi,10LL)*

          (-10710LL*Power(xi,12LL)*Power(xj,12LL)*

             (-3555LL - 127008LL*r*xj + 138384LL*Power(r,2LL)*Power(xj,2LL) - 

               74556LL*Power(r,3LL)*Power(xj,3LL) - 22284LL*Power(r,4LL)*Power(xj,4LL) + 

               408LL*Power(r,5LL)*Power(xj,5LL) + 576LL*Power(r,6LL)*Power(xj,6LL) + 

               60LL*Power(r,7LL)*Power(xj,7LL) + 2LL*Power(r,8LL)*Power(xj,8LL)) + 

            2LL*Power(xi,20LL)*Power(xj,4LL)*

             (963900LL + 1735020LL*r*xj + 1542240LL*Power(r,2LL)*Power(xj,2LL) + 

               899640LL*Power(r,3LL)*Power(xj,3LL) + 385560LL*Power(r,4LL)*Power(xj,4LL) + 

               128520LL*Power(r,5LL)*Power(xj,5LL) + 34272LL*Power(r,6LL)*Power(xj,6LL) + 

               9126LL*Power(r,7LL)*Power(xj,7LL) + 333LL*Power(r,8LL)*Power(xj,8LL) - 

               20LL*Power(r,9LL)*Power(xj,9LL)) - 

            2LL*Power(xj,24LL)*(119041650LL + 107137485LL*r*xj + 

               45110520LL*Power(r,2LL)*Power(xj,2LL) + 

               11695320LL*Power(r,3LL)*Power(xj,3LL) + 

               2063880LL*Power(r,4LL)*Power(xj,4LL) + 

               257985LL*Power(r,5LL)*Power(xj,5LL) + 22932LL*Power(r,6LL)*Power(xj,6LL) + 

               1404LL*Power(r,7LL)*Power(xj,7LL) + 54LL*Power(r,8LL)*Power(xj,8LL) + 

               Power(r,9LL)*Power(xj,9LL)) + 

            2LL*Power(xi,2LL)*Power(xj,22LL)*

             (-3264488325LL - 2505368880LL*r*xj - 

               881390160LL*Power(r,2LL)*Power(xj,2LL) - 

               185775660LL*Power(r,3LL)*Power(xj,3LL) - 

               25639740LL*Power(r,4LL)*Power(xj,4LL) - 

               2361555LL*Power(r,5LL)*Power(xj,5LL) - 

               139356LL*Power(r,6LL)*Power(xj,6LL) - 4482LL*Power(r,7LL)*Power(xj,7LL) - 

               27LL*Power(r,8LL)*Power(xj,8LL) + 2LL*Power(r,9LL)*Power(xj,9LL)) + 

            Power(xi,24LL)*(14175LL + 25515LL*r*xj + 22680LL*Power(r,2LL)*Power(xj,2LL) + 

               13230LL*Power(r,3LL)*Power(xj,3LL) + 5670LL*Power(r,4LL)*Power(xj,4LL) + 

               1890LL*Power(r,5LL)*Power(xj,5LL) + 504LL*Power(r,6LL)*Power(xj,6LL) + 

               108LL*Power(r,7LL)*Power(xj,7LL) + 18LL*Power(r,8LL)*Power(xj,8LL) + 

               2LL*Power(r,9LL)*Power(xj,9LL)) - 

            102LL*Power(xi,10LL)*Power(xj,14LL)*

             (44986725LL - 97433280LL*r*xj + 44467920LL*Power(r,2LL)*Power(xj,2LL) + 

               15857100LL*Power(r,3LL)*Power(xj,3LL) - 

               457380LL*Power(r,4LL)*Power(xj,4LL) - 620550LL*Power(r,5LL)*Power(xj,5LL) - 

               83160LL*Power(r,6LL)*Power(xj,6LL) - 4068LL*Power(r,7LL)*Power(xj,7LL) - 

               6LL*Power(r,8LL)*Power(xj,8LL) + 4LL*Power(r,9LL)*Power(xj,9LL)) + 

            102LL*Power(xi,14LL)*Power(xj,10LL)*

             (-859950LL - 1437345LL*r*xj - 2260440LL*Power(r,2LL)*Power(xj,2LL) + 

               810810LL*Power(r,3LL)*Power(xj,3LL) - 

               1056510LL*Power(r,4LL)*Power(xj,4LL) - 

               217854LL*Power(r,5LL)*Power(xj,5LL) + 6552LL*Power(r,6LL)*Power(xj,6LL) + 

               3852LL*Power(r,7LL)*Power(xj,7LL) + 258LL*Power(r,8LL)*Power(xj,8LL) + 

               4LL*Power(r,9LL)*Power(xj,9LL)) - 

            Power(xi,22LL)*Power(xj,2LL)*

             (240975LL + 433755LL*r*xj + 385560LL*Power(r,2LL)*Power(xj,2LL) + 

               224910LL*Power(r,3LL)*Power(xj,3LL) + 96390LL*Power(r,4LL)*Power(xj,4LL) + 

               32130LL*Power(r,5LL)*Power(xj,5LL) + 8568LL*Power(r,6LL)*Power(xj,6LL) + 

               1836LL*Power(r,7LL)*Power(xj,7LL) + 306LL*Power(r,8LL)*Power(xj,8LL) + 

               4LL*Power(r,9LL)*Power(xj,9LL)) + 

            2LL*Power(xi,4LL)*Power(xj,20LL)*

             (-18032978565LL - 9823683240LL*r*xj - 

               2047323600LL*Power(r,2LL)*Power(xj,2LL) - 

               129098340LL*Power(r,3LL)*Power(xj,3LL) + 

               26410860LL*Power(r,4LL)*Power(xj,4LL) + 

               7094304LL*Power(r,5LL)*Power(xj,5LL) + 

               788256LL*Power(r,6LL)*Power(xj,6LL) + 48654LL*Power(r,7LL)*Power(xj,7LL) + 

               1593LL*Power(r,8LL)*Power(xj,8LL) + 20LL*Power(r,9LL)*Power(xj,9LL)) - 

            6LL*Power(xi,16LL)*Power(xj,8LL)*

             (-5622750LL - 10120950LL*r*xj - 8996400LL*Power(r,2LL)*Power(xj,2LL) - 

               5698350LL*Power(r,3LL)*Power(xj,3LL) - 

               897750LL*Power(r,4LL)*Power(xj,4LL) - 

               1641591LL*Power(r,5LL)*Power(xj,5LL) - 

               211932LL*Power(r,6LL)*Power(xj,6LL) + 10224LL*Power(r,7LL)*Power(xj,7LL) + 

               2364LL*Power(r,8LL)*Power(xj,8LL) + 73LL*Power(r,9LL)*Power(xj,9LL)) + 

            2LL*Power(xi,18LL)*Power(xj,6LL)*

             (-4819500LL - 8675100LL*r*xj - 7711200LL*Power(r,2LL)*Power(xj,2LL) - 

               4498200LL*Power(r,3LL)*Power(xj,3LL) - 

               1927800LL*Power(r,4LL)*Power(xj,4LL) - 

               561519LL*Power(r,5LL)*Power(xj,5LL) - 279468LL*Power(r,6LL)*Power(xj,6LL) - 

               20682LL*Power(r,7LL)*Power(xj,7LL) + 1305LL*Power(r,8LL)*Power(xj,8LL) + 

               106LL*Power(r,9LL)*Power(xj,9LL)) + 

            3LL*Power(xi,8LL)*Power(xj,16LL)*

             (-9364244085LL + 6940428705LL*r*xj + 

               2117684520LL*Power(r,2LL)*Power(xj,2LL) - 

               230268150LL*Power(r,3LL)*Power(xj,3LL) - 

               149610510LL*Power(r,4LL)*Power(xj,4LL) - 

               21824334LL*Power(r,5LL)*Power(xj,5LL) - 

               1223208LL*Power(r,6LL)*Power(xj,6LL) + 12708LL*Power(r,7LL)*Power(xj,7LL) + 

               4470LL*Power(r,8LL)*Power(xj,8LL) + 146LL*Power(r,9LL)*Power(xj,9LL)) - 

            Power(xi,6LL)*Power(xj,18LL)*

             (57304872765LL + 7147185255LL*r*xj - 

               5801702760LL*Power(r,2LL)*Power(xj,2LL) - 

               2053388610LL*Power(r,3LL)*Power(xj,3LL) - 

               271655370LL*Power(r,4LL)*Power(xj,4LL) - 

               10864854LL*Power(r,5LL)*Power(xj,5LL) + 

               1337112LL*Power(r,6LL)*Power(xj,6LL) + 202716LL*Power(r,7LL)*Power(xj,7LL) + 

               10746LL*Power(r,8LL)*Power(xj,8LL) + 212LL*Power(r,9LL)*Power(xj,9LL))))/

       (56700LL*exp(2LL*r*(xi + xj))*Power(r,2LL)*Power(xi - xj,17LL)*

         Power(xi + xj,17LL)) + (56700LL*exp(2LL*r*(xi + xj))*

          Power(Power(xi,2LL) - Power(xj,2LL),17LL) + 

         9LL*exp(2LL*r*xj)*Power(xj,12LL)*

          (-980LL*Power(r,6LL)*Power(xi,28LL) - 20LL*Power(r,7LL)*Power(xi,29LL) + 

            6300LL*Power(xj,22LL) + 11025LL*r*xi*Power(xj,22LL) - 

            50LL*Power(r,5LL)*Power(xi,27LL)*(441LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            3150LL*Power(xi,2LL)*Power(xj,20LL)*(-34LL + 3LL*Power(r,2LL)*Power(xj,2LL)) + 

            525LL*r*Power(xi,3LL)*Power(xj,20LL)*

             (-357LL + 10LL*Power(r,2LL)*Power(xj,2LL)) - 

            420LL*Power(r,4LL)*Power(xi,26LL)*(700LL + 19LL*Power(r,2LL)*Power(xj,2LL)) + 

            1050LL*Power(xi,4LL)*Power(xj,18LL)*

             (816LL - 153LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            210LL*r*Power(xi,5LL)*Power(xj,18LL)*

             (7140LL - 425LL*Power(r,2LL)*Power(xj,2LL) + 3LL*Power(r,4LL)*Power(xj,4LL)) + 

            42LL*Power(r,3LL)*Power(xi,25LL)*

             (-59500LL - 6035LL*Power(r,2LL)*Power(xj,2LL) + 

               18LL*Power(r,4LL)*Power(xj,4LL)) + 

            84LL*Power(r,2LL)*Power(xi,24LL)*

             (-160650LL - 52700LL*Power(r,2LL)*Power(xj,2LL) + 

               397LL*Power(r,4LL)*Power(xj,4LL)) - 

            28LL*Power(xi,12LL)*Power(xj,10LL)*

             (100849950LL + 27100125LL*Power(r,2LL)*Power(xj,2LL) + 

               186150LL*Power(r,4LL)*Power(xj,4LL) - 2177LL*Power(r,6LL)*Power(xj,6LL)) + 

            140LL*Power(xi,6LL)*Power(xj,16LL)*

             (-30600LL + 9180LL*Power(r,2LL)*Power(xj,2LL) - 

               255LL*Power(r,4LL)*Power(xj,4LL) + Power(r,6LL)*Power(xj,6LL)) - 

            2380LL*Power(xi,8LL)*Power(xj,14LL)*

             (-6300LL + 2700LL*Power(r,2LL)*Power(xj,2LL) - 

               120LL*Power(r,4LL)*Power(xj,4LL) + Power(r,6LL)*Power(xj,6LL)) + 

            10LL*r*Power(xi,7LL)*Power(xj,16LL)*

             (-749700LL + 71400LL*Power(r,2LL)*Power(xj,2LL) - 

               1071LL*Power(r,4LL)*Power(xj,4LL) + 2LL*Power(r,6LL)*Power(xj,6LL)) + 

            204LL*r*Power(xi,15LL)*Power(xj,8LL)*

             (28962255LL - 1744750LL*Power(r,2LL)*Power(xj,2LL) + 

               9555LL*Power(r,4LL)*Power(xj,4LL) + 6LL*Power(r,6LL)*Power(xj,6LL)) - 

            42LL*r*Power(xi,11LL)*Power(xj,12LL)*

             (-12911925LL - 1634550LL*Power(r,2LL)*Power(xj,2LL) - 

               7103LL*Power(r,4LL)*Power(xj,4LL) + 18LL*Power(r,6LL)*Power(xj,6LL)) + 

            2LL*r*Power(xi,9LL)*Power(xj,14LL)*

             (16948575LL - 1184400LL*Power(r,2LL)*Power(xj,2LL) + 

               63861LL*Power(r,4LL)*Power(xj,4LL) + 50LL*Power(r,6LL)*Power(xj,6LL)) + 

            28LL*Power(xi,22LL)*(-2180250LL - 10993050LL*Power(r,2LL)*Power(xj,2LL) + 

               14925LL*Power(r,4LL)*Power(xj,4LL) + 73LL*Power(r,6LL)*Power(xj,6LL)) - 

            952LL*Power(xi,14LL)*Power(xj,8LL)*

             (16966215LL + 725175LL*Power(r,2LL)*Power(xj,2LL) - 

               36075LL*Power(r,4LL)*Power(xj,4LL) + 79LL*Power(r,6LL)*Power(xj,6LL)) - 

            84LL*Power(xi,10LL)*Power(xj,12LL)*

             (1723800LL + 279225LL*Power(r,2LL)*Power(xj,2LL) + 

               45600LL*Power(r,4LL)*Power(xj,4LL) + 107LL*Power(r,6LL)*Power(xj,6LL)) - 

            35LL*r*Power(xi,17LL)*Power(xj,6LL)*

             (132637869LL - 2205240LL*Power(r,2LL)*Power(xj,2LL) - 

               48348LL*Power(r,4LL)*Power(xj,4LL) + 136LL*Power(r,6LL)*Power(xj,6LL)) - 

            6LL*r*Power(xi,21LL)*Power(xj,2LL)*

             (192298050LL + 12644275LL*Power(r,2LL)*Power(xj,2LL) - 

               218029LL*Power(r,4LL)*Power(xj,4LL) + 204LL*Power(r,6LL)*Power(xj,6LL)) + 

            4LL*r*Power(xi,13LL)*Power(xj,10LL)*

             (1259522775LL + 15895425LL*Power(r,2LL)*Power(xj,2LL) - 

               493017LL*Power(r,4LL)*Power(xj,4LL) + 263LL*Power(r,6LL)*Power(xj,6LL)) - 

            140LL*Power(xi,16LL)*Power(xj,6LL)*

             (180826281LL - 15101406LL*Power(r,2LL)*Power(xj,2LL) + 

               160140LL*Power(r,4LL)*Power(xj,4LL) + 442LL*Power(r,6LL)*Power(xj,6LL)) - 

            2LL*r*Power(xi,23LL)*(21366450LL + 23526300LL*Power(r,2LL)*Power(xj,2LL) - 

               246729LL*Power(r,4LL)*Power(xj,4LL) + 526LL*Power(r,6LL)*Power(xj,6LL)) + 

            7LL*r*Power(xi,19LL)*Power(xj,4LL)*

             (-811081215LL + 39095550LL*Power(r,2LL)*Power(xj,2LL) - 

               515916LL*Power(r,4LL)*Power(xj,4LL) + 680LL*Power(r,6LL)*Power(xj,6LL)) + 

            70LL*Power(xi,18LL)*Power(xj,4LL)*

             (-180554454LL + 9873711LL*Power(r,2LL)*Power(xj,2LL) - 

               414120LL*Power(r,4LL)*Power(xj,4LL) + 2924LL*Power(r,6LL)*Power(xj,6LL)) - 

            14LL*Power(xi,20LL)*Power(xj,2LL)*

             (136919700LL + 71867115LL*Power(r,2LL)*Power(xj,2LL) - 

               2154150LL*Power(r,4LL)*Power(xj,4LL) + 10268LL*Power(r,6LL)*Power(xj,6LL))) 

    - 4LL*exp(2LL*r*xi)*Power(xi,10LL)*

          (-10710LL*Power(xi,12LL)*Power(xj,12LL)*

             (-3555LL - 127008LL*r*xj + 138384LL*Power(r,2LL)*Power(xj,2LL) - 

               74556LL*Power(r,3LL)*Power(xj,3LL) - 22284LL*Power(r,4LL)*Power(xj,4LL) + 

               408LL*Power(r,5LL)*Power(xj,5LL) + 576LL*Power(r,6LL)*Power(xj,6LL) + 

               60LL*Power(r,7LL)*Power(xj,7LL) + 2LL*Power(r,8LL)*Power(xj,8LL)) + 

            2LL*Power(xi,20LL)*Power(xj,4LL)*

             (963900LL + 1735020LL*r*xj + 1542240LL*Power(r,2LL)*Power(xj,2LL) + 

               899640LL*Power(r,3LL)*Power(xj,3LL) + 385560LL*Power(r,4LL)*Power(xj,4LL) + 

               128520LL*Power(r,5LL)*Power(xj,5LL) + 34272LL*Power(r,6LL)*Power(xj,6LL) + 

               9126LL*Power(r,7LL)*Power(xj,7LL) + 333LL*Power(r,8LL)*Power(xj,8LL) - 

               20LL*Power(r,9LL)*Power(xj,9LL)) - 

            2LL*Power(xj,24LL)*(119041650LL + 107137485LL*r*xj + 

               45110520LL*Power(r,2LL)*Power(xj,2LL) + 

               11695320LL*Power(r,3LL)*Power(xj,3LL) + 

               2063880LL*Power(r,4LL)*Power(xj,4LL) + 

               257985LL*Power(r,5LL)*Power(xj,5LL) + 22932LL*Power(r,6LL)*Power(xj,6LL) + 

               1404LL*Power(r,7LL)*Power(xj,7LL) + 54LL*Power(r,8LL)*Power(xj,8LL) + 

               Power(r,9LL)*Power(xj,9LL)) + 

            2LL*Power(xi,2LL)*Power(xj,22LL)*

             (-3264488325LL - 2505368880LL*r*xj - 

               881390160LL*Power(r,2LL)*Power(xj,2LL) - 

               185775660LL*Power(r,3LL)*Power(xj,3LL) - 

               25639740LL*Power(r,4LL)*Power(xj,4LL) - 

               2361555LL*Power(r,5LL)*Power(xj,5LL) - 

               139356LL*Power(r,6LL)*Power(xj,6LL) - 4482LL*Power(r,7LL)*Power(xj,7LL) - 

               27LL*Power(r,8LL)*Power(xj,8LL) + 2LL*Power(r,9LL)*Power(xj,9LL)) + 

            Power(xi,24LL)*(14175LL + 25515LL*r*xj + 22680LL*Power(r,2LL)*Power(xj,2LL) + 

               13230LL*Power(r,3LL)*Power(xj,3LL) + 5670LL*Power(r,4LL)*Power(xj,4LL) + 

               1890LL*Power(r,5LL)*Power(xj,5LL) + 504LL*Power(r,6LL)*Power(xj,6LL) + 

               108LL*Power(r,7LL)*Power(xj,7LL) + 18LL*Power(r,8LL)*Power(xj,8LL) + 

               2LL*Power(r,9LL)*Power(xj,9LL)) - 

            102LL*Power(xi,10LL)*Power(xj,14LL)*

             (44986725LL - 97433280LL*r*xj + 44467920LL*Power(r,2LL)*Power(xj,2LL) + 

               15857100LL*Power(r,3LL)*Power(xj,3LL) - 

               457380LL*Power(r,4LL)*Power(xj,4LL) - 620550LL*Power(r,5LL)*Power(xj,5LL) - 

               83160LL*Power(r,6LL)*Power(xj,6LL) - 4068LL*Power(r,7LL)*Power(xj,7LL) - 

               6LL*Power(r,8LL)*Power(xj,8LL) + 4LL*Power(r,9LL)*Power(xj,9LL)) + 

            102LL*Power(xi,14LL)*Power(xj,10LL)*

             (-859950LL - 1437345LL*r*xj - 2260440LL*Power(r,2LL)*Power(xj,2LL) + 

               810810LL*Power(r,3LL)*Power(xj,3LL) - 

               1056510LL*Power(r,4LL)*Power(xj,4LL) - 

               217854LL*Power(r,5LL)*Power(xj,5LL) + 6552LL*Power(r,6LL)*Power(xj,6LL) + 

               3852LL*Power(r,7LL)*Power(xj,7LL) + 258LL*Power(r,8LL)*Power(xj,8LL) + 

               4LL*Power(r,9LL)*Power(xj,9LL)) - 

            Power(xi,22LL)*Power(xj,2LL)*

             (240975LL + 433755LL*r*xj + 385560LL*Power(r,2LL)*Power(xj,2LL) + 

               224910LL*Power(r,3LL)*Power(xj,3LL) + 96390LL*Power(r,4LL)*Power(xj,4LL) + 

               32130LL*Power(r,5LL)*Power(xj,5LL) + 8568LL*Power(r,6LL)*Power(xj,6LL) + 

               1836LL*Power(r,7LL)*Power(xj,7LL) + 306LL*Power(r,8LL)*Power(xj,8LL) + 

               4LL*Power(r,9LL)*Power(xj,9LL)) + 

            2LL*Power(xi,4LL)*Power(xj,20LL)*

             (-18032978565LL - 9823683240LL*r*xj - 

               2047323600LL*Power(r,2LL)*Power(xj,2LL) - 

               129098340LL*Power(r,3LL)*Power(xj,3LL) + 

               26410860LL*Power(r,4LL)*Power(xj,4LL) + 

               7094304LL*Power(r,5LL)*Power(xj,5LL) + 

               788256LL*Power(r,6LL)*Power(xj,6LL) + 48654LL*Power(r,7LL)*Power(xj,7LL) + 

               1593LL*Power(r,8LL)*Power(xj,8LL) + 20LL*Power(r,9LL)*Power(xj,9LL)) - 

            6LL*Power(xi,16LL)*Power(xj,8LL)*

             (-5622750LL - 10120950LL*r*xj - 8996400LL*Power(r,2LL)*Power(xj,2LL) - 

               5698350LL*Power(r,3LL)*Power(xj,3LL) - 

               897750LL*Power(r,4LL)*Power(xj,4LL) - 

               1641591LL*Power(r,5LL)*Power(xj,5LL) - 

               211932LL*Power(r,6LL)*Power(xj,6LL) + 10224LL*Power(r,7LL)*Power(xj,7LL) + 

               2364LL*Power(r,8LL)*Power(xj,8LL) + 73LL*Power(r,9LL)*Power(xj,9LL)) + 

            2LL*Power(xi,18LL)*Power(xj,6LL)*

             (-4819500LL - 8675100LL*r*xj - 7711200LL*Power(r,2LL)*Power(xj,2LL) - 

               4498200LL*Power(r,3LL)*Power(xj,3LL) - 

               1927800LL*Power(r,4LL)*Power(xj,4LL) - 

               561519LL*Power(r,5LL)*Power(xj,5LL) - 279468LL*Power(r,6LL)*Power(xj,6LL) - 

               20682LL*Power(r,7LL)*Power(xj,7LL) + 1305LL*Power(r,8LL)*Power(xj,8LL) + 

               106LL*Power(r,9LL)*Power(xj,9LL)) + 

            3LL*Power(xi,8LL)*Power(xj,16LL)*

             (-9364244085LL + 6940428705LL*r*xj + 

               2117684520LL*Power(r,2LL)*Power(xj,2LL) - 

               230268150LL*Power(r,3LL)*Power(xj,3LL) - 

               149610510LL*Power(r,4LL)*Power(xj,4LL) - 

               21824334LL*Power(r,5LL)*Power(xj,5LL) - 

               1223208LL*Power(r,6LL)*Power(xj,6LL) + 12708LL*Power(r,7LL)*Power(xj,7LL) + 

               4470LL*Power(r,8LL)*Power(xj,8LL) + 146LL*Power(r,9LL)*Power(xj,9LL)) - 

            Power(xi,6LL)*Power(xj,18LL)*

             (57304872765LL + 7147185255LL*r*xj - 

               5801702760LL*Power(r,2LL)*Power(xj,2LL) - 

               2053388610LL*Power(r,3LL)*Power(xj,3LL) - 

               271655370LL*Power(r,4LL)*Power(xj,4LL) - 

               10864854LL*Power(r,5LL)*Power(xj,5LL) + 

               1337112LL*Power(r,6LL)*Power(xj,6LL) + 202716LL*Power(r,7LL)*Power(xj,7LL) + 

               10746LL*Power(r,8LL)*Power(xj,8LL) + 212LL*Power(r,9LL)*Power(xj,9LL))))/

       (28350LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,17LL)*Power(xi + xj,16LL)) - 

      (113400LL*exp(2LL*r*(xi + xj))*(xi + xj)*

          Power(Power(xi,2LL) - Power(xj,2LL),17LL) + 

         9LL*exp(2LL*r*xj)*Power(xj,12LL)*

          (-5880LL*Power(r,5LL)*Power(xi,28LL) - 140LL*Power(r,6LL)*Power(xi,29LL) - 

            15960LL*Power(r,5LL)*Power(xi,26LL)*Power(xj,2LL) - 

            200LL*Power(r,6LL)*Power(xi,27LL)*Power(xj,2LL) + 11025LL*xi*Power(xj,22LL) + 

            18900LL*r*Power(xi,2LL)*Power(xj,22LL) + 

            10500LL*Power(r,2LL)*Power(xi,3LL)*Power(xj,22LL) - 

            250LL*Power(r,4LL)*Power(xi,27LL)*(441LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            525LL*Power(xi,3LL)*Power(xj,20LL)*(-357LL + 10LL*Power(r,2LL)*Power(xj,2LL)) - 

            1680LL*Power(r,3LL)*Power(xi,26LL)*(700LL + 19LL*Power(r,2LL)*Power(xj,2LL)) + 

            1050LL*Power(xi,4LL)*Power(xj,18LL)*

             (-306LL*r*Power(xj,2LL) + 8LL*Power(r,3LL)*Power(xj,4LL)) + 

            210LL*r*Power(xi,5LL)*Power(xj,18LL)*

             (-850LL*r*Power(xj,2LL) + 12LL*Power(r,3LL)*Power(xj,4LL)) + 

            42LL*Power(r,3LL)*Power(xi,25LL)*

             (-12070LL*r*Power(xj,2LL) + 72LL*Power(r,3LL)*Power(xj,4LL)) + 

            84LL*Power(r,2LL)*Power(xi,24LL)*

             (-105400LL*r*Power(xj,2LL) + 1588LL*Power(r,3LL)*Power(xj,4LL)) + 

            210LL*Power(xi,5LL)*Power(xj,18LL)*

             (7140LL - 425LL*Power(r,2LL)*Power(xj,2LL) + 3LL*Power(r,4LL)*Power(xj,4LL)) + 

            126LL*Power(r,2LL)*Power(xi,25LL)*

             (-59500LL - 6035LL*Power(r,2LL)*Power(xj,2LL) + 18LL*Power(r,4LL)*Power(xj,4LL)) 

    + 168LL*r*Power(xi,24LL)*(-160650LL - 52700LL*Power(r,2LL)*Power(xj,2LL) + 

               397LL*Power(r,4LL)*Power(xj,4LL)) - 

            28LL*Power(xi,12LL)*Power(xj,10LL)*

             (54200250LL*r*Power(xj,2LL) + 744600LL*Power(r,3LL)*Power(xj,4LL) - 

               13062LL*Power(r,5LL)*Power(xj,6LL)) + 

            140LL*Power(xi,6LL)*Power(xj,16LL)*

             (18360LL*r*Power(xj,2LL) - 1020LL*Power(r,3LL)*Power(xj,4LL) + 

               6LL*Power(r,5LL)*Power(xj,6LL)) - 

            2380LL*Power(xi,8LL)*Power(xj,14LL)*

             (5400LL*r*Power(xj,2LL) - 480LL*Power(r,3LL)*Power(xj,4LL) + 

               6LL*Power(r,5LL)*Power(xj,6LL)) + 

            10LL*r*Power(xi,7LL)*Power(xj,16LL)*

             (142800LL*r*Power(xj,2LL) - 4284LL*Power(r,3LL)*Power(xj,4LL) + 

               12LL*Power(r,5LL)*Power(xj,6LL)) + 

            204LL*r*Power(xi,15LL)*Power(xj,8LL)*

             (-3489500LL*r*Power(xj,2LL) + 38220LL*Power(r,3LL)*Power(xj,4LL) + 

               36LL*Power(r,5LL)*Power(xj,6LL)) - 

            42LL*r*Power(xi,11LL)*Power(xj,12LL)*

             (-3269100LL*r*Power(xj,2LL) - 28412LL*Power(r,3LL)*Power(xj,4LL) + 

               108LL*Power(r,5LL)*Power(xj,6LL)) + 

            2LL*r*Power(xi,9LL)*Power(xj,14LL)*

             (-2368800LL*r*Power(xj,2LL) + 255444LL*Power(r,3LL)*Power(xj,4LL) + 

               300LL*Power(r,5LL)*Power(xj,6LL)) + 

            28LL*Power(xi,22LL)*(-21986100LL*r*Power(xj,2LL) + 

               59700LL*Power(r,3LL)*Power(xj,4LL) + 438LL*Power(r,5LL)*Power(xj,6LL)) - 

            952LL*Power(xi,14LL)*Power(xj,8LL)*

             (1450350LL*r*Power(xj,2LL) - 144300LL*Power(r,3LL)*Power(xj,4LL) + 

               474LL*Power(r,5LL)*Power(xj,6LL)) - 

            84LL*Power(xi,10LL)*Power(xj,12LL)*

             (558450LL*r*Power(xj,2LL) + 182400LL*Power(r,3LL)*Power(xj,4LL) + 

               642LL*Power(r,5LL)*Power(xj,6LL)) - 

            35LL*r*Power(xi,17LL)*Power(xj,6LL)*

             (-4410480LL*r*Power(xj,2LL) - 193392LL*Power(r,3LL)*Power(xj,4LL) + 

               816LL*Power(r,5LL)*Power(xj,6LL)) - 

            6LL*r*Power(xi,21LL)*Power(xj,2LL)*

             (25288550LL*r*Power(xj,2LL) - 872116LL*Power(r,3LL)*Power(xj,4LL) + 

               1224LL*Power(r,5LL)*Power(xj,6LL)) + 

            4LL*r*Power(xi,13LL)*Power(xj,10LL)*

             (31790850LL*r*Power(xj,2LL) - 1972068LL*Power(r,3LL)*Power(xj,4LL) + 

               1578LL*Power(r,5LL)*Power(xj,6LL)) - 

            140LL*Power(xi,16LL)*Power(xj,6LL)*

             (-30202812LL*r*Power(xj,2LL) + 640560LL*Power(r,3LL)*Power(xj,4LL) + 

               2652LL*Power(r,5LL)*Power(xj,6LL)) - 

            2LL*r*Power(xi,23LL)*(47052600LL*r*Power(xj,2LL) - 

               986916LL*Power(r,3LL)*Power(xj,4LL) + 3156LL*Power(r,5LL)*Power(xj,6LL)) + 

            7LL*r*Power(xi,19LL)*Power(xj,4LL)*

             (78191100LL*r*Power(xj,2LL) - 2063664LL*Power(r,3LL)*Power(xj,4LL) + 

               4080LL*Power(r,5LL)*Power(xj,6LL)) + 

            70LL*Power(xi,18LL)*Power(xj,4LL)*

             (19747422LL*r*Power(xj,2LL) - 1656480LL*Power(r,3LL)*Power(xj,4LL) + 

               17544LL*Power(r,5LL)*Power(xj,6LL)) - 

            14LL*Power(xi,20LL)*Power(xj,2LL)*

             (143734230LL*r*Power(xj,2LL) - 8616600LL*Power(r,3LL)*Power(xj,4LL) + 

               61608LL*Power(r,5LL)*Power(xj,6LL)) + 

            10LL*Power(xi,7LL)*Power(xj,16LL)*

             (-749700LL + 71400LL*Power(r,2LL)*Power(xj,2LL) - 

               1071LL*Power(r,4LL)*Power(xj,4LL) + 2LL*Power(r,6LL)*Power(xj,6LL)) + 

            204LL*Power(xi,15LL)*Power(xj,8LL)*

             (28962255LL - 1744750LL*Power(r,2LL)*Power(xj,2LL) + 

               9555LL*Power(r,4LL)*Power(xj,4LL) + 6LL*Power(r,6LL)*Power(xj,6LL)) - 

            42LL*Power(xi,11LL)*Power(xj,12LL)*

             (-12911925LL - 1634550LL*Power(r,2LL)*Power(xj,2LL) - 

               7103LL*Power(r,4LL)*Power(xj,4LL) + 18LL*Power(r,6LL)*Power(xj,6LL)) + 

            2LL*Power(xi,9LL)*Power(xj,14LL)*

             (16948575LL - 1184400LL*Power(r,2LL)*Power(xj,2LL) + 

               63861LL*Power(r,4LL)*Power(xj,4LL) + 50LL*Power(r,6LL)*Power(xj,6LL)) - 

            35LL*Power(xi,17LL)*Power(xj,6LL)*

             (132637869LL - 2205240LL*Power(r,2LL)*Power(xj,2LL) - 

               48348LL*Power(r,4LL)*Power(xj,4LL) + 136LL*Power(r,6LL)*Power(xj,6LL)) - 

            6LL*Power(xi,21LL)*Power(xj,2LL)*

             (192298050LL + 12644275LL*Power(r,2LL)*Power(xj,2LL) - 

               218029LL*Power(r,4LL)*Power(xj,4LL) + 204LL*Power(r,6LL)*Power(xj,6LL)) + 

            4LL*Power(xi,13LL)*Power(xj,10LL)*

             (1259522775LL + 15895425LL*Power(r,2LL)*Power(xj,2LL) - 

               493017LL*Power(r,4LL)*Power(xj,4LL) + 263LL*Power(r,6LL)*Power(xj,6LL)) - 

            2LL*Power(xi,23LL)*(21366450LL + 23526300LL*Power(r,2LL)*Power(xj,2LL) - 

               246729LL*Power(r,4LL)*Power(xj,4LL) + 526LL*Power(r,6LL)*Power(xj,6LL)) + 

            7LL*Power(xi,19LL)*Power(xj,4LL)*

             (-811081215LL + 39095550LL*Power(r,2LL)*Power(xj,2LL) - 

               515916LL*Power(r,4LL)*Power(xj,4LL) + 680LL*Power(r,6LL)*Power(xj,6LL))) + 

         18LL*exp(2LL*r*xj)*Power(xj,13LL)*

          (-980LL*Power(r,6LL)*Power(xi,28LL) - 20LL*Power(r,7LL)*Power(xi,29LL) + 

            6300LL*Power(xj,22LL) + 11025LL*r*xi*Power(xj,22LL) - 

            50LL*Power(r,5LL)*Power(xi,27LL)*(441LL + 2LL*Power(r,2LL)*Power(xj,2LL)) + 

            3150LL*Power(xi,2LL)*Power(xj,20LL)*(-34LL + 3LL*Power(r,2LL)*Power(xj,2LL)) + 

            525LL*r*Power(xi,3LL)*Power(xj,20LL)*(-357LL + 10LL*Power(r,2LL)*Power(xj,2LL)) - 

            420LL*Power(r,4LL)*Power(xi,26LL)*(700LL + 19LL*Power(r,2LL)*Power(xj,2LL)) + 

            1050LL*Power(xi,4LL)*Power(xj,18LL)*

             (816LL - 153LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            210LL*r*Power(xi,5LL)*Power(xj,18LL)*

             (7140LL - 425LL*Power(r,2LL)*Power(xj,2LL) + 3LL*Power(r,4LL)*Power(xj,4LL)) + 

            42LL*Power(r,3LL)*Power(xi,25LL)*

             (-59500LL - 6035LL*Power(r,2LL)*Power(xj,2LL) + 18LL*Power(r,4LL)*Power(xj,4LL)) 

    + 84LL*Power(r,2LL)*Power(xi,24LL)*(-160650LL - 52700LL*Power(r,2LL)*Power(xj,2LL) + 

               397LL*Power(r,4LL)*Power(xj,4LL)) - 

            28LL*Power(xi,12LL)*Power(xj,10LL)*

             (100849950LL + 27100125LL*Power(r,2LL)*Power(xj,2LL) + 

               186150LL*Power(r,4LL)*Power(xj,4LL) - 2177LL*Power(r,6LL)*Power(xj,6LL)) + 

            140LL*Power(xi,6LL)*Power(xj,16LL)*

             (-30600LL + 9180LL*Power(r,2LL)*Power(xj,2LL) - 

               255LL*Power(r,4LL)*Power(xj,4LL) + Power(r,6LL)*Power(xj,6LL)) - 

            2380LL*Power(xi,8LL)*Power(xj,14LL)*

             (-6300LL + 2700LL*Power(r,2LL)*Power(xj,2LL) - 

               120LL*Power(r,4LL)*Power(xj,4LL) + Power(r,6LL)*Power(xj,6LL)) + 

            10LL*r*Power(xi,7LL)*Power(xj,16LL)*

             (-749700LL + 71400LL*Power(r,2LL)*Power(xj,2LL) - 

               1071LL*Power(r,4LL)*Power(xj,4LL) + 2LL*Power(r,6LL)*Power(xj,6LL)) + 

            204LL*r*Power(xi,15LL)*Power(xj,8LL)*

             (28962255LL - 1744750LL*Power(r,2LL)*Power(xj,2LL) + 

               9555LL*Power(r,4LL)*Power(xj,4LL) + 6LL*Power(r,6LL)*Power(xj,6LL)) - 

            42LL*r*Power(xi,11LL)*Power(xj,12LL)*

             (-12911925LL - 1634550LL*Power(r,2LL)*Power(xj,2LL) - 

               7103LL*Power(r,4LL)*Power(xj,4LL) + 18LL*Power(r,6LL)*Power(xj,6LL)) + 

            2LL*r*Power(xi,9LL)*Power(xj,14LL)*

             (16948575LL - 1184400LL*Power(r,2LL)*Power(xj,2LL) + 

               63861LL*Power(r,4LL)*Power(xj,4LL) + 50LL*Power(r,6LL)*Power(xj,6LL)) + 

            28LL*Power(xi,22LL)*(-2180250LL - 10993050LL*Power(r,2LL)*Power(xj,2LL) + 

               14925LL*Power(r,4LL)*Power(xj,4LL) + 73LL*Power(r,6LL)*Power(xj,6LL)) - 

            952LL*Power(xi,14LL)*Power(xj,8LL)*

             (16966215LL + 725175LL*Power(r,2LL)*Power(xj,2LL) - 

               36075LL*Power(r,4LL)*Power(xj,4LL) + 79LL*Power(r,6LL)*Power(xj,6LL)) - 

            84LL*Power(xi,10LL)*Power(xj,12LL)*

             (1723800LL + 279225LL*Power(r,2LL)*Power(xj,2LL) + 

               45600LL*Power(r,4LL)*Power(xj,4LL) + 107LL*Power(r,6LL)*Power(xj,6LL)) - 

            35LL*r*Power(xi,17LL)*Power(xj,6LL)*

             (132637869LL - 2205240LL*Power(r,2LL)*Power(xj,2LL) - 

               48348LL*Power(r,4LL)*Power(xj,4LL) + 136LL*Power(r,6LL)*Power(xj,6LL)) - 

            6LL*r*Power(xi,21LL)*Power(xj,2LL)*

             (192298050LL + 12644275LL*Power(r,2LL)*Power(xj,2LL) - 

               218029LL*Power(r,4LL)*Power(xj,4LL) + 204LL*Power(r,6LL)*Power(xj,6LL)) + 

            4LL*r*Power(xi,13LL)*Power(xj,10LL)*

             (1259522775LL + 15895425LL*Power(r,2LL)*Power(xj,2LL) - 

               493017LL*Power(r,4LL)*Power(xj,4LL) + 263LL*Power(r,6LL)*Power(xj,6LL)) - 

            140LL*Power(xi,16LL)*Power(xj,6LL)*

             (180826281LL - 15101406LL*Power(r,2LL)*Power(xj,2LL) + 

               160140LL*Power(r,4LL)*Power(xj,4LL) + 442LL*Power(r,6LL)*Power(xj,6LL)) - 

            2LL*r*Power(xi,23LL)*(21366450LL + 23526300LL*Power(r,2LL)*Power(xj,2LL) - 

               246729LL*Power(r,4LL)*Power(xj,4LL) + 526LL*Power(r,6LL)*Power(xj,6LL)) + 

            7LL*r*Power(xi,19LL)*Power(xj,4LL)*

             (-811081215LL + 39095550LL*Power(r,2LL)*Power(xj,2LL) - 

               515916LL*Power(r,4LL)*Power(xj,4LL) + 680LL*Power(r,6LL)*Power(xj,6LL)) + 

            70LL*Power(xi,18LL)*Power(xj,4LL)*

             (-180554454LL + 9873711LL*Power(r,2LL)*Power(xj,2LL) - 

               414120LL*Power(r,4LL)*Power(xj,4LL) + 2924LL*Power(r,6LL)*Power(xj,6LL)) - 

            14LL*Power(xi,20LL)*Power(xj,2LL)*

             (136919700LL + 71867115LL*Power(r,2LL)*Power(xj,2LL) - 

               2154150LL*Power(r,4LL)*Power(xj,4LL) + 10268LL*Power(r,6LL)*Power(xj,6LL))) - 

         4LL*exp(2LL*r*xi)*Power(xi,10LL)*

          (-10710LL*Power(xi,12LL)*Power(xj,12LL)*

             (-127008LL*xj + 276768LL*r*Power(xj,2LL) - 

               223668LL*Power(r,2LL)*Power(xj,3LL) - 89136LL*Power(r,3LL)*Power(xj,4LL) + 

               2040LL*Power(r,4LL)*Power(xj,5LL) + 3456LL*Power(r,5LL)*Power(xj,6LL) + 

               420LL*Power(r,6LL)*Power(xj,7LL) + 16LL*Power(r,7LL)*Power(xj,8LL)) + 

            2LL*Power(xi,20LL)*Power(xj,4LL)*

             (1735020LL*xj + 3084480LL*r*Power(xj,2LL) + 

               2698920LL*Power(r,2LL)*Power(xj,3LL) + 

               1542240LL*Power(r,3LL)*Power(xj,4LL) + 

               642600LL*Power(r,4LL)*Power(xj,5LL) + 205632LL*Power(r,5LL)*Power(xj,6LL) + 

               63882LL*Power(r,6LL)*Power(xj,7LL) + 2664LL*Power(r,7LL)*Power(xj,8LL) - 

               180LL*Power(r,8LL)*Power(xj,9LL)) - 

            2LL*Power(xj,24LL)*(107137485LL*xj + 90221040LL*r*Power(xj,2LL) + 

               35085960LL*Power(r,2LL)*Power(xj,3LL) + 

               8255520LL*Power(r,3LL)*Power(xj,4LL) + 

               1289925LL*Power(r,4LL)*Power(xj,5LL) + 

               137592LL*Power(r,5LL)*Power(xj,6LL) + 9828LL*Power(r,6LL)*Power(xj,7LL) + 

               432LL*Power(r,7LL)*Power(xj,8LL) + 9LL*Power(r,8LL)*Power(xj,9LL)) + 

            2LL*Power(xi,2LL)*Power(xj,22LL)*

             (-2505368880LL*xj - 1762780320LL*r*Power(xj,2LL) - 

               557326980LL*Power(r,2LL)*Power(xj,3LL) - 

               102558960LL*Power(r,3LL)*Power(xj,4LL) - 

               11807775LL*Power(r,4LL)*Power(xj,5LL) - 

               836136LL*Power(r,5LL)*Power(xj,6LL) - 31374LL*Power(r,6LL)*Power(xj,7LL) - 

               216LL*Power(r,7LL)*Power(xj,8LL) + 18LL*Power(r,8LL)*Power(xj,9LL)) + 

            Power(xi,24LL)*(25515LL*xj + 45360LL*r*Power(xj,2LL) + 

               39690LL*Power(r,2LL)*Power(xj,3LL) + 22680LL*Power(r,3LL)*Power(xj,4LL) + 

               9450LL*Power(r,4LL)*Power(xj,5LL) + 3024LL*Power(r,5LL)*Power(xj,6LL) + 

               756LL*Power(r,6LL)*Power(xj,7LL) + 144LL*Power(r,7LL)*Power(xj,8LL) + 

               18LL*Power(r,8LL)*Power(xj,9LL)) - 

            102LL*Power(xi,10LL)*Power(xj,14LL)*

             (-97433280LL*xj + 88935840LL*r*Power(xj,2LL) + 

               47571300LL*Power(r,2LL)*Power(xj,3LL) - 

               1829520LL*Power(r,3LL)*Power(xj,4LL) - 

               3102750LL*Power(r,4LL)*Power(xj,5LL) - 

               498960LL*Power(r,5LL)*Power(xj,6LL) - 28476LL*Power(r,6LL)*Power(xj,7LL) - 

               48LL*Power(r,7LL)*Power(xj,8LL) + 36LL*Power(r,8LL)*Power(xj,9LL)) + 

            102LL*Power(xi,14LL)*Power(xj,10LL)*

             (-1437345LL*xj - 4520880LL*r*Power(xj,2LL) + 

               2432430LL*Power(r,2LL)*Power(xj,3LL) - 

               4226040LL*Power(r,3LL)*Power(xj,4LL) - 

               1089270LL*Power(r,4LL)*Power(xj,5LL) + 39312LL*Power(r,5LL)*Power(xj,6LL) + 

               26964LL*Power(r,6LL)*Power(xj,7LL) + 2064LL*Power(r,7LL)*Power(xj,8LL) + 

               36LL*Power(r,8LL)*Power(xj,9LL)) - 

            Power(xi,22LL)*Power(xj,2LL)*

             (433755LL*xj + 771120LL*r*Power(xj,2LL) + 

               674730LL*Power(r,2LL)*Power(xj,3LL) + 385560LL*Power(r,3LL)*Power(xj,4LL) + 

               160650LL*Power(r,4LL)*Power(xj,5LL) + 51408LL*Power(r,5LL)*Power(xj,6LL) + 

               12852LL*Power(r,6LL)*Power(xj,7LL) + 2448LL*Power(r,7LL)*Power(xj,8LL) + 

               36LL*Power(r,8LL)*Power(xj,9LL)) + 

            2LL*Power(xi,4LL)*Power(xj,20LL)*

             (-9823683240LL*xj - 4094647200LL*r*Power(xj,2LL) - 

               387295020LL*Power(r,2LL)*Power(xj,3LL) + 

               105643440LL*Power(r,3LL)*Power(xj,4LL) + 

               35471520LL*Power(r,4LL)*Power(xj,5LL) + 

               4729536LL*Power(r,5LL)*Power(xj,6LL) + 

               340578LL*Power(r,6LL)*Power(xj,7LL) + 12744LL*Power(r,7LL)*Power(xj,8LL) + 

               180LL*Power(r,8LL)*Power(xj,9LL)) - 

            6LL*Power(xi,16LL)*Power(xj,8LL)*

             (-10120950LL*xj - 17992800LL*r*Power(xj,2LL) - 

               17095050LL*Power(r,2LL)*Power(xj,3LL) - 

               3591000LL*Power(r,3LL)*Power(xj,4LL) - 

               8207955LL*Power(r,4LL)*Power(xj,5LL) - 

               1271592LL*Power(r,5LL)*Power(xj,6LL) + 71568LL*Power(r,6LL)*Power(xj,7LL) + 

               18912LL*Power(r,7LL)*Power(xj,8LL) + 657LL*Power(r,8LL)*Power(xj,9LL)) + 

            2LL*Power(xi,18LL)*Power(xj,6LL)*

             (-8675100LL*xj - 15422400LL*r*Power(xj,2LL) - 

               13494600LL*Power(r,2LL)*Power(xj,3LL) - 

               7711200LL*Power(r,3LL)*Power(xj,4LL) - 

               2807595LL*Power(r,4LL)*Power(xj,5LL) - 

               1676808LL*Power(r,5LL)*Power(xj,6LL) - 

               144774LL*Power(r,6LL)*Power(xj,7LL) + 10440LL*Power(r,7LL)*Power(xj,8LL) + 

               954LL*Power(r,8LL)*Power(xj,9LL)) + 

            3LL*Power(xi,8LL)*Power(xj,16LL)*

             (6940428705LL*xj + 4235369040LL*r*Power(xj,2LL) - 

               690804450LL*Power(r,2LL)*Power(xj,3LL) - 

               598442040LL*Power(r,3LL)*Power(xj,4LL) - 

               109121670LL*Power(r,4LL)*Power(xj,5LL) - 

               7339248LL*Power(r,5LL)*Power(xj,6LL) + 88956LL*Power(r,6LL)*Power(xj,7LL) + 

               35760LL*Power(r,7LL)*Power(xj,8LL) + 1314LL*Power(r,8LL)*Power(xj,9LL)) - 

            Power(xi,6LL)*Power(xj,18LL)*

             (7147185255LL*xj - 11603405520LL*r*Power(xj,2LL) - 

               6160165830LL*Power(r,2LL)*Power(xj,3LL) - 

               1086621480LL*Power(r,3LL)*Power(xj,4LL) - 

               54324270LL*Power(r,4LL)*Power(xj,5LL) + 

               8022672LL*Power(r,5LL)*Power(xj,6LL) + 

               1419012LL*Power(r,6LL)*Power(xj,7LL) + 85968LL*Power(r,7LL)*Power(xj,8LL) + 

               1908LL*Power(r,8LL)*Power(xj,9LL))) - 

         8LL*exp(2LL*r*xi)*Power(xi,11LL)*

          (-10710LL*Power(xi,12LL)*Power(xj,12LL)*

             (-3555LL - 127008LL*r*xj + 138384LL*Power(r,2LL)*Power(xj,2LL) - 

               74556LL*Power(r,3LL)*Power(xj,3LL) - 22284LL*Power(r,4LL)*Power(xj,4LL) + 

               408LL*Power(r,5LL)*Power(xj,5LL) + 576LL*Power(r,6LL)*Power(xj,6LL) + 

               60LL*Power(r,7LL)*Power(xj,7LL) + 2LL*Power(r,8LL)*Power(xj,8LL)) + 

            2LL*Power(xi,20LL)*Power(xj,4LL)*

             (963900LL + 1735020LL*r*xj + 1542240LL*Power(r,2LL)*Power(xj,2LL) + 

               899640LL*Power(r,3LL)*Power(xj,3LL) + 385560LL*Power(r,4LL)*Power(xj,4LL) + 

               128520LL*Power(r,5LL)*Power(xj,5LL) + 34272LL*Power(r,6LL)*Power(xj,6LL) + 

               9126LL*Power(r,7LL)*Power(xj,7LL) + 333LL*Power(r,8LL)*Power(xj,8LL) - 

               20LL*Power(r,9LL)*Power(xj,9LL)) - 

            2LL*Power(xj,24LL)*(119041650LL + 107137485LL*r*xj + 

               45110520LL*Power(r,2LL)*Power(xj,2LL) + 

               11695320LL*Power(r,3LL)*Power(xj,3LL) + 

               2063880LL*Power(r,4LL)*Power(xj,4LL) + 257985LL*Power(r,5LL)*Power(xj,5LL) + 

               22932LL*Power(r,6LL)*Power(xj,6LL) + 1404LL*Power(r,7LL)*Power(xj,7LL) + 

               54LL*Power(r,8LL)*Power(xj,8LL) + Power(r,9LL)*Power(xj,9LL)) + 

            2LL*Power(xi,2LL)*Power(xj,22LL)*

             (-3264488325LL - 2505368880LL*r*xj - 

               881390160LL*Power(r,2LL)*Power(xj,2LL) - 

               185775660LL*Power(r,3LL)*Power(xj,3LL) - 

               25639740LL*Power(r,4LL)*Power(xj,4LL) - 

               2361555LL*Power(r,5LL)*Power(xj,5LL) - 139356LL*Power(r,6LL)*Power(xj,6LL) - 

               4482LL*Power(r,7LL)*Power(xj,7LL) - 27LL*Power(r,8LL)*Power(xj,8LL) + 

               2LL*Power(r,9LL)*Power(xj,9LL)) + 

            Power(xi,24LL)*(14175LL + 25515LL*r*xj + 22680LL*Power(r,2LL)*Power(xj,2LL) + 

               13230LL*Power(r,3LL)*Power(xj,3LL) + 5670LL*Power(r,4LL)*Power(xj,4LL) + 

               1890LL*Power(r,5LL)*Power(xj,5LL) + 504LL*Power(r,6LL)*Power(xj,6LL) + 

               108LL*Power(r,7LL)*Power(xj,7LL) + 18LL*Power(r,8LL)*Power(xj,8LL) + 

               2LL*Power(r,9LL)*Power(xj,9LL)) - 

            102LL*Power(xi,10LL)*Power(xj,14LL)*

             (44986725LL - 97433280LL*r*xj + 44467920LL*Power(r,2LL)*Power(xj,2LL) + 

               15857100LL*Power(r,3LL)*Power(xj,3LL) - 

               457380LL*Power(r,4LL)*Power(xj,4LL) - 620550LL*Power(r,5LL)*Power(xj,5LL) - 

               83160LL*Power(r,6LL)*Power(xj,6LL) - 4068LL*Power(r,7LL)*Power(xj,7LL) - 

               6LL*Power(r,8LL)*Power(xj,8LL) + 4LL*Power(r,9LL)*Power(xj,9LL)) + 

            102LL*Power(xi,14LL)*Power(xj,10LL)*

             (-859950LL - 1437345LL*r*xj - 2260440LL*Power(r,2LL)*Power(xj,2LL) + 

               810810LL*Power(r,3LL)*Power(xj,3LL) - 1056510LL*Power(r,4LL)*Power(xj,4LL) - 

               217854LL*Power(r,5LL)*Power(xj,5LL) + 6552LL*Power(r,6LL)*Power(xj,6LL) + 

               3852LL*Power(r,7LL)*Power(xj,7LL) + 258LL*Power(r,8LL)*Power(xj,8LL) + 

               4LL*Power(r,9LL)*Power(xj,9LL)) - 

            Power(xi,22LL)*Power(xj,2LL)*

             (240975LL + 433755LL*r*xj + 385560LL*Power(r,2LL)*Power(xj,2LL) + 

               224910LL*Power(r,3LL)*Power(xj,3LL) + 96390LL*Power(r,4LL)*Power(xj,4LL) + 

               32130LL*Power(r,5LL)*Power(xj,5LL) + 8568LL*Power(r,6LL)*Power(xj,6LL) + 

               1836LL*Power(r,7LL)*Power(xj,7LL) + 306LL*Power(r,8LL)*Power(xj,8LL) + 

               4LL*Power(r,9LL)*Power(xj,9LL)) + 

            2LL*Power(xi,4LL)*Power(xj,20LL)*

             (-18032978565LL - 9823683240LL*r*xj - 

               2047323600LL*Power(r,2LL)*Power(xj,2LL) - 

               129098340LL*Power(r,3LL)*Power(xj,3LL) + 

               26410860LL*Power(r,4LL)*Power(xj,4LL) + 

               7094304LL*Power(r,5LL)*Power(xj,5LL) + 788256LL*Power(r,6LL)*Power(xj,6LL) + 

               48654LL*Power(r,7LL)*Power(xj,7LL) + 1593LL*Power(r,8LL)*Power(xj,8LL) + 

               20LL*Power(r,9LL)*Power(xj,9LL)) - 

            6LL*Power(xi,16LL)*Power(xj,8LL)*

             (-5622750LL - 10120950LL*r*xj - 8996400LL*Power(r,2LL)*Power(xj,2LL) - 

               5698350LL*Power(r,3LL)*Power(xj,3LL) - 897750LL*Power(r,4LL)*Power(xj,4LL) - 

               1641591LL*Power(r,5LL)*Power(xj,5LL) - 211932LL*Power(r,6LL)*Power(xj,6LL) + 

               10224LL*Power(r,7LL)*Power(xj,7LL) + 2364LL*Power(r,8LL)*Power(xj,8LL) + 

               73LL*Power(r,9LL)*Power(xj,9LL)) + 

            2LL*Power(xi,18LL)*Power(xj,6LL)*

             (-4819500LL - 8675100LL*r*xj - 7711200LL*Power(r,2LL)*Power(xj,2LL) - 

               4498200LL*Power(r,3LL)*Power(xj,3LL) - 

               1927800LL*Power(r,4LL)*Power(xj,4LL) - 561519LL*Power(r,5LL)*Power(xj,5LL) - 

               279468LL*Power(r,6LL)*Power(xj,6LL) - 20682LL*Power(r,7LL)*Power(xj,7LL) + 

               1305LL*Power(r,8LL)*Power(xj,8LL) + 106LL*Power(r,9LL)*Power(xj,9LL)) + 

            3LL*Power(xi,8LL)*Power(xj,16LL)*

             (-9364244085LL + 6940428705LL*r*xj + 

               2117684520LL*Power(r,2LL)*Power(xj,2LL) - 

               230268150LL*Power(r,3LL)*Power(xj,3LL) - 

               149610510LL*Power(r,4LL)*Power(xj,4LL) - 

               21824334LL*Power(r,5LL)*Power(xj,5LL) - 

               1223208LL*Power(r,6LL)*Power(xj,6LL) + 12708LL*Power(r,7LL)*Power(xj,7LL) + 

               4470LL*Power(r,8LL)*Power(xj,8LL) + 146LL*Power(r,9LL)*Power(xj,9LL)) - 

            Power(xi,6LL)*Power(xj,18LL)*

             (57304872765LL + 7147185255LL*r*xj - 

               5801702760LL*Power(r,2LL)*Power(xj,2LL) - 

               2053388610LL*Power(r,3LL)*Power(xj,3LL) - 

               271655370LL*Power(r,4LL)*Power(xj,4LL) - 

               10864854LL*Power(r,5LL)*Power(xj,5LL) + 

               1337112LL*Power(r,6LL)*Power(xj,6LL) + 202716LL*Power(r,7LL)*Power(xj,7LL) + 

               10746LL*Power(r,8LL)*Power(xj,8LL) + 212LL*Power(r,9LL)*Power(xj,9LL))))/

       (56700LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,17LL)*Power(xi + xj,17LL))

    ; }
   
  }
  return S;
}

cl_R DSlater_4S_1S(cl_R r,cl_R xi,cl_R xj)
{
  return DSlater_1S_4S(r,xj,xi);
}

cl_R DSlater_4S_2S(cl_R r,cl_R xi,cl_R xj)
{
  return DSlater_2S_4S(r,xj,xi);
}

cl_R DSlater_4S_3S(cl_R r,cl_R xi,cl_R xj)
{
  return DSlater_3S_4S(r,xj,xi);
}

cl_R DSlater_5S_5S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = 0LL

    ; } else {  S = -(-(1e9*cl_float(223247882,precision)+3524011562LL)*xi + (1e9*cl_float(243290200,precision)+8176640000LL)*exp(2LL*r*xi)*xi - 

          (1e9*cl_float(406411127,precision)+7742766250LL)*r*Power(xi,2LL) - 

          (1e9*cl_float(366777722,precision)+5928565750LL)*Power(r,2LL)*Power(xi,3LL) - 

          (1e9*cl_float(218697853,precision)+9397652000LL)*Power(r,3LL)*Power(xi,4LL) - 

          (1e9*cl_float(968809971,precision)+432303000LL)*Power(r,4LL)*Power(xi,5LL) - 

          (1e9*cl_float(339917275,precision)+693195200LL)*Power(r,5LL)*Power(xi,6LL) - 

          983239817883523200LL*Power(r,6LL)*Power(xi,7LL) - 

          240924879420825600LL*Power(r,7LL)*Power(xi,8LL) - 

          50973581199340800LL*Power(r,8LL)*Power(xi,9LL) - 

          9439831425024000LL*Power(r,9LL)*Power(xi,10LL) - 

          1544699687731200LL*Power(r,10LL)*Power(xi,11LL) - 

          224683590942720LL*Power(r,11LL)*Power(xi,12LL) - 

          29125650677760LL*Power(r,12LL)*Power(xi,13LL) - 

          3360652001280LL*Power(r,13LL)*Power(xi,14LL) - 

          342923673600LL*Power(r,14LL)*Power(xi,15LL) - 

          30482104320LL*Power(r,15LL)*Power(xi,16LL) - 

          2286157824LL*Power(r,16LL)*Power(xi,17LL) - 

          134479872LL*Power(r,17LL)*Power(xi,18LL) - 4980736LL*Power(r,18LL)*Power(xi,19LL))/

       (1.21645100408832e19*exp(2LL*r*xi)*r) + 

      (-(1e9*cl_float(121645100,precision)+4088320000LL) + (1e9*cl_float(121645100,precision)+4088320000LL)*exp(2LL*r*xi) - 

         (1e9*cl_float(223247882,precision)+3524011562LL)*r*xi - 

         (1e9*cl_float(203205563,precision)+8871383125LL)*Power(r,2LL)*Power(xi,2LL) - 

         (1e9*cl_float(122259240,precision)+8642855250LL)*Power(r,3LL)*Power(xi,3LL) - 

         (1e9*cl_float(546744634,precision)+849413000LL)*Power(r,4LL)*Power(xi,4LL) - 

         (1e9*cl_float(193761994,precision)+286460600LL)*Power(r,5LL)*Power(xi,5LL) - 

         566528792821992000LL*Power(r,6LL)*Power(xi,6LL) - 

         140462831126217600LL*Power(r,7LL)*Power(xi,7LL) - 

         30115609927603200LL*Power(r,8LL)*Power(xi,8LL) - 

         5663731244371200LL*Power(r,9LL)*Power(xi,9LL) - 

         943983142502400LL*Power(r,10LL)*Power(xi,10LL) - 

         140427244339200LL*Power(r,11LL)*Power(xi,11LL) - 

         18723632578560LL*Power(r,12LL)*Power(xi,12LL) - 

         2240434667520LL*Power(r,13LL)*Power(xi,13LL) - 

         240046571520LL*Power(r,14LL)*Power(xi,14LL) - 

         22861578240LL*Power(r,15LL)*Power(xi,15LL) - 

         1905131520LL*Power(r,16LL)*Power(xi,16LL) - 

         134479872LL*Power(r,17LL)*Power(xi,17LL) - 

         7471104LL*Power(r,18LL)*Power(xi,18LL) - 262144LL*Power(r,19LL)*Power(xi,19LL))/

       (1.21645100408832e19*exp(2LL*r*xi)*Power(r,2LL)) + 

      (xi*(-(1e9*cl_float(121645100,precision)+4088320000LL) + (1e9*cl_float(121645100,precision)+4088320000LL)*exp(2LL*r*xi) - 

           (1e9*cl_float(223247882,precision)+3524011562LL)*r*xi - 

           (1e9*cl_float(203205563,precision)+8871383125LL)*Power(r,2LL)*Power(xi,2LL) - 

           (1e9*cl_float(122259240,precision)+8642855250LL)*Power(r,3LL)*Power(xi,3LL) - 

           (1e9*cl_float(546744634,precision)+849413000LL)*Power(r,4LL)*Power(xi,4LL) - 

           (1e9*cl_float(193761994,precision)+286460600LL)*Power(r,5LL)*Power(xi,5LL) - 

           566528792821992000LL*Power(r,6LL)*Power(xi,6LL) - 

           140462831126217600LL*Power(r,7LL)*Power(xi,7LL) - 

           30115609927603200LL*Power(r,8LL)*Power(xi,8LL) - 

           5663731244371200LL*Power(r,9LL)*Power(xi,9LL) - 

           943983142502400LL*Power(r,10LL)*Power(xi,10LL) - 

           140427244339200LL*Power(r,11LL)*Power(xi,11LL) - 

           18723632578560LL*Power(r,12LL)*Power(xi,12LL) - 

           2240434667520LL*Power(r,13LL)*Power(xi,13LL) - 

           240046571520LL*Power(r,14LL)*Power(xi,14LL) - 

           22861578240LL*Power(r,15LL)*Power(xi,15LL) - 

           1905131520LL*Power(r,16LL)*Power(xi,16LL) - 

           134479872LL*Power(r,17LL)*Power(xi,17LL) - 

           7471104LL*Power(r,18LL)*Power(xi,18LL) - 262144LL*Power(r,19LL)*Power(xi,19LL)))/

       (6.0822550204416e18*exp(2LL*r*xi)*r)

    ; }
 
  }
  else {
      if (r == 0LL) {  S = 0LL

    ; } else {  S = (70875LL*exp(2LL*r*(xi + xj))*Power(Power(xi,2LL) - Power(xj,2LL),19LL) + 

         exp(2LL*r*xj)*Power(xj,12LL)*

          (-630LL*Power(r,8LL)*Power(xi,34LL) - 10LL*Power(r,9LL)*Power(xi,35LL) + 

            70875LL*Power(xj,26LL) + 127575LL*r*xi*Power(xj,26LL) - 

            30LL*Power(r,7LL)*Power(xi,33LL)*(630LL + Power(r,2LL)*Power(xj,2LL)) + 

            14175LL*Power(xi,2LL)*Power(xj,24LL)*(-95LL + 8LL*Power(r,2LL)*Power(xj,2LL)) + 

            4725LL*r*Power(xi,3LL)*Power(xj,24LL)*

             (-513LL + 14LL*Power(r,2LL)*Power(xj,2LL)) - 

            90LL*Power(r,6LL)*Power(xi,32LL)*(3920LL + 43LL*Power(r,2LL)*Power(xj,2LL)) + 

            4725LL*r*Power(xi,5LL)*Power(xj,22LL)*

             (4617LL - 266LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            14175LL*Power(xi,4LL)*Power(xj,22LL)*

             (855LL - 152LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            36LL*Power(r,5LL)*Power(xi,31LL)*

             (-124950LL - 4985LL*Power(r,2LL)*Power(xj,2LL) + 

               13LL*Power(r,4LL)*Power(xj,4LL)) + 

            36LL*Power(r,4LL)*Power(xi,30LL)*

             (-1124550LL - 127960LL*Power(r,2LL)*Power(xj,2LL) + 

               863LL*Power(r,4LL)*Power(xj,4LL)) + 

            135LL*r*Power(xi,7LL)*Power(xj,20LL)*

             (-915705LL + 83790LL*Power(r,2LL)*Power(xj,2LL) - 

               1330LL*Power(r,4LL)*Power(xj,4LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) + 

            315LL*Power(xi,6LL)*Power(xj,20LL)*

             (-218025LL + 61560LL*Power(r,2LL)*Power(xj,2LL) - 

               1710LL*Power(r,4LL)*Power(xj,4LL) + 8LL*Power(r,6LL)*Power(xj,6LL)) - 

            36LL*Power(r,3LL)*Power(xi,29LL)*

             (7122150LL + 2102730LL*Power(r,2LL)*Power(xj,2LL) - 

               23294LL*Power(r,4LL)*Power(xj,4LL) + 37LL*Power(r,6LL)*Power(xj,6LL)) - 

            36LL*Power(r,2LL)*Power(xi,28LL)*

             (30523500LL + 23401350LL*Power(r,2LL)*Power(xj,2LL) - 

               299250LL*Power(r,4LL)*Power(xj,4LL) + 1297LL*Power(r,6LL)*Power(xj,6LL)) + 

            r*Power(xi,17LL)*Power(xj,10LL)*

             (1073961177975LL - 21753487980LL*Power(r,2LL)*Power(xj,2LL) - 

               745994340LL*Power(r,4LL)*Power(xj,4LL) + 

               5307156LL*Power(r,6LL)*Power(xj,6LL) - 818LL*Power(r,8LL)*Power(xj,8LL)) + 

            10LL*r*Power(xi,9LL)*Power(xj,18LL)*

             (49448070LL - 6409935LL*Power(r,2LL)*Power(xj,2LL) + 

               161595LL*Power(r,4LL)*Power(xj,4LL) - 1026LL*Power(r,6LL)*Power(xj,6LL) + 

               Power(r,8LL)*Power(xj,8LL)) + 

            90LL*Power(xi,8LL)*Power(xj,18LL)*

             (3052350LL - 1220940LL*Power(r,2LL)*Power(xj,2LL) + 

               53865LL*Power(r,4LL)*Power(xj,4LL) - 532LL*Power(r,6LL)*Power(xj,6LL) + 

               Power(r,8LL)*Power(xj,8LL)) - 

            1710LL*Power(xi,10LL)*Power(xj,16LL)*

             (481950LL - 257040LL*Power(r,2LL)*Power(xj,2LL) + 

               16065LL*Power(r,4LL)*Power(xj,4LL) - 252LL*Power(r,6LL)*Power(xj,6LL) + 

               Power(r,8LL)*Power(xj,8LL)) + 

            6LL*r*Power(xi,11LL)*Power(xj,16LL)*

             (-207559800LL + 50390550LL*Power(r,2LL)*Power(xj,2LL) - 

               1165815LL*Power(r,4LL)*Power(xj,4LL) + 

               21396LL*Power(r,6LL)*Power(xj,6LL) + 5LL*Power(r,8LL)*Power(xj,8LL)) - 

            18LL*r*Power(xi,13LL)*Power(xj,14LL)*

             (-1703720025LL - 155669850LL*Power(r,2LL)*Power(xj,2LL) - 

               7410270LL*Power(r,4LL)*Power(xj,4LL) - 1532LL*Power(r,6LL)*Power(xj,6LL) + 

               26LL*Power(r,8LL)*Power(xj,8LL)) + 

            18LL*r*Power(xi,15LL)*Power(xj,12LL)*

             (19380896325LL + 1329128850LL*Power(r,2LL)*Power(xj,2LL) - 

               7608930LL*Power(r,4LL)*Power(xj,4LL) - 

               116238LL*Power(r,6LL)*Power(xj,6LL) + 74LL*Power(r,8LL)*Power(xj,8LL)) - 

            18LL*Power(xi,12LL)*Power(xj,14LL)*

             (89026875LL + 179071200LL*Power(r,2LL)*Power(xj,2LL) + 

               1552950LL*Power(r,4LL)*Power(xj,4LL) + 

               295820LL*Power(r,6LL)*Power(xj,6LL) + 146LL*Power(r,8LL)*Power(xj,8LL)) + 

            18LL*r*Power(xi,25LL)*Power(xj,2LL)*

             (-5449970925LL - 1137574935LL*Power(r,2LL)*Power(xj,2LL) + 

               37834755LL*Power(r,4LL)*Power(xj,4LL) - 

               273062LL*Power(r,6LL)*Power(xj,6LL) + 171LL*Power(r,8LL)*Power(xj,8LL)) - 

            9LL*r*Power(xi,19LL)*Power(xj,8LL)*

             (-37914907275LL + 7613889570LL*Power(r,2LL)*Power(xj,2LL) - 

               170524620LL*Power(r,4LL)*Power(xj,4LL) + 

               397936LL*Power(r,6LL)*Power(xj,6LL) + 342LL*Power(r,8LL)*Power(xj,8LL)) - 

            3LL*r*Power(xi,23LL)*Power(xj,4LL)*

             (219130630425LL - 11118046590LL*Power(r,2LL)*Power(xj,2LL) + 

               327611970LL*Power(r,4LL)*Power(xj,4LL) - 

               2920908LL*Power(r,6LL)*Power(xj,6LL) + 2584LL*Power(r,8LL)*Power(xj,8LL)) + 

            3LL*r*Power(xi,21LL)*Power(xj,6LL)*

             (-345162539925LL + 19030764690LL*Power(r,2LL)*Power(xj,2LL) - 

               141976170LL*Power(r,4LL)*Power(xj,4LL) - 

               1441872LL*Power(r,6LL)*Power(xj,6LL) + 2584LL*Power(r,8LL)*Power(xj,8LL)) + 

            63LL*Power(xi,20LL)*Power(xj,6LL)*

             (-50980542525LL + 6240202920LL*Power(r,2LL)*Power(xj,2LL) - 

               201314310LL*Power(r,4LL)*Power(xj,4LL) + 

               956080LL*Power(r,6LL)*Power(xj,6LL) + 2584LL*Power(r,8LL)*Power(xj,8LL)) + 

            18LL*Power(xi,14LL)*Power(xj,12LL)*

             (-7803332775LL - 2519206200LL*Power(r,2LL)*Power(xj,2LL) - 

               119719950LL*Power(r,4LL)*Power(xj,4LL) + 

               182280LL*Power(r,6LL)*Power(xj,6LL) + 2734LL*Power(r,8LL)*Power(xj,8LL)) - 

            18LL*Power(xi,26LL)*(195859125LL + 1794781800LL*Power(r,2LL)*Power(xj,2LL) + 

               67337235LL*Power(r,4LL)*Power(xj,4LL) - 

               1659700LL*Power(r,6LL)*Power(xj,6LL) + 4089LL*Power(r,8LL)*Power(xj,8LL)) + 

            9LL*Power(xi,18LL)*Power(xj,8LL)*

             (-357591274425LL + 8328390840LL*Power(r,2LL)*Power(xj,2LL) + 

               912042180LL*Power(r,4LL)*Power(xj,4LL) - 

               12842480LL*Power(r,6LL)*Power(xj,6LL) + 10678LL*Power(r,8LL)*Power(xj,8LL)) 

    - 9LL*Power(xi,16LL)*Power(xj,10LL)*

             (128599724925LL + 21298077360LL*Power(r,2LL)*Power(xj,2LL) - 

               267928500LL*Power(r,4LL)*Power(xj,4LL) - 

               5458320LL*Power(r,6LL)*Power(xj,6LL) + 14722LL*Power(r,8LL)*Power(xj,8LL)) 

    + 18LL*Power(xi,24LL)*Power(xj,2LL)*

             (-7604930025LL - 8866107180LL*Power(r,2LL)*Power(xj,2LL) + 

               399272265LL*Power(r,4LL)*Power(xj,4LL) - 

               5925780LL*Power(r,6LL)*Power(xj,6LL) + 17651LL*Power(r,8LL)*Power(xj,8LL)) 

    - 9LL*Power(xi,22LL)*Power(xj,4LL)*(129194933175LL + 

               3909863160LL*Power(r,2LL)*Power(xj,2LL) + 

               91420770LL*Power(r,4LL)*Power(xj,4LL) - 

               8762040LL*Power(r,6LL)*Power(xj,6LL) + 43928LL*Power(r,8LL)*Power(xj,8LL)) 

    + Power(xi,27LL)*(-2884470750LL*r - 6409935000LL*Power(r,3LL)*Power(xj,2LL) + 

               28332990LL*Power(r,5LL)*Power(xj,4LL) + 

               58104LL*Power(r,7LL)*Power(xj,6LL) + 818LL*Power(r,9LL)*Power(xj,8LL))) + 

         exp(2LL*r*xi)*Power(xi,12LL)*

          (Power(xi,8LL)*Power(xj,18LL)*

             (3218321469825LL - 341234165475LL*r*xj - 

               393132783960LL*Power(r,2LL)*Power(xj,2LL) - 

               57092294070LL*Power(r,3LL)*Power(xj,3LL) + 

               822786930LL*Power(r,4LL)*Power(xj,4LL) + 

               982835910LL*Power(r,5LL)*Power(xj,5LL) + 

               106664040LL*Power(r,6LL)*Power(xj,6LL) + 

               4915116LL*Power(r,7LL)*Power(xj,7LL) + 73602LL*Power(r,8LL)*Power(xj,8LL) - 

               818LL*Power(r,9LL)*Power(xj,9LL)) + 

            10LL*Power(xj,26LL)*(352546425LL + 288447075LL*r*xj + 

               109884600LL*Power(r,2LL)*Power(xj,2LL) + 

               25639740LL*Power(r,3LL)*Power(xj,3LL) + 

               4048380LL*Power(r,4LL)*Power(xj,4LL) + 

               449820LL*Power(r,5LL)*Power(xj,5LL) + 35280LL*Power(r,6LL)*Power(xj,6LL) + 

               1890LL*Power(r,7LL)*Power(xj,7LL) + 63LL*Power(r,8LL)*Power(xj,8LL) + 

               Power(r,9LL)*Power(xj,9LL)) + 

            30LL*Power(xi,2LL)*Power(xj,24LL)*

             (4562958015LL + 3269982555LL*r*xj + 

               1076869080LL*Power(r,2LL)*Power(xj,2LL) + 

               213664500LL*Power(r,3LL)*Power(xj,3LL) + 

               28081620LL*Power(r,4LL)*Power(xj,4LL) + 

               2523276LL*Power(r,5LL)*Power(xj,5LL) + 

               153552LL*Power(r,6LL)*Power(xj,6LL) + 5982LL*Power(r,7LL)*Power(xj,7LL) + 

               129LL*Power(r,8LL)*Power(xj,8LL) + Power(r,9LL)*Power(xj,9LL)) - 

            15LL*Power(xi,24LL)*Power(xj,2LL)*

             (-89775LL - 161595LL*r*xj - 143640LL*Power(r,2LL)*Power(xj,2LL) - 

               83790LL*Power(r,3LL)*Power(xj,3LL) - 35910LL*Power(r,4LL)*Power(xj,4LL) - 

               11970LL*Power(r,5LL)*Power(xj,5LL) - 3192LL*Power(r,6LL)*Power(xj,6LL) - 

               684LL*Power(r,7LL)*Power(xj,7LL) - 114LL*Power(r,8LL)*Power(xj,8LL) + 

               2LL*Power(r,9LL)*Power(xj,9LL)) - 

            5LL*Power(xi,26LL)*(14175LL + 25515LL*r*xj + 

               22680LL*Power(r,2LL)*Power(xj,2LL) + 13230LL*Power(r,3LL)*Power(xj,3LL) + 

               5670LL*Power(r,4LL)*Power(xj,4LL) + 1890LL*Power(r,5LL)*Power(xj,5LL) + 

               504LL*Power(r,6LL)*Power(xj,6LL) + 108LL*Power(r,7LL)*Power(xj,7LL) + 

               18LL*Power(r,8LL)*Power(xj,8LL) + 2LL*Power(r,9LL)*Power(xj,9LL)) - 

            1938LL*Power(xi,14LL)*Power(xj,12LL)*

             (-826875LL + 15824025LL*r*xj - 23398200LL*Power(r,2LL)*Power(xj,2LL) + 

               12344850LL*Power(r,3LL)*Power(xj,3LL) + 

               1244250LL*Power(r,4LL)*Power(xj,4LL) - 

               384930LL*Power(r,5LL)*Power(xj,5LL) - 59640LL*Power(r,6LL)*Power(xj,6LL) - 

               1848LL*Power(r,7LL)*Power(xj,7LL) + 84LL*Power(r,8LL)*Power(xj,8LL) + 

               4LL*Power(r,9LL)*Power(xj,9LL)) + 

            1938LL*Power(xi,12LL)*Power(xj,14LL)*

             (72476775LL - 180008325LL*r*xj + 98907480LL*Power(r,2LL)*Power(xj,2LL) + 

               11224710LL*Power(r,3LL)*Power(xj,3LL) - 

               4235490LL*Power(r,4LL)*Power(xj,4LL) - 

               791910LL*Power(r,5LL)*Power(xj,5LL) - 31080LL*Power(r,6LL)*Power(xj,6LL) + 

               2232LL*Power(r,7LL)*Power(xj,7LL) + 204LL*Power(r,8LL)*Power(xj,8LL) + 

               4LL*Power(r,9LL)*Power(xj,9LL)) + 

            342LL*Power(xi,16LL)*Power(xj,10LL)*

             (2409750LL + 3641400LL*r*xj + 9424800LL*Power(r,2LL)*Power(xj,2LL) - 

               8193150LL*Power(r,3LL)*Power(xj,3LL) + 

               6301050LL*Power(r,4LL)*Power(xj,4LL) + 

               400470LL*Power(r,5LL)*Power(xj,5LL) - 143640LL*Power(r,6LL)*Power(xj,6LL) - 

               15518LL*Power(r,7LL)*Power(xj,7LL) - 281LL*Power(r,8LL)*Power(xj,8LL) + 

               9LL*Power(r,9LL)*Power(xj,9LL)) - 

            171LL*Power(xi,10LL)*Power(xj,16LL)*

             (-6768406575LL + 6280474725LL*r*xj + 

               438336360LL*Power(r,2LL)*Power(xj,2LL) - 

               400731030LL*Power(r,3LL)*Power(xj,3LL) - 

               74168430LL*Power(r,4LL)*Power(xj,4LL) - 

               2490810LL*Power(r,5LL)*Power(xj,5LL) + 

               461160LL*Power(r,6LL)*Power(xj,6LL) + 51244LL*Power(r,7LL)*Power(xj,7LL) + 

               1858LL*Power(r,8LL)*Power(xj,8LL) + 18LL*Power(r,9LL)*Power(xj,9LL)) + 

            9LL*Power(xi,22LL)*Power(xj,4LL)*

             (-1346625LL - 2423925LL*r*xj - 2154600LL*Power(r,2LL)*Power(xj,2LL) - 

               1256850LL*Power(r,3LL)*Power(xj,3LL) - 

               538650LL*Power(r,4LL)*Power(xj,4LL) - 179550LL*Power(r,5LL)*Power(xj,5LL) - 

               47880LL*Power(r,6LL)*Power(xj,6LL) - 14264LL*Power(r,7LL)*Power(xj,7LL) + 

               292LL*Power(r,8LL)*Power(xj,8LL) + 52LL*Power(r,9LL)*Power(xj,9LL)) - 

            9LL*Power(xi,4LL)*Power(xj,22LL)*

             (-129194933175LL - 73043543475LL*r*xj - 

               17732214360LL*Power(r,2LL)*Power(xj,2LL) - 

               2275149870LL*Power(r,3LL)*Power(xj,3LL) - 

               134674470LL*Power(r,4LL)*Power(xj,4LL) + 

               3148110LL*Power(r,5LL)*Power(xj,5LL) + 

               1197000LL*Power(r,6LL)*Power(xj,6LL) + 93176LL*Power(r,7LL)*Power(xj,7LL) + 

               3452LL*Power(r,8LL)*Power(xj,8LL) + 52LL*Power(r,9LL)*Power(xj,9LL)) + 

            9LL*Power(xi,6LL)*Power(xj,20LL)*

             (356863797675LL + 115054179975LL*r*xj + 

               3909863160LL*Power(r,2LL)*Power(xj,2LL) - 

               3706015530LL*Power(r,3LL)*Power(xj,3LL) - 

               798544530LL*Power(r,4LL)*Power(xj,4LL) - 

               75669510LL*Power(r,5LL)*Power(xj,5LL) - 

               3319400LL*Power(r,6LL)*Power(xj,6LL) - 6456LL*Power(r,7LL)*Power(xj,7LL) + 

               5188LL*Power(r,8LL)*Power(xj,8LL) + 148LL*Power(r,9LL)*Power(xj,9LL)) - 

            9LL*Power(xi,20LL)*Power(xj,6LL)*

             (-7630875LL - 13735575LL*r*xj - 12209400LL*Power(r,2LL)*Power(xj,2LL) - 

               7122150LL*Power(r,3LL)*Power(xj,3LL) - 

               3052350LL*Power(r,4LL)*Power(xj,4LL) - 

               777210LL*Power(r,5LL)*Power(xj,5LL) - 591640LL*Power(r,6LL)*Power(xj,6LL) + 

               3064LL*Power(r,7LL)*Power(xj,7LL) + 5468LL*Power(r,8LL)*Power(xj,8LL) + 

               148LL*Power(r,9LL)*Power(xj,9LL)) + 

            2LL*Power(xi,18LL)*Power(xj,8LL)*

             (-137355750LL - 247240350LL*r*xj - 219769200LL*Power(r,2LL)*Power(xj,2LL) - 

               151171650LL*Power(r,3LL)*Power(xj,3LL) + 

               13976550LL*Power(r,4LL)*Power(xj,4LL) - 

               66692430LL*Power(r,5LL)*Power(xj,5LL) - 

               1640520LL*Power(r,6LL)*Power(xj,6LL) + 

               1046142LL*Power(r,7LL)*Power(xj,7LL) + 66249LL*Power(r,8LL)*Power(xj,8LL) + 

               409LL*Power(r,9LL)*Power(xj,9LL))))/

       (70875LL*exp(2LL*r*(xi + xj))*Power(r,2LL)*Power(xi - xj,19LL)*

         Power(xi + xj,19LL)) + (2LL*(70875LL*exp(2LL*r*(xi + xj))*

            Power(Power(xi,2LL) - Power(xj,2LL),19LL) + 

           exp(2LL*r*xj)*Power(xj,12LL)*

            (-630LL*Power(r,8LL)*Power(xi,34LL) - 10LL*Power(r,9LL)*Power(xi,35LL) + 

              70875LL*Power(xj,26LL) + 127575LL*r*xi*Power(xj,26LL) - 

              30LL*Power(r,7LL)*Power(xi,33LL)*(630LL + Power(r,2LL)*Power(xj,2LL)) + 

              14175LL*Power(xi,2LL)*Power(xj,24LL)*

               (-95LL + 8LL*Power(r,2LL)*Power(xj,2LL)) + 

              4725LL*r*Power(xi,3LL)*Power(xj,24LL)*

               (-513LL + 14LL*Power(r,2LL)*Power(xj,2LL)) - 

              90LL*Power(r,6LL)*Power(xi,32LL)*(3920LL + 43LL*Power(r,2LL)*Power(xj,2LL)) + 

              4725LL*r*Power(xi,5LL)*Power(xj,22LL)*

               (4617LL - 266LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) 

    + 14175LL*Power(xi,4LL)*Power(xj,22LL)*

               (855LL - 152LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

              36LL*Power(r,5LL)*Power(xi,31LL)*

               (-124950LL - 4985LL*Power(r,2LL)*Power(xj,2LL) + 

                 13LL*Power(r,4LL)*Power(xj,4LL)) + 

              36LL*Power(r,4LL)*Power(xi,30LL)*

               (-1124550LL - 127960LL*Power(r,2LL)*Power(xj,2LL) + 

                 863LL*Power(r,4LL)*Power(xj,4LL)) + 

              135LL*r*Power(xi,7LL)*Power(xj,20LL)*

               (-915705LL + 83790LL*Power(r,2LL)*Power(xj,2LL) - 

                 1330LL*Power(r,4LL)*Power(xj,4LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) + 

              315LL*Power(xi,6LL)*Power(xj,20LL)*

               (-218025LL + 61560LL*Power(r,2LL)*Power(xj,2LL) - 

                 1710LL*Power(r,4LL)*Power(xj,4LL) + 8LL*Power(r,6LL)*Power(xj,6LL)) - 

              36LL*Power(r,3LL)*Power(xi,29LL)*

               (7122150LL + 2102730LL*Power(r,2LL)*Power(xj,2LL) - 

                 23294LL*Power(r,4LL)*Power(xj,4LL) + 37LL*Power(r,6LL)*Power(xj,6LL)) - 

              36LL*Power(r,2LL)*Power(xi,28LL)*

               (30523500LL + 23401350LL*Power(r,2LL)*Power(xj,2LL) - 

                 299250LL*Power(r,4LL)*Power(xj,4LL) + 1297LL*Power(r,6LL)*Power(xj,6LL)) 

    + r*Power(xi,17LL)*Power(xj,10LL)*

               (1073961177975LL - 21753487980LL*Power(r,2LL)*Power(xj,2LL) - 

                 745994340LL*Power(r,4LL)*Power(xj,4LL) + 

                 5307156LL*Power(r,6LL)*Power(xj,6LL) - 818LL*Power(r,8LL)*Power(xj,8LL)) 

    + 10LL*r*Power(xi,9LL)*Power(xj,18LL)*

               (49448070LL - 6409935LL*Power(r,2LL)*Power(xj,2LL) + 

                 161595LL*Power(r,4LL)*Power(xj,4LL) - 

                 1026LL*Power(r,6LL)*Power(xj,6LL) + Power(r,8LL)*Power(xj,8LL)) + 

              90LL*Power(xi,8LL)*Power(xj,18LL)*

               (3052350LL - 1220940LL*Power(r,2LL)*Power(xj,2LL) + 

                 53865LL*Power(r,4LL)*Power(xj,4LL) - 532LL*Power(r,6LL)*Power(xj,6LL) + 

                 Power(r,8LL)*Power(xj,8LL)) - 

              1710LL*Power(xi,10LL)*Power(xj,16LL)*

               (481950LL - 257040LL*Power(r,2LL)*Power(xj,2LL) + 

                 16065LL*Power(r,4LL)*Power(xj,4LL) - 252LL*Power(r,6LL)*Power(xj,6LL) + 

                 Power(r,8LL)*Power(xj,8LL)) + 

              6LL*r*Power(xi,11LL)*Power(xj,16LL)*

               (-207559800LL + 50390550LL*Power(r,2LL)*Power(xj,2LL) - 

                 1165815LL*Power(r,4LL)*Power(xj,4LL) + 

                 21396LL*Power(r,6LL)*Power(xj,6LL) + 5LL*Power(r,8LL)*Power(xj,8LL)) - 

              18LL*r*Power(xi,13LL)*Power(xj,14LL)*

               (-1703720025LL - 155669850LL*Power(r,2LL)*Power(xj,2LL) - 

                 7410270LL*Power(r,4LL)*Power(xj,4LL) - 

                 1532LL*Power(r,6LL)*Power(xj,6LL) + 26LL*Power(r,8LL)*Power(xj,8LL)) + 

              18LL*r*Power(xi,15LL)*Power(xj,12LL)*

               (19380896325LL + 1329128850LL*Power(r,2LL)*Power(xj,2LL) - 

                 7608930LL*Power(r,4LL)*Power(xj,4LL) - 

                 116238LL*Power(r,6LL)*Power(xj,6LL) + 74LL*Power(r,8LL)*Power(xj,8LL)) - 

              18LL*Power(xi,12LL)*Power(xj,14LL)*

               (89026875LL + 179071200LL*Power(r,2LL)*Power(xj,2LL) + 

                 1552950LL*Power(r,4LL)*Power(xj,4LL) + 

                 295820LL*Power(r,6LL)*Power(xj,6LL) + 146LL*Power(r,8LL)*Power(xj,8LL)) + 

              18LL*r*Power(xi,25LL)*Power(xj,2LL)*

               (-5449970925LL - 1137574935LL*Power(r,2LL)*Power(xj,2LL) + 

                 37834755LL*Power(r,4LL)*Power(xj,4LL) - 

                 273062LL*Power(r,6LL)*Power(xj,6LL) + 171LL*Power(r,8LL)*Power(xj,8LL)) - 

              9LL*r*Power(xi,19LL)*Power(xj,8LL)*

               (-37914907275LL + 7613889570LL*Power(r,2LL)*Power(xj,2LL) - 

                 170524620LL*Power(r,4LL)*Power(xj,4LL) + 

                 397936LL*Power(r,6LL)*Power(xj,6LL) + 342LL*Power(r,8LL)*Power(xj,8LL)) - 

              3LL*r*Power(xi,23LL)*Power(xj,4LL)*

               (219130630425LL - 11118046590LL*Power(r,2LL)*Power(xj,2LL) + 

                 327611970LL*Power(r,4LL)*Power(xj,4LL) - 

                 2920908LL*Power(r,6LL)*Power(xj,6LL) + 2584LL*Power(r,8LL)*Power(xj,8LL)) 

    + 3LL*r*Power(xi,21LL)*Power(xj,6LL)*

               (-345162539925LL + 19030764690LL*Power(r,2LL)*Power(xj,2LL) - 

                 141976170LL*Power(r,4LL)*Power(xj,4LL) - 

                 1441872LL*Power(r,6LL)*Power(xj,6LL) + 2584LL*Power(r,8LL)*Power(xj,8LL)) 

    + 63LL*Power(xi,20LL)*Power(xj,6LL)*

               (-50980542525LL + 6240202920LL*Power(r,2LL)*Power(xj,2LL) - 

                 201314310LL*Power(r,4LL)*Power(xj,4LL) + 

                 956080LL*Power(r,6LL)*Power(xj,6LL) + 2584LL*Power(r,8LL)*Power(xj,8LL)) 

    + 18LL*Power(xi,14LL)*Power(xj,12LL)*

               (-7803332775LL - 2519206200LL*Power(r,2LL)*Power(xj,2LL) - 

                 119719950LL*Power(r,4LL)*Power(xj,4LL) + 

                 182280LL*Power(r,6LL)*Power(xj,6LL) + 2734LL*Power(r,8LL)*Power(xj,8LL)) 

    - 18LL*Power(xi,26LL)*(195859125LL + 1794781800LL*Power(r,2LL)*Power(xj,2LL) + 

                 67337235LL*Power(r,4LL)*Power(xj,4LL) - 

                 1659700LL*Power(r,6LL)*Power(xj,6LL) + 4089LL*Power(r,8LL)*Power(xj,8LL)) 

    + 9LL*Power(xi,18LL)*Power(xj,8LL)*(-357591274425LL + 

                 8328390840LL*Power(r,2LL)*Power(xj,2LL) + 

                 912042180LL*Power(r,4LL)*Power(xj,4LL) - 

                 12842480LL*Power(r,6LL)*Power(xj,6LL) + 

                 10678LL*Power(r,8LL)*Power(xj,8LL)) - 

              9LL*Power(xi,16LL)*Power(xj,10LL)*

               (128599724925LL + 21298077360LL*Power(r,2LL)*Power(xj,2LL) - 

                 267928500LL*Power(r,4LL)*Power(xj,4LL) - 

                 5458320LL*Power(r,6LL)*Power(xj,6LL) + 14722LL*Power(r,8LL)*Power(xj,8LL)

    ) + 18LL*Power(xi,24LL)*Power(xj,2LL)*

               (-7604930025LL - 8866107180LL*Power(r,2LL)*Power(xj,2LL) + 

                 399272265LL*Power(r,4LL)*Power(xj,4LL) - 

                 5925780LL*Power(r,6LL)*Power(xj,6LL) + 17651LL*Power(r,8LL)*Power(xj,8LL)

    ) - 9LL*Power(xi,22LL)*Power(xj,4LL)*

               (129194933175LL + 3909863160LL*Power(r,2LL)*Power(xj,2LL) + 

                 91420770LL*Power(r,4LL)*Power(xj,4LL) - 

                 8762040LL*Power(r,6LL)*Power(xj,6LL) + 43928LL*Power(r,8LL)*Power(xj,8LL)

    ) + Power(xi,27LL)*(-2884470750LL*r - 6409935000LL*Power(r,3LL)*Power(xj,2LL) + 

                 28332990LL*Power(r,5LL)*Power(xj,4LL) + 

                 58104LL*Power(r,7LL)*Power(xj,6LL) + 818LL*Power(r,9LL)*Power(xj,8LL))) + 

           exp(2LL*r*xi)*Power(xi,12LL)*

            (Power(xi,8LL)*Power(xj,18LL)*

               (3218321469825LL - 341234165475LL*r*xj - 

                 393132783960LL*Power(r,2LL)*Power(xj,2LL) - 

                 57092294070LL*Power(r,3LL)*Power(xj,3LL) + 

                 822786930LL*Power(r,4LL)*Power(xj,4LL) + 

                 982835910LL*Power(r,5LL)*Power(xj,5LL) + 

                 106664040LL*Power(r,6LL)*Power(xj,6LL) + 

                 4915116LL*Power(r,7LL)*Power(xj,7LL) + 

                 73602LL*Power(r,8LL)*Power(xj,8LL) - 818LL*Power(r,9LL)*Power(xj,9LL)) + 

              10LL*Power(xj,26LL)*(352546425LL + 288447075LL*r*xj + 

                 109884600LL*Power(r,2LL)*Power(xj,2LL) + 

                 25639740LL*Power(r,3LL)*Power(xj,3LL) + 

                 4048380LL*Power(r,4LL)*Power(xj,4LL) + 

                 449820LL*Power(r,5LL)*Power(xj,5LL) + 

                 35280LL*Power(r,6LL)*Power(xj,6LL) + 1890LL*Power(r,7LL)*Power(xj,7LL) + 

                 63LL*Power(r,8LL)*Power(xj,8LL) + Power(r,9LL)*Power(xj,9LL)) + 

              30LL*Power(xi,2LL)*Power(xj,24LL)*

               (4562958015LL + 3269982555LL*r*xj + 

                 1076869080LL*Power(r,2LL)*Power(xj,2LL) + 

                 213664500LL*Power(r,3LL)*Power(xj,3LL) + 

                 28081620LL*Power(r,4LL)*Power(xj,4LL) + 

                 2523276LL*Power(r,5LL)*Power(xj,5LL) + 

                 153552LL*Power(r,6LL)*Power(xj,6LL) + 5982LL*Power(r,7LL)*Power(xj,7LL) + 

                 129LL*Power(r,8LL)*Power(xj,8LL) + Power(r,9LL)*Power(xj,9LL)) - 

              15LL*Power(xi,24LL)*Power(xj,2LL)*

               (-89775LL - 161595LL*r*xj - 143640LL*Power(r,2LL)*Power(xj,2LL) - 

                 83790LL*Power(r,3LL)*Power(xj,3LL) - 35910LL*Power(r,4LL)*Power(xj,4LL) - 

                 11970LL*Power(r,5LL)*Power(xj,5LL) - 3192LL*Power(r,6LL)*Power(xj,6LL) - 

                 684LL*Power(r,7LL)*Power(xj,7LL) - 114LL*Power(r,8LL)*Power(xj,8LL) + 

                 2LL*Power(r,9LL)*Power(xj,9LL)) - 

              5LL*Power(xi,26LL)*(14175LL + 25515LL*r*xj + 

                 22680LL*Power(r,2LL)*Power(xj,2LL) + 13230LL*Power(r,3LL)*Power(xj,3LL) + 

                 5670LL*Power(r,4LL)*Power(xj,4LL) + 1890LL*Power(r,5LL)*Power(xj,5LL) + 

                 504LL*Power(r,6LL)*Power(xj,6LL) + 108LL*Power(r,7LL)*Power(xj,7LL) + 

                 18LL*Power(r,8LL)*Power(xj,8LL) + 2LL*Power(r,9LL)*Power(xj,9LL)) - 

              1938LL*Power(xi,14LL)*Power(xj,12LL)*

               (-826875LL + 15824025LL*r*xj - 23398200LL*Power(r,2LL)*Power(xj,2LL) + 

                 12344850LL*Power(r,3LL)*Power(xj,3LL) + 

                 1244250LL*Power(r,4LL)*Power(xj,4LL) - 

                 384930LL*Power(r,5LL)*Power(xj,5LL) - 

                 59640LL*Power(r,6LL)*Power(xj,6LL) - 1848LL*Power(r,7LL)*Power(xj,7LL) + 

                 84LL*Power(r,8LL)*Power(xj,8LL) + 4LL*Power(r,9LL)*Power(xj,9LL)) + 

              1938LL*Power(xi,12LL)*Power(xj,14LL)*

               (72476775LL - 180008325LL*r*xj + 98907480LL*Power(r,2LL)*Power(xj,2LL) + 

                 11224710LL*Power(r,3LL)*Power(xj,3LL) - 

                 4235490LL*Power(r,4LL)*Power(xj,4LL) - 

                 791910LL*Power(r,5LL)*Power(xj,5LL) - 

                 31080LL*Power(r,6LL)*Power(xj,6LL) + 2232LL*Power(r,7LL)*Power(xj,7LL) + 

                 204LL*Power(r,8LL)*Power(xj,8LL) + 4LL*Power(r,9LL)*Power(xj,9LL)) + 

              342LL*Power(xi,16LL)*Power(xj,10LL)*

               (2409750LL + 3641400LL*r*xj + 9424800LL*Power(r,2LL)*Power(xj,2LL) - 

                 8193150LL*Power(r,3LL)*Power(xj,3LL) + 

                 6301050LL*Power(r,4LL)*Power(xj,4LL) + 

                 400470LL*Power(r,5LL)*Power(xj,5LL) - 

                 143640LL*Power(r,6LL)*Power(xj,6LL) - 

                 15518LL*Power(r,7LL)*Power(xj,7LL) - 281LL*Power(r,8LL)*Power(xj,8LL) + 

                 9LL*Power(r,9LL)*Power(xj,9LL)) - 

              171LL*Power(xi,10LL)*Power(xj,16LL)*

               (-6768406575LL + 6280474725LL*r*xj + 

                 438336360LL*Power(r,2LL)*Power(xj,2LL) - 

                 400731030LL*Power(r,3LL)*Power(xj,3LL) - 

                 74168430LL*Power(r,4LL)*Power(xj,4LL) - 

                 2490810LL*Power(r,5LL)*Power(xj,5LL) + 

                 461160LL*Power(r,6LL)*Power(xj,6LL) + 

                 51244LL*Power(r,7LL)*Power(xj,7LL) + 1858LL*Power(r,8LL)*Power(xj,8LL) + 

                 18LL*Power(r,9LL)*Power(xj,9LL)) + 

              9LL*Power(xi,22LL)*Power(xj,4LL)*

               (-1346625LL - 2423925LL*r*xj - 2154600LL*Power(r,2LL)*Power(xj,2LL) - 

                 1256850LL*Power(r,3LL)*Power(xj,3LL) - 

                 538650LL*Power(r,4LL)*Power(xj,4LL) - 

                 179550LL*Power(r,5LL)*Power(xj,5LL) - 

                 47880LL*Power(r,6LL)*Power(xj,6LL) - 14264LL*Power(r,7LL)*Power(xj,7LL) + 

                 292LL*Power(r,8LL)*Power(xj,8LL) + 52LL*Power(r,9LL)*Power(xj,9LL)) - 

              9LL*Power(xi,4LL)*Power(xj,22LL)*

               (-129194933175LL - 73043543475LL*r*xj - 

                 17732214360LL*Power(r,2LL)*Power(xj,2LL) - 

                 2275149870LL*Power(r,3LL)*Power(xj,3LL) - 

                 134674470LL*Power(r,4LL)*Power(xj,4LL) + 

                 3148110LL*Power(r,5LL)*Power(xj,5LL) + 

                 1197000LL*Power(r,6LL)*Power(xj,6LL) + 

                 93176LL*Power(r,7LL)*Power(xj,7LL) + 3452LL*Power(r,8LL)*Power(xj,8LL) + 

                 52LL*Power(r,9LL)*Power(xj,9LL)) + 

              9LL*Power(xi,6LL)*Power(xj,20LL)*

               (356863797675LL + 115054179975LL*r*xj + 

                 3909863160LL*Power(r,2LL)*Power(xj,2LL) - 

                 3706015530LL*Power(r,3LL)*Power(xj,3LL) - 

                 798544530LL*Power(r,4LL)*Power(xj,4LL) - 

                 75669510LL*Power(r,5LL)*Power(xj,5LL) - 

                 3319400LL*Power(r,6LL)*Power(xj,6LL) - 

                 6456LL*Power(r,7LL)*Power(xj,7LL) + 5188LL*Power(r,8LL)*Power(xj,8LL) + 

                 148LL*Power(r,9LL)*Power(xj,9LL)) - 

              9LL*Power(xi,20LL)*Power(xj,6LL)*

               (-7630875LL - 13735575LL*r*xj - 12209400LL*Power(r,2LL)*Power(xj,2LL) - 

                 7122150LL*Power(r,3LL)*Power(xj,3LL) - 

                 3052350LL*Power(r,4LL)*Power(xj,4LL) - 

                 777210LL*Power(r,5LL)*Power(xj,5LL) - 

                 591640LL*Power(r,6LL)*Power(xj,6LL) + 3064LL*Power(r,7LL)*Power(xj,7LL) + 

                 5468LL*Power(r,8LL)*Power(xj,8LL) + 148LL*Power(r,9LL)*Power(xj,9LL)) + 

              2LL*Power(xi,18LL)*Power(xj,8LL)*

               (-137355750LL - 247240350LL*r*xj - 

                 219769200LL*Power(r,2LL)*Power(xj,2LL) - 

                 151171650LL*Power(r,3LL)*Power(xj,3LL) + 

                 13976550LL*Power(r,4LL)*Power(xj,4LL) - 

                 66692430LL*Power(r,5LL)*Power(xj,5LL) - 

                 1640520LL*Power(r,6LL)*Power(xj,6LL) + 

                 1046142LL*Power(r,7LL)*Power(xj,7LL) + 

                 66249LL*Power(r,8LL)*Power(xj,8LL) + 409LL*Power(r,9LL)*Power(xj,9LL)))))/

       (70875LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,19LL)*Power(xi + xj,18LL)) - 

      (141750LL*exp(2LL*r*(xi + xj))*(xi + xj)*

          Power(Power(xi,2LL) - Power(xj,2LL),19LL) + 

         exp(2LL*r*xj)*Power(xj,12LL)*

          (-5040LL*Power(r,7LL)*Power(xi,34LL) - 90LL*Power(r,8LL)*Power(xi,35LL) - 

            7740LL*Power(r,7LL)*Power(xi,32LL)*Power(xj,2LL) - 

            60LL*Power(r,8LL)*Power(xi,33LL)*Power(xj,2LL) + 127575LL*xi*Power(xj,26LL) + 

            226800LL*r*Power(xi,2LL)*Power(xj,26LL) + 

            132300LL*Power(r,2LL)*Power(xi,3LL)*Power(xj,26LL) - 

            210LL*Power(r,6LL)*Power(xi,33LL)*(630LL + Power(r,2LL)*Power(xj,2LL)) + 

            4725LL*Power(xi,3LL)*Power(xj,24LL)*(-513LL + 14LL*Power(r,2LL)*Power(xj,2LL)) - 

            540LL*Power(r,5LL)*Power(xi,32LL)*(3920LL + 43LL*Power(r,2LL)*Power(xj,2LL)) + 

            4725LL*r*Power(xi,5LL)*Power(xj,22LL)*

             (-532LL*r*Power(xj,2LL) + 8LL*Power(r,3LL)*Power(xj,4LL)) + 

            14175LL*Power(xi,4LL)*Power(xj,22LL)*

             (-304LL*r*Power(xj,2LL) + 8LL*Power(r,3LL)*Power(xj,4LL)) + 

            36LL*Power(r,5LL)*Power(xi,31LL)*

             (-9970LL*r*Power(xj,2LL) + 52LL*Power(r,3LL)*Power(xj,4LL)) + 

            36LL*Power(r,4LL)*Power(xi,30LL)*

             (-255920LL*r*Power(xj,2LL) + 3452LL*Power(r,3LL)*Power(xj,4LL)) + 

            4725LL*Power(xi,5LL)*Power(xj,22LL)*

             (4617LL - 266LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            180LL*Power(r,4LL)*Power(xi,31LL)*

             (-124950LL - 4985LL*Power(r,2LL)*Power(xj,2LL) + 

               13LL*Power(r,4LL)*Power(xj,4LL)) + 

            144LL*Power(r,3LL)*Power(xi,30LL)*

             (-1124550LL - 127960LL*Power(r,2LL)*Power(xj,2LL) + 

               863LL*Power(r,4LL)*Power(xj,4LL)) + 

            135LL*r*Power(xi,7LL)*Power(xj,20LL)*

             (167580LL*r*Power(xj,2LL) - 5320LL*Power(r,3LL)*Power(xj,4LL) + 

               24LL*Power(r,5LL)*Power(xj,6LL)) + 

            315LL*Power(xi,6LL)*Power(xj,20LL)*

             (123120LL*r*Power(xj,2LL) - 6840LL*Power(r,3LL)*Power(xj,4LL) + 

               48LL*Power(r,5LL)*Power(xj,6LL)) - 

            36LL*Power(r,3LL)*Power(xi,29LL)*

             (4205460LL*r*Power(xj,2LL) - 93176LL*Power(r,3LL)*Power(xj,4LL) + 

               222LL*Power(r,5LL)*Power(xj,6LL)) - 

            36LL*Power(r,2LL)*Power(xi,28LL)*

             (46802700LL*r*Power(xj,2LL) - 1197000LL*Power(r,3LL)*Power(xj,4LL) + 

               7782LL*Power(r,5LL)*Power(xj,6LL)) + 

            135LL*Power(xi,7LL)*Power(xj,20LL)*

             (-915705LL + 83790LL*Power(r,2LL)*Power(xj,2LL) - 

               1330LL*Power(r,4LL)*Power(xj,4LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) - 

            108LL*Power(r,2LL)*Power(xi,29LL)*

             (7122150LL + 2102730LL*Power(r,2LL)*Power(xj,2LL) - 

               23294LL*Power(r,4LL)*Power(xj,4LL) + 37LL*Power(r,6LL)*Power(xj,6LL)) - 

            72LL*r*Power(xi,28LL)*(30523500LL + 23401350LL*Power(r,2LL)*Power(xj,2LL) - 

               299250LL*Power(r,4LL)*Power(xj,4LL) + 1297LL*Power(r,6LL)*Power(xj,6LL)) + 

            r*Power(xi,17LL)*Power(xj,10LL)*

             (-43506975960LL*r*Power(xj,2LL) - 2983977360LL*Power(r,3LL)*Power(xj,4LL) + 

               31842936LL*Power(r,5LL)*Power(xj,6LL) - 6544LL*Power(r,7LL)*Power(xj,8LL)) + 

            10LL*r*Power(xi,9LL)*Power(xj,18LL)*

             (-12819870LL*r*Power(xj,2LL) + 646380LL*Power(r,3LL)*Power(xj,4LL) - 

               6156LL*Power(r,5LL)*Power(xj,6LL) + 8LL*Power(r,7LL)*Power(xj,8LL)) + 

            90LL*Power(xi,8LL)*Power(xj,18LL)*

             (-2441880LL*r*Power(xj,2LL) + 215460LL*Power(r,3LL)*Power(xj,4LL) - 

               3192LL*Power(r,5LL)*Power(xj,6LL) + 8LL*Power(r,7LL)*Power(xj,8LL)) - 

            1710LL*Power(xi,10LL)*Power(xj,16LL)*

             (-514080LL*r*Power(xj,2LL) + 64260LL*Power(r,3LL)*Power(xj,4LL) - 

               1512LL*Power(r,5LL)*Power(xj,6LL) + 8LL*Power(r,7LL)*Power(xj,8LL)) + 

            6LL*r*Power(xi,11LL)*Power(xj,16LL)*

             (100781100LL*r*Power(xj,2LL) - 4663260LL*Power(r,3LL)*Power(xj,4LL) + 

               128376LL*Power(r,5LL)*Power(xj,6LL) + 40LL*Power(r,7LL)*Power(xj,8LL)) - 

            18LL*r*Power(xi,13LL)*Power(xj,14LL)*

             (-311339700LL*r*Power(xj,2LL) - 29641080LL*Power(r,3LL)*Power(xj,4LL) - 

               9192LL*Power(r,5LL)*Power(xj,6LL) + 208LL*Power(r,7LL)*Power(xj,8LL)) + 

            18LL*r*Power(xi,15LL)*Power(xj,12LL)*

             (2658257700LL*r*Power(xj,2LL) - 30435720LL*Power(r,3LL)*Power(xj,4LL) - 

               697428LL*Power(r,5LL)*Power(xj,6LL) + 592LL*Power(r,7LL)*Power(xj,8LL)) - 

            18LL*Power(xi,12LL)*Power(xj,14LL)*

             (358142400LL*r*Power(xj,2LL) + 6211800LL*Power(r,3LL)*Power(xj,4LL) + 

               1774920LL*Power(r,5LL)*Power(xj,6LL) + 1168LL*Power(r,7LL)*Power(xj,8LL)) + 

            18LL*r*Power(xi,25LL)*Power(xj,2LL)*

             (-2275149870LL*r*Power(xj,2LL) + 151339020LL*Power(r,3LL)*Power(xj,4LL) - 

               1638372LL*Power(r,5LL)*Power(xj,6LL) + 1368LL*Power(r,7LL)*Power(xj,8LL)) - 

            9LL*r*Power(xi,19LL)*Power(xj,8LL)*

             (15227779140LL*r*Power(xj,2LL) - 682098480LL*Power(r,3LL)*Power(xj,4LL) + 

               2387616LL*Power(r,5LL)*Power(xj,6LL) + 2736LL*Power(r,7LL)*Power(xj,8LL)) - 

            3LL*r*Power(xi,23LL)*Power(xj,4LL)*

             (-22236093180LL*r*Power(xj,2LL) + 1310447880LL*Power(r,3LL)*Power(xj,4LL) - 

               17525448LL*Power(r,5LL)*Power(xj,6LL) + 20672LL*Power(r,7LL)*Power(xj,8LL)) 

    + 3LL*r*Power(xi,21LL)*Power(xj,6LL)*

             (38061529380LL*r*Power(xj,2LL) - 567904680LL*Power(r,3LL)*Power(xj,4LL) - 

               8651232LL*Power(r,5LL)*Power(xj,6LL) + 20672LL*Power(r,7LL)*Power(xj,8LL)) + 

            63LL*Power(xi,20LL)*Power(xj,6LL)*

             (12480405840LL*r*Power(xj,2LL) - 805257240LL*Power(r,3LL)*Power(xj,4LL) + 

               5736480LL*Power(r,5LL)*Power(xj,6LL) + 20672LL*Power(r,7LL)*Power(xj,8LL)) + 

            18LL*Power(xi,14LL)*Power(xj,12LL)*

             (-5038412400LL*r*Power(xj,2LL) - 478879800LL*Power(r,3LL)*Power(xj,4LL) + 

               1093680LL*Power(r,5LL)*Power(xj,6LL) + 21872LL*Power(r,7LL)*Power(xj,8LL)) - 

            18LL*Power(xi,26LL)*(3589563600LL*r*Power(xj,2LL) + 

               269348940LL*Power(r,3LL)*Power(xj,4LL) - 

               9958200LL*Power(r,5LL)*Power(xj,6LL) + 32712LL*Power(r,7LL)*Power(xj,8LL)) + 

            9LL*Power(xi,18LL)*Power(xj,8LL)*

             (16656781680LL*r*Power(xj,2LL) + 3648168720LL*Power(r,3LL)*Power(xj,4LL) - 

               77054880LL*Power(r,5LL)*Power(xj,6LL) + 85424LL*Power(r,7LL)*Power(xj,8LL)) 

    - 9LL*Power(xi,16LL)*Power(xj,10LL)*(42596154720LL*r*Power(xj,2LL) - 

               1071714000LL*Power(r,3LL)*Power(xj,4LL) - 

               32749920LL*Power(r,5LL)*Power(xj,6LL) + 117776LL*Power(r,7LL)*Power(xj,8LL)) 

    + 18LL*Power(xi,24LL)*Power(xj,2LL)*(-17732214360LL*r*Power(xj,2LL) + 

               1597089060LL*Power(r,3LL)*Power(xj,4LL) - 

               35554680LL*Power(r,5LL)*Power(xj,6LL) + 141208LL*Power(r,7LL)*Power(xj,8LL)) 

    - 9LL*Power(xi,22LL)*Power(xj,4LL)*(7819726320LL*r*Power(xj,2LL) + 

               365683080LL*Power(r,3LL)*Power(xj,4LL) - 

               52572240LL*Power(r,5LL)*Power(xj,6LL) + 351424LL*Power(r,7LL)*Power(xj,8LL)) 

    + Power(xi,17LL)*Power(xj,10LL)*(1073961177975LL - 

               21753487980LL*Power(r,2LL)*Power(xj,2LL) - 

               745994340LL*Power(r,4LL)*Power(xj,4LL) + 

               5307156LL*Power(r,6LL)*Power(xj,6LL) - 818LL*Power(r,8LL)*Power(xj,8LL)) + 

            10LL*Power(xi,9LL)*Power(xj,18LL)*

             (49448070LL - 6409935LL*Power(r,2LL)*Power(xj,2LL) + 

               161595LL*Power(r,4LL)*Power(xj,4LL) - 1026LL*Power(r,6LL)*Power(xj,6LL) + 

               Power(r,8LL)*Power(xj,8LL)) + 

            6LL*Power(xi,11LL)*Power(xj,16LL)*

             (-207559800LL + 50390550LL*Power(r,2LL)*Power(xj,2LL) - 

               1165815LL*Power(r,4LL)*Power(xj,4LL) + 21396LL*Power(r,6LL)*Power(xj,6LL) + 

               5LL*Power(r,8LL)*Power(xj,8LL)) - 

            18LL*Power(xi,13LL)*Power(xj,14LL)*

             (-1703720025LL - 155669850LL*Power(r,2LL)*Power(xj,2LL) - 

               7410270LL*Power(r,4LL)*Power(xj,4LL) - 1532LL*Power(r,6LL)*Power(xj,6LL) + 

               26LL*Power(r,8LL)*Power(xj,8LL)) + 

            18LL*Power(xi,15LL)*Power(xj,12LL)*

             (19380896325LL + 1329128850LL*Power(r,2LL)*Power(xj,2LL) - 

               7608930LL*Power(r,4LL)*Power(xj,4LL) - 

               116238LL*Power(r,6LL)*Power(xj,6LL) + 74LL*Power(r,8LL)*Power(xj,8LL)) + 

            18LL*Power(xi,25LL)*Power(xj,2LL)*

             (-5449970925LL - 1137574935LL*Power(r,2LL)*Power(xj,2LL) + 

               37834755LL*Power(r,4LL)*Power(xj,4LL) - 

               273062LL*Power(r,6LL)*Power(xj,6LL) + 171LL*Power(r,8LL)*Power(xj,8LL)) - 

            9LL*Power(xi,19LL)*Power(xj,8LL)*

             (-37914907275LL + 7613889570LL*Power(r,2LL)*Power(xj,2LL) - 

               170524620LL*Power(r,4LL)*Power(xj,4LL) + 

               397936LL*Power(r,6LL)*Power(xj,6LL) + 342LL*Power(r,8LL)*Power(xj,8LL)) - 

            3LL*Power(xi,23LL)*Power(xj,4LL)*

             (219130630425LL - 11118046590LL*Power(r,2LL)*Power(xj,2LL) + 

               327611970LL*Power(r,4LL)*Power(xj,4LL) - 

               2920908LL*Power(r,6LL)*Power(xj,6LL) + 2584LL*Power(r,8LL)*Power(xj,8LL)) + 

            3LL*Power(xi,21LL)*Power(xj,6LL)*

             (-345162539925LL + 19030764690LL*Power(r,2LL)*Power(xj,2LL) - 

               141976170LL*Power(r,4LL)*Power(xj,4LL) - 

               1441872LL*Power(r,6LL)*Power(xj,6LL) + 2584LL*Power(r,8LL)*Power(xj,8LL)) + 

            Power(xi,27LL)*(-2884470750LL - 19229805000LL*Power(r,2LL)*Power(xj,2LL) + 

               141664950LL*Power(r,4LL)*Power(xj,4LL) + 

               406728LL*Power(r,6LL)*Power(xj,6LL) + 7362LL*Power(r,8LL)*Power(xj,8LL))) + 

         2LL*exp(2LL*r*xj)*Power(xj,13LL)*

          (-630LL*Power(r,8LL)*Power(xi,34LL) - 10LL*Power(r,9LL)*Power(xi,35LL) + 

            70875LL*Power(xj,26LL) + 127575LL*r*xi*Power(xj,26LL) - 

            30LL*Power(r,7LL)*Power(xi,33LL)*(630LL + Power(r,2LL)*Power(xj,2LL)) + 

            14175LL*Power(xi,2LL)*Power(xj,24LL)*(-95LL + 8LL*Power(r,2LL)*Power(xj,2LL)) + 

            4725LL*r*Power(xi,3LL)*Power(xj,24LL)*

             (-513LL + 14LL*Power(r,2LL)*Power(xj,2LL)) - 

            90LL*Power(r,6LL)*Power(xi,32LL)*(3920LL + 43LL*Power(r,2LL)*Power(xj,2LL)) + 

            4725LL*r*Power(xi,5LL)*Power(xj,22LL)*

             (4617LL - 266LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            14175LL*Power(xi,4LL)*Power(xj,22LL)*

             (855LL - 152LL*Power(r,2LL)*Power(xj,2LL) + 2LL*Power(r,4LL)*Power(xj,4LL)) + 

            36LL*Power(r,5LL)*Power(xi,31LL)*

             (-124950LL - 4985LL*Power(r,2LL)*Power(xj,2LL) + 

               13LL*Power(r,4LL)*Power(xj,4LL)) + 

            36LL*Power(r,4LL)*Power(xi,30LL)*

             (-1124550LL - 127960LL*Power(r,2LL)*Power(xj,2LL) + 

               863LL*Power(r,4LL)*Power(xj,4LL)) + 

            135LL*r*Power(xi,7LL)*Power(xj,20LL)*

             (-915705LL + 83790LL*Power(r,2LL)*Power(xj,2LL) - 

               1330LL*Power(r,4LL)*Power(xj,4LL) + 4LL*Power(r,6LL)*Power(xj,6LL)) + 

            315LL*Power(xi,6LL)*Power(xj,20LL)*

             (-218025LL + 61560LL*Power(r,2LL)*Power(xj,2LL) - 

               1710LL*Power(r,4LL)*Power(xj,4LL) + 8LL*Power(r,6LL)*Power(xj,6LL)) - 

            36LL*Power(r,3LL)*Power(xi,29LL)*

             (7122150LL + 2102730LL*Power(r,2LL)*Power(xj,2LL) - 

               23294LL*Power(r,4LL)*Power(xj,4LL) + 37LL*Power(r,6LL)*Power(xj,6LL)) - 

            36LL*Power(r,2LL)*Power(xi,28LL)*

             (30523500LL + 23401350LL*Power(r,2LL)*Power(xj,2LL) - 

               299250LL*Power(r,4LL)*Power(xj,4LL) + 1297LL*Power(r,6LL)*Power(xj,6LL)) + 

            r*Power(xi,17LL)*Power(xj,10LL)*

             (1073961177975LL - 21753487980LL*Power(r,2LL)*Power(xj,2LL) - 

               745994340LL*Power(r,4LL)*Power(xj,4LL) + 

               5307156LL*Power(r,6LL)*Power(xj,6LL) - 818LL*Power(r,8LL)*Power(xj,8LL)) + 

            10LL*r*Power(xi,9LL)*Power(xj,18LL)*

             (49448070LL - 6409935LL*Power(r,2LL)*Power(xj,2LL) + 

               161595LL*Power(r,4LL)*Power(xj,4LL) - 1026LL*Power(r,6LL)*Power(xj,6LL) + 

               Power(r,8LL)*Power(xj,8LL)) + 

            90LL*Power(xi,8LL)*Power(xj,18LL)*

             (3052350LL - 1220940LL*Power(r,2LL)*Power(xj,2LL) + 

               53865LL*Power(r,4LL)*Power(xj,4LL) - 532LL*Power(r,6LL)*Power(xj,6LL) + 

               Power(r,8LL)*Power(xj,8LL)) - 

            1710LL*Power(xi,10LL)*Power(xj,16LL)*

             (481950LL - 257040LL*Power(r,2LL)*Power(xj,2LL) + 

               16065LL*Power(r,4LL)*Power(xj,4LL) - 252LL*Power(r,6LL)*Power(xj,6LL) + 

               Power(r,8LL)*Power(xj,8LL)) + 

            6LL*r*Power(xi,11LL)*Power(xj,16LL)*

             (-207559800LL + 50390550LL*Power(r,2LL)*Power(xj,2LL) - 

               1165815LL*Power(r,4LL)*Power(xj,4LL) + 21396LL*Power(r,6LL)*Power(xj,6LL) + 

               5LL*Power(r,8LL)*Power(xj,8LL)) - 

            18LL*r*Power(xi,13LL)*Power(xj,14LL)*

             (-1703720025LL - 155669850LL*Power(r,2LL)*Power(xj,2LL) - 

               7410270LL*Power(r,4LL)*Power(xj,4LL) - 1532LL*Power(r,6LL)*Power(xj,6LL) + 

               26LL*Power(r,8LL)*Power(xj,8LL)) + 

            18LL*r*Power(xi,15LL)*Power(xj,12LL)*

             (19380896325LL + 1329128850LL*Power(r,2LL)*Power(xj,2LL) - 

               7608930LL*Power(r,4LL)*Power(xj,4LL) - 

               116238LL*Power(r,6LL)*Power(xj,6LL) + 74LL*Power(r,8LL)*Power(xj,8LL)) - 

            18LL*Power(xi,12LL)*Power(xj,14LL)*

             (89026875LL + 179071200LL*Power(r,2LL)*Power(xj,2LL) + 

               1552950LL*Power(r,4LL)*Power(xj,4LL) + 

               295820LL*Power(r,6LL)*Power(xj,6LL) + 146LL*Power(r,8LL)*Power(xj,8LL)) + 

            18LL*r*Power(xi,25LL)*Power(xj,2LL)*

             (-5449970925LL - 1137574935LL*Power(r,2LL)*Power(xj,2LL) + 

               37834755LL*Power(r,4LL)*Power(xj,4LL) - 

               273062LL*Power(r,6LL)*Power(xj,6LL) + 171LL*Power(r,8LL)*Power(xj,8LL)) - 

            9LL*r*Power(xi,19LL)*Power(xj,8LL)*

             (-37914907275LL + 7613889570LL*Power(r,2LL)*Power(xj,2LL) - 

               170524620LL*Power(r,4LL)*Power(xj,4LL) + 

               397936LL*Power(r,6LL)*Power(xj,6LL) + 342LL*Power(r,8LL)*Power(xj,8LL)) - 

            3LL*r*Power(xi,23LL)*Power(xj,4LL)*

             (219130630425LL - 11118046590LL*Power(r,2LL)*Power(xj,2LL) + 

               327611970LL*Power(r,4LL)*Power(xj,4LL) - 

               2920908LL*Power(r,6LL)*Power(xj,6LL) + 2584LL*Power(r,8LL)*Power(xj,8LL)) + 

            3LL*r*Power(xi,21LL)*Power(xj,6LL)*

             (-345162539925LL + 19030764690LL*Power(r,2LL)*Power(xj,2LL) - 

               141976170LL*Power(r,4LL)*Power(xj,4LL) - 

               1441872LL*Power(r,6LL)*Power(xj,6LL) + 2584LL*Power(r,8LL)*Power(xj,8LL)) + 

            63LL*Power(xi,20LL)*Power(xj,6LL)*

             (-50980542525LL + 6240202920LL*Power(r,2LL)*Power(xj,2LL) - 

               201314310LL*Power(r,4LL)*Power(xj,4LL) + 

               956080LL*Power(r,6LL)*Power(xj,6LL) + 2584LL*Power(r,8LL)*Power(xj,8LL)) + 

            18LL*Power(xi,14LL)*Power(xj,12LL)*

             (-7803332775LL - 2519206200LL*Power(r,2LL)*Power(xj,2LL) - 

               119719950LL*Power(r,4LL)*Power(xj,4LL) + 

               182280LL*Power(r,6LL)*Power(xj,6LL) + 2734LL*Power(r,8LL)*Power(xj,8LL)) - 

            18LL*Power(xi,26LL)*(195859125LL + 1794781800LL*Power(r,2LL)*Power(xj,2LL) + 

               67337235LL*Power(r,4LL)*Power(xj,4LL) - 

               1659700LL*Power(r,6LL)*Power(xj,6LL) + 4089LL*Power(r,8LL)*Power(xj,8LL)) + 

            9LL*Power(xi,18LL)*Power(xj,8LL)*

             (-357591274425LL + 8328390840LL*Power(r,2LL)*Power(xj,2LL) + 

               912042180LL*Power(r,4LL)*Power(xj,4LL) - 

               12842480LL*Power(r,6LL)*Power(xj,6LL) + 10678LL*Power(r,8LL)*Power(xj,8LL)) 

    - 9LL*Power(xi,16LL)*Power(xj,10LL)*(128599724925LL + 

               21298077360LL*Power(r,2LL)*Power(xj,2LL) - 

               267928500LL*Power(r,4LL)*Power(xj,4LL) - 

               5458320LL*Power(r,6LL)*Power(xj,6LL) + 14722LL*Power(r,8LL)*Power(xj,8LL)) + 

            18LL*Power(xi,24LL)*Power(xj,2LL)*

             (-7604930025LL - 8866107180LL*Power(r,2LL)*Power(xj,2LL) + 

               399272265LL*Power(r,4LL)*Power(xj,4LL) - 

               5925780LL*Power(r,6LL)*Power(xj,6LL) + 17651LL*Power(r,8LL)*Power(xj,8LL)) - 

            9LL*Power(xi,22LL)*Power(xj,4LL)*

             (129194933175LL + 3909863160LL*Power(r,2LL)*Power(xj,2LL) + 

               91420770LL*Power(r,4LL)*Power(xj,4LL) - 

               8762040LL*Power(r,6LL)*Power(xj,6LL) + 43928LL*Power(r,8LL)*Power(xj,8LL)) + 

            Power(xi,27LL)*(-2884470750LL*r - 6409935000LL*Power(r,3LL)*Power(xj,2LL) + 

               28332990LL*Power(r,5LL)*Power(xj,4LL) + 58104LL*Power(r,7LL)*Power(xj,6LL) + 

               818LL*Power(r,9LL)*Power(xj,8LL))) + 

         exp(2LL*r*xi)*Power(xi,12LL)*

          (Power(xi,8LL)*Power(xj,18LL)*

             (-341234165475LL*xj - 786265567920LL*r*Power(xj,2LL) - 

               171276882210LL*Power(r,2LL)*Power(xj,3LL) + 

               3291147720LL*Power(r,3LL)*Power(xj,4LL) + 

               4914179550LL*Power(r,4LL)*Power(xj,5LL) + 

               639984240LL*Power(r,5LL)*Power(xj,6LL) + 

               34405812LL*Power(r,6LL)*Power(xj,7LL) + 

               588816LL*Power(r,7LL)*Power(xj,8LL) - 7362LL*Power(r,8LL)*Power(xj,9LL)) + 

            10LL*Power(xj,26LL)*(288447075LL*xj + 219769200LL*r*Power(xj,2LL) + 

               76919220LL*Power(r,2LL)*Power(xj,3LL) + 

               16193520LL*Power(r,3LL)*Power(xj,4LL) + 

               2249100LL*Power(r,4LL)*Power(xj,5LL) + 

               211680LL*Power(r,5LL)*Power(xj,6LL) + 13230LL*Power(r,6LL)*Power(xj,7LL) + 

               504LL*Power(r,7LL)*Power(xj,8LL) + 9LL*Power(r,8LL)*Power(xj,9LL)) + 

            30LL*Power(xi,2LL)*Power(xj,24LL)*

             (3269982555LL*xj + 2153738160LL*r*Power(xj,2LL) + 

               640993500LL*Power(r,2LL)*Power(xj,3LL) + 

               112326480LL*Power(r,3LL)*Power(xj,4LL) + 

               12616380LL*Power(r,4LL)*Power(xj,5LL) + 

               921312LL*Power(r,5LL)*Power(xj,6LL) + 41874LL*Power(r,6LL)*Power(xj,7LL) + 

               1032LL*Power(r,7LL)*Power(xj,8LL) + 9LL*Power(r,8LL)*Power(xj,9LL)) - 

            15LL*Power(xi,24LL)*Power(xj,2LL)*

             (-161595LL*xj - 287280LL*r*Power(xj,2LL) - 

               251370LL*Power(r,2LL)*Power(xj,3LL) - 143640LL*Power(r,3LL)*Power(xj,4LL) - 

               59850LL*Power(r,4LL)*Power(xj,5LL) - 19152LL*Power(r,5LL)*Power(xj,6LL) - 

               4788LL*Power(r,6LL)*Power(xj,7LL) - 912LL*Power(r,7LL)*Power(xj,8LL) + 

               18LL*Power(r,8LL)*Power(xj,9LL)) - 

            5LL*Power(xi,26LL)*(25515LL*xj + 45360LL*r*Power(xj,2LL) + 

               39690LL*Power(r,2LL)*Power(xj,3LL) + 22680LL*Power(r,3LL)*Power(xj,4LL) + 

               9450LL*Power(r,4LL)*Power(xj,5LL) + 3024LL*Power(r,5LL)*Power(xj,6LL) + 

               756LL*Power(r,6LL)*Power(xj,7LL) + 144LL*Power(r,7LL)*Power(xj,8LL) + 

               18LL*Power(r,8LL)*Power(xj,9LL)) - 

            1938LL*Power(xi,14LL)*Power(xj,12LL)*

             (15824025LL*xj - 46796400LL*r*Power(xj,2LL) + 

               37034550LL*Power(r,2LL)*Power(xj,3LL) + 

               4977000LL*Power(r,3LL)*Power(xj,4LL) - 

               1924650LL*Power(r,4LL)*Power(xj,5LL) - 

               357840LL*Power(r,5LL)*Power(xj,6LL) - 12936LL*Power(r,6LL)*Power(xj,7LL) + 

               672LL*Power(r,7LL)*Power(xj,8LL) + 36LL*Power(r,8LL)*Power(xj,9LL)) + 

            1938LL*Power(xi,12LL)*Power(xj,14LL)*

             (-180008325LL*xj + 197814960LL*r*Power(xj,2LL) + 

               33674130LL*Power(r,2LL)*Power(xj,3LL) - 

               16941960LL*Power(r,3LL)*Power(xj,4LL) - 

               3959550LL*Power(r,4LL)*Power(xj,5LL) - 

               186480LL*Power(r,5LL)*Power(xj,6LL) + 15624LL*Power(r,6LL)*Power(xj,7LL) + 

               1632LL*Power(r,7LL)*Power(xj,8LL) + 36LL*Power(r,8LL)*Power(xj,9LL)) + 

            342LL*Power(xi,16LL)*Power(xj,10LL)*

             (3641400LL*xj + 18849600LL*r*Power(xj,2LL) - 

               24579450LL*Power(r,2LL)*Power(xj,3LL) + 

               25204200LL*Power(r,3LL)*Power(xj,4LL) + 

               2002350LL*Power(r,4LL)*Power(xj,5LL) - 

               861840LL*Power(r,5LL)*Power(xj,6LL) - 108626LL*Power(r,6LL)*Power(xj,7LL) - 

               2248LL*Power(r,7LL)*Power(xj,8LL) + 81LL*Power(r,8LL)*Power(xj,9LL)) - 

            171LL*Power(xi,10LL)*Power(xj,16LL)*

             (6280474725LL*xj + 876672720LL*r*Power(xj,2LL) - 

               1202193090LL*Power(r,2LL)*Power(xj,3LL) - 

               296673720LL*Power(r,3LL)*Power(xj,4LL) - 

               12454050LL*Power(r,4LL)*Power(xj,5LL) + 

               2766960LL*Power(r,5LL)*Power(xj,6LL) + 

               358708LL*Power(r,6LL)*Power(xj,7LL) + 14864LL*Power(r,7LL)*Power(xj,8LL) + 

               162LL*Power(r,8LL)*Power(xj,9LL)) + 

            9LL*Power(xi,22LL)*Power(xj,4LL)*

             (-2423925LL*xj - 4309200LL*r*Power(xj,2LL) - 

               3770550LL*Power(r,2LL)*Power(xj,3LL) - 

               2154600LL*Power(r,3LL)*Power(xj,4LL) - 

               897750LL*Power(r,4LL)*Power(xj,5LL) - 287280LL*Power(r,5LL)*Power(xj,6LL) - 

               99848LL*Power(r,6LL)*Power(xj,7LL) + 2336LL*Power(r,7LL)*Power(xj,8LL) + 

               468LL*Power(r,8LL)*Power(xj,9LL)) - 

            9LL*Power(xi,4LL)*Power(xj,22LL)*

             (-73043543475LL*xj - 35464428720LL*r*Power(xj,2LL) - 

               6825449610LL*Power(r,2LL)*Power(xj,3LL) - 

               538697880LL*Power(r,3LL)*Power(xj,4LL) + 

               15740550LL*Power(r,4LL)*Power(xj,5LL) + 

               7182000LL*Power(r,5LL)*Power(xj,6LL) + 

               652232LL*Power(r,6LL)*Power(xj,7LL) + 27616LL*Power(r,7LL)*Power(xj,8LL) + 

               468LL*Power(r,8LL)*Power(xj,9LL)) + 

            9LL*Power(xi,6LL)*Power(xj,20LL)*

             (115054179975LL*xj + 7819726320LL*r*Power(xj,2LL) - 

               11118046590LL*Power(r,2LL)*Power(xj,3LL) - 

               3194178120LL*Power(r,3LL)*Power(xj,4LL) - 

               378347550LL*Power(r,4LL)*Power(xj,5LL) - 

               19916400LL*Power(r,5LL)*Power(xj,6LL) - 

               45192LL*Power(r,6LL)*Power(xj,7LL) + 41504LL*Power(r,7LL)*Power(xj,8LL) + 

               1332LL*Power(r,8LL)*Power(xj,9LL)) - 

            9LL*Power(xi,20LL)*Power(xj,6LL)*

             (-13735575LL*xj - 24418800LL*r*Power(xj,2LL) - 

               21366450LL*Power(r,2LL)*Power(xj,3LL) - 

               12209400LL*Power(r,3LL)*Power(xj,4LL) - 

               3886050LL*Power(r,4LL)*Power(xj,5LL) - 

               3549840LL*Power(r,5LL)*Power(xj,6LL) + 21448LL*Power(r,6LL)*Power(xj,7LL) + 

               43744LL*Power(r,7LL)*Power(xj,8LL) + 1332LL*Power(r,8LL)*Power(xj,9LL)) + 

            2LL*Power(xi,18LL)*Power(xj,8LL)*

             (-247240350LL*xj - 439538400LL*r*Power(xj,2LL) - 

               453514950LL*Power(r,2LL)*Power(xj,3LL) + 

               55906200LL*Power(r,3LL)*Power(xj,4LL) - 

               333462150LL*Power(r,4LL)*Power(xj,5LL) - 

               9843120LL*Power(r,5LL)*Power(xj,6LL) + 

               7322994LL*Power(r,6LL)*Power(xj,7LL) + 529992LL*Power(r,7LL)*Power(xj,8LL) + 

               3681LL*Power(r,8LL)*Power(xj,9LL))) + 

         2LL*exp(2LL*r*xi)*Power(xi,13LL)*

          (Power(xi,8LL)*Power(xj,18LL)*

             (3218321469825LL - 341234165475LL*r*xj - 

               393132783960LL*Power(r,2LL)*Power(xj,2LL) - 

               57092294070LL*Power(r,3LL)*Power(xj,3LL) + 

               822786930LL*Power(r,4LL)*Power(xj,4LL) + 

               982835910LL*Power(r,5LL)*Power(xj,5LL) + 

               106664040LL*Power(r,6LL)*Power(xj,6LL) + 

               4915116LL*Power(r,7LL)*Power(xj,7LL) + 73602LL*Power(r,8LL)*Power(xj,8LL) - 

               818LL*Power(r,9LL)*Power(xj,9LL)) + 

            10LL*Power(xj,26LL)*(352546425LL + 288447075LL*r*xj + 

               109884600LL*Power(r,2LL)*Power(xj,2LL) + 

               25639740LL*Power(r,3LL)*Power(xj,3LL) + 

               4048380LL*Power(r,4LL)*Power(xj,4LL) + 449820LL*Power(r,5LL)*Power(xj,5LL) + 

               35280LL*Power(r,6LL)*Power(xj,6LL) + 1890LL*Power(r,7LL)*Power(xj,7LL) + 

               63LL*Power(r,8LL)*Power(xj,8LL) + Power(r,9LL)*Power(xj,9LL)) + 

            30LL*Power(xi,2LL)*Power(xj,24LL)*

             (4562958015LL + 3269982555LL*r*xj + 

               1076869080LL*Power(r,2LL)*Power(xj,2LL) + 

               213664500LL*Power(r,3LL)*Power(xj,3LL) + 

               28081620LL*Power(r,4LL)*Power(xj,4LL) + 

               2523276LL*Power(r,5LL)*Power(xj,5LL) + 153552LL*Power(r,6LL)*Power(xj,6LL) + 

               5982LL*Power(r,7LL)*Power(xj,7LL) + 129LL*Power(r,8LL)*Power(xj,8LL) + 

               Power(r,9LL)*Power(xj,9LL)) - 

            15LL*Power(xi,24LL)*Power(xj,2LL)*

             (-89775LL - 161595LL*r*xj - 143640LL*Power(r,2LL)*Power(xj,2LL) - 

               83790LL*Power(r,3LL)*Power(xj,3LL) - 35910LL*Power(r,4LL)*Power(xj,4LL) - 

               11970LL*Power(r,5LL)*Power(xj,5LL) - 3192LL*Power(r,6LL)*Power(xj,6LL) - 

               684LL*Power(r,7LL)*Power(xj,7LL) - 114LL*Power(r,8LL)*Power(xj,8LL) + 

               2LL*Power(r,9LL)*Power(xj,9LL)) - 

            5LL*Power(xi,26LL)*(14175LL + 25515LL*r*xj + 22680LL*Power(r,2LL)*Power(xj,2LL) + 

               13230LL*Power(r,3LL)*Power(xj,3LL) + 5670LL*Power(r,4LL)*Power(xj,4LL) + 

               1890LL*Power(r,5LL)*Power(xj,5LL) + 504LL*Power(r,6LL)*Power(xj,6LL) + 

               108LL*Power(r,7LL)*Power(xj,7LL) + 18LL*Power(r,8LL)*Power(xj,8LL) + 

               2LL*Power(r,9LL)*Power(xj,9LL)) - 

            1938LL*Power(xi,14LL)*Power(xj,12LL)*

             (-826875LL + 15824025LL*r*xj - 23398200LL*Power(r,2LL)*Power(xj,2LL) + 

               12344850LL*Power(r,3LL)*Power(xj,3LL) + 

               1244250LL*Power(r,4LL)*Power(xj,4LL) - 384930LL*Power(r,5LL)*Power(xj,5LL) - 

               59640LL*Power(r,6LL)*Power(xj,6LL) - 1848LL*Power(r,7LL)*Power(xj,7LL) + 

               84LL*Power(r,8LL)*Power(xj,8LL) + 4LL*Power(r,9LL)*Power(xj,9LL)) + 

            1938LL*Power(xi,12LL)*Power(xj,14LL)*

             (72476775LL - 180008325LL*r*xj + 98907480LL*Power(r,2LL)*Power(xj,2LL) + 

               11224710LL*Power(r,3LL)*Power(xj,3LL) - 

               4235490LL*Power(r,4LL)*Power(xj,4LL) - 791910LL*Power(r,5LL)*Power(xj,5LL) - 

               31080LL*Power(r,6LL)*Power(xj,6LL) + 2232LL*Power(r,7LL)*Power(xj,7LL) + 

               204LL*Power(r,8LL)*Power(xj,8LL) + 4LL*Power(r,9LL)*Power(xj,9LL)) + 

            342LL*Power(xi,16LL)*Power(xj,10LL)*

             (2409750LL + 3641400LL*r*xj + 9424800LL*Power(r,2LL)*Power(xj,2LL) - 

               8193150LL*Power(r,3LL)*Power(xj,3LL) + 

               6301050LL*Power(r,4LL)*Power(xj,4LL) + 400470LL*Power(r,5LL)*Power(xj,5LL) - 

               143640LL*Power(r,6LL)*Power(xj,6LL) - 15518LL*Power(r,7LL)*Power(xj,7LL) - 

               281LL*Power(r,8LL)*Power(xj,8LL) + 9LL*Power(r,9LL)*Power(xj,9LL)) - 

            171LL*Power(xi,10LL)*Power(xj,16LL)*

             (-6768406575LL + 6280474725LL*r*xj + 

               438336360LL*Power(r,2LL)*Power(xj,2LL) - 

               400731030LL*Power(r,3LL)*Power(xj,3LL) - 

               74168430LL*Power(r,4LL)*Power(xj,4LL) - 

               2490810LL*Power(r,5LL)*Power(xj,5LL) + 461160LL*Power(r,6LL)*Power(xj,6LL) + 

               51244LL*Power(r,7LL)*Power(xj,7LL) + 1858LL*Power(r,8LL)*Power(xj,8LL) + 

               18LL*Power(r,9LL)*Power(xj,9LL)) + 

            9LL*Power(xi,22LL)*Power(xj,4LL)*

             (-1346625LL - 2423925LL*r*xj - 2154600LL*Power(r,2LL)*Power(xj,2LL) - 

               1256850LL*Power(r,3LL)*Power(xj,3LL) - 538650LL*Power(r,4LL)*Power(xj,4LL) - 

               179550LL*Power(r,5LL)*Power(xj,5LL) - 47880LL*Power(r,6LL)*Power(xj,6LL) - 

               14264LL*Power(r,7LL)*Power(xj,7LL) + 292LL*Power(r,8LL)*Power(xj,8LL) + 

               52LL*Power(r,9LL)*Power(xj,9LL)) - 

            9LL*Power(xi,4LL)*Power(xj,22LL)*

             (-129194933175LL - 73043543475LL*r*xj - 

               17732214360LL*Power(r,2LL)*Power(xj,2LL) - 

               2275149870LL*Power(r,3LL)*Power(xj,3LL) - 

               134674470LL*Power(r,4LL)*Power(xj,4LL) + 

               3148110LL*Power(r,5LL)*Power(xj,5LL) + 

               1197000LL*Power(r,6LL)*Power(xj,6LL) + 93176LL*Power(r,7LL)*Power(xj,7LL) + 

               3452LL*Power(r,8LL)*Power(xj,8LL) + 52LL*Power(r,9LL)*Power(xj,9LL)) + 

            9LL*Power(xi,6LL)*Power(xj,20LL)*

             (356863797675LL + 115054179975LL*r*xj + 

               3909863160LL*Power(r,2LL)*Power(xj,2LL) - 

               3706015530LL*Power(r,3LL)*Power(xj,3LL) - 

               798544530LL*Power(r,4LL)*Power(xj,4LL) - 

               75669510LL*Power(r,5LL)*Power(xj,5LL) - 

               3319400LL*Power(r,6LL)*Power(xj,6LL) - 6456LL*Power(r,7LL)*Power(xj,7LL) + 

               5188LL*Power(r,8LL)*Power(xj,8LL) + 148LL*Power(r,9LL)*Power(xj,9LL)) - 

            9LL*Power(xi,20LL)*Power(xj,6LL)*

             (-7630875LL - 13735575LL*r*xj - 12209400LL*Power(r,2LL)*Power(xj,2LL) - 

               7122150LL*Power(r,3LL)*Power(xj,3LL) - 

               3052350LL*Power(r,4LL)*Power(xj,4LL) - 777210LL*Power(r,5LL)*Power(xj,5LL) - 

               591640LL*Power(r,6LL)*Power(xj,6LL) + 3064LL*Power(r,7LL)*Power(xj,7LL) + 

               5468LL*Power(r,8LL)*Power(xj,8LL) + 148LL*Power(r,9LL)*Power(xj,9LL)) + 

            2LL*Power(xi,18LL)*Power(xj,8LL)*

             (-137355750LL - 247240350LL*r*xj - 219769200LL*Power(r,2LL)*Power(xj,2LL) - 

               151171650LL*Power(r,3LL)*Power(xj,3LL) + 

               13976550LL*Power(r,4LL)*Power(xj,4LL) - 

               66692430LL*Power(r,5LL)*Power(xj,5LL) - 

               1640520LL*Power(r,6LL)*Power(xj,6LL) + 1046142LL*Power(r,7LL)*Power(xj,7LL) + 

               66249LL*Power(r,8LL)*Power(xj,8LL) + 409LL*Power(r,9LL)*Power(xj,9LL))))/

       (70875LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj,19LL)*Power(xi + xj,19LL))

    ; }
   
  }
  return S;
}

cl_R DSlater_5S_1S(cl_R r,cl_R xi,cl_R xj)
{
  return DSlater_1S_5S(r,xj,xi);
}

cl_R DSlater_5S_2S(cl_R r,cl_R xi,cl_R xj)
{
  return DSlater_2S_5S(r,xj,xi);
}

cl_R DSlater_5S_3S(cl_R r,cl_R xi,cl_R xj)
{
  return DSlater_3S_5S(r,xj,xi);
}

cl_R DSlater_5S_4S(cl_R r,cl_R xi,cl_R xj)
{
  return DSlater_4S_5S(r,xj,xi);
}

cl_R Nuclear_1S(cl_R r,cl_R xi)
{
  cl_R S = ZERO;
    S = 1LL/r - (1LL + r*xi)/(exp(2LL*r*xi)*r)

    ;

  return S;
}

cl_R Nuclear_2S(cl_R r,cl_R xi)
{
  cl_R S = ZERO;
    S = 1LL/r - (6LL + 9LL*r*xi + 6LL*Power(r,2LL)*Power(xi,2LL) + 2LL*Power(r,3LL)*Power(xi,3LL))/

     (6LL*exp(2LL*r*xi)*r)

    ;

  return S;
}

cl_R Nuclear_3S(cl_R r,cl_R xi)
{
  cl_R S = ZERO;
    S = 1LL/r - (45LL + 75LL*r*xi + 60LL*Power(r,2LL)*Power(xi,2LL) + 

       30LL*Power(r,3LL)*Power(xi,3LL) + 10LL*Power(r,4LL)*Power(xi,4LL) + 

       2LL*Power(r,5LL)*Power(xi,5LL))/(45LL*exp(2LL*r*xi)*r)

    ;

  return S;
}

cl_R Nuclear_4S(cl_R r,cl_R xi)
{
  cl_R S = ZERO;
    S = 1LL/r - (1260LL + 2205LL*r*xi + 1890LL*Power(r,2LL)*Power(xi,2LL) + 

       1050LL*Power(r,3LL)*Power(xi,3LL) + 420LL*Power(r,4LL)*Power(xi,4LL) + 

       126LL*Power(r,5LL)*Power(xi,5LL) + 28LL*Power(r,6LL)*Power(xi,6LL) + 

       4LL*Power(r,7LL)*Power(xi,7LL))/(1260LL*exp(2LL*r*xi)*r)

    ;

  return S;
}

cl_R Nuclear_5S(cl_R r,cl_R xi)
{
  cl_R S = ZERO;
    S = 1LL/r - (14175LL + 25515LL*r*xi + 22680LL*Power(r,2LL)*Power(xi,2LL) + 

       13230LL*Power(r,3LL)*Power(xi,3LL) + 5670LL*Power(r,4LL)*Power(xi,4LL) + 

       1890LL*Power(r,5LL)*Power(xi,5LL) + 504LL*Power(r,6LL)*Power(xi,6LL) + 

       108LL*Power(r,7LL)*Power(xi,7LL) + 18LL*Power(r,8LL)*Power(xi,8LL) + 

       2LL*Power(r,9LL)*Power(xi,9LL))/(14175LL*exp(2LL*r*xi)*r)

    ;

  return S;
}

cl_R DNuclear_1S(cl_R r,cl_R xi)
{
  cl_R S = ZERO;
    S = Power(r,-2LL) - (1LL + 2LL*r*xi + 2LL*Power(r,2LL)*Power(xi,2LL))/

     (exp(2LL*r*xi)*Power(r,2LL))

    ;

  return S;
}

cl_R DNuclear_2S(cl_R r,cl_R xi)
{
  cl_R S = ZERO;
    S = Power(r,-2LL) - (3LL + 6LL*r*xi + 6LL*Power(r,2LL)*Power(xi,2LL) + 

       4LL*Power(r,3LL)*Power(xi,3LL) + 2LL*Power(r,4LL)*Power(xi,4LL))/

     (3LL*exp(2LL*r*xi)*Power(r,2LL))

    ;

  return S;
}

cl_R DNuclear_3S(cl_R r,cl_R xi)
{
  cl_R S = ZERO;
    S = Power(r,-2LL) - (45LL + 90LL*r*xi + 90LL*Power(r,2LL)*Power(xi,2LL) + 

       60LL*Power(r,3LL)*Power(xi,3LL) + 30LL*Power(r,4LL)*Power(xi,4LL) + 

       12LL*Power(r,5LL)*Power(xi,5LL) + 4LL*Power(r,6LL)*Power(xi,6LL))/

     (45LL*exp(2LL*r*xi)*Power(r,2LL))

    ;

  return S;
}

cl_R DNuclear_4S(cl_R r,cl_R xi)
{
  cl_R S = ZERO;
    S = Power(r,-2LL) - (315LL + 630LL*r*xi + 630LL*Power(r,2LL)*Power(xi,2LL) + 

       420LL*Power(r,3LL)*Power(xi,3LL) + 210LL*Power(r,4LL)*Power(xi,4LL) + 

       84LL*Power(r,5LL)*Power(xi,5LL) + 28LL*Power(r,6LL)*Power(xi,6LL) + 

       8LL*Power(r,7LL)*Power(xi,7LL) + 2LL*Power(r,8LL)*Power(xi,8LL))/

     (315LL*exp(2LL*r*xi)*Power(r,2LL))

    ;

  return S;
}

cl_R DNuclear_5S(cl_R r,cl_R xi)
{
  cl_R S = ZERO;
    S = Power(r,-2LL) - (14175LL + 28350LL*r*xi + 28350LL*Power(r,2LL)*Power(xi,2LL) + 

       18900LL*Power(r,3LL)*Power(xi,3LL) + 9450LL*Power(r,4LL)*Power(xi,4LL) + 

       3780LL*Power(r,5LL)*Power(xi,5LL) + 1260LL*Power(r,6LL)*Power(xi,6LL) + 

       360LL*Power(r,7LL)*Power(xi,7LL) + 90LL*Power(r,8LL)*Power(xi,8LL) + 

       20LL*Power(r,9LL)*Power(xi,9LL) + 4LL*Power(r,10LL)*Power(xi,10LL))/

     (14175LL*exp(2LL*r*xi)*Power(r,2LL))

    ;

  return S;
}

typedef cl_R t_slater_SS_func(cl_R r,cl_R xi,cl_R xj);
typedef cl_R t_slater_NS_func(cl_R r,cl_R xi);
t_slater_SS_func (*Slater_SS[SLATER_MAX][SLATER_MAX]) = {
  {  Slater_1S_1S,  Slater_1S_2S,  Slater_1S_3S,  Slater_1S_4S,  Slater_1S_5S},
  {  Slater_2S_1S,  Slater_2S_2S,  Slater_2S_3S,  Slater_2S_4S,  Slater_2S_5S},
  {  Slater_3S_1S,  Slater_3S_2S,  Slater_3S_3S,  Slater_3S_4S,  Slater_3S_5S},
  {  Slater_4S_1S,  Slater_4S_2S,  Slater_4S_3S,  Slater_4S_4S,  Slater_4S_5S},
  {  Slater_5S_1S,  Slater_5S_2S,  Slater_5S_3S,  Slater_5S_4S,  Slater_5S_5S}
};

t_slater_SS_func (*DSlater_SS[SLATER_MAX][SLATER_MAX]) = {
  {  DSlater_1S_1S,  DSlater_1S_2S,  DSlater_1S_3S,  DSlater_1S_4S,  DSlater_1S_5S},
  {  DSlater_2S_1S,  DSlater_2S_2S,  DSlater_2S_3S,  DSlater_2S_4S,  DSlater_2S_5S},
  {  DSlater_3S_1S,  DSlater_3S_2S,  DSlater_3S_3S,  DSlater_3S_4S,  DSlater_3S_5S},
  {  DSlater_4S_1S,  DSlater_4S_2S,  DSlater_4S_3S,  DSlater_4S_4S,  DSlater_4S_5S},
  {  DSlater_5S_1S,  DSlater_5S_2S,  DSlater_5S_3S,  DSlater_5S_4S,  DSlater_5S_5S}
};

t_slater_NS_func (*Slater_NS[SLATER_MAX]) = {
  Nuclear_1S,  Nuclear_2S,  Nuclear_3S,  Nuclear_4S,  Nuclear_5S
};

t_slater_NS_func (*DSlater_NS[SLATER_MAX]) = {
  DNuclear_1S,  DNuclear_2S,  DNuclear_3S,  DNuclear_4S,  DNuclear_5S
};

static char *my_ftoa(double d)
{
  static char buf[256];
  sprintf(buf,"%g",d);
  if (strchr(buf,'.') == NULL) strcat(buf,".0");
  strcat(buf,"_80");
  return buf;
}

#endif
/* HAVE_CLN_CLN_H */

extern "C" double Coulomb_SS(double r,int i,int j,double xi,double xj)
{
  char buf[256];
  double S;
#ifdef HAVE_CLN_CLN_H
  cl_R cr,cxi,cxj,cS;

  if ((i < 1) || (i > SLATER_MAX) || (j < 1) || (j > SLATER_MAX)) {
    cerr << "Slater-Slater integral " << i << "  " << j << " not supported." << endl;
    exit(1);  }
  cxi = my_ftoa(xi);
  cxj = my_ftoa(xj);
  cr = my_ftoa(r);
  cS = Slater_SS[i-1][j-1](cr,cxi,cxj);
  return double_approx(cS);
#else
  cerr << "Can not compute Slater integrals without the CLN library" << endl;
  return 0.0;
#endif
}

extern "C" double Nuclear_SS(double r,int i,double xi)
{
  char buf[256];
  double S;
#ifdef HAVE_CLN_CLN_H
  cl_R cr,cxi,cxj,cS;

  if ((i < 1) || (i > SLATER_MAX)) {
    cerr << "Slater-Nuclear integral " << i << " not supported." << endl;
    exit(1);
  }
  cxi = my_ftoa(xi);
  cr = my_ftoa(r);
  cS = Slater_NS[i-1](cr,cxi);
  return double_approx(cS);
#else
  cerr << "Can not compute Slater integrals without the CLN library" << endl;
  return 0.0;
#endif
}

extern "C" double DCoulomb_SS(double r,int i,int j,double xi,double xj)
{
  char buf[256];
  double S;
#ifdef HAVE_CLN_CLN_H
  cl_R cr,cxi,cxj,cS;

  if ((i < 1) || (i > SLATER_MAX) || (j < 1) || (j > SLATER_MAX)) {
    cerr << "Slater-Slater integral " << i << "  " << j << " not supported." << endl;
    exit(1);  }
  cxi = my_ftoa(xi);
  cxj = my_ftoa(xj);
  cr = my_ftoa(r);
  cS = DSlater_SS[i-1][j-1](cr,cxi,cxj);
  return double_approx(cS);
#else
  cerr << "Can not compute Slater integrals without the CLN library" << endl;
  return 0.0;
#endif
}

extern "C" double DNuclear_SS(double r,int i,double xi)
{
  char buf[256];
  double S;
#ifdef HAVE_CLN_CLN_H
  cl_R cr,cxi,cxj,cS;

  if ((i < 1) || (i > SLATER_MAX)) {
    cerr << "Slater-Nuclear integral " << i << " not supported." << endl;
    exit(1);
  }
  cxi = my_ftoa(xi);
  cr = my_ftoa(r);
  cS = DSlater_NS[i-1](cr,cxi);
  return double_approx(cS);
#else
  cerr << "Can not compute Slater integrals without the CLN library" << endl;
  return 0.0;
#endif
}

