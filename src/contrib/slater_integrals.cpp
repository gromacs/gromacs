// slater_integrals.cpp (c) 2008 Paul J. van Maaren and David van der Spoel
#include <iostream>
using namespace std;
#include <config.h>
#ifdef HAVE_CLN_CLN_H
#include <cln/cln.h>
#include "slater_integrals.h"
using namespace cln;

#define PRECISION 80

static cl_F     ZERO = "0.0_80";
static cl_F      ONE = "1.0_80";
static cl_F      TWO = "2.0_80";
static cl_F    THREE = "3.0_80";
static cl_F     FOUR = "4.0_80";
static cl_F     FIVE = "5.0_80";
static cl_F      SIX = "6.0_80";
static cl_F    SEVEN = "7.0_80";
static cl_F    EIGHT = "8.0_80";
static cl_F     NINE = "9.0_80";
static float_format_t precision = float_format(80);

#define Power(x,y) (exp(y*ln(x)))

cl_F Slater_1S_1S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = (FIVE*xi)/cl_float(8.0,precision)

    ; } else {  S = (ONE/r)*(-cl_float(24.0,precision) + cl_float(24.0,precision)*exp(TWO*r*xi) - cl_float(33.0,precision)*r*xi - cl_float(18.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

        FOUR*Power(r,THREE)*Power(xi,THREE))/(cl_float(24.0,precision)*exp(TWO*r*xi))

    ; }
 
  }
  else {
      if (r == ZERO) { S = (xi*xj*(Power(xi,TWO) + THREE*xi*xj + Power(xj,TWO)))/Power(xi + xj,THREE)

    ; } else { S = (ONE/r)*(exp(TWO*(rxi + rxj))*Power(Power(rxi,TWO) - Power(rxj,TWO),THREE) + 

        exp(TWO*rxj)*Power(rxj,FOUR)*

         (-THREE*Power(rxi,TWO) - Power(rxi,THREE) + Power(rxj,TWO) + rxi*Power(rxj,TWO)) - 

        exp(TWO*rxi)*Power(rxi,FOUR)*

         (Power(rxi,TWO)*(ONE + rxj) - Power(rxj,TWO)*(THREE + rxj)))/

      (exp(TWO*(rxi + rxj))*Power(rxi - rxj,THREE)*Power(rxi + rxj,THREE))

     ; }
   
  }
  return S;
}

cl_F Slater_1S_2S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = (SEVEN*xi)/cl_float(16.0,precision)

    ; } else {  S = (ONE/r)*(-cl_float(240.0,precision) + cl_float(240.0,precision)*exp(TWO*r*xi) - cl_float(375.0,precision)*r*xi - cl_float(270.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

        cl_float(115.0,precision)*Power(r,THREE)*Power(xi,THREE) - cl_float(30.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - 

        FOUR*Power(r,FIVE)*Power(xi,FIVE))/(cl_float(240.0,precision)*exp(TWO*r*xi))

    ; }
 
  }
  else {
      if (r == ZERO) { S = (xi*xj*(Power(xi,FOUR) + FIVE*Power(xi,THREE)*xj + cl_float(10.0,precision)*Power(xi,TWO)*Power(xj,TWO) + 

          cl_float(10.0,precision)*xi*Power(xj,THREE) + TWO*Power(xj,FOUR)))/(cl_float(2.0,precision)*Power(xi + xj,FIVE))

    ; } else { S = (ONE/r)*(SIX*exp(TWO*(rxi + rxj))*Power(Power(rxi,TWO) - Power(rxj,TWO),FIVE) + 

        SIX*exp(TWO*rxj)*Power(rxj,SIX)*

         (-FOUR*Power(rxi,FOUR) - Power(rxi,FIVE) - FIVE*Power(rxi,TWO)*Power(rxj,TWO) + 

           Power(rxj,FOUR) + rxi*Power(rxj,FOUR)) - 

        exp(TWO*rxi)*Power(rxi,FOUR)*

         (Power(rxi,SIX)*(SIX + NINE*rxj + SIX*Power(rxj,TWO) + TWO*Power(rxj,THREE)) - 

           THREE*Power(rxi,FOUR)*Power(rxj,TWO)*

            (cl_float(10.0,precision) + cl_float(15.0,precision)*rxj + cl_float(10.0,precision)*Power(rxj,TWO) + TWO*Power(rxj,THREE)) + 

           THREE*Power(rxi,TWO)*Power(rxj,FOUR)*

            (cl_float(20.0,precision) + cl_float(33.0,precision)*rxj + cl_float(14.0,precision)*Power(rxj,TWO) + TWO*Power(rxj,THREE)) - 

           Power(rxj,SIX)*(cl_float(84.0,precision) + cl_float(63.0,precision)*rxj + cl_float(18.0,precision)*Power(rxj,TWO) + TWO*Power(rxj,THREE))))/

      (cl_float(6.0,precision)*exp(TWO*(rxi + rxj))*Power(rxi - rxj,FIVE)*Power(rxi + rxj,FIVE))

     ; }
   
  }
  return S;
}

cl_F Slater_1S_3S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = (cl_float(41.0,precision)*xi)/cl_float(128.0,precision)

    ; } else {  S = (ONE/r)*(-cl_float(120960.0,precision) + cl_float(120960.0,precision)*exp(TWO*r*xi) - cl_float(203175.0,precision)*r*xi - 

        cl_float(164430.0,precision)*Power(r,TWO)*Power(xi,TWO) - cl_float(84420.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

        cl_float(30240.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - cl_float(7728.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

        cl_float(1344.0,precision)*Power(r,SIX)*Power(xi,SIX) - cl_float(128.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN))/

      (cl_float(120960.0,precision)*exp(TWO*r*xi))

    ; }
 
  }
  else {
      if (r == ZERO) { S = (xi*xj*(Power(xi,SIX) + SEVEN*Power(xi,FIVE)*xj + cl_float(21.0,precision)*Power(xi,FOUR)*Power(xj,TWO) + 

          cl_float(35.0,precision)*Power(xi,THREE)*Power(xj,THREE) + cl_float(35.0,precision)*Power(xi,TWO)*Power(xj,FOUR) + 

          cl_float(21.0,precision)*xi*Power(xj,FIVE) + THREE*Power(xj,SIX)))/(cl_float(3.0,precision)*Power(xi + xj,SEVEN))

    ; } else { S = (ONE/r)*(cl_float(45.0,precision)*exp(TWO*(rxi + rxj))*Power(Power(rxi,TWO) - Power(rxj,TWO),SEVEN) + 

        cl_float(15.0,precision)*exp(TWO*rxj)*Power(rxj,EIGHT)*

         (-cl_float(15.0,precision)*Power(rxi,SIX) - THREE*Power(rxi,SEVEN) - cl_float(63.0,precision)*Power(rxi,FOUR)*Power(rxj,TWO) - 

           SEVEN*Power(rxi,FIVE)*Power(rxj,TWO) - cl_float(21.0,precision)*Power(rxi,TWO)*Power(rxj,FOUR) + 

           SEVEN*Power(rxi,THREE)*Power(rxj,FOUR) + THREE*Power(rxj,SIX) + THREE*rxi*Power(rxj,SIX)) + 

        exp(TWO*rxi)*Power(rxi,FOUR)*

         (-cl_float(10.0,precision)*Power(rxi,TWO)*Power(rxj,EIGHT)*

            (cl_float(135.0,precision) + cl_float(333.0,precision)*rxj + cl_float(228.0,precision)*Power(rxj,TWO) + cl_float(75.0,precision)*Power(rxj,THREE) + 

              cl_float(13.0,precision)*Power(rxj,FOUR) + Power(rxj,FIVE)) + 

           TWO*Power(rxj,cl_float(10.0,precision))*(cl_float(945.0,precision) + cl_float(945.0,precision)*rxj + cl_float(420.0,precision)*Power(rxj,TWO) + 

              cl_float(105.0,precision)*Power(rxj,THREE) + cl_float(15.0,precision)*Power(rxj,FOUR) + Power(rxj,FIVE)) - 

           Power(rxi,cl_float(10.0,precision))*(cl_float(45.0,precision) + cl_float(75.0,precision)*rxj + cl_float(60.0,precision)*Power(rxj,TWO) + cl_float(30.0,precision)*Power(rxj,THREE) + 

              cl_float(10.0,precision)*Power(rxj,FOUR) + TWO*Power(rxj,FIVE)) + 

           FIVE*Power(rxi,EIGHT)*Power(rxj,TWO)*

            (cl_float(63.0,precision) + cl_float(105.0,precision)*rxj + cl_float(84.0,precision)*Power(rxj,TWO) + cl_float(42.0,precision)*Power(rxj,THREE) + 

              cl_float(14.0,precision)*Power(rxj,FOUR) + TWO*Power(rxj,FIVE)) - 

           FIVE*Power(rxi,SIX)*Power(rxj,FOUR)*

            (cl_float(189.0,precision) + cl_float(315.0,precision)*rxj + cl_float(252.0,precision)*Power(rxj,TWO) + cl_float(132.0,precision)*Power(rxj,THREE) + 

              cl_float(36.0,precision)*Power(rxj,FOUR) + FOUR*Power(rxj,FIVE)) + 

           FIVE*Power(rxi,FOUR)*Power(rxj,SIX)*

            (cl_float(315.0,precision) + cl_float(513.0,precision)*rxj + cl_float(468.0,precision)*Power(rxj,TWO) + cl_float(204.0,precision)*Power(rxj,THREE) + 

              cl_float(44.0,precision)*Power(rxj,FOUR) + FOUR*Power(rxj,FIVE))))/

      (cl_float(45.0,precision)*exp(TWO*(rxi + rxj))*Power(rxi - rxj,SEVEN)*Power(rxi + rxj,SEVEN))

     ; }
   
  }
  return S;
}

cl_F Slater_1S_4S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = (cl_float(253.0,precision)*xi)/cl_float(1024.0,precision)

    ; } else {  S = (ONE/r)*(-cl_float(2903040.0,precision) + cl_float(2903040.0,precision)*exp(TWO*r*xi) - cl_float(5088825.0,precision)*r*xi - 

        cl_float(4371570.0,precision)*Power(r,TWO)*Power(xi,TWO) - cl_float(2439990.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

        cl_float(986580.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - cl_float(303912.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

        cl_float(72576.0,precision)*Power(r,SIX)*Power(xi,SIX) - cl_float(13248.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

        cl_float(1728.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - cl_float(128.0,precision)*Power(r,NINE)*Power(xi,NINE))/

      (cl_float(2.90304e6,precision)*exp(TWO*r*xi))

    ; }
 
  }
  else {
      if (r == ZERO) { S = (xi*xj*(Power(xi,EIGHT) + NINE*Power(xi,SEVEN)*xj + cl_float(36.0,precision)*Power(xi,SIX)*Power(xj,TWO) + 

          cl_float(84.0,precision)*Power(xi,FIVE)*Power(xj,THREE) + cl_float(126.0,precision)*Power(xi,FOUR)*Power(xj,FOUR) + 

          cl_float(126.0,precision)*Power(xi,THREE)*Power(xj,FIVE) + cl_float(84.0,precision)*Power(xi,TWO)*Power(xj,SIX) + 

          cl_float(36.0,precision)*xi*Power(xj,SEVEN) + FOUR*Power(xj,EIGHT)))/(cl_float(4.0,precision)*Power(xi + xj,NINE))

    ; } else { S = (ONE/r)*(cl_float(1260.0,precision)*exp(TWO*(rxi + rxj))*Power(Power(rxi,TWO) - Power(rxj,TWO),NINE) + 

        cl_float(1260.0,precision)*exp(TWO*rxj)*Power(rxj,cl_float(10.0,precision))*

         (-SIX*Power(rxi,EIGHT) - Power(rxi,NINE) - cl_float(51.0,precision)*Power(rxi,SIX)*Power(rxj,TWO) - 

           SIX*Power(rxi,SEVEN)*Power(rxj,TWO) - cl_float(63.0,precision)*Power(rxi,FOUR)*Power(rxj,FOUR) - 

           NINE*Power(rxi,TWO)*Power(rxj,SIX) + SIX*Power(rxi,THREE)*Power(rxj,SIX) + 

           Power(rxj,EIGHT) + rxi*Power(rxj,EIGHT)) - 

        exp(TWO*rxi)*Power(rxi,FOUR)*

         (cl_float(42.0,precision)*Power(rxi,cl_float(10.0,precision))*Power(rxj,FOUR)*

            (cl_float(1080.0,precision) + cl_float(1890.0,precision)*rxj + cl_float(1620.0,precision)*Power(rxj,TWO) + cl_float(900.0,precision)*Power(rxj,THREE) + 

              cl_float(360.0,precision)*Power(rxj,FOUR) + cl_float(111.0,precision)*Power(rxj,FIVE) + cl_float(22.0,precision)*Power(rxj,SIX) + 

              TWO*Power(rxj,SEVEN)) - cl_float(70.0,precision)*Power(rxi,EIGHT)*Power(rxj,SIX)*

            (cl_float(1512.0,precision) + cl_float(2646.0,precision)*rxj + cl_float(2268.0,precision)*Power(rxj,TWO) + cl_float(1248.0,precision)*Power(rxj,THREE) + 

              cl_float(528.0,precision)*Power(rxj,FOUR) + cl_float(153.0,precision)*Power(rxj,FIVE) + cl_float(26.0,precision)*Power(rxj,SIX) + 

              TWO*Power(rxj,SEVEN)) + cl_float(14.0,precision)*Power(rxi,TWO)*Power(rxj,cl_float(12.0,precision))*

            (cl_float(2970.0,precision) + cl_float(16335.0,precision)*rxj + cl_float(15390.0,precision)*Power(rxj,TWO) + cl_float(7110.0,precision)*Power(rxj,THREE) + 

              cl_float(1980.0,precision)*Power(rxj,FOUR) + cl_float(351.0,precision)*Power(rxj,FIVE) + cl_float(38.0,precision)*Power(rxj,SIX) + 

              TWO*Power(rxj,SEVEN)) - TWO*Power(rxj,cl_float(14.0,precision))*

            (cl_float(62370.0,precision) + cl_float(72765.0,precision)*rxj + cl_float(39690.0,precision)*Power(rxj,TWO) + cl_float(13230.0,precision)*Power(rxj,THREE) + 

              cl_float(2940.0,precision)*Power(rxj,FOUR) + cl_float(441.0,precision)*Power(rxj,FIVE) + cl_float(42.0,precision)*Power(rxj,SIX) + 

              TWO*Power(rxj,SEVEN)) + Power(rxi,cl_float(14.0,precision))*

            (cl_float(1260.0,precision) + cl_float(2205.0,precision)*rxj + cl_float(1890.0,precision)*Power(rxj,TWO) + cl_float(1050.0,precision)*Power(rxj,THREE) + 

              cl_float(420.0,precision)*Power(rxj,FOUR) + cl_float(126.0,precision)*Power(rxj,FIVE) + cl_float(28.0,precision)*Power(rxj,SIX) + 

              FOUR*Power(rxj,SEVEN)) - SEVEN*Power(rxi,cl_float(12.0,precision))*Power(rxj,TWO)*

            (cl_float(1620.0,precision) + cl_float(2835.0,precision)*rxj + cl_float(2430.0,precision)*Power(rxj,TWO) + cl_float(1350.0,precision)*Power(rxj,THREE) + 

              cl_float(540.0,precision)*Power(rxj,FOUR) + cl_float(162.0,precision)*Power(rxj,FIVE) + cl_float(36.0,precision)*Power(rxj,SIX) + 

              FOUR*Power(rxj,SEVEN)) + cl_float(35.0,precision)*Power(rxi,SIX)*Power(rxj,EIGHT)*

            (cl_float(4536.0,precision) + cl_float(7983.0,precision)*rxj + cl_float(6534.0,precision)*Power(rxj,TWO) + cl_float(4014.0,precision)*Power(rxj,THREE) + 

              cl_float(1644.0,precision)*Power(rxj,FOUR) + cl_float(414.0,precision)*Power(rxj,FIVE) + cl_float(60.0,precision)*Power(rxj,SIX) + 

              FOUR*Power(rxj,SEVEN)) - cl_float(21.0,precision)*Power(rxi,FOUR)*Power(rxj,cl_float(10.0,precision))*

            (cl_float(7920.0,precision) + cl_float(11385.0,precision)*rxj + cl_float(12330.0,precision)*Power(rxj,TWO) + cl_float(7410.0,precision)*Power(rxj,THREE) + 

              cl_float(2580.0,precision)*Power(rxj,FOUR) + cl_float(546.0,precision)*Power(rxj,FIVE) + cl_float(68.0,precision)*Power(rxj,SIX) + 

              FOUR*Power(rxj,SEVEN))))/

      (cl_float(1260.0,precision)*exp(TWO*(rxi + rxj))*Power(rxi - rxj,NINE)*Power(rxi + rxj,NINE))

     ; }
   
  }
  return S;
}

cl_F Slater_1S_5S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = (cl_float(2041.0,precision)*xi)/cl_float(10240.0,precision)

    ; } else {  S = (ONE/r)*(-cl_float(1596672000.0,precision) + cl_float(1596672000.0,precision)*exp(TWO*r*xi) - cl_float(2875101075.0,precision)*r*xi - 

        cl_float(2556858150.0,precision)*Power(r,TWO)*Power(xi,TWO) - cl_float(1492929900.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

        cl_float(641163600.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - cl_float(214719120.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

        cl_float(57879360.0,precision)*Power(r,SIX)*Power(xi,SIX) - cl_float(12735360.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

        cl_float(2280960.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - cl_float(323840.0,precision)*Power(r,NINE)*Power(xi,NINE) - 

        cl_float(33792.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - cl_float(2048.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)))/

      (cl_float(1.596672e9,precision)*exp(TWO*r*xi))

    ; }
 
  }
  else {
      if (r == ZERO) { S = (xi*xj*(Power(xi,cl_float(10.0,precision)) + cl_float(11.0,precision)*Power(xi,NINE)*xj + cl_float(55.0,precision)*Power(xi,EIGHT)*Power(xj,TWO) + 

          cl_float(165.0,precision)*Power(xi,SEVEN)*Power(xj,THREE) + cl_float(330.0,precision)*Power(xi,SIX)*Power(xj,FOUR) + 

          cl_float(462.0,precision)*Power(xi,FIVE)*Power(xj,FIVE) + cl_float(462.0,precision)*Power(xi,FOUR)*Power(xj,SIX) + 

          cl_float(330.0,precision)*Power(xi,THREE)*Power(xj,SEVEN) + cl_float(165.0,precision)*Power(xi,TWO)*Power(xj,EIGHT) + 

          cl_float(55.0,precision)*xi*Power(xj,NINE) + FIVE*Power(xj,cl_float(10.0,precision))))/(cl_float(5.0,precision)*Power(xi + xj,cl_float(11.0,precision)))

    ; } else { S = (ONE/r)*(cl_float(14175.0,precision)*exp(TWO*(rxi + rxj))*Power(Power(rxi,TWO) - Power(rxj,TWO),cl_float(11.0,precision)) + 

        cl_float(2835.0,precision)*exp(TWO*rxj)*Power(rxj,cl_float(12.0,precision))*

         (-cl_float(35.0,precision)*Power(rxi,cl_float(10.0,precision)) - FIVE*Power(rxi,cl_float(11.0,precision)) - cl_float(495.0,precision)*Power(rxi,EIGHT)*Power(rxj,TWO) - 

           cl_float(55.0,precision)*Power(rxi,NINE)*Power(rxj,TWO) - cl_float(1254.0,precision)*Power(rxi,SIX)*Power(rxj,FOUR) - 

           cl_float(66.0,precision)*Power(rxi,SEVEN)*Power(rxj,FOUR) - cl_float(726.0,precision)*Power(rxi,FOUR)*Power(rxj,SIX) + 

           cl_float(66.0,precision)*Power(rxi,FIVE)*Power(rxj,SIX) - cl_float(55.0,precision)*Power(rxi,TWO)*Power(rxj,EIGHT) + 

           cl_float(55.0,precision)*Power(rxi,THREE)*Power(rxj,EIGHT) + FIVE*Power(rxj,cl_float(10.0,precision)) + FIVE*rxi*Power(rxj,cl_float(10.0,precision))) 

    - exp(TWO*rxi)*Power(rxi,FOUR)*(Power(rxi,cl_float(18.0,precision))*

            (cl_float(14175.0,precision) + cl_float(25515.0,precision)*rxj + cl_float(22680.0,precision)*Power(rxj,TWO) + cl_float(13230.0,precision)*Power(rxj,THREE) + 

              cl_float(5670.0,precision)*Power(rxj,FOUR) + cl_float(1890.0,precision)*Power(rxj,FIVE) + cl_float(504.0,precision)*Power(rxj,SIX) + 

              cl_float(108.0,precision)*Power(rxj,SEVEN) + cl_float(18.0,precision)*Power(rxj,EIGHT) + TWO*Power(rxj,NINE)) - 

           NINE*Power(rxi,cl_float(16.0,precision))*Power(rxj,TWO)*

            (cl_float(17325.0,precision) + cl_float(31185.0,precision)*rxj + cl_float(27720.0,precision)*Power(rxj,TWO) + cl_float(16170.0,precision)*Power(rxj,THREE) + 

              cl_float(6930.0,precision)*Power(rxj,FOUR) + cl_float(2310.0,precision)*Power(rxj,FIVE) + cl_float(616.0,precision)*Power(rxj,SIX) + 

              cl_float(132.0,precision)*Power(rxj,SEVEN) + cl_float(22.0,precision)*Power(rxj,EIGHT) + TWO*Power(rxj,NINE)) + 

           cl_float(126.0,precision)*Power(rxi,cl_float(10.0,precision))*Power(rxj,EIGHT)*

            (cl_float(37125.0,precision) + cl_float(66825.0,precision)*rxj + cl_float(59400.0,precision)*Power(rxj,TWO) + cl_float(34725.0,precision)*Power(rxj,THREE) + 

              cl_float(14625.0,precision)*Power(rxj,FOUR) + cl_float(5043.0,precision)*Power(rxj,FIVE) + cl_float(1396.0,precision)*Power(rxj,SIX) + 

              cl_float(276.0,precision)*Power(rxj,SEVEN) + cl_float(34.0,precision)*Power(rxj,EIGHT) + TWO*Power(rxj,NINE)) - 

           cl_float(126.0,precision)*Power(rxi,EIGHT)*Power(rxj,cl_float(10.0,precision))*

            (cl_float(51975.0,precision) + cl_float(93420.0,precision)*rxj + cl_float(84240.0,precision)*Power(rxj,TWO) + cl_float(46815.0,precision)*Power(rxj,THREE) + 

              cl_float(20835.0,precision)*Power(rxj,FOUR) + cl_float(7485.0,precision)*Power(rxj,FIVE) + cl_float(1964.0,precision)*Power(rxj,SIX) + 

              cl_float(348.0,precision)*Power(rxj,SEVEN) + cl_float(38.0,precision)*Power(rxj,EIGHT) + TWO*Power(rxj,NINE)) + 

           NINE*Power(rxi,TWO)*Power(rxj,cl_float(16.0,precision))*

            (-cl_float(135135.0,precision) + cl_float(405405.0,precision)*rxj + cl_float(582120.0,precision)*Power(rxj,TWO) + cl_float(346500.0,precision)*Power(rxj,THREE) + 

              cl_float(124740.0,precision)*Power(rxj,FOUR) + cl_float(30492.0,precision)*Power(rxj,FIVE) + cl_float(5264.0,precision)*Power(rxj,SIX) + 

              cl_float(636.0,precision)*Power(rxj,SEVEN) + cl_float(50.0,precision)*Power(rxj,EIGHT) + TWO*Power(rxj,NINE)) - 

           Power(rxj,cl_float(18.0,precision))*(cl_float(2837835.0,precision) + cl_float(3648645.0,precision)*rxj + cl_float(2245320.0,precision)*Power(rxj,TWO) + 

              cl_float(873180.0,precision)*Power(rxj,THREE) + cl_float(238140.0,precision)*Power(rxj,FOUR) + cl_float(47628.0,precision)*Power(rxj,FIVE) + 

              cl_float(7056.0,precision)*Power(rxj,SIX) + cl_float(756.0,precision)*Power(rxj,SEVEN) + cl_float(54.0,precision)*Power(rxj,EIGHT) + 

              TWO*Power(rxj,NINE)) + NINE*Power(rxi,cl_float(14.0,precision))*Power(rxj,FOUR)*

            (cl_float(86625.0,precision) + cl_float(155925.0,precision)*rxj + cl_float(138600.0,precision)*Power(rxj,TWO) + cl_float(80850.0,precision)*Power(rxj,THREE) + 

              cl_float(34650.0,precision)*Power(rxj,FOUR) + cl_float(11550.0,precision)*Power(rxj,FIVE) + cl_float(3080.0,precision)*Power(rxj,SIX) + 

              cl_float(672.0,precision)*Power(rxj,SEVEN) + cl_float(104.0,precision)*Power(rxj,EIGHT) + EIGHT*Power(rxj,NINE)) - 

           cl_float(21.0,precision)*Power(rxi,cl_float(12.0,precision))*Power(rxj,SIX)*

            (cl_float(111375.0,precision) + cl_float(200475.0,precision)*rxj + cl_float(178200.0,precision)*Power(rxj,TWO) + cl_float(103950.0,precision)*Power(rxj,THREE) + 

              cl_float(44550.0,precision)*Power(rxj,FOUR) + cl_float(14778.0,precision)*Power(rxj,FIVE) + cl_float(4056.0,precision)*Power(rxj,SIX) + 

              cl_float(864.0,precision)*Power(rxj,SEVEN) + cl_float(120.0,precision)*Power(rxj,EIGHT) + EIGHT*Power(rxj,NINE)) + 

           cl_float(21.0,precision)*Power(rxi,SIX)*Power(rxj,cl_float(12.0,precision))*

            (cl_float(307125.0,precision) + cl_float(594945.0,precision)*rxj + cl_float(456840.0,precision)*Power(rxj,TWO) + cl_float(281790.0,precision)*Power(rxj,THREE) + 

              cl_float(137430.0,precision)*Power(rxj,FOUR) + cl_float(47250.0,precision)*Power(rxj,FIVE) + cl_float(11064.0,precision)*Power(rxj,SIX) + 

              cl_float(1728.0,precision)*Power(rxj,SEVEN) + cl_float(168.0,precision)*Power(rxj,EIGHT) + EIGHT*Power(rxj,NINE)) - 

           NINE*Power(rxi,FOUR)*Power(rxj,cl_float(14.0,precision))*

            (cl_float(675675.0,precision) + cl_float(675675.0,precision)*rxj + cl_float(748440.0,precision)*Power(rxj,TWO) + cl_float(561330.0,precision)*Power(rxj,THREE) + 

              cl_float(256410.0,precision)*Power(rxj,FOUR) + cl_float(76230.0,precision)*Power(rxj,FIVE) + cl_float(15400.0,precision)*Power(rxj,SIX) + 

              cl_float(2112.0,precision)*Power(rxj,SEVEN) + cl_float(184.0,precision)*Power(rxj,EIGHT) + EIGHT*Power(rxj,NINE))))/

      (cl_float(14175.0,precision)*exp(TWO*(rxi + rxj))*Power(rxi - rxj,cl_float(11.0,precision))*Power(rxi + rxj,cl_float(11.0,precision)))

     ; }
   
  }
  return S;
}

cl_F Slater_2S_2S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = (cl_float(93.0,precision)*xi)/cl_float(256.0,precision)

    ; } else {  S = (ONE/r)*(-cl_float(80640.0,precision) + cl_float(80640.0,precision)*exp(TWO*r*xi) - cl_float(131985.0,precision)*r*xi - 

        cl_float(102690.0,precision)*Power(r,TWO)*Power(xi,TWO) - cl_float(49980.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

        cl_float(16800.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - cl_float(4032.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

        cl_float(672.0,precision)*Power(r,SIX)*Power(xi,SIX) - cl_float(64.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN))/

      (cl_float(80640.0,precision)*exp(TWO*r*xi))

    ; }
 
  }
  else {
      if (r == ZERO) { S = (xi*xj*(Power(xi,SIX) + SEVEN*Power(xi,FIVE)*xj + cl_float(21.0,precision)*Power(xi,FOUR)*Power(xj,TWO) + 

          cl_float(35.0,precision)*Power(xi,THREE)*Power(xj,THREE) + cl_float(21.0,precision)*Power(xi,TWO)*Power(xj,FOUR) + 

          SEVEN*xi*Power(xj,FIVE) + Power(xj,SIX)))/(cl_float(2.0,precision)*Power(xi + xj,SEVEN))

    ; } else { S = (ONE/r)*(SIX*exp(TWO*(rxi + rxj))*Power(Power(rxi,TWO) - Power(rxj,TWO),SEVEN) - 

        exp(TWO*rxi)*Power(rxi,SIX)*

         (cl_float(21.0,precision)*Power(rxi,FOUR)*Power(rxj,FOUR)*(SIX + cl_float(11.0,precision)*rxj + TWO*Power(rxj,TWO)) - 

           TWO*Power(rxj,EIGHT)*(cl_float(90.0,precision) + cl_float(54.0,precision)*rxj + cl_float(12.0,precision)*Power(rxj,TWO) + Power(rxj,THREE)) + 

           Power(rxi,EIGHT)*(SIX + NINE*rxj + SIX*Power(rxj,TWO) + TWO*Power(rxj,THREE)) + 

           Power(rxi,TWO)*Power(rxj,SIX)*

            (-cl_float(390.0,precision) - cl_float(69.0,precision)*rxj + cl_float(18.0,precision)*Power(rxj,TWO) + FOUR*Power(rxj,THREE)) - 

           Power(rxi,SIX)*Power(rxj,TWO)*

            (cl_float(42.0,precision) + cl_float(63.0,precision)*rxj + cl_float(42.0,precision)*Power(rxj,TWO) + FOUR*Power(rxj,THREE))) + 

        exp(TWO*rxj)*Power(rxj,SIX)*

         (-cl_float(24.0,precision)*Power(rxi,cl_float(10.0,precision)) - TWO*Power(rxi,cl_float(11.0,precision)) - cl_float(69.0,precision)*Power(rxi,SEVEN)*Power(rxj,TWO) + 

           SIX*Power(rxj,EIGHT) + NINE*rxi*Power(rxj,EIGHT) + 

           FOUR*Power(rxi,NINE)*(-cl_float(27.0,precision) + Power(rxj,TWO)) + 

           cl_float(18.0,precision)*Power(rxi,EIGHT)*(-cl_float(10.0,precision) + Power(rxj,TWO)) + 

           SIX*Power(rxi,TWO)*Power(rxj,SIX)*(-SEVEN + Power(rxj,TWO)) - 

           cl_float(42.0,precision)*Power(rxi,FOUR)*Power(rxj,FOUR)*(-THREE + Power(rxj,TWO)) + 

           Power(rxi,THREE)*Power(rxj,SIX)*(-cl_float(63.0,precision) + TWO*Power(rxj,TWO)) + 

           SIX*Power(rxi,SIX)*Power(rxj,TWO)*(-cl_float(65.0,precision) + SEVEN*Power(rxj,TWO)) + 

           Power(rxi,FIVE)*(cl_float(231.0,precision)*Power(rxj,FOUR) - FOUR*Power(rxj,SIX))))/

      (cl_float(6.0,precision)*exp(TWO*(rxi + rxj))*Power(rxi - rxj,SEVEN)*Power(rxi + rxj,SEVEN))

     ; }
   
  }
  return S;
}

cl_F Slater_2S_3S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = (cl_float(451.0,precision)*xi)/cl_float(1536.0,precision)

    ; } else {  S = (ONE/r)*(-cl_float(4354560.0,precision) + cl_float(4354560.0,precision)*exp(TWO*r*xi) - cl_float(7430535.0,precision)*r*xi - 

        cl_float(6151950.0,precision)*Power(r,TWO)*Power(xi,TWO) - cl_float(3275370.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

        cl_float(1251180.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - cl_float(361368.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

        cl_float(80640.0,precision)*Power(r,SIX)*Power(xi,SIX) - cl_float(13824.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

        cl_float(1728.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - cl_float(128.0,precision)*Power(r,NINE)*Power(xi,NINE))/

      (cl_float(4.35456e6,precision)*exp(TWO*r*xi))

    ; }
 
  }
  else {
      if (r == ZERO) { S = (xi*xj*(TWO*Power(xi,EIGHT) + cl_float(18.0,precision)*Power(xi,SEVEN)*xj + cl_float(72.0,precision)*Power(xi,SIX)*Power(xj,TWO) + 

          cl_float(168.0,precision)*Power(xi,FIVE)*Power(xj,THREE) + cl_float(252.0,precision)*Power(xi,FOUR)*Power(xj,FOUR) + 

          cl_float(252.0,precision)*Power(xi,THREE)*Power(xj,FIVE) + cl_float(108.0,precision)*Power(xi,TWO)*Power(xj,SIX) + 

          cl_float(27.0,precision)*xi*Power(xj,SEVEN) + THREE*Power(xj,EIGHT)))/(cl_float(6.0,precision)*Power(xi + xj,NINE))

    ; } else { S = (ONE/r)*(cl_float(90.0,precision)*exp(TWO*(rxi + rxj))*Power(Power(rxi,TWO) - Power(rxj,TWO),NINE) + 

        FIVE*exp(TWO*rxj)*Power(rxj,EIGHT)*

         (-cl_float(90.0,precision)*Power(rxi,cl_float(12.0,precision)) - SIX*Power(rxi,cl_float(13.0,precision)) + cl_float(18.0,precision)*Power(rxj,cl_float(10.0,precision)) + 

           cl_float(27.0,precision)*rxi*Power(rxj,cl_float(10.0,precision)) + cl_float(18.0,precision)*Power(rxi,TWO)*Power(rxj,EIGHT)*

            (-NINE + Power(rxj,TWO)) - cl_float(162.0,precision)*Power(rxi,FOUR)*Power(rxj,SIX)*

            (-FOUR + Power(rxj,TWO)) - cl_float(198.0,precision)*Power(rxi,cl_float(10.0,precision))*(FIVE + Power(rxj,TWO)) - 

           cl_float(108.0,precision)*Power(rxi,SIX)*Power(rxj,FOUR)*(cl_float(36.0,precision) + Power(rxj,TWO)) + 

           TWO*Power(rxi,FIVE)*Power(rxj,SIX)*(cl_float(675.0,precision) + Power(rxj,TWO)) - 

           cl_float(18.0,precision)*Power(rxi,SEVEN)*Power(rxj,FOUR)*(-cl_float(81.0,precision) + TWO*Power(rxj,TWO)) + 

           THREE*Power(rxi,THREE)*Power(rxj,EIGHT)*(-cl_float(81.0,precision) + TWO*Power(rxj,TWO)) - 

           Power(rxi,cl_float(11.0,precision))*(cl_float(495.0,precision) + TWO*Power(rxj,TWO)) + 

           NINE*Power(rxi,NINE)*Power(rxj,TWO)*(-cl_float(233.0,precision) + FOUR*Power(rxj,TWO)) + 

           SIX*Power(rxi,EIGHT)*Power(rxj,TWO)*(-cl_float(1063.0,precision) + cl_float(90.0,precision)*Power(rxj,TWO))) - 

        TWO*exp(TWO*rxi)*Power(rxi,SIX)*

         (-cl_float(90.0,precision)*Power(rxi,SIX)*Power(rxj,SIX)*

            (cl_float(42.0,precision) + cl_float(65.0,precision)*rxj + cl_float(76.0,precision)*Power(rxj,TWO) + cl_float(22.0,precision)*Power(rxj,THREE) + TWO*Power(rxj,FOUR)) - 

           TWO*Power(rxj,cl_float(12.0,precision))*(cl_float(2970.0,precision) + cl_float(2475.0,precision)*rxj + cl_float(900.0,precision)*Power(rxj,TWO) + 

              cl_float(180.0,precision)*Power(rxj,THREE) + cl_float(20.0,precision)*Power(rxj,FOUR) + Power(rxj,FIVE)) + 

           cl_float(10.0,precision)*Power(rxi,EIGHT)*Power(rxj,FOUR)*

            (cl_float(162.0,precision) + cl_float(270.0,precision)*rxj + cl_float(216.0,precision)*Power(rxj,TWO) + cl_float(122.0,precision)*Power(rxj,THREE) + 

              cl_float(22.0,precision)*Power(rxj,FOUR) + Power(rxj,FIVE)) - 

           FIVE*Power(rxi,FOUR)*Power(rxj,EIGHT)*

            (-cl_float(639.0,precision) - cl_float(3555.0,precision)*rxj - cl_float(1452.0,precision)*Power(rxj,TWO) - cl_float(174.0,precision)*Power(rxj,THREE) + 

              SIX*Power(rxj,FOUR) + TWO*Power(rxj,FIVE)) + 

           Power(rxi,cl_float(12.0,precision))*(cl_float(45.0,precision) + cl_float(75.0,precision)*rxj + cl_float(60.0,precision)*Power(rxj,TWO) + cl_float(30.0,precision)*Power(rxj,THREE) + 

              cl_float(10.0,precision)*Power(rxj,FOUR) + TWO*Power(rxj,FIVE)) - 

           Power(rxi,cl_float(10.0,precision))*Power(rxj,TWO)*

            (cl_float(405.0,precision) + cl_float(675.0,precision)*rxj + cl_float(540.0,precision)*Power(rxj,TWO) + cl_float(270.0,precision)*Power(rxj,THREE) + 

              cl_float(90.0,precision)*Power(rxj,FOUR) + EIGHT*Power(rxj,FIVE)) + 

           Power(rxi,TWO)*Power(rxj,cl_float(10.0,precision))*

            (-cl_float(21615.0,precision) - cl_float(9075.0,precision)*rxj - cl_float(300.0,precision)*Power(rxj,TWO) + cl_float(490.0,precision)*Power(rxj,THREE) + 

              cl_float(110.0,precision)*Power(rxj,FOUR) + EIGHT*Power(rxj,FIVE))))/

      (cl_float(90.0,precision)*exp(TWO*(rxi + rxj))*Power(rxi - rxj,NINE)*Power(rxi + rxj,NINE))

     ; }
   
  }
  return S;
}

cl_F Slater_2S_4S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = (cl_float(975.0,precision)*xi)/cl_float(4096.0,precision)

    ; } else {  S = (ONE/r)*(-cl_float(638668800.0,precision) + cl_float(638668800.0,precision)*exp(TWO*r*xi) - cl_float(1125310725.0,precision)*r*xi - 

        cl_float(973283850.0,precision)*Power(r,TWO)*Power(xi,TWO) - cl_float(549063900.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

        cl_float(226195200.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - cl_float(72099720.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

        cl_float(18350640.0,precision)*Power(r,SIX)*Power(xi,SIX) - cl_float(3785760.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

        cl_float(633600.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - cl_float(84480.0,precision)*Power(r,NINE)*Power(xi,NINE) - 

        cl_float(8448.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - cl_float(512.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)))/

      (cl_float(6.386688e8,precision)*exp(TWO*r*xi))

    ; }
 
  }
  else {
      if (r == ZERO) { S = (xi*xj*(Power(xi,cl_float(10.0,precision)) + cl_float(11.0,precision)*Power(xi,NINE)*xj + cl_float(55.0,precision)*Power(xi,EIGHT)*Power(xj,TWO) + 

          cl_float(165.0,precision)*Power(xi,SEVEN)*Power(xj,THREE) + cl_float(330.0,precision)*Power(xi,SIX)*Power(xj,FOUR) + 

          cl_float(462.0,precision)*Power(xi,FIVE)*Power(xj,FIVE) + cl_float(462.0,precision)*Power(xi,FOUR)*Power(xj,SIX) + 

          cl_float(330.0,precision)*Power(xi,THREE)*Power(xj,SEVEN) + cl_float(110.0,precision)*Power(xi,TWO)*Power(xj,EIGHT) + 

          cl_float(22.0,precision)*xi*Power(xj,NINE) + TWO*Power(xj,cl_float(10.0,precision))))/(cl_float(4.0,precision)*Power(xi + xj,cl_float(11.0,precision)))

    ; } else { S = (ONE/r)*(cl_float(1260.0,precision)*exp(TWO*(rxi + rxj))*Power(Power(rxi,TWO) - Power(rxj,TWO),cl_float(11.0,precision)) + 

        cl_float(210.0,precision)*exp(TWO*rxj)*Power(rxj,cl_float(10.0,precision))*

         (-cl_float(36.0,precision)*Power(rxi,cl_float(14.0,precision)) - TWO*Power(rxi,cl_float(15.0,precision)) - 

           cl_float(1287.0,precision)*Power(rxi,NINE)*Power(rxj,FOUR) + SIX*Power(rxj,cl_float(12.0,precision)) + 

           NINE*rxi*Power(rxj,cl_float(12.0,precision)) - cl_float(22.0,precision)*Power(rxi,SEVEN)*Power(rxj,SIX)*

            (-cl_float(135.0,precision) + Power(rxj,TWO)) + 

           SIX*Power(rxi,TWO)*Power(rxj,cl_float(10.0,precision))*(-cl_float(11.0,precision) + Power(rxj,TWO)) - 

           cl_float(66.0,precision)*Power(rxi,FOUR)*Power(rxj,EIGHT)*(-FIVE + Power(rxj,TWO)) + 

           EIGHT*Power(rxi,FIVE)*Power(rxj,EIGHT)*(cl_float(99.0,precision) + Power(rxj,TWO)) + 

           Power(rxi,THREE)*Power(rxj,cl_float(10.0,precision))*(-cl_float(99.0,precision) + TWO*Power(rxj,TWO)) - 

           cl_float(132.0,precision)*Power(rxi,SIX)*Power(rxj,SIX)*(cl_float(27.0,precision) + TWO*Power(rxj,TWO)) - 

           cl_float(78.0,precision)*Power(rxi,cl_float(12.0,precision))*(SEVEN + THREE*Power(rxj,TWO)) - 

           TWO*Power(rxi,cl_float(13.0,precision))*(cl_float(117.0,precision) + FOUR*Power(rxj,TWO)) + 

           cl_float(66.0,precision)*Power(rxi,EIGHT)*Power(rxj,FOUR)*(-cl_float(191.0,precision) + SIX*Power(rxj,TWO)) + 

           Power(rxi,cl_float(11.0,precision))*Power(rxj,TWO)*(-cl_float(2151.0,precision) + cl_float(22.0,precision)*Power(rxj,TWO)) + 

           SIX*Power(rxi,cl_float(10.0,precision))*Power(rxj,TWO)*(-cl_float(1099.0,precision) + cl_float(33.0,precision)*Power(rxj,TWO))) - 

        exp(TWO*rxi)*Power(rxi,SIX)*

         (cl_float(385.0,precision)*Power(rxi,EIGHT)*Power(rxj,EIGHT)*

            (cl_float(1080.0,precision) + cl_float(1935.0,precision)*rxj + cl_float(1350.0,precision)*Power(rxj,TWO) + cl_float(1170.0,precision)*Power(rxj,THREE) + 

              cl_float(420.0,precision)*Power(rxj,FOUR) + cl_float(66.0,precision)*Power(rxj,FIVE) + FOUR*Power(rxj,SIX)) - 

           FOUR*Power(rxj,cl_float(16.0,precision))*(cl_float(135135.0,precision) + cl_float(135135.0,precision)*rxj + cl_float(62370.0,precision)*Power(rxj,TWO) + 

              cl_float(17325.0,precision)*Power(rxj,THREE) + cl_float(3150.0,precision)*Power(rxj,FOUR) + cl_float(378.0,precision)*Power(rxj,FIVE) + 

              cl_float(28.0,precision)*Power(rxj,SIX) + Power(rxj,SEVEN)) + 

           Power(rxi,cl_float(16.0,precision))*(cl_float(1260.0,precision) + cl_float(2205.0,precision)*rxj + cl_float(1890.0,precision)*Power(rxj,TWO) + 

              cl_float(1050.0,precision)*Power(rxj,THREE) + cl_float(420.0,precision)*Power(rxj,FOUR) + cl_float(126.0,precision)*Power(rxj,FIVE) + 

              cl_float(28.0,precision)*Power(rxj,SIX) + FOUR*Power(rxj,SEVEN)) + 

           SEVEN*Power(rxi,SIX)*Power(rxj,cl_float(10.0,precision))*

            (-cl_float(99540.0,precision) - cl_float(58095.0,precision)*rxj - cl_float(190710.0,precision)*Power(rxj,TWO) - cl_float(100950.0,precision)*Power(rxj,THREE) - 

              cl_float(21660.0,precision)*Power(rxj,FOUR) - cl_float(1938.0,precision)*Power(rxj,FIVE) - FOUR*Power(rxj,SIX) + 

              EIGHT*Power(rxj,SEVEN)) - SEVEN*Power(rxi,FOUR)*Power(rxj,cl_float(12.0,precision))*

            (cl_float(114660.0,precision) - cl_float(343395.0,precision)*rxj - cl_float(242910.0,precision)*Power(rxj,TWO) - cl_float(61950.0,precision)*Power(rxj,THREE) - 

              cl_float(6060.0,precision)*Power(rxj,FOUR) + cl_float(282.0,precision)*Power(rxj,FIVE) + cl_float(116.0,precision)*Power(rxj,SIX) + 

              EIGHT*Power(rxj,SEVEN)) + SEVEN*Power(rxi,cl_float(12.0,precision))*Power(rxj,FOUR)*

            (cl_float(9900.0,precision) + cl_float(17325.0,precision)*rxj + cl_float(14850.0,precision)*Power(rxj,TWO) + cl_float(8250.0,precision)*Power(rxj,THREE) + 

              cl_float(3300.0,precision)*Power(rxj,FOUR) + cl_float(1074.0,precision)*Power(rxj,FIVE) + cl_float(164.0,precision)*Power(rxj,SIX) + 

              EIGHT*Power(rxj,SEVEN)) - SEVEN*Power(rxi,cl_float(10.0,precision))*Power(rxj,SIX)*

            (cl_float(29700.0,precision) + cl_float(51975.0,precision)*rxj + cl_float(44550.0,precision)*Power(rxj,TWO) + cl_float(23850.0,precision)*Power(rxj,THREE) + 

              cl_float(11700.0,precision)*Power(rxj,FOUR) + cl_float(2814.0,precision)*Power(rxj,FIVE) + cl_float(284.0,precision)*Power(rxj,SIX) + 

              EIGHT*Power(rxj,SEVEN)) - Power(rxi,cl_float(14.0,precision))*Power(rxj,TWO)*

            (cl_float(13860.0,precision) + cl_float(24255.0,precision)*rxj + cl_float(20790.0,precision)*Power(rxj,TWO) + cl_float(11550.0,precision)*Power(rxj,THREE) + 

              cl_float(4620.0,precision)*Power(rxj,FOUR) + cl_float(1386.0,precision)*Power(rxj,FIVE) + cl_float(308.0,precision)*Power(rxj,SIX) + 

              cl_float(24.0,precision)*Power(rxj,SEVEN)) + Power(rxi,TWO)*Power(rxj,cl_float(14.0,precision))*

            (-cl_float(3063060.0,precision) - cl_float(1936935.0,precision)*rxj - cl_float(408870.0,precision)*Power(rxj,TWO) + cl_float(11550.0,precision)*Power(rxj,THREE) + 

              cl_float(23100.0,precision)*Power(rxj,FOUR) + cl_float(5082.0,precision)*Power(rxj,FIVE) + cl_float(532.0,precision)*Power(rxj,SIX) + 

              cl_float(24.0,precision)*Power(rxj,SEVEN))))/

      (cl_float(1260.0,precision)*exp(TWO*(rxi + rxj))*Power(rxi - rxj,cl_float(11.0,precision))*Power(rxi + rxj,cl_float(11.0,precision)))

     ; }
   
  }
  return S;
}

cl_F Slater_2S_5S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = (cl_float(2011.0,precision)*xi)/cl_float(10240.0,precision)

    ; } else {  S = (ONE/r)*(-cl_float(124540416000.0,precision) + cl_float(124540416000.0,precision)*exp(TWO*r*xi) - cl_float(224622748350.0,precision)*r*xi - 

        cl_float(200164664700.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

        cl_float(117249207075.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

        cl_float(50639138550.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - 

        cl_float(17132415300.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - cl_float(4704860160.0,precision)*Power(r,SIX)*Power(xi,SIX) - 

        cl_float(1071195840.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - cl_float(204478560.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - 

        cl_float(32809920.0,precision)*Power(r,NINE)*Power(xi,NINE) - cl_float(4392960.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - 

        cl_float(479232.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)) - cl_float(39936.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(12.0,precision)) - 

        cl_float(2048.0,precision)*Power(r,cl_float(13.0,precision))*Power(xi,cl_float(13.0,precision)))/(cl_float(1.24540416e11,precision)*exp(TWO*r*xi))

    ; }
 
  }
  else {
      if (r == ZERO) { S = (xi*xj*(TWO*Power(xi,cl_float(12.0,precision)) + cl_float(26.0,precision)*Power(xi,cl_float(11.0,precision))*xj + cl_float(156.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,TWO) + 

          cl_float(572.0,precision)*Power(xi,NINE)*Power(xj,THREE) + cl_float(1430.0,precision)*Power(xi,EIGHT)*Power(xj,FOUR) + 

          cl_float(2574.0,precision)*Power(xi,SEVEN)*Power(xj,FIVE) + cl_float(3432.0,precision)*Power(xi,SIX)*Power(xj,SIX) + 

          cl_float(3432.0,precision)*Power(xi,FIVE)*Power(xj,SEVEN) + cl_float(2574.0,precision)*Power(xi,FOUR)*Power(xj,EIGHT) + 

          cl_float(1430.0,precision)*Power(xi,THREE)*Power(xj,NINE) + cl_float(390.0,precision)*Power(xi,TWO)*Power(xj,cl_float(10.0,precision)) + 

          cl_float(65.0,precision)*xi*Power(xj,cl_float(11.0,precision)) + FIVE*Power(xj,cl_float(12.0,precision))))/(cl_float(10.0,precision)*Power(xi + xj,cl_float(13.0,precision)))

    ; } else { S = (ONE/r)*(cl_float(28350.0,precision)*exp(TWO*(rxi + rxj))*Power(Power(rxi,TWO) - Power(rxj,TWO),cl_float(13.0,precision)) + 

        cl_float(945.0,precision)*exp(TWO*rxj)*Power(rxj,cl_float(12.0,precision))*

         (-cl_float(210.0,precision)*Power(rxi,cl_float(16.0,precision)) - cl_float(10.0,precision)*Power(rxi,cl_float(17.0,precision)) + cl_float(30.0,precision)*Power(rxj,cl_float(14.0,precision)) + 

           cl_float(45.0,precision)*rxi*Power(rxj,cl_float(14.0,precision)) + cl_float(39.0,precision)*Power(rxi,SEVEN)*Power(rxj,EIGHT)*

            (cl_float(1309.0,precision) - TWO*Power(rxj,TWO)) + 

           cl_float(858.0,precision)*Power(rxi,EIGHT)*Power(rxj,SIX)*(-cl_float(305.0,precision) + Power(rxj,TWO)) + 

           cl_float(30.0,precision)*Power(rxi,TWO)*Power(rxj,cl_float(12.0,precision))*(-cl_float(13.0,precision) + Power(rxj,TWO)) - 

           cl_float(390.0,precision)*Power(rxi,FOUR)*Power(rxj,cl_float(10.0,precision))*(-SIX + Power(rxj,TWO)) - 

           cl_float(143.0,precision)*Power(rxi,NINE)*Power(rxj,SIX)*(-cl_float(153.0,precision) + TWO*Power(rxj,TWO)) + 

           FIVE*Power(rxi,THREE)*Power(rxj,cl_float(12.0,precision))*(-cl_float(117.0,precision) + TWO*Power(rxj,TWO)) - 

           cl_float(45.0,precision)*Power(rxi,cl_float(15.0,precision))*(cl_float(35.0,precision) + TWO*Power(rxj,TWO)) - 

           cl_float(138.0,precision)*Power(rxi,cl_float(12.0,precision))*Power(rxj,TWO)*(cl_float(580.0,precision) + cl_float(13.0,precision)*Power(rxj,TWO)) - 

           cl_float(150.0,precision)*Power(rxi,cl_float(14.0,precision))*(cl_float(28.0,precision) + cl_float(17.0,precision)*Power(rxj,TWO)) + 

           cl_float(13.0,precision)*Power(rxi,cl_float(11.0,precision))*Power(rxj,FOUR)*(-cl_float(4071.0,precision) + cl_float(22.0,precision)*Power(rxj,TWO)) + 

           THREE*Power(rxi,cl_float(13.0,precision))*Power(rxj,TWO)*(-cl_float(8135.0,precision) + cl_float(26.0,precision)*Power(rxj,TWO)) + 

           THREE*Power(rxi,FIVE)*Power(rxj,cl_float(10.0,precision))*(cl_float(2171.0,precision) + cl_float(30.0,precision)*Power(rxj,TWO)) + 

           cl_float(234.0,precision)*Power(rxi,cl_float(10.0,precision))*Power(rxj,FOUR)*(-cl_float(1235.0,precision) + cl_float(33.0,precision)*Power(rxj,TWO)) - 

           cl_float(78.0,precision)*Power(rxi,SIX)*Power(rxj,EIGHT)*(cl_float(550.0,precision) + cl_float(47.0,precision)*Power(rxj,TWO))) - 

        TWO*exp(TWO*rxi)*Power(rxi,SIX)*

         (-cl_float(819.0,precision)*Power(rxi,cl_float(10.0,precision))*Power(rxj,cl_float(10.0,precision))*

            (cl_float(22275.0,precision) + cl_float(39780.0,precision)*rxj + cl_float(38160.0,precision)*Power(rxj,TWO) + cl_float(16560.0,precision)*Power(rxj,THREE) + 

              cl_float(9840.0,precision)*Power(rxj,FOUR) + cl_float(3900.0,precision)*Power(rxj,FIVE) + cl_float(816.0,precision)*Power(rxj,SIX) + 

              cl_float(88.0,precision)*Power(rxj,SEVEN) + FOUR*Power(rxj,EIGHT)) + 

           Power(rxi,cl_float(20.0,precision))*(cl_float(14175.0,precision) + cl_float(25515.0,precision)*rxj + cl_float(22680.0,precision)*Power(rxj,TWO) + 

              cl_float(13230.0,precision)*Power(rxj,THREE) + cl_float(5670.0,precision)*Power(rxj,FOUR) + cl_float(1890.0,precision)*Power(rxj,FIVE) + 

              cl_float(504.0,precision)*Power(rxj,SIX) + cl_float(108.0,precision)*Power(rxj,SEVEN) + cl_float(18.0,precision)*Power(rxj,EIGHT) + 

              TWO*Power(rxj,NINE)) - Power(rxj,cl_float(20.0,precision))*

            (cl_float(16216200.0,precision) + cl_float(18243225.0,precision)*rxj + cl_float(9729720.0,precision)*Power(rxj,TWO) + 

              cl_float(3243240.0,precision)*Power(rxj,THREE) + cl_float(748440.0,precision)*Power(rxj,FOUR) + 

              cl_float(124740.0,precision)*Power(rxj,FIVE) + cl_float(15120.0,precision)*Power(rxj,SIX) + cl_float(1296.0,precision)*Power(rxj,SEVEN) + 

              cl_float(72.0,precision)*Power(rxj,EIGHT) + TWO*Power(rxj,NINE)) + 

           cl_float(18.0,precision)*Power(rxi,cl_float(16.0,precision))*Power(rxj,FOUR)*

            (cl_float(61425.0,precision) + cl_float(110565.0,precision)*rxj + cl_float(98280.0,precision)*Power(rxj,TWO) + cl_float(57330.0,precision)*Power(rxj,THREE) + 

              cl_float(24570.0,precision)*Power(rxj,FOUR) + cl_float(8190.0,precision)*Power(rxj,FIVE) + cl_float(2184.0,precision)*Power(rxj,SIX) + 

              cl_float(496.0,precision)*Power(rxj,SEVEN) + cl_float(64.0,precision)*Power(rxj,EIGHT) + THREE*Power(rxj,NINE)) - 

           cl_float(18.0,precision)*Power(rxi,FOUR)*Power(rxj,cl_float(16.0,precision))*

            (cl_float(6572475.0,precision) - cl_float(3161340.0,precision)*rxj - cl_float(4782960.0,precision)*Power(rxj,TWO) - 

              cl_float(1912365.0,precision)*Power(rxj,THREE) - cl_float(378105.0,precision)*Power(rxj,FOUR) - cl_float(34125.0,precision)*Power(rxj,FIVE) + 

              cl_float(1092.0,precision)*Power(rxj,SIX) + cl_float(650.0,precision)*Power(rxj,SEVEN) + cl_float(71.0,precision)*Power(rxj,EIGHT) + 

              THREE*Power(rxj,NINE)) - cl_float(21.0,precision)*Power(rxi,EIGHT)*Power(rxj,cl_float(12.0,precision))*

            (-cl_float(1063800.0,precision) - cl_float(2775735.0,precision)*rxj - cl_float(862920.0,precision)*Power(rxj,TWO) - 

              cl_float(1132020.0,precision)*Power(rxj,THREE) - cl_float(698580.0,precision)*Power(rxj,FOUR) - 

              cl_float(196920.0,precision)*Power(rxj,FIVE) - cl_float(28992.0,precision)*Power(rxj,SIX) - cl_float(2064.0,precision)*Power(rxj,SEVEN) - 

              cl_float(24.0,precision)*Power(rxj,EIGHT) + FOUR*Power(rxj,NINE)) + 

           cl_float(21.0,precision)*Power(rxi,cl_float(12.0,precision))*Power(rxj,EIGHT)*

            (cl_float(482625.0,precision) + cl_float(868725.0,precision)*rxj + cl_float(772200.0,precision)*Power(rxj,TWO) + cl_float(455400.0,precision)*Power(rxj,THREE) + 

              cl_float(178200.0,precision)*Power(rxj,FOUR) + cl_float(72180.0,precision)*Power(rxj,FIVE) + cl_float(19920.0,precision)*Power(rxj,SIX) + 

              cl_float(2952.0,precision)*Power(rxj,SEVEN) + cl_float(204.0,precision)*Power(rxj,EIGHT) + FOUR*Power(rxj,NINE)) + 

           SIX*Power(rxi,SIX)*Power(rxj,cl_float(14.0,precision))*

            (-cl_float(10357200.0,precision) + cl_float(5071815.0,precision)*rxj - cl_float(6463800.0,precision)*Power(rxj,TWO) - 

              cl_float(7151130.0,precision)*Power(rxj,THREE) - cl_float(2572290.0,precision)*Power(rxj,FOUR) - 

              cl_float(468720.0,precision)*Power(rxj,FIVE) - cl_float(42672.0,precision)*Power(rxj,SIX) - cl_float(648.0,precision)*Power(rxj,SEVEN) + 

              cl_float(228.0,precision)*Power(rxj,EIGHT) + cl_float(16.0,precision)*Power(rxj,NINE)) - 

           Power(rxi,cl_float(18.0,precision))*Power(rxj,TWO)*

            (cl_float(184275.0,precision) + cl_float(331695.0,precision)*rxj + cl_float(294840.0,precision)*Power(rxj,TWO) + cl_float(171990.0,precision)*Power(rxj,THREE) + 

              cl_float(73710.0,precision)*Power(rxj,FOUR) + cl_float(24570.0,precision)*Power(rxj,FIVE) + cl_float(6552.0,precision)*Power(rxj,SIX) + 

              cl_float(1404.0,precision)*Power(rxj,SEVEN) + cl_float(234.0,precision)*Power(rxj,EIGHT) + cl_float(16.0,precision)*Power(rxj,NINE)) + 

           Power(rxi,TWO)*Power(rxj,cl_float(18.0,precision))*

            (-cl_float(133783650.0,precision) - cl_float(107432325.0,precision)*rxj - cl_float(35675640.0,precision)*Power(rxj,TWO) - 

              cl_float(5135130.0,precision)*Power(rxj,THREE) + cl_float(270270.0,precision)*Power(rxj,FOUR) + 

              cl_float(270270.0,precision)*Power(rxj,FIVE) + cl_float(57960.0,precision)*Power(rxj,SIX) + cl_float(6948.0,precision)*Power(rxj,SEVEN) + 

              cl_float(486.0,precision)*Power(rxj,EIGHT) + cl_float(16.0,precision)*Power(rxj,NINE)) - 

           SIX*Power(rxi,cl_float(14.0,precision))*Power(rxj,SIX)*

            (cl_float(675675.0,precision) + cl_float(1216215.0,precision)*rxj + cl_float(1081080.0,precision)*Power(rxj,TWO) + cl_float(630630.0,precision)*Power(rxj,THREE) + 

              cl_float(270270.0,precision)*Power(rxj,FOUR) + cl_float(88200.0,precision)*Power(rxj,FIVE) + cl_float(26544.0,precision)*Power(rxj,SIX) + 

              cl_float(5160.0,precision)*Power(rxj,SEVEN) + cl_float(492.0,precision)*Power(rxj,EIGHT) + cl_float(16.0,precision)*Power(rxj,NINE))))/

      (cl_float(28350.0,precision)*exp(TWO*(rxi + rxj))*Power(rxi - rxj,cl_float(13.0,precision))*Power(rxi + rxj,cl_float(13.0,precision)))

     ; }
   
  }
  return S;
}

cl_F Slater_2S_1S(cl_F r,cl_F xi,cl_F xj)
{
  return Slater_1S_2S(r,xj,xi);
}

cl_F Slater_3S_3S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = (cl_float(793.0,precision)*xi)/cl_float(3072.0,precision)

    ; } else {  S = (ONE/r)*(-cl_float(1437004800.0,precision) + cl_float(1437004800.0,precision)*exp(TWO*r*xi) - cl_float(2503064025.0,precision)*r*xi - 

        cl_float(2132118450.0,precision)*Power(r,TWO)*Power(xi,TWO) - cl_float(1180664100.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

        cl_float(476506800.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - cl_float(148856400.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

        cl_float(37255680.0,precision)*Power(r,SIX)*Power(xi,SIX) - cl_float(7603200.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

        cl_float(1267200.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - cl_float(168960.0,precision)*Power(r,NINE)*Power(xi,NINE) - 

        cl_float(16896.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - cl_float(1024.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)))/

      (cl_float(1.4370048e9,precision)*exp(TWO*r*xi))

    ; }
 
  }
  else {
      if (r == ZERO) { S = (xi*xj*(Power(xi,cl_float(10.0,precision)) + cl_float(11.0,precision)*Power(xi,NINE)*xj + cl_float(55.0,precision)*Power(xi,EIGHT)*Power(xj,TWO) + 

          cl_float(165.0,precision)*Power(xi,SEVEN)*Power(xj,THREE) + cl_float(330.0,precision)*Power(xi,SIX)*Power(xj,FOUR) + 

          cl_float(462.0,precision)*Power(xi,FIVE)*Power(xj,FIVE) + cl_float(330.0,precision)*Power(xi,FOUR)*Power(xj,SIX) + 

          cl_float(165.0,precision)*Power(xi,THREE)*Power(xj,SEVEN) + cl_float(55.0,precision)*Power(xi,TWO)*Power(xj,EIGHT) + 

          cl_float(11.0,precision)*xi*Power(xj,NINE) + Power(xj,cl_float(10.0,precision))))/(cl_float(3.0,precision)*Power(xi + xj,cl_float(11.0,precision)))

    ; } else { S = (ONE/r)*(cl_float(135.0,precision)*exp(TWO*(rxi + rxj))*Power(Power(rxi,TWO) - Power(rxj,TWO),cl_float(11.0,precision)) + 

        exp(TWO*rxj)*Power(rxj,EIGHT)*

         (-cl_float(150.0,precision)*Power(rxi,cl_float(18.0,precision)) - SIX*Power(rxi,cl_float(19.0,precision)) + cl_float(135.0,precision)*Power(rxj,cl_float(14.0,precision)) + 

           cl_float(225.0,precision)*rxi*Power(rxj,cl_float(14.0,precision)) + cl_float(10.0,precision)*Power(rxi,cl_float(17.0,precision))*(-cl_float(165.0,precision) + Power(rxj,TWO)) - 

           cl_float(30.0,precision)*Power(rxi,cl_float(16.0,precision))*(cl_float(330.0,precision) + Power(rxj,TWO)) + 

           cl_float(45.0,precision)*Power(rxi,THREE)*Power(rxj,cl_float(12.0,precision))*(-cl_float(55.0,precision) + TWO*Power(rxj,TWO)) + 

           cl_float(45.0,precision)*Power(rxi,TWO)*Power(rxj,cl_float(12.0,precision))*(-cl_float(33.0,precision) + FOUR*Power(rxj,TWO)) + 

           Power(rxi,NINE)*Power(rxj,SIX)*

            (cl_float(234135.0,precision) - cl_float(4950.0,precision)*Power(rxj,TWO) - cl_float(34.0,precision)*Power(rxj,FOUR)) - 

           FIVE*Power(rxi,SEVEN)*Power(rxj,EIGHT)*

            (cl_float(6237.0,precision) - cl_float(1242.0,precision)*Power(rxj,TWO) + TWO*Power(rxj,FOUR)) + 

           THREE*Power(rxi,FIVE)*Power(rxj,cl_float(10.0,precision))*

            (cl_float(4125.0,precision) - cl_float(330.0,precision)*Power(rxj,TWO) + TWO*Power(rxj,FOUR)) + 

           cl_float(15.0,precision)*Power(rxi,FOUR)*Power(rxj,cl_float(10.0,precision))*

            (cl_float(495.0,precision) - cl_float(132.0,precision)*Power(rxj,TWO) + TWO*Power(rxj,FOUR)) - 

           cl_float(165.0,precision)*Power(rxi,SIX)*Power(rxj,EIGHT)*

            (cl_float(135.0,precision) - cl_float(60.0,precision)*Power(rxj,TWO) + TWO*Power(rxj,FOUR)) - 

           FIVE*Power(rxi,cl_float(13.0,precision))*Power(rxj,TWO)*

            (cl_float(43875.0,precision) - cl_float(3438.0,precision)*Power(rxj,TWO) + cl_float(22.0,precision)*Power(rxj,FOUR)) + 

           FIVE*Power(rxi,cl_float(11.0,precision))*Power(rxj,FOUR)*

            (cl_float(7695.0,precision) - cl_float(2442.0,precision)*Power(rxj,TWO) + cl_float(22.0,precision)*Power(rxj,FOUR)) + 

           cl_float(15.0,precision)*Power(rxi,EIGHT)*Power(rxj,SIX)*

            (-cl_float(33.0,precision) - cl_float(3564.0,precision)*Power(rxj,TWO) + cl_float(26.0,precision)*Power(rxj,FOUR)) + 

           Power(rxi,cl_float(15.0,precision))*(-cl_float(32175.0,precision) - cl_float(3690.0,precision)*Power(rxj,TWO) + cl_float(34.0,precision)*Power(rxj,FOUR)) + 

           cl_float(15.0,precision)*Power(rxi,cl_float(10.0,precision))*Power(rxj,FOUR)*

            (-cl_float(32277.0,precision) + cl_float(1364.0,precision)*Power(rxj,TWO) + cl_float(66.0,precision)*Power(rxj,FOUR)) + 

           cl_float(15.0,precision)*Power(rxi,cl_float(14.0,precision))*(-cl_float(3003.0,precision) - cl_float(2932.0,precision)*Power(rxj,TWO) + cl_float(94.0,precision)*Power(rxj,FOUR)) - 

           cl_float(15.0,precision)*Power(rxi,cl_float(12.0,precision))*Power(rxj,TWO)*

            (cl_float(28119.0,precision) - cl_float(5252.0,precision)*Power(rxj,TWO) + cl_float(154.0,precision)*Power(rxj,FOUR))) + 

        exp(TWO*rxi)*Power(rxi,EIGHT)*

         (-FIVE*Power(rxi,TWO)*Power(rxj,cl_float(12.0,precision))*

            (-cl_float(84357.0,precision) - cl_float(43875.0,precision)*rxj - cl_float(8796.0,precision)*Power(rxj,TWO) - cl_float(738.0,precision)*Power(rxj,THREE) - 

              SIX*Power(rxj,FOUR) + TWO*Power(rxj,FIVE)) - 

           THREE*Power(rxi,cl_float(14.0,precision))*(cl_float(45.0,precision) + cl_float(75.0,precision)*rxj + cl_float(60.0,precision)*Power(rxj,TWO) + cl_float(30.0,precision)*Power(rxj,THREE) + 

              cl_float(10.0,precision)*Power(rxj,FOUR) + TWO*Power(rxj,FIVE)) - 

           cl_float(55.0,precision)*Power(rxi,EIGHT)*Power(rxj,SIX)*

            (-cl_float(405.0,precision) - cl_float(567.0,precision)*rxj - cl_float(972.0,precision)*Power(rxj,TWO) - cl_float(90.0,precision)*Power(rxj,THREE) + 

              cl_float(18.0,precision)*Power(rxj,FOUR) + TWO*Power(rxj,FIVE)) + 

           cl_float(55.0,precision)*Power(rxi,SIX)*Power(rxj,EIGHT)*

            (NINE - cl_float(4257.0,precision)*rxj - cl_float(372.0,precision)*Power(rxj,TWO) + cl_float(222.0,precision)*Power(rxj,THREE) + 

              cl_float(42.0,precision)*Power(rxj,FOUR) + TWO*Power(rxj,FIVE)) + 

           THREE*Power(rxj,cl_float(14.0,precision))*(cl_float(15015.0,precision) + cl_float(10725.0,precision)*rxj + cl_float(3300.0,precision)*Power(rxj,TWO) + 

              cl_float(550.0,precision)*Power(rxj,THREE) + cl_float(50.0,precision)*Power(rxj,FOUR) + TWO*Power(rxj,FIVE)) + 

           FIVE*Power(rxi,cl_float(12.0,precision))*Power(rxj,TWO)*

            (cl_float(297.0,precision) + cl_float(495.0,precision)*rxj + cl_float(396.0,precision)*Power(rxj,TWO) + cl_float(198.0,precision)*Power(rxj,THREE) + 

              cl_float(66.0,precision)*Power(rxj,FOUR) + TWO*Power(rxj,FIVE)) + 

           Power(rxi,cl_float(10.0,precision))*Power(rxj,FOUR)*

            (-cl_float(7425.0,precision) - cl_float(12375.0,precision)*rxj - cl_float(9900.0,precision)*Power(rxj,TWO) - cl_float(6210.0,precision)*Power(rxj,THREE) - 

              cl_float(390.0,precision)*Power(rxj,FOUR) + cl_float(34.0,precision)*Power(rxj,FIVE)) - 

           Power(rxi,FOUR)*Power(rxj,cl_float(10.0,precision))*

            (-cl_float(484155.0,precision) + cl_float(38475.0,precision)*rxj + cl_float(78780.0,precision)*Power(rxj,TWO) + cl_float(17190.0,precision)*Power(rxj,THREE) + 

              cl_float(1410.0,precision)*Power(rxj,FOUR) + cl_float(34.0,precision)*Power(rxj,FIVE))))/

      (cl_float(135.0,precision)*exp(TWO*(rxi + rxj))*Power(rxi - rxj,cl_float(11.0,precision))*Power(rxi + rxj,cl_float(11.0,precision)))

     ; }
   
  }
  return S;
}

cl_F Slater_3S_4S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = (cl_float(1363.0,precision)*xi)/cl_float(6144.0,precision)

    ; } else {  S = (ONE/r)*(-cl_float(74724249600.0,precision) + cl_float(74724249600.0,precision)*exp(TWO*r*xi) - cl_float(132871488750.0,precision)*r*xi - 

        cl_float(116294478300.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

        cl_float(66678987375.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

        cl_float(28114836750.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - cl_float(9274044780.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

        cl_float(2484321840.0,precision)*Power(r,SIX)*Power(xi,SIX) - cl_float(553204080.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

        cl_float(103783680.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - cl_float(16473600.0,precision)*Power(r,NINE)*Power(xi,NINE) - 

        cl_float(2196480.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - cl_float(239616.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)) - 

        cl_float(19968.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(12.0,precision)) - cl_float(1024.0,precision)*Power(r,cl_float(13.0,precision))*Power(xi,cl_float(13.0,precision)))/

      (cl_float(7.47242496e10,precision)*exp(TWO*r*xi))

    ; }
 
  }
  else {
      if (r == ZERO) { S = (xi*xj*(THREE*Power(xi,cl_float(12.0,precision)) + cl_float(39.0,precision)*Power(xi,cl_float(11.0,precision))*xj + cl_float(234.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,TWO) + 

          cl_float(858.0,precision)*Power(xi,NINE)*Power(xj,THREE) + cl_float(2145.0,precision)*Power(xi,EIGHT)*Power(xj,FOUR) + 

          cl_float(3861.0,precision)*Power(xi,SEVEN)*Power(xj,FIVE) + cl_float(5148.0,precision)*Power(xi,SIX)*Power(xj,SIX) + 

          cl_float(5148.0,precision)*Power(xi,FIVE)*Power(xj,SEVEN) + cl_float(2860.0,precision)*Power(xi,FOUR)*Power(xj,EIGHT) + 

          cl_float(1144.0,precision)*Power(xi,THREE)*Power(xj,NINE) + cl_float(312.0,precision)*Power(xi,TWO)*Power(xj,cl_float(10.0,precision)) + 

          cl_float(52.0,precision)*xi*Power(xj,cl_float(11.0,precision)) + FOUR*Power(xj,cl_float(12.0,precision))))/(cl_float(12.0,precision)*Power(xi + xj,cl_float(13.0,precision)))

    ; } else { S = (ONE/r)*(cl_float(3780.0,precision)*exp(TWO*(rxi + rxj))*Power(Power(rxi,TWO) - Power(rxj,TWO),cl_float(13.0,precision)) + 

        cl_float(84.0,precision)*exp(TWO*rxj)*Power(rxj,cl_float(10.0,precision))*

         (-cl_float(60.0,precision)*Power(rxi,cl_float(20.0,precision)) - TWO*Power(rxi,cl_float(21.0,precision)) + cl_float(45.0,precision)*Power(rxj,cl_float(16.0,precision)) + 

           cl_float(75.0,precision)*rxi*Power(rxj,cl_float(16.0,precision)) - FOUR*Power(rxi,cl_float(19.0,precision))*(cl_float(195.0,precision) + Power(rxj,TWO)) + 

           cl_float(15.0,precision)*Power(rxi,THREE)*Power(rxj,cl_float(14.0,precision))*(-cl_float(65.0,precision) + TWO*Power(rxj,TWO)) + 

           cl_float(15.0,precision)*Power(rxi,TWO)*Power(rxj,cl_float(14.0,precision))*(-cl_float(39.0,precision) + FOUR*Power(rxj,TWO)) - 

           cl_float(30.0,precision)*Power(rxi,cl_float(18.0,precision))*(cl_float(182.0,precision) + NINE*Power(rxj,TWO)) + 

           cl_float(30.0,precision)*Power(rxi,cl_float(13.0,precision))*Power(rxj,FOUR)*(-cl_float(13047.0,precision) + cl_float(377.0,precision)*Power(rxj,TWO)) + 

           TWO*Power(rxi,FIVE)*Power(rxj,cl_float(12.0,precision))*

            (cl_float(2925.0,precision) - cl_float(195.0,precision)*Power(rxj,TWO) + Power(rxj,FOUR)) + 

           cl_float(10.0,precision)*Power(rxi,FOUR)*Power(rxj,cl_float(12.0,precision))*

            (cl_float(351.0,precision) - cl_float(78.0,precision)*Power(rxj,TWO) + Power(rxj,FOUR)) - 

           cl_float(130.0,precision)*Power(rxi,SIX)*Power(rxj,cl_float(10.0,precision))*

            (cl_float(99.0,precision) - cl_float(36.0,precision)*Power(rxj,TWO) + Power(rxj,FOUR)) + 

           cl_float(13.0,precision)*Power(rxi,cl_float(11.0,precision))*Power(rxj,SIX)*

            (cl_float(30735.0,precision) - cl_float(1650.0,precision)*Power(rxj,TWO) + FOUR*Power(rxj,FOUR)) + 

           Power(rxi,SEVEN)*Power(rxj,cl_float(10.0,precision))*

            (-cl_float(15015.0,precision) + cl_float(3330.0,precision)*Power(rxj,TWO) + FOUR*Power(rxj,FOUR)) + 

           cl_float(210.0,precision)*Power(rxi,cl_float(16.0,precision))*(-cl_float(156.0,precision) - cl_float(262.0,precision)*Power(rxj,TWO) + FIVE*Power(rxj,FOUR)) - 

           SIX*Power(rxi,NINE)*Power(rxj,EIGHT)*

            (-cl_float(48620.0,precision) - cl_float(715.0,precision)*Power(rxj,TWO) + SIX*Power(rxj,FOUR)) + 

           THREE*Power(rxi,cl_float(17.0,precision))*(-cl_float(6825.0,precision) - cl_float(1870.0,precision)*Power(rxj,TWO) + cl_float(12.0,precision)*Power(rxj,FOUR)) - 

           cl_float(30.0,precision)*Power(rxi,cl_float(14.0,precision))*Power(rxj,TWO)*

            (cl_float(17934.0,precision) - cl_float(12.0,precision)*Power(rxj,TWO) + cl_float(13.0,precision)*Power(rxj,FOUR)) - 

           cl_float(15.0,precision)*Power(rxi,EIGHT)*Power(rxj,EIGHT)*

            (cl_float(2145.0,precision) + cl_float(2860.0,precision)*Power(rxj,TWO) + cl_float(14.0,precision)*Power(rxj,FOUR)) + 

           cl_float(65.0,precision)*Power(rxi,cl_float(10.0,precision))*Power(rxj,SIX)*

            (-cl_float(13725.0,precision) - cl_float(792.0,precision)*Power(rxj,TWO) + cl_float(22.0,precision)*Power(rxj,FOUR)) - 

           cl_float(10.0,precision)*Power(rxi,cl_float(12.0,precision))*Power(rxj,FOUR)*

            (cl_float(153630.0,precision) - cl_float(15054.0,precision)*Power(rxj,TWO) + cl_float(143.0,precision)*Power(rxj,FOUR)) + 

           Power(rxi,cl_float(15.0,precision))*(-cl_float(269325.0,precision)*Power(rxj,TWO) + cl_float(9270.0,precision)*Power(rxj,FOUR) - 

              cl_float(52.0,precision)*Power(rxj,SIX))) + exp(TWO*rxi)*Power(rxi,EIGHT)*

         (Power(rxi,TWO)*Power(rxj,cl_float(16.0,precision))*

            (cl_float(70073640.0,precision) + cl_float(47669895.0,precision)*rxj + cl_float(13931190.0,precision)*Power(rxj,TWO) + 

              cl_float(2170350.0,precision)*Power(rxj,THREE) + cl_float(169260.0,precision)*Power(rxj,FOUR) + cl_float(1638.0,precision)*Power(rxj,FIVE) - 

              cl_float(756.0,precision)*Power(rxj,SIX) - cl_float(44.0,precision)*Power(rxj,SEVEN)) + 

           cl_float(364.0,precision)*Power(rxi,cl_float(10.0,precision))*Power(rxj,EIGHT)*

            (-cl_float(7425.0,precision) - cl_float(13860.0,precision)*rxj - cl_float(5940.0,precision)*Power(rxj,TWO) - cl_float(11880.0,precision)*Power(rxj,THREE) - 

              cl_float(2640.0,precision)*Power(rxj,FOUR) - cl_float(45.0,precision)*Power(rxj,FIVE) + cl_float(30.0,precision)*Power(rxj,SIX) + 

              TWO*Power(rxj,SEVEN)) - cl_float(364.0,precision)*Power(rxi,EIGHT)*Power(rxj,cl_float(10.0,precision))*

            (-cl_float(20925.0,precision) + cl_float(18270.0,precision)*rxj - cl_float(58320.0,precision)*Power(rxj,TWO) - cl_float(17730.0,precision)*Power(rxj,THREE) - 

              cl_float(300.0,precision)*Power(rxj,FOUR) + cl_float(423.0,precision)*Power(rxj,FIVE) + cl_float(54.0,precision)*Power(rxj,SIX) + 

              TWO*Power(rxj,SEVEN)) - THREE*Power(rxi,cl_float(18.0,precision))*

            (cl_float(1260.0,precision) + cl_float(2205.0,precision)*rxj + cl_float(1890.0,precision)*Power(rxj,TWO) + cl_float(1050.0,precision)*Power(rxj,THREE) + 

              cl_float(420.0,precision)*Power(rxj,FOUR) + cl_float(126.0,precision)*Power(rxj,FIVE) + cl_float(28.0,precision)*Power(rxj,SIX) + 

              FOUR*Power(rxj,SEVEN)) + THREE*Power(rxj,cl_float(18.0,precision))*

            (cl_float(1801800.0,precision) + cl_float(1576575.0,precision)*rxj + cl_float(630630.0,precision)*Power(rxj,TWO) + 

              cl_float(150150.0,precision)*Power(rxj,THREE) + cl_float(23100.0,precision)*Power(rxj,FOUR) + cl_float(2310.0,precision)*Power(rxj,FIVE) + 

              cl_float(140.0,precision)*Power(rxj,SIX) + FOUR*Power(rxj,SEVEN)) + 

           TWO*Power(rxi,cl_float(14.0,precision))*Power(rxj,FOUR)*

            (-cl_float(147420.0,precision) - cl_float(257985.0,precision)*rxj - cl_float(221130.0,precision)*Power(rxj,TWO) - cl_float(122850.0,precision)*Power(rxj,THREE) - 

              cl_float(49140.0,precision)*Power(rxj,FOUR) - cl_float(17388.0,precision)*Power(rxj,FIVE) - cl_float(1512.0,precision)*Power(rxj,SIX) + 

              EIGHT*Power(rxj,SEVEN)) - cl_float(42.0,precision)*Power(rxi,cl_float(12.0,precision))*Power(rxj,SIX)*

            (-cl_float(25740.0,precision) - cl_float(45045.0,precision)*rxj - cl_float(38610.0,precision)*Power(rxj,TWO) - cl_float(19470.0,precision)*Power(rxj,THREE) - 

              cl_float(12540.0,precision)*Power(rxj,FOUR) - cl_float(1836.0,precision)*Power(rxj,FIVE) - EIGHT*Power(rxj,SIX) + 

              EIGHT*Power(rxj,SEVEN)) + cl_float(42.0,precision)*Power(rxi,SIX)*Power(rxj,cl_float(12.0,precision))*

            (cl_float(921600.0,precision) - cl_float(1640835.0,precision)*rxj - cl_float(546030.0,precision)*Power(rxj,TWO) + cl_float(20730.0,precision)*Power(rxj,THREE) + 

              cl_float(30180.0,precision)*Power(rxj,FOUR) + cl_float(5028.0,precision)*Power(rxj,FIVE) + cl_float(344.0,precision)*Power(rxj,SIX) + 

              EIGHT*Power(rxj,SEVEN)) - TWO*Power(rxi,FOUR)*Power(rxj,cl_float(14.0,precision))*

            (-cl_float(67767840.0,precision) - cl_float(13377735.0,precision)*rxj + cl_float(6601770.0,precision)*Power(rxj,TWO) + 

              cl_float(3115350.0,precision)*Power(rxj,THREE) + cl_float(548940.0,precision)*Power(rxj,FOUR) + cl_float(48132.0,precision)*Power(rxj,FIVE) + 

              cl_float(1848.0,precision)*Power(rxj,SIX) + EIGHT*Power(rxj,SEVEN)) + 

           Power(rxi,cl_float(16.0,precision))*Power(rxj,TWO)*

            (cl_float(49140.0,precision) + cl_float(85995.0,precision)*rxj + cl_float(73710.0,precision)*Power(rxj,TWO) + cl_float(40950.0,precision)*Power(rxj,THREE) + 

              cl_float(16380.0,precision)*Power(rxj,FOUR) + cl_float(4914.0,precision)*Power(rxj,FIVE) + cl_float(1092.0,precision)*Power(rxj,SIX) + 

              cl_float(44.0,precision)*Power(rxj,SEVEN))))/

      (cl_float(3780.0,precision)*exp(TWO*(rxi + rxj))*Power(rxi - rxj,cl_float(13.0,precision))*Power(rxi + rxj,cl_float(13.0,precision)))

     ; }
   
  }
  return S;
}

cl_F Slater_3S_5S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = (cl_float(31059.0,precision)*xi)/cl_float(163840.0,precision)

    ; } else {  S = (ONE/r)*(-cl_float(313841848320000.0,precision) + cl_float(313841848320000.0,precision)*exp(TWO*r*xi) - cl_float(568188982486125.0,precision)*r*xi - 

        cl_float(508694268332250.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

        cl_float(299892470377500.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

        cl_float(130753815192000.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - 

        cl_float(44881155118800.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

        cl_float(12601803614400.0,precision)*Power(r,SIX)*Power(xi,SIX) - 

        cl_float(2967953788800.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

        cl_float(596237241600.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - 

        cl_float(103264761600.0,precision)*Power(r,NINE)*Power(xi,NINE) - 

        cl_float(15498362880.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - 

        cl_float(2012774400.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)) - 

        cl_float(223641600.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(12.0,precision)) - cl_float(20643840.0,precision)*Power(r,cl_float(13.0,precision))*Power(xi,cl_float(13.0,precision)) - 

        cl_float(1474560.0,precision)*Power(r,cl_float(14.0,precision))*Power(xi,cl_float(14.0,precision)) - cl_float(65536.0,precision)*Power(r,cl_float(15.0,precision))*Power(xi,cl_float(15.0,precision)))/

      (cl_float(3.1384184832e14,precision)*exp(TWO*r*xi))

    ; }
 
  }
  else {
      if (r == ZERO) { S = (xi*xj*(THREE*Power(xi,cl_float(14.0,precision)) + cl_float(45.0,precision)*Power(xi,cl_float(13.0,precision))*xj + cl_float(315.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,TWO) + 

          cl_float(1365.0,precision)*Power(xi,cl_float(11.0,precision))*Power(xj,THREE) + cl_float(4095.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,FOUR) + 

          cl_float(9009.0,precision)*Power(xi,NINE)*Power(xj,FIVE) + cl_float(15015.0,precision)*Power(xi,EIGHT)*Power(xj,SIX) + 

          cl_float(19305.0,precision)*Power(xi,SEVEN)*Power(xj,SEVEN) + cl_float(19305.0,precision)*Power(xi,SIX)*Power(xj,EIGHT) + 

          cl_float(15015.0,precision)*Power(xi,FIVE)*Power(xj,NINE) + cl_float(6825.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(10.0,precision)) + 

          cl_float(2275.0,precision)*Power(xi,THREE)*Power(xj,cl_float(11.0,precision)) + cl_float(525.0,precision)*Power(xi,TWO)*Power(xj,cl_float(12.0,precision)) + 

          cl_float(75.0,precision)*xi*Power(xj,cl_float(13.0,precision)) + FIVE*Power(xj,cl_float(14.0,precision))))/(cl_float(15.0,precision)*Power(xi + xj,cl_float(15.0,precision)))

    ; } else { S = (ONE/r)*(cl_float(42525.0,precision)*exp(TWO*(rxi + rxj))*Power(Power(rxi,TWO) - Power(rxj,TWO),cl_float(15.0,precision)) + 

        cl_float(189.0,precision)*exp(TWO*rxj)*Power(rxj,cl_float(12.0,precision))*

         (-cl_float(350.0,precision)*Power(rxi,cl_float(22.0,precision)) - cl_float(10.0,precision)*Power(rxi,cl_float(23.0,precision)) + cl_float(225.0,precision)*Power(rxj,cl_float(18.0,precision)) + 

           cl_float(375.0,precision)*rxi*Power(rxj,cl_float(18.0,precision)) - cl_float(70.0,precision)*Power(rxi,cl_float(21.0,precision))*(cl_float(75.0,precision) + Power(rxj,TWO)) + 

           cl_float(75.0,precision)*Power(rxi,THREE)*Power(rxj,cl_float(16.0,precision))*(-cl_float(75.0,precision) + TWO*Power(rxj,TWO)) + 

           cl_float(75.0,precision)*Power(rxi,TWO)*Power(rxj,cl_float(16.0,precision))*(-cl_float(45.0,precision) + FOUR*Power(rxj,TWO)) - 

           cl_float(50.0,precision)*Power(rxi,cl_float(20.0,precision))*(cl_float(840.0,precision) + cl_float(71.0,precision)*Power(rxj,TWO)) + 

           Power(rxi,NINE)*Power(rxj,cl_float(10.0,precision))*

            (cl_float(4694625.0,precision) + cl_float(124800.0,precision)*Power(rxj,TWO) - cl_float(248.0,precision)*Power(rxj,FOUR)) + 

           cl_float(20.0,precision)*Power(rxi,cl_float(17.0,precision))*Power(rxj,TWO)*

            (-cl_float(185895.0,precision) - cl_float(948.0,precision)*Power(rxj,TWO) + TWO*Power(rxj,FOUR)) + 

           FIVE*Power(rxi,FIVE)*Power(rxj,cl_float(14.0,precision))*

            (cl_float(7875.0,precision) - cl_float(450.0,precision)*Power(rxj,TWO) + TWO*Power(rxj,FOUR)) + 

           cl_float(25.0,precision)*Power(rxi,FOUR)*Power(rxj,cl_float(14.0,precision))*

            (cl_float(945.0,precision) - cl_float(180.0,precision)*Power(rxj,TWO) + TWO*Power(rxj,FOUR)) - 

           cl_float(375.0,precision)*Power(rxi,SIX)*Power(rxj,cl_float(12.0,precision))*

            (cl_float(273.0,precision) - cl_float(84.0,precision)*Power(rxj,TWO) + TWO*Power(rxj,FOUR)) - 

           FIVE*Power(rxi,cl_float(11.0,precision))*Power(rxj,EIGHT)*

            (-cl_float(2803125.0,precision) + cl_float(49140.0,precision)*Power(rxj,TWO) + EIGHT*Power(rxj,FOUR)) + 

           FIVE*Power(rxi,SEVEN)*Power(rxj,cl_float(12.0,precision))*

            (-cl_float(16965.0,precision) + cl_float(5152.0,precision)*Power(rxj,TWO) + cl_float(14.0,precision)*Power(rxj,FOUR)) + 

           cl_float(325.0,precision)*Power(rxi,cl_float(10.0,precision))*Power(rxj,EIGHT)*

            (-cl_float(60117.0,precision) - cl_float(5340.0,precision)*Power(rxj,TWO) + cl_float(40.0,precision)*Power(rxj,FOUR)) - 

           cl_float(15.0,precision)*Power(rxi,cl_float(15.0,precision))*Power(rxj,FOUR)*

            (cl_float(845085.0,precision) - cl_float(22960.0,precision)*Power(rxj,TWO) + cl_float(52.0,precision)*Power(rxj,FOUR)) + 

           cl_float(15.0,precision)*Power(rxi,cl_float(13.0,precision))*Power(rxj,SIX)*

            (-cl_float(139125.0,precision) - cl_float(10140.0,precision)*Power(rxj,TWO) + cl_float(52.0,precision)*Power(rxj,FOUR)) + 

           cl_float(75.0,precision)*Power(rxi,cl_float(12.0,precision))*Power(rxj,SIX)*

            (-cl_float(729687.0,precision) + cl_float(25532.0,precision)*Power(rxj,TWO) + cl_float(52.0,precision)*Power(rxj,FOUR)) + 

           cl_float(60.0,precision)*Power(rxi,cl_float(18.0,precision))*(-cl_float(5355.0,precision) - cl_float(11940.0,precision)*Power(rxj,TWO) + cl_float(86.0,precision)*Power(rxj,FOUR)) + 

           TWO*Power(rxi,cl_float(19.0,precision))*(-cl_float(89250.0,precision) - cl_float(35425.0,precision)*Power(rxj,TWO) + cl_float(124.0,precision)*Power(rxj,FOUR)) + 

           cl_float(100.0,precision)*Power(rxi,cl_float(16.0,precision))*Power(rxj,TWO)*

            (-cl_float(79713.0,precision) - cl_float(13311.0,precision)*Power(rxj,TWO) + cl_float(146.0,precision)*Power(rxj,FOUR)) - 

           FIVE*Power(rxi,EIGHT)*Power(rxj,cl_float(10.0,precision))*

            (cl_float(157365.0,precision) + cl_float(95940.0,precision)*Power(rxj,TWO) + cl_float(952.0,precision)*Power(rxj,FOUR)) - 

           cl_float(15.0,precision)*Power(rxi,cl_float(14.0,precision))*Power(rxj,FOUR)*

            (cl_float(2638467.0,precision) - cl_float(157500.0,precision)*Power(rxj,TWO) + cl_float(1820.0,precision)*Power(rxj,FOUR))) - 

        exp(TWO*rxi)*Power(rxi,EIGHT)*

         (THREE*Power(rxi,cl_float(14.0,precision))*Power(rxj,EIGHT)*

            (cl_float(19348875.0,precision) + cl_float(34827975.0,precision)*rxj + cl_float(30958200.0,precision)*Power(rxj,TWO) + 

              cl_float(18689580.0,precision)*Power(rxj,THREE) + cl_float(5847660.0,precision)*Power(rxj,FOUR) + 

              cl_float(3723300.0,precision)*Power(rxj,FIVE) + cl_float(845040.0,precision)*Power(rxj,SIX) + cl_float(58680.0,precision)*Power(rxj,SEVEN) - 

              cl_float(1548.0,precision)*Power(rxj,EIGHT) - cl_float(236.0,precision)*Power(rxj,NINE)) - 

           cl_float(42.0,precision)*Power(rxi,FOUR)*Power(rxj,cl_float(18.0,precision))*

            (cl_float(251336925.0,precision) + cl_float(104824125.0,precision)*rxj + cl_float(340200.0,precision)*Power(rxj,TWO) - 

              cl_float(9122085.0,precision)*Power(rxj,THREE) - cl_float(2798145.0,precision)*Power(rxj,FOUR) - 

              cl_float(433755.0,precision)*Power(rxj,FIVE) - cl_float(39060.0,precision)*Power(rxj,SIX) - cl_float(1890.0,precision)*Power(rxj,SEVEN) - 

              cl_float(27.0,precision)*Power(rxj,EIGHT) + Power(rxj,NINE)) - 

           SIX*Power(rxj,cl_float(22.0,precision))*(cl_float(34459425.0,precision) + cl_float(34459425.0,precision)*rxj + cl_float(16216200.0,precision)*Power(rxj,TWO) + 

              cl_float(4729725.0,precision)*Power(rxj,THREE) + cl_float(945945.0,precision)*Power(rxj,FOUR) + 

              cl_float(135135.0,precision)*Power(rxj,FIVE) + cl_float(13860.0,precision)*Power(rxj,SIX) + cl_float(990.0,precision)*Power(rxj,SEVEN) + 

              cl_float(45.0,precision)*Power(rxj,EIGHT) + Power(rxj,NINE)) + 

           THREE*Power(rxi,cl_float(22.0,precision))*(cl_float(14175.0,precision) + cl_float(25515.0,precision)*rxj + cl_float(22680.0,precision)*Power(rxj,TWO) + 

              cl_float(13230.0,precision)*Power(rxj,THREE) + cl_float(5670.0,precision)*Power(rxj,FOUR) + cl_float(1890.0,precision)*Power(rxj,FIVE) + 

              cl_float(504.0,precision)*Power(rxj,SIX) + cl_float(108.0,precision)*Power(rxj,SEVEN) + cl_float(18.0,precision)*Power(rxj,EIGHT) + 

              TWO*Power(rxj,NINE)) + cl_float(21.0,precision)*Power(rxi,cl_float(18.0,precision))*Power(rxj,FOUR)*

            (cl_float(212625.0,precision) + cl_float(382725.0,precision)*rxj + cl_float(340200.0,precision)*Power(rxj,TWO) + cl_float(198450.0,precision)*Power(rxj,THREE) + 

              cl_float(85050.0,precision)*Power(rxj,FOUR) + cl_float(28350.0,precision)*Power(rxj,FIVE) + cl_float(7560.0,precision)*Power(rxj,SIX) + 

              cl_float(1836.0,precision)*Power(rxj,SEVEN) + cl_float(162.0,precision)*Power(rxj,EIGHT) + TWO*Power(rxj,NINE)) - 

           cl_float(54.0,precision)*Power(rxi,SIX)*Power(rxj,cl_float(16.0,precision))*

            (cl_float(133451955.0,precision) - cl_float(73700865.0,precision)*rxj - cl_float(54096840.0,precision)*Power(rxj,TWO) - 

              cl_float(8306235.0,precision)*Power(rxj,THREE) + cl_float(966945.0,precision)*Power(rxj,FOUR) + 

              cl_float(516747.0,precision)*Power(rxj,FIVE) + cl_float(80724.0,precision)*Power(rxj,SIX) + cl_float(6434.0,precision)*Power(rxj,SEVEN) + 

              cl_float(251.0,precision)*Power(rxj,EIGHT) + THREE*Power(rxj,NINE)) + 

           cl_float(315.0,precision)*Power(rxi,cl_float(12.0,precision))*Power(rxj,cl_float(10.0,precision))*

            (-cl_float(405405.0,precision) - cl_float(710073.0,precision)*rxj - cl_float(805896.0,precision)*Power(rxj,TWO) - cl_float(101556.0,precision)*Power(rxj,THREE) - 

              cl_float(258804.0,precision)*Power(rxj,FOUR) - cl_float(90972.0,precision)*Power(rxj,FIVE) - cl_float(9744.0,precision)*Power(rxj,SIX) + 

              cl_float(120.0,precision)*Power(rxj,SEVEN) + cl_float(84.0,precision)*Power(rxj,EIGHT) + FOUR*Power(rxj,NINE)) - 

           cl_float(315.0,precision)*Power(rxi,cl_float(10.0,precision))*Power(rxj,cl_float(12.0,precision))*

            (-cl_float(482895.0,precision) - cl_float(2656395.0,precision)*rxj + cl_float(1186920.0,precision)*Power(rxj,TWO) - 

              cl_float(1155420.0,precision)*Power(rxj,THREE) - cl_float(643356.0,precision)*Power(rxj,FOUR) - cl_float(93492.0,precision)*Power(rxj,FIVE) + 

              cl_float(336.0,precision)*Power(rxj,SIX) + cl_float(1368.0,precision)*Power(rxj,SEVEN) + cl_float(132.0,precision)*Power(rxj,EIGHT) + 

              FOUR*Power(rxj,NINE)) + cl_float(27.0,precision)*Power(rxi,cl_float(16.0,precision))*Power(rxj,SIX)*

            (-cl_float(716625.0,precision) - cl_float(1289925.0,precision)*rxj - cl_float(1146600.0,precision)*Power(rxj,TWO) - 

              cl_float(668850.0,precision)*Power(rxj,THREE) - cl_float(286650.0,precision)*Power(rxj,FOUR) - cl_float(90006.0,precision)*Power(rxj,FIVE) - 

              cl_float(32872.0,precision)*Power(rxj,SIX) - cl_float(4812.0,precision)*Power(rxj,SEVEN) - cl_float(178.0,precision)*Power(rxj,EIGHT) + 

              SIX*Power(rxj,NINE)) + TWO*Power(rxi,TWO)*Power(rxj,cl_float(20.0,precision))*

            (-cl_float(1782492075.0,precision) - cl_float(1449175455.0,precision)*rxj - cl_float(533365560.0,precision)*Power(rxj,TWO) - 

              cl_float(114631335.0,precision)*Power(rxj,THREE) - cl_float(15221115.0,precision)*Power(rxj,FOUR) - 

              cl_float(1142505.0,precision)*Power(rxj,FIVE) - cl_float(18396.0,precision)*Power(rxj,SIX) + cl_float(5238.0,precision)*Power(rxj,SEVEN) + 

              cl_float(513.0,precision)*Power(rxj,EIGHT) + cl_float(17.0,precision)*Power(rxj,NINE)) - 

           Power(rxi,cl_float(20.0,precision))*Power(rxj,TWO)*

            (cl_float(637875.0,precision) + cl_float(1148175.0,precision)*rxj + cl_float(1020600.0,precision)*Power(rxj,TWO) + 

              cl_float(595350.0,precision)*Power(rxj,THREE) + cl_float(255150.0,precision)*Power(rxj,FOUR) + cl_float(85050.0,precision)*Power(rxj,FIVE) + 

              cl_float(22680.0,precision)*Power(rxj,SIX) + cl_float(4860.0,precision)*Power(rxj,SEVEN) + cl_float(810.0,precision)*Power(rxj,EIGHT) + 

              cl_float(34.0,precision)*Power(rxj,NINE)) + THREE*Power(rxi,EIGHT)*Power(rxj,cl_float(14.0,precision))*

            (-cl_float(593408025.0,precision) + cl_float(946053675.0,precision)*rxj - cl_float(394427880.0,precision)*Power(rxj,TWO) - 

              cl_float(315870660.0,precision)*Power(rxj,THREE) - cl_float(53891460.0,precision)*Power(rxj,FOUR) + 

              cl_float(910980.0,precision)*Power(rxj,FIVE) + cl_float(1409520.0,precision)*Power(rxj,SIX) + cl_float(192168.0,precision)*Power(rxj,SEVEN) + 

              cl_float(11196.0,precision)*Power(rxj,EIGHT) + cl_float(236.0,precision)*Power(rxj,NINE))))/

      (cl_float(42525.0,precision)*exp(TWO*(rxi + rxj))*Power(rxi - rxj,cl_float(15.0,precision))*Power(rxi + rxj,cl_float(15.0,precision)))

     ; }
   
  }
  return S;
}

cl_F Slater_3S_1S(cl_F r,cl_F xi,cl_F xj)
{
  return Slater_1S_3S(r,xj,xi);
}

cl_F Slater_3S_2S(cl_F r,cl_F xi,cl_F xj)
{
  return Slater_2S_3S(r,xj,xi);
}

cl_F Slater_4S_4S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = (cl_float(26333.0,precision)*xi)/cl_float(131072.0,precision)

    ; } else {  S = (ONE/r)*(-cl_float(83691159552000.0,precision) + cl_float(83691159552000.0,precision)*exp(TWO*r*xi) - cl_float(150568359566625.0,precision)*r*xi - 

        cl_float(133754400029250.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

        cl_float(78142908343500.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

        cl_float(33740723016000.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - 

        cl_float(11470756096800.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

        cl_float(3193358968800.0,precision)*Power(r,SIX)*Power(xi,SIX) - 

        cl_float(747112766400.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

        cl_float(149448499200.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - 

        cl_float(25830604800.0,precision)*Power(r,NINE)*Power(xi,NINE) - 

        cl_float(3874590720.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - 

        cl_float(503193600.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)) - cl_float(55910400.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(12.0,precision)) - 

        cl_float(5160960.0,precision)*Power(r,cl_float(13.0,precision))*Power(xi,cl_float(13.0,precision)) - cl_float(368640.0,precision)*Power(r,cl_float(14.0,precision))*Power(xi,cl_float(14.0,precision)) - 

        cl_float(16384.0,precision)*Power(r,cl_float(15.0,precision))*Power(xi,cl_float(15.0,precision)))/(cl_float(8.3691159552e13,precision)*exp(TWO*r*xi))

    ; }
 
  }
  else {
      if (r == ZERO) { S = (xi*xj*(Power(xi,cl_float(14.0,precision)) + cl_float(15.0,precision)*Power(xi,cl_float(13.0,precision))*xj + cl_float(105.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,TWO) + 

          cl_float(455.0,precision)*Power(xi,cl_float(11.0,precision))*Power(xj,THREE) + cl_float(1365.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,FOUR) + 

          cl_float(3003.0,precision)*Power(xi,NINE)*Power(xj,FIVE) + cl_float(5005.0,precision)*Power(xi,EIGHT)*Power(xj,SIX) + 

          cl_float(6435.0,precision)*Power(xi,SEVEN)*Power(xj,SEVEN) + cl_float(5005.0,precision)*Power(xi,SIX)*Power(xj,EIGHT) + 

          cl_float(3003.0,precision)*Power(xi,FIVE)*Power(xj,NINE) + cl_float(1365.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(10.0,precision)) + 

          cl_float(455.0,precision)*Power(xi,THREE)*Power(xj,cl_float(11.0,precision)) + cl_float(105.0,precision)*Power(xi,TWO)*Power(xj,cl_float(12.0,precision)) + 

          cl_float(15.0,precision)*xi*Power(xj,cl_float(13.0,precision)) + Power(xj,cl_float(14.0,precision))))/(cl_float(4.0,precision)*Power(xi + xj,cl_float(15.0,precision)))

    ; } else { S = (ONE/r)*(cl_float(1260.0,precision)*exp(TWO*(rxi + rxj))*Power(Power(rxi,TWO) - Power(rxj,TWO),cl_float(15.0,precision)) + 

        exp(TWO*rxj)*Power(rxj,cl_float(10.0,precision))*

         (-cl_float(3276.0,precision)*Power(rxi,cl_float(25.0,precision)) - cl_float(168.0,precision)*Power(rxi,cl_float(26.0,precision)) - FOUR*Power(rxi,cl_float(27.0,precision)) + 

           cl_float(1260.0,precision)*Power(rxj,cl_float(20.0,precision)) + cl_float(2205.0,precision)*rxi*Power(rxj,cl_float(20.0,precision)) + 

           cl_float(1890.0,precision)*Power(rxi,TWO)*Power(rxj,cl_float(18.0,precision))*(-cl_float(10.0,precision) + Power(rxj,TWO)) - 

           cl_float(420.0,precision)*Power(rxi,cl_float(24.0,precision))*(cl_float(91.0,precision) + Power(rxj,TWO)) + 

           cl_float(525.0,precision)*Power(rxi,THREE)*Power(rxj,cl_float(18.0,precision))*(-cl_float(63.0,precision) + TWO*Power(rxj,TWO)) + 

           cl_float(42.0,precision)*Power(rxi,cl_float(23.0,precision))*(-cl_float(6825.0,precision) - cl_float(405.0,precision)*Power(rxj,TWO) + TWO*Power(rxj,FOUR)) + 

           cl_float(63.0,precision)*Power(rxi,FIVE)*Power(rxj,cl_float(16.0,precision))*

            (cl_float(3675.0,precision) - cl_float(250.0,precision)*Power(rxj,TWO) + TWO*Power(rxj,FOUR)) + 

           cl_float(210.0,precision)*Power(rxi,FOUR)*Power(rxj,cl_float(16.0,precision))*

            (cl_float(630.0,precision) - cl_float(135.0,precision)*Power(rxj,TWO) + TWO*Power(rxj,FOUR)) + 

           cl_float(252.0,precision)*Power(rxi,cl_float(22.0,precision))*(-cl_float(5460.0,precision) - cl_float(1225.0,precision)*Power(rxj,TWO) + cl_float(17.0,precision)*Power(rxj,FOUR)) - 

           cl_float(1260.0,precision)*Power(rxi,cl_float(17.0,precision))*Power(rxj,FOUR)*

            (cl_float(141729.0,precision) - cl_float(10145.0,precision)*Power(rxj,TWO) + cl_float(116.0,precision)*Power(rxj,FOUR)) + 

           cl_float(21.0,precision)*Power(rxi,NINE)*Power(rxj,cl_float(12.0,precision))*

            (cl_float(164775.0,precision) - cl_float(18460.0,precision)*Power(rxj,TWO) + cl_float(828.0,precision)*Power(rxj,FOUR)) + 

           cl_float(14.0,precision)*Power(rxi,SIX)*Power(rxj,cl_float(14.0,precision))*

            (-cl_float(40950.0,precision) + cl_float(14175.0,precision)*Power(rxj,TWO) - cl_float(450.0,precision)*Power(rxj,FOUR) + TWO*Power(rxj,SIX)) - 

           cl_float(210.0,precision)*Power(rxi,EIGHT)*Power(rxj,cl_float(12.0,precision))*

            (-cl_float(8190.0,precision) + cl_float(4095.0,precision)*Power(rxj,TWO) - cl_float(210.0,precision)*Power(rxj,FOUR) + TWO*Power(rxj,SIX)) + 

           cl_float(42.0,precision)*Power(rxi,cl_float(10.0,precision))*Power(rxj,cl_float(10.0,precision))*

            (-cl_float(209430.0,precision) - cl_float(2925.0,precision)*Power(rxj,TWO) - cl_float(8840.0,precision)*Power(rxj,FOUR) + FOUR*Power(rxj,SIX)) 

    + Power(rxi,SEVEN)*Power(rxj,cl_float(14.0,precision))*(-cl_float(1003275.0,precision) + cl_float(110250.0,precision)*Power(rxj,TWO) - 

              cl_float(1890.0,precision)*Power(rxj,FOUR) + FOUR*Power(rxj,SIX)) - 

           cl_float(21.0,precision)*Power(rxi,cl_float(11.0,precision))*Power(rxj,cl_float(10.0,precision))*

            (-cl_float(1033695.0,precision) - cl_float(218400.0,precision)*Power(rxj,TWO) + cl_float(552.0,precision)*Power(rxj,FOUR) + 

              FOUR*Power(rxj,SIX)) + cl_float(280.0,precision)*Power(rxi,cl_float(18.0,precision))*Power(rxj,TWO)*

            (-cl_float(385560.0,precision) - cl_float(73953.0,precision)*Power(rxj,TWO) + cl_float(2370.0,precision)*Power(rxj,FOUR) + FOUR*Power(rxj,SIX)) 

    - cl_float(35.0,precision)*Power(rxi,cl_float(15.0,precision))*Power(rxj,SIX)*

            (-cl_float(1565613.0,precision) + cl_float(359520.0,precision)*Power(rxj,TWO) - cl_float(7020.0,precision)*Power(rxj,FOUR) + 

              EIGHT*Power(rxj,SIX)) + cl_float(14.0,precision)*Power(rxi,cl_float(19.0,precision))*Power(rxj,TWO)*

            (-cl_float(4980150.0,precision) + cl_float(126765.0,precision)*Power(rxj,TWO) - cl_float(3852.0,precision)*Power(rxj,FOUR) + 

              cl_float(20.0,precision)*Power(rxj,SIX)) - cl_float(630.0,precision)*Power(rxi,cl_float(14.0,precision))*Power(rxj,SIX)*

            (cl_float(708714.0,precision) - cl_float(14385.0,precision)*Power(rxj,TWO) - cl_float(2340.0,precision)*Power(rxj,FOUR) + cl_float(20.0,precision)*Power(rxj,SIX)) 

    + cl_float(210.0,precision)*Power(rxi,cl_float(16.0,precision))*Power(rxj,FOUR)*

            (-cl_float(2087532.0,precision) + cl_float(328491.0,precision)*Power(rxj,TWO) - cl_float(11740.0,precision)*Power(rxj,FOUR) + 

              cl_float(52.0,precision)*Power(rxj,SIX)) - cl_float(84.0,precision)*Power(rxi,cl_float(20.0,precision))*

            (cl_float(59670.0,precision) + cl_float(236250.0,precision)*Power(rxj,TWO) - cl_float(8745.0,precision)*Power(rxj,FOUR) + cl_float(92.0,precision)*Power(rxj,SIX)) 

    - TWO*Power(rxi,cl_float(21.0,precision))*(cl_float(1949220.0,precision) + cl_float(1598625.0,precision)*Power(rxj,TWO) - cl_float(41391.0,precision)*Power(rxj,FOUR) + 

              cl_float(128.0,precision)*Power(rxj,SIX)) + Power(rxi,cl_float(13.0,precision))*Power(rxj,EIGHT)*

            (cl_float(173037375.0,precision) - cl_float(2784600.0,precision)*Power(rxj,TWO) - cl_float(112140.0,precision)*Power(rxj,FOUR) + 

              cl_float(256.0,precision)*Power(rxj,SIX)) + cl_float(14.0,precision)*Power(rxi,cl_float(12.0,precision))*Power(rxj,EIGHT)*

            (-cl_float(7260750.0,precision) - cl_float(2521935.0,precision)*Power(rxj,TWO) + cl_float(19500.0,precision)*Power(rxj,FOUR) + 

              cl_float(344.0,precision)*Power(rxj,SIX))) + 

        exp(TWO*rxi)*Power(rxi,cl_float(10.0,precision))*

         (cl_float(210.0,precision)*Power(rxi,TWO)*Power(rxj,cl_float(18.0,precision))*

            (cl_float(514080.0,precision) + cl_float(332010.0,precision)*rxj + cl_float(94500.0,precision)*Power(rxj,TWO) + cl_float(15225.0,precision)*Power(rxj,THREE) + 

              cl_float(1470.0,precision)*Power(rxj,FOUR) + cl_float(81.0,precision)*Power(rxj,FIVE) + TWO*Power(rxj,SIX)) + 

           cl_float(105.0,precision)*Power(rxi,cl_float(18.0,precision))*Power(rxj,TWO)*

            (cl_float(180.0,precision) + cl_float(315.0,precision)*rxj + cl_float(270.0,precision)*Power(rxj,TWO) + cl_float(150.0,precision)*Power(rxj,THREE) + 

              cl_float(60.0,precision)*Power(rxj,FOUR) + cl_float(18.0,precision)*Power(rxj,FIVE) + FOUR*Power(rxj,SIX)) - 

           cl_float(1365.0,precision)*Power(rxi,cl_float(10.0,precision))*Power(rxj,cl_float(10.0,precision))*

            (-cl_float(6444.0,precision) + cl_float(15903.0,precision)*rxj - cl_float(25866.0,precision)*Power(rxj,TWO) - cl_float(2040.0,precision)*Power(rxj,THREE) + 

              cl_float(1080.0,precision)*Power(rxj,FOUR) + cl_float(180.0,precision)*Power(rxj,FIVE) + EIGHT*Power(rxj,SIX)) + 

           Power(rxi,cl_float(14.0,precision))*Power(rxj,SIX)*

            (cl_float(573300.0,precision) + cl_float(1003275.0,precision)*rxj + cl_float(859950.0,precision)*Power(rxj,TWO) + cl_float(387660.0,precision)*Power(rxj,THREE) + 

              cl_float(371280.0,precision)*Power(rxj,FOUR) + cl_float(11592.0,precision)*Power(rxj,FIVE) - cl_float(4816.0,precision)*Power(rxj,SIX) - 

              cl_float(256.0,precision)*Power(rxj,SEVEN)) + TWO*Power(rxj,cl_float(20.0,precision))*

            (cl_float(2506140.0,precision) + cl_float(1949220.0,precision)*rxj + cl_float(687960.0,precision)*Power(rxj,TWO) + 

              cl_float(143325.0,precision)*Power(rxj,THREE) + cl_float(19110.0,precision)*Power(rxj,FOUR) + cl_float(1638.0,precision)*Power(rxj,FIVE) + 

              cl_float(84.0,precision)*Power(rxj,SIX) + TWO*Power(rxj,SEVEN)) - 

           cl_float(42.0,precision)*Power(rxi,FOUR)*Power(rxj,cl_float(16.0,precision))*

            (-cl_float(10437660.0,precision) - cl_float(4251870.0,precision)*rxj - cl_float(493020.0,precision)*Power(rxj,TWO) + 

              cl_float(42255.0,precision)*Power(rxj,THREE) + cl_float(17490.0,precision)*Power(rxj,FOUR) + cl_float(1971.0,precision)*Power(rxj,FIVE) + 

              cl_float(102.0,precision)*Power(rxj,SIX) + TWO*Power(rxj,SEVEN)) + 

           cl_float(21.0,precision)*Power(rxi,cl_float(16.0,precision))*Power(rxj,FOUR)*

            (-cl_float(6300.0,precision) - cl_float(11025.0,precision)*rxj - cl_float(9450.0,precision)*Power(rxj,TWO) - cl_float(5250.0,precision)*Power(rxj,THREE) - 

              cl_float(2100.0,precision)*Power(rxj,FOUR) - cl_float(828.0,precision)*Power(rxj,FIVE) - EIGHT*Power(rxj,SIX) + 

              FOUR*Power(rxj,SEVEN)) - Power(rxi,cl_float(20.0,precision))*

            (cl_float(1260.0,precision) + cl_float(2205.0,precision)*rxj + cl_float(1890.0,precision)*Power(rxj,TWO) + cl_float(1050.0,precision)*Power(rxj,THREE) + 

              cl_float(420.0,precision)*Power(rxj,FOUR) + cl_float(126.0,precision)*Power(rxj,FIVE) + cl_float(28.0,precision)*Power(rxj,SIX) + 

              FOUR*Power(rxj,SEVEN)) - cl_float(35.0,precision)*Power(rxi,EIGHT)*Power(rxj,cl_float(12.0,precision))*

            (-cl_float(2904300.0,precision) + cl_float(4943925.0,precision)*rxj + cl_float(258930.0,precision)*Power(rxj,TWO) - 

              cl_float(359520.0,precision)*Power(rxj,THREE) - cl_float(70440.0,precision)*Power(rxj,FOUR) - cl_float(4176.0,precision)*Power(rxj,FIVE) + 

              cl_float(32.0,precision)*Power(rxj,SIX) + EIGHT*Power(rxj,SEVEN)) + 

           cl_float(35.0,precision)*Power(rxi,cl_float(12.0,precision))*Power(rxj,EIGHT)*

            (-cl_float(49140.0,precision) - cl_float(98865.0,precision)*rxj + cl_float(3510.0,precision)*Power(rxj,TWO) - cl_float(131040.0,precision)*Power(rxj,THREE) - 

              cl_float(7800.0,precision)*Power(rxj,FOUR) + cl_float(3204.0,precision)*Power(rxj,FIVE) + cl_float(360.0,precision)*Power(rxj,SIX) + 

              EIGHT*Power(rxj,SEVEN)) + Power(rxi,SIX)*Power(rxj,cl_float(14.0,precision))*

            (cl_float(446489820.0,precision) - cl_float(54796455.0,precision)*rxj - cl_float(68983110.0,precision)*Power(rxj,TWO) - 

              cl_float(12782700.0,precision)*Power(rxj,THREE) - cl_float(663600.0,precision)*Power(rxj,FOUR) + cl_float(53928.0,precision)*Power(rxj,FIVE) + 

              cl_float(7728.0,precision)*Power(rxj,SIX) + cl_float(256.0,precision)*Power(rxj,SEVEN))))/

      (cl_float(1260.0,precision)*exp(TWO*(rxi + rxj))*Power(rxi - rxj,cl_float(15.0,precision))*Power(rxi + rxj,cl_float(15.0,precision)))

     ; }
   
  }
  return S;
}

cl_F Slater_4S_5S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = (cl_float(234137.0,precision)*xi)/cl_float(1.31072e6,precision)

    ; } else {  S = (ONE/r)*(-cl_float(14227497123840000.0,precision) + cl_float(14227497123840000.0,precision)*exp(TWO*r*xi) - 

        cl_float(25913502934444125.0,precision)*r*xi - cl_float(23372011621208250.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

        cl_float(13907709869303250.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

        cl_float(6137735659555500.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - 

        cl_float(2140857388870200.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

        cl_float(614116575072000.0,precision)*Power(r,SIX)*Power(xi,SIX) - 

        cl_float(148809580920000.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

        cl_float(31036639233600.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - 

        cl_float(5645342102400.0,precision)*Power(r,NINE)*Power(xi,NINE) - 

        cl_float(903333150720.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - 

        cl_float(127744081920.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)) - 

        cl_float(15968010240.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(12.0,precision)) - 

        cl_float(1754726400.0,precision)*Power(r,cl_float(13.0,precision))*Power(xi,cl_float(13.0,precision)) - 

        cl_float(167116800.0,precision)*Power(r,cl_float(14.0,precision))*Power(xi,cl_float(14.0,precision)) - cl_float(13369344.0,precision)*Power(r,cl_float(15.0,precision))*Power(xi,cl_float(15.0,precision)) - 

        cl_float(835584.0,precision)*Power(r,cl_float(16.0,precision))*Power(xi,cl_float(16.0,precision)) - cl_float(32768.0,precision)*Power(r,cl_float(17.0,precision))*Power(xi,cl_float(17.0,precision)))/

      (cl_float(1.422749712384e16,precision)*exp(TWO*r*xi))

    ; }
 
  }
  else {
      if (r == ZERO) { S = (xi*xj*(FOUR*Power(xi,cl_float(16.0,precision)) + cl_float(68.0,precision)*Power(xi,cl_float(15.0,precision))*xj + cl_float(544.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,TWO) + 

          cl_float(2720.0,precision)*Power(xi,cl_float(13.0,precision))*Power(xj,THREE) + cl_float(9520.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,FOUR) + 

          cl_float(24752.0,precision)*Power(xi,cl_float(11.0,precision))*Power(xj,FIVE) + cl_float(49504.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,SIX) + 

          cl_float(77792.0,precision)*Power(xi,NINE)*Power(xj,SEVEN) + cl_float(97240.0,precision)*Power(xi,EIGHT)*Power(xj,EIGHT) + 

          cl_float(97240.0,precision)*Power(xi,SEVEN)*Power(xj,NINE) + cl_float(61880.0,precision)*Power(xi,SIX)*Power(xj,cl_float(10.0,precision)) + 

          cl_float(30940.0,precision)*Power(xi,FIVE)*Power(xj,cl_float(11.0,precision)) + cl_float(11900.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(12.0,precision)) + 

          cl_float(3400.0,precision)*Power(xi,THREE)*Power(xj,cl_float(13.0,precision)) + cl_float(680.0,precision)*Power(xi,TWO)*Power(xj,cl_float(14.0,precision)) + 

          cl_float(85.0,precision)*xi*Power(xj,cl_float(15.0,precision)) + FIVE*Power(xj,cl_float(16.0,precision))))/(cl_float(20.0,precision)*Power(xi + xj,cl_float(17.0,precision)))

    ; } else { S = (ONE/r)*(cl_float(56700.0,precision)*exp(TWO*(rxi + rxj))*Power(Power(rxi,TWO) - Power(rxj,TWO),cl_float(17.0,precision)) + 

        NINE*exp(TWO*rxj)*Power(rxj,cl_float(12.0,precision))*

         (-cl_float(980.0,precision)*Power(rxi,cl_float(28.0,precision)) - cl_float(20.0,precision)*Power(rxi,cl_float(29.0,precision)) + cl_float(6300.0,precision)*Power(rxj,cl_float(22.0,precision)) + 

           cl_float(11025.0,precision)*rxi*Power(rxj,cl_float(22.0,precision)) - cl_float(50.0,precision)*Power(rxi,cl_float(27.0,precision))*(cl_float(441.0,precision) + TWO*Power(rxj,TWO)) + 

           cl_float(3150.0,precision)*Power(rxi,TWO)*Power(rxj,cl_float(20.0,precision))*(-cl_float(34.0,precision) + THREE*Power(rxj,TWO)) + 

           cl_float(525.0,precision)*Power(rxi,THREE)*Power(rxj,cl_float(20.0,precision))*(-cl_float(357.0,precision) + cl_float(10.0,precision)*Power(rxj,TWO)) - 

           cl_float(420.0,precision)*Power(rxi,cl_float(26.0,precision))*(cl_float(700.0,precision) + cl_float(19.0,precision)*Power(rxj,TWO)) + 

           cl_float(1050.0,precision)*Power(rxi,FOUR)*Power(rxj,cl_float(18.0,precision))*

            (cl_float(816.0,precision) - cl_float(153.0,precision)*Power(rxj,TWO) + TWO*Power(rxj,FOUR)) + 

           cl_float(210.0,precision)*Power(rxi,FIVE)*Power(rxj,cl_float(18.0,precision))*

            (cl_float(7140.0,precision) - cl_float(425.0,precision)*Power(rxj,TWO) + THREE*Power(rxj,FOUR)) + 

           cl_float(42.0,precision)*Power(rxi,cl_float(25.0,precision))*(-cl_float(59500.0,precision) - cl_float(6035.0,precision)*Power(rxj,TWO) + cl_float(18.0,precision)*Power(rxj,FOUR)) + 

           cl_float(84.0,precision)*Power(rxi,cl_float(24.0,precision))*(-cl_float(160650.0,precision) - cl_float(52700.0,precision)*Power(rxj,TWO) + cl_float(397.0,precision)*Power(rxj,FOUR)) - 

           cl_float(28.0,precision)*Power(rxi,cl_float(12.0,precision))*Power(rxj,cl_float(10.0,precision))*

            (cl_float(100849950.0,precision) + cl_float(27100125.0,precision)*Power(rxj,TWO) + cl_float(186150.0,precision)*Power(rxj,FOUR) - 

              cl_float(2177.0,precision)*Power(rxj,SIX)) + 

           cl_float(140.0,precision)*Power(rxi,SIX)*Power(rxj,cl_float(16.0,precision))*

            (-cl_float(30600.0,precision) + cl_float(9180.0,precision)*Power(rxj,TWO) - cl_float(255.0,precision)*Power(rxj,FOUR) + Power(rxj,SIX)) - 

           cl_float(2380.0,precision)*Power(rxi,EIGHT)*Power(rxj,cl_float(14.0,precision))*

            (-cl_float(6300.0,precision) + cl_float(2700.0,precision)*Power(rxj,TWO) - cl_float(120.0,precision)*Power(rxj,FOUR) + Power(rxj,SIX)) + 

           cl_float(10.0,precision)*Power(rxi,SEVEN)*Power(rxj,cl_float(16.0,precision))*

            (-cl_float(749700.0,precision) + cl_float(71400.0,precision)*Power(rxj,TWO) - cl_float(1071.0,precision)*Power(rxj,FOUR) + TWO*Power(rxj,SIX)) 

    + cl_float(204.0,precision)*Power(rxi,cl_float(15.0,precision))*Power(rxj,EIGHT)*

            (cl_float(28962255.0,precision) - cl_float(1744750.0,precision)*Power(rxj,TWO) + cl_float(9555.0,precision)*Power(rxj,FOUR) + 

              SIX*Power(rxj,SIX)) - cl_float(42.0,precision)*Power(rxi,cl_float(11.0,precision))*Power(rxj,cl_float(12.0,precision))*

            (-cl_float(12911925.0,precision) - cl_float(1634550.0,precision)*Power(rxj,TWO) - cl_float(7103.0,precision)*Power(rxj,FOUR) + 

              cl_float(18.0,precision)*Power(rxj,SIX)) + TWO*Power(rxi,NINE)*Power(rxj,cl_float(14.0,precision))*

            (cl_float(16948575.0,precision) - cl_float(1184400.0,precision)*Power(rxj,TWO) + cl_float(63861.0,precision)*Power(rxj,FOUR) + 

              cl_float(50.0,precision)*Power(rxj,SIX)) + cl_float(28.0,precision)*Power(rxi,cl_float(22.0,precision))*

            (-cl_float(2180250.0,precision) - cl_float(10993050.0,precision)*Power(rxj,TWO) + cl_float(14925.0,precision)*Power(rxj,FOUR) + 

              cl_float(73.0,precision)*Power(rxj,SIX)) - cl_float(952.0,precision)*Power(rxi,cl_float(14.0,precision))*Power(rxj,EIGHT)*

            (cl_float(16966215.0,precision) + cl_float(725175.0,precision)*Power(rxj,TWO) - cl_float(36075.0,precision)*Power(rxj,FOUR) + 

              cl_float(79.0,precision)*Power(rxj,SIX)) - cl_float(84.0,precision)*Power(rxi,cl_float(10.0,precision))*Power(rxj,cl_float(12.0,precision))*

            (cl_float(1723800.0,precision) + cl_float(279225.0,precision)*Power(rxj,TWO) + cl_float(45600.0,precision)*Power(rxj,FOUR) + 

              cl_float(107.0,precision)*Power(rxj,SIX)) - cl_float(35.0,precision)*Power(rxi,cl_float(17.0,precision))*Power(rxj,SIX)*

            (cl_float(132637869.0,precision) - cl_float(2205240.0,precision)*Power(rxj,TWO) - cl_float(48348.0,precision)*Power(rxj,FOUR) + 

              cl_float(136.0,precision)*Power(rxj,SIX)) - SIX*Power(rxi,cl_float(21.0,precision))*Power(rxj,TWO)*

            (cl_float(192298050.0,precision) + cl_float(12644275.0,precision)*Power(rxj,TWO) - cl_float(218029.0,precision)*Power(rxj,FOUR) + 

              cl_float(204.0,precision)*Power(rxj,SIX)) + FOUR*Power(rxi,cl_float(13.0,precision))*Power(rxj,cl_float(10.0,precision))*

            (cl_float(1259522775.0,precision) + cl_float(15895425.0,precision)*Power(rxj,TWO) - cl_float(493017.0,precision)*Power(rxj,FOUR) + 

              cl_float(263.0,precision)*Power(rxj,SIX)) - cl_float(140.0,precision)*Power(rxi,cl_float(16.0,precision))*Power(rxj,SIX)*

            (cl_float(180826281.0,precision) - cl_float(15101406.0,precision)*Power(rxj,TWO) + cl_float(160140.0,precision)*Power(rxj,FOUR) + 

              cl_float(442.0,precision)*Power(rxj,SIX)) - TWO*Power(rxi,cl_float(23.0,precision))*

            (cl_float(21366450.0,precision) + cl_float(23526300.0,precision)*Power(rxj,TWO) - cl_float(246729.0,precision)*Power(rxj,FOUR) + 

              cl_float(526.0,precision)*Power(rxj,SIX)) + SEVEN*Power(rxi,cl_float(19.0,precision))*Power(rxj,FOUR)*

            (-cl_float(811081215.0,precision) + cl_float(39095550.0,precision)*Power(rxj,TWO) - cl_float(515916.0,precision)*Power(rxj,FOUR) + 

              cl_float(680.0,precision)*Power(rxj,SIX)) + cl_float(70.0,precision)*Power(rxi,cl_float(18.0,precision))*Power(rxj,FOUR)*

            (-cl_float(180554454.0,precision) + cl_float(9873711.0,precision)*Power(rxj,TWO) - cl_float(414120.0,precision)*Power(rxj,FOUR) + 

              cl_float(2924.0,precision)*Power(rxj,SIX)) - 

           cl_float(14.0,precision)*Power(rxi,cl_float(20.0,precision))*Power(rxj,TWO)*

            (cl_float(136919700.0,precision) + cl_float(71867115.0,precision)*Power(rxj,TWO) - cl_float(2154150.0,precision)*Power(rxj,FOUR) + 

              cl_float(10268.0,precision)*Power(rxj,SIX))) - 

        FOUR*exp(TWO*rxi)*Power(rxi,cl_float(10.0,precision))*

         (-cl_float(10710.0,precision)*Power(rxi,cl_float(12.0,precision))*Power(rxj,cl_float(12.0,precision))*

            (-cl_float(3555.0,precision) - cl_float(127008.0,precision)*rxj + cl_float(138384.0,precision)*Power(rxj,TWO) - cl_float(74556.0,precision)*Power(rxj,THREE) - 

              cl_float(22284.0,precision)*Power(rxj,FOUR) + cl_float(408.0,precision)*Power(rxj,FIVE) + cl_float(576.0,precision)*Power(rxj,SIX) + 

              cl_float(60.0,precision)*Power(rxj,SEVEN) + TWO*Power(rxj,EIGHT)) + 

           TWO*Power(rxi,cl_float(20.0,precision))*Power(rxj,FOUR)*

            (cl_float(963900.0,precision) + cl_float(1735020.0,precision)*rxj + cl_float(1542240.0,precision)*Power(rxj,TWO) + 

              cl_float(899640.0,precision)*Power(rxj,THREE) + cl_float(385560.0,precision)*Power(rxj,FOUR) + cl_float(128520.0,precision)*Power(rxj,FIVE) + 

              cl_float(34272.0,precision)*Power(rxj,SIX) + cl_float(9126.0,precision)*Power(rxj,SEVEN) + cl_float(333.0,precision)*Power(rxj,EIGHT) - 

              cl_float(20.0,precision)*Power(rxj,NINE)) - TWO*Power(rxj,cl_float(24.0,precision))*

            (cl_float(119041650.0,precision) + cl_float(107137485.0,precision)*rxj + cl_float(45110520.0,precision)*Power(rxj,TWO) + 

              cl_float(11695320.0,precision)*Power(rxj,THREE) + cl_float(2063880.0,precision)*Power(rxj,FOUR) + 

              cl_float(257985.0,precision)*Power(rxj,FIVE) + cl_float(22932.0,precision)*Power(rxj,SIX) + cl_float(1404.0,precision)*Power(rxj,SEVEN) + 

              cl_float(54.0,precision)*Power(rxj,EIGHT) + Power(rxj,NINE)) + 

           TWO*Power(rxi,TWO)*Power(rxj,cl_float(22.0,precision))*

            (-cl_float(3264488325.0,precision) - cl_float(2505368880.0,precision)*rxj - cl_float(881390160.0,precision)*Power(rxj,TWO) - 

              cl_float(185775660.0,precision)*Power(rxj,THREE) - cl_float(25639740.0,precision)*Power(rxj,FOUR) - 

              cl_float(2361555.0,precision)*Power(rxj,FIVE) - cl_float(139356.0,precision)*Power(rxj,SIX) - cl_float(4482.0,precision)*Power(rxj,SEVEN) - 

              cl_float(27.0,precision)*Power(rxj,EIGHT) + TWO*Power(rxj,NINE)) + 

           Power(rxi,cl_float(24.0,precision))*(cl_float(14175.0,precision) + cl_float(25515.0,precision)*rxj + cl_float(22680.0,precision)*Power(rxj,TWO) + 

              cl_float(13230.0,precision)*Power(rxj,THREE) + cl_float(5670.0,precision)*Power(rxj,FOUR) + cl_float(1890.0,precision)*Power(rxj,FIVE) + 

              cl_float(504.0,precision)*Power(rxj,SIX) + cl_float(108.0,precision)*Power(rxj,SEVEN) + cl_float(18.0,precision)*Power(rxj,EIGHT) + 

              TWO*Power(rxj,NINE)) - cl_float(102.0,precision)*Power(rxi,cl_float(10.0,precision))*Power(rxj,cl_float(14.0,precision))*

            (cl_float(44986725.0,precision) - cl_float(97433280.0,precision)*rxj + cl_float(44467920.0,precision)*Power(rxj,TWO) + 

              cl_float(15857100.0,precision)*Power(rxj,THREE) - cl_float(457380.0,precision)*Power(rxj,FOUR) - 

              cl_float(620550.0,precision)*Power(rxj,FIVE) - cl_float(83160.0,precision)*Power(rxj,SIX) - cl_float(4068.0,precision)*Power(rxj,SEVEN) - 

              SIX*Power(rxj,EIGHT) + FOUR*Power(rxj,NINE)) + 

           cl_float(102.0,precision)*Power(rxi,cl_float(14.0,precision))*Power(rxj,cl_float(10.0,precision))*

            (-cl_float(859950.0,precision) - cl_float(1437345.0,precision)*rxj - cl_float(2260440.0,precision)*Power(rxj,TWO) + 

              cl_float(810810.0,precision)*Power(rxj,THREE) - cl_float(1056510.0,precision)*Power(rxj,FOUR) - 

              cl_float(217854.0,precision)*Power(rxj,FIVE) + cl_float(6552.0,precision)*Power(rxj,SIX) + cl_float(3852.0,precision)*Power(rxj,SEVEN) + 

              cl_float(258.0,precision)*Power(rxj,EIGHT) + FOUR*Power(rxj,NINE)) - 

           Power(rxi,cl_float(22.0,precision))*Power(rxj,TWO)*

            (cl_float(240975.0,precision) + cl_float(433755.0,precision)*rxj + cl_float(385560.0,precision)*Power(rxj,TWO) + cl_float(224910.0,precision)*Power(rxj,THREE) + 

              cl_float(96390.0,precision)*Power(rxj,FOUR) + cl_float(32130.0,precision)*Power(rxj,FIVE) + cl_float(8568.0,precision)*Power(rxj,SIX) + 

              cl_float(1836.0,precision)*Power(rxj,SEVEN) + cl_float(306.0,precision)*Power(rxj,EIGHT) + FOUR*Power(rxj,NINE)) + 

           TWO*Power(rxi,FOUR)*Power(rxj,cl_float(20.0,precision))*

            (-cl_float(18032978565.0,precision) - cl_float(9823683240.0,precision)*rxj - cl_float(2047323600.0,precision)*Power(rxj,TWO) - 

              cl_float(129098340.0,precision)*Power(rxj,THREE) + cl_float(26410860.0,precision)*Power(rxj,FOUR) + 

              cl_float(7094304.0,precision)*Power(rxj,FIVE) + cl_float(788256.0,precision)*Power(rxj,SIX) + cl_float(48654.0,precision)*Power(rxj,SEVEN) + 

              cl_float(1593.0,precision)*Power(rxj,EIGHT) + cl_float(20.0,precision)*Power(rxj,NINE)) - 

           SIX*Power(rxi,cl_float(16.0,precision))*Power(rxj,EIGHT)*

            (-cl_float(5622750.0,precision) - cl_float(10120950.0,precision)*rxj - cl_float(8996400.0,precision)*Power(rxj,TWO) - 

              cl_float(5698350.0,precision)*Power(rxj,THREE) - cl_float(897750.0,precision)*Power(rxj,FOUR) - 

              cl_float(1641591.0,precision)*Power(rxj,FIVE) - cl_float(211932.0,precision)*Power(rxj,SIX) + cl_float(10224.0,precision)*Power(rxj,SEVEN) + 

              cl_float(2364.0,precision)*Power(rxj,EIGHT) + cl_float(73.0,precision)*Power(rxj,NINE)) + 

           TWO*Power(rxi,cl_float(18.0,precision))*Power(rxj,SIX)*

            (-cl_float(4819500.0,precision) - cl_float(8675100.0,precision)*rxj - cl_float(7711200.0,precision)*Power(rxj,TWO) - 

              cl_float(4498200.0,precision)*Power(rxj,THREE) - cl_float(1927800.0,precision)*Power(rxj,FOUR) - 

              cl_float(561519.0,precision)*Power(rxj,FIVE) - cl_float(279468.0,precision)*Power(rxj,SIX) - cl_float(20682.0,precision)*Power(rxj,SEVEN) + 

              cl_float(1305.0,precision)*Power(rxj,EIGHT) + cl_float(106.0,precision)*Power(rxj,NINE)) + 

           THREE*Power(rxi,EIGHT)*Power(rxj,cl_float(16.0,precision))*

            (-cl_float(9364244085.0,precision) + cl_float(6940428705.0,precision)*rxj + cl_float(2117684520.0,precision)*Power(rxj,TWO) - 

              cl_float(230268150.0,precision)*Power(rxj,THREE) - cl_float(149610510.0,precision)*Power(rxj,FOUR) - 

              cl_float(21824334.0,precision)*Power(rxj,FIVE) - cl_float(1223208.0,precision)*Power(rxj,SIX) + 

              cl_float(12708.0,precision)*Power(rxj,SEVEN) + cl_float(4470.0,precision)*Power(rxj,EIGHT) + cl_float(146.0,precision)*Power(rxj,NINE)) - 

           Power(rxi,SIX)*Power(rxj,cl_float(18.0,precision))*

            (cl_float(57304872765.0,precision) + cl_float(7147185255.0,precision)*rxj - cl_float(5801702760.0,precision)*Power(rxj,TWO) - 

              cl_float(2053388610.0,precision)*Power(rxj,THREE) - cl_float(271655370.0,precision)*Power(rxj,FOUR) - 

              cl_float(10864854.0,precision)*Power(rxj,FIVE) + cl_float(1337112.0,precision)*Power(rxj,SIX) + 

              cl_float(202716.0,precision)*Power(rxj,SEVEN) + cl_float(10746.0,precision)*Power(rxj,EIGHT) + cl_float(212.0,precision)*Power(rxj,NINE))))/

      (cl_float(56700.0,precision)*exp(TWO*(rxi + rxj))*Power(rxi - rxj,cl_float(17.0,precision))*Power(rxi + rxj,cl_float(17.0,precision)))

     ; }
   
  }
  return S;
}

cl_F Slater_4S_1S(cl_F r,cl_F xi,cl_F xj)
{
  return Slater_1S_4S(r,xj,xi);
}

cl_F Slater_4S_2S(cl_F r,cl_F xi,cl_F xj)
{
  return Slater_2S_4S(r,xj,xi);
}

cl_F Slater_4S_3S(cl_F r,cl_F xi,cl_F xj)
{
  return Slater_3S_4S(r,xj,xi);
}

cl_F Slater_5S_5S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = (cl_float(43191.0,precision)*xi)/cl_float(262144.0,precision)

    ; } else {  S = (ONE/r)*(-cl_float(12164510040883200000.0,precision) + cl_float(12164510040883200000.0,precision)*exp(TWO*r*xi) - 

        cl_float(22324788235240115625.0,precision)*r*xi - 

        cl_float(20320556388713831250.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

        cl_float(12225924086428552500.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

        cl_float(5467446348494130000.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - 

        cl_float(1937619942864606000.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

        cl_float(566528792821992000.0,precision)*Power(r,SIX)*Power(xi,SIX) - 

        cl_float(140462831126217600.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

        cl_float(30115609927603200.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - 

        cl_float(5663731244371200.0,precision)*Power(r,NINE)*Power(xi,NINE) - 

        cl_float(943983142502400.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - 

        cl_float(140427244339200.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)) - 

        cl_float(18723632578560.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(12.0,precision)) - 

        cl_float(2240434667520.0,precision)*Power(r,cl_float(13.0,precision))*Power(xi,cl_float(13.0,precision)) - 

        cl_float(240046571520.0,precision)*Power(r,cl_float(14.0,precision))*Power(xi,cl_float(14.0,precision)) - 

        cl_float(22861578240.0,precision)*Power(r,cl_float(15.0,precision))*Power(xi,cl_float(15.0,precision)) - 

        cl_float(1905131520.0,precision)*Power(r,cl_float(16.0,precision))*Power(xi,cl_float(16.0,precision)) - 

        cl_float(134479872.0,precision)*Power(r,cl_float(17.0,precision))*Power(xi,cl_float(17.0,precision)) - cl_float(7471104.0,precision)*Power(r,cl_float(18.0,precision))*Power(xi,cl_float(18.0,precision)) - 

        cl_float(262144.0,precision)*Power(r,cl_float(19.0,precision))*Power(xi,cl_float(19.0,precision)))/(cl_float(1.21645100408832e19,precision)*exp(TWO*r*xi))

    ; }
 
  }
  else {
      if (r == ZERO) { S = (xi*xj*(Power(xi,cl_float(18.0,precision)) + cl_float(19.0,precision)*Power(xi,cl_float(17.0,precision))*xj + cl_float(171.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,TWO) + 

          cl_float(969.0,precision)*Power(xi,cl_float(15.0,precision))*Power(xj,THREE) + cl_float(3876.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,FOUR) + 

          cl_float(11628.0,precision)*Power(xi,cl_float(13.0,precision))*Power(xj,FIVE) + cl_float(27132.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,SIX) + 

          cl_float(50388.0,precision)*Power(xi,cl_float(11.0,precision))*Power(xj,SEVEN) + cl_float(75582.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,EIGHT) + 

          cl_float(92378.0,precision)*Power(xi,NINE)*Power(xj,NINE) + cl_float(75582.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(10.0,precision)) + 

          cl_float(50388.0,precision)*Power(xi,SEVEN)*Power(xj,cl_float(11.0,precision)) + cl_float(27132.0,precision)*Power(xi,SIX)*Power(xj,cl_float(12.0,precision)) + 

          cl_float(11628.0,precision)*Power(xi,FIVE)*Power(xj,cl_float(13.0,precision)) + cl_float(3876.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(14.0,precision)) + 

          cl_float(969.0,precision)*Power(xi,THREE)*Power(xj,cl_float(15.0,precision)) + cl_float(171.0,precision)*Power(xi,TWO)*Power(xj,cl_float(16.0,precision)) + 

          cl_float(19.0,precision)*xi*Power(xj,cl_float(17.0,precision)) + Power(xj,cl_float(18.0,precision))))/(cl_float(5.0,precision)*Power(xi + xj,cl_float(19.0,precision)))

    ; } else { S = (ONE/r)*(cl_float(70875.0,precision)*exp(TWO*(rxi + rxj))*Power(Power(rxi,TWO) - Power(rxj,TWO),cl_float(19.0,precision)) + 

        exp(TWO*rxj)*Power(rxj,cl_float(12.0,precision))*

         (-cl_float(630.0,precision)*Power(rxi,cl_float(34.0,precision)) - cl_float(10.0,precision)*Power(rxi,cl_float(35.0,precision)) + cl_float(70875.0,precision)*Power(rxj,cl_float(26.0,precision)) + 

           cl_float(127575.0,precision)*rxi*Power(rxj,cl_float(26.0,precision)) - cl_float(30.0,precision)*Power(rxi,cl_float(33.0,precision))*(cl_float(630.0,precision) + Power(rxj,TWO)) + 

           cl_float(14175.0,precision)*Power(rxi,TWO)*Power(rxj,cl_float(24.0,precision))*(-cl_float(95.0,precision) + EIGHT*Power(rxj,TWO)) + 

           cl_float(4725.0,precision)*Power(rxi,THREE)*Power(rxj,cl_float(24.0,precision))*(-cl_float(513.0,precision) + cl_float(14.0,precision)*Power(rxj,TWO)) - 

           cl_float(90.0,precision)*Power(rxi,cl_float(32.0,precision))*(cl_float(3920.0,precision) + cl_float(43.0,precision)*Power(rxj,TWO)) + 

           cl_float(4725.0,precision)*Power(rxi,FIVE)*Power(rxj,cl_float(22.0,precision))*

            (cl_float(4617.0,precision) - cl_float(266.0,precision)*Power(rxj,TWO) + TWO*Power(rxj,FOUR)) + 

           cl_float(14175.0,precision)*Power(rxi,FOUR)*Power(rxj,cl_float(22.0,precision))*

            (cl_float(855.0,precision) - cl_float(152.0,precision)*Power(rxj,TWO) + TWO*Power(rxj,FOUR)) + 

           cl_float(36.0,precision)*Power(rxi,cl_float(31.0,precision))*(-cl_float(124950.0,precision) - cl_float(4985.0,precision)*Power(rxj,TWO) + cl_float(13.0,precision)*Power(rxj,FOUR)) + 

           cl_float(36.0,precision)*Power(rxi,cl_float(30.0,precision))*(-cl_float(1124550.0,precision) - cl_float(127960.0,precision)*Power(rxj,TWO) + 

              cl_float(863.0,precision)*Power(rxj,FOUR)) + cl_float(135.0,precision)*Power(rxi,SEVEN)*Power(rxj,cl_float(20.0,precision))*

            (-cl_float(915705.0,precision) + cl_float(83790.0,precision)*Power(rxj,TWO) - cl_float(1330.0,precision)*Power(rxj,FOUR) + FOUR*Power(rxj,SIX)) 

    + cl_float(315.0,precision)*Power(rxi,SIX)*Power(rxj,cl_float(20.0,precision))*

            (-cl_float(218025.0,precision) + cl_float(61560.0,precision)*Power(rxj,TWO) - cl_float(1710.0,precision)*Power(rxj,FOUR) + EIGHT*Power(rxj,SIX)) 

    - cl_float(36.0,precision)*Power(rxi,cl_float(29.0,precision))*(cl_float(7122150.0,precision) + cl_float(2102730.0,precision)*Power(rxj,TWO) - cl_float(23294.0,precision)*Power(rxj,FOUR) + 

              cl_float(37.0,precision)*Power(rxj,SIX)) - cl_float(36.0,precision)*Power(rxi,cl_float(28.0,precision))*

            (cl_float(30523500.0,precision) + cl_float(23401350.0,precision)*Power(rxj,TWO) - cl_float(299250.0,precision)*Power(rxj,FOUR) + 

              cl_float(1297.0,precision)*Power(rxj,SIX)) + 

           Power(rxi,cl_float(17.0,precision))*Power(rxj,cl_float(10.0,precision))*

            (cl_float(1073961177975.0,precision) - cl_float(21753487980.0,precision)*Power(rxj,TWO) - 

              cl_float(745994340.0,precision)*Power(rxj,FOUR) + cl_float(5307156.0,precision)*Power(rxj,SIX) - cl_float(818.0,precision)*Power(rxj,EIGHT)) 

    + cl_float(10.0,precision)*Power(rxi,NINE)*Power(rxj,cl_float(18.0,precision))*

            (cl_float(49448070.0,precision) - cl_float(6409935.0,precision)*Power(rxj,TWO) + cl_float(161595.0,precision)*Power(rxj,FOUR) - 

              cl_float(1026.0,precision)*Power(rxj,SIX) + Power(rxj,EIGHT)) + 

           cl_float(90.0,precision)*Power(rxi,EIGHT)*Power(rxj,cl_float(18.0,precision))*

            (cl_float(3052350.0,precision) - cl_float(1220940.0,precision)*Power(rxj,TWO) + cl_float(53865.0,precision)*Power(rxj,FOUR) - 

              cl_float(532.0,precision)*Power(rxj,SIX) + Power(rxj,EIGHT)) - 

           cl_float(1710.0,precision)*Power(rxi,cl_float(10.0,precision))*Power(rxj,cl_float(16.0,precision))*

            (cl_float(481950.0,precision) - cl_float(257040.0,precision)*Power(rxj,TWO) + cl_float(16065.0,precision)*Power(rxj,FOUR) - 

              cl_float(252.0,precision)*Power(rxj,SIX) + Power(rxj,EIGHT)) + 

           SIX*Power(rxi,cl_float(11.0,precision))*Power(rxj,cl_float(16.0,precision))*

            (-cl_float(207559800.0,precision) + cl_float(50390550.0,precision)*Power(rxj,TWO) - cl_float(1165815.0,precision)*Power(rxj,FOUR) + 

              cl_float(21396.0,precision)*Power(rxj,SIX) + FIVE*Power(rxj,EIGHT)) - 

           cl_float(18.0,precision)*Power(rxi,cl_float(13.0,precision))*Power(rxj,cl_float(14.0,precision))*

            (-cl_float(1703720025.0,precision) - cl_float(155669850.0,precision)*Power(rxj,TWO) - cl_float(7410270.0,precision)*Power(rxj,FOUR) - 

              cl_float(1532.0,precision)*Power(rxj,SIX) + cl_float(26.0,precision)*Power(rxj,EIGHT)) + 

           cl_float(18.0,precision)*Power(rxi,cl_float(15.0,precision))*Power(rxj,cl_float(12.0,precision))*

            (cl_float(19380896325.0,precision) + cl_float(1329128850.0,precision)*Power(rxj,TWO) - cl_float(7608930.0,precision)*Power(rxj,FOUR) - 

              cl_float(116238.0,precision)*Power(rxj,SIX) + cl_float(74.0,precision)*Power(rxj,EIGHT)) - 

           cl_float(18.0,precision)*Power(rxi,cl_float(12.0,precision))*Power(rxj,cl_float(14.0,precision))*

            (cl_float(89026875.0,precision) + cl_float(179071200.0,precision)*Power(rxj,TWO) + cl_float(1552950.0,precision)*Power(rxj,FOUR) + 

              cl_float(295820.0,precision)*Power(rxj,SIX) + cl_float(146.0,precision)*Power(rxj,EIGHT)) + 

           cl_float(18.0,precision)*Power(rxi,cl_float(25.0,precision))*Power(rxj,TWO)*

            (-cl_float(5449970925.0,precision) - cl_float(1137574935.0,precision)*Power(rxj,TWO) + cl_float(37834755.0,precision)*Power(rxj,FOUR) - 

              cl_float(273062.0,precision)*Power(rxj,SIX) + cl_float(171.0,precision)*Power(rxj,EIGHT)) - 

           NINE*Power(rxi,cl_float(19.0,precision))*Power(rxj,EIGHT)*

            (-cl_float(37914907275.0,precision) + cl_float(7613889570.0,precision)*Power(rxj,TWO) - cl_float(170524620.0,precision)*Power(rxj,FOUR) + 

              cl_float(397936.0,precision)*Power(rxj,SIX) + cl_float(342.0,precision)*Power(rxj,EIGHT)) + 

           Power(rxi,cl_float(27.0,precision))*(-cl_float(2884470750.0,precision) - cl_float(6409935000.0,precision)*Power(rxj,TWO) + 

              cl_float(28332990.0,precision)*Power(rxj,FOUR) + cl_float(58104.0,precision)*Power(rxj,SIX) + cl_float(818.0,precision)*Power(rxj,EIGHT)) - 

           THREE*Power(rxi,cl_float(23.0,precision))*Power(rxj,FOUR)*

            (cl_float(219130630425.0,precision) - cl_float(11118046590.0,precision)*Power(rxj,TWO) + 

              cl_float(327611970.0,precision)*Power(rxj,FOUR) - cl_float(2920908.0,precision)*Power(rxj,SIX) + cl_float(2584.0,precision)*Power(rxj,EIGHT)

    ) + THREE*Power(rxi,cl_float(21.0,precision))*Power(rxj,SIX)*

            (-cl_float(345162539925.0,precision) + cl_float(19030764690.0,precision)*Power(rxj,TWO) - 

              cl_float(141976170.0,precision)*Power(rxj,FOUR) - cl_float(1441872.0,precision)*Power(rxj,SIX) + cl_float(2584.0,precision)*Power(rxj,EIGHT)

    ) + cl_float(63.0,precision)*Power(rxi,cl_float(20.0,precision))*Power(rxj,SIX)*

            (-cl_float(50980542525.0,precision) + cl_float(6240202920.0,precision)*Power(rxj,TWO) - cl_float(201314310.0,precision)*Power(rxj,FOUR) + 

              cl_float(956080.0,precision)*Power(rxj,SIX) + cl_float(2584.0,precision)*Power(rxj,EIGHT)) + 

           cl_float(18.0,precision)*Power(rxi,cl_float(14.0,precision))*Power(rxj,cl_float(12.0,precision))*

            (-cl_float(7803332775.0,precision) - cl_float(2519206200.0,precision)*Power(rxj,TWO) - cl_float(119719950.0,precision)*Power(rxj,FOUR) + 

              cl_float(182280.0,precision)*Power(rxj,SIX) + cl_float(2734.0,precision)*Power(rxj,EIGHT)) - 

           cl_float(18.0,precision)*Power(rxi,cl_float(26.0,precision))*(cl_float(195859125.0,precision) + cl_float(1794781800.0,precision)*Power(rxj,TWO) + 

              cl_float(67337235.0,precision)*Power(rxj,FOUR) - cl_float(1659700.0,precision)*Power(rxj,SIX) + cl_float(4089.0,precision)*Power(rxj,EIGHT)) 

    + NINE*Power(rxi,cl_float(18.0,precision))*Power(rxj,EIGHT)*

            (-cl_float(357591274425.0,precision) + cl_float(8328390840.0,precision)*Power(rxj,TWO) + 

              cl_float(912042180.0,precision)*Power(rxj,FOUR) - cl_float(12842480.0,precision)*Power(rxj,SIX) + 

              cl_float(10678.0,precision)*Power(rxj,EIGHT)) - 

           NINE*Power(rxi,cl_float(16.0,precision))*Power(rxj,cl_float(10.0,precision))*

            (cl_float(128599724925.0,precision) + cl_float(21298077360.0,precision)*Power(rxj,TWO) - 

              cl_float(267928500.0,precision)*Power(rxj,FOUR) - cl_float(5458320.0,precision)*Power(rxj,SIX) + 

              cl_float(14722.0,precision)*Power(rxj,EIGHT)) + 

           cl_float(18.0,precision)*Power(rxi,cl_float(24.0,precision))*Power(rxj,TWO)*

            (-cl_float(7604930025.0,precision) - cl_float(8866107180.0,precision)*Power(rxj,TWO) + cl_float(399272265.0,precision)*Power(rxj,FOUR) - 

              cl_float(5925780.0,precision)*Power(rxj,SIX) + cl_float(17651.0,precision)*Power(rxj,EIGHT)) - 

           NINE*Power(rxi,cl_float(22.0,precision))*Power(rxj,FOUR)*

            (cl_float(129194933175.0,precision) + cl_float(3909863160.0,precision)*Power(rxj,TWO) + cl_float(91420770.0,precision)*Power(rxj,FOUR) - 

              cl_float(8762040.0,precision)*Power(rxj,SIX) + cl_float(43928.0,precision)*Power(rxj,EIGHT))) + 

        exp(TWO*rxi)*Power(rxi,cl_float(12.0,precision))*

         (Power(rxi,EIGHT)*Power(rxj,cl_float(18.0,precision))*

            (cl_float(3218321469825.0,precision) - cl_float(341234165475.0,precision)*rxj - cl_float(393132783960.0,precision)*Power(rxj,TWO) - 

              cl_float(57092294070.0,precision)*Power(rxj,THREE) + cl_float(822786930.0,precision)*Power(rxj,FOUR) + 

              cl_float(982835910.0,precision)*Power(rxj,FIVE) + cl_float(106664040.0,precision)*Power(rxj,SIX) + 

              cl_float(4915116.0,precision)*Power(rxj,SEVEN) + cl_float(73602.0,precision)*Power(rxj,EIGHT) - cl_float(818.0,precision)*Power(rxj,NINE)) + 

           cl_float(10.0,precision)*Power(rxj,cl_float(26.0,precision))*(cl_float(352546425.0,precision) + cl_float(288447075.0,precision)*rxj + 

              cl_float(109884600.0,precision)*Power(rxj,TWO) + cl_float(25639740.0,precision)*Power(rxj,THREE) + 

              cl_float(4048380.0,precision)*Power(rxj,FOUR) + cl_float(449820.0,precision)*Power(rxj,FIVE) + cl_float(35280.0,precision)*Power(rxj,SIX) + 

              cl_float(1890.0,precision)*Power(rxj,SEVEN) + cl_float(63.0,precision)*Power(rxj,EIGHT) + Power(rxj,NINE)) + 

           cl_float(30.0,precision)*Power(rxi,TWO)*Power(rxj,cl_float(24.0,precision))*

            (cl_float(4562958015.0,precision) + cl_float(3269982555.0,precision)*rxj + cl_float(1076869080.0,precision)*Power(rxj,TWO) + 

              cl_float(213664500.0,precision)*Power(rxj,THREE) + cl_float(28081620.0,precision)*Power(rxj,FOUR) + 

              cl_float(2523276.0,precision)*Power(rxj,FIVE) + cl_float(153552.0,precision)*Power(rxj,SIX) + cl_float(5982.0,precision)*Power(rxj,SEVEN) + 

              cl_float(129.0,precision)*Power(rxj,EIGHT) + Power(rxj,NINE)) - 

           cl_float(15.0,precision)*Power(rxi,cl_float(24.0,precision))*Power(rxj,TWO)*

            (-cl_float(89775.0,precision) - cl_float(161595.0,precision)*rxj - cl_float(143640.0,precision)*Power(rxj,TWO) - cl_float(83790.0,precision)*Power(rxj,THREE) - 

              cl_float(35910.0,precision)*Power(rxj,FOUR) - cl_float(11970.0,precision)*Power(rxj,FIVE) - cl_float(3192.0,precision)*Power(rxj,SIX) - 

              cl_float(684.0,precision)*Power(rxj,SEVEN) - cl_float(114.0,precision)*Power(rxj,EIGHT) + TWO*Power(rxj,NINE)) - 

           FIVE*Power(rxi,cl_float(26.0,precision))*(cl_float(14175.0,precision) + cl_float(25515.0,precision)*rxj + cl_float(22680.0,precision)*Power(rxj,TWO) + 

              cl_float(13230.0,precision)*Power(rxj,THREE) + cl_float(5670.0,precision)*Power(rxj,FOUR) + cl_float(1890.0,precision)*Power(rxj,FIVE) + 

              cl_float(504.0,precision)*Power(rxj,SIX) + cl_float(108.0,precision)*Power(rxj,SEVEN) + cl_float(18.0,precision)*Power(rxj,EIGHT) + 

              TWO*Power(rxj,NINE)) - cl_float(1938.0,precision)*Power(rxi,cl_float(14.0,precision))*Power(rxj,cl_float(12.0,precision))*

            (-cl_float(826875.0,precision) + cl_float(15824025.0,precision)*rxj - cl_float(23398200.0,precision)*Power(rxj,TWO) + 

              cl_float(12344850.0,precision)*Power(rxj,THREE) + cl_float(1244250.0,precision)*Power(rxj,FOUR) - 

              cl_float(384930.0,precision)*Power(rxj,FIVE) - cl_float(59640.0,precision)*Power(rxj,SIX) - cl_float(1848.0,precision)*Power(rxj,SEVEN) + 

              cl_float(84.0,precision)*Power(rxj,EIGHT) + FOUR*Power(rxj,NINE)) + 

           cl_float(1938.0,precision)*Power(rxi,cl_float(12.0,precision))*Power(rxj,cl_float(14.0,precision))*

            (cl_float(72476775.0,precision) - cl_float(180008325.0,precision)*rxj + cl_float(98907480.0,precision)*Power(rxj,TWO) + 

              cl_float(11224710.0,precision)*Power(rxj,THREE) - cl_float(4235490.0,precision)*Power(rxj,FOUR) - 

              cl_float(791910.0,precision)*Power(rxj,FIVE) - cl_float(31080.0,precision)*Power(rxj,SIX) + cl_float(2232.0,precision)*Power(rxj,SEVEN) + 

              cl_float(204.0,precision)*Power(rxj,EIGHT) + FOUR*Power(rxj,NINE)) + 

           cl_float(342.0,precision)*Power(rxi,cl_float(16.0,precision))*Power(rxj,cl_float(10.0,precision))*

            (cl_float(2409750.0,precision) + cl_float(3641400.0,precision)*rxj + cl_float(9424800.0,precision)*Power(rxj,TWO) - 

              cl_float(8193150.0,precision)*Power(rxj,THREE) + cl_float(6301050.0,precision)*Power(rxj,FOUR) + 

              cl_float(400470.0,precision)*Power(rxj,FIVE) - cl_float(143640.0,precision)*Power(rxj,SIX) - cl_float(15518.0,precision)*Power(rxj,SEVEN) - 

              cl_float(281.0,precision)*Power(rxj,EIGHT) + NINE*Power(rxj,NINE)) - 

           cl_float(171.0,precision)*Power(rxi,cl_float(10.0,precision))*Power(rxj,cl_float(16.0,precision))*

            (-cl_float(6768406575.0,precision) + cl_float(6280474725.0,precision)*rxj + cl_float(438336360.0,precision)*Power(rxj,TWO) - 

              cl_float(400731030.0,precision)*Power(rxj,THREE) - cl_float(74168430.0,precision)*Power(rxj,FOUR) - 

              cl_float(2490810.0,precision)*Power(rxj,FIVE) + cl_float(461160.0,precision)*Power(rxj,SIX) + cl_float(51244.0,precision)*Power(rxj,SEVEN) + 

              cl_float(1858.0,precision)*Power(rxj,EIGHT) + cl_float(18.0,precision)*Power(rxj,NINE)) + 

           NINE*Power(rxi,cl_float(22.0,precision))*Power(rxj,FOUR)*

            (-cl_float(1346625.0,precision) - cl_float(2423925.0,precision)*rxj - cl_float(2154600.0,precision)*Power(rxj,TWO) - 

              cl_float(1256850.0,precision)*Power(rxj,THREE) - cl_float(538650.0,precision)*Power(rxj,FOUR) - 

              cl_float(179550.0,precision)*Power(rxj,FIVE) - cl_float(47880.0,precision)*Power(rxj,SIX) - cl_float(14264.0,precision)*Power(rxj,SEVEN) + 

              cl_float(292.0,precision)*Power(rxj,EIGHT) + cl_float(52.0,precision)*Power(rxj,NINE)) - 

           NINE*Power(rxi,FOUR)*Power(rxj,cl_float(22.0,precision))*

            (-cl_float(129194933175.0,precision) - cl_float(73043543475.0,precision)*rxj - cl_float(17732214360.0,precision)*Power(rxj,TWO) - 

              cl_float(2275149870.0,precision)*Power(rxj,THREE) - cl_float(134674470.0,precision)*Power(rxj,FOUR) + 

              cl_float(3148110.0,precision)*Power(rxj,FIVE) + cl_float(1197000.0,precision)*Power(rxj,SIX) + 

              cl_float(93176.0,precision)*Power(rxj,SEVEN) + cl_float(3452.0,precision)*Power(rxj,EIGHT) + cl_float(52.0,precision)*Power(rxj,NINE)) + 

           NINE*Power(rxi,SIX)*Power(rxj,cl_float(20.0,precision))*

            (cl_float(356863797675.0,precision) + cl_float(115054179975.0,precision)*rxj + cl_float(3909863160.0,precision)*Power(rxj,TWO) - 

              cl_float(3706015530.0,precision)*Power(rxj,THREE) - cl_float(798544530.0,precision)*Power(rxj,FOUR) - 

              cl_float(75669510.0,precision)*Power(rxj,FIVE) - cl_float(3319400.0,precision)*Power(rxj,SIX) - 

              cl_float(6456.0,precision)*Power(rxj,SEVEN) + cl_float(5188.0,precision)*Power(rxj,EIGHT) + cl_float(148.0,precision)*Power(rxj,NINE)) - 

           NINE*Power(rxi,cl_float(20.0,precision))*Power(rxj,SIX)*

            (-cl_float(7630875.0,precision) - cl_float(13735575.0,precision)*rxj - cl_float(12209400.0,precision)*Power(rxj,TWO) - 

              cl_float(7122150.0,precision)*Power(rxj,THREE) - cl_float(3052350.0,precision)*Power(rxj,FOUR) - 

              cl_float(777210.0,precision)*Power(rxj,FIVE) - cl_float(591640.0,precision)*Power(rxj,SIX) + cl_float(3064.0,precision)*Power(rxj,SEVEN) + 

              cl_float(5468.0,precision)*Power(rxj,EIGHT) + cl_float(148.0,precision)*Power(rxj,NINE)) + 

           TWO*Power(rxi,cl_float(18.0,precision))*Power(rxj,EIGHT)*

            (-cl_float(137355750.0,precision) - cl_float(247240350.0,precision)*rxj - cl_float(219769200.0,precision)*Power(rxj,TWO) - 

              cl_float(151171650.0,precision)*Power(rxj,THREE) + cl_float(13976550.0,precision)*Power(rxj,FOUR) - 

              cl_float(66692430.0,precision)*Power(rxj,FIVE) - cl_float(1640520.0,precision)*Power(rxj,SIX) + 

              cl_float(1046142.0,precision)*Power(rxj,SEVEN) + cl_float(66249.0,precision)*Power(rxj,EIGHT) + cl_float(409.0,precision)*Power(rxj,NINE))))/

      (cl_float(70875.0,precision)*exp(TWO*(rxi + rxj))*Power(rxi - rxj,cl_float(19.0,precision))*Power(rxi + rxj,cl_float(19.0,precision)))

     ; }
   
  }
  return S;
}

cl_F Slater_5S_1S(cl_F r,cl_F xi,cl_F xj)
{
  return Slater_1S_5S(r,xj,xi);
}

cl_F Slater_5S_2S(cl_F r,cl_F xi,cl_F xj)
{
  return Slater_2S_5S(r,xj,xi);
}

cl_F Slater_5S_3S(cl_F r,cl_F xi,cl_F xj)
{
  return Slater_3S_5S(r,xj,xi);
}

cl_F Slater_5S_4S(cl_F r,cl_F xi,cl_F xj)
{
  return Slater_4S_5S(r,xj,xi);
}

cl_F DSlater_1S_1S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = ZERO

    ; } else {  S = -(-cl_float(33.0,precision)*xi + cl_float(48.0,precision)*exp(TWO*r*xi)*xi - cl_float(36.0,precision)*r*Power(xi,TWO) - 

          cl_float(12.0,precision)*Power(r,TWO)*Power(xi,THREE))/(cl_float(24.0,precision)*exp(TWO*r*xi)*r) + 

      (-cl_float(24.0,precision) + cl_float(24.0,precision)*exp(TWO*r*xi) - cl_float(33.0,precision)*r*xi - cl_float(18.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

         FOUR*Power(r,THREE)*Power(xi,THREE))/(cl_float(24.0,precision)*exp(TWO*r*xi)*Power(r,TWO)) + 

      (xi*(-cl_float(24.0,precision) + cl_float(24.0,precision)*exp(TWO*r*xi) - cl_float(33.0,precision)*r*xi - cl_float(18.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

           FOUR*Power(r,THREE)*Power(xi,THREE)))/(cl_float(12.0,precision)*exp(TWO*r*xi)*r)

    ; }
 
  }
  else {
      if (r == ZERO) {  S = ZERO

    ; } else {  S = (exp(TWO*r*(xi + xj))*Power(Power(xi,TWO) - Power(xj,TWO),THREE) + 

         exp(TWO*r*xj)*Power(xj,FOUR)*

          (-THREE*Power(xi,TWO) - r*Power(xi,THREE) + Power(xj,TWO) + r*xi*Power(xj,TWO)) - 

         exp(TWO*r*xi)*Power(xi,FOUR)*

          (Power(xi,TWO)*(ONE + r*xj) - Power(xj,TWO)*(THREE + r*xj)))/

       (exp(TWO*r*(xi + xj))*Power(r,TWO)*Power(xi - xj,THREE)*Power(xi + xj,THREE)) + 

      (TWO*(exp(TWO*r*(xi + xj))*Power(Power(xi,TWO) - Power(xj,TWO),THREE) + 

           exp(TWO*r*xj)*Power(xj,FOUR)*

            (-THREE*Power(xi,TWO) - r*Power(xi,THREE) + Power(xj,TWO) + r*xi*Power(xj,TWO)) - 

           exp(TWO*r*xi)*Power(xi,FOUR)*

            (Power(xi,TWO)*(ONE + r*xj) - Power(xj,TWO)*(THREE + r*xj))))/

       (exp(TWO*r*(xi + xj))*r*Power(xi - xj,THREE)*Power(xi + xj,TWO)) - 

      (TWO*exp(TWO*r*(xi + xj))*(xi + xj)*Power(Power(xi,TWO) - Power(xj,TWO),THREE) + 

         exp(TWO*r*xj)*Power(xj,FOUR)*(-Power(xi,THREE) + xi*Power(xj,TWO)) + 

         TWO*exp(TWO*r*xj)*Power(xj,FIVE)*

          (-THREE*Power(xi,TWO) - r*Power(xi,THREE) + Power(xj,TWO) + r*xi*Power(xj,TWO)) - 

         exp(TWO*r*xi)*Power(xi,FOUR)*(Power(xi,TWO)*xj - Power(xj,THREE)) - 

         TWO*exp(TWO*r*xi)*Power(xi,FIVE)*

          (Power(xi,TWO)*(ONE + r*xj) - Power(xj,TWO)*(THREE + r*xj)))/

       (exp(TWO*r*(xi + xj))*r*Power(xi - xj,THREE)*Power(xi + xj,THREE))

    ; }
   
  }
  return S;
}

cl_F DSlater_1S_2S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = ZERO

    ; } else {  S = -(-cl_float(375.0,precision)*xi + cl_float(480.0,precision)*exp(TWO*r*xi)*xi - cl_float(540.0,precision)*r*Power(xi,TWO) - 

          cl_float(345.0,precision)*Power(r,TWO)*Power(xi,THREE) - cl_float(120.0,precision)*Power(r,THREE)*Power(xi,FOUR) - 

          cl_float(20.0,precision)*Power(r,FOUR)*Power(xi,FIVE))/(cl_float(240.0,precision)*exp(TWO*r*xi)*r) + 

      (-cl_float(240.0,precision) + cl_float(240.0,precision)*exp(TWO*r*xi) - cl_float(375.0,precision)*r*xi - cl_float(270.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

         cl_float(115.0,precision)*Power(r,THREE)*Power(xi,THREE) - cl_float(30.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - 

         FOUR*Power(r,FIVE)*Power(xi,FIVE))/(cl_float(240.0,precision)*exp(TWO*r*xi)*Power(r,TWO)) + 

      (xi*(-cl_float(240.0,precision) + cl_float(240.0,precision)*exp(TWO*r*xi) - cl_float(375.0,precision)*r*xi - cl_float(270.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

           cl_float(115.0,precision)*Power(r,THREE)*Power(xi,THREE) - cl_float(30.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - 

           FOUR*Power(r,FIVE)*Power(xi,FIVE)))/(cl_float(120.0,precision)*exp(TWO*r*xi)*r)

    ; }
 
  }
  else {
      if (r == ZERO) {  S = ZERO

    ; } else {  S = (SIX*exp(TWO*r*(xi + xj))*Power(Power(xi,TWO) - Power(xj,TWO),FIVE) + 

         SIX*exp(TWO*r*xj)*Power(xj,SIX)*

          (-FOUR*Power(xi,FOUR) - r*Power(xi,FIVE) - FIVE*Power(xi,TWO)*Power(xj,TWO) + 

            Power(xj,FOUR) + r*xi*Power(xj,FOUR)) - 

         exp(TWO*r*xi)*Power(xi,FOUR)*

          (Power(xi,SIX)*(SIX + NINE*r*xj + SIX*Power(r,TWO)*Power(xj,TWO) + 

               TWO*Power(r,THREE)*Power(xj,THREE)) - 

            THREE*Power(xi,FOUR)*Power(xj,TWO)*

             (cl_float(10.0,precision) + cl_float(15.0,precision)*r*xj + cl_float(10.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               TWO*Power(r,THREE)*Power(xj,THREE)) + 

            THREE*Power(xi,TWO)*Power(xj,FOUR)*

             (cl_float(20.0,precision) + cl_float(33.0,precision)*r*xj + cl_float(14.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               TWO*Power(r,THREE)*Power(xj,THREE)) - 

            Power(xj,SIX)*(cl_float(84.0,precision) + cl_float(63.0,precision)*r*xj + cl_float(18.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               TWO*Power(r,THREE)*Power(xj,THREE))))/

       (cl_float(6.0,precision)*exp(TWO*r*(xi + xj))*Power(r,TWO)*Power(xi - xj,FIVE)*Power(xi + xj,FIVE)) + 

      (SIX*exp(TWO*r*(xi + xj))*Power(Power(xi,TWO) - Power(xj,TWO),FIVE) + 

         SIX*exp(TWO*r*xj)*Power(xj,SIX)*

          (-FOUR*Power(xi,FOUR) - r*Power(xi,FIVE) - FIVE*Power(xi,TWO)*Power(xj,TWO) + 

            Power(xj,FOUR) + r*xi*Power(xj,FOUR)) - 

         exp(TWO*r*xi)*Power(xi,FOUR)*

          (Power(xi,SIX)*(SIX + NINE*r*xj + SIX*Power(r,TWO)*Power(xj,TWO) + 

               TWO*Power(r,THREE)*Power(xj,THREE)) - 

            THREE*Power(xi,FOUR)*Power(xj,TWO)*

             (cl_float(10.0,precision) + cl_float(15.0,precision)*r*xj + cl_float(10.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               TWO*Power(r,THREE)*Power(xj,THREE)) + 

            THREE*Power(xi,TWO)*Power(xj,FOUR)*

             (cl_float(20.0,precision) + cl_float(33.0,precision)*r*xj + cl_float(14.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               TWO*Power(r,THREE)*Power(xj,THREE)) - 

            Power(xj,SIX)*(cl_float(84.0,precision) + cl_float(63.0,precision)*r*xj + cl_float(18.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               TWO*Power(r,THREE)*Power(xj,THREE))))/

       (cl_float(3.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,FIVE)*Power(xi + xj,FOUR)) - 

      (cl_float(12.0,precision)*exp(TWO*r*(xi + xj))*(xi + xj)*Power(Power(xi,TWO) - Power(xj,TWO),FIVE) + 

         SIX*exp(TWO*r*xj)*Power(xj,SIX)*(-Power(xi,FIVE) + xi*Power(xj,FOUR)) + 

         cl_float(12.0,precision)*exp(TWO*r*xj)*Power(xj,SEVEN)*

          (-FOUR*Power(xi,FOUR) - r*Power(xi,FIVE) - FIVE*Power(xi,TWO)*Power(xj,TWO) + 

            Power(xj,FOUR) + r*xi*Power(xj,FOUR)) - 

         exp(TWO*r*xi)*Power(xi,FOUR)*

          (Power(xi,SIX)*(NINE*xj + cl_float(12.0,precision)*r*Power(xj,TWO) + SIX*Power(r,TWO)*Power(xj,THREE)) - 

            THREE*Power(xi,FOUR)*Power(xj,TWO)*

             (cl_float(15.0,precision)*xj + cl_float(20.0,precision)*r*Power(xj,TWO) + SIX*Power(r,TWO)*Power(xj,THREE)) + 

            THREE*Power(xi,TWO)*Power(xj,FOUR)*

             (cl_float(33.0,precision)*xj + cl_float(28.0,precision)*r*Power(xj,TWO) + SIX*Power(r,TWO)*Power(xj,THREE)) - 

            Power(xj,SIX)*(cl_float(63.0,precision)*xj + cl_float(36.0,precision)*r*Power(xj,TWO) + SIX*Power(r,TWO)*Power(xj,THREE))) - 

         TWO*exp(TWO*r*xi)*Power(xi,FIVE)*

          (Power(xi,SIX)*(SIX + NINE*r*xj + SIX*Power(r,TWO)*Power(xj,TWO) + 

               TWO*Power(r,THREE)*Power(xj,THREE)) - 

            THREE*Power(xi,FOUR)*Power(xj,TWO)*

             (cl_float(10.0,precision) + cl_float(15.0,precision)*r*xj + cl_float(10.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               TWO*Power(r,THREE)*Power(xj,THREE)) + 

            THREE*Power(xi,TWO)*Power(xj,FOUR)*

             (cl_float(20.0,precision) + cl_float(33.0,precision)*r*xj + cl_float(14.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               TWO*Power(r,THREE)*Power(xj,THREE)) - 

            Power(xj,SIX)*(cl_float(84.0,precision) + cl_float(63.0,precision)*r*xj + cl_float(18.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               TWO*Power(r,THREE)*Power(xj,THREE))))/

       (cl_float(6.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,FIVE)*Power(xi + xj,FIVE))

    ; }
   
  }
  return S;
}

cl_F DSlater_1S_3S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = ZERO

    ; } else {  S = -(-cl_float(203175.0,precision)*xi + cl_float(241920.0,precision)*exp(TWO*r*xi)*xi - cl_float(328860.0,precision)*r*Power(xi,TWO) - 

          cl_float(253260.0,precision)*Power(r,TWO)*Power(xi,THREE) - cl_float(120960.0,precision)*Power(r,THREE)*Power(xi,FOUR) - 

          cl_float(38640.0,precision)*Power(r,FOUR)*Power(xi,FIVE) - cl_float(8064.0,precision)*Power(r,FIVE)*Power(xi,SIX) - 

          cl_float(896.0,precision)*Power(r,SIX)*Power(xi,SEVEN))/(cl_float(120960.0,precision)*exp(TWO*r*xi)*r) + 

      (-cl_float(120960.0,precision) + cl_float(120960.0,precision)*exp(TWO*r*xi) - cl_float(203175.0,precision)*r*xi - 

         cl_float(164430.0,precision)*Power(r,TWO)*Power(xi,TWO) - cl_float(84420.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

         cl_float(30240.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - cl_float(7728.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

         cl_float(1344.0,precision)*Power(r,SIX)*Power(xi,SIX) - cl_float(128.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN))/

       (cl_float(120960.0,precision)*exp(TWO*r*xi)*Power(r,TWO)) + 

      (xi*(-cl_float(120960.0,precision) + cl_float(120960.0,precision)*exp(TWO*r*xi) - cl_float(203175.0,precision)*r*xi - 

           cl_float(164430.0,precision)*Power(r,TWO)*Power(xi,TWO) - cl_float(84420.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

           cl_float(30240.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - cl_float(7728.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

           cl_float(1344.0,precision)*Power(r,SIX)*Power(xi,SIX) - cl_float(128.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN)))/

       (cl_float(60480.0,precision)*exp(TWO*r*xi)*r)

    ; }
 
  }
  else {
      if (r == ZERO) {  S = ZERO

    ; } else {  S = (cl_float(45.0,precision)*exp(TWO*r*(xi + xj))*Power(Power(xi,TWO) - Power(xj,TWO),SEVEN) + 

         cl_float(15.0,precision)*exp(TWO*r*xj)*Power(xj,EIGHT)*

          (-cl_float(15.0,precision)*Power(xi,SIX) - THREE*r*Power(xi,SEVEN) - cl_float(63.0,precision)*Power(xi,FOUR)*Power(xj,TWO) - 

            SEVEN*r*Power(xi,FIVE)*Power(xj,TWO) - cl_float(21.0,precision)*Power(xi,TWO)*Power(xj,FOUR) + 

            SEVEN*r*Power(xi,THREE)*Power(xj,FOUR) + THREE*Power(xj,SIX) + THREE*r*xi*Power(xj,SIX)) + 

         exp(TWO*r*xi)*Power(xi,FOUR)*

          (-cl_float(10.0,precision)*Power(xi,TWO)*Power(xj,EIGHT)*

             (cl_float(135.0,precision) + cl_float(333.0,precision)*r*xj + cl_float(228.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(75.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(13.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               Power(r,FIVE)*Power(xj,FIVE)) + 

            TWO*Power(xj,cl_float(10.0,precision))*(cl_float(945.0,precision) + cl_float(945.0,precision)*r*xj + cl_float(420.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(105.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(15.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               Power(r,FIVE)*Power(xj,FIVE)) - 

            Power(xi,cl_float(10.0,precision))*(cl_float(45.0,precision) + cl_float(75.0,precision)*r*xj + cl_float(60.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(30.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(10.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               TWO*Power(r,FIVE)*Power(xj,FIVE)) + 

            FIVE*Power(xi,EIGHT)*Power(xj,TWO)*

             (cl_float(63.0,precision) + cl_float(105.0,precision)*r*xj + cl_float(84.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(42.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(14.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               TWO*Power(r,FIVE)*Power(xj,FIVE)) - 

            FIVE*Power(xi,SIX)*Power(xj,FOUR)*

             (cl_float(189.0,precision) + cl_float(315.0,precision)*r*xj + cl_float(252.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(132.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(36.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               FOUR*Power(r,FIVE)*Power(xj,FIVE)) + 

            FIVE*Power(xi,FOUR)*Power(xj,SIX)*

             (cl_float(315.0,precision) + cl_float(513.0,precision)*r*xj + cl_float(468.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(204.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(44.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               FOUR*Power(r,FIVE)*Power(xj,FIVE))))/

       (cl_float(45.0,precision)*exp(TWO*r*(xi + xj))*Power(r,TWO)*Power(xi - xj,SEVEN)*Power(xi + xj,SEVEN)) 

    + (TWO*(cl_float(45.0,precision)*exp(TWO*r*(xi + xj))*Power(Power(xi,TWO) - Power(xj,TWO),SEVEN) + 

           cl_float(15.0,precision)*exp(TWO*r*xj)*Power(xj,EIGHT)*

            (-cl_float(15.0,precision)*Power(xi,SIX) - THREE*r*Power(xi,SEVEN) - cl_float(63.0,precision)*Power(xi,FOUR)*Power(xj,TWO) - 

              SEVEN*r*Power(xi,FIVE)*Power(xj,TWO) - cl_float(21.0,precision)*Power(xi,TWO)*Power(xj,FOUR) + 

              SEVEN*r*Power(xi,THREE)*Power(xj,FOUR) + THREE*Power(xj,SIX) + THREE*r*xi*Power(xj,SIX)) 

    + exp(TWO*r*xi)*Power(xi,FOUR)*(-cl_float(10.0,precision)*Power(xi,TWO)*Power(xj,EIGHT)*

               (cl_float(135.0,precision) + cl_float(333.0,precision)*r*xj + cl_float(228.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(75.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(13.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 Power(r,FIVE)*Power(xj,FIVE)) + 

              TWO*Power(xj,cl_float(10.0,precision))*(cl_float(945.0,precision) + cl_float(945.0,precision)*r*xj + cl_float(420.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(105.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(15.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 Power(r,FIVE)*Power(xj,FIVE)) - 

              Power(xi,cl_float(10.0,precision))*(cl_float(45.0,precision) + cl_float(75.0,precision)*r*xj + cl_float(60.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(30.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(10.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 TWO*Power(r,FIVE)*Power(xj,FIVE)) + 

              FIVE*Power(xi,EIGHT)*Power(xj,TWO)*

               (cl_float(63.0,precision) + cl_float(105.0,precision)*r*xj + cl_float(84.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(42.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(14.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 TWO*Power(r,FIVE)*Power(xj,FIVE)) - 

              FIVE*Power(xi,SIX)*Power(xj,FOUR)*

               (cl_float(189.0,precision) + cl_float(315.0,precision)*r*xj + cl_float(252.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(132.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(36.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 FOUR*Power(r,FIVE)*Power(xj,FIVE)) + 

              FIVE*Power(xi,FOUR)*Power(xj,SIX)*

               (cl_float(315.0,precision) + cl_float(513.0,precision)*r*xj + cl_float(468.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(204.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(44.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 FOUR*Power(r,FIVE)*Power(xj,FIVE)))))/

       (cl_float(45.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,SEVEN)*Power(xi + xj,SIX)) - 

      (cl_float(90.0,precision)*exp(TWO*r*(xi + xj))*(xi + xj)*Power(Power(xi,TWO) - Power(xj,TWO),SEVEN) + 

         cl_float(15.0,precision)*exp(TWO*r*xj)*Power(xj,EIGHT)*

          (-THREE*Power(xi,SEVEN) - SEVEN*Power(xi,FIVE)*Power(xj,TWO) + 

            SEVEN*Power(xi,THREE)*Power(xj,FOUR) + THREE*xi*Power(xj,SIX)) + 

         cl_float(30.0,precision)*exp(TWO*r*xj)*Power(xj,NINE)*

          (-cl_float(15.0,precision)*Power(xi,SIX) - THREE*r*Power(xi,SEVEN) - cl_float(63.0,precision)*Power(xi,FOUR)*Power(xj,TWO) - 

            SEVEN*r*Power(xi,FIVE)*Power(xj,TWO) - cl_float(21.0,precision)*Power(xi,TWO)*Power(xj,FOUR) + 

            SEVEN*r*Power(xi,THREE)*Power(xj,FOUR) + THREE*Power(xj,SIX) + THREE*r*xi*Power(xj,SIX)) + 

         exp(TWO*r*xi)*Power(xi,FOUR)*

          (-cl_float(10.0,precision)*Power(xi,TWO)*Power(xj,EIGHT)*

             (cl_float(333.0,precision)*xj + cl_float(456.0,precision)*r*Power(xj,TWO) + cl_float(225.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(52.0,precision)*Power(r,THREE)*Power(xj,FOUR) + FIVE*Power(r,FOUR)*Power(xj,FIVE)) + 

            TWO*Power(xj,cl_float(10.0,precision))*(cl_float(945.0,precision)*xj + cl_float(840.0,precision)*r*Power(xj,TWO) + 

               cl_float(315.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(60.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               FIVE*Power(r,FOUR)*Power(xj,FIVE)) - 

            Power(xi,cl_float(10.0,precision))*(cl_float(75.0,precision)*xj + cl_float(120.0,precision)*r*Power(xj,TWO) + 

               cl_float(90.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(40.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(10.0,precision)*Power(r,FOUR)*Power(xj,FIVE)) + 

            FIVE*Power(xi,EIGHT)*Power(xj,TWO)*

             (cl_float(105.0,precision)*xj + cl_float(168.0,precision)*r*Power(xj,TWO) + cl_float(126.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(56.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(10.0,precision)*Power(r,FOUR)*Power(xj,FIVE)) - 

            FIVE*Power(xi,SIX)*Power(xj,FOUR)*

             (cl_float(315.0,precision)*xj + cl_float(504.0,precision)*r*Power(xj,TWO) + cl_float(396.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(144.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(20.0,precision)*Power(r,FOUR)*Power(xj,FIVE)) + 

            FIVE*Power(xi,FOUR)*Power(xj,SIX)*

             (cl_float(513.0,precision)*xj + cl_float(936.0,precision)*r*Power(xj,TWO) + cl_float(612.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(176.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(20.0,precision)*Power(r,FOUR)*Power(xj,FIVE))) + 

         TWO*exp(TWO*r*xi)*Power(xi,FIVE)*

          (-cl_float(10.0,precision)*Power(xi,TWO)*Power(xj,EIGHT)*

             (cl_float(135.0,precision) + cl_float(333.0,precision)*r*xj + cl_float(228.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(75.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(13.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               Power(r,FIVE)*Power(xj,FIVE)) + 

            TWO*Power(xj,cl_float(10.0,precision))*(cl_float(945.0,precision) + cl_float(945.0,precision)*r*xj + cl_float(420.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(105.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(15.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               Power(r,FIVE)*Power(xj,FIVE)) - 

            Power(xi,cl_float(10.0,precision))*(cl_float(45.0,precision) + cl_float(75.0,precision)*r*xj + cl_float(60.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(30.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(10.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               TWO*Power(r,FIVE)*Power(xj,FIVE)) + 

            FIVE*Power(xi,EIGHT)*Power(xj,TWO)*

             (cl_float(63.0,precision) + cl_float(105.0,precision)*r*xj + cl_float(84.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(42.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(14.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               TWO*Power(r,FIVE)*Power(xj,FIVE)) - 

            FIVE*Power(xi,SIX)*Power(xj,FOUR)*

             (cl_float(189.0,precision) + cl_float(315.0,precision)*r*xj + cl_float(252.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(132.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(36.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               FOUR*Power(r,FIVE)*Power(xj,FIVE)) + 

            FIVE*Power(xi,FOUR)*Power(xj,SIX)*

             (cl_float(315.0,precision) + cl_float(513.0,precision)*r*xj + cl_float(468.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(204.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(44.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               FOUR*Power(r,FIVE)*Power(xj,FIVE))))/

       (cl_float(45.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,SEVEN)*Power(xi + xj,SEVEN))

    ; }
   
  }
  return S;
}

cl_F DSlater_1S_4S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = ZERO

    ; } else {  S = -(-cl_float(5088825.0,precision)*xi + cl_float(5806080.0,precision)*exp(TWO*r*xi)*xi - cl_float(8743140.0,precision)*r*Power(xi,TWO) - 

          cl_float(7319970.0,precision)*Power(r,TWO)*Power(xi,THREE) - cl_float(3946320.0,precision)*Power(r,THREE)*Power(xi,FOUR) - 

          cl_float(1519560.0,precision)*Power(r,FOUR)*Power(xi,FIVE) - cl_float(435456.0,precision)*Power(r,FIVE)*Power(xi,SIX) - 

          cl_float(92736.0,precision)*Power(r,SIX)*Power(xi,SEVEN) - cl_float(13824.0,precision)*Power(r,SEVEN)*Power(xi,EIGHT) - 

          cl_float(1152.0,precision)*Power(r,EIGHT)*Power(xi,NINE))/(cl_float(2.90304e6,precision)*exp(TWO*r*xi)*r) + 

      (-cl_float(2903040.0,precision) + cl_float(2903040.0,precision)*exp(TWO*r*xi) - cl_float(5088825.0,precision)*r*xi - 

         cl_float(4371570.0,precision)*Power(r,TWO)*Power(xi,TWO) - cl_float(2439990.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

         cl_float(986580.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - cl_float(303912.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

         cl_float(72576.0,precision)*Power(r,SIX)*Power(xi,SIX) - cl_float(13248.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

         cl_float(1728.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - cl_float(128.0,precision)*Power(r,NINE)*Power(xi,NINE))/

       (cl_float(2.90304e6,precision)*exp(TWO*r*xi)*Power(r,TWO)) + 

      (xi*(-cl_float(2903040.0,precision) + cl_float(2903040.0,precision)*exp(TWO*r*xi) - cl_float(5088825.0,precision)*r*xi - 

           cl_float(4371570.0,precision)*Power(r,TWO)*Power(xi,TWO) - cl_float(2439990.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

           cl_float(986580.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - cl_float(303912.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

           cl_float(72576.0,precision)*Power(r,SIX)*Power(xi,SIX) - cl_float(13248.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

           cl_float(1728.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - cl_float(128.0,precision)*Power(r,NINE)*Power(xi,NINE)))/

       (cl_float(1.45152e6,precision)*exp(TWO*r*xi)*r)

    ; }
 
  }
  else {
      if (r == ZERO) {  S = ZERO

    ; } else {  S = (cl_float(1260.0,precision)*exp(TWO*r*(xi + xj))*Power(Power(xi,TWO) - Power(xj,TWO),NINE) + 

         cl_float(1260.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(10.0,precision))*

          (-SIX*Power(xi,EIGHT) - r*Power(xi,NINE) - cl_float(51.0,precision)*Power(xi,SIX)*Power(xj,TWO) - 

            SIX*r*Power(xi,SEVEN)*Power(xj,TWO) - cl_float(63.0,precision)*Power(xi,FOUR)*Power(xj,FOUR) - 

            NINE*Power(xi,TWO)*Power(xj,SIX) + SIX*r*Power(xi,THREE)*Power(xj,SIX) + 

            Power(xj,EIGHT) + r*xi*Power(xj,EIGHT)) + 

         exp(TWO*r*xi)*Power(xi,FOUR)*

          (-cl_float(42.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,FOUR)*

             (cl_float(1080.0,precision) + cl_float(1890.0,precision)*r*xj + cl_float(1620.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(900.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(360.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(111.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(22.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               TWO*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            cl_float(70.0,precision)*Power(xi,EIGHT)*Power(xj,SIX)*

             (cl_float(1512.0,precision) + cl_float(2646.0,precision)*r*xj + cl_float(2268.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1248.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(528.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(153.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(26.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               TWO*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            cl_float(14.0,precision)*Power(xi,TWO)*Power(xj,cl_float(12.0,precision))*

             (cl_float(2970.0,precision) + cl_float(16335.0,precision)*r*xj + cl_float(15390.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(7110.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(1980.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(351.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(38.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               TWO*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            TWO*Power(xj,cl_float(14.0,precision))*(cl_float(62370.0,precision) + cl_float(72765.0,precision)*r*xj + 

               cl_float(39690.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(13230.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(2940.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(441.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(42.0,precision)*Power(r,SIX)*Power(xj,SIX) + TWO*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            Power(xi,cl_float(14.0,precision))*(cl_float(1260.0,precision) + cl_float(2205.0,precision)*r*xj + cl_float(1890.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1050.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(420.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(126.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            SEVEN*Power(xi,cl_float(12.0,precision))*Power(xj,TWO)*

             (cl_float(1620.0,precision) + cl_float(2835.0,precision)*r*xj + cl_float(2430.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1350.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(540.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(162.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(36.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            cl_float(35.0,precision)*Power(xi,SIX)*Power(xj,EIGHT)*

             (cl_float(4536.0,precision) + cl_float(7983.0,precision)*r*xj + cl_float(6534.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(4014.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(1644.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(414.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(60.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            cl_float(21.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(10.0,precision))*

             (cl_float(7920.0,precision) + cl_float(11385.0,precision)*r*xj + cl_float(12330.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(7410.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(2580.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(546.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(68.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN))))/

       (cl_float(1260.0,precision)*exp(TWO*r*(xi + xj))*Power(r,TWO)*Power(xi - xj,NINE)*

         Power(xi + xj,NINE)) + (cl_float(1260.0,precision)*exp(TWO*r*(xi + xj))*

          Power(Power(xi,TWO) - Power(xj,TWO),NINE) + 

         cl_float(1260.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(10.0,precision))*

          (-SIX*Power(xi,EIGHT) - r*Power(xi,NINE) - cl_float(51.0,precision)*Power(xi,SIX)*Power(xj,TWO) - 

            SIX*r*Power(xi,SEVEN)*Power(xj,TWO) - cl_float(63.0,precision)*Power(xi,FOUR)*Power(xj,FOUR) - 

            NINE*Power(xi,TWO)*Power(xj,SIX) + SIX*r*Power(xi,THREE)*Power(xj,SIX) + 

            Power(xj,EIGHT) + r*xi*Power(xj,EIGHT)) + 

         exp(TWO*r*xi)*Power(xi,FOUR)*

          (-cl_float(42.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,FOUR)*

             (cl_float(1080.0,precision) + cl_float(1890.0,precision)*r*xj + cl_float(1620.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(900.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(360.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(111.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(22.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               TWO*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            cl_float(70.0,precision)*Power(xi,EIGHT)*Power(xj,SIX)*

             (cl_float(1512.0,precision) + cl_float(2646.0,precision)*r*xj + cl_float(2268.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1248.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(528.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(153.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(26.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               TWO*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            cl_float(14.0,precision)*Power(xi,TWO)*Power(xj,cl_float(12.0,precision))*

             (cl_float(2970.0,precision) + cl_float(16335.0,precision)*r*xj + cl_float(15390.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(7110.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(1980.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(351.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(38.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               TWO*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            TWO*Power(xj,cl_float(14.0,precision))*(cl_float(62370.0,precision) + cl_float(72765.0,precision)*r*xj + 

               cl_float(39690.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(13230.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(2940.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(441.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(42.0,precision)*Power(r,SIX)*Power(xj,SIX) + TWO*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            Power(xi,cl_float(14.0,precision))*(cl_float(1260.0,precision) + cl_float(2205.0,precision)*r*xj + cl_float(1890.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1050.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(420.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(126.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            SEVEN*Power(xi,cl_float(12.0,precision))*Power(xj,TWO)*

             (cl_float(1620.0,precision) + cl_float(2835.0,precision)*r*xj + cl_float(2430.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1350.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(540.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(162.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(36.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            cl_float(35.0,precision)*Power(xi,SIX)*Power(xj,EIGHT)*

             (cl_float(4536.0,precision) + cl_float(7983.0,precision)*r*xj + cl_float(6534.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(4014.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(1644.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(414.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(60.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            cl_float(21.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(10.0,precision))*

             (cl_float(7920.0,precision) + cl_float(11385.0,precision)*r*xj + cl_float(12330.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(7410.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(2580.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(546.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(68.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN))))/

       (cl_float(630.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,NINE)*Power(xi + xj,EIGHT)) - 

      (cl_float(2520.0,precision)*exp(TWO*r*(xi + xj))*(xi + xj)*

          Power(Power(xi,TWO) - Power(xj,TWO),NINE) + 

         cl_float(1260.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(10.0,precision))*

          (-Power(xi,NINE) - SIX*Power(xi,SEVEN)*Power(xj,TWO) + 

            SIX*Power(xi,THREE)*Power(xj,SIX) + xi*Power(xj,EIGHT)) + 

         cl_float(2520.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(11.0,precision))*

          (-SIX*Power(xi,EIGHT) - r*Power(xi,NINE) - cl_float(51.0,precision)*Power(xi,SIX)*Power(xj,TWO) - 

            SIX*r*Power(xi,SEVEN)*Power(xj,TWO) - cl_float(63.0,precision)*Power(xi,FOUR)*Power(xj,FOUR) - 

            NINE*Power(xi,TWO)*Power(xj,SIX) + SIX*r*Power(xi,THREE)*Power(xj,SIX) + 

            Power(xj,EIGHT) + r*xi*Power(xj,EIGHT)) + 

         exp(TWO*r*xi)*Power(xi,FOUR)*

          (-cl_float(42.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,FOUR)*

             (cl_float(1890.0,precision)*xj + cl_float(3240.0,precision)*r*Power(xj,TWO) + cl_float(2700.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(1440.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(555.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(132.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(14.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) + 

            cl_float(70.0,precision)*Power(xi,EIGHT)*Power(xj,SIX)*

             (cl_float(2646.0,precision)*xj + cl_float(4536.0,precision)*r*Power(xj,TWO) + cl_float(3744.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(2112.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(765.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(156.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(14.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) - 

            cl_float(14.0,precision)*Power(xi,TWO)*Power(xj,cl_float(12.0,precision))*

             (cl_float(16335.0,precision)*xj + cl_float(30780.0,precision)*r*Power(xj,TWO) + cl_float(21330.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(7920.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(1755.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(228.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(14.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) + 

            TWO*Power(xj,cl_float(14.0,precision))*(cl_float(72765.0,precision)*xj + cl_float(79380.0,precision)*r*Power(xj,TWO) + 

               cl_float(39690.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(11760.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(2205.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(252.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(14.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) - 

            Power(xi,cl_float(14.0,precision))*(cl_float(2205.0,precision)*xj + cl_float(3780.0,precision)*r*Power(xj,TWO) + 

               cl_float(3150.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(1680.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(630.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(168.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) + 

            SEVEN*Power(xi,cl_float(12.0,precision))*Power(xj,TWO)*

             (cl_float(2835.0,precision)*xj + cl_float(4860.0,precision)*r*Power(xj,TWO) + cl_float(4050.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(2160.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(810.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(216.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) - 

            cl_float(35.0,precision)*Power(xi,SIX)*Power(xj,EIGHT)*

             (cl_float(7983.0,precision)*xj + cl_float(13068.0,precision)*r*Power(xj,TWO) + cl_float(12042.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(6576.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(2070.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(360.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) + 

            cl_float(21.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(10.0,precision))*

             (cl_float(11385.0,precision)*xj + cl_float(24660.0,precision)*r*Power(xj,TWO) + cl_float(22230.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(10320.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(2730.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(408.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SEVEN))) + 

         TWO*exp(TWO*r*xi)*Power(xi,FIVE)*

          (-cl_float(42.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,FOUR)*

             (cl_float(1080.0,precision) + cl_float(1890.0,precision)*r*xj + cl_float(1620.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(900.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(360.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(111.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(22.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               TWO*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            cl_float(70.0,precision)*Power(xi,EIGHT)*Power(xj,SIX)*

             (cl_float(1512.0,precision) + cl_float(2646.0,precision)*r*xj + cl_float(2268.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1248.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(528.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(153.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(26.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               TWO*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            cl_float(14.0,precision)*Power(xi,TWO)*Power(xj,cl_float(12.0,precision))*

             (cl_float(2970.0,precision) + cl_float(16335.0,precision)*r*xj + cl_float(15390.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(7110.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(1980.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(351.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(38.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               TWO*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            TWO*Power(xj,cl_float(14.0,precision))*(cl_float(62370.0,precision) + cl_float(72765.0,precision)*r*xj + cl_float(39690.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(13230.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(2940.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(441.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(42.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               TWO*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            Power(xi,cl_float(14.0,precision))*(cl_float(1260.0,precision) + cl_float(2205.0,precision)*r*xj + cl_float(1890.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1050.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(420.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(126.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            SEVEN*Power(xi,cl_float(12.0,precision))*Power(xj,TWO)*

             (cl_float(1620.0,precision) + cl_float(2835.0,precision)*r*xj + cl_float(2430.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1350.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(540.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(162.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(36.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            cl_float(35.0,precision)*Power(xi,SIX)*Power(xj,EIGHT)*

             (cl_float(4536.0,precision) + cl_float(7983.0,precision)*r*xj + cl_float(6534.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(4014.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(1644.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(414.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(60.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            cl_float(21.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(10.0,precision))*

             (cl_float(7920.0,precision) + cl_float(11385.0,precision)*r*xj + cl_float(12330.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(7410.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(2580.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(546.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(68.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN))))/

       (cl_float(1260.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,NINE)*Power(xi + xj,NINE))

    ; }
   
  }
  return S;
}

cl_F DSlater_1S_5S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = ZERO

    ; } else {  S = -(-cl_float(2875101075.0,precision)*xi + cl_float(3193344000.0,precision)*exp(TWO*r*xi)*xi - 

          cl_float(5113716300.0,precision)*r*Power(xi,TWO) - cl_float(4478789700.0,precision)*Power(r,TWO)*Power(xi,THREE) - 

          cl_float(2564654400.0,precision)*Power(r,THREE)*Power(xi,FOUR) - 

          cl_float(1073595600.0,precision)*Power(r,FOUR)*Power(xi,FIVE) - 

          cl_float(347276160.0,precision)*Power(r,FIVE)*Power(xi,SIX) - cl_float(89147520.0,precision)*Power(r,SIX)*Power(xi,SEVEN) - 

          cl_float(18247680.0,precision)*Power(r,SEVEN)*Power(xi,EIGHT) - cl_float(2914560.0,precision)*Power(r,EIGHT)*Power(xi,NINE) - 

          cl_float(337920.0,precision)*Power(r,NINE)*Power(xi,cl_float(10.0,precision)) - cl_float(22528.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(11.0,precision)))/

       (cl_float(1.596672e9,precision)*exp(TWO*r*xi)*r) + 

      (-cl_float(1596672000.0,precision) + cl_float(1596672000.0,precision)*exp(TWO*r*xi) - cl_float(2875101075.0,precision)*r*xi - 

         cl_float(2556858150.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

         cl_float(1492929900.0,precision)*Power(r,THREE)*Power(xi,THREE) - cl_float(641163600.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - 

         cl_float(214719120.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - cl_float(57879360.0,precision)*Power(r,SIX)*Power(xi,SIX) - 

         cl_float(12735360.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - cl_float(2280960.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - 

         cl_float(323840.0,precision)*Power(r,NINE)*Power(xi,NINE) - cl_float(33792.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - 

         cl_float(2048.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)))/(cl_float(1.596672e9,precision)*exp(TWO*r*xi)*Power(r,TWO)) + 

      (xi*(-cl_float(1596672000.0,precision) + cl_float(1596672000.0,precision)*exp(TWO*r*xi) - cl_float(2875101075.0,precision)*r*xi - 

           cl_float(2556858150.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

           cl_float(1492929900.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

           cl_float(641163600.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - cl_float(214719120.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

           cl_float(57879360.0,precision)*Power(r,SIX)*Power(xi,SIX) - cl_float(12735360.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

           cl_float(2280960.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - cl_float(323840.0,precision)*Power(r,NINE)*Power(xi,NINE) - 

           cl_float(33792.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - cl_float(2048.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision))))/

       (cl_float(7.98336e8,precision)*exp(TWO*r*xi)*r)

    ; }
 
  }
  else {
      if (r == ZERO) {  S = ZERO

    ; } else {  S = (cl_float(14175.0,precision)*exp(TWO*r*(xi + xj))*Power(Power(xi,TWO) - Power(xj,TWO),cl_float(11.0,precision)) + 

         cl_float(2835.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(12.0,precision))*

          (-cl_float(35.0,precision)*Power(xi,cl_float(10.0,precision)) - FIVE*r*Power(xi,cl_float(11.0,precision)) - cl_float(495.0,precision)*Power(xi,EIGHT)*Power(xj,TWO) - 

            cl_float(55.0,precision)*r*Power(xi,NINE)*Power(xj,TWO) - cl_float(1254.0,precision)*Power(xi,SIX)*Power(xj,FOUR) - 

            cl_float(66.0,precision)*r*Power(xi,SEVEN)*Power(xj,FOUR) - cl_float(726.0,precision)*Power(xi,FOUR)*Power(xj,SIX) + 

            cl_float(66.0,precision)*r*Power(xi,FIVE)*Power(xj,SIX) - cl_float(55.0,precision)*Power(xi,TWO)*Power(xj,EIGHT) + 

            cl_float(55.0,precision)*r*Power(xi,THREE)*Power(xj,EIGHT) + FIVE*Power(xj,cl_float(10.0,precision)) + FIVE*r*xi*Power(xj,cl_float(10.0,precision))

    ) + exp(TWO*r*xi)*Power(xi,FOUR)*

          (-(Power(xi,cl_float(18.0,precision))*(cl_float(14175.0,precision) + cl_float(25515.0,precision)*r*xj + cl_float(22680.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(13230.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(5670.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(1890.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(504.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

                 cl_float(108.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

                 TWO*Power(r,NINE)*Power(xj,NINE))) + 

            NINE*Power(xi,cl_float(16.0,precision))*Power(xj,TWO)*

             (cl_float(17325.0,precision) + cl_float(31185.0,precision)*r*xj + cl_float(27720.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(16170.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(6930.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(2310.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(616.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(132.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(22.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               TWO*Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(126.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,EIGHT)*

             (cl_float(37125.0,precision) + cl_float(66825.0,precision)*r*xj + cl_float(59400.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(34725.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(14625.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(5043.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(1396.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(276.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(34.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               TWO*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(126.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(10.0,precision))*

             (cl_float(51975.0,precision) + cl_float(93420.0,precision)*r*xj + cl_float(84240.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(46815.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(20835.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(7485.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(1964.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(348.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(38.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               TWO*Power(r,NINE)*Power(xj,NINE)) - 

            NINE*Power(xi,TWO)*Power(xj,cl_float(16.0,precision))*

             (-cl_float(135135.0,precision) + cl_float(405405.0,precision)*r*xj + cl_float(582120.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(346500.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(124740.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(30492.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(5264.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(636.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(50.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               TWO*Power(r,NINE)*Power(xj,NINE)) + 

            Power(xj,cl_float(18.0,precision))*(cl_float(2837835.0,precision) + cl_float(3648645.0,precision)*r*xj + 

               cl_float(2245320.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(873180.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(238140.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(47628.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(7056.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(756.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(54.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               TWO*Power(r,NINE)*Power(xj,NINE)) - 

            NINE*Power(xi,cl_float(14.0,precision))*Power(xj,FOUR)*

             (cl_float(86625.0,precision) + cl_float(155925.0,precision)*r*xj + cl_float(138600.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(80850.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(34650.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(11550.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(3080.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(672.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(104.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               EIGHT*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(21.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,SIX)*

             (cl_float(111375.0,precision) + cl_float(200475.0,precision)*r*xj + cl_float(178200.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(103950.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(44550.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(14778.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(4056.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(864.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(120.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               EIGHT*Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(21.0,precision)*Power(xi,SIX)*Power(xj,cl_float(12.0,precision))*

             (cl_float(307125.0,precision) + cl_float(594945.0,precision)*r*xj + cl_float(456840.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(281790.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(137430.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(47250.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(11064.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(1728.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(168.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               EIGHT*Power(r,NINE)*Power(xj,NINE)) + 

            NINE*Power(xi,FOUR)*Power(xj,cl_float(14.0,precision))*

             (cl_float(675675.0,precision) + cl_float(675675.0,precision)*r*xj + cl_float(748440.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(561330.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(256410.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(76230.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(15400.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(2112.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(184.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               EIGHT*Power(r,NINE)*Power(xj,NINE))))/

       (cl_float(14175.0,precision)*exp(TWO*r*(xi + xj))*Power(r,TWO)*Power(xi - xj,cl_float(11.0,precision))*

         Power(xi + xj,cl_float(11.0,precision))) + (TWO*(cl_float(14175.0,precision)*exp(TWO*r*(xi + xj))*

            Power(Power(xi,TWO) - Power(xj,TWO),cl_float(11.0,precision)) + 

           cl_float(2835.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(12.0,precision))*

            (-cl_float(35.0,precision)*Power(xi,cl_float(10.0,precision)) - FIVE*r*Power(xi,cl_float(11.0,precision)) - 

              cl_float(495.0,precision)*Power(xi,EIGHT)*Power(xj,TWO) - cl_float(55.0,precision)*r*Power(xi,NINE)*Power(xj,TWO) - 

              cl_float(1254.0,precision)*Power(xi,SIX)*Power(xj,FOUR) - cl_float(66.0,precision)*r*Power(xi,SEVEN)*Power(xj,FOUR) - 

              cl_float(726.0,precision)*Power(xi,FOUR)*Power(xj,SIX) + cl_float(66.0,precision)*r*Power(xi,FIVE)*Power(xj,SIX) - 

              cl_float(55.0,precision)*Power(xi,TWO)*Power(xj,EIGHT) + cl_float(55.0,precision)*r*Power(xi,THREE)*Power(xj,EIGHT) + 

              FIVE*Power(xj,cl_float(10.0,precision)) + FIVE*r*xi*Power(xj,cl_float(10.0,precision))) + 

           exp(TWO*r*xi)*Power(xi,FOUR)*

            (-(Power(xi,cl_float(18.0,precision))*(cl_float(14175.0,precision) + cl_float(25515.0,precision)*r*xj + 

                   cl_float(22680.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                   cl_float(13230.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

                   cl_float(5670.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(1890.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

                   cl_float(504.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(108.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

                   cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + TWO*Power(r,NINE)*Power(xj,NINE))) + 

              NINE*Power(xi,cl_float(16.0,precision))*Power(xj,TWO)*

               (cl_float(17325.0,precision) + cl_float(31185.0,precision)*r*xj + cl_float(27720.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(16170.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(6930.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(2310.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(616.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

                 cl_float(132.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(22.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

                 TWO*Power(r,NINE)*Power(xj,NINE)) - 

              cl_float(126.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,EIGHT)*

               (cl_float(37125.0,precision) + cl_float(66825.0,precision)*r*xj + cl_float(59400.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(34725.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(14625.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(5043.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(1396.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

                 cl_float(276.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(34.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

                 TWO*Power(r,NINE)*Power(xj,NINE)) + 

              cl_float(126.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(10.0,precision))*

               (cl_float(51975.0,precision) + cl_float(93420.0,precision)*r*xj + cl_float(84240.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(46815.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(20835.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(7485.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(1964.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

                 cl_float(348.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(38.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

                 TWO*Power(r,NINE)*Power(xj,NINE)) - 

              NINE*Power(xi,TWO)*Power(xj,cl_float(16.0,precision))*

               (-cl_float(135135.0,precision) + cl_float(405405.0,precision)*r*xj + cl_float(582120.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(346500.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

                 cl_float(124740.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(30492.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(5264.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

                 cl_float(636.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(50.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

                 TWO*Power(r,NINE)*Power(xj,NINE)) + 

              Power(xj,cl_float(18.0,precision))*(cl_float(2837835.0,precision) + cl_float(3648645.0,precision)*r*xj + 

                 cl_float(2245320.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(873180.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

                 cl_float(238140.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(47628.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(7056.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

                 cl_float(756.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(54.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

                 TWO*Power(r,NINE)*Power(xj,NINE)) - 

              NINE*Power(xi,cl_float(14.0,precision))*Power(xj,FOUR)*

               (cl_float(86625.0,precision) + cl_float(155925.0,precision)*r*xj + cl_float(138600.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(80850.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(34650.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(11550.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(3080.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

                 cl_float(672.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(104.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

                 EIGHT*Power(r,NINE)*Power(xj,NINE)) + 

              cl_float(21.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,SIX)*

               (cl_float(111375.0,precision) + cl_float(200475.0,precision)*r*xj + cl_float(178200.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(103950.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

                 cl_float(44550.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(14778.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

                 cl_float(4056.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(864.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

                 cl_float(120.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + EIGHT*Power(r,NINE)*Power(xj,NINE)) - 

              cl_float(21.0,precision)*Power(xi,SIX)*Power(xj,cl_float(12.0,precision))*

               (cl_float(307125.0,precision) + cl_float(594945.0,precision)*r*xj + cl_float(456840.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(281790.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

                 cl_float(137430.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(47250.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(11064.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

                 cl_float(1728.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(168.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

                 EIGHT*Power(r,NINE)*Power(xj,NINE)) + 

              NINE*Power(xi,FOUR)*Power(xj,cl_float(14.0,precision))*

               (cl_float(675675.0,precision) + cl_float(675675.0,precision)*r*xj + cl_float(748440.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(561330.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

                 cl_float(256410.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(76230.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

                 cl_float(15400.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(2112.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

                 cl_float(184.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + EIGHT*Power(r,NINE)*Power(xj,NINE)))))/

       (cl_float(14175.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,cl_float(11.0,precision))*Power(xi + xj,cl_float(10.0,precision))) - 

      (cl_float(28350.0,precision)*exp(TWO*r*(xi + xj))*(xi + xj)*

          Power(Power(xi,TWO) - Power(xj,TWO),cl_float(11.0,precision)) + 

         cl_float(2835.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(12.0,precision))*

          (-FIVE*Power(xi,cl_float(11.0,precision)) - cl_float(55.0,precision)*Power(xi,NINE)*Power(xj,TWO) - 

            cl_float(66.0,precision)*Power(xi,SEVEN)*Power(xj,FOUR) + cl_float(66.0,precision)*Power(xi,FIVE)*Power(xj,SIX) + 

            cl_float(55.0,precision)*Power(xi,THREE)*Power(xj,EIGHT) + FIVE*xi*Power(xj,cl_float(10.0,precision))) + 

         cl_float(5670.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(13.0,precision))*

          (-cl_float(35.0,precision)*Power(xi,cl_float(10.0,precision)) - FIVE*r*Power(xi,cl_float(11.0,precision)) - cl_float(495.0,precision)*Power(xi,EIGHT)*Power(xj,TWO) - 

            cl_float(55.0,precision)*r*Power(xi,NINE)*Power(xj,TWO) - cl_float(1254.0,precision)*Power(xi,SIX)*Power(xj,FOUR) - 

            cl_float(66.0,precision)*r*Power(xi,SEVEN)*Power(xj,FOUR) - cl_float(726.0,precision)*Power(xi,FOUR)*Power(xj,SIX) + 

            cl_float(66.0,precision)*r*Power(xi,FIVE)*Power(xj,SIX) - cl_float(55.0,precision)*Power(xi,TWO)*Power(xj,EIGHT) + 

            cl_float(55.0,precision)*r*Power(xi,THREE)*Power(xj,EIGHT) + FIVE*Power(xj,cl_float(10.0,precision)) + FIVE*r*xi*Power(xj,cl_float(10.0,precision))) 

    + exp(TWO*r*xi)*Power(xi,FOUR)*(-(Power(xi,cl_float(18.0,precision))*

               (cl_float(25515.0,precision)*xj + cl_float(45360.0,precision)*r*Power(xj,TWO) + 

                 cl_float(39690.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(22680.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

                 cl_float(9450.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(3024.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

                 cl_float(756.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(144.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

                 cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,NINE))) + 

            NINE*Power(xi,cl_float(16.0,precision))*Power(xj,TWO)*

             (cl_float(31185.0,precision)*xj + cl_float(55440.0,precision)*r*Power(xj,TWO) + cl_float(48510.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(27720.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(11550.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(3696.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(924.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + 

               cl_float(176.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) - 

            cl_float(126.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,EIGHT)*

             (cl_float(66825.0,precision)*xj + cl_float(118800.0,precision)*r*Power(xj,TWO) + 

               cl_float(104175.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(58500.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(25215.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(8376.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(1932.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(272.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            cl_float(126.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(10.0,precision))*

             (cl_float(93420.0,precision)*xj + cl_float(168480.0,precision)*r*Power(xj,TWO) + 

               cl_float(140445.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(83340.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(37425.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(11784.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(2436.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(304.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) - 

            NINE*Power(xi,TWO)*Power(xj,cl_float(16.0,precision))*

             (cl_float(405405.0,precision)*xj + cl_float(1164240.0,precision)*r*Power(xj,TWO) + 

               cl_float(1039500.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(498960.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(152460.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(31584.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(4452.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + 

               cl_float(400.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            Power(xj,cl_float(18.0,precision))*(cl_float(3648645.0,precision)*xj + cl_float(4490640.0,precision)*r*Power(xj,TWO) + 

               cl_float(2619540.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(952560.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(238140.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(42336.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(5292.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + 

               cl_float(432.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) - 

            NINE*Power(xi,cl_float(14.0,precision))*Power(xj,FOUR)*

             (cl_float(155925.0,precision)*xj + cl_float(277200.0,precision)*r*Power(xj,TWO) + 

               cl_float(242550.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(138600.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(57750.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(18480.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(4704.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(832.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(72.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            cl_float(21.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,SIX)*

             (cl_float(200475.0,precision)*xj + cl_float(356400.0,precision)*r*Power(xj,TWO) + 

               cl_float(311850.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(178200.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(73890.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(24336.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(6048.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(960.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(72.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) - 

            cl_float(21.0,precision)*Power(xi,SIX)*Power(xj,cl_float(12.0,precision))*

             (cl_float(594945.0,precision)*xj + cl_float(913680.0,precision)*r*Power(xj,TWO) + 

               cl_float(845370.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(549720.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(236250.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(66384.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(12096.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(1344.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(72.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            NINE*Power(xi,FOUR)*Power(xj,cl_float(14.0,precision))*

             (cl_float(675675.0,precision)*xj + cl_float(1496880.0,precision)*r*Power(xj,TWO) + 

               cl_float(1683990.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(1025640.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(381150.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(92400.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(14784.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + 

               cl_float(1472.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + cl_float(72.0,precision)*Power(r,EIGHT)*Power(xj,NINE))) + 

         TWO*exp(TWO*r*xi)*Power(xi,FIVE)*

          (-(Power(xi,cl_float(18.0,precision))*(cl_float(14175.0,precision) + cl_float(25515.0,precision)*r*xj + cl_float(22680.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(13230.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(5670.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(1890.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(504.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

                 cl_float(108.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

                 TWO*Power(r,NINE)*Power(xj,NINE))) + 

            NINE*Power(xi,cl_float(16.0,precision))*Power(xj,TWO)*

             (cl_float(17325.0,precision) + cl_float(31185.0,precision)*r*xj + cl_float(27720.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(16170.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(6930.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(2310.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(616.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(132.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(22.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               TWO*Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(126.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,EIGHT)*

             (cl_float(37125.0,precision) + cl_float(66825.0,precision)*r*xj + cl_float(59400.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(34725.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(14625.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(5043.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(1396.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(276.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(34.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               TWO*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(126.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(10.0,precision))*

             (cl_float(51975.0,precision) + cl_float(93420.0,precision)*r*xj + cl_float(84240.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(46815.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(20835.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(7485.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(1964.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(348.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(38.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               TWO*Power(r,NINE)*Power(xj,NINE)) - 

            NINE*Power(xi,TWO)*Power(xj,cl_float(16.0,precision))*

             (-cl_float(135135.0,precision) + cl_float(405405.0,precision)*r*xj + cl_float(582120.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(346500.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(124740.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(30492.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(5264.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(636.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(50.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               TWO*Power(r,NINE)*Power(xj,NINE)) + 

            Power(xj,cl_float(18.0,precision))*(cl_float(2837835.0,precision) + cl_float(3648645.0,precision)*r*xj + 

               cl_float(2245320.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(873180.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(238140.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(47628.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(7056.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(756.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(54.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + TWO*Power(r,NINE)*Power(xj,NINE)) - 

            NINE*Power(xi,cl_float(14.0,precision))*Power(xj,FOUR)*

             (cl_float(86625.0,precision) + cl_float(155925.0,precision)*r*xj + cl_float(138600.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(80850.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(34650.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(11550.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(3080.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(672.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(104.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               EIGHT*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(21.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,SIX)*

             (cl_float(111375.0,precision) + cl_float(200475.0,precision)*r*xj + cl_float(178200.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(103950.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(44550.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(14778.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(4056.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(864.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(120.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               EIGHT*Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(21.0,precision)*Power(xi,SIX)*Power(xj,cl_float(12.0,precision))*

             (cl_float(307125.0,precision) + cl_float(594945.0,precision)*r*xj + cl_float(456840.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(281790.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(137430.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(47250.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(11064.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(1728.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(168.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               EIGHT*Power(r,NINE)*Power(xj,NINE)) + 

            NINE*Power(xi,FOUR)*Power(xj,cl_float(14.0,precision))*

             (cl_float(675675.0,precision) + cl_float(675675.0,precision)*r*xj + cl_float(748440.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(561330.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(256410.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(76230.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(15400.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(2112.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(184.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               EIGHT*Power(r,NINE)*Power(xj,NINE))))/

       (cl_float(14175.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,cl_float(11.0,precision))*Power(xi + xj,cl_float(11.0,precision)))

    ; }
   
  }
  return S;
}

cl_F DSlater_2S_2S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = ZERO

    ; } else {  S = -(-cl_float(131985.0,precision)*xi + cl_float(161280.0,precision)*exp(TWO*r*xi)*xi - cl_float(205380.0,precision)*r*Power(xi,TWO) - 

          cl_float(149940.0,precision)*Power(r,TWO)*Power(xi,THREE) - cl_float(67200.0,precision)*Power(r,THREE)*Power(xi,FOUR) - 

          cl_float(20160.0,precision)*Power(r,FOUR)*Power(xi,FIVE) - cl_float(4032.0,precision)*Power(r,FIVE)*Power(xi,SIX) - 

          cl_float(448.0,precision)*Power(r,SIX)*Power(xi,SEVEN))/(cl_float(80640.0,precision)*exp(TWO*r*xi)*r) + 

      (-cl_float(80640.0,precision) + cl_float(80640.0,precision)*exp(TWO*r*xi) - cl_float(131985.0,precision)*r*xi - 

         cl_float(102690.0,precision)*Power(r,TWO)*Power(xi,TWO) - cl_float(49980.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

         cl_float(16800.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - cl_float(4032.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

         cl_float(672.0,precision)*Power(r,SIX)*Power(xi,SIX) - cl_float(64.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN))/

       (cl_float(80640.0,precision)*exp(TWO*r*xi)*Power(r,TWO)) + 

      (xi*(-cl_float(80640.0,precision) + cl_float(80640.0,precision)*exp(TWO*r*xi) - cl_float(131985.0,precision)*r*xi - 

           cl_float(102690.0,precision)*Power(r,TWO)*Power(xi,TWO) - cl_float(49980.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

           cl_float(16800.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - cl_float(4032.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

           cl_float(672.0,precision)*Power(r,SIX)*Power(xi,SIX) - cl_float(64.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN)))/

       (cl_float(40320.0,precision)*exp(TWO*r*xi)*r)

    ; }
 
  }
  else {
      if (r == ZERO) {  S = ZERO

    ; } else {  S = (SIX*exp(TWO*r*(xi + xj))*Power(Power(xi,TWO) - Power(xj,TWO),SEVEN) - 

         exp(TWO*r*xi)*Power(xi,SIX)*

          (cl_float(21.0,precision)*Power(xi,FOUR)*Power(xj,FOUR)*

             (SIX + cl_float(11.0,precision)*r*xj + TWO*Power(r,TWO)*Power(xj,TWO)) - 

            TWO*Power(xj,EIGHT)*(cl_float(90.0,precision) + cl_float(54.0,precision)*r*xj + cl_float(12.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               Power(r,THREE)*Power(xj,THREE)) + 

            Power(xi,EIGHT)*(SIX + NINE*r*xj + SIX*Power(r,TWO)*Power(xj,TWO) + 

               TWO*Power(r,THREE)*Power(xj,THREE)) + 

            Power(xi,TWO)*Power(xj,SIX)*

             (-cl_float(390.0,precision) - cl_float(69.0,precision)*r*xj + cl_float(18.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               FOUR*Power(r,THREE)*Power(xj,THREE)) - 

            Power(xi,SIX)*Power(xj,TWO)*

             (cl_float(42.0,precision) + cl_float(63.0,precision)*r*xj + cl_float(42.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               FOUR*Power(r,THREE)*Power(xj,THREE))) + 

         exp(TWO*r*xj)*Power(xj,SIX)*

          (-cl_float(24.0,precision)*Power(r,TWO)*Power(xi,cl_float(10.0,precision)) - TWO*Power(r,THREE)*Power(xi,cl_float(11.0,precision)) - 

            cl_float(69.0,precision)*r*Power(xi,SEVEN)*Power(xj,TWO) + SIX*Power(xj,EIGHT) + NINE*r*xi*Power(xj,EIGHT) + 

            FOUR*r*Power(xi,NINE)*(-cl_float(27.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(18.0,precision)*Power(xi,EIGHT)*(-cl_float(10.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            SIX*Power(xi,TWO)*Power(xj,SIX)*(-SEVEN + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(42.0,precision)*Power(xi,FOUR)*Power(xj,FOUR)*(-THREE + Power(r,TWO)*Power(xj,TWO)) + 

            r*Power(xi,THREE)*Power(xj,SIX)*(-cl_float(63.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            SIX*Power(xi,SIX)*Power(xj,TWO)*(-cl_float(65.0,precision) + SEVEN*Power(r,TWO)*Power(xj,TWO)) + 

            Power(xi,FIVE)*(cl_float(231.0,precision)*r*Power(xj,FOUR) - FOUR*Power(r,THREE)*Power(xj,SIX))))/

       (cl_float(6.0,precision)*exp(TWO*r*(xi + xj))*Power(r,TWO)*Power(xi - xj,SEVEN)*Power(xi + xj,SEVEN)) + 

      (SIX*exp(TWO*r*(xi + xj))*Power(Power(xi,TWO) - Power(xj,TWO),SEVEN) - 

         exp(TWO*r*xi)*Power(xi,SIX)*

          (cl_float(21.0,precision)*Power(xi,FOUR)*Power(xj,FOUR)*

             (SIX + cl_float(11.0,precision)*r*xj + TWO*Power(r,TWO)*Power(xj,TWO)) - 

            TWO*Power(xj,EIGHT)*(cl_float(90.0,precision) + cl_float(54.0,precision)*r*xj + cl_float(12.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               Power(r,THREE)*Power(xj,THREE)) + 

            Power(xi,EIGHT)*(SIX + NINE*r*xj + SIX*Power(r,TWO)*Power(xj,TWO) + 

               TWO*Power(r,THREE)*Power(xj,THREE)) + 

            Power(xi,TWO)*Power(xj,SIX)*

             (-cl_float(390.0,precision) - cl_float(69.0,precision)*r*xj + cl_float(18.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               FOUR*Power(r,THREE)*Power(xj,THREE)) - 

            Power(xi,SIX)*Power(xj,TWO)*

             (cl_float(42.0,precision) + cl_float(63.0,precision)*r*xj + cl_float(42.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               FOUR*Power(r,THREE)*Power(xj,THREE))) + 

         exp(TWO*r*xj)*Power(xj,SIX)*

          (-cl_float(24.0,precision)*Power(r,TWO)*Power(xi,cl_float(10.0,precision)) - TWO*Power(r,THREE)*Power(xi,cl_float(11.0,precision)) - 

            cl_float(69.0,precision)*r*Power(xi,SEVEN)*Power(xj,TWO) + SIX*Power(xj,EIGHT) + NINE*r*xi*Power(xj,EIGHT) + 

            FOUR*r*Power(xi,NINE)*(-cl_float(27.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(18.0,precision)*Power(xi,EIGHT)*(-cl_float(10.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            SIX*Power(xi,TWO)*Power(xj,SIX)*(-SEVEN + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(42.0,precision)*Power(xi,FOUR)*Power(xj,FOUR)*(-THREE + Power(r,TWO)*Power(xj,TWO)) + 

            r*Power(xi,THREE)*Power(xj,SIX)*(-cl_float(63.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            SIX*Power(xi,SIX)*Power(xj,TWO)*(-cl_float(65.0,precision) + SEVEN*Power(r,TWO)*Power(xj,TWO)) + 

            Power(xi,FIVE)*(cl_float(231.0,precision)*r*Power(xj,FOUR) - FOUR*Power(r,THREE)*Power(xj,SIX))))/

       (cl_float(3.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,SEVEN)*Power(xi + xj,SIX)) - 

      (cl_float(12.0,precision)*exp(TWO*r*(xi + xj))*(xi + xj)*Power(Power(xi,TWO) - Power(xj,TWO),SEVEN) - 

         exp(TWO*r*xi)*Power(xi,SIX)*

          (cl_float(21.0,precision)*Power(xi,FOUR)*Power(xj,FOUR)*(cl_float(11.0,precision)*xj + FOUR*r*Power(xj,TWO)) - 

            TWO*Power(xj,EIGHT)*(cl_float(54.0,precision)*xj + cl_float(24.0,precision)*r*Power(xj,TWO) + 

               THREE*Power(r,TWO)*Power(xj,THREE)) + 

            Power(xi,EIGHT)*(NINE*xj + cl_float(12.0,precision)*r*Power(xj,TWO) + SIX*Power(r,TWO)*Power(xj,THREE)) + 

            Power(xi,TWO)*Power(xj,SIX)*

             (-cl_float(69.0,precision)*xj + cl_float(36.0,precision)*r*Power(xj,TWO) + cl_float(12.0,precision)*Power(r,TWO)*Power(xj,THREE)) - 

            Power(xi,SIX)*Power(xj,TWO)*

             (cl_float(63.0,precision)*xj + cl_float(84.0,precision)*r*Power(xj,TWO) + cl_float(12.0,precision)*Power(r,TWO)*Power(xj,THREE))) - 

         TWO*exp(TWO*r*xi)*Power(xi,SEVEN)*

          (cl_float(21.0,precision)*Power(xi,FOUR)*Power(xj,FOUR)*

             (SIX + cl_float(11.0,precision)*r*xj + TWO*Power(r,TWO)*Power(xj,TWO)) - 

            TWO*Power(xj,EIGHT)*(cl_float(90.0,precision) + cl_float(54.0,precision)*r*xj + cl_float(12.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               Power(r,THREE)*Power(xj,THREE)) + 

            Power(xi,EIGHT)*(SIX + NINE*r*xj + SIX*Power(r,TWO)*Power(xj,TWO) + 

               TWO*Power(r,THREE)*Power(xj,THREE)) + 

            Power(xi,TWO)*Power(xj,SIX)*

             (-cl_float(390.0,precision) - cl_float(69.0,precision)*r*xj + cl_float(18.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               FOUR*Power(r,THREE)*Power(xj,THREE)) - 

            Power(xi,SIX)*Power(xj,TWO)*

             (cl_float(42.0,precision) + cl_float(63.0,precision)*r*xj + cl_float(42.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               FOUR*Power(r,THREE)*Power(xj,THREE))) + 

         exp(TWO*r*xj)*Power(xj,SIX)*

          (-cl_float(48.0,precision)*r*Power(xi,cl_float(10.0,precision)) - SIX*Power(r,TWO)*Power(xi,cl_float(11.0,precision)) - 

            cl_float(69.0,precision)*Power(xi,SEVEN)*Power(xj,TWO) + cl_float(36.0,precision)*r*Power(xi,EIGHT)*Power(xj,TWO) + 

            EIGHT*Power(r,TWO)*Power(xi,NINE)*Power(xj,TWO) + 

            cl_float(84.0,precision)*r*Power(xi,SIX)*Power(xj,FOUR) - cl_float(84.0,precision)*r*Power(xi,FOUR)*Power(xj,SIX) + 

            NINE*xi*Power(xj,EIGHT) + cl_float(12.0,precision)*r*Power(xi,TWO)*Power(xj,EIGHT) + 

            FOUR*Power(r,TWO)*Power(xi,THREE)*Power(xj,EIGHT) + 

            FOUR*Power(xi,NINE)*(-cl_float(27.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            Power(xi,THREE)*Power(xj,SIX)*(-cl_float(63.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            Power(xi,FIVE)*(cl_float(231.0,precision)*Power(xj,FOUR) - cl_float(12.0,precision)*Power(r,TWO)*Power(xj,SIX))) + 

         TWO*exp(TWO*r*xj)*Power(xj,SEVEN)*

          (-cl_float(24.0,precision)*Power(r,TWO)*Power(xi,cl_float(10.0,precision)) - TWO*Power(r,THREE)*Power(xi,cl_float(11.0,precision)) - 

            cl_float(69.0,precision)*r*Power(xi,SEVEN)*Power(xj,TWO) + SIX*Power(xj,EIGHT) + NINE*r*xi*Power(xj,EIGHT) + 

            FOUR*r*Power(xi,NINE)*(-cl_float(27.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(18.0,precision)*Power(xi,EIGHT)*(-cl_float(10.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            SIX*Power(xi,TWO)*Power(xj,SIX)*(-SEVEN + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(42.0,precision)*Power(xi,FOUR)*Power(xj,FOUR)*(-THREE + Power(r,TWO)*Power(xj,TWO)) + 

            r*Power(xi,THREE)*Power(xj,SIX)*(-cl_float(63.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            SIX*Power(xi,SIX)*Power(xj,TWO)*(-cl_float(65.0,precision) + SEVEN*Power(r,TWO)*Power(xj,TWO)) + 

            Power(xi,FIVE)*(cl_float(231.0,precision)*r*Power(xj,FOUR) - FOUR*Power(r,THREE)*Power(xj,SIX))))/

       (cl_float(6.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,SEVEN)*Power(xi + xj,SEVEN))

    ; }
   
  }
  return S;
}

cl_F DSlater_2S_3S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = ZERO

    ; } else {  S = -(-cl_float(7430535.0,precision)*xi + cl_float(8709120.0,precision)*exp(TWO*r*xi)*xi - cl_float(12303900.0,precision)*r*Power(xi,TWO) - 

          cl_float(9826110.0,precision)*Power(r,TWO)*Power(xi,THREE) - cl_float(5004720.0,precision)*Power(r,THREE)*Power(xi,FOUR) - 

          cl_float(1806840.0,precision)*Power(r,FOUR)*Power(xi,FIVE) - cl_float(483840.0,precision)*Power(r,FIVE)*Power(xi,SIX) - 

          cl_float(96768.0,precision)*Power(r,SIX)*Power(xi,SEVEN) - cl_float(13824.0,precision)*Power(r,SEVEN)*Power(xi,EIGHT) - 

          cl_float(1152.0,precision)*Power(r,EIGHT)*Power(xi,NINE))/(cl_float(4.35456e6,precision)*exp(TWO*r*xi)*r) + 

      (-cl_float(4354560.0,precision) + cl_float(4354560.0,precision)*exp(TWO*r*xi) - cl_float(7430535.0,precision)*r*xi - 

         cl_float(6151950.0,precision)*Power(r,TWO)*Power(xi,TWO) - cl_float(3275370.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

         cl_float(1251180.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - cl_float(361368.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

         cl_float(80640.0,precision)*Power(r,SIX)*Power(xi,SIX) - cl_float(13824.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

         cl_float(1728.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - cl_float(128.0,precision)*Power(r,NINE)*Power(xi,NINE))/

       (cl_float(4.35456e6,precision)*exp(TWO*r*xi)*Power(r,TWO)) + 

      (xi*(-cl_float(4354560.0,precision) + cl_float(4354560.0,precision)*exp(TWO*r*xi) - cl_float(7430535.0,precision)*r*xi - 

           cl_float(6151950.0,precision)*Power(r,TWO)*Power(xi,TWO) - cl_float(3275370.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

           cl_float(1251180.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - cl_float(361368.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

           cl_float(80640.0,precision)*Power(r,SIX)*Power(xi,SIX) - cl_float(13824.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

           cl_float(1728.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - cl_float(128.0,precision)*Power(r,NINE)*Power(xi,NINE)))/

       (cl_float(2.17728e6,precision)*exp(TWO*r*xi)*r)

    ; }
 
  }
  else {
      if (r == ZERO) {  S = ZERO

    ; } else {  S = (cl_float(90.0,precision)*exp(TWO*r*(xi + xj))*Power(Power(xi,TWO) - Power(xj,TWO),NINE) + 

         FIVE*exp(TWO*r*xj)*Power(xj,EIGHT)*

          (-cl_float(90.0,precision)*Power(r,TWO)*Power(xi,cl_float(12.0,precision)) - SIX*Power(r,THREE)*Power(xi,cl_float(13.0,precision)) + 

            cl_float(18.0,precision)*Power(xj,cl_float(10.0,precision)) + cl_float(27.0,precision)*r*xi*Power(xj,cl_float(10.0,precision)) + 

            cl_float(18.0,precision)*Power(xi,TWO)*Power(xj,EIGHT)*(-NINE + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(162.0,precision)*Power(xi,FOUR)*Power(xj,SIX)*(-FOUR + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(198.0,precision)*Power(xi,cl_float(10.0,precision))*(FIVE + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(108.0,precision)*Power(xi,SIX)*Power(xj,FOUR)*(cl_float(36.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            TWO*r*Power(xi,FIVE)*Power(xj,SIX)*(cl_float(675.0,precision) + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(18.0,precision)*r*Power(xi,SEVEN)*Power(xj,FOUR)*(-cl_float(81.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            THREE*r*Power(xi,THREE)*Power(xj,EIGHT)*(-cl_float(81.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) - 

            r*Power(xi,cl_float(11.0,precision))*(cl_float(495.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            NINE*r*Power(xi,NINE)*Power(xj,TWO)*(-cl_float(233.0,precision) + FOUR*Power(r,TWO)*Power(xj,TWO)) + 

            SIX*Power(xi,EIGHT)*Power(xj,TWO)*(-cl_float(1063.0,precision) + cl_float(90.0,precision)*Power(r,TWO)*Power(xj,TWO))) - 

         TWO*exp(TWO*r*xi)*Power(xi,SIX)*

          (-cl_float(90.0,precision)*Power(xi,SIX)*Power(xj,SIX)*

             (cl_float(42.0,precision) + cl_float(65.0,precision)*r*xj + cl_float(76.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(22.0,precision)*Power(r,THREE)*Power(xj,THREE) + TWO*Power(r,FOUR)*Power(xj,FOUR)) - 

            TWO*Power(xj,cl_float(12.0,precision))*(cl_float(2970.0,precision) + cl_float(2475.0,precision)*r*xj + cl_float(900.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(180.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(20.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               Power(r,FIVE)*Power(xj,FIVE)) + 

            cl_float(10.0,precision)*Power(xi,EIGHT)*Power(xj,FOUR)*

             (cl_float(162.0,precision) + cl_float(270.0,precision)*r*xj + cl_float(216.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(122.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(22.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               Power(r,FIVE)*Power(xj,FIVE)) - 

            FIVE*Power(xi,FOUR)*Power(xj,EIGHT)*

             (-cl_float(639.0,precision) - cl_float(3555.0,precision)*r*xj - cl_float(1452.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(174.0,precision)*Power(r,THREE)*Power(xj,THREE) + SIX*Power(r,FOUR)*Power(xj,FOUR) + 

               TWO*Power(r,FIVE)*Power(xj,FIVE)) + 

            Power(xi,cl_float(12.0,precision))*(cl_float(45.0,precision) + cl_float(75.0,precision)*r*xj + cl_float(60.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(30.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(10.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               TWO*Power(r,FIVE)*Power(xj,FIVE)) - 

            Power(xi,cl_float(10.0,precision))*Power(xj,TWO)*

             (cl_float(405.0,precision) + cl_float(675.0,precision)*r*xj + cl_float(540.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(270.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(90.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               EIGHT*Power(r,FIVE)*Power(xj,FIVE)) + 

            Power(xi,TWO)*Power(xj,cl_float(10.0,precision))*

             (-cl_float(21615.0,precision) - cl_float(9075.0,precision)*r*xj - cl_float(300.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(490.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(110.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               EIGHT*Power(r,FIVE)*Power(xj,FIVE))))/

       (cl_float(90.0,precision)*exp(TWO*r*(xi + xj))*Power(r,TWO)*Power(xi - xj,NINE)*Power(xi + xj,NINE)) 

    + (cl_float(90.0,precision)*exp(TWO*r*(xi + xj))*Power(Power(xi,TWO) - Power(xj,TWO),NINE) + 

         FIVE*exp(TWO*r*xj)*Power(xj,EIGHT)*

          (-cl_float(90.0,precision)*Power(r,TWO)*Power(xi,cl_float(12.0,precision)) - SIX*Power(r,THREE)*Power(xi,cl_float(13.0,precision)) + 

            cl_float(18.0,precision)*Power(xj,cl_float(10.0,precision)) + cl_float(27.0,precision)*r*xi*Power(xj,cl_float(10.0,precision)) + 

            cl_float(18.0,precision)*Power(xi,TWO)*Power(xj,EIGHT)*(-NINE + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(162.0,precision)*Power(xi,FOUR)*Power(xj,SIX)*(-FOUR + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(198.0,precision)*Power(xi,cl_float(10.0,precision))*(FIVE + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(108.0,precision)*Power(xi,SIX)*Power(xj,FOUR)*(cl_float(36.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            TWO*r*Power(xi,FIVE)*Power(xj,SIX)*(cl_float(675.0,precision) + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(18.0,precision)*r*Power(xi,SEVEN)*Power(xj,FOUR)*(-cl_float(81.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            THREE*r*Power(xi,THREE)*Power(xj,EIGHT)*(-cl_float(81.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) - 

            r*Power(xi,cl_float(11.0,precision))*(cl_float(495.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            NINE*r*Power(xi,NINE)*Power(xj,TWO)*(-cl_float(233.0,precision) + FOUR*Power(r,TWO)*Power(xj,TWO)) + 

            SIX*Power(xi,EIGHT)*Power(xj,TWO)*(-cl_float(1063.0,precision) + cl_float(90.0,precision)*Power(r,TWO)*Power(xj,TWO))) - 

         TWO*exp(TWO*r*xi)*Power(xi,SIX)*

          (-cl_float(90.0,precision)*Power(xi,SIX)*Power(xj,SIX)*

             (cl_float(42.0,precision) + cl_float(65.0,precision)*r*xj + cl_float(76.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(22.0,precision)*Power(r,THREE)*Power(xj,THREE) + TWO*Power(r,FOUR)*Power(xj,FOUR)) - 

            TWO*Power(xj,cl_float(12.0,precision))*(cl_float(2970.0,precision) + cl_float(2475.0,precision)*r*xj + cl_float(900.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(180.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(20.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               Power(r,FIVE)*Power(xj,FIVE)) + 

            cl_float(10.0,precision)*Power(xi,EIGHT)*Power(xj,FOUR)*

             (cl_float(162.0,precision) + cl_float(270.0,precision)*r*xj + cl_float(216.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(122.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(22.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               Power(r,FIVE)*Power(xj,FIVE)) - 

            FIVE*Power(xi,FOUR)*Power(xj,EIGHT)*

             (-cl_float(639.0,precision) - cl_float(3555.0,precision)*r*xj - cl_float(1452.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(174.0,precision)*Power(r,THREE)*Power(xj,THREE) + SIX*Power(r,FOUR)*Power(xj,FOUR) + 

               TWO*Power(r,FIVE)*Power(xj,FIVE)) + 

            Power(xi,cl_float(12.0,precision))*(cl_float(45.0,precision) + cl_float(75.0,precision)*r*xj + cl_float(60.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(30.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(10.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               TWO*Power(r,FIVE)*Power(xj,FIVE)) - 

            Power(xi,cl_float(10.0,precision))*Power(xj,TWO)*

             (cl_float(405.0,precision) + cl_float(675.0,precision)*r*xj + cl_float(540.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(270.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(90.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               EIGHT*Power(r,FIVE)*Power(xj,FIVE)) + 

            Power(xi,TWO)*Power(xj,cl_float(10.0,precision))*

             (-cl_float(21615.0,precision) - cl_float(9075.0,precision)*r*xj - cl_float(300.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(490.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(110.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               EIGHT*Power(r,FIVE)*Power(xj,FIVE))))/

       (cl_float(45.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,NINE)*Power(xi + xj,EIGHT)) - 

      (cl_float(180.0,precision)*exp(TWO*r*(xi + xj))*(xi + xj)*Power(Power(xi,TWO) - Power(xj,TWO),NINE) + 

         FIVE*exp(TWO*r*xj)*Power(xj,EIGHT)*

          (-cl_float(180.0,precision)*r*Power(xi,cl_float(12.0,precision)) - cl_float(18.0,precision)*Power(r,TWO)*Power(xi,cl_float(13.0,precision)) - 

            cl_float(396.0,precision)*r*Power(xi,cl_float(10.0,precision))*Power(xj,TWO) - 

            FOUR*Power(r,TWO)*Power(xi,cl_float(11.0,precision))*Power(xj,TWO) + 

            cl_float(1080.0,precision)*r*Power(xi,EIGHT)*Power(xj,FOUR) + 

            cl_float(72.0,precision)*Power(r,TWO)*Power(xi,NINE)*Power(xj,FOUR) - 

            cl_float(216.0,precision)*r*Power(xi,SIX)*Power(xj,SIX) - 

            cl_float(72.0,precision)*Power(r,TWO)*Power(xi,SEVEN)*Power(xj,SIX) - 

            cl_float(324.0,precision)*r*Power(xi,FOUR)*Power(xj,EIGHT) + 

            FOUR*Power(r,TWO)*Power(xi,FIVE)*Power(xj,EIGHT) + cl_float(27.0,precision)*xi*Power(xj,cl_float(10.0,precision)) + 

            cl_float(36.0,precision)*r*Power(xi,TWO)*Power(xj,cl_float(10.0,precision)) + 

            cl_float(12.0,precision)*Power(r,TWO)*Power(xi,THREE)*Power(xj,cl_float(10.0,precision)) + 

            TWO*Power(xi,FIVE)*Power(xj,SIX)*(cl_float(675.0,precision) + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(18.0,precision)*Power(xi,SEVEN)*Power(xj,FOUR)*(-cl_float(81.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            THREE*Power(xi,THREE)*Power(xj,EIGHT)*(-cl_float(81.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) - 

            Power(xi,cl_float(11.0,precision))*(cl_float(495.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            NINE*Power(xi,NINE)*Power(xj,TWO)*(-cl_float(233.0,precision) + FOUR*Power(r,TWO)*Power(xj,TWO))) + 

         cl_float(10.0,precision)*exp(TWO*r*xj)*Power(xj,NINE)*

          (-cl_float(90.0,precision)*Power(r,TWO)*Power(xi,cl_float(12.0,precision)) - SIX*Power(r,THREE)*Power(xi,cl_float(13.0,precision)) + 

            cl_float(18.0,precision)*Power(xj,cl_float(10.0,precision)) + cl_float(27.0,precision)*r*xi*Power(xj,cl_float(10.0,precision)) + 

            cl_float(18.0,precision)*Power(xi,TWO)*Power(xj,EIGHT)*(-NINE + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(162.0,precision)*Power(xi,FOUR)*Power(xj,SIX)*(-FOUR + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(198.0,precision)*Power(xi,cl_float(10.0,precision))*(FIVE + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(108.0,precision)*Power(xi,SIX)*Power(xj,FOUR)*(cl_float(36.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            TWO*r*Power(xi,FIVE)*Power(xj,SIX)*(cl_float(675.0,precision) + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(18.0,precision)*r*Power(xi,SEVEN)*Power(xj,FOUR)*(-cl_float(81.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            THREE*r*Power(xi,THREE)*Power(xj,EIGHT)*(-cl_float(81.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) - 

            r*Power(xi,cl_float(11.0,precision))*(cl_float(495.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            NINE*r*Power(xi,NINE)*Power(xj,TWO)*(-cl_float(233.0,precision) + FOUR*Power(r,TWO)*Power(xj,TWO)) + 

            SIX*Power(xi,EIGHT)*Power(xj,TWO)*(-cl_float(1063.0,precision) + cl_float(90.0,precision)*Power(r,TWO)*Power(xj,TWO))) - 

         TWO*exp(TWO*r*xi)*Power(xi,SIX)*

          (-cl_float(90.0,precision)*Power(xi,SIX)*Power(xj,SIX)*

             (cl_float(65.0,precision)*xj + cl_float(152.0,precision)*r*Power(xj,TWO) + cl_float(66.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               EIGHT*Power(r,THREE)*Power(xj,FOUR)) - 

            TWO*Power(xj,cl_float(12.0,precision))*(cl_float(2475.0,precision)*xj + cl_float(1800.0,precision)*r*Power(xj,TWO) + 

               cl_float(540.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(80.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               FIVE*Power(r,FOUR)*Power(xj,FIVE)) + 

            cl_float(10.0,precision)*Power(xi,EIGHT)*Power(xj,FOUR)*

             (cl_float(270.0,precision)*xj + cl_float(432.0,precision)*r*Power(xj,TWO) + cl_float(366.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(88.0,precision)*Power(r,THREE)*Power(xj,FOUR) + FIVE*Power(r,FOUR)*Power(xj,FIVE)) - 

            FIVE*Power(xi,FOUR)*Power(xj,EIGHT)*

             (-cl_float(3555.0,precision)*xj - cl_float(2904.0,precision)*r*Power(xj,TWO) - cl_float(522.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(24.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(10.0,precision)*Power(r,FOUR)*Power(xj,FIVE)) + 

            Power(xi,cl_float(12.0,precision))*(cl_float(75.0,precision)*xj + cl_float(120.0,precision)*r*Power(xj,TWO) + 

               cl_float(90.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(40.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(10.0,precision)*Power(r,FOUR)*Power(xj,FIVE)) - 

            Power(xi,cl_float(10.0,precision))*Power(xj,TWO)*

             (cl_float(675.0,precision)*xj + cl_float(1080.0,precision)*r*Power(xj,TWO) + cl_float(810.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(360.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(40.0,precision)*Power(r,FOUR)*Power(xj,FIVE)) + 

            Power(xi,TWO)*Power(xj,cl_float(10.0,precision))*

             (-cl_float(9075.0,precision)*xj - cl_float(600.0,precision)*r*Power(xj,TWO) + cl_float(1470.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(440.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(40.0,precision)*Power(r,FOUR)*Power(xj,FIVE))) - 

         FOUR*exp(TWO*r*xi)*Power(xi,SEVEN)*

          (-cl_float(90.0,precision)*Power(xi,SIX)*Power(xj,SIX)*

             (cl_float(42.0,precision) + cl_float(65.0,precision)*r*xj + cl_float(76.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(22.0,precision)*Power(r,THREE)*Power(xj,THREE) + TWO*Power(r,FOUR)*Power(xj,FOUR)) - 

            TWO*Power(xj,cl_float(12.0,precision))*(cl_float(2970.0,precision) + cl_float(2475.0,precision)*r*xj + cl_float(900.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(180.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(20.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               Power(r,FIVE)*Power(xj,FIVE)) + 

            cl_float(10.0,precision)*Power(xi,EIGHT)*Power(xj,FOUR)*

             (cl_float(162.0,precision) + cl_float(270.0,precision)*r*xj + cl_float(216.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(122.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(22.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               Power(r,FIVE)*Power(xj,FIVE)) - 

            FIVE*Power(xi,FOUR)*Power(xj,EIGHT)*

             (-cl_float(639.0,precision) - cl_float(3555.0,precision)*r*xj - cl_float(1452.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(174.0,precision)*Power(r,THREE)*Power(xj,THREE) + SIX*Power(r,FOUR)*Power(xj,FOUR) + 

               TWO*Power(r,FIVE)*Power(xj,FIVE)) + 

            Power(xi,cl_float(12.0,precision))*(cl_float(45.0,precision) + cl_float(75.0,precision)*r*xj + cl_float(60.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(30.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(10.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               TWO*Power(r,FIVE)*Power(xj,FIVE)) - 

            Power(xi,cl_float(10.0,precision))*Power(xj,TWO)*

             (cl_float(405.0,precision) + cl_float(675.0,precision)*r*xj + cl_float(540.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(270.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(90.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               EIGHT*Power(r,FIVE)*Power(xj,FIVE)) + 

            Power(xi,TWO)*Power(xj,cl_float(10.0,precision))*

             (-cl_float(21615.0,precision) - cl_float(9075.0,precision)*r*xj - cl_float(300.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(490.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(110.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               EIGHT*Power(r,FIVE)*Power(xj,FIVE))))/

       (cl_float(90.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,NINE)*Power(xi + xj,NINE))

    ; }
   
  }
  return S;
}

cl_F DSlater_2S_4S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = ZERO

    ; } else {  S = -(-cl_float(1125310725.0,precision)*xi + cl_float(1277337600.0,precision)*exp(TWO*r*xi)*xi - 

          cl_float(1946567700.0,precision)*r*Power(xi,TWO) - cl_float(1647191700.0,precision)*Power(r,TWO)*Power(xi,THREE) - 

          cl_float(904780800.0,precision)*Power(r,THREE)*Power(xi,FOUR) - cl_float(360498600.0,precision)*Power(r,FOUR)*Power(xi,FIVE) - 

          cl_float(110103840.0,precision)*Power(r,FIVE)*Power(xi,SIX) - cl_float(26500320.0,precision)*Power(r,SIX)*Power(xi,SEVEN) - 

          cl_float(5068800.0,precision)*Power(r,SEVEN)*Power(xi,EIGHT) - cl_float(760320.0,precision)*Power(r,EIGHT)*Power(xi,NINE) - 

          cl_float(84480.0,precision)*Power(r,NINE)*Power(xi,cl_float(10.0,precision)) - cl_float(5632.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(11.0,precision)))/

       (cl_float(6.386688e8,precision)*exp(TWO*r*xi)*r) + 

      (-cl_float(638668800.0,precision) + cl_float(638668800.0,precision)*exp(TWO*r*xi) - cl_float(1125310725.0,precision)*r*xi - 

         cl_float(973283850.0,precision)*Power(r,TWO)*Power(xi,TWO) - cl_float(549063900.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

         cl_float(226195200.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - cl_float(72099720.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

         cl_float(18350640.0,precision)*Power(r,SIX)*Power(xi,SIX) - cl_float(3785760.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

         cl_float(633600.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - cl_float(84480.0,precision)*Power(r,NINE)*Power(xi,NINE) - 

         cl_float(8448.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - cl_float(512.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)))/

       (cl_float(6.386688e8,precision)*exp(TWO*r*xi)*Power(r,TWO)) + 

      (xi*(-cl_float(638668800.0,precision) + cl_float(638668800.0,precision)*exp(TWO*r*xi) - cl_float(1125310725.0,precision)*r*xi - 

           cl_float(973283850.0,precision)*Power(r,TWO)*Power(xi,TWO) - cl_float(549063900.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

           cl_float(226195200.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - cl_float(72099720.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

           cl_float(18350640.0,precision)*Power(r,SIX)*Power(xi,SIX) - cl_float(3785760.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

           cl_float(633600.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - cl_float(84480.0,precision)*Power(r,NINE)*Power(xi,NINE) - 

           cl_float(8448.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - cl_float(512.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision))))/

       (cl_float(3.193344e8,precision)*exp(TWO*r*xi)*r)

    ; }
 
  }
  else {
      if (r == ZERO) {  S = ZERO

    ; } else {  S = (cl_float(1260.0,precision)*exp(TWO*r*(xi + xj))*Power(Power(xi,TWO) - Power(xj,TWO),cl_float(11.0,precision)) + 

         cl_float(210.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(10.0,precision))*

          (-cl_float(36.0,precision)*Power(r,TWO)*Power(xi,cl_float(14.0,precision)) - TWO*Power(r,THREE)*Power(xi,cl_float(15.0,precision)) - 

            cl_float(1287.0,precision)*r*Power(xi,NINE)*Power(xj,FOUR) + SIX*Power(xj,cl_float(12.0,precision)) + 

            NINE*r*xi*Power(xj,cl_float(12.0,precision)) - cl_float(22.0,precision)*r*Power(xi,SEVEN)*Power(xj,SIX)*

             (-cl_float(135.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            SIX*Power(xi,TWO)*Power(xj,cl_float(10.0,precision))*(-cl_float(11.0,precision) + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(66.0,precision)*Power(xi,FOUR)*Power(xj,EIGHT)*(-FIVE + Power(r,TWO)*Power(xj,TWO)) + 

            EIGHT*r*Power(xi,FIVE)*Power(xj,EIGHT)*(cl_float(99.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            r*Power(xi,THREE)*Power(xj,cl_float(10.0,precision))*(-cl_float(99.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(132.0,precision)*Power(xi,SIX)*Power(xj,SIX)*(cl_float(27.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(78.0,precision)*Power(xi,cl_float(12.0,precision))*(SEVEN + THREE*Power(r,TWO)*Power(xj,TWO)) - 

            TWO*r*Power(xi,cl_float(13.0,precision))*(cl_float(117.0,precision) + FOUR*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(66.0,precision)*Power(xi,EIGHT)*Power(xj,FOUR)*(-cl_float(191.0,precision) + SIX*Power(r,TWO)*Power(xj,TWO)) + 

            r*Power(xi,cl_float(11.0,precision))*Power(xj,TWO)*(-cl_float(2151.0,precision) + cl_float(22.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            SIX*Power(xi,cl_float(10.0,precision))*Power(xj,TWO)*(-cl_float(1099.0,precision) + cl_float(33.0,precision)*Power(r,TWO)*Power(xj,TWO))) + 

         exp(TWO*r*xi)*Power(xi,SIX)*

          (-cl_float(385.0,precision)*Power(xi,EIGHT)*Power(xj,EIGHT)*

             (cl_float(1080.0,precision) + cl_float(1935.0,precision)*r*xj + cl_float(1350.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1170.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(420.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(66.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + FOUR*Power(r,SIX)*Power(xj,SIX)) + 

            SEVEN*Power(xi,SIX)*Power(xj,cl_float(10.0,precision))*

             (cl_float(99540.0,precision) + cl_float(58095.0,precision)*r*xj + cl_float(190710.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(100950.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(21660.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(1938.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + FOUR*Power(r,SIX)*Power(xj,SIX) - 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            FOUR*Power(xj,cl_float(16.0,precision))*(cl_float(135135.0,precision) + cl_float(135135.0,precision)*r*xj + 

               cl_float(62370.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(17325.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(3150.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(378.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SIX) + Power(r,SEVEN)*Power(xj,SEVEN)) - 

            Power(xi,cl_float(16.0,precision))*(cl_float(1260.0,precision) + cl_float(2205.0,precision)*r*xj + cl_float(1890.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1050.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(420.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(126.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            SEVEN*Power(xi,FOUR)*Power(xj,cl_float(12.0,precision))*

             (cl_float(114660.0,precision) - cl_float(343395.0,precision)*r*xj - cl_float(242910.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(61950.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(6060.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(282.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(116.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            SEVEN*Power(xi,cl_float(12.0,precision))*Power(xj,FOUR)*

             (cl_float(9900.0,precision) + cl_float(17325.0,precision)*r*xj + cl_float(14850.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(8250.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(3300.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(1074.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(164.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            SEVEN*Power(xi,cl_float(10.0,precision))*Power(xj,SIX)*

             (cl_float(29700.0,precision) + cl_float(51975.0,precision)*r*xj + cl_float(44550.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(23850.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(11700.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(2814.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(284.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            Power(xi,cl_float(14.0,precision))*Power(xj,TWO)*

             (cl_float(13860.0,precision) + cl_float(24255.0,precision)*r*xj + cl_float(20790.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(11550.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(4620.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(1386.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(308.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(24.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            Power(xi,TWO)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(3063060.0,precision) - cl_float(1936935.0,precision)*r*xj - cl_float(408870.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(11550.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(23100.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(5082.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(532.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(24.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN))))/

       (cl_float(1260.0,precision)*exp(TWO*r*(xi + xj))*Power(r,TWO)*Power(xi - xj,cl_float(11.0,precision))*

         Power(xi + xj,cl_float(11.0,precision))) + (cl_float(1260.0,precision)*exp(TWO*r*(xi + xj))*

          Power(Power(xi,TWO) - Power(xj,TWO),cl_float(11.0,precision)) + 

         cl_float(210.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(10.0,precision))*

          (-cl_float(36.0,precision)*Power(r,TWO)*Power(xi,cl_float(14.0,precision)) - TWO*Power(r,THREE)*Power(xi,cl_float(15.0,precision)) - 

            cl_float(1287.0,precision)*r*Power(xi,NINE)*Power(xj,FOUR) + SIX*Power(xj,cl_float(12.0,precision)) + 

            NINE*r*xi*Power(xj,cl_float(12.0,precision)) - cl_float(22.0,precision)*r*Power(xi,SEVEN)*Power(xj,SIX)*

             (-cl_float(135.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            SIX*Power(xi,TWO)*Power(xj,cl_float(10.0,precision))*(-cl_float(11.0,precision) + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(66.0,precision)*Power(xi,FOUR)*Power(xj,EIGHT)*(-FIVE + Power(r,TWO)*Power(xj,TWO)) + 

            EIGHT*r*Power(xi,FIVE)*Power(xj,EIGHT)*(cl_float(99.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            r*Power(xi,THREE)*Power(xj,cl_float(10.0,precision))*(-cl_float(99.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(132.0,precision)*Power(xi,SIX)*Power(xj,SIX)*(cl_float(27.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(78.0,precision)*Power(xi,cl_float(12.0,precision))*(SEVEN + THREE*Power(r,TWO)*Power(xj,TWO)) - 

            TWO*r*Power(xi,cl_float(13.0,precision))*(cl_float(117.0,precision) + FOUR*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(66.0,precision)*Power(xi,EIGHT)*Power(xj,FOUR)*(-cl_float(191.0,precision) + SIX*Power(r,TWO)*Power(xj,TWO)) + 

            r*Power(xi,cl_float(11.0,precision))*Power(xj,TWO)*(-cl_float(2151.0,precision) + cl_float(22.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            SIX*Power(xi,cl_float(10.0,precision))*Power(xj,TWO)*(-cl_float(1099.0,precision) + cl_float(33.0,precision)*Power(r,TWO)*Power(xj,TWO))) + 

         exp(TWO*r*xi)*Power(xi,SIX)*

          (-cl_float(385.0,precision)*Power(xi,EIGHT)*Power(xj,EIGHT)*

             (cl_float(1080.0,precision) + cl_float(1935.0,precision)*r*xj + cl_float(1350.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1170.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(420.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(66.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + FOUR*Power(r,SIX)*Power(xj,SIX)) + 

            SEVEN*Power(xi,SIX)*Power(xj,cl_float(10.0,precision))*

             (cl_float(99540.0,precision) + cl_float(58095.0,precision)*r*xj + cl_float(190710.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(100950.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(21660.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(1938.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + FOUR*Power(r,SIX)*Power(xj,SIX) - 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            FOUR*Power(xj,cl_float(16.0,precision))*(cl_float(135135.0,precision) + cl_float(135135.0,precision)*r*xj + 

               cl_float(62370.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(17325.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(3150.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(378.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SIX) + Power(r,SEVEN)*Power(xj,SEVEN)) - 

            Power(xi,cl_float(16.0,precision))*(cl_float(1260.0,precision) + cl_float(2205.0,precision)*r*xj + cl_float(1890.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1050.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(420.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(126.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            SEVEN*Power(xi,FOUR)*Power(xj,cl_float(12.0,precision))*

             (cl_float(114660.0,precision) - cl_float(343395.0,precision)*r*xj - cl_float(242910.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(61950.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(6060.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(282.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(116.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            SEVEN*Power(xi,cl_float(12.0,precision))*Power(xj,FOUR)*

             (cl_float(9900.0,precision) + cl_float(17325.0,precision)*r*xj + cl_float(14850.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(8250.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(3300.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(1074.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(164.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            SEVEN*Power(xi,cl_float(10.0,precision))*Power(xj,SIX)*

             (cl_float(29700.0,precision) + cl_float(51975.0,precision)*r*xj + cl_float(44550.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(23850.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(11700.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(2814.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(284.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            Power(xi,cl_float(14.0,precision))*Power(xj,TWO)*

             (cl_float(13860.0,precision) + cl_float(24255.0,precision)*r*xj + cl_float(20790.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(11550.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(4620.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(1386.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(308.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(24.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            Power(xi,TWO)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(3063060.0,precision) - cl_float(1936935.0,precision)*r*xj - cl_float(408870.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(11550.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(23100.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(5082.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(532.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(24.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN))))/

       (cl_float(630.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,cl_float(11.0,precision))*Power(xi + xj,cl_float(10.0,precision))) - 

      (cl_float(2520.0,precision)*exp(TWO*r*(xi + xj))*(xi + xj)*

          Power(Power(xi,TWO) - Power(xj,TWO),cl_float(11.0,precision)) + 

         cl_float(210.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(10.0,precision))*

          (-cl_float(72.0,precision)*r*Power(xi,cl_float(14.0,precision)) - SIX*Power(r,TWO)*Power(xi,cl_float(15.0,precision)) - 

            cl_float(468.0,precision)*r*Power(xi,cl_float(12.0,precision))*Power(xj,TWO) - 

            cl_float(16.0,precision)*Power(r,TWO)*Power(xi,cl_float(13.0,precision))*Power(xj,TWO) - 

            cl_float(1287.0,precision)*Power(xi,NINE)*Power(xj,FOUR) + cl_float(396.0,precision)*r*Power(xi,cl_float(10.0,precision))*Power(xj,FOUR) + 

            cl_float(44.0,precision)*Power(r,TWO)*Power(xi,cl_float(11.0,precision))*Power(xj,FOUR) + 

            cl_float(792.0,precision)*r*Power(xi,EIGHT)*Power(xj,SIX) - cl_float(528.0,precision)*r*Power(xi,SIX)*Power(xj,EIGHT) - 

            cl_float(44.0,precision)*Power(r,TWO)*Power(xi,SEVEN)*Power(xj,EIGHT) - 

            cl_float(132.0,precision)*r*Power(xi,FOUR)*Power(xj,cl_float(10.0,precision)) + 

            cl_float(16.0,precision)*Power(r,TWO)*Power(xi,FIVE)*Power(xj,cl_float(10.0,precision)) + NINE*xi*Power(xj,cl_float(12.0,precision)) + 

            cl_float(12.0,precision)*r*Power(xi,TWO)*Power(xj,cl_float(12.0,precision)) + 

            FOUR*Power(r,TWO)*Power(xi,THREE)*Power(xj,cl_float(12.0,precision)) - 

            cl_float(22.0,precision)*Power(xi,SEVEN)*Power(xj,SIX)*(-cl_float(135.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            EIGHT*Power(xi,FIVE)*Power(xj,EIGHT)*(cl_float(99.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            Power(xi,THREE)*Power(xj,cl_float(10.0,precision))*(-cl_float(99.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) - 

            TWO*Power(xi,cl_float(13.0,precision))*(cl_float(117.0,precision) + FOUR*Power(r,TWO)*Power(xj,TWO)) + 

            Power(xi,cl_float(11.0,precision))*Power(xj,TWO)*(-cl_float(2151.0,precision) + cl_float(22.0,precision)*Power(r,TWO)*Power(xj,TWO))) + 

         cl_float(420.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(11.0,precision))*

          (-cl_float(36.0,precision)*Power(r,TWO)*Power(xi,cl_float(14.0,precision)) - TWO*Power(r,THREE)*Power(xi,cl_float(15.0,precision)) - 

            cl_float(1287.0,precision)*r*Power(xi,NINE)*Power(xj,FOUR) + SIX*Power(xj,cl_float(12.0,precision)) + 

            NINE*r*xi*Power(xj,cl_float(12.0,precision)) - cl_float(22.0,precision)*r*Power(xi,SEVEN)*Power(xj,SIX)*

             (-cl_float(135.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            SIX*Power(xi,TWO)*Power(xj,cl_float(10.0,precision))*(-cl_float(11.0,precision) + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(66.0,precision)*Power(xi,FOUR)*Power(xj,EIGHT)*(-FIVE + Power(r,TWO)*Power(xj,TWO)) + 

            EIGHT*r*Power(xi,FIVE)*Power(xj,EIGHT)*(cl_float(99.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            r*Power(xi,THREE)*Power(xj,cl_float(10.0,precision))*(-cl_float(99.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(132.0,precision)*Power(xi,SIX)*Power(xj,SIX)*(cl_float(27.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(78.0,precision)*Power(xi,cl_float(12.0,precision))*(SEVEN + THREE*Power(r,TWO)*Power(xj,TWO)) - 

            TWO*r*Power(xi,cl_float(13.0,precision))*(cl_float(117.0,precision) + FOUR*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(66.0,precision)*Power(xi,EIGHT)*Power(xj,FOUR)*(-cl_float(191.0,precision) + SIX*Power(r,TWO)*Power(xj,TWO)) + 

            r*Power(xi,cl_float(11.0,precision))*Power(xj,TWO)*(-cl_float(2151.0,precision) + cl_float(22.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            SIX*Power(xi,cl_float(10.0,precision))*Power(xj,TWO)*(-cl_float(1099.0,precision) + cl_float(33.0,precision)*Power(r,TWO)*Power(xj,TWO))) + 

         exp(TWO*r*xi)*Power(xi,SIX)*

          (-cl_float(385.0,precision)*Power(xi,EIGHT)*Power(xj,EIGHT)*

             (cl_float(1935.0,precision)*xj + cl_float(2700.0,precision)*r*Power(xj,TWO) + cl_float(3510.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(1680.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(330.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(24.0,precision)*Power(r,FIVE)*Power(xj,SIX)) + 

            SEVEN*Power(xi,SIX)*Power(xj,cl_float(10.0,precision))*

             (cl_float(58095.0,precision)*xj + cl_float(381420.0,precision)*r*Power(xj,TWO) + 

               cl_float(302850.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(86640.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(9690.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(24.0,precision)*Power(r,FIVE)*Power(xj,SIX) - 

               cl_float(56.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) + 

            FOUR*Power(xj,cl_float(16.0,precision))*(cl_float(135135.0,precision)*xj + cl_float(124740.0,precision)*r*Power(xj,TWO) + 

               cl_float(51975.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(12600.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(1890.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(168.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               SEVEN*Power(r,SIX)*Power(xj,SEVEN)) - 

            Power(xi,cl_float(16.0,precision))*(cl_float(2205.0,precision)*xj + cl_float(3780.0,precision)*r*Power(xj,TWO) + 

               cl_float(3150.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(1680.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(630.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(168.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) + 

            SEVEN*Power(xi,FOUR)*Power(xj,cl_float(12.0,precision))*

             (-cl_float(343395.0,precision)*xj - cl_float(485820.0,precision)*r*Power(xj,TWO) - 

               cl_float(185850.0,precision)*Power(r,TWO)*Power(xj,THREE) - cl_float(24240.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(1410.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(696.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(56.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) - 

            SEVEN*Power(xi,cl_float(12.0,precision))*Power(xj,FOUR)*

             (cl_float(17325.0,precision)*xj + cl_float(29700.0,precision)*r*Power(xj,TWO) + cl_float(24750.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(13200.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(5370.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(984.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(56.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) + 

            SEVEN*Power(xi,cl_float(10.0,precision))*Power(xj,SIX)*

             (cl_float(51975.0,precision)*xj + cl_float(89100.0,precision)*r*Power(xj,TWO) + cl_float(71550.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(46800.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(14070.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(1704.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(56.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) + 

            Power(xi,cl_float(14.0,precision))*Power(xj,TWO)*

             (cl_float(24255.0,precision)*xj + cl_float(41580.0,precision)*r*Power(xj,TWO) + cl_float(34650.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(18480.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(6930.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(1848.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(168.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) - 

            Power(xi,TWO)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(1936935.0,precision)*xj - cl_float(817740.0,precision)*r*Power(xj,TWO) + 

               cl_float(34650.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(92400.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(25410.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(3192.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(168.0,precision)*Power(r,SIX)*Power(xj,SEVEN))) + 

         TWO*exp(TWO*r*xi)*Power(xi,SEVEN)*

          (-cl_float(385.0,precision)*Power(xi,EIGHT)*Power(xj,EIGHT)*

             (cl_float(1080.0,precision) + cl_float(1935.0,precision)*r*xj + cl_float(1350.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1170.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(420.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(66.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + FOUR*Power(r,SIX)*Power(xj,SIX)) + 

            SEVEN*Power(xi,SIX)*Power(xj,cl_float(10.0,precision))*

             (cl_float(99540.0,precision) + cl_float(58095.0,precision)*r*xj + cl_float(190710.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(100950.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(21660.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(1938.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + FOUR*Power(r,SIX)*Power(xj,SIX) - 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            FOUR*Power(xj,cl_float(16.0,precision))*(cl_float(135135.0,precision) + cl_float(135135.0,precision)*r*xj + 

               cl_float(62370.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(17325.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(3150.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(378.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SIX) + Power(r,SEVEN)*Power(xj,SEVEN)) - 

            Power(xi,cl_float(16.0,precision))*(cl_float(1260.0,precision) + cl_float(2205.0,precision)*r*xj + cl_float(1890.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1050.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(420.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(126.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            SEVEN*Power(xi,FOUR)*Power(xj,cl_float(12.0,precision))*

             (cl_float(114660.0,precision) - cl_float(343395.0,precision)*r*xj - cl_float(242910.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(61950.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(6060.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(282.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(116.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            SEVEN*Power(xi,cl_float(12.0,precision))*Power(xj,FOUR)*

             (cl_float(9900.0,precision) + cl_float(17325.0,precision)*r*xj + cl_float(14850.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(8250.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(3300.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(1074.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(164.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            SEVEN*Power(xi,cl_float(10.0,precision))*Power(xj,SIX)*

             (cl_float(29700.0,precision) + cl_float(51975.0,precision)*r*xj + cl_float(44550.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(23850.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(11700.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(2814.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(284.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            Power(xi,cl_float(14.0,precision))*Power(xj,TWO)*

             (cl_float(13860.0,precision) + cl_float(24255.0,precision)*r*xj + cl_float(20790.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(11550.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(4620.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(1386.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(308.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(24.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            Power(xi,TWO)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(3063060.0,precision) - cl_float(1936935.0,precision)*r*xj - cl_float(408870.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(11550.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(23100.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(5082.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(532.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(24.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN))))/

       (cl_float(1260.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,cl_float(11.0,precision))*Power(xi + xj,cl_float(11.0,precision)))

    ; }
   
  }
  return S;
}

cl_F DSlater_2S_5S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = ZERO

    ; } else {  S = -(-cl_float(224622748350.0,precision)*xi + cl_float(249080832000.0,precision)*exp(TWO*r*xi)*xi - 

          cl_float(400329329400.0,precision)*r*Power(xi,TWO) - cl_float(351747621225.0,precision)*Power(r,TWO)*Power(xi,THREE) - 

          cl_float(202556554200.0,precision)*Power(r,THREE)*Power(xi,FOUR) - 

          cl_float(85662076500.0,precision)*Power(r,FOUR)*Power(xi,FIVE) - 

          cl_float(28229160960.0,precision)*Power(r,FIVE)*Power(xi,SIX) - 

          cl_float(7498370880.0,precision)*Power(r,SIX)*Power(xi,SEVEN) - 

          cl_float(1635828480.0,precision)*Power(r,SEVEN)*Power(xi,EIGHT) - 

          cl_float(295289280.0,precision)*Power(r,EIGHT)*Power(xi,NINE) - cl_float(43929600.0,precision)*Power(r,NINE)*Power(xi,cl_float(10.0,precision)) - 

          cl_float(5271552.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(11.0,precision)) - cl_float(479232.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(12.0,precision)) - 

          cl_float(26624.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(13.0,precision)))/(cl_float(1.24540416e11,precision)*exp(TWO*r*xi)*r) + 

      (-cl_float(124540416000.0,precision) + cl_float(124540416000.0,precision)*exp(TWO*r*xi) - cl_float(224622748350.0,precision)*r*xi - 

         cl_float(200164664700.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

         cl_float(117249207075.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

         cl_float(50639138550.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - 

         cl_float(17132415300.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

         cl_float(4704860160.0,precision)*Power(r,SIX)*Power(xi,SIX) - 

         cl_float(1071195840.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - cl_float(204478560.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - 

         cl_float(32809920.0,precision)*Power(r,NINE)*Power(xi,NINE) - cl_float(4392960.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - 

         cl_float(479232.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)) - cl_float(39936.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(12.0,precision)) - 

         cl_float(2048.0,precision)*Power(r,cl_float(13.0,precision))*Power(xi,cl_float(13.0,precision)))/

       (cl_float(1.24540416e11,precision)*exp(TWO*r*xi)*Power(r,TWO)) + 

      (xi*(-cl_float(124540416000.0,precision) + cl_float(124540416000.0,precision)*exp(TWO*r*xi) - cl_float(224622748350.0,precision)*r*xi - 

           cl_float(200164664700.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

           cl_float(117249207075.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

           cl_float(50639138550.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - 

           cl_float(17132415300.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

           cl_float(4704860160.0,precision)*Power(r,SIX)*Power(xi,SIX) - 

           cl_float(1071195840.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

           cl_float(204478560.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - cl_float(32809920.0,precision)*Power(r,NINE)*Power(xi,NINE) - 

           cl_float(4392960.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - cl_float(479232.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)) - 

           cl_float(39936.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(12.0,precision)) - cl_float(2048.0,precision)*Power(r,cl_float(13.0,precision))*Power(xi,cl_float(13.0,precision))))/

       (cl_float(6.2270208e10,precision)*exp(TWO*r*xi)*r)

    ; }
 
  }
  else {
      if (r == ZERO) {  S = ZERO

    ; } else {  S = (cl_float(28350.0,precision)*exp(TWO*r*(xi + xj))*Power(Power(xi,TWO) - Power(xj,TWO),cl_float(13.0,precision)) + 

         cl_float(945.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(12.0,precision))*

          (-cl_float(210.0,precision)*Power(r,TWO)*Power(xi,cl_float(16.0,precision)) - cl_float(10.0,precision)*Power(r,THREE)*Power(xi,cl_float(17.0,precision)) + 

            cl_float(30.0,precision)*Power(xj,cl_float(14.0,precision)) + cl_float(45.0,precision)*r*xi*Power(xj,cl_float(14.0,precision)) + 

            cl_float(39.0,precision)*r*Power(xi,SEVEN)*Power(xj,EIGHT)*(cl_float(1309.0,precision) - TWO*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(858.0,precision)*Power(xi,EIGHT)*Power(xj,SIX)*(-cl_float(305.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(30.0,precision)*Power(xi,TWO)*Power(xj,cl_float(12.0,precision))*(-cl_float(13.0,precision) + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(390.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(10.0,precision))*(-SIX + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(143.0,precision)*r*Power(xi,NINE)*Power(xj,SIX)*(-cl_float(153.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            FIVE*r*Power(xi,THREE)*Power(xj,cl_float(12.0,precision))*(-cl_float(117.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(45.0,precision)*r*Power(xi,cl_float(15.0,precision))*(cl_float(35.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(138.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,TWO)*(cl_float(580.0,precision) + cl_float(13.0,precision)*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(150.0,precision)*Power(xi,cl_float(14.0,precision))*(cl_float(28.0,precision) + cl_float(17.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(13.0,precision)*r*Power(xi,cl_float(11.0,precision))*Power(xj,FOUR)*

             (-cl_float(4071.0,precision) + cl_float(22.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            THREE*r*Power(xi,cl_float(13.0,precision))*Power(xj,TWO)*(-cl_float(8135.0,precision) + cl_float(26.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            THREE*r*Power(xi,FIVE)*Power(xj,cl_float(10.0,precision))*(cl_float(2171.0,precision) + cl_float(30.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(234.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,FOUR)*(-cl_float(1235.0,precision) + cl_float(33.0,precision)*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(78.0,precision)*Power(xi,SIX)*Power(xj,EIGHT)*(cl_float(550.0,precision) + cl_float(47.0,precision)*Power(r,TWO)*Power(xj,TWO))) - 

         TWO*exp(TWO*r*xi)*Power(xi,SIX)*

          (-cl_float(819.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(10.0,precision))*

             (cl_float(22275.0,precision) + cl_float(39780.0,precision)*r*xj + cl_float(38160.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(16560.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(9840.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(3900.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(816.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(88.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + FOUR*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            Power(xi,cl_float(20.0,precision))*(cl_float(14175.0,precision) + cl_float(25515.0,precision)*r*xj + cl_float(22680.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(13230.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(5670.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(1890.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(504.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(108.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               TWO*Power(r,NINE)*Power(xj,NINE)) - 

            Power(xj,cl_float(20.0,precision))*(cl_float(16216200.0,precision) + cl_float(18243225.0,precision)*r*xj + 

               cl_float(9729720.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(3243240.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(748440.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(124740.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(15120.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(1296.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(72.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + TWO*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(18.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,FOUR)*

             (cl_float(61425.0,precision) + cl_float(110565.0,precision)*r*xj + cl_float(98280.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(57330.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(24570.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(8190.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(2184.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(496.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(64.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               THREE*Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(18.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(16.0,precision))*

             (cl_float(6572475.0,precision) - cl_float(3161340.0,precision)*r*xj - cl_float(4782960.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1912365.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(378105.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(34125.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(1092.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(650.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(71.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + THREE*Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(21.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(12.0,precision))*

             (-cl_float(1063800.0,precision) - cl_float(2775735.0,precision)*r*xj - cl_float(862920.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1132020.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(698580.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(196920.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(28992.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(2064.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - 

               cl_float(24.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + FOUR*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(21.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,EIGHT)*

             (cl_float(482625.0,precision) + cl_float(868725.0,precision)*r*xj + cl_float(772200.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(455400.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(178200.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(72180.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(19920.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(2952.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(204.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               FOUR*Power(r,NINE)*Power(xj,NINE)) + 

            SIX*Power(xi,SIX)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(10357200.0,precision) + cl_float(5071815.0,precision)*r*xj - cl_float(6463800.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(7151130.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(2572290.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(468720.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(42672.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

               cl_float(648.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(228.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               cl_float(16.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

            Power(xi,cl_float(18.0,precision))*Power(xj,TWO)*

             (cl_float(184275.0,precision) + cl_float(331695.0,precision)*r*xj + cl_float(294840.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(171990.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(73710.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(24570.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(6552.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(1404.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(234.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               cl_float(16.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

            Power(xi,TWO)*Power(xj,cl_float(18.0,precision))*

             (-cl_float(133783650.0,precision) - cl_float(107432325.0,precision)*r*xj - cl_float(35675640.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(5135130.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(270270.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(270270.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(57960.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(6948.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(486.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(16.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

            SIX*Power(xi,cl_float(14.0,precision))*Power(xj,SIX)*

             (cl_float(675675.0,precision) + cl_float(1216215.0,precision)*r*xj + cl_float(1081080.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(630630.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(270270.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(88200.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(26544.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(5160.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(492.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               cl_float(16.0,precision)*Power(r,NINE)*Power(xj,NINE))))/

       (cl_float(28350.0,precision)*exp(TWO*r*(xi + xj))*Power(r,TWO)*Power(xi - xj,cl_float(13.0,precision))*

         Power(xi + xj,cl_float(13.0,precision))) + (cl_float(28350.0,precision)*exp(TWO*r*(xi + xj))*

          Power(Power(xi,TWO) - Power(xj,TWO),cl_float(13.0,precision)) + 

         cl_float(945.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(12.0,precision))*

          (-cl_float(210.0,precision)*Power(r,TWO)*Power(xi,cl_float(16.0,precision)) - cl_float(10.0,precision)*Power(r,THREE)*Power(xi,cl_float(17.0,precision)) + 

            cl_float(30.0,precision)*Power(xj,cl_float(14.0,precision)) + cl_float(45.0,precision)*r*xi*Power(xj,cl_float(14.0,precision)) + 

            cl_float(39.0,precision)*r*Power(xi,SEVEN)*Power(xj,EIGHT)*(cl_float(1309.0,precision) - TWO*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(858.0,precision)*Power(xi,EIGHT)*Power(xj,SIX)*(-cl_float(305.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(30.0,precision)*Power(xi,TWO)*Power(xj,cl_float(12.0,precision))*(-cl_float(13.0,precision) + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(390.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(10.0,precision))*(-SIX + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(143.0,precision)*r*Power(xi,NINE)*Power(xj,SIX)*(-cl_float(153.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            FIVE*r*Power(xi,THREE)*Power(xj,cl_float(12.0,precision))*(-cl_float(117.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(45.0,precision)*r*Power(xi,cl_float(15.0,precision))*(cl_float(35.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(138.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,TWO)*(cl_float(580.0,precision) + cl_float(13.0,precision)*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(150.0,precision)*Power(xi,cl_float(14.0,precision))*(cl_float(28.0,precision) + cl_float(17.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(13.0,precision)*r*Power(xi,cl_float(11.0,precision))*Power(xj,FOUR)*

             (-cl_float(4071.0,precision) + cl_float(22.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            THREE*r*Power(xi,cl_float(13.0,precision))*Power(xj,TWO)*(-cl_float(8135.0,precision) + cl_float(26.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            THREE*r*Power(xi,FIVE)*Power(xj,cl_float(10.0,precision))*(cl_float(2171.0,precision) + cl_float(30.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(234.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,FOUR)*(-cl_float(1235.0,precision) + cl_float(33.0,precision)*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(78.0,precision)*Power(xi,SIX)*Power(xj,EIGHT)*(cl_float(550.0,precision) + cl_float(47.0,precision)*Power(r,TWO)*Power(xj,TWO))) - 

         TWO*exp(TWO*r*xi)*Power(xi,SIX)*

          (-cl_float(819.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(10.0,precision))*

             (cl_float(22275.0,precision) + cl_float(39780.0,precision)*r*xj + cl_float(38160.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(16560.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(9840.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(3900.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(816.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(88.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + FOUR*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            Power(xi,cl_float(20.0,precision))*(cl_float(14175.0,precision) + cl_float(25515.0,precision)*r*xj + cl_float(22680.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(13230.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(5670.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(1890.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(504.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(108.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               TWO*Power(r,NINE)*Power(xj,NINE)) - 

            Power(xj,cl_float(20.0,precision))*(cl_float(16216200.0,precision) + cl_float(18243225.0,precision)*r*xj + 

               cl_float(9729720.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(3243240.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(748440.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(124740.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(15120.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(1296.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(72.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + TWO*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(18.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,FOUR)*

             (cl_float(61425.0,precision) + cl_float(110565.0,precision)*r*xj + cl_float(98280.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(57330.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(24570.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(8190.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(2184.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(496.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(64.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               THREE*Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(18.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(16.0,precision))*

             (cl_float(6572475.0,precision) - cl_float(3161340.0,precision)*r*xj - cl_float(4782960.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1912365.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(378105.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(34125.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(1092.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(650.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(71.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + THREE*Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(21.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(12.0,precision))*

             (-cl_float(1063800.0,precision) - cl_float(2775735.0,precision)*r*xj - cl_float(862920.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1132020.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(698580.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(196920.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(28992.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(2064.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - 

               cl_float(24.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + FOUR*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(21.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,EIGHT)*

             (cl_float(482625.0,precision) + cl_float(868725.0,precision)*r*xj + cl_float(772200.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(455400.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(178200.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(72180.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(19920.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(2952.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(204.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               FOUR*Power(r,NINE)*Power(xj,NINE)) + 

            SIX*Power(xi,SIX)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(10357200.0,precision) + cl_float(5071815.0,precision)*r*xj - cl_float(6463800.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(7151130.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(2572290.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(468720.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(42672.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

               cl_float(648.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(228.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               cl_float(16.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

            Power(xi,cl_float(18.0,precision))*Power(xj,TWO)*

             (cl_float(184275.0,precision) + cl_float(331695.0,precision)*r*xj + cl_float(294840.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(171990.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(73710.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(24570.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(6552.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(1404.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(234.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               cl_float(16.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

            Power(xi,TWO)*Power(xj,cl_float(18.0,precision))*

             (-cl_float(133783650.0,precision) - cl_float(107432325.0,precision)*r*xj - cl_float(35675640.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(5135130.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(270270.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(270270.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(57960.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(6948.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(486.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(16.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

            SIX*Power(xi,cl_float(14.0,precision))*Power(xj,SIX)*

             (cl_float(675675.0,precision) + cl_float(1216215.0,precision)*r*xj + cl_float(1081080.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(630630.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(270270.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(88200.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(26544.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(5160.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(492.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               cl_float(16.0,precision)*Power(r,NINE)*Power(xj,NINE))))/

       (cl_float(14175.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,cl_float(13.0,precision))*Power(xi + xj,cl_float(12.0,precision))) - 

      (cl_float(56700.0,precision)*exp(TWO*r*(xi + xj))*(xi + xj)*

          Power(Power(xi,TWO) - Power(xj,TWO),cl_float(13.0,precision)) + 

         cl_float(945.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(12.0,precision))*

          (-cl_float(420.0,precision)*r*Power(xi,cl_float(16.0,precision)) - cl_float(30.0,precision)*Power(r,TWO)*Power(xi,cl_float(17.0,precision)) - 

            cl_float(5100.0,precision)*r*Power(xi,cl_float(14.0,precision))*Power(xj,TWO) - 

            cl_float(180.0,precision)*Power(r,TWO)*Power(xi,cl_float(15.0,precision))*Power(xj,TWO) - 

            cl_float(3588.0,precision)*r*Power(xi,cl_float(12.0,precision))*Power(xj,FOUR) + 

            cl_float(156.0,precision)*Power(r,TWO)*Power(xi,cl_float(13.0,precision))*Power(xj,FOUR) + 

            cl_float(15444.0,precision)*r*Power(xi,cl_float(10.0,precision))*Power(xj,SIX) + 

            cl_float(572.0,precision)*Power(r,TWO)*Power(xi,cl_float(11.0,precision))*Power(xj,SIX) + 

            cl_float(1716.0,precision)*r*Power(xi,EIGHT)*Power(xj,EIGHT) - 

            cl_float(572.0,precision)*Power(r,TWO)*Power(xi,NINE)*Power(xj,EIGHT) - 

            cl_float(7332.0,precision)*r*Power(xi,SIX)*Power(xj,cl_float(10.0,precision)) - 

            cl_float(156.0,precision)*Power(r,TWO)*Power(xi,SEVEN)*Power(xj,cl_float(10.0,precision)) - 

            cl_float(780.0,precision)*r*Power(xi,FOUR)*Power(xj,cl_float(12.0,precision)) + 

            cl_float(180.0,precision)*Power(r,TWO)*Power(xi,FIVE)*Power(xj,cl_float(12.0,precision)) + cl_float(45.0,precision)*xi*Power(xj,cl_float(14.0,precision)) + 

            cl_float(60.0,precision)*r*Power(xi,TWO)*Power(xj,cl_float(14.0,precision)) + 

            cl_float(20.0,precision)*Power(r,TWO)*Power(xi,THREE)*Power(xj,cl_float(14.0,precision)) + 

            cl_float(39.0,precision)*Power(xi,SEVEN)*Power(xj,EIGHT)*(cl_float(1309.0,precision) - TWO*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(143.0,precision)*Power(xi,NINE)*Power(xj,SIX)*(-cl_float(153.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            FIVE*Power(xi,THREE)*Power(xj,cl_float(12.0,precision))*(-cl_float(117.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(45.0,precision)*Power(xi,cl_float(15.0,precision))*(cl_float(35.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(13.0,precision)*Power(xi,cl_float(11.0,precision))*Power(xj,FOUR)*(-cl_float(4071.0,precision) + cl_float(22.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            THREE*Power(xi,cl_float(13.0,precision))*Power(xj,TWO)*(-cl_float(8135.0,precision) + cl_float(26.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            THREE*Power(xi,FIVE)*Power(xj,cl_float(10.0,precision))*(cl_float(2171.0,precision) + cl_float(30.0,precision)*Power(r,TWO)*Power(xj,TWO))) + 

         cl_float(1890.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(13.0,precision))*

          (-cl_float(210.0,precision)*Power(r,TWO)*Power(xi,cl_float(16.0,precision)) - cl_float(10.0,precision)*Power(r,THREE)*Power(xi,cl_float(17.0,precision)) + 

            cl_float(30.0,precision)*Power(xj,cl_float(14.0,precision)) + cl_float(45.0,precision)*r*xi*Power(xj,cl_float(14.0,precision)) + 

            cl_float(39.0,precision)*r*Power(xi,SEVEN)*Power(xj,EIGHT)*(cl_float(1309.0,precision) - TWO*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(858.0,precision)*Power(xi,EIGHT)*Power(xj,SIX)*(-cl_float(305.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(30.0,precision)*Power(xi,TWO)*Power(xj,cl_float(12.0,precision))*(-cl_float(13.0,precision) + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(390.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(10.0,precision))*(-SIX + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(143.0,precision)*r*Power(xi,NINE)*Power(xj,SIX)*(-cl_float(153.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            FIVE*r*Power(xi,THREE)*Power(xj,cl_float(12.0,precision))*(-cl_float(117.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(45.0,precision)*r*Power(xi,cl_float(15.0,precision))*(cl_float(35.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(138.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,TWO)*(cl_float(580.0,precision) + cl_float(13.0,precision)*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(150.0,precision)*Power(xi,cl_float(14.0,precision))*(cl_float(28.0,precision) + cl_float(17.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(13.0,precision)*r*Power(xi,cl_float(11.0,precision))*Power(xj,FOUR)*(-cl_float(4071.0,precision) + cl_float(22.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            THREE*r*Power(xi,cl_float(13.0,precision))*Power(xj,TWO)*(-cl_float(8135.0,precision) + cl_float(26.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            THREE*r*Power(xi,FIVE)*Power(xj,cl_float(10.0,precision))*(cl_float(2171.0,precision) + cl_float(30.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(234.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,FOUR)*(-cl_float(1235.0,precision) + cl_float(33.0,precision)*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(78.0,precision)*Power(xi,SIX)*Power(xj,EIGHT)*(cl_float(550.0,precision) + cl_float(47.0,precision)*Power(r,TWO)*Power(xj,TWO))) - 

         TWO*exp(TWO*r*xi)*Power(xi,SIX)*

          (-cl_float(819.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(10.0,precision))*

             (cl_float(39780.0,precision)*xj + cl_float(76320.0,precision)*r*Power(xj,TWO) + cl_float(49680.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(39360.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(19500.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(4896.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(616.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + 

               cl_float(32.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT)) + 

            Power(xi,cl_float(20.0,precision))*(cl_float(25515.0,precision)*xj + cl_float(45360.0,precision)*r*Power(xj,TWO) + 

               cl_float(39690.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(22680.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(9450.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(3024.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(756.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(144.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) - 

            Power(xj,cl_float(20.0,precision))*(cl_float(18243225.0,precision)*xj + cl_float(19459440.0,precision)*r*Power(xj,TWO) + 

               cl_float(9729720.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(2993760.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(623700.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(90720.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(9072.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(576.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            cl_float(18.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,FOUR)*

             (cl_float(110565.0,precision)*xj + cl_float(196560.0,precision)*r*Power(xj,TWO) + 

               cl_float(171990.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(98280.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(40950.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(13104.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(3472.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(512.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(27.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) - 

            cl_float(18.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(16.0,precision))*

             (-cl_float(3161340.0,precision)*xj - cl_float(9565920.0,precision)*r*Power(xj,TWO) - 

               cl_float(5737095.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(1512420.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(170625.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(6552.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(4550.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(568.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(27.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) - 

            cl_float(21.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(12.0,precision))*

             (-cl_float(2775735.0,precision)*xj - cl_float(1725840.0,precision)*r*Power(xj,TWO) - 

               cl_float(3396060.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(2794320.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(984600.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - cl_float(173952.0,precision)*Power(r,FIVE)*Power(xj,SIX) - 

               cl_float(14448.0,precision)*Power(r,SIX)*Power(xj,SEVEN) - cl_float(192.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(36.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            cl_float(21.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,EIGHT)*

             (cl_float(868725.0,precision)*xj + cl_float(1544400.0,precision)*r*Power(xj,TWO) + 

               cl_float(1366200.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(712800.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(360900.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(119520.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(20664.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + 

               cl_float(1632.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + cl_float(36.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            SIX*Power(xi,SIX)*Power(xj,cl_float(14.0,precision))*

             (cl_float(5071815.0,precision)*xj - cl_float(12927600.0,precision)*r*Power(xj,TWO) - 

               cl_float(21453390.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(10289160.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(2343600.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - 

               cl_float(256032.0,precision)*Power(r,FIVE)*Power(xj,SIX) - cl_float(4536.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + 

               cl_float(1824.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + cl_float(144.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) - 

            Power(xi,cl_float(18.0,precision))*Power(xj,TWO)*

             (cl_float(331695.0,precision)*xj + cl_float(589680.0,precision)*r*Power(xj,TWO) + 

               cl_float(515970.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(294840.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(122850.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(39312.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(9828.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(1872.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(144.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            Power(xi,TWO)*Power(xj,cl_float(18.0,precision))*

             (-cl_float(107432325.0,precision)*xj - cl_float(71351280.0,precision)*r*Power(xj,TWO) - 

               cl_float(15405390.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(1081080.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(1351350.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(347760.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(48636.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + 

               cl_float(3888.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + cl_float(144.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) - 

            SIX*Power(xi,cl_float(14.0,precision))*Power(xj,SIX)*

             (cl_float(1216215.0,precision)*xj + cl_float(2162160.0,precision)*r*Power(xj,TWO) + 

               cl_float(1891890.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(1081080.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(441000.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(159264.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(36120.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + 

               cl_float(3936.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + cl_float(144.0,precision)*Power(r,EIGHT)*Power(xj,NINE))) - 

         FOUR*exp(TWO*r*xi)*Power(xi,SEVEN)*

          (-cl_float(819.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(10.0,precision))*

             (cl_float(22275.0,precision) + cl_float(39780.0,precision)*r*xj + cl_float(38160.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(16560.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(9840.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(3900.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(816.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(88.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + FOUR*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            Power(xi,cl_float(20.0,precision))*(cl_float(14175.0,precision) + cl_float(25515.0,precision)*r*xj + cl_float(22680.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(13230.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(5670.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(1890.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(504.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(108.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               TWO*Power(r,NINE)*Power(xj,NINE)) - 

            Power(xj,cl_float(20.0,precision))*(cl_float(16216200.0,precision) + cl_float(18243225.0,precision)*r*xj + 

               cl_float(9729720.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(3243240.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(748440.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(124740.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(15120.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(1296.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(72.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               TWO*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(18.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,FOUR)*

             (cl_float(61425.0,precision) + cl_float(110565.0,precision)*r*xj + cl_float(98280.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(57330.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(24570.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(8190.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(2184.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(496.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(64.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               THREE*Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(18.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(16.0,precision))*

             (cl_float(6572475.0,precision) - cl_float(3161340.0,precision)*r*xj - cl_float(4782960.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1912365.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(378105.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(34125.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(1092.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(650.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(71.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               THREE*Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(21.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(12.0,precision))*

             (-cl_float(1063800.0,precision) - cl_float(2775735.0,precision)*r*xj - cl_float(862920.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1132020.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(698580.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(196920.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(28992.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

               cl_float(2064.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - cl_float(24.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               FOUR*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(21.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,EIGHT)*

             (cl_float(482625.0,precision) + cl_float(868725.0,precision)*r*xj + cl_float(772200.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(455400.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(178200.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(72180.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(19920.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(2952.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(204.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               FOUR*Power(r,NINE)*Power(xj,NINE)) + 

            SIX*Power(xi,SIX)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(10357200.0,precision) + cl_float(5071815.0,precision)*r*xj - cl_float(6463800.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(7151130.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(2572290.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(468720.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(42672.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(648.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(228.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(16.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

            Power(xi,cl_float(18.0,precision))*Power(xj,TWO)*

             (cl_float(184275.0,precision) + cl_float(331695.0,precision)*r*xj + cl_float(294840.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(171990.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(73710.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(24570.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(6552.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(1404.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(234.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               cl_float(16.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

            Power(xi,TWO)*Power(xj,cl_float(18.0,precision))*

             (-cl_float(133783650.0,precision) - cl_float(107432325.0,precision)*r*xj - cl_float(35675640.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(5135130.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(270270.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(270270.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(57960.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(6948.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(486.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               cl_float(16.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

            SIX*Power(xi,cl_float(14.0,precision))*Power(xj,SIX)*

             (cl_float(675675.0,precision) + cl_float(1216215.0,precision)*r*xj + cl_float(1081080.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(630630.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(270270.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(88200.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(26544.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(5160.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(492.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               cl_float(16.0,precision)*Power(r,NINE)*Power(xj,NINE))))/

       (cl_float(28350.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,cl_float(13.0,precision))*Power(xi + xj,cl_float(13.0,precision)))

    ; }
   
  }
  return S;
}

cl_F DSlater_2S_1S(cl_F r,cl_F xi,cl_F xj)
{
  return DSlater_1S_2S(r,xj,xi);
}

cl_F DSlater_3S_3S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = ZERO

    ; } else {  S = -(-cl_float(2503064025.0,precision)*xi + cl_float(2874009600.0,precision)*exp(TWO*r*xi)*xi - 

          cl_float(4264236900.0,precision)*r*Power(xi,TWO) - cl_float(3541992300.0,precision)*Power(r,TWO)*Power(xi,THREE) - 

          cl_float(1906027200.0,precision)*Power(r,THREE)*Power(xi,FOUR) - 

          cl_float(744282000.0,precision)*Power(r,FOUR)*Power(xi,FIVE) - cl_float(223534080.0,precision)*Power(r,FIVE)*Power(xi,SIX) - 

          cl_float(53222400.0,precision)*Power(r,SIX)*Power(xi,SEVEN) - cl_float(10137600.0,precision)*Power(r,SEVEN)*Power(xi,EIGHT) - 

          cl_float(1520640.0,precision)*Power(r,EIGHT)*Power(xi,NINE) - cl_float(168960.0,precision)*Power(r,NINE)*Power(xi,cl_float(10.0,precision)) - 

          cl_float(11264.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(11.0,precision)))/(cl_float(1.4370048e9,precision)*exp(TWO*r*xi)*r) + 

      (-cl_float(1437004800.0,precision) + cl_float(1437004800.0,precision)*exp(TWO*r*xi) - cl_float(2503064025.0,precision)*r*xi - 

         cl_float(2132118450.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

         cl_float(1180664100.0,precision)*Power(r,THREE)*Power(xi,THREE) - cl_float(476506800.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - 

         cl_float(148856400.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - cl_float(37255680.0,precision)*Power(r,SIX)*Power(xi,SIX) - 

         cl_float(7603200.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - cl_float(1267200.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - 

         cl_float(168960.0,precision)*Power(r,NINE)*Power(xi,NINE) - cl_float(16896.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - 

         cl_float(1024.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)))/(cl_float(1.4370048e9,precision)*exp(TWO*r*xi)*Power(r,TWO)) 

    + (xi*(-cl_float(1437004800.0,precision) + cl_float(1437004800.0,precision)*exp(TWO*r*xi) - cl_float(2503064025.0,precision)*r*xi - 

           cl_float(2132118450.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

           cl_float(1180664100.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

           cl_float(476506800.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - cl_float(148856400.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

           cl_float(37255680.0,precision)*Power(r,SIX)*Power(xi,SIX) - cl_float(7603200.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

           cl_float(1267200.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - cl_float(168960.0,precision)*Power(r,NINE)*Power(xi,NINE) - 

           cl_float(16896.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - cl_float(1024.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision))))/

       (cl_float(7.185024e8,precision)*exp(TWO*r*xi)*r)

    ; }
 
  }
  else {
      if (r == ZERO) {  S = ZERO

    ; } else {  S = (cl_float(135.0,precision)*exp(TWO*r*(xi + xj))*Power(Power(xi,TWO) - Power(xj,TWO),cl_float(11.0,precision)) + 

         exp(TWO*r*xj)*Power(xj,EIGHT)*

          (-cl_float(150.0,precision)*Power(r,FOUR)*Power(xi,cl_float(18.0,precision)) - SIX*Power(r,FIVE)*Power(xi,cl_float(19.0,precision)) + 

            cl_float(135.0,precision)*Power(xj,cl_float(14.0,precision)) + cl_float(225.0,precision)*r*xi*Power(xj,cl_float(14.0,precision)) + 

            cl_float(10.0,precision)*Power(r,THREE)*Power(xi,cl_float(17.0,precision))*(-cl_float(165.0,precision) + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(30.0,precision)*Power(r,TWO)*Power(xi,cl_float(16.0,precision))*(cl_float(330.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(45.0,precision)*r*Power(xi,THREE)*Power(xj,cl_float(12.0,precision))*(-cl_float(55.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(45.0,precision)*Power(xi,TWO)*Power(xj,cl_float(12.0,precision))*(-cl_float(33.0,precision) + FOUR*Power(r,TWO)*Power(xj,TWO)) + 

            r*Power(xi,NINE)*Power(xj,SIX)*

             (cl_float(234135.0,precision) - cl_float(4950.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(34.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

            FIVE*r*Power(xi,SEVEN)*Power(xj,EIGHT)*

             (cl_float(6237.0,precision) - cl_float(1242.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            THREE*r*Power(xi,FIVE)*Power(xj,cl_float(10.0,precision))*

             (cl_float(4125.0,precision) - cl_float(330.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(15.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(10.0,precision))*

             (cl_float(495.0,precision) - cl_float(132.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) - 

            cl_float(165.0,precision)*Power(xi,SIX)*Power(xj,EIGHT)*

             (cl_float(135.0,precision) - cl_float(60.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) - 

            FIVE*r*Power(xi,cl_float(13.0,precision))*Power(xj,TWO)*

             (cl_float(43875.0,precision) - cl_float(3438.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(22.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) 

    + FIVE*r*Power(xi,cl_float(11.0,precision))*Power(xj,FOUR)*

             (cl_float(7695.0,precision) - cl_float(2442.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(22.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) 

    + cl_float(15.0,precision)*Power(xi,EIGHT)*Power(xj,SIX)*(-cl_float(33.0,precision) - cl_float(3564.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(26.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            r*Power(xi,cl_float(15.0,precision))*(-cl_float(32175.0,precision) - cl_float(3690.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(34.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(15.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,FOUR)*

             (-cl_float(32277.0,precision) + cl_float(1364.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(66.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(15.0,precision)*Power(xi,cl_float(14.0,precision))*(-cl_float(3003.0,precision) - cl_float(2932.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(94.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

            cl_float(15.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,TWO)*

             (cl_float(28119.0,precision) - cl_float(5252.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(154.0,precision)*Power(r,FOUR)*Power(xj,FOUR))

    ) + exp(TWO*r*xi)*Power(xi,EIGHT)*

          (-FIVE*Power(xi,TWO)*Power(xj,cl_float(12.0,precision))*

             (-cl_float(84357.0,precision) - cl_float(43875.0,precision)*r*xj - cl_float(8796.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(738.0,precision)*Power(r,THREE)*Power(xj,THREE) - SIX*Power(r,FOUR)*Power(xj,FOUR) + 

               TWO*Power(r,FIVE)*Power(xj,FIVE)) - 

            THREE*Power(xi,cl_float(14.0,precision))*(cl_float(45.0,precision) + cl_float(75.0,precision)*r*xj + cl_float(60.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(30.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(10.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               TWO*Power(r,FIVE)*Power(xj,FIVE)) - 

            cl_float(55.0,precision)*Power(xi,EIGHT)*Power(xj,SIX)*

             (-cl_float(405.0,precision) - cl_float(567.0,precision)*r*xj - cl_float(972.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(90.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(18.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               TWO*Power(r,FIVE)*Power(xj,FIVE)) + 

            cl_float(55.0,precision)*Power(xi,SIX)*Power(xj,EIGHT)*

             (NINE - cl_float(4257.0,precision)*r*xj - cl_float(372.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(222.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(42.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               TWO*Power(r,FIVE)*Power(xj,FIVE)) + 

            THREE*Power(xj,cl_float(14.0,precision))*(cl_float(15015.0,precision) + cl_float(10725.0,precision)*r*xj + cl_float(3300.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(550.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(50.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               TWO*Power(r,FIVE)*Power(xj,FIVE)) + 

            FIVE*Power(xi,cl_float(12.0,precision))*Power(xj,TWO)*

             (cl_float(297.0,precision) + cl_float(495.0,precision)*r*xj + cl_float(396.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(198.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(66.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               TWO*Power(r,FIVE)*Power(xj,FIVE)) + 

            Power(xi,cl_float(10.0,precision))*Power(xj,FOUR)*

             (-cl_float(7425.0,precision) - cl_float(12375.0,precision)*r*xj - cl_float(9900.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(6210.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(390.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(34.0,precision)*Power(r,FIVE)*Power(xj,FIVE)) - 

            Power(xi,FOUR)*Power(xj,cl_float(10.0,precision))*

             (-cl_float(484155.0,precision) + cl_float(38475.0,precision)*r*xj + cl_float(78780.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(17190.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(1410.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(34.0,precision)*Power(r,FIVE)*Power(xj,FIVE))))/

       (cl_float(135.0,precision)*exp(TWO*r*(xi + xj))*Power(r,TWO)*Power(xi - xj,cl_float(11.0,precision))*

         Power(xi + xj,cl_float(11.0,precision))) + (TWO*(cl_float(135.0,precision)*exp(TWO*r*(xi + xj))*

            Power(Power(xi,TWO) - Power(xj,TWO),cl_float(11.0,precision)) + 

           exp(TWO*r*xj)*Power(xj,EIGHT)*

            (-cl_float(150.0,precision)*Power(r,FOUR)*Power(xi,cl_float(18.0,precision)) - SIX*Power(r,FIVE)*Power(xi,cl_float(19.0,precision)) + 

              cl_float(135.0,precision)*Power(xj,cl_float(14.0,precision)) + cl_float(225.0,precision)*r*xi*Power(xj,cl_float(14.0,precision)) + 

              cl_float(10.0,precision)*Power(r,THREE)*Power(xi,cl_float(17.0,precision))*(-cl_float(165.0,precision) + Power(r,TWO)*Power(xj,TWO)) - 

              cl_float(30.0,precision)*Power(r,TWO)*Power(xi,cl_float(16.0,precision))*(cl_float(330.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

              cl_float(45.0,precision)*r*Power(xi,THREE)*Power(xj,cl_float(12.0,precision))*(-cl_float(55.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

              cl_float(45.0,precision)*Power(xi,TWO)*Power(xj,cl_float(12.0,precision))*(-cl_float(33.0,precision) + FOUR*Power(r,TWO)*Power(xj,TWO)) + 

              r*Power(xi,NINE)*Power(xj,SIX)*

               (cl_float(234135.0,precision) - cl_float(4950.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(34.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

              FIVE*r*Power(xi,SEVEN)*Power(xj,EIGHT)*

               (cl_float(6237.0,precision) - cl_float(1242.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) 

    + THREE*r*Power(xi,FIVE)*Power(xj,cl_float(10.0,precision))*

               (cl_float(4125.0,precision) - cl_float(330.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) 

    + cl_float(15.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(10.0,precision))*

               (cl_float(495.0,precision) - cl_float(132.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) - 

              cl_float(165.0,precision)*Power(xi,SIX)*Power(xj,EIGHT)*

               (cl_float(135.0,precision) - cl_float(60.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) - 

              FIVE*r*Power(xi,cl_float(13.0,precision))*Power(xj,TWO)*

               (cl_float(43875.0,precision) - cl_float(3438.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(22.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

              FIVE*r*Power(xi,cl_float(11.0,precision))*Power(xj,FOUR)*

               (cl_float(7695.0,precision) - cl_float(2442.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(22.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

              cl_float(15.0,precision)*Power(xi,EIGHT)*Power(xj,SIX)*

               (-cl_float(33.0,precision) - cl_float(3564.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(26.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) 

    + r*Power(xi,cl_float(15.0,precision))*(-cl_float(32175.0,precision) - cl_float(3690.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(34.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

              cl_float(15.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,FOUR)*

               (-cl_float(32277.0,precision) + cl_float(1364.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(66.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

              cl_float(15.0,precision)*Power(xi,cl_float(14.0,precision))*(-cl_float(3003.0,precision) - cl_float(2932.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(94.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

              cl_float(15.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,TWO)*

               (cl_float(28119.0,precision) - cl_float(5252.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(154.0,precision)*Power(r,FOUR)*Power(xj,FOUR))) + 

           exp(TWO*r*xi)*Power(xi,EIGHT)*

            (-FIVE*Power(xi,TWO)*Power(xj,cl_float(12.0,precision))*

               (-cl_float(84357.0,precision) - cl_float(43875.0,precision)*r*xj - cl_float(8796.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(738.0,precision)*Power(r,THREE)*Power(xj,THREE) - SIX*Power(r,FOUR)*Power(xj,FOUR) + 

                 TWO*Power(r,FIVE)*Power(xj,FIVE)) - 

              THREE*Power(xi,cl_float(14.0,precision))*(cl_float(45.0,precision) + cl_float(75.0,precision)*r*xj + cl_float(60.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(30.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(10.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 TWO*Power(r,FIVE)*Power(xj,FIVE)) - 

              cl_float(55.0,precision)*Power(xi,EIGHT)*Power(xj,SIX)*

               (-cl_float(405.0,precision) - cl_float(567.0,precision)*r*xj - cl_float(972.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(90.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(18.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 TWO*Power(r,FIVE)*Power(xj,FIVE)) + 

              cl_float(55.0,precision)*Power(xi,SIX)*Power(xj,EIGHT)*

               (NINE - cl_float(4257.0,precision)*r*xj - cl_float(372.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(222.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(42.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 TWO*Power(r,FIVE)*Power(xj,FIVE)) + 

              THREE*Power(xj,cl_float(14.0,precision))*(cl_float(15015.0,precision) + cl_float(10725.0,precision)*r*xj + 

                 cl_float(3300.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(550.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

                 cl_float(50.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + TWO*Power(r,FIVE)*Power(xj,FIVE)) + 

              FIVE*Power(xi,cl_float(12.0,precision))*Power(xj,TWO)*

               (cl_float(297.0,precision) + cl_float(495.0,precision)*r*xj + cl_float(396.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(198.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(66.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 TWO*Power(r,FIVE)*Power(xj,FIVE)) + 

              Power(xi,cl_float(10.0,precision))*Power(xj,FOUR)*

               (-cl_float(7425.0,precision) - cl_float(12375.0,precision)*r*xj - cl_float(9900.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(6210.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(390.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(34.0,precision)*Power(r,FIVE)*Power(xj,FIVE)) - 

              Power(xi,FOUR)*Power(xj,cl_float(10.0,precision))*

               (-cl_float(484155.0,precision) + cl_float(38475.0,precision)*r*xj + cl_float(78780.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(17190.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(1410.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(34.0,precision)*Power(r,FIVE)*Power(xj,FIVE)))))/

       (cl_float(135.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,cl_float(11.0,precision))*Power(xi + xj,cl_float(10.0,precision))) - 

      (cl_float(270.0,precision)*exp(TWO*r*(xi + xj))*(xi + xj)*

          Power(Power(xi,TWO) - Power(xj,TWO),cl_float(11.0,precision)) + 

         exp(TWO*r*xj)*Power(xj,EIGHT)*

          (-cl_float(600.0,precision)*Power(r,THREE)*Power(xi,cl_float(18.0,precision)) - cl_float(30.0,precision)*Power(r,FOUR)*Power(xi,cl_float(19.0,precision)) - 

            cl_float(60.0,precision)*Power(r,THREE)*Power(xi,cl_float(16.0,precision))*Power(xj,TWO) + 

            cl_float(20.0,precision)*Power(r,FOUR)*Power(xi,cl_float(17.0,precision))*Power(xj,TWO) + cl_float(225.0,precision)*xi*Power(xj,cl_float(14.0,precision)) + 

            cl_float(360.0,precision)*r*Power(xi,TWO)*Power(xj,cl_float(14.0,precision)) + 

            cl_float(180.0,precision)*Power(r,TWO)*Power(xi,THREE)*Power(xj,cl_float(14.0,precision)) + 

            cl_float(30.0,precision)*Power(r,TWO)*Power(xi,cl_float(17.0,precision))*(-cl_float(165.0,precision) + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(60.0,precision)*r*Power(xi,cl_float(16.0,precision))*(cl_float(330.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(45.0,precision)*Power(xi,THREE)*Power(xj,cl_float(12.0,precision))*(-cl_float(55.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            r*Power(xi,NINE)*Power(xj,SIX)*

             (-cl_float(9900.0,precision)*r*Power(xj,TWO) - cl_float(136.0,precision)*Power(r,THREE)*Power(xj,FOUR)) - 

            FIVE*r*Power(xi,SEVEN)*Power(xj,EIGHT)*

             (-cl_float(2484.0,precision)*r*Power(xj,TWO) + EIGHT*Power(r,THREE)*Power(xj,FOUR)) + 

            THREE*r*Power(xi,FIVE)*Power(xj,cl_float(10.0,precision))*

             (-cl_float(660.0,precision)*r*Power(xj,TWO) + EIGHT*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(15.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(10.0,precision))*

             (-cl_float(264.0,precision)*r*Power(xj,TWO) + EIGHT*Power(r,THREE)*Power(xj,FOUR)) - 

            cl_float(165.0,precision)*Power(xi,SIX)*Power(xj,EIGHT)*

             (-cl_float(120.0,precision)*r*Power(xj,TWO) + EIGHT*Power(r,THREE)*Power(xj,FOUR)) - 

            FIVE*r*Power(xi,cl_float(13.0,precision))*Power(xj,TWO)*

             (-cl_float(6876.0,precision)*r*Power(xj,TWO) + cl_float(88.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            FIVE*r*Power(xi,cl_float(11.0,precision))*Power(xj,FOUR)*

             (-cl_float(4884.0,precision)*r*Power(xj,TWO) + cl_float(88.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(15.0,precision)*Power(xi,EIGHT)*Power(xj,SIX)*

             (-cl_float(7128.0,precision)*r*Power(xj,TWO) + cl_float(104.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            r*Power(xi,cl_float(15.0,precision))*(-cl_float(7380.0,precision)*r*Power(xj,TWO) + cl_float(136.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(15.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,FOUR)*

             (cl_float(2728.0,precision)*r*Power(xj,TWO) + cl_float(264.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(15.0,precision)*Power(xi,cl_float(14.0,precision))*(-cl_float(5864.0,precision)*r*Power(xj,TWO) + 

               cl_float(376.0,precision)*Power(r,THREE)*Power(xj,FOUR)) - 

            cl_float(15.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,TWO)*

             (-cl_float(10504.0,precision)*r*Power(xj,TWO) + cl_float(616.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            Power(xi,NINE)*Power(xj,SIX)*

             (cl_float(234135.0,precision) - cl_float(4950.0,precision)*Power(r,TWO)*Power(xj,TWO) - cl_float(34.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) 

    - FIVE*Power(xi,SEVEN)*Power(xj,EIGHT)*(cl_float(6237.0,precision) - cl_float(1242.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            THREE*Power(xi,FIVE)*Power(xj,cl_float(10.0,precision))*

             (cl_float(4125.0,precision) - cl_float(330.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) - 

            FIVE*Power(xi,cl_float(13.0,precision))*Power(xj,TWO)*

             (cl_float(43875.0,precision) - cl_float(3438.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(22.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) 

    + FIVE*Power(xi,cl_float(11.0,precision))*Power(xj,FOUR)*(cl_float(7695.0,precision) - cl_float(2442.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(22.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            Power(xi,cl_float(15.0,precision))*(-cl_float(32175.0,precision) - cl_float(3690.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(34.0,precision)*Power(r,FOUR)*Power(xj,FOUR))) + 

         TWO*exp(TWO*r*xj)*Power(xj,NINE)*

          (-cl_float(150.0,precision)*Power(r,FOUR)*Power(xi,cl_float(18.0,precision)) - SIX*Power(r,FIVE)*Power(xi,cl_float(19.0,precision)) + 

            cl_float(135.0,precision)*Power(xj,cl_float(14.0,precision)) + cl_float(225.0,precision)*r*xi*Power(xj,cl_float(14.0,precision)) + 

            cl_float(10.0,precision)*Power(r,THREE)*Power(xi,cl_float(17.0,precision))*(-cl_float(165.0,precision) + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(30.0,precision)*Power(r,TWO)*Power(xi,cl_float(16.0,precision))*(cl_float(330.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(45.0,precision)*r*Power(xi,THREE)*Power(xj,cl_float(12.0,precision))*(-cl_float(55.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(45.0,precision)*Power(xi,TWO)*Power(xj,cl_float(12.0,precision))*(-cl_float(33.0,precision) + FOUR*Power(r,TWO)*Power(xj,TWO)) + 

            r*Power(xi,NINE)*Power(xj,SIX)*

             (cl_float(234135.0,precision) - cl_float(4950.0,precision)*Power(r,TWO)*Power(xj,TWO) - cl_float(34.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) 

    - FIVE*r*Power(xi,SEVEN)*Power(xj,EIGHT)*(cl_float(6237.0,precision) - cl_float(1242.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            THREE*r*Power(xi,FIVE)*Power(xj,cl_float(10.0,precision))*

             (cl_float(4125.0,precision) - cl_float(330.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(15.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(10.0,precision))*

             (cl_float(495.0,precision) - cl_float(132.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) - 

            cl_float(165.0,precision)*Power(xi,SIX)*Power(xj,EIGHT)*

             (cl_float(135.0,precision) - cl_float(60.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) - 

            FIVE*r*Power(xi,cl_float(13.0,precision))*Power(xj,TWO)*

             (cl_float(43875.0,precision) - cl_float(3438.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(22.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) 

    + FIVE*r*Power(xi,cl_float(11.0,precision))*Power(xj,FOUR)*

             (cl_float(7695.0,precision) - cl_float(2442.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(22.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(15.0,precision)*Power(xi,EIGHT)*Power(xj,SIX)*

             (-cl_float(33.0,precision) - cl_float(3564.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(26.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            r*Power(xi,cl_float(15.0,precision))*(-cl_float(32175.0,precision) - cl_float(3690.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(34.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(15.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,FOUR)*

             (-cl_float(32277.0,precision) + cl_float(1364.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(66.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) 

    + cl_float(15.0,precision)*Power(xi,cl_float(14.0,precision))*(-cl_float(3003.0,precision) - cl_float(2932.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(94.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

            cl_float(15.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,TWO)*

             (cl_float(28119.0,precision) - cl_float(5252.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(154.0,precision)*Power(r,FOUR)*Power(xj,FOUR))) 

    + exp(TWO*r*xi)*Power(xi,EIGHT)*(-FIVE*Power(xi,TWO)*Power(xj,cl_float(12.0,precision))*

             (-cl_float(43875.0,precision)*xj - cl_float(17592.0,precision)*r*Power(xj,TWO) - cl_float(2214.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(24.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(10.0,precision)*Power(r,FOUR)*Power(xj,FIVE)) - 

            THREE*Power(xi,cl_float(14.0,precision))*(cl_float(75.0,precision)*xj + cl_float(120.0,precision)*r*Power(xj,TWO) + 

               cl_float(90.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(40.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(10.0,precision)*Power(r,FOUR)*Power(xj,FIVE)) - 

            cl_float(55.0,precision)*Power(xi,EIGHT)*Power(xj,SIX)*

             (-cl_float(567.0,precision)*xj - cl_float(1944.0,precision)*r*Power(xj,TWO) - cl_float(270.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(72.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(10.0,precision)*Power(r,FOUR)*Power(xj,FIVE)) + 

            cl_float(55.0,precision)*Power(xi,SIX)*Power(xj,EIGHT)*

             (-cl_float(4257.0,precision)*xj - cl_float(744.0,precision)*r*Power(xj,TWO) + cl_float(666.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(168.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(10.0,precision)*Power(r,FOUR)*Power(xj,FIVE)) + 

            THREE*Power(xj,cl_float(14.0,precision))*(cl_float(10725.0,precision)*xj + cl_float(6600.0,precision)*r*Power(xj,TWO) + 

               cl_float(1650.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(200.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(10.0,precision)*Power(r,FOUR)*Power(xj,FIVE)) + 

            FIVE*Power(xi,cl_float(12.0,precision))*Power(xj,TWO)*

             (cl_float(495.0,precision)*xj + cl_float(792.0,precision)*r*Power(xj,TWO) + cl_float(594.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(264.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(10.0,precision)*Power(r,FOUR)*Power(xj,FIVE)) + 

            Power(xi,cl_float(10.0,precision))*Power(xj,FOUR)*

             (-cl_float(12375.0,precision)*xj - cl_float(19800.0,precision)*r*Power(xj,TWO) - cl_float(18630.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(1560.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(170.0,precision)*Power(r,FOUR)*Power(xj,FIVE)) - 

            Power(xi,FOUR)*Power(xj,cl_float(10.0,precision))*

             (cl_float(38475.0,precision)*xj + cl_float(157560.0,precision)*r*Power(xj,TWO) + cl_float(51570.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(5640.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(170.0,precision)*Power(r,FOUR)*Power(xj,FIVE))) + 

         TWO*exp(TWO*r*xi)*Power(xi,NINE)*

          (-FIVE*Power(xi,TWO)*Power(xj,cl_float(12.0,precision))*

             (-cl_float(84357.0,precision) - cl_float(43875.0,precision)*r*xj - cl_float(8796.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(738.0,precision)*Power(r,THREE)*Power(xj,THREE) - SIX*Power(r,FOUR)*Power(xj,FOUR) + 

               TWO*Power(r,FIVE)*Power(xj,FIVE)) - 

            THREE*Power(xi,cl_float(14.0,precision))*(cl_float(45.0,precision) + cl_float(75.0,precision)*r*xj + cl_float(60.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(30.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(10.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               TWO*Power(r,FIVE)*Power(xj,FIVE)) - 

            cl_float(55.0,precision)*Power(xi,EIGHT)*Power(xj,SIX)*

             (-cl_float(405.0,precision) - cl_float(567.0,precision)*r*xj - cl_float(972.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(90.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(18.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               TWO*Power(r,FIVE)*Power(xj,FIVE)) + 

            cl_float(55.0,precision)*Power(xi,SIX)*Power(xj,EIGHT)*

             (NINE - cl_float(4257.0,precision)*r*xj - cl_float(372.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(222.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(42.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               TWO*Power(r,FIVE)*Power(xj,FIVE)) + 

            THREE*Power(xj,cl_float(14.0,precision))*(cl_float(15015.0,precision) + cl_float(10725.0,precision)*r*xj + cl_float(3300.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(550.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(50.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               TWO*Power(r,FIVE)*Power(xj,FIVE)) + 

            FIVE*Power(xi,cl_float(12.0,precision))*Power(xj,TWO)*

             (cl_float(297.0,precision) + cl_float(495.0,precision)*r*xj + cl_float(396.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(198.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(66.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               TWO*Power(r,FIVE)*Power(xj,FIVE)) + 

            Power(xi,cl_float(10.0,precision))*Power(xj,FOUR)*

             (-cl_float(7425.0,precision) - cl_float(12375.0,precision)*r*xj - cl_float(9900.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(6210.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(390.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(34.0,precision)*Power(r,FIVE)*Power(xj,FIVE)) - 

            Power(xi,FOUR)*Power(xj,cl_float(10.0,precision))*

             (-cl_float(484155.0,precision) + cl_float(38475.0,precision)*r*xj + cl_float(78780.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(17190.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(1410.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(34.0,precision)*Power(r,FIVE)*Power(xj,FIVE))))/

       (cl_float(135.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,cl_float(11.0,precision))*Power(xi + xj,cl_float(11.0,precision)))

    ; }
   
  }
  return S;
}

cl_F DSlater_3S_4S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = ZERO

    ; } else {  S = -(-cl_float(132871488750.0,precision)*xi + cl_float(149448499200.0,precision)*exp(TWO*r*xi)*xi - 

          cl_float(232588956600.0,precision)*r*Power(xi,TWO) - cl_float(200036962125.0,precision)*Power(r,TWO)*Power(xi,THREE) - 

          cl_float(112459347000.0,precision)*Power(r,THREE)*Power(xi,FOUR) - 

          cl_float(46370223900.0,precision)*Power(r,FOUR)*Power(xi,FIVE) - 

          cl_float(14905931040.0,precision)*Power(r,FIVE)*Power(xi,SIX) - 

          cl_float(3872428560.0,precision)*Power(r,SIX)*Power(xi,SEVEN) - 

          cl_float(830269440.0,precision)*Power(r,SEVEN)*Power(xi,EIGHT) - cl_float(148262400.0,precision)*Power(r,EIGHT)*Power(xi,NINE) - 

          cl_float(21964800.0,precision)*Power(r,NINE)*Power(xi,cl_float(10.0,precision)) - cl_float(2635776.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(11.0,precision)) - 

          cl_float(239616.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(12.0,precision)) - cl_float(13312.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(13.0,precision)))/

       (cl_float(7.47242496e10,precision)*exp(TWO*r*xi)*r) + 

      (-cl_float(74724249600.0,precision) + cl_float(74724249600.0,precision)*exp(TWO*r*xi) - cl_float(132871488750.0,precision)*r*xi - 

         cl_float(116294478300.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

         cl_float(66678987375.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

         cl_float(28114836750.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - 

         cl_float(9274044780.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

         cl_float(2484321840.0,precision)*Power(r,SIX)*Power(xi,SIX) - cl_float(553204080.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

         cl_float(103783680.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - cl_float(16473600.0,precision)*Power(r,NINE)*Power(xi,NINE) - 

         cl_float(2196480.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - cl_float(239616.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)) - 

         cl_float(19968.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(12.0,precision)) - cl_float(1024.0,precision)*Power(r,cl_float(13.0,precision))*Power(xi,cl_float(13.0,precision)))/

       (cl_float(7.47242496e10,precision)*exp(TWO*r*xi)*Power(r,TWO)) + 

      (xi*(-cl_float(74724249600.0,precision) + cl_float(74724249600.0,precision)*exp(TWO*r*xi) - cl_float(132871488750.0,precision)*r*xi - 

           cl_float(116294478300.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

           cl_float(66678987375.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

           cl_float(28114836750.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - 

           cl_float(9274044780.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

           cl_float(2484321840.0,precision)*Power(r,SIX)*Power(xi,SIX) - 

           cl_float(553204080.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - cl_float(103783680.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - 

           cl_float(16473600.0,precision)*Power(r,NINE)*Power(xi,NINE) - cl_float(2196480.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - 

           cl_float(239616.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)) - cl_float(19968.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(12.0,precision)) - 

           cl_float(1024.0,precision)*Power(r,cl_float(13.0,precision))*Power(xi,cl_float(13.0,precision))))/(cl_float(3.73621248e10,precision)*exp(TWO*r*xi)*r)

    ; }
 
  }
  else {
      if (r == ZERO) {  S = ZERO

    ; } else {  S = (cl_float(3780.0,precision)*exp(TWO*r*(xi + xj))*Power(Power(xi,TWO) - Power(xj,TWO),cl_float(13.0,precision)) + 

         cl_float(84.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(10.0,precision))*

          (-cl_float(60.0,precision)*Power(r,FOUR)*Power(xi,cl_float(20.0,precision)) - TWO*Power(r,FIVE)*Power(xi,cl_float(21.0,precision)) + 

            cl_float(45.0,precision)*Power(xj,cl_float(16.0,precision)) + cl_float(75.0,precision)*r*xi*Power(xj,cl_float(16.0,precision)) - 

            FOUR*Power(r,THREE)*Power(xi,cl_float(19.0,precision))*(cl_float(195.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(15.0,precision)*r*Power(xi,THREE)*Power(xj,cl_float(14.0,precision))*(-cl_float(65.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(15.0,precision)*Power(xi,TWO)*Power(xj,cl_float(14.0,precision))*(-cl_float(39.0,precision) + FOUR*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(30.0,precision)*Power(r,TWO)*Power(xi,cl_float(18.0,precision))*(cl_float(182.0,precision) + NINE*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(30.0,precision)*r*Power(xi,cl_float(13.0,precision))*Power(xj,FOUR)*

             (-cl_float(13047.0,precision) + cl_float(377.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            TWO*r*Power(xi,FIVE)*Power(xj,cl_float(12.0,precision))*

             (cl_float(2925.0,precision) - cl_float(195.0,precision)*Power(r,TWO)*Power(xj,TWO) + Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(10.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(12.0,precision))*

             (cl_float(351.0,precision) - cl_float(78.0,precision)*Power(r,TWO)*Power(xj,TWO) + Power(r,FOUR)*Power(xj,FOUR)) - 

            cl_float(130.0,precision)*Power(xi,SIX)*Power(xj,cl_float(10.0,precision))*

             (cl_float(99.0,precision) - cl_float(36.0,precision)*Power(r,TWO)*Power(xj,TWO) + Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(13.0,precision)*r*Power(xi,cl_float(11.0,precision))*Power(xj,SIX)*

             (cl_float(30735.0,precision) - cl_float(1650.0,precision)*Power(r,TWO)*Power(xj,TWO) + FOUR*Power(r,FOUR)*Power(xj,FOUR)) 

    + r*Power(xi,SEVEN)*Power(xj,cl_float(10.0,precision))*(-cl_float(15015.0,precision) + cl_float(3330.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               FOUR*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(210.0,precision)*Power(xi,cl_float(16.0,precision))*(-cl_float(156.0,precision) - cl_float(262.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               FIVE*Power(r,FOUR)*Power(xj,FOUR)) - 

            SIX*r*Power(xi,NINE)*Power(xj,EIGHT)*

             (-cl_float(48620.0,precision) - cl_float(715.0,precision)*Power(r,TWO)*Power(xj,TWO) + SIX*Power(r,FOUR)*Power(xj,FOUR)) 

    + THREE*r*Power(xi,cl_float(17.0,precision))*(-cl_float(6825.0,precision) - cl_float(1870.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(12.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

            cl_float(30.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,TWO)*

             (cl_float(17934.0,precision) - cl_float(12.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(13.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

            cl_float(15.0,precision)*Power(xi,EIGHT)*Power(xj,EIGHT)*

             (cl_float(2145.0,precision) + cl_float(2860.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(14.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) 

    + cl_float(65.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,SIX)*

             (-cl_float(13725.0,precision) - cl_float(792.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(22.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) 

    - cl_float(10.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,FOUR)*

             (cl_float(153630.0,precision) - cl_float(15054.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(143.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            Power(xi,cl_float(15.0,precision))*(-cl_float(269325.0,precision)*r*Power(xj,TWO) + 

               cl_float(9270.0,precision)*Power(r,THREE)*Power(xj,FOUR) - cl_float(52.0,precision)*Power(r,FIVE)*Power(xj,SIX))) + 

         exp(TWO*r*xi)*Power(xi,EIGHT)*

          (Power(xi,TWO)*Power(xj,cl_float(16.0,precision))*

             (cl_float(70073640.0,precision) + cl_float(47669895.0,precision)*r*xj + cl_float(13931190.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(2170350.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(169260.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(1638.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(756.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(44.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            cl_float(364.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,EIGHT)*

             (-cl_float(7425.0,precision) - cl_float(13860.0,precision)*r*xj - cl_float(5940.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(11880.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(2640.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(45.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(30.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               TWO*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            cl_float(364.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(10.0,precision))*

             (-cl_float(20925.0,precision) + cl_float(18270.0,precision)*r*xj - cl_float(58320.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(17730.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(300.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(423.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(54.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               TWO*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            THREE*Power(xi,cl_float(18.0,precision))*(cl_float(1260.0,precision) + cl_float(2205.0,precision)*r*xj + cl_float(1890.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1050.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(420.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(126.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            THREE*Power(xj,cl_float(18.0,precision))*(cl_float(1801800.0,precision) + cl_float(1576575.0,precision)*r*xj + 

               cl_float(630630.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(150150.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(23100.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(2310.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(140.0,precision)*Power(r,SIX)*Power(xj,SIX) + FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            TWO*Power(xi,cl_float(14.0,precision))*Power(xj,FOUR)*

             (-cl_float(147420.0,precision) - cl_float(257985.0,precision)*r*xj - cl_float(221130.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(122850.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(49140.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(17388.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(1512.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            cl_float(42.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,SIX)*

             (-cl_float(25740.0,precision) - cl_float(45045.0,precision)*r*xj - cl_float(38610.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(19470.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(12540.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(1836.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - EIGHT*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            cl_float(42.0,precision)*Power(xi,SIX)*Power(xj,cl_float(12.0,precision))*

             (cl_float(921600.0,precision) - cl_float(1640835.0,precision)*r*xj - cl_float(546030.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(20730.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(30180.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(5028.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(344.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            TWO*Power(xi,FOUR)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(67767840.0,precision) - cl_float(13377735.0,precision)*r*xj + cl_float(6601770.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(3115350.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(548940.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(48132.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(1848.0,precision)*Power(r,SIX)*Power(xj,SIX) + EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            Power(xi,cl_float(16.0,precision))*Power(xj,TWO)*

             (cl_float(49140.0,precision) + cl_float(85995.0,precision)*r*xj + cl_float(73710.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(40950.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(16380.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(4914.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(1092.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(44.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN))))/

       (cl_float(3780.0,precision)*exp(TWO*r*(xi + xj))*Power(r,TWO)*Power(xi - xj,cl_float(13.0,precision))*

         Power(xi + xj,cl_float(13.0,precision))) + (cl_float(3780.0,precision)*exp(TWO*r*(xi + xj))*

          Power(Power(xi,TWO) - Power(xj,TWO),cl_float(13.0,precision)) + 

         cl_float(84.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(10.0,precision))*

          (-cl_float(60.0,precision)*Power(r,FOUR)*Power(xi,cl_float(20.0,precision)) - TWO*Power(r,FIVE)*Power(xi,cl_float(21.0,precision)) + 

            cl_float(45.0,precision)*Power(xj,cl_float(16.0,precision)) + cl_float(75.0,precision)*r*xi*Power(xj,cl_float(16.0,precision)) - 

            FOUR*Power(r,THREE)*Power(xi,cl_float(19.0,precision))*(cl_float(195.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(15.0,precision)*r*Power(xi,THREE)*Power(xj,cl_float(14.0,precision))*(-cl_float(65.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(15.0,precision)*Power(xi,TWO)*Power(xj,cl_float(14.0,precision))*(-cl_float(39.0,precision) + FOUR*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(30.0,precision)*Power(r,TWO)*Power(xi,cl_float(18.0,precision))*(cl_float(182.0,precision) + NINE*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(30.0,precision)*r*Power(xi,cl_float(13.0,precision))*Power(xj,FOUR)*

             (-cl_float(13047.0,precision) + cl_float(377.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            TWO*r*Power(xi,FIVE)*Power(xj,cl_float(12.0,precision))*

             (cl_float(2925.0,precision) - cl_float(195.0,precision)*Power(r,TWO)*Power(xj,TWO) + Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(10.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(12.0,precision))*

             (cl_float(351.0,precision) - cl_float(78.0,precision)*Power(r,TWO)*Power(xj,TWO) + Power(r,FOUR)*Power(xj,FOUR)) - 

            cl_float(130.0,precision)*Power(xi,SIX)*Power(xj,cl_float(10.0,precision))*

             (cl_float(99.0,precision) - cl_float(36.0,precision)*Power(r,TWO)*Power(xj,TWO) + Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(13.0,precision)*r*Power(xi,cl_float(11.0,precision))*Power(xj,SIX)*

             (cl_float(30735.0,precision) - cl_float(1650.0,precision)*Power(r,TWO)*Power(xj,TWO) + FOUR*Power(r,FOUR)*Power(xj,FOUR)) 

    + r*Power(xi,SEVEN)*Power(xj,cl_float(10.0,precision))*(-cl_float(15015.0,precision) + cl_float(3330.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               FOUR*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(210.0,precision)*Power(xi,cl_float(16.0,precision))*(-cl_float(156.0,precision) - cl_float(262.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               FIVE*Power(r,FOUR)*Power(xj,FOUR)) - 

            SIX*r*Power(xi,NINE)*Power(xj,EIGHT)*

             (-cl_float(48620.0,precision) - cl_float(715.0,precision)*Power(r,TWO)*Power(xj,TWO) + SIX*Power(r,FOUR)*Power(xj,FOUR)) 

    + THREE*r*Power(xi,cl_float(17.0,precision))*(-cl_float(6825.0,precision) - cl_float(1870.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(12.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

            cl_float(30.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,TWO)*

             (cl_float(17934.0,precision) - cl_float(12.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(13.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

            cl_float(15.0,precision)*Power(xi,EIGHT)*Power(xj,EIGHT)*

             (cl_float(2145.0,precision) + cl_float(2860.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(14.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) 

    + cl_float(65.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,SIX)*

             (-cl_float(13725.0,precision) - cl_float(792.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(22.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) 

    - cl_float(10.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,FOUR)*

             (cl_float(153630.0,precision) - cl_float(15054.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(143.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            Power(xi,cl_float(15.0,precision))*(-cl_float(269325.0,precision)*r*Power(xj,TWO) + 

               cl_float(9270.0,precision)*Power(r,THREE)*Power(xj,FOUR) - cl_float(52.0,precision)*Power(r,FIVE)*Power(xj,SIX))) + 

         exp(TWO*r*xi)*Power(xi,EIGHT)*

          (Power(xi,TWO)*Power(xj,cl_float(16.0,precision))*

             (cl_float(70073640.0,precision) + cl_float(47669895.0,precision)*r*xj + cl_float(13931190.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(2170350.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(169260.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(1638.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(756.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(44.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            cl_float(364.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,EIGHT)*

             (-cl_float(7425.0,precision) - cl_float(13860.0,precision)*r*xj - cl_float(5940.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(11880.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(2640.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(45.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(30.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               TWO*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            cl_float(364.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(10.0,precision))*

             (-cl_float(20925.0,precision) + cl_float(18270.0,precision)*r*xj - cl_float(58320.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(17730.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(300.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(423.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(54.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               TWO*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            THREE*Power(xi,cl_float(18.0,precision))*(cl_float(1260.0,precision) + cl_float(2205.0,precision)*r*xj + cl_float(1890.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1050.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(420.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(126.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            THREE*Power(xj,cl_float(18.0,precision))*(cl_float(1801800.0,precision) + cl_float(1576575.0,precision)*r*xj + 

               cl_float(630630.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(150150.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(23100.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(2310.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(140.0,precision)*Power(r,SIX)*Power(xj,SIX) + FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            TWO*Power(xi,cl_float(14.0,precision))*Power(xj,FOUR)*

             (-cl_float(147420.0,precision) - cl_float(257985.0,precision)*r*xj - cl_float(221130.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(122850.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(49140.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(17388.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(1512.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            cl_float(42.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,SIX)*

             (-cl_float(25740.0,precision) - cl_float(45045.0,precision)*r*xj - cl_float(38610.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(19470.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(12540.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(1836.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - EIGHT*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            cl_float(42.0,precision)*Power(xi,SIX)*Power(xj,cl_float(12.0,precision))*

             (cl_float(921600.0,precision) - cl_float(1640835.0,precision)*r*xj - cl_float(546030.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(20730.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(30180.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(5028.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(344.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            TWO*Power(xi,FOUR)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(67767840.0,precision) - cl_float(13377735.0,precision)*r*xj + cl_float(6601770.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(3115350.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(548940.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(48132.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(1848.0,precision)*Power(r,SIX)*Power(xj,SIX) + EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            Power(xi,cl_float(16.0,precision))*Power(xj,TWO)*

             (cl_float(49140.0,precision) + cl_float(85995.0,precision)*r*xj + cl_float(73710.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(40950.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(16380.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(4914.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(1092.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(44.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN))))/

       (cl_float(1890.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,cl_float(13.0,precision))*Power(xi + xj,cl_float(12.0,precision))) - 

      (cl_float(7560.0,precision)*exp(TWO*r*(xi + xj))*(xi + xj)*

          Power(Power(xi,TWO) - Power(xj,TWO),cl_float(13.0,precision)) + 

         cl_float(84.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(10.0,precision))*

          (-cl_float(240.0,precision)*Power(r,THREE)*Power(xi,cl_float(20.0,precision)) - cl_float(10.0,precision)*Power(r,FOUR)*Power(xi,cl_float(21.0,precision)) - 

            cl_float(540.0,precision)*Power(r,THREE)*Power(xi,cl_float(18.0,precision))*Power(xj,TWO) - 

            EIGHT*Power(r,FOUR)*Power(xi,cl_float(19.0,precision))*Power(xj,TWO) + 

            cl_float(22620.0,precision)*Power(r,TWO)*Power(xi,cl_float(13.0,precision))*Power(xj,SIX) + cl_float(75.0,precision)*xi*Power(xj,cl_float(16.0,precision)) + 

            cl_float(120.0,precision)*r*Power(xi,TWO)*Power(xj,cl_float(16.0,precision)) + 

            cl_float(60.0,precision)*Power(r,TWO)*Power(xi,THREE)*Power(xj,cl_float(16.0,precision)) - 

            cl_float(12.0,precision)*Power(r,TWO)*Power(xi,cl_float(19.0,precision))*(cl_float(195.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(15.0,precision)*Power(xi,THREE)*Power(xj,cl_float(14.0,precision))*(-cl_float(65.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(60.0,precision)*r*Power(xi,cl_float(18.0,precision))*(cl_float(182.0,precision) + NINE*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(30.0,precision)*Power(xi,cl_float(13.0,precision))*Power(xj,FOUR)*(-cl_float(13047.0,precision) + cl_float(377.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            TWO*r*Power(xi,FIVE)*Power(xj,cl_float(12.0,precision))*

             (-cl_float(390.0,precision)*r*Power(xj,TWO) + FOUR*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(10.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(12.0,precision))*

             (-cl_float(156.0,precision)*r*Power(xj,TWO) + FOUR*Power(r,THREE)*Power(xj,FOUR)) - 

            cl_float(130.0,precision)*Power(xi,SIX)*Power(xj,cl_float(10.0,precision))*

             (-cl_float(72.0,precision)*r*Power(xj,TWO) + FOUR*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(13.0,precision)*r*Power(xi,cl_float(11.0,precision))*Power(xj,SIX)*

             (-cl_float(3300.0,precision)*r*Power(xj,TWO) + cl_float(16.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            r*Power(xi,SEVEN)*Power(xj,cl_float(10.0,precision))*

             (cl_float(6660.0,precision)*r*Power(xj,TWO) + cl_float(16.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(210.0,precision)*Power(xi,cl_float(16.0,precision))*(-cl_float(524.0,precision)*r*Power(xj,TWO) + cl_float(20.0,precision)*Power(r,THREE)*Power(xj,FOUR)) - 

            SIX*r*Power(xi,NINE)*Power(xj,EIGHT)*

             (-cl_float(1430.0,precision)*r*Power(xj,TWO) + cl_float(24.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            THREE*r*Power(xi,cl_float(17.0,precision))*(-cl_float(3740.0,precision)*r*Power(xj,TWO) + 

               cl_float(48.0,precision)*Power(r,THREE)*Power(xj,FOUR)) - 

            cl_float(30.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,TWO)*

             (-cl_float(24.0,precision)*r*Power(xj,TWO) + cl_float(52.0,precision)*Power(r,THREE)*Power(xj,FOUR)) - 

            cl_float(15.0,precision)*Power(xi,EIGHT)*Power(xj,EIGHT)*

             (cl_float(5720.0,precision)*r*Power(xj,TWO) + cl_float(56.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(65.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,SIX)*

             (-cl_float(1584.0,precision)*r*Power(xj,TWO) + cl_float(88.0,precision)*Power(r,THREE)*Power(xj,FOUR)) - 

            cl_float(10.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,FOUR)*

             (-cl_float(30108.0,precision)*r*Power(xj,TWO) + cl_float(572.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            TWO*Power(xi,FIVE)*Power(xj,cl_float(12.0,precision))*

             (cl_float(2925.0,precision) - cl_float(195.0,precision)*Power(r,TWO)*Power(xj,TWO) + Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(13.0,precision)*Power(xi,cl_float(11.0,precision))*Power(xj,SIX)*

             (cl_float(30735.0,precision) - cl_float(1650.0,precision)*Power(r,TWO)*Power(xj,TWO) + FOUR*Power(r,FOUR)*Power(xj,FOUR)) + 

            Power(xi,SEVEN)*Power(xj,cl_float(10.0,precision))*

             (-cl_float(15015.0,precision) + cl_float(3330.0,precision)*Power(r,TWO)*Power(xj,TWO) + FOUR*Power(r,FOUR)*Power(xj,FOUR)) 

    - SIX*Power(xi,NINE)*Power(xj,EIGHT)*(-cl_float(48620.0,precision) - cl_float(715.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               SIX*Power(r,FOUR)*Power(xj,FOUR)) + 

            THREE*Power(xi,cl_float(17.0,precision))*(-cl_float(6825.0,precision) - cl_float(1870.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(12.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            Power(xi,cl_float(15.0,precision))*(-cl_float(269325.0,precision)*Power(xj,TWO) + cl_float(27810.0,precision)*Power(r,TWO)*Power(xj,FOUR) - 

               cl_float(260.0,precision)*Power(r,FOUR)*Power(xj,SIX))) + 

         cl_float(168.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(11.0,precision))*

          (-cl_float(60.0,precision)*Power(r,FOUR)*Power(xi,cl_float(20.0,precision)) - TWO*Power(r,FIVE)*Power(xi,cl_float(21.0,precision)) + 

            cl_float(45.0,precision)*Power(xj,cl_float(16.0,precision)) + cl_float(75.0,precision)*r*xi*Power(xj,cl_float(16.0,precision)) - 

            FOUR*Power(r,THREE)*Power(xi,cl_float(19.0,precision))*(cl_float(195.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(15.0,precision)*r*Power(xi,THREE)*Power(xj,cl_float(14.0,precision))*(-cl_float(65.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(15.0,precision)*Power(xi,TWO)*Power(xj,cl_float(14.0,precision))*(-cl_float(39.0,precision) + FOUR*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(30.0,precision)*Power(r,TWO)*Power(xi,cl_float(18.0,precision))*(cl_float(182.0,precision) + NINE*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(30.0,precision)*r*Power(xi,cl_float(13.0,precision))*Power(xj,FOUR)*

             (-cl_float(13047.0,precision) + cl_float(377.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            TWO*r*Power(xi,FIVE)*Power(xj,cl_float(12.0,precision))*

             (cl_float(2925.0,precision) - cl_float(195.0,precision)*Power(r,TWO)*Power(xj,TWO) + Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(10.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(12.0,precision))*

             (cl_float(351.0,precision) - cl_float(78.0,precision)*Power(r,TWO)*Power(xj,TWO) + Power(r,FOUR)*Power(xj,FOUR)) - 

            cl_float(130.0,precision)*Power(xi,SIX)*Power(xj,cl_float(10.0,precision))*

             (cl_float(99.0,precision) - cl_float(36.0,precision)*Power(r,TWO)*Power(xj,TWO) + Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(13.0,precision)*r*Power(xi,cl_float(11.0,precision))*Power(xj,SIX)*

             (cl_float(30735.0,precision) - cl_float(1650.0,precision)*Power(r,TWO)*Power(xj,TWO) + FOUR*Power(r,FOUR)*Power(xj,FOUR)) + 

            r*Power(xi,SEVEN)*Power(xj,cl_float(10.0,precision))*

             (-cl_float(15015.0,precision) + cl_float(3330.0,precision)*Power(r,TWO)*Power(xj,TWO) + FOUR*Power(r,FOUR)*Power(xj,FOUR)) 

    + cl_float(210.0,precision)*Power(xi,cl_float(16.0,precision))*(-cl_float(156.0,precision) - cl_float(262.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               FIVE*Power(r,FOUR)*Power(xj,FOUR)) - 

            SIX*r*Power(xi,NINE)*Power(xj,EIGHT)*

             (-cl_float(48620.0,precision) - cl_float(715.0,precision)*Power(r,TWO)*Power(xj,TWO) + SIX*Power(r,FOUR)*Power(xj,FOUR)) + 

            THREE*r*Power(xi,cl_float(17.0,precision))*(-cl_float(6825.0,precision) - cl_float(1870.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(12.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

            cl_float(30.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,TWO)*

             (cl_float(17934.0,precision) - cl_float(12.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(13.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

            cl_float(15.0,precision)*Power(xi,EIGHT)*Power(xj,EIGHT)*

             (cl_float(2145.0,precision) + cl_float(2860.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(14.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(65.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,SIX)*

             (-cl_float(13725.0,precision) - cl_float(792.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(22.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) 

    - cl_float(10.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,FOUR)*(cl_float(153630.0,precision) - cl_float(15054.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(143.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            Power(xi,cl_float(15.0,precision))*(-cl_float(269325.0,precision)*r*Power(xj,TWO) + cl_float(9270.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(52.0,precision)*Power(r,FIVE)*Power(xj,SIX))) + 

         exp(TWO*r*xi)*Power(xi,EIGHT)*

          (Power(xi,TWO)*Power(xj,cl_float(16.0,precision))*

             (cl_float(47669895.0,precision)*xj + cl_float(27862380.0,precision)*r*Power(xj,TWO) + 

               cl_float(6511050.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(677040.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(8190.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - 

               cl_float(4536.0,precision)*Power(r,FIVE)*Power(xj,SIX) - cl_float(308.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) + 

            cl_float(364.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,EIGHT)*

             (-cl_float(13860.0,precision)*xj - cl_float(11880.0,precision)*r*Power(xj,TWO) - cl_float(35640.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(10560.0,precision)*Power(r,THREE)*Power(xj,FOUR) - cl_float(225.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(180.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(14.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) - 

            cl_float(364.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(10.0,precision))*

             (cl_float(18270.0,precision)*xj - cl_float(116640.0,precision)*r*Power(xj,TWO) - cl_float(53190.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(1200.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(2115.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(324.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(14.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) - 

            THREE*Power(xi,cl_float(18.0,precision))*(cl_float(2205.0,precision)*xj + cl_float(3780.0,precision)*r*Power(xj,TWO) + 

               cl_float(3150.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(1680.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(630.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(168.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) + 

            THREE*Power(xj,cl_float(18.0,precision))*(cl_float(1576575.0,precision)*xj + cl_float(1261260.0,precision)*r*Power(xj,TWO) + 

               cl_float(450450.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(92400.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(11550.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(840.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) + 

            TWO*Power(xi,cl_float(14.0,precision))*Power(xj,FOUR)*

             (-cl_float(257985.0,precision)*xj - cl_float(442260.0,precision)*r*Power(xj,TWO) - 

               cl_float(368550.0,precision)*Power(r,TWO)*Power(xj,THREE) - cl_float(196560.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(86940.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - cl_float(9072.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(56.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) - 

            cl_float(42.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,SIX)*

             (-cl_float(45045.0,precision)*xj - cl_float(77220.0,precision)*r*Power(xj,TWO) - cl_float(58410.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(50160.0,precision)*Power(r,THREE)*Power(xj,FOUR) - cl_float(9180.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - 

               cl_float(48.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(56.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) + 

            cl_float(42.0,precision)*Power(xi,SIX)*Power(xj,cl_float(12.0,precision))*

             (-cl_float(1640835.0,precision)*xj - cl_float(1092060.0,precision)*r*Power(xj,TWO) + 

               cl_float(62190.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(120720.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(25140.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(2064.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(56.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) - 

            TWO*Power(xi,FOUR)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(13377735.0,precision)*xj + cl_float(13203540.0,precision)*r*Power(xj,TWO) + 

               cl_float(9346050.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(2195760.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(240660.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(11088.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(56.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) + 

            Power(xi,cl_float(16.0,precision))*Power(xj,TWO)*

             (cl_float(85995.0,precision)*xj + cl_float(147420.0,precision)*r*Power(xj,TWO) + cl_float(122850.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(65520.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(24570.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(6552.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(308.0,precision)*Power(r,SIX)*Power(xj,SEVEN))) + 

         TWO*exp(TWO*r*xi)*Power(xi,NINE)*

          (Power(xi,TWO)*Power(xj,cl_float(16.0,precision))*

             (cl_float(70073640.0,precision) + cl_float(47669895.0,precision)*r*xj + cl_float(13931190.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(2170350.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(169260.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(1638.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(756.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

               cl_float(44.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            cl_float(364.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,EIGHT)*

             (-cl_float(7425.0,precision) - cl_float(13860.0,precision)*r*xj - cl_float(5940.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(11880.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(2640.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(45.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(30.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               TWO*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            cl_float(364.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(10.0,precision))*

             (-cl_float(20925.0,precision) + cl_float(18270.0,precision)*r*xj - cl_float(58320.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(17730.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(300.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(423.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(54.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               TWO*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            THREE*Power(xi,cl_float(18.0,precision))*(cl_float(1260.0,precision) + cl_float(2205.0,precision)*r*xj + cl_float(1890.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1050.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(420.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(126.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            THREE*Power(xj,cl_float(18.0,precision))*(cl_float(1801800.0,precision) + cl_float(1576575.0,precision)*r*xj + 

               cl_float(630630.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(150150.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(23100.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(2310.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(140.0,precision)*Power(r,SIX)*Power(xj,SIX) + FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            TWO*Power(xi,cl_float(14.0,precision))*Power(xj,FOUR)*

             (-cl_float(147420.0,precision) - cl_float(257985.0,precision)*r*xj - cl_float(221130.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(122850.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(49140.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(17388.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(1512.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            cl_float(42.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,SIX)*

             (-cl_float(25740.0,precision) - cl_float(45045.0,precision)*r*xj - cl_float(38610.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(19470.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(12540.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(1836.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - EIGHT*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            cl_float(42.0,precision)*Power(xi,SIX)*Power(xj,cl_float(12.0,precision))*

             (cl_float(921600.0,precision) - cl_float(1640835.0,precision)*r*xj - cl_float(546030.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(20730.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(30180.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(5028.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(344.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            TWO*Power(xi,FOUR)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(67767840.0,precision) - cl_float(13377735.0,precision)*r*xj + cl_float(6601770.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(3115350.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(548940.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(48132.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(1848.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            Power(xi,cl_float(16.0,precision))*Power(xj,TWO)*

             (cl_float(49140.0,precision) + cl_float(85995.0,precision)*r*xj + cl_float(73710.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(40950.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(16380.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(4914.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(1092.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(44.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN))))/

       (cl_float(3780.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,cl_float(13.0,precision))*Power(xi + xj,cl_float(13.0,precision)))

    ; }
   
  }
  return S;
}

cl_F DSlater_3S_5S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = ZERO

    ; } else {  S = -(-cl_float(568188982486125.0,precision)*xi + cl_float(627683696640000.0,precision)*exp(TWO*r*xi)*xi - 

          cl_float(1017388536664500.0,precision)*r*Power(xi,TWO) - 

          cl_float(899677411132500.0,precision)*Power(r,TWO)*Power(xi,THREE) - 

          cl_float(523015260768000.0,precision)*Power(r,THREE)*Power(xi,FOUR) - 

          cl_float(224405775594000.0,precision)*Power(r,FOUR)*Power(xi,FIVE) - 

          cl_float(75610821686400.0,precision)*Power(r,FIVE)*Power(xi,SIX) - 

          cl_float(20775676521600.0,precision)*Power(r,SIX)*Power(xi,SEVEN) - 

          cl_float(4769897932800.0,precision)*Power(r,SEVEN)*Power(xi,EIGHT) - 

          cl_float(929382854400.0,precision)*Power(r,EIGHT)*Power(xi,NINE) - 

          cl_float(154983628800.0,precision)*Power(r,NINE)*Power(xi,cl_float(10.0,precision)) - 

          cl_float(22140518400.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(11.0,precision)) - 

          cl_float(2683699200.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(12.0,precision)) - 

          cl_float(268369920.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(13.0,precision)) - 

          cl_float(20643840.0,precision)*Power(r,cl_float(13.0,precision))*Power(xi,cl_float(14.0,precision)) - cl_float(983040.0,precision)*Power(r,cl_float(14.0,precision))*Power(xi,cl_float(15.0,precision)))/

       (cl_float(3.1384184832e14,precision)*exp(TWO*r*xi)*r) + 

      (-cl_float(313841848320000.0,precision) + cl_float(313841848320000.0,precision)*exp(TWO*r*xi) - 

         cl_float(568188982486125.0,precision)*r*xi - cl_float(508694268332250.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

         cl_float(299892470377500.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

         cl_float(130753815192000.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - 

         cl_float(44881155118800.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

         cl_float(12601803614400.0,precision)*Power(r,SIX)*Power(xi,SIX) - 

         cl_float(2967953788800.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

         cl_float(596237241600.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - 

         cl_float(103264761600.0,precision)*Power(r,NINE)*Power(xi,NINE) - 

         cl_float(15498362880.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - 

         cl_float(2012774400.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)) - 

         cl_float(223641600.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(12.0,precision)) - 

         cl_float(20643840.0,precision)*Power(r,cl_float(13.0,precision))*Power(xi,cl_float(13.0,precision)) - cl_float(1474560.0,precision)*Power(r,cl_float(14.0,precision))*Power(xi,cl_float(14.0,precision)) - 

         cl_float(65536.0,precision)*Power(r,cl_float(15.0,precision))*Power(xi,cl_float(15.0,precision)))/

       (cl_float(3.1384184832e14,precision)*exp(TWO*r*xi)*Power(r,TWO)) + 

      (xi*(-cl_float(313841848320000.0,precision) + cl_float(313841848320000.0,precision)*exp(TWO*r*xi) - 

           cl_float(568188982486125.0,precision)*r*xi - cl_float(508694268332250.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

           cl_float(299892470377500.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

           cl_float(130753815192000.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - 

           cl_float(44881155118800.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

           cl_float(12601803614400.0,precision)*Power(r,SIX)*Power(xi,SIX) - 

           cl_float(2967953788800.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

           cl_float(596237241600.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - 

           cl_float(103264761600.0,precision)*Power(r,NINE)*Power(xi,NINE) - 

           cl_float(15498362880.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - 

           cl_float(2012774400.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)) - 

           cl_float(223641600.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(12.0,precision)) - 

           cl_float(20643840.0,precision)*Power(r,cl_float(13.0,precision))*Power(xi,cl_float(13.0,precision)) - 

           cl_float(1474560.0,precision)*Power(r,cl_float(14.0,precision))*Power(xi,cl_float(14.0,precision)) - cl_float(65536.0,precision)*Power(r,cl_float(15.0,precision))*Power(xi,cl_float(15.0,precision))))/

       (cl_float(1.5692092416e14,precision)*exp(TWO*r*xi)*r)

    ; }
 
  }
  else {
      if (r == ZERO) {  S = ZERO

    ; } else {  S = (cl_float(42525.0,precision)*exp(TWO*r*(xi + xj))*Power(Power(xi,TWO) - Power(xj,TWO),cl_float(15.0,precision)) + 

         cl_float(189.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(12.0,precision))*

          (-cl_float(350.0,precision)*Power(r,FOUR)*Power(xi,cl_float(22.0,precision)) - cl_float(10.0,precision)*Power(r,FIVE)*Power(xi,cl_float(23.0,precision)) + 

            cl_float(225.0,precision)*Power(xj,cl_float(18.0,precision)) + cl_float(375.0,precision)*r*xi*Power(xj,cl_float(18.0,precision)) - 

            cl_float(70.0,precision)*Power(r,THREE)*Power(xi,cl_float(21.0,precision))*(cl_float(75.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(75.0,precision)*r*Power(xi,THREE)*Power(xj,cl_float(16.0,precision))*(-cl_float(75.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(75.0,precision)*Power(xi,TWO)*Power(xj,cl_float(16.0,precision))*(-cl_float(45.0,precision) + FOUR*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(50.0,precision)*Power(r,TWO)*Power(xi,cl_float(20.0,precision))*(cl_float(840.0,precision) + cl_float(71.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            r*Power(xi,NINE)*Power(xj,cl_float(10.0,precision))*

             (cl_float(4694625.0,precision) + cl_float(124800.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(248.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(20.0,precision)*r*Power(xi,cl_float(17.0,precision))*Power(xj,TWO)*

             (-cl_float(185895.0,precision) - cl_float(948.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) 

    + FIVE*r*Power(xi,FIVE)*Power(xj,cl_float(14.0,precision))*

             (cl_float(7875.0,precision) - cl_float(450.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(25.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(14.0,precision))*

             (cl_float(945.0,precision) - cl_float(180.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) - 

            cl_float(375.0,precision)*Power(xi,SIX)*Power(xj,cl_float(12.0,precision))*

             (cl_float(273.0,precision) - cl_float(84.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) - 

            FIVE*r*Power(xi,cl_float(11.0,precision))*Power(xj,EIGHT)*

             (-cl_float(2803125.0,precision) + cl_float(49140.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               EIGHT*Power(r,FOUR)*Power(xj,FOUR)) + 

            FIVE*r*Power(xi,SEVEN)*Power(xj,cl_float(12.0,precision))*

             (-cl_float(16965.0,precision) + cl_float(5152.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(14.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(325.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,EIGHT)*

             (-cl_float(60117.0,precision) - cl_float(5340.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(40.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

            cl_float(15.0,precision)*r*Power(xi,cl_float(15.0,precision))*Power(xj,FOUR)*

             (cl_float(845085.0,precision) - cl_float(22960.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(52.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(15.0,precision)*r*Power(xi,cl_float(13.0,precision))*Power(xj,SIX)*

             (-cl_float(139125.0,precision) - cl_float(10140.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(52.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(75.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,SIX)*

             (-cl_float(729687.0,precision) + cl_float(25532.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(52.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(60.0,precision)*Power(xi,cl_float(18.0,precision))*(-cl_float(5355.0,precision) - cl_float(11940.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(86.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            TWO*r*Power(xi,cl_float(19.0,precision))*(-cl_float(89250.0,precision) - cl_float(35425.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(124.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(100.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,TWO)*

             (-cl_float(79713.0,precision) - cl_float(13311.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(146.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

            FIVE*Power(xi,EIGHT)*Power(xj,cl_float(10.0,precision))*

             (cl_float(157365.0,precision) + cl_float(95940.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(952.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

            cl_float(15.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,FOUR)*

             (cl_float(2638467.0,precision) - cl_float(157500.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1820.0,precision)*Power(r,FOUR)*Power(xj,FOUR))) + 

         exp(TWO*r*xi)*Power(xi,EIGHT)*

          (TWO*Power(xi,TWO)*Power(xj,cl_float(20.0,precision))*

             (cl_float(1782492075.0,precision) + cl_float(1449175455.0,precision)*r*xj + 

               cl_float(533365560.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(114631335.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(15221115.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(1142505.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(18396.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

               cl_float(5238.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - cl_float(513.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) - 

               cl_float(17.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(42.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(18.0,precision))*

             (cl_float(251336925.0,precision) + cl_float(104824125.0,precision)*r*xj + cl_float(340200.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(9122085.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(2798145.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(433755.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(39060.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

               cl_float(1890.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - cl_float(27.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               Power(r,NINE)*Power(xj,NINE)) + 

            SIX*Power(xj,cl_float(22.0,precision))*(cl_float(34459425.0,precision) + cl_float(34459425.0,precision)*r*xj + 

               cl_float(16216200.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(4729725.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(945945.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(135135.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(13860.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(990.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(45.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + Power(r,NINE)*Power(xj,NINE)) - 

            THREE*Power(xi,cl_float(22.0,precision))*(cl_float(14175.0,precision) + cl_float(25515.0,precision)*r*xj + 

               cl_float(22680.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(13230.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(5670.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(1890.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(504.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(108.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + TWO*Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(21.0,precision)*Power(xi,cl_float(18.0,precision))*Power(xj,FOUR)*

             (cl_float(212625.0,precision) + cl_float(382725.0,precision)*r*xj + cl_float(340200.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(198450.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(85050.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(28350.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(7560.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(1836.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(162.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               TWO*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(54.0,precision)*Power(xi,SIX)*Power(xj,cl_float(16.0,precision))*

             (cl_float(133451955.0,precision) - cl_float(73700865.0,precision)*r*xj - cl_float(54096840.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(8306235.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(966945.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(516747.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(80724.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(6434.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(251.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + THREE*Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(315.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,cl_float(10.0,precision))*

             (-cl_float(405405.0,precision) - cl_float(710073.0,precision)*r*xj - cl_float(805896.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(101556.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(258804.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(90972.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(9744.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(120.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(84.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               FOUR*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(315.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(12.0,precision))*

             (-cl_float(482895.0,precision) - cl_float(2656395.0,precision)*r*xj + cl_float(1186920.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1155420.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(643356.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(93492.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(336.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(1368.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(132.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + FOUR*Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(27.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,SIX)*

             (-cl_float(716625.0,precision) - cl_float(1289925.0,precision)*r*xj - cl_float(1146600.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(668850.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(286650.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(90006.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(32872.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

               cl_float(4812.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - cl_float(178.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               SIX*Power(r,NINE)*Power(xj,NINE)) + 

            Power(xi,cl_float(20.0,precision))*Power(xj,TWO)*

             (cl_float(637875.0,precision) + cl_float(1148175.0,precision)*r*xj + cl_float(1020600.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(595350.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(255150.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(85050.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(22680.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(4860.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(810.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               cl_float(34.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

            THREE*Power(xi,cl_float(14.0,precision))*Power(xj,EIGHT)*

             (-cl_float(19348875.0,precision) - cl_float(34827975.0,precision)*r*xj - cl_float(30958200.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(18689580.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(5847660.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(3723300.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(845040.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(58680.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(1548.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(236.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

            THREE*Power(xi,EIGHT)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(593408025.0,precision) + cl_float(946053675.0,precision)*r*xj - cl_float(394427880.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(315870660.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(53891460.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(910980.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(1409520.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(192168.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(11196.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               cl_float(236.0,precision)*Power(r,NINE)*Power(xj,NINE))))/

       (cl_float(42525.0,precision)*exp(TWO*r*(xi + xj))*Power(r,TWO)*Power(xi - xj,cl_float(15.0,precision))*

         Power(xi + xj,cl_float(15.0,precision))) + (TWO*(cl_float(42525.0,precision)*exp(TWO*r*(xi + xj))*

            Power(Power(xi,TWO) - Power(xj,TWO),cl_float(15.0,precision)) + 

           cl_float(189.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(12.0,precision))*

            (-cl_float(350.0,precision)*Power(r,FOUR)*Power(xi,cl_float(22.0,precision)) - cl_float(10.0,precision)*Power(r,FIVE)*Power(xi,cl_float(23.0,precision)) + 

              cl_float(225.0,precision)*Power(xj,cl_float(18.0,precision)) + cl_float(375.0,precision)*r*xi*Power(xj,cl_float(18.0,precision)) - 

              cl_float(70.0,precision)*Power(r,THREE)*Power(xi,cl_float(21.0,precision))*(cl_float(75.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

              cl_float(75.0,precision)*r*Power(xi,THREE)*Power(xj,cl_float(16.0,precision))*(-cl_float(75.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

              cl_float(75.0,precision)*Power(xi,TWO)*Power(xj,cl_float(16.0,precision))*(-cl_float(45.0,precision) + FOUR*Power(r,TWO)*Power(xj,TWO)) - 

              cl_float(50.0,precision)*Power(r,TWO)*Power(xi,cl_float(20.0,precision))*(cl_float(840.0,precision) + cl_float(71.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

              r*Power(xi,NINE)*Power(xj,cl_float(10.0,precision))*

               (cl_float(4694625.0,precision) + cl_float(124800.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(248.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

              cl_float(20.0,precision)*r*Power(xi,cl_float(17.0,precision))*Power(xj,TWO)*

               (-cl_float(185895.0,precision) - cl_float(948.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

              FIVE*r*Power(xi,FIVE)*Power(xj,cl_float(14.0,precision))*

               (cl_float(7875.0,precision) - cl_float(450.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) 

    + cl_float(25.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(14.0,precision))*

               (cl_float(945.0,precision) - cl_float(180.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) - 

              cl_float(375.0,precision)*Power(xi,SIX)*Power(xj,cl_float(12.0,precision))*

               (cl_float(273.0,precision) - cl_float(84.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) - 

              FIVE*r*Power(xi,cl_float(11.0,precision))*Power(xj,EIGHT)*

               (-cl_float(2803125.0,precision) + cl_float(49140.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 EIGHT*Power(r,FOUR)*Power(xj,FOUR)) + 

              FIVE*r*Power(xi,SEVEN)*Power(xj,cl_float(12.0,precision))*

               (-cl_float(16965.0,precision) + cl_float(5152.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(14.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

              cl_float(325.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,EIGHT)*

               (-cl_float(60117.0,precision) - cl_float(5340.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(40.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

              cl_float(15.0,precision)*r*Power(xi,cl_float(15.0,precision))*Power(xj,FOUR)*

               (cl_float(845085.0,precision) - cl_float(22960.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(52.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

              cl_float(15.0,precision)*r*Power(xi,cl_float(13.0,precision))*Power(xj,SIX)*

               (-cl_float(139125.0,precision) - cl_float(10140.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(52.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

              cl_float(75.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,SIX)*

               (-cl_float(729687.0,precision) + cl_float(25532.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(52.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

              cl_float(60.0,precision)*Power(xi,cl_float(18.0,precision))*(-cl_float(5355.0,precision) - cl_float(11940.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(86.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

              TWO*r*Power(xi,cl_float(19.0,precision))*(-cl_float(89250.0,precision) - cl_float(35425.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(124.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

              cl_float(100.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,TWO)*

               (-cl_float(79713.0,precision) - cl_float(13311.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(146.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

              FIVE*Power(xi,EIGHT)*Power(xj,cl_float(10.0,precision))*

               (cl_float(157365.0,precision) + cl_float(95940.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(952.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

              cl_float(15.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,FOUR)*

               (cl_float(2638467.0,precision) - cl_float(157500.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(1820.0,precision)*Power(r,FOUR)*Power(xj,FOUR))) + 

           exp(TWO*r*xi)*Power(xi,EIGHT)*

            (TWO*Power(xi,TWO)*Power(xj,cl_float(20.0,precision))*

               (cl_float(1782492075.0,precision) + cl_float(1449175455.0,precision)*r*xj + 

                 cl_float(533365560.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(114631335.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

                 cl_float(15221115.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(1142505.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

                 cl_float(18396.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(5238.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - 

                 cl_float(513.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) - cl_float(17.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

              cl_float(42.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(18.0,precision))*

               (cl_float(251336925.0,precision) + cl_float(104824125.0,precision)*r*xj + cl_float(340200.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(9122085.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

                 cl_float(2798145.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(433755.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

                 cl_float(39060.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(1890.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - 

                 cl_float(27.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + Power(r,NINE)*Power(xj,NINE)) + 

              SIX*Power(xj,cl_float(22.0,precision))*(cl_float(34459425.0,precision) + cl_float(34459425.0,precision)*r*xj + 

                 cl_float(16216200.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(4729725.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

                 cl_float(945945.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(135135.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

                 cl_float(13860.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(990.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

                 cl_float(45.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + Power(r,NINE)*Power(xj,NINE)) - 

              THREE*Power(xi,cl_float(22.0,precision))*(cl_float(14175.0,precision) + cl_float(25515.0,precision)*r*xj + 

                 cl_float(22680.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(13230.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

                 cl_float(5670.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(1890.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

                 cl_float(504.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(108.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

                 cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + TWO*Power(r,NINE)*Power(xj,NINE)) - 

              cl_float(21.0,precision)*Power(xi,cl_float(18.0,precision))*Power(xj,FOUR)*

               (cl_float(212625.0,precision) + cl_float(382725.0,precision)*r*xj + cl_float(340200.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(198450.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

                 cl_float(85050.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(28350.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

                 cl_float(7560.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(1836.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

                 cl_float(162.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + TWO*Power(r,NINE)*Power(xj,NINE)) + 

              cl_float(54.0,precision)*Power(xi,SIX)*Power(xj,cl_float(16.0,precision))*

               (cl_float(133451955.0,precision) - cl_float(73700865.0,precision)*r*xj - cl_float(54096840.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(8306235.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

                 cl_float(966945.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(516747.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

                 cl_float(80724.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(6434.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

                 cl_float(251.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + THREE*Power(r,NINE)*Power(xj,NINE)) - 

              cl_float(315.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,cl_float(10.0,precision))*

               (-cl_float(405405.0,precision) - cl_float(710073.0,precision)*r*xj - cl_float(805896.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(101556.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

                 cl_float(258804.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(90972.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(9744.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

                 cl_float(120.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(84.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

                 FOUR*Power(r,NINE)*Power(xj,NINE)) + 

              cl_float(315.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(12.0,precision))*

               (-cl_float(482895.0,precision) - cl_float(2656395.0,precision)*r*xj + cl_float(1186920.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(1155420.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

                 cl_float(643356.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(93492.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(336.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

                 cl_float(1368.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(132.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

                 FOUR*Power(r,NINE)*Power(xj,NINE)) - 

              cl_float(27.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,SIX)*

               (-cl_float(716625.0,precision) - cl_float(1289925.0,precision)*r*xj - cl_float(1146600.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(668850.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

                 cl_float(286650.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(90006.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(32872.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

                 cl_float(4812.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - cl_float(178.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

                 SIX*Power(r,NINE)*Power(xj,NINE)) + 

              Power(xi,cl_float(20.0,precision))*Power(xj,TWO)*

               (cl_float(637875.0,precision) + cl_float(1148175.0,precision)*r*xj + cl_float(1020600.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(595350.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

                 cl_float(255150.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(85050.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(22680.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

                 cl_float(4860.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(810.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

                 cl_float(34.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

              THREE*Power(xi,cl_float(14.0,precision))*Power(xj,EIGHT)*

               (-cl_float(19348875.0,precision) - cl_float(34827975.0,precision)*r*xj - cl_float(30958200.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(18689580.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

                 cl_float(5847660.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(3723300.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

                 cl_float(845040.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

                 cl_float(58680.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(1548.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

                 cl_float(236.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

              THREE*Power(xi,EIGHT)*Power(xj,cl_float(14.0,precision))*

               (-cl_float(593408025.0,precision) + cl_float(946053675.0,precision)*r*xj - 

                 cl_float(394427880.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(315870660.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

                 cl_float(53891460.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(910980.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

                 cl_float(1409520.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

                 cl_float(192168.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(11196.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

                 cl_float(236.0,precision)*Power(r,NINE)*Power(xj,NINE)))))/

       (cl_float(42525.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,cl_float(15.0,precision))*Power(xi + xj,cl_float(14.0,precision))) - 

      (cl_float(85050.0,precision)*exp(TWO*r*(xi + xj))*(xi + xj)*

          Power(Power(xi,TWO) - Power(xj,TWO),cl_float(15.0,precision)) + 

         cl_float(189.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(12.0,precision))*

          (-cl_float(1400.0,precision)*Power(r,THREE)*Power(xi,cl_float(22.0,precision)) - cl_float(50.0,precision)*Power(r,FOUR)*Power(xi,cl_float(23.0,precision)) - 

            cl_float(7100.0,precision)*Power(r,THREE)*Power(xi,cl_float(20.0,precision))*Power(xj,TWO) - 

            cl_float(140.0,precision)*Power(r,FOUR)*Power(xi,cl_float(21.0,precision))*Power(xj,TWO) + cl_float(375.0,precision)*xi*Power(xj,cl_float(18.0,precision)) + 

            cl_float(600.0,precision)*r*Power(xi,TWO)*Power(xj,cl_float(18.0,precision)) + 

            cl_float(300.0,precision)*Power(r,TWO)*Power(xi,THREE)*Power(xj,cl_float(18.0,precision)) - 

            cl_float(210.0,precision)*Power(r,TWO)*Power(xi,cl_float(21.0,precision))*(cl_float(75.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(75.0,precision)*Power(xi,THREE)*Power(xj,cl_float(16.0,precision))*(-cl_float(75.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(100.0,precision)*r*Power(xi,cl_float(20.0,precision))*(cl_float(840.0,precision) + cl_float(71.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            r*Power(xi,NINE)*Power(xj,cl_float(10.0,precision))*

             (cl_float(249600.0,precision)*r*Power(xj,TWO) - cl_float(992.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(20.0,precision)*r*Power(xi,cl_float(17.0,precision))*Power(xj,TWO)*

             (-cl_float(1896.0,precision)*r*Power(xj,TWO) + EIGHT*Power(r,THREE)*Power(xj,FOUR)) + 

            FIVE*r*Power(xi,FIVE)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(900.0,precision)*r*Power(xj,TWO) + EIGHT*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(25.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(360.0,precision)*r*Power(xj,TWO) + EIGHT*Power(r,THREE)*Power(xj,FOUR)) - 

            cl_float(375.0,precision)*Power(xi,SIX)*Power(xj,cl_float(12.0,precision))*

             (-cl_float(168.0,precision)*r*Power(xj,TWO) + EIGHT*Power(r,THREE)*Power(xj,FOUR)) - 

            FIVE*r*Power(xi,cl_float(11.0,precision))*Power(xj,EIGHT)*

             (cl_float(98280.0,precision)*r*Power(xj,TWO) + cl_float(32.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            FIVE*r*Power(xi,SEVEN)*Power(xj,cl_float(12.0,precision))*

             (cl_float(10304.0,precision)*r*Power(xj,TWO) + cl_float(56.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(325.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,EIGHT)*

             (-cl_float(10680.0,precision)*r*Power(xj,TWO) + cl_float(160.0,precision)*Power(r,THREE)*Power(xj,FOUR)) - 

            cl_float(15.0,precision)*r*Power(xi,cl_float(15.0,precision))*Power(xj,FOUR)*

             (-cl_float(45920.0,precision)*r*Power(xj,TWO) + cl_float(208.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(15.0,precision)*r*Power(xi,cl_float(13.0,precision))*Power(xj,SIX)*

             (-cl_float(20280.0,precision)*r*Power(xj,TWO) + cl_float(208.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(75.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,SIX)*

             (cl_float(51064.0,precision)*r*Power(xj,TWO) + cl_float(208.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(60.0,precision)*Power(xi,cl_float(18.0,precision))*(-cl_float(23880.0,precision)*r*Power(xj,TWO) + 

               cl_float(344.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            TWO*r*Power(xi,cl_float(19.0,precision))*(-cl_float(70850.0,precision)*r*Power(xj,TWO) + 

               cl_float(496.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(100.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,TWO)*

             (-cl_float(26622.0,precision)*r*Power(xj,TWO) + cl_float(584.0,precision)*Power(r,THREE)*Power(xj,FOUR)) - 

            FIVE*Power(xi,EIGHT)*Power(xj,cl_float(10.0,precision))*

             (cl_float(191880.0,precision)*r*Power(xj,TWO) + cl_float(3808.0,precision)*Power(r,THREE)*Power(xj,FOUR)) - 

            cl_float(15.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,FOUR)*

             (-cl_float(315000.0,precision)*r*Power(xj,TWO) + cl_float(7280.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            Power(xi,NINE)*Power(xj,cl_float(10.0,precision))*

             (cl_float(4694625.0,precision) + cl_float(124800.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(248.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(20.0,precision)*Power(xi,cl_float(17.0,precision))*Power(xj,TWO)*

             (-cl_float(185895.0,precision) - cl_float(948.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) 

    + FIVE*Power(xi,FIVE)*Power(xj,cl_float(14.0,precision))*(cl_float(7875.0,precision) - cl_float(450.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               TWO*Power(r,FOUR)*Power(xj,FOUR)) - 

            FIVE*Power(xi,cl_float(11.0,precision))*Power(xj,EIGHT)*

             (-cl_float(2803125.0,precision) + cl_float(49140.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               EIGHT*Power(r,FOUR)*Power(xj,FOUR)) + 

            FIVE*Power(xi,SEVEN)*Power(xj,cl_float(12.0,precision))*

             (-cl_float(16965.0,precision) + cl_float(5152.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(14.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) 

    - cl_float(15.0,precision)*Power(xi,cl_float(15.0,precision))*Power(xj,FOUR)*(cl_float(845085.0,precision) - cl_float(22960.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(52.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(15.0,precision)*Power(xi,cl_float(13.0,precision))*Power(xj,SIX)*

             (-cl_float(139125.0,precision) - cl_float(10140.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(52.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            TWO*Power(xi,cl_float(19.0,precision))*(-cl_float(89250.0,precision) - cl_float(35425.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(124.0,precision)*Power(r,FOUR)*Power(xj,FOUR))) + 

         cl_float(378.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(13.0,precision))*

          (-cl_float(350.0,precision)*Power(r,FOUR)*Power(xi,cl_float(22.0,precision)) - cl_float(10.0,precision)*Power(r,FIVE)*Power(xi,cl_float(23.0,precision)) + 

            cl_float(225.0,precision)*Power(xj,cl_float(18.0,precision)) + cl_float(375.0,precision)*r*xi*Power(xj,cl_float(18.0,precision)) - 

            cl_float(70.0,precision)*Power(r,THREE)*Power(xi,cl_float(21.0,precision))*(cl_float(75.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(75.0,precision)*r*Power(xi,THREE)*Power(xj,cl_float(16.0,precision))*(-cl_float(75.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(75.0,precision)*Power(xi,TWO)*Power(xj,cl_float(16.0,precision))*(-cl_float(45.0,precision) + FOUR*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(50.0,precision)*Power(r,TWO)*Power(xi,cl_float(20.0,precision))*(cl_float(840.0,precision) + cl_float(71.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            r*Power(xi,NINE)*Power(xj,cl_float(10.0,precision))*

             (cl_float(4694625.0,precision) + cl_float(124800.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(248.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(20.0,precision)*r*Power(xi,cl_float(17.0,precision))*Power(xj,TWO)*

             (-cl_float(185895.0,precision) - cl_float(948.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) 

    + FIVE*r*Power(xi,FIVE)*Power(xj,cl_float(14.0,precision))*

             (cl_float(7875.0,precision) - cl_float(450.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(25.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(14.0,precision))*

             (cl_float(945.0,precision) - cl_float(180.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) - 

            cl_float(375.0,precision)*Power(xi,SIX)*Power(xj,cl_float(12.0,precision))*

             (cl_float(273.0,precision) - cl_float(84.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) - 

            FIVE*r*Power(xi,cl_float(11.0,precision))*Power(xj,EIGHT)*

             (-cl_float(2803125.0,precision) + cl_float(49140.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               EIGHT*Power(r,FOUR)*Power(xj,FOUR)) + 

            FIVE*r*Power(xi,SEVEN)*Power(xj,cl_float(12.0,precision))*

             (-cl_float(16965.0,precision) + cl_float(5152.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(14.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) 

    + cl_float(325.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,EIGHT)*

             (-cl_float(60117.0,precision) - cl_float(5340.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(40.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) 

    - cl_float(15.0,precision)*r*Power(xi,cl_float(15.0,precision))*Power(xj,FOUR)*

             (cl_float(845085.0,precision) - cl_float(22960.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(52.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(15.0,precision)*r*Power(xi,cl_float(13.0,precision))*Power(xj,SIX)*

             (-cl_float(139125.0,precision) - cl_float(10140.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(52.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(75.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,SIX)*

             (-cl_float(729687.0,precision) + cl_float(25532.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(52.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(60.0,precision)*Power(xi,cl_float(18.0,precision))*(-cl_float(5355.0,precision) - cl_float(11940.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(86.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            TWO*r*Power(xi,cl_float(19.0,precision))*(-cl_float(89250.0,precision) - cl_float(35425.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(124.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(100.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,TWO)*

             (-cl_float(79713.0,precision) - cl_float(13311.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(146.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

            FIVE*Power(xi,EIGHT)*Power(xj,cl_float(10.0,precision))*

             (cl_float(157365.0,precision) + cl_float(95940.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(952.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

            cl_float(15.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,FOUR)*

             (cl_float(2638467.0,precision) - cl_float(157500.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1820.0,precision)*Power(r,FOUR)*Power(xj,FOUR))) + 

         exp(TWO*r*xi)*Power(xi,EIGHT)*

          (TWO*Power(xi,TWO)*Power(xj,cl_float(20.0,precision))*

             (cl_float(1449175455.0,precision)*xj + cl_float(1066731120.0,precision)*r*Power(xj,TWO) + 

               cl_float(343894005.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(60884460.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(5712525.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(110376.0,precision)*Power(r,FIVE)*Power(xj,SIX) - cl_float(36666.0,precision)*Power(r,SIX)*Power(xj,SEVEN) - 

               cl_float(4104.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) - cl_float(153.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            cl_float(42.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(18.0,precision))*

             (cl_float(104824125.0,precision)*xj + cl_float(680400.0,precision)*r*Power(xj,TWO) - 

               cl_float(27366255.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(11192580.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(2168775.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - 

               cl_float(234360.0,precision)*Power(r,FIVE)*Power(xj,SIX) - cl_float(13230.0,precision)*Power(r,SIX)*Power(xj,SEVEN) - 

               cl_float(216.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + NINE*Power(r,EIGHT)*Power(xj,NINE)) + 

            SIX*Power(xj,cl_float(22.0,precision))*(cl_float(34459425.0,precision)*xj + cl_float(32432400.0,precision)*r*Power(xj,TWO) + 

               cl_float(14189175.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(3783780.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(675675.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(83160.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(6930.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(360.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               NINE*Power(r,EIGHT)*Power(xj,NINE)) - 

            THREE*Power(xi,cl_float(22.0,precision))*(cl_float(25515.0,precision)*xj + cl_float(45360.0,precision)*r*Power(xj,TWO) + 

               cl_float(39690.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(22680.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(9450.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(3024.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(756.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(144.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) - 

            cl_float(21.0,precision)*Power(xi,cl_float(18.0,precision))*Power(xj,FOUR)*

             (cl_float(382725.0,precision)*xj + cl_float(680400.0,precision)*r*Power(xj,TWO) + 

               cl_float(595350.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(340200.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(141750.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(45360.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(12852.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(1296.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            cl_float(54.0,precision)*Power(xi,SIX)*Power(xj,cl_float(16.0,precision))*

             (-cl_float(73700865.0,precision)*xj - cl_float(108193680.0,precision)*r*Power(xj,TWO) - 

               cl_float(24918705.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(3867780.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(2583735.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(484344.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(45038.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + 

               cl_float(2008.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + cl_float(27.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) - 

            cl_float(315.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,cl_float(10.0,precision))*

             (-cl_float(710073.0,precision)*xj - cl_float(1611792.0,precision)*r*Power(xj,TWO) - 

               cl_float(304668.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(1035216.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(454860.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - cl_float(58464.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(840.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(672.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(36.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            cl_float(315.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(12.0,precision))*

             (-cl_float(2656395.0,precision)*xj + cl_float(2373840.0,precision)*r*Power(xj,TWO) - 

               cl_float(3466260.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(2573424.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(467460.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(2016.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(9576.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(1056.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(36.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) - 

            cl_float(27.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,SIX)*

             (-cl_float(1289925.0,precision)*xj - cl_float(2293200.0,precision)*r*Power(xj,TWO) - 

               cl_float(2006550.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(1146600.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(450030.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - cl_float(197232.0,precision)*Power(r,FIVE)*Power(xj,SIX) - 

               cl_float(33684.0,precision)*Power(r,SIX)*Power(xj,SEVEN) - cl_float(1424.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(54.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            Power(xi,cl_float(20.0,precision))*Power(xj,TWO)*

             (cl_float(1148175.0,precision)*xj + cl_float(2041200.0,precision)*r*Power(xj,TWO) + 

               cl_float(1786050.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(1020600.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(425250.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(136080.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(34020.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(6480.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(306.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            THREE*Power(xi,cl_float(14.0,precision))*Power(xj,EIGHT)*

             (-cl_float(34827975.0,precision)*xj - cl_float(61916400.0,precision)*r*Power(xj,TWO) - 

               cl_float(56068740.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(23390640.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(18616500.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - 

               cl_float(5070240.0,precision)*Power(r,FIVE)*Power(xj,SIX) - 

               cl_float(410760.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(12384.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(2124.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) - 

            THREE*Power(xi,EIGHT)*Power(xj,cl_float(14.0,precision))*

             (cl_float(946053675.0,precision)*xj - cl_float(788855760.0,precision)*r*Power(xj,TWO) - 

               cl_float(947611980.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(215565840.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(4554900.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(8457120.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(1345176.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(89568.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(2124.0,precision)*Power(r,EIGHT)*Power(xj,NINE))) + 

         TWO*exp(TWO*r*xi)*Power(xi,NINE)*

          (TWO*Power(xi,TWO)*Power(xj,cl_float(20.0,precision))*

             (cl_float(1782492075.0,precision) + cl_float(1449175455.0,precision)*r*xj + cl_float(533365560.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(114631335.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(15221115.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(1142505.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(18396.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

               cl_float(5238.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - cl_float(513.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) - 

               cl_float(17.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(42.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(18.0,precision))*

             (cl_float(251336925.0,precision) + cl_float(104824125.0,precision)*r*xj + cl_float(340200.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(9122085.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(2798145.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(433755.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(39060.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(1890.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - 

               cl_float(27.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + Power(r,NINE)*Power(xj,NINE)) + 

            SIX*Power(xj,cl_float(22.0,precision))*(cl_float(34459425.0,precision) + cl_float(34459425.0,precision)*r*xj + 

               cl_float(16216200.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(4729725.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(945945.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(135135.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(13860.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(990.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(45.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               Power(r,NINE)*Power(xj,NINE)) - 

            THREE*Power(xi,cl_float(22.0,precision))*(cl_float(14175.0,precision) + cl_float(25515.0,precision)*r*xj + cl_float(22680.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(13230.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(5670.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(1890.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(504.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(108.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               TWO*Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(21.0,precision)*Power(xi,cl_float(18.0,precision))*Power(xj,FOUR)*

             (cl_float(212625.0,precision) + cl_float(382725.0,precision)*r*xj + cl_float(340200.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(198450.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(85050.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(28350.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(7560.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(1836.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(162.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               TWO*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(54.0,precision)*Power(xi,SIX)*Power(xj,cl_float(16.0,precision))*

             (cl_float(133451955.0,precision) - cl_float(73700865.0,precision)*r*xj - cl_float(54096840.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(8306235.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(966945.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(516747.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(80724.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(6434.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(251.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               THREE*Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(315.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,cl_float(10.0,precision))*

             (-cl_float(405405.0,precision) - cl_float(710073.0,precision)*r*xj - cl_float(805896.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(101556.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(258804.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(90972.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(9744.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(120.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(84.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               FOUR*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(315.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(12.0,precision))*

             (-cl_float(482895.0,precision) - cl_float(2656395.0,precision)*r*xj + cl_float(1186920.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1155420.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(643356.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(93492.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(336.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(1368.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(132.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               FOUR*Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(27.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,SIX)*

             (-cl_float(716625.0,precision) - cl_float(1289925.0,precision)*r*xj - cl_float(1146600.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(668850.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(286650.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(90006.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(32872.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

               cl_float(4812.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - cl_float(178.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               SIX*Power(r,NINE)*Power(xj,NINE)) + 

            Power(xi,cl_float(20.0,precision))*Power(xj,TWO)*

             (cl_float(637875.0,precision) + cl_float(1148175.0,precision)*r*xj + cl_float(1020600.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(595350.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(255150.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(85050.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(22680.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(4860.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(810.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               cl_float(34.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

            THREE*Power(xi,cl_float(14.0,precision))*Power(xj,EIGHT)*

             (-cl_float(19348875.0,precision) - cl_float(34827975.0,precision)*r*xj - cl_float(30958200.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(18689580.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(5847660.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(3723300.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(845040.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

               cl_float(58680.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(1548.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               cl_float(236.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

            THREE*Power(xi,EIGHT)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(593408025.0,precision) + cl_float(946053675.0,precision)*r*xj - cl_float(394427880.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(315870660.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(53891460.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(910980.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(1409520.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(192168.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(11196.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(236.0,precision)*Power(r,NINE)*Power(xj,NINE))))/

       (cl_float(42525.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,cl_float(15.0,precision))*Power(xi + xj,cl_float(15.0,precision)))

    ; }
   
  }
  return S;
}

cl_F DSlater_3S_1S(cl_F r,cl_F xi,cl_F xj)
{
  return DSlater_1S_3S(r,xj,xi);
}

cl_F DSlater_3S_2S(cl_F r,cl_F xi,cl_F xj)
{
  return DSlater_2S_3S(r,xj,xi);
}

cl_F DSlater_4S_4S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = ZERO

    ; } else {  S = -(-cl_float(150568359566625.0,precision)*xi + cl_float(167382319104000.0,precision)*exp(TWO*r*xi)*xi - 

          cl_float(267508800058500.0,precision)*r*Power(xi,TWO) - 

          cl_float(234428725030500.0,precision)*Power(r,TWO)*Power(xi,THREE) - 

          cl_float(134962892064000.0,precision)*Power(r,THREE)*Power(xi,FOUR) - 

          cl_float(57353780484000.0,precision)*Power(r,FOUR)*Power(xi,FIVE) - 

          cl_float(19160153812800.0,precision)*Power(r,FIVE)*Power(xi,SIX) - 

          cl_float(5229789364800.0,precision)*Power(r,SIX)*Power(xi,SEVEN) - 

          cl_float(1195587993600.0,precision)*Power(r,SEVEN)*Power(xi,EIGHT) - 

          cl_float(232475443200.0,precision)*Power(r,EIGHT)*Power(xi,NINE) - 

          cl_float(38745907200.0,precision)*Power(r,NINE)*Power(xi,cl_float(10.0,precision)) - 

          cl_float(5535129600.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(11.0,precision)) - 

          cl_float(670924800.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(12.0,precision)) - 

          cl_float(67092480.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(13.0,precision)) - 

          cl_float(5160960.0,precision)*Power(r,cl_float(13.0,precision))*Power(xi,cl_float(14.0,precision)) - cl_float(245760.0,precision)*Power(r,cl_float(14.0,precision))*Power(xi,cl_float(15.0,precision)))/

       (cl_float(8.3691159552e13,precision)*exp(TWO*r*xi)*r) + 

      (-cl_float(83691159552000.0,precision) + cl_float(83691159552000.0,precision)*exp(TWO*r*xi) - 

         cl_float(150568359566625.0,precision)*r*xi - cl_float(133754400029250.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

         cl_float(78142908343500.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

         cl_float(33740723016000.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - 

         cl_float(11470756096800.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

         cl_float(3193358968800.0,precision)*Power(r,SIX)*Power(xi,SIX) - 

         cl_float(747112766400.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

         cl_float(149448499200.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - 

         cl_float(25830604800.0,precision)*Power(r,NINE)*Power(xi,NINE) - 

         cl_float(3874590720.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - 

         cl_float(503193600.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)) - 

         cl_float(55910400.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(12.0,precision)) - cl_float(5160960.0,precision)*Power(r,cl_float(13.0,precision))*Power(xi,cl_float(13.0,precision)) - 

         cl_float(368640.0,precision)*Power(r,cl_float(14.0,precision))*Power(xi,cl_float(14.0,precision)) - cl_float(16384.0,precision)*Power(r,cl_float(15.0,precision))*Power(xi,cl_float(15.0,precision)))/

       (cl_float(8.3691159552e13,precision)*exp(TWO*r*xi)*Power(r,TWO)) + 

      (xi*(-cl_float(83691159552000.0,precision) + cl_float(83691159552000.0,precision)*exp(TWO*r*xi) - 

           cl_float(150568359566625.0,precision)*r*xi - cl_float(133754400029250.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

           cl_float(78142908343500.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

           cl_float(33740723016000.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - 

           cl_float(11470756096800.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

           cl_float(3193358968800.0,precision)*Power(r,SIX)*Power(xi,SIX) - 

           cl_float(747112766400.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

           cl_float(149448499200.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - 

           cl_float(25830604800.0,precision)*Power(r,NINE)*Power(xi,NINE) - 

           cl_float(3874590720.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - 

           cl_float(503193600.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)) - 

           cl_float(55910400.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(12.0,precision)) - 

           cl_float(5160960.0,precision)*Power(r,cl_float(13.0,precision))*Power(xi,cl_float(13.0,precision)) - cl_float(368640.0,precision)*Power(r,cl_float(14.0,precision))*Power(xi,cl_float(14.0,precision)) - 

           cl_float(16384.0,precision)*Power(r,cl_float(15.0,precision))*Power(xi,cl_float(15.0,precision))))/(cl_float(4.1845579776e13,precision)*exp(TWO*r*xi)*r)

    ; }
 
  }
  else {
      if (r == ZERO) {  S = ZERO

    ; } else {  S = (cl_float(1260.0,precision)*exp(TWO*r*(xi + xj))*Power(Power(xi,TWO) - Power(xj,TWO),cl_float(15.0,precision)) + 

         exp(TWO*r*xj)*Power(xj,cl_float(10.0,precision))*

          (-cl_float(3276.0,precision)*Power(r,FIVE)*Power(xi,cl_float(25.0,precision)) - cl_float(168.0,precision)*Power(r,SIX)*Power(xi,cl_float(26.0,precision)) - 

            FOUR*Power(r,SEVEN)*Power(xi,cl_float(27.0,precision)) + cl_float(1260.0,precision)*Power(xj,cl_float(20.0,precision)) + 

            cl_float(2205.0,precision)*r*xi*Power(xj,cl_float(20.0,precision)) + 

            cl_float(1890.0,precision)*Power(xi,TWO)*Power(xj,cl_float(18.0,precision))*(-cl_float(10.0,precision) + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(420.0,precision)*Power(r,FOUR)*Power(xi,cl_float(24.0,precision))*(cl_float(91.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(525.0,precision)*r*Power(xi,THREE)*Power(xj,cl_float(18.0,precision))*(-cl_float(63.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(42.0,precision)*Power(r,THREE)*Power(xi,cl_float(23.0,precision))*

             (-cl_float(6825.0,precision) - cl_float(405.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(63.0,precision)*r*Power(xi,FIVE)*Power(xj,cl_float(16.0,precision))*

             (cl_float(3675.0,precision) - cl_float(250.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(210.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(16.0,precision))*

             (cl_float(630.0,precision) - cl_float(135.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(252.0,precision)*Power(r,TWO)*Power(xi,cl_float(22.0,precision))*

             (-cl_float(5460.0,precision) - cl_float(1225.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(17.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) 

    - cl_float(1260.0,precision)*r*Power(xi,cl_float(17.0,precision))*Power(xj,FOUR)*

             (cl_float(141729.0,precision) - cl_float(10145.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(116.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(21.0,precision)*r*Power(xi,NINE)*Power(xj,cl_float(12.0,precision))*

             (cl_float(164775.0,precision) - cl_float(18460.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(828.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(14.0,precision)*Power(xi,SIX)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(40950.0,precision) + cl_float(14175.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(450.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + TWO*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(210.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(12.0,precision))*

             (-cl_float(8190.0,precision) + cl_float(4095.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(210.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + TWO*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(42.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(10.0,precision))*

             (-cl_float(209430.0,precision) - cl_float(2925.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(8840.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + FOUR*Power(r,SIX)*Power(xj,SIX)) + 

            r*Power(xi,SEVEN)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(1003275.0,precision) + cl_float(110250.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1890.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + FOUR*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(21.0,precision)*r*Power(xi,cl_float(11.0,precision))*Power(xj,cl_float(10.0,precision))*

             (-cl_float(1033695.0,precision) - cl_float(218400.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(552.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + FOUR*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(280.0,precision)*Power(xi,cl_float(18.0,precision))*Power(xj,TWO)*

             (-cl_float(385560.0,precision) - cl_float(73953.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(2370.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + FOUR*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(35.0,precision)*r*Power(xi,cl_float(15.0,precision))*Power(xj,SIX)*

             (-cl_float(1565613.0,precision) + cl_float(359520.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(7020.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + EIGHT*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(14.0,precision)*r*Power(xi,cl_float(19.0,precision))*Power(xj,TWO)*

             (-cl_float(4980150.0,precision) + cl_float(126765.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(3852.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(20.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(630.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,SIX)*

             (cl_float(708714.0,precision) - cl_float(14385.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(2340.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(20.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(210.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,FOUR)*

             (-cl_float(2087532.0,precision) + cl_float(328491.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(11740.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(52.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(84.0,precision)*Power(xi,cl_float(20.0,precision))*(cl_float(59670.0,precision) + cl_float(236250.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(8745.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(92.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            TWO*r*Power(xi,cl_float(21.0,precision))*(cl_float(1949220.0,precision) + cl_float(1598625.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(41391.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(128.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            r*Power(xi,cl_float(13.0,precision))*Power(xj,EIGHT)*

             (cl_float(173037375.0,precision) - cl_float(2784600.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(112140.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(256.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(14.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,EIGHT)*

             (-cl_float(7260750.0,precision) - cl_float(2521935.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(19500.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(344.0,precision)*Power(r,SIX)*Power(xj,SIX))) + 

         exp(TWO*r*xi)*Power(xi,cl_float(10.0,precision))*

          (cl_float(210.0,precision)*Power(xi,TWO)*Power(xj,cl_float(18.0,precision))*

             (cl_float(514080.0,precision) + cl_float(332010.0,precision)*r*xj + cl_float(94500.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(15225.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(1470.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(81.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + TWO*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(105.0,precision)*Power(xi,cl_float(18.0,precision))*Power(xj,TWO)*

             (cl_float(180.0,precision) + cl_float(315.0,precision)*r*xj + cl_float(270.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(150.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(60.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(18.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + FOUR*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(1365.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(10.0,precision))*

             (-cl_float(6444.0,precision) + cl_float(15903.0,precision)*r*xj - cl_float(25866.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(2040.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(1080.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(180.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + EIGHT*Power(r,SIX)*Power(xj,SIX)) + 

            Power(xi,cl_float(14.0,precision))*Power(xj,SIX)*

             (cl_float(573300.0,precision) + cl_float(1003275.0,precision)*r*xj + cl_float(859950.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(387660.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(371280.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(11592.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(4816.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

               cl_float(256.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            TWO*Power(xj,cl_float(20.0,precision))*(cl_float(2506140.0,precision) + cl_float(1949220.0,precision)*r*xj + 

               cl_float(687960.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(143325.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(19110.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(1638.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(84.0,precision)*Power(r,SIX)*Power(xj,SIX) + TWO*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            cl_float(42.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(16.0,precision))*

             (-cl_float(10437660.0,precision) - cl_float(4251870.0,precision)*r*xj - cl_float(493020.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(42255.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(17490.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(1971.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(102.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               TWO*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            cl_float(21.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,FOUR)*

             (-cl_float(6300.0,precision) - cl_float(11025.0,precision)*r*xj - cl_float(9450.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(5250.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(2100.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(828.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - EIGHT*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            Power(xi,cl_float(20.0,precision))*(cl_float(1260.0,precision) + cl_float(2205.0,precision)*r*xj + cl_float(1890.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1050.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(420.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(126.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            cl_float(35.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(12.0,precision))*

             (-cl_float(2904300.0,precision) + cl_float(4943925.0,precision)*r*xj + cl_float(258930.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(359520.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(70440.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(4176.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(32.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            cl_float(35.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,EIGHT)*

             (-cl_float(49140.0,precision) - cl_float(98865.0,precision)*r*xj + cl_float(3510.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(131040.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(7800.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(3204.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(360.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            Power(xi,SIX)*Power(xj,cl_float(14.0,precision))*

             (cl_float(446489820.0,precision) - cl_float(54796455.0,precision)*r*xj - cl_float(68983110.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(12782700.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(663600.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(53928.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(7728.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(256.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN))))/

       (cl_float(1260.0,precision)*exp(TWO*r*(xi + xj))*Power(r,TWO)*Power(xi - xj,cl_float(15.0,precision))*

         Power(xi + xj,cl_float(15.0,precision))) + (cl_float(1260.0,precision)*exp(TWO*r*(xi + xj))*

          Power(Power(xi,TWO) - Power(xj,TWO),cl_float(15.0,precision)) + 

         exp(TWO*r*xj)*Power(xj,cl_float(10.0,precision))*

          (-cl_float(3276.0,precision)*Power(r,FIVE)*Power(xi,cl_float(25.0,precision)) - cl_float(168.0,precision)*Power(r,SIX)*Power(xi,cl_float(26.0,precision)) - 

            FOUR*Power(r,SEVEN)*Power(xi,cl_float(27.0,precision)) + cl_float(1260.0,precision)*Power(xj,cl_float(20.0,precision)) + 

            cl_float(2205.0,precision)*r*xi*Power(xj,cl_float(20.0,precision)) + 

            cl_float(1890.0,precision)*Power(xi,TWO)*Power(xj,cl_float(18.0,precision))*(-cl_float(10.0,precision) + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(420.0,precision)*Power(r,FOUR)*Power(xi,cl_float(24.0,precision))*(cl_float(91.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(525.0,precision)*r*Power(xi,THREE)*Power(xj,cl_float(18.0,precision))*(-cl_float(63.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(42.0,precision)*Power(r,THREE)*Power(xi,cl_float(23.0,precision))*

             (-cl_float(6825.0,precision) - cl_float(405.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(63.0,precision)*r*Power(xi,FIVE)*Power(xj,cl_float(16.0,precision))*

             (cl_float(3675.0,precision) - cl_float(250.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(210.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(16.0,precision))*

             (cl_float(630.0,precision) - cl_float(135.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(252.0,precision)*Power(r,TWO)*Power(xi,cl_float(22.0,precision))*

             (-cl_float(5460.0,precision) - cl_float(1225.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(17.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) 

    - cl_float(1260.0,precision)*r*Power(xi,cl_float(17.0,precision))*Power(xj,FOUR)*

             (cl_float(141729.0,precision) - cl_float(10145.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(116.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(21.0,precision)*r*Power(xi,NINE)*Power(xj,cl_float(12.0,precision))*

             (cl_float(164775.0,precision) - cl_float(18460.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(828.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(14.0,precision)*Power(xi,SIX)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(40950.0,precision) + cl_float(14175.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(450.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + TWO*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(210.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(12.0,precision))*

             (-cl_float(8190.0,precision) + cl_float(4095.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(210.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + TWO*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(42.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(10.0,precision))*

             (-cl_float(209430.0,precision) - cl_float(2925.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(8840.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + FOUR*Power(r,SIX)*Power(xj,SIX)) + 

            r*Power(xi,SEVEN)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(1003275.0,precision) + cl_float(110250.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1890.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + FOUR*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(21.0,precision)*r*Power(xi,cl_float(11.0,precision))*Power(xj,cl_float(10.0,precision))*

             (-cl_float(1033695.0,precision) - cl_float(218400.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(552.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + FOUR*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(280.0,precision)*Power(xi,cl_float(18.0,precision))*Power(xj,TWO)*

             (-cl_float(385560.0,precision) - cl_float(73953.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(2370.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + FOUR*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(35.0,precision)*r*Power(xi,cl_float(15.0,precision))*Power(xj,SIX)*

             (-cl_float(1565613.0,precision) + cl_float(359520.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(7020.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + EIGHT*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(14.0,precision)*r*Power(xi,cl_float(19.0,precision))*Power(xj,TWO)*

             (-cl_float(4980150.0,precision) + cl_float(126765.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(3852.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(20.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(630.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,SIX)*

             (cl_float(708714.0,precision) - cl_float(14385.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(2340.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(20.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(210.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,FOUR)*

             (-cl_float(2087532.0,precision) + cl_float(328491.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(11740.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(52.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(84.0,precision)*Power(xi,cl_float(20.0,precision))*(cl_float(59670.0,precision) + cl_float(236250.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(8745.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(92.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            TWO*r*Power(xi,cl_float(21.0,precision))*(cl_float(1949220.0,precision) + cl_float(1598625.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(41391.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(128.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            r*Power(xi,cl_float(13.0,precision))*Power(xj,EIGHT)*

             (cl_float(173037375.0,precision) - cl_float(2784600.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(112140.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(256.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(14.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,EIGHT)*

             (-cl_float(7260750.0,precision) - cl_float(2521935.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(19500.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(344.0,precision)*Power(r,SIX)*Power(xj,SIX))) + 

         exp(TWO*r*xi)*Power(xi,cl_float(10.0,precision))*

          (cl_float(210.0,precision)*Power(xi,TWO)*Power(xj,cl_float(18.0,precision))*

             (cl_float(514080.0,precision) + cl_float(332010.0,precision)*r*xj + cl_float(94500.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(15225.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(1470.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(81.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + TWO*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(105.0,precision)*Power(xi,cl_float(18.0,precision))*Power(xj,TWO)*

             (cl_float(180.0,precision) + cl_float(315.0,precision)*r*xj + cl_float(270.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(150.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(60.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(18.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + FOUR*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(1365.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(10.0,precision))*

             (-cl_float(6444.0,precision) + cl_float(15903.0,precision)*r*xj - cl_float(25866.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(2040.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(1080.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(180.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + EIGHT*Power(r,SIX)*Power(xj,SIX)) + 

            Power(xi,cl_float(14.0,precision))*Power(xj,SIX)*

             (cl_float(573300.0,precision) + cl_float(1003275.0,precision)*r*xj + cl_float(859950.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(387660.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(371280.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(11592.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(4816.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

               cl_float(256.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            TWO*Power(xj,cl_float(20.0,precision))*(cl_float(2506140.0,precision) + cl_float(1949220.0,precision)*r*xj + 

               cl_float(687960.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(143325.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(19110.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(1638.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(84.0,precision)*Power(r,SIX)*Power(xj,SIX) + TWO*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            cl_float(42.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(16.0,precision))*

             (-cl_float(10437660.0,precision) - cl_float(4251870.0,precision)*r*xj - cl_float(493020.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(42255.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(17490.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(1971.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(102.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               TWO*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            cl_float(21.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,FOUR)*

             (-cl_float(6300.0,precision) - cl_float(11025.0,precision)*r*xj - cl_float(9450.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(5250.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(2100.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(828.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - EIGHT*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            Power(xi,cl_float(20.0,precision))*(cl_float(1260.0,precision) + cl_float(2205.0,precision)*r*xj + cl_float(1890.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1050.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(420.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(126.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            cl_float(35.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(12.0,precision))*

             (-cl_float(2904300.0,precision) + cl_float(4943925.0,precision)*r*xj + cl_float(258930.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(359520.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(70440.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(4176.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(32.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            cl_float(35.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,EIGHT)*

             (-cl_float(49140.0,precision) - cl_float(98865.0,precision)*r*xj + cl_float(3510.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(131040.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(7800.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(3204.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(360.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            Power(xi,SIX)*Power(xj,cl_float(14.0,precision))*

             (cl_float(446489820.0,precision) - cl_float(54796455.0,precision)*r*xj - cl_float(68983110.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(12782700.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(663600.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(53928.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(7728.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(256.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN))))/

       (cl_float(630.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,cl_float(15.0,precision))*Power(xi + xj,cl_float(14.0,precision))) - 

      (cl_float(2520.0,precision)*exp(TWO*r*(xi + xj))*(xi + xj)*

          Power(Power(xi,TWO) - Power(xj,TWO),cl_float(15.0,precision)) + 

         exp(TWO*r*xj)*Power(xj,cl_float(10.0,precision))*

          (-cl_float(16380.0,precision)*Power(r,FOUR)*Power(xi,cl_float(25.0,precision)) - cl_float(1008.0,precision)*Power(r,FIVE)*Power(xi,cl_float(26.0,precision)) - 

            cl_float(28.0,precision)*Power(r,SIX)*Power(xi,cl_float(27.0,precision)) - 

            cl_float(840.0,precision)*Power(r,FIVE)*Power(xi,cl_float(24.0,precision))*Power(xj,TWO) + cl_float(2205.0,precision)*xi*Power(xj,cl_float(20.0,precision)) + 

            cl_float(3780.0,precision)*r*Power(xi,TWO)*Power(xj,cl_float(20.0,precision)) + 

            cl_float(2100.0,precision)*Power(r,TWO)*Power(xi,THREE)*Power(xj,cl_float(20.0,precision)) - 

            cl_float(1680.0,precision)*Power(r,THREE)*Power(xi,cl_float(24.0,precision))*(cl_float(91.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(525.0,precision)*Power(xi,THREE)*Power(xj,cl_float(18.0,precision))*(-cl_float(63.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(42.0,precision)*Power(r,THREE)*Power(xi,cl_float(23.0,precision))*

             (-cl_float(810.0,precision)*r*Power(xj,TWO) + EIGHT*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(63.0,precision)*r*Power(xi,FIVE)*Power(xj,cl_float(16.0,precision))*

             (-cl_float(500.0,precision)*r*Power(xj,TWO) + EIGHT*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(210.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(16.0,precision))*

             (-cl_float(270.0,precision)*r*Power(xj,TWO) + EIGHT*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(252.0,precision)*Power(r,TWO)*Power(xi,cl_float(22.0,precision))*

             (-cl_float(2450.0,precision)*r*Power(xj,TWO) + cl_float(68.0,precision)*Power(r,THREE)*Power(xj,FOUR)) - 

            cl_float(1260.0,precision)*r*Power(xi,cl_float(17.0,precision))*Power(xj,FOUR)*

             (-cl_float(20290.0,precision)*r*Power(xj,TWO) + cl_float(464.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(21.0,precision)*r*Power(xi,NINE)*Power(xj,cl_float(12.0,precision))*

             (-cl_float(36920.0,precision)*r*Power(xj,TWO) + cl_float(3312.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(126.0,precision)*Power(r,TWO)*Power(xi,cl_float(23.0,precision))*

             (-cl_float(6825.0,precision) - cl_float(405.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(63.0,precision)*Power(xi,FIVE)*Power(xj,cl_float(16.0,precision))*

             (cl_float(3675.0,precision) - cl_float(250.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(504.0,precision)*r*Power(xi,cl_float(22.0,precision))*(-cl_float(5460.0,precision) - cl_float(1225.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(17.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

            cl_float(1260.0,precision)*Power(xi,cl_float(17.0,precision))*Power(xj,FOUR)*

             (cl_float(141729.0,precision) - cl_float(10145.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(116.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(21.0,precision)*Power(xi,NINE)*Power(xj,cl_float(12.0,precision))*

             (cl_float(164775.0,precision) - cl_float(18460.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(828.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(14.0,precision)*Power(xi,SIX)*Power(xj,cl_float(14.0,precision))*

             (cl_float(28350.0,precision)*r*Power(xj,TWO) - cl_float(1800.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(12.0,precision)*Power(r,FIVE)*Power(xj,SIX)) - 

            cl_float(210.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(12.0,precision))*

             (cl_float(8190.0,precision)*r*Power(xj,TWO) - cl_float(840.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(12.0,precision)*Power(r,FIVE)*Power(xj,SIX)) + 

            cl_float(42.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(10.0,precision))*

             (-cl_float(5850.0,precision)*r*Power(xj,TWO) - cl_float(35360.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(24.0,precision)*Power(r,FIVE)*Power(xj,SIX)) + 

            r*Power(xi,SEVEN)*Power(xj,cl_float(14.0,precision))*

             (cl_float(220500.0,precision)*r*Power(xj,TWO) - cl_float(7560.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(24.0,precision)*Power(r,FIVE)*Power(xj,SIX)) - 

            cl_float(21.0,precision)*r*Power(xi,cl_float(11.0,precision))*Power(xj,cl_float(10.0,precision))*

             (-cl_float(436800.0,precision)*r*Power(xj,TWO) + cl_float(2208.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(24.0,precision)*Power(r,FIVE)*Power(xj,SIX)) + 

            cl_float(280.0,precision)*Power(xi,cl_float(18.0,precision))*Power(xj,TWO)*

             (-cl_float(147906.0,precision)*r*Power(xj,TWO) + cl_float(9480.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(24.0,precision)*Power(r,FIVE)*Power(xj,SIX)) - 

            cl_float(35.0,precision)*r*Power(xi,cl_float(15.0,precision))*Power(xj,SIX)*

             (cl_float(719040.0,precision)*r*Power(xj,TWO) - cl_float(28080.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(48.0,precision)*Power(r,FIVE)*Power(xj,SIX)) + 

            cl_float(14.0,precision)*r*Power(xi,cl_float(19.0,precision))*Power(xj,TWO)*

             (cl_float(253530.0,precision)*r*Power(xj,TWO) - cl_float(15408.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(120.0,precision)*Power(r,FIVE)*Power(xj,SIX)) - 

            cl_float(630.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,SIX)*

             (-cl_float(28770.0,precision)*r*Power(xj,TWO) - cl_float(9360.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(120.0,precision)*Power(r,FIVE)*Power(xj,SIX)) + 

            cl_float(210.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,FOUR)*

             (cl_float(656982.0,precision)*r*Power(xj,TWO) - cl_float(46960.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(312.0,precision)*Power(r,FIVE)*Power(xj,SIX)) - 

            cl_float(84.0,precision)*Power(xi,cl_float(20.0,precision))*(cl_float(472500.0,precision)*r*Power(xj,TWO) - 

               cl_float(34980.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(552.0,precision)*Power(r,FIVE)*Power(xj,SIX)) - 

            TWO*r*Power(xi,cl_float(21.0,precision))*(cl_float(3197250.0,precision)*r*Power(xj,TWO) - 

               cl_float(165564.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(768.0,precision)*Power(r,FIVE)*Power(xj,SIX)) + 

            r*Power(xi,cl_float(13.0,precision))*Power(xj,EIGHT)*

             (-cl_float(5569200.0,precision)*r*Power(xj,TWO) - cl_float(448560.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(1536.0,precision)*Power(r,FIVE)*Power(xj,SIX)) + 

            cl_float(14.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,EIGHT)*

             (-cl_float(5043870.0,precision)*r*Power(xj,TWO) + cl_float(78000.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(2064.0,precision)*Power(r,FIVE)*Power(xj,SIX)) + 

            Power(xi,SEVEN)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(1003275.0,precision) + cl_float(110250.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1890.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + FOUR*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(21.0,precision)*Power(xi,cl_float(11.0,precision))*Power(xj,cl_float(10.0,precision))*

             (-cl_float(1033695.0,precision) - cl_float(218400.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(552.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + FOUR*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(35.0,precision)*Power(xi,cl_float(15.0,precision))*Power(xj,SIX)*

             (-cl_float(1565613.0,precision) + cl_float(359520.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(7020.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + EIGHT*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(14.0,precision)*Power(xi,cl_float(19.0,precision))*Power(xj,TWO)*

             (-cl_float(4980150.0,precision) + cl_float(126765.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(3852.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(20.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            TWO*Power(xi,cl_float(21.0,precision))*(cl_float(1949220.0,precision) + cl_float(1598625.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(41391.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(128.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            Power(xi,cl_float(13.0,precision))*Power(xj,EIGHT)*

             (cl_float(173037375.0,precision) - cl_float(2784600.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(112140.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(256.0,precision)*Power(r,SIX)*Power(xj,SIX))) + 

         TWO*exp(TWO*r*xj)*Power(xj,cl_float(11.0,precision))*

          (-cl_float(3276.0,precision)*Power(r,FIVE)*Power(xi,cl_float(25.0,precision)) - cl_float(168.0,precision)*Power(r,SIX)*Power(xi,cl_float(26.0,precision)) - 

            FOUR*Power(r,SEVEN)*Power(xi,cl_float(27.0,precision)) + cl_float(1260.0,precision)*Power(xj,cl_float(20.0,precision)) + 

            cl_float(2205.0,precision)*r*xi*Power(xj,cl_float(20.0,precision)) + 

            cl_float(1890.0,precision)*Power(xi,TWO)*Power(xj,cl_float(18.0,precision))*(-cl_float(10.0,precision) + Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(420.0,precision)*Power(r,FOUR)*Power(xi,cl_float(24.0,precision))*(cl_float(91.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(525.0,precision)*r*Power(xi,THREE)*Power(xj,cl_float(18.0,precision))*(-cl_float(63.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(42.0,precision)*Power(r,THREE)*Power(xi,cl_float(23.0,precision))*

             (-cl_float(6825.0,precision) - cl_float(405.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(63.0,precision)*r*Power(xi,FIVE)*Power(xj,cl_float(16.0,precision))*

             (cl_float(3675.0,precision) - cl_float(250.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(210.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(16.0,precision))*

             (cl_float(630.0,precision) - cl_float(135.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(252.0,precision)*Power(r,TWO)*Power(xi,cl_float(22.0,precision))*

             (-cl_float(5460.0,precision) - cl_float(1225.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(17.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) 

    - cl_float(1260.0,precision)*r*Power(xi,cl_float(17.0,precision))*Power(xj,FOUR)*

             (cl_float(141729.0,precision) - cl_float(10145.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(116.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(21.0,precision)*r*Power(xi,NINE)*Power(xj,cl_float(12.0,precision))*

             (cl_float(164775.0,precision) - cl_float(18460.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(828.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(14.0,precision)*Power(xi,SIX)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(40950.0,precision) + cl_float(14175.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(450.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + TWO*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(210.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(12.0,precision))*

             (-cl_float(8190.0,precision) + cl_float(4095.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(210.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + TWO*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(42.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(10.0,precision))*

             (-cl_float(209430.0,precision) - cl_float(2925.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(8840.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + FOUR*Power(r,SIX)*Power(xj,SIX)) + 

            r*Power(xi,SEVEN)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(1003275.0,precision) + cl_float(110250.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1890.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + FOUR*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(21.0,precision)*r*Power(xi,cl_float(11.0,precision))*Power(xj,cl_float(10.0,precision))*

             (-cl_float(1033695.0,precision) - cl_float(218400.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(552.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + FOUR*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(280.0,precision)*Power(xi,cl_float(18.0,precision))*Power(xj,TWO)*

             (-cl_float(385560.0,precision) - cl_float(73953.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(2370.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + FOUR*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(35.0,precision)*r*Power(xi,cl_float(15.0,precision))*Power(xj,SIX)*

             (-cl_float(1565613.0,precision) + cl_float(359520.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(7020.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + EIGHT*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(14.0,precision)*r*Power(xi,cl_float(19.0,precision))*Power(xj,TWO)*

             (-cl_float(4980150.0,precision) + cl_float(126765.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(3852.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(20.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(630.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,SIX)*

             (cl_float(708714.0,precision) - cl_float(14385.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(2340.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(20.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(210.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,FOUR)*

             (-cl_float(2087532.0,precision) + cl_float(328491.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(11740.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(52.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(84.0,precision)*Power(xi,cl_float(20.0,precision))*(cl_float(59670.0,precision) + cl_float(236250.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(8745.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(92.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            TWO*r*Power(xi,cl_float(21.0,precision))*(cl_float(1949220.0,precision) + cl_float(1598625.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(41391.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(128.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            r*Power(xi,cl_float(13.0,precision))*Power(xj,EIGHT)*

             (cl_float(173037375.0,precision) - cl_float(2784600.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(112140.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(256.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(14.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,EIGHT)*

             (-cl_float(7260750.0,precision) - cl_float(2521935.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(19500.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(344.0,precision)*Power(r,SIX)*Power(xj,SIX))) + 

         exp(TWO*r*xi)*Power(xi,cl_float(10.0,precision))*

          (cl_float(210.0,precision)*Power(xi,TWO)*Power(xj,cl_float(18.0,precision))*

             (cl_float(332010.0,precision)*xj + cl_float(189000.0,precision)*r*Power(xj,TWO) + 

               cl_float(45675.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(5880.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(405.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(12.0,precision)*Power(r,FIVE)*Power(xj,SIX)) + 

            cl_float(105.0,precision)*Power(xi,cl_float(18.0,precision))*Power(xj,TWO)*

             (cl_float(315.0,precision)*xj + cl_float(540.0,precision)*r*Power(xj,TWO) + cl_float(450.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(240.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(90.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(24.0,precision)*Power(r,FIVE)*Power(xj,SIX)) - 

            cl_float(1365.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(10.0,precision))*

             (cl_float(15903.0,precision)*xj - cl_float(51732.0,precision)*r*Power(xj,TWO) - cl_float(6120.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(4320.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(900.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(48.0,precision)*Power(r,FIVE)*Power(xj,SIX)) + 

            Power(xi,cl_float(14.0,precision))*Power(xj,SIX)*

             (cl_float(1003275.0,precision)*xj + cl_float(1719900.0,precision)*r*Power(xj,TWO) + 

               cl_float(1162980.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(1485120.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(57960.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - 

               cl_float(28896.0,precision)*Power(r,FIVE)*Power(xj,SIX) - cl_float(1792.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) + 

            TWO*Power(xj,cl_float(20.0,precision))*(cl_float(1949220.0,precision)*xj + cl_float(1375920.0,precision)*r*Power(xj,TWO) + 

               cl_float(429975.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(76440.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(8190.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(504.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(14.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) - 

            cl_float(42.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(16.0,precision))*

             (-cl_float(4251870.0,precision)*xj - cl_float(986040.0,precision)*r*Power(xj,TWO) + 

               cl_float(126765.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(69960.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(9855.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(612.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(14.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) + 

            cl_float(21.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,FOUR)*

             (-cl_float(11025.0,precision)*xj - cl_float(18900.0,precision)*r*Power(xj,TWO) - cl_float(15750.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(8400.0,precision)*Power(r,THREE)*Power(xj,FOUR) - cl_float(4140.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - 

               cl_float(48.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) - 

            Power(xi,cl_float(20.0,precision))*(cl_float(2205.0,precision)*xj + cl_float(3780.0,precision)*r*Power(xj,TWO) + 

               cl_float(3150.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(1680.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(630.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(168.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) - 

            cl_float(35.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(12.0,precision))*

             (cl_float(4943925.0,precision)*xj + cl_float(517860.0,precision)*r*Power(xj,TWO) - 

               cl_float(1078560.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(281760.0,precision)*Power(r,THREE)*Power(xj,FOUR) - cl_float(20880.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(192.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(56.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) + 

            cl_float(35.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,EIGHT)*

             (-cl_float(98865.0,precision)*xj + cl_float(7020.0,precision)*r*Power(xj,TWO) - cl_float(393120.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(31200.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(16020.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(2160.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(56.0,precision)*Power(r,SIX)*Power(xj,SEVEN)) + 

            Power(xi,SIX)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(54796455.0,precision)*xj - cl_float(137966220.0,precision)*r*Power(xj,TWO) - 

               cl_float(38348100.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(2654400.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(269640.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(46368.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(1792.0,precision)*Power(r,SIX)*Power(xj,SEVEN))) + 

         TWO*exp(TWO*r*xi)*Power(xi,cl_float(11.0,precision))*

          (cl_float(210.0,precision)*Power(xi,TWO)*Power(xj,cl_float(18.0,precision))*

             (cl_float(514080.0,precision) + cl_float(332010.0,precision)*r*xj + cl_float(94500.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(15225.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(1470.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(81.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + TWO*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(105.0,precision)*Power(xi,cl_float(18.0,precision))*Power(xj,TWO)*

             (cl_float(180.0,precision) + cl_float(315.0,precision)*r*xj + cl_float(270.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(150.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(60.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(18.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + FOUR*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(1365.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(10.0,precision))*

             (-cl_float(6444.0,precision) + cl_float(15903.0,precision)*r*xj - cl_float(25866.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(2040.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(1080.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(180.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + EIGHT*Power(r,SIX)*Power(xj,SIX)) + 

            Power(xi,cl_float(14.0,precision))*Power(xj,SIX)*

             (cl_float(573300.0,precision) + cl_float(1003275.0,precision)*r*xj + cl_float(859950.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(387660.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(371280.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(11592.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(4816.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

               cl_float(256.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            TWO*Power(xj,cl_float(20.0,precision))*(cl_float(2506140.0,precision) + cl_float(1949220.0,precision)*r*xj + 

               cl_float(687960.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(143325.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(19110.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(1638.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(84.0,precision)*Power(r,SIX)*Power(xj,SIX) + TWO*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            cl_float(42.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(16.0,precision))*

             (-cl_float(10437660.0,precision) - cl_float(4251870.0,precision)*r*xj - cl_float(493020.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(42255.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(17490.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(1971.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(102.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               TWO*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            cl_float(21.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,FOUR)*

             (-cl_float(6300.0,precision) - cl_float(11025.0,precision)*r*xj - cl_float(9450.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(5250.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(2100.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(828.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - EIGHT*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            Power(xi,cl_float(20.0,precision))*(cl_float(1260.0,precision) + cl_float(2205.0,precision)*r*xj + cl_float(1890.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1050.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(420.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(126.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(28.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               FOUR*Power(r,SEVEN)*Power(xj,SEVEN)) - 

            cl_float(35.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(12.0,precision))*

             (-cl_float(2904300.0,precision) + cl_float(4943925.0,precision)*r*xj + cl_float(258930.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(359520.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(70440.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(4176.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(32.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            cl_float(35.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,EIGHT)*

             (-cl_float(49140.0,precision) - cl_float(98865.0,precision)*r*xj + cl_float(3510.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(131040.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(7800.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(3204.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(360.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               EIGHT*Power(r,SEVEN)*Power(xj,SEVEN)) + 

            Power(xi,SIX)*Power(xj,cl_float(14.0,precision))*

             (cl_float(446489820.0,precision) - cl_float(54796455.0,precision)*r*xj - cl_float(68983110.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(12782700.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(663600.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(53928.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(7728.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(256.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN))))/

       (cl_float(1260.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,cl_float(15.0,precision))*Power(xi + xj,cl_float(15.0,precision)))

    ; }
   
  }
  return S;
}

cl_F DSlater_4S_5S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = ZERO

    ; } else {  S = -(-cl_float(25913502934444125.0,precision)*xi + cl_float(28454994247680000.0,precision)*exp(TWO*r*xi)*xi - 

          cl_float(46744023242416500.0,precision)*r*Power(xi,TWO) - 

          cl_float(41723129607909750.0,precision)*Power(r,TWO)*Power(xi,THREE) - 

          cl_float(24550942638222000.0,precision)*Power(r,THREE)*Power(xi,FOUR) - 

          cl_float(10704286944351000.0,precision)*Power(r,FOUR)*Power(xi,FIVE) - 

          cl_float(3684699450432000.0,precision)*Power(r,FIVE)*Power(xi,SIX) - 

          cl_float(1041667066440000.0,precision)*Power(r,SIX)*Power(xi,SEVEN) - 

          cl_float(248293113868800.0,precision)*Power(r,SEVEN)*Power(xi,EIGHT) - 

          cl_float(50808078921600.0,precision)*Power(r,EIGHT)*Power(xi,NINE) - 

          cl_float(9033331507200.0,precision)*Power(r,NINE)*Power(xi,cl_float(10.0,precision)) - 

          cl_float(1405184901120.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(11.0,precision)) - 

          cl_float(191616122880.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(12.0,precision)) - 

          cl_float(22811443200.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(13.0,precision)) - 

          cl_float(2339635200.0,precision)*Power(r,cl_float(13.0,precision))*Power(xi,cl_float(14.0,precision)) - 

          cl_float(200540160.0,precision)*Power(r,cl_float(14.0,precision))*Power(xi,cl_float(15.0,precision)) - 

          cl_float(13369344.0,precision)*Power(r,cl_float(15.0,precision))*Power(xi,cl_float(16.0,precision)) - cl_float(557056.0,precision)*Power(r,cl_float(16.0,precision))*Power(xi,cl_float(17.0,precision)))/

       (cl_float(1.422749712384e16,precision)*exp(TWO*r*xi)*r) + 

      (-cl_float(14227497123840000.0,precision) + cl_float(14227497123840000.0,precision)*exp(TWO*r*xi) - 

         cl_float(25913502934444125.0,precision)*r*xi - cl_float(23372011621208250.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

         cl_float(13907709869303250.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

         cl_float(6137735659555500.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - 

         cl_float(2140857388870200.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

         cl_float(614116575072000.0,precision)*Power(r,SIX)*Power(xi,SIX) - 

         cl_float(148809580920000.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

         cl_float(31036639233600.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - 

         cl_float(5645342102400.0,precision)*Power(r,NINE)*Power(xi,NINE) - 

         cl_float(903333150720.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - 

         cl_float(127744081920.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)) - 

         cl_float(15968010240.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(12.0,precision)) - 

         cl_float(1754726400.0,precision)*Power(r,cl_float(13.0,precision))*Power(xi,cl_float(13.0,precision)) - 

         cl_float(167116800.0,precision)*Power(r,cl_float(14.0,precision))*Power(xi,cl_float(14.0,precision)) - 

         cl_float(13369344.0,precision)*Power(r,cl_float(15.0,precision))*Power(xi,cl_float(15.0,precision)) - cl_float(835584.0,precision)*Power(r,cl_float(16.0,precision))*Power(xi,cl_float(16.0,precision)) - 

         cl_float(32768.0,precision)*Power(r,cl_float(17.0,precision))*Power(xi,cl_float(17.0,precision)))/

       (cl_float(1.422749712384e16,precision)*exp(TWO*r*xi)*Power(r,TWO)) + 

      (xi*(-cl_float(14227497123840000.0,precision) + cl_float(14227497123840000.0,precision)*exp(TWO*r*xi) - 

           cl_float(25913502934444125.0,precision)*r*xi - cl_float(23372011621208250.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

           cl_float(13907709869303250.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

           cl_float(6137735659555500.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - 

           cl_float(2140857388870200.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

           cl_float(614116575072000.0,precision)*Power(r,SIX)*Power(xi,SIX) - 

           cl_float(148809580920000.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

           cl_float(31036639233600.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - 

           cl_float(5645342102400.0,precision)*Power(r,NINE)*Power(xi,NINE) - 

           cl_float(903333150720.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - 

           cl_float(127744081920.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)) - 

           cl_float(15968010240.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(12.0,precision)) - 

           cl_float(1754726400.0,precision)*Power(r,cl_float(13.0,precision))*Power(xi,cl_float(13.0,precision)) - 

           cl_float(167116800.0,precision)*Power(r,cl_float(14.0,precision))*Power(xi,cl_float(14.0,precision)) - 

           cl_float(13369344.0,precision)*Power(r,cl_float(15.0,precision))*Power(xi,cl_float(15.0,precision)) - cl_float(835584.0,precision)*Power(r,cl_float(16.0,precision))*Power(xi,cl_float(16.0,precision)) - 

           cl_float(32768.0,precision)*Power(r,cl_float(17.0,precision))*Power(xi,cl_float(17.0,precision))))/(cl_float(7.11374856192e15,precision)*exp(TWO*r*xi)*r)

    ; }
 
  }
  else {
      if (r == ZERO) {  S = ZERO

    ; } else {  S = (cl_float(56700.0,precision)*exp(TWO*r*(xi + xj))*Power(Power(xi,TWO) - Power(xj,TWO),cl_float(17.0,precision)) + 

         NINE*exp(TWO*r*xj)*Power(xj,cl_float(12.0,precision))*

          (-cl_float(980.0,precision)*Power(r,SIX)*Power(xi,cl_float(28.0,precision)) - cl_float(20.0,precision)*Power(r,SEVEN)*Power(xi,cl_float(29.0,precision)) + 

            cl_float(6300.0,precision)*Power(xj,cl_float(22.0,precision)) + cl_float(11025.0,precision)*r*xi*Power(xj,cl_float(22.0,precision)) - 

            cl_float(50.0,precision)*Power(r,FIVE)*Power(xi,cl_float(27.0,precision))*(cl_float(441.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(3150.0,precision)*Power(xi,TWO)*Power(xj,cl_float(20.0,precision))*(-cl_float(34.0,precision) + THREE*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(525.0,precision)*r*Power(xi,THREE)*Power(xj,cl_float(20.0,precision))*

             (-cl_float(357.0,precision) + cl_float(10.0,precision)*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(420.0,precision)*Power(r,FOUR)*Power(xi,cl_float(26.0,precision))*(cl_float(700.0,precision) + cl_float(19.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(1050.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(18.0,precision))*

             (cl_float(816.0,precision) - cl_float(153.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(210.0,precision)*r*Power(xi,FIVE)*Power(xj,cl_float(18.0,precision))*

             (cl_float(7140.0,precision) - cl_float(425.0,precision)*Power(r,TWO)*Power(xj,TWO) + THREE*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(42.0,precision)*Power(r,THREE)*Power(xi,cl_float(25.0,precision))*

             (-cl_float(59500.0,precision) - cl_float(6035.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(18.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(84.0,precision)*Power(r,TWO)*Power(xi,cl_float(24.0,precision))*

             (-cl_float(160650.0,precision) - cl_float(52700.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(397.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

            cl_float(28.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,cl_float(10.0,precision))*

             (cl_float(100849950.0,precision) + cl_float(27100125.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(186150.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(2177.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(140.0,precision)*Power(xi,SIX)*Power(xj,cl_float(16.0,precision))*

             (-cl_float(30600.0,precision) + cl_float(9180.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(255.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(2380.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(6300.0,precision) + cl_float(2700.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(120.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(10.0,precision)*r*Power(xi,SEVEN)*Power(xj,cl_float(16.0,precision))*

             (-cl_float(749700.0,precision) + cl_float(71400.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1071.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + TWO*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(204.0,precision)*r*Power(xi,cl_float(15.0,precision))*Power(xj,EIGHT)*

             (cl_float(28962255.0,precision) - cl_float(1744750.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(9555.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + SIX*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(42.0,precision)*r*Power(xi,cl_float(11.0,precision))*Power(xj,cl_float(12.0,precision))*

             (-cl_float(12911925.0,precision) - cl_float(1634550.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(7103.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(18.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            TWO*r*Power(xi,NINE)*Power(xj,cl_float(14.0,precision))*

             (cl_float(16948575.0,precision) - cl_float(1184400.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(63861.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(50.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(28.0,precision)*Power(xi,cl_float(22.0,precision))*(-cl_float(2180250.0,precision) - cl_float(10993050.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(14925.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(73.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(952.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,EIGHT)*

             (cl_float(16966215.0,precision) + cl_float(725175.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(36075.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(79.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(84.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(12.0,precision))*

             (cl_float(1723800.0,precision) + cl_float(279225.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(45600.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(107.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(35.0,precision)*r*Power(xi,cl_float(17.0,precision))*Power(xj,SIX)*

             (cl_float(132637869.0,precision) - cl_float(2205240.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(48348.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(136.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            SIX*r*Power(xi,cl_float(21.0,precision))*Power(xj,TWO)*

             (cl_float(192298050.0,precision) + cl_float(12644275.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(218029.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(204.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            FOUR*r*Power(xi,cl_float(13.0,precision))*Power(xj,cl_float(10.0,precision))*

             (cl_float(1259522775.0,precision) + cl_float(15895425.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(493017.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(263.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(140.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,SIX)*

             (cl_float(180826281.0,precision) - cl_float(15101406.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(160140.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(442.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            TWO*r*Power(xi,cl_float(23.0,precision))*(cl_float(21366450.0,precision) + cl_float(23526300.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(246729.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(526.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            SEVEN*r*Power(xi,cl_float(19.0,precision))*Power(xj,FOUR)*

             (-cl_float(811081215.0,precision) + cl_float(39095550.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(515916.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(680.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(70.0,precision)*Power(xi,cl_float(18.0,precision))*Power(xj,FOUR)*

             (-cl_float(180554454.0,precision) + cl_float(9873711.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(414120.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(2924.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(14.0,precision)*Power(xi,cl_float(20.0,precision))*Power(xj,TWO)*

             (cl_float(136919700.0,precision) + cl_float(71867115.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(2154150.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(10268.0,precision)*Power(r,SIX)*Power(xj,SIX))) 

    - FOUR*exp(TWO*r*xi)*Power(xi,cl_float(10.0,precision))*

          (-cl_float(10710.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,cl_float(12.0,precision))*

             (-cl_float(3555.0,precision) - cl_float(127008.0,precision)*r*xj + cl_float(138384.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(74556.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(22284.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(408.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(576.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(60.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + TWO*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            TWO*Power(xi,cl_float(20.0,precision))*Power(xj,FOUR)*

             (cl_float(963900.0,precision) + cl_float(1735020.0,precision)*r*xj + cl_float(1542240.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(899640.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(385560.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(128520.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(34272.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(9126.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(333.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) - 

               cl_float(20.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

            TWO*Power(xj,cl_float(24.0,precision))*(cl_float(119041650.0,precision) + cl_float(107137485.0,precision)*r*xj + 

               cl_float(45110520.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(11695320.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(2063880.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(257985.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(22932.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(1404.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(54.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               Power(r,NINE)*Power(xj,NINE)) + 

            TWO*Power(xi,TWO)*Power(xj,cl_float(22.0,precision))*

             (-cl_float(3264488325.0,precision) - cl_float(2505368880.0,precision)*r*xj - 

               cl_float(881390160.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(185775660.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(25639740.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(2361555.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(139356.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(4482.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - 

               cl_float(27.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + TWO*Power(r,NINE)*Power(xj,NINE)) + 

            Power(xi,cl_float(24.0,precision))*(cl_float(14175.0,precision) + cl_float(25515.0,precision)*r*xj + cl_float(22680.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(13230.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(5670.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(1890.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(504.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(108.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               TWO*Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(102.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(14.0,precision))*

             (cl_float(44986725.0,precision) - cl_float(97433280.0,precision)*r*xj + cl_float(44467920.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(15857100.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(457380.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(620550.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(83160.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(4068.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - 

               SIX*Power(r,EIGHT)*Power(xj,EIGHT) + FOUR*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(102.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,cl_float(10.0,precision))*

             (-cl_float(859950.0,precision) - cl_float(1437345.0,precision)*r*xj - cl_float(2260440.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(810810.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(1056510.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(217854.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(6552.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(3852.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(258.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               FOUR*Power(r,NINE)*Power(xj,NINE)) - 

            Power(xi,cl_float(22.0,precision))*Power(xj,TWO)*

             (cl_float(240975.0,precision) + cl_float(433755.0,precision)*r*xj + cl_float(385560.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(224910.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(96390.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(32130.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(8568.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(1836.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(306.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               FOUR*Power(r,NINE)*Power(xj,NINE)) + 

            TWO*Power(xi,FOUR)*Power(xj,cl_float(20.0,precision))*

             (-cl_float(18032978565.0,precision) - cl_float(9823683240.0,precision)*r*xj - 

               cl_float(2047323600.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(129098340.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(26410860.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(7094304.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(788256.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(48654.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(1593.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(20.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

            SIX*Power(xi,cl_float(16.0,precision))*Power(xj,EIGHT)*

             (-cl_float(5622750.0,precision) - cl_float(10120950.0,precision)*r*xj - cl_float(8996400.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(5698350.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(897750.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(1641591.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(211932.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(10224.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(2364.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(73.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

            TWO*Power(xi,cl_float(18.0,precision))*Power(xj,SIX)*

             (-cl_float(4819500.0,precision) - cl_float(8675100.0,precision)*r*xj - cl_float(7711200.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(4498200.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(1927800.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(561519.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(279468.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

               cl_float(20682.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(1305.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               cl_float(106.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

            THREE*Power(xi,EIGHT)*Power(xj,cl_float(16.0,precision))*

             (-cl_float(9364244085.0,precision) + cl_float(6940428705.0,precision)*r*xj + 

               cl_float(2117684520.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(230268150.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(149610510.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(21824334.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(1223208.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(12708.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(4470.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(146.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

            Power(xi,SIX)*Power(xj,cl_float(18.0,precision))*

             (cl_float(57304872765.0,precision) + cl_float(7147185255.0,precision)*r*xj - 

               cl_float(5801702760.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(2053388610.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(271655370.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(10864854.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(1337112.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(202716.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(10746.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(212.0,precision)*Power(r,NINE)*Power(xj,NINE))))/

       (cl_float(56700.0,precision)*exp(TWO*r*(xi + xj))*Power(r,TWO)*Power(xi - xj,cl_float(17.0,precision))*

         Power(xi + xj,cl_float(17.0,precision))) + (cl_float(56700.0,precision)*exp(TWO*r*(xi + xj))*

          Power(Power(xi,TWO) - Power(xj,TWO),cl_float(17.0,precision)) + 

         NINE*exp(TWO*r*xj)*Power(xj,cl_float(12.0,precision))*

          (-cl_float(980.0,precision)*Power(r,SIX)*Power(xi,cl_float(28.0,precision)) - cl_float(20.0,precision)*Power(r,SEVEN)*Power(xi,cl_float(29.0,precision)) + 

            cl_float(6300.0,precision)*Power(xj,cl_float(22.0,precision)) + cl_float(11025.0,precision)*r*xi*Power(xj,cl_float(22.0,precision)) - 

            cl_float(50.0,precision)*Power(r,FIVE)*Power(xi,cl_float(27.0,precision))*(cl_float(441.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(3150.0,precision)*Power(xi,TWO)*Power(xj,cl_float(20.0,precision))*(-cl_float(34.0,precision) + THREE*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(525.0,precision)*r*Power(xi,THREE)*Power(xj,cl_float(20.0,precision))*

             (-cl_float(357.0,precision) + cl_float(10.0,precision)*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(420.0,precision)*Power(r,FOUR)*Power(xi,cl_float(26.0,precision))*(cl_float(700.0,precision) + cl_float(19.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(1050.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(18.0,precision))*

             (cl_float(816.0,precision) - cl_float(153.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(210.0,precision)*r*Power(xi,FIVE)*Power(xj,cl_float(18.0,precision))*

             (cl_float(7140.0,precision) - cl_float(425.0,precision)*Power(r,TWO)*Power(xj,TWO) + THREE*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(42.0,precision)*Power(r,THREE)*Power(xi,cl_float(25.0,precision))*

             (-cl_float(59500.0,precision) - cl_float(6035.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(18.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(84.0,precision)*Power(r,TWO)*Power(xi,cl_float(24.0,precision))*

             (-cl_float(160650.0,precision) - cl_float(52700.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(397.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

            cl_float(28.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,cl_float(10.0,precision))*

             (cl_float(100849950.0,precision) + cl_float(27100125.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(186150.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(2177.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(140.0,precision)*Power(xi,SIX)*Power(xj,cl_float(16.0,precision))*

             (-cl_float(30600.0,precision) + cl_float(9180.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(255.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(2380.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(6300.0,precision) + cl_float(2700.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(120.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(10.0,precision)*r*Power(xi,SEVEN)*Power(xj,cl_float(16.0,precision))*

             (-cl_float(749700.0,precision) + cl_float(71400.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1071.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + TWO*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(204.0,precision)*r*Power(xi,cl_float(15.0,precision))*Power(xj,EIGHT)*

             (cl_float(28962255.0,precision) - cl_float(1744750.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(9555.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + SIX*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(42.0,precision)*r*Power(xi,cl_float(11.0,precision))*Power(xj,cl_float(12.0,precision))*

             (-cl_float(12911925.0,precision) - cl_float(1634550.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(7103.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(18.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            TWO*r*Power(xi,NINE)*Power(xj,cl_float(14.0,precision))*

             (cl_float(16948575.0,precision) - cl_float(1184400.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(63861.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(50.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(28.0,precision)*Power(xi,cl_float(22.0,precision))*(-cl_float(2180250.0,precision) - cl_float(10993050.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(14925.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(73.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(952.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,EIGHT)*

             (cl_float(16966215.0,precision) + cl_float(725175.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(36075.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(79.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(84.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(12.0,precision))*

             (cl_float(1723800.0,precision) + cl_float(279225.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(45600.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(107.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(35.0,precision)*r*Power(xi,cl_float(17.0,precision))*Power(xj,SIX)*

             (cl_float(132637869.0,precision) - cl_float(2205240.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(48348.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(136.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            SIX*r*Power(xi,cl_float(21.0,precision))*Power(xj,TWO)*

             (cl_float(192298050.0,precision) + cl_float(12644275.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(218029.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(204.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            FOUR*r*Power(xi,cl_float(13.0,precision))*Power(xj,cl_float(10.0,precision))*

             (cl_float(1259522775.0,precision) + cl_float(15895425.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(493017.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(263.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(140.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,SIX)*

             (cl_float(180826281.0,precision) - cl_float(15101406.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(160140.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(442.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            TWO*r*Power(xi,cl_float(23.0,precision))*(cl_float(21366450.0,precision) + cl_float(23526300.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(246729.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(526.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            SEVEN*r*Power(xi,cl_float(19.0,precision))*Power(xj,FOUR)*

             (-cl_float(811081215.0,precision) + cl_float(39095550.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(515916.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(680.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(70.0,precision)*Power(xi,cl_float(18.0,precision))*Power(xj,FOUR)*

             (-cl_float(180554454.0,precision) + cl_float(9873711.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(414120.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(2924.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(14.0,precision)*Power(xi,cl_float(20.0,precision))*Power(xj,TWO)*

             (cl_float(136919700.0,precision) + cl_float(71867115.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(2154150.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(10268.0,precision)*Power(r,SIX)*Power(xj,SIX))) 

    - FOUR*exp(TWO*r*xi)*Power(xi,cl_float(10.0,precision))*

          (-cl_float(10710.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,cl_float(12.0,precision))*

             (-cl_float(3555.0,precision) - cl_float(127008.0,precision)*r*xj + cl_float(138384.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(74556.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(22284.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(408.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(576.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(60.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + TWO*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            TWO*Power(xi,cl_float(20.0,precision))*Power(xj,FOUR)*

             (cl_float(963900.0,precision) + cl_float(1735020.0,precision)*r*xj + cl_float(1542240.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(899640.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(385560.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(128520.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(34272.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(9126.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(333.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) - 

               cl_float(20.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

            TWO*Power(xj,cl_float(24.0,precision))*(cl_float(119041650.0,precision) + cl_float(107137485.0,precision)*r*xj + 

               cl_float(45110520.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(11695320.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(2063880.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(257985.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(22932.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(1404.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(54.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               Power(r,NINE)*Power(xj,NINE)) + 

            TWO*Power(xi,TWO)*Power(xj,cl_float(22.0,precision))*

             (-cl_float(3264488325.0,precision) - cl_float(2505368880.0,precision)*r*xj - 

               cl_float(881390160.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(185775660.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(25639740.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(2361555.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(139356.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(4482.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - 

               cl_float(27.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + TWO*Power(r,NINE)*Power(xj,NINE)) + 

            Power(xi,cl_float(24.0,precision))*(cl_float(14175.0,precision) + cl_float(25515.0,precision)*r*xj + cl_float(22680.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(13230.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(5670.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(1890.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(504.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(108.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               TWO*Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(102.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(14.0,precision))*

             (cl_float(44986725.0,precision) - cl_float(97433280.0,precision)*r*xj + cl_float(44467920.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(15857100.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(457380.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(620550.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(83160.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(4068.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - 

               SIX*Power(r,EIGHT)*Power(xj,EIGHT) + FOUR*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(102.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,cl_float(10.0,precision))*

             (-cl_float(859950.0,precision) - cl_float(1437345.0,precision)*r*xj - cl_float(2260440.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(810810.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(1056510.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(217854.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(6552.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(3852.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(258.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               FOUR*Power(r,NINE)*Power(xj,NINE)) - 

            Power(xi,cl_float(22.0,precision))*Power(xj,TWO)*

             (cl_float(240975.0,precision) + cl_float(433755.0,precision)*r*xj + cl_float(385560.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(224910.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(96390.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(32130.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(8568.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(1836.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(306.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               FOUR*Power(r,NINE)*Power(xj,NINE)) + 

            TWO*Power(xi,FOUR)*Power(xj,cl_float(20.0,precision))*

             (-cl_float(18032978565.0,precision) - cl_float(9823683240.0,precision)*r*xj - 

               cl_float(2047323600.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(129098340.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(26410860.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(7094304.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(788256.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(48654.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(1593.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(20.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

            SIX*Power(xi,cl_float(16.0,precision))*Power(xj,EIGHT)*

             (-cl_float(5622750.0,precision) - cl_float(10120950.0,precision)*r*xj - cl_float(8996400.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(5698350.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(897750.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(1641591.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(211932.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(10224.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(2364.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(73.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

            TWO*Power(xi,cl_float(18.0,precision))*Power(xj,SIX)*

             (-cl_float(4819500.0,precision) - cl_float(8675100.0,precision)*r*xj - cl_float(7711200.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(4498200.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(1927800.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(561519.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(279468.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

               cl_float(20682.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(1305.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               cl_float(106.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

            THREE*Power(xi,EIGHT)*Power(xj,cl_float(16.0,precision))*

             (-cl_float(9364244085.0,precision) + cl_float(6940428705.0,precision)*r*xj + 

               cl_float(2117684520.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(230268150.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(149610510.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(21824334.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(1223208.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(12708.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(4470.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(146.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

            Power(xi,SIX)*Power(xj,cl_float(18.0,precision))*

             (cl_float(57304872765.0,precision) + cl_float(7147185255.0,precision)*r*xj - 

               cl_float(5801702760.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(2053388610.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(271655370.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(10864854.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(1337112.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(202716.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(10746.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(212.0,precision)*Power(r,NINE)*Power(xj,NINE))))/

       (cl_float(28350.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,cl_float(17.0,precision))*Power(xi + xj,cl_float(16.0,precision))) - 

      (cl_float(113400.0,precision)*exp(TWO*r*(xi + xj))*(xi + xj)*

          Power(Power(xi,TWO) - Power(xj,TWO),cl_float(17.0,precision)) + 

         NINE*exp(TWO*r*xj)*Power(xj,cl_float(12.0,precision))*

          (-cl_float(5880.0,precision)*Power(r,FIVE)*Power(xi,cl_float(28.0,precision)) - cl_float(140.0,precision)*Power(r,SIX)*Power(xi,cl_float(29.0,precision)) - 

            cl_float(15960.0,precision)*Power(r,FIVE)*Power(xi,cl_float(26.0,precision))*Power(xj,TWO) - 

            cl_float(200.0,precision)*Power(r,SIX)*Power(xi,cl_float(27.0,precision))*Power(xj,TWO) + cl_float(11025.0,precision)*xi*Power(xj,cl_float(22.0,precision)) + 

            cl_float(18900.0,precision)*r*Power(xi,TWO)*Power(xj,cl_float(22.0,precision)) + 

            cl_float(10500.0,precision)*Power(r,TWO)*Power(xi,THREE)*Power(xj,cl_float(22.0,precision)) - 

            cl_float(250.0,precision)*Power(r,FOUR)*Power(xi,cl_float(27.0,precision))*(cl_float(441.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(525.0,precision)*Power(xi,THREE)*Power(xj,cl_float(20.0,precision))*(-cl_float(357.0,precision) + cl_float(10.0,precision)*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(1680.0,precision)*Power(r,THREE)*Power(xi,cl_float(26.0,precision))*(cl_float(700.0,precision) + cl_float(19.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(1050.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(18.0,precision))*

             (-cl_float(306.0,precision)*r*Power(xj,TWO) + EIGHT*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(210.0,precision)*r*Power(xi,FIVE)*Power(xj,cl_float(18.0,precision))*

             (-cl_float(850.0,precision)*r*Power(xj,TWO) + cl_float(12.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(42.0,precision)*Power(r,THREE)*Power(xi,cl_float(25.0,precision))*

             (-cl_float(12070.0,precision)*r*Power(xj,TWO) + cl_float(72.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(84.0,precision)*Power(r,TWO)*Power(xi,cl_float(24.0,precision))*

             (-cl_float(105400.0,precision)*r*Power(xj,TWO) + cl_float(1588.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(210.0,precision)*Power(xi,FIVE)*Power(xj,cl_float(18.0,precision))*

             (cl_float(7140.0,precision) - cl_float(425.0,precision)*Power(r,TWO)*Power(xj,TWO) + THREE*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(126.0,precision)*Power(r,TWO)*Power(xi,cl_float(25.0,precision))*

             (-cl_float(59500.0,precision) - cl_float(6035.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(18.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) 

    + cl_float(168.0,precision)*r*Power(xi,cl_float(24.0,precision))*(-cl_float(160650.0,precision) - cl_float(52700.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(397.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

            cl_float(28.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,cl_float(10.0,precision))*

             (cl_float(54200250.0,precision)*r*Power(xj,TWO) + cl_float(744600.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(13062.0,precision)*Power(r,FIVE)*Power(xj,SIX)) + 

            cl_float(140.0,precision)*Power(xi,SIX)*Power(xj,cl_float(16.0,precision))*

             (cl_float(18360.0,precision)*r*Power(xj,TWO) - cl_float(1020.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               SIX*Power(r,FIVE)*Power(xj,SIX)) - 

            cl_float(2380.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(14.0,precision))*

             (cl_float(5400.0,precision)*r*Power(xj,TWO) - cl_float(480.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               SIX*Power(r,FIVE)*Power(xj,SIX)) + 

            cl_float(10.0,precision)*r*Power(xi,SEVEN)*Power(xj,cl_float(16.0,precision))*

             (cl_float(142800.0,precision)*r*Power(xj,TWO) - cl_float(4284.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(12.0,precision)*Power(r,FIVE)*Power(xj,SIX)) + 

            cl_float(204.0,precision)*r*Power(xi,cl_float(15.0,precision))*Power(xj,EIGHT)*

             (-cl_float(3489500.0,precision)*r*Power(xj,TWO) + cl_float(38220.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(36.0,precision)*Power(r,FIVE)*Power(xj,SIX)) - 

            cl_float(42.0,precision)*r*Power(xi,cl_float(11.0,precision))*Power(xj,cl_float(12.0,precision))*

             (-cl_float(3269100.0,precision)*r*Power(xj,TWO) - cl_float(28412.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(108.0,precision)*Power(r,FIVE)*Power(xj,SIX)) + 

            TWO*r*Power(xi,NINE)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(2368800.0,precision)*r*Power(xj,TWO) + cl_float(255444.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(300.0,precision)*Power(r,FIVE)*Power(xj,SIX)) + 

            cl_float(28.0,precision)*Power(xi,cl_float(22.0,precision))*(-cl_float(21986100.0,precision)*r*Power(xj,TWO) + 

               cl_float(59700.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(438.0,precision)*Power(r,FIVE)*Power(xj,SIX)) - 

            cl_float(952.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,EIGHT)*

             (cl_float(1450350.0,precision)*r*Power(xj,TWO) - cl_float(144300.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(474.0,precision)*Power(r,FIVE)*Power(xj,SIX)) - 

            cl_float(84.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(12.0,precision))*

             (cl_float(558450.0,precision)*r*Power(xj,TWO) + cl_float(182400.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(642.0,precision)*Power(r,FIVE)*Power(xj,SIX)) - 

            cl_float(35.0,precision)*r*Power(xi,cl_float(17.0,precision))*Power(xj,SIX)*

             (-cl_float(4410480.0,precision)*r*Power(xj,TWO) - cl_float(193392.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(816.0,precision)*Power(r,FIVE)*Power(xj,SIX)) - 

            SIX*r*Power(xi,cl_float(21.0,precision))*Power(xj,TWO)*

             (cl_float(25288550.0,precision)*r*Power(xj,TWO) - cl_float(872116.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(1224.0,precision)*Power(r,FIVE)*Power(xj,SIX)) + 

            FOUR*r*Power(xi,cl_float(13.0,precision))*Power(xj,cl_float(10.0,precision))*

             (cl_float(31790850.0,precision)*r*Power(xj,TWO) - cl_float(1972068.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(1578.0,precision)*Power(r,FIVE)*Power(xj,SIX)) - 

            cl_float(140.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,SIX)*

             (-cl_float(30202812.0,precision)*r*Power(xj,TWO) + cl_float(640560.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(2652.0,precision)*Power(r,FIVE)*Power(xj,SIX)) - 

            TWO*r*Power(xi,cl_float(23.0,precision))*(cl_float(47052600.0,precision)*r*Power(xj,TWO) - 

               cl_float(986916.0,precision)*Power(r,THREE)*Power(xj,FOUR) + cl_float(3156.0,precision)*Power(r,FIVE)*Power(xj,SIX)) + 

            SEVEN*r*Power(xi,cl_float(19.0,precision))*Power(xj,FOUR)*

             (cl_float(78191100.0,precision)*r*Power(xj,TWO) - cl_float(2063664.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(4080.0,precision)*Power(r,FIVE)*Power(xj,SIX)) + 

            cl_float(70.0,precision)*Power(xi,cl_float(18.0,precision))*Power(xj,FOUR)*

             (cl_float(19747422.0,precision)*r*Power(xj,TWO) - cl_float(1656480.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(17544.0,precision)*Power(r,FIVE)*Power(xj,SIX)) - 

            cl_float(14.0,precision)*Power(xi,cl_float(20.0,precision))*Power(xj,TWO)*

             (cl_float(143734230.0,precision)*r*Power(xj,TWO) - cl_float(8616600.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(61608.0,precision)*Power(r,FIVE)*Power(xj,SIX)) + 

            cl_float(10.0,precision)*Power(xi,SEVEN)*Power(xj,cl_float(16.0,precision))*

             (-cl_float(749700.0,precision) + cl_float(71400.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1071.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + TWO*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(204.0,precision)*Power(xi,cl_float(15.0,precision))*Power(xj,EIGHT)*

             (cl_float(28962255.0,precision) - cl_float(1744750.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(9555.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + SIX*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(42.0,precision)*Power(xi,cl_float(11.0,precision))*Power(xj,cl_float(12.0,precision))*

             (-cl_float(12911925.0,precision) - cl_float(1634550.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(7103.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(18.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            TWO*Power(xi,NINE)*Power(xj,cl_float(14.0,precision))*

             (cl_float(16948575.0,precision) - cl_float(1184400.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(63861.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(50.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(35.0,precision)*Power(xi,cl_float(17.0,precision))*Power(xj,SIX)*

             (cl_float(132637869.0,precision) - cl_float(2205240.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(48348.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(136.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            SIX*Power(xi,cl_float(21.0,precision))*Power(xj,TWO)*

             (cl_float(192298050.0,precision) + cl_float(12644275.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(218029.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(204.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            FOUR*Power(xi,cl_float(13.0,precision))*Power(xj,cl_float(10.0,precision))*

             (cl_float(1259522775.0,precision) + cl_float(15895425.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(493017.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(263.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            TWO*Power(xi,cl_float(23.0,precision))*(cl_float(21366450.0,precision) + cl_float(23526300.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(246729.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(526.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            SEVEN*Power(xi,cl_float(19.0,precision))*Power(xj,FOUR)*

             (-cl_float(811081215.0,precision) + cl_float(39095550.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(515916.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(680.0,precision)*Power(r,SIX)*Power(xj,SIX))) + 

         cl_float(18.0,precision)*exp(TWO*r*xj)*Power(xj,cl_float(13.0,precision))*

          (-cl_float(980.0,precision)*Power(r,SIX)*Power(xi,cl_float(28.0,precision)) - cl_float(20.0,precision)*Power(r,SEVEN)*Power(xi,cl_float(29.0,precision)) + 

            cl_float(6300.0,precision)*Power(xj,cl_float(22.0,precision)) + cl_float(11025.0,precision)*r*xi*Power(xj,cl_float(22.0,precision)) - 

            cl_float(50.0,precision)*Power(r,FIVE)*Power(xi,cl_float(27.0,precision))*(cl_float(441.0,precision) + TWO*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(3150.0,precision)*Power(xi,TWO)*Power(xj,cl_float(20.0,precision))*(-cl_float(34.0,precision) + THREE*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(525.0,precision)*r*Power(xi,THREE)*Power(xj,cl_float(20.0,precision))*(-cl_float(357.0,precision) + cl_float(10.0,precision)*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(420.0,precision)*Power(r,FOUR)*Power(xi,cl_float(26.0,precision))*(cl_float(700.0,precision) + cl_float(19.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(1050.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(18.0,precision))*

             (cl_float(816.0,precision) - cl_float(153.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(210.0,precision)*r*Power(xi,FIVE)*Power(xj,cl_float(18.0,precision))*

             (cl_float(7140.0,precision) - cl_float(425.0,precision)*Power(r,TWO)*Power(xj,TWO) + THREE*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(42.0,precision)*Power(r,THREE)*Power(xi,cl_float(25.0,precision))*

             (-cl_float(59500.0,precision) - cl_float(6035.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(18.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) 

    + cl_float(84.0,precision)*Power(r,TWO)*Power(xi,cl_float(24.0,precision))*(-cl_float(160650.0,precision) - cl_float(52700.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(397.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) - 

            cl_float(28.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,cl_float(10.0,precision))*

             (cl_float(100849950.0,precision) + cl_float(27100125.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(186150.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(2177.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(140.0,precision)*Power(xi,SIX)*Power(xj,cl_float(16.0,precision))*

             (-cl_float(30600.0,precision) + cl_float(9180.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(255.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(2380.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(14.0,precision))*

             (-cl_float(6300.0,precision) + cl_float(2700.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(120.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(10.0,precision)*r*Power(xi,SEVEN)*Power(xj,cl_float(16.0,precision))*

             (-cl_float(749700.0,precision) + cl_float(71400.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1071.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + TWO*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(204.0,precision)*r*Power(xi,cl_float(15.0,precision))*Power(xj,EIGHT)*

             (cl_float(28962255.0,precision) - cl_float(1744750.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(9555.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + SIX*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(42.0,precision)*r*Power(xi,cl_float(11.0,precision))*Power(xj,cl_float(12.0,precision))*

             (-cl_float(12911925.0,precision) - cl_float(1634550.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(7103.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(18.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            TWO*r*Power(xi,NINE)*Power(xj,cl_float(14.0,precision))*

             (cl_float(16948575.0,precision) - cl_float(1184400.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(63861.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(50.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(28.0,precision)*Power(xi,cl_float(22.0,precision))*(-cl_float(2180250.0,precision) - cl_float(10993050.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(14925.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(73.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(952.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,EIGHT)*

             (cl_float(16966215.0,precision) + cl_float(725175.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(36075.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(79.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(84.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(12.0,precision))*

             (cl_float(1723800.0,precision) + cl_float(279225.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(45600.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(107.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(35.0,precision)*r*Power(xi,cl_float(17.0,precision))*Power(xj,SIX)*

             (cl_float(132637869.0,precision) - cl_float(2205240.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(48348.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(136.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            SIX*r*Power(xi,cl_float(21.0,precision))*Power(xj,TWO)*

             (cl_float(192298050.0,precision) + cl_float(12644275.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(218029.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(204.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            FOUR*r*Power(xi,cl_float(13.0,precision))*Power(xj,cl_float(10.0,precision))*

             (cl_float(1259522775.0,precision) + cl_float(15895425.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(493017.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(263.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(140.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,SIX)*

             (cl_float(180826281.0,precision) - cl_float(15101406.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(160140.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(442.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            TWO*r*Power(xi,cl_float(23.0,precision))*(cl_float(21366450.0,precision) + cl_float(23526300.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(246729.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(526.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            SEVEN*r*Power(xi,cl_float(19.0,precision))*Power(xj,FOUR)*

             (-cl_float(811081215.0,precision) + cl_float(39095550.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(515916.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(680.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(70.0,precision)*Power(xi,cl_float(18.0,precision))*Power(xj,FOUR)*

             (-cl_float(180554454.0,precision) + cl_float(9873711.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(414120.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(2924.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(14.0,precision)*Power(xi,cl_float(20.0,precision))*Power(xj,TWO)*

             (cl_float(136919700.0,precision) + cl_float(71867115.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(2154150.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(10268.0,precision)*Power(r,SIX)*Power(xj,SIX))) - 

         FOUR*exp(TWO*r*xi)*Power(xi,cl_float(10.0,precision))*

          (-cl_float(10710.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,cl_float(12.0,precision))*

             (-cl_float(127008.0,precision)*xj + cl_float(276768.0,precision)*r*Power(xj,TWO) - 

               cl_float(223668.0,precision)*Power(r,TWO)*Power(xj,THREE) - cl_float(89136.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(2040.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(3456.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(420.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(16.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT)) + 

            TWO*Power(xi,cl_float(20.0,precision))*Power(xj,FOUR)*

             (cl_float(1735020.0,precision)*xj + cl_float(3084480.0,precision)*r*Power(xj,TWO) + 

               cl_float(2698920.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(1542240.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(642600.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(205632.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(63882.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(2664.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) - 

               cl_float(180.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) - 

            TWO*Power(xj,cl_float(24.0,precision))*(cl_float(107137485.0,precision)*xj + cl_float(90221040.0,precision)*r*Power(xj,TWO) + 

               cl_float(35085960.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(8255520.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(1289925.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(137592.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(9828.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + 

               cl_float(432.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + NINE*Power(r,EIGHT)*Power(xj,NINE)) + 

            TWO*Power(xi,TWO)*Power(xj,cl_float(22.0,precision))*

             (-cl_float(2505368880.0,precision)*xj - cl_float(1762780320.0,precision)*r*Power(xj,TWO) - 

               cl_float(557326980.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(102558960.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(11807775.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - 

               cl_float(836136.0,precision)*Power(r,FIVE)*Power(xj,SIX) - cl_float(31374.0,precision)*Power(r,SIX)*Power(xj,SEVEN) - 

               cl_float(216.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            Power(xi,cl_float(24.0,precision))*(cl_float(25515.0,precision)*xj + cl_float(45360.0,precision)*r*Power(xj,TWO) + 

               cl_float(39690.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(22680.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(9450.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(3024.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(756.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(144.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) - 

            cl_float(102.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(14.0,precision))*

             (-cl_float(97433280.0,precision)*xj + cl_float(88935840.0,precision)*r*Power(xj,TWO) + 

               cl_float(47571300.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(1829520.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(3102750.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - 

               cl_float(498960.0,precision)*Power(r,FIVE)*Power(xj,SIX) - cl_float(28476.0,precision)*Power(r,SIX)*Power(xj,SEVEN) - 

               cl_float(48.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + cl_float(36.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            cl_float(102.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,cl_float(10.0,precision))*

             (-cl_float(1437345.0,precision)*xj - cl_float(4520880.0,precision)*r*Power(xj,TWO) + 

               cl_float(2432430.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(4226040.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(1089270.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(39312.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(26964.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(2064.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(36.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) - 

            Power(xi,cl_float(22.0,precision))*Power(xj,TWO)*

             (cl_float(433755.0,precision)*xj + cl_float(771120.0,precision)*r*Power(xj,TWO) + 

               cl_float(674730.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(385560.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(160650.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(51408.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(12852.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(2448.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(36.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            TWO*Power(xi,FOUR)*Power(xj,cl_float(20.0,precision))*

             (-cl_float(9823683240.0,precision)*xj - cl_float(4094647200.0,precision)*r*Power(xj,TWO) - 

               cl_float(387295020.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(105643440.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(35471520.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(4729536.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(340578.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(12744.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(180.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) - 

            SIX*Power(xi,cl_float(16.0,precision))*Power(xj,EIGHT)*

             (-cl_float(10120950.0,precision)*xj - cl_float(17992800.0,precision)*r*Power(xj,TWO) - 

               cl_float(17095050.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(3591000.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(8207955.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - 

               cl_float(1271592.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(71568.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + 

               cl_float(18912.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + cl_float(657.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            TWO*Power(xi,cl_float(18.0,precision))*Power(xj,SIX)*

             (-cl_float(8675100.0,precision)*xj - cl_float(15422400.0,precision)*r*Power(xj,TWO) - 

               cl_float(13494600.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(7711200.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(2807595.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - 

               cl_float(1676808.0,precision)*Power(r,FIVE)*Power(xj,SIX) - 

               cl_float(144774.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(10440.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(954.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            THREE*Power(xi,EIGHT)*Power(xj,cl_float(16.0,precision))*

             (cl_float(6940428705.0,precision)*xj + cl_float(4235369040.0,precision)*r*Power(xj,TWO) - 

               cl_float(690804450.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(598442040.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(109121670.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - 

               cl_float(7339248.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(88956.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + 

               cl_float(35760.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + cl_float(1314.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) - 

            Power(xi,SIX)*Power(xj,cl_float(18.0,precision))*

             (cl_float(7147185255.0,precision)*xj - cl_float(11603405520.0,precision)*r*Power(xj,TWO) - 

               cl_float(6160165830.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(1086621480.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(54324270.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(8022672.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(1419012.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(85968.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(1908.0,precision)*Power(r,EIGHT)*Power(xj,NINE))) - 

         EIGHT*exp(TWO*r*xi)*Power(xi,cl_float(11.0,precision))*

          (-cl_float(10710.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,cl_float(12.0,precision))*

             (-cl_float(3555.0,precision) - cl_float(127008.0,precision)*r*xj + cl_float(138384.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(74556.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(22284.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(408.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(576.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(60.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + TWO*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            TWO*Power(xi,cl_float(20.0,precision))*Power(xj,FOUR)*

             (cl_float(963900.0,precision) + cl_float(1735020.0,precision)*r*xj + cl_float(1542240.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(899640.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(385560.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(128520.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(34272.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(9126.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(333.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) - 

               cl_float(20.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

            TWO*Power(xj,cl_float(24.0,precision))*(cl_float(119041650.0,precision) + cl_float(107137485.0,precision)*r*xj + 

               cl_float(45110520.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(11695320.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(2063880.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(257985.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(22932.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(1404.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(54.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + Power(r,NINE)*Power(xj,NINE)) + 

            TWO*Power(xi,TWO)*Power(xj,cl_float(22.0,precision))*

             (-cl_float(3264488325.0,precision) - cl_float(2505368880.0,precision)*r*xj - 

               cl_float(881390160.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(185775660.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(25639740.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(2361555.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(139356.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

               cl_float(4482.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - cl_float(27.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               TWO*Power(r,NINE)*Power(xj,NINE)) + 

            Power(xi,cl_float(24.0,precision))*(cl_float(14175.0,precision) + cl_float(25515.0,precision)*r*xj + cl_float(22680.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(13230.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(5670.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(1890.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(504.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(108.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               TWO*Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(102.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(14.0,precision))*

             (cl_float(44986725.0,precision) - cl_float(97433280.0,precision)*r*xj + cl_float(44467920.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(15857100.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(457380.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(620550.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(83160.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(4068.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - 

               SIX*Power(r,EIGHT)*Power(xj,EIGHT) + FOUR*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(102.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,cl_float(10.0,precision))*

             (-cl_float(859950.0,precision) - cl_float(1437345.0,precision)*r*xj - cl_float(2260440.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(810810.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(1056510.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(217854.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(6552.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(3852.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(258.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               FOUR*Power(r,NINE)*Power(xj,NINE)) - 

            Power(xi,cl_float(22.0,precision))*Power(xj,TWO)*

             (cl_float(240975.0,precision) + cl_float(433755.0,precision)*r*xj + cl_float(385560.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(224910.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(96390.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(32130.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(8568.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(1836.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(306.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               FOUR*Power(r,NINE)*Power(xj,NINE)) + 

            TWO*Power(xi,FOUR)*Power(xj,cl_float(20.0,precision))*

             (-cl_float(18032978565.0,precision) - cl_float(9823683240.0,precision)*r*xj - 

               cl_float(2047323600.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(129098340.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(26410860.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(7094304.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(788256.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(48654.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(1593.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               cl_float(20.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

            SIX*Power(xi,cl_float(16.0,precision))*Power(xj,EIGHT)*

             (-cl_float(5622750.0,precision) - cl_float(10120950.0,precision)*r*xj - cl_float(8996400.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(5698350.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(897750.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(1641591.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(211932.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(10224.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(2364.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               cl_float(73.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

            TWO*Power(xi,cl_float(18.0,precision))*Power(xj,SIX)*

             (-cl_float(4819500.0,precision) - cl_float(8675100.0,precision)*r*xj - cl_float(7711200.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(4498200.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(1927800.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(561519.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(279468.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(20682.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(1305.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(106.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

            THREE*Power(xi,EIGHT)*Power(xj,cl_float(16.0,precision))*

             (-cl_float(9364244085.0,precision) + cl_float(6940428705.0,precision)*r*xj + 

               cl_float(2117684520.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(230268150.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(149610510.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(21824334.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(1223208.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(12708.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(4470.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(146.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

            Power(xi,SIX)*Power(xj,cl_float(18.0,precision))*

             (cl_float(57304872765.0,precision) + cl_float(7147185255.0,precision)*r*xj - 

               cl_float(5801702760.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(2053388610.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(271655370.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(10864854.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(1337112.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(202716.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(10746.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(212.0,precision)*Power(r,NINE)*Power(xj,NINE))))/

       (cl_float(56700.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,cl_float(17.0,precision))*Power(xi + xj,cl_float(17.0,precision)))

    ; }
   
  }
  return S;
}

cl_F DSlater_4S_1S(cl_F r,cl_F xi,cl_F xj)
{
  return DSlater_1S_4S(r,xj,xi);
}

cl_F DSlater_4S_2S(cl_F r,cl_F xi,cl_F xj)
{
  return DSlater_2S_4S(r,xj,xi);
}

cl_F DSlater_4S_3S(cl_F r,cl_F xi,cl_F xj)
{
  return DSlater_3S_4S(r,xj,xi);
}

cl_F DSlater_5S_5S(cl_F r,cl_F xi,cl_F xj)
{
  cl_F S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == ZERO) {  S = ZERO

    ; } else {  S = -(-cl_float(22324788235240115625.0,precision)*xi + cl_float(24329020081766400000.0,precision)*exp(TWO*r*xi)*xi - 

          cl_float(40641112777427662500.0,precision)*r*Power(xi,TWO) - 

          cl_float(36677772259285657500.0,precision)*Power(r,TWO)*Power(xi,THREE) - 

          cl_float(21869785393976520000.0,precision)*Power(r,THREE)*Power(xi,FOUR) - 

          cl_float(9688099714323030000.0,precision)*Power(r,FOUR)*Power(xi,FIVE) - 

          cl_float(3399172756931952000.0,precision)*Power(r,FIVE)*Power(xi,SIX) - 

          cl_float(983239817883523200.0,precision)*Power(r,SIX)*Power(xi,SEVEN) - 

          cl_float(240924879420825600.0,precision)*Power(r,SEVEN)*Power(xi,EIGHT) - 

          cl_float(50973581199340800.0,precision)*Power(r,EIGHT)*Power(xi,NINE) - 

          cl_float(9439831425024000.0,precision)*Power(r,NINE)*Power(xi,cl_float(10.0,precision)) - 

          cl_float(1544699687731200.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(11.0,precision)) - 

          cl_float(224683590942720.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(12.0,precision)) - 

          cl_float(29125650677760.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(13.0,precision)) - 

          cl_float(3360652001280.0,precision)*Power(r,cl_float(13.0,precision))*Power(xi,cl_float(14.0,precision)) - 

          cl_float(342923673600.0,precision)*Power(r,cl_float(14.0,precision))*Power(xi,cl_float(15.0,precision)) - 

          cl_float(30482104320.0,precision)*Power(r,cl_float(15.0,precision))*Power(xi,cl_float(16.0,precision)) - 

          cl_float(2286157824.0,precision)*Power(r,cl_float(16.0,precision))*Power(xi,cl_float(17.0,precision)) - 

          cl_float(134479872.0,precision)*Power(r,cl_float(17.0,precision))*Power(xi,cl_float(18.0,precision)) - cl_float(4980736.0,precision)*Power(r,cl_float(18.0,precision))*Power(xi,cl_float(19.0,precision)))/

       (cl_float(1.21645100408832e19,precision)*exp(TWO*r*xi)*r) + 

      (-cl_float(12164510040883200000.0,precision) + cl_float(12164510040883200000.0,precision)*exp(TWO*r*xi) - 

         cl_float(22324788235240115625.0,precision)*r*xi - 

         cl_float(20320556388713831250.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

         cl_float(12225924086428552500.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

         cl_float(5467446348494130000.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - 

         cl_float(1937619942864606000.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

         cl_float(566528792821992000.0,precision)*Power(r,SIX)*Power(xi,SIX) - 

         cl_float(140462831126217600.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

         cl_float(30115609927603200.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - 

         cl_float(5663731244371200.0,precision)*Power(r,NINE)*Power(xi,NINE) - 

         cl_float(943983142502400.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - 

         cl_float(140427244339200.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)) - 

         cl_float(18723632578560.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(12.0,precision)) - 

         cl_float(2240434667520.0,precision)*Power(r,cl_float(13.0,precision))*Power(xi,cl_float(13.0,precision)) - 

         cl_float(240046571520.0,precision)*Power(r,cl_float(14.0,precision))*Power(xi,cl_float(14.0,precision)) - 

         cl_float(22861578240.0,precision)*Power(r,cl_float(15.0,precision))*Power(xi,cl_float(15.0,precision)) - 

         cl_float(1905131520.0,precision)*Power(r,cl_float(16.0,precision))*Power(xi,cl_float(16.0,precision)) - 

         cl_float(134479872.0,precision)*Power(r,cl_float(17.0,precision))*Power(xi,cl_float(17.0,precision)) - 

         cl_float(7471104.0,precision)*Power(r,cl_float(18.0,precision))*Power(xi,cl_float(18.0,precision)) - cl_float(262144.0,precision)*Power(r,cl_float(19.0,precision))*Power(xi,cl_float(19.0,precision)))/

       (cl_float(1.21645100408832e19,precision)*exp(TWO*r*xi)*Power(r,TWO)) + 

      (xi*(-cl_float(12164510040883200000.0,precision) + cl_float(12164510040883200000.0,precision)*exp(TWO*r*xi) - 

           cl_float(22324788235240115625.0,precision)*r*xi - 

           cl_float(20320556388713831250.0,precision)*Power(r,TWO)*Power(xi,TWO) - 

           cl_float(12225924086428552500.0,precision)*Power(r,THREE)*Power(xi,THREE) - 

           cl_float(5467446348494130000.0,precision)*Power(r,FOUR)*Power(xi,FOUR) - 

           cl_float(1937619942864606000.0,precision)*Power(r,FIVE)*Power(xi,FIVE) - 

           cl_float(566528792821992000.0,precision)*Power(r,SIX)*Power(xi,SIX) - 

           cl_float(140462831126217600.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) - 

           cl_float(30115609927603200.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) - 

           cl_float(5663731244371200.0,precision)*Power(r,NINE)*Power(xi,NINE) - 

           cl_float(943983142502400.0,precision)*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)) - 

           cl_float(140427244339200.0,precision)*Power(r,cl_float(11.0,precision))*Power(xi,cl_float(11.0,precision)) - 

           cl_float(18723632578560.0,precision)*Power(r,cl_float(12.0,precision))*Power(xi,cl_float(12.0,precision)) - 

           cl_float(2240434667520.0,precision)*Power(r,cl_float(13.0,precision))*Power(xi,cl_float(13.0,precision)) - 

           cl_float(240046571520.0,precision)*Power(r,cl_float(14.0,precision))*Power(xi,cl_float(14.0,precision)) - 

           cl_float(22861578240.0,precision)*Power(r,cl_float(15.0,precision))*Power(xi,cl_float(15.0,precision)) - 

           cl_float(1905131520.0,precision)*Power(r,cl_float(16.0,precision))*Power(xi,cl_float(16.0,precision)) - 

           cl_float(134479872.0,precision)*Power(r,cl_float(17.0,precision))*Power(xi,cl_float(17.0,precision)) - 

           cl_float(7471104.0,precision)*Power(r,cl_float(18.0,precision))*Power(xi,cl_float(18.0,precision)) - cl_float(262144.0,precision)*Power(r,cl_float(19.0,precision))*Power(xi,cl_float(19.0,precision))))/

       (cl_float(6.0822550204416e18,precision)*exp(TWO*r*xi)*r)

    ; }
 
  }
  else {
      if (r == ZERO) {  S = ZERO

    ; } else {  S = (cl_float(70875.0,precision)*exp(TWO*r*(xi + xj))*Power(Power(xi,TWO) - Power(xj,TWO),cl_float(19.0,precision)) + 

         exp(TWO*r*xj)*Power(xj,cl_float(12.0,precision))*

          (-cl_float(630.0,precision)*Power(r,EIGHT)*Power(xi,cl_float(34.0,precision)) - cl_float(10.0,precision)*Power(r,NINE)*Power(xi,cl_float(35.0,precision)) + 

            cl_float(70875.0,precision)*Power(xj,cl_float(26.0,precision)) + cl_float(127575.0,precision)*r*xi*Power(xj,cl_float(26.0,precision)) - 

            cl_float(30.0,precision)*Power(r,SEVEN)*Power(xi,cl_float(33.0,precision))*(cl_float(630.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(14175.0,precision)*Power(xi,TWO)*Power(xj,cl_float(24.0,precision))*(-cl_float(95.0,precision) + EIGHT*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(4725.0,precision)*r*Power(xi,THREE)*Power(xj,cl_float(24.0,precision))*

             (-cl_float(513.0,precision) + cl_float(14.0,precision)*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(90.0,precision)*Power(r,SIX)*Power(xi,cl_float(32.0,precision))*(cl_float(3920.0,precision) + cl_float(43.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(4725.0,precision)*r*Power(xi,FIVE)*Power(xj,cl_float(22.0,precision))*

             (cl_float(4617.0,precision) - cl_float(266.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(14175.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(22.0,precision))*

             (cl_float(855.0,precision) - cl_float(152.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(36.0,precision)*Power(r,FIVE)*Power(xi,cl_float(31.0,precision))*

             (-cl_float(124950.0,precision) - cl_float(4985.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(13.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(36.0,precision)*Power(r,FOUR)*Power(xi,cl_float(30.0,precision))*

             (-cl_float(1124550.0,precision) - cl_float(127960.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(863.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(135.0,precision)*r*Power(xi,SEVEN)*Power(xj,cl_float(20.0,precision))*

             (-cl_float(915705.0,precision) + cl_float(83790.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1330.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + FOUR*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(315.0,precision)*Power(xi,SIX)*Power(xj,cl_float(20.0,precision))*

             (-cl_float(218025.0,precision) + cl_float(61560.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1710.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + EIGHT*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(36.0,precision)*Power(r,THREE)*Power(xi,cl_float(29.0,precision))*

             (cl_float(7122150.0,precision) + cl_float(2102730.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(23294.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(37.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(36.0,precision)*Power(r,TWO)*Power(xi,cl_float(28.0,precision))*

             (cl_float(30523500.0,precision) + cl_float(23401350.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(299250.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(1297.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            r*Power(xi,cl_float(17.0,precision))*Power(xj,cl_float(10.0,precision))*

             (cl_float(1073961177975.0,precision) - cl_float(21753487980.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(745994340.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(5307156.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(818.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            cl_float(10.0,precision)*r*Power(xi,NINE)*Power(xj,cl_float(18.0,precision))*

             (cl_float(49448070.0,precision) - cl_float(6409935.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(161595.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(1026.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               Power(r,EIGHT)*Power(xj,EIGHT)) + 

            cl_float(90.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(18.0,precision))*

             (cl_float(3052350.0,precision) - cl_float(1220940.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(53865.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(532.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               Power(r,EIGHT)*Power(xj,EIGHT)) - 

            cl_float(1710.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(16.0,precision))*

             (cl_float(481950.0,precision) - cl_float(257040.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(16065.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(252.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               Power(r,EIGHT)*Power(xj,EIGHT)) + 

            SIX*r*Power(xi,cl_float(11.0,precision))*Power(xj,cl_float(16.0,precision))*

             (-cl_float(207559800.0,precision) + cl_float(50390550.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1165815.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(21396.0,precision)*Power(r,SIX)*Power(xj,SIX) + FIVE*Power(r,EIGHT)*Power(xj,EIGHT)) - 

            cl_float(18.0,precision)*r*Power(xi,cl_float(13.0,precision))*Power(xj,cl_float(14.0,precision))*

             (-cl_float(1703720025.0,precision) - cl_float(155669850.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(7410270.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(1532.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(26.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            cl_float(18.0,precision)*r*Power(xi,cl_float(15.0,precision))*Power(xj,cl_float(12.0,precision))*

             (cl_float(19380896325.0,precision) + cl_float(1329128850.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(7608930.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(116238.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(74.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) - 

            cl_float(18.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,cl_float(14.0,precision))*

             (cl_float(89026875.0,precision) + cl_float(179071200.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1552950.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(295820.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(146.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            cl_float(18.0,precision)*r*Power(xi,cl_float(25.0,precision))*Power(xj,TWO)*

             (-cl_float(5449970925.0,precision) - cl_float(1137574935.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(37834755.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(273062.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(171.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) - 

            NINE*r*Power(xi,cl_float(19.0,precision))*Power(xj,EIGHT)*

             (-cl_float(37914907275.0,precision) + cl_float(7613889570.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(170524620.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(397936.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(342.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) - 

            THREE*r*Power(xi,cl_float(23.0,precision))*Power(xj,FOUR)*

             (cl_float(219130630425.0,precision) - cl_float(11118046590.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(327611970.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(2920908.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(2584.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            THREE*r*Power(xi,cl_float(21.0,precision))*Power(xj,SIX)*

             (-cl_float(345162539925.0,precision) + cl_float(19030764690.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(141976170.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(1441872.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(2584.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            cl_float(63.0,precision)*Power(xi,cl_float(20.0,precision))*Power(xj,SIX)*

             (-cl_float(50980542525.0,precision) + cl_float(6240202920.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(201314310.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(956080.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(2584.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            cl_float(18.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,cl_float(12.0,precision))*

             (-cl_float(7803332775.0,precision) - cl_float(2519206200.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(119719950.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(182280.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(2734.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) - 

            cl_float(18.0,precision)*Power(xi,cl_float(26.0,precision))*(cl_float(195859125.0,precision) + cl_float(1794781800.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(67337235.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(1659700.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(4089.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            NINE*Power(xi,cl_float(18.0,precision))*Power(xj,EIGHT)*

             (-cl_float(357591274425.0,precision) + cl_float(8328390840.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(912042180.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(12842480.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(10678.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) 

    - NINE*Power(xi,cl_float(16.0,precision))*Power(xj,cl_float(10.0,precision))*

             (cl_float(128599724925.0,precision) + cl_float(21298077360.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(267928500.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(5458320.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(14722.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) 

    + cl_float(18.0,precision)*Power(xi,cl_float(24.0,precision))*Power(xj,TWO)*

             (-cl_float(7604930025.0,precision) - cl_float(8866107180.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(399272265.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(5925780.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(17651.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) 

    - NINE*Power(xi,cl_float(22.0,precision))*Power(xj,FOUR)*(cl_float(129194933175.0,precision) + 

               cl_float(3909863160.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(91420770.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(8762040.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(43928.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) 

    + Power(xi,cl_float(27.0,precision))*(-cl_float(2884470750.0,precision)*r - cl_float(6409935000.0,precision)*Power(r,THREE)*Power(xj,TWO) + 

               cl_float(28332990.0,precision)*Power(r,FIVE)*Power(xj,FOUR) + 

               cl_float(58104.0,precision)*Power(r,SEVEN)*Power(xj,SIX) + cl_float(818.0,precision)*Power(r,NINE)*Power(xj,EIGHT))) + 

         exp(TWO*r*xi)*Power(xi,cl_float(12.0,precision))*

          (Power(xi,EIGHT)*Power(xj,cl_float(18.0,precision))*

             (cl_float(3218321469825.0,precision) - cl_float(341234165475.0,precision)*r*xj - 

               cl_float(393132783960.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(57092294070.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(822786930.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(982835910.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(106664040.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(4915116.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(73602.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) - 

               cl_float(818.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(10.0,precision)*Power(xj,cl_float(26.0,precision))*(cl_float(352546425.0,precision) + cl_float(288447075.0,precision)*r*xj + 

               cl_float(109884600.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(25639740.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(4048380.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(449820.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(35280.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(1890.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(63.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(30.0,precision)*Power(xi,TWO)*Power(xj,cl_float(24.0,precision))*

             (cl_float(4562958015.0,precision) + cl_float(3269982555.0,precision)*r*xj + 

               cl_float(1076869080.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(213664500.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(28081620.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(2523276.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(153552.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(5982.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(129.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(15.0,precision)*Power(xi,cl_float(24.0,precision))*Power(xj,TWO)*

             (-cl_float(89775.0,precision) - cl_float(161595.0,precision)*r*xj - cl_float(143640.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(83790.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(35910.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(11970.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(3192.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

               cl_float(684.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - cl_float(114.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               TWO*Power(r,NINE)*Power(xj,NINE)) - 

            FIVE*Power(xi,cl_float(26.0,precision))*(cl_float(14175.0,precision) + cl_float(25515.0,precision)*r*xj + 

               cl_float(22680.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(13230.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(5670.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(1890.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(504.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(108.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + TWO*Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(1938.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,cl_float(12.0,precision))*

             (-cl_float(826875.0,precision) + cl_float(15824025.0,precision)*r*xj - cl_float(23398200.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(12344850.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(1244250.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(384930.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(59640.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

               cl_float(1848.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(84.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               FOUR*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(1938.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,cl_float(14.0,precision))*

             (cl_float(72476775.0,precision) - cl_float(180008325.0,precision)*r*xj + cl_float(98907480.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(11224710.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(4235490.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(791910.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(31080.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(2232.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(204.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               FOUR*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(342.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,cl_float(10.0,precision))*

             (cl_float(2409750.0,precision) + cl_float(3641400.0,precision)*r*xj + cl_float(9424800.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(8193150.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(6301050.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(400470.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(143640.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

               cl_float(15518.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - cl_float(281.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               NINE*Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(171.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(16.0,precision))*

             (-cl_float(6768406575.0,precision) + cl_float(6280474725.0,precision)*r*xj + 

               cl_float(438336360.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(400731030.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(74168430.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(2490810.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(461160.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(51244.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(1858.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(18.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

            NINE*Power(xi,cl_float(22.0,precision))*Power(xj,FOUR)*

             (-cl_float(1346625.0,precision) - cl_float(2423925.0,precision)*r*xj - cl_float(2154600.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1256850.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(538650.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(179550.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(47880.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(14264.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(292.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(52.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

            NINE*Power(xi,FOUR)*Power(xj,cl_float(22.0,precision))*

             (-cl_float(129194933175.0,precision) - cl_float(73043543475.0,precision)*r*xj - 

               cl_float(17732214360.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(2275149870.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(134674470.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(3148110.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(1197000.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(93176.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(3452.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(52.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

            NINE*Power(xi,SIX)*Power(xj,cl_float(20.0,precision))*

             (cl_float(356863797675.0,precision) + cl_float(115054179975.0,precision)*r*xj + 

               cl_float(3909863160.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(3706015530.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(798544530.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(75669510.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(3319400.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(6456.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(5188.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(148.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

            NINE*Power(xi,cl_float(20.0,precision))*Power(xj,SIX)*

             (-cl_float(7630875.0,precision) - cl_float(13735575.0,precision)*r*xj - cl_float(12209400.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(7122150.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(3052350.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(777210.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(591640.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(3064.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(5468.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               cl_float(148.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

            TWO*Power(xi,cl_float(18.0,precision))*Power(xj,EIGHT)*

             (-cl_float(137355750.0,precision) - cl_float(247240350.0,precision)*r*xj - cl_float(219769200.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(151171650.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(13976550.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(66692430.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(1640520.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(1046142.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(66249.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               cl_float(409.0,precision)*Power(r,NINE)*Power(xj,NINE))))/

       (cl_float(70875.0,precision)*exp(TWO*r*(xi + xj))*Power(r,TWO)*Power(xi - xj,cl_float(19.0,precision))*

         Power(xi + xj,cl_float(19.0,precision))) + (TWO*(cl_float(70875.0,precision)*exp(TWO*r*(xi + xj))*

            Power(Power(xi,TWO) - Power(xj,TWO),cl_float(19.0,precision)) + 

           exp(TWO*r*xj)*Power(xj,cl_float(12.0,precision))*

            (-cl_float(630.0,precision)*Power(r,EIGHT)*Power(xi,cl_float(34.0,precision)) - cl_float(10.0,precision)*Power(r,NINE)*Power(xi,cl_float(35.0,precision)) + 

              cl_float(70875.0,precision)*Power(xj,cl_float(26.0,precision)) + cl_float(127575.0,precision)*r*xi*Power(xj,cl_float(26.0,precision)) - 

              cl_float(30.0,precision)*Power(r,SEVEN)*Power(xi,cl_float(33.0,precision))*(cl_float(630.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

              cl_float(14175.0,precision)*Power(xi,TWO)*Power(xj,cl_float(24.0,precision))*

               (-cl_float(95.0,precision) + EIGHT*Power(r,TWO)*Power(xj,TWO)) + 

              cl_float(4725.0,precision)*r*Power(xi,THREE)*Power(xj,cl_float(24.0,precision))*

               (-cl_float(513.0,precision) + cl_float(14.0,precision)*Power(r,TWO)*Power(xj,TWO)) - 

              cl_float(90.0,precision)*Power(r,SIX)*Power(xi,cl_float(32.0,precision))*(cl_float(3920.0,precision) + cl_float(43.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

              cl_float(4725.0,precision)*r*Power(xi,FIVE)*Power(xj,cl_float(22.0,precision))*

               (cl_float(4617.0,precision) - cl_float(266.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) 

    + cl_float(14175.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(22.0,precision))*

               (cl_float(855.0,precision) - cl_float(152.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

              cl_float(36.0,precision)*Power(r,FIVE)*Power(xi,cl_float(31.0,precision))*

               (-cl_float(124950.0,precision) - cl_float(4985.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(13.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

              cl_float(36.0,precision)*Power(r,FOUR)*Power(xi,cl_float(30.0,precision))*

               (-cl_float(1124550.0,precision) - cl_float(127960.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(863.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

              cl_float(135.0,precision)*r*Power(xi,SEVEN)*Power(xj,cl_float(20.0,precision))*

               (-cl_float(915705.0,precision) + cl_float(83790.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(1330.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + FOUR*Power(r,SIX)*Power(xj,SIX)) + 

              cl_float(315.0,precision)*Power(xi,SIX)*Power(xj,cl_float(20.0,precision))*

               (-cl_float(218025.0,precision) + cl_float(61560.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(1710.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + EIGHT*Power(r,SIX)*Power(xj,SIX)) - 

              cl_float(36.0,precision)*Power(r,THREE)*Power(xi,cl_float(29.0,precision))*

               (cl_float(7122150.0,precision) + cl_float(2102730.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(23294.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(37.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

              cl_float(36.0,precision)*Power(r,TWO)*Power(xi,cl_float(28.0,precision))*

               (cl_float(30523500.0,precision) + cl_float(23401350.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(299250.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(1297.0,precision)*Power(r,SIX)*Power(xj,SIX)) 

    + r*Power(xi,cl_float(17.0,precision))*Power(xj,cl_float(10.0,precision))*

               (cl_float(1073961177975.0,precision) - cl_float(21753487980.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(745994340.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(5307156.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(818.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) 

    + cl_float(10.0,precision)*r*Power(xi,NINE)*Power(xj,cl_float(18.0,precision))*

               (cl_float(49448070.0,precision) - cl_float(6409935.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(161595.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(1026.0,precision)*Power(r,SIX)*Power(xj,SIX) + Power(r,EIGHT)*Power(xj,EIGHT)) + 

              cl_float(90.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(18.0,precision))*

               (cl_float(3052350.0,precision) - cl_float(1220940.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(53865.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(532.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

                 Power(r,EIGHT)*Power(xj,EIGHT)) - 

              cl_float(1710.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(16.0,precision))*

               (cl_float(481950.0,precision) - cl_float(257040.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(16065.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(252.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

                 Power(r,EIGHT)*Power(xj,EIGHT)) + 

              SIX*r*Power(xi,cl_float(11.0,precision))*Power(xj,cl_float(16.0,precision))*

               (-cl_float(207559800.0,precision) + cl_float(50390550.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(1165815.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(21396.0,precision)*Power(r,SIX)*Power(xj,SIX) + FIVE*Power(r,EIGHT)*Power(xj,EIGHT)) - 

              cl_float(18.0,precision)*r*Power(xi,cl_float(13.0,precision))*Power(xj,cl_float(14.0,precision))*

               (-cl_float(1703720025.0,precision) - cl_float(155669850.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(7410270.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(1532.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(26.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) + 

              cl_float(18.0,precision)*r*Power(xi,cl_float(15.0,precision))*Power(xj,cl_float(12.0,precision))*

               (cl_float(19380896325.0,precision) + cl_float(1329128850.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(7608930.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(116238.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(74.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) - 

              cl_float(18.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,cl_float(14.0,precision))*

               (cl_float(89026875.0,precision) + cl_float(179071200.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(1552950.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(295820.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(146.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) + 

              cl_float(18.0,precision)*r*Power(xi,cl_float(25.0,precision))*Power(xj,TWO)*

               (-cl_float(5449970925.0,precision) - cl_float(1137574935.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(37834755.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(273062.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(171.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) - 

              NINE*r*Power(xi,cl_float(19.0,precision))*Power(xj,EIGHT)*

               (-cl_float(37914907275.0,precision) + cl_float(7613889570.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(170524620.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(397936.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(342.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) - 

              THREE*r*Power(xi,cl_float(23.0,precision))*Power(xj,FOUR)*

               (cl_float(219130630425.0,precision) - cl_float(11118046590.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(327611970.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(2920908.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(2584.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) 

    + THREE*r*Power(xi,cl_float(21.0,precision))*Power(xj,SIX)*

               (-cl_float(345162539925.0,precision) + cl_float(19030764690.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(141976170.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(1441872.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(2584.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) 

    + cl_float(63.0,precision)*Power(xi,cl_float(20.0,precision))*Power(xj,SIX)*

               (-cl_float(50980542525.0,precision) + cl_float(6240202920.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(201314310.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(956080.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(2584.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) 

    + cl_float(18.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,cl_float(12.0,precision))*

               (-cl_float(7803332775.0,precision) - cl_float(2519206200.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(119719950.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(182280.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(2734.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) 

    - cl_float(18.0,precision)*Power(xi,cl_float(26.0,precision))*(cl_float(195859125.0,precision) + cl_float(1794781800.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(67337235.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(1659700.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(4089.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) 

    + NINE*Power(xi,cl_float(18.0,precision))*Power(xj,EIGHT)*(-cl_float(357591274425.0,precision) + 

                 cl_float(8328390840.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(912042180.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(12842480.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

                 cl_float(10678.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) - 

              NINE*Power(xi,cl_float(16.0,precision))*Power(xj,cl_float(10.0,precision))*

               (cl_float(128599724925.0,precision) + cl_float(21298077360.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(267928500.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(5458320.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(14722.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)

    ) + cl_float(18.0,precision)*Power(xi,cl_float(24.0,precision))*Power(xj,TWO)*

               (-cl_float(7604930025.0,precision) - cl_float(8866107180.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(399272265.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(5925780.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(17651.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)

    ) - NINE*Power(xi,cl_float(22.0,precision))*Power(xj,FOUR)*

               (cl_float(129194933175.0,precision) + cl_float(3909863160.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(91420770.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(8762040.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(43928.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)

    ) + Power(xi,cl_float(27.0,precision))*(-cl_float(2884470750.0,precision)*r - cl_float(6409935000.0,precision)*Power(r,THREE)*Power(xj,TWO) + 

                 cl_float(28332990.0,precision)*Power(r,FIVE)*Power(xj,FOUR) + 

                 cl_float(58104.0,precision)*Power(r,SEVEN)*Power(xj,SIX) + cl_float(818.0,precision)*Power(r,NINE)*Power(xj,EIGHT))) + 

           exp(TWO*r*xi)*Power(xi,cl_float(12.0,precision))*

            (Power(xi,EIGHT)*Power(xj,cl_float(18.0,precision))*

               (cl_float(3218321469825.0,precision) - cl_float(341234165475.0,precision)*r*xj - 

                 cl_float(393132783960.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(57092294070.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

                 cl_float(822786930.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(982835910.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

                 cl_float(106664040.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

                 cl_float(4915116.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

                 cl_float(73602.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) - cl_float(818.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

              cl_float(10.0,precision)*Power(xj,cl_float(26.0,precision))*(cl_float(352546425.0,precision) + cl_float(288447075.0,precision)*r*xj + 

                 cl_float(109884600.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(25639740.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

                 cl_float(4048380.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(449820.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

                 cl_float(35280.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(1890.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

                 cl_float(63.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + Power(r,NINE)*Power(xj,NINE)) + 

              cl_float(30.0,precision)*Power(xi,TWO)*Power(xj,cl_float(24.0,precision))*

               (cl_float(4562958015.0,precision) + cl_float(3269982555.0,precision)*r*xj + 

                 cl_float(1076869080.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(213664500.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

                 cl_float(28081620.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(2523276.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

                 cl_float(153552.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(5982.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

                 cl_float(129.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + Power(r,NINE)*Power(xj,NINE)) - 

              cl_float(15.0,precision)*Power(xi,cl_float(24.0,precision))*Power(xj,TWO)*

               (-cl_float(89775.0,precision) - cl_float(161595.0,precision)*r*xj - cl_float(143640.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(83790.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(35910.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(11970.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(3192.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

                 cl_float(684.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - cl_float(114.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

                 TWO*Power(r,NINE)*Power(xj,NINE)) - 

              FIVE*Power(xi,cl_float(26.0,precision))*(cl_float(14175.0,precision) + cl_float(25515.0,precision)*r*xj + 

                 cl_float(22680.0,precision)*Power(r,TWO)*Power(xj,TWO) + cl_float(13230.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

                 cl_float(5670.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(1890.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

                 cl_float(504.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(108.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

                 cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + TWO*Power(r,NINE)*Power(xj,NINE)) - 

              cl_float(1938.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,cl_float(12.0,precision))*

               (-cl_float(826875.0,precision) + cl_float(15824025.0,precision)*r*xj - cl_float(23398200.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(12344850.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

                 cl_float(1244250.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(384930.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

                 cl_float(59640.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(1848.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

                 cl_float(84.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + FOUR*Power(r,NINE)*Power(xj,NINE)) + 

              cl_float(1938.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,cl_float(14.0,precision))*

               (cl_float(72476775.0,precision) - cl_float(180008325.0,precision)*r*xj + cl_float(98907480.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

                 cl_float(11224710.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

                 cl_float(4235490.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(791910.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

                 cl_float(31080.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(2232.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

                 cl_float(204.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + FOUR*Power(r,NINE)*Power(xj,NINE)) + 

              cl_float(342.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,cl_float(10.0,precision))*

               (cl_float(2409750.0,precision) + cl_float(3641400.0,precision)*r*xj + cl_float(9424800.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(8193150.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

                 cl_float(6301050.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(400470.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

                 cl_float(143640.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

                 cl_float(15518.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - cl_float(281.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

                 NINE*Power(r,NINE)*Power(xj,NINE)) - 

              cl_float(171.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(16.0,precision))*

               (-cl_float(6768406575.0,precision) + cl_float(6280474725.0,precision)*r*xj + 

                 cl_float(438336360.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(400731030.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

                 cl_float(74168430.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(2490810.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

                 cl_float(461160.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

                 cl_float(51244.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(1858.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

                 cl_float(18.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

              NINE*Power(xi,cl_float(22.0,precision))*Power(xj,FOUR)*

               (-cl_float(1346625.0,precision) - cl_float(2423925.0,precision)*r*xj - cl_float(2154600.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(1256850.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

                 cl_float(538650.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(179550.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

                 cl_float(47880.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(14264.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

                 cl_float(292.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(52.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

              NINE*Power(xi,FOUR)*Power(xj,cl_float(22.0,precision))*

               (-cl_float(129194933175.0,precision) - cl_float(73043543475.0,precision)*r*xj - 

                 cl_float(17732214360.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(2275149870.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

                 cl_float(134674470.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

                 cl_float(3148110.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

                 cl_float(1197000.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

                 cl_float(93176.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(3452.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

                 cl_float(52.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

              NINE*Power(xi,SIX)*Power(xj,cl_float(20.0,precision))*

               (cl_float(356863797675.0,precision) + cl_float(115054179975.0,precision)*r*xj + 

                 cl_float(3909863160.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(3706015530.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

                 cl_float(798544530.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(75669510.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

                 cl_float(3319400.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

                 cl_float(6456.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(5188.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

                 cl_float(148.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

              NINE*Power(xi,cl_float(20.0,precision))*Power(xj,SIX)*

               (-cl_float(7630875.0,precision) - cl_float(13735575.0,precision)*r*xj - cl_float(12209400.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(7122150.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

                 cl_float(3052350.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(777210.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

                 cl_float(591640.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(3064.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

                 cl_float(5468.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(148.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

              TWO*Power(xi,cl_float(18.0,precision))*Power(xj,EIGHT)*

               (-cl_float(137355750.0,precision) - cl_float(247240350.0,precision)*r*xj - 

                 cl_float(219769200.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

                 cl_float(151171650.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

                 cl_float(13976550.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

                 cl_float(66692430.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

                 cl_float(1640520.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

                 cl_float(1046142.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

                 cl_float(66249.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(409.0,precision)*Power(r,NINE)*Power(xj,NINE)))))/

       (cl_float(70875.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,cl_float(19.0,precision))*Power(xi + xj,cl_float(18.0,precision))) - 

      (cl_float(141750.0,precision)*exp(TWO*r*(xi + xj))*(xi + xj)*

          Power(Power(xi,TWO) - Power(xj,TWO),cl_float(19.0,precision)) + 

         exp(TWO*r*xj)*Power(xj,cl_float(12.0,precision))*

          (-cl_float(5040.0,precision)*Power(r,SEVEN)*Power(xi,cl_float(34.0,precision)) - cl_float(90.0,precision)*Power(r,EIGHT)*Power(xi,cl_float(35.0,precision)) - 

            cl_float(7740.0,precision)*Power(r,SEVEN)*Power(xi,cl_float(32.0,precision))*Power(xj,TWO) - 

            cl_float(60.0,precision)*Power(r,EIGHT)*Power(xi,cl_float(33.0,precision))*Power(xj,TWO) + cl_float(127575.0,precision)*xi*Power(xj,cl_float(26.0,precision)) + 

            cl_float(226800.0,precision)*r*Power(xi,TWO)*Power(xj,cl_float(26.0,precision)) + 

            cl_float(132300.0,precision)*Power(r,TWO)*Power(xi,THREE)*Power(xj,cl_float(26.0,precision)) - 

            cl_float(210.0,precision)*Power(r,SIX)*Power(xi,cl_float(33.0,precision))*(cl_float(630.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(4725.0,precision)*Power(xi,THREE)*Power(xj,cl_float(24.0,precision))*(-cl_float(513.0,precision) + cl_float(14.0,precision)*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(540.0,precision)*Power(r,FIVE)*Power(xi,cl_float(32.0,precision))*(cl_float(3920.0,precision) + cl_float(43.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(4725.0,precision)*r*Power(xi,FIVE)*Power(xj,cl_float(22.0,precision))*

             (-cl_float(532.0,precision)*r*Power(xj,TWO) + EIGHT*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(14175.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(22.0,precision))*

             (-cl_float(304.0,precision)*r*Power(xj,TWO) + EIGHT*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(36.0,precision)*Power(r,FIVE)*Power(xi,cl_float(31.0,precision))*

             (-cl_float(9970.0,precision)*r*Power(xj,TWO) + cl_float(52.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(36.0,precision)*Power(r,FOUR)*Power(xi,cl_float(30.0,precision))*

             (-cl_float(255920.0,precision)*r*Power(xj,TWO) + cl_float(3452.0,precision)*Power(r,THREE)*Power(xj,FOUR)) + 

            cl_float(4725.0,precision)*Power(xi,FIVE)*Power(xj,cl_float(22.0,precision))*

             (cl_float(4617.0,precision) - cl_float(266.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(180.0,precision)*Power(r,FOUR)*Power(xi,cl_float(31.0,precision))*

             (-cl_float(124950.0,precision) - cl_float(4985.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(13.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(144.0,precision)*Power(r,THREE)*Power(xi,cl_float(30.0,precision))*

             (-cl_float(1124550.0,precision) - cl_float(127960.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(863.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(135.0,precision)*r*Power(xi,SEVEN)*Power(xj,cl_float(20.0,precision))*

             (cl_float(167580.0,precision)*r*Power(xj,TWO) - cl_float(5320.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(24.0,precision)*Power(r,FIVE)*Power(xj,SIX)) + 

            cl_float(315.0,precision)*Power(xi,SIX)*Power(xj,cl_float(20.0,precision))*

             (cl_float(123120.0,precision)*r*Power(xj,TWO) - cl_float(6840.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(48.0,precision)*Power(r,FIVE)*Power(xj,SIX)) - 

            cl_float(36.0,precision)*Power(r,THREE)*Power(xi,cl_float(29.0,precision))*

             (cl_float(4205460.0,precision)*r*Power(xj,TWO) - cl_float(93176.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(222.0,precision)*Power(r,FIVE)*Power(xj,SIX)) - 

            cl_float(36.0,precision)*Power(r,TWO)*Power(xi,cl_float(28.0,precision))*

             (cl_float(46802700.0,precision)*r*Power(xj,TWO) - cl_float(1197000.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(7782.0,precision)*Power(r,FIVE)*Power(xj,SIX)) + 

            cl_float(135.0,precision)*Power(xi,SEVEN)*Power(xj,cl_float(20.0,precision))*

             (-cl_float(915705.0,precision) + cl_float(83790.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1330.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + FOUR*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(108.0,precision)*Power(r,TWO)*Power(xi,cl_float(29.0,precision))*

             (cl_float(7122150.0,precision) + cl_float(2102730.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(23294.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(37.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(72.0,precision)*r*Power(xi,cl_float(28.0,precision))*(cl_float(30523500.0,precision) + cl_float(23401350.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(299250.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(1297.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            r*Power(xi,cl_float(17.0,precision))*Power(xj,cl_float(10.0,precision))*

             (-cl_float(43506975960.0,precision)*r*Power(xj,TWO) - cl_float(2983977360.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(31842936.0,precision)*Power(r,FIVE)*Power(xj,SIX) - cl_float(6544.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT)) + 

            cl_float(10.0,precision)*r*Power(xi,NINE)*Power(xj,cl_float(18.0,precision))*

             (-cl_float(12819870.0,precision)*r*Power(xj,TWO) + cl_float(646380.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(6156.0,precision)*Power(r,FIVE)*Power(xj,SIX) + EIGHT*Power(r,SEVEN)*Power(xj,EIGHT)) + 

            cl_float(90.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(18.0,precision))*

             (-cl_float(2441880.0,precision)*r*Power(xj,TWO) + cl_float(215460.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(3192.0,precision)*Power(r,FIVE)*Power(xj,SIX) + EIGHT*Power(r,SEVEN)*Power(xj,EIGHT)) - 

            cl_float(1710.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(16.0,precision))*

             (-cl_float(514080.0,precision)*r*Power(xj,TWO) + cl_float(64260.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(1512.0,precision)*Power(r,FIVE)*Power(xj,SIX) + EIGHT*Power(r,SEVEN)*Power(xj,EIGHT)) + 

            SIX*r*Power(xi,cl_float(11.0,precision))*Power(xj,cl_float(16.0,precision))*

             (cl_float(100781100.0,precision)*r*Power(xj,TWO) - cl_float(4663260.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(128376.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(40.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT)) - 

            cl_float(18.0,precision)*r*Power(xi,cl_float(13.0,precision))*Power(xj,cl_float(14.0,precision))*

             (-cl_float(311339700.0,precision)*r*Power(xj,TWO) - cl_float(29641080.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(9192.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(208.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT)) + 

            cl_float(18.0,precision)*r*Power(xi,cl_float(15.0,precision))*Power(xj,cl_float(12.0,precision))*

             (cl_float(2658257700.0,precision)*r*Power(xj,TWO) - cl_float(30435720.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(697428.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(592.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT)) - 

            cl_float(18.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,cl_float(14.0,precision))*

             (cl_float(358142400.0,precision)*r*Power(xj,TWO) + cl_float(6211800.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(1774920.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(1168.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT)) + 

            cl_float(18.0,precision)*r*Power(xi,cl_float(25.0,precision))*Power(xj,TWO)*

             (-cl_float(2275149870.0,precision)*r*Power(xj,TWO) + cl_float(151339020.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(1638372.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(1368.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT)) - 

            NINE*r*Power(xi,cl_float(19.0,precision))*Power(xj,EIGHT)*

             (cl_float(15227779140.0,precision)*r*Power(xj,TWO) - cl_float(682098480.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(2387616.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(2736.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT)) - 

            THREE*r*Power(xi,cl_float(23.0,precision))*Power(xj,FOUR)*

             (-cl_float(22236093180.0,precision)*r*Power(xj,TWO) + cl_float(1310447880.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(17525448.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(20672.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT)) 

    + THREE*r*Power(xi,cl_float(21.0,precision))*Power(xj,SIX)*

             (cl_float(38061529380.0,precision)*r*Power(xj,TWO) - cl_float(567904680.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(8651232.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(20672.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT)) + 

            cl_float(63.0,precision)*Power(xi,cl_float(20.0,precision))*Power(xj,SIX)*

             (cl_float(12480405840.0,precision)*r*Power(xj,TWO) - cl_float(805257240.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(5736480.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(20672.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT)) + 

            cl_float(18.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,cl_float(12.0,precision))*

             (-cl_float(5038412400.0,precision)*r*Power(xj,TWO) - cl_float(478879800.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(1093680.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(21872.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT)) - 

            cl_float(18.0,precision)*Power(xi,cl_float(26.0,precision))*(cl_float(3589563600.0,precision)*r*Power(xj,TWO) + 

               cl_float(269348940.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(9958200.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(32712.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT)) + 

            NINE*Power(xi,cl_float(18.0,precision))*Power(xj,EIGHT)*

             (cl_float(16656781680.0,precision)*r*Power(xj,TWO) + cl_float(3648168720.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(77054880.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(85424.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT)) 

    - NINE*Power(xi,cl_float(16.0,precision))*Power(xj,cl_float(10.0,precision))*(cl_float(42596154720.0,precision)*r*Power(xj,TWO) - 

               cl_float(1071714000.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(32749920.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(117776.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT)) 

    + cl_float(18.0,precision)*Power(xi,cl_float(24.0,precision))*Power(xj,TWO)*(-cl_float(17732214360.0,precision)*r*Power(xj,TWO) + 

               cl_float(1597089060.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(35554680.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(141208.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT)) 

    - NINE*Power(xi,cl_float(22.0,precision))*Power(xj,FOUR)*(cl_float(7819726320.0,precision)*r*Power(xj,TWO) + 

               cl_float(365683080.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(52572240.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(351424.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT)) 

    + Power(xi,cl_float(17.0,precision))*Power(xj,cl_float(10.0,precision))*(cl_float(1073961177975.0,precision) - 

               cl_float(21753487980.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(745994340.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(5307156.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(818.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            cl_float(10.0,precision)*Power(xi,NINE)*Power(xj,cl_float(18.0,precision))*

             (cl_float(49448070.0,precision) - cl_float(6409935.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(161595.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(1026.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               Power(r,EIGHT)*Power(xj,EIGHT)) + 

            SIX*Power(xi,cl_float(11.0,precision))*Power(xj,cl_float(16.0,precision))*

             (-cl_float(207559800.0,precision) + cl_float(50390550.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1165815.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(21396.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               FIVE*Power(r,EIGHT)*Power(xj,EIGHT)) - 

            cl_float(18.0,precision)*Power(xi,cl_float(13.0,precision))*Power(xj,cl_float(14.0,precision))*

             (-cl_float(1703720025.0,precision) - cl_float(155669850.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(7410270.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(1532.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(26.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            cl_float(18.0,precision)*Power(xi,cl_float(15.0,precision))*Power(xj,cl_float(12.0,precision))*

             (cl_float(19380896325.0,precision) + cl_float(1329128850.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(7608930.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(116238.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(74.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            cl_float(18.0,precision)*Power(xi,cl_float(25.0,precision))*Power(xj,TWO)*

             (-cl_float(5449970925.0,precision) - cl_float(1137574935.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(37834755.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(273062.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(171.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) - 

            NINE*Power(xi,cl_float(19.0,precision))*Power(xj,EIGHT)*

             (-cl_float(37914907275.0,precision) + cl_float(7613889570.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(170524620.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(397936.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(342.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) - 

            THREE*Power(xi,cl_float(23.0,precision))*Power(xj,FOUR)*

             (cl_float(219130630425.0,precision) - cl_float(11118046590.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(327611970.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(2920908.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(2584.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            THREE*Power(xi,cl_float(21.0,precision))*Power(xj,SIX)*

             (-cl_float(345162539925.0,precision) + cl_float(19030764690.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(141976170.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(1441872.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(2584.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            Power(xi,cl_float(27.0,precision))*(-cl_float(2884470750.0,precision) - cl_float(19229805000.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(141664950.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(406728.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(7362.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT))) + 

         TWO*exp(TWO*r*xj)*Power(xj,cl_float(13.0,precision))*

          (-cl_float(630.0,precision)*Power(r,EIGHT)*Power(xi,cl_float(34.0,precision)) - cl_float(10.0,precision)*Power(r,NINE)*Power(xi,cl_float(35.0,precision)) + 

            cl_float(70875.0,precision)*Power(xj,cl_float(26.0,precision)) + cl_float(127575.0,precision)*r*xi*Power(xj,cl_float(26.0,precision)) - 

            cl_float(30.0,precision)*Power(r,SEVEN)*Power(xi,cl_float(33.0,precision))*(cl_float(630.0,precision) + Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(14175.0,precision)*Power(xi,TWO)*Power(xj,cl_float(24.0,precision))*(-cl_float(95.0,precision) + EIGHT*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(4725.0,precision)*r*Power(xi,THREE)*Power(xj,cl_float(24.0,precision))*

             (-cl_float(513.0,precision) + cl_float(14.0,precision)*Power(r,TWO)*Power(xj,TWO)) - 

            cl_float(90.0,precision)*Power(r,SIX)*Power(xi,cl_float(32.0,precision))*(cl_float(3920.0,precision) + cl_float(43.0,precision)*Power(r,TWO)*Power(xj,TWO)) + 

            cl_float(4725.0,precision)*r*Power(xi,FIVE)*Power(xj,cl_float(22.0,precision))*

             (cl_float(4617.0,precision) - cl_float(266.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(14175.0,precision)*Power(xi,FOUR)*Power(xj,cl_float(22.0,precision))*

             (cl_float(855.0,precision) - cl_float(152.0,precision)*Power(r,TWO)*Power(xj,TWO) + TWO*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(36.0,precision)*Power(r,FIVE)*Power(xi,cl_float(31.0,precision))*

             (-cl_float(124950.0,precision) - cl_float(4985.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(13.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(36.0,precision)*Power(r,FOUR)*Power(xi,cl_float(30.0,precision))*

             (-cl_float(1124550.0,precision) - cl_float(127960.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(863.0,precision)*Power(r,FOUR)*Power(xj,FOUR)) + 

            cl_float(135.0,precision)*r*Power(xi,SEVEN)*Power(xj,cl_float(20.0,precision))*

             (-cl_float(915705.0,precision) + cl_float(83790.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1330.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + FOUR*Power(r,SIX)*Power(xj,SIX)) + 

            cl_float(315.0,precision)*Power(xi,SIX)*Power(xj,cl_float(20.0,precision))*

             (-cl_float(218025.0,precision) + cl_float(61560.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1710.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + EIGHT*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(36.0,precision)*Power(r,THREE)*Power(xi,cl_float(29.0,precision))*

             (cl_float(7122150.0,precision) + cl_float(2102730.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(23294.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(37.0,precision)*Power(r,SIX)*Power(xj,SIX)) - 

            cl_float(36.0,precision)*Power(r,TWO)*Power(xi,cl_float(28.0,precision))*

             (cl_float(30523500.0,precision) + cl_float(23401350.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(299250.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(1297.0,precision)*Power(r,SIX)*Power(xj,SIX)) + 

            r*Power(xi,cl_float(17.0,precision))*Power(xj,cl_float(10.0,precision))*

             (cl_float(1073961177975.0,precision) - cl_float(21753487980.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(745994340.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(5307156.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(818.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            cl_float(10.0,precision)*r*Power(xi,NINE)*Power(xj,cl_float(18.0,precision))*

             (cl_float(49448070.0,precision) - cl_float(6409935.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(161595.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(1026.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               Power(r,EIGHT)*Power(xj,EIGHT)) + 

            cl_float(90.0,precision)*Power(xi,EIGHT)*Power(xj,cl_float(18.0,precision))*

             (cl_float(3052350.0,precision) - cl_float(1220940.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(53865.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(532.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               Power(r,EIGHT)*Power(xj,EIGHT)) - 

            cl_float(1710.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(16.0,precision))*

             (cl_float(481950.0,precision) - cl_float(257040.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(16065.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(252.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               Power(r,EIGHT)*Power(xj,EIGHT)) + 

            SIX*r*Power(xi,cl_float(11.0,precision))*Power(xj,cl_float(16.0,precision))*

             (-cl_float(207559800.0,precision) + cl_float(50390550.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1165815.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(21396.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               FIVE*Power(r,EIGHT)*Power(xj,EIGHT)) - 

            cl_float(18.0,precision)*r*Power(xi,cl_float(13.0,precision))*Power(xj,cl_float(14.0,precision))*

             (-cl_float(1703720025.0,precision) - cl_float(155669850.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(7410270.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(1532.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(26.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            cl_float(18.0,precision)*r*Power(xi,cl_float(15.0,precision))*Power(xj,cl_float(12.0,precision))*

             (cl_float(19380896325.0,precision) + cl_float(1329128850.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(7608930.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(116238.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(74.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) - 

            cl_float(18.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,cl_float(14.0,precision))*

             (cl_float(89026875.0,precision) + cl_float(179071200.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(1552950.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(295820.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(146.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            cl_float(18.0,precision)*r*Power(xi,cl_float(25.0,precision))*Power(xj,TWO)*

             (-cl_float(5449970925.0,precision) - cl_float(1137574935.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(37834755.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(273062.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(171.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) - 

            NINE*r*Power(xi,cl_float(19.0,precision))*Power(xj,EIGHT)*

             (-cl_float(37914907275.0,precision) + cl_float(7613889570.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(170524620.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(397936.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(342.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) - 

            THREE*r*Power(xi,cl_float(23.0,precision))*Power(xj,FOUR)*

             (cl_float(219130630425.0,precision) - cl_float(11118046590.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(327611970.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(2920908.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(2584.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            THREE*r*Power(xi,cl_float(21.0,precision))*Power(xj,SIX)*

             (-cl_float(345162539925.0,precision) + cl_float(19030764690.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(141976170.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(1441872.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(2584.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            cl_float(63.0,precision)*Power(xi,cl_float(20.0,precision))*Power(xj,SIX)*

             (-cl_float(50980542525.0,precision) + cl_float(6240202920.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(201314310.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(956080.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(2584.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            cl_float(18.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,cl_float(12.0,precision))*

             (-cl_float(7803332775.0,precision) - cl_float(2519206200.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(119719950.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(182280.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(2734.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) - 

            cl_float(18.0,precision)*Power(xi,cl_float(26.0,precision))*(cl_float(195859125.0,precision) + cl_float(1794781800.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(67337235.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(1659700.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(4089.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            NINE*Power(xi,cl_float(18.0,precision))*Power(xj,EIGHT)*

             (-cl_float(357591274425.0,precision) + cl_float(8328390840.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(912042180.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(12842480.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(10678.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) 

    - NINE*Power(xi,cl_float(16.0,precision))*Power(xj,cl_float(10.0,precision))*(cl_float(128599724925.0,precision) + 

               cl_float(21298077360.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(267928500.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(5458320.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(14722.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            cl_float(18.0,precision)*Power(xi,cl_float(24.0,precision))*Power(xj,TWO)*

             (-cl_float(7604930025.0,precision) - cl_float(8866107180.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(399272265.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(5925780.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(17651.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) - 

            NINE*Power(xi,cl_float(22.0,precision))*Power(xj,FOUR)*

             (cl_float(129194933175.0,precision) + cl_float(3909863160.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(91420770.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(8762040.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(43928.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT)) + 

            Power(xi,cl_float(27.0,precision))*(-cl_float(2884470750.0,precision)*r - cl_float(6409935000.0,precision)*Power(r,THREE)*Power(xj,TWO) + 

               cl_float(28332990.0,precision)*Power(r,FIVE)*Power(xj,FOUR) + cl_float(58104.0,precision)*Power(r,SEVEN)*Power(xj,SIX) + 

               cl_float(818.0,precision)*Power(r,NINE)*Power(xj,EIGHT))) + 

         exp(TWO*r*xi)*Power(xi,cl_float(12.0,precision))*

          (Power(xi,EIGHT)*Power(xj,cl_float(18.0,precision))*

             (-cl_float(341234165475.0,precision)*xj - cl_float(786265567920.0,precision)*r*Power(xj,TWO) - 

               cl_float(171276882210.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(3291147720.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(4914179550.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(639984240.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(34405812.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + 

               cl_float(588816.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) - cl_float(7362.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            cl_float(10.0,precision)*Power(xj,cl_float(26.0,precision))*(cl_float(288447075.0,precision)*xj + cl_float(219769200.0,precision)*r*Power(xj,TWO) + 

               cl_float(76919220.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(16193520.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(2249100.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(211680.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(13230.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + 

               cl_float(504.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + NINE*Power(r,EIGHT)*Power(xj,NINE)) + 

            cl_float(30.0,precision)*Power(xi,TWO)*Power(xj,cl_float(24.0,precision))*

             (cl_float(3269982555.0,precision)*xj + cl_float(2153738160.0,precision)*r*Power(xj,TWO) + 

               cl_float(640993500.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(112326480.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(12616380.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(921312.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(41874.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + 

               cl_float(1032.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + NINE*Power(r,EIGHT)*Power(xj,NINE)) - 

            cl_float(15.0,precision)*Power(xi,cl_float(24.0,precision))*Power(xj,TWO)*

             (-cl_float(161595.0,precision)*xj - cl_float(287280.0,precision)*r*Power(xj,TWO) - 

               cl_float(251370.0,precision)*Power(r,TWO)*Power(xj,THREE) - cl_float(143640.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(59850.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - cl_float(19152.0,precision)*Power(r,FIVE)*Power(xj,SIX) - 

               cl_float(4788.0,precision)*Power(r,SIX)*Power(xj,SEVEN) - cl_float(912.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) - 

            FIVE*Power(xi,cl_float(26.0,precision))*(cl_float(25515.0,precision)*xj + cl_float(45360.0,precision)*r*Power(xj,TWO) + 

               cl_float(39690.0,precision)*Power(r,TWO)*Power(xj,THREE) + cl_float(22680.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(9450.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + cl_float(3024.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(756.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(144.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) - 

            cl_float(1938.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,cl_float(12.0,precision))*

             (cl_float(15824025.0,precision)*xj - cl_float(46796400.0,precision)*r*Power(xj,TWO) + 

               cl_float(37034550.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(4977000.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(1924650.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - 

               cl_float(357840.0,precision)*Power(r,FIVE)*Power(xj,SIX) - cl_float(12936.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + 

               cl_float(672.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + cl_float(36.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            cl_float(1938.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,cl_float(14.0,precision))*

             (-cl_float(180008325.0,precision)*xj + cl_float(197814960.0,precision)*r*Power(xj,TWO) + 

               cl_float(33674130.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(16941960.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(3959550.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - 

               cl_float(186480.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(15624.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + 

               cl_float(1632.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + cl_float(36.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            cl_float(342.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,cl_float(10.0,precision))*

             (cl_float(3641400.0,precision)*xj + cl_float(18849600.0,precision)*r*Power(xj,TWO) - 

               cl_float(24579450.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(25204200.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(2002350.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - 

               cl_float(861840.0,precision)*Power(r,FIVE)*Power(xj,SIX) - cl_float(108626.0,precision)*Power(r,SIX)*Power(xj,SEVEN) - 

               cl_float(2248.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + cl_float(81.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) - 

            cl_float(171.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(16.0,precision))*

             (cl_float(6280474725.0,precision)*xj + cl_float(876672720.0,precision)*r*Power(xj,TWO) - 

               cl_float(1202193090.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(296673720.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(12454050.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(2766960.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(358708.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(14864.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(162.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            NINE*Power(xi,cl_float(22.0,precision))*Power(xj,FOUR)*

             (-cl_float(2423925.0,precision)*xj - cl_float(4309200.0,precision)*r*Power(xj,TWO) - 

               cl_float(3770550.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(2154600.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(897750.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - cl_float(287280.0,precision)*Power(r,FIVE)*Power(xj,SIX) - 

               cl_float(99848.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(2336.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(468.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) - 

            NINE*Power(xi,FOUR)*Power(xj,cl_float(22.0,precision))*

             (-cl_float(73043543475.0,precision)*xj - cl_float(35464428720.0,precision)*r*Power(xj,TWO) - 

               cl_float(6825449610.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(538697880.0,precision)*Power(r,THREE)*Power(xj,FOUR) + 

               cl_float(15740550.0,precision)*Power(r,FOUR)*Power(xj,FIVE) + 

               cl_float(7182000.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(652232.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(27616.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(468.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            NINE*Power(xi,SIX)*Power(xj,cl_float(20.0,precision))*

             (cl_float(115054179975.0,precision)*xj + cl_float(7819726320.0,precision)*r*Power(xj,TWO) - 

               cl_float(11118046590.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(3194178120.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(378347550.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - 

               cl_float(19916400.0,precision)*Power(r,FIVE)*Power(xj,SIX) - 

               cl_float(45192.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(41504.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(1332.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) - 

            NINE*Power(xi,cl_float(20.0,precision))*Power(xj,SIX)*

             (-cl_float(13735575.0,precision)*xj - cl_float(24418800.0,precision)*r*Power(xj,TWO) - 

               cl_float(21366450.0,precision)*Power(r,TWO)*Power(xj,THREE) - 

               cl_float(12209400.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(3886050.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - 

               cl_float(3549840.0,precision)*Power(r,FIVE)*Power(xj,SIX) + cl_float(21448.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + 

               cl_float(43744.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + cl_float(1332.0,precision)*Power(r,EIGHT)*Power(xj,NINE)) + 

            TWO*Power(xi,cl_float(18.0,precision))*Power(xj,EIGHT)*

             (-cl_float(247240350.0,precision)*xj - cl_float(439538400.0,precision)*r*Power(xj,TWO) - 

               cl_float(453514950.0,precision)*Power(r,TWO)*Power(xj,THREE) + 

               cl_float(55906200.0,precision)*Power(r,THREE)*Power(xj,FOUR) - 

               cl_float(333462150.0,precision)*Power(r,FOUR)*Power(xj,FIVE) - 

               cl_float(9843120.0,precision)*Power(r,FIVE)*Power(xj,SIX) + 

               cl_float(7322994.0,precision)*Power(r,SIX)*Power(xj,SEVEN) + cl_float(529992.0,precision)*Power(r,SEVEN)*Power(xj,EIGHT) + 

               cl_float(3681.0,precision)*Power(r,EIGHT)*Power(xj,NINE))) + 

         TWO*exp(TWO*r*xi)*Power(xi,cl_float(13.0,precision))*

          (Power(xi,EIGHT)*Power(xj,cl_float(18.0,precision))*

             (cl_float(3218321469825.0,precision) - cl_float(341234165475.0,precision)*r*xj - 

               cl_float(393132783960.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(57092294070.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(822786930.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(982835910.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(106664040.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(4915116.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(73602.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) - 

               cl_float(818.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(10.0,precision)*Power(xj,cl_float(26.0,precision))*(cl_float(352546425.0,precision) + cl_float(288447075.0,precision)*r*xj + 

               cl_float(109884600.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(25639740.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(4048380.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(449820.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(35280.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(1890.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(63.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(30.0,precision)*Power(xi,TWO)*Power(xj,cl_float(24.0,precision))*

             (cl_float(4562958015.0,precision) + cl_float(3269982555.0,precision)*r*xj + 

               cl_float(1076869080.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(213664500.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(28081620.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(2523276.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(153552.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(5982.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(129.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(15.0,precision)*Power(xi,cl_float(24.0,precision))*Power(xj,TWO)*

             (-cl_float(89775.0,precision) - cl_float(161595.0,precision)*r*xj - cl_float(143640.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(83790.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(35910.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(11970.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(3192.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

               cl_float(684.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - cl_float(114.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               TWO*Power(r,NINE)*Power(xj,NINE)) - 

            FIVE*Power(xi,cl_float(26.0,precision))*(cl_float(14175.0,precision) + cl_float(25515.0,precision)*r*xj + cl_float(22680.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(13230.0,precision)*Power(r,THREE)*Power(xj,THREE) + cl_float(5670.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(1890.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(504.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(108.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(18.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               TWO*Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(1938.0,precision)*Power(xi,cl_float(14.0,precision))*Power(xj,cl_float(12.0,precision))*

             (-cl_float(826875.0,precision) + cl_float(15824025.0,precision)*r*xj - cl_float(23398200.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(12344850.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(1244250.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(384930.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(59640.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(1848.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(84.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + FOUR*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(1938.0,precision)*Power(xi,cl_float(12.0,precision))*Power(xj,cl_float(14.0,precision))*

             (cl_float(72476775.0,precision) - cl_float(180008325.0,precision)*r*xj + cl_float(98907480.0,precision)*Power(r,TWO)*Power(xj,TWO) + 

               cl_float(11224710.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(4235490.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(791910.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(31080.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(2232.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(204.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + FOUR*Power(r,NINE)*Power(xj,NINE)) + 

            cl_float(342.0,precision)*Power(xi,cl_float(16.0,precision))*Power(xj,cl_float(10.0,precision))*

             (cl_float(2409750.0,precision) + cl_float(3641400.0,precision)*r*xj + cl_float(9424800.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(8193150.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(6301050.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + cl_float(400470.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(143640.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(15518.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) - 

               cl_float(281.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + NINE*Power(r,NINE)*Power(xj,NINE)) - 

            cl_float(171.0,precision)*Power(xi,cl_float(10.0,precision))*Power(xj,cl_float(16.0,precision))*

             (-cl_float(6768406575.0,precision) + cl_float(6280474725.0,precision)*r*xj + 

               cl_float(438336360.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(400731030.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(74168430.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(2490810.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + cl_float(461160.0,precision)*Power(r,SIX)*Power(xj,SIX) + 

               cl_float(51244.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(1858.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               cl_float(18.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

            NINE*Power(xi,cl_float(22.0,precision))*Power(xj,FOUR)*

             (-cl_float(1346625.0,precision) - cl_float(2423925.0,precision)*r*xj - cl_float(2154600.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(1256850.0,precision)*Power(r,THREE)*Power(xj,THREE) - cl_float(538650.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(179550.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - cl_float(47880.0,precision)*Power(r,SIX)*Power(xj,SIX) - 

               cl_float(14264.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + cl_float(292.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + 

               cl_float(52.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

            NINE*Power(xi,FOUR)*Power(xj,cl_float(22.0,precision))*

             (-cl_float(129194933175.0,precision) - cl_float(73043543475.0,precision)*r*xj - 

               cl_float(17732214360.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(2275149870.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(134674470.0,precision)*Power(r,FOUR)*Power(xj,FOUR) + 

               cl_float(3148110.0,precision)*Power(r,FIVE)*Power(xj,FIVE) + 

               cl_float(1197000.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(93176.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(3452.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(52.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

            NINE*Power(xi,SIX)*Power(xj,cl_float(20.0,precision))*

             (cl_float(356863797675.0,precision) + cl_float(115054179975.0,precision)*r*xj + 

               cl_float(3909863160.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(3706015530.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(798544530.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(75669510.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(3319400.0,precision)*Power(r,SIX)*Power(xj,SIX) - cl_float(6456.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(5188.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(148.0,precision)*Power(r,NINE)*Power(xj,NINE)) - 

            NINE*Power(xi,cl_float(20.0,precision))*Power(xj,SIX)*

             (-cl_float(7630875.0,precision) - cl_float(13735575.0,precision)*r*xj - cl_float(12209400.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(7122150.0,precision)*Power(r,THREE)*Power(xj,THREE) - 

               cl_float(3052350.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - cl_float(777210.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(591640.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(3064.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(5468.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(148.0,precision)*Power(r,NINE)*Power(xj,NINE)) + 

            TWO*Power(xi,cl_float(18.0,precision))*Power(xj,EIGHT)*

             (-cl_float(137355750.0,precision) - cl_float(247240350.0,precision)*r*xj - cl_float(219769200.0,precision)*Power(r,TWO)*Power(xj,TWO) - 

               cl_float(151171650.0,precision)*Power(r,THREE)*Power(xj,THREE) + 

               cl_float(13976550.0,precision)*Power(r,FOUR)*Power(xj,FOUR) - 

               cl_float(66692430.0,precision)*Power(r,FIVE)*Power(xj,FIVE) - 

               cl_float(1640520.0,precision)*Power(r,SIX)*Power(xj,SIX) + cl_float(1046142.0,precision)*Power(r,SEVEN)*Power(xj,SEVEN) + 

               cl_float(66249.0,precision)*Power(r,EIGHT)*Power(xj,EIGHT) + cl_float(409.0,precision)*Power(r,NINE)*Power(xj,NINE))))/

       (cl_float(70875.0,precision)*exp(TWO*r*(xi + xj))*r*Power(xi - xj,cl_float(19.0,precision))*Power(xi + xj,cl_float(19.0,precision)))

    ; }
   
  }
  return S;
}

cl_F DSlater_5S_1S(cl_F r,cl_F xi,cl_F xj)
{
  return DSlater_1S_5S(r,xj,xi);
}

cl_F DSlater_5S_2S(cl_F r,cl_F xi,cl_F xj)
{
  return DSlater_2S_5S(r,xj,xi);
}

cl_F DSlater_5S_3S(cl_F r,cl_F xi,cl_F xj)
{
  return DSlater_3S_5S(r,xj,xi);
}

cl_F DSlater_5S_4S(cl_F r,cl_F xi,cl_F xj)
{
  return DSlater_4S_5S(r,xj,xi);
}

cl_F Nuclear_1S(cl_F r,cl_F xi)
{
  cl_F S = ZERO;
    S = ONE/r - (ONE + r*xi)/(exp(TWO*r*xi)*r)

    ;

  return S;
}

cl_F Nuclear_2S(cl_F r,cl_F xi)
{
  cl_F S = ZERO;
    S = ONE/r - (SIX + NINE*r*xi + SIX*Power(r,TWO)*Power(xi,TWO) + TWO*Power(r,THREE)*Power(xi,THREE))/

     (cl_float(6.0,precision)*exp(TWO*r*xi)*r)

    ;

  return S;
}

cl_F Nuclear_3S(cl_F r,cl_F xi)
{
  cl_F S = ZERO;
    S = ONE/r - (cl_float(45.0,precision) + cl_float(75.0,precision)*r*xi + cl_float(60.0,precision)*Power(r,TWO)*Power(xi,TWO) + 

       cl_float(30.0,precision)*Power(r,THREE)*Power(xi,THREE) + cl_float(10.0,precision)*Power(r,FOUR)*Power(xi,FOUR) + 

       TWO*Power(r,FIVE)*Power(xi,FIVE))/(cl_float(45.0,precision)*exp(TWO*r*xi)*r)

    ;

  return S;
}

cl_F Nuclear_4S(cl_F r,cl_F xi)
{
  cl_F S = ZERO;
    S = ONE/r - (cl_float(1260.0,precision) + cl_float(2205.0,precision)*r*xi + cl_float(1890.0,precision)*Power(r,TWO)*Power(xi,TWO) + 

       cl_float(1050.0,precision)*Power(r,THREE)*Power(xi,THREE) + cl_float(420.0,precision)*Power(r,FOUR)*Power(xi,FOUR) + 

       cl_float(126.0,precision)*Power(r,FIVE)*Power(xi,FIVE) + cl_float(28.0,precision)*Power(r,SIX)*Power(xi,SIX) + 

       FOUR*Power(r,SEVEN)*Power(xi,SEVEN))/(cl_float(1260.0,precision)*exp(TWO*r*xi)*r)

    ;

  return S;
}

cl_F Nuclear_5S(cl_F r,cl_F xi)
{
  cl_F S = ZERO;
    S = ONE/r - (cl_float(14175.0,precision) + cl_float(25515.0,precision)*r*xi + cl_float(22680.0,precision)*Power(r,TWO)*Power(xi,TWO) + 

       cl_float(13230.0,precision)*Power(r,THREE)*Power(xi,THREE) + cl_float(5670.0,precision)*Power(r,FOUR)*Power(xi,FOUR) + 

       cl_float(1890.0,precision)*Power(r,FIVE)*Power(xi,FIVE) + cl_float(504.0,precision)*Power(r,SIX)*Power(xi,SIX) + 

       cl_float(108.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) + cl_float(18.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) + 

       TWO*Power(r,NINE)*Power(xi,NINE))/(cl_float(14175.0,precision)*exp(TWO*r*xi)*r)

    ;

  return S;
}

cl_F DNuclear_1S(cl_F r,cl_F xi)
{
  cl_F S = ZERO;
    S = Power(r,-TWO) - (ONE + TWO*r*xi + TWO*Power(r,TWO)*Power(xi,TWO))/

     (exp(TWO*r*xi)*Power(r,TWO))

    ;

  return S;
}

cl_F DNuclear_2S(cl_F r,cl_F xi)
{
  cl_F S = ZERO;
    S = Power(r,-TWO) - (THREE + SIX*r*xi + SIX*Power(r,TWO)*Power(xi,TWO) + 

       FOUR*Power(r,THREE)*Power(xi,THREE) + TWO*Power(r,FOUR)*Power(xi,FOUR))/

     (cl_float(3.0,precision)*exp(TWO*r*xi)*Power(r,TWO))

    ;

  return S;
}

cl_F DNuclear_3S(cl_F r,cl_F xi)
{
  cl_F S = ZERO;
    S = Power(r,-TWO) - (cl_float(45.0,precision) + cl_float(90.0,precision)*r*xi + cl_float(90.0,precision)*Power(r,TWO)*Power(xi,TWO) + 

       cl_float(60.0,precision)*Power(r,THREE)*Power(xi,THREE) + cl_float(30.0,precision)*Power(r,FOUR)*Power(xi,FOUR) + 

       cl_float(12.0,precision)*Power(r,FIVE)*Power(xi,FIVE) + FOUR*Power(r,SIX)*Power(xi,SIX))/

     (cl_float(45.0,precision)*exp(TWO*r*xi)*Power(r,TWO))

    ;

  return S;
}

cl_F DNuclear_4S(cl_F r,cl_F xi)
{
  cl_F S = ZERO;
    S = Power(r,-TWO) - (cl_float(315.0,precision) + cl_float(630.0,precision)*r*xi + cl_float(630.0,precision)*Power(r,TWO)*Power(xi,TWO) + 

       cl_float(420.0,precision)*Power(r,THREE)*Power(xi,THREE) + cl_float(210.0,precision)*Power(r,FOUR)*Power(xi,FOUR) + 

       cl_float(84.0,precision)*Power(r,FIVE)*Power(xi,FIVE) + cl_float(28.0,precision)*Power(r,SIX)*Power(xi,SIX) + 

       EIGHT*Power(r,SEVEN)*Power(xi,SEVEN) + TWO*Power(r,EIGHT)*Power(xi,EIGHT))/

     (cl_float(315.0,precision)*exp(TWO*r*xi)*Power(r,TWO))

    ;

  return S;
}

cl_F DNuclear_5S(cl_F r,cl_F xi)
{
  cl_F S = ZERO;
    S = Power(r,-TWO) - (cl_float(14175.0,precision) + cl_float(28350.0,precision)*r*xi + cl_float(28350.0,precision)*Power(r,TWO)*Power(xi,TWO) + 

       cl_float(18900.0,precision)*Power(r,THREE)*Power(xi,THREE) + cl_float(9450.0,precision)*Power(r,FOUR)*Power(xi,FOUR) + 

       cl_float(3780.0,precision)*Power(r,FIVE)*Power(xi,FIVE) + cl_float(1260.0,precision)*Power(r,SIX)*Power(xi,SIX) + 

       cl_float(360.0,precision)*Power(r,SEVEN)*Power(xi,SEVEN) + cl_float(90.0,precision)*Power(r,EIGHT)*Power(xi,EIGHT) + 

       cl_float(20.0,precision)*Power(r,NINE)*Power(xi,NINE) + FOUR*Power(r,cl_float(10.0,precision))*Power(xi,cl_float(10.0,precision)))/

     (cl_float(14175.0,precision)*exp(TWO*r*xi)*Power(r,TWO))

    ;

  return S;
}

typedef cl_F t_slater_SS_func(cl_F r,cl_F xi,cl_F xj);
typedef cl_F t_slater_NS_func(cl_F r,cl_F xi);
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
  cl_F cr,cxi,cxj,cS;

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
  cl_F cr,cxi,cxj,cS;

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
  cl_F cr,cxi,cxj,cS;

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
  cl_F cr,cxi,cxj,cS;

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

