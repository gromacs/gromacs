/* slater_S_integrals.c (c) 2008 Paul J. van Maaren and David van der Spoel */
#include <stdio.h>
#include <math.h>
#include "slater_S_integrals.h"

#define Power(x, y)     (pow((double)(x), (double)(y)))
#define Pi              3.14159265358979323846264
#define E               2.71828182845904523536029

static double Slater_1S_1S(double rij,double xii,double xij)
{
  double S;

  if (xii == xij) {
    S =     (-24 + 24*Power(E,2*rij*xii) - 33*rij*xii - 18*Power(rij,2)*Power(xii,2) - 
        4*Power(rij,3)*Power(xii,3))/(24.*Power(E,2*rij*xii)*rij)
    ;
  }
  else {
    S =     (Power(E,2*rij*(xii + xij))*Power(Power(xii,2) - Power(xij,2),3) + 
        Power(E,2*rij*xij)*Power(xij,4)*
         (-3*Power(xii,2) - rij*Power(xii,3) + Power(xij,2) + 
           rij*xii*Power(xij,2)) - 
        Power(E,2*rij*xii)*Power(xii,4)*
         (Power(xii,2)*(1 + rij*xij) - Power(xij,2)*(3 + rij*xij)))/
      (Power(E,2*rij*(xii + xij))*rij*Power(Power(xii,2) - Power(xij,2),3))
    ;
  }
  return S;
}

static double Slater_2S_1S(double rij,double xii,double xij)
{
  double S;

  if (xii == xij) {
    S =     (-240 + 240*Power(E,2*rij*xii) - 375*rij*xii - 
        270*Power(rij,2)*Power(xii,2) - 115*Power(rij,3)*Power(xii,3) - 
        30*Power(rij,4)*Power(xii,4) - 4*Power(rij,5)*Power(xii,5))/
      (240.*Power(E,2*rij*xii)*rij)
    ;
  }
  else {
    S =     (6*Power(E,2*rij*(xii + xij))*Power(Power(xii,2) - Power(xij,2),5) - 
        6*Power(E,2*rij*xii)*Power(xii,6)*
         (-5*Power(xii,2)*Power(xij,2) + Power(xii,4)*(1 + rij*xij) - 
           Power(xij,4)*(4 + rij*xij)) + 
        Power(E,2*rij*xij)*Power(xij,4)*
         (-18*Power(rij,2)*Power(xii,8) - 2*Power(rij,3)*Power(xii,9) + 
           6*Power(xij,6) + 9*rij*xii*Power(xij,6) + 
           6*Power(xii,2)*Power(xij,4)*(-5 + Power(rij,2)*Power(xij,2)) + 
           42*Power(xii,6)*(-2 + Power(rij,2)*Power(xij,2)) - 
           30*Power(xii,4)*Power(xij,2)*(-2 + Power(rij,2)*Power(xij,2)) + 
           rij*Power(xii,3)*Power(xij,4)*(-45 + 2*Power(rij,2)*Power(xij,2)) + 
           Power(xii,7)*(-63*rij + 6*Power(rij,3)*Power(xij,2)) + 
           Power(xii,5)*(99*rij*Power(xij,2) - 6*Power(rij,3)*Power(xij,4))))/
      (6.*Power(E,2*rij*(xii + xij))*rij*Power(Power(xii,2) - Power(xij,2),5))
    ;
  }
  return S;
}

static double Slater_3S_1S(double rij,double xii,double xij)
{
  double S;

  if (xii == xij) {
    S =     (-120960 + 120960*Power(E,2*rij*xii) - 203175*rij*xii - 
        164430*Power(rij,2)*Power(xii,2) - 84420*Power(rij,3)*Power(xii,3) - 
        30240*Power(rij,4)*Power(xii,4) - 7728*Power(rij,5)*Power(xii,5) - 
        1344*Power(rij,6)*Power(xii,6) - 128*Power(rij,7)*Power(xii,7))/
      (120960.*Power(E,2*rij*xii)*rij)
    ;
  }
  else {
    S =     (45*Power(E,2*rij*(xii + xij))*Power(Power(xii,2) - Power(xij,2),7) - 
        15*Power(E,2*rij*xii)*Power(xii,8)*
         (7*Power(xii,4)*Power(xij,2)*(-3 + rij*xij) + 
           3*Power(xii,6)*(1 + rij*xij) - 3*Power(xij,6)*(5 + rij*xij) - 
           7*Power(xii,2)*Power(xij,4)*(9 + rij*xij)) + 
        Power(E,2*rij*xij)*Power(xij,4)*
         (-30*Power(rij,4)*Power(xii,14) - 2*Power(rij,5)*Power(xii,15) + 
           45*Power(xij,10) + 75*rij*xii*Power(xij,10) + 
           10*Power(rij,3)*Power(xii,13)*(-21 + Power(rij,2)*Power(xij,2)) + 
           15*rij*Power(xii,3)*Power(xij,8)*
            (-35 + 2*Power(rij,2)*Power(xij,2)) + 
           15*Power(xii,2)*Power(xij,8)*(-21 + 4*Power(rij,2)*Power(xij,2)) + 
           10*Power(rij,2)*Power(xii,12)*(-84 + 13*Power(rij,2)*Power(xij,2)) + 
           rij*Power(xii,5)*Power(xij,6)*
            (1575 - 210*Power(rij,2)*Power(xij,2) + 2*Power(rij,4)*Power(xij,4)) \
    - 5*rij*Power(xii,7)*Power(xij,4)*
            (513 - 132*Power(rij,2)*Power(xij,2) + 2*Power(rij,4)*Power(xij,4)) \
    + 10*rij*Power(xii,9)*Power(xij,2)*
            (333 - 102*Power(rij,2)*Power(xij,2) + 2*Power(rij,4)*Power(xij,4)) \
    + 5*Power(xii,4)*Power(xij,6)*(189 - 84*Power(rij,2)*Power(xij,2) + 
              2*Power(rij,4)*Power(xij,4)) - 
           10*rij*Power(xii,11)*(189 - 75*Power(rij,2)*Power(xij,2) + 
              2*Power(rij,4)*Power(xij,4)) - 
           35*Power(xii,6)*Power(xij,4)*
            (45 - 36*Power(rij,2)*Power(xij,2) + 2*Power(rij,4)*Power(xij,4)) + 
           90*Power(xii,8)*Power(xij,2)*
            (15 - 26*Power(rij,2)*Power(xij,2) + 2*Power(rij,4)*Power(xij,4)) - 
           10*Power(xii,10)*(189 - 228*Power(rij,2)*Power(xij,2) + 
              22*Power(rij,4)*Power(xij,4))))/
      (45.*Power(E,2*rij*(xii + xij))*rij*Power(Power(xii,2) - Power(xij,2),7))
    ;
  }
  return S;
}

static double Slater_4S_1S(double rij,double xii,double xij)
{
  double S;

  if (xii == xij) {
    S =     (-2903040 + 2903040*Power(E,2*rij*xii) - 5088825*rij*xii - 
        4371570*Power(rij,2)*Power(xii,2) - 2439990*Power(rij,3)*Power(xii,3) - 
        986580*Power(rij,4)*Power(xii,4) - 303912*Power(rij,5)*Power(xii,5) - 
        72576*Power(rij,6)*Power(xii,6) - 13248*Power(rij,7)*Power(xii,7) - 
        1728*Power(rij,8)*Power(xii,8) - 128*Power(rij,9)*Power(xii,9))/
      (2.90304e6*Power(E,2*rij*xii)*rij)
    ;
  }
  else {
    S =     (1260*Power(E,2*rij*(xii + xij))*Power(Power(xii,2) - Power(xij,2),9) - 
        1260*Power(E,2*rij*xii)*Power(xii,10)*
         (-63*Power(xii,4)*Power(xij,4) + Power(xii,8)*(1 + rij*xij) - 
           Power(xij,8)*(6 + rij*xij) + 
           3*Power(xii,6)*Power(xij,2)*(-3 + 2*rij*xij) - 
           3*Power(xii,2)*Power(xij,6)*(17 + 2*rij*xij)) + 
        Power(E,2*rij*xij)*Power(xij,4)*
         (-84*Power(rij,6)*Power(xii,20) - 4*Power(rij,7)*Power(xii,21) + 
           1260*Power(xij,14) + 2205*rij*xii*Power(xij,14) + 
           1890*Power(xii,2)*Power(xij,12)*(-6 + Power(rij,2)*Power(xij,2)) + 
           14*Power(rij,5)*Power(xii,19)*(-63 + 2*Power(rij,2)*Power(xij,2)) + 
           105*rij*Power(xii,3)*Power(xij,12)*
            (-189 + 10*Power(rij,2)*Power(xij,2)) + 
           28*Power(rij,4)*Power(xii,18)*(-210 + 19*Power(rij,2)*Power(xij,2)) + 
           126*rij*Power(xii,5)*Power(xij,10)*
            (630 - 75*Power(rij,2)*Power(xij,2) + Power(rij,4)*Power(xij,4)) - 
           42*Power(rij,3)*Power(xii,17)*
            (630 - 117*Power(rij,2)*Power(xij,2) + 2*Power(rij,4)*Power(xij,4)) \
    + 210*Power(xii,4)*Power(xij,10)*
            (216 - 81*Power(rij,2)*Power(xij,2) + 2*Power(rij,4)*Power(xij,4)) - 
           84*Power(rij,2)*Power(xii,16)*
            (945 - 330*Power(rij,2)*Power(xij,2) + 17*Power(rij,4)*Power(xij,4)) \
    + 28*Power(xii,6)*Power(xij,8)*
            (-3780 + 2430*Power(rij,2)*Power(xij,2) - 
              135*Power(rij,4)*Power(xij,4) + Power(rij,6)*Power(xij,6)) - 
           252*Power(xii,8)*Power(xij,6)*
            (-630 + 630*Power(rij,2)*Power(xij,2) - 
              60*Power(rij,4)*Power(xij,4) + Power(rij,6)*Power(xij,6)) + 
           2*rij*Power(xii,7)*Power(xij,8)*
            (-92610 + 18900*Power(rij,2)*Power(xij,2) - 
              567*Power(rij,4)*Power(xij,4) + 2*Power(rij,6)*Power(xij,6)) - 
           70*rij*Power(xii,13)*Power(xij,2)*
            (-3267 + 2223*Power(rij,2)*Power(xij,2) - 
              207*Power(rij,4)*Power(xij,4) + 2*Power(rij,6)*Power(xij,6)) + 
           462*Power(xii,10)*Power(xij,4)*
            (-360 + 495*Power(rij,2)*Power(xij,2) - 
              80*Power(rij,4)*Power(xij,4) + 2*Power(rij,6)*Power(xij,6)) - 
           7*rij*Power(xii,9)*Power(xij,6)*
            (-39915 + 12480*Power(rij,2)*Power(xij,2) - 
              666*Power(rij,4)*Power(xij,4) + 4*Power(rij,6)*Power(xij,6)) + 
           21*rij*Power(xii,11)*Power(xij,4)*
            (-11385 + 6690*Power(rij,2)*Power(xij,2) - 
              510*Power(rij,4)*Power(xij,4) + 4*Power(rij,6)*Power(xij,6)) + 
           420*Power(xii,14)*(-297 + 513*Power(rij,2)*Power(xij,2) - 
              129*Power(rij,4)*Power(xij,4) + 5*Power(rij,6)*Power(xij,6)) + 
           14*rij*Power(xii,15)*(-10395 + 7110*Power(rij,2)*Power(xij,2) - 
              819*Power(rij,4)*Power(xij,4) + 10*Power(rij,6)*Power(xij,6)) - 
           70*Power(xii,12)*Power(xij,2)*
            (-594 + 3699*Power(rij,2)*Power(xij,2) - 
              822*Power(rij,4)*Power(xij,4) + 26*Power(rij,6)*Power(xij,6))))/
      (1260.*Power(E,2*rij*(xii + xij))*rij*Power(Power(xii,2) - Power(xij,2),9))
    ;
  }
  return S;
}

static double Slater_2S_2S(double rij,double xii,double xij)
{
  double S;

  if (xii == xij) {
    S =     (-80640 + 80640*Power(E,2*rij*xii) - 131985*rij*xii - 
        102690*Power(rij,2)*Power(xii,2) - 49980*Power(rij,3)*Power(xii,3) - 
        16800*Power(rij,4)*Power(xii,4) - 4032*Power(rij,5)*Power(xii,5) - 
        672*Power(rij,6)*Power(xii,6) - 64*Power(rij,7)*Power(xii,7))/
      (80640.*Power(E,2*rij*xii)*rij)
    ;
  }
  else {
    S =     (6*Power(E,2*rij*(xii + xij))*Power(Power(xii,2) - Power(xij,2),7) - 
        Power(E,2*rij*xii)*Power(xii,6)*
         (21*Power(xii,4)*Power(xij,4)*
            (6 + 11*rij*xij + 2*Power(rij,2)*Power(xij,2)) - 
           2*Power(xij,8)*(90 + 54*rij*xij + 12*Power(rij,2)*Power(xij,2) + 
              Power(rij,3)*Power(xij,3)) + 
           Power(xii,8)*(6 + 9*rij*xij + 6*Power(rij,2)*Power(xij,2) + 
              2*Power(rij,3)*Power(xij,3)) + 
           Power(xii,2)*Power(xij,6)*
            (-390 - 69*rij*xij + 18*Power(rij,2)*Power(xij,2) + 
              4*Power(rij,3)*Power(xij,3)) - 
           Power(xii,6)*Power(xij,2)*
            (42 + 63*rij*xij + 42*Power(rij,2)*Power(xij,2) + 
              4*Power(rij,3)*Power(xij,3))) + 
        Power(E,2*rij*xij)*Power(xij,6)*
         (-24*Power(rij,2)*Power(xii,10) - 2*Power(rij,3)*Power(xii,11) - 
           69*rij*Power(xii,7)*Power(xij,2) + 6*Power(xij,8) + 
           9*rij*xii*Power(xij,8) + 
           4*rij*Power(xii,9)*(-27 + Power(rij,2)*Power(xij,2)) + 
           18*Power(xii,8)*(-10 + Power(rij,2)*Power(xij,2)) + 
           6*Power(xii,2)*Power(xij,6)*(-7 + Power(rij,2)*Power(xij,2)) - 
           42*Power(xii,4)*Power(xij,4)*(-3 + Power(rij,2)*Power(xij,2)) + 
           rij*Power(xii,3)*Power(xij,6)*(-63 + 2*Power(rij,2)*Power(xij,2)) + 
           6*Power(xii,6)*Power(xij,2)*(-65 + 7*Power(rij,2)*Power(xij,2)) + 
           Power(xii,5)*(231*rij*Power(xij,4) - 4*Power(rij,3)*Power(xij,6))))/
      (6.*Power(E,2*rij*(xii + xij))*rij*Power(Power(xii,2) - Power(xij,2),7))
    ;
  }
  return S;
}

static double Slater_3S_2S(double rij,double xii,double xij)
{
  double S;

  if (xii == xij) {
    S =     (-4354560 + 4354560*Power(E,2*rij*xii) - 7430535*rij*xii - 
        6151950*Power(rij,2)*Power(xii,2) - 3275370*Power(rij,3)*Power(xii,3) - 
        1251180*Power(rij,4)*Power(xii,4) - 361368*Power(rij,5)*Power(xii,5) - 
        80640*Power(rij,6)*Power(xii,6) - 13824*Power(rij,7)*Power(xii,7) - 
        1728*Power(rij,8)*Power(xii,8) - 128*Power(rij,9)*Power(xii,9))/
      (4.35456e6*Power(E,2*rij*xii)*rij)
    ;
  }
  else {
    S =     (90*Power(E,2*rij*(xii + xij))*Power(Power(xii,2) - Power(xij,2),9) - 
        5*Power(E,2*rij*xii)*Power(xii,8)*
         (Power(xii,8)*Power(xij,2)*
            (-162 - 243*rij*xij - 162*Power(rij,2)*Power(xij,2) + 
              2*Power(rij,3)*Power(xij,3)) - 
           18*Power(xii,6)*Power(xij,4)*
            (-36 - 75*rij*xij + 6*Power(rij,2)*Power(xij,2) + 
              2*Power(rij,3)*Power(xij,3)) + 
           3*Power(xii,10)*(6 + 9*rij*xij + 6*Power(rij,2)*Power(xij,2) + 
              2*Power(rij,3)*Power(xij,3)) + 
           18*Power(xii,4)*Power(xij,6)*
            (-216 + 81*rij*xij + 30*Power(rij,2)*Power(xij,2) + 
              2*Power(rij,3)*Power(xij,3)) - 
           3*Power(xij,10)*(330 + 165*rij*xij + 30*Power(rij,2)*Power(xij,2) + 
              2*Power(rij,3)*Power(xij,3)) - 
           Power(xii,2)*Power(xij,8)*
            (6378 + 2097*rij*xij + 198*Power(rij,2)*Power(xij,2) + 
              2*Power(rij,3)*Power(xij,3))) + 
        2*Power(E,2*rij*xij)*Power(xij,6)*
         (-40*Power(rij,4)*Power(xii,16) - 2*Power(rij,5)*Power(xii,17) + 
           45*Power(xij,12) + 75*rij*xii*Power(xij,12) + 
           8*Power(rij,3)*Power(xii,15)*(-45 + Power(rij,2)*Power(xij,2)) + 
           15*rij*Power(xii,3)*Power(xij,10)*
            (-45 + 2*Power(rij,2)*Power(xij,2)) + 
           15*Power(xii,2)*Power(xij,10)*(-27 + 4*Power(rij,2)*Power(xij,2)) + 
           10*Power(rij,2)*Power(xii,14)*(-180 + 11*Power(rij,2)*Power(xij,2)) + 
           15*rij*Power(xii,11)*Power(xij,2)*
            (-605 + 58*Power(rij,2)*Power(xij,2)) + 
           2*rij*Power(xii,5)*Power(xij,8)*
            (1350 - 135*Power(rij,2)*Power(xij,2) + Power(rij,4)*Power(xij,4)) + 
           10*Power(xii,4)*Power(xij,8)*
            (162 - 54*Power(rij,2)*Power(xij,2) + Power(rij,4)*Power(xij,4)) - 
           10*rij*Power(xii,13)*(495 - 49*Power(rij,2)*Power(xij,2) + 
              Power(rij,4)*Power(xij,4)) - 
           90*Power(xii,6)*Power(xij,6)*
            (42 - 24*Power(rij,2)*Power(xij,2) + Power(rij,4)*Power(xij,4)) - 
           30*Power(xii,12)*(198 + 10*Power(rij,2)*Power(xij,2) + 
              Power(rij,4)*Power(xij,4)) + 
           5*rij*Power(xii,9)*Power(xij,4)*
            (3555 - 396*Power(rij,2)*Power(xij,2) + 2*Power(rij,4)*Power(xij,4)) \
    - 2*rij*Power(xii,7)*Power(xij,6)*
            (2925 - 610*Power(rij,2)*Power(xij,2) + 4*Power(rij,4)*Power(xij,4)) \
    - 15*Power(xii,10)*Power(xij,2)*
            (1441 - 484*Power(rij,2)*Power(xij,2) + 
              12*Power(rij,4)*Power(xij,4)) + 
           5*Power(xii,8)*Power(xij,4)*
            (639 - 1368*Power(rij,2)*Power(xij,2) + 44*Power(rij,4)*Power(xij,4))\
    ))/(90.*Power(E,2*rij*(xii + xij))*rij*Power(Power(xii,2) - Power(xij,2),9))
    ;
  }
  return S;
}

static double Slater_4S_2S(double rij,double xii,double xij)
{
  double S;

  if (xii == xij) {
    S =     (-638668800 + 638668800*Power(E,2*rij*xii) - 1125310725.0*rij*xii - 
        973283850*Power(rij,2)*Power(xii,2) - 
        549063900*Power(rij,3)*Power(xii,3) - 
        226195200*Power(rij,4)*Power(xii,4) - 
        72099720*Power(rij,5)*Power(xii,5) - 
        18350640*Power(rij,6)*Power(xii,6) - 3785760*Power(rij,7)*Power(xii,7) - 
        633600*Power(rij,8)*Power(xii,8) - 84480*Power(rij,9)*Power(xii,9) - 
        8448*Power(rij,10)*Power(xii,10) - 512*Power(rij,11)*Power(xii,11))/
      (6.386688e8*Power(E,2*rij*xii)*rij)
    ;
  }
  else {
    S =     (1260*Power(E,2*rij*(xii + xij))*Power(Power(xii,2) - Power(xij,2),11) - 
        210*Power(E,2*rij*xii)*Power(xii,10)*
         (198*Power(xii,6)*Power(xij,6)*
            (-18 + 15*rij*xij + 2*Power(rij,2)*Power(xij,2)) - 
           22*Power(xii,8)*Power(xij,4)*
            (-15 - 36*rij*xij + 12*Power(rij,2)*Power(xij,2) + 
              Power(rij,3)*Power(xij,3)) - 
           2*Power(xij,12)*(273 + 117*rij*xij + 18*Power(rij,2)*Power(xij,2) + 
              Power(rij,3)*Power(xij,3)) + 
           Power(xii,12)*(6 + 9*rij*xij + 6*Power(rij,2)*Power(xij,2) + 
              2*Power(rij,3)*Power(xij,3)) + 
           11*Power(xii,4)*Power(xij,8)*
            (-1146 - 117*rij*xij + 18*Power(rij,2)*Power(xij,2) + 
              2*Power(rij,3)*Power(xij,3)) + 
           Power(xii,10)*Power(xij,2)*
            (-66 - 99*rij*xij - 66*Power(rij,2)*Power(xij,2) + 
              8*Power(rij,3)*Power(xij,3)) - 
           Power(xii,2)*Power(xij,10)*
            (6594 + 2151*rij*xij + 234*Power(rij,2)*Power(xij,2) + 
              8*Power(rij,3)*Power(xij,3))) + 
        Power(E,2*rij*xij)*Power(xij,6)*
         (-112*Power(rij,6)*Power(xii,22) - 4*Power(rij,7)*Power(xii,23) + 
           1260*Power(xij,16) + 2205*rij*xii*Power(xij,16) + 
           24*Power(rij,5)*Power(xii,21)*(-63 + Power(rij,2)*Power(xij,2)) + 
           630*Power(xii,2)*Power(xij,14)*(-22 + 3*Power(rij,2)*Power(xij,2)) + 
           105*rij*Power(xii,3)*Power(xij,14)*
            (-231 + 10*Power(rij,2)*Power(xij,2)) + 
           28*Power(rij,4)*Power(xii,20)*(-450 + 19*Power(rij,2)*Power(xij,2)) + 
           210*Power(xii,4)*Power(xij,12)*
            (330 - 99*Power(rij,2)*Power(xij,2) + 2*Power(rij,4)*Power(xij,4)) - 
           14*Power(rij,3)*Power(xii,19)*
            (4950 - 363*Power(rij,2)*Power(xij,2) + 4*Power(rij,4)*Power(xij,4)) \
    + 21*rij*Power(xii,5)*Power(xij,12)*
            (5775 - 550*Power(rij,2)*Power(xij,2) + 6*Power(rij,4)*Power(xij,4)) \
    - 28*Power(rij,2)*Power(xii,18)*
            (8910 - 825*Power(rij,2)*Power(xij,2) + 
              29*Power(rij,4)*Power(xij,4)) - 
           21*rij*Power(xii,15)*Power(xij,2)*
            (92235 - 20650*Power(rij,2)*Power(xij,2) + 
              646*Power(rij,4)*Power(xij,4)) - 
           14*Power(xii,16)*(38610 + 29205*Power(rij,2)*Power(xij,2) - 
              3030*Power(rij,4)*Power(xij,4) + 2*Power(rij,6)*Power(xij,6)) + 
           14*Power(xii,6)*Power(xij,10)*
            (-14850 + 7425*Power(rij,2)*Power(xij,2) - 
              330*Power(rij,4)*Power(xij,4) + 2*Power(rij,6)*Power(xij,6)) - 
           154*Power(xii,8)*Power(xij,8)*
            (-2700 + 2025*Power(rij,2)*Power(xij,2) - 
              150*Power(rij,4)*Power(xij,4) + 2*Power(rij,6)*Power(xij,6)) + 
           rij*Power(xii,7)*Power(xij,10)*
            (-363825 + 57750*Power(rij,2)*Power(xij,2) - 
              1386*Power(rij,4)*Power(xij,4) + 4*Power(rij,6)*Power(xij,6)) + 
           14*rij*Power(xii,17)*(-38610 + 825*Power(rij,2)*Power(xij,2) - 
              141*Power(rij,4)*Power(xij,4) + 4*Power(rij,6)*Power(xij,6)) - 
           7*rij*Power(xii,13)*Power(xij,4)*
            (-343395 + 100950*Power(rij,2)*Power(xij,2) - 
              3630*Power(rij,4)*Power(xij,4) + 8*Power(rij,6)*Power(xij,6)) + 
           7*rij*Power(xii,11)*Power(xij,6)*
            (-58095 + 64350*Power(rij,2)*Power(xij,2) - 
              2814*Power(rij,4)*Power(xij,4) + 8*Power(rij,6)*Power(xij,6)) - 
           3*rij*Power(xii,9)*Power(xij,8)*
            (-248325 + 55650*Power(rij,2)*Power(xij,2) - 
              2506*Power(rij,4)*Power(xij,4) + 8*Power(rij,6)*Power(xij,6)) + 
           70*Power(xii,14)*Power(xij,2)*
            (-43758 + 24291*Power(rij,2)*Power(xij,2) - 
              2166*Power(rij,4)*Power(xij,4) + 22*Power(rij,6)*Power(xij,6)) + 
           14*Power(xii,10)*Power(xij,6)*
            (-49770 + 37125*Power(rij,2)*Power(xij,2) - 
              5850*Power(rij,4)*Power(xij,4) + 82*Power(rij,6)*Power(xij,6)) - 
           14*Power(xii,12)*Power(xij,4)*
            (57330 + 95355*Power(rij,2)*Power(xij,2) - 
              11550*Power(rij,4)*Power(xij,4) + 142*Power(rij,6)*Power(xij,6))))/
      (1260.*Power(E,2*rij*(xii + xij))*rij*Power(Power(xii,2) - Power(xij,2),11))
    ;
  }
  return S;
}

double Slater_1S_2S(double rij,double xii,double xij)
{
  return Slater_2S_1S(rij,xij,xii);
}

static double Slater_3S_3S(double rij,double xii,double xij)
{
  double S;

  if (xii == xij) {
    S =     (-1437004800.0 + 1437004800.0*Power(E,2*rij*xii) - 2503064025.0*rij*xii - 
        2132118450.0*Power(rij,2)*Power(xii,2) - 
        1180664100.0*Power(rij,3)*Power(xii,3) - 
        476506800*Power(rij,4)*Power(xii,4) - 
        148856400*Power(rij,5)*Power(xii,5) - 
        37255680*Power(rij,6)*Power(xii,6) - 7603200*Power(rij,7)*Power(xii,7) - 
        1267200*Power(rij,8)*Power(xii,8) - 168960*Power(rij,9)*Power(xii,9) - 
        16896*Power(rij,10)*Power(xii,10) - 1024*Power(rij,11)*Power(xii,11))/
      (1.4370048e9*Power(E,2*rij*xii)*rij)
    ;
  }
  else {
    S =     (135*Power(E,2*rij*(xii + xij))*Power(Power(xii,2) - Power(xij,2),11) + 
        Power(E,2*rij*xij)*Power(xij,8)*
         (-150*Power(rij,4)*Power(xii,18) - 6*Power(rij,5)*Power(xii,19) + 
           135*Power(xij,14) + 225*rij*xii*Power(xij,14) + 
           10*Power(rij,3)*Power(xii,17)*(-165 + Power(rij,2)*Power(xij,2)) - 
           30*Power(rij,2)*Power(xii,16)*(330 + Power(rij,2)*Power(xij,2)) + 
           45*rij*Power(xii,3)*Power(xij,12)*
            (-55 + 2*Power(rij,2)*Power(xij,2)) + 
           45*Power(xii,2)*Power(xij,12)*(-33 + 4*Power(rij,2)*Power(xij,2)) + 
           rij*Power(xii,9)*Power(xij,6)*
            (234135 - 4950*Power(rij,2)*Power(xij,2) - 
              34*Power(rij,4)*Power(xij,4)) - 
           5*rij*Power(xii,7)*Power(xij,8)*
            (6237 - 1242*Power(rij,2)*Power(xij,2) + 
              2*Power(rij,4)*Power(xij,4)) + 
           3*rij*Power(xii,5)*Power(xij,10)*
            (4125 - 330*Power(rij,2)*Power(xij,2) + 
              2*Power(rij,4)*Power(xij,4)) + 
           15*Power(xii,4)*Power(xij,10)*
            (495 - 132*Power(rij,2)*Power(xij,2) + 2*Power(rij,4)*Power(xij,4)) \
    - 165*Power(xii,6)*Power(xij,8)*
            (135 - 60*Power(rij,2)*Power(xij,2) + 2*Power(rij,4)*Power(xij,4)) \
    - 5*rij*Power(xii,13)*Power(xij,2)*
            (43875 - 3438*Power(rij,2)*Power(xij,2) + 
              22*Power(rij,4)*Power(xij,4)) + 
           5*rij*Power(xii,11)*Power(xij,4)*
            (7695 - 2442*Power(rij,2)*Power(xij,2) + 
              22*Power(rij,4)*Power(xij,4)) + 
           15*Power(xii,8)*Power(xij,6)*
            (-33 - 3564*Power(rij,2)*Power(xij,2) + 
              26*Power(rij,4)*Power(xij,4)) + 
           rij*Power(xii,15)*(-32175 - 3690*Power(rij,2)*Power(xij,2) + 
              34*Power(rij,4)*Power(xij,4)) + 
           15*Power(xii,10)*Power(xij,4)*
            (-32277 + 1364*Power(rij,2)*Power(xij,2) + 
              66*Power(rij,4)*Power(xij,4)) + 
           15*Power(xii,14)*(-3003 - 2932*Power(rij,2)*Power(xij,2) + 
              94*Power(rij,4)*Power(xij,4)) - 
           15*Power(xii,12)*Power(xij,2)*
            (28119 - 5252*Power(rij,2)*Power(xij,2) + 
              154*Power(rij,4)*Power(xij,4))) + 
        Power(E,2*rij*xii)*Power(xii,8)*
         (-5*Power(xii,2)*Power(xij,12)*
            (-84357 - 43875*rij*xij - 8796*Power(rij,2)*Power(xij,2) - 
              738*Power(rij,3)*Power(xij,3) - 6*Power(rij,4)*Power(xij,4) + 
              2*Power(rij,5)*Power(xij,5)) - 
           3*Power(xii,14)*(45 + 75*rij*xij + 60*Power(rij,2)*Power(xij,2) + 
              30*Power(rij,3)*Power(xij,3) + 10*Power(rij,4)*Power(xij,4) + 
              2*Power(rij,5)*Power(xij,5)) - 
           55*Power(xii,8)*Power(xij,6)*
            (-405 - 567*rij*xij - 972*Power(rij,2)*Power(xij,2) - 
              90*Power(rij,3)*Power(xij,3) + 18*Power(rij,4)*Power(xij,4) + 
              2*Power(rij,5)*Power(xij,5)) + 
           55*Power(xii,6)*Power(xij,8)*
            (9 - 4257*rij*xij - 372*Power(rij,2)*Power(xij,2) + 
              222*Power(rij,3)*Power(xij,3) + 42*Power(rij,4)*Power(xij,4) + 
              2*Power(rij,5)*Power(xij,5)) + 
           3*Power(xij,14)*(15015 + 10725*rij*xij + 
              3300*Power(rij,2)*Power(xij,2) + 550*Power(rij,3)*Power(xij,3) + 
              50*Power(rij,4)*Power(xij,4) + 2*Power(rij,5)*Power(xij,5)) + 
           5*Power(xii,12)*Power(xij,2)*
            (297 + 495*rij*xij + 396*Power(rij,2)*Power(xij,2) + 
              198*Power(rij,3)*Power(xij,3) + 66*Power(rij,4)*Power(xij,4) + 
              2*Power(rij,5)*Power(xij,5)) + 
           Power(xii,10)*Power(xij,4)*
            (-7425 - 12375*rij*xij - 9900*Power(rij,2)*Power(xij,2) - 
              6210*Power(rij,3)*Power(xij,3) - 390*Power(rij,4)*Power(xij,4) + 
              34*Power(rij,5)*Power(xij,5)) - 
           Power(xii,4)*Power(xij,10)*
            (-484155 + 38475*rij*xij + 78780*Power(rij,2)*Power(xij,2) + 
              17190*Power(rij,3)*Power(xij,3) + 1410*Power(rij,4)*Power(xij,4) + 
              34*Power(rij,5)*Power(xij,5))))/
      (135.*Power(E,2*rij*(xii + xij))*rij*Power(Power(xii,2) - Power(xij,2),11))
    ;
  }
  return S;
}

static double Slater_4S_3S(double rij,double xii,double xij)
{
  double S;

  if (xii == xij) {
    S =     (-74724249600.0 + 74724249600.0*Power(E,2*rij*xii) - 132871488750.0*rij*xii - 
        116294478300.0*Power(rij,2)*Power(xii,2) - 
        66678987375.0*Power(rij,3)*Power(xii,3) - 
        28114836750.0*Power(rij,4)*Power(xii,4) - 
        9274044780.0*Power(rij,5)*Power(xii,5) - 
        2484321840.0*Power(rij,6)*Power(xii,6) - 
        553204080*Power(rij,7)*Power(xii,7) - 
        103783680*Power(rij,8)*Power(xii,8) - 
        16473600*Power(rij,9)*Power(xii,9) - 
        2196480*Power(rij,10)*Power(xii,10) - 
        239616*Power(rij,11)*Power(xii,11) - 19968*Power(rij,12)*Power(xii,12) - 
        1024*Power(rij,13)*Power(xii,13))/(7.47242496e10*Power(E,2*rij*xii)*rij)
    ;
  }
  else {
    S =     (3780*Power(E,2*rij*(xii + xij))*Power(Power(xii,2) - Power(xij,2),13) - 
        84*Power(E,2*rij*xii)*Power(xii,10)*
         (-715*Power(xii,8)*Power(xij,8)*
            (45 - 408*rij*xij + 72*Power(rij,2)*Power(xij,2) + 
              30*Power(rij,3)*Power(xij,3) + 2*Power(rij,4)*Power(xij,4)) + 
           Power(xii,16)*(45 + 75*rij*xij + 60*Power(rij,2)*Power(xij,2) + 
              30*Power(rij,3)*Power(xij,3) + 10*Power(rij,4)*Power(xij,4) + 
              2*Power(rij,5)*Power(xij,5)) - 
           Power(xij,16)*(32760 + 20475*rij*xij + 
              5460*Power(rij,2)*Power(xij,2) + 780*Power(rij,3)*Power(xij,3) + 
              60*Power(rij,4)*Power(xij,4) + 2*Power(rij,5)*Power(xij,5)) + 
           Power(xii,14)*Power(xij,2)*
            (-585 - 975*rij*xij - 780*Power(rij,2)*Power(xij,2) - 
              390*Power(rij,3)*Power(xij,3) - 130*Power(rij,4)*Power(xij,4) + 
              4*Power(rij,5)*Power(xij,5)) - 
           13*Power(xii,6)*Power(xij,10)*
            (68625 - 30735*rij*xij - 11580*Power(rij,2)*Power(xij,2) - 
              870*Power(rij,3)*Power(xij,3) + 30*Power(rij,4)*Power(xij,4) + 
              4*Power(rij,5)*Power(xij,5)) + 
           13*Power(xii,10)*Power(xij,6)*
            (-990 - 1155*rij*xij - 3300*Power(rij,2)*Power(xij,2) + 
              330*Power(rij,3)*Power(xij,3) + 110*Power(rij,4)*Power(xij,4) + 
              4*Power(rij,5)*Power(xij,5)) - 
           Power(xii,2)*Power(xij,14)*
            (538020 + 269325*rij*xij + 55020*Power(rij,2)*Power(xij,2) + 
              5610*Power(rij,3)*Power(xij,3) + 270*Power(rij,4)*Power(xij,4) + 
              4*Power(rij,5)*Power(xij,5)) - 
           6*Power(xii,12)*Power(xij,4)*
            (-585 - 975*rij*xij - 780*Power(rij,2)*Power(xij,2) - 
              555*Power(rij,3)*Power(xij,3) + 35*Power(rij,4)*Power(xij,4) + 
              6*Power(rij,5)*Power(xij,5)) + 
           6*Power(xii,4)*Power(xij,12)*
            (-256050 - 65235*rij*xij + 60*Power(rij,2)*Power(xij,2) + 
              1545*Power(rij,3)*Power(xij,3) + 175*Power(rij,4)*Power(xij,4) + 
              6*Power(rij,5)*Power(xij,5))) + 
        Power(E,2*rij*xij)*Power(xij,8)*
         (-420*Power(rij,6)*Power(xii,24) - 12*Power(rij,7)*Power(xii,25) + 
           3780*Power(xij,18) + 6615*rij*xii*Power(xij,18) + 
           22*Power(rij,5)*Power(xii,23)*(-315 + 2*Power(rij,2)*Power(xij,2)) + 
           252*Power(rij,4)*Power(xii,22)*(-275 + 3*Power(rij,2)*Power(xij,2)) + 
           1890*Power(xii,2)*Power(xij,16)*(-26 + 3*Power(rij,2)*Power(xij,2)) + 
           315*rij*Power(xii,3)*Power(xij,16)*
            (-273 + 10*Power(rij,2)*Power(xij,2)) + 
           630*Power(xii,4)*Power(xij,14)*
            (468 - 117*Power(rij,2)*Power(xij,2) + 2*Power(rij,4)*Power(xij,4)) \
    + 126*rij*Power(xii,5)*Power(xij,14)*
            (4095 - 325*Power(rij,2)*Power(xij,2) + 3*Power(rij,4)*Power(xij,4)) \
    + 2*Power(rij,3)*Power(xii,21)*
            (-225225 - 819*Power(rij,2)*Power(xij,2) + 
              8*Power(rij,4)*Power(xij,4)) + 
           42*Power(rij,2)*Power(xii,20)*
            (-45045 - 4030*Power(rij,2)*Power(xij,2) + 
              88*Power(rij,4)*Power(xij,4)) + 
           4*rij*Power(xii,9)*Power(xij,10)*
            (1261260 - 204435*Power(rij,2)*Power(xij,2) + 
              8694*Power(rij,4)*Power(xij,4) - 11*Power(rij,6)*Power(xij,6)) - 
           336*Power(xii,12)*Power(xij,6)*
            (115200 + 63180*Power(rij,2)*Power(xij,2) - 
              2860*Power(rij,4)*Power(xij,4) + Power(rij,6)*Power(xij,6)) + 
           84*Power(xii,6)*Power(xij,12)*
            (-12870 + 5265*Power(rij,2)*Power(xij,2) - 
              195*Power(rij,4)*Power(xij,4) + Power(rij,6)*Power(xij,6)) - 
           1092*Power(xii,8)*Power(xij,10)*
            (-2475 + 1485*Power(rij,2)*Power(xij,2) - 
              90*Power(rij,4)*Power(xij,4) + Power(rij,6)*Power(xij,6)) + 
           6*rij*Power(xii,7)*Power(xij,12)*
            (-315315 + 40950*Power(rij,2)*Power(xij,2) - 
              819*Power(rij,4)*Power(xij,4) + 2*Power(rij,6)*Power(xij,6)) - 
           8*rij*Power(xii,11)*Power(xij,8)*
            (-831285 - 540540*Power(rij,2)*Power(xij,2) + 
              9639*Power(rij,4)*Power(xij,4) + 2*Power(rij,6)*Power(xij,6)) + 
           42*rij*Power(xii,13)*Power(xij,6)*
            (1640835 - 153660*Power(rij,2)*Power(xij,2) + 
              390*Power(rij,4)*Power(xij,4) + 8*Power(rij,6)*Power(xij,6)) + 
           252*Power(xii,10)*Power(xij,8)*
            (-30225 + 8580*Power(rij,2)*Power(xij,2) - 
              2090*Power(rij,4)*Power(xij,4) + 12*Power(rij,6)*Power(xij,6)) - 
           21*rij*Power(xii,19)*(225225 + 103350*Power(rij,2)*Power(xij,2) - 
              4584*Power(rij,4)*Power(xij,4) + 16*Power(rij,6)*Power(xij,6)) - 
           420*Power(xii,14)*Power(xij,4)*
            (322704 - 54603*Power(rij,2)*Power(xij,2) + 
              260*Power(rij,4)*Power(xij,4) + 26*Power(rij,6)*Power(xij,6)) - 
           14*rij*Power(xii,15)*Power(xij,4)*
            (1911105 + 62190*Power(rij,2)*Power(xij,2) - 
              10998*Power(rij,4)*Power(xij,4) + 52*Power(rij,6)*Power(xij,6)) + 
           252*Power(xii,16)*Power(xij,2)*
            (-278070 + 52395*Power(rij,2)*Power(xij,2) - 
              5030*Power(rij,4)*Power(xij,4) + 78*Power(rij,6)*Power(xij,6)) + 
           7*rij*Power(xii,17)*Power(xij,2)*
            (-6809985 + 890100*Power(rij,2)*Power(xij,2) - 
              30168*Power(rij,4)*Power(xij,4) + 104*Power(rij,6)*Power(xij,6)) - 
           42*Power(xii,18)*(128700 + 331695*Power(rij,2)*Power(xij,2) - 
              26140*Power(rij,4)*Power(xij,4) + 344*Power(rij,6)*Power(xij,6))))/
      (3780.*Power(E,2*rij*(xii + xij))*rij*Power(Power(xii,2) - Power(xij,2),13))
    ;
  }
  return S;
}

double Slater_1S_3S(double rij,double xii,double xij)
{
  return Slater_3S_1S(rij,xij,xii);
}

double Slater_2S_3S(double rij,double xii,double xij)
{
  return Slater_3S_2S(rij,xij,xii);
}

static double Slater_4S_4S(double rij,double xii,double xij)
{
  double S;

  if (xii == xij) {
    S =     (-83691159552000.0 + 83691159552000.0*Power(E,2*rij*xii) - 
        150568359566625.0*rij*xii - 133754400029250.0*Power(rij,2)*Power(xii,2) - 
        78142908343500.0*Power(rij,3)*Power(xii,3) - 
        33740723016000.0*Power(rij,4)*Power(xii,4) - 
        11470756096800.0*Power(rij,5)*Power(xii,5) - 
        3193358968800.0*Power(rij,6)*Power(xii,6) - 
        747112766400.0*Power(rij,7)*Power(xii,7) - 
        149448499200.0*Power(rij,8)*Power(xii,8) - 
        25830604800.0*Power(rij,9)*Power(xii,9) - 
        3874590720.0*Power(rij,10)*Power(xii,10) - 
        503193600*Power(rij,11)*Power(xii,11) - 
        55910400*Power(rij,12)*Power(xii,12) - 
        5160960*Power(rij,13)*Power(xii,13) - 
        368640*Power(rij,14)*Power(xii,14) - 16384*Power(rij,15)*Power(xii,15))/
      (8.3691159552e13*Power(E,2*rij*xii)*rij)
    ;
  }
  else {
    S =     (1260*Power(E,2*rij*(xii + xij))*Power(Power(xii,2) - Power(xij,2),15) + 
        Power(E,2*rij*xij)*Power(xij,10)*
         (-3276*Power(rij,5)*Power(xii,25) - 168*Power(rij,6)*Power(xii,26) - 
           4*Power(rij,7)*Power(xii,27) + 1260*Power(xij,20) + 
           2205*rij*xii*Power(xij,20) + 
           1890*Power(xii,2)*Power(xij,18)*(-10 + Power(rij,2)*Power(xij,2)) - 
           420*Power(rij,4)*Power(xii,24)*(91 + Power(rij,2)*Power(xij,2)) + 
           525*rij*Power(xii,3)*Power(xij,18)*
            (-63 + 2*Power(rij,2)*Power(xij,2)) + 
           42*Power(rij,3)*Power(xii,23)*
            (-6825 - 405*Power(rij,2)*Power(xij,2) + 
              2*Power(rij,4)*Power(xij,4)) + 
           63*rij*Power(xii,5)*Power(xij,16)*
            (3675 - 250*Power(rij,2)*Power(xij,2) + 
              2*Power(rij,4)*Power(xij,4)) + 
           210*Power(xii,4)*Power(xij,16)*
            (630 - 135*Power(rij,2)*Power(xij,2) + 2*Power(rij,4)*Power(xij,4)) \
    + 252*Power(rij,2)*Power(xii,22)*
            (-5460 - 1225*Power(rij,2)*Power(xij,2) + 
              17*Power(rij,4)*Power(xij,4)) - 
           1260*rij*Power(xii,17)*Power(xij,4)*
            (141729 - 10145*Power(rij,2)*Power(xij,2) + 
              116*Power(rij,4)*Power(xij,4)) + 
           21*rij*Power(xii,9)*Power(xij,12)*
            (164775 - 18460*Power(rij,2)*Power(xij,2) + 
              828*Power(rij,4)*Power(xij,4)) + 
           14*Power(xii,6)*Power(xij,14)*
            (-40950 + 14175*Power(rij,2)*Power(xij,2) - 
              450*Power(rij,4)*Power(xij,4) + 2*Power(rij,6)*Power(xij,6)) - 
           210*Power(xii,8)*Power(xij,12)*
            (-8190 + 4095*Power(rij,2)*Power(xij,2) - 
              210*Power(rij,4)*Power(xij,4) + 2*Power(rij,6)*Power(xij,6)) + 
           42*Power(xii,10)*Power(xij,10)*
            (-209430 - 2925*Power(rij,2)*Power(xij,2) - 
              8840*Power(rij,4)*Power(xij,4) + 4*Power(rij,6)*Power(xij,6)) + 
           rij*Power(xii,7)*Power(xij,14)*
            (-1003275 + 110250*Power(rij,2)*Power(xij,2) - 
              1890*Power(rij,4)*Power(xij,4) + 4*Power(rij,6)*Power(xij,6)) - 
           21*rij*Power(xii,11)*Power(xij,10)*
            (-1033695 - 218400*Power(rij,2)*Power(xij,2) + 
              552*Power(rij,4)*Power(xij,4) + 4*Power(rij,6)*Power(xij,6)) + 
           280*Power(xii,18)*Power(xij,2)*
            (-385560 - 73953*Power(rij,2)*Power(xij,2) + 
              2370*Power(rij,4)*Power(xij,4) + 4*Power(rij,6)*Power(xij,6)) - 
           35*rij*Power(xii,15)*Power(xij,6)*
            (-1565613 + 359520*Power(rij,2)*Power(xij,2) - 
              7020*Power(rij,4)*Power(xij,4) + 8*Power(rij,6)*Power(xij,6)) + 
           14*rij*Power(xii,19)*Power(xij,2)*
            (-4980150 + 126765*Power(rij,2)*Power(xij,2) - 
              3852*Power(rij,4)*Power(xij,4) + 20*Power(rij,6)*Power(xij,6)) - 
           630*Power(xii,14)*Power(xij,6)*
            (708714 - 14385*Power(rij,2)*Power(xij,2) - 
              2340*Power(rij,4)*Power(xij,4) + 20*Power(rij,6)*Power(xij,6)) + 
           210*Power(xii,16)*Power(xij,4)*
            (-2087532 + 328491*Power(rij,2)*Power(xij,2) - 
              11740*Power(rij,4)*Power(xij,4) + 52*Power(rij,6)*Power(xij,6)) - 
           84*Power(xii,20)*(59670 + 236250*Power(rij,2)*Power(xij,2) - 
              8745*Power(rij,4)*Power(xij,4) + 92*Power(rij,6)*Power(xij,6)) - 
           2*rij*Power(xii,21)*(1949220 + 1598625*Power(rij,2)*Power(xij,2) - 
              41391*Power(rij,4)*Power(xij,4) + 128*Power(rij,6)*Power(xij,6)) \
    + rij*Power(xii,13)*Power(xij,8)*
            (173037375 - 2784600*Power(rij,2)*Power(xij,2) - 
              112140*Power(rij,4)*Power(xij,4) + 256*Power(rij,6)*Power(xij,6)) \
    + 14*Power(xii,12)*Power(xij,8)*
            (-7260750 - 2521935*Power(rij,2)*Power(xij,2) + 
              19500*Power(rij,4)*Power(xij,4) + 344*Power(rij,6)*Power(xij,6))) \
    + Power(E,2*rij*xii)*Power(xii,10)*
         (210*Power(xii,2)*Power(xij,18)*
            (514080 + 332010*rij*xij + 94500*Power(rij,2)*Power(xij,2) + 
              15225*Power(rij,3)*Power(xij,3) + 
              1470*Power(rij,4)*Power(xij,4) + 81*Power(rij,5)*Power(xij,5) + 
              2*Power(rij,6)*Power(xij,6)) + 
           105*Power(xii,18)*Power(xij,2)*
            (180 + 315*rij*xij + 270*Power(rij,2)*Power(xij,2) + 
              150*Power(rij,3)*Power(xij,3) + 60*Power(rij,4)*Power(xij,4) + 
              18*Power(rij,5)*Power(xij,5) + 4*Power(rij,6)*Power(xij,6)) - 
           1365*Power(xii,10)*Power(xij,10)*
            (-6444 + 15903*rij*xij - 25866*Power(rij,2)*Power(xij,2) - 
              2040*Power(rij,3)*Power(xij,3) + 1080*Power(rij,4)*Power(xij,4) + 
              180*Power(rij,5)*Power(xij,5) + 8*Power(rij,6)*Power(xij,6)) + 
           Power(xii,14)*Power(xij,6)*
            (573300 + 1003275*rij*xij + 859950*Power(rij,2)*Power(xij,2) + 
              387660*Power(rij,3)*Power(xij,3) + 
              371280*Power(rij,4)*Power(xij,4) + 
              11592*Power(rij,5)*Power(xij,5) - 
              4816*Power(rij,6)*Power(xij,6) - 256*Power(rij,7)*Power(xij,7)) + 
           2*Power(xij,20)*(2506140 + 1949220*rij*xij + 
              687960*Power(rij,2)*Power(xij,2) + 
              143325*Power(rij,3)*Power(xij,3) + 
              19110*Power(rij,4)*Power(xij,4) + 
              1638*Power(rij,5)*Power(xij,5) + 84*Power(rij,6)*Power(xij,6) + 
              2*Power(rij,7)*Power(xij,7)) - 
           42*Power(xii,4)*Power(xij,16)*
            (-10437660 - 4251870*rij*xij - 493020*Power(rij,2)*Power(xij,2) + 
              42255*Power(rij,3)*Power(xij,3) + 
              17490*Power(rij,4)*Power(xij,4) + 
              1971*Power(rij,5)*Power(xij,5) + 102*Power(rij,6)*Power(xij,6) + 
              2*Power(rij,7)*Power(xij,7)) + 
           21*Power(xii,16)*Power(xij,4)*
            (-6300 - 11025*rij*xij - 9450*Power(rij,2)*Power(xij,2) - 
              5250*Power(rij,3)*Power(xij,3) - 2100*Power(rij,4)*Power(xij,4) - 
              828*Power(rij,5)*Power(xij,5) - 8*Power(rij,6)*Power(xij,6) + 
              4*Power(rij,7)*Power(xij,7)) - 
           Power(xii,20)*(1260 + 2205*rij*xij + 
              1890*Power(rij,2)*Power(xij,2) + 1050*Power(rij,3)*Power(xij,3) + 
              420*Power(rij,4)*Power(xij,4) + 126*Power(rij,5)*Power(xij,5) + 
              28*Power(rij,6)*Power(xij,6) + 4*Power(rij,7)*Power(xij,7)) - 
           35*Power(xii,8)*Power(xij,12)*
            (-2904300 + 4943925*rij*xij + 258930*Power(rij,2)*Power(xij,2) - 
              359520*Power(rij,3)*Power(xij,3) - 
              70440*Power(rij,4)*Power(xij,4) - 
              4176*Power(rij,5)*Power(xij,5) + 32*Power(rij,6)*Power(xij,6) + 
              8*Power(rij,7)*Power(xij,7)) + 
           35*Power(xii,12)*Power(xij,8)*
            (-49140 - 98865*rij*xij + 3510*Power(rij,2)*Power(xij,2) - 
              131040*Power(rij,3)*Power(xij,3) - 
              7800*Power(rij,4)*Power(xij,4) + 3204*Power(rij,5)*Power(xij,5) + 
              360*Power(rij,6)*Power(xij,6) + 8*Power(rij,7)*Power(xij,7)) + 
           Power(xii,6)*Power(xij,14)*
            (446489820 - 54796455*rij*xij - 68983110*Power(rij,2)*Power(xij,2) - 
              12782700*Power(rij,3)*Power(xij,3) - 
              663600*Power(rij,4)*Power(xij,4) + 
              53928*Power(rij,5)*Power(xij,5) + 7728*Power(rij,6)*Power(xij,6) + 
              256*Power(rij,7)*Power(xij,7))))/
      (1260.*Power(E,2*rij*(xii + xij))*rij*Power(Power(xii,2) - Power(xij,2),15))
    ;
  }
  return S;
}

double Slater_1S_4S(double rij,double xii,double xij)
{
  return Slater_4S_1S(rij,xij,xii);
}

double Slater_2S_4S(double rij,double xii,double xij)
{
  return Slater_4S_2S(rij,xij,xii);
}

double Slater_3S_4S(double rij,double xii,double xij)
{
  return Slater_4S_3S(rij,xij,xii);
}

static double DSlater_1S_1S(double rij,double xii,double xij)
{
  double S;

  if (xii == xij) {
    S =     (-6 + 6*Power(E,2*rij*xii) - 12*rij*xii - 12*Power(rij,2)*Power(xii,2) - 
        7*Power(rij,3)*Power(xii,3) - 2*Power(rij,4)*Power(xii,4))/
      (6.*Power(E,2*rij*xii)*Power(rij,2))
    ;
  }
  else {
    S =     (Power(E,2*rij*(xii + xij))*Power(Power(xii,2) - Power(xij,2),3) + 
        Power(E,2*rij*xij)*Power(xij,4)*
         (-6*rij*Power(xii,3) - 2*Power(rij,2)*Power(xii,4) + Power(xij,2) + 
           2*rij*xii*Power(xij,2) + 
           Power(xii,2)*(-3 + 2*Power(rij,2)*Power(xij,2))) - 
        Power(E,2*rij*xii)*Power(xii,4)*
         (Power(xii,2)*(1 + 2*rij*xij + 2*Power(rij,2)*Power(xij,2)) - 
           Power(xij,2)*(3 + 6*rij*xij + 2*Power(rij,2)*Power(xij,2))))/
      (Power(E,2*rij*(xii + xij))*Power(rij,2)*
        Power(Power(xii,2) - Power(xij,2),3))
    ;
  }
  return S;
}

static double DSlater_2S_1S(double rij,double xii,double xij)
{
  double S;

  if (xii == xij) {
    S =     (-120 + 120*Power(E,2*rij*xii) - 240*rij*xii - 
        240*Power(rij,2)*Power(xii,2) - 155*Power(rij,3)*Power(xii,3) - 
        70*Power(rij,4)*Power(xii,4) - 22*Power(rij,5)*Power(xii,5) - 
        4*Power(rij,6)*Power(xii,6))/(120.*Power(E,2*rij*xii)*Power(rij,2))
    ;
  }
  else {
    S =     (3*Power(E,2*rij*(xii + xij))*Power(Power(xii,2) - Power(xij,2),5) - 
        3*Power(E,2*rij*xii)*Power(xii,6)*
         (-5*Power(xii,2)*Power(xij,2)*(1 + 2*rij*xij) - 
           2*Power(xij,4)*(2 + 4*rij*xij + Power(rij,2)*Power(xij,2)) + 
           Power(xii,4)*(1 + 2*rij*xij + 2*Power(rij,2)*Power(xij,2))) + 
        Power(E,2*rij*xij)*Power(xij,4)*
         (-16*Power(rij,3)*Power(xii,9) - 2*Power(rij,4)*Power(xii,10) + 
           3*Power(xij,6) + 6*rij*xii*Power(xij,6) + 
           6*Power(rij,2)*Power(xii,8)*(-9 + Power(rij,2)*Power(xij,2)) + 
           3*Power(xii,2)*Power(xij,4)*(-5 + 2*Power(rij,2)*Power(xij,2)) + 
           12*rij*Power(xii,7)*(-7 + 3*Power(rij,2)*Power(xij,2)) + 
           Power(xii,5)*(60*rij*Power(xij,2) - 24*Power(rij,3)*Power(xij,4)) + 
           2*Power(xii,4)*Power(xij,2)*
            (15 - 15*Power(rij,2)*Power(xij,2) + Power(rij,4)*Power(xij,4)) - 
           6*Power(xii,6)*(7 - 13*Power(rij,2)*Power(xij,2) + 
              Power(rij,4)*Power(xij,4)) + 
           Power(xii,3)*(-30*rij*Power(xij,4) + 4*Power(rij,3)*Power(xij,6))))/
      (3.*Power(E,2*rij*(xii + xij))*Power(rij,2)*
        Power(Power(xii,2) - Power(xij,2),5))
    ;
  }
  return S;
}

static double DSlater_3S_1S(double rij,double xii,double xij)
{
  double S;

  if (xii == xij) {
    S =     (-30240 + 30240*Power(E,2*rij*xii) - 60480*rij*xii - 
        60480*Power(rij,2)*Power(xii,2) - 40005*Power(rij,3)*Power(xii,3) - 
        19530*Power(rij,4)*Power(xii,4) - 7392*Power(rij,5)*Power(xii,5) - 
        2184*Power(rij,6)*Power(xii,6) - 480*Power(rij,7)*Power(xii,7) - 
        64*Power(rij,8)*Power(xii,8))/(30240.*Power(E,2*rij*xii)*Power(rij,2))
    ;
  }
  else {
    S =     (45*Power(E,2*rij*(xii + xij))*Power(Power(xii,2) - Power(xij,2),7) - 
        15*Power(E,2*rij*xii)*Power(xii,8)*
         (7*Power(xii,4)*Power(xij,2)*
            (-3 - 6*rij*xij + 2*Power(rij,2)*Power(xij,2)) - 
           3*Power(xij,6)*(5 + 10*rij*xij + 2*Power(rij,2)*Power(xij,2)) - 
           7*Power(xii,2)*Power(xij,4)*
            (9 + 18*rij*xij + 2*Power(rij,2)*Power(xij,2)) + 
           Power(xii,6)*(3 + 6*rij*xij + 6*Power(rij,2)*Power(xij,2))) + 
        Power(E,2*rij*xij)*Power(xij,4)*
         (-52*Power(rij,5)*Power(xii,15) - 4*Power(rij,6)*Power(xii,16) + 
           45*Power(xij,10) + 90*rij*xii*Power(xij,10) + 
           10*Power(rij,4)*Power(xii,14)*(-33 + 2*Power(rij,2)*Power(xij,2)) + 
           30*rij*Power(xii,3)*Power(xij,8)*
            (-21 + 2*Power(rij,2)*Power(xij,2)) + 
           45*Power(xii,2)*Power(xij,8)*(-7 + 2*Power(rij,2)*Power(xij,2)) + 
           20*Power(rij,3)*Power(xii,13)*(-63 + 11*Power(rij,2)*Power(xij,2)) + 
           6*rij*Power(xii,5)*Power(xij,6)*
            (315 - 70*Power(rij,2)*Power(xij,2) + 2*Power(rij,4)*Power(xij,4)) + 
           15*Power(xii,4)*Power(xij,6)*
            (63 - 42*Power(rij,2)*Power(xij,2) + 2*Power(rij,4)*Power(xij,4)) - 
           50*rij*Power(xii,7)*Power(xij,4)*
            (63 - 24*Power(rij,2)*Power(xij,2) + 2*Power(rij,4)*Power(xij,4)) - 
           180*rij*Power(xii,11)*(21 - 17*Power(rij,2)*Power(xij,2) + 
              2*Power(rij,4)*Power(xij,4)) - 
           10*Power(rij,2)*Power(xii,12)*
            (294 - 111*Power(rij,2)*Power(xij,2) + 4*Power(rij,4)*Power(xij,4)) \
    + 20*rij*Power(xii,9)*Power(xij,2)*
            (135 - 132*Power(rij,2)*Power(xij,2) + 14*Power(rij,4)*Power(xij,4)) \
    - 10*Power(xii,8)*Power(xij,2)*
            (-135 + 279*Power(rij,2)*Power(xij,2) - 
              78*Power(rij,4)*Power(xij,4) + 2*Power(rij,6)*Power(xij,6)) + 
           Power(xii,6)*Power(xij,4)*
            (-1575 + 1890*Power(rij,2)*Power(xij,2) - 
              210*Power(rij,4)*Power(xij,4) + 4*Power(rij,6)*Power(xij,6)) + 
           10*Power(xii,10)*(-189 + 438*Power(rij,2)*Power(xij,2) - 
              138*Power(rij,4)*Power(xij,4) + 4*Power(rij,6)*Power(xij,6))))/
      (45.*Power(E,2*rij*(xii + xij))*Power(rij,2)*
        Power(Power(xii,2) - Power(xij,2),7))
    ;
  }
  return S;
}

static double DSlater_4S_1S(double rij,double xii,double xij)
{
  double S;

  if (xii == xij) {
    S =     (-362880 + 362880*Power(E,2*rij*xii) - 725760*rij*xii - 
        725760*Power(rij,2)*Power(xii,2) - 482895*Power(rij,3)*Power(xii,3) - 
        240030*Power(rij,4)*Power(xii,4) - 94689*Power(rij,5)*Power(xii,5) - 
        30618*Power(rij,6)*Power(xii,6) - 8208*Power(rij,7)*Power(xii,7) - 
        1800*Power(rij,8)*Power(xii,8) - 304*Power(rij,9)*Power(xii,9) - 
        32*Power(rij,10)*Power(xii,10))/(362880.*Power(E,2*rij*xii)*Power(rij,2))
    ;
  }
  else {
    S =     (315*Power(E,2*rij*(xii + xij))*Power(Power(xii,2) - Power(xij,2),9) - 
        315*Power(E,2*rij*xii)*Power(xii,10)*
         (-63*Power(xii,4)*Power(xij,4)*(1 + 2*rij*xij) - 
           2*Power(xij,8)*(3 + 6*rij*xij + Power(rij,2)*Power(xij,2)) + 
           Power(xii,8)*(1 + 2*rij*xij + 2*Power(rij,2)*Power(xij,2)) + 
           3*Power(xii,6)*Power(xij,2)*
            (-3 - 6*rij*xij + 4*Power(rij,2)*Power(xij,2)) - 
           3*Power(xii,2)*Power(xij,6)*
            (17 + 34*rij*xij + 4*Power(rij,2)*Power(xij,2))) + 
        Power(E,2*rij*xij)*Power(xij,4)*
         (-36*Power(rij,7)*Power(xii,21) - 2*Power(rij,8)*Power(xii,22) + 
           315*Power(xij,14) + 630*rij*xii*Power(xij,14) + 
           14*Power(rij,6)*Power(xii,20)*(-24 + Power(rij,2)*Power(xij,2)) + 
           210*rij*Power(xii,3)*Power(xij,12)*
            (-27 + 2*Power(rij,2)*Power(xij,2)) + 
           315*Power(xii,2)*Power(xij,12)*(-9 + 2*Power(rij,2)*Power(xij,2)) + 
           14*Power(rij,5)*Power(xii,19)*(-147 + 16*Power(rij,2)*Power(xij,2)) + 
           84*rij*Power(xii,5)*Power(xij,10)*
            (270 - 45*Power(rij,2)*Power(xij,2) + Power(rij,4)*Power(xij,4)) + 
           210*Power(xii,4)*Power(xij,10)*
            (54 - 27*Power(rij,2)*Power(xij,2) + Power(rij,4)*Power(xij,4)) - 
           14*Power(rij,4)*Power(xii,18)*
            (630 - 128*Power(rij,2)*Power(xij,2) + 3*Power(rij,4)*Power(xij,4)) \
    - 42*Power(rij,3)*Power(xii,17)*
            (630 - 213*Power(rij,2)*Power(xij,2) + 14*Power(rij,4)*Power(xij,4)) \
    + 4*rij*Power(xii,7)*Power(xij,8)*
            (-13230 + 3780*Power(rij,2)*Power(xij,2) - 
              189*Power(rij,4)*Power(xij,4) + 2*Power(rij,6)*Power(xij,6)) + 
           14*Power(xii,6)*Power(xij,8)*
            (-1890 + 1620*Power(rij,2)*Power(xij,2) - 
              135*Power(rij,4)*Power(xij,4) + 2*Power(rij,6)*Power(xij,6)) - 
           42*rij*Power(xii,9)*Power(xij,6)*
            (-1890 + 850*Power(rij,2)*Power(xij,2) - 
              69*Power(rij,4)*Power(xij,4) + 2*Power(rij,6)*Power(xij,6)) + 
           14*Power(rij,2)*Power(xii,16)*
            (-3780 + 2070*Power(rij,2)*Power(xij,2) - 
              282*Power(rij,4)*Power(xij,4) + 5*Power(rij,6)*Power(xij,6)) + 
           42*rij*Power(xii,11)*Power(xij,4)*
            (-1980 + 1050*Power(rij,2)*Power(xij,2) - 
              185*Power(rij,4)*Power(xij,4) + 8*Power(rij,6)*Power(xij,6)) - 
           70*rij*Power(xii,13)*Power(xij,2)*
            (-297 + 738*Power(rij,2)*Power(xij,2) - 
              204*Power(rij,4)*Power(xij,4) + 10*Power(rij,6)*Power(xij,6)) + 
           42*rij*Power(xii,15)*(-1485 + 1380*Power(rij,2)*Power(xij,2) - 
              372*Power(rij,4)*Power(xij,4) + 20*Power(rij,6)*Power(xij,6)) + 
           2*Power(xii,8)*Power(xij,6)*
            (19845 - 26460*Power(rij,2)*Power(xij,2) + 
              3780*Power(rij,4)*Power(xij,4) - 126*Power(rij,6)*Power(xij,6) + 
              Power(rij,8)*Power(xij,8)) - 
           14*Power(xii,10)*Power(xij,4)*
            (2970 - 5895*Power(rij,2)*Power(xij,2) + 
              1140*Power(rij,4)*Power(xij,4) - 84*Power(rij,6)*Power(xij,6) + 
              Power(rij,8)*Power(xij,8)) - 
           35*Power(xii,14)*(891 - 1728*Power(rij,2)*Power(xij,2) + 
              1062*Power(rij,4)*Power(xij,4) - 132*Power(rij,6)*Power(xij,6) + 
              2*Power(rij,8)*Power(xij,8)) + 
           7*Power(xii,12)*Power(xij,2)*
            (1485 - 7830*Power(rij,2)*Power(xij,2) + 
              3870*Power(rij,4)*Power(xij,4) - 440*Power(rij,6)*Power(xij,6) + 
              6*Power(rij,8)*Power(xij,8))))/
      (315.*Power(E,2*rij*(xii + xij))*Power(rij,2)*
        Power(Power(xii,2) - Power(xij,2),9))
    ;
  }
  return S;
}

static double DSlater_2S_2S(double rij,double xii,double xij)
{
  double S;

  if (xii == xij) {
    S =     (-20160 + 20160*Power(E,2*rij*xii) - 40320*rij*xii - 
        40320*Power(rij,2)*Power(xii,2) - 26355*Power(rij,3)*Power(xii,3) - 
        12390*Power(rij,4)*Power(xii,4) - 4368*Power(rij,5)*Power(xii,5) - 
        1176*Power(rij,6)*Power(xii,6) - 240*Power(rij,7)*Power(xii,7) - 
        32*Power(rij,8)*Power(xii,8))/(20160.*Power(E,2*rij*xii)*Power(rij,2))
    ;
  }
  else {
    S =     (3*Power(E,2*rij*(xii + xij))*Power(Power(xii,2) - Power(xij,2),7) + 
        Power(E,2*rij*xii)*Power(xii,6)*
         (-21*Power(xii,4)*Power(xij,4)*
            (3 + 6*rij*xij + 10*Power(rij,2)*Power(xij,2) + 
              2*Power(rij,3)*Power(xij,3)) + 
           Power(xii,2)*Power(xij,6)*
            (195 + 390*rij*xij + 78*Power(rij,2)*Power(xij,2) - 
              14*Power(rij,3)*Power(xij,3) - 4*Power(rij,4)*Power(xij,4)) + 
           2*Power(xij,8)*(45 + 90*rij*xij + 48*Power(rij,2)*Power(xij,2) + 
              11*Power(rij,3)*Power(xij,3) + Power(rij,4)*Power(xij,4)) - 
           Power(xii,8)*(3 + 6*rij*xij + 6*Power(rij,2)*Power(xij,2) + 
              4*Power(rij,3)*Power(xij,3) + 2*Power(rij,4)*Power(xij,4)) + 
           Power(xii,6)*Power(xij,2)*
            (21 + 42*rij*xij + 42*Power(rij,2)*Power(xij,2) + 
              38*Power(rij,3)*Power(xij,3) + 4*Power(rij,4)*Power(xij,4))) + 
        Power(E,2*rij*xij)*Power(xij,6)*
         (-22*Power(rij,3)*Power(xii,11) - 2*Power(rij,4)*Power(xii,12) + 
           3*Power(xij,8) + 6*rij*xii*Power(xij,8) + 
           4*Power(rij,2)*Power(xii,10)*(-24 + Power(rij,2)*Power(xij,2)) + 
           3*Power(xii,2)*Power(xij,6)*(-7 + 2*Power(rij,2)*Power(xij,2)) + 
           2*rij*Power(xii,9)*(-90 + 7*Power(rij,2)*Power(xij,2)) + 
           6*rij*Power(xii,7)*Power(xij,2)*(-65 + 7*Power(rij,2)*Power(xij,2)) - 
           6*Power(xii,8)*(15 + 13*Power(rij,2)*Power(xij,2)) + 
           Power(xii,4)*Power(xij,4)*
            (63 - 42*Power(rij,2)*Power(xij,2) + 2*Power(rij,4)*Power(xij,4)) + 
           Power(xii,5)*(126*rij*Power(xij,4) - 38*Power(rij,3)*Power(xij,6)) + 
           Power(xii,6)*(-195*Power(xij,2) + 210*Power(rij,2)*Power(xij,4) - 
              4*Power(rij,4)*Power(xij,6)) + 
           Power(xii,3)*(-42*rij*Power(xij,6) + 4*Power(rij,3)*Power(xij,8))))/
      (3.*Power(E,2*rij*(xii + xij))*Power(rij,2)*
        Power(Power(xii,2) - Power(xij,2),7))
    ;
  }
  return S;
}

static double DSlater_3S_2S(double rij,double xii,double xij)
{
  double S;

  if (xii == xij) {
    S =     (-544320 + 544320*Power(E,2*rij*xii) - 1088640*rij*xii - 
        1088640*Power(rij,2)*Power(xii,2) - 719145*Power(rij,3)*Power(xii,3) - 
        349650*Power(rij,4)*Power(xii,4) - 132111*Power(rij,5)*Power(xii,5) - 
        39942*Power(rij,6)*Power(xii,6) - 9792*Power(rij,7)*Power(xii,7) - 
        1944*Power(rij,8)*Power(xii,8) - 304*Power(rij,9)*Power(xii,9) - 
        32*Power(rij,10)*Power(xii,10))/(544320.*Power(E,2*rij*xii)*Power(rij,2))
    ;
  }
  else {
    S =     (45*Power(E,2*rij*(xii + xij))*Power(Power(xii,2) - Power(xij,2),9) - 
        5*Power(E,2*rij*xii)*Power(xii,8)*
         (-36*Power(xii,6)*Power(xij,4)*
            (-9 - 18*rij*xij - 39*Power(rij,2)*Power(xij,2) + 
              2*Power(rij,3)*Power(xij,3) + Power(rij,4)*Power(xij,4)) + 
           36*Power(xii,4)*Power(xij,6)*
            (-54 - 108*rij*xij + 33*Power(rij,2)*Power(xij,2) + 
              14*Power(rij,3)*Power(xij,3) + Power(rij,4)*Power(xij,4)) + 
           Power(xii,8)*Power(xij,2)*
            (-81 - 162*rij*xij - 162*Power(rij,2)*Power(xij,2) - 
              164*Power(rij,3)*Power(xij,3) + 2*Power(rij,4)*Power(xij,4)) + 
           3*Power(xii,10)*(3 + 6*rij*xij + 6*Power(rij,2)*Power(xij,2) + 
              4*Power(rij,3)*Power(xij,3) + 2*Power(rij,4)*Power(xij,4)) - 
           3*Power(xij,10)*(165 + 330*rij*xij + 
              150*Power(rij,2)*Power(xij,2) + 28*Power(rij,3)*Power(xij,3) + 
              2*Power(rij,4)*Power(xij,4)) - 
           Power(xii,2)*Power(xij,8)*
            (3189 + 6378*rij*xij + 1998*Power(rij,2)*Power(xij,2) + 
              196*Power(rij,3)*Power(xij,3) + 2*Power(rij,4)*Power(xij,4))) + 
        Power(E,2*rij*xij)*Power(xij,6)*
         (-72*Power(rij,5)*Power(xii,17) - 4*Power(rij,6)*Power(xii,18) + 
           45*Power(xij,12) + 90*rij*xii*Power(xij,12) + 
           8*Power(rij,4)*Power(xii,16)*(-75 + 2*Power(rij,2)*Power(xij,2)) + 
           30*rij*Power(xii,3)*Power(xij,10)*
            (-27 + 2*Power(rij,2)*Power(xij,2)) + 
           45*Power(xii,2)*Power(xij,10)*(-9 + 2*Power(rij,2)*Power(xij,2)) + 
           4*Power(rij,3)*Power(xii,15)*(-720 + 47*Power(rij,2)*Power(xij,2)) + 
           12*rij*Power(xii,5)*Power(xij,8)*
            (270 - 45*Power(rij,2)*Power(xij,2) + Power(rij,4)*Power(xij,4)) + 
           30*Power(xii,4)*Power(xij,8)*
            (54 - 27*Power(rij,2)*Power(xij,2) + Power(rij,4)*Power(xij,4)) - 
           20*rij*Power(xii,13)*(594 + 79*Power(rij,2)*Power(xij,2) + 
              Power(rij,4)*Power(xij,4)) - 
           10*Power(rij,2)*Power(xii,14)*
            (810 - 65*Power(rij,2)*Power(xij,2) + 2*Power(rij,4)*Power(xij,4)) - 
           30*rij*Power(xii,11)*Power(xij,2)*
            (1441 - 426*Power(rij,2)*Power(xij,2) + 
              12*Power(rij,4)*Power(xij,4)) - 
           4*rij*Power(xii,7)*Power(xij,6)*
            (1890 - 470*Power(rij,2)*Power(xij,2) + 
              37*Power(rij,4)*Power(xij,4)) + 
           10*rij*Power(xii,9)*Power(xij,4)*
            (639 - 972*Power(rij,2)*Power(xij,2) + 40*Power(rij,4)*Power(xij,4)) \
    + 30*Power(xii,12)*(-198 - 595*Power(rij,2)*Power(xij,2) + 
              61*Power(rij,4)*Power(xij,4)) + 
           Power(xii,8)*Power(xij,4)*
            (3195 - 4860*Power(rij,2)*Power(xij,2) + 
              1780*Power(rij,4)*Power(xij,4) - 16*Power(rij,6)*Power(xij,6)) + 
           2*Power(xii,6)*Power(xij,6)*
            (-1890 + 1620*Power(rij,2)*Power(xij,2) - 
              135*Power(rij,4)*Power(xij,4) + 2*Power(rij,6)*Power(xij,6)) + 
           5*Power(xii,10)*Power(xij,2)*
            (-4323 + 5658*Power(rij,2)*Power(xij,2) - 
              684*Power(rij,4)*Power(xij,4) + 4*Power(rij,6)*Power(xij,6))))/
      (45.*Power(E,2*rij*(xii + xij))*Power(rij,2)*
        Power(Power(xii,2) - Power(xij,2),9))
    ;
  }
  return S;
}

static double DSlater_4S_2S(double rij,double xii,double xij)
{
  double S;

  if (xii == xij) {
    S =     (-159667200 + 159667200*Power(E,2*rij*xii) - 319334400*rij*xii - 
        319334400*Power(rij,2)*Power(xii,2) - 
        212109975*Power(rij,3)*Power(xii,3) - 
        104885550*Power(rij,4)*Power(xii,4) - 
        40997880*Power(rij,5)*Power(xii,5) - 
        13111560*Power(rij,6)*Power(xii,6) - 3496680*Power(rij,7)*Power(xii,7) - 
        784080*Power(rij,8)*Power(xii,8) - 147840*Power(rij,9)*Power(xii,9) - 
        23232*Power(rij,10)*Power(xii,10) - 2944*Power(rij,11)*Power(xii,11) - 
        256*Power(rij,12)*Power(xii,12))/
      (1.596672e8*Power(E,2*rij*xii)*Power(rij,2))
    ;
  }
  else {
    S =     (315*Power(E,2*rij*(xii + xij))*Power(Power(xii,2) - Power(xij,2),11) - 
        105*Power(E,2*rij*xii)*Power(xii,10)*
         (198*Power(xii,6)*Power(xij,6)*
            (-9 - 18*rij*xij + 14*Power(rij,2)*Power(xij,2) + 
              2*Power(rij,3)*Power(xij,3)) + 
           Power(xii,12)*(3 + 6*rij*xij + 6*Power(rij,2)*Power(xij,2) + 
              4*Power(rij,3)*Power(xij,3) + 2*Power(rij,4)*Power(xij,4)) + 
           11*Power(xii,4)*Power(xij,8)*
            (-573 - 1146*rij*xij - 126*Power(rij,2)*Power(xij,2) + 
              16*Power(rij,3)*Power(xij,3) + 2*Power(rij,4)*Power(xij,4)) - 
           11*Power(xii,8)*Power(xij,4)*
            (-15 - 30*rij*xij - 84*Power(rij,2)*Power(xij,2) + 
              22*Power(rij,3)*Power(xij,3) + 2*Power(rij,4)*Power(xij,4)) - 
           Power(xij,12)*(273 + 546*rij*xij + 216*Power(rij,2)*Power(xij,2) + 
              34*Power(rij,3)*Power(xij,3) + 2*Power(rij,4)*Power(xij,4)) + 
           Power(xii,10)*Power(xij,2)*
            (-33 - 66*rij*xij - 66*Power(rij,2)*Power(xij,2) - 
              74*Power(rij,3)*Power(xij,3) + 8*Power(rij,4)*Power(xij,4)) - 
           Power(xii,2)*Power(xij,10)*
            (3297 + 6594*rij*xij + 2034*Power(rij,2)*Power(xij,2) + 
              226*Power(rij,3)*Power(xij,3) + 8*Power(rij,4)*Power(xij,4))) + 
        Power(E,2*rij*xij)*Power(xij,6)*
         (-50*Power(rij,7)*Power(xii,23) - 2*Power(rij,8)*Power(xii,24) + 
           315*Power(xij,16) + 630*rij*xii*Power(xij,16) + 
           210*rij*Power(xii,3)*Power(xij,14)*
            (-33 + 2*Power(rij,2)*Power(xij,2)) + 
           315*Power(xii,2)*Power(xij,14)*(-11 + 2*Power(rij,2)*Power(xij,2)) + 
           4*Power(rij,6)*Power(xii,22)*(-154 + 3*Power(rij,2)*Power(xij,2)) + 
           Power(xii,21)*(-4788*Power(rij,5) + 230*Power(rij,7)*Power(xij,2)) - 
           28*Power(rij,4)*Power(xii,20)*
            (900 - 67*Power(rij,2)*Power(xij,2) + Power(rij,4)*Power(xij,4)) + 
           42*rij*Power(xii,5)*Power(xij,12)*
            (825 - 110*Power(rij,2)*Power(xij,2) + 2*Power(rij,4)*Power(xij,4)) \
    + 105*Power(xii,4)*Power(xij,12)*
            (165 - 66*Power(rij,2)*Power(xij,2) + 2*Power(rij,4)*Power(xij,4)) - 
           14*Power(rij,3)*Power(xii,19)*
            (6435 - 462*Power(rij,2)*Power(xij,2) + 
              23*Power(rij,4)*Power(xij,4)) + 
           2*rij*Power(xii,9)*Power(xij,8)*
            (103950 - 36225*Power(rij,2)*Power(xij,2) + 
              2016*Power(rij,4)*Power(xij,4) - 59*Power(rij,6)*Power(xij,6)) + 
           14*Power(rij,2)*Power(xii,18)*
            (-14850 - 825*Power(rij,2)*Power(xij,2) + 
              2*Power(rij,4)*Power(xij,4) + 2*Power(rij,6)*Power(xij,6)) + 
           2*rij*Power(xii,7)*Power(xij,10)*
            (-51975 + 11550*Power(rij,2)*Power(xij,2) - 
              462*Power(rij,4)*Power(xij,4) + 4*Power(rij,6)*Power(xij,6)) + 
           7*Power(xii,6)*Power(xij,10)*
            (-7425 + 4950*Power(rij,2)*Power(xij,2) - 
              330*Power(rij,4)*Power(xij,4) + 4*Power(rij,6)*Power(xij,6)) - 
           14*rij*Power(xii,17)*(19305 + 15015*Power(rij,2)*Power(xij,2) - 
              1656*Power(rij,4)*Power(xij,4) + 7*Power(rij,6)*Power(xij,6)) - 
           70*rij*Power(xii,13)*Power(xij,4)*
            (5733 + 4488*Power(rij,2)*Power(xij,2) - 
              792*Power(rij,4)*Power(xij,4) + 13*Power(rij,6)*Power(xij,6)) + 
           14*rij*Power(xii,11)*Power(xij,6)*
            (-24885 + 2475*Power(rij,2)*Power(xij,2) - 
              1518*Power(rij,4)*Power(xij,4) + 35*Power(rij,6)*Power(xij,6)) + 
           14*rij*Power(xii,15)*Power(xij,2)*
            (-109395 + 45240*Power(rij,2)*Power(xij,2) - 
              4446*Power(rij,4)*Power(xij,4) + 55*Power(rij,6)*Power(xij,6)) - 
           7*Power(xii,16)*(19305 + 123750*Power(rij,2)*Power(xij,2) - 
              26430*Power(rij,4)*Power(xij,4) + 964*Power(rij,6)*Power(xij,6)) + 
           Power(xii,10)*Power(xij,6)*
            (-174195 + 242550*Power(rij,2)*Power(xij,2) - 
              22050*Power(rij,4)*Power(xij,4) + 
              2324*Power(rij,6)*Power(xij,6) - 12*Power(rij,8)*Power(xij,8)) + 
           2*Power(xii,8)*Power(xij,8)*
            (51975 - 51975*Power(rij,2)*Power(xij,2) + 
              5775*Power(rij,4)*Power(xij,4) - 154*Power(rij,6)*Power(xij,6) + 
              Power(rij,8)*Power(xij,8)) - 
           7*Power(xii,14)*Power(xij,2)*
            (109395 - 110970*Power(rij,2)*Power(xij,2) + 
              34230*Power(rij,4)*Power(xij,4) - 
              1540*Power(rij,6)*Power(xij,6) + 4*Power(rij,8)*Power(xij,8)) + 
           7*Power(xii,12)*Power(xij,4)*
            (-28665 + 18630*Power(rij,2)*Power(xij,2) + 
              14850*Power(rij,4)*Power(xij,4) - 1052*Power(rij,6)*Power(xij,6) + 
              4*Power(rij,8)*Power(xij,8))))/
      (315.*Power(E,2*rij*(xii + xij))*Power(rij,2)*
        Power(Power(xii,2) - Power(xij,2),11))
    ;
  }
  return S;
}

double DSlater_1S_2S(double rij,double xii,double xij)
{
  return DSlater_2S_1S(rij,xij,xii);
}

static double DSlater_3S_3S(double rij,double xii,double xij)
{
  double S;

  if (xii == xij) {
    S =     (-359251200 + 359251200*Power(E,2*rij*xii) - 718502400*rij*xii - 
        718502400*Power(rij,2)*Power(xii,2) - 
        475727175*Power(rij,3)*Power(xii,3) - 
        232951950*Power(rij,4)*Power(xii,4) - 
        89397000*Power(rij,5)*Power(xii,5) - 
        27858600*Power(rij,6)*Power(xii,6) - 7223040*Power(rij,7)*Power(xii,7) - 
        1584000*Power(rij,8)*Power(xii,8) - 295680*Power(rij,9)*Power(xii,9) - 
        46464*Power(rij,10)*Power(xii,10) - 5888*Power(rij,11)*Power(xii,11) - 
        512*Power(rij,12)*Power(xii,12))/
      (3.592512e8*Power(E,2*rij*xii)*Power(rij,2))
    ;
  }
  else {
    S =     (135*Power(E,2*rij*(xii + xij))*Power(Power(xii,2) - Power(xij,2),11) + 
        Power(E,2*rij*xij)*Power(xij,8)*
         (-276*Power(rij,5)*Power(xii,19) - 12*Power(rij,6)*Power(xii,20) + 
           135*Power(xij,14) + 270*rij*xii*Power(xij,14) - 
           100*Power(rij,3)*Power(xii,17)*(165 + Power(rij,2)*Power(xij,2)) + 
           10*Power(rij,4)*Power(xii,18)*(-285 + 2*Power(rij,2)*Power(xij,2)) + 
           90*rij*Power(xii,3)*Power(xij,12)*
            (-33 + 2*Power(rij,2)*Power(xij,2)) + 
           135*Power(xii,2)*Power(xij,12)*(-11 + 2*Power(rij,2)*Power(xij,2)) + 
           18*rij*Power(xii,5)*Power(xij,10)*
            (825 - 110*Power(rij,2)*Power(xij,2) + 2*Power(rij,4)*Power(xij,4)) \
    + 45*Power(xii,4)*Power(xij,10)*
            (165 - 66*Power(rij,2)*Power(xij,2) + 2*Power(rij,4)*Power(xij,4)) \
    - 10*rij*Power(xii,7)*Power(xij,8)*
            (4455 - 738*Power(rij,2)*Power(xij,2) + 
              62*Power(rij,4)*Power(xij,4)) + 
           10*rij*Power(xii,11)*Power(xij,4)*
            (-96831 + 6534*Power(rij,2)*Power(xij,2) + 
              154*Power(rij,4)*Power(xij,4)) - 
           10*rij*Power(xii,13)*Power(xij,2)*
            (84357 - 12318*Power(rij,2)*Power(xij,2) + 
              418*Power(rij,4)*Power(xij,4)) + 
           2*rij*Power(xii,9)*Power(xij,6)*
            (-495 - 48510*Power(rij,2)*Power(xij,2) + 
              458*Power(rij,4)*Power(xij,4)) + 
           Power(xii,15)*(-90090*rij - 80580*Power(rij,3)*Power(xij,2) + 
              2684*Power(rij,5)*Power(xij,4)) + 
           Power(xii,16)*(-54450*Power(rij,2) - 
              7290*Power(rij,4)*Power(xij,2) + 68*Power(rij,6)*Power(xij,4)) - 
           5*Power(xii,8)*Power(xij,6)*
            (99 + 1782*Power(rij,2)*Power(xij,2) - 
              2250*Power(rij,4)*Power(xij,4) + 4*Power(rij,6)*Power(xij,6)) + 
           3*Power(xii,6)*Power(xij,8)*
            (-7425 + 4950*Power(rij,2)*Power(xij,2) - 
              330*Power(rij,4)*Power(xij,4) + 4*Power(rij,6)*Power(xij,6)) - 
           5*Power(xii,14)*(9009 + 78954*Power(rij,2)*Power(xij,2) - 
              6030*Power(rij,4)*Power(xij,4) + 44*Power(rij,6)*Power(xij,6)) + 
           5*Power(xii,12)*Power(xij,2)*
            (-84357 - 366*Power(rij,2)*Power(xij,2) - 
              3498*Power(rij,4)*Power(xij,4) + 44*Power(rij,6)*Power(xij,6)) - 
           Power(xii,10)*Power(xij,4)*
            (484155 - 447810*Power(rij,2)*Power(xij,2) + 
              12870*Power(rij,4)*Power(xij,4) + 68*Power(rij,6)*Power(xij,6))) + 
        Power(E,2*rij*xii)*Power(xii,8)*
         (Power(xii,4)*Power(xij,10)*
            (484155 + 968310*rij*xij + 1830*Power(rij,2)*Power(xij,2) - 
              123180*Power(rij,3)*Power(xij,3) - 
              30150*Power(rij,4)*Power(xij,4) - 
              2684*Power(rij,5)*Power(xij,5) - 68*Power(rij,6)*Power(xij,6)) + 
           5*Power(xii,2)*Power(xij,12)*
            (84357 + 168714*rij*xij + 78954*Power(rij,2)*Power(xij,2) + 
              16116*Power(rij,3)*Power(xij,3) + 
              1458*Power(rij,4)*Power(xij,4) + 20*Power(rij,5)*Power(xij,5) - 
              4*Power(rij,6)*Power(xij,6)) - 
           3*Power(xii,14)*(45 + 90*rij*xij + 90*Power(rij,2)*Power(xij,2) + 
              60*Power(rij,3)*Power(xij,3) + 30*Power(rij,4)*Power(xij,4) + 
              12*Power(rij,5)*Power(xij,5) + 4*Power(rij,6)*Power(xij,6)) - 
           55*Power(xii,8)*Power(xij,6)*
            (-405 - 810*rij*xij - 162*Power(rij,2)*Power(xij,2) - 
              1764*Power(rij,3)*Power(xij,3) - 234*Power(rij,4)*Power(xij,4) + 
              28*Power(rij,5)*Power(xij,5) + 4*Power(rij,6)*Power(xij,6)) + 
           55*Power(xii,6)*Power(xij,8)*
            (9 + 18*rij*xij - 8142*Power(rij,2)*Power(xij,2) - 
              1188*Power(rij,3)*Power(xij,3) + 318*Power(rij,4)*Power(xij,4) + 
              76*Power(rij,5)*Power(xij,5) + 4*Power(rij,6)*Power(xij,6)) + 
           3*Power(xij,14)*(15015 + 30030*rij*xij + 
              18150*Power(rij,2)*Power(xij,2) + 
              5500*Power(rij,3)*Power(xij,3) + 950*Power(rij,4)*Power(xij,4) + 
              92*Power(rij,5)*Power(xij,5) + 4*Power(rij,6)*Power(xij,6)) + 
           5*Power(xii,12)*Power(xij,2)*
            (297 + 594*rij*xij + 594*Power(rij,2)*Power(xij,2) + 
              396*Power(rij,3)*Power(xij,3) + 198*Power(rij,4)*Power(xij,4) + 
              124*Power(rij,5)*Power(xij,5) + 4*Power(rij,6)*Power(xij,6)) + 
           Power(xii,10)*Power(xij,4)*
            (-7425 - 14850*rij*xij - 14850*Power(rij,2)*Power(xij,2) - 
              7380*Power(rij,3)*Power(xij,3) - 11250*Power(rij,4)*Power(xij,4) - 
              916*Power(rij,5)*Power(xij,5) + 68*Power(rij,6)*Power(xij,6))))/
      (135.*Power(E,2*rij*(xii + xij))*Power(rij,2)*
        Power(Power(xii,2) - Power(xij,2),11))
    ;
  }
  return S;
}

static double DSlater_4S_3S(double rij,double xii,double xij)
{
  double S;

  if (xii == xij) {
    S =     (-37362124800.0 + 37362124800.0*Power(E,2*rij*xii) - 74724249600.0*rij*xii - 
        74724249600.0*Power(rij,2)*Power(xii,2) - 
        49615490925.0*Power(rij,3)*Power(xii,3) - 
        24506732250.0*Power(rij,4)*Power(xii,4) - 
        9566747190.0*Power(rij,5)*Power(xii,5) - 
        3063240180.0*Power(rij,6)*Power(xii,6) - 
        824709600*Power(rij,7)*Power(xii,7) - 
        189961200*Power(rij,8)*Power(xii,8) - 
        37889280*Power(rij,9)*Power(xii,9) - 
        6589440*Power(rij,10)*Power(xii,10) - 
        998400*Power(rij,11)*Power(xii,11) - 
        129792*Power(rij,12)*Power(xii,12) - 13824*Power(rij,13)*Power(xii,13) - 
        1024*Power(rij,14)*Power(xii,14))/
      (3.73621248e10*Power(E,2*rij*xii)*Power(rij,2))
    ;
  }
  else {
    S =     (945*Power(E,2*rij*(xii + xij))*Power(Power(xii,2) - Power(xij,2),13) - 
        21*Power(E,2*rij*xii)*Power(xii,10)*
         (-715*Power(xii,8)*Power(xij,8)*
            (45 + 90*rij*xij - 888*Power(rij,2)*Power(xij,2) + 
              84*Power(rij,3)*Power(xij,3) + 54*Power(rij,4)*Power(xij,4) + 
              4*Power(rij,5)*Power(xij,5)) + 
           6*Power(xii,12)*Power(xij,4)*
            (585 + 1170*rij*xij + 1170*Power(rij,2)*Power(xij,2) + 
              450*Power(rij,3)*Power(xij,3) + 1215*Power(rij,4)*Power(xij,4) - 
              46*Power(rij,5)*Power(xij,5) - 12*Power(rij,6)*Power(xij,6)) - 
           2*Power(xij,16)*(16380 + 32760*rij*xij + 
              17745*Power(rij,2)*Power(xij,2) + 
              4680*Power(rij,3)*Power(xij,3) + 690*Power(rij,4)*Power(xij,4) + 
              56*Power(rij,5)*Power(xij,5) + 2*Power(rij,6)*Power(xij,6)) + 
           Power(xii,16)*(45 + 90*rij*xij + 90*Power(rij,2)*Power(xij,2) + 
              60*Power(rij,3)*Power(xij,3) + 30*Power(rij,4)*Power(xij,4) + 
              12*Power(rij,5)*Power(xij,5) + 4*Power(rij,6)*Power(xij,6)) + 
           26*Power(xii,10)*Power(xij,6)*
            (-495 - 990*rij*xij + 495*Power(rij,2)*Power(xij,2) - 
              3630*Power(rij,3)*Power(xij,3) + 165*Power(rij,4)*Power(xij,4) + 
              102*Power(rij,5)*Power(xij,5) + 4*Power(rij,6)*Power(xij,6)) - 
           2*Power(xii,2)*Power(xij,14)*
            (269010 + 538020*rij*xij + 241815*Power(rij,2)*Power(xij,2) + 
              49410*Power(rij,3)*Power(xij,3) + 
              5205*Power(rij,4)*Power(xij,4) + 262*Power(rij,5)*Power(xij,5) + 
              4*Power(rij,6)*Power(xij,6)) + 
           Power(xii,14)*Power(xij,2)*
            (-585 - 1170*rij*xij - 1170*Power(rij,2)*Power(xij,2) - 
              780*Power(rij,3)*Power(xij,3) - 390*Power(rij,4)*Power(xij,4) - 
              276*Power(rij,5)*Power(xij,5) + 8*Power(rij,6)*Power(xij,6)) - 
           13*Power(xii,6)*Power(xij,10)*
            (68625 + 137250*rij*xij - 49890*Power(rij,2)*Power(xij,2) - 
              21420*Power(rij,3)*Power(xij,3) - 
              1830*Power(rij,4)*Power(xij,4) + 44*Power(rij,5)*Power(xij,5) + 
              8*Power(rij,6)*Power(xij,6)) + 
           6*Power(xii,4)*Power(xij,12)*
            (-256050 - 512100*rij*xij - 130530*Power(rij,2)*Power(xij,2) - 
              2970*Power(rij,3)*Power(xij,3) + 2565*Power(rij,4)*Power(xij,4) + 
              326*Power(rij,5)*Power(xij,5) + 12*Power(rij,6)*Power(xij,6))) + 
        Power(E,2*rij*xij)*Power(xij,8)*
         (-192*Power(rij,7)*Power(xii,25) - 6*Power(rij,8)*Power(xii,26) + 
           945*Power(xij,18) + 1890*rij*xii*Power(xij,18) + 
           630*rij*Power(xii,3)*Power(xij,16)*
            (-39 + 2*Power(rij,2)*Power(xij,2)) + 
           945*Power(xii,2)*Power(xij,16)*(-13 + 2*Power(rij,2)*Power(xij,2)) + 
           24*Power(rij,5)*Power(xii,23)*
            (-1155 + 13*Power(rij,2)*Power(xij,2)) + 
           Power(xii,24)*(-2940*Power(rij,6) + 22*Power(rij,8)*Power(xij,2)) + 
           252*rij*Power(xii,5)*Power(xij,14)*
            (585 - 65*Power(rij,2)*Power(xij,2) + Power(rij,4)*Power(xij,4)) + 
           630*Power(xii,4)*Power(xij,14)*
            (117 - 39*Power(rij,2)*Power(xij,2) + Power(rij,4)*Power(xij,4)) + 
           2*Power(rij,4)*Power(xii,22)*
            (-86625 - 882*Power(rij,2)*Power(xij,2) + 
              4*Power(rij,4)*Power(xij,4)) + 
           48*Power(rij,3)*Power(xii,21)*
            (-15015 - 1729*Power(rij,2)*Power(xij,2) + 
              38*Power(rij,4)*Power(xij,4)) + 
           12*rij*Power(xii,7)*Power(xij,12)*
            (-45045 + 8190*Power(rij,2)*Power(xij,2) - 
              273*Power(rij,4)*Power(xij,4) + 2*Power(rij,6)*Power(xij,6)) + 
           42*Power(xii,6)*Power(xij,12)*
            (-6435 + 3510*Power(rij,2)*Power(xij,2) - 
              195*Power(rij,4)*Power(xij,4) + 2*Power(rij,6)*Power(xij,6)) - 
           42*Power(rij,2)*Power(xii,20)*
            (45045 + 22815*Power(rij,2)*Power(xij,2) - 
              1036*Power(rij,4)*Power(xij,4) + 4*Power(rij,6)*Power(xij,6)) - 
           84*rij*Power(xii,13)*Power(xij,6)*
            (230400 + 87945*Power(rij,2)*Power(xij,2) - 
              5525*Power(rij,4)*Power(xij,4) + 8*Power(rij,6)*Power(xij,6)) - 
           84*rij*Power(xii,15)*Power(xij,4)*
            (806760 - 141690*Power(rij,2)*Power(xij,2) + 
              2483*Power(rij,4)*Power(xij,4) + 52*Power(rij,6)*Power(xij,6)) - 
           84*rij*Power(xii,19)*(32175 + 70005*Power(rij,2)*Power(xij,2) - 
              5389*Power(rij,4)*Power(xij,4) + 80*Power(rij,6)*Power(xij,6)) - 
           6*rij*Power(xii,9)*Power(xij,10)*
            (-225225 + 66990*Power(rij,2)*Power(xij,2) - 
              2394*Power(rij,4)*Power(xij,4) + 80*Power(rij,6)*Power(xij,6)) + 
           84*rij*Power(xii,17)*Power(xij,2)*
            (-417105 + 41505*Power(rij,2)*Power(xij,2) - 
              5031*Power(rij,4)*Power(xij,4) + 104*Power(rij,6)*Power(xij,6)) + 
           6*rij*Power(xii,11)*Power(xij,8)*
            (-634725 - 180180*Power(rij,2)*Power(xij,2) - 
              31038*Power(rij,4)*Power(xij,4) + 256*Power(rij,6)*Power(xij,6)) + 
           Power(xii,10)*Power(xij,8)*
            (-1904175 + 1981980*Power(rij,2)*Power(xij,2) - 
              13860*Power(rij,4)*Power(xij,4) + 
              13608*Power(rij,6)*Power(xij,6) - 22*Power(rij,8)*Power(xij,8)) + 
           3*Power(xii,8)*Power(xij,10)*
            (225225 - 180180*Power(rij,2)*Power(xij,2) + 
              16380*Power(rij,4)*Power(xij,4) - 364*Power(rij,6)*Power(xij,6) + 
              2*Power(rij,8)*Power(xij,8)) + 
           84*Power(xii,14)*Power(xij,4)*
            (-403380 + 341955*Power(rij,2)*Power(xij,2) - 
              37440*Power(rij,4)*Power(xij,4) + 260*Power(rij,6)*Power(xij,6) + 
              2*Power(rij,8)*Power(xij,8)) - 
           4*Power(xii,12)*Power(xij,6)*
            (2419200 - 2158065*Power(rij,2)*Power(xij,2) - 
              360360*Power(rij,4)*Power(xij,4) + 
              9534*Power(rij,6)*Power(xij,6) + 2*Power(rij,8)*Power(xij,8)) + 
           14*Power(xii,18)*(-96525 - 1453725*Power(rij,2)*Power(xij,2) + 
              163710*Power(rij,4)*Power(xij,4) - 
              6252*Power(rij,6)*Power(xij,6) + 26*Power(rij,8)*Power(xij,8)) - 
           14*Power(xii,16)*Power(xij,2)*
            (1251315 + 1191330*Power(rij,2)*Power(xij,2) - 
              36810*Power(rij,4)*Power(xij,4) - 3744*Power(rij,6)*Power(xij,6) + 
              26*Power(rij,8)*Power(xij,8))))/
      (945.*Power(E,2*rij*(xii + xij))*Power(rij,2)*
        Power(Power(xii,2) - Power(xij,2),13))
    ;
  }
  return S;
}

double DSlater_1S_3S(double rij,double xii,double xij)
{
  return DSlater_3S_1S(rij,xij,xii);
}

double DSlater_2S_3S(double rij,double xii,double xij)
{
  return DSlater_3S_2S(rij,xij,xii);
}

static double DSlater_4S_4S(double rij,double xii,double xij)
{
  double S;

  if (xii == xij) {
    S =     (-20922789888000.0 + 20922789888000.0*Power(E,2*rij*xii) - 
        41845579776000.0*rij*xii - 41845579776000.0*Power(rij,2)*Power(xii,2) - 
        27805745842875.0*Power(rij,3)*Power(xii,3) - 
        13765911909750.0*Power(rij,4)*Power(xii,4) - 
        5399605411200.0*Power(rij,5)*Power(xii,5) - 
        1743679337400.0*Power(rij,6)*Power(xii,6) - 
        476010334800.0*Power(rij,7)*Power(xii,7) - 
        112021509600.0*Power(rij,8)*Power(xii,8) - 
        23063040000.0*Power(rij,9)*Power(xii,9) - 
        4197473280.0*Power(rij,10)*Power(xii,10) - 
        679311360*Power(rij,11)*Power(xii,11) - 
        97843200*Power(rij,12)*Power(xii,12) - 
        12472320*Power(rij,13)*Power(xii,13) - 
        1382400*Power(rij,14)*Power(xii,14) - 
        126976*Power(rij,15)*Power(xii,15) - 8192*Power(rij,16)*Power(xii,16))/
      (2.0922789888e13*Power(E,2*rij*xii)*Power(rij,2))
    ;
  }
  else {
    S =     (315*Power(E,2*rij*(xii + xij))*Power(Power(xii,2) - Power(xij,2),15) + 
        Power(E,2*rij*xij)*Power(xij,10)*
         (-1428*Power(rij,6)*Power(xii,26) - 78*Power(rij,7)*Power(xii,27) - 
           2*Power(rij,8)*Power(xii,28) + 315*Power(xij,20) + 
           630*rij*xii*Power(xij,20) + 
           210*rij*Power(xii,3)*Power(xij,18)*
            (-45 + 2*Power(rij,2)*Power(xij,2)) + 
           315*Power(xii,2)*Power(xij,18)*(-15 + 2*Power(rij,2)*Power(xij,2)) - 
           42*Power(rij,5)*Power(xii,25)*(377 + 5*Power(rij,2)*Power(xij,2)) + 
           42*Power(rij,4)*Power(xii,24)*
            (-2730 - 190*Power(rij,2)*Power(xij,2) + Power(rij,4)*Power(xij,4)) \
    + 42*rij*Power(xii,5)*Power(xij,16)*
            (1575 - 150*Power(rij,2)*Power(xij,2) + 
              2*Power(rij,4)*Power(xij,4)) + 
           105*Power(xii,4)*Power(xij,16)*
            (315 - 90*Power(rij,2)*Power(xij,2) + 2*Power(rij,4)*Power(xij,4)) \
    + 63*Power(rij,3)*Power(xii,23)*
            (-8645 - 2180*Power(rij,2)*Power(xij,2) + 
              32*Power(rij,4)*Power(xij,4)) + 
           2*rij*Power(xii,7)*Power(xij,14)*
            (-143325 + 22050*Power(rij,2)*Power(xij,2) - 
              630*Power(rij,4)*Power(xij,4) + 4*Power(rij,6)*Power(xij,6)) + 
           7*Power(xii,6)*Power(xij,14)*
            (-20475 + 9450*Power(rij,2)*Power(xij,2) - 
              450*Power(rij,4)*Power(xij,4) + 4*Power(rij,6)*Power(xij,6)) + 
           21*rij*Power(xii,11)*Power(xij,10)*
            (-209430 - 112125*Power(rij,2)*Power(xij,2) - 
              8288*Power(rij,4)*Power(xij,4) + 10*Power(rij,6)*Power(xij,6)) - 
           21*rij*Power(xii,9)*Power(xij,12)*
            (-40950 + 11245*Power(rij,2)*Power(xij,2) - 
              222*Power(rij,4)*Power(xij,4) + 10*Power(rij,6)*Power(xij,6)) + 
           7*rij*Power(xii,19)*Power(xij,2)*
            (-7711200 - 1605825*Power(rij,2)*Power(xij,2) + 
              55104*Power(rij,4)*Power(xij,4) + 20*Power(rij,6)*Power(xij,6)) - 
           4*Power(rij,2)*Power(xii,22)*
            (401310 + 341775*Power(rij,2)*Power(xij,2) - 
              9009*Power(rij,4)*Power(xij,4) + 32*Power(rij,6)*Power(xij,6)) + 
           105*rij*Power(xii,17)*Power(xij,4)*
            (-2087532 + 267621*Power(rij,2)*Power(xij,2) - 
              10348*Power(rij,4)*Power(xij,4) + 52*Power(rij,6)*Power(xij,6)) - 
           105*rij*Power(xii,15)*Power(xij,6)*
            (2126142 - 103075*Power(rij,2)*Power(xij,2) - 
              4680*Power(rij,4)*Power(xij,4) + 56*Power(rij,6)*Power(xij,6)) + 
           21*Power(xii,10)*Power(xij,10)*
            (-104715 + 83850*Power(rij,2)*Power(xij,2) + 
              4030*Power(rij,4)*Power(xij,4) + 404*Power(rij,6)*Power(xij,6)) - 
           70*Power(xii,18)*Power(xij,2)*
            (385560 + 1201608*Power(rij,2)*Power(xij,2) - 
              84195*Power(rij,4)*Power(xij,4) + 1064*Power(rij,6)*Power(xij,6)) \
    - 3*rij*Power(xii,21)*(835380 + 2774625*Power(rij,2)*Power(xij,2) - 
              94836*Power(rij,4)*Power(xij,4) + 1160*Power(rij,6)*Power(xij,6)) \
    + rij*Power(xii,13)*Power(xij,8)*
            (-50825250 - 16261245*Power(rij,2)*Power(xij,2) + 
              248640*Power(rij,4)*Power(xij,4) + 2024*Power(rij,6)*Power(xij,6)\
    ) - 70*Power(xii,16)*Power(xij,4)*
            (1565649 - 145035*Power(rij,2)*Power(xij,2) + 
              63465*Power(rij,4)*Power(xij,4) - 
              1560*Power(rij,6)*Power(xij,6) + 2*Power(rij,8)*Power(xij,8)) + 
           Power(xii,8)*Power(xij,12)*
            (429975 - 286650*Power(rij,2)*Power(xij,2) + 
              22050*Power(rij,4)*Power(xij,4) - 
              420*Power(rij,6)*Power(xij,6) + 2*Power(rij,8)*Power(xij,8)) - 
           7*Power(xii,12)*Power(xij,8)*
            (3630375 - 2811510*Power(rij,2)*Power(xij,2) - 
              298350*Power(rij,4)*Power(xij,4) + 
              1688*Power(rij,6)*Power(xij,6) + 6*Power(rij,8)*Power(xij,8)) + 
           14*Power(xii,20)*(-89505 - 2135700*Power(rij,2)*Power(xij,2) + 
              24030*Power(rij,4)*Power(xij,4) - 
              1236*Power(rij,6)*Power(xij,6) + 10*Power(rij,8)*Power(xij,8)) + 
           Power(xii,14)*Power(xij,6)*
            (-111622455 + 84253050*Power(rij,2)*Power(xij,2) - 
              2497950*Power(rij,4)*Power(xij,4) - 
              40320*Power(rij,6)*Power(xij,6) + 128*Power(rij,8)*Power(xij,8))) \
    + Power(E,2*rij*xii)*Power(xii,10)*
         (105*Power(xii,18)*Power(xij,2)*
            (45 + 90*rij*xij + 90*Power(rij,2)*Power(xij,2) + 
              60*Power(rij,3)*Power(xij,3) + 30*Power(rij,4)*Power(xij,4) + 
              12*Power(rij,5)*Power(xij,5) + 4*Power(rij,6)*Power(xij,6) + 
              2*Power(rij,7)*Power(xij,7)) + 
           105*Power(xii,2)*Power(xij,18)*
            (257040 + 514080*rij*xij + 284760*Power(rij,2)*Power(xij,2) + 
              79275*Power(rij,3)*Power(xij,3) + 
              13020*Power(rij,4)*Power(xij,4) + 
              1308*Power(rij,5)*Power(xij,5) + 76*Power(rij,6)*Power(xij,6) + 
              2*Power(rij,7)*Power(xij,7)) - 
           1365*Power(xii,10)*Power(xij,10)*
            (-1611 - 3222*rij*xij + 14418*Power(rij,2)*Power(xij,2) - 
              11913*Power(rij,3)*Power(xij,3) - 
              1830*Power(rij,4)*Power(xij,4) + 360*Power(rij,5)*Power(xij,5) + 
              80*Power(rij,6)*Power(xij,6) + 4*Power(rij,7)*Power(xij,7)) + 
           Power(xii,14)*Power(xij,6)*
            (143325 + 286650*rij*xij + 286650*Power(rij,2)*Power(xij,2) + 
              236145*Power(rij,3)*Power(xij,3) - 
              84630*Power(rij,4)*Power(xij,4) + 
              174048*Power(rij,5)*Power(xij,5) + 
              11816*Power(rij,6)*Power(xij,6) - 
              2024*Power(rij,7)*Power(xij,7) - 128*Power(rij,8)*Power(xij,8)) + 
           21*Power(xii,16)*Power(xij,4)*
            (-1575 - 3150*rij*xij - 3150*Power(rij,2)*Power(xij,2) - 
              2100*Power(rij,3)*Power(xij,3) - 1050*Power(rij,4)*Power(xij,4) - 
              222*Power(rij,5)*Power(xij,5) - 404*Power(rij,6)*Power(xij,6) - 
              10*Power(rij,7)*Power(xij,7) + 2*Power(rij,8)*Power(xij,8)) - 
           Power(xii,20)*(315 + 630*rij*xij + 630*Power(rij,2)*Power(xij,2) + 
              420*Power(rij,3)*Power(xij,3) + 210*Power(rij,4)*Power(xij,4) + 
              84*Power(rij,5)*Power(xij,5) + 28*Power(rij,6)*Power(xij,6) + 
              8*Power(rij,7)*Power(xij,7) + 2*Power(rij,8)*Power(xij,8)) + 
           Power(xij,20)*(1253070 + 2506140*rij*xij + 
              1605240*Power(rij,2)*Power(xij,2) + 
              544635*Power(rij,3)*Power(xij,3) + 
              114660*Power(rij,4)*Power(xij,4) + 
              15834*Power(rij,5)*Power(xij,5) + 
              1428*Power(rij,6)*Power(xij,6) + 78*Power(rij,7)*Power(xij,7) + 
              2*Power(rij,8)*Power(xij,8)) - 
           21*Power(xii,4)*Power(xij,16)*
            (-5218830 - 10437660*rij*xij - 4005360*Power(rij,2)*Power(xij,2) - 
              535275*Power(rij,3)*Power(xij,3) + 
              16020*Power(rij,4)*Power(xij,4) + 
              13548*Power(rij,5)*Power(xij,5) + 
              1716*Power(rij,6)*Power(xij,6) + 96*Power(rij,7)*Power(xij,7) + 
              2*Power(rij,8)*Power(xij,8)) - 
           35*Power(xii,8)*Power(xij,12)*
            (-726075 - 1452150*rij*xij + 2407230*Power(rij,2)*Power(xij,2) + 
              309225*Power(rij,3)*Power(xij,3) - 
              126930*Power(rij,4)*Power(xij,4) - 
              31044*Power(rij,5)*Power(xij,5) - 
              2128*Power(rij,6)*Power(xij,6) + 4*Power(rij,7)*Power(xij,7) + 
              4*Power(rij,8)*Power(xij,8)) + 
           35*Power(xii,12)*Power(xij,8)*
            (-12285 - 24570*rij*xij - 50310*Power(rij,2)*Power(xij,2) + 
              67275*Power(rij,3)*Power(xij,3) - 
              59670*Power(rij,4)*Power(xij,4) - 
              7104*Power(rij,5)*Power(xij,5) + 1152*Power(rij,6)*Power(xij,6) + 
              168*Power(rij,7)*Power(xij,7) + 4*Power(rij,8)*Power(xij,8)) + 
           Power(xii,6)*Power(xij,14)*
            (111622455 + 223244910*rij*xij - 
              10152450*Power(rij,2)*Power(xij,2) - 
              28100205*Power(rij,3)*Power(xij,3) - 
              5893650*Power(rij,4)*Power(xij,4) - 
              385728*Power(rij,5)*Power(xij,5) + 
              17304*Power(rij,6)*Power(xij,6) + 3480*Power(rij,7)*Power(xij,7) + 
              128*Power(rij,8)*Power(xij,8))))/
      (315.*Power(E,2*rij*(xii + xij))*Power(rij,2)*
        Power(Power(xii,2) - Power(xij,2),15))
    ;
  }
  return S;
}

double DSlater_1S_4S(double rij,double xii,double xij)
{
  return DSlater_4S_1S(rij,xij,xii);
}

double DSlater_2S_4S(double rij,double xii,double xij)
{
  return DSlater_4S_2S(rij,xij,xii);
}

double DSlater_3S_4S(double rij,double xii,double xij)
{
  return DSlater_4S_3S(rij,xij,xii);
}

double Nuclear_1S(double rij,double xii)
{
  double S;

  S = 
  1/rij - (1 + rij*xii)/(Power(E,2*rij*xii)*rij)
    ;
  return S;
}

double Nuclear_2S(double rij,double xii)
{
  double S;

  S = 
  1/rij - (6 + 9*rij*xii + 6*Power(rij,2)*Power(xii,2) + 
       2*Power(rij,3)*Power(xii,3))/(6.*Power(E,2*rij*xii)*rij)
    ;
  return S;
}

double Nuclear_3S(double rij,double xii)
{
  double S;

  S = 
  1/rij - (45 + 75*rij*xii + 60*Power(rij,2)*Power(xii,2) + 
       30*Power(rij,3)*Power(xii,3) + 10*Power(rij,4)*Power(xii,4) + 
       2*Power(rij,5)*Power(xii,5))/(45.*Power(E,2*rij*xii)*rij)
    ;
  return S;
}

double Nuclear_4S(double rij,double xii)
{
  double S;

  S = 
  1/rij - (1260 + 2205*rij*xii + 1890*Power(rij,2)*Power(xii,2) + 
       1050*Power(rij,3)*Power(xii,3) + 420*Power(rij,4)*Power(xii,4) + 
       126*Power(rij,5)*Power(xii,5) + 28*Power(rij,6)*Power(xii,6) + 
       4*Power(rij,7)*Power(xii,7))/(1260.*Power(E,2*rij*xii)*rij)
    ;
  return S;
}

double DNuclear_1S(double rij,double xii)
{
  double S;

  S = 
  Power(rij,-2) - (1 + 2*rij*xii + 2*Power(rij,2)*Power(xii,2))/
     (Power(E,2*rij*xii)*Power(rij,2))
    ;
  return S;
}

double DNuclear_2S(double rij,double xii)
{
  double S;

  S = 
  Power(rij,-2) - (3 + 6*rij*xii + 6*Power(rij,2)*Power(xii,2) + 
       4*Power(rij,3)*Power(xii,3) + 2*Power(rij,4)*Power(xii,4))/
     (3.*Power(E,2*rij*xii)*Power(rij,2))
    ;
  return S;
}

double DNuclear_3S(double rij,double xii)
{
  double S;

  S = 
  Power(rij,-2) - (45 + 90*rij*xii + 90*Power(rij,2)*Power(xii,2) + 
       60*Power(rij,3)*Power(xii,3) + 30*Power(rij,4)*Power(xii,4) + 
       12*Power(rij,5)*Power(xii,5) + 4*Power(rij,6)*Power(xii,6))/
     (45.*Power(E,2*rij*xii)*Power(rij,2))
    ;
  return S;
}

double DNuclear_4S(double rij,double xii)
{
  double S;

  S = 
  Power(rij,-2) - (315 + 630*rij*xii + 630*Power(rij,2)*Power(xii,2) + 
       420*Power(rij,3)*Power(xii,3) + 210*Power(rij,4)*Power(xii,4) + 
       84*Power(rij,5)*Power(xii,5) + 28*Power(rij,6)*Power(xii,6) + 
       8*Power(rij,7)*Power(xii,7) + 2*Power(rij,8)*Power(xii,8))/
     (315.*Power(E,2*rij*xii)*Power(rij,2))
    ;
  return S;
}

t_slater_SS_func (*Slater_SS[SLATER_MAX][SLATER_MAX]) = {
  {  Slater_1S_1S,  Slater_2S_1S,  Slater_3S_1S,  Slater_4S_1S},
  {  Slater_1S_2S,  Slater_2S_2S,  Slater_3S_2S,  Slater_4S_2S},
  {  Slater_1S_3S,  Slater_2S_3S,  Slater_3S_3S,  Slater_4S_3S},
  {  Slater_1S_4S,  Slater_2S_4S,  Slater_3S_4S,  Slater_4S_4S}
};

t_slater_SS_func (*DSlater_SS[SLATER_MAX][SLATER_MAX]) = {
  {  DSlater_1S_1S,  DSlater_2S_1S,  DSlater_3S_1S,  DSlater_4S_1S},
  {  DSlater_1S_2S,  DSlater_2S_2S,  DSlater_3S_2S,  DSlater_4S_2S},
  {  DSlater_1S_3S,  DSlater_2S_3S,  DSlater_3S_3S,  DSlater_4S_3S},
  {  DSlater_1S_4S,  DSlater_2S_4S,  DSlater_3S_4S,  DSlater_4S_4S}
};

t_slater_NS_func (*Slater_NS[SLATER_MAX]) = {
  Nuclear_1S,  Nuclear_2S,  Nuclear_3S,  Nuclear_4S
};

t_slater_NS_func (*DSlater_NS[SLATER_MAX]) = {
  DNuclear_1S,  DNuclear_2S,  DNuclear_3S,  DNuclear_4S
};

