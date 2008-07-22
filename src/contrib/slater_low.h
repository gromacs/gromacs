/* slater_integrals.cpp (c) 2008 Paul J. van Maaren and David van der Spoel */
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

extern cl_R Power(cl_R a,int b);

extern cl_R Slater_1S_1S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_1S_2S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_2S_1S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_1S_3S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_3S_1S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_1S_4S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_4S_1S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_1S_5S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_5S_1S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_1S_6S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_6S_1S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_2S_2S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_2S_3S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_3S_2S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_2S_4S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_4S_2S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_2S_5S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_5S_2S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_2S_6S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_6S_2S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_3S_3S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_3S_4S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_4S_3S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_3S_5S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_5S_3S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_3S_6S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_6S_3S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_4S_4S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_4S_5S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_5S_4S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_4S_6S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_6S_4S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_5S_5S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_5S_6S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_6S_5S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Slater_6S_6S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_1S_1S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_1S_2S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_2S_1S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_1S_3S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_3S_1S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_1S_4S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_4S_1S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_1S_5S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_5S_1S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_1S_6S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_6S_1S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_2S_2S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_2S_3S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_3S_2S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_2S_4S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_4S_2S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_2S_5S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_5S_2S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_2S_6S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_6S_2S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_3S_3S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_3S_4S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_4S_3S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_3S_5S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_5S_3S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_3S_6S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_6S_3S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_4S_4S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_4S_5S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_5S_4S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_4S_6S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_6S_4S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_5S_5S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_5S_6S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_6S_5S(cl_R r,cl_R xi,cl_R xj);

extern cl_R DSlater_6S_6S(cl_R r,cl_R xi,cl_R xj);

extern cl_R Nuclear_1S(cl_R r,cl_R xi);

extern cl_R Nuclear_2S(cl_R r,cl_R xi);

extern cl_R Nuclear_3S(cl_R r,cl_R xi);

extern cl_R Nuclear_4S(cl_R r,cl_R xi);

extern cl_R Nuclear_5S(cl_R r,cl_R xi);

extern cl_R Nuclear_6S(cl_R r,cl_R xi);

extern cl_R DNuclear_1S(cl_R r,cl_R xi);

extern cl_R DNuclear_2S(cl_R r,cl_R xi);

extern cl_R DNuclear_3S(cl_R r,cl_R xi);

extern cl_R DNuclear_4S(cl_R r,cl_R xi);

extern cl_R DNuclear_5S(cl_R r,cl_R xi);

extern cl_R DNuclear_6S(cl_R r,cl_R xi);

