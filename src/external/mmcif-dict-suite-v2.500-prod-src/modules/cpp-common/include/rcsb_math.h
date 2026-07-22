//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$

#ifndef RCSB_MATH_H
#define RCSB_MATH_H

// IRIX does not support cmath
#ifdef IRIX_OS
#include <math.h>
#else
#include <cmath>
#endif

#endif // RCSB_MATH_H

