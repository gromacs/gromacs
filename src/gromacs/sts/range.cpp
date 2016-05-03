#include "range.h"

Ratio operator*(Ratio r, int n) { return Ratio(r.nom_*n,r.denom_); }
Ratio operator*(int n, Ratio r) { return r*n; }
