#include <math.h>
#include "../gmx_blas.h"
#include "../gmx_lapack.h"

void 
F77_FUNC(slasd8,SLASD8)(int *icompq, 
	int *k, 
	float *d__, 
     	float *z__, 
	float *vf, 
	float *vl, 
	float *difl, 
	float *difr, 
	int *lddifr, 
	float *dsigma, 
	float *work, 
	int *info)
{
    int difr_dim1, difr_offset, i__1, i__2;
    float d__2;
    float *p1,*p2,t1,t2;

    int i__, j;
    float dj, rho;
    int iwk1, iwk2, iwk3;
    float temp;
    int iwk2i, iwk3i;
    float diflj, difrj, dsigj;
    float dsigjp;
    int c__1 = 1;
    int c__0 = 0;
    float one = 1.;

    /* avoid warnings on high gcc optimization levels */
    difrj = dsigjp = 0;

     --d__;
    --z__;
    --vf;
    --vl;
    --difl;
    difr_dim1 = *lddifr;
    difr_offset = 1 + difr_dim1;
    difr -= difr_offset;
    --dsigma;
    --work;

    *info = 0;

    p1 = &t1;
    p2 = &t2;

    if (*k == 1) {
	d__[1] = fabs(z__[1]);
	difl[1] = d__[1];
	if (*icompq == 1) {
	    difl[2] = 1.;
	    difr[(difr_dim1 << 1) + 1] = 1.;
	}
	return;
    }

    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
      t1 = dsigma[i__];
      t2 = dsigma[i__];
      /* force store and reload from memory */
      d__2 = (*p1) + (*p2) - dsigma[i__];
    }

    iwk1 = 1;
    iwk2 = iwk1 + *k;
    iwk3 = iwk2 + *k;
    iwk2i = iwk2 - 1;
    iwk3i = iwk3 - 1;

    rho = F77_FUNC(snrm2,SNRM2)(k, &z__[1], &c__1);
    F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &rho, &one, k, &c__1, &z__[1], k, info);
    rho *= rho;

    F77_FUNC(slaset,SLASET)("A", k, &c__1, &one, &one, &work[iwk3], k);

    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	F77_FUNC(slasd4,SLASD4)(k, &j, &dsigma[1], &z__[1], &work[iwk1], &rho, &d__[j], &work[
		iwk2], info);

	if (*info != 0) {
	    return;
	}
	work[iwk3i + j] = work[iwk3i + j] * work[j] * work[iwk2i + j];
	difl[j] = -work[j];
	difr[j + difr_dim1] = -work[j + 1];
	i__2 = j - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[iwk3i + i__] = work[iwk3i + i__] * work[i__] * work[iwk2i + 
		    i__] / (dsigma[i__] - dsigma[j]) / (dsigma[i__] + dsigma[
		    j]);
	}
	i__2 = *k;
	for (i__ = j + 1; i__ <= i__2; ++i__) {
	    work[iwk3i + i__] = work[iwk3i + i__] * work[i__] * work[iwk2i + 
		    i__] / (dsigma[i__] - dsigma[j]) / (dsigma[i__] + dsigma[
		    j]);
	}
    }

    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__2 = sqrt(fabs(work[iwk3i + i__]));
	z__[i__] = (z__[i__] > 0) ? d__2 : -d__2;
    }

    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	diflj = difl[j];
	dj = d__[j];
	dsigj = -dsigma[j];
	if (j < *k) {
	    difrj = -difr[j + difr_dim1];
	    dsigjp = -dsigma[j + 1];
	}
	work[j] = -z__[j] / diflj / (dsigma[j] + dj);
	i__2 = j - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	  t1 = dsigma[i__];
	  t2 = dsigj;
	  /* force store and reload from memory */
	  t1 = (*p1) + (*p2) - diflj;
	  work[i__] = z__[i__] / t1 / ( dsigma[i__] + dj);
	}
	i__2 = *k;
	for (i__ = j + 1; i__ <= i__2; ++i__) {
	  t1 = dsigma[i__];
	  t2 = dsigjp;
	  /* force store and reload from memory */
	  t1 = (*p1) + (*p2) - difrj;
	    work[i__] = z__[i__] / t1 / (dsigma[i__] + dj);
	}
	temp = F77_FUNC(snrm2,SNRM2)(k, &work[1], &c__1);
	work[iwk2i + j] = F77_FUNC(sdot,SDOT)(k, &work[1], &c__1, &vf[1], &c__1) / temp;
	work[iwk3i + j] = F77_FUNC(sdot,SDOT)(k, &work[1], &c__1, &vl[1], &c__1) / temp;
	if (*icompq == 1) {
	    difr[j + (difr_dim1 << 1)] = temp;
	}
    }

    F77_FUNC(scopy,SCOPY)(k, &work[iwk2], &c__1, &vf[1], &c__1);
    F77_FUNC(scopy,SCOPY)(k, &work[iwk3], &c__1, &vl[1], &c__1);

    return;

} 
