/* ilaenv.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "../gmx_blas.h"
#include "../gmx_lapack.h"
#include <string.h>

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))


#ifdef KR_headers
VOID s_copy(a, b, la, lb) register char *a, *b; int la, lb;
#else
void s_copy(register char *a, register char *b, int la, int lb)
#endif
{
	register char *aend, *bend;

	aend = a + la;

	if(la <= lb)
#ifndef NO_OVERWRITE
		if (a <= b || a >= b + la)
#endif
			while(a < aend)
				*a++ = *b++;
#ifndef NO_OVERWRITE
		else
			for(b += la; a < aend; )
				*--aend = *--b;
#endif

	else {
		bend = b + lb;
#ifndef NO_OVERWRITE
		if (a <= b || a >= bend)
#endif
			while(b < bend)
				*a++ = *b++;
#ifndef NO_OVERWRITE
		else {
			a += lb;
			while(b < bend)
				*--a = *--bend;
			a += lb;
			}
#endif
		while(a < aend)
			*a++ = ' ';
		}
	}
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* compare two strings */

#ifdef KR_headers
int s_cmp(a0, b0, la, lb) char *a0, *b0; int la, lb;
#else
int s_cmp(char *a0, char *b0, int la, int lb)
#endif
{
register unsigned char *a, *aend, *b, *bend;
a = (unsigned char *)a0;
b = (unsigned char *)b0;
aend = a + la;
bend = b + lb;

if(la <= lb)
	{
	while(a < aend)
		if(*a != *b)
			return( *a - *b );
		else
			{ ++a; ++b; }

	while(b < bend)
		if(*b != ' ')
			return( ' ' - *b );
		else	++b;
	}

else
	{
	while(b < bend)
		if(*a == *b)
			{ ++a; ++b; }
		else
			return( *a - *b );
	while(a < aend)
		if(*a != ' ')
			return(*a - ' ');
		else	++a;
	}
return(0);
}
#ifdef __cplusplus
}
#endif

int ieeeck_(int *ispec, float *zero, float *one)
{
    /* System generated locals */
    int ret_val;

    /* Local variables */
    float nan1, nan2, nan3, nan4, nan5, nan6, neginf, posinf, negzro, newzro;


/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  IEEECK is called from the ILAENV to verify that Infinity and */
/*  possibly NaN arithmetic is safe (i.e. will not trap). */

/*  Arguments */
/*  ========= */

/*  ISPEC   (input) INTEGER */
/*          Specifies whether to test just for inifinity arithmetic */
/*          or whether to test for infinity and NaN arithmetic. */
/*          = 0: Verify infinity arithmetic only. */
/*          = 1: Verify infinity and NaN arithmetic. */

/*  ZERO    (input) REAL */
/*          Must contain the value 0.0 */
/*          This is passed to prevent the compiler from optimizing */
/*          away this code. */

/*  ONE     (input) REAL */
/*          Must contain the value 1.0 */
/*          This is passed to prevent the compiler from optimizing */
/*          away this code. */

/*  RETURN VALUE:  INTEGER */
/*          = 0:  Arithmetic failed to produce the correct answers */
/*          = 1:  Arithmetic produced the correct answers */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */
    ret_val = 1;

    posinf = *one / *zero;
    if (posinf <= *one) {
	ret_val = 0;
	return ret_val;
    }

    neginf = -(*one) / *zero;
    if (neginf >= *zero) {
	ret_val = 0;
	return ret_val;
    }

    negzro = *one / (neginf + *one);
    if (negzro != *zero) {
	ret_val = 0;
	return ret_val;
    }

    neginf = *one / negzro;
    if (neginf >= *zero) {
	ret_val = 0;
	return ret_val;
    }

    newzro = negzro + *zero;
    if (newzro != *zero) {
	ret_val = 0;
	return ret_val;
    }

    posinf = *one / newzro;
    if (posinf <= *one) {
	ret_val = 0;
	return ret_val;
    }

    neginf *= posinf;
    if (neginf >= *zero) {
	ret_val = 0;
	return ret_val;
    }

    posinf *= posinf;
    if (posinf <= *one) {
	ret_val = 0;
	return ret_val;
    }




/*     Return if we were only asked to check infinity arithmetic */

    if (*ispec == 0) {
	return ret_val;
    }

    nan1 = posinf + neginf;

    nan2 = posinf / neginf;

    nan3 = posinf / posinf;

    nan4 = posinf * *zero;

    nan5 = neginf * negzro;

    nan6 = nan5 * 0.f;

    if (nan1 == nan1) {
	ret_val = 0;
	return ret_val;
    }

    if (nan2 == nan2) {
	ret_val = 0;
	return ret_val;
    }

    if (nan3 == nan3) {
	ret_val = 0;
	return ret_val;
    }

    if (nan4 == nan4) {
	ret_val = 0;
	return ret_val;
    }

    if (nan5 == nan5) {
	ret_val = 0;
	return ret_val;
    }

    if (nan6 == nan6) {
	ret_val = 0;
	return ret_val;
    }

    return ret_val;
} /* ieeeck_ */



/* Table of constant values */

static int c__1 = 1;
static float c_b163 = 0.f;
static float c_b164 = 1.f;
static int c__0 = 0;

int F77_FUNC(ilaenv,ILAENV) (int *ispec, char *name__, char *opts, int *n1,
	int *n2, int *n3, int *n4)
{
    /* System generated locals */
    int ret_val;

    /* Builtin functions */
    /* Subroutine */ 
    int s_cmp(char *, char *, int, int);

    /* Local variables */
    int i__;
    char c1[1], c2[1], c3[1], c4[1];
    int ic, nb, iz, nx;
    int cname;
    int nbmin;
    int sname;
    extern int ieeeck_(int *, float *, float *);
    char subnam[5];
    extern int iparmq_(int *, char *, char *, int *, int *, 
	    int *, int *);

    int name_len, opts_len;

    name_len = strlen (name__);
    opts_len = strlen (opts);

/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     January 2007 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  ILAENV is called from the LAPACK routines to choose problem-dependent */
/*  parameters for the local environment.  See ISPEC for a description of */
/*  the parameters. */

/*  ILAENV returns an INTEGER */
/*  if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC */
/*  if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value. */

/*  This version provides a set of parameters which should give good, */
/*  but not optimal, performance on many of the currently available */
/*  computers.  Users are encouraged to modify this subroutine to set */
/*  the tuning parameters for their particular machine using the option */
/*  and problem size information in the arguments. */

/*  This routine will not function correctly if it is converted to all */
/*  lower case.  Converting it to all upper case is allowed. */

/*  Arguments */
/*  ========= */

/*  ISPEC   (input) INTEGER */
/*          Specifies the parameter to be returned as the value of */
/*          ILAENV. */
/*          = 1: the optimal blocksize; if this value is 1, an unblocked */
/*               algorithm will give the best performance. */
/*          = 2: the minimum block size for which the block routine */
/*               should be used; if the usable block size is less than */
/*               this value, an unblocked routine should be used. */
/*          = 3: the crossover point (in a block routine, for N less */
/*               than this value, an unblocked routine should be used) */
/*          = 4: the number of shifts, used in the nonsymmetric */
/*               eigenvalue routines (DEPRECATED) */
/*          = 5: the minimum column dimension for blocking to be used; */
/*               rectangular blocks must have dimension at least k by m, */
/*               where k is given by ILAENV(2,...) and m by ILAENV(5,...) */
/*          = 6: the crossover point for the SVD (when reducing an m by n */
/*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds */
/*               this value, a QR factorization is used first to reduce */
/*               the matrix to a triangular form.) */
/*          = 7: the number of processors */
/*          = 8: the crossover point for the multishift QR method */
/*               for nonsymmetric eigenvalue problems (DEPRECATED) */
/*          = 9: maximum size of the subproblems at the bottom of the */
/*               computation tree in the divide-and-conquer algorithm */
/*               (used by xGELSD and xGESDD) */
/*          =10: ieee NaN arithmetic can be trusted not to trap */
/*          =11: infinity arithmetic can be trusted not to trap */
/*          12 <= ISPEC <= 16: */
/*               xHSEQR or one of its subroutines, */
/*               see IPARMQ for detailed explanation */

/*  NAME    (input) CHARACTER*(*) */
/*          The name of the calling subroutine, in either upper case or */
/*          lower case. */

/*  OPTS    (input) CHARACTER*(*) */
/*          The character options to the subroutine NAME, concatenated */
/*          into a single character string.  For example, UPLO = 'U', */
/*          TRANS = 'T', and DIAG = 'N' for a triangular routine would */
/*          be specified as OPTS = 'UTN'. */

/*  N1      (input) INTEGER */
/*  N2      (input) INTEGER */
/*  N3      (input) INTEGER */
/*  N4      (input) INTEGER */
/*          Problem dimensions for the subroutine NAME; these may not all */
/*          be required. */

/*  Further Details */
/*  =============== */

/*  The following conventions have been used when calling ILAENV from the */
/*  LAPACK routines: */
/*  1)  OPTS is a concatenation of all of the character options to */
/*      subroutine NAME, in the same order that they appear in the */
/*      argument list for NAME, even if they are not used in determining */
/*      the value of the parameter specified by ISPEC. */
/*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order */
/*      that they appear in the argument list for NAME.  N1 is used */
/*      first, N2 second, and so on, and unused problem dimensions are */
/*      passed a value of -1. */
/*  3)  The parameter value returned by ILAENV is checked for validity in */
/*      the calling subroutine.  For example, ILAENV is used to retrieve */
/*      the optimal blocksize for STRTRI as follows: */

/*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 ) */
/*      IF( NB.LE.1 ) NB = MAX( 1, N ) */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    switch (*ispec) {
	case 1:  goto L10;
	case 2:  goto L10;
	case 3:  goto L10;
	case 4:  goto L80;
	case 5:  goto L90;
	case 6:  goto L100;
	case 7:  goto L110;
	case 8:  goto L120;
	case 9:  goto L130;
	case 10:  goto L140;
	case 11:  goto L150;
	case 12:  goto L160;
	case 13:  goto L160;
	case 14:  goto L160;
	case 15:  goto L160;
	case 16:  goto L160;
    }

/*     Invalid value for ISPEC */

    ret_val = -1;
    return ret_val;

L10:

/*     Convert NAME to upper case if the first character is lower case. */

    ret_val = 1;
    s_copy(subnam, name__, (int)1, name_len);
    ic = *(unsigned char *)subnam;
    iz = 'Z';
    if (iz == 90 || iz == 122) {

/*        ASCII character set */

	if (ic >= 97 && ic <= 122) {
	    *(unsigned char *)subnam = (char) (ic - 32);
	    for (i__ = 2; i__ <= 6; ++i__) {
		ic = *(unsigned char *)&subnam[i__ - 1];
		if (ic >= 97 && ic <= 122) {
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
		}
/* L20: */
	    }
	}

    } else if (iz == 233 || iz == 169) {

/*        EBCDIC character set */

	if (((ic >= 129) && (ic <= 137)) || ((ic >= 145) && (ic <= 153)) || 
	    ((ic >= 162) && (ic <= 169))) {
	    *(unsigned char *)subnam = (char) (ic + 64);
	    for (i__ = 2; i__ <= 6; ++i__) {
		ic = *(unsigned char *)&subnam[i__ - 1];
		if (((ic >= 129) && (ic <= 137)) || ((ic >= 145) && (ic <= 153)) || ((ic >= 162) && (ic <= 169))) {
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic + 64);
		}
/* L30: */
	    }
	}

    } else if (iz == 218 || iz == 250) {

/*        Prime machines:  ASCII+128 */

	if (ic >= 225 && ic <= 250) {
	    *(unsigned char *)subnam = (char) (ic - 32);
	    for (i__ = 2; i__ <= 6; ++i__) {
		ic = *(unsigned char *)&subnam[i__ - 1];
		if (ic >= 225 && ic <= 250) {
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
		}
/* L40: */
	    }
	}
    }

    *(unsigned char *)c1 = *(unsigned char *)subnam;
    sname = *(unsigned char *)c1 == 'S' || *(unsigned char *)c1 == 'D';
    cname = *(unsigned char *)c1 == 'C' || *(unsigned char *)c1 == 'Z';
    if (! (cname || sname)) {
	return ret_val;
    }
    s_copy(c2, subnam + 1, (int)1, (int)2);
    s_copy(c3, subnam + 3, (int)1, (int)3);
    s_copy(c4, c3 + 1, (int)1, (int)2);

    switch (*ispec) {
	case 1:  goto L50;
	case 2:  goto L60;
	case 3:  goto L70;
    }

L50:

/*     ISPEC = 1:  block size */

/*     In these examples, separate code is provided for setting NB for */
/*     real and complex.  We assume that NB will take the same value in */
/*     single or double precision. */

    nb = 1;

    if (s_cmp(c2, "GE", 1, 2) == 0) {
	if (s_cmp(c3, "TRF", 1, 3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	} else if (s_cmp(c3, "QRF", 1, 3) == 0 || s_cmp(c3,
		"RQF", 1, 3) == 0 || s_cmp(c3, "LQF",
		1, 3) == 0 || s_cmp(c3, "QLF", 1, 3)
		== 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (s_cmp(c3, "HRD", (int)1, (int)3) == 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (s_cmp(c3, "BRD", (int)1, (int)3) == 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (s_cmp(c3, "TRI", (int)1, (int)3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (s_cmp(c2, "PO", (int)1, (int)2) == 0) {
	if (s_cmp(c3, "TRF", (int)1, (int)3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (s_cmp(c2, "SY", (int)1, (int)2) == 0) {
	if (s_cmp(c3, "TRF", (int)1, (int)3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	} else if (sname && s_cmp(c3, "TRD", (int)1, (int)3) == 0) {
	    nb = 32;
	} else if (sname && s_cmp(c3, "GST", (int)1, (int)3) == 0) {
	    nb = 64;
	}
    } else if (cname && s_cmp(c2, "HE", (int)1, (int)2) == 0) {
	if (s_cmp(c3, "TRF", (int)1, (int)3) == 0) {
	    nb = 64;
	} else if (s_cmp(c3, "TRD", (int)1, (int)3) == 0) {
	    nb = 32;
	} else if (s_cmp(c3, "GST", (int)1, (int)3) == 0) {
	    nb = 64;
	}
    } else if (sname && s_cmp(c2, "OR", (int)1, (int)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (int)1, (int)2) == 0 || s_cmp(c4, "RQ", 
		    (int)1, (int)2) == 0 || s_cmp(c4, "LQ", (int)1, (
		    int)2) == 0 || s_cmp(c4, "QL", (int)1, (int)2) ==
		     0 || s_cmp(c4, "HR", (int)1, (int)2) == 0 || s_cmp(
		    c4, "TR", (int)1, (int)2) == 0 || s_cmp(c4, "BR", (
		    int)1, (int)2) == 0) {
		nb = 32;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", (int)1, (int)2) == 0 || s_cmp(c4, "RQ", 
		    (int)1, (int)2) == 0 || s_cmp(c4, "LQ", (int)1, (
		    int)2) == 0 || s_cmp(c4, "QL", (int)1, (int)2) ==
		     0 || s_cmp(c4, "HR", (int)1, (int)2) == 0 || s_cmp(
		    c4, "TR", (int)1, (int)2) == 0 || s_cmp(c4, "BR", (
		    int)1, (int)2) == 0) {
		nb = 32;
	    }
	}
    } else if (cname && s_cmp(c2, "UN", (int)1, (int)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (int)1, (int)2) == 0 || s_cmp(c4, "RQ", 
		    (int)1, (int)2) == 0 || s_cmp(c4, "LQ", (int)1, (
		    int)2) == 0 || s_cmp(c4, "QL", (int)1, (int)2) ==
		     0 || s_cmp(c4, "HR", (int)1, (int)2) == 0 || s_cmp(
		    c4, "TR", (int)1, (int)2) == 0 || s_cmp(c4, "BR", (
		    int)1, (int)2) == 0) {
		nb = 32;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", (int)1, (int)2) == 0 || s_cmp(c4, "RQ", 
		    (int)1, (int)2) == 0 || s_cmp(c4, "LQ", (int)1, (
		    int)2) == 0 || s_cmp(c4, "QL", (int)1, (int)2) ==
		     0 || s_cmp(c4, "HR", (int)1, (int)2) == 0 || s_cmp(
		    c4, "TR", (int)1, (int)2) == 0 || s_cmp(c4, "BR", (
		    int)1, (int)2) == 0) {
		nb = 32;
	    }
	}
    } else if (s_cmp(c2, "GB", (int)1, (int)2) == 0) {
	if (s_cmp(c3, "TRF", (int)1, (int)3) == 0) {
	    if (sname) {
		if (*n4 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    } else {
		if (*n4 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    }
	}
    } else if (s_cmp(c2, "PB", (int)1, (int)2) == 0) {
	if (s_cmp(c3, "TRF", (int)1, (int)3) == 0) {
	    if (sname) {
		if (*n2 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    } else {
		if (*n2 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    }
	}
    } else if (s_cmp(c2, "TR", (int)1, (int)2) == 0) {
	if (s_cmp(c3, "TRI", (int)1, (int)3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (s_cmp(c2, "LA", (int)1, (int)2) == 0) {
	if (s_cmp(c3, "UUM", (int)1, (int)3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (sname && s_cmp(c2, "ST", (int)1, (int)2) == 0) {
	if (s_cmp(c3, "EBZ", (int)1, (int)3) == 0) {
	    nb = 1;
	}
    }
    ret_val = nb;
    return ret_val;

L60:

/*     ISPEC = 2:  minimum block size */

    nbmin = 2;
    if (s_cmp(c2, "GE", (int)1, (int)2) == 0) {
	if (s_cmp(c3, "QRF", (int)1, (int)3) == 0 || s_cmp(c3, "RQF", (
		int)1, (int)3) == 0 || s_cmp(c3, "LQF", (int)1, (
		int)3) == 0 || s_cmp(c3, "QLF", (int)1, (int)3) == 0)
		 {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (s_cmp(c3, "HRD", (int)1, (int)3) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (s_cmp(c3, "BRD", (int)1, (int)3) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (s_cmp(c3, "TRI", (int)1, (int)3) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	}
    } else if (s_cmp(c2, "SY", (int)1, (int)2) == 0) {
	if (s_cmp(c3, "TRF", (int)1, (int)3) == 0) {
	    if (sname) {
		nbmin = 8;
	    } else {
		nbmin = 8;
	    }
	} else if (sname && s_cmp(c3, "TRD", (int)1, (int)3) == 0) {
	    nbmin = 2;
	}
    } else if (cname && s_cmp(c2, "HE", (int)1, (int)2) == 0) {
	if (s_cmp(c3, "TRD", (int)1, (int)3) == 0) {
	    nbmin = 2;
	}
    } else if (sname && s_cmp(c2, "OR", (int)1, (int)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (int)1, (int)2) == 0 || s_cmp(c4, "RQ", 
		    (int)1, (int)2) == 0 || s_cmp(c4, "LQ", (int)1, (
		    int)2) == 0 || s_cmp(c4, "QL", (int)1, (int)2) ==
		     0 || s_cmp(c4, "HR", (int)1, (int)2) == 0 || s_cmp(
		    c4, "TR", (int)1, (int)2) == 0 || s_cmp(c4, "BR", (
		    int)1, (int)2) == 0) {
		nbmin = 2;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", (int)1, (int)2) == 0 || s_cmp(c4, "RQ", 
		    (int)1, (int)2) == 0 || s_cmp(c4, "LQ", (int)1, (
		    int)2) == 0 || s_cmp(c4, "QL", (int)1, (int)2) ==
		     0 || s_cmp(c4, "HR", (int)1, (int)2) == 0 || s_cmp(
		    c4, "TR", (int)1, (int)2) == 0 || s_cmp(c4, "BR", (
		    int)1, (int)2) == 0) {
		nbmin = 2;
	    }
	}
    } else if (cname && s_cmp(c2, "UN", (int)1, (int)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (int)1, (int)2) == 0 || s_cmp(c4, "RQ", 
		    (int)1, (int)2) == 0 || s_cmp(c4, "LQ", (int)1, (
		    int)2) == 0 || s_cmp(c4, "QL", (int)1, (int)2) ==
		     0 || s_cmp(c4, "HR", (int)1, (int)2) == 0 || s_cmp(
		    c4, "TR", (int)1, (int)2) == 0 || s_cmp(c4, "BR", (
		    int)1, (int)2) == 0) {
		nbmin = 2;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", (int)1, (int)2) == 0 || s_cmp(c4, "RQ", 
		    (int)1, (int)2) == 0 || s_cmp(c4, "LQ", (int)1, (
		    int)2) == 0 || s_cmp(c4, "QL", (int)1, (int)2) ==
		     0 || s_cmp(c4, "HR", (int)1, (int)2) == 0 || s_cmp(
		    c4, "TR", (int)1, (int)2) == 0 || s_cmp(c4, "BR", (
		    int)1, (int)2) == 0) {
		nbmin = 2;
	    }
	}
    }
    ret_val = nbmin;
    return ret_val;

L70:

/*     ISPEC = 3:  crossover point */

    nx = 0;
    if (s_cmp(c2, "GE", (int)1, (int)2) == 0) {
	if (s_cmp(c3, "QRF", (int)1, (int)3) == 0 || s_cmp(c3, "RQF", (
		int)1, (int)3) == 0 || s_cmp(c3, "LQF", (int)1, (
		int)3) == 0 || s_cmp(c3, "QLF", (int)1, (int)3) == 0)
		 {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	} else if (s_cmp(c3, "HRD", (int)1, (int)3) == 0) {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	} else if (s_cmp(c3, "BRD", (int)1, (int)3) == 0) {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	}
    } else if (s_cmp(c2, "SY", (int)1, (int)2) == 0) {
	if (sname && s_cmp(c3, "TRD", (int)1, (int)3) == 0) {
	    nx = 32;
	}
    } else if (cname && s_cmp(c2, "HE", (int)1, (int)2) == 0) {
	if (s_cmp(c3, "TRD", (int)1, (int)3) == 0) {
	    nx = 32;
	}
    } else if (sname && s_cmp(c2, "OR", (int)1, (int)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (int)1, (int)2) == 0 || s_cmp(c4, "RQ", 
		    (int)1, (int)2) == 0 || s_cmp(c4, "LQ", (int)1, (
		    int)2) == 0 || s_cmp(c4, "QL", (int)1, (int)2) ==
		     0 || s_cmp(c4, "HR", (int)1, (int)2) == 0 || s_cmp(
		    c4, "TR", (int)1, (int)2) == 0 || s_cmp(c4, "BR", (
		    int)1, (int)2) == 0) {
		nx = 128;
	    }
	}
    } else if (cname && s_cmp(c2, "UN", (int)1, (int)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (int)1, (int)2) == 0 || s_cmp(c4, "RQ", 
		    (int)1, (int)2) == 0 || s_cmp(c4, "LQ", (int)1, (
		    int)2) == 0 || s_cmp(c4, "QL", (int)1, (int)2) ==
		     0 || s_cmp(c4, "HR", (int)1, (int)2) == 0 || s_cmp(
		    c4, "TR", (int)1, (int)2) == 0 || s_cmp(c4, "BR", (
		    int)1, (int)2) == 0) {
		nx = 128;
	    }
	}
    }
    ret_val = nx;
    return ret_val;

L80:

/*     ISPEC = 4:  number of shifts (used by xHSEQR) */

    ret_val = 6;
    return ret_val;

L90:

/*     ISPEC = 5:  minimum column dimension (not used) */

    ret_val = 2;
    return ret_val;

L100:

/*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD) */

    ret_val = (int) ((float) min(*n1,*n2) * 1.6f);
    return ret_val;

L110:

/*     ISPEC = 7:  number of processors (not used) */

    ret_val = 1;
    return ret_val;

L120:

/*     ISPEC = 8:  crossover point for multishift (used by xHSEQR) */

    ret_val = 50;
    return ret_val;

L130:

/*     ISPEC = 9:  maximum size of the subproblems at the bottom of the */
/*                 computation tree in the divide-and-conquer algorithm */
/*                 (used by xGELSD and xGESDD) */

    ret_val = 25;
    return ret_val;

L140:

/*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap */

/*     ILAENV = 0 */
    ret_val = 1;
    if (ret_val == 1) {
	ret_val = ieeeck_(&c__1, &c_b163, &c_b164);
    }
    return ret_val;

L150:

/*     ISPEC = 11: infinity arithmetic can be trusted not to trap */

/*     ILAENV = 0 */
    ret_val = 1;
    if (ret_val == 1) {
	ret_val = ieeeck_(&c__0, &c_b163, &c_b164);
    }
    return ret_val;

L160:

/*     12 <= ISPEC <= 16: xHSEQR or one of its subroutines. */

    ret_val = F77_FUNC(iparmq, IPARMQ) (ispec, name__, opts, n1, n2, n3, n4)
	    ;
    return ret_val;

/*     End of ILAENV */

} /* ilaenv_ */
