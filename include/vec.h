/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
 * University of Groningen, The Netherlands
 *
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Good gRace! Old Maple Actually Chews Slate
 */
#include "typedefs.h"

#ifndef	_vec_h
#define	_vec_h

#ifdef HAVE_IDENT
#ident	"@(#) vec.h 1.8 12/16/92"
#endif /* HAVE_IDENT */

#ifdef CPLUSPLUS
extern "C" {
#endif

#ifndef DOUBLE
#ifdef C_INVSQRT
extern float invsqrt(float x);
#define INVSQRT_DONE
#endif
#endif

#ifndef INVSQRT_DONE
#define invsqrt(x) (1.0f/sqrt(x))
#endif

extern void rvec_add(rvec a,rvec b,rvec c);
extern void rvec_inc(rvec a,rvec b);
extern void rvec_sub(rvec a,rvec b,rvec c);
extern void rvec_dec(rvec a,rvec b);
extern void copy_rvec(rvec a,rvec b);
extern void copy_mat(matrix a,matrix b);
extern void svmul(real a,rvec v1,rvec v2);
extern void clear_rvec(rvec a);
extern void clear_mat(matrix a);
extern real iprod(rvec a,rvec b);
extern real norm2(rvec a);
extern real norm(rvec a);
extern real cos_angle(rvec a,rvec b);
extern void oprod(rvec a,rvec b,rvec c);
extern void mmul(matrix a,matrix b,matrix dest);
extern real det(matrix a);
extern void m_add(matrix a,matrix b,matrix dest);
extern void m_sub(matrix a,matrix b,matrix dest);
extern void msmul(matrix m1,real r1,matrix dest);
extern void m_inv(matrix src,matrix dest);
extern void mvmul(matrix a,rvec src,rvec dest);
extern void unitv(rvec src,rvec dest);
extern real trace(matrix m);
extern real _divide(real a,real b,char *file,int line);
extern int  _mod(int a,int b,char *file,int line);
#define divide(a,b) _divide((a),(b),__FILE__,__LINE__)
#define mod(a,b)    _mod((a),(b),__FILE__,__LINE__)

#ifdef CPLUSPLUS
}
#endif

#endif	/* _vec_h */
