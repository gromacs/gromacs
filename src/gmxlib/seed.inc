C     #define	EXP_LSB		frac0x00800000
C     #define	EXP_SEED_SIZE	256
C     #define	EXP_MASK	0x7f800000
C     #define	EXP_SHIFT	23
C     #define	MAX_FRACT	0x007fffff
C     #define	FRACT_MASK	0x007fffff
C     #define	FRACT_SIZE	11              significant part of fraction
C     #define	FRACT_SHIFT	(EXP_SHIFT-FRACT_SIZE)
C     #define	FRACT_SEED_SIZE	(1<<(FRACT_SIZE+1))   one bit for even/odd
C     #define	FRACT_FIRST	(0x3f800000>>FRACT_SHIFT)
C     #define	NOT_INITIALISED	~0
C     #define	EXP_ADDR(val)	(((val)&EXP_MASK)>>EXP_SHIFT)
C     #define	FRACT_ADDR(val)	(((val)&(FRACT_MASK|EXP_LSB))>>FRACT_SHIFT)
      
      implicit none
      integer*4 nexp,nfract,maxfract
      integer*4 fractshift,fractmask,fractf
      integer*4 expshift,expmask,explsb
      integer*4 expseed,fracseed
      parameter (nexp=256,nfract=4096,maxfract=8388607)
      parameter (fractshift=12,fractmask=8388607,fractf=1065353216)
      parameter (expshift=23,expmask=2139095040,explsb=8388608)
      common /seed/ expseed(nexp),fracseed(nfract)
      integer*4 fl2i,fractaddr,expaddr
      real*4    i2fl
