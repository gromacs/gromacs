/*
 *
 * Gromacs 4.0                         Copyright (c) 1991-2003
 * David van der Spoel, Erik Lindahl, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


/*
 * The ia64-assembly Gromacs inner loops would not have been
 * possible without a lot of support and optimization suggestions 
 * from John Worley at Hewlett-Packard.
 */

/* Each thread locks a counter and grabs a couple of neighborlists.
 * Available sizes for this chunk: 1,2,4,8, or 16 
 */
#define THREAD_CHUNK_SIZE       8
#define JJNR_PREFETCH_DISTANCE  128

//	ia64 General Register definitions:	
#define	zero	r0	/* permanent zero					*/
#define	gp		r1	/* global data pointer				*/

#define	at0		r2	/* temp, target of addi				*/
#define	at1		r3	/* temp, target of addi				*/

#define	S0		r4	/* callee saves register			*/
#define	S1		r5	/* callee saves register			*/
#define	S2		r6	/* callee saves register			*/
#define	S3		r7	/* callee saves register			*/

#define	v0		r8	/* 1st fixed point return value/ptr	*/
#define	v1		r9	/* 2nd fixed return value/ptr		*/
#define	v2		r10	/* 3rd fixed return value/ptr		*/
#define	v3		r11	/* 4th fixed return value/ptr		*/

#define	sp		r12	/* memory stack pointer				*/
#define	tp		r13	/* thread pointer					*/

#define	t0		r14	/* caller saves register			*/
#define	t1		r15	/* caller saves register			*/
#define	t2		r16	/* caller saves register			*/
#define	t3		r17	/* caller saves register			*/
#define	t4		r18	/* caller saves register			*/
#define	t5		r19	/* caller saves register			*/

#define	t6		r20	/* caller saves register			*/
#define	t7		r21	/* caller saves register			*/
#define	t8		r22	/* caller saves register			*/
#define	t9		r23	/* caller saves register			*/
#define	t10		r24	/* caller saves register			*/
#define	t11		r25	/* caller saves register			*/
#define	t12		r26	/* caller saves register			*/
#define	t13		r27	/* caller saves register			*/
#define	t14		r28	/* caller saves register			*/
#define	t15		r29	/* caller saves register			*/
#define	t16		r30	/* caller saves register			*/
#define	t17		r31	/* caller saves register			*/


//	ia64 Floating-point register definitions
#define	fZero	f0	/* permanent floating point 0.0		*/
#define	fOne	f1	/* permanent floating point 1.0		*/

#define	fs0		f2	/* callee saves register			*/
#define	fs1		f3	/* callee saves register			*/
#define	fs2		f4	/* callee saves register			*/
#define	fs3		f5	/* callee saves register			*/
	
#define	ft0		f6	/* caller saves register			*/
#define	ft1		f7	/* caller saves register			*/

#define	fa0		f8	/* argument register 0				*/
#define	fa1		f9	/* argument register 1				*/
#define	fa2		f10	/* argument register 2				*/
#define	fa3		f11	/* argument register 3				*/
#define	fa4		f12	/* argument register 4				*/
#define	fa5		f13	/* argument register 5				*/
#define	fa6		f14	/* argument register 6				*/
#define	fa7		f15	/* argument register 7				*/

#define	fv0		f8	/* return value register 0			*/
#define	fv1		f9	/* return value register 1			*/
#define	fv2		f10	/* return value register 2			*/
#define	fv3		f11	/* return value register 3			*/
#define	fv4		f12	/* return value register 4			*/
#define	fv5		f13	/* return value register 5			*/
#define	fv6		f14	/* return value register 6			*/
#define	fv7		f15	/* return value register 7			*/

#define	fs4		f16	/* callee saves register			*/
#define	fs5		f17	/* callee saves register			*/
#define	fs6		f18	/* callee saves register			*/
#define	fs7		f19	/* callee saves register			*/
#define	fs8		f20	/* callee saves register			*/
#define	fs9		f21	/* callee saves register			*/
#define	fs10	f22	/* callee saves register			*/
#define	fs11	f23	/* callee saves register			*/

#define	fs12	f24	/* callee saves register			*/
#define	fs13	f25	/* callee saves register			*/
#define	fs14	f26	/* callee saves register			*/
#define	fs15	f27	/* callee saves register			*/
#define	fs16	f28	/* callee saves register			*/
#define	fs17	f29	/* callee saves register			*/
#define	fs18	f30	/* callee saves register			*/
#define	fs19	f31	/* callee saves register			*/

// ia64 predicate register definitions
#define	pone	p0	/* permanent one predicate			*/
#define	pTrue	p0	/* permanent one predicate			*/

#define	ps0		p1	/* callee saves predicate			*/
#define	ps1		p2	/* callee saves predicate			*/
#define	ps2		p3	/* callee saves predicate			*/
#define	ps3		p4	/* callee saves predicate			*/
#define	ps4		p5	/* callee saves predicate			*/

#define	pt0		p6	/* caller saves predicate			*/
#define	pt1		p7	/* caller saves predicate			*/
#define	pt2		p8	/* caller saves predicate			*/
#define	pt3		p9	/* caller saves predicate			*/
#define	pt4		p10	/* caller saves predicate			*/
#define	pt5		p11	/* caller saves predicate			*/
#define	pt6		p12	/* caller saves predicate			*/
#define	pt7		p13	/* caller saves predicate			*/
#define	pt8		p14	/* caller saves predicate			*/
#define	pt9		p15	/* caller saves predicate			*/

// ia64 branch register definitions
#define	rb		b0	/* return link						*/

#define	bs0		b1	/* callee saves branch register		*/
#define	bs1		b2	/* callee saves branch register		*/
#define	bs2		b3	/* callee saves branch register		*/
#define	bs3		b4	/* callee saves branch register		*/
#define	bs4		b5	/* callee saves branch register		*/
	
#define	bt0		b6	/* caller saves branch register		*/
#define	bt1		b7	/* caller saves branch register		*/
	
		
.text

#define	CHARGE		t10
#define	FACTION		t9
#define	FActII		loc0
#define	FActIX		fs1
#define	FActIY		fs2
#define	FActIZ		fs3
#define	FIX			fs6
#define	FIY			fs7
#define	FIZ			fs8
#define	FSHIFT		t6
#define	FShiftIS	loc1
#define	FShiftX		fa7
#define	FShiftY		fs4
#define	FShiftZ		fs5
#define	ICharge		fa6
#define	InnerCnt	t17
#define	II			t13
#define	II3			in7
#define	IQ			fa6
#define	IS			t12
#define	IS3			in6
#define	IX			fa3
#define	IY			fa4
#define	IZ			fa5
#define	In_FSHIFT	in6
#define	In_GID		in7
#define	In_IINR		in1
#define	In_JINDEX	in2
#define	In_JJNR		in3
#define	In_NRI		in0
#define	In_SHIFT	in4
#define	In_SHIFTVEC	in5
#define NRI			loc4
#define IINR		loc5
#define JINDEX		loc6
#define JJNR		loc7
#define SHIFT		loc8
#define GID			loc9
#define COUNT		loc10
#define	JX			DX[0]
#define	JY			DY[0]
#define	JZ			DZ[0]
#define	LCSave		at0
#define	NJ0			t14
#define	NJ1			t15
#define	POSITION	t8
#define	PRSave		at1
#define	PosX		f88
#define	PosY		f89
#define	PosZ		f90
#define	QCharge		Charge[0]
#define	QQ			Charge[3]
#define	SHIFTVEC	t5
#define	VC			t11
#define	VCPtr		ggid
#define	VCTotal		fs0
#define VNBTotal    fs6
#define	VCoul		Charge[3]
#define	VCoul2		Charge[4]
#define R           RT[0]
#define Vvdw12		C12[2]
#define RInv12      RInv6[1]
#define	argPtr		t4
#define	chargePtr	v2
#define	Tmp1		t0
#define	Tmp2		t17
#define	Tmp3		loc3
#define	Tmp4		t2
#define Tmp5		t3
#define	fHALF		ft0
#define	f3_8		ft1
#define f5_16		fa0
#define fSIX        fa1
#define fTWELVE     fa2
#define	fillP0		v0
#define	fillP1		v1
#define NN0			t0
#define NN1			loc3
#define	ggid		loc2
#define	gidPtr		t7
#define	iinrPtr		t1
#define	jindexPtr	t2
#define	jjnrPtr		t3
#define	jnr			t16
#define	jnr3		v0
#define	nriCount	t0
#define	pCont		pt0
#define	pDone		pt1
#define	pJJNR		pt2
#define	pMore		pt3
#define	pLast		pt4
#define	posPtr		v3
#define fTWO        fa3
#define Facel       fa1
#define Tabscale    fa2
#define	shX			fa3
#define	shY			fa4
#define	shZ			fa5
#define	shiftPtr	t4
#define	shiftVPtr	v1
#define	spillPtr	v0
#define	spillPtr2   t0
#define	xPFS		at0
#define VFTab       loc11
#define nnn         loc12
#define Nouter      loc13
#define Ninner      loc14
#define OuterIter   loc15
#define InnerIter   loc16

#define	_NINPUTS	8
#define	_NLOCALS	17
#define	_NOUTPUT	0
#define	_NROTATE	8

	.regstk	8, 17, 0, 8
	.rotr	FActPtr[8]
	.rotf	DX[8], DY[8], DZ[8], FActX[5], FActY[5], FActZ[5], Charge[8], RSqr[3], RInv[7], RInvT[2], RInvU[2], RInvErr[2], FScalar[2], Y[3], F[3], G[3], H[2], RT[3], n0[2], eps[2]
	.rotp	pPipe[8]

#define	PIPE_DEPTH	8

#define	EXP(n)					(0xffff + (n))

#define	POS_STK_OFFSET			0x10
#define	FACTION_STK_OFFSET		0x18
#define	CHARGE_STK_OFFSET		0x20
#define	FACEL_STK_OFFSET		0x28
#define	KRF_STK_OFFSET			0x30
#define	CRF_STK_OFFSET			0x38
#define	VC_STK_OFFSET			0x40
#define	TYPE_STK_OFFSET			0x48
#define	NTYPE_STK_OFFSET		0x50
#define	NBFP_STK_OFFSET			0x58
#define	VNB_STK_OFFSET			0x60
#define	TABSCALE_STK_OFFSET		0x68
#define	VFTAB_STK_OFFSET		0x70
#define	INVSQRTA_STK_OFFSET		0x78
#define	DVDA_STK_OFFSET			0x80
#define	GBTABSCALE_STK_OFFSET		0x88
#define	GBTAB_STK_OFFSET		0x90
#define	NTHREADS_STK_OFFSET		0x98
#define	COUNT_STK_OFFSET		0xA0
#define MTX_STK_OFFSET			0xA8
#define OUTERITER_STK_OFFSET		0xB0
#define INNERITER_STK_OFFSET		0xB8
#define WORK_STK_OFFSET     		0xC0

	.global nb_kernel300_ia64_double
	.proc	nb_kernel300_ia64_double
	.align	32

nb_kernel300_ia64_double:
//	INIT 1
	{	.mmi
		alloc			xPFS = ar.pfs, _NINPUTS, _NLOCALS, _NOUTPUT, _NROTATE
		mov			spillPtr = sp
		mov			Tmp1 = EXP(-1)
	}
	{	.mfi
		nop			0x0
		nop			0x0
		add			argPtr = FACEL_STK_OFFSET, sp
	} ;;
//	INIT 2
	{	.mfi									
		ld8			Ninner = [argPtr], COUNT_STK_OFFSET - FACEL_STK_OFFSET
		nop			0x0
		nop			0x0
	}
	{	.mii
		stf.spill		[spillPtr] = fs0, -16
		nop			0x0
		add			sp = -8 * 16, sp
	} ;;
//	INIT 3	
	{	.mlx									
		stf.spill		[spillPtr] = fs1, -16
		movl			Tmp3 = 0x3ec00000
	} 
	{	.mmi								
		setf.exp		fHALF = Tmp1
		ld8			COUNT = [argPtr], POS_STK_OFFSET - COUNT_STK_OFFSET
		mov			PRSave	= pr
	} ;;
//	INIT 4	
	{	.mlx
		ld8			POSITION = [argPtr], CHARGE_STK_OFFSET - POS_STK_OFFSET
		movl			Tmp4 = 0x3ea00000
	}
	{	.mmi									
		ldfd			Facel = [Ninner]
		stf.spill		[spillPtr] = fs2, -16
		mov			pr.rot	= 0x0
	} ;;
//  INIT 5	
	{	.mmi								
		ld8			CHARGE = [argPtr], FACTION_STK_OFFSET - CHARGE_STK_OFFSET
		setf.d			f3_8     = Tmp3
		mov			SHIFTVEC   = In_SHIFTVEC
	}
   	{	.mmf
		stf.spill		[spillPtr] = fs3, -16
		ld4			NRI = [In_NRI]
		nop			0x0
	} ;;
//  INIT 6	
	{	.mmi
		ld8			FACTION = [argPtr], VC_STK_OFFSET - FACTION_STK_OFFSET
		setf.d			f5_16 = Tmp4
		mov			FSHIFT   = In_FSHIFT
	} 
	{	.mmi
		stf.spill		[spillPtr] = fs4, -16
		nop			0x0
		mov			GID = In_GID
	} ;;
//  INIT 7	
	{	.mmf
		ld8			VC = [argPtr], TABSCALE_STK_OFFSET - VC_STK_OFFSET
		mov			SHIFT = In_SHIFT
		nop			0x0
	}
	{	.mfi
		stf.spill		[spillPtr] = fs5, -16
		nop			0x0
		nop			0x0

	} ;;
//  INIT 8	
	{	.mmf
		ld8			Ninner = [argPtr], VFTAB_STK_OFFSET - TABSCALE_STK_OFFSET
		mov			JJNR = In_JJNR
		fnorm			f3_8 = f3_8
	}	
	{	.mii
		stf.spill		[spillPtr] = fs6, -16
		mov			IINR = In_IINR
		mov			JINDEX = In_JINDEX
	} ;;
//  INIT 9		
	{	.mmf
		ld8			VFTab = [argPtr], OUTERITER_STK_OFFSET - VFTAB_STK_OFFSET
		stf.spill		[spillPtr] = fs7, -16
		fnorm			f5_16 = f5_16
	} ;; 
//  INIT 10
	{	.mfi
		stf.spill		[spillPtr] = fs8, -16
		fnorm			fHALF = fHALF
		mov			LCSave = ar.lc
	} ;;
//  INIT 11
	{	.mmi
		ld8			OuterIter = [argPtr], INNERITER_STK_OFFSET - OUTERITER_STK_OFFSET
		ldfd			Tabscale = [Ninner]
		mov			Nouter = 0
	} ;;
//  INIT 12
	{	.mfi
		ld8			InnerIter = [argPtr]
		fnorm			Facel = Facel
		mov			Ninner = 0
	} ;;
//  20 bundles used for init - still aligned.




	
threadLoop:
//  THREAD PROLOGUE 1	
	{	.mfi		
		fetchadd4.rel	NN0 = [COUNT], THREAD_CHUNK_SIZE
		nop				0x0
		nop				0x0
	} 
	{	.mfi		// alignment bundle
		nop				0x0
		nop				0x0
		nop				0x0
	} ;;

//  THREAD PROLOGUE 2 - at least 12 cycle latency hole before this bundle (fetchadd4)
	{	.mmi		
		cmp.lt			pCont, pDone = NN0, NRI
		shladd			gidPtr = NN0, 2, GID
		adds			NN1 = THREAD_CHUNK_SIZE, NN0
	}
	{	.mmi
		shladd			jindexPtr = NN0, 2, JINDEX
		shladd   		shiftPtr  = NN0, 2, SHIFT
		shladd			iinrPtr   = NN0, 2, IINR
	} ;; 
//  THREAD PROLOGUE 3 	
	{ .mmi				
	(pCont) ld4			II = [iinrPtr], 4
	(pCont) ld4			IS = [shiftPtr], 4
		cmp.ge			pLast, pMore = NN1, NRI
	}
	{ .mib
	(pCont) ld4			NJ0 = [jindexPtr], 4
	(pCont) adds		Tmp2 = 1, NN0
	(pDone) br.cond.spnt.few finish
	} ;; 		
//  THREAD PROLOGUE 4	
	{ .mmi				
		ld4				ggid = [gidPtr], 4
		shladd			II3 = II, 1, II
		shladd			IS3 = IS, 1, IS
	}
	{ .mmi
		ld4				NJ1 = [jindexPtr], 4
		shladd			chargePtr = II, 3, CHARGE
		shladd			jjnrPtr = NJ0, 2, JJNR
	} ;;
//  THREAD PROLOGUE 5	
	{ .mmi
		cmp.lt			pCont, pDone = Tmp2, NRI						
		shladd			FShiftIS  = IS3, 3, FSHIFT
		nop				0x0
	}	
	{ .mmi
		shladd			posPtr    = II3, 3, POSITION
		shladd			FActII    = II3, 3, FACTION
		shladd			shiftVPtr = IS3, 3, SHIFTVEC	
	} ;;
//  THREAD PROLOGUE 6	
	{ .mmi				
		ld4				jnr = [jjnrPtr], 4
	(pCont)	ld4			IS = [shiftPtr], 4
		nop				0x0
	}	
	{ .mmi
	(pCont) ld4			II = [iinrPtr], 4
			nop			0x0
	(pLast)	mov			NN1 = NRI
	} ;;
//  12 bundles in thread prologue - still aligned





outerLoop:
	//	At this point in the outer loop, the following values are ready
	//
	//		FActII		Pointer to FACTION XYZ for II
	//		FShiftIS	Pointer to FSHIFT XYZ for IS
	//		shiftVPtr	Pointer to current shift XYZ values
	//		posPtr		Pointer to current XYZ position
	//		chargePtr	Pointer to current atom charge
	//		ggid		Index for Vc array
	//		jjnr		Pointer to next neighbor index
	//		jnr			Current jnr value
	//		NJ0, NJ1	Bounds of current neighbor list
	//
	//	Load up all the floating-point values (yes, McKinley can do 4 FP loads
	//	per cycle) and initialize the loop counters and predicates. Compute
	//	the initial position <x, y, z> and charge. If this isn't the last time
	//	through the loop, start loading the next value for NJ1 - we already
	//	moved the previous NJ1 -> NJ0.
//	OUTER PROLOGUE 1
	{	.mfi						
		nop 		0x0
		mov			FIX = f0
		add		Nouter = 1, Nouter
	}
	{	.mmf
		ldfd		shX = [shiftVPtr], 8
		ldfd		PosX = [posPtr], 8
		mov			FIY = f0
	} ;;
//	OUTER PROLOGUE 2
	{	.mmf
		nop			0x0
		ldfd		shY = [shiftVPtr], 8
		nop			0x0
	}
	{	.mfi
		ldfd		PosY = [posPtr], 8
		nop			0x0
		nop			0x0		
	} ;;
//	OUTER PROLOGUE 3
	{	.mmf						
		ldfd		shZ = [shiftVPtr]
		ldfd		PosZ = [posPtr]
		mov			FIZ = f0
	}
	{	.mmi
		ldfd		FShiftX = [FShiftIS], 8
		ldfd		FActIX = [FActII], 8
		nop			0x0
	} ;;
//	OUTER PROLOGUE 4
	{	.mmf	
		ldfd		FShiftY = [FShiftIS], 8
		ldfd		FActIY = [FActII], 8
		nop			0x0
	}
	{ 	.mii
		shladd		VCPtr = ggid, 3, VC
		sub			InnerCnt = NJ1, NJ0, 1
		mov			NJ0 = NJ1
	} ;;
//	OUTER PROLOGUE 5
	{	.mmi
		ldfd		FActIZ = [FActII], -16
		ldfd		FShiftZ = [FShiftIS], -16
		nop			0x0
	} ;;
//	OUTER PROLOGUE 6
	{	.mmf		
		ldfd		ICharge = [chargePtr]
		nop			0x0
		fadd		IX = shX, PosX
	} ;;
//	OUTER PROLOGUE 7
	{	.mfi
		ldfd		VCTotal = [VCPtr]
		fadd		IY = shY, PosY
		add			NN0 = 1, NN0
	}
	{	.mfi
	(pCont)	ld4		NJ1 = [jindexPtr], 4
		nop			0x0

		//	This may seem strange, but we set the first stage of the
		//	pipe to execute this way because setting pr.rot doesn't take
		//	into account how much the predicates have rotated. If this is
		//	the first time through, we cleared all the pipeline predicates
		//	in the initialization. If not, flushing the pipeline set all
		//	the pipeline predicates to 0

		cmp.eq		pPipe[0], p0 = zero, zero
	} ;;
//	OUTER PROLOGUE 8
	{	.mfi		
		cmp.lt		pCont, pDone = NN0, NN1
		fadd		IZ = shZ, PosZ
		mov		    ar.lc = InnerCnt
	} ;;
//	OUTER PROLOGUE 9
	{	.mfi		
		nop			0x0
		fmpy		IQ = ICharge, Facel
		mov			ar.ec = PIPE_DEPTH
	} ;;
// 14 bundles in outer loop - still aligned.

	//	The inner loop is a 6-stage pipeline. The serial sequence of float ops
	//	is folded into a 17-cycle loop (17 * 2 = 34 float ops, one empty), 
    //  then divided
	//	into 5 stages. 



innerLoop:
//	INNER LOOP 1
	{	.mfi	
	(pPipe[5])	ldfpd	Y[0], F[0] = [nnn], 16
	(pPipe[1])	fsub	DX[1] = IX, DX[1]
	(pPipe[0])	shladd	jnr3 = jnr, 1, jnr
	}
	//	We march through jjnr[] sequentially, so it's usually a good idea
	//	to preload the next value. However, we don't want to do this if
	//	(1) we're in the epilogue or (2) this is the last time through and
	//	there are no more atoms to inspect. Thus, we keep track of the loop
	//	trip and use the logic below to see if we should load ahead

	.pred.rel "mutex", pCont, pDone

	{	.mfi
	(pCont)		cmp.ge	pJJNR, p0 = InnerCnt, zero
	(pPipe[3])	fmpy	RInvErr[1] = RSqr[2], RInv[2]
	(pDone)		cmp.gt	pJJNR, p0 = InnerCnt, zero
	} ;;
//	INNER LOOP 2
	{	.mfi	
	(pPipe[5])	ldfpd	G[0], H[0] = [nnn]
	(pPipe[1])	fsub	DY[1] = IY, DY[1]
	(pPipe[0])	add		InnerCnt = -1, InnerCnt
	}
	{	.mfi
	(pPipe[0])	shladd	posPtr = jnr3, 3, POSITION
	(pPipe[7])	fnma	F[2] = F[2], RInv[6], fZero
	(pPipe[0])	shladd	FActPtr[0] = jnr3, 3, FACTION
	} ;;
//	INNER LOOP 3
	{	.mfi									
	(pPipe[0])	ldfd	JX = [posPtr], 8
	(pPipe[1])	fsub	DZ[1] = IZ, DZ[1]
	(pPipe[0])	shladd	chargePtr = jnr, 3, CHARGE
	}
	{  	.mfi
				nop		0x0
	(pPipe[2])	fmpy	RInvErr[0] = RInv[1], RSqr[1]
				nop		0x0
	} ;;
//	INNER LOOP 4
	{	.mfi	
	(pPipe[0])	ldfd	JY = [posPtr], 8
	(pPipe[3])	fmpy	RInvT[1] = fHALF, RInv[2]
				nop		0x0
	}
	{	.mfi
	(pJJNR)		ld4		jnr = [jjnrPtr], 4
	(pPipe[4])	fcvt.fx.trunc n0[0] = RT[1]
	(pPipe[0])	add	Ninner = 1, Ninner
	} ;;
//	INNER LOOP 5
	{	.mfi									
	(pPipe[0])	ldfd	JZ = [posPtr], 8
	(pPipe[1])	fmpy	RSqr[0] = DX[1], DX[1]
				nop		0x0
	}
	{	.mfi
				nop		0x0
	(pPipe[3])	fnma	RInvErr[1] = RInvErr[1], RInv[2], fOne
				nop		0x0
	} ;;
//	INNER LOOP 6
	{	.mfi	
	(pPipe[3])	ldfd	FActX[0] = [FActPtr[3]], 8
	(pPipe[7])	fnma	FActX[4] = F[2], DX[7], FActX[4]	
	(pJJNR)     add     jjnrPtr = JJNR_PREFETCH_DISTANCE, jjnrPtr
	}
	{	.mfi
				nop		0x0
	(pPipe[7])	fnma	FActY[4] = F[2], DY[7], FActY[4]	
				nop		0x0
	} ;;
//	INNER LOOP 7
	{	.mfi	
	(pPipe[3])	ldfd	FActY[0] = [FActPtr[3]], 8
	(pPipe[2])	fnma	RInvErr[0] = RInvErr[0], RInv[1], fOne
				nop		0x0
	}
	{	.mfi
				nop		0x0
	(pPipe[6])	fma		G[1] = eps[1], H[1], G[1]	
				nop		0x0
	} ;;

//	INNER LOOP 8
	{	.mfi	
	(pPipe[3])	ldfd	FActZ[0] = [FActPtr[3]], -16
	(pPipe[7])	fnma	FActZ[4] = F[2], DZ[7], FActZ[4]	
				nop		0x0
	}
	{	.mfi
	(pPipe[4])	getf.sig	nnn = n0[0]
	(pPipe[4])	fcvt.xf 	n0[0] = n0[0]
				nop			0x0
	} ;;
//	INNER LOOP 9
	{	.mfi	
	(pPipe[0])	ldfd	Charge[0] = [chargePtr]	
	(pPipe[1])	fma		RSqr[0] = DY[1], DY[1], RSqr[0]
				nop		0x0
	}
	{	.mfi
	(pJJNR)     lfetch.nta  [jjnrPtr]
	(pPipe[3])	fma		RInv[2] = RInvErr[1], RInvT[1], RInv[2]
				nop		0x0
	} ;;
//	INNER LOOP 10
	{	.mfi	
				nop		0x0
	(pPipe[7])	fma 	FIX = DX[7], F[2], FIX
				nop		0x0
	}
	{	.mfi
				nop		0x0
	(pPipe[7])	fma 	FIY = DY[7], F[2], FIY
				nop		0x0
	} ;;
//	INNER LOOP 11
	{	.mfi	
				nop		0x0
	(pPipe[2])	fma		RInvT[0] = RInvErr[0], f3_8, fHALF
	(pJJNR)     add     jjnrPtr = -JJNR_PREFETCH_DISTANCE, jjnrPtr
	}
	{	.mfi
				nop		0x0
	(pPipe[2])	fmpy	RInvU[0] = RInv[1], RInvErr[0]
				nop		0x0
	} ;;
//	INNER LOOP 12
	{	.mfi	
				nop		0x0	
	(pPipe[6])	fma		F[1] = eps[1], G[1], F[1]
				nop		0x0
	}
	{	.mfi
				nop		0x0
	(pPipe[6])	fma		G[1] = eps[1], H[1], G[1]	
				nop		0x0
	} ;;
//	INNER LOOP 13
	{	.mfi	
				nop		0x0	
	(pPipe[1])	fma		RSqr[0] = DZ[1], DZ[1], RSqr[0]
				nop		0x0
	}
	{	.mfi
				nop		0x0
	(pPipe[3])	fmpy	RInv[2] = RInv[2], Tabscale
				nop		0x0
	} ;;
//	INNER LOOP 14
	{	.mfi	
				nop		0x0
	(pPipe[5])	fmpy	RInv[4] = RInv[4], Charge[5]
				nop		0x0
	}
	{	.mfi
				nop		0x0
	(pPipe[7])	fma 	FIZ = DZ[7], F[2], FIZ
				nop		0x0
	} ;;
//	INNER LOOP 15
	{	.mfi	
				nop		0x0
	(pPipe[2])	fma		RInv[1] = RInvU[0], RInvT[0], RInv[1]
				nop		0x0
	}
	{	.mfi
				nop		0x0
	(pPipe[5])	fsub  	eps[0] = RT[2], n0[1]
				nop		0x0
	} ;;
//	INNER LOOP 16
	{	.mfi	
	(pPipe[7])	stfd	[FActPtr[7]] = FActX[4], 8
	(pPipe[6])	fma     Y[1] = eps[1], F[1], Y[1]
				nop		0x0
	}
	{	.mfi
				nop		0x0
	(pPipe[6])	fma		F[1] = eps[1], G[1], F[1]
	(pPipe[4])	shladd 	nnn = nnn, 3, zero
	} ;;
//	INNER LOOP 17
	{	.mfi	
	(pPipe[7])	stfd	[FActPtr[7]] = FActY[4], 8
	(pPipe[1])	frsqrta RInv[0], p0 = RSqr[0]
				nop		0x0
	}
	{	.mfi
				nop		0x0
	(pPipe[3])	fmpy	RT[0] = RInv[2], RSqr[2]
	(pPipe[4])	shladd 	nnn = nnn, 2, VFTab
	} ;;
//	INNER LOOP 18
	{	.mfi	
	(pPipe[7])	stfd	[FActPtr[7]] = FActZ[4]
	(pPipe[7])	fma		VCTotal = Charge[7], Y[2], VCTotal	
			nop		0x0
	}
	{	.mfb
				nop		0x0
	(pPipe[2])	fmpy	Charge[2] = Charge[2], IQ
			br.ctop.sptk.many	innerLoop
	} ;;


// 	End of modulo-scheduled inner loop

	//	Having finshed the loop, we now compute various quantities to
	//	store. In paralllel, start computing computing some of the values
	//	for the next loop trip, if we're going there.
//	OUTER EPILOGUE 1
    {   .mfi
	  	nop	0x0
		nop	0x0
	(pCont)	shladd	II3 = II, 1, II
    }
	{	.mfi								
	(pCont)	shladd	chargePtr = II, 3, CHARGE
		nop 0x0
	(pCont)	shladd	IS3 = IS, 1, IS
    } ;;
//	OUTER EPILOGUE 2
    {   .mfi
	(pCont)	ld4	II = [iinrPtr] ,4
		fadd	FActIX = FActIX, FIX
		nop 0x0
	}
    {   .mfi
	(pCont)	ld4	IS = [shiftPtr], 4
		fadd	FShiftX = FShiftX, FIX
		nop 0x0
	} ;;
// 	OUTER EPILOGUE 3
    {   .mfi
		nop 0x0
		fadd	FActIY = FActIY, FIY
	(pCont)	shladd	shiftVPtr = IS3, 3, SHIFTVEC						
	} 
    {   .mfi
		nop 0x0
		fadd	FShiftY = FShiftY, FIY
	(pCont)	shladd	posPtr = II3, 3, POSITION
	} ;;
//	OUTER EPILOGUE 4
    {   .mfi
		nop 	0x0
		fadd	FActIZ = FActIZ, FIZ
		nop 	0x0
	} 
    {   .mfi
		nop 	0x0
		fadd	FShiftZ = FShiftZ, FIZ		
		nop 	0x0
	} ;;
//	OUTER EPILOGUE 5
	{	.mmi
		stfd	[FActII] = FActIX, 8
		stfd	[FShiftIS] = FShiftX, 8
		nop 	0x0
	}
    {   .mmi
		stfd    [VCPtr] = VCTotal
	(pCont)		ld4     ggid = [gidPtr], 4 
		nop 	0x0
	} ;;
//	OUTER EPILOGUE 6
	{	.mmi
		stfd	[FActII] = FActIY, 8
		stfd	[FShiftIS] = FShiftY, 8
		nop 0x0
	} ;;
//	OUTER EPILOGUE 7
	{	.mfi
		stfd	[FActII] = FActIZ
		nop		0x0
	(pCont)	shladd	FActII = II3, 3, FACTION
	}
	{	.mib
		stfd	[FShiftIS] = FShiftZ
	(pCont)	shladd	FShiftIS = IS3, 3, FSHIFT
	(pCont)	br.cond.sptk.many	outerLoop
	} ;;




	// Finish if this was the last chunk, or do another thread-loop iteration
//  THREAD EPILOGUE 1
	{ .mib				
		nop				0x0
		nop				0x0
	(pMore) br.cond.sptk.many threadLoop
	} ;;
	



	//	Ready to exit - restore the floating-point registers we saved, the
	//	loop counter, and the predicates, then we're done. Note that the
	//	stack pointer has the address of the last saved FP register.

finish:
//  EXIT 1
	{	.mmi
		mov			fillP0 = sp
		add			fillP1 = 16, sp
		nop			0x0
	}  
	{	.mmi
		st4			[OuterIter] = Nouter
		st4			[InnerIter] = Ninner
		nop			0x0
	} ;;
//  EXIT 2
	{	.mmi
		ldf.fill		fs8 = [fillP0], 32
		ldf.fill		fs7 = [fillP1], 32
		nop				0x0
	} ;;
//  EXIT 3
	{	.mmi
		ldf.fill		fs6 = [fillP0], 32
		ldf.fill		fs5 = [fillP1], 32
		mov				ar.lc = LCSave
	} ;;
//  EXIT 4
	{	.mmi
		ldf.fill		fs4 = [fillP0], 32
		ldf.fill		fs3 = [fillP1], 32
		mov				pr = PRSave, 0x1ffff
	} ;;
//  EXIT 5
	{	.mmi
		ldf.fill		fs2 = [fillP0], 32
		ldf.fill		fs1 = [fillP1], 32
		add				sp = 8 * 16, sp
	} ;;
//  EXIT 6
	{	.mmb
		ldf.fill		fs0 = [fillP0]
		nop				0x0
		br.ret.sptk.few	rp
	} ;;

	.endp	 nb_kernel300_ia64_double



