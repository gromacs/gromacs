;#
;#
;# Gromacs 4.0                         Copyright (c) 1991-2003 
;# David van der Spoel, Erik Lindahl
;#
;# This program is free software; you can redistribute it and/or
;# modify it under the terms of the GNU General Public License
;# as published by the Free Software Foundation; either version 2
;# of the License, or (at your option) any later version.
;#
;# To help us fund GROMACS development, we humbly ask that you cite
;# the research papers on the package. Check out http://www.gromacs.org
;# 
;# And Hey:
;# Gnomes, ROck Monsters And Chili Sauce
;#

;# These files require GNU binutils 2.10 or later, since we
;# use intel syntax for portability, or a recent version 
;# of NASM that understands Extended 3DNow and SSE2 instructions.
;# (NASM is normally only used with MS Visual C++).
;# Since NASM and gnu as disagree on some definitions and use 
;# completely different preprocessing options I have to introduce a
;# trick: NASM uses ';' for comments, while gnu as uses '#' on x86.
;# Gnu as treats ';' as a line break, i.e. ignores it. This is the
;# reason why all comments need both symbols...
;# The source is written for GNU as, with intel syntax. When you use
;# NASM we redefine a couple of things. The false if-statement around 
;# the following code is seen by GNU as, but NASM doesn't see it, so 
;# the code inside is read by NASM but not gcc.

; .if 0    # block below only read by NASM
%define .section	section
%define .long		dd
%define .align		align
%define .globl		global
;# NASM only wants 'dword', not 'dword ptr'.
%define ptr
%macro .equiv                  2
   %1 equ %2
%endmacro
; .endif                   # End of NASM-specific block
; .intel_syntax noprefix   # Line only read by gnu as

.section .text

	
.globl nb_kernel334_ia32_sse2
.globl _nb_kernel334_ia32_sse2
nb_kernel334_ia32_sse2:	
_nb_kernel334_ia32_sse2:	
.equiv          nb334_p_nri,            8
.equiv          nb334_iinr,             12
.equiv          nb334_jindex,           16
.equiv          nb334_jjnr,             20
.equiv          nb334_shift,            24
.equiv          nb334_shiftvec,         28
.equiv          nb334_fshift,           32
.equiv          nb334_gid,              36
.equiv          nb334_pos,              40
.equiv          nb334_faction,          44
.equiv          nb334_charge,           48
.equiv          nb334_p_facel,          52
.equiv          nb334_argkrf,           56
.equiv          nb334_argcrf,           60
.equiv          nb334_Vc,               64
.equiv          nb334_type,             68
.equiv          nb334_p_ntype,          72
.equiv          nb334_vdwparam,         76
.equiv          nb334_Vvdw,             80
.equiv          nb334_p_tabscale,       84
.equiv          nb334_VFtab,            88
.equiv          nb334_invsqrta,         92
.equiv          nb334_dvda,             96
.equiv          nb334_p_gbtabscale,     100
.equiv          nb334_GBtab,            104
.equiv          nb334_p_nthreads,       108
.equiv          nb334_count,            112
.equiv          nb334_mtx,              116
.equiv          nb334_outeriter,        120
.equiv          nb334_inneriter,        124
.equiv          nb334_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb334_ixO,              0
.equiv          nb334_iyO,              16
.equiv          nb334_izO,              32
.equiv          nb334_ixH1,             48
.equiv          nb334_iyH1,             64
.equiv          nb334_izH1,             80
.equiv          nb334_ixH2,             96
.equiv          nb334_iyH2,             112
.equiv          nb334_izH2,             128
.equiv          nb334_ixM,              144
.equiv          nb334_iyM,              160
.equiv          nb334_izM,              176
.equiv          nb334_jxO,              192
.equiv          nb334_jyO,              208
.equiv          nb334_jzO,              224
.equiv          nb334_jxH1,             240
.equiv          nb334_jyH1,             256
.equiv          nb334_jzH1,             272
.equiv          nb334_jxH2,             288
.equiv          nb334_jyH2,             304
.equiv          nb334_jzH2,             320
.equiv          nb334_jxM,              336
.equiv          nb334_jyM,              352
.equiv          nb334_jzM,              368
.equiv          nb334_dxOO,             384
.equiv          nb334_dyOO,             400
.equiv          nb334_dzOO,             416
.equiv          nb334_dxH1H1,           432
.equiv          nb334_dyH1H1,           448
.equiv          nb334_dzH1H1,           464
.equiv          nb334_dxH1H2,           480
.equiv          nb334_dyH1H2,           496
.equiv          nb334_dzH1H2,           512
.equiv          nb334_dxH1M,            528
.equiv          nb334_dyH1M,            544
.equiv          nb334_dzH1M,            560
.equiv          nb334_dxH2H1,           576
.equiv          nb334_dyH2H1,           592
.equiv          nb334_dzH2H1,           608
.equiv          nb334_dxH2H2,           624
.equiv          nb334_dyH2H2,           640
.equiv          nb334_dzH2H2,           656
.equiv          nb334_dxH2M,            672
.equiv          nb334_dyH2M,            688
.equiv          nb334_dzH2M,            704
.equiv          nb334_dxMH1,            720
.equiv          nb334_dyMH1,            736
.equiv          nb334_dzMH1,            752
.equiv          nb334_dxMH2,            768
.equiv          nb334_dyMH2,            784
.equiv          nb334_dzMH2,            800
.equiv          nb334_dxMM,             816
.equiv          nb334_dyMM,             832
.equiv          nb334_dzMM,             848
.equiv          nb334_qqMM,             864
.equiv          nb334_qqMH,             880
.equiv          nb334_qqHH,             896
.equiv          nb334_two,              912
.equiv          nb334_tsc,              928
.equiv          nb334_c6,               944
.equiv          nb334_c12,              960
.equiv          nb334_vctot,            976
.equiv          nb334_Vvdwtot,          992
.equiv          nb334_fixO,             1008
.equiv          nb334_fiyO,             1024
.equiv          nb334_fizO,             1040
.equiv          nb334_fixH1,            1056
.equiv          nb334_fiyH1,            1072
.equiv          nb334_fizH1,            1088
.equiv          nb334_fixH2,            1104
.equiv          nb334_fiyH2,            1120
.equiv          nb334_fizH2,            1136
.equiv          nb334_fixM,             1152
.equiv          nb334_fiyM,             1168
.equiv          nb334_fizM,             1184
.equiv          nb334_fjxO,             1200
.equiv          nb334_fjyO,             1216
.equiv          nb334_fjzO,             1232
.equiv          nb334_fjxH1,            1248
.equiv          nb334_fjyH1,            1264
.equiv          nb334_fjzH1,            1280
.equiv          nb334_fjxH2,            1296
.equiv          nb334_fjyH2,            1312
.equiv          nb334_fjzH2,            1328
.equiv          nb334_fjxM,             1344
.equiv          nb334_fjyM,             1360
.equiv          nb334_fjzM,             1376
.equiv          nb334_half,             1392
.equiv          nb334_three,            1408
.equiv          nb334_rsqOO,            1424
.equiv          nb334_rsqH1H1,          1440
.equiv          nb334_rsqH1H2,          1456
.equiv          nb334_rsqH1M,           1472
.equiv          nb334_rsqH2H1,          1488
.equiv          nb334_rsqH2H2,          1504
.equiv          nb334_rsqH2M,           1520
.equiv          nb334_rsqMH1,           1536
.equiv          nb334_rsqMH2,           1552
.equiv          nb334_rsqMM,            1568
.equiv          nb334_rinvOO,           1584
.equiv          nb334_rinvH1H1,         1600
.equiv          nb334_rinvH1H2,         1616
.equiv          nb334_rinvH1M,          1632
.equiv          nb334_rinvH2H1,         1648
.equiv          nb334_rinvH2H2,         1664
.equiv          nb334_rinvH2M,          1680
.equiv          nb334_rinvMH1,          1696
.equiv          nb334_rinvMH2,          1712
.equiv          nb334_rinvMM,           1728
.equiv          nb334_fscal,            1744
.equiv          nb334_is3,              1760
.equiv          nb334_ii3,              1764
.equiv          nb334_innerjjnr,        1768
.equiv          nb334_innerk,           1772
.equiv          nb334_n,                1776
.equiv          nb334_nn1,              1780
.equiv          nb334_nri,              1784
.equiv          nb334_nouter,           1788
.equiv          nb334_ninner,           1792
.equiv          nb334_salign,           1796
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 1800		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb334_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb334_p_nri]
	mov ecx, [ecx]
	mov [esp + nb334_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb334_nouter], eax
	mov [esp + nb334_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb334_half], eax
	mov [esp + nb334_half+4], ebx
	movsd xmm1, [esp + nb334_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb334_half], xmm1
	movapd [esp + nb334_two], xmm2
	movapd [esp + nb334_three], xmm3
	mov eax, [ebp + nb334_p_tabscale]
	movsd xmm3, [eax]

	shufpd xmm3, xmm3, 0
	movapd [esp + nb334_tsc],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb334_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb334_charge]
	movsd xmm3, [edx + ebx*8 + 24]	
	movsd xmm4, xmm3	
	movsd xmm5, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb334_p_facel]
	movsd xmm6, [esi]
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm5
	mulsd  xmm5, xmm5
	mulsd  xmm3, xmm6
	mulsd  xmm4, xmm6
	mulsd  xmm5, xmm6
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb334_qqMM], xmm3
	movapd [esp + nb334_qqMH], xmm4
	movapd [esp + nb334_qqHH], xmm5
		
	xorpd xmm0, xmm0
	mov   edx, [ebp + nb334_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov edi, [ebp + nb334_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   eax, [ebp + nb334_vdwparam]
	movlpd xmm0, [eax + edx*8]
	movhpd xmm0, [eax + edx*8 + 8]
	movhlps xmm1, xmm0
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	movapd [esp + nb334_c6], xmm0
	movapd [esp + nb334_c12], xmm1

.nb334_threadloop:
        mov   esi, [ebp + nb334_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb334_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb334_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb334_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb334_n], eax
        mov [esp + nb334_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb334_outerstart
        jmp .nb334_end

.nb334_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb334_nouter]
	mov [esp + nb334_nouter], ebx

.nb334_outer:
	mov   eax, [ebp + nb334_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb334_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb334_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb334_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb334_pos]    ;# eax = base of pos[]  
	mov   [esp + nb334_ii3], ebx		

	addsd xmm3, [eax + ebx*8] 	;# ox
	addsd xmm4, [eax + ebx*8 + 8] 	;# oy
	addsd xmm5, [eax + ebx*8 + 16]	;# oz	
	addsd xmm6, [eax + ebx*8 + 24] 	;# h1x
	addsd xmm7, [eax + ebx*8 + 32] 	;# h1y
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	shufpd xmm7, xmm7, 0
	movapd [esp + nb334_ixO], xmm3
	movapd [esp + nb334_iyO], xmm4
	movapd [esp + nb334_izO], xmm5
	movapd [esp + nb334_ixH1], xmm6
	movapd [esp + nb334_iyH1], xmm7

	movsd xmm6, xmm2
	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm6, [eax + ebx*8 + 40] ;# h1z
	addsd xmm0, [eax + ebx*8 + 48] ;# h2x
	addsd xmm1, [eax + ebx*8 + 56] ;# h2y
	addsd xmm2, [eax + ebx*8 + 64] ;# h2z
	addsd xmm3, [eax + ebx*8 + 72] ;# mx
	addsd xmm4, [eax + ebx*8 + 80] ;# my
	addsd xmm5, [eax + ebx*8 + 88] ;# mz

	shufpd xmm6, xmm6, 0
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb334_izH1], xmm6
	movapd [esp + nb334_ixH2], xmm0
	movapd [esp + nb334_iyH2], xmm1
	movapd [esp + nb334_izH2], xmm2
	movapd [esp + nb334_ixM], xmm3
	movapd [esp + nb334_iyM], xmm4
	movapd [esp + nb334_izM], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb334_vctot], xmm4
	movapd [esp + nb334_Vvdwtot], xmm4
	movapd [esp + nb334_fixO], xmm4
	movapd [esp + nb334_fiyO], xmm4
	movapd [esp + nb334_fizO], xmm4
	movapd [esp + nb334_fixH1], xmm4
	movapd [esp + nb334_fiyH1], xmm4
	movapd [esp + nb334_fizH1], xmm4
	movapd [esp + nb334_fixH2], xmm4
	movapd [esp + nb334_fiyH2], xmm4
	movapd [esp + nb334_fizH2], xmm4
	movapd [esp + nb334_fixM], xmm4
	movapd [esp + nb334_fiyM], xmm4
	movapd [esp + nb334_fizM], xmm4
	
	mov   eax, [ebp + nb334_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb334_pos] 
	mov   edi, [ebp + nb334_faction]	
	mov   eax, [ebp + nb334_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb334_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb334_ninner]
	mov   [esp + nb334_ninner], ecx
	add   edx, 0
	mov   [esp + nb334_innerk], edx    ;# number of innerloop atoms 
	jge   .nb334_unroll_loop
	jmp   .nb334_checksingle
.nb334_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb334_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	
	add dword ptr [esp + nb334_innerjjnr], 8 ;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb334_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	
	;# move j coordinates to local temp variables 
	;# load ox, oy, oz, h1x
	movlpd xmm0, [esi + eax*8]
	movlpd xmm2, [esi + ebx*8]
	movhpd xmm0, [esi + eax*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 8]
	movlpd xmm3, [esi + eax*8 + 16]
	movlpd xmm5, [esi + ebx*8 + 16]
	movhpd xmm3, [esi + eax*8 + 24]
	movhpd xmm5, [esi + ebx*8 + 24]
	movapd xmm1, xmm0 
	movapd xmm4, xmm3
	unpcklpd xmm0, xmm2 ;# ox 
	unpckhpd xmm1, xmm2 ;# oy
	unpcklpd xmm3, xmm5 ;# ox 
	unpckhpd xmm4, xmm5 ;# oy
	movapd 	[esp + nb334_jxO], xmm0
	movapd 	[esp + nb334_jyO], xmm1
	movapd 	[esp + nb334_jzO], xmm3
	movapd 	[esp + nb334_jxH1], xmm4
	
	;# load h1y, h1z, h2x, h2y 
	movlpd xmm0, [esi + eax*8 + 32]
	movlpd xmm2, [esi + ebx*8 + 32]
	movhpd xmm0, [esi + eax*8 + 40]
	movhpd xmm2, [esi + ebx*8 + 40]
	movlpd xmm3, [esi + eax*8 + 48]
	movlpd xmm5, [esi + ebx*8 + 48]
	movhpd xmm3, [esi + eax*8 + 56]
	movhpd xmm5, [esi + ebx*8 + 56]
	movapd xmm1, xmm0 
	movapd xmm4, xmm3
	unpcklpd xmm0, xmm2 ;# h1y
	unpckhpd xmm1, xmm2 ;# h1z
	unpcklpd xmm3, xmm5 ;# h2x
	unpckhpd xmm4, xmm5 ;# h2y
	movapd 	[esp + nb334_jyH1], xmm0
	movapd 	[esp + nb334_jzH1], xmm1
	movapd 	[esp + nb334_jxH2], xmm3
	movapd 	[esp + nb334_jyH2], xmm4
	
	;# load h2z, mx, my, mz
	movlpd xmm0, [esi + eax*8 + 64]
	movlpd xmm2, [esi + ebx*8 + 64]
	movhpd xmm0, [esi + eax*8 + 72]
	movhpd xmm2, [esi + ebx*8 + 72]
	movlpd xmm3, [esi + eax*8 + 80]
	movlpd xmm5, [esi + ebx*8 + 80]
	movhpd xmm3, [esi + eax*8 + 88]
	movhpd xmm5, [esi + ebx*8 + 88]
	movapd xmm1, xmm0 
	movapd xmm4, xmm3
	unpcklpd xmm0, xmm2 ;# h2z
	unpckhpd xmm1, xmm2 ;# mx
	unpcklpd xmm3, xmm5 ;# my
	unpckhpd xmm4, xmm5 ;# mz
	movapd 	[esp + nb334_jzH2], xmm0
	movapd 	[esp + nb334_jxM], xmm1
	movapd 	[esp + nb334_jyM], xmm3
	movapd 	[esp + nb334_jzM], xmm4
	
	;# start calculating pairwise distances
	movapd xmm0, [esp + nb334_ixO]
	movapd xmm1, [esp + nb334_iyO]
	movapd xmm2, [esp + nb334_izO]
	movapd xmm3, [esp + nb334_ixH1]
	movapd xmm4, [esp + nb334_iyH1]
	movapd xmm5, [esp + nb334_izH1]
	subpd  xmm0, [esp + nb334_jxO]
	subpd  xmm1, [esp + nb334_jyO]
	subpd  xmm2, [esp + nb334_jzO]
	subpd  xmm3, [esp + nb334_jxH1]
	subpd  xmm4, [esp + nb334_jyH1]
	subpd  xmm5, [esp + nb334_jzH1]
	movapd [esp + nb334_dxOO], xmm0
	movapd [esp + nb334_dyOO], xmm1
	movapd [esp + nb334_dzOO], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb334_dxH1H1], xmm3
	movapd [esp + nb334_dyH1H1], xmm4
	movapd [esp + nb334_dzH1H1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb334_rsqOO], xmm0
	movapd [esp + nb334_rsqH1H1], xmm3

	movapd xmm0, [esp + nb334_ixH1]
	movapd xmm1, [esp + nb334_iyH1]
	movapd xmm2, [esp + nb334_izH1]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [esp + nb334_jxH2]
	subpd  xmm1, [esp + nb334_jyH2]
	subpd  xmm2, [esp + nb334_jzH2]
	subpd  xmm3, [esp + nb334_jxM]
	subpd  xmm4, [esp + nb334_jyM]
	subpd  xmm5, [esp + nb334_jzM]
	movapd [esp + nb334_dxH1H2], xmm0
	movapd [esp + nb334_dyH1H2], xmm1
	movapd [esp + nb334_dzH1H2], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb334_dxH1M], xmm3
	movapd [esp + nb334_dyH1M], xmm4
	movapd [esp + nb334_dzH1M], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb334_rsqH1H2], xmm0
	movapd [esp + nb334_rsqH1M], xmm3

	movapd xmm0, [esp + nb334_ixH2]
	movapd xmm1, [esp + nb334_iyH2]
	movapd xmm2, [esp + nb334_izH2]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [esp + nb334_jxH1]
	subpd  xmm1, [esp + nb334_jyH1]
	subpd  xmm2, [esp + nb334_jzH1]
	subpd  xmm3, [esp + nb334_jxH2]
	subpd  xmm4, [esp + nb334_jyH2]
	subpd  xmm5, [esp + nb334_jzH2]
	movapd [esp + nb334_dxH2H1], xmm0
	movapd [esp + nb334_dyH2H1], xmm1
	movapd [esp + nb334_dzH2H1], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb334_dxH2H2], xmm3
	movapd [esp + nb334_dyH2H2], xmm4
	movapd [esp + nb334_dzH2H2], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb334_rsqH2H1], xmm0
	movapd [esp + nb334_rsqH2H2], xmm3

	movapd xmm0, [esp + nb334_ixH2]
	movapd xmm1, [esp + nb334_iyH2]
	movapd xmm2, [esp + nb334_izH2]
	movapd xmm3, [esp + nb334_ixM]
	movapd xmm4, [esp + nb334_iyM]
	movapd xmm5, [esp + nb334_izM]
	subpd  xmm0, [esp + nb334_jxM]
	subpd  xmm1, [esp + nb334_jyM]
	subpd  xmm2, [esp + nb334_jzM]
	subpd  xmm3, [esp + nb334_jxH1]
	subpd  xmm4, [esp + nb334_jyH1]
	subpd  xmm5, [esp + nb334_jzH1]
	movapd [esp + nb334_dxH2M], xmm0
	movapd [esp + nb334_dyH2M], xmm1
	movapd [esp + nb334_dzH2M], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb334_dxMH1], xmm3
	movapd [esp + nb334_dyMH1], xmm4
	movapd [esp + nb334_dzMH1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm4, xmm3
	addpd  xmm4, xmm5
	movapd [esp + nb334_rsqH2M], xmm0
	movapd [esp + nb334_rsqMH1], xmm4

	movapd xmm0, [esp + nb334_ixM]
	movapd xmm1, [esp + nb334_iyM]
	movapd xmm2, [esp + nb334_izM]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [esp + nb334_jxH2]
	subpd  xmm1, [esp + nb334_jyH2]
	subpd  xmm2, [esp + nb334_jzH2]
	subpd  xmm3, [esp + nb334_jxM]
	subpd  xmm4, [esp + nb334_jyM]
	subpd  xmm5, [esp + nb334_jzM]
	movapd [esp + nb334_dxMH2], xmm0
	movapd [esp + nb334_dyMH2], xmm1
	movapd [esp + nb334_dzMH2], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb334_dxMM], xmm3
	movapd [esp + nb334_dyMM], xmm4
	movapd [esp + nb334_dzMM], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm4, xmm3
	addpd  xmm4, xmm5
	movapd [esp + nb334_rsqMH2], xmm0
	movapd [esp + nb334_rsqMM], xmm4
	
	;# Invsqrt form rsq M-H2 (rsq in xmm0) and MM (rsq in xmm4) 
	cvtpd2ps xmm1, xmm0
	cvtpd2ps xmm5, xmm4
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1   ;# luA
	cvtps2pd xmm5, xmm5   ;# luB
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb334_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb334_half] ;# iter1 
	mulpd   xmm7, [esp + nb334_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb334_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb334_half] ;# rinv 
	mulpd   xmm5, [esp + nb334_half] ;# rinv 
	movapd [esp + nb334_rinvMH2], xmm1
	movapd [esp + nb334_rinvMM], xmm5

	movapd xmm0, [esp + nb334_rsqOO]
	movapd xmm4, [esp + nb334_rsqH1H1]	
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb334_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb334_half] ;# iter1 of  
	mulpd   xmm7, [esp + nb334_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb334_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb334_half] ;# rinv 
	mulpd   xmm5, [esp + nb334_half] ;# rinv
	movapd [esp + nb334_rinvOO], xmm1
	movapd [esp + nb334_rinvH1H1], xmm5

	movapd xmm0, [esp + nb334_rsqH1H2]
	movapd xmm4, [esp + nb334_rsqH1M]	
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb334_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb334_half] ;# iter1 
	mulpd   xmm7, [esp + nb334_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb334_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb334_half] ;# rinv 
	mulpd   xmm5, [esp + nb334_half] ;# rinv 
	movapd [esp + nb334_rinvH1H2], xmm1
	movapd [esp + nb334_rinvH1M], xmm5

	movapd xmm0, [esp + nb334_rsqH2H1]
	movapd xmm4, [esp + nb334_rsqH2H2]	
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb334_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb334_half] ;# iter1a 
	mulpd   xmm7, [esp + nb334_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb334_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb334_half] ;# rinv 
	mulpd   xmm5, [esp + nb334_half] ;# rinv 
	movapd [esp + nb334_rinvH2H1], xmm1
	movapd [esp + nb334_rinvH2H2], xmm5

	movapd xmm0, [esp + nb334_rsqMH1]
	movapd xmm4, [esp + nb334_rsqH2M]	
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb334_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb334_half] ;# iter1a 
	mulpd   xmm7, [esp + nb334_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb334_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb334_half] ;# rinv 
	mulpd   xmm5, [esp + nb334_half] ;# rinv 
	movapd [esp + nb334_rinvMH1], xmm1
	movapd [esp + nb334_rinvH2M], xmm5

	;# start with OO interaction 
	movapd xmm0, [esp + nb334_rinvOO]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb334_rsqOO] ;# xmm1=r 
	mulpd  xmm1, [esp + nb334_tsc]

	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	movd mm0, eax	
	movd mm1, ebx
	mov  esi, [ebp + nb334_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 
	lea   ebx, [ebx + ebx*2]	;# idx*=3 (12 total now) 

	;# Dispersion 
	movlpd xmm4, [esi + eax*8 + 32]	;# Y1
	movlpd xmm3, [esi + ebx*8 + 32]	;# Y2
	movhpd xmm4, [esi + eax*8 + 40]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 40]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 48]	;# G1
	movlpd xmm3, [esi + ebx*8 + 48]	;# G2
	movhpd xmm6, [esi + eax*8 + 56]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 56]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# Dispersion table ready, in xmm4-xmm7  		

		mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb334_two]	;# two*Heps2 
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb334_c6]
	mulpd  xmm7, xmm4	 ;# fijD 
	mulpd  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack. Update Vvdwtot directly 
	addpd  xmm5, [esp + nb334_Vvdwtot]
	movapd [esp + nb334_fscal], xmm7
	movapd [esp + nb334_Vvdwtot], xmm5

	;# Repulsion 
	movlpd xmm4, [esi + eax*8 + 64]	;# Y1
	movlpd xmm3, [esi + ebx*8 + 64]	;# Y2
	movhpd xmm4, [esi + eax*8 + 72]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 72]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 80]	;# G1
	movlpd xmm3, [esi + ebx*8 + 80]	;# G2
	movhpd xmm6, [esi + eax*8 + 88]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 88]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# Repulsion table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb334_two]	;# two*Heps2 
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb334_c12]
	mulpd  xmm7, xmm4 ;# fijR 
	mulpd  xmm5, xmm4 ;# Vvdw12 
	addpd  xmm7, [esp + nb334_fscal]

	addpd  xmm5, [esp + nb334_Vvdwtot]
	movapd [esp + nb334_Vvdwtot], xmm5
	xorpd  xmm4, xmm4

	mulpd xmm7, [esp + nb334_tsc]
	mulpd xmm7, [esp + nb334_rinvOO]
	subpd xmm4, xmm7

	movapd xmm0, xmm4
	movapd xmm1, xmm4
	movapd xmm2, xmm4

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb334_dxOO]
	mulpd xmm1, [esp + nb334_dyOO]
	mulpd xmm2, [esp + nb334_dzOO]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb334_fixO]
	addpd xmm1, [esp + nb334_fiyO]
	addpd xmm2, [esp + nb334_fizO]
	movapd [esp + nb334_fjxO], xmm3
	movapd [esp + nb334_fjyO], xmm4
	movapd [esp + nb334_fjzO], xmm5
	movapd [esp + nb334_fixO], xmm0
	movapd [esp + nb334_fiyO], xmm1
	movapd [esp + nb334_fizO], xmm2

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb334_rinvH1H1]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb334_rsqH1H1] ;# xmm1=r 
	mulpd  xmm1, [esp + nb334_tsc]

	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb334_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 
	lea   ebx, [ebx + ebx*2]	;# idx*=3 (12 total now) 

	movlpd xmm4, [esi + eax*8]	;# Y1	
	movlpd xmm3, [esi + ebx*8]	;# Y2
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb334_two]	;# two*Heps2 
	movapd xmm3, [esp + nb334_qqHH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm3 fijC 

    	addpd  xmm5, [esp + nb334_vctot]
    	movapd [esp + nb334_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb334_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb334_dxH1H1]
	mulpd xmm1, [esp + nb334_dyH1H1]
	mulpd xmm2, [esp + nb334_dzH1H1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb334_fixH1]
	addpd xmm1, [esp + nb334_fiyH1]
	addpd xmm2, [esp + nb334_fizH1]
	movapd [esp + nb334_fjxH1], xmm3
	movapd [esp + nb334_fjyH1], xmm4
	movapd [esp + nb334_fjzH1], xmm5
	movapd [esp + nb334_fixH1], xmm0
	movapd [esp + nb334_fiyH1], xmm1
	movapd [esp + nb334_fizH1], xmm2

	;# H1-H2 interaction  
	movapd xmm0, [esp + nb334_rinvH1H2]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb334_rsqH1H2] ;# xmm1=r 
	mulpd  xmm1, [esp + nb334_tsc]
	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb334_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 
	lea   ebx, [ebx + ebx*2]	;# idx*=3 (12 total now) 

	movlpd xmm4, [esi + eax*8]	;# Y1	
	movlpd xmm3, [esi + ebx*8]	;# Y2
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb334_two]	;# two*Heps2 
	movapd xmm3, [esp + nb334_qqHH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm3 fijC 

    	addpd  xmm5, [esp + nb334_vctot]
    	movapd [esp + nb334_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb334_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb334_dxH1H2]
	mulpd xmm1, [esp + nb334_dyH1H2]
	mulpd xmm2, [esp + nb334_dzH1H2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb334_fixH1]
	addpd xmm1, [esp + nb334_fiyH1]
	addpd xmm2, [esp + nb334_fizH1]
	movapd [esp + nb334_fjxH2], xmm3
	movapd [esp + nb334_fjyH2], xmm4
	movapd [esp + nb334_fjzH2], xmm5
	movapd [esp + nb334_fixH1], xmm0
	movapd [esp + nb334_fiyH1], xmm1
	movapd [esp + nb334_fizH1], xmm2

	;# H1-M interaction 
	movapd xmm0, [esp + nb334_rinvH1M]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb334_rsqH1M] ;# xmm1=r 
	mulpd  xmm1, [esp + nb334_tsc]
	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb334_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 
	lea   ebx, [ebx + ebx*2]	;# idx*=3 (12 total now) 

	movlpd xmm4, [esi + eax*8]	;# Y1	
	movlpd xmm3, [esi + ebx*8]	;# Y2
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb334_two]	;# two*Heps2 
	movapd xmm3, [esp + nb334_qqMH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm3 fijC 

    	addpd  xmm5, [esp + nb334_vctot]
    	movapd [esp + nb334_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb334_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb334_dxH1M]
	mulpd xmm1, [esp + nb334_dyH1M]
	mulpd xmm2, [esp + nb334_dzH1M]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb334_fixH1]
	addpd xmm1, [esp + nb334_fiyH1]
	addpd xmm2, [esp + nb334_fizH1]
	movapd [esp + nb334_fjxM], xmm3
	movapd [esp + nb334_fjyM], xmm4
	movapd [esp + nb334_fjzM], xmm5
	movapd [esp + nb334_fixH1], xmm0
	movapd [esp + nb334_fiyH1], xmm1
	movapd [esp + nb334_fizH1], xmm2

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb334_rinvH2H1]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb334_rsqH2H1] ;# xmm1=r 
	mulpd  xmm1, [esp + nb334_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb334_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 
	lea   ebx, [ebx + ebx*2]	;# idx*=3 (12 total now) 

	movlpd xmm4, [esi + eax*8]	;# Y1	
	movlpd xmm3, [esi + ebx*8]	;# Y2
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb334_two]	;# two*Heps2 
	movapd xmm3, [esp + nb334_qqHH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm3 fijC 

    	addpd  xmm5, [esp + nb334_vctot]
    	movapd [esp + nb334_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb334_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb334_fjxH1]
	movapd xmm4, [esp + nb334_fjyH1]
	movapd xmm5, [esp + nb334_fjzH1]
	mulpd xmm0, [esp + nb334_dxH2H1]
	mulpd xmm1, [esp + nb334_dyH2H1]
	mulpd xmm2, [esp + nb334_dzH2H1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb334_fixH2]
	addpd xmm1, [esp + nb334_fiyH2]
	addpd xmm2, [esp + nb334_fizH2]
	movapd [esp + nb334_fjxH1], xmm3
	movapd [esp + nb334_fjyH1], xmm4
	movapd [esp + nb334_fjzH1], xmm5
	movapd [esp + nb334_fixH2], xmm0
	movapd [esp + nb334_fiyH2], xmm1
	movapd [esp + nb334_fizH2], xmm2

	;# H2-H2 interaction 
	movapd xmm0, [esp + nb334_rinvH2H2]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb334_rsqH2H2] ;# xmm1=r 
	mulpd  xmm1, [esp + nb334_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb334_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 
	lea   ebx, [ebx + ebx*2]	;# idx*=3 (12 total now) 

	movlpd xmm4, [esi + eax*8]	;# Y1	
	movlpd xmm3, [esi + ebx*8]	;# Y2
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb334_two]	;# two*Heps2 
	movapd xmm3, [esp + nb334_qqHH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm3 fijC 

    	addpd  xmm5, [esp + nb334_vctot]
    	movapd [esp + nb334_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb334_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb334_fjxH2]
	movapd xmm4, [esp + nb334_fjyH2]
	movapd xmm5, [esp + nb334_fjzH2]
	mulpd xmm0, [esp + nb334_dxH2H2]
	mulpd xmm1, [esp + nb334_dyH2H2]
	mulpd xmm2, [esp + nb334_dzH2H2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb334_fixH2]
	addpd xmm1, [esp + nb334_fiyH2]
	addpd xmm2, [esp + nb334_fizH2]
	movapd [esp + nb334_fjxH2], xmm3
	movapd [esp + nb334_fjyH2], xmm4
	movapd [esp + nb334_fjzH2], xmm5
	movapd [esp + nb334_fixH2], xmm0
	movapd [esp + nb334_fiyH2], xmm1
	movapd [esp + nb334_fizH2], xmm2

	;# H2-M interaction 
	movapd xmm0, [esp + nb334_rinvH2M]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb334_rsqH2M] ;# xmm1=r 
	mulpd  xmm1, [esp + nb334_tsc]
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb334_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 
	lea   ebx, [ebx + ebx*2]	;# idx*=3 (12 total now) 

	movlpd xmm4, [esi + eax*8]	;# Y1	
	movlpd xmm3, [esi + ebx*8]	;# Y2
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb334_two]	;# two*Heps2 
	movapd xmm3, [esp + nb334_qqMH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm3 fijC 

    	addpd  xmm5, [esp + nb334_vctot]
    	movapd [esp + nb334_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb334_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb334_fjxM]
	movapd xmm4, [esp + nb334_fjyM]
	movapd xmm5, [esp + nb334_fjzM]
	mulpd xmm0, [esp + nb334_dxH2M]
	mulpd xmm1, [esp + nb334_dyH2M]
	mulpd xmm2, [esp + nb334_dzH2M]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb334_fixH2]
	addpd xmm1, [esp + nb334_fiyH2]
	addpd xmm2, [esp + nb334_fizH2]
	movapd [esp + nb334_fjxM], xmm3
	movapd [esp + nb334_fjyM], xmm4
	movapd [esp + nb334_fjzM], xmm5
	movapd [esp + nb334_fixH2], xmm0
	movapd [esp + nb334_fiyH2], xmm1
	movapd [esp + nb334_fizH2], xmm2

	;# M-H1 interaction 
	movapd xmm0, [esp + nb334_rinvMH1]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb334_rsqMH1] ;# xmm1=r 
	mulpd  xmm1, [esp + nb334_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb334_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 
	lea   ebx, [ebx + ebx*2]	;# idx*=3 (12 total now) 

	movlpd xmm4, [esi + eax*8]	;# Y1	
	movlpd xmm3, [esi + ebx*8]	;# Y2
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb334_two]	;# two*Heps2 
	movapd xmm3, [esp + nb334_qqMH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm3 fijC 

    	addpd  xmm5, [esp + nb334_vctot]
    	movapd [esp + nb334_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb334_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1

	movapd xmm3, [esp + nb334_fjxH1]
	movapd xmm4, [esp + nb334_fjyH1]
	movapd xmm5, [esp + nb334_fjzH1]
	mulpd xmm0, [esp + nb334_dxMH1]
	mulpd xmm1, [esp + nb334_dyMH1]
	mulpd xmm2, [esp + nb334_dzMH1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb334_fixM]
	addpd xmm1, [esp + nb334_fiyM]
	addpd xmm2, [esp + nb334_fizM]
	movapd [esp + nb334_fjxH1], xmm3
	movapd [esp + nb334_fjyH1], xmm4
	movapd [esp + nb334_fjzH1], xmm5
	movapd [esp + nb334_fixM], xmm0
	movapd [esp + nb334_fiyM], xmm1
	movapd [esp + nb334_fizM], xmm2

	;# M-H2 interaction 
	movapd xmm0, [esp + nb334_rinvMH2]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb334_rsqMH2] ;# xmm1=r 
	mulpd  xmm1, [esp + nb334_tsc]
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb334_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 
	lea   ebx, [ebx + ebx*2]	;# idx*=3 (12 total now) 

	movlpd xmm4, [esi + eax*8]	;# Y1	
	movlpd xmm3, [esi + ebx*8]	;# Y2
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb334_two]	;# two*Heps2 
	movapd xmm3, [esp + nb334_qqMH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm3 fijC 

    	addpd  xmm5, [esp + nb334_vctot]
    	movapd [esp + nb334_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb334_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb334_fjxH2]
	movapd xmm4, [esp + nb334_fjyH2]
	movapd xmm5, [esp + nb334_fjzH2]
	mulpd xmm0, [esp + nb334_dxMH2]
	mulpd xmm1, [esp + nb334_dyMH2]
	mulpd xmm2, [esp + nb334_dzMH2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb334_fixM]
	addpd xmm1, [esp + nb334_fiyM]
	addpd xmm2, [esp + nb334_fizM]
	movapd [esp + nb334_fjxH2], xmm3
	movapd [esp + nb334_fjyH2], xmm4
	movapd [esp + nb334_fjzH2], xmm5
	movapd [esp + nb334_fixM], xmm0
	movapd [esp + nb334_fiyM], xmm1
	movapd [esp + nb334_fizM], xmm2

	;# M-M interaction 
	movapd xmm0, [esp + nb334_rinvMM]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb334_rsqMM] ;# xmm1=r 
	mulpd  xmm1, [esp + nb334_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb334_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 
	lea   ebx, [ebx + ebx*2]	;# idx*=3 (12 total now) 

	movlpd xmm4, [esi + eax*8]	;# Y1	
	movlpd xmm3, [esi + ebx*8]	;# Y2
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb334_two]	;# two*Heps2 
	movapd xmm3, [esp + nb334_qqMM]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm3 fijC 
	
    	addpd  xmm5, [esp + nb334_vctot]
    	movapd [esp + nb334_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb334_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb334_fjxM]
	movapd xmm4, [esp + nb334_fjyM] 
	movapd xmm5, [esp + nb334_fjzM]
	mulpd xmm0, [esp + nb334_dxMM]
	mulpd xmm1, [esp + nb334_dyMM]
	mulpd xmm2, [esp + nb334_dzMM]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb334_fixM]
	addpd xmm1, [esp + nb334_fiyM]
	addpd xmm2, [esp + nb334_fizM]
	movapd [esp + nb334_fjxM], xmm3
	movapd [esp + nb334_fjyM], xmm4
	movapd [esp + nb334_fjzM], xmm5
	movapd [esp + nb334_fixM], xmm0
	movapd [esp + nb334_fiyM], xmm1
	movapd [esp + nb334_fizM], xmm2

	mov edi, [ebp + nb334_faction]

	movd eax, mm0
	movd ebx, mm1
	
	;# Did all interactions - now update j forces 
	;# Step1 - transpose fjxO, fjyO and fjzO, fjxH1
	movapd xmm0, [esp + nb334_fjxO]
	movapd xmm1, [esp + nb334_fjyO]
	movapd xmm2, [esp + nb334_fjzO]
	movapd xmm3, [esp + nb334_fjxH1]
	movapd xmm4, xmm0
	movapd xmm5, xmm2
	unpcklpd xmm0, xmm1   ;# fjOxA fjOyA
	unpckhpd xmm4, xmm1   ;# fjOxB fjOyB
	unpcklpd xmm2, xmm3   ;# fjOzA fjH1xA
	unpckhpd xmm5, xmm3   ;# fjOzB fjH1xB
	movlpd xmm1, [edi + eax*8]
	movlpd xmm3, [edi + ebx*8]
	movhpd xmm1, [edi + eax*8 + 8]
	movhpd xmm3, [edi + ebx*8 + 8]
	movlpd xmm6, [edi + eax*8 + 16]
	movlpd xmm7, [edi + ebx*8 + 16]
	movhpd xmm6, [edi + eax*8 + 24]
	movhpd xmm7, [edi + ebx*8 + 24]
	addpd  xmm1, xmm0
	addpd  xmm3, xmm4
	addpd  xmm6, xmm2
	addpd  xmm7, xmm5
	movlpd [edi + eax*8],      xmm1
	movlpd [edi + ebx*8],      xmm3
	movhpd [edi + eax*8 + 8],  xmm1
	movhpd [edi + ebx*8 + 8],  xmm3
	movlpd [edi + eax*8 + 16], xmm6
	movlpd [edi + ebx*8 + 16], xmm7
	movhpd [edi + eax*8 + 24], xmm6
	movhpd [edi + ebx*8 + 24], xmm7

	;# Step2 - transpose fjyH1, fjzH1 and fjxH2, fjyH2
	movapd xmm0, [esp + nb334_fjyH1]
	movapd xmm1, [esp + nb334_fjzH1]
	movapd xmm2, [esp + nb334_fjxH2]
	movapd xmm3, [esp + nb334_fjyH2]
	movapd xmm4, xmm0
	movapd xmm5, xmm2
	unpcklpd xmm0, xmm1   ;# fjOxA fjOyA
	unpckhpd xmm4, xmm1   ;# fjOxB fjOyB
	unpcklpd xmm2, xmm3   ;# fjOzA fjH1xA
	unpckhpd xmm5, xmm3   ;# fjOzB fjH1xB
	movlpd xmm1, [edi + eax*8 + 32]
	movlpd xmm3, [edi + ebx*8 + 32]
	movhpd xmm1, [edi + eax*8 + 40]
	movhpd xmm3, [edi + ebx*8 + 40]
	movlpd xmm6, [edi + eax*8 + 48]
	movlpd xmm7, [edi + ebx*8 + 48]
	movhpd xmm6, [edi + eax*8 + 56]
	movhpd xmm7, [edi + ebx*8 + 56]
	addpd  xmm1, xmm0
	addpd  xmm3, xmm4
	addpd  xmm6, xmm2
	addpd  xmm7, xmm5
	movlpd [edi + eax*8 + 32], xmm1
	movlpd [edi + ebx*8 + 32], xmm3
	movhpd [edi + eax*8 + 40], xmm1
	movhpd [edi + ebx*8 + 40], xmm3
	movlpd [edi + eax*8 + 48], xmm6
	movlpd [edi + ebx*8 + 48], xmm7
	movhpd [edi + eax*8 + 56], xmm6
	movhpd [edi + ebx*8 + 56], xmm7

	;# Step3 - transpose fjzH2, fjxM and fjyM, fjzM
	movapd xmm0, [esp + nb334_fjzH2]
	movapd xmm1, [esp + nb334_fjxM]
	movapd xmm2, [esp + nb334_fjyM]
	movapd xmm3, [esp + nb334_fjzM]
	movapd xmm4, xmm0
	movapd xmm5, xmm2
	unpcklpd xmm0, xmm1   ;# fjOxA fjOyA
	unpckhpd xmm4, xmm1   ;# fjOxB fjOyB
	unpcklpd xmm2, xmm3   ;# fjOzA fjH1xA
	unpckhpd xmm5, xmm3   ;# fjOzB fjH1xB
	movlpd xmm1, [edi + eax*8 + 64]
	movlpd xmm3, [edi + ebx*8 + 64]
	movhpd xmm1, [edi + eax*8 + 72]
	movhpd xmm3, [edi + ebx*8 + 72]
	movlpd xmm6, [edi + eax*8 + 80]
	movlpd xmm7, [edi + ebx*8 + 80]
	movhpd xmm6, [edi + eax*8 + 88]
	movhpd xmm7, [edi + ebx*8 + 88]
	addpd  xmm1, xmm0
	addpd  xmm3, xmm4
	addpd  xmm6, xmm2
	addpd  xmm7, xmm5
	movlpd [edi + eax*8 + 64], xmm1
	movlpd [edi + ebx*8 + 64], xmm3
	movhpd [edi + eax*8 + 72], xmm1
	movhpd [edi + ebx*8 + 72], xmm3
	movlpd [edi + eax*8 + 80], xmm6
	movlpd [edi + ebx*8 + 80], xmm7
	movhpd [edi + eax*8 + 88], xmm6
	movhpd [edi + ebx*8 + 88], xmm7

	;# should we do one more iteration? 
	sub dword ptr [esp + nb334_innerk],  2
	jl    .nb334_checksingle
	jmp   .nb334_unroll_loop
.nb334_checksingle:
	mov   edx, [esp + nb334_innerk]
	and   edx, 1
	jnz   .nb334_dosingle
	jmp   .nb334_updateouterdata
.nb334_dosingle:
	mov   edx, [esp + nb334_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]

	mov esi, [ebp + nb334_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	
	;# move j coordinates to local temp variables 
	;# load ox, oy, oz, h1x
	movlpd xmm0, [esi + eax*8] 
	movhpd xmm0, [esi + eax*8 + 8]
	movlpd xmm1, [esi + eax*8 + 16]
	movhpd xmm1, [esi + eax*8 + 24]
	movlpd xmm2, [esi + eax*8 + 32]
	movhpd xmm2, [esi + eax*8 + 40]
	movlpd xmm3, [esi + eax*8 + 48]
	movhpd xmm3, [esi + eax*8 + 56]
	movlpd xmm4, [esi + eax*8 + 64]
	movhpd xmm4, [esi + eax*8 + 72]
	movlpd xmm5, [esi + eax*8 + 80]
	movhpd xmm5, [esi + eax*8 + 88]
	movsd  [esp + nb334_jxO], xmm0
	movsd  [esp + nb334_jzO], xmm1
	movsd  [esp + nb334_jyH1], xmm2
	movsd  [esp + nb334_jxH2], xmm3
	movsd  [esp + nb334_jzH2], xmm4
	movsd  [esp + nb334_jyM], xmm5
	unpckhpd xmm0, xmm0
	unpckhpd xmm1, xmm1
	unpckhpd xmm2, xmm2
	unpckhpd xmm3, xmm3
	unpckhpd xmm4, xmm4
	unpckhpd xmm5, xmm5
	movsd  [esp + nb334_jyO], xmm0
	movsd  [esp + nb334_jxH1], xmm1
	movsd  [esp + nb334_jzH1], xmm2
	movsd  [esp + nb334_jyH2], xmm3
	movsd  [esp + nb334_jxM], xmm4
	movsd  [esp + nb334_jzM], xmm5

	;# start calculating pairwise distances
	movapd xmm0, [esp + nb334_ixO]
	movapd xmm1, [esp + nb334_iyO]
	movapd xmm2, [esp + nb334_izO]
	movapd xmm3, [esp + nb334_ixH1]
	movapd xmm4, [esp + nb334_iyH1]
	movapd xmm5, [esp + nb334_izH1]
	subsd  xmm0, [esp + nb334_jxO]
	subsd  xmm1, [esp + nb334_jyO]
	subsd  xmm2, [esp + nb334_jzO]
	subsd  xmm3, [esp + nb334_jxH1]
	subsd  xmm4, [esp + nb334_jyH1]
	subsd  xmm5, [esp + nb334_jzH1]
	movapd [esp + nb334_dxOO], xmm0
	movapd [esp + nb334_dyOO], xmm1
	movapd [esp + nb334_dzOO], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb334_dxH1H1], xmm3
	movapd [esp + nb334_dyH1H1], xmm4
	movapd [esp + nb334_dzH1H1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb334_rsqOO], xmm0
	movapd [esp + nb334_rsqH1H1], xmm3

	movapd xmm0, [esp + nb334_ixH1]
	movapd xmm1, [esp + nb334_iyH1]
	movapd xmm2, [esp + nb334_izH1]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [esp + nb334_jxH2]
	subsd  xmm1, [esp + nb334_jyH2]
	subsd  xmm2, [esp + nb334_jzH2]
	subsd  xmm3, [esp + nb334_jxM]
	subsd  xmm4, [esp + nb334_jyM]
	subsd  xmm5, [esp + nb334_jzM]
	movapd [esp + nb334_dxH1H2], xmm0
	movapd [esp + nb334_dyH1H2], xmm1
	movapd [esp + nb334_dzH1H2], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb334_dxH1M], xmm3
	movapd [esp + nb334_dyH1M], xmm4
	movapd [esp + nb334_dzH1M], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb334_rsqH1H2], xmm0
	movapd [esp + nb334_rsqH1M], xmm3

	movapd xmm0, [esp + nb334_ixH2]
	movapd xmm1, [esp + nb334_iyH2]
	movapd xmm2, [esp + nb334_izH2]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [esp + nb334_jxH1]
	subsd  xmm1, [esp + nb334_jyH1]
	subsd  xmm2, [esp + nb334_jzH1]
	subsd  xmm3, [esp + nb334_jxH2]
	subsd  xmm4, [esp + nb334_jyH2]
	subsd  xmm5, [esp + nb334_jzH2]
	movapd [esp + nb334_dxH2H1], xmm0
	movapd [esp + nb334_dyH2H1], xmm1
	movapd [esp + nb334_dzH2H1], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb334_dxH2H2], xmm3
	movapd [esp + nb334_dyH2H2], xmm4
	movapd [esp + nb334_dzH2H2], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb334_rsqH2H1], xmm0
	movapd [esp + nb334_rsqH2H2], xmm3

	movapd xmm0, [esp + nb334_ixH2]
	movapd xmm1, [esp + nb334_iyH2]
	movapd xmm2, [esp + nb334_izH2]
	movapd xmm3, [esp + nb334_ixM]
	movapd xmm4, [esp + nb334_iyM]
	movapd xmm5, [esp + nb334_izM]
	subsd  xmm0, [esp + nb334_jxM]
	subsd  xmm1, [esp + nb334_jyM]
	subsd  xmm2, [esp + nb334_jzM]
	subsd  xmm3, [esp + nb334_jxH1]
	subsd  xmm4, [esp + nb334_jyH1]
	subsd  xmm5, [esp + nb334_jzH1]
	movapd [esp + nb334_dxH2M], xmm0
	movapd [esp + nb334_dyH2M], xmm1
	movapd [esp + nb334_dzH2M], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb334_dxMH1], xmm3
	movapd [esp + nb334_dyMH1], xmm4
	movapd [esp + nb334_dzMH1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm4, xmm3
	addsd  xmm4, xmm5
	movapd [esp + nb334_rsqH2M], xmm0
	movapd [esp + nb334_rsqMH1], xmm4

	movapd xmm0, [esp + nb334_ixM]
	movapd xmm1, [esp + nb334_iyM]
	movapd xmm2, [esp + nb334_izM]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [esp + nb334_jxH2]
	subsd  xmm1, [esp + nb334_jyH2]
	subsd  xmm2, [esp + nb334_jzH2]
	subsd  xmm3, [esp + nb334_jxM]
	subsd  xmm4, [esp + nb334_jyM]
	subsd  xmm5, [esp + nb334_jzM]
	movapd [esp + nb334_dxMH2], xmm0
	movapd [esp + nb334_dyMH2], xmm1
	movapd [esp + nb334_dzMH2], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb334_dxMM], xmm3
	movapd [esp + nb334_dyMM], xmm4
	movapd [esp + nb334_dzMM], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm4, xmm3
	addsd  xmm4, xmm5
	movapd [esp + nb334_rsqMH2], xmm0
	movapd [esp + nb334_rsqMM], xmm4

	;# Invsqrt form rsq M-H2 (rsq in xmm0) and MM (rsq in xmm4) 
	cvtsd2ss xmm1, xmm0
	cvtsd2ss xmm5, xmm4
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1   ;# luA
	cvtss2sd xmm5, xmm5   ;# luB
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb334_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb334_half] ;# iter1 
	mulsd   xmm7, [esp + nb334_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb334_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb334_half] ;# rinv 
	mulsd   xmm5, [esp + nb334_half] ;# rinv 
	movapd [esp + nb334_rinvMH2], xmm1
	movapd [esp + nb334_rinvMM], xmm5

	movapd xmm0, [esp + nb334_rsqOO]
	movapd xmm4, [esp + nb334_rsqH1H1]	
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb334_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb334_half] ;# iter1 of  
	mulsd   xmm7, [esp + nb334_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb334_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb334_half] ;# rinv 
	mulsd   xmm5, [esp + nb334_half] ;# rinv
	movapd [esp + nb334_rinvOO], xmm1
	movapd [esp + nb334_rinvH1H1], xmm5

	movapd xmm0, [esp + nb334_rsqH1H2]
	movapd xmm4, [esp + nb334_rsqH1M]	
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb334_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb334_half] ;# iter1 
	mulsd   xmm7, [esp + nb334_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb334_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb334_half] ;# rinv 
	mulsd   xmm5, [esp + nb334_half] ;# rinv 
	movapd [esp + nb334_rinvH1H2], xmm1
	movapd [esp + nb334_rinvH1M], xmm5

	movapd xmm0, [esp + nb334_rsqH2H1]
	movapd xmm4, [esp + nb334_rsqH2H2]	
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb334_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb334_half] ;# iter1a 
	mulsd   xmm7, [esp + nb334_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb334_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb334_half] ;# rinv 
	mulsd   xmm5, [esp + nb334_half] ;# rinv 
	movapd [esp + nb334_rinvH2H1], xmm1
	movapd [esp + nb334_rinvH2H2], xmm5

	movapd xmm0, [esp + nb334_rsqMH1]
	movapd xmm4, [esp + nb334_rsqH2M]	
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb334_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb334_half] ;# iter1a 
	mulsd   xmm7, [esp + nb334_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb334_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb334_half] ;# rinv 
	mulsd   xmm5, [esp + nb334_half] ;# rinv 
	movapd [esp + nb334_rinvMH1], xmm1
	movapd [esp + nb334_rinvH2M], xmm5

	;# start with OO interaction 
	movapd xmm0, [esp + nb334_rinvOO]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb334_rsqOO] ;# xmm1=r 
	mulsd  xmm1, [esp + nb334_tsc]

	movd mm0, eax	


	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 

		
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb334_VFtab]
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 

	;# Dispersion 
	movsd xmm4, [esi + eax*8 + 32]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 40]	;# F1 	
	movsd xmm6, [esi + eax*8 + 48]	;# G1	
	movsd xmm7, [esi + eax*8 + 56]	;# H1 	

	;# Dispersion table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb334_two]	;# two*Heps2 
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb334_c6]
	mulsd  xmm7, xmm4	 ;# fijD 
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack. Update Vvdwtot directly 
	addsd  xmm5, [esp + nb334_Vvdwtot]
	movsd [esp + nb334_fscal], xmm7
	movsd [esp + nb334_Vvdwtot], xmm5

	;# Repulsion 
	movsd xmm4, [esi + eax*8 + 64]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 72]	;# F1 	
	movsd xmm6, [esi + eax*8 + 80]	;# G1
	movsd xmm7, [esi + eax*8 + 88]	;# H1 	
	;# Repulsion table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb334_two]	;# two*Heps2 
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb334_c12]
	mulsd  xmm7, xmm4 ;# fijR 
	mulsd  xmm5, xmm4 ;# Vvdw12 
	addsd  xmm7, [esp + nb334_fscal]

	addsd  xmm5, [esp + nb334_Vvdwtot]
	movsd [esp + nb334_Vvdwtot], xmm5
	xorpd  xmm4, xmm4

	mulsd xmm7, [esp + nb334_tsc]
	mulsd xmm7, [esp + nb334_rinvOO]
	subsd xmm4, xmm7

	movapd xmm0, xmm4
	movapd xmm1, xmm4
	movapd xmm2, xmm4

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb334_dxOO]
	mulsd xmm1, [esp + nb334_dyOO]
	mulsd xmm2, [esp + nb334_dzOO]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb334_fixO]
	addsd xmm1, [esp + nb334_fiyO]
	addsd xmm2, [esp + nb334_fizO]
	movsd [esp + nb334_fjxO], xmm3
	movsd [esp + nb334_fjyO], xmm4
	movsd [esp + nb334_fjzO], xmm5
	movsd [esp + nb334_fixO], xmm0
	movsd [esp + nb334_fiyO], xmm1
	movsd [esp + nb334_fizO], xmm2

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb334_rinvH1H1]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb334_rsqH1H1] ;# xmm1=r 
	mulsd  xmm1, [esp + nb334_tsc]

	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb334_VFtab]
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	

	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb334_two]	;# two*Heps2 
	movapd xmm3, [esp + nb334_qqHH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm3 fijC 

    	addsd  xmm5, [esp + nb334_vctot]
    	movsd [esp + nb334_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb334_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb334_dxH1H1]
	mulsd xmm1, [esp + nb334_dyH1H1]
	mulsd xmm2, [esp + nb334_dzH1H1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb334_fixH1]
	addsd xmm1, [esp + nb334_fiyH1]
	addsd xmm2, [esp + nb334_fizH1]
	movsd [esp + nb334_fjxH1], xmm3
	movsd [esp + nb334_fjyH1], xmm4
	movsd [esp + nb334_fjzH1], xmm5
	movsd [esp + nb334_fixH1], xmm0
	movsd [esp + nb334_fiyH1], xmm1
	movsd [esp + nb334_fizH1], xmm2

	;# H1-H2 interaction  
	movapd xmm0, [esp + nb334_rinvH1H2]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb334_rsqH1H2] ;# xmm1=r 
	mulsd  xmm1, [esp + nb334_tsc]
	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb334_VFtab]
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 

	movsd xmm4, [esi + eax*8]	;# Y1  	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb334_two]	;# two*Heps2 
	movapd xmm3, [esp + nb334_qqHH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm3 fijC 

    	addsd  xmm5, [esp + nb334_vctot]
    	movsd [esp + nb334_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb334_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb334_dxH1H2]
	mulsd xmm1, [esp + nb334_dyH1H2]
	mulsd xmm2, [esp + nb334_dzH1H2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb334_fixH1]
	addsd xmm1, [esp + nb334_fiyH1]
	addsd xmm2, [esp + nb334_fizH1]
	movsd [esp + nb334_fjxH2], xmm3
	movsd [esp + nb334_fjyH2], xmm4
	movsd [esp + nb334_fjzH2], xmm5
	movsd [esp + nb334_fixH1], xmm0
	movsd [esp + nb334_fiyH1], xmm1
	movsd [esp + nb334_fizH1], xmm2

	;# H1-M interaction 
	movapd xmm0, [esp + nb334_rinvH1M]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb334_rsqH1M] ;# xmm1=r 
	mulsd  xmm1, [esp + nb334_tsc]
	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb334_VFtab]
	lea  eax, [eax + eax*2]	;# idx*=3 (12 total now) 

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb334_two]	;# two*Heps2 
	movapd xmm3, [esp + nb334_qqMH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm3 fijC 

    	addsd  xmm5, [esp + nb334_vctot]
    	movsd [esp + nb334_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb334_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb334_dxH1M]
	mulsd xmm1, [esp + nb334_dyH1M]
	mulsd xmm2, [esp + nb334_dzH1M]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb334_fixH1]
	addsd xmm1, [esp + nb334_fiyH1]
	addsd xmm2, [esp + nb334_fizH1]
	movsd [esp + nb334_fjxM], xmm3
	movsd [esp + nb334_fjyM], xmm4
	movsd [esp + nb334_fjzM], xmm5
	movsd [esp + nb334_fixH1], xmm0
	movsd [esp + nb334_fiyH1], xmm1
	movsd [esp + nb334_fizH1], xmm2

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb334_rinvH2H1]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb334_rsqH2H1] ;# xmm1=r 
	mulsd  xmm1, [esp + nb334_tsc]	

	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb334_VFtab]
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb334_two]	;# two*Heps2 
	movapd xmm3, [esp + nb334_qqHH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm3 fijC 

    	addsd  xmm5, [esp + nb334_vctot]
    	movsd [esp + nb334_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb334_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb334_fjxH1]
	movapd xmm4, [esp + nb334_fjyH1]
	movapd xmm5, [esp + nb334_fjzH1]
	mulsd xmm0, [esp + nb334_dxH2H1]
	mulsd xmm1, [esp + nb334_dyH2H1]
	mulsd xmm2, [esp + nb334_dzH2H1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb334_fixH2]
	addsd xmm1, [esp + nb334_fiyH2]
	addsd xmm2, [esp + nb334_fizH2]
	movsd [esp + nb334_fjxH1], xmm3
	movsd [esp + nb334_fjyH1], xmm4
	movsd [esp + nb334_fjzH1], xmm5
	movsd [esp + nb334_fixH2], xmm0
	movsd [esp + nb334_fiyH2], xmm1
	movsd [esp + nb334_fizH2], xmm2

	;# H2-H2 interaction 
	movapd xmm0, [esp + nb334_rinvH2H2]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb334_rsqH2H2] ;# xmm1=r 
	mulsd  xmm1, [esp + nb334_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb334_VFtab]
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb334_two]	;# two*Heps2 
	movapd xmm3, [esp + nb334_qqHH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm3 fijC 

    	addsd  xmm5, [esp + nb334_vctot]
    	movsd [esp + nb334_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb334_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb334_fjxH2]
	movapd xmm4, [esp + nb334_fjyH2]
	movapd xmm5, [esp + nb334_fjzH2]
	mulsd xmm0, [esp + nb334_dxH2H2]
	mulsd xmm1, [esp + nb334_dyH2H2]
	mulsd xmm2, [esp + nb334_dzH2H2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb334_fixH2]
	addsd xmm1, [esp + nb334_fiyH2]
	addsd xmm2, [esp + nb334_fizH2]
	movsd [esp + nb334_fjxH2], xmm3
	movsd [esp + nb334_fjyH2], xmm4
	movsd [esp + nb334_fjzH2], xmm5
	movsd [esp + nb334_fixH2], xmm0
	movsd [esp + nb334_fiyH2], xmm1
	movsd [esp + nb334_fizH2], xmm2

	;# H2-M interaction 
	movapd xmm0, [esp + nb334_rinvH2M]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb334_rsqH2M] ;# xmm1=r 
	mulsd  xmm1, [esp + nb334_tsc]
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb334_VFtab]
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb334_two]	;# two*Heps2 
	movapd xmm3, [esp + nb334_qqMH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm3 fijC 

    	addsd  xmm5, [esp + nb334_vctot]
    	movsd [esp + nb334_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb334_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb334_fjxM]
	movapd xmm4, [esp + nb334_fjyM]
	movapd xmm5, [esp + nb334_fjzM]
	mulsd xmm0, [esp + nb334_dxH2M]
	mulsd xmm1, [esp + nb334_dyH2M]
	mulsd xmm2, [esp + nb334_dzH2M]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb334_fixH2]
	addsd xmm1, [esp + nb334_fiyH2]
	addsd xmm2, [esp + nb334_fizH2]
	movsd [esp + nb334_fjxM], xmm3
	movsd [esp + nb334_fjyM], xmm4
	movsd [esp + nb334_fjzM], xmm5
	movsd [esp + nb334_fixH2], xmm0
	movsd [esp + nb334_fiyH2], xmm1
	movsd [esp + nb334_fizH2], xmm2

	;# M-H1 interaction 
	movapd xmm0, [esp + nb334_rinvMH1]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb334_rsqMH1] ;# xmm1=r 
	mulsd  xmm1, [esp + nb334_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb334_VFtab]
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb334_two]	;# two*Heps2 
	movapd xmm3, [esp + nb334_qqMH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm3 fijC 

    	addsd  xmm5, [esp + nb334_vctot]
    	movsd [esp + nb334_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb334_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1

	movapd xmm3, [esp + nb334_fjxH1]
	movapd xmm4, [esp + nb334_fjyH1]
	movapd xmm5, [esp + nb334_fjzH1]
	mulsd xmm0, [esp + nb334_dxMH1]
	mulsd xmm1, [esp + nb334_dyMH1]
	mulsd xmm2, [esp + nb334_dzMH1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb334_fixM]
	addsd xmm1, [esp + nb334_fiyM]
	addsd xmm2, [esp + nb334_fizM]
	movsd [esp + nb334_fjxH1], xmm3
	movsd [esp + nb334_fjyH1], xmm4
	movsd [esp + nb334_fjzH1], xmm5
	movsd [esp + nb334_fixM], xmm0
	movsd [esp + nb334_fiyM], xmm1
	movsd [esp + nb334_fizM], xmm2

	;# M-H2 interaction 
	movapd xmm0, [esp + nb334_rinvMH2]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb334_rsqMH2] ;# xmm1=r 
	mulsd  xmm1, [esp + nb334_tsc]
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb334_VFtab]
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb334_two]	;# two*Heps2 
	movapd xmm3, [esp + nb334_qqMH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm3 fijC 

    	addsd  xmm5, [esp + nb334_vctot]
    	movsd [esp + nb334_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb334_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb334_fjxH2]
	movapd xmm4, [esp + nb334_fjyH2]
	movapd xmm5, [esp + nb334_fjzH2]
	mulsd xmm0, [esp + nb334_dxMH2]
	mulsd xmm1, [esp + nb334_dyMH2]
	mulsd xmm2, [esp + nb334_dzMH2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb334_fixM]
	addsd xmm1, [esp + nb334_fiyM]
	addsd xmm2, [esp + nb334_fizM]
	movsd [esp + nb334_fjxH2], xmm3
	movsd [esp + nb334_fjyH2], xmm4
	movsd [esp + nb334_fjzH2], xmm5
	movsd [esp + nb334_fixM], xmm0
	movsd [esp + nb334_fiyM], xmm1
	movsd [esp + nb334_fizM], xmm2

	;# M-M interaction 
	movapd xmm0, [esp + nb334_rinvMM]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb334_rsqMM] ;# xmm1=r 
	mulsd  xmm1, [esp + nb334_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb334_VFtab]
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb334_two]	;# two*Heps2 
	movapd xmm3, [esp + nb334_qqMM]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm3 fijC 
	
    	addsd  xmm5, [esp + nb334_vctot]
    	movsd [esp + nb334_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb334_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb334_fjxM]
	movapd xmm4, [esp + nb334_fjyM]
	movapd xmm5, [esp + nb334_fjzM]
	mulsd xmm0, [esp + nb334_dxMM]
	mulsd xmm1, [esp + nb334_dyMM]
	mulsd xmm2, [esp + nb334_dzMM]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb334_fixM]
	addsd xmm1, [esp + nb334_fiyM]
	addsd xmm2, [esp + nb334_fizM]
	movsd [esp + nb334_fjxM], xmm3
	movsd [esp + nb334_fjyM], xmm4
	movsd [esp + nb334_fjzM], xmm5
	movsd [esp + nb334_fixM], xmm0
	movsd [esp + nb334_fiyM], xmm1
	movsd [esp + nb334_fizM], xmm2

	mov edi, [ebp + nb334_faction]

	movd eax, mm0

	;# Did all interactions - now update j forces 
	;# Step1 - merge forces
	movlpd xmm0, [esp + nb334_fjxO]
	movlpd xmm1, [esp + nb334_fjzO]
	movlpd xmm2, [esp + nb334_fjyH1]
	movlpd xmm3, [esp + nb334_fjxH2]
	movlpd xmm4, [esp + nb334_fjzH2]
	movlpd xmm5, [esp + nb334_fjyM]

	movhpd xmm0, [esp + nb334_fjyO]
	movhpd xmm1, [esp + nb334_fjxH1]
	movhpd xmm2, [esp + nb334_fjzH1]
	movhpd xmm3, [esp + nb334_fjyH2]
	movhpd xmm4, [esp + nb334_fjxM]
	movhpd xmm5, [esp + nb334_fjzM]

	movlps xmm6, [edi + eax*8]
	movhps xmm6, [edi + eax*8 + 8]
	movlps xmm7, [edi + eax*8 + 16]
	movhps xmm7, [edi + eax*8 + 24]
	addpd  xmm0, xmm6
	addpd  xmm1, xmm7
	movlps xmm6, [edi + eax*8 + 32]
	movhps xmm6, [edi + eax*8 + 40]
	movlps xmm7, [edi + eax*8 + 48]
	movhps xmm7, [edi + eax*8 + 56]
	addpd  xmm2, xmm6
	addpd  xmm3, xmm7
	movlps xmm6, [edi + eax*8 + 64]
	movhps xmm6, [edi + eax*8 + 72]
	movlps xmm7, [edi + eax*8 + 80]
	movhps xmm7, [edi + eax*8 + 88]
	addpd  xmm4, xmm6
	addpd  xmm5, xmm7
	
	movlpd [edi + eax*8],      xmm0
	movhpd [edi + eax*8 + 8],  xmm0
	movlpd [edi + eax*8 + 16], xmm1
	movhpd [edi + eax*8 + 24], xmm1
	movlpd [edi + eax*8 + 32], xmm2
	movhpd [edi + eax*8 + 40], xmm2
	movlpd [edi + eax*8 + 48], xmm3
	movhpd [edi + eax*8 + 56], xmm3
	movlpd [edi + eax*8 + 64], xmm4
	movhpd [edi + eax*8 + 72], xmm4
	movlpd [edi + eax*8 + 80], xmm5
	movhpd [edi + eax*8 + 88], xmm5

	
.nb334_updateouterdata:
	mov   ecx, [esp + nb334_ii3]
	mov   edi, [ebp + nb334_faction]
	mov   esi, [ebp + nb334_fshift]
	mov   edx, [esp + nb334_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb334_fixO]
	movapd xmm1, [esp + nb334_fiyO]
	movapd xmm2, [esp + nb334_fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	movapd xmm3, xmm0	
	movapd xmm4, xmm1	
	movapd xmm5, xmm2	

	;# increment i force 
	movsd  xmm3, [edi + ecx*8]
	movsd  xmm4, [edi + ecx*8 + 8]
	movsd  xmm5, [edi + ecx*8 + 16]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movsd  [edi + ecx*8],     xmm3
	movsd  [edi + ecx*8 + 8], xmm4
	movsd  [edi + ecx*8 + 16], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movapd xmm6, xmm0
	movsd xmm7, xmm2
	unpcklpd xmm6, xmm1

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb334_fixH1]
	movapd xmm1, [esp + nb334_fiyH1]
	movapd xmm2, [esp + nb334_fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	;# increment i force 
	movsd  xmm3, [edi + ecx*8 + 24]
	movsd  xmm4, [edi + ecx*8 + 32]
	movsd  xmm5, [edi + ecx*8 + 40]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movsd  [edi + ecx*8 + 24], xmm3
	movsd  [edi + ecx*8 + 32], xmm4
	movsd  [edi + ecx*8 + 40], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addsd xmm7, xmm2
	unpcklpd xmm0, xmm1
	addpd xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb334_fixH2]
	movapd xmm1, [esp + nb334_fiyH2]
	movapd xmm2, [esp + nb334_fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	movapd xmm3, xmm0	
	movapd xmm4, xmm1	
	movapd xmm5, xmm2	

	;# increment i force 
	movsd  xmm3, [edi + ecx*8 + 48]
	movsd  xmm4, [edi + ecx*8 + 56]
	movsd  xmm5, [edi + ecx*8 + 64]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movsd  [edi + ecx*8 + 48], xmm3
	movsd  [edi + ecx*8 + 56], xmm4
	movsd  [edi + ecx*8 + 64], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addsd xmm7, xmm2
	unpcklpd xmm0, xmm1
	addpd xmm6, xmm0

	;# accumulate Mi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb334_fixM]
	movapd xmm1, [esp + nb334_fiyM]
	movapd xmm2, [esp + nb334_fizM]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	movapd xmm3, xmm0	
	movapd xmm4, xmm1	
	movapd xmm5, xmm2	

	;# increment i force 
	movsd  xmm3, [edi + ecx*8 + 72]
	movsd  xmm4, [edi + ecx*8 + 80]
	movsd  xmm5, [edi + ecx*8 + 88]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movsd  [edi + ecx*8 + 72], xmm3
	movsd  [edi + ecx*8 + 80], xmm4
	movsd  [edi + ecx*8 + 88], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addsd xmm7, xmm2
	unpcklpd xmm0, xmm1
	addpd xmm6, xmm0

	;# increment fshift force 
	movlpd xmm3, [esi + edx*8]
	movhpd xmm3, [esi + edx*8 + 8]
	movsd  xmm4, [esi + edx*8 + 16]
	addpd  xmm3, xmm6
	addsd  xmm4, xmm7
	movlpd [esi + edx*8],      xmm3
	movhpd [esi + edx*8 + 8],  xmm3
	movsd  [esi + edx*8 + 16], xmm4

	;# get n from stack
	mov esi, [esp + nb334_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb334_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb334_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb334_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb334_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb334_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb334_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb334_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb334_n], esi
        jmp .nb334_outer
.nb334_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb334_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb334_end
        ;# non-zero, do one more workunit
        jmp   .nb334_threadloop
.nb334_end:
	emms

	mov eax, [esp + nb334_nouter]
	mov ebx, [esp + nb334_ninner]
	mov ecx, [ebp + nb334_outeriter]
	mov edx, [ebp + nb334_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb334_salign]
	add esp, eax
	add esp, 1800
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret


	
.globl nb_kernel334nf_ia32_sse2
.globl _nb_kernel334nf_ia32_sse2
nb_kernel334nf_ia32_sse2:	
_nb_kernel334nf_ia32_sse2:	
.equiv          nb334nf_p_nri,          8
.equiv          nb334nf_iinr,           12
.equiv          nb334nf_jindex,         16
.equiv          nb334nf_jjnr,           20
.equiv          nb334nf_shift,          24
.equiv          nb334nf_shiftvec,       28
.equiv          nb334nf_fshift,         32
.equiv          nb334nf_gid,            36
.equiv          nb334nf_pos,            40
.equiv          nb334nf_faction,        44
.equiv          nb334nf_charge,         48
.equiv          nb334nf_p_facel,        52
.equiv          nb334nf_argkrf,         56
.equiv          nb334nf_argcrf,         60
.equiv          nb334nf_Vc,             64
.equiv          nb334nf_type,           68
.equiv          nb334nf_p_ntype,        72
.equiv          nb334nf_vdwparam,       76
.equiv          nb334nf_Vvdw,           80
.equiv          nb334nf_p_tabscale,     84
.equiv          nb334nf_VFtab,          88
.equiv          nb334nf_invsqrta,       92
.equiv          nb334nf_dvda,           96
.equiv          nb334nf_p_gbtabscale,   100
.equiv          nb334nf_GBtab,          104
.equiv          nb334nf_p_nthreads,     108
.equiv          nb334nf_count,          112
.equiv          nb334nf_mtx,            116
.equiv          nb334nf_outeriter,      120
.equiv          nb334nf_inneriter,      124
.equiv          nb334nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb334nf_ixO,            0
.equiv          nb334nf_iyO,            16
.equiv          nb334nf_izO,            32
.equiv          nb334nf_ixH1,           48
.equiv          nb334nf_iyH1,           64
.equiv          nb334nf_izH1,           80
.equiv          nb334nf_ixH2,           96
.equiv          nb334nf_iyH2,           112
.equiv          nb334nf_izH2,           128
.equiv          nb334nf_ixM,            144
.equiv          nb334nf_iyM,            160
.equiv          nb334nf_izM,            176
.equiv          nb334nf_jxO,            192
.equiv          nb334nf_jyO,            208
.equiv          nb334nf_jzO,            224
.equiv          nb334nf_jxH1,           240
.equiv          nb334nf_jyH1,           256
.equiv          nb334nf_jzH1,           272
.equiv          nb334nf_jxH2,           288
.equiv          nb334nf_jyH2,           304
.equiv          nb334nf_jzH2,           320
.equiv          nb334nf_jxM,            336
.equiv          nb334nf_jyM,            352
.equiv          nb334nf_jzM,            368
.equiv          nb334nf_qqMM,           384
.equiv          nb334nf_qqMH,           400
.equiv          nb334nf_qqHH,           416
.equiv          nb334nf_tsc,            432
.equiv          nb334nf_c6,             448
.equiv          nb334nf_c12,            464
.equiv          nb334nf_vctot,          480
.equiv          nb334nf_Vvdwtot,        496
.equiv          nb334nf_half,           512
.equiv          nb334nf_three,          528
.equiv          nb334nf_rsqOO,          544
.equiv          nb334nf_rsqH1H1,        560
.equiv          nb334nf_rsqH1H2,        576
.equiv          nb334nf_rsqH1M,         592
.equiv          nb334nf_rsqH2H1,        608
.equiv          nb334nf_rsqH2H2,        624
.equiv          nb334nf_rsqH2M,         640
.equiv          nb334nf_rsqMH1,         656
.equiv          nb334nf_rsqMH2,         672
.equiv          nb334nf_rsqMM,          688
.equiv          nb334nf_rinvOO,         704
.equiv          nb334nf_rinvH1H1,       720
.equiv          nb334nf_rinvH1H2,       736
.equiv          nb334nf_rinvH1M,        752
.equiv          nb334nf_rinvH2H1,       768
.equiv          nb334nf_rinvH2H2,       784
.equiv          nb334nf_rinvH2M,        800
.equiv          nb334nf_rinvMH1,        816
.equiv          nb334nf_rinvMH2,        832
.equiv          nb334nf_rinvMM,         848
.equiv          nb334nf_is3,            864
.equiv          nb334nf_ii3,            868
.equiv          nb334nf_innerjjnr,      872
.equiv          nb334nf_innerk,         876
.equiv          nb334nf_n,              880
.equiv          nb334nf_nn1,            884
.equiv          nb334nf_nri,            888
.equiv          nb334nf_nouter,         892
.equiv          nb334nf_ninner,         896
.equiv          nb334nf_salign,         900
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 904		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb334nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb334nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb334nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb334nf_nouter], eax
	mov [esp + nb334nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb334nf_half], eax
	mov [esp + nb334nf_half+4], ebx
	movsd xmm1, [esp + nb334nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb334nf_half], xmm1
	movapd [esp + nb334nf_three], xmm3
	mov eax, [ebp + nb334nf_p_tabscale]
	movsd xmm3, [eax]

	shufpd xmm3, xmm3, 0
	movapd [esp + nb334nf_tsc],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb334nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb334nf_charge]
	movsd xmm3, [edx + ebx*8 + 24]	
	movsd xmm4, xmm3	
	movsd xmm5, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb334nf_p_facel]
	movsd xmm6, [esi]
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm5
	mulsd  xmm5, xmm5
	mulsd  xmm3, xmm6
	mulsd  xmm4, xmm6
	mulsd  xmm5, xmm6
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb334nf_qqMM], xmm3
	movapd [esp + nb334nf_qqMH], xmm4
	movapd [esp + nb334nf_qqHH], xmm5
		
	xorpd xmm0, xmm0
	mov   edx, [ebp + nb334nf_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov edi, [ebp + nb334nf_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   eax, [ebp + nb334nf_vdwparam]
	movlpd xmm0, [eax + edx*8]
	movlpd xmm1, [eax + edx*8 + 8]
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	movapd [esp + nb334nf_c6], xmm0
	movapd [esp + nb334nf_c12], xmm1

.nb334nf_threadloop:
        mov   esi, [ebp + nb334nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb334nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb334nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb334nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb334nf_n], eax
        mov [esp + nb334nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb334nf_outerstart
        jmp .nb334nf_end

.nb334nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb334nf_nouter]
	mov [esp + nb334nf_nouter], ebx

.nb334nf_outer:
	mov   eax, [ebp + nb334nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb334nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb334nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb334nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb334nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb334nf_ii3], ebx		

	addsd xmm3, [eax + ebx*8] 	;# ox
	addsd xmm4, [eax + ebx*8 + 8] 	;# oy
	addsd xmm5, [eax + ebx*8 + 16]	;# oz	
	addsd xmm6, [eax + ebx*8 + 24] 	;# h1x
	addsd xmm7, [eax + ebx*8 + 32] 	;# h1y
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	shufpd xmm7, xmm7, 0
	movapd [esp + nb334nf_ixO], xmm3
	movapd [esp + nb334nf_iyO], xmm4
	movapd [esp + nb334nf_izO], xmm5
	movapd [esp + nb334nf_ixH1], xmm6
	movapd [esp + nb334nf_iyH1], xmm7

	movsd xmm6, xmm2
	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm6, [eax + ebx*8 + 40] ;# h1z
	addsd xmm0, [eax + ebx*8 + 48] ;# h2x
	addsd xmm1, [eax + ebx*8 + 56] ;# h2y
	addsd xmm2, [eax + ebx*8 + 64] ;# h2z
	addsd xmm3, [eax + ebx*8 + 72] ;# mx
	addsd xmm4, [eax + ebx*8 + 80] ;# my
	addsd xmm5, [eax + ebx*8 + 88] ;# mz

	shufpd xmm6, xmm6, 0
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb334nf_izH1], xmm6
	movapd [esp + nb334nf_ixH2], xmm0
	movapd [esp + nb334nf_iyH2], xmm1
	movapd [esp + nb334nf_izH2], xmm2
	movapd [esp + nb334nf_ixM], xmm3
	movapd [esp + nb334nf_iyM], xmm4
	movapd [esp + nb334nf_izM], xmm5

	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [esp + nb334nf_vctot], xmm4
	movapd [esp + nb334nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb334nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb334nf_pos] 
	mov   edi, [ebp + nb334nf_faction]	
	mov   eax, [ebp + nb334nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb334nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb334nf_ninner]
	mov   [esp + nb334nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb334nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb334nf_unroll_loop
	jmp   .nb334nf_checksingle
.nb334nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb334nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	
	add dword ptr [esp + nb334nf_innerjjnr], 8 ;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb334nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	
	;# move j coordinates to local temp variables 
	;# load ox, oy, oz, h1x
	movlpd xmm0, [esi + eax*8]
	movlpd xmm2, [esi + ebx*8]
	movhpd xmm0, [esi + eax*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 8]
	movlpd xmm3, [esi + eax*8 + 16]
	movlpd xmm5, [esi + ebx*8 + 16]
	movhpd xmm3, [esi + eax*8 + 24]
	movhpd xmm5, [esi + ebx*8 + 24]
	movapd xmm1, xmm0 
	movapd xmm4, xmm3
	unpcklpd xmm0, xmm2 ;# ox 
	unpckhpd xmm1, xmm2 ;# oy
	unpcklpd xmm3, xmm5 ;# ox 
	unpckhpd xmm4, xmm5 ;# oy
	movapd 	[esp + nb334nf_jxO], xmm0
	movapd 	[esp + nb334nf_jyO], xmm1
	movapd 	[esp + nb334nf_jzO], xmm3
	movapd 	[esp + nb334nf_jxH1], xmm4
	
	;# load h1y, h1z, h2x, h2y 
	movlpd xmm0, [esi + eax*8 + 32]
	movlpd xmm2, [esi + ebx*8 + 32]
	movhpd xmm0, [esi + eax*8 + 40]
	movhpd xmm2, [esi + ebx*8 + 40]
	movlpd xmm3, [esi + eax*8 + 48]
	movlpd xmm5, [esi + ebx*8 + 48]
	movhpd xmm3, [esi + eax*8 + 56]
	movhpd xmm5, [esi + ebx*8 + 56]
	movapd xmm1, xmm0 
	movapd xmm4, xmm3
	unpcklpd xmm0, xmm2 ;# h1y
	unpckhpd xmm1, xmm2 ;# h1z
	unpcklpd xmm3, xmm5 ;# h2x
	unpckhpd xmm4, xmm5 ;# h2y
	movapd 	[esp + nb334nf_jyH1], xmm0
	movapd 	[esp + nb334nf_jzH1], xmm1
	movapd 	[esp + nb334nf_jxH2], xmm3
	movapd 	[esp + nb334nf_jyH2], xmm4
	
	;# load h2z, mx, my, mz
	movlpd xmm0, [esi + eax*8 + 64]
	movlpd xmm2, [esi + ebx*8 + 64]
	movhpd xmm0, [esi + eax*8 + 72]
	movhpd xmm2, [esi + ebx*8 + 72]
	movlpd xmm3, [esi + eax*8 + 80]
	movlpd xmm5, [esi + ebx*8 + 80]
	movhpd xmm3, [esi + eax*8 + 88]
	movhpd xmm5, [esi + ebx*8 + 88]
	movapd xmm1, xmm0 
	movapd xmm4, xmm3
	unpcklpd xmm0, xmm2 ;# h2z
	unpckhpd xmm1, xmm2 ;# mx
	unpcklpd xmm3, xmm5 ;# my
	unpckhpd xmm4, xmm5 ;# mz
	movapd 	[esp + nb334nf_jzH2], xmm0
	movapd 	[esp + nb334nf_jxM], xmm1
	movapd 	[esp + nb334nf_jyM], xmm3
	movapd 	[esp + nb334nf_jzM], xmm4
	
	;# start calculating pairwise distances
	movapd xmm0, [esp + nb334nf_ixO]
	movapd xmm1, [esp + nb334nf_iyO]
	movapd xmm2, [esp + nb334nf_izO]
	movapd xmm3, [esp + nb334nf_ixH1]
	movapd xmm4, [esp + nb334nf_iyH1]
	movapd xmm5, [esp + nb334nf_izH1]
	subpd  xmm0, [esp + nb334nf_jxO]
	subpd  xmm1, [esp + nb334nf_jyO]
	subpd  xmm2, [esp + nb334nf_jzO]
	subpd  xmm3, [esp + nb334nf_jxH1]
	subpd  xmm4, [esp + nb334nf_jyH1]
	subpd  xmm5, [esp + nb334nf_jzH1]
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb334nf_rsqOO], xmm0
	movapd [esp + nb334nf_rsqH1H1], xmm3

	movapd xmm0, [esp + nb334nf_ixH1]
	movapd xmm1, [esp + nb334nf_iyH1]
	movapd xmm2, [esp + nb334nf_izH1]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [esp + nb334nf_jxH2]
	subpd  xmm1, [esp + nb334nf_jyH2]
	subpd  xmm2, [esp + nb334nf_jzH2]
	subpd  xmm3, [esp + nb334nf_jxM]
	subpd  xmm4, [esp + nb334nf_jyM]
	subpd  xmm5, [esp + nb334nf_jzM]
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb334nf_rsqH1H2], xmm0
	movapd [esp + nb334nf_rsqH1M], xmm3

	movapd xmm0, [esp + nb334nf_ixH2]
	movapd xmm1, [esp + nb334nf_iyH2]
	movapd xmm2, [esp + nb334nf_izH2]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [esp + nb334nf_jxH1]
	subpd  xmm1, [esp + nb334nf_jyH1]
	subpd  xmm2, [esp + nb334nf_jzH1]
	subpd  xmm3, [esp + nb334nf_jxH2]
	subpd  xmm4, [esp + nb334nf_jyH2]
	subpd  xmm5, [esp + nb334nf_jzH2]
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb334nf_rsqH2H1], xmm0
	movapd [esp + nb334nf_rsqH2H2], xmm3

	movapd xmm0, [esp + nb334nf_ixH2]
	movapd xmm1, [esp + nb334nf_iyH2]
	movapd xmm2, [esp + nb334nf_izH2]
	movapd xmm3, [esp + nb334nf_ixM]
	movapd xmm4, [esp + nb334nf_iyM]
	movapd xmm5, [esp + nb334nf_izM]
	subpd  xmm0, [esp + nb334nf_jxM]
	subpd  xmm1, [esp + nb334nf_jyM]
	subpd  xmm2, [esp + nb334nf_jzM]
	subpd  xmm3, [esp + nb334nf_jxH1]
	subpd  xmm4, [esp + nb334nf_jyH1]
	subpd  xmm5, [esp + nb334nf_jzH1]
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm4, xmm3
	addpd  xmm4, xmm5
	movapd [esp + nb334nf_rsqH2M], xmm0
	movapd [esp + nb334nf_rsqMH1], xmm4

	movapd xmm0, [esp + nb334nf_ixM]
	movapd xmm1, [esp + nb334nf_iyM]
	movapd xmm2, [esp + nb334nf_izM]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [esp + nb334nf_jxH2]
	subpd  xmm1, [esp + nb334nf_jyH2]
	subpd  xmm2, [esp + nb334nf_jzH2]
	subpd  xmm3, [esp + nb334nf_jxM]
	subpd  xmm4, [esp + nb334nf_jyM]
	subpd  xmm5, [esp + nb334nf_jzM]
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm4, xmm3
	addpd  xmm4, xmm5
	movapd [esp + nb334nf_rsqMH2], xmm0
	movapd [esp + nb334nf_rsqMM], xmm4
	
	;# Invsqrt form rsq M-H2 (rsq in xmm0) and MM (rsq in xmm4) 
	cvtpd2ps xmm1, xmm0
	cvtpd2ps xmm5, xmm4
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1   ;# luA
	cvtps2pd xmm5, xmm5   ;# luB
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb334nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb334nf_half] ;# iter1 
	mulpd   xmm7, [esp + nb334nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb334nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb334nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb334nf_half] ;# rinv 
	movapd [esp + nb334nf_rinvMH2], xmm1
	movapd [esp + nb334nf_rinvMM], xmm5

	movapd xmm0, [esp + nb334nf_rsqOO]
	movapd xmm4, [esp + nb334nf_rsqH1H1]	
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb334nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb334nf_half] ;# iter1 of  
	mulpd   xmm7, [esp + nb334nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb334nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb334nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb334nf_half] ;# rinv
	movapd [esp + nb334nf_rinvOO], xmm1
	movapd [esp + nb334nf_rinvH1H1], xmm5

	movapd xmm0, [esp + nb334nf_rsqH1H2]
	movapd xmm4, [esp + nb334nf_rsqH1M]	
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb334nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb334nf_half] ;# iter1 
	mulpd   xmm7, [esp + nb334nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb334nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb334nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb334nf_half] ;# rinv 
	movapd [esp + nb334nf_rinvH1H2], xmm1
	movapd [esp + nb334nf_rinvH1M], xmm5

	movapd xmm0, [esp + nb334nf_rsqH2H1]
	movapd xmm4, [esp + nb334nf_rsqH2H2]	
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb334nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb334nf_half] ;# iter1a 
	mulpd   xmm7, [esp + nb334nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb334nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb334nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb334nf_half] ;# rinv 
	movapd [esp + nb334nf_rinvH2H1], xmm1
	movapd [esp + nb334nf_rinvH2H2], xmm5

	movapd xmm0, [esp + nb334nf_rsqMH1]
	movapd xmm4, [esp + nb334nf_rsqH2M]	
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb334nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb334nf_half] ;# iter1a 
	mulpd   xmm7, [esp + nb334nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb334nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb334nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb334nf_half] ;# rinv 
	movapd [esp + nb334nf_rinvMH1], xmm1
	movapd [esp + nb334nf_rinvH2M], xmm5

	;# start with OO interaction 
	movapd xmm0, [esp + nb334nf_rinvOO]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb334nf_rsqOO] ;# xmm1=r 
	mulpd  xmm1, [esp + nb334nf_tsc]

	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	movd mm0, eax	
	movd mm1, ebx
	mov  esi, [ebp + nb334nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 
	lea   ebx, [ebx + ebx*2]	;# idx*=3 (12 total now) 

	;# Dispersion 
	movlpd xmm4, [esi + eax*8 + 32]	;# Y1
	movlpd xmm3, [esi + ebx*8 + 32]	;# Y2
	movhpd xmm4, [esi + eax*8 + 40]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 40]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 48]	;# G1
	movlpd xmm3, [esi + ebx*8 + 48]	;# G2
	movhpd xmm6, [esi + eax*8 + 56]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 56]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# Dispersion table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb334nf_c6]
	mulpd  xmm5, xmm4	 ;# Vvdw6 
	addpd  xmm5, [esp + nb334nf_Vvdwtot]
	movapd [esp + nb334nf_Vvdwtot], xmm5

	;# Repulsion 
	movlpd xmm4, [esi + eax*8 + 64]	;# Y1
	movlpd xmm3, [esi + ebx*8 + 64]	;# Y2
	movhpd xmm4, [esi + eax*8 + 72]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 72]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 80]	;# G1
	movlpd xmm3, [esi + ebx*8 + 80]	;# G2
	movhpd xmm6, [esi + eax*8 + 88]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 88]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# Dispersion table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb334nf_c12]
	mulpd  xmm5, xmm4 ;# Vvdw12 
	addpd  xmm5, [esp + nb334nf_Vvdwtot]
	movapd [esp + nb334nf_Vvdwtot], xmm5

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb334nf_rinvH1H1]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb334nf_rsqH1H1] ;# xmm1=r 
	mulpd  xmm1, [esp + nb334nf_tsc]

	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb334nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 
	lea   ebx, [ebx + ebx*2]	;# idx*=3 (12 total now) 

	movlpd xmm4, [esi + eax*8]	;# Y1	
	movlpd xmm3, [esi + ebx*8]	;# Y2
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb334nf_qqHH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    	addpd  xmm5, [esp + nb334nf_vctot]
    	movapd [esp + nb334nf_vctot], xmm5

	;# H1-H2 interaction  
	movapd xmm0, [esp + nb334nf_rinvH1H2]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb334nf_rsqH1H2] ;# xmm1=r 
	mulpd  xmm1, [esp + nb334nf_tsc]
	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb334nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 
	lea   ebx, [ebx + ebx*2]	;# idx*=3 (12 total now) 

	movlpd xmm4, [esi + eax*8]	;# Y1	
	movlpd xmm3, [esi + ebx*8]	;# Y2
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb334nf_qqHH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    	addpd  xmm5, [esp + nb334nf_vctot]
    	movapd [esp + nb334nf_vctot], xmm5

	;# H1-M interaction 
	movapd xmm0, [esp + nb334nf_rinvH1M]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb334nf_rsqH1M] ;# xmm1=r 
	mulpd  xmm1, [esp + nb334nf_tsc]
	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb334nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 
	lea   ebx, [ebx + ebx*2]	;# idx*=3 (12 total now) 

	movlpd xmm4, [esi + eax*8]	;# Y1	
	movlpd xmm3, [esi + ebx*8]	;# Y2
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb334nf_qqMH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addpd  xmm5, [esp + nb334nf_vctot]
    	movapd [esp + nb334nf_vctot], xmm5

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb334nf_rinvH2H1]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb334nf_rsqH2H1] ;# xmm1=r 
	mulpd  xmm1, [esp + nb334nf_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb334nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 
	lea   ebx, [ebx + ebx*2]	;# idx*=3 (12 total now) 

	movlpd xmm4, [esi + eax*8]	;# Y1	
	movlpd xmm3, [esi + ebx*8]	;# Y2
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb334nf_qqHH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addpd  xmm5, [esp + nb334nf_vctot]
    	movapd [esp + nb334nf_vctot], xmm5

	;# H2-H2 interaction 
	movapd xmm0, [esp + nb334nf_rinvH2H2]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb334nf_rsqH2H2] ;# xmm1=r 
	mulpd  xmm1, [esp + nb334nf_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb334nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 
	lea   ebx, [ebx + ebx*2]	;# idx*=3 (12 total now) 

	movlpd xmm4, [esi + eax*8]	;# Y1	
	movlpd xmm3, [esi + ebx*8]	;# Y2
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb334nf_qqHH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addpd  xmm5, [esp + nb334nf_vctot]
    	movapd [esp + nb334nf_vctot], xmm5

	;# H2-M interaction 
	movapd xmm0, [esp + nb334nf_rinvH2M]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb334nf_rsqH2M] ;# xmm1=r 
	mulpd  xmm1, [esp + nb334nf_tsc]
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb334nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 
	lea   ebx, [ebx + ebx*2]	;# idx*=3 (12 total now) 

	movlpd xmm4, [esi + eax*8]	;# Y1	
	movlpd xmm3, [esi + ebx*8]	;# Y2
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb334nf_qqMH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addpd  xmm5, [esp + nb334nf_vctot]
    	movapd [esp + nb334nf_vctot], xmm5

	;# M-H1 interaction 
	movapd xmm0, [esp + nb334nf_rinvMH1]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb334nf_rsqMH1] ;# xmm1=r 
	mulpd  xmm1, [esp + nb334nf_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb334nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 
	lea   ebx, [ebx + ebx*2]	;# idx*=3 (12 total now) 

	movlpd xmm4, [esi + eax*8]	;# Y1	
	movlpd xmm3, [esi + ebx*8]	;# Y2
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb334nf_qqMH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addpd  xmm5, [esp + nb334nf_vctot]
    	movapd [esp + nb334nf_vctot], xmm5

	;# M-H2 interaction 
	movapd xmm0, [esp + nb334nf_rinvMH2]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb334nf_rsqMH2] ;# xmm1=r 
	mulpd  xmm1, [esp + nb334nf_tsc]
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb334nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 
	lea   ebx, [ebx + ebx*2]	;# idx*=3 (12 total now) 

	movlpd xmm4, [esi + eax*8]	;# Y1	
	movlpd xmm3, [esi + ebx*8]	;# Y2
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb334nf_qqMH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addpd  xmm5, [esp + nb334nf_vctot]
    	movapd [esp + nb334nf_vctot], xmm5

	;# M-M interaction 
	movapd xmm0, [esp + nb334nf_rinvMM]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb334nf_rsqMM] ;# xmm1=r 
	mulpd  xmm1, [esp + nb334nf_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb334nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 
	lea   ebx, [ebx + ebx*2]	;# idx*=3 (12 total now) 

	movlpd xmm4, [esi + eax*8]	;# Y1	
	movlpd xmm3, [esi + ebx*8]	;# Y2
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb334nf_qqMM]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	
    	addpd  xmm5, [esp + nb334nf_vctot]
    	movapd [esp + nb334nf_vctot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [esp + nb334nf_innerk],  2
	jl    .nb334nf_checksingle
	jmp   .nb334nf_unroll_loop
.nb334nf_checksingle:
	mov   edx, [esp + nb334nf_innerk]
	and   edx, 1
	jnz   .nb334nf_dosingle
	jmp   .nb334nf_updateouterdata
.nb334nf_dosingle:
	mov   edx, [esp + nb334nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]

	mov esi, [ebp + nb334nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	
	;# move j coordinates to local temp variables 
	;# load ox, oy, oz, h1x
	movlpd xmm0, [esi + eax*8] 
	movhpd xmm0, [esi + eax*8 + 8]
	movlpd xmm1, [esi + eax*8 + 16]
	movhpd xmm1, [esi + eax*8 + 24]
	movlpd xmm2, [esi + eax*8 + 32]
	movhpd xmm2, [esi + eax*8 + 40]
	movlpd xmm3, [esi + eax*8 + 48]
	movhpd xmm3, [esi + eax*8 + 56]
	movlpd xmm4, [esi + eax*8 + 64]
	movhpd xmm4, [esi + eax*8 + 72]
	movlpd xmm5, [esi + eax*8 + 80]
	movhpd xmm5, [esi + eax*8 + 88]
	movsd  [esp + nb334nf_jxO], xmm0
	movsd  [esp + nb334nf_jzO], xmm1
	movsd  [esp + nb334nf_jyH1], xmm2
	movsd  [esp + nb334nf_jxH2], xmm3
	movsd  [esp + nb334nf_jzH2], xmm4
	movsd  [esp + nb334nf_jyM], xmm5
	unpckhpd xmm0, xmm0
	unpckhpd xmm1, xmm1
	unpckhpd xmm2, xmm2
	unpckhpd xmm3, xmm3
	unpckhpd xmm4, xmm4
	unpckhpd xmm5, xmm5
	movsd  [esp + nb334nf_jyO], xmm0
	movsd  [esp + nb334nf_jxH1], xmm1
	movsd  [esp + nb334nf_jzH1], xmm2
	movsd  [esp + nb334nf_jyH2], xmm3
	movsd  [esp + nb334nf_jxM], xmm4
	movsd  [esp + nb334nf_jzM], xmm5

	;# start calculating pairwise distances
	movapd xmm0, [esp + nb334nf_ixO]
	movapd xmm1, [esp + nb334nf_iyO]
	movapd xmm2, [esp + nb334nf_izO]
	movapd xmm3, [esp + nb334nf_ixH1]
	movapd xmm4, [esp + nb334nf_iyH1]
	movapd xmm5, [esp + nb334nf_izH1]
	subsd  xmm0, [esp + nb334nf_jxO]
	subsd  xmm1, [esp + nb334nf_jyO]
	subsd  xmm2, [esp + nb334nf_jzO]
	subsd  xmm3, [esp + nb334nf_jxH1]
	subsd  xmm4, [esp + nb334nf_jyH1]
	subsd  xmm5, [esp + nb334nf_jzH1]
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb334nf_rsqOO], xmm0
	movapd [esp + nb334nf_rsqH1H1], xmm3

	movapd xmm0, [esp + nb334nf_ixH1]
	movapd xmm1, [esp + nb334nf_iyH1]
	movapd xmm2, [esp + nb334nf_izH1]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [esp + nb334nf_jxH2]
	subsd  xmm1, [esp + nb334nf_jyH2]
	subsd  xmm2, [esp + nb334nf_jzH2]
	subsd  xmm3, [esp + nb334nf_jxM]
	subsd  xmm4, [esp + nb334nf_jyM]
	subsd  xmm5, [esp + nb334nf_jzM]
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb334nf_rsqH1H2], xmm0
	movapd [esp + nb334nf_rsqH1M], xmm3

	movapd xmm0, [esp + nb334nf_ixH2]
	movapd xmm1, [esp + nb334nf_iyH2]
	movapd xmm2, [esp + nb334nf_izH2]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [esp + nb334nf_jxH1]
	subsd  xmm1, [esp + nb334nf_jyH1]
	subsd  xmm2, [esp + nb334nf_jzH1]
	subsd  xmm3, [esp + nb334nf_jxH2]
	subsd  xmm4, [esp + nb334nf_jyH2]
	subsd  xmm5, [esp + nb334nf_jzH2]
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb334nf_rsqH2H1], xmm0
	movapd [esp + nb334nf_rsqH2H2], xmm3

	movapd xmm0, [esp + nb334nf_ixH2]
	movapd xmm1, [esp + nb334nf_iyH2]
	movapd xmm2, [esp + nb334nf_izH2]
	movapd xmm3, [esp + nb334nf_ixM]
	movapd xmm4, [esp + nb334nf_iyM]
	movapd xmm5, [esp + nb334nf_izM]
	subsd  xmm0, [esp + nb334nf_jxM]
	subsd  xmm1, [esp + nb334nf_jyM]
	subsd  xmm2, [esp + nb334nf_jzM]
	subsd  xmm3, [esp + nb334nf_jxH1]
	subsd  xmm4, [esp + nb334nf_jyH1]
	subsd  xmm5, [esp + nb334nf_jzH1]
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm4, xmm3
	addsd  xmm4, xmm5
	movapd [esp + nb334nf_rsqH2M], xmm0
	movapd [esp + nb334nf_rsqMH1], xmm4

	movapd xmm0, [esp + nb334nf_ixM]
	movapd xmm1, [esp + nb334nf_iyM]
	movapd xmm2, [esp + nb334nf_izM]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [esp + nb334nf_jxH2]
	subsd  xmm1, [esp + nb334nf_jyH2]
	subsd  xmm2, [esp + nb334nf_jzH2]
	subsd  xmm3, [esp + nb334nf_jxM]
	subsd  xmm4, [esp + nb334nf_jyM]
	subsd  xmm5, [esp + nb334nf_jzM]
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm4, xmm3
	addsd  xmm4, xmm5
	movapd [esp + nb334nf_rsqMH2], xmm0
	movapd [esp + nb334nf_rsqMM], xmm4

	;# Invsqrt form rsq M-H2 (rsq in xmm0) and MM (rsq in xmm4) 
	cvtsd2ss xmm1, xmm0
	cvtsd2ss xmm5, xmm4
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1   ;# luA
	cvtss2sd xmm5, xmm5   ;# luB
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb334nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb334nf_half] ;# iter1 
	mulsd   xmm7, [esp + nb334nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb334nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb334nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb334nf_half] ;# rinv 
	movapd [esp + nb334nf_rinvMH2], xmm1
	movapd [esp + nb334nf_rinvMM], xmm5

	movapd xmm0, [esp + nb334nf_rsqOO]
	movapd xmm4, [esp + nb334nf_rsqH1H1]	
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb334nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb334nf_half] ;# iter1 of  
	mulsd   xmm7, [esp + nb334nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb334nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb334nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb334nf_half] ;# rinv
	movapd [esp + nb334nf_rinvOO], xmm1
	movapd [esp + nb334nf_rinvH1H1], xmm5

	movapd xmm0, [esp + nb334nf_rsqH1H2]
	movapd xmm4, [esp + nb334nf_rsqH1M]	
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb334nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb334nf_half] ;# iter1 
	mulsd   xmm7, [esp + nb334nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb334nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb334nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb334nf_half] ;# rinv 
	movapd [esp + nb334nf_rinvH1H2], xmm1
	movapd [esp + nb334nf_rinvH1M], xmm5

	movapd xmm0, [esp + nb334nf_rsqH2H1]
	movapd xmm4, [esp + nb334nf_rsqH2H2]	
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb334nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb334nf_half] ;# iter1a 
	mulsd   xmm7, [esp + nb334nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb334nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb334nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb334nf_half] ;# rinv 
	movapd [esp + nb334nf_rinvH2H1], xmm1
	movapd [esp + nb334nf_rinvH2H2], xmm5

	movapd xmm0, [esp + nb334nf_rsqMH1]
	movapd xmm4, [esp + nb334nf_rsqH2M]	
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb334nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb334nf_half] ;# iter1a 
	mulsd   xmm7, [esp + nb334nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb334nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb334nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb334nf_half] ;# rinv 
	movapd [esp + nb334nf_rinvMH1], xmm1
	movapd [esp + nb334nf_rinvH2M], xmm5

	;# start with OO interaction 
	movapd xmm0, [esp + nb334nf_rinvOO]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb334nf_rsqOO] ;# xmm1=r 
	mulsd  xmm1, [esp + nb334nf_tsc]

	movd mm0, eax	


	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 

		
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb334nf_VFtab]
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 

	;# Dispersion 
	movsd xmm4, [esi + eax*8 + 32]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 40]	;# F1 	
	movsd xmm6, [esi + eax*8 + 48]	;# G1 	
	movsd xmm7, [esi + eax*8 + 56]	;# H1 	
	;# Dispersion table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb334nf_c6]
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	addsd  xmm5, [esp + nb334nf_Vvdwtot]
	movsd [esp + nb334nf_Vvdwtot], xmm5

	;# Repulsion 
	movsd xmm4, [esi + eax*8 + 64]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 72]	;# F1 	
	movsd xmm6, [esi + eax*8 + 80]	;# G1 	
	movsd xmm7, [esi + eax*8 + 88]	;# H1 	

	;# Dispersion table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb334nf_c12]
	mulsd  xmm5, xmm4 ;# Vvdw12 


	addsd  xmm5, [esp + nb334nf_Vvdwtot]
	movsd [esp + nb334nf_Vvdwtot], xmm5

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb334nf_rinvH1H1]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb334nf_rsqH1H1] ;# xmm1=r 
	mulsd  xmm1, [esp + nb334nf_tsc]

	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb334nf_VFtab]
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb334nf_qqHH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  


    	addsd  xmm5, [esp + nb334nf_vctot]
    	movsd [esp + nb334nf_vctot], xmm5

	;# H1-H2 interaction  
	movapd xmm0, [esp + nb334nf_rinvH1H2]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb334nf_rsqH1H2] ;# xmm1=r 
	mulsd  xmm1, [esp + nb334nf_tsc]
	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb334nf_VFtab]
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	

	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb334nf_qqHH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    	addsd  xmm5, [esp + nb334nf_vctot]
    	movsd [esp + nb334nf_vctot], xmm5

	;# H1-M interaction 
	movapd xmm0, [esp + nb334nf_rinvH1M]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb334nf_rsqH1M] ;# xmm1=r 
	mulsd  xmm1, [esp + nb334nf_tsc]
	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb334nf_VFtab]
	lea  eax, [eax + eax*2]	;# idx*=3 (12 total now) 
	
	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb334nf_qqMH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addsd  xmm5, [esp + nb334nf_vctot]
    	movsd [esp + nb334nf_vctot], xmm5

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb334nf_rinvH2H1]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb334nf_rsqH2H1] ;# xmm1=r 
	mulsd  xmm1, [esp + nb334nf_tsc]	

	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb334nf_VFtab]
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb334nf_qqHH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addsd  xmm5, [esp + nb334nf_vctot]
    	movsd [esp + nb334nf_vctot], xmm5

	;# H2-H2 interaction 
	movapd xmm0, [esp + nb334nf_rinvH2H2]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb334nf_rsqH2H2] ;# xmm1=r 
	mulsd  xmm1, [esp + nb334nf_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb334nf_VFtab]
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	

	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb334nf_qqHH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addsd  xmm5, [esp + nb334nf_vctot]
    	movsd [esp + nb334nf_vctot], xmm5

	;# H2-M interaction 
	movapd xmm0, [esp + nb334nf_rinvH2M]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb334nf_rsqH2M] ;# xmm1=r 
	mulsd  xmm1, [esp + nb334nf_tsc]
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb334nf_VFtab]
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	

	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb334nf_qqMH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addsd  xmm5, [esp + nb334nf_vctot]
    	movsd [esp + nb334nf_vctot], xmm5

	;# M-H1 interaction 
	movapd xmm0, [esp + nb334nf_rinvMH1]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb334nf_rsqMH1] ;# xmm1=r 
	mulsd  xmm1, [esp + nb334nf_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb334nf_VFtab]
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb334nf_qqMH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addsd  xmm5, [esp + nb334nf_vctot]
    	movsd [esp + nb334nf_vctot], xmm5

	;# M-H2 interaction 
	movapd xmm0, [esp + nb334nf_rinvMH2]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb334nf_rsqMH2] ;# xmm1=r 
	mulsd  xmm1, [esp + nb334nf_tsc]
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb334nf_VFtab]
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb334nf_qqMH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addsd  xmm5, [esp + nb334nf_vctot]
    	movsd [esp + nb334nf_vctot], xmm5

	;# M-M interaction 
	movapd xmm0, [esp + nb334nf_rinvMM]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb334nf_rsqMM] ;# xmm1=r 
	mulsd  xmm1, [esp + nb334nf_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb334nf_VFtab]
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb334nf_qqMM]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	
    	addsd  xmm5, [esp + nb334nf_vctot]
    	movsd [esp + nb334nf_vctot], xmm5
	
.nb334nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb334nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb334nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb334nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb334nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb334nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb334nf_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb334nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb334nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb334nf_n], esi
        jmp .nb334nf_outer
.nb334nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb334nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb334nf_end
        ;# non-zero, do one more workunit
        jmp   .nb334nf_threadloop
.nb334nf_end:
	emms

	mov eax, [esp + nb334nf_nouter]
	mov ebx, [esp + nb334nf_ninner]
	mov ecx, [ebp + nb334nf_outeriter]
	mov edx, [ebp + nb334nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb334nf_salign]
	add esp, eax
	add esp, 904
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret
