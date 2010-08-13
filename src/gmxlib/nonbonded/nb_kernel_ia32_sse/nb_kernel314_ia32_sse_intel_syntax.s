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


	
.globl nb_kernel314_ia32_sse
.globl _nb_kernel314_ia32_sse
nb_kernel314_ia32_sse:	
_nb_kernel314_ia32_sse:	
.equiv          nb314_p_nri,            8
.equiv          nb314_iinr,             12
.equiv          nb314_jindex,           16
.equiv          nb314_jjnr,             20
.equiv          nb314_shift,            24
.equiv          nb314_shiftvec,         28
.equiv          nb314_fshift,           32
.equiv          nb314_gid,              36
.equiv          nb314_pos,              40
.equiv          nb314_faction,          44
.equiv          nb314_charge,           48
.equiv          nb314_p_facel,          52
.equiv          nb314_argkrf,           56
.equiv          nb314_argcrf,           60
.equiv          nb314_Vc,               64
.equiv          nb314_type,             68
.equiv          nb314_p_ntype,          72
.equiv          nb314_vdwparam,         76
.equiv          nb314_Vvdw,             80
.equiv          nb314_p_tabscale,       84
.equiv          nb314_VFtab,            88
.equiv          nb314_invsqrta,         92
.equiv          nb314_dvda,             96
.equiv          nb314_p_gbtabscale,     100
.equiv          nb314_GBtab,            104
.equiv          nb314_p_nthreads,       108
.equiv          nb314_count,            112
.equiv          nb314_mtx,              116
.equiv          nb314_outeriter,        120
.equiv          nb314_inneriter,        124
.equiv          nb314_work,             128
;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb314_ixO,              0
.equiv          nb314_iyO,              16
.equiv          nb314_izO,              32
.equiv          nb314_ixH1,             48
.equiv          nb314_iyH1,             64
.equiv          nb314_izH1,             80
.equiv          nb314_ixH2,             96
.equiv          nb314_iyH2,             112
.equiv          nb314_izH2,             128
.equiv          nb314_ixM,              144
.equiv          nb314_iyM,              160
.equiv          nb314_izM,              176
.equiv          nb314_jxO,              192
.equiv          nb314_jyO,              208
.equiv          nb314_jzO,              224
.equiv          nb314_jxH1,             240
.equiv          nb314_jyH1,             256
.equiv          nb314_jzH1,             272
.equiv          nb314_jxH2,             288
.equiv          nb314_jyH2,             304
.equiv          nb314_jzH2,             320
.equiv          nb314_jxM,              336
.equiv          nb314_jyM,              352
.equiv          nb314_jzM,              368
.equiv          nb314_dxOO,             384
.equiv          nb314_dyOO,             400
.equiv          nb314_dzOO,             416
.equiv          nb314_dxH1H1,           432
.equiv          nb314_dyH1H1,           448
.equiv          nb314_dzH1H1,           464
.equiv          nb314_dxH1H2,           480
.equiv          nb314_dyH1H2,           496
.equiv          nb314_dzH1H2,           512
.equiv          nb314_dxH1M,            528
.equiv          nb314_dyH1M,            544
.equiv          nb314_dzH1M,            560
.equiv          nb314_dxH2H1,           576
.equiv          nb314_dyH2H1,           592
.equiv          nb314_dzH2H1,           608
.equiv          nb314_dxH2H2,           624
.equiv          nb314_dyH2H2,           640
.equiv          nb314_dzH2H2,           656
.equiv          nb314_dxH2M,            672
.equiv          nb314_dyH2M,            688
.equiv          nb314_dzH2M,            704
.equiv          nb314_dxMH1,            720
.equiv          nb314_dyMH1,            736
.equiv          nb314_dzMH1,            752
.equiv          nb314_dxMH2,            768
.equiv          nb314_dyMH2,            784
.equiv          nb314_dzMH2,            800
.equiv          nb314_dxMM,             816
.equiv          nb314_dyMM,             832
.equiv          nb314_dzMM,             848
.equiv          nb314_qqMM,             864
.equiv          nb314_qqMH,             880
.equiv          nb314_qqHH,             896
.equiv          nb314_two,              912
.equiv          nb314_tsc,              928
.equiv          nb314_c6,               944
.equiv          nb314_c12,              960
.equiv          nb314_six,              976
.equiv          nb314_twelve,           992
.equiv          nb314_vctot,            1008
.equiv          nb314_Vvdwtot,          1024
.equiv          nb314_fixO,             1040
.equiv          nb314_fiyO,             1056
.equiv          nb314_fizO,             1072
.equiv          nb314_fixH1,            1088
.equiv          nb314_fiyH1,            1104
.equiv          nb314_fizH1,            1120
.equiv          nb314_fixH2,            1136
.equiv          nb314_fiyH2,            1152
.equiv          nb314_fizH2,            1168
.equiv          nb314_fixM,             1184
.equiv          nb314_fiyM,             1200
.equiv          nb314_fizM,             1216
.equiv          nb314_fjxO,             1232
.equiv          nb314_fjyO,             1248
.equiv          nb314_fjzO,             1264
.equiv          nb314_fjxH1,            1280
.equiv          nb314_fjyH1,            1296
.equiv          nb314_fjzH1,            1312
.equiv          nb314_fjxH2,            1328
.equiv          nb314_fjyH2,            1344
.equiv          nb314_fjzH2,            1360
.equiv          nb314_fjxM,             1376
.equiv          nb314_fjyM,             1392
.equiv          nb314_fjzM,             1408
.equiv          nb314_half,             1424
.equiv          nb314_three,            1440
.equiv          nb314_rsqOO,            1456
.equiv          nb314_rsqH1H1,          1472
.equiv          nb314_rsqH1H2,          1488
.equiv          nb314_rsqH1M,           1504
.equiv          nb314_rsqH2H1,          1520
.equiv          nb314_rsqH2H2,          1536
.equiv          nb314_rsqH2M,           1552
.equiv          nb314_rsqMH1,           1568
.equiv          nb314_rsqMH2,           1584
.equiv          nb314_rsqMM,            1600
.equiv          nb314_rinvsqOO,         1616
.equiv          nb314_rinvH1H1,         1632
.equiv          nb314_rinvH1H2,         1648
.equiv          nb314_rinvH1M,          1664
.equiv          nb314_rinvH2H1,         1680
.equiv          nb314_rinvH2H2,         1696
.equiv          nb314_rinvH2M,          1712
.equiv          nb314_rinvMH1,          1728
.equiv          nb314_rinvMH2,          1744
.equiv          nb314_rinvMM,           1760
.equiv          nb314_fstmp,            1776
.equiv          nb314_is3,              1792
.equiv          nb314_ii3,              1796
.equiv          nb314_innerjjnr,        1800
.equiv          nb314_innerk,           1804
.equiv          nb314_n,                1808
.equiv          nb314_nn1,              1812
.equiv          nb314_nri,              1816
.equiv          nb314_nouter,           1820
.equiv          nb314_ninner,           1824
.equiv          nb314_salign,           1828
	push ebp
	mov ebp,esp	
	push eax
	push ebx
	push ecx
	push edx
	push esi
	push edi
	sub esp, 1832		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb314_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb314_p_nri]
	mov ecx, [ecx]
	mov [esp + nb314_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb314_nouter], eax
	mov [esp + nb314_ninner], eax


	mov eax, [ebp + nb314_p_tabscale]
	movss xmm5, [eax]
	shufps xmm5, xmm5, 0
	movaps [esp + nb314_tsc],  xmm5
	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb314_half], eax
	movss xmm1, [esp + nb314_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps xmm4, xmm3
	addps  xmm4, xmm4	;# 6.0
	movaps xmm5, xmm4
	addps  xmm5, xmm5	;# constant 12.0
	movaps [esp + nb314_half],  xmm1
	movaps [esp + nb314_two],  xmm2
	movaps [esp + nb314_three],  xmm3
	movaps [esp + nb314_six],  xmm4
	movaps [esp + nb314_twelve],  xmm5

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb314_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]		;# ebx =ii 

	mov   edx, [ebp + nb314_charge]
	movss xmm5, [edx + ebx*4 + 4]	
	movss xmm3, [edx + ebx*4 + 12]	
	movss xmm4, xmm3	
	mov esi, [ebp + nb314_p_facel]
	movss xmm6, [esi]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb314_qqMM], xmm3
	movaps [esp + nb314_qqMH], xmm4
	movaps [esp + nb314_qqHH], xmm5
		
	xorps xmm0, xmm0
	mov   edx, [ebp + nb314_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov edi, [ebp + nb314_p_ntype]
	imul  ecx, [edi]  	;# ecx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   eax, [ebp + nb314_vdwparam]
	movlps xmm0, [eax + edx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0x55
	movaps [esp + nb314_c6], xmm0
	movaps [esp + nb314_c12], xmm1

.nb314_threadloop:
        mov   esi, [ebp + nb314_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb314_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb314_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb314_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb314_n], eax
        mov [esp + nb314_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb314_outerstart
        jmp .nb314_end
	
.nb314_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb314_nouter]
	mov [esp + nb314_nouter], ebx

.nb314_outer:
	mov   eax, [ebp + nb314_shift]  	;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]	;# ebx=3*is 
	mov   [esp + nb314_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb314_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb314_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]		;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb314_pos]	;# eax = base of pos[]  
	mov   [esp + nb314_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	movaps xmm6, xmm0
	movaps xmm7, xmm1

	addss xmm3, [eax + ebx*4]  	;# ox
	addss xmm4, [eax + ebx*4 + 4]  ;# oy
	addss xmm5, [eax + ebx*4 + 8]  ;# oz
	addss xmm6, [eax + ebx*4 + 12] ;# h1x
	addss xmm7, [eax + ebx*4 + 16] ;# h1y
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	shufps xmm7, xmm7, 0
	movaps [esp + nb314_ixO], xmm3
	movaps [esp + nb314_iyO], xmm4
	movaps [esp + nb314_izO], xmm5
	movaps [esp + nb314_ixH1], xmm6
	movaps [esp + nb314_iyH1], xmm7

	movss xmm6, xmm2
	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm6, [eax + ebx*4 + 20] ;# h1z
	addss xmm0, [eax + ebx*4 + 24] ;# h2x
	addss xmm1, [eax + ebx*4 + 28] ;# h2y
	addss xmm2, [eax + ebx*4 + 32] ;# h2z
	addss xmm3, [eax + ebx*4 + 36] ;# mx
	addss xmm4, [eax + ebx*4 + 40] ;# my
	addss xmm5, [eax + ebx*4 + 44] ;# mz

	shufps xmm6, xmm6, 0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb314_izH1], xmm6
	movaps [esp + nb314_ixH2], xmm0
	movaps [esp + nb314_iyH2], xmm1
	movaps [esp + nb314_izH2], xmm2
	movaps [esp + nb314_ixM], xmm3
	movaps [esp + nb314_iyM], xmm4
	movaps [esp + nb314_izM], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb314_vctot], xmm4
	movaps [esp + nb314_Vvdwtot], xmm4
	movaps [esp + nb314_fixO], xmm4
	movaps [esp + nb314_fiyO], xmm4
	movaps [esp + nb314_fizO], xmm4
	movaps [esp + nb314_fixH1], xmm4
	movaps [esp + nb314_fiyH1], xmm4
	movaps [esp + nb314_fizH1], xmm4
	movaps [esp + nb314_fixH2], xmm4
	movaps [esp + nb314_fiyH2], xmm4
	movaps [esp + nb314_fizH2], xmm4
	movaps [esp + nb314_fixM], xmm4
	movaps [esp + nb314_fiyM], xmm4
	movaps [esp + nb314_fizM], xmm4
	
	mov   eax, [ebp + nb314_jindex]
	mov   ecx, [eax + esi*4]	 	;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   esi, [ebp + nb314_pos]
	mov   edi, [ebp + nb314_faction]	
	mov   eax, [ebp + nb314_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb314_innerjjnr], eax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb314_ninner]
	mov   [esp + nb314_ninner], ecx
	add   edx, 0
	mov   [esp + nb314_innerk], edx	;# number of innerloop atoms 
	jge   .nb314_unroll_loop
	jmp   .nb314_single_check
.nb314_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb314_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	mov   ecx, [edx + 8]
	mov   edx, [edx + 12]     	;# eax-edx=jnr1-4 
	
	add dword ptr [esp + nb314_innerjjnr],  16 ;# advance pointer (unroll 4) 

	mov esi, [ebp + nb314_pos]   	;# base of pos[] 

	lea   eax, [eax + eax*2] 	;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2] 	;# replace jnr with j3 
	lea   edx, [edx + edx*2]	
	
	;# move j coordinates to local temp variables
	;# Load Ox, Oy, Oz, H1x 
	movlps xmm1, [esi + eax*4]	;#  Oxa   Oya    -    -
	movlps xmm4, [esi + ecx*4]	;#  Oxc   Oyc    -    -
	movhps xmm1, [esi + ebx*4]	;#  Oxa   Oya   Oxb   Oyb 
	movhps xmm4, [esi + edx*4]	;#  Oxc   Oyc   Oxd   Oyd 
	movaps xmm0, xmm1		;#  Oxa   Oya   Oxb   Oyb 
	shufps xmm0, xmm4, 0x88		;#  Oxa   Oxb   Oxc   Oxd
	shufps xmm1, xmm4, 0xDD		;#  Oya   Oyb   Oyc   Oyd
	movlps xmm3, [esi + eax*4 + 8]	;#  Oza  H1xa    -    -
	movlps xmm5, [esi + ecx*4 + 8]	;#  Ozc  H1xc    -    -
	movhps xmm3, [esi + ebx*4 + 8]	;#  Oza  H1xa   Ozb  H1xb 
	movhps xmm5, [esi + edx*4 + 8]	;#  Ozc  H1xc   Ozd  H1xd 
	movaps xmm2, xmm3		;#  Oza  H1xa   Ozb  H1xb 
	shufps xmm2, xmm5, 0x88		;#  Oza   Ozb   Ozc   Ozd
	shufps xmm3, xmm5, 0xDD		;# H1xa  H1xb  H1xc  H1xd
	;# coordinates in xmm0-xmm3	
	;# store
	movaps [esp + nb314_jxO], xmm0
	movaps [esp + nb314_jyO], xmm1
	movaps [esp + nb314_jzO], xmm2
	movaps [esp + nb314_jxH1], xmm3

	;# Load H1y H1z H2x H2y 
	movlps xmm1, [esi + eax*4 + 16]	
	movlps xmm4, [esi + ecx*4 + 16]	
	movhps xmm1, [esi + ebx*4 + 16]	
	movhps xmm4, [esi + edx*4 + 16]	
	movaps xmm0, xmm1		
	shufps xmm0, xmm4, 0x88		
	shufps xmm1, xmm4, 0xDD		
	movlps xmm3, [esi + eax*4 + 24]	
	movlps xmm5, [esi + ecx*4 + 24]	
	movhps xmm3, [esi + ebx*4 + 24]	
	movhps xmm5, [esi + edx*4 + 24]	
	movaps xmm2, xmm3		
	shufps xmm2, xmm5, 0x88		
	shufps xmm3, xmm5, 0xDD		
	;# coordinates in xmm0-xmm3	
	;# store
	movaps [esp + nb314_jyH1], xmm0
	movaps [esp + nb314_jzH1], xmm1
	movaps [esp + nb314_jxH2], xmm2
	movaps [esp + nb314_jyH2], xmm3

	;# Load H2z Mx My Mz 
	movlps xmm1, [esi + eax*4 + 32]	
	movlps xmm4, [esi + ecx*4 + 32]	
	movhps xmm1, [esi + ebx*4 + 32]	
	movhps xmm4, [esi + edx*4 + 32]	
	movaps xmm0, xmm1		
	shufps xmm0, xmm4, 0x88		
	shufps xmm1, xmm4, 0xDD		
	movlps xmm3, [esi + eax*4 + 40]	
	movlps xmm5, [esi + ecx*4 + 40]	
	movhps xmm3, [esi + ebx*4 + 40]	
	movhps xmm5, [esi + edx*4 + 40]	
	movaps xmm2, xmm3		
	shufps xmm2, xmm5, 0x88		
	shufps xmm3, xmm5, 0xDD		
	;# coordinates in xmm0-xmm3	
	;# store
	movaps [esp + nb314_jzH2], xmm0
	movaps [esp + nb314_jxM], xmm1
	movaps [esp + nb314_jyM], xmm2
	movaps [esp + nb314_jzM], xmm3
	
	;# start calculating pairwise distances
	movaps xmm0, [esp + nb314_ixO]
	movaps xmm1, [esp + nb314_iyO]
	movaps xmm2, [esp + nb314_izO]
	movaps xmm3, [esp + nb314_ixH1]
	movaps xmm4, [esp + nb314_iyH1]
	movaps xmm5, [esp + nb314_izH1]
	subps  xmm0, [esp + nb314_jxO]
	subps  xmm1, [esp + nb314_jyO]
	subps  xmm2, [esp + nb314_jzO]
	subps  xmm3, [esp + nb314_jxH1]
	subps  xmm4, [esp + nb314_jyH1]
	subps  xmm5, [esp + nb314_jzH1]
	movaps [esp + nb314_dxOO], xmm0
	movaps [esp + nb314_dyOO], xmm1
	movaps [esp + nb314_dzOO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb314_dxH1H1], xmm3
	movaps [esp + nb314_dyH1H1], xmm4
	movaps [esp + nb314_dzH1H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb314_rsqOO], xmm0
	movaps [esp + nb314_rsqH1H1], xmm3

	movaps xmm0, [esp + nb314_ixH1]
	movaps xmm1, [esp + nb314_iyH1]
	movaps xmm2, [esp + nb314_izH1]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [esp + nb314_jxH2]
	subps  xmm1, [esp + nb314_jyH2]
	subps  xmm2, [esp + nb314_jzH2]
	subps  xmm3, [esp + nb314_jxM]
	subps  xmm4, [esp + nb314_jyM]
	subps  xmm5, [esp + nb314_jzM]
	movaps [esp + nb314_dxH1H2], xmm0
	movaps [esp + nb314_dyH1H2], xmm1
	movaps [esp + nb314_dzH1H2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb314_dxH1M], xmm3
	movaps [esp + nb314_dyH1M], xmm4
	movaps [esp + nb314_dzH1M], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb314_rsqH1H2], xmm0
	movaps [esp + nb314_rsqH1M], xmm3

	movaps xmm0, [esp + nb314_ixH2]
	movaps xmm1, [esp + nb314_iyH2]
	movaps xmm2, [esp + nb314_izH2]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [esp + nb314_jxH1]
	subps  xmm1, [esp + nb314_jyH1]
	subps  xmm2, [esp + nb314_jzH1]
	subps  xmm3, [esp + nb314_jxH2]
	subps  xmm4, [esp + nb314_jyH2]
	subps  xmm5, [esp + nb314_jzH2]
	movaps [esp + nb314_dxH2H1], xmm0
	movaps [esp + nb314_dyH2H1], xmm1
	movaps [esp + nb314_dzH2H1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb314_dxH2H2], xmm3
	movaps [esp + nb314_dyH2H2], xmm4
	movaps [esp + nb314_dzH2H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb314_rsqH2H1], xmm0
	movaps [esp + nb314_rsqH2H2], xmm3

	movaps xmm0, [esp + nb314_ixH2]
	movaps xmm1, [esp + nb314_iyH2]
	movaps xmm2, [esp + nb314_izH2]
	movaps xmm3, [esp + nb314_ixM]
	movaps xmm4, [esp + nb314_iyM]
	movaps xmm5, [esp + nb314_izM]
	subps  xmm0, [esp + nb314_jxM]
	subps  xmm1, [esp + nb314_jyM]
	subps  xmm2, [esp + nb314_jzM]
	subps  xmm3, [esp + nb314_jxH1]
	subps  xmm4, [esp + nb314_jyH1]
	subps  xmm5, [esp + nb314_jzH1]
	movaps [esp + nb314_dxH2M], xmm0
	movaps [esp + nb314_dyH2M], xmm1
	movaps [esp + nb314_dzH2M], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb314_dxMH1], xmm3
	movaps [esp + nb314_dyMH1], xmm4
	movaps [esp + nb314_dzMH1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [esp + nb314_rsqH2M], xmm0
	movaps [esp + nb314_rsqMH1], xmm4

	movaps xmm0, [esp + nb314_ixM]
	movaps xmm1, [esp + nb314_iyM]
	movaps xmm2, [esp + nb314_izM]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [esp + nb314_jxH2]
	subps  xmm1, [esp + nb314_jyH2]
	subps  xmm2, [esp + nb314_jzH2]
	subps  xmm3, [esp + nb314_jxM]
	subps  xmm4, [esp + nb314_jyM]
	subps  xmm5, [esp + nb314_jzM]
	movaps [esp + nb314_dxMH2], xmm0
	movaps [esp + nb314_dyMH2], xmm1
	movaps [esp + nb314_dzMH2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb314_dxMM], xmm3
	movaps [esp + nb314_dyMM], xmm4
	movaps [esp + nb314_dzMM], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [esp + nb314_rsqMH2], xmm0
	movaps [esp + nb314_rsqMM], xmm4

	;# start by doing reciprocal for OO
	movaps  xmm7, [esp + nb314_rsqOO]
	rcpps   xmm2, xmm7
	movaps  xmm1, [esp + nb314_two]
	mulps   xmm7, xmm2
	subps   xmm1, xmm7
	mulps   xmm2, xmm1 ;# rinvsq 
	movaps [esp + nb314_rinvsqOO], xmm2
	
	;# next step is invsqrt - do two at a time.
	rsqrtps xmm1, [esp + nb314_rsqH1H1]
	rsqrtps xmm5, [esp + nb314_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb314_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb314_rsqH1H1]
	mulps   xmm5, [esp + nb314_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb314_half] ;# rinvH1H1 
	mulps   xmm7, [esp + nb314_half] ;# rinvH1H2 
	movaps  [esp + nb314_rinvH1H1], xmm3
	movaps  [esp + nb314_rinvH1H2], xmm7
			
	rsqrtps xmm1, [esp + nb314_rsqH1M]
	rsqrtps xmm5, [esp + nb314_rsqH2H1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb314_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb314_rsqH1M]
	mulps   xmm5, [esp + nb314_rsqH2H1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb314_half] 
	mulps   xmm7, [esp + nb314_half]
	movaps  [esp + nb314_rinvH1M], xmm3
	movaps  [esp + nb314_rinvH2H1], xmm7
			
	rsqrtps xmm1, [esp + nb314_rsqH2H2]
	rsqrtps xmm5, [esp + nb314_rsqH2M]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb314_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb314_rsqH2H2]
	mulps   xmm5, [esp + nb314_rsqH2M]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb314_half] 
	mulps   xmm7, [esp + nb314_half]
	movaps  [esp + nb314_rinvH2H2], xmm3
	movaps  [esp + nb314_rinvH2M], xmm7
	
	rsqrtps xmm1, [esp + nb314_rsqMH1]
	rsqrtps xmm5, [esp + nb314_rsqMH2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb314_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb314_rsqMH1]
	mulps   xmm5, [esp + nb314_rsqMH2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb314_half] 
	mulps   xmm7, [esp + nb314_half]
	movaps  [esp + nb314_rinvMH1], xmm3
	movaps  [esp + nb314_rinvMH2], xmm7
        		
	rsqrtps xmm1, [esp + nb314_rsqMM]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + nb314_three]
	mulps   xmm1, [esp + nb314_rsqMM]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [esp + nb314_half] 
	movaps  [esp + nb314_rinvMM], xmm3
	
	;# start with OO LJ interaction
	movaps xmm0, [esp + nb314_rinvsqOO]
	movaps xmm1, xmm0
	mulps  xmm1, xmm1	;# rinv4
	mulps  xmm1, xmm0	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm1, [esp + nb314_c6]
	mulps  xmm2, [esp + nb314_c12]
	movaps xmm4, xmm2
	subps  xmm4, xmm1
	addps  xmm4, [esp + nb314_Vvdwtot]
	mulps  xmm1, [esp + nb314_six]
	mulps  xmm2, [esp + nb314_twelve]
	movaps [esp + nb314_Vvdwtot], xmm4
	subps  xmm2, xmm1
	mulps  xmm0, xmm2 	;# fscal 
	movaps xmm1, xmm0
	movaps xmm2, xmm0		

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + nb314_dxOO]
	mulps xmm1, [esp + nb314_dyOO]
	mulps xmm2, [esp + nb314_dzOO]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb314_fixO]
	addps xmm1, [esp + nb314_fiyO]
	addps xmm2, [esp + nb314_fizO]
	movaps [esp + nb314_fjxO], xmm3
	movaps [esp + nb314_fjyO], xmm4
	movaps [esp + nb314_fjzO], xmm5
	movaps [esp + nb314_fixO], xmm0
	movaps [esp + nb314_fiyO], xmm1
	movaps [esp + nb314_fizO], xmm2

	;# Coulomb interactions - first H1H1
	movaps xmm0, [esp + nb314_rinvH1H1]
	
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb314_rsqH1H1] ;# xmm1=r 
	mulps  xmm1, [esp + nb314_tsc]
		
	movhlps xmm2, xmm1
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 
	pslld   mm6, 2
	pslld   mm7, 2
	
	movd mm0, eax
	movd mm1, ebx
	movd mm2, ecx
	movd mm3, edx

	mov  esi, [ebp + nb314_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [esp + nb314_two]   	;# two*Heps2 
	movaps xmm3, [esp + nb314_qqHH]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 
	;# update vctot 
        addps  xmm5, [esp + nb314_vctot]
        movaps [esp + nb314_vctot], xmm5
		
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + nb314_tsc]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	
	mulps xmm0, [esp + nb314_dxH1H1]
	mulps xmm1, [esp + nb314_dyH1H1]
	mulps xmm2, [esp + nb314_dzH1H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb314_fixH1]
	addps xmm1, [esp + nb314_fiyH1]
	addps xmm2, [esp + nb314_fizH1]
	movaps [esp + nb314_fjxH1], xmm3
	movaps [esp + nb314_fjyH1], xmm4
	movaps [esp + nb314_fjzH1], xmm5
	movaps [esp + nb314_fixH1], xmm0
	movaps [esp + nb314_fiyH1], xmm1
	movaps [esp + nb314_fizH1], xmm2

	;# H1-H2 interaction 
	movaps xmm0, [esp + nb314_rinvH1H2]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb314_rsqH1H2] ;# xmm1=r 
	mulps  xmm1, [esp + nb314_tsc]	
	movhlps xmm2, xmm1
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2
	
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [esp + nb314_two]   	;# two*Heps2 
	movaps xmm3, [esp + nb314_qqHH]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 

	addps  xmm5, [esp + nb314_vctot]
	movaps [esp + nb314_vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + nb314_tsc]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + nb314_dxH1H2]
	mulps xmm1, [esp + nb314_dyH1H2]
	mulps xmm2, [esp + nb314_dzH1H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb314_fixH1]
	addps xmm1, [esp + nb314_fiyH1]
	addps xmm2, [esp + nb314_fizH1]
	movaps [esp + nb314_fjxH2], xmm3
	movaps [esp + nb314_fjyH2], xmm4
	movaps [esp + nb314_fjzH2], xmm5
	movaps [esp + nb314_fixH1], xmm0
	movaps [esp + nb314_fiyH1], xmm1
	movaps [esp + nb314_fizH1], xmm2

	;# H1-M interaction  
	movaps xmm0, [esp + nb314_rinvH1M]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb314_rsqH1M] ;# xmm1=r 
	mulps  xmm1, [esp + nb314_tsc]	
	movhlps xmm2, xmm1
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2

	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [esp + nb314_two]   	;# two*Heps2 
	movaps xmm3, [esp + nb314_qqMH]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 

	addps  xmm5, [esp + nb314_vctot]
	movaps [esp + nb314_vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + nb314_tsc]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + nb314_dxH1M]
	mulps xmm1, [esp + nb314_dyH1M]
	mulps xmm2, [esp + nb314_dzH1M]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb314_fixH1]
	addps xmm1, [esp + nb314_fiyH1]
	addps xmm2, [esp + nb314_fizH1]
	movaps [esp + nb314_fjxM], xmm3
	movaps [esp + nb314_fjyM], xmm4
	movaps [esp + nb314_fjzM], xmm5
	movaps [esp + nb314_fixH1], xmm0
	movaps [esp + nb314_fiyH1], xmm1
	movaps [esp + nb314_fizH1], xmm2

	;# H2-H1 interaction 
	movaps xmm0, [esp + nb314_rinvH2H1]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb314_rsqH2H1] ;# xmm1=r 
	mulps  xmm1, [esp + nb314_tsc]	
	movhlps xmm2, xmm1
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2

	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [esp + nb314_two]   	;# two*Heps2 
	movaps xmm3, [esp + nb314_qqHH]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 

	addps  xmm5, [esp + nb314_vctot]
	movaps [esp + nb314_vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + nb314_tsc]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [esp + nb314_fjxH1]
	movaps xmm4, [esp + nb314_fjyH1]
	movaps xmm5, [esp + nb314_fjzH1]
	mulps xmm0, [esp + nb314_dxH2H1]
	mulps xmm1, [esp + nb314_dyH2H1]
	mulps xmm2, [esp + nb314_dzH2H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb314_fixH2]
	addps xmm1, [esp + nb314_fiyH2]
	addps xmm2, [esp + nb314_fizH2]
	movaps [esp + nb314_fjxH1], xmm3
	movaps [esp + nb314_fjyH1], xmm4
	movaps [esp + nb314_fjzH1], xmm5
	movaps [esp + nb314_fixH2], xmm0
	movaps [esp + nb314_fiyH2], xmm1
	movaps [esp + nb314_fizH2], xmm2

	;# H2-H2 interaction 
	movaps xmm0, [esp + nb314_rinvH2H2]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb314_rsqH2H2] ;# xmm1=r 
	mulps  xmm1, [esp + nb314_tsc]	
	movhlps xmm2, xmm1
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2

	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [esp + nb314_two]   	;# two*Heps2 
	movaps xmm3, [esp + nb314_qqHH]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 

	addps  xmm5, [esp + nb314_vctot]
	movaps [esp + nb314_vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + nb314_tsc]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [esp + nb314_fjxH2]
	movaps xmm4, [esp + nb314_fjyH2]
	movaps xmm5, [esp + nb314_fjzH2]
	mulps xmm0, [esp + nb314_dxH2H2]
	mulps xmm1, [esp + nb314_dyH2H2]
	mulps xmm2, [esp + nb314_dzH2H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb314_fixH2]
	addps xmm1, [esp + nb314_fiyH2]
	addps xmm2, [esp + nb314_fizH2]
	movaps [esp + nb314_fjxH2], xmm3
	movaps [esp + nb314_fjyH2], xmm4
	movaps [esp + nb314_fjzH2], xmm5
	movaps [esp + nb314_fixH2], xmm0
	movaps [esp + nb314_fiyH2], xmm1
	movaps [esp + nb314_fizH2], xmm2

	;# H2-M interaction 
	movaps xmm0, [esp + nb314_rinvH2M]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb314_rsqH2M] ;# xmm1=r 
	mulps  xmm1, [esp + nb314_tsc]
	movhlps xmm2, xmm1
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2

	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [esp + nb314_two]   	;# two*Heps2 
	movaps xmm3, [esp + nb314_qqMH]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 

	addps  xmm5, [esp + nb314_vctot]
	movaps [esp + nb314_vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + nb314_tsc]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [esp + nb314_fjxM]
	movaps xmm4, [esp + nb314_fjyM]
	movaps xmm5, [esp + nb314_fjzM]
	mulps xmm0, [esp + nb314_dxH2M]
	mulps xmm1, [esp + nb314_dyH2M]
	mulps xmm2, [esp + nb314_dzH2M]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb314_fixH2]
	addps xmm1, [esp + nb314_fiyH2]
	addps xmm2, [esp + nb314_fizH2]
	movaps [esp + nb314_fjxM], xmm3
	movaps [esp + nb314_fjyM], xmm4
	movaps [esp + nb314_fjzM], xmm5
	movaps [esp + nb314_fixH2], xmm0
	movaps [esp + nb314_fiyH2], xmm1
	movaps [esp + nb314_fizH2], xmm2

	;# M-H1 interaction 
	movaps xmm0, [esp + nb314_rinvMH1]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb314_rsqMH1] ;# xmm1=r 
	mulps  xmm1, [esp + nb314_tsc]	
	movhlps xmm2, xmm1
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2

	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [esp + nb314_two]   	;# two*Heps2 
	movaps xmm3, [esp + nb314_qqMH]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 

	addps  xmm5, [esp + nb314_vctot]
	movaps [esp + nb314_vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + nb314_tsc]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1

	movaps xmm3, [esp + nb314_fjxH1]
	movaps xmm4, [esp + nb314_fjyH1]
	movaps xmm5, [esp + nb314_fjzH1]
	mulps xmm0, [esp + nb314_dxMH1]
	mulps xmm1, [esp + nb314_dyMH1]
	mulps xmm2, [esp + nb314_dzMH1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb314_fixM]
	addps xmm1, [esp + nb314_fiyM]
	addps xmm2, [esp + nb314_fizM]
	movaps [esp + nb314_fjxH1], xmm3
	movaps [esp + nb314_fjyH1], xmm4
	movaps [esp + nb314_fjzH1], xmm5
	movaps [esp + nb314_fixM], xmm0
	movaps [esp + nb314_fiyM], xmm1
	movaps [esp + nb314_fizM], xmm2

	;# M-H2 interaction 
	movaps xmm0, [esp + nb314_rinvMH2]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb314_rsqMH2] ;# xmm1=r 
	mulps  xmm1, [esp + nb314_tsc]
	movhlps xmm2, xmm1
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2

	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [esp + nb314_two]   	;# two*Heps2 
	movaps xmm3, [esp + nb314_qqMH]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 

	addps  xmm5, [esp + nb314_vctot]
	movaps [esp + nb314_vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + nb314_tsc]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [esp + nb314_fjxH2]
	movaps xmm4, [esp + nb314_fjyH2]
	movaps xmm5, [esp + nb314_fjzH2]
	mulps xmm0, [esp + nb314_dxMH2]
	mulps xmm1, [esp + nb314_dyMH2]
	mulps xmm2, [esp + nb314_dzMH2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb314_fixM]
	addps xmm1, [esp + nb314_fiyM]
	addps xmm2, [esp + nb314_fizM]
	movaps [esp + nb314_fjxH2], xmm3
	movaps [esp + nb314_fjyH2], xmm4
	movaps [esp + nb314_fjzH2], xmm5
	movaps [esp + nb314_fixM], xmm0
	movaps [esp + nb314_fiyM], xmm1
	movaps [esp + nb314_fizM], xmm2

	;# M-M interaction 
	movaps xmm0, [esp + nb314_rinvMM]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb314_rsqMM] ;# xmm1=r 
	mulps  xmm1, [esp + nb314_tsc]	
	movhlps xmm2, xmm1
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2

	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [esp + nb314_two]   	;# two*Heps2 
	movaps xmm3, [esp + nb314_qqMM]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 

	addps  xmm5, [esp + nb314_vctot]
	movaps [esp + nb314_vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + nb314_tsc]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [esp + nb314_fjxM]
	movaps xmm4, [esp + nb314_fjyM]
	movaps xmm5, [esp + nb314_fjzM]
	mulps xmm0, [esp + nb314_dxMM]
	mulps xmm1, [esp + nb314_dyMM]
	mulps xmm2, [esp + nb314_dzMM]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb314_fixM]
	addps xmm1, [esp + nb314_fiyM]
	addps xmm2, [esp + nb314_fizM]
	movaps [esp + nb314_fjxM], xmm3
	movaps [esp + nb314_fjyM], xmm4
	movaps [esp + nb314_fjzM], xmm5
	movaps [esp + nb314_fixM], xmm0
	movaps [esp + nb314_fiyM], xmm1
	movaps [esp + nb314_fizM], xmm2

	mov edi, [ebp + nb314_faction]

	movd eax, mm0
	movd ebx, mm1
	movd ecx, mm2
	movd edx, mm3
	
	;# Did all interactions - now update j forces 
	;# 4 j waters with four atoms each.
	;# step 1 : transpose fjxO, fjyO, fjzO, fjxH1
	movaps xmm0, [esp + nb314_fjxO]
	movaps xmm1, [esp + nb314_fjyO]
	movaps xmm2, [esp + nb314_fjzO]
	movaps xmm3, [esp + nb314_fjxH1]
	movaps xmm4, xmm0
	movaps xmm5, xmm1
	unpcklps xmm4, xmm2 	
	unpcklps xmm5, xmm3	
	unpckhps xmm0, xmm2 
	unpckhps xmm1, xmm3 
	movaps xmm2, xmm4
	movaps xmm3, xmm0
	unpcklps xmm4, xmm5 
	unpckhps xmm2, xmm5 

	unpcklps xmm0, xmm1 
	unpckhps xmm3, xmm1 
	;# results are now in xmm4, xmm2, xmm0, xmm3
	;# load the corresponding j forces from memory
	movlps   xmm1, [edi + eax*4]
	movlps   xmm5, [edi + ebx*4]
	movlps   xmm6, [edi + ecx*4]
	movlps   xmm7, [edi + edx*4]
	movhps   xmm1, [edi + eax*4 + 8]
	movhps   xmm5, [edi + ebx*4 + 8]
	movhps   xmm6, [edi + ecx*4 + 8]
	movhps   xmm7, [edi + edx*4 + 8]
	;# add
	addps    xmm1, xmm4
	addps    xmm5, xmm2
	addps    xmm6, xmm0
	addps    xmm7, xmm3
	;# store back
	movlps   [edi + eax*4], xmm1
	movlps   [edi + ebx*4], xmm5
	movlps   [edi + ecx*4], xmm6
	movlps   [edi + edx*4], xmm7
	movhps   [edi + eax*4 + 8], xmm1
	movhps   [edi + ebx*4 + 8], xmm5
	movhps   [edi + ecx*4 + 8], xmm6
	movhps   [edi + edx*4 + 8], xmm7

	;# step 2 : transpose fjyH1, fjzH1, fjxH2, fjyH2
	movaps xmm0, [esp + nb314_fjyH1]
	movaps xmm1, [esp + nb314_fjzH1]
	movaps xmm2, [esp + nb314_fjxH2]
	movaps xmm3, [esp + nb314_fjyH2]
	movaps xmm4, xmm0
	movaps xmm5, xmm1
	unpcklps xmm4, xmm2 	
	unpcklps xmm5, xmm3	
	unpckhps xmm0, xmm2 
	unpckhps xmm1, xmm3 
	movaps xmm2, xmm4
	movaps xmm3, xmm0
	unpcklps xmm4, xmm5 
	unpckhps xmm2, xmm5 

	unpcklps xmm0, xmm1 
	unpckhps xmm3, xmm1 
	;# results are now in xmm4, xmm2, xmm0, xmm3
	;# load the corresponding j forces from memory
	movlps   xmm1, [edi + eax*4 + 16]
	movlps   xmm5, [edi + ebx*4 + 16]
	movlps   xmm6, [edi + ecx*4 + 16]
	movlps   xmm7, [edi + edx*4 + 16]
	movhps   xmm1, [edi + eax*4 + 24]
	movhps   xmm5, [edi + ebx*4 + 24]
	movhps   xmm6, [edi + ecx*4 + 24]
	movhps   xmm7, [edi + edx*4 + 24]
	;# add
	addps    xmm1, xmm4
	addps    xmm5, xmm2
	addps    xmm6, xmm0
	addps    xmm7, xmm3
	;# store back
	movlps   [edi + eax*4 + 16], xmm1
	movlps   [edi + ebx*4 + 16], xmm5
	movlps   [edi + ecx*4 + 16], xmm6
	movlps   [edi + edx*4 + 16], xmm7
	movhps   [edi + eax*4 + 24], xmm1
	movhps   [edi + ebx*4 + 24], xmm5
	movhps   [edi + ecx*4 + 24], xmm6
	movhps   [edi + edx*4 + 24], xmm7

	;# step 3 : transpose fjzH2, fjxM, fjyM, fjzM. xmm4 is scratch
	movaps xmm0, [esp + nb314_fjzH2]
	movaps xmm1, [esp + nb314_fjxM]
	movaps xmm2, [esp + nb314_fjyM]
	movaps xmm3, [esp + nb314_fjzM]
	
        movaps xmm4, xmm0
        movaps xmm5, xmm1
        unpcklps xmm4, xmm2
        unpcklps xmm5, xmm3
        unpckhps xmm0, xmm2
        unpckhps xmm1, xmm3
        movaps xmm2, xmm4
        movaps xmm3, xmm0
        unpcklps xmm4, xmm5
        unpckhps xmm2, xmm5
	
        unpcklps xmm0, xmm1
        unpckhps xmm3, xmm1
	
	;# results are now in xmm0, xmm1, xmm2, xmm3
	;# load the corresponding j forces from memory
	movlps   xmm1, [edi + eax*4 + 32]
	movlps   xmm5, [edi + ebx*4 + 32]
	movlps   xmm6, [edi + ecx*4 + 32]
	movlps   xmm7, [edi + edx*4 + 32]
	movhps   xmm1, [edi + eax*4 + 40]
	movhps   xmm5, [edi + ebx*4 + 40]
	movhps   xmm6, [edi + ecx*4 + 40]
	movhps   xmm7, [edi + edx*4 + 40]
	;# add
	addps    xmm1, xmm4
	addps    xmm5, xmm2
	addps    xmm6, xmm0
	addps    xmm7, xmm3
	;# store back
	movlps   [edi + eax*4 + 32], xmm1
	movlps   [edi + ebx*4 + 32], xmm5
	movlps   [edi + ecx*4 + 32], xmm6
	movlps   [edi + edx*4 + 32], xmm7
	movhps   [edi + eax*4 + 40], xmm1
	movhps   [edi + ebx*4 + 40], xmm5
	movhps   [edi + ecx*4 + 40], xmm6
	movhps   [edi + edx*4 + 40], xmm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb314_innerk],  4
	jl    .nb314_single_check
	jmp   .nb314_unroll_loop
.nb314_single_check:
	add dword ptr [esp + nb314_innerk],  4
	jnz   .nb314_single_loop
	jmp   .nb314_updateouterdata
.nb314_single_loop:
	mov   edx, [esp + nb314_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb314_innerjjnr],  4	

	mov esi, [ebp + nb314_pos]
	lea   eax, [eax + eax*2]  

	;# fetch j coordinates
	movlps xmm3,  [esi + eax*4]		;#  Ox  Oy  
	movlps xmm4,  [esi + eax*4 + 16]	;# H1y H1z 
	movlps xmm5,  [esi + eax*4 + 32]	;# H2z  Mx 
	movhps xmm3,  [esi + eax*4 + 8]   	;#  Ox  Oy  Oz H1x
	movhps xmm4,  [esi + eax*4 + 24]	;# H1y H1z H2x H2y
	movhps xmm5,  [esi + eax*4 + 40]	;# H2z  Mx  My  Mz
	;# transpose
	movaps xmm0, xmm4
	movaps xmm1, xmm3
	movaps xmm2, xmm4
	movaps xmm6, xmm3
	shufps xmm4, xmm5, 18  ;# (00010010)  h2x - Mx  - 
	shufps xmm3, xmm0, 193 ;# (11000001)  Oy  - H1y - 
	shufps xmm2, xmm5, 35  ;# (00100011) H2y - My  - 
 	shufps xmm1, xmm0, 18  ;# (00010010)  Oz  - H1z - 
	;#  xmm6: Ox - - H1x   xmm5: H2z - - Mz 
	shufps xmm6, xmm4, 140 ;# (10001100) Ox H1x H2x Mx 
	shufps xmm3, xmm2, 136 ;# (10001000) Oy H1y H2y My 
	shufps xmm1, xmm5, 200 ;# (11001000) Oz H1z H2z Mz

	;# store all j coordinates in jO  
	movaps [esp + nb314_jxO], xmm6
	movaps [esp + nb314_jyO], xmm3
	movaps [esp + nb314_jzO], xmm1

	;# do O and H1 in parallel
	movaps xmm0, [esp + nb314_ixO]
	movaps xmm1, [esp + nb314_iyO]
	movaps xmm2, [esp + nb314_izO]
	movaps xmm3, [esp + nb314_ixH1]
	movaps xmm4, [esp + nb314_iyH1]
	movaps xmm5, [esp + nb314_izH1]
	subps  xmm0, [esp + nb314_jxO]
	subps  xmm1, [esp + nb314_jyO]
	subps  xmm2, [esp + nb314_jzO]
	subps  xmm3, [esp + nb314_jxO]
	subps  xmm4, [esp + nb314_jyO]
	subps  xmm5, [esp + nb314_jzO]
	
	movaps [esp + nb314_dxOO], xmm0
	movaps [esp + nb314_dyOO], xmm1
	movaps [esp + nb314_dzOO], xmm2
	movaps [esp + nb314_dxH1H1], xmm3
	movaps [esp + nb314_dyH1H1], xmm4
	movaps [esp + nb314_dzH1H1], xmm5
	
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm4, xmm3
	addps xmm4, xmm5	;# have rsq in xmm4
	;# Save H1 data in H1H1 
	movaps [esp + nb314_rsqH1H1], xmm4
	
	;# do 1/x for O and 1/sqrt(x) for H1
	rcpss  xmm1, xmm0
	rsqrtps xmm5, xmm4
	movss  xmm2, [esp + nb314_two]
	movaps  xmm6, xmm5	
	mulss  xmm0, xmm1
	mulps   xmm5, xmm5
	subss  xmm2, xmm0
	movaps  xmm7, [esp + nb314_three]
	mulss  xmm2, xmm1 	;# constant 1/r2
	
	
	mulps   xmm5, xmm4
	movss  xmm0, xmm2
	subps   xmm7, xmm5
	mulss  xmm2, xmm2
	mulps   xmm7, xmm6							
	mulss  xmm2, xmm0 	;# constant 1/r6
	mulps   xmm7, [esp + nb314_half] ;# rinv iH1 - j water 
	movss  xmm1, xmm2
	movaps [esp + nb314_rinvH1H1], xmm7

	mulss  xmm2, xmm2 	;# constant 1/r12
	mulss  xmm1, [esp + nb314_c6]
	mulss  xmm2, [esp + nb314_c12]
	movss  xmm3, xmm2
	subss  xmm3, xmm1
	addss  xmm3, [esp + nb314_Vvdwtot]
	movss  [esp + nb314_Vvdwtot], xmm3
	mulss  xmm1, [esp + nb314_six]
	mulss  xmm2, [esp + nb314_twelve]
	subss  xmm2, xmm1
	mulss  xmm0, xmm2 	;# fscal
	movss  xmm1, xmm0
	movss  xmm2, xmm0
	mulss  xmm0, [esp + nb314_dxOO]
	mulss  xmm1, [esp + nb314_dyOO]
	mulss  xmm2, [esp + nb314_dzOO]
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2
	movaps  [esp + nb314_fjxO], xmm3
	movaps  [esp + nb314_fjyO], xmm4
	movaps  [esp + nb314_fjzO], xmm5
	addss   xmm0, [esp + nb314_fixO]
	addss   xmm1, [esp + nb314_fiyO]
	addss   xmm2, [esp + nb314_fizO]
	movss  [esp + nb314_fixO], xmm0
	movss  [esp + nb314_fiyO], xmm1
	movss  [esp + nb314_fizO], xmm2

	;# do  H1 coulomb interaction
	movaps xmm0, [esp + nb314_rinvH1H1] ;# rinv 
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb314_rsqH1H1] 	;# r
	mulps xmm1, [esp + nb314_tsc]
	
	movhlps xmm2, xmm1	
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 
	pslld   mm6, 2
	pslld   mm7, 2
	
	psrlq mm6, 32
	movd ebx, mm6
	movd ecx, mm7
	psrlq mm7, 32
	movd edx, mm7		;# table indices in ebx,ecx,edx 

	mov esi, [ebp + nb314_VFtab]

	movlps xmm4, [esi + ebx*4]
	movlps xmm3, [esi + ecx*4]
	movlps xmm7, [esi + edx*4]
	movhps xmm4, [esi + ebx*4 + 8]
	movhps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + edx*4 + 8]
	movaps xmm6, xmm3
	unpcklps xmm6, xmm7
	unpckhps xmm3, xmm7
	movaps xmm5, xmm4
	movaps xmm7, xmm4
	shufps xmm4, xmm6, 0x40
	shufps xmm5, xmm6, 0xE4
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 0x48
	shufps xmm7, xmm3, 0xEC
	;# coulomb table ready, in xmm4-xmm7
	
	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [esp + nb314_two]   	;# two*Heps2 

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [esp + nb314_qqHH]
	movhps  xmm3, [esp + nb314_qqMH]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point xmm5 contains vcoul and xmm3 fijC 
	
	addps  xmm5, [esp + nb314_vctot]
	movaps [esp + nb314_vctot], xmm5

	mulps  xmm3, [esp + nb314_tsc]
	xorps  xmm2, xmm2
	subps  xmm2, xmm3
	mulps  xmm0, xmm2
	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [esp + nb314_dxH1H1]
	mulps   xmm1, [esp + nb314_dyH1H1]
	mulps   xmm2, [esp + nb314_dzH1H1]
	;# update forces H1 - j water 
	movaps  xmm3, [esp + nb314_fjxO]
	movaps  xmm4, [esp + nb314_fjyO]
	movaps  xmm5, [esp + nb314_fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + nb314_fjxO], xmm3
	movaps  [esp + nb314_fjyO], xmm4
	movaps  [esp + nb314_fjzO], xmm5
	addps   xmm0, [esp + nb314_fixH1]
	addps   xmm1, [esp + nb314_fiyH1]
	addps   xmm2, [esp + nb314_fizH1]
	movaps  [esp + nb314_fixH1], xmm0
	movaps  [esp + nb314_fiyH1], xmm1
	movaps  [esp + nb314_fizH1], xmm2	
	
	;# i H2 & M simultaneously first get i particle coords: 
	movaps  xmm0, [esp + nb314_ixH2]
	movaps  xmm1, [esp + nb314_iyH2]
	movaps  xmm2, [esp + nb314_izH2]	
	movaps  xmm3, [esp + nb314_ixM] 
	movaps  xmm4, [esp + nb314_iyM] 
	movaps  xmm5, [esp + nb314_izM] 
	subps   xmm0, [esp + nb314_jxO]
	subps   xmm1, [esp + nb314_jyO]
	subps   xmm2, [esp + nb314_jzO]
	subps   xmm3, [esp + nb314_jxO]
	subps   xmm4, [esp + nb314_jyO]
	subps   xmm5, [esp + nb314_jzO]
	movaps [esp + nb314_dxH2H2], xmm0
	movaps [esp + nb314_dyH2H2], xmm1
	movaps [esp + nb314_dzH2H2], xmm2
	movaps [esp + nb314_dxMM], xmm3
	movaps [esp + nb314_dyMM], xmm4
	movaps [esp + nb314_dzMM], xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;# have rsqH2 in xmm0 
	addps xmm4, xmm5	;# have rsqM in xmm4 

	;# start with H2, save data 
	movaps [esp + nb314_rsqH2H2], xmm0
	movaps [esp + nb314_rsqMM], xmm4	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb314_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb314_half] ;# rinv H2 - j water 
	mulps   xmm7, [esp + nb314_half] ;# rinv M - j water  

	movaps [esp + nb314_rinvH2H2], xmm3
	movaps [esp + nb314_rinvMM], xmm7

	movaps xmm1, xmm3
	mulps  xmm1, [esp + nb314_rsqH2H2]	;# xmm1=r 
	movaps xmm0, xmm3	;# xmm0=rinv 
	mulps  xmm1, [esp + nb314_tsc]
	
	movhlps xmm2, xmm1	
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 
	pslld   mm6, 2
	pslld   mm7, 2

	psrlq mm6, 32
	movd ebx, mm6
	movd ecx, mm7
	psrlq mm7, 32
	movd edx, mm7		;# table indices in ebx,ecx,edx 

	movlps xmm4, [esi + ebx*4]
	movlps xmm3, [esi + ecx*4]
	movlps xmm7, [esi + edx*4]
	movhps xmm4, [esi + ebx*4 + 8]
	movhps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + edx*4 + 8]
	movaps xmm6, xmm3
	unpcklps xmm6, xmm7
	unpckhps xmm3, xmm7
	movaps xmm5, xmm4
	movaps xmm7, xmm4
	shufps xmm4, xmm6, 0x40
	shufps xmm5, xmm6, 0xE4
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 0x48
	shufps xmm7, xmm3, 0xEC
	;# coulomb table ready, in xmm4-xmm7
	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [esp + nb314_two]   	;# two*Heps2 

	xorps  xmm3, xmm3
	
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [esp + nb314_qqHH]
	movhps  xmm3, [esp + nb314_qqMH]
	shufps  xmm3, xmm3, 193 ;# constant 11000001
		
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point xmm5 contains vcoul and xmm3 fijC 
	addps  xmm5, [esp + nb314_vctot]
	movaps [esp + nb314_vctot], xmm5	

	xorps  xmm1, xmm1

	mulps xmm3, [esp + nb314_tsc]
	mulps xmm3, xmm0
	subps  xmm1, xmm3
	
	movaps  xmm0, xmm1
	movaps  xmm2, xmm1
	mulps   xmm0, [esp + nb314_dxH2H2]
	mulps   xmm1, [esp + nb314_dyH2H2]
	mulps   xmm2, [esp + nb314_dzH2H2]
	;# update forces H1 - j water 
	movaps  xmm3, [esp + nb314_fjxO]
	movaps  xmm4, [esp + nb314_fjyO]
	movaps  xmm5, [esp + nb314_fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + nb314_fjxO], xmm3
	movaps  [esp + nb314_fjyO], xmm4
	movaps  [esp + nb314_fjzO], xmm5
	addps   xmm0, [esp + nb314_fixH2]
	addps   xmm1, [esp + nb314_fiyH2]
	addps   xmm2, [esp + nb314_fizH2]
	movaps  [esp + nb314_fixH2], xmm0
	movaps  [esp + nb314_fiyH2], xmm1
	movaps  [esp + nb314_fizH2], xmm2
	
	;# do table for i M - j water interaction 
	movaps xmm0, [esp + nb314_rinvMM]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb314_rsqMM]	;# xmm0=rinv, xmm1=r 
	mulps  xmm1, [esp + nb314_tsc]
	
	movhlps xmm2, xmm1	
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 
	pslld   mm6, 2
	pslld   mm7, 2
 
	psrlq mm6, 32
	movd ebx, mm6
	movd ecx, mm7
	psrlq mm7, 32
	movd edx, mm7		;# table indices in ebx,ecx,edx 

	movlps xmm4, [esi + ebx*4]
	movlps xmm3, [esi + ecx*4]
	movlps xmm7, [esi + edx*4]
	movhps xmm4, [esi + ebx*4 + 8]
	movhps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + edx*4 + 8]
	movaps xmm6, xmm3
	unpcklps xmm6, xmm7
	unpckhps xmm3, xmm7
	movaps xmm5, xmm4
	movaps xmm7, xmm4
	shufps xmm4, xmm6, 0x40
	shufps xmm5, xmm6, 0xE4
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 0x48
	shufps xmm7, xmm3, 0xEC
	;# coulomb table ready, in xmm4-xmm7
	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [esp + nb314_two]   	;# two*Heps2 

	xorps  xmm3, xmm3
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [esp + nb314_qqMH]
	movhps  xmm3, [esp + nb314_qqMM]
	shufps  xmm3, xmm3, 193 ;# constant 11000001
		
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point xmm5 contains vcoul and xmm3 fijC 
	addps  xmm5, [esp + nb314_vctot]
	movaps [esp + nb314_vctot], xmm5	

	xorps  xmm1, xmm1

	mulps xmm3, [esp + nb314_tsc]
	mulps xmm3, xmm0
	subps  xmm1, xmm3
	
	movaps  xmm0, xmm1
	movaps  xmm2, xmm1
	
	mulps   xmm0, [esp + nb314_dxMM]
	mulps   xmm1, [esp + nb314_dyMM]
	mulps   xmm2, [esp + nb314_dzMM]
	movaps  xmm3, [esp + nb314_fjxO]
	movaps  xmm4, [esp + nb314_fjyO]
	movaps  xmm5, [esp + nb314_fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	mov     esi, [ebp + nb314_faction]
	movaps  [esp + nb314_fjxO], xmm3
	movaps  [esp + nb314_fjyO], xmm4
	movaps  [esp + nb314_fjzO], xmm5
	addps   xmm0, [esp + nb314_fixM]
	addps   xmm1, [esp + nb314_fiyM]
	addps   xmm2, [esp + nb314_fizM]
	movaps  [esp + nb314_fixM], xmm0
	movaps  [esp + nb314_fiyM], xmm1
	movaps  [esp + nb314_fizM], xmm2

	;# update j water forces from local variables.
	;# transpose back first
	movaps  xmm0, [esp + nb314_fjxO] ;# Ox H1x H2x Mx 
	movaps  xmm1, [esp + nb314_fjyO] ;# Oy H1y H2y My
	movaps  xmm2, [esp + nb314_fjzO] ;# Oz H1z H2z Mz
	 
	movaps  xmm3, xmm0
	movaps  xmm4, xmm0
	unpcklps xmm3, xmm1       	;# Ox Oy - -
	shufps  xmm4, xmm2, 0x1	  	;# h1x - Oz -
	movaps  xmm5, xmm1
	movaps  xmm6, xmm0
	unpcklps xmm5, xmm2 	  	;# - - H1y H1z
	unpckhps xmm6, xmm1	  	;# h2x h2y - - 
	unpckhps xmm1, xmm2	  	;# - - My Mz
	
	shufps   xmm2, xmm0, 0x32  ;# (00110010) h2z - Mx -
	shufps   xmm3, xmm4, 36  ;# constant 00100100 ;# Ox Oy Oz H1x 
	shufps   xmm5, xmm6, 78  ;# constant 01001110 ;# h1y h1z h2x h2y
	shufps   xmm2, xmm1, 232  ;# constant 11101000 ;# h2z mx my mz

	movlps  xmm0, [esi + eax*4]
	movlps  xmm1, [esi + eax*4 + 16]
	movlps  xmm4, [esi + eax*4 + 32]
	movhps  xmm0, [esi + eax*4 + 8]
	movhps  xmm1, [esi + eax*4 + 24]
	movhps  xmm4, [esi + eax*4 + 40]
	addps   xmm0, xmm3
	addps   xmm1, xmm5
	addps   xmm4, xmm2
	movlps   [esi + eax*4], xmm0
	movlps   [esi + eax*4 + 16], xmm1
	movlps   [esi + eax*4 + 32], xmm4
	movhps   [esi + eax*4 + 8], xmm0
	movhps   [esi + eax*4 + 24], xmm1
	movhps   [esi + eax*4 + 40], xmm4

	dec dword ptr [esp + nb314_innerk]
	jz    .nb314_updateouterdata
	jmp   .nb314_single_loop
.nb314_updateouterdata:
	mov   ecx, [esp + nb314_ii3]
	mov   edi, [ebp + nb314_faction]
	mov   esi, [ebp + nb314_fshift]
	mov   edx, [esp + nb314_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb314_fixO]
	movaps xmm1, [esp + nb314_fiyO] 
	movaps xmm2, [esp + nb314_fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# constant 00001000	

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb314_fixH1]
	movaps xmm1, [esp + nb314_fiyH1]
	movaps xmm2, [esp + nb314_fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [edi + ecx*4 + 12]
	movss  xmm4, [edi + ecx*4 + 16]
	movss  xmm5, [edi + ecx*4 + 20]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 12], xmm3
	movss  [edi + ecx*4 + 16], xmm4
	movss  [edi + ecx*4 + 20], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb314_fixH2]
	movaps xmm1, [esp + nb314_fiyH2]
	movaps xmm2, [esp + nb314_fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [edi + ecx*4 + 24]
	movss  xmm4, [edi + ecx*4 + 28]
	movss  xmm5, [edi + ecx*4 + 32]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 24], xmm3
	movss  [edi + ecx*4 + 28], xmm4
	movss  [edi + ecx*4 + 32], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# accumulate Mi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb314_fixM]
	movaps xmm1, [esp + nb314_fiyM]
	movaps xmm2, [esp + nb314_fizM]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [edi + ecx*4 + 36]
	movss  xmm4, [edi + ecx*4 + 40]
	movss  xmm5, [edi + ecx*4 + 44]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 36], xmm3
	movss  [edi + ecx*4 + 40], xmm4
	movss  [edi + ecx*4 + 44], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# increment fshift force  
	movlps  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [esi + edx*4],    xmm3
	movss  [esi + edx*4 + 8], xmm4

	;# get n from stack
	mov esi, [esp + nb314_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb314_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb314_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb314_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [esp + nb314_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb314_Vvdw]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb314_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb314_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb314_n], esi
        jmp .nb314_outer
.nb314_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb314_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb314_end
        ;# non-zero, do one more workunit
        jmp   .nb314_threadloop
.nb314_end:
	emms

	mov eax, [esp + nb314_nouter]
	mov ebx, [esp + nb314_ninner]
	mov ecx, [ebp + nb314_outeriter]
	mov edx, [ebp + nb314_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb314_salign]
	add esp, eax
	add esp, 1832
	pop edi
	pop esi
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret

	


.globl nb_kernel314nf_ia32_sse
.globl _nb_kernel314nf_ia32_sse
nb_kernel314nf_ia32_sse:	
_nb_kernel314nf_ia32_sse:	
.equiv          nb314nf_p_nri,          8
.equiv          nb314nf_iinr,           12
.equiv          nb314nf_jindex,         16
.equiv          nb314nf_jjnr,           20
.equiv          nb314nf_shift,          24
.equiv          nb314nf_shiftvec,       28
.equiv          nb314nf_fshift,         32
.equiv          nb314nf_gid,            36
.equiv          nb314nf_pos,            40
.equiv          nb314nf_faction,        44
.equiv          nb314nf_charge,         48
.equiv          nb314nf_p_facel,        52
.equiv          nb314nf_argkrf,         56
.equiv          nb314nf_argcrf,         60
.equiv          nb314nf_Vc,             64
.equiv          nb314nf_type,           68
.equiv          nb314nf_p_ntype,        72
.equiv          nb314nf_vdwparam,       76
.equiv          nb314nf_Vvdw,           80
.equiv          nb314nf_p_tabscale,     84
.equiv          nb314nf_VFtab,          88
.equiv          nb314nf_invsqrta,       92
.equiv          nb314nf_dvda,           96
.equiv          nb314nf_p_gbtabscale,   100
.equiv          nb314nf_GBtab,          104
.equiv          nb314nf_p_nthreads,     108
.equiv          nb314nf_count,          112
.equiv          nb314nf_mtx,            116
.equiv          nb314nf_outeriter,      120
.equiv          nb314nf_inneriter,      124
.equiv          nb314nf_work,           128
;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb314nf_ixO,            0
.equiv          nb314nf_iyO,            16
.equiv          nb314nf_izO,            32
.equiv          nb314nf_ixH1,           48
.equiv          nb314nf_iyH1,           64
.equiv          nb314nf_izH1,           80
.equiv          nb314nf_ixH2,           96
.equiv          nb314nf_iyH2,           112
.equiv          nb314nf_izH2,           128
.equiv          nb314nf_ixM,            144
.equiv          nb314nf_iyM,            160
.equiv          nb314nf_izM,            176
.equiv          nb314nf_jxO,            192
.equiv          nb314nf_jyO,            208
.equiv          nb314nf_jzO,            224
.equiv          nb314nf_jxH1,           240
.equiv          nb314nf_jyH1,           256
.equiv          nb314nf_jzH1,           272
.equiv          nb314nf_jxH2,           288
.equiv          nb314nf_jyH2,           304
.equiv          nb314nf_jzH2,           320
.equiv          nb314nf_jxM,            336
.equiv          nb314nf_jyM,            352
.equiv          nb314nf_jzM,            368
.equiv          nb314nf_qqMM,           384
.equiv          nb314nf_qqMH,           400
.equiv          nb314nf_qqHH,           416
.equiv          nb314nf_two,            432
.equiv          nb314nf_tsc,            448
.equiv          nb314nf_c6,             464
.equiv          nb314nf_c12,            480
.equiv          nb314nf_vctot,          496
.equiv          nb314nf_Vvdwtot,        512
.equiv          nb314nf_half,           528
.equiv          nb314nf_three,          544
.equiv          nb314nf_rsqOO,          560
.equiv          nb314nf_rsqH1H1,        576
.equiv          nb314nf_rsqH1H2,        592
.equiv          nb314nf_rsqH1M,         608
.equiv          nb314nf_rsqH2H1,        624
.equiv          nb314nf_rsqH2H2,        640
.equiv          nb314nf_rsqH2M,         656
.equiv          nb314nf_rsqMH1,         672
.equiv          nb314nf_rsqMH2,         688
.equiv          nb314nf_rsqMM,          704
.equiv          nb314nf_rinvsqOO,       720
.equiv          nb314nf_rinvH1H1,       736
.equiv          nb314nf_rinvH1H2,       752
.equiv          nb314nf_rinvH1M,        768
.equiv          nb314nf_rinvH2H1,       784
.equiv          nb314nf_rinvH2H2,       800
.equiv          nb314nf_rinvH2M,        816
.equiv          nb314nf_rinvMH1,        832
.equiv          nb314nf_rinvMH2,        848
.equiv          nb314nf_rinvMM,         864
.equiv          nb314nf_is3,            880
.equiv          nb314nf_ii3,            884
.equiv          nb314nf_innerjjnr,      888
.equiv          nb314nf_innerk,         892
.equiv          nb314nf_n,              896
.equiv          nb314nf_nn1,            900
.equiv          nb314nf_nri,            904
.equiv          nb314nf_nouter,         908
.equiv          nb314nf_ninner,         912
.equiv          nb314nf_salign,         916
	push ebp
	mov ebp,esp	
	push eax
	push ebx
	push ecx
	push edx
	push esi
	push edi
	sub esp, 920		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb314nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb314nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb314nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb314nf_nouter], eax
	mov [esp + nb314nf_ninner], eax


	mov eax, [ebp + nb314nf_p_tabscale]
	movss xmm5, [eax]
	shufps xmm5, xmm5, 0
	movaps [esp + nb314nf_tsc],  xmm5
	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb314nf_half], eax
	movss xmm1, [esp + nb314nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb314nf_half],  xmm1
	movaps [esp + nb314nf_two],  xmm2
	movaps [esp + nb314nf_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb314nf_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]		;# ebx =ii 

	mov   edx, [ebp + nb314nf_charge]
	movss xmm5, [edx + ebx*4 + 4]	
	movss xmm3, [edx + ebx*4 + 12]	
	movss xmm4, xmm3	
	mov esi, [ebp + nb314nf_p_facel]
	movss xmm6, [esi]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb314nf_qqMM], xmm3
	movaps [esp + nb314nf_qqMH], xmm4
	movaps [esp + nb314nf_qqHH], xmm5
		
	xorps xmm0, xmm0
	mov   edx, [ebp + nb314nf_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov edi, [ebp + nb314nf_p_ntype]
	imul  ecx, [edi]  	;# ecx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   eax, [ebp + nb314nf_vdwparam]
	movlps xmm0, [eax + edx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0x55
	movaps [esp + nb314nf_c6], xmm0
	movaps [esp + nb314nf_c12], xmm1

.nb314nf_threadloop:
        mov   esi, [ebp + nb314nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb314nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb314nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb314nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb314nf_n], eax
        mov [esp + nb314nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb314nf_outerstart
        jmp .nb314nf_end
	
.nb314nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb314nf_nouter]
	mov [esp + nb314nf_nouter], ebx

.nb314nf_outer:
	mov   eax, [ebp + nb314nf_shift]  	;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]	;# ebx=3*is 
	mov   [esp + nb314nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb314nf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb314nf_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]		;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb314nf_pos]	;# eax = base of pos[]  
	mov   [esp + nb314nf_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	movaps xmm6, xmm0
	movaps xmm7, xmm1

	addss xmm3, [eax + ebx*4]  	;# ox
	addss xmm4, [eax + ebx*4 + 4]  ;# oy
	addss xmm5, [eax + ebx*4 + 8]  ;# oz
	addss xmm6, [eax + ebx*4 + 12] ;# h1x
	addss xmm7, [eax + ebx*4 + 16] ;# h1y
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	shufps xmm7, xmm7, 0
	movaps [esp + nb314nf_ixO], xmm3
	movaps [esp + nb314nf_iyO], xmm4
	movaps [esp + nb314nf_izO], xmm5
	movaps [esp + nb314nf_ixH1], xmm6
	movaps [esp + nb314nf_iyH1], xmm7

	movss xmm6, xmm2
	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm6, [eax + ebx*4 + 20] ;# h1z
	addss xmm0, [eax + ebx*4 + 24] ;# h2x
	addss xmm1, [eax + ebx*4 + 28] ;# h2y
	addss xmm2, [eax + ebx*4 + 32] ;# h2z
	addss xmm3, [eax + ebx*4 + 36] ;# mx
	addss xmm4, [eax + ebx*4 + 40] ;# my
	addss xmm5, [eax + ebx*4 + 44] ;# mz

	shufps xmm6, xmm6, 0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb314nf_izH1], xmm6
	movaps [esp + nb314nf_ixH2], xmm0
	movaps [esp + nb314nf_iyH2], xmm1
	movaps [esp + nb314nf_izH2], xmm2
	movaps [esp + nb314nf_ixM], xmm3
	movaps [esp + nb314nf_iyM], xmm4
	movaps [esp + nb314nf_izM], xmm5

	;# clear vctot 	
	xorps xmm4, xmm4
	movaps [esp + nb314nf_vctot], xmm4
	movaps [esp + nb314nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb314nf_jindex]
	mov   ecx, [eax + esi*4]	 	;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   esi, [ebp + nb314nf_pos]
	mov   edi, [ebp + nb314nf_faction]	
	mov   eax, [ebp + nb314nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb314nf_innerjjnr], eax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb314nf_ninner]
	mov   [esp + nb314nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb314nf_innerk], edx	;# number of innerloop atoms 
	jge   .nb314nf_unroll_loop
	jmp   .nb314nf_single_check
.nb314nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb314nf_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	mov   ecx, [edx + 8]
	mov   edx, [edx + 12]     	;# eax-edx=jnr1-4 
	
	add dword ptr [esp + nb314nf_innerjjnr],  16 ;# advance pointer (unroll 4) 

	mov esi, [ebp + nb314nf_pos]   	;# base of pos[] 

	lea   eax, [eax + eax*2] 	;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2] 	;# replace jnr with j3 
	lea   edx, [edx + edx*2]	
	
	;# move j coordinates to local temp variables
	;# Load Ox, Oy, Oz, H1x 
	movlps xmm1, [esi + eax*4]	;#  Oxa   Oya    -    -
	movlps xmm4, [esi + ecx*4]	;#  Oxc   Oyc    -    -
	movhps xmm1, [esi + ebx*4]	;#  Oxa   Oya   Oxb   Oyb 
	movhps xmm4, [esi + edx*4]	;#  Oxc   Oyc   Oxd   Oyd 
	movaps xmm0, xmm1		;#  Oxa   Oya   Oxb   Oyb 
	shufps xmm0, xmm4, 0x88		;#  Oxa   Oxb   Oxc   Oxd
	shufps xmm1, xmm4, 0xDD		;#  Oya   Oyb   Oyc   Oyd
	movlps xmm3, [esi + eax*4 + 8]	;#  Oza  H1xa    -    -
	movlps xmm5, [esi + ecx*4 + 8]	;#  Ozc  H1xc    -    -
	movhps xmm3, [esi + ebx*4 + 8]	;#  Oza  H1xa   Ozb  H1xb 
	movhps xmm5, [esi + edx*4 + 8]	;#  Ozc  H1xc   Ozd  H1xd 
	movaps xmm2, xmm3		;#  Oza  H1xa   Ozb  H1xb 
	shufps xmm2, xmm5, 0x88		;#  Oza   Ozb   Ozc   Ozd
	shufps xmm3, xmm5, 0xDD		;# H1xa  H1xb  H1xc  H1xd
	;# coordinates in xmm0-xmm3	
	;# store
	movaps [esp + nb314nf_jxO], xmm0
	movaps [esp + nb314nf_jyO], xmm1
	movaps [esp + nb314nf_jzO], xmm2
	movaps [esp + nb314nf_jxH1], xmm3

	;# Load H1y H1z H2x H2y 
	movlps xmm1, [esi + eax*4 + 16]	
	movlps xmm4, [esi + ecx*4 + 16]	
	movhps xmm1, [esi + ebx*4 + 16]	
	movhps xmm4, [esi + edx*4 + 16]	
	movaps xmm0, xmm1		
	shufps xmm0, xmm4, 0x88		
	shufps xmm1, xmm4, 0xDD		
	movlps xmm3, [esi + eax*4 + 24]	
	movlps xmm5, [esi + ecx*4 + 24]	
	movhps xmm3, [esi + ebx*4 + 24]	
	movhps xmm5, [esi + edx*4 + 24]	
	movaps xmm2, xmm3		
	shufps xmm2, xmm5, 0x88		
	shufps xmm3, xmm5, 0xDD		
	;# coordinates in xmm0-xmm3	
	;# store
	movaps [esp + nb314nf_jyH1], xmm0
	movaps [esp + nb314nf_jzH1], xmm1
	movaps [esp + nb314nf_jxH2], xmm2
	movaps [esp + nb314nf_jyH2], xmm3

	;# Load H2z Mx My Mz 
	movlps xmm1, [esi + eax*4 + 32]	
	movlps xmm4, [esi + ecx*4 + 32]	
	movhps xmm1, [esi + ebx*4 + 32]	
	movhps xmm4, [esi + edx*4 + 32]	
	movaps xmm0, xmm1		
	shufps xmm0, xmm4, 0x88		
	shufps xmm1, xmm4, 0xDD		
	movlps xmm3, [esi + eax*4 + 40]	
	movlps xmm5, [esi + ecx*4 + 40]	
	movhps xmm3, [esi + ebx*4 + 40]	
	movhps xmm5, [esi + edx*4 + 40]	
	movaps xmm2, xmm3		
	shufps xmm2, xmm5, 0x88		
	shufps xmm3, xmm5, 0xDD		
	;# coordinates in xmm0-xmm3	
	;# store
	movaps [esp + nb314nf_jzH2], xmm0
	movaps [esp + nb314nf_jxM], xmm1
	movaps [esp + nb314nf_jyM], xmm2
	movaps [esp + nb314nf_jzM], xmm3
	
	;# start calculating pairwise distances
	movaps xmm0, [esp + nb314nf_ixO]
	movaps xmm1, [esp + nb314nf_iyO]
	movaps xmm2, [esp + nb314nf_izO]
	movaps xmm3, [esp + nb314nf_ixH1]
	movaps xmm4, [esp + nb314nf_iyH1]
	movaps xmm5, [esp + nb314nf_izH1]
	subps  xmm0, [esp + nb314nf_jxO]
	subps  xmm1, [esp + nb314nf_jyO]
	subps  xmm2, [esp + nb314nf_jzO]
	subps  xmm3, [esp + nb314nf_jxH1]
	subps  xmm4, [esp + nb314nf_jyH1]
	subps  xmm5, [esp + nb314nf_jzH1]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb314nf_rsqOO], xmm0
	movaps [esp + nb314nf_rsqH1H1], xmm3

	movaps xmm0, [esp + nb314nf_ixH1]
	movaps xmm1, [esp + nb314nf_iyH1]
	movaps xmm2, [esp + nb314nf_izH1]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [esp + nb314nf_jxH2]
	subps  xmm1, [esp + nb314nf_jyH2]
	subps  xmm2, [esp + nb314nf_jzH2]
	subps  xmm3, [esp + nb314nf_jxM]
	subps  xmm4, [esp + nb314nf_jyM]
	subps  xmm5, [esp + nb314nf_jzM]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb314nf_rsqH1H2], xmm0
	movaps [esp + nb314nf_rsqH1M], xmm3

	movaps xmm0, [esp + nb314nf_ixH2]
	movaps xmm1, [esp + nb314nf_iyH2]
	movaps xmm2, [esp + nb314nf_izH2]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [esp + nb314nf_jxH1]
	subps  xmm1, [esp + nb314nf_jyH1]
	subps  xmm2, [esp + nb314nf_jzH1]
	subps  xmm3, [esp + nb314nf_jxH2]
	subps  xmm4, [esp + nb314nf_jyH2]
	subps  xmm5, [esp + nb314nf_jzH2]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb314nf_rsqH2H1], xmm0
	movaps [esp + nb314nf_rsqH2H2], xmm3

	movaps xmm0, [esp + nb314nf_ixH2]
	movaps xmm1, [esp + nb314nf_iyH2]
	movaps xmm2, [esp + nb314nf_izH2]
	movaps xmm3, [esp + nb314nf_ixM]
	movaps xmm4, [esp + nb314nf_iyM]
	movaps xmm5, [esp + nb314nf_izM]
	subps  xmm0, [esp + nb314nf_jxM]
	subps  xmm1, [esp + nb314nf_jyM]
	subps  xmm2, [esp + nb314nf_jzM]
	subps  xmm3, [esp + nb314nf_jxH1]
	subps  xmm4, [esp + nb314nf_jyH1]
	subps  xmm5, [esp + nb314nf_jzH1]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [esp + nb314nf_rsqH2M], xmm0
	movaps [esp + nb314nf_rsqMH1], xmm4

	movaps xmm0, [esp + nb314nf_ixM]
	movaps xmm1, [esp + nb314nf_iyM]
	movaps xmm2, [esp + nb314nf_izM]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [esp + nb314nf_jxH2]
	subps  xmm1, [esp + nb314nf_jyH2]
	subps  xmm2, [esp + nb314nf_jzH2]
	subps  xmm3, [esp + nb314nf_jxM]
	subps  xmm4, [esp + nb314nf_jyM]
	subps  xmm5, [esp + nb314nf_jzM]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [esp + nb314nf_rsqMH2], xmm0
	movaps [esp + nb314nf_rsqMM], xmm4

	;# start by doing reciprocal for OO
	movaps  xmm7, [esp + nb314nf_rsqOO]
	rcpps   xmm2, xmm7
	movaps  xmm1, [esp + nb314nf_two]
	mulps   xmm7, xmm2
	subps   xmm1, xmm7
	mulps   xmm2, xmm1 ;# rinvsq 
	movaps [esp + nb314nf_rinvsqOO], xmm2
	
	;# next step is invsqrt - do two at a time.
	rsqrtps xmm1, [esp + nb314nf_rsqH1H1]
	rsqrtps xmm5, [esp + nb314nf_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb314nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb314nf_rsqH1H1]
	mulps   xmm5, [esp + nb314nf_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb314nf_half] ;# rinvH1H1 
	mulps   xmm7, [esp + nb314nf_half] ;# rinvH1H2 
	movaps  [esp + nb314nf_rinvH1H1], xmm3
	movaps  [esp + nb314nf_rinvH1H2], xmm7
			
	rsqrtps xmm1, [esp + nb314nf_rsqH1M]
	rsqrtps xmm5, [esp + nb314nf_rsqH2H1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb314nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb314nf_rsqH1M]
	mulps   xmm5, [esp + nb314nf_rsqH2H1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb314nf_half] 
	mulps   xmm7, [esp + nb314nf_half]
	movaps  [esp + nb314nf_rinvH1M], xmm3
	movaps  [esp + nb314nf_rinvH2H1], xmm7
			
	rsqrtps xmm1, [esp + nb314nf_rsqH2H2]
	rsqrtps xmm5, [esp + nb314nf_rsqH2M]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb314nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb314nf_rsqH2H2]
	mulps   xmm5, [esp + nb314nf_rsqH2M]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb314nf_half] 
	mulps   xmm7, [esp + nb314nf_half]
	movaps  [esp + nb314nf_rinvH2H2], xmm3
	movaps  [esp + nb314nf_rinvH2M], xmm7
	
	rsqrtps xmm1, [esp + nb314nf_rsqMH1]
	rsqrtps xmm5, [esp + nb314nf_rsqMH2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb314nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb314nf_rsqMH1]
	mulps   xmm5, [esp + nb314nf_rsqMH2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb314nf_half] 
	mulps   xmm7, [esp + nb314nf_half]
	movaps  [esp + nb314nf_rinvMH1], xmm3
	movaps  [esp + nb314nf_rinvMH2], xmm7
        		
	rsqrtps xmm1, [esp + nb314nf_rsqMM]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + nb314nf_three]
	mulps   xmm1, [esp + nb314nf_rsqMM]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [esp + nb314nf_half] 
	movaps  [esp + nb314nf_rinvMM], xmm3
	
	;# start with OO LJ interaction
	movaps xmm0, [esp + nb314nf_rinvsqOO]
	movaps xmm1, xmm0
	mulps  xmm1, xmm1	;# rinv4
	mulps  xmm1, xmm0	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm1, [esp + nb314nf_c6]
	mulps  xmm2, [esp + nb314nf_c12]
	movaps xmm4, xmm2
	subps  xmm4, xmm1
	addps  xmm4, [esp + nb314nf_Vvdwtot]
	movaps [esp + nb314nf_Vvdwtot], xmm4

	;# Coulomb interactions - first H1H1
	movaps xmm0, [esp + nb314nf_rinvH1H1]
	
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb314nf_rsqH1H1] ;# xmm1=r 
	mulps  xmm1, [esp + nb314nf_tsc]
		
	movhlps xmm2, xmm1
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 
	pslld   mm6, 2
	pslld   mm7, 2
	
	movd mm0, eax
	movd mm1, ebx
	movd mm2, ecx
	movd mm3, edx

	mov  esi, [ebp + nb314nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	movaps xmm3, [esp + nb314nf_qqHH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 
	;# update vctot 
        addps  xmm5, [esp + nb314nf_vctot]
        movaps [esp + nb314nf_vctot], xmm5
		
	;# H1-H2 interaction 
	movaps xmm0, [esp + nb314nf_rinvH1H2]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb314nf_rsqH1H2] ;# xmm1=r 
	mulps  xmm1, [esp + nb314nf_tsc]	
	movhlps xmm2, xmm1
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2
	
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	movaps xmm3, [esp + nb314nf_qqHH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 

	addps  xmm5, [esp + nb314nf_vctot]
	movaps [esp + nb314nf_vctot], xmm5

	;# H1-M interaction  
	movaps xmm0, [esp + nb314nf_rinvH1M]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb314nf_rsqH1M] ;# xmm1=r 
	mulps  xmm1, [esp + nb314nf_tsc]	
	movhlps xmm2, xmm1
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2

	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	movaps xmm3, [esp + nb314nf_qqMH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 

	addps  xmm5, [esp + nb314nf_vctot]
	movaps [esp + nb314nf_vctot], xmm5

	;# H2-H1 interaction 
	movaps xmm0, [esp + nb314nf_rinvH2H1]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb314nf_rsqH2H1] ;# xmm1=r 
	mulps  xmm1, [esp + nb314nf_tsc]	
	movhlps xmm2, xmm1
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2

	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	movaps xmm3, [esp + nb314nf_qqHH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 

	addps  xmm5, [esp + nb314nf_vctot]
	movaps [esp + nb314nf_vctot], xmm5

	;# H2-H2 interaction 
	movaps xmm0, [esp + nb314nf_rinvH2H2]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb314nf_rsqH2H2] ;# xmm1=r 
	mulps  xmm1, [esp + nb314nf_tsc]	
	movhlps xmm2, xmm1
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2

	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	movaps xmm3, [esp + nb314nf_qqHH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 

	addps  xmm5, [esp + nb314nf_vctot]
	movaps [esp + nb314nf_vctot], xmm5

	;# H2-M interaction 
	movaps xmm0, [esp + nb314nf_rinvH2M]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb314nf_rsqH2M] ;# xmm1=r 
	mulps  xmm1, [esp + nb314nf_tsc]
	movhlps xmm2, xmm1
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2

	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	movaps xmm3, [esp + nb314nf_qqMH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 

	addps  xmm5, [esp + nb314nf_vctot]
	movaps [esp + nb314nf_vctot], xmm5

	;# M-H1 interaction 
	movaps xmm0, [esp + nb314nf_rinvMH1]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb314nf_rsqMH1] ;# xmm1=r 
	mulps  xmm1, [esp + nb314nf_tsc]	
	movhlps xmm2, xmm1
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2

	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	movaps xmm3, [esp + nb314nf_qqMH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 

	addps  xmm5, [esp + nb314nf_vctot]
	movaps [esp + nb314nf_vctot], xmm5

	;# M-H2 interaction 
	movaps xmm0, [esp + nb314nf_rinvMH2]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb314nf_rsqMH2] ;# xmm1=r 
	mulps  xmm1, [esp + nb314nf_tsc]
	movhlps xmm2, xmm1
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2

	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	movaps xmm3, [esp + nb314nf_qqMH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 

	addps  xmm5, [esp + nb314nf_vctot]
	movaps [esp + nb314nf_vctot], xmm5

	;# M-M interaction 
	movaps xmm0, [esp + nb314nf_rinvMM]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb314nf_rsqMM] ;# xmm1=r 
	mulps  xmm1, [esp + nb314nf_tsc]	
	movhlps xmm2, xmm1
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2

	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	movaps xmm3, [esp + nb314nf_qqMM]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 

	addps  xmm5, [esp + nb314nf_vctot]
	movaps [esp + nb314nf_vctot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [esp + nb314nf_innerk],  4
	jl    .nb314nf_single_check
	jmp   .nb314nf_unroll_loop
.nb314nf_single_check:
	add dword ptr [esp + nb314nf_innerk],  4
	jnz   .nb314nf_single_loop
	jmp   .nb314nf_updateouterdata
.nb314nf_single_loop:
	mov   edx, [esp + nb314nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb314nf_innerjjnr],  4	

	mov esi, [ebp + nb314nf_pos]
	lea   eax, [eax + eax*2]  

	;# fetch j coordinates
	movlps xmm3,  [esi + eax*4]		;#  Ox  Oy  
	movlps xmm4,  [esi + eax*4 + 16]	;# H1y H1z 
	movlps xmm5,  [esi + eax*4 + 32]	;# H2z  Mx 
	movhps xmm3,  [esi + eax*4 + 8]   	;#  Ox  Oy  Oz H1x
	movhps xmm4,  [esi + eax*4 + 24]	;# H1y H1z H2x H2y
	movhps xmm5,  [esi + eax*4 + 40]	;# H2z  Mx  My  Mz
	;# transpose
	movaps xmm0, xmm4
	movaps xmm1, xmm3
	movaps xmm2, xmm4
	movaps xmm6, xmm3
	shufps xmm4, xmm5, 18  ;# (00010010)  h2x - Mx  - 
	shufps xmm3, xmm0, 193 ;# (11000001)  Oy  - H1y - 
	shufps xmm2, xmm5, 35  ;# (00100011) H2y - My  - 
 	shufps xmm1, xmm0, 18  ;# (00010010)  Oz  - H1z - 
	;#  xmm6: Ox - - H1x   xmm5: H2z - - Mz 
	shufps xmm6, xmm4, 140 ;# (10001100) Ox H1x H2x Mx 
	shufps xmm3, xmm2, 136 ;# (10001000) Oy H1y H2y My 
	shufps xmm1, xmm5, 200 ;# (11001000) Oz H1z H2z Mz

	;# store all j coordinates in jO  
	movaps [esp + nb314nf_jxO], xmm6
	movaps [esp + nb314nf_jyO], xmm3
	movaps [esp + nb314nf_jzO], xmm1

	;# do O and H1 in parallel
	movaps xmm0, [esp + nb314nf_ixO]
	movaps xmm1, [esp + nb314nf_iyO]
	movaps xmm2, [esp + nb314nf_izO]
	movaps xmm3, [esp + nb314nf_ixH1]
	movaps xmm4, [esp + nb314nf_iyH1]
	movaps xmm5, [esp + nb314nf_izH1]
	subps  xmm0, [esp + nb314nf_jxO]
	subps  xmm1, [esp + nb314nf_jyO]
	subps  xmm2, [esp + nb314nf_jzO]
	subps  xmm3, [esp + nb314nf_jxO]
	subps  xmm4, [esp + nb314nf_jyO]
	subps  xmm5, [esp + nb314nf_jzO]
	
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm4, xmm3
	addps xmm4, xmm5	;# have rsq in xmm4
	;# Save H1 data in H1H1 
	movaps [esp + nb314nf_rsqH1H1], xmm4
	
	;# do 1/x for O and 1/sqrt(x) for H1
	rcpss  xmm1, xmm0
	rsqrtps xmm5, xmm4
	movss  xmm2, [esp + nb314nf_two]
	movaps  xmm6, xmm5	
	mulss  xmm0, xmm1
	mulps   xmm5, xmm5
	subss  xmm2, xmm0
	movaps  xmm7, [esp + nb314nf_three]
	mulss  xmm2, xmm1 	;# constant 1/r2
	
	
	mulps   xmm5, xmm4
	movss  xmm0, xmm2
	subps   xmm7, xmm5
	mulss  xmm2, xmm2
	mulps   xmm7, xmm6							
	mulss  xmm2, xmm0 	;# constant 1/r6
	mulps   xmm7, [esp + nb314nf_half] ;# rinv iH1 - j water 
	movss  xmm1, xmm2
	movaps [esp + nb314nf_rinvH1H1], xmm7

	mulss  xmm2, xmm2 	;# constant 1/r12
	mulss  xmm1, [esp + nb314nf_c6]
	mulss  xmm2, [esp + nb314nf_c12]
	movss  xmm3, xmm2
	subss  xmm3, xmm1
	addss  xmm3, [esp + nb314nf_Vvdwtot]
	movss  [esp + nb314nf_Vvdwtot], xmm3

	;# do  H1 coulomb interaction
	movaps xmm0, [esp + nb314nf_rinvH1H1] ;# rinv 
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb314nf_rsqH1H1] 	;# r
	mulps xmm1, [esp + nb314nf_tsc]
	
	movhlps xmm2, xmm1	
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 
	pslld   mm6, 2
	pslld   mm7, 2
	
	psrlq mm6, 32
	movd ebx, mm6
	movd ecx, mm7
	psrlq mm7, 32
	movd edx, mm7		;# table indices in ebx,ecx,edx 

	mov esi, [ebp + nb314nf_VFtab]

	movlps xmm4, [esi + ebx*4]
	movlps xmm3, [esi + ecx*4]
	movlps xmm7, [esi + edx*4]
	movhps xmm4, [esi + ebx*4 + 8]
	movhps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + edx*4 + 8]
	movaps xmm6, xmm3
	unpcklps xmm6, xmm7
	unpckhps xmm3, xmm7
	movaps xmm5, xmm4
	movaps xmm7, xmm4
	shufps xmm4, xmm6, 0x40
	shufps xmm5, xmm6, 0xE4
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 0x48
	shufps xmm7, xmm3, 0xEC
	;# coulomb table ready, in xmm4-xmm7
	
	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [esp + nb314nf_qqHH]
	movhps  xmm3, [esp + nb314nf_qqMH]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point xmm5 contains vcoul 
	
	addps  xmm5, [esp + nb314nf_vctot]
	movaps [esp + nb314nf_vctot], xmm5
	
	;# i H2 & M simultaneously first get i particle coords: 
	movaps  xmm0, [esp + nb314nf_ixH2]
	movaps  xmm1, [esp + nb314nf_iyH2]
	movaps  xmm2, [esp + nb314nf_izH2]	
	movaps  xmm3, [esp + nb314nf_ixM] 
	movaps  xmm4, [esp + nb314nf_iyM] 
	movaps  xmm5, [esp + nb314nf_izM] 
	subps   xmm0, [esp + nb314nf_jxO]
	subps   xmm1, [esp + nb314nf_jyO]
	subps   xmm2, [esp + nb314nf_jzO]
	subps   xmm3, [esp + nb314nf_jxO]
	subps   xmm4, [esp + nb314nf_jyO]
	subps   xmm5, [esp + nb314nf_jzO]
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;# have rsqH2 in xmm0 
	addps xmm4, xmm5	;# have rsqM in xmm4 

	;# start with H2, save data 
	movaps [esp + nb314nf_rsqH2H2], xmm0
	movaps [esp + nb314nf_rsqMM], xmm4	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb314nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb314nf_half] ;# rinv H2 - j water 
	mulps   xmm7, [esp + nb314nf_half] ;# rinv M - j water  

	movaps [esp + nb314nf_rinvH2H2], xmm3
	movaps [esp + nb314nf_rinvMM], xmm7

	movaps xmm1, xmm3
	mulps  xmm1, [esp + nb314nf_rsqH2H2]	;# xmm1=r 
	movaps xmm0, xmm3	;# xmm0=rinv 
	mulps  xmm1, [esp + nb314nf_tsc]
	
	movhlps xmm2, xmm1	
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 
	pslld   mm6, 2
	pslld   mm7, 2

	psrlq mm6, 32
	movd ebx, mm6
	movd ecx, mm7
	psrlq mm7, 32
	movd edx, mm7		;# table indices in ebx,ecx,edx 

	movlps xmm4, [esi + ebx*4]
	movlps xmm3, [esi + ecx*4]
	movlps xmm7, [esi + edx*4]
	movhps xmm4, [esi + ebx*4 + 8]
	movhps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + edx*4 + 8]
	movaps xmm6, xmm3
	unpcklps xmm6, xmm7
	unpckhps xmm3, xmm7
	movaps xmm5, xmm4
	movaps xmm7, xmm4
	shufps xmm4, xmm6, 0x40
	shufps xmm5, xmm6, 0xE4
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 0x48
	shufps xmm7, xmm3, 0xEC
	;# coulomb table ready, in xmm4-xmm7
	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 

	xorps  xmm3, xmm3
	
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [esp + nb314nf_qqHH]
	movhps  xmm3, [esp + nb314nf_qqMH]
	shufps  xmm3, xmm3, 193 ;# constant 11000001
		
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point xmm5 contains vcoul 
	addps  xmm5, [esp + nb314nf_vctot]
	movaps [esp + nb314nf_vctot], xmm5	

	;# do table for i M - j water interaction 
	movaps xmm0, [esp + nb314nf_rinvMM]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb314nf_rsqMM]	;# xmm0=rinv, xmm1=r 
	mulps  xmm1, [esp + nb314nf_tsc]
	
	movhlps xmm2, xmm1	
	cvttps2pi mm6, xmm1
	cvttps2pi mm7, xmm2 	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2   	;# xmm2=eps2 
	pslld   mm6, 2
	pslld   mm7, 2
 
	psrlq mm6, 32
	movd ebx, mm6
	movd ecx, mm7
	psrlq mm7, 32
	movd edx, mm7		;# table indices in ebx,ecx,edx 

	movlps xmm4, [esi + ebx*4]
	movlps xmm3, [esi + ecx*4]
	movlps xmm7, [esi + edx*4]
	movhps xmm4, [esi + ebx*4 + 8]
	movhps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + edx*4 + 8]
	movaps xmm6, xmm3
	unpcklps xmm6, xmm7
	unpckhps xmm3, xmm7
	movaps xmm5, xmm4
	movaps xmm7, xmm4
	shufps xmm4, xmm6, 0x40
	shufps xmm5, xmm6, 0xE4
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 0x48
	shufps xmm7, xmm3, 0xEC
	;# coulomb table ready, in xmm4-xmm7
	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 

	xorps  xmm3, xmm3
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [esp + nb314nf_qqMH]
	movhps  xmm3, [esp + nb314nf_qqMM]
	shufps  xmm3, xmm3, 193 ;# constant 11000001
		
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point xmm5 contains vcoul 
	addps  xmm5, [esp + nb314nf_vctot]
	movaps [esp + nb314nf_vctot], xmm5	

	dec dword ptr [esp + nb314nf_innerk]
	jz    .nb314nf_updateouterdata
	jmp   .nb314nf_single_loop
.nb314nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb314nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb314nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb314nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb314nf_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [esp + nb314nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb314nf_Vvdw]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb314nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb314nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb314nf_n], esi
        jmp .nb314nf_outer
.nb314nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb314nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb314nf_end
        ;# non-zero, do one more workunit
        jmp   .nb314nf_threadloop
.nb314nf_end:
	emms

	mov eax, [esp + nb314nf_nouter]
	mov ebx, [esp + nb314nf_ninner]
	mov ecx, [ebp + nb314nf_outeriter]
	mov edx, [ebp + nb314nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb314nf_salign]
	add esp, eax
	add esp, 920
	pop edi
	pop esi
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret
