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

	
.globl nb_kernel214_ia32_sse
.globl _nb_kernel214_ia32_sse
nb_kernel214_ia32_sse:	
_nb_kernel214_ia32_sse:	
.equiv          nb214_p_nri,            8
.equiv          nb214_iinr,             12
.equiv          nb214_jindex,           16
.equiv          nb214_jjnr,             20
.equiv          nb214_shift,            24
.equiv          nb214_shiftvec,         28
.equiv          nb214_fshift,           32
.equiv          nb214_gid,              36
.equiv          nb214_pos,              40
.equiv          nb214_faction,          44
.equiv          nb214_charge,           48
.equiv          nb214_p_facel,          52
.equiv          nb214_argkrf,           56
.equiv          nb214_argcrf,           60
.equiv          nb214_Vc,               64
.equiv          nb214_type,             68
.equiv          nb214_p_ntype,          72
.equiv          nb214_vdwparam,         76
.equiv          nb214_Vvdw,             80
.equiv          nb214_p_tabscale,       84
.equiv          nb214_VFtab,            88
.equiv          nb214_invsqrta,         92
.equiv          nb214_dvda,             96
.equiv          nb214_p_gbtabscale,     100
.equiv          nb214_GBtab,            104
.equiv          nb214_p_nthreads,       108
.equiv          nb214_count,            112
.equiv          nb214_mtx,              116
.equiv          nb214_outeriter,        120
.equiv          nb214_inneriter,        124
.equiv          nb214_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb214_ixO,              0
.equiv          nb214_iyO,              16
.equiv          nb214_izO,              32
.equiv          nb214_ixH1,             48
.equiv          nb214_iyH1,             64
.equiv          nb214_izH1,             80
.equiv          nb214_ixH2,             96
.equiv          nb214_iyH2,             112
.equiv          nb214_izH2,             128
.equiv          nb214_ixM,              144
.equiv          nb214_iyM,              160
.equiv          nb214_izM,              176
.equiv          nb214_jxO,              192
.equiv          nb214_jyO,              208
.equiv          nb214_jzO,              224
.equiv          nb214_jxH1,             240
.equiv          nb214_jyH1,             256
.equiv          nb214_jzH1,             272
.equiv          nb214_jxH2,             288
.equiv          nb214_jyH2,             304
.equiv          nb214_jzH2,             320
.equiv          nb214_jxM,              336
.equiv          nb214_jyM,              352
.equiv          nb214_jzM,              368
.equiv          nb214_dxOO,             384
.equiv          nb214_dyOO,             400
.equiv          nb214_dzOO,             416
.equiv          nb214_dxH1H1,           432
.equiv          nb214_dyH1H1,           448
.equiv          nb214_dzH1H1,           464
.equiv          nb214_dxH1H2,           480
.equiv          nb214_dyH1H2,           496
.equiv          nb214_dzH1H2,           512
.equiv          nb214_dxH1M,            528
.equiv          nb214_dyH1M,            544
.equiv          nb214_dzH1M,            560
.equiv          nb214_dxH2H1,           576
.equiv          nb214_dyH2H1,           592
.equiv          nb214_dzH2H1,           608
.equiv          nb214_dxH2H2,           624
.equiv          nb214_dyH2H2,           640
.equiv          nb214_dzH2H2,           656
.equiv          nb214_dxH2M,            672
.equiv          nb214_dyH2M,            688
.equiv          nb214_dzH2M,            704
.equiv          nb214_dxMH1,            720
.equiv          nb214_dyMH1,            736
.equiv          nb214_dzMH1,            752
.equiv          nb214_dxMH2,            768
.equiv          nb214_dyMH2,            784
.equiv          nb214_dzMH2,            800
.equiv          nb214_dxMM,             816
.equiv          nb214_dyMM,             832
.equiv          nb214_dzMM,             848
.equiv          nb214_qqMM,             864
.equiv          nb214_qqMH,             880
.equiv          nb214_qqHH,             896
.equiv          nb214_two,              912
.equiv          nb214_c6,               928
.equiv          nb214_c12,              944
.equiv          nb214_six,              960
.equiv          nb214_twelve,           976
.equiv          nb214_vctot,            992
.equiv          nb214_Vvdwtot,          1008
.equiv          nb214_fixO,             1024
.equiv          nb214_fiyO,             1040
.equiv          nb214_fizO,             1056
.equiv          nb214_fixH1,            1072
.equiv          nb214_fiyH1,            1088
.equiv          nb214_fizH1,            1104
.equiv          nb214_fixH2,            1120
.equiv          nb214_fiyH2,            1136
.equiv          nb214_fizH2,            1152
.equiv          nb214_fixM,             1168
.equiv          nb214_fiyM,             1184
.equiv          nb214_fizM,             1200
.equiv          nb214_fjxO,             1216
.equiv          nb214_fjyO,             1232
.equiv          nb214_fjzO,             1248
.equiv          nb214_fjxH1,            1264
.equiv          nb214_fjyH1,            1280
.equiv          nb214_fjzH1,            1296
.equiv          nb214_fjxH2,            1312
.equiv          nb214_fjyH2,            1328
.equiv          nb214_fjzH2,            1344
.equiv          nb214_fjxM,             1360
.equiv          nb214_fjyM,             1376
.equiv          nb214_fjzM,             1392
.equiv          nb214_half,             1408
.equiv          nb214_three,            1424
.equiv          nb214_rsqOO,            1440
.equiv          nb214_rsqH1H1,          1456
.equiv          nb214_rsqH1H2,          1472
.equiv          nb214_rsqH1M,           1488
.equiv          nb214_rsqH2H1,          1504
.equiv          nb214_rsqH2H2,          1520
.equiv          nb214_rsqH2M,           1536
.equiv          nb214_rsqMH1,           1552
.equiv          nb214_rsqMH2,           1568
.equiv          nb214_rsqMM,            1584
.equiv          nb214_rinvsqOO,         1600
.equiv          nb214_rinvH1H1,         1616
.equiv          nb214_rinvH1H2,         1632
.equiv          nb214_rinvH1M,          1648
.equiv          nb214_rinvH2H1,         1664
.equiv          nb214_rinvH2H2,         1680
.equiv          nb214_rinvH2M,          1696
.equiv          nb214_rinvMH1,          1712
.equiv          nb214_rinvMH2,          1728
.equiv          nb214_rinvMM,           1744
.equiv          nb214_fstmp,            1760
.equiv          nb214_krf,              1776
.equiv          nb214_crf,              1792
.equiv          nb214_is3,              1808
.equiv          nb214_ii3,              1812
.equiv          nb214_innerjjnr,        1816
.equiv          nb214_innerk,           1820
.equiv          nb214_n,                1824
.equiv          nb214_nn1,              1828
.equiv          nb214_nri,              1832
.equiv          nb214_nouter,           1836
.equiv          nb214_ninner,           1840
.equiv          nb214_salign,           1844
	push ebp
	mov ebp,esp	
	push eax
	push ebx
	push ecx
	push edx
	push esi
	push edi
	sub esp, 1848		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb214_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb214_p_nri]
	mov ecx, [ecx]
	mov [esp + nb214_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb214_nouter], eax
	mov [esp + nb214_ninner], eax


	mov esi, [ebp + nb214_argkrf]
	mov edi, [ebp + nb214_argcrf]
	movss xmm5, [esi]
	movss xmm6, [edi]
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	movaps [esp + nb214_krf], xmm5
	movaps [esp + nb214_crf], xmm6
	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb214_half], eax
	movss xmm1, [esp + nb214_half]
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
	movaps [esp + nb214_half],  xmm1
	movaps [esp + nb214_two],  xmm2
	movaps [esp + nb214_three],  xmm3
	movaps [esp + nb214_six],  xmm4
	movaps [esp + nb214_twelve],  xmm5

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb214_iinr]   ;# ecx = pointer into iinr[]
	mov   ebx, [ecx]		;# ebx =ii 

	mov   edx, [ebp + nb214_charge]
	movss xmm5, [edx + ebx*4 + 4]	
	movss xmm3, [edx + ebx*4 + 12]	
	movss xmm4, xmm3	
	mov esi, [ebp + nb214_p_facel]
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
	movaps [esp + nb214_qqMM], xmm3
	movaps [esp + nb214_qqMH], xmm4
	movaps [esp + nb214_qqHH], xmm5
		
	xorps xmm0, xmm0
	mov   edx, [ebp + nb214_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov edi, [ebp + nb214_p_ntype]
	imul  ecx, [edi]  ;# ecx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   eax, [ebp + nb214_vdwparam]
	movlps xmm0, [eax + edx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0x55
	movaps [esp + nb214_c6], xmm0
	movaps [esp + nb214_c12], xmm1

.nb214_threadloop:
        mov   esi, [ebp + nb214_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb214_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb214_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb214_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb214_n], eax
        mov [esp + nb214_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb214_outerstart
        jmp .nb214_end
	
.nb214_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb214_nouter]
	mov [esp + nb214_nouter], ebx

.nb214_outer:
	mov   eax, [ebp + nb214_shift]  	;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]	;# ebx=3*is 
	mov   [esp + nb214_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb214_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb214_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]		;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb214_pos]	;# eax = base of pos[]  
	mov   [esp + nb214_ii3], ebx	
	
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
	movaps [esp + nb214_ixO], xmm3
	movaps [esp + nb214_iyO], xmm4
	movaps [esp + nb214_izO], xmm5
	movaps [esp + nb214_ixH1], xmm6
	movaps [esp + nb214_iyH1], xmm7

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
	movaps [esp + nb214_izH1], xmm6
	movaps [esp + nb214_ixH2], xmm0
	movaps [esp + nb214_iyH2], xmm1
	movaps [esp + nb214_izH2], xmm2
	movaps [esp + nb214_ixM], xmm3
	movaps [esp + nb214_iyM], xmm4
	movaps [esp + nb214_izM], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb214_vctot], xmm4
	movaps [esp + nb214_Vvdwtot], xmm4
	movaps [esp + nb214_fixO], xmm4
	movaps [esp + nb214_fiyO], xmm4
	movaps [esp + nb214_fizO], xmm4
	movaps [esp + nb214_fixH1], xmm4
	movaps [esp + nb214_fiyH1], xmm4
	movaps [esp + nb214_fizH1], xmm4
	movaps [esp + nb214_fixH2], xmm4
	movaps [esp + nb214_fiyH2], xmm4
	movaps [esp + nb214_fizH2], xmm4
	movaps [esp + nb214_fixM], xmm4
	movaps [esp + nb214_fiyM], xmm4
	movaps [esp + nb214_fizM], xmm4
	
	mov   eax, [ebp + nb214_jindex]
	mov   ecx, [eax + esi*4]	 	;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   esi, [ebp + nb214_pos]
	mov   edi, [ebp + nb214_faction]	
	mov   eax, [ebp + nb214_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb214_innerjjnr], eax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb214_ninner]
	mov   [esp + nb214_ninner], ecx
	add   edx, 0
	mov   [esp + nb214_innerk], edx	;# number of innerloop atoms 
	jge   .nb214_unroll_loop
	jmp   .nb214_single_check
.nb214_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb214_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	mov   ecx, [edx + 8]
	mov   edx, [edx + 12]     	;# eax-edx=jnr1-4 
	
	add dword ptr [esp + nb214_innerjjnr],  16 ;# advance pointer (unroll 4) 

	mov esi, [ebp + nb214_pos]   	;# base of pos[] 

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
	movaps [esp + nb214_jxO], xmm0
	movaps [esp + nb214_jyO], xmm1
	movaps [esp + nb214_jzO], xmm2
	movaps [esp + nb214_jxH1], xmm3

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
	movaps [esp + nb214_jyH1], xmm0
	movaps [esp + nb214_jzH1], xmm1
	movaps [esp + nb214_jxH2], xmm2
	movaps [esp + nb214_jyH2], xmm3

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
	movaps [esp + nb214_jzH2], xmm0
	movaps [esp + nb214_jxM], xmm1
	movaps [esp + nb214_jyM], xmm2
	movaps [esp + nb214_jzM], xmm3
	
	;# start calculating pairwise distances
	movaps xmm0, [esp + nb214_ixO]
	movaps xmm1, [esp + nb214_iyO]
	movaps xmm2, [esp + nb214_izO]
	movaps xmm3, [esp + nb214_ixH1]
	movaps xmm4, [esp + nb214_iyH1]
	movaps xmm5, [esp + nb214_izH1]
	subps  xmm0, [esp + nb214_jxO]
	subps  xmm1, [esp + nb214_jyO]
	subps  xmm2, [esp + nb214_jzO]
	subps  xmm3, [esp + nb214_jxH1]
	subps  xmm4, [esp + nb214_jyH1]
	subps  xmm5, [esp + nb214_jzH1]
	movaps [esp + nb214_dxOO], xmm0
	movaps [esp + nb214_dyOO], xmm1
	movaps [esp + nb214_dzOO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb214_dxH1H1], xmm3
	movaps [esp + nb214_dyH1H1], xmm4
	movaps [esp + nb214_dzH1H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb214_rsqOO], xmm0
	movaps [esp + nb214_rsqH1H1], xmm3

	movaps xmm0, [esp + nb214_ixH1]
	movaps xmm1, [esp + nb214_iyH1]
	movaps xmm2, [esp + nb214_izH1]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [esp + nb214_jxH2]
	subps  xmm1, [esp + nb214_jyH2]
	subps  xmm2, [esp + nb214_jzH2]
	subps  xmm3, [esp + nb214_jxM]
	subps  xmm4, [esp + nb214_jyM]
	subps  xmm5, [esp + nb214_jzM]
	movaps [esp + nb214_dxH1H2], xmm0
	movaps [esp + nb214_dyH1H2], xmm1
	movaps [esp + nb214_dzH1H2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb214_dxH1M], xmm3
	movaps [esp + nb214_dyH1M], xmm4
	movaps [esp + nb214_dzH1M], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb214_rsqH1H2], xmm0
	movaps [esp + nb214_rsqH1M], xmm3

	movaps xmm0, [esp + nb214_ixH2]
	movaps xmm1, [esp + nb214_iyH2]
	movaps xmm2, [esp + nb214_izH2]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [esp + nb214_jxH1]
	subps  xmm1, [esp + nb214_jyH1]
	subps  xmm2, [esp + nb214_jzH1]
	subps  xmm3, [esp + nb214_jxH2]
	subps  xmm4, [esp + nb214_jyH2]
	subps  xmm5, [esp + nb214_jzH2]
	movaps [esp + nb214_dxH2H1], xmm0
	movaps [esp + nb214_dyH2H1], xmm1
	movaps [esp + nb214_dzH2H1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb214_dxH2H2], xmm3
	movaps [esp + nb214_dyH2H2], xmm4
	movaps [esp + nb214_dzH2H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb214_rsqH2H1], xmm0
	movaps [esp + nb214_rsqH2H2], xmm3

	movaps xmm0, [esp + nb214_ixH2]
	movaps xmm1, [esp + nb214_iyH2]
	movaps xmm2, [esp + nb214_izH2]
	movaps xmm3, [esp + nb214_ixM]
	movaps xmm4, [esp + nb214_iyM]
	movaps xmm5, [esp + nb214_izM]
	subps  xmm0, [esp + nb214_jxM]
	subps  xmm1, [esp + nb214_jyM]
	subps  xmm2, [esp + nb214_jzM]
	subps  xmm3, [esp + nb214_jxH1]
	subps  xmm4, [esp + nb214_jyH1]
	subps  xmm5, [esp + nb214_jzH1]
	movaps [esp + nb214_dxH2M], xmm0
	movaps [esp + nb214_dyH2M], xmm1
	movaps [esp + nb214_dzH2M], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb214_dxMH1], xmm3
	movaps [esp + nb214_dyMH1], xmm4
	movaps [esp + nb214_dzMH1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [esp + nb214_rsqH2M], xmm0
	movaps [esp + nb214_rsqMH1], xmm4

	movaps xmm0, [esp + nb214_ixM]
	movaps xmm1, [esp + nb214_iyM]
	movaps xmm2, [esp + nb214_izM]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [esp + nb214_jxH2]
	subps  xmm1, [esp + nb214_jyH2]
	subps  xmm2, [esp + nb214_jzH2]
	subps  xmm3, [esp + nb214_jxM]
	subps  xmm4, [esp + nb214_jyM]
	subps  xmm5, [esp + nb214_jzM]
	movaps [esp + nb214_dxMH2], xmm0
	movaps [esp + nb214_dyMH2], xmm1
	movaps [esp + nb214_dzMH2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb214_dxMM], xmm3
	movaps [esp + nb214_dyMM], xmm4
	movaps [esp + nb214_dzMM], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [esp + nb214_rsqMH2], xmm0
	movaps [esp + nb214_rsqMM], xmm4

	;# start by doing reciprocal for OO
	movaps  xmm7, [esp + nb214_rsqOO]
	rcpps   xmm2, xmm7
	movaps  xmm1, [esp + nb214_two]
	mulps   xmm7, xmm2
	subps   xmm1, xmm7
	mulps   xmm2, xmm1 ;# rinvsq 
	movaps [esp + nb214_rinvsqOO], xmm2
	
	;# next step is invsqrt - do two at a time.
	rsqrtps xmm1, [esp + nb214_rsqH1H1]
	rsqrtps xmm5, [esp + nb214_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb214_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb214_rsqH1H1]
	mulps   xmm5, [esp + nb214_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb214_half] ;# rinvH1H1 
	mulps   xmm7, [esp + nb214_half] ;# rinvH1H2 
	movaps  [esp + nb214_rinvH1H1], xmm3
	movaps  [esp + nb214_rinvH1H2], xmm7
			
	rsqrtps xmm1, [esp + nb214_rsqH1M]
	rsqrtps xmm5, [esp + nb214_rsqH2H1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb214_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb214_rsqH1M]
	mulps   xmm5, [esp + nb214_rsqH2H1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb214_half] 
	mulps   xmm7, [esp + nb214_half]
	movaps  [esp + nb214_rinvH1M], xmm3
	movaps  [esp + nb214_rinvH2H1], xmm7
			
	rsqrtps xmm1, [esp + nb214_rsqH2H2]
	rsqrtps xmm5, [esp + nb214_rsqH2M]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb214_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb214_rsqH2H2]
	mulps   xmm5, [esp + nb214_rsqH2M]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb214_half] 
	mulps   xmm7, [esp + nb214_half]
	movaps  [esp + nb214_rinvH2H2], xmm3
	movaps  [esp + nb214_rinvH2M], xmm7
	
	rsqrtps xmm1, [esp + nb214_rsqMH1]
	rsqrtps xmm5, [esp + nb214_rsqMH2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb214_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb214_rsqMH1]
	mulps   xmm5, [esp + nb214_rsqMH2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb214_half] 
	mulps   xmm7, [esp + nb214_half]
	movaps  [esp + nb214_rinvMH1], xmm3
	movaps  [esp + nb214_rinvMH2], xmm7
        		
	rsqrtps xmm1, [esp + nb214_rsqMM]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + nb214_three]
	mulps   xmm1, [esp + nb214_rsqMM]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [esp + nb214_half] 
	movaps  [esp + nb214_rinvMM], xmm3
	
	;# start with OO LJ interaction
	movaps xmm0, [esp + nb214_rinvsqOO]
	movaps xmm1, xmm0
	mulps  xmm1, xmm1	;# rinv4
	mulps  xmm1, xmm0	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm1, [esp + nb214_c6]
	mulps  xmm2, [esp + nb214_c12]
	movaps xmm4, xmm2
	subps  xmm4, xmm1
	addps  xmm4, [esp + nb214_Vvdwtot]
	mulps  xmm1, [esp + nb214_six]
	mulps  xmm2, [esp + nb214_twelve]
	movaps [esp + nb214_Vvdwtot], xmm4
	subps  xmm2, xmm1
	mulps  xmm0, xmm2 	;# fscal 
	movaps xmm1, xmm0
	movaps xmm2, xmm0		

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + nb214_dxOO]
	mulps xmm1, [esp + nb214_dyOO]
	mulps xmm2, [esp + nb214_dzOO]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb214_fixO]
	addps xmm1, [esp + nb214_fiyO]
	addps xmm2, [esp + nb214_fizO]
	movaps [esp + nb214_fjxO], xmm3
	movaps [esp + nb214_fjyO], xmm4
	movaps [esp + nb214_fjzO], xmm5
	movaps [esp + nb214_fixO], xmm0
	movaps [esp + nb214_fiyO], xmm1
	movaps [esp + nb214_fizO], xmm2

	;# Coulomb interactions 
	;# start with H1-H1 interaction 
	movaps xmm0, [esp + nb214_rinvH1H1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb214_krf]
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	mulps  xmm5, [esp + nb214_rsqH1H1] ;# xmm5=krsq 
	movaps xmm6, xmm5
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [esp + nb214_crf]
	mulps  xmm6, [esp + nb214_qqHH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	mulps xmm5, [esp + nb214_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb214_qqHH] ;# xmm7 = coul part of fscal 
	
	addps  xmm6, [esp + nb214_vctot] ;# local vctot summation variable 
	mulps  xmm0, xmm7
	
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + nb214_dxH1H1]
	mulps xmm1, [esp + nb214_dyH1H1]
	mulps xmm2, [esp + nb214_dzH1H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb214_fixH1]
	addps xmm1, [esp + nb214_fiyH1]
	addps xmm2, [esp + nb214_fizH1]
	movaps [esp + nb214_fjxH1], xmm3
	movaps [esp + nb214_fjyH1], xmm4
	movaps [esp + nb214_fjzH1], xmm5
	movaps [esp + nb214_fixH1], xmm0
	movaps [esp + nb214_fiyH1], xmm1
	movaps [esp + nb214_fizH1], xmm2

	;# H1-H2 interaction 
	movaps xmm0, [esp + nb214_rinvH1H2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb214_krf]
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb214_rsqH1H2] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb214_crf]
	mulps  xmm0, xmm0
	mulps  xmm4, [esp + nb214_qqHH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb214_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb214_qqHH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH1  
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + nb214_dxH1H2]
	mulps xmm1, [esp + nb214_dyH1H2]
	mulps xmm2, [esp + nb214_dzH1H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb214_fixH1]
	addps xmm1, [esp + nb214_fiyH1]
	addps xmm2, [esp + nb214_fizH1]
	movaps [esp + nb214_fjxH2], xmm3
	movaps [esp + nb214_fjyH2], xmm4
	movaps [esp + nb214_fjzH2], xmm5
	movaps [esp + nb214_fixH1], xmm0
	movaps [esp + nb214_fiyH1], xmm1
	movaps [esp + nb214_fizH1], xmm2

	;# H1-M interaction  
	movaps xmm0, [esp + nb214_rinvH1M]
	movaps xmm7, xmm0	;# xmm7=Rinv 
	movaps xmm5, [esp + nb214_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb214_rsqH1M] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb214_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + nb214_qqMH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb214_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb214_qqMH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + nb214_dxH1M]
	mulps xmm1, [esp + nb214_dyH1M]
	mulps xmm2, [esp + nb214_dzH1M]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb214_fixH1]
	addps xmm1, [esp + nb214_fiyH1]
	addps xmm2, [esp + nb214_fizH1]
	movaps [esp + nb214_fjxM], xmm3
	movaps [esp + nb214_fjyM], xmm4
	movaps [esp + nb214_fjzM], xmm5
	movaps [esp + nb214_fixH1], xmm0
	movaps [esp + nb214_fiyH1], xmm1
	movaps [esp + nb214_fizH1], xmm2

	;# H2-H1 interaction 
	movaps xmm0, [esp + nb214_rinvH2H1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb214_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb214_rsqH2H1] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb214_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + nb214_qqHH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb214_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb214_qqHH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [esp + nb214_fjxH1]
	movaps xmm4, [esp + nb214_fjyH1]
	movaps xmm5, [esp + nb214_fjzH1]
	mulps xmm0, [esp + nb214_dxH2H1]
	mulps xmm1, [esp + nb214_dyH2H1]
	mulps xmm2, [esp + nb214_dzH2H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb214_fixH2]
	addps xmm1, [esp + nb214_fiyH2]
	addps xmm2, [esp + nb214_fizH2]
	movaps [esp + nb214_fjxH1], xmm3
	movaps [esp + nb214_fjyH1], xmm4
	movaps [esp + nb214_fjzH1], xmm5
	movaps [esp + nb214_fixH2], xmm0
	movaps [esp + nb214_fiyH2], xmm1
	movaps [esp + nb214_fizH2], xmm2

	;# H2-H2 interaction 
	movaps xmm0, [esp + nb214_rinvH2H2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb214_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb214_rsqH2H2] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb214_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + nb214_qqHH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb214_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb214_qqHH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [esp + nb214_fjxH2]
	movaps xmm4, [esp + nb214_fjyH2]
	movaps xmm5, [esp + nb214_fjzH2]
	mulps xmm0, [esp + nb214_dxH2H2]
	mulps xmm1, [esp + nb214_dyH2H2]
	mulps xmm2, [esp + nb214_dzH2H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb214_fixH2]
	addps xmm1, [esp + nb214_fiyH2]
	addps xmm2, [esp + nb214_fizH2]
	movaps [esp + nb214_fjxH2], xmm3
	movaps [esp + nb214_fjyH2], xmm4
	movaps [esp + nb214_fjzH2], xmm5
	movaps [esp + nb214_fixH2], xmm0
	movaps [esp + nb214_fiyH2], xmm1
	movaps [esp + nb214_fizH2], xmm2

	;# H2-M interaction 
	movaps xmm0, [esp + nb214_rinvH2M]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb214_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb214_rsqH2M] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb214_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + nb214_qqMH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb214_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb214_qqMH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	movaps xmm3, [esp + nb214_fjxM]
	movaps xmm4, [esp + nb214_fjyM]
	movaps xmm5, [esp + nb214_fjzM]
	mulps xmm0, [esp + nb214_dxH2M]
	mulps xmm1, [esp + nb214_dyH2M]
	mulps xmm2, [esp + nb214_dzH2M]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb214_fixH2]
	addps xmm1, [esp + nb214_fiyH2]
	addps xmm2, [esp + nb214_fizH2]
	movaps [esp + nb214_fjxM], xmm3
	movaps [esp + nb214_fjyM], xmm4
	movaps [esp + nb214_fjzM], xmm5
	movaps [esp + nb214_fixH2], xmm0
	movaps [esp + nb214_fiyH2], xmm1
	movaps [esp + nb214_fizH2], xmm2

	;# M-H1 interaction 
	movaps xmm0, [esp + nb214_rinvMH1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb214_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb214_rsqMH1] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb214_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + nb214_qqMH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb214_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb214_qqMH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [esp + nb214_fjxH1]
	movaps xmm4, [esp + nb214_fjyH1]
	movaps xmm5, [esp + nb214_fjzH1]
	mulps xmm0, [esp + nb214_dxMH1]
	mulps xmm1, [esp + nb214_dyMH1]
	mulps xmm2, [esp + nb214_dzMH1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb214_fixM]
	addps xmm1, [esp + nb214_fiyM]
	addps xmm2, [esp + nb214_fizM]
	movaps [esp + nb214_fjxH1], xmm3
	movaps [esp + nb214_fjyH1], xmm4
	movaps [esp + nb214_fjzH1], xmm5
	movaps [esp + nb214_fixM], xmm0
	movaps [esp + nb214_fiyM], xmm1
	movaps [esp + nb214_fizM], xmm2

	;# M-H2 interaction 
	movaps xmm0, [esp + nb214_rinvMH2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb214_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb214_rsqMH2] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb214_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + nb214_qqMH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb214_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb214_qqMH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [esp + nb214_fjxH2]
	movaps xmm4, [esp + nb214_fjyH2]
	movaps xmm5, [esp + nb214_fjzH2]
	mulps xmm0, [esp + nb214_dxMH2]
	mulps xmm1, [esp + nb214_dyMH2]
	mulps xmm2, [esp + nb214_dzMH2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb214_fixM]
	addps xmm1, [esp + nb214_fiyM]
	addps xmm2, [esp + nb214_fizM]
	movaps [esp + nb214_fjxH2], xmm3
	movaps [esp + nb214_fjyH2], xmm4
	movaps [esp + nb214_fjzH2], xmm5
	movaps [esp + nb214_fixM], xmm0
	movaps [esp + nb214_fiyM], xmm1
	movaps [esp + nb214_fizM], xmm2

	;# M-M interaction 
	movaps xmm0, [esp + nb214_rinvMM]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb214_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb214_rsqMM] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb214_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + nb214_qqMM] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb214_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb214_qqMM] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm1, xmm0
	movaps [esp + nb214_vctot], xmm6
	movaps xmm2, xmm0
	
	movaps xmm3, [esp + nb214_fjxM]
	movaps xmm4, [esp + nb214_fjyM]
	movaps xmm5, [esp + nb214_fjzM]
	mulps xmm0, [esp + nb214_dxMM]
	mulps xmm1, [esp + nb214_dyMM]
	mulps xmm2, [esp + nb214_dzMM]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb214_fixM]
	addps xmm1, [esp + nb214_fiyM]
	addps xmm2, [esp + nb214_fizM]
	movaps [esp + nb214_fjxM], xmm3
	movaps [esp + nb214_fjyM], xmm4
	movaps [esp + nb214_fjzM], xmm5
	movaps [esp + nb214_fixM], xmm0
	movaps [esp + nb214_fiyM], xmm1
	movaps [esp + nb214_fizM], xmm2

	mov edi, [ebp + nb214_faction]
	;# update j forces 
	;# 4 j waters with four atoms each.
	;# step 1 : transpose fjxO, fjyO, fjzO, fjxH1
	movaps xmm0, [esp + nb214_fjxO]
	movaps xmm1, [esp + nb214_fjyO]
	movaps xmm2, [esp + nb214_fjzO]
	movaps xmm3, [esp + nb214_fjxH1]
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
	movaps xmm0, [esp + nb214_fjyH1]
	movaps xmm1, [esp + nb214_fjzH1]
	movaps xmm2, [esp + nb214_fjxH2]
	movaps xmm3, [esp + nb214_fjyH2]
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
	movaps xmm0, [esp + nb214_fjzH2]
	movaps xmm1, [esp + nb214_fjxM]
	movaps xmm2, [esp + nb214_fjyM]
	movaps xmm3, [esp + nb214_fjzM]
	
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
	sub dword ptr [esp + nb214_innerk],  4
	jl    .nb214_single_check
	jmp   .nb214_unroll_loop
.nb214_single_check:
	add dword ptr [esp + nb214_innerk],  4
	jnz   .nb214_single_loop
	jmp   .nb214_updateouterdata
.nb214_single_loop:
	mov   edx, [esp + nb214_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb214_innerjjnr],  4	

	mov esi, [ebp + nb214_pos]
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
	movaps [esp + nb214_jxO], xmm6
	movaps [esp + nb214_jyO], xmm3
	movaps [esp + nb214_jzO], xmm1

	;# do O and M in parallel
	movaps xmm0, [esp + nb214_ixO]
	movaps xmm1, [esp + nb214_iyO]
	movaps xmm2, [esp + nb214_izO]
	movaps xmm3, [esp + nb214_ixM]
	movaps xmm4, [esp + nb214_iyM]
	movaps xmm5, [esp + nb214_izM]
	subps  xmm0, [esp + nb214_jxO]
	subps  xmm1, [esp + nb214_jyO]
	subps  xmm2, [esp + nb214_jzO]
	subps  xmm3, [esp + nb214_jxO]
	subps  xmm4, [esp + nb214_jyO]
	subps  xmm5, [esp + nb214_jzO]
	
	movaps [esp + nb214_dxOO], xmm0
	movaps [esp + nb214_dyOO], xmm1
	movaps [esp + nb214_dzOO], xmm2
	movaps [esp + nb214_dxMM], xmm3
	movaps [esp + nb214_dyMM], xmm4
	movaps [esp + nb214_dzMM], xmm5
	
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
	;# Save M data 
	movaps [esp + nb214_rsqMM], xmm4
	
	;# do 1/x for O and 1/sqrt(x) for M
	rcpss  xmm1, xmm0
	rsqrtps xmm5, xmm4
	movss  xmm2, [esp + nb214_two]
	movaps  xmm6, xmm5	
	mulss  xmm0, xmm1
	mulps   xmm5, xmm5
	subss  xmm2, xmm0
	movaps  xmm7, [esp + nb214_three]
	mulss  xmm2, xmm1 	;# constant 1/r2
	
	mulps   xmm5, xmm4
	movss  xmm0, xmm2
	subps   xmm7, xmm5
	mulss  xmm2, xmm2
	mulps   xmm7, xmm6
	mulss  xmm2, xmm0 	;# constant 1/r6
	mulps   xmm7, [esp + nb214_half] ;# rinv iH1 - j water 
	movss  xmm1, xmm2
	movaps [esp + nb214_rinvMM], xmm7

	mulss  xmm2, xmm2 	;# constant 1/r12
	mulss  xmm1, [esp + nb214_c6]
	mulss  xmm2, [esp + nb214_c12]
	movss  xmm3, xmm2
	subss  xmm3, xmm1
	addss  xmm3, [esp + nb214_Vvdwtot]
	movss  [esp + nb214_Vvdwtot], xmm3
	mulss  xmm1, [esp + nb214_six]
	mulss  xmm2, [esp + nb214_twelve]
	subss  xmm2, xmm1
	mulss  xmm0, xmm2 	;# fscal
	movss  xmm1, xmm0
	movss  xmm2, xmm0
	mulss  xmm0, [esp + nb214_dxOO]
	mulss  xmm1, [esp + nb214_dyOO]
	mulss  xmm2, [esp + nb214_dzOO]
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2
	movaps  [esp + nb214_fjxO], xmm3
	movaps  [esp + nb214_fjyO], xmm4
	movaps  [esp + nb214_fjzO], xmm5
	addss   xmm0, [esp + nb214_fixO]
	addss   xmm1, [esp + nb214_fiyO]
	addss   xmm2, [esp + nb214_fizO]
	movss  [esp + nb214_fixO], xmm0
	movss  [esp + nb214_fiyO], xmm1
	movss  [esp + nb214_fizO], xmm2

	;# do  M coulomb interaction
	movaps xmm0, [esp + nb214_rinvMM]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb214_krf]
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [esp + nb214_qqMH]
	movhps  xmm3, [esp + nb214_qqMM]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	mulps  xmm5, [esp + nb214_rsqMM] ;# xmm5=krsq 
	movaps xmm6, xmm5
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [esp + nb214_crf]
	mulps  xmm6, xmm3 ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	mulps xmm5, [esp + nb214_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, xmm3 ;# xmm7 = coul part of fscal 
	
	addps  xmm6, [esp + nb214_vctot] 
	movaps [esp + nb214_vctot], xmm6
	mulps  xmm0, xmm7

	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [esp + nb214_dxMM]
	mulps   xmm1, [esp + nb214_dyMM]
	mulps   xmm2, [esp + nb214_dzMM]
	;# update forces M - j water 
	movaps  xmm3, [esp + nb214_fjxO]
	movaps  xmm4, [esp + nb214_fjyO]
	movaps  xmm5, [esp + nb214_fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + nb214_fjxO], xmm3
	movaps  [esp + nb214_fjyO], xmm4
	movaps  [esp + nb214_fjzO], xmm5
	addps   xmm0, [esp + nb214_fixM]
	addps   xmm1, [esp + nb214_fiyM]
	addps   xmm2, [esp + nb214_fizM]
	movaps  [esp + nb214_fixM], xmm0
	movaps  [esp + nb214_fiyM], xmm1
	movaps  [esp + nb214_fizM], xmm2	
	
	;# i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [esp + nb214_ixH1]
	movaps  xmm1, [esp + nb214_iyH1]
	movaps  xmm2, [esp + nb214_izH1]	
	movaps  xmm3, [esp + nb214_ixH2] 
	movaps  xmm4, [esp + nb214_iyH2] 
	movaps  xmm5, [esp + nb214_izH2] 
	subps   xmm0, [esp + nb214_jxO]
	subps   xmm1, [esp + nb214_jyO]
	subps   xmm2, [esp + nb214_jzO]
	subps   xmm3, [esp + nb214_jxO]
	subps   xmm4, [esp + nb214_jyO]
	subps   xmm5, [esp + nb214_jzO]
	movaps [esp + nb214_dxH1H1], xmm0
	movaps [esp + nb214_dyH1H1], xmm1
	movaps [esp + nb214_dzH1H1], xmm2
	movaps [esp + nb214_dxH2H2], xmm3
	movaps [esp + nb214_dyH2H2], xmm4
	movaps [esp + nb214_dzH2H2], xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;# have rsqH1 in xmm0 
	addps xmm4, xmm5	;# have rsqH2 in xmm4 
	movaps  [esp + nb214_rsqH1H1], xmm0
	movaps  [esp + nb214_rsqH2H2], xmm4
		
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb214_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb214_half] ;# rinvH1H1
	mulps   xmm7, [esp + nb214_half] ;# rinvH2H2
	movaps  [esp + nb214_rinvH1H1], xmm3
	movaps  [esp + nb214_rinvH2H2], xmm7
	
	;# Do H1 coulomb interaction
	movaps xmm0, [esp + nb214_rinvH1H1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb214_krf]
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [esp + nb214_qqHH]
	movhps  xmm3, [esp + nb214_qqMH]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	mulps  xmm5, [esp + nb214_rsqH1H1] ;# xmm5=krsq 
	movaps xmm6, xmm5
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [esp + nb214_crf]
	mulps  xmm6, xmm3 ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	mulps xmm5, [esp + nb214_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, xmm3 ;# xmm7 = coul part of fscal 
	
	addps  xmm6, [esp + nb214_vctot] 
	movaps [esp + nb214_vctot], xmm6

	mulps  xmm0, xmm7

	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [esp + nb214_dxH1H1]
	mulps   xmm1, [esp + nb214_dyH1H1]
	mulps   xmm2, [esp + nb214_dzH1H1]
	;# update forces H1 - j water 
	movaps  xmm3, [esp + nb214_fjxO]
	movaps  xmm4, [esp + nb214_fjyO]
	movaps  xmm5, [esp + nb214_fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + nb214_fjxO], xmm3
	movaps  [esp + nb214_fjyO], xmm4
	movaps  [esp + nb214_fjzO], xmm5
	addps   xmm0, [esp + nb214_fixH1]
	addps   xmm1, [esp + nb214_fiyH1]
	addps   xmm2, [esp + nb214_fizH1]
	movaps  [esp + nb214_fixH1], xmm0
	movaps  [esp + nb214_fiyH1], xmm1
	movaps  [esp + nb214_fizH1], xmm2	

	;# H2 Coulomb
	movaps xmm0, [esp + nb214_rinvH2H2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb214_krf]
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [esp + nb214_qqHH]
	movhps  xmm3, [esp + nb214_qqMH]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	mulps  xmm5, [esp + nb214_rsqH2H2] ;# xmm5=krsq 
	movaps xmm6, xmm5
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [esp + nb214_crf]
	mulps  xmm6, xmm3 ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	mulps xmm5, [esp + nb214_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, xmm3 ;# xmm7 = coul part of fscal 
	
	addps  xmm6, [esp + nb214_vctot] ;# local vctot summation variable
	movaps [esp + nb214_vctot], xmm6
	mulps  xmm0, xmm7

	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [esp + nb214_dxH2H2]
	mulps   xmm1, [esp + nb214_dyH2H2]
	mulps   xmm2, [esp + nb214_dzH2H2]
	;# update forces H2 - j water 
	movaps  xmm3, [esp + nb214_fjxO]
	movaps  xmm4, [esp + nb214_fjyO]
	movaps  xmm5, [esp + nb214_fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + nb214_fjxO], xmm3
	movaps  [esp + nb214_fjyO], xmm4
	movaps  [esp + nb214_fjzO], xmm5
	addps   xmm0, [esp + nb214_fixH2]
	addps   xmm1, [esp + nb214_fiyH2]
	addps   xmm2, [esp + nb214_fizH2]
	movaps  [esp + nb214_fixH2], xmm0
	movaps  [esp + nb214_fiyH2], xmm1
	movaps  [esp + nb214_fizH2], xmm2			
	
	mov     esi, [ebp + nb214_faction]
	;# update j water forces from local variables.
	;# transpose back first
	movaps  xmm0, [esp + nb214_fjxO] ;# Ox H1x H2x Mx 
	movaps  xmm1, [esp + nb214_fjyO] ;# Oy H1y H2y My
	movaps  xmm2, [esp + nb214_fjzO] ;# Oz H1z H2z Mz
	 
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
	
	dec dword ptr [esp + nb214_innerk]
	jz    .nb214_updateouterdata
	jmp   .nb214_single_loop
.nb214_updateouterdata:
	mov   ecx, [esp + nb214_ii3]
	mov   edi, [ebp + nb214_faction]
	mov   esi, [ebp + nb214_fshift]
	mov   edx, [esp + nb214_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb214_fixO]
	movaps xmm1, [esp + nb214_fiyO] 
	movaps xmm2, [esp + nb214_fizO]

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
	movaps xmm0, [esp + nb214_fixH1]
	movaps xmm1, [esp + nb214_fiyH1]
	movaps xmm2, [esp + nb214_fizH1]

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
	movaps xmm0, [esp + nb214_fixH2]
	movaps xmm1, [esp + nb214_fiyH2]
	movaps xmm2, [esp + nb214_fizH2]

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
	movaps xmm0, [esp + nb214_fixM]
	movaps xmm1, [esp + nb214_fiyM]
	movaps xmm2, [esp + nb214_fizM]

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
	mov esi, [esp + nb214_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb214_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb214_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb214_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [esp + nb214_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb214_Vvdw]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb214_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb214_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb214_n], esi
        jmp .nb214_outer
.nb214_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb214_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb214_end
        ;# non-zero, do one more workunit
        jmp   .nb214_threadloop
.nb214_end:
	emms

	mov eax, [esp + nb214_nouter]
	mov ebx, [esp + nb214_ninner]
	mov ecx, [ebp + nb214_outeriter]
	mov edx, [ebp + nb214_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb214_salign]
	add esp, eax
	add esp, 1848
	pop edi
	pop esi
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret


	
.globl nb_kernel214nf_ia32_sse
.globl _nb_kernel214nf_ia32_sse
nb_kernel214nf_ia32_sse:	
_nb_kernel214nf_ia32_sse:	
.equiv          nb214nf_p_nri,          8
.equiv          nb214nf_iinr,           12
.equiv          nb214nf_jindex,         16
.equiv          nb214nf_jjnr,           20
.equiv          nb214nf_shift,          24
.equiv          nb214nf_shiftvec,       28
.equiv          nb214nf_fshift,         32
.equiv          nb214nf_gid,            36
.equiv          nb214nf_pos,            40
.equiv          nb214nf_faction,        44
.equiv          nb214nf_charge,         48
.equiv          nb214nf_p_facel,        52
.equiv          nb214nf_argkrf,         56
.equiv          nb214nf_argcrf,         60
.equiv          nb214nf_Vc,             64
.equiv          nb214nf_type,           68
.equiv          nb214nf_p_ntype,        72
.equiv          nb214nf_vdwparam,       76
.equiv          nb214nf_Vvdw,           80
.equiv          nb214nf_p_tabscale,     84
.equiv          nb214nf_VFtab,          88
.equiv          nb214nf_invsqrta,       92
.equiv          nb214nf_dvda,           96
.equiv          nb214nf_p_gbtabscale,   100
.equiv          nb214nf_GBtab,          104
.equiv          nb214nf_p_nthreads,     108
.equiv          nb214nf_count,          112
.equiv          nb214nf_mtx,            116
.equiv          nb214nf_outeriter,      120
.equiv          nb214nf_inneriter,      124
.equiv          nb214nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb214nf_ixO,            0
.equiv          nb214nf_iyO,            16
.equiv          nb214nf_izO,            32
.equiv          nb214nf_ixH1,           48
.equiv          nb214nf_iyH1,           64
.equiv          nb214nf_izH1,           80
.equiv          nb214nf_ixH2,           96
.equiv          nb214nf_iyH2,           112
.equiv          nb214nf_izH2,           128
.equiv          nb214nf_ixM,            144
.equiv          nb214nf_iyM,            160
.equiv          nb214nf_izM,            176
.equiv          nb214nf_jxO,            192
.equiv          nb214nf_jyO,            208
.equiv          nb214nf_jzO,            224
.equiv          nb214nf_jxH1,           240
.equiv          nb214nf_jyH1,           256
.equiv          nb214nf_jzH1,           272
.equiv          nb214nf_jxH2,           288
.equiv          nb214nf_jyH2,           304
.equiv          nb214nf_jzH2,           320
.equiv          nb214nf_jxM,            336
.equiv          nb214nf_jyM,            352
.equiv          nb214nf_jzM,            368
.equiv          nb214nf_qqMM,           384
.equiv          nb214nf_qqMH,           400
.equiv          nb214nf_qqHH,           416
.equiv          nb214nf_two,            432
.equiv          nb214nf_c6,             448
.equiv          nb214nf_c12,            464
.equiv          nb214nf_vctot,          480
.equiv          nb214nf_Vvdwtot,        496
.equiv          nb214nf_half,           512
.equiv          nb214nf_three,          528
.equiv          nb214nf_rsqOO,          544
.equiv          nb214nf_rsqH1H1,        560
.equiv          nb214nf_rsqH1H2,        576
.equiv          nb214nf_rsqH1M,         592
.equiv          nb214nf_rsqH2H1,        608
.equiv          nb214nf_rsqH2H2,        624
.equiv          nb214nf_rsqH2M,         640
.equiv          nb214nf_rsqMH1,         656
.equiv          nb214nf_rsqMH2,         672
.equiv          nb214nf_rsqMM,          688
.equiv          nb214nf_rinvsqOO,       704
.equiv          nb214nf_rinvH1H1,       720
.equiv          nb214nf_rinvH1H2,       736
.equiv          nb214nf_rinvH1M,        752
.equiv          nb214nf_rinvH2H1,       768
.equiv          nb214nf_rinvH2H2,       784
.equiv          nb214nf_rinvH2M,        800
.equiv          nb214nf_rinvMH1,        816
.equiv          nb214nf_rinvMH2,        832
.equiv          nb214nf_rinvMM,         848
.equiv          nb214nf_krf,            864
.equiv          nb214nf_crf,            880
.equiv          nb214nf_is3,            896
.equiv          nb214nf_ii3,            900
.equiv          nb214nf_innerjjnr,      904
.equiv          nb214nf_innerk,         908
.equiv          nb214nf_n,              912
.equiv          nb214nf_nn1,            916
.equiv          nb214nf_nri,            920
.equiv          nb214nf_nouter,         924
.equiv          nb214nf_ninner,         928
.equiv          nb214nf_salign,         932
	push ebp
	mov ebp,esp	
	push eax
	push ebx
	push ecx
	push edx
	push esi
	push edi
	sub esp, 936		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb214nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb214nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb214nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb214nf_nouter], eax
	mov [esp + nb214nf_ninner], eax


	mov esi, [ebp + nb214nf_argkrf]
	mov edi, [ebp + nb214nf_argcrf]
	movss xmm5, [esi]
	movss xmm6, [edi]
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	movaps [esp + nb214nf_krf], xmm5
	movaps [esp + nb214nf_crf], xmm6
	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb214nf_half], eax
	movss xmm1, [esp + nb214nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb214nf_half],  xmm1
	movaps [esp + nb214nf_two],  xmm2
	movaps [esp + nb214nf_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb214nf_iinr]   ;# ecx = pointer into iinr[]
	mov   ebx, [ecx]		;# ebx =ii 

	mov   edx, [ebp + nb214nf_charge]
	movss xmm5, [edx + ebx*4 + 4]	
	movss xmm3, [edx + ebx*4 + 12]	
	movss xmm4, xmm3	
	mov esi, [ebp + nb214nf_p_facel]
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
	movaps [esp + nb214nf_qqMM], xmm3
	movaps [esp + nb214nf_qqMH], xmm4
	movaps [esp + nb214nf_qqHH], xmm5
		
	xorps xmm0, xmm0
	mov   edx, [ebp + nb214nf_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov edi, [ebp + nb214nf_p_ntype]
	imul  ecx, [edi]  ;# ecx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   eax, [ebp + nb214nf_vdwparam]
	movlps xmm0, [eax + edx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0x55
	movaps [esp + nb214nf_c6], xmm0
	movaps [esp + nb214nf_c12], xmm1

.nb214nf_threadloop:
        mov   esi, [ebp + nb214nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb214nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb214nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb214nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb214nf_n], eax
        mov [esp + nb214nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb214nf_outerstart
        jmp .nb214nf_end
	
.nb214nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb214nf_nouter]
	mov [esp + nb214nf_nouter], ebx

.nb214nf_outer:
	mov   eax, [ebp + nb214nf_shift]  	;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]	;# ebx=3*is 
	mov   [esp + nb214nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb214nf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb214nf_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]		;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb214nf_pos]	;# eax = base of pos[]  
	mov   [esp + nb214nf_ii3], ebx	
	
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
	movaps [esp + nb214nf_ixO], xmm3
	movaps [esp + nb214nf_iyO], xmm4
	movaps [esp + nb214nf_izO], xmm5
	movaps [esp + nb214nf_ixH1], xmm6
	movaps [esp + nb214nf_iyH1], xmm7

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
	movaps [esp + nb214nf_izH1], xmm6
	movaps [esp + nb214nf_ixH2], xmm0
	movaps [esp + nb214nf_iyH2], xmm1
	movaps [esp + nb214nf_izH2], xmm2
	movaps [esp + nb214nf_ixM], xmm3
	movaps [esp + nb214nf_iyM], xmm4
	movaps [esp + nb214nf_izM], xmm5

	;# clear vctot 
	xorps xmm4, xmm4
	movaps [esp + nb214nf_vctot], xmm4
	movaps [esp + nb214nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb214nf_jindex]
	mov   ecx, [eax + esi*4]	 	;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   esi, [ebp + nb214nf_pos]
	mov   edi, [ebp + nb214nf_faction]	
	mov   eax, [ebp + nb214nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb214nf_innerjjnr], eax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb214nf_ninner]
	mov   [esp + nb214nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb214nf_innerk], edx	;# number of innerloop atoms 
	jge   .nb214nf_unroll_loop
	jmp   .nb214nf_single_check
.nb214nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb214nf_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	mov   ecx, [edx + 8]
	mov   edx, [edx + 12]     	;# eax-edx=jnr1-4 
	
	add dword ptr [esp + nb214nf_innerjjnr],  16 ;# advance pointer (unroll 4) 

	mov esi, [ebp + nb214nf_pos]   	;# base of pos[] 

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
	movaps [esp + nb214nf_jxO], xmm0
	movaps [esp + nb214nf_jyO], xmm1
	movaps [esp + nb214nf_jzO], xmm2
	movaps [esp + nb214nf_jxH1], xmm3

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
	movaps [esp + nb214nf_jyH1], xmm0
	movaps [esp + nb214nf_jzH1], xmm1
	movaps [esp + nb214nf_jxH2], xmm2
	movaps [esp + nb214nf_jyH2], xmm3

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
	movaps [esp + nb214nf_jzH2], xmm0
	movaps [esp + nb214nf_jxM], xmm1
	movaps [esp + nb214nf_jyM], xmm2
	movaps [esp + nb214nf_jzM], xmm3
	
	;# start calculating pairwise distances
	movaps xmm0, [esp + nb214nf_ixO]
	movaps xmm1, [esp + nb214nf_iyO]
	movaps xmm2, [esp + nb214nf_izO]
	movaps xmm3, [esp + nb214nf_ixH1]
	movaps xmm4, [esp + nb214nf_iyH1]
	movaps xmm5, [esp + nb214nf_izH1]
	subps  xmm0, [esp + nb214nf_jxO]
	subps  xmm1, [esp + nb214nf_jyO]
	subps  xmm2, [esp + nb214nf_jzO]
	subps  xmm3, [esp + nb214nf_jxH1]
	subps  xmm4, [esp + nb214nf_jyH1]
	subps  xmm5, [esp + nb214nf_jzH1]
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
	movaps [esp + nb214nf_rsqOO], xmm0
	movaps [esp + nb214nf_rsqH1H1], xmm3

	movaps xmm0, [esp + nb214nf_ixH1]
	movaps xmm1, [esp + nb214nf_iyH1]
	movaps xmm2, [esp + nb214nf_izH1]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [esp + nb214nf_jxH2]
	subps  xmm1, [esp + nb214nf_jyH2]
	subps  xmm2, [esp + nb214nf_jzH2]
	subps  xmm3, [esp + nb214nf_jxM]
	subps  xmm4, [esp + nb214nf_jyM]
	subps  xmm5, [esp + nb214nf_jzM]
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
	movaps [esp + nb214nf_rsqH1H2], xmm0
	movaps [esp + nb214nf_rsqH1M], xmm3

	movaps xmm0, [esp + nb214nf_ixH2]
	movaps xmm1, [esp + nb214nf_iyH2]
	movaps xmm2, [esp + nb214nf_izH2]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [esp + nb214nf_jxH1]
	subps  xmm1, [esp + nb214nf_jyH1]
	subps  xmm2, [esp + nb214nf_jzH1]
	subps  xmm3, [esp + nb214nf_jxH2]
	subps  xmm4, [esp + nb214nf_jyH2]
	subps  xmm5, [esp + nb214nf_jzH2]
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
	movaps [esp + nb214nf_rsqH2H1], xmm0
	movaps [esp + nb214nf_rsqH2H2], xmm3

	movaps xmm0, [esp + nb214nf_ixH2]
	movaps xmm1, [esp + nb214nf_iyH2]
	movaps xmm2, [esp + nb214nf_izH2]
	movaps xmm3, [esp + nb214nf_ixM]
	movaps xmm4, [esp + nb214nf_iyM]
	movaps xmm5, [esp + nb214nf_izM]
	subps  xmm0, [esp + nb214nf_jxM]
	subps  xmm1, [esp + nb214nf_jyM]
	subps  xmm2, [esp + nb214nf_jzM]
	subps  xmm3, [esp + nb214nf_jxH1]
	subps  xmm4, [esp + nb214nf_jyH1]
	subps  xmm5, [esp + nb214nf_jzH1]
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
	movaps [esp + nb214nf_rsqH2M], xmm0
	movaps [esp + nb214nf_rsqMH1], xmm4

	movaps xmm0, [esp + nb214nf_ixM]
	movaps xmm1, [esp + nb214nf_iyM]
	movaps xmm2, [esp + nb214nf_izM]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [esp + nb214nf_jxH2]
	subps  xmm1, [esp + nb214nf_jyH2]
	subps  xmm2, [esp + nb214nf_jzH2]
	subps  xmm3, [esp + nb214nf_jxM]
	subps  xmm4, [esp + nb214nf_jyM]
	subps  xmm5, [esp + nb214nf_jzM]
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
	movaps [esp + nb214nf_rsqMH2], xmm0
	movaps [esp + nb214nf_rsqMM], xmm4

	;# start by doing reciprocal for OO
	movaps  xmm7, [esp + nb214nf_rsqOO]
	rcpps   xmm2, xmm7
	movaps  xmm1, [esp + nb214nf_two]
	mulps   xmm7, xmm2
	subps   xmm1, xmm7
	mulps   xmm2, xmm1 ;# rinvsq 
	movaps [esp + nb214nf_rinvsqOO], xmm2
	
	;# next step is invsqrt - do two at a time.
	rsqrtps xmm1, [esp + nb214nf_rsqH1H1]
	rsqrtps xmm5, [esp + nb214nf_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb214nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb214nf_rsqH1H1]
	mulps   xmm5, [esp + nb214nf_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb214nf_half] ;# rinvH1H1 
	mulps   xmm7, [esp + nb214nf_half] ;# rinvH1H2 
	movaps  [esp + nb214nf_rinvH1H1], xmm3
	movaps  [esp + nb214nf_rinvH1H2], xmm7
			
	rsqrtps xmm1, [esp + nb214nf_rsqH1M]
	rsqrtps xmm5, [esp + nb214nf_rsqH2H1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb214nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb214nf_rsqH1M]
	mulps   xmm5, [esp + nb214nf_rsqH2H1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb214nf_half] 
	mulps   xmm7, [esp + nb214nf_half]
	movaps  [esp + nb214nf_rinvH1M], xmm3
	movaps  [esp + nb214nf_rinvH2H1], xmm7
			
	rsqrtps xmm1, [esp + nb214nf_rsqH2H2]
	rsqrtps xmm5, [esp + nb214nf_rsqH2M]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb214nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb214nf_rsqH2H2]
	mulps   xmm5, [esp + nb214nf_rsqH2M]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb214nf_half] 
	mulps   xmm7, [esp + nb214nf_half]
	movaps  [esp + nb214nf_rinvH2H2], xmm3
	movaps  [esp + nb214nf_rinvH2M], xmm7
	
	rsqrtps xmm1, [esp + nb214nf_rsqMH1]
	rsqrtps xmm5, [esp + nb214nf_rsqMH2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb214nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb214nf_rsqMH1]
	mulps   xmm5, [esp + nb214nf_rsqMH2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb214nf_half] 
	mulps   xmm7, [esp + nb214nf_half]
	movaps  [esp + nb214nf_rinvMH1], xmm3
	movaps  [esp + nb214nf_rinvMH2], xmm7
        		
	rsqrtps xmm1, [esp + nb214nf_rsqMM]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + nb214nf_three]
	mulps   xmm1, [esp + nb214nf_rsqMM]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [esp + nb214nf_half] 
	movaps  [esp + nb214nf_rinvMM], xmm3
	
	;# start with OO LJ interaction
	movaps xmm0, [esp + nb214nf_rinvsqOO]
	movaps xmm1, xmm0
	mulps  xmm1, xmm1	;# rinv4
	mulps  xmm1, xmm0	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm1, [esp + nb214nf_c6]
	mulps  xmm2, [esp + nb214nf_c12]
	movaps xmm4, xmm2
	subps  xmm4, xmm1
	addps  xmm4, [esp + nb214nf_Vvdwtot]
	movaps [esp + nb214nf_Vvdwtot], xmm4

	;# Coulomb interactions 
	;# add all H-H rsq in xmm2, H-M rsq xmm4
	;# H-H rinv in xmm0, H-M in xmm1
	movaps xmm0, [esp + nb214nf_rinvH1H1]
	movaps xmm1, [esp + nb214nf_rinvH1M]
	movaps xmm2, [esp + nb214nf_rsqH1H1]
	movaps xmm4, [esp + nb214nf_rsqH1M]
	addps  xmm0, [esp + nb214nf_rinvH1H2]
	addps  xmm1, [esp + nb214nf_rinvH2M]
	addps  xmm2, [esp + nb214nf_rsqH1H2]
	addps  xmm4, [esp + nb214nf_rsqH2M]
	addps  xmm0, [esp + nb214nf_rinvH2H1]
	addps  xmm1, [esp + nb214nf_rinvMH1]
	addps  xmm2, [esp + nb214nf_rsqH2H1]
	addps  xmm4, [esp + nb214nf_rsqMH1]
	addps  xmm0, [esp + nb214nf_rinvH2H2]
	addps  xmm1, [esp + nb214nf_rinvMH2]
	addps  xmm2, [esp + nb214nf_rsqH2H2]
	addps  xmm4, [esp + nb214nf_rsqMH2]
	movaps xmm5, [esp + nb214nf_krf]
	movaps xmm6, [esp + nb214nf_crf]

	;# calc 4*crf in xmm7
	movaps xmm7, xmm6
	addps  xmm7, xmm7
	addps  xmm7, xmm7
	mulps  xmm2, xmm5 ;# H-H krsq
	mulps  xmm4, xmm5 ;# H-M krsq
	addps  xmm0, xmm2 ;# H-H rinv+krsq
	addps  xmm1, xmm4 ;# H-M rinv+krsq
	subps  xmm0, xmm7 ;# H-H rinv+krsq-crf
	subps  xmm1, xmm7 ;# H-M rinv+krsq-crf
	mulps  xmm0, [esp + nb214nf_qqHH]
	mulps  xmm1, [esp + nb214nf_qqMH]
	addps  xmm0, xmm1 
	addps  xmm0, [esp + nb214nf_vctot]
	;# M-M interaction
	movaps xmm4, [esp + nb214nf_rinvMM]
	mulps  xmm5, [esp + nb214nf_rsqMM] ;# krsq
	addps  xmm5, xmm4                  ;# rinv+krsq
	subps xmm5, [esp + nb214nf_crf] ;# xmm5=rinv+ krsq-crf 
	mulps xmm5, [esp + nb214nf_qqMM]
	addps xmm5, xmm0
	movaps [esp + nb214nf_vctot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [esp + nb214nf_innerk],  4
	jl    .nb214nf_single_check
	jmp   .nb214nf_unroll_loop
.nb214nf_single_check:
	add dword ptr [esp + nb214nf_innerk],  4
	jnz   .nb214nf_single_loop
	jmp   .nb214nf_updateouterdata
.nb214nf_single_loop:
	mov   edx, [esp + nb214nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb214nf_innerjjnr],  4	

	mov esi, [ebp + nb214nf_pos]
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
	movaps [esp + nb214nf_jxO], xmm6
	movaps [esp + nb214nf_jyO], xmm3
	movaps [esp + nb214nf_jzO], xmm1

	;# do O and M in parallel
	movaps xmm0, [esp + nb214nf_ixO]
	movaps xmm1, [esp + nb214nf_iyO]
	movaps xmm2, [esp + nb214nf_izO]
	movaps xmm3, [esp + nb214nf_ixM]
	movaps xmm4, [esp + nb214nf_iyM]
	movaps xmm5, [esp + nb214nf_izM]
	subps  xmm0, [esp + nb214nf_jxO]
	subps  xmm1, [esp + nb214nf_jyO]
	subps  xmm2, [esp + nb214nf_jzO]
	subps  xmm3, [esp + nb214nf_jxO]
	subps  xmm4, [esp + nb214nf_jyO]
	subps  xmm5, [esp + nb214nf_jzO]
		
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
	;# Save M data 
	movaps [esp + nb214nf_rsqMM], xmm4
	
	;# do 1/x for O and 1/sqrt(x) for M
	rcpss  xmm1, xmm0
	rsqrtps xmm5, xmm4
	movss  xmm2, [esp + nb214nf_two]
	movaps  xmm6, xmm5	
	mulss  xmm0, xmm1
	mulps   xmm5, xmm5
	subss  xmm2, xmm0
	movaps  xmm7, [esp + nb214nf_three]
	mulss  xmm2, xmm1 	;# constant 1/r2
	
	mulps   xmm5, xmm4
	movss  xmm0, xmm2
	subps   xmm7, xmm5
	mulss  xmm2, xmm2
	mulps   xmm7, xmm6
	mulss  xmm2, xmm0 	;# constant 1/r6
	mulps   xmm7, [esp + nb214nf_half] ;# rinv iH1 - j water 
	movss  xmm1, xmm2
	movaps [esp + nb214nf_rinvMM], xmm7

	mulss  xmm2, xmm2 	;# constant 1/r12
	mulss  xmm1, [esp + nb214nf_c6]
	mulss  xmm2, [esp + nb214nf_c12]
	movss  xmm3, xmm2
	subss  xmm3, xmm1
	addss  xmm3, [esp + nb214nf_Vvdwtot]
	movss  [esp + nb214nf_Vvdwtot], xmm3

	;# do  M coulomb interaction
	movaps xmm7, [esp + nb214nf_rinvMM]
	movaps xmm6, [esp + nb214nf_krf]

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [esp + nb214nf_qqMH]
	movhps  xmm3, [esp + nb214nf_qqMM]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	mulps  xmm6, [esp + nb214nf_rsqMM] ;# xmm6=krsq 
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [esp + nb214nf_crf]
	mulps  xmm6, xmm3 ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addps  xmm6, [esp + nb214nf_vctot] 
	movaps [esp + nb214nf_vctot], xmm6
	
	;# i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [esp + nb214nf_ixH1]
	movaps  xmm1, [esp + nb214nf_iyH1]
	movaps  xmm2, [esp + nb214nf_izH1]	
	movaps  xmm3, [esp + nb214nf_ixH2] 
	movaps  xmm4, [esp + nb214nf_iyH2] 
	movaps  xmm5, [esp + nb214nf_izH2] 
	subps   xmm0, [esp + nb214nf_jxO]
	subps   xmm1, [esp + nb214nf_jyO]
	subps   xmm2, [esp + nb214nf_jzO]
	subps   xmm3, [esp + nb214nf_jxO]
	subps   xmm4, [esp + nb214nf_jyO]
	subps   xmm5, [esp + nb214nf_jzO]
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;# have rsqH1 in xmm0 
	addps xmm4, xmm5	;# have rsqH2 in xmm4 
	movaps  [esp + nb214nf_rsqH1H1], xmm0
	movaps  [esp + nb214nf_rsqH2H2], xmm4
		
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb214nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb214nf_half] ;# rinvH1H1
	mulps   xmm7, [esp + nb214nf_half] ;# rinvH2H2
	movaps  [esp + nb214nf_rinvH1H1], xmm3
	movaps  [esp + nb214nf_rinvH2H2], xmm7
	
	;# Do H1 coulomb interaction
	movaps xmm7, [esp + nb214nf_rinvH1H1]
	movaps xmm6, [esp + nb214nf_krf]

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [esp + nb214nf_qqHH]
	movhps  xmm3, [esp + nb214nf_qqMH]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	mulps  xmm6, [esp + nb214nf_rsqH1H1] ;# xmm6=krsq 
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [esp + nb214nf_crf]
	mulps  xmm6, xmm3 ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addps  xmm6, [esp + nb214nf_vctot] 
	movaps [esp + nb214nf_vctot], xmm6

	;# H2 Coulomb
	movaps xmm7, [esp + nb214nf_rinvH2H2]
	movaps xmm6, [esp + nb214nf_krf]

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [esp + nb214nf_qqHH]
	movhps  xmm3, [esp + nb214nf_qqMH]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	mulps  xmm6, [esp + nb214nf_rsqH2H2] ;# xmm6=krsq 
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [esp + nb214nf_crf]
	mulps  xmm6, xmm3 ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addps  xmm6, [esp + nb214nf_vctot] ;# local vctot summation variable
	movaps [esp + nb214nf_vctot], xmm6
	dec dword ptr [esp + nb214nf_innerk]
	jz    .nb214nf_updateouterdata
	jmp   .nb214nf_single_loop
.nb214nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb214nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb214nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb214nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb214nf_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [esp + nb214nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb214nf_Vvdw]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb214nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb214nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb214nf_n], esi
        jmp .nb214nf_outer
.nb214nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb214nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb214nf_end
        ;# non-zero, do one more workunit
        jmp   .nb214nf_threadloop
.nb214nf_end:
	emms

	mov eax, [esp + nb214nf_nouter]
	mov ebx, [esp + nb214nf_ninner]
	mov ecx, [ebp + nb214nf_outeriter]
	mov edx, [ebp + nb214nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb214nf_salign]
	add esp, eax
	add esp, 936
	pop edi
	pop esi
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret
