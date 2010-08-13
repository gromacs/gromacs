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


	
.globl nb_kernel204_ia32_sse
.globl _nb_kernel204_ia32_sse
nb_kernel204_ia32_sse:	
_nb_kernel204_ia32_sse:	
.equiv          nb204_p_nri,            8
.equiv          nb204_iinr,             12
.equiv          nb204_jindex,           16
.equiv          nb204_jjnr,             20
.equiv          nb204_shift,            24
.equiv          nb204_shiftvec,         28
.equiv          nb204_fshift,           32
.equiv          nb204_gid,              36
.equiv          nb204_pos,              40
.equiv          nb204_faction,          44
.equiv          nb204_charge,           48
.equiv          nb204_p_facel,          52
.equiv          nb204_argkrf,           56
.equiv          nb204_argcrf,           60
.equiv          nb204_Vc,               64
.equiv          nb204_type,             68
.equiv          nb204_p_ntype,          72
.equiv          nb204_vdwparam,         76
.equiv          nb204_Vvdw,             80
.equiv          nb204_p_tabscale,       84
.equiv          nb204_VFtab,            88
.equiv          nb204_invsqrta,         92
.equiv          nb204_dvda,             96
.equiv          nb204_p_gbtabscale,     100
.equiv          nb204_GBtab,            104
.equiv          nb204_p_nthreads,       108
.equiv          nb204_count,            112
.equiv          nb204_mtx,              116
.equiv          nb204_outeriter,        120
.equiv          nb204_inneriter,        124
.equiv          nb204_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb204_ixH1,             0
.equiv          nb204_iyH1,             16
.equiv          nb204_izH1,             32
.equiv          nb204_ixH2,             48
.equiv          nb204_iyH2,             64
.equiv          nb204_izH2,             80
.equiv          nb204_ixM,              96
.equiv          nb204_iyM,              112
.equiv          nb204_izM,              128
.equiv          nb204_jxH1,             144
.equiv          nb204_jyH1,             160
.equiv          nb204_jzH1,             176
.equiv          nb204_jxH2,             192
.equiv          nb204_jyH2,             208
.equiv          nb204_jzH2,             224
.equiv          nb204_jxM,              240
.equiv          nb204_jyM,              256
.equiv          nb204_jzM,              272
.equiv          nb204_dxH1H1,           288
.equiv          nb204_dyH1H1,           304
.equiv          nb204_dzH1H1,           320
.equiv          nb204_dxH1H2,           336
.equiv          nb204_dyH1H2,           352
.equiv          nb204_dzH1H2,           368
.equiv          nb204_dxH1M,            384
.equiv          nb204_dyH1M,            400
.equiv          nb204_dzH1M,            416
.equiv          nb204_dxH2H1,           432
.equiv          nb204_dyH2H1,           448
.equiv          nb204_dzH2H1,           464
.equiv          nb204_dxH2H2,           480
.equiv          nb204_dyH2H2,           496
.equiv          nb204_dzH2H2,           512
.equiv          nb204_dxH2M,            528
.equiv          nb204_dyH2M,            544
.equiv          nb204_dzH2M,            560
.equiv          nb204_dxMH1,            576
.equiv          nb204_dyMH1,            592
.equiv          nb204_dzMH1,            608
.equiv          nb204_dxMH2,            624
.equiv          nb204_dyMH2,            640
.equiv          nb204_dzMH2,            656
.equiv          nb204_dxMM,             672
.equiv          nb204_dyMM,             688
.equiv          nb204_dzMM,             704
.equiv          nb204_qqHH,             720
.equiv          nb204_qqMH,             736
.equiv          nb204_qqMM,             752
.equiv          nb204_vctot,            768
.equiv          nb204_fixH1,            784
.equiv          nb204_fiyH1,            800
.equiv          nb204_fizH1,            816
.equiv          nb204_fixH2,            832
.equiv          nb204_fiyH2,            848
.equiv          nb204_fizH2,            864
.equiv          nb204_fixM,             880
.equiv          nb204_fiyM,             896
.equiv          nb204_fizM,             912
.equiv          nb204_fjxH1,            928
.equiv          nb204_fjyH1,            944
.equiv          nb204_fjzH1,            960
.equiv          nb204_fjxH2,            976
.equiv          nb204_fjyH2,            992
.equiv          nb204_fjzH2,            1008
.equiv          nb204_fjxM,             1024
.equiv          nb204_fjyM,             1040
.equiv          nb204_fjzM,             1056
.equiv          nb204_fjzMb,            1060
.equiv          nb204_fjzMc,            1064
.equiv          nb204_fjzMd,            1068
.equiv          nb204_half,             1072
.equiv          nb204_three,            1088
.equiv          nb204_rsqH1H1,          1104
.equiv          nb204_rsqH1H2,          1120
.equiv          nb204_rsqH1M,           1136
.equiv          nb204_rsqH2H1,          1152
.equiv          nb204_rsqH2H2,          1168
.equiv          nb204_rsqH2M,           1184
.equiv          nb204_rsqMH1,           1200
.equiv          nb204_rsqMH2,           1216
.equiv          nb204_rsqMM,            1232
.equiv          nb204_rinvH1H1,         1248
.equiv          nb204_rinvH1H2,         1264
.equiv          nb204_rinvH1M,          1280
.equiv          nb204_rinvH2H1,         1296
.equiv          nb204_rinvH2H2,         1312
.equiv          nb204_rinvH2M,          1328
.equiv          nb204_rinvMH1,          1344
.equiv          nb204_rinvMH2,          1360
.equiv          nb204_rinvMM,           1376
.equiv          nb204_two,              1392
.equiv          nb204_krf,              1408
.equiv          nb204_crf,              1424
.equiv          nb204_is3,              1440
.equiv          nb204_ii3,              1444
.equiv          nb204_innerjjnr,        1448
.equiv          nb204_innerk,           1452
.equiv          nb204_n,                1456
.equiv          nb204_nn1,              1460
.equiv          nb204_nri,              1464
.equiv          nb204_nouter,           1468
.equiv          nb204_ninner,           1472
.equiv          nb204_salign,           1476
	push ebp
	mov ebp,esp	
	push eax
	push ebx
	push ecx
	push edx
	push esi
	push edi
	sub esp, 1480		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb204_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb204_p_nri]
	mov ecx, [ecx]
	mov [esp + nb204_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb204_nouter], eax
	mov [esp + nb204_ninner], eax


	mov esi, [ebp + nb204_argkrf]
	mov edi, [ebp + nb204_argcrf]
	movss xmm5, [esi]
	movss xmm6, [edi]
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	movaps [esp + nb204_krf], xmm5
	movaps [esp + nb204_crf], xmm6
	
	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb204_half], eax
	movss xmm1, [esp + nb204_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb204_half],  xmm1
	movaps [esp + nb204_two],  xmm2
	movaps [esp + nb204_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb204_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]		;# ebx =ii 

	mov   edx, [ebp + nb204_charge]
	movss xmm3, [edx + ebx*4 + 4]	
	movss xmm4, xmm3	
	movss xmm5, [edx + ebx*4 + 12]	
	mov esi, [ebp + nb204_p_facel]
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
	movaps [esp + nb204_qqHH], xmm3
	movaps [esp + nb204_qqMH], xmm4
	movaps [esp + nb204_qqMM], xmm5
	
.nb204_threadloop:
        mov   esi, [ebp + nb204_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb204_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb204_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb204_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb204_n], eax
        mov [esp + nb204_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb204_outerstart
        jmp .nb204_end
	
.nb204_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb204_nouter]
	mov [esp + nb204_nouter], ebx

.nb204_outer:
	mov   eax, [ebp + nb204_shift]  	;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]	;# ebx=3*is 
	mov   [esp + nb204_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb204_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb204_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]		;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb204_pos]	;# eax = base of pos[]  
	mov   [esp + nb204_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [eax + ebx*4 + 12]
	addss xmm4, [eax + ebx*4 + 16]
	addss xmm5, [eax + ebx*4 + 20]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb204_ixH1], xmm3
	movaps [esp + nb204_iyH1], xmm4
	movaps [esp + nb204_izH1], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 24]
	addss xmm1, [eax + ebx*4 + 28]
	addss xmm2, [eax + ebx*4 + 32]		
	addss xmm3, [eax + ebx*4 + 36]
	addss xmm4, [eax + ebx*4 + 40]
	addss xmm5, [eax + ebx*4 + 44]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb204_ixH2], xmm0
	movaps [esp + nb204_iyH2], xmm1
	movaps [esp + nb204_izH2], xmm2
	movaps [esp + nb204_ixM], xmm3
	movaps [esp + nb204_iyM], xmm4
	movaps [esp + nb204_izM], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb204_vctot], xmm4
	movaps [esp + nb204_fixH1], xmm4
	movaps [esp + nb204_fiyH1], xmm4
	movaps [esp + nb204_fizH1], xmm4
	movaps [esp + nb204_fixH2], xmm4
	movaps [esp + nb204_fiyH2], xmm4
	movaps [esp + nb204_fizH2], xmm4
	movaps [esp + nb204_fixM], xmm4
	movaps [esp + nb204_fiyM], xmm4
	movaps [esp + nb204_fizM], xmm4
	
	mov   eax, [ebp + nb204_jindex]
	mov   ecx, [eax + esi*4]	 	;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   esi, [ebp + nb204_pos]
	mov   edi, [ebp + nb204_faction]	
	mov   eax, [ebp + nb204_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb204_innerjjnr], eax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb204_ninner]
	mov   [esp + nb204_ninner], ecx
	add   edx, 0
	mov   [esp + nb204_innerk], edx	;# number of innerloop atoms 
	jge   .nb204_unroll_loop
	jmp   .nb204_single_check
.nb204_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb204_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	mov   ecx, [edx + 8]
	mov   edx, [edx + 12]     	;# eax-edx=jnr1-4 
	
	add dword ptr [esp + nb204_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb204_pos]   	;# base of pos[] 

	lea   eax, [eax + eax*2] 	;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2] 	;# replace jnr with j3 
	lea   edx, [edx + edx*2]	
	
	;# move j coordinates to local temp variables 
	movlps xmm2, [esi + eax*4 + 12]
	movlps xmm3, [esi + eax*4 + 24]
	movlps xmm4, [esi + eax*4 + 36]

	movlps xmm5, [esi + ebx*4 + 12]
	movlps xmm6, [esi + ebx*4 + 24]
	movlps xmm7, [esi + ebx*4 + 36]

	movhps xmm2, [esi + ecx*4 + 12]
	movhps xmm3, [esi + ecx*4 + 24]
	movhps xmm4, [esi + ecx*4 + 36]

	movhps xmm5, [esi + edx*4 + 12]
	movhps xmm6, [esi + edx*4 + 24]
	movhps xmm7, [esi + edx*4 + 36]
	
	movaps xmm0, xmm2
	movaps xmm1, xmm3
	unpcklps xmm0, xmm5
	unpcklps xmm1, xmm6
	unpckhps xmm2, xmm5
	unpckhps xmm3, xmm6
	movaps xmm5, xmm4
	movaps   xmm6, xmm0
	unpcklps xmm4, xmm7
	unpckhps xmm5, xmm7
	movaps   xmm7, xmm1
	movlhps  xmm0, xmm2
	movaps [esp + nb204_jxH1], xmm0
	movhlps  xmm2, xmm6
	movaps [esp + nb204_jyH1], xmm2
	movlhps  xmm1, xmm3
	movaps [esp + nb204_jxH2], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [esp + nb204_jyH2], xmm3
	movlhps  xmm4, xmm5
	movaps [esp + nb204_jxM], xmm4
	movhlps  xmm5, xmm6
	movaps [esp + nb204_jyM], xmm5

	movss  xmm0, [esi + eax*4 + 20]
	movss  xmm1, [esi + eax*4 + 32]
	movss  xmm2, [esi + eax*4 + 44]

	movss  xmm3, [esi + ecx*4 + 20]
	movss  xmm4, [esi + ecx*4 + 32]
	movss  xmm5, [esi + ecx*4 + 44]

	movhps xmm0, [esi + ebx*4 + 16]
	movhps xmm1, [esi + ebx*4 + 28]
	movhps xmm2, [esi + ebx*4 + 40]
	
	movhps xmm3, [esi + edx*4 + 16]
	movhps xmm4, [esi + edx*4 + 28]
	movhps xmm5, [esi + edx*4 + 40]
	
	shufps xmm0, xmm3, 204  ;# constant 11001100
	shufps xmm1, xmm4, 204  ;# constant 11001100
	shufps xmm2, xmm5, 204  ;# constant 11001100
	movaps [esp + nb204_jzH1],  xmm0
	movaps [esp + nb204_jzH2],  xmm1
	movaps [esp + nb204_jzM],  xmm2

	movaps xmm0, [esp + nb204_ixH1]
	movaps xmm1, [esp + nb204_iyH1]
	movaps xmm2, [esp + nb204_izH1]
	movaps xmm3, [esp + nb204_ixH1]
	movaps xmm4, [esp + nb204_iyH1]
	movaps xmm5, [esp + nb204_izH1]
	subps  xmm0, [esp + nb204_jxH1]
	subps  xmm1, [esp + nb204_jyH1]
	subps  xmm2, [esp + nb204_jzH1]
	subps  xmm3, [esp + nb204_jxH2]
	subps  xmm4, [esp + nb204_jyH2]
	subps  xmm5, [esp + nb204_jzH2]
	movaps [esp + nb204_dxH1H1], xmm0
	movaps [esp + nb204_dyH1H1], xmm1
	movaps [esp + nb204_dzH1H1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb204_dxH1H2], xmm3
	movaps [esp + nb204_dyH1H2], xmm4
	movaps [esp + nb204_dzH1H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb204_rsqH1H1], xmm0
	movaps [esp + nb204_rsqH1H2], xmm3

	movaps xmm0, [esp + nb204_ixH1]
	movaps xmm1, [esp + nb204_iyH1]
	movaps xmm2, [esp + nb204_izH1]
	movaps xmm3, [esp + nb204_ixH2]
	movaps xmm4, [esp + nb204_iyH2]
	movaps xmm5, [esp + nb204_izH2]
	subps  xmm0, [esp + nb204_jxM]
	subps  xmm1, [esp + nb204_jyM]
	subps  xmm2, [esp + nb204_jzM]
	subps  xmm3, [esp + nb204_jxH1]
	subps  xmm4, [esp + nb204_jyH1]
	subps  xmm5, [esp + nb204_jzH1]
	movaps [esp + nb204_dxH1M], xmm0
	movaps [esp + nb204_dyH1M], xmm1
	movaps [esp + nb204_dzH1M], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb204_dxH2H1], xmm3
	movaps [esp + nb204_dyH2H1], xmm4
	movaps [esp + nb204_dzH2H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb204_rsqH1M], xmm0
	movaps [esp + nb204_rsqH2H1], xmm3

	movaps xmm0, [esp + nb204_ixH2]
	movaps xmm1, [esp + nb204_iyH2]
	movaps xmm2, [esp + nb204_izH2]
	movaps xmm3, [esp + nb204_ixH2]
	movaps xmm4, [esp + nb204_iyH2]
	movaps xmm5, [esp + nb204_izH2]
	subps  xmm0, [esp + nb204_jxH2]
	subps  xmm1, [esp + nb204_jyH2]
	subps  xmm2, [esp + nb204_jzH2]
	subps  xmm3, [esp + nb204_jxM]
	subps  xmm4, [esp + nb204_jyM]
	subps  xmm5, [esp + nb204_jzM]
	movaps [esp + nb204_dxH2H2], xmm0
	movaps [esp + nb204_dyH2H2], xmm1
	movaps [esp + nb204_dzH2H2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb204_dxH2M], xmm3
	movaps [esp + nb204_dyH2M], xmm4
	movaps [esp + nb204_dzH2M], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb204_rsqH2H2], xmm0
	movaps [esp + nb204_rsqH2M], xmm3

	movaps xmm0, [esp + nb204_ixM]
	movaps xmm1, [esp + nb204_iyM]
	movaps xmm2, [esp + nb204_izM]
	movaps xmm3, [esp + nb204_ixM]
	movaps xmm4, [esp + nb204_iyM]
	movaps xmm5, [esp + nb204_izM]
	subps  xmm0, [esp + nb204_jxH1]
	subps  xmm1, [esp + nb204_jyH1]
	subps  xmm2, [esp + nb204_jzH1]
	subps  xmm3, [esp + nb204_jxH2]
	subps  xmm4, [esp + nb204_jyH2]
	subps  xmm5, [esp + nb204_jzH2]
	movaps [esp + nb204_dxMH1], xmm0
	movaps [esp + nb204_dyMH1], xmm1
	movaps [esp + nb204_dzMH1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb204_dxMH2], xmm3
	movaps [esp + nb204_dyMH2], xmm4
	movaps [esp + nb204_dzMH2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [esp + nb204_rsqMH1], xmm0
	movaps [esp + nb204_rsqMH2], xmm4

	movaps xmm0, [esp + nb204_ixM]
	movaps xmm1, [esp + nb204_iyM]
	movaps xmm2, [esp + nb204_izM]
	subps  xmm0, [esp + nb204_jxM]
	subps  xmm1, [esp + nb204_jyM]
	subps  xmm2, [esp + nb204_jzM]
	movaps [esp + nb204_dxMM], xmm0
	movaps [esp + nb204_dyMM], xmm1
	movaps [esp + nb204_dzMM], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [esp + nb204_rsqMM], xmm0
		
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb204_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb204_half] ;# rinvH2H2 
	mulps   xmm7, [esp + nb204_half] ;# rinvH2H1 
	movaps  [esp + nb204_rinvMM], xmm3
	movaps  [esp + nb204_rinvMH2], xmm7
	
	rsqrtps xmm1, [esp + nb204_rsqH1H1]
	rsqrtps xmm5, [esp + nb204_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb204_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb204_rsqH1H1]
	mulps   xmm5, [esp + nb204_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb204_half] 
	mulps   xmm7, [esp + nb204_half]
	movaps  [esp + nb204_rinvH1H1], xmm3
	movaps  [esp + nb204_rinvH1H2], xmm7
	
	rsqrtps xmm1, [esp + nb204_rsqH1M]
	rsqrtps xmm5, [esp + nb204_rsqH2H1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb204_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb204_rsqH1M]
	mulps   xmm5, [esp + nb204_rsqH2H1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb204_half] 
	mulps   xmm7, [esp + nb204_half]
	movaps  [esp + nb204_rinvH1M], xmm3
	movaps  [esp + nb204_rinvH2H1], xmm7
	
	rsqrtps xmm1, [esp + nb204_rsqH2H2]
	rsqrtps xmm5, [esp + nb204_rsqH2M]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb204_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb204_rsqH2H2]
	mulps   xmm5, [esp + nb204_rsqH2M]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb204_half] 
	mulps   xmm7, [esp + nb204_half]
	movaps  [esp + nb204_rinvH2H2], xmm3
	movaps  [esp + nb204_rinvH2M], xmm7
	
	rsqrtps xmm1, [esp + nb204_rsqMH1]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + nb204_three]
	mulps   xmm1, [esp + nb204_rsqMH1]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [esp + nb204_half] 
	movaps  [esp + nb204_rinvMH1], xmm3

	;# start with H1-H1 interaction 
	movaps xmm0, [esp + nb204_rinvH1H1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb204_krf]
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	mulps  xmm5, [esp + nb204_rsqH1H1] ;# xmm5=krsq 
	movaps xmm6, xmm5
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [esp + nb204_crf]
	mulps  xmm6, [esp + nb204_qqHH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	mulps xmm5, [esp + nb204_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb204_qqHH] ;# xmm7 = coul part of fscal 
	
	addps  xmm6, [esp + nb204_vctot] ;# local vctot summation variable 
	mulps  xmm0, xmm7
	
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + nb204_dxH1H1]
	mulps xmm1, [esp + nb204_dyH1H1]
	mulps xmm2, [esp + nb204_dzH1H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb204_fixH1]
	addps xmm1, [esp + nb204_fiyH1]
	addps xmm2, [esp + nb204_fizH1]
	movaps [esp + nb204_fjxH1], xmm3
	movaps [esp + nb204_fjyH1], xmm4
	movaps [esp + nb204_fjzH1], xmm5
	movaps [esp + nb204_fixH1], xmm0
	movaps [esp + nb204_fiyH1], xmm1
	movaps [esp + nb204_fizH1], xmm2

	;# H1-H2 interaction 
	movaps xmm0, [esp + nb204_rinvH1H2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb204_krf]
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb204_rsqH1H2] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb204_crf]
	mulps  xmm0, xmm0
	mulps  xmm4, [esp + nb204_qqHH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb204_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb204_qqHH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH1  
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + nb204_dxH1H2]
	mulps xmm1, [esp + nb204_dyH1H2]
	mulps xmm2, [esp + nb204_dzH1H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb204_fixH1]
	addps xmm1, [esp + nb204_fiyH1]
	addps xmm2, [esp + nb204_fizH1]
	movaps [esp + nb204_fjxH2], xmm3
	movaps [esp + nb204_fjyH2], xmm4
	movaps [esp + nb204_fjzH2], xmm5
	movaps [esp + nb204_fixH1], xmm0
	movaps [esp + nb204_fiyH1], xmm1
	movaps [esp + nb204_fizH1], xmm2

	;# H1-M interaction  
	movaps xmm0, [esp + nb204_rinvH1M]
	movaps xmm7, xmm0	;# xmm7=Rinv 
	movaps xmm5, [esp + nb204_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb204_rsqH1M] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb204_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + nb204_qqMH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb204_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb204_qqMH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + nb204_dxH1M]
	mulps xmm1, [esp + nb204_dyH1M]
	mulps xmm2, [esp + nb204_dzH1M]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb204_fixH1]
	addps xmm1, [esp + nb204_fiyH1]
	addps xmm2, [esp + nb204_fizH1]
	movaps [esp + nb204_fjxM], xmm3
	movaps [esp + nb204_fjyM], xmm4
	movaps [esp + nb204_fjzM], xmm5
	movaps [esp + nb204_fixH1], xmm0
	movaps [esp + nb204_fiyH1], xmm1
	movaps [esp + nb204_fizH1], xmm2

	;# H2-H1 interaction 
	movaps xmm0, [esp + nb204_rinvH2H1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb204_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb204_rsqH2H1] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb204_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + nb204_qqHH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb204_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb204_qqHH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [esp + nb204_fjxH1]
	movaps xmm4, [esp + nb204_fjyH1]
	movaps xmm5, [esp + nb204_fjzH1]
	mulps xmm0, [esp + nb204_dxH2H1]
	mulps xmm1, [esp + nb204_dyH2H1]
	mulps xmm2, [esp + nb204_dzH2H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb204_fixH2]
	addps xmm1, [esp + nb204_fiyH2]
	addps xmm2, [esp + nb204_fizH2]
	movaps [esp + nb204_fjxH1], xmm3
	movaps [esp + nb204_fjyH1], xmm4
	movaps [esp + nb204_fjzH1], xmm5
	movaps [esp + nb204_fixH2], xmm0
	movaps [esp + nb204_fiyH2], xmm1
	movaps [esp + nb204_fizH2], xmm2

	;# H2-H2 interaction 
	movaps xmm0, [esp + nb204_rinvH2H2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb204_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb204_rsqH2H2] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb204_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + nb204_qqHH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb204_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb204_qqHH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [esp + nb204_fjxH2]
	movaps xmm4, [esp + nb204_fjyH2]
	movaps xmm5, [esp + nb204_fjzH2]
	mulps xmm0, [esp + nb204_dxH2H2]
	mulps xmm1, [esp + nb204_dyH2H2]
	mulps xmm2, [esp + nb204_dzH2H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb204_fixH2]
	addps xmm1, [esp + nb204_fiyH2]
	addps xmm2, [esp + nb204_fizH2]
	movaps [esp + nb204_fjxH2], xmm3
	movaps [esp + nb204_fjyH2], xmm4
	movaps [esp + nb204_fjzH2], xmm5
	movaps [esp + nb204_fixH2], xmm0
	movaps [esp + nb204_fiyH2], xmm1
	movaps [esp + nb204_fizH2], xmm2

	;# H2-M interaction 
	movaps xmm0, [esp + nb204_rinvH2M]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb204_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb204_rsqH2M] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb204_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + nb204_qqMH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb204_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb204_qqMH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	movaps xmm3, [esp + nb204_fjxM]
	movaps xmm4, [esp + nb204_fjyM]
	movaps xmm5, [esp + nb204_fjzM]
	mulps xmm0, [esp + nb204_dxH2M]
	mulps xmm1, [esp + nb204_dyH2M]
	mulps xmm2, [esp + nb204_dzH2M]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb204_fixH2]
	addps xmm1, [esp + nb204_fiyH2]
	addps xmm2, [esp + nb204_fizH2]
	movaps [esp + nb204_fjxM], xmm3
	movaps [esp + nb204_fjyM], xmm4
	movaps [esp + nb204_fjzM], xmm5
	movaps [esp + nb204_fixH2], xmm0
	movaps [esp + nb204_fiyH2], xmm1
	movaps [esp + nb204_fizH2], xmm2

	;# M-H1 interaction 
	movaps xmm0, [esp + nb204_rinvMH1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb204_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb204_rsqMH1] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb204_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + nb204_qqMH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb204_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb204_qqMH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [esp + nb204_fjxH1]
	movaps xmm4, [esp + nb204_fjyH1]
	movaps xmm5, [esp + nb204_fjzH1]
	mulps xmm0, [esp + nb204_dxMH1]
	mulps xmm1, [esp + nb204_dyMH1]
	mulps xmm2, [esp + nb204_dzMH1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb204_fixM]
	addps xmm1, [esp + nb204_fiyM]
	addps xmm2, [esp + nb204_fizM]
	movaps [esp + nb204_fjxH1], xmm3
	movaps [esp + nb204_fjyH1], xmm4
	movaps [esp + nb204_fjzH1], xmm5
	movaps [esp + nb204_fixM], xmm0
	movaps [esp + nb204_fiyM], xmm1
	movaps [esp + nb204_fizM], xmm2

	;# M-H2 interaction 
	movaps xmm0, [esp + nb204_rinvMH2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb204_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb204_rsqMH2] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb204_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + nb204_qqMH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb204_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb204_qqMH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [esp + nb204_fjxH2]
	movaps xmm4, [esp + nb204_fjyH2]
	movaps xmm5, [esp + nb204_fjzH2]
	mulps xmm0, [esp + nb204_dxMH2]
	mulps xmm1, [esp + nb204_dyMH2]
	mulps xmm2, [esp + nb204_dzMH2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb204_fixM]
	addps xmm1, [esp + nb204_fiyM]
	addps xmm2, [esp + nb204_fizM]
	movaps [esp + nb204_fjxH2], xmm3
	movaps [esp + nb204_fjyH2], xmm4
	movaps [esp + nb204_fjzH2], xmm5
	movaps [esp + nb204_fixM], xmm0
	movaps [esp + nb204_fiyM], xmm1
	movaps [esp + nb204_fizM], xmm2

	;# M-M interaction 
	movaps xmm0, [esp + nb204_rinvMM]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb204_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb204_rsqMM] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb204_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + nb204_qqMM] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb204_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb204_qqMM] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm1, xmm0
	movaps [esp + nb204_vctot], xmm6
	movaps xmm2, xmm0
	
	movaps xmm3, [esp + nb204_fjxM]
	movaps xmm4, [esp + nb204_fjyM]
	movaps xmm5, [esp + nb204_fjzM]
	mulps xmm0, [esp + nb204_dxMM]
	mulps xmm1, [esp + nb204_dyMM]
	mulps xmm2, [esp + nb204_dzMM]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb204_fixM]
	addps xmm1, [esp + nb204_fiyM]
	addps xmm2, [esp + nb204_fizM]
	movaps [esp + nb204_fjxM], xmm3
	movaps [esp + nb204_fjyM], xmm4
	movaps [esp + nb204_fjzM], xmm5
	movaps [esp + nb204_fixM], xmm0
	movaps [esp + nb204_fiyM], xmm1
	movaps [esp + nb204_fizM], xmm2

	mov edi, [ebp + nb204_faction]
		
	;# Did all interactions - now update j forces 
	;# At this stage forces are still on the stack, in positions:
	;# fjxH1, fjyH1, fjzH1, ... , fjzM.
	;# Each position is a quadruplet of forces for the four 
	;# corresponding j waters, so we need to transpose them before
	;# adding to the memory positions.
	;# 
	;# This _used_ to be a simple transpose, but the resulting high number
	;# of unaligned 128-bit load/stores might trigger a possible hardware 
	;# bug on Athlon and Opteron chips, so I have worked around it
	;# to use 64-bit load/stores instead. The performance hit should be
	;# very modest, since the 128-bit unaligned memory instructions were
	;# slow anyway. 

	;# 4 j waters with three atoms each - first do 1st Hydrogen X & Y forces for 4 j particles 
	movaps xmm0, [esp + nb204_fjxH1] ;# xmm0= fjxH1a  fjxH1b  fjxH1c  fjxH1d 
	movaps xmm2, [esp + nb204_fjyH1] ;# xmm1= fjyH1a  fjyH1b  fjyH1c  fjyH1d
	movlps xmm3, [edi + eax*4 + 12]
	movlps xmm4, [edi + ecx*4 + 12]
	movaps xmm1, xmm0
	unpcklps xmm0, xmm2    	   ;# xmm0= fjxH1a  fjyH1a  fjxH1b  fjyH1b
	unpckhps xmm1, xmm2        ;# xmm1= fjxH1c  fjyH1c  fjxH1d  fjyH1d
	movhps xmm3, [edi + ebx*4 + 12]
	movhps xmm4, [edi + edx*4 + 12]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	movlps [edi + eax*4 + 12], xmm3
	movlps [edi + ecx*4 + 12], xmm4
	movhps [edi + ebx*4 + 12], xmm3
	movhps [edi + edx*4 + 12], xmm4

	;# constant 1st Hydrogen Z & 2nd hydrogen X forces for 4 j particles 
	movaps xmm0, [esp + nb204_fjzH1]  ;# xmm0= fjzH1a   fjzH1b   fjzH1c   fjzH1d 
	movaps xmm2, [esp + nb204_fjxH2] ;# xmm1= fjxH2a  fjxH2b  fjxH2c  fjxH2d
	movlps xmm3, [edi + eax*4 + 20]
	movlps xmm4, [edi + ecx*4 + 20]
	movaps xmm1, xmm0
	unpcklps xmm0, xmm2    	   ;# xmm0= fjzH1a  fjxH2a  fjzH1b  fjxH2b
	unpckhps xmm1, xmm2        ;# xmm1= fjzH1c  fjxH2c  fjzH1d  fjxH2d
	movhps xmm3, [edi + ebx*4 + 20]
	movhps xmm4, [edi + edx*4 + 20]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	movlps [edi + eax*4 + 20], xmm3
	movlps [edi + ecx*4 + 20], xmm4
	movhps [edi + ebx*4 + 20], xmm3
	movhps [edi + edx*4 + 20], xmm4
	
	;# constant 2nd hydrogen Y & Z forces for 4 j particles 
	movaps xmm0, [esp + nb204_fjyH2] ;# xmm0= fjyH2a  fjyH2b  fjyH2c  fjyH2d 
	movaps xmm2, [esp + nb204_fjzH2] ;# xmm1= fjzH2a  fjzH2b  fjzH2c  fjzH2d
	movlps xmm3, [edi + eax*4 + 28]
	movlps xmm4, [edi + ecx*4 + 28]
	movaps xmm1, xmm0
	unpcklps xmm0, xmm2		;# xmm0= fjyH2a  fjzH2a  fjyH2b  fjzH2b
	unpckhps xmm1, xmm2		;# xmm1= fjyH2c  fjzH2c  fjyH2d  fjzH2d
	movhps xmm3, [edi + ebx*4 + 28]
	movhps xmm4, [edi + edx*4 + 28]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	movlps [edi + eax*4 + 28], xmm3
	movlps [edi + ecx*4 + 28], xmm4
	movhps [edi + ebx*4 + 28], xmm3
	movhps [edi + edx*4 + 28], xmm4

	;# Dummy (M) X & Y forces for 4 j particles 
	movaps xmm0, [esp + nb204_fjxM] ;# xmm0= fjxMa  fjxMb  fjxMc  fjxMd 
	movaps xmm2, [esp + nb204_fjyM] ;# xmm1= fjyMa  fjyMb  fjyMc  fjyMd
	movlps xmm3, [edi + eax*4 + 36]
	movlps xmm4, [edi + ecx*4 + 36]
	movaps xmm1, xmm0
	unpcklps xmm0, xmm2		;# xmm0= fjxMa  fjyMa  fjxMb  fjyMb
	unpckhps xmm1, xmm2		;# xmm1= fjxMc  fjyMc  fjxMd  fjyMd
	movhps xmm3, [edi + ebx*4 + 36]
	movhps xmm4, [edi + edx*4 + 36]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	movlps [edi + eax*4 + 36], xmm3
	movlps [edi + ecx*4 + 36], xmm4
	movhps [edi + ebx*4 + 36], xmm3
	movhps [edi + edx*4 + 36], xmm4

	
	;# Dummy (M) Z forces for 4 j particles 
	;# Just load the four Z coords into one reg. each
	movss xmm4, [edi + eax*4 + 44]
	movss xmm5, [edi + ebx*4 + 44]
	movss xmm6, [edi + ecx*4 + 44]
	movss xmm7, [edi + edx*4 + 44]
	;# add what we have on the stack
	addss xmm4, [esp + nb204_fjzM] 
	addss xmm5, [esp + nb204_fjzMb] 
	addss xmm6, [esp + nb204_fjzMc] 
	addss xmm7, [esp + nb204_fjzMd]
	;# store back
	movss [edi + eax*4 + 44], xmm4
	movss [edi + ebx*4 + 44], xmm5
	movss [edi + ecx*4 + 44], xmm6
	movss [edi + edx*4 + 44], xmm7
		
	;# should we do one more iteration? 
	sub dword ptr [esp + nb204_innerk],  4
	jl    .nb204_single_check
	jmp   .nb204_unroll_loop
.nb204_single_check:
	add dword ptr [esp + nb204_innerk],  4
	jnz   .nb204_single_loop
	jmp   .nb204_updateouterdata
.nb204_single_loop:
	mov   edx, [esp + nb204_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb204_innerjjnr],  4	

	mov esi, [ebp + nb204_pos]
	lea   eax, [eax + eax*2]  

	;# fetch j coordinates 
	xorps xmm3, xmm3
	xorps xmm4, xmm4
	xorps xmm5, xmm5
	movss xmm3, [esi + eax*4 + 36]		;# jxM  -  -  -
	movss xmm4, [esi + eax*4 + 40]		;# jyM  -  -  -
	movss xmm5, [esi + eax*4 + 44]		;# jzM  -  -  -  

	movlps xmm6, [esi + eax*4 + 12]		;# xmm6 = jxH1 jyH1   -    -
	movss  xmm7, [esi + eax*4 + 20]		;# xmm7 = jzH1   -    -    - 
	movhps xmm6, [esi + eax*4 + 24]		;# xmm6 = jxH1 jyH1 jxH2 jyH2
	movss  xmm2, [esi + eax*4 + 32]		;# xmm2 = jzH2   -    -    -
	
	;# have all coords, time for some shuffling.

	shufps xmm6, xmm6, 216 ;# constant 11011000	;# xmm6 = jxH1 jxH2 jyH1 jyH2 
	unpcklps xmm7, xmm2			;# xmm7 = jzH1 jzH2   -    -
	movaps  xmm0, [esp + nb204_ixM]     
	movaps  xmm1, [esp + nb204_iyM]
	movaps  xmm2, [esp + nb204_izM]	
	movlhps xmm3, xmm6			;# xmm3 = jxM   0   jxH1 jxH2 
	shufps  xmm4, xmm6, 228 ;# constant 11100100	;# xmm4 = jyM   0   jyH1 jyH2 
	shufps  xmm5, xmm7, 68  ;# constant 01000100	;# xmm5 = jzM   0   jzH1 jzH2
	
	;# store all j coordinates in jM
	movaps [esp + nb204_jxM], xmm3
	movaps [esp + nb204_jyM], xmm4
	movaps [esp + nb204_jzM], xmm5
	subps  xmm0, xmm3
	subps  xmm1, xmm4
	subps  xmm2, xmm5
	movaps [esp + nb204_dxMM], xmm0
	movaps [esp + nb204_dyMM], xmm1
	movaps [esp + nb204_dzMM], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 

	movaps xmm6, xmm0
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	mulps   xmm6, [esp + nb204_krf] ;# xmm6=krsq 
	movaps  xmm2, xmm1
	movaps  xmm7, xmm6 	;# xmm7=krsq 
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + nb204_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [esp + nb204_half] ;# rinv iO - j water 
	
	addps   xmm6, xmm3	;# xmm6=rinv+ krsq 
	mulps   xmm7, [esp + nb204_two]
	subps   xmm6, [esp + nb204_crf] ;# xmm6=rinv+ krsq-crf 
	
	xorps   xmm1, xmm1
	movaps  xmm0, xmm3
	subps   xmm3, xmm7	;# xmm3=rinv-2*krsq 
	xorps   xmm4, xmm4
	mulps   xmm0, xmm0	;# xmm0=rinvsq 
	;# fetch charges to xmm4 (temporary) 
	movss   xmm4, [esp + nb204_qqMM]
	movhps  xmm4, [esp + nb204_qqMH]

	mulps xmm6, xmm4	;# vcoul  
	mulps xmm3, xmm4	;# coul part of fs  


	addps   xmm6, [esp + nb204_vctot]
	mulps   xmm0, xmm3	;# total fscal 
	movaps  [esp + nb204_vctot], xmm6	

	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	mulps   xmm0, [esp + nb204_dxMM]
	mulps   xmm1, [esp + nb204_dyMM]
	mulps   xmm2, [esp + nb204_dzMM]
	
	;# initial update for j forces 
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + nb204_fjxM], xmm3
	movaps  [esp + nb204_fjyM], xmm4
	movaps  [esp + nb204_fjzM], xmm5
	addps   xmm0, [esp + nb204_fixM]
	addps   xmm1, [esp + nb204_fiyM]
	addps   xmm2, [esp + nb204_fizM]
	movaps  [esp + nb204_fixM], xmm0
	movaps  [esp + nb204_fiyM], xmm1
	movaps  [esp + nb204_fizM], xmm2

	
	;# done with i M Now do i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [esp + nb204_ixH1]
	movaps  xmm1, [esp + nb204_iyH1]
	movaps  xmm2, [esp + nb204_izH1]	
	movaps  xmm3, [esp + nb204_ixH2] 
	movaps  xmm4, [esp + nb204_iyH2] 
	movaps  xmm5, [esp + nb204_izH2] 
	subps   xmm0, [esp + nb204_jxM]
	subps   xmm1, [esp + nb204_jyM]
	subps   xmm2, [esp + nb204_jzM]
	subps   xmm3, [esp + nb204_jxM]
	subps   xmm4, [esp + nb204_jyM]
	subps   xmm5, [esp + nb204_jzM]
	movaps [esp + nb204_dxH1M], xmm0
	movaps [esp + nb204_dyH1M], xmm1
	movaps [esp + nb204_dzH1M], xmm2
	movaps [esp + nb204_dxH2M], xmm3
	movaps [esp + nb204_dyH2M], xmm4
	movaps [esp + nb204_dzH2M], xmm5
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
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb204_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb204_half] ;# rinv H1 - j water 
	mulps   xmm7, [esp + nb204_half] ;# rinv H2 - j water  

	mulps xmm0, [esp + nb204_krf] ;# krsq 
	mulps xmm4, [esp + nb204_krf] ;# krsq  

	;# assemble charges in xmm6 
	xorps   xmm6, xmm6
	movss   xmm6, [esp + nb204_qqMH]
	movhps  xmm6, [esp + nb204_qqHH]
	movaps  xmm1, xmm0
	movaps  xmm5, xmm4
	addps   xmm0, xmm3	;# krsq+ rinv 
	addps   xmm4, xmm7	;# krsq+ rinv 
	subps   xmm0, [esp + nb204_crf]
	subps   xmm4, [esp + nb204_crf]
	mulps   xmm1, [esp + nb204_two]
	mulps   xmm5, [esp + nb204_two]
	mulps   xmm0, xmm6	;# vcoul 
	mulps   xmm4, xmm6	;# vcoul 
	addps   xmm4, xmm0		
	addps   xmm4, [esp + nb204_vctot]
	movaps  [esp + nb204_vctot], xmm4
	movaps  xmm0, xmm3
	movaps  xmm4, xmm7
	mulps   xmm3, xmm3
	mulps   xmm7, xmm7
	subps   xmm0, xmm1
	subps   xmm4, xmm5
	mulps   xmm0, xmm6
	mulps   xmm4, xmm6
	mulps   xmm0, xmm3	;# fscal 
	mulps   xmm7, xmm4	;# fscal 
	
	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	mulps   xmm0, [esp + nb204_dxH1M]
	mulps   xmm1, [esp + nb204_dyH1M]
	mulps   xmm2, [esp + nb204_dzH1M]
	;# update forces H1 - j water 
	movaps  xmm3, [esp + nb204_fjxM]
	movaps  xmm4, [esp + nb204_fjyM]
	movaps  xmm5, [esp + nb204_fjzM]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + nb204_fjxM], xmm3
	movaps  [esp + nb204_fjyM], xmm4
	movaps  [esp + nb204_fjzM], xmm5
	addps   xmm0, [esp + nb204_fixH1]
	addps   xmm1, [esp + nb204_fiyH1]
	addps   xmm2, [esp + nb204_fizH1]
	movaps  [esp + nb204_fixH1], xmm0
	movaps  [esp + nb204_fiyH1], xmm1
	movaps  [esp + nb204_fizH1], xmm2
	;# do forces H2 - j water 
	movaps xmm0, xmm7
	movaps xmm1, xmm7
	movaps xmm2, xmm7
	mulps   xmm0, [esp + nb204_dxH2M]
	mulps   xmm1, [esp + nb204_dyH2M]
	mulps   xmm2, [esp + nb204_dzH2M]
	movaps  xmm3, [esp + nb204_fjxM]
	movaps  xmm4, [esp + nb204_fjyM]
	movaps  xmm5, [esp + nb204_fjzM]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	mov     esi, [ebp + nb204_faction]
	movaps  [esp + nb204_fjxM], xmm3
	movaps  [esp + nb204_fjyM], xmm4
	movaps  [esp + nb204_fjzM], xmm5
	addps   xmm0, [esp + nb204_fixH2]
	addps   xmm1, [esp + nb204_fiyH2]
	addps   xmm2, [esp + nb204_fizH2]
	movaps  [esp + nb204_fixH2], xmm0
	movaps  [esp + nb204_fiyH2], xmm1
	movaps  [esp + nb204_fizH2], xmm2

	;# update j water forces from local variables 
	movlps  xmm0, [esi + eax*4 + 36]
	movlps  xmm1, [esi + eax*4 + 12]
	movhps  xmm1, [esi + eax*4 + 24]
	movaps  xmm3, [esp + nb204_fjxM]
	movaps  xmm4, [esp + nb204_fjyM]
	movaps  xmm5, [esp + nb204_fjzM]
	movaps  xmm6, xmm5
	movaps  xmm7, xmm5
	shufps  xmm6, xmm6, 2 ;# constant 00000010
	shufps  xmm7, xmm7, 3 ;# constant 00000011
	addss   xmm5, [esi + eax*4 + 44]
	addss   xmm6, [esi + eax*4 + 20]
	addss   xmm7, [esi + eax*4 + 32]
	movss   [esi + eax*4 + 44], xmm5
	movss   [esi + eax*4 + 20], xmm6
	movss   [esi + eax*4 + 32], xmm7
	movaps   xmm5, xmm3
	unpcklps xmm3, xmm4
	unpckhps xmm5, xmm4
	addps    xmm0, xmm3
	addps    xmm1, xmm5
	movlps  [esi + eax*4 + 36], xmm0 
	movlps  [esi + eax*4 + 12], xmm1 
	movhps  [esi + eax*4 + 24], xmm1 
	
	dec dword ptr [esp + nb204_innerk]
	jz    .nb204_updateouterdata
	jmp   .nb204_single_loop
.nb204_updateouterdata:
	mov   ecx, [esp + nb204_ii3]
	mov   edi, [ebp + nb204_faction]
	mov   esi, [ebp + nb204_fshift]
	mov   edx, [esp + nb204_is3]

	;# accumulate  H1 i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb204_fixH1]
	movaps xmm1, [esp + nb204_fiyH1] 
	movaps xmm2, [esp + nb204_fizH1]

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
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# constant 00001000	
 
	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb204_fixH2]
	movaps xmm1, [esp + nb204_fiyH2]
	movaps xmm2, [esp + nb204_fizH2]

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

	;# accumulate M i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb204_fixM]
	movaps xmm1, [esp + nb204_fiyM]
	movaps xmm2, [esp + nb204_fizM]

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
	mov esi, [esp + nb204_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb204_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb204_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb204_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 	
	
        ;# finish if last 
        mov ecx, [esp + nb204_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb204_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb204_n], esi
        jmp .nb204_outer
.nb204_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb204_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb204_end
        ;# non-zero, do one more workunit
        jmp   .nb204_threadloop
.nb204_end:
	emms

	mov eax, [esp + nb204_nouter]
	mov ebx, [esp + nb204_ninner]
	mov ecx, [ebp + nb204_outeriter]
	mov edx, [ebp + nb204_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb204_salign]
	add esp, eax
	add esp, 1480
	pop edi
	pop esi
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret



	
.globl nb_kernel204nf_ia32_sse
.globl _nb_kernel204nf_ia32_sse
nb_kernel204nf_ia32_sse:	
_nb_kernel204nf_ia32_sse:	
.equiv          nb204nf_p_nri,          8
.equiv          nb204nf_iinr,           12
.equiv          nb204nf_jindex,         16
.equiv          nb204nf_jjnr,           20
.equiv          nb204nf_shift,          24
.equiv          nb204nf_shiftvec,       28
.equiv          nb204nf_fshift,         32
.equiv          nb204nf_gid,            36
.equiv          nb204nf_pos,            40
.equiv          nb204nf_faction,        44
.equiv          nb204nf_charge,         48
.equiv          nb204nf_p_facel,        52
.equiv          nb204nf_argkrf,         56
.equiv          nb204nf_argcrf,         60
.equiv          nb204nf_Vc,             64
.equiv          nb204nf_type,           68
.equiv          nb204nf_p_ntype,        72
.equiv          nb204nf_vdwparam,       76
.equiv          nb204nf_Vvdw,           80
.equiv          nb204nf_p_tabscale,     84
.equiv          nb204nf_VFtab,          88
.equiv          nb204nf_invsqrta,       92
.equiv          nb204nf_dvda,           96
.equiv          nb204nf_p_gbtabscale,   100
.equiv          nb204nf_GBtab,          104
.equiv          nb204nf_p_nthreads,     108
.equiv          nb204nf_count,          112
.equiv          nb204nf_mtx,            116
.equiv          nb204nf_outeriter,      120
.equiv          nb204nf_inneriter,      124
.equiv          nb204nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb204nf_ixH1,           0
.equiv          nb204nf_iyH1,           16
.equiv          nb204nf_izH1,           32
.equiv          nb204nf_ixH2,           48
.equiv          nb204nf_iyH2,           64
.equiv          nb204nf_izH2,           80
.equiv          nb204nf_ixM,            96
.equiv          nb204nf_iyM,            112
.equiv          nb204nf_izM,            128
.equiv          nb204nf_jxH1,           144
.equiv          nb204nf_jyH1,           160
.equiv          nb204nf_jzH1,           176
.equiv          nb204nf_jxH2,           192
.equiv          nb204nf_jyH2,           208
.equiv          nb204nf_jzH2,           224
.equiv          nb204nf_jxM,            240
.equiv          nb204nf_jyM,            256
.equiv          nb204nf_jzM,            272
.equiv          nb204nf_qqHH,           288
.equiv          nb204nf_qqMH,           304
.equiv          nb204nf_qqMM,           320
.equiv          nb204nf_vctot,          336
.equiv          nb204nf_half,           352
.equiv          nb204nf_three,          368
.equiv          nb204nf_rsqH1H1,        384
.equiv          nb204nf_rsqH1H2,        400
.equiv          nb204nf_rsqH1M,         416
.equiv          nb204nf_rsqH2H1,        432
.equiv          nb204nf_rsqH2H2,        448
.equiv          nb204nf_rsqH2M,         464
.equiv          nb204nf_rsqMH1,         480
.equiv          nb204nf_rsqMH2,         496
.equiv          nb204nf_rsqMM,          512
.equiv          nb204nf_rinvH1H1,       528
.equiv          nb204nf_rinvH1H2,       544
.equiv          nb204nf_rinvH1M,        560
.equiv          nb204nf_rinvH2H1,       576
.equiv          nb204nf_rinvH2H2,       592
.equiv          nb204nf_rinvH2M,        608
.equiv          nb204nf_rinvMH1,        624
.equiv          nb204nf_rinvMH2,        640
.equiv          nb204nf_rinvMM,         656
.equiv          nb204nf_krf,            672
.equiv          nb204nf_crf,            688
.equiv          nb204nf_is3,            704
.equiv          nb204nf_ii3,            708
.equiv          nb204nf_innerjjnr,      712
.equiv          nb204nf_innerk,         716
.equiv          nb204nf_n,              720
.equiv          nb204nf_nn1,            724
.equiv          nb204nf_nri,            728
.equiv          nb204nf_nouter,         732
.equiv          nb204nf_ninner,         736
.equiv          nb204nf_salign,         740
	push ebp
	mov ebp,esp
	push eax
	push ebx
	push ecx
	push edx
	push esi
	push edi
	sub esp, 744		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb204nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb204nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb204nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb204nf_nouter], eax
	mov [esp + nb204nf_ninner], eax


	mov esi, [ebp + nb204nf_argkrf]
	mov edi, [ebp + nb204nf_argcrf]
	movss xmm5, [esi]
	movss xmm6, [edi]
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	movaps [esp + nb204nf_krf], xmm5
	movaps [esp + nb204nf_crf], xmm6
	
	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb204nf_half], eax
	movss xmm1, [esp + nb204nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb204nf_half],  xmm1
	movaps [esp + nb204nf_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb204nf_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]		;# ebx =ii 

	mov   edx, [ebp + nb204nf_charge]
	movss xmm3, [edx + ebx*4 + 4]	
	movss xmm4, xmm3	
	movss xmm5, [edx + ebx*4 + 12]	
	mov esi, [ebp + nb204nf_p_facel]
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
	movaps [esp + nb204nf_qqHH], xmm3
	movaps [esp + nb204nf_qqMH], xmm4
	movaps [esp + nb204nf_qqMM], xmm5
	
.nb204nf_threadloop:
        mov   esi, [ebp + nb204nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb204nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb204nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb204nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb204nf_n], eax
        mov [esp + nb204nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb204nf_outerstart
        jmp .nb204nf_end
	
.nb204nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb204nf_nouter]
	mov [esp + nb204nf_nouter], ebx

.nb204nf_outer:
	mov   eax, [ebp + nb204nf_shift]  	;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]	;# ebx=3*is 
	mov   [esp + nb204nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb204nf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb204nf_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]		;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb204nf_pos]	;# eax = base of pos[]  
	mov   [esp + nb204nf_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [eax + ebx*4 + 12]
	addss xmm4, [eax + ebx*4 + 16]
	addss xmm5, [eax + ebx*4 + 20]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb204nf_ixH1], xmm3
	movaps [esp + nb204nf_iyH1], xmm4
	movaps [esp + nb204nf_izH1], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 24]
	addss xmm1, [eax + ebx*4 + 28]
	addss xmm2, [eax + ebx*4 + 32]		
	addss xmm3, [eax + ebx*4 + 36]
	addss xmm4, [eax + ebx*4 + 40]
	addss xmm5, [eax + ebx*4 + 44]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb204nf_ixH2], xmm0
	movaps [esp + nb204nf_iyH2], xmm1
	movaps [esp + nb204nf_izH2], xmm2
	movaps [esp + nb204nf_ixM], xmm3
	movaps [esp + nb204nf_iyM], xmm4
	movaps [esp + nb204nf_izM], xmm5

	;# clear vctot 
	xorps xmm4, xmm4
	movaps [esp + nb204nf_vctot], xmm4
	
	mov   eax, [ebp + nb204nf_jindex]
	mov   ecx, [eax + esi*4]	 	;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   esi, [ebp + nb204nf_pos]
	mov   edi, [ebp + nb204nf_faction]	
	mov   eax, [ebp + nb204nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb204nf_innerjjnr], eax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb204nf_ninner]
	mov   [esp + nb204nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb204nf_innerk], edx	;# number of innerloop atoms 
	jge   .nb204nf_unroll_loop
	jmp   .nb204nf_single_check
.nb204nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb204nf_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	mov   ecx, [edx + 8]
	mov   edx, [edx + 12]     	;# eax-edx=jnr1-4 
	
	add dword ptr [esp + nb204nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb204nf_pos]   	;# base of pos[] 

	lea   eax, [eax + eax*2] 	;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2] 	;# replace jnr with j3 
	lea   edx, [edx + edx*2]	
	
	;# move j coordinates to local temp variables 
	movlps xmm2, [esi + eax*4 + 12]
	movlps xmm3, [esi + eax*4 + 24]
	movlps xmm4, [esi + eax*4 + 36]

	movlps xmm5, [esi + ebx*4 + 12]
	movlps xmm6, [esi + ebx*4 + 24]
	movlps xmm7, [esi + ebx*4 + 36]

	movhps xmm2, [esi + ecx*4 + 12]
	movhps xmm3, [esi + ecx*4 + 24]
	movhps xmm4, [esi + ecx*4 + 36]

	movhps xmm5, [esi + edx*4 + 12]
	movhps xmm6, [esi + edx*4 + 24]
	movhps xmm7, [esi + edx*4 + 36]
	
	movaps xmm0, xmm2
	movaps xmm1, xmm3
	unpcklps xmm0, xmm5
	unpcklps xmm1, xmm6
	unpckhps xmm2, xmm5
	unpckhps xmm3, xmm6
	movaps xmm5, xmm4
	movaps   xmm6, xmm0
	unpcklps xmm4, xmm7
	unpckhps xmm5, xmm7
	movaps   xmm7, xmm1
	movlhps  xmm0, xmm2
	movaps [esp + nb204nf_jxH1], xmm0
	movhlps  xmm2, xmm6
	movaps [esp + nb204nf_jyH1], xmm2
	movlhps  xmm1, xmm3
	movaps [esp + nb204nf_jxH2], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [esp + nb204nf_jyH2], xmm3
	movlhps  xmm4, xmm5
	movaps [esp + nb204nf_jxM], xmm4
	movhlps  xmm5, xmm6
	movaps [esp + nb204nf_jyM], xmm5

	movss  xmm0, [esi + eax*4 + 20]
	movss  xmm1, [esi + eax*4 + 32]
	movss  xmm2, [esi + eax*4 + 44]

	movss  xmm3, [esi + ecx*4 + 20]
	movss  xmm4, [esi + ecx*4 + 32]
	movss  xmm5, [esi + ecx*4 + 44]

	movhps xmm0, [esi + ebx*4 + 16]
	movhps xmm1, [esi + ebx*4 + 28]
	movhps xmm2, [esi + ebx*4 + 40]
	
	movhps xmm3, [esi + edx*4 + 16]
	movhps xmm4, [esi + edx*4 + 28]
	movhps xmm5, [esi + edx*4 + 40]
	
	shufps xmm0, xmm3, 204  ;# constant 11001100
	shufps xmm1, xmm4, 204  ;# constant 11001100
	shufps xmm2, xmm5, 204  ;# constant 11001100
	movaps [esp + nb204nf_jzH1],  xmm0
	movaps [esp + nb204nf_jzH2],  xmm1
	movaps [esp + nb204nf_jzM],  xmm2

	movaps xmm0, [esp + nb204nf_ixH1]
	movaps xmm1, [esp + nb204nf_iyH1]
	movaps xmm2, [esp + nb204nf_izH1]
	movaps xmm3, [esp + nb204nf_ixH1]
	movaps xmm4, [esp + nb204nf_iyH1]
	movaps xmm5, [esp + nb204nf_izH1]
	subps  xmm0, [esp + nb204nf_jxH1]
	subps  xmm1, [esp + nb204nf_jyH1]
	subps  xmm2, [esp + nb204nf_jzH1]
	subps  xmm3, [esp + nb204nf_jxH2]
	subps  xmm4, [esp + nb204nf_jyH2]
	subps  xmm5, [esp + nb204nf_jzH2]
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
	movaps [esp + nb204nf_rsqH1H1], xmm0
	movaps [esp + nb204nf_rsqH1H2], xmm3

	movaps xmm0, [esp + nb204nf_ixH1]
	movaps xmm1, [esp + nb204nf_iyH1]
	movaps xmm2, [esp + nb204nf_izH1]
	movaps xmm3, [esp + nb204nf_ixH2]
	movaps xmm4, [esp + nb204nf_iyH2]
	movaps xmm5, [esp + nb204nf_izH2]
	subps  xmm0, [esp + nb204nf_jxM]
	subps  xmm1, [esp + nb204nf_jyM]
	subps  xmm2, [esp + nb204nf_jzM]
	subps  xmm3, [esp + nb204nf_jxH1]
	subps  xmm4, [esp + nb204nf_jyH1]
	subps  xmm5, [esp + nb204nf_jzH1]
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
	movaps [esp + nb204nf_rsqH1M], xmm0
	movaps [esp + nb204nf_rsqH2H1], xmm3

	movaps xmm0, [esp + nb204nf_ixH2]
	movaps xmm1, [esp + nb204nf_iyH2]
	movaps xmm2, [esp + nb204nf_izH2]
	movaps xmm3, [esp + nb204nf_ixH2]
	movaps xmm4, [esp + nb204nf_iyH2]
	movaps xmm5, [esp + nb204nf_izH2]
	subps  xmm0, [esp + nb204nf_jxH2]
	subps  xmm1, [esp + nb204nf_jyH2]
	subps  xmm2, [esp + nb204nf_jzH2]
	subps  xmm3, [esp + nb204nf_jxM]
	subps  xmm4, [esp + nb204nf_jyM]
	subps  xmm5, [esp + nb204nf_jzM]
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
	movaps [esp + nb204nf_rsqH2H2], xmm0
	movaps [esp + nb204nf_rsqH2M], xmm3

	movaps xmm0, [esp + nb204nf_ixM]
	movaps xmm1, [esp + nb204nf_iyM]
	movaps xmm2, [esp + nb204nf_izM]
	movaps xmm3, [esp + nb204nf_ixM]
	movaps xmm4, [esp + nb204nf_iyM]
	movaps xmm5, [esp + nb204nf_izM]
	subps  xmm0, [esp + nb204nf_jxH1]
	subps  xmm1, [esp + nb204nf_jyH1]
	subps  xmm2, [esp + nb204nf_jzH1]
	subps  xmm3, [esp + nb204nf_jxH2]
	subps  xmm4, [esp + nb204nf_jyH2]
	subps  xmm5, [esp + nb204nf_jzH2]
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
	movaps [esp + nb204nf_rsqMH1], xmm0
	movaps [esp + nb204nf_rsqMH2], xmm4

	movaps xmm0, [esp + nb204nf_ixM]
	movaps xmm1, [esp + nb204nf_iyM]
	movaps xmm2, [esp + nb204nf_izM]
	subps  xmm0, [esp + nb204nf_jxM]
	subps  xmm1, [esp + nb204nf_jyM]
	subps  xmm2, [esp + nb204nf_jzM]
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [esp + nb204nf_rsqMM], xmm0
		
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb204nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb204nf_half] ;# rinvH2H2 
	mulps   xmm7, [esp + nb204nf_half] ;# rinvH2H1 
	movaps  [esp + nb204nf_rinvMM], xmm3
	movaps  [esp + nb204nf_rinvMH2], xmm7
	
	rsqrtps xmm1, [esp + nb204nf_rsqH1H1]
	rsqrtps xmm5, [esp + nb204nf_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb204nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb204nf_rsqH1H1]
	mulps   xmm5, [esp + nb204nf_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb204nf_half] 
	mulps   xmm7, [esp + nb204nf_half]
	movaps  [esp + nb204nf_rinvH1H1], xmm3
	movaps  [esp + nb204nf_rinvH1H2], xmm7
	
	rsqrtps xmm1, [esp + nb204nf_rsqH1M]
	rsqrtps xmm5, [esp + nb204nf_rsqH2H1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb204nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb204nf_rsqH1M]
	mulps   xmm5, [esp + nb204nf_rsqH2H1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb204nf_half] 
	mulps   xmm7, [esp + nb204nf_half]
	movaps  [esp + nb204nf_rinvH1M], xmm3
	movaps  [esp + nb204nf_rinvH2H1], xmm7
	
	rsqrtps xmm1, [esp + nb204nf_rsqH2H2]
	rsqrtps xmm5, [esp + nb204nf_rsqH2M]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb204nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb204nf_rsqH2H2]
	mulps   xmm5, [esp + nb204nf_rsqH2M]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb204nf_half] 
	mulps   xmm7, [esp + nb204nf_half]
	movaps  [esp + nb204nf_rinvH2H2], xmm3
	movaps  [esp + nb204nf_rinvH2M], xmm7
	
	rsqrtps xmm1, [esp + nb204nf_rsqMH1]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + nb204nf_three]
	mulps   xmm1, [esp + nb204nf_rsqMH1]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [esp + nb204nf_half] 
	movaps  [esp + nb204nf_rinvMH1], xmm3

	;# Coulomb interactions 
	;# add all H-H rsq in xmm2, H-M rsq xmm4
	;# H-H rinv in xmm0, H-M in xmm1
	movaps xmm0, [esp + nb204nf_rinvH1H1]
	movaps xmm1, [esp + nb204nf_rinvH1M]
	movaps xmm2, [esp + nb204nf_rsqH1H1]
	movaps xmm4, [esp + nb204nf_rsqH1M]
	addps  xmm0, [esp + nb204nf_rinvH1H2]
	addps  xmm1, [esp + nb204nf_rinvH2M]
	addps  xmm2, [esp + nb204nf_rsqH1H2]
	addps  xmm4, [esp + nb204nf_rsqH2M]
	addps  xmm0, [esp + nb204nf_rinvH2H1]
	addps  xmm1, [esp + nb204nf_rinvMH1]
	addps  xmm2, [esp + nb204nf_rsqH2H1]
	addps  xmm4, [esp + nb204nf_rsqMH1]
	addps  xmm0, [esp + nb204nf_rinvH2H2]
	addps  xmm1, [esp + nb204nf_rinvMH2]
	addps  xmm2, [esp + nb204nf_rsqH2H2]
	addps  xmm4, [esp + nb204nf_rsqMH2]
	movaps xmm5, [esp + nb204nf_krf]
	movaps xmm6, [esp + nb204nf_crf]

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
	mulps  xmm0, [esp + nb204nf_qqHH]
	mulps  xmm1, [esp + nb204nf_qqMH]
	addps  xmm0, xmm1 
	addps  xmm0, [esp + nb204nf_vctot]
	;# M-M interaction
	movaps xmm4, [esp + nb204nf_rinvMM]
	mulps  xmm5, [esp + nb204nf_rsqMM] ;# krsq
	addps  xmm5, xmm4                  ;# rinv+krsq
	subps xmm5, [esp + nb204nf_crf] ;# xmm5=rinv+ krsq-crf 
	mulps xmm5, [esp + nb204nf_qqMM]
	addps xmm5, xmm0
	movaps [esp + nb204nf_vctot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [esp + nb204nf_innerk],  4
	jl    .nb204nf_single_check
	jmp   .nb204nf_unroll_loop
.nb204nf_single_check:
	add dword ptr [esp + nb204nf_innerk],  4
	jnz   .nb204nf_single_loop
	jmp   .nb204nf_updateouterdata
.nb204nf_single_loop:
	mov   edx, [esp + nb204nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb204nf_innerjjnr],  4	

	mov esi, [ebp + nb204nf_pos]
	lea   eax, [eax + eax*2]  

	;# fetch j coordinates 
	xorps xmm3, xmm3
	xorps xmm4, xmm4
	xorps xmm5, xmm5
	movss xmm3, [esi + eax*4 + 36]		;# jxM  -  -  -
	movss xmm4, [esi + eax*4 + 40]		;# jyM  -  -  -
	movss xmm5, [esi + eax*4 + 44]		;# jzM  -  -  -  

	movlps xmm6, [esi + eax*4 + 12]		;# xmm6 = jxH1 jyH1   -    -
	movss  xmm7, [esi + eax*4 + 20]		;# xmm7 = jzH1   -    -    - 
	movhps xmm6, [esi + eax*4 + 24]		;# xmm6 = jxH1 jyH1 jxH2 jyH2
	movss  xmm2, [esi + eax*4 + 32]		;# xmm2 = jzH2   -    -    -
	
	;# have all coords, time for some shuffling.

	shufps xmm6, xmm6, 216 ;# constant 11011000	;# xmm6 = jxH1 jxH2 jyH1 jyH2 
	unpcklps xmm7, xmm2			;# xmm7 = jzH1 jzH2   -    -
	movaps  xmm0, [esp + nb204nf_ixM]     
	movaps  xmm1, [esp + nb204nf_iyM]
	movaps  xmm2, [esp + nb204nf_izM]	
	movlhps xmm3, xmm6			;# xmm3 = jxM   0   jxH1 jxH2 
	shufps  xmm4, xmm6, 228 ;# constant 11100100	;# xmm4 = jyM   0   jyH1 jyH2 
	shufps  xmm5, xmm7, 68  ;# constant 01000100	;# xmm5 = jzM   0   jzH1 jzH2
	
	;# store all j coordinates in jM
	movaps [esp + nb204nf_jxM], xmm3
	movaps [esp + nb204nf_jyM], xmm4
	movaps [esp + nb204nf_jzM], xmm5
	subps  xmm0, xmm3
	subps  xmm1, xmm4
	subps  xmm2, xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 

	movaps xmm6, xmm0
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	mulps   xmm6, [esp + nb204nf_krf] ;# xmm6=krsq 
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + nb204nf_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [esp + nb204nf_half] ;# rinv iO - j water 
	
	addps   xmm6, xmm3	;# xmm6=rinv+ krsq 
	subps   xmm6, [esp + nb204nf_crf] ;# xmm6=rinv+ krsq-crf 
	
	xorps   xmm4, xmm4
	;# fetch charges to xmm4 (temporary) 
	movss   xmm4, [esp + nb204nf_qqMM]
	movhps  xmm4, [esp + nb204nf_qqMH]
	mulps xmm6, xmm4	;# vcoul  
	addps   xmm6, [esp + nb204nf_vctot]
	movaps  [esp + nb204nf_vctot], xmm6	
	
	;# done with i M Now do i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [esp + nb204nf_ixH1]
	movaps  xmm1, [esp + nb204nf_iyH1]
	movaps  xmm2, [esp + nb204nf_izH1]	
	movaps  xmm3, [esp + nb204nf_ixH2] 
	movaps  xmm4, [esp + nb204nf_iyH2] 
	movaps  xmm5, [esp + nb204nf_izH2] 
	subps   xmm0, [esp + nb204nf_jxM]
	subps   xmm1, [esp + nb204nf_jyM]
	subps   xmm2, [esp + nb204nf_jzM]
	subps   xmm3, [esp + nb204nf_jxM]
	subps   xmm4, [esp + nb204nf_jyM]
	subps   xmm5, [esp + nb204nf_jzM]
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
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb204nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb204nf_half] ;# rinv H1 - j water 
	mulps   xmm7, [esp + nb204nf_half] ;# rinv H2 - j water  

	mulps xmm0, [esp + nb204nf_krf] ;# krsq 
	mulps xmm4, [esp + nb204nf_krf] ;# krsq  

	;# assemble charges in xmm6 
	xorps   xmm6, xmm6
	movss   xmm6, [esp + nb204nf_qqMH]
	movhps  xmm6, [esp + nb204nf_qqHH]
	addps   xmm0, xmm3	;# krsq+ rinv 
	addps   xmm4, xmm7	;# krsq+ rinv 
	subps   xmm0, [esp + nb204nf_crf]
	subps   xmm4, [esp + nb204nf_crf]
	mulps   xmm0, xmm6	;# vcoul 
	mulps   xmm4, xmm6	;# vcoul 
	addps   xmm4, xmm0		
	addps   xmm4, [esp + nb204nf_vctot]
	movaps  [esp + nb204nf_vctot], xmm4
	dec dword ptr [esp + nb204nf_innerk]
	jz    .nb204nf_updateouterdata
	jmp   .nb204nf_single_loop
.nb204nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb204nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb204nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb204nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb204nf_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 	
	
        ;# finish if last 
        mov ecx, [esp + nb204nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb204nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb204nf_n], esi
        jmp .nb204nf_outer
.nb204nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb204nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb204nf_end
        ;# non-zero, do one more workunit
        jmp   .nb204nf_threadloop
.nb204nf_end:
	emms

	mov eax, [esp + nb204nf_nouter]
	mov ebx, [esp + nb204nf_ninner]
	mov ecx, [ebp + nb204nf_outeriter]
	mov edx, [ebp + nb204nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb204nf_salign]
	add esp, eax
	add esp, 744
	pop edi
	pop esi
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret
