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


.globl nb_kernel132_ia32_sse
.globl _nb_kernel132_ia32_sse
nb_kernel132_ia32_sse:	
_nb_kernel132_ia32_sse:	
.equiv          nb132_p_nri,            8
.equiv          nb132_iinr,             12
.equiv          nb132_jindex,           16
.equiv          nb132_jjnr,             20
.equiv          nb132_shift,            24
.equiv          nb132_shiftvec,         28
.equiv          nb132_fshift,           32
.equiv          nb132_gid,              36
.equiv          nb132_pos,              40
.equiv          nb132_faction,          44
.equiv          nb132_charge,           48
.equiv          nb132_p_facel,          52
.equiv          nb132_argkrf,           56
.equiv          nb132_argcrf,           60
.equiv          nb132_Vc,               64
.equiv          nb132_type,             68
.equiv          nb132_p_ntype,          72
.equiv          nb132_vdwparam,         76
.equiv          nb132_Vvdw,             80
.equiv          nb132_p_tabscale,       84
.equiv          nb132_VFtab,            88
.equiv          nb132_invsqrta,         92
.equiv          nb132_dvda,             96
.equiv          nb132_p_gbtabscale,     100
.equiv          nb132_GBtab,            104
.equiv          nb132_p_nthreads,       108
.equiv          nb132_count,            112
.equiv          nb132_mtx,              116
.equiv          nb132_outeriter,        120
.equiv          nb132_inneriter,        124
.equiv          nb132_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb132_ixO,              0
.equiv          nb132_iyO,              16
.equiv          nb132_izO,              32
.equiv          nb132_ixH1,             48
.equiv          nb132_iyH1,             64
.equiv          nb132_izH1,             80
.equiv          nb132_ixH2,             96
.equiv          nb132_iyH2,             112
.equiv          nb132_izH2,             128
.equiv          nb132_jxO,              144
.equiv          nb132_jyO,              160
.equiv          nb132_jzO,              176
.equiv          nb132_jxH1,             192
.equiv          nb132_jyH1,             208
.equiv          nb132_jzH1,             224
.equiv          nb132_jxH2,             240
.equiv          nb132_jyH2,             256
.equiv          nb132_jzH2,             272
.equiv          nb132_dxOO,             288
.equiv          nb132_dyOO,             304
.equiv          nb132_dzOO,             320
.equiv          nb132_dxOH1,            336
.equiv          nb132_dyOH1,            352
.equiv          nb132_dzOH1,            368
.equiv          nb132_dxOH2,            384
.equiv          nb132_dyOH2,            400
.equiv          nb132_dzOH2,            416
.equiv          nb132_dxH1O,            432
.equiv          nb132_dyH1O,            448
.equiv          nb132_dzH1O,            464
.equiv          nb132_dxH1H1,           480
.equiv          nb132_dyH1H1,           496
.equiv          nb132_dzH1H1,           512
.equiv          nb132_dxH1H2,           528
.equiv          nb132_dyH1H2,           544
.equiv          nb132_dzH1H2,           560
.equiv          nb132_dxH2O,            576
.equiv          nb132_dyH2O,            592
.equiv          nb132_dzH2O,            608
.equiv          nb132_dxH2H1,           624
.equiv          nb132_dyH2H1,           640
.equiv          nb132_dzH2H1,           656
.equiv          nb132_dxH2H2,           672
.equiv          nb132_dyH2H2,           688
.equiv          nb132_dzH2H2,           704
.equiv          nb132_qqOO,             720
.equiv          nb132_qqOH,             736
.equiv          nb132_qqHH,             752
.equiv          nb132_c6,               768
.equiv          nb132_c12,              784
.equiv          nb132_tsc,              800
.equiv          nb132_fstmp,            816
.equiv          nb132_vctot,            832
.equiv          nb132_Vvdwtot,          848
.equiv          nb132_fixO,             864
.equiv          nb132_fiyO,             880
.equiv          nb132_fizO,             896
.equiv          nb132_fixH1,            912
.equiv          nb132_fiyH1,            928
.equiv          nb132_fizH1,            944
.equiv          nb132_fixH2,            960
.equiv          nb132_fiyH2,            976
.equiv          nb132_fizH2,            992
.equiv          nb132_fjxO,             1008
.equiv          nb132_fjyO,             1024
.equiv          nb132_fjzO,             1040
.equiv          nb132_fjxH1,            1056
.equiv          nb132_fjyH1,            1072
.equiv          nb132_fjzH1,            1088
.equiv          nb132_fjxH2,            1104
.equiv          nb132_fjyH2,            1120
.equiv          nb132_fjzH2,            1136
.equiv          nb132_fjzH2b,           1140
.equiv          nb132_fjzH2c,           1144
.equiv          nb132_fjzH2d,           1148
.equiv          nb132_half,             1152
.equiv          nb132_three,            1168
.equiv          nb132_rsqOO,            1184
.equiv          nb132_rsqOH1,           1200
.equiv          nb132_rsqOH2,           1216
.equiv          nb132_rsqH1O,           1232
.equiv          nb132_rsqH1H1,          1248
.equiv          nb132_rsqH1H2,          1264
.equiv          nb132_rsqH2O,           1280
.equiv          nb132_rsqH2H1,          1296
.equiv          nb132_rsqH2H2,          1312
.equiv          nb132_rinvOO,           1328
.equiv          nb132_rinvOH1,          1344
.equiv          nb132_rinvOH2,          1360
.equiv          nb132_rinvH1O,          1376
.equiv          nb132_rinvH1H1,         1392
.equiv          nb132_rinvH1H2,         1408
.equiv          nb132_rinvH2O,          1424
.equiv          nb132_rinvH2H1,         1440
.equiv          nb132_rinvH2H2,         1456
.equiv          nb132_two,              1472
.equiv          nb132_is3,              1520
.equiv          nb132_ii3,              1524
.equiv          nb132_innerjjnr,        1528
.equiv          nb132_innerk,           1532
.equiv          nb132_n,                1536
.equiv          nb132_nn1,              1540
.equiv          nb132_nri,              1544
.equiv          nb132_nouter,           1548
.equiv          nb132_ninner,           1552
.equiv          nb132_salign,           1556
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 1560		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb132_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb132_p_nri]
	mov ecx, [ecx]
	mov [esp + nb132_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb132_nouter], eax
	mov [esp + nb132_ninner], eax

	mov eax, [ebp + nb132_p_tabscale]
	movss xmm3, [eax]
	shufps xmm3, xmm3, 0
	movaps [esp + nb132_tsc], xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb132_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb132_charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, xmm3	
	movss xmm5, [edx + ebx*4 + 4]	
	mov esi, [ebp + nb132_p_facel]
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
	movaps [esp + nb132_qqOO], xmm3
	movaps [esp + nb132_qqOH], xmm4
	movaps [esp + nb132_qqHH], xmm5
		
	xorps xmm0, xmm0
	mov   edx, [ebp + nb132_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov edi, [ebp + nb132_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   eax, [ebp + nb132_vdwparam]
	movlps xmm0, [eax + edx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 85  ;# constant 01010101
	movaps [esp + nb132_c6], xmm0
	movaps [esp + nb132_c12], xmm1

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb132_half], eax
	movss xmm1, [esp + nb132_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb132_half],  xmm1
	movaps [esp + nb132_two],  xmm2
	movaps [esp + nb132_three],  xmm3

.nb132_threadloop:
        mov   esi, [ebp + nb132_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb132_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb132_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb132_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb132_n], eax
        mov [esp + nb132_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb132_outerstart
        jmp .nb132_end

.nb132_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb132_nouter]
	mov [esp + nb132_nouter], ebx

.nb132_outer:
	mov   eax, [ebp + nb132_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb132_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb132_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb132_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb132_pos]    ;# eax = base of pos[]  
	mov   [esp + nb132_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb132_ixO], xmm3
	movaps [esp + nb132_iyO], xmm4
	movaps [esp + nb132_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb132_ixH1], xmm0
	movaps [esp + nb132_iyH1], xmm1
	movaps [esp + nb132_izH1], xmm2
	movaps [esp + nb132_ixH2], xmm3
	movaps [esp + nb132_iyH2], xmm4
	movaps [esp + nb132_izH2], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb132_vctot], xmm4
	movaps [esp + nb132_Vvdwtot], xmm4
	movaps [esp + nb132_fixO], xmm4
	movaps [esp + nb132_fiyO], xmm4
	movaps [esp + nb132_fizO], xmm4
	movaps [esp + nb132_fixH1], xmm4
	movaps [esp + nb132_fiyH1], xmm4
	movaps [esp + nb132_fizH1], xmm4
	movaps [esp + nb132_fixH2], xmm4
	movaps [esp + nb132_fiyH2], xmm4
	movaps [esp + nb132_fizH2], xmm4
	
	mov   eax, [ebp + nb132_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb132_pos]
	mov   edi, [ebp + nb132_faction]	
	mov   eax, [ebp + nb132_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb132_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb132_ninner]
	mov   [esp + nb132_ninner], ecx
	add   edx, 0
	mov   [esp + nb132_innerk], edx    ;# number of innerloop atoms 
	jge   .nb132_unroll_loop
	jmp   .nb132_single_check
.nb132_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb132_innerjjnr]     ;# pointer to jjnr[k] 

	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	mov   ecx, [edx + 8]
	mov   edx, [edx + 12]         ;# eax-edx=jnr1-4 
	
	add dword ptr [esp + nb132_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb132_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]     ;# replace jnr with j3 
	lea   edx, [edx + edx*2]	
	
	;# move j coordinates to local temp variables 
	movlps xmm2, [esi + eax*4]
	movlps xmm3, [esi + eax*4 + 12]
	movlps xmm4, [esi + eax*4 + 24]

	movlps xmm5, [esi + ebx*4]
	movlps xmm6, [esi + ebx*4 + 12]
	movlps xmm7, [esi + ebx*4 + 24]

	movhps xmm2, [esi + ecx*4]
	movhps xmm3, [esi + ecx*4 + 12]
	movhps xmm4, [esi + ecx*4 + 24]

	movhps xmm5, [esi + edx*4]
	movhps xmm6, [esi + edx*4 + 12]
	movhps xmm7, [esi + edx*4 + 24]

	;# current state: 	
	;# xmm2= jxOa  jyOa  jxOc  jyOc 
	;# xmm3= jxH1a jyH1a jxH1c jyH1c 
	;# xmm4= jxH2a jyH2a jxH2c jyH2c 
	;# xmm5= jxOb  jyOb  jxOd  jyOd 
	;# xmm6= jxH1b jyH1b jxH1d jyH1d 
	;# xmm7= jxH2b jyH2b jxH2d jyH2d 
	
	movaps xmm0, xmm2
	movaps xmm1, xmm3
	unpcklps xmm0, xmm5	;# xmm0= jxOa  jxOb  jyOa  jyOb 
	unpcklps xmm1, xmm6	;# xmm1= jxH1a jxH1b jyH1a jyH1b 
	unpckhps xmm2, xmm5	;# xmm2= jxOc  jxOd  jyOc  jyOd 
	unpckhps xmm3, xmm6	;# xmm3= jxH1c jxH1d jyH1c jyH1d 
	movaps xmm5, xmm4
	movaps   xmm6, xmm0
	unpcklps xmm4, xmm7	;# xmm4= jxH2a jxH2b jyH2a jyH2b 		
	unpckhps xmm5, xmm7	;# xmm5= jxH2c jxH2d jyH2c jyH2d 
	movaps   xmm7, xmm1
	movlhps  xmm0, xmm2	;# xmm0= jxOa  jxOb  jxOc  jxOd 
	movaps [esp + nb132_jxO], xmm0
	movhlps  xmm2, xmm6	;# xmm2= jyOa  jyOb  jyOc  jyOd 
	movaps [esp + nb132_jyO], xmm2
	movlhps  xmm1, xmm3
	movaps [esp + nb132_jxH1], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [esp + nb132_jyH1], xmm3
	movlhps  xmm4, xmm5
	movaps [esp + nb132_jxH2], xmm4
	movhlps  xmm5, xmm6
	movaps [esp + nb132_jyH2], xmm5

	movss  xmm0, [esi + eax*4 + 8]
	movss  xmm1, [esi + eax*4 + 20]
	movss  xmm2, [esi + eax*4 + 32]

	movss  xmm3, [esi + ecx*4 + 8]
	movss  xmm4, [esi + ecx*4 + 20]
	movss  xmm5, [esi + ecx*4 + 32]

	movhps xmm0, [esi + ebx*4 + 4]
	movhps xmm1, [esi + ebx*4 + 16]
	movhps xmm2, [esi + ebx*4 + 28]
	
	movhps xmm3, [esi + edx*4 + 4]
	movhps xmm4, [esi + edx*4 + 16]
	movhps xmm5, [esi + edx*4 + 28]
	
	shufps xmm0, xmm3, 204  ;# constant 11001100
	shufps xmm1, xmm4, 204  ;# constant 11001100
	shufps xmm2, xmm5, 204  ;# constant 11001100
	movaps [esp + nb132_jzO],  xmm0
	movaps [esp + nb132_jzH1],  xmm1
	movaps [esp + nb132_jzH2],  xmm2

	movaps xmm0, [esp + nb132_ixO]
	movaps xmm1, [esp + nb132_iyO]
	movaps xmm2, [esp + nb132_izO]
	movaps xmm3, [esp + nb132_ixO]
	movaps xmm4, [esp + nb132_iyO]
	movaps xmm5, [esp + nb132_izO]
	subps  xmm0, [esp + nb132_jxO]
	subps  xmm1, [esp + nb132_jyO]
	subps  xmm2, [esp + nb132_jzO]
	subps  xmm3, [esp + nb132_jxH1]
	subps  xmm4, [esp + nb132_jyH1]
	subps  xmm5, [esp + nb132_jzH1]
	movaps [esp + nb132_dxOO], xmm0
	movaps [esp + nb132_dyOO], xmm1
	movaps [esp + nb132_dzOO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb132_dxOH1], xmm3
	movaps [esp + nb132_dyOH1], xmm4
	movaps [esp + nb132_dzOH1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb132_rsqOO], xmm0
	movaps [esp + nb132_rsqOH1], xmm3

	movaps xmm0, [esp + nb132_ixO]
	movaps xmm1, [esp + nb132_iyO]
	movaps xmm2, [esp + nb132_izO]
	movaps xmm3, [esp + nb132_ixH1]
	movaps xmm4, [esp + nb132_iyH1]
	movaps xmm5, [esp + nb132_izH1]
	subps  xmm0, [esp + nb132_jxH2]
	subps  xmm1, [esp + nb132_jyH2]
	subps  xmm2, [esp + nb132_jzH2]
	subps  xmm3, [esp + nb132_jxO]
	subps  xmm4, [esp + nb132_jyO]
	subps  xmm5, [esp + nb132_jzO]
	movaps [esp + nb132_dxOH2], xmm0
	movaps [esp + nb132_dyOH2], xmm1
	movaps [esp + nb132_dzOH2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb132_dxH1O], xmm3
	movaps [esp + nb132_dyH1O], xmm4
	movaps [esp + nb132_dzH1O], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb132_rsqOH2], xmm0
	movaps [esp + nb132_rsqH1O], xmm3

	movaps xmm0, [esp + nb132_ixH1]
	movaps xmm1, [esp + nb132_iyH1]
	movaps xmm2, [esp + nb132_izH1]
	movaps xmm3, [esp + nb132_ixH1]
	movaps xmm4, [esp + nb132_iyH1]
	movaps xmm5, [esp + nb132_izH1]
	subps  xmm0, [esp + nb132_jxH1]
	subps  xmm1, [esp + nb132_jyH1]
	subps  xmm2, [esp + nb132_jzH1]
	subps  xmm3, [esp + nb132_jxH2]
	subps  xmm4, [esp + nb132_jyH2]
	subps  xmm5, [esp + nb132_jzH2]
	movaps [esp + nb132_dxH1H1], xmm0
	movaps [esp + nb132_dyH1H1], xmm1
	movaps [esp + nb132_dzH1H1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb132_dxH1H2], xmm3
	movaps [esp + nb132_dyH1H2], xmm4
	movaps [esp + nb132_dzH1H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb132_rsqH1H1], xmm0
	movaps [esp + nb132_rsqH1H2], xmm3

	movaps xmm0, [esp + nb132_ixH2]
	movaps xmm1, [esp + nb132_iyH2]
	movaps xmm2, [esp + nb132_izH2]
	movaps xmm3, [esp + nb132_ixH2]
	movaps xmm4, [esp + nb132_iyH2]
	movaps xmm5, [esp + nb132_izH2]
	subps  xmm0, [esp + nb132_jxO]
	subps  xmm1, [esp + nb132_jyO]
	subps  xmm2, [esp + nb132_jzO]
	subps  xmm3, [esp + nb132_jxH1]
	subps  xmm4, [esp + nb132_jyH1]
	subps  xmm5, [esp + nb132_jzH1]
	movaps [esp + nb132_dxH2O], xmm0
	movaps [esp + nb132_dyH2O], xmm1
	movaps [esp + nb132_dzH2O], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb132_dxH2H1], xmm3
	movaps [esp + nb132_dyH2H1], xmm4
	movaps [esp + nb132_dzH2H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [esp + nb132_rsqH2O], xmm0
	movaps [esp + nb132_rsqH2H1], xmm4

	movaps xmm0, [esp + nb132_ixH2]
	movaps xmm1, [esp + nb132_iyH2]
	movaps xmm2, [esp + nb132_izH2]
	subps  xmm0, [esp + nb132_jxH2]
	subps  xmm1, [esp + nb132_jyH2]
	subps  xmm2, [esp + nb132_jzH2]
	movaps [esp + nb132_dxH2H2], xmm0
	movaps [esp + nb132_dyH2H2], xmm1
	movaps [esp + nb132_dzH2H2], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [esp + nb132_rsqH2H2], xmm0
		
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb132_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb132_half] ;# rinvH2H2 
	mulps   xmm7, [esp + nb132_half] ;# rinvH2H1 
	movaps  [esp + nb132_rinvH2H2], xmm3
	movaps  [esp + nb132_rinvH2H1], xmm7
	
	rsqrtps xmm1, [esp + nb132_rsqOO]
	rsqrtps xmm5, [esp + nb132_rsqOH1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb132_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb132_rsqOO]
	mulps   xmm5, [esp + nb132_rsqOH1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb132_half] 
	mulps   xmm7, [esp + nb132_half]
	movaps  [esp + nb132_rinvOO], xmm3
	movaps  [esp + nb132_rinvOH1], xmm7
	
	rsqrtps xmm1, [esp + nb132_rsqOH2]
	rsqrtps xmm5, [esp + nb132_rsqH1O]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb132_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb132_rsqOH2]
	mulps   xmm5, [esp + nb132_rsqH1O]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb132_half] 
	mulps   xmm7, [esp + nb132_half]
	movaps  [esp + nb132_rinvOH2], xmm3
	movaps  [esp + nb132_rinvH1O], xmm7
	
	rsqrtps xmm1, [esp + nb132_rsqH1H1]
	rsqrtps xmm5, [esp + nb132_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb132_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb132_rsqH1H1]
	mulps   xmm5, [esp + nb132_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb132_half] 
	mulps   xmm7, [esp + nb132_half]
	movaps  [esp + nb132_rinvH1H1], xmm3
	movaps  [esp + nb132_rinvH1H2], xmm7
	
	rsqrtps xmm1, [esp + nb132_rsqH2O]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + nb132_three]
	mulps   xmm1, [esp + nb132_rsqH2O]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [esp + nb132_half] 
	movaps  [esp + nb132_rinvH2O], xmm3

	;# start with OO interaction - first the table LJ part
	movaps xmm0, [esp + nb132_rinvOO]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb132_rsqOO] ;# xmm1=r 
	mulps  xmm1, [esp + nb132_tsc]
		
	movhlps xmm2, xmm1
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 3
    pslld mm7, 3

    movd mm0, eax
    movd mm1, ebx
    movd mm2, ecx
    movd mm3, edx

    movd eax, mm6
    psrlq mm6, 32
    movd ecx, mm7
    psrlq mm7, 32
    movd ebx, mm6
    movd edx, mm7

    mov  esi, [ebp + nb132_VFtab]
	
    ;# dispersion 
    movlps xmm5, [esi + eax*4]
    movlps xmm7, [esi + ecx*4]
    movhps xmm5, [esi + ebx*4]
    movhps xmm7, [esi + edx*4] ;# got half table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# constant 10001000
    shufps xmm5, xmm7, 221  ;# constant 11011101

    movlps xmm7, [esi + eax*4 + 8]
    movlps xmm3, [esi + ecx*4 + 8]
    movhps xmm7, [esi + ebx*4 + 8]
    movhps xmm3, [esi + edx*4 + 8] ;# other half of table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# constant 10001000
    shufps xmm7, xmm3, 221  ;# constant 11011101
    ;# dispersion table ready, in xmm4-xmm7 
    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    mulps  xmm7, [esp + nb132_two]       ;# two*Heps2 
    addps  xmm7, xmm6
    addps  xmm7, xmm5 ;# xmm7=FF 
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 

    movaps xmm4, [esp + nb132_c6]
    mulps  xmm7, xmm4    ;# fijD 
    mulps  xmm5, xmm4    ;# Vvdw6 
    movaps  [esp + nb132_fstmp], xmm7

    addps  xmm5, [esp + nb132_Vvdwtot]
    movaps [esp + nb132_Vvdwtot], xmm5

    ;# repulsion 
    movlps xmm5, [esi + eax*4 + 16]
    movlps xmm7, [esi + ecx*4 + 16]
    movhps xmm5, [esi + ebx*4 + 16]
    movhps xmm7, [esi + edx*4 + 16] ;# got half table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# constant 10001000
    shufps xmm5, xmm7, 221  ;# constant 11011101

    movlps xmm7, [esi + eax*4 + 24]
    movlps xmm3, [esi + ecx*4 + 24]
    movhps xmm7, [esi + ebx*4 + 24]
    movhps xmm3, [esi + edx*4 + 24] ;# other half of table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# constant 10001000
    shufps xmm7, xmm3, 221  ;# constant 11011101
    ;# repulsion table ready, in xmm4-xmm7 
    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    mulps  xmm7, [esp + nb132_two]       ;# two*Heps2 
    addps  xmm7, xmm6
    addps  xmm7, xmm5 ;# xmm7=FF 
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
 
    movaps xmm4, [esp + nb132_c12]
    mulps  xmm7, xmm4 ;# fijR 
    mulps  xmm5, xmm4 ;# Vvdw12 
    addps  xmm7, [esp + nb132_fstmp] 
    mulps xmm7, [esp + nb132_tsc]

    addps  xmm5, [esp + nb132_Vvdwtot]
    movaps [esp + nb132_Vvdwtot], xmm5

	movaps xmm0, [esp + nb132_rinvOO]
	movaps xmm6, xmm0	
	mulps  xmm6, [esp + nb132_qqOO] 
	movaps xmm2, xmm6
	
	addps  xmm6, [esp + nb132_vctot] ;# local vctot summation variable 
	mulps  xmm2, xmm0
	subps  xmm2, xmm7
	mulps  xmm2, xmm0
	
	movaps xmm0, xmm2
	movaps xmm1, xmm2
	
    movd  eax, mm0
    movd  ebx, mm1
    movd  ecx, mm2
    movd  edx, mm3

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + nb132_dxOO]
	mulps xmm1, [esp + nb132_dyOO]
	mulps xmm2, [esp + nb132_dzOO]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb132_fixO]
	addps xmm1, [esp + nb132_fiyO]
	addps xmm2, [esp + nb132_fizO]
	movaps [esp + nb132_fjxO], xmm3
	movaps [esp + nb132_fjyO], xmm4
	movaps [esp + nb132_fjzO], xmm5
	movaps [esp + nb132_fixO], xmm0
	movaps [esp + nb132_fiyO], xmm1
	movaps [esp + nb132_fizO], xmm2

	;# O-H1 interaction 
	movaps xmm0, [esp + nb132_rinvOH1]
	movaps xmm4, xmm0	;# xmm7=rinv 
	mulps  xmm0, xmm0   ;# rinvsq
	mulps  xmm4, [esp + nb132_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulps  xmm0, xmm4

	addps  xmm6, xmm4	;# add to local vctot 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + nb132_dxOH1]
	mulps xmm1, [esp + nb132_dyOH1]
	mulps xmm2, [esp + nb132_dzOH1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb132_fixO]
	addps xmm1, [esp + nb132_fiyO]
	addps xmm2, [esp + nb132_fizO]
	movaps [esp + nb132_fjxH1], xmm3
	movaps [esp + nb132_fjyH1], xmm4
	movaps [esp + nb132_fjzH1], xmm5
	movaps [esp + nb132_fixO], xmm0
	movaps [esp + nb132_fiyO], xmm1
	movaps [esp + nb132_fizO], xmm2

	;# O-H2 interaction  
	movaps xmm0, [esp + nb132_rinvOH2]
	movaps xmm4, xmm0	;# xmm7=rinv 
	mulps  xmm0, xmm0   ;# rinvsq
	mulps  xmm4, [esp + nb132_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulps  xmm0, xmm4

	addps  xmm6, xmm4	;# add to local vctot 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + nb132_dxOH2]
	mulps xmm1, [esp + nb132_dyOH2]
	mulps xmm2, [esp + nb132_dzOH2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb132_fixO]
	addps xmm1, [esp + nb132_fiyO]
	addps xmm2, [esp + nb132_fizO]
	movaps [esp + nb132_fjxH2], xmm3
	movaps [esp + nb132_fjyH2], xmm4
	movaps [esp + nb132_fjzH2], xmm5
	movaps [esp + nb132_fixO], xmm0
	movaps [esp + nb132_fiyO], xmm1
	movaps [esp + nb132_fizO], xmm2

	;# H1-O interaction 
	movaps xmm0, [esp + nb132_rinvH1O]
	movaps xmm4, xmm0	;# xmm7=rinv 
	mulps  xmm0, xmm0   ;# rinvsq
	mulps  xmm4, [esp + nb132_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulps  xmm0, xmm4

	addps  xmm6, xmm4	;# add to local vctot 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	movaps xmm3, [esp + nb132_fjxO]
	movaps xmm4, [esp + nb132_fjyO]
	movaps xmm5, [esp + nb132_fjzO]
	mulps xmm0, [esp + nb132_dxH1O]
	mulps xmm1, [esp + nb132_dyH1O]
	mulps xmm2, [esp + nb132_dzH1O]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb132_fixH1]
	addps xmm1, [esp + nb132_fiyH1]
	addps xmm2, [esp + nb132_fizH1]
	movaps [esp + nb132_fjxO], xmm3
	movaps [esp + nb132_fjyO], xmm4
	movaps [esp + nb132_fjzO], xmm5
	movaps [esp + nb132_fixH1], xmm0
	movaps [esp + nb132_fiyH1], xmm1
	movaps [esp + nb132_fizH1], xmm2

	;# H1-H1 interaction 
	movaps xmm0, [esp + nb132_rinvH1H1]
	movaps xmm4, xmm0	;# xmm7=rinv 
	mulps  xmm0, xmm0   ;# rinvsq
	mulps  xmm4, [esp + nb132_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulps  xmm0, xmm4

	addps  xmm6, xmm4	;# add to local vctot 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	movaps xmm3, [esp + nb132_fjxH1]
	movaps xmm4, [esp + nb132_fjyH1]
	movaps xmm5, [esp + nb132_fjzH1]
	mulps xmm0, [esp + nb132_dxH1H1]
	mulps xmm1, [esp + nb132_dyH1H1]
	mulps xmm2, [esp + nb132_dzH1H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb132_fixH1]
	addps xmm1, [esp + nb132_fiyH1]
	addps xmm2, [esp + nb132_fizH1]
	movaps [esp + nb132_fjxH1], xmm3
	movaps [esp + nb132_fjyH1], xmm4
	movaps [esp + nb132_fjzH1], xmm5
	movaps [esp + nb132_fixH1], xmm0
	movaps [esp + nb132_fiyH1], xmm1
	movaps [esp + nb132_fizH1], xmm2

	;# H1-H2 interaction 
	movaps xmm0, [esp + nb132_rinvH1H2]
	movaps xmm4, xmm0	;# xmm7=rinv 
	mulps  xmm0, xmm0   ;# rinvsq
	mulps  xmm4, [esp + nb132_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulps  xmm0, xmm4

	addps  xmm6, xmm4	;# add to local vctot 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	movaps xmm3, [esp + nb132_fjxH2]
	movaps xmm4, [esp + nb132_fjyH2]
	movaps xmm5, [esp + nb132_fjzH2]
	mulps xmm0, [esp + nb132_dxH1H2]
	mulps xmm1, [esp + nb132_dyH1H2]
	mulps xmm2, [esp + nb132_dzH1H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb132_fixH1]
	addps xmm1, [esp + nb132_fiyH1]
	addps xmm2, [esp + nb132_fizH1]
	movaps [esp + nb132_fjxH2], xmm3
	movaps [esp + nb132_fjyH2], xmm4
	movaps [esp + nb132_fjzH2], xmm5
	movaps [esp + nb132_fixH1], xmm0
	movaps [esp + nb132_fiyH1], xmm1
	movaps [esp + nb132_fizH1], xmm2

	;# H2-O interaction 
	movaps xmm0, [esp + nb132_rinvH2O]
	movaps xmm4, xmm0	;# xmm7=rinv 
	mulps  xmm0, xmm0   ;# rinvsq
	mulps  xmm4, [esp + nb132_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulps  xmm0, xmm4

	addps  xmm6, xmm4	;# add to local vctot 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	movaps xmm3, [esp + nb132_fjxO]
	movaps xmm4, [esp + nb132_fjyO]
	movaps xmm5, [esp + nb132_fjzO]
	mulps xmm0, [esp + nb132_dxH2O]
	mulps xmm1, [esp + nb132_dyH2O]
	mulps xmm2, [esp + nb132_dzH2O]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb132_fixH2]
	addps xmm1, [esp + nb132_fiyH2]
	addps xmm2, [esp + nb132_fizH2]
	movaps [esp + nb132_fjxO], xmm3
	movaps [esp + nb132_fjyO], xmm4
	movaps [esp + nb132_fjzO], xmm5
	movaps [esp + nb132_fixH2], xmm0
	movaps [esp + nb132_fiyH2], xmm1
	movaps [esp + nb132_fizH2], xmm2

	;# H2-H1 interaction 
	movaps xmm0, [esp + nb132_rinvH2H1]
	movaps xmm4, xmm0	;# xmm7=rinv 
	mulps  xmm0, xmm0   ;# rinvsq
	mulps  xmm4, [esp + nb132_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulps  xmm0, xmm4

	addps  xmm6, xmm4	;# add to local vctot 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	movaps xmm3, [esp + nb132_fjxH1]
	movaps xmm4, [esp + nb132_fjyH1]
	movaps xmm5, [esp + nb132_fjzH1]
	mulps xmm0, [esp + nb132_dxH2H1]
	mulps xmm1, [esp + nb132_dyH2H1]
	mulps xmm2, [esp + nb132_dzH2H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb132_fixH2]
	addps xmm1, [esp + nb132_fiyH2]
	addps xmm2, [esp + nb132_fizH2]
	movaps [esp + nb132_fjxH1], xmm3
	movaps [esp + nb132_fjyH1], xmm4
	movaps [esp + nb132_fjzH1], xmm5
	movaps [esp + nb132_fixH2], xmm0
	movaps [esp + nb132_fiyH2], xmm1
	movaps [esp + nb132_fizH2], xmm2

	;# H2-H2 interaction 
	movaps xmm0, [esp + nb132_rinvH2H2]
	movaps xmm4, xmm0	;# xmm7=rinv 
	mulps  xmm0, xmm0   ;# rinvsq
	mulps  xmm4, [esp + nb132_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulps  xmm0, xmm4

	addps  xmm6, xmm4	;# add to local vctot 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	movaps xmm1, xmm0
	movaps [esp + nb132_vctot], xmm6
	movaps xmm2, xmm0
	
	movaps xmm3, [esp + nb132_fjxH2]
	movaps xmm4, [esp + nb132_fjyH2]
	movaps xmm5, [esp + nb132_fjzH2]
	mulps xmm0, [esp + nb132_dxH2H2]
	mulps xmm1, [esp + nb132_dyH2H2]
	mulps xmm2, [esp + nb132_dzH2H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb132_fixH2]
	addps xmm1, [esp + nb132_fiyH2]
	addps xmm2, [esp + nb132_fizH2]
	movaps [esp + nb132_fjxH2], xmm3
	movaps [esp + nb132_fjyH2], xmm4
	movaps [esp + nb132_fjzH2], xmm5
	movaps [esp + nb132_fixH2], xmm0
	movaps [esp + nb132_fiyH2], xmm1
	movaps [esp + nb132_fizH2], xmm2

	mov edi, [ebp + nb132_faction]
		
	;# Did all interactions - now update j forces 
	;# At this stage forces are still on the stack, in positions:
	;# fjxO, fjyO, fjzO, ... , fjzH2.
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
	
		
	;# 4 j waters with three atoms each - first do Oxygen X & Y forces for 4 j particles 
	movaps xmm0, [esp + nb132_fjxO] ;# xmm0= fjxOa  fjxOb  fjxOc  fjxOd 
	movaps xmm2, [esp + nb132_fjyO] ;# xmm1= fjyOa  fjyOb  fjyOc  fjyOd
	movlps xmm3, [edi + eax*4]
	movlps xmm4, [edi + ecx*4]
	movaps xmm1, xmm0
	unpcklps xmm0, xmm2    	   ;# xmm0= fjxOa  fjyOa  fjxOb  fjyOb
	unpckhps xmm1, xmm2        ;# xmm1= fjxOc  fjyOc  fjxOd  fjyOd
	movhps xmm3, [edi + ebx*4]
	movhps xmm4, [edi + edx*4]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	movlps [edi + eax*4], xmm3
	movlps [edi + ecx*4], xmm4
	movhps [edi + ebx*4], xmm3
	movhps [edi + edx*4], xmm4
	
	;# Oxygen Z & first hydrogen X forces for 4 j particles 
	movaps xmm0, [esp + nb132_fjzO]  ;# xmm0= fjzOa   fjzOb   fjzOc   fjzOd 
	movaps xmm2, [esp + nb132_fjxH1] ;# xmm1= fjxH1a  fjxH1b  fjxH1c  fjxH1d
	movlps xmm3, [edi + eax*4 + 8]
	movlps xmm4, [edi + ecx*4 + 8]
	movaps xmm1, xmm0
	unpcklps xmm0, xmm2    	   ;# xmm0= fjzOa  fjxH1a  fjzOb  fjxH1b
	unpckhps xmm1, xmm2        ;# xmm1= fjzOc  fjxH1c  fjzOd  fjxH1d
	movhps xmm3, [edi + ebx*4 + 8]
	movhps xmm4, [edi + edx*4 + 8]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	movlps [edi + eax*4 + 8], xmm3
	movlps [edi + ecx*4 + 8], xmm4
	movhps [edi + ebx*4 + 8], xmm3
	movhps [edi + edx*4 + 8], xmm4

	
	;# First hydrogen Y & Z forces for 4 j particles 
	movaps xmm0, [esp + nb132_fjyH1]  ;# xmm0= fjyH1a  fjyH1b  fjyH1c  fjyH1d 
	movaps xmm2, [esp + nb132_fjzH1] ;# xmm1= fjzH1a  fjzH1b  fjzH1c  fjzH1d
	movlps xmm3, [edi + eax*4 + 16]
	movlps xmm4, [edi + ecx*4 + 16]
	movaps xmm1, xmm0
	unpcklps xmm0, xmm2		;# xmm0= fjyH1a  fjzH1a  fjyH1b  fjzH1b
	unpckhps xmm1, xmm2		;# xmm1= fjyH1c  fjzH1c  fjyH1d  fjzH1d
	movhps xmm3, [edi + ebx*4 + 16]
	movhps xmm4, [edi + edx*4 + 16]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	movlps [edi + eax*4 + 16], xmm3
	movlps [edi + ecx*4 + 16], xmm4
	movhps [edi + ebx*4 + 16], xmm3
	movhps [edi + edx*4 + 16], xmm4

	
	;# Second hydrogen X & Y forces for 4 j particles 
	movaps xmm0, [esp + nb132_fjxH2]  ;# xmm0= fjxH2a  fjxH2b  fjxH2c  fjxH2d 
	movaps xmm2, [esp + nb132_fjyH2] ;# xmm1= fjyH2a  fjyH2b  fjyH2c  fjyH2d
	movlps xmm3, [edi + eax*4 + 24]
	movlps xmm4, [edi + ecx*4 + 24]
	movaps xmm1, xmm0
	unpcklps xmm0, xmm2		;# xmm0= fjxH2a  fjyH2a  fjxH2b  fjyH2b
	unpckhps xmm1, xmm2		;# xmm1= fjxH2c  fjyH2c  fjxH2d  fjyH2d
	movhps xmm3, [edi + ebx*4 + 24]
	movhps xmm4, [edi + edx*4 + 24]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	movlps [edi + eax*4 + 24], xmm3
	movlps [edi + ecx*4 + 24], xmm4
	movhps [edi + ebx*4 + 24], xmm3
	movhps [edi + edx*4 + 24], xmm4

	
	;# Second hydrogen Z forces for 4 j particles 
	;# Just load the four Z coords into one reg. each
	movss xmm4, [edi + eax*4 + 32]
	movss xmm5, [edi + ebx*4 + 32]
	movss xmm6, [edi + ecx*4 + 32]
	movss xmm7, [edi + edx*4 + 32]
	;# add what we have on the stack
	addss xmm4, [esp + nb132_fjzH2] 
	addss xmm5, [esp + nb132_fjzH2b] 
	addss xmm6, [esp + nb132_fjzH2c] 
	addss xmm7, [esp + nb132_fjzH2d]
	;# store back
	movss [edi + eax*4 + 32], xmm4
	movss [edi + ebx*4 + 32], xmm5
	movss [edi + ecx*4 + 32], xmm6
	movss [edi + edx*4 + 32], xmm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb132_innerk],  4
	jl    .nb132_single_check
	jmp   .nb132_unroll_loop
.nb132_single_check:
	add dword ptr [esp + nb132_innerk],  4
	jnz   .nb132_single_loop
	jmp   .nb132_updateouterdata
.nb132_single_loop:
	mov   edx, [esp + nb132_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb132_innerjjnr],  4	

	mov esi, [ebp + nb132_pos]
	lea   eax, [eax + eax*2]  

	;# fetch j coordinates 
	xorps xmm3, xmm3
	xorps xmm4, xmm4
	xorps xmm5, xmm5
	
	movss xmm3, [esi + eax*4]		;# jxO  -  -  -
	movss xmm4, [esi + eax*4 + 4]		;# jyO  -  -  -
	movss xmm5, [esi + eax*4 + 8]		;# jzO  -  -  -  

	movlps xmm6, [esi + eax*4 + 12]		;# xmm6 = jxH1 jyH1   -    -
	movss  xmm7, [esi + eax*4 + 20]		;# xmm7 = jzH1   -    -    - 
	movhps xmm6, [esi + eax*4 + 24]		;# xmm6 = jxH1 jyH1 jxH2 jyH2
	movss  xmm2, [esi + eax*4 + 32]		;# xmm2 = jzH2   -    -    -
	
	;# have all coords, time for some shuffling.

	shufps xmm6, xmm6, 216 ;# constant 11011000	;# xmm6 = jxH1 jxH2 jyH1 jyH2 
	unpcklps xmm7, xmm2			;# xmm7 = jzH1 jzH2   -    -
	movaps  xmm0, [esp + nb132_ixO]     
	movaps  xmm1, [esp + nb132_iyO]
	movaps  xmm2, [esp + nb132_izO]	
	movlhps xmm3, xmm6			;# xmm3 = jxO   0   jxH1 jxH2 
	shufps  xmm4, xmm6, 228 ;# constant 11100100	;# xmm4 = jyO   0   jyH1 jyH2 
	shufps  xmm5, xmm7, 68  ;# constant 01000100	;# xmm5 = jzO   0   jzH1 jzH2
	
	;# store all j coordinates in jO  
	movaps [esp + nb132_jxO], xmm3
	movaps [esp + nb132_jyO], xmm4
	movaps [esp + nb132_jzO], xmm5
	subps  xmm0, xmm3
	subps  xmm1, xmm4
	subps  xmm2, xmm5
	movaps [esp + nb132_dxOO], xmm0
	movaps [esp + nb132_dyOO], xmm1
	movaps [esp + nb132_dzOO], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 
	movaps [esp + nb132_rsqOO], xmm0
	
	movaps xmm6, xmm0
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + nb132_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [esp + nb132_half] ;# rinv iO - j water in xmm3
	movaps  [esp + nb132_rinvOO], xmm3

	movaps  xmm0, xmm3  ;# rinv

	xorps   xmm4, xmm4
	;# fetch charges to xmm4 (temporary) 
	movss   xmm4, [esp + nb132_qqOO]
	movhps  xmm4, [esp + nb132_qqOH]

	mulps  xmm3, xmm4   ;# vcoul
	movaps xmm6, xmm3   
	
	mulps xmm3, xmm0
	movaps [esp + nb132_fstmp], xmm3 ;# save it
	addps  xmm6, [esp + nb132_vctot]
    movaps [esp + nb132_vctot], xmm6

	movaps xmm0, [esp + nb132_rinvOO]
	movss xmm1, xmm0
	mulss  xmm1, [esp + nb132_rsqOO] ;# xmm1=r 
	mulss  xmm1, [esp + nb132_tsc]
		
    cvttps2pi mm6, xmm1
    cvtpi2ps xmm3, mm6
	subss    xmm1, xmm3	;# xmm1=eps 
    movss xmm2, xmm1
    mulss  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 3
	
    movd mm0, eax
	
    mov  esi, [ebp + nb132_VFtab]
    movd eax, mm6
	
    ;# dispersion 
    movlps xmm5, [esi + eax*4]
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# constant 10001000
    shufps xmm5, xmm7, 221  ;# constant 11011101

    movlps xmm7, [esi + eax*4 + 8]
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# constant 10001000
    shufps xmm7, xmm3, 221  ;# constant 11011101
    ;# dispersion table ready, in xmm4-xmm7 
    mulss  xmm6, xmm1       ;# xmm6=Geps 
    mulss  xmm7, xmm2       ;# xmm7=Heps2 
    addss  xmm5, xmm6
    addss  xmm5, xmm7       ;# xmm5=Fp 
    mulss  xmm7, [esp + nb132_two]       ;# two*Heps2 
    addss  xmm7, xmm6
    addss  xmm7, xmm5 ;# xmm7=FF 
    mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
    addss  xmm5, xmm4 ;# xmm5=VV 

    movss xmm4, [esp + nb132_c6]
    mulss  xmm7, xmm4    ;# fijD 
    mulss  xmm5, xmm4    ;# Vvdw6 
	movss  xmm3, [esp + nb132_fstmp]
	mulps  xmm7, [esp + nb132_tsc]
	subss  xmm3, xmm7
	movss  [esp + nb132_fstmp], xmm3 

    addss  xmm5, [esp + nb132_Vvdwtot]
    movss [esp + nb132_Vvdwtot], xmm5

    ;# repulsion 
    movlps xmm5, [esi + eax*4 + 16]
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# constant 10001000
    shufps xmm5, xmm7, 221  ;# constant 11011101

    movlps xmm7, [esi + eax*4 + 24]
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# constant 10001000
    shufps xmm7, xmm3, 221  ;# constant 11011101
    ;# table ready, in xmm4-xmm7 
    mulss  xmm6, xmm1       ;# xmm6=Geps 
    mulss  xmm7, xmm2       ;# xmm7=Heps2 
    addss  xmm5, xmm6
    addss  xmm5, xmm7       ;# xmm5=Fp 
    mulss  xmm7, [esp + nb132_two]       ;# two*Heps2 
    addss  xmm7, xmm6
    addss  xmm7, xmm5 ;# xmm7=FF 
    mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
    addss  xmm5, xmm4 ;# xmm5=VV 

    movss xmm4, [esp + nb132_c12]
    mulss  xmm7, xmm4 ;# fijR 
    mulss  xmm5, xmm4 ;# Vvdw12 
	movaps xmm3, [esp + nb132_fstmp]
	mulss  xmm7, [esp + nb132_tsc]
	subss  xmm3, xmm7

    addss  xmm5, [esp + nb132_Vvdwtot]
    movss [esp + nb132_Vvdwtot], xmm5

	mulps  xmm0, xmm3
	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	
	mulps   xmm0, [esp + nb132_dxOO]
	mulps   xmm1, [esp + nb132_dyOO]
	mulps   xmm2, [esp + nb132_dzOO]
	
	movd eax, mm0

	;# initial update for j forces 
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + nb132_fjxO], xmm3
	movaps  [esp + nb132_fjyO], xmm4
	movaps  [esp + nb132_fjzO], xmm5
	addps   xmm0, [esp + nb132_fixO]
	addps   xmm1, [esp + nb132_fiyO]
	addps   xmm2, [esp + nb132_fizO]
	movaps  [esp + nb132_fixO], xmm0
	movaps  [esp + nb132_fiyO], xmm1
	movaps  [esp + nb132_fizO], xmm2

	
	;# done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [esp + nb132_ixH1]
	movaps  xmm1, [esp + nb132_iyH1]
	movaps  xmm2, [esp + nb132_izH1]	
	movaps  xmm3, [esp + nb132_ixH2] 
	movaps  xmm4, [esp + nb132_iyH2] 
	movaps  xmm5, [esp + nb132_izH2] 
	subps   xmm0, [esp + nb132_jxO]
	subps   xmm1, [esp + nb132_jyO]
	subps   xmm2, [esp + nb132_jzO]
	subps   xmm3, [esp + nb132_jxO]
	subps   xmm4, [esp + nb132_jyO]
	subps   xmm5, [esp + nb132_jzO]
	movaps [esp + nb132_dxH1O], xmm0
	movaps [esp + nb132_dyH1O], xmm1
	movaps [esp + nb132_dzH1O], xmm2
	movaps [esp + nb132_dxH2O], xmm3
	movaps [esp + nb132_dyH2O], xmm4
	movaps [esp + nb132_dzH2O], xmm5
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
	movaps  xmm3, [esp + nb132_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb132_half] ;# rinv H1 - j water 
	mulps   xmm7, [esp + nb132_half] ;# rinv H2 - j water  

	;# assemble charges in xmm6 
	xorps   xmm6, xmm6
	movss   xmm6, [esp + nb132_qqOH]
	movhps  xmm6, [esp + nb132_qqHH]
	
	movaps  xmm0, xmm3
	movaps  xmm1, xmm7
	mulps   xmm3, xmm6 ;# vcoul
	mulps   xmm7, xmm6 ;# vcoul
	mulps   xmm0, xmm0
	mulps   xmm1, xmm1
	mulps   xmm0, xmm3
	addps   xmm3, xmm7
	mulps   xmm7, xmm1
	addps   xmm3, [esp + nb132_vctot]
	movaps  [esp + nb132_vctot], xmm3

	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	mulps   xmm0, [esp + nb132_dxH1O]
	mulps   xmm1, [esp + nb132_dyH1O]
	mulps   xmm2, [esp + nb132_dzH1O]
	;# update forces H1 - j water 
	movaps  xmm3, [esp + nb132_fjxO]
	movaps  xmm4, [esp + nb132_fjyO]
	movaps  xmm5, [esp + nb132_fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + nb132_fjxO], xmm3
	movaps  [esp + nb132_fjyO], xmm4
	movaps  [esp + nb132_fjzO], xmm5
	addps   xmm0, [esp + nb132_fixH1]
	addps   xmm1, [esp + nb132_fiyH1]
	addps   xmm2, [esp + nb132_fizH1]
	movaps  [esp + nb132_fixH1], xmm0
	movaps  [esp + nb132_fiyH1], xmm1
	movaps  [esp + nb132_fizH1], xmm2
	;# do forces H2 - j water 
	movaps xmm0, xmm7
	movaps xmm1, xmm7
	movaps xmm2, xmm7
	mulps   xmm0, [esp + nb132_dxH2O]
	mulps   xmm1, [esp + nb132_dyH2O]
	mulps   xmm2, [esp + nb132_dzH2O]
	movaps  xmm3, [esp + nb132_fjxO]
	movaps  xmm4, [esp + nb132_fjyO]
	movaps  xmm5, [esp + nb132_fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	mov     esi, [ebp + nb132_faction]
	movaps  [esp + nb132_fjxO], xmm3
	movaps  [esp + nb132_fjyO], xmm4
	movaps  [esp + nb132_fjzO], xmm5
	addps   xmm0, [esp + nb132_fixH2]
	addps   xmm1, [esp + nb132_fiyH2]
	addps   xmm2, [esp + nb132_fizH2]
	movaps  [esp + nb132_fixH2], xmm0
	movaps  [esp + nb132_fiyH2], xmm1
	movaps  [esp + nb132_fizH2], xmm2

	;# update j water forces from local variables 
	movlps  xmm0, [esi + eax*4]
	movlps  xmm1, [esi + eax*4 + 12]
	movhps  xmm1, [esi + eax*4 + 24]
	movaps  xmm3, [esp + nb132_fjxO]
	movaps  xmm4, [esp + nb132_fjyO]
	movaps  xmm5, [esp + nb132_fjzO]
	movaps  xmm6, xmm5
	movaps  xmm7, xmm5
	shufps  xmm6, xmm6, 2 ;# constant 00000010
	shufps  xmm7, xmm7, 3 ;# constant 00000011
	addss   xmm5, [esi + eax*4 + 8]
	addss   xmm6, [esi + eax*4 + 20]
	addss   xmm7, [esi + eax*4 + 32]
	movss   [esi + eax*4 + 8], xmm5
	movss   [esi + eax*4 + 20], xmm6
	movss   [esi + eax*4 + 32], xmm7
	movaps   xmm5, xmm3
	unpcklps xmm3, xmm4
	unpckhps xmm5, xmm4
	addps    xmm0, xmm3
	addps    xmm1, xmm5
	movlps  [esi + eax*4], xmm0 
	movlps  [esi + eax*4 + 12], xmm1 
	movhps  [esi + eax*4 + 24], xmm1 
	
	dec dword ptr [esp + nb132_innerk]
	jz    .nb132_updateouterdata
	jmp   .nb132_single_loop
.nb132_updateouterdata:
	mov   ecx, [esp + nb132_ii3]
	mov   edi, [ebp + nb132_faction]
	mov   esi, [ebp + nb132_fshift]
	mov   edx, [esp + nb132_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb132_fixO]
	movaps xmm1, [esp + nb132_fiyO] 
	movaps xmm2, [esp + nb132_fizO]

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
	movaps xmm0, [esp + nb132_fixH1]
	movaps xmm1, [esp + nb132_fiyH1]
	movaps xmm2, [esp + nb132_fizH1]

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
	movaps xmm0, [esp + nb132_fixH2]
	movaps xmm1, [esp + nb132_fiyH2]
	movaps xmm2, [esp + nb132_fizH2]

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

	;# increment fshift force  
	movlps  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [esi + edx*4],    xmm3
	movss  [esi + edx*4 + 8], xmm4

	;# get n from stack
	mov esi, [esp + nb132_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb132_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb132_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb132_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [esp + nb132_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb132_Vvdw]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb132_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb132_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb132_n], esi
        jmp .nb132_outer
.nb132_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb132_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb132_end
        ;# non-zero, do one more workunit
        jmp   .nb132_threadloop
.nb132_end:
	emms

	mov eax, [esp + nb132_nouter]
	mov ebx, [esp + nb132_ninner]
	mov ecx, [ebp + nb132_outeriter]
	mov edx, [ebp + nb132_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb132_salign]
	add esp, eax
	add esp, 1560
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret


.globl nb_kernel132nf_ia32_sse
.globl _nb_kernel132nf_ia32_sse
nb_kernel132nf_ia32_sse:	
_nb_kernel132nf_ia32_sse:	
.equiv          nb132nf_p_nri,            8
.equiv          nb132nf_iinr,             12
.equiv          nb132nf_jindex,           16
.equiv          nb132nf_jjnr,             20
.equiv          nb132nf_shift,            24
.equiv          nb132nf_shiftvec,         28
.equiv          nb132nf_fshift,           32
.equiv          nb132nf_gid,              36
.equiv          nb132nf_pos,              40
.equiv          nb132nf_faction,          44
.equiv          nb132nf_charge,           48
.equiv          nb132nf_p_facel,          52
.equiv          nb132nf_argkrf,           56
.equiv          nb132nf_argcrf,           60
.equiv          nb132nf_Vc,               64
.equiv          nb132nf_type,             68
.equiv          nb132nf_p_ntype,          72
.equiv          nb132nf_vdwparam,         76
.equiv          nb132nf_Vvdw,             80
.equiv          nb132nf_p_tabscale,       84
.equiv          nb132nf_VFtab,            88
.equiv          nb132nf_invsqrta,         92
.equiv          nb132nf_dvda,             96
.equiv          nb132nf_p_gbtabscale,     100
.equiv          nb132nf_GBtab,            104
.equiv          nb132nf_p_nthreads,       108
.equiv          nb132nf_count,            112
.equiv          nb132nf_mtx,              116
.equiv          nb132nf_outeriter,        120
.equiv          nb132nf_inneriter,        124
.equiv          nb132nf_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb132nf_ixO,              0
.equiv          nb132nf_iyO,              16
.equiv          nb132nf_izO,              32
.equiv          nb132nf_ixH1,             48
.equiv          nb132nf_iyH1,             64
.equiv          nb132nf_izH1,             80
.equiv          nb132nf_ixH2,             96
.equiv          nb132nf_iyH2,             112
.equiv          nb132nf_izH2,             128
.equiv          nb132nf_jxO,              144
.equiv          nb132nf_jyO,              160
.equiv          nb132nf_jzO,              176
.equiv          nb132nf_jxH1,             192
.equiv          nb132nf_jyH1,             208
.equiv          nb132nf_jzH1,             224
.equiv          nb132nf_jxH2,             240
.equiv          nb132nf_jyH2,             256
.equiv          nb132nf_jzH2,             272
.equiv          nb132nf_qqOO,             288
.equiv          nb132nf_qqOH,             304
.equiv          nb132nf_qqHH,             320
.equiv          nb132nf_c6,               336
.equiv          nb132nf_c12,              352
.equiv          nb132nf_tsc,              368
.equiv          nb132nf_vctot,            384
.equiv          nb132nf_Vvdwtot,          400
.equiv          nb132nf_half,             416
.equiv          nb132nf_three,            432
.equiv          nb132nf_rsqOO,            448
.equiv          nb132nf_rsqOH1,           464
.equiv          nb132nf_rsqOH2,           480
.equiv          nb132nf_rsqH1O,           496
.equiv          nb132nf_rsqH1H1,          512
.equiv          nb132nf_rsqH1H2,          528
.equiv          nb132nf_rsqH2O,           544
.equiv          nb132nf_rsqH2H1,          560
.equiv          nb132nf_rsqH2H2,          576
.equiv          nb132nf_rinvOO,           592
.equiv          nb132nf_rinvOH1,          608
.equiv          nb132nf_rinvOH2,          624
.equiv          nb132nf_rinvH1O,          640
.equiv          nb132nf_rinvH1H1,         656
.equiv          nb132nf_rinvH1H2,         672
.equiv          nb132nf_rinvH2O,          688
.equiv          nb132nf_rinvH2H1,         704
.equiv          nb132nf_rinvH2H2,         720
.equiv          nb132nf_is3,              768
.equiv          nb132nf_ii3,              772
.equiv          nb132nf_innerjjnr,        776
.equiv          nb132nf_innerk,           780
.equiv          nb132nf_n,                784
.equiv          nb132nf_nn1,              788
.equiv          nb132nf_nri,              792
.equiv          nb132nf_nouter,           796
.equiv          nb132nf_ninner,           800
.equiv          nb132nf_salign,           804
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 808		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb132nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb132nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb132nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb132nf_nouter], eax
	mov [esp + nb132nf_ninner], eax

	mov eax, [ebp + nb132nf_p_tabscale]
	movss xmm3, [eax]
	shufps xmm3, xmm3, 0
	movaps [esp + nb132nf_tsc], xmm3

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb132nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb132nf_charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, xmm3	
	movss xmm5, [edx + ebx*4 + 4]	
	mov esi, [ebp + nb132nf_p_facel]
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
	movaps [esp + nb132nf_qqOO], xmm3
	movaps [esp + nb132nf_qqOH], xmm4
	movaps [esp + nb132nf_qqHH], xmm5
		
	xorps xmm0, xmm0
	mov   edx, [ebp + nb132nf_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov edi, [ebp + nb132nf_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   eax, [ebp + nb132nf_vdwparam]
	movlps xmm0, [eax + edx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 85  ;# constant 01010101
	movaps [esp + nb132nf_c6], xmm0
	movaps [esp + nb132nf_c12], xmm1

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb132nf_half], eax
	movss xmm1, [esp + nb132nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb132nf_half],  xmm1
	movaps [esp + nb132nf_three],  xmm3

.nb132nf_threadloop:
        mov   esi, [ebp + nb132nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb132nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb132nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb132nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb132nf_n], eax
        mov [esp + nb132nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb132nf_outerstart
        jmp .nb132nf_end

.nb132nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb132nf_nouter]
	mov [esp + nb132nf_nouter], ebx

.nb132nf_outer:
	mov   eax, [ebp + nb132nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb132nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb132nf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb132nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb132nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb132nf_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb132nf_ixO], xmm3
	movaps [esp + nb132nf_iyO], xmm4
	movaps [esp + nb132nf_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb132nf_ixH1], xmm0
	movaps [esp + nb132nf_iyH1], xmm1
	movaps [esp + nb132nf_izH1], xmm2
	movaps [esp + nb132nf_ixH2], xmm3
	movaps [esp + nb132nf_iyH2], xmm4
	movaps [esp + nb132nf_izH2], xmm5

	;# clear vctot
	xorps xmm4, xmm4
	movaps [esp + nb132nf_vctot], xmm4
	movaps [esp + nb132nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb132nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb132nf_pos]
	mov   eax, [ebp + nb132nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb132nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb132nf_ninner]
	mov   [esp + nb132nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb132nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb132nf_unroll_loop
	jmp   .nb132nf_single_check
.nb132nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb132nf_innerjjnr]     ;# pointer to jjnr[k] 

	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	mov   ecx, [edx + 8]
	mov   edx, [edx + 12]         ;# eax-edx=jnr1-4 
	
	add dword ptr [esp + nb132nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb132nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]     ;# replace jnr with j3 
	lea   edx, [edx + edx*2]	
	
	;# move j coordinates to local temp variables 
	movlps xmm2, [esi + eax*4]
	movlps xmm3, [esi + eax*4 + 12]
	movlps xmm4, [esi + eax*4 + 24]

	movlps xmm5, [esi + ebx*4]
	movlps xmm6, [esi + ebx*4 + 12]
	movlps xmm7, [esi + ebx*4 + 24]

	movhps xmm2, [esi + ecx*4]
	movhps xmm3, [esi + ecx*4 + 12]
	movhps xmm4, [esi + ecx*4 + 24]

	movhps xmm5, [esi + edx*4]
	movhps xmm6, [esi + edx*4 + 12]
	movhps xmm7, [esi + edx*4 + 24]

	;# current state: 	
	;# xmm2= jxOa  jyOa  jxOc  jyOc 
	;# xmm3= jxH1a jyH1a jxH1c jyH1c 
	;# xmm4= jxH2a jyH2a jxH2c jyH2c 
	;# xmm5= jxOb  jyOb  jxOd  jyOd 
	;# xmm6= jxH1b jyH1b jxH1d jyH1d 
	;# xmm7= jxH2b jyH2b jxH2d jyH2d 
	
	movaps xmm0, xmm2
	movaps xmm1, xmm3
	unpcklps xmm0, xmm5	;# xmm0= jxOa  jxOb  jyOa  jyOb 
	unpcklps xmm1, xmm6	;# xmm1= jxH1a jxH1b jyH1a jyH1b 
	unpckhps xmm2, xmm5	;# xmm2= jxOc  jxOd  jyOc  jyOd 
	unpckhps xmm3, xmm6	;# xmm3= jxH1c jxH1d jyH1c jyH1d 
	movaps xmm5, xmm4
	movaps   xmm6, xmm0
	unpcklps xmm4, xmm7	;# xmm4= jxH2a jxH2b jyH2a jyH2b 		
	unpckhps xmm5, xmm7	;# xmm5= jxH2c jxH2d jyH2c jyH2d 
	movaps   xmm7, xmm1
	movlhps  xmm0, xmm2	;# xmm0= jxOa  jxOb  jxOc  jxOd 
	movaps [esp + nb132nf_jxO], xmm0
	movhlps  xmm2, xmm6	;# xmm2= jyOa  jyOb  jyOc  jyOd 
	movaps [esp + nb132nf_jyO], xmm2
	movlhps  xmm1, xmm3
	movaps [esp + nb132nf_jxH1], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [esp + nb132nf_jyH1], xmm3
	movlhps  xmm4, xmm5
	movaps [esp + nb132nf_jxH2], xmm4
	movhlps  xmm5, xmm6
	movaps [esp + nb132nf_jyH2], xmm5

	movss  xmm0, [esi + eax*4 + 8]
	movss  xmm1, [esi + eax*4 + 20]
	movss  xmm2, [esi + eax*4 + 32]

	movss  xmm3, [esi + ecx*4 + 8]
	movss  xmm4, [esi + ecx*4 + 20]
	movss  xmm5, [esi + ecx*4 + 32]

	movhps xmm0, [esi + ebx*4 + 4]
	movhps xmm1, [esi + ebx*4 + 16]
	movhps xmm2, [esi + ebx*4 + 28]
	
	movhps xmm3, [esi + edx*4 + 4]
	movhps xmm4, [esi + edx*4 + 16]
	movhps xmm5, [esi + edx*4 + 28]
	
	shufps xmm0, xmm3, 204  ;# constant 11001100
	shufps xmm1, xmm4, 204  ;# constant 11001100
	shufps xmm2, xmm5, 204  ;# constant 11001100
	movaps [esp + nb132nf_jzO],  xmm0
	movaps [esp + nb132nf_jzH1],  xmm1
	movaps [esp + nb132nf_jzH2],  xmm2

	movaps xmm0, [esp + nb132nf_ixO]
	movaps xmm1, [esp + nb132nf_iyO]
	movaps xmm2, [esp + nb132nf_izO]
	movaps xmm3, [esp + nb132nf_ixO]
	movaps xmm4, [esp + nb132nf_iyO]
	movaps xmm5, [esp + nb132nf_izO]
	subps  xmm0, [esp + nb132nf_jxO]
	subps  xmm1, [esp + nb132nf_jyO]
	subps  xmm2, [esp + nb132nf_jzO]
	subps  xmm3, [esp + nb132nf_jxH1]
	subps  xmm4, [esp + nb132nf_jyH1]
	subps  xmm5, [esp + nb132nf_jzH1]
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
	movaps [esp + nb132nf_rsqOO], xmm0
	movaps [esp + nb132nf_rsqOH1], xmm3

	movaps xmm0, [esp + nb132nf_ixO]
	movaps xmm1, [esp + nb132nf_iyO]
	movaps xmm2, [esp + nb132nf_izO]
	movaps xmm3, [esp + nb132nf_ixH1]
	movaps xmm4, [esp + nb132nf_iyH1]
	movaps xmm5, [esp + nb132nf_izH1]
	subps  xmm0, [esp + nb132nf_jxH2]
	subps  xmm1, [esp + nb132nf_jyH2]
	subps  xmm2, [esp + nb132nf_jzH2]
	subps  xmm3, [esp + nb132nf_jxO]
	subps  xmm4, [esp + nb132nf_jyO]
	subps  xmm5, [esp + nb132nf_jzO]
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
	movaps [esp + nb132nf_rsqOH2], xmm0
	movaps [esp + nb132nf_rsqH1O], xmm3

	movaps xmm0, [esp + nb132nf_ixH1]
	movaps xmm1, [esp + nb132nf_iyH1]
	movaps xmm2, [esp + nb132nf_izH1]
	movaps xmm3, [esp + nb132nf_ixH1]
	movaps xmm4, [esp + nb132nf_iyH1]
	movaps xmm5, [esp + nb132nf_izH1]
	subps  xmm0, [esp + nb132nf_jxH1]
	subps  xmm1, [esp + nb132nf_jyH1]
	subps  xmm2, [esp + nb132nf_jzH1]
	subps  xmm3, [esp + nb132nf_jxH2]
	subps  xmm4, [esp + nb132nf_jyH2]
	subps  xmm5, [esp + nb132nf_jzH2]
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
	movaps [esp + nb132nf_rsqH1H1], xmm0
	movaps [esp + nb132nf_rsqH1H2], xmm3

	movaps xmm0, [esp + nb132nf_ixH2]
	movaps xmm1, [esp + nb132nf_iyH2]
	movaps xmm2, [esp + nb132nf_izH2]
	movaps xmm3, [esp + nb132nf_ixH2]
	movaps xmm4, [esp + nb132nf_iyH2]
	movaps xmm5, [esp + nb132nf_izH2]
	subps  xmm0, [esp + nb132nf_jxO]
	subps  xmm1, [esp + nb132nf_jyO]
	subps  xmm2, [esp + nb132nf_jzO]
	subps  xmm3, [esp + nb132nf_jxH1]
	subps  xmm4, [esp + nb132nf_jyH1]
	subps  xmm5, [esp + nb132nf_jzH1]
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
	movaps [esp + nb132nf_rsqH2O], xmm0
	movaps [esp + nb132nf_rsqH2H1], xmm4

	movaps xmm0, [esp + nb132nf_ixH2]
	movaps xmm1, [esp + nb132nf_iyH2]
	movaps xmm2, [esp + nb132nf_izH2]
	subps  xmm0, [esp + nb132nf_jxH2]
	subps  xmm1, [esp + nb132nf_jyH2]
	subps  xmm2, [esp + nb132nf_jzH2]
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [esp + nb132nf_rsqH2H2], xmm0
		
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb132nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb132nf_half] ;# rinvH2H2 
	mulps   xmm7, [esp + nb132nf_half] ;# rinvH2H1 
	movaps  [esp + nb132nf_rinvH2H2], xmm3
	movaps  [esp + nb132nf_rinvH2H1], xmm7
	
	rsqrtps xmm1, [esp + nb132nf_rsqOO]
	rsqrtps xmm5, [esp + nb132nf_rsqOH1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb132nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb132nf_rsqOO]
	mulps   xmm5, [esp + nb132nf_rsqOH1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb132nf_half] 
	mulps   xmm7, [esp + nb132nf_half]
	movaps  [esp + nb132nf_rinvOO], xmm3
	movaps  [esp + nb132nf_rinvOH1], xmm7
	
	rsqrtps xmm1, [esp + nb132nf_rsqOH2]
	rsqrtps xmm5, [esp + nb132nf_rsqH1O]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb132nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb132nf_rsqOH2]
	mulps   xmm5, [esp + nb132nf_rsqH1O]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb132nf_half] 
	mulps   xmm7, [esp + nb132nf_half]
	movaps  [esp + nb132nf_rinvOH2], xmm3
	movaps  [esp + nb132nf_rinvH1O], xmm7
	
	rsqrtps xmm1, [esp + nb132nf_rsqH1H1]
	rsqrtps xmm5, [esp + nb132nf_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb132nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb132nf_rsqH1H1]
	mulps   xmm5, [esp + nb132nf_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb132nf_half] 
	mulps   xmm7, [esp + nb132nf_half]
	movaps  [esp + nb132nf_rinvH1H1], xmm3
	movaps  [esp + nb132nf_rinvH1H2], xmm7
	
	rsqrtps xmm1, [esp + nb132nf_rsqH2O]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + nb132nf_three]
	mulps   xmm1, [esp + nb132nf_rsqH2O]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [esp + nb132nf_half] 
	movaps  [esp + nb132nf_rinvH2O], xmm3

	;# start with OO interaction - first the table LJ part
	movaps xmm0, [esp + nb132nf_rinvOO]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + nb132nf_rsqOO] ;# xmm1=r 
	mulps  xmm1, [esp + nb132nf_tsc]
		
	movhlps xmm2, xmm1
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 3
    pslld mm7, 3

    movd eax, mm6
    psrlq mm6, 32
    movd ecx, mm7
    psrlq mm7, 32
    movd ebx, mm6
    movd edx, mm7

    mov  esi, [ebp + nb132nf_VFtab]
	
    ;# dispersion 
    movlps xmm5, [esi + eax*4]
    movlps xmm7, [esi + ecx*4]
    movhps xmm5, [esi + ebx*4]
    movhps xmm7, [esi + edx*4] ;# got half table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# constant 10001000
    shufps xmm5, xmm7, 221  ;# constant 11011101

    movlps xmm7, [esi + eax*4 + 8]
    movlps xmm3, [esi + ecx*4 + 8]
    movhps xmm7, [esi + ebx*4 + 8]
    movhps xmm3, [esi + edx*4 + 8] ;# other half of table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# constant 10001000
    shufps xmm7, xmm3, 221  ;# constant 11011101
    ;# dispersion table ready, in xmm4-xmm7 
    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 

    movaps xmm4, [esp + nb132nf_c6]
    mulps  xmm5, xmm4    ;# Vvdw6 

    addps  xmm5, [esp + nb132nf_Vvdwtot]
    movaps [esp + nb132nf_Vvdwtot], xmm5

    ;# repulsion 
    movlps xmm5, [esi + eax*4 + 16]
    movlps xmm7, [esi + ecx*4 + 16]
    movhps xmm5, [esi + ebx*4 + 16]
    movhps xmm7, [esi + edx*4 + 16] ;# got half table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# constant 10001000
    shufps xmm5, xmm7, 221  ;# constant 11011101

    movlps xmm7, [esi + eax*4 + 24]
    movlps xmm3, [esi + ecx*4 + 24]
    movhps xmm7, [esi + ebx*4 + 24]
    movhps xmm3, [esi + edx*4 + 24] ;# other half of table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# constant 10001000
    shufps xmm7, xmm3, 221  ;# constant 11011101
    ;# repulsion table ready, in xmm4-xmm7 
    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
 
    movaps xmm4, [esp + nb132nf_c12]
    mulps  xmm5, xmm4 ;# Vvdw12 

    addps  xmm5, [esp + nb132nf_Vvdwtot]
    movaps [esp + nb132nf_Vvdwtot], xmm5

	;# Coulomb interactions 
	movaps xmm0, [esp + nb132nf_rinvOH1]
	movaps xmm1, [esp + nb132nf_rinvH1H1]
	movaps xmm2, [esp + nb132nf_rinvOO]
	addps  xmm0, [esp + nb132nf_rinvOH2]
	addps  xmm1, [esp + nb132nf_rinvH1H2]
	addps  xmm0, [esp + nb132nf_rinvH1O]
	addps  xmm1, [esp + nb132nf_rinvH2H1]
	addps  xmm0, [esp + nb132nf_rinvH2O]
	addps  xmm1, [esp + nb132nf_rinvH2H2]
	
	mulps  xmm0, [esp + nb132nf_qqOH]
	mulps  xmm1, [esp + nb132nf_qqHH]
	mulps  xmm2, [esp + nb132nf_qqOO]
	
	addps  xmm0, xmm1
	addps  xmm2, [esp + nb132nf_vctot]
	addps  xmm0, xmm2
	movaps [esp + nb132nf_vctot], xmm0

	;# should we do one more iteration? 
	sub dword ptr [esp + nb132nf_innerk],  4
	jl    .nb132nf_single_check
	jmp   .nb132nf_unroll_loop
.nb132nf_single_check:
	add dword ptr [esp + nb132nf_innerk],  4
	jnz   .nb132nf_single_loop
	jmp   .nb132nf_updateouterdata
.nb132nf_single_loop:
	mov   edx, [esp + nb132nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb132nf_innerjjnr],  4	

	mov esi, [ebp + nb132nf_pos]
	lea   eax, [eax + eax*2]  

	;# fetch j coordinates 
	xorps xmm3, xmm3
	xorps xmm4, xmm4
	xorps xmm5, xmm5
	
	movss xmm3, [esi + eax*4]		;# jxO  -  -  -
	movss xmm4, [esi + eax*4 + 4]		;# jyO  -  -  -
	movss xmm5, [esi + eax*4 + 8]		;# jzO  -  -  -  

	movlps xmm6, [esi + eax*4 + 12]		;# xmm6 = jxH1 jyH1   -    -
	movss  xmm7, [esi + eax*4 + 20]		;# xmm7 = jzH1   -    -    - 
	movhps xmm6, [esi + eax*4 + 24]		;# xmm6 = jxH1 jyH1 jxH2 jyH2
	movss  xmm2, [esi + eax*4 + 32]		;# xmm2 = jzH2   -    -    -
	
	;# have all coords, time for some shuffling.

	shufps xmm6, xmm6, 216 ;# constant 11011000	;# xmm6 = jxH1 jxH2 jyH1 jyH2 
	unpcklps xmm7, xmm2			;# xmm7 = jzH1 jzH2   -    -
	movaps  xmm0, [esp + nb132nf_ixO]     
	movaps  xmm1, [esp + nb132nf_iyO]
	movaps  xmm2, [esp + nb132nf_izO]	
	movlhps xmm3, xmm6			;# xmm3 = jxO   0   jxH1 jxH2 
	shufps  xmm4, xmm6, 228 ;# constant 11100100	;# xmm4 = jyO   0   jyH1 jyH2 
	shufps  xmm5, xmm7, 68  ;# constant 01000100	;# xmm5 = jzO   0   jzH1 jzH2
	
	;# store all j coordinates in jO  
	movaps [esp + nb132nf_jxO], xmm3
	movaps [esp + nb132nf_jyO], xmm4
	movaps [esp + nb132nf_jzO], xmm5
	subps  xmm0, xmm3
	subps  xmm1, xmm4
	subps  xmm2, xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 
	movaps [esp + nb132nf_rsqOO], xmm0
	
	movaps xmm6, xmm0
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + nb132nf_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [esp + nb132nf_half] ;# rinv iO - j water in xmm3
	movaps  [esp + nb132nf_rinvOO], xmm3
	
	
	xorps   xmm4, xmm4
	;# fetch charges to xmm4 (temporary) 
	movss   xmm4, [esp + nb132nf_qqOO]
	movhps  xmm4, [esp + nb132nf_qqOH]

	mulps xmm3, xmm4	;# vcoul  
	addps  xmm3, [esp + nb132nf_vctot]
    movaps [esp + nb132nf_vctot], xmm3


	movaps xmm0, [esp + nb132nf_rinvOO]
	movss xmm1, xmm0
	mulss  xmm1, [esp + nb132nf_rsqOO] ;# xmm1=r 
	mulss  xmm1, [esp + nb132nf_tsc]
		
    cvttps2pi mm6, xmm1
    cvtpi2ps xmm3, mm6
	subss    xmm1, xmm3	;# xmm1=eps 
    movss xmm2, xmm1
    mulss  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 3
	
    mov  esi, [ebp + nb132nf_VFtab]
    movd eax, mm6
	
    ;# dispersion 
    movlps xmm5, [esi + eax*4]
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# constant 10001000
    shufps xmm5, xmm7, 221  ;# constant 11011101

    movlps xmm7, [esi + eax*4 + 8]
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# constant 10001000
    shufps xmm7, xmm3, 221  ;# constant 11011101
    ;# dispersion table ready, in xmm4-xmm7 
    mulss  xmm6, xmm1       ;# xmm6=Geps 
    mulss  xmm7, xmm2       ;# xmm7=Heps2 
    addss  xmm5, xmm6
    addss  xmm5, xmm7       ;# xmm5=Fp 
    mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
    addss  xmm5, xmm4 ;# xmm5=VV 

    movss xmm4, [esp + nb132nf_c6]
    mulss  xmm5, xmm4    ;# Vvdw6 

    addss  xmm5, [esp + nb132nf_Vvdwtot]
    movss [esp + nb132nf_Vvdwtot], xmm5

    ;# repulsion 
    movlps xmm5, [esi + eax*4 + 16]
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# constant 10001000
    shufps xmm5, xmm7, 221  ;# constant 11011101

    movlps xmm7, [esi + eax*4 + 24]
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# constant 10001000
    shufps xmm7, xmm3, 221  ;# constant 11011101
    ;# table ready, in xmm4-xmm7 
    mulss  xmm6, xmm1       ;# xmm6=Geps 
    mulss  xmm7, xmm2       ;# xmm7=Heps2 
    addss  xmm5, xmm6
    addss  xmm5, xmm7       ;# xmm5=Fp 
    mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
    addss  xmm5, xmm4 ;# xmm5=VV 

    movss xmm4, [esp + nb132nf_c12]
    mulss  xmm5, xmm4 ;# Vvdw12 
    addss  xmm5, [esp + nb132nf_Vvdwtot]
    movss [esp + nb132nf_Vvdwtot], xmm5

	;# done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [esp + nb132nf_ixH1]
	movaps  xmm1, [esp + nb132nf_iyH1]
	movaps  xmm2, [esp + nb132nf_izH1]	
	movaps  xmm3, [esp + nb132nf_ixH2] 
	movaps  xmm4, [esp + nb132nf_iyH2] 
	movaps  xmm5, [esp + nb132nf_izH2] 
	subps   xmm0, [esp + nb132nf_jxO]
	subps   xmm1, [esp + nb132nf_jyO]
	subps   xmm2, [esp + nb132nf_jzO]
	subps   xmm3, [esp + nb132nf_jxO]
	subps   xmm4, [esp + nb132nf_jyO]
	subps   xmm5, [esp + nb132nf_jzO]
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
	movaps  xmm3, [esp + nb132nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb132nf_half] ;# rinv H1 - j water 
	mulps   xmm7, [esp + nb132nf_half] ;# rinv H2 - j water  

	;# assemble charges in xmm6 
	xorps   xmm6, xmm6
	movss   xmm6, [esp + nb132nf_qqOH]
	movhps  xmm6, [esp + nb132nf_qqHH]

	addps   xmm3, xmm7
	mulps   xmm3, xmm6

	addps   xmm3, [esp + nb132nf_vctot]
	movaps  [esp + nb132nf_vctot], xmm3
	
	dec dword ptr [esp + nb132nf_innerk]
	jz    .nb132nf_updateouterdata
	jmp   .nb132nf_single_loop
.nb132nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb132nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb132nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb132nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb132nf_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [esp + nb132nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb132nf_Vvdw]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb132nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb132nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb132nf_n], esi
        jmp .nb132nf_outer
.nb132nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb132nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb132nf_end
        ;# non-zero, do one more workunit
        jmp   .nb132nf_threadloop
.nb132nf_end:
	emms

	mov eax, [esp + nb132nf_nouter]
	mov ebx, [esp + nb132nf_ninner]
	mov ecx, [ebp + nb132nf_outeriter]
	mov edx, [ebp + nb132nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb132nf_salign]
	add esp, eax
	add esp, 808
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret

