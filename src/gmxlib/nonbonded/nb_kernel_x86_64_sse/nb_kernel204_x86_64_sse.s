;#
;# $Id$
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
.equiv          .equiv                  2
   %1 equ %2
%endmacro
; .endif                   # End of NASM-specific block
; .intel_syntax noprefix   # Line only read by gnu as



	
.globl nb_kernel204_x86_64_sse
.globl _nb_kernel204_x86_64_sse
nb_kernel204_x86_64_sse:	
_nb_kernel204_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb204_fshift,           16
.equiv          nb204_gid,              24
.equiv          nb204_pos,              32
.equiv          nb204_faction,          40
.equiv          nb204_charge,           48
.equiv          nb204_p_facel,          56
.equiv          nb204_argkrf,           64
.equiv          nb204_argcrf,           72
.equiv          nb204_Vc,               80
.equiv          nb204_type,             88
.equiv          nb204_p_ntype,          96
.equiv          nb204_vdwparam,         104
.equiv          nb204_Vvdw,             112
.equiv          nb204_p_tabscale,       120
.equiv          nb204_VFtab,            128
.equiv          nb204_invsqrta,         136
.equiv          nb204_dvda,             144
.equiv          nb204_p_gbtabscale,     152
.equiv          nb204_GBtab,            160
.equiv          nb204_p_nthreads,       168
.equiv          nb204_count,            176
.equiv          nb204_mtx,              184
.equiv          nb204_outeriter,        192
.equiv          nb204_inneriter,        200
.equiv          nb204_work,             208
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
.equiv          nb204_nri,              1456
.equiv          nb204_iinr,             1464
.equiv          nb204_jindex,           1472
.equiv          nb204_jjnr,             1480
.equiv          nb204_shift,            1488
.equiv          nb204_shiftvec,         1496
.equiv          nb204_facel,            1504
.equiv          nb204_innerk,           1512
.equiv          nb204_n,                1516
.equiv          nb204_nn1,              1520
.equiv          nb204_nouter,           1524
.equiv          nb204_ninner,           1528

	push rbp
	mov  rbp, rsp
	push rbx

	femms
	sub rsp, 1544		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb204_nouter], eax
	mov [rsp + nb204_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb204_nri], edi
	mov [rsp + nb204_iinr], rsi
	mov [rsp + nb204_jindex], rdx
	mov [rsp + nb204_jjnr], rcx
	mov [rsp + nb204_shift], r8
	mov [rsp + nb204_shiftvec], r9
	mov rsi, [rbp + nb204_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb204_facel], xmm0


	mov rsi, [rbp + nb204_argkrf]
	mov rdi, [rbp + nb204_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb204_krf], xmm1
	movaps [rsp + nb204_crf], xmm2
	
	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb204_half], eax
	movss xmm1, [rsp + nb204_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb204_half],  xmm1
	movaps [rsp + nb204_two],  xmm2
	movaps [rsp + nb204_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb204_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]		;# ebx =ii 

	mov   rdx, [rbp + nb204_charge]
	movss xmm3, [rdx + rbx*4 + 4]	
	movss xmm4, xmm3	
	movss xmm5, [rdx + rbx*4 + 12]	
	mov rsi, [rbp + nb204_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb204_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb204_qqHH], xmm3
	movaps [rsp + nb204_qqMH], xmm4
	movaps [rsp + nb204_qqMM], xmm5
	
.nb204_threadloop:
        mov   rsi, [rbp + nb204_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb204_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   rbx, 1                           ;# rbx=nn1=nn0+10
        lock cmpxchg [rsi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb204_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb204_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb204_n], eax
        mov [rsp + nb204_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb204_outerstart
        jmp .nb204_end
	
.nb204_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb204_nouter]
	mov [rsp + nb204_nouter], ebx

.nb204_outer:
	mov   rax, [rsp + nb204_shift]  	;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]	;# ebx=3*is 
	mov   [rsp + nb204_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb204_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb204_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   rax, [rbp + nb204_pos]	;# rax = base of pos[]  
	mov   [rsp + nb204_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [rax + rbx*4 + 12]
	addss xmm4, [rax + rbx*4 + 16]
	addss xmm5, [rax + rbx*4 + 20]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb204_ixH1], xmm3
	movaps [rsp + nb204_iyH1], xmm4
	movaps [rsp + nb204_izH1], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [rax + rbx*4 + 24]
	addss xmm1, [rax + rbx*4 + 28]
	addss xmm2, [rax + rbx*4 + 32]		
	addss xmm3, [rax + rbx*4 + 36]
	addss xmm4, [rax + rbx*4 + 40]
	addss xmm5, [rax + rbx*4 + 44]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb204_ixH2], xmm0
	movaps [rsp + nb204_iyH2], xmm1
	movaps [rsp + nb204_izH2], xmm2
	movaps [rsp + nb204_ixM], xmm3
	movaps [rsp + nb204_iyM], xmm4
	movaps [rsp + nb204_izM], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb204_vctot], xmm4
	movaps [rsp + nb204_fixH1], xmm4
	movaps [rsp + nb204_fiyH1], xmm4
	movaps [rsp + nb204_fizH1], xmm4
	movaps [rsp + nb204_fixH2], xmm4
	movaps [rsp + nb204_fiyH2], xmm4
	movaps [rsp + nb204_fizH2], xmm4
	movaps [rsp + nb204_fixM], xmm4
	movaps [rsp + nb204_fiyM], xmm4
	movaps [rsp + nb204_fizM], xmm4
	
	mov   rax, [rsp + nb204_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb204_pos]
	mov   rdi, [rbp + nb204_faction]	
	mov   rax, [rsp + nb204_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb204_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb204_ninner]
	mov   [rsp + nb204_ninner], ecx
	add   edx, 0
	mov   [rsp + nb204_innerk], edx	;# number of innerloop atoms 
	jge   .nb204_unroll_loop
	jmp   .nb204_single_check
.nb204_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb204_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb204_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb204_pos]   	;# base of pos[] 

	lea   eax, [eax + eax*2] 	;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2] 	;# replace jnr with j3 
	lea   edx, [edx + edx*2]	
	
	;# move j coordinates to local temp variables 
	movlps xmm2, [rsi + rax*4 + 12]
	movlps xmm3, [rsi + rax*4 + 24]
	movlps xmm4, [rsi + rax*4 + 36]

	movlps xmm5, [rsi + rbx*4 + 12]
	movlps xmm6, [rsi + rbx*4 + 24]
	movlps xmm7, [rsi + rbx*4 + 36]

	movhps xmm2, [rsi + rcx*4 + 12]
	movhps xmm3, [rsi + rcx*4 + 24]
	movhps xmm4, [rsi + rcx*4 + 36]

	movhps xmm5, [rsi + rdx*4 + 12]
	movhps xmm6, [rsi + rdx*4 + 24]
	movhps xmm7, [rsi + rdx*4 + 36]
	
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
	movaps [rsp + nb204_jxH1], xmm0
	movhlps  xmm2, xmm6
	movaps [rsp + nb204_jyH1], xmm2
	movlhps  xmm1, xmm3
	movaps [rsp + nb204_jxH2], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [rsp + nb204_jyH2], xmm3
	movlhps  xmm4, xmm5
	movaps [rsp + nb204_jxM], xmm4
	movhlps  xmm5, xmm6
	movaps [rsp + nb204_jyM], xmm5

	movss  xmm0, [rsi + rax*4 + 20]
	movss  xmm1, [rsi + rax*4 + 32]
	movss  xmm2, [rsi + rax*4 + 44]

	movss  xmm3, [rsi + rcx*4 + 20]
	movss  xmm4, [rsi + rcx*4 + 32]
	movss  xmm5, [rsi + rcx*4 + 44]

	movhps xmm0, [rsi + rbx*4 + 16]
	movhps xmm1, [rsi + rbx*4 + 28]
	movhps xmm2, [rsi + rbx*4 + 40]
	
	movhps xmm3, [rsi + rdx*4 + 16]
	movhps xmm4, [rsi + rdx*4 + 28]
	movhps xmm5, [rsi + rdx*4 + 40]
	
	shufps xmm0, xmm3, 204  ;# 11001100
	shufps xmm1, xmm4, 204  ;# 11001100
	shufps xmm2, xmm5, 204  ;# 11001100
	movaps [rsp + nb204_jzH1],  xmm0
	movaps [rsp + nb204_jzH2],  xmm1
	movaps [rsp + nb204_jzM],  xmm2

	movaps xmm0, [rsp + nb204_ixH1]
	movaps xmm1, [rsp + nb204_iyH1]
	movaps xmm2, [rsp + nb204_izH1]
	movaps xmm3, [rsp + nb204_ixH1]
	movaps xmm4, [rsp + nb204_iyH1]
	movaps xmm5, [rsp + nb204_izH1]
	subps  xmm0, [rsp + nb204_jxH1]
	subps  xmm1, [rsp + nb204_jyH1]
	subps  xmm2, [rsp + nb204_jzH1]
	subps  xmm3, [rsp + nb204_jxH2]
	subps  xmm4, [rsp + nb204_jyH2]
	subps  xmm5, [rsp + nb204_jzH2]
	movaps [rsp + nb204_dxH1H1], xmm0
	movaps [rsp + nb204_dyH1H1], xmm1
	movaps [rsp + nb204_dzH1H1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb204_dxH1H2], xmm3
	movaps [rsp + nb204_dyH1H2], xmm4
	movaps [rsp + nb204_dzH1H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + nb204_rsqH1H1], xmm0
	movaps [rsp + nb204_rsqH1H2], xmm3

	movaps xmm0, [rsp + nb204_ixH1]
	movaps xmm1, [rsp + nb204_iyH1]
	movaps xmm2, [rsp + nb204_izH1]
	movaps xmm3, [rsp + nb204_ixH2]
	movaps xmm4, [rsp + nb204_iyH2]
	movaps xmm5, [rsp + nb204_izH2]
	subps  xmm0, [rsp + nb204_jxM]
	subps  xmm1, [rsp + nb204_jyM]
	subps  xmm2, [rsp + nb204_jzM]
	subps  xmm3, [rsp + nb204_jxH1]
	subps  xmm4, [rsp + nb204_jyH1]
	subps  xmm5, [rsp + nb204_jzH1]
	movaps [rsp + nb204_dxH1M], xmm0
	movaps [rsp + nb204_dyH1M], xmm1
	movaps [rsp + nb204_dzH1M], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb204_dxH2H1], xmm3
	movaps [rsp + nb204_dyH2H1], xmm4
	movaps [rsp + nb204_dzH2H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + nb204_rsqH1M], xmm0
	movaps [rsp + nb204_rsqH2H1], xmm3

	movaps xmm0, [rsp + nb204_ixH2]
	movaps xmm1, [rsp + nb204_iyH2]
	movaps xmm2, [rsp + nb204_izH2]
	movaps xmm3, [rsp + nb204_ixH2]
	movaps xmm4, [rsp + nb204_iyH2]
	movaps xmm5, [rsp + nb204_izH2]
	subps  xmm0, [rsp + nb204_jxH2]
	subps  xmm1, [rsp + nb204_jyH2]
	subps  xmm2, [rsp + nb204_jzH2]
	subps  xmm3, [rsp + nb204_jxM]
	subps  xmm4, [rsp + nb204_jyM]
	subps  xmm5, [rsp + nb204_jzM]
	movaps [rsp + nb204_dxH2H2], xmm0
	movaps [rsp + nb204_dyH2H2], xmm1
	movaps [rsp + nb204_dzH2H2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb204_dxH2M], xmm3
	movaps [rsp + nb204_dyH2M], xmm4
	movaps [rsp + nb204_dzH2M], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + nb204_rsqH2H2], xmm0
	movaps [rsp + nb204_rsqH2M], xmm3

	movaps xmm0, [rsp + nb204_ixM]
	movaps xmm1, [rsp + nb204_iyM]
	movaps xmm2, [rsp + nb204_izM]
	movaps xmm3, [rsp + nb204_ixM]
	movaps xmm4, [rsp + nb204_iyM]
	movaps xmm5, [rsp + nb204_izM]
	subps  xmm0, [rsp + nb204_jxH1]
	subps  xmm1, [rsp + nb204_jyH1]
	subps  xmm2, [rsp + nb204_jzH1]
	subps  xmm3, [rsp + nb204_jxH2]
	subps  xmm4, [rsp + nb204_jyH2]
	subps  xmm5, [rsp + nb204_jzH2]
	movaps [rsp + nb204_dxMH1], xmm0
	movaps [rsp + nb204_dyMH1], xmm1
	movaps [rsp + nb204_dzMH1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb204_dxMH2], xmm3
	movaps [rsp + nb204_dyMH2], xmm4
	movaps [rsp + nb204_dzMH2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [rsp + nb204_rsqMH1], xmm0
	movaps [rsp + nb204_rsqMH2], xmm4

	movaps xmm0, [rsp + nb204_ixM]
	movaps xmm1, [rsp + nb204_iyM]
	movaps xmm2, [rsp + nb204_izM]
	subps  xmm0, [rsp + nb204_jxM]
	subps  xmm1, [rsp + nb204_jyM]
	subps  xmm2, [rsp + nb204_jzM]
	movaps [rsp + nb204_dxMM], xmm0
	movaps [rsp + nb204_dyMM], xmm1
	movaps [rsp + nb204_dzMM], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [rsp + nb204_rsqMM], xmm0
	
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb204_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb204_half] ;# rinvH2H2 
	mulps   xmm7, [rsp + nb204_half] ;# rinvH2H1 
	movaps  [rsp + nb204_rinvMM], xmm3
	movaps  [rsp + nb204_rinvMH2], xmm7
	
	rsqrtps xmm1, [rsp + nb204_rsqH1H1]
	rsqrtps xmm5, [rsp + nb204_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb204_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb204_rsqH1H1]
	mulps   xmm5, [rsp + nb204_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb204_half] 
	mulps   xmm7, [rsp + nb204_half]
	movaps  [rsp + nb204_rinvH1H1], xmm3
	movaps  [rsp + nb204_rinvH1H2], xmm7
	
	rsqrtps xmm1, [rsp + nb204_rsqH1M]
	rsqrtps xmm5, [rsp + nb204_rsqH2H1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb204_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb204_rsqH1M]
	mulps   xmm5, [rsp + nb204_rsqH2H1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb204_half] 
	mulps   xmm7, [rsp + nb204_half]
	movaps  [rsp + nb204_rinvH1M], xmm3
	movaps  [rsp + nb204_rinvH2H1], xmm7
	
	rsqrtps xmm1, [rsp + nb204_rsqH2H2]
	rsqrtps xmm5, [rsp + nb204_rsqH2M]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb204_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb204_rsqH2H2]
	mulps   xmm5, [rsp + nb204_rsqH2M]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb204_half] 
	mulps   xmm7, [rsp + nb204_half]
	movaps  [rsp + nb204_rinvH2H2], xmm3
	movaps  [rsp + nb204_rinvH2M], xmm7
	
	rsqrtps xmm1, [rsp + nb204_rsqMH1]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb204_three]
	mulps   xmm1, [rsp + nb204_rsqMH1]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + nb204_half] 
	movaps  [rsp + nb204_rinvMH1], xmm3

	;# start with H1-H1 interaction 
	movaps xmm0, [rsp + nb204_rinvH1H1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb204_krf]
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	mulps  xmm5, [rsp + nb204_rsqH1H1] ;# xmm5=krsq 
	movaps xmm6, xmm5
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [rsp + nb204_crf]
	mulps  xmm6, [rsp + nb204_qqHH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	mulps xmm5, [rsp + nb204_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [rsp + nb204_qqHH] ;# xmm7 = coul part of fscal 
	
	addps  xmm6, [rsp + nb204_vctot] ;# local vctot summation variable 
	mulps  xmm0, xmm7
	
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [rsp + nb204_dxH1H1]
	mulps xmm1, [rsp + nb204_dyH1H1]
	mulps xmm2, [rsp + nb204_dzH1H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb204_fixH1]
	addps xmm1, [rsp + nb204_fiyH1]
	addps xmm2, [rsp + nb204_fizH1]
	movaps [rsp + nb204_fjxH1], xmm3
	movaps [rsp + nb204_fjyH1], xmm4
	movaps [rsp + nb204_fjzH1], xmm5
	movaps [rsp + nb204_fixH1], xmm0
	movaps [rsp + nb204_fiyH1], xmm1
	movaps [rsp + nb204_fizH1], xmm2

	;# H1-H2 interaction 
	movaps xmm0, [rsp + nb204_rinvH1H2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb204_krf]
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb204_rsqH1H2] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb204_crf]
	mulps  xmm0, xmm0
	mulps  xmm4, [rsp + nb204_qqHH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [rsp + nb204_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [rsp + nb204_qqHH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH1  
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [rsp + nb204_dxH1H2]
	mulps xmm1, [rsp + nb204_dyH1H2]
	mulps xmm2, [rsp + nb204_dzH1H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb204_fixH1]
	addps xmm1, [rsp + nb204_fiyH1]
	addps xmm2, [rsp + nb204_fizH1]
	movaps [rsp + nb204_fjxH2], xmm3
	movaps [rsp + nb204_fjyH2], xmm4
	movaps [rsp + nb204_fjzH2], xmm5
	movaps [rsp + nb204_fixH1], xmm0
	movaps [rsp + nb204_fiyH1], xmm1
	movaps [rsp + nb204_fizH1], xmm2

	;# H1-M interaction  
	movaps xmm0, [rsp + nb204_rinvH1M]
	movaps xmm7, xmm0	;# xmm7=Rinv 
	movaps xmm5, [rsp + nb204_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb204_rsqH1M] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb204_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb204_qqMH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [rsp + nb204_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [rsp + nb204_qqMH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [rsp + nb204_dxH1M]
	mulps xmm1, [rsp + nb204_dyH1M]
	mulps xmm2, [rsp + nb204_dzH1M]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb204_fixH1]
	addps xmm1, [rsp + nb204_fiyH1]
	addps xmm2, [rsp + nb204_fizH1]
	movaps [rsp + nb204_fjxM], xmm3
	movaps [rsp + nb204_fjyM], xmm4
	movaps [rsp + nb204_fjzM], xmm5
	movaps [rsp + nb204_fixH1], xmm0
	movaps [rsp + nb204_fiyH1], xmm1
	movaps [rsp + nb204_fizH1], xmm2

	;# H2-H1 interaction 
	movaps xmm0, [rsp + nb204_rinvH2H1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb204_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb204_rsqH2H1] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb204_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb204_qqHH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [rsp + nb204_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [rsp + nb204_qqHH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [rsp + nb204_fjxH1]
	movaps xmm4, [rsp + nb204_fjyH1]
	movaps xmm5, [rsp + nb204_fjzH1]
	mulps xmm0, [rsp + nb204_dxH2H1]
	mulps xmm1, [rsp + nb204_dyH2H1]
	mulps xmm2, [rsp + nb204_dzH2H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb204_fixH2]
	addps xmm1, [rsp + nb204_fiyH2]
	addps xmm2, [rsp + nb204_fizH2]
	movaps [rsp + nb204_fjxH1], xmm3
	movaps [rsp + nb204_fjyH1], xmm4
	movaps [rsp + nb204_fjzH1], xmm5
	movaps [rsp + nb204_fixH2], xmm0
	movaps [rsp + nb204_fiyH2], xmm1
	movaps [rsp + nb204_fizH2], xmm2

	;# H2-H2 interaction 
	movaps xmm0, [rsp + nb204_rinvH2H2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb204_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb204_rsqH2H2] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb204_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb204_qqHH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [rsp + nb204_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [rsp + nb204_qqHH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [rsp + nb204_fjxH2]
	movaps xmm4, [rsp + nb204_fjyH2]
	movaps xmm5, [rsp + nb204_fjzH2]
	mulps xmm0, [rsp + nb204_dxH2H2]
	mulps xmm1, [rsp + nb204_dyH2H2]
	mulps xmm2, [rsp + nb204_dzH2H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb204_fixH2]
	addps xmm1, [rsp + nb204_fiyH2]
	addps xmm2, [rsp + nb204_fizH2]
	movaps [rsp + nb204_fjxH2], xmm3
	movaps [rsp + nb204_fjyH2], xmm4
	movaps [rsp + nb204_fjzH2], xmm5
	movaps [rsp + nb204_fixH2], xmm0
	movaps [rsp + nb204_fiyH2], xmm1
	movaps [rsp + nb204_fizH2], xmm2

	;# H2-M interaction 
	movaps xmm0, [rsp + nb204_rinvH2M]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb204_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb204_rsqH2M] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb204_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb204_qqMH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [rsp + nb204_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [rsp + nb204_qqMH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	movaps xmm3, [rsp + nb204_fjxM]
	movaps xmm4, [rsp + nb204_fjyM]
	movaps xmm5, [rsp + nb204_fjzM]
	mulps xmm0, [rsp + nb204_dxH2M]
	mulps xmm1, [rsp + nb204_dyH2M]
	mulps xmm2, [rsp + nb204_dzH2M]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb204_fixH2]
	addps xmm1, [rsp + nb204_fiyH2]
	addps xmm2, [rsp + nb204_fizH2]
	movaps [rsp + nb204_fjxM], xmm3
	movaps [rsp + nb204_fjyM], xmm4
	movaps [rsp + nb204_fjzM], xmm5
	movaps [rsp + nb204_fixH2], xmm0
	movaps [rsp + nb204_fiyH2], xmm1
	movaps [rsp + nb204_fizH2], xmm2

	;# M-H1 interaction 
	movaps xmm0, [rsp + nb204_rinvMH1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb204_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb204_rsqMH1] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb204_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb204_qqMH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [rsp + nb204_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [rsp + nb204_qqMH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [rsp + nb204_fjxH1]
	movaps xmm4, [rsp + nb204_fjyH1]
	movaps xmm5, [rsp + nb204_fjzH1]
	mulps xmm0, [rsp + nb204_dxMH1]
	mulps xmm1, [rsp + nb204_dyMH1]
	mulps xmm2, [rsp + nb204_dzMH1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb204_fixM]
	addps xmm1, [rsp + nb204_fiyM]
	addps xmm2, [rsp + nb204_fizM]
	movaps [rsp + nb204_fjxH1], xmm3
	movaps [rsp + nb204_fjyH1], xmm4
	movaps [rsp + nb204_fjzH1], xmm5
	movaps [rsp + nb204_fixM], xmm0
	movaps [rsp + nb204_fiyM], xmm1
	movaps [rsp + nb204_fizM], xmm2

	;# M-H2 interaction 
	movaps xmm0, [rsp + nb204_rinvMH2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb204_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb204_rsqMH2] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb204_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb204_qqMH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [rsp + nb204_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [rsp + nb204_qqMH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [rsp + nb204_fjxH2]
	movaps xmm4, [rsp + nb204_fjyH2]
	movaps xmm5, [rsp + nb204_fjzH2]
	mulps xmm0, [rsp + nb204_dxMH2]
	mulps xmm1, [rsp + nb204_dyMH2]
	mulps xmm2, [rsp + nb204_dzMH2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb204_fixM]
	addps xmm1, [rsp + nb204_fiyM]
	addps xmm2, [rsp + nb204_fizM]
	movaps [rsp + nb204_fjxH2], xmm3
	movaps [rsp + nb204_fjyH2], xmm4
	movaps [rsp + nb204_fjzH2], xmm5
	movaps [rsp + nb204_fixM], xmm0
	movaps [rsp + nb204_fiyM], xmm1
	movaps [rsp + nb204_fizM], xmm2

	;# M-M interaction 
	movaps xmm0, [rsp + nb204_rinvMM]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb204_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb204_rsqMM] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb204_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb204_qqMM] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [rsp + nb204_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [rsp + nb204_qqMM] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm1, xmm0
	movaps [rsp + nb204_vctot], xmm6
	movaps xmm2, xmm0
	
	movaps xmm3, [rsp + nb204_fjxM]
	movaps xmm4, [rsp + nb204_fjyM]
	movaps xmm5, [rsp + nb204_fjzM]
	mulps xmm0, [rsp + nb204_dxMM]
	mulps xmm1, [rsp + nb204_dyMM]
	mulps xmm2, [rsp + nb204_dzMM]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb204_fixM]
	addps xmm1, [rsp + nb204_fiyM]
	addps xmm2, [rsp + nb204_fizM]
	movaps [rsp + nb204_fjxM], xmm3
	movaps [rsp + nb204_fjyM], xmm4
	movaps [rsp + nb204_fjzM], xmm5
	movaps [rsp + nb204_fixM], xmm0
	movaps [rsp + nb204_fiyM], xmm1
	movaps [rsp + nb204_fizM], xmm2

	mov rdi, [rbp + nb204_faction]
	
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
	movaps xmm0, [rsp + nb204_fjxH1] ;# xmm0= fjxH1a  fjxH1b  fjxH1c  fjxH1d 
	movaps xmm2, [rsp + nb204_fjyH1] ;# xmm1= fjyH1a  fjyH1b  fjyH1c  fjyH1d
	movlps xmm3, [rdi + rax*4 + 12]
	movlps xmm4, [rdi + rcx*4 + 12]
	movaps xmm1, xmm0
	unpcklps xmm0, xmm2    	   ;# xmm0= fjxH1a  fjyH1a  fjxH1b  fjyH1b
	unpckhps xmm1, xmm2        ;# xmm1= fjxH1c  fjyH1c  fjxH1d  fjyH1d
	movhps xmm3, [rdi + rbx*4 + 12]
	movhps xmm4, [rdi + rdx*4 + 12]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	movlps [rdi + rax*4 + 12], xmm3
	movlps [rdi + rcx*4 + 12], xmm4
	movhps [rdi + rbx*4 + 12], xmm3
	movhps [rdi + rdx*4 + 12], xmm4

	;# 1st Hydrogen Z & 2nd hydrogen X forces for 4 j particles 
	movaps xmm0, [rsp + nb204_fjzH1]  ;# xmm0= fjzH1a   fjzH1b   fjzH1c   fjzH1d 
	movaps xmm2, [rsp + nb204_fjxH2] ;# xmm1= fjxH2a  fjxH2b  fjxH2c  fjxH2d
	movlps xmm3, [rdi + rax*4 + 20]
	movlps xmm4, [rdi + rcx*4 + 20]
	movaps xmm1, xmm0
	unpcklps xmm0, xmm2    	   ;# xmm0= fjzH1a  fjxH2a  fjzH1b  fjxH2b
	unpckhps xmm1, xmm2        ;# xmm1= fjzH1c  fjxH2c  fjzH1d  fjxH2d
	movhps xmm3, [rdi + rbx*4 + 20]
	movhps xmm4, [rdi + rdx*4 + 20]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	movlps [rdi + rax*4 + 20], xmm3
	movlps [rdi + rcx*4 + 20], xmm4
	movhps [rdi + rbx*4 + 20], xmm3
	movhps [rdi + rdx*4 + 20], xmm4
	
	;# 2nd hydrogen Y & Z forces for 4 j particles 
	movaps xmm0, [rsp + nb204_fjyH2] ;# xmm0= fjyH2a  fjyH2b  fjyH2c  fjyH2d 
	movaps xmm2, [rsp + nb204_fjzH2] ;# xmm1= fjzH2a  fjzH2b  fjzH2c  fjzH2d
	movlps xmm3, [rdi + rax*4 + 28]
	movlps xmm4, [rdi + rcx*4 + 28]
	movaps xmm1, xmm0
	unpcklps xmm0, xmm2		;# xmm0= fjyH2a  fjzH2a  fjyH2b  fjzH2b
	unpckhps xmm1, xmm2		;# xmm1= fjyH2c  fjzH2c  fjyH2d  fjzH2d
	movhps xmm3, [rdi + rbx*4 + 28]
	movhps xmm4, [rdi + rdx*4 + 28]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	movlps [rdi + rax*4 + 28], xmm3
	movlps [rdi + rcx*4 + 28], xmm4
	movhps [rdi + rbx*4 + 28], xmm3
	movhps [rdi + rdx*4 + 28], xmm4

	;# Dummy (M) X & Y forces for 4 j particles 
	movaps xmm0, [rsp + nb204_fjxM] ;# xmm0= fjxMa  fjxMb  fjxMc  fjxMd 
	movaps xmm2, [rsp + nb204_fjyM] ;# xmm1= fjyMa  fjyMb  fjyMc  fjyMd
	movlps xmm3, [rdi + rax*4 + 36]
	movlps xmm4, [rdi + rcx*4 + 36]
	movaps xmm1, xmm0
	unpcklps xmm0, xmm2		;# xmm0= fjxMa  fjyMa  fjxMb  fjyMb
	unpckhps xmm1, xmm2		;# xmm1= fjxMc  fjyMc  fjxMd  fjyMd
	movhps xmm3, [rdi + rbx*4 + 36]
	movhps xmm4, [rdi + rdx*4 + 36]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	movlps [rdi + rax*4 + 36], xmm3
	movlps [rdi + rcx*4 + 36], xmm4
	movhps [rdi + rbx*4 + 36], xmm3
	movhps [rdi + rdx*4 + 36], xmm4

	
	;# Dummy (M) Z forces for 4 j particles 
	;# Just load the four Z coords into one reg. each
	movss xmm4, [rdi + rax*4 + 44]
	movss xmm5, [rdi + rbx*4 + 44]
	movss xmm6, [rdi + rcx*4 + 44]
	movss xmm7, [rdi + rdx*4 + 44]
	;# add what we have on the stack
	addss xmm4, [rsp + nb204_fjzM] 
	addss xmm5, [rsp + nb204_fjzM + 4] 
	addss xmm6, [rsp + nb204_fjzM + 8] 
	addss xmm7, [rsp + nb204_fjzM + 12]
	;# store back
	movss [rdi + rax*4 + 44], xmm4
	movss [rdi + rbx*4 + 44], xmm5
	movss [rdi + rcx*4 + 44], xmm6
	movss [rdi + rdx*4 + 44], xmm7
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb204_innerk],  4
	jl    .nb204_single_check
	jmp   .nb204_unroll_loop
.nb204_single_check:
	add dword ptr [rsp + nb204_innerk],  4
	jnz   .nb204_single_loop
	jmp   .nb204_updateouterdata
.nb204_single_loop:
	mov   rdx, [rsp + nb204_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb204_innerjjnr],  4	

	mov rsi, [rbp + nb204_pos]
	lea   eax, [eax + eax*2]  

	;# fetch j coordinates 
	xorps xmm3, xmm3
	xorps xmm4, xmm4
	xorps xmm5, xmm5
	movss xmm3, [rsi + rax*4 + 36]		;# jxM  -  -  -
	movss xmm4, [rsi + rax*4 + 40]		;# jyM  -  -  -
	movss xmm5, [rsi + rax*4 + 44]		;# jzM  -  -  -  

	movlps xmm6, [rsi + rax*4 + 12]		;# xmm6 = jxH1 jyH1   -    -
	movss  xmm7, [rsi + rax*4 + 20]		;# xmm7 = jzH1   -    -    - 
	movhps xmm6, [rsi + rax*4 + 24]		;# xmm6 = jxH1 jyH1 jxH2 jyH2
	movss  xmm2, [rsi + rax*4 + 32]		;# xmm2 = jzH2   -    -    -
	
	;# have all coords, time for some shuffling.

	shufps xmm6, xmm6, 216 ;# 11011000	;# xmm6 = jxH1 jxH2 jyH1 jyH2 
	unpcklps xmm7, xmm2			;# xmm7 = jzH1 jzH2   -    -
	movaps  xmm0, [rsp + nb204_ixM]     
	movaps  xmm1, [rsp + nb204_iyM]
	movaps  xmm2, [rsp + nb204_izM]	
	movlhps xmm3, xmm6			;# xmm3 = jxM   0   jxH1 jxH2 
	shufps  xmm4, xmm6, 228 ;# 11100100	;# xmm4 = jyM   0   jyH1 jyH2 
	shufps  xmm5, xmm7, 68  ;# 01000100	;# xmm5 = jzM   0   jzH1 jzH2

	;# store all j coordinates in jM 
	movaps [rsp + nb204_jxM], xmm3
	movaps [rsp + nb204_jyM], xmm4
	movaps [rsp + nb204_jzM], xmm5
	subps  xmm0, xmm3
	subps  xmm1, xmm4
	subps  xmm2, xmm5
	movaps [rsp + nb204_dxMM], xmm0
	movaps [rsp + nb204_dyMM], xmm1
	movaps [rsp + nb204_dzMM], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 

	movaps xmm6, xmm0
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	mulps   xmm6, [rsp + nb204_krf] ;# xmm6=krsq 
	movaps  xmm2, xmm1
	movaps  xmm7, xmm6 	;# xmm7=krsq 
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb204_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + nb204_half] ;# rinv iO - j water 
	
	addps   xmm6, xmm3	;# xmm6=rinv+ krsq 
	mulps   xmm7, [rsp + nb204_two]
	subps   xmm6, [rsp + nb204_crf] ;# xmm6=rinv+ krsq-crf 
	
	xorps   xmm1, xmm1
	movaps  xmm0, xmm3
	subps   xmm3, xmm7	;# xmm3=rinv-2*krsq 
	xorps   xmm4, xmm4
	mulps   xmm0, xmm0	;# xmm0=rinvsq 
	;# fetch charges to xmm4 (temporary) 
	movss   xmm4, [rsp + nb204_qqMM]
	movhps  xmm4, [rsp + nb204_qqMH]

	mulps xmm6, xmm4	;# vcoul  
	mulps xmm3, xmm4	;# coul part of fs  


	addps   xmm6, [rsp + nb204_vctot]
	mulps   xmm0, xmm3	;# total fscal 
	movaps  [rsp + nb204_vctot], xmm6	

	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	mulps   xmm0, [rsp + nb204_dxMM]
	mulps   xmm1, [rsp + nb204_dyMM]
	mulps   xmm2, [rsp + nb204_dzMM]
	
	;# initial update for j forces 
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [rsp + nb204_fjxM], xmm3
	movaps  [rsp + nb204_fjyM], xmm4
	movaps  [rsp + nb204_fjzM], xmm5
	addps   xmm0, [rsp + nb204_fixM]
	addps   xmm1, [rsp + nb204_fiyM]
	addps   xmm2, [rsp + nb204_fizM]
	movaps  [rsp + nb204_fixM], xmm0
	movaps  [rsp + nb204_fiyM], xmm1
	movaps  [rsp + nb204_fizM], xmm2

	
	;# done with i M Now do i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [rsp + nb204_ixH1]
	movaps  xmm1, [rsp + nb204_iyH1]
	movaps  xmm2, [rsp + nb204_izH1]	
	movaps  xmm3, [rsp + nb204_ixH2] 
	movaps  xmm4, [rsp + nb204_iyH2] 
	movaps  xmm5, [rsp + nb204_izH2] 
	subps   xmm0, [rsp + nb204_jxM]
	subps   xmm1, [rsp + nb204_jyM]
	subps   xmm2, [rsp + nb204_jzM]
	subps   xmm3, [rsp + nb204_jxM]
	subps   xmm4, [rsp + nb204_jyM]
	subps   xmm5, [rsp + nb204_jzM]
	movaps [rsp + nb204_dxH1M], xmm0
	movaps [rsp + nb204_dyH1M], xmm1
	movaps [rsp + nb204_dzH1M], xmm2
	movaps [rsp + nb204_dxH2M], xmm3
	movaps [rsp + nb204_dyH2M], xmm4
	movaps [rsp + nb204_dzH2M], xmm5
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
	movaps  xmm3, [rsp + nb204_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb204_half] ;# rinv H1 - j water 
	mulps   xmm7, [rsp + nb204_half] ;# rinv H2 - j water  

	mulps xmm0, [rsp + nb204_krf] ;# krsq 
	mulps xmm4, [rsp + nb204_krf] ;# krsq  

	;# assemble charges in xmm6 
	xorps   xmm6, xmm6
	movss   xmm6, [rsp + nb204_qqMH]
	movhps  xmm6, [rsp + nb204_qqHH]
	movaps  xmm1, xmm0
	movaps  xmm5, xmm4
	addps   xmm0, xmm3	;# krsq+ rinv 
	addps   xmm4, xmm7	;# krsq+ rinv 
	subps   xmm0, [rsp + nb204_crf]
	subps   xmm4, [rsp + nb204_crf]
	mulps   xmm1, [rsp + nb204_two]
	mulps   xmm5, [rsp + nb204_two]
	mulps   xmm0, xmm6	;# vcoul 
	mulps   xmm4, xmm6	;# vcoul 
	addps   xmm4, xmm0		
	addps   xmm4, [rsp + nb204_vctot]
	movaps  [rsp + nb204_vctot], xmm4
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
	mulps   xmm0, [rsp + nb204_dxH1M]
	mulps   xmm1, [rsp + nb204_dyH1M]
	mulps   xmm2, [rsp + nb204_dzH1M]
	;# update forces H1 - j water 
	movaps  xmm3, [rsp + nb204_fjxM]
	movaps  xmm4, [rsp + nb204_fjyM]
	movaps  xmm5, [rsp + nb204_fjzM]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [rsp + nb204_fjxM], xmm3
	movaps  [rsp + nb204_fjyM], xmm4
	movaps  [rsp + nb204_fjzM], xmm5
	addps   xmm0, [rsp + nb204_fixH1]
	addps   xmm1, [rsp + nb204_fiyH1]
	addps   xmm2, [rsp + nb204_fizH1]
	movaps  [rsp + nb204_fixH1], xmm0
	movaps  [rsp + nb204_fiyH1], xmm1
	movaps  [rsp + nb204_fizH1], xmm2
	;# do forces H2 - j water 
	movaps xmm0, xmm7
	movaps xmm1, xmm7
	movaps xmm2, xmm7
	mulps   xmm0, [rsp + nb204_dxH2M]
	mulps   xmm1, [rsp + nb204_dyH2M]
	mulps   xmm2, [rsp + nb204_dzH2M]
	movaps  xmm3, [rsp + nb204_fjxM]
	movaps  xmm4, [rsp + nb204_fjyM]
	movaps  xmm5, [rsp + nb204_fjzM]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	mov     rsi, [rbp + nb204_faction]
	movaps  [rsp + nb204_fjxM], xmm3
	movaps  [rsp + nb204_fjyM], xmm4
	movaps  [rsp + nb204_fjzM], xmm5
	addps   xmm0, [rsp + nb204_fixH2]
	addps   xmm1, [rsp + nb204_fiyH2]
	addps   xmm2, [rsp + nb204_fizH2]
	movaps  [rsp + nb204_fixH2], xmm0
	movaps  [rsp + nb204_fiyH2], xmm1
	movaps  [rsp + nb204_fizH2], xmm2

	;# update j water forces from local variables 
	movlps  xmm0, [rsi + rax*4 + 36]
	movlps  xmm1, [rsi + rax*4 + 12]
	movhps  xmm1, [rsi + rax*4 + 24]
	movaps  xmm3, [rsp + nb204_fjxM]
	movaps  xmm4, [rsp + nb204_fjyM]
	movaps  xmm5, [rsp + nb204_fjzM]
	movaps  xmm6, xmm5
	movaps  xmm7, xmm5
	shufps  xmm6, xmm6, 2 ;# 00000010
	shufps  xmm7, xmm7, 3 ;# 00000011
	addss   xmm5, [rsi + rax*4 + 44]
	addss   xmm6, [rsi + rax*4 + 20]
	addss   xmm7, [rsi + rax*4 + 32]
	movss   [rsi + rax*4 + 44], xmm5
	movss   [rsi + rax*4 + 20], xmm6
	movss   [rsi + rax*4 + 32], xmm7
	movaps   xmm5, xmm3
	unpcklps xmm3, xmm4
	unpckhps xmm5, xmm4
	addps    xmm0, xmm3
	addps    xmm1, xmm5
	movlps  [rsi + rax*4 + 36], xmm0 
	movlps  [rsi + rax*4 + 12], xmm1 
	movhps  [rsi + rax*4 + 24], xmm1 
	
	dec dword ptr [rsp + nb204_innerk]
	jz    .nb204_updateouterdata
	jmp   .nb204_single_loop
.nb204_updateouterdata:
	mov   ecx, [rsp + nb204_ii3]
	mov   rdi, [rbp + nb204_faction]
	mov   rsi, [rbp + nb204_fshift]
	mov   edx, [rsp + nb204_is3]

	;# accumulate  H1 i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb204_fixH1]
	movaps xmm1, [rsp + nb204_fiyH1] 
	movaps xmm2, [rsp + nb204_fizH1]

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
	movss  xmm3, [rdi + rcx*4 + 12]
	movss  xmm4, [rdi + rcx*4 + 16]
	movss  xmm5, [rdi + rcx*4 + 20]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [rdi + rcx*4 + 12], xmm3
	movss  [rdi + rcx*4 + 16], xmm4
	movss  [rdi + rcx*4 + 20], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# 00001000	
 
	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb204_fixH2]
	movaps xmm1, [rsp + nb204_fiyH2]
	movaps xmm2, [rsp + nb204_fizH2]

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
	movss  xmm3, [rdi + rcx*4 + 24]
	movss  xmm4, [rdi + rcx*4 + 28]
	movss  xmm5, [rdi + rcx*4 + 32]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [rdi + rcx*4 + 24], xmm3
	movss  [rdi + rcx*4 + 28], xmm4
	movss  [rdi + rcx*4 + 32], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# 00001000	
	addps   xmm6, xmm0

	;# accumulate M i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb204_fixM]
	movaps xmm1, [rsp + nb204_fiyM]
	movaps xmm2, [rsp + nb204_fizM]

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
	movss  xmm3, [rdi + rcx*4 + 36]
	movss  xmm4, [rdi + rcx*4 + 40]
	movss  xmm5, [rdi + rcx*4 + 44]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [rdi + rcx*4 + 36], xmm3
	movss  [rdi + rcx*4 + 40], xmm4
	movss  [rdi + rcx*4 + 44], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# 00001000	
	addps   xmm6, xmm0

	;# increment fshift force  
	movlps  xmm3, [rsi + rdx*4]
	movss  xmm4, [rsi + rdx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [rsi + rdx*4],    xmm3
	movss  [rsi + rdx*4 + 8], xmm4

	;# get n from stack
	mov esi, [rsp + nb204_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb204_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb204_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb204_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 	
	
        ;# finish if last 
        mov ecx, [rsp + nb204_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb204_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb204_n], esi
        jmp .nb204_outer
.nb204_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb204_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb204_end
        ;# non-zero, do one more workunit
        jmp   .nb204_threadloop
.nb204_end:

	mov eax, [rsp + nb204_nouter]
	mov ebx, [rsp + nb204_ninner]
	mov rcx, [rbp + nb204_outeriter]
	mov rdx, [rbp + nb204_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 1544
	femms

	pop rbx
	pop	rbp
	ret







	
.globl nb_kernel204nf_x86_64_sse
.globl _nb_kernel204nf_x86_64_sse
nb_kernel204nf_x86_64_sse:	
_nb_kernel204nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb204nf_fshift,         16
.equiv          nb204nf_gid,            24
.equiv          nb204nf_pos,            32
.equiv          nb204nf_faction,        40
.equiv          nb204nf_charge,         48
.equiv          nb204nf_p_facel,        56
.equiv          nb204nf_argkrf,         64
.equiv          nb204nf_argcrf,         72
.equiv          nb204nf_Vc,             80
.equiv          nb204nf_type,           88
.equiv          nb204nf_p_ntype,        96
.equiv          nb204nf_vdwparam,       104
.equiv          nb204nf_Vvdw,           112
.equiv          nb204nf_p_tabscale,     120
.equiv          nb204nf_VFtab,          128
.equiv          nb204nf_invsqrta,       136
.equiv          nb204nf_dvda,           144
.equiv          nb204nf_p_gbtabscale,   152
.equiv          nb204nf_GBtab,          160
.equiv          nb204nf_p_nthreads,     168
.equiv          nb204nf_count,          176
.equiv          nb204nf_mtx,            184
.equiv          nb204nf_outeriter,      192
.equiv          nb204nf_inneriter,      200
.equiv          nb204nf_work,           208
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
.equiv          nb204nf_nri,            720
.equiv          nb204nf_iinr,           728
.equiv          nb204nf_jindex,         736
.equiv          nb204nf_jjnr,           744
.equiv          nb204nf_shift,          752
.equiv          nb204nf_shiftvec,       760
.equiv          nb204nf_facel,          768
.equiv          nb204nf_innerk,         776
.equiv          nb204nf_n,              780
.equiv          nb204nf_nn1,            784
.equiv          nb204nf_nouter,         788
.equiv          nb204nf_ninner,         792

	push rbp
	mov  rbp, rsp
	push rbx

	femms
	sub rsp, 808		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb204nf_nouter], eax
	mov [rsp + nb204nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb204nf_nri], edi
	mov [rsp + nb204nf_iinr], rsi
	mov [rsp + nb204nf_jindex], rdx
	mov [rsp + nb204nf_jjnr], rcx
	mov [rsp + nb204nf_shift], r8
	mov [rsp + nb204nf_shiftvec], r9
	mov rsi, [rbp + nb204nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb204nf_facel], xmm0

	mov rsi, [rbp + nb204nf_argkrf]
	mov rdi, [rbp + nb204nf_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb204nf_krf], xmm1
	movaps [rsp + nb204nf_crf], xmm2

	
	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb204nf_half], eax
	movss xmm1, [rsp + nb204nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb204nf_half],  xmm1
	movaps [rsp + nb204nf_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb204nf_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]		;# ebx =ii 

	mov   rdx, [rbp + nb204nf_charge]
	movss xmm3, [rdx + rbx*4 + 4]	
	movss xmm4, xmm3	
	movss xmm5, [rdx + rbx*4 + 12]	
	mov rsi, [rbp + nb204nf_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb204nf_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb204nf_qqHH], xmm3
	movaps [rsp + nb204nf_qqMH], xmm4
	movaps [rsp + nb204nf_qqMM], xmm5
	
.nb204nf_threadloop:
        mov   rsi, [rbp + nb204nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb204nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   rbx, 1                           ;# rbx=nn1=nn0+10
        lock cmpxchg [rsi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb204nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb204nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb204nf_n], eax
        mov [rsp + nb204nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb204nf_outerstart
        jmp .nb204nf_end
	
.nb204nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb204nf_nouter]
	mov [rsp + nb204nf_nouter], ebx

.nb204nf_outer:
	mov   rax, [rsp + nb204nf_shift]  	;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]	;# ebx=3*is 
	mov   [rsp + nb204nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb204nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb204nf_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   rax, [rbp + nb204nf_pos]	;# rax = base of pos[]  
	mov   [rsp + nb204nf_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [rax + rbx*4 + 12]
	addss xmm4, [rax + rbx*4 + 16]
	addss xmm5, [rax + rbx*4 + 20]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb204nf_ixH1], xmm3
	movaps [rsp + nb204nf_iyH1], xmm4
	movaps [rsp + nb204nf_izH1], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [rax + rbx*4 + 24]
	addss xmm1, [rax + rbx*4 + 28]
	addss xmm2, [rax + rbx*4 + 32]		
	addss xmm3, [rax + rbx*4 + 36]
	addss xmm4, [rax + rbx*4 + 40]
	addss xmm5, [rax + rbx*4 + 44]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb204nf_ixH2], xmm0
	movaps [rsp + nb204nf_iyH2], xmm1
	movaps [rsp + nb204nf_izH2], xmm2
	movaps [rsp + nb204nf_ixM], xmm3
	movaps [rsp + nb204nf_iyM], xmm4
	movaps [rsp + nb204nf_izM], xmm5

	;# clear vctot 
	xorps xmm4, xmm4
	movaps [rsp + nb204nf_vctot], xmm4
	
	mov   rax, [rsp + nb204nf_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb204nf_pos]
	mov   rdi, [rbp + nb204nf_faction]	
	mov   rax, [rsp + nb204nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb204nf_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb204nf_ninner]
	mov   [rsp + nb204nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb204nf_innerk], edx	;# number of innerloop atoms 
	jge   .nb204nf_unroll_loop
	jmp   .nb204nf_single_check
.nb204nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb204nf_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb204nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb204nf_pos]   	;# base of pos[] 

	lea   eax, [eax + eax*2] 	;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2] 	;# replace jnr with j3 
	lea   edx, [edx + edx*2]	
	
	;# move j coordinates to local temp variables 
	movlps xmm2, [rsi + rax*4 + 12]
	movlps xmm3, [rsi + rax*4 + 24]
	movlps xmm4, [rsi + rax*4 + 36]

	movlps xmm5, [rsi + rbx*4 + 12]
	movlps xmm6, [rsi + rbx*4 + 24]
	movlps xmm7, [rsi + rbx*4 + 36]

	movhps xmm2, [rsi + rcx*4 + 12]
	movhps xmm3, [rsi + rcx*4 + 24]
	movhps xmm4, [rsi + rcx*4 + 36]

	movhps xmm5, [rsi + rdx*4 + 12]
	movhps xmm6, [rsi + rdx*4 + 24]
	movhps xmm7, [rsi + rdx*4 + 36]
	
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
	movaps [rsp + nb204nf_jxH1], xmm0
	movhlps  xmm2, xmm6
	movaps [rsp + nb204nf_jyH1], xmm2
	movlhps  xmm1, xmm3
	movaps [rsp + nb204nf_jxH2], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [rsp + nb204nf_jyH2], xmm3
	movlhps  xmm4, xmm5
	movaps [rsp + nb204nf_jxM], xmm4
	movhlps  xmm5, xmm6
	movaps [rsp + nb204nf_jyM], xmm5

	movss  xmm0, [rsi + rax*4 + 20]
	movss  xmm1, [rsi + rax*4 + 32]
	movss  xmm2, [rsi + rax*4 + 44]

	movss  xmm3, [rsi + rcx*4 + 20]
	movss  xmm4, [rsi + rcx*4 + 32]
	movss  xmm5, [rsi + rcx*4 + 44]

	movhps xmm0, [rsi + rbx*4 + 16]
	movhps xmm1, [rsi + rbx*4 + 28]
	movhps xmm2, [rsi + rbx*4 + 40]
	
	movhps xmm3, [rsi + rdx*4 + 16]
	movhps xmm4, [rsi + rdx*4 + 28]
	movhps xmm5, [rsi + rdx*4 + 40]
	
	shufps xmm0, xmm3, 204  ;# 11001100
	shufps xmm1, xmm4, 204  ;# 11001100
	shufps xmm2, xmm5, 204  ;# 11001100
	movaps [rsp + nb204nf_jzH1],  xmm0
	movaps [rsp + nb204nf_jzH2],  xmm1
	movaps [rsp + nb204nf_jzM],  xmm2

	movaps xmm0, [rsp + nb204nf_ixH1]
	movaps xmm1, [rsp + nb204nf_iyH1]
	movaps xmm2, [rsp + nb204nf_izH1]
	movaps xmm3, [rsp + nb204nf_ixH1]
	movaps xmm4, [rsp + nb204nf_iyH1]
	movaps xmm5, [rsp + nb204nf_izH1]
	subps  xmm0, [rsp + nb204nf_jxH1]
	subps  xmm1, [rsp + nb204nf_jyH1]
	subps  xmm2, [rsp + nb204nf_jzH1]
	subps  xmm3, [rsp + nb204nf_jxH2]
	subps  xmm4, [rsp + nb204nf_jyH2]
	subps  xmm5, [rsp + nb204nf_jzH2]
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
	movaps [rsp + nb204nf_rsqH1H1], xmm0
	movaps [rsp + nb204nf_rsqH1H2], xmm3

	movaps xmm0, [rsp + nb204nf_ixH1]
	movaps xmm1, [rsp + nb204nf_iyH1]
	movaps xmm2, [rsp + nb204nf_izH1]
	movaps xmm3, [rsp + nb204nf_ixH2]
	movaps xmm4, [rsp + nb204nf_iyH2]
	movaps xmm5, [rsp + nb204nf_izH2]
	subps  xmm0, [rsp + nb204nf_jxM]
	subps  xmm1, [rsp + nb204nf_jyM]
	subps  xmm2, [rsp + nb204nf_jzM]
	subps  xmm3, [rsp + nb204nf_jxH1]
	subps  xmm4, [rsp + nb204nf_jyH1]
	subps  xmm5, [rsp + nb204nf_jzH1]
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
	movaps [rsp + nb204nf_rsqH1M], xmm0
	movaps [rsp + nb204nf_rsqH2H1], xmm3

	movaps xmm0, [rsp + nb204nf_ixH2]
	movaps xmm1, [rsp + nb204nf_iyH2]
	movaps xmm2, [rsp + nb204nf_izH2]
	movaps xmm3, [rsp + nb204nf_ixH2]
	movaps xmm4, [rsp + nb204nf_iyH2]
	movaps xmm5, [rsp + nb204nf_izH2]
	subps  xmm0, [rsp + nb204nf_jxH2]
	subps  xmm1, [rsp + nb204nf_jyH2]
	subps  xmm2, [rsp + nb204nf_jzH2]
	subps  xmm3, [rsp + nb204nf_jxM]
	subps  xmm4, [rsp + nb204nf_jyM]
	subps  xmm5, [rsp + nb204nf_jzM]
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
	movaps [rsp + nb204nf_rsqH2H2], xmm0
	movaps [rsp + nb204nf_rsqH2M], xmm3

	movaps xmm0, [rsp + nb204nf_ixM]
	movaps xmm1, [rsp + nb204nf_iyM]
	movaps xmm2, [rsp + nb204nf_izM]
	movaps xmm3, [rsp + nb204nf_ixM]
	movaps xmm4, [rsp + nb204nf_iyM]
	movaps xmm5, [rsp + nb204nf_izM]
	subps  xmm0, [rsp + nb204nf_jxH1]
	subps  xmm1, [rsp + nb204nf_jyH1]
	subps  xmm2, [rsp + nb204nf_jzH1]
	subps  xmm3, [rsp + nb204nf_jxH2]
	subps  xmm4, [rsp + nb204nf_jyH2]
	subps  xmm5, [rsp + nb204nf_jzH2]
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
	movaps [rsp + nb204nf_rsqMH1], xmm0
	movaps [rsp + nb204nf_rsqMH2], xmm4

	movaps xmm0, [rsp + nb204nf_ixM]
	movaps xmm1, [rsp + nb204nf_iyM]
	movaps xmm2, [rsp + nb204nf_izM]
	subps  xmm0, [rsp + nb204nf_jxM]
	subps  xmm1, [rsp + nb204nf_jyM]
	subps  xmm2, [rsp + nb204nf_jzM]
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [rsp + nb204nf_rsqMM], xmm0
	
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb204nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb204nf_half] ;# rinvH2H2 
	mulps   xmm7, [rsp + nb204nf_half] ;# rinvH2H1 
	movaps  [rsp + nb204nf_rinvMM], xmm3
	movaps  [rsp + nb204nf_rinvMH2], xmm7
	
	rsqrtps xmm1, [rsp + nb204nf_rsqH1H1]
	rsqrtps xmm5, [rsp + nb204nf_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb204nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb204nf_rsqH1H1]
	mulps   xmm5, [rsp + nb204nf_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb204nf_half] 
	mulps   xmm7, [rsp + nb204nf_half]
	movaps  [rsp + nb204nf_rinvH1H1], xmm3
	movaps  [rsp + nb204nf_rinvH1H2], xmm7
	
	rsqrtps xmm1, [rsp + nb204nf_rsqH1M]
	rsqrtps xmm5, [rsp + nb204nf_rsqH2H1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb204nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb204nf_rsqH1M]
	mulps   xmm5, [rsp + nb204nf_rsqH2H1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb204nf_half] 
	mulps   xmm7, [rsp + nb204nf_half]
	movaps  [rsp + nb204nf_rinvH1M], xmm3
	movaps  [rsp + nb204nf_rinvH2H1], xmm7
	
	rsqrtps xmm1, [rsp + nb204nf_rsqH2H2]
	rsqrtps xmm5, [rsp + nb204nf_rsqH2M]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb204nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb204nf_rsqH2H2]
	mulps   xmm5, [rsp + nb204nf_rsqH2M]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb204nf_half] 
	mulps   xmm7, [rsp + nb204nf_half]
	movaps  [rsp + nb204nf_rinvH2H2], xmm3
	movaps  [rsp + nb204nf_rinvH2M], xmm7
	
	rsqrtps xmm1, [rsp + nb204nf_rsqMH1]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb204nf_three]
	mulps   xmm1, [rsp + nb204nf_rsqMH1]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + nb204nf_half] 
	movaps  [rsp + nb204nf_rinvMH1], xmm3

	;# Coulomb interactions 
	;# add all H-H rsq in xmm2, H-M rsq xmm4
	;# H-H rinv in xmm0, H-M in xmm1
	movaps xmm0, [rsp + nb204nf_rinvH1H1]
	movaps xmm1, [rsp + nb204nf_rinvH1M]
	movaps xmm2, [rsp + nb204nf_rsqH1H1]
	movaps xmm4, [rsp + nb204nf_rsqH1M]
	addps  xmm0, [rsp + nb204nf_rinvH1H2]
	addps  xmm1, [rsp + nb204nf_rinvH2M]
	addps  xmm2, [rsp + nb204nf_rsqH1H2]
	addps  xmm4, [rsp + nb204nf_rsqH2M]
	addps  xmm0, [rsp + nb204nf_rinvH2H1]
	addps  xmm1, [rsp + nb204nf_rinvMH1]
	addps  xmm2, [rsp + nb204nf_rsqH2H1]
	addps  xmm4, [rsp + nb204nf_rsqMH1]
	addps  xmm0, [rsp + nb204nf_rinvH2H2]
	addps  xmm1, [rsp + nb204nf_rinvMH2]
	addps  xmm2, [rsp + nb204nf_rsqH2H2]
	addps  xmm4, [rsp + nb204nf_rsqMH2]
	movaps xmm5, [rsp + nb204nf_krf]
	movaps xmm6, [rsp + nb204nf_crf]

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
	mulps  xmm0, [rsp + nb204nf_qqHH]
	mulps  xmm1, [rsp + nb204nf_qqMH]
	addps  xmm0, xmm1 
	addps  xmm0, [rsp + nb204nf_vctot]
	;# M-M interaction
	movaps xmm4, [rsp + nb204nf_rinvMM]
	mulps  xmm5, [rsp + nb204nf_rsqMM] ;# krsq
	addps  xmm5, xmm4                  ;# rinv+krsq
	subps xmm5, [rsp + nb204nf_crf] ;# xmm5=rinv+ krsq-crf 
	mulps xmm5, [rsp + nb204nf_qqMM]
	addps xmm5, xmm0
	movaps [rsp + nb204nf_vctot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb204nf_innerk],  4
	jl    .nb204nf_single_check
	jmp   .nb204nf_unroll_loop
.nb204nf_single_check:
	add dword ptr [rsp + nb204nf_innerk],  4
	jnz   .nb204nf_single_loop
	jmp   .nb204nf_updateouterdata
.nb204nf_single_loop:
	mov   rdx, [rsp + nb204nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb204nf_innerjjnr],  4	

	mov rsi, [rbp + nb204nf_pos]
	lea   eax, [eax + eax*2]  

	;# fetch j coordinates 
	xorps xmm3, xmm3
	xorps xmm4, xmm4
	xorps xmm5, xmm5
	movss xmm3, [rsi + rax*4 + 36]		;# jxM  -  -  -
	movss xmm4, [rsi + rax*4 + 40]		;# jyM  -  -  -
	movss xmm5, [rsi + rax*4 + 44]		;# jzM  -  -  -  

	movlps xmm6, [rsi + rax*4 + 12]		;# xmm6 = jxH1 jyH1   -    -
	movss  xmm7, [rsi + rax*4 + 20]		;# xmm7 = jzH1   -    -    - 
	movhps xmm6, [rsi + rax*4 + 24]		;# xmm6 = jxH1 jyH1 jxH2 jyH2
	movss  xmm2, [rsi + rax*4 + 32]		;# xmm2 = jzH2   -    -    -
	
	;# have all coords, time for some shuffling.

	shufps xmm6, xmm6, 216 ;# 11011000	;# xmm6 = jxH1 jxH2 jyH1 jyH2 
	unpcklps xmm7, xmm2			;# xmm7 = jzH1 jzH2   -    -
	movaps  xmm0, [rsp + nb204nf_ixM]     
	movaps  xmm1, [rsp + nb204nf_iyM]
	movaps  xmm2, [rsp + nb204nf_izM]	
	movlhps xmm3, xmm6			;# xmm3 = jxM   0   jxH1 jxH2 
	shufps  xmm4, xmm6, 228 ;# 11100100	;# xmm4 = jyM   0   jyH1 jyH2 
	shufps  xmm5, xmm7, 68  ;# 01000100	;# xmm5 = jzM   0   jzH1 jzH2

	;# store all j coordinates in jM
	movaps [rsp + nb204nf_jxM], xmm3
	movaps [rsp + nb204nf_jyM], xmm4
	movaps [rsp + nb204nf_jzM], xmm5
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
	mulps   xmm6, [rsp + nb204nf_krf] ;# xmm6=krsq 
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb204nf_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + nb204nf_half] ;# rinv iO - j water 
	
	addps   xmm6, xmm3	;# xmm6=rinv+ krsq 
	subps   xmm6, [rsp + nb204nf_crf] ;# xmm6=rinv+ krsq-crf 
	
	xorps   xmm4, xmm4
	;# fetch charges to xmm4 (temporary) 
	movss   xmm4, [rsp + nb204nf_qqMM]
	movhps  xmm4, [rsp + nb204nf_qqMH]
	mulps xmm6, xmm4	;# vcoul  
	addps   xmm6, [rsp + nb204nf_vctot]
	movaps  [rsp + nb204nf_vctot], xmm6	
	
	;# done with i M Now do i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [rsp + nb204nf_ixH1]
	movaps  xmm1, [rsp + nb204nf_iyH1]
	movaps  xmm2, [rsp + nb204nf_izH1]	
	movaps  xmm3, [rsp + nb204nf_ixH2] 
	movaps  xmm4, [rsp + nb204nf_iyH2] 
	movaps  xmm5, [rsp + nb204nf_izH2] 
	subps   xmm0, [rsp + nb204nf_jxM]
	subps   xmm1, [rsp + nb204nf_jyM]
	subps   xmm2, [rsp + nb204nf_jzM]
	subps   xmm3, [rsp + nb204nf_jxM]
	subps   xmm4, [rsp + nb204nf_jyM]
	subps   xmm5, [rsp + nb204nf_jzM]
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
	movaps  xmm3, [rsp + nb204nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb204nf_half] ;# rinv H1 - j water 
	mulps   xmm7, [rsp + nb204nf_half] ;# rinv H2 - j water  

	mulps xmm0, [rsp + nb204nf_krf] ;# krsq 
	mulps xmm4, [rsp + nb204nf_krf] ;# krsq  

	;# assemble charges in xmm6 
	xorps   xmm6, xmm6
	movss   xmm6, [rsp + nb204nf_qqMH]
	movhps  xmm6, [rsp + nb204nf_qqHH]
	addps   xmm0, xmm3	;# krsq+ rinv 
	addps   xmm4, xmm7	;# krsq+ rinv 
	subps   xmm0, [rsp + nb204nf_crf]
	subps   xmm4, [rsp + nb204nf_crf]
	mulps   xmm0, xmm6	;# vcoul 
	mulps   xmm4, xmm6	;# vcoul 
	addps   xmm4, xmm0		
	addps   xmm4, [rsp + nb204nf_vctot]
	movaps  [rsp + nb204nf_vctot], xmm4
	dec dword ptr [rsp + nb204nf_innerk]
	jz    .nb204nf_updateouterdata
	jmp   .nb204nf_single_loop
.nb204nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb204nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb204nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb204nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb204nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 	
	
        ;# finish if last 
        mov ecx, [rsp + nb204nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb204nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb204nf_n], esi
        jmp .nb204nf_outer
.nb204nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb204nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb204nf_end
        ;# non-zero, do one more workunit
        jmp   .nb204nf_threadloop
.nb204nf_end:

	mov eax, [rsp + nb204nf_nouter]
	mov ebx, [rsp + nb204nf_ninner]
	mov rcx, [rbp + nb204nf_outeriter]
	mov rdx, [rbp + nb204nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 808
	femms

	pop rbx
	pop	rbp
	ret
