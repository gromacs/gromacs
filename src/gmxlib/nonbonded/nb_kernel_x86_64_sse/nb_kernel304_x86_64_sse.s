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


	

	
.globl nb_kernel304_x86_64_sse
.globl _nb_kernel304_x86_64_sse
nb_kernel304_x86_64_sse:	
_nb_kernel304_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb304_fshift,           16
.equiv          nb304_gid,              24
.equiv          nb304_pos,              32
.equiv          nb304_faction,          40
.equiv          nb304_charge,           48
.equiv          nb304_p_facel,          56
.equiv          nb304_argkrf,           64
.equiv          nb304_argcrf,           72
.equiv          nb304_Vc,               80
.equiv          nb304_type,             88
.equiv          nb304_p_ntype,          96
.equiv          nb304_vdwparam,         104
.equiv          nb304_Vvdw,             112
.equiv          nb304_p_tabscale,       120
.equiv          nb304_VFtab,            128
.equiv          nb304_invsqrta,         136
.equiv          nb304_dvda,             144
.equiv          nb304_p_gbtabscale,     152
.equiv          nb304_GBtab,            160
.equiv          nb304_p_nthreads,       168
.equiv          nb304_count,            176
.equiv          nb304_mtx,              184
.equiv          nb304_outeriter,        192
.equiv          nb304_inneriter,        200
.equiv          nb304_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb304_ixH1,             0
.equiv          nb304_iyH1,             16
.equiv          nb304_izH1,             32
.equiv          nb304_ixH2,             48
.equiv          nb304_iyH2,             64
.equiv          nb304_izH2,             80
.equiv          nb304_ixM,              96
.equiv          nb304_iyM,              112
.equiv          nb304_izM,              128
.equiv          nb304_jxH1,             144
.equiv          nb304_jyH1,             160
.equiv          nb304_jzH1,             176
.equiv          nb304_jxH2,             192
.equiv          nb304_jyH2,             208
.equiv          nb304_jzH2,             224
.equiv          nb304_jxM,              240
.equiv          nb304_jyM,              256
.equiv          nb304_jzM,              272
.equiv          nb304_dxH1H1,           288
.equiv          nb304_dyH1H1,           304
.equiv          nb304_dzH1H1,           320
.equiv          nb304_dxH1H2,           336
.equiv          nb304_dyH1H2,           352
.equiv          nb304_dzH1H2,           368
.equiv          nb304_dxH1M,            384
.equiv          nb304_dyH1M,            400
.equiv          nb304_dzH1M,            416
.equiv          nb304_dxH2H1,           432
.equiv          nb304_dyH2H1,           448
.equiv          nb304_dzH2H1,           464
.equiv          nb304_dxH2H2,           480
.equiv          nb304_dyH2H2,           496
.equiv          nb304_dzH2H2,           512
.equiv          nb304_dxH2M,            528
.equiv          nb304_dyH2M,            544
.equiv          nb304_dzH2M,            560
.equiv          nb304_dxMH1,            576
.equiv          nb304_dyMH1,            592
.equiv          nb304_dzMH1,            608
.equiv          nb304_dxMH2,            624
.equiv          nb304_dyMH2,            640
.equiv          nb304_dzMH2,            656
.equiv          nb304_dxMM,             672
.equiv          nb304_dyMM,             688
.equiv          nb304_dzMM,             704
.equiv          nb304_qqHH,             720
.equiv          nb304_qqMH,             736
.equiv          nb304_qqMM,             752
.equiv          nb304_two,              768
.equiv          nb304_tsc,              784
.equiv          nb304_vctot,            800
.equiv          nb304_fixH1,            816
.equiv          nb304_fiyH1,            832
.equiv          nb304_fizH1,            848
.equiv          nb304_fixH2,            864
.equiv          nb304_fiyH2,            880
.equiv          nb304_fizH2,            896
.equiv          nb304_fixM,             912
.equiv          nb304_fiyM,             928
.equiv          nb304_fizM,             944
.equiv          nb304_fjxH1,            960
.equiv          nb304_fjyH1,            976
.equiv          nb304_fjzH1,            992
.equiv          nb304_fjxH2,            1008
.equiv          nb304_fjyH2,            1024
.equiv          nb304_fjzH2,            1040
.equiv          nb304_fjxM,             1056
.equiv          nb304_fjyM,             1072
.equiv          nb304_fjzM,             1088
.equiv          nb304_half,             1104
.equiv          nb304_three,            1120
.equiv          nb304_rsqH1H1,          1136
.equiv          nb304_rsqH1H2,          1152
.equiv          nb304_rsqH1M,           1168
.equiv          nb304_rsqH2H1,          1184
.equiv          nb304_rsqH2H2,          1200
.equiv          nb304_rsqH2M,           1216
.equiv          nb304_rsqMH1,           1232
.equiv          nb304_rsqMH2,           1248
.equiv          nb304_rsqMM,            1264
.equiv          nb304_rinvH1H1,         1280
.equiv          nb304_rinvH1H2,         1296
.equiv          nb304_rinvH1M,          1312
.equiv          nb304_rinvH2H1,         1328
.equiv          nb304_rinvH2H2,         1344
.equiv          nb304_rinvH2M,          1360
.equiv          nb304_rinvMH1,          1376
.equiv          nb304_rinvMH2,          1392
.equiv          nb304_rinvMM,           1408
.equiv          nb304_is3,              1424
.equiv          nb304_ii3,              1428
.equiv          nb304_nri,              1432
.equiv          nb304_iinr,             1440
.equiv          nb304_jindex,           1448
.equiv          nb304_jjnr,             1456
.equiv          nb304_shift,            1464
.equiv          nb304_shiftvec,         1472
.equiv          nb304_facel,            1480
.equiv          nb304_innerjjnr,        1488
.equiv          nb304_innerk,           1496
.equiv          nb304_n,                1500
.equiv          nb304_nn1,              1504
.equiv          nb304_nouter,           1508
.equiv          nb304_ninner,           1512
	push rbp
	mov  rbp, rsp
	push rbx

	femms
	sub rsp, 1528		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb304_nouter], eax
	mov [rsp + nb304_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb304_nri], edi
	mov [rsp + nb304_iinr], rsi
	mov [rsp + nb304_jindex], rdx
	mov [rsp + nb304_jjnr], rcx
	mov [rsp + nb304_shift], r8
	mov [rsp + nb304_shiftvec], r9
	mov rsi, [rbp + nb304_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb304_facel], xmm0

	mov rax, [rbp + nb304_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb304_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb304_half], eax
	movss xmm1, [rsp + nb304_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb304_half],  xmm1
	movaps [rsp + nb304_two],  xmm2
	movaps [rsp + nb304_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb304_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]		;# ebx =ii 

	mov   rdx, [rbp + nb304_charge]
	movss xmm3, [rdx + rbx*4 + 4]	
	movss xmm4, xmm3	
	movss xmm5, [rdx + rbx*4 + 12]	
	mov rsi, [rbp + nb304_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb304_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb304_qqHH], xmm3
	movaps [rsp + nb304_qqMH], xmm4
	movaps [rsp + nb304_qqMM], xmm5		

.nb304_threadloop:
        mov   rsi, [rbp + nb304_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb304_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [rsi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb304_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb304_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb304_n], eax
        mov [rsp + nb304_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb304_outerstart
        jmp .nb304_end
	
.nb304_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb304_nouter]
	mov [rsp + nb304_nouter], ebx

.nb304_outer:
	mov   rax, [rsp + nb304_shift]  	;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb304_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb304_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb304_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb304_pos]	;# rax = base of pos[]  
	mov   [rsp + nb304_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [rax + rbx*4 + 12]
	addss xmm4, [rax + rbx*4 + 16]
	addss xmm5, [rax + rbx*4 + 20]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb304_ixH1], xmm3
	movaps [rsp + nb304_iyH1], xmm4
	movaps [rsp + nb304_izH1], xmm5

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
	movaps [rsp + nb304_ixH2], xmm0
	movaps [rsp + nb304_iyH2], xmm1
	movaps [rsp + nb304_izH2], xmm2
	movaps [rsp + nb304_ixM], xmm3
	movaps [rsp + nb304_iyM], xmm4
	movaps [rsp + nb304_izM], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb304_vctot], xmm4
	movaps [rsp + nb304_fixH1], xmm4
	movaps [rsp + nb304_fiyH1], xmm4
	movaps [rsp + nb304_fizH1], xmm4
	movaps [rsp + nb304_fixH2], xmm4
	movaps [rsp + nb304_fiyH2], xmm4
	movaps [rsp + nb304_fizH2], xmm4
	movaps [rsp + nb304_fixM], xmm4
	movaps [rsp + nb304_fiyM], xmm4
	movaps [rsp + nb304_fizM], xmm4
	
	mov   rax, [rsp + nb304_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb304_pos]
	mov   rdi, [rbp + nb304_faction]	
	mov   rax, [rsp + nb304_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb304_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb304_ninner]
	mov   [rsp + nb304_ninner], ecx
	add   edx, 0
	mov   [rsp + nb304_innerk], edx	;# number of innerloop atoms 
	jge   .nb304_unroll_loop
	jmp   .nb304_single_check
.nb304_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb304_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb304_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb304_pos]   	;# base of pos[] 

	lea   rax, [rax + rax*2] 	;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2] 	;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	
	
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
	movaps [rsp + nb304_jxH1], xmm0
	movhlps  xmm2, xmm6
	movaps [rsp + nb304_jyH1], xmm2
	movlhps  xmm1, xmm3
	movaps [rsp + nb304_jxH2], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [rsp + nb304_jyH2], xmm3
	movlhps  xmm4, xmm5
	movaps [rsp + nb304_jxM], xmm4
	movhlps  xmm5, xmm6
	movaps [rsp + nb304_jyM], xmm5

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
	movaps [rsp + nb304_jzH1],  xmm0
	movaps [rsp + nb304_jzH2],  xmm1
	movaps [rsp + nb304_jzM],  xmm2

	movaps xmm0, [rsp + nb304_ixH1]
	movaps xmm1, [rsp + nb304_iyH1]
	movaps xmm2, [rsp + nb304_izH1]
	movaps xmm3, [rsp + nb304_ixH1]
	movaps xmm4, [rsp + nb304_iyH1]
	movaps xmm5, [rsp + nb304_izH1]
	subps  xmm0, [rsp + nb304_jxH1]
	subps  xmm1, [rsp + nb304_jyH1]
	subps  xmm2, [rsp + nb304_jzH1]
	subps  xmm3, [rsp + nb304_jxH2]
	subps  xmm4, [rsp + nb304_jyH2]
	subps  xmm5, [rsp + nb304_jzH2]
	movaps [rsp + nb304_dxH1H1], xmm0
	movaps [rsp + nb304_dyH1H1], xmm1
	movaps [rsp + nb304_dzH1H1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb304_dxH1H2], xmm3
	movaps [rsp + nb304_dyH1H2], xmm4
	movaps [rsp + nb304_dzH1H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + nb304_rsqH1H1], xmm0
	movaps [rsp + nb304_rsqH1H2], xmm3

	movaps xmm0, [rsp + nb304_ixH1]
	movaps xmm1, [rsp + nb304_iyH1]
	movaps xmm2, [rsp + nb304_izH1]
	movaps xmm3, [rsp + nb304_ixH2]
	movaps xmm4, [rsp + nb304_iyH2]
	movaps xmm5, [rsp + nb304_izH2]
	subps  xmm0, [rsp + nb304_jxM]
	subps  xmm1, [rsp + nb304_jyM]
	subps  xmm2, [rsp + nb304_jzM]
	subps  xmm3, [rsp + nb304_jxH1]
	subps  xmm4, [rsp + nb304_jyH1]
	subps  xmm5, [rsp + nb304_jzH1]
	movaps [rsp + nb304_dxH1M], xmm0
	movaps [rsp + nb304_dyH1M], xmm1
	movaps [rsp + nb304_dzH1M], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb304_dxH2H1], xmm3
	movaps [rsp + nb304_dyH2H1], xmm4
	movaps [rsp + nb304_dzH2H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + nb304_rsqH1M], xmm0
	movaps [rsp + nb304_rsqH2H1], xmm3

	movaps xmm0, [rsp + nb304_ixH2]
	movaps xmm1, [rsp + nb304_iyH2]
	movaps xmm2, [rsp + nb304_izH2]
	movaps xmm3, [rsp + nb304_ixH2]
	movaps xmm4, [rsp + nb304_iyH2]
	movaps xmm5, [rsp + nb304_izH2]
	subps  xmm0, [rsp + nb304_jxH2]
	subps  xmm1, [rsp + nb304_jyH2]
	subps  xmm2, [rsp + nb304_jzH2]
	subps  xmm3, [rsp + nb304_jxM]
	subps  xmm4, [rsp + nb304_jyM]
	subps  xmm5, [rsp + nb304_jzM]
	movaps [rsp + nb304_dxH2H2], xmm0
	movaps [rsp + nb304_dyH2H2], xmm1
	movaps [rsp + nb304_dzH2H2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb304_dxH2M], xmm3
	movaps [rsp + nb304_dyH2M], xmm4
	movaps [rsp + nb304_dzH2M], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + nb304_rsqH2H2], xmm0
	movaps [rsp + nb304_rsqH2M], xmm3

	movaps xmm0, [rsp + nb304_ixM]
	movaps xmm1, [rsp + nb304_iyM]
	movaps xmm2, [rsp + nb304_izM]
	movaps xmm3, [rsp + nb304_ixM]
	movaps xmm4, [rsp + nb304_iyM]
	movaps xmm5, [rsp + nb304_izM]
	subps  xmm0, [rsp + nb304_jxH1]
	subps  xmm1, [rsp + nb304_jyH1]
	subps  xmm2, [rsp + nb304_jzH1]
	subps  xmm3, [rsp + nb304_jxH2]
	subps  xmm4, [rsp + nb304_jyH2]
	subps  xmm5, [rsp + nb304_jzH2]
	movaps [rsp + nb304_dxMH1], xmm0
	movaps [rsp + nb304_dyMH1], xmm1
	movaps [rsp + nb304_dzMH1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb304_dxMH2], xmm3
	movaps [rsp + nb304_dyMH2], xmm4
	movaps [rsp + nb304_dzMH2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [rsp + nb304_rsqMH1], xmm0
	movaps [rsp + nb304_rsqMH2], xmm4

	movaps xmm0, [rsp + nb304_ixM]
	movaps xmm1, [rsp + nb304_iyM]
	movaps xmm2, [rsp + nb304_izM]
	subps  xmm0, [rsp + nb304_jxM]
	subps  xmm1, [rsp + nb304_jyM]
	subps  xmm2, [rsp + nb304_jzM]
	movaps [rsp + nb304_dxMM], xmm0
	movaps [rsp + nb304_dyMM], xmm1
	movaps [rsp + nb304_dzMM], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [rsp + nb304_rsqMM], xmm0
	
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb304_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb304_half] ;# rinvMM
	mulps   xmm7, [rsp + nb304_half] ;# rinvMH2 
	movaps  [rsp + nb304_rinvMM], xmm3
	movaps  [rsp + nb304_rinvMH2], xmm7
	
	rsqrtps xmm1, [rsp + nb304_rsqH1H1]
	rsqrtps xmm5, [rsp + nb304_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb304_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb304_rsqH1H1]
	mulps   xmm5, [rsp + nb304_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb304_half] 
	mulps   xmm7, [rsp + nb304_half]
	movaps  [rsp + nb304_rinvH1H1], xmm3
	movaps  [rsp + nb304_rinvH1H2], xmm7
	
	rsqrtps xmm1, [rsp + nb304_rsqH1M]
	rsqrtps xmm5, [rsp + nb304_rsqH2H1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb304_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb304_rsqH1M]
	mulps   xmm5, [rsp + nb304_rsqH2H1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb304_half] 
	mulps   xmm7, [rsp + nb304_half]
	movaps  [rsp + nb304_rinvH1M], xmm3
	movaps  [rsp + nb304_rinvH2H1], xmm7
	
	rsqrtps xmm1, [rsp + nb304_rsqH2H2]
	rsqrtps xmm5, [rsp + nb304_rsqH2M]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb304_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb304_rsqH2H2]
	mulps   xmm5, [rsp + nb304_rsqH2M]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb304_half] 
	mulps   xmm7, [rsp + nb304_half]
	movaps  [rsp + nb304_rinvH2H2], xmm3
	movaps  [rsp + nb304_rinvH2M], xmm7
	
	rsqrtps xmm1, [rsp + nb304_rsqMH1]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb304_three]
	mulps   xmm1, [rsp + nb304_rsqMH1]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + nb304_half] 
	movaps  [rsp + nb304_rinvMH1], xmm3

	;# start with H1-H1 interaction 
	movaps xmm0, [rsp + nb304_rinvH1H1]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb304_rsqH1H1] ;# xmm1=r 
	mulps  xmm1, [rsp + nb304_tsc]
	
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

	mov  rsi, [rbp + nb304_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [rsp + nb304_two]   	;# two*Heps2 
	movaps xmm3, [rsp + nb304_qqHH]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addps  xmm5, [rsp + nb304_vctot]
	xorps  xmm2, xmm2
	movaps [rsp + nb304_vctot], xmm5
	mulps  xmm3, [rsp + nb304_tsc]
	
	subps  xmm2, xmm3
	mulps  xmm0, xmm2
	
	movaps xmm1, xmm0
	movaps xmm2, xmm0		

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [rsp + nb304_dxH1H1]
	mulps xmm1, [rsp + nb304_dyH1H1]
	mulps xmm2, [rsp + nb304_dzH1H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb304_fixH1]
	addps xmm1, [rsp + nb304_fiyH1]
	addps xmm2, [rsp + nb304_fizH1]
	movaps [rsp + nb304_fjxH1], xmm3
	movaps [rsp + nb304_fjyH1], xmm4
	movaps [rsp + nb304_fjzH1], xmm5
	movaps [rsp + nb304_fixH1], xmm0
	movaps [rsp + nb304_fiyH1], xmm1
	movaps [rsp + nb304_fizH1], xmm2

	;# H1-H2 interaction 
	movaps xmm0, [rsp + nb304_rinvH1H2]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb304_rsqH1H2] ;# xmm1=r 
	mulps  xmm1, [rsp + nb304_tsc]	
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

	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [rsp + nb304_two]   	;# two*Heps2 
	movaps xmm3, [rsp + nb304_qqHH]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 

	addps  xmm5, [rsp + nb304_vctot]
	movaps [rsp + nb304_vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [rsp + nb304_tsc]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [rsp + nb304_dxH1H2]
	mulps xmm1, [rsp + nb304_dyH1H2]
	mulps xmm2, [rsp + nb304_dzH1H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb304_fixH1]
	addps xmm1, [rsp + nb304_fiyH1]
	addps xmm2, [rsp + nb304_fizH1]
	movaps [rsp + nb304_fjxH2], xmm3
	movaps [rsp + nb304_fjyH2], xmm4
	movaps [rsp + nb304_fjzH2], xmm5
	movaps [rsp + nb304_fixH1], xmm0
	movaps [rsp + nb304_fiyH1], xmm1
	movaps [rsp + nb304_fizH1], xmm2

	;# H1-M interaction  
	movaps xmm0, [rsp + nb304_rinvH1M]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb304_rsqH1M] ;# xmm1=r 
	mulps  xmm1, [rsp + nb304_tsc]	
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

	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [rsp + nb304_two]   	;# two*Heps2 
	movaps xmm3, [rsp + nb304_qqMH]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 

	addps  xmm5, [rsp + nb304_vctot]
	movaps [rsp + nb304_vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [rsp + nb304_tsc]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [rsp + nb304_dxH1M]
	mulps xmm1, [rsp + nb304_dyH1M]
	mulps xmm2, [rsp + nb304_dzH1M]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb304_fixH1]
	addps xmm1, [rsp + nb304_fiyH1]
	addps xmm2, [rsp + nb304_fizH1]
	movaps [rsp + nb304_fjxM], xmm3
	movaps [rsp + nb304_fjyM], xmm4
	movaps [rsp + nb304_fjzM], xmm5
	movaps [rsp + nb304_fixH1], xmm0
	movaps [rsp + nb304_fiyH1], xmm1
	movaps [rsp + nb304_fizH1], xmm2

	;# H2-H1 interaction 
	movaps xmm0, [rsp + nb304_rinvH2H1]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb304_rsqH2H1] ;# xmm1=r 
	mulps  xmm1, [rsp + nb304_tsc]	
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

	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [rsp + nb304_two]   	;# two*Heps2 
	movaps xmm3, [rsp + nb304_qqHH]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 

	addps  xmm5, [rsp + nb304_vctot]
	movaps [rsp + nb304_vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [rsp + nb304_tsc]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [rsp + nb304_fjxH1]
	movaps xmm4, [rsp + nb304_fjyH1]
	movaps xmm5, [rsp + nb304_fjzH1]
	mulps xmm0, [rsp + nb304_dxH2H1]
	mulps xmm1, [rsp + nb304_dyH2H1]
	mulps xmm2, [rsp + nb304_dzH2H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb304_fixH2]
	addps xmm1, [rsp + nb304_fiyH2]
	addps xmm2, [rsp + nb304_fizH2]
	movaps [rsp + nb304_fjxH1], xmm3
	movaps [rsp + nb304_fjyH1], xmm4
	movaps [rsp + nb304_fjzH1], xmm5
	movaps [rsp + nb304_fixH2], xmm0
	movaps [rsp + nb304_fiyH2], xmm1
	movaps [rsp + nb304_fizH2], xmm2

	;# H2-H2 interaction 
	movaps xmm0, [rsp + nb304_rinvH2H2]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb304_rsqH2H2] ;# xmm1=r 
	mulps  xmm1, [rsp + nb304_tsc]	
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

	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [rsp + nb304_two]   	;# two*Heps2 
	movaps xmm3, [rsp + nb304_qqHH]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 

	addps  xmm5, [rsp + nb304_vctot]
	movaps [rsp + nb304_vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [rsp + nb304_tsc]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [rsp + nb304_fjxH2]
	movaps xmm4, [rsp + nb304_fjyH2]
	movaps xmm5, [rsp + nb304_fjzH2]
	mulps xmm0, [rsp + nb304_dxH2H2]
	mulps xmm1, [rsp + nb304_dyH2H2]
	mulps xmm2, [rsp + nb304_dzH2H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb304_fixH2]
	addps xmm1, [rsp + nb304_fiyH2]
	addps xmm2, [rsp + nb304_fizH2]
	movaps [rsp + nb304_fjxH2], xmm3
	movaps [rsp + nb304_fjyH2], xmm4
	movaps [rsp + nb304_fjzH2], xmm5
	movaps [rsp + nb304_fixH2], xmm0
	movaps [rsp + nb304_fiyH2], xmm1
	movaps [rsp + nb304_fizH2], xmm2

	;# H2-M interaction 
	movaps xmm0, [rsp + nb304_rinvH2M]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb304_rsqH2M] ;# xmm1=r 
	mulps  xmm1, [rsp + nb304_tsc]
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

	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [rsp + nb304_two]   	;# two*Heps2 
	movaps xmm3, [rsp + nb304_qqMH]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 

	addps  xmm5, [rsp + nb304_vctot]
	movaps [rsp + nb304_vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [rsp + nb304_tsc]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [rsp + nb304_fjxM]
	movaps xmm4, [rsp + nb304_fjyM]
	movaps xmm5, [rsp + nb304_fjzM]
	mulps xmm0, [rsp + nb304_dxH2M]
	mulps xmm1, [rsp + nb304_dyH2M]
	mulps xmm2, [rsp + nb304_dzH2M]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb304_fixH2]
	addps xmm1, [rsp + nb304_fiyH2]
	addps xmm2, [rsp + nb304_fizH2]
	movaps [rsp + nb304_fjxM], xmm3
	movaps [rsp + nb304_fjyM], xmm4
	movaps [rsp + nb304_fjzM], xmm5
	movaps [rsp + nb304_fixH2], xmm0
	movaps [rsp + nb304_fiyH2], xmm1
	movaps [rsp + nb304_fizH2], xmm2

	;# M-H1 interaction 
	movaps xmm0, [rsp + nb304_rinvMH1]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb304_rsqMH1] ;# xmm1=r 
	mulps  xmm1, [rsp + nb304_tsc]	
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

	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [rsp + nb304_two]   	;# two*Heps2 
	movaps xmm3, [rsp + nb304_qqMH]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 

	addps  xmm5, [rsp + nb304_vctot]
	movaps [rsp + nb304_vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [rsp + nb304_tsc]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1

	movaps xmm3, [rsp + nb304_fjxH1]
	movaps xmm4, [rsp + nb304_fjyH1]
	movaps xmm5, [rsp + nb304_fjzH1]
	mulps xmm0, [rsp + nb304_dxMH1]
	mulps xmm1, [rsp + nb304_dyMH1]
	mulps xmm2, [rsp + nb304_dzMH1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb304_fixM]
	addps xmm1, [rsp + nb304_fiyM]
	addps xmm2, [rsp + nb304_fizM]
	movaps [rsp + nb304_fjxH1], xmm3
	movaps [rsp + nb304_fjyH1], xmm4
	movaps [rsp + nb304_fjzH1], xmm5
	movaps [rsp + nb304_fixM], xmm0
	movaps [rsp + nb304_fiyM], xmm1
	movaps [rsp + nb304_fizM], xmm2

	;# M-H2 interaction 
	movaps xmm0, [rsp + nb304_rinvMH2]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb304_rsqMH2] ;# xmm1=r 
	mulps  xmm1, [rsp + nb304_tsc]
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

	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [rsp + nb304_two]   	;# two*Heps2 
	movaps xmm3, [rsp + nb304_qqMH]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 

	addps  xmm5, [rsp + nb304_vctot]
	movaps [rsp + nb304_vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [rsp + nb304_tsc]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [rsp + nb304_fjxH2]
	movaps xmm4, [rsp + nb304_fjyH2]
	movaps xmm5, [rsp + nb304_fjzH2]
	mulps xmm0, [rsp + nb304_dxMH2]
	mulps xmm1, [rsp + nb304_dyMH2]
	mulps xmm2, [rsp + nb304_dzMH2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb304_fixM]
	addps xmm1, [rsp + nb304_fiyM]
	addps xmm2, [rsp + nb304_fizM]
	movaps [rsp + nb304_fjxH2], xmm3
	movaps [rsp + nb304_fjyH2], xmm4
	movaps [rsp + nb304_fjzH2], xmm5
	movaps [rsp + nb304_fixM], xmm0
	movaps [rsp + nb304_fiyM], xmm1
	movaps [rsp + nb304_fizM], xmm2

	;# M-M interaction 
	movaps xmm0, [rsp + nb304_rinvMM]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb304_rsqMM] ;# xmm1=r 
	mulps  xmm1, [rsp + nb304_tsc]	
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

	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [rsp + nb304_two]   	;# two*Heps2 
	movaps xmm3, [rsp + nb304_qqMM]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 

	addps  xmm5, [rsp + nb304_vctot]
	movaps [rsp + nb304_vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [rsp + nb304_tsc]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [rsp + nb304_fjxM]
	movaps xmm4, [rsp + nb304_fjyM]
	movaps xmm5, [rsp + nb304_fjzM]
	mulps xmm0, [rsp + nb304_dxMM]
	mulps xmm1, [rsp + nb304_dyMM]
	mulps xmm2, [rsp + nb304_dzMM]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb304_fixM]
	addps xmm1, [rsp + nb304_fiyM]
	addps xmm2, [rsp + nb304_fizM]
	movaps [rsp + nb304_fjxM], xmm3
	movaps [rsp + nb304_fjyM], xmm4
	movaps [rsp + nb304_fjzM], xmm5
	movaps [rsp + nb304_fixM], xmm0
	movaps [rsp + nb304_fiyM], xmm1
	movaps [rsp + nb304_fizM], xmm2

	mov rdi, [rbp + nb304_faction]

	movd eax, mm0
	movd ebx, mm1
	movd ecx, mm2
	movd edx, mm3
	
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
	movaps xmm0, [rsp + nb304_fjxH1] ;# xmm0= fjxH1a  fjxH1b  fjxH1c  fjxH1d 
	movaps xmm2, [rsp + nb304_fjyH1] ;# xmm1= fjyH1a  fjyH1b  fjyH1c  fjyH1d
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
	movaps xmm0, [rsp + nb304_fjzH1]  ;# xmm0= fjzH1a   fjzH1b   fjzH1c   fjzH1d 
	movaps xmm2, [rsp + nb304_fjxH2] ;# xmm1= fjxH2a  fjxH2b  fjxH2c  fjxH2d
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
	movaps xmm0, [rsp + nb304_fjyH2] ;# xmm0= fjyH2a  fjyH2b  fjyH2c  fjyH2d 
	movaps xmm2, [rsp + nb304_fjzH2] ;# xmm1= fjzH2a  fjzH2b  fjzH2c  fjzH2d
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
	movaps xmm0, [rsp + nb304_fjxM] ;# xmm0= fjxMa  fjxMb  fjxMc  fjxMd 
	movaps xmm2, [rsp + nb304_fjyM] ;# xmm1= fjyMa  fjyMb  fjyMc  fjyMd
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
	addss xmm4, [rsp + nb304_fjzM] 
	addss xmm5, [rsp + nb304_fjzM + 4] 
	addss xmm6, [rsp + nb304_fjzM + 8] 
	addss xmm7, [rsp + nb304_fjzM + 12]
	;# store back
	movss [rdi + rax*4 + 44], xmm4
	movss [rdi + rbx*4 + 44], xmm5
	movss [rdi + rcx*4 + 44], xmm6
	movss [rdi + rdx*4 + 44], xmm7
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb304_innerk],  4
	jl    .nb304_single_check
	jmp   .nb304_unroll_loop
.nb304_single_check:
	add dword ptr [rsp + nb304_innerk],  4
	jnz   .nb304_single_loop
	jmp   .nb304_updateouterdata
.nb304_single_loop:
	mov   rdx, [rsp + nb304_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb304_innerjjnr],  4	

	mov rsi, [rbp + nb304_pos]
	lea   rax, [rax + rax*2]  

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
	movaps  xmm0, [rsp + nb304_ixM]     
	movaps  xmm1, [rsp + nb304_iyM]
	movaps  xmm2, [rsp + nb304_izM]	
	movlhps xmm3, xmm6			;# xmm3 = jxM   0   jxH1 jxH2 
	shufps  xmm4, xmm6, 228 ;# 11100100	;# xmm4 = jyM   0   jyH1 jyH2 
	shufps  xmm5, xmm7, 68  ;# 01000100	;# xmm5 = jzM   0   jzH1 jzH2
	
	;# store all j coordinates in jM
	movaps [rsp + nb304_jxM], xmm3
	movaps [rsp + nb304_jyM], xmm4
	movaps [rsp + nb304_jzM], xmm5
	subps  xmm0, xmm3
	subps  xmm1, xmm4
	subps  xmm2, xmm5
	movaps [rsp + nb304_dxMM], xmm0
	movaps [rsp + nb304_dyMM], xmm1
	movaps [rsp + nb304_dzMM], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	movaps  xmm2, xmm1	
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb304_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [rsp + nb304_half] ;# rinv iO - j water 

	movaps  xmm1, xmm3
	mulps   xmm1, xmm0	;# xmm1=r 
	movaps  xmm0, xmm3	;# xmm0=rinv 
	mulps  xmm1, [rsp + nb304_tsc]
	
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
	
	movd ebx, mm6
	movd ecx, mm7
	psrlq mm7, 32
	movd edx, mm7		;# table indices in ebx,ecx,edx 

	mov rsi, [rbp + nb304_VFtab]
	
	movlps xmm5, [rsi + rbx*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rbx*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  
	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [rsp + nb304_two]   	;# two*Heps2 

	xorps  xmm3, xmm3

	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [rsp + nb304_qqMM]
	movhps  xmm3, [rsp + nb304_qqMH]
	
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point xmm5 contains vcoul and xmm3 fijC 
	
	addps  xmm5, [rsp + nb304_vctot]
	movaps [rsp + nb304_vctot], xmm5
	xorps  xmm2, xmm2
	mulps  xmm3, [rsp + nb304_tsc]

	subps  xmm2, xmm3
	mulps  xmm0, xmm2
	
	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [rsp + nb304_dxMM]
	mulps   xmm1, [rsp + nb304_dyMM]
	mulps   xmm2, [rsp + nb304_dzMM]
	;# initial update for j forces 
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [rsp + nb304_fjxM], xmm3
	movaps  [rsp + nb304_fjyM], xmm4
	movaps  [rsp + nb304_fjzM], xmm5
	addps   xmm0, [rsp + nb304_fixM]
	addps   xmm1, [rsp + nb304_fiyM]
	addps   xmm2, [rsp + nb304_fizM]
	movaps  [rsp + nb304_fixM], xmm0
	movaps  [rsp + nb304_fiyM], xmm1
	movaps  [rsp + nb304_fizM], xmm2

	
	;# done with i M Now do i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [rsp + nb304_ixH1]
	movaps  xmm1, [rsp + nb304_iyH1]
	movaps  xmm2, [rsp + nb304_izH1]	
	movaps  xmm3, [rsp + nb304_ixH2] 
	movaps  xmm4, [rsp + nb304_iyH2] 
	movaps  xmm5, [rsp + nb304_izH2] 
	subps   xmm0, [rsp + nb304_jxM]
	subps   xmm1, [rsp + nb304_jyM]
	subps   xmm2, [rsp + nb304_jzM]
	subps   xmm3, [rsp + nb304_jxM]
	subps   xmm4, [rsp + nb304_jyM]
	subps   xmm5, [rsp + nb304_jzM]
	movaps [rsp + nb304_dxH1M], xmm0
	movaps [rsp + nb304_dyH1M], xmm1
	movaps [rsp + nb304_dzH1M], xmm2
	movaps [rsp + nb304_dxH2M], xmm3
	movaps [rsp + nb304_dyH2M], xmm4
	movaps [rsp + nb304_dzH2M], xmm5
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

	;# start with H1, save H2 data 
	movaps [rsp + nb304_rsqH2M], xmm4
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb304_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb304_half] ;# rinv H1 - j water 
	mulps   xmm7, [rsp + nb304_half] ;# rinv H2 - j water  

	;# start with H1, save H2 data 
	movaps [rsp + nb304_rinvH2M], xmm7

	movaps xmm1, xmm3
	mulps  xmm1, xmm0	;# xmm1=r 
	movaps xmm0, xmm3	;# xmm0=rinv 
	mulps  xmm1, [rsp + nb304_tsc]
	
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

	movd ebx, mm6
	movd ecx, mm7
	psrlq mm7, 32
	movd edx, mm7		;# table indices in ebx,ecx,edx 

	movlps xmm5, [rsi + rbx*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rbx*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  
	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [rsp + nb304_two]   	;# two*Heps2 

	xorps  xmm3, xmm3
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [rsp + nb304_qqMH]
	movhps  xmm3, [rsp + nb304_qqHH]
	
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point xmm5 contains vcoul and xmm3 fijC 
	addps  xmm5, [rsp + nb304_vctot]
	movaps [rsp + nb304_vctot], xmm5	

	xorps  xmm1, xmm1

	mulps xmm3, [rsp + nb304_tsc]
	mulps xmm3, xmm0
	subps  xmm1, xmm3
	
	movaps  xmm0, xmm1
	movaps  xmm2, xmm1
	mulps   xmm0, [rsp + nb304_dxH1M]
	mulps   xmm1, [rsp + nb304_dyH1M]
	mulps   xmm2, [rsp + nb304_dzH1M]
	;# update forces H1 - j water 
	movaps  xmm3, [rsp + nb304_fjxM]
	movaps  xmm4, [rsp + nb304_fjyM]
	movaps  xmm5, [rsp + nb304_fjzM]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [rsp + nb304_fjxM], xmm3
	movaps  [rsp + nb304_fjyM], xmm4
	movaps  [rsp + nb304_fjzM], xmm5
	addps   xmm0, [rsp + nb304_fixH1]
	addps   xmm1, [rsp + nb304_fiyH1]
	addps   xmm2, [rsp + nb304_fizH1]
	movaps  [rsp + nb304_fixH1], xmm0
	movaps  [rsp + nb304_fiyH1], xmm1
	movaps  [rsp + nb304_fizH1], xmm2
	;# do table for H2 - j water interaction 
	movaps xmm0, [rsp + nb304_rinvH2M]
	movaps xmm1, [rsp + nb304_rsqH2M]
	mulps  xmm1, xmm0	;# xmm0=rinv, xmm1=r 
	mulps  xmm1, [rsp + nb304_tsc]
	
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

	movd ebx, mm6
	movd ecx, mm7
	psrlq mm7, 32
	movd edx, mm7		;# table indices in ebx,ecx,edx 

	movlps xmm5, [rsi + rbx*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rbx*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  
	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [rsp + nb304_two]   	;# two*Heps2 

	xorps  xmm3, xmm3
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [rsp + nb304_qqMH]
	movhps  xmm3, [rsp + nb304_qqHH]
	
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point xmm5 contains vcoul and xmm3 fijC 
	addps  xmm5, [rsp + nb304_vctot]
	movaps [rsp + nb304_vctot], xmm5	

	xorps  xmm1, xmm1

	mulps xmm3, [rsp + nb304_tsc]
	mulps xmm3, xmm0
	subps  xmm1, xmm3
	
	movaps  xmm0, xmm1
	movaps  xmm2, xmm1
	
	mulps   xmm0, [rsp + nb304_dxH2M]
	mulps   xmm1, [rsp + nb304_dyH2M]
	mulps   xmm2, [rsp + nb304_dzH2M]
	movaps  xmm3, [rsp + nb304_fjxM]
	movaps  xmm4, [rsp + nb304_fjyM]
	movaps  xmm5, [rsp + nb304_fjzM]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	mov     rsi, [rbp + nb304_faction]
	movaps  [rsp + nb304_fjxM], xmm3
	movaps  [rsp + nb304_fjyM], xmm4
	movaps  [rsp + nb304_fjzM], xmm5
	addps   xmm0, [rsp + nb304_fixH2]
	addps   xmm1, [rsp + nb304_fiyH2]
	addps   xmm2, [rsp + nb304_fizH2]
	movaps  [rsp + nb304_fixH2], xmm0
	movaps  [rsp + nb304_fiyH2], xmm1
	movaps  [rsp + nb304_fizH2], xmm2

	;# update j water forces from local variables 
	movlps  xmm0, [rsi + rax*4 + 36]
	movlps  xmm1, [rsi + rax*4 + 12]
	movhps  xmm1, [rsi + rax*4 + 24]
	movaps  xmm3, [rsp + nb304_fjxM]
	movaps  xmm4, [rsp + nb304_fjyM]
	movaps  xmm5, [rsp + nb304_fjzM]
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
	
	dec dword ptr [rsp + nb304_innerk]
	jz    .nb304_updateouterdata
	jmp   .nb304_single_loop
.nb304_updateouterdata:
	mov   ecx, [rsp + nb304_ii3]
	mov   rdi, [rbp + nb304_faction]
	mov   rsi, [rbp + nb304_fshift]
	mov   edx, [rsp + nb304_is3]

	;# accumulate  H1 i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb304_fixH1]
	movaps xmm1, [rsp + nb304_fiyH1] 
	movaps xmm2, [rsp + nb304_fizH1]

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

	;# accumulate H2 i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb304_fixH2]
	movaps xmm1, [rsp + nb304_fiyH2]
	movaps xmm2, [rsp + nb304_fizH2]

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
	movaps xmm0, [rsp + nb304_fixM]
	movaps xmm1, [rsp + nb304_fiyM]
	movaps xmm2, [rsp + nb304_fizM]

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
	mov esi, [rsp + nb304_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb304_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb304_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb304_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb304_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb304_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb304_n], esi
        jmp .nb304_outer
.nb304_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb304_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb304_end
        ;# non-zero, do one more workunit
        jmp   .nb304_threadloop
.nb304_end:
	mov eax, [rsp + nb304_nouter]
	mov ebx, [rsp + nb304_ninner]
	mov rcx, [rbp + nb304_outeriter]
	mov rdx, [rbp + nb304_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 1528
	femms

	pop rbx
	pop	rbp
	ret

	

	

	
.globl nb_kernel304nf_x86_64_sse
.globl _nb_kernel304nf_x86_64_sse
nb_kernel304nf_x86_64_sse:	
_nb_kernel304nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb304nf_fshift,         16
.equiv          nb304nf_gid,            24
.equiv          nb304nf_pos,            32
.equiv          nb304nf_faction,        40
.equiv          nb304nf_charge,         48
.equiv          nb304nf_p_facel,        56
.equiv          nb304nf_argkrf,         64
.equiv          nb304nf_argcrf,         72
.equiv          nb304nf_Vc,             80
.equiv          nb304nf_type,           88
.equiv          nb304nf_p_ntype,        96
.equiv          nb304nf_vdwparam,       104
.equiv          nb304nf_Vvdw,           112
.equiv          nb304nf_p_tabscale,     120
.equiv          nb304nf_VFtab,          128
.equiv          nb304nf_invsqrta,       136
.equiv          nb304nf_dvda,           144
.equiv          nb304nf_p_gbtabscale,   152
.equiv          nb304nf_GBtab,          160
.equiv          nb304nf_p_nthreads,     168
.equiv          nb304nf_count,          176
.equiv          nb304nf_mtx,            184
.equiv          nb304nf_outeriter,      192
.equiv          nb304nf_inneriter,      200
.equiv          nb304nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb304nf_ixH1,           0
.equiv          nb304nf_iyH1,           16
.equiv          nb304nf_izH1,           32
.equiv          nb304nf_ixH2,           48
.equiv          nb304nf_iyH2,           64
.equiv          nb304nf_izH2,           80
.equiv          nb304nf_ixM,            96
.equiv          nb304nf_iyM,            112
.equiv          nb304nf_izM,            128
.equiv          nb304nf_jxH1,           144
.equiv          nb304nf_jyH1,           160
.equiv          nb304nf_jzH1,           176
.equiv          nb304nf_jxH2,           192
.equiv          nb304nf_jyH2,           208
.equiv          nb304nf_jzH2,           224
.equiv          nb304nf_jxM,            240
.equiv          nb304nf_jyM,            256
.equiv          nb304nf_jzM,            272
.equiv          nb304nf_qqHH,           288
.equiv          nb304nf_qqMH,           304
.equiv          nb304nf_qqMM,           320
.equiv          nb304nf_tsc,            336
.equiv          nb304nf_vctot,          352
.equiv          nb304nf_half,           368
.equiv          nb304nf_three,          384
.equiv          nb304nf_rsqH1H1,        400
.equiv          nb304nf_rsqH1H2,        416
.equiv          nb304nf_rsqH1M,         432
.equiv          nb304nf_rsqH2H1,        448
.equiv          nb304nf_rsqH2H2,        464
.equiv          nb304nf_rsqH2M,         480
.equiv          nb304nf_rsqMH1,         496
.equiv          nb304nf_rsqMH2,         512
.equiv          nb304nf_rsqMM,          528
.equiv          nb304nf_rinvH1H1,       544
.equiv          nb304nf_rinvH1H2,       560
.equiv          nb304nf_rinvH1M,        576
.equiv          nb304nf_rinvH2H1,       592
.equiv          nb304nf_rinvH2H2,       608
.equiv          nb304nf_rinvH2M,        624
.equiv          nb304nf_rinvMH1,        640
.equiv          nb304nf_rinvMH2,        656
.equiv          nb304nf_rinvMM,         672
.equiv          nb304nf_is3,            688
.equiv          nb304nf_ii3,            692
.equiv          nb304nf_nri,            696
.equiv          nb304nf_iinr,           704
.equiv          nb304nf_jindex,         712
.equiv          nb304nf_jjnr,           720
.equiv          nb304nf_shift,          728
.equiv          nb304nf_shiftvec,       736
.equiv          nb304nf_facel,          744
.equiv          nb304nf_innerjjnr,      752
.equiv          nb304nf_innerk,         760
.equiv          nb304nf_n,              764
.equiv          nb304nf_nn1,            768
.equiv          nb304nf_nouter,         772
.equiv          nb304nf_ninner,         776
	push rbp
	mov  rbp, rsp
	push rbx

	femms
	sub rsp, 792		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb304nf_nouter], eax
	mov [rsp + nb304nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb304nf_nri], edi
	mov [rsp + nb304nf_iinr], rsi
	mov [rsp + nb304nf_jindex], rdx
	mov [rsp + nb304nf_jjnr], rcx
	mov [rsp + nb304nf_shift], r8
	mov [rsp + nb304nf_shiftvec], r9
	mov rsi, [rbp + nb304nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb304nf_facel], xmm0

	mov rax, [rbp + nb304nf_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb304nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb304nf_half], eax
	movss xmm1, [rsp + nb304nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb304nf_half],  xmm1
	movaps [rsp + nb304nf_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb304nf_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]		;# ebx =ii 

	mov   rdx, [rbp + nb304nf_charge]
	movss xmm3, [rdx + rbx*4 + 4]	
	movss xmm4, xmm3	
	movss xmm5, [rdx + rbx*4 + 12]	
	mov rsi, [rbp + nb304nf_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb304nf_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb304nf_qqHH], xmm3
	movaps [rsp + nb304nf_qqMH], xmm4
	movaps [rsp + nb304nf_qqMM], xmm5		

.nb304nf_threadloop:
        mov   rsi, [rbp + nb304nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb304nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [rsi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb304nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb304nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb304nf_n], eax
        mov [rsp + nb304nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb304nf_outerstart
        jmp .nb304nf_end
	
.nb304nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb304nf_nouter]
	mov [rsp + nb304nf_nouter], ebx

.nb304nf_outer:
	mov   rax, [rsp + nb304nf_shift]  	;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb304nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb304nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb304nf_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb304nf_pos]	;# rax = base of pos[]  
	mov   [rsp + nb304nf_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [rax + rbx*4 + 12]
	addss xmm4, [rax + rbx*4 + 16]
	addss xmm5, [rax + rbx*4 + 20]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb304nf_ixH1], xmm3
	movaps [rsp + nb304nf_iyH1], xmm4
	movaps [rsp + nb304nf_izH1], xmm5

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
	movaps [rsp + nb304nf_ixH2], xmm0
	movaps [rsp + nb304nf_iyH2], xmm1
	movaps [rsp + nb304nf_izH2], xmm2
	movaps [rsp + nb304nf_ixM], xmm3
	movaps [rsp + nb304nf_iyM], xmm4
	movaps [rsp + nb304nf_izM], xmm5

	;# clear vctot 
	xorps xmm4, xmm4
	movaps [rsp + nb304nf_vctot], xmm4
	
	mov   rax, [rsp + nb304nf_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb304nf_pos]
	mov   rdi, [rbp + nb304nf_faction]	
	mov   rax, [rsp + nb304nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb304nf_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb304nf_ninner]
	mov   [rsp + nb304nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb304nf_innerk], edx	;# number of innerloop atoms 
	jge   .nb304nf_unroll_loop
	jmp   .nb304nf_single_check
.nb304nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb304nf_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb304nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb304nf_pos]   	;# base of pos[] 

	lea   rax, [rax + rax*2] 	;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2] 	;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	
	
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
	movaps [rsp + nb304nf_jxH1], xmm0
	movhlps  xmm2, xmm6
	movaps [rsp + nb304nf_jyH1], xmm2
	movlhps  xmm1, xmm3
	movaps [rsp + nb304nf_jxH2], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [rsp + nb304nf_jyH2], xmm3
	movlhps  xmm4, xmm5
	movaps [rsp + nb304nf_jxM], xmm4
	movhlps  xmm5, xmm6
	movaps [rsp + nb304nf_jyM], xmm5

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
	movaps [rsp + nb304nf_jzH1],  xmm0
	movaps [rsp + nb304nf_jzH2],  xmm1
	movaps [rsp + nb304nf_jzM],  xmm2

	movaps xmm0, [rsp + nb304nf_ixH1]
	movaps xmm1, [rsp + nb304nf_iyH1]
	movaps xmm2, [rsp + nb304nf_izH1]
	movaps xmm3, [rsp + nb304nf_ixH1]
	movaps xmm4, [rsp + nb304nf_iyH1]
	movaps xmm5, [rsp + nb304nf_izH1]
	subps  xmm0, [rsp + nb304nf_jxH1]
	subps  xmm1, [rsp + nb304nf_jyH1]
	subps  xmm2, [rsp + nb304nf_jzH1]
	subps  xmm3, [rsp + nb304nf_jxH2]
	subps  xmm4, [rsp + nb304nf_jyH2]
	subps  xmm5, [rsp + nb304nf_jzH2]
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
	movaps [rsp + nb304nf_rsqH1H1], xmm0
	movaps [rsp + nb304nf_rsqH1H2], xmm3

	movaps xmm0, [rsp + nb304nf_ixH1]
	movaps xmm1, [rsp + nb304nf_iyH1]
	movaps xmm2, [rsp + nb304nf_izH1]
	movaps xmm3, [rsp + nb304nf_ixH2]
	movaps xmm4, [rsp + nb304nf_iyH2]
	movaps xmm5, [rsp + nb304nf_izH2]
	subps  xmm0, [rsp + nb304nf_jxM]
	subps  xmm1, [rsp + nb304nf_jyM]
	subps  xmm2, [rsp + nb304nf_jzM]
	subps  xmm3, [rsp + nb304nf_jxH1]
	subps  xmm4, [rsp + nb304nf_jyH1]
	subps  xmm5, [rsp + nb304nf_jzH1]
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
	movaps [rsp + nb304nf_rsqH1M], xmm0
	movaps [rsp + nb304nf_rsqH2H1], xmm3

	movaps xmm0, [rsp + nb304nf_ixH2]
	movaps xmm1, [rsp + nb304nf_iyH2]
	movaps xmm2, [rsp + nb304nf_izH2]
	movaps xmm3, [rsp + nb304nf_ixH2]
	movaps xmm4, [rsp + nb304nf_iyH2]
	movaps xmm5, [rsp + nb304nf_izH2]
	subps  xmm0, [rsp + nb304nf_jxH2]
	subps  xmm1, [rsp + nb304nf_jyH2]
	subps  xmm2, [rsp + nb304nf_jzH2]
	subps  xmm3, [rsp + nb304nf_jxM]
	subps  xmm4, [rsp + nb304nf_jyM]
	subps  xmm5, [rsp + nb304nf_jzM]
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
	movaps [rsp + nb304nf_rsqH2H2], xmm0
	movaps [rsp + nb304nf_rsqH2M], xmm3

	movaps xmm0, [rsp + nb304nf_ixM]
	movaps xmm1, [rsp + nb304nf_iyM]
	movaps xmm2, [rsp + nb304nf_izM]
	movaps xmm3, [rsp + nb304nf_ixM]
	movaps xmm4, [rsp + nb304nf_iyM]
	movaps xmm5, [rsp + nb304nf_izM]
	subps  xmm0, [rsp + nb304nf_jxH1]
	subps  xmm1, [rsp + nb304nf_jyH1]
	subps  xmm2, [rsp + nb304nf_jzH1]
	subps  xmm3, [rsp + nb304nf_jxH2]
	subps  xmm4, [rsp + nb304nf_jyH2]
	subps  xmm5, [rsp + nb304nf_jzH2]
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
	movaps [rsp + nb304nf_rsqMH1], xmm0
	movaps [rsp + nb304nf_rsqMH2], xmm4

	movaps xmm0, [rsp + nb304nf_ixM]
	movaps xmm1, [rsp + nb304nf_iyM]
	movaps xmm2, [rsp + nb304nf_izM]
	subps  xmm0, [rsp + nb304nf_jxM]
	subps  xmm1, [rsp + nb304nf_jyM]
	subps  xmm2, [rsp + nb304nf_jzM]
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [rsp + nb304nf_rsqMM], xmm0
	
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb304nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb304nf_half] ;# rinvMM
	mulps   xmm7, [rsp + nb304nf_half] ;# rinvMH2 
	movaps  [rsp + nb304nf_rinvMM], xmm3
	movaps  [rsp + nb304nf_rinvMH2], xmm7
	
	rsqrtps xmm1, [rsp + nb304nf_rsqH1H1]
	rsqrtps xmm5, [rsp + nb304nf_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb304nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb304nf_rsqH1H1]
	mulps   xmm5, [rsp + nb304nf_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb304nf_half] 
	mulps   xmm7, [rsp + nb304nf_half]
	movaps  [rsp + nb304nf_rinvH1H1], xmm3
	movaps  [rsp + nb304nf_rinvH1H2], xmm7
	
	rsqrtps xmm1, [rsp + nb304nf_rsqH1M]
	rsqrtps xmm5, [rsp + nb304nf_rsqH2H1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb304nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb304nf_rsqH1M]
	mulps   xmm5, [rsp + nb304nf_rsqH2H1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb304nf_half] 
	mulps   xmm7, [rsp + nb304nf_half]
	movaps  [rsp + nb304nf_rinvH1M], xmm3
	movaps  [rsp + nb304nf_rinvH2H1], xmm7
	
	rsqrtps xmm1, [rsp + nb304nf_rsqH2H2]
	rsqrtps xmm5, [rsp + nb304nf_rsqH2M]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb304nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb304nf_rsqH2H2]
	mulps   xmm5, [rsp + nb304nf_rsqH2M]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb304nf_half] 
	mulps   xmm7, [rsp + nb304nf_half]
	movaps  [rsp + nb304nf_rinvH2H2], xmm3
	movaps  [rsp + nb304nf_rinvH2M], xmm7
	
	rsqrtps xmm1, [rsp + nb304nf_rsqMH1]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb304nf_three]
	mulps   xmm1, [rsp + nb304nf_rsqMH1]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + nb304nf_half] 
	movaps  [rsp + nb304nf_rinvMH1], xmm3

	;# start with H1-H1 interaction 
	movaps xmm0, [rsp + nb304nf_rinvH1H1]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb304nf_rsqH1H1] ;# xmm1=r 
	mulps  xmm1, [rsp + nb304nf_tsc]
	
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

	mov  rsi, [rbp + nb304nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	movaps xmm3, [rsp + nb304nf_qqHH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 
	;# update vctot 
	addps  xmm5, [rsp + nb304nf_vctot]
	movaps [rsp + nb304nf_vctot], xmm5

	;# H1-H2 interaction 
	movaps xmm0, [rsp + nb304nf_rinvH1H2]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb304nf_rsqH1H2] ;# xmm1=r 
	mulps  xmm1, [rsp + nb304nf_tsc]	
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

	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	movaps xmm3, [rsp + nb304nf_qqHH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 

	addps  xmm5, [rsp + nb304nf_vctot]
	movaps [rsp + nb304nf_vctot], xmm5

	;# H1-M interaction  
	movaps xmm0, [rsp + nb304nf_rinvH1M]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb304nf_rsqH1M] ;# xmm1=r 
	mulps  xmm1, [rsp + nb304nf_tsc]	
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

	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	movaps xmm3, [rsp + nb304nf_qqMH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 

	addps  xmm5, [rsp + nb304nf_vctot]
	movaps [rsp + nb304nf_vctot], xmm5

	;# H2-H1 interaction 
	movaps xmm0, [rsp + nb304nf_rinvH2H1]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb304nf_rsqH2H1] ;# xmm1=r 
	mulps  xmm1, [rsp + nb304nf_tsc]	
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

	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	movaps xmm3, [rsp + nb304nf_qqHH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 

	addps  xmm5, [rsp + nb304nf_vctot]
	movaps [rsp + nb304nf_vctot], xmm5

	;# H2-H2 interaction 
	movaps xmm0, [rsp + nb304nf_rinvH2H2]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb304nf_rsqH2H2] ;# xmm1=r 
	mulps  xmm1, [rsp + nb304nf_tsc]	
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

	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	movaps xmm3, [rsp + nb304nf_qqHH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 

	addps  xmm5, [rsp + nb304nf_vctot]
	movaps [rsp + nb304nf_vctot], xmm5

	;# H2-M interaction 
	movaps xmm0, [rsp + nb304nf_rinvH2M]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb304nf_rsqH2M] ;# xmm1=r 
	mulps  xmm1, [rsp + nb304nf_tsc]
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

	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	movaps xmm3, [rsp + nb304nf_qqMH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 
	addps  xmm5, [rsp + nb304nf_vctot]
	movaps [rsp + nb304nf_vctot], xmm5

	;# M-H1 interaction 
	movaps xmm0, [rsp + nb304nf_rinvMH1]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb304nf_rsqMH1] ;# xmm1=r 
	mulps  xmm1, [rsp + nb304nf_tsc]	
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

	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	movaps xmm3, [rsp + nb304nf_qqMH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 

	addps  xmm5, [rsp + nb304nf_vctot]
	movaps [rsp + nb304nf_vctot], xmm5

	;# M-H2 interaction 
	movaps xmm0, [rsp + nb304nf_rinvMH2]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb304nf_rsqMH2] ;# xmm1=r 
	mulps  xmm1, [rsp + nb304nf_tsc]
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

	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	movaps xmm3, [rsp + nb304nf_qqMH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 

	addps  xmm5, [rsp + nb304nf_vctot]
	movaps [rsp + nb304nf_vctot], xmm5

	;# M-M interaction 
	movaps xmm0, [rsp + nb304nf_rinvMM]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb304nf_rsqMM] ;# xmm1=r 
	mulps  xmm1, [rsp + nb304nf_tsc]	
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

	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	movaps xmm3, [rsp + nb304nf_qqMM]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 

	addps  xmm5, [rsp + nb304nf_vctot]
	movaps [rsp + nb304nf_vctot], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb304nf_innerk],  4
	jl    .nb304nf_single_check
	jmp   .nb304nf_unroll_loop
.nb304nf_single_check:
	add dword ptr [rsp + nb304nf_innerk],  4
	jnz   .nb304nf_single_loop
	jmp   .nb304nf_updateouterdata
.nb304nf_single_loop:
	mov   rdx, [rsp + nb304nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb304nf_innerjjnr],  4	

	mov rsi, [rbp + nb304nf_pos]
	lea   rax, [rax + rax*2]  

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
	movaps  xmm0, [rsp + nb304nf_ixM]     
	movaps  xmm1, [rsp + nb304nf_iyM]
	movaps  xmm2, [rsp + nb304nf_izM]	
	movlhps xmm3, xmm6			;# xmm3 = jxM   0   jxH1 jxH2 
	shufps  xmm4, xmm6, 228 ;# 11100100	;# xmm4 = jyM   0   jyH1 jyH2 
	shufps  xmm5, xmm7, 68  ;# 01000100	;# xmm5 = jzM   0   jzH1 jzH2
	
	;# store all j coordinates in jM 
	movaps [rsp + nb304nf_jxM], xmm3
	movaps [rsp + nb304nf_jyM], xmm4
	movaps [rsp + nb304nf_jzM], xmm5
	subps  xmm0, xmm3
	subps  xmm1, xmm4
	subps  xmm2, xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	movaps  xmm2, xmm1	
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb304nf_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [rsp + nb304nf_half] ;# rinv iO - j water 

	movaps  xmm1, xmm3
	mulps   xmm1, xmm0	;# xmm1=r 
	movaps  xmm0, xmm3	;# xmm0=rinv 
	mulps  xmm1, [rsp + nb304nf_tsc]
	
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
	
	movd ebx, mm6
	movd ecx, mm7
	psrlq mm7, 32
	movd edx, mm7		;# table indices in ebx,ecx,edx 

	mov rsi, [rbp + nb304nf_VFtab]
	
	movlps xmm5, [rsi + rbx*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rbx*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  
	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 

	xorps  xmm3, xmm3

	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [rsp + nb304nf_qqMM]
	movhps  xmm3, [rsp + nb304nf_qqMH]
	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point xmm5 contains vcoul 
	
	addps  xmm5, [rsp + nb304nf_vctot]
	movaps [rsp + nb304nf_vctot], xmm5
	
	;# done with i M Now do i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [rsp + nb304nf_ixH1]
	movaps  xmm1, [rsp + nb304nf_iyH1]
	movaps  xmm2, [rsp + nb304nf_izH1]	
	movaps  xmm3, [rsp + nb304nf_ixH2] 
	movaps  xmm4, [rsp + nb304nf_iyH2] 
	movaps  xmm5, [rsp + nb304nf_izH2] 
	subps   xmm0, [rsp + nb304nf_jxM]
	subps   xmm1, [rsp + nb304nf_jyM]
	subps   xmm2, [rsp + nb304nf_jzM]
	subps   xmm3, [rsp + nb304nf_jxM]
	subps   xmm4, [rsp + nb304nf_jyM]
	subps   xmm5, [rsp + nb304nf_jzM]
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

	;# start with H1, save H2 data 
	movaps [rsp + nb304nf_rsqH2M], xmm4
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb304nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb304nf_half] ;# rinv H1 - j water 
	mulps   xmm7, [rsp + nb304nf_half] ;# rinv H2 - j water  

	;# start with H1, save H2 data 
	movaps [rsp + nb304nf_rinvH2M], xmm7

	movaps xmm1, xmm3
	mulps  xmm1, xmm0	;# xmm1=r 
	movaps xmm0, xmm3	;# xmm0=rinv 
	mulps  xmm1, [rsp + nb304nf_tsc]
	
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

	movd ebx, mm6
	movd ecx, mm7
	psrlq mm7, 32
	movd edx, mm7		;# table indices in ebx,ecx,edx 

	movlps xmm5, [rsi + rbx*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rbx*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  
	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 

	xorps  xmm3, xmm3
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [rsp + nb304nf_qqMH]
	movhps  xmm3, [rsp + nb304nf_qqHH]
	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point xmm5 contains vcoul 
	addps  xmm5, [rsp + nb304nf_vctot]
	movaps [rsp + nb304nf_vctot], xmm5	

	;# do table for H2 - j water interaction 
	movaps xmm0, [rsp + nb304nf_rinvH2M]
	movaps xmm1, [rsp + nb304nf_rsqH2M]
	mulps  xmm1, xmm0	;# xmm0=rinv, xmm1=r 
	mulps  xmm1, [rsp + nb304nf_tsc]
	
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

	movd ebx, mm6
	movd ecx, mm7
	psrlq mm7, 32
	movd edx, mm7		;# table indices in ebx,ecx,edx 

	movlps xmm5, [rsi + rbx*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rbx*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  
	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 

	xorps  xmm3, xmm3
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [rsp + nb304nf_qqMH]
	movhps  xmm3, [rsp + nb304nf_qqHH]
	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point xmm5 contains vcoul 
	addps  xmm5, [rsp + nb304nf_vctot]
	movaps [rsp + nb304nf_vctot], xmm5	
	
	dec dword ptr [rsp + nb304nf_innerk]
	jz    .nb304nf_updateouterdata
	jmp   .nb304nf_single_loop
.nb304nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb304nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb304nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb304nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb304nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb304nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb304nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb304nf_n], esi
        jmp .nb304nf_outer
.nb304nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb304nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb304nf_end
        ;# non-zero, do one more workunit
        jmp   .nb304nf_threadloop
.nb304nf_end:
	mov eax, [rsp + nb304nf_nouter]
	mov ebx, [rsp + nb304nf_ninner]
	mov rcx, [rbp + nb304nf_outeriter]
	mov rdx, [rbp + nb304nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 792
	femms

	pop rbx
	pop	rbp
	ret



	
