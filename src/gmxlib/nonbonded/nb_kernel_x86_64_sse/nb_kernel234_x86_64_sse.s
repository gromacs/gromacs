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


	
.globl nb_kernel234_x86_64_sse
.globl _nb_kernel234_x86_64_sse
nb_kernel234_x86_64_sse:	
_nb_kernel234_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb234_fshift,           16
.equiv          nb234_gid,              24
.equiv          nb234_pos,              32
.equiv          nb234_faction,          40
.equiv          nb234_charge,           48
.equiv          nb234_p_facel,          56
.equiv          nb234_argkrf,           64
.equiv          nb234_argcrf,           72
.equiv          nb234_Vc,               80
.equiv          nb234_type,             88
.equiv          nb234_p_ntype,          96
.equiv          nb234_vdwparam,         104
.equiv          nb234_Vvdw,             112
.equiv          nb234_p_tabscale,       120
.equiv          nb234_VFtab,            128
.equiv          nb234_invsqrta,         136
.equiv          nb234_dvda,             144
.equiv          nb234_p_gbtabscale,     152
.equiv          nb234_GBtab,            160
.equiv          nb234_p_nthreads,       168
.equiv          nb234_count,            176
.equiv          nb234_mtx,              184
.equiv          nb234_outeriter,        192
.equiv          nb234_inneriter,        200
.equiv          nb234_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb234_ixO,              0
.equiv          nb234_iyO,              16
.equiv          nb234_izO,              32
.equiv          nb234_ixH1,             48
.equiv          nb234_iyH1,             64
.equiv          nb234_izH1,             80
.equiv          nb234_ixH2,             96
.equiv          nb234_iyH2,             112
.equiv          nb234_izH2,             128
.equiv          nb234_ixM,              144
.equiv          nb234_iyM,              160
.equiv          nb234_izM,              176
.equiv          nb234_jxO,              192
.equiv          nb234_jyO,              208
.equiv          nb234_jzO,              224
.equiv          nb234_jxH1,             240
.equiv          nb234_jyH1,             256
.equiv          nb234_jzH1,             272
.equiv          nb234_jxH2,             288
.equiv          nb234_jyH2,             304
.equiv          nb234_jzH2,             320
.equiv          nb234_jxM,              336
.equiv          nb234_jyM,              352
.equiv          nb234_jzM,              368
.equiv          nb234_dxOO,             384
.equiv          nb234_dyOO,             400
.equiv          nb234_dzOO,             416
.equiv          nb234_dxH1H1,           432
.equiv          nb234_dyH1H1,           448
.equiv          nb234_dzH1H1,           464
.equiv          nb234_dxH1H2,           480
.equiv          nb234_dyH1H2,           496
.equiv          nb234_dzH1H2,           512
.equiv          nb234_dxH1M,            528
.equiv          nb234_dyH1M,            544
.equiv          nb234_dzH1M,            560
.equiv          nb234_dxH2H1,           576
.equiv          nb234_dyH2H1,           592
.equiv          nb234_dzH2H1,           608
.equiv          nb234_dxH2H2,           624
.equiv          nb234_dyH2H2,           640
.equiv          nb234_dzH2H2,           656
.equiv          nb234_dxH2M,            672
.equiv          nb234_dyH2M,            688
.equiv          nb234_dzH2M,            704
.equiv          nb234_dxMH1,            720
.equiv          nb234_dyMH1,            736
.equiv          nb234_dzMH1,            752
.equiv          nb234_dxMH2,            768
.equiv          nb234_dyMH2,            784
.equiv          nb234_dzMH2,            800
.equiv          nb234_dxMM,             816
.equiv          nb234_dyMM,             832
.equiv          nb234_dzMM,             848
.equiv          nb234_qqMM,             864
.equiv          nb234_qqMH,             880
.equiv          nb234_qqHH,             896
.equiv          nb234_two,              912
.equiv          nb234_c6,               928
.equiv          nb234_c12,              944
.equiv          nb234_tsc,              960
.equiv          nb234_fstmp ,           976
.equiv          nb234_vctot,            992
.equiv          nb234_Vvdwtot,          1008
.equiv          nb234_fixO,             1024
.equiv          nb234_fiyO,             1040
.equiv          nb234_fizO,             1056
.equiv          nb234_fixH1,            1072
.equiv          nb234_fiyH1,            1088
.equiv          nb234_fizH1,            1104
.equiv          nb234_fixH2,            1120
.equiv          nb234_fiyH2,            1136
.equiv          nb234_fizH2,            1152
.equiv          nb234_fixM,             1168
.equiv          nb234_fiyM,             1184
.equiv          nb234_fizM,             1200
.equiv          nb234_fjxO,             1216
.equiv          nb234_fjyO,             1232
.equiv          nb234_fjzO,             1248
.equiv          nb234_fjxH1,            1264
.equiv          nb234_fjyH1,            1280
.equiv          nb234_fjzH1,            1296
.equiv          nb234_fjxH2,            1312
.equiv          nb234_fjyH2,            1328
.equiv          nb234_fjzH2,            1344
.equiv          nb234_fjxM,             1360
.equiv          nb234_fjyM,             1376
.equiv          nb234_fjzM,             1392
.equiv          nb234_half,             1408
.equiv          nb234_three,            1424
.equiv          nb234_rsqOO,            1440
.equiv          nb234_rsqH1H1,          1456
.equiv          nb234_rsqH1H2,          1472
.equiv          nb234_rsqH1M,           1488
.equiv          nb234_rsqH2H1,          1504
.equiv          nb234_rsqH2H2,          1520
.equiv          nb234_rsqH2M,           1536
.equiv          nb234_rsqMH1,           1552
.equiv          nb234_rsqMH2,           1568
.equiv          nb234_rsqMM,            1584
.equiv          nb234_rinvOO,           1600
.equiv          nb234_rinvH1H1,         1616
.equiv          nb234_rinvH1H2,         1632
.equiv          nb234_rinvH1M,          1648
.equiv          nb234_rinvH2H1,         1664
.equiv          nb234_rinvH2H2,         1680
.equiv          nb234_rinvH2M,          1696
.equiv          nb234_rinvMH1,          1712
.equiv          nb234_rinvMH2,          1728
.equiv          nb234_rinvMM,           1744
.equiv          nb234_krf,              1776
.equiv          nb234_crf,              1792
.equiv          nb234_is3,              1808
.equiv          nb234_ii3,              1812
.equiv          nb234_innerjjnr,        1816
.equiv          nb234_nri,              1824
.equiv          nb234_iinr,             1832
.equiv          nb234_jindex,           1840
.equiv          nb234_jjnr,             1848
.equiv          nb234_shift,            1856
.equiv          nb234_shiftvec,         1864
.equiv          nb234_facel,            1872
.equiv          nb234_innerk,           1880
.equiv          nb234_n,                1884
.equiv          nb234_nn1,              1888
.equiv          nb234_nouter,           1892
.equiv          nb234_ninner,           1896

	push rbp
	mov  rbp, rsp
	push rbx

	
	femms
	sub rsp, 1912		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb234_nouter], eax
	mov [rsp + nb234_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb234_nri], edi
	mov [rsp + nb234_iinr], rsi
	mov [rsp + nb234_jindex], rdx
	mov [rsp + nb234_jjnr], rcx
	mov [rsp + nb234_shift], r8
	mov [rsp + nb234_shiftvec], r9
	mov rsi, [rbp + nb234_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb234_facel], xmm0

	mov rax, [rbp + nb234_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb234_tsc], xmm3

	mov rsi, [rbp + nb234_argkrf]
	mov rdi, [rbp + nb234_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb234_krf], xmm1
	movaps [rsp + nb234_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb234_half], eax
	movss xmm1, [rsp + nb234_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb234_half],  xmm1
	movaps [rsp + nb234_two],  xmm2
	movaps [rsp + nb234_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb234_iinr]   ;# rcx = pointer into iinr[]
	mov   ebx, [rcx]		;# ebx =ii 

	mov   rdx, [rbp + nb234_charge]
	movss xmm5, [rdx + rbx*4 + 4]	
	movss xmm3, [rdx + rbx*4 + 12]	
	movss xmm4, xmm3	
	mov rsi, [rbp + nb234_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb234_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb234_qqMM], xmm3
	movaps [rsp + nb234_qqMH], xmm4
	movaps [rsp + nb234_qqHH], xmm5
	
	xorps xmm0, xmm0
	mov   rdx, [rbp + nb234_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb234_p_ntype]
	imul  ecx, [rdi]  ;# rcx = ntia = 2*ntype*type[ii0] 
	add   rdx, rcx
	mov   rax, [rbp + nb234_vdwparam]
	movlps xmm0, [rax + rdx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0x55
	movaps [rsp + nb234_c6], xmm0
	movaps [rsp + nb234_c12], xmm1

.nb234_threadloop:
        mov   rsi, [rbp + nb234_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb234_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [rsi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb234_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb234_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb234_n], eax
        mov [rsp + nb234_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb234_outerstart
        jmp .nb234_end
	
.nb234_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb234_nouter]
	mov [rsp + nb234_nouter], ebx

.nb234_outer:
	mov   rax, [rsp + nb234_shift]  	;# eax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb234_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb234_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb234_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb234_pos]	;# eax = base of pos[]  
	mov   [rsp + nb234_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	movaps xmm6, xmm0
	movaps xmm7, xmm1

	addss xmm3, [rax + rbx*4]  	;# ox
	addss xmm4, [rax + rbx*4 + 4]  ;# oy
	addss xmm5, [rax + rbx*4 + 8]  ;# oz
	addss xmm6, [rax + rbx*4 + 12] ;# h1x
	addss xmm7, [rax + rbx*4 + 16] ;# h1y
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	shufps xmm7, xmm7, 0
	movaps [rsp + nb234_ixO], xmm3
	movaps [rsp + nb234_iyO], xmm4
	movaps [rsp + nb234_izO], xmm5
	movaps [rsp + nb234_ixH1], xmm6
	movaps [rsp + nb234_iyH1], xmm7

	movss xmm6, xmm2
	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm6, [rax + rbx*4 + 20] ;# h1z
	addss xmm0, [rax + rbx*4 + 24] ;# h2x
	addss xmm1, [rax + rbx*4 + 28] ;# h2y
	addss xmm2, [rax + rbx*4 + 32] ;# h2z
	addss xmm3, [rax + rbx*4 + 36] ;# mx
	addss xmm4, [rax + rbx*4 + 40] ;# my
	addss xmm5, [rax + rbx*4 + 44] ;# mz

	shufps xmm6, xmm6, 0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb234_izH1], xmm6
	movaps [rsp + nb234_ixH2], xmm0
	movaps [rsp + nb234_iyH2], xmm1
	movaps [rsp + nb234_izH2], xmm2
	movaps [rsp + nb234_ixM], xmm3
	movaps [rsp + nb234_iyM], xmm4
	movaps [rsp + nb234_izM], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb234_vctot], xmm4
	movaps [rsp + nb234_Vvdwtot], xmm4
	movaps [rsp + nb234_fixO], xmm4
	movaps [rsp + nb234_fiyO], xmm4
	movaps [rsp + nb234_fizO], xmm4
	movaps [rsp + nb234_fixH1], xmm4
	movaps [rsp + nb234_fiyH1], xmm4
	movaps [rsp + nb234_fizH1], xmm4
	movaps [rsp + nb234_fixH2], xmm4
	movaps [rsp + nb234_fiyH2], xmm4
	movaps [rsp + nb234_fizH2], xmm4
	movaps [rsp + nb234_fixM], xmm4
	movaps [rsp + nb234_fiyM], xmm4
	movaps [rsp + nb234_fizM], xmm4
	
	mov   rax, [rsp + nb234_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb234_pos]
	mov   rdi, [rbp + nb234_faction]	
	mov   rax, [rsp + nb234_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb234_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb234_ninner]
	mov   [rsp + nb234_ninner], ecx
	add   edx, 0
	mov   [rsp + nb234_innerk], edx	;# number of innerloop atoms 
	jge   .nb234_unroll_loop
	jmp   .nb234_single_check
.nb234_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb234_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb234_innerjjnr],  16 ;# advance pointer (unroll 4) 

	mov rsi, [rbp + nb234_pos]   	;# base of pos[] 

	lea   rax, [rax + rax*2] 	;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2] 	;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	
	
	;# move j coordinates to local temp variables
	;# Load Ox, Oy, Oz, H1x 
	movlps xmm1, [rsi + rax*4]	;#  Oxa   Oya    -    -
	movlps xmm4, [rsi + rcx*4]	;#  Oxc   Oyc    -    -
	movhps xmm1, [rsi + rbx*4]	;#  Oxa   Oya   Oxb   Oyb 
	movhps xmm4, [rsi + rdx*4]	;#  Oxc   Oyc   Oxd   Oyd 
	movaps xmm0, xmm1		;#  Oxa   Oya   Oxb   Oyb 
	shufps xmm0, xmm4, 0x88		;#  Oxa   Oxb   Oxc   Oxd
	shufps xmm1, xmm4, 0xDD		;#  Oya   Oyb   Oyc   Oyd
	movlps xmm3, [rsi + rax*4 + 8]	;#  Oza  H1xa    -    -
	movlps xmm5, [rsi + rcx*4 + 8]	;#  Ozc  H1xc    -    -
	movhps xmm3, [rsi + rbx*4 + 8]	;#  Oza  H1xa   Ozb  H1xb 
	movhps xmm5, [rsi + rdx*4 + 8]	;#  Ozc  H1xc   Ozd  H1xd 
	movaps xmm2, xmm3		;#  Oza  H1xa   Ozb  H1xb 
	shufps xmm2, xmm5, 0x88		;#  Oza   Ozb   Ozc   Ozd
	shufps xmm3, xmm5, 0xDD		;# H1xa  H1xb  H1xc  H1xd
	;# coordinates in xmm0-xmm3	
	;# store
	movaps [rsp + nb234_jxO], xmm0
	movaps [rsp + nb234_jyO], xmm1
	movaps [rsp + nb234_jzO], xmm2
	movaps [rsp + nb234_jxH1], xmm3

	;# Load H1y H1z H2x H2y 
	movlps xmm1, [rsi + rax*4 + 16]	
	movlps xmm4, [rsi + rcx*4 + 16]	
	movhps xmm1, [rsi + rbx*4 + 16]	
	movhps xmm4, [rsi + rdx*4 + 16]	
	movaps xmm0, xmm1		
	shufps xmm0, xmm4, 0x88		
	shufps xmm1, xmm4, 0xDD		
	movlps xmm3, [rsi + rax*4 + 24]	
	movlps xmm5, [rsi + rcx*4 + 24]	
	movhps xmm3, [rsi + rbx*4 + 24]	
	movhps xmm5, [rsi + rdx*4 + 24]	
	movaps xmm2, xmm3		
	shufps xmm2, xmm5, 0x88		
	shufps xmm3, xmm5, 0xDD		
	;# coordinates in xmm0-xmm3	
	;# store
	movaps [rsp + nb234_jyH1], xmm0
	movaps [rsp + nb234_jzH1], xmm1
	movaps [rsp + nb234_jxH2], xmm2
	movaps [rsp + nb234_jyH2], xmm3

	;# Load H2z Mx My Mz 
	movlps xmm1, [rsi + rax*4 + 32]	
	movlps xmm4, [rsi + rcx*4 + 32]	
	movhps xmm1, [rsi + rbx*4 + 32]	
	movhps xmm4, [rsi + rdx*4 + 32]	
	movaps xmm0, xmm1		
	shufps xmm0, xmm4, 0x88		
	shufps xmm1, xmm4, 0xDD		
	movlps xmm3, [rsi + rax*4 + 40]	
	movlps xmm5, [rsi + rcx*4 + 40]	
	movhps xmm3, [rsi + rbx*4 + 40]	
	movhps xmm5, [rsi + rdx*4 + 40]	
	movaps xmm2, xmm3		
	shufps xmm2, xmm5, 0x88		
	shufps xmm3, xmm5, 0xDD		
	;# coordinates in xmm0-xmm3	
	;# store
	movaps [rsp + nb234_jzH2], xmm0
	movaps [rsp + nb234_jxM], xmm1
	movaps [rsp + nb234_jyM], xmm2
	movaps [rsp + nb234_jzM], xmm3
	
	;# start calculating pairwise distances
	movaps xmm0, [rsp + nb234_ixO]
	movaps xmm1, [rsp + nb234_iyO]
	movaps xmm2, [rsp + nb234_izO]
	movaps xmm3, [rsp + nb234_ixH1]
	movaps xmm4, [rsp + nb234_iyH1]
	movaps xmm5, [rsp + nb234_izH1]
	subps  xmm0, [rsp + nb234_jxO]
	subps  xmm1, [rsp + nb234_jyO]
	subps  xmm2, [rsp + nb234_jzO]
	subps  xmm3, [rsp + nb234_jxH1]
	subps  xmm4, [rsp + nb234_jyH1]
	subps  xmm5, [rsp + nb234_jzH1]
	movaps [rsp + nb234_dxOO], xmm0
	movaps [rsp + nb234_dyOO], xmm1
	movaps [rsp + nb234_dzOO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb234_dxH1H1], xmm3
	movaps [rsp + nb234_dyH1H1], xmm4
	movaps [rsp + nb234_dzH1H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + nb234_rsqOO], xmm0
	movaps [rsp + nb234_rsqH1H1], xmm3

	movaps xmm0, [rsp + nb234_ixH1]
	movaps xmm1, [rsp + nb234_iyH1]
	movaps xmm2, [rsp + nb234_izH1]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [rsp + nb234_jxH2]
	subps  xmm1, [rsp + nb234_jyH2]
	subps  xmm2, [rsp + nb234_jzH2]
	subps  xmm3, [rsp + nb234_jxM]
	subps  xmm4, [rsp + nb234_jyM]
	subps  xmm5, [rsp + nb234_jzM]
	movaps [rsp + nb234_dxH1H2], xmm0
	movaps [rsp + nb234_dyH1H2], xmm1
	movaps [rsp + nb234_dzH1H2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb234_dxH1M], xmm3
	movaps [rsp + nb234_dyH1M], xmm4
	movaps [rsp + nb234_dzH1M], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + nb234_rsqH1H2], xmm0
	movaps [rsp + nb234_rsqH1M], xmm3

	movaps xmm0, [rsp + nb234_ixH2]
	movaps xmm1, [rsp + nb234_iyH2]
	movaps xmm2, [rsp + nb234_izH2]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [rsp + nb234_jxH1]
	subps  xmm1, [rsp + nb234_jyH1]
	subps  xmm2, [rsp + nb234_jzH1]
	subps  xmm3, [rsp + nb234_jxH2]
	subps  xmm4, [rsp + nb234_jyH2]
	subps  xmm5, [rsp + nb234_jzH2]
	movaps [rsp + nb234_dxH2H1], xmm0
	movaps [rsp + nb234_dyH2H1], xmm1
	movaps [rsp + nb234_dzH2H1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb234_dxH2H2], xmm3
	movaps [rsp + nb234_dyH2H2], xmm4
	movaps [rsp + nb234_dzH2H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + nb234_rsqH2H1], xmm0
	movaps [rsp + nb234_rsqH2H2], xmm3

	movaps xmm0, [rsp + nb234_ixH2]
	movaps xmm1, [rsp + nb234_iyH2]
	movaps xmm2, [rsp + nb234_izH2]
	movaps xmm3, [rsp + nb234_ixM]
	movaps xmm4, [rsp + nb234_iyM]
	movaps xmm5, [rsp + nb234_izM]
	subps  xmm0, [rsp + nb234_jxM]
	subps  xmm1, [rsp + nb234_jyM]
	subps  xmm2, [rsp + nb234_jzM]
	subps  xmm3, [rsp + nb234_jxH1]
	subps  xmm4, [rsp + nb234_jyH1]
	subps  xmm5, [rsp + nb234_jzH1]
	movaps [rsp + nb234_dxH2M], xmm0
	movaps [rsp + nb234_dyH2M], xmm1
	movaps [rsp + nb234_dzH2M], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb234_dxMH1], xmm3
	movaps [rsp + nb234_dyMH1], xmm4
	movaps [rsp + nb234_dzMH1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [rsp + nb234_rsqH2M], xmm0
	movaps [rsp + nb234_rsqMH1], xmm4

	movaps xmm0, [rsp + nb234_ixM]
	movaps xmm1, [rsp + nb234_iyM]
	movaps xmm2, [rsp + nb234_izM]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [rsp + nb234_jxH2]
	subps  xmm1, [rsp + nb234_jyH2]
	subps  xmm2, [rsp + nb234_jzH2]
	subps  xmm3, [rsp + nb234_jxM]
	subps  xmm4, [rsp + nb234_jyM]
	subps  xmm5, [rsp + nb234_jzM]
	movaps [rsp + nb234_dxMH2], xmm0
	movaps [rsp + nb234_dyMH2], xmm1
	movaps [rsp + nb234_dzMH2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb234_dxMM], xmm3
	movaps [rsp + nb234_dyMM], xmm4
	movaps [rsp + nb234_dzMM], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [rsp + nb234_rsqMH2], xmm0
	movaps [rsp + nb234_rsqMM], xmm4

	;# start by doing invsqrt for OO
	rsqrtps xmm1, [rsp + nb234_rsqOO]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb234_three]
	mulps   xmm1, [rsp + nb234_rsqOO]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + nb234_half] 
	movaps  [rsp + nb234_rinvOO], xmm3
	
	;# more invsqrt ops - do two at a time.
	rsqrtps xmm1, [rsp + nb234_rsqH1H1]
	rsqrtps xmm5, [rsp + nb234_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb234_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb234_rsqH1H1]
	mulps   xmm5, [rsp + nb234_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb234_half] ;# rinvH1H1 
	mulps   xmm7, [rsp + nb234_half] ;# rinvH1H2 
	movaps  [rsp + nb234_rinvH1H1], xmm3
	movaps  [rsp + nb234_rinvH1H2], xmm7
			
	rsqrtps xmm1, [rsp + nb234_rsqH1M]
	rsqrtps xmm5, [rsp + nb234_rsqH2H1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb234_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb234_rsqH1M]
	mulps   xmm5, [rsp + nb234_rsqH2H1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb234_half] 
	mulps   xmm7, [rsp + nb234_half]
	movaps  [rsp + nb234_rinvH1M], xmm3
	movaps  [rsp + nb234_rinvH2H1], xmm7
			
	rsqrtps xmm1, [rsp + nb234_rsqH2H2]
	rsqrtps xmm5, [rsp + nb234_rsqH2M]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb234_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb234_rsqH2H2]
	mulps   xmm5, [rsp + nb234_rsqH2M]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb234_half] 
	mulps   xmm7, [rsp + nb234_half]
	movaps  [rsp + nb234_rinvH2H2], xmm3
	movaps  [rsp + nb234_rinvH2M], xmm7
	
	rsqrtps xmm1, [rsp + nb234_rsqMH1]
	rsqrtps xmm5, [rsp + nb234_rsqMH2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb234_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb234_rsqMH1]
	mulps   xmm5, [rsp + nb234_rsqMH2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb234_half] 
	mulps   xmm7, [rsp + nb234_half]
	movaps  [rsp + nb234_rinvMH1], xmm3
	movaps  [rsp + nb234_rinvMH2], xmm7
        		
	rsqrtps xmm1, [rsp + nb234_rsqMM]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb234_three]
	mulps   xmm1, [rsp + nb234_rsqMM]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + nb234_half] 
	movaps  [rsp + nb234_rinvMM], xmm3
	
	;# start with OO LJ interaction
	movaps xmm0, [rsp + nb234_rinvOO]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb234_rsqOO] ;# xmm1=r 
	mulps  xmm1, [rsp + nb234_tsc]
		
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

    mov  rsi, [rbp + nb234_VFtab]
	
    ;# dispersion 
    movlps xmm5, [rsi + rax*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm5, [rsi + rbx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# constant 10001000
    shufps xmm5, xmm7, 221  ;# constant 11011101

    movlps xmm7, [rsi + rax*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm7, [rsi + rbx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# constant 10001000
    shufps xmm7, xmm3, 221  ;# constant 11011101
    ;# dispersion table ready, in xmm4-xmm7 
    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    mulps  xmm7, [rsp + nb234_two]       ;# two*Heps2 
    addps  xmm7, xmm6
    addps  xmm7, xmm5 ;# xmm7=FF 
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 

    movaps xmm4, [rsp + nb234_c6]
    mulps  xmm7, xmm4    ;# fijD 
    mulps  xmm5, xmm4    ;# Vvdw6 
    movaps  [rsp + nb234_fstmp], xmm7

    addps  xmm5, [rsp + nb234_Vvdwtot]
    movaps [rsp + nb234_Vvdwtot], xmm5

    ;# repulsion 
    movlps xmm5, [rsi + rax*4 + 16]
    movlps xmm7, [rsi + rcx*4 + 16]
    movhps xmm5, [rsi + rbx*4 + 16]
    movhps xmm7, [rsi + rdx*4 + 16] ;# got half table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# constant 10001000
    shufps xmm5, xmm7, 221  ;# constant 11011101

    movlps xmm7, [rsi + rax*4 + 24]
    movlps xmm3, [rsi + rcx*4 + 24]
    movhps xmm7, [rsi + rbx*4 + 24]
    movhps xmm3, [rsi + rdx*4 + 24] ;# other half of table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# constant 10001000
    shufps xmm7, xmm3, 221  ;# constant 11011101
    ;# repulsion table ready, in xmm4-xmm7 
    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    mulps  xmm7, [rsp + nb234_two]       ;# two*Heps2 
    addps  xmm7, xmm6
    addps  xmm7, xmm5 ;# xmm7=FF 
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
 
    movaps xmm4, [rsp + nb234_c12]
    mulps  xmm7, xmm4 ;# fijR 
    mulps  xmm5, xmm4 ;# Vvdw12 
    addps  xmm7, [rsp + nb234_fstmp] 
    mulps xmm7, [rsp + nb234_tsc]

    addps  xmm5, [rsp + nb234_Vvdwtot]
    movaps [rsp + nb234_Vvdwtot], xmm5

	xorps xmm2, xmm2
	mulps  xmm7, xmm0
	subps xmm2, xmm7

	movaps xmm0, xmm2
	movaps xmm1, xmm2
	
    movd  eax, mm0
    movd  ebx, mm1
    movd  ecx, mm2
    movd  edx, mm3

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [rsp + nb234_dxOO]
	mulps xmm1, [rsp + nb234_dyOO]
	mulps xmm2, [rsp + nb234_dzOO]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb234_fixO]
	addps xmm1, [rsp + nb234_fiyO]
	addps xmm2, [rsp + nb234_fizO]
	movaps [rsp + nb234_fjxO], xmm3
	movaps [rsp + nb234_fjyO], xmm4
	movaps [rsp + nb234_fjzO], xmm5
	movaps [rsp + nb234_fixO], xmm0
	movaps [rsp + nb234_fiyO], xmm1
	movaps [rsp + nb234_fizO], xmm2

	;# Coulomb interactions 
	;# start with H1-H1 interaction 
	movaps xmm0, [rsp + nb234_rinvH1H1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb234_krf]
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	mulps  xmm5, [rsp + nb234_rsqH1H1] ;# xmm5=krsq 
	movaps xmm6, xmm5
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [rsp + nb234_crf]
	mulps  xmm6, [rsp + nb234_qqHH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	mulps xmm5, [rsp + nb234_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [rsp + nb234_qqHH] ;# xmm7 = coul part of fscal 
	
	addps  xmm6, [rsp + nb234_vctot] ;# local vctot summation variable 
	mulps  xmm0, xmm7
	
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [rsp + nb234_dxH1H1]
	mulps xmm1, [rsp + nb234_dyH1H1]
	mulps xmm2, [rsp + nb234_dzH1H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb234_fixH1]
	addps xmm1, [rsp + nb234_fiyH1]
	addps xmm2, [rsp + nb234_fizH1]
	movaps [rsp + nb234_fjxH1], xmm3
	movaps [rsp + nb234_fjyH1], xmm4
	movaps [rsp + nb234_fjzH1], xmm5
	movaps [rsp + nb234_fixH1], xmm0
	movaps [rsp + nb234_fiyH1], xmm1
	movaps [rsp + nb234_fizH1], xmm2

	;# H1-H2 interaction 
	movaps xmm0, [rsp + nb234_rinvH1H2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb234_krf]
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb234_rsqH1H2] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb234_crf]
	mulps  xmm0, xmm0
	mulps  xmm4, [rsp + nb234_qqHH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [rsp + nb234_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [rsp + nb234_qqHH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH1  
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [rsp + nb234_dxH1H2]
	mulps xmm1, [rsp + nb234_dyH1H2]
	mulps xmm2, [rsp + nb234_dzH1H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb234_fixH1]
	addps xmm1, [rsp + nb234_fiyH1]
	addps xmm2, [rsp + nb234_fizH1]
	movaps [rsp + nb234_fjxH2], xmm3
	movaps [rsp + nb234_fjyH2], xmm4
	movaps [rsp + nb234_fjzH2], xmm5
	movaps [rsp + nb234_fixH1], xmm0
	movaps [rsp + nb234_fiyH1], xmm1
	movaps [rsp + nb234_fizH1], xmm2

	;# H1-M interaction  
	movaps xmm0, [rsp + nb234_rinvH1M]
	movaps xmm7, xmm0	;# xmm7=Rinv 
	movaps xmm5, [rsp + nb234_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb234_rsqH1M] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb234_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb234_qqMH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [rsp + nb234_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [rsp + nb234_qqMH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [rsp + nb234_dxH1M]
	mulps xmm1, [rsp + nb234_dyH1M]
	mulps xmm2, [rsp + nb234_dzH1M]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb234_fixH1]
	addps xmm1, [rsp + nb234_fiyH1]
	addps xmm2, [rsp + nb234_fizH1]
	movaps [rsp + nb234_fjxM], xmm3
	movaps [rsp + nb234_fjyM], xmm4
	movaps [rsp + nb234_fjzM], xmm5
	movaps [rsp + nb234_fixH1], xmm0
	movaps [rsp + nb234_fiyH1], xmm1
	movaps [rsp + nb234_fizH1], xmm2

	;# H2-H1 interaction 
	movaps xmm0, [rsp + nb234_rinvH2H1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb234_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb234_rsqH2H1] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb234_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb234_qqHH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [rsp + nb234_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [rsp + nb234_qqHH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [rsp + nb234_fjxH1]
	movaps xmm4, [rsp + nb234_fjyH1]
	movaps xmm5, [rsp + nb234_fjzH1]
	mulps xmm0, [rsp + nb234_dxH2H1]
	mulps xmm1, [rsp + nb234_dyH2H1]
	mulps xmm2, [rsp + nb234_dzH2H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb234_fixH2]
	addps xmm1, [rsp + nb234_fiyH2]
	addps xmm2, [rsp + nb234_fizH2]
	movaps [rsp + nb234_fjxH1], xmm3
	movaps [rsp + nb234_fjyH1], xmm4
	movaps [rsp + nb234_fjzH1], xmm5
	movaps [rsp + nb234_fixH2], xmm0
	movaps [rsp + nb234_fiyH2], xmm1
	movaps [rsp + nb234_fizH2], xmm2

	;# H2-H2 interaction 
	movaps xmm0, [rsp + nb234_rinvH2H2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb234_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb234_rsqH2H2] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb234_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb234_qqHH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [rsp + nb234_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [rsp + nb234_qqHH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [rsp + nb234_fjxH2]
	movaps xmm4, [rsp + nb234_fjyH2]
	movaps xmm5, [rsp + nb234_fjzH2]
	mulps xmm0, [rsp + nb234_dxH2H2]
	mulps xmm1, [rsp + nb234_dyH2H2]
	mulps xmm2, [rsp + nb234_dzH2H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb234_fixH2]
	addps xmm1, [rsp + nb234_fiyH2]
	addps xmm2, [rsp + nb234_fizH2]
	movaps [rsp + nb234_fjxH2], xmm3
	movaps [rsp + nb234_fjyH2], xmm4
	movaps [rsp + nb234_fjzH2], xmm5
	movaps [rsp + nb234_fixH2], xmm0
	movaps [rsp + nb234_fiyH2], xmm1
	movaps [rsp + nb234_fizH2], xmm2

	;# H2-M interaction 
	movaps xmm0, [rsp + nb234_rinvH2M]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb234_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb234_rsqH2M] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb234_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb234_qqMH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [rsp + nb234_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [rsp + nb234_qqMH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	movaps xmm3, [rsp + nb234_fjxM]
	movaps xmm4, [rsp + nb234_fjyM]
	movaps xmm5, [rsp + nb234_fjzM]
	mulps xmm0, [rsp + nb234_dxH2M]
	mulps xmm1, [rsp + nb234_dyH2M]
	mulps xmm2, [rsp + nb234_dzH2M]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb234_fixH2]
	addps xmm1, [rsp + nb234_fiyH2]
	addps xmm2, [rsp + nb234_fizH2]
	movaps [rsp + nb234_fjxM], xmm3
	movaps [rsp + nb234_fjyM], xmm4
	movaps [rsp + nb234_fjzM], xmm5
	movaps [rsp + nb234_fixH2], xmm0
	movaps [rsp + nb234_fiyH2], xmm1
	movaps [rsp + nb234_fizH2], xmm2

	;# M-H1 interaction 
	movaps xmm0, [rsp + nb234_rinvMH1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb234_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb234_rsqMH1] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb234_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb234_qqMH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [rsp + nb234_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [rsp + nb234_qqMH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [rsp + nb234_fjxH1]
	movaps xmm4, [rsp + nb234_fjyH1]
	movaps xmm5, [rsp + nb234_fjzH1]
	mulps xmm0, [rsp + nb234_dxMH1]
	mulps xmm1, [rsp + nb234_dyMH1]
	mulps xmm2, [rsp + nb234_dzMH1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb234_fixM]
	addps xmm1, [rsp + nb234_fiyM]
	addps xmm2, [rsp + nb234_fizM]
	movaps [rsp + nb234_fjxH1], xmm3
	movaps [rsp + nb234_fjyH1], xmm4
	movaps [rsp + nb234_fjzH1], xmm5
	movaps [rsp + nb234_fixM], xmm0
	movaps [rsp + nb234_fiyM], xmm1
	movaps [rsp + nb234_fizM], xmm2

	;# M-H2 interaction 
	movaps xmm0, [rsp + nb234_rinvMH2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb234_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb234_rsqMH2] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb234_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb234_qqMH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [rsp + nb234_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [rsp + nb234_qqMH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [rsp + nb234_fjxH2]
	movaps xmm4, [rsp + nb234_fjyH2]
	movaps xmm5, [rsp + nb234_fjzH2]
	mulps xmm0, [rsp + nb234_dxMH2]
	mulps xmm1, [rsp + nb234_dyMH2]
	mulps xmm2, [rsp + nb234_dzMH2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb234_fixM]
	addps xmm1, [rsp + nb234_fiyM]
	addps xmm2, [rsp + nb234_fizM]
	movaps [rsp + nb234_fjxH2], xmm3
	movaps [rsp + nb234_fjyH2], xmm4
	movaps [rsp + nb234_fjzH2], xmm5
	movaps [rsp + nb234_fixM], xmm0
	movaps [rsp + nb234_fiyM], xmm1
	movaps [rsp + nb234_fizM], xmm2

	;# M-M interaction 
	movaps xmm0, [rsp + nb234_rinvMM]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb234_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb234_rsqMM] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb234_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb234_qqMM] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [rsp + nb234_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [rsp + nb234_qqMM] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm1, xmm0
	movaps [rsp + nb234_vctot], xmm6
	movaps xmm2, xmm0
	
	movaps xmm3, [rsp + nb234_fjxM]
	movaps xmm4, [rsp + nb234_fjyM]
	movaps xmm5, [rsp + nb234_fjzM]
	mulps xmm0, [rsp + nb234_dxMM]
	mulps xmm1, [rsp + nb234_dyMM]
	mulps xmm2, [rsp + nb234_dzMM]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [rsp + nb234_fixM]
	addps xmm1, [rsp + nb234_fiyM]
	addps xmm2, [rsp + nb234_fizM]
	movaps [rsp + nb234_fjxM], xmm3
	movaps [rsp + nb234_fjyM], xmm4
	movaps [rsp + nb234_fjzM], xmm5
	movaps [rsp + nb234_fixM], xmm0
	movaps [rsp + nb234_fiyM], xmm1
	movaps [rsp + nb234_fizM], xmm2

	mov rdi, [rbp + nb234_faction]
	;# update j forces 
	;# 4 j waters with four atoms each.
	;# step 1 : transpose fjxO, fjyO, fjzO, fjxH1
	movaps xmm0, [rsp + nb234_fjxO]
	movaps xmm1, [rsp + nb234_fjyO]
	movaps xmm2, [rsp + nb234_fjzO]
	movaps xmm3, [rsp + nb234_fjxH1]
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
	movlps   xmm1, [rdi + rax*4]
	movlps   xmm5, [rdi + rbx*4]
	movlps   xmm6, [rdi + rcx*4]
	movlps   xmm7, [rdi + rdx*4]
	movhps   xmm1, [rdi + rax*4 + 8]
	movhps   xmm5, [rdi + rbx*4 + 8]
	movhps   xmm6, [rdi + rcx*4 + 8]
	movhps   xmm7, [rdi + rdx*4 + 8]
	;# add
	addps    xmm1, xmm4
	addps    xmm5, xmm2
	addps    xmm6, xmm0
	addps    xmm7, xmm3
	;# store back
	movlps   [rdi + rax*4], xmm1
	movlps   [rdi + rbx*4], xmm5
	movlps   [rdi + rcx*4], xmm6
	movlps   [rdi + rdx*4], xmm7
	movhps   [rdi + rax*4 + 8], xmm1
	movhps   [rdi + rbx*4 + 8], xmm5
	movhps   [rdi + rcx*4 + 8], xmm6
	movhps   [rdi + rdx*4 + 8], xmm7

	;# step 2 : transpose fjyH1, fjzH1, fjxH2, fjyH2
	movaps xmm0, [rsp + nb234_fjyH1]
	movaps xmm1, [rsp + nb234_fjzH1]
	movaps xmm2, [rsp + nb234_fjxH2]
	movaps xmm3, [rsp + nb234_fjyH2]
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
	movlps   xmm1, [rdi + rax*4 + 16]
	movlps   xmm5, [rdi + rbx*4 + 16]
	movlps   xmm6, [rdi + rcx*4 + 16]
	movlps   xmm7, [rdi + rdx*4 + 16]
	movhps   xmm1, [rdi + rax*4 + 24]
	movhps   xmm5, [rdi + rbx*4 + 24]
	movhps   xmm6, [rdi + rcx*4 + 24]
	movhps   xmm7, [rdi + rdx*4 + 24]
	;# add
	addps    xmm1, xmm4
	addps    xmm5, xmm2
	addps    xmm6, xmm0
	addps    xmm7, xmm3
	;# store back
	movlps   [rdi + rax*4 + 16], xmm1
	movlps   [rdi + rbx*4 + 16], xmm5
	movlps   [rdi + rcx*4 + 16], xmm6
	movlps   [rdi + rdx*4 + 16], xmm7
	movhps   [rdi + rax*4 + 24], xmm1
	movhps   [rdi + rbx*4 + 24], xmm5
	movhps   [rdi + rcx*4 + 24], xmm6
	movhps   [rdi + rdx*4 + 24], xmm7

	;# step 3 : transpose fjzH2, fjxM, fjyM, fjzM. xmm4 is scratch
	movaps xmm0, [rsp + nb234_fjzH2]
	movaps xmm1, [rsp + nb234_fjxM]
	movaps xmm2, [rsp + nb234_fjyM]
	movaps xmm3, [rsp + nb234_fjzM]
	
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
	movlps   xmm1, [rdi + rax*4 + 32]
	movlps   xmm5, [rdi + rbx*4 + 32]
	movlps   xmm6, [rdi + rcx*4 + 32]
	movlps   xmm7, [rdi + rdx*4 + 32]
	movhps   xmm1, [rdi + rax*4 + 40]
	movhps   xmm5, [rdi + rbx*4 + 40]
	movhps   xmm6, [rdi + rcx*4 + 40]
	movhps   xmm7, [rdi + rdx*4 + 40]
	;# add
	addps    xmm1, xmm4
	addps    xmm5, xmm2
	addps    xmm6, xmm0
	addps    xmm7, xmm3
	;# store back
	movlps   [rdi + rax*4 + 32], xmm1
	movlps   [rdi + rbx*4 + 32], xmm5
	movlps   [rdi + rcx*4 + 32], xmm6
	movlps   [rdi + rdx*4 + 32], xmm7
	movhps   [rdi + rax*4 + 40], xmm1
	movhps   [rdi + rbx*4 + 40], xmm5
	movhps   [rdi + rcx*4 + 40], xmm6
	movhps   [rdi + rdx*4 + 40], xmm7
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb234_innerk],  4
	jl    .nb234_single_check
	jmp   .nb234_unroll_loop
.nb234_single_check:
	add dword ptr [rsp + nb234_innerk],  4
	jnz   .nb234_single_loop
	jmp   .nb234_updateouterdata
.nb234_single_loop:
	mov   rdx, [rsp + nb234_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb234_innerjjnr],  4	

	mov rsi, [rbp + nb234_pos]
	lea   rax, [rax + rax*2]  

	;# fetch j coordinates
	movlps xmm3,  [rsi + rax*4]		;#  Ox  Oy  
	movlps xmm4,  [rsi + rax*4 + 16]	;# H1y H1z 
	movlps xmm5,  [rsi + rax*4 + 32]	;# H2z  Mx 
	movhps xmm3,  [rsi + rax*4 + 8]   	;#  Ox  Oy  Oz H1x
	movhps xmm4,  [rsi + rax*4 + 24]	;# H1y H1z H2x H2y
	movhps xmm5,  [rsi + rax*4 + 40]	;# H2z  Mx  My  Mz
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
	movaps [rsp + nb234_jxO], xmm6
	movaps [rsp + nb234_jyO], xmm3
	movaps [rsp + nb234_jzO], xmm1

	;# do O and M in parallel
	movaps xmm0, [rsp + nb234_ixO]
	movaps xmm1, [rsp + nb234_iyO]
	movaps xmm2, [rsp + nb234_izO]
	movaps xmm3, [rsp + nb234_ixM]
	movaps xmm4, [rsp + nb234_iyM]
	movaps xmm5, [rsp + nb234_izM]
	subps  xmm0, [rsp + nb234_jxO]
	subps  xmm1, [rsp + nb234_jyO]
	subps  xmm2, [rsp + nb234_jzO]
	subps  xmm3, [rsp + nb234_jxO]
	subps  xmm4, [rsp + nb234_jyO]
	subps  xmm5, [rsp + nb234_jzO]
	
	movaps [rsp + nb234_dxOO], xmm0
	movaps [rsp + nb234_dyOO], xmm1
	movaps [rsp + nb234_dzOO], xmm2
	movaps [rsp + nb234_dxMM], xmm3
	movaps [rsp + nb234_dyMM], xmm4
	movaps [rsp + nb234_dzMM], xmm5
	
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
	;# Save data 
	movaps [rsp + nb234_rsqOO], xmm0
	movaps [rsp + nb234_rsqMM], xmm4
	
	;# do 1/x for O
	rsqrtps xmm1, xmm0
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb234_three]
	mulps   xmm1, xmm0	;# rsq*lu*lu 
	subps   xmm3, xmm1	;# constant 30-rsq*lu*lu 
	mulps   xmm3, xmm2	;# lu*(3-rsq*lu*lu) 
	mulps   xmm3, [rsp + nb234_half]
	movaps  [rsp + nb234_rinvOO], xmm3	;# rinvH2 
	
	;# 1/sqrt(x) for M
	rsqrtps xmm5, xmm4
	movaps  xmm6, xmm5	
	mulps   xmm5, xmm5
	movaps  xmm7, [rsp + nb234_three]
	mulps   xmm5, xmm4
	subps   xmm7, xmm5
	mulps   xmm7, xmm6
	mulps   xmm7, [rsp + nb234_half] ;# rinv iH1 - j water 
	movaps [rsp + nb234_rinvMM], xmm7


	;# LJ table interaction
	movaps xmm0, [rsp + nb234_rinvOO]
	movss xmm1, xmm0
	mulss  xmm1, [rsp + nb234_rsqOO] ;# xmm1=r 
	mulss  xmm1, [rsp + nb234_tsc]
		
    cvttps2pi mm6, xmm1
    cvtpi2ps xmm3, mm6
	subss    xmm1, xmm3	;# xmm1=eps 
    movss xmm2, xmm1
    mulss  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 3
	
    movd mm0, eax
	
    mov  rsi, [rbp + nb234_VFtab]
    movd eax, mm6
	
    ;# dispersion 
    movlps xmm5, [rsi + rax*4]
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# constant 10001000
    shufps xmm5, xmm7, 221  ;# constant 11011101

    movlps xmm7, [rsi + rax*4 + 8]
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# constant 10001000
    shufps xmm7, xmm3, 221  ;# constant 11011101
    ;# dispersion table ready, in xmm4-xmm7 
    mulss  xmm6, xmm1       ;# xmm6=Geps 
    mulss  xmm7, xmm2       ;# xmm7=Heps2 
    addss  xmm5, xmm6
    addss  xmm5, xmm7       ;# xmm5=Fp 
    mulss  xmm7, [rsp + nb234_two]       ;# two*Heps2 
    addss  xmm7, xmm6
    addss  xmm7, xmm5 ;# xmm7=FF 
    mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
    addss  xmm5, xmm4 ;# xmm5=VV 

    movss xmm4, [rsp + nb234_c6]
    mulss  xmm7, xmm4    ;# fijD 
    mulss  xmm5, xmm4    ;# Vvdw6 
	xorps  xmm3, xmm3
	
	mulps  xmm7, [rsp + nb234_tsc]
	subss  xmm3, xmm7
	movss  [rsp + nb234_fstmp], xmm3 

    addss  xmm5, [rsp + nb234_Vvdwtot]
    movss [rsp + nb234_Vvdwtot], xmm5

    ;# repulsion 
    movlps xmm5, [rsi + rax*4 + 16]
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# constant 10001000
    shufps xmm5, xmm7, 221  ;# constant 11011101

    movlps xmm7, [rsi + rax*4 + 24]
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# constant 10001000
    shufps xmm7, xmm3, 221  ;# constant 11011101
    ;# table ready, in xmm4-xmm7 
    mulss  xmm6, xmm1       ;# xmm6=Geps 
    mulss  xmm7, xmm2       ;# xmm7=Heps2 
    addss  xmm5, xmm6
    addss  xmm5, xmm7       ;# xmm5=Fp 
    mulss  xmm7, [rsp + nb234_two]       ;# two*Heps2 
    addss  xmm7, xmm6
    addss  xmm7, xmm5 ;# xmm7=FF 
    mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
    addss  xmm5, xmm4 ;# xmm5=VV 

    movss xmm4, [rsp + nb234_c12]
    mulss  xmm7, xmm4 ;# fijR 
    mulss  xmm5, xmm4 ;# Vvdw12 
	movaps xmm3, [rsp + nb234_fstmp]
	mulss  xmm7, [rsp + nb234_tsc]
	subss  xmm3, xmm7

    addss  xmm5, [rsp + nb234_Vvdwtot]
    movss [rsp + nb234_Vvdwtot], xmm5

	mulss xmm0, xmm3
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	movd  eax, mm0
	
	mulss  xmm0, [rsp + nb234_dxOO]
	mulss  xmm1, [rsp + nb234_dyOO]
	mulss  xmm2, [rsp + nb234_dzOO]
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2
	movaps  [rsp + nb234_fjxO], xmm3
	movaps  [rsp + nb234_fjyO], xmm4
	movaps  [rsp + nb234_fjzO], xmm5
	addss   xmm0, [rsp + nb234_fixO]
	addss   xmm1, [rsp + nb234_fiyO]
	addss   xmm2, [rsp + nb234_fizO]
	movss  [rsp + nb234_fixO], xmm0
	movss  [rsp + nb234_fiyO], xmm1
	movss  [rsp + nb234_fizO], xmm2

	;# do  M coulomb interaction
	movaps xmm0, [rsp + nb234_rinvMM]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb234_krf]
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [rsp + nb234_qqMH]
	movhps  xmm3, [rsp + nb234_qqMM]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	mulps  xmm5, [rsp + nb234_rsqMM] ;# xmm5=krsq 
	movaps xmm6, xmm5
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [rsp + nb234_crf]
	mulps  xmm6, xmm3 ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	mulps xmm5, [rsp + nb234_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, xmm3 ;# xmm7 = coul part of fscal 
	
	addps  xmm6, [rsp + nb234_vctot] 
	movaps [rsp + nb234_vctot], xmm6
	mulps  xmm0, xmm7

	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [rsp + nb234_dxMM]
	mulps   xmm1, [rsp + nb234_dyMM]
	mulps   xmm2, [rsp + nb234_dzMM]
	;# update forces M - j water 
	movaps  xmm3, [rsp + nb234_fjxO]
	movaps  xmm4, [rsp + nb234_fjyO]
	movaps  xmm5, [rsp + nb234_fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [rsp + nb234_fjxO], xmm3
	movaps  [rsp + nb234_fjyO], xmm4
	movaps  [rsp + nb234_fjzO], xmm5
	addps   xmm0, [rsp + nb234_fixM]
	addps   xmm1, [rsp + nb234_fiyM]
	addps   xmm2, [rsp + nb234_fizM]
	movaps  [rsp + nb234_fixM], xmm0
	movaps  [rsp + nb234_fiyM], xmm1
	movaps  [rsp + nb234_fizM], xmm2	
	
	;# i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [rsp + nb234_ixH1]
	movaps  xmm1, [rsp + nb234_iyH1]
	movaps  xmm2, [rsp + nb234_izH1]	
	movaps  xmm3, [rsp + nb234_ixH2] 
	movaps  xmm4, [rsp + nb234_iyH2] 
	movaps  xmm5, [rsp + nb234_izH2] 
	subps   xmm0, [rsp + nb234_jxO]
	subps   xmm1, [rsp + nb234_jyO]
	subps   xmm2, [rsp + nb234_jzO]
	subps   xmm3, [rsp + nb234_jxO]
	subps   xmm4, [rsp + nb234_jyO]
	subps   xmm5, [rsp + nb234_jzO]
	movaps [rsp + nb234_dxH1H1], xmm0
	movaps [rsp + nb234_dyH1H1], xmm1
	movaps [rsp + nb234_dzH1H1], xmm2
	movaps [rsp + nb234_dxH2H2], xmm3
	movaps [rsp + nb234_dyH2H2], xmm4
	movaps [rsp + nb234_dzH2H2], xmm5
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
	movaps  [rsp + nb234_rsqH1H1], xmm0
	movaps  [rsp + nb234_rsqH2H2], xmm4
		
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb234_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb234_half] ;# rinvH1H1
	mulps   xmm7, [rsp + nb234_half] ;# rinvH2H2
	movaps  [rsp + nb234_rinvH1H1], xmm3
	movaps  [rsp + nb234_rinvH2H2], xmm7
	
	;# Do H1 coulomb interaction
	movaps xmm0, [rsp + nb234_rinvH1H1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb234_krf]
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [rsp + nb234_qqHH]
	movhps  xmm3, [rsp + nb234_qqMH]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	mulps  xmm5, [rsp + nb234_rsqH1H1] ;# xmm5=krsq 
	movaps xmm6, xmm5
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [rsp + nb234_crf]
	mulps  xmm6, xmm3 ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	mulps xmm5, [rsp + nb234_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, xmm3 ;# xmm7 = coul part of fscal 
	
	addps  xmm6, [rsp + nb234_vctot] 
	movaps [rsp + nb234_vctot], xmm6

	mulps  xmm0, xmm7

	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [rsp + nb234_dxH1H1]
	mulps   xmm1, [rsp + nb234_dyH1H1]
	mulps   xmm2, [rsp + nb234_dzH1H1]
	;# update forces H1 - j water 
	movaps  xmm3, [rsp + nb234_fjxO]
	movaps  xmm4, [rsp + nb234_fjyO]
	movaps  xmm5, [rsp + nb234_fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [rsp + nb234_fjxO], xmm3
	movaps  [rsp + nb234_fjyO], xmm4
	movaps  [rsp + nb234_fjzO], xmm5
	addps   xmm0, [rsp + nb234_fixH1]
	addps   xmm1, [rsp + nb234_fiyH1]
	addps   xmm2, [rsp + nb234_fizH1]
	movaps  [rsp + nb234_fixH1], xmm0
	movaps  [rsp + nb234_fiyH1], xmm1
	movaps  [rsp + nb234_fizH1], xmm2	

	;# H2 Coulomb
	movaps xmm0, [rsp + nb234_rinvH2H2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb234_krf]
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [rsp + nb234_qqHH]
	movhps  xmm3, [rsp + nb234_qqMH]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	mulps  xmm5, [rsp + nb234_rsqH2H2] ;# xmm5=krsq 
	movaps xmm6, xmm5
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [rsp + nb234_crf]
	mulps  xmm6, xmm3 ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	mulps xmm5, [rsp + nb234_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, xmm3 ;# xmm7 = coul part of fscal 
	
	addps  xmm6, [rsp + nb234_vctot] ;# local vctot summation variable
	movaps [rsp + nb234_vctot], xmm6
	mulps  xmm0, xmm7

	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [rsp + nb234_dxH2H2]
	mulps   xmm1, [rsp + nb234_dyH2H2]
	mulps   xmm2, [rsp + nb234_dzH2H2]
	;# update forces H2 - j water 
	movaps  xmm3, [rsp + nb234_fjxO]
	movaps  xmm4, [rsp + nb234_fjyO]
	movaps  xmm5, [rsp + nb234_fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [rsp + nb234_fjxO], xmm3
	movaps  [rsp + nb234_fjyO], xmm4
	movaps  [rsp + nb234_fjzO], xmm5
	addps   xmm0, [rsp + nb234_fixH2]
	addps   xmm1, [rsp + nb234_fiyH2]
	addps   xmm2, [rsp + nb234_fizH2]
	movaps  [rsp + nb234_fixH2], xmm0
	movaps  [rsp + nb234_fiyH2], xmm1
	movaps  [rsp + nb234_fizH2], xmm2			
	
	mov     rsi, [rbp + nb234_faction]
	;# update j water forces from local variables.
	;# transpose back first
	movaps  xmm0, [rsp + nb234_fjxO] ;# Ox H1x H2x Mx 
	movaps  xmm1, [rsp + nb234_fjyO] ;# Oy H1y H2y My
	movaps  xmm2, [rsp + nb234_fjzO] ;# Oz H1z H2z Mz

	movaps  xmm3, xmm0
	movaps  xmm4, xmm0
	unpcklps xmm3, xmm1       	;# Ox Oy - -
	shufps  xmm4, xmm2, 0x1	  	;# h1x - Oz -
	movaps  xmm5, xmm1
	movaps  xmm6, xmm0
	unpcklps xmm5, xmm2 	  	;# - - H1y H1z
	unpckhps xmm6, xmm1	  	;# h2x h2y - - 
	unpckhps xmm1, xmm2	  	;# - - My Mz
	
	shufps   xmm2, xmm0, 50  ;# (00110010) h2z - Mx -
	shufps   xmm3, xmm4, 36  ;# constant 00100100 ;# Ox Oy Oz H1x 
	shufps   xmm5, xmm6, 78  ;# constant 01001110 ;# h1y h1z h2x h2y
	shufps   xmm2, xmm1, 232  ;# constant 11101000 ;# h2z mx my mz

	movlps  xmm0, [rsi + rax*4]
	movlps  xmm1, [rsi + rax*4 + 16]
	movlps  xmm4, [rsi + rax*4 + 32]
	movhps  xmm0, [rsi + rax*4 + 8]
	movhps  xmm1, [rsi + rax*4 + 24]
	movhps  xmm4, [rsi + rax*4 + 40]
	addps   xmm0, xmm3
	addps   xmm1, xmm5
	addps   xmm4, xmm2
	movlps   [rsi + rax*4], xmm0
	movlps   [rsi + rax*4 + 16], xmm1
	movlps   [rsi + rax*4 + 32], xmm4
	movhps   [rsi + rax*4 + 8], xmm0
	movhps   [rsi + rax*4 + 24], xmm1
	movhps   [rsi + rax*4 + 40], xmm4
	
	dec dword ptr [rsp + nb234_innerk]
	jz    .nb234_updateouterdata
	jmp   .nb234_single_loop
.nb234_updateouterdata:
	mov   ecx, [rsp + nb234_ii3]
	mov   rdi, [rbp + nb234_faction]
	mov   rsi, [rbp + nb234_fshift]
	mov   edx, [rsp + nb234_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb234_fixO]
	movaps xmm1, [rsp + nb234_fiyO] 
	movaps xmm2, [rsp + nb234_fizO]

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
	movss  xmm3, [rdi + rcx*4]
	movss  xmm4, [rdi + rcx*4 + 4]
	movss  xmm5, [rdi + rcx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [rdi + rcx*4],     xmm3
	movss  [rdi + rcx*4 + 4], xmm4
	movss  [rdi + rcx*4 + 8], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# constant 00001000	

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb234_fixH1]
	movaps xmm1, [rsp + nb234_fiyH1]
	movaps xmm2, [rsp + nb234_fizH1]

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
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb234_fixH2]
	movaps xmm1, [rsp + nb234_fiyH2]
	movaps xmm2, [rsp + nb234_fizH2]

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
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# accumulate Mi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb234_fixM]
	movaps xmm1, [rsp + nb234_fiyM]
	movaps xmm2, [rsp + nb234_fizM]

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
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# increment fshift force  
	movlps  xmm3, [rsi + rdx*4]
	movss  xmm4, [rsi + rdx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [rsi + rdx*4],    xmm3
	movss  [rsi + rdx*4 + 8], xmm4

	;# get n from stack
	mov esi, [rsp + nb234_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb234_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb234_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb234_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb234_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb234_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb234_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb234_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb234_n], esi
        jmp .nb234_outer
.nb234_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb234_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb234_end
        ;# non-zero, do one more workunit
        jmp   .nb234_threadloop
.nb234_end:
	mov eax, [rsp + nb234_nouter]
	mov ebx, [rsp + nb234_ninner]
	mov rcx, [rbp + nb234_outeriter]
	mov rdx, [rbp + nb234_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 1912
	femms

	pop rbx
	pop	rbp
	ret


	
.globl nb_kernel234nf_x86_64_sse
.globl _nb_kernel234nf_x86_64_sse
nb_kernel234nf_x86_64_sse:	
_nb_kernel234nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb234nf_fshift,         16
.equiv          nb234nf_gid,            24
.equiv          nb234nf_pos,            32
.equiv          nb234nf_faction,        40
.equiv          nb234nf_charge,         48
.equiv          nb234nf_p_facel,        56
.equiv          nb234nf_argkrf,         64
.equiv          nb234nf_argcrf,         72
.equiv          nb234nf_Vc,             80
.equiv          nb234nf_type,           88
.equiv          nb234nf_p_ntype,        96
.equiv          nb234nf_vdwparam,       104
.equiv          nb234nf_Vvdw,           112
.equiv          nb234nf_p_tabscale,     120
.equiv          nb234nf_VFtab,          128
.equiv          nb234nf_invsqrta,       136
.equiv          nb234nf_dvda,           144
.equiv          nb234nf_p_gbtabscale,   152
.equiv          nb234nf_GBtab,          160
.equiv          nb234nf_p_nthreads,     168
.equiv          nb234nf_count,          176
.equiv          nb234nf_mtx,            184
.equiv          nb234nf_outeriter,      192
.equiv          nb234nf_inneriter,      200
.equiv          nb234nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb234nf_ixO,            0
.equiv          nb234nf_iyO,            16
.equiv          nb234nf_izO,            32
.equiv          nb234nf_ixH1,           48
.equiv          nb234nf_iyH1,           64
.equiv          nb234nf_izH1,           80
.equiv          nb234nf_ixH2,           96
.equiv          nb234nf_iyH2,           112
.equiv          nb234nf_izH2,           128
.equiv          nb234nf_ixM,            144
.equiv          nb234nf_iyM,            160
.equiv          nb234nf_izM,            176
.equiv          nb234nf_jxO,            192
.equiv          nb234nf_jyO,            208
.equiv          nb234nf_jzO,            224
.equiv          nb234nf_jxH1,           240
.equiv          nb234nf_jyH1,           256
.equiv          nb234nf_jzH1,           272
.equiv          nb234nf_jxH2,           288
.equiv          nb234nf_jyH2,           304
.equiv          nb234nf_jzH2,           320
.equiv          nb234nf_jxM,            336
.equiv          nb234nf_jyM,            352
.equiv          nb234nf_jzM,            368
.equiv          nb234nf_qqMM,           384
.equiv          nb234nf_qqMH,           400
.equiv          nb234nf_qqHH,           416
.equiv          nb234nf_tsc,            432
.equiv          nb234nf_c6,             448
.equiv          nb234nf_c12,            464
.equiv          nb234nf_vctot,          480
.equiv          nb234nf_Vvdwtot,        496
.equiv          nb234nf_half,           512
.equiv          nb234nf_three,          528
.equiv          nb234nf_rsqOO,          544
.equiv          nb234nf_rsqH1H1,        560
.equiv          nb234nf_rsqH1H2,        576
.equiv          nb234nf_rsqH1M,         592
.equiv          nb234nf_rsqH2H1,        608
.equiv          nb234nf_rsqH2H2,        624
.equiv          nb234nf_rsqH2M,         640
.equiv          nb234nf_rsqMH1,         656
.equiv          nb234nf_rsqMH2,         672
.equiv          nb234nf_rsqMM,          688
.equiv          nb234nf_rinvOO,         704
.equiv          nb234nf_rinvH1H1,       720
.equiv          nb234nf_rinvH1H2,       736
.equiv          nb234nf_rinvH1M,        752
.equiv          nb234nf_rinvH2H1,       768
.equiv          nb234nf_rinvH2H2,       784
.equiv          nb234nf_rinvH2M,        800
.equiv          nb234nf_rinvMH1,        816
.equiv          nb234nf_rinvMH2,        832
.equiv          nb234nf_rinvMM,         848
.equiv          nb234nf_krf,            864
.equiv          nb234nf_crf,            880
.equiv          nb234nf_is3,            896
.equiv          nb234nf_ii3,            900
.equiv          nb234nf_innerjjnr,      904
.equiv          nb234nf_nri,            912
.equiv          nb234nf_iinr,           920
.equiv          nb234nf_jindex,         928
.equiv          nb234nf_jjnr,           936
.equiv          nb234nf_shift,          944
.equiv          nb234nf_shiftvec,       952
.equiv          nb234nf_facel,          960
.equiv          nb234nf_innerk,         968
.equiv          nb234nf_n,              972
.equiv          nb234nf_nn1,            976
.equiv          nb234nf_nouter,         980
.equiv          nb234nf_ninner,         984


	push rbp
	mov  rbp, rsp
	push rbx
	
	femms
	sub rsp, 1000		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb234nf_nouter], eax
	mov [rsp + nb234nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb234nf_nri], edi
	mov [rsp + nb234nf_iinr], rsi
	mov [rsp + nb234nf_jindex], rdx
	mov [rsp + nb234nf_jjnr], rcx
	mov [rsp + nb234nf_shift], r8
	mov [rsp + nb234nf_shiftvec], r9
	mov rsi, [rbp + nb234nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb234nf_facel], xmm0

	mov rax, [rbp + nb234nf_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb234nf_tsc], xmm3

	mov rsi, [rbp + nb234nf_argkrf]
	mov rdi, [rbp + nb234nf_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb234nf_krf], xmm1
	movaps [rsp + nb234nf_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb234nf_half], eax
	movss xmm1, [rsp + nb234nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb234nf_half],  xmm1
	movaps [rsp + nb234nf_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb234nf_iinr]   ;# rcx = pointer into iinr[]
	mov   ebx, [rcx]		;# ebx =ii 

	mov   rdx, [rbp + nb234nf_charge]
	movss xmm5, [rdx + rbx*4 + 4]	
	movss xmm3, [rdx + rbx*4 + 12]	
	movss xmm4, xmm3	
	mov rsi, [rbp + nb234nf_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb234nf_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb234nf_qqMM], xmm3
	movaps [rsp + nb234nf_qqMH], xmm4
	movaps [rsp + nb234nf_qqHH], xmm5
	
	xorps xmm0, xmm0
	mov   rdx, [rbp + nb234nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb234nf_p_ntype]
	imul  ecx, [rdi]  ;# rcx = ntia = 2*ntype*type[ii0] 
	add   rdx, rcx
	mov   rax, [rbp + nb234nf_vdwparam]
	movlps xmm0, [rax + rdx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0x55
	movaps [rsp + nb234nf_c6], xmm0
	movaps [rsp + nb234nf_c12], xmm1

.nb234nf_threadloop:
        mov   rsi, [rbp + nb234nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb234nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [rsi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb234nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb234nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb234nf_n], eax
        mov [rsp + nb234nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb234nf_outerstart
        jmp .nb234nf_end
	
.nb234nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb234nf_nouter]
	mov [rsp + nb234nf_nouter], ebx

.nb234nf_outer:
	mov   rax, [rsp + nb234nf_shift]  	;# eax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb234nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb234nf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb234nf_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb234nf_pos]	;# eax = base of pos[]  
	mov   [rsp + nb234nf_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	movaps xmm6, xmm0
	movaps xmm7, xmm1

	addss xmm3, [rax + rbx*4]  	;# ox
	addss xmm4, [rax + rbx*4 + 4]  ;# oy
	addss xmm5, [rax + rbx*4 + 8]  ;# oz
	addss xmm6, [rax + rbx*4 + 12] ;# h1x
	addss xmm7, [rax + rbx*4 + 16] ;# h1y
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	shufps xmm7, xmm7, 0
	movaps [rsp + nb234nf_ixO], xmm3
	movaps [rsp + nb234nf_iyO], xmm4
	movaps [rsp + nb234nf_izO], xmm5
	movaps [rsp + nb234nf_ixH1], xmm6
	movaps [rsp + nb234nf_iyH1], xmm7

	movss xmm6, xmm2
	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm6, [rax + rbx*4 + 20] ;# h1z
	addss xmm0, [rax + rbx*4 + 24] ;# h2x
	addss xmm1, [rax + rbx*4 + 28] ;# h2y
	addss xmm2, [rax + rbx*4 + 32] ;# h2z
	addss xmm3, [rax + rbx*4 + 36] ;# mx
	addss xmm4, [rax + rbx*4 + 40] ;# my
	addss xmm5, [rax + rbx*4 + 44] ;# mz

	shufps xmm6, xmm6, 0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb234nf_izH1], xmm6
	movaps [rsp + nb234nf_ixH2], xmm0
	movaps [rsp + nb234nf_iyH2], xmm1
	movaps [rsp + nb234nf_izH2], xmm2
	movaps [rsp + nb234nf_ixM], xmm3
	movaps [rsp + nb234nf_iyM], xmm4
	movaps [rsp + nb234nf_izM], xmm5

	;# clear vctot 
	xorps xmm4, xmm4
	movaps [rsp + nb234nf_vctot], xmm4
	movaps [rsp + nb234nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb234nf_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb234nf_pos]
	mov   rax, [rsp + nb234nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb234nf_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb234nf_ninner]
	mov   [rsp + nb234nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb234nf_innerk], edx	;# number of innerloop atoms 
	jge   .nb234nf_unroll_loop
	jmp   .nb234nf_single_check
.nb234nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb234nf_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb234nf_innerjjnr],  16 ;# advance pointer (unroll 4) 

	mov rsi, [rbp + nb234nf_pos]   	;# base of pos[] 

	lea   rax, [rax + rax*2] 	;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2] 	;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	
	
	;# move j coordinates to local temp variables
	;# Load Ox, Oy, Oz, H1x 
	movlps xmm1, [rsi + rax*4]	;#  Oxa   Oya    -    -
	movlps xmm4, [rsi + rcx*4]	;#  Oxc   Oyc    -    -
	movhps xmm1, [rsi + rbx*4]	;#  Oxa   Oya   Oxb   Oyb 
	movhps xmm4, [rsi + rdx*4]	;#  Oxc   Oyc   Oxd   Oyd 
	movaps xmm0, xmm1		;#  Oxa   Oya   Oxb   Oyb 
	shufps xmm0, xmm4, 0x88		;#  Oxa   Oxb   Oxc   Oxd
	shufps xmm1, xmm4, 0xDD		;#  Oya   Oyb   Oyc   Oyd
	movlps xmm3, [rsi + rax*4 + 8]	;#  Oza  H1xa    -    -
	movlps xmm5, [rsi + rcx*4 + 8]	;#  Ozc  H1xc    -    -
	movhps xmm3, [rsi + rbx*4 + 8]	;#  Oza  H1xa   Ozb  H1xb 
	movhps xmm5, [rsi + rdx*4 + 8]	;#  Ozc  H1xc   Ozd  H1xd 
	movaps xmm2, xmm3		;#  Oza  H1xa   Ozb  H1xb 
	shufps xmm2, xmm5, 0x88		;#  Oza   Ozb   Ozc   Ozd
	shufps xmm3, xmm5, 0xDD		;# H1xa  H1xb  H1xc  H1xd
	;# coordinates in xmm0-xmm3	
	;# store
	movaps [rsp + nb234nf_jxO], xmm0
	movaps [rsp + nb234nf_jyO], xmm1
	movaps [rsp + nb234nf_jzO], xmm2
	movaps [rsp + nb234nf_jxH1], xmm3

	;# Load H1y H1z H2x H2y 
	movlps xmm1, [rsi + rax*4 + 16]	
	movlps xmm4, [rsi + rcx*4 + 16]	
	movhps xmm1, [rsi + rbx*4 + 16]	
	movhps xmm4, [rsi + rdx*4 + 16]	
	movaps xmm0, xmm1		
	shufps xmm0, xmm4, 0x88		
	shufps xmm1, xmm4, 0xDD		
	movlps xmm3, [rsi + rax*4 + 24]	
	movlps xmm5, [rsi + rcx*4 + 24]	
	movhps xmm3, [rsi + rbx*4 + 24]	
	movhps xmm5, [rsi + rdx*4 + 24]	
	movaps xmm2, xmm3		
	shufps xmm2, xmm5, 0x88		
	shufps xmm3, xmm5, 0xDD		
	;# coordinates in xmm0-xmm3	
	;# store
	movaps [rsp + nb234nf_jyH1], xmm0
	movaps [rsp + nb234nf_jzH1], xmm1
	movaps [rsp + nb234nf_jxH2], xmm2
	movaps [rsp + nb234nf_jyH2], xmm3

	;# Load H2z Mx My Mz 
	movlps xmm1, [rsi + rax*4 + 32]	
	movlps xmm4, [rsi + rcx*4 + 32]	
	movhps xmm1, [rsi + rbx*4 + 32]	
	movhps xmm4, [rsi + rdx*4 + 32]	
	movaps xmm0, xmm1		
	shufps xmm0, xmm4, 0x88		
	shufps xmm1, xmm4, 0xDD		
	movlps xmm3, [rsi + rax*4 + 40]	
	movlps xmm5, [rsi + rcx*4 + 40]	
	movhps xmm3, [rsi + rbx*4 + 40]	
	movhps xmm5, [rsi + rdx*4 + 40]	
	movaps xmm2, xmm3		
	shufps xmm2, xmm5, 0x88		
	shufps xmm3, xmm5, 0xDD		
	;# coordinates in xmm0-xmm3	
	;# store
	movaps [rsp + nb234nf_jzH2], xmm0
	movaps [rsp + nb234nf_jxM], xmm1
	movaps [rsp + nb234nf_jyM], xmm2
	movaps [rsp + nb234nf_jzM], xmm3
	
	;# start calculating pairwise distances
	movaps xmm0, [rsp + nb234nf_ixO]
	movaps xmm1, [rsp + nb234nf_iyO]
	movaps xmm2, [rsp + nb234nf_izO]
	movaps xmm3, [rsp + nb234nf_ixH1]
	movaps xmm4, [rsp + nb234nf_iyH1]
	movaps xmm5, [rsp + nb234nf_izH1]
	subps  xmm0, [rsp + nb234nf_jxO]
	subps  xmm1, [rsp + nb234nf_jyO]
	subps  xmm2, [rsp + nb234nf_jzO]
	subps  xmm3, [rsp + nb234nf_jxH1]
	subps  xmm4, [rsp + nb234nf_jyH1]
	subps  xmm5, [rsp + nb234nf_jzH1]
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
	movaps [rsp + nb234nf_rsqOO], xmm0
	movaps [rsp + nb234nf_rsqH1H1], xmm3

	movaps xmm0, [rsp + nb234nf_ixH1]
	movaps xmm1, [rsp + nb234nf_iyH1]
	movaps xmm2, [rsp + nb234nf_izH1]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [rsp + nb234nf_jxH2]
	subps  xmm1, [rsp + nb234nf_jyH2]
	subps  xmm2, [rsp + nb234nf_jzH2]
	subps  xmm3, [rsp + nb234nf_jxM]
	subps  xmm4, [rsp + nb234nf_jyM]
	subps  xmm5, [rsp + nb234nf_jzM]
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
	movaps [rsp + nb234nf_rsqH1H2], xmm0
	movaps [rsp + nb234nf_rsqH1M], xmm3

	movaps xmm0, [rsp + nb234nf_ixH2]
	movaps xmm1, [rsp + nb234nf_iyH2]
	movaps xmm2, [rsp + nb234nf_izH2]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [rsp + nb234nf_jxH1]
	subps  xmm1, [rsp + nb234nf_jyH1]
	subps  xmm2, [rsp + nb234nf_jzH1]
	subps  xmm3, [rsp + nb234nf_jxH2]
	subps  xmm4, [rsp + nb234nf_jyH2]
	subps  xmm5, [rsp + nb234nf_jzH2]
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
	movaps [rsp + nb234nf_rsqH2H1], xmm0
	movaps [rsp + nb234nf_rsqH2H2], xmm3

	movaps xmm0, [rsp + nb234nf_ixH2]
	movaps xmm1, [rsp + nb234nf_iyH2]
	movaps xmm2, [rsp + nb234nf_izH2]
	movaps xmm3, [rsp + nb234nf_ixM]
	movaps xmm4, [rsp + nb234nf_iyM]
	movaps xmm5, [rsp + nb234nf_izM]
	subps  xmm0, [rsp + nb234nf_jxM]
	subps  xmm1, [rsp + nb234nf_jyM]
	subps  xmm2, [rsp + nb234nf_jzM]
	subps  xmm3, [rsp + nb234nf_jxH1]
	subps  xmm4, [rsp + nb234nf_jyH1]
	subps  xmm5, [rsp + nb234nf_jzH1]
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
	movaps [rsp + nb234nf_rsqH2M], xmm0
	movaps [rsp + nb234nf_rsqMH1], xmm4

	movaps xmm0, [rsp + nb234nf_ixM]
	movaps xmm1, [rsp + nb234nf_iyM]
	movaps xmm2, [rsp + nb234nf_izM]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [rsp + nb234nf_jxH2]
	subps  xmm1, [rsp + nb234nf_jyH2]
	subps  xmm2, [rsp + nb234nf_jzH2]
	subps  xmm3, [rsp + nb234nf_jxM]
	subps  xmm4, [rsp + nb234nf_jyM]
	subps  xmm5, [rsp + nb234nf_jzM]
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
	movaps [rsp + nb234nf_rsqMH2], xmm0
	movaps [rsp + nb234nf_rsqMM], xmm4

	;# start by doing invsqrt for OO
	rsqrtps xmm1, [rsp + nb234nf_rsqOO]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb234nf_three]
	mulps   xmm1, [rsp + nb234nf_rsqOO]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + nb234nf_half] 
	movaps  [rsp + nb234nf_rinvOO], xmm3
	
	;# more invsqrt ops - do two at a time.
	rsqrtps xmm1, [rsp + nb234nf_rsqH1H1]
	rsqrtps xmm5, [rsp + nb234nf_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb234nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb234nf_rsqH1H1]
	mulps   xmm5, [rsp + nb234nf_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb234nf_half] ;# rinvH1H1 
	mulps   xmm7, [rsp + nb234nf_half] ;# rinvH1H2 
	movaps  [rsp + nb234nf_rinvH1H1], xmm3
	movaps  [rsp + nb234nf_rinvH1H2], xmm7
			
	rsqrtps xmm1, [rsp + nb234nf_rsqH1M]
	rsqrtps xmm5, [rsp + nb234nf_rsqH2H1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb234nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb234nf_rsqH1M]
	mulps   xmm5, [rsp + nb234nf_rsqH2H1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb234nf_half] 
	mulps   xmm7, [rsp + nb234nf_half]
	movaps  [rsp + nb234nf_rinvH1M], xmm3
	movaps  [rsp + nb234nf_rinvH2H1], xmm7
			
	rsqrtps xmm1, [rsp + nb234nf_rsqH2H2]
	rsqrtps xmm5, [rsp + nb234nf_rsqH2M]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb234nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb234nf_rsqH2H2]
	mulps   xmm5, [rsp + nb234nf_rsqH2M]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb234nf_half] 
	mulps   xmm7, [rsp + nb234nf_half]
	movaps  [rsp + nb234nf_rinvH2H2], xmm3
	movaps  [rsp + nb234nf_rinvH2M], xmm7
	
	rsqrtps xmm1, [rsp + nb234nf_rsqMH1]
	rsqrtps xmm5, [rsp + nb234nf_rsqMH2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb234nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb234nf_rsqMH1]
	mulps   xmm5, [rsp + nb234nf_rsqMH2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb234nf_half] 
	mulps   xmm7, [rsp + nb234nf_half]
	movaps  [rsp + nb234nf_rinvMH1], xmm3
	movaps  [rsp + nb234nf_rinvMH2], xmm7
        		
	rsqrtps xmm1, [rsp + nb234nf_rsqMM]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb234nf_three]
	mulps   xmm1, [rsp + nb234nf_rsqMM]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + nb234nf_half] 
	movaps  [rsp + nb234nf_rinvMM], xmm3
	
	;# start with OO LJ interaction
	movaps xmm0, [rsp + nb234nf_rinvOO]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb234nf_rsqOO] ;# xmm1=r 
	mulps  xmm1, [rsp + nb234nf_tsc]
		
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

    mov  rsi, [rbp + nb234nf_VFtab]
	
    ;# dispersion 
    movlps xmm5, [rsi + rax*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm5, [rsi + rbx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# constant 10001000
    shufps xmm5, xmm7, 221  ;# constant 11011101

    movlps xmm7, [rsi + rax*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm7, [rsi + rbx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of table  
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

    movaps xmm4, [rsp + nb234nf_c6]
    mulps  xmm5, xmm4    ;# Vvdw6 

    addps  xmm5, [rsp + nb234nf_Vvdwtot]
    movaps [rsp + nb234nf_Vvdwtot], xmm5

    ;# repulsion 
    movlps xmm5, [rsi + rax*4 + 16]
    movlps xmm7, [rsi + rcx*4 + 16]
    movhps xmm5, [rsi + rbx*4 + 16]
    movhps xmm7, [rsi + rdx*4 + 16] ;# got half table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# constant 10001000
    shufps xmm5, xmm7, 221  ;# constant 11011101

    movlps xmm7, [rsi + rax*4 + 24]
    movlps xmm3, [rsi + rcx*4 + 24]
    movhps xmm7, [rsi + rbx*4 + 24]
    movhps xmm3, [rsi + rdx*4 + 24] ;# other half of table  
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
 
    movaps xmm4, [rsp + nb234nf_c12]
    mulps  xmm5, xmm4 ;# Vvdw12 

    addps  xmm5, [rsp + nb234nf_Vvdwtot]
    movaps [rsp + nb234nf_Vvdwtot], xmm5

	;# Coulomb interactions 
	;# start with H1-H1 interaction 
	movaps xmm0, [rsp + nb234nf_rinvH1H1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb234nf_krf]
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	mulps  xmm5, [rsp + nb234nf_rsqH1H1] ;# xmm5=krsq 
	movaps xmm6, xmm5
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [rsp + nb234nf_crf]
	mulps  xmm6, [rsp + nb234nf_qqHH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addps  xmm6, [rsp + nb234nf_vctot] ;# local vctot summation variable 

	;# H1-H2 interaction 
	movaps xmm0, [rsp + nb234nf_rinvH1H2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb234nf_krf]
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb234nf_rsqH1H2] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb234nf_crf]
	mulps  xmm0, xmm0
	mulps  xmm4, [rsp + nb234nf_qqHH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	addps  xmm6, xmm4	;# add to local vctot 

	;# H1-M interaction  
	movaps xmm0, [rsp + nb234nf_rinvH1M]
	movaps xmm7, xmm0	;# xmm7=Rinv 
	movaps xmm5, [rsp + nb234nf_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb234nf_rsqH1M] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb234nf_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb234nf_qqMH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	addps  xmm6, xmm4	;# add to local vctot 

	;# H2-H1 interaction 
	movaps xmm0, [rsp + nb234nf_rinvH2H1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb234nf_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb234nf_rsqH2H1] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb234nf_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb234nf_qqHH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	addps  xmm6, xmm4	;# add to local vctot 

	;# H2-H2 interaction 
	movaps xmm0, [rsp + nb234nf_rinvH2H2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb234nf_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb234nf_rsqH2H2] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb234nf_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb234nf_qqHH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	addps  xmm6, xmm4	;# add to local vctot 

	;# H2-M interaction 
	movaps xmm0, [rsp + nb234nf_rinvH2M]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb234nf_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb234nf_rsqH2M] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb234nf_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb234nf_qqMH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	addps  xmm6, xmm4	;# add to local vctot 

	;# M-H1 interaction 
	movaps xmm0, [rsp + nb234nf_rinvMH1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb234nf_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb234nf_rsqMH1] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb234nf_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb234nf_qqMH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	addps  xmm6, xmm4	;# add to local vctot 

	;# M-H2 interaction 
	movaps xmm0, [rsp + nb234nf_rinvMH2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb234nf_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb234nf_rsqMH2] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb234nf_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb234nf_qqMH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	addps  xmm6, xmm4	;# add to local vctot 

	;# M-M interaction 
	movaps xmm0, [rsp + nb234nf_rinvMM]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb234nf_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb234nf_rsqMM] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb234nf_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb234nf_qqMM] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	addps  xmm6, xmm4	;# add to local vctot 
	movaps [rsp + nb234nf_vctot], xmm6

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb234nf_innerk],  4
	jl    .nb234nf_single_check
	jmp   .nb234nf_unroll_loop
.nb234nf_single_check:
	add dword ptr [rsp + nb234nf_innerk],  4
	jnz   .nb234nf_single_loop
	jmp   .nb234nf_updateouterdata
.nb234nf_single_loop:
	mov   rdx, [rsp + nb234nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb234nf_innerjjnr],  4

	mov rsi, [rbp + nb234nf_pos]
	lea   rax, [rax + rax*2]  

	;# fetch j coordinates
	movlps xmm3,  [rsi + rax*4]		;#  Ox  Oy  
	movlps xmm4,  [rsi + rax*4 + 16]	;# H1y H1z 
	movlps xmm5,  [rsi + rax*4 + 32]	;# H2z  Mx 
	movhps xmm3,  [rsi + rax*4 + 8]   	;#  Ox  Oy  Oz H1x
	movhps xmm4,  [rsi + rax*4 + 24]	;# H1y H1z H2x H2y
	movhps xmm5,  [rsi + rax*4 + 40]	;# H2z  Mx  My  Mz
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
	movaps [rsp + nb234nf_jxO], xmm6
	movaps [rsp + nb234nf_jyO], xmm3
	movaps [rsp + nb234nf_jzO], xmm1

	;# do O and M in parallel
	movaps xmm0, [rsp + nb234nf_ixO]
	movaps xmm1, [rsp + nb234nf_iyO]
	movaps xmm2, [rsp + nb234nf_izO]
	movaps xmm3, [rsp + nb234nf_ixM]
	movaps xmm4, [rsp + nb234nf_iyM]
	movaps xmm5, [rsp + nb234nf_izM]
	subps  xmm0, [rsp + nb234nf_jxO]
	subps  xmm1, [rsp + nb234nf_jyO]
	subps  xmm2, [rsp + nb234nf_jzO]
	subps  xmm3, [rsp + nb234nf_jxO]
	subps  xmm4, [rsp + nb234nf_jyO]
	subps  xmm5, [rsp + nb234nf_jzO]
	
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
	;# Save data 
	movaps [rsp + nb234nf_rsqOO], xmm0
	movaps [rsp + nb234nf_rsqMM], xmm4
	
	;# do 1/x for O
	rsqrtps xmm1, xmm0
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb234nf_three]
	mulps   xmm1, xmm0	;# rsq*lu*lu 
	subps   xmm3, xmm1	;# constant 30-rsq*lu*lu 
	mulps   xmm3, xmm2	;# lu*(3-rsq*lu*lu) 
	mulps   xmm3, [rsp + nb234nf_half]
	movaps  [rsp + nb234nf_rinvOO], xmm3	;# rinvH2 
	
	;# 1/sqrt(x) for M
	rsqrtps xmm5, xmm4
	movaps  xmm6, xmm5	
	mulps   xmm5, xmm5
	movaps  xmm7, [rsp + nb234nf_three]
	mulps   xmm5, xmm4
	subps   xmm7, xmm5
	mulps   xmm7, xmm6
	mulps   xmm7, [rsp + nb234nf_half] ;# rinv iH1 - j water 
	movaps [rsp + nb234nf_rinvMM], xmm7


	;# LJ table interaction
	movaps xmm0, [rsp + nb234nf_rinvOO]
	movss xmm1, xmm0
	mulss  xmm1, [rsp + nb234nf_rsqOO] ;# xmm1=r 
	mulss  xmm1, [rsp + nb234nf_tsc]
		
    cvttps2pi mm6, xmm1
    cvtpi2ps xmm3, mm6
	subss    xmm1, xmm3	;# xmm1=eps 
    movss xmm2, xmm1
    mulss  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 3
	
    mov  rsi, [rbp + nb234nf_VFtab]
    movd eax, mm6
	
    ;# dispersion 
    movlps xmm5, [rsi + rax*4]
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# constant 10001000
    shufps xmm5, xmm7, 221  ;# constant 11011101

    movlps xmm7, [rsi + rax*4 + 8]
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

    movss xmm4, [rsp + nb234nf_c6]
    mulss  xmm5, xmm4    ;# Vvdw6 

    addss  xmm5, [rsp + nb234nf_Vvdwtot]
    movss [rsp + nb234nf_Vvdwtot], xmm5

    ;# repulsion 
    movlps xmm5, [rsi + rax*4 + 16]
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# constant 10001000
    shufps xmm5, xmm7, 221  ;# constant 11011101

    movlps xmm7, [rsi + rax*4 + 24]
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

    movss xmm4, [rsp + nb234nf_c12]
    mulss  xmm5, xmm4 ;# Vvdw12 
    addss  xmm5, [rsp + nb234nf_Vvdwtot]
    movss [rsp + nb234nf_Vvdwtot], xmm5

	;# do  M coulomb interaction
	movaps xmm0, [rsp + nb234nf_rinvMM]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb234nf_krf]
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [rsp + nb234nf_qqMH]
	movhps  xmm3, [rsp + nb234nf_qqMM]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	mulps  xmm5, [rsp + nb234nf_rsqMM] ;# xmm5=krsq 
	movaps xmm6, xmm5
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [rsp + nb234nf_crf]
	mulps  xmm6, xmm3 ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addps  xmm6, [rsp + nb234nf_vctot] 
	movaps [rsp + nb234nf_vctot], xmm6
	
	;# i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [rsp + nb234nf_ixH1]
	movaps  xmm1, [rsp + nb234nf_iyH1]
	movaps  xmm2, [rsp + nb234nf_izH1]	
	movaps  xmm3, [rsp + nb234nf_ixH2] 
	movaps  xmm4, [rsp + nb234nf_iyH2] 
	movaps  xmm5, [rsp + nb234nf_izH2] 
	subps   xmm0, [rsp + nb234nf_jxO]
	subps   xmm1, [rsp + nb234nf_jyO]
	subps   xmm2, [rsp + nb234nf_jzO]
	subps   xmm3, [rsp + nb234nf_jxO]
	subps   xmm4, [rsp + nb234nf_jyO]
	subps   xmm5, [rsp + nb234nf_jzO]
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
	movaps  [rsp + nb234nf_rsqH1H1], xmm0
	movaps  [rsp + nb234nf_rsqH2H2], xmm4
		
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb234nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb234nf_half] ;# rinvH1H1
	mulps   xmm7, [rsp + nb234nf_half] ;# rinvH2H2
	movaps  [rsp + nb234nf_rinvH1H1], xmm3
	movaps  [rsp + nb234nf_rinvH2H2], xmm7
	
	;# Do H1 coulomb interaction
	movaps xmm0, [rsp + nb234nf_rinvH1H1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb234nf_krf]
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [rsp + nb234nf_qqHH]
	movhps  xmm3, [rsp + nb234nf_qqMH]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	mulps  xmm5, [rsp + nb234nf_rsqH1H1] ;# xmm5=krsq 
	movaps xmm6, xmm5
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [rsp + nb234nf_crf]
	mulps  xmm6, xmm3 ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addps  xmm6, [rsp + nb234nf_vctot] 
	movaps [rsp + nb234nf_vctot], xmm6

	;# H2 Coulomb
	movaps xmm0, [rsp + nb234nf_rinvH2H2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb234nf_krf]
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [rsp + nb234nf_qqHH]
	movhps  xmm3, [rsp + nb234nf_qqMH]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	mulps  xmm5, [rsp + nb234nf_rsqH2H2] ;# xmm5=krsq 
	movaps xmm6, xmm5
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [rsp + nb234nf_crf]
	mulps  xmm6, xmm3 ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	
	addps  xmm6, [rsp + nb234nf_vctot] ;# local vctot summation variable
	movaps [rsp + nb234nf_vctot], xmm6
	
	dec dword ptr [rsp + nb234nf_innerk]
	jz    .nb234nf_updateouterdata
	jmp   .nb234nf_single_loop
.nb234nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb234nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb234nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb234nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb234nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb234nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb234nf_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb234nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb234nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb234nf_n], esi
        jmp .nb234nf_outer
.nb234nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb234nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb234nf_end
        ;# non-zero, do one more workunit
        jmp   .nb234nf_threadloop
.nb234nf_end:
	mov eax, [rsp + nb234nf_nouter]
	mov ebx, [rsp + nb234nf_ninner]
	mov rcx, [rbp + nb234nf_outeriter]
	mov rdx, [rbp + nb234nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 1000
	femms

	pop rbx
	pop	rbp
	ret



