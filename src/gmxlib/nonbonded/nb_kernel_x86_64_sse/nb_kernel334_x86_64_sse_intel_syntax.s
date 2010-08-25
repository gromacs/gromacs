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


.globl nb_kernel334_x86_64_sse
.globl _nb_kernel334_x86_64_sse
nb_kernel334_x86_64_sse:	
_nb_kernel334_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb334_fshift,           16
.equiv          nb334_gid,              24
.equiv          nb334_pos,              32
.equiv          nb334_faction,          40
.equiv          nb334_charge,           48
.equiv          nb334_p_facel,          56
.equiv          nb334_argkrf,           64
.equiv          nb334_argcrf,           72
.equiv          nb334_Vc,               80
.equiv          nb334_type,             88
.equiv          nb334_p_ntype,          96
.equiv          nb334_vdwparam,         104
.equiv          nb334_Vvdw,             112
.equiv          nb334_p_tabscale,       120
.equiv          nb334_VFtab,            128
.equiv          nb334_invsqrta,         136
.equiv          nb334_dvda,             144
.equiv          nb334_p_gbtabscale,     152
.equiv          nb334_GBtab,            160
.equiv          nb334_p_nthreads,       168
.equiv          nb334_count,            176
.equiv          nb334_mtx,              184
.equiv          nb334_outeriter,        192
.equiv          nb334_inneriter,        200
.equiv          nb334_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
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
.equiv          nb334_epsH1,            1248
.equiv          nb334_epsH2,            1264
.equiv          nb334_epsM,             1280
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
.equiv          nb334_fstmp,            1744
.equiv          nb334_is3,              1760
.equiv          nb334_ii3,              1764
.equiv          nb334_nri,              1768
.equiv          nb334_iinr,             1776
.equiv          nb334_jindex,           1784
.equiv          nb334_jjnr,             1792
.equiv          nb334_shift,            1800
.equiv          nb334_shiftvec,         1808
.equiv          nb334_facel,            1816
.equiv          nb334_innerjjnr,        1824
.equiv          nb334_innerk,           1832
.equiv          nb334_n,                1836
.equiv          nb334_nn1,              1840
.equiv          nb334_nouter,           1844
.equiv          nb334_ninner,           1848
	push rbp
	mov  rbp, rsp
    
    ;# Push integer registers on stack
	push rbx
    push rsi
    push rdi
    push r12
    push r13
    push r14
    push r15

    ;# Make room for registers xmm6-xmm15 (10 registers=160 bytes)
    sub rsp, 168
    
    ;# Save xmm registers to stack
    movaps [rsp      ], xmm6
    movaps [rsp + 16 ], xmm7
    movaps [rsp + 32 ], xmm8
    movaps [rsp + 48 ], xmm9
    movaps [rsp + 64 ], xmm10
    movaps [rsp + 80 ], xmm11
    movaps [rsp + 96 ], xmm12
    movaps [rsp + 112], xmm13
    movaps [rsp + 128], xmm14
    movaps [rsp + 144], xmm15
    
	emms

	sub rsp, 1856		;# local variable stack space (n*16+8)
; .if 0    # block below only read by NASM - special calling convention on win64
%ifidn __OUTPUT_FORMAT__, win64
    ;# Adjust rbp to account for shadow space (32) & two extra args (2*8) on stack
    add rbp, 48
    ;# Adjust stack pointer for different alignment
    ;# Move around arguments to fit AMD64 convention below
    ;# AMD64 passes args in: rdi,rsi,rdx,rcx,r8,r9 + stack
    ;# win64 passes args in: rcx,rdx,r8,r9         + stack
    mov rdi, rcx
    mov rsi, rdx
    mov rdx, r8
    mov rcx, r9
    mov r8,  [rbp]
    mov r9,  [rbp + 8]
%endif
; .endif   # end NASM- and win64-specific block

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb334_nouter], eax
	mov [rsp + nb334_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb334_nri], edi
	mov [rsp + nb334_iinr], rsi
	mov [rsp + nb334_jindex], rdx
	mov [rsp + nb334_jjnr], rcx
	mov [rsp + nb334_shift], r8
	mov [rsp + nb334_shiftvec], r9
	mov rsi, [rbp + nb334_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb334_facel], xmm0

	mov rax, [rbp + nb334_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb334_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb334_half], eax
	movss xmm1, [rsp + nb334_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb334_half],  xmm1
	movaps [rsp + nb334_two],  xmm2
	movaps [rsp + nb334_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb334_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]		;# ebx =ii 

	mov   rdx, [rbp + nb334_charge]
	movss xmm5, [rdx + rbx*4 + 4]	
	movss xmm3, [rdx + rbx*4 + 12]	
	movss xmm4, xmm3	
	mov rsi, [rbp + nb334_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb334_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb334_qqMM], xmm3
	movaps [rsp + nb334_qqMH], xmm4
	movaps [rsp + nb334_qqHH], xmm5
	
	xorps xmm0, xmm0
	mov   rdx, [rbp + nb334_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb334_p_ntype]
	imul  ecx, [rdi]  	;# rcx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   rax, [rbp + nb334_vdwparam]
	movlps xmm0, [rax + rdx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0x55
	movaps [rsp + nb334_c6], xmm0
	movaps [rsp + nb334_c12], xmm1

.nb334_threadloop:
        mov   rsi, [rbp + nb334_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb334_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb334_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb334_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb334_n], eax
        mov [rsp + nb334_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb334_outerstart
        jmp .nb334_end
	
.nb334_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb334_nouter]
	mov [rsp + nb334_nouter], ebx

.nb334_outer:
	mov   rax, [rsp + nb334_shift]  	;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb334_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb334_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb334_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb334_pos]	;# rax = base of pos[]  
	mov   [rsp + nb334_ii3], ebx	
	
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
	movaps [rsp + nb334_ixO], xmm3
	movaps [rsp + nb334_iyO], xmm4
	movaps [rsp + nb334_izO], xmm5
	movaps [rsp + nb334_ixH1], xmm6
	movaps [rsp + nb334_iyH1], xmm7

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
	movaps [rsp + nb334_izH1], xmm6
	movaps [rsp + nb334_ixH2], xmm0
	movaps [rsp + nb334_iyH2], xmm1
	movaps [rsp + nb334_izH2], xmm2
	movaps [rsp + nb334_ixM], xmm3
	movaps [rsp + nb334_iyM], xmm4
	movaps [rsp + nb334_izM], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb334_vctot], xmm4
	movaps [rsp + nb334_Vvdwtot], xmm4
	movaps [rsp + nb334_fixO], xmm4
	movaps [rsp + nb334_fiyO], xmm4
	movaps [rsp + nb334_fizO], xmm4
	movaps [rsp + nb334_fixH1], xmm4
	movaps [rsp + nb334_fiyH1], xmm4
	movaps [rsp + nb334_fizH1], xmm4
	movaps [rsp + nb334_fixH2], xmm4
	movaps [rsp + nb334_fiyH2], xmm4
	movaps [rsp + nb334_fizH2], xmm4
	movaps [rsp + nb334_fixM], xmm4
	movaps [rsp + nb334_fiyM], xmm4
	movaps [rsp + nb334_fizM], xmm4
	
	mov   rax, [rsp + nb334_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb334_pos]
	mov   rdi, [rbp + nb334_faction]	
	mov   rax, [rsp + nb334_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb334_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb334_ninner]
	mov   [rsp + nb334_ninner], ecx
	add   edx, 0
	mov   [rsp + nb334_innerk], edx	;# number of innerloop atoms 
	jge   .nb334_unroll_loop
	jmp   .nb334_single_check
.nb334_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb334_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb334_innerjjnr],  16 ;# advance pointer (unroll 4) 

	mov rsi, [rbp + nb334_pos]   	;# base of pos[] 

	lea   rax, [rax + rax*2] 	;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2] 	;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	
	
	;# load j O coordinates 
    movlps xmm0, [rsi + rax*4] ;# jxOa jyOa  -   -
    movlps xmm1, [rsi + rcx*4] ;# jxOc jyOc  -   -
    movhps xmm0, [rsi + rbx*4] ;# jxOa jyOa jxOb jyOb 
    movhps xmm1, [rsi + rdx*4] ;# jxOc jyOc jxOd jyOd 

    movss  xmm2, [rsi + rax*4 + 8] ;# jzOa  -  -  -
    movss  xmm3, [rsi + rcx*4 + 8] ;# jzOc  -  -  -
    movhps xmm2, [rsi + rbx*4 + 8] ;# jzOa  -  jzOb  -
    movhps xmm3, [rsi + rdx*4 + 8] ;# jzOc  -  jzOd -
    
    movaps xmm4, xmm0
    unpcklps xmm0, xmm1  ;# jxOa jxOc jyOa jyOc        
    unpckhps xmm4, xmm1  ;# jxOb jxOd jyOb jyOd
    movaps xmm1, xmm0
    unpcklps xmm0, xmm4 ;# x
    unpckhps xmm1, xmm4 ;# y

    shufps xmm2, xmm3,  136  ;# 10001000 => jzOa jzOb jzOc jzOd

    ;# xmm0 = Ox
    ;# xmm1 = Oy
    ;# xmm2 = Oz
        
    subps xmm0, [rsp + nb334_ixO]
    subps xmm1, [rsp + nb334_iyO]
    subps xmm2, [rsp + nb334_izO]

    ;# store dx/dy/dz
    movaps xmm13, xmm0
    movaps xmm14, xmm1
    movaps xmm15, xmm2
    
    ;# square it
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
    
   	addps  xmm1, xmm0
	addps  xmm1, xmm2
    ;# rsq in xmm1

    ;# calculate rinv=1/sqrt(rsq)
	rsqrtps xmm5, xmm1
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm4, [rsp + nb334_three]
	mulps xmm5, xmm1	;# rsq*lu*lu 	
    subps xmm4, xmm5	;# 30-rsq*lu*lu 
	mulps xmm4, xmm2	
	mulps xmm4, [rsp + nb334_half]	
	movaps xmm2, xmm4
	mulps  xmm1, xmm4	
    ;# xmm2=rinv
    ;# xmm1=r

    mulps xmm1, [rsp + nb334_tsc] ;# rtab

    ;# truncate and convert to integers
    cvttps2dq xmm5, xmm1
    
    ;# convert back to float
    cvtdq2ps  xmm4, xmm5
    
    ;# multiply by 4
    pslld   xmm5, 2

    ;# multiply by three (copy, mult. by two, add back)
    movaps  xmm6, xmm5
    pslld   xmm5, 1
    paddd   xmm5, xmm6

    ;# calculate eps
    subps     xmm1, xmm4

    ;# move to integer registers
    movhlps xmm6, xmm5
    movd    r8d, xmm5
    movd    r10d, xmm6
    pshufd  xmm5, xmm5, 1
    pshufd  xmm6, xmm6, 1
    movd    r9d, xmm5
    movd    r11d, xmm6
    ;# table indices in r8-r11
    
    ;# xmm1=eps
    ;# xmm2=rinv
    
	mov rsi, [rbp + nb334_VFtab]
    ;# load LJ dispersion and repulsion in parallel 
    ;# NB: We are using a combined (LJ+coul) table, 
    ;# so the LJ table data is offset 4*4 = 16 bytes
    movlps xmm5, [rsi + r8*4 + 16]
   	movlps xmm9, [rsi + r8*4 + 32]
	movlps xmm7,  [rsi + r10*4 + 16]
	movlps xmm11, [rsi + r10*4 + 32]
	movhps xmm5, [rsi + r9*4 + 16]
	movhps xmm9, [rsi + r9*4 + 32]
	movhps xmm7,  [rsi + r11*4 + 16]
	movhps xmm11, [rsi + r11*4 + 32]

    movaps xmm4, xmm5
    movaps xmm8, xmm9
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm8, xmm11, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101
	shufps xmm9, xmm11, 221  ;# 11011101

	movlps xmm7,  [rsi + r8*4 + 24]
	movlps xmm11, [rsi + r8*4 + 40]
	movlps xmm0,  [rsi + r10*4 + 24]
	movlps xmm3,  [rsi + r10*4 + 40]
	movhps xmm7,  [rsi + r9*4 + 24]
	movhps xmm11, [rsi + r9*4 + 40]
	movhps xmm0,  [rsi + r11*4 + 24]
	movhps xmm3,  [rsi + r11*4 + 40]

    movaps xmm6, xmm7
    movaps xmm10, xmm11
    
	shufps xmm6, xmm0, 136  ;# 10001000
	shufps xmm10, xmm3, 136  ;# 10001000
	shufps xmm7, xmm0, 221  ;# 11011101
	shufps xmm11, xmm3, 221  ;# 11011101
    ;# dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11
    
    mulps  xmm7, xmm1    ;# Heps
    mulps  xmm11, xmm1 
    mulps  xmm6, xmm1   ;# Geps
    mulps  xmm10, xmm1 
    mulps  xmm7, xmm1   ;# Heps2
    mulps  xmm11, xmm1 
    addps  xmm5, xmm6  ;# F+Geps
    addps  xmm9, xmm10 
    addps  xmm5, xmm7   ;# F+Geps+Heps2 = Fp
    addps  xmm9, xmm11 
    addps  xmm7, xmm7    ;# 2*Heps2
    addps  xmm11, xmm11
    addps  xmm7, xmm6   ;# 2*Heps2+Geps
    addps  xmm11, xmm10
    
    addps  xmm7, xmm5  ;# FF = Fp + 2*Heps2 + Geps
    addps  xmm11, xmm9
    mulps  xmm5, xmm1  ;# eps*Fp
    mulps  xmm9, xmm1
    addps  xmm5, xmm4 ;# VV
    addps  xmm9, xmm8

    mulps  xmm5, [rsp + nb334_c6]  ;# VV*c6 = vnb6
    mulps  xmm9, [rsp + nb334_c12]  ;# VV*c12 = vnb12
    addps  xmm5, xmm9
    addps  xmm5, [rsp + nb334_Vvdwtot]
    movaps [rsp + nb334_Vvdwtot], xmm5
        
    mulps  xmm7, [rsp + nb334_c6]   ;# FF*c6 = fnb6
    mulps  xmm11, [rsp + nb334_c12]   ;# FF*c12  = fnb12
    addps  xmm7, xmm11
    
    mulps  xmm7, [rsp + nb334_tsc]
    mulps  xmm7, xmm2
    xorps  xmm9, xmm9
    
    subps  xmm9, xmm7

    ;# fx/fy/fz
    mulps  xmm13, xmm9
    mulps  xmm14, xmm9
    mulps  xmm15, xmm9

    ;# increment i force
    movaps xmm0, [rsp + nb334_fixO]
    movaps xmm1, [rsp + nb334_fiyO]
    movaps xmm2, [rsp + nb334_fizO]
    addps  xmm0, xmm13
    addps  xmm1, xmm14
    addps  xmm2, xmm15
    movaps [rsp + nb334_fixO], xmm0
    movaps [rsp + nb334_fiyO], xmm1
    movaps [rsp + nb334_fizO], xmm2

	;# move j O forces to local temp variables 
    mov rdi, [rbp + nb334_faction]
    movlps xmm4, [rdi + rax*4] ;# jxOa jyOa  -   -
    movlps xmm5, [rdi + rcx*4] ;# jxOc jyOc  -   -
    movhps xmm4, [rdi + rbx*4] ;# jxOa jyOa jxOb jyOb 
    movhps xmm5, [rdi + rdx*4] ;# jxOc jyOc jxOd jyOd 

    movss  xmm6, [rdi + rax*4 + 8] ;# jzOa  -  -  -
    movss  xmm7, [rdi + rcx*4 + 8] ;# jzOc  -  -  -
    movhps xmm6, [rdi + rbx*4 + 8] ;# jzOa  -  jzOb  -
    movhps xmm7, [rdi + rdx*4 + 8] ;# jzOc  -  jzOd -

    shufps xmm6, xmm7,  136  ;# 10001000 => jzOa jzOb jzOc jzOd

    ;# xmm4: jxOa jyOa jxOb jyOb 
    ;# xmm5: jxOc jyOc jxOd jyOd
    ;# xmm6: jzOa jzOb jzOc jzOd

    ;# update O forces
    movaps xmm12, xmm13
    unpcklps xmm13, xmm14   ;# (local) fjx1 fjx1 fjy1 fjy2
    unpckhps xmm12, xmm14   ;# (local) fjx3 fjx4 fjy3 fjy4
    
    addps xmm4, xmm13
    addps xmm5, xmm12
    addps xmm6, xmm15
    
    movhlps  xmm7, xmm6 ;# fH1zc fH1zd

    movlps [rdi + rax*4], xmm4
    movhps [rdi + rbx*4], xmm4
    movlps [rdi + rcx*4], xmm5
    movhps [rdi + rdx*4], xmm5
    movss  [rdi + rax*4 + 8], xmm6
    movss  [rdi + rcx*4 + 8], xmm7
    shufps xmm6, xmm6, 1
    shufps xmm7, xmm7, 1
    movss  [rdi + rbx*4 + 8], xmm6
    movss  [rdi + rdx*4 + 8], xmm7
    ;# done with OO interaction
    
    ;# move j H1 coordinates to local temp variables 
    mov rsi, [rbp + nb334_pos]
    movlps xmm0, [rsi + rax*4 + 12] ;# jxH1a jyH1a  -   -
    movlps xmm1, [rsi + rcx*4 + 12] ;# jxH1c jyH1c  -   -
    movhps xmm0, [rsi + rbx*4 + 12] ;# jxH1a jyH1a jxH1b jyH1b 
    movhps xmm1, [rsi + rdx*4 + 12] ;# jxH1c jyH1c jxH1d jyH1d 

    movss  xmm2, [rsi + rax*4 + 20] ;# jzH1a  -  -  -
    movss  xmm3, [rsi + rcx*4 + 20] ;# jzH1c  -  -  -
    movhps xmm2, [rsi + rbx*4 + 20] ;# jzH1a  -  jzH1b  -
    movhps xmm3, [rsi + rdx*4 + 20] ;# jzH1c  -  jzH1d -

    movd mm0, eax ;# save j3 in mm0-mm3
    movd mm1, ebx
    movd mm2, ecx
    movd mm3, edx
    
    movaps xmm4, xmm0
    unpcklps xmm0, xmm1  ;# jxH1a jxH1c jyH1a jyH1c        
    unpckhps xmm4, xmm1  ;# jxH1b jxH1d jyH1b jyH1d
    movaps xmm1, xmm0
    unpcklps xmm0, xmm4 ;# x
    unpckhps xmm1, xmm4 ;# y

    shufps   xmm2, xmm3,  136  ;# 10001000 => jzH1a jzH1b jzH1c jzH1d

    ;# xmm0 = H1x
    ;# xmm1 = H1y
    ;# xmm2 = H1z
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + nb334_ixH1]
    subps xmm1, [rsp + nb334_iyH1]
    subps xmm2, [rsp + nb334_izH1]
    subps xmm3, [rsp + nb334_ixH2]
    subps xmm4, [rsp + nb334_iyH2]
    subps xmm5, [rsp + nb334_izH2]
    subps xmm6, [rsp + nb334_ixM]
    subps xmm7, [rsp + nb334_iyM]
    subps xmm8, [rsp + nb334_izM]
    
	movaps [rsp + nb334_dxH1H1], xmm0
	movaps [rsp + nb334_dyH1H1], xmm1
	movaps [rsp + nb334_dzH1H1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb334_dxH2H1], xmm3
	movaps [rsp + nb334_dyH2H1], xmm4
	movaps [rsp + nb334_dzH2H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb334_dxMH1], xmm6
	movaps [rsp + nb334_dyMH1], xmm7
	movaps [rsp + nb334_dzMH1], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt for jH1 atoms
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + nb334_three]
	movaps  xmm10, xmm9
    movaps  xmm11, xmm9

	mulps   xmm1, xmm0 ;# rsq*lu*lu
	mulps   xmm4, xmm3 ;# rsq*lu*lu 
    mulps   xmm7, xmm6 ;# rsq*lu*lu
	
	subps   xmm9, xmm1
	subps   xmm10, xmm4
    subps   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulps   xmm9, xmm2
	mulps   xmm10, xmm5
    mulps   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movaps  xmm4, [rsp + nb334_half]
	mulps   xmm9, xmm4  ;# rinvH1H1 
	mulps   xmm10, xmm4 ;# rinvH2H1
    mulps   xmm11, xmm4 ;# rinvMH1

	movaps  [rsp + nb334_rinvH1H1], xmm9
	movaps  [rsp + nb334_rinvH2H1], xmm10
	movaps  [rsp + nb334_rinvMH1], xmm11
	
	;# H1 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movaps xmm1, [rsp + nb334_tsc]
    mulps  xmm0, xmm9  ;# r
    mulps  xmm3, xmm10
    mulps  xmm6, xmm11
    mulps  xmm0, xmm1 ;# rtab
    mulps  xmm3, xmm1
    mulps  xmm6, xmm1
    
    ;# truncate and convert to integers
    cvttps2dq xmm1, xmm0
    cvttps2dq xmm4, xmm3
    cvttps2dq xmm7, xmm6        

    ;# convert back to float
    cvtdq2ps  xmm2, xmm1
    cvtdq2ps  xmm5, xmm4
    cvtdq2ps  xmm8, xmm7
    
    ;# multiply by 4
    pslld   xmm1, 2
    pslld   xmm4, 2
    pslld   xmm7, 2
    
    ;# multiply by three (copy, mult. by two, add back)
    movaps  xmm10, xmm1
    movaps  xmm11, xmm4
    movaps  xmm12, xmm7
    pslld   xmm1, 1
    pslld   xmm4, 1
    pslld   xmm7, 1    
    paddd   xmm1, xmm10
    paddd   xmm4, xmm11
    paddd   xmm7, xmm12    
    
    ;# move to integer registers
    movhlps xmm13, xmm1
    movhlps xmm14, xmm4
    movhlps xmm15, xmm7
    movd    eax, xmm1
    movd    r8d, xmm4
    movd    r12d, xmm7
    movd    ecx, xmm13
    movd    r10d, xmm14
    movd    r14d, xmm15
    pshufd  xmm1, xmm1, 1
    pshufd  xmm4, xmm4, 1
    pshufd  xmm7, xmm7, 1
    pshufd  xmm13, xmm13, 1
    pshufd  xmm14, xmm14, 1
    pshufd  xmm15, xmm15, 1
    movd    ebx, xmm1
    movd    r9d, xmm4
    movd    r13d, xmm7    
    movd    edx, xmm13
    movd    r11d, xmm14
    movd    r15d, xmm15   
        
    mov  rsi, [rbp + nb334_VFtab]

    ;# calculate eps
    subps     xmm0, xmm2
    subps     xmm3, xmm5
    subps     xmm6, xmm8

    movaps    [rsp + nb334_epsH1], xmm0
    movaps    [rsp + nb334_epsH2], xmm3
    movaps    [rsp + nb334_epsM], xmm6

    ;# Load LOTS of table data
   	movlps xmm1, [rsi + rax*4]
   	movlps xmm5, [rsi + r8*4]
   	movlps xmm9, [rsi + r12*4]

	movlps xmm3, [rsi + rcx*4]
	movlps xmm7, [rsi + r10*4]
	movlps xmm11, [rsi + r14*4]

	movhps xmm1, [rsi + rbx*4]
	movhps xmm5, [rsi + r9*4]
	movhps xmm9, [rsi + r13*4]

	movhps xmm3, [rsi + rdx*4]
	movhps xmm7, [rsi + r11*4]
	movhps xmm11, [rsi + r15*4]

    movaps xmm0, xmm1
    movaps xmm4, xmm5
    movaps xmm8, xmm9
	shufps xmm0, xmm3, 136  ;# 10001000
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm8, xmm11, 136  ;# 10001000
	shufps xmm1, xmm3, 221  ;# 11011101
	shufps xmm5, xmm7, 221  ;# 11011101
	shufps xmm9, xmm11, 221  ;# 11011101
    
	movlps xmm3, [rsi + rax*4 + 8]
	movlps xmm7, [rsi + r8*4 + 8]
	movlps xmm11, [rsi + r12*4 + 8]
    
	movlps xmm12, [rsi + rcx*4 + 8]
	movlps xmm13, [rsi + r10*4 + 8]
	movlps xmm14, [rsi + r14*4 + 8]

	movhps xmm3, [rsi + rbx*4 + 8]
	movhps xmm7, [rsi + r9*4 + 8]
	movhps xmm11, [rsi + r13*4 + 8]
    
	movhps xmm12, [rsi + rdx*4 + 8]
	movhps xmm13, [rsi + r11*4 + 8]
	movhps xmm14, [rsi + r15*4 + 8]

    movaps xmm2, xmm3
    movaps xmm6, xmm7
    movaps xmm10, xmm11
    
	shufps xmm2, xmm12, 136  ;# 10001000
	shufps xmm6, xmm13, 136  ;# 10001000
	shufps xmm10, xmm14, 136  ;# 10001000
	shufps xmm3, xmm12, 221  ;# 11011101
	shufps xmm7, xmm13, 221  ;# 11011101
	shufps xmm11, xmm14, 221  ;# 11011101
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    movaps xmm12, [rsp + nb334_epsH1]
    movaps xmm13, [rsp + nb334_epsH2]
    movaps xmm14, [rsp + nb334_epsM]
    
    mulps  xmm3, xmm12   ;# Heps
    mulps  xmm7, xmm13
    mulps  xmm11, xmm14 
    mulps  xmm2, xmm12   ;# Geps
    mulps  xmm6, xmm13
    mulps  xmm10, xmm14 
    mulps  xmm3, xmm12   ;# Heps2
    mulps  xmm7, xmm13
    mulps  xmm11, xmm14 

    addps  xmm1, xmm2   ;# F+Geps
    addps  xmm5, xmm6
    addps  xmm9, xmm10 
    addps  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addps  xmm5, xmm7
    addps  xmm9, xmm11 
    addps  xmm3, xmm3    ;# 2*Heps2
    addps  xmm7, xmm7
    addps  xmm11, xmm11
    addps  xmm3, xmm2    ;# 2*Heps2+Geps
    addps  xmm7, xmm6  
    addps  xmm11, xmm10
    addps  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addps  xmm7, xmm5
    addps  xmm11, xmm9
    mulps  xmm1, xmm12   ;# eps*Fp
    mulps  xmm5, xmm13
    mulps  xmm9, xmm14
    movaps xmm12, [rsp + nb334_qqHH]
    movaps xmm13, [rsp + nb334_qqMH]
    addps  xmm1, xmm0     ;# VV
    addps  xmm5, xmm4
    addps  xmm9, xmm8
    mulps  xmm1, xmm12   ;# VV*qq = vcoul
    mulps  xmm5, xmm12
    mulps  xmm9, xmm13
    mulps  xmm3, xmm12    ;# FF*qq = fij
    mulps  xmm7, xmm12
    mulps  xmm11, xmm13
    
    ;# accumulate vctot
    addps  xmm1, [rsp + nb334_vctot]
    addps  xmm5, xmm9
    addps  xmm1, xmm5
    movaps [rsp + nb334_vctot], xmm1
    
    movaps xmm10, [rsp + nb334_tsc]
    mulps  xmm3, xmm10  ;# fscal
    mulps  xmm7, xmm10
    mulps  xmm10, xmm11
    
    movd eax, mm0 ;# restore j3 from mm0-mm3
    movd ebx, mm1
    movd ecx, mm2
    movd edx, mm3

	;# move j H1 forces to local temp variables 
    mov rdi, [rbp + nb334_faction]
    movlps xmm11, [rdi + rax*4 + 12] ;# jxH1a jyH1a  -   -
    movlps xmm12, [rdi + rcx*4 + 12] ;# jxH1c jyH1c  -   -
    movhps xmm11, [rdi + rbx*4 + 12] ;# jxH1a jyH1a jxH1b jyH1b 
    movhps xmm12, [rdi + rdx*4 + 12] ;# jxH1c jyH1c jxH1d jyH1d 

    movss  xmm13, [rdi + rax*4 + 20] ;# jzH1a  -  -  -
    movss  xmm14, [rdi + rcx*4 + 20] ;# jzH1c  -  -  -
    movhps xmm13, [rdi + rbx*4 + 20] ;# jzH1a  -  jzH1b  -
    movhps xmm14, [rdi + rdx*4 + 20] ;# jzH1c  -  jzH1d -
    
    shufps xmm13, xmm14,  136  ;# 10001000 => jzH1a jzH1b jzH1c jzH1d

    ;# xmm11: jxH1a jyH1a jxH1b jyH1b 
    ;# xmm12: jxH1c jyH1c jxH1d jyH1d
    ;# xmm13: jzH1a jzH1b jzH1c jzH1d

    xorps  xmm0, xmm0
    xorps  xmm4, xmm4
    xorps  xmm8, xmm8

    mulps  xmm3, [rsp + nb334_rinvH1H1]
    mulps  xmm7, [rsp + nb334_rinvH2H1]
    mulps  xmm10, [rsp + nb334_rinvMH1]
    
    subps  xmm0, xmm3
    subps  xmm4, xmm7
    subps  xmm8, xmm10
    
    movaps xmm1, xmm0
    movaps xmm2, xmm0
    movaps xmm3, xmm4
    movaps xmm5, xmm4
    movaps xmm6, xmm8
    movaps xmm7, xmm8

	mulps xmm0, [rsp + nb334_dxH1H1]
	mulps xmm1, [rsp + nb334_dyH1H1]
	mulps xmm2, [rsp + nb334_dzH1H1]
	mulps xmm3, [rsp + nb334_dxH2H1]
	mulps xmm4, [rsp + nb334_dyH2H1]
	mulps xmm5, [rsp + nb334_dzH2H1]
	mulps xmm6, [rsp + nb334_dxMH1]
	mulps xmm7, [rsp + nb334_dyMH1]
	mulps xmm8, [rsp + nb334_dzMH1]

    movaps xmm14, xmm0
    movaps xmm15, xmm1
    addps xmm13,  xmm2
    addps xmm0, [rsp + nb334_fixH1]
    addps xmm1, [rsp + nb334_fiyH1]
    addps xmm2, [rsp + nb334_fizH1]

    addps xmm14, xmm3
    addps xmm15, xmm4
    addps xmm13, xmm5
    addps xmm3, [rsp + nb334_fixH2]
    addps xmm4, [rsp + nb334_fiyH2]
    addps xmm5, [rsp + nb334_fizH2]

    addps xmm14, xmm6
    addps xmm15, xmm7
    addps xmm13, xmm8
    addps xmm6, [rsp + nb334_fixM]
    addps xmm7, [rsp + nb334_fiyM]
    addps xmm8, [rsp + nb334_fizM]

    movaps [rsp + nb334_fixH1], xmm0
    movaps [rsp + nb334_fiyH1], xmm1
    movaps [rsp + nb334_fizH1], xmm2
    movaps [rsp + nb334_fixH2], xmm3
    movaps [rsp + nb334_fiyH2], xmm4
    movaps [rsp + nb334_fizH2], xmm5
    movaps [rsp + nb334_fixM], xmm6
    movaps [rsp + nb334_fiyM], xmm7
    movaps [rsp + nb334_fizM], xmm8
    
    ;# xmm14 = fH1x
    ;# xmm15 = fH1y
    ;# xmm13 = fH1z
    movaps xmm0, xmm14
    unpcklps xmm14, xmm15
    unpckhps xmm0,  xmm15
    
    addps  xmm11, xmm14
    addps  xmm12, xmm0
    
    movhlps  xmm14, xmm13 ;# fH1zc fH1zd
    
    movlps [rdi + rax*4 + 12], xmm11
    movhps [rdi + rbx*4 + 12], xmm11
    movlps [rdi + rcx*4 + 12], xmm12
    movhps [rdi + rdx*4 + 12], xmm12
    movss  [rdi + rax*4 + 20], xmm13
    movss  [rdi + rcx*4 + 20], xmm14
    shufps xmm13, xmm13, 1
    shufps xmm14, xmm14, 1
    movss  [rdi + rbx*4 + 20], xmm13
    movss  [rdi + rdx*4 + 20], xmm14
    
	;# move j H2 coordinates to local temp variables 
	mov   rsi, [rbp + nb334_pos]
    movlps xmm0, [rsi + rax*4 + 24] ;# jxH2a jyH2a  -   -
    movlps xmm1, [rsi + rcx*4 + 24] ;# jxH2c jyH2c  -   -
    movhps xmm0, [rsi + rbx*4 + 24] ;# jxH2a jyH2a jxH2b jyH2b 
    movhps xmm1, [rsi + rdx*4 + 24] ;# jxH2c jyH2c jxH2d jyH2d 

    movss  xmm2, [rsi + rax*4 + 32] ;# jzH2a  -  -  -
    movss  xmm3, [rsi + rcx*4 + 32] ;# jzH2c  -  -  -
    movhps xmm2, [rsi + rbx*4 + 32] ;# jzH2a  -  jzH2b  -
    movhps xmm3, [rsi + rdx*4 + 32] ;# jzH2c  -  jzH2d -

    movaps xmm4, xmm0
    unpcklps xmm0, xmm1  ;# jxH2a jxH2c jyH2a jyH2c        
    unpckhps xmm4, xmm1  ;# jxH2b jxH2d jyH2b jyH2d
    movaps xmm1, xmm0
    unpcklps xmm0, xmm4 ;# x
    unpckhps xmm1, xmm4 ;# y

    shufps   xmm2, xmm3,  136  ;# 10001000 => jzH2a jzH2b jzH2c jzH2d

    ;# xmm0 = H2x
    ;# xmm1 = H2y
    ;# xmm2 = H2z
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + nb334_ixH1]
    subps xmm1, [rsp + nb334_iyH1]
    subps xmm2, [rsp + nb334_izH1]
    subps xmm3, [rsp + nb334_ixH2]
    subps xmm4, [rsp + nb334_iyH2]
    subps xmm5, [rsp + nb334_izH2]
    subps xmm6, [rsp + nb334_ixM]
    subps xmm7, [rsp + nb334_iyM]
    subps xmm8, [rsp + nb334_izM]
    
	movaps [rsp + nb334_dxH1H2], xmm0
	movaps [rsp + nb334_dyH1H2], xmm1
	movaps [rsp + nb334_dzH1H2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb334_dxH2H2], xmm3
	movaps [rsp + nb334_dyH2H2], xmm4
	movaps [rsp + nb334_dzH2H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb334_dxMH2], xmm6
	movaps [rsp + nb334_dyMH2], xmm7
	movaps [rsp + nb334_dzMH2], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt for jH2 atoms
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + nb334_three]
	movaps  xmm10, xmm9
    movaps  xmm11, xmm9

	mulps   xmm1, xmm0 ;# rsq*lu*lu
	mulps   xmm4, xmm3 ;# rsq*lu*lu 
    mulps   xmm7, xmm6 ;# rsq*lu*lu
	
	subps   xmm9, xmm1
	subps   xmm10, xmm4
    subps   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulps   xmm9, xmm2
	mulps   xmm10, xmm5
    mulps   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movaps  xmm4, [rsp + nb334_half]
	mulps   xmm9, xmm4  ;# rinvH1H2
	mulps   xmm10, xmm4 ;# rinvH2H2
    mulps   xmm11, xmm4 ;# rinvMH2

	movaps  [rsp + nb334_rinvH1H2], xmm9
	movaps  [rsp + nb334_rinvH2H2], xmm10
	movaps  [rsp + nb334_rinvMH2], xmm11
	
	;# H2 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movaps xmm1, [rsp + nb334_tsc]
    mulps  xmm0, xmm9  ;# r
    mulps  xmm3, xmm10
    mulps  xmm6, xmm11
    mulps  xmm0, xmm1 ;# rtab
    mulps  xmm3, xmm1
    mulps  xmm6, xmm1

    mov  rsi, [rbp + nb334_VFtab]
    
    ;# truncate and convert to integers
    cvttps2dq xmm1, xmm0
    cvttps2dq xmm4, xmm3
    cvttps2dq xmm7, xmm6        
    
    ;# convert back to float
    cvtdq2ps  xmm2, xmm1
    cvtdq2ps  xmm5, xmm4
    cvtdq2ps  xmm8, xmm7
         
    ;# multiply by 4
    pslld   xmm1, 2
    pslld   xmm4, 2
    pslld   xmm7, 2
    
    ;# multiply by three (copy, mult. by two, add back)
    movaps  xmm10, xmm1
    movaps  xmm11, xmm4
    movaps  xmm12, xmm7
    pslld   xmm1, 1
    pslld   xmm4, 1
    pslld   xmm7, 1    
    paddd   xmm1, xmm10
    paddd   xmm4, xmm11
    paddd   xmm7, xmm12    
    
    ;# move to integer registers
    movhlps xmm13, xmm1
    movhlps xmm14, xmm4
    movhlps xmm15, xmm7
    movd    eax, xmm1
    movd    r8d, xmm4
    movd    r12d, xmm7
    movd    ecx, xmm13
    movd    r10d, xmm14
    movd    r14d, xmm15
    pshufd  xmm1, xmm1, 1
    pshufd  xmm4, xmm4, 1
    pshufd  xmm7, xmm7, 1
    pshufd  xmm13, xmm13, 1
    pshufd  xmm14, xmm14, 1
    pshufd  xmm15, xmm15, 1
    movd    ebx, xmm1
    movd    r9d, xmm4
    movd    r13d, xmm7    
    movd    edx, xmm13
    movd    r11d, xmm14
    movd    r15d, xmm15   

    ;# calculate eps
    subps     xmm0, xmm2
    subps     xmm3, xmm5
    subps     xmm6, xmm8

    movaps    [rsp + nb334_epsH1], xmm0
    movaps    [rsp + nb334_epsH2], xmm3
    movaps    [rsp + nb334_epsM], xmm6


    ;# Load LOTS of table data
   	movlps xmm1, [rsi + rax*4]
   	movlps xmm5, [rsi + r8*4]
   	movlps xmm9, [rsi + r12*4]

	movlps xmm3, [rsi + rcx*4]
	movlps xmm7, [rsi + r10*4]
	movlps xmm11, [rsi + r14*4]

	movhps xmm1, [rsi + rbx*4]
	movhps xmm5, [rsi + r9*4]
	movhps xmm9, [rsi + r13*4]

	movhps xmm3, [rsi + rdx*4]
	movhps xmm7, [rsi + r11*4]
	movhps xmm11, [rsi + r15*4]

    movaps xmm0, xmm1
    movaps xmm4, xmm5
    movaps xmm8, xmm9
	shufps xmm0, xmm3, 136  ;# 10001000
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm8, xmm11, 136  ;# 10001000
	shufps xmm1, xmm3, 221  ;# 11011101
	shufps xmm5, xmm7, 221  ;# 11011101
	shufps xmm9, xmm11, 221  ;# 11011101
    
	movlps xmm3, [rsi + rax*4 + 8]
	movlps xmm7, [rsi + r8*4 + 8]
	movlps xmm11, [rsi + r12*4 + 8]
    
	movlps xmm12, [rsi + rcx*4 + 8]
	movlps xmm13, [rsi + r10*4 + 8]
	movlps xmm14, [rsi + r14*4 + 8]

	movhps xmm3, [rsi + rbx*4 + 8]
	movhps xmm7, [rsi + r9*4 + 8]
	movhps xmm11, [rsi + r13*4 + 8]
    
	movhps xmm12, [rsi + rdx*4 + 8]
	movhps xmm13, [rsi + r11*4 + 8]
	movhps xmm14, [rsi + r15*4 + 8]

    movaps xmm2, xmm3
    movaps xmm6, xmm7
    movaps xmm10, xmm11
    
	shufps xmm2, xmm12, 136  ;# 10001000
	shufps xmm6, xmm13, 136  ;# 10001000
	shufps xmm10, xmm14, 136  ;# 10001000
	shufps xmm3, xmm12, 221  ;# 11011101
	shufps xmm7, xmm13, 221  ;# 11011101
	shufps xmm11, xmm14, 221  ;# 11011101
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    movaps xmm12, [rsp + nb334_epsH1]
    movaps xmm13, [rsp + nb334_epsH2]
    movaps xmm14, [rsp + nb334_epsM]
    
    mulps  xmm3, xmm12   ;# Heps
    mulps  xmm7, xmm13
    mulps  xmm11, xmm14 
    mulps  xmm2, xmm12   ;# Geps
    mulps  xmm6, xmm13
    mulps  xmm10, xmm14 
    mulps  xmm3, xmm12   ;# Heps2
    mulps  xmm7, xmm13
    mulps  xmm11, xmm14 

    addps  xmm1, xmm2   ;# F+Geps
    addps  xmm5, xmm6
    addps  xmm9, xmm10 
    addps  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addps  xmm5, xmm7
    addps  xmm9, xmm11 
    addps  xmm3, xmm3    ;# 2*Heps2
    addps  xmm7, xmm7
    addps  xmm11, xmm11
    addps  xmm3, xmm2    ;# 2*Heps2+Geps
    addps  xmm7, xmm6  
    addps  xmm11, xmm10
    addps  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addps  xmm7, xmm5
    addps  xmm11, xmm9
    mulps  xmm1, xmm12   ;# eps*Fp
    mulps  xmm5, xmm13
    mulps  xmm9, xmm14
    movaps xmm12, [rsp + nb334_qqHH]
    movaps xmm13, [rsp + nb334_qqMH]
    addps  xmm1, xmm0     ;# VV
    addps  xmm5, xmm4
    addps  xmm9, xmm8
    mulps  xmm1, xmm12   ;# VV*qq = vcoul
    mulps  xmm5, xmm12
    mulps  xmm9, xmm13
    mulps  xmm3, xmm12    ;# FF*qq = fij
    mulps  xmm7, xmm12
    mulps  xmm11, xmm13
    
    ;# accumulate vctot
    addps  xmm1, [rsp + nb334_vctot]
    addps  xmm5, xmm9
    addps  xmm1, xmm5
    movaps [rsp + nb334_vctot], xmm1
    
    movaps xmm10, [rsp + nb334_tsc]
    mulps  xmm3, xmm10  ;# fscal
    mulps  xmm7, xmm10
    mulps  xmm10, xmm11
        
    movd eax, mm0 ;# restore j3 from mm0-mm3
    movd ebx, mm1
    movd ecx, mm2
    movd edx, mm3
    
	;# move j H2 forces to local temp variables 
    mov rdi, [rbp + nb334_faction]
    movlps xmm11, [rdi + rax*4 + 24] ;# jxH2a jyH2a  -   -
    movlps xmm12, [rdi + rcx*4 + 24] ;# jxH2c jyH2c  -   -
    movhps xmm11, [rdi + rbx*4 + 24] ;# jxH2a jyH2a jxH2b jyH2b 
    movhps xmm12, [rdi + rdx*4 + 24] ;# jxH2c jyH2c jxH2d jyH2d 

    movss  xmm13, [rdi + rax*4 + 32] ;# jzH2a  -  -  -
    movss  xmm14, [rdi + rcx*4 + 32] ;# jzH2c  -  -  -
    movhps xmm13, [rdi + rbx*4 + 32] ;# jzH2a  -  jzH2b  -
    movhps xmm14, [rdi + rdx*4 + 32] ;# jzH2c  -  jzH2d -
    
    shufps xmm13, xmm14,  136  ;# 10001000 => jzH2a jzH2b jzH2c jzH2d

    ;# xmm11: jxH2a jyH2a jxH2b jyH2b 
    ;# xmm12: jxH2c jyH2c jxH2d jyH2d
    ;# xmm13: jzH2a jzH2b jzH2c jzH2d

    xorps  xmm0, xmm0
    xorps  xmm4, xmm4
    xorps  xmm8, xmm8

    mulps  xmm3, [rsp + nb334_rinvH1H2]
    mulps  xmm7, [rsp + nb334_rinvH2H2]
    mulps  xmm10, [rsp + nb334_rinvMH2]
    
    subps  xmm0, xmm3
    subps  xmm4, xmm7
    subps  xmm8, xmm10
    
    movaps xmm1, xmm0
    movaps xmm2, xmm0
    movaps xmm3, xmm4
    movaps xmm5, xmm4
    movaps xmm6, xmm8
    movaps xmm7, xmm8

	mulps xmm0, [rsp + nb334_dxH1H2]
	mulps xmm1, [rsp + nb334_dyH1H2]
	mulps xmm2, [rsp + nb334_dzH1H2]
	mulps xmm3, [rsp + nb334_dxH2H2]
	mulps xmm4, [rsp + nb334_dyH2H2]
	mulps xmm5, [rsp + nb334_dzH2H2]
	mulps xmm6, [rsp + nb334_dxMH2]
	mulps xmm7, [rsp + nb334_dyMH2]
	mulps xmm8, [rsp + nb334_dzMH2]

    movaps xmm14, xmm0
    movaps xmm15, xmm1
    addps xmm13, xmm2
    addps xmm0, [rsp + nb334_fixH1]
    addps xmm1, [rsp + nb334_fiyH1]
    addps xmm2, [rsp + nb334_fizH1]

    addps xmm14, xmm3
    addps xmm15, xmm4
    addps xmm13, xmm5
    addps xmm3, [rsp + nb334_fixH2]
    addps xmm4, [rsp + nb334_fiyH2]
    addps xmm5, [rsp + nb334_fizH2]

    addps xmm14, xmm6
    addps xmm15, xmm7
    addps xmm13, xmm8
    addps xmm6, [rsp + nb334_fixM]
    addps xmm7, [rsp + nb334_fiyM]
    addps xmm8, [rsp + nb334_fizM]

    movaps [rsp + nb334_fixH1], xmm0
    movaps [rsp + nb334_fiyH1], xmm1
    movaps [rsp + nb334_fizH1], xmm2
    movaps [rsp + nb334_fixH2], xmm3
    movaps [rsp + nb334_fiyH2], xmm4
    movaps [rsp + nb334_fizH2], xmm5
    movaps [rsp + nb334_fixM], xmm6
    movaps [rsp + nb334_fiyM], xmm7
    movaps [rsp + nb334_fizM], xmm8
    
    ;# xmm14 = fH2x
    ;# xmm15 = fH2y
    ;# xmm13 = fH2z
    movaps xmm0, xmm14
    unpcklps xmm14, xmm15
    unpckhps xmm0, xmm15
    
    addps  xmm11, xmm14
    addps  xmm12, xmm0
    
    movhlps  xmm14, xmm13 ;# fH2zc fH2zd
    
    movlps [rdi + rax*4 + 24], xmm11
    movhps [rdi + rbx*4 + 24], xmm11
    movlps [rdi + rcx*4 + 24], xmm12
    movhps [rdi + rdx*4 + 24], xmm12
    movss  [rdi + rax*4 + 32], xmm13
    movss  [rdi + rcx*4 + 32], xmm14
    shufps xmm13, xmm13, 1
    shufps xmm14, xmm14, 1
    movss  [rdi + rbx*4 + 32], xmm13
    movss  [rdi + rdx*4 + 32], xmm14
    
   	mov   rsi, [rbp + nb334_pos]
	;# move j M coordinates to local temp variables 
    movlps xmm0, [rsi + rax*4 + 36] ;# jxMa jyMa  -   -
    movlps xmm1, [rsi + rcx*4 + 36] ;# jxMc jyMc  -   -
    movhps xmm0, [rsi + rbx*4 + 36] ;# jxMa jyMa jxMb jyMb 
    movhps xmm1, [rsi + rdx*4 + 36] ;# jxMc jyMc jxMd jyMd 

    movss  xmm2, [rsi + rax*4 + 44] ;# jzMa  -  -  -
    movss  xmm3, [rsi + rcx*4 + 44] ;# jzMc  -  -  -
    movss  xmm5, [rsi + rbx*4 + 44] ;# jzMb  -  -  -
    movss  xmm6, [rsi + rdx*4 + 44] ;# jzMd  -  -  -
    movlhps xmm2, xmm5 ;# jzMa  -  jzMb  -
    movlhps xmm3, xmm6 ;# jzMc  -  jzMd -

    movaps xmm4, xmm0
    unpcklps xmm0, xmm1  ;# jxMa jxMc jyMa jyMc        
    unpckhps xmm4, xmm1  ;# jxMb jxMd jyMb jyMd
    movaps xmm1, xmm0
    unpcklps xmm0, xmm4 ;# x
    unpckhps xmm1, xmm4 ;# y

    shufps   xmm2, xmm3,  136  ;# 10001000 => jzMa jzMb jzMc jzMd

    ;# xmm0 = Mx
    ;# xmm1 = My
    ;# xmm2 = Mz
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + nb334_ixH1]
    subps xmm1, [rsp + nb334_iyH1]
    subps xmm2, [rsp + nb334_izH1]
    subps xmm3, [rsp + nb334_ixH2]
    subps xmm4, [rsp + nb334_iyH2]
    subps xmm5, [rsp + nb334_izH2]
    subps xmm6, [rsp + nb334_ixM]
    subps xmm7, [rsp + nb334_iyM]
    subps xmm8, [rsp + nb334_izM]
    
	movaps [rsp + nb334_dxH1M], xmm0
	movaps [rsp + nb334_dyH1M], xmm1
	movaps [rsp + nb334_dzH1M], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb334_dxH2M], xmm3
	movaps [rsp + nb334_dyH2M], xmm4
	movaps [rsp + nb334_dzH2M], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb334_dxMM], xmm6
	movaps [rsp + nb334_dyMM], xmm7
	movaps [rsp + nb334_dzMM], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt for jM atoms
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + nb334_three]
	movaps  xmm10, xmm9
    movaps  xmm11, xmm9

	mulps   xmm1, xmm0 ;# rsq*lu*lu
	mulps   xmm4, xmm3 ;# rsq*lu*lu 
    mulps   xmm7, xmm6 ;# rsq*lu*lu
	
	subps   xmm9, xmm1
	subps   xmm10, xmm4
    subps   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulps   xmm9, xmm2
	mulps   xmm10, xmm5
    mulps   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movaps  xmm4, [rsp + nb334_half]
	mulps   xmm9, xmm4  ;# rinvH1M
	mulps   xmm10, xmm4 ;# rinvH2M
    mulps   xmm11, xmm4 ;# rinvMM

	movaps  [rsp + nb334_rinvH1M], xmm9
	movaps  [rsp + nb334_rinvH2M], xmm10
	movaps  [rsp + nb334_rinvMM], xmm11
	
	;# M interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movaps xmm1, [rsp + nb334_tsc]
    mulps  xmm0, xmm9  ;# r
    mulps  xmm3, xmm10
    mulps  xmm6, xmm11
    mulps  xmm0, xmm1 ;# rtab
    mulps  xmm3, xmm1
    mulps  xmm6, xmm1
    
    ;# truncate and convert to integers
    cvttps2dq xmm1, xmm0
    cvttps2dq xmm4, xmm3
    cvttps2dq xmm7, xmm6        
    
    ;# convert back to float
    cvtdq2ps  xmm2, xmm1
    cvtdq2ps  xmm5, xmm4
    cvtdq2ps  xmm8, xmm7
     
    ;# multiply by 4
    pslld   xmm1, 2
    pslld   xmm4, 2
    pslld   xmm7, 2
    
    ;# multiply by three (copy, mult. by two, add back)
    movaps  xmm10, xmm1
    movaps  xmm11, xmm4
    movaps  xmm12, xmm7
    pslld   xmm1, 1
    pslld   xmm4, 1
    pslld   xmm7, 1    
    paddd   xmm1, xmm10
    paddd   xmm4, xmm11
    paddd   xmm7, xmm12    
    
    ;# move to integer registers
    movhlps xmm13, xmm1
    movhlps xmm14, xmm4
    movhlps xmm15, xmm7
    movd    eax, xmm1
    movd    r8d, xmm4
    movd    r12d, xmm7
    movd    ecx, xmm13
    movd    r10d, xmm14
    movd    r14d, xmm15
    pshufd  xmm1, xmm1, 1
    pshufd  xmm4, xmm4, 1
    pshufd  xmm7, xmm7, 1
    pshufd  xmm13, xmm13, 1
    pshufd  xmm14, xmm14, 1
    pshufd  xmm15, xmm15, 1
    movd    ebx, xmm1
    movd    r9d, xmm4
    movd    r13d, xmm7    
    movd    edx, xmm13
    movd    r11d, xmm14
    movd    r15d, xmm15   
    
    mov  rsi, [rbp + nb334_VFtab]

    ;# calculate eps
    subps     xmm0, xmm2
    subps     xmm3, xmm5
    subps     xmm6, xmm8

    movaps    [rsp + nb334_epsH1], xmm0
    movaps    [rsp + nb334_epsH2], xmm3
    movaps    [rsp + nb334_epsM], xmm6

    ;# Load LOTS of table data
   	movlps xmm1, [rsi + rax*4]
   	movlps xmm5, [rsi + r8*4]
   	movlps xmm9, [rsi + r12*4]

	movlps xmm3, [rsi + rcx*4]
	movlps xmm7, [rsi + r10*4]
	movlps xmm11, [rsi + r14*4]

	movhps xmm1, [rsi + rbx*4]
	movhps xmm5, [rsi + r9*4]
	movhps xmm9, [rsi + r13*4]

	movhps xmm3, [rsi + rdx*4]
	movhps xmm7, [rsi + r11*4]
	movhps xmm11, [rsi + r15*4]

    movaps xmm0, xmm1
    movaps xmm4, xmm5
    movaps xmm8, xmm9
	shufps xmm0, xmm3, 136  ;# 10001000
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm8, xmm11, 136  ;# 10001000
	shufps xmm1, xmm3, 221  ;# 11011101
	shufps xmm5, xmm7, 221  ;# 11011101
	shufps xmm9, xmm11, 221  ;# 11011101
    
	movlps xmm3, [rsi + rax*4 + 8]
	movlps xmm7, [rsi + r8*4 + 8]
	movlps xmm11, [rsi + r12*4 + 8]
    
	movlps xmm12, [rsi + rcx*4 + 8]
	movlps xmm13, [rsi + r10*4 + 8]
	movlps xmm14, [rsi + r14*4 + 8]

	movhps xmm3, [rsi + rbx*4 + 8]
	movhps xmm7, [rsi + r9*4 + 8]
	movhps xmm11, [rsi + r13*4 + 8]
    
	movhps xmm12, [rsi + rdx*4 + 8]
	movhps xmm13, [rsi + r11*4 + 8]
	movhps xmm14, [rsi + r15*4 + 8]

    movaps xmm2, xmm3
    movaps xmm6, xmm7
    movaps xmm10, xmm11
    
	shufps xmm2, xmm12, 136  ;# 10001000
	shufps xmm6, xmm13, 136  ;# 10001000
	shufps xmm10, xmm14, 136  ;# 10001000
	shufps xmm3, xmm12, 221  ;# 11011101
	shufps xmm7, xmm13, 221  ;# 11011101
	shufps xmm11, xmm14, 221  ;# 11011101
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    movaps xmm12, [rsp + nb334_epsH1]
    movaps xmm13, [rsp + nb334_epsH2]
    movaps xmm14, [rsp + nb334_epsM]
    
    mulps  xmm3, xmm12   ;# Heps
    mulps  xmm7, xmm13
    mulps  xmm11, xmm14 
    mulps  xmm2, xmm12   ;# Geps
    mulps  xmm6, xmm13
    mulps  xmm10, xmm14 
    mulps  xmm3, xmm12   ;# Heps2
    mulps  xmm7, xmm13
    mulps  xmm11, xmm14 

    addps  xmm1, xmm2   ;# F+Geps
    addps  xmm5, xmm6
    addps  xmm9, xmm10 
    addps  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addps  xmm5, xmm7
    addps  xmm9, xmm11 
    addps  xmm3, xmm3    ;# 2*Heps2
    addps  xmm7, xmm7
    addps  xmm11, xmm11
    addps  xmm3, xmm2    ;# 2*Heps2+Geps
    addps  xmm7, xmm6  
    addps  xmm11, xmm10
    addps  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addps  xmm7, xmm5
    addps  xmm11, xmm9
    mulps  xmm1, xmm12   ;# eps*Fp
    mulps  xmm5, xmm13
    mulps  xmm9, xmm14
    movaps xmm12, [rsp + nb334_qqMH]
    movaps xmm13, [rsp + nb334_qqMM]
    addps  xmm1, xmm0     ;# VV
    addps  xmm5, xmm4
    addps  xmm9, xmm8
    mulps  xmm1, xmm12   ;# VV*qq = vcoul
    mulps  xmm5, xmm12
    mulps  xmm9, xmm13
    mulps  xmm3, xmm12    ;# FF*qq = fij
    mulps  xmm7, xmm12
    mulps  xmm11, xmm13
    
    ;# accumulate vctot
    addps  xmm1, [rsp + nb334_vctot]
    addps  xmm5, xmm9
    addps  xmm1, xmm5
    movaps [rsp + nb334_vctot], xmm1
    
    movaps xmm10, [rsp + nb334_tsc]
    mulps  xmm3, xmm10  ;# fscal
    mulps  xmm7, xmm10
    mulps  xmm10, xmm11
        
    movd eax, mm0 ;# restore j3 from mm0-mm3
    movd ebx, mm1
    movd ecx, mm2
    movd edx, mm3
    
	;# move j M forces to local temp variables 
    mov rdi, [rbp + nb334_faction]
    movlps xmm11, [rdi + rax*4 + 36] ;# jxMa jyMa  -   -
    movlps xmm12, [rdi + rcx*4 + 36] ;# jxMc jyMc  -   -
    movhps xmm11, [rdi + rbx*4 + 36] ;# jxMa jyMa jxMb jyMb 
    movhps xmm12, [rdi + rdx*4 + 36] ;# jxMc jyMc jxMd jyMd 

    movss  xmm13, [rdi + rax*4 + 44] ;# jzMa  -  -  -
    movss  xmm14, [rdi + rcx*4 + 44] ;# jzMc  -  -  -
    movss  xmm1,  [rdi + rbx*4 + 44] ;# jzMb  -  -  -
    movss  xmm2,  [rdi + rdx*4 + 44] ;# jzMd  -  -  -
    movlhps xmm13, xmm1 ;# jzMa  -  jzMb  -
    movlhps xmm14, xmm2 ;# jzMc  -  jzMd -
    
    shufps xmm13, xmm14,  136  ;# 10001000 => jzMa jzMb jzMc jzMd

    ;# xmm11: jxMa jyMa jxMb jyMb 
    ;# xmm12: jxMc jyMc jxMd jyMd
    ;# xmm13: jzMa jzMb jzMc jzMd

    xorps  xmm0, xmm0
    xorps  xmm4, xmm4
    xorps  xmm8, xmm8

    mulps  xmm3, [rsp + nb334_rinvH1M]
    mulps  xmm7, [rsp + nb334_rinvH2M]
    mulps  xmm10, [rsp + nb334_rinvMM]
    
    subps  xmm0, xmm3
    subps  xmm4, xmm7
    subps  xmm8, xmm10
    
    movaps xmm1, xmm0
    movaps xmm2, xmm0
    movaps xmm3, xmm4
    movaps xmm5, xmm4
    movaps xmm6, xmm8
    movaps xmm7, xmm8

	mulps xmm0, [rsp + nb334_dxH1M]
	mulps xmm1, [rsp + nb334_dyH1M]
	mulps xmm2, [rsp + nb334_dzH1M]
	mulps xmm3, [rsp + nb334_dxH2M]
	mulps xmm4, [rsp + nb334_dyH2M]
	mulps xmm5, [rsp + nb334_dzH2M]
	mulps xmm6, [rsp + nb334_dxMM]
	mulps xmm7, [rsp + nb334_dyMM]
	mulps xmm8, [rsp + nb334_dzMM]

    movaps xmm14, xmm0
    movaps xmm15, xmm1
    addps xmm13, xmm2
    addps xmm0, [rsp + nb334_fixH1]
    addps xmm1, [rsp + nb334_fiyH1]
    addps xmm2, [rsp + nb334_fizH1]

    addps xmm14, xmm3
    addps xmm15, xmm4
    addps xmm13, xmm5
    addps xmm3, [rsp + nb334_fixH2]
    addps xmm4, [rsp + nb334_fiyH2]
    addps xmm5, [rsp + nb334_fizH2]

    addps xmm14, xmm6
    addps xmm15, xmm7
    addps xmm13, xmm8
    addps xmm6, [rsp + nb334_fixM]
    addps xmm7, [rsp + nb334_fiyM]
    addps xmm8, [rsp + nb334_fizM]

    movaps [rsp + nb334_fixH1], xmm0
    movaps [rsp + nb334_fiyH1], xmm1
    movaps [rsp + nb334_fizH1], xmm2
    movaps [rsp + nb334_fixH2], xmm3
    movaps [rsp + nb334_fiyH2], xmm4
    movaps [rsp + nb334_fizH2], xmm5
    movaps [rsp + nb334_fixM], xmm6
    movaps [rsp + nb334_fiyM], xmm7
    movaps [rsp + nb334_fizM], xmm8
    
    ;# xmm14 = fMx
    ;# xmm15 = fMy
    ;# xmm13 = fMz
    movaps xmm0, xmm14
    unpcklps xmm14, xmm15
    unpckhps xmm0,  xmm15
    
    addps  xmm11, xmm14
    addps  xmm12, xmm0
    
    movhlps  xmm14, xmm13 ;# fMzc fMzd
    
    movlps [rdi + rax*4 + 36], xmm11
    movhps [rdi + rbx*4 + 36], xmm11
    movlps [rdi + rcx*4 + 36], xmm12
    movhps [rdi + rdx*4 + 36], xmm12
    movss  [rdi + rax*4 + 44], xmm13
    movss  [rdi + rcx*4 + 44], xmm14
    shufps xmm13, xmm13, 1
    shufps xmm14, xmm14, 1
    movss  [rdi + rbx*4 + 44], xmm13
    movss  [rdi + rdx*4 + 44], xmm14
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb334_innerk],  4
	jl    .nb334_single_check
	jmp   .nb334_unroll_loop
.nb334_single_check:
	add dword ptr [rsp + nb334_innerk],  4
	jnz   .nb334_single_loop
	jmp   .nb334_updateouterdata
.nb334_single_loop:
	mov   rdx, [rsp + nb334_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb334_innerjjnr],  4	

	mov rsi, [rbp + nb334_pos]
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
	movaps [rsp + nb334_jxO], xmm6
	movaps [rsp + nb334_jyO], xmm3
	movaps [rsp + nb334_jzO], xmm1

    movaps xmm0, xmm6   ;# jxO
    movaps xmm2, xmm1   ;# jzO
    movaps xmm1, xmm3   ;# jyO
    movaps xmm4, xmm3   ;# jyO
    movaps xmm3, xmm6   ;# jxO
    movaps xmm5, xmm2   ;# jzO
        
	;# do O and H1 in parallel
	subps  xmm0, [rsp + nb334_ixO] 
	subps  xmm1, [rsp + nb334_iyO] 
	subps  xmm2, [rsp + nb334_izO] 
	subps  xmm3, [rsp + nb334_ixH1] 
	subps  xmm4, [rsp + nb334_iyH1]
	subps  xmm5, [rsp + nb334_izH1]
	
	movaps [rsp + nb334_dxOO], xmm0
	movaps [rsp + nb334_dyOO], xmm1
	movaps [rsp + nb334_dzOO], xmm2
	movaps [rsp + nb334_dxH1H1], xmm3
	movaps [rsp + nb334_dyH1H1], xmm4
	movaps [rsp + nb334_dzH1H1], xmm5
	
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
	movaps [rsp + nb334_rsqOO], xmm0
	movaps [rsp + nb334_rsqH1H1], xmm4
	
	;# do 1/sqrt(x) for O and  H1
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb334_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb334_half] ;# rinv O - j water 
	mulps   xmm7, [rsp + nb334_half] ;# rinv H1 - j water  

	movaps [rsp + nb334_rinvOO], xmm3
	movaps [rsp + nb334_rinvH1H1], xmm7

	mov rsi, [rbp + nb334_VFtab]
	
	;# do O table LJ interaction
	movaps xmm0, xmm3
	movaps xmm1, xmm0
	mulss  xmm1, [rsp + nb334_rsqOO] ;# xmm1=r 
	mulss  xmm1, [rsp + nb334_tsc]

	cvttps2pi mm6, xmm1
	cvtpi2ps xmm3, mm6
	subss    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulss  xmm2, xmm2   	;# xmm2=eps2 
	pslld   mm6, 2

	movd ebx, mm6
	lea   rbx, [rbx + rbx*2]

	;# load dispersion table data into xmm4
	movlps xmm4, [rsi + rbx*4 + 16]
	movlps xmm6, [rsi + rbx*4 + 24]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 0x1
	shufps xmm7, xmm7, 0x1
	;# dispersion table YFGH ready in xmm4-xmm7
	mulss  xmm6, xmm1   	;# xmm6=Geps 
	mulss  xmm7, xmm2   	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7   	;# xmm5=Fp 
	mulss  xmm7, [rsp + nb334_two]   	;# two*Heps2 
	addss  xmm7, xmm6
	addss  xmm7, xmm5 ;# xmm7=FF 
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 

	movaps xmm4, [rsp + nb334_c6]
	mulss  xmm7, xmm4	;# fijD 
	mulss  xmm5, xmm4	;# Vvdw6 

	;# save scalar force in xmm3. Update Vvdwtot directly 
	addss  xmm5, [rsp + nb334_Vvdwtot]
	movaps xmm3, xmm7 ;# fscal 
	movss [rsp + nb334_Vvdwtot], xmm5
	
	;# load repulsion table data into xmm4
	movlps xmm4, [rsi + rbx*4 + 32]
	movlps xmm6, [rsi + rbx*4 + 40]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 0x1
	shufps xmm7, xmm7, 0x1
	;# repulsion table YFGH ready in xmm4-xmm7
	
	mulss  xmm6, xmm1   	;# xmm6=Geps 
	mulss  xmm7, xmm2   	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7   	;# xmm5=Fp 
	mulss  xmm7, [rsp + nb334_two]   	;# two*Heps2 
	addss  xmm7, xmm6
	addss  xmm7, xmm5 ;# xmm7=FF 
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 
 
	movaps xmm4, [rsp + nb334_c12]
	mulss  xmm7, xmm4 ;# fijR 
	mulss  xmm5, xmm4 ;# Vvdw12 
	addss  xmm7, xmm3

	addss  xmm5, [rsp + nb334_Vvdwtot]
	movss [rsp + nb334_Vvdwtot], xmm5

	xorps  xmm1, xmm1
	mulss xmm7, [rsp + nb334_tsc]
	mulss xmm7, xmm0
	subss  xmm1, xmm7

	movaps xmm0, xmm1
	movaps xmm2, xmm1		
	
	mulss  xmm0, [rsp + nb334_dxOO]
	mulss  xmm1, [rsp + nb334_dyOO]
	mulss  xmm2, [rsp + nb334_dzOO]
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	addss   xmm3, xmm0
	addss   xmm4, xmm1
	addss   xmm5, xmm2
	movaps  [rsp + nb334_fjxO], xmm3
	movaps  [rsp + nb334_fjyO], xmm4
	movaps  [rsp + nb334_fjzO], xmm5
	addss   xmm0, [rsp + nb334_fixO]
	addss   xmm1, [rsp + nb334_fiyO]
	addss   xmm2, [rsp + nb334_fizO]
	movss  [rsp + nb334_fixO], xmm0
	movss  [rsp + nb334_fiyO], xmm1
	movss  [rsp + nb334_fizO], xmm2	

	;# do  H1 coulomb interaction
	movaps xmm0, [rsp + nb334_rinvH1H1] ;# rinv 
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb334_rsqH1H1] 	;# r
	mulps xmm1, [rsp + nb334_tsc]
	
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

	lea   rbx, [rbx + rbx*2]
	lea   rcx, [rcx + rcx*2]
	lea   rdx, [rdx + rdx*2]
	
	mov rsi, [rbp + nb334_VFtab]

	movlps xmm4, [rsi + rbx*4]
	movlps xmm3, [rsi + rcx*4]
	movlps xmm7, [rsi + rdx*4]
	movhps xmm4, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rdx*4 + 8]
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
	mulps  xmm7, [rsp + nb334_two]   	;# two*Heps2 

	xorps  xmm3, xmm3
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [rsp + nb334_qqHH]
	movhps  xmm3, [rsp + nb334_qqMH]
	shufps  xmm3, xmm3, 193 ;# 11000001 

	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point xmm5 contains vcoul and xmm3 fijC 
	
	addps  xmm5, [rsp + nb334_vctot]
	movaps [rsp + nb334_vctot], xmm5

	mulps  xmm3, [rsp + nb334_tsc]
	xorps  xmm2, xmm2
	subps  xmm2, xmm3
	mulps  xmm0, xmm2
	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [rsp + nb334_dxH1H1]
	mulps   xmm1, [rsp + nb334_dyH1H1]
	mulps   xmm2, [rsp + nb334_dzH1H1]
	;# update forces H1 - j water 
	movaps  xmm3, [rsp + nb334_fjxO]
	movaps  xmm4, [rsp + nb334_fjyO]
	movaps  xmm5, [rsp + nb334_fjzO]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  [rsp + nb334_fjxO], xmm3
	movaps  [rsp + nb334_fjyO], xmm4
	movaps  [rsp + nb334_fjzO], xmm5
	addps   xmm0, [rsp + nb334_fixH1]
	addps   xmm1, [rsp + nb334_fiyH1]
	addps   xmm2, [rsp + nb334_fizH1]
	movaps  [rsp + nb334_fixH1], xmm0
	movaps  [rsp + nb334_fiyH1], xmm1
	movaps  [rsp + nb334_fizH1], xmm2	
	
	;# i H2 & M simultaneously first get i particle coords: 
    movaps  xmm0, [rsp + nb334_jxO]
    movaps  xmm1, [rsp + nb334_jyO]
    movaps  xmm2, [rsp + nb334_jzO]
    movaps  xmm3, xmm0
    movaps  xmm4, xmm1
    movaps  xmm5, xmm2
	subps   xmm0, [rsp + nb334_ixH2]
	subps   xmm1, [rsp + nb334_iyH2]
	subps   xmm2, [rsp + nb334_izH2]
	subps   xmm3, [rsp + nb334_ixM] 
	subps   xmm4, [rsp + nb334_iyM] 
	subps   xmm5, [rsp + nb334_izM] 
	movaps [rsp + nb334_dxH2H2], xmm0
	movaps [rsp + nb334_dyH2H2], xmm1
	movaps [rsp + nb334_dzH2H2], xmm2
	movaps [rsp + nb334_dxMM], xmm3
	movaps [rsp + nb334_dyMM], xmm4
	movaps [rsp + nb334_dzMM], xmm5
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
	movaps [rsp + nb334_rsqH2H2], xmm0
	movaps [rsp + nb334_rsqMM], xmm4	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb334_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb334_half] ;# rinv H2 - j water 
	mulps   xmm7, [rsp + nb334_half] ;# rinv M - j water  

	movaps [rsp + nb334_rinvH2H2], xmm3
	movaps [rsp + nb334_rinvMM], xmm7

	movaps xmm1, xmm3
	mulps  xmm1, [rsp + nb334_rsqH2H2]	;# xmm1=r 
	movaps xmm0, xmm3	;# xmm0=rinv 
	mulps  xmm1, [rsp + nb334_tsc]
	
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

	lea   rbx, [rbx + rbx*2]
	lea   rcx, [rcx + rcx*2]
	lea   rdx, [rdx + rdx*2]

	movlps xmm4, [rsi + rbx*4]
	movlps xmm3, [rsi + rcx*4]
	movlps xmm7, [rsi + rdx*4]
	movhps xmm4, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rdx*4 + 8]
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
	mulps  xmm7, [rsp + nb334_two]   	;# two*Heps2 

	xorps  xmm3, xmm3
	
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [rsp + nb334_qqHH]
	movhps  xmm3, [rsp + nb334_qqMH]
	shufps  xmm3, xmm3, 193 ;# 11000001
	
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point xmm5 contains vcoul and xmm3 fijC 
	addps  xmm5, [rsp + nb334_vctot]
	movaps [rsp + nb334_vctot], xmm5	

	xorps  xmm1, xmm1

	mulps xmm3, [rsp + nb334_tsc]
	mulps xmm3, xmm0
	subps  xmm1, xmm3
	
	movaps  xmm0, xmm1
	movaps  xmm2, xmm1
	mulps   xmm0, [rsp + nb334_dxH2H2]
	mulps   xmm1, [rsp + nb334_dyH2H2]
	mulps   xmm2, [rsp + nb334_dzH2H2]
	;# update forces H1 - j water 
	movaps  xmm3, [rsp + nb334_fjxO]
	movaps  xmm4, [rsp + nb334_fjyO]
	movaps  xmm5, [rsp + nb334_fjzO]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  [rsp + nb334_fjxO], xmm3
	movaps  [rsp + nb334_fjyO], xmm4
	movaps  [rsp + nb334_fjzO], xmm5
	addps   xmm0, [rsp + nb334_fixH2]
	addps   xmm1, [rsp + nb334_fiyH2]
	addps   xmm2, [rsp + nb334_fizH2]
	movaps  [rsp + nb334_fixH2], xmm0
	movaps  [rsp + nb334_fiyH2], xmm1
	movaps  [rsp + nb334_fizH2], xmm2
	
	;# do table for i M - j water interaction 
	movaps xmm0, [rsp + nb334_rinvMM]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb334_rsqMM]	;# xmm0=rinv, xmm1=r 
	mulps  xmm1, [rsp + nb334_tsc]
	
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

	lea   rbx, [rbx + rbx*2]
	lea   rcx, [rcx + rcx*2]
	lea   rdx, [rdx + rdx*2]

        mov rsi, [rbp + nb334_VFtab]

        movlps xmm4, [rsi + rbx*4]
        movlps xmm3, [rsi + rcx*4]
        movlps xmm7, [rsi + rdx*4]
        movhps xmm4, [rsi + rbx*4 + 8]
        movhps xmm3, [rsi + rcx*4 + 8]
        movhps xmm7, [rsi + rdx*4 + 8]
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
	;; # coulomb table ready, in xmm4-xmm7

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm7, [rsp + nb334_two]   	;# two*Heps2 

	xorps  xmm3, xmm3
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [rsp + nb334_qqMH]
	movhps  xmm3, [rsp + nb334_qqMM]
	shufps  xmm3, xmm3, 193 ;# 11000001
	
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point xmm5 contains vcoul and xmm3 fijC 
	addps  xmm5, [rsp + nb334_vctot]
	movaps [rsp + nb334_vctot], xmm5	

	xorps  xmm1, xmm1

	mulps xmm3, [rsp + nb334_tsc]
	mulps xmm3, xmm0
	subps  xmm1, xmm3
	
	movaps  xmm0, xmm1
	movaps  xmm2, xmm1
	
	mulps   xmm0, [rsp + nb334_dxMM]
	mulps   xmm1, [rsp + nb334_dyMM]
	mulps   xmm2, [rsp + nb334_dzMM]
	movaps  xmm3, [rsp + nb334_fjxO]
	movaps  xmm4, [rsp + nb334_fjyO]
	movaps  xmm5, [rsp + nb334_fjzO]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	mov     rsi, [rbp + nb334_faction]
	movaps  [rsp + nb334_fjxO], xmm3
	movaps  [rsp + nb334_fjyO], xmm4
	movaps  [rsp + nb334_fjzO], xmm5
	addps   xmm0, [rsp + nb334_fixM]
	addps   xmm1, [rsp + nb334_fiyM]
	addps   xmm2, [rsp + nb334_fizM]
	movaps  [rsp + nb334_fixM], xmm0
	movaps  [rsp + nb334_fiyM], xmm1
	movaps  [rsp + nb334_fizM], xmm2

	;# update j water forces from local variables.
	;# transpose back first
	movaps  xmm0, [rsp + nb334_fjxO] ;# Ox H1x H2x Mx 
	movaps  xmm1, [rsp + nb334_fjyO] ;# Oy H1y H2y My
	movaps  xmm2, [rsp + nb334_fjzO] ;# Oz H1z H2z Mz
	 
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
	shufps   xmm3, xmm4, 36  ;# 00100100 ;# Ox Oy Oz H1x 
	shufps   xmm5, xmm6, 78  ;# 01001110 ;# h1y h1z h2x h2y
	shufps   xmm2, xmm1, 232  ;# 11101000 ;# h2z mx my mz

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

	dec dword ptr [rsp + nb334_innerk]
	jz    .nb334_updateouterdata
	jmp   .nb334_single_loop
.nb334_updateouterdata:
	mov   ecx, [rsp + nb334_ii3]
	mov   rdi, [rbp + nb334_faction]
	mov   rsi, [rbp + nb334_fshift]
	mov   edx, [rsp + nb334_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb334_fixO]
	movaps xmm1, [rsp + nb334_fiyO] 
	movaps xmm2, [rsp + nb334_fizO]

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
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4],     xmm3
	movss  [rdi + rcx*4 + 4], xmm4
	movss  [rdi + rcx*4 + 8], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# 00001000	

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb334_fixH1]
	movaps xmm1, [rsp + nb334_fiyH1]
	movaps xmm2, [rsp + nb334_fizH1]

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
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4 + 12], xmm3
	movss  [rdi + rcx*4 + 16], xmm4
	movss  [rdi + rcx*4 + 20], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# 00001000	
	addps   xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb334_fixH2]
	movaps xmm1, [rsp + nb334_fiyH2]
	movaps xmm2, [rsp + nb334_fizH2]

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
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4 + 24], xmm3
	movss  [rdi + rcx*4 + 28], xmm4
	movss  [rdi + rcx*4 + 32], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# 00001000	
	addps   xmm6, xmm0

	;# accumulate Mi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb334_fixM]
	movaps xmm1, [rsp + nb334_fiyM]
	movaps xmm2, [rsp + nb334_fizM]

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
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
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
	subps  xmm3, xmm6
	subss  xmm4, xmm7
	movlps  [rsi + rdx*4],    xmm3
	movss  [rsi + rdx*4 + 8], xmm4

	;# get n from stack
	mov esi, [rsp + nb334_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb334_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb334_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb334_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb334_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb334_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb334_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb334_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb334_n], esi
        jmp .nb334_outer
.nb334_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb334_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb334_end
        ;# non-zero, do one more workunit
        jmp   .nb334_threadloop
.nb334_end:
	mov eax, [rsp + nb334_nouter]
	mov ebx, [rsp + nb334_ninner]
	mov rcx, [rbp + nb334_outeriter]
	mov rdx, [rbp + nb334_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 1856
	emms

    ;# Save xmm registers to stack
    movaps xmm6,  [rsp      ]
    movaps xmm7,  [rsp + 16 ]
    movaps xmm8,  [rsp + 32 ]
    movaps xmm9,  [rsp + 48 ]
    movaps xmm10, [rsp + 64 ]
    movaps xmm11, [rsp + 80 ]
    movaps xmm12, [rsp + 96 ]
    movaps xmm13, [rsp + 112]
    movaps xmm14, [rsp + 128]
    movaps xmm15, [rsp + 144]

    ;# Reset pointers after restoring xmm6-15
    add rsp, 168

    pop r15
    pop r14
    pop r13
    pop r12
    pop rdi
    pop rsi
    pop rbx
    
	pop	rbp
	ret




.globl nb_kernel334nf_x86_64_sse
.globl _nb_kernel334nf_x86_64_sse
nb_kernel334nf_x86_64_sse:	
_nb_kernel334nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb334nf_fshift,         16
.equiv          nb334nf_gid,            24
.equiv          nb334nf_pos,            32
.equiv          nb334nf_faction,        40
.equiv          nb334nf_charge,         48
.equiv          nb334nf_p_facel,        56
.equiv          nb334nf_argkrf,         64
.equiv          nb334nf_argcrf,         72
.equiv          nb334nf_Vc,             80
.equiv          nb334nf_type,           88
.equiv          nb334nf_p_ntype,        96
.equiv          nb334nf_vdwparam,       104
.equiv          nb334nf_Vvdw,           112
.equiv          nb334nf_p_tabscale,     120
.equiv          nb334nf_VFtab,          128
.equiv          nb334nf_invsqrta,       136
.equiv          nb334nf_dvda,           144
.equiv          nb334nf_p_gbtabscale,   152
.equiv          nb334nf_GBtab,          160
.equiv          nb334nf_p_nthreads,     168
.equiv          nb334nf_count,          176
.equiv          nb334nf_mtx,            184
.equiv          nb334nf_outeriter,      192
.equiv          nb334nf_inneriter,      200
.equiv          nb334nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
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
.equiv          nb334nf_rsqH2H2,        704
.equiv          nb334nf_rsqH2M,         720
.equiv          nb334nf_rsqMH1,         736
.equiv          nb334nf_rsqMH2,         752
.equiv          nb334nf_rsqMM,          768
.equiv          nb334nf_rinvOO,         784
.equiv          nb334nf_rinvH1H1,       800
.equiv          nb334nf_rinvH1H2,       816
.equiv          nb334nf_rinvH1M,        832
.equiv          nb334nf_rinvH2H1,       848
.equiv          nb334nf_rinvH2H2,       864
.equiv          nb334nf_rinvH2M,        880
.equiv          nb334nf_rinvMH1,        896
.equiv          nb334nf_rinvMH2,        912
.equiv          nb334nf_rinvMM,         928
.equiv          nb334nf_is3,            944
.equiv          nb334nf_ii3,            948
.equiv          nb334nf_nri,            952
.equiv          nb334nf_iinr,           960
.equiv          nb334nf_jindex,         968
.equiv          nb334nf_jjnr,           976
.equiv          nb334nf_shift,          984
.equiv          nb334nf_shiftvec,       992
.equiv          nb334nf_facel,          1000
.equiv          nb334nf_innerjjnr,      1008
.equiv          nb334nf_innerk,         1016
.equiv          nb334nf_n,              1020
.equiv          nb334nf_nn1,            1024
.equiv          nb334nf_nouter,         1028
.equiv          nb334nf_ninner,         1032

	push rbp
	mov  rbp, rsp
    
    ;# Push integer registers on stack
	push rbx
    push rsi
    push rdi
    push r12
    push r13
    push r14
    push r15

    ;# Make room for registers xmm6-xmm15 (10 registers=160 bytes)
    sub rsp, 168
    
    ;# Save xmm registers to stack
    movaps [rsp      ], xmm6
    movaps [rsp + 16 ], xmm7
    movaps [rsp + 32 ], xmm8
    movaps [rsp + 48 ], xmm9
    movaps [rsp + 64 ], xmm10
    movaps [rsp + 80 ], xmm11
    movaps [rsp + 96 ], xmm12
    movaps [rsp + 112], xmm13
    movaps [rsp + 128], xmm14
    movaps [rsp + 144], xmm15
    
	emms
	sub rsp, 1040		;# local variable stack space (n*16+8)
; .if 0    # block below only read by NASM - special calling convention on win64
%ifidn __OUTPUT_FORMAT__, win64
    ;# Adjust rbp to account for shadow space (32) & two extra args (2*8) on stack
    add rbp, 48
    ;# Adjust stack pointer for different alignment
    ;# Move around arguments to fit AMD64 convention below
    ;# AMD64 passes args in: rdi,rsi,rdx,rcx,r8,r9 + stack
    ;# win64 passes args in: rcx,rdx,r8,r9         + stack
    mov rdi, rcx
    mov rsi, rdx
    mov rdx, r8
    mov rcx, r9
    mov r8,  [rbp]
    mov r9,  [rbp + 8]
%endif
; .endif   # end NASM- and win64-specific block

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb334nf_nouter], eax
	mov [rsp + nb334nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb334nf_nri], edi
	mov [rsp + nb334nf_iinr], rsi
	mov [rsp + nb334nf_jindex], rdx
	mov [rsp + nb334nf_jjnr], rcx
	mov [rsp + nb334nf_shift], r8
	mov [rsp + nb334nf_shiftvec], r9
	mov rsi, [rbp + nb334nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb334nf_facel], xmm0

	mov rax, [rbp + nb334nf_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb334nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb334nf_half], eax
	movss xmm1, [rsp + nb334nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb334nf_half],  xmm1
	movaps [rsp + nb334nf_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb334nf_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]		;# ebx =ii 

	mov   rdx, [rbp + nb334nf_charge]
	movss xmm5, [rdx + rbx*4 + 4]	
	movss xmm3, [rdx + rbx*4 + 12]	
	movss xmm4, xmm3	
	mov rsi, [rbp + nb334nf_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb334nf_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb334nf_qqMM], xmm3
	movaps [rsp + nb334nf_qqMH], xmm4
	movaps [rsp + nb334nf_qqHH], xmm5
	
	xorps xmm0, xmm0
	mov   rdx, [rbp + nb334nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb334nf_p_ntype]
	imul  ecx, [rdi]  	;# rcx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   rax, [rbp + nb334nf_vdwparam]
	movlps xmm0, [rax + rdx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0x55
	movaps [rsp + nb334nf_c6], xmm0
	movaps [rsp + nb334nf_c12], xmm1

.nb334nf_threadloop:
        mov   rsi, [rbp + nb334nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb334nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb334nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb334nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb334nf_n], eax
        mov [rsp + nb334nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb334nf_outerstart
        jmp .nb334nf_end
	
.nb334nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb334nf_nouter]
	mov [rsp + nb334nf_nouter], ebx

.nb334nf_outer:
	mov   rax, [rsp + nb334nf_shift]  	;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb334nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb334nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb334nf_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb334nf_pos]	;# rax = base of pos[]  
	mov   [rsp + nb334nf_ii3], ebx	
	
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
	movaps [rsp + nb334nf_ixO], xmm3
	movaps [rsp + nb334nf_iyO], xmm4
	movaps [rsp + nb334nf_izO], xmm5
	movaps [rsp + nb334nf_ixH1], xmm6
	movaps [rsp + nb334nf_iyH1], xmm7

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
	movaps [rsp + nb334nf_izH1], xmm6
	movaps [rsp + nb334nf_ixH2], xmm0
	movaps [rsp + nb334nf_iyH2], xmm1
	movaps [rsp + nb334nf_izH2], xmm2
	movaps [rsp + nb334nf_ixM], xmm3
	movaps [rsp + nb334nf_iyM], xmm4
	movaps [rsp + nb334nf_izM], xmm5

	;# clear vctot 
	xorps xmm4, xmm4
	movaps [rsp + nb334nf_vctot], xmm4
	movaps [rsp + nb334nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb334nf_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb334nf_pos]
	mov   rdi, [rbp + nb334nf_faction]	
	mov   rax, [rsp + nb334nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb334nf_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb334nf_ninner]
	mov   [rsp + nb334nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb334nf_innerk], edx	;# number of innerloop atoms 
	jge   .nb334nf_unroll_loop
	jmp   .nb334nf_single_check
.nb334nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb334nf_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb334nf_innerjjnr],  16 ;# advance pointer (unroll 4) 

	mov rsi, [rbp + nb334nf_pos]   	;# base of pos[] 

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
	movaps [rsp + nb334nf_jxO], xmm0
	movaps [rsp + nb334nf_jyO], xmm1
	movaps [rsp + nb334nf_jzO], xmm2
	movaps [rsp + nb334nf_jxH1], xmm3

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
	movaps [rsp + nb334nf_jyH1], xmm0
	movaps [rsp + nb334nf_jzH1], xmm1
	movaps [rsp + nb334nf_jxH2], xmm2
	movaps [rsp + nb334nf_jyH2], xmm3

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
	movaps [rsp + nb334nf_jzH2], xmm0
	movaps [rsp + nb334nf_jxM], xmm1
	movaps [rsp + nb334nf_jyM], xmm2
	movaps [rsp + nb334nf_jzM], xmm3
	
	;# start calculating pairwise distances
	movaps xmm0, [rsp + nb334nf_ixO]
	movaps xmm1, [rsp + nb334nf_iyO]
	movaps xmm2, [rsp + nb334nf_izO]
	movaps xmm3, [rsp + nb334nf_ixH1]
	movaps xmm4, [rsp + nb334nf_iyH1]
	movaps xmm5, [rsp + nb334nf_izH1]
	subps  xmm0, [rsp + nb334nf_jxO]
	subps  xmm1, [rsp + nb334nf_jyO]
	subps  xmm2, [rsp + nb334nf_jzO]
	subps  xmm3, [rsp + nb334nf_jxH1]
	subps  xmm4, [rsp + nb334nf_jyH1]
	subps  xmm5, [rsp + nb334nf_jzH1]
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
	movaps [rsp + nb334nf_rsqOO], xmm0
	movaps [rsp + nb334nf_rsqH1H1], xmm3

	movaps xmm0, [rsp + nb334nf_ixH1]
	movaps xmm1, [rsp + nb334nf_iyH1]
	movaps xmm2, [rsp + nb334nf_izH1]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [rsp + nb334nf_jxH2]
	subps  xmm1, [rsp + nb334nf_jyH2]
	subps  xmm2, [rsp + nb334nf_jzH2]
	subps  xmm3, [rsp + nb334nf_jxM]
	subps  xmm4, [rsp + nb334nf_jyM]
	subps  xmm5, [rsp + nb334nf_jzM]
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
	movaps [rsp + nb334nf_rsqH1H2], xmm0
	movaps [rsp + nb334nf_rsqH1M], xmm3

	movaps xmm0, [rsp + nb334nf_ixH2]
	movaps xmm1, [rsp + nb334nf_iyH2]
	movaps xmm2, [rsp + nb334nf_izH2]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [rsp + nb334nf_jxH1]
	subps  xmm1, [rsp + nb334nf_jyH1]
	subps  xmm2, [rsp + nb334nf_jzH1]
	subps  xmm3, [rsp + nb334nf_jxH2]
	subps  xmm4, [rsp + nb334nf_jyH2]
	subps  xmm5, [rsp + nb334nf_jzH2]
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
	movaps [rsp + nb334nf_rsqH2H1], xmm0
	movaps [rsp + nb334nf_rsqH2H2], xmm3

	movaps xmm0, [rsp + nb334nf_ixH2]
	movaps xmm1, [rsp + nb334nf_iyH2]
	movaps xmm2, [rsp + nb334nf_izH2]
	movaps xmm3, [rsp + nb334nf_ixM]
	movaps xmm4, [rsp + nb334nf_iyM]
	movaps xmm5, [rsp + nb334nf_izM]
	subps  xmm0, [rsp + nb334nf_jxM]
	subps  xmm1, [rsp + nb334nf_jyM]
	subps  xmm2, [rsp + nb334nf_jzM]
	subps  xmm3, [rsp + nb334nf_jxH1]
	subps  xmm4, [rsp + nb334nf_jyH1]
	subps  xmm5, [rsp + nb334nf_jzH1]
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
	movaps [rsp + nb334nf_rsqH2M], xmm0
	movaps [rsp + nb334nf_rsqMH1], xmm4

	movaps xmm0, [rsp + nb334nf_ixM]
	movaps xmm1, [rsp + nb334nf_iyM]
	movaps xmm2, [rsp + nb334nf_izM]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [rsp + nb334nf_jxH2]
	subps  xmm1, [rsp + nb334nf_jyH2]
	subps  xmm2, [rsp + nb334nf_jzH2]
	subps  xmm3, [rsp + nb334nf_jxM]
	subps  xmm4, [rsp + nb334nf_jyM]
	subps  xmm5, [rsp + nb334nf_jzM]
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
	movaps [rsp + nb334nf_rsqMH2], xmm0
	movaps [rsp + nb334nf_rsqMM], xmm4
	
	;# Invsqrt for O-O
	rsqrtps  xmm1, [rsp + nb334nf_rsqOO]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb334nf_three]	
	mulps   xmm1, [rsp + nb334nf_rsqOO]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + nb334nf_half] ;# rinvOO
	movaps [rsp + nb334nf_rinvOO], xmm3
	
	;# Invsqrt for H1-H1 and H1-H2
	rsqrtps xmm1, [rsp + nb334nf_rsqH1H1]
	rsqrtps xmm5, [rsp + nb334nf_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb334nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb334nf_rsqH1H1]
	mulps   xmm5, [rsp + nb334nf_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb334nf_half] ;# rinvH1H1 
	mulps   xmm7, [rsp + nb334nf_half] ;# rinvH1H2 
	movaps  [rsp + nb334nf_rinvH1H1], xmm3
	movaps  [rsp + nb334nf_rinvH1H2], xmm7
	
	rsqrtps xmm1, [rsp + nb334nf_rsqH1M]
	rsqrtps xmm5, [rsp + nb334nf_rsqH2H1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb334nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb334nf_rsqH1M]
	mulps   xmm5, [rsp + nb334nf_rsqH2H1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb334nf_half] 
	mulps   xmm7, [rsp + nb334nf_half]
	movaps  [rsp + nb334nf_rinvH1M], xmm3
	movaps  [rsp + nb334nf_rinvH2H1], xmm7
	
	rsqrtps xmm1, [rsp + nb334nf_rsqH2H2]
	rsqrtps xmm5, [rsp + nb334nf_rsqH2M]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb334nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb334nf_rsqH2H2]
	mulps   xmm5, [rsp + nb334nf_rsqH2M]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb334nf_half] 
	mulps   xmm7, [rsp + nb334nf_half]
	movaps  [rsp + nb334nf_rinvH2H2], xmm3
	movaps  [rsp + nb334nf_rinvH2M], xmm7
	
	rsqrtps xmm1, [rsp + nb334nf_rsqMH1]
	rsqrtps xmm5, [rsp + nb334nf_rsqMH2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb334nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb334nf_rsqMH1]
	mulps   xmm5, [rsp + nb334nf_rsqMH2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb334nf_half] 
	mulps   xmm7, [rsp + nb334nf_half]
	movaps  [rsp + nb334nf_rinvMH1], xmm3
	movaps  [rsp + nb334nf_rinvMH2], xmm7
        		
	rsqrtps xmm1, [rsp + nb334nf_rsqMM]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb334nf_three]
	mulps   xmm1, [rsp + nb334nf_rsqMM]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + nb334nf_half] 
	movaps  [rsp + nb334nf_rinvMM], xmm3

	;# start with OO table interaction
	movaps xmm0, [rsp + nb334nf_rinvOO]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb334nf_rsqOO] ;# xmm1=r
	mulps  xmm1, [rsp + nb334nf_tsc]
	
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

	mov  rsi, [rbp + nb334nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	lea   rax, [rax + rax*2]
	lea   rbx, [rbx + rbx*2]
	lea   rcx, [rcx + rcx*2]
	lea   rdx, [rdx + rdx*2]
	
	;# load dispersion table data into xmm4-xmm7
	movlps xmm5, [rsi + rax*4 + 16]
	movlps xmm7, [rsi + rcx*4 + 16]
	movhps xmm5, [rsi + rbx*4 + 16]
	movhps xmm7, [rsi + rdx*4 + 16] ;# got half table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 24]
	movlps xmm3, [rsi + rcx*4 + 24]
	movhps xmm7, [rsi + rbx*4 + 24]
	movhps xmm3, [rsi + rdx*4 + 24] ;# other half of table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# dispersion table YFGH ready in xmm4-xmm7
	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 

	movaps xmm4, [rsp + nb334nf_c6]
	mulps  xmm5, xmm4	;# Vvdw6 

	;# Update Vvdwtot directly 
	addps  xmm5, [rsp + nb334nf_Vvdwtot]
	movaps [rsp + nb334nf_Vvdwtot], xmm5

	;# load repulsion table data into xmm4-xmm7
	movlps xmm5, [rsi + rax*4 + 32]
	movlps xmm7, [rsi + rcx*4 + 32]
	movhps xmm5, [rsi + rbx*4 + 32]
	movhps xmm7, [rsi + rdx*4 + 32] ;# got half table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 40]
	movlps xmm3, [rsi + rcx*4 + 40]
	movhps xmm7, [rsi + rbx*4 + 40]
	movhps xmm3, [rsi + rdx*4 + 40] ;# other half of table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# repulsion table YFGH ready in xmm4-xmm7
	
	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
 
	movaps xmm4, [rsp + nb334nf_c12]
	mulps  xmm5, xmm4 ;# Vvdw12 

	addps  xmm5, [rsp + nb334nf_Vvdwtot]
	movaps [rsp + nb334nf_Vvdwtot], xmm5

	;# Coulomb interactions - first H1H1
	movaps xmm0, [rsp + nb334nf_rinvH1H1]
	
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb334nf_rsqH1H1] ;# xmm1=r 
	mulps  xmm1, [rsp + nb334nf_tsc]
	
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

	mov  rsi, [rbp + nb334nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	lea   rax, [rax + rax*2]
	lea   rbx, [rbx + rbx*2]
	lea   rcx, [rcx + rcx*2]
	lea   rdx, [rdx + rdx*2]

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
	movaps xmm3, [rsp + nb334nf_qqHH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 
	;# update vctot 
    	addps  xmm5, [rsp + nb334nf_vctot]
        movaps [rsp + nb334nf_vctot], xmm5
	
	;# H1-H2 interaction 
	movaps xmm0, [rsp + nb334nf_rinvH1H2]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb334nf_rsqH1H2] ;# xmm1=r 
	mulps  xmm1, [rsp + nb334nf_tsc]	
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

	lea   rax, [rax + rax*2]
	lea   rbx, [rbx + rbx*2]
	lea   rcx, [rcx + rcx*2]
	lea   rdx, [rdx + rdx*2]

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
	movaps xmm3, [rsp + nb334nf_qqHH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 

	addps  xmm5, [rsp + nb334nf_vctot]
	movaps [rsp + nb334nf_vctot], xmm5

	;# H1-M interaction  
	movaps xmm0, [rsp + nb334nf_rinvH1M]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb334nf_rsqH1M] ;# xmm1=r 
	mulps  xmm1, [rsp + nb334nf_tsc]	
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
	
	lea   rax, [rax + rax*2]
	lea   rbx, [rbx + rbx*2]
	lea   rcx, [rcx + rcx*2]
	lea   rdx, [rdx + rdx*2]


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
	movaps xmm3, [rsp + nb334nf_qqMH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 

	addps  xmm5, [rsp + nb334nf_vctot]
	movaps [rsp + nb334nf_vctot], xmm5

	;# H2-H1 interaction 
	movaps xmm0, [rsp + nb334nf_rinvH2H1]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb334nf_rsqH2H1] ;# xmm1=r 
	mulps  xmm1, [rsp + nb334nf_tsc]	
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

	lea   rax, [rax + rax*2]
	lea   rbx, [rbx + rbx*2]
	lea   rcx, [rcx + rcx*2]
	lea   rdx, [rdx + rdx*2]

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
	movaps xmm3, [rsp + nb334nf_qqHH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 

	addps  xmm5, [rsp + nb334nf_vctot]
	movaps [rsp + nb334nf_vctot], xmm5

	;# H2-H2 interaction 
	movaps xmm0, [rsp + nb334nf_rinvH2H2]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb334nf_rsqH2H2] ;# xmm1=r 
	mulps  xmm1, [rsp + nb334nf_tsc]	
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

	lea   rax, [rax + rax*2]
	lea   rbx, [rbx + rbx*2]
	lea   rcx, [rcx + rcx*2]
	lea   rdx, [rdx + rdx*2]

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
	movaps xmm3, [rsp + nb334nf_qqHH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 

	addps  xmm5, [rsp + nb334nf_vctot]
	movaps [rsp + nb334nf_vctot], xmm5

	;# H2-M interaction 
	movaps xmm0, [rsp + nb334nf_rinvH2M]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb334nf_rsqH2M] ;# xmm1=r 
	mulps  xmm1, [rsp + nb334nf_tsc]
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

	lea   rax, [rax + rax*2]
	lea   rbx, [rbx + rbx*2]
	lea   rcx, [rcx + rcx*2]
	lea   rdx, [rdx + rdx*2]

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
	movaps xmm3, [rsp + nb334nf_qqMH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 

	addps  xmm5, [rsp + nb334nf_vctot]
	movaps [rsp + nb334nf_vctot], xmm5

	;# M-H1 interaction 
	movaps xmm0, [rsp + nb334nf_rinvMH1]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb334nf_rsqMH1] ;# xmm1=r 
	mulps  xmm1, [rsp + nb334nf_tsc]	
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

	lea   rax, [rax + rax*2]
	lea   rbx, [rbx + rbx*2]
	lea   rcx, [rcx + rcx*2]
	lea   rdx, [rdx + rdx*2]

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
	movaps xmm3, [rsp + nb334nf_qqMH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 

	addps  xmm5, [rsp + nb334nf_vctot]
	movaps [rsp + nb334nf_vctot], xmm5

	;# M-H2 interaction 
	movaps xmm0, [rsp + nb334nf_rinvMH2]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb334nf_rsqMH2] ;# xmm1=r 
	mulps  xmm1, [rsp + nb334nf_tsc]
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

	lea   rax, [rax + rax*2]
	lea   rbx, [rbx + rbx*2]
	lea   rcx, [rcx + rcx*2]
	lea   rdx, [rdx + rdx*2]

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
	movaps xmm3, [rsp + nb334nf_qqMH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 

	addps  xmm5, [rsp + nb334nf_vctot]
	movaps [rsp + nb334nf_vctot], xmm5

	;# M-M interaction 
	movaps xmm0, [rsp + nb334nf_rinvMM]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb334nf_rsqMM] ;# xmm1=r 
	mulps  xmm1, [rsp + nb334nf_tsc]	
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

	lea   rax, [rax + rax*2]
	lea   rbx, [rbx + rbx*2]
	lea   rcx, [rcx + rcx*2]
	lea   rdx, [rdx + rdx*2]

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
	movaps xmm3, [rsp + nb334nf_qqMM]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 

	addps  xmm5, [rsp + nb334nf_vctot]
	movaps [rsp + nb334nf_vctot], xmm5
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb334nf_innerk],  4
	jl    .nb334nf_single_check
	jmp   .nb334nf_unroll_loop
.nb334nf_single_check:
	add dword ptr [rsp + nb334nf_innerk],  4
	jnz   .nb334nf_single_loop
	jmp   .nb334nf_updateouterdata
.nb334nf_single_loop:
	mov   rdx, [rsp + nb334nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb334nf_innerjjnr],  4	

	mov rsi, [rbp + nb334nf_pos]
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
	movaps [rsp + nb334nf_jxO], xmm6
	movaps [rsp + nb334nf_jyO], xmm3
	movaps [rsp + nb334nf_jzO], xmm1

	;# do O and H1 in parallel
	movaps xmm0, [rsp + nb334nf_ixO]
	movaps xmm1, [rsp + nb334nf_iyO]
	movaps xmm2, [rsp + nb334nf_izO]
	movaps xmm3, [rsp + nb334nf_ixH1]
	movaps xmm4, [rsp + nb334nf_iyH1]
	movaps xmm5, [rsp + nb334nf_izH1]
	subps  xmm0, [rsp + nb334nf_jxO]
	subps  xmm1, [rsp + nb334nf_jyO]
	subps  xmm2, [rsp + nb334nf_jzO]
	subps  xmm3, [rsp + nb334nf_jxO]
	subps  xmm4, [rsp + nb334nf_jyO]
	subps  xmm5, [rsp + nb334nf_jzO]
	
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
	movaps [rsp + nb334nf_rsqOO], xmm0
	movaps [rsp + nb334nf_rsqH1H1], xmm4
	
	;# do 1/sqrt(x) for O and  H1
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb334nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb334nf_half] ;# rinv O - j water 
	mulps   xmm7, [rsp + nb334nf_half] ;# rinv H1 - j water  

	movaps [rsp + nb334nf_rinvOO], xmm3
	movaps [rsp + nb334nf_rinvH1H1], xmm7

	mov rsi, [rbp + nb334nf_VFtab]
	
	;# do O table LJ interaction
	movaps xmm0, xmm3
	movaps xmm1, xmm0
	mulss  xmm1, [rsp + nb334nf_rsqOO] ;# xmm1=r 
	mulss  xmm1, [rsp + nb334nf_tsc]

	cvttps2pi mm6, xmm1
	cvtpi2ps xmm3, mm6
	subss    xmm1, xmm3	;# xmm1=eps 
	movaps xmm2, xmm1
	mulss  xmm2, xmm2   	;# xmm2=eps2 
	pslld   mm6, 2

	movd ebx, mm6
	lea   rbx, [rbx + rbx*2]

	;# load dispersion table data into xmm4
	movlps xmm4, [rsi + rbx*4 + 16]
	movlps xmm6, [rsi + rbx*4 + 24]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 0x1
	shufps xmm7, xmm7, 0x1
	;# dispersion table YFGH ready in xmm4-xmm7
	mulss  xmm6, xmm1   	;# xmm6=Geps 
	mulss  xmm7, xmm2   	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7   	;# xmm5=Fp 
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 

	movaps xmm4, [rsp + nb334nf_c6]
	mulss  xmm5, xmm4	;# Vvdw6 

	;# Update Vvdwtot directly 
	addss  xmm5, [rsp + nb334nf_Vvdwtot]
	movss [rsp + nb334nf_Vvdwtot], xmm5
	
	;# load repulsion table data into xmm4
	movlps xmm4, [rsi + rbx*4 + 32]
	movlps xmm6, [rsi + rbx*4 + 40]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 0x1
	shufps xmm7, xmm7, 0x1
	;# repulsion table YFGH ready in xmm4-xmm7
	
	mulss  xmm6, xmm1   	;# xmm6=Geps 
	mulss  xmm7, xmm2   	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7   	;# xmm5=Fp 
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 
 
	movaps xmm4, [rsp + nb334nf_c12]
	mulss  xmm5, xmm4 ;# Vvdw12 

	addss  xmm5, [rsp + nb334nf_Vvdwtot]
	movss [rsp + nb334nf_Vvdwtot], xmm5

	;# do  H1 coulomb interaction
	movaps xmm0, [rsp + nb334nf_rinvH1H1] ;# rinv 
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb334nf_rsqH1H1] 	;# r
	mulps xmm1, [rsp + nb334nf_tsc]
	
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

	lea   rbx, [rbx + rbx*2]
	lea   rcx, [rcx + rcx*2]
	lea   rdx, [rdx + rdx*2]
	
	mov rsi, [rbp + nb334_VFtab]

	movlps xmm4, [rsi + rbx*4]
	movlps xmm3, [rsi + rcx*4]
	movlps xmm7, [rsi + rdx*4]
	movhps xmm4, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rdx*4 + 8]
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
	movss   xmm3, [rsp + nb334nf_qqHH]
	movhps  xmm3, [rsp + nb334nf_qqMH]
	shufps  xmm3, xmm3, 193 ;# 11000001 

	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point xmm5 contains vcoul 
	
	addps  xmm5, [rsp + nb334nf_vctot]
	movaps [rsp + nb334nf_vctot], xmm5

	;# i H2 & M simultaneously first get i particle coords: 
	movaps  xmm0, [rsp + nb334nf_ixH2]
	movaps  xmm1, [rsp + nb334nf_iyH2]
	movaps  xmm2, [rsp + nb334nf_izH2]	
	movaps  xmm3, [rsp + nb334nf_ixM] 
	movaps  xmm4, [rsp + nb334nf_iyM] 
	movaps  xmm5, [rsp + nb334nf_izM] 
	subps   xmm0, [rsp + nb334nf_jxO]
	subps   xmm1, [rsp + nb334nf_jyO]
	subps   xmm2, [rsp + nb334nf_jzO]
	subps   xmm3, [rsp + nb334nf_jxO]
	subps   xmm4, [rsp + nb334nf_jyO]
	subps   xmm5, [rsp + nb334nf_jzO]
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
	movaps [rsp + nb334nf_rsqH2H2], xmm0
	movaps [rsp + nb334nf_rsqMM], xmm4	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb334nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb334nf_half] ;# rinv H2 - j water 
	mulps   xmm7, [rsp + nb334nf_half] ;# rinv M - j water  

	movaps [rsp + nb334nf_rinvH2H2], xmm3
	movaps [rsp + nb334nf_rinvMM], xmm7

	movaps xmm1, xmm3
	mulps  xmm1, [rsp + nb334nf_rsqH2H2]	;# xmm1=r 
	movaps xmm0, xmm3	;# xmm0=rinv 
	mulps  xmm1, [rsp + nb334nf_tsc]
	
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

	lea   rbx, [rbx + rbx*2]
	lea   rcx, [rcx + rcx*2]
	lea   rdx, [rdx + rdx*2]

	movlps xmm4, [rsi + rbx*4]
	movlps xmm3, [rsi + rcx*4]
	movlps xmm7, [rsi + rdx*4]
	movhps xmm4, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rdx*4 + 8]
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
	movss   xmm3, [rsp + nb334nf_qqHH]
	movhps  xmm3, [rsp + nb334nf_qqMH]
	shufps  xmm3, xmm3, 193 ;# 11000001
	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point xmm5 contains vcoul 
	addps  xmm5, [rsp + nb334nf_vctot]
	movaps [rsp + nb334nf_vctot], xmm5	

	;# do table for i M - j water interaction 
	movaps xmm0, [rsp + nb334nf_rinvMM]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb334nf_rsqMM]	;# xmm0=rinv, xmm1=r 
	mulps  xmm1, [rsp + nb334nf_tsc]
	
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

	lea   rbx, [rbx + rbx*2]
	lea   rcx, [rcx + rcx*2]
	lea   rdx, [rdx + rdx*2]

        mov rsi, [rbp + nb334_VFtab]

        movlps xmm4, [rsi + rbx*4]
        movlps xmm3, [rsi + rcx*4]
        movlps xmm7, [rsi + rdx*4]
        movhps xmm4, [rsi + rbx*4 + 8]
        movhps xmm3, [rsi + rcx*4 + 8]
        movhps xmm7, [rsi + rdx*4 + 8]
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
	;; # coulomb table ready, in xmm4-xmm7

	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp 

	xorps  xmm3, xmm3
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [rsp + nb334nf_qqMH]
	movhps  xmm3, [rsp + nb334nf_qqMM]
	shufps  xmm3, xmm3, 193 ;# 11000001
	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point xmm5 contains vcoul
	addps  xmm5, [rsp + nb334nf_vctot]
	movaps [rsp + nb334nf_vctot], xmm5	

	dec dword ptr [rsp + nb334nf_innerk]
	jz    .nb334nf_updateouterdata
	jmp   .nb334nf_single_loop
.nb334nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb334nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb334nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb334nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb334nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb334nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb334nf_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb334nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb334nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb334nf_n], esi
        jmp .nb334nf_outer
.nb334nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb334nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb334nf_end
        ;# non-zero, do one more workunit
        jmp   .nb334nf_threadloop
.nb334nf_end:
	mov eax, [rsp + nb334nf_nouter]
	mov ebx, [rsp + nb334nf_ninner]
	mov rcx, [rbp + nb334nf_outeriter]
	mov rdx, [rbp + nb334nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 1040
	emms

    ;# Save xmm registers to stack
    movaps xmm6,  [rsp      ]
    movaps xmm7,  [rsp + 16 ]
    movaps xmm8,  [rsp + 32 ]
    movaps xmm9,  [rsp + 48 ]
    movaps xmm10, [rsp + 64 ]
    movaps xmm11, [rsp + 80 ]
    movaps xmm12, [rsp + 96 ]
    movaps xmm13, [rsp + 112]
    movaps xmm14, [rsp + 128]
    movaps xmm15, [rsp + 144]

    ;# Reset pointers after restoring xmm6-15
    add rsp, 168

    pop r15
    pop r14
    pop r13
    pop r12
    pop rdi
    pop rsi
    pop rbx
    
	pop	rbp
	ret
