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


	
.globl nb_kernel114_x86_64_sse 
.globl _nb_kernel114_x86_64_sse
nb_kernel114_x86_64_sse:	
_nb_kernel114_x86_64_sse:	
.equiv          nb114_fshift,           16
.equiv          nb114_gid,              24
.equiv          nb114_pos,              32
.equiv          nb114_faction,          40
.equiv          nb114_charge,           48
.equiv          nb114_p_facel,          56
.equiv          nb114_argkrf,           64
.equiv          nb114_argcrf,           72
.equiv          nb114_Vc,               80
.equiv          nb114_type,             88
.equiv          nb114_p_ntype,          96
.equiv          nb114_vdwparam,         104
.equiv          nb114_Vvdw,             112
.equiv          nb114_p_tabscale,       120
.equiv          nb114_VFtab,            128
.equiv          nb114_invsqrta,         136
.equiv          nb114_dvda,             144
.equiv          nb114_p_gbtabscale,     152
.equiv          nb114_GBtab,            160
.equiv          nb114_p_nthreads,       168
.equiv          nb114_count,            176
.equiv          nb114_mtx,              184
.equiv          nb114_outeriter,        192
.equiv          nb114_inneriter,        200
.equiv          nb114_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb114_ixO,              0
.equiv          nb114_iyO,              16
.equiv          nb114_izO,              32
.equiv          nb114_ixH1,             48
.equiv          nb114_iyH1,             64
.equiv          nb114_izH1,             80
.equiv          nb114_ixH2,             96
.equiv          nb114_iyH2,             112
.equiv          nb114_izH2,             128
.equiv          nb114_ixM,              144
.equiv          nb114_iyM,              160
.equiv          nb114_izM,              176
.equiv          nb114_jxO,              192
.equiv          nb114_jyO,              208
.equiv          nb114_jzO,              224
.equiv          nb114_jxH1,             240
.equiv          nb114_jyH1,             256
.equiv          nb114_jzH1,             272
.equiv          nb114_jxH2,             288
.equiv          nb114_jyH2,             304
.equiv          nb114_jzH2,             320
.equiv          nb114_jxM,              336
.equiv          nb114_jyM,              352
.equiv          nb114_jzM,              368
.equiv          nb114_dxOO,             384
.equiv          nb114_dyOO,             400
.equiv          nb114_dzOO,             416
.equiv          nb114_dxH1H1,           432
.equiv          nb114_dyH1H1,           448
.equiv          nb114_dzH1H1,           464
.equiv          nb114_dxH1H2,           480
.equiv          nb114_dyH1H2,           496
.equiv          nb114_dzH1H2,           512
.equiv          nb114_dxH1M,            528
.equiv          nb114_dyH1M,            544
.equiv          nb114_dzH1M,            560
.equiv          nb114_dxH2H1,           576
.equiv          nb114_dyH2H1,           592
.equiv          nb114_dzH2H1,           608
.equiv          nb114_dxH2H2,           624
.equiv          nb114_dyH2H2,           640
.equiv          nb114_dzH2H2,           656
.equiv          nb114_dxH2M,            672
.equiv          nb114_dyH2M,            688
.equiv          nb114_dzH2M,            704
.equiv          nb114_dxMH1,            720
.equiv          nb114_dyMH1,            736
.equiv          nb114_dzMH1,            752
.equiv          nb114_dxMH2,            768
.equiv          nb114_dyMH2,            784
.equiv          nb114_dzMH2,            800
.equiv          nb114_dxMM,             816
.equiv          nb114_dyMM,             832
.equiv          nb114_dzMM,             848
.equiv          nb114_qqMM,             864
.equiv          nb114_qqMH,             880
.equiv          nb114_qqHH,             896
.equiv          nb114_two,              912
.equiv          nb114_c6,               928
.equiv          nb114_c12,              944
.equiv          nb114_six,              960
.equiv          nb114_twelve,           976
.equiv          nb114_vctot,            992
.equiv          nb114_Vvdwtot,          1008
.equiv          nb114_fixO,             1024
.equiv          nb114_fiyO,             1040
.equiv          nb114_fizO,             1056
.equiv          nb114_fixH1,            1072
.equiv          nb114_fiyH1,            1088
.equiv          nb114_fizH1,            1104
.equiv          nb114_fixH2,            1120
.equiv          nb114_fiyH2,            1136
.equiv          nb114_fizH2,            1152
.equiv          nb114_fixM,             1168
.equiv          nb114_fiyM,             1184
.equiv          nb114_fizM,             1200
.equiv          nb114_fjxO,             1216
.equiv          nb114_fjyO,             1232
.equiv          nb114_fjzO,             1248
.equiv          nb114_fjxH1,            1264
.equiv          nb114_fjyH1,            1280
.equiv          nb114_fjzH1,            1296
.equiv          nb114_fjxH2,            1312
.equiv          nb114_fjyH2,            1328
.equiv          nb114_fjzH2,            1344
.equiv          nb114_fjxM,             1360
.equiv          nb114_fjyM,             1376
.equiv          nb114_fjzM,             1392
.equiv          nb114_half,             1408
.equiv          nb114_three,            1424
.equiv          nb114_rsqOO,            1440
.equiv          nb114_rsqH1H1,          1456
.equiv          nb114_rsqH1H2,          1472
.equiv          nb114_rsqH1M,           1488
.equiv          nb114_rsqH2H1,          1504
.equiv          nb114_rsqH2H2,          1520
.equiv          nb114_rsqH2M,           1536
.equiv          nb114_rsqMH1,           1552
.equiv          nb114_rsqMH2,           1568
.equiv          nb114_rsqMM,            1584
.equiv          nb114_rinvsqOO,         1600
.equiv          nb114_rinvH1H1,         1616
.equiv          nb114_rinvH1H2,         1632
.equiv          nb114_rinvH1M,          1648
.equiv          nb114_rinvH2H1,         1664
.equiv          nb114_rinvH2H2,         1680
.equiv          nb114_rinvH2M,          1696
.equiv          nb114_rinvMH1,          1712
.equiv          nb114_rinvMH2,          1728
.equiv          nb114_rinvMM,           1744
.equiv          nb114_fstmp,            1760
.equiv          nb114_nri,              1776
.equiv          nb114_iinr,             1784
.equiv          nb114_jindex,           1792
.equiv          nb114_jjnr,             1800
.equiv          nb114_shift,            1808
.equiv          nb114_shiftvec,         1816
.equiv          nb114_facel,            1824
.equiv          nb114_innerjjnr,        1832
.equiv          nb114_is3,              1840
.equiv          nb114_ii3,              1844
.equiv          nb114_innerk,           1848
.equiv          nb114_n,                1856
.equiv          nb114_nn1,              1860
.equiv          nb114_nouter,           1864
.equiv          nb114_ninner,           1868

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
	sub rsp, 1872
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
	mov [rsp + nb114_nouter], eax
	mov [rsp + nb114_ninner], eax


	mov edi, [rdi]
	mov [rsp + nb114_nri], edi
	mov [rsp + nb114_iinr], rsi
	mov [rsp + nb114_jindex], rdx
	mov [rsp + nb114_jjnr], rcx
	mov [rsp + nb114_shift], r8
	mov [rsp + nb114_shiftvec], r9
	mov rsi, [rbp + nb114_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb114_facel], xmm0

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb114_half], eax
	movss xmm1, [rsp + nb114_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps xmm4, xmm3
	addps  xmm4, xmm4	;# six
	movaps xmm5, xmm4
	addps  xmm5, xmm5	;# twelve
	movaps [rsp + nb114_half],  xmm1
	movaps [rsp + nb114_two],  xmm2
	movaps [rsp + nb114_three],  xmm3
	movaps [rsp + nb114_six],  xmm4
	movaps [rsp + nb114_twelve],  xmm5

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb114_iinr]   ;# rcx = pointer into iinr[]
	mov   ebx, [rcx]		;# ebx =ii 

	mov   rdx, [rbp + nb114_charge]
	movss xmm5, [rdx + rbx*4 + 4]	
	movss xmm3, [rdx + rbx*4 + 12]	
	movss xmm4, xmm3	
	mov rsi, [rbp + nb114_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb114_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb114_qqMM], xmm3
	movaps [rsp + nb114_qqMH], xmm4
	movaps [rsp + nb114_qqHH], xmm5
	
	xorps xmm0, xmm0
	mov   rdx, [rbp + nb114_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb114_p_ntype]
	imul  ecx, [rdi]  ;# rcx = ntia = 2*ntype*type[ii0] 
	add   rdx, rcx
	mov   rax, [rbp + nb114_vdwparam]
	movlps xmm0, [rax + rdx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0x55
	movaps [rsp + nb114_c6], xmm0
	movaps [rsp + nb114_c12], xmm1

.nb114_threadloop:
        mov   rsi, [rbp + nb114_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb114_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   rbx, 1                           ;# rbx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb114_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb114_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb114_n], eax
        mov [rsp + nb114_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb114_outerstart
        jmp .nb114_end
	
.nb114_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb114_nouter]
	mov [rsp + nb114_nouter], ebx

.nb114_outer:
	mov   rax, [rsp + nb114_shift]  	;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb114_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb114_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb114_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb114_pos]	;# rax = base of pos[]  
	mov   [rsp + nb114_ii3], ebx	
	
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
	movaps [rsp + nb114_ixO], xmm3
	movaps [rsp + nb114_iyO], xmm4
	movaps [rsp + nb114_izO], xmm5
	movaps [rsp + nb114_ixH1], xmm6
	movaps [rsp + nb114_iyH1], xmm7

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
	movaps [rsp + nb114_izH1], xmm6
	movaps [rsp + nb114_ixH2], xmm0
	movaps [rsp + nb114_iyH2], xmm1
	movaps [rsp + nb114_izH2], xmm2
	movaps [rsp + nb114_ixM], xmm3
	movaps [rsp + nb114_iyM], xmm4
	movaps [rsp + nb114_izM], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb114_vctot], xmm4
	movaps [rsp + nb114_Vvdwtot], xmm4
	movaps [rsp + nb114_fixO], xmm4
	movaps [rsp + nb114_fiyO], xmm4
	movaps [rsp + nb114_fizO], xmm4
	movaps [rsp + nb114_fixH1], xmm4
	movaps [rsp + nb114_fiyH1], xmm4
	movaps [rsp + nb114_fizH1], xmm4
	movaps [rsp + nb114_fixH2], xmm4
	movaps [rsp + nb114_fiyH2], xmm4
	movaps [rsp + nb114_fizH2], xmm4
	movaps [rsp + nb114_fixM], xmm4
	movaps [rsp + nb114_fiyM], xmm4
	movaps [rsp + nb114_fizM], xmm4
	
	mov   rax, [rsp + nb114_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb114_pos]
	mov   rdi, [rbp + nb114_faction]	
	mov   rax, [rsp + nb114_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb114_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb114_ninner]
	mov   [rsp + nb114_ninner], ecx
	add   edx, 0
	mov   [rsp + nb114_innerk], edx	;# number of innerloop atoms 
	jge   .nb114_unroll_loop
	jmp   .nb114_single_check
.nb114_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb114_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb114_innerjjnr],  16 ;# advance pointer (unroll 4) 

	mov rsi, [rbp + nb114_pos]   	;# base of pos[] 

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
        
    subps xmm0, [rsp + nb114_ixO]
    subps xmm1, [rsp + nb114_iyO]
    subps xmm2, [rsp + nb114_izO]

    movaps xmm4, xmm0
    movaps xmm5, xmm1
    movaps xmm6, xmm2
    
    ;# square it
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
    
   	addps  xmm1, xmm0
	addps  xmm1, xmm2
    ;# rsq in xmm1

	;# move j O forces to local temp variables 
    movlps xmm10, [rdi + rax*4] ;# jxOa jyOa  -   -
    movlps xmm11, [rdi + rcx*4] ;# jxOc jyOc  -   -
    movhps xmm10, [rdi + rbx*4] ;# jxOa jyOa jxOb jyOb 
    movhps xmm11, [rdi + rdx*4] ;# jxOc jyOc jxOd jyOd 

    movss  xmm12, [rdi + rax*4 + 8] ;# jzOa  -  -  -
    movss  xmm13, [rdi + rcx*4 + 8] ;# jzOc  -  -  -
    movhps xmm12, [rdi + rbx*4 + 8] ;# jzOa  -  jzOb  -
    movhps xmm13, [rdi + rdx*4 + 8] ;# jzOc  -  jzOd -

    shufps xmm12, xmm13,  136  ;# 10001000 => jzOa jzOb jzOc jzOd

    ;# xmm10: jxOa jyOa jxOb jyOb 
    ;# xmm11: jxOc jyOc jxOd jyOd
    ;# xmm12: jzOa jzOb jzOc jzOd

    ;# calc rinvsq=1/rsq
	rcpps xmm2, xmm1
	movaps xmm0, [rsp + nb114_two]
	mulps xmm1, xmm2
	subps xmm0, xmm1
	mulps xmm0, xmm2	;# xmm0=rinvsq

	movaps xmm1, xmm0   ;# rinvsq

	mulps  xmm0, xmm0   ;# rinv4
	mulps  xmm0, xmm1	;# rinv6
	movaps xmm2, xmm0
	mulps  xmm2, xmm2	;# xmm2=rinv12

	mulps  xmm0, [rsp + nb114_c6]
	mulps  xmm2, [rsp + nb114_c12]
	movaps xmm8, xmm2
    subps  xmm2, xmm0	;# Vvdw=Vvdw12-Vvdw6 
	mulps  xmm0, [rsp + nb114_six]
	mulps  xmm8, [rsp + nb114_twelve]
	subps  xmm8, xmm0
	mulps  xmm8, xmm1	;# xmm8=total fscal 
        
    ;# add potential to Vvdwtot
	addps  xmm2, [rsp + nb114_Vvdwtot]
    movaps [rsp + nb114_Vvdwtot], xmm2

    ;# calculate scalar force by multiplying dx/dy/dz with fscal
	mulps  xmm4, xmm8
	mulps  xmm5, xmm8
	mulps  xmm6, xmm8

    ;# increment i force
    movaps xmm0, [rsp + nb114_fixO]
    movaps xmm1, [rsp + nb114_fiyO]
    movaps xmm2, [rsp + nb114_fizO]
    addps  xmm0, xmm4
    addps  xmm1, xmm5
    addps  xmm2, xmm6
    movaps [rsp + nb114_fixO], xmm0
    movaps [rsp + nb114_fiyO], xmm1
    movaps [rsp + nb114_fizO], xmm2

    ;# update O forces
    ;# xmm3 = fH1x , xmm4 = fH1y
    movaps xmm3, xmm4
    unpcklps xmm4, xmm5   ;# fjx1 fjx1 fjy1 fjy2
    unpckhps xmm3, xmm5   ;# fjx3 fjx4 fjy3 fjy4
    
    addps xmm10, xmm4
    addps xmm11, xmm3
    addps xmm12, xmm6
    
    movhlps  xmm13, xmm12 ;# fH1zc fH1zd

    movlps [rdi + rax*4], xmm10
    movhps [rdi + rbx*4], xmm10
    movlps [rdi + rcx*4], xmm11
    movhps [rdi + rdx*4], xmm11
    movss  [rdi + rax*4 + 8], xmm12
    movss  [rdi + rcx*4 + 8], xmm13
    shufps xmm12, xmm12, 1
    shufps xmm13, xmm13, 1
    movss  [rdi + rbx*4 + 8], xmm12
    movss  [rdi + rdx*4 + 8], xmm13
    ;# done with OO interaction.
    
	;# move j H1 coordinates to local temp variables 
    movlps xmm0, [rsi + rax*4 + 12] ;# jxH1a jyH1a  -   -
    movlps xmm1, [rsi + rcx*4 + 12] ;# jxH1c jyH1c  -   -
    movhps xmm0, [rsi + rbx*4 + 12] ;# jxH1a jyH1a jxH1b jyH1b 
    movhps xmm1, [rsi + rdx*4 + 12] ;# jxH1c jyH1c jxH1d jyH1d 

    movss  xmm2, [rsi + rax*4 + 20] ;# jzH1a  -  -  -
    movss  xmm3, [rsi + rcx*4 + 20] ;# jzH1c  -  -  -
    movhps xmm2, [rsi + rbx*4 + 20] ;# jzH1a  -  jzH1b  -
    movhps xmm3, [rsi + rdx*4 + 20] ;# jzH1c  -  jzH1d -
    
    movaps xmm4, xmm0
    unpcklps xmm0, xmm1  ;# jxH1a jxH1c jyH1a jyH1c        
    unpckhps xmm4, xmm1  ;# jxH1b jxH1d jyH1b jyH1d
    movaps xmm1, xmm0
    unpcklps xmm0, xmm4 ;# x
    unpckhps xmm1, xmm4 ;# y

    shufps xmm2, xmm3,  136  ;# 10001000 => jzH1a jzH1b jzH1c jzH1d

    ;# xmm0 = H1x
    ;# xmm1 = H1y
    ;# xmm2 = H1z
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2    
    
    subps xmm0, [rsp + nb114_ixH1]
    subps xmm1, [rsp + nb114_iyH1]
    subps xmm2, [rsp + nb114_izH1]
    subps xmm3, [rsp + nb114_ixH2]
    subps xmm4, [rsp + nb114_iyH2]
    subps xmm5, [rsp + nb114_izH2]
    subps xmm6, [rsp + nb114_ixM]
    subps xmm7, [rsp + nb114_iyM]
    subps xmm8, [rsp + nb114_izM]
    
	movaps [rsp + nb114_dxH1H1], xmm0
	movaps [rsp + nb114_dyH1H1], xmm1
	movaps [rsp + nb114_dzH1H1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb114_dxH2H1], xmm3
	movaps [rsp + nb114_dyH2H1], xmm4
	movaps [rsp + nb114_dzH2H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb114_dxMH1], xmm6
	movaps [rsp + nb114_dyMH1], xmm7
	movaps [rsp + nb114_dzMH1], xmm8
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
		
	movaps  xmm9, [rsp + nb114_three]
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

	movaps  xmm0, [rsp + nb114_half]
	mulps   xmm9, xmm0  ;# rinvH1H1
	mulps   xmm10, xmm0 ;# rinvH2H1
    mulps   xmm11, xmm0 ;# rinvMH1
	
	;# H1 interactions 
    movaps xmm0, xmm9
    movaps xmm1, xmm10
    movaps xmm2, xmm11
    mulps  xmm9, xmm9
    mulps  xmm10, xmm10
    mulps  xmm11, xmm11
    mulps  xmm0, [rsp + nb114_qqHH] 
    mulps  xmm1, [rsp + nb114_qqHH] 
    mulps  xmm2, [rsp + nb114_qqMH] 
    mulps  xmm9, xmm0
    mulps  xmm10, xmm1
    mulps  xmm11, xmm2
    
    addps xmm0, [rsp + nb114_vctot] 
    addps xmm1, xmm2
    addps xmm0, xmm1
    movaps [rsp + nb114_vctot], xmm0
    
	;# move j H1 forces to local temp variables 
    movlps xmm0, [rdi + rax*4 + 12] ;# jxH1a jyH1a  -   -
    movlps xmm1, [rdi + rcx*4 + 12] ;# jxH1c jyH1c  -   -
    movhps xmm0, [rdi + rbx*4 + 12] ;# jxH1a jyH1a jxH1b jyH1b 
    movhps xmm1, [rdi + rdx*4 + 12] ;# jxH1c jyH1c jxH1d jyH1d 

    movss  xmm2, [rdi + rax*4 + 20] ;# jzH1a  -  -  -
    movss  xmm3, [rdi + rcx*4 + 20] ;# jzH1c  -  -  -
    movhps xmm2, [rdi + rbx*4 + 20] ;# jzH1a  -  jzH1b  -
    movhps xmm3, [rdi + rdx*4 + 20] ;# jzH1c  -  jzH1d -

    shufps xmm2, xmm3,  136  ;# 10001000 => jzH1a jzH1b jzH1c jzH1d

    ;# xmm0: jxH1a jyH1a jxH1b jyH1b 
    ;# xmm1: jxH1c jyH1c jxH1d jyH1d
    ;# xmm2: jzH1a jzH1b jzH1c jzH1d

    movaps xmm7, xmm9
    movaps xmm8, xmm9
    movaps xmm13, xmm11
    movaps xmm14, xmm11
    movaps xmm15, xmm11
    movaps xmm11, xmm10
    movaps xmm12, xmm10

	mulps xmm7, [rsp + nb114_dxH1H1]
	mulps xmm8, [rsp + nb114_dyH1H1]
	mulps xmm9, [rsp + nb114_dzH1H1]
	mulps xmm10, [rsp + nb114_dxH2H1]
	mulps xmm11, [rsp + nb114_dyH2H1]
	mulps xmm12, [rsp + nb114_dzH2H1]
	mulps xmm13, [rsp + nb114_dxMH1]
	mulps xmm14, [rsp + nb114_dyMH1]
	mulps xmm15, [rsp + nb114_dzMH1]

    movaps xmm3, xmm7
    movaps xmm4, xmm8
    addps xmm2, xmm9
    addps xmm7, [rsp + nb114_fixH1]
    addps xmm8, [rsp + nb114_fiyH1]
    addps xmm9, [rsp + nb114_fizH1]

    addps xmm3, xmm10
    addps xmm4, xmm11
    addps xmm2, xmm12
    addps xmm10, [rsp + nb114_fixH2]
    addps xmm11, [rsp + nb114_fiyH2]
    addps xmm12, [rsp + nb114_fizH2]

    addps xmm3, xmm13
    addps xmm4, xmm14
    addps xmm2, xmm15
    addps xmm13, [rsp + nb114_fixM]
    addps xmm14, [rsp + nb114_fiyM]
    addps xmm15, [rsp + nb114_fizM]

    movaps [rsp + nb114_fixH1], xmm7
    movaps [rsp + nb114_fiyH1], xmm8
    movaps [rsp + nb114_fizH1], xmm9
    movaps [rsp + nb114_fixH2], xmm10
    movaps [rsp + nb114_fiyH2], xmm11
    movaps [rsp + nb114_fizH2], xmm12
    movaps [rsp + nb114_fixM], xmm13
    movaps [rsp + nb114_fiyM], xmm14
    movaps [rsp + nb114_fizM], xmm15
    
    ;# xmm3 = fH1x , xmm4 = fH1y
    movaps xmm5, xmm3
    unpcklps xmm3, xmm4
    unpckhps xmm5, xmm4
    
    addps xmm0, xmm3
    addps xmm1, xmm5

    movhlps  xmm3, xmm2 ;# fH1zc fH1zd
    
    
    movlps [rdi + rax*4 + 12], xmm0
    movhps [rdi + rbx*4 + 12], xmm0
    movlps [rdi + rcx*4 + 12], xmm1
    movhps [rdi + rdx*4 + 12], xmm1
    movss  [rdi + rax*4 + 20], xmm2
    movss  [rdi + rcx*4 + 20], xmm3
    shufps xmm2, xmm2, 1
    shufps xmm3, xmm3, 1
    movss  [rdi + rbx*4 + 20], xmm2
    movss  [rdi + rdx*4 + 20], xmm3

	;# move j H2 coordinates to local temp variables 
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
    
    subps xmm0, [rsp + nb114_ixH1]
    subps xmm1, [rsp + nb114_iyH1]
    subps xmm2, [rsp + nb114_izH1]
    subps xmm3, [rsp + nb114_ixH2]
    subps xmm4, [rsp + nb114_iyH2]
    subps xmm5, [rsp + nb114_izH2]
    subps xmm6, [rsp + nb114_ixM]
    subps xmm7, [rsp + nb114_iyM]
    subps xmm8, [rsp + nb114_izM]
    
	movaps [rsp + nb114_dxH1H2], xmm0
	movaps [rsp + nb114_dyH1H2], xmm1
	movaps [rsp + nb114_dzH1H2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb114_dxH2H2], xmm3
	movaps [rsp + nb114_dyH2H2], xmm4
	movaps [rsp + nb114_dzH2H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb114_dxMH2], xmm6
	movaps [rsp + nb114_dyMH2], xmm7
	movaps [rsp + nb114_dzMH2], xmm8
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
		
	movaps  xmm9, [rsp + nb114_three]
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

	movaps  xmm0, [rsp + nb114_half]
	mulps   xmm9, xmm0  ;# rinvH1H2
	mulps   xmm10, xmm0 ;# rinvH2H2
    mulps   xmm11, xmm0 ;# rinvMH2
	
	;# H2 interactions 
    movaps xmm0, xmm9
    movaps xmm1, xmm10
    movaps xmm2, xmm11
    mulps  xmm9, xmm9
    mulps  xmm10, xmm10
    mulps  xmm11, xmm11
    mulps  xmm0, [rsp + nb114_qqHH] 
    mulps  xmm1, [rsp + nb114_qqHH] 
    mulps  xmm2, [rsp + nb114_qqMH] 
    mulps  xmm9, xmm0
    mulps  xmm10, xmm1
    mulps  xmm11, xmm2
    
    addps xmm0, [rsp + nb114_vctot] 
    addps xmm1, xmm2
    addps xmm0, xmm1
    movaps [rsp + nb114_vctot], xmm0
    
	;# move j H2 forces to local temp variables 
    movlps xmm0, [rdi + rax*4 + 24] ;# jxH2a jyH2a  -   -
    movlps xmm1, [rdi + rcx*4 + 24] ;# jxH2c jyH2c  -   -
    movhps xmm0, [rdi + rbx*4 + 24] ;# jxH2a jyH2a jxH2b jyH2b 
    movhps xmm1, [rdi + rdx*4 + 24] ;# jxH2c jyH2c jxH2d jyH2d 

    movss  xmm2, [rdi + rax*4 + 32] ;# jzH2a  -  -  -
    movss  xmm3, [rdi + rcx*4 + 32] ;# jzH2c  -  -  -
    movhps xmm2, [rdi + rbx*4 + 32] ;# jzH2a  -  jzH2b  -
    movhps xmm3, [rdi + rdx*4 + 32] ;# jzH2c  -  jzH2d -
    
    shufps xmm2, xmm3,  136  ;# 10001000 => jzH2a jzH2b jzH2c jzH2d

    ;# xmm0: jxH2a jyH2a jxH2b jyH2b 
    ;# xmm1: jxH2c jyH2c jxH2d jyH2d
    ;# xmm2: jzH2a jzH2b jzH2c jzH2d

    movaps xmm7, xmm9
    movaps xmm8, xmm9
    movaps xmm13, xmm11
    movaps xmm14, xmm11
    movaps xmm15, xmm11
    movaps xmm11, xmm10
    movaps xmm12, xmm10

	mulps xmm7, [rsp + nb114_dxH1H2]
	mulps xmm8, [rsp + nb114_dyH1H2]
	mulps xmm9, [rsp + nb114_dzH1H2]
	mulps xmm10, [rsp + nb114_dxH2H2]
	mulps xmm11, [rsp + nb114_dyH2H2]
	mulps xmm12, [rsp + nb114_dzH2H2]
	mulps xmm13, [rsp + nb114_dxMH2]
	mulps xmm14, [rsp + nb114_dyMH2]
	mulps xmm15, [rsp + nb114_dzMH2]

    movaps xmm3, xmm7
    movaps xmm4, xmm8
    addps xmm2, xmm9
    addps xmm7, [rsp + nb114_fixH1]
    addps xmm8, [rsp + nb114_fiyH1]
    addps xmm9, [rsp + nb114_fizH1]

    addps xmm3, xmm10
    addps xmm4, xmm11
    addps xmm2, xmm12
    addps xmm10, [rsp + nb114_fixH2]
    addps xmm11, [rsp + nb114_fiyH2]
    addps xmm12, [rsp + nb114_fizH2]

    addps xmm3, xmm13
    addps xmm4, xmm14
    addps xmm2, xmm15
    addps xmm13, [rsp + nb114_fixM]
    addps xmm14, [rsp + nb114_fiyM]
    addps xmm15, [rsp + nb114_fizM]
    
    movaps [rsp + nb114_fixH1], xmm7
    movaps [rsp + nb114_fiyH1], xmm8
    movaps [rsp + nb114_fizH1], xmm9
    movaps [rsp + nb114_fixH2], xmm10
    movaps [rsp + nb114_fiyH2], xmm11
    movaps [rsp + nb114_fizH2], xmm12
    movaps [rsp + nb114_fixM], xmm13
    movaps [rsp + nb114_fiyM], xmm14
    movaps [rsp + nb114_fizM], xmm15

    ;# xmm3 = fH2x , xmm4 = fH2y
    movaps xmm5, xmm3
    unpcklps xmm3, xmm4
    unpckhps xmm5, xmm4
    
    addps xmm0, xmm3
    addps xmm1, xmm5

    movhlps  xmm3, xmm2 ;# fH2zc fH2zd
    
    movlps [rdi + rax*4 + 24], xmm0
    movhps [rdi + rbx*4 + 24], xmm0
    movlps [rdi + rcx*4 + 24], xmm1
    movhps [rdi + rdx*4 + 24], xmm1
    movss  [rdi + rax*4 + 32], xmm2
    movss  [rdi + rcx*4 + 32], xmm3
    shufps xmm2, xmm2, 1 
    shufps xmm3, xmm3, 1
    movss  [rdi + rbx*4 + 32], xmm2
    movss  [rdi + rdx*4 + 32], xmm3

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
    
    subps xmm0, [rsp + nb114_ixH1]
    subps xmm1, [rsp + nb114_iyH1]
    subps xmm2, [rsp + nb114_izH1]
    subps xmm3, [rsp + nb114_ixH2]
    subps xmm4, [rsp + nb114_iyH2]
    subps xmm5, [rsp + nb114_izH2]
    subps xmm6, [rsp + nb114_ixM]
    subps xmm7, [rsp + nb114_iyM]
    subps xmm8, [rsp + nb114_izM]
    
	movaps [rsp + nb114_dxH1M], xmm0
	movaps [rsp + nb114_dyH1M], xmm1
	movaps [rsp + nb114_dzH1M], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb114_dxH2M], xmm3
	movaps [rsp + nb114_dyH2M], xmm4
	movaps [rsp + nb114_dzH2M], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb114_dxMM], xmm6
	movaps [rsp + nb114_dyMM], xmm7
	movaps [rsp + nb114_dzMM], xmm8
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
		
	movaps  xmm9, [rsp + nb114_three]
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

	movaps  xmm0, [rsp + nb114_half]
	mulps   xmm9, xmm0  ;# rinvH1M 
	mulps   xmm10, xmm0 ;# rinvH2M
    mulps   xmm11, xmm0 ;# rinvMM
	
	;# M interactions 
    movaps xmm0, xmm9
    movaps xmm1, xmm10
    movaps xmm2, xmm11
    mulps  xmm9, xmm9
    mulps  xmm10, xmm10
    mulps  xmm11, xmm11
    mulps  xmm0, [rsp + nb114_qqMH] 
    mulps  xmm1, [rsp + nb114_qqMH] 
    mulps  xmm2, [rsp + nb114_qqMM] 
    mulps  xmm9, xmm0
    mulps  xmm10, xmm1
    mulps  xmm11, xmm2
    
    addps xmm0, [rsp + nb114_vctot] 
    addps xmm1, xmm2
    addps xmm0, xmm1
    movaps [rsp + nb114_vctot], xmm0
    
	;# move j M forces to local temp variables 
    movlps xmm0, [rdi + rax*4 + 36] ;# jxMa jyMa  -   -
    movlps xmm1, [rdi + rcx*4 + 36] ;# jxMc jyMc  -   -
    movhps xmm0, [rdi + rbx*4 + 36] ;# jxMa jyMa jxMb jyMb 
    movhps xmm1, [rdi + rdx*4 + 36] ;# jxMc jyMc jxMd jyMd 

    movss  xmm2, [rdi + rax*4 + 44] ;# jzMa  -  -  -
    movss  xmm3, [rdi + rcx*4 + 44] ;# jzMc  -  -  -
    movss  xmm7, [rdi + rbx*4 + 44] ;# jzMb  -  -  -
    movss  xmm8, [rdi + rdx*4 + 44] ;# jzMd  -  -  -
    movlhps xmm2, xmm7 ;# jzMa  -  jzMb  -
    movlhps xmm3, xmm8 ;# jzMc  -  jzMd -
    
    shufps xmm2, xmm3,  136  ;# 10001000 => jzMa jzMb jzMc jzMd

    ;# xmm0: jxMa jyMa jxMb jyMb 
    ;# xmm1: jxMc jyMc jxMd jyMd
    ;# xmm2: jzMa jzMb jzMc jzMd

    movaps xmm7, xmm9
    movaps xmm8, xmm9
    movaps xmm13, xmm11
    movaps xmm14, xmm11
    movaps xmm15, xmm11
    movaps xmm11, xmm10
    movaps xmm12, xmm10

	mulps xmm7, [rsp + nb114_dxH1M]
	mulps xmm8, [rsp + nb114_dyH1M]
	mulps xmm9, [rsp + nb114_dzH1M]
	mulps xmm10, [rsp + nb114_dxH2M]
	mulps xmm11, [rsp + nb114_dyH2M]
	mulps xmm12, [rsp + nb114_dzH2M]
	mulps xmm13, [rsp + nb114_dxMM]
	mulps xmm14, [rsp + nb114_dyMM]
	mulps xmm15, [rsp + nb114_dzMM]

    movaps xmm3, xmm7
    movaps xmm4, xmm8
    addps xmm2, xmm9
    addps xmm7, [rsp + nb114_fixH1]
    addps xmm8, [rsp + nb114_fiyH1]
    addps xmm9, [rsp + nb114_fizH1]

    addps xmm3, xmm10
    addps xmm4, xmm11
    addps xmm2, xmm12
    addps xmm10, [rsp + nb114_fixH2]
    addps xmm11, [rsp + nb114_fiyH2]
    addps xmm12, [rsp + nb114_fizH2]

    addps xmm3, xmm13
    addps xmm4, xmm14
    addps xmm2, xmm15
    addps xmm13, [rsp + nb114_fixM]
    addps xmm14, [rsp + nb114_fiyM]
    addps xmm15, [rsp + nb114_fizM]
    
    movaps [rsp + nb114_fixH1], xmm7
    movaps [rsp + nb114_fiyH1], xmm8
    movaps [rsp + nb114_fizH1], xmm9
    movaps [rsp + nb114_fixH2], xmm10
    movaps [rsp + nb114_fiyH2], xmm11
    movaps [rsp + nb114_fizH2], xmm12
    movaps [rsp + nb114_fixM], xmm13
    movaps [rsp + nb114_fiyM], xmm14
    movaps [rsp + nb114_fizM], xmm15

    ;# xmm3 = fMx , xmm4 = fMy
    movaps xmm5, xmm3
    unpcklps xmm3, xmm4
    unpckhps xmm5, xmm4
    
    addps xmm0, xmm3
    addps xmm1, xmm5

    movhlps  xmm3, xmm2 ;# fMzc fMzd
    
    movlps [rdi + rax*4 + 36], xmm0
    movhps [rdi + rbx*4 + 36], xmm0
    movlps [rdi + rcx*4 + 36], xmm1
    movhps [rdi + rdx*4 + 36], xmm1
    movss  [rdi + rax*4 + 44], xmm2
    movss  [rdi + rcx*4 + 44], xmm3
    shufps xmm2, xmm2, 1
    shufps xmm3, xmm3, 1
    movss  [rdi + rbx*4 + 44], xmm2
    movss  [rdi + rdx*4 + 44], xmm3

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb114_innerk],  4
	jl    .nb114_single_check
	jmp   .nb114_unroll_loop
.nb114_single_check:
	add dword ptr [rsp + nb114_innerk],  4
	jnz   .nb114_single_loop
	jmp   .nb114_updateouterdata
.nb114_single_loop:
	mov   rdx, [rsp + nb114_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb114_innerjjnr],  4	

	mov rsi, [rbp + nb114_pos]
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
	movaps [rsp + nb114_jxO], xmm6
	movaps [rsp + nb114_jyO], xmm3
	movaps [rsp + nb114_jzO], xmm1

    movaps xmm0, xmm6   ;# jxO
    movaps xmm2, xmm1   ;# jzO
    movaps xmm1, xmm3   ;# jyO
    movaps xmm4, xmm3   ;# jyO
    movaps xmm3, xmm6   ;# jxO
    movaps xmm5, xmm2   ;# jzO
        
	;# do O and M in parallel
	subps  xmm0, [rsp + nb114_ixO] 
	subps  xmm1, [rsp + nb114_iyO] 
	subps  xmm2, [rsp + nb114_izO] 
	subps  xmm3, [rsp + nb114_ixM] 
	subps  xmm4, [rsp + nb114_iyM]
	subps  xmm5, [rsp + nb114_izM]
	
	movaps [rsp + nb114_dxOO], xmm0
	movaps [rsp + nb114_dyOO], xmm1
	movaps [rsp + nb114_dzOO], xmm2
	movaps [rsp + nb114_dxMM], xmm3
	movaps [rsp + nb114_dyMM], xmm4
	movaps [rsp + nb114_dzMM], xmm5
	
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
	movaps [rsp + nb114_rsqMM], xmm4
	
	;# do 1/x for O and 1/sqrt(x) for M
	rcpss  xmm1, xmm0
	rsqrtps xmm5, xmm4
	movss  xmm2, [rsp + nb114_two]
	movaps  xmm6, xmm5	
	mulss  xmm0, xmm1
	mulps   xmm5, xmm5
	subss  xmm2, xmm0
	movaps  xmm7, [rsp + nb114_three]
	mulss  xmm2, xmm1 	;# 1/r2
	
	mulps   xmm5, xmm4
	movss  xmm0, xmm2
	subps   xmm7, xmm5
	mulss  xmm2, xmm2
	mulps   xmm7, xmm6
	mulss  xmm2, xmm0 	;# 1/r6
	mulps   xmm7, [rsp + nb114_half] ;# rinv iH1 - j water 
	movss  xmm1, xmm2
	movaps [rsp + nb114_rinvMM], xmm7

	mulss  xmm2, xmm2 	;# 1/r12
	mulss  xmm1, [rsp + nb114_c6]
	mulss  xmm2, [rsp + nb114_c12]
	movss  xmm3, xmm2
	subss  xmm3, xmm1
	addss  xmm3, [rsp + nb114_Vvdwtot]
	movss  [rsp + nb114_Vvdwtot], xmm3
	mulss  xmm1, [rsp + nb114_six]
	mulss  xmm2, [rsp + nb114_twelve]
	subss  xmm2, xmm1
	mulss  xmm0, xmm2 	;# fscal
	movss  xmm1, xmm0
	movss  xmm2, xmm0
	mulss  xmm0, [rsp + nb114_dxOO]
	mulss  xmm1, [rsp + nb114_dyOO]
	mulss  xmm2, [rsp + nb114_dzOO]
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	addss   xmm3, xmm0
	addss   xmm4, xmm1
	addss   xmm5, xmm2
	movaps  [rsp + nb114_fjxO], xmm3
	movaps  [rsp + nb114_fjyO], xmm4
	movaps  [rsp + nb114_fjzO], xmm5
	addss   xmm0, [rsp + nb114_fixO]
	addss   xmm1, [rsp + nb114_fiyO]
	addss   xmm2, [rsp + nb114_fizO]
	movss  [rsp + nb114_fixO], xmm0
	movss  [rsp + nb114_fiyO], xmm1
	movss  [rsp + nb114_fizO], xmm2

	;# do  M coulomb interaction
	movaps xmm0, [rsp + nb114_rinvMM]
	movaps xmm4, xmm0	;# xmm4=rinv
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [rsp + nb114_qqMH]
	movhps  xmm3, [rsp + nb114_qqMM]
	shufps  xmm3, xmm3, 193 ;# 11000001 

	mulps  xmm4, xmm3 	;# xmm4=voul	
	mulps  xmm0, xmm4
	addps  xmm4, [rsp + nb114_vctot] 
	movaps [rsp + nb114_vctot], xmm4

	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [rsp + nb114_dxMM]
	mulps   xmm1, [rsp + nb114_dyMM]
	mulps   xmm2, [rsp + nb114_dzMM]
	;# update forces M - j water 
	movaps  xmm3, [rsp + nb114_fjxO]
	movaps  xmm4, [rsp + nb114_fjyO]
	movaps  xmm5, [rsp + nb114_fjzO]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  [rsp + nb114_fjxO], xmm3
	movaps  [rsp + nb114_fjyO], xmm4
	movaps  [rsp + nb114_fjzO], xmm5
	addps   xmm0, [rsp + nb114_fixM]
	addps   xmm1, [rsp + nb114_fiyM]
	addps   xmm2, [rsp + nb114_fizM]
	movaps  [rsp + nb114_fixM], xmm0
	movaps  [rsp + nb114_fiyM], xmm1
	movaps  [rsp + nb114_fizM], xmm2	
	
	;# i H1 & H2 simultaneously first get i particle coords: 
    movaps  xmm0, [rsp + nb114_jxO]
    movaps  xmm1, [rsp + nb114_jyO]
    movaps  xmm2, [rsp + nb114_jzO]
    movaps  xmm3, xmm0
    movaps  xmm4, xmm1
    movaps  xmm5, xmm2
	subps   xmm0, [rsp + nb114_ixH1]
	subps   xmm1, [rsp + nb114_iyH1]
	subps   xmm2, [rsp + nb114_izH1]
	subps   xmm3, [rsp + nb114_ixH2] 
	subps   xmm4, [rsp + nb114_iyH2] 
	subps   xmm5, [rsp + nb114_izH2] 
	movaps [rsp + nb114_dxH1H1], xmm0
	movaps [rsp + nb114_dyH1H1], xmm1
	movaps [rsp + nb114_dzH1H1], xmm2
	movaps [rsp + nb114_dxH2H2], xmm3
	movaps [rsp + nb114_dyH2H2], xmm4
	movaps [rsp + nb114_dzH2H2], xmm5
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
	movaps  [rsp + nb114_rsqH1H1], xmm0
	movaps  [rsp + nb114_rsqH2H2], xmm4
	
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb114_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb114_half] ;# rinvH1H1
	mulps   xmm7, [rsp + nb114_half] ;# rinvH2H2
	movaps  [rsp + nb114_rinvH1H1], xmm3
	movaps  [rsp + nb114_rinvH2H2], xmm7
	
	;# Do H1 coulomb interaction
	movaps xmm0, [rsp + nb114_rinvH1H1]
	movaps xmm4, xmm0	;# xmm4=rinv 
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [rsp + nb114_qqHH]
	movhps  xmm3, [rsp + nb114_qqMH]
	shufps  xmm3, xmm3, 193 ;# 11000001 

	mulps  xmm4, xmm3 	;# xmm4=voul
	mulps  xmm0, xmm4
	addps  xmm4, [rsp + nb114_vctot] 
	movaps [rsp + nb114_vctot], xmm4


	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [rsp + nb114_dxH1H1]
	mulps   xmm1, [rsp + nb114_dyH1H1]
	mulps   xmm2, [rsp + nb114_dzH1H1]
	;# update forces H1 - j water 
	movaps  xmm3, [rsp + nb114_fjxO]
	movaps  xmm4, [rsp + nb114_fjyO]
	movaps  xmm5, [rsp + nb114_fjzO]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  [rsp + nb114_fjxO], xmm3
	movaps  [rsp + nb114_fjyO], xmm4
	movaps  [rsp + nb114_fjzO], xmm5
	addps   xmm0, [rsp + nb114_fixH1]
	addps   xmm1, [rsp + nb114_fiyH1]
	addps   xmm2, [rsp + nb114_fizH1]
	movaps  [rsp + nb114_fixH1], xmm0
	movaps  [rsp + nb114_fiyH1], xmm1
	movaps  [rsp + nb114_fizH1], xmm2	

	;# H2 Coulomb
	movaps xmm0, [rsp + nb114_rinvH2H2]
	movaps xmm4, xmm0	;# xmm4=rinv 
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [rsp + nb114_qqHH]
	movhps  xmm3, [rsp + nb114_qqMH]
	shufps  xmm3, xmm3, 193 ;# 11000001 

	mulps  xmm4, xmm3       ;# xmm4=voul
	mulps  xmm0, xmm4
	addps  xmm4, [rsp + nb114_vctot] ;# local vctot summation variable
	movaps [rsp + nb114_vctot], xmm4

	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [rsp + nb114_dxH2H2]
	mulps   xmm1, [rsp + nb114_dyH2H2]
	mulps   xmm2, [rsp + nb114_dzH2H2]
	;# update forces H2 - j water 
	movaps  xmm3, [rsp + nb114_fjxO]
	movaps  xmm4, [rsp + nb114_fjyO]
	movaps  xmm5, [rsp + nb114_fjzO]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  [rsp + nb114_fjxO], xmm3
	movaps  [rsp + nb114_fjyO], xmm4
	movaps  [rsp + nb114_fjzO], xmm5
	addps   xmm0, [rsp + nb114_fixH2]
	addps   xmm1, [rsp + nb114_fiyH2]
	addps   xmm2, [rsp + nb114_fizH2]
	movaps  [rsp + nb114_fixH2], xmm0
	movaps  [rsp + nb114_fiyH2], xmm1
	movaps  [rsp + nb114_fizH2], xmm2			
	
	mov     rsi, [rbp + nb114_faction]
	;# update j water forces from local variables.
	;# transpose back first
	movaps  xmm0, [rsp + nb114_fjxO] ;# Ox H1x H2x Mx 
	movaps  xmm1, [rsp + nb114_fjyO] ;# Oy H1y H2y My
	movaps  xmm2, [rsp + nb114_fjzO] ;# Oz H1z H2z Mz
	 
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

	dec dword ptr [rsp + nb114_innerk]
	jz    .nb114_updateouterdata
	jmp   .nb114_single_loop
.nb114_updateouterdata:
	mov   ecx, [rsp + nb114_ii3]
	mov   rdi, [rbp + nb114_faction]
	mov   rsi, [rbp + nb114_fshift]
	mov   edx, [rsp + nb114_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb114_fixO]
	movaps xmm1, [rsp + nb114_fiyO] 
	movaps xmm2, [rsp + nb114_fizO]

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
	movaps xmm0, [rsp + nb114_fixH1]
	movaps xmm1, [rsp + nb114_fiyH1]
	movaps xmm2, [rsp + nb114_fizH1]

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
	movaps xmm0, [rsp + nb114_fixH2]
	movaps xmm1, [rsp + nb114_fiyH2]
	movaps xmm2, [rsp + nb114_fizH2]

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
	movaps xmm0, [rsp + nb114_fixM]
	movaps xmm1, [rsp + nb114_fiyM]
	movaps xmm2, [rsp + nb114_fizM]

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
	mov esi, [rsp + nb114_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb114_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb114_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb114_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb114_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb114_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb114_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb114_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb114_n], esi
        jmp .nb114_outer
.nb114_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb114_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb114_end
        ;# non-zero, do one more workunit
        jmp   .nb114_threadloop
.nb114_end:

	mov eax, [rsp + nb114_nouter]
	mov ebx, [rsp + nb114_ninner]
	mov rcx, [rbp + nb114_outeriter]
	mov rdx, [rbp + nb114_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 1872
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



	
.globl nb_kernel114nf_x86_64_sse 
.globl _nb_kernel114nf_x86_64_sse
nb_kernel114nf_x86_64_sse:	
_nb_kernel114nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb114nf_fshift,         16
.equiv          nb114nf_gid,            24
.equiv          nb114nf_pos,            32
.equiv          nb114nf_faction,        40
.equiv          nb114nf_charge,         48
.equiv          nb114nf_p_facel,        56
.equiv          nb114nf_argkrf,         64
.equiv          nb114nf_argcrf,         72
.equiv          nb114nf_Vc,             80
.equiv          nb114nf_type,           88
.equiv          nb114nf_p_ntype,        96
.equiv          nb114nf_vdwparam,       104
.equiv          nb114nf_Vvdw,           112
.equiv          nb114nf_p_tabscale,     120
.equiv          nb114nf_VFtab,          128
.equiv          nb114nf_invsqrta,       136
.equiv          nb114nf_dvda,           144
.equiv          nb114nf_p_gbtabscale,   152
.equiv          nb114nf_GBtab,          160
.equiv          nb114nf_p_nthreads,     168
.equiv          nb114nf_count,          176
.equiv          nb114nf_mtx,            184
.equiv          nb114nf_outeriter,      192
.equiv          nb114nf_inneriter,      200
.equiv          nb114nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb114nf_ixO,            0
.equiv          nb114nf_iyO,            16
.equiv          nb114nf_izO,            32
.equiv          nb114nf_ixH1,           48
.equiv          nb114nf_iyH1,           64
.equiv          nb114nf_izH1,           80
.equiv          nb114nf_ixH2,           96
.equiv          nb114nf_iyH2,           112
.equiv          nb114nf_izH2,           128
.equiv          nb114nf_ixM,            144
.equiv          nb114nf_iyM,            160
.equiv          nb114nf_izM,            176
.equiv          nb114nf_jxO,            192
.equiv          nb114nf_jyO,            208
.equiv          nb114nf_jzO,            224
.equiv          nb114nf_jxH1,           240
.equiv          nb114nf_jyH1,           256
.equiv          nb114nf_jzH1,           272
.equiv          nb114nf_jxH2,           288
.equiv          nb114nf_jyH2,           304
.equiv          nb114nf_jzH2,           320
.equiv          nb114nf_jxM,            336
.equiv          nb114nf_jyM,            352
.equiv          nb114nf_jzM,            368
.equiv          nb114nf_qqMM,           384
.equiv          nb114nf_qqMH,           400
.equiv          nb114nf_qqHH,           416
.equiv          nb114nf_two,            432
.equiv          nb114nf_c6,             448
.equiv          nb114nf_c12,            464
.equiv          nb114nf_vctot,          480
.equiv          nb114nf_Vvdwtot,        496
.equiv          nb114nf_half,           512
.equiv          nb114nf_three,          528
.equiv          nb114nf_rsqOO,          544
.equiv          nb114nf_rsqH1H1,        560
.equiv          nb114nf_rsqH1H2,        576
.equiv          nb114nf_rsqH1M,         592
.equiv          nb114nf_rsqH2H1,        608
.equiv          nb114nf_rsqH2H2,        624
.equiv          nb114nf_rsqH2M,         640
.equiv          nb114nf_rsqMH1,         656
.equiv          nb114nf_rsqMH2,         672
.equiv          nb114nf_rsqMM,          688
.equiv          nb114nf_rinvsqOO,       704
.equiv          nb114nf_rinvH1H1,       720
.equiv          nb114nf_rinvH1H2,       736
.equiv          nb114nf_rinvH1M,        752
.equiv          nb114nf_rinvH2H1,       768
.equiv          nb114nf_rinvH2H2,       784
.equiv          nb114nf_rinvH2M,        800
.equiv          nb114nf_rinvMH1,        816
.equiv          nb114nf_rinvMH2,        832
.equiv          nb114nf_rinvMM,         848
.equiv          nb114nf_nri,            864
.equiv          nb114nf_iinr,           872
.equiv          nb114nf_jindex,         880
.equiv          nb114nf_jjnr,           888
.equiv          nb114nf_shift,          896
.equiv          nb114nf_shiftvec,       904
.equiv          nb114nf_facel,          912
.equiv          nb114nf_innerjjnr,      920
.equiv          nb114nf_innerk,         928
.equiv          nb114nf_is3,            932
.equiv          nb114nf_ii3,            936
.equiv          nb114nf_n,              940
.equiv          nb114nf_nn1,            944
.equiv          nb114nf_nouter,         948
.equiv          nb114nf_ninner,         952
.equiv          nb114nf_salign,         956

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
	sub rsp, 960
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
	mov [rsp + nb114nf_nouter], eax
	mov [rsp + nb114nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb114nf_nri], edi
	mov [rsp + nb114nf_iinr], rsi
	mov [rsp + nb114nf_jindex], rdx
	mov [rsp + nb114nf_jjnr], rcx
	mov [rsp + nb114nf_shift], r8
	mov [rsp + nb114nf_shiftvec], r9
	mov rsi, [rbp + nb114nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb114nf_facel], xmm0

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb114nf_half], eax
	movss xmm1, [rsp + nb114nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb114nf_half],  xmm1
	movaps [rsp + nb114nf_two],  xmm2
	movaps [rsp + nb114nf_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb114nf_iinr]   ;# rcx = pointer into iinr[]
	mov   ebx, [rcx]		;# ebx =ii 

	mov   rdx, [rbp + nb114nf_charge]
	movss xmm5, [rdx + rbx*4 + 4]	
	movss xmm3, [rdx + rbx*4 + 12]	
	movss xmm4, xmm3	
	mov rsi, [rbp + nb114nf_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb114nf_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb114nf_qqMM], xmm3
	movaps [rsp + nb114nf_qqMH], xmm4
	movaps [rsp + nb114nf_qqHH], xmm5
	
	xorps xmm0, xmm0
	mov   rdx, [rbp + nb114nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb114nf_p_ntype]
	imul  ecx, [rdi]  ;# rcx = ntia = 2*ntype*type[ii0] 
	add   rdx, rcx
	mov   rax, [rbp + nb114nf_vdwparam]
	movlps xmm0, [rax + rdx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0x55
	movaps [rsp + nb114nf_c6], xmm0
	movaps [rsp + nb114nf_c12], xmm1

.nb114nf_threadloop:
        mov   rsi, [rbp + nb114nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb114nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   rbx, 1                           ;# rbx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb114nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb114nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb114nf_n], eax
        mov [rsp + nb114nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb114nf_outerstart
        jmp .nb114nf_end
	
.nb114nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb114nf_nouter]
	mov [rsp + nb114nf_nouter], ebx

.nb114nf_outer:
	mov   rax, [rsp + nb114nf_shift]  	;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb114nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb114nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb114nf_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb114nf_pos]	;# rax = base of pos[]  
	mov   [rsp + nb114nf_ii3], ebx	
	
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
	movaps [rsp + nb114nf_ixO], xmm3
	movaps [rsp + nb114nf_iyO], xmm4
	movaps [rsp + nb114nf_izO], xmm5
	movaps [rsp + nb114nf_ixH1], xmm6
	movaps [rsp + nb114nf_iyH1], xmm7

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
	movaps [rsp + nb114nf_izH1], xmm6
	movaps [rsp + nb114nf_ixH2], xmm0
	movaps [rsp + nb114nf_iyH2], xmm1
	movaps [rsp + nb114nf_izH2], xmm2
	movaps [rsp + nb114nf_ixM], xmm3
	movaps [rsp + nb114nf_iyM], xmm4
	movaps [rsp + nb114nf_izM], xmm5

	;# clear vctot 
	xorps xmm4, xmm4
	movaps [rsp + nb114nf_vctot], xmm4
	movaps [rsp + nb114nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb114nf_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb114nf_pos]
	mov   rdi, [rbp + nb114nf_faction]	
	mov   rax, [rsp + nb114nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb114nf_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb114nf_ninner]
	mov   [rsp + nb114nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb114nf_innerk], edx	;# number of innerloop atoms 
	jge   .nb114nf_unroll_loop
	jmp   .nb114nf_single_check
.nb114nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb114nf_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb114nf_innerjjnr],  16 ;# advance pointer (unroll 4) 

	mov rsi, [rbp + nb114nf_pos]   	;# base of pos[] 

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
	movaps [rsp + nb114nf_jxO], xmm0
	movaps [rsp + nb114nf_jyO], xmm1
	movaps [rsp + nb114nf_jzO], xmm2
	movaps [rsp + nb114nf_jxH1], xmm3

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
	movaps [rsp + nb114nf_jyH1], xmm0
	movaps [rsp + nb114nf_jzH1], xmm1
	movaps [rsp + nb114nf_jxH2], xmm2
	movaps [rsp + nb114nf_jyH2], xmm3

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
	movaps [rsp + nb114nf_jzH2], xmm0
	movaps [rsp + nb114nf_jxM], xmm1
	movaps [rsp + nb114nf_jyM], xmm2
	movaps [rsp + nb114nf_jzM], xmm3

	;# start calculating pairwise distances
	movaps xmm0, [rsp + nb114nf_ixO]
	movaps xmm1, [rsp + nb114nf_iyO]
	movaps xmm2, [rsp + nb114nf_izO]
	movaps xmm3, [rsp + nb114nf_ixH1]
	movaps xmm4, [rsp + nb114nf_iyH1]
	movaps xmm5, [rsp + nb114nf_izH1]
	subps  xmm0, [rsp + nb114nf_jxO]
	subps  xmm1, [rsp + nb114nf_jyO]
	subps  xmm2, [rsp + nb114nf_jzO]
	subps  xmm3, [rsp + nb114nf_jxH1]
	subps  xmm4, [rsp + nb114nf_jyH1]
	subps  xmm5, [rsp + nb114nf_jzH1]
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
	movaps [rsp + nb114nf_rsqOO], xmm0
	movaps [rsp + nb114nf_rsqH1H1], xmm3

	movaps xmm0, [rsp + nb114nf_ixH1]
	movaps xmm1, [rsp + nb114nf_iyH1]
	movaps xmm2, [rsp + nb114nf_izH1]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [rsp + nb114nf_jxH2]
	subps  xmm1, [rsp + nb114nf_jyH2]
	subps  xmm2, [rsp + nb114nf_jzH2]
	subps  xmm3, [rsp + nb114nf_jxM]
	subps  xmm4, [rsp + nb114nf_jyM]
	subps  xmm5, [rsp + nb114nf_jzM]
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
	movaps [rsp + nb114nf_rsqH1H2], xmm0
	movaps [rsp + nb114nf_rsqH1M], xmm3

	movaps xmm0, [rsp + nb114nf_ixH2]
	movaps xmm1, [rsp + nb114nf_iyH2]
	movaps xmm2, [rsp + nb114nf_izH2]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [rsp + nb114nf_jxH1]
	subps  xmm1, [rsp + nb114nf_jyH1]
	subps  xmm2, [rsp + nb114nf_jzH1]
	subps  xmm3, [rsp + nb114nf_jxH2]
	subps  xmm4, [rsp + nb114nf_jyH2]
	subps  xmm5, [rsp + nb114nf_jzH2]
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
	movaps [rsp + nb114nf_rsqH2H1], xmm0
	movaps [rsp + nb114nf_rsqH2H2], xmm3

	movaps xmm0, [rsp + nb114nf_ixH2]
	movaps xmm1, [rsp + nb114nf_iyH2]
	movaps xmm2, [rsp + nb114nf_izH2]
	movaps xmm3, [rsp + nb114nf_ixM]
	movaps xmm4, [rsp + nb114nf_iyM]
	movaps xmm5, [rsp + nb114nf_izM]
	subps  xmm0, [rsp + nb114nf_jxM]
	subps  xmm1, [rsp + nb114nf_jyM]
	subps  xmm2, [rsp + nb114nf_jzM]
	subps  xmm3, [rsp + nb114nf_jxH1]
	subps  xmm4, [rsp + nb114nf_jyH1]
	subps  xmm5, [rsp + nb114nf_jzH1]
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
	movaps [rsp + nb114nf_rsqH2M], xmm0
	movaps [rsp + nb114nf_rsqMH1], xmm4

	movaps xmm0, [rsp + nb114nf_ixM]
	movaps xmm1, [rsp + nb114nf_iyM]
	movaps xmm2, [rsp + nb114nf_izM]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [rsp + nb114nf_jxH2]
	subps  xmm1, [rsp + nb114nf_jyH2]
	subps  xmm2, [rsp + nb114nf_jzH2]
	subps  xmm3, [rsp + nb114nf_jxM]
	subps  xmm4, [rsp + nb114nf_jyM]
	subps  xmm5, [rsp + nb114nf_jzM]
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
	movaps [rsp + nb114nf_rsqMH2], xmm0
	movaps [rsp + nb114nf_rsqMM], xmm4

	;# start by doing reciprocal for OO
	movaps  xmm7, [rsp + nb114nf_rsqOO]
	rcpps   xmm2, xmm7
	movaps  xmm1, [rsp + nb114nf_two]
	mulps   xmm7, xmm2
	subps   xmm1, xmm7
	mulps   xmm2, xmm1 ;# rinvsq 
	movaps [rsp + nb114nf_rinvsqOO], xmm2
	
	;# next step is invsqrt - do two at a time.
	rsqrtps xmm1, [rsp + nb114nf_rsqH1H1]
	rsqrtps xmm5, [rsp + nb114nf_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb114nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb114nf_rsqH1H1]
	mulps   xmm5, [rsp + nb114nf_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb114nf_half] ;# rinvH1H1 
	mulps   xmm7, [rsp + nb114nf_half] ;# rinvH1H2 
	movaps  [rsp + nb114nf_rinvH1H1], xmm3
	movaps  [rsp + nb114nf_rinvH1H2], xmm7
	
	rsqrtps xmm1, [rsp + nb114nf_rsqH1M]
	rsqrtps xmm5, [rsp + nb114nf_rsqH2H1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb114nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb114nf_rsqH1M]
	mulps   xmm5, [rsp + nb114nf_rsqH2H1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb114nf_half] 
	mulps   xmm7, [rsp + nb114nf_half]
	movaps  [rsp + nb114nf_rinvH1M], xmm3
	movaps  [rsp + nb114nf_rinvH2H1], xmm7
	
	rsqrtps xmm1, [rsp + nb114nf_rsqH2H2]
	rsqrtps xmm5, [rsp + nb114nf_rsqH2M]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb114nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb114nf_rsqH2H2]
	mulps   xmm5, [rsp + nb114nf_rsqH2M]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb114nf_half] 
	mulps   xmm7, [rsp + nb114nf_half]
	movaps  [rsp + nb114nf_rinvH2H2], xmm3
	movaps  [rsp + nb114nf_rinvH2M], xmm7
	
	rsqrtps xmm1, [rsp + nb114nf_rsqMH1]
	rsqrtps xmm5, [rsp + nb114nf_rsqMH2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb114nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb114nf_rsqMH1]
	mulps   xmm5, [rsp + nb114nf_rsqMH2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb114nf_half] 
	mulps   xmm7, [rsp + nb114nf_half]
	movaps  [rsp + nb114nf_rinvMH1], xmm3
	movaps  [rsp + nb114nf_rinvMH2], xmm7
        		
	rsqrtps xmm1, [rsp + nb114nf_rsqMM]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb114nf_three]
	mulps   xmm1, [rsp + nb114nf_rsqMM]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + nb114nf_half] 
	movaps  [rsp + nb114nf_rinvMM], xmm3
	
	;# start with OO LJ interaction
	movaps xmm0, [rsp + nb114nf_rinvsqOO]
	movaps xmm1, xmm0
	mulps  xmm1, xmm1	;# rinv4
	mulps  xmm1, xmm0	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm1, [rsp + nb114nf_c6]
	mulps  xmm2, [rsp + nb114nf_c12]
	movaps xmm4, xmm2
	subps  xmm4, xmm1
	addps  xmm4, [rsp + nb114nf_Vvdwtot]
	movaps [rsp + nb114nf_Vvdwtot], xmm4

	;# Coulomb interactions 
	;# all H-H interactions
	movaps xmm0, [rsp + nb114nf_rinvH1H1]
	addps  xmm0, [rsp + nb114nf_rinvH1H2]
	addps  xmm0, [rsp + nb114nf_rinvH2H1]
	addps  xmm0, [rsp + nb114nf_rinvH2H2]
	mulps  xmm0, [rsp + nb114nf_qqHH]
	;# all M-H interactions
	movaps xmm1, [rsp + nb114nf_rinvH1M]
	addps  xmm1, [rsp + nb114nf_rinvH2M]
	addps  xmm1, [rsp + nb114nf_rinvMH1]
	addps  xmm1, [rsp + nb114nf_rinvMH2]
	mulps  xmm1, [rsp + nb114nf_qqMH]
	;# The M-M interaction
	movaps xmm2, [rsp + nb114nf_rinvMM]
	mulps  xmm2, [rsp + nb114nf_qqMM]
	addps  xmm0, xmm1
	addps  xmm2, [rsp + nb114nf_vctot] 
	addps  xmm0, xmm2
	movaps [rsp + nb114nf_vctot], xmm0 
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb114nf_innerk],  4
	jl    .nb114nf_single_check
	jmp   .nb114nf_unroll_loop
.nb114nf_single_check:
	add dword ptr [rsp + nb114nf_innerk],  4
	jnz   .nb114nf_single_loop
	jmp   .nb114nf_updateouterdata
.nb114nf_single_loop:
	mov   rdx, [rsp + nb114nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb114nf_innerjjnr],  4	

	mov rsi, [rbp + nb114nf_pos]
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
	movaps [rsp + nb114nf_jxO], xmm6
	movaps [rsp + nb114nf_jyO], xmm3
	movaps [rsp + nb114nf_jzO], xmm1

	;# do O and M in parallel
	movaps xmm0, [rsp + nb114nf_ixO]
	movaps xmm1, [rsp + nb114nf_iyO]
	movaps xmm2, [rsp + nb114nf_izO]
	movaps xmm3, [rsp + nb114nf_ixM]
	movaps xmm4, [rsp + nb114nf_iyM]
	movaps xmm5, [rsp + nb114nf_izM]
	subps  xmm0, [rsp + nb114nf_jxO]
	subps  xmm1, [rsp + nb114nf_jyO]
	subps  xmm2, [rsp + nb114nf_jzO]
	subps  xmm3, [rsp + nb114nf_jxO]
	subps  xmm4, [rsp + nb114nf_jyO]
	subps  xmm5, [rsp + nb114nf_jzO]
	
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
	movaps [rsp + nb114nf_rsqMM], xmm4
	
	;# do 1/x for O and 1/sqrt(x) for M
	rcpss  xmm1, xmm0
	rsqrtps xmm5, xmm4
	movss  xmm2, [rsp + nb114nf_two]
	movaps  xmm6, xmm5	
	mulss  xmm0, xmm1
	mulps   xmm5, xmm5
	subss  xmm2, xmm0
	movaps  xmm7, [rsp + nb114nf_three]
	mulss  xmm2, xmm1 	;# 1/r2
	
	mulps   xmm5, xmm4
	movss  xmm0, xmm2
	subps   xmm7, xmm5
	mulss  xmm2, xmm2
	mulps   xmm7, xmm6
	mulss  xmm2, xmm0 	;# 1/r6
	mulps   xmm7, [rsp + nb114nf_half] ;# rinv iH1 - j water 
	movss  xmm1, xmm2
	movaps [rsp + nb114nf_rinvMM], xmm7

	mulss  xmm2, xmm2 	;# 1/r12
	mulss  xmm1, [rsp + nb114nf_c6]
	mulss  xmm2, [rsp + nb114nf_c12]
	movss  xmm3, xmm2
	subss  xmm3, xmm1
	addss  xmm3, [rsp + nb114nf_Vvdwtot]
	movss  [rsp + nb114nf_Vvdwtot], xmm3

	;# do  M coulomb interaction
	movaps xmm0, [rsp + nb114nf_rinvMM]
	movaps xmm4, xmm0	;# xmm4=rinv

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [rsp + nb114nf_qqMH]
	movhps  xmm3, [rsp + nb114nf_qqMM]
	shufps  xmm3, xmm3, 193 ;# 11000001 

	mulps  xmm4, xmm3 	;# xmm4=voul	
	addps  xmm4, [rsp + nb114nf_vctot] 
	movaps [rsp + nb114nf_vctot], xmm4

	;# i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [rsp + nb114nf_ixH1]
	movaps  xmm1, [rsp + nb114nf_iyH1]
	movaps  xmm2, [rsp + nb114nf_izH1]	
	movaps  xmm3, [rsp + nb114nf_ixH2] 
	movaps  xmm4, [rsp + nb114nf_iyH2] 
	movaps  xmm5, [rsp + nb114nf_izH2] 
	subps   xmm0, [rsp + nb114nf_jxO]
	subps   xmm1, [rsp + nb114nf_jyO]
	subps   xmm2, [rsp + nb114nf_jzO]
	subps   xmm3, [rsp + nb114nf_jxO]
	subps   xmm4, [rsp + nb114nf_jyO]
	subps   xmm5, [rsp + nb114nf_jzO]
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
	movaps  [rsp + nb114nf_rsqH1H1], xmm0
	movaps  [rsp + nb114nf_rsqH2H2], xmm4
	
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb114nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb114nf_half] ;# rinvH1H1
	mulps   xmm7, [rsp + nb114nf_half] ;# rinvH2H2
	movaps  [rsp + nb114nf_rinvH1H1], xmm3
	movaps  [rsp + nb114nf_rinvH2H2], xmm7
	
	;# Do H1 & H2 coulomb interaction
	movaps xmm0, [rsp + nb114nf_rinvH1H1]
        addps  xmm0, [rsp + nb114nf_rinvH2H2]

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [rsp + nb114nf_qqHH]
	movhps  xmm3, [rsp + nb114nf_qqMH]
	shufps  xmm3, xmm3, 193 ;# 11000001 

	mulps  xmm0, xmm3 	;# xmm4=voul
	addps  xmm0, [rsp + nb114nf_vctot] 
	movaps [rsp + nb114nf_vctot], xmm0

	dec dword ptr [rsp + nb114nf_innerk]
	jz    .nb114nf_updateouterdata
	jmp   .nb114nf_single_loop
.nb114nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb114nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb114nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb114nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb114nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb114nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb114nf_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb114nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb114nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb114nf_n], esi
        jmp .nb114nf_outer
.nb114nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb114nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb114nf_end
        ;# non-zero, do one more workunit
        jmp   .nb114nf_threadloop
.nb114nf_end:

	mov eax, [rsp + nb114nf_nouter]
	mov ebx, [rsp + nb114nf_ninner]
	mov rcx, [rbp + nb114nf_outeriter]
	mov rdx, [rbp + nb114nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 960
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
