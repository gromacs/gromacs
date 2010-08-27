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

	
.globl nb_kernel314_x86_64_sse2
.globl _nb_kernel314_x86_64_sse2
nb_kernel314_x86_64_sse2:	
_nb_kernel314_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb314_fshift,           16
.equiv          nb314_gid,              24
.equiv          nb314_pos,              32
.equiv          nb314_faction,          40
.equiv          nb314_charge,           48
.equiv          nb314_p_facel,          56
.equiv          nb314_argkrf,           64
.equiv          nb314_argcrf,           72
.equiv          nb314_Vc,               80
.equiv          nb314_type,             88
.equiv          nb314_p_ntype,          96
.equiv          nb314_vdwparam,         104
.equiv          nb314_Vvdw,             112
.equiv          nb314_p_tabscale,       120
.equiv          nb314_VFtab,            128
.equiv          nb314_invsqrta,         136
.equiv          nb314_dvda,             144
.equiv          nb314_p_gbtabscale,     152
.equiv          nb314_GBtab,            160
.equiv          nb314_p_nthreads,       168
.equiv          nb314_count,            176
.equiv          nb314_mtx,              184
.equiv          nb314_outeriter,        192
.equiv          nb314_inneriter,        200
.equiv          nb314_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
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
.equiv          nb314_vctot,            976
.equiv          nb314_Vvdwtot,          992
.equiv          nb314_fixO,             1008
.equiv          nb314_fiyO,             1024
.equiv          nb314_fizO,             1040
.equiv          nb314_fixH1,            1056
.equiv          nb314_fiyH1,            1072
.equiv          nb314_fizH1,            1088
.equiv          nb314_fixH2,            1104
.equiv          nb314_fiyH2,            1120
.equiv          nb314_fizH2,            1136
.equiv          nb314_fixM,             1152
.equiv          nb314_fiyM,             1168
.equiv          nb314_fizM,             1184
.equiv          nb314_fjxO,             1200
.equiv          nb314_fjyO,             1216
.equiv          nb314_fjzO,             1232
.equiv          nb314_fjxH1,            1248
.equiv          nb314_fjyH1,            1264
.equiv          nb314_fjzH1,            1280
.equiv          nb314_fjxH2,            1296
.equiv          nb314_fjyH2,            1312
.equiv          nb314_fjzH2,            1328
.equiv          nb314_epsH1,            1344
.equiv          nb314_epsH2,            1360
.equiv          nb314_epsM,             1376
.equiv          nb314_half,             1392
.equiv          nb314_three,            1408
.equiv          nb314_six,              1424
.equiv          nb314_twelve,           1440
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
.equiv          nb314_is3,              1776
.equiv          nb314_ii3,              1780
.equiv          nb314_nri,              1784
.equiv          nb314_iinr,             1792
.equiv          nb314_jindex,           1800
.equiv          nb314_jjnr,             1808
.equiv          nb314_shift,            1816
.equiv          nb314_shiftvec,         1824
.equiv          nb314_facel,            1832
.equiv          nb314_innerjjnr,        1840
.equiv          nb314_innerk,           1848
.equiv          nb314_n,                1852
.equiv          nb314_nn1,              1856
.equiv          nb314_nouter,           1860
.equiv          nb314_ninner,           1864

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

	emms
	sub rsp, 1872		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb314_nouter], eax
	mov [rsp + nb314_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb314_nri], edi
	mov [rsp + nb314_iinr], rsi
	mov [rsp + nb314_jindex], rdx
	mov [rsp + nb314_jjnr], rcx
	mov [rsp + nb314_shift], r8
	mov [rsp + nb314_shiftvec], r9
	mov rsi, [rbp + nb314_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb314_facel], xmm0

	mov rax, [rbp + nb314_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb314_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb314_half], eax
	mov [rsp + nb314_half+4], ebx
	movsd xmm1, [rsp + nb314_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd xmm4, xmm3
	addpd  xmm4, xmm4       ;# six
	movapd xmm5, xmm4
	addpd  xmm5, xmm5       ;# twelve
	movapd [rsp + nb314_half], xmm1
	movapd [rsp + nb314_two], xmm2
	movapd [rsp + nb314_three], xmm3
	movapd [rsp + nb314_six], xmm4
	movapd [rsp + nb314_twelve], xmm5

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb314_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb314_charge]
	movsd xmm3, [rdx + rbx*8 + 24]	
	movsd xmm4, xmm3	
	movsd xmm5, [rdx + rbx*8 + 8]	
	mov rsi, [rbp + nb314_p_facel]
	movsd xmm0, [rsi]
	movsd xmm6, [rsp + nb314_facel]
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm5
	mulsd  xmm5, xmm5
	mulsd  xmm3, xmm6
	mulsd  xmm4, xmm6
	mulsd  xmm5, xmm6
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb314_qqMM], xmm3
	movapd [rsp + nb314_qqMH], xmm4
	movapd [rsp + nb314_qqHH], xmm5
	
	xorpd xmm0, xmm0
	mov   rdx, [rbp + nb314_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb314_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   rax, [rbp + nb314_vdwparam]
	movlpd xmm0, [rax + rdx*8]
	movlpd xmm1, [rax + rdx*8 + 8]
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	movapd [rsp + nb314_c6], xmm0
	movapd [rsp + nb314_c12], xmm1

.nb314_threadloop:
        mov   rsi, [rbp + nb314_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb314_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb314_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb314_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb314_n], eax
        mov [rsp + nb314_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb314_outerstart
        jmp .nb314_end

.nb314_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb314_nouter]
	mov [rsp + nb314_nouter], ebx

.nb314_outer:
	mov   rax, [rsp + nb314_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb314_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb314_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb314_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb314_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb314_ii3], ebx		

	addsd xmm3, [rax + rbx*8] 	;# ox
	addsd xmm4, [rax + rbx*8 + 8] 	;# oy
	addsd xmm5, [rax + rbx*8 + 16]	;# oz	
	addsd xmm6, [rax + rbx*8 + 24] 	;# h1x
	addsd xmm7, [rax + rbx*8 + 32] 	;# h1y
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	shufpd xmm7, xmm7, 0
	movapd [rsp + nb314_ixO], xmm3
	movapd [rsp + nb314_iyO], xmm4
	movapd [rsp + nb314_izO], xmm5
	movapd [rsp + nb314_ixH1], xmm6
	movapd [rsp + nb314_iyH1], xmm7

	movsd xmm6, xmm2
	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm6, [rax + rbx*8 + 40] ;# h1z
	addsd xmm0, [rax + rbx*8 + 48] ;# h2x
	addsd xmm1, [rax + rbx*8 + 56] ;# h2y
	addsd xmm2, [rax + rbx*8 + 64] ;# h2z
	addsd xmm3, [rax + rbx*8 + 72] ;# mx
	addsd xmm4, [rax + rbx*8 + 80] ;# my
	addsd xmm5, [rax + rbx*8 + 88] ;# mz

	shufpd xmm6, xmm6, 0
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb314_izH1], xmm6
	movapd [rsp + nb314_ixH2], xmm0
	movapd [rsp + nb314_iyH2], xmm1
	movapd [rsp + nb314_izH2], xmm2
	movapd [rsp + nb314_ixM], xmm3
	movapd [rsp + nb314_iyM], xmm4
	movapd [rsp + nb314_izM], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb314_vctot], xmm4
	movapd [rsp + nb314_Vvdwtot], xmm4
	movapd [rsp + nb314_fixO], xmm4
	movapd [rsp + nb314_fiyO], xmm4
	movapd [rsp + nb314_fizO], xmm4
	movapd [rsp + nb314_fixH1], xmm4
	movapd [rsp + nb314_fiyH1], xmm4
	movapd [rsp + nb314_fizH1], xmm4
	movapd [rsp + nb314_fixH2], xmm4
	movapd [rsp + nb314_fiyH2], xmm4
	movapd [rsp + nb314_fizH2], xmm4
	movapd [rsp + nb314_fixM], xmm4
	movapd [rsp + nb314_fiyM], xmm4
	movapd [rsp + nb314_fizM], xmm4
	
	mov   rax, [rsp + nb314_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb314_pos] 
	mov   rdi, [rbp + nb314_faction]	
	mov   rax, [rsp + nb314_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb314_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb314_ninner]
	mov   [rsp + nb314_ninner], ecx
	add   edx, 0
	mov   [rsp + nb314_innerk], edx    ;# number of innerloop atoms 
	jge   .nb314_unroll_loop
	jmp   .nb314_checksingle
.nb314_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb314_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	
	add qword ptr [rsp + nb314_innerjjnr], 8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb314_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	
	;# load j O coordinates
    movlpd xmm4, [rsi + rax*8] 
    movlpd xmm5, [rsi + rax*8 + 8] 
    movlpd xmm6, [rsi + rax*8 + 16] 
    movhpd xmm4, [rsi + rbx*8] 
    movhpd xmm5, [rsi + rbx*8 + 8] 
    movhpd xmm6, [rsi + rbx*8 + 16] 

    ;# xmm4 = Ox
    ;# xmm5 = Oy
    ;# xmm6 = Oz
        
    subpd xmm4, [rsp + nb314_ixO]
    subpd xmm5, [rsp + nb314_iyO]
    subpd xmm6, [rsp + nb314_izO]

    ;# store dx/dy/dz
    movapd xmm9, xmm4
    movapd xmm10, xmm5
    movapd xmm11, xmm6
    
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6	;# rsq in xmm4 

	cvtpd2ps xmm6, xmm4	
	rcpps xmm6, xmm6
	cvtps2pd xmm6, xmm6	;# lu in low xmm6 
	
	;# 1/x lookup seed in xmm6 
	movapd xmm0, [rsp + nb314_two]
	movapd xmm5, xmm4
	mulpd xmm4, xmm6	;# lu*rsq 
	subpd xmm0, xmm4	;# 2-lu*rsq 
	mulpd xmm6, xmm0	;# (new lu) 
	
	movapd xmm0, [rsp + nb314_two]
	mulpd xmm5, xmm6	;# lu*rsq 
	subpd xmm0, xmm5	;# 2-lu*rsq 
	mulpd xmm0, xmm6	;# xmm0=rinvsq 

	movapd xmm1, xmm0
	mulpd  xmm1, xmm0
	mulpd  xmm1, xmm0	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulpd  xmm2, xmm2	;# xmm2=rinvtwelve 

	mulpd  xmm1, [rsp + nb314_c6]   ;# mult by c6
	mulpd  xmm2, [rsp + nb314_c12]   ;# mult by c12
	movapd xmm5, xmm2
	subpd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	mulpd  xmm1, [rsp + nb314_six]
	mulpd  xmm2, [rsp + nb314_twelve]
	subpd  xmm2, xmm1
	mulpd  xmm0, xmm2	;# xmm0=total fscal 

    ;# increment potential
    addpd  xmm5, [rsp + nb314_Vvdwtot]
    movapd [rsp + nb314_Vvdwtot], xmm5

	mov    rdi, [rbp + nb314_faction]
	mulpd  xmm9,  xmm0
	mulpd  xmm10, xmm0
	mulpd  xmm11, xmm0

    movapd xmm0, [rsp + nb314_fixO]
    movapd xmm1, [rsp + nb314_fiyO]
    movapd xmm2, [rsp + nb314_fizO]
    
    ;# accumulate i forces
    addpd xmm0, xmm9
    addpd xmm1, xmm10
    addpd xmm2, xmm11
    movapd [rsp + nb314_fixO], xmm0
    movapd [rsp + nb314_fiyO], xmm1
    movapd [rsp + nb314_fizO], xmm2
    
	;# the fj's - start by accumulating forces from memory 
    mov rdi, [rbp + nb314_faction]
	movlpd xmm3, [rdi + rax*8]
	movlpd xmm4, [rdi + rax*8 + 8]
	movlpd xmm5, [rdi + rax*8 + 16]
	movhpd xmm3, [rdi + rbx*8]
	movhpd xmm4, [rdi + rbx*8 + 8]
	movhpd xmm5, [rdi + rbx*8 + 16]
	addpd xmm3, xmm9
	addpd xmm4, xmm10
	addpd xmm5, xmm11
	movlpd [rdi + rax*8], xmm3
	movlpd [rdi + rax*8 + 8], xmm4
	movlpd [rdi + rax*8 + 16], xmm5
	movhpd [rdi + rbx*8], xmm3
	movhpd [rdi + rbx*8 + 8], xmm4
	movhpd [rdi + rbx*8 + 16], xmm5
    ;# done with OO interaction

    mov rsi, [rbp + nb314_pos]
    ;# move j H1 coordinates to local temp variables 
    movlpd xmm0, [rsi + rax*8 + 24] 
    movlpd xmm1, [rsi + rax*8 + 32] 
    movlpd xmm2, [rsi + rax*8 + 40] 
    movhpd xmm0, [rsi + rbx*8 + 24] 
    movhpd xmm1, [rsi + rbx*8 + 32] 
    movhpd xmm2, [rsi + rbx*8 + 40] 

    ;# xmm0 = H1x
    ;# xmm1 = H1y
    ;# xmm2 = H1z
        
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subpd xmm0, [rsp + nb314_ixH1]
    subpd xmm1, [rsp + nb314_iyH1]
    subpd xmm2, [rsp + nb314_izH1]
    subpd xmm3, [rsp + nb314_ixH2]
    subpd xmm4, [rsp + nb314_iyH2]
    subpd xmm5, [rsp + nb314_izH2]
    subpd xmm6, [rsp + nb314_ixM]
    subpd xmm7, [rsp + nb314_iyM]
    subpd xmm8, [rsp + nb314_izM]
    
	movapd [rsp + nb314_dxH1H1], xmm0
	movapd [rsp + nb314_dyH1H1], xmm1
	movapd [rsp + nb314_dzH1H1], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb314_dxH2H1], xmm3
	movapd [rsp + nb314_dyH2H1], xmm4
	movapd [rsp + nb314_dzH2H1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb314_dxMH1], xmm6
	movapd [rsp + nb314_dyMH1], xmm7
	movapd [rsp + nb314_dzMH1], xmm8
	mulpd  xmm6, xmm6
	mulpd  xmm7, xmm7
	mulpd  xmm8, xmm8
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
    addpd  xmm6, xmm7
    addpd  xmm6, xmm8

	;# start doing invsqrt for jH1 atoms
    cvtpd2ps xmm1, xmm0
    cvtpd2ps xmm4, xmm3
    cvtpd2ps xmm7, xmm6
	rsqrtps xmm1, xmm1
	rsqrtps xmm4, xmm4
    rsqrtps xmm7, xmm7
    cvtps2pd xmm1, xmm1
    cvtps2pd xmm4, xmm4
    cvtps2pd xmm7, xmm7
	
	movapd  xmm2, xmm1
	movapd  xmm5, xmm4
    movapd  xmm8, xmm7
    
	mulpd   xmm1, xmm1 ;# lu*lu
	mulpd   xmm4, xmm4 ;# lu*lu
    mulpd   xmm7, xmm7 ;# lu*lu
		
	movapd  xmm9, [rsp + nb314_three]
	movapd  xmm10, xmm9
    movapd  xmm11, xmm9

	mulpd   xmm1, xmm0 ;# rsq*lu*lu
	mulpd   xmm4, xmm3 ;# rsq*lu*lu 
    mulpd   xmm7, xmm6 ;# rsq*lu*lu
	
	subpd   xmm9, xmm1
	subpd   xmm10, xmm4
    subpd   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulpd   xmm9, xmm2
	mulpd   xmm10, xmm5
    mulpd   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movapd  xmm15, [rsp + nb314_half]
	mulpd   xmm9, xmm15  ;# first iteration for rinvH1H1 
	mulpd   xmm10, xmm15 ;# first iteration for rinvH2H1
    mulpd   xmm11, xmm15 ;# first iteration for rinvMH1	

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulpd   xmm2, xmm2 ;# lu*lu
	mulpd   xmm5, xmm5 ;# lu*lu
    mulpd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb314_three]
	movapd  xmm4, xmm1
    movapd  xmm7, xmm1

	mulpd   xmm2, xmm0 ;# rsq*lu*lu
	mulpd   xmm5, xmm3 ;# rsq*lu*lu 
    mulpd   xmm8, xmm6 ;# rsq*lu*lu
	
	subpd   xmm1, xmm2
	subpd   xmm4, xmm5
    subpd   xmm7, xmm8 ;# 3-rsq*lu*lu

	mulpd   xmm9, xmm1
	mulpd   xmm10, xmm4
    mulpd   xmm11, xmm7 ;# lu*(3-rsq*lu*lu)

	movapd  xmm15, [rsp + nb314_half]
	mulpd   xmm9, xmm15  ;#  rinvH1H1 
	mulpd   xmm10, xmm15 ;#   rinvH2H1
    mulpd   xmm11, xmm15 ;#   rinvMH1
	
	movapd  [rsp + nb314_rinvH1H1], xmm9
	movapd  [rsp + nb314_rinvH2H1], xmm10
	movapd  [rsp + nb314_rinvMH1], xmm11
	
	;# H1 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movapd xmm1, [rsp + nb314_tsc]
    mulpd  xmm0, xmm9  ;# r
    mulpd  xmm3, xmm10
    mulpd  xmm6, xmm11
    mulpd  xmm0, xmm1 ;# rtab
    mulpd  xmm3, xmm1
    mulpd  xmm6, xmm1
    
    ;# truncate and convert to integers
    cvttpd2dq xmm1, xmm0
    cvttpd2dq xmm4, xmm3
    cvttpd2dq xmm7, xmm6        

    ;# convert back to float
    cvtdq2pd  xmm2, xmm1
    cvtdq2pd  xmm5, xmm4
    cvtdq2pd  xmm8, xmm7
    
    ;# multiply by 4
    pslld   xmm1, 2
    pslld   xmm4, 2
    pslld   xmm7, 2
    
    ;# move to integer registers
    pshufd xmm13, xmm1, 1
    pshufd xmm14, xmm4, 1
    pshufd xmm15, xmm7, 1
    movd    r8d, xmm1
    movd    r10d, xmm4
    movd    r12d, xmm7
    movd    r9d, xmm13
    movd    r11d, xmm14
    movd    r13d, xmm15
        
    mov  rsi, [rbp + nb314_VFtab]

    ;# calculate eps
    subpd     xmm0, xmm2
    subpd     xmm3, xmm5
    subpd     xmm6, xmm8

    movapd    [rsp + nb314_epsH1], xmm0
    movapd    [rsp + nb314_epsH2], xmm3
    movapd    [rsp + nb314_epsM], xmm6

    ;# Load LOTS of table data
    movlpd xmm0,  [rsi + r8*8]
    movlpd xmm1,  [rsi + r8*8 + 8]
    movlpd xmm2,  [rsi + r8*8 + 16]
    movlpd xmm3,  [rsi + r8*8 + 24]
    movlpd xmm4,  [rsi + r10*8]
    movlpd xmm5,  [rsi + r10*8 + 8]
    movlpd xmm6,  [rsi + r10*8 + 16]
    movlpd xmm7,  [rsi + r10*8 + 24]
    movlpd xmm8,  [rsi + r12*8]
    movlpd xmm9,  [rsi + r12*8 + 8]
    movlpd xmm10, [rsi + r12*8 + 16]
    movlpd xmm11, [rsi + r12*8 + 24]
    movhpd xmm0,  [rsi + r9*8]
    movhpd xmm1,  [rsi + r9*8 + 8]
    movhpd xmm2,  [rsi + r9*8 + 16]
    movhpd xmm3,  [rsi + r9*8 + 24]
    movhpd xmm4,  [rsi + r11*8]
    movhpd xmm5,  [rsi + r11*8 + 8]
    movhpd xmm6,  [rsi + r11*8 + 16]
    movhpd xmm7,  [rsi + r11*8 + 24]
    movhpd xmm8,  [rsi + r13*8]
    movhpd xmm9,  [rsi + r13*8 + 8]
    movhpd xmm10, [rsi + r13*8 + 16]
    movhpd xmm11, [rsi + r13*8 + 24]
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    movapd xmm12, [rsp + nb314_epsH1]
    movapd xmm13, [rsp + nb314_epsH2]
    movapd xmm14, [rsp + nb314_epsM]
    
    mulpd  xmm3, xmm12   ;# Heps
    mulpd  xmm7, xmm13
    mulpd  xmm11, xmm14 
    mulpd  xmm2, xmm12   ;# Geps
    mulpd  xmm6, xmm13
    mulpd  xmm10, xmm14 
    mulpd  xmm3, xmm12   ;# Heps2
    mulpd  xmm7, xmm13
    mulpd  xmm11, xmm14 

    addpd  xmm1, xmm2   ;# F+Geps
    addpd  xmm5, xmm6
    addpd  xmm9, xmm10 
    addpd  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addpd  xmm5, xmm7
    addpd  xmm9, xmm11 
    addpd  xmm3, xmm3    ;# 2*Heps2
    addpd  xmm7, xmm7
    addpd  xmm11, xmm11
    addpd  xmm3, xmm2    ;# 2*Heps2+Geps
    addpd  xmm7, xmm6  
    addpd  xmm11, xmm10
    addpd  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addpd  xmm7, xmm5
    addpd  xmm11, xmm9
    mulpd  xmm1, xmm12   ;# eps*Fp
    mulpd  xmm5, xmm13
    mulpd  xmm9, xmm14
    movapd xmm12, [rsp + nb314_qqHH]
    movapd xmm13, [rsp + nb314_qqMH]
    addpd  xmm1, xmm0     ;# VV
    addpd  xmm5, xmm4
    addpd  xmm9, xmm8
    mulpd  xmm1, xmm12   ;# VV*qq = vcoul
    mulpd  xmm5, xmm12
    mulpd  xmm9, xmm13
    mulpd  xmm3, xmm12    ;# FF*qq = fij
    mulpd  xmm7, xmm12
    mulpd  xmm11, xmm13
    
    ;# accumulate vctot
    addpd  xmm1, [rsp + nb314_vctot]
    addpd  xmm5, xmm9
    addpd  xmm1, xmm5
    movapd [rsp + nb314_vctot], xmm1
    
    movapd xmm10, [rsp + nb314_tsc]
    mulpd  xmm3, xmm10  ;# fscal
    mulpd  xmm7, xmm10
    mulpd  xmm10, xmm11
    
    xorpd xmm4, xmm4
    xorpd xmm8, xmm8
    xorpd xmm11, xmm11
    
    subpd xmm4, xmm3
    subpd xmm8, xmm7
    subpd xmm11, xmm10

    mulpd xmm4, [rsp + nb314_rinvH1H1]
    mulpd xmm8, [rsp + nb314_rinvH2H1]
    mulpd xmm11, [rsp + nb314_rinvMH1]
    
    ;# move j H1 forces to xmm0-xmm2
    mov rdi, [rbp + nb314_faction]
	movlpd xmm0, [rdi + rax*8 + 24]
	movlpd xmm1, [rdi + rax*8 + 32]
	movlpd xmm2, [rdi + rax*8 + 40]
	movhpd xmm0, [rdi + rbx*8 + 24]
	movhpd xmm1, [rdi + rbx*8 + 32]
	movhpd xmm2, [rdi + rbx*8 + 40]

    movapd xmm3, xmm4
    movapd xmm5, xmm4
    movapd xmm7, xmm8
    movapd xmm9, xmm8
    movapd xmm10, xmm11
    movapd xmm12, xmm11

	mulpd xmm3, [rsp + nb314_dxH1H1]
	mulpd xmm4, [rsp + nb314_dyH1H1]
	mulpd xmm5, [rsp + nb314_dzH1H1]
	mulpd xmm7, [rsp + nb314_dxH2H1]
	mulpd xmm8, [rsp + nb314_dyH2H1]
	mulpd xmm9, [rsp + nb314_dzH2H1]
	mulpd xmm10, [rsp + nb314_dxMH1]
	mulpd xmm11, [rsp + nb314_dyMH1]
	mulpd xmm12, [rsp + nb314_dzMH1]

    addpd xmm0, xmm3
    addpd xmm1, xmm4
    addpd xmm2, xmm5
    addpd xmm3, [rsp + nb314_fixH1]
    addpd xmm4, [rsp + nb314_fiyH1]
    addpd xmm5, [rsp + nb314_fizH1]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb314_fixH2]
    addpd xmm8, [rsp + nb314_fiyH2]
    addpd xmm9, [rsp + nb314_fizH2]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb314_fixM]
    addpd xmm11, [rsp + nb314_fiyM]
    addpd xmm12, [rsp + nb314_fizM]

    movapd [rsp + nb314_fixH1], xmm3
    movapd [rsp + nb314_fiyH1], xmm4
    movapd [rsp + nb314_fizH1], xmm5
    movapd [rsp + nb314_fixH2], xmm7
    movapd [rsp + nb314_fiyH2], xmm8
    movapd [rsp + nb314_fizH2], xmm9
    movapd [rsp + nb314_fixM], xmm10
    movapd [rsp + nb314_fiyM], xmm11
    movapd [rsp + nb314_fizM], xmm12
   
    ;# store back j H1 forces from xmm0-xmm2
	movlpd [rdi + rax*8 + 24], xmm0
	movlpd [rdi + rax*8 + 32], xmm1
	movlpd [rdi + rax*8 + 40], xmm2
	movhpd [rdi + rbx*8 + 24], xmm0
	movhpd [rdi + rbx*8 + 32], xmm1
	movhpd [rdi + rbx*8 + 40], xmm2

	;# move j H2 coordinates to local temp variables 
    mov rsi, [rbp + nb314_pos]
    movlpd xmm0, [rsi + rax*8 + 48] 
    movlpd xmm1, [rsi + rax*8 + 56] 
    movlpd xmm2, [rsi + rax*8 + 64] 
    movhpd xmm0, [rsi + rbx*8 + 48] 
    movhpd xmm1, [rsi + rbx*8 + 56] 
    movhpd xmm2, [rsi + rbx*8 + 64] 

    ;# xmm0 = H2x
    ;# xmm1 = H2y
    ;# xmm2 = H2z
        
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subpd xmm0, [rsp + nb314_ixH1]
    subpd xmm1, [rsp + nb314_iyH1]
    subpd xmm2, [rsp + nb314_izH1]
    subpd xmm3, [rsp + nb314_ixH2]
    subpd xmm4, [rsp + nb314_iyH2]
    subpd xmm5, [rsp + nb314_izH2]
    subpd xmm6, [rsp + nb314_ixM]
    subpd xmm7, [rsp + nb314_iyM]
    subpd xmm8, [rsp + nb314_izM]
    
	movapd [rsp + nb314_dxH1H2], xmm0
	movapd [rsp + nb314_dyH1H2], xmm1
	movapd [rsp + nb314_dzH1H2], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb314_dxH2H2], xmm3
	movapd [rsp + nb314_dyH2H2], xmm4
	movapd [rsp + nb314_dzH2H2], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb314_dxMH2], xmm6
	movapd [rsp + nb314_dyMH2], xmm7
	movapd [rsp + nb314_dzMH2], xmm8
	mulpd  xmm6, xmm6
	mulpd  xmm7, xmm7
	mulpd  xmm8, xmm8
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
    addpd  xmm6, xmm7
    addpd  xmm6, xmm8

	;# start doing invsqrt for jH2 atoms
    cvtpd2ps xmm1, xmm0
    cvtpd2ps xmm4, xmm3
    cvtpd2ps xmm7, xmm6
	rsqrtps xmm1, xmm1
	rsqrtps xmm4, xmm4
    rsqrtps xmm7, xmm7
    cvtps2pd xmm1, xmm1
    cvtps2pd xmm4, xmm4
    cvtps2pd xmm7, xmm7
	
	movapd  xmm2, xmm1
	movapd  xmm5, xmm4
    movapd  xmm8, xmm7
    
	mulpd   xmm1, xmm1 ;# lu*lu
	mulpd   xmm4, xmm4 ;# lu*lu
    mulpd   xmm7, xmm7 ;# lu*lu
		
	movapd  xmm9, [rsp + nb314_three]
	movapd  xmm10, xmm9
    movapd  xmm11, xmm9

	mulpd   xmm1, xmm0 ;# rsq*lu*lu
	mulpd   xmm4, xmm3 ;# rsq*lu*lu 
    mulpd   xmm7, xmm6 ;# rsq*lu*lu
	
	subpd   xmm9, xmm1
	subpd   xmm10, xmm4
    subpd   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulpd   xmm9, xmm2
	mulpd   xmm10, xmm5
    mulpd   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movapd  xmm15, [rsp + nb314_half]
	mulpd   xmm9, xmm15  ;# first iteration for rinvH1H2 
	mulpd   xmm10, xmm15 ;# first iteration for rinvH2H2
    mulpd   xmm11, xmm15 ;# first iteration for rinvMH2

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulpd   xmm2, xmm2 ;# lu*lu
	mulpd   xmm5, xmm5 ;# lu*lu
    mulpd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb314_three]
	movapd  xmm4, xmm1
    movapd  xmm7, xmm1

	mulpd   xmm2, xmm0 ;# rsq*lu*lu
	mulpd   xmm5, xmm3 ;# rsq*lu*lu 
    mulpd   xmm8, xmm6 ;# rsq*lu*lu
	
	subpd   xmm1, xmm2
	subpd   xmm4, xmm5
    subpd   xmm7, xmm8 ;# 3-rsq*lu*lu

	mulpd   xmm9, xmm1
	mulpd   xmm10, xmm4
    mulpd   xmm11, xmm7 ;# lu*(3-rsq*lu*lu)

	movapd  xmm15, [rsp + nb314_half]
	mulpd   xmm9, xmm15  ;#  rinvH1H2
	mulpd   xmm10, xmm15 ;#   rinvH2H2
    mulpd   xmm11, xmm15 ;#   rinvMH2
	
	movapd  [rsp + nb314_rinvH1H2], xmm9
	movapd  [rsp + nb314_rinvH2H2], xmm10
	movapd  [rsp + nb314_rinvMH2], xmm11
	
	;# H2 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movapd xmm1, [rsp + nb314_tsc]
    mulpd  xmm0, xmm9  ;# r
    mulpd  xmm3, xmm10
    mulpd  xmm6, xmm11
    mulpd  xmm0, xmm1 ;# rtab
    mulpd  xmm3, xmm1
    mulpd  xmm6, xmm1
    
    ;# truncate and convert to integers
    cvttpd2dq xmm1, xmm0
    cvttpd2dq xmm4, xmm3
    cvttpd2dq xmm7, xmm6        

    ;# convert back to float
    cvtdq2pd  xmm2, xmm1
    cvtdq2pd  xmm5, xmm4
    cvtdq2pd  xmm8, xmm7
    
    ;# multiply by 4
    pslld   xmm1, 2
    pslld   xmm4, 2
    pslld   xmm7, 2
    
    ;# move to integer registers
    pshufd xmm13, xmm1, 1
    pshufd xmm14, xmm4, 1
    pshufd xmm15, xmm7, 1
    movd    r8d, xmm1
    movd    r10d, xmm4
    movd    r12d, xmm7
    movd    r9d, xmm13
    movd    r11d, xmm14
    movd    r13d, xmm15
        
    mov  rsi, [rbp + nb314_VFtab]

    ;# calculate eps
    subpd     xmm0, xmm2
    subpd     xmm3, xmm5
    subpd     xmm6, xmm8

    movapd    [rsp + nb314_epsH1], xmm0
    movapd    [rsp + nb314_epsH2], xmm3
    movapd    [rsp + nb314_epsM], xmm6

    ;# Load LOTS of table data
    movlpd xmm0,  [rsi + r8*8]
    movlpd xmm1,  [rsi + r8*8 + 8]
    movlpd xmm2,  [rsi + r8*8 + 16]
    movlpd xmm3,  [rsi + r8*8 + 24]
    movlpd xmm4,  [rsi + r10*8]
    movlpd xmm5,  [rsi + r10*8 + 8]
    movlpd xmm6,  [rsi + r10*8 + 16]
    movlpd xmm7,  [rsi + r10*8 + 24]
    movlpd xmm8,  [rsi + r12*8]
    movlpd xmm9,  [rsi + r12*8 + 8]
    movlpd xmm10, [rsi + r12*8 + 16]
    movlpd xmm11, [rsi + r12*8 + 24]
    movhpd xmm0,  [rsi + r9*8]
    movhpd xmm1,  [rsi + r9*8 + 8]
    movhpd xmm2,  [rsi + r9*8 + 16]
    movhpd xmm3,  [rsi + r9*8 + 24]
    movhpd xmm4,  [rsi + r11*8]
    movhpd xmm5,  [rsi + r11*8 + 8]
    movhpd xmm6,  [rsi + r11*8 + 16]
    movhpd xmm7,  [rsi + r11*8 + 24]
    movhpd xmm8,  [rsi + r13*8]
    movhpd xmm9,  [rsi + r13*8 + 8]
    movhpd xmm10, [rsi + r13*8 + 16]
    movhpd xmm11, [rsi + r13*8 + 24]
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    movapd xmm12, [rsp + nb314_epsH1]
    movapd xmm13, [rsp + nb314_epsH2]
    movapd xmm14, [rsp + nb314_epsM]
    
    mulpd  xmm3, xmm12   ;# Heps
    mulpd  xmm7, xmm13
    mulpd  xmm11, xmm14 
    mulpd  xmm2, xmm12   ;# Geps
    mulpd  xmm6, xmm13
    mulpd  xmm10, xmm14 
    mulpd  xmm3, xmm12   ;# Heps2
    mulpd  xmm7, xmm13
    mulpd  xmm11, xmm14 

    addpd  xmm1, xmm2   ;# F+Geps
    addpd  xmm5, xmm6
    addpd  xmm9, xmm10 
    addpd  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addpd  xmm5, xmm7
    addpd  xmm9, xmm11 
    addpd  xmm3, xmm3    ;# 2*Heps2
    addpd  xmm7, xmm7
    addpd  xmm11, xmm11
    addpd  xmm3, xmm2    ;# 2*Heps2+Geps
    addpd  xmm7, xmm6  
    addpd  xmm11, xmm10
    addpd  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addpd  xmm7, xmm5
    addpd  xmm11, xmm9
    mulpd  xmm1, xmm12   ;# eps*Fp
    mulpd  xmm5, xmm13
    mulpd  xmm9, xmm14
    movapd xmm12, [rsp + nb314_qqHH]
    movapd xmm13, [rsp + nb314_qqMH]
    addpd  xmm1, xmm0     ;# VV
    addpd  xmm5, xmm4
    addpd  xmm9, xmm8
    mulpd  xmm1, xmm12   ;# VV*qq = vcoul
    mulpd  xmm5, xmm12
    mulpd  xmm9, xmm13
    mulpd  xmm3, xmm12    ;# FF*qq = fij
    mulpd  xmm7, xmm12
    mulpd  xmm11, xmm13
    
    ;# accumulate vctot
    addpd  xmm1, [rsp + nb314_vctot]
    addpd  xmm5, xmm9
    addpd  xmm1, xmm5
    movapd [rsp + nb314_vctot], xmm1
    
    movapd xmm10, [rsp + nb314_tsc]
    mulpd  xmm3, xmm10  ;# fscal
    mulpd  xmm7, xmm10
    mulpd  xmm10, xmm11
    
    xorpd xmm4, xmm4
    xorpd xmm8, xmm8
    xorpd xmm11, xmm11
    
    subpd xmm4, xmm3
    subpd xmm8, xmm7
    subpd xmm11, xmm10

    mulpd xmm4, [rsp + nb314_rinvH1H2]
    mulpd xmm8, [rsp + nb314_rinvH2H2]
    mulpd xmm11, [rsp + nb314_rinvMH2]
    
    ;# move j H2 forces to xmm0-xmm2
    mov rdi, [rbp + nb314_faction]
	movlpd xmm0, [rdi + rax*8 + 48]
	movlpd xmm1, [rdi + rax*8 + 56]
	movlpd xmm2, [rdi + rax*8 + 64]
	movhpd xmm0, [rdi + rbx*8 + 48]
	movhpd xmm1, [rdi + rbx*8 + 56]
	movhpd xmm2, [rdi + rbx*8 + 64]

    movapd xmm3, xmm4
    movapd xmm5, xmm4
    movapd xmm7, xmm8
    movapd xmm9, xmm8
    movapd xmm10, xmm11
    movapd xmm12, xmm11

	mulpd xmm3, [rsp + nb314_dxH1H2]
	mulpd xmm4, [rsp + nb314_dyH1H2]
	mulpd xmm5, [rsp + nb314_dzH1H2]
	mulpd xmm7, [rsp + nb314_dxH2H2]
	mulpd xmm8, [rsp + nb314_dyH2H2]
	mulpd xmm9, [rsp + nb314_dzH2H2]
	mulpd xmm10, [rsp + nb314_dxMH2]
	mulpd xmm11, [rsp + nb314_dyMH2]
	mulpd xmm12, [rsp + nb314_dzMH2]

    addpd xmm0, xmm3
    addpd xmm1, xmm4
    addpd xmm2, xmm5
    addpd xmm3, [rsp + nb314_fixH1]
    addpd xmm4, [rsp + nb314_fiyH1]
    addpd xmm5, [rsp + nb314_fizH1]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb314_fixH2]
    addpd xmm8, [rsp + nb314_fiyH2]
    addpd xmm9, [rsp + nb314_fizH2]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb314_fixM]
    addpd xmm11, [rsp + nb314_fiyM]
    addpd xmm12, [rsp + nb314_fizM]

    movapd [rsp + nb314_fixH1], xmm3
    movapd [rsp + nb314_fiyH1], xmm4
    movapd [rsp + nb314_fizH1], xmm5
    movapd [rsp + nb314_fixH2], xmm7
    movapd [rsp + nb314_fiyH2], xmm8
    movapd [rsp + nb314_fizH2], xmm9
    movapd [rsp + nb314_fixM], xmm10
    movapd [rsp + nb314_fiyM], xmm11
    movapd [rsp + nb314_fizM], xmm12
   
    ;# store back j H2 forces from xmm0-xmm2
	movlpd [rdi + rax*8 + 48], xmm0
	movlpd [rdi + rax*8 + 56], xmm1
	movlpd [rdi + rax*8 + 64], xmm2
	movhpd [rdi + rbx*8 + 48], xmm0
	movhpd [rdi + rbx*8 + 56], xmm1
	movhpd [rdi + rbx*8 + 64], xmm2
       
	;# move j M coordinates to local temp variables 
    mov rsi, [rbp + nb314_pos]
    movlpd xmm0, [rsi + rax*8 + 72] 
    movlpd xmm1, [rsi + rax*8 + 80] 
    movlpd xmm2, [rsi + rax*8 + 88] 
    movhpd xmm0, [rsi + rbx*8 + 72] 
    movhpd xmm1, [rsi + rbx*8 + 80] 
    movhpd xmm2, [rsi + rbx*8 + 88] 

    ;# xmm0 = Mx
    ;# xmm1 = My
    ;# xmm2 = Mz
        
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subpd xmm0, [rsp + nb314_ixH1]
    subpd xmm1, [rsp + nb314_iyH1]
    subpd xmm2, [rsp + nb314_izH1]
    subpd xmm3, [rsp + nb314_ixH2]
    subpd xmm4, [rsp + nb314_iyH2]
    subpd xmm5, [rsp + nb314_izH2]
    subpd xmm6, [rsp + nb314_ixM]
    subpd xmm7, [rsp + nb314_iyM]
    subpd xmm8, [rsp + nb314_izM]
    
	movapd [rsp + nb314_dxH1M], xmm0
	movapd [rsp + nb314_dyH1M], xmm1
	movapd [rsp + nb314_dzH1M], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb314_dxH2M], xmm3
	movapd [rsp + nb314_dyH2M], xmm4
	movapd [rsp + nb314_dzH2M], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb314_dxMM], xmm6
	movapd [rsp + nb314_dyMM], xmm7
	movapd [rsp + nb314_dzMM], xmm8
	mulpd  xmm6, xmm6
	mulpd  xmm7, xmm7
	mulpd  xmm8, xmm8
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
    addpd  xmm6, xmm7
    addpd  xmm6, xmm8

	;# start doing invsqrt for jM atoms
    cvtpd2ps xmm1, xmm0
    cvtpd2ps xmm4, xmm3
    cvtpd2ps xmm7, xmm6
	rsqrtps xmm1, xmm1
	rsqrtps xmm4, xmm4
    rsqrtps xmm7, xmm7
    cvtps2pd xmm1, xmm1
    cvtps2pd xmm4, xmm4
    cvtps2pd xmm7, xmm7
	
	movapd  xmm2, xmm1
	movapd  xmm5, xmm4
    movapd  xmm8, xmm7
    
	mulpd   xmm1, xmm1 ;# lu*lu
	mulpd   xmm4, xmm4 ;# lu*lu
    mulpd   xmm7, xmm7 ;# lu*lu
		
	movapd  xmm9, [rsp + nb314_three]
	movapd  xmm10, xmm9
    movapd  xmm11, xmm9

	mulpd   xmm1, xmm0 ;# rsq*lu*lu
	mulpd   xmm4, xmm3 ;# rsq*lu*lu 
    mulpd   xmm7, xmm6 ;# rsq*lu*lu
	
	subpd   xmm9, xmm1
	subpd   xmm10, xmm4
    subpd   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulpd   xmm9, xmm2
	mulpd   xmm10, xmm5
    mulpd   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movapd  xmm15, [rsp + nb314_half]
	mulpd   xmm9, xmm15  ;# first iteration for rinvH1M 
	mulpd   xmm10, xmm15 ;# first iteration for rinvH2M
    mulpd   xmm11, xmm15 ;# first iteration for rinvMM

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulpd   xmm2, xmm2 ;# lu*lu
	mulpd   xmm5, xmm5 ;# lu*lu
    mulpd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb314_three]
	movapd  xmm4, xmm1
    movapd  xmm7, xmm1

	mulpd   xmm2, xmm0 ;# rsq*lu*lu
	mulpd   xmm5, xmm3 ;# rsq*lu*lu 
    mulpd   xmm8, xmm6 ;# rsq*lu*lu
	
	subpd   xmm1, xmm2
	subpd   xmm4, xmm5
    subpd   xmm7, xmm8 ;# 3-rsq*lu*lu

	mulpd   xmm9, xmm1
	mulpd   xmm10, xmm4
    mulpd   xmm11, xmm7 ;# lu*(3-rsq*lu*lu)

	movapd  xmm15, [rsp + nb314_half]
	mulpd   xmm9, xmm15  ;#  rinvH1M
	mulpd   xmm10, xmm15 ;#   rinvH2M
    mulpd   xmm11, xmm15 ;#   rinvMM
    
	
	movapd  [rsp + nb314_rinvH1M], xmm9
	movapd  [rsp + nb314_rinvH2M], xmm10
	movapd  [rsp + nb314_rinvMM], xmm11
	
	;# M interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movapd xmm1, [rsp + nb314_tsc]
    mulpd  xmm0, xmm9  ;# r
    mulpd  xmm3, xmm10
    mulpd  xmm6, xmm11
    mulpd  xmm0, xmm1 ;# rtab
    mulpd  xmm3, xmm1
    mulpd  xmm6, xmm1
    
    ;# truncate and convert to integers
    cvttpd2dq xmm1, xmm0
    cvttpd2dq xmm4, xmm3
    cvttpd2dq xmm7, xmm6        

    ;# convert back to float
    cvtdq2pd  xmm2, xmm1
    cvtdq2pd  xmm5, xmm4
    cvtdq2pd  xmm8, xmm7
    
    ;# multiply by 4
    pslld   xmm1, 2
    pslld   xmm4, 2
    pslld   xmm7, 2
    
    ;# move to integer registers
    pshufd xmm13, xmm1, 1
    pshufd xmm14, xmm4, 1
    pshufd xmm15, xmm7, 1
    movd    r8d, xmm1
    movd    r10d, xmm4
    movd    r12d, xmm7
    movd    r9d, xmm13
    movd    r11d, xmm14
    movd    r13d, xmm15
        
    mov  rsi, [rbp + nb314_VFtab]

    ;# calculate eps
    subpd     xmm0, xmm2
    subpd     xmm3, xmm5
    subpd     xmm6, xmm8

    movapd    [rsp + nb314_epsH1], xmm0
    movapd    [rsp + nb314_epsH2], xmm3
    movapd    [rsp + nb314_epsM], xmm6

    ;# Load LOTS of table data
    movlpd xmm0,  [rsi + r8*8]
    movlpd xmm1,  [rsi + r8*8 + 8]
    movlpd xmm2,  [rsi + r8*8 + 16]
    movlpd xmm3,  [rsi + r8*8 + 24]
    movlpd xmm4,  [rsi + r10*8]
    movlpd xmm5,  [rsi + r10*8 + 8]
    movlpd xmm6,  [rsi + r10*8 + 16]
    movlpd xmm7,  [rsi + r10*8 + 24]
    movlpd xmm8,  [rsi + r12*8]
    movlpd xmm9,  [rsi + r12*8 + 8]
    movlpd xmm10, [rsi + r12*8 + 16]
    movlpd xmm11, [rsi + r12*8 + 24]
    movhpd xmm0,  [rsi + r9*8]
    movhpd xmm1,  [rsi + r9*8 + 8]
    movhpd xmm2,  [rsi + r9*8 + 16]
    movhpd xmm3,  [rsi + r9*8 + 24]
    movhpd xmm4,  [rsi + r11*8]
    movhpd xmm5,  [rsi + r11*8 + 8]
    movhpd xmm6,  [rsi + r11*8 + 16]
    movhpd xmm7,  [rsi + r11*8 + 24]
    movhpd xmm8,  [rsi + r13*8]
    movhpd xmm9,  [rsi + r13*8 + 8]
    movhpd xmm10, [rsi + r13*8 + 16]
    movhpd xmm11, [rsi + r13*8 + 24]
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    movapd xmm12, [rsp + nb314_epsH1]
    movapd xmm13, [rsp + nb314_epsH2]
    movapd xmm14, [rsp + nb314_epsM]
    
    mulpd  xmm3, xmm12   ;# Heps
    mulpd  xmm7, xmm13
    mulpd  xmm11, xmm14 
    mulpd  xmm2, xmm12   ;# Geps
    mulpd  xmm6, xmm13
    mulpd  xmm10, xmm14 
    mulpd  xmm3, xmm12   ;# Heps2
    mulpd  xmm7, xmm13
    mulpd  xmm11, xmm14 

    addpd  xmm1, xmm2   ;# F+Geps
    addpd  xmm5, xmm6
    addpd  xmm9, xmm10 
    addpd  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addpd  xmm5, xmm7
    addpd  xmm9, xmm11 
    addpd  xmm3, xmm3    ;# 2*Heps2
    addpd  xmm7, xmm7
    addpd  xmm11, xmm11
    addpd  xmm3, xmm2    ;# 2*Heps2+Geps
    addpd  xmm7, xmm6  
    addpd  xmm11, xmm10
    addpd  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addpd  xmm7, xmm5
    addpd  xmm11, xmm9
    mulpd  xmm1, xmm12   ;# eps*Fp
    mulpd  xmm5, xmm13
    mulpd  xmm9, xmm14
    movapd xmm12, [rsp + nb314_qqMH]
    movapd xmm13, [rsp + nb314_qqMM]
    addpd  xmm1, xmm0     ;# VV
    addpd  xmm5, xmm4
    addpd  xmm9, xmm8
    mulpd  xmm1, xmm12   ;# VV*qq = vcoul
    mulpd  xmm5, xmm12
    mulpd  xmm9, xmm13
    mulpd  xmm3, xmm12    ;# FF*qq = fij
    mulpd  xmm7, xmm12
    mulpd  xmm11, xmm13
    
    ;# accumulate vctot
    addpd  xmm1, [rsp + nb314_vctot]
    addpd  xmm5, xmm9
    addpd  xmm1, xmm5
    movapd [rsp + nb314_vctot], xmm1
    
    movapd xmm10, [rsp + nb314_tsc]
    mulpd  xmm3, xmm10  ;# fscal
    mulpd  xmm7, xmm10
    mulpd  xmm10, xmm11
        
    xorpd xmm4, xmm4
    xorpd xmm8, xmm8
    xorpd xmm11, xmm11
    
    subpd xmm4, xmm3
    subpd xmm8, xmm7
    subpd xmm11, xmm10

    mulpd xmm4, [rsp + nb314_rinvH1M]
    mulpd xmm8, [rsp + nb314_rinvH2M]
    mulpd xmm11, [rsp + nb314_rinvMM]
    
    ;# move j M forces to xmm0-xmm2
    mov rdi, [rbp + nb314_faction]
	movlpd xmm0, [rdi + rax*8 + 72]
	movlpd xmm1, [rdi + rax*8 + 80]
	movlpd xmm2, [rdi + rax*8 + 88]
	movhpd xmm0, [rdi + rbx*8 + 72]
	movhpd xmm1, [rdi + rbx*8 + 80]
	movhpd xmm2, [rdi + rbx*8 + 88]

    movapd xmm3, xmm4
    movapd xmm5, xmm4
    movapd xmm7, xmm8
    movapd xmm9, xmm8
    movapd xmm10, xmm11
    movapd xmm12, xmm11

	mulpd xmm3, [rsp + nb314_dxH1M]
	mulpd xmm4, [rsp + nb314_dyH1M]
	mulpd xmm5, [rsp + nb314_dzH1M]
	mulpd xmm7, [rsp + nb314_dxH2M]
	mulpd xmm8, [rsp + nb314_dyH2M]
	mulpd xmm9, [rsp + nb314_dzH2M]
	mulpd xmm10, [rsp + nb314_dxMM]
	mulpd xmm11, [rsp + nb314_dyMM]
	mulpd xmm12, [rsp + nb314_dzMM]

    addpd xmm0, xmm3
    addpd xmm1, xmm4
    addpd xmm2, xmm5
    addpd xmm3, [rsp + nb314_fixH1]
    addpd xmm4, [rsp + nb314_fiyH1]
    addpd xmm5, [rsp + nb314_fizH1]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb314_fixH2]
    addpd xmm8, [rsp + nb314_fiyH2]
    addpd xmm9, [rsp + nb314_fizH2]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb314_fixM]
    addpd xmm11, [rsp + nb314_fiyM]
    addpd xmm12, [rsp + nb314_fizM]

    movapd [rsp + nb314_fixH1], xmm3
    movapd [rsp + nb314_fiyH1], xmm4
    movapd [rsp + nb314_fizH1], xmm5
    movapd [rsp + nb314_fixH2], xmm7
    movapd [rsp + nb314_fiyH2], xmm8
    movapd [rsp + nb314_fizH2], xmm9
    movapd [rsp + nb314_fixM], xmm10
    movapd [rsp + nb314_fiyM], xmm11
    movapd [rsp + nb314_fizM], xmm12
   
    ;# store back j M forces from xmm0-xmm2
	movlpd [rdi + rax*8 + 72], xmm0
	movlpd [rdi + rax*8 + 80], xmm1
	movlpd [rdi + rax*8 + 88], xmm2
	movhpd [rdi + rbx*8 + 72], xmm0
	movhpd [rdi + rbx*8 + 80], xmm1
	movhpd [rdi + rbx*8 + 88], xmm2
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb314_innerk],  2
	jl    .nb314_checksingle
	jmp   .nb314_unroll_loop
.nb314_checksingle:
	mov   edx, [rsp + nb314_innerk]
	and   edx, 1
	jnz   .nb314_dosingle
	jmp   .nb314_updateouterdata
.nb314_dosingle:
	mov   rdx, [rsp + nb314_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	

	mov rsi, [rbp + nb314_pos]
	lea   rax, [rax + rax*2]  

	;# load j O coordinates
    movsd xmm4, [rsi + rax*8] 
    movsd xmm5, [rsi + rax*8 + 8] 
    movsd xmm6, [rsi + rax*8 + 16] 

    ;# xmm4 = Ox
    ;# xmm5 = Oy
    ;# xmm6 = Oz
        
    subsd xmm4, [rsp + nb314_ixO]
    subsd xmm5, [rsp + nb314_iyO]
    subsd xmm6, [rsp + nb314_izO]

    ;# store dx/dy/dz
    movapd xmm9, xmm4
    movapd xmm10, xmm5
    movapd xmm11, xmm6
    
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6	;# rsq in xmm4 

	cvtsd2ss xmm6, xmm4	
	rcpss xmm6, xmm6
	cvtss2sd xmm6, xmm6	;# lu in low xmm6 
	
	;# 1/x lookup seed in xmm6 
	movapd xmm0, [rsp + nb314_two]
	movapd xmm5, xmm4
	mulsd xmm4, xmm6	;# lu*rsq 
	subsd xmm0, xmm4	;# 2-lu*rsq 
	mulsd xmm6, xmm0	;# (new lu) 
	
	movapd xmm0, [rsp + nb314_two]
	mulsd xmm5, xmm6	;# lu*rsq 
	subsd xmm0, xmm5	;# 2-lu*rsq 
	mulsd xmm0, xmm6	;# xmm0=rinvsq 

	movapd xmm1, xmm0
	mulsd  xmm1, xmm0
	mulsd  xmm1, xmm0	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 

	mulsd  xmm1, [rsp + nb314_c6]   ;# mult by c6
	mulsd  xmm2, [rsp + nb314_c12]   ;# mult by c12
	movapd xmm5, xmm2
	subsd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	mulsd  xmm1, [rsp + nb314_six]
	mulsd  xmm2, [rsp + nb314_twelve]
	subsd  xmm2, xmm1
	mulsd  xmm0, xmm2	;# xmm0=total fscal 

    ;# increment potential
    addsd  xmm5, [rsp + nb314_Vvdwtot]
    movsd [rsp + nb314_Vvdwtot], xmm5

	mov    rdi, [rbp + nb314_faction]
	mulsd  xmm9,  xmm0
	mulsd  xmm10, xmm0
	mulsd  xmm11, xmm0

    movapd xmm0, [rsp + nb314_fixO]
    movapd xmm1, [rsp + nb314_fiyO]
    movapd xmm2, [rsp + nb314_fizO]
    
    ;# accumulate i forces
    addsd xmm0, xmm9
    addsd xmm1, xmm10
    addsd xmm2, xmm11
    movsd [rsp + nb314_fixO], xmm0
    movsd [rsp + nb314_fiyO], xmm1
    movsd [rsp + nb314_fizO], xmm2
    
	;# the fj's - start by accumulating forces from memory 
    mov rdi, [rbp + nb314_faction]
	addsd xmm9,  [rdi + rax*8]
	addsd xmm10, [rdi + rax*8 + 8]
	addsd xmm11, [rdi + rax*8 + 16]
	movsd [rdi + rax*8], xmm9
	movsd [rdi + rax*8 + 8], xmm10
	movsd [rdi + rax*8 + 16], xmm11
    ;# done with OO interaction

	;# move j H1 coordinates to local temp variables 
    mov rsi, [rbp + nb314_pos]
    movsd xmm0, [rsi + rax*8 + 24] 
    movsd xmm1, [rsi + rax*8 + 32] 
    movsd xmm2, [rsi + rax*8 + 40] 

    ;# xmm0 = H1x
    ;# xmm1 = H1y
    ;# xmm2 = H1z
        
    movsd xmm3, xmm0
    movsd xmm4, xmm1
    movsd xmm5, xmm2
    movsd xmm6, xmm0
    movsd xmm7, xmm1
    movsd xmm8, xmm2
    
    subsd xmm0, [rsp + nb314_ixH1]
    subsd xmm1, [rsp + nb314_iyH1]
    subsd xmm2, [rsp + nb314_izH1]
    subsd xmm3, [rsp + nb314_ixH2]
    subsd xmm4, [rsp + nb314_iyH2]
    subsd xmm5, [rsp + nb314_izH2]
    subsd xmm6, [rsp + nb314_ixM]
    subsd xmm7, [rsp + nb314_iyM]
    subsd xmm8, [rsp + nb314_izM]
    
	movsd [rsp + nb314_dxH1H1], xmm0
	movsd [rsp + nb314_dyH1H1], xmm1
	movsd [rsp + nb314_dzH1H1], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movsd [rsp + nb314_dxH2H1], xmm3
	movsd [rsp + nb314_dyH2H1], xmm4
	movsd [rsp + nb314_dzH2H1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	movsd [rsp + nb314_dxMH1], xmm6
	movsd [rsp + nb314_dyMH1], xmm7
	movsd [rsp + nb314_dzMH1], xmm8
	mulsd  xmm6, xmm6
	mulsd  xmm7, xmm7
	mulsd  xmm8, xmm8
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
    addsd  xmm6, xmm7
    addsd  xmm6, xmm8

	;# start doing invsqrt for jH1 atoms
    cvtsd2ss xmm1, xmm0
    cvtsd2ss xmm4, xmm3
    cvtsd2ss xmm7, xmm6
	rsqrtss xmm1, xmm1
	rsqrtss xmm4, xmm4
    rsqrtss xmm7, xmm7
    cvtss2sd xmm1, xmm1
    cvtss2sd xmm4, xmm4
    cvtss2sd xmm7, xmm7
	
	movsd  xmm2, xmm1
	movsd  xmm5, xmm4
    movsd  xmm8, xmm7
    
	mulsd   xmm1, xmm1 ;# lu*lu
	mulsd   xmm4, xmm4 ;# lu*lu
    mulsd   xmm7, xmm7 ;# lu*lu
		
	movsd  xmm9, [rsp + nb314_three]
	movsd  xmm10, xmm9
    movsd  xmm11, xmm9

	mulsd   xmm1, xmm0 ;# rsq*lu*lu
	mulsd   xmm4, xmm3 ;# rsq*lu*lu 
    mulsd   xmm7, xmm6 ;# rsq*lu*lu
	
	subsd   xmm9, xmm1
	subsd   xmm10, xmm4
    subsd   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm2
	mulsd   xmm10, xmm5
    mulsd   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movsd  xmm15, [rsp + nb314_half]
	mulsd   xmm9, xmm15  ;# first iteration for rinvH1H1 
	mulsd   xmm10, xmm15 ;# first iteration for rinvH2H1
    mulsd   xmm11, xmm15 ;# first iteration for rinvMH1	

    ;# second iteration step    
	movsd  xmm2, xmm9
	movsd  xmm5, xmm10
    movsd  xmm8, xmm11
    
	mulsd   xmm2, xmm2 ;# lu*lu
	mulsd   xmm5, xmm5 ;# lu*lu
    mulsd   xmm8, xmm8 ;# lu*lu
		
	movsd  xmm1, [rsp + nb314_three]
	movsd  xmm4, xmm1
    movsd  xmm7, xmm1

	mulsd   xmm2, xmm0 ;# rsq*lu*lu
	mulsd   xmm5, xmm3 ;# rsq*lu*lu 
    mulsd   xmm8, xmm6 ;# rsq*lu*lu
	
	subsd   xmm1, xmm2
	subsd   xmm4, xmm5
    subsd   xmm7, xmm8 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm1
	mulsd   xmm10, xmm4
    mulsd   xmm11, xmm7 ;# lu*(3-rsq*lu*lu)

	movsd  xmm15, [rsp + nb314_half]
	mulsd   xmm9, xmm15  ;#  rinvH1H1 
	mulsd   xmm10, xmm15 ;#   rinvH2H1
    mulsd   xmm11, xmm15 ;#   rinvMH1
		
	movsd  [rsp + nb314_rinvH1H1], xmm9
	movsd  [rsp + nb314_rinvH2H1], xmm10
	movsd  [rsp + nb314_rinvMH1], xmm11
	
	;# H1 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movsd  xmm1, [rsp + nb314_tsc]
    mulsd  xmm0, xmm9  ;# r
    mulsd  xmm3, xmm10
    mulsd  xmm6, xmm11
    mulsd  xmm0, xmm1 ;# rtab
    mulsd  xmm3, xmm1
    mulsd  xmm6, xmm1
    
    ;# truncate and convert to integers
    cvttsd2si r8d, xmm0
    cvttsd2si r10d, xmm3
    cvttsd2si r12d, xmm6        

    ;# convert back to float
    cvtsi2sd  xmm2, r8d
    cvtsi2sd  xmm5, r10d
    cvtsi2sd  xmm8, r12d
    
    ;# multiply by 4
    shl   r8, 2
    shl   r10d, 2
    shl   r12d, 2
    
    mov  rsi, [rbp + nb314_VFtab]

    ;# calculate eps
    subsd     xmm0, xmm2
    subsd     xmm3, xmm5
    subsd     xmm6, xmm8

    movsd    [rsp + nb314_epsH1], xmm0
    movsd    [rsp + nb314_epsH2], xmm3
    movsd    [rsp + nb314_epsM], xmm6

    ;# Load LOTS of table data
    movsd xmm0,  [rsi + r8*8]
    movsd xmm1,  [rsi + r8*8 + 8]
    movsd xmm2,  [rsi + r8*8 + 16]
    movsd xmm3,  [rsi + r8*8 + 24]
    movsd xmm4,  [rsi + r10*8]
    movsd xmm5,  [rsi + r10*8 + 8]
    movsd xmm6,  [rsi + r10*8 + 16]
    movsd xmm7,  [rsi + r10*8 + 24]
    movsd xmm8,  [rsi + r12*8]
    movsd xmm9,  [rsi + r12*8 + 8]
    movsd xmm10, [rsi + r12*8 + 16]
    movsd xmm11, [rsi + r12*8 + 24]
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    movsd xmm12, [rsp + nb314_epsH1]
    movsd xmm13, [rsp + nb314_epsH2]
    movsd xmm14, [rsp + nb314_epsM]
    
    mulsd  xmm3, xmm12   ;# Heps
    mulsd  xmm7, xmm13
    mulsd  xmm11, xmm14 
    mulsd  xmm2, xmm12   ;# Geps
    mulsd  xmm6, xmm13
    mulsd  xmm10, xmm14 
    mulsd  xmm3, xmm12   ;# Heps2
    mulsd  xmm7, xmm13
    mulsd  xmm11, xmm14 

    addsd  xmm1, xmm2   ;# F+Geps
    addsd  xmm5, xmm6
    addsd  xmm9, xmm10 
    addsd  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addsd  xmm5, xmm7
    addsd  xmm9, xmm11 
    addsd  xmm3, xmm3    ;# 2*Heps2
    addsd  xmm7, xmm7
    addsd  xmm11, xmm11
    addsd  xmm3, xmm2    ;# 2*Heps2+Geps
    addsd  xmm7, xmm6  
    addsd  xmm11, xmm10
    addsd  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addsd  xmm7, xmm5
    addsd  xmm11, xmm9
    mulsd  xmm1, xmm12   ;# eps*Fp
    mulsd  xmm5, xmm13
    mulsd  xmm9, xmm14
    movsd xmm12, [rsp + nb314_qqHH]
    movsd xmm13, [rsp + nb314_qqMH]
    addsd  xmm1, xmm0     ;# VV
    addsd  xmm5, xmm4
    addsd  xmm9, xmm8
    mulsd  xmm1, xmm12   ;# VV*qq = vcoul
    mulsd  xmm5, xmm12
    mulsd  xmm9, xmm13
    mulsd  xmm3, xmm12    ;# FF*qq = fij
    mulsd  xmm7, xmm12
    mulsd  xmm11, xmm13
    
    ;# accumulate vctot
    addsd  xmm1, [rsp + nb314_vctot]
    addsd  xmm5, xmm9
    addsd  xmm1, xmm5
    movsd [rsp + nb314_vctot], xmm1
    
    movsd xmm10, [rsp + nb314_tsc]
    mulsd  xmm3, xmm10  ;# fscal
    mulsd  xmm7, xmm10
    mulsd  xmm10, xmm11
    
    xorpd xmm4, xmm4
    xorpd xmm8, xmm8
    xorpd xmm11, xmm11
    
    subsd xmm4, xmm3
    subsd xmm8, xmm7
    subsd xmm11, xmm10

    mulsd xmm4, [rsp + nb314_rinvH1H1]
    mulsd xmm8, [rsp + nb314_rinvH2H1]
    mulsd xmm11, [rsp + nb314_rinvMH1]
    
    ;# move j H1 forces to xmm0-xmm2
    mov rdi, [rbp + nb314_faction]
	movsd xmm0, [rdi + rax*8 + 24]
	movsd xmm1, [rdi + rax*8 + 32]
	movsd xmm2, [rdi + rax*8 + 40]

    movapd xmm3, xmm4
    movapd xmm5, xmm4
    movapd xmm7, xmm8
    movapd xmm9, xmm8
    movapd xmm10, xmm11
    movapd xmm12, xmm11

	mulsd xmm3, [rsp + nb314_dxH1H1]
	mulsd xmm4, [rsp + nb314_dyH1H1]
	mulsd xmm5, [rsp + nb314_dzH1H1]
	mulsd xmm7, [rsp + nb314_dxH2H1]
	mulsd xmm8, [rsp + nb314_dyH2H1]
	mulsd xmm9, [rsp + nb314_dzH2H1]
	mulsd xmm10, [rsp + nb314_dxMH1]
	mulsd xmm11, [rsp + nb314_dyMH1]
	mulsd xmm12, [rsp + nb314_dzMH1]

    addsd xmm0, xmm3
    addsd xmm1, xmm4
    addsd xmm2, xmm5
    addsd xmm3, [rsp + nb314_fixH1]
    addsd xmm4, [rsp + nb314_fiyH1]
    addsd xmm5, [rsp + nb314_fizH1]

    addsd xmm0, xmm7
    addsd xmm1, xmm8
    addsd xmm2, xmm9
    addsd xmm7, [rsp + nb314_fixH2]
    addsd xmm8, [rsp + nb314_fiyH2]
    addsd xmm9, [rsp + nb314_fizH2]

    addsd xmm0, xmm10
    addsd xmm1, xmm11
    addsd xmm2, xmm12
    addsd xmm10, [rsp + nb314_fixM]
    addsd xmm11, [rsp + nb314_fiyM]
    addsd xmm12, [rsp + nb314_fizM]

    movsd [rsp + nb314_fixH1], xmm3
    movsd [rsp + nb314_fiyH1], xmm4
    movsd [rsp + nb314_fizH1], xmm5
    movsd [rsp + nb314_fixH2], xmm7
    movsd [rsp + nb314_fiyH2], xmm8
    movsd [rsp + nb314_fizH2], xmm9
    movsd [rsp + nb314_fixM], xmm10
    movsd [rsp + nb314_fiyM], xmm11
    movsd [rsp + nb314_fizM], xmm12
   
    ;# store back j H1 forces from xmm0-xmm2
	movsd [rdi + rax*8 + 24], xmm0
	movsd [rdi + rax*8 + 32], xmm1
	movsd [rdi + rax*8 + 40], xmm2

	;# move j H2 coordinates to local temp variables 
    mov rsi, [rbp + nb314_pos]
    movsd xmm0, [rsi + rax*8 + 48] 
    movsd xmm1, [rsi + rax*8 + 56] 
    movsd xmm2, [rsi + rax*8 + 64] 

    ;# xmm0 = H2x
    ;# xmm1 = H2y
    ;# xmm2 = H2z
        
    movsd xmm3, xmm0
    movsd xmm4, xmm1
    movsd xmm5, xmm2
    movsd xmm6, xmm0
    movsd xmm7, xmm1
    movsd xmm8, xmm2
    
    subsd xmm0, [rsp + nb314_ixH1]
    subsd xmm1, [rsp + nb314_iyH1]
    subsd xmm2, [rsp + nb314_izH1]
    subsd xmm3, [rsp + nb314_ixH2]
    subsd xmm4, [rsp + nb314_iyH2]
    subsd xmm5, [rsp + nb314_izH2]
    subsd xmm6, [rsp + nb314_ixM]
    subsd xmm7, [rsp + nb314_iyM]
    subsd xmm8, [rsp + nb314_izM]
    
	movsd [rsp + nb314_dxH1H2], xmm0
	movsd [rsp + nb314_dyH1H2], xmm1
	movsd [rsp + nb314_dzH1H2], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movsd [rsp + nb314_dxH2H2], xmm3
	movsd [rsp + nb314_dyH2H2], xmm4
	movsd [rsp + nb314_dzH2H2], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	movsd [rsp + nb314_dxMH2], xmm6
	movsd [rsp + nb314_dyMH2], xmm7
	movsd [rsp + nb314_dzMH2], xmm8
	mulsd  xmm6, xmm6
	mulsd  xmm7, xmm7
	mulsd  xmm8, xmm8
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
    addsd  xmm6, xmm7
    addsd  xmm6, xmm8

	;# start doing invsqrt for jH2 atoms
    cvtsd2ss xmm1, xmm0
    cvtsd2ss xmm4, xmm3
    cvtsd2ss xmm7, xmm6
	rsqrtss xmm1, xmm1
	rsqrtss xmm4, xmm4
    rsqrtss xmm7, xmm7
    cvtss2sd xmm1, xmm1
    cvtss2sd xmm4, xmm4
    cvtss2sd xmm7, xmm7
	
	movsd  xmm2, xmm1
	movsd  xmm5, xmm4
    movsd  xmm8, xmm7
    
	mulsd   xmm1, xmm1 ;# lu*lu
	mulsd   xmm4, xmm4 ;# lu*lu
    mulsd   xmm7, xmm7 ;# lu*lu
		
	movsd  xmm9, [rsp + nb314_three]
	movsd  xmm10, xmm9
    movsd  xmm11, xmm9

	mulsd   xmm1, xmm0 ;# rsq*lu*lu
	mulsd   xmm4, xmm3 ;# rsq*lu*lu 
    mulsd   xmm7, xmm6 ;# rsq*lu*lu
	
	subsd   xmm9, xmm1
	subsd   xmm10, xmm4
    subsd   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm2
	mulsd   xmm10, xmm5
    mulsd   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movsd  xmm15, [rsp + nb314_half]
	mulsd   xmm9, xmm15  ;# first iteration for rinvH1H2 
	mulsd   xmm10, xmm15 ;# first iteration for rinvH2H2
    mulsd   xmm11, xmm15 ;# first iteration for rinvMH2

    ;# second iteration step    
	movsd  xmm2, xmm9
	movsd  xmm5, xmm10
    movsd  xmm8, xmm11
    
	mulsd   xmm2, xmm2 ;# lu*lu
	mulsd   xmm5, xmm5 ;# lu*lu
    mulsd   xmm8, xmm8 ;# lu*lu
		
	movsd  xmm1, [rsp + nb314_three]
	movsd  xmm4, xmm1
    movsd  xmm7, xmm1

	mulsd   xmm2, xmm0 ;# rsq*lu*lu
	mulsd   xmm5, xmm3 ;# rsq*lu*lu 
    mulsd   xmm8, xmm6 ;# rsq*lu*lu
	
	subsd   xmm1, xmm2
	subsd   xmm4, xmm5
    subsd   xmm7, xmm8 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm1
	mulsd   xmm10, xmm4
    mulsd   xmm11, xmm7 ;# lu*(3-rsq*lu*lu)

	movsd  xmm15, [rsp + nb314_half]
	mulsd   xmm9, xmm15  ;#  rinvH1H2
	mulsd   xmm10, xmm15 ;#   rinvH2H2
    mulsd   xmm11, xmm15 ;#   rinvMH2
	
	movsd  [rsp + nb314_rinvH1H2], xmm9
	movsd  [rsp + nb314_rinvH2H2], xmm10
	movsd  [rsp + nb314_rinvMH2], xmm11
	
	;# H2 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movsd  xmm1, [rsp + nb314_tsc]
    mulsd  xmm0, xmm9  ;# r
    mulsd  xmm3, xmm10
    mulsd  xmm6, xmm11
    mulsd  xmm0, xmm1 ;# rtab
    mulsd  xmm3, xmm1
    mulsd  xmm6, xmm1
    
    ;# truncate and convert to integers
    cvttsd2si r8d, xmm0
    cvttsd2si r10d, xmm3
    cvttsd2si r12d, xmm6        

    ;# convert back to float
    cvtsi2sd  xmm2, r8d
    cvtsi2sd  xmm5, r10d
    cvtsi2sd  xmm8, r12d
    
    ;# multiply by 4
    shl   r8, 2
    shl   r10d, 2
    shl   r12d, 2
    
    mov  rsi, [rbp + nb314_VFtab]

    ;# calculate eps
    subsd     xmm0, xmm2
    subsd     xmm3, xmm5
    subsd     xmm6, xmm8

    movsd    [rsp + nb314_epsH1], xmm0
    movsd    [rsp + nb314_epsH2], xmm3
    movsd    [rsp + nb314_epsM], xmm6

    ;# Load LOTS of table data
    movsd xmm0,  [rsi + r8*8]
    movsd xmm1,  [rsi + r8*8 + 8]
    movsd xmm2,  [rsi + r8*8 + 16]
    movsd xmm3,  [rsi + r8*8 + 24]
    movsd xmm4,  [rsi + r10*8]
    movsd xmm5,  [rsi + r10*8 + 8]
    movsd xmm6,  [rsi + r10*8 + 16]
    movsd xmm7,  [rsi + r10*8 + 24]
    movsd xmm8,  [rsi + r12*8]
    movsd xmm9,  [rsi + r12*8 + 8]
    movsd xmm10, [rsi + r12*8 + 16]
    movsd xmm11, [rsi + r12*8 + 24]
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    movsd xmm12, [rsp + nb314_epsH1]
    movsd xmm13, [rsp + nb314_epsH2]
    movsd xmm14, [rsp + nb314_epsM]
    
    mulsd  xmm3, xmm12   ;# Heps
    mulsd  xmm7, xmm13
    mulsd  xmm11, xmm14 
    mulsd  xmm2, xmm12   ;# Geps
    mulsd  xmm6, xmm13
    mulsd  xmm10, xmm14 
    mulsd  xmm3, xmm12   ;# Heps2
    mulsd  xmm7, xmm13
    mulsd  xmm11, xmm14 

    addsd  xmm1, xmm2   ;# F+Geps
    addsd  xmm5, xmm6
    addsd  xmm9, xmm10 
    addsd  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addsd  xmm5, xmm7
    addsd  xmm9, xmm11 
    addsd  xmm3, xmm3    ;# 2*Heps2
    addsd  xmm7, xmm7
    addsd  xmm11, xmm11
    addsd  xmm3, xmm2    ;# 2*Heps2+Geps
    addsd  xmm7, xmm6  
    addsd  xmm11, xmm10
    addsd  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addsd  xmm7, xmm5
    addsd  xmm11, xmm9
    mulsd  xmm1, xmm12   ;# eps*Fp
    mulsd  xmm5, xmm13
    mulsd  xmm9, xmm14
    movsd xmm12, [rsp + nb314_qqHH]
    movsd xmm13, [rsp + nb314_qqMH]
    addsd  xmm1, xmm0     ;# VV
    addsd  xmm5, xmm4
    addsd  xmm9, xmm8
    mulsd  xmm1, xmm12   ;# VV*qq = vcoul
    mulsd  xmm5, xmm12
    mulsd  xmm9, xmm13
    mulsd  xmm3, xmm12    ;# FF*qq = fij
    mulsd  xmm7, xmm12
    mulsd  xmm11, xmm13
    
    ;# accumulate vctot
    addsd  xmm1, [rsp + nb314_vctot]
    addsd  xmm5, xmm9
    addsd  xmm1, xmm5
    movsd [rsp + nb314_vctot], xmm1
    
    movsd xmm10, [rsp + nb314_tsc]
    mulsd  xmm3, xmm10  ;# fscal
    mulsd  xmm7, xmm10
    mulsd  xmm10, xmm11
    
    xorpd xmm4, xmm4
    xorpd xmm8, xmm8
    xorpd xmm11, xmm11
    
    subsd xmm4, xmm3
    subsd xmm8, xmm7
    subsd xmm11, xmm10

    mulsd xmm4, [rsp + nb314_rinvH1H2]
    mulsd xmm8, [rsp + nb314_rinvH2H2]
    mulsd xmm11, [rsp + nb314_rinvMH2]
    
    ;# move j H2 forces to xmm0-xmm2
    mov rdi, [rbp + nb314_faction]
	movsd xmm0, [rdi + rax*8 + 48]
	movsd xmm1, [rdi + rax*8 + 56]
	movsd xmm2, [rdi + rax*8 + 64]

    movapd xmm3, xmm4
    movapd xmm5, xmm4
    movapd xmm7, xmm8
    movapd xmm9, xmm8
    movapd xmm10, xmm11
    movapd xmm12, xmm11

	mulsd xmm3, [rsp + nb314_dxH1H2]
	mulsd xmm4, [rsp + nb314_dyH1H2]
	mulsd xmm5, [rsp + nb314_dzH1H2]
	mulsd xmm7, [rsp + nb314_dxH2H2]
	mulsd xmm8, [rsp + nb314_dyH2H2]
	mulsd xmm9, [rsp + nb314_dzH2H2]
	mulsd xmm10, [rsp + nb314_dxMH2]
	mulsd xmm11, [rsp + nb314_dyMH2]
	mulsd xmm12, [rsp + nb314_dzMH2]

    addsd xmm0, xmm3
    addsd xmm1, xmm4
    addsd xmm2, xmm5
    addsd xmm3, [rsp + nb314_fixH1]
    addsd xmm4, [rsp + nb314_fiyH1]
    addsd xmm5, [rsp + nb314_fizH1]

    addsd xmm0, xmm7
    addsd xmm1, xmm8
    addsd xmm2, xmm9
    addsd xmm7, [rsp + nb314_fixH2]
    addsd xmm8, [rsp + nb314_fiyH2]
    addsd xmm9, [rsp + nb314_fizH2]

    addsd xmm0, xmm10
    addsd xmm1, xmm11
    addsd xmm2, xmm12
    addsd xmm10, [rsp + nb314_fixM]
    addsd xmm11, [rsp + nb314_fiyM]
    addsd xmm12, [rsp + nb314_fizM]

    movsd [rsp + nb314_fixH1], xmm3
    movsd [rsp + nb314_fiyH1], xmm4
    movsd [rsp + nb314_fizH1], xmm5
    movsd [rsp + nb314_fixH2], xmm7
    movsd [rsp + nb314_fiyH2], xmm8
    movsd [rsp + nb314_fizH2], xmm9
    movsd [rsp + nb314_fixM], xmm10
    movsd [rsp + nb314_fiyM], xmm11
    movsd [rsp + nb314_fizM], xmm12
   
    ;# store back j H2 forces from xmm0-xmm2
	movsd [rdi + rax*8 + 48], xmm0
	movsd [rdi + rax*8 + 56], xmm1
	movsd [rdi + rax*8 + 64], xmm2
       
	;# move j M coordinates to local temp variables 
    mov rsi, [rbp + nb314_pos]
    movsd xmm0, [rsi + rax*8 + 72] 
    movsd xmm1, [rsi + rax*8 + 80] 
    movsd xmm2, [rsi + rax*8 + 88] 

    ;# xmm0 = Mx
    ;# xmm1 = My
    ;# xmm2 = Mz
        
    movsd xmm3, xmm0
    movsd xmm4, xmm1
    movsd xmm5, xmm2
    movsd xmm6, xmm0
    movsd xmm7, xmm1
    movsd xmm8, xmm2
    
    subsd xmm0, [rsp + nb314_ixH1]
    subsd xmm1, [rsp + nb314_iyH1]
    subsd xmm2, [rsp + nb314_izH1]
    subsd xmm3, [rsp + nb314_ixH2]
    subsd xmm4, [rsp + nb314_iyH2]
    subsd xmm5, [rsp + nb314_izH2]
    subsd xmm6, [rsp + nb314_ixM]
    subsd xmm7, [rsp + nb314_iyM]
    subsd xmm8, [rsp + nb314_izM]
    
	movsd [rsp + nb314_dxH1M], xmm0
	movsd [rsp + nb314_dyH1M], xmm1
	movsd [rsp + nb314_dzH1M], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movsd [rsp + nb314_dxH2M], xmm3
	movsd [rsp + nb314_dyH2M], xmm4
	movsd [rsp + nb314_dzH2M], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	movsd [rsp + nb314_dxMM], xmm6
	movsd [rsp + nb314_dyMM], xmm7
	movsd [rsp + nb314_dzMM], xmm8
	mulsd  xmm6, xmm6
	mulsd  xmm7, xmm7
	mulsd  xmm8, xmm8
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
    addsd  xmm6, xmm7
    addsd  xmm6, xmm8

	;# start doing invsqrt for jM atoms
    cvtsd2ss xmm1, xmm0
    cvtsd2ss xmm4, xmm3
    cvtsd2ss xmm7, xmm6
	rsqrtss xmm1, xmm1
	rsqrtss xmm4, xmm4
    rsqrtss xmm7, xmm7
    cvtss2sd xmm1, xmm1
    cvtss2sd xmm4, xmm4
    cvtss2sd xmm7, xmm7
	
	movsd  xmm2, xmm1
	movsd  xmm5, xmm4
    movsd  xmm8, xmm7
    
	mulsd   xmm1, xmm1 ;# lu*lu
	mulsd   xmm4, xmm4 ;# lu*lu
    mulsd   xmm7, xmm7 ;# lu*lu
		
	movsd  xmm9, [rsp + nb314_three]
	movsd  xmm10, xmm9
    movsd  xmm11, xmm9

	mulsd   xmm1, xmm0 ;# rsq*lu*lu
	mulsd   xmm4, xmm3 ;# rsq*lu*lu 
    mulsd   xmm7, xmm6 ;# rsq*lu*lu
	
	subsd   xmm9, xmm1
	subsd   xmm10, xmm4
    subsd   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm2
	mulsd   xmm10, xmm5
    mulsd   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movsd  xmm15, [rsp + nb314_half]
	mulsd   xmm9, xmm15  ;# first iteration for rinvH1M 
	mulsd   xmm10, xmm15 ;# first iteration for rinvH2M
    mulsd   xmm11, xmm15 ;# first iteration for rinvMM

    ;# second iteration step    
	movsd  xmm2, xmm9
	movsd  xmm5, xmm10
    movsd  xmm8, xmm11
    
	mulsd   xmm2, xmm2 ;# lu*lu
	mulsd   xmm5, xmm5 ;# lu*lu
    mulsd   xmm8, xmm8 ;# lu*lu
		
	movsd  xmm1, [rsp + nb314_three]
	movsd  xmm4, xmm1
    movsd  xmm7, xmm1

	mulsd   xmm2, xmm0 ;# rsq*lu*lu
	mulsd   xmm5, xmm3 ;# rsq*lu*lu 
    mulsd   xmm8, xmm6 ;# rsq*lu*lu
	
	subsd   xmm1, xmm2
	subsd   xmm4, xmm5
    subsd   xmm7, xmm8 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm1
	mulsd   xmm10, xmm4
    mulsd   xmm11, xmm7 ;# lu*(3-rsq*lu*lu)

	movsd  xmm15, [rsp + nb314_half]
	mulsd   xmm9, xmm15  ;#  rinvH1M
	mulsd   xmm10, xmm15 ;#   rinvH2M
    mulsd   xmm11, xmm15 ;#   rinvMM
	
	movsd  [rsp + nb314_rinvH1M], xmm9
	movsd  [rsp + nb314_rinvH2M], xmm10
	movsd  [rsp + nb314_rinvMM], xmm11
	
	;# M interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movsd  xmm1, [rsp + nb314_tsc]
    mulsd  xmm0, xmm9  ;# r
    mulsd  xmm3, xmm10
    mulsd  xmm6, xmm11
    mulsd  xmm0, xmm1 ;# rtab
    mulsd  xmm3, xmm1
    mulsd  xmm6, xmm1
    
    ;# truncate and convert to integers
    cvttsd2si r8d, xmm0
    cvttsd2si r10d, xmm3
    cvttsd2si r12d, xmm6        

    ;# convert back to float
    cvtsi2sd  xmm2, r8d
    cvtsi2sd  xmm5, r10d
    cvtsi2sd  xmm8, r12d
    
    ;# multiply by 4
    shl   r8, 2
    shl   r10d, 2
    shl   r12d, 2
    
    mov  rsi, [rbp + nb314_VFtab]

    ;# calculate eps
    subsd     xmm0, xmm2
    subsd     xmm3, xmm5
    subsd     xmm6, xmm8

    movsd    [rsp + nb314_epsH1], xmm0
    movsd    [rsp + nb314_epsH2], xmm3
    movsd    [rsp + nb314_epsM], xmm6

    ;# Load LOTS of table data
    movsd xmm0,  [rsi + r8*8]
    movsd xmm1,  [rsi + r8*8 + 8]
    movsd xmm2,  [rsi + r8*8 + 16]
    movsd xmm3,  [rsi + r8*8 + 24]
    movsd xmm4,  [rsi + r10*8]
    movsd xmm5,  [rsi + r10*8 + 8]
    movsd xmm6,  [rsi + r10*8 + 16]
    movsd xmm7,  [rsi + r10*8 + 24]
    movsd xmm8,  [rsi + r12*8]
    movsd xmm9,  [rsi + r12*8 + 8]
    movsd xmm10, [rsi + r12*8 + 16]
    movsd xmm11, [rsi + r12*8 + 24]
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    movsd xmm12, [rsp + nb314_epsH1]
    movsd xmm13, [rsp + nb314_epsH2]
    movsd xmm14, [rsp + nb314_epsM]
    
    mulsd  xmm3, xmm12   ;# Heps
    mulsd  xmm7, xmm13
    mulsd  xmm11, xmm14 
    mulsd  xmm2, xmm12   ;# Geps
    mulsd  xmm6, xmm13
    mulsd  xmm10, xmm14 
    mulsd  xmm3, xmm12   ;# Heps2
    mulsd  xmm7, xmm13
    mulsd  xmm11, xmm14 

    addsd  xmm1, xmm2   ;# F+Geps
    addsd  xmm5, xmm6
    addsd  xmm9, xmm10 
    addsd  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addsd  xmm5, xmm7
    addsd  xmm9, xmm11 
    addsd  xmm3, xmm3    ;# 2*Heps2
    addsd  xmm7, xmm7
    addsd  xmm11, xmm11
    addsd  xmm3, xmm2    ;# 2*Heps2+Geps
    addsd  xmm7, xmm6  
    addsd  xmm11, xmm10
    addsd  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addsd  xmm7, xmm5
    addsd  xmm11, xmm9
    mulsd  xmm1, xmm12   ;# eps*Fp
    mulsd  xmm5, xmm13
    mulsd  xmm9, xmm14
    movsd xmm12, [rsp + nb314_qqMH]
    movsd xmm13, [rsp + nb314_qqMM]
    addsd  xmm1, xmm0     ;# VV
    addsd  xmm5, xmm4
    addsd  xmm9, xmm8
    mulsd  xmm1, xmm12   ;# VV*qq = vcoul
    mulsd  xmm5, xmm12
    mulsd  xmm9, xmm13
    mulsd  xmm3, xmm12    ;# FF*qq = fij
    mulsd  xmm7, xmm12
    mulsd  xmm11, xmm13
    
    ;# accumulate vctot
    addsd  xmm1, [rsp + nb314_vctot]
    addsd  xmm5, xmm9
    addsd  xmm1, xmm5
    movsd [rsp + nb314_vctot], xmm1
    
    movsd xmm10, [rsp + nb314_tsc]
    mulsd  xmm3, xmm10  ;# fscal
    mulsd  xmm7, xmm10
    mulsd  xmm10, xmm11
    
    xorpd xmm4, xmm4
    xorpd xmm8, xmm8
    xorpd xmm11, xmm11
    
    subsd xmm4, xmm3
    subsd xmm8, xmm7
    subsd xmm11, xmm10

    mulsd xmm4, [rsp + nb314_rinvH1M]
    mulsd xmm8, [rsp + nb314_rinvH2M]
    mulsd xmm11, [rsp + nb314_rinvMM]
    
    ;# move j M forces to xmm0-xmm2
    mov rdi, [rbp + nb314_faction]
	movsd xmm0, [rdi + rax*8 + 72]
	movsd xmm1, [rdi + rax*8 + 80]
	movsd xmm2, [rdi + rax*8 + 88]

    movapd xmm3, xmm4
    movapd xmm5, xmm4
    movapd xmm7, xmm8
    movapd xmm9, xmm8
    movapd xmm10, xmm11
    movapd xmm12, xmm11

	mulsd xmm3, [rsp + nb314_dxH1M]
	mulsd xmm4, [rsp + nb314_dyH1M]
	mulsd xmm5, [rsp + nb314_dzH1M]
	mulsd xmm7, [rsp + nb314_dxH2M]
	mulsd xmm8, [rsp + nb314_dyH2M]
	mulsd xmm9, [rsp + nb314_dzH2M]
	mulsd xmm10, [rsp + nb314_dxMM]
	mulsd xmm11, [rsp + nb314_dyMM]
	mulsd xmm12, [rsp + nb314_dzMM]

    addsd xmm0, xmm3
    addsd xmm1, xmm4
    addsd xmm2, xmm5
    addsd xmm3, [rsp + nb314_fixH1]
    addsd xmm4, [rsp + nb314_fiyH1]
    addsd xmm5, [rsp + nb314_fizH1]

    addsd xmm0, xmm7
    addsd xmm1, xmm8
    addsd xmm2, xmm9
    addsd xmm7, [rsp + nb314_fixH2]
    addsd xmm8, [rsp + nb314_fiyH2]
    addsd xmm9, [rsp + nb314_fizH2]

    addsd xmm0, xmm10
    addsd xmm1, xmm11
    addsd xmm2, xmm12
    addsd xmm10, [rsp + nb314_fixM]
    addsd xmm11, [rsp + nb314_fiyM]
    addsd xmm12, [rsp + nb314_fizM]

    movsd [rsp + nb314_fixH1], xmm3
    movsd [rsp + nb314_fiyH1], xmm4
    movsd [rsp + nb314_fizH1], xmm5
    movsd [rsp + nb314_fixH2], xmm7
    movsd [rsp + nb314_fiyH2], xmm8
    movsd [rsp + nb314_fizH2], xmm9
    movsd [rsp + nb314_fixM], xmm10
    movsd [rsp + nb314_fiyM], xmm11
    movsd [rsp + nb314_fizM], xmm12
   
    ;# store back j M forces from xmm0-xmm2
	movsd [rdi + rax*8 + 72], xmm0
	movsd [rdi + rax*8 + 80], xmm1
	movsd [rdi + rax*8 + 88], xmm2
	
.nb314_updateouterdata:
	mov   ecx, [rsp + nb314_ii3]
	mov   rdi, [rbp + nb314_faction]
	mov   rsi, [rbp + nb314_fshift]
	mov   edx, [rsp + nb314_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb314_fixO]
	movapd xmm1, [rsp + nb314_fiyO]
	movapd xmm2, [rsp + nb314_fizO]

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
	movsd  xmm3, [rdi + rcx*8]
	movsd  xmm4, [rdi + rcx*8 + 8]
	movsd  xmm5, [rdi + rcx*8 + 16]
	subsd  xmm3, xmm0
	subsd  xmm4, xmm1
	subsd  xmm5, xmm2
	movsd  [rdi + rcx*8],     xmm3
	movsd  [rdi + rcx*8 + 8], xmm4
	movsd  [rdi + rcx*8 + 16], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movapd xmm6, xmm0
	movsd xmm7, xmm2
	unpcklpd xmm6, xmm1

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb314_fixH1]
	movapd xmm1, [rsp + nb314_fiyH1]
	movapd xmm2, [rsp + nb314_fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	;# increment i force 
	movsd  xmm3, [rdi + rcx*8 + 24]
	movsd  xmm4, [rdi + rcx*8 + 32]
	movsd  xmm5, [rdi + rcx*8 + 40]
	subsd  xmm3, xmm0
	subsd  xmm4, xmm1
	subsd  xmm5, xmm2
	movsd  [rdi + rcx*8 + 24], xmm3
	movsd  [rdi + rcx*8 + 32], xmm4
	movsd  [rdi + rcx*8 + 40], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addsd xmm7, xmm2
	unpcklpd xmm0, xmm1
	addpd xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb314_fixH2]
	movapd xmm1, [rsp + nb314_fiyH2]
	movapd xmm2, [rsp + nb314_fizH2]

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
	movsd  xmm3, [rdi + rcx*8 + 48]
	movsd  xmm4, [rdi + rcx*8 + 56]
	movsd  xmm5, [rdi + rcx*8 + 64]
	subsd  xmm3, xmm0
	subsd  xmm4, xmm1
	subsd  xmm5, xmm2
	movsd  [rdi + rcx*8 + 48], xmm3
	movsd  [rdi + rcx*8 + 56], xmm4
	movsd  [rdi + rcx*8 + 64], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addsd xmm7, xmm2
	unpcklpd xmm0, xmm1
	addpd xmm6, xmm0

	;# accumulate Mi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb314_fixM]
	movapd xmm1, [rsp + nb314_fiyM]
	movapd xmm2, [rsp + nb314_fizM]

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
	movsd  xmm3, [rdi + rcx*8 + 72]
	movsd  xmm4, [rdi + rcx*8 + 80]
	movsd  xmm5, [rdi + rcx*8 + 88]
	subsd  xmm3, xmm0
	subsd  xmm4, xmm1
	subsd  xmm5, xmm2
	movsd  [rdi + rcx*8 + 72], xmm3
	movsd  [rdi + rcx*8 + 80], xmm4
	movsd  [rdi + rcx*8 + 88], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addsd xmm7, xmm2
	unpcklpd xmm0, xmm1
	addpd xmm6, xmm0

	;# increment fshift force 
	movlpd xmm3, [rsi + rdx*8]
	movhpd xmm3, [rsi + rdx*8 + 8]
	movsd  xmm4, [rsi + rdx*8 + 16]
	subpd  xmm3, xmm6
	subsd  xmm4, xmm7
	movlpd [rsi + rdx*8],      xmm3
	movhpd [rsi + rdx*8 + 8],  xmm3
	movsd  [rsi + rdx*8 + 16], xmm4

	;# get n from stack
	mov esi, [rsp + nb314_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb314_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb314_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb314_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb314_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb314_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb314_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb314_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb314_n], esi
        jmp .nb314_outer
.nb314_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb314_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb314_end
        ;# non-zero, do one more workunit
        jmp   .nb314_threadloop
.nb314_end:
	mov eax, [rsp + nb314_nouter]
	mov ebx, [rsp + nb314_ninner]
	mov rcx, [rbp + nb314_outeriter]
	mov rdx, [rbp + nb314_inneriter]
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




	
.globl nb_kernel314nf_x86_64_sse2
.globl _nb_kernel314nf_x86_64_sse2
nb_kernel314nf_x86_64_sse2:	
_nb_kernel314nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb314nf_fshift,         16
.equiv          nb314nf_gid,            24
.equiv          nb314nf_pos,            32
.equiv          nb314nf_faction,        40
.equiv          nb314nf_charge,         48
.equiv          nb314nf_p_facel,        56
.equiv          nb314nf_argkrf,         64
.equiv          nb314nf_argcrf,         72
.equiv          nb314nf_Vc,             80
.equiv          nb314nf_type,           88
.equiv          nb314nf_p_ntype,        96
.equiv          nb314nf_vdwparam,       104
.equiv          nb314nf_Vvdw,           112
.equiv          nb314nf_p_tabscale,     120
.equiv          nb314nf_VFtab,          128
.equiv          nb314nf_invsqrta,       136
.equiv          nb314nf_dvda,           144
.equiv          nb314nf_p_gbtabscale,   152
.equiv          nb314nf_GBtab,          160
.equiv          nb314nf_p_nthreads,     168
.equiv          nb314nf_count,          176
.equiv          nb314nf_mtx,            184
.equiv          nb314nf_outeriter,      192
.equiv          nb314nf_inneriter,      200
.equiv          nb314nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
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
.equiv          nb314nf_tsc,            432
.equiv          nb314nf_c6,             448
.equiv          nb314nf_c12,            464
.equiv          nb314nf_vctot,          480
.equiv          nb314nf_Vvdwtot,        496
.equiv          nb314nf_half,           512
.equiv          nb314nf_three,          528
.equiv          nb314nf_rsqOO,          544
.equiv          nb314nf_rsqH1H1,        560
.equiv          nb314nf_rsqH1H2,        576
.equiv          nb314nf_rsqH1M,         592
.equiv          nb314nf_rsqH2H1,        608
.equiv          nb314nf_rsqH2H2,        624
.equiv          nb314nf_rsqH2M,         640
.equiv          nb314nf_rsqMH1,         656
.equiv          nb314nf_rsqMH2,         672
.equiv          nb314nf_rsqMM,          688
.equiv          nb314nf_rinvsqOO,       704
.equiv          nb314nf_rinvH1H1,       720
.equiv          nb314nf_rinvH1H2,       736
.equiv          nb314nf_rinvH1M,        752
.equiv          nb314nf_rinvH2H1,       768
.equiv          nb314nf_rinvH2H2,       784
.equiv          nb314nf_rinvH2M,        800
.equiv          nb314nf_rinvMH1,        816
.equiv          nb314nf_rinvMH2,        832
.equiv          nb314nf_rinvMM,         848
.equiv          nb314nf_is3,            864
.equiv          nb314nf_ii3,            868
.equiv          nb314nf_nri,            872
.equiv          nb314nf_iinr,           880
.equiv          nb314nf_jindex,         888
.equiv          nb314nf_jjnr,           896
.equiv          nb314nf_shift,          904
.equiv          nb314nf_shiftvec,       912
.equiv          nb314nf_facel,          920
.equiv          nb314nf_innerjjnr,      928
.equiv          nb314nf_innerk,         936
.equiv          nb314nf_n,              940
.equiv          nb314nf_nn1,            944
.equiv          nb314nf_nouter,         948
.equiv          nb314nf_ninner,         952

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

	emms
	sub rsp, 960		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb314nf_nouter], eax
	mov [rsp + nb314nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb314nf_nri], edi
	mov [rsp + nb314nf_iinr], rsi
	mov [rsp + nb314nf_jindex], rdx
	mov [rsp + nb314nf_jjnr], rcx
	mov [rsp + nb314nf_shift], r8
	mov [rsp + nb314nf_shiftvec], r9
	mov rsi, [rbp + nb314nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb314nf_facel], xmm0

	mov rax, [rbp + nb314nf_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb314nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb314nf_half], eax
	mov [rsp + nb314nf_half+4], ebx
	movsd xmm1, [rsp + nb314nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb314nf_half], xmm1
	movapd [rsp + nb314nf_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb314nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb314nf_charge]
	movsd xmm3, [rdx + rbx*8 + 24]	
	movsd xmm4, xmm3	
	movsd xmm5, [rdx + rbx*8 + 8]	
	mov rsi, [rbp + nb314nf_p_facel]
	movsd xmm0, [rsi]
	movsd xmm6, [rsp + nb314nf_facel]
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm5
	mulsd  xmm5, xmm5
	mulsd  xmm3, xmm6
	mulsd  xmm4, xmm6
	mulsd  xmm5, xmm6
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb314nf_qqMM], xmm3
	movapd [rsp + nb314nf_qqMH], xmm4
	movapd [rsp + nb314nf_qqHH], xmm5
	
	xorpd xmm0, xmm0
	mov   rdx, [rbp + nb314nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb314nf_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   rax, [rbp + nb314nf_vdwparam]
	movlpd xmm0, [rax + rdx*8]
	movlpd xmm1, [rax + rdx*8 + 8]
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	movapd [rsp + nb314nf_c6], xmm0
	movapd [rsp + nb314nf_c12], xmm1

.nb314nf_threadloop:
        mov   rsi, [rbp + nb314nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb314nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb314nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb314nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb314nf_n], eax
        mov [rsp + nb314nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb314nf_outerstart
        jmp .nb314nf_end

.nb314nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb314nf_nouter]
	mov [rsp + nb314nf_nouter], ebx

.nb314nf_outer:
	mov   rax, [rsp + nb314nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb314nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb314nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb314nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb314nf_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb314nf_ii3], ebx		

	addsd xmm3, [rax + rbx*8] 	;# ox
	addsd xmm4, [rax + rbx*8 + 8] 	;# oy
	addsd xmm5, [rax + rbx*8 + 16]	;# oz	
	addsd xmm6, [rax + rbx*8 + 24] 	;# h1x
	addsd xmm7, [rax + rbx*8 + 32] 	;# h1y
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	shufpd xmm7, xmm7, 0
	movapd [rsp + nb314nf_ixO], xmm3
	movapd [rsp + nb314nf_iyO], xmm4
	movapd [rsp + nb314nf_izO], xmm5
	movapd [rsp + nb314nf_ixH1], xmm6
	movapd [rsp + nb314nf_iyH1], xmm7

	movsd xmm6, xmm2
	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm6, [rax + rbx*8 + 40] ;# h1z
	addsd xmm0, [rax + rbx*8 + 48] ;# h2x
	addsd xmm1, [rax + rbx*8 + 56] ;# h2y
	addsd xmm2, [rax + rbx*8 + 64] ;# h2z
	addsd xmm3, [rax + rbx*8 + 72] ;# mx
	addsd xmm4, [rax + rbx*8 + 80] ;# my
	addsd xmm5, [rax + rbx*8 + 88] ;# mz

	shufpd xmm6, xmm6, 0
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb314nf_izH1], xmm6
	movapd [rsp + nb314nf_ixH2], xmm0
	movapd [rsp + nb314nf_iyH2], xmm1
	movapd [rsp + nb314nf_izH2], xmm2
	movapd [rsp + nb314nf_ixM], xmm3
	movapd [rsp + nb314nf_iyM], xmm4
	movapd [rsp + nb314nf_izM], xmm5

	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [rsp + nb314nf_vctot], xmm4
	movapd [rsp + nb314nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb314nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb314nf_pos] 
	mov   rdi, [rbp + nb314nf_faction]	
	mov   rax, [rsp + nb314nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb314nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb314nf_ninner]
	mov   [rsp + nb314nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb314nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb314nf_unroll_loop
	jmp   .nb314nf_checksingle
.nb314nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb314nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	
	add qword ptr [rsp + nb314nf_innerjjnr], 8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb314nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	
	;# move j coordinates to local temp variables 
	;# load ox, oy, oz, h1x
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm2, [rsi + rbx*8]
	movhpd xmm0, [rsi + rax*8 + 8]
	movhpd xmm2, [rsi + rbx*8 + 8]
	movlpd xmm3, [rsi + rax*8 + 16]
	movlpd xmm5, [rsi + rbx*8 + 16]
	movhpd xmm3, [rsi + rax*8 + 24]
	movhpd xmm5, [rsi + rbx*8 + 24]
	movapd xmm1, xmm0 
	movapd xmm4, xmm3
	unpcklpd xmm0, xmm2 ;# ox 
	unpckhpd xmm1, xmm2 ;# oy
	unpcklpd xmm3, xmm5 ;# ox 
	unpckhpd xmm4, xmm5 ;# oy
	movapd 	[rsp + nb314nf_jxO], xmm0
	movapd 	[rsp + nb314nf_jyO], xmm1
	movapd 	[rsp + nb314nf_jzO], xmm3
	movapd 	[rsp + nb314nf_jxH1], xmm4
	
	;# load h1y, h1z, h2x, h2y 
	movlpd xmm0, [rsi + rax*8 + 32]
	movlpd xmm2, [rsi + rbx*8 + 32]
	movhpd xmm0, [rsi + rax*8 + 40]
	movhpd xmm2, [rsi + rbx*8 + 40]
	movlpd xmm3, [rsi + rax*8 + 48]
	movlpd xmm5, [rsi + rbx*8 + 48]
	movhpd xmm3, [rsi + rax*8 + 56]
	movhpd xmm5, [rsi + rbx*8 + 56]
	movapd xmm1, xmm0 
	movapd xmm4, xmm3
	unpcklpd xmm0, xmm2 ;# h1y
	unpckhpd xmm1, xmm2 ;# h1z
	unpcklpd xmm3, xmm5 ;# h2x
	unpckhpd xmm4, xmm5 ;# h2y
	movapd 	[rsp + nb314nf_jyH1], xmm0
	movapd 	[rsp + nb314nf_jzH1], xmm1
	movapd 	[rsp + nb314nf_jxH2], xmm3
	movapd 	[rsp + nb314nf_jyH2], xmm4
	
	;# load h2z, mx, my, mz
	movlpd xmm0, [rsi + rax*8 + 64]
	movlpd xmm2, [rsi + rbx*8 + 64]
	movhpd xmm0, [rsi + rax*8 + 72]
	movhpd xmm2, [rsi + rbx*8 + 72]
	movlpd xmm3, [rsi + rax*8 + 80]
	movlpd xmm5, [rsi + rbx*8 + 80]
	movhpd xmm3, [rsi + rax*8 + 88]
	movhpd xmm5, [rsi + rbx*8 + 88]
	movapd xmm1, xmm0 
	movapd xmm4, xmm3
	unpcklpd xmm0, xmm2 ;# h2z
	unpckhpd xmm1, xmm2 ;# mx
	unpcklpd xmm3, xmm5 ;# my
	unpckhpd xmm4, xmm5 ;# mz
	movapd 	[rsp + nb314nf_jzH2], xmm0
	movapd 	[rsp + nb314nf_jxM], xmm1
	movapd 	[rsp + nb314nf_jyM], xmm3
	movapd 	[rsp + nb314nf_jzM], xmm4
	
	;# start calculating pairwise distances
	movapd xmm0, [rsp + nb314nf_ixO]
	movapd xmm1, [rsp + nb314nf_iyO]
	movapd xmm2, [rsp + nb314nf_izO]
	movapd xmm3, [rsp + nb314nf_ixH1]
	movapd xmm4, [rsp + nb314nf_iyH1]
	movapd xmm5, [rsp + nb314nf_izH1]
	subpd  xmm0, [rsp + nb314nf_jxO]
	subpd  xmm1, [rsp + nb314nf_jyO]
	subpd  xmm2, [rsp + nb314nf_jzO]
	subpd  xmm3, [rsp + nb314nf_jxH1]
	subpd  xmm4, [rsp + nb314nf_jyH1]
	subpd  xmm5, [rsp + nb314nf_jzH1]
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
	movapd [rsp + nb314nf_rsqOO], xmm0
	movapd [rsp + nb314nf_rsqH1H1], xmm3

	movapd xmm0, [rsp + nb314nf_ixH1]
	movapd xmm1, [rsp + nb314nf_iyH1]
	movapd xmm2, [rsp + nb314nf_izH1]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [rsp + nb314nf_jxH2]
	subpd  xmm1, [rsp + nb314nf_jyH2]
	subpd  xmm2, [rsp + nb314nf_jzH2]
	subpd  xmm3, [rsp + nb314nf_jxM]
	subpd  xmm4, [rsp + nb314nf_jyM]
	subpd  xmm5, [rsp + nb314nf_jzM]
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
	movapd [rsp + nb314nf_rsqH1H2], xmm0
	movapd [rsp + nb314nf_rsqH1M], xmm3

	movapd xmm0, [rsp + nb314nf_ixH2]
	movapd xmm1, [rsp + nb314nf_iyH2]
	movapd xmm2, [rsp + nb314nf_izH2]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [rsp + nb314nf_jxH1]
	subpd  xmm1, [rsp + nb314nf_jyH1]
	subpd  xmm2, [rsp + nb314nf_jzH1]
	subpd  xmm3, [rsp + nb314nf_jxH2]
	subpd  xmm4, [rsp + nb314nf_jyH2]
	subpd  xmm5, [rsp + nb314nf_jzH2]
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
	movapd [rsp + nb314nf_rsqH2H1], xmm0
	movapd [rsp + nb314nf_rsqH2H2], xmm3

	movapd xmm0, [rsp + nb314nf_ixH2]
	movapd xmm1, [rsp + nb314nf_iyH2]
	movapd xmm2, [rsp + nb314nf_izH2]
	movapd xmm3, [rsp + nb314nf_ixM]
	movapd xmm4, [rsp + nb314nf_iyM]
	movapd xmm5, [rsp + nb314nf_izM]
	subpd  xmm0, [rsp + nb314nf_jxM]
	subpd  xmm1, [rsp + nb314nf_jyM]
	subpd  xmm2, [rsp + nb314nf_jzM]
	subpd  xmm3, [rsp + nb314nf_jxH1]
	subpd  xmm4, [rsp + nb314nf_jyH1]
	subpd  xmm5, [rsp + nb314nf_jzH1]
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
	movapd [rsp + nb314nf_rsqH2M], xmm0
	movapd [rsp + nb314nf_rsqMH1], xmm4

	movapd xmm0, [rsp + nb314nf_ixM]
	movapd xmm1, [rsp + nb314nf_iyM]
	movapd xmm2, [rsp + nb314nf_izM]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [rsp + nb314nf_jxH2]
	subpd  xmm1, [rsp + nb314nf_jyH2]
	subpd  xmm2, [rsp + nb314nf_jzH2]
	subpd  xmm3, [rsp + nb314nf_jxM]
	subpd  xmm4, [rsp + nb314nf_jyM]
	subpd  xmm5, [rsp + nb314nf_jzM]
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
	movapd [rsp + nb314nf_rsqMH2], xmm0
	movapd [rsp + nb314nf_rsqMM], xmm4
	
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
	movapd  xmm3, [rsp + nb314nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb314nf_half] ;# iter1 
	mulpd   xmm7, [rsp + nb314nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb314nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb314nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb314nf_half] ;# rinv 
	movapd [rsp + nb314nf_rinvMH2], xmm1
	movapd [rsp + nb314nf_rinvMM], xmm5

	movapd xmm0, [rsp + nb314nf_rsqOO]
	movapd xmm4, [rsp + nb314nf_rsqH1H1]	
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
	movapd  xmm3, [rsp + nb314nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb314nf_half] ;# iter1 of  
	mulpd   xmm7, [rsp + nb314nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb314nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb314nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb314nf_half] ;# rinv
	mulpd   xmm1, xmm1
	movapd [rsp + nb314nf_rinvsqOO], xmm1
	movapd [rsp + nb314nf_rinvH1H1], xmm5

	movapd xmm0, [rsp + nb314nf_rsqH1H2]
	movapd xmm4, [rsp + nb314nf_rsqH1M]	
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
	movapd  xmm3, [rsp + nb314nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb314nf_half] ;# iter1 
	mulpd   xmm7, [rsp + nb314nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb314nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb314nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb314nf_half] ;# rinv 
	movapd [rsp + nb314nf_rinvH1H2], xmm1
	movapd [rsp + nb314nf_rinvH1M], xmm5

	movapd xmm0, [rsp + nb314nf_rsqH2H1]
	movapd xmm4, [rsp + nb314nf_rsqH2H2]	
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
	movapd  xmm3, [rsp + nb314nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb314nf_half] ;# iter1a 
	mulpd   xmm7, [rsp + nb314nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb314nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb314nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb314nf_half] ;# rinv 
	movapd [rsp + nb314nf_rinvH2H1], xmm1
	movapd [rsp + nb314nf_rinvH2H2], xmm5

	movapd xmm0, [rsp + nb314nf_rsqMH1]
	movapd xmm4, [rsp + nb314nf_rsqH2M]	
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
	movapd  xmm3, [rsp + nb314nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb314nf_half] ;# iter1a 
	mulpd   xmm7, [rsp + nb314nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb314nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb314nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb314nf_half] ;# rinv 
	movapd [rsp + nb314nf_rinvMH1], xmm1
	movapd [rsp + nb314nf_rinvH2M], xmm5

	;# start with OO interaction 
	movapd xmm0, [rsp + nb314nf_rinvsqOO] ;# xmm0=rinvsq
	movapd  xmm1, xmm0	
	mulpd   xmm1, xmm1 ;# rinv4
	mulpd   xmm1, xmm0 ;#rinvsix
	movapd  xmm2, xmm1
	mulpd	xmm2, xmm2 ;# rinvtwelve
	mulpd  xmm1, [rsp + nb314nf_c6]
	mulpd  xmm2, [rsp + nb314nf_c12]
	movapd xmm3, xmm2
	subpd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addpd  xmm3, [rsp + nb314nf_Vvdwtot]
	movapd [rsp + nb314nf_Vvdwtot], xmm3

	;# H1-H1 interaction 
	movapd xmm0, [rsp + nb314nf_rinvH1H1]
	movapd xmm1, xmm0
	mulpd  xmm1, [rsp + nb314nf_rsqH1H1] ;# xmm1=r 
	mulpd  xmm1, [rsp + nb314nf_tsc]

	movd mm0, eax
	movd mm1, ebx

	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb314nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb314nf_qqHH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    	addpd  xmm5, [rsp + nb314nf_vctot]
    	movapd [rsp + nb314nf_vctot], xmm5

	;# H1-H2 interaction  
	movapd xmm0, [rsp + nb314nf_rinvH1H2]
	movapd xmm1, xmm0
	mulpd  xmm1, [rsp + nb314nf_rsqH1H2] ;# xmm1=r 
	mulpd  xmm1, [rsp + nb314nf_tsc]
	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb314nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb314nf_qqHH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    	addpd  xmm5, [rsp + nb314nf_vctot]
    	movapd [rsp + nb314nf_vctot], xmm5

	;# H1-M interaction 
	movapd xmm0, [rsp + nb314nf_rinvH1M]
	movapd xmm1, xmm0
	mulpd  xmm1, [rsp + nb314nf_rsqH1M] ;# xmm1=r 
	mulpd  xmm1, [rsp + nb314nf_tsc]
	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb314nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb314nf_qqMH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addpd  xmm5, [rsp + nb314nf_vctot]
    	movapd [rsp + nb314nf_vctot], xmm5

	;# H2-H1 interaction 
	movapd xmm0, [rsp + nb314nf_rinvH2H1]
	movapd xmm1, xmm0
	mulpd  xmm1, [rsp + nb314nf_rsqH2H1] ;# xmm1=r 
	mulpd  xmm1, [rsp + nb314nf_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb314nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb314nf_qqHH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addpd  xmm5, [rsp + nb314nf_vctot]
    	movapd [rsp + nb314nf_vctot], xmm5

	;# H2-H2 interaction 
	movapd xmm0, [rsp + nb314nf_rinvH2H2]
	movapd xmm1, xmm0
	mulpd  xmm1, [rsp + nb314nf_rsqH2H2] ;# xmm1=r 
	mulpd  xmm1, [rsp + nb314nf_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb314nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb314nf_qqHH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addpd  xmm5, [rsp + nb314nf_vctot]
    	movapd [rsp + nb314nf_vctot], xmm5

	;# H2-M interaction 
	movapd xmm0, [rsp + nb314nf_rinvH2M]
	movapd xmm1, xmm0
	mulpd  xmm1, [rsp + nb314nf_rsqH2M] ;# xmm1=r 
	mulpd  xmm1, [rsp + nb314nf_tsc]
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb314nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb314nf_qqMH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addpd  xmm5, [rsp + nb314nf_vctot]
    	movapd [rsp + nb314nf_vctot], xmm5

	;# M-H1 interaction 
	movapd xmm0, [rsp + nb314nf_rinvMH1]
	movapd xmm1, xmm0
	mulpd  xmm1, [rsp + nb314nf_rsqMH1] ;# xmm1=r 
	mulpd  xmm1, [rsp + nb314nf_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb314nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb314nf_qqMH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addpd  xmm5, [rsp + nb314nf_vctot]
    	movapd [rsp + nb314nf_vctot], xmm5

	;# M-H2 interaction 
	movapd xmm0, [rsp + nb314nf_rinvMH2]
	movapd xmm1, xmm0
	mulpd  xmm1, [rsp + nb314nf_rsqMH2] ;# xmm1=r 
	mulpd  xmm1, [rsp + nb314nf_tsc]
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb314nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb314nf_qqMH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addpd  xmm5, [rsp + nb314nf_vctot]
    	movapd [rsp + nb314nf_vctot], xmm5

	;# M-M interaction 
	movapd xmm0, [rsp + nb314nf_rinvMM]
	movapd xmm1, xmm0
	mulpd  xmm1, [rsp + nb314nf_rsqMM] ;# xmm1=r 
	mulpd  xmm1, [rsp + nb314nf_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb314nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb314nf_qqMM]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	
    	addpd  xmm5, [rsp + nb314nf_vctot]
    	movapd [rsp + nb314nf_vctot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb314nf_innerk],  2
	jl    .nb314nf_checksingle
	jmp   .nb314nf_unroll_loop
.nb314nf_checksingle:
	mov   edx, [rsp + nb314nf_innerk]
	and   edx, 1
	jnz   .nb314nf_dosingle
	jmp   .nb314nf_updateouterdata
.nb314nf_dosingle:
	mov   rdx, [rsp + nb314nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]

	mov rsi, [rbp + nb314nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	
	;# move j coordinates to local temp variables 
	;# load ox, oy, oz, h1x
	movlpd xmm0, [rsi + rax*8] 
	movhpd xmm0, [rsi + rax*8 + 8]
	movlpd xmm1, [rsi + rax*8 + 16]
	movhpd xmm1, [rsi + rax*8 + 24]
	movlpd xmm2, [rsi + rax*8 + 32]
	movhpd xmm2, [rsi + rax*8 + 40]
	movlpd xmm3, [rsi + rax*8 + 48]
	movhpd xmm3, [rsi + rax*8 + 56]
	movlpd xmm4, [rsi + rax*8 + 64]
	movhpd xmm4, [rsi + rax*8 + 72]
	movlpd xmm5, [rsi + rax*8 + 80]
	movhpd xmm5, [rsi + rax*8 + 88]
	movsd  [rsp + nb314nf_jxO], xmm0
	movsd  [rsp + nb314nf_jzO], xmm1
	movsd  [rsp + nb314nf_jyH1], xmm2
	movsd  [rsp + nb314nf_jxH2], xmm3
	movsd  [rsp + nb314nf_jzH2], xmm4
	movsd  [rsp + nb314nf_jyM], xmm5
	unpckhpd xmm0, xmm0
	unpckhpd xmm1, xmm1
	unpckhpd xmm2, xmm2
	unpckhpd xmm3, xmm3
	unpckhpd xmm4, xmm4
	unpckhpd xmm5, xmm5
	movsd  [rsp + nb314nf_jyO], xmm0
	movsd  [rsp + nb314nf_jxH1], xmm1
	movsd  [rsp + nb314nf_jzH1], xmm2
	movsd  [rsp + nb314nf_jyH2], xmm3
	movsd  [rsp + nb314nf_jxM], xmm4
	movsd  [rsp + nb314nf_jzM], xmm5

	;# start calculating pairwise distances
	movapd xmm0, [rsp + nb314nf_ixO]
	movapd xmm1, [rsp + nb314nf_iyO]
	movapd xmm2, [rsp + nb314nf_izO]
	movapd xmm3, [rsp + nb314nf_ixH1]
	movapd xmm4, [rsp + nb314nf_iyH1]
	movapd xmm5, [rsp + nb314nf_izH1]
	subsd  xmm0, [rsp + nb314nf_jxO]
	subsd  xmm1, [rsp + nb314nf_jyO]
	subsd  xmm2, [rsp + nb314nf_jzO]
	subsd  xmm3, [rsp + nb314nf_jxH1]
	subsd  xmm4, [rsp + nb314nf_jyH1]
	subsd  xmm5, [rsp + nb314nf_jzH1]
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
	movapd [rsp + nb314nf_rsqOO], xmm0
	movapd [rsp + nb314nf_rsqH1H1], xmm3

	movapd xmm0, [rsp + nb314nf_ixH1]
	movapd xmm1, [rsp + nb314nf_iyH1]
	movapd xmm2, [rsp + nb314nf_izH1]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [rsp + nb314nf_jxH2]
	subsd  xmm1, [rsp + nb314nf_jyH2]
	subsd  xmm2, [rsp + nb314nf_jzH2]
	subsd  xmm3, [rsp + nb314nf_jxM]
	subsd  xmm4, [rsp + nb314nf_jyM]
	subsd  xmm5, [rsp + nb314nf_jzM]
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
	movapd [rsp + nb314nf_rsqH1H2], xmm0
	movapd [rsp + nb314nf_rsqH1M], xmm3

	movapd xmm0, [rsp + nb314nf_ixH2]
	movapd xmm1, [rsp + nb314nf_iyH2]
	movapd xmm2, [rsp + nb314nf_izH2]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [rsp + nb314nf_jxH1]
	subsd  xmm1, [rsp + nb314nf_jyH1]
	subsd  xmm2, [rsp + nb314nf_jzH1]
	subsd  xmm3, [rsp + nb314nf_jxH2]
	subsd  xmm4, [rsp + nb314nf_jyH2]
	subsd  xmm5, [rsp + nb314nf_jzH2]
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
	movapd [rsp + nb314nf_rsqH2H1], xmm0
	movapd [rsp + nb314nf_rsqH2H2], xmm3

	movapd xmm0, [rsp + nb314nf_ixH2]
	movapd xmm1, [rsp + nb314nf_iyH2]
	movapd xmm2, [rsp + nb314nf_izH2]
	movapd xmm3, [rsp + nb314nf_ixM]
	movapd xmm4, [rsp + nb314nf_iyM]
	movapd xmm5, [rsp + nb314nf_izM]
	subsd  xmm0, [rsp + nb314nf_jxM]
	subsd  xmm1, [rsp + nb314nf_jyM]
	subsd  xmm2, [rsp + nb314nf_jzM]
	subsd  xmm3, [rsp + nb314nf_jxH1]
	subsd  xmm4, [rsp + nb314nf_jyH1]
	subsd  xmm5, [rsp + nb314nf_jzH1]
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
	movapd [rsp + nb314nf_rsqH2M], xmm0
	movapd [rsp + nb314nf_rsqMH1], xmm4

	movapd xmm0, [rsp + nb314nf_ixM]
	movapd xmm1, [rsp + nb314nf_iyM]
	movapd xmm2, [rsp + nb314nf_izM]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [rsp + nb314nf_jxH2]
	subsd  xmm1, [rsp + nb314nf_jyH2]
	subsd  xmm2, [rsp + nb314nf_jzH2]
	subsd  xmm3, [rsp + nb314nf_jxM]
	subsd  xmm4, [rsp + nb314nf_jyM]
	subsd  xmm5, [rsp + nb314nf_jzM]
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
	movapd [rsp + nb314nf_rsqMH2], xmm0
	movapd [rsp + nb314nf_rsqMM], xmm4

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
	movapd  xmm3, [rsp + nb314nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb314nf_half] ;# iter1 
	mulsd   xmm7, [rsp + nb314nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb314nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb314nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb314nf_half] ;# rinv 
	movapd [rsp + nb314nf_rinvMH2], xmm1
	movapd [rsp + nb314nf_rinvMM], xmm5

	movapd xmm0, [rsp + nb314nf_rsqOO]
	movapd xmm4, [rsp + nb314nf_rsqH1H1]	
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
	movapd  xmm3, [rsp + nb314nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb314nf_half] ;# iter1 of  
	mulsd   xmm7, [rsp + nb314nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb314nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb314nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb314nf_half] ;# rinv
	mulsd   xmm1, xmm1
	movapd [rsp + nb314nf_rinvsqOO], xmm1
	movapd [rsp + nb314nf_rinvH1H1], xmm5

	movapd xmm0, [rsp + nb314nf_rsqH1H2]
	movapd xmm4, [rsp + nb314nf_rsqH1M]	
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
	movapd  xmm3, [rsp + nb314nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb314nf_half] ;# iter1 
	mulsd   xmm7, [rsp + nb314nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb314nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb314nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb314nf_half] ;# rinv 
	movapd [rsp + nb314nf_rinvH1H2], xmm1
	movapd [rsp + nb314nf_rinvH1M], xmm5

	movapd xmm0, [rsp + nb314nf_rsqH2H1]
	movapd xmm4, [rsp + nb314nf_rsqH2H2]	
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
	movapd  xmm3, [rsp + nb314nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb314nf_half] ;# iter1a 
	mulsd   xmm7, [rsp + nb314nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb314nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb314nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb314nf_half] ;# rinv 
	movapd [rsp + nb314nf_rinvH2H1], xmm1
	movapd [rsp + nb314nf_rinvH2H2], xmm5

	movapd xmm0, [rsp + nb314nf_rsqMH1]
	movapd xmm4, [rsp + nb314nf_rsqH2M]	
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
	movapd  xmm3, [rsp + nb314nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb314nf_half] ;# iter1a 
	mulsd   xmm7, [rsp + nb314nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb314nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb314nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb314nf_half] ;# rinv 
	movapd [rsp + nb314nf_rinvMH1], xmm1
	movapd [rsp + nb314nf_rinvH2M], xmm5


	;# start with OO interaction 
	movsd xmm0, [rsp + nb314nf_rinvsqOO] ;# xmm0=rinvsq
	movsd  xmm1, xmm0	
	mulsd   xmm1, xmm1 ;# rinv4
	mulsd   xmm1, xmm0 ;#rinvsix
	movsd  xmm2, xmm1
	mulsd	xmm2, xmm2 ;# rinvtwelve
	mulsd  xmm1, [rsp + nb314nf_c6]
	mulsd  xmm2, [rsp + nb314nf_c12]
	movsd xmm3, xmm2
	subsd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addsd  xmm3, [rsp + nb314nf_Vvdwtot]
	movsd [rsp + nb314nf_Vvdwtot], xmm3

	;# H1-H1 interaction 
	movapd xmm0, [rsp + nb314nf_rinvH1H1]
	movapd xmm1, xmm0
	mulsd  xmm1, [rsp + nb314nf_rsqH1H1] ;# xmm1=r 
	mulsd  xmm1, [rsp + nb314nf_tsc]

	movd mm0, eax

	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, rax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb314nf_VFtab]

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm4
	unpckhpd xmm5, xmm5

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm6
	unpckhpd xmm7, xmm7
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb314nf_qqHH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  


    	addsd  xmm5, [rsp + nb314nf_vctot]
    	movsd [rsp + nb314nf_vctot], xmm5

	;# H1-H2 interaction  
	movapd xmm0, [rsp + nb314nf_rinvH1H2]
	movapd xmm1, xmm0
	mulsd  xmm1, [rsp + nb314nf_rsqH1H2] ;# xmm1=r 
	mulsd  xmm1, [rsp + nb314nf_tsc]
	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb314nf_VFtab]

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm4
	unpckhpd xmm5, xmm5

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm6
	unpckhpd xmm7, xmm7
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb314nf_qqHH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    	addsd  xmm5, [rsp + nb314nf_vctot]
    	movsd [rsp + nb314nf_vctot], xmm5

	;# H1-M interaction 
	movapd xmm0, [rsp + nb314nf_rinvH1M]
	movapd xmm1, xmm0
	mulsd  xmm1, [rsp + nb314nf_rsqH1M] ;# xmm1=r 
	mulsd  xmm1, [rsp + nb314nf_tsc]
	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb314nf_VFtab]
	
	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm4
	unpckhpd xmm5, xmm5

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm6
	unpckhpd xmm7, xmm7
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb314nf_qqMH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addsd  xmm5, [rsp + nb314nf_vctot]
    	movsd [rsp + nb314nf_vctot], xmm5

	;# H2-H1 interaction 
	movapd xmm0, [rsp + nb314nf_rinvH2H1]
	movapd xmm1, xmm0
	mulsd  xmm1, [rsp + nb314nf_rsqH2H1] ;# xmm1=r 
	mulsd  xmm1, [rsp + nb314nf_tsc]	

	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb314nf_VFtab]

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm4
	unpckhpd xmm5, xmm5

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm6
	unpckhpd xmm7, xmm7
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb314nf_qqHH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addsd  xmm5, [rsp + nb314nf_vctot]
    	movsd [rsp + nb314nf_vctot], xmm5

	;# H2-H2 interaction 
	movapd xmm0, [rsp + nb314nf_rinvH2H2]
	movapd xmm1, xmm0
	mulsd  xmm1, [rsp + nb314nf_rsqH2H2] ;# xmm1=r 
	mulsd  xmm1, [rsp + nb314nf_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb314nf_VFtab]

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm4
	unpckhpd xmm5, xmm5

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm6
	unpckhpd xmm7, xmm7
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb314nf_qqHH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addsd  xmm5, [rsp + nb314nf_vctot]
    	movsd [rsp + nb314nf_vctot], xmm5

	;# H2-M interaction 
	movapd xmm0, [rsp + nb314nf_rinvH2M]
	movapd xmm1, xmm0
	mulsd  xmm1, [rsp + nb314nf_rsqH2M] ;# xmm1=r 
	mulsd  xmm1, [rsp + nb314nf_tsc]
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb314nf_VFtab]

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm4
	unpckhpd xmm5, xmm5

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm6
	unpckhpd xmm7, xmm7
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb314nf_qqMH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addsd  xmm5, [rsp + nb314nf_vctot]
    	movsd [rsp + nb314nf_vctot], xmm5

	;# M-H1 interaction 
	movapd xmm0, [rsp + nb314nf_rinvMH1]
	movapd xmm1, xmm0
	mulsd  xmm1, [rsp + nb314nf_rsqMH1] ;# xmm1=r 
	mulsd  xmm1, [rsp + nb314nf_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb314nf_VFtab]

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm4
	unpckhpd xmm5, xmm5

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm6
	unpckhpd xmm7, xmm7
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb314nf_qqMH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addsd  xmm5, [rsp + nb314nf_vctot]
    	movsd [rsp + nb314nf_vctot], xmm5

	;# M-H2 interaction 
	movapd xmm0, [rsp + nb314nf_rinvMH2]
	movapd xmm1, xmm0
	mulsd  xmm1, [rsp + nb314nf_rsqMH2] ;# xmm1=r 
	mulsd  xmm1, [rsp + nb314nf_tsc]
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb314nf_VFtab]

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm4
	unpckhpd xmm5, xmm5

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm6
	unpckhpd xmm7, xmm7
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb314nf_qqMH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  

    	addsd  xmm5, [rsp + nb314nf_vctot]
    	movsd [rsp + nb314nf_vctot], xmm5

	;# M-M interaction 
	movapd xmm0, [rsp + nb314nf_rinvMM]
	movapd xmm1, xmm0
	mulsd  xmm1, [rsp + nb314nf_rsqMM] ;# xmm1=r 
	mulsd  xmm1, [rsp + nb314nf_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb314nf_VFtab]

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm4
	unpckhpd xmm5, xmm5

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm6
	unpckhpd xmm7, xmm7
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb314nf_qqMM]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	
    	addsd  xmm5, [rsp + nb314nf_vctot]
    	movsd [rsp + nb314nf_vctot], xmm5
	
.nb314nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb314nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb314nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb314nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb314nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb314nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb314nf_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb314nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb314nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb314nf_n], esi
        jmp .nb314nf_outer
.nb314nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb314nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb314nf_end
        ;# non-zero, do one more workunit
        jmp   .nb314nf_threadloop
.nb314nf_end:
	mov eax, [rsp + nb314nf_nouter]
	mov ebx, [rsp + nb314nf_ninner]
	mov rcx, [rbp + nb314nf_outeriter]
	mov rdx, [rbp + nb314nf_inneriter]
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
