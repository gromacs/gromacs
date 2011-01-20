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

	
.globl nb_kernel214_x86_64_sse
.globl _nb_kernel214_x86_64_sse
nb_kernel214_x86_64_sse:	
_nb_kernel214_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb214_fshift,           16
.equiv          nb214_gid,              24
.equiv          nb214_pos,              32
.equiv          nb214_faction,          40
.equiv          nb214_charge,           48
.equiv          nb214_p_facel,          56
.equiv          nb214_argkrf,           64
.equiv          nb214_argcrf,           72
.equiv          nb214_Vc,               80
.equiv          nb214_type,             88
.equiv          nb214_p_ntype,          96
.equiv          nb214_vdwparam,         104
.equiv          nb214_Vvdw,             112
.equiv          nb214_p_tabscale,       120
.equiv          nb214_VFtab,            128
.equiv          nb214_invsqrta,         136
.equiv          nb214_dvda,             144
.equiv          nb214_p_gbtabscale,     152
.equiv          nb214_GBtab,            160
.equiv          nb214_p_nthreads,       168
.equiv          nb214_count,            176
.equiv          nb214_mtx,              184
.equiv          nb214_outeriter,        192
.equiv          nb214_inneriter,        200
.equiv          nb214_work,             208
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
.equiv          nb214_nri,              1824
.equiv          nb214_iinr,             1832
.equiv          nb214_jindex,           1840
.equiv          nb214_jjnr,             1848
.equiv          nb214_shift,            1856
.equiv          nb214_shiftvec,         1864
.equiv          nb214_facel,            1872
.equiv          nb214_innerk,           1880
.equiv          nb214_n,                1884
.equiv          nb214_nn1,              1888
.equiv          nb214_nouter,           1892
.equiv          nb214_ninner,           1896

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

	sub rsp, 1904		;# local variable stack space (n*16+8)
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
	mov [rsp + nb214_nouter], eax
	mov [rsp + nb214_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb214_nri], edi
	mov [rsp + nb214_iinr], rsi
	mov [rsp + nb214_jindex], rdx
	mov [rsp + nb214_jjnr], rcx
	mov [rsp + nb214_shift], r8
	mov [rsp + nb214_shiftvec], r9
	mov rsi, [rbp + nb214_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb214_facel], xmm0


	mov rsi, [rbp + nb214_argkrf]
	mov rdi, [rbp + nb214_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb214_krf], xmm1
	movaps [rsp + nb214_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb214_half], eax
	movss xmm1, [rsp + nb214_half]
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
	movaps [rsp + nb214_half],  xmm1
	movaps [rsp + nb214_two],  xmm2
	movaps [rsp + nb214_three],  xmm3
	movaps [rsp + nb214_six],  xmm4
	movaps [rsp + nb214_twelve],  xmm5

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb214_iinr]   ;# rcx = pointer into iinr[]
	mov   ebx, [rcx]		;# ebx =ii 

	mov   rdx, [rbp + nb214_charge]
	movss xmm5, [rdx + rbx*4 + 4]	
	movss xmm3, [rdx + rbx*4 + 12]	
	movss xmm4, xmm3	
	mov rsi, [rbp + nb214_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb214_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb214_qqMM], xmm3
	movaps [rsp + nb214_qqMH], xmm4
	movaps [rsp + nb214_qqHH], xmm5
	
	xorps xmm0, xmm0
	mov   rdx, [rbp + nb214_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb214_p_ntype]
	imul  ecx, [rdi]  ;# rcx = ntia = 2*ntype*type[ii0] 
	add   rdx, rcx
	mov   rax, [rbp + nb214_vdwparam]
	movlps xmm0, [rax + rdx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0x55
	movaps [rsp + nb214_c6], xmm0
	movaps [rsp + nb214_c12], xmm1

.nb214_threadloop:
        mov   rsi, [rbp + nb214_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb214_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   rbx, 1                           ;# rbx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb214_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb214_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb214_n], eax
        mov [rsp + nb214_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb214_outerstart
        jmp .nb214_end
	
.nb214_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb214_nouter]
	mov [rsp + nb214_nouter], ebx

.nb214_outer:
	mov   rax, [rsp + nb214_shift]  	;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb214_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb214_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb214_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb214_pos]	;# rax = base of pos[]  
	mov   [rsp + nb214_ii3], ebx	
	
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
	movaps [rsp + nb214_ixO], xmm3
	movaps [rsp + nb214_iyO], xmm4
	movaps [rsp + nb214_izO], xmm5
	movaps [rsp + nb214_ixH1], xmm6
	movaps [rsp + nb214_iyH1], xmm7

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
	movaps [rsp + nb214_izH1], xmm6
	movaps [rsp + nb214_ixH2], xmm0
	movaps [rsp + nb214_iyH2], xmm1
	movaps [rsp + nb214_izH2], xmm2
	movaps [rsp + nb214_ixM], xmm3
	movaps [rsp + nb214_iyM], xmm4
	movaps [rsp + nb214_izM], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb214_vctot], xmm4
	movaps [rsp + nb214_Vvdwtot], xmm4
	movaps [rsp + nb214_fixO], xmm4
	movaps [rsp + nb214_fiyO], xmm4
	movaps [rsp + nb214_fizO], xmm4
	movaps [rsp + nb214_fixH1], xmm4
	movaps [rsp + nb214_fiyH1], xmm4
	movaps [rsp + nb214_fizH1], xmm4
	movaps [rsp + nb214_fixH2], xmm4
	movaps [rsp + nb214_fiyH2], xmm4
	movaps [rsp + nb214_fizH2], xmm4
	movaps [rsp + nb214_fixM], xmm4
	movaps [rsp + nb214_fiyM], xmm4
	movaps [rsp + nb214_fizM], xmm4
	
	mov   rax, [rsp + nb214_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb214_pos]
	mov   rdi, [rbp + nb214_faction]	
	mov   rax, [rsp + nb214_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb214_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb214_ninner]
	mov   [rsp + nb214_ninner], ecx
	add   edx, 0
	mov   [rsp + nb214_innerk], edx	;# number of innerloop atoms 
	jge   .nb214_unroll_loop
	jmp   .nb214_single_check
.nb214_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb214_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb214_innerjjnr],  16 ;# advance pointer (unroll 4) 

	mov rsi, [rbp + nb214_pos]   	;# base of pos[] 

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
        
    subps xmm0, [rsp + nb214_ixO]
    subps xmm1, [rsp + nb214_iyO]
    subps xmm2, [rsp + nb214_izO]

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
    movss  xmm0,  [rdi + rbx*4 + 8] ;# jzb
    movss  xmm8,  [rdi + rdx*4 + 8] ;# jzd
    movlhps xmm12, xmm0 ;# jza  -  jzb  -
    movlhps xmm13, xmm8 ;# jzc  -  jzd -

    shufps xmm12, xmm13,  136  ;# 10001000 => jzOa jzOb jzOc jzOd

    ;# xmm10: jxOa jyOa jxOb jyOb 
    ;# xmm11: jxOc jyOc jxOd jyOd
    ;# xmm12: jzOa jzOb jzOc jzOd

    ;# calc rinvsq=1/rsq
	rcpps xmm2, xmm1
	movaps xmm0, [rsp + nb214_two]
	mulps xmm1, xmm2
	subps xmm0, xmm1
	mulps xmm0, xmm2	;# xmm0=rinvsq

	movaps xmm1, xmm0   ;# rinvsq

	mulps  xmm0, xmm0   ;# rinv4
	mulps  xmm0, xmm1	;# rinv6
	movaps xmm2, xmm0
	mulps  xmm2, xmm2	;# xmm2=rinv12

	mulps  xmm0, [rsp + nb214_c6]
	mulps  xmm2, [rsp + nb214_c12]
	movaps xmm8, xmm2
    subps  xmm2, xmm0	;# Vvdw=Vvdw12-Vvdw6 
	mulps  xmm0, [rsp + nb214_six]
	mulps  xmm8, [rsp + nb214_twelve]
	subps  xmm8, xmm0
	mulps  xmm8, xmm1	;# xmm8=total fscal 
        
    ;# add potential to Vvdwtot
	addps  xmm2, [rsp + nb214_Vvdwtot]
    movaps [rsp + nb214_Vvdwtot], xmm2

    ;# calculate scalar force by multiplying dx/dy/dz with fscal
	mulps  xmm4, xmm8
	mulps  xmm5, xmm8
	mulps  xmm6, xmm8

    ;# increment i force
    movaps xmm0, [rsp + nb214_fixO]
    movaps xmm1, [rsp + nb214_fiyO]
    movaps xmm2, [rsp + nb214_fizO]
    addps  xmm0, xmm4
    addps  xmm1, xmm5
    addps  xmm2, xmm6
    movaps [rsp + nb214_fixO], xmm0
    movaps [rsp + nb214_fiyO], xmm1
    movaps [rsp + nb214_fizO], xmm2

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
    
    subps xmm0, [rsp + nb214_ixH1]
    subps xmm1, [rsp + nb214_iyH1]
    subps xmm2, [rsp + nb214_izH1]
    subps xmm3, [rsp + nb214_ixH2]
    subps xmm4, [rsp + nb214_iyH2]
    subps xmm5, [rsp + nb214_izH2]
    subps xmm6, [rsp + nb214_ixM]
    subps xmm7, [rsp + nb214_iyM]
    subps xmm8, [rsp + nb214_izM]
    
	movaps [rsp + nb214_dxH1H1], xmm0
	movaps [rsp + nb214_dyH1H1], xmm1
	movaps [rsp + nb214_dzH1H1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb214_dxH2H1], xmm3
	movaps [rsp + nb214_dyH2H1], xmm4
	movaps [rsp + nb214_dzH2H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb214_dxMH1], xmm6
	movaps [rsp + nb214_dyMH1], xmm7
	movaps [rsp + nb214_dzMH1], xmm8
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
		
	movaps  xmm9, [rsp + nb214_three]
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

	movaps  xmm4, [rsp + nb214_half]
	mulps   xmm9, xmm4  ;# rinvH1H1 
	mulps   xmm10, xmm4 ;# rinvH2H1
    mulps   xmm11, xmm4 ;# rinvMH1
	
	;# H1 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movaps xmm1, xmm9 ;# copy of rinv
    movaps xmm4, xmm10
    movaps xmm7, xmm11
    movaps xmm2, [rsp + nb214_krf]    
    mulps  xmm9, xmm9   ;# rinvsq
    mulps  xmm10, xmm10
    mulps  xmm11, xmm11
    mulps  xmm0, xmm2  ;# k*rsq
    mulps  xmm3, xmm2
    mulps  xmm6, xmm2
    movaps xmm2, xmm0 ;# copy of k*rsq
    movaps xmm5, xmm3
    movaps xmm8, xmm6
    addps  xmm2, xmm1  ;# rinv+krsq
    addps  xmm5, xmm4
    addps  xmm8, xmm7
    movaps xmm14, [rsp + nb214_crf]
    subps  xmm2, xmm14   ;# rinv+krsq-crf
    subps  xmm5, xmm14
    subps  xmm8, xmm14
    movaps xmm12, [rsp + nb214_qqHH]
    movaps xmm13, [rsp + nb214_qqMH]    
    mulps  xmm2, xmm12 ;# voul=qq*(rinv+ krsq-crf)
    mulps  xmm5, xmm12 ;# voul=qq*(rinv+ krsq-crf)
    mulps  xmm8, xmm13 ;# voul=qq*(rinv+ krsq-crf)
    addps  xmm0, xmm0 ;# 2*krsq
    addps  xmm3, xmm3 
    addps  xmm6, xmm6 
    subps  xmm1, xmm0 ;# rinv-2*krsq
    subps  xmm4, xmm3
    subps  xmm7, xmm6
    mulps  xmm1, xmm12   ;# (rinv-2*krsq)*qq
    mulps  xmm4, xmm12
    mulps  xmm7, xmm13
    addps  xmm2, [rsp + nb214_vctot]
    addps  xmm5, xmm8
    addps  xmm2, xmm5
    movaps xmm15, xmm2
    
    mulps  xmm1, xmm9   ;# fscal
    mulps  xmm4, xmm10
    mulps  xmm7, xmm11

	;# move j H1 forces to local temp variables 
    movlps xmm9, [rdi + rax*4 + 12] ;# jxH1a jyH1a  -   -
    movlps xmm10, [rdi + rcx*4 + 12] ;# jxH1c jyH1c  -   -
    movhps xmm9, [rdi + rbx*4 + 12] ;# jxH1a jyH1a jxH1b jyH1b 
    movhps xmm10, [rdi + rdx*4 + 12] ;# jxH1c jyH1c jxH1d jyH1d 

    movss  xmm11, [rdi + rax*4 + 20] ;# jzH1a  -  -  -
    movss  xmm12, [rdi + rcx*4 + 20] ;# jzH1c  -  -  -
    movhps xmm11, [rdi + rbx*4 + 20] ;# jzH1a  -  jzH1b  -
    movhps xmm12, [rdi + rdx*4 + 20] ;# jzH1c  -  jzH1d -
    
    shufps xmm11, xmm12,  136  ;# 10001000 => jzH1a jzH1b jzH1c jzH1d

    ;# xmm9: jxH1a jyH1a jxH1b jyH1b 
    ;# xmm10: jxH1c jyH1c jxH1d jyH1d
    ;# xmm11: jzH1a jzH1b jzH1c jzH1d

    movaps xmm0, xmm1
    movaps xmm2, xmm1
    movaps xmm3, xmm4
    movaps xmm5, xmm4
    movaps xmm6, xmm7
    movaps xmm8, xmm7

	mulps xmm0, [rsp + nb214_dxH1H1]
	mulps xmm1, [rsp + nb214_dyH1H1]
	mulps xmm2, [rsp + nb214_dzH1H1]
	mulps xmm3, [rsp + nb214_dxH2H1]
	mulps xmm4, [rsp + nb214_dyH2H1]
	mulps xmm5, [rsp + nb214_dzH2H1]
	mulps xmm6, [rsp + nb214_dxMH1]
	mulps xmm7, [rsp + nb214_dyMH1]
	mulps xmm8, [rsp + nb214_dzMH1]

    movaps xmm13, xmm0
    movaps xmm14, xmm1
    addps xmm11, xmm2
    addps xmm0, [rsp + nb214_fixH1]
    addps xmm1, [rsp + nb214_fiyH1]
    addps xmm2, [rsp + nb214_fizH1]

    addps xmm13,  xmm3
    addps xmm14, xmm4
    addps xmm11, xmm5
    addps xmm3, [rsp + nb214_fixH2]
    addps xmm4, [rsp + nb214_fiyH2]
    addps xmm5, [rsp + nb214_fizH2]

    addps xmm13, xmm6
    addps xmm14, xmm7
    addps xmm11, xmm8
    addps xmm6, [rsp + nb214_fixM]
    addps xmm7, [rsp + nb214_fiyM]
    addps xmm8, [rsp + nb214_fizM]

    movaps [rsp + nb214_fixH1], xmm0
    movaps [rsp + nb214_fiyH1], xmm1
    movaps [rsp + nb214_fizH1], xmm2
    movaps [rsp + nb214_fixH2], xmm3
    movaps [rsp + nb214_fiyH2], xmm4
    movaps [rsp + nb214_fizH2], xmm5
    movaps [rsp + nb214_fixM], xmm6
    movaps [rsp + nb214_fiyM], xmm7
    movaps [rsp + nb214_fizM], xmm8
    
    ;# xmm9 = fH1x
    ;# xmm10 = fH1y
    ;# xmm11 = fH1z
    movaps xmm0, xmm13
    unpcklps xmm13, xmm14
    unpckhps xmm0, xmm14
    
    addps xmm9, xmm13
    addps xmm10, xmm0

    movhlps  xmm12, xmm11 ;# fH1zc fH1zd
    
    movlps [rdi + rax*4 + 12], xmm9
    movhps [rdi + rbx*4 + 12], xmm9
    movlps [rdi + rcx*4 + 12], xmm10
    movhps [rdi + rdx*4 + 12], xmm10
    movss  [rdi + rax*4 + 20], xmm11
    movss  [rdi + rcx*4 + 20], xmm12
    shufps xmm11, xmm11, 1
    shufps xmm12, xmm12, 1
    movss  [rdi + rbx*4 + 20], xmm11
    movss  [rdi + rdx*4 + 20], xmm12

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
    
    
    subps xmm0, [rsp + nb214_ixH1]
    subps xmm1, [rsp + nb214_iyH1]
    subps xmm2, [rsp + nb214_izH1]
    subps xmm3, [rsp + nb214_ixH2]
    subps xmm4, [rsp + nb214_iyH2]
    subps xmm5, [rsp + nb214_izH2]
    subps xmm6, [rsp + nb214_ixM]
    subps xmm7, [rsp + nb214_iyM]
    subps xmm8, [rsp + nb214_izM]
    
	movaps [rsp + nb214_dxH1H2], xmm0
	movaps [rsp + nb214_dyH1H2], xmm1
	movaps [rsp + nb214_dzH1H2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb214_dxH2H2], xmm3
	movaps [rsp + nb214_dyH2H2], xmm4
	movaps [rsp + nb214_dzH2H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb214_dxMH2], xmm6
	movaps [rsp + nb214_dyMH2], xmm7
	movaps [rsp + nb214_dzMH2], xmm8
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
		
	movaps  xmm9, [rsp + nb214_three]
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

	movaps  xmm4, [rsp + nb214_half]
	mulps   xmm9, xmm4  ;# rinvH1H2
	mulps   xmm10, xmm4 ;# rinvH2H2
    mulps   xmm11, xmm4 ;# rinvMH2
	
	;# H2 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movaps xmm1, xmm9 ;# copy of rinv
    movaps xmm4, xmm10
    movaps xmm7, xmm11
    movaps xmm2, [rsp + nb214_krf]    
    mulps  xmm9, xmm9   ;# rinvsq
    mulps  xmm10, xmm10
    mulps  xmm11, xmm11
    mulps  xmm0, xmm2  ;# k*rsq
    mulps  xmm3, xmm2
    mulps  xmm6, xmm2
    movaps xmm2, xmm0 ;# copy of k*rsq
    movaps xmm5, xmm3
    movaps xmm8, xmm6
    addps  xmm2, xmm1  ;# rinv+krsq
    addps  xmm5, xmm4
    addps  xmm8, xmm7
    movaps xmm14, [rsp + nb214_crf]
    subps  xmm2, xmm14   ;# rinv+krsq-crf
    subps  xmm5, xmm14
    subps  xmm8, xmm14
    movaps xmm12, [rsp + nb214_qqHH]
    movaps xmm13, [rsp + nb214_qqMH]    
    mulps  xmm2, xmm12 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    mulps  xmm5, xmm12 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    mulps  xmm8, xmm13 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    addps  xmm0, xmm0 ;# 2*krsq
    addps  xmm3, xmm3 
    addps  xmm6, xmm6 
    subps  xmm1, xmm0 ;# rinv-2*krsq
    subps  xmm4, xmm3
    subps  xmm7, xmm6
    mulps  xmm1, xmm12   ;# (rinv-2*krsq)*qq
    mulps  xmm4, xmm12
    mulps  xmm7, xmm13
    addps  xmm15, xmm2
    addps  xmm5, xmm8
    addps  xmm15, xmm5
    
    mulps  xmm1, xmm9   ;# fscal
    mulps  xmm4, xmm10
    mulps  xmm7, xmm11

	;# move j H2 forces to local temp variables 
    movlps xmm9, [rdi + rax*4 + 24] ;# jxH2a jyH2a  -   -
    movlps xmm10, [rdi + rcx*4 + 24] ;# jxH2c jyH2c  -   -
    movhps xmm9, [rdi + rbx*4 + 24] ;# jxH2a jyH2a jxH2b jyH2b 
    movhps xmm10, [rdi + rdx*4 + 24] ;# jxH2c jyH2c jxH2d jyH2d 

    movss  xmm11, [rdi + rax*4 + 32] ;# jzH2a  -  -  -
    movss  xmm12, [rdi + rcx*4 + 32] ;# jzH2c  -  -  -
    movhps xmm11, [rdi + rbx*4 + 32] ;# jzH2a  -  jzH2b  -
    movhps xmm12, [rdi + rdx*4 + 32] ;# jzH2c  -  jzH2d -
    
    shufps xmm11, xmm12,  136  ;# 10001000 => jzH2a jzH2b jzH2c jzH2d

    ;# xmm9: jxH2a jyH2a jxH2b jyH2b 
    ;# xmm10: jxH2c jyH2c jxH2d jyH2d
    ;# xmm11: jzH2a jzH2b jzH2c jzH2d

    movaps xmm0, xmm1
    movaps xmm2, xmm1
    movaps xmm3, xmm4
    movaps xmm5, xmm4
    movaps xmm6, xmm7
    movaps xmm8, xmm7

	mulps xmm0, [rsp + nb214_dxH1H2]
	mulps xmm1, [rsp + nb214_dyH1H2]
	mulps xmm2, [rsp + nb214_dzH1H2]
	mulps xmm3, [rsp + nb214_dxH2H2]
	mulps xmm4, [rsp + nb214_dyH2H2]
	mulps xmm5, [rsp + nb214_dzH2H2]
	mulps xmm6, [rsp + nb214_dxMH2]
	mulps xmm7, [rsp + nb214_dyMH2]
	mulps xmm8, [rsp + nb214_dzMH2]

    movaps xmm13, xmm0
    movaps xmm14, xmm1
    addps xmm11, xmm2
    addps xmm0, [rsp + nb214_fixH1]
    addps xmm1, [rsp + nb214_fiyH1]
    addps xmm2, [rsp + nb214_fizH1]

    addps xmm13, xmm3
    addps xmm14, xmm4
    addps xmm11, xmm5
    addps xmm3, [rsp + nb214_fixH2]
    addps xmm4, [rsp + nb214_fiyH2]
    addps xmm5, [rsp + nb214_fizH2]

    addps xmm13, xmm6
    addps xmm14, xmm7
    addps xmm11, xmm8
    addps xmm6, [rsp + nb214_fixM]
    addps xmm7, [rsp + nb214_fiyM]
    addps xmm8, [rsp + nb214_fizM]

    movaps [rsp + nb214_fixH1], xmm0
    movaps [rsp + nb214_fiyH1], xmm1
    movaps [rsp + nb214_fizH1], xmm2
    movaps [rsp + nb214_fixH2], xmm3
    movaps [rsp + nb214_fiyH2], xmm4
    movaps [rsp + nb214_fizH2], xmm5
    movaps [rsp + nb214_fixM], xmm6
    movaps [rsp + nb214_fiyM], xmm7
    movaps [rsp + nb214_fizM], xmm8
    
    ;# xmm9  = fH2x
    ;# xmm10 = fH2y
    ;# xmm11 = fH2z
    movaps xmm0, xmm13
    unpcklps xmm13, xmm14
    unpckhps xmm0, xmm14
    
    addps xmm9, xmm13
    addps xmm10, xmm0

    movhlps  xmm12, xmm11 ;# fH2zc fH2zd
    
    movlps [rdi + rax*4 + 24], xmm9
    movhps [rdi + rbx*4 + 24], xmm9
    movlps [rdi + rcx*4 + 24], xmm10
    movhps [rdi + rdx*4 + 24], xmm10
    movss  [rdi + rax*4 + 32], xmm11
    movss  [rdi + rcx*4 + 32], xmm12
    shufps xmm11, xmm11, 1
    shufps xmm12, xmm12, 1
    movss  [rdi + rbx*4 + 32], xmm11
    movss  [rdi + rdx*4 + 32], xmm12

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

    shufps xmm2, xmm3,  136  ;# 10001000 => jzMa jzMb jzMc jzMd

    ;# xmm0 = Mx
    ;# xmm1 = My
    ;# xmm2 = Mz
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + nb214_ixH1]
    subps xmm1, [rsp + nb214_iyH1]
    subps xmm2, [rsp + nb214_izH1]
    subps xmm3, [rsp + nb214_ixH2]
    subps xmm4, [rsp + nb214_iyH2]
    subps xmm5, [rsp + nb214_izH2]
    subps xmm6, [rsp + nb214_ixM]
    subps xmm7, [rsp + nb214_iyM]
    subps xmm8, [rsp + nb214_izM]
    
	movaps [rsp + nb214_dxH1M], xmm0
	movaps [rsp + nb214_dyH1M], xmm1
	movaps [rsp + nb214_dzH1M], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb214_dxH2M], xmm3
	movaps [rsp + nb214_dyH2M], xmm4
	movaps [rsp + nb214_dzH2M], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb214_dxMM], xmm6
	movaps [rsp + nb214_dyMM], xmm7
	movaps [rsp + nb214_dzMM], xmm8
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
		
	movaps  xmm9, [rsp + nb214_three]
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

	movaps  xmm4, [rsp + nb214_half]
	mulps   xmm9, xmm4  ;# rinvH1M
	mulps   xmm10, xmm4 ;# rinvH2M
    mulps   xmm11, xmm4 ;# rinvMM
	
	;# M interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movaps xmm1, xmm9 ;# copy of rinv
    movaps xmm4, xmm10
    movaps xmm7, xmm11
    movaps xmm2, [rsp + nb214_krf]    
    mulps  xmm9, xmm9   ;# rinvsq
    mulps  xmm10, xmm10
    mulps  xmm11, xmm11
    mulps  xmm0, xmm2  ;# k*rsq
    mulps  xmm3, xmm2
    mulps  xmm6, xmm2
    movaps xmm2, xmm0 ;# copy of k*rsq
    movaps xmm5, xmm3
    movaps xmm8, xmm6
    addps  xmm2, xmm1  ;# rinv+krsq
    addps  xmm5, xmm4
    addps  xmm8, xmm7
    movaps xmm14, [rsp + nb214_crf]
    subps  xmm2, xmm14   ;# rinv+krsq-crf
    subps  xmm5, xmm14
    subps  xmm8, xmm14
    movaps xmm12, [rsp + nb214_qqMH]
    movaps xmm13, [rsp + nb214_qqMM]    
    mulps  xmm2, xmm12 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    mulps  xmm5, xmm12 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    mulps  xmm8, xmm13 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    addps  xmm0, xmm0 ;# 2*krsq
    addps  xmm3, xmm3 
    addps  xmm6, xmm6 
    subps  xmm1, xmm0 ;# rinv-2*krsq
    subps  xmm4, xmm3
    subps  xmm7, xmm6
    mulps  xmm1, xmm12   ;# (rinv-2*krsq)*qq
    mulps  xmm4, xmm12
    mulps  xmm7, xmm13
    addps  xmm5, xmm8
    addps  xmm2, xmm15
    addps  xmm2, xmm5
    movaps [rsp + nb214_vctot], xmm2
    
    mulps  xmm1, xmm9   ;# fscal
    mulps  xmm4, xmm10
    mulps  xmm7, xmm11

	;# move j M forces to local temp variables 
    movlps xmm9, [rdi + rax*4 + 36] ;# jxMa jyMa  -   -
    movlps xmm10, [rdi + rcx*4 + 36] ;# jxMc jyMc  -   -
    movhps xmm9, [rdi + rbx*4 + 36] ;# jxMa jyMa jxMb jyMb 
    movhps xmm10, [rdi + rdx*4 + 36] ;# jxMc jyMc jxMd jyMd 

    movss  xmm11, [rdi + rax*4 + 44] ;# jzMa  -  -  -
    movss  xmm12, [rdi + rcx*4 + 44] ;# jzMc  -  -  -
    movss  xmm2,  [rdi + rbx*4 + 44] ;# jzMb  -  -  -
    movss  xmm3,  [rdi + rdx*4 + 44] ;# jzMd  -  -  -
    movlhps xmm11, xmm2 ;# jzMa  -  jzMb  -
    movlhps xmm12, xmm3 ;# jzMc  -  jzMd -
    
    shufps xmm11, xmm12,  136  ;# 10001000 => jzMa jzMb jzMc jzMd

    ;# xmm9: jxMa jyMa jxMb jyMb 
    ;# xmm10: jxMc jyMc jxMd jyMd
    ;# xmm11: jzMa jzMb jzMc jzMd

    movaps xmm0, xmm1
    movaps xmm2, xmm1
    movaps xmm3, xmm4
    movaps xmm5, xmm4
    movaps xmm6, xmm7
    movaps xmm8, xmm7

	mulps xmm0, [rsp + nb214_dxH1M]
	mulps xmm1, [rsp + nb214_dyH1M]
	mulps xmm2, [rsp + nb214_dzH1M]
	mulps xmm3, [rsp + nb214_dxH2M]
	mulps xmm4, [rsp + nb214_dyH2M]
	mulps xmm5, [rsp + nb214_dzH2M]
	mulps xmm6, [rsp + nb214_dxMM]
	mulps xmm7, [rsp + nb214_dyMM]
	mulps xmm8, [rsp + nb214_dzMM]

    movaps xmm13, xmm0
    movaps xmm14, xmm1
    addps xmm11, xmm2
    addps xmm0, [rsp + nb214_fixH1]
    addps xmm1, [rsp + nb214_fiyH1]
    addps xmm2, [rsp + nb214_fizH1]

    addps xmm13, xmm3
    addps xmm14, xmm4
    addps xmm11, xmm5
    addps xmm3, [rsp + nb214_fixH2]
    addps xmm4, [rsp + nb214_fiyH2]
    addps xmm5, [rsp + nb214_fizH2]

    addps xmm13, xmm6
    addps xmm14, xmm7
    addps xmm11, xmm8
    addps xmm6, [rsp + nb214_fixM]
    addps xmm7, [rsp + nb214_fiyM]
    addps xmm8, [rsp + nb214_fizM]

    movaps [rsp + nb214_fixH1], xmm0
    movaps [rsp + nb214_fiyH1], xmm1
    movaps [rsp + nb214_fizH1], xmm2
    movaps [rsp + nb214_fixH2], xmm3
    movaps [rsp + nb214_fiyH2], xmm4
    movaps [rsp + nb214_fizH2], xmm5
    movaps [rsp + nb214_fixM], xmm6
    movaps [rsp + nb214_fiyM], xmm7
    movaps [rsp + nb214_fizM], xmm8
    
    ;# xmm0 = fMx
    ;# xmm1 = fMy
    ;# xmm2 = fMz
    movaps xmm0, xmm13
    unpcklps xmm13, xmm14
    unpckhps xmm0, xmm14
    
    addps xmm9, xmm13
    addps xmm10, xmm0

    movhlps  xmm12, xmm11 ;# fMzc fMzd
    
    movlps [rdi + rax*4 + 36], xmm9
    movhps [rdi + rbx*4 + 36], xmm9
    movlps [rdi + rcx*4 + 36], xmm10
    movhps [rdi + rdx*4 + 36], xmm10
    movss  [rdi + rax*4 + 44], xmm11
    movss  [rdi + rcx*4 + 44], xmm12
    shufps xmm11, xmm11, 1
    shufps xmm12, xmm12, 1
    movss  [rdi + rbx*4 + 44], xmm11
    movss  [rdi + rdx*4 + 44], xmm12
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb214_innerk],  4
	jl    .nb214_single_check
	jmp   .nb214_unroll_loop
.nb214_single_check:
	add dword ptr [rsp + nb214_innerk],  4
	jnz   .nb214_single_loop
	jmp   .nb214_updateouterdata
.nb214_single_loop:
	mov   rdx, [rsp + nb214_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb214_innerjjnr],  4	

	mov rsi, [rbp + nb214_pos]
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
	movaps [rsp + nb214_jxO], xmm6
	movaps [rsp + nb214_jyO], xmm3
	movaps [rsp + nb214_jzO], xmm1

    movaps xmm0, xmm6   ;# jxO
    movaps xmm2, xmm1   ;# jzO
    movaps xmm1, xmm3   ;# jyO
    movaps xmm4, xmm3   ;# jyO
    movaps xmm3, xmm6   ;# jxO
    movaps xmm5, xmm2   ;# jzO
        
	;# do O and M in parallel
	subps  xmm0, [rsp + nb214_ixO] 
	subps  xmm1, [rsp + nb214_iyO] 
	subps  xmm2, [rsp + nb214_izO] 
	subps  xmm3, [rsp + nb214_ixM] 
	subps  xmm4, [rsp + nb214_iyM]
	subps  xmm5, [rsp + nb214_izM]
	
	movaps [rsp + nb214_dxOO], xmm0
	movaps [rsp + nb214_dyOO], xmm1
	movaps [rsp + nb214_dzOO], xmm2
	movaps [rsp + nb214_dxMM], xmm3
	movaps [rsp + nb214_dyMM], xmm4
	movaps [rsp + nb214_dzMM], xmm5
	
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
	movaps [rsp + nb214_rsqMM], xmm4
	
	;# do 1/x for O and 1/sqrt(x) for M
	rcpss  xmm1, xmm0
	rsqrtps xmm5, xmm4
	movss  xmm2, [rsp + nb214_two]
	movaps  xmm6, xmm5	
	mulss  xmm0, xmm1
	mulps   xmm5, xmm5
	subss  xmm2, xmm0
	movaps  xmm7, [rsp + nb214_three]
	mulss  xmm2, xmm1 	;# 1/r2
	
	mulps   xmm5, xmm4
	movss  xmm0, xmm2
	subps   xmm7, xmm5
	mulss  xmm2, xmm2
	mulps   xmm7, xmm6
	mulss  xmm2, xmm0 	;# 1/r6
	mulps   xmm7, [rsp + nb214_half] ;# rinv iH1 - j water 
	movss  xmm1, xmm2
	movaps [rsp + nb214_rinvMM], xmm7

	mulss  xmm2, xmm2 	;# 1/r12
	mulss  xmm1, [rsp + nb214_c6]
	mulss  xmm2, [rsp + nb214_c12]
	movss  xmm3, xmm2
	subss  xmm3, xmm1
	addss  xmm3, [rsp + nb214_Vvdwtot]
	movss  [rsp + nb214_Vvdwtot], xmm3
	mulss  xmm1, [rsp + nb214_six]
	mulss  xmm2, [rsp + nb214_twelve]
	subss  xmm2, xmm1
	mulss  xmm0, xmm2 	;# fscal
	movss  xmm1, xmm0
	movss  xmm2, xmm0
	mulss  xmm0, [rsp + nb214_dxOO]
	mulss  xmm1, [rsp + nb214_dyOO]
	mulss  xmm2, [rsp + nb214_dzOO]
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	addss   xmm3, xmm0
	addss   xmm4, xmm1
	addss   xmm5, xmm2
	movaps  [rsp + nb214_fjxO], xmm3
	movaps  [rsp + nb214_fjyO], xmm4
	movaps  [rsp + nb214_fjzO], xmm5
	addss   xmm0, [rsp + nb214_fixO]
	addss   xmm1, [rsp + nb214_fiyO]
	addss   xmm2, [rsp + nb214_fizO]
	movss  [rsp + nb214_fixO], xmm0
	movss  [rsp + nb214_fiyO], xmm1
	movss  [rsp + nb214_fizO], xmm2

	;# do  M coulomb interaction
	movaps xmm0, [rsp + nb214_rinvMM]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb214_krf]
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [rsp + nb214_qqMH]
	movhps  xmm3, [rsp + nb214_qqMM]
	shufps  xmm3, xmm3, 193 ;# 11000001 

	mulps  xmm5, [rsp + nb214_rsqMM] ;# xmm5=krsq 
	movaps xmm6, xmm5
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [rsp + nb214_crf]
	mulps  xmm6, xmm3 ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	mulps xmm5, [rsp + nb214_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, xmm3 ;# xmm7 = coul part of fscal 
	
	addps  xmm6, [rsp + nb214_vctot] 
	movaps [rsp + nb214_vctot], xmm6
	mulps  xmm0, xmm7

	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [rsp + nb214_dxMM]
	mulps   xmm1, [rsp + nb214_dyMM]
	mulps   xmm2, [rsp + nb214_dzMM]
	;# update forces M - j water 
	movaps  xmm3, [rsp + nb214_fjxO]
	movaps  xmm4, [rsp + nb214_fjyO]
	movaps  xmm5, [rsp + nb214_fjzO]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  [rsp + nb214_fjxO], xmm3
	movaps  [rsp + nb214_fjyO], xmm4
	movaps  [rsp + nb214_fjzO], xmm5
	addps   xmm0, [rsp + nb214_fixM]
	addps   xmm1, [rsp + nb214_fiyM]
	addps   xmm2, [rsp + nb214_fizM]
	movaps  [rsp + nb214_fixM], xmm0
	movaps  [rsp + nb214_fiyM], xmm1
	movaps  [rsp + nb214_fizM], xmm2	
	
	;# i H1 & H2 simultaneously first get i particle coords: 
    movaps  xmm0, [rsp + nb214_jxO]
    movaps  xmm1, [rsp + nb214_jyO]
    movaps  xmm2, [rsp + nb214_jzO]
    movaps  xmm3, xmm0
    movaps  xmm4, xmm1
    movaps  xmm5, xmm2
	subps   xmm0, [rsp + nb214_ixH1]
	subps   xmm1, [rsp + nb214_iyH1]
	subps   xmm2, [rsp + nb214_izH1]
	subps   xmm3, [rsp + nb214_ixH2] 
	subps   xmm4, [rsp + nb214_iyH2] 
	subps   xmm5, [rsp + nb214_izH2] 
	movaps [rsp + nb214_dxH1H1], xmm0
	movaps [rsp + nb214_dyH1H1], xmm1
	movaps [rsp + nb214_dzH1H1], xmm2
	movaps [rsp + nb214_dxH2H2], xmm3
	movaps [rsp + nb214_dyH2H2], xmm4
	movaps [rsp + nb214_dzH2H2], xmm5
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
	movaps  [rsp + nb214_rsqH1H1], xmm0
	movaps  [rsp + nb214_rsqH2H2], xmm4
	
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb214_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb214_half] ;# rinvH1H1
	mulps   xmm7, [rsp + nb214_half] ;# rinvH2H2
	movaps  [rsp + nb214_rinvH1H1], xmm3
	movaps  [rsp + nb214_rinvH2H2], xmm7
	
	;# Do H1 coulomb interaction
	movaps xmm0, [rsp + nb214_rinvH1H1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb214_krf]
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [rsp + nb214_qqHH]
	movhps  xmm3, [rsp + nb214_qqMH]
	shufps  xmm3, xmm3, 193 ;# 11000001 

	mulps  xmm5, [rsp + nb214_rsqH1H1] ;# xmm5=krsq 
	movaps xmm6, xmm5
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [rsp + nb214_crf]
	mulps  xmm6, xmm3 ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	mulps xmm5, [rsp + nb214_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, xmm3 ;# xmm7 = coul part of fscal 
	
	addps  xmm6, [rsp + nb214_vctot] 
	movaps [rsp + nb214_vctot], xmm6

	mulps  xmm0, xmm7

	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [rsp + nb214_dxH1H1]
	mulps   xmm1, [rsp + nb214_dyH1H1]
	mulps   xmm2, [rsp + nb214_dzH1H1]
	;# update forces H1 - j water 
	movaps  xmm3, [rsp + nb214_fjxO]
	movaps  xmm4, [rsp + nb214_fjyO]
	movaps  xmm5, [rsp + nb214_fjzO]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  [rsp + nb214_fjxO], xmm3
	movaps  [rsp + nb214_fjyO], xmm4
	movaps  [rsp + nb214_fjzO], xmm5
	addps   xmm0, [rsp + nb214_fixH1]
	addps   xmm1, [rsp + nb214_fiyH1]
	addps   xmm2, [rsp + nb214_fizH1]
	movaps  [rsp + nb214_fixH1], xmm0
	movaps  [rsp + nb214_fiyH1], xmm1
	movaps  [rsp + nb214_fizH1], xmm2	

	;# H2 Coulomb
	movaps xmm0, [rsp + nb214_rinvH2H2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb214_krf]
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [rsp + nb214_qqHH]
	movhps  xmm3, [rsp + nb214_qqMH]
	shufps  xmm3, xmm3, 193 ;# 11000001 

	mulps  xmm5, [rsp + nb214_rsqH2H2] ;# xmm5=krsq 
	movaps xmm6, xmm5
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [rsp + nb214_crf]
	mulps  xmm6, xmm3 ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	mulps xmm5, [rsp + nb214_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, xmm3 ;# xmm7 = coul part of fscal 
	
	addps  xmm6, [rsp + nb214_vctot] ;# local vctot summation variable
	movaps [rsp + nb214_vctot], xmm6
	mulps  xmm0, xmm7

	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [rsp + nb214_dxH2H2]
	mulps   xmm1, [rsp + nb214_dyH2H2]
	mulps   xmm2, [rsp + nb214_dzH2H2]
	;# update forces H2 - j water 
	movaps  xmm3, [rsp + nb214_fjxO]
	movaps  xmm4, [rsp + nb214_fjyO]
	movaps  xmm5, [rsp + nb214_fjzO]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  [rsp + nb214_fjxO], xmm3
	movaps  [rsp + nb214_fjyO], xmm4
	movaps  [rsp + nb214_fjzO], xmm5
	addps   xmm0, [rsp + nb214_fixH2]
	addps   xmm1, [rsp + nb214_fiyH2]
	addps   xmm2, [rsp + nb214_fizH2]
	movaps  [rsp + nb214_fixH2], xmm0
	movaps  [rsp + nb214_fiyH2], xmm1
	movaps  [rsp + nb214_fizH2], xmm2			
	
	mov     rsi, [rbp + nb214_faction]
	;# update j water forces from local variables.
	;# transpose back first
	movaps  xmm0, [rsp + nb214_fjxO] ;# Ox H1x H2x Mx 
	movaps  xmm1, [rsp + nb214_fjyO] ;# Oy H1y H2y My
	movaps  xmm2, [rsp + nb214_fjzO] ;# Oz H1z H2z Mz
	 
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
	
	dec dword ptr [rsp + nb214_innerk]
	jz    .nb214_updateouterdata
	jmp   .nb214_single_loop
.nb214_updateouterdata:
	mov   ecx, [rsp + nb214_ii3]
	mov   rdi, [rbp + nb214_faction]
	mov   rsi, [rbp + nb214_fshift]
	mov   edx, [rsp + nb214_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb214_fixO]
	movaps xmm1, [rsp + nb214_fiyO] 
	movaps xmm2, [rsp + nb214_fizO]

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
	movaps xmm0, [rsp + nb214_fixH1]
	movaps xmm1, [rsp + nb214_fiyH1]
	movaps xmm2, [rsp + nb214_fizH1]

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
	movaps xmm0, [rsp + nb214_fixH2]
	movaps xmm1, [rsp + nb214_fiyH2]
	movaps xmm2, [rsp + nb214_fizH2]

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
	movaps xmm0, [rsp + nb214_fixM]
	movaps xmm1, [rsp + nb214_fiyM]
	movaps xmm2, [rsp + nb214_fizM]

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
	mov esi, [rsp + nb214_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb214_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb214_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb214_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb214_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb214_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb214_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb214_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb214_n], esi
        jmp .nb214_outer
.nb214_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb214_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb214_end
        ;# non-zero, do one more workunit
        jmp   .nb214_threadloop
.nb214_end:
	mov eax, [rsp + nb214_nouter]
	mov ebx, [rsp + nb214_ninner]
	mov rcx, [rbp + nb214_outeriter]
	mov rdx, [rbp + nb214_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 1904
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


	
.globl nb_kernel214nf_x86_64_sse
.globl _nb_kernel214nf_x86_64_sse
nb_kernel214nf_x86_64_sse:	
_nb_kernel214nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb214nf_fshift,         16
.equiv          nb214nf_gid,            24
.equiv          nb214nf_pos,            32
.equiv          nb214nf_faction,        40
.equiv          nb214nf_charge,         48
.equiv          nb214nf_p_facel,        56
.equiv          nb214nf_argkrf,         64
.equiv          nb214nf_argcrf,         72
.equiv          nb214nf_Vc,             80
.equiv          nb214nf_type,           88
.equiv          nb214nf_p_ntype,        96
.equiv          nb214nf_vdwparam,       104
.equiv          nb214nf_Vvdw,           112
.equiv          nb214nf_p_tabscale,     120
.equiv          nb214nf_VFtab,          128
.equiv          nb214nf_invsqrta,       136
.equiv          nb214nf_dvda,           144
.equiv          nb214nf_p_gbtabscale,   152
.equiv          nb214nf_GBtab,          160
.equiv          nb214nf_p_nthreads,     168
.equiv          nb214nf_count,          176
.equiv          nb214nf_mtx,            184
.equiv          nb214nf_outeriter,      192
.equiv          nb214nf_inneriter,      200
.equiv          nb214nf_work,           208
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
.equiv          nb214nf_nri,            912
.equiv          nb214nf_iinr,           920
.equiv          nb214nf_jindex,         928
.equiv          nb214nf_jjnr,           936
.equiv          nb214nf_shift,          944
.equiv          nb214nf_shiftvec,       952
.equiv          nb214nf_facel,          960
.equiv          nb214nf_innerk,         968
.equiv          nb214nf_n,              972
.equiv          nb214nf_nn1,            976
.equiv          nb214nf_nouter,         980
.equiv          nb214nf_ninner,         984

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
	sub rsp, 992		;# local variable stack space (n*16+8)
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
	mov [rsp + nb214nf_nouter], eax
	mov [rsp + nb214nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb214nf_nri], edi
	mov [rsp + nb214nf_iinr], rsi
	mov [rsp + nb214nf_jindex], rdx
	mov [rsp + nb214nf_jjnr], rcx
	mov [rsp + nb214nf_shift], r8
	mov [rsp + nb214nf_shiftvec], r9
	mov rsi, [rbp + nb214nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb214nf_facel], xmm0


	mov rsi, [rbp + nb214nf_argkrf]
	mov rdi, [rbp + nb214nf_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb214nf_krf], xmm1
	movaps [rsp + nb214nf_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb214nf_half], eax
	movss xmm1, [rsp + nb214nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb214nf_half],  xmm1
	movaps [rsp + nb214nf_two],  xmm2
	movaps [rsp + nb214nf_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb214nf_iinr]   ;# rcx = pointer into iinr[]
	mov   ebx, [rcx]		;# ebx =ii 

	mov   rdx, [rbp + nb214nf_charge]
	movss xmm5, [rdx + rbx*4 + 4]	
	movss xmm3, [rdx + rbx*4 + 12]	
	movss xmm4, xmm3	
	mov rsi, [rbp + nb214nf_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb214nf_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb214nf_qqMM], xmm3
	movaps [rsp + nb214nf_qqMH], xmm4
	movaps [rsp + nb214nf_qqHH], xmm5
	
	xorps xmm0, xmm0
	mov   rdx, [rbp + nb214nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb214nf_p_ntype]
	imul  ecx, [rdi]  ;# rcx = ntia = 2*ntype*type[ii0] 
	add   rdx, rcx
	mov   rax, [rbp + nb214nf_vdwparam]
	movlps xmm0, [rax + rdx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0x55
	movaps [rsp + nb214nf_c6], xmm0
	movaps [rsp + nb214nf_c12], xmm1

.nb214nf_threadloop:
        mov   rsi, [rbp + nb214nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb214nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   rbx, 1                           ;# rbx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb214nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb214nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb214nf_n], eax
        mov [rsp + nb214nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb214nf_outerstart
        jmp .nb214nf_end
	
.nb214nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb214nf_nouter]
	mov [rsp + nb214nf_nouter], ebx

.nb214nf_outer:
	mov   rax, [rsp + nb214nf_shift]  	;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb214nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb214nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb214nf_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb214nf_pos]	;# rax = base of pos[]  
	mov   [rsp + nb214nf_ii3], ebx	
	
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
	movaps [rsp + nb214nf_ixO], xmm3
	movaps [rsp + nb214nf_iyO], xmm4
	movaps [rsp + nb214nf_izO], xmm5
	movaps [rsp + nb214nf_ixH1], xmm6
	movaps [rsp + nb214nf_iyH1], xmm7

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
	movaps [rsp + nb214nf_izH1], xmm6
	movaps [rsp + nb214nf_ixH2], xmm0
	movaps [rsp + nb214nf_iyH2], xmm1
	movaps [rsp + nb214nf_izH2], xmm2
	movaps [rsp + nb214nf_ixM], xmm3
	movaps [rsp + nb214nf_iyM], xmm4
	movaps [rsp + nb214nf_izM], xmm5

	;# clear vctot 
	xorps xmm4, xmm4
	movaps [rsp + nb214nf_vctot], xmm4
	movaps [rsp + nb214nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb214nf_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb214nf_pos]
	mov   rdi, [rbp + nb214nf_faction]	
	mov   rax, [rsp + nb214nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb214nf_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb214nf_ninner]
	mov   [rsp + nb214nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb214nf_innerk], edx	;# number of innerloop atoms 
	jge   .nb214nf_unroll_loop
	jmp   .nb214nf_single_check
.nb214nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb214nf_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb214nf_innerjjnr],  16 ;# advance pointer (unroll 4) 

	mov rsi, [rbp + nb214nf_pos]   	;# base of pos[] 

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
	movaps [rsp + nb214nf_jxO], xmm0
	movaps [rsp + nb214nf_jyO], xmm1
	movaps [rsp + nb214nf_jzO], xmm2
	movaps [rsp + nb214nf_jxH1], xmm3

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
	movaps [rsp + nb214nf_jyH1], xmm0
	movaps [rsp + nb214nf_jzH1], xmm1
	movaps [rsp + nb214nf_jxH2], xmm2
	movaps [rsp + nb214nf_jyH2], xmm3

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
	movaps [rsp + nb214nf_jzH2], xmm0
	movaps [rsp + nb214nf_jxM], xmm1
	movaps [rsp + nb214nf_jyM], xmm2
	movaps [rsp + nb214nf_jzM], xmm3
	
	;# start calculating pairwise distances
	movaps xmm0, [rsp + nb214nf_ixO]
	movaps xmm1, [rsp + nb214nf_iyO]
	movaps xmm2, [rsp + nb214nf_izO]
	movaps xmm3, [rsp + nb214nf_ixH1]
	movaps xmm4, [rsp + nb214nf_iyH1]
	movaps xmm5, [rsp + nb214nf_izH1]
	subps  xmm0, [rsp + nb214nf_jxO]
	subps  xmm1, [rsp + nb214nf_jyO]
	subps  xmm2, [rsp + nb214nf_jzO]
	subps  xmm3, [rsp + nb214nf_jxH1]
	subps  xmm4, [rsp + nb214nf_jyH1]
	subps  xmm5, [rsp + nb214nf_jzH1]
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
	movaps [rsp + nb214nf_rsqOO], xmm0
	movaps [rsp + nb214nf_rsqH1H1], xmm3

	movaps xmm0, [rsp + nb214nf_ixH1]
	movaps xmm1, [rsp + nb214nf_iyH1]
	movaps xmm2, [rsp + nb214nf_izH1]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [rsp + nb214nf_jxH2]
	subps  xmm1, [rsp + nb214nf_jyH2]
	subps  xmm2, [rsp + nb214nf_jzH2]
	subps  xmm3, [rsp + nb214nf_jxM]
	subps  xmm4, [rsp + nb214nf_jyM]
	subps  xmm5, [rsp + nb214nf_jzM]
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
	movaps [rsp + nb214nf_rsqH1H2], xmm0
	movaps [rsp + nb214nf_rsqH1M], xmm3

	movaps xmm0, [rsp + nb214nf_ixH2]
	movaps xmm1, [rsp + nb214nf_iyH2]
	movaps xmm2, [rsp + nb214nf_izH2]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [rsp + nb214nf_jxH1]
	subps  xmm1, [rsp + nb214nf_jyH1]
	subps  xmm2, [rsp + nb214nf_jzH1]
	subps  xmm3, [rsp + nb214nf_jxH2]
	subps  xmm4, [rsp + nb214nf_jyH2]
	subps  xmm5, [rsp + nb214nf_jzH2]
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
	movaps [rsp + nb214nf_rsqH2H1], xmm0
	movaps [rsp + nb214nf_rsqH2H2], xmm3

	movaps xmm0, [rsp + nb214nf_ixH2]
	movaps xmm1, [rsp + nb214nf_iyH2]
	movaps xmm2, [rsp + nb214nf_izH2]
	movaps xmm3, [rsp + nb214nf_ixM]
	movaps xmm4, [rsp + nb214nf_iyM]
	movaps xmm5, [rsp + nb214nf_izM]
	subps  xmm0, [rsp + nb214nf_jxM]
	subps  xmm1, [rsp + nb214nf_jyM]
	subps  xmm2, [rsp + nb214nf_jzM]
	subps  xmm3, [rsp + nb214nf_jxH1]
	subps  xmm4, [rsp + nb214nf_jyH1]
	subps  xmm5, [rsp + nb214nf_jzH1]
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
	movaps [rsp + nb214nf_rsqH2M], xmm0
	movaps [rsp + nb214nf_rsqMH1], xmm4

	movaps xmm0, [rsp + nb214nf_ixM]
	movaps xmm1, [rsp + nb214nf_iyM]
	movaps xmm2, [rsp + nb214nf_izM]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [rsp + nb214nf_jxH2]
	subps  xmm1, [rsp + nb214nf_jyH2]
	subps  xmm2, [rsp + nb214nf_jzH2]
	subps  xmm3, [rsp + nb214nf_jxM]
	subps  xmm4, [rsp + nb214nf_jyM]
	subps  xmm5, [rsp + nb214nf_jzM]
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
	movaps [rsp + nb214nf_rsqMH2], xmm0
	movaps [rsp + nb214nf_rsqMM], xmm4

	;# start by doing reciprocal for OO
	movaps  xmm7, [rsp + nb214nf_rsqOO]
	rcpps   xmm2, xmm7
	movaps  xmm1, [rsp + nb214nf_two]
	mulps   xmm7, xmm2
	subps   xmm1, xmm7
	mulps   xmm2, xmm1 ;# rinvsq 
	movaps [rsp + nb214nf_rinvsqOO], xmm2
	
	;# next step is invsqrt - do two at a time.
	rsqrtps xmm1, [rsp + nb214nf_rsqH1H1]
	rsqrtps xmm5, [rsp + nb214nf_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb214nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb214nf_rsqH1H1]
	mulps   xmm5, [rsp + nb214nf_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb214nf_half] ;# rinvH1H1 
	mulps   xmm7, [rsp + nb214nf_half] ;# rinvH1H2 
	movaps  [rsp + nb214nf_rinvH1H1], xmm3
	movaps  [rsp + nb214nf_rinvH1H2], xmm7
	
	rsqrtps xmm1, [rsp + nb214nf_rsqH1M]
	rsqrtps xmm5, [rsp + nb214nf_rsqH2H1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb214nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb214nf_rsqH1M]
	mulps   xmm5, [rsp + nb214nf_rsqH2H1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb214nf_half] 
	mulps   xmm7, [rsp + nb214nf_half]
	movaps  [rsp + nb214nf_rinvH1M], xmm3
	movaps  [rsp + nb214nf_rinvH2H1], xmm7
	
	rsqrtps xmm1, [rsp + nb214nf_rsqH2H2]
	rsqrtps xmm5, [rsp + nb214nf_rsqH2M]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb214nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb214nf_rsqH2H2]
	mulps   xmm5, [rsp + nb214nf_rsqH2M]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb214nf_half] 
	mulps   xmm7, [rsp + nb214nf_half]
	movaps  [rsp + nb214nf_rinvH2H2], xmm3
	movaps  [rsp + nb214nf_rinvH2M], xmm7
	
	rsqrtps xmm1, [rsp + nb214nf_rsqMH1]
	rsqrtps xmm5, [rsp + nb214nf_rsqMH2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb214nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb214nf_rsqMH1]
	mulps   xmm5, [rsp + nb214nf_rsqMH2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb214nf_half] 
	mulps   xmm7, [rsp + nb214nf_half]
	movaps  [rsp + nb214nf_rinvMH1], xmm3
	movaps  [rsp + nb214nf_rinvMH2], xmm7
        		
	rsqrtps xmm1, [rsp + nb214nf_rsqMM]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb214nf_three]
	mulps   xmm1, [rsp + nb214nf_rsqMM]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + nb214nf_half] 
	movaps  [rsp + nb214nf_rinvMM], xmm3
	
	;# start with OO LJ interaction
	movaps xmm0, [rsp + nb214nf_rinvsqOO]
	movaps xmm1, xmm0
	mulps  xmm1, xmm1	;# rinv4
	mulps  xmm1, xmm0	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm1, [rsp + nb214nf_c6]
	mulps  xmm2, [rsp + nb214nf_c12]
	movaps xmm4, xmm2
	subps  xmm4, xmm1
	addps  xmm4, [rsp + nb214nf_Vvdwtot]
	movaps [rsp + nb214nf_Vvdwtot], xmm4

	;# Coulomb interactions 
	;# add all H-H rsq in xmm2, H-M rsq xmm4
	;# H-H rinv in xmm0, H-M in xmm1
	movaps xmm0, [rsp + nb214nf_rinvH1H1]
	movaps xmm1, [rsp + nb214nf_rinvH1M]
	movaps xmm2, [rsp + nb214nf_rsqH1H1]
	movaps xmm4, [rsp + nb214nf_rsqH1M]
	addps  xmm0, [rsp + nb214nf_rinvH1H2]
	addps  xmm1, [rsp + nb214nf_rinvH2M]
	addps  xmm2, [rsp + nb214nf_rsqH1H2]
	addps  xmm4, [rsp + nb214nf_rsqH2M]
	addps  xmm0, [rsp + nb214nf_rinvH2H1]
	addps  xmm1, [rsp + nb214nf_rinvMH1]
	addps  xmm2, [rsp + nb214nf_rsqH2H1]
	addps  xmm4, [rsp + nb214nf_rsqMH1]
	addps  xmm0, [rsp + nb214nf_rinvH2H2]
	addps  xmm1, [rsp + nb214nf_rinvMH2]
	addps  xmm2, [rsp + nb214nf_rsqH2H2]
	addps  xmm4, [rsp + nb214nf_rsqMH2]
	movaps xmm5, [rsp + nb214nf_krf]
	movaps xmm6, [rsp + nb214nf_crf]

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
	mulps  xmm0, [rsp + nb214nf_qqHH]
	mulps  xmm1, [rsp + nb214nf_qqMH]
	addps  xmm0, xmm1 
	addps  xmm0, [rsp + nb214nf_vctot]
	;# M-M interaction
	movaps xmm4, [rsp + nb214nf_rinvMM]
	mulps  xmm5, [rsp + nb214nf_rsqMM] ;# krsq
	addps  xmm5, xmm4                  ;# rinv+krsq
	subps xmm5, [rsp + nb214nf_crf] ;# xmm5=rinv+ krsq-crf 
	mulps xmm5, [rsp + nb214nf_qqMM]
	addps xmm5, xmm0
	movaps [rsp + nb214nf_vctot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb214nf_innerk],  4
	jl    .nb214nf_single_check
	jmp   .nb214nf_unroll_loop
.nb214nf_single_check:
	add dword ptr [rsp + nb214nf_innerk],  4
	jnz   .nb214nf_single_loop
	jmp   .nb214nf_updateouterdata
.nb214nf_single_loop:
	mov   rdx, [rsp + nb214nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb214nf_innerjjnr],  4	

	mov rsi, [rbp + nb214nf_pos]
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
	movaps [rsp + nb214nf_jxO], xmm6
	movaps [rsp + nb214nf_jyO], xmm3
	movaps [rsp + nb214nf_jzO], xmm1

	;# do O and M in parallel
	movaps xmm0, [rsp + nb214nf_ixO]
	movaps xmm1, [rsp + nb214nf_iyO]
	movaps xmm2, [rsp + nb214nf_izO]
	movaps xmm3, [rsp + nb214nf_ixM]
	movaps xmm4, [rsp + nb214nf_iyM]
	movaps xmm5, [rsp + nb214nf_izM]
	subps  xmm0, [rsp + nb214nf_jxO]
	subps  xmm1, [rsp + nb214nf_jyO]
	subps  xmm2, [rsp + nb214nf_jzO]
	subps  xmm3, [rsp + nb214nf_jxO]
	subps  xmm4, [rsp + nb214nf_jyO]
	subps  xmm5, [rsp + nb214nf_jzO]
	
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
	movaps [rsp + nb214nf_rsqMM], xmm4
	
	;# do 1/x for O and 1/sqrt(x) for M
	rcpss  xmm1, xmm0
	rsqrtps xmm5, xmm4
	movss  xmm2, [rsp + nb214nf_two]
	movaps  xmm6, xmm5	
	mulss  xmm0, xmm1
	mulps   xmm5, xmm5
	subss  xmm2, xmm0
	movaps  xmm7, [rsp + nb214nf_three]
	mulss  xmm2, xmm1 	;# 1/r2
	
	mulps   xmm5, xmm4
	movss  xmm0, xmm2
	subps   xmm7, xmm5
	mulss  xmm2, xmm2
	mulps   xmm7, xmm6
	mulss  xmm2, xmm0 	;# 1/r6
	mulps   xmm7, [rsp + nb214nf_half] ;# rinv iH1 - j water 
	movss  xmm1, xmm2
	movaps [rsp + nb214nf_rinvMM], xmm7

	mulss  xmm2, xmm2 	;# 1/r12
	mulss  xmm1, [rsp + nb214nf_c6]
	mulss  xmm2, [rsp + nb214nf_c12]
	movss  xmm3, xmm2
	subss  xmm3, xmm1
	addss  xmm3, [rsp + nb214nf_Vvdwtot]
	movss  [rsp + nb214nf_Vvdwtot], xmm3

	;# do  M coulomb interaction
	movaps xmm7, [rsp + nb214nf_rinvMM]
	movaps xmm6, [rsp + nb214nf_krf]

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [rsp + nb214nf_qqMH]
	movhps  xmm3, [rsp + nb214nf_qqMM]
	shufps  xmm3, xmm3, 193 ;# 11000001 

	mulps  xmm6, [rsp + nb214nf_rsqMM] ;# xmm6=krsq 
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [rsp + nb214nf_crf]
	mulps  xmm6, xmm3 ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addps  xmm6, [rsp + nb214nf_vctot] 
	movaps [rsp + nb214nf_vctot], xmm6
	
	;# i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [rsp + nb214nf_ixH1]
	movaps  xmm1, [rsp + nb214nf_iyH1]
	movaps  xmm2, [rsp + nb214nf_izH1]	
	movaps  xmm3, [rsp + nb214nf_ixH2] 
	movaps  xmm4, [rsp + nb214nf_iyH2] 
	movaps  xmm5, [rsp + nb214nf_izH2] 
	subps   xmm0, [rsp + nb214nf_jxO]
	subps   xmm1, [rsp + nb214nf_jyO]
	subps   xmm2, [rsp + nb214nf_jzO]
	subps   xmm3, [rsp + nb214nf_jxO]
	subps   xmm4, [rsp + nb214nf_jyO]
	subps   xmm5, [rsp + nb214nf_jzO]
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
	movaps  [rsp + nb214nf_rsqH1H1], xmm0
	movaps  [rsp + nb214nf_rsqH2H2], xmm4
	
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb214nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb214nf_half] ;# rinvH1H1
	mulps   xmm7, [rsp + nb214nf_half] ;# rinvH2H2
	movaps  [rsp + nb214nf_rinvH1H1], xmm3
	movaps  [rsp + nb214nf_rinvH2H2], xmm7
	
	;# Do H1 coulomb interaction
	movaps xmm7, [rsp + nb214nf_rinvH1H1]
	movaps xmm6, [rsp + nb214nf_krf]

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [rsp + nb214nf_qqHH]
	movhps  xmm3, [rsp + nb214nf_qqMH]
	shufps  xmm3, xmm3, 193 ;# 11000001 

	mulps  xmm6, [rsp + nb214nf_rsqH1H1] ;# xmm6=krsq 
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [rsp + nb214nf_crf]
	mulps  xmm6, xmm3 ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addps  xmm6, [rsp + nb214nf_vctot] 
	movaps [rsp + nb214nf_vctot], xmm6

	;# H2 Coulomb
	movaps xmm7, [rsp + nb214nf_rinvH2H2]
	movaps xmm6, [rsp + nb214nf_krf]

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [rsp + nb214nf_qqHH]
	movhps  xmm3, [rsp + nb214nf_qqMH]
	shufps  xmm3, xmm3, 193 ;# 11000001 

	mulps  xmm6, [rsp + nb214nf_rsqH2H2] ;# xmm6=krsq 
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [rsp + nb214nf_crf]
	mulps  xmm6, xmm3 ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addps  xmm6, [rsp + nb214nf_vctot] ;# local vctot summation variable
	movaps [rsp + nb214nf_vctot], xmm6
	dec dword ptr [rsp + nb214nf_innerk]
	jz    .nb214nf_updateouterdata
	jmp   .nb214nf_single_loop
.nb214nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb214nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb214nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb214nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb214nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb214nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb214nf_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb214nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb214nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb214nf_n], esi
        jmp .nb214nf_outer
.nb214nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb214nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb214nf_end
        ;# non-zero, do one more workunit
        jmp   .nb214nf_threadloop
.nb214nf_end:
	mov eax, [rsp + nb214nf_nouter]
	mov ebx, [rsp + nb214nf_ninner]
	mov rcx, [rbp + nb214nf_outeriter]
	mov rdx, [rbp + nb214nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 992
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
