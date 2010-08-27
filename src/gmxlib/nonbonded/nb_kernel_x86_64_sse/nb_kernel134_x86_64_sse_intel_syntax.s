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

	
.globl nb_kernel134_x86_64_sse
.globl _nb_kernel134_x86_64_sse
nb_kernel134_x86_64_sse:	
_nb_kernel134_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb134_fshift,           16
.equiv          nb134_gid,              24
.equiv          nb134_pos,              32
.equiv          nb134_faction,          40
.equiv          nb134_charge,           48
.equiv          nb134_p_facel,          56
.equiv          nb134_argkrf,           64
.equiv          nb134_argcrf,           72
.equiv          nb134_Vc,               80
.equiv          nb134_type,             88
.equiv          nb134_p_ntype,          96
.equiv          nb134_vdwparam,         104
.equiv          nb134_Vvdw,             112
.equiv          nb134_p_tabscale,       120
.equiv          nb134_VFtab,            128
.equiv          nb134_invsqrta,         136
.equiv          nb134_dvda,             144
.equiv          nb134_p_gbtabscale,     152
.equiv          nb134_GBtab,            160
.equiv          nb134_p_nthreads,       168
.equiv          nb134_count,            176
.equiv          nb134_mtx,              184
.equiv          nb134_outeriter,        192
.equiv          nb134_inneriter,        200
.equiv          nb134_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb134_ixO,              0
.equiv          nb134_iyO,              16
.equiv          nb134_izO,              32
.equiv          nb134_ixH1,             48
.equiv          nb134_iyH1,             64
.equiv          nb134_izH1,             80
.equiv          nb134_ixH2,             96
.equiv          nb134_iyH2,             112
.equiv          nb134_izH2,             128
.equiv          nb134_ixM,              144
.equiv          nb134_iyM,              160
.equiv          nb134_izM,              176
.equiv          nb134_jxO,              192
.equiv          nb134_jyO,              208
.equiv          nb134_jzO,              224
.equiv          nb134_jxH1,             240
.equiv          nb134_jyH1,             256
.equiv          nb134_jzH1,             272
.equiv          nb134_jxH2,             288
.equiv          nb134_jyH2,             304
.equiv          nb134_jzH2,             320
.equiv          nb134_jxM,              336
.equiv          nb134_jyM,              352
.equiv          nb134_jzM,              368
.equiv          nb134_dxOO,             384
.equiv          nb134_dyOO,             400
.equiv          nb134_dzOO,             416
.equiv          nb134_dxH1H1,           432
.equiv          nb134_dyH1H1,           448
.equiv          nb134_dzH1H1,           464
.equiv          nb134_dxH1H2,           480
.equiv          nb134_dyH1H2,           496
.equiv          nb134_dzH1H2,           512
.equiv          nb134_dxH1M,            528
.equiv          nb134_dyH1M,            544
.equiv          nb134_dzH1M,            560
.equiv          nb134_dxH2H1,           576
.equiv          nb134_dyH2H1,           592
.equiv          nb134_dzH2H1,           608
.equiv          nb134_dxH2H2,           624
.equiv          nb134_dyH2H2,           640
.equiv          nb134_dzH2H2,           656
.equiv          nb134_dxH2M,            672
.equiv          nb134_dyH2M,            688
.equiv          nb134_dzH2M,            704
.equiv          nb134_dxMH1,            720
.equiv          nb134_dyMH1,            736
.equiv          nb134_dzMH1,            752
.equiv          nb134_dxMH2,            768
.equiv          nb134_dyMH2,            784
.equiv          nb134_dzMH2,            800
.equiv          nb134_dxMM,             816
.equiv          nb134_dyMM,             832
.equiv          nb134_dzMM,             848
.equiv          nb134_qqMM,             864
.equiv          nb134_qqMH,             880
.equiv          nb134_qqHH,             896
.equiv          nb134_two,              912
.equiv          nb134_c6,               928
.equiv          nb134_c12,              944
.equiv          nb134_tsc,              960
.equiv          nb134_fstmp,            976
.equiv          nb134_vctot,            992
.equiv          nb134_Vvdwtot,          1008
.equiv          nb134_fixO,             1024
.equiv          nb134_fiyO,             1040
.equiv          nb134_fizO,             1056
.equiv          nb134_fixH1,            1072
.equiv          nb134_fiyH1,            1088
.equiv          nb134_fizH1,            1104
.equiv          nb134_fixH2,            1120
.equiv          nb134_fiyH2,            1136
.equiv          nb134_fizH2,            1152
.equiv          nb134_fixM,             1168
.equiv          nb134_fiyM,             1184
.equiv          nb134_fizM,             1200
.equiv          nb134_fjxO,             1216
.equiv          nb134_fjyO,             1232
.equiv          nb134_fjzO,             1248
.equiv          nb134_fjxH1,            1264
.equiv          nb134_fjyH1,            1280
.equiv          nb134_fjzH1,            1296
.equiv          nb134_fjxH2,            1312
.equiv          nb134_fjyH2,            1328
.equiv          nb134_fjzH2,            1344
.equiv          nb134_fjxM,             1360
.equiv          nb134_fjyM,             1376
.equiv          nb134_fjzM,             1392
.equiv          nb134_half,             1408
.equiv          nb134_three,            1424
.equiv          nb134_rsqOO,            1440
.equiv          nb134_rsqH1H1,          1456
.equiv          nb134_rsqH1H2,          1472
.equiv          nb134_rsqH1M,           1488
.equiv          nb134_rsqH2H1,          1504
.equiv          nb134_rsqH2H2,          1520
.equiv          nb134_rsqH2M,           1536
.equiv          nb134_rsqMH1,           1552
.equiv          nb134_rsqMH2,           1568
.equiv          nb134_rsqMM,            1584
.equiv          nb134_rinvOO,           1600
.equiv          nb134_rinvH1H1,         1616
.equiv          nb134_rinvH1H2,         1632
.equiv          nb134_rinvH1M,          1648
.equiv          nb134_rinvH2H1,         1664
.equiv          nb134_rinvH2H2,         1680
.equiv          nb134_rinvH2M,          1696
.equiv          nb134_rinvMH1,          1712
.equiv          nb134_rinvMH2,          1728
.equiv          nb134_rinvMM,           1744
.equiv          nb134_krf,              1776
.equiv          nb134_crf,              1792
.equiv          nb134_is3,              1808
.equiv          nb134_ii3,              1812
.equiv          nb134_innerjjnr,        1816
.equiv          nb134_nri,              1824
.equiv          nb134_iinr,             1832
.equiv          nb134_jindex,           1840
.equiv          nb134_jjnr,             1848
.equiv          nb134_shift,            1856
.equiv          nb134_shiftvec,         1864
.equiv          nb134_facel,            1872
.equiv          nb134_innerk,           1880
.equiv          nb134_n,                1884
.equiv          nb134_nn1,              1888
.equiv          nb134_nouter,           1892
.equiv          nb134_ninner,           1896

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
	mov [rsp + nb134_nouter], eax
	mov [rsp + nb134_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb134_nri], edi
	mov [rsp + nb134_iinr], rsi
	mov [rsp + nb134_jindex], rdx
	mov [rsp + nb134_jjnr], rcx
	mov [rsp + nb134_shift], r8
	mov [rsp + nb134_shiftvec], r9
	mov rsi, [rbp + nb134_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb134_facel], xmm0

	mov rax, [rbp + nb134_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb134_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb134_half], eax
	movss xmm1, [rsp + nb134_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb134_half],  xmm1
	movaps [rsp + nb134_two],  xmm2
	movaps [rsp + nb134_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb134_iinr]   ;# rcx = pointer into iinr[]
	mov   ebx, [rcx]		;# ebx =ii 

	mov   rdx, [rbp + nb134_charge]
	movss xmm5, [rdx + rbx*4 + 4]	
	movss xmm3, [rdx + rbx*4 + 12]	
	movss xmm4, xmm3	
	mov rsi, [rbp + nb134_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb134_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb134_qqMM], xmm3
	movaps [rsp + nb134_qqMH], xmm4
	movaps [rsp + nb134_qqHH], xmm5
	
	xorps xmm0, xmm0
	mov   rdx, [rbp + nb134_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb134_p_ntype]
	imul  ecx, [rdi]  ;# rcx = ntia = 2*ntype*type[ii0] 
	add   rdx, rcx
	mov   rax, [rbp + nb134_vdwparam]
	movlps xmm0, [rax + rdx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0x55
	movaps [rsp + nb134_c6], xmm0
	movaps [rsp + nb134_c12], xmm1

.nb134_threadloop:
        mov   rsi, [rbp + nb134_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb134_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb134_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb134_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb134_n], eax
        mov [rsp + nb134_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb134_outerstart
        jmp .nb134_end
	
.nb134_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb134_nouter]
	mov [rsp + nb134_nouter], ebx

.nb134_outer:
	mov   rax, [rsp + nb134_shift]  	;# eax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb134_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb134_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb134_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb134_pos]	;# eax = base of pos[]  
	mov   [rsp + nb134_ii3], ebx	
	
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
	movaps [rsp + nb134_ixO], xmm3
	movaps [rsp + nb134_iyO], xmm4
	movaps [rsp + nb134_izO], xmm5
	movaps [rsp + nb134_ixH1], xmm6
	movaps [rsp + nb134_iyH1], xmm7

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
	movaps [rsp + nb134_izH1], xmm6
	movaps [rsp + nb134_ixH2], xmm0
	movaps [rsp + nb134_iyH2], xmm1
	movaps [rsp + nb134_izH2], xmm2
	movaps [rsp + nb134_ixM], xmm3
	movaps [rsp + nb134_iyM], xmm4
	movaps [rsp + nb134_izM], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb134_vctot], xmm4
	movaps [rsp + nb134_Vvdwtot], xmm4
	movaps [rsp + nb134_fixO], xmm4
	movaps [rsp + nb134_fiyO], xmm4
	movaps [rsp + nb134_fizO], xmm4
	movaps [rsp + nb134_fixH1], xmm4
	movaps [rsp + nb134_fiyH1], xmm4
	movaps [rsp + nb134_fizH1], xmm4
	movaps [rsp + nb134_fixH2], xmm4
	movaps [rsp + nb134_fiyH2], xmm4
	movaps [rsp + nb134_fizH2], xmm4
	movaps [rsp + nb134_fixM], xmm4
	movaps [rsp + nb134_fiyM], xmm4
	movaps [rsp + nb134_fizM], xmm4
	
	mov   rax, [rsp + nb134_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb134_pos]
	mov   rdi, [rbp + nb134_faction]	
	mov   rax, [rsp + nb134_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb134_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb134_ninner]
	mov   [rsp + nb134_ninner], ecx
	add   edx, 0
	mov   [rsp + nb134_innerk], edx	;# number of innerloop atoms 
	jge   .nb134_unroll_loop
	jmp   .nb134_single_check
.nb134_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb134_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb134_innerjjnr],  16 ;# advance pointer (unroll 4) 

	mov rsi, [rbp + nb134_pos]   	;# base of pos[] 

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
        
    subps xmm0, [rsp + nb134_ixO]
    subps xmm1, [rsp + nb134_iyO]
    subps xmm2, [rsp + nb134_izO]

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
	movaps xmm4, [rsp + nb134_three]
	mulps xmm5, xmm1	;# rsq*lu*lu 	
    subps xmm4, xmm5	;# 30-rsq*lu*lu 
	mulps xmm4, xmm2	
	mulps xmm4, [rsp + nb134_half]	
	movaps xmm2, xmm4
	mulps  xmm1, xmm4	
    ;# xmm2=rinv
    ;# xmm1=r

    mulps xmm1, [rsp + nb134_tsc] ;# rtab

    ;# truncate and convert to integers
    cvttps2dq xmm5, xmm1
    
    ;# convert back to float
    cvtdq2ps  xmm4, xmm5
    
    ;# multiply by 8
    pslld   xmm5, 3

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
    
	mov rsi, [rbp + nb134_VFtab]
    ;# load LJ dispersion and repulsion in parallel
    movlps xmm5, [rsi + r8*4]
   	movlps xmm9, [rsi + r8*4 + 16]
	movlps xmm7,  [rsi + r10*4]
	movlps xmm11, [rsi + r10*4 + 16]
	movhps xmm5, [rsi + r9*4]
	movhps xmm9, [rsi + r9*4 + 16]
	movhps xmm7,  [rsi + r11*4]
	movhps xmm11, [rsi + r11*4 + 16]

    movaps xmm4, xmm5
    movaps xmm8, xmm9
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm8, xmm11, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101
	shufps xmm9, xmm11, 221  ;# 11011101

	movlps xmm7,  [rsi + r8*4 + 8]
	movlps xmm11, [rsi + r8*4 + 24]
	movlps xmm0,  [rsi + r10*4 + 8]
	movlps xmm3,  [rsi + r10*4 + 24]
	movhps xmm7,  [rsi + r9*4 + 8]
	movhps xmm11, [rsi + r9*4 + 24]
	movhps xmm0,  [rsi + r11*4 + 8]
	movhps xmm3,  [rsi + r11*4 + 24]

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

    mulps  xmm5, [rsp + nb134_c6]  ;# VV*c6 = vnb6
    mulps  xmm9, [rsp + nb134_c12]  ;# VV*c12 = vnb12
    addps  xmm5, xmm9
    addps  xmm5, [rsp + nb134_Vvdwtot]
    movaps [rsp + nb134_Vvdwtot], xmm5
        
    mulps  xmm7, [rsp + nb134_c6]   ;# FF*c6 = fnb6
    mulps  xmm11, [rsp + nb134_c12]   ;# FF*c12  = fnb12
    addps  xmm7, xmm11
    
    mulps  xmm7, [rsp + nb134_tsc]
    mulps  xmm7, xmm2   ;# fscal
    xorps  xmm9, xmm9
    
    subps  xmm9, xmm7

    ;# fx/fy/fz
    mulps  xmm13, xmm9
    mulps  xmm14, xmm9
    mulps  xmm15, xmm9

    ;# increment i force
    movaps xmm0, [rsp + nb134_fixO]
    movaps xmm1, [rsp + nb134_fiyO]
    movaps xmm2, [rsp + nb134_fizO]
    addps  xmm0, xmm13
    addps  xmm1, xmm14
    addps  xmm2, xmm15
    movaps [rsp + nb134_fixO], xmm0
    movaps [rsp + nb134_fiyO], xmm1
    movaps [rsp + nb134_fizO], xmm2

	;# move j O forces to local temp variables 
	mov rdi, [rbp + nb134_faction]
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
	mov rsi, [rbp + nb134_pos]       ;# base of pos[] 
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
    
    subps xmm0, [rsp + nb134_ixH1]
    subps xmm1, [rsp + nb134_iyH1]
    subps xmm2, [rsp + nb134_izH1]
    subps xmm3, [rsp + nb134_ixH2]
    subps xmm4, [rsp + nb134_iyH2]
    subps xmm5, [rsp + nb134_izH2]
    subps xmm6, [rsp + nb134_ixM]
    subps xmm7, [rsp + nb134_iyM]
    subps xmm8, [rsp + nb134_izM]
    
	movaps [rsp + nb134_dxH1H1], xmm0
	movaps [rsp + nb134_dyH1H1], xmm1
	movaps [rsp + nb134_dzH1H1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb134_dxH2H1], xmm3
	movaps [rsp + nb134_dyH2H1], xmm4
	movaps [rsp + nb134_dzH2H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb134_dxMH1], xmm6
	movaps [rsp + nb134_dyMH1], xmm7
	movaps [rsp + nb134_dzMH1], xmm8
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
		
	movaps  xmm9, [rsp + nb134_three]
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

	movaps  xmm0, [rsp + nb134_half]
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
    mulps  xmm0, [rsp + nb134_qqHH] 
    mulps  xmm1, [rsp + nb134_qqHH] 
    mulps  xmm2, [rsp + nb134_qqMH] 
    mulps  xmm9, xmm0
    mulps  xmm10, xmm1
    mulps  xmm11, xmm2
    
    addps xmm0, [rsp + nb134_vctot] 
    addps xmm1, xmm2
    addps xmm0, xmm1
    movaps [rsp + nb134_vctot], xmm0
    
	;# move j H1 forces to local temp variables 
	mov rdi, [rbp + nb134_faction]
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

	mulps xmm7, [rsp + nb134_dxH1H1]
	mulps xmm8, [rsp + nb134_dyH1H1]
	mulps xmm9, [rsp + nb134_dzH1H1]
	mulps xmm10, [rsp + nb134_dxH2H1]
	mulps xmm11, [rsp + nb134_dyH2H1]
	mulps xmm12, [rsp + nb134_dzH2H1]
	mulps xmm13, [rsp + nb134_dxMH1]
	mulps xmm14, [rsp + nb134_dyMH1]
	mulps xmm15, [rsp + nb134_dzMH1]

    movaps xmm3, xmm7
    movaps xmm4, xmm8
    addps xmm2, xmm9
    addps xmm7, [rsp + nb134_fixH1]
    addps xmm8, [rsp + nb134_fiyH1]
    addps xmm9, [rsp + nb134_fizH1]

    addps xmm3, xmm10
    addps xmm4, xmm11
    addps xmm2, xmm12
    addps xmm10, [rsp + nb134_fixH2]
    addps xmm11, [rsp + nb134_fiyH2]
    addps xmm12, [rsp + nb134_fizH2]

    addps xmm3, xmm13
    addps xmm4, xmm14
    addps xmm2, xmm15
    addps xmm13, [rsp + nb134_fixM]
    addps xmm14, [rsp + nb134_fiyM]
    addps xmm15, [rsp + nb134_fizM]

    movaps [rsp + nb134_fixH1], xmm7
    movaps [rsp + nb134_fiyH1], xmm8
    movaps [rsp + nb134_fizH1], xmm9
    movaps [rsp + nb134_fixH2], xmm10
    movaps [rsp + nb134_fiyH2], xmm11
    movaps [rsp + nb134_fizH2], xmm12
    movaps [rsp + nb134_fixM], xmm13
    movaps [rsp + nb134_fiyM], xmm14
    movaps [rsp + nb134_fizM], xmm15
    
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
	mov rsi, [rbp + nb134_pos]       ;# base of pos[] 
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
    
    subps xmm0, [rsp + nb134_ixH1]
    subps xmm1, [rsp + nb134_iyH1]
    subps xmm2, [rsp + nb134_izH1]
    subps xmm3, [rsp + nb134_ixH2]
    subps xmm4, [rsp + nb134_iyH2]
    subps xmm5, [rsp + nb134_izH2]
    subps xmm6, [rsp + nb134_ixM]
    subps xmm7, [rsp + nb134_iyM]
    subps xmm8, [rsp + nb134_izM]
    
	movaps [rsp + nb134_dxH1H2], xmm0
	movaps [rsp + nb134_dyH1H2], xmm1
	movaps [rsp + nb134_dzH1H2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb134_dxH2H2], xmm3
	movaps [rsp + nb134_dyH2H2], xmm4
	movaps [rsp + nb134_dzH2H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb134_dxMH2], xmm6
	movaps [rsp + nb134_dyMH2], xmm7
	movaps [rsp + nb134_dzMH2], xmm8
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
		
	movaps  xmm9, [rsp + nb134_three]
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

	movaps  xmm0, [rsp + nb134_half]
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
    mulps  xmm0, [rsp + nb134_qqHH] 
    mulps  xmm1, [rsp + nb134_qqHH] 
    mulps  xmm2, [rsp + nb134_qqMH] 
    mulps  xmm9, xmm0
    mulps  xmm10, xmm1
    mulps  xmm11, xmm2
    
    addps xmm0, [rsp + nb134_vctot] 
    addps xmm1, xmm2
    addps xmm0, xmm1
    movaps [rsp + nb134_vctot], xmm0
    
	;# move j H2 forces to local temp variables 
	mov rdi, [rbp + nb134_faction]
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

	mulps xmm7, [rsp + nb134_dxH1H2]
	mulps xmm8, [rsp + nb134_dyH1H2]
	mulps xmm9, [rsp + nb134_dzH1H2]
	mulps xmm10, [rsp + nb134_dxH2H2]
	mulps xmm11, [rsp + nb134_dyH2H2]
	mulps xmm12, [rsp + nb134_dzH2H2]
	mulps xmm13, [rsp + nb134_dxMH2]
	mulps xmm14, [rsp + nb134_dyMH2]
	mulps xmm15, [rsp + nb134_dzMH2]

    movaps xmm3, xmm7
    movaps xmm4, xmm8
    addps xmm2, xmm9
    addps xmm7, [rsp + nb134_fixH1]
    addps xmm8, [rsp + nb134_fiyH1]
    addps xmm9, [rsp + nb134_fizH1]

    addps xmm3, xmm10
    addps xmm4, xmm11
    addps xmm2, xmm12
    addps xmm10, [rsp + nb134_fixH2]
    addps xmm11, [rsp + nb134_fiyH2]
    addps xmm12, [rsp + nb134_fizH2]

    addps xmm3, xmm13
    addps xmm4, xmm14
    addps xmm2, xmm15
    addps xmm13, [rsp + nb134_fixM]
    addps xmm14, [rsp + nb134_fiyM]
    addps xmm15, [rsp + nb134_fizM]
    
    movaps [rsp + nb134_fixH1], xmm7
    movaps [rsp + nb134_fiyH1], xmm8
    movaps [rsp + nb134_fizH1], xmm9
    movaps [rsp + nb134_fixH2], xmm10
    movaps [rsp + nb134_fiyH2], xmm11
    movaps [rsp + nb134_fizH2], xmm12
    movaps [rsp + nb134_fixM], xmm13
    movaps [rsp + nb134_fiyM], xmm14
    movaps [rsp + nb134_fizM], xmm15

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
	mov rsi, [rbp + nb134_pos]       ;# base of pos[] 
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
    
    subps xmm0, [rsp + nb134_ixH1]
    subps xmm1, [rsp + nb134_iyH1]
    subps xmm2, [rsp + nb134_izH1]
    subps xmm3, [rsp + nb134_ixH2]
    subps xmm4, [rsp + nb134_iyH2]
    subps xmm5, [rsp + nb134_izH2]
    subps xmm6, [rsp + nb134_ixM]
    subps xmm7, [rsp + nb134_iyM]
    subps xmm8, [rsp + nb134_izM]
    
	movaps [rsp + nb134_dxH1M], xmm0
	movaps [rsp + nb134_dyH1M], xmm1
	movaps [rsp + nb134_dzH1M], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb134_dxH2M], xmm3
	movaps [rsp + nb134_dyH2M], xmm4
	movaps [rsp + nb134_dzH2M], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb134_dxMM], xmm6
	movaps [rsp + nb134_dyMM], xmm7
	movaps [rsp + nb134_dzMM], xmm8
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
		
	movaps  xmm9, [rsp + nb134_three]
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

	movaps  xmm0, [rsp + nb134_half]
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
    mulps  xmm0, [rsp + nb134_qqMH] 
    mulps  xmm1, [rsp + nb134_qqMH] 
    mulps  xmm2, [rsp + nb134_qqMM] 
    mulps  xmm9, xmm0
    mulps  xmm10, xmm1
    mulps  xmm11, xmm2
    
    addps xmm0, [rsp + nb134_vctot] 
    addps xmm1, xmm2
    addps xmm0, xmm1
    movaps [rsp + nb134_vctot], xmm0
    
	;# move j M forces to local temp variables 
	mov rdi, [rbp + nb134_faction]
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

	mulps xmm7, [rsp + nb134_dxH1M]
	mulps xmm8, [rsp + nb134_dyH1M]
	mulps xmm9, [rsp + nb134_dzH1M]
	mulps xmm10, [rsp + nb134_dxH2M]
	mulps xmm11, [rsp + nb134_dyH2M]
	mulps xmm12, [rsp + nb134_dzH2M]
	mulps xmm13, [rsp + nb134_dxMM]
	mulps xmm14, [rsp + nb134_dyMM]
	mulps xmm15, [rsp + nb134_dzMM]

    movaps xmm3, xmm7
    movaps xmm4, xmm8
    addps xmm2, xmm9
    addps xmm7, [rsp + nb134_fixH1]
    addps xmm8, [rsp + nb134_fiyH1]
    addps xmm9, [rsp + nb134_fizH1]

    addps xmm3, xmm10
    addps xmm4, xmm11
    addps xmm2, xmm12
    addps xmm10, [rsp + nb134_fixH2]
    addps xmm11, [rsp + nb134_fiyH2]
    addps xmm12, [rsp + nb134_fizH2]

    addps xmm3, xmm13
    addps xmm4, xmm14
    addps xmm2, xmm15
    addps xmm13, [rsp + nb134_fixM]
    addps xmm14, [rsp + nb134_fiyM]
    addps xmm15, [rsp + nb134_fizM]
    
    movaps [rsp + nb134_fixH1], xmm7
    movaps [rsp + nb134_fiyH1], xmm8
    movaps [rsp + nb134_fizH1], xmm9
    movaps [rsp + nb134_fixH2], xmm10
    movaps [rsp + nb134_fiyH2], xmm11
    movaps [rsp + nb134_fizH2], xmm12
    movaps [rsp + nb134_fixM], xmm13
    movaps [rsp + nb134_fiyM], xmm14
    movaps [rsp + nb134_fizM], xmm15

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
	sub dword ptr [rsp + nb134_innerk],  4
	jl    .nb134_single_check
	jmp   .nb134_unroll_loop
.nb134_single_check:
	add dword ptr [rsp + nb134_innerk],  4
	jnz   .nb134_single_loop
	jmp   .nb134_updateouterdata
.nb134_single_loop:
	mov   rdx, [rsp + nb134_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb134_innerjjnr],  4	

	mov rsi, [rbp + nb134_pos]
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
	movaps [rsp + nb134_jxO], xmm6
	movaps [rsp + nb134_jyO], xmm3
	movaps [rsp + nb134_jzO], xmm1

    movaps xmm0, xmm6   ;# jxO
    movaps xmm2, xmm1   ;# jzO
    movaps xmm1, xmm3   ;# jyO
    movaps xmm4, xmm3   ;# jyO
    movaps xmm3, xmm6   ;# jxO
    movaps xmm5, xmm2   ;# jzO
        
	;# do O and M in parallel
	subps  xmm0, [rsp + nb134_ixO] 
	subps  xmm1, [rsp + nb134_iyO] 
	subps  xmm2, [rsp + nb134_izO] 
	subps  xmm3, [rsp + nb134_ixM] 
	subps  xmm4, [rsp + nb134_iyM]
	subps  xmm5, [rsp + nb134_izM]
	
	movaps [rsp + nb134_dxOO], xmm0
	movaps [rsp + nb134_dyOO], xmm1
	movaps [rsp + nb134_dzOO], xmm2
	movaps [rsp + nb134_dxMM], xmm3
	movaps [rsp + nb134_dyMM], xmm4
	movaps [rsp + nb134_dzMM], xmm5
	
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
	movaps [rsp + nb134_rsqOO], xmm0
	movaps [rsp + nb134_rsqMM], xmm4
	
	;# do 1/x for O
	rsqrtps xmm1, xmm0
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb134_three]
	mulps   xmm1, xmm0	;# rsq*lu*lu 
	subps   xmm3, xmm1	;# constant 30-rsq*lu*lu 
	mulps   xmm3, xmm2	;# lu*(3-rsq*lu*lu) 
	mulps   xmm3, [rsp + nb134_half]
	movaps  [rsp + nb134_rinvOO], xmm3	;# rinvH2 
	
	;# 1/sqrt(x) for M
	rsqrtps xmm5, xmm4
	movaps  xmm6, xmm5	
	mulps   xmm5, xmm5
	movaps  xmm7, [rsp + nb134_three]
	mulps   xmm5, xmm4
	subps   xmm7, xmm5
	mulps   xmm7, xmm6
	mulps   xmm7, [rsp + nb134_half] ;# rinv iH1 - j water 
	movaps [rsp + nb134_rinvMM], xmm7


	;# LJ table interaction
	movaps xmm0, [rsp + nb134_rinvOO]
	movss xmm1, xmm0
	mulss  xmm1, [rsp + nb134_rsqOO] ;# xmm1=r 
	mulss  xmm1, [rsp + nb134_tsc]
		
    cvttps2pi mm6, xmm1
    cvtpi2ps xmm3, mm6
	subss    xmm1, xmm3	;# xmm1=eps 
    movss xmm2, xmm1
    mulss  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 3
		
    mov  rsi, [rbp + nb134_VFtab]
    movd r8d, mm6
	
    ;# dispersion 
    movlps xmm5, [rsi + r8*4]
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# constant 10001000
    shufps xmm5, xmm7, 221  ;# constant 11011101

    movlps xmm7, [rsi + r8*4 + 8]
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# constant 10001000
    shufps xmm7, xmm3, 221  ;# constant 11011101
    ;# dispersion table ready, in xmm4-xmm7 
    mulss  xmm6, xmm1       ;# xmm6=Geps 
    mulss  xmm7, xmm2       ;# xmm7=Heps2 
    addss  xmm5, xmm6
    addss  xmm5, xmm7       ;# xmm5=Fp 
    mulss  xmm7, [rsp + nb134_two]       ;# two*Heps2 
    addss  xmm7, xmm6
    addss  xmm7, xmm5 ;# xmm7=FF 
    mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
    addss  xmm5, xmm4 ;# xmm5=VV 

    movss xmm4, [rsp + nb134_c6]
    mulss  xmm7, xmm4    ;# fijD 
    mulss  xmm5, xmm4    ;# Vvdw6 
	xorps  xmm3, xmm3
	
	mulps  xmm7, [rsp + nb134_tsc]
	subss  xmm3, xmm7
	movss  [rsp + nb134_fstmp], xmm3 

    addss  xmm5, [rsp + nb134_Vvdwtot]
    movss [rsp + nb134_Vvdwtot], xmm5

    ;# repulsion 
    movlps xmm5, [rsi + r8*4 + 16]
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# constant 10001000
    shufps xmm5, xmm7, 221  ;# constant 11011101

    movlps xmm7, [rsi + r8*4 + 24]
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# constant 10001000
    shufps xmm7, xmm3, 221  ;# constant 11011101
    ;# table ready, in xmm4-xmm7 
    mulss  xmm6, xmm1       ;# xmm6=Geps 
    mulss  xmm7, xmm2       ;# xmm7=Heps2 
    addss  xmm5, xmm6
    addss  xmm5, xmm7       ;# xmm5=Fp 
    mulss  xmm7, [rsp + nb134_two]       ;# two*Heps2 
    addss  xmm7, xmm6
    addss  xmm7, xmm5 ;# xmm7=FF 
    mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
    addss  xmm5, xmm4 ;# xmm5=VV 

    movss xmm4, [rsp + nb134_c12]
    mulss  xmm7, xmm4 ;# fijR 
    mulss  xmm5, xmm4 ;# Vvdw12 
	movaps xmm3, [rsp + nb134_fstmp]
	mulss  xmm7, [rsp + nb134_tsc]
	subss  xmm3, xmm7

    addss  xmm5, [rsp + nb134_Vvdwtot]
    movss [rsp + nb134_Vvdwtot], xmm5

	mulss xmm0, xmm3
	movaps xmm1, xmm0
	movaps xmm2, xmm0
		
	mulss  xmm0, [rsp + nb134_dxOO]
	mulss  xmm1, [rsp + nb134_dyOO]
	mulss  xmm2, [rsp + nb134_dzOO]
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	addss   xmm3, xmm0
	addss   xmm4, xmm1
	addss   xmm5, xmm2
	movaps  [rsp + nb134_fjxO], xmm3
	movaps  [rsp + nb134_fjyO], xmm4
	movaps  [rsp + nb134_fjzO], xmm5
	addss   xmm0, [rsp + nb134_fixO]
	addss   xmm1, [rsp + nb134_fiyO]
	addss   xmm2, [rsp + nb134_fizO]
	movss  [rsp + nb134_fixO], xmm0
	movss  [rsp + nb134_fiyO], xmm1
	movss  [rsp + nb134_fizO], xmm2

	;# do  M coulomb interaction
	movaps xmm0, [rsp + nb134_rinvMM]
	movaps xmm7, xmm0	;# xmm7=rinv 
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [rsp + nb134_qqMH]
	movhps  xmm3, [rsp + nb134_qqMM]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	mulps  xmm7, xmm3 ;# vcoul=rinv*qq
	movaps xmm6, xmm7
	addps  xmm7, [rsp + nb134_vctot] 
	movaps [rsp + nb134_vctot], xmm7
	mulps  xmm0, xmm6

	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [rsp + nb134_dxMM]
	mulps   xmm1, [rsp + nb134_dyMM]
	mulps   xmm2, [rsp + nb134_dzMM]
	;# update forces M - j water 
	movaps  xmm3, [rsp + nb134_fjxO]
	movaps  xmm4, [rsp + nb134_fjyO]
	movaps  xmm5, [rsp + nb134_fjzO]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  [rsp + nb134_fjxO], xmm3
	movaps  [rsp + nb134_fjyO], xmm4
	movaps  [rsp + nb134_fjzO], xmm5
	addps   xmm0, [rsp + nb134_fixM]
	addps   xmm1, [rsp + nb134_fiyM]
	addps   xmm2, [rsp + nb134_fizM]
	movaps  [rsp + nb134_fixM], xmm0
	movaps  [rsp + nb134_fiyM], xmm1
	movaps  [rsp + nb134_fizM], xmm2	
	
	;# i H1 & H2 simultaneously first get i particle coords: 
    movaps  xmm0, [rsp + nb134_jxO]
    movaps  xmm1, [rsp + nb134_jyO]
    movaps  xmm2, [rsp + nb134_jzO]
    movaps  xmm3, xmm0
    movaps  xmm4, xmm1
    movaps  xmm5, xmm2
	subps   xmm0, [rsp + nb134_ixH1]
	subps   xmm1, [rsp + nb134_iyH1]
	subps   xmm2, [rsp + nb134_izH1]
	subps   xmm3, [rsp + nb134_ixH2] 
	subps   xmm4, [rsp + nb134_iyH2] 
	subps   xmm5, [rsp + nb134_izH2] 
	movaps [rsp + nb134_dxH1H1], xmm0
	movaps [rsp + nb134_dyH1H1], xmm1
	movaps [rsp + nb134_dzH1H1], xmm2
	movaps [rsp + nb134_dxH2H2], xmm3
	movaps [rsp + nb134_dyH2H2], xmm4
	movaps [rsp + nb134_dzH2H2], xmm5
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
	movaps  [rsp + nb134_rsqH1H1], xmm0
	movaps  [rsp + nb134_rsqH2H2], xmm4
		
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb134_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb134_half] ;# rinvH1H1
	mulps   xmm7, [rsp + nb134_half] ;# rinvH2H2
	movaps  [rsp + nb134_rinvH1H1], xmm3
	movaps  [rsp + nb134_rinvH2H2], xmm7
	
	;# Do H1 coulomb interaction
	movaps xmm0, xmm3
	movaps xmm7, xmm0	;# xmm7=rinv 
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [rsp + nb134_qqHH]
	movhps  xmm3, [rsp + nb134_qqMH]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	mulps  xmm7, xmm3 ;#vcoul 
	mulps  xmm0, xmm7
	
	addps  xmm7, [rsp + nb134_vctot] 
	movaps [rsp + nb134_vctot], xmm7

	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [rsp + nb134_dxH1H1]
	mulps   xmm1, [rsp + nb134_dyH1H1]
	mulps   xmm2, [rsp + nb134_dzH1H1]
	;# update forces H1 - j water 
	movaps  xmm3, [rsp + nb134_fjxO]
	movaps  xmm4, [rsp + nb134_fjyO]
	movaps  xmm5, [rsp + nb134_fjzO]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  [rsp + nb134_fjxO], xmm3
	movaps  [rsp + nb134_fjyO], xmm4
	movaps  [rsp + nb134_fjzO], xmm5
	addps   xmm0, [rsp + nb134_fixH1]
	addps   xmm1, [rsp + nb134_fiyH1]
	addps   xmm2, [rsp + nb134_fizH1]
	movaps  [rsp + nb134_fixH1], xmm0
	movaps  [rsp + nb134_fiyH1], xmm1
	movaps  [rsp + nb134_fizH1], xmm2	

	;# H2 Coulomb
	movaps xmm0, [rsp + nb134_rinvH2H2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [rsp + nb134_qqHH]
	movhps  xmm3, [rsp + nb134_qqMH]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	mulps  xmm7, xmm3
	mulps  xmm0, xmm7
	
	addps  xmm7, [rsp + nb134_vctot] ;# local vctot summation variable
	movaps [rsp + nb134_vctot], xmm7

	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [rsp + nb134_dxH2H2]
	mulps   xmm1, [rsp + nb134_dyH2H2]
	mulps   xmm2, [rsp + nb134_dzH2H2]
	;# update forces H2 - j water 
	movaps  xmm3, [rsp + nb134_fjxO]
	movaps  xmm4, [rsp + nb134_fjyO]
	movaps  xmm5, [rsp + nb134_fjzO]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  [rsp + nb134_fjxO], xmm3
	movaps  [rsp + nb134_fjyO], xmm4
	movaps  [rsp + nb134_fjzO], xmm5
	addps   xmm0, [rsp + nb134_fixH2]
	addps   xmm1, [rsp + nb134_fiyH2]
	addps   xmm2, [rsp + nb134_fizH2]
	movaps  [rsp + nb134_fixH2], xmm0
	movaps  [rsp + nb134_fiyH2], xmm1
	movaps  [rsp + nb134_fizH2], xmm2			
	
	mov     rsi, [rbp + nb134_faction]
	;# update j water forces from local variables.
	;# transpose back first
	movaps  xmm0, [rsp + nb134_fjxO] ;# Ox H1x H2x Mx 
	movaps  xmm1, [rsp + nb134_fjyO] ;# Oy H1y H2y My
	movaps  xmm2, [rsp + nb134_fjzO] ;# Oz H1z H2z Mz

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
	
	dec dword ptr [rsp + nb134_innerk]
	jz    .nb134_updateouterdata
	jmp   .nb134_single_loop
.nb134_updateouterdata:
	mov   ecx, [rsp + nb134_ii3]
	mov   rdi, [rbp + nb134_faction]
	mov   rsi, [rbp + nb134_fshift]
	mov   edx, [rsp + nb134_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb134_fixO]
	movaps xmm1, [rsp + nb134_fiyO] 
	movaps xmm2, [rsp + nb134_fizO]

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
	shufps  xmm6, xmm6, 8 ;# constant 00001000	

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb134_fixH1]
	movaps xmm1, [rsp + nb134_fiyH1]
	movaps xmm2, [rsp + nb134_fizH1]

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
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb134_fixH2]
	movaps xmm1, [rsp + nb134_fiyH2]
	movaps xmm2, [rsp + nb134_fizH2]

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
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# accumulate Mi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb134_fixM]
	movaps xmm1, [rsp + nb134_fiyM]
	movaps xmm2, [rsp + nb134_fizM]

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
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# increment fshift force  
	movlps  xmm3, [rsi + rdx*4]
	movss  xmm4, [rsi + rdx*4 + 8]
	subps  xmm3, xmm6
	subss  xmm4, xmm7
	movlps  [rsi + rdx*4],    xmm3
	movss  [rsi + rdx*4 + 8], xmm4

	;# get n from stack
	mov esi, [rsp + nb134_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb134_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb134_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb134_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb134_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb134_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb134_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb134_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb134_n], esi
        jmp .nb134_outer
.nb134_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb134_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb134_end
        ;# non-zero, do one more workunit
        jmp   .nb134_threadloop
.nb134_end:
	mov eax, [rsp + nb134_nouter]
	mov ebx, [rsp + nb134_ninner]
	mov rcx, [rbp + nb134_outeriter]
	mov rdx, [rbp + nb134_inneriter]
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

	
.globl nb_kernel134nf_x86_64_sse
.globl _nb_kernel134nf_x86_64_sse
nb_kernel134nf_x86_64_sse:	
_nb_kernel134nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb134nf_fshift,         16
.equiv          nb134nf_gid,            24
.equiv          nb134nf_pos,            32
.equiv          nb134nf_faction,        40
.equiv          nb134nf_charge,         48
.equiv          nb134nf_p_facel,        56
.equiv          nb134nf_argkrf,         64
.equiv          nb134nf_argcrf,         72
.equiv          nb134nf_Vc,             80
.equiv          nb134nf_type,           88
.equiv          nb134nf_p_ntype,        96
.equiv          nb134nf_vdwparam,       104
.equiv          nb134nf_Vvdw,           112
.equiv          nb134nf_p_tabscale,     120
.equiv          nb134nf_VFtab,          128
.equiv          nb134nf_invsqrta,       136
.equiv          nb134nf_dvda,           144
.equiv          nb134nf_p_gbtabscale,   152
.equiv          nb134nf_GBtab,          160
.equiv          nb134nf_p_nthreads,     168
.equiv          nb134nf_count,          176
.equiv          nb134nf_mtx,            184
.equiv          nb134nf_outeriter,      192
.equiv          nb134nf_inneriter,      200
.equiv          nb134nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb134nf_ixO,            0
.equiv          nb134nf_iyO,            16
.equiv          nb134nf_izO,            32
.equiv          nb134nf_ixH1,           48
.equiv          nb134nf_iyH1,           64
.equiv          nb134nf_izH1,           80
.equiv          nb134nf_ixH2,           96
.equiv          nb134nf_iyH2,           112
.equiv          nb134nf_izH2,           128
.equiv          nb134nf_ixM,            144
.equiv          nb134nf_iyM,            160
.equiv          nb134nf_izM,            176
.equiv          nb134nf_jxO,            192
.equiv          nb134nf_jyO,            208
.equiv          nb134nf_jzO,            224
.equiv          nb134nf_jxH1,           240
.equiv          nb134nf_jyH1,           256
.equiv          nb134nf_jzH1,           272
.equiv          nb134nf_jxH2,           288
.equiv          nb134nf_jyH2,           304
.equiv          nb134nf_jzH2,           320
.equiv          nb134nf_jxM,            336
.equiv          nb134nf_jyM,            352
.equiv          nb134nf_jzM,            368
.equiv          nb134nf_qqMM,           384
.equiv          nb134nf_qqMH,           400
.equiv          nb134nf_qqHH,           416
.equiv          nb134nf_tsc,            432
.equiv          nb134nf_c6,             448
.equiv          nb134nf_c12,            464
.equiv          nb134nf_vctot,          480
.equiv          nb134nf_Vvdwtot,        496
.equiv          nb134nf_half,           512
.equiv          nb134nf_three,          528
.equiv          nb134nf_rsqOO,          544
.equiv          nb134nf_rsqH1H1,        560
.equiv          nb134nf_rsqH1H2,        576
.equiv          nb134nf_rsqH1M,         592
.equiv          nb134nf_rsqH2H1,        608
.equiv          nb134nf_rsqH2H2,        624
.equiv          nb134nf_rsqH2M,         640
.equiv          nb134nf_rsqMH1,         656
.equiv          nb134nf_rsqMH2,         672
.equiv          nb134nf_rsqMM,          688
.equiv          nb134nf_rinvOO,         704
.equiv          nb134nf_rinvH1H1,       720
.equiv          nb134nf_rinvH1H2,       736
.equiv          nb134nf_rinvH1M,        752
.equiv          nb134nf_rinvH2H1,       768
.equiv          nb134nf_rinvH2H2,       784
.equiv          nb134nf_rinvH2M,        800
.equiv          nb134nf_rinvMH1,        816
.equiv          nb134nf_rinvMH2,        832
.equiv          nb134nf_rinvMM,         848
.equiv          nb134nf_krf,            864
.equiv          nb134nf_crf,            880
.equiv          nb134nf_is3,            896
.equiv          nb134nf_ii3,            900
.equiv          nb134nf_innerjjnr,      904
.equiv          nb134nf_nri,            912
.equiv          nb134nf_iinr,           920
.equiv          nb134nf_jindex,         928
.equiv          nb134nf_jjnr,           936
.equiv          nb134nf_shift,          944
.equiv          nb134nf_shiftvec,       952
.equiv          nb134nf_facel,          960
.equiv          nb134nf_innerk,         968
.equiv          nb134nf_n,              972
.equiv          nb134nf_nn1,            976
.equiv          nb134nf_nouter,         980
.equiv          nb134nf_ninner,         984

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
	mov [rsp + nb134nf_nouter], eax
	mov [rsp + nb134nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb134nf_nri], edi
	mov [rsp + nb134nf_iinr], rsi
	mov [rsp + nb134nf_jindex], rdx
	mov [rsp + nb134nf_jjnr], rcx
	mov [rsp + nb134nf_shift], r8
	mov [rsp + nb134nf_shiftvec], r9
	mov rsi, [rbp + nb134nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb134nf_facel], xmm0

	mov rax, [rbp + nb134nf_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb134nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb134nf_half], eax
	movss xmm1, [rsp + nb134nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb134nf_half],  xmm1
	movaps [rsp + nb134nf_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb134nf_iinr]   ;# rcx = pointer into iinr[]
	mov   ebx, [rcx]		;# ebx =ii 

	mov   rdx, [rbp + nb134nf_charge]
	movss xmm5, [rdx + rbx*4 + 4]	
	movss xmm3, [rdx + rbx*4 + 12]	
	movss xmm4, xmm3	
	mov rsi, [rbp + nb134nf_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb134nf_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb134nf_qqMM], xmm3
	movaps [rsp + nb134nf_qqMH], xmm4
	movaps [rsp + nb134nf_qqHH], xmm5
	
	xorps xmm0, xmm0
	mov   rdx, [rbp + nb134nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb134nf_p_ntype]
	imul  ecx, [rdi]  ;# rcx = ntia = 2*ntype*type[ii0] 
	add   rdx, rcx
	mov   rax, [rbp + nb134nf_vdwparam]
	movlps xmm0, [rax + rdx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0x55
	movaps [rsp + nb134nf_c6], xmm0
	movaps [rsp + nb134nf_c12], xmm1

.nb134nf_threadloop:
        mov   rsi, [rbp + nb134nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb134nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb134nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb134nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb134nf_n], eax
        mov [rsp + nb134nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb134nf_outerstart
        jmp .nb134nf_end
	
.nb134nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb134nf_nouter]
	mov [rsp + nb134nf_nouter], ebx

.nb134nf_outer:
	mov   rax, [rsp + nb134nf_shift]  	;# eax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb134nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb134nf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb134nf_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb134nf_pos]	;# eax = base of pos[]  
	mov   [rsp + nb134nf_ii3], ebx	
	
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
	movaps [rsp + nb134nf_ixO], xmm3
	movaps [rsp + nb134nf_iyO], xmm4
	movaps [rsp + nb134nf_izO], xmm5
	movaps [rsp + nb134nf_ixH1], xmm6
	movaps [rsp + nb134nf_iyH1], xmm7

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
	movaps [rsp + nb134nf_izH1], xmm6
	movaps [rsp + nb134nf_ixH2], xmm0
	movaps [rsp + nb134nf_iyH2], xmm1
	movaps [rsp + nb134nf_izH2], xmm2
	movaps [rsp + nb134nf_ixM], xmm3
	movaps [rsp + nb134nf_iyM], xmm4
	movaps [rsp + nb134nf_izM], xmm5

	;# clear vctot 
	xorps xmm4, xmm4
	movaps [rsp + nb134nf_vctot], xmm4
	movaps [rsp + nb134nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb134nf_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb134nf_pos]
	mov   rax, [rsp + nb134nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb134nf_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb134nf_ninner]
	mov   [rsp + nb134nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb134nf_innerk], edx	;# number of innerloop atoms 
	jge   .nb134nf_unroll_loop
	jmp   .nb134nf_single_check
.nb134nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb134nf_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb134nf_innerjjnr],  16 ;# advance pointer (unroll 4) 

	mov rsi, [rbp + nb134nf_pos]   	;# base of pos[] 

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
	movaps [rsp + nb134nf_jxO], xmm0
	movaps [rsp + nb134nf_jyO], xmm1
	movaps [rsp + nb134nf_jzO], xmm2
	movaps [rsp + nb134nf_jxH1], xmm3

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
	movaps [rsp + nb134nf_jyH1], xmm0
	movaps [rsp + nb134nf_jzH1], xmm1
	movaps [rsp + nb134nf_jxH2], xmm2
	movaps [rsp + nb134nf_jyH2], xmm3

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
	movaps [rsp + nb134nf_jzH2], xmm0
	movaps [rsp + nb134nf_jxM], xmm1
	movaps [rsp + nb134nf_jyM], xmm2
	movaps [rsp + nb134nf_jzM], xmm3
	
	;# start calculating pairwise distances
	movaps xmm0, [rsp + nb134nf_ixO]
	movaps xmm1, [rsp + nb134nf_iyO]
	movaps xmm2, [rsp + nb134nf_izO]
	movaps xmm3, [rsp + nb134nf_ixH1]
	movaps xmm4, [rsp + nb134nf_iyH1]
	movaps xmm5, [rsp + nb134nf_izH1]
	subps  xmm0, [rsp + nb134nf_jxO]
	subps  xmm1, [rsp + nb134nf_jyO]
	subps  xmm2, [rsp + nb134nf_jzO]
	subps  xmm3, [rsp + nb134nf_jxH1]
	subps  xmm4, [rsp + nb134nf_jyH1]
	subps  xmm5, [rsp + nb134nf_jzH1]
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
	movaps [rsp + nb134nf_rsqOO], xmm0
	movaps [rsp + nb134nf_rsqH1H1], xmm3

	movaps xmm0, [rsp + nb134nf_ixH1]
	movaps xmm1, [rsp + nb134nf_iyH1]
	movaps xmm2, [rsp + nb134nf_izH1]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [rsp + nb134nf_jxH2]
	subps  xmm1, [rsp + nb134nf_jyH2]
	subps  xmm2, [rsp + nb134nf_jzH2]
	subps  xmm3, [rsp + nb134nf_jxM]
	subps  xmm4, [rsp + nb134nf_jyM]
	subps  xmm5, [rsp + nb134nf_jzM]
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
	movaps [rsp + nb134nf_rsqH1H2], xmm0
	movaps [rsp + nb134nf_rsqH1M], xmm3

	movaps xmm0, [rsp + nb134nf_ixH2]
	movaps xmm1, [rsp + nb134nf_iyH2]
	movaps xmm2, [rsp + nb134nf_izH2]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [rsp + nb134nf_jxH1]
	subps  xmm1, [rsp + nb134nf_jyH1]
	subps  xmm2, [rsp + nb134nf_jzH1]
	subps  xmm3, [rsp + nb134nf_jxH2]
	subps  xmm4, [rsp + nb134nf_jyH2]
	subps  xmm5, [rsp + nb134nf_jzH2]
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
	movaps [rsp + nb134nf_rsqH2H1], xmm0
	movaps [rsp + nb134nf_rsqH2H2], xmm3

	movaps xmm0, [rsp + nb134nf_ixH2]
	movaps xmm1, [rsp + nb134nf_iyH2]
	movaps xmm2, [rsp + nb134nf_izH2]
	movaps xmm3, [rsp + nb134nf_ixM]
	movaps xmm4, [rsp + nb134nf_iyM]
	movaps xmm5, [rsp + nb134nf_izM]
	subps  xmm0, [rsp + nb134nf_jxM]
	subps  xmm1, [rsp + nb134nf_jyM]
	subps  xmm2, [rsp + nb134nf_jzM]
	subps  xmm3, [rsp + nb134nf_jxH1]
	subps  xmm4, [rsp + nb134nf_jyH1]
	subps  xmm5, [rsp + nb134nf_jzH1]
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
	movaps [rsp + nb134nf_rsqH2M], xmm0
	movaps [rsp + nb134nf_rsqMH1], xmm4

	movaps xmm0, [rsp + nb134nf_ixM]
	movaps xmm1, [rsp + nb134nf_iyM]
	movaps xmm2, [rsp + nb134nf_izM]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [rsp + nb134nf_jxH2]
	subps  xmm1, [rsp + nb134nf_jyH2]
	subps  xmm2, [rsp + nb134nf_jzH2]
	subps  xmm3, [rsp + nb134nf_jxM]
	subps  xmm4, [rsp + nb134nf_jyM]
	subps  xmm5, [rsp + nb134nf_jzM]
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
	movaps [rsp + nb134nf_rsqMH2], xmm0
	movaps [rsp + nb134nf_rsqMM], xmm4

	;# start by doing invsqrt for OO
	rsqrtps xmm1, [rsp + nb134nf_rsqOO]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb134nf_three]
	mulps   xmm1, [rsp + nb134nf_rsqOO]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + nb134nf_half] 
	movaps  [rsp + nb134nf_rinvOO], xmm3
	
	;# more invsqrt ops - do two at a time.
	rsqrtps xmm1, [rsp + nb134nf_rsqH1H1]
	rsqrtps xmm5, [rsp + nb134nf_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb134nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb134nf_rsqH1H1]
	mulps   xmm5, [rsp + nb134nf_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb134nf_half] ;# rinvH1H1 
	mulps   xmm7, [rsp + nb134nf_half] ;# rinvH1H2 
	movaps  [rsp + nb134nf_rinvH1H1], xmm3
	movaps  [rsp + nb134nf_rinvH1H2], xmm7
			
	rsqrtps xmm1, [rsp + nb134nf_rsqH1M]
	rsqrtps xmm5, [rsp + nb134nf_rsqH2H1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb134nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb134nf_rsqH1M]
	mulps   xmm5, [rsp + nb134nf_rsqH2H1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb134nf_half] 
	mulps   xmm7, [rsp + nb134nf_half]
	movaps  [rsp + nb134nf_rinvH1M], xmm3
	movaps  [rsp + nb134nf_rinvH2H1], xmm7
			
	rsqrtps xmm1, [rsp + nb134nf_rsqH2H2]
	rsqrtps xmm5, [rsp + nb134nf_rsqH2M]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb134nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb134nf_rsqH2H2]
	mulps   xmm5, [rsp + nb134nf_rsqH2M]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb134nf_half] 
	mulps   xmm7, [rsp + nb134nf_half]
	movaps  [rsp + nb134nf_rinvH2H2], xmm3
	movaps  [rsp + nb134nf_rinvH2M], xmm7
	
	rsqrtps xmm1, [rsp + nb134nf_rsqMH1]
	rsqrtps xmm5, [rsp + nb134nf_rsqMH2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb134nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb134nf_rsqMH1]
	mulps   xmm5, [rsp + nb134nf_rsqMH2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb134nf_half] 
	mulps   xmm7, [rsp + nb134nf_half]
	movaps  [rsp + nb134nf_rinvMH1], xmm3
	movaps  [rsp + nb134nf_rinvMH2], xmm7
        		
	rsqrtps xmm1, [rsp + nb134nf_rsqMM]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb134nf_three]
	mulps   xmm1, [rsp + nb134nf_rsqMM]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + nb134nf_half] 
	movaps  [rsp + nb134nf_rinvMM], xmm3
	
	;# start with OO LJ interaction
	movaps xmm0, [rsp + nb134nf_rinvOO]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb134nf_rsqOO] ;# xmm1=r 
	mulps  xmm1, [rsp + nb134nf_tsc]
		
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

    mov  rsi, [rbp + nb134nf_VFtab]
	
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

    movaps xmm4, [rsp + nb134nf_c6]
    mulps  xmm5, xmm4    ;# Vvdw6 

    addps  xmm5, [rsp + nb134nf_Vvdwtot]
    movaps [rsp + nb134nf_Vvdwtot], xmm5

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
 
    movaps xmm4, [rsp + nb134nf_c12]
    mulps  xmm5, xmm4 ;# Vvdw12 

    addps  xmm5, [rsp + nb134nf_Vvdwtot]
    movaps [rsp + nb134nf_Vvdwtot], xmm5

	;# Coulomb interactions 
	movaps xmm0, [rsp + nb134nf_rinvH1H1]
	movaps xmm1, [rsp + nb134nf_rinvH1M]
	movaps xmm2, [rsp + nb134nf_rinvMM]
	addps  xmm0, [rsp + nb134nf_rinvH1H2]
	addps  xmm1, [rsp + nb134nf_rinvH2M]
	addps  xmm0, [rsp + nb134nf_rinvH2H1]
	addps  xmm1, [rsp + nb134nf_rinvMH1]
	addps  xmm0, [rsp + nb134nf_rinvH2H2]
	addps  xmm1, [rsp + nb134nf_rinvMH2]
	
	mulps  xmm0, [rsp + nb134nf_qqHH]
	mulps  xmm1, [rsp + nb134nf_qqMH]
	mulps  xmm2, [rsp + nb134nf_qqMM]
	
	addps  xmm0, xmm1
	addps  xmm2, [rsp + nb134nf_vctot]
	addps  xmm0, xmm2
	movaps [rsp + nb134nf_vctot], xmm0

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb134nf_innerk],  4
	jl    .nb134nf_single_check
	jmp   .nb134nf_unroll_loop
.nb134nf_single_check:
	add dword ptr [rsp + nb134nf_innerk],  4
	jnz   .nb134nf_single_loop
	jmp   .nb134nf_updateouterdata
.nb134nf_single_loop:
	mov   rdx, [rsp + nb134nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb134nf_innerjjnr],  4

	mov rsi, [rbp + nb134nf_pos]
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
	movaps [rsp + nb134nf_jxO], xmm6
	movaps [rsp + nb134nf_jyO], xmm3
	movaps [rsp + nb134nf_jzO], xmm1

	;# do O and M in parallel
	movaps xmm0, [rsp + nb134nf_ixO]
	movaps xmm1, [rsp + nb134nf_iyO]
	movaps xmm2, [rsp + nb134nf_izO]
	movaps xmm3, [rsp + nb134nf_ixM]
	movaps xmm4, [rsp + nb134nf_iyM]
	movaps xmm5, [rsp + nb134nf_izM]
	subps  xmm0, [rsp + nb134nf_jxO]
	subps  xmm1, [rsp + nb134nf_jyO]
	subps  xmm2, [rsp + nb134nf_jzO]
	subps  xmm3, [rsp + nb134nf_jxO]
	subps  xmm4, [rsp + nb134nf_jyO]
	subps  xmm5, [rsp + nb134nf_jzO]
	
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
	movaps [rsp + nb134nf_rsqOO], xmm0
	movaps [rsp + nb134nf_rsqMM], xmm4
	
	;# do 1/x for O
	rsqrtps xmm1, xmm0
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb134nf_three]
	mulps   xmm1, xmm0	;# rsq*lu*lu 
	subps   xmm3, xmm1	;# constant 30-rsq*lu*lu 
	mulps   xmm3, xmm2	;# lu*(3-rsq*lu*lu) 
	mulps   xmm3, [rsp + nb134nf_half]
	movaps  [rsp + nb134nf_rinvOO], xmm3	;# rinvH2 
	
	;# 1/sqrt(x) for M
	rsqrtps xmm5, xmm4
	movaps  xmm6, xmm5	
	mulps   xmm5, xmm5
	movaps  xmm7, [rsp + nb134nf_three]
	mulps   xmm5, xmm4
	subps   xmm7, xmm5
	mulps   xmm7, xmm6
	mulps   xmm7, [rsp + nb134nf_half] ;# rinv iH1 - j water 
	movaps [rsp + nb134nf_rinvMM], xmm7


	;# LJ table interaction
	movaps xmm0, [rsp + nb134nf_rinvOO]
	movss xmm1, xmm0
	mulss  xmm1, [rsp + nb134nf_rsqOO] ;# xmm1=r 
	mulss  xmm1, [rsp + nb134nf_tsc]
		
    cvttps2pi mm6, xmm1
    cvtpi2ps xmm3, mm6
	subss    xmm1, xmm3	;# xmm1=eps 
    movss xmm2, xmm1
    mulss  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 3
	
    mov  rsi, [rbp + nb134nf_VFtab]
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

    movss xmm4, [rsp + nb134nf_c6]
    mulss  xmm5, xmm4    ;# Vvdw6 

    addss  xmm5, [rsp + nb134nf_Vvdwtot]
    movss [rsp + nb134nf_Vvdwtot], xmm5

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

    movss xmm4, [rsp + nb134nf_c12]
    mulss  xmm5, xmm4 ;# Vvdw12 
    addss  xmm5, [rsp + nb134nf_Vvdwtot]
    movss [rsp + nb134nf_Vvdwtot], xmm5

	;# do  M coulomb interaction
	movaps xmm0, [rsp + nb134nf_rinvMM]

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [rsp + nb134nf_qqMH]
	movhps  xmm3, [rsp + nb134nf_qqMM]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	mulps  xmm0, xmm3
	addps  xmm0, [rsp + nb134nf_vctot] 
	movaps [rsp + nb134nf_vctot], xmm0
	
	;# i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [rsp + nb134nf_ixH1]
	movaps  xmm1, [rsp + nb134nf_iyH1]
	movaps  xmm2, [rsp + nb134nf_izH1]	
	movaps  xmm3, [rsp + nb134nf_ixH2] 
	movaps  xmm4, [rsp + nb134nf_iyH2] 
	movaps  xmm5, [rsp + nb134nf_izH2] 
	subps   xmm0, [rsp + nb134nf_jxO]
	subps   xmm1, [rsp + nb134nf_jyO]
	subps   xmm2, [rsp + nb134nf_jzO]
	subps   xmm3, [rsp + nb134nf_jxO]
	subps   xmm4, [rsp + nb134nf_jyO]
	subps   xmm5, [rsp + nb134nf_jzO]
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
	movaps  [rsp + nb134nf_rsqH1H1], xmm0
	movaps  [rsp + nb134nf_rsqH2H2], xmm4
		
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb134nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb134nf_half] ;# rinvH1H1
	mulps   xmm7, [rsp + nb134nf_half] ;# rinvH2H2
	movaps  [rsp + nb134nf_rinvH1H1], xmm3
	movaps  [rsp + nb134nf_rinvH2H2], xmm7
	

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [rsp + nb134nf_qqHH]
	movhps  xmm3, [rsp + nb134nf_qqMH]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	;# Do coulomb interactions	
	movaps xmm0, [rsp + nb134nf_rinvH1H1]
	addps  xmm0, [rsp + nb134nf_rinvH2H2]
	mulps  xmm0, xmm3
	addps  xmm0, [rsp + nb134nf_vctot] 
	movaps [rsp + nb134nf_vctot], xmm0
	
	dec dword ptr [rsp + nb134nf_innerk]
	jz    .nb134nf_updateouterdata
	jmp   .nb134nf_single_loop
.nb134nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb134nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb134nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb134nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb134nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb134nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb134nf_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb134nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb134nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb134nf_n], esi
        jmp .nb134nf_outer
.nb134nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb134nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb134nf_end
        ;# non-zero, do one more workunit
        jmp   .nb134nf_threadloop
.nb134nf_end:
	mov eax, [rsp + nb134nf_nouter]
	mov ebx, [rsp + nb134nf_ninner]
	mov rcx, [rbp + nb134nf_outeriter]
	mov rdx, [rbp + nb134nf_inneriter]
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

