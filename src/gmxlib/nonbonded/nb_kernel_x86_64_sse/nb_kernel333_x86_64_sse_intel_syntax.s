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


.globl nb_kernel333_x86_64_sse
.globl _nb_kernel333_x86_64_sse
nb_kernel333_x86_64_sse:	
_nb_kernel333_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb333_fshift,           16
.equiv          nb333_gid,              24
.equiv          nb333_pos,              32
.equiv          nb333_faction,          40
.equiv          nb333_charge,           48
.equiv          nb333_p_facel,          56
.equiv          nb333_argkrf,           64
.equiv          nb333_argcrf,           72
.equiv          nb333_Vc,               80
.equiv          nb333_type,             88
.equiv          nb333_p_ntype,          96
.equiv          nb333_vdwparam,         104
.equiv          nb333_Vvdw,             112
.equiv          nb333_p_tabscale,       120
.equiv          nb333_VFtab,            128
.equiv          nb333_invsqrta,         136
.equiv          nb333_dvda,             144
.equiv          nb333_p_gbtabscale,     152
.equiv          nb333_GBtab,            160
.equiv          nb333_p_nthreads,       168
.equiv          nb333_count,            176
.equiv          nb333_mtx,              184
.equiv          nb333_outeriter,        192
.equiv          nb333_inneriter,        200
.equiv          nb333_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb333_ixO,              0
.equiv          nb333_iyO,              16
.equiv          nb333_izO,              32
.equiv          nb333_ixH1,             48
.equiv          nb333_iyH1,             64
.equiv          nb333_izH1,             80
.equiv          nb333_ixH2,             96
.equiv          nb333_iyH2,             112
.equiv          nb333_izH2,             128
.equiv          nb333_ixM,              144
.equiv          nb333_iyM,              160
.equiv          nb333_izM,              176
.equiv          nb333_iqM,              192
.equiv          nb333_iqH,              208
.equiv          nb333_dxO,              224
.equiv          nb333_dyO,              240
.equiv          nb333_dzO,              256
.equiv          nb333_dxH1,             272
.equiv          nb333_dyH1,             288
.equiv          nb333_dzH1,             304
.equiv          nb333_dxH2,             320
.equiv          nb333_dyH2,             336
.equiv          nb333_dzH2,             352
.equiv          nb333_dxM,              368
.equiv          nb333_dyM,              384
.equiv          nb333_dzM,              400
.equiv          nb333_qqM,              416
.equiv          nb333_qqH,              432
.equiv          nb333_rinvO,            448
.equiv          nb333_rinvH1,           464
.equiv          nb333_rinvH2,           480
.equiv          nb333_rinvM,            496
.equiv          nb333_rO,               512
.equiv          nb333_rH1,              528
.equiv          nb333_rH2,              544
.equiv          nb333_rM,               560
.equiv          nb333_tsc,              576
.equiv          nb333_two,              592
.equiv          nb333_c6,               608
.equiv          nb333_c12,              624
.equiv          nb333_vctot,            640
.equiv          nb333_Vvdwtot,          656
.equiv          nb333_fixO,             672
.equiv          nb333_fiyO,             688
.equiv          nb333_fizO,             704
.equiv          nb333_fixH1,            720
.equiv          nb333_fiyH1,            736
.equiv          nb333_fizH1,            752
.equiv          nb333_fixH2,            768
.equiv          nb333_fiyH2,            784
.equiv          nb333_fizH2,            800
.equiv          nb333_fixM,             816
.equiv          nb333_fiyM,             832
.equiv          nb333_fizM,             848
.equiv          nb333_fjx,              864
.equiv          nb333_fjy,              880
.equiv          nb333_fjz,              896
.equiv          nb333_epsH1,            912
.equiv          nb333_epsH2,            928
.equiv          nb333_epsM,             944
.equiv          nb333_half,             960
.equiv          nb333_three,            976
.equiv          nb333_is3,              992
.equiv          nb333_ii3,              996
.equiv          nb333_nri,              1000
.equiv          nb333_iinr,             1008
.equiv          nb333_jindex,           1016
.equiv          nb333_jjnr,             1024
.equiv          nb333_shift,            1032
.equiv          nb333_shiftvec,         1040
.equiv          nb333_facel,            1048
.equiv          nb333_innerjjnr,        1056
.equiv          nb333_ntia,             1064
.equiv          nb333_innerk,           1068
.equiv          nb333_n,                1072
.equiv          nb333_nn1,              1076
.equiv          nb333_nouter,           1080
.equiv          nb333_ninner,           1084
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

	sub rsp, 1104		;# local variable stack space (n*16+8)
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
	mov [rsp + nb333_nouter], eax
	mov [rsp + nb333_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb333_nri], edi
	mov [rsp + nb333_iinr], rsi
	mov [rsp + nb333_jindex], rdx
	mov [rsp + nb333_jjnr], rcx
	mov [rsp + nb333_shift], r8
	mov [rsp + nb333_shiftvec], r9
	mov rsi, [rbp + nb333_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb333_facel], xmm0

	mov rax, [rbp + nb333_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb333_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb333_half], eax
	movss xmm1, [rsp + nb333_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb333_half],  xmm1
	movaps [rsp + nb333_two],  xmm2
	movaps [rsp + nb333_three],  xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb333_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]		;# ebx =ii 

	mov   rdx, [rbp + nb333_charge]
	movss xmm4, [rdx + rbx*4 + 4]	
	movss xmm3, [rdx + rbx*4 + 12]	
	mov rsi, [rbp + nb333_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb333_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb333_iqM], xmm3
	movaps [rsp + nb333_iqH], xmm4
	
	mov   rdx, [rbp + nb333_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb333_p_ntype]
	imul  ecx, [rdi]  	;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb333_ntia], ecx		

.nb333_threadloop:
        mov   rsi, [rbp + nb333_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb333_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb333_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb333_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb333_n], eax
        mov [rsp + nb333_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb333_outerstart
        jmp .nb333_end
	
.nb333_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb333_nouter]
	mov [rsp + nb333_nouter], ebx

.nb333_outer:
	mov   rax, [rsp + nb333_shift]  	;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb333_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb333_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb333_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2	
	movaps xmm6, xmm0
	movaps xmm7, xmm1
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb333_pos]	;# rax = base of pos[]  
	mov   [rsp + nb333_ii3], ebx

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
	movaps [rsp + nb333_ixO], xmm3
	movaps [rsp + nb333_iyO], xmm4
	movaps [rsp + nb333_izO], xmm5
	movaps [rsp + nb333_ixH1], xmm6
	movaps [rsp + nb333_iyH1], xmm7

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
	movaps [rsp + nb333_izH1], xmm6
	movaps [rsp + nb333_ixH2], xmm0
	movaps [rsp + nb333_iyH2], xmm1
	movaps [rsp + nb333_izH2], xmm2
	movaps [rsp + nb333_ixM], xmm3
	movaps [rsp + nb333_iyM], xmm4
	movaps [rsp + nb333_izM], xmm5
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb333_vctot], xmm4
	movaps [rsp + nb333_Vvdwtot], xmm4
	movaps [rsp + nb333_fixO], xmm4
	movaps [rsp + nb333_fiyO], xmm4
	movaps [rsp + nb333_fizO], xmm4
	movaps [rsp + nb333_fixH1], xmm4
	movaps [rsp + nb333_fiyH1], xmm4
	movaps [rsp + nb333_fizH1], xmm4
	movaps [rsp + nb333_fixH2], xmm4
	movaps [rsp + nb333_fiyH2], xmm4
	movaps [rsp + nb333_fizH2], xmm4
	movaps [rsp + nb333_fixM], xmm4
	movaps [rsp + nb333_fiyM], xmm4
	movaps [rsp + nb333_fizM], xmm4
	
	mov   rax, [rsp + nb333_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb333_pos]
	mov   rdi, [rbp + nb333_faction]	
	mov   rax, [rsp + nb333_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb333_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb333_ninner]
	mov   [rsp + nb333_ninner], ecx
	add   edx, 0
	mov   [rsp + nb333_innerk], edx	;# number of innerloop atoms 
	jge   .nb333_unroll_loop
	jmp   .nb333_odd_inner
.nb333_unroll_loop:
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb333_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4
	
	add qword ptr [rsp + nb333_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb333_charge]	;# base of charge[] 

	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]
	
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	 	;# and in xmm4 
	mulps  xmm3, [rsp + nb333_iqM]
	mulps  xmm4, [rsp + nb333_iqH]

	movaps  [rsp + nb333_qqM], xmm3
	movaps  [rsp + nb333_qqH], xmm4
	
	mov rsi, [rbp + nb333_type]
	mov r8d, [rsi + rax*4]
	mov r9d, [rsi + rbx*4]
	mov r10d, [rsi + rcx*4]
	mov r11d, [rsi + rdx*4]
	mov rsi, [rbp + nb333_vdwparam]
	shl r8d, 1	
	shl r9d, 1	
	shl r10d, 1	
	shl r11d, 1	
	mov edi, [rsp + nb333_ntia]
	add r8d, edi
	add r9d, edi
	add r10d, edi
	add r11d, edi

	movlps xmm6, [rsi + r8*4]
	movlps xmm7, [rsi + r10*4]
	movhps xmm6, [rsi + r9*4]
	movhps xmm7, [rsi + r11*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm6, xmm7, 221  ;# 11011101

	movaps [rsp + nb333_c6], xmm4
	movaps [rsp + nb333_c12], xmm6

	mov rsi, [rbp + nb333_pos]   	;# base of pos[] 

	lea   rax, [rax + rax*2] 	;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2] 	;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	

	;# move four coordinates to xmm0-xmm2 	
	movlps xmm4, [rsi + rax*4]
	movlps xmm5, [rsi + rcx*4]
	movss xmm2, [rsi + rax*4 + 8]
	movss xmm6, [rsi + rcx*4 + 8]

	movhps xmm4, [rsi + rbx*4]
	movhps xmm5, [rsi + rdx*4]

	movss xmm0, [rsi + rbx*4 + 8]
	movss xmm1, [rsi + rdx*4 + 8]

	shufps xmm2, xmm0, 0
	shufps xmm6, xmm1, 0
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 136  ;# 10001000
	
	shufps xmm0, xmm5, 136  ;# 10001000
	shufps xmm1, xmm5, 221  ;# 11011101		

    ;# xmm0 = jx
    ;# xmm1 = jy
    ;# xmm2 = jz
    
    ;# O interaction
    ;# copy to xmm3-xmm5
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    
    subps xmm3, [rsp + nb333_ixO]
    subps xmm4, [rsp + nb333_iyO]
    subps xmm5, [rsp + nb333_izO]
    
    movaps [rsp + nb333_dxO], xmm3
    movaps [rsp + nb333_dyO], xmm4
    movaps [rsp + nb333_dzO], xmm5
    
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm3, xmm4
	addps  xmm3, xmm5
    ;# xmm3=rsq
    
    ;# calculate rinv=1/sqrt(rsq)
	rsqrtps xmm5, xmm3
	movaps xmm15, xmm5
	mulps xmm5, xmm5
	movaps xmm4, [rsp + nb333_three]
	mulps xmm5, xmm3	;# rsq*lu*lu 	
    subps xmm4, xmm5	;# 30-rsq*lu*lu 
	mulps xmm4, xmm15	
	mulps xmm4, [rsp + nb333_half]	
	movaps xmm15, xmm4
	mulps  xmm3, xmm4	
    ;# xmm15=rinv
    ;# xmm3=r

    mulps xmm3, [rsp + nb333_tsc] ;# rtab

    ;# truncate and convert to integers
    cvttps2dq xmm5, xmm3
    
    ;# convert back to float
    cvtdq2ps  xmm4, xmm5
    
    ;# multiply by 4
    pslld   xmm5, 2

    ;# multiply by three (copy, mult. by two, add back)
    movaps  xmm10, xmm5
    pslld   xmm10, 1
    paddd   xmm5, xmm10    

    ;# calculate eps
    subps     xmm3, xmm4    ;# xmm3=eps
    
    ;# move to integer registers
    movhlps xmm6, xmm5
    movd    r8d, xmm5
    movd    r10d, xmm6
    pshufd  xmm5, xmm5, 1
    pshufd  xmm6, xmm6, 1
    movd    r9d, xmm5
    movd    r11d, xmm6
    ;# xmm3=eps
    ;# xmm15=rinv

	mov rsi, [rbp + nb333_VFtab]
    ;# calculate LJ table
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
    
	movlps xmm13, [rsi + r10*4 + 24]
	movlps xmm14, [rsi + r10*4 + 40]

	movhps xmm7,  [rsi + r9*4 + 24]
	movhps xmm11, [rsi + r9*4 + 40]
    
	movhps xmm13, [rsi + r11*4 + 24]
	movhps xmm14, [rsi + r11*4 + 40]

    movaps xmm6, xmm7
    movaps xmm10, xmm11
    
	shufps xmm6, xmm13, 136  ;# 10001000
	shufps xmm10, xmm14, 136  ;# 10001000
	shufps xmm7, xmm13, 221  ;# 11011101
	shufps xmm11, xmm14, 221  ;# 11011101
    ;# dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11
    
    mulps  xmm7, xmm3    ;# Heps
    mulps  xmm11, xmm3 
    mulps  xmm6, xmm3   ;# Geps
    mulps  xmm10, xmm3 
    mulps  xmm7, xmm3   ;# Heps2
    mulps  xmm11, xmm3 
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
    mulps  xmm5, xmm3  ;# eps*Fp
    mulps  xmm9, xmm3
    movaps xmm12, [rsp + nb333_c6]
    movaps xmm13, [rsp + nb333_c12]
    addps  xmm5, xmm4 ;# VV
    addps  xmm9, xmm8

    mulps  xmm5, xmm12  ;# VV*c6 = vnb6
    mulps  xmm9, xmm13  ;# VV*c12 = vnb12
    addps  xmm5, xmm9
    addps  xmm5, [rsp + nb333_Vvdwtot]
    movaps [rsp + nb333_Vvdwtot], xmm5
        
    mulps  xmm7, xmm12   ;# FF*c6 = fnb6
    mulps  xmm11, xmm13   ;# FF*c12  = fnb12
    addps  xmm7, xmm11
    
    mulps  xmm7, [rsp + nb333_tsc]
    mulps  xmm7, xmm15   ;# -fscal
    xorps  xmm9, xmm9
    
    subps  xmm9, xmm7     ;# fscal
    movaps xmm10, xmm9
    movaps xmm11, xmm9

    mulps  xmm9,  [rsp + nb333_dxO] ;# fx/fy/fz
    mulps  xmm10, [rsp + nb333_dyO]
    mulps  xmm11, [rsp + nb333_dzO]

    ;# save j force temporarily
    movaps [rsp + nb333_fjx], xmm9
    movaps [rsp + nb333_fjy], xmm10
    movaps [rsp + nb333_fjz], xmm11
    
    ;# increment i O force
    addps xmm9, [rsp + nb333_fixO]
    addps xmm10, [rsp + nb333_fiyO]
    addps xmm11, [rsp + nb333_fizO]
    movaps [rsp + nb333_fixO], xmm9
    movaps [rsp + nb333_fiyO], xmm10
    movaps [rsp + nb333_fizO], xmm11
    ;# finished O LJ interaction.


    ;# do H1, H2, and M interactions in parallel.
    ;# xmm0-xmm2 still contain j coordinates.                
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + nb333_ixH1]
    subps xmm1, [rsp + nb333_iyH1]
    subps xmm2, [rsp + nb333_izH1]
    subps xmm3, [rsp + nb333_ixH2]
    subps xmm4, [rsp + nb333_iyH2]
    subps xmm5, [rsp + nb333_izH2]
    subps xmm6, [rsp + nb333_ixM]
    subps xmm7, [rsp + nb333_iyM]
    subps xmm8, [rsp + nb333_izM]    

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps [rsp + nb333_dxH1], xmm0
	movaps [rsp + nb333_dyH1], xmm1
	movaps [rsp + nb333_dzH1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb333_dxH2], xmm3
	movaps [rsp + nb333_dyH2], xmm4
	movaps [rsp + nb333_dzH2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb333_dxM], xmm6
	movaps [rsp + nb333_dyM], xmm7
	movaps [rsp + nb333_dzM], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt for j atoms
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + nb333_three]
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

	movaps  xmm4, [rsp + nb333_half]
	mulps   xmm9, xmm4  ;# rinvH1 
	mulps   xmm10, xmm4 ;# rinvH2
    mulps   xmm11, xmm4 ;# rinvM

	movaps  [rsp + nb333_rinvH1], xmm9
	movaps  [rsp + nb333_rinvH2], xmm10
	movaps  [rsp + nb333_rinvM], xmm11
	
	;# interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movaps xmm1, [rsp + nb333_tsc]
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
        
    mov  rsi, [rbp + nb333_VFtab]

    ;# calculate eps
    subps     xmm0, xmm2
    subps     xmm3, xmm5
    subps     xmm6, xmm8

    movaps    [rsp + nb333_epsH1], xmm0
    movaps    [rsp + nb333_epsH2], xmm3
    movaps    [rsp + nb333_epsM], xmm6

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
    
    movaps xmm12, [rsp + nb333_epsH1]
    movaps xmm13, [rsp + nb333_epsH2]
    movaps xmm14, [rsp + nb333_epsM]
    
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
    movaps xmm12, [rsp + nb333_qqH]
    movaps xmm13, [rsp + nb333_qqM]
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
    addps  xmm1, [rsp + nb333_vctot]
    addps  xmm5, xmm9
    addps  xmm1, xmm5
    movaps [rsp + nb333_vctot], xmm1
    
    movaps xmm10, [rsp + nb333_tsc]
    mulps  xmm3, xmm10  ;# fscal
    mulps  xmm7, xmm10
    mulps  xmm10, xmm11
    
    movd eax, mm0
    movd ebx, mm1
    movd ecx, mm2
    movd edx, mm3

	;# move j forces to local temp variables 
    mov rdi, [rbp + nb333_faction]
    movlps xmm11, [rdi + rax*4] ;# jxa jya  -   -
    movlps xmm12, [rdi + rcx*4] ;# jxc jyc  -   -
    movhps xmm11, [rdi + rbx*4] ;# jxa jya jxb jyb 
    movhps xmm12, [rdi + rdx*4] ;# jxc jyc jxd jyd 

    movss  xmm13, [rdi + rax*4 + 8] ;# jza  -  -  -
    movss  xmm14, [rdi + rcx*4 + 8] ;# jzc  -  -  -
    movss  xmm2,  [rdi + rbx*4 + 8] ;# jzb
    movss  xmm5,  [rdi + rdx*4 + 8] ;# jzd
    movlhps xmm13, xmm2 ;# jza  -  jzb  -
    movlhps xmm14, xmm5 ;# jzc  -  jzd -
    
    shufps xmm13, xmm14,  136  ;# 10001000 => jza jzb jzc jzd

    ;# xmm11: jxa jya jxb jyb 
    ;# xmm12: jxc jyc jxd jyd
    ;# xmm13: jza jzb jzc jzd

    xorps  xmm0, xmm0
    xorps  xmm4, xmm4
    xorps  xmm8, xmm8

    mulps  xmm3, [rsp + nb333_rinvH1]
    mulps  xmm7, [rsp + nb333_rinvH2]
    mulps  xmm10, [rsp + nb333_rinvM]
    
    subps  xmm0, xmm3
    subps  xmm4, xmm7
    subps  xmm8, xmm10
    
    movaps xmm1, xmm0
    movaps xmm2, xmm0
    movaps xmm3, xmm4
    movaps xmm5, xmm4
    movaps xmm6, xmm8
    movaps xmm7, xmm8

	mulps xmm0, [rsp + nb333_dxH1]
	mulps xmm1, [rsp + nb333_dyH1]
	mulps xmm2, [rsp + nb333_dzH1]
	mulps xmm3, [rsp + nb333_dxH2]
	mulps xmm4, [rsp + nb333_dyH2]
	mulps xmm5, [rsp + nb333_dzH2]
	mulps xmm6, [rsp + nb333_dxM]
	mulps xmm7, [rsp + nb333_dyM]
	mulps xmm8, [rsp + nb333_dzM]

    ;# fetch forces from O interaction
    movaps xmm14, [rsp + nb333_fjx]
    movaps xmm15, [rsp + nb333_fjy]
    addps  xmm13, [rsp + nb333_fjz]

    addps xmm14, xmm0
    addps xmm15, xmm1
    addps xmm13,  xmm2
    addps xmm0, [rsp + nb333_fixH1]
    addps xmm1, [rsp + nb333_fiyH1]
    addps xmm2, [rsp + nb333_fizH1]

    addps xmm14, xmm3
    addps xmm15, xmm4
    addps xmm13, xmm5
    addps xmm3, [rsp + nb333_fixH2]
    addps xmm4, [rsp + nb333_fiyH2]
    addps xmm5, [rsp + nb333_fizH2]

    addps xmm14, xmm6
    addps xmm15, xmm7
    addps xmm13, xmm8
    addps xmm6, [rsp + nb333_fixM]
    addps xmm7, [rsp + nb333_fiyM]
    addps xmm8, [rsp + nb333_fizM]

    movaps [rsp + nb333_fixH1], xmm0
    movaps [rsp + nb333_fiyH1], xmm1
    movaps [rsp + nb333_fizH1], xmm2
    movaps [rsp + nb333_fixH2], xmm3
    movaps [rsp + nb333_fiyH2], xmm4
    movaps [rsp + nb333_fizH2], xmm5
    movaps [rsp + nb333_fixM], xmm6
    movaps [rsp + nb333_fiyM], xmm7
    movaps [rsp + nb333_fizM], xmm8
    
    ;# xmm14 = fjx
    ;# xmm15 = fjy
    ;# xmm13 = fjz
    movaps xmm0, xmm14
    unpcklps xmm14, xmm15
    unpckhps xmm0,  xmm15
    
    addps  xmm11, xmm14
    addps  xmm12, xmm0
    
    movhlps  xmm14, xmm13 ;# fjzc fjzd
    
    movlps [rdi + rax*4], xmm11
    movhps [rdi + rbx*4], xmm11
    movlps [rdi + rcx*4], xmm12
    movhps [rdi + rdx*4], xmm12
    movss  [rdi + rax*4 + 8], xmm13
    movss  [rdi + rcx*4 + 8], xmm14
    shufps xmm13, xmm13, 1
    shufps xmm14, xmm14, 1
    movss  [rdi + rbx*4 + 8], xmm13
    movss  [rdi + rdx*4 + 8], xmm14
    
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb333_innerk],  4
	jl    .nb333_odd_inner
	jmp   .nb333_unroll_loop
.nb333_odd_inner:	
	add dword ptr [rsp + nb333_innerk],  4
	jnz   .nb333_odd_loop
	jmp   .nb333_updateouterdata
.nb333_odd_loop:
	mov   rdx, [rsp + nb333_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb333_innerjjnr],  4	

 	xorps xmm4, xmm4  	;# clear reg.
	movss xmm4, [rsp + nb333_iqM]
	mov rsi, [rbp + nb333_charge] 
	movhps xmm4, [rsp + nb333_iqH]  ;# [qM  0  qH  qH] 
	shufps xmm4, xmm4, 41	;# [0 qH qH qM]

	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb333_qqM], xmm3	;# use dummy qq for storage 
	
	xorps xmm6, xmm6
	mov rsi, [rbp + nb333_type]
	mov ebx, [rsi + rax*4]
	mov rsi, [rbp + nb333_vdwparam]
	shl ebx, 1	
	add ebx, [rsp + nb333_ntia]
	movlps xmm6, [rsi + rbx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 252  ;# 11111100
	shufps xmm7, xmm7, 253  ;# 11111101
	movaps [rsp + nb333_c6], xmm6
	movaps [rsp + nb333_c12], xmm7
	
	mov rsi, [rbp + nb333_pos]
	lea rax, [rax + rax*2]  

	movss xmm0, [rsp + nb333_ixO]
	movss xmm1, [rsp + nb333_iyO]
	movss xmm2, [rsp + nb333_izO]
	movss xmm3, [rsp + nb333_ixH1]
	movss xmm4, [rsp + nb333_iyH1]
	movss xmm5, [rsp + nb333_izH1]
	unpcklps xmm0, [rsp + nb333_ixH2] 	;# ixO ixH2 - -
	unpcklps xmm1, [rsp + nb333_iyH2]  	;# iyO iyH2 - -
	unpcklps xmm2, [rsp + nb333_izH2]	;# izO izH2 - -
	unpcklps xmm3, [rsp + nb333_ixM] 	;# ixH1 ixM - -
	unpcklps xmm4, [rsp + nb333_iyM]  	;# iyH1 iyM - -
	unpcklps xmm5, [rsp + nb333_izM]	;# izH1 izM - -
	unpcklps xmm0, xmm3  	;# ixO ixH1 ixH2 ixM
	unpcklps xmm1, xmm4 	;# same for y
	unpcklps xmm2, xmm5 	;# same for z
	
	;# move j coords to xmm0-xmm2 
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rax*4 + 4]
	movss xmm5, [rsi + rax*4 + 8]
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	;# use O distances for storage
	movaps [rsp + nb333_dxO], xmm3
	movaps [rsp + nb333_dyO], xmm4
	movaps [rsp + nb333_dzO], xmm5

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	;# rsq in xmm4 
	
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb333_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb333_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv

	movaps [rsp + nb333_rinvM], xmm0
	mulps  xmm4, xmm0  	;# r	 
	mulps xmm4, [rsp + nb333_tsc]
	
	movhlps xmm7, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm7	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm7, mm7
	movlhps xmm3, xmm7

	subps   xmm4, xmm3	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=eps2
	
	pslld mm6, 2
	pslld mm7, 2

	movd mm0, eax

	mov  rsi, [rbp + nb333_VFtab]
	movd eax, mm6
    	psrlq mm6, 32
	movd ebx, mm6
	movd ecx, mm7
	psrlq mm7, 32
	movd edx, mm7

	lea   rax, [rax + rax*2]
	lea   rbx, [rbx + rbx*2]
	lea   rcx, [rcx + rcx*2]
	lea   rdx, [rdx + rdx*2]

	;# first do LJ table for O
	;# load dispersion table data into xmm4
	movlps xmm4, [rsi + rax*4 + 16]
	movlps xmm6, [rsi + rax*4 + 24]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 0x1
	shufps xmm7, xmm7, 0x1

	;# dispersion table YFGH ready in xmm4-xmm7
	mulss  xmm6, xmm1   	;# xmm6=Geps 
	mulss  xmm7, xmm2   	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7   	;# xmm5=Fp 
	mulss  xmm7, [rsp + nb333_two]   	;# two*Heps2 
	addss  xmm7, xmm6
	addss  xmm7, xmm5 ;# xmm7=FF 
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 

	movaps xmm4, [rsp + nb333_c6]
	mulss  xmm7, xmm4	;# fijD 
	mulss  xmm5, xmm4	;# Vvdw6 

	;# save scalar force in xmm3. Update Vvdwtot directly 
	addss  xmm5, [rsp + nb333_Vvdwtot]
	xorps xmm3, xmm3
	movss xmm3, xmm7 ;# fscal 
	movss [rsp + nb333_Vvdwtot], xmm5
	
	;# load repulsion table data into xmm4
	movlps xmm4, [rsi + rax*4 + 32]
	movlps xmm6, [rsi + rax*4 + 40]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 0x1
	shufps xmm7, xmm7, 0x1
	;# repulsion table YFGH ready in xmm4-xmm7
	
	mulss  xmm6, xmm1   	;# xmm6=Geps 
	mulss  xmm7, xmm2   	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7   	;# xmm5=Fp 
	mulss  xmm7, [rsp + nb333_two]   	;# two*Heps2 
	addss  xmm7, xmm6
	addss  xmm7, xmm5 ;# xmm7=FF 
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 
 
	movaps xmm4, [rsp + nb333_c12]
	mulss  xmm7, xmm4 ;# fijR 
	mulss  xmm5, xmm4 ;# Vvdw12 
	addss  xmm3, xmm7
	
	addss  xmm5, [rsp + nb333_Vvdwtot]
	movss [rsp + nb333_Vvdwtot], xmm5

	movaps [rsp+nb333_rinvO], xmm3 ;# save fscal temp. in rinvO

	;# do the Coulomb interaction for H1,H2,M
	xorps  xmm5, xmm5
	movlps xmm3, [rsi + rcx*4]	;# data: Y3 F3  -  - 
	movhps xmm5, [rsi + rbx*4]	;# data:  0  0 Y2 F2
	movhps xmm3, [rsi + rdx*4]      ;# data: Y3 F3 Y4 F4 

	movaps xmm4, xmm5		;# data:  0  0 Y2 F2 
	shufps xmm4, xmm3, 0x88		;# data:  0 Y2 Y3 Y3
	shufps xmm5, xmm3, 0xDD	        ;# data:  0 F2 F3 F4 

	xorps  xmm7, xmm7
	movlps xmm3, [rsi + rcx*4 + 8]	;# data: G3 H3  -  - 
	movhps xmm7, [rsi + rbx*4 + 8]	;# data:  0  0 G2 H2
	movhps xmm3, [rsi + rdx*4 + 8]  ;# data: G3 H3 G4 H4 

	movaps xmm6, xmm7		;# data:  0  0 G2 H2 
	shufps xmm6, xmm3, 0x88		;# data:  0 G2 G3 G3
	shufps xmm7, xmm3, 0xDD	        ;# data:  0 H2 H3 H4 

	;# xmm4 =  0  Y2 Y3 Y4
	;# xmm5 =  0  F2 F3 F4
	;# xmm6 =  0  G2 G3 G4
	;# xmm7 =  0  H2 H3 H4
	;# coulomb table ready, in xmm4-xmm7      
	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp        
	mulps  xmm7, [rsp + nb333_two]   	;# two*Heps2 
	movaps xmm0, [rsp + nb333_qqM]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm0 ;# vcoul=qq*VV  
	mulps  xmm0, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and xmm0 fijC 
	;# increment vcoul - then we can get rid of mm5 
	addps  xmm5, [rsp + nb333_vctot]
	movaps [rsp + nb333_vctot], xmm5
	
	addps xmm0, [rsp+nb333_rinvO] ;# total fscal (temp. storage in rinvO)

	xorps xmm4, xmm4
	mulps  xmm0, [rsp + nb333_rinvM]
	mulps  xmm0, [rsp + nb333_tsc]
	subps  xmm4, xmm0
	
	movaps xmm0, [rsp + nb333_dxO]
	movaps xmm1, [rsp + nb333_dyO]
	movaps xmm2, [rsp + nb333_dzO]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4 ;# xmm0-xmm2 now contains tx-tz (partial force)
	
	movss  xmm3, [rsp + nb333_fixO]	
	movss  xmm4, [rsp + nb333_fiyO]	
	movss  xmm5, [rsp + nb333_fizO]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [rsp + nb333_fixO], xmm3	
	movss  [rsp + nb333_fiyO], xmm4	
	movss  [rsp + nb333_fizO], xmm5	;# updated the O force now do the H's
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2      
	shufps xmm3, xmm3, 0x39	;# shift right 
	shufps xmm4, xmm4, 0x39
	shufps xmm5, xmm5, 0x39
	addss  xmm3, [rsp + nb333_fixH1]
	addss  xmm4, [rsp + nb333_fiyH1]
	addss  xmm5, [rsp + nb333_fizH1]
	movss  [rsp + nb333_fixH1], xmm3	
	movss  [rsp + nb333_fiyH1], xmm4	
	movss  [rsp + nb333_fizH1], xmm5	;# updated the H1 force 

	shufps xmm3, xmm3, 0x39
	shufps xmm4, xmm4, 0x39
	shufps xmm5, xmm5, 0x39
	addss  xmm3, [rsp + nb333_fixH2]
	addss  xmm4, [rsp + nb333_fiyH2]
	addss  xmm5, [rsp + nb333_fizH2]
	movss  [rsp + nb333_fixH2], xmm3	
	movss  [rsp + nb333_fiyH2], xmm4	
	movss  [rsp + nb333_fizH2], xmm5	;# updated the H2 force 

	mov rdi, [rbp + nb333_faction]
	shufps xmm3, xmm3, 0x39
	shufps xmm4, xmm4, 0x39
	shufps xmm5, xmm5, 0x39
	addss  xmm3, [rsp + nb333_fixM]
	addss  xmm4, [rsp + nb333_fiyM]
	addss  xmm5, [rsp + nb333_fizM]
	movss  [rsp + nb333_fixM], xmm3	
	movss  [rsp + nb333_fiyM], xmm4	
	movss  [rsp + nb333_fizM], xmm5	;# updated the M force 

	movd eax, mm0
	;# the fj's - move in from mem start by acc. tx/ty/tz in xmm0, xmm1
	movlps xmm6, [rdi + rax*4]
	movss  xmm7, [rdi + rax*4 + 8]
	
	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  xmm0, xmm3
	movaps  xmm1, xmm4
	movaps  xmm2, xmm5
	
	shufps xmm3, xmm3, 0x39	;# shift right 
	shufps xmm4, xmm4, 0x39
	shufps xmm5, xmm5, 0x39
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5
	unpcklps xmm0, xmm1 	;# x,y sum in xmm0, z sum in xmm2
	
	addps    xmm6, xmm0
	addss    xmm7, xmm2
	
	movlps [rdi + rax*4],     xmm6
	movss  [rdi + rax*4 + 8], xmm7

	dec dword ptr [rsp + nb333_innerk]
	jz    .nb333_updateouterdata
	jmp   .nb333_odd_loop
.nb333_updateouterdata:
	mov   ecx, [rsp + nb333_ii3]
	mov   rdi, [rbp + nb333_faction]
	mov   rsi, [rbp + nb333_fshift]
	mov   edx, [rsp + nb333_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb333_fixO]
	movaps xmm1, [rsp + nb333_fiyO]
	movaps xmm2, [rsp + nb333_fizO]

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
	movaps xmm0, [rsp + nb333_fixH1]
	movaps xmm1, [rsp + nb333_fiyH1]
	movaps xmm2, [rsp + nb333_fizH1]

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
	movaps xmm0, [rsp + nb333_fixH2]
	movaps xmm1, [rsp + nb333_fiyH2]
	movaps xmm2, [rsp + nb333_fizH2]
	
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
	movaps xmm0, [rsp + nb333_fixM]
	movaps xmm1, [rsp + nb333_fiyM]
	movaps xmm2, [rsp + nb333_fizM]

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
	mov esi, [rsp + nb333_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb333_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb333_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb333_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb333_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb333_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb333_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb333_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb333_n], esi
        jmp .nb333_outer
.nb333_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb333_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb333_end
        ;# non-zero, do one more workunit
        jmp   .nb333_threadloop
.nb333_end:
	mov eax, [rsp + nb333_nouter]
	mov ebx, [rsp + nb333_ninner]
	mov rcx, [rbp + nb333_outeriter]
	mov rdx, [rbp + nb333_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 1104
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
	


.globl nb_kernel333nf_x86_64_sse
.globl _nb_kernel333nf_x86_64_sse
nb_kernel333nf_x86_64_sse:	
_nb_kernel333nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb333nf_fshift,         16
.equiv          nb333nf_gid,            24
.equiv          nb333nf_pos,            32
.equiv          nb333nf_faction,        40
.equiv          nb333nf_charge,         48
.equiv          nb333nf_p_facel,        56
.equiv          nb333nf_argkrf,         64
.equiv          nb333nf_argcrf,         72
.equiv          nb333nf_Vc,             80
.equiv          nb333nf_type,           88
.equiv          nb333nf_p_ntype,        96
.equiv          nb333nf_vdwparam,       104
.equiv          nb333nf_Vvdw,           112
.equiv          nb333nf_p_tabscale,     120
.equiv          nb333nf_VFtab,          128
.equiv          nb333nf_invsqrta,       136
.equiv          nb333nf_dvda,           144
.equiv          nb333nf_p_gbtabscale,   152
.equiv          nb333nf_GBtab,          160
.equiv          nb333nf_p_nthreads,     168
.equiv          nb333nf_count,          176
.equiv          nb333nf_mtx,            184
.equiv          nb333nf_outeriter,      192
.equiv          nb333nf_inneriter,      200
.equiv          nb333nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb333nf_ixO,            0
.equiv          nb333nf_iyO,            16
.equiv          nb333nf_izO,            32
.equiv          nb333nf_ixH1,           48
.equiv          nb333nf_iyH1,           64
.equiv          nb333nf_izH1,           80
.equiv          nb333nf_ixH2,           96
.equiv          nb333nf_iyH2,           112
.equiv          nb333nf_izH2,           128
.equiv          nb333nf_ixM,            144
.equiv          nb333nf_iyM,            160
.equiv          nb333nf_izM,            176
.equiv          nb333nf_iqM,            192
.equiv          nb333nf_iqH,            208
.equiv          nb333nf_qqM,            224
.equiv          nb333nf_qqH,            240
.equiv          nb333nf_rinvO,          256
.equiv          nb333nf_rinvH1,         272
.equiv          nb333nf_rinvH2,         288
.equiv          nb333nf_rinvM,          304
.equiv          nb333nf_rO,             320
.equiv          nb333nf_rH1,            336
.equiv          nb333nf_rH2,            352
.equiv          nb333nf_rM,             368
.equiv          nb333nf_tsc,            384
.equiv          nb333nf_c6,             400
.equiv          nb333nf_c12,            416
.equiv          nb333nf_vctot,          432
.equiv          nb333nf_Vvdwtot,        448
.equiv          nb333nf_half,           464
.equiv          nb333nf_three,          480
.equiv          nb333nf_is3,            496
.equiv          nb333nf_ii3,            500
.equiv          nb333nf_nri,            504
.equiv          nb333nf_iinr,           512
.equiv          nb333nf_jindex,         520
.equiv          nb333nf_jjnr,           528
.equiv          nb333nf_shift,          536
.equiv          nb333nf_shiftvec,       544
.equiv          nb333nf_facel,          552
.equiv          nb333nf_innerjjnr,      560
.equiv          nb333nf_ntia,           568
.equiv          nb333nf_innerk,         572
.equiv          nb333nf_n,              576
.equiv          nb333nf_nn1,            580
.equiv          nb333nf_nouter,         584
.equiv          nb333nf_ninner,         588
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

	sub rsp, 592		;# local variable stack space (n*16+8)
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
	mov [rsp + nb333nf_nouter], eax
	mov [rsp + nb333nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb333nf_nri], edi
	mov [rsp + nb333nf_iinr], rsi
	mov [rsp + nb333nf_jindex], rdx
	mov [rsp + nb333nf_jjnr], rcx
	mov [rsp + nb333nf_shift], r8
	mov [rsp + nb333nf_shiftvec], r9
	mov rsi, [rbp + nb333nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb333nf_facel], xmm0

	mov rax, [rbp + nb333nf_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb333nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb333nf_half], eax
	movss xmm1, [rsp + nb333nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb333nf_half],  xmm1
	movaps [rsp + nb333nf_three],  xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb333nf_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]		;# ebx =ii 

	mov   rdx, [rbp + nb333nf_charge]
	movss xmm4, [rdx + rbx*4 + 4]	
	movss xmm3, [rdx + rbx*4 + 12]	
	mov rsi, [rbp + nb333nf_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb333nf_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb333nf_iqM], xmm3
	movaps [rsp + nb333nf_iqH], xmm4
	
	mov   rdx, [rbp + nb333nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb333nf_p_ntype]
	imul  ecx, [rdi]  	;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb333nf_ntia], ecx		

.nb333nf_threadloop:
        mov   rsi, [rbp + nb333nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb333nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb333nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb333nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb333nf_n], eax
        mov [rsp + nb333nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb333nf_outerstart
        jmp .nb333nf_end
	
.nb333nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb333nf_nouter]
	mov [rsp + nb333nf_nouter], ebx

.nb333nf_outer:
	mov   rax, [rsp + nb333nf_shift]  	;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb333nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb333nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb333nf_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2	
	movaps xmm6, xmm0
	movaps xmm7, xmm1
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb333nf_pos]	;# rax = base of pos[]  
	mov   [rsp + nb333nf_ii3], ebx

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
	movaps [rsp + nb333nf_ixO], xmm3
	movaps [rsp + nb333nf_iyO], xmm4
	movaps [rsp + nb333nf_izO], xmm5
	movaps [rsp + nb333nf_ixH1], xmm6
	movaps [rsp + nb333nf_iyH1], xmm7

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
	movaps [rsp + nb333nf_izH1], xmm6
	movaps [rsp + nb333nf_ixH2], xmm0
	movaps [rsp + nb333nf_iyH2], xmm1
	movaps [rsp + nb333nf_izH2], xmm2
	movaps [rsp + nb333nf_ixM], xmm3
	movaps [rsp + nb333nf_iyM], xmm4
	movaps [rsp + nb333nf_izM], xmm5
	
	;# clear vctot 
	xorps xmm4, xmm4
	movaps [rsp + nb333nf_vctot], xmm4
	movaps [rsp + nb333nf_Vvdwtot], xmm4

	mov   rax, [rsp + nb333nf_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb333nf_pos]
	mov   rdi, [rbp + nb333nf_faction]	
	mov   rax, [rsp + nb333nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb333nf_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb333nf_ninner]
	mov   [rsp + nb333nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb333nf_innerk], edx	;# number of innerloop atoms 
	jge   .nb333nf_unroll_loop
	jmp   .nb333nf_odd_inner
.nb333nf_unroll_loop:
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb333nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4
	
	add qword ptr [rsp + nb333nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb333nf_charge]	;# base of charge[] 

	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]
	
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	 	;# and in xmm4 
	mulps  xmm3, [rsp + nb333nf_iqM]
	mulps  xmm4, [rsp + nb333nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps  [rsp + nb333nf_qqM], xmm3
	movaps  [rsp + nb333nf_qqH], xmm4
	
	mov rsi, [rbp + nb333nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]
	mov rsi, [rbp + nb333nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb333nf_ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm6, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm6, xmm7, 221  ;# 11011101
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [rsp + nb333nf_c6], xmm4
	movaps [rsp + nb333nf_c12], xmm6

	mov rsi, [rbp + nb333nf_pos]   	;# base of pos[] 

	lea   rax, [rax + rax*2] 	;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2] 	;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	

	;# move four coordinates to xmm0-xmm2 	
	movlps xmm4, [rsi + rax*4]
	movlps xmm5, [rsi + rcx*4]
	movss xmm2, [rsi + rax*4 + 8]
	movss xmm6, [rsi + rcx*4 + 8]

	movhps xmm4, [rsi + rbx*4]
	movhps xmm5, [rsi + rdx*4]

	movss xmm0, [rsi + rbx*4 + 8]
	movss xmm1, [rsi + rdx*4 + 8]

	shufps xmm2, xmm0, 0
	shufps xmm6, xmm1, 0
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 136  ;# 10001000
	
	shufps xmm0, xmm5, 136  ;# 10001000
	shufps xmm1, xmm5, 221  ;# 11011101		

	;# move ixO-izO to xmm4-xmm6 
	movaps xmm4, [rsp + nb333nf_ixO]
	movaps xmm5, [rsp + nb333nf_iyO]
	movaps xmm6, [rsp + nb333nf_izO]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	movaps xmm7, xmm4
	;# rsqO in xmm7

	;# move ixH1-izH1 to xmm4-xmm6 
	movaps xmm4, [rsp + nb333nf_ixH1]
	movaps xmm5, [rsp + nb333nf_iyH1]
	movaps xmm6, [rsp + nb333nf_izH1]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm6, xmm5
	addps xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movaps xmm3, [rsp + nb333nf_ixH2]
	movaps xmm4, [rsp + nb333nf_iyH2]
	movaps xmm5, [rsp + nb333nf_izH2]

	;# calc dr 
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	;# square it 
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	addps xmm5, xmm4
	addps xmm5, xmm3
	
	;# move ixM-izM to xmm2-xmm4  
	movaps xmm3, [rsp + nb333nf_iyM]
	movaps xmm4, [rsp + nb333nf_izM]
	subps  xmm3, xmm1
	subps  xmm4, xmm2
	movaps xmm2, [rsp + nb333nf_ixM]
	subps  xmm2, xmm0	

	;# square it 
	mulps xmm2,xmm2
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	addps xmm4, xmm3
	addps xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7

	;# rsqH1 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm0, [rsp + nb333nf_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm0, [rsp + nb333nf_half]
	movaps  [rsp + nb333nf_rinvH1], xmm0	;# rinvH1  
	mulps   xmm6, xmm0
	movaps  [rsp + nb333nf_rH1], xmm6

	;# rsqH2 - seed to xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm0, [rsp + nb333nf_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm0, [rsp + nb333nf_half]
	movaps  [rsp + nb333nf_rinvH2], xmm0	;# rinvH2 
	mulps   xmm5, xmm0
	movaps  [rsp + nb333nf_rH2], xmm5

	;# rsqM - seed to xmm2 
	rsqrtps xmm2, xmm4
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm0, [rsp + nb333nf_three]
	mulps   xmm2, xmm4	;# rsq*lu*lu 
	subps   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm0, [rsp + nb333nf_half]
	movaps  [rsp + nb333nf_rinvM], xmm0	;# rinvM 
	mulps   xmm4, xmm0
	movaps  [rsp + nb333nf_rM], xmm4	
	
	;# Do the O LJ table interaction directly.
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm0, [rsp + nb333nf_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm0, [rsp + nb333nf_half] ;# rinv
	
	movaps xmm1, xmm0
	mulps  xmm1, xmm7	;# xmm1=r
	mulps  xmm1, [rsp + nb333nf_tsc] ;# r*tabscale
	
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

	mov  rsi, [rbp + nb333nf_VFtab]
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
	movhps xmm7, [rsi + rdx*4 + 16] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 24]
	movlps xmm3, [rsi + rcx*4 + 24]
	movhps xmm7, [rsi + rbx*4 + 24]
	movhps xmm3, [rsi + rdx*4 + 24] ;# other half of coulomb table  
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

	movaps xmm4, [rsp + nb333nf_c6]
	mulps  xmm5, xmm4	;# Vvdw6 

	;# Update Vvdwtot directly	
	addps  xmm5, [rsp + nb333nf_Vvdwtot]
	movaps [rsp + nb333nf_Vvdwtot], xmm5

	;# load repulsion table data into xmm4-xmm7
	movlps xmm5, [rsi + rax*4 + 32]
	movlps xmm7, [rsi + rcx*4 + 32]
	movhps xmm5, [rsi + rbx*4 + 32]
	movhps xmm7, [rsi + rdx*4 + 32] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 40]
	movlps xmm3, [rsi + rcx*4 + 40]
	movhps xmm7, [rsi + rbx*4 + 40]
	movhps xmm3, [rsi + rdx*4 + 40] ;# other half of coulomb table  
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
 
	movaps xmm4, [rsp + nb333nf_c12]
	mulps  xmm5, xmm4 ;# Vvdw12 
	addps  xmm5, [rsp + nb333nf_Vvdwtot]
	movaps [rsp + nb333nf_Vvdwtot], xmm5

	;# Do H1 interaction
	mov  rsi, [rbp + nb333nf_VFtab]
	
	movaps xmm7, [rsp + nb333nf_rH1]
	mulps   xmm7, [rsp + nb333nf_tsc]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4	;# mm6/mm7 contain lu indices 
	
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm4, mm7
	movlhps xmm3, xmm4
	
	subps xmm7, xmm3
	movaps xmm1, xmm7	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 2
	pslld mm7, 2
	
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
	movaps xmm0, [rsp + nb333nf_qqH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm0 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 
	;# increment vcoul 
	addps  xmm5, [rsp + nb333nf_vctot]
	movaps [rsp + nb333nf_vctot], xmm5 

	;# Done with H1, do H2 interactions 
	movaps xmm7, [rsp + nb333nf_rH2]
	mulps   xmm7, [rsp + nb333nf_tsc]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4	;# mm6/mm7 contain lu indices 
	
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm4, mm7
	movlhps xmm3, xmm4
	
	subps xmm7, xmm3
	movaps xmm1, xmm7	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 2
	pslld mm7, 2
	
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
	movaps xmm0, [rsp + nb333nf_qqH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm0 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 
	;# increment vcoul 
	addps  xmm5, [rsp + nb333nf_vctot]
	movaps [rsp + nb333nf_vctot], xmm5 

	;# Done with H2, do M interactions 
	movaps xmm7, [rsp + nb333nf_rM]
	mulps   xmm7, [rsp + nb333nf_tsc]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4	;# mm6/mm7 contain lu indices 
	
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm4, mm7
	movlhps xmm3, xmm4
	
	subps xmm7, xmm3
	movaps xmm1, xmm7	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 2
	pslld mm7, 2
	
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
	movaps xmm0, [rsp + nb333nf_qqM]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm0 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 
	;# increment vcoul 
	addps  xmm5, [rsp + nb333nf_vctot]
	movaps [rsp + nb333nf_vctot], xmm5 
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb333nf_innerk],  4
	jl    .nb333nf_odd_inner
	jmp   .nb333nf_unroll_loop
.nb333nf_odd_inner:	
	add dword ptr [rsp + nb333nf_innerk],  4
	jnz   .nb333nf_odd_loop
	jmp   .nb333nf_updateouterdata
.nb333nf_odd_loop:
	mov   rdx, [rsp + nb333nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb333nf_innerjjnr],  4	

 	xorps xmm4, xmm4  	;# clear reg.
	movss xmm4, [rsp + nb333nf_iqM]
	mov rsi, [rbp + nb333nf_charge] 
	movhps xmm4, [rsp + nb333nf_iqH]  ;# [qM  0  qH  qH] 
	shufps xmm4, xmm4, 41	;# [0 qH qH qM]

	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb333nf_qqM], xmm3	;# use dummy qq for storage 
	
	xorps xmm6, xmm6
	mov rsi, [rbp + nb333nf_type]
	mov ebx, [rsi + rax*4]
	mov rsi, [rbp + nb333nf_vdwparam]
	shl ebx, 1	
	add ebx, [rsp + nb333nf_ntia]
	movlps xmm6, [rsi + rbx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 252  ;# 11111100
	shufps xmm7, xmm7, 253  ;# 11111101
	movaps [rsp + nb333nf_c6], xmm6
	movaps [rsp + nb333nf_c12], xmm7
	
	mov rsi, [rbp + nb333nf_pos]
	lea rax, [rax + rax*2]  

	movss xmm3, [rsp + nb333nf_ixO]
	movss xmm4, [rsp + nb333nf_iyO]
	movss xmm5, [rsp + nb333nf_izO]
	movss xmm0, [rsp + nb333nf_ixH1]
	movss xmm1, [rsp + nb333nf_iyH1]
	movss xmm2, [rsp + nb333nf_izH1]
	unpcklps xmm3, [rsp + nb333nf_ixH2] 	;# ixO ixH2 - -
	unpcklps xmm4, [rsp + nb333nf_iyH2]  	;# iyO iyH2 - -
	unpcklps xmm5, [rsp + nb333nf_izH2]	;# izO izH2 - -
	unpcklps xmm0, [rsp + nb333nf_ixM] 	;# ixH1 ixM - -
	unpcklps xmm1, [rsp + nb333nf_iyM]  	;# iyH1 iyM - -
	unpcklps xmm2, [rsp + nb333nf_izM]	;# izH1 izM - -
	unpcklps xmm3, xmm0  	;# ixO ixH1 ixH2 ixM
	unpcklps xmm4, xmm1 	;# same for y
	unpcklps xmm5, xmm2 	;# same for z
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [rsi + rax*4]
	movss xmm1, [rsi + rax*4 + 4]
	movss xmm2, [rsi + rax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	;# rsq in xmm4 
	
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb333nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb333nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv

	movaps [rsp + nb333nf_rinvM], xmm0
	mulps  xmm4, xmm0  	;# r	 
	mulps xmm4, [rsp + nb333nf_tsc]
	
	movhlps xmm7, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm7	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm7, mm7
	movlhps xmm3, xmm7

	subps   xmm4, xmm3	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=eps2
	
	pslld mm6, 2
	pslld mm7, 2

	movd mm0, eax

	mov  rsi, [rbp + nb333nf_VFtab]
	movd eax, mm6
    	psrlq mm6, 32
	movd ebx, mm6
	movd ecx, mm7
	psrlq mm7, 32
	movd edx, mm7

	lea   rax, [rax + rax*2]
	lea   rbx, [rbx + rbx*2]
	lea   rcx, [rcx + rcx*2]
	lea   rdx, [rdx + rdx*2]

	;# first do LJ table for O
	;# load dispersion table data into xmm4
	movlps xmm4, [rsi + rax*4 + 16]
	movlps xmm6, [rsi + rax*4 + 24]
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

	movaps xmm4, [rsp + nb333nf_c6]
	mulss  xmm5, xmm4	;# Vvdw6 
	addss  xmm5, [rsp + nb333nf_Vvdwtot]
	movss [rsp + nb333nf_Vvdwtot], xmm5

	;# load repulsion table data into xmm4
	movlps xmm4, [rsi + rax*4 + 32]
	movlps xmm6, [rsi + rax*4 + 40]
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
 
	movaps xmm4, [rsp + nb333nf_c12]
	mulss  xmm5, xmm4 ;# Vvdw12 
	addss  xmm5, [rsp + nb333nf_Vvdwtot]
	movss [rsp + nb333nf_Vvdwtot], xmm5
	;# do the Coulomb interaction for H1,H2,M
	xorps  xmm5, xmm5
	movlps xmm3, [rsi + rcx*4]	;# data: Y3 F3  -  - 
	movhps xmm5, [rsi + rbx*4]	;# data:  0  0 Y2 F2
	movhps xmm3, [rsi + rdx*4]      ;# data: Y3 F3 Y4 F4 

	movaps xmm4, xmm5		;# data:  0  0 Y2 F2 
	shufps xmm4, xmm3, 0x88		;# data:  0 Y2 Y3 Y3
	shufps xmm5, xmm3, 0xDD	        ;# data:  0 F2 F3 F4 

	xorps  xmm7, xmm7
	movlps xmm3, [rsi + rcx*4 + 8]	;# data: G3 H3  -  - 
	movhps xmm7, [rsi + rbx*4 + 8]	;# data:  0  0 G2 H2
	movhps xmm3, [rsi + rdx*4 + 8]  ;# data: G3 H3 G4 H4 

	movaps xmm6, xmm7		;# data:  0  0 G2 H2 
	shufps xmm6, xmm3, 0x88		;# data:  0 G2 G3 G3
	shufps xmm7, xmm3, 0xDD	        ;# data:  0 H2 H3 H4 

	;# xmm4 =  0  Y2 Y3 Y4
	;# xmm5 =  0  F2 F3 F4
	;# xmm6 =  0  G2 G3 G4
	;# xmm7 =  0  H2 H3 H4	
	;# coulomb table ready, in xmm4-xmm7      
	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp        
	movaps xmm0, [rsp + nb333nf_qqM]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm0 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 
	;# increment vcoul
	addps  xmm5, [rsp + nb333nf_vctot]
	movaps [rsp + nb333nf_vctot], xmm5

	dec dword ptr [rsp + nb333nf_innerk]
	jz    .nb333nf_updateouterdata
	jmp   .nb333nf_odd_loop
.nb333nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb333nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb333nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb333nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb333nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb333nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb333nf_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb333nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb333nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb333nf_n], esi
        jmp .nb333nf_outer
.nb333nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb333nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb333nf_end
        ;# non-zero, do one more workunit
        jmp   .nb333nf_threadloop
.nb333nf_end:
	mov eax, [rsp + nb333nf_nouter]
	mov ebx, [rsp + nb333nf_ninner]
	mov rcx, [rbp + nb333nf_outeriter]
	mov rdx, [rbp + nb333nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 592
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
