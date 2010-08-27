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

	
.globl nb_kernel133_x86_64_sse
.globl _nb_kernel133_x86_64_sse
nb_kernel133_x86_64_sse:	
_nb_kernel133_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb133_fshift,           16
.equiv          nb133_gid,              24
.equiv          nb133_pos,              32
.equiv          nb133_faction,          40
.equiv          nb133_charge,           48
.equiv          nb133_p_facel,          56
.equiv          nb133_argkrf,           64
.equiv          nb133_argcrf,           72
.equiv          nb133_Vc,               80
.equiv          nb133_type,             88
.equiv          nb133_p_ntype,          96
.equiv          nb133_vdwparam,         104
.equiv          nb133_Vvdw,             112
.equiv          nb133_p_tabscale,       120
.equiv          nb133_VFtab,            128
.equiv          nb133_invsqrta,         136
.equiv          nb133_dvda,             144
.equiv          nb133_p_gbtabscale,     152
.equiv          nb133_GBtab,            160
.equiv          nb133_p_nthreads,       168
.equiv          nb133_count,            176
.equiv          nb133_mtx,              184
.equiv          nb133_outeriter,        192
.equiv          nb133_inneriter,        200
.equiv          nb133_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb133_ixO,              0
.equiv          nb133_iyO,              16
.equiv          nb133_izO,              32
.equiv          nb133_ixH1,             48
.equiv          nb133_iyH1,             64
.equiv          nb133_izH1,             80
.equiv          nb133_ixH2,             96
.equiv          nb133_iyH2,             112
.equiv          nb133_izH2,             128
.equiv          nb133_ixM,              144
.equiv          nb133_iyM,              160
.equiv          nb133_izM,              176
.equiv          nb133_iqM,              192
.equiv          nb133_iqH,              208
.equiv          nb133_dxO,              224
.equiv          nb133_dyO,              240
.equiv          nb133_dzO,              256
.equiv          nb133_dxH1,             272
.equiv          nb133_dyH1,             288
.equiv          nb133_dzH1,             304
.equiv          nb133_dxH2,             320
.equiv          nb133_dyH2,             336
.equiv          nb133_dzH2,             352
.equiv          nb133_dxM,              368
.equiv          nb133_dyM,              384
.equiv          nb133_dzM,              400
.equiv          nb133_qqM,              416
.equiv          nb133_qqH,              432
.equiv          nb133_rinvH1,           448
.equiv          nb133_rinvH2,           464
.equiv          nb133_rinvM,            480
.equiv          nb133_two,              496
.equiv          nb133_c6,               512
.equiv          nb133_c12,              528
.equiv          nb133_tsc,              544
.equiv          nb133_fstmp,            560
.equiv          nb133_krf,              576
.equiv          nb133_crf,              592
.equiv          nb133_krsqH1,           608
.equiv          nb133_krsqH2,           624
.equiv          nb133_krsqM,            640
.equiv          nb133_vctot,            656
.equiv          nb133_Vvdwtot,          672
.equiv          nb133_fixO,             688
.equiv          nb133_fiyO,             704
.equiv          nb133_fizO,             720
.equiv          nb133_fixH1,            736
.equiv          nb133_fiyH1,            752
.equiv          nb133_fizH1,            768
.equiv          nb133_fixH2,            784
.equiv          nb133_fiyH2,            800
.equiv          nb133_fizH2,            816
.equiv          nb133_fixM,             832
.equiv          nb133_fiyM,             848
.equiv          nb133_fizM,             864
.equiv          nb133_fjx,              880
.equiv          nb133_fjy,              896
.equiv          nb133_fjz,              912
.equiv          nb133_half,             928
.equiv          nb133_three,            944
.equiv          nb133_rsqOO,            960
.equiv          nb133_facel,            976
.equiv          nb133_iinr,             984
.equiv          nb133_jindex,           992
.equiv          nb133_jjnr,             1000
.equiv          nb133_shift,            1008
.equiv          nb133_shiftvec,         1016
.equiv          nb133_innerjjnr,        1024
.equiv          nb133_is3,              1032
.equiv          nb133_ii3,              1036
.equiv          nb133_nri,              1040
.equiv          nb133_ntia,             1044
.equiv          nb133_innerk,           1048
.equiv          nb133_n,                1052
.equiv          nb133_nn1,              1056
.equiv          nb133_nouter,           1060
.equiv          nb133_ninner,           1064
	
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
	sub rsp, 1072		;# local variable stack space (n*16+8)
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
	mov [rsp + nb133_nouter], eax
	mov [rsp + nb133_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb133_nri], edi
	mov [rsp + nb133_iinr], rsi
	mov [rsp + nb133_jindex], rdx
	mov [rsp + nb133_jjnr], rcx
	mov [rsp + nb133_shift], r8
	mov [rsp + nb133_shiftvec], r9
	mov rsi, [rbp + nb133_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb133_facel], xmm0

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb133_half], eax
	movss xmm1, [rsp + nb133_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb133_half],  xmm1
	movaps [rsp + nb133_two],  xmm2
	movaps [rsp + nb133_three],  xmm3
	
	mov rax, [rbp + nb133_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb133_tsc], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb133_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb133_charge]
	movss xmm4, [rdx + rbx*4 + 4]	
	movss xmm3, [rdx + rbx*4 + 12]	
	mov rsi, [rbp + nb133_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb133_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb133_iqM], xmm3
	movaps [rsp + nb133_iqH], xmm4
	
	mov   rdx, [rbp + nb133_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb133_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb133_ntia], ecx		
.nb133_threadloop:
        mov   rsi, [rbp + nb133_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb133_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb133_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb133_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb133_n], eax
        mov [rsp + nb133_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb133_outerstart
        jmp .nb133_end
	
.nb133_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb133_nouter]
	mov [rsp + nb133_nouter], ebx

.nb133_outer:
	mov   rax, [rsp + nb133_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb133_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb133_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb133_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2	
	movaps xmm6, xmm0
	movaps xmm7, xmm1
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb133_pos]	;# eax = base of pos[]  
	mov   [rsp + nb133_ii3], ebx

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
	movaps [rsp + nb133_ixO], xmm3
	movaps [rsp + nb133_iyO], xmm4
	movaps [rsp + nb133_izO], xmm5
	movaps [rsp + nb133_ixH1], xmm6
	movaps [rsp + nb133_iyH1], xmm7

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
	movaps [rsp + nb133_izH1], xmm6
	movaps [rsp + nb133_ixH2], xmm0
	movaps [rsp + nb133_iyH2], xmm1
	movaps [rsp + nb133_izH2], xmm2
	movaps [rsp + nb133_ixM], xmm3
	movaps [rsp + nb133_iyM], xmm4
	movaps [rsp + nb133_izM], xmm5
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb133_vctot], xmm4
	movaps [rsp + nb133_Vvdwtot], xmm4
	movaps [rsp + nb133_fixO], xmm4
	movaps [rsp + nb133_fiyO], xmm4
	movaps [rsp + nb133_fizO], xmm4
	movaps [rsp + nb133_fixH1], xmm4
	movaps [rsp + nb133_fiyH1], xmm4
	movaps [rsp + nb133_fizH1], xmm4
	movaps [rsp + nb133_fixH2], xmm4
	movaps [rsp + nb133_fiyH2], xmm4
	movaps [rsp + nb133_fizH2], xmm4
	movaps [rsp + nb133_fixM], xmm4
	movaps [rsp + nb133_fiyM], xmm4
	movaps [rsp + nb133_fizM], xmm4
	
	mov   rax, [rsp + nb133_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb133_pos]
	mov   rdi, [rbp + nb133_faction]	
	mov   rax, [rsp + nb133_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb133_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb133_ninner]
	mov   [rsp + nb133_ninner], ecx
	add   edx, 0
	mov   [rsp + nb133_innerk], edx	;# number of innerloop atoms 
	jge   .nb133_unroll_loop
	jmp   .nb133_odd_inner
.nb133_unroll_loop:
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb133_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb133_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb133_charge]	;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	 	;# and in xmm4 
	mulps  xmm3, [rsp + nb133_iqM]
	mulps  xmm4, [rsp + nb133_iqH]

	movaps  [rsp + nb133_qqM], xmm3
	movaps  [rsp + nb133_qqH], xmm4
	
	mov rsi, [rbp + nb133_type]
	mov r8d, [rsi + rax*4]
	mov r9d, [rsi + rbx*4]
	mov r10d, [rsi + rcx*4]
	mov r11d, [rsi + rdx*4]
	mov rsi, [rbp + nb133_vdwparam]
	shl r8d, 1	
	shl r9d, 1	
	shl r10d, 1	
	shl r11d, 1	
	mov edi, [rsp + nb133_ntia]
	add r8d, edi
	add r9d, edi
	add r10d, edi
	add r11d, edi

	movlps xmm6, [rsi + r8*4]
	movlps xmm7, [rsi + r10*4]
	movhps xmm6, [rsi + r9*4]
	movhps xmm7, [rsi + r11*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm6, xmm7, 221  ;# constant 11011101

	movaps [rsp + nb133_c6], xmm4
	movaps [rsp + nb133_c12], xmm6

	mov rsi, [rbp + nb133_pos]   	;# base of pos[] 

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

	shufps xmm2, xmm6, 136  ;# constant 10001000	
	shufps xmm0, xmm5, 136  ;# constant 10001000
	shufps xmm1, xmm5, 221  ;# constant 11011101		

    ;# xmm0 = jx
    ;# xmm1 = jy
    ;# xmm2 = jz
    
    ;# O interaction
    ;# copy to xmm3-xmm5
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    
    subps xmm3, [rsp + nb133_ixO]
    subps xmm4, [rsp + nb133_iyO]
    subps xmm5, [rsp + nb133_izO]
    
    movaps [rsp + nb133_dxO], xmm3
    movaps [rsp + nb133_dyO], xmm4
    movaps [rsp + nb133_dzO], xmm5
    
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
	movaps xmm4, [rsp + nb133_three]
	mulps xmm5, xmm3	;# rsq*lu*lu 	
    subps xmm4, xmm5	;# 30-rsq*lu*lu 
	mulps xmm4, xmm15	
	mulps xmm4, [rsp + nb133_half]	
	movaps xmm15, xmm4
	mulps  xmm3, xmm4	
    ;# xmm15=rinv
    ;# xmm3=r

    mulps xmm3, [rsp + nb133_tsc] ;# rtab

    ;# truncate and convert to integers
    cvttps2dq xmm5, xmm3
    
    ;# convert back to float
    cvtdq2ps  xmm4, xmm5
    
    ;# multiply by 8
    pslld   xmm5, 3

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

	mov rsi, [rbp + nb133_VFtab]
    ;# calculate LJ table
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
    
	movlps xmm13, [rsi + r10*4 + 8]
	movlps xmm14, [rsi + r10*4 + 24]

	movhps xmm7,  [rsi + r9*4 + 8]
	movhps xmm11, [rsi + r9*4 + 24]
    
	movhps xmm13, [rsi + r11*4 + 8]
	movhps xmm14, [rsi + r11*4 + 24]

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
    movaps xmm12, [rsp + nb133_c6]
    movaps xmm13, [rsp + nb133_c12]
    addps  xmm5, xmm4 ;# VV
    addps  xmm9, xmm8

    mulps  xmm5, xmm12  ;# VV*c6 = vnb6
    mulps  xmm9, xmm13  ;# VV*c12 = vnb12
    addps  xmm5, xmm9
    addps  xmm5, [rsp + nb133_Vvdwtot]
    movaps [rsp + nb133_Vvdwtot], xmm5
        
    mulps  xmm7, xmm12   ;# FF*c6 = fnb6
    mulps  xmm11, xmm13   ;# FF*c12  = fnb12
    addps  xmm7, xmm11
    
    mulps  xmm7, [rsp + nb133_tsc]
    mulps  xmm7, xmm15   ;# -fscal
    xorps  xmm9, xmm9
    
    subps  xmm9, xmm7     ;# fscal
    movaps xmm10, xmm9
    movaps xmm11, xmm9

    mulps  xmm9,  [rsp + nb133_dxO] ;# fx/fy/fz
    mulps  xmm10, [rsp + nb133_dyO]
    mulps  xmm11, [rsp + nb133_dzO]

    ;# save j force temporarily
    movaps [rsp + nb133_fjx], xmm9
    movaps [rsp + nb133_fjy], xmm10
    movaps [rsp + nb133_fjz], xmm11
    
    ;# increment i O force
    addps xmm9, [rsp + nb133_fixO]
    addps xmm10, [rsp + nb133_fiyO]
    addps xmm11, [rsp + nb133_fizO]
    movaps [rsp + nb133_fixO], xmm9
    movaps [rsp + nb133_fiyO], xmm10
    movaps [rsp + nb133_fizO], xmm11
    ;# finished O LJ interaction.

    ;# do H1, H2, and M interactions in parallel.
    ;# xmm0-xmm2 still contain j coordinates.        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + nb133_ixH1]
    subps xmm1, [rsp + nb133_iyH1]
    subps xmm2, [rsp + nb133_izH1]
    subps xmm3, [rsp + nb133_ixH2]
    subps xmm4, [rsp + nb133_iyH2]
    subps xmm5, [rsp + nb133_izH2]
    subps xmm6, [rsp + nb133_ixM]
    subps xmm7, [rsp + nb133_iyM]
    subps xmm8, [rsp + nb133_izM]
    
	movaps [rsp + nb133_dxH1], xmm0
	movaps [rsp + nb133_dyH1], xmm1
	movaps [rsp + nb133_dzH1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb133_dxH2], xmm3
	movaps [rsp + nb133_dyH2], xmm4
	movaps [rsp + nb133_dzH2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb133_dxM], xmm6
	movaps [rsp + nb133_dyM], xmm7
	movaps [rsp + nb133_dzM], xmm8
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
		
	movaps  xmm9, [rsp + nb133_three]
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

	movaps  xmm0, [rsp + nb133_half]
	mulps   xmm9, xmm0  ;# rinvH1
	mulps   xmm10, xmm0 ;# rinvH2
    mulps   xmm11, xmm0 ;# rinvM
	
	;# interactions 
    movaps xmm0, xmm9    ;# rinv
    movaps xmm1, xmm10
    movaps xmm2, xmm11
    mulps  xmm9, xmm9    ;# rinvsq
    mulps  xmm10, xmm10
    mulps  xmm11, xmm11
    mulps  xmm0, [rsp + nb133_qqH] 
    mulps  xmm1, [rsp + nb133_qqH] 
    mulps  xmm2, [rsp + nb133_qqM] 
    mulps  xmm9, xmm0
    mulps  xmm10, xmm1
    mulps  xmm11, xmm2
    
    addps xmm0, [rsp + nb133_vctot] 
    addps xmm1, xmm2
    addps xmm0, xmm1
    movaps [rsp + nb133_vctot], xmm0
    
	;# move j forces to local temp variables 
	mov rdi, [rbp + nb133_faction]
    movlps xmm0, [rdi + rax*4] ;# jxa jya  -   -
    movlps xmm1, [rdi + rcx*4] ;# jxc jyc  -   -
    movhps xmm0, [rdi + rbx*4] ;# jxa jya jxb jyb 
    movhps xmm1, [rdi + rdx*4] ;# jxc jyc jxd jyd 

    movss  xmm2, [rdi + rax*4 + 8] ;# jza  -  -  -
    movss  xmm3, [rdi + rcx*4 + 8] ;# jzc  -  -  -
    movss  xmm5, [rdi + rbx*4 + 8] ;# jzb  - - -
    movss  xmm6, [rdi + rdx*4 + 8] ;# jzd - - -
    movlhps xmm2, xmm5
    movlhps xmm3, xmm6
    
    shufps xmm2, xmm3,  136  ;# 10001000 => jza jzb jzc jzd

    ;# xmm0: jxa jya jxb jyb 
    ;# xmm1: jxc jyc jxd jyd
    ;# xmm2: jza jzb jzc jzd

    movaps xmm7, xmm9
    movaps xmm8, xmm9
    movaps xmm13, xmm11
    movaps xmm14, xmm11
    movaps xmm15, xmm11
    movaps xmm11, xmm10
    movaps xmm12, xmm10

	mulps xmm7, [rsp + nb133_dxH1]
	mulps xmm8, [rsp + nb133_dyH1]
	mulps xmm9, [rsp + nb133_dzH1]
	mulps xmm10, [rsp + nb133_dxH2]
	mulps xmm11, [rsp + nb133_dyH2]
	mulps xmm12, [rsp + nb133_dzH2]
	mulps xmm13, [rsp + nb133_dxM]
	mulps xmm14, [rsp + nb133_dyM]
	mulps xmm15, [rsp + nb133_dzM]

    ;# fetch forces from O interaction
    movaps xmm3, [rsp + nb133_fjx]
    movaps xmm4, [rsp + nb133_fjy]
    addps  xmm2, [rsp + nb133_fjz]

    addps xmm3, xmm7
    addps xmm4, xmm8
    addps xmm2, xmm9
    addps xmm7, [rsp + nb133_fixH1]
    addps xmm8, [rsp + nb133_fiyH1]
    addps xmm9, [rsp + nb133_fizH1]

    addps xmm3, xmm10
    addps xmm4, xmm11
    addps xmm2, xmm12
    addps xmm10, [rsp + nb133_fixH2]
    addps xmm11, [rsp + nb133_fiyH2]
    addps xmm12, [rsp + nb133_fizH2]

    addps xmm3, xmm13
    addps xmm4, xmm14
    addps xmm2, xmm15
    addps xmm13, [rsp + nb133_fixM]
    addps xmm14, [rsp + nb133_fiyM]
    addps xmm15, [rsp + nb133_fizM]

    movaps [rsp + nb133_fixH1], xmm7
    movaps [rsp + nb133_fiyH1], xmm8
    movaps [rsp + nb133_fizH1], xmm9
    movaps [rsp + nb133_fixH2], xmm10
    movaps [rsp + nb133_fiyH2], xmm11
    movaps [rsp + nb133_fizH2], xmm12
    movaps [rsp + nb133_fixM], xmm13
    movaps [rsp + nb133_fiyM], xmm14
    movaps [rsp + nb133_fizM], xmm15
    
    ;# xmm3 = fjx , xmm4 = fjy  , xmm2=fjz, already updated.
    movaps xmm5, xmm3   
    unpcklps xmm3, xmm4   ;# fjx1 fjy1 fjx2 fjy2
    unpckhps xmm5, xmm4   ;# fjx3 fjy3 fjx4 fjy4
    
    addps xmm0, xmm3 
    addps xmm1, xmm5
    movhlps  xmm3, xmm2 ;# fjzc fjzd
    
    movlps [rdi + rax*4], xmm0
    movhps [rdi + rbx*4], xmm0
    movlps [rdi + rcx*4], xmm1
    movhps [rdi + rdx*4], xmm1
    movss  [rdi + rax*4 + 8], xmm2
    movss  [rdi + rcx*4 + 8], xmm3
    shufps xmm2, xmm2, 1
    shufps xmm3, xmm3, 1
    movss  [rdi + rbx*4 + 8], xmm2
    movss  [rdi + rdx*4 + 8], xmm3

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb133_innerk],  4
	jl    .nb133_odd_inner
	jmp   .nb133_unroll_loop
.nb133_odd_inner:	
	add dword ptr [rsp + nb133_innerk],  4
	jnz   .nb133_odd_loop
	jmp   .nb133_updateouterdata
.nb133_odd_loop:
	mov   rdx, [rsp + nb133_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb133_innerjjnr],  4

 	xorps xmm4, xmm4  	;# clear reg.
	movss xmm4, [rsp + nb133_iqM]
	mov rsi, [rbp + nb133_charge] 
	movhps xmm4, [rsp + nb133_iqH]  ;# [qM  0  qH  qH] 
	shufps xmm4, xmm4, 41	;# [0 qH qH qM]

	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb133_qqM], xmm3	;# use dummy qq for storage 
	
	xorps xmm6, xmm6
	mov rsi, [rbp + nb133_type]
	mov ebx, [rsi + rax*4]
	mov rsi, [rbp + nb133_vdwparam]
	shl ebx, 1	
	add ebx, [rsp + nb133_ntia]
	movlps xmm6, [rsi + rbx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 252  ;# constant 11111100
	shufps xmm7, xmm7, 253  ;# constant 11111101
	movaps [rsp + nb133_c6], xmm6
	movaps [rsp + nb133_c12], xmm7
	
	mov rsi, [rbp + nb133_pos]
	lea rax, [rax + rax*2]  

	movss xmm0, [rsp + nb133_ixO]
	movss xmm1, [rsp + nb133_iyO]
	movss xmm2, [rsp + nb133_izO]
	movss xmm3, [rsp + nb133_ixH1]
	movss xmm4, [rsp + nb133_iyH1]
	movss xmm5, [rsp + nb133_izH1]
	unpcklps xmm0, [rsp + nb133_ixH2] 	;# ixO ixH2 - -
	unpcklps xmm1, [rsp + nb133_iyH2]  	;# iyO iyH2 - -
	unpcklps xmm2, [rsp + nb133_izH2]	;# izO izH2 - -
	unpcklps xmm3, [rsp + nb133_ixM] 	;# ixH1 ixM - -
	unpcklps xmm4, [rsp + nb133_iyM]  	;# iyH1 iyM - -
	unpcklps xmm5, [rsp + nb133_izM]	;# izH1 izM - -
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
	movaps [rsp + nb133_dxO], xmm3
	movaps [rsp + nb133_dyO], xmm4
	movaps [rsp + nb133_dzO], xmm5

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
	movaps xmm1, [rsp + nb133_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb133_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv, xmm4=rsq

	;# LJ table interaction
	mulps xmm4, xmm0
	mulps  xmm4, [rsp + nb133_tsc] ;# rtab
	
	cvttps2pi mm6, xmm4
	cvtpi2ps xmm6, mm6
	subss  xmm4, xmm6	
	movss xmm1, xmm4	;# xmm1=eps 
	movss xmm2, xmm1	
	mulss  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3

	movd mm0, eax	

	mov  rsi, [rbp + nb133_VFtab]
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

	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm7, [rsp + nb133_two]	;# two*Heps2 
	addss  xmm7, xmm6
	addss  xmm7, xmm5 ;# xmm7=FF 
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 

	movss xmm4, [rsp + nb133_c6]
	mulss  xmm7, xmm4	 ;# fijD 
	mulss  xmm5, xmm4	 ;# Vvdw6 
	mulss  xmm7, [rsp + nb133_tsc]
	;# put scalar force on stack Update Vvdwtot directly 
	addss  xmm5, [rsp + nb133_Vvdwtot]
	movss [rsp + nb133_fstmp], xmm7
	movss [rsp + nb133_Vvdwtot], xmm5

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
	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm7, [rsp + nb133_two]	;# two*Heps2 
	addss  xmm7, xmm6
	addss  xmm7, xmm5 ;# xmm7=FF 
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 
 	
	movss xmm4, [rsp + nb133_c12]
	mulss  xmm7, xmm4 ;# fijR 
	mulss  xmm5, xmm4 ;# Vvdw12 
	mulss  xmm7, [rsp + nb133_tsc]
	addss  xmm7, [rsp + nb133_fstmp]
	movss [rsp + nb133_fstmp], xmm7
	addss  xmm5, [rsp + nb133_Vvdwtot]
	movss [rsp + nb133_Vvdwtot], xmm5	

	movd eax, mm0	

	movaps xmm4, xmm0
	mulps  xmm4, [rsp + nb133_qqM] 
	movaps xmm2, xmm4
	mulps  xmm4, xmm0
	subss  xmm4, [rsp + nb133_fstmp]
	mulps  xmm4, xmm0
	
	addps  xmm2, [rsp + nb133_vctot]	
	movaps [rsp + nb133_vctot], xmm2
		
	movaps xmm0, [rsp + nb133_dxO]
	movaps xmm1, [rsp + nb133_dyO]
	movaps xmm2, [rsp + nb133_dzO]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4 ;# xmm0-xmm2 now contains tx-tz (partial force)
	
	movss  xmm3, [rsp + nb133_fixO]	
	movss  xmm4, [rsp + nb133_fiyO]	
	movss  xmm5, [rsp + nb133_fizO]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [rsp + nb133_fixO], xmm3	
	movss  [rsp + nb133_fiyO], xmm4	
	movss  [rsp + nb133_fizO], xmm5	;# updated the O force now do the H's
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2      
	shufps xmm3, xmm3, 0x39	;# shift right 
	shufps xmm4, xmm4, 0x39
	shufps xmm5, xmm5, 0x39
	addss  xmm3, [rsp + nb133_fixH1]
	addss  xmm4, [rsp + nb133_fiyH1]
	addss  xmm5, [rsp + nb133_fizH1]
	movss  [rsp + nb133_fixH1], xmm3	
	movss  [rsp + nb133_fiyH1], xmm4	
	movss  [rsp + nb133_fizH1], xmm5	;# updated the H1 force 

	shufps xmm3, xmm3, 0x39
	shufps xmm4, xmm4, 0x39
	shufps xmm5, xmm5, 0x39
	addss  xmm3, [rsp + nb133_fixH2]
	addss  xmm4, [rsp + nb133_fiyH2]
	addss  xmm5, [rsp + nb133_fizH2]
	movss  [rsp + nb133_fixH2], xmm3	
	movss  [rsp + nb133_fiyH2], xmm4	
	movss  [rsp + nb133_fizH2], xmm5	;# updated the H2 force 

	mov rdi, [rbp + nb133_faction]
	shufps xmm3, xmm3, 0x39
	shufps xmm4, xmm4, 0x39
	shufps xmm5, xmm5, 0x39
	addss  xmm3, [rsp + nb133_fixM]
	addss  xmm4, [rsp + nb133_fiyM]
	addss  xmm5, [rsp + nb133_fizM]
	movss  [rsp + nb133_fixM], xmm3	
	movss  [rsp + nb133_fiyM], xmm4	
	movss  [rsp + nb133_fizM], xmm5	;# updated the M force 

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

	dec dword ptr [rsp + nb133_innerk]
	jz    .nb133_updateouterdata
	jmp   .nb133_odd_loop
.nb133_updateouterdata:
	mov   ecx, [rsp + nb133_ii3]
	mov   rdi, [rbp + nb133_faction]
	mov   rsi, [rbp + nb133_fshift]
	mov   edx, [rsp + nb133_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb133_fixO]
	movaps xmm1, [rsp + nb133_fiyO]
	movaps xmm2, [rsp + nb133_fizO]

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
	movaps xmm0, [rsp + nb133_fixH1]
	movaps xmm1, [rsp + nb133_fiyH1]
	movaps xmm2, [rsp + nb133_fizH1]

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
	movaps xmm0, [rsp + nb133_fixH2]
	movaps xmm1, [rsp + nb133_fiyH2]
	movaps xmm2, [rsp + nb133_fizH2]

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
	movaps xmm0, [rsp + nb133_fixM]
	movaps xmm1, [rsp + nb133_fiyM]
	movaps xmm2, [rsp + nb133_fizM]

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
	mov esi, [rsp + nb133_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb133_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb133_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb133_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb133_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb133_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb133_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb133_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb133_n], esi
        jmp .nb133_outer
.nb133_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb133_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb133_end
        ;# non-zero, do one more workunit
        jmp   .nb133_threadloop
.nb133_end:
	mov eax, [rsp + nb133_nouter]
	mov ebx, [rsp + nb133_ninner]
	mov rcx, [rbp + nb133_outeriter]
	mov rdx, [rbp + nb133_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 1072
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





	
.globl nb_kernel133nf_x86_64_sse
.globl _nb_kernel133nf_x86_64_sse
nb_kernel133nf_x86_64_sse:	
_nb_kernel133nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb133nf_fshift,         16
.equiv          nb133nf_gid,            24
.equiv          nb133nf_pos,            32
.equiv          nb133nf_faction,        40
.equiv          nb133nf_charge,         48
.equiv          nb133nf_p_facel,        56
.equiv          nb133nf_argkrf,         64
.equiv          nb133nf_argcrf,         72
.equiv          nb133nf_Vc,             80
.equiv          nb133nf_type,           88
.equiv          nb133nf_p_ntype,        96
.equiv          nb133nf_vdwparam,       104
.equiv          nb133nf_Vvdw,           112
.equiv          nb133nf_p_tabscale,     120
.equiv          nb133nf_VFtab,          128
.equiv          nb133nf_invsqrta,       136
.equiv          nb133nf_dvda,           144
.equiv          nb133nf_p_gbtabscale,   152
.equiv          nb133nf_GBtab,          160
.equiv          nb133nf_p_nthreads,     168
.equiv          nb133nf_count,          176
.equiv          nb133nf_mtx,            184
.equiv          nb133nf_outeriter,      192
.equiv          nb133nf_inneriter,      200
.equiv          nb133nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb133nf_ixO,            0
.equiv          nb133nf_iyO,            16
.equiv          nb133nf_izO,            32
.equiv          nb133nf_ixH1,           48
.equiv          nb133nf_iyH1,           64
.equiv          nb133nf_izH1,           80
.equiv          nb133nf_ixH2,           96
.equiv          nb133nf_iyH2,           112
.equiv          nb133nf_izH2,           128
.equiv          nb133nf_ixM,            144
.equiv          nb133nf_iyM,            160
.equiv          nb133nf_izM,            176
.equiv          nb133nf_iqM,            192
.equiv          nb133nf_iqH,            208
.equiv          nb133nf_qqM,            224
.equiv          nb133nf_qqH,            240
.equiv          nb133nf_rinvH1,         256
.equiv          nb133nf_rinvH2,         272
.equiv          nb133nf_rinvM,          288
.equiv          nb133nf_tsc,            304
.equiv          nb133nf_c6,             320
.equiv          nb133nf_c12,            336
.equiv          nb133nf_krf,            352
.equiv          nb133nf_crf,            368
.equiv          nb133nf_krsqH1,         384
.equiv          nb133nf_krsqH2,         400
.equiv          nb133nf_krsqM,          416
.equiv          nb133nf_vctot,          432
.equiv          nb133nf_Vvdwtot,        448
.equiv          nb133nf_half,           464
.equiv          nb133nf_three,          480
.equiv          nb133nf_nri,            496
.equiv          nb133nf_iinr,           504
.equiv          nb133nf_jindex,         512
.equiv          nb133nf_jjnr,           520
.equiv          nb133nf_shift,          528
.equiv          nb133nf_shiftvec,       536
.equiv          nb133nf_facel,          544
.equiv          nb133nf_innerjjnr,      552
.equiv          nb133nf_is3,            560
.equiv          nb133nf_ii3,            564
.equiv          nb133nf_ntia,           568
.equiv          nb133nf_innerk,         572
.equiv          nb133nf_n,              576
.equiv          nb133nf_nn1,            580
.equiv          nb133nf_nouter,         584
.equiv          nb133nf_ninner,         588

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
	mov [rsp + nb133nf_nouter], eax
	mov [rsp + nb133nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb133nf_nri], edi
	mov [rsp + nb133nf_iinr], rsi
	mov [rsp + nb133nf_jindex], rdx
	mov [rsp + nb133nf_jjnr], rcx
	mov [rsp + nb133nf_shift], r8
	mov [rsp + nb133nf_shiftvec], r9
	mov rsi, [rbp + nb133nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb133nf_facel], xmm0

	mov rax, [rbp + nb133nf_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb133nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb133nf_half], eax
	movss xmm1, [rsp + nb133nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb133nf_half],  xmm1
	movaps [rsp + nb133nf_three],  xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb133nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb133nf_charge]
	movss xmm4, [rdx + rbx*4 + 4]	
	movss xmm3, [rdx + rbx*4 + 12]	
	mov rsi, [rbp + nb133nf_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb133nf_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb133nf_iqM], xmm3
	movaps [rsp + nb133nf_iqH], xmm4
	
	mov   rdx, [rbp + nb133nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb133nf_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb133nf_ntia], ecx		

.nb133nf_threadloop:
        mov   rsi, [rbp + nb133nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb133nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb133nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb133nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb133nf_n], eax
        mov [rsp + nb133nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb133nf_outerstart
        jmp .nb133nf_end
	
.nb133nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb133nf_nouter]
	mov [rsp + nb133nf_nouter], ebx

.nb133nf_outer:
	mov   rax, [rsp + nb133nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb133nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb133nf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb133nf_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2	
	movaps xmm6, xmm0
	movaps xmm7, xmm1
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb133nf_pos]	;# eax = base of pos[]  
	mov   [rsp + nb133nf_ii3], ebx

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
	movaps [rsp + nb133nf_ixO], xmm3
	movaps [rsp + nb133nf_iyO], xmm4
	movaps [rsp + nb133nf_izO], xmm5
	movaps [rsp + nb133nf_ixH1], xmm6
	movaps [rsp + nb133nf_iyH1], xmm7

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
	movaps [rsp + nb133nf_izH1], xmm6
	movaps [rsp + nb133nf_ixH2], xmm0
	movaps [rsp + nb133nf_iyH2], xmm1
	movaps [rsp + nb133nf_izH2], xmm2
	movaps [rsp + nb133nf_ixM], xmm3
	movaps [rsp + nb133nf_iyM], xmm4
	movaps [rsp + nb133nf_izM], xmm5
	
	;# clear vctot 
	xorps xmm4, xmm4
	movaps [rsp + nb133nf_vctot], xmm4
	movaps [rsp + nb133nf_Vvdwtot], xmm4

	mov   rax, [rsp + nb133nf_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb133nf_pos]
	mov   rax, [rsp + nb133nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb133nf_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb133nf_ninner]
	mov   [rsp + nb133nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb133nf_innerk], edx	;# number of innerloop atoms 
	jge   .nb133nf_unroll_loop
	jmp   .nb133nf_odd_inner
.nb133nf_unroll_loop:
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb133nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb133nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb133nf_charge]	;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	 	;# and in xmm4 
	mulps  xmm3, [rsp + nb133nf_iqM]
	mulps  xmm4, [rsp + nb133nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps  [rsp + nb133nf_qqM], xmm3
	movaps  [rsp + nb133nf_qqH], xmm4
	
	mov rsi, [rbp + nb133nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]
	mov rsi, [rbp + nb133nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb133nf_ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm6, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm6, xmm7, 221  ;# constant 11011101
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [rsp + nb133nf_c6], xmm4
	movaps [rsp + nb133nf_c12], xmm6

	mov rsi, [rbp + nb133nf_pos]   	;# base of pos[] 

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

	shufps xmm2, xmm6, 136  ;# constant 10001000
	
	shufps xmm0, xmm5, 136  ;# constant 10001000
	shufps xmm1, xmm5, 221  ;# constant 11011101		

	;# move ixO-izO to xmm4-xmm6 
	movaps xmm4, [rsp + nb133nf_ixO]
	movaps xmm5, [rsp + nb133nf_iyO]
	movaps xmm6, [rsp + nb133nf_izO]

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
	movaps xmm4, [rsp + nb133nf_ixH1]
	movaps xmm5, [rsp + nb133nf_iyH1]
	movaps xmm6, [rsp + nb133nf_izH1]

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
	movaps xmm3, [rsp + nb133nf_ixH2]
	movaps xmm4, [rsp + nb133nf_iyH2]
	movaps xmm5, [rsp + nb133nf_izH2]

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
	movaps xmm3, [rsp + nb133nf_iyM]
	movaps xmm4, [rsp + nb133nf_izM]
	subps  xmm3, xmm1
	subps  xmm4, xmm2
	movaps xmm2, [rsp + nb133nf_ixM]
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
	movaps  xmm0, [rsp + nb133nf_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm0, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm0, [rsp + nb133nf_half]
	movaps  [rsp + nb133nf_rinvH1], xmm0	;# rinvH1 

	;# rsqH2 - seed to xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm0, [rsp + nb133nf_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm0, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm0, [rsp + nb133nf_half]
	movaps  [rsp + nb133nf_rinvH2], xmm0	;# rinvH2 

	;# rsqM - seed to xmm2 
	rsqrtps xmm2, xmm4
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm0, [rsp + nb133nf_three]
	mulps   xmm2, xmm4	;# rsq*lu*lu 
	subps   xmm0, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm0, [rsp + nb133nf_half]
	movaps  [rsp + nb133nf_rinvM], xmm0
	
	;# Do the O LJ-only interaction directly.	
	;# rsqO is in xmm7
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb133nf_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb133nf_half]
	movaps  xmm0, xmm4
	;# xmm0=rinvO
	
	mulps xmm7, xmm0
	mulps xmm7, [rsp + nb133nf_tsc] ;# rtab
	
	movhlps xmm5, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm5	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps  xmm7, xmm6	
	movaps xmm1, xmm7	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3
	pslld mm7, 3

	mov  rsi, [rbp + nb133nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7
	
	;# dispersion 
	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half dispersion table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101
	
	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of dispersion table 
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# dispersion table ready, in xmm4-xmm7 	

	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 

	movaps xmm4, [rsp + nb133nf_c6]
	mulps  xmm5, xmm4	 ;# Vvdw6 

	addps  xmm5, [rsp + nb133nf_Vvdwtot]
	movaps [rsp + nb133nf_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm5, [rsi + rax*4 + 16]
	movlps xmm7, [rsi + rcx*4 + 16]
	movhps xmm5, [rsi + rbx*4 + 16]
	movhps xmm7, [rsi + rdx*4 + 16] ;# got half repulsion table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [rsi + rax*4 + 24]
	movlps xmm3, [rsi + rcx*4 + 24]
	movhps xmm7, [rsi + rbx*4 + 24]
	movhps xmm3, [rsi + rdx*4 + 24] ;# other half of repulsion table 
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
 	
	movaps xmm4, [rsp + nb133nf_c12]
	mulps  xmm5, xmm4 ;# Vvdw12 

	addps  xmm5, [rsp + nb133nf_Vvdwtot]
	movaps [rsp + nb133nf_Vvdwtot], xmm5

	;# Do H1-H2-M interactions	
	movaps  xmm7, [rsp + nb133nf_rinvH1]
	addps   xmm7, [rsp + nb133nf_rinvH2]
	movaps  xmm6, [rsp + nb133nf_rinvM]
	
	mulps   xmm7, [rsp + nb133nf_qqH]
	mulps   xmm6, [rsp + nb133nf_qqM]
	addps   xmm7, xmm6

	addps  xmm7, [rsp + nb133nf_vctot]	
	movaps [rsp + nb133nf_vctot], xmm7

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb133nf_innerk],  4
	jl    .nb133nf_odd_inner
	jmp   .nb133nf_unroll_loop
.nb133nf_odd_inner:	
	add dword ptr [rsp + nb133nf_innerk],  4
	jnz   .nb133nf_odd_loop
	jmp   .nb133nf_updateouterdata
.nb133nf_odd_loop:
	mov   rdx, [rsp + nb133nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb133nf_innerjjnr],  4

 	xorps xmm4, xmm4  	;# clear reg.
	movss xmm4, [rsp + nb133nf_iqM]
	mov rsi, [rbp + nb133nf_charge] 
	movhps xmm4, [rsp + nb133nf_iqH]  ;# [qM  0  qH  qH] 
	shufps xmm4, xmm4, 41	;# [0 qH qH qM]

	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb133nf_qqM], xmm3	;# use dummy qq for storage 
	
	xorps xmm6, xmm6
	mov rsi, [rbp + nb133nf_type]
	mov ebx, [rsi + rax*4]
	mov rsi, [rbp + nb133nf_vdwparam]
	shl ebx, 1	
	add ebx, [rsp + nb133nf_ntia]
	movlps xmm6, [rsi + rbx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 252  ;# constant 11111100
	shufps xmm7, xmm7, 253  ;# constant 11111101
	movaps [rsp + nb133nf_c6], xmm6
	movaps [rsp + nb133nf_c12], xmm7
	
	mov rsi, [rbp + nb133nf_pos]
	lea rax, [rax + rax*2]  

	movss xmm3, [rsp + nb133nf_ixO]
	movss xmm4, [rsp + nb133nf_iyO]
	movss xmm5, [rsp + nb133nf_izO]
	movss xmm0, [rsp + nb133nf_ixH1]
	movss xmm1, [rsp + nb133nf_iyH1]
	movss xmm2, [rsp + nb133nf_izH1]
	unpcklps xmm3, [rsp + nb133nf_ixH2] 	;# ixO ixH2 - -
	unpcklps xmm4, [rsp + nb133nf_iyH2]  	;# iyO iyH2 - -
	unpcklps xmm5, [rsp + nb133nf_izH2]	;# izO izH2 - -
	unpcklps xmm0, [rsp + nb133nf_ixM] 	;# ixH1 ixM - -
	unpcklps xmm1, [rsp + nb133nf_iyM]  	;# iyH1 iyM - -
	unpcklps xmm2, [rsp + nb133nf_izM]	;# izH1 izM - -
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
	movaps xmm1, [rsp + nb133nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb133nf_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv, xmm4=rsq

	;# LJ table interaction
	mulps xmm4, xmm0
	mulps  xmm4, [rsp + nb133nf_tsc] ;# rtab
	
	cvttps2pi mm6, xmm4
	cvtpi2ps xmm6, mm6
	subss  xmm4, xmm6	
	movss xmm1, xmm4	;# xmm1=eps 
	movss xmm2, xmm1	
	mulss  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3

	mov  rsi, [rbp + nb133nf_VFtab]
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

	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 

	movss xmm4, [rsp + nb133nf_c6]
	mulss  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack Update Vvdwtot directly 
	addss  xmm5, [rsp + nb133nf_Vvdwtot]
	movss [rsp + nb133nf_Vvdwtot], xmm5

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
	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 
 	
	movss xmm4, [rsp + nb133nf_c12]
	mulss  xmm5, xmm4 ;# Vvdw12 

	addss  xmm5, [rsp + nb133nf_Vvdwtot]
	movss [rsp + nb133nf_Vvdwtot], xmm5	

	mulps  xmm0, [rsp + nb133nf_qqM]	;# xmm0=vcoul 
	
	addps  xmm0, [rsp + nb133nf_vctot]	
	movaps [rsp + nb133nf_vctot], xmm0
	
	dec dword ptr [rsp + nb133nf_innerk]
	jz    .nb133nf_updateouterdata
	jmp   .nb133nf_odd_loop
.nb133nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb133nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb133nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb133nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb133nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb133nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb133nf_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb133nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb133nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb133nf_n], esi
        jmp .nb133nf_outer
.nb133nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb133nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb133nf_end
        ;# non-zero, do one more workunit
        jmp   .nb133nf_threadloop
.nb133nf_end:

	mov eax, [rsp + nb133nf_nouter]
	mov ebx, [rsp + nb133nf_ninner]
	mov rcx, [rbp + nb133nf_outeriter]
	mov rdx, [rbp + nb133nf_inneriter]
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
