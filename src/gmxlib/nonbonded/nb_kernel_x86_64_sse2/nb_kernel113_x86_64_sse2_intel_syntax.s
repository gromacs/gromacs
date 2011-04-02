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

.globl nb_kernel113_x86_64_sse2
.globl _nb_kernel113_x86_64_sse2
nb_kernel113_x86_64_sse2:	
_nb_kernel113_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb113_fshift,           16
.equiv          nb113_gid,              24
.equiv          nb113_pos,              32
.equiv          nb113_faction,          40
.equiv          nb113_charge,           48
.equiv          nb113_p_facel,          56
.equiv          nb113_argkrf,           64
.equiv          nb113_argcrf,           72
.equiv          nb113_Vc,               80
.equiv          nb113_type,             88
.equiv          nb113_p_ntype,          96
.equiv          nb113_vdwparam,         104
.equiv          nb113_Vvdw,             112
.equiv          nb113_p_tabscale,       120
.equiv          nb113_VFtab,            128
.equiv          nb113_invsqrta,         136
.equiv          nb113_dvda,             144
.equiv          nb113_p_gbtabscale,     152
.equiv          nb113_GBtab,            160
.equiv          nb113_p_nthreads,       168
.equiv          nb113_count,            176
.equiv          nb113_mtx,              184
.equiv          nb113_outeriter,        192
.equiv          nb113_inneriter,        200
.equiv          nb113_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb113_ixO,              0
.equiv          nb113_iyO,              16
.equiv          nb113_izO,              32
.equiv          nb113_ixH1,             48
.equiv          nb113_iyH1,             64
.equiv          nb113_izH1,             80
.equiv          nb113_ixH2,             96
.equiv          nb113_iyH2,             112
.equiv          nb113_izH2,             128
.equiv          nb113_ixM,              144
.equiv          nb113_iyM,              160
.equiv          nb113_izM,              176
.equiv          nb113_iqH,              192
.equiv          nb113_iqM,              208
.equiv          nb113_dxO,              224
.equiv          nb113_dyO,              240
.equiv          nb113_dzO,              256
.equiv          nb113_dxH1,             272
.equiv          nb113_dyH1,             288
.equiv          nb113_dzH1,             304
.equiv          nb113_dxH2,             320
.equiv          nb113_dyH2,             336
.equiv          nb113_dzH2,             352
.equiv          nb113_dxM,              368
.equiv          nb113_dyM,              384
.equiv          nb113_dzM,              400
.equiv          nb113_qqH,              416
.equiv          nb113_qqM,              432
.equiv          nb113_c6,               448
.equiv          nb113_c12,              464
.equiv          nb113_six,              480
.equiv          nb113_twelve,           496
.equiv          nb113_vctot,            512
.equiv          nb113_Vvdwtot,          528
.equiv          nb113_fixO,             544
.equiv          nb113_fiyO,             560
.equiv          nb113_fizO,             576
.equiv          nb113_fixH1,            592
.equiv          nb113_fiyH1,            608
.equiv          nb113_fizH1,            624
.equiv          nb113_fixH2,            640
.equiv          nb113_fiyH2,            656
.equiv          nb113_fizH2,            672
.equiv          nb113_fixM,             688
.equiv          nb113_fiyM,             704
.equiv          nb113_fizM,             720
.equiv          nb113_fjx,              736
.equiv          nb113_fjy,              752
.equiv          nb113_fjz,              768
.equiv          nb113_half,             784
.equiv          nb113_three,            800
.equiv          nb113_two,              816
.equiv          nb113_rinvH1,           832
.equiv          nb113_rinvH2,           848
.equiv          nb113_rinvM,            864
.equiv          nb113_nri,              880
.equiv          nb113_iinr,             888
.equiv          nb113_jindex,           896
.equiv          nb113_jjnr,             904
.equiv          nb113_shift,            912
.equiv          nb113_shiftvec,         920
.equiv          nb113_facel,            928
.equiv          nb113_innerjjnr,        936
.equiv          nb113_is3,              944
.equiv          nb113_ii3,              948
.equiv          nb113_ntia,             952
.equiv          nb113_innerk,           956
.equiv          nb113_n,                960
.equiv          nb113_nn1,              964
.equiv          nb113_nouter,           968
.equiv          nb113_ninner,           972

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
	sub rsp, 976		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb113_nouter], eax
	mov [rsp + nb113_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb113_nri], edi
	mov [rsp + nb113_iinr], rsi
	mov [rsp + nb113_jindex], rdx
	mov [rsp + nb113_jjnr], rcx
	mov [rsp + nb113_shift], r8
	mov [rsp + nb113_shiftvec], r9
	mov rsi, [rbp + nb113_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb113_facel], xmm0

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000     ;# upper half of half
	mov [rsp + nb113_half], eax
	mov [rsp + nb113_half+4], ebx
	movsd xmm1, [rsp + nb113_half]
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
	movapd [rsp + nb113_half], xmm1
	movapd [rsp + nb113_two], xmm2
	movapd [rsp + nb113_three], xmm3
	movapd [rsp + nb113_six], xmm4
	movapd [rsp + nb113_twelve], xmm5

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb113_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb113_charge]
	movsd xmm3, [rdx + rbx*8 + 8]	
	movsd xmm4, [rdx + rbx*8 + 24]	
	mov rsi, [rbp + nb113_p_facel]
	movsd xmm0, [rsi]
	movsd xmm5, [rsp + nb113_facel]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb113_iqH], xmm3
	movapd [rsp + nb113_iqM], xmm4
	
	mov   rdx, [rbp + nb113_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb113_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb113_ntia], ecx		
.nb113_threadloop:
        mov   rsi, [rbp + nb113_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb113_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb113_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb113_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb113_n], eax
        mov [rsp + nb113_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb113_outerstart
        jmp .nb113_end

.nb113_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb113_nouter]
	mov [rsp + nb113_nouter], ebx

.nb113_outer:
	mov   rax, [rsp + nb113_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb113_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb113_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb113_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb113_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb113_ii3], ebx

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
	movapd [rsp + nb113_ixO], xmm3
	movapd [rsp + nb113_iyO], xmm4
	movapd [rsp + nb113_izO], xmm5
	movapd [rsp + nb113_ixH1], xmm6
	movapd [rsp + nb113_iyH1], xmm7

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
	movapd [rsp + nb113_izH1], xmm6
	movapd [rsp + nb113_ixH2], xmm0
	movapd [rsp + nb113_iyH2], xmm1
	movapd [rsp + nb113_izH2], xmm2
	movapd [rsp + nb113_ixM], xmm3
	movapd [rsp + nb113_iyM], xmm4
	movapd [rsp + nb113_izM], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb113_vctot], xmm4
	movapd [rsp + nb113_Vvdwtot], xmm4
	movapd [rsp + nb113_fixO], xmm4
	movapd [rsp + nb113_fiyO], xmm4
	movapd [rsp + nb113_fizO], xmm4
	movapd [rsp + nb113_fixH1], xmm4
	movapd [rsp + nb113_fiyH1], xmm4
	movapd [rsp + nb113_fizH1], xmm4
	movapd [rsp + nb113_fixH2], xmm4
	movapd [rsp + nb113_fiyH2], xmm4
	movapd [rsp + nb113_fizH2], xmm4
	movapd [rsp + nb113_fixM], xmm4
	movapd [rsp + nb113_fiyM], xmm4
	movapd [rsp + nb113_fizM], xmm4
	
	mov   rax, [rsp + nb113_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb113_pos]
	mov   rdi, [rbp + nb113_faction]	
	mov   rax, [rsp + nb113_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb113_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb113_ninner]
	mov   [rsp + nb113_ninner], ecx
	add   edx, 0
	mov   [rsp + nb113_innerk], edx    ;# number of innerloop atoms 
	jge   .nb113_unroll_loop
	jmp   .nb113_checksingle
.nb113_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb113_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]

	add qword ptr [rsp + nb113_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb113_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [rsp + nb113_iqM]
	mulpd  xmm4, [rsp + nb113_iqH]

	movapd  [rsp + nb113_qqM], xmm3
	movapd  [rsp + nb113_qqH], xmm4
	
	mov rsi, [rbp + nb113_type]
	mov r8d, [rsi + rax*4]
	mov r9d, [rsi + rbx*4]
	mov rsi, [rbp + nb113_vdwparam]
	shl r8d, 1	
	shl r9d, 1	
	mov edi, [rsp + nb113_ntia]
	add r8d, edi
	add r9d, edi

	movlpd xmm6, [rsi + r8*8]	;# c6a c6b
	movhpd xmm6, [rsi + r9*8]	;# 
	movlpd xmm7, [rsi + r8*8 + 8]	;# c12a c12b 
	movhpd xmm7, [rsi + r9*8 + 8]	;# 
	
	movapd [rsp + nb113_c6], xmm6
	movapd [rsp + nb113_c12], xmm7
	
	mov rsi, [rbp + nb113_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# move j coordinates to local temp variables 
    movlpd xmm0, [rsi + rax*8] 
    movlpd xmm1, [rsi + rax*8 + 8] 
    movlpd xmm2, [rsi + rax*8 + 16] 
    movhpd xmm0, [rsi + rbx*8] 
    movhpd xmm1, [rsi + rbx*8 + 8] 
    movhpd xmm2, [rsi + rbx*8 + 16] 

    ;# xmm0 = jx
    ;# xmm1 = jy
    ;# xmm2 = jz
    
    ;# O interaction
    ;# copy to xmm3-xmm5
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    
    subpd xmm3, [rsp + nb113_ixO]
    subpd xmm4, [rsp + nb113_iyO]
    subpd xmm5, [rsp + nb113_izO]
    
    movapd xmm13, xmm3
    movapd xmm14, xmm4
    movapd xmm15, xmm5
    
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5

	addpd  xmm3, xmm4
	addpd  xmm3, xmm5

    ;# calc 1/rsq
    cvtpd2ps xmm6, xmm3    
    rcpps xmm6, xmm6
    cvtps2pd xmm6, xmm6     ;# lu in low xmm6 
        
    ;# 1/x lookup seed in xmm6 
    movapd xmm4, [rsp + nb113_two]
    movapd xmm5, xmm3       ;# rsq
    mulpd xmm3, xmm6        ;# lu*rsq 
    subpd xmm4, xmm3        ;# 2-lu*rsq 
    mulpd xmm6, xmm4        ;# (new lu) 
        
    movapd xmm4, [rsp + nb113_two]
    mulpd xmm5, xmm6        ;# lu*rsq 
    subpd xmm4, xmm5        ;# 2-lu*rsq 
    mulpd xmm4, xmm6        ;# xmm4=rinvsq 

    movapd xmm3, xmm4       ;# rinvsq
    mulpd  xmm4, xmm4       ;# rinv4
    mulpd  xmm4, xmm3       ;# rinv6
    movapd xmm5, xmm4     
    mulpd  xmm5, xmm5       ;# rinv12
    mulpd  xmm4, [rsp + nb113_c6]
    mulpd  xmm5, [rsp + nb113_c12]
    movapd xmm6, xmm5
    subpd  xmm6, xmm4  ;# Vvdw=vvdw12-vvdw6
    mulpd  xmm4, [rsp + nb113_six]
    mulpd  xmm5, [rsp + nb113_twelve]
    subpd  xmm5, xmm4
    mulpd  xmm3, xmm5   ;# fscal
    
    addpd  xmm6, [rsp + nb113_Vvdwtot]
    movapd [rsp + nb113_Vvdwtot], xmm6
    
    mulpd  xmm13, xmm3 ;# fx
    mulpd  xmm14, xmm3 ;# fy
    mulpd  xmm15, xmm3 ;# fz

    ;# save j force temporarily
    movapd [rsp + nb113_fjx], xmm13
    movapd [rsp + nb113_fjy], xmm14
    movapd [rsp + nb113_fjz], xmm15
    
    ;# increment i O force
    addpd xmm13, [rsp + nb113_fixO]
    addpd xmm14, [rsp + nb113_fiyO]
    addpd xmm15, [rsp + nb113_fizO]
    movapd [rsp + nb113_fixO], xmm13
    movapd [rsp + nb113_fiyO], xmm14
    movapd [rsp + nb113_fizO], xmm15    
    ;# finished O LJ interaction.


    ;# do H1, H2, and M interactions in parallel.
    ;# xmm0-xmm2 still contain j coordinates.        
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subpd xmm0, [rsp + nb113_ixH1]
    subpd xmm1, [rsp + nb113_iyH1]
    subpd xmm2, [rsp + nb113_izH1]
    subpd xmm3, [rsp + nb113_ixH2]
    subpd xmm4, [rsp + nb113_iyH2]
    subpd xmm5, [rsp + nb113_izH2]
    subpd xmm6, [rsp + nb113_ixM]
    subpd xmm7, [rsp + nb113_iyM]
    subpd xmm8, [rsp + nb113_izM]
    
	movapd [rsp + nb113_dxH1], xmm0
	movapd [rsp + nb113_dyH1], xmm1
	movapd [rsp + nb113_dzH1], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb113_dxH2], xmm3
	movapd [rsp + nb113_dyH2], xmm4
	movapd [rsp + nb113_dzH2], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb113_dxM], xmm6
	movapd [rsp + nb113_dyM], xmm7
	movapd [rsp + nb113_dzM], xmm8
	mulpd  xmm6, xmm6
	mulpd  xmm7, xmm7
	mulpd  xmm8, xmm8
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
    addpd  xmm6, xmm7
    addpd  xmm6, xmm8

	;# start doing invsqrt for j atoms
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
		
	movapd  xmm9, [rsp + nb113_three]
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

	movapd  xmm15, [rsp + nb113_half]
	mulpd   xmm9, xmm15  ;# first iteration for rinvH1
	mulpd   xmm10, xmm15 ;# first iteration for rinvH2
    mulpd   xmm11, xmm15 ;# first iteration for rinvM

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulpd   xmm2, xmm2 ;# lu*lu
	mulpd   xmm5, xmm5 ;# lu*lu
    mulpd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb113_three]
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

	movapd  xmm15, [rsp + nb113_half]
	mulpd   xmm9, xmm15  ;#  rinvH1
	mulpd   xmm10, xmm15 ;#   rinvH2
    mulpd   xmm11, xmm15 ;#   rinvM
	
	;# interactions 
    movapd xmm0, xmm9
    movapd xmm1, xmm10
    movapd xmm2, xmm11
    mulpd  xmm9, xmm9
    mulpd  xmm10, xmm10
    mulpd  xmm11, xmm11
    mulpd  xmm0, [rsp + nb113_qqH] 
    mulpd  xmm1, [rsp + nb113_qqH] 
    mulpd  xmm2, [rsp + nb113_qqM] 
    mulpd  xmm9, xmm0
    mulpd  xmm10, xmm1
    mulpd  xmm11, xmm2
    
    addpd xmm0, [rsp + nb113_vctot] 
    addpd xmm1, xmm2
    addpd xmm0, xmm1
    movapd [rsp + nb113_vctot], xmm0
    
    ;# move j forces to xmm0-xmm2
    mov rdi, [rbp + nb113_faction]
	movlpd xmm0, [rdi + rax*8]
	movlpd xmm1, [rdi + rax*8 + 8]
	movlpd xmm2, [rdi + rax*8 + 16]
	movhpd xmm0, [rdi + rbx*8]
	movhpd xmm1, [rdi + rbx*8 + 8]
	movhpd xmm2, [rdi + rbx*8 + 16]

    movapd xmm7, xmm9
    movapd xmm8, xmm9
    movapd xmm13, xmm11
    movapd xmm14, xmm11
    movapd xmm15, xmm11
    movapd xmm11, xmm10
    movapd xmm12, xmm10

    ;# add forces from O interaction
    addpd xmm0, [rsp + nb113_fjx]
    addpd xmm1, [rsp + nb113_fjy]
    addpd xmm2, [rsp + nb113_fjz]

	mulpd xmm7, [rsp + nb113_dxH1]
	mulpd xmm8, [rsp + nb113_dyH1]
	mulpd xmm9, [rsp + nb113_dzH1]
	mulpd xmm10, [rsp + nb113_dxH2]
	mulpd xmm11, [rsp + nb113_dyH2]
	mulpd xmm12, [rsp + nb113_dzH2]
	mulpd xmm13, [rsp + nb113_dxM]
	mulpd xmm14, [rsp + nb113_dyM]
	mulpd xmm15, [rsp + nb113_dzM]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb113_fixH1]
    addpd xmm8, [rsp + nb113_fiyH1]
    addpd xmm9, [rsp + nb113_fizH1]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb113_fixH2]
    addpd xmm11, [rsp + nb113_fiyH2]
    addpd xmm12, [rsp + nb113_fizH2]

    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    addpd xmm13, [rsp + nb113_fixM]
    addpd xmm14, [rsp + nb113_fiyM]
    addpd xmm15, [rsp + nb113_fizM]

    movapd [rsp + nb113_fixH1], xmm7
    movapd [rsp + nb113_fiyH1], xmm8
    movapd [rsp + nb113_fizH1], xmm9
    movapd [rsp + nb113_fixH2], xmm10
    movapd [rsp + nb113_fiyH2], xmm11
    movapd [rsp + nb113_fizH2], xmm12
    movapd [rsp + nb113_fixM], xmm13
    movapd [rsp + nb113_fiyM], xmm14
    movapd [rsp + nb113_fizM], xmm15
   
    ;# store back j forces from xmm0-xmm2
	movlpd [rdi + rax*8], xmm0
	movlpd [rdi + rax*8 + 8], xmm1
	movlpd [rdi + rax*8 + 16], xmm2
	movhpd [rdi + rbx*8], xmm0
	movhpd [rdi + rbx*8 + 8], xmm1
	movhpd [rdi + rbx*8 + 16], xmm2

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb113_innerk],  2
	jl   .nb113_checksingle
	jmp  .nb113_unroll_loop
.nb113_checksingle:	
	mov   edx, [rsp + nb113_innerk]
	and   edx, 1
	jnz  .nb113_dosingle
	jmp  .nb113_updateouterdata
.nb113_dosingle:
	mov   rdx, [rsp + nb113_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb113_innerjjnr],  4	

	mov rsi, [rbp + nb113_charge]    ;# base of charge[] 

	xorpd xmm3, xmm3
	movlpd xmm3, [rsi + rax*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [rsp + nb113_iqM]
	mulpd  xmm4, [rsp + nb113_iqH]

	movapd  [rsp + nb113_qqM], xmm3
	movapd  [rsp + nb113_qqH], xmm4
	
	mov rsi, [rbp + nb113_type]
	mov r8d, [rsi + rax*4]
	mov rsi, [rbp + nb113_vdwparam]
	shl r8d, 1	
	mov edi, [rsp + nb113_ntia]
	add r8d, edi

	movsd xmm6, [rsi + r8*8]	;# c6a
    movsd xmm7, [rsi + r8*8 + 8]	;# c12a 	
	movapd [rsp + nb113_c6], xmm6
	movapd [rsp + nb113_c12], xmm7
	
	mov rsi, [rbp + nb113_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2  and xmm4-xmm6
	movlpd xmm4, [rsi + rax*8]
	movlpd xmm5, [rsi + rax*8 + 8]
	movlpd xmm6, [rsi + rax*8 + 16]
    movapd xmm0, xmm4
    movapd xmm1, xmm5
    movapd xmm2, xmm6

	;# calc dr 
	subsd xmm4, [rsp + nb113_ixO]
	subsd xmm5, [rsp + nb113_iyO]
	subsd xmm6, [rsp + nb113_izO]

	;# store dr 
	movapd [rsp + nb113_dxO], xmm4
	movapd [rsp + nb113_dyO], xmm5
	movapd [rsp + nb113_dzO], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 

	;# move j coords to xmm4-xmm6 
	movapd xmm4, xmm0
	movapd xmm5, xmm1
	movapd xmm6, xmm2

	;# calc dr 
	subsd xmm4, [rsp + nb113_ixH1]
	subsd xmm5, [rsp + nb113_iyH1]
	subsd xmm6, [rsp + nb113_izH1]

	;# store dr 
	movapd [rsp + nb113_dxH1], xmm4
	movapd [rsp + nb113_dyH1], xmm5
	movapd [rsp + nb113_dzH1], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm6, xmm5
	addsd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move j coords to xmm3-xmm5
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	;# calc dr 
	subsd xmm3, [rsp + nb113_ixH2]
	subsd xmm4, [rsp + nb113_iyH2]
	subsd xmm5, [rsp + nb113_izH2]


	;# store dr 
	movapd [rsp + nb113_dxH2], xmm3
	movapd [rsp + nb113_dyH2], xmm4
	movapd [rsp + nb113_dzH2], xmm5
	;# square it 
	mulsd xmm3,xmm3
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	addsd xmm5, xmm4
	addsd xmm5, xmm3

	;# move j coords to xmm4-xmm2
	movapd xmm4, xmm0
	movapd xmm3, xmm1
    ;# xmm2 already contains z

	;# calc dr 
	subsd xmm4, [rsp + nb113_ixM]
	subsd xmm3, [rsp + nb113_iyM]
	subsd xmm2, [rsp + nb113_izM]

	;# store dr 
	movapd [rsp + nb113_dxM], xmm4
	movapd [rsp + nb113_dyM], xmm3
	movapd [rsp + nb113_dzM], xmm2
    
	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# start with rsqH1 - put seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [rsp + nb113_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb113_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb113_three]
	subsd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb113_half] ;# rinv 
	movapd [rsp + nb113_rinvH1], xmm1
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [rsp + nb113_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb113_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb113_three]
	subsd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb113_half] ;# rinv 
	movapd [rsp + nb113_rinvH2], xmm1
	
	;# rsqM - seed in xmm2 
	cvtsd2ss xmm2, xmm4
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [rsp + nb113_three]
	mulsd   xmm2, xmm4	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb113_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb113_three]
	subsd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb113_half] ;# rinv 
	movapd [rsp + nb113_rinvM], xmm1

	;# do O interactions directly. xmm7=rsq
	cvtsd2ss xmm2, xmm7
	movapd   xmm6, xmm7
	rcpps    xmm2, xmm2
	cvtss2sd xmm2, xmm2
	movapd   xmm1, [rsp + nb113_two]
	movapd   xmm0, xmm1
	mulsd   xmm7, xmm2
	subsd   xmm1, xmm7
	mulsd   xmm2, xmm1 ;# iter1 
	mulsd   xmm6, xmm2
	subsd   xmm0, xmm6
	mulsd   xmm0, xmm2 ;# xmm0=rinvsq
	movapd  xmm1, xmm0	
	mulsd   xmm1, xmm1 ;# rinv4
	mulsd   xmm1, xmm0 ;#rinvsix
	movapd  xmm2, xmm1
	mulsd	xmm2, xmm2 ;# rinvtwelve
	mulsd  xmm1, [rsp + nb113_c6]
	mulsd  xmm2, [rsp + nb113_c12]
	movapd xmm3, xmm2
	subsd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addsd  xmm3, [rsp + nb113_Vvdwtot]
	mulsd  xmm1, [rsp + nb113_six]
	mulsd  xmm2, [rsp + nb113_twelve]
	subsd  xmm2, xmm1
	mulsd  xmm2, xmm0
	movapd xmm4, xmm2 ;# total fsO 
	movsd [rsp + nb113_Vvdwtot], xmm3

	movapd xmm0, [rsp + nb113_dxO]
	movapd xmm1, [rsp + nb113_dyO]
	movapd xmm2, [rsp + nb113_dzO]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update O forces 
	movapd xmm3, [rsp + nb113_fixO]
	movapd xmm4, [rsp + nb113_fiyO]
	movapd xmm7, [rsp + nb113_fizO]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [rsp + nb113_fixO], xmm3
	movsd [rsp + nb113_fiyO], xmm4
	movsd [rsp + nb113_fizO], xmm7
	;# update j forces with water O 
	movsd [rsp + nb113_fjx], xmm0
	movsd [rsp + nb113_fjy], xmm1
	movsd [rsp + nb113_fjz], xmm2

	;# H1 interactions
	movapd  xmm6, [rsp + nb113_rinvH1] 
	movapd  xmm4, xmm6	
	mulsd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	mulsd  xmm6, [rsp + nb113_qqH]	;# xmm6=vcoul 
	mulsd  xmm4, xmm6		;# total fsH1 in xmm4 
	
	addsd  xmm6, [rsp + nb113_vctot]

	movapd xmm0, [rsp + nb113_dxH1]
	movapd xmm1, [rsp + nb113_dyH1]
	movapd xmm2, [rsp + nb113_dzH1]
	movsd [rsp + nb113_vctot], xmm6
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [rsp + nb113_fixH1]
	movapd xmm4, [rsp + nb113_fiyH1]
	movapd xmm7, [rsp + nb113_fizH1]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [rsp + nb113_fixH1], xmm3
	movsd [rsp + nb113_fiyH1], xmm4
	movsd [rsp + nb113_fizH1], xmm7
	;# update j forces with water H1 
	addsd  xmm0, [rsp + nb113_fjx]
	addsd  xmm1, [rsp + nb113_fjy]
	addsd  xmm2, [rsp + nb113_fjz]
	movsd [rsp + nb113_fjx], xmm0
	movsd [rsp + nb113_fjy], xmm1
	movsd [rsp + nb113_fjz], xmm2

	;# H2 interactions 
	movapd  xmm5, [rsp + nb113_rinvH2] 
	movapd  xmm4, xmm5	
	mulsd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	mulsd  xmm5, [rsp + nb113_qqH]	;# xmm5=vcoul 
	mulsd  xmm4, xmm5		;# total fsH1 in xmm4 
	
	addsd  xmm5, [rsp + nb113_vctot]

	movapd xmm0, [rsp + nb113_dxH2]
	movapd xmm1, [rsp + nb113_dyH2]
	movapd xmm2, [rsp + nb113_dzH2]
	movsd [rsp + nb113_vctot], xmm5
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [rsp + nb113_fixH2]
	movapd xmm4, [rsp + nb113_fiyH2]
	movapd xmm7, [rsp + nb113_fizH2]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [rsp + nb113_fixH2], xmm3
	movsd [rsp + nb113_fiyH2], xmm4
	movsd [rsp + nb113_fizH2], xmm7
	;# update j forces with water H2 
	addsd  xmm0, [rsp + nb113_fjx]
	addsd  xmm1, [rsp + nb113_fjy]
	addsd  xmm2, [rsp + nb113_fjz]
	movsd [rsp + nb113_fjx], xmm0
	movsd [rsp + nb113_fjy], xmm1
	movsd [rsp + nb113_fjz], xmm2

	;# M interactions 
	movapd  xmm5, [rsp + nb113_rinvM] 
	movapd  xmm4, xmm5	
	mulsd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	mulsd  xmm5, [rsp + nb113_qqM]	;# xmm5=vcoul 
	mulsd  xmm4, xmm5		;# total fsH1 in xmm4 
	
	addsd  xmm5, [rsp + nb113_vctot]

	movapd xmm0, [rsp + nb113_dxM]
	movapd xmm1, [rsp + nb113_dyM]
	movapd xmm2, [rsp + nb113_dzM]
	movsd [rsp + nb113_vctot], xmm5
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update M forces 
	movapd xmm3, [rsp + nb113_fixM]
	movapd xmm4, [rsp + nb113_fiyM]
	movapd xmm7, [rsp + nb113_fizM]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [rsp + nb113_fixM], xmm3
	movsd [rsp + nb113_fiyM], xmm4
	movsd [rsp + nb113_fizM], xmm7

	mov rdi, [rbp + nb113_faction]
	;# update j forces 
	addsd  xmm0, [rsp + nb113_fjx]
	addsd  xmm1, [rsp + nb113_fjy]
	addsd  xmm2, [rsp + nb113_fjz]
	movlpd xmm3, [rdi + rax*8]
	movlpd xmm4, [rdi + rax*8 + 8]
	movlpd xmm5, [rdi + rax*8 + 16]
	addsd xmm3, xmm0
	addsd xmm4, xmm1
	addsd xmm5, xmm2
	movlpd [rdi + rax*8], xmm3
	movlpd [rdi + rax*8 + 8], xmm4
	movlpd [rdi + rax*8 + 16], xmm5

.nb113_updateouterdata:
	mov   ecx, [rsp + nb113_ii3]
	mov   rdi, [rbp + nb113_faction]
	mov   rsi, [rbp + nb113_fshift]
	mov   edx, [rsp + nb113_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb113_fixO]
	movapd xmm1, [rsp + nb113_fiyO]
	movapd xmm2, [rsp + nb113_fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

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
	unpcklpd xmm6,xmm1 

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb113_fixH1]
	movapd xmm1, [rsp + nb113_fiyH1]
	movapd xmm2, [rsp + nb113_fizH1]

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
	movapd xmm0, [rsp + nb113_fixH2]
	movapd xmm1, [rsp + nb113_fiyH2]
	movapd xmm2, [rsp + nb113_fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

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
	movapd xmm0, [rsp + nb113_fixM]
	movapd xmm1, [rsp + nb113_fiyM]
	movapd xmm2, [rsp + nb113_fizM]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

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
	mov esi, [rsp + nb113_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb113_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb113_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb113_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb113_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb113_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
       ;# finish if last 
        mov ecx, [rsp + nb113_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb113_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb113_n], esi
        jmp .nb113_outer
.nb113_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb113_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb113_end
        ;# non-zero, do one more workunit
        jmp   .nb113_threadloop
.nb113_end:
	mov eax, [rsp + nb113_nouter]
	mov ebx, [rsp + nb113_ninner]
	mov rcx, [rbp + nb113_outeriter]
	mov rdx, [rbp + nb113_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 976
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





.globl nb_kernel113nf_x86_64_sse2
.globl _nb_kernel113nf_x86_64_sse2
nb_kernel113nf_x86_64_sse2:	
_nb_kernel113nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb113nf_fshift,         16
.equiv          nb113nf_gid,            24
.equiv          nb113nf_pos,            32
.equiv          nb113nf_faction,        40
.equiv          nb113nf_charge,         48
.equiv          nb113nf_p_facel,        56
.equiv          nb113nf_argkrf,         64
.equiv          nb113nf_argcrf,         72
.equiv          nb113nf_Vc,             80
.equiv          nb113nf_type,           88
.equiv          nb113nf_p_ntype,        96
.equiv          nb113nf_vdwparam,       104
.equiv          nb113nf_Vvdw,           112
.equiv          nb113nf_p_tabscale,     120
.equiv          nb113nf_VFtab,          128
.equiv          nb113nf_invsqrta,       136
.equiv          nb113nf_dvda,           144
.equiv          nb113nf_p_gbtabscale,   152
.equiv          nb113nf_GBtab,          160
.equiv          nb113nf_p_nthreads,     168
.equiv          nb113nf_count,          176
.equiv          nb113nf_mtx,            184
.equiv          nb113nf_outeriter,      192
.equiv          nb113nf_inneriter,      200
.equiv          nb113nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb113nf_ixO,            0
.equiv          nb113nf_iyO,            16
.equiv          nb113nf_izO,            32
.equiv          nb113nf_ixH1,           48
.equiv          nb113nf_iyH1,           64
.equiv          nb113nf_izH1,           80
.equiv          nb113nf_ixH2,           96
.equiv          nb113nf_iyH2,           112
.equiv          nb113nf_izH2,           128
.equiv          nb113nf_ixM,            144
.equiv          nb113nf_iyM,            160
.equiv          nb113nf_izM,            176
.equiv          nb113nf_iqH,            192
.equiv          nb113nf_iqM,            208
.equiv          nb113nf_qqH,            224
.equiv          nb113nf_qqM,            240
.equiv          nb113nf_c6,             256
.equiv          nb113nf_c12,            272
.equiv          nb113nf_vctot,          288
.equiv          nb113nf_Vvdwtot,        304
.equiv          nb113nf_half,           320
.equiv          nb113nf_three,          336
.equiv          nb113nf_two,            352
.equiv          nb113nf_is3,            368
.equiv          nb113nf_ii3,            372
.equiv          nb113nf_nri,            376
.equiv          nb113nf_iinr,           384
.equiv          nb113nf_jindex,         392
.equiv          nb113nf_jjnr,           400
.equiv          nb113nf_shift,          408
.equiv          nb113nf_shiftvec,       416
.equiv          nb113nf_facel,          424
.equiv          nb113nf_innerjjnr,      432
.equiv          nb113nf_ntia,           440
.equiv          nb113nf_innerk,         444
.equiv          nb113nf_n,              448
.equiv          nb113nf_nn1,            452
.equiv          nb113nf_nouter,         456
.equiv          nb113nf_ninner,         460

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
	sub rsp, 464		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb113nf_nouter], eax
	mov [rsp + nb113nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb113nf_nri], edi
	mov [rsp + nb113nf_iinr], rsi
	mov [rsp + nb113nf_jindex], rdx
	mov [rsp + nb113nf_jjnr], rcx
	mov [rsp + nb113nf_shift], r8
	mov [rsp + nb113nf_shiftvec], r9
	mov rsi, [rbp + nb113nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb113nf_facel], xmm0

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb113nf_half], eax
	mov [rsp + nb113nf_half+4], ebx
	movsd xmm1, [rsp + nb113nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb113nf_half], xmm1
	movapd [rsp + nb113nf_two], xmm2
	movapd [rsp + nb113nf_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb113nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb113nf_charge]
	movsd xmm3, [rdx + rbx*8 + 8]	
	movsd xmm4, [rdx + rbx*8 + 24]	
	mov rsi, [rbp + nb113nf_p_facel]
	movsd xmm0, [rsi]
	movsd xmm5, [rsp + nb113nf_facel]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb113nf_iqH], xmm3
	movapd [rsp + nb113nf_iqM], xmm4
	
	mov   rdx, [rbp + nb113nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb113nf_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb113nf_ntia], ecx		

.nb113nf_threadloop:
        mov   rsi, [rbp + nb113nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb113nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb113nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb113nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb113nf_n], eax
        mov [rsp + nb113nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb113nf_outerstart
        jmp .nb113nf_end

.nb113nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb113nf_nouter]
	mov [rsp + nb113nf_nouter], ebx

.nb113nf_outer:
	mov   rax, [rsp + nb113nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb113nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb113nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb113nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb113nf_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb113nf_ii3], ebx

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
	movapd [rsp + nb113nf_ixO], xmm3
	movapd [rsp + nb113nf_iyO], xmm4
	movapd [rsp + nb113nf_izO], xmm5
	movapd [rsp + nb113nf_ixH1], xmm6
	movapd [rsp + nb113nf_iyH1], xmm7

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
	movapd [rsp + nb113nf_izH1], xmm6
	movapd [rsp + nb113nf_ixH2], xmm0
	movapd [rsp + nb113nf_iyH2], xmm1
	movapd [rsp + nb113nf_izH2], xmm2
	movapd [rsp + nb113nf_ixM], xmm3
	movapd [rsp + nb113nf_iyM], xmm4
	movapd [rsp + nb113nf_izM], xmm5

	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [rsp + nb113nf_vctot], xmm4
	movapd [rsp + nb113nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb113nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb113nf_pos]
	mov   rdi, [rbp + nb113nf_faction]	
	mov   rax, [rsp + nb113nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb113nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb113nf_ninner]
	mov   [rsp + nb113nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb113nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb113nf_unroll_loop
	jmp   .nb113nf_checksingle
.nb113nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb113nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]

	add qword ptr [rsp + nb113nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb113nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [rsp + nb113nf_iqM]
	mulpd  xmm4, [rsp + nb113nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx

	movapd  [rsp + nb113nf_qqM], xmm3
	movapd  [rsp + nb113nf_qqH], xmm4
	
	mov rsi, [rbp + nb113nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov rsi, [rbp + nb113nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	mov edi, [rsp + nb113nf_ntia]
	add eax, edi
	add ebx, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movlpd xmm7, [rsi + rbx*8]	;# c6b
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 
	movhpd xmm7, [rsi + rbx*8 + 8]	;# c6b c12b 

	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [rsp + nb113nf_c6], xmm4
	movapd [rsp + nb113nf_c12], xmm6
	
	mov rsi, [rbp + nb113nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	movhpd xmm0, [rsi + rbx*8]
	movhpd xmm1, [rsi + rbx*8 + 8]
	movhpd xmm2, [rsi + rbx*8 + 16]		

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [rsp + nb113nf_ixO]
	movapd xmm5, [rsp + nb113nf_iyO]
	movapd xmm6, [rsp + nb113nf_izO]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [rsp + nb113nf_ixH1]
	movapd xmm5, [rsp + nb113nf_iyH1]
	movapd xmm6, [rsp + nb113nf_izH1]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm6, xmm5
	addpd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [rsp + nb113nf_ixH2]
	movapd xmm4, [rsp + nb113nf_iyH2]
	movapd xmm5, [rsp + nb113nf_izH2]

	;# calc dr 
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2

	;# square it 
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	addpd xmm5, xmm4
	addpd xmm5, xmm3

	;# move ixM-izM to xmm2-xmm4  
	movapd xmm3, [rsp + nb113nf_iyM]
	movapd xmm4, [rsp + nb113nf_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [rsp + nb113nf_ixM]
	subpd  xmm2, xmm0	

	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# start with rsqH1 - put seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2
	
	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [rsp + nb113nf_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [rsp + nb113nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb113nf_three]
	subpd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [rsp + nb113nf_half] ;# rinv 
	movapd  xmm6, xmm1	;# rinvH1

	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [rsp + nb113nf_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [rsp + nb113nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb113nf_three]
	subpd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [rsp + nb113nf_half] ;# rinv 
	movapd  xmm5, xmm1	;# rinvH2
	
	;# rsqM - seed in xmm2 
	cvtpd2ps xmm2, xmm4	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [rsp + nb113nf_three]
	mulpd   xmm2, xmm4	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [rsp + nb113nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb113nf_three]
	subpd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [rsp + nb113nf_half] ;# rinv 
	movapd  xmm4, xmm1	;# rinvM

	;# calculate coulomb potentials from rinv.
	addpd   xmm6, xmm5	;# rinvH1+rinvH2
	mulpd	xmm4, [rsp + nb113nf_qqM]
	mulpd	xmm6, [rsp + nb113nf_qqH]
	addpd   xmm4, xmm6
	addpd   xmm4, [rsp + nb113nf_vctot]
	movapd  [rsp + nb113nf_vctot], xmm4

	;# do O interactions - rsqO is in xmm7
	cvtpd2ps xmm2, xmm7
	movapd   xmm6, xmm7
	rcpps    xmm2, xmm2
	cvtps2pd xmm2, xmm2
	movapd   xmm1, [rsp + nb113nf_two]
	movapd   xmm0, xmm1
	mulpd   xmm7, xmm2
	subpd   xmm1, xmm7
	mulpd   xmm2, xmm1 ;# iter1 
	mulpd   xmm6, xmm2
	subpd   xmm0, xmm6
	mulpd   xmm0, xmm2 ;# xmm0=rinvsq
	movapd  xmm1, xmm0	
	mulpd   xmm1, xmm1 ;# rinv4
	mulpd   xmm1, xmm0 ;#rinvsix
	movapd  xmm2, xmm1
	mulpd	xmm2, xmm2 ;# rinvtwelve
	mulpd  xmm1, [rsp + nb113nf_c6]
	mulpd  xmm2, [rsp + nb113nf_c12]
	movapd xmm3, xmm2
	subpd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addpd  xmm3, [rsp + nb113nf_Vvdwtot]
	movapd [rsp + nb113nf_Vvdwtot], xmm3
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb113nf_innerk],  2
	jl   .nb113nf_checksingle
	jmp  .nb113nf_unroll_loop
.nb113nf_checksingle:	
	add dword ptr [rsp + nb113nf_innerk],  2
	jnz  .nb113nf_dosingle
	jmp  .nb113nf_updateouterdata
.nb113nf_dosingle:
	mov   rdx, [rsp + nb113nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb113nf_innerjjnr],  4	

	mov rsi, [rbp + nb113nf_charge]    ;# base of charge[] 

	xorpd xmm3, xmm3
	movlpd xmm3, [rsi + rax*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [rsp + nb113nf_iqM]
	mulpd  xmm4, [rsp + nb113nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 

	movapd  [rsp + nb113nf_qqM], xmm3
	movapd  [rsp + nb113nf_qqH], xmm4
	
	mov rsi, [rbp + nb113nf_type]
	mov eax, [rsi + rax*4]
	mov rsi, [rbp + nb113nf_vdwparam]
	shl eax, 1	
	mov edi, [rsp + nb113nf_ntia]
	add eax, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [rsp + nb113nf_c6], xmm4
	movapd [rsp + nb113nf_c12], xmm6
	
	mov rsi, [rbp + nb113nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [rsp + nb113nf_ixO]
	movapd xmm5, [rsp + nb113nf_iyO]
	movapd xmm6, [rsp + nb113nf_izO]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [rsp + nb113nf_ixH1]
	movapd xmm5, [rsp + nb113nf_iyH1]
	movapd xmm6, [rsp + nb113nf_izH1]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm6, xmm5
	addsd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [rsp + nb113nf_ixH2]
	movapd xmm4, [rsp + nb113nf_iyH2]
	movapd xmm5, [rsp + nb113nf_izH2]

	;# calc dr 
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2

	;# square it 
	mulsd xmm3,xmm3
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	addsd xmm5, xmm4
	addsd xmm5, xmm3
	;# move ixM-izM to xmm2-xmm4  
	movapd xmm3, [rsp + nb113nf_iyM]
	movapd xmm4, [rsp + nb113nf_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [rsp + nb113nf_ixM]
	subpd  xmm2, xmm0	

	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# start with rsqH1 - put seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [rsp + nb113nf_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb113nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb113nf_three]
	subsd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb113nf_half] ;# rinv 
	movapd xmm6, xmm1
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [rsp + nb113nf_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb113nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb113nf_three]
	subsd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb113nf_half] ;# rinv 
	movapd xmm5, xmm1
	
	;# rsqM - seed in xmm2 
	cvtsd2ss xmm2, xmm4
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [rsp + nb113nf_three]
	mulsd   xmm2, xmm4	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb113nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb113nf_three]
	subsd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb113nf_half] ;# rinv 
	movapd xmm4, xmm1

	;# Calculate coulomb potential
	addsd  xmm6, xmm5 	;# rinvH1+rinvH2
	mulsd  xmm4, [rsp + nb113nf_qqM]
	mulsd  xmm6, [rsp + nb113nf_qqH]
	addsd  xmm4, xmm6
	addsd xmm4, [rsp + nb113nf_vctot]
	movsd [rsp + nb113nf_vctot], xmm4

	;# do O interactions directly. xmm7=rsq
	cvtsd2ss xmm2, xmm7
	movapd   xmm6, xmm7
	rcpps    xmm2, xmm2
	cvtss2sd xmm2, xmm2
	movapd   xmm1, [rsp + nb113nf_two]
	movapd   xmm0, xmm1
	mulsd   xmm7, xmm2
	subsd   xmm1, xmm7
	mulsd   xmm2, xmm1 ;# iter1 
	mulsd   xmm6, xmm2
	subsd   xmm0, xmm6
	mulsd   xmm0, xmm2 ;# xmm0=rinvsq
	movapd  xmm1, xmm0	
	mulsd   xmm1, xmm1 ;# rinv4
	mulsd   xmm1, xmm0 ;#rinvsix
	movapd  xmm2, xmm1
	mulsd	xmm2, xmm2 ;# rinvtwelve
	mulsd  xmm1, [rsp + nb113nf_c6]
	mulsd  xmm2, [rsp + nb113nf_c12]
	movapd xmm3, xmm2
	subsd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addsd  xmm3, [rsp + nb113nf_Vvdwtot]
	movsd [rsp + nb113nf_Vvdwtot], xmm3

.nb113nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb113nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb113nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb113nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb113nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb113nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb113nf_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
       ;# finish if last 
        mov ecx, [rsp + nb113nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb113nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb113nf_n], esi
        jmp .nb113nf_outer
.nb113nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb113nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb113nf_end
        ;# non-zero, do one more workunit
        jmp   .nb113nf_threadloop
.nb113nf_end:
	mov eax, [rsp + nb113nf_nouter]
	mov ebx, [rsp + nb113nf_ninner]
	mov rcx, [rbp + nb113nf_outeriter]
	mov rdx, [rbp + nb113nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 464
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
