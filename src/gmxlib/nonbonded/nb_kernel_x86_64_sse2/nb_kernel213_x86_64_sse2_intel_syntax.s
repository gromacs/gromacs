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


.globl nb_kernel213_x86_64_sse2
.globl _nb_kernel213_x86_64_sse2
nb_kernel213_x86_64_sse2:	
_nb_kernel213_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb213_fshift,           16
.equiv          nb213_gid,              24
.equiv          nb213_pos,              32
.equiv          nb213_faction,          40
.equiv          nb213_charge,           48
.equiv          nb213_p_facel,          56
.equiv          nb213_argkrf,           64
.equiv          nb213_argcrf,           72
.equiv          nb213_Vc,               80
.equiv          nb213_type,             88
.equiv          nb213_p_ntype,          96
.equiv          nb213_vdwparam,         104
.equiv          nb213_Vvdw,             112
.equiv          nb213_p_tabscale,       120
.equiv          nb213_VFtab,            128
.equiv          nb213_invsqrta,         136
.equiv          nb213_dvda,             144
.equiv          nb213_p_gbtabscale,     152
.equiv          nb213_GBtab,            160
.equiv          nb213_p_nthreads,       168
.equiv          nb213_count,            176
.equiv          nb213_mtx,              184
.equiv          nb213_outeriter,        192
.equiv          nb213_inneriter,        200
.equiv          nb213_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb213_ixO,              0
.equiv          nb213_iyO,              16
.equiv          nb213_izO,              32
.equiv          nb213_ixH1,             48
.equiv          nb213_iyH1,             64
.equiv          nb213_izH1,             80
.equiv          nb213_ixH2,             96
.equiv          nb213_iyH2,             112
.equiv          nb213_izH2,             128
.equiv          nb213_ixM,              144
.equiv          nb213_iyM,              160
.equiv          nb213_izM,              176
.equiv          nb213_iqH,              192
.equiv          nb213_iqM,              208
.equiv          nb213_dxO,              224
.equiv          nb213_dyO,              240
.equiv          nb213_dzO,              256
.equiv          nb213_dxH1,             272
.equiv          nb213_dyH1,             288
.equiv          nb213_dzH1,             304
.equiv          nb213_dxH2,             320
.equiv          nb213_dyH2,             336
.equiv          nb213_dzH2,             352
.equiv          nb213_dxM,              368
.equiv          nb213_dyM,              384
.equiv          nb213_dzM,              400
.equiv          nb213_qqH,              416
.equiv          nb213_qqM,              432
.equiv          nb213_c6,               448
.equiv          nb213_c12,              464
.equiv          nb213_six,              480
.equiv          nb213_twelve,           496
.equiv          nb213_vctot,            512
.equiv          nb213_Vvdwtot,          528
.equiv          nb213_fixO,             544
.equiv          nb213_fiyO,             560
.equiv          nb213_fizO,             576
.equiv          nb213_fixH1,            592
.equiv          nb213_fiyH1,            608
.equiv          nb213_fizH1,            624
.equiv          nb213_fixH2,            640
.equiv          nb213_fiyH2,            656
.equiv          nb213_fizH2,            672
.equiv          nb213_fixM,             688
.equiv          nb213_fiyM,             704
.equiv          nb213_fizM,             720
.equiv          nb213_fjx,              736
.equiv          nb213_fjy,              752
.equiv          nb213_fjz,              768
.equiv          nb213_half,             784
.equiv          nb213_three,            800
.equiv          nb213_two,              816
.equiv          nb213_rinvH1,           832
.equiv          nb213_rinvH2,           848
.equiv          nb213_rinvM,            864
.equiv          nb213_krsqH1,           880
.equiv          nb213_krsqH2,           896
.equiv          nb213_krsqM,            912
.equiv          nb213_krf,              928
.equiv          nb213_crf,              944
.equiv          nb213_is3,              960
.equiv          nb213_ii3,              964
.equiv          nb213_nri,              968
.equiv          nb213_iinr,             976
.equiv          nb213_jindex,           984
.equiv          nb213_jjnr,             992
.equiv          nb213_shift,            1000
.equiv          nb213_shiftvec,         1008
.equiv          nb213_facel,            1016
.equiv          nb213_innerjjnr,        1024
.equiv          nb213_ntia,             1032
.equiv          nb213_innerk,           1036
.equiv          nb213_n,                1040
.equiv          nb213_nn1,              1044
.equiv          nb213_nouter,           1048
.equiv          nb213_ninner,           1052

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
	sub rsp, 1056		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb213_nouter], eax
	mov [rsp + nb213_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb213_nri], edi
	mov [rsp + nb213_iinr], rsi
	mov [rsp + nb213_jindex], rdx
	mov [rsp + nb213_jjnr], rcx
	mov [rsp + nb213_shift], r8
	mov [rsp + nb213_shiftvec], r9
	mov rsi, [rbp + nb213_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb213_facel], xmm0

	mov rsi, [rbp + nb213_argkrf]
	mov rdi, [rbp + nb213_argcrf]
	movsd xmm1, [rsi]
	movsd xmm2, [rdi]
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	movapd [rsp + nb213_krf], xmm1
	movapd [rsp + nb213_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb213_half], eax
	mov [rsp + nb213_half+4], ebx
	movsd xmm1, [rsp + nb213_half]
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
	movapd [rsp + nb213_half], xmm1
	movapd [rsp + nb213_two], xmm2
	movapd [rsp + nb213_three], xmm3
	movapd [rsp + nb213_six], xmm4
	movapd [rsp + nb213_twelve], xmm5

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb213_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb213_charge]
	movsd xmm3, [rdx + rbx*8 + 8]	
	movsd xmm4, [rdx + rbx*8 + 24]	

	movsd xmm5, [rsp + nb213_facel]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb213_iqH], xmm3
	movapd [rsp + nb213_iqM], xmm4
	
	mov   rdx, [rbp + nb213_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb213_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb213_ntia], ecx		
.nb213_threadloop:
        mov   rsi, [rbp + nb213_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb213_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb213_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb213_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb213_n], eax
        mov [rsp + nb213_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb213_outerstart
        jmp .nb213_end

.nb213_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb213_nouter]
	mov [rsp + nb213_nouter], ebx

.nb213_outer:
	mov   rax, [rsp + nb213_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb213_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb213_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb213_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb213_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb213_ii3], ebx

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
	movapd [rsp + nb213_ixO], xmm3
	movapd [rsp + nb213_iyO], xmm4
	movapd [rsp + nb213_izO], xmm5
	movapd [rsp + nb213_ixH1], xmm6
	movapd [rsp + nb213_iyH1], xmm7

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
	movapd [rsp + nb213_izH1], xmm6
	movapd [rsp + nb213_ixH2], xmm0
	movapd [rsp + nb213_iyH2], xmm1
	movapd [rsp + nb213_izH2], xmm2
	movapd [rsp + nb213_ixM], xmm3
	movapd [rsp + nb213_iyM], xmm4
	movapd [rsp + nb213_izM], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb213_vctot], xmm4
	movapd [rsp + nb213_Vvdwtot], xmm4
	movapd [rsp + nb213_fixO], xmm4
	movapd [rsp + nb213_fiyO], xmm4
	movapd [rsp + nb213_fizO], xmm4
	movapd [rsp + nb213_fixH1], xmm4
	movapd [rsp + nb213_fiyH1], xmm4
	movapd [rsp + nb213_fizH1], xmm4
	movapd [rsp + nb213_fixH2], xmm4
	movapd [rsp + nb213_fiyH2], xmm4
	movapd [rsp + nb213_fizH2], xmm4
	movapd [rsp + nb213_fixM], xmm4
	movapd [rsp + nb213_fiyM], xmm4
	movapd [rsp + nb213_fizM], xmm4
	
	mov   rax, [rsp + nb213_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb213_pos]
	mov   rdi, [rbp + nb213_faction]	
	mov   rax, [rsp + nb213_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb213_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb213_ninner]
	mov   [rsp + nb213_ninner], ecx
	add   edx, 0
	mov   [rsp + nb213_innerk], edx    ;# number of innerloop atoms 
	jge   .nb213_unroll_loop
	jmp   .nb213_checksingle
.nb213_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb213_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]

	add qword ptr [rsp + nb213_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb213_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [rsp + nb213_iqM]
	mulpd  xmm4, [rsp + nb213_iqH]

	movapd  [rsp + nb213_qqM], xmm3
	movapd  [rsp + nb213_qqH], xmm4
	
	mov rsi, [rbp + nb213_type]
	mov r8d, [rsi + rax*4]
	mov r9d, [rsi + rbx*4]
	mov rsi, [rbp + nb213_vdwparam]
	shl r8d, 1	
	shl r9d, 1	
	mov edi, [rsp + nb213_ntia]
	add r8d, edi
	add r9d, edi

	movlpd xmm6, [rsi + r8*8]	;# c6a
	movlpd xmm7, [rsi + r9*8]	;# c6b
	movhpd xmm6, [rsi + r8*8 + 8]	;# c6a c12a 
	movhpd xmm7, [rsi + r9*8 + 8]	;# c6b c12b 
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7

	movapd [rsp + nb213_c6], xmm4
	movapd [rsp + nb213_c12], xmm6
	
	mov rsi, [rbp + nb213_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# load j coordinates 
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
    
    subpd xmm3, [rsp + nb213_ixO]
    subpd xmm4, [rsp + nb213_iyO]
    subpd xmm5, [rsp + nb213_izO]
    
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
    movapd xmm4, [rsp + nb213_two]
    movapd xmm5, xmm3
    mulpd xmm3, xmm6        ;# lu*rsq 
    subpd xmm4, xmm3        ;# 2-lu*rsq 
    mulpd xmm6, xmm4        ;# (new lu) 
        
    movapd xmm4, [rsp + nb213_two]
    mulpd xmm5, xmm6        ;# lu*rsq 
    subpd xmm4, xmm5        ;# 2-lu*rsq 
    mulpd xmm4, xmm6        ;# xmm4=rinvsq 

    movapd xmm3, xmm4       ;# rinvsq
    mulpd  xmm4, xmm4       ;# rinv4
    mulpd  xmm4, xmm3       ;# rinv6
    movapd xmm5, xmm4     
    mulpd  xmm5, xmm5       ;# rinv12
    mulpd  xmm4, [rsp + nb213_c6]
    mulpd  xmm5, [rsp + nb213_c12]
    movapd xmm6, xmm5
    subpd  xmm6, xmm4  ;# Vvdw=vvdw12-vvdw6
    mulpd  xmm4, [rsp + nb213_six]
    mulpd  xmm5, [rsp + nb213_twelve]
    subpd  xmm5, xmm4
    mulpd  xmm3, xmm5   ;# fscal
    
    addpd  xmm6, [rsp + nb213_Vvdwtot]
    movapd [rsp + nb213_Vvdwtot], xmm6
    
    mulpd  xmm13, xmm3 ;# fx
    mulpd  xmm14, xmm3 ;# fy
    mulpd  xmm15, xmm3 ;# fz

    ;# save j force temporarily
    movapd [rsp + nb213_fjx], xmm13
    movapd [rsp + nb213_fjy], xmm14
    movapd [rsp + nb213_fjz], xmm15
    
    ;# increment i O force
    addpd xmm13, [rsp + nb213_fixO]
    addpd xmm14, [rsp + nb213_fiyO]
    addpd xmm15, [rsp + nb213_fizO]
    movapd [rsp + nb213_fixO], xmm13
    movapd [rsp + nb213_fiyO], xmm14
    movapd [rsp + nb213_fizO], xmm15    
    ;# finished O LJ interaction.


    ;# do H1, H2, and M interactions in parallel.
    ;# xmm0-xmm2 still contain j coordinates.                
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subpd xmm0, [rsp + nb213_ixH1]
    subpd xmm1, [rsp + nb213_iyH1]
    subpd xmm2, [rsp + nb213_izH1]
    subpd xmm3, [rsp + nb213_ixH2]
    subpd xmm4, [rsp + nb213_iyH2]
    subpd xmm5, [rsp + nb213_izH2]
    subpd xmm6, [rsp + nb213_ixM]
    subpd xmm7, [rsp + nb213_iyM]
    subpd xmm8, [rsp + nb213_izM]
    
	movapd [rsp + nb213_dxH1], xmm0
	movapd [rsp + nb213_dyH1], xmm1
	movapd [rsp + nb213_dzH1], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb213_dxH2], xmm3
	movapd [rsp + nb213_dyH2], xmm4
	movapd [rsp + nb213_dzH2], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb213_dxM], xmm6
	movapd [rsp + nb213_dyM], xmm7
	movapd [rsp + nb213_dzM], xmm8
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
		
	movapd  xmm9, [rsp + nb213_three]
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

	movapd  xmm15, [rsp + nb213_half]
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
		
	movapd  xmm1, [rsp + nb213_three]
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

	movapd  xmm15, [rsp + nb213_half]
	mulpd   xmm9, xmm15  ;#  rinvH1
	mulpd   xmm10, xmm15 ;#   rinvH2
    mulpd   xmm11, xmm15 ;#   rinvM
	
	;# interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movapd xmm1, xmm9 ;# copy of rinv
    movapd xmm4, xmm10
    movapd xmm7, xmm11
    movapd xmm2, [rsp + nb213_krf]    
    mulpd  xmm9, xmm9   ;# rinvsq
    mulpd  xmm10, xmm10
    mulpd  xmm11, xmm11
    mulpd  xmm0, xmm2  ;# k*rsq
    mulpd  xmm3, xmm2
    mulpd  xmm6, xmm2
    movapd xmm2, xmm0 ;# copy of k*rsq
    movapd xmm5, xmm3
    movapd xmm8, xmm6
    addpd  xmm2, xmm1  ;# rinv+krsq
    addpd  xmm5, xmm4
    addpd  xmm8, xmm7
    movapd xmm14, [rsp + nb213_crf]
    subpd  xmm2, xmm14   ;# rinv+krsq-crf
    subpd  xmm5, xmm14
    subpd  xmm8, xmm14
    movapd xmm12, [rsp + nb213_qqH]
    movapd xmm13, [rsp + nb213_qqM]    
    mulpd  xmm2, xmm12 ;# voul=qq*(rinv+ krsq-crf)
    mulpd  xmm5, xmm12 ;# voul=qq*(rinv+ krsq-crf)
    mulpd  xmm8, xmm13 ;# voul=qq*(rinv+ krsq-crf)
    addpd  xmm0, xmm0 ;# 2*krsq
    addpd  xmm3, xmm3 
    addpd  xmm6, xmm6 
    subpd  xmm1, xmm0 ;# rinv-2*krsq
    subpd  xmm4, xmm3
    subpd  xmm7, xmm6
    mulpd  xmm1, xmm12   ;# (rinv-2*krsq)*qq
    mulpd  xmm4, xmm12
    mulpd  xmm7, xmm13
    addpd  xmm2, [rsp + nb213_vctot]
    addpd  xmm5, xmm8
    addpd  xmm2, xmm5
    movapd [rsp + nb213_vctot], xmm2
    
    mulpd  xmm9, xmm1   ;# fscal
    mulpd  xmm10, xmm4
    mulpd  xmm11, xmm7

    ;# move j forces to xmm0-xmm2
    mov rdi, [rbp + nb213_faction]
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
    addpd xmm0, [rsp + nb213_fjx]
    addpd xmm1, [rsp + nb213_fjy]
    addpd xmm2, [rsp + nb213_fjz]

	mulpd xmm7, [rsp + nb213_dxH1]
	mulpd xmm8, [rsp + nb213_dyH1]
	mulpd xmm9, [rsp + nb213_dzH1]
	mulpd xmm10, [rsp + nb213_dxH2]
	mulpd xmm11, [rsp + nb213_dyH2]
	mulpd xmm12, [rsp + nb213_dzH2]
	mulpd xmm13, [rsp + nb213_dxM]
	mulpd xmm14, [rsp + nb213_dyM]
	mulpd xmm15, [rsp + nb213_dzM]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb213_fixH1]
    addpd xmm8, [rsp + nb213_fiyH1]
    addpd xmm9, [rsp + nb213_fizH1]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb213_fixH2]
    addpd xmm11, [rsp + nb213_fiyH2]
    addpd xmm12, [rsp + nb213_fizH2]

    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    addpd xmm13, [rsp + nb213_fixM]
    addpd xmm14, [rsp + nb213_fiyM]
    addpd xmm15, [rsp + nb213_fizM]

    movapd [rsp + nb213_fixH1], xmm7
    movapd [rsp + nb213_fiyH1], xmm8
    movapd [rsp + nb213_fizH1], xmm9
    movapd [rsp + nb213_fixH2], xmm10
    movapd [rsp + nb213_fiyH2], xmm11
    movapd [rsp + nb213_fizH2], xmm12
    movapd [rsp + nb213_fixM], xmm13
    movapd [rsp + nb213_fiyM], xmm14
    movapd [rsp + nb213_fizM], xmm15
   
    ;# store back j forces from xmm0-xmm2
	movlpd [rdi + rax*8], xmm0
	movlpd [rdi + rax*8 + 8], xmm1
	movlpd [rdi + rax*8 + 16], xmm2
	movhpd [rdi + rbx*8], xmm0
	movhpd [rdi + rbx*8 + 8], xmm1
	movhpd [rdi + rbx*8 + 16], xmm2

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb213_innerk],  2
	jl   .nb213_checksingle
	jmp  .nb213_unroll_loop
.nb213_checksingle:	
	mov   edx, [rsp + nb213_innerk]
	and   edx, 1
	jnz  .nb213_dosingle
	jmp  .nb213_updateouterdata
.nb213_dosingle:
	mov   rdx, [rsp + nb213_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb213_innerjjnr],  4	

	mov rsi, [rbp + nb213_charge]    ;# base of charge[] 

	xorpd xmm3, xmm3
	movlpd xmm3, [rsi + rax*8]
	movapd xmm4, xmm3
	mulsd  xmm3, [rsp + nb213_iqM]
	mulsd  xmm4, [rsp + nb213_iqH]

	movapd  [rsp + nb213_qqM], xmm3
	movapd  [rsp + nb213_qqH], xmm4
	
	mov rsi, [rbp + nb213_type]
	mov r8d, [rsi + rax*4]
	mov rsi, [rbp + nb213_vdwparam]
	shl r8d, 1	
	mov edi, [rsp + nb213_ntia]
	add r8d, edi

	movlpd xmm6, [rsi + r8*8]	;# c6a
	movhpd xmm6, [rsi + r8*8 + 8]	;# c6a c12a 
	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movapd [rsp + nb213_c6], xmm4
	movapd [rsp + nb213_c12], xmm6
	
	mov rsi, [rbp + nb213_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2  and xmm4-xmm6
	movlpd xmm4, [rsi + rax*8]
	movlpd xmm5, [rsi + rax*8 + 8]
	movlpd xmm6, [rsi + rax*8 + 16]
    movapd xmm0, xmm4
    movapd xmm1, xmm5
    movapd xmm2, xmm6

	;# calc dr 
	subsd xmm4, [rsp + nb213_ixO]
	subsd xmm5, [rsp + nb213_iyO]
	subsd xmm6, [rsp + nb213_izO]

	;# store dr 
	movapd [rsp + nb213_dxO], xmm4
	movapd [rsp + nb213_dyO], xmm5
	movapd [rsp + nb213_dzO], xmm6
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
	subsd xmm4, [rsp + nb213_ixH1]
	subsd xmm5, [rsp + nb213_iyH1]
	subsd xmm6, [rsp + nb213_izH1]

	;# store dr 
	movapd [rsp + nb213_dxH1], xmm4
	movapd [rsp + nb213_dyH1], xmm5
	movapd [rsp + nb213_dzH1], xmm6
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
	subsd xmm3, [rsp + nb213_ixH2]
	subsd xmm4, [rsp + nb213_iyH2]
	subsd xmm5, [rsp + nb213_izH2]

	;# store dr 
	movapd [rsp + nb213_dxH2], xmm3
	movapd [rsp + nb213_dyH2], xmm4
	movapd [rsp + nb213_dzH2], xmm5
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
	subsd xmm4, [rsp + nb213_ixM]
	subsd xmm3, [rsp + nb213_iyM]
	subsd xmm2, [rsp + nb213_izM]

	;# store dr 
	movapd [rsp + nb213_dxM], xmm4
	movapd [rsp + nb213_dyM], xmm3
	movapd [rsp + nb213_dzM], xmm2

	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# calculate krsq
	movsd xmm0, [rsp + nb213_krf]
	movsd xmm1, xmm0
	movsd xmm2, xmm0
	mulsd xmm0, xmm4  
	mulsd xmm1, xmm5
	mulsd xmm2, xmm6
	movsd [rsp + nb213_krsqM], xmm0
	movsd [rsp + nb213_krsqH2], xmm1
	movsd [rsp + nb213_krsqH1], xmm2

	;# start with rsqH1 - put seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [rsp + nb213_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb213_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb213_three]
	subsd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb213_half] ;# rinv 
	movapd [rsp + nb213_rinvH1], xmm1
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [rsp + nb213_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb213_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb213_three]
	subsd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb213_half] ;# rinv 
	movapd [rsp + nb213_rinvH2], xmm1
	
	;# rsqM - seed in xmm2 
	cvtsd2ss xmm2, xmm4
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [rsp + nb213_three]
	mulsd   xmm2, xmm4	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb213_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb213_three]
	subsd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb213_half] ;# rinv 
	movapd [rsp + nb213_rinvM], xmm1

	;# do O interactions directly. xmm7=rsq
	cvtsd2ss xmm2, xmm7
	movapd   xmm6, xmm7
	rcpps    xmm2, xmm2
	cvtss2sd xmm2, xmm2
	movapd   xmm1, [rsp + nb213_two]
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
	mulsd  xmm1, [rsp + nb213_c6]
	mulsd  xmm2, [rsp + nb213_c12]
	movapd xmm3, xmm2
	subsd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addsd  xmm3, [rsp + nb213_Vvdwtot]
	mulsd  xmm1, [rsp + nb213_six]
	mulsd  xmm2, [rsp + nb213_twelve]
	subsd  xmm2, xmm1
	mulsd  xmm2, xmm0
	movapd xmm4, xmm2 ;# total fsO 
	movsd [rsp + nb213_Vvdwtot], xmm3

	movapd xmm0, [rsp + nb213_dxO]
	movapd xmm1, [rsp + nb213_dyO]
	movapd xmm2, [rsp + nb213_dzO]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update O forces 
	movapd xmm3, [rsp + nb213_fixO]
	movapd xmm4, [rsp + nb213_fiyO]
	movapd xmm7, [rsp + nb213_fizO]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [rsp + nb213_fixO], xmm3
	movsd [rsp + nb213_fiyO], xmm4
	movsd [rsp + nb213_fizO], xmm7
	;# update j forces with water O 
	movsd [rsp + nb213_fjx], xmm0
	movsd [rsp + nb213_fjy], xmm1
	movsd [rsp + nb213_fjz], xmm2

	;# H1 interactions
	movsd  xmm6, [rsp + nb213_rinvH1] 
	movsd  xmm4, xmm6
	mulsd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movsd  xmm7, xmm6
	movsd  xmm0, [rsp + nb213_krsqH1]
	addsd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	mulsd   xmm0, [rsp + nb213_two]
	subsd   xmm6, [rsp + nb213_crf]
	subsd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulsd   xmm6, [rsp + nb213_qqH] ;# vcoul 
	mulsd   xmm7, [rsp + nb213_qqH]
	mulsd  xmm4, xmm7		;# total fsH1 in xmm4 
	
	addsd  xmm6, [rsp + nb213_vctot]

	movapd xmm0, [rsp + nb213_dxH1]
	movapd xmm1, [rsp + nb213_dyH1]
	movapd xmm2, [rsp + nb213_dzH1]
	movsd [rsp + nb213_vctot], xmm6
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [rsp + nb213_fixH1]
	movapd xmm4, [rsp + nb213_fiyH1]
	movapd xmm7, [rsp + nb213_fizH1]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [rsp + nb213_fixH1], xmm3
	movsd [rsp + nb213_fiyH1], xmm4
	movsd [rsp + nb213_fizH1], xmm7
	;# update j forces with water H1 
	addsd  xmm0, [rsp + nb213_fjx]
	addsd  xmm1, [rsp + nb213_fjy]
	addsd  xmm2, [rsp + nb213_fjz]
	movsd [rsp + nb213_fjx], xmm0
	movsd [rsp + nb213_fjy], xmm1
	movsd [rsp + nb213_fjz], xmm2

	;# H2 interactions 
	movsd  xmm5, [rsp + nb213_rinvH2] 
	movsd  xmm4, xmm5	
	mulsd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movsd  xmm7, xmm5
	movsd  xmm0, [rsp + nb213_krsqH2]
	addsd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	mulsd   xmm0, [rsp + nb213_two]
	subsd   xmm5, [rsp + nb213_crf]
	subsd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulsd   xmm5, [rsp + nb213_qqH] ;# vcoul 
	mulsd   xmm7, [rsp + nb213_qqH]
	mulsd  xmm4, xmm7		;# total fsH2 in xmm4 
	
	addsd  xmm5, [rsp + nb213_vctot]

	movapd xmm0, [rsp + nb213_dxH2]
	movapd xmm1, [rsp + nb213_dyH2]
	movapd xmm2, [rsp + nb213_dzH2]
	movsd [rsp + nb213_vctot], xmm5
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [rsp + nb213_fixH2]
	movapd xmm4, [rsp + nb213_fiyH2]
	movapd xmm7, [rsp + nb213_fizH2]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [rsp + nb213_fixH2], xmm3
	movsd [rsp + nb213_fiyH2], xmm4
	movsd [rsp + nb213_fizH2], xmm7
	;# update j forces with water H2 
	addsd  xmm0, [rsp + nb213_fjx]
	addsd  xmm1, [rsp + nb213_fjy]
	addsd  xmm2, [rsp + nb213_fjz]
	movsd [rsp + nb213_fjx], xmm0
	movsd [rsp + nb213_fjy], xmm1
	movsd [rsp + nb213_fjz], xmm2

	;# M interactions 
	movsd  xmm5, [rsp + nb213_rinvM] 
	movsd  xmm4, xmm5	
	mulsd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movsd  xmm7, xmm5
	movsd  xmm0, [rsp + nb213_krsqM]
	addsd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	mulsd   xmm0, [rsp + nb213_two]
	subsd   xmm5, [rsp + nb213_crf]
	subsd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulsd   xmm5, [rsp + nb213_qqM] ;# vcoul 
	mulsd   xmm7, [rsp + nb213_qqM]
	mulsd  xmm4, xmm7		;# total fsH2 in xmm4 
	
	addsd  xmm5, [rsp + nb213_vctot]

	movapd xmm0, [rsp + nb213_dxM]
	movapd xmm1, [rsp + nb213_dyM]
	movapd xmm2, [rsp + nb213_dzM]
	movsd [rsp + nb213_vctot], xmm5
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update M forces 
	movapd xmm3, [rsp + nb213_fixM]
	movapd xmm4, [rsp + nb213_fiyM]
	movapd xmm7, [rsp + nb213_fizM]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [rsp + nb213_fixM], xmm3
	movsd [rsp + nb213_fiyM], xmm4
	movsd [rsp + nb213_fizM], xmm7

	mov rdi, [rbp + nb213_faction]
	;# update j forces 
	addsd  xmm0, [rsp + nb213_fjx]
	addsd  xmm1, [rsp + nb213_fjy]
	addsd  xmm2, [rsp + nb213_fjz]
	movlpd xmm3, [rdi + rax*8]
	movlpd xmm4, [rdi + rax*8 + 8]
	movlpd xmm5, [rdi + rax*8 + 16]
	addsd xmm3, xmm0
	addsd xmm4, xmm1
	addsd xmm5, xmm2
	movlpd [rdi + rax*8], xmm3
	movlpd [rdi + rax*8 + 8], xmm4
	movlpd [rdi + rax*8 + 16], xmm5

.nb213_updateouterdata:
	mov   ecx, [rsp + nb213_ii3]
	mov   rdi, [rbp + nb213_faction]
	mov   rsi, [rbp + nb213_fshift]
	mov   edx, [rsp + nb213_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb213_fixO]
	movapd xmm1, [rsp + nb213_fiyO]
	movapd xmm2, [rsp + nb213_fizO]

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
	movapd xmm0, [rsp + nb213_fixH1]
	movapd xmm1, [rsp + nb213_fiyH1]
	movapd xmm2, [rsp + nb213_fizH1]

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
	movapd xmm0, [rsp + nb213_fixH2]
	movapd xmm1, [rsp + nb213_fiyH2]
	movapd xmm2, [rsp + nb213_fizH2]

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
	movapd xmm0, [rsp + nb213_fixM]
	movapd xmm1, [rsp + nb213_fiyM]
	movapd xmm2, [rsp + nb213_fizM]

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
	mov esi, [rsp + nb213_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb213_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb213_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb213_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb213_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb213_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
       ;# finish if last 
        mov ecx, [rsp + nb213_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb213_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb213_n], esi
        jmp .nb213_outer
.nb213_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb213_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb213_end
        ;# non-zero, do one more workunit
        jmp   .nb213_threadloop
.nb213_end:
	mov eax, [rsp + nb213_nouter]
	mov ebx, [rsp + nb213_ninner]
	mov rcx, [rbp + nb213_outeriter]
	mov rdx, [rbp + nb213_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 1056
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





.globl nb_kernel213nf_x86_64_sse2
.globl _nb_kernel213nf_x86_64_sse2
nb_kernel213nf_x86_64_sse2:	
_nb_kernel213nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb213nf_fshift,         16
.equiv          nb213nf_gid,            24
.equiv          nb213nf_pos,            32
.equiv          nb213nf_faction,        40
.equiv          nb213nf_charge,         48
.equiv          nb213nf_p_facel,        56
.equiv          nb213nf_argkrf,         64
.equiv          nb213nf_argcrf,         72
.equiv          nb213nf_Vc,             80
.equiv          nb213nf_type,           88
.equiv          nb213nf_p_ntype,        96
.equiv          nb213nf_vdwparam,       104
.equiv          nb213nf_Vvdw,           112
.equiv          nb213nf_p_tabscale,     120
.equiv          nb213nf_VFtab,          128
.equiv          nb213nf_invsqrta,       136
.equiv          nb213nf_dvda,           144
.equiv          nb213nf_p_gbtabscale,   152
.equiv          nb213nf_GBtab,          160
.equiv          nb213nf_p_nthreads,     168
.equiv          nb213nf_count,          176
.equiv          nb213nf_mtx,            184
.equiv          nb213nf_outeriter,      192
.equiv          nb213nf_inneriter,      200
.equiv          nb213nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb213nf_ixO,            0
.equiv          nb213nf_iyO,            16
.equiv          nb213nf_izO,            32
.equiv          nb213nf_ixH1,           48
.equiv          nb213nf_iyH1,           64
.equiv          nb213nf_izH1,           80
.equiv          nb213nf_ixH2,           96
.equiv          nb213nf_iyH2,           112
.equiv          nb213nf_izH2,           128
.equiv          nb213nf_ixM,            144
.equiv          nb213nf_iyM,            160
.equiv          nb213nf_izM,            176
.equiv          nb213nf_iqH,            192
.equiv          nb213nf_iqM,            208
.equiv          nb213nf_qqH,            224
.equiv          nb213nf_qqM,            240
.equiv          nb213nf_c6,             256
.equiv          nb213nf_c12,            272
.equiv          nb213nf_vctot,          288
.equiv          nb213nf_Vvdwtot,        304
.equiv          nb213nf_half,           320
.equiv          nb213nf_three,          336
.equiv          nb213nf_two,            352
.equiv          nb213nf_rinvH1,         368
.equiv          nb213nf_rinvH2,         384
.equiv          nb213nf_rinvM,          400
.equiv          nb213nf_krsqH1,         416
.equiv          nb213nf_krsqH2,         432
.equiv          nb213nf_krsqM,          448
.equiv          nb213nf_krf,            464
.equiv          nb213nf_crf,            480
.equiv          nb213nf_nri,            496
.equiv          nb213nf_iinr,           504
.equiv          nb213nf_jindex,         512
.equiv          nb213nf_jjnr,           520
.equiv          nb213nf_shift,          528
.equiv          nb213nf_shiftvec,       536
.equiv          nb213nf_facel,          544
.equiv          nb213nf_innerjjnr,      552
.equiv          nb213nf_is3,            560
.equiv          nb213nf_ii3,            564
.equiv          nb213nf_ntia,           568
.equiv          nb213nf_innerk,         572
.equiv          nb213nf_n,              576
.equiv          nb213nf_nn1,            580
.equiv          nb213nf_nouter,         584
.equiv          nb213nf_ninner,         588

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
	sub rsp, 592		;# local variable stack space (n*16+8)
	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb213nf_nouter], eax
	mov [rsp + nb213nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb213nf_nri], edi
	mov [rsp + nb213nf_iinr], rsi
	mov [rsp + nb213nf_jindex], rdx
	mov [rsp + nb213nf_jjnr], rcx
	mov [rsp + nb213nf_shift], r8
	mov [rsp + nb213nf_shiftvec], r9
	mov rsi, [rbp + nb213nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb213nf_facel], xmm0

	mov rsi, [rbp + nb213nf_argkrf]
	mov rdi, [rbp + nb213nf_argcrf]
	movsd xmm1, [rsi]
	movsd xmm2, [rdi]
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	movapd [rsp + nb213nf_krf], xmm1
	movapd [rsp + nb213nf_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb213nf_half], eax
	mov [rsp + nb213nf_half+4], ebx
	movsd xmm1, [rsp + nb213nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb213nf_half], xmm1
	movapd [rsp + nb213nf_two], xmm2
	movapd [rsp + nb213nf_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb213nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb213nf_charge]
	movsd xmm3, [rdx + rbx*8 + 8]	
	movsd xmm4, [rdx + rbx*8 + 24]	

	movsd xmm5, [rsp + nb213nf_facel]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb213nf_iqH], xmm3
	movapd [rsp + nb213nf_iqM], xmm4
	
	mov   rdx, [rbp + nb213nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb213nf_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb213nf_ntia], ecx		
.nb213nf_threadloop:
        mov   rsi, [rbp + nb213nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb213nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb213nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb213nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb213nf_n], eax
        mov [rsp + nb213nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb213nf_outerstart
        jmp .nb213nf_end

.nb213nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb213nf_nouter]
	mov [rsp + nb213nf_nouter], ebx

.nb213nf_outer:
	mov   rax, [rsp + nb213nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb213nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb213nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb213nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb213nf_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb213nf_ii3], ebx

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
	movapd [rsp + nb213nf_ixO], xmm3
	movapd [rsp + nb213nf_iyO], xmm4
	movapd [rsp + nb213nf_izO], xmm5
	movapd [rsp + nb213nf_ixH1], xmm6
	movapd [rsp + nb213nf_iyH1], xmm7

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
	movapd [rsp + nb213nf_izH1], xmm6
	movapd [rsp + nb213nf_ixH2], xmm0
	movapd [rsp + nb213nf_iyH2], xmm1
	movapd [rsp + nb213nf_izH2], xmm2
	movapd [rsp + nb213nf_ixM], xmm3
	movapd [rsp + nb213nf_iyM], xmm4
	movapd [rsp + nb213nf_izM], xmm5

	;# clear vctot
	xorpd xmm4, xmm4
	movapd [rsp + nb213nf_vctot], xmm4
	movapd [rsp + nb213nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb213nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb213nf_pos]
	mov   rdi, [rbp + nb213nf_faction]	
	mov   rax, [rsp + nb213nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb213nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb213nf_ninner]
	mov   [rsp + nb213nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb213nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb213nf_unroll_loop
	jmp   .nb213nf_checksingle
.nb213nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb213nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]

	add qword ptr [rsp + nb213nf_innerjjnr],  8	
	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb213nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [rsp + nb213nf_iqM]
	mulpd  xmm4, [rsp + nb213nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx

	movapd  [rsp + nb213nf_qqM], xmm3
	movapd  [rsp + nb213nf_qqH], xmm4
	
	mov rsi, [rbp + nb213nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov rsi, [rbp + nb213nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	mov edi, [rsp + nb213nf_ntia]
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
	movapd [rsp + nb213nf_c6], xmm4
	movapd [rsp + nb213nf_c12], xmm6
	
	mov rsi, [rbp + nb213nf_pos]       ;# base of pos[] 

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
	movapd xmm4, [rsp + nb213nf_ixO]
	movapd xmm5, [rsp + nb213nf_iyO]
	movapd xmm6, [rsp + nb213nf_izO]

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
	movapd xmm4, [rsp + nb213nf_ixH1]
	movapd xmm5, [rsp + nb213nf_iyH1]
	movapd xmm6, [rsp + nb213nf_izH1]

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
	movapd xmm3, [rsp + nb213nf_ixH2]
	movapd xmm4, [rsp + nb213nf_iyH2]
	movapd xmm5, [rsp + nb213nf_izH2]

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
	movapd xmm3, [rsp + nb213nf_iyM]
	movapd xmm4, [rsp + nb213nf_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [rsp + nb213nf_ixM]
	subpd  xmm2, xmm0	

	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# calculate krsq
	movapd xmm0, [rsp + nb213nf_krf]
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	mulpd xmm0, xmm4  
	mulpd xmm1, xmm5
	mulpd xmm2, xmm6
	movapd [rsp + nb213nf_krsqM], xmm0
	movapd [rsp + nb213nf_krsqH2], xmm1
	movapd [rsp + nb213nf_krsqH1], xmm2

	;# start with rsqH1 - put seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2
	
	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [rsp + nb213nf_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [rsp + nb213nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb213nf_three]
	subpd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [rsp + nb213nf_half] ;# rinv 
	movapd  [rsp + nb213nf_rinvH1], xmm1	

	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [rsp + nb213nf_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [rsp + nb213nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb213nf_three]
	subpd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [rsp + nb213nf_half] ;# rinv 
	movapd  [rsp + nb213nf_rinvH2], xmm1	
	
	;# rsqM - seed in xmm2 
	cvtpd2ps xmm2, xmm4	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [rsp + nb213nf_three]
	mulpd   xmm2, xmm4	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [rsp + nb213nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb213nf_three]
	subpd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [rsp + nb213nf_half] ;# rinv 
	movapd  [rsp + nb213nf_rinvM], xmm1	

	;# do O interactions directly - rsqO is in xmm7
	cvtpd2ps xmm2, xmm7
	movapd   xmm6, xmm7
	rcpps    xmm2, xmm2
	cvtps2pd xmm2, xmm2
	movapd   xmm1, [rsp + nb213nf_two]
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
	mulpd  xmm1, [rsp + nb213nf_c6]
	mulpd  xmm2, [rsp + nb213nf_c12]
	movapd xmm3, xmm2
	subpd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addpd  xmm3, [rsp + nb213nf_Vvdwtot]
	movapd [rsp + nb213nf_Vvdwtot], xmm3

	;# H1 interactions 
	movapd  xmm6, [rsp + nb213nf_rinvH1]
	addpd   xmm6, [rsp + nb213nf_krsqH1]
 	subpd   xmm6, [rsp + nb213nf_crf]
	mulpd   xmm6, [rsp + nb213nf_qqH] ;# vcoul 
	addpd   xmm6, [rsp + nb213nf_vctot]
	movapd [rsp + nb213nf_vctot], xmm6
	
	;# H2 interactions 
	movapd  xmm6, [rsp + nb213nf_rinvH2]
	addpd   xmm6, [rsp + nb213nf_krsqH2]
 	subpd   xmm6, [rsp + nb213nf_crf]
	mulpd   xmm6, [rsp + nb213nf_qqH] ;# vcoul 
	addpd   xmm6, [rsp + nb213nf_vctot]
	movapd [rsp + nb213nf_vctot], xmm6

	;# M interactions 
	movapd  xmm6, [rsp + nb213nf_rinvM]
	addpd   xmm6, [rsp + nb213nf_krsqM]
 	subpd   xmm6, [rsp + nb213nf_crf]
	mulpd   xmm6, [rsp + nb213nf_qqM] ;# vcoul 
	addpd   xmm6, [rsp + nb213nf_vctot]
	movapd [rsp + nb213nf_vctot], xmm6
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb213nf_innerk],  2
	jl   .nb213nf_checksingle
	jmp  .nb213nf_unroll_loop
.nb213nf_checksingle:	
	mov   edx, [rsp + nb213nf_innerk]
	and   edx, 1
	jnz  .nb213nf_dosingle
	jmp  .nb213nf_updateouterdata
.nb213nf_dosingle:
	mov   rdx, [rsp + nb213nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb213nf_innerjjnr],  4	

	mov rsi, [rbp + nb213nf_charge]    ;# base of charge[] 

	xorpd xmm3, xmm3
	movlpd xmm3, [rsi + rax*8]
	movapd xmm4, xmm3
	mulsd  xmm3, [rsp + nb213nf_iqM]
	mulsd  xmm4, [rsp + nb213nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 

	movapd  [rsp + nb213nf_qqM], xmm3
	movapd  [rsp + nb213nf_qqH], xmm4
	
	mov rsi, [rbp + nb213nf_type]
	mov eax, [rsi + rax*4]
	mov rsi, [rbp + nb213nf_vdwparam]
	shl eax, 1	
	mov edi, [rsp + nb213nf_ntia]
	add eax, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 
	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [rsp + nb213nf_c6], xmm4
	movapd [rsp + nb213nf_c12], xmm6
	
	mov rsi, [rbp + nb213nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [rsp + nb213nf_ixO]
	movapd xmm5, [rsp + nb213nf_iyO]
	movapd xmm6, [rsp + nb213nf_izO]

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
	movapd xmm4, [rsp + nb213nf_ixH1]
	movapd xmm5, [rsp + nb213nf_iyH1]
	movapd xmm6, [rsp + nb213nf_izH1]

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
	movapd xmm3, [rsp + nb213nf_ixH2]
	movapd xmm4, [rsp + nb213nf_iyH2]
	movapd xmm5, [rsp + nb213nf_izH2]

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
	movapd xmm3, [rsp + nb213nf_iyM]
	movapd xmm4, [rsp + nb213nf_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [rsp + nb213nf_ixM]
	subpd  xmm2, xmm0	

	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# calculate krsq
	movsd xmm0, [rsp + nb213nf_krf]
	movsd xmm1, xmm0
	movsd xmm2, xmm0
	mulsd xmm0, xmm4  
	mulsd xmm1, xmm5
	mulsd xmm2, xmm6
	movsd [rsp + nb213nf_krsqM], xmm0
	movsd [rsp + nb213nf_krsqH2], xmm1
	movsd [rsp + nb213nf_krsqH1], xmm2

	;# start with rsqH1 - put seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [rsp + nb213nf_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb213nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb213nf_three]
	subsd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb213nf_half] ;# rinv 
	movapd [rsp + nb213nf_rinvH1], xmm1
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [rsp + nb213nf_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb213nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb213nf_three]
	subsd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb213nf_half] ;# rinv 
	movapd [rsp + nb213nf_rinvH2], xmm1
	
	;# rsqM - seed in xmm2 
	cvtsd2ss xmm2, xmm4
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [rsp + nb213nf_three]
	mulsd   xmm2, xmm4	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb213nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb213nf_three]
	subsd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb213nf_half] ;# rinv 
	movapd [rsp + nb213nf_rinvM], xmm1

	;# do O interactions directly. xmm7=rsq
	cvtsd2ss xmm2, xmm7
	movapd   xmm6, xmm7
	rcpps    xmm2, xmm2
	cvtss2sd xmm2, xmm2
	movapd   xmm1, [rsp + nb213nf_two]
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
	mulsd  xmm1, [rsp + nb213nf_c6]
	mulsd  xmm2, [rsp + nb213nf_c12]
	movapd xmm3, xmm2
	subsd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addsd  xmm3, [rsp + nb213nf_Vvdwtot]
	movsd [rsp + nb213nf_Vvdwtot], xmm3

	;# H1 interactions 
	movsd  xmm6, [rsp + nb213nf_rinvH1]
	addsd   xmm6, [rsp + nb213nf_krsqH1]
 	subsd   xmm6, [rsp + nb213nf_crf]
	mulsd   xmm6, [rsp + nb213nf_qqH] ;# vcoul 
	addsd   xmm6, [rsp + nb213nf_vctot]
	movsd [rsp + nb213nf_vctot], xmm6
	
	;# H2 interactions 
	movsd  xmm6, [rsp + nb213nf_rinvH2]
	addsd   xmm6, [rsp + nb213nf_krsqH2]
 	subsd   xmm6, [rsp + nb213nf_crf]
	mulsd   xmm6, [rsp + nb213nf_qqH] ;# vcoul 
	addsd   xmm6, [rsp + nb213nf_vctot]
	movsd [rsp + nb213nf_vctot], xmm6

	;# M interactions 
	movsd  xmm6, [rsp + nb213nf_rinvM]
	addsd   xmm6, [rsp + nb213nf_krsqM]
 	subsd   xmm6, [rsp + nb213nf_crf]
	mulsd   xmm6, [rsp + nb213nf_qqM] ;# vcoul 
	addsd   xmm6, [rsp + nb213nf_vctot]
	movsd [rsp + nb213nf_vctot], xmm6
	
.nb213nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb213nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb213nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb213nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb213nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb213nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb213nf_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
       ;# finish if last 
        mov ecx, [rsp + nb213nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb213nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb213nf_n], esi
        jmp .nb213nf_outer
.nb213nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb213nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb213nf_end
        ;# non-zero, do one more workunit
        jmp   .nb213nf_threadloop
.nb213nf_end:
	mov eax, [rsp + nb213nf_nouter]
	mov ebx, [rsp + nb213nf_ninner]
	mov rcx, [rbp + nb213nf_outeriter]
	mov rdx, [rbp + nb213nf_inneriter]
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
