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





.globl nb_kernel313_x86_64_sse2
.globl _nb_kernel313_x86_64_sse2
nb_kernel313_x86_64_sse2:	
_nb_kernel313_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb313_fshift,           16
.equiv          nb313_gid,              24
.equiv          nb313_pos,              32
.equiv          nb313_faction,          40
.equiv          nb313_charge,           48
.equiv          nb313_p_facel,          56
.equiv          nb313_argkrf,           64
.equiv          nb313_argcrf,           72
.equiv          nb313_Vc,               80
.equiv          nb313_type,             88
.equiv          nb313_p_ntype,          96
.equiv          nb313_vdwparam,         104
.equiv          nb313_Vvdw,             112
.equiv          nb313_p_tabscale,       120
.equiv          nb313_VFtab,            128
.equiv          nb313_invsqrta,         136
.equiv          nb313_dvda,             144
.equiv          nb313_p_gbtabscale,     152
.equiv          nb313_GBtab,            160
.equiv          nb313_p_nthreads,       168
.equiv          nb313_count,            176
.equiv          nb313_mtx,              184
.equiv          nb313_outeriter,        192
.equiv          nb313_inneriter,        200
.equiv          nb313_work,             208
	;# stack offsets for local variables 
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb313_ixO,              0
.equiv          nb313_iyO,              16
.equiv          nb313_izO,              32
.equiv          nb313_ixH1,             48
.equiv          nb313_iyH1,             64
.equiv          nb313_izH1,             80
.equiv          nb313_ixH2,             96
.equiv          nb313_iyH2,             112
.equiv          nb313_izH2,             128
.equiv          nb313_ixM,              144
.equiv          nb313_iyM,              160
.equiv          nb313_izM,              176
.equiv          nb313_iqM,              192
.equiv          nb313_iqH,              208
.equiv          nb313_dxO,              224
.equiv          nb313_dyO,              240
.equiv          nb313_dzO,              256
.equiv          nb313_dxH1,             272
.equiv          nb313_dyH1,             288
.equiv          nb313_dzH1,             304
.equiv          nb313_dxH2,             320
.equiv          nb313_dyH2,             336
.equiv          nb313_dzH2,             352
.equiv          nb313_dxM,              368
.equiv          nb313_dyM,              384
.equiv          nb313_dzM,              400
.equiv          nb313_qqM,              416
.equiv          nb313_qqH,              432
.equiv          nb313_rinvsqO,          448
.equiv          nb313_rinvH1,           464
.equiv          nb313_rinvH2,           480
.equiv          nb313_rinvM,            496
.equiv          nb313_rO,               512
.equiv          nb313_rH1,              528
.equiv          nb313_rH2,              544
.equiv          nb313_rM,               560
.equiv          nb313_tsc,              576
.equiv          nb313_two,              592
.equiv          nb313_c6,               608
.equiv          nb313_c12,              624
.equiv          nb313_six,              640
.equiv          nb313_twelve,           656
.equiv          nb313_vctot,            672
.equiv          nb313_Vvdwtot,          688
.equiv          nb313_fixO,             704
.equiv          nb313_fiyO,             720
.equiv          nb313_fizO,             736
.equiv          nb313_fixH1,            752
.equiv          nb313_fiyH1,            768
.equiv          nb313_fizH1,            784
.equiv          nb313_fixH2,            800
.equiv          nb313_fiyH2,            816
.equiv          nb313_fizH2,            832
.equiv          nb313_fixM,             848
.equiv          nb313_fiyM,             864
.equiv          nb313_fizM,             880
.equiv          nb313_fjx,              896
.equiv          nb313_fjy,              912
.equiv          nb313_fjz,              928
.equiv          nb313_half,             944
.equiv          nb313_three,            960
.equiv          nb313_is3,              976
.equiv          nb313_ii3,              980
.equiv          nb313_nri,              984
.equiv          nb313_iinr,             992
.equiv          nb313_jindex,           1000
.equiv          nb313_jjnr,             1008
.equiv          nb313_shift,            1016
.equiv          nb313_shiftvec,         1024
.equiv          nb313_facel,            1032
.equiv          nb313_innerjjnr,        1040
.equiv          nb313_ntia,             1048
.equiv          nb313_innerk,           1052
.equiv          nb313_n,                1056
.equiv          nb313_nn1,              1060
.equiv          nb313_nouter,           1064
.equiv          nb313_ninner,           1068

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
	sub rsp, 1072		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb313_nouter], eax
	mov [rsp + nb313_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb313_nri], edi
	mov [rsp + nb313_iinr], rsi
	mov [rsp + nb313_jindex], rdx
	mov [rsp + nb313_jjnr], rcx
	mov [rsp + nb313_shift], r8
	mov [rsp + nb313_shiftvec], r9
	mov rsi, [rbp + nb313_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb313_facel], xmm0

	mov rax, [rbp + nb313_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb313_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb313_half], eax
	mov [rsp + nb313_half+4], ebx
	movsd xmm1, [rsp + nb313_half]
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
	movapd [rsp + nb313_half], xmm1
	movapd [rsp + nb313_two], xmm2
	movapd [rsp + nb313_three], xmm3
	movapd [rsp + nb313_six], xmm4
	movapd [rsp + nb313_twelve], xmm5

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb313_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb313_charge]
	movsd xmm3, [rdx + rbx*8 + 8]	
	movsd xmm4, [rdx + rbx*8 + 24]	
	mov rsi, [rbp + nb313_p_facel]
	movsd xmm0, [rsi]
	movsd xmm5, [rsp + nb313_facel]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb313_iqH], xmm3
	movapd [rsp + nb313_iqM], xmm4
	
	mov   rdx, [rbp + nb313_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb313_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb313_ntia], ecx		
.nb313_threadloop:
        mov   rsi, [rbp + nb313_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb313_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb313_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb313_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb313_n], eax
        mov [rsp + nb313_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb313_outerstart
        jmp .nb313_end

.nb313_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb313_nouter]
	mov [rsp + nb313_nouter], ebx

.nb313_outer:
	mov   rax, [rsp + nb313_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb313_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb313_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb313_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb313_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb313_ii3], ebx

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
	movapd [rsp + nb313_ixO], xmm3
	movapd [rsp + nb313_iyO], xmm4
	movapd [rsp + nb313_izO], xmm5
	movapd [rsp + nb313_ixH1], xmm6
	movapd [rsp + nb313_iyH1], xmm7

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
	movapd [rsp + nb313_izH1], xmm6
	movapd [rsp + nb313_ixH2], xmm0
	movapd [rsp + nb313_iyH2], xmm1
	movapd [rsp + nb313_izH2], xmm2
	movapd [rsp + nb313_ixM], xmm3
	movapd [rsp + nb313_iyM], xmm4
	movapd [rsp + nb313_izM], xmm5
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb313_vctot], xmm4
	movapd [rsp + nb313_Vvdwtot], xmm4
	movapd [rsp + nb313_fixO], xmm4
	movapd [rsp + nb313_fiyO], xmm4
	movapd [rsp + nb313_fizO], xmm4
	movapd [rsp + nb313_fixH1], xmm4
	movapd [rsp + nb313_fiyH1], xmm4
	movapd [rsp + nb313_fizH1], xmm4
	movapd [rsp + nb313_fixH2], xmm4
	movapd [rsp + nb313_fiyH2], xmm4
	movapd [rsp + nb313_fizH2], xmm4
	movapd [rsp + nb313_fixM], xmm4
	movapd [rsp + nb313_fiyM], xmm4
	movapd [rsp + nb313_fizM], xmm4
	
	mov   rax, [rsp + nb313_jindex]
	mov   ecx, [rax+rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb313_pos]
	mov   rdi, [rbp + nb313_faction]	
	mov   rax, [rsp + nb313_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb313_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb313_ninner]
	mov   [rsp + nb313_ninner], ecx
	add   edx, 0
	mov   [rsp + nb313_innerk], edx    ;# number of innerloop atoms 
	jge   .nb313_unroll_loop
	jmp   .nb313_checksingle
.nb313_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb313_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]
	
	add qword ptr [rsp + nb313_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb313_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]
	movapd xmm4, xmm3	     
	mulpd  xmm3, [rsp + nb313_iqM]
	mulpd  xmm4, [rsp + nb313_iqH]
	movapd  [rsp + nb313_qqM], xmm3
	movapd  [rsp + nb313_qqH], xmm4	

	mov rsi, [rbp + nb313_type]
	mov r8d, [rsi + rax*4]
	mov r9d, [rsi + rbx*4]
	mov rsi, [rbp + nb313_vdwparam]
	shl r8d, 1	
	shl r9d, 1	
	mov edi, [rsp + nb313_ntia]
	add r8d, edi
	add r9d, edi

	movlpd xmm6, [rsi + r8*8]	;# c6a
	movlpd xmm7, [rsi + r9*8]	;# c6b
	movhpd xmm6, [rsi + r8*8 + 8]	;# c6a c12a 
	movhpd xmm7, [rsi + r9*8 + 8]	;# c6b c12b 
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movapd [rsp + nb313_c6], xmm4
	movapd [rsp + nb313_c12], xmm6

	mov rsi, [rbp + nb313_pos]       ;# base of pos[] 

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
    
    subpd xmm3, [rsp + nb313_ixO]
    subpd xmm4, [rsp + nb313_iyO]
    subpd xmm5, [rsp + nb313_izO]
    
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
    movapd xmm4, [rsp + nb313_two]
    movapd xmm5, xmm3
    mulpd xmm3, xmm6        ;# lu*rsq 
    subpd xmm4, xmm3        ;# 2-lu*rsq 
    mulpd xmm6, xmm4        ;# (new lu) 
        
    movapd xmm4, [rsp + nb313_two]
    mulpd xmm5, xmm6        ;# lu*rsq 
    subpd xmm4, xmm5        ;# 2-lu*rsq 
    mulpd xmm4, xmm6        ;# xmm4=rinvsq 

    movapd xmm3, xmm4       ;# rinvsq
    mulpd  xmm4, xmm4       ;# rinv4
    mulpd  xmm4, xmm3       ;# rinv6
    movapd xmm5, xmm4     
    mulpd  xmm5, xmm5       ;# rinv12
    mulpd  xmm4, [rsp + nb313_c6]
    mulpd  xmm5, [rsp + nb313_c12]
    movapd xmm6, xmm5
    subpd  xmm6, xmm4  ;# Vvdw=vvdw12-vvdw6
    mulpd  xmm4, [rsp + nb313_six]
    mulpd  xmm5, [rsp + nb313_twelve]
    subpd  xmm5, xmm4
    mulpd  xmm3, xmm5   ;# fscal
    
    addpd  xmm6, [rsp + nb313_Vvdwtot]
    movapd [rsp + nb313_Vvdwtot], xmm6
    
    mulpd  xmm13, xmm3 ;# fx
    mulpd  xmm14, xmm3 ;# fy
    mulpd  xmm15, xmm3 ;# fz

    ;# save j force temporarily
    movapd [rsp + nb313_fjx], xmm13
    movapd [rsp + nb313_fjy], xmm14
    movapd [rsp + nb313_fjz], xmm15
    
    ;# increment i O force
    addpd xmm13, [rsp + nb313_fixO]
    addpd xmm14, [rsp + nb313_fiyO]
    addpd xmm15, [rsp + nb313_fizO]
    movapd [rsp + nb313_fixO], xmm13
    movapd [rsp + nb313_fiyO], xmm14
    movapd [rsp + nb313_fizO], xmm15    
    ;# finished O LJ interaction.


    ;# do H1, H2, and M interactions in parallel.
    ;# xmm0-xmm2 still contain j coordinates.                
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subpd xmm0, [rsp + nb313_ixH1]
    subpd xmm1, [rsp + nb313_iyH1]
    subpd xmm2, [rsp + nb313_izH1]
    subpd xmm3, [rsp + nb313_ixH2]
    subpd xmm4, [rsp + nb313_iyH2]
    subpd xmm5, [rsp + nb313_izH2]
    subpd xmm6, [rsp + nb313_ixM]
    subpd xmm7, [rsp + nb313_iyM]
    subpd xmm8, [rsp + nb313_izM]
    
	movapd [rsp + nb313_dxH1], xmm0
	movapd [rsp + nb313_dyH1], xmm1
	movapd [rsp + nb313_dzH1], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb313_dxH2], xmm3
	movapd [rsp + nb313_dyH2], xmm4
	movapd [rsp + nb313_dzH2], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb313_dxM], xmm6
	movapd [rsp + nb313_dyM], xmm7
	movapd [rsp + nb313_dzM], xmm8
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
		
	movapd  xmm9, [rsp + nb313_three]
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

	movapd  xmm15, [rsp + nb313_half]
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
		
	movapd  xmm1, [rsp + nb313_three]
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

	movapd  xmm15, [rsp + nb313_half]
	mulpd   xmm9, xmm15  ;#  rinvH1
	mulpd   xmm10, xmm15 ;#   rinvH2
    mulpd   xmm11, xmm15 ;#   rinvM
	
	movapd  [rsp + nb313_rinvH1], xmm9
	movapd  [rsp + nb313_rinvH2], xmm10
	movapd  [rsp + nb313_rinvM], xmm11
	
	;# interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movapd xmm1, [rsp + nb313_tsc]
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
        
    mov  rsi, [rbp + nb313_VFtab]

    ;# calculate eps
    subpd     xmm0, xmm2
    subpd     xmm3, xmm5
    subpd     xmm6, xmm8

    movapd xmm12, xmm0  ;# epsH1
    movapd xmm13, xmm3  ;# epsH2
    movapd xmm14, xmm6  ;# epsM

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
    movapd xmm12, [rsp + nb313_qqH]
    movapd xmm13, [rsp + nb313_qqM]
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
    addpd  xmm1, [rsp + nb313_vctot]
    addpd  xmm5, xmm9
    addpd  xmm1, xmm5
    movapd [rsp + nb313_vctot], xmm1
    
    movapd xmm10, [rsp + nb313_tsc]
    mulpd  xmm3, xmm10  ;# fscal
    mulpd  xmm7, xmm10
    mulpd  xmm10, xmm11
    
    xorpd xmm4, xmm4
    xorpd xmm8, xmm8
    xorpd xmm11, xmm11
    
    subpd xmm4, xmm3
    subpd xmm8, xmm7
    subpd xmm11, xmm10

    mulpd xmm4, [rsp + nb313_rinvH1]
    mulpd xmm8, [rsp + nb313_rinvH2]
    mulpd xmm11, [rsp + nb313_rinvM]
    
    ;# move j forces to xmm0-xmm2
    mov rdi, [rbp + nb313_faction]
	movlpd xmm0, [rdi + rax*8]
	movlpd xmm1, [rdi + rax*8 + 8]
	movlpd xmm2, [rdi + rax*8 + 16]
	movhpd xmm0, [rdi + rbx*8]
	movhpd xmm1, [rdi + rbx*8 + 8]
	movhpd xmm2, [rdi + rbx*8 + 16]

    movapd xmm3, xmm4
    movapd xmm5, xmm4
    movapd xmm7, xmm8
    movapd xmm9, xmm8
    movapd xmm10, xmm11
    movapd xmm12, xmm11

    ;# add forces from O interaction
    addpd xmm0, [rsp + nb313_fjx]
    addpd xmm1, [rsp + nb313_fjy]
    addpd xmm2, [rsp + nb313_fjz]

	mulpd xmm3, [rsp + nb313_dxH1]
	mulpd xmm4, [rsp + nb313_dyH1]
	mulpd xmm5, [rsp + nb313_dzH1]
	mulpd xmm7, [rsp + nb313_dxH2]
	mulpd xmm8, [rsp + nb313_dyH2]
	mulpd xmm9, [rsp + nb313_dzH2]
	mulpd xmm10, [rsp + nb313_dxM]
	mulpd xmm11, [rsp + nb313_dyM]
	mulpd xmm12, [rsp + nb313_dzM]

    addpd xmm0, xmm3
    addpd xmm1, xmm4
    addpd xmm2, xmm5
    addpd xmm3, [rsp + nb313_fixH1]
    addpd xmm4, [rsp + nb313_fiyH1]
    addpd xmm5, [rsp + nb313_fizH1]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb313_fixH2]
    addpd xmm8, [rsp + nb313_fiyH2]
    addpd xmm9, [rsp + nb313_fizH2]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb313_fixM]
    addpd xmm11, [rsp + nb313_fiyM]
    addpd xmm12, [rsp + nb313_fizM]

    movapd [rsp + nb313_fixH1], xmm3
    movapd [rsp + nb313_fiyH1], xmm4
    movapd [rsp + nb313_fizH1], xmm5
    movapd [rsp + nb313_fixH2], xmm7
    movapd [rsp + nb313_fiyH2], xmm8
    movapd [rsp + nb313_fizH2], xmm9
    movapd [rsp + nb313_fixM], xmm10
    movapd [rsp + nb313_fiyM], xmm11
    movapd [rsp + nb313_fizM], xmm12
   
    ;# store back j forces from xmm0-xmm2
	movlpd [rdi + rax*8], xmm0
	movlpd [rdi + rax*8 + 8], xmm1
	movlpd [rdi + rax*8 + 16], xmm2
	movhpd [rdi + rbx*8], xmm0
	movhpd [rdi + rbx*8 + 8], xmm1
	movhpd [rdi + rbx*8 + 16], xmm2

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb313_innerk],  2
	jl    .nb313_checksingle
	jmp   .nb313_unroll_loop
.nb313_checksingle:	
	mov   edx, [rsp + nb313_innerk]
	and   edx, 1
	jnz   .nb313_dosingle
	jmp   .nb313_updateouterdata
.nb313_dosingle:
	mov   rdx, [rsp + nb313_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]

	mov rsi, [rbp + nb313_charge]    ;# base of charge[] 
	xorpd xmm3, xmm3	
	movlpd xmm3, [rsi + rax*8]
	movapd xmm4, xmm3	     
	mulpd  xmm3, [rsp + nb313_iqM]
	mulpd  xmm4, [rsp + nb313_iqH]

	movapd  [rsp + nb313_qqM], xmm3
	movapd  [rsp + nb313_qqH], xmm4	
	
	mov rsi, [rbp + nb313_type]
	mov r8d, [rsi + rax*4]
	mov rsi, [rbp + nb313_vdwparam]
	shl r8d, 1	
	mov edi, [rsp + nb313_ntia]
	add r8d, edi

	movlpd xmm6, [rsi + r8*8]	;# c6a
	movhpd xmm6, [rsi + r8*8 + 8]	;# c6a c12a 
	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movapd [rsp + nb313_c6], xmm4
	movapd [rsp + nb313_c12], xmm6
	
	mov rsi, [rbp + nb313_pos]       ;# base of pos[] 
	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	
	;# move coordinates to xmm0-xmm2  and xmm4-xmm6
	movlpd xmm4, [rsi + rax*8]
	movlpd xmm5, [rsi + rax*8 + 8]
	movlpd xmm6, [rsi + rax*8 + 16]
    movapd xmm0, xmm4
    movapd xmm1, xmm5
    movapd xmm2, xmm6

	;# calc dr 
	subsd xmm4, [rsp + nb313_ixO]
	subsd xmm5, [rsp + nb313_iyO]
	subsd xmm6, [rsp + nb313_izO]

	;# store dr 
	movapd [rsp + nb313_dxO], xmm4
	movapd [rsp + nb313_dyO], xmm5
	movapd [rsp + nb313_dzO], xmm6
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
	subsd xmm4, [rsp + nb313_ixH1]
	subsd xmm5, [rsp + nb313_iyH1]
	subsd xmm6, [rsp + nb313_izH1]

	;# store dr 
	movapd [rsp + nb313_dxH1], xmm4
	movapd [rsp + nb313_dyH1], xmm5
	movapd [rsp + nb313_dzH1], xmm6
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
	subsd xmm3, [rsp + nb313_ixH2]
	subsd xmm4, [rsp + nb313_iyH2]
	subsd xmm5, [rsp + nb313_izH2]

	;# store dr 
	movapd [rsp + nb313_dxH2], xmm3
	movapd [rsp + nb313_dyH2], xmm4
	movapd [rsp + nb313_dzH2], xmm5
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
	subsd xmm4, [rsp + nb313_ixM]
	subsd xmm3, [rsp + nb313_iyM]
	subsd xmm2, [rsp + nb313_izM]

	;# store dr 
	movapd [rsp + nb313_dxM], xmm4
	movapd [rsp + nb313_dyM], xmm3
	movapd [rsp + nb313_dzM], xmm2

	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# 1/x for O - rsqO is in xmm7
	cvtsd2ss xmm2, xmm7
	movsd   xmm3, xmm7
	rcpps    xmm2, xmm2
	cvtss2sd xmm2, xmm2
	movsd   xmm1, [rsp + nb313_two]
	movsd   xmm0, xmm1
	mulsd   xmm7, xmm2
	subsd   xmm1, xmm7
	mulsd   xmm2, xmm1 ;# iter1 
	mulsd   xmm3, xmm2
	subsd   xmm0, xmm3
	mulsd   xmm0, xmm2 ;# xmm0=rinvsq
	movsd  [rsp + nb313_rinvsqO], xmm0

	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm0, [rsp + nb313_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm0, [rsp + nb313_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm6
	movapd xmm3, xmm0
	mulsd xmm0, xmm0	;# lu*lu 
	mulsd xmm2, xmm0	;# rsq*lu*lu 
	movapd xmm0, [rsp + nb313_three]
	subsd xmm0, xmm2	;# 3-rsq*lu*lu 
	mulsd xmm0, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm0, [rsp + nb313_half] ;# rinv 
	movapd [rsp + nb313_rinvH1], xmm0	;# rinvH1 
	mulsd  xmm6, xmm0
	movapd [rsp + nb313_rH1], xmm6	;# rH1 
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm0, [rsp + nb313_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm0, [rsp + nb313_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm5
	movapd xmm3, xmm0
	mulsd xmm0, xmm0	;# lu*lu 
	mulsd xmm2, xmm0	;# rsq*lu*lu 
	movapd xmm0, [rsp + nb313_three]
	subsd xmm0, xmm2	;# 3-rsq*lu*lu 
	mulsd xmm0, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm0, [rsp + nb313_half] ;# rinv 
	movapd [rsp + nb313_rinvH2], xmm0 ;# rinv 
	mulsd xmm5, xmm0
	movapd [rsp + nb313_rH2], xmm5 ;# r 

	;# rsqM - seed in xmm2 
	cvtsd2ss xmm2, xmm4	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm0, [rsp + nb313_three]
	mulsd   xmm2, xmm4	;# rsq*lu*lu 
	subsd   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm0, [rsp + nb313_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm4
	movapd xmm3, xmm0
	mulsd xmm0, xmm0	;# lu*lu 
	mulsd xmm2, xmm0	;# rsq*lu*lu 
	movapd xmm0, [rsp + nb313_three]
	subsd xmm0, xmm2	;# 3-rsq*lu*lu 
	mulsd xmm0, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm0, [rsp + nb313_half] ;# rinv 
	movapd [rsp + nb313_rinvM], xmm0 ;# rinv 
	mulsd xmm4, xmm0
	movapd [rsp + nb313_rM], xmm4 ;# r 

	;# do O interactions 
	movapd  xmm0, [rsp + nb313_rinvsqO]
	movapd  xmm1, xmm0	
	mulsd   xmm1, xmm1 ;# rinv4
	mulsd   xmm1, xmm0 ;#rinvsix
	movapd  xmm2, xmm1
	mulsd	xmm2, xmm2 ;# rinvtwelve
	mulsd  xmm1, [rsp + nb313_c6]
	mulsd  xmm2, [rsp + nb313_c12]
	movapd xmm3, xmm2
	subsd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addsd  xmm3, [rsp + nb313_Vvdwtot]
	mulsd  xmm1, [rsp + nb313_six]
	mulsd  xmm2, [rsp + nb313_twelve]
	subsd  xmm2, xmm1
	mulsd  xmm2, xmm0
	movapd xmm4, xmm2 ;# total fsO 
	movsd [rsp + nb313_Vvdwtot], xmm3

	movapd xmm0, [rsp + nb313_dxO]
	movapd xmm1, [rsp + nb313_dyO]
	movapd xmm2, [rsp + nb313_dzO]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4
	;# tx in xmm0-xmm2 

	;# update O forces 
	movapd xmm3, [rsp + nb313_fixO]
	movapd xmm4, [rsp + nb313_fiyO]
	movapd xmm7, [rsp + nb313_fizO]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [rsp + nb313_fixO], xmm3
	movlpd [rsp + nb313_fiyO], xmm4
	movlpd [rsp + nb313_fizO], xmm7
	;# update j forces with water O 
	movlpd [rsp + nb313_fjx], xmm0
	movlpd [rsp + nb313_fjy], xmm1
	movlpd [rsp + nb313_fjz], xmm2

	;# Done with O interactions - now H1! 
	movapd xmm7, [rsp + nb313_rH1]
	mulpd xmm7, [rsp + nb313_tsc]
	cvttsd2si r8d, xmm7	;# mm6 = lu idx 
	cvtsi2sd xmm6, r8d
	subpd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	shl r8d, 2		;# idx *= 4 
	mov  rsi, [rbp + nb313_VFtab]

	movapd xmm4, [rsi + r8*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + r8*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [rsp + nb313_two]	;# two*Heps2 
	movapd xmm3, [rsp + nb313_qqH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 
    ;# increment vcoul 
	xorpd  xmm4, xmm4
    addsd  xmm5, [rsp + nb313_vctot]
	mulsd  xmm3, [rsp + nb313_rinvH1]
    movlpd [rsp + nb313_vctot], xmm5 
	mulsd  xmm3, [rsp + nb313_tsc]
	subsd xmm4, xmm3

	movapd xmm0, [rsp + nb313_dxH1]
	movapd xmm1, [rsp + nb313_dyH1]
	movapd xmm2, [rsp + nb313_dzH1]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [rsp + nb313_fixH1]
	movapd xmm4, [rsp + nb313_fiyH1]
	movapd xmm7, [rsp + nb313_fizH1]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [rsp + nb313_fixH1], xmm3
	movlpd [rsp + nb313_fiyH1], xmm4
	movlpd [rsp + nb313_fizH1], xmm7
	;# update j forces with water H1 
	addsd  xmm0, [rsp + nb313_fjx]
	addsd  xmm1, [rsp + nb313_fjy]
	addsd  xmm2, [rsp + nb313_fjz]
	movlpd [rsp + nb313_fjx], xmm0
	movlpd [rsp + nb313_fjy], xmm1
	movlpd [rsp + nb313_fjz], xmm2

	;#  H2 interactions 
	movapd xmm7, [rsp + nb313_rH2]
	mulsd   xmm7, [rsp + nb313_tsc]
	cvttsd2si r8d, xmm7	;# mm6 = lu idx 
	cvtsi2sd xmm6, r8d
	subsd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl r8d, 2		;# idx *= 4 
	mov  rsi, [rbp + nb313_VFtab]

	movapd xmm4, [rsi + r8*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + r8*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [rsp + nb313_two]	;# two*Heps2 
	movapd xmm3, [rsp + nb313_qqH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm3 fijC 
    	;# increment vcoul 
	xorpd  xmm4, xmm4
    	addsd  xmm5, [rsp + nb313_vctot]
	mulsd  xmm3, [rsp + nb313_rinvH2]
    	movlpd [rsp + nb313_vctot], xmm5 
	mulsd  xmm3, [rsp + nb313_tsc]
	subsd  xmm4, xmm3

	movapd xmm0, [rsp + nb313_dxH2]
	movapd xmm1, [rsp + nb313_dyH2]
	movapd xmm2, [rsp + nb313_dzH2]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4
	
	;# update H2 forces 
	movapd xmm3, [rsp + nb313_fixH2]
	movapd xmm4, [rsp + nb313_fiyH2]
	movapd xmm7, [rsp + nb313_fizH2]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [rsp + nb313_fixH2], xmm3
	movlpd [rsp + nb313_fiyH2], xmm4
	movlpd [rsp + nb313_fizH2], xmm7
	;# update j forces with water H1 
	addsd  xmm0, [rsp + nb313_fjx]
	addsd  xmm1, [rsp + nb313_fjy]
	addsd  xmm2, [rsp + nb313_fjz]
	movlpd [rsp + nb313_fjx], xmm0
	movlpd [rsp + nb313_fjy], xmm1
	movlpd [rsp + nb313_fjz], xmm2

	;# M interactions 
	movapd xmm7, [rsp + nb313_rM]
	mulsd   xmm7, [rsp + nb313_tsc]
	cvttsd2si r8d, xmm7	;# mm6 = lu idx 
	cvtsi2sd xmm6, r8d
	subsd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl r8d, 2		;# idx *= 4 
	mov  rsi, [rbp + nb313_VFtab]

	movapd xmm4, [rsi + r8*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + r8*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [rsp + nb313_two]	;# two*Heps2 
	movapd xmm3, [rsp + nb313_qqM]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm3 fijC 
    	;# increment vcoul 
	xorpd  xmm4, xmm4
    	addsd  xmm5, [rsp + nb313_vctot]
	mulsd  xmm3, [rsp + nb313_rinvM]
    	movlpd [rsp + nb313_vctot], xmm5 
	mulsd  xmm3, [rsp + nb313_tsc]
	subsd  xmm4, xmm3

	movapd xmm0, [rsp + nb313_dxM]
	movapd xmm1, [rsp + nb313_dyM]
	movapd xmm2, [rsp + nb313_dzM]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4
	
	;# update M forces 
	movapd xmm3, [rsp + nb313_fixM]
	movapd xmm4, [rsp + nb313_fiyM]
	movapd xmm7, [rsp + nb313_fizM]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [rsp + nb313_fixM], xmm3
	movlpd [rsp + nb313_fiyM], xmm4
	movlpd [rsp + nb313_fizM], xmm7

	mov rdi, [rbp + nb313_faction]
	;# update j forces 
	addsd  xmm0, [rsp + nb313_fjx]
	addsd  xmm1, [rsp + nb313_fjy]
	addsd  xmm2, [rsp + nb313_fjz]

	;# the fj's - start by accumulating forces from memory 
	movlpd xmm3, [rdi + rax*8]
	movlpd xmm4, [rdi + rax*8 + 8]
	movlpd xmm5, [rdi + rax*8 + 16]
	addsd xmm3, xmm0
	addsd xmm4, xmm1
	addsd xmm5, xmm2
	movlpd [rdi + rax*8], xmm3
	movlpd [rdi + rax*8 + 8], xmm4
	movlpd [rdi + rax*8 + 16], xmm5

.nb313_updateouterdata:
	mov   ecx, [rsp + nb313_ii3]
	mov   rdi, [rbp + nb313_faction]
	mov   rsi, [rbp + nb313_fshift]
	mov   edx, [rsp + nb313_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb313_fixO]
	movapd xmm1, [rsp + nb313_fiyO]
	movapd xmm2, [rsp + nb313_fizO]

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
	movapd xmm0, [rsp + nb313_fixH1]
	movapd xmm1, [rsp + nb313_fiyH1]
	movapd xmm2, [rsp + nb313_fizH1]

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
	movapd xmm0, [rsp + nb313_fixH2]
	movapd xmm1, [rsp + nb313_fiyH2]
	movapd xmm2, [rsp + nb313_fizH2]

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
	movapd xmm0, [rsp + nb313_fixM]
	movapd xmm1, [rsp + nb313_fiyM]
	movapd xmm2, [rsp + nb313_fizM]

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
	mov esi, [rsp + nb313_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb313_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb313_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb313_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb313_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb313_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb313_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb313_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb313_n], esi
        jmp .nb313_outer
.nb313_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb313_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb313_end
        ;# non-zero, do one more workunit
        jmp   .nb313_threadloop
.nb313_end:
	mov eax, [rsp + nb313_nouter]
	mov ebx, [rsp + nb313_ninner]
	mov rcx, [rbp + nb313_outeriter]
	mov rdx, [rbp + nb313_inneriter]
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
	



.globl nb_kernel313nf_x86_64_sse2
.globl _nb_kernel313nf_x86_64_sse2
nb_kernel313nf_x86_64_sse2:	
_nb_kernel313nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb313nf_fshift,         16
.equiv          nb313nf_gid,            24
.equiv          nb313nf_pos,            32
.equiv          nb313nf_faction,        40
.equiv          nb313nf_charge,         48
.equiv          nb313nf_p_facel,        56
.equiv          nb313nf_argkrf,         64
.equiv          nb313nf_argcrf,         72
.equiv          nb313nf_Vc,             80
.equiv          nb313nf_type,           88
.equiv          nb313nf_p_ntype,        96
.equiv          nb313nf_vdwparam,       104
.equiv          nb313nf_Vvdw,           112
.equiv          nb313nf_p_tabscale,     120
.equiv          nb313nf_VFtab,          128
.equiv          nb313nf_invsqrta,       136
.equiv          nb313nf_dvda,           144
.equiv          nb313nf_p_gbtabscale,   152
.equiv          nb313nf_GBtab,          160
.equiv          nb313nf_p_nthreads,     168
.equiv          nb313nf_count,          176
.equiv          nb313nf_mtx,            184
.equiv          nb313nf_outeriter,      192
.equiv          nb313nf_inneriter,      200
.equiv          nb313nf_work,           208
	;# stack offsets for local variables 
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb313nf_ixO,            0
.equiv          nb313nf_iyO,            16
.equiv          nb313nf_izO,            32
.equiv          nb313nf_ixH1,           48
.equiv          nb313nf_iyH1,           64
.equiv          nb313nf_izH1,           80
.equiv          nb313nf_ixH2,           96
.equiv          nb313nf_iyH2,           112
.equiv          nb313nf_izH2,           128
.equiv          nb313nf_ixM,            144
.equiv          nb313nf_iyM,            160
.equiv          nb313nf_izM,            176
.equiv          nb313nf_iqM,            192
.equiv          nb313nf_iqH,            208
.equiv          nb313nf_qqM,            224
.equiv          nb313nf_qqH,            240
.equiv          nb313nf_rinvsqO,        256
.equiv          nb313nf_rinvH1,         272
.equiv          nb313nf_rinvH2,         288
.equiv          nb313nf_rinvM,          304
.equiv          nb313nf_rO,             320
.equiv          nb313nf_rH1,            336
.equiv          nb313nf_rH2,            352
.equiv          nb313nf_rM,             368
.equiv          nb313nf_tsc,            384
.equiv          nb313nf_two,            400
.equiv          nb313nf_c6,             416
.equiv          nb313nf_c12,            432
.equiv          nb313nf_vctot,          448
.equiv          nb313nf_Vvdwtot,        464
.equiv          nb313nf_half,           480
.equiv          nb313nf_three,          496
.equiv          nb313nf_is3,            512
.equiv          nb313nf_ii3,            516
.equiv          nb313nf_nri,            520
.equiv          nb313nf_iinr,           528
.equiv          nb313nf_jindex,         536
.equiv          nb313nf_jjnr,           544
.equiv          nb313nf_shift,          552
.equiv          nb313nf_shiftvec,       560
.equiv          nb313nf_facel,          568
.equiv          nb313nf_innerjjnr,      576
.equiv          nb313nf_ntia,           584
.equiv          nb313nf_innerk,         588
.equiv          nb313nf_n,              592
.equiv          nb313nf_nn1,            596
.equiv          nb313nf_nouter,         600
.equiv          nb313nf_ninner,         604

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
	sub rsp, 608		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb313nf_nouter], eax
	mov [rsp + nb313nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb313nf_nri], edi
	mov [rsp + nb313nf_iinr], rsi
	mov [rsp + nb313nf_jindex], rdx
	mov [rsp + nb313nf_jjnr], rcx
	mov [rsp + nb313nf_shift], r8
	mov [rsp + nb313nf_shiftvec], r9
	mov rsi, [rbp + nb313nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb313nf_facel], xmm0

	mov rax, [rbp + nb313nf_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb313nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb313nf_half], eax
	mov [rsp + nb313nf_half+4], ebx
	movsd xmm1, [rsp + nb313nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb313nf_half], xmm1
	movapd [rsp + nb313nf_two], xmm2
	movapd [rsp + nb313nf_three], xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb313nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb313nf_charge]
	movsd xmm3, [rdx + rbx*8 + 8]	
	movsd xmm4, [rdx + rbx*8 + 24]	
	mov rsi, [rbp + nb313nf_p_facel]
	movsd xmm0, [rsi]
	movsd xmm5, [rsp + nb313nf_facel]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb313nf_iqH], xmm3
	movapd [rsp + nb313nf_iqM], xmm4
	
	mov   rdx, [rbp + nb313nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb313nf_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb313nf_ntia], ecx		
.nb313nf_threadloop:
        mov   rsi, [rbp + nb313nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb313nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb313nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb313nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb313nf_n], eax
        mov [rsp + nb313nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb313nf_outerstart
        jmp .nb313nf_end

.nb313nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb313nf_nouter]
	mov [rsp + nb313nf_nouter], ebx

.nb313nf_outer:
	mov   rax, [rsp + nb313nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb313nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb313nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb313nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb313nf_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb313nf_ii3], ebx

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
	movapd [rsp + nb313nf_ixO], xmm3
	movapd [rsp + nb313nf_iyO], xmm4
	movapd [rsp + nb313nf_izO], xmm5
	movapd [rsp + nb313nf_ixH1], xmm6
	movapd [rsp + nb313nf_iyH1], xmm7

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
	movapd [rsp + nb313nf_izH1], xmm6
	movapd [rsp + nb313nf_ixH2], xmm0
	movapd [rsp + nb313nf_iyH2], xmm1
	movapd [rsp + nb313nf_izH2], xmm2
	movapd [rsp + nb313nf_ixM], xmm3
	movapd [rsp + nb313nf_iyM], xmm4
	movapd [rsp + nb313nf_izM], xmm5
	
	;# clear vctot
	xorpd xmm4, xmm4
	movapd [rsp + nb313nf_vctot], xmm4
	movapd [rsp + nb313nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb313nf_jindex]
	mov   ecx, [rax+rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb313nf_pos]
	mov   rdi, [rbp + nb313nf_faction]	
	mov   rax, [rsp + nb313nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb313nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb313nf_ninner]
	mov   [rsp + nb313nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb313nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb313nf_unroll_loop
	jmp   .nb313nf_checksingle
.nb313nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb313nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]
	
	add qword ptr [rsp + nb313nf_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb313nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]
	movapd xmm4, xmm3	     
	mulpd  xmm3, [rsp + nb313nf_iqM]
	mulpd  xmm4, [rsp + nb313nf_iqH]
	movapd  [rsp + nb313nf_qqM], xmm3
	movapd  [rsp + nb313nf_qqH], xmm4	

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	mov rsi, [rbp + nb313nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov rsi, [rbp + nb313nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	mov edi, [rsp + nb313nf_ntia]
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
	movapd [rsp + nb313nf_c6], xmm4
	movapd [rsp + nb313nf_c12], xmm6

	mov rsi, [rbp + nb313nf_pos]       ;# base of pos[] 

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
	movapd xmm4, [rsp + nb313nf_ixO]
	movapd xmm5, [rsp + nb313nf_iyO]
	movapd xmm6, [rsp + nb313nf_izO]

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
	movapd xmm4, [rsp + nb313nf_ixH1]
	movapd xmm5, [rsp + nb313nf_iyH1]
	movapd xmm6, [rsp + nb313nf_izH1]

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
	movapd xmm3, [rsp + nb313nf_ixH2]
	movapd xmm4, [rsp + nb313nf_iyH2]
	movapd xmm5, [rsp + nb313nf_izH2]

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
	movapd xmm3, [rsp + nb313nf_iyM]
	movapd xmm4, [rsp + nb313nf_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [rsp + nb313nf_ixM]
	subpd  xmm2, xmm0	

	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# 1/x for O - rsqO is in xmm7
	cvtpd2ps xmm2, xmm7
	movapd   xmm3, xmm7
	rcpps    xmm2, xmm2
	cvtps2pd xmm2, xmm2
	movapd   xmm1, [rsp + nb313nf_two]
	movapd   xmm0, xmm1
	mulpd   xmm7, xmm2
	subpd   xmm1, xmm7
	mulpd   xmm2, xmm1 ;# iter1 
	mulpd   xmm3, xmm2
	subpd   xmm0, xmm3
	mulpd   xmm0, xmm2 ;# xmm0=rinvsq
	movapd  [rsp + nb313nf_rinvsqO], xmm0
	
	;# rsqH1 - seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm0, [rsp + nb313nf_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm0, [rsp + nb313nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm6
	movapd xmm3, xmm0
	mulpd xmm0, xmm0	;# lu*lu 
	mulpd xmm2, xmm0	;# rsq*lu*lu 
	movapd xmm0, [rsp + nb313nf_three]
	subpd xmm0, xmm2	;# 3-rsq*lu*lu 
	mulpd xmm0, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm0, [rsp + nb313nf_half] ;# rinv 
	movapd [rsp + nb313nf_rinvH1], xmm0	;# rinvH1 
	mulpd  xmm6, xmm0
	movapd [rsp + nb313nf_rH1], xmm6	;# rH1 
	
	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm0, [rsp + nb313nf_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm0, [rsp + nb313nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm5
	movapd xmm3, xmm0
	mulpd xmm0, xmm0	;# lu*lu 
	mulpd xmm2, xmm0	;# rsq*lu*lu 
	movapd xmm0, [rsp + nb313nf_three]
	subpd xmm0, xmm2	;# 3-rsq*lu*lu 
	mulpd xmm0, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm0, [rsp + nb313nf_half] ;# rinv 
	movapd [rsp + nb313nf_rinvH2], xmm0 ;# rinv 
	mulpd xmm5, xmm0
	movapd [rsp + nb313nf_rH2], xmm5 ;# r 

	;# rsqM - seed in xmm2 
	cvtpd2ps xmm2, xmm4	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm0, [rsp + nb313nf_three]
	mulpd   xmm2, xmm4	;# rsq*lu*lu 
	subpd   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm0, [rsp + nb313nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm4
	movapd xmm3, xmm0
	mulpd xmm0, xmm0	;# lu*lu 
	mulpd xmm2, xmm0	;# rsq*lu*lu 
	movapd xmm0, [rsp + nb313nf_three]
	subpd xmm0, xmm2	;# 3-rsq*lu*lu 
	mulpd xmm0, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm0, [rsp + nb313nf_half] ;# rinv 
	movapd [rsp + nb313nf_rinvM], xmm0 ;# rinv 
	mulpd xmm4, xmm0
	movapd [rsp + nb313nf_rM], xmm4 ;# r 

	;# do O interactions 
	movapd xmm0, [rsp + nb313nf_rinvsqO]
	movapd  xmm1, xmm0	
	mulpd   xmm1, xmm1 ;# rinv4
	mulpd   xmm1, xmm0 ;#rinvsix
	movapd  xmm2, xmm1
	mulpd	xmm2, xmm2 ;# rinvtwelve
	mulpd  xmm1, [rsp + nb313nf_c6]
	mulpd  xmm2, [rsp + nb313nf_c12]
	movapd xmm3, xmm2
	subpd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addpd  xmm3, [rsp + nb313nf_Vvdwtot]
	movapd [rsp + nb313nf_Vvdwtot], xmm3

	;# Done with O interactions - now H1! 
	movapd xmm7, [rsp + nb313nf_rH1]
	mulpd xmm7, [rsp + nb313nf_tsc]
	cvttpd2pi mm6, xmm7	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	movd mm0, eax	
	movd mm1, ebx

	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb313nf_VFtab]
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
	movapd xmm3, [rsp + nb313nf_qqH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    	;# increment vcoul 
    	addpd  xmm5, [rsp + nb313nf_vctot]
    	movapd [rsp + nb313nf_vctot], xmm5 

	;# H2 interactions 
	movapd xmm7, [rsp + nb313nf_rH2]
	mulpd   xmm7, [rsp + nb313nf_tsc]
	cvttpd2pi mm6, xmm7	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb313nf_VFtab]
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
	movapd xmm3, [rsp + nb313nf_qqH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    	;# increment vcoul 
    	addpd  xmm5, [rsp + nb313nf_vctot]
    	movapd [rsp + nb313nf_vctot], xmm5 

	;# M interactions 
	movapd xmm7, [rsp + nb313nf_rM]
	mulpd   xmm7, [rsp + nb313nf_tsc]
	cvttpd2pi mm6, xmm7	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	movd mm0, eax	
	movd mm1, ebx
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb313nf_VFtab]
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
	movapd xmm3, [rsp + nb313nf_qqM]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    	;# increment vcoul 
    	addpd  xmm5, [rsp + nb313nf_vctot]
    	movapd [rsp + nb313nf_vctot], xmm5 

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb313nf_innerk],  2
	jl    .nb313nf_checksingle
	jmp   .nb313nf_unroll_loop
.nb313nf_checksingle:	
	mov   edx, [rsp + nb313nf_innerk]
	and   edx, 1
	jnz   .nb313nf_dosingle
	jmp   .nb313nf_updateouterdata
.nb313nf_dosingle:
	mov   rdx, [rsp + nb313nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]

	mov rsi, [rbp + nb313nf_charge]    ;# base of charge[] 
	xorpd xmm3, xmm3	
	movlpd xmm3, [rsi + rax*8]
	movapd xmm4, xmm3	     
	mulpd  xmm3, [rsp + nb313nf_iqM]
	mulpd  xmm4, [rsp + nb313nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movapd  [rsp + nb313nf_qqM], xmm3
	movapd  [rsp + nb313nf_qqH], xmm4	
	
	mov rsi, [rbp + nb313nf_type]
	mov eax, [rsi + rax*4]
	mov rsi, [rbp + nb313nf_vdwparam]
	shl eax, 1	
	mov edi, [rsp + nb313nf_ntia]
	add eax, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 
	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movapd [rsp + nb313nf_c6], xmm4
	movapd [rsp + nb313nf_c12], xmm6
	
	mov rsi, [rbp + nb313nf_pos]       ;# base of pos[] 
	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	
	;# move coords to xmm0-xmm2 
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [rsp + nb313nf_ixO]
	movapd xmm5, [rsp + nb313nf_iyO]
	movapd xmm6, [rsp + nb313nf_izO]

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
	movapd xmm4, [rsp + nb313nf_ixH1]
	movapd xmm5, [rsp + nb313nf_iyH1]
	movapd xmm6, [rsp + nb313nf_izH1]

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
	movapd xmm3, [rsp + nb313nf_ixH2]
	movapd xmm4, [rsp + nb313nf_iyH2]
	movapd xmm5, [rsp + nb313nf_izH2]

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
	movapd xmm3, [rsp + nb313nf_iyM]
	movapd xmm4, [rsp + nb313nf_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [rsp + nb313nf_ixM]
	subpd  xmm2, xmm0	

	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# 1/x for O - rsqO is in xmm7
	cvtsd2ss xmm2, xmm7
	movsd   xmm3, xmm7
	rcpps    xmm2, xmm2
	cvtss2sd xmm2, xmm2
	movsd   xmm1, [rsp + nb313nf_two]
	movsd   xmm0, xmm1
	mulsd   xmm7, xmm2
	subsd   xmm1, xmm7
	mulsd   xmm2, xmm1 ;# iter1 
	mulsd   xmm3, xmm2
	subsd   xmm0, xmm3
	mulsd   xmm0, xmm2 ;# xmm0=rinvsq
	movsd  [rsp + nb313nf_rinvsqO], xmm0

	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm0, [rsp + nb313nf_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm0, [rsp + nb313nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm6
	movapd xmm3, xmm0
	mulsd xmm0, xmm0	;# lu*lu 
	mulsd xmm2, xmm0	;# rsq*lu*lu 
	movapd xmm0, [rsp + nb313nf_three]
	subsd xmm0, xmm2	;# 3-rsq*lu*lu 
	mulsd xmm0, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm0, [rsp + nb313nf_half] ;# rinv 
	movapd [rsp + nb313nf_rinvH1], xmm0	;# rinvH1 
	mulsd  xmm6, xmm0
	movapd [rsp + nb313nf_rH1], xmm6	;# rH1 
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm0, [rsp + nb313nf_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm0, [rsp + nb313nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm5
	movapd xmm3, xmm0
	mulsd xmm0, xmm0	;# lu*lu 
	mulsd xmm2, xmm0	;# rsq*lu*lu 
	movapd xmm0, [rsp + nb313nf_three]
	subsd xmm0, xmm2	;# 3-rsq*lu*lu 
	mulsd xmm0, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm0, [rsp + nb313nf_half] ;# rinv 
	movapd [rsp + nb313nf_rinvH2], xmm0 ;# rinv 
	mulsd xmm5, xmm0
	movapd [rsp + nb313nf_rH2], xmm5 ;# r 

	;# rsqM - seed in xmm2 
	cvtsd2ss xmm2, xmm4	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm0, [rsp + nb313nf_three]
	mulsd   xmm2, xmm4	;# rsq*lu*lu 
	subsd   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm0, [rsp + nb313nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm4
	movapd xmm3, xmm0
	mulsd xmm0, xmm0	;# lu*lu 
	mulsd xmm2, xmm0	;# rsq*lu*lu 
	movapd xmm0, [rsp + nb313nf_three]
	subsd xmm0, xmm2	;# 3-rsq*lu*lu 
	mulsd xmm0, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm0, [rsp + nb313nf_half] ;# rinv 
	movapd [rsp + nb313nf_rinvM], xmm0 ;# rinv 
	mulsd xmm4, xmm0
	movapd [rsp + nb313nf_rM], xmm4 ;# r 

	;# do O interactions 
	movapd  xmm0, [rsp + nb313nf_rinvsqO]
	movapd  xmm1, xmm0	
	mulsd   xmm1, xmm1 ;# rinv4
	mulsd   xmm1, xmm0 ;#rinvsix
	movapd  xmm2, xmm1
	mulsd	xmm2, xmm2 ;# rinvtwelve
	mulsd  xmm1, [rsp + nb313nf_c6]
	mulsd  xmm2, [rsp + nb313nf_c12]
	movapd xmm3, xmm2
	subsd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addsd  xmm3, [rsp + nb313nf_Vvdwtot]
	movsd [rsp + nb313nf_Vvdwtot], xmm3

	;# Done with O interactions - now H1! 
	movapd xmm7, [rsp + nb313nf_rH1]
	mulpd xmm7, [rsp + nb313nf_tsc]
	cvttsd2si eax, xmm7	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subpd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb313nf_VFtab]

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb313nf_qqH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    	;# increment vcoul 
    	addsd  xmm5, [rsp + nb313nf_vctot]
    	movlpd [rsp + nb313nf_vctot], xmm5 

	;#  H2 interactions 
	movapd xmm7, [rsp + nb313nf_rH2]
	mulsd   xmm7, [rsp + nb313nf_tsc]
	cvttsd2si eax, xmm7	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb313nf_VFtab]

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb313nf_qqH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    	;# increment vcoul 
    	addsd  xmm5, [rsp + nb313nf_vctot]
    	movlpd [rsp + nb313nf_vctot], xmm5 

	;# M interactions 
	movapd xmm7, [rsp + nb313nf_rM]
	mulsd   xmm7, [rsp + nb313nf_tsc]
	cvttsd2si eax, xmm7	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb313nf_VFtab]

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb313nf_qqM]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    	;# increment vcoul 
    	addsd  xmm5, [rsp + nb313nf_vctot]
    	movlpd [rsp + nb313nf_vctot], xmm5 

.nb313nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb313nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb313nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb313nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb313nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb313nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb313nf_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb313nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb313nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb313nf_n], esi
        jmp .nb313nf_outer
.nb313nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb313nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb313nf_end
        ;# non-zero, do one more workunit
        jmp   .nb313nf_threadloop
.nb313nf_end:
	mov eax, [rsp + nb313nf_nouter]
	mov ebx, [rsp + nb313nf_ninner]
	mov rcx, [rbp + nb313nf_outeriter]
	mov rdx, [rbp + nb313nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 608
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
