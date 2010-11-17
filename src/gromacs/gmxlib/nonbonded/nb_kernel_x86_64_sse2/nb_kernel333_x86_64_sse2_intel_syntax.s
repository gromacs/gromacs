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





.globl nb_kernel333_x86_64_sse2
.globl _nb_kernel333_x86_64_sse2
nb_kernel333_x86_64_sse2:	
_nb_kernel333_x86_64_sse2:	
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
	;# bottom of stack is cache-aligned for sse2 use 
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
.equiv          nb333_half,             912
.equiv          nb333_three,            928
.equiv          nb333_is3,              944
.equiv          nb333_ii3,              948
.equiv          nb333_nri,              952
.equiv          nb333_iinr,             960
.equiv          nb333_jindex,           968
.equiv          nb333_jjnr,             976
.equiv          nb333_shift,            984
.equiv          nb333_shiftvec,         992
.equiv          nb333_facel,            1000
.equiv          nb333_innerjjnr,        1008
.equiv          nb333_ntia,             1016
.equiv          nb333_innerk,           1020
.equiv          nb333_n,                1024
.equiv          nb333_nn1,              1028
.equiv          nb333_nouter,           1032
.equiv          nb333_ninner,           1036

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
	sub rsp, 1040		;# local variable stack space (n*16+8)

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
	movsd xmm0, [rsi]
	movsd [rsp + nb333_facel], xmm0

	mov rax, [rbp + nb333_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb333_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb333_half], eax
	mov [rsp + nb333_half+4], ebx
	movsd xmm1, [rsp + nb333_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb333_half], xmm1
	movapd [rsp + nb333_two], xmm2
	movapd [rsp + nb333_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb333_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb333_charge]
	movsd xmm3, [rdx + rbx*8 + 8]	
	movsd xmm4, [rdx + rbx*8 + 24]	
	mov rsi, [rbp + nb333_p_facel]
	movsd xmm0, [rsi]
	movsd xmm5, [rsp + nb333_facel]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb333_iqH], xmm3
	movapd [rsp + nb333_iqM], xmm4
	
	mov   rdx, [rbp + nb333_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb333_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
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
	mov   rax, [rsp + nb333_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb333_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb333_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb333_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb333_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb333_ii3], ebx

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
	movapd [rsp + nb333_ixO], xmm3
	movapd [rsp + nb333_iyO], xmm4
	movapd [rsp + nb333_izO], xmm5
	movapd [rsp + nb333_ixH1], xmm6
	movapd [rsp + nb333_iyH1], xmm7

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
	movapd [rsp + nb333_izH1], xmm6
	movapd [rsp + nb333_ixH2], xmm0
	movapd [rsp + nb333_iyH2], xmm1
	movapd [rsp + nb333_izH2], xmm2
	movapd [rsp + nb333_ixM], xmm3
	movapd [rsp + nb333_iyM], xmm4
	movapd [rsp + nb333_izM], xmm5
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb333_vctot], xmm4
	movapd [rsp + nb333_Vvdwtot], xmm4
	movapd [rsp + nb333_fixO], xmm4
	movapd [rsp + nb333_fiyO], xmm4
	movapd [rsp + nb333_fizO], xmm4
	movapd [rsp + nb333_fixH1], xmm4
	movapd [rsp + nb333_fiyH1], xmm4
	movapd [rsp + nb333_fizH1], xmm4
	movapd [rsp + nb333_fixH2], xmm4
	movapd [rsp + nb333_fiyH2], xmm4
	movapd [rsp + nb333_fizH2], xmm4
	movapd [rsp + nb333_fixM], xmm4
	movapd [rsp + nb333_fiyM], xmm4
	movapd [rsp + nb333_fizM], xmm4
	
	mov   rax, [rsp + nb333_jindex]
	mov   ecx, [rax+rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb333_pos]
	mov   rdi, [rbp + nb333_faction]	
	mov   rax, [rsp + nb333_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb333_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb333_ninner]
	mov   [rsp + nb333_ninner], ecx
	add   edx, 0
	mov   [rsp + nb333_innerk], edx    ;# number of innerloop atoms 
	jge   .nb333_unroll_loop
	jmp   .nb333_checksingle
.nb333_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb333_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]
	
	add qword ptr [rsp + nb333_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb333_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]
	movapd xmm4, xmm3	     
	mulpd  xmm3, [rsp + nb333_iqM]
	mulpd  xmm4, [rsp + nb333_iqH]
	movapd  [rsp + nb333_qqM], xmm3
	movapd  [rsp + nb333_qqH], xmm4	

	mov rsi, [rbp + nb333_type]
	mov r8d, [rsi + rax*4]
	mov r9d, [rsi + rbx*4]
	mov rsi, [rbp + nb333_vdwparam]
	shl r8d, 1	
	shl r9d, 1	
	mov edi, [rsp + nb333_ntia]
	add r8d, edi
	add r9d, edi

	movlpd xmm6, [rsi + r8*8]	;# c6a
	movlpd xmm7, [rsi + r9*8]	;# c6b
	movhpd xmm6, [rsi + r8*8 + 8]	;# c6a c12a 
	movhpd xmm7, [rsi + r9*8 + 8]	;# c6b c12b 
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movapd [rsp + nb333_c6], xmm4
	movapd [rsp + nb333_c12], xmm6

	mov rsi, [rbp + nb333_pos]       ;# base of pos[] 

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
    
    subpd xmm3, [rsp + nb333_ixO]
    subpd xmm4, [rsp + nb333_iyO]
    subpd xmm5, [rsp + nb333_izO]
    
    movapd [rsp + nb333_dxO], xmm3
    movapd [rsp + nb333_dyO], xmm4
    movapd [rsp + nb333_dzO], xmm5
    
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5

	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
    ;# xmm3=rsq

    cvtpd2ps xmm5, xmm3     
    rsqrtps xmm5, xmm5
    cvtps2pd xmm15, xmm5     ;# lu in low xmm2 

    ;# lookup seed in xmm2 
    movapd xmm5, xmm15       ;# copy of lu 
    mulpd xmm15, xmm15        ;# lu*lu 
    movapd xmm7, [rsp + nb333_three]
    mulpd xmm15, xmm3        ;# rsq*lu*lu                    
    movapd xmm6, [rsp + nb333_half]
    subpd xmm7, xmm15        ;# 30-rsq*lu*lu 
    mulpd xmm7, xmm5        
    mulpd xmm7, xmm6        ;# xmm0=iter1 of rinv (new lu) 

    movapd xmm5, xmm7       ;# copy of lu 
    mulpd xmm7, xmm7        ;# lu*lu 
    movapd xmm15, [rsp + nb333_three]
    mulpd xmm7, xmm3        ;# rsq*lu*lu                    
    movapd xmm6, [rsp + nb333_half]
    subpd xmm15, xmm7        ;# 30-rsq*lu*lu 
    mulpd xmm15, xmm5        
    mulpd xmm15, xmm6        ;# xmm15=rinv
        
    mulpd xmm3, xmm15        ;# xmm3=r 

    ;# xmm15=rinv
    ;# xmm3=r

    mulpd xmm3, [rsp + nb333_tsc] ;# rtab

    ;# truncate and convert to integers
    cvttpd2pi mm6, xmm3
    
    ;# convert back to float
    cvtpi2pd  xmm4, mm6
    
    ;# multiply by 4
    pslld   mm6, 2

    ;# calculate eps
    subpd     xmm3, xmm4    ;# xmm3=eps
    
    ;# move to integer registers
    movd r10d, mm6
    psrlq mm6, 32
    movd r11d, mm6
    
    ;# multiply by 3
    lea   r10, [r10 + r10*2]
	lea   r11, [r11 + r11*2]

    ;# xmm3=eps
    ;# xmm15=rinv

	mov rsi, [rbp + nb333_VFtab]
    ;# indices in r10, r11. Load dispersion and repulsion tables in parallel.
    movapd xmm4, [rsi + r10*8 + 32]     ;# Y1d F1d  
    movapd xmm12, [rsi + r11*8 + 32]    ;# Y2d F2d 
    movapd xmm8, [rsi + r10*8 + 64]     ;# Y1r F1r  
    movapd xmm13, [rsi + r11*8 + 64]    ;# Y2r F2r 
    movapd xmm5, xmm4
    movapd xmm9, xmm8
    unpcklpd xmm4, xmm12    ;# Y1d Y2d 
    unpckhpd xmm5, xmm12    ;# F1d F2d 
    unpcklpd xmm8, xmm13    ;# Y1r Y2r 
    unpckhpd xmm9, xmm13    ;# F1r F2r 

    movapd xmm6, [rsi + r10*8 + 48]     ;# G1d H1d  
    movapd xmm12, [rsi + r11*8 + 48]    ;# G2d H2d 
    movapd xmm10, [rsi + r10*8 + 80]    ;# G1r H1r      
    movapd xmm13, [rsi + r11*8 + 80]    ;# G2r H2r 
    movapd xmm7, xmm6
    movapd xmm11, xmm10
    unpcklpd xmm6, xmm12    ;# G1d G2d 
    unpckhpd xmm7, xmm12    ;# H1d H2d 
    unpcklpd xmm10, xmm13   ;# G1r G2r 
    unpckhpd xmm11, xmm13   ;# H1r H2r 
    ;# dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11
    
    mulpd  xmm7, xmm3    ;# Heps
    mulpd  xmm11, xmm3 
    mulpd  xmm6, xmm3   ;# Geps
    mulpd  xmm10, xmm3 
    mulpd  xmm7, xmm3   ;# Heps2
    mulpd  xmm11, xmm3 
    addpd  xmm5, xmm6  ;# F+Geps
    addpd  xmm9, xmm10 
    addpd  xmm5, xmm7   ;# F+Geps+Heps2 = Fp
    addpd  xmm9, xmm11 
    addpd  xmm7, xmm7    ;# 2*Heps2
    addpd  xmm11, xmm11
    addpd  xmm7, xmm6   ;# 2*Heps2+Geps
    addpd  xmm11, xmm10
    
    addpd  xmm7, xmm5  ;# FF = Fp + 2*Heps2 + Geps
    addpd  xmm11, xmm9
    mulpd  xmm5, xmm3  ;# eps*Fp
    mulpd  xmm9, xmm3
    movapd xmm12, [rsp + nb333_c6]
    movapd xmm13, [rsp + nb333_c12]
    addpd  xmm5, xmm4 ;# VV
    addpd  xmm9, xmm8

    mulpd  xmm5, xmm12  ;# VV*c6 = vnb6
    mulpd  xmm9, xmm13  ;# VV*c12 = vnb12
    addpd  xmm5, xmm9
    addpd  xmm5, [rsp + nb333_Vvdwtot]
    movapd [rsp + nb333_Vvdwtot], xmm5
        
    mulpd  xmm7, xmm12   ;# FF*c6 = fnb6
    mulpd  xmm11, xmm13   ;# FF*c12  = fnb12
    addpd  xmm7, xmm11
    
    mulpd  xmm7, [rsp + nb333_tsc]
    mulpd  xmm7, xmm15   ;# -fscal
    xorpd  xmm9, xmm9
    
    subpd  xmm9, xmm7     ;# fscal
    movapd xmm10, xmm9
    movapd xmm11, xmm9

    mulpd  xmm9,  [rsp + nb333_dxO] ;# fx/fy/fz
    mulpd  xmm10, [rsp + nb333_dyO]
    mulpd  xmm11, [rsp + nb333_dzO]

    ;# save j force temporarily
    movapd [rsp + nb333_fjx], xmm9
    movapd [rsp + nb333_fjy], xmm10
    movapd [rsp + nb333_fjz], xmm11
    
    ;# increment i O force
    addpd xmm9, [rsp + nb333_fixO]
    addpd xmm10, [rsp + nb333_fiyO]
    addpd xmm11, [rsp + nb333_fizO]
    movapd [rsp + nb333_fixO], xmm9
    movapd [rsp + nb333_fiyO], xmm10
    movapd [rsp + nb333_fizO], xmm11
    ;# finished O LJ interaction.


    ;# do H1, H2, and M interactions in parallel.
    ;# xmm0-xmm2 still contain j coordinates.                
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subpd xmm0, [rsp + nb333_ixH1]
    subpd xmm1, [rsp + nb333_iyH1]
    subpd xmm2, [rsp + nb333_izH1]
    subpd xmm3, [rsp + nb333_ixH2]
    subpd xmm4, [rsp + nb333_iyH2]
    subpd xmm5, [rsp + nb333_izH2]
    subpd xmm6, [rsp + nb333_ixM]
    subpd xmm7, [rsp + nb333_iyM]
    subpd xmm8, [rsp + nb333_izM]
    
	movapd [rsp + nb333_dxH1], xmm0
	movapd [rsp + nb333_dyH1], xmm1
	movapd [rsp + nb333_dzH1], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb333_dxH2], xmm3
	movapd [rsp + nb333_dyH2], xmm4
	movapd [rsp + nb333_dzH2], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb333_dxM], xmm6
	movapd [rsp + nb333_dyM], xmm7
	movapd [rsp + nb333_dzM], xmm8
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
		
	movapd  xmm9, [rsp + nb333_three]
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

	movapd  xmm15, [rsp + nb333_half]
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
		
	movapd  xmm1, [rsp + nb333_three]
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

	movapd  xmm15, [rsp + nb333_half]
	mulpd   xmm9, xmm15  ;#  rinvH1
	mulpd   xmm10, xmm15 ;#   rinvH2
    mulpd   xmm11, xmm15 ;#   rinvM
	
	movapd  [rsp + nb333_rinvH1], xmm9
	movapd  [rsp + nb333_rinvH2], xmm10
	movapd  [rsp + nb333_rinvM], xmm11
	
	;# interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movapd xmm1, [rsp + nb333_tsc]
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
    
    ;# multiply by three (copy, mult. by two, add back)
    movapd  xmm10, xmm1
    movapd  xmm11, xmm4
    movapd  xmm12, xmm7
    pslld   xmm1, 1
    pslld   xmm4, 1
    pslld   xmm7, 1    
    paddd   xmm1, xmm10
    paddd   xmm4, xmm11
    paddd   xmm7, xmm12    
    
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
        
    mov  rsi, [rbp + nb333_VFtab]

    ;# calculate eps
    subpd     xmm0, xmm2
    subpd     xmm3, xmm5
    subpd     xmm6, xmm8

    movapd    xmm12, xmm0  ;# epsH1
    movapd    xmm13, xmm3  ;# epsH2
    movapd    xmm14, xmm6  ;# epsM

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
    movapd xmm12, [rsp + nb333_qqH]
    movapd xmm13, [rsp + nb333_qqM]
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
    addpd  xmm1, [rsp + nb333_vctot]
    addpd  xmm5, xmm9
    addpd  xmm1, xmm5
    movapd [rsp + nb333_vctot], xmm1
    
    movapd xmm10, [rsp + nb333_tsc]
    mulpd  xmm3, xmm10  ;# fscal
    mulpd  xmm7, xmm10
    mulpd  xmm10, xmm11
    
    xorpd  xmm4, xmm4
    xorpd  xmm8, xmm8
    xorpd  xmm11, xmm11
    
    subpd  xmm4, xmm3
    subpd  xmm8, xmm7
    subpd  xmm11, xmm10

    mulpd  xmm4, [rsp + nb333_rinvH1]
    mulpd  xmm8, [rsp + nb333_rinvH2]
    mulpd  xmm11, [rsp + nb333_rinvM]
    
    ;# move j forces to xmm0-xmm2
    mov rdi, [rbp + nb333_faction]
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
    addpd xmm0, [rsp + nb333_fjx]
    addpd xmm1, [rsp + nb333_fjy]
    addpd xmm2, [rsp + nb333_fjz]

	mulpd xmm3, [rsp + nb333_dxH1]
	mulpd xmm4, [rsp + nb333_dyH1]
	mulpd xmm5, [rsp + nb333_dzH1]
	mulpd xmm7, [rsp + nb333_dxH2]
	mulpd xmm8, [rsp + nb333_dyH2]
	mulpd xmm9, [rsp + nb333_dzH2]
	mulpd xmm10, [rsp + nb333_dxM]
	mulpd xmm11, [rsp + nb333_dyM]
	mulpd xmm12, [rsp + nb333_dzM]

    addpd xmm0, xmm3
    addpd xmm1, xmm4
    addpd xmm2, xmm5
    addpd xmm3, [rsp + nb333_fixH1]
    addpd xmm4, [rsp + nb333_fiyH1]
    addpd xmm5, [rsp + nb333_fizH1]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb333_fixH2]
    addpd xmm8, [rsp + nb333_fiyH2]
    addpd xmm9, [rsp + nb333_fizH2]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb333_fixM]
    addpd xmm11, [rsp + nb333_fiyM]
    addpd xmm12, [rsp + nb333_fizM]

    movapd [rsp + nb333_fixH1], xmm3
    movapd [rsp + nb333_fiyH1], xmm4
    movapd [rsp + nb333_fizH1], xmm5
    movapd [rsp + nb333_fixH2], xmm7
    movapd [rsp + nb333_fiyH2], xmm8
    movapd [rsp + nb333_fizH2], xmm9
    movapd [rsp + nb333_fixM], xmm10
    movapd [rsp + nb333_fiyM], xmm11
    movapd [rsp + nb333_fizM], xmm12
   
    ;# store back j forces from xmm0-xmm2
	movlpd [rdi + rax*8], xmm0
	movlpd [rdi + rax*8 + 8], xmm1
	movlpd [rdi + rax*8 + 16], xmm2
	movhpd [rdi + rbx*8], xmm0
	movhpd [rdi + rbx*8 + 8], xmm1
	movhpd [rdi + rbx*8 + 16], xmm2

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb333_innerk],  2
	jl    .nb333_checksingle
	jmp   .nb333_unroll_loop
.nb333_checksingle:	
	mov   edx, [rsp + nb333_innerk]
	and   edx, 1
	jnz   .nb333_dosingle
	jmp   .nb333_updateouterdata
.nb333_dosingle:
	mov   rdx, [rsp + nb333_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]
	
	add qword ptr [rsp + nb333_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb333_charge]    ;# base of charge[] 
	
	movsd xmm3, [rsi + rax*8]
	movapd xmm4, xmm3	     
	mulsd  xmm3, [rsp + nb333_iqM]
	mulsd  xmm4, [rsp + nb333_iqH]
	movapd  [rsp + nb333_qqM], xmm3
	movapd  [rsp + nb333_qqH], xmm4	

	mov rsi, [rbp + nb333_type]
	mov r8d, [rsi + rax*4]
	mov rsi, [rbp + nb333_vdwparam]
	shl r8d, 1	
	mov edi, [rsp + nb333_ntia]
	add r8d, edi

	movlpd xmm6, [rsi + r8*8]	;# c6a
	movlpd xmm7, [rsi + r8*8 + 8]	;# c12a
	
	movapd [rsp + nb333_c6], xmm6
	movapd [rsp + nb333_c12], xmm7

	mov rsi, [rbp + nb333_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move j coordinates to local temp variables 
    movsd xmm0, [rsi + rax*8] 
    movsd xmm1, [rsi + rax*8 + 8] 
    movsd xmm2, [rsi + rax*8 + 16] 

    ;# xmm0 = jx
    ;# xmm1 = jy
    ;# xmm2 = jz
        
    ;# O interaction
    ;# copy to xmm3-xmm5
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    
    subsd xmm3, [rsp + nb333_ixO]
    subsd xmm4, [rsp + nb333_iyO]
    subsd xmm5, [rsp + nb333_izO]
    
    movapd [rsp + nb333_dxO], xmm3
    movapd [rsp + nb333_dyO], xmm4
    movapd [rsp + nb333_dzO], xmm5
    
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5

	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
    ;# xmm3=rsq

    cvtsd2ss xmm5, xmm3     
    rsqrtss xmm5, xmm5
    cvtss2sd xmm15, xmm5     ;# lu in low xmm2 

    ;# lookup seed in xmm2 
    movapd xmm5, xmm15       ;# copy of lu 
    mulsd xmm15, xmm15        ;# lu*lu 
    movapd xmm7, [rsp + nb333_three]
    mulsd xmm15, xmm3        ;# rsq*lu*lu                    
    movapd xmm6, [rsp + nb333_half]
    subsd xmm7, xmm15        ;# 30-rsq*lu*lu 
    mulsd xmm7, xmm5        
    mulsd xmm7, xmm6        ;# xmm0=iter1 of rinv (new lu) 

    movapd xmm5, xmm7       ;# copy of lu 
    mulsd xmm7, xmm7        ;# lu*lu 
    movapd xmm15, [rsp + nb333_three]
    mulsd xmm7, xmm3        ;# rsq*lu*lu                    
    movapd xmm6, [rsp + nb333_half]
    subsd xmm15, xmm7        ;# 30-rsq*lu*lu 
    mulsd xmm15, xmm5        
    mulsd xmm15, xmm6        ;# xmm15=rinv
        
    mulsd xmm3, xmm15        ;# xmm3=r 

    ;# xmm15=rinv
    ;# xmm3=r

    mulsd xmm3, [rsp + nb333_tsc] ;# rtab

    ;# truncate and convert to integers
    cvttsd2si r10d, xmm3
    
    ;# convert back to float
    cvtsi2sd  xmm4, r10d
    
    ;# multiply by 4
    shl    r10d, 2
    
    ;# calculate eps
    subsd     xmm3, xmm4    ;# xmm3=eps
    
    ;# multiply by 3
    lea   r10, [r10 + r10*2]

    ;# xmm3=eps
    ;# xmm15=rinv

	mov rsi, [rbp + nb333_VFtab]
    movsd  xmm4, [rsi + r10*8 + 32]
    movsd  xmm5, [rsi + r10*8 + 40]
    movsd  xmm6, [rsi + r10*8 + 48]
    movsd  xmm7, [rsi + r10*8 + 56]
    movsd  xmm8, [rsi + r10*8 + 64]
    movsd  xmm9, [rsi + r10*8 + 72]
    movsd  xmm10, [rsi + r10*8 + 80]
    movsd  xmm11, [rsi + r10*8 + 88]
    ;# dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11
    
    mulsd  xmm7, xmm3    ;# Heps
    mulsd  xmm11, xmm3 
    mulsd  xmm6, xmm3   ;# Geps
    mulsd  xmm10, xmm3 
    mulsd  xmm7, xmm3   ;# Heps2
    mulsd  xmm11, xmm3 
    addsd  xmm5, xmm6  ;# F+Geps
    addsd  xmm9, xmm10 
    addsd  xmm5, xmm7   ;# F+Geps+Heps2 = Fp
    addsd  xmm9, xmm11 
    addsd  xmm7, xmm7    ;# 2*Heps2
    addsd  xmm11, xmm11
    addsd  xmm7, xmm6   ;# 2*Heps2+Geps
    addsd  xmm11, xmm10
    
    addsd  xmm7, xmm5  ;# FF = Fp + 2*Heps2 + Geps
    addsd  xmm11, xmm9
    mulsd  xmm5, xmm3  ;# eps*Fp
    mulsd  xmm9, xmm3
    movapd xmm12, [rsp + nb333_c6]
    movapd xmm13, [rsp + nb333_c12]
    addsd  xmm5, xmm4 ;# VV
    addsd  xmm9, xmm8

    mulsd  xmm5, xmm12  ;# VV*c6 = vnb6
    mulsd  xmm9, xmm13  ;# VV*c12 = vnb12
    addsd  xmm5, xmm9
    addsd  xmm5, [rsp + nb333_Vvdwtot]
    movsd [rsp + nb333_Vvdwtot], xmm5
        
    mulsd  xmm7, xmm12   ;# FF*c6 = fnb6
    mulsd  xmm11, xmm13   ;# FF*c12  = fnb12
    addsd  xmm7, xmm11
    
    mulsd  xmm7, [rsp + nb333_tsc]
    mulsd  xmm7, xmm15   ;# -fscal
    xorpd  xmm9, xmm9
    
    subsd  xmm9, xmm7     ;# fscal
    movapd xmm10, xmm9
    movapd xmm11, xmm9

    mulsd  xmm9,  [rsp + nb333_dxO] ;# fx/fy/fz
    mulsd  xmm10, [rsp + nb333_dyO]
    mulsd  xmm11, [rsp + nb333_dzO]

    ;# save j force temporarily
    movapd [rsp + nb333_fjx], xmm9
    movapd [rsp + nb333_fjy], xmm10
    movapd [rsp + nb333_fjz], xmm11
    
    ;# increment i O force
    addsd xmm9, [rsp + nb333_fixO]
    addsd xmm10, [rsp + nb333_fiyO]
    addsd xmm11, [rsp + nb333_fizO]
    movsd [rsp + nb333_fixO], xmm9
    movsd [rsp + nb333_fiyO], xmm10
    movsd [rsp + nb333_fizO], xmm11
    ;# finished O LJ interaction.


    ;# do H1, H2, and M interactions in parallel.
    ;# xmm0-xmm2 still contain j coordinates.                
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subsd xmm0, [rsp + nb333_ixH1]
    subsd xmm1, [rsp + nb333_iyH1]
    subsd xmm2, [rsp + nb333_izH1]
    subsd xmm3, [rsp + nb333_ixH2]
    subsd xmm4, [rsp + nb333_iyH2]
    subsd xmm5, [rsp + nb333_izH2]
    subsd xmm6, [rsp + nb333_ixM]
    subsd xmm7, [rsp + nb333_iyM]
    subsd xmm8, [rsp + nb333_izM]
    
	movapd [rsp + nb333_dxH1], xmm0
	movapd [rsp + nb333_dyH1], xmm1
	movapd [rsp + nb333_dzH1], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [rsp + nb333_dxH2], xmm3
	movapd [rsp + nb333_dyH2], xmm4
	movapd [rsp + nb333_dzH2], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	movapd [rsp + nb333_dxM], xmm6
	movapd [rsp + nb333_dyM], xmm7
	movapd [rsp + nb333_dzM], xmm8
	mulsd  xmm6, xmm6
	mulsd  xmm7, xmm7
	mulsd  xmm8, xmm8
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
    addsd  xmm6, xmm7
    addsd  xmm6, xmm8

	;# start doing invsqrt for j atoms
    cvtsd2ss xmm1, xmm0
    cvtsd2ss xmm4, xmm3
    cvtsd2ss xmm7, xmm6
	rsqrtss xmm1, xmm1
	rsqrtss xmm4, xmm4
    rsqrtss xmm7, xmm7
    cvtss2sd xmm1, xmm1
    cvtss2sd xmm4, xmm4
    cvtss2sd xmm7, xmm7
	
	movapd  xmm2, xmm1
	movapd  xmm5, xmm4
    movapd  xmm8, xmm7
    
	mulsd   xmm1, xmm1 ;# lu*lu
	mulsd   xmm4, xmm4 ;# lu*lu
    mulsd   xmm7, xmm7 ;# lu*lu
		
	movapd  xmm9, [rsp + nb333_three]
	movapd  xmm10, xmm9
    movapd  xmm11, xmm9

	mulsd   xmm1, xmm0 ;# rsq*lu*lu
	mulsd   xmm4, xmm3 ;# rsq*lu*lu 
    mulsd   xmm7, xmm6 ;# rsq*lu*lu
	
	subsd   xmm9, xmm1
	subsd   xmm10, xmm4
    subsd   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm2
	mulsd   xmm10, xmm5
    mulsd   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movapd  xmm15, [rsp + nb333_half]
	mulsd   xmm9, xmm15  ;# first iteration for rinvH1
	mulsd   xmm10, xmm15 ;# first iteration for rinvH2
    mulsd   xmm11, xmm15 ;# first iteration for rinvM

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulsd   xmm2, xmm2 ;# lu*lu
	mulsd   xmm5, xmm5 ;# lu*lu
    mulsd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb333_three]
	movapd  xmm4, xmm1
    movapd  xmm7, xmm1

	mulsd   xmm2, xmm0 ;# rsq*lu*lu
	mulsd   xmm5, xmm3 ;# rsq*lu*lu 
    mulsd   xmm8, xmm6 ;# rsq*lu*lu
	
	subsd   xmm1, xmm2
	subsd   xmm4, xmm5
    subsd   xmm7, xmm8 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm1
	mulsd   xmm10, xmm4
    mulsd   xmm11, xmm7 ;# lu*(3-rsq*lu*lu)

	movapd  xmm15, [rsp + nb333_half]
	mulsd   xmm9, xmm15  ;#  rinvH1
	mulsd   xmm10, xmm15 ;#   rinvH2
    mulsd   xmm11, xmm15 ;#   rinvM
	
	movapd  [rsp + nb333_rinvH1], xmm9
	movapd  [rsp + nb333_rinvH2], xmm10
	movapd  [rsp + nb333_rinvM], xmm11
	
	;# interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movapd xmm1, [rsp + nb333_tsc]
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
    shl   r8d, 2
    shl   r10d, 2
    shl   r12d, 2
    
    mov  rsi, [rbp + nb333_VFtab]

	lea   r8, [r8 + r8*2]
    lea   r10, [r10 + r10*2]
    lea   r12, [r12 + r12*2]
            
    ;# calculate eps
    subsd     xmm0, xmm2
    subsd     xmm3, xmm5
    subsd     xmm6, xmm8

    movapd    xmm12, xmm0  ;# epsH1
    movapd    xmm13, xmm3  ;# epsH2
    movapd    xmm14, xmm6  ;# epsM

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
    movapd xmm12, [rsp + nb333_qqH]
    movapd xmm13, [rsp + nb333_qqM]
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
    addsd  xmm1, [rsp + nb333_vctot]
    addsd  xmm5, xmm9
    addsd  xmm1, xmm5
    movsd [rsp + nb333_vctot], xmm1

    movapd xmm10, [rsp + nb333_tsc]
    mulsd  xmm3, xmm10  ;# fscal
    mulsd  xmm7, xmm10
    mulsd  xmm10, xmm11
    
    xorpd  xmm4, xmm4
    xorpd  xmm8, xmm8
    xorpd  xmm11, xmm11
    
    subsd  xmm4, xmm3
    subsd  xmm8, xmm7
    subsd  xmm11, xmm10

    mulsd  xmm4, [rsp + nb333_rinvH1]
    mulsd  xmm8, [rsp + nb333_rinvH2]
    mulsd  xmm11, [rsp + nb333_rinvM]
    
    ;# move j forces to xmm0-xmm2
    mov rdi, [rbp + nb333_faction]
	movsd xmm0, [rdi + rax*8]
	movsd xmm1, [rdi + rax*8 + 8]
	movsd xmm2, [rdi + rax*8 + 16]

    movapd xmm3, xmm4
    movapd xmm5, xmm4
    movapd xmm7, xmm8
    movapd xmm9, xmm8
    movapd xmm10, xmm11
    movapd xmm12, xmm11

    ;# add forces from O interaction
    addsd xmm0, [rsp + nb333_fjx]
    addsd xmm1, [rsp + nb333_fjy]
    addsd xmm2, [rsp + nb333_fjz]

	mulsd xmm3, [rsp + nb333_dxH1]
	mulsd xmm4, [rsp + nb333_dyH1]
	mulsd xmm5, [rsp + nb333_dzH1]
	mulsd xmm7, [rsp + nb333_dxH2]
	mulsd xmm8, [rsp + nb333_dyH2]
	mulsd xmm9, [rsp + nb333_dzH2]
	mulsd xmm10, [rsp + nb333_dxM]
	mulsd xmm11, [rsp + nb333_dyM]
	mulsd xmm12, [rsp + nb333_dzM]

    addsd xmm0, xmm3
    addsd xmm1, xmm4
    addsd xmm2, xmm5
    addsd xmm3, [rsp + nb333_fixH1]
    addsd xmm4, [rsp + nb333_fiyH1]
    addsd xmm5, [rsp + nb333_fizH1]

    addsd xmm0, xmm7
    addsd xmm1, xmm8
    addsd xmm2, xmm9
    addsd xmm7, [rsp + nb333_fixH2]
    addsd xmm8, [rsp + nb333_fiyH2]
    addsd xmm9, [rsp + nb333_fizH2]

    addsd xmm0, xmm10
    addsd xmm1, xmm11
    addsd xmm2, xmm12
    addsd xmm10, [rsp + nb333_fixM]
    addsd xmm11, [rsp + nb333_fiyM]
    addsd xmm12, [rsp + nb333_fizM]

    movsd [rsp + nb333_fixH1], xmm3
    movsd [rsp + nb333_fiyH1], xmm4
    movsd [rsp + nb333_fizH1], xmm5
    movsd [rsp + nb333_fixH2], xmm7
    movsd [rsp + nb333_fiyH2], xmm8
    movsd [rsp + nb333_fizH2], xmm9
    movsd [rsp + nb333_fixM], xmm10
    movsd [rsp + nb333_fiyM], xmm11
    movsd [rsp + nb333_fizM], xmm12
   
    ;# store back j forces from xmm0-xmm2
	movsd [rdi + rax*8], xmm0
	movsd [rdi + rax*8 + 8], xmm1
	movsd [rdi + rax*8 + 16], xmm2

.nb333_updateouterdata:
	mov   ecx, [rsp + nb333_ii3]
	mov   rdi, [rbp + nb333_faction]
	mov   rsi, [rbp + nb333_fshift]
	mov   edx, [rsp + nb333_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb333_fixO]
	movapd xmm1, [rsp + nb333_fiyO]
	movapd xmm2, [rsp + nb333_fizO]

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
	movapd xmm0, [rsp + nb333_fixH1]
	movapd xmm1, [rsp + nb333_fiyH1]
	movapd xmm2, [rsp + nb333_fizH1]

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
	movapd xmm0, [rsp + nb333_fixH2]
	movapd xmm1, [rsp + nb333_fiyH2]
	movapd xmm2, [rsp + nb333_fizH2]

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
	movapd xmm0, [rsp + nb333_fixM]
	movapd xmm1, [rsp + nb333_fiyM]
	movapd xmm2, [rsp + nb333_fizM]

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
	mov esi, [rsp + nb333_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb333_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb333_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb333_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb333_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb333_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
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
	





.globl nb_kernel333nf_x86_64_sse2
.globl _nb_kernel333nf_x86_64_sse2
nb_kernel333nf_x86_64_sse2:	
_nb_kernel333nf_x86_64_sse2:	
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
	;# bottom of stack is cache-aligned for sse2 use 
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
.equiv          nb333nf_rO,             256
.equiv          nb333nf_rH1,            272
.equiv          nb333nf_rH2,            288
.equiv          nb333nf_rM,             304
.equiv          nb333nf_tsc,            320
.equiv          nb333nf_c6,             336
.equiv          nb333nf_c12,            352
.equiv          nb333nf_vctot,          368
.equiv          nb333nf_Vvdwtot,        384
.equiv          nb333nf_half,           400
.equiv          nb333nf_three,          416
.equiv          nb333nf_is3,            432
.equiv          nb333nf_ii3,            436
.equiv          nb333nf_nri,            440
.equiv          nb333nf_iinr,           448
.equiv          nb333nf_jindex,         456
.equiv          nb333nf_jjnr,           464
.equiv          nb333nf_shift,          472
.equiv          nb333nf_shiftvec,       480
.equiv          nb333nf_facel,          488
.equiv          nb333nf_innerjjnr,      496
.equiv          nb333nf_ntia,           504
.equiv          nb333nf_innerk,         508
.equiv          nb333nf_n,              512
.equiv          nb333nf_nn1,            516
.equiv          nb333nf_nouter,         520
.equiv          nb333nf_ninner,         524

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
	sub rsp, 528		;# local variable stack space (n*16+8)

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
	movsd xmm0, [rsi]
	movsd [rsp + nb333nf_facel], xmm0

	mov rax, [rbp + nb333nf_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb333nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb333nf_half], eax
	mov [rsp + nb333nf_half+4], ebx
	movsd xmm1, [rsp + nb333nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb333nf_half], xmm1
	movapd [rsp + nb333nf_three], xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb333nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb333nf_charge]
	movsd xmm3, [rdx + rbx*8 + 8]	
	movsd xmm4, [rdx + rbx*8 + 24]	
	mov rsi, [rbp + nb333nf_p_facel]
	movsd xmm0, [rsi]
	movsd xmm5, [rsp + nb333nf_facel]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb333nf_iqH], xmm3
	movapd [rsp + nb333nf_iqM], xmm4
	
	mov   rdx, [rbp + nb333nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb333nf_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
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
	mov   rax, [rsp + nb333nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb333nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb333nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb333nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb333nf_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb333nf_ii3], ebx

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
	movapd [rsp + nb333nf_ixO], xmm3
	movapd [rsp + nb333nf_iyO], xmm4
	movapd [rsp + nb333nf_izO], xmm5
	movapd [rsp + nb333nf_ixH1], xmm6
	movapd [rsp + nb333nf_iyH1], xmm7

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
	movapd [rsp + nb333nf_izH1], xmm6
	movapd [rsp + nb333nf_ixH2], xmm0
	movapd [rsp + nb333nf_iyH2], xmm1
	movapd [rsp + nb333nf_izH2], xmm2
	movapd [rsp + nb333nf_ixM], xmm3
	movapd [rsp + nb333nf_iyM], xmm4
	movapd [rsp + nb333nf_izM], xmm5
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [rsp + nb333nf_vctot], xmm4
	movapd [rsp + nb333nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb333nf_jindex]
	mov   ecx, [rax+rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb333nf_pos]
	mov   rdi, [rbp + nb333nf_faction]	
	mov   rax, [rsp + nb333nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb333nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb333nf_ninner]
	mov   [rsp + nb333nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb333nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb333nf_unroll_loop
	jmp   .nb333nf_checksingle
.nb333nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb333nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]
	
	add qword ptr [rsp + nb333nf_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb333nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]
	movapd xmm4, xmm3	     
	mulpd  xmm3, [rsp + nb333nf_iqM]
	mulpd  xmm4, [rsp + nb333nf_iqH]
	movapd  [rsp + nb333nf_qqM], xmm3
	movapd  [rsp + nb333nf_qqH], xmm4	

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	mov rsi, [rbp + nb333nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov rsi, [rbp + nb333nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	mov edi, [rsp + nb333nf_ntia]
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
	movapd [rsp + nb333nf_c6], xmm4
	movapd [rsp + nb333nf_c12], xmm6

	mov rsi, [rbp + nb333nf_pos]       ;# base of pos[] 

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
	movapd xmm4, [rsp + nb333nf_ixO]
	movapd xmm5, [rsp + nb333nf_iyO]
	movapd xmm6, [rsp + nb333nf_izO]

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
	movapd xmm4, [rsp + nb333nf_ixH1]
	movapd xmm5, [rsp + nb333nf_iyH1]
	movapd xmm6, [rsp + nb333nf_izH1]

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
	movapd xmm3, [rsp + nb333nf_ixH2]
	movapd xmm4, [rsp + nb333nf_iyH2]
	movapd xmm5, [rsp + nb333nf_izH2]

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
	movapd xmm3, [rsp + nb333nf_iyM]
	movapd xmm4, [rsp + nb333nf_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [rsp + nb333nf_ixM]
	subpd  xmm2, xmm0	

	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# rsqO - put seed in xmm2 
	cvtpd2ps xmm2, xmm7	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm0, [rsp + nb333nf_three]
	mulpd   xmm2, xmm7	;# rsq*lu*lu 
	subpd   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm0, [rsp + nb333nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm7
	movapd xmm3, xmm0
	mulpd xmm0, xmm0	;# lu*lu 
	mulpd xmm2, xmm0	;# rsq*lu*lu 
	movapd xmm0, [rsp + nb333nf_three]
	subpd xmm0, xmm2	;# 3-rsq*lu*lu 
	mulpd xmm0, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm0, [rsp + nb333nf_half] ;# rinv 
	mulpd   xmm7, xmm0
	movapd  [rsp + nb333nf_rO], xmm7	;# r in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm0, [rsp + nb333nf_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm0, [rsp + nb333nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm6
	movapd xmm3, xmm0
	mulpd xmm0, xmm0	;# lu*lu 
	mulpd xmm2, xmm0	;# rsq*lu*lu 
	movapd xmm0, [rsp + nb333nf_three]
	subpd xmm0, xmm2	;# 3-rsq*lu*lu 
	mulpd xmm0, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm0, [rsp + nb333nf_half] ;# rinv 
	mulpd  xmm6, xmm0
	movapd [rsp + nb333nf_rH1], xmm6	;# rH1 
	
	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm0, [rsp + nb333nf_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm0, [rsp + nb333nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm5
	movapd xmm3, xmm0
	mulpd xmm0, xmm0	;# lu*lu 
	mulpd xmm2, xmm0	;# rsq*lu*lu 
	movapd xmm0, [rsp + nb333nf_three]
	subpd xmm0, xmm2	;# 3-rsq*lu*lu 
	mulpd xmm0, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm0, [rsp + nb333nf_half] ;# rinv 
	mulpd xmm5, xmm0
	movapd [rsp + nb333nf_rH2], xmm5 ;# r 

	;# rsqM - seed in xmm2 
	cvtpd2ps xmm2, xmm4	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm0, [rsp + nb333nf_three]
	mulpd   xmm2, xmm4	;# rsq*lu*lu 
	subpd   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm0, [rsp + nb333nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm4
	movapd xmm3, xmm0
	mulpd xmm0, xmm0	;# lu*lu 
	mulpd xmm2, xmm0	;# rsq*lu*lu 
	movapd xmm0, [rsp + nb333nf_three]
	subpd xmm0, xmm2	;# 3-rsq*lu*lu 
	mulpd xmm0, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm0, [rsp + nb333nf_half] ;# rinv 
	mulpd xmm4, xmm0
	movapd [rsp + nb333nf_rM], xmm4 ;# r 

	;# do O interactions 
	;# rO is still in xmm7 
	mulpd   xmm7, [rsp + nb333nf_tsc]
	cvttpd2pi mm6, xmm7	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	movd mm0, eax	
	movd mm1, ebx
	mov  rsi, [rbp + nb333nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   rax, [rax + rax*2] ;# idx *= 3 (total *=12 now) 
	lea   rbx, [rbx + rbx*2]
	
	;# Dispersion 
	movapd xmm4, [rsi + rax*8 + 32]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8 + 32]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 48]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 48]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# Dispersion table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [rsp + nb333nf_c6]
	mulpd  xmm5, xmm4	 ;# Vvdw6 		
	addpd  xmm5, [rsp + nb333nf_Vvdwtot]
	movapd [rsp + nb333nf_Vvdwtot], xmm5
	
	;# Repulsion 
	movapd xmm4, [rsi + rax*8 + 64]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8 + 64]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 80]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 80]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# Dispersion table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [rsp + nb333nf_c12]
	mulpd  xmm5, xmm4 ;# Vvdw12 
	addpd  xmm5, [rsp + nb333nf_Vvdwtot]
	movapd [rsp + nb333nf_Vvdwtot], xmm5

	;# Done with O interactions - now H1! 
	movapd xmm7, [rsp + nb333nf_rH1]
	mulpd xmm7, [rsp + nb333nf_tsc]
	cvttpd2pi mm6, xmm7	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb333nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   rax, [rax + rax*2] ;# idx *= 3 (total *=12 now) 	
	lea   rbx, [rbx + rbx*2]
	
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
	movapd xmm3, [rsp + nb333nf_qqH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    	;# increment vcoul 
	xorpd  xmm4, xmm4
    	addpd  xmm5, [rsp + nb333nf_vctot]
    	movapd [rsp + nb333nf_vctot], xmm5 

	;# H2 interactions 
	movapd xmm7, [rsp + nb333nf_rH2]
	mulpd   xmm7, [rsp + nb333nf_tsc]
	cvttpd2pi mm6, xmm7	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb333nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   rax, [rax + rax*2] ;# idx *= 3 (total *=12 now)
	lea   rbx, [rbx + rbx*2]

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
	movapd xmm3, [rsp + nb333nf_qqH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    	;# increment vcoul 
	xorpd  xmm4, xmm4
    	addpd  xmm5, [rsp + nb333nf_vctot]
    	movapd [rsp + nb333nf_vctot], xmm5 

	;# M interactions 
	movapd xmm7, [rsp + nb333nf_rM]
	mulpd   xmm7, [rsp + nb333nf_tsc]
	cvttpd2pi mm6, xmm7	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb333nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   rax, [rax + rax*2] ;# idx *= 3 (total *=12 now)
	lea   rbx, [rbx + rbx*2]

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
	movapd xmm3, [rsp + nb333nf_qqM]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    	;# increment vcoul 
	xorpd  xmm4, xmm4
    	addpd  xmm5, [rsp + nb333nf_vctot]
    	movapd [rsp + nb333nf_vctot], xmm5 

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb333nf_innerk],  2
	jl    .nb333nf_checksingle
	jmp   .nb333nf_unroll_loop
.nb333nf_checksingle:	
	mov   edx, [rsp + nb333nf_innerk]
	and   edx, 1
	jnz   .nb333nf_dosingle
	jmp   .nb333nf_updateouterdata
.nb333nf_dosingle:
	mov   rdx, [rsp + nb333nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]

	mov rsi, [rbp + nb333nf_charge]    ;# base of charge[] 
	xorpd xmm3, xmm3	
	movlpd xmm3, [rsi + rax*8]
	movapd xmm4, xmm3	     
	mulpd  xmm3, [rsp + nb333nf_iqM]
	mulpd  xmm4, [rsp + nb333nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movapd  [rsp + nb333nf_qqM], xmm3
	movapd  [rsp + nb333nf_qqH], xmm4	
	
	mov rsi, [rbp + nb333nf_type]
	mov eax, [rsi + rax*4]
	mov rsi, [rbp + nb333nf_vdwparam]
	shl eax, 1	
	mov edi, [rsp + nb333nf_ntia]
	add eax, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movapd [rsp + nb333nf_c6], xmm4
	movapd [rsp + nb333nf_c12], xmm6
	
	mov rsi, [rbp + nb333nf_pos]       ;# base of pos[] 
	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	
	;# move coords to xmm0-xmm2 
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [rsp + nb333nf_ixO]
	movapd xmm5, [rsp + nb333nf_iyO]
	movapd xmm6, [rsp + nb333nf_izO]

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
	movapd xmm4, [rsp + nb333nf_ixH1]
	movapd xmm5, [rsp + nb333nf_iyH1]
	movapd xmm6, [rsp + nb333nf_izH1]

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
	movapd xmm3, [rsp + nb333nf_ixH2]
	movapd xmm4, [rsp + nb333nf_iyH2]
	movapd xmm5, [rsp + nb333nf_izH2]

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
	movapd xmm3, [rsp + nb333nf_iyM]
	movapd xmm4, [rsp + nb333nf_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [rsp + nb333nf_ixM]
	subpd  xmm2, xmm0	

	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# start with rsqO - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm0, [rsp + nb333nf_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm0, [rsp + nb333nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm7
	movapd xmm3, xmm0
	mulsd xmm0, xmm0	;# lu*lu 
	mulsd xmm2, xmm0	;# rsq*lu*lu 
	movapd xmm0, [rsp + nb333nf_three]
	subsd xmm0, xmm2	;# 3-rsq*lu*lu 
	mulsd xmm0, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm0, [rsp + nb333nf_half] ;# rinv 
	mulsd   xmm7, xmm0
	movapd  [rsp + nb333nf_rO], xmm7	;# r in xmm7 

	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm0, [rsp + nb333nf_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm0, [rsp + nb333nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm6
	movapd xmm3, xmm0
	mulsd xmm0, xmm0	;# lu*lu 
	mulsd xmm2, xmm0	;# rsq*lu*lu 
	movapd xmm0, [rsp + nb333nf_three]
	subsd xmm0, xmm2	;# 3-rsq*lu*lu 
	mulsd xmm0, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm0, [rsp + nb333nf_half] ;# rinv 
	mulsd  xmm6, xmm0
	movapd [rsp + nb333nf_rH1], xmm6	;# rH1 
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm0, [rsp + nb333nf_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm0, [rsp + nb333nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm5
	movapd xmm3, xmm0
	mulsd xmm0, xmm0	;# lu*lu 
	mulsd xmm2, xmm0	;# rsq*lu*lu 
	movapd xmm0, [rsp + nb333nf_three]
	subsd xmm0, xmm2	;# 3-rsq*lu*lu 
	mulsd xmm0, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm0, [rsp + nb333nf_half] ;# rinv 
	mulsd xmm5, xmm0
	movapd [rsp + nb333nf_rH2], xmm5 ;# r 

	;# rsqM - seed in xmm2 
	cvtsd2ss xmm2, xmm4	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm0, [rsp + nb333nf_three]
	mulsd   xmm2, xmm4	;# rsq*lu*lu 
	subsd   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm0, [rsp + nb333nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm4
	movapd xmm3, xmm0
	mulsd xmm0, xmm0	;# lu*lu 
	mulsd xmm2, xmm0	;# rsq*lu*lu 
	movapd xmm0, [rsp + nb333nf_three]
	subsd xmm0, xmm2	;# 3-rsq*lu*lu 
	mulsd xmm0, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm0, [rsp + nb333nf_half] ;# rinv 
	mulsd xmm4, xmm0
	movapd [rsp + nb333nf_rM], xmm4 ;# r 

	;# do O interactions 
	movd mm0, eax	
	;# rO is still in xmm7 
	mulsd   xmm7, [rsp + nb333nf_tsc]
	cvttsd2si eax, xmm7	;# lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm7
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	 	
	shl eax, 2		
	mov  rsi, [rbp + nb333nf_VFtab]
	lea   rax, [rax + rax*2] ;# idx *= 3 (total *=12 now) 	
	
	;# Dispersion 
	movapd xmm4, [rsi + rax*8 + 32]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + rax*8 + 48]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# Dispersion table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [rsp + nb333nf_c6]
	mulsd  xmm5, xmm4	 ;# Vvdw6 
	
	addsd  xmm5, [rsp + nb333nf_Vvdwtot]
	movsd [rsp + nb333nf_Vvdwtot], xmm5
	
	;# Repulsion 
	movapd xmm4, [rsi + rax*8 + 64]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + rax*8 + 80]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# Dispersion table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [rsp + nb333nf_c12]
	mulsd  xmm5, xmm4 ;# Vvdw12 
	addsd  xmm5, [rsp + nb333nf_Vvdwtot]
	movsd [rsp + nb333nf_Vvdwtot], xmm5

	;# Done with O interactions - now H1! 
	movapd xmm7, [rsp + nb333nf_rH1]
	mulpd xmm7, [rsp + nb333nf_tsc]
	cvttsd2si eax, xmm7	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subpd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb333nf_VFtab]
	lea   rax, [rax + rax*2] ;# idx *= 3 (total *=12 now)	

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
	movapd xmm3, [rsp + nb333nf_qqH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    	;# increment vcoul 
	xorpd  xmm4, xmm4
    	addsd  xmm5, [rsp + nb333nf_vctot]
    	movlpd [rsp + nb333nf_vctot], xmm5 

	;#  H2 interactions 
	movapd xmm7, [rsp + nb333nf_rH2]
	mulsd   xmm7, [rsp + nb333nf_tsc]
	cvttsd2si eax, xmm7	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb333nf_VFtab]
	lea   rax, [rax + rax*2] ;# idx *= 3 (total *=12 now)	

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
	movapd xmm3, [rsp + nb333nf_qqH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    	;# increment vcoul 
	xorpd  xmm4, xmm4
    	addsd  xmm5, [rsp + nb333nf_vctot]
    	movlpd [rsp + nb333nf_vctot], xmm5 

	;# M interactions 
	movapd xmm7, [rsp + nb333nf_rM]
	mulsd   xmm7, [rsp + nb333nf_tsc]
	cvttsd2si eax, xmm7	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb333nf_VFtab]
	lea   rax, [rax + rax*2] ;# idx *= 3 (total *=12 now)	

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
	movapd xmm3, [rsp + nb333nf_qqM]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    	;# increment vcoul 
	xorpd  xmm4, xmm4
    	addsd  xmm5, [rsp + nb333nf_vctot]
    	movlpd [rsp + nb333nf_vctot], xmm5 

.nb333nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb333nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb333nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb333nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb333nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb333nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb333nf_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
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

	add rsp, 528
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
