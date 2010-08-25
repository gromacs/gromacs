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

.globl nb_kernel133_x86_64_sse2
.globl _nb_kernel133_x86_64_sse2
nb_kernel133_x86_64_sse2:	
_nb_kernel133_x86_64_sse2:	
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
	;# bottom of stack is cache-aligned for sse2 use 
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
.equiv          nb133_iqH,              192
.equiv          nb133_iqM,              208
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
.equiv          nb133_qqH,              416
.equiv          nb133_qqM,              432
.equiv          nb133_c6,               448
.equiv          nb133_c12,              464
.equiv          nb133_tsc,              480
.equiv          nb133_fstmp,            496
.equiv          nb133_vctot,            512
.equiv          nb133_Vvdwtot,          528
.equiv          nb133_fixO,             544
.equiv          nb133_fiyO,             560
.equiv          nb133_fizO,             576
.equiv          nb133_fixH1,            592
.equiv          nb133_fiyH1,            608
.equiv          nb133_fizH1,            624
.equiv          nb133_fixH2,            640
.equiv          nb133_fiyH2,            656
.equiv          nb133_fizH2,            672
.equiv          nb133_fixM,             688
.equiv          nb133_fiyM,             704
.equiv          nb133_fizM,             720
.equiv          nb133_fjx,              736
.equiv          nb133_fjy,              752
.equiv          nb133_fjz,              768
.equiv          nb133_half,             784
.equiv          nb133_three,            800
.equiv          nb133_two,              816
.equiv          nb133_rinvH1,           832
.equiv          nb133_rinvH2,           848
.equiv          nb133_rinvM,            864
.equiv          nb133_krsqH1,           880
.equiv          nb133_krsqH2,           896
.equiv          nb133_krsqM,            912
.equiv          nb133_krf,              928
.equiv          nb133_crf,              944
.equiv          nb133_rsqO,             960
.equiv          nb133_facel,            976	
.equiv          nb133_iinr,             992
.equiv          nb133_jindex,           1000
.equiv          nb133_jjnr,             1008
.equiv          nb133_shift,            1016
.equiv          nb133_shiftvec,         1024
.equiv          nb133_innerjjnr,        1032
.equiv          nb133_is3,              1040
.equiv          nb133_ii3,              1044
.equiv          nb133_nri,              1048
.equiv          nb133_ntia,             1052
.equiv          nb133_innerk,           1056
.equiv          nb133_n,                1060
.equiv          nb133_nn1,              1064
.equiv          nb133_nouter,           1068
.equiv          nb133_ninner,           1072

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
	movsd xmm0, [rsi]
	movsd [rsp + nb133_facel], xmm0

	mov rax, [rbp + nb133_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb133_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb133_half], eax
	mov [rsp + nb133_half+4], ebx
	movsd xmm1, [rsp + nb133_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb133_half], xmm1
	movapd [rsp + nb133_two], xmm2
	movapd [rsp + nb133_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb133_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb133_charge]
	movsd xmm3, [rdx + rbx*8 + 8]	
	movsd xmm4, [rdx + rbx*8 + 24]	

	movsd xmm5, [rsp + nb133_facel]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb133_iqH], xmm3
	movapd [rsp + nb133_iqM], xmm4
	
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
	mov   ebx, [rax+rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb133_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb133_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb133_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb133_pos]    ;# eax = base of pos[]  
	mov   [rsp + nb133_ii3], ebx

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
	movapd [rsp + nb133_ixO], xmm3
	movapd [rsp + nb133_iyO], xmm4
	movapd [rsp + nb133_izO], xmm5
	movapd [rsp + nb133_ixH1], xmm6
	movapd [rsp + nb133_iyH1], xmm7

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
	movapd [rsp + nb133_izH1], xmm6
	movapd [rsp + nb133_ixH2], xmm0
	movapd [rsp + nb133_iyH2], xmm1
	movapd [rsp + nb133_izH2], xmm2
	movapd [rsp + nb133_ixM], xmm3
	movapd [rsp + nb133_iyM], xmm4
	movapd [rsp + nb133_izM], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb133_vctot], xmm4
	movapd [rsp + nb133_Vvdwtot], xmm4
	movapd [rsp + nb133_fixO], xmm4
	movapd [rsp + nb133_fiyO], xmm4
	movapd [rsp + nb133_fizO], xmm4
	movapd [rsp + nb133_fixH1], xmm4
	movapd [rsp + nb133_fiyH1], xmm4
	movapd [rsp + nb133_fizH1], xmm4
	movapd [rsp + nb133_fixH2], xmm4
	movapd [rsp + nb133_fiyH2], xmm4
	movapd [rsp + nb133_fizH2], xmm4
	movapd [rsp + nb133_fixM], xmm4
	movapd [rsp + nb133_fiyM], xmm4
	movapd [rsp + nb133_fizM], xmm4
	
	mov   rax, [rsp + nb133_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb133_pos]
	mov   rdi, [rbp + nb133_faction]	
	mov   rax, [rsp + nb133_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb133_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb133_ninner]
	mov   [rsp + nb133_ninner], ecx
	add   edx, 0
	mov   [rsp + nb133_innerk], edx    ;# number of innerloop atoms 
	jge   .nb133_unroll_loop
	jmp   .nb133_checksingle
.nb133_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb133_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]

	add qword ptr [rsp + nb133_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb133_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [rsp + nb133_iqM]
	mulpd  xmm4, [rsp + nb133_iqH]

	movapd  [rsp + nb133_qqM], xmm3
	movapd  [rsp + nb133_qqH], xmm4
	
	mov rsi, [rbp + nb133_type]
	mov r8d, [rsi + rax*4]
	mov r9d, [rsi + rbx*4]
	mov rsi, [rbp + nb133_vdwparam]
	shl r8d, 1	
	shl r9d, 1	
	mov edi, [rsp + nb133_ntia]
	add r8d, edi
	add r9d, edi

	movlpd xmm6, [rsi + r8*8]	;# c6a
	movlpd xmm7, [rsi + r9*8]	;# c6b
	movhpd xmm6, [rsi + r8*8 + 8]	;# c6a c12a 
	movhpd xmm7, [rsi + r9*8 + 8]	;# c6b c12b 
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movapd [rsp + nb133_c6], xmm4
	movapd [rsp + nb133_c12], xmm6
	
	mov rsi, [rbp + nb133_pos]       ;# base of pos[] 

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
    
    subpd xmm3, [rsp + nb133_ixO]
    subpd xmm4, [rsp + nb133_iyO]
    subpd xmm5, [rsp + nb133_izO]
    
    movapd [rsp + nb133_dxO], xmm3
    movapd [rsp + nb133_dyO], xmm4
    movapd [rsp + nb133_dzO], xmm5
    
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5

	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
    ;# xmm3=rsq

    cvtpd2ps xmm15, xmm3     
    rsqrtps xmm15, xmm15
    cvtps2pd xmm15, xmm15     ;# lu in low xmm2 

    ;# lookup seed in xmm15
    movapd xmm5, xmm15       ;# copy of lu 
    mulpd xmm15, xmm15        ;# lu*lu 
    movapd xmm7, [rsp + nb133_three]
    mulpd xmm15, xmm3        ;# rsq*lu*lu                    
    movapd xmm6, [rsp + nb133_half]
    subpd xmm7, xmm15        ;# 30-rsq*lu*lu 
    mulpd xmm7, xmm5        
    mulpd xmm7, xmm6        ;# xmm0=iter1 of rinv (new lu) 

    movapd xmm5, xmm7       ;# copy of lu 
    mulpd xmm7, xmm7        ;# lu*lu 
    movapd xmm15, [rsp + nb133_three]
    mulpd xmm7, xmm3        ;# rsq*lu*lu                    
    movapd xmm6, [rsp + nb133_half]
    subpd xmm15, xmm7        ;# 30-rsq*lu*lu 
    mulpd xmm15, xmm5        
    mulpd xmm15, xmm6        ;# xmm15=rinv
        
    mulpd xmm3, xmm15        ;# xmm3=r 

    ;# xmm15=rinv
    ;# xmm3=r

    mulpd xmm3, [rsp + nb133_tsc] ;# rtab

    ;# truncate and convert to integers
    cvttpd2pi mm6, xmm3
    
    ;# convert back to float
    cvtpi2pd  xmm4, mm6
    
    ;# multiply by 8
    pslld   mm6, 3

    ;# calculate eps
    subpd     xmm3, xmm4    ;# xmm3=eps
    
    ;# move to integer registers
    movd r10d, mm6
    psrlq mm6, 32
    movd r11d, mm6
    
    ;# xmm3=eps
    ;# xmm15=rinv

	mov rsi, [rbp + nb133_VFtab]
    ;# indices in r10, r11. Load dispersion and repulsion tables in parallel.
    movapd xmm4, [rsi + r10*8]          ;# Y1d F1d  
    movapd xmm12, [rsi + r11*8]         ;# Y2d F2d 
    movapd xmm8, [rsi + r10*8 + 32]     ;# Y1r F1r  
    movapd xmm13, [rsi + r11*8 + 32]    ;# Y2r F2r 
    movapd xmm5, xmm4
    movapd xmm9, xmm8
    unpcklpd xmm4, xmm12    ;# Y1d Y2d 
    unpckhpd xmm5, xmm12    ;# F1d F2d 
    unpcklpd xmm8, xmm13    ;# Y1r Y2r 
    unpckhpd xmm9, xmm13    ;# F1r F2r 

    movapd xmm6, [rsi + r10*8 + 16]     ;# G1d H1d  
    movapd xmm12, [rsi + r11*8 + 16]        ;# G2d H2d 
    movapd xmm10, [rsi + r10*8 + 48]        ;# G1r H1r      
    movapd xmm13, [rsi + r11*8 + 48]        ;# G2r H2r 
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
    movapd xmm12, [rsp + nb133_c6]
    movapd xmm13, [rsp + nb133_c12]
    addpd  xmm5, xmm4 ;# VV
    addpd  xmm9, xmm8

    mulpd  xmm5, xmm12  ;# VV*c6 = vnb6
    mulpd  xmm9, xmm13  ;# VV*c12 = vnb12
    addpd  xmm5, xmm9
    addpd  xmm5, [rsp + nb133_Vvdwtot]
    movapd [rsp + nb133_Vvdwtot], xmm5
        
    mulpd  xmm7, xmm12   ;# FF*c6 = fnb6
    mulpd  xmm11, xmm13   ;# FF*c12  = fnb12
    addpd  xmm7, xmm11
    
    mulpd  xmm7, [rsp + nb133_tsc]
    mulpd  xmm7, xmm15   ;# -fscal
    xorpd  xmm9, xmm9
    
    subpd  xmm9, xmm7     ;# fscal
    movapd xmm10, xmm9
    movapd xmm11, xmm9

    mulpd  xmm9,  [rsp + nb133_dxO] ;# fx/fy/fz
    mulpd  xmm10, [rsp + nb133_dyO]
    mulpd  xmm11, [rsp + nb133_dzO]

    ;# save j force temporarily
    movapd [rsp + nb133_fjx], xmm9
    movapd [rsp + nb133_fjy], xmm10
    movapd [rsp + nb133_fjz], xmm11
    
    ;# increment i O force
    addpd xmm9, [rsp + nb133_fixO]
    addpd xmm10, [rsp + nb133_fiyO]
    addpd xmm11, [rsp + nb133_fizO]
    movapd [rsp + nb133_fixO], xmm9
    movapd [rsp + nb133_fiyO], xmm10
    movapd [rsp + nb133_fizO], xmm11
    ;# finished O LJ interaction.


    ;# do H1, H2, and M interactions in parallel.
    ;# xmm0-xmm2 still contain j coordinates.        
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subpd xmm0, [rsp + nb133_ixH1]
    subpd xmm1, [rsp + nb133_iyH1]
    subpd xmm2, [rsp + nb133_izH1]
    subpd xmm3, [rsp + nb133_ixH2]
    subpd xmm4, [rsp + nb133_iyH2]
    subpd xmm5, [rsp + nb133_izH2]
    subpd xmm6, [rsp + nb133_ixM]
    subpd xmm7, [rsp + nb133_iyM]
    subpd xmm8, [rsp + nb133_izM]
    
	movapd [rsp + nb133_dxH1], xmm0
	movapd [rsp + nb133_dyH1], xmm1
	movapd [rsp + nb133_dzH1], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb133_dxH2], xmm3
	movapd [rsp + nb133_dyH2], xmm4
	movapd [rsp + nb133_dzH2], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb133_dxM], xmm6
	movapd [rsp + nb133_dyM], xmm7
	movapd [rsp + nb133_dzM], xmm8
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
		
	movapd  xmm9, [rsp + nb133_three]
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

	movapd  xmm15, [rsp + nb133_half]
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
		
	movapd  xmm1, [rsp + nb133_three]
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

	movapd  xmm15, [rsp + nb133_half]
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
    mulpd  xmm0, [rsp + nb133_qqH] 
    mulpd  xmm1, [rsp + nb133_qqH] 
    mulpd  xmm2, [rsp + nb133_qqM] 
    mulpd  xmm9, xmm0
    mulpd  xmm10, xmm1
    mulpd  xmm11, xmm2
    
    addpd xmm0, [rsp + nb133_vctot] 
    addpd xmm1, xmm2
    addpd xmm0, xmm1
    movapd [rsp + nb133_vctot], xmm0
    
    ;# move j forces to xmm0-xmm2
	mov   rdi, [rbp + nb133_faction]	
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
    addpd xmm0, [rsp + nb133_fjx]
    addpd xmm1, [rsp + nb133_fjy]
    addpd xmm2, [rsp + nb133_fjz]

	mulpd xmm7, [rsp + nb133_dxH1]
	mulpd xmm8, [rsp + nb133_dyH1]
	mulpd xmm9, [rsp + nb133_dzH1]
	mulpd xmm10, [rsp + nb133_dxH2]
	mulpd xmm11, [rsp + nb133_dyH2]
	mulpd xmm12, [rsp + nb133_dzH2]
	mulpd xmm13, [rsp + nb133_dxM]
	mulpd xmm14, [rsp + nb133_dyM]
	mulpd xmm15, [rsp + nb133_dzM]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb133_fixH1]
    addpd xmm8, [rsp + nb133_fiyH1]
    addpd xmm9, [rsp + nb133_fizH1]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb133_fixH2]
    addpd xmm11, [rsp + nb133_fiyH2]
    addpd xmm12, [rsp + nb133_fizH2]

    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    addpd xmm13, [rsp + nb133_fixM]
    addpd xmm14, [rsp + nb133_fiyM]
    addpd xmm15, [rsp + nb133_fizM]

    movapd [rsp + nb133_fixH1], xmm7
    movapd [rsp + nb133_fiyH1], xmm8
    movapd [rsp + nb133_fizH1], xmm9
    movapd [rsp + nb133_fixH2], xmm10
    movapd [rsp + nb133_fiyH2], xmm11
    movapd [rsp + nb133_fizH2], xmm12
    movapd [rsp + nb133_fixM], xmm13
    movapd [rsp + nb133_fiyM], xmm14
    movapd [rsp + nb133_fizM], xmm15

    ;# store back j forces from xmm0-xmm2
	movlpd [rdi + rax*8], xmm0
	movlpd [rdi + rax*8 + 8], xmm1
	movlpd [rdi + rax*8 + 16], xmm2
	movhpd [rdi + rbx*8], xmm0
	movhpd [rdi + rbx*8 + 8], xmm1
	movhpd [rdi + rbx*8 + 16], xmm2

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb133_innerk],  2
	jl   .nb133_checksingle
	jmp  .nb133_unroll_loop
.nb133_checksingle:	
	mov   edx, [rsp + nb133_innerk]
	and   edx, 1
	jnz  .nb133_dosingle
	jmp  .nb133_updateouterdata
.nb133_dosingle:
	mov   rdx, [rsp + nb133_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb133_innerjjnr],  4

	mov rsi, [rbp + nb133_charge]    ;# base of charge[] 

	xorpd xmm3, xmm3
	movlpd xmm3, [rsi + rax*8]
	movapd xmm4, xmm3
	mulsd  xmm3, [rsp + nb133_iqM]
	mulsd  xmm4, [rsp + nb133_iqH]

	movapd  [rsp + nb133_qqM], xmm3
	movapd  [rsp + nb133_qqH], xmm4
	
	mov rsi, [rbp + nb133_type]
	mov r8d, [rsi + rax*4]
	mov rsi, [rbp + nb133_vdwparam]
	shl r8d, 1	
	mov edi, [rsp + nb133_ntia]
	add r8d, edi

	movlpd xmm6, [rsi + r8*8]	;# c6a
	movhpd xmm6, [rsi + r8*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movapd [rsp + nb133_c6], xmm4
	movapd [rsp + nb133_c12], xmm6
	
	mov rsi, [rbp + nb133_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2  and xmm4-xmm6
	movlpd xmm4, [rsi + rax*8]
	movlpd xmm5, [rsi + rax*8 + 8]
	movlpd xmm6, [rsi + rax*8 + 16]
    movapd xmm0, xmm4
    movapd xmm1, xmm5
    movapd xmm2, xmm6

	;# calc dr 
	subsd xmm4, [rsp + nb133_ixO]
	subsd xmm5, [rsp + nb133_iyO]
	subsd xmm6, [rsp + nb133_izO]

	;# store dr 
	movapd [rsp + nb133_dxO], xmm4
	movapd [rsp + nb133_dyO], xmm5
	movapd [rsp + nb133_dzO], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 
	movapd [rsp + nb133_rsqO], xmm7
	
	;# move j coords to xmm4-xmm6 
	movapd xmm4, xmm0
	movapd xmm5, xmm1
	movapd xmm6, xmm2

	;# calc dr 
	subsd xmm4, [rsp + nb133_ixH1]
	subsd xmm5, [rsp + nb133_iyH1]
	subsd xmm6, [rsp + nb133_izH1]

	;# store dr 
	movapd [rsp + nb133_dxH1], xmm4
	movapd [rsp + nb133_dyH1], xmm5
	movapd [rsp + nb133_dzH1], xmm6
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
	subsd xmm3, [rsp + nb133_ixH2]
	subsd xmm4, [rsp + nb133_iyH2]
	subsd xmm5, [rsp + nb133_izH2]

	;# store dr 
	movapd [rsp + nb133_dxH2], xmm3
	movapd [rsp + nb133_dyH2], xmm4
	movapd [rsp + nb133_dzH2], xmm5
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
	subsd xmm4, [rsp + nb133_ixM]
	subsd xmm3, [rsp + nb133_iyM]
	subsd xmm2, [rsp + nb133_izM]

	;# store dr 
	movapd [rsp + nb133_dxM], xmm4
	movapd [rsp + nb133_dyM], xmm3
	movapd [rsp + nb133_dzM], xmm2

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
	movapd  xmm1, [rsp + nb133_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb133_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb133_three]
	subsd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb133_half] ;# rinv 
	movapd [rsp + nb133_rinvH1], xmm1
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [rsp + nb133_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb133_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb133_three]
	subsd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb133_half] ;# rinv 
	movapd [rsp + nb133_rinvH2], xmm1
	
	;# rsqM - seed in xmm2 
	cvtsd2ss xmm2, xmm4
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [rsp + nb133_three]
	mulsd   xmm2, xmm4	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb133_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb133_three]
	subsd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb133_half] ;# rinv 
	movapd [rsp + nb133_rinvM], xmm1

	;# rsqO - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movsd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movsd  xmm4, [rsp + nb133_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb133_half] ;# iter1 ( new lu) 

	movsd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movsd xmm4, [rsp + nb133_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb133_half] ;# rinv 
	movsd  xmm7, xmm4	;# rinvO in xmm7 
	
	movsd xmm4, [rsp + nb133_rsqO]
	movapd xmm0, xmm7
	;# LJ table interaction.
	mulsd xmm4, xmm7	;# xmm4=r 
	mulsd xmm4, [rsp + nb133_tsc]
	
	cvttsd2si ebx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, ebx
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	shl ebx, 3

	mov  rsi, [rbp + nb133_VFtab]

	;# dispersion 
	movlpd xmm4, [rsi + rbx*8]	;# Y1 	
	movhpd xmm4, [rsi + rbx*8 + 8]	;# Y1 F1 	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rbx*8 + 16]	;# G1
	movhpd xmm6, [rsi + rbx*8 + 24]	;# G1 H1 	
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# dispersion table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [rsp + nb133_two]	;# two*Heps2 
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movsd xmm4, [rsp + nb133_c6]
	mulsd  xmm7, xmm4	 ;# fijD 
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack Update Vvdwtot directly 
	addsd  xmm5, [rsp + nb133_Vvdwtot]
	xorpd  xmm3, xmm3
	mulsd  xmm7, [rsp + nb133_tsc]
	subsd  xmm3, xmm7
	movsd [rsp + nb133_fstmp], xmm3
	movsd [rsp + nb133_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [rsi + rbx*8 + 32]	;# Y1 	
	movhpd xmm4, [rsi + rbx*8 + 40]	;# Y1 F1 	

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rbx*8 + 48]	;# G1
	movhpd xmm6, [rsi + rbx*8 + 56]	;# G1 H1 	

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	
	;# table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [rsp + nb133_two]	;# two*Heps2 
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	
	movsd xmm4, [rsp + nb133_c12]
	mulsd  xmm7, xmm4 
	mulsd  xmm5, xmm4  
	
	addsd  xmm5, [rsp + nb133_Vvdwtot]
	movsd xmm3, [rsp + nb133_fstmp]
	mulsd  xmm7, [rsp + nb133_tsc]
	subsd  xmm3, xmm7
	movsd [rsp + nb133_Vvdwtot], xmm5

	mulsd  xmm3, xmm0
		
		
	movsd xmm0, [rsp + nb133_dxO]
	movsd xmm1, [rsp + nb133_dyO]
	movsd xmm2, [rsp + nb133_dzO]

	mov    rdi, [rbp + nb133_faction]
	mulsd  xmm0, xmm3
	mulsd  xmm1, xmm3
	mulsd  xmm2, xmm3

	;# update O forces 
	movapd xmm3, [rsp + nb133_fixO]
	movapd xmm4, [rsp + nb133_fiyO]
	movapd xmm7, [rsp + nb133_fizO]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [rsp + nb133_fixO], xmm3
	movsd [rsp + nb133_fiyO], xmm4
	movsd [rsp + nb133_fizO], xmm7
	;# update j forces with water O 
	movsd [rsp + nb133_fjx], xmm0
	movsd [rsp + nb133_fjy], xmm1
	movsd [rsp + nb133_fjz], xmm2

	;# H1 interactions
	movsd  xmm6, [rsp + nb133_rinvH1] 
	movsd  xmm4, xmm6
	mulsd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	mulsd   xmm6, [rsp + nb133_qqH] ;# vcoul 
	mulsd   xmm4, xmm6    ;# fscal
	addsd  xmm6, [rsp + nb133_vctot]
	movsd [rsp + nb133_vctot], xmm6

	movapd xmm0, [rsp + nb133_dxH1]
	movapd xmm1, [rsp + nb133_dyH1]
	movapd xmm2, [rsp + nb133_dzH1]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [rsp + nb133_fixH1]
	movapd xmm4, [rsp + nb133_fiyH1]
	movapd xmm7, [rsp + nb133_fizH1]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [rsp + nb133_fixH1], xmm3
	movsd [rsp + nb133_fiyH1], xmm4
	movsd [rsp + nb133_fizH1], xmm7
	;# update j forces with water H1 
	addsd  xmm0, [rsp + nb133_fjx]
	addsd  xmm1, [rsp + nb133_fjy]
	addsd  xmm2, [rsp + nb133_fjz]
	movsd [rsp + nb133_fjx], xmm0
	movsd [rsp + nb133_fjy], xmm1
	movsd [rsp + nb133_fjz], xmm2

	;# H2 interactions 
	movsd  xmm6, [rsp + nb133_rinvH2] 
	movsd  xmm4, xmm6
	mulsd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	mulsd   xmm6, [rsp + nb133_qqH] ;# vcoul 
	mulsd   xmm4, xmm6    ;# fscal
	addsd  xmm6, [rsp + nb133_vctot]
	movsd [rsp + nb133_vctot], xmm6

	movapd xmm0, [rsp + nb133_dxH2]
	movapd xmm1, [rsp + nb133_dyH2]
	movapd xmm2, [rsp + nb133_dzH2]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [rsp + nb133_fixH2]
	movapd xmm4, [rsp + nb133_fiyH2]
	movapd xmm7, [rsp + nb133_fizH2]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [rsp + nb133_fixH2], xmm3
	movsd [rsp + nb133_fiyH2], xmm4
	movsd [rsp + nb133_fizH2], xmm7
	;# update j forces with water H2 
	addsd  xmm0, [rsp + nb133_fjx]
	addsd  xmm1, [rsp + nb133_fjy]
	addsd  xmm2, [rsp + nb133_fjz]
	movsd [rsp + nb133_fjx], xmm0
	movsd [rsp + nb133_fjy], xmm1
	movsd [rsp + nb133_fjz], xmm2

	;# M interactions 
	movsd  xmm6, [rsp + nb133_rinvM] 
	movsd  xmm4, xmm6
	mulsd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	mulsd   xmm6, [rsp + nb133_qqM] ;# vcoul 
	mulsd   xmm4, xmm6    ;# fscal
	addsd  xmm6, [rsp + nb133_vctot]
	movsd [rsp + nb133_vctot], xmm6

	movapd xmm0, [rsp + nb133_dxM]
	movapd xmm1, [rsp + nb133_dyM]
	movapd xmm2, [rsp + nb133_dzM]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update M forces 
	movapd xmm3, [rsp + nb133_fixM]
	movapd xmm4, [rsp + nb133_fiyM]
	movapd xmm7, [rsp + nb133_fizM]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [rsp + nb133_fixM], xmm3
	movsd [rsp + nb133_fiyM], xmm4
	movsd [rsp + nb133_fizM], xmm7

	mov rdi, [rbp + nb133_faction]
	;# update j forces 
	addsd  xmm0, [rsp + nb133_fjx]
	addsd  xmm1, [rsp + nb133_fjy]
	addsd  xmm2, [rsp + nb133_fjz]
	movlpd xmm3, [rdi + rax*8]
	movlpd xmm4, [rdi + rax*8 + 8]
	movlpd xmm5, [rdi + rax*8 + 16]
	addsd xmm3, xmm0
	addsd xmm4, xmm1
	addsd xmm5, xmm2
	movlpd [rdi + rax*8], xmm3
	movlpd [rdi + rax*8 + 8], xmm4
	movlpd [rdi + rax*8 + 16], xmm5

.nb133_updateouterdata:
	mov   ecx, [rsp + nb133_ii3]
	mov   rdi, [rbp + nb133_faction]
	mov   rsi, [rbp + nb133_fshift]
	mov   edx, [rsp + nb133_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb133_fixO]
	movapd xmm1, [rsp + nb133_fiyO]
	movapd xmm2, [rsp + nb133_fizO]

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
	movapd xmm0, [rsp + nb133_fixH1]
	movapd xmm1, [rsp + nb133_fiyH1]
	movapd xmm2, [rsp + nb133_fizH1]

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
	movapd xmm0, [rsp + nb133_fixH2]
	movapd xmm1, [rsp + nb133_fiyH2]
	movapd xmm2, [rsp + nb133_fizH2]

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
	movapd xmm0, [rsp + nb133_fixM]
	movapd xmm1, [rsp + nb133_fiyM]
	movapd xmm2, [rsp + nb133_fizM]

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
	mov esi, [rsp + nb133_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb133_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb133_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb133_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb133_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb133_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
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



.globl nb_kernel133nf_x86_64_sse2
.globl _nb_kernel133nf_x86_64_sse2
nb_kernel133nf_x86_64_sse2:	
_nb_kernel133nf_x86_64_sse2:	
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
	;# bottom of stack is cache-aligned for sse2 use 
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
.equiv          nb133nf_iqH,            192
.equiv          nb133nf_iqM,            208
.equiv          nb133nf_qqH,            224
.equiv          nb133nf_qqM,            240
.equiv          nb133nf_c6,             256
.equiv          nb133nf_c12,            272
.equiv          nb133nf_vctot,          288
.equiv          nb133nf_Vvdwtot,        304
.equiv          nb133nf_half,           320
.equiv          nb133nf_three,          336
.equiv          nb133nf_tsc,            352
.equiv          nb133nf_rinvH1,         368
.equiv          nb133nf_rinvH2,         384
.equiv          nb133nf_rinvM,          400
.equiv          nb133nf_krsqH1,         416
.equiv          nb133nf_krsqH2,         432
.equiv          nb133nf_krsqM,          448
.equiv          nb133nf_rsqO,           464
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
	movsd xmm0, [rsi]
	movsd [rsp + nb133nf_facel], xmm0

	mov rax, [rbp + nb133nf_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb133nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb133nf_half], eax
	mov [rsp + nb133nf_half+4], ebx
	movsd xmm1, [rsp + nb133nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb133nf_half], xmm1
	movapd [rsp + nb133nf_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb133nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb133nf_charge]
	movsd xmm3, [rdx + rbx*8 + 8]	
	movsd xmm4, [rdx + rbx*8 + 24]	

	movsd xmm5, [rsp + nb133nf_facel]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb133nf_iqH], xmm3
	movapd [rsp + nb133nf_iqM], xmm4
	
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
	mov   ebx, [rax+rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb133nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb133nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb133nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb133nf_pos]    ;# eax = base of pos[]  
	mov   [rsp + nb133nf_ii3], ebx

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
	movapd [rsp + nb133nf_ixO], xmm3
	movapd [rsp + nb133nf_iyO], xmm4
	movapd [rsp + nb133nf_izO], xmm5
	movapd [rsp + nb133nf_ixH1], xmm6
	movapd [rsp + nb133nf_iyH1], xmm7

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
	movapd [rsp + nb133nf_izH1], xmm6
	movapd [rsp + nb133nf_ixH2], xmm0
	movapd [rsp + nb133nf_iyH2], xmm1
	movapd [rsp + nb133nf_izH2], xmm2
	movapd [rsp + nb133nf_ixM], xmm3
	movapd [rsp + nb133nf_iyM], xmm4
	movapd [rsp + nb133nf_izM], xmm5

	;# clear vctot
	xorpd xmm4, xmm4
	movapd [rsp + nb133nf_vctot], xmm4
	movapd [rsp + nb133nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb133nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb133nf_pos]
	mov   rdi, [rbp + nb133nf_faction]	
	mov   rax, [rsp + nb133nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb133nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb133nf_ninner]
	mov   [rsp + nb133nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb133nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb133nf_unroll_loop
	jmp   .nb133nf_checksingle
.nb133nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb133nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]

	add qword ptr [rsp + nb133nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb133nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [rsp + nb133nf_iqM]
	mulpd  xmm4, [rsp + nb133nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx

	movapd  [rsp + nb133nf_qqM], xmm3
	movapd  [rsp + nb133nf_qqH], xmm4
	
	mov rsi, [rbp + nb133nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov rsi, [rbp + nb133nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	mov edi, [rsp + nb133nf_ntia]
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
	movapd [rsp + nb133nf_c6], xmm4
	movapd [rsp + nb133nf_c12], xmm6
	
	mov rsi, [rbp + nb133nf_pos]       ;# base of pos[] 

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
	movapd xmm4, [rsp + nb133nf_ixO]
	movapd xmm5, [rsp + nb133nf_iyO]
	movapd xmm6, [rsp + nb133nf_izO]

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
	movapd xmm4, [rsp + nb133nf_ixH1]
	movapd xmm5, [rsp + nb133nf_iyH1]
	movapd xmm6, [rsp + nb133nf_izH1]

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
	movapd xmm3, [rsp + nb133nf_ixH2]
	movapd xmm4, [rsp + nb133nf_iyH2]
	movapd xmm5, [rsp + nb133nf_izH2]

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
	movapd xmm3, [rsp + nb133nf_iyM]
	movapd xmm4, [rsp + nb133nf_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [rsp + nb133nf_ixM]
	subpd  xmm2, xmm0	

	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 
	movapd [rsp + nb133nf_rsqO], xmm7
	
	;# start with rsqH1 - put seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2
	
	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [rsp + nb133nf_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [rsp + nb133nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb133nf_three]
	subpd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [rsp + nb133nf_half] ;# rinv 
	movapd  [rsp + nb133nf_rinvH1], xmm1	

	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [rsp + nb133nf_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [rsp + nb133nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb133nf_three]
	subpd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [rsp + nb133nf_half] ;# rinv 
	movapd  [rsp + nb133nf_rinvH2], xmm1	
	
	;# rsqM - seed in xmm2 
	cvtpd2ps xmm2, xmm4	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [rsp + nb133nf_three]
	mulpd   xmm2, xmm4	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [rsp + nb133nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb133nf_three]
	subpd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [rsp + nb133nf_half] ;# rinv 
	movapd  [rsp + nb133nf_rinvM], xmm1	

		
	;# rsqO - put seed in xmm2 
	cvtpd2ps xmm2, xmm7	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [rsp + nb133nf_three]
	mulpd   xmm2, xmm7	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [rsp + nb133nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb133nf_three]
	subpd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [rsp + nb133nf_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	
	
	movapd xmm4, [rsp + nb133nf_rsqO]
	movapd xmm0, xmm7
	;# LJ table interaction.
	mulpd xmm4, xmm7	;# xmm4=r 
	mulpd xmm4, [rsp + nb133nf_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 3		;# idx *= 8 
	
	mov  rsi, [rbp + nb133nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6

	;# dispersion 
	movlpd xmm4, [rsi + rax*8]	;# Y1 	
	movlpd xmm3, [rsi + rbx*8]	;# Y2 
	movhpd xmm4, [rsi + rax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [rsi + rbx*8 + 8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rax*8 + 16]	;# G1
	movlpd xmm3, [rsi + rbx*8 + 16]	;# G2
	movhpd xmm6, [rsi + rax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [rsi + rbx*8 + 24]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# dispersion table ready, in xmm4-xmm7 	
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [rsp + nb133nf_c6]
	mulpd  xmm5, xmm4	 ;# Vvdw6 

	;# Update Vvdwtot directly 
	addpd  xmm5, [rsp + nb133nf_Vvdwtot]
	movapd [rsp + nb133nf_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [rsi + rax*8 + 32]	;# Y1 	
	movlpd xmm3, [rsi + rbx*8 + 32]	;# Y2 
	movhpd xmm4, [rsi + rax*8 + 40]	;# Y1 F1 	
	movhpd xmm3, [rsi + rbx*8 + 40]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rax*8 + 48]	;# G1
	movlpd xmm3, [rsi + rbx*8 + 48]	;# G2
	movhpd xmm6, [rsi + rax*8 + 56]	;# G1 H1 	
	movhpd xmm3, [rsi + rbx*8 + 56]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	
	;# table ready, in xmm4-xmm7 	
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	
	movapd xmm4, [rsp + nb133nf_c12]
	mulpd  xmm5, xmm4  
	
	addpd  xmm5, [rsp + nb133nf_Vvdwtot]
	movapd [rsp + nb133nf_Vvdwtot], xmm5

	;# H1/H2/M interactions 
	movapd  xmm6, [rsp + nb133nf_rinvH1] 
	addpd   xmm6, [rsp + nb133nf_rinvH2] 
	movapd  xmm7, [rsp + nb133nf_rinvM] 
	mulpd   xmm6, [rsp + nb133nf_qqH]
	mulpd   xmm7, [rsp + nb133nf_qqM]
	addpd   xmm6, xmm7
	addpd   xmm6, [rsp + nb133nf_vctot]
	movapd  [rsp + nb133nf_vctot], xmm6
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb133nf_innerk],  2
	jl   .nb133nf_checksingle
	jmp  .nb133nf_unroll_loop
.nb133nf_checksingle:	
	mov   edx, [rsp + nb133nf_innerk]
	and   edx, 1
	jnz  .nb133nf_dosingle
	jmp  .nb133nf_updateouterdata
.nb133nf_dosingle:
	mov   rdx, [rsp + nb133nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb133nf_innerjjnr],  4

	mov rsi, [rbp + nb133nf_charge]    ;# base of charge[] 

	xorpd xmm3, xmm3
	movlpd xmm3, [rsi + rax*8]
	movapd xmm4, xmm3
	mulsd  xmm3, [rsp + nb133nf_iqM]
	mulsd  xmm4, [rsp + nb133nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 

	movapd  [rsp + nb133nf_qqM], xmm3
	movapd  [rsp + nb133nf_qqH], xmm4
	
	mov rsi, [rbp + nb133nf_type]
	mov eax, [rsi + rax*4]
	mov rsi, [rbp + nb133nf_vdwparam]
	shl eax, 1	
	mov edi, [rsp + nb133nf_ntia]
	add eax, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [rsp + nb133nf_c6], xmm4
	movapd [rsp + nb133nf_c12], xmm6
	
	mov rsi, [rbp + nb133nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [rsp + nb133nf_ixO]
	movapd xmm5, [rsp + nb133nf_iyO]
	movapd xmm6, [rsp + nb133nf_izO]

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
	movapd [rsp + nb133nf_rsqO], xmm7
	
	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [rsp + nb133nf_ixH1]
	movapd xmm5, [rsp + nb133nf_iyH1]
	movapd xmm6, [rsp + nb133nf_izH1]

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
	movapd xmm3, [rsp + nb133nf_ixH2]
	movapd xmm4, [rsp + nb133nf_iyH2]
	movapd xmm5, [rsp + nb133nf_izH2]

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
	movapd xmm3, [rsp + nb133nf_iyM]
	movapd xmm4, [rsp + nb133nf_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [rsp + nb133nf_ixM]
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
	movapd  xmm1, [rsp + nb133nf_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb133nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb133nf_three]
	subsd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb133nf_half] ;# rinv 
	movapd [rsp + nb133nf_rinvH1], xmm1
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [rsp + nb133nf_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb133nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb133nf_three]
	subsd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb133nf_half] ;# rinv 
	movapd [rsp + nb133nf_rinvH2], xmm1
	
	;# rsqM - seed in xmm2 
	cvtsd2ss xmm2, xmm4
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [rsp + nb133nf_three]
	mulsd   xmm2, xmm4	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb133nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb133nf_three]
	subsd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb133nf_half] ;# rinv 
	movapd [rsp + nb133nf_rinvM], xmm1

	;# rsqO - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movsd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movsd  xmm4, [rsp + nb133nf_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb133nf_half] ;# iter1 ( new lu) 

	movsd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movsd xmm4, [rsp + nb133nf_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb133nf_half] ;# rinv 
	movsd  xmm7, xmm4	;# rinvO in xmm7 
	
	movsd xmm4, [rsp + nb133nf_rsqO]
	movapd xmm0, xmm7
	;# LJ table interaction.
	mulsd xmm4, xmm7	;# xmm4=r 
	mulsd xmm4, [rsp + nb133nf_tsc]
	
	cvttsd2si ebx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, ebx
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	shl ebx, 3

	mov  rsi, [rbp + nb133nf_VFtab]

	;# dispersion 
	movlpd xmm4, [rsi + rbx*8]	;# Y1 	
	movhpd xmm4, [rsi + rbx*8 + 8]	;# Y1 F1 	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rbx*8 + 16]	;# G1
	movhpd xmm6, [rsi + rbx*8 + 24]	;# G1 H1 	
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# dispersion table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movsd xmm4, [rsp + nb133nf_c6]
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack Update Vvdwtot directly 
	addsd  xmm5, [rsp + nb133nf_Vvdwtot]
	movsd [rsp + nb133nf_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [rsi + rbx*8 + 32]	;# Y1 	
	movhpd xmm4, [rsi + rbx*8 + 40]	;# Y1 F1 	

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rbx*8 + 48]	;# G1
	movhpd xmm6, [rsi + rbx*8 + 56]	;# G1 H1 	

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	
	;# table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	
	movsd xmm4, [rsp + nb133nf_c12]
	mulsd  xmm5, xmm4  
	
	addsd  xmm5, [rsp + nb133nf_Vvdwtot]
	movsd [rsp + nb133nf_Vvdwtot], xmm5

	;# H1/H2/M interactions 
	movsd  xmm6, [rsp + nb133nf_rinvH1] 
	addsd  xmm6, [rsp + nb133nf_rinvH2] 
	movsd  xmm7, [rsp + nb133nf_rinvM] 
	mulsd  xmm6, [rsp + nb133nf_qqH]
	mulsd  xmm7, [rsp + nb133nf_qqM]
	addsd  xmm6, xmm7
	addsd  xmm6, [rsp + nb133nf_vctot]
	movsd  [rsp + nb133nf_vctot], xmm6
	
.nb133nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb133nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb133nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb133nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb133nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb133nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb133nf_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
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

