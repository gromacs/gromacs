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

	
.globl nb_kernel213_x86_64_sse
.globl _nb_kernel213_x86_64_sse
nb_kernel213_x86_64_sse:	
_nb_kernel213_x86_64_sse:	
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
	;# bottom of stack is cache-aligned for sse use 
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
.equiv          nb213_iqM,              192
.equiv          nb213_iqH,              208
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
.equiv          nb213_qqM,              416
.equiv          nb213_qqH,              432
.equiv          nb213_rinvH1,           448
.equiv          nb213_rinvH2,           464
.equiv          nb213_rinvM,            480
.equiv          nb213_two,              496
.equiv          nb213_c6,               512
.equiv          nb213_c12,              528
.equiv          nb213_six,              544
.equiv          nb213_twelve,           560
.equiv          nb213_krf,              576
.equiv          nb213_crf,              592
.equiv          nb213_krsqH1,           608
.equiv          nb213_krsqH2,           624
.equiv          nb213_krsqM,            640
.equiv          nb213_vctot,            656
.equiv          nb213_Vvdwtot,          672
.equiv          nb213_fixO,             688
.equiv          nb213_fiyO,             704
.equiv          nb213_fizO,             720
.equiv          nb213_fixH1,            736
.equiv          nb213_fiyH1,            752
.equiv          nb213_fizH1,            768
.equiv          nb213_fixH2,            784
.equiv          nb213_fiyH2,            800
.equiv          nb213_fizH2,            816
.equiv          nb213_fixM,             832
.equiv          nb213_fiyM,             848
.equiv          nb213_fizM,             864
.equiv          nb213_fjx,              880
.equiv          nb213_fjy,              896
.equiv          nb213_fjz,              912
.equiv          nb213_half,             928
.equiv          nb213_three,            944
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
    
	emms
	sub rsp, 1056		;# local variable stack space (n*16+8)
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
	movss xmm0, [rsi]
	movss [rsp + nb213_facel], xmm0


	mov rsi, [rbp + nb213_argkrf]
	mov rdi, [rbp + nb213_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb213_krf], xmm1
	movaps [rsp + nb213_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb213_half], eax
	movss xmm1, [rsp + nb213_half]
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
	movaps [rsp + nb213_half],  xmm1
	movaps [rsp + nb213_two],  xmm2
	movaps [rsp + nb213_three],  xmm3
	movaps [rsp + nb213_six],  xmm4
	movaps [rsp + nb213_twelve],  xmm5
	
	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb213_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb213_charge]
	movss xmm4, [rdx + rbx*4 + 4]	
	movss xmm3, [rdx + rbx*4 + 12]	
	mov rsi, [rbp + nb213_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb213_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb213_iqM], xmm3
	movaps [rsp + nb213_iqH], xmm4
	
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
        add   rbx, 1                           ;# rbx=nn1=nn0+10
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
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb213_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb213_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb213_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2	
	movaps xmm6, xmm0
	movaps xmm7, xmm1
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb213_pos]	;# rax = base of pos[]  
	mov   [rsp + nb213_ii3], ebx

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
	movaps [rsp + nb213_ixO], xmm3
	movaps [rsp + nb213_iyO], xmm4
	movaps [rsp + nb213_izO], xmm5
	movaps [rsp + nb213_ixH1], xmm6
	movaps [rsp + nb213_iyH1], xmm7

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
	movaps [rsp + nb213_izH1], xmm6
	movaps [rsp + nb213_ixH2], xmm0
	movaps [rsp + nb213_iyH2], xmm1
	movaps [rsp + nb213_izH2], xmm2
	movaps [rsp + nb213_ixM], xmm3
	movaps [rsp + nb213_iyM], xmm4
	movaps [rsp + nb213_izM], xmm5
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb213_vctot], xmm4
	movaps [rsp + nb213_Vvdwtot], xmm4
	movaps [rsp + nb213_fixO], xmm4
	movaps [rsp + nb213_fiyO], xmm4
	movaps [rsp + nb213_fizO], xmm4
	movaps [rsp + nb213_fixH1], xmm4
	movaps [rsp + nb213_fiyH1], xmm4
	movaps [rsp + nb213_fizH1], xmm4
	movaps [rsp + nb213_fixH2], xmm4
	movaps [rsp + nb213_fiyH2], xmm4
	movaps [rsp + nb213_fizH2], xmm4
	movaps [rsp + nb213_fixM], xmm4
	movaps [rsp + nb213_fiyM], xmm4
	movaps [rsp + nb213_fizM], xmm4
	
	mov   rax, [rsp + nb213_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb213_pos]
	mov   rdi, [rbp + nb213_faction]	
	mov   rax, [rsp + nb213_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb213_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb213_ninner]
	mov   [rsp + nb213_ninner], ecx
	add   edx, 0
	mov   [rsp + nb213_innerk], edx	;# number of innerloop atoms 
	jge   .nb213_unroll_loop
	jmp   .nb213_odd_inner
.nb213_unroll_loop:
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb213_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb213_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb213_charge]	;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	 	;# and in xmm4 
	mulps  xmm3, [rsp + nb213_iqM]
	mulps  xmm4, [rsp + nb213_iqH]

	movaps  [rsp + nb213_qqM], xmm3
	movaps  [rsp + nb213_qqH], xmm4
	
	mov rsi, [rbp + nb213_type]
	mov r8d, [rsi + rax*4]
	mov r9d, [rsi + rbx*4]
	mov r10d, [rsi + rcx*4]
	mov r11d, [rsi + rdx*4]
	mov rsi, [rbp + nb213_vdwparam]
	shl r8d, 1	
	shl r9d, 1	
	shl r10d, 1	
	shl r11d, 1	
	mov edi, [rsp + nb213_ntia]
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
	
	movaps [rsp + nb213_c6], xmm4
	movaps [rsp + nb213_c12], xmm6

	mov rsi, [rbp + nb213_pos]   	;# base of pos[] 

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
    
    subps xmm3, [rsp + nb213_ixO]
    subps xmm4, [rsp + nb213_iyO]
    subps xmm5, [rsp + nb213_izO]
    
    movaps xmm13, xmm3
    movaps xmm14, xmm4
    movaps xmm15, xmm5
    
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm3, xmm4
	addps  xmm3, xmm5

    ;# calc 1/rsq
    rcpps xmm5, xmm3
    movaps xmm4, [rsp + nb213_two]
    mulps xmm3, xmm5
    subps xmm4, xmm3
    mulps xmm4, xmm5        ;# xmm4=rinvsq

    movaps xmm3, xmm4       ;# rinvsq
    mulps  xmm4, xmm4       ;# rinv4
    mulps  xmm4, xmm3       ;# rinv6
    movaps xmm5, xmm4     
    mulps  xmm5, xmm5       ;# rinv12
    mulps  xmm4, [rsp + nb213_c6]
    mulps  xmm5, [rsp + nb213_c12]
    movaps xmm6, xmm5
    subps  xmm6, xmm4  ;# Vvdw=vvdw12-vvdw6
    mulps  xmm4, [rsp + nb213_six]
    mulps  xmm5, [rsp + nb213_twelve]
    subps  xmm5, xmm4
    mulps  xmm3, xmm5   ;# fscal
    
    addps  xmm6, [rsp + nb213_Vvdwtot]
    movaps [rsp + nb213_Vvdwtot], xmm6
    
    mulps  xmm13, xmm3 ;# fx
    mulps  xmm14, xmm3 ;# fy
    mulps  xmm15, xmm3 ;# fz

    ;# save j force temporarily
    movaps [rsp + nb213_fjx], xmm13
    movaps [rsp + nb213_fjy], xmm14
    movaps [rsp + nb213_fjz], xmm15
    
    ;# increment i O force
    addps xmm13, [rsp + nb213_fixO]
    addps xmm14, [rsp + nb213_fiyO]
    addps xmm15, [rsp + nb213_fizO]
    movaps [rsp + nb213_fixO], xmm13
    movaps [rsp + nb213_fiyO], xmm14
    movaps [rsp + nb213_fizO], xmm15    
    ;# finished O LJ interaction.


    ;# do H1, H2, and M interactions in parallel.
    ;# xmm0-xmm2 still contain j coordinates.        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + nb213_ixH1]
    subps xmm1, [rsp + nb213_iyH1]
    subps xmm2, [rsp + nb213_izH1]
    subps xmm3, [rsp + nb213_ixH2]
    subps xmm4, [rsp + nb213_iyH2]
    subps xmm5, [rsp + nb213_izH2]
    subps xmm6, [rsp + nb213_ixM]
    subps xmm7, [rsp + nb213_iyM]
    subps xmm8, [rsp + nb213_izM]
    
	movaps [rsp + nb213_dxH1], xmm0
	movaps [rsp + nb213_dyH1], xmm1
	movaps [rsp + nb213_dzH1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb213_dxH2], xmm3
	movaps [rsp + nb213_dyH2], xmm4
	movaps [rsp + nb213_dzH2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb213_dxM], xmm6
	movaps [rsp + nb213_dyM], xmm7
	movaps [rsp + nb213_dzM], xmm8
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
		
	movaps  xmm9, [rsp + nb213_three]
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

	movaps  xmm4, [rsp + nb213_half]
	mulps   xmm9, xmm4  ;# rinvH1
	mulps   xmm10, xmm4 ;# rinvH2
    mulps   xmm11, xmm4 ;# rinvM
	
	;# interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movaps xmm1, xmm9 ;# copy of rinv
    movaps xmm4, xmm10
    movaps xmm7, xmm11
    movaps xmm2, [rsp + nb213_krf]    
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
    movaps xmm14, [rsp + nb213_crf]
    subps  xmm2, xmm14   ;# rinv+krsq-crf
    subps  xmm5, xmm14
    subps  xmm8, xmm14
    movaps xmm12, [rsp + nb213_qqH]
    movaps xmm13, [rsp + nb213_qqM]    
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
    addps  xmm2, [rsp + nb213_vctot]
    addps  xmm5, xmm8
    addps  xmm2, xmm5
    movaps [rsp + nb213_vctot], xmm2
    
    mulps  xmm1, xmm9   ;# fscal
    mulps  xmm4, xmm10
    mulps  xmm7, xmm11

	mov   rdi, [rbp + nb213_faction]	
	;# move j  forces to local temp variables 
    movlps xmm9, [rdi + rax*4] ;# jxa jya  -   -
    movlps xmm10, [rdi + rcx*4] ;# jxc jyc  -   -
    movhps xmm9, [rdi + rbx*4] ;# jxa jya jxb jyb 
    movhps xmm10, [rdi + rdx*4] ;# jxc jyc jxd jyd 

    movss  xmm11, [rdi + rax*4 + 8] ;# jza  -  -  -
    movss  xmm12, [rdi + rcx*4 + 8] ;# jzc  -  -  -
    movss  xmm6,  [rdi + rbx*4 + 8] ;# jzb
    movss  xmm8,  [rdi + rdx*4 + 8] ;# jzd
    movlhps xmm11, xmm6 ;# jza  -  jzb  -
    movlhps xmm12, xmm8 ;# jzc  -  jzd -
    
    shufps xmm11, xmm12,  136  ;# 10001000 => jza jzb jzc jzd

    ;# xmm9: jxa jya jxb jyb 
    ;# xmm10: jxc jyc jxd jyd
    ;# xmm11: jza jzb jzc jzd

    movaps xmm0, xmm1
    movaps xmm2, xmm1
    movaps xmm3, xmm4
    movaps xmm5, xmm4
    movaps xmm6, xmm7
    movaps xmm8, xmm7

	mulps xmm0, [rsp + nb213_dxH1]
	mulps xmm1, [rsp + nb213_dyH1]
	mulps xmm2, [rsp + nb213_dzH1]
	mulps xmm3, [rsp + nb213_dxH2]
	mulps xmm4, [rsp + nb213_dyH2]
	mulps xmm5, [rsp + nb213_dzH2]
	mulps xmm6, [rsp + nb213_dxM]
	mulps xmm7, [rsp + nb213_dyM]
	mulps xmm8, [rsp + nb213_dzM]

    ;# fetch forces from O interaction
    movaps xmm13, [rsp + nb213_fjx]
    movaps xmm14, [rsp + nb213_fjy]
    addps  xmm11, [rsp + nb213_fjz]

    addps xmm13, xmm0
    addps xmm14, xmm1
    addps xmm11, xmm2
    addps xmm0, [rsp + nb213_fixH1]
    addps xmm1, [rsp + nb213_fiyH1]
    addps xmm2, [rsp + nb213_fizH1]

    addps xmm13,  xmm3
    addps xmm14, xmm4
    addps xmm11, xmm5
    addps xmm3, [rsp + nb213_fixH2]
    addps xmm4, [rsp + nb213_fiyH2]
    addps xmm5, [rsp + nb213_fizH2]

    addps xmm13, xmm6
    addps xmm14, xmm7
    addps xmm11, xmm8
    addps xmm6, [rsp + nb213_fixM]
    addps xmm7, [rsp + nb213_fiyM]
    addps xmm8, [rsp + nb213_fizM]

    movaps [rsp + nb213_fixH1], xmm0
    movaps [rsp + nb213_fiyH1], xmm1
    movaps [rsp + nb213_fizH1], xmm2
    movaps [rsp + nb213_fixH2], xmm3
    movaps [rsp + nb213_fiyH2], xmm4
    movaps [rsp + nb213_fizH2], xmm5
    movaps [rsp + nb213_fixM], xmm6
    movaps [rsp + nb213_fiyM], xmm7
    movaps [rsp + nb213_fizM], xmm8
    
    ;# xmm9 = fjx
    ;# xmm10 = fjy
    ;# xmm11 = fjz
    movaps xmm15, xmm13
    unpcklps xmm13, xmm14
    unpckhps xmm15, xmm14
    
    addps xmm9, xmm13
    addps xmm10, xmm15

    movhlps  xmm12, xmm11 ;# fjzc fjzd
    
    movlps [rdi + rax*4], xmm9
    movhps [rdi + rbx*4], xmm9
    movlps [rdi + rcx*4], xmm10
    movhps [rdi + rdx*4], xmm10
    movss  [rdi + rax*4 + 8], xmm11
    movss  [rdi + rcx*4 + 8], xmm12
    shufps xmm11, xmm11, 1
    shufps xmm12, xmm12, 1
    movss  [rdi + rbx*4 + 8], xmm11
    movss  [rdi + rdx*4 + 8], xmm12

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb213_innerk],  4
	jl    .nb213_odd_inner
	jmp   .nb213_unroll_loop
.nb213_odd_inner:	
	add dword ptr [rsp + nb213_innerk],  4
	jnz   .nb213_odd_loop
	jmp   .nb213_updateouterdata
.nb213_odd_loop:
	mov   rdx, [rsp + nb213_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb213_innerjjnr],  4	

 	xorps xmm4, xmm4  	;# clear reg.
	movss xmm4, [rsp + nb213_iqM]
	mov rsi, [rbp + nb213_charge] 
	movhps xmm4, [rsp + nb213_iqH]  ;# [qM  0  qH  qH] 
	shufps xmm4, xmm4, 41	;# [0 qH qH qM]

	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb213_qqM], xmm3	;# use dummy qq for storage 
	
	xorps xmm6, xmm6
	mov rsi, [rbp + nb213_type]
	mov ebx, [rsi + rax*4]
	mov rsi, [rbp + nb213_vdwparam]
	shl ebx, 1	
	add ebx, [rsp + nb213_ntia]
	movlps xmm6, [rsi + rbx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 252  ;# 11111100
	shufps xmm7, xmm7, 253  ;# 11111101
	movaps [rsp + nb213_c6], xmm6
	movaps [rsp + nb213_c12], xmm7
	
	mov rsi, [rbp + nb213_pos]
	lea rax, [rax + rax*2]  

	movss xmm0, [rsp + nb213_ixO]
	movss xmm1, [rsp + nb213_iyO]
	movss xmm2, [rsp + nb213_izO]
	movss xmm3, [rsp + nb213_ixH1]
	movss xmm4, [rsp + nb213_iyH1]
	movss xmm5, [rsp + nb213_izH1]
	unpcklps xmm0, [rsp + nb213_ixH2] 	;# ixO ixH2 - -
	unpcklps xmm1, [rsp + nb213_iyH2]  	;# iyO iyH2 - -
	unpcklps xmm2, [rsp + nb213_izH2]	;# izO izH2 - -
	unpcklps xmm3, [rsp + nb213_ixM] 	;# ixH1 ixM - -
	unpcklps xmm4, [rsp + nb213_iyM]  	;# iyH1 iyM - -
	unpcklps xmm5, [rsp + nb213_izM]	;# izH1 izM - -
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
	movaps [rsp + nb213_dxO], xmm3
	movaps [rsp + nb213_dyO], xmm4
	movaps [rsp + nb213_dzO], xmm5

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	;# rsq in xmm4 
	movaps xmm0, xmm4
	mulps xmm0, [rsp + nb213_krf]
	movaps [rsp + nb213_krsqM], xmm0
	
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb213_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb213_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv

	
	movaps xmm4, xmm0
	movaps xmm1, xmm0
	movaps xmm5, xmm0
	mulps  xmm4, xmm4	;# xmm1=rinv, xmm4=rinvsq
	movaps xmm3, [rsp + nb213_krsqM]
	addps  xmm5, xmm3	;# xmm0=rinv+ krsq 
	subps  xmm5, [rsp + nb213_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm3, [rsp + nb213_two]
	subps  xmm1, xmm3	;# xmm1=rinv-2*krsq
	movaps xmm7, xmm5
	mulps  xmm7, [rsp + nb213_qqM]	;# xmm0=vcoul 
	mulps  xmm1, [rsp + nb213_qqM] 	;# xmm1=coul part of fs 
	movaps xmm6, xmm1
	
	addps  xmm7, [rsp + nb213_vctot]	
	movaps [rsp + nb213_vctot], xmm7

	movaps xmm1, xmm0
	mulps  xmm1, xmm1
	movaps xmm2, xmm1
	mulss  xmm1, xmm1
	mulss  xmm1, xmm2	;# xmm1=rinvsix
	xorps  xmm4, xmm4
	movss  xmm4, xmm1
	mulss  xmm4, xmm4	;# xmm4=rinvtwelve 
	mulss  xmm1, [rsp + nb213_c6]
	mulss  xmm4, [rsp + nb213_c12]
	movaps xmm3, xmm4
	subss  xmm3, xmm1	;# xmm3=Vvdw12-Vvdw6 
	mulss  xmm1, [rsp + nb213_six]
	mulss  xmm4, [rsp + nb213_twelve]
	subss  xmm4, xmm1
	addss  xmm3, [rsp + nb213_Vvdwtot]
	movss  [rsp + nb213_Vvdwtot], xmm3
	addps  xmm4, xmm6
	mulps  xmm4, xmm0
	mulps  xmm4, xmm0  	;# fscal
	
	movaps xmm0, [rsp + nb213_dxO]
	movaps xmm1, [rsp + nb213_dyO]
	movaps xmm2, [rsp + nb213_dzO]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4 ;# xmm0-xmm2 now contains tx-tz (partial force)
	
	movss  xmm3, [rsp + nb213_fixO]	
	movss  xmm4, [rsp + nb213_fiyO]	
	movss  xmm5, [rsp + nb213_fizO]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [rsp + nb213_fixO], xmm3	
	movss  [rsp + nb213_fiyO], xmm4	
	movss  [rsp + nb213_fizO], xmm5	;# updated the O force now do the H's
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2      
	shufps xmm3, xmm3, 0x39	;# shift right 
	shufps xmm4, xmm4, 0x39
	shufps xmm5, xmm5, 0x39
	addss  xmm3, [rsp + nb213_fixH1]
	addss  xmm4, [rsp + nb213_fiyH1]
	addss  xmm5, [rsp + nb213_fizH1]
	movss  [rsp + nb213_fixH1], xmm3	
	movss  [rsp + nb213_fiyH1], xmm4	
	movss  [rsp + nb213_fizH1], xmm5	;# updated the H1 force 

	shufps xmm3, xmm3, 0x39
	shufps xmm4, xmm4, 0x39
	shufps xmm5, xmm5, 0x39
	addss  xmm3, [rsp + nb213_fixH2]
	addss  xmm4, [rsp + nb213_fiyH2]
	addss  xmm5, [rsp + nb213_fizH2]
	movss  [rsp + nb213_fixH2], xmm3	
	movss  [rsp + nb213_fiyH2], xmm4	
	movss  [rsp + nb213_fizH2], xmm5	;# updated the H2 force 

	mov rdi, [rbp + nb213_faction]
	shufps xmm3, xmm3, 0x39
	shufps xmm4, xmm4, 0x39
	shufps xmm5, xmm5, 0x39
	addss  xmm3, [rsp + nb213_fixM]
	addss  xmm4, [rsp + nb213_fiyM]
	addss  xmm5, [rsp + nb213_fizM]
	movss  [rsp + nb213_fixM], xmm3	
	movss  [rsp + nb213_fiyM], xmm4	
	movss  [rsp + nb213_fizM], xmm5	;# updated the M force 

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

	dec dword ptr [rsp + nb213_innerk]
	jz    .nb213_updateouterdata
	jmp   .nb213_odd_loop
.nb213_updateouterdata:
	mov   ecx, [rsp + nb213_ii3]
	mov   rdi, [rbp + nb213_faction]
	mov   rsi, [rbp + nb213_fshift]
	mov   edx, [rsp + nb213_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb213_fixO]
	movaps xmm1, [rsp + nb213_fiyO]
	movaps xmm2, [rsp + nb213_fizO]

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
	movaps xmm0, [rsp + nb213_fixH1]
	movaps xmm1, [rsp + nb213_fiyH1]
	movaps xmm2, [rsp + nb213_fizH1]

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
	movaps xmm0, [rsp + nb213_fixH2]
	movaps xmm1, [rsp + nb213_fiyH2]
	movaps xmm2, [rsp + nb213_fizH2]

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
	movaps xmm0, [rsp + nb213_fixM]
	movaps xmm1, [rsp + nb213_fiyM]
	movaps xmm2, [rsp + nb213_fizM]

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
	mov esi, [rsp + nb213_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb213_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb213_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb213_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb213_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb213_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
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






	
.globl nb_kernel213nf_x86_64_sse
.globl _nb_kernel213nf_x86_64_sse
nb_kernel213nf_x86_64_sse:	
_nb_kernel213nf_x86_64_sse:	
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
	;# bottom of stack is cache-aligned for sse use 
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
.equiv          nb213nf_iqM,            192
.equiv          nb213nf_iqH,            208
.equiv          nb213nf_qqM,            224
.equiv          nb213nf_qqH,            240
.equiv          nb213nf_rinvH1,         256
.equiv          nb213nf_rinvH2,         272
.equiv          nb213nf_rinvM,          288
.equiv          nb213nf_two,            304
.equiv          nb213nf_c6,             320
.equiv          nb213nf_c12,            336
.equiv          nb213nf_krf,            352
.equiv          nb213nf_crf,            368
.equiv          nb213nf_krsqH1,         384
.equiv          nb213nf_krsqH2,         400
.equiv          nb213nf_krsqM,          416
.equiv          nb213nf_vctot,          432
.equiv          nb213nf_Vvdwtot,        448
.equiv          nb213nf_half,           464
.equiv          nb213nf_three,          480
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
	movss xmm0, [rsi]
	movss [rsp + nb213nf_facel], xmm0


	mov rsi, [rbp + nb213nf_argkrf]
	mov rdi, [rbp + nb213nf_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb213nf_krf], xmm1
	movaps [rsp + nb213nf_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb213nf_half], eax
	movss xmm1, [rsp + nb213nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb213nf_half],  xmm1
	movaps [rsp + nb213nf_two],  xmm2
	movaps [rsp + nb213nf_three],  xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb213nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb213nf_charge]
	movss xmm4, [rdx + rbx*4 + 4]	
	movss xmm3, [rdx + rbx*4 + 12]	
	mov rsi, [rbp + nb213nf_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb213nf_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb213nf_iqM], xmm3
	movaps [rsp + nb213nf_iqH], xmm4
	
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
        add   rbx, 1                           ;# rbx=nn1=nn0+10
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
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb213nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb213nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb213nf_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2	
	movaps xmm6, xmm0
	movaps xmm7, xmm1
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb213nf_pos]	;# rax = base of pos[]  
	mov   [rsp + nb213nf_ii3], ebx

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
	movaps [rsp + nb213nf_ixO], xmm3
	movaps [rsp + nb213nf_iyO], xmm4
	movaps [rsp + nb213nf_izO], xmm5
	movaps [rsp + nb213nf_ixH1], xmm6
	movaps [rsp + nb213nf_iyH1], xmm7

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
	movaps [rsp + nb213nf_izH1], xmm6
	movaps [rsp + nb213nf_ixH2], xmm0
	movaps [rsp + nb213nf_iyH2], xmm1
	movaps [rsp + nb213nf_izH2], xmm2
	movaps [rsp + nb213nf_ixM], xmm3
	movaps [rsp + nb213nf_iyM], xmm4
	movaps [rsp + nb213nf_izM], xmm5
	
	;# clear vctot 
	xorps xmm4, xmm4
	movaps [rsp + nb213nf_vctot], xmm4
	movaps [rsp + nb213nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb213nf_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb213nf_pos]
	mov   rdi, [rbp + nb213nf_faction]	
	mov   rax, [rsp + nb213nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb213nf_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb213nf_ninner]
	mov   [rsp + nb213nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb213nf_innerk], edx	;# number of innerloop atoms 


	jge   .nb213nf_unroll_loop
	jmp   .nb213nf_odd_inner
.nb213nf_unroll_loop:
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb213nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb213nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb213nf_charge]	;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	 	;# and in xmm4 
	mulps  xmm3, [rsp + nb213nf_iqM]
	mulps  xmm4, [rsp + nb213nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps  [rsp + nb213nf_qqM], xmm3
	movaps  [rsp + nb213nf_qqH], xmm4
	
	mov rsi, [rbp + nb213nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]
	mov rsi, [rbp + nb213nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb213nf_ntia]
	add rax, rdi
	add rbx, rdi
	add rcx, rdi
	add rdx, rdi

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

	movaps [rsp + nb213nf_c6], xmm4
	movaps [rsp + nb213nf_c12], xmm6

	mov rsi, [rbp + nb213nf_pos]   	;# base of pos[] 

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
	movaps xmm4, [rsp + nb213nf_ixO]
	movaps xmm5, [rsp + nb213nf_iyO]
	movaps xmm6, [rsp + nb213nf_izO]

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
	movaps xmm4, [rsp + nb213nf_ixH1]
	movaps xmm5, [rsp + nb213nf_iyH1]
	movaps xmm6, [rsp + nb213nf_izH1]

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
	movaps xmm3, [rsp + nb213nf_ixH2]
	movaps xmm4, [rsp + nb213nf_iyH2]
	movaps xmm5, [rsp + nb213nf_izH2]

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
	movaps xmm3, [rsp + nb213nf_iyM]
	movaps xmm4, [rsp + nb213nf_izM]
	subps  xmm3, xmm1
	subps  xmm4, xmm2
	movaps xmm2, [rsp + nb213nf_ixM]
	subps  xmm2, xmm0	

	;# square it 
	mulps xmm2,xmm2
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	addps xmm4, xmm3
	addps xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 
	movaps xmm0, xmm4
	movaps xmm1, xmm5
	movaps xmm2, xmm6
	mulps  xmm0, [rsp + nb213nf_krf]	
	mulps  xmm1, [rsp + nb213nf_krf]	
	mulps  xmm2, [rsp + nb213nf_krf]	
	movaps [rsp + nb213nf_krsqM], xmm0
	movaps [rsp + nb213nf_krsqH2], xmm1
	movaps [rsp + nb213nf_krsqH1], xmm2

	;# rsqH1 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm0, [rsp + nb213nf_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm0, [rsp + nb213nf_half]
	movaps  [rsp + nb213nf_rinvH1], xmm0	;# rinvH1 

	;# rsqH2 - seed to xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm0, [rsp + nb213nf_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm0, [rsp + nb213nf_half]
	movaps  [rsp + nb213nf_rinvH2], xmm0	;# rinvH2 

	;# rsqM - seed to xmm2 
	rsqrtps xmm2, xmm4
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm0, [rsp + nb213nf_three]
	mulps   xmm2, xmm4	;# rsq*lu*lu 
	subps   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm0, [rsp + nb213nf_half]
	movaps  [rsp + nb213nf_rinvM], xmm0
	
	;# Do the O LJ-only interaction directly.	
	rcpps   xmm2, xmm7
	movaps  xmm1, [rsp + nb213nf_two]
	mulps   xmm7, xmm2
	subps   xmm1, xmm7
	mulps   xmm2, xmm1 ;# rinvsq 
	movaps  xmm0, xmm2
	mulps   xmm0, xmm2  	;# r4
	mulps   xmm0, xmm2 	;# r6
	movaps  xmm1, xmm0
	mulps   xmm1, xmm1  	;# r12
	mulps   xmm0, [rsp + nb213nf_c6]
	mulps   xmm1, [rsp + nb213nf_c12]
	movaps  xmm3, xmm1
	subps   xmm3, xmm0  	;# Vvdw12-Vvdw6
	addps   xmm3, [rsp + nb213nf_Vvdwtot]
	movaps  [rsp + nb213nf_Vvdwtot], xmm3

	;# do H1 interactions
	movaps  xmm7, [rsp + nb213nf_rinvH1]
	addps xmm7, [rsp + nb213nf_krsqH1]
	subps xmm7, [rsp + nb213nf_crf] ;# xmm7=rinv+ krsq-crf 
	mulps xmm7, [rsp + nb213nf_qqH]
	addps xmm7, [rsp + nb213nf_vctot]

	;# H2 interactions 
	movaps  xmm6, [rsp + nb213nf_rinvH2]
	addps xmm6, [rsp + nb213nf_krsqH2]
	subps xmm6, [rsp + nb213nf_crf] ;# xmm6=rinv+ krsq-crf 
	mulps xmm6, [rsp + nb213nf_qqH]
	addps xmm6, xmm7 

	;# M interactions 
	movaps xmm5, [rsp + nb213nf_rinvM]
	addps xmm5, [rsp + nb213nf_krsqM]
	subps xmm5, [rsp + nb213nf_crf] ;# xmm5=rinv+ krsq-crf 
	mulps xmm5, [rsp + nb213nf_qqM]
	addps xmm5, xmm6
	movaps [rsp + nb213nf_vctot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb213nf_innerk],  4
	jl    .nb213nf_odd_inner
	jmp   .nb213nf_unroll_loop
.nb213nf_odd_inner:	
	add dword ptr [rsp + nb213nf_innerk],  4
	jnz   .nb213nf_odd_loop
	jmp   .nb213nf_updateouterdata
.nb213nf_odd_loop:

	mov   rdx, [rsp + nb213nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb213nf_innerjjnr],  4	

 	xorps xmm4, xmm4  	;# clear reg.
	movss xmm4, [rsp + nb213nf_iqM]
	mov rsi, [rbp + nb213nf_charge] 
	movhps xmm4, [rsp + nb213nf_iqH]  ;# [qM  0  qH  qH] 
	shufps xmm4, xmm4, 41	;# [0 qH qH qM]

	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb213nf_qqM], xmm3	;# use dummy qq for storage 
	
	xorps xmm6, xmm6
	mov rsi, [rbp + nb213nf_type]
	mov ebx, [rsi + rax*4]
	mov rsi, [rbp + nb213nf_vdwparam]
	shl ebx, 1	
	add ebx, [rsp + nb213nf_ntia]
	movlps xmm6, [rsi + rbx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 252  ;# 11111100
	shufps xmm7, xmm7, 253  ;# 11111101
	movaps [rsp + nb213nf_c6], xmm6
	movaps [rsp + nb213nf_c12], xmm7
	
	mov rsi, [rbp + nb213nf_pos]
	lea rax, [rax + rax*2]  

	movss xmm3, [rsp + nb213nf_ixO]
	movss xmm4, [rsp + nb213nf_iyO]
	movss xmm5, [rsp + nb213nf_izO]
	movss xmm0, [rsp + nb213nf_ixH1]
	movss xmm1, [rsp + nb213nf_iyH1]
	movss xmm2, [rsp + nb213nf_izH1]
	unpcklps xmm3, [rsp + nb213nf_ixH2] 	;# ixO ixH2 - -
	unpcklps xmm4, [rsp + nb213nf_iyH2]  	;# iyO iyH2 - -
	unpcklps xmm5, [rsp + nb213nf_izH2]	;# izO izH2 - -
	unpcklps xmm0, [rsp + nb213nf_ixM] 	;# ixH1 ixM - -
	unpcklps xmm1, [rsp + nb213nf_iyM]  	;# iyH1 iyM - -
	unpcklps xmm2, [rsp + nb213nf_izM]	;# izH1 izM - -
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
	movaps xmm0, xmm4
	mulps xmm0, [rsp + nb213nf_krf]
	movaps [rsp + nb213nf_krsqM], xmm0
	
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb213nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb213nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv

	movaps xmm1, xmm0	;# xmm1=rinv
	addps  xmm1, [rsp + nb213nf_krsqM]
	subps  xmm1, [rsp + nb213nf_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm1, [rsp + nb213nf_qqM]
	addps  xmm1, [rsp + nb213nf_vctot]	
	movaps [rsp + nb213nf_vctot], xmm1

	movaps xmm1, xmm0
	mulps  xmm1, xmm1
	movaps xmm2, xmm1
	mulss  xmm1, xmm1
	mulss  xmm1, xmm2	;# xmm1=rinvsix
	xorps  xmm4, xmm4
	movss  xmm4, xmm1
	mulss  xmm4, xmm4	;# xmm4=rinvtwelve 
	mulss  xmm1, [rsp + nb213nf_c6]
	mulss  xmm4, [rsp + nb213nf_c12]
	movaps xmm3, xmm4 
	subss  xmm3, xmm1	;# xmm3=Vvdw12-Vvdw6 
	addss  xmm3, [rsp + nb213nf_Vvdwtot]
	movss  [rsp + nb213nf_Vvdwtot], xmm3

	dec dword ptr [rsp + nb213nf_innerk]
	jz    .nb213nf_updateouterdata
	jmp   .nb213nf_odd_loop
.nb213nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb213nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb213nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb213nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb213nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb213nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb213nf_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
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
