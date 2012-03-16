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


.globl nb_kernel232_x86_64_sse
.globl _nb_kernel232_x86_64_sse
nb_kernel232_x86_64_sse:	
_nb_kernel232_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb232_fshift,           16
.equiv          nb232_gid,              24
.equiv          nb232_pos,              32
.equiv          nb232_faction,          40
.equiv          nb232_charge,           48
.equiv          nb232_p_facel,          56
.equiv          nb232_argkrf,           64
.equiv          nb232_argcrf,           72
.equiv          nb232_Vc,               80
.equiv          nb232_type,             88
.equiv          nb232_p_ntype,          96
.equiv          nb232_vdwparam,         104
.equiv          nb232_Vvdw,             112
.equiv          nb232_p_tabscale,       120
.equiv          nb232_VFtab,            128
.equiv          nb232_invsqrta,         136
.equiv          nb232_dvda,             144
.equiv          nb232_p_gbtabscale,     152
.equiv          nb232_GBtab,            160
.equiv          nb232_p_nthreads,       168
.equiv          nb232_count,            176
.equiv          nb232_mtx,              184
.equiv          nb232_outeriter,        192
.equiv          nb232_inneriter,        200
.equiv          nb232_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb232_ixO,              0
.equiv          nb232_iyO,              16
.equiv          nb232_izO,              32
.equiv          nb232_ixH1,             48
.equiv          nb232_iyH1,             64
.equiv          nb232_izH1,             80
.equiv          nb232_ixH2,             96
.equiv          nb232_iyH2,             112
.equiv          nb232_izH2,             128
.equiv          nb232_jxO,              144
.equiv          nb232_jyO,              160
.equiv          nb232_jzO,              176
.equiv          nb232_jxH1,             192
.equiv          nb232_jyH1,             208
.equiv          nb232_jzH1,             224
.equiv          nb232_jxH2,             240
.equiv          nb232_jyH2,             256
.equiv          nb232_jzH2,             272
.equiv          nb232_dxOO,             288
.equiv          nb232_dyOO,             304
.equiv          nb232_dzOO,             320
.equiv          nb232_dxOH1,            336
.equiv          nb232_dyOH1,            352
.equiv          nb232_dzOH1,            368
.equiv          nb232_dxOH2,            384
.equiv          nb232_dyOH2,            400
.equiv          nb232_dzOH2,            416
.equiv          nb232_dxH1O,            432
.equiv          nb232_dyH1O,            448
.equiv          nb232_dzH1O,            464
.equiv          nb232_dxH1H1,           480
.equiv          nb232_dyH1H1,           496
.equiv          nb232_dzH1H1,           512
.equiv          nb232_dxH1H2,           528
.equiv          nb232_dyH1H2,           544
.equiv          nb232_dzH1H2,           560
.equiv          nb232_dxH2O,            576
.equiv          nb232_dyH2O,            592
.equiv          nb232_dzH2O,            608
.equiv          nb232_dxH2H1,           624
.equiv          nb232_dyH2H1,           640
.equiv          nb232_dzH2H1,           656
.equiv          nb232_dxH2H2,           672
.equiv          nb232_dyH2H2,           688
.equiv          nb232_dzH2H2,           704
.equiv          nb232_qqOO,             720
.equiv          nb232_qqOH,             736
.equiv          nb232_qqHH,             752
.equiv          nb232_c6,               768
.equiv          nb232_c12,              784
.equiv          nb232_tsc,              800
.equiv          nb232_fstmp,            816
.equiv          nb232_vctot,            832
.equiv          nb232_Vvdwtot,          848
.equiv          nb232_fixO,             864
.equiv          nb232_fiyO,             880
.equiv          nb232_fizO,             896
.equiv          nb232_fixH1,            912
.equiv          nb232_fiyH1,            928
.equiv          nb232_fizH1,            944
.equiv          nb232_fixH2,            960
.equiv          nb232_fiyH2,            976
.equiv          nb232_fizH2,            992
.equiv          nb232_fjxO,             1008
.equiv          nb232_fjyO,             1024
.equiv          nb232_fjzO,             1040
.equiv          nb232_fjxH1,            1056
.equiv          nb232_fjyH1,            1072
.equiv          nb232_fjzH1,            1088
.equiv          nb232_fjxH2,            1104
.equiv          nb232_fjyH2,            1120
.equiv          nb232_fjzH2,            1136
.equiv          nb232_half,             1152
.equiv          nb232_three,            1168
.equiv          nb232_rsqOO,            1184
.equiv          nb232_rsqOH1,           1200
.equiv          nb232_rsqOH2,           1216
.equiv          nb232_rsqH1O,           1232
.equiv          nb232_rsqH1H1,          1248
.equiv          nb232_rsqH1H2,          1264
.equiv          nb232_rsqH2O,           1280
.equiv          nb232_rsqH2H1,          1296
.equiv          nb232_rsqH2H2,          1312
.equiv          nb232_rinvOO,           1328
.equiv          nb232_rinvOH1,          1344
.equiv          nb232_rinvOH2,          1360
.equiv          nb232_rinvH1O,          1376
.equiv          nb232_rinvH1H1,         1392
.equiv          nb232_rinvH1H2,         1408
.equiv          nb232_rinvH2O,          1424
.equiv          nb232_rinvH2H1,         1440
.equiv          nb232_rinvH2H2,         1456
.equiv          nb232_two,              1472
.equiv          nb232_krf,              1488
.equiv          nb232_crf,              1504
.equiv          nb232_nri,              1520
.equiv          nb232_iinr,             1528
.equiv          nb232_jindex,           1536
.equiv          nb232_jjnr,             1544
.equiv          nb232_shift,            1552
.equiv          nb232_shiftvec,         1560
.equiv          nb232_facel,            1568
.equiv          nb232_innerjjnr,        1576
.equiv          nb232_is3,              1584
.equiv          nb232_ii3,              1588
.equiv          nb232_innerk,           1592
.equiv          nb232_n,                1596
.equiv          nb232_nn1,              1600
.equiv          nb232_nouter,           1604
.equiv          nb232_ninner,           1608

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
	sub rsp, 1616		;# local variable stack space (n*16+8)
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
	mov [rsp + nb232_nouter], eax
	mov [rsp + nb232_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb232_nri], edi
	mov [rsp + nb232_iinr], rsi
	mov [rsp + nb232_jindex], rdx
	mov [rsp + nb232_jjnr], rcx
	mov [rsp + nb232_shift], r8
	mov [rsp + nb232_shiftvec], r9
	mov rsi, [rbp + nb232_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb232_facel], xmm0

	mov rax, [rbp + nb232_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb232_tsc], xmm3

	mov rsi, [rbp + nb232_argkrf]
	mov rdi, [rbp + nb232_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb232_krf], xmm1
	movaps [rsp + nb232_crf], xmm2
	
	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb232_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb232_charge]
	movss xmm3, [rdx + rbx*4]	
	movss xmm4, xmm3	
	movss xmm5, [rdx + rbx*4 + 4]	
	mov rsi, [rbp + nb232_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb232_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb232_qqOO], xmm3
	movaps [rsp + nb232_qqOH], xmm4
	movaps [rsp + nb232_qqHH], xmm5
	
	xorps xmm0, xmm0
	mov   rdx, [rbp + nb232_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb232_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	add   rdx, rcx
	mov   rax, [rbp + nb232_vdwparam]
	movlps xmm0, [rax + rdx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 85  ;# 01010101
	movaps [rsp + nb232_c6], xmm0
	movaps [rsp + nb232_c12], xmm1

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb232_half], eax
	movss xmm1, [rsp + nb232_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb232_half],  xmm1
	movaps [rsp + nb232_two],  xmm2
	movaps [rsp + nb232_three],  xmm3

.nb232_threadloop:
        mov   rsi, [rbp + nb232_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb232_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb232_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb232_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb232_n], eax
        mov [rsp + nb232_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb232_outerstart
        jmp .nb232_end

.nb232_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb232_nouter]
	mov [rsp + nb232_nouter], ebx

.nb232_outer:
	mov   rax, [rsp + nb232_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb232_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb232_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb232_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb232_pos]    ;# eax = base of pos[]  
	mov   [rsp + nb232_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [rax + rbx*4]
	addss xmm4, [rax + rbx*4 + 4]
	addss xmm5, [rax + rbx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb232_ixO], xmm3
	movaps [rsp + nb232_iyO], xmm4
	movaps [rsp + nb232_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [rax + rbx*4 + 12]
	addss xmm1, [rax + rbx*4 + 16]
	addss xmm2, [rax + rbx*4 + 20]		
	addss xmm3, [rax + rbx*4 + 24]
	addss xmm4, [rax + rbx*4 + 28]
	addss xmm5, [rax + rbx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb232_ixH1], xmm0
	movaps [rsp + nb232_iyH1], xmm1
	movaps [rsp + nb232_izH1], xmm2
	movaps [rsp + nb232_ixH2], xmm3
	movaps [rsp + nb232_iyH2], xmm4
	movaps [rsp + nb232_izH2], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb232_vctot], xmm4
	movaps [rsp + nb232_Vvdwtot], xmm4
	movaps [rsp + nb232_fixO], xmm4
	movaps [rsp + nb232_fiyO], xmm4
	movaps [rsp + nb232_fizO], xmm4
	movaps [rsp + nb232_fixH1], xmm4
	movaps [rsp + nb232_fiyH1], xmm4
	movaps [rsp + nb232_fizH1], xmm4
	movaps [rsp + nb232_fixH2], xmm4
	movaps [rsp + nb232_fiyH2], xmm4
	movaps [rsp + nb232_fizH2], xmm4
	
	mov   rax, [rsp + nb232_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb232_pos]
	mov   rdi, [rbp + nb232_faction]	
	mov   rax, [rsp + nb232_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb232_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb232_ninner]
	mov   [rsp + nb232_ninner], ecx
	add   edx, 0
	mov   [rsp + nb232_innerk], edx    ;# number of innerloop atoms 
	jge   .nb232_unroll_loop
	jmp   .nb232_single_check
.nb232_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb232_innerjjnr]     ;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb232_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb232_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2]     ;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	

	;# move j O coordinates to local temp variables 
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
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + nb232_ixO]
    subps xmm1, [rsp + nb232_iyO]
    subps xmm2, [rsp + nb232_izO]
    subps xmm3, [rsp + nb232_ixH1]
    subps xmm4, [rsp + nb232_iyH1]
    subps xmm5, [rsp + nb232_izH1]
    subps xmm6, [rsp + nb232_ixH2]
    subps xmm7, [rsp + nb232_iyH2]
    subps xmm8, [rsp + nb232_izH2]
    
	movaps [rsp + nb232_dxOO], xmm0
	movaps [rsp + nb232_dyOO], xmm1
	movaps [rsp + nb232_dzOO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb232_dxH1O], xmm3
	movaps [rsp + nb232_dyH1O], xmm4
	movaps [rsp + nb232_dzH1O], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb232_dxH2O], xmm6
	movaps [rsp + nb232_dyH2O], xmm7
	movaps [rsp + nb232_dzH2O], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt for jO atoms
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + nb232_three]
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

	movaps  xmm4, [rsp + nb232_half]
	mulps   xmm9, xmm4  ;# rinvOO 
	mulps   xmm10, xmm4 ;# rinvH1O
    mulps   xmm11, xmm4 ;# rinvH2O
	
	;# O interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11
    movaps [rsp + nb232_rsqOO], xmm0
    movaps [rsp + nb232_rsqOH1], xmm3
    movaps [rsp + nb232_rsqOH2], xmm6
    movaps [rsp + nb232_rinvOO], xmm9
    movaps [rsp + nb232_rinvOH1], xmm10
    movaps [rsp + nb232_rinvOH2], xmm11
    
    

    ;# table LJ interaction
    mulps  xmm0, xmm9
    mulps  xmm0, [rsp + nb232_tsc] ;# rtab

    ;# truncate and convert to integers
    cvttps2dq xmm1, xmm0

    ;# convert back to float
    cvtdq2ps  xmm2, xmm1
         
    ;# multiply by 8
    pslld   xmm1, 3

    ;# move to integer registers
    movhlps xmm13, xmm1
    movd    r8d, xmm1
    movd    r10d, xmm13
    shufps  xmm1, xmm1, 1
    shufps  xmm13, xmm13, 1
    movd    r9d, xmm1
    movd    r11d, xmm13
    
    ;# calculate eps
    subps     xmm0, xmm2
    mov  rsi, [rbp + nb232_VFtab]
            
    movlps xmm5, [rsi + r8*4]
   	movlps xmm9, [rsi + r8*4 + 16]

	movlps xmm7, [rsi + r10*4]
	movlps xmm11, [rsi + r10*4 + 16]

	movhps xmm5, [rsi + r9*4]
	movhps xmm9, [rsi + r9*4 + 16]

	movhps xmm7, [rsi + r11*4]
	movhps xmm11, [rsi + r11*4 + 16]

    movaps xmm4, xmm5
    movaps xmm8, xmm9
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm8, xmm11, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101
	shufps xmm9, xmm11, 221  ;# 11011101

	movlps xmm7, [rsi + r8*4 + 8]
	movlps xmm11, [rsi + r8*4 + 24]
    
	movlps xmm13, [rsi + r10*4 + 8]
	movlps xmm14, [rsi + r10*4 + 24]

	movhps xmm7, [rsi + r9*4 + 8]
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

    mulps  xmm7, xmm0    ;# Heps
    mulps  xmm11, xmm0 
    mulps  xmm6, xmm0   ;# Geps
    mulps  xmm10, xmm0 
    mulps  xmm7, xmm0   ;# Heps2
    mulps  xmm11, xmm0 
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
    mulps  xmm5, xmm0  ;# eps*Fp
    mulps  xmm9, xmm0
    movaps xmm12, [rsp + nb232_c6]
    movaps xmm13, [rsp + nb232_c12]
    addps  xmm5, xmm4 ;# VV
    addps  xmm9, xmm8

    mulps  xmm5, xmm12  ;# VV*c6 = vnb6
    mulps  xmm9, xmm13  ;# VV*c12 = vnb12
    addps  xmm5, xmm9
    addps  xmm5, [rsp + nb232_Vvdwtot]
    movaps [rsp + nb232_Vvdwtot], xmm5
        
    mulps  xmm7, xmm12   ;# FF*c6 = fnb6
    mulps  xmm11, xmm13   ;# FF*c12  = fnb12
    addps  xmm7, xmm11
    mulps  xmm7, [rsp + nb232_tsc]
    movaps [rsp + nb232_fstmp], xmm7 
           
    ;# Coulomb reaction-field interaction
    movaps xmm0, [rsp + nb232_rsqOO]
    movaps xmm3, [rsp + nb232_rsqOH1]
    movaps xmm6, [rsp + nb232_rsqOH2]
    movaps xmm9, [rsp + nb232_rinvOO]
    movaps xmm10, [rsp + nb232_rinvOH1]
    movaps xmm11, [rsp + nb232_rinvOH2]
    
    movaps xmm1, xmm9 ;# copy of rinv
    movaps xmm4, xmm10
    movaps xmm7, xmm11
    movaps xmm2, [rsp + nb232_krf]    
    mulps  xmm10, xmm10 ;# rinvsq
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
    subps  xmm2, [rsp + nb232_crf]   ;# rinv+krsq-crf
    subps  xmm5, [rsp + nb232_crf]
    subps  xmm8, [rsp + nb232_crf]   
    mulps  xmm2, [rsp + nb232_qqOO] ;# voul=qq*(rinv+ krsq-crf)
    mulps  xmm5, [rsp + nb232_qqOH] ;# voul=qq*(rinv+ krsq-crf)
    mulps  xmm8, [rsp + nb232_qqOH] ;# voul=qq*(rinv+ krsq-crf)
    addps  xmm0, xmm0 ;# 2*krsq
    addps  xmm3, xmm3 
    addps  xmm6, xmm6 
    subps  xmm1, xmm0 ;# rinv-2*krsq
    subps  xmm4, xmm3
    subps  xmm7, xmm6
    mulps  xmm1, [rsp + nb232_qqOO]   ;# (rinv-2*krsq)*qq
    mulps  xmm4, [rsp + nb232_qqOH] 
    mulps  xmm7, [rsp + nb232_qqOH] 
    addps  xmm2, [rsp + nb232_vctot]
    addps  xmm5, xmm8
    addps  xmm2, xmm5
    movaps xmm15, xmm2
    mulps  xmm1, xmm9   ;# fijC
    mulps  xmm4, xmm10
    mulps  xmm7, xmm11    
    ;# xmm1, xmm4, xmm7 contains fscal coul
    ;# xmm15 contains vctot
    subps xmm1, [rsp + nb232_fstmp]
    mulps xmm1, xmm9
    
	;# move j O forces to local temp variables 
	mov rdi, [rbp + nb232_faction]
    movlps xmm9, [rdi + rax*4] ;# jxOa jyOa  -   -
    movlps xmm10, [rdi + rcx*4] ;# jxOc jyOc  -   -
    movhps xmm9, [rdi + rbx*4] ;# jxOa jyOa jxOb jyOb 
    movhps xmm10, [rdi + rdx*4] ;# jxOc jyOc jxOd jyOd 

    movss  xmm11, [rdi + rax*4 + 8] ;# jzOa  -  -  -
    movss  xmm12, [rdi + rcx*4 + 8] ;# jzOc  -  -  -
    movhps xmm11, [rdi + rbx*4 + 8] ;# jzOa  -  jzOb  -
    movhps xmm12, [rdi + rdx*4 + 8] ;# jzOc  -  jzOd -
    
    shufps xmm11, xmm12,  136  ;# 10001000 => jzOa jzOb jzOc jzOd

    ;# xmm9: jxOa jyOa jxOb jyOb 
    ;# xmm10: jxOc jyOc jxOd jyOd
    ;# xmm11: jzOa jzOb jzOc jzOd


    movaps xmm0, xmm1
    movaps xmm2, xmm1
    movaps xmm3, xmm4
    movaps xmm5, xmm4
    movaps xmm6, xmm7
    movaps xmm8, xmm7

	mulps xmm0, [rsp + nb232_dxOO]
	mulps xmm1, [rsp + nb232_dyOO]
	mulps xmm2, [rsp + nb232_dzOO]
	mulps xmm3, [rsp + nb232_dxH1O]
	mulps xmm4, [rsp + nb232_dyH1O]
	mulps xmm5, [rsp + nb232_dzH1O]
	mulps xmm6, [rsp + nb232_dxH2O]
	mulps xmm7, [rsp + nb232_dyH2O]
	mulps xmm8, [rsp + nb232_dzH2O]

    movaps xmm13, xmm0
    movaps xmm14, xmm1
    addps xmm11, xmm2
    addps xmm0, [rsp + nb232_fixO]
    addps xmm1, [rsp + nb232_fiyO]
    addps xmm2, [rsp + nb232_fizO]

    addps xmm13,  xmm3
    addps xmm14, xmm4
    addps xmm11, xmm5
    addps xmm3, [rsp + nb232_fixH1]
    addps xmm4, [rsp + nb232_fiyH1]
    addps xmm5, [rsp + nb232_fizH1]

    addps xmm13,  xmm6
    addps xmm14, xmm7
    addps xmm11, xmm8
    addps xmm6, [rsp + nb232_fixH2]
    addps xmm7, [rsp + nb232_fiyH2]
    addps xmm8, [rsp + nb232_fizH2]

    movaps [rsp + nb232_fixO], xmm0
    movaps [rsp + nb232_fiyO], xmm1
    movaps [rsp + nb232_fizO], xmm2
    movaps [rsp + nb232_fixH1], xmm3
    movaps [rsp + nb232_fiyH1], xmm4
    movaps [rsp + nb232_fizH1], xmm5
    movaps [rsp + nb232_fixH2], xmm6
    movaps [rsp + nb232_fiyH2], xmm7
    movaps [rsp + nb232_fizH2], xmm8
    
    ;# xmm9 = fOx
    ;# xmm10 = fOy
    ;# xmm11 = fOz
    movaps xmm0, xmm13
    unpcklps xmm13, xmm14
    unpckhps xmm0, xmm14
    
    addps xmm9, xmm13
    addps xmm10, xmm0

    movhlps  xmm12, xmm11 ;# fOzc fOzd
    
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
    
    mov rsi, [ rbp + nb232_pos ]

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

    shufps   xmm2, xmm3,  136  ;# 10001000 => jzH1a jzH1b jzH1c jzH1d

    ;# xmm0 = H1x
    ;# xmm1 = H1y
    ;# xmm2 = H1z
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    
    subps xmm0, [rsp + nb232_ixO]
    subps xmm1, [rsp + nb232_iyO]
    subps xmm2, [rsp + nb232_izO]
    subps xmm3, [rsp + nb232_ixH1]
    subps xmm4, [rsp + nb232_iyH1]
    subps xmm5, [rsp + nb232_izH1]
    subps xmm6, [rsp + nb232_ixH2]
    subps xmm7, [rsp + nb232_iyH2]
    subps xmm8, [rsp + nb232_izH2]
    
	movaps [rsp + nb232_dxOH1], xmm0
	movaps [rsp + nb232_dyOH1], xmm1
	movaps [rsp + nb232_dzOH1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb232_dxH1H1], xmm3
	movaps [rsp + nb232_dyH1H1], xmm4
	movaps [rsp + nb232_dzH1H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb232_dxH2H1], xmm6
	movaps [rsp + nb232_dyH2H1], xmm7
	movaps [rsp + nb232_dzH2H1], xmm8
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
		
	movaps  xmm9, [rsp + nb232_three]
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

	movaps  xmm4, [rsp + nb232_half]
	mulps   xmm9, xmm4  ;# rinvOH1
	mulps   xmm10, xmm4 ;# rinvH1H1
    mulps   xmm11, xmm4 ;# rinvH2H1
	
	;# H1 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movaps xmm1, xmm9 ;# copy of rinv
    movaps xmm4, xmm10
    movaps xmm7, xmm11
    movaps xmm2, [rsp + nb232_krf]    
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
    movaps xmm14, [rsp + nb232_crf]
    subps  xmm2, xmm14   ;# rinv+krsq-crf
    subps  xmm5, xmm14
    subps  xmm8, xmm14
    movaps xmm12, [rsp + nb232_qqOH]
    movaps xmm13, [rsp + nb232_qqHH]    
    mulps  xmm2, xmm12 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    mulps  xmm5, xmm13 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    mulps  xmm8, xmm13 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    addps  xmm0, xmm0 ;# 2*krsq
    addps  xmm3, xmm3 
    addps  xmm6, xmm6 
    subps  xmm1, xmm0 ;# rinv-2*krsq
    subps  xmm4, xmm3
    subps  xmm7, xmm6
    mulps  xmm1, xmm12   ;# (rinv-2*krsq)*qq
    mulps  xmm4, xmm13
    mulps  xmm7, xmm13
    addps  xmm15, xmm2
    addps  xmm5, xmm8
    addps  xmm15, xmm5
    
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

	mulps xmm0, [rsp + nb232_dxOH1]
	mulps xmm1, [rsp + nb232_dyOH1]
	mulps xmm2, [rsp + nb232_dzOH1]
	mulps xmm3, [rsp + nb232_dxH1H1]
	mulps xmm4, [rsp + nb232_dyH1H1]
	mulps xmm5, [rsp + nb232_dzH1H1]
	mulps xmm6, [rsp + nb232_dxH2H1]
	mulps xmm7, [rsp + nb232_dyH2H1]
	mulps xmm8, [rsp + nb232_dzH2H1]

    movaps xmm13,  xmm0
    movaps xmm14, xmm1
    addps xmm11, xmm2
    addps xmm0, [rsp + nb232_fixO]
    addps xmm1, [rsp + nb232_fiyO]
    addps xmm2, [rsp + nb232_fizO]

    addps xmm13,  xmm3
    addps xmm14, xmm4
    addps xmm11, xmm5
    addps xmm3, [rsp + nb232_fixH1]
    addps xmm4, [rsp + nb232_fiyH1]
    addps xmm5, [rsp + nb232_fizH1]

    addps xmm13,  xmm6
    addps xmm14, xmm7
    addps xmm11, xmm8
    addps xmm6, [rsp + nb232_fixH2]
    addps xmm7, [rsp + nb232_fiyH2]
    addps xmm8, [rsp + nb232_fizH2]

    movaps [rsp + nb232_fixO], xmm0
    movaps [rsp + nb232_fiyO], xmm1
    movaps [rsp + nb232_fizO], xmm2
    movaps [rsp + nb232_fixH1], xmm3
    movaps [rsp + nb232_fiyH1], xmm4
    movaps [rsp + nb232_fizH1], xmm5
    movaps [rsp + nb232_fixH2], xmm6
    movaps [rsp + nb232_fiyH2], xmm7
    movaps [rsp + nb232_fizH2], xmm8
    
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
    movss  xmm5, [rsi + rbx*4 + 32] ;# jzH2b  -  -  -
    movss  xmm6, [rsi + rdx*4 + 32] ;# jzH2d  -  -  -
    movlhps xmm2, xmm5 ;# jzH2a  -  jzH2b  -
    movlhps xmm3, xmm6 ;# jzH2c  -  jzH2d -
    
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
    
    subps xmm0, [rsp + nb232_ixO]
    subps xmm1, [rsp + nb232_iyO]
    subps xmm2, [rsp + nb232_izO]
    subps xmm3, [rsp + nb232_ixH1]
    subps xmm4, [rsp + nb232_iyH1]
    subps xmm5, [rsp + nb232_izH1]
    subps xmm6, [rsp + nb232_ixH2]
    subps xmm7, [rsp + nb232_iyH2]
    subps xmm8, [rsp + nb232_izH2]
    
	movaps [rsp + nb232_dxOH2], xmm0
	movaps [rsp + nb232_dyOH2], xmm1
	movaps [rsp + nb232_dzOH2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb232_dxH1H2], xmm3
	movaps [rsp + nb232_dyH1H2], xmm4
	movaps [rsp + nb232_dzH1H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb232_dxH2H2], xmm6
	movaps [rsp + nb232_dyH2H2], xmm7
	movaps [rsp + nb232_dzH2H2], xmm8
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
		
	movaps  xmm9, [rsp + nb232_three]
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

	movaps  xmm4, [rsp + nb232_half]
	mulps   xmm9, xmm4  ;# rinvOH2
	mulps   xmm10, xmm4 ;# rinvH1H2
    mulps   xmm11, xmm4 ;# rinvH2H2
	
	;# H2 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movaps xmm1, xmm9 ;# copy of rinv
    movaps xmm4, xmm10
    movaps xmm7, xmm11
    movaps xmm2, [rsp + nb232_krf]    
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
    movaps xmm14, [rsp + nb232_crf]
    subps  xmm2, xmm14   ;# rinv+krsq-crf
    subps  xmm5, xmm14
    subps  xmm8, xmm14
    movaps xmm12, [rsp + nb232_qqOH]
    movaps xmm13, [rsp + nb232_qqHH]    
    mulps  xmm2, xmm12 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    mulps  xmm5, xmm13 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    mulps  xmm8, xmm13 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    addps  xmm0, xmm0 ;# 2*krsq
    addps  xmm3, xmm3 
    addps  xmm6, xmm6 
    subps  xmm1, xmm0 ;# rinv-2*krsq
    subps  xmm4, xmm3
    subps  xmm7, xmm6
    mulps  xmm1, xmm12   ;# (rinv-2*krsq)*qq
    mulps  xmm4, xmm13
    mulps  xmm7, xmm13
    addps  xmm5, xmm8
    addps  xmm2, xmm15
    addps  xmm2, xmm5
    movaps [rsp + nb232_vctot], xmm2
    
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
    movss  xmm2, [rdi + rbx*4 + 32] ;# jzH2b  -  -  -
    movss  xmm3, [rdi + rdx*4 + 32] ;# jzH2d  -  -  -
    movlhps xmm11, xmm2 ;# jzH2a  -  jzH2b  -
    movlhps xmm12, xmm3 ;# jzH2c  -  jzH2d -
    
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

	mulps xmm0, [rsp + nb232_dxOH2]
	mulps xmm1, [rsp + nb232_dyOH2]
	mulps xmm2, [rsp + nb232_dzOH2]
	mulps xmm3, [rsp + nb232_dxH1H2]
	mulps xmm4, [rsp + nb232_dyH1H2]
	mulps xmm5, [rsp + nb232_dzH1H2]
	mulps xmm6, [rsp + nb232_dxH2H2]
	mulps xmm7, [rsp + nb232_dyH2H2]
	mulps xmm8, [rsp + nb232_dzH2H2]

    movaps xmm13, xmm0
    movaps xmm14, xmm1
    addps xmm11, xmm2
    addps xmm0, [rsp + nb232_fixO]
    addps xmm1, [rsp + nb232_fiyO]
    addps xmm2, [rsp + nb232_fizO]

    addps xmm13,  xmm3
    addps xmm14, xmm4
    addps xmm11, xmm5
    addps xmm3, [rsp + nb232_fixH1]
    addps xmm4, [rsp + nb232_fiyH1]
    addps xmm5, [rsp + nb232_fizH1]

    addps xmm13,  xmm6
    addps xmm14, xmm7
    addps xmm11, xmm8
    addps xmm6, [rsp + nb232_fixH2]
    addps xmm7, [rsp + nb232_fiyH2]
    addps xmm8, [rsp + nb232_fizH2]

    movaps [rsp + nb232_fixO], xmm0
    movaps [rsp + nb232_fiyO], xmm1
    movaps [rsp + nb232_fizO], xmm2
    movaps [rsp + nb232_fixH1], xmm3
    movaps [rsp + nb232_fiyH1], xmm4
    movaps [rsp + nb232_fizH1], xmm5
    movaps [rsp + nb232_fixH2], xmm6
    movaps [rsp + nb232_fiyH2], xmm7
    movaps [rsp + nb232_fizH2], xmm8
    
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
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb232_innerk],  4
	jl    .nb232_single_check
	jmp   .nb232_unroll_loop
.nb232_single_check:
	add dword ptr [rsp + nb232_innerk],  4
	jnz   .nb232_single_loop
	jmp   .nb232_updateouterdata
.nb232_single_loop:
	mov   rdx, [rsp + nb232_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]
	add qword ptr [rsp + nb232_innerjjnr],  4

	mov rsi, [rbp + nb232_pos]
	lea   rax, [rax + rax*2]

	;# fetch j coordinates 
	xorps xmm0, xmm0
	xorps xmm1, xmm1
	xorps xmm2, xmm2
	
	movss xmm0, [rsi + rax*4]		;# jxO  -  -  -
	movss xmm1, [rsi + rax*4 + 4]		;# jyO  -  -  -
	movss xmm2, [rsi + rax*4 + 8]		;# jzO  -  -  -  

	movlps xmm6, [rsi + rax*4 + 12]		;# xmm6 = jxH1 jyH1   -    -
	movss  xmm7, [rsi + rax*4 + 20]		;# xmm7 = jzH1   -    -    - 
	movhps xmm6, [rsi + rax*4 + 24]		;# xmm6 = jxH1 jyH1 jxH2 jyH2
	movss  xmm5, [rsi + rax*4 + 32]		;# xmm5 = jzH2   -    -    -
	
	;# have all coords, time for some shuffling.

	shufps xmm6, xmm6, 216 ;# 11011000	;# xmm6 = jxH1 jxH2 jyH1 jyH2 
	unpcklps xmm7, xmm5			;# xmm7 = jzH1 jzH2   -    -

	movlhps xmm0, xmm6			;# xmm0 = jxO   0   jxH1 jxH2 
	shufps  xmm1, xmm6, 228 ;# 11100100	;# xmm1 = jyO   0   jyH1 jyH2 
	shufps  xmm2, xmm7, 68  ;# 01000100	;# xmm2 = jzO   0   jzH1 jzH2

	movaps [rsp + nb232_jxO], xmm0
	movaps [rsp + nb232_jyO], xmm1
	movaps [rsp + nb232_jzO], xmm2
	subps  xmm0, [rsp + nb232_ixO]
	subps  xmm1, [rsp + nb232_iyO]
	subps  xmm2, [rsp + nb232_izO]
	movaps [rsp + nb232_dxOO], xmm0
	movaps [rsp + nb232_dyOO], xmm1
	movaps [rsp + nb232_dzOO], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 
	movaps [rsp + nb232_rsqOO], xmm0
	
	movaps xmm6, xmm0
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	mulps   xmm6, [rsp + nb232_krf] ;# xmm6=krsq 
	movaps  xmm2, xmm1
	movaps  xmm7, xmm6
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb232_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [rsp + nb232_half] ;# rinv iO - j water in xmm3
	movaps  [rsp + nb232_rinvOO], xmm3
	
	addps   xmm6, xmm3	;# xmm6=rinv+ krsq 
	mulps   xmm7, [rsp + nb232_two]
	subps  xmm6, [rsp + nb232_crf]	;# xmm6=rinv+ krsq-crf 
	
	movaps  xmm0, xmm3  ;# rinv
	subps   xmm3, xmm7	;# xmm3=rinv-2*krsq 
	xorps   xmm4, xmm4
	;# fetch charges to xmm4 (temporary) 
	movss   xmm4, [rsp + nb232_qqOO]
	movhps  xmm4, [rsp + nb232_qqOH]

	mulps xmm6, xmm4	;# vcoul  
	mulps xmm3, xmm4	;# coul part of fs  
	mulps xmm3, xmm0
	movaps [rsp + nb232_fstmp], xmm3 ;# save it
	addps  xmm6, [rsp + nb232_vctot]
    movaps [rsp + nb232_vctot], xmm6

	movaps xmm0, [rsp + nb232_rinvOO]
	movss xmm1, xmm0
	mulss  xmm1, [rsp + nb232_rsqOO] ;# xmm1=r 
	mulss  xmm1, [rsp + nb232_tsc]
		
    cvttps2pi mm6, xmm1
    cvtpi2ps xmm3, mm6
	subss    xmm1, xmm3	;# xmm1=eps 
    movss xmm2, xmm1
    mulss  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 3
	
    mov  rsi, [rbp + nb232_VFtab]
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
    mulss  xmm7, [rsp + nb232_two]       ;# two*Heps2 
    addss  xmm7, xmm6
    addss  xmm7, xmm5 ;# xmm7=FF 
    mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
    addss  xmm5, xmm4 ;# xmm5=VV 

    movss xmm4, [rsp + nb232_c6]
    mulss  xmm7, xmm4    ;# fijD 
    mulss  xmm5, xmm4    ;# Vvdw6 
	movss  xmm3, [rsp + nb232_fstmp]
	mulps  xmm7, [rsp + nb232_tsc]
	subss  xmm3, xmm7
	movss  [rsp + nb232_fstmp], xmm3 

    addss  xmm5, [rsp + nb232_Vvdwtot]
    movss [rsp + nb232_Vvdwtot], xmm5

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
    mulss  xmm7, [rsp + nb232_two]       ;# two*Heps2 
    addss  xmm7, xmm6
    addss  xmm7, xmm5 ;# xmm7=FF 
    mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
    addss  xmm5, xmm4 ;# xmm5=VV 

    movss xmm4, [rsp + nb232_c12]
    mulss  xmm7, xmm4 ;# fijR 
    mulss  xmm5, xmm4 ;# Vvdw12 
	movaps xmm3, [rsp + nb232_fstmp]
	mulss  xmm7, [rsp + nb232_tsc]
	subss  xmm3, xmm7

    addss  xmm5, [rsp + nb232_Vvdwtot]
    movss [rsp + nb232_Vvdwtot], xmm5

	mulps  xmm0, xmm3
	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	
	mulps   xmm0, [rsp + nb232_dxOO]
	mulps   xmm1, [rsp + nb232_dyOO]
	mulps   xmm2, [rsp + nb232_dzOO]

	;# initial update for j forces 
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  [rsp + nb232_fjxO], xmm3
	movaps  [rsp + nb232_fjyO], xmm4
	movaps  [rsp + nb232_fjzO], xmm5
	addps   xmm0, [rsp + nb232_fixO]
	addps   xmm1, [rsp + nb232_fiyO]
	addps   xmm2, [rsp + nb232_fizO]
	movaps  [rsp + nb232_fixO], xmm0
	movaps  [rsp + nb232_fiyO], xmm1
	movaps  [rsp + nb232_fizO], xmm2

	
	;# done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
    movaps  xmm0, [rsp + nb232_jxO]
    movaps  xmm1, [rsp + nb232_jyO]
    movaps  xmm2, [rsp + nb232_jzO]
    movaps  xmm3, xmm0
    movaps  xmm4, xmm1
    movaps  xmm5, xmm2
    subps   xmm0, [rsp + nb232_ixH1]
    subps   xmm1, [rsp + nb232_iyH1]
    subps   xmm2, [rsp + nb232_izH1]	
    subps   xmm3, [rsp + nb232_ixH2]
    subps   xmm4, [rsp + nb232_iyH2]
    subps   xmm5, [rsp + nb232_izH2]	

	movaps [rsp + nb232_dxH1O], xmm0
	movaps [rsp + nb232_dyH1O], xmm1
	movaps [rsp + nb232_dzH1O], xmm2
	movaps [rsp + nb232_dxH2O], xmm3
	movaps [rsp + nb232_dyH2O], xmm4
	movaps [rsp + nb232_dzH2O], xmm5
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
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb232_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb232_half] ;# rinv H1 - j water 
	mulps   xmm7, [rsp + nb232_half] ;# rinv H2 - j water  

	mulps xmm0, [rsp + nb232_krf] ;# krsq 
	mulps xmm4, [rsp + nb232_krf] ;# krsq  


	;# assemble charges in xmm6 
	xorps   xmm6, xmm6
	movss   xmm6, [rsp + nb232_qqOH]
	movhps  xmm6, [rsp + nb232_qqHH]
	movaps  xmm1, xmm0
	movaps  xmm5, xmm4
	addps   xmm0, xmm3	;# krsq+ rinv 
	addps   xmm4, xmm7	;# krsq+ rinv 
	subps xmm0, [rsp + nb232_crf]
	subps xmm4, [rsp + nb232_crf]
	mulps   xmm1, [rsp + nb232_two]
	mulps   xmm5, [rsp + nb232_two]
	mulps   xmm0, xmm6	;# vcoul 
	mulps   xmm4, xmm6	;# vcoul 
	addps   xmm4, xmm0		
	addps   xmm4, [rsp + nb232_vctot]
	movaps  [rsp + nb232_vctot], xmm4
	movaps  xmm0, xmm3
	movaps  xmm4, xmm7
	mulps   xmm3, xmm3
	mulps   xmm7, xmm7
	subps   xmm0, xmm1
	subps   xmm4, xmm5
	mulps   xmm0, xmm6
	mulps   xmm4, xmm6
	mulps   xmm0, xmm3	;# fscal 
	mulps   xmm7, xmm4	;# fscal 
	
	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	mulps   xmm0, [rsp + nb232_dxH1O]
	mulps   xmm1, [rsp + nb232_dyH1O]
	mulps   xmm2, [rsp + nb232_dzH1O]
	;# update forces H1 - j water 
	movaps  xmm3, [rsp + nb232_fjxO]
	movaps  xmm4, [rsp + nb232_fjyO]
	movaps  xmm5, [rsp + nb232_fjzO]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  [rsp + nb232_fjxO], xmm3
	movaps  [rsp + nb232_fjyO], xmm4
	movaps  [rsp + nb232_fjzO], xmm5
	addps   xmm0, [rsp + nb232_fixH1]
	addps   xmm1, [rsp + nb232_fiyH1]
	addps   xmm2, [rsp + nb232_fizH1]
	movaps  [rsp + nb232_fixH1], xmm0
	movaps  [rsp + nb232_fiyH1], xmm1
	movaps  [rsp + nb232_fizH1], xmm2
	;# do forces H2 - j water 
	movaps xmm0, xmm7
	movaps xmm1, xmm7
	movaps xmm2, xmm7
	mulps   xmm0, [rsp + nb232_dxH2O]
	mulps   xmm1, [rsp + nb232_dyH2O]
	mulps   xmm2, [rsp + nb232_dzH2O]
	movaps  xmm3, [rsp + nb232_fjxO]
	movaps  xmm4, [rsp + nb232_fjyO]
	movaps  xmm5, [rsp + nb232_fjzO]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	mov     rsi, [rbp + nb232_faction]
	movaps  [rsp + nb232_fjxO], xmm3
	movaps  [rsp + nb232_fjyO], xmm4
	movaps  [rsp + nb232_fjzO], xmm5
	addps   xmm0, [rsp + nb232_fixH2]
	addps   xmm1, [rsp + nb232_fiyH2]
	addps   xmm2, [rsp + nb232_fizH2]
	movaps  [rsp + nb232_fixH2], xmm0
	movaps  [rsp + nb232_fiyH2], xmm1
	movaps  [rsp + nb232_fizH2], xmm2

	;# update j water forces from local variables 
	movlps  xmm0, [rsi + rax*4]
	movlps  xmm1, [rsi + rax*4 + 12]
	movhps  xmm1, [rsi + rax*4 + 24]
	movaps  xmm3, [rsp + nb232_fjxO]
	movaps  xmm4, [rsp + nb232_fjyO]
	movaps  xmm5, [rsp + nb232_fjzO]
	movaps  xmm6, xmm5
	movaps  xmm7, xmm5
	shufps  xmm6, xmm6, 2 ;# constant 00000010
	shufps  xmm7, xmm7, 3 ;# constant 00000011
	addss   xmm5, [rsi + rax*4 + 8]
	addss   xmm6, [rsi + rax*4 + 20]
	addss   xmm7, [rsi + rax*4 + 32]
	movss   [rsi + rax*4 + 8], xmm5
	movss   [rsi + rax*4 + 20], xmm6
	movss   [rsi + rax*4 + 32], xmm7
	movaps   xmm5, xmm3
	unpcklps xmm3, xmm4
	unpckhps xmm5, xmm4
	addps    xmm0, xmm3
	addps    xmm1, xmm5
	movlps  [rsi + rax*4], xmm0 
	movlps  [rsi + rax*4 + 12], xmm1 
	movhps  [rsi + rax*4 + 24], xmm1 
	
	dec dword ptr [rsp + nb232_innerk]
	jz    .nb232_updateouterdata
	jmp   .nb232_single_loop
.nb232_updateouterdata:
	mov   ecx, [rsp + nb232_ii3]
	mov   rdi, [rbp + nb232_faction]
	mov   rsi, [rbp + nb232_fshift]
	mov   edx, [rsp + nb232_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb232_fixO]
	movaps xmm1, [rsp + nb232_fiyO] 
	movaps xmm2, [rsp + nb232_fizO]

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
	movaps xmm0, [rsp + nb232_fixH1]
	movaps xmm1, [rsp + nb232_fiyH1]
	movaps xmm2, [rsp + nb232_fizH1]

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
	movaps xmm0, [rsp + nb232_fixH2]
	movaps xmm1, [rsp + nb232_fiyH2]
	movaps xmm2, [rsp + nb232_fizH2]

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

	;# increment fshift force  
	movlps  xmm3, [rsi + rdx*4]
	movss  xmm4, [rsi + rdx*4 + 8]
	subps  xmm3, xmm6
	subss  xmm4, xmm7
	movlps  [rsi + rdx*4],    xmm3
	movss  [rsi + rdx*4 + 8], xmm4

	;# get n from stack
	mov esi, [rsp + nb232_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb232_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb232_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb232_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb232_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb232_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb232_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb232_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb232_n], esi
        jmp .nb232_outer
.nb232_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb232_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb232_end
        ;# non-zero, do one more workunit
        jmp   .nb232_threadloop
.nb232_end:
	mov eax, [rsp + nb232_nouter]
	mov ebx, [rsp + nb232_ninner]
	mov rcx, [rbp + nb232_outeriter]
	mov rdx, [rbp + nb232_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 1616
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





.globl nb_kernel232nf_x86_64_sse
.globl _nb_kernel232nf_x86_64_sse
nb_kernel232nf_x86_64_sse:	
_nb_kernel232nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb232nf_fshift,         16
.equiv          nb232nf_gid,            24
.equiv          nb232nf_pos,            32
.equiv          nb232nf_faction,        40
.equiv          nb232nf_charge,         48
.equiv          nb232nf_p_facel,        56
.equiv          nb232nf_argkrf,         64
.equiv          nb232nf_argcrf,         72
.equiv          nb232nf_Vc,             80
.equiv          nb232nf_type,           88
.equiv          nb232nf_p_ntype,        96
.equiv          nb232nf_vdwparam,       104
.equiv          nb232nf_Vvdw,           112
.equiv          nb232nf_p_tabscale,     120
.equiv          nb232nf_VFtab,          128
.equiv          nb232nf_invsqrta,       136
.equiv          nb232nf_dvda,           144
.equiv          nb232nf_p_gbtabscale,   152
.equiv          nb232nf_GBtab,          160
.equiv          nb232nf_p_nthreads,     168
.equiv          nb232nf_count,          176
.equiv          nb232nf_mtx,            184
.equiv          nb232nf_outeriter,      192
.equiv          nb232nf_inneriter,      200
.equiv          nb232nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb232nf_ixO,            0
.equiv          nb232nf_iyO,            16
.equiv          nb232nf_izO,            32
.equiv          nb232nf_ixH1,           48
.equiv          nb232nf_iyH1,           64
.equiv          nb232nf_izH1,           80
.equiv          nb232nf_ixH2,           96
.equiv          nb232nf_iyH2,           112
.equiv          nb232nf_izH2,           128
.equiv          nb232nf_jxO,            144
.equiv          nb232nf_jyO,            160
.equiv          nb232nf_jzO,            176
.equiv          nb232nf_jxH1,           192
.equiv          nb232nf_jyH1,           208
.equiv          nb232nf_jzH1,           224
.equiv          nb232nf_jxH2,           240
.equiv          nb232nf_jyH2,           256
.equiv          nb232nf_jzH2,           272
.equiv          nb232nf_qqOO,           288
.equiv          nb232nf_qqOH,           304
.equiv          nb232nf_qqHH,           320
.equiv          nb232nf_c6,             336
.equiv          nb232nf_c12,            352
.equiv          nb232nf_vctot,          368
.equiv          nb232nf_Vvdwtot,        384
.equiv          nb232nf_half,           400
.equiv          nb232nf_three,          416
.equiv          nb232nf_rsqOO,          432
.equiv          nb232nf_rsqOH1,         448
.equiv          nb232nf_rsqOH2,         464
.equiv          nb232nf_rsqH1O,         480
.equiv          nb232nf_rsqH1H1,        496
.equiv          nb232nf_rsqH1H2,        512
.equiv          nb232nf_rsqH2O,         528
.equiv          nb232nf_rsqH2H1,        544
.equiv          nb232nf_rsqH2H2,        560
.equiv          nb232nf_rinvOO,         576
.equiv          nb232nf_rinvOH1,        592
.equiv          nb232nf_rinvOH2,        608
.equiv          nb232nf_rinvH1O,        624
.equiv          nb232nf_rinvH1H1,       640
.equiv          nb232nf_rinvH1H2,       656
.equiv          nb232nf_rinvH2O,        672
.equiv          nb232nf_rinvH2H1,       688
.equiv          nb232nf_rinvH2H2,       704
.equiv          nb232nf_krf,            720
.equiv          nb232nf_crf,            736
.equiv          nb232nf_tsc,            752
.equiv          nb232nf_nri,            768
.equiv          nb232nf_iinr,           776
.equiv          nb232nf_jindex,         784
.equiv          nb232nf_jjnr,           792
.equiv          nb232nf_shift,          800
.equiv          nb232nf_shiftvec,       808
.equiv          nb232nf_facel,          816
.equiv          nb232nf_innerjjnr,      824
.equiv          nb232nf_is3,            832
.equiv          nb232nf_ii3,            836
.equiv          nb232nf_innerk,         840
.equiv          nb232nf_n,              844
.equiv          nb232nf_nn1,            848
.equiv          nb232nf_nouter,         852
.equiv          nb232nf_ninner,         856
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
	sub rsp, 864		;# local variable stack space (n*16+8)
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
	mov [rsp + nb232nf_nouter], eax
	mov [rsp + nb232nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb232nf_nri], edi
	mov [rsp + nb232nf_iinr], rsi
	mov [rsp + nb232nf_jindex], rdx
	mov [rsp + nb232nf_jjnr], rcx
	mov [rsp + nb232nf_shift], r8
	mov [rsp + nb232nf_shiftvec], r9
	mov rsi, [rbp + nb232nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb232nf_facel], xmm0

	mov rax, [rbp + nb232nf_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb232nf_tsc], xmm3

	mov rsi, [rbp + nb232nf_argkrf]
	mov rdi, [rbp + nb232nf_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb232nf_krf], xmm1
	movaps [rsp + nb232nf_crf], xmm2

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb232nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb232nf_charge]
	movss xmm3, [rdx + rbx*4]	
	movss xmm4, xmm3	
	movss xmm5, [rdx + rbx*4 + 4]	
	mov rsi, [rbp + nb232nf_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb232nf_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb232nf_qqOO], xmm3
	movaps [rsp + nb232nf_qqOH], xmm4
	movaps [rsp + nb232nf_qqHH], xmm5
	
	xorps xmm0, xmm0
	mov   rdx, [rbp + nb232nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb232nf_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	add   rdx, rcx
	mov   rax, [rbp + nb232nf_vdwparam]
	movlps xmm0, [rax + rdx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 85  ;# 01010101
	movaps [rsp + nb232nf_c6], xmm0
	movaps [rsp + nb232nf_c12], xmm1

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb232nf_half], eax
	movss xmm1, [rsp + nb232nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb232nf_half],  xmm1
	movaps [rsp + nb232nf_three],  xmm3

.nb232nf_threadloop:
        mov   rsi, [rbp + nb232nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb232nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb232nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb232nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb232nf_n], eax
        mov [rsp + nb232nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb232nf_outerstart
        jmp .nb232nf_end

.nb232nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb232nf_nouter]
	mov [rsp + nb232nf_nouter], ebx

.nb232nf_outer:
	mov   rax, [rsp + nb232nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb232nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb232nf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb232nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb232nf_pos]    ;# eax = base of pos[]  
	mov   [rsp + nb232nf_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [rax + rbx*4]
	addss xmm4, [rax + rbx*4 + 4]
	addss xmm5, [rax + rbx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb232nf_ixO], xmm3
	movaps [rsp + nb232nf_iyO], xmm4
	movaps [rsp + nb232nf_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [rax + rbx*4 + 12]
	addss xmm1, [rax + rbx*4 + 16]
	addss xmm2, [rax + rbx*4 + 20]		
	addss xmm3, [rax + rbx*4 + 24]
	addss xmm4, [rax + rbx*4 + 28]
	addss xmm5, [rax + rbx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb232nf_ixH1], xmm0
	movaps [rsp + nb232nf_iyH1], xmm1
	movaps [rsp + nb232nf_izH1], xmm2
	movaps [rsp + nb232nf_ixH2], xmm3
	movaps [rsp + nb232nf_iyH2], xmm4
	movaps [rsp + nb232nf_izH2], xmm5

	;# clear vctot
	xorps xmm4, xmm4
	movaps [rsp + nb232nf_vctot], xmm4
	movaps [rsp + nb232nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb232nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb232nf_pos]
	mov   rax, [rsp + nb232nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb232nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb232nf_ninner]
	mov   [rsp + nb232nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb232nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb232nf_unroll_loop
	jmp   .nb232nf_single_check
.nb232nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb232nf_innerjjnr]     ;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb232nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb232nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2]     ;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	
	
	;# move j coordinates to local temp variables 
	movlps xmm2, [rsi + rax*4]
	movlps xmm3, [rsi + rax*4 + 12]
	movlps xmm4, [rsi + rax*4 + 24]

	movlps xmm5, [rsi + rbx*4]
	movlps xmm6, [rsi + rbx*4 + 12]
	movlps xmm7, [rsi + rbx*4 + 24]

	movhps xmm2, [rsi + rcx*4]
	movhps xmm3, [rsi + rcx*4 + 12]
	movhps xmm4, [rsi + rcx*4 + 24]

	movhps xmm5, [rsi + rdx*4]
	movhps xmm6, [rsi + rdx*4 + 12]
	movhps xmm7, [rsi + rdx*4 + 24]

	;# current state: 	
	;# xmm2= jxOa  jyOa  jxOc  jyOc 
	;# xmm3= jxH1a jyH1a jxH1c jyH1c 
	;# xmm4= jxH2a jyH2a jxH2c jyH2c 
	;# xmm5= jxOb  jyOb  jxOd  jyOd 
	;# xmm6= jxH1b jyH1b jxH1d jyH1d 
	;# xmm7= jxH2b jyH2b jxH2d jyH2d 
	
	movaps xmm0, xmm2
	movaps xmm1, xmm3
	unpcklps xmm0, xmm5	;# xmm0= jxOa  jxOb  jyOa  jyOb 
	unpcklps xmm1, xmm6	;# xmm1= jxH1a jxH1b jyH1a jyH1b 
	unpckhps xmm2, xmm5	;# xmm2= jxOc  jxOd  jyOc  jyOd 
	unpckhps xmm3, xmm6	;# xmm3= jxH1c jxH1d jyH1c jyH1d 
	movaps xmm5, xmm4
	movaps   xmm6, xmm0
	unpcklps xmm4, xmm7	;# xmm4= jxH2a jxH2b jyH2a jyH2b 		
	unpckhps xmm5, xmm7	;# xmm5= jxH2c jxH2d jyH2c jyH2d 
	movaps   xmm7, xmm1
	movlhps  xmm0, xmm2	;# xmm0= jxOa  jxOb  jxOc  jxOd 
	movaps [rsp + nb232nf_jxO], xmm0
	movhlps  xmm2, xmm6	;# xmm2= jyOa  jyOb  jyOc  jyOd 
	movaps [rsp + nb232nf_jyO], xmm2
	movlhps  xmm1, xmm3
	movaps [rsp + nb232nf_jxH1], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [rsp + nb232nf_jyH1], xmm3
	movlhps  xmm4, xmm5
	movaps [rsp + nb232nf_jxH2], xmm4
	movhlps  xmm5, xmm6
	movaps [rsp + nb232nf_jyH2], xmm5

	movss  xmm0, [rsi + rax*4 + 8]
	movss  xmm1, [rsi + rax*4 + 20]
	movss  xmm2, [rsi + rax*4 + 32]

	movss  xmm3, [rsi + rcx*4 + 8]
	movss  xmm4, [rsi + rcx*4 + 20]
	movss  xmm5, [rsi + rcx*4 + 32]

	movhps xmm0, [rsi + rbx*4 + 4]
	movhps xmm1, [rsi + rbx*4 + 16]
	movhps xmm2, [rsi + rbx*4 + 28]
	
	movhps xmm3, [rsi + rdx*4 + 4]
	movhps xmm4, [rsi + rdx*4 + 16]
	movhps xmm5, [rsi + rdx*4 + 28]
	
	shufps xmm0, xmm3, 204  ;# constant 11001100
	shufps xmm1, xmm4, 204  ;# constant 11001100
	shufps xmm2, xmm5, 204  ;# constant 11001100
	movaps [rsp + nb232nf_jzO],  xmm0
	movaps [rsp + nb232nf_jzH1],  xmm1
	movaps [rsp + nb232nf_jzH2],  xmm2

	movaps xmm0, [rsp + nb232nf_ixO]
	movaps xmm1, [rsp + nb232nf_iyO]
	movaps xmm2, [rsp + nb232nf_izO]
	movaps xmm3, [rsp + nb232nf_ixO]
	movaps xmm4, [rsp + nb232nf_iyO]
	movaps xmm5, [rsp + nb232nf_izO]
	subps  xmm0, [rsp + nb232nf_jxO]
	subps  xmm1, [rsp + nb232nf_jyO]
	subps  xmm2, [rsp + nb232nf_jzO]
	subps  xmm3, [rsp + nb232nf_jxH1]
	subps  xmm4, [rsp + nb232nf_jyH1]
	subps  xmm5, [rsp + nb232nf_jzH1]
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
	movaps [rsp + nb232nf_rsqOO], xmm0
	movaps [rsp + nb232nf_rsqOH1], xmm3

	movaps xmm0, [rsp + nb232nf_ixO]
	movaps xmm1, [rsp + nb232nf_iyO]
	movaps xmm2, [rsp + nb232nf_izO]
	movaps xmm3, [rsp + nb232nf_ixH1]
	movaps xmm4, [rsp + nb232nf_iyH1]
	movaps xmm5, [rsp + nb232nf_izH1]
	subps  xmm0, [rsp + nb232nf_jxH2]
	subps  xmm1, [rsp + nb232nf_jyH2]
	subps  xmm2, [rsp + nb232nf_jzH2]
	subps  xmm3, [rsp + nb232nf_jxO]
	subps  xmm4, [rsp + nb232nf_jyO]
	subps  xmm5, [rsp + nb232nf_jzO]
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
	movaps [rsp + nb232nf_rsqOH2], xmm0
	movaps [rsp + nb232nf_rsqH1O], xmm3

	movaps xmm0, [rsp + nb232nf_ixH1]
	movaps xmm1, [rsp + nb232nf_iyH1]
	movaps xmm2, [rsp + nb232nf_izH1]
	movaps xmm3, [rsp + nb232nf_ixH1]
	movaps xmm4, [rsp + nb232nf_iyH1]
	movaps xmm5, [rsp + nb232nf_izH1]
	subps  xmm0, [rsp + nb232nf_jxH1]
	subps  xmm1, [rsp + nb232nf_jyH1]
	subps  xmm2, [rsp + nb232nf_jzH1]
	subps  xmm3, [rsp + nb232nf_jxH2]
	subps  xmm4, [rsp + nb232nf_jyH2]
	subps  xmm5, [rsp + nb232nf_jzH2]
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
	movaps [rsp + nb232nf_rsqH1H1], xmm0
	movaps [rsp + nb232nf_rsqH1H2], xmm3

	movaps xmm0, [rsp + nb232nf_ixH2]
	movaps xmm1, [rsp + nb232nf_iyH2]
	movaps xmm2, [rsp + nb232nf_izH2]
	movaps xmm3, [rsp + nb232nf_ixH2]
	movaps xmm4, [rsp + nb232nf_iyH2]
	movaps xmm5, [rsp + nb232nf_izH2]
	subps  xmm0, [rsp + nb232nf_jxO]
	subps  xmm1, [rsp + nb232nf_jyO]
	subps  xmm2, [rsp + nb232nf_jzO]
	subps  xmm3, [rsp + nb232nf_jxH1]
	subps  xmm4, [rsp + nb232nf_jyH1]
	subps  xmm5, [rsp + nb232nf_jzH1]
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
	movaps [rsp + nb232nf_rsqH2O], xmm0
	movaps [rsp + nb232nf_rsqH2H1], xmm4

	movaps xmm0, [rsp + nb232nf_ixH2]
	movaps xmm1, [rsp + nb232nf_iyH2]
	movaps xmm2, [rsp + nb232nf_izH2]
	subps  xmm0, [rsp + nb232nf_jxH2]
	subps  xmm1, [rsp + nb232nf_jyH2]
	subps  xmm2, [rsp + nb232nf_jzH2]
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [rsp + nb232nf_rsqH2H2], xmm0
		
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb232nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb232nf_half] ;# rinvH2H2 
	mulps   xmm7, [rsp + nb232nf_half] ;# rinvH2H1 
	movaps  [rsp + nb232nf_rinvH2H2], xmm3
	movaps  [rsp + nb232nf_rinvH2H1], xmm7
	
	rsqrtps xmm1, [rsp + nb232nf_rsqOO]
	rsqrtps xmm5, [rsp + nb232nf_rsqOH1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb232nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb232nf_rsqOO]
	mulps   xmm5, [rsp + nb232nf_rsqOH1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb232nf_half] 
	mulps   xmm7, [rsp + nb232nf_half]
	movaps  [rsp + nb232nf_rinvOO], xmm3
	movaps  [rsp + nb232nf_rinvOH1], xmm7
	
	rsqrtps xmm1, [rsp + nb232nf_rsqOH2]
	rsqrtps xmm5, [rsp + nb232nf_rsqH1O]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb232nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb232nf_rsqOH2]
	mulps   xmm5, [rsp + nb232nf_rsqH1O]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb232nf_half] 
	mulps   xmm7, [rsp + nb232nf_half]
	movaps  [rsp + nb232nf_rinvOH2], xmm3
	movaps  [rsp + nb232nf_rinvH1O], xmm7
	
	rsqrtps xmm1, [rsp + nb232nf_rsqH1H1]
	rsqrtps xmm5, [rsp + nb232nf_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb232nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb232nf_rsqH1H1]
	mulps   xmm5, [rsp + nb232nf_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb232nf_half] 
	mulps   xmm7, [rsp + nb232nf_half]
	movaps  [rsp + nb232nf_rinvH1H1], xmm3
	movaps  [rsp + nb232nf_rinvH1H2], xmm7
	
	rsqrtps xmm1, [rsp + nb232nf_rsqH2O]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb232nf_three]
	mulps   xmm1, [rsp + nb232nf_rsqH2O]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + nb232nf_half] 
	movaps  [rsp + nb232nf_rinvH2O], xmm3

	;# start with OO interaction - first the table LJ part
	movaps xmm0, [rsp + nb232nf_rinvOO]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb232nf_rsqOO] ;# xmm1=r 
	mulps  xmm1, [rsp + nb232nf_tsc]
		
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

    mov  rsi, [rbp + nb232nf_VFtab]
	
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

    movaps xmm4, [rsp + nb232nf_c6]
    mulps  xmm5, xmm4    ;# Vvdw6 

    addps  xmm5, [rsp + nb232nf_Vvdwtot]
    movaps [rsp + nb232nf_Vvdwtot], xmm5

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
 
    movaps xmm4, [rsp + nb232nf_c12]
    mulps  xmm5, xmm4 ;# Vvdw12 

    addps  xmm5, [rsp + nb232nf_Vvdwtot]
    movaps [rsp + nb232nf_Vvdwtot], xmm5

	movaps xmm0, [rsp + nb232nf_rinvOO]
	movaps xmm2, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb232nf_krf]
	mulps  xmm5, [rsp + nb232nf_rsqOO] ;# xmm5=krsq 
	movaps xmm6, xmm5
	addps  xmm6, xmm2	;# xmm6=rinv+ krsq 
	subps  xmm6, [rsp + nb232nf_crf]
	
	mulps  xmm6, [rsp + nb232nf_qqOO] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addps  xmm6, [rsp + nb232nf_vctot] ;# local vctot summation variable 

	;# O-H1 interaction 
	movaps xmm0, [rsp + nb232nf_rinvOH1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb232nf_krf]
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb232nf_rsqOH1] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=rinv+ krsq 
	mulps  xmm0, xmm0
	subps  xmm4, [rsp + nb232nf_crf]
	mulps  xmm4, [rsp + nb232nf_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 
	addps  xmm6, xmm4	;# add to local vctot 

	;# O-H2 interaction  
	movaps xmm0, [rsp + nb232nf_rinvOH2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb232nf_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb232nf_rsqOH2] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	mulps xmm0, xmm0
	subps  xmm4, [rsp + nb232nf_crf]
	mulps  xmm4, [rsp + nb232nf_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 
	addps  xmm6, xmm4	;# add to local vctot 

	;# H1-O interaction 
	movaps xmm0, [rsp + nb232nf_rinvH1O]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb232nf_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb232nf_rsqH1O] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=rinv+ krsq 
	mulps xmm0, xmm0
	subps  xmm4, [rsp + nb232nf_crf]
	mulps  xmm4, [rsp + nb232nf_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 
	addps  xmm6, xmm4	;# add to local vctot 

	;# H1-H1 interaction 
	movaps xmm0, [rsp + nb232nf_rinvH1H1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb232nf_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb232nf_rsqH1H1] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb232nf_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb232nf_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	addps  xmm6, xmm4	;# add to local vctot 

	;# H1-H2 interaction 
	movaps xmm0, [rsp + nb232nf_rinvH1H2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb232nf_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb232nf_rsqH1H2] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	mulps xmm0, xmm0
	subps  xmm4, [rsp + nb232nf_crf]
	mulps  xmm4, [rsp + nb232nf_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	addps  xmm6, xmm4	;# add to local vctot 

	;# H2-O interaction 
	movaps xmm0, [rsp + nb232nf_rinvH2O]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb232nf_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb232nf_rsqH2O] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb232nf_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb232nf_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 
	addps  xmm6, xmm4	;# add to local vctot 

	;# H2-H1 interaction 
	movaps xmm0, [rsp + nb232nf_rinvH2H1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb232nf_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb232nf_rsqH2H1] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb232nf_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb232nf_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	addps  xmm6, xmm4	;# add to local vctot 

	;# H2-H2 interaction 
	movaps xmm0, [rsp + nb232nf_rinvH2H2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [rsp + nb232nf_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [rsp + nb232nf_rsqH2H2] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [rsp + nb232nf_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [rsp + nb232nf_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	addps  xmm6, xmm4	;# add to local vctot 
	movaps [rsp + nb232nf_vctot], xmm6
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb232nf_innerk],  4
	jl    .nb232nf_single_check
	jmp   .nb232nf_unroll_loop
.nb232nf_single_check:
	add dword ptr [rsp + nb232nf_innerk],  4
	jnz   .nb232nf_single_loop
	jmp   .nb232nf_updateouterdata
.nb232nf_single_loop:
	mov   rdx, [rsp + nb232nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb232nf_innerjjnr],  4

	mov rsi, [rbp + nb232nf_pos]
	lea   rax, [rax + rax*2]  

	;# fetch j coordinates 
	xorps xmm3, xmm3
	xorps xmm4, xmm4
	xorps xmm5, xmm5
	
	movss xmm3, [rsi + rax*4]		;# jxO  -  -  -
	movss xmm4, [rsi + rax*4 + 4]		;# jyO  -  -  -
	movss xmm5, [rsi + rax*4 + 8]		;# jzO  -  -  -  

	movlps xmm6, [rsi + rax*4 + 12]		;# xmm6 = jxH1 jyH1   -    -
	movss  xmm7, [rsi + rax*4 + 20]		;# xmm7 = jzH1   -    -    - 
	movhps xmm6, [rsi + rax*4 + 24]		;# xmm6 = jxH1 jyH1 jxH2 jyH2
	movss  xmm2, [rsi + rax*4 + 32]		;# xmm2 = jzH2   -    -    -
	
	;# have all coords, time for some shuffling.

	shufps xmm6, xmm6, 216 ;# constant 11011000	;# xmm6 = jxH1 jxH2 jyH1 jyH2 
	unpcklps xmm7, xmm2			;# xmm7 = jzH1 jzH2   -    -
	movaps  xmm0, [rsp + nb232nf_ixO]     
	movaps  xmm1, [rsp + nb232nf_iyO]
	movaps  xmm2, [rsp + nb232nf_izO]	
	movlhps xmm3, xmm6			;# xmm3 = jxO   0   jxH1 jxH2 
	shufps  xmm4, xmm6, 228 ;# constant 11100100	;# xmm4 = jyO   0   jyH1 jyH2 
	shufps  xmm5, xmm7, 68  ;# constant 01000100	;# xmm5 = jzO   0   jzH1 jzH2
	
	;# store all j coordinates in jO  
	movaps [rsp + nb232nf_jxO], xmm3
	movaps [rsp + nb232nf_jyO], xmm4
	movaps [rsp + nb232nf_jzO], xmm5
	subps  xmm0, xmm3
	subps  xmm1, xmm4
	subps  xmm2, xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 
	movaps [rsp + nb232nf_rsqOO], xmm0
	
	movaps xmm6, xmm0
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	mulps   xmm6, [rsp + nb232nf_krf] ;# xmm6=krsq 
	movaps  xmm2, xmm1
	movaps  xmm7, xmm6
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb232nf_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [rsp + nb232nf_half] ;# rinv iO - j water in xmm3
	movaps  [rsp + nb232nf_rinvOO], xmm3
	
	addps   xmm6, xmm3	;# xmm6=rinv+ krsq 
	subps  xmm6, [rsp + nb232nf_crf]	;# xmm6=rinv+ krsq-crf 
	
	movaps  xmm0, xmm3  ;# rinv
	xorps   xmm4, xmm4
	;# fetch charges to xmm4 (temporary) 
	movss   xmm4, [rsp + nb232nf_qqOO]
	movhps  xmm4, [rsp + nb232nf_qqOH]

	mulps xmm6, xmm4	;# vcoul  
	addps  xmm6, [rsp + nb232nf_vctot]
    movaps [rsp + nb232nf_vctot], xmm6

	movaps xmm0, [rsp + nb232nf_rinvOO]
	movss xmm1, xmm0
	mulss  xmm1, [rsp + nb232nf_rsqOO] ;# xmm1=r 
	mulss  xmm1, [rsp + nb232nf_tsc]
		
    cvttps2pi mm6, xmm1
    cvtpi2ps xmm3, mm6
	subss    xmm1, xmm3	;# xmm1=eps 
    movss xmm2, xmm1
    mulss  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 3
	
    mov  rsi, [rbp + nb232nf_VFtab]
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

    movss xmm4, [rsp + nb232nf_c6]
    mulss  xmm5, xmm4    ;# Vvdw6 

    addss  xmm5, [rsp + nb232nf_Vvdwtot]
    movss [rsp + nb232nf_Vvdwtot], xmm5

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

    movss xmm4, [rsp + nb232nf_c12]
    mulss  xmm5, xmm4 ;# Vvdw12 
    addss  xmm5, [rsp + nb232nf_Vvdwtot]
    movss [rsp + nb232nf_Vvdwtot], xmm5

	;# done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [rsp + nb232nf_ixH1]
	movaps  xmm1, [rsp + nb232nf_iyH1]
	movaps  xmm2, [rsp + nb232nf_izH1]	
	movaps  xmm3, [rsp + nb232nf_ixH2] 
	movaps  xmm4, [rsp + nb232nf_iyH2] 
	movaps  xmm5, [rsp + nb232nf_izH2] 
	subps   xmm0, [rsp + nb232nf_jxO]
	subps   xmm1, [rsp + nb232nf_jyO]
	subps   xmm2, [rsp + nb232nf_jzO]
	subps   xmm3, [rsp + nb232nf_jxO]
	subps   xmm4, [rsp + nb232nf_jyO]
	subps   xmm5, [rsp + nb232nf_jzO]
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
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb232nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb232nf_half] ;# rinv H1 - j water 
	mulps   xmm7, [rsp + nb232nf_half] ;# rinv H2 - j water  

	mulps xmm0, [rsp + nb232nf_krf] ;# krsq 
	mulps xmm4, [rsp + nb232nf_krf] ;# krsq  


	;# assemble charges in xmm6 
	xorps   xmm6, xmm6
	movss   xmm6, [rsp + nb232nf_qqOH]
	movhps  xmm6, [rsp + nb232nf_qqHH]
	movaps  xmm1, xmm0
	movaps  xmm5, xmm4
	addps   xmm0, xmm3	;# krsq+ rinv 
	addps   xmm4, xmm7	;# krsq+ rinv 
	subps xmm0, [rsp + nb232nf_crf]
	subps xmm4, [rsp + nb232nf_crf]
	mulps   xmm0, xmm6	;# vcoul 
	mulps   xmm4, xmm6	;# vcoul 
	addps   xmm4, xmm0		
	addps   xmm4, [rsp + nb232nf_vctot]
	movaps  [rsp + nb232nf_vctot], xmm4
	
	dec dword ptr [rsp + nb232nf_innerk]
	jz    .nb232nf_updateouterdata
	jmp   .nb232nf_single_loop
.nb232nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb232nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb232nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb232nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb232nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb232nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb232nf_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb232nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb232nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb232nf_n], esi
        jmp .nb232nf_outer
.nb232nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb232nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb232nf_end
        ;# non-zero, do one more workunit
        jmp   .nb232nf_threadloop
.nb232nf_end:
	mov eax, [rsp + nb232nf_nouter]
	mov ebx, [rsp + nb232nf_ninner]
	mov rcx, [rbp + nb232nf_outeriter]
	mov rdx, [rbp + nb232nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 864
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
