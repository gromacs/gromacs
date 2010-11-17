##
##
## Gromacs 4.0                         Copyright (c) 1991-2003 
## David van der Spoel, Erik Lindahl
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## To help us fund GROMACS development, we humbly ask that you cite
## the research papers on the package. Check out http://www.gromacs.org
## 
## And Hey:
## Gnomes, ROck Monsters And Chili Sauce
##





.globl nb_kernel303_x86_64_sse2
.globl _nb_kernel303_x86_64_sse2
nb_kernel303_x86_64_sse2:       
_nb_kernel303_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb303_fshift, 16
.set nb303_gid, 24
.set nb303_pos, 32
.set nb303_faction, 40
.set nb303_charge, 48
.set nb303_p_facel, 56
.set nb303_argkrf, 64
.set nb303_argcrf, 72
.set nb303_Vc, 80
.set nb303_type, 88
.set nb303_p_ntype, 96
.set nb303_vdwparam, 104
.set nb303_Vvdw, 112
.set nb303_p_tabscale, 120
.set nb303_VFtab, 128
.set nb303_invsqrta, 136
.set nb303_dvda, 144
.set nb303_p_gbtabscale, 152
.set nb303_GBtab, 160
.set nb303_p_nthreads, 168
.set nb303_count, 176
.set nb303_mtx, 184
.set nb303_outeriter, 192
.set nb303_inneriter, 200
.set nb303_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb303_ixM, 0
.set nb303_iyM, 16
.set nb303_izM, 32
.set nb303_ixH1, 48
.set nb303_iyH1, 64
.set nb303_izH1, 80
.set nb303_ixH2, 96
.set nb303_iyH2, 112
.set nb303_izH2, 128
.set nb303_iqM, 144
.set nb303_iqH, 160
.set nb303_dxM, 176
.set nb303_dyM, 192
.set nb303_dzM, 208
.set nb303_dxH1, 224
.set nb303_dyH1, 240
.set nb303_dzH1, 256
.set nb303_dxH2, 272
.set nb303_dyH2, 288
.set nb303_dzH2, 304
.set nb303_qqM, 320
.set nb303_qqH, 336
.set nb303_rinvM, 352
.set nb303_rinvH1, 368
.set nb303_rinvH2, 384
.set nb303_rM, 400
.set nb303_rH1, 416
.set nb303_rH2, 432
.set nb303_tsc, 448
.set nb303_two, 464
.set nb303_vctot, 480
.set nb303_fixM, 496
.set nb303_fiyM, 512
.set nb303_fizM, 528
.set nb303_fixH1, 544
.set nb303_fiyH1, 560
.set nb303_fizH1, 576
.set nb303_fixH2, 592
.set nb303_fiyH2, 608
.set nb303_fizH2, 624
.set nb303_fjx, 640
.set nb303_fjy, 656
.set nb303_fjz, 672
.set nb303_half, 688
.set nb303_three, 704
.set nb303_is3, 720
.set nb303_ii3, 724
.set nb303_nri, 728
.set nb303_iinr, 736
.set nb303_jindex, 744
.set nb303_jjnr, 752
.set nb303_shift, 760
.set nb303_shiftvec, 768
.set nb303_facel, 776
.set nb303_innerjjnr, 784
.set nb303_innerk, 792
.set nb303_n, 796
.set nb303_nn1, 800
.set nb303_nouter, 804
.set nb303_ninner, 808
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $824,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb303_nouter(%rsp)
        movl %eax,nb303_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb303_nri(%rsp)
        movq %rsi,nb303_iinr(%rsp)
        movq %rdx,nb303_jindex(%rsp)
        movq %rcx,nb303_jjnr(%rsp)
        movq %r8,nb303_shift(%rsp)
        movq %r9,nb303_shiftvec(%rsp)
        movq nb303_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb303_facel(%rsp)

        movq nb303_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb303_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb303_half(%rsp)
        movl %ebx,nb303_half+4(%rsp)
        movsd nb303_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb303_half(%rsp)
        movapd %xmm2,nb303_two(%rsp)
        movapd %xmm3,nb303_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb303_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb303_charge(%rbp),%rdx
        movsd 8(%rdx,%rbx,8),%xmm3
        movsd 24(%rdx,%rbx,8),%xmm4
        movq nb303_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb303_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb303_iqH(%rsp)
        movapd %xmm4,nb303_iqM(%rsp)

_nb_kernel303_x86_64_sse2.nb303_threadloop: 
        movq  nb303_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel303_x86_64_sse2.nb303_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel303_x86_64_sse2.nb303_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb303_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb303_n(%rsp)
        movl %ebx,nb303_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel303_x86_64_sse2.nb303_outerstart
        jmp _nb_kernel303_x86_64_sse2.nb303_end

_nb_kernel303_x86_64_sse2.nb303_outerstart: 
        ## ebx contains number of outer iterations
        addl nb303_nouter(%rsp),%ebx
        movl %ebx,nb303_nouter(%rsp)

_nb_kernel303_x86_64_sse2.nb303_outer: 
        movq  nb303_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb303_is3(%rsp)      ## store is3 

        movq  nb303_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb303_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb303_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb303_ii3(%rsp)

        addsd 24(%rax,%rbx,8),%xmm3
        addsd 32(%rax,%rbx,8),%xmm4
        addsd 40(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb303_ixH1(%rsp)
        movapd %xmm4,nb303_iyH1(%rsp)
        movapd %xmm5,nb303_izH1(%rsp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 48(%rax,%rbx,8),%xmm0
        addsd 56(%rax,%rbx,8),%xmm1
        addsd 64(%rax,%rbx,8),%xmm2
        addsd 72(%rax,%rbx,8),%xmm3
        addsd 80(%rax,%rbx,8),%xmm4
        addsd 88(%rax,%rbx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb303_ixH2(%rsp)
        movapd %xmm1,nb303_iyH2(%rsp)
        movapd %xmm2,nb303_izH2(%rsp)
        movapd %xmm3,nb303_ixM(%rsp)
        movapd %xmm4,nb303_iyM(%rsp)
        movapd %xmm5,nb303_izM(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb303_vctot(%rsp)
        movapd %xmm4,nb303_fixM(%rsp)
        movapd %xmm4,nb303_fiyM(%rsp)
        movapd %xmm4,nb303_fizM(%rsp)
        movapd %xmm4,nb303_fixH1(%rsp)
        movapd %xmm4,nb303_fiyH1(%rsp)
        movapd %xmm4,nb303_fizH1(%rsp)
        movapd %xmm4,nb303_fixH2(%rsp)
        movapd %xmm4,nb303_fiyH2(%rsp)
        movapd %xmm4,nb303_fizH2(%rsp)

        movq  nb303_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb303_pos(%rbp),%rsi
        movq  nb303_faction(%rbp),%rdi
        movq  nb303_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb303_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb303_ninner(%rsp),%ecx
        movl  %ecx,nb303_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb303_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel303_x86_64_sse2.nb303_unroll_loop
        jmp   _nb_kernel303_x86_64_sse2.nb303_checksingle
_nb_kernel303_x86_64_sse2.nb303_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb303_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb303_innerjjnr(%rsp)                   ## advance pointer (unrolled 2) 
        movq nb303_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb303_iqM(%rsp),%xmm3
        mulpd  nb303_iqH(%rsp),%xmm4

        movapd  %xmm3,nb303_qqM(%rsp)
        movapd  %xmm4,nb303_qqH(%rsp)

        movq nb303_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move j coordinates to local temp variables 
    movlpd (%rsi,%rax,8),%xmm0
    movlpd 8(%rsi,%rax,8),%xmm1
    movlpd 16(%rsi,%rax,8),%xmm2
    movhpd (%rsi,%rbx,8),%xmm0
    movhpd 8(%rsi,%rbx,8),%xmm1
    movhpd 16(%rsi,%rbx,8),%xmm2

    ## xmm0 = jx
    ## xmm1 = jy
    ## xmm2 = jz

    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb303_ixH1(%rsp),%xmm0
    subpd nb303_iyH1(%rsp),%xmm1
    subpd nb303_izH1(%rsp),%xmm2
    subpd nb303_ixH2(%rsp),%xmm3
    subpd nb303_iyH2(%rsp),%xmm4
    subpd nb303_izH2(%rsp),%xmm5
    subpd nb303_ixM(%rsp),%xmm6
    subpd nb303_iyM(%rsp),%xmm7
    subpd nb303_izM(%rsp),%xmm8

        movapd %xmm0,nb303_dxH1(%rsp)
        movapd %xmm1,nb303_dyH1(%rsp)
        movapd %xmm2,nb303_dzH1(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb303_dxH2(%rsp)
        movapd %xmm4,nb303_dyH2(%rsp)
        movapd %xmm5,nb303_dzH2(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb303_dxM(%rsp)
        movapd %xmm7,nb303_dyM(%rsp)
        movapd %xmm8,nb303_dzM(%rsp)
        mulpd  %xmm6,%xmm6
        mulpd  %xmm7,%xmm7
        mulpd  %xmm8,%xmm8
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
    addpd  %xmm7,%xmm6
    addpd  %xmm8,%xmm6

        ## start doing invsqrt for j atoms
    cvtpd2ps %xmm0,%xmm1
    cvtpd2ps %xmm3,%xmm4
    cvtpd2ps %xmm6,%xmm7
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm4,%xmm4
    rsqrtps %xmm7,%xmm7
    cvtps2pd %xmm1,%xmm1
    cvtps2pd %xmm4,%xmm4
    cvtps2pd %xmm7,%xmm7

        movapd  %xmm1,%xmm2
        movapd  %xmm4,%xmm5
    movapd  %xmm7,%xmm8

        mulpd   %xmm1,%xmm1 ## lu*lu
        mulpd   %xmm4,%xmm4 ## lu*lu
    mulpd   %xmm7,%xmm7 ## lu*lu

        movapd  nb303_three(%rsp),%xmm9
        movapd  %xmm9,%xmm10
    movapd  %xmm9,%xmm11

        mulpd   %xmm0,%xmm1 ## rsq*lu*lu
        mulpd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm7 ## rsq*lu*lu

        subpd   %xmm1,%xmm9
        subpd   %xmm4,%xmm10
    subpd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulpd   %xmm2,%xmm9
        mulpd   %xmm5,%xmm10
    mulpd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb303_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvH1
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH2
    mulpd   %xmm15,%xmm11 ## first iteration for rinvM

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb303_three(%rsp),%xmm1
        movapd  %xmm1,%xmm4
    movapd  %xmm1,%xmm7

        mulpd   %xmm0,%xmm2 ## rsq*lu*lu
        mulpd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm8 ## rsq*lu*lu

        subpd   %xmm2,%xmm1
        subpd   %xmm5,%xmm4
    subpd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulpd   %xmm1,%xmm9
        mulpd   %xmm4,%xmm10
    mulpd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb303_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvH1
        mulpd   %xmm15,%xmm10 ##   rinvH2
    mulpd   %xmm15,%xmm11 ##   rinvM

        movapd  %xmm9,nb303_rinvH1(%rsp)
        movapd  %xmm10,nb303_rinvH2(%rsp)
        movapd  %xmm11,nb303_rinvM(%rsp)

        ## interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd nb303_tsc(%rsp),%xmm1
    mulpd  %xmm9,%xmm0 ## r
    mulpd  %xmm10,%xmm3
    mulpd  %xmm11,%xmm6
    mulpd  %xmm1,%xmm0 ## rtab
    mulpd  %xmm1,%xmm3
    mulpd  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttpd2dq %xmm0,%xmm1
    cvttpd2dq %xmm3,%xmm4
    cvttpd2dq %xmm6,%xmm7

    ## convert back to float
    cvtdq2pd  %xmm1,%xmm2
    cvtdq2pd  %xmm4,%xmm5
    cvtdq2pd  %xmm7,%xmm8

    ## multiply by 4
    pslld   $2,%xmm1
    pslld   $2,%xmm4
    pslld   $2,%xmm7

    ## move to integer registers
    pshufd $1,%xmm1,%xmm13
    pshufd $1,%xmm4,%xmm14
    pshufd $1,%xmm7,%xmm15
    movd    %xmm1,%r8d
    movd    %xmm4,%r10d
    movd    %xmm7,%r12d
    movd    %xmm13,%r9d
    movd    %xmm14,%r11d
    movd    %xmm15,%r13d

    movq nb303_VFtab(%rbp),%rsi

    ## calculate eps
    subpd     %xmm2,%xmm0
    subpd     %xmm5,%xmm3
    subpd     %xmm8,%xmm6

    movapd    %xmm0,%xmm12 ## epsH1
    movapd    %xmm3,%xmm13 ## epsH2
    movapd    %xmm6,%xmm14 ## epsM

    ## Load LOTS of table data
    movlpd (%rsi,%r8,8),%xmm0
    movlpd 8(%rsi,%r8,8),%xmm1
    movlpd 16(%rsi,%r8,8),%xmm2
    movlpd 24(%rsi,%r8,8),%xmm3
    movlpd (%rsi,%r10,8),%xmm4
    movlpd 8(%rsi,%r10,8),%xmm5
    movlpd 16(%rsi,%r10,8),%xmm6
    movlpd 24(%rsi,%r10,8),%xmm7
    movlpd (%rsi,%r12,8),%xmm8
    movlpd 8(%rsi,%r12,8),%xmm9
    movlpd 16(%rsi,%r12,8),%xmm10
    movlpd 24(%rsi,%r12,8),%xmm11
    movhpd (%rsi,%r9,8),%xmm0
    movhpd 8(%rsi,%r9,8),%xmm1
    movhpd 16(%rsi,%r9,8),%xmm2
    movhpd 24(%rsi,%r9,8),%xmm3
    movhpd (%rsi,%r11,8),%xmm4
    movhpd 8(%rsi,%r11,8),%xmm5
    movhpd 16(%rsi,%r11,8),%xmm6
    movhpd 24(%rsi,%r11,8),%xmm7
    movhpd (%rsi,%r13,8),%xmm8
    movhpd 8(%rsi,%r13,8),%xmm9
    movhpd 16(%rsi,%r13,8),%xmm10
    movhpd 24(%rsi,%r13,8),%xmm11
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    mulpd  %xmm12,%xmm3  ## Heps
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11
    mulpd  %xmm12,%xmm2  ## Geps
    mulpd  %xmm13,%xmm6
    mulpd  %xmm14,%xmm10
    mulpd  %xmm12,%xmm3  ## Heps2
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11

    addpd  %xmm2,%xmm1  ## F+Geps
    addpd  %xmm6,%xmm5
    addpd  %xmm10,%xmm9
    addpd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addpd  %xmm7,%xmm5
    addpd  %xmm11,%xmm9
    addpd  %xmm3,%xmm3   ## 2*Heps2
    addpd  %xmm7,%xmm7
    addpd  %xmm11,%xmm11
    addpd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addpd  %xmm6,%xmm7
    addpd  %xmm10,%xmm11
    addpd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addpd  %xmm5,%xmm7
    addpd  %xmm9,%xmm11
    mulpd  %xmm12,%xmm1  ## eps*Fp
    mulpd  %xmm13,%xmm5
    mulpd  %xmm14,%xmm9
    movapd nb303_qqH(%rsp),%xmm12
    movapd nb303_qqM(%rsp),%xmm13
    addpd  %xmm0,%xmm1    ## VV
    addpd  %xmm4,%xmm5
    addpd  %xmm8,%xmm9
    mulpd  %xmm12,%xmm1  ## VV*qq = vcoul
    mulpd  %xmm12,%xmm5
    mulpd  %xmm13,%xmm9
    mulpd  %xmm12,%xmm3   ## FF*qq = fij
    mulpd  %xmm12,%xmm7
    mulpd  %xmm13,%xmm11

    ## accumulate vctot
    addpd  nb303_vctot(%rsp),%xmm1
    addpd  %xmm9,%xmm5
    addpd  %xmm5,%xmm1
    movapd %xmm1,nb303_vctot(%rsp)

    movapd nb303_tsc(%rsp),%xmm10
    mulpd  %xmm10,%xmm3 ## fscal
    mulpd  %xmm10,%xmm7
    mulpd  %xmm11,%xmm10

    xorpd %xmm4,%xmm4
    xorpd %xmm8,%xmm8
    xorpd %xmm11,%xmm11

    subpd %xmm3,%xmm4
    subpd %xmm7,%xmm8
    subpd %xmm10,%xmm11

    mulpd nb303_rinvH1(%rsp),%xmm4
    mulpd nb303_rinvH2(%rsp),%xmm8
    mulpd nb303_rinvM(%rsp),%xmm11

    ## move j forces to xmm0-xmm2
    movq nb303_faction(%rbp),%rdi
        movlpd (%rdi,%rax,8),%xmm0
        movlpd 8(%rdi,%rax,8),%xmm1
        movlpd 16(%rdi,%rax,8),%xmm2
        movhpd (%rdi,%rbx,8),%xmm0
        movhpd 8(%rdi,%rbx,8),%xmm1
        movhpd 16(%rdi,%rbx,8),%xmm2

    movapd %xmm4,%xmm3
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm7
    movapd %xmm8,%xmm9
    movapd %xmm11,%xmm10
    movapd %xmm11,%xmm12

        mulpd nb303_dxH1(%rsp),%xmm3
        mulpd nb303_dyH1(%rsp),%xmm4
        mulpd nb303_dzH1(%rsp),%xmm5
        mulpd nb303_dxH2(%rsp),%xmm7
        mulpd nb303_dyH2(%rsp),%xmm8
        mulpd nb303_dzH2(%rsp),%xmm9
        mulpd nb303_dxM(%rsp),%xmm10
        mulpd nb303_dyM(%rsp),%xmm11
        mulpd nb303_dzM(%rsp),%xmm12

    addpd %xmm3,%xmm0
    addpd %xmm4,%xmm1
    addpd %xmm5,%xmm2
    addpd nb303_fixH1(%rsp),%xmm3
    addpd nb303_fiyH1(%rsp),%xmm4
    addpd nb303_fizH1(%rsp),%xmm5

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb303_fixH2(%rsp),%xmm7
    addpd nb303_fiyH2(%rsp),%xmm8
    addpd nb303_fizH2(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb303_fixM(%rsp),%xmm10
    addpd nb303_fiyM(%rsp),%xmm11
    addpd nb303_fizM(%rsp),%xmm12

    movapd %xmm3,nb303_fixH1(%rsp)
    movapd %xmm4,nb303_fiyH1(%rsp)
    movapd %xmm5,nb303_fizH1(%rsp)
    movapd %xmm7,nb303_fixH2(%rsp)
    movapd %xmm8,nb303_fiyH2(%rsp)
    movapd %xmm9,nb303_fizH2(%rsp)
    movapd %xmm10,nb303_fixM(%rsp)
    movapd %xmm11,nb303_fiyM(%rsp)
    movapd %xmm12,nb303_fizM(%rsp)

    ## store back j forces from xmm0-xmm2
        movlpd %xmm0,(%rdi,%rax,8)
        movlpd %xmm1,8(%rdi,%rax,8)
        movlpd %xmm2,16(%rdi,%rax,8)
        movhpd %xmm0,(%rdi,%rbx,8)
        movhpd %xmm1,8(%rdi,%rbx,8)
        movhpd %xmm2,16(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb303_innerk(%rsp)
        jl    _nb_kernel303_x86_64_sse2.nb303_checksingle
        jmp   _nb_kernel303_x86_64_sse2.nb303_unroll_loop
_nb_kernel303_x86_64_sse2.nb303_checksingle: 
        movl  nb303_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel303_x86_64_sse2.nb303_dosingle
        jmp   _nb_kernel303_x86_64_sse2.nb303_updateouterdata
_nb_kernel303_x86_64_sse2.nb303_dosingle: 
        movq  nb303_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb303_charge(%rbp),%rsi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb303_iqM(%rsp),%xmm3
        mulpd  nb303_iqH(%rsp),%xmm4

        movapd  %xmm3,nb303_qqM(%rsp)
        movapd  %xmm4,nb303_qqH(%rsp)

        movq nb303_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        ## move coordinates to xmm4-xmm6 & xmm0-xmm2    
        movlpd (%rsi,%rax,8),%xmm4
        movlpd 8(%rsi,%rax,8),%xmm5
        movlpd 16(%rsi,%rax,8),%xmm6
    movapd %xmm4,%xmm0
    movapd %xmm5,%xmm1
    movapd %xmm6,%xmm2

        ## calc dr 
        subsd nb303_ixM(%rsp),%xmm4
        subsd nb303_iyM(%rsp),%xmm5
        subsd nb303_izM(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb303_dxM(%rsp)
        movapd %xmm5,nb303_dyM(%rsp)
        movapd %xmm6,nb303_dzM(%rsp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqM in xmm7 

        ## move j coords to xmm4-xmm6 
        movapd %xmm0,%xmm4
        movapd %xmm1,%xmm5
        movapd %xmm2,%xmm6

        ## calc dr 
        subsd nb303_ixH1(%rsp),%xmm4
        subsd nb303_iyH1(%rsp),%xmm5
        subsd nb303_izH1(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb303_dxH1(%rsp)
        movapd %xmm5,nb303_dyH1(%rsp)
        movapd %xmm6,nb303_dzH1(%rsp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move j coords to xmm3-xmm5 
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        ## calc dr 
        subsd nb303_ixH2(%rsp),%xmm3
        subsd nb303_iyH2(%rsp),%xmm4
        subsd nb303_izH2(%rsp),%xmm5

        ## store dr 
        movapd %xmm3,nb303_dxH2(%rsp)
        movapd %xmm4,nb303_dyH2(%rsp)
        movapd %xmm5,nb303_dzH2(%rsp)
        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqM in xmm7 

        ## start with rsqM - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb303_three(%rsp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb303_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm7,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb303_three(%rsp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb303_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,nb303_rinvM(%rsp)         ## rinvM in xmm4 
        mulsd   %xmm4,%xmm7
        movapd  %xmm7,nb303_rM(%rsp)    ## r in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb303_three(%rsp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb303_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb303_three(%rsp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb303_half(%rsp),%xmm4   ## rinv 
        movapd %xmm4,nb303_rinvH1(%rsp)         ## rinvH1 
        mulsd  %xmm4,%xmm6
        movapd %xmm6,nb303_rH1(%rsp)    ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb303_three(%rsp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb303_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb303_three(%rsp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb303_half(%rsp),%xmm4   ## rinv 
        movapd %xmm4,nb303_rinvH2(%rsp)   ## rinv 
        mulsd %xmm4,%xmm5
        movapd %xmm5,nb303_rH2(%rsp)   ## r 

        ## do M interactions 
        ## rM is still in xmm7 
        mulsd   nb303_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%r8d    ## mm6 = lu idx 
        cvtsi2sd %r8d,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%r8d            ## idx *= 4 
        movq nb303_VFtab(%rbp),%rsi

        movapd (%rsi,%r8,8),%xmm4       ## Y1 F1 
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1  

        movapd 16(%rsi,%r8,8),%xmm6     ## G1 H1 
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb303_two(%rsp),%xmm7    ## two*Heps2 
        movapd nb303_qqM(%rsp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 
    ## increment vcoul - then we can get rid of mm5 
    addsd  nb303_vctot(%rsp),%xmm5
    movlpd %xmm5,nb303_vctot(%rsp)
        xorpd  %xmm4,%xmm4

        mulsd  nb303_tsc(%rsp),%xmm3
        mulsd  nb303_rinvM(%rsp),%xmm3
        subsd  %xmm3,%xmm4

        movapd nb303_dxM(%rsp),%xmm0
        movapd nb303_dyM(%rsp),%xmm1
        movapd nb303_dzM(%rsp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2      ## tx in xmm0-xmm2 

        ## update M forces 
        movapd nb303_fixM(%rsp),%xmm3
        movapd nb303_fiyM(%rsp),%xmm4
        movapd nb303_fizM(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb303_fixM(%rsp)
        movlpd %xmm4,nb303_fiyM(%rsp)
        movlpd %xmm7,nb303_fizM(%rsp)
        ## update j forces with water M 
        movlpd %xmm0,nb303_fjx(%rsp)
        movlpd %xmm1,nb303_fjy(%rsp)
        movlpd %xmm2,nb303_fjz(%rsp)

        ## Done with M interactions - now H1! 
        movapd nb303_rH1(%rsp),%xmm7
        mulsd nb303_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%r8d    ## mm6 = lu idx 
        cvtsi2sd %r8d,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%r8d            ## idx *= 4 
        movq nb303_VFtab(%rbp),%rsi

        movapd (%rsi,%r8,8),%xmm4       ## Y1 F1 
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movapd 16(%rsi,%r8,8),%xmm6     ## G1 H1 
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb303_two(%rsp),%xmm7    ## two*Heps2 
        movapd nb303_qqH(%rsp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 
    ## increment vcoul 
        xorpd  %xmm4,%xmm4
    addsd  nb303_vctot(%rsp),%xmm5
        mulsd  nb303_rinvH1(%rsp),%xmm3
    movlpd %xmm5,nb303_vctot(%rsp)
        mulsd  nb303_tsc(%rsp),%xmm3
        subsd %xmm3,%xmm4

        movapd nb303_dxH1(%rsp),%xmm0
        movapd nb303_dyH1(%rsp),%xmm1
        movapd nb303_dzH1(%rsp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb303_fixH1(%rsp),%xmm3
        movapd nb303_fiyH1(%rsp),%xmm4
        movapd nb303_fizH1(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb303_fixH1(%rsp)
        movlpd %xmm4,nb303_fiyH1(%rsp)
        movlpd %xmm7,nb303_fizH1(%rsp)
        ## update j forces with water H1 
        addsd  nb303_fjx(%rsp),%xmm0
        addsd  nb303_fjy(%rsp),%xmm1
        addsd  nb303_fjz(%rsp),%xmm2
        movlpd %xmm0,nb303_fjx(%rsp)
        movlpd %xmm1,nb303_fjy(%rsp)
        movlpd %xmm2,nb303_fjz(%rsp)

        ## Done with H1, finally we do H2 interactions 
        movapd nb303_rH2(%rsp),%xmm7
        mulsd   nb303_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%r8d    ## mm6 = lu idx 
        cvtsi2sd %r8d,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%r8d            ## idx *= 4 
        movq nb303_VFtab(%rbp),%rsi

        movapd (%rsi,%r8,8),%xmm4       ## Y1 F1 
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%r8,8),%xmm6     ## G1 H1 
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb303_two(%rsp),%xmm7    ## two*Heps2 
        movapd nb303_qqH(%rsp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 
    ## increment vcoul 
        xorpd  %xmm4,%xmm4
    addsd  nb303_vctot(%rsp),%xmm5
        mulsd  nb303_rinvH2(%rsp),%xmm3
    movlpd %xmm5,nb303_vctot(%rsp)
        mulsd  nb303_tsc(%rsp),%xmm3
        subsd  %xmm3,%xmm4

        movapd nb303_dxH2(%rsp),%xmm0
        movapd nb303_dyH2(%rsp),%xmm1
        movapd nb303_dzH2(%rsp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb303_fixH2(%rsp),%xmm3
        movapd nb303_fiyH2(%rsp),%xmm4
        movapd nb303_fizH2(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb303_fixH2(%rsp)
        movlpd %xmm4,nb303_fiyH2(%rsp)
        movlpd %xmm7,nb303_fizH2(%rsp)

        movq nb303_faction(%rbp),%rdi
        ## update j forces 
        ## update j forces with water H1 
        addsd  nb303_fjx(%rsp),%xmm0
        addsd  nb303_fjy(%rsp),%xmm1
        addsd  nb303_fjz(%rsp),%xmm2

        ## the fj's - start by accumulating forces from memory 
        movlpd (%rdi,%rax,8),%xmm3
        movlpd 8(%rdi,%rax,8),%xmm4
        movlpd 16(%rdi,%rax,8),%xmm5
        addsd %xmm0,%xmm3
        addsd %xmm1,%xmm4
        addsd %xmm2,%xmm5
        movlpd %xmm3,(%rdi,%rax,8)
        movlpd %xmm4,8(%rdi,%rax,8)
        movlpd %xmm5,16(%rdi,%rax,8)

_nb_kernel303_x86_64_sse2.nb303_updateouterdata: 
        movl  nb303_ii3(%rsp),%ecx
        movq  nb303_faction(%rbp),%rdi
        movq  nb303_fshift(%rbp),%rsi
        movl  nb303_is3(%rsp),%edx

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movapd nb303_fixH1(%rsp),%xmm0
        movapd nb303_fiyH1(%rsp),%xmm1
        movapd nb303_fizH1(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  24(%rdi,%rcx,8),%xmm3
        movsd  32(%rdi,%rcx,8),%xmm4
        movsd  40(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,24(%rdi,%rcx,8)
        movsd  %xmm4,32(%rdi,%rcx,8)
        movsd  %xmm5,40(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        movapd %xmm0,%xmm6
        movsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movapd nb303_fixH2(%rsp),%xmm0
        movapd nb303_fiyH2(%rsp),%xmm1
        movapd nb303_fizH2(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  48(%rdi,%rcx,8),%xmm3
        movsd  56(%rdi,%rcx,8),%xmm4
        movsd  64(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,48(%rdi,%rcx,8)
        movsd  %xmm4,56(%rdi,%rcx,8)
        movsd  %xmm5,64(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movapd nb303_fixM(%rsp),%xmm0
        movapd nb303_fiyM(%rsp),%xmm1
        movapd nb303_fizM(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        ## increment i force 
        movsd  72(%rdi,%rcx,8),%xmm3
        movsd  80(%rdi,%rcx,8),%xmm4
        movsd  88(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,72(%rdi,%rcx,8)
        movsd  %xmm4,80(%rdi,%rcx,8)
        movsd  %xmm5,88(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## increment fshift force 
        movlpd (%rsi,%rdx,8),%xmm3
        movhpd 8(%rsi,%rdx,8),%xmm3
        movsd  16(%rsi,%rdx,8),%xmm4
        subpd  %xmm6,%xmm3
        subsd  %xmm7,%xmm4
        movlpd %xmm3,(%rsi,%rdx,8)
        movhpd %xmm3,8(%rsi,%rdx,8)
        movsd  %xmm4,16(%rsi,%rdx,8)

        ## get n from stack
        movl nb303_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb303_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb303_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb303_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb303_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel303_x86_64_sse2.nb303_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb303_n(%rsp)
        jmp _nb_kernel303_x86_64_sse2.nb303_outer
_nb_kernel303_x86_64_sse2.nb303_outerend: 
        ## check if more outer neighborlists remain
        movl  nb303_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel303_x86_64_sse2.nb303_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel303_x86_64_sse2.nb303_threadloop
_nb_kernel303_x86_64_sse2.nb303_end: 
        movl nb303_nouter(%rsp),%eax
        movl nb303_ninner(%rsp),%ebx
        movq nb303_outeriter(%rbp),%rcx
        movq nb303_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $824,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret




.globl nb_kernel303nf_x86_64_sse2
.globl _nb_kernel303nf_x86_64_sse2
nb_kernel303nf_x86_64_sse2:     
_nb_kernel303nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb303nf_fshift, 16
.set nb303nf_gid, 24
.set nb303nf_pos, 32
.set nb303nf_faction, 40
.set nb303nf_charge, 48
.set nb303nf_p_facel, 56
.set nb303nf_argkrf, 64
.set nb303nf_argcrf, 72
.set nb303nf_Vc, 80
.set nb303nf_type, 88
.set nb303nf_p_ntype, 96
.set nb303nf_vdwparam, 104
.set nb303nf_Vvdw, 112
.set nb303nf_p_tabscale, 120
.set nb303nf_VFtab, 128
.set nb303nf_invsqrta, 136
.set nb303nf_dvda, 144
.set nb303nf_p_gbtabscale, 152
.set nb303nf_GBtab, 160
.set nb303nf_p_nthreads, 168
.set nb303nf_count, 176
.set nb303nf_mtx, 184
.set nb303nf_outeriter, 192
.set nb303nf_inneriter, 200
.set nb303nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb303nf_ixM, 0
.set nb303nf_iyM, 16
.set nb303nf_izM, 32
.set nb303nf_ixH1, 48
.set nb303nf_iyH1, 64
.set nb303nf_izH1, 80
.set nb303nf_ixH2, 96
.set nb303nf_iyH2, 112
.set nb303nf_izH2, 128
.set nb303nf_iqM, 144
.set nb303nf_iqH, 160
.set nb303nf_qqM, 176
.set nb303nf_qqH, 192
.set nb303nf_rinvM, 208
.set nb303nf_rinvH1, 224
.set nb303nf_rinvH2, 240
.set nb303nf_rM, 256
.set nb303nf_rH1, 272
.set nb303nf_rH2, 288
.set nb303nf_tsc, 304
.set nb303nf_vctot, 320
.set nb303nf_half, 336
.set nb303nf_three, 352
.set nb303nf_is3, 368
.set nb303nf_ii3, 372
.set nb303nf_nri, 376
.set nb303nf_iinr, 384
.set nb303nf_jindex, 392
.set nb303nf_jjnr, 400
.set nb303nf_shift, 408
.set nb303nf_shiftvec, 416
.set nb303nf_facel, 424
.set nb303nf_innerjjnr, 432
.set nb303nf_innerk, 440
.set nb303nf_n, 444
.set nb303nf_nn1, 448
.set nb303nf_nouter, 452
.set nb303nf_ninner, 456
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $472,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb303nf_nouter(%rsp)
        movl %eax,nb303nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb303nf_nri(%rsp)
        movq %rsi,nb303nf_iinr(%rsp)
        movq %rdx,nb303nf_jindex(%rsp)
        movq %rcx,nb303nf_jjnr(%rsp)
        movq %r8,nb303nf_shift(%rsp)
        movq %r9,nb303nf_shiftvec(%rsp)
        movq nb303nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb303nf_facel(%rsp)

        movq nb303nf_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb303nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb303nf_half(%rsp)
        movl %ebx,nb303nf_half+4(%rsp)
        movsd nb303nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb303nf_half(%rsp)
        movapd %xmm3,nb303nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb303nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb303nf_charge(%rbp),%rdx
        movsd 8(%rdx,%rbx,8),%xmm3
        movsd 24(%rdx,%rbx,8),%xmm4
        movq nb303nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb303nf_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb303nf_iqH(%rsp)
        movapd %xmm4,nb303nf_iqM(%rsp)

_nb_kernel303nf_x86_64_sse2.nb303nf_threadloop: 
        movq  nb303nf_count(%rbp),%rsi          ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel303nf_x86_64_sse2.nb303nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel303nf_x86_64_sse2.nb303nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb303nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb303nf_n(%rsp)
        movl %ebx,nb303nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel303nf_x86_64_sse2.nb303nf_outerstart
        jmp _nb_kernel303nf_x86_64_sse2.nb303nf_end

_nb_kernel303nf_x86_64_sse2.nb303nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb303nf_nouter(%rsp),%ebx
        movl %ebx,nb303nf_nouter(%rsp)

_nb_kernel303nf_x86_64_sse2.nb303nf_outer: 
        movq  nb303nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 

        movq  nb303nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb303nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb303nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb303nf_ii3(%rsp)

        addsd 24(%rax,%rbx,8),%xmm3
        addsd 32(%rax,%rbx,8),%xmm4
        addsd 40(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb303nf_ixH1(%rsp)
        movapd %xmm4,nb303nf_iyH1(%rsp)
        movapd %xmm5,nb303nf_izH1(%rsp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 48(%rax,%rbx,8),%xmm0
        addsd 56(%rax,%rbx,8),%xmm1
        addsd 64(%rax,%rbx,8),%xmm2
        addsd 72(%rax,%rbx,8),%xmm3
        addsd 80(%rax,%rbx,8),%xmm4
        addsd 88(%rax,%rbx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb303nf_ixH2(%rsp)
        movapd %xmm1,nb303nf_iyH2(%rsp)
        movapd %xmm2,nb303nf_izH2(%rsp)
        movapd %xmm3,nb303nf_ixM(%rsp)
        movapd %xmm4,nb303nf_iyM(%rsp)
        movapd %xmm5,nb303nf_izM(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb303nf_vctot(%rsp)

        movq  nb303nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb303nf_pos(%rbp),%rsi
        movq  nb303nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb303nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb303nf_ninner(%rsp),%ecx
        movl  %ecx,nb303nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb303nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel303nf_x86_64_sse2.nb303nf_unroll_loop
        jmp   _nb_kernel303nf_x86_64_sse2.nb303nf_checksingle
_nb_kernel303nf_x86_64_sse2.nb303nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb303nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb303nf_innerjjnr(%rsp)                 ## advance pointer (unrolled 2) 
        movq nb303nf_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb303nf_iqM(%rsp),%xmm3
        mulpd  nb303nf_iqH(%rsp),%xmm4

        movapd  %xmm3,nb303nf_qqM(%rsp)
        movapd  %xmm4,nb303nf_qqH(%rsp)

        movq nb303nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move two coordinates to xmm0-xmm2    
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2
        movhpd (%rsi,%rbx,8),%xmm0
        movhpd 8(%rsi,%rbx,8),%xmm1
        movhpd 16(%rsi,%rbx,8),%xmm2

        ## move ixM-izM to xmm4-xmm6 
        movapd nb303nf_ixM(%rsp),%xmm4
        movapd nb303nf_iyM(%rsp),%xmm5
        movapd nb303nf_izM(%rsp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqM in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb303nf_ixH1(%rsp),%xmm4
        movapd nb303nf_iyH1(%rsp),%xmm5
        movapd nb303nf_izH1(%rsp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm6
        addpd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb303nf_ixH2(%rsp),%xmm3
        movapd nb303nf_iyH2(%rsp),%xmm4
        movapd nb303nf_izH2(%rsp),%xmm5

        ## calc dr 
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5

        ## square it 
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        addpd %xmm4,%xmm5
        addpd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqM in xmm7 

        ## start with rsqM - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb303nf_three(%rsp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb303nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm7,%xmm2
        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb303nf_three(%rsp),%xmm4
        subpd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb303nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,nb303nf_rinvM(%rsp)       ## rinvM in xmm4 
        mulpd   %xmm4,%xmm7
        movapd  %xmm7,nb303nf_rM(%rsp)          ## r in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb303nf_three(%rsp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb303nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb303nf_three(%rsp),%xmm4
        subpd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb303nf_half(%rsp),%xmm4   ## rinv 
        movapd %xmm4,nb303nf_rinvH1(%rsp)       ## rinvH1 
        mulpd  %xmm4,%xmm6
        movapd %xmm6,nb303nf_rH1(%rsp)          ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb303nf_three(%rsp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb303nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb303nf_three(%rsp),%xmm4
        subpd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb303nf_half(%rsp),%xmm4   ## rinv 
        movapd %xmm4,nb303nf_rinvH2(%rsp)   ## rinv 
        mulpd %xmm4,%xmm5
        movapd %xmm5,nb303nf_rH2(%rsp)   ## r 

        ## do M interactions 
        ## rM is still in xmm7 
        mulpd nb303nf_tsc(%rsp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb303nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb303nf_qqM(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul - then we can get rid of mm5 
    addpd  nb303nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb303nf_vctot(%rsp)

        ## Done with M interactions - now H1! 
        movapd nb303nf_rH1(%rsp),%xmm7
        mulpd nb303nf_tsc(%rsp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb303nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb303nf_qqH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul 
    addpd  nb303nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb303nf_vctot(%rsp)

        ## Done with H1, finally we do H2 interactions 
        movapd nb303nf_rH2(%rsp),%xmm7
        mulpd   nb303nf_tsc(%rsp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb303nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb303nf_qqH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul 
    addpd  nb303nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb303nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb303nf_innerk(%rsp)
        jl    _nb_kernel303nf_x86_64_sse2.nb303nf_checksingle
        jmp   _nb_kernel303nf_x86_64_sse2.nb303nf_unroll_loop
_nb_kernel303nf_x86_64_sse2.nb303nf_checksingle: 
        movl  nb303nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel303nf_x86_64_sse2.nb303nf_dosingle
        jmp   _nb_kernel303nf_x86_64_sse2.nb303nf_updateouterdata
_nb_kernel303nf_x86_64_sse2.nb303nf_dosingle: 
        movq  nb303nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb303nf_charge(%rbp),%rsi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb303nf_iqM(%rsp),%xmm3
        mulpd  nb303nf_iqH(%rsp),%xmm4

        movapd  %xmm3,nb303nf_qqM(%rsp)
        movapd  %xmm4,nb303nf_qqH(%rsp)

        movq nb303nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        ## move coordinates to xmm0-xmm2        
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move ixM-izM to xmm4-xmm6 
        movapd nb303nf_ixM(%rsp),%xmm4
        movapd nb303nf_iyM(%rsp),%xmm5
        movapd nb303nf_izM(%rsp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqM in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb303nf_ixH1(%rsp),%xmm4
        movapd nb303nf_iyH1(%rsp),%xmm5
        movapd nb303nf_izH1(%rsp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb303nf_ixH2(%rsp),%xmm3
        movapd nb303nf_iyH2(%rsp),%xmm4
        movapd nb303nf_izH2(%rsp),%xmm5

        ## calc dr 
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5

        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqM in xmm7 

        ## start with rsqM - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb303nf_three(%rsp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb303nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm7,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb303nf_three(%rsp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb303nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,nb303nf_rinvM(%rsp)       ## rinvM in xmm4 
        mulsd   %xmm4,%xmm7
        movapd  %xmm7,nb303nf_rM(%rsp)          ## r in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb303nf_three(%rsp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb303nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb303nf_three(%rsp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb303nf_half(%rsp),%xmm4   ## rinv 
        movapd %xmm4,nb303nf_rinvH1(%rsp)       ## rinvH1 
        mulsd  %xmm4,%xmm6
        movapd %xmm6,nb303nf_rH1(%rsp)          ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb303nf_three(%rsp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb303nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb303nf_three(%rsp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb303nf_half(%rsp),%xmm4   ## rinv 
        movapd %xmm4,nb303nf_rinvH2(%rsp)   ## rinv 
        mulsd %xmm4,%xmm5
        movapd %xmm5,nb303nf_rH2(%rsp)   ## r 

        ## do M interactions 
        movd %eax,%mm0
        ## rM is still in xmm7 
        mulsd   nb303nf_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb303nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1 
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1  

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1 
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb303nf_qqM(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul - then we can get rid of mm5 
    addsd  nb303nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb303nf_vctot(%rsp)

        ## Done with M interactions - now H1! 
        movapd nb303nf_rH1(%rsp),%xmm7
        mulsd nb303nf_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb303nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1 
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1 
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb303nf_qqH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul 
    addsd  nb303nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb303nf_vctot(%rsp)


        ## Done with H1, finally we do H2 interactions 
        movapd nb303nf_rH2(%rsp),%xmm7
        mulsd   nb303nf_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb303nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1 
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1 
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb303nf_qqH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul 
    addsd  nb303nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb303nf_vctot(%rsp)

_nb_kernel303nf_x86_64_sse2.nb303nf_updateouterdata: 
        ## get group index for i particle 
        ## get n from stack
        movl nb303nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb303nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb303nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb303nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb303nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel303nf_x86_64_sse2.nb303nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb303nf_n(%rsp)
        jmp _nb_kernel303nf_x86_64_sse2.nb303nf_outer
_nb_kernel303nf_x86_64_sse2.nb303nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb303nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel303nf_x86_64_sse2.nb303nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel303nf_x86_64_sse2.nb303nf_threadloop
_nb_kernel303nf_x86_64_sse2.nb303nf_end: 
        movl nb303nf_nouter(%rsp),%eax
        movl nb303nf_ninner(%rsp),%ebx
        movq nb303nf_outeriter(%rbp),%rcx
        movq nb303nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $472,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret


