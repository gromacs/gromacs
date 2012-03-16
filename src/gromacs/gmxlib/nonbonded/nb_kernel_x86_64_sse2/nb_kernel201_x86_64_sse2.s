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







.globl nb_kernel201_x86_64_sse2
.globl _nb_kernel201_x86_64_sse2
nb_kernel201_x86_64_sse2:       
_nb_kernel201_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb201_fshift, 16
.set nb201_gid, 24
.set nb201_pos, 32
.set nb201_faction, 40
.set nb201_charge, 48
.set nb201_p_facel, 56
.set nb201_argkrf, 64
.set nb201_argcrf, 72
.set nb201_Vc, 80
.set nb201_type, 88
.set nb201_p_ntype, 96
.set nb201_vdwparam, 104
.set nb201_Vvdw, 112
.set nb201_p_tabscale, 120
.set nb201_VFtab, 128
.set nb201_invsqrta, 136
.set nb201_dvda, 144
.set nb201_p_gbtabscale, 152
.set nb201_GBtab, 160
.set nb201_p_nthreads, 168
.set nb201_count, 176
.set nb201_mtx, 184
.set nb201_outeriter, 192
.set nb201_inneriter, 200
.set nb201_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb201_ixO, 0
.set nb201_iyO, 16
.set nb201_izO, 32
.set nb201_ixH1, 48
.set nb201_iyH1, 64
.set nb201_izH1, 80
.set nb201_ixH2, 96
.set nb201_iyH2, 112
.set nb201_izH2, 128
.set nb201_iqO, 144
.set nb201_iqH, 160
.set nb201_dxO, 176
.set nb201_dyO, 192
.set nb201_dzO, 208
.set nb201_dxH1, 224
.set nb201_dyH1, 240
.set nb201_dzH1, 256
.set nb201_dxH2, 272
.set nb201_dyH2, 288
.set nb201_dzH2, 304
.set nb201_qqO, 320
.set nb201_qqH, 336
.set nb201_vctot, 352
.set nb201_fixO, 384
.set nb201_fiyO, 400
.set nb201_fizO, 416
.set nb201_fixH1, 432
.set nb201_fiyH1, 448
.set nb201_fizH1, 464
.set nb201_fixH2, 480
.set nb201_fiyH2, 496
.set nb201_fizH2, 512
.set nb201_fjx, 528
.set nb201_fjy, 544
.set nb201_fjz, 560
.set nb201_half, 576
.set nb201_three, 592
.set nb201_two, 608
.set nb201_krf, 624
.set nb201_crf, 640
.set nb201_krsqO, 656
.set nb201_krsqH1, 672
.set nb201_krsqH2, 688
.set nb201_nri, 704
.set nb201_iinr, 712
.set nb201_jindex, 720
.set nb201_jjnr, 728
.set nb201_shift, 736
.set nb201_shiftvec, 744
.set nb201_facel, 752
.set nb201_innerjjnr, 760
.set nb201_is3, 768
.set nb201_ii3, 772
.set nb201_innerk, 776
.set nb201_n, 780
.set nb201_nn1, 784
.set nb201_nouter, 788
.set nb201_ninner, 792
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $808,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb201_nouter(%rsp)
        movl %eax,nb201_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb201_nri(%rsp)
        movq %rsi,nb201_iinr(%rsp)
        movq %rdx,nb201_jindex(%rsp)
        movq %rcx,nb201_jjnr(%rsp)
        movq %r8,nb201_shift(%rsp)
        movq %r9,nb201_shiftvec(%rsp)
        movq nb201_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb201_facel(%rsp)

        movq nb201_argkrf(%rbp),%rsi
        movq nb201_argcrf(%rbp),%rdi
        movsd (%rsi),%xmm1
        movsd (%rdi),%xmm2
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        movapd %xmm1,nb201_krf(%rsp)
        movapd %xmm2,nb201_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb201_half(%rsp)
        movl %ebx,nb201_half+4(%rsp)
        movsd nb201_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb201_half(%rsp)
        movapd %xmm2,nb201_two(%rsp)
        movapd %xmm3,nb201_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb201_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb201_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        movsd 8(%rdx,%rbx,8),%xmm4

        movsd nb201_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb201_iqO(%rsp)
        movapd %xmm4,nb201_iqH(%rsp)

_nb_kernel201_x86_64_sse2.nb201_threadloop: 
        movq  nb201_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel201_x86_64_sse2.nb201_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel201_x86_64_sse2.nb201_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb201_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb201_n(%rsp)
        movl %ebx,nb201_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel201_x86_64_sse2.nb201_outerstart
        jmp _nb_kernel201_x86_64_sse2.nb201_end

_nb_kernel201_x86_64_sse2.nb201_outerstart: 
        ## ebx contains number of outer iterations
        addl nb201_nouter(%rsp),%ebx
        movl %ebx,nb201_nouter(%rsp)

_nb_kernel201_x86_64_sse2.nb201_outer: 
        movq  nb201_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb201_is3(%rsp)      ## store is3 

        movq  nb201_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb201_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb201_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb201_ii3(%rsp)

        addsd (%rax,%rbx,8),%xmm3
        addsd 8(%rax,%rbx,8),%xmm4
        addsd 16(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb201_ixO(%rsp)
        movapd %xmm4,nb201_iyO(%rsp)
        movapd %xmm5,nb201_izO(%rsp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 24(%rax,%rbx,8),%xmm0
        addsd 32(%rax,%rbx,8),%xmm1
        addsd 40(%rax,%rbx,8),%xmm2
        addsd 48(%rax,%rbx,8),%xmm3
        addsd 56(%rax,%rbx,8),%xmm4
        addsd 64(%rax,%rbx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb201_ixH1(%rsp)
        movapd %xmm1,nb201_iyH1(%rsp)
        movapd %xmm2,nb201_izH1(%rsp)
        movapd %xmm3,nb201_ixH2(%rsp)
        movapd %xmm4,nb201_iyH2(%rsp)
        movapd %xmm5,nb201_izH2(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb201_vctot(%rsp)
        movapd %xmm4,nb201_fixO(%rsp)
        movapd %xmm4,nb201_fiyO(%rsp)
        movapd %xmm4,nb201_fizO(%rsp)
        movapd %xmm4,nb201_fixH1(%rsp)
        movapd %xmm4,nb201_fiyH1(%rsp)
        movapd %xmm4,nb201_fizH1(%rsp)
        movapd %xmm4,nb201_fixH2(%rsp)
        movapd %xmm4,nb201_fiyH2(%rsp)
        movapd %xmm4,nb201_fizH2(%rsp)

        movq  nb201_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb201_pos(%rbp),%rsi
        movq  nb201_faction(%rbp),%rdi
        movq  nb201_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb201_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb201_ninner(%rsp),%ecx
        movl  %ecx,nb201_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb201_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel201_x86_64_sse2.nb201_unroll_loop
        jmp   _nb_kernel201_x86_64_sse2.nb201_checksingle
_nb_kernel201_x86_64_sse2.nb201_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb201_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb201_innerjjnr(%rsp)                   ## advance pointer (unrolled 2) 

        movq nb201_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb201_iqO(%rsp),%xmm3
        mulpd  nb201_iqH(%rsp),%xmm4
        movapd  %xmm3,nb201_qqO(%rsp)
        movapd  %xmm4,nb201_qqH(%rsp)

        movq nb201_pos(%rbp),%rsi        ## base of pos[] 

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

    subpd nb201_ixO(%rsp),%xmm0
    subpd nb201_iyO(%rsp),%xmm1
    subpd nb201_izO(%rsp),%xmm2
    subpd nb201_ixH1(%rsp),%xmm3
    subpd nb201_iyH1(%rsp),%xmm4
    subpd nb201_izH1(%rsp),%xmm5
    subpd nb201_ixH2(%rsp),%xmm6
    subpd nb201_iyH2(%rsp),%xmm7
    subpd nb201_izH2(%rsp),%xmm8

        movapd %xmm0,nb201_dxO(%rsp)
        movapd %xmm1,nb201_dyO(%rsp)
        movapd %xmm2,nb201_dzO(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb201_dxH1(%rsp)
        movapd %xmm4,nb201_dyH1(%rsp)
        movapd %xmm5,nb201_dzH1(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb201_dxH2(%rsp)
        movapd %xmm7,nb201_dyH2(%rsp)
        movapd %xmm8,nb201_dzH2(%rsp)
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

        movapd  nb201_three(%rsp),%xmm9
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

        movapd  nb201_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvO
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH1
    mulpd   %xmm15,%xmm11 ## first iteration for rinvH2

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb201_three(%rsp),%xmm1
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

        movapd  nb201_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvO
        mulpd   %xmm15,%xmm10 ##   rinvH1
    mulpd   %xmm15,%xmm11 ##   rinvH2

        ## interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd %xmm9,%xmm1 ## copy of rinv
    movapd %xmm10,%xmm4
    movapd %xmm11,%xmm7
    movapd nb201_krf(%rsp),%xmm2
    mulpd  %xmm9,%xmm9  ## rinvsq
    mulpd  %xmm10,%xmm10
    mulpd  %xmm11,%xmm11
    mulpd  %xmm2,%xmm0 ## k*rsq
    mulpd  %xmm2,%xmm3
    mulpd  %xmm2,%xmm6
    movapd %xmm0,%xmm2 ## copy of k*rsq
    movapd %xmm3,%xmm5
    movapd %xmm6,%xmm8
    addpd  %xmm1,%xmm2 ## rinv+krsq
    addpd  %xmm4,%xmm5
    addpd  %xmm7,%xmm8
    movapd nb201_crf(%rsp),%xmm14
    subpd  %xmm14,%xmm2  ## rinv+krsq-crf
    subpd  %xmm14,%xmm5
    subpd  %xmm14,%xmm8
    movapd nb201_qqO(%rsp),%xmm12
    movapd nb201_qqH(%rsp),%xmm13
    mulpd  %xmm12,%xmm2 ## voul=qq*(rinv+ krsq-crf)
    mulpd  %xmm13,%xmm5 ## voul=qq*(rinv+ krsq-crf)
    mulpd  %xmm13,%xmm8 ## voul=qq*(rinv+ krsq-crf)
    addpd  %xmm0,%xmm0 ## 2*krsq
    addpd  %xmm3,%xmm3
    addpd  %xmm6,%xmm6
    subpd  %xmm0,%xmm1 ## rinv-2*krsq
    subpd  %xmm3,%xmm4
    subpd  %xmm6,%xmm7
    mulpd  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulpd  %xmm13,%xmm4
    mulpd  %xmm13,%xmm7
    addpd  nb201_vctot(%rsp),%xmm2
    addpd  %xmm8,%xmm5
    addpd  %xmm5,%xmm2
    movapd %xmm2,nb201_vctot(%rsp)

    mulpd  %xmm1,%xmm9  ## fscal
    mulpd  %xmm4,%xmm10
    mulpd  %xmm7,%xmm11

    ## move j forces to xmm0-xmm2
        movlpd (%rdi,%rax,8),%xmm0
        movlpd 8(%rdi,%rax,8),%xmm1
        movlpd 16(%rdi,%rax,8),%xmm2
        movhpd (%rdi,%rbx,8),%xmm0
        movhpd 8(%rdi,%rbx,8),%xmm1
        movhpd 16(%rdi,%rbx,8),%xmm2

    movapd %xmm9,%xmm7
    movapd %xmm9,%xmm8
    movapd %xmm11,%xmm13
    movapd %xmm11,%xmm14
    movapd %xmm11,%xmm15
    movapd %xmm10,%xmm11
    movapd %xmm10,%xmm12

        mulpd nb201_dxO(%rsp),%xmm7
        mulpd nb201_dyO(%rsp),%xmm8
        mulpd nb201_dzO(%rsp),%xmm9
        mulpd nb201_dxH1(%rsp),%xmm10
        mulpd nb201_dyH1(%rsp),%xmm11
        mulpd nb201_dzH1(%rsp),%xmm12
        mulpd nb201_dxH2(%rsp),%xmm13
        mulpd nb201_dyH2(%rsp),%xmm14
        mulpd nb201_dzH2(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb201_fixO(%rsp),%xmm7
    addpd nb201_fiyO(%rsp),%xmm8
    addpd nb201_fizO(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb201_fixH1(%rsp),%xmm10
    addpd nb201_fiyH1(%rsp),%xmm11
    addpd nb201_fizH1(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb201_fixH2(%rsp),%xmm13
    addpd nb201_fiyH2(%rsp),%xmm14
    addpd nb201_fizH2(%rsp),%xmm15

    movapd %xmm7,nb201_fixO(%rsp)
    movapd %xmm8,nb201_fiyO(%rsp)
    movapd %xmm9,nb201_fizO(%rsp)
    movapd %xmm10,nb201_fixH1(%rsp)
    movapd %xmm11,nb201_fiyH1(%rsp)
    movapd %xmm12,nb201_fizH1(%rsp)
    movapd %xmm13,nb201_fixH2(%rsp)
    movapd %xmm14,nb201_fiyH2(%rsp)
    movapd %xmm15,nb201_fizH2(%rsp)

    ## store back j O forces from xmm0-xmm2
        movlpd %xmm0,(%rdi,%rax,8)
        movlpd %xmm1,8(%rdi,%rax,8)
        movlpd %xmm2,16(%rdi,%rax,8)
        movhpd %xmm0,(%rdi,%rbx,8)
        movhpd %xmm1,8(%rdi,%rbx,8)
        movhpd %xmm2,16(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb201_innerk(%rsp)
        jl    _nb_kernel201_x86_64_sse2.nb201_checksingle
        jmp   _nb_kernel201_x86_64_sse2.nb201_unroll_loop
_nb_kernel201_x86_64_sse2.nb201_checksingle: 
        movl  nb201_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel201_x86_64_sse2.nb201_dosingle
        jmp   _nb_kernel201_x86_64_sse2.nb201_updateouterdata
_nb_kernel201_x86_64_sse2.nb201_dosingle: 
        movq  nb201_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb201_innerjjnr(%rsp)

        movq nb201_charge(%rbp),%rsi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb201_iqO(%rsp),%xmm3
        mulpd  nb201_iqH(%rsp),%xmm4
        movapd  %xmm3,nb201_qqO(%rsp)
        movapd  %xmm4,nb201_qqH(%rsp)

        movq nb201_pos(%rbp),%rsi        ## base of pos[] 
        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2        
        movlpd (%rsi,%rax,8),%xmm4
        movlpd 8(%rsi,%rax,8),%xmm5
        movlpd 16(%rsi,%rax,8),%xmm6
    movapd %xmm4,%xmm0
    movapd %xmm5,%xmm1
    movapd %xmm6,%xmm2

        ## calc dr 
        subsd nb201_ixO(%rsp),%xmm4
        subsd nb201_iyO(%rsp),%xmm5
        subsd nb201_izO(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb201_dxO(%rsp)
        movapd %xmm5,nb201_dyO(%rsp)
        movapd %xmm6,nb201_dzO(%rsp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move j coords to xmm4-xmm6 
        movapd %xmm0,%xmm4
        movapd %xmm1,%xmm5
        movapd %xmm2,%xmm6

        ## calc dr 
        subsd nb201_ixH1(%rsp),%xmm4
        subsd nb201_iyH1(%rsp),%xmm5
        subsd nb201_izH1(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb201_dxH1(%rsp)
        movapd %xmm5,nb201_dyH1(%rsp)
        movapd %xmm6,nb201_dzH1(%rsp)
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
        subsd nb201_ixH2(%rsp),%xmm3
        subsd nb201_iyH2(%rsp),%xmm4
        subsd nb201_izH2(%rsp),%xmm5

        ## store dr 
        movapd %xmm3,nb201_dxH2(%rsp)
        movapd %xmm4,nb201_dyH2(%rsp)
        movapd %xmm5,nb201_dzH2(%rsp)
        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        movapd %xmm5,%xmm0
        movapd %xmm6,%xmm1
        movapd %xmm7,%xmm2

        mulsd  nb201_krf(%rsp),%xmm0
        mulsd  nb201_krf(%rsp),%xmm1
        mulsd  nb201_krf(%rsp),%xmm2

        movapd %xmm0,nb201_krsqH2(%rsp)
        movapd %xmm1,nb201_krsqH1(%rsp)
        movapd %xmm2,nb201_krsqO(%rsp)

        ## start with rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb201_three(%rsp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb201_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb201_three(%rsp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb201_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb201_three(%rsp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb201_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb201_three(%rsp),%xmm4
        subsd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb201_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb201_three(%rsp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb201_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb201_three(%rsp),%xmm4
        subsd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb201_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movapd  %xmm7,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movapd  %xmm7,%xmm3
        movapd  nb201_krsqO(%rsp),%xmm0
        addsd   %xmm0,%xmm7     ## xmm6=rinv+ krsq 
        mulsd   nb201_two(%rsp),%xmm0
        subsd   nb201_crf(%rsp),%xmm7
        subsd   %xmm0,%xmm3     ## xmm7=rinv-2*krsq 
        mulsd   nb201_qqO(%rsp),%xmm7   ## vcoul 
        mulsd   nb201_qqO(%rsp),%xmm3
        mulsd  %xmm3,%xmm4      ## total fsH1 in xmm4 

        addsd  nb201_vctot(%rsp),%xmm7

        movapd nb201_dxO(%rsp),%xmm0
        movapd nb201_dyO(%rsp),%xmm1
        movapd nb201_dzO(%rsp),%xmm2
        movlpd %xmm7,nb201_vctot(%rsp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update O forces 
        movapd nb201_fixO(%rsp),%xmm3
        movapd nb201_fiyO(%rsp),%xmm4
        movapd nb201_fizO(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb201_fixO(%rsp)
        movlpd %xmm4,nb201_fiyO(%rsp)
        movlpd %xmm7,nb201_fizO(%rsp)
        ## update j forces with water O 
        movlpd %xmm0,nb201_fjx(%rsp)
        movlpd %xmm1,nb201_fjy(%rsp)
        movlpd %xmm2,nb201_fjz(%rsp)

        ## H1 interactions 
        movapd  %xmm6,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movapd  %xmm6,%xmm7
        movapd  nb201_krsqH1(%rsp),%xmm0
        addsd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        mulsd   nb201_two(%rsp),%xmm0
        subsd   nb201_crf(%rsp),%xmm6
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb201_qqH(%rsp),%xmm6   ## vcoul 
        mulsd   nb201_qqH(%rsp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH1 in xmm4 

        addsd  nb201_vctot(%rsp),%xmm6

        movapd nb201_dxH1(%rsp),%xmm0
        movapd nb201_dyH1(%rsp),%xmm1
        movapd nb201_dzH1(%rsp),%xmm2
        movlpd %xmm6,nb201_vctot(%rsp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb201_fixH1(%rsp),%xmm3
        movapd nb201_fiyH1(%rsp),%xmm4
        movapd nb201_fizH1(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb201_fixH1(%rsp)
        movlpd %xmm4,nb201_fiyH1(%rsp)
        movlpd %xmm7,nb201_fizH1(%rsp)
        ## update j forces with water H1 
        addsd  nb201_fjx(%rsp),%xmm0
        addsd  nb201_fjy(%rsp),%xmm1
        addsd  nb201_fjz(%rsp),%xmm2
        movlpd %xmm0,nb201_fjx(%rsp)
        movlpd %xmm1,nb201_fjy(%rsp)
        movlpd %xmm2,nb201_fjz(%rsp)

        ## H2 interactions 
        movapd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movapd  %xmm5,%xmm7
        movapd  nb201_krsqH2(%rsp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        mulsd   nb201_two(%rsp),%xmm0
        subsd   nb201_crf(%rsp),%xmm5
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb201_qqH(%rsp),%xmm5   ## vcoul 
        mulsd   nb201_qqH(%rsp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addsd  nb201_vctot(%rsp),%xmm5

        movapd nb201_dxH2(%rsp),%xmm0
        movapd nb201_dyH2(%rsp),%xmm1
        movapd nb201_dzH2(%rsp),%xmm2
        movlpd %xmm5,nb201_vctot(%rsp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb201_fixH2(%rsp),%xmm3
        movapd nb201_fiyH2(%rsp),%xmm4
        movapd nb201_fizH2(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb201_fixH2(%rsp)
        movlpd %xmm4,nb201_fiyH2(%rsp)
        movlpd %xmm7,nb201_fizH2(%rsp)

        movq nb201_faction(%rbp),%rdi
        ## update j forces 
        addsd  nb201_fjx(%rsp),%xmm0
        addsd  nb201_fjy(%rsp),%xmm1
        addsd  nb201_fjz(%rsp),%xmm2
        movlpd (%rdi,%rax,8),%xmm3
        movlpd 8(%rdi,%rax,8),%xmm4
        movlpd 16(%rdi,%rax,8),%xmm5
        addsd %xmm0,%xmm3
        addsd %xmm1,%xmm4
        addsd %xmm2,%xmm5
        movlpd %xmm3,(%rdi,%rax,8)
        movlpd %xmm4,8(%rdi,%rax,8)
        movlpd %xmm5,16(%rdi,%rax,8)

_nb_kernel201_x86_64_sse2.nb201_updateouterdata: 
        movl  nb201_ii3(%rsp),%ecx
        movq  nb201_faction(%rbp),%rdi
        movq  nb201_fshift(%rbp),%rsi
        movl  nb201_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb201_fixO(%rsp),%xmm0
        movapd nb201_fiyO(%rsp),%xmm1
        movapd nb201_fizO(%rsp),%xmm2

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
        movsd  (%rdi,%rcx,8),%xmm3
        movsd  8(%rdi,%rcx,8),%xmm4
        movsd  16(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,(%rdi,%rcx,8)
        movsd  %xmm4,8(%rdi,%rcx,8)
        movsd  %xmm5,16(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        movapd %xmm0,%xmm6
        movsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm6

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movapd nb201_fixH1(%rsp),%xmm0
        movapd nb201_fiyH1(%rsp),%xmm1
        movapd nb201_fizH1(%rsp),%xmm2

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
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movapd nb201_fixH2(%rsp),%xmm0
        movapd nb201_fiyH2(%rsp),%xmm1
        movapd nb201_fizH2(%rsp),%xmm2

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
        movl nb201_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb201_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb201_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb201_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb201_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel201_x86_64_sse2.nb201_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb201_n(%rsp)
        jmp _nb_kernel201_x86_64_sse2.nb201_outer
_nb_kernel201_x86_64_sse2.nb201_outerend: 
        ## check if more outer neighborlists remain
        movl  nb201_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel201_x86_64_sse2.nb201_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel201_x86_64_sse2.nb201_threadloop
_nb_kernel201_x86_64_sse2.nb201_end: 
        movl nb201_nouter(%rsp),%eax
        movl nb201_ninner(%rsp),%ebx
        movq nb201_outeriter(%rbp),%rcx
        movq nb201_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $808,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret





.globl nb_kernel201nf_x86_64_sse2
.globl _nb_kernel201nf_x86_64_sse2
nb_kernel201nf_x86_64_sse2:     
_nb_kernel201nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb201nf_fshift, 16
.set nb201nf_gid, 24
.set nb201nf_pos, 32
.set nb201nf_faction, 40
.set nb201nf_charge, 48
.set nb201nf_p_facel, 56
.set nb201nf_argkrf, 64
.set nb201nf_argcrf, 72
.set nb201nf_Vc, 80
.set nb201nf_type, 88
.set nb201nf_p_ntype, 96
.set nb201nf_vdwparam, 104
.set nb201nf_Vvdw, 112
.set nb201nf_p_tabscale, 120
.set nb201nf_VFtab, 128
.set nb201nf_invsqrta, 136
.set nb201nf_dvda, 144
.set nb201nf_p_gbtabscale, 152
.set nb201nf_GBtab, 160
.set nb201nf_p_nthreads, 168
.set nb201nf_count, 176
.set nb201nf_mtx, 184
.set nb201nf_outeriter, 192
.set nb201nf_inneriter, 200
.set nb201nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb201nf_ixO, 0
.set nb201nf_iyO, 16
.set nb201nf_izO, 32
.set nb201nf_ixH1, 48
.set nb201nf_iyH1, 64
.set nb201nf_izH1, 80
.set nb201nf_ixH2, 96
.set nb201nf_iyH2, 112
.set nb201nf_izH2, 128
.set nb201nf_iqO, 144
.set nb201nf_iqH, 160
.set nb201nf_qqO, 176
.set nb201nf_qqH, 192
.set nb201nf_vctot, 208
.set nb201nf_half, 224
.set nb201nf_three, 240
.set nb201nf_krf, 256
.set nb201nf_crf, 272
.set nb201nf_krsqO, 288
.set nb201nf_krsqH1, 304
.set nb201nf_krsqH2, 320
.set nb201nf_nri, 336
.set nb201nf_iinr, 344
.set nb201nf_jindex, 352
.set nb201nf_jjnr, 360
.set nb201nf_shift, 368
.set nb201nf_shiftvec, 376
.set nb201nf_facel, 384
.set nb201nf_innerjjnr, 392
.set nb201nf_is3, 400
.set nb201nf_ii3, 404
.set nb201nf_innerk, 408
.set nb201nf_n, 412
.set nb201nf_nn1, 416
.set nb201nf_nouter, 420
.set nb201nf_ninner, 424
        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $440,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb201nf_nouter(%rsp)
        movl %eax,nb201nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb201nf_nri(%rsp)
        movq %rsi,nb201nf_iinr(%rsp)
        movq %rdx,nb201nf_jindex(%rsp)
        movq %rcx,nb201nf_jjnr(%rsp)
        movq %r8,nb201nf_shift(%rsp)
        movq %r9,nb201nf_shiftvec(%rsp)
        movq nb201nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb201nf_facel(%rsp)

        movq nb201nf_argkrf(%rbp),%rsi
        movq nb201nf_argcrf(%rbp),%rdi
        movsd (%rsi),%xmm1
        movsd (%rdi),%xmm2
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        movapd %xmm1,nb201nf_krf(%rsp)
        movapd %xmm2,nb201nf_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb201nf_half(%rsp)
        movl %ebx,nb201nf_half+4(%rsp)
        movsd nb201nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb201nf_half(%rsp)
        movapd %xmm3,nb201nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb201nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb201nf_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        movsd 8(%rdx,%rbx,8),%xmm4

        movsd nb201nf_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb201nf_iqO(%rsp)
        movapd %xmm4,nb201nf_iqH(%rsp)

_nb_kernel201nf_x86_64_sse2.nb201nf_threadloop: 
        movq  nb201nf_count(%rbp),%rsi          ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel201nf_x86_64_sse2.nb201nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel201nf_x86_64_sse2.nb201nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb201nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb201nf_n(%rsp)
        movl %ebx,nb201nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel201nf_x86_64_sse2.nb201nf_outerstart
        jmp _nb_kernel201nf_x86_64_sse2.nb201nf_end

_nb_kernel201nf_x86_64_sse2.nb201nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb201nf_nouter(%rsp),%ebx
        movl %ebx,nb201nf_nouter(%rsp)

_nb_kernel201nf_x86_64_sse2.nb201nf_outer: 
        movq  nb201nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 

        movq  nb201nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb201nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb201nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb201nf_ii3(%rsp)

        addsd (%rax,%rbx,8),%xmm3
        addsd 8(%rax,%rbx,8),%xmm4
        addsd 16(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb201nf_ixO(%rsp)
        movapd %xmm4,nb201nf_iyO(%rsp)
        movapd %xmm5,nb201nf_izO(%rsp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 24(%rax,%rbx,8),%xmm0
        addsd 32(%rax,%rbx,8),%xmm1
        addsd 40(%rax,%rbx,8),%xmm2
        addsd 48(%rax,%rbx,8),%xmm3
        addsd 56(%rax,%rbx,8),%xmm4
        addsd 64(%rax,%rbx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb201nf_ixH1(%rsp)
        movapd %xmm1,nb201nf_iyH1(%rsp)
        movapd %xmm2,nb201nf_izH1(%rsp)
        movapd %xmm3,nb201nf_ixH2(%rsp)
        movapd %xmm4,nb201nf_iyH2(%rsp)
        movapd %xmm5,nb201nf_izH2(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb201nf_vctot(%rsp)

        movq  nb201nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx     ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb201nf_pos(%rbp),%rsi
        movq  nb201nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb201nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb201nf_ninner(%rsp),%ecx
        movl  %ecx,nb201nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb201nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel201nf_x86_64_sse2.nb201nf_unroll_loop
        jmp   _nb_kernel201nf_x86_64_sse2.nb201nf_checksingle
_nb_kernel201nf_x86_64_sse2.nb201nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb201nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb201nf_innerjjnr(%rsp)                 ## advance pointer (unrolled 2) 

        movq nb201nf_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb201nf_iqO(%rsp),%xmm3
        mulpd  nb201nf_iqH(%rsp),%xmm4
        movapd  %xmm3,nb201nf_qqO(%rsp)
        movapd  %xmm4,nb201nf_qqH(%rsp)

        movq nb201nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move two coordinates to xmm0-xmm2 
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2
        movhpd (%rsi,%rbx,8),%xmm0
        movhpd 8(%rsi,%rbx,8),%xmm1
        movhpd 16(%rsi,%rbx,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb201nf_ixO(%rsp),%xmm4
        movapd nb201nf_iyO(%rsp),%xmm5
        movapd nb201nf_izO(%rsp),%xmm6

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
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb201nf_ixH1(%rsp),%xmm4
        movapd nb201nf_iyH1(%rsp),%xmm5
        movapd nb201nf_izH1(%rsp),%xmm6

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
        movapd nb201nf_ixH2(%rsp),%xmm3
        movapd nb201nf_iyH2(%rsp),%xmm4
        movapd nb201nf_izH2(%rsp),%xmm5

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
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        movapd %xmm5,%xmm0
        movapd %xmm6,%xmm1
        movapd %xmm7,%xmm2

        mulpd  nb201nf_krf(%rsp),%xmm0
        mulpd  nb201nf_krf(%rsp),%xmm1
        mulpd  nb201nf_krf(%rsp),%xmm2

        movapd %xmm0,nb201nf_krsqH2(%rsp)
        movapd %xmm1,nb201nf_krsqH1(%rsp)
        movapd %xmm2,nb201nf_krsqO(%rsp)

        ## start with rsqO - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb201nf_three(%rsp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb201nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb201nf_three(%rsp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb201nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb201nf_three(%rsp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb201nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb201nf_three(%rsp),%xmm4
        subpd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb201nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb201nf_three(%rsp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb201nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb201nf_three(%rsp),%xmm4
        subpd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb201nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movapd  nb201nf_krsqO(%rsp),%xmm0
        addpd   %xmm0,%xmm7     ## xmm7=rinv+ krsq 
        subpd   nb201nf_crf(%rsp),%xmm7
        mulpd   nb201nf_qqO(%rsp),%xmm7   ## vcoul      
        addpd  nb201nf_vctot(%rsp),%xmm7

        ## H1 interactions 
        movapd  nb201nf_krsqH1(%rsp),%xmm0
        addpd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        subpd   nb201nf_crf(%rsp),%xmm6
        mulpd   nb201nf_qqH(%rsp),%xmm6   ## vcoul 
        addpd  %xmm7,%xmm6

        ## H2 interactions 
        movapd  nb201nf_krsqH2(%rsp),%xmm0
        addpd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        subpd   nb201nf_crf(%rsp),%xmm5
        mulpd   nb201nf_qqH(%rsp),%xmm5   ## vcoul 
        addpd  %xmm6,%xmm5
        movapd %xmm5,nb201nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb201nf_innerk(%rsp)
        jl    _nb_kernel201nf_x86_64_sse2.nb201nf_checksingle
        jmp   _nb_kernel201nf_x86_64_sse2.nb201nf_unroll_loop
_nb_kernel201nf_x86_64_sse2.nb201nf_checksingle: 
        movl  nb201nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel201nf_x86_64_sse2.nb201nf_dosingle
        jmp   _nb_kernel201nf_x86_64_sse2.nb201nf_updateouterdata
_nb_kernel201nf_x86_64_sse2.nb201nf_dosingle: 
        movq  nb201nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb201nf_innerjjnr(%rsp)

        movq nb201nf_charge(%rbp),%rsi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb201nf_iqO(%rsp),%xmm3
        mulpd  nb201nf_iqH(%rsp),%xmm4
        movapd  %xmm3,nb201nf_qqO(%rsp)
        movapd  %xmm4,nb201nf_qqH(%rsp)

        movq nb201nf_pos(%rbp),%rsi        ## base of pos[] 
        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb201nf_ixO(%rsp),%xmm4
        movapd nb201nf_iyO(%rsp),%xmm5
        movapd nb201nf_izO(%rsp),%xmm6

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
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb201nf_ixH1(%rsp),%xmm4
        movapd nb201nf_iyH1(%rsp),%xmm5
        movapd nb201nf_izH1(%rsp),%xmm6

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
        movapd nb201nf_ixH2(%rsp),%xmm3
        movapd nb201nf_iyH2(%rsp),%xmm4
        movapd nb201nf_izH2(%rsp),%xmm5

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
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        movapd %xmm5,%xmm0
        movapd %xmm6,%xmm1
        movapd %xmm7,%xmm2

        mulsd  nb201nf_krf(%rsp),%xmm0
        mulsd  nb201nf_krf(%rsp),%xmm1
        mulsd  nb201nf_krf(%rsp),%xmm2

        movapd %xmm0,nb201nf_krsqH2(%rsp)
        movapd %xmm1,nb201nf_krsqH1(%rsp)
        movapd %xmm2,nb201nf_krsqO(%rsp)

        ## start with rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb201nf_three(%rsp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb201nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb201nf_three(%rsp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb201nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb201nf_three(%rsp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb201nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb201nf_three(%rsp),%xmm4
        subsd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb201nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb201nf_three(%rsp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb201nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb201nf_three(%rsp),%xmm4
        subsd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb201nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movapd  nb201nf_krsqO(%rsp),%xmm0
        addsd   %xmm0,%xmm7     ## xmm7=rinv+ krsq 
        subsd   nb201nf_crf(%rsp),%xmm7
        mulsd   nb201nf_qqO(%rsp),%xmm7   ## vcoul      
        addsd  nb201nf_vctot(%rsp),%xmm7

        ## H1 interactions 
        movapd  nb201nf_krsqH1(%rsp),%xmm0
        addsd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        subsd   nb201nf_crf(%rsp),%xmm6
        mulsd   nb201nf_qqH(%rsp),%xmm6   ## vcoul 
        addsd  %xmm7,%xmm6

        ## H2 interactions 
        movapd  nb201nf_krsqH2(%rsp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        subsd   nb201nf_crf(%rsp),%xmm5
        mulsd   nb201nf_qqH(%rsp),%xmm5   ## vcoul 
        addsd  %xmm6,%xmm5
        movlpd %xmm5,nb201nf_vctot(%rsp)

_nb_kernel201nf_x86_64_sse2.nb201nf_updateouterdata: 
        ## get n from stack
        movl nb201nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb201nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb201nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb201nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb201nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel201nf_x86_64_sse2.nb201nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb201nf_n(%rsp)
        jmp _nb_kernel201nf_x86_64_sse2.nb201nf_outer
_nb_kernel201nf_x86_64_sse2.nb201nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb201nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel201nf_x86_64_sse2.nb201nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel201nf_x86_64_sse2.nb201nf_threadloop
_nb_kernel201nf_x86_64_sse2.nb201nf_end: 
        movl nb201nf_nouter(%rsp),%eax
        movl nb201nf_ninner(%rsp),%ebx
        movq nb201nf_outeriter(%rbp),%rcx
        movq nb201nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $440,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret

