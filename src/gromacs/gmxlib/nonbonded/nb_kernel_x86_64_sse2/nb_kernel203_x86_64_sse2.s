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







.globl nb_kernel203_x86_64_sse2
.globl _nb_kernel203_x86_64_sse2
nb_kernel203_x86_64_sse2:       
_nb_kernel203_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb203_fshift, 16
.set nb203_gid, 24
.set nb203_pos, 32
.set nb203_faction, 40
.set nb203_charge, 48
.set nb203_p_facel, 56
.set nb203_argkrf, 64
.set nb203_argcrf, 72
.set nb203_Vc, 80
.set nb203_type, 88
.set nb203_p_ntype, 96
.set nb203_vdwparam, 104
.set nb203_Vvdw, 112
.set nb203_p_tabscale, 120
.set nb203_VFtab, 128
.set nb203_invsqrta, 136
.set nb203_dvda, 144
.set nb203_p_gbtabscale, 152
.set nb203_GBtab, 160
.set nb203_p_nthreads, 168
.set nb203_count, 176
.set nb203_mtx, 184
.set nb203_outeriter, 192
.set nb203_inneriter, 200
.set nb203_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb203_ixM, 0
.set nb203_iyM, 16
.set nb203_izM, 32
.set nb203_ixH1, 48
.set nb203_iyH1, 64
.set nb203_izH1, 80
.set nb203_ixH2, 96
.set nb203_iyH2, 112
.set nb203_izH2, 128
.set nb203_iqM, 144
.set nb203_iqH, 160
.set nb203_dxM, 176
.set nb203_dyM, 192
.set nb203_dzM, 208
.set nb203_dxH1, 224
.set nb203_dyH1, 240
.set nb203_dzH1, 256
.set nb203_dxH2, 272
.set nb203_dyH2, 288
.set nb203_dzH2, 304
.set nb203_qqM, 320
.set nb203_qqH, 336
.set nb203_vctot, 352
.set nb203_fixM, 384
.set nb203_fiyM, 400
.set nb203_fizM, 416
.set nb203_fixH1, 432
.set nb203_fiyH1, 448
.set nb203_fizH1, 464
.set nb203_fixH2, 480
.set nb203_fiyH2, 496
.set nb203_fizH2, 512
.set nb203_fjx, 528
.set nb203_fjy, 544
.set nb203_fjz, 560
.set nb203_half, 576
.set nb203_three, 592
.set nb203_two, 608
.set nb203_krf, 624
.set nb203_crf, 640
.set nb203_krsqM, 656
.set nb203_krsqH1, 672
.set nb203_krsqH2, 688
.set nb203_is3, 704
.set nb203_ii3, 708
.set nb203_nri, 712
.set nb203_iinr, 720
.set nb203_jindex, 728
.set nb203_jjnr, 736
.set nb203_shift, 744
.set nb203_shiftvec, 752
.set nb203_facel, 760
.set nb203_innerjjnr, 768
.set nb203_innerk, 776
.set nb203_n, 780
.set nb203_nn1, 784
.set nb203_nouter, 788
.set nb203_ninner, 792
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
        movl %eax,nb203_nouter(%rsp)
        movl %eax,nb203_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb203_nri(%rsp)
        movq %rsi,nb203_iinr(%rsp)
        movq %rdx,nb203_jindex(%rsp)
        movq %rcx,nb203_jjnr(%rsp)
        movq %r8,nb203_shift(%rsp)
        movq %r9,nb203_shiftvec(%rsp)
        movq nb203_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb203_facel(%rsp)

        movq nb203_argkrf(%rbp),%rsi
        movq nb203_argcrf(%rbp),%rdi
        movsd (%rsi),%xmm1
        movsd (%rdi),%xmm2
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        movapd %xmm1,nb203_krf(%rsp)
        movapd %xmm2,nb203_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb203_half(%rsp)
        movl %ebx,nb203_half+4(%rsp)
        movsd nb203_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb203_half(%rsp)
        movapd %xmm2,nb203_two(%rsp)
        movapd %xmm3,nb203_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb203_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb203_charge(%rbp),%rdx
        movsd 8(%rdx,%rbx,8),%xmm3
        movsd 24(%rdx,%rbx,8),%xmm4

        movsd nb203_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb203_iqH(%rsp)
        movapd %xmm4,nb203_iqM(%rsp)

_nb_kernel203_x86_64_sse2.nb203_threadloop: 
        movq  nb203_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel203_x86_64_sse2.nb203_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel203_x86_64_sse2.nb203_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb203_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb203_n(%rsp)
        movl %ebx,nb203_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel203_x86_64_sse2.nb203_outerstart
        jmp _nb_kernel203_x86_64_sse2.nb203_end

_nb_kernel203_x86_64_sse2.nb203_outerstart: 
        ## ebx contains number of outer iterations
        addl nb203_nouter(%rsp),%ebx
        movl %ebx,nb203_nouter(%rsp)

_nb_kernel203_x86_64_sse2.nb203_outer: 
        movq  nb203_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb203_is3(%rsp)      ## store is3 

        movq  nb203_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb203_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb203_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb203_ii3(%rsp)

        addsd 24(%rax,%rbx,8),%xmm3
        addsd 32(%rax,%rbx,8),%xmm4
        addsd 40(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb203_ixH1(%rsp)
        movapd %xmm4,nb203_iyH1(%rsp)
        movapd %xmm5,nb203_izH1(%rsp)

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
        movapd %xmm0,nb203_ixH2(%rsp)
        movapd %xmm1,nb203_iyH2(%rsp)
        movapd %xmm2,nb203_izH2(%rsp)
        movapd %xmm3,nb203_ixM(%rsp)
        movapd %xmm4,nb203_iyM(%rsp)
        movapd %xmm5,nb203_izM(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb203_vctot(%rsp)
        movapd %xmm4,nb203_fixM(%rsp)
        movapd %xmm4,nb203_fiyM(%rsp)
        movapd %xmm4,nb203_fizM(%rsp)
        movapd %xmm4,nb203_fixH1(%rsp)
        movapd %xmm4,nb203_fiyH1(%rsp)
        movapd %xmm4,nb203_fizH1(%rsp)
        movapd %xmm4,nb203_fixH2(%rsp)
        movapd %xmm4,nb203_fiyH2(%rsp)
        movapd %xmm4,nb203_fizH2(%rsp)

        movq  nb203_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb203_pos(%rbp),%rsi
        movq  nb203_faction(%rbp),%rdi
        movq  nb203_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb203_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb203_ninner(%rsp),%ecx
        movl  %ecx,nb203_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb203_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel203_x86_64_sse2.nb203_unroll_loop
        jmp   _nb_kernel203_x86_64_sse2.nb203_checksingle
_nb_kernel203_x86_64_sse2.nb203_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb203_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb203_innerjjnr(%rsp)                   ## advance pointer (unrolled 2) 

        movq nb203_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb203_iqM(%rsp),%xmm3
        mulpd  nb203_iqH(%rsp),%xmm4
        movapd  %xmm3,nb203_qqM(%rsp)
        movapd  %xmm4,nb203_qqH(%rsp)

        movq nb203_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## load j coordinates 
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

    subpd nb203_ixH1(%rsp),%xmm0
    subpd nb203_iyH1(%rsp),%xmm1
    subpd nb203_izH1(%rsp),%xmm2
    subpd nb203_ixH2(%rsp),%xmm3
    subpd nb203_iyH2(%rsp),%xmm4
    subpd nb203_izH2(%rsp),%xmm5
    subpd nb203_ixM(%rsp),%xmm6
    subpd nb203_iyM(%rsp),%xmm7
    subpd nb203_izM(%rsp),%xmm8

        movapd %xmm0,nb203_dxH1(%rsp)
        movapd %xmm1,nb203_dyH1(%rsp)
        movapd %xmm2,nb203_dzH1(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb203_dxH2(%rsp)
        movapd %xmm4,nb203_dyH2(%rsp)
        movapd %xmm5,nb203_dzH2(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb203_dxM(%rsp)
        movapd %xmm7,nb203_dyM(%rsp)
        movapd %xmm8,nb203_dzM(%rsp)
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

        movapd  nb203_three(%rsp),%xmm9
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

        movapd  nb203_half(%rsp),%xmm15
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

        movapd  nb203_three(%rsp),%xmm1
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

        movapd  nb203_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvH1
        mulpd   %xmm15,%xmm10 ##   rinvH2
    mulpd   %xmm15,%xmm11 ##   rinvM

        ## interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd %xmm9,%xmm1 ## copy of rinv
    movapd %xmm10,%xmm4
    movapd %xmm11,%xmm7
    movapd nb203_krf(%rsp),%xmm2
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
    movapd nb203_crf(%rsp),%xmm14
    subpd  %xmm14,%xmm2  ## rinv+krsq-crf
    subpd  %xmm14,%xmm5
    subpd  %xmm14,%xmm8
    movapd nb203_qqH(%rsp),%xmm12
    movapd nb203_qqM(%rsp),%xmm13
    mulpd  %xmm12,%xmm2 ## voul=qq*(rinv+ krsq-crf)
    mulpd  %xmm12,%xmm5 ## voul=qq*(rinv+ krsq-crf)
    mulpd  %xmm13,%xmm8 ## voul=qq*(rinv+ krsq-crf)
    addpd  %xmm0,%xmm0 ## 2*krsq
    addpd  %xmm3,%xmm3
    addpd  %xmm6,%xmm6
    subpd  %xmm0,%xmm1 ## rinv-2*krsq
    subpd  %xmm3,%xmm4
    subpd  %xmm6,%xmm7
    mulpd  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulpd  %xmm12,%xmm4
    mulpd  %xmm13,%xmm7
    addpd  nb203_vctot(%rsp),%xmm2
    addpd  %xmm8,%xmm5
    addpd  %xmm5,%xmm2
    movapd %xmm2,nb203_vctot(%rsp)

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

        mulpd nb203_dxH1(%rsp),%xmm7
        mulpd nb203_dyH1(%rsp),%xmm8
        mulpd nb203_dzH1(%rsp),%xmm9
        mulpd nb203_dxH2(%rsp),%xmm10
        mulpd nb203_dyH2(%rsp),%xmm11
        mulpd nb203_dzH2(%rsp),%xmm12
        mulpd nb203_dxM(%rsp),%xmm13
        mulpd nb203_dyM(%rsp),%xmm14
        mulpd nb203_dzM(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb203_fixH1(%rsp),%xmm7
    addpd nb203_fiyH1(%rsp),%xmm8
    addpd nb203_fizH1(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb203_fixH2(%rsp),%xmm10
    addpd nb203_fiyH2(%rsp),%xmm11
    addpd nb203_fizH2(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb203_fixM(%rsp),%xmm13
    addpd nb203_fiyM(%rsp),%xmm14
    addpd nb203_fizM(%rsp),%xmm15

    movapd %xmm7,nb203_fixH1(%rsp)
    movapd %xmm8,nb203_fiyH1(%rsp)
    movapd %xmm9,nb203_fizH1(%rsp)
    movapd %xmm10,nb203_fixH2(%rsp)
    movapd %xmm11,nb203_fiyH2(%rsp)
    movapd %xmm12,nb203_fizH2(%rsp)
    movapd %xmm13,nb203_fixM(%rsp)
    movapd %xmm14,nb203_fiyM(%rsp)
    movapd %xmm15,nb203_fizM(%rsp)

    ## store back j forces from xmm0-xmm2
        movlpd %xmm0,(%rdi,%rax,8)
        movlpd %xmm1,8(%rdi,%rax,8)
        movlpd %xmm2,16(%rdi,%rax,8)
        movhpd %xmm0,(%rdi,%rbx,8)
        movhpd %xmm1,8(%rdi,%rbx,8)
        movhpd %xmm2,16(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb203_innerk(%rsp)
        jl    _nb_kernel203_x86_64_sse2.nb203_checksingle
        jmp   _nb_kernel203_x86_64_sse2.nb203_unroll_loop
_nb_kernel203_x86_64_sse2.nb203_checksingle: 
        movl  nb203_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel203_x86_64_sse2.nb203_dosingle
        jmp   _nb_kernel203_x86_64_sse2.nb203_updateouterdata
_nb_kernel203_x86_64_sse2.nb203_dosingle: 
        movq  nb203_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb203_innerjjnr(%rsp)

        movq nb203_charge(%rbp),%rsi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb203_iqM(%rsp),%xmm3
        mulpd  nb203_iqH(%rsp),%xmm4
        movapd  %xmm3,nb203_qqM(%rsp)
        movapd  %xmm4,nb203_qqH(%rsp)

        movq nb203_pos(%rbp),%rsi        ## base of pos[] 
        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coordinates to xmm4-xmm6 & xmm0-xmm2    
        movlpd (%rsi,%rax,8),%xmm4
        movlpd 8(%rsi,%rax,8),%xmm5
        movlpd 16(%rsi,%rax,8),%xmm6
    movapd %xmm4,%xmm0
    movapd %xmm5,%xmm1
    movapd %xmm6,%xmm2

        ## calc dr 
        subsd nb203_ixM(%rsp),%xmm4
        subsd nb203_iyM(%rsp),%xmm5
        subsd nb203_izM(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb203_dxM(%rsp)
        movapd %xmm5,nb203_dyM(%rsp)
        movapd %xmm6,nb203_dzM(%rsp)
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
        subsd nb203_ixH1(%rsp),%xmm4
        subsd nb203_iyH1(%rsp),%xmm5
        subsd nb203_izH1(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb203_dxH1(%rsp)
        movapd %xmm5,nb203_dyH1(%rsp)
        movapd %xmm6,nb203_dzH1(%rsp)
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
        subsd nb203_ixH2(%rsp),%xmm3
        subsd nb203_iyH2(%rsp),%xmm4
        subsd nb203_izH2(%rsp),%xmm5

        ## store dr 
        movapd %xmm3,nb203_dxH2(%rsp)
        movapd %xmm4,nb203_dyH2(%rsp)
        movapd %xmm5,nb203_dzH2(%rsp)
        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqM in xmm7 

        movapd %xmm5,%xmm0
        movapd %xmm6,%xmm1
        movapd %xmm7,%xmm2

        mulsd  nb203_krf(%rsp),%xmm0
        mulsd  nb203_krf(%rsp),%xmm1
        mulsd  nb203_krf(%rsp),%xmm2

        movapd %xmm0,nb203_krsqH2(%rsp)
        movapd %xmm1,nb203_krsqH1(%rsp)
        movapd %xmm2,nb203_krsqM(%rsp)

        ## start with rsqM - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb203_three(%rsp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb203_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb203_three(%rsp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb203_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvM in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb203_three(%rsp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb203_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb203_three(%rsp),%xmm4
        subsd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb203_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb203_three(%rsp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb203_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb203_three(%rsp),%xmm4
        subsd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb203_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do M interactions 
        movapd  %xmm7,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movapd  %xmm7,%xmm3
        movapd  nb203_krsqM(%rsp),%xmm0
        addsd   %xmm0,%xmm7     ## xmm6=rinv+ krsq 
        mulsd   nb203_two(%rsp),%xmm0
        subsd   nb203_crf(%rsp),%xmm7
        subsd   %xmm0,%xmm3     ## xmm7=rinv-2*krsq 
        mulsd   nb203_qqM(%rsp),%xmm7   ## vcoul 
        mulsd   nb203_qqM(%rsp),%xmm3
        mulsd  %xmm3,%xmm4      ## total fsH1 in xmm4 

        addsd  nb203_vctot(%rsp),%xmm7

        movapd nb203_dxM(%rsp),%xmm0
        movapd nb203_dyM(%rsp),%xmm1
        movapd nb203_dzM(%rsp),%xmm2
        movlpd %xmm7,nb203_vctot(%rsp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update M forces 
        movapd nb203_fixM(%rsp),%xmm3
        movapd nb203_fiyM(%rsp),%xmm4
        movapd nb203_fizM(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb203_fixM(%rsp)
        movlpd %xmm4,nb203_fiyM(%rsp)
        movlpd %xmm7,nb203_fizM(%rsp)
        ## update j forces with water M 
        movlpd %xmm0,nb203_fjx(%rsp)
        movlpd %xmm1,nb203_fjy(%rsp)
        movlpd %xmm2,nb203_fjz(%rsp)

        ## H1 interactions 
        movapd  %xmm6,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movapd  %xmm6,%xmm7
        movapd  nb203_krsqH1(%rsp),%xmm0
        addsd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        mulsd   nb203_two(%rsp),%xmm0
        subsd   nb203_crf(%rsp),%xmm6
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb203_qqH(%rsp),%xmm6   ## vcoul 
        mulsd   nb203_qqH(%rsp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH1 in xmm4 

        addsd  nb203_vctot(%rsp),%xmm6

        movapd nb203_dxH1(%rsp),%xmm0
        movapd nb203_dyH1(%rsp),%xmm1
        movapd nb203_dzH1(%rsp),%xmm2
        movlpd %xmm6,nb203_vctot(%rsp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb203_fixH1(%rsp),%xmm3
        movapd nb203_fiyH1(%rsp),%xmm4
        movapd nb203_fizH1(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb203_fixH1(%rsp)
        movlpd %xmm4,nb203_fiyH1(%rsp)
        movlpd %xmm7,nb203_fizH1(%rsp)
        ## update j forces with water H1 
        addsd  nb203_fjx(%rsp),%xmm0
        addsd  nb203_fjy(%rsp),%xmm1
        addsd  nb203_fjz(%rsp),%xmm2
        movlpd %xmm0,nb203_fjx(%rsp)
        movlpd %xmm1,nb203_fjy(%rsp)
        movlpd %xmm2,nb203_fjz(%rsp)

        ## H2 interactions 
        movapd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movapd  %xmm5,%xmm7
        movapd  nb203_krsqH2(%rsp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        mulsd   nb203_two(%rsp),%xmm0
        subsd   nb203_crf(%rsp),%xmm5
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb203_qqH(%rsp),%xmm5   ## vcoul 
        mulsd   nb203_qqH(%rsp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addsd  nb203_vctot(%rsp),%xmm5

        movapd nb203_dxH2(%rsp),%xmm0
        movapd nb203_dyH2(%rsp),%xmm1
        movapd nb203_dzH2(%rsp),%xmm2
        movlpd %xmm5,nb203_vctot(%rsp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb203_fixH2(%rsp),%xmm3
        movapd nb203_fiyH2(%rsp),%xmm4
        movapd nb203_fizH2(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb203_fixH2(%rsp)
        movlpd %xmm4,nb203_fiyH2(%rsp)
        movlpd %xmm7,nb203_fizH2(%rsp)

        movq nb203_faction(%rbp),%rdi
        ## update j forces 
        addsd  nb203_fjx(%rsp),%xmm0
        addsd  nb203_fjy(%rsp),%xmm1
        addsd  nb203_fjz(%rsp),%xmm2
        movlpd (%rdi,%rax,8),%xmm3
        movlpd 8(%rdi,%rax,8),%xmm4
        movlpd 16(%rdi,%rax,8),%xmm5
        addsd %xmm0,%xmm3
        addsd %xmm1,%xmm4
        addsd %xmm2,%xmm5
        movlpd %xmm3,(%rdi,%rax,8)
        movlpd %xmm4,8(%rdi,%rax,8)
        movlpd %xmm5,16(%rdi,%rax,8)

_nb_kernel203_x86_64_sse2.nb203_updateouterdata: 
        movl  nb203_ii3(%rsp),%ecx
        movq  nb203_faction(%rbp),%rdi
        movq  nb203_fshift(%rbp),%rsi
        movl  nb203_is3(%rsp),%edx

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movapd nb203_fixH1(%rsp),%xmm0
        movapd nb203_fiyH1(%rsp),%xmm1
        movapd nb203_fizH1(%rsp),%xmm2

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
        movapd nb203_fixH2(%rsp),%xmm0
        movapd nb203_fiyH2(%rsp),%xmm1
        movapd nb203_fizH2(%rsp),%xmm2

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
        movapd nb203_fixM(%rsp),%xmm0
        movapd nb203_fiyM(%rsp),%xmm1
        movapd nb203_fizM(%rsp),%xmm2

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
        movl nb203_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb203_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb203_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb203_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb203_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel203_x86_64_sse2.nb203_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb203_n(%rsp)
        jmp _nb_kernel203_x86_64_sse2.nb203_outer
_nb_kernel203_x86_64_sse2.nb203_outerend: 
        ## check if more outer neighborlists remain
        movl  nb203_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel203_x86_64_sse2.nb203_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel203_x86_64_sse2.nb203_threadloop
_nb_kernel203_x86_64_sse2.nb203_end: 
        movl nb203_nouter(%rsp),%eax
        movl nb203_ninner(%rsp),%ebx
        movq nb203_outeriter(%rbp),%rcx
        movq nb203_inneriter(%rbp),%rdx
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




.globl nb_kernel203nf_x86_64_sse2
.globl _nb_kernel203nf_x86_64_sse2
nb_kernel203nf_x86_64_sse2:     
_nb_kernel203nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb203nf_fshift, 16
.set nb203nf_gid, 24
.set nb203nf_pos, 32
.set nb203nf_faction, 40
.set nb203nf_charge, 48
.set nb203nf_p_facel, 56
.set nb203nf_argkrf, 64
.set nb203nf_argcrf, 72
.set nb203nf_Vc, 80
.set nb203nf_type, 88
.set nb203nf_p_ntype, 96
.set nb203nf_vdwparam, 104
.set nb203nf_Vvdw, 112
.set nb203nf_p_tabscale, 120
.set nb203nf_VFtab, 128
.set nb203nf_invsqrta, 136
.set nb203nf_dvda, 144
.set nb203nf_p_gbtabscale, 152
.set nb203nf_GBtab, 160
.set nb203nf_p_nthreads, 168
.set nb203nf_count, 176
.set nb203nf_mtx, 184
.set nb203nf_outeriter, 192
.set nb203nf_inneriter, 200
.set nb203nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb203nf_ixM, 0
.set nb203nf_iyM, 16
.set nb203nf_izM, 32
.set nb203nf_ixH1, 48
.set nb203nf_iyH1, 64
.set nb203nf_izH1, 80
.set nb203nf_ixH2, 96
.set nb203nf_iyH2, 112
.set nb203nf_izH2, 128
.set nb203nf_iqM, 144
.set nb203nf_iqH, 160
.set nb203nf_qqM, 176
.set nb203nf_qqH, 192
.set nb203nf_vctot, 208
.set nb203nf_half, 224
.set nb203nf_three, 240
.set nb203nf_krf, 256
.set nb203nf_crf, 272
.set nb203nf_krsqM, 288
.set nb203nf_krsqH1, 304
.set nb203nf_krsqH2, 320
.set nb203nf_is3, 336
.set nb203nf_ii3, 340
.set nb203nf_nri, 344
.set nb203nf_iinr, 352
.set nb203nf_jindex, 360
.set nb203nf_jjnr, 368
.set nb203nf_shift, 376
.set nb203nf_shiftvec, 384
.set nb203nf_facel, 392
.set nb203nf_innerjjnr, 400
.set nb203nf_innerk, 408
.set nb203nf_n, 412
.set nb203nf_nn1, 416
.set nb203nf_nouter, 420
.set nb203nf_ninner, 424
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
        movl %eax,nb203nf_nouter(%rsp)
        movl %eax,nb203nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb203nf_nri(%rsp)
        movq %rsi,nb203nf_iinr(%rsp)
        movq %rdx,nb203nf_jindex(%rsp)
        movq %rcx,nb203nf_jjnr(%rsp)
        movq %r8,nb203nf_shift(%rsp)
        movq %r9,nb203nf_shiftvec(%rsp)
        movq nb203nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb203nf_facel(%rsp)

        movq nb203nf_argkrf(%rbp),%rsi
        movq nb203nf_argcrf(%rbp),%rdi
        movsd (%rsi),%xmm1
        movsd (%rdi),%xmm2
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        movapd %xmm1,nb203nf_krf(%rsp)
        movapd %xmm2,nb203nf_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb203nf_half(%rsp)
        movl %ebx,nb203nf_half+4(%rsp)
        movsd nb203nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb203nf_half(%rsp)
        movapd %xmm3,nb203nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb203nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb203nf_charge(%rbp),%rdx
        movsd 8(%rdx,%rbx,8),%xmm3
        movsd 24(%rdx,%rbx,8),%xmm4

        movsd nb203nf_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb203nf_iqH(%rsp)
        movapd %xmm4,nb203nf_iqM(%rsp)

_nb_kernel203nf_x86_64_sse2.nb203nf_threadloop: 
        movq  nb203nf_count(%rbp),%rsi          ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel203nf_x86_64_sse2.nb203nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel203nf_x86_64_sse2.nb203nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb203nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb203nf_n(%rsp)
        movl %ebx,nb203nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel203nf_x86_64_sse2.nb203nf_outerstart
        jmp _nb_kernel203nf_x86_64_sse2.nb203nf_end

_nb_kernel203nf_x86_64_sse2.nb203nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb203nf_nouter(%rsp),%ebx
        movl %ebx,nb203nf_nouter(%rsp)

_nb_kernel203nf_x86_64_sse2.nb203nf_outer: 
        movq  nb203nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 

        movq  nb203nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb203nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb203nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb203nf_ii3(%rsp)

        addsd 24(%rax,%rbx,8),%xmm3
        addsd 32(%rax,%rbx,8),%xmm4
        addsd 40(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb203nf_ixH1(%rsp)
        movapd %xmm4,nb203nf_iyH1(%rsp)
        movapd %xmm5,nb203nf_izH1(%rsp)

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
        movapd %xmm0,nb203nf_ixH2(%rsp)
        movapd %xmm1,nb203nf_iyH2(%rsp)
        movapd %xmm2,nb203nf_izH2(%rsp)
        movapd %xmm3,nb203nf_ixM(%rsp)
        movapd %xmm4,nb203nf_iyM(%rsp)
        movapd %xmm5,nb203nf_izM(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb203nf_vctot(%rsp)

        movq  nb203nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx     ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb203nf_pos(%rbp),%rsi
        movq  nb203nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb203nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb203nf_ninner(%rsp),%ecx
        movl  %ecx,nb203nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb203nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel203nf_x86_64_sse2.nb203nf_unroll_loop
        jmp   _nb_kernel203nf_x86_64_sse2.nb203nf_checksingle
_nb_kernel203nf_x86_64_sse2.nb203nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb203nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb203nf_innerjjnr(%rsp)                 ## advance pointer (unrolled 2) 

        movq nb203nf_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb203nf_iqM(%rsp),%xmm3
        mulpd  nb203nf_iqH(%rsp),%xmm4
        movapd  %xmm3,nb203nf_qqM(%rsp)
        movapd  %xmm4,nb203nf_qqH(%rsp)

        movq nb203nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movapd nb203nf_ixM(%rsp),%xmm4
        movapd nb203nf_iyM(%rsp),%xmm5
        movapd nb203nf_izM(%rsp),%xmm6

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
        movapd nb203nf_ixH1(%rsp),%xmm4
        movapd nb203nf_iyH1(%rsp),%xmm5
        movapd nb203nf_izH1(%rsp),%xmm6

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
        movapd nb203nf_ixH2(%rsp),%xmm3
        movapd nb203nf_iyH2(%rsp),%xmm4
        movapd nb203nf_izH2(%rsp),%xmm5

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

        movapd %xmm5,%xmm0
        movapd %xmm6,%xmm1
        movapd %xmm7,%xmm2

        mulpd  nb203nf_krf(%rsp),%xmm0
        mulpd  nb203nf_krf(%rsp),%xmm1
        mulpd  nb203nf_krf(%rsp),%xmm2

        movapd %xmm0,nb203nf_krsqH2(%rsp)
        movapd %xmm1,nb203nf_krsqH1(%rsp)
        movapd %xmm2,nb203nf_krsqM(%rsp)

        ## start with rsqM - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb203nf_three(%rsp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb203nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb203nf_three(%rsp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb203nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvM in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb203nf_three(%rsp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb203nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb203nf_three(%rsp),%xmm4
        subpd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb203nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb203nf_three(%rsp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb203nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb203nf_three(%rsp),%xmm4
        subpd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb203nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do M interactions 
        movapd  nb203nf_krsqM(%rsp),%xmm0
        addpd   %xmm0,%xmm7     ## xmm7=rinv+ krsq 
        subpd   nb203nf_crf(%rsp),%xmm7
        mulpd   nb203nf_qqM(%rsp),%xmm7   ## vcoul      
        addpd  nb203nf_vctot(%rsp),%xmm7

        ## H1 interactions 
        movapd  nb203nf_krsqH1(%rsp),%xmm0
        addpd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        subpd   nb203nf_crf(%rsp),%xmm6
        mulpd   nb203nf_qqH(%rsp),%xmm6   ## vcoul 
        addpd  %xmm7,%xmm6

        ## H2 interactions 
        movapd  nb203nf_krsqH2(%rsp),%xmm0
        addpd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        subpd   nb203nf_crf(%rsp),%xmm5
        mulpd   nb203nf_qqH(%rsp),%xmm5   ## vcoul 
        addpd  %xmm6,%xmm5
        movapd %xmm5,nb203nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb203nf_innerk(%rsp)
        jl    _nb_kernel203nf_x86_64_sse2.nb203nf_checksingle
        jmp   _nb_kernel203nf_x86_64_sse2.nb203nf_unroll_loop
_nb_kernel203nf_x86_64_sse2.nb203nf_checksingle: 
        movl  nb203nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel203nf_x86_64_sse2.nb203nf_dosingle
        jmp   _nb_kernel203nf_x86_64_sse2.nb203nf_updateouterdata
_nb_kernel203nf_x86_64_sse2.nb203nf_dosingle: 
        movq  nb203nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb203nf_innerjjnr(%rsp)

        movq nb203nf_charge(%rbp),%rsi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb203nf_iqM(%rsp),%xmm3
        mulpd  nb203nf_iqH(%rsp),%xmm4
        movapd  %xmm3,nb203nf_qqM(%rsp)
        movapd  %xmm4,nb203nf_qqH(%rsp)

        movq nb203nf_pos(%rbp),%rsi        ## base of pos[] 
        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move ixM-izM to xmm4-xmm6 
        movapd nb203nf_ixM(%rsp),%xmm4
        movapd nb203nf_iyM(%rsp),%xmm5
        movapd nb203nf_izM(%rsp),%xmm6

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
        movapd nb203nf_ixH1(%rsp),%xmm4
        movapd nb203nf_iyH1(%rsp),%xmm5
        movapd nb203nf_izH1(%rsp),%xmm6

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
        movapd nb203nf_ixH2(%rsp),%xmm3
        movapd nb203nf_iyH2(%rsp),%xmm4
        movapd nb203nf_izH2(%rsp),%xmm5

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

        movapd %xmm5,%xmm0
        movapd %xmm6,%xmm1
        movapd %xmm7,%xmm2

        mulsd  nb203nf_krf(%rsp),%xmm0
        mulsd  nb203nf_krf(%rsp),%xmm1
        mulsd  nb203nf_krf(%rsp),%xmm2

        movapd %xmm0,nb203nf_krsqH2(%rsp)
        movapd %xmm1,nb203nf_krsqH1(%rsp)
        movapd %xmm2,nb203nf_krsqM(%rsp)

        ## start with rsqM - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb203nf_three(%rsp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb203nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb203nf_three(%rsp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb203nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvM in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb203nf_three(%rsp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb203nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb203nf_three(%rsp),%xmm4
        subsd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb203nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb203nf_three(%rsp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb203nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb203nf_three(%rsp),%xmm4
        subsd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb203nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do M interactions 
        movapd  nb203nf_krsqM(%rsp),%xmm0
        addsd   %xmm0,%xmm7     ## xmm7=rinv+ krsq 
        subsd   nb203nf_crf(%rsp),%xmm7
        mulsd   nb203nf_qqM(%rsp),%xmm7   ## vcoul      
        addsd  nb203nf_vctot(%rsp),%xmm7

        ## H1 interactions 
        movapd  nb203nf_krsqH1(%rsp),%xmm0
        addsd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        subsd   nb203nf_crf(%rsp),%xmm6
        mulsd   nb203nf_qqH(%rsp),%xmm6   ## vcoul 
        addsd  %xmm7,%xmm6

        ## H2 interactions 
        movapd  nb203nf_krsqH2(%rsp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        subsd   nb203nf_crf(%rsp),%xmm5
        mulsd   nb203nf_qqH(%rsp),%xmm5   ## vcoul 
        addsd  %xmm6,%xmm5
        movlpd %xmm5,nb203nf_vctot(%rsp)

_nb_kernel203nf_x86_64_sse2.nb203nf_updateouterdata: 
        ## get n from stack
        movl nb203nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb203nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb203nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb203nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb203nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel203nf_x86_64_sse2.nb203nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb203nf_n(%rsp)
        jmp _nb_kernel203nf_x86_64_sse2.nb203nf_outer
_nb_kernel203nf_x86_64_sse2.nb203nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb203nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel203nf_x86_64_sse2.nb203nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel203nf_x86_64_sse2.nb203nf_threadloop
_nb_kernel203nf_x86_64_sse2.nb203nf_end: 
        movl nb203nf_nouter(%rsp),%eax
        movl nb203nf_ninner(%rsp),%ebx
        movq nb203nf_outeriter(%rbp),%rcx
        movq nb203nf_inneriter(%rbp),%rdx
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

