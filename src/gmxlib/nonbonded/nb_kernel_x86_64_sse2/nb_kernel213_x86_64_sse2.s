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





.globl nb_kernel213_x86_64_sse2
.globl _nb_kernel213_x86_64_sse2
nb_kernel213_x86_64_sse2:       
_nb_kernel213_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb213_fshift, 16
.set nb213_gid, 24
.set nb213_pos, 32
.set nb213_faction, 40
.set nb213_charge, 48
.set nb213_p_facel, 56
.set nb213_argkrf, 64
.set nb213_argcrf, 72
.set nb213_Vc, 80
.set nb213_type, 88
.set nb213_p_ntype, 96
.set nb213_vdwparam, 104
.set nb213_Vvdw, 112
.set nb213_p_tabscale, 120
.set nb213_VFtab, 128
.set nb213_invsqrta, 136
.set nb213_dvda, 144
.set nb213_p_gbtabscale, 152
.set nb213_GBtab, 160
.set nb213_p_nthreads, 168
.set nb213_count, 176
.set nb213_mtx, 184
.set nb213_outeriter, 192
.set nb213_inneriter, 200
.set nb213_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb213_ixO, 0
.set nb213_iyO, 16
.set nb213_izO, 32
.set nb213_ixH1, 48
.set nb213_iyH1, 64
.set nb213_izH1, 80
.set nb213_ixH2, 96
.set nb213_iyH2, 112
.set nb213_izH2, 128
.set nb213_ixM, 144
.set nb213_iyM, 160
.set nb213_izM, 176
.set nb213_iqH, 192
.set nb213_iqM, 208
.set nb213_dxO, 224
.set nb213_dyO, 240
.set nb213_dzO, 256
.set nb213_dxH1, 272
.set nb213_dyH1, 288
.set nb213_dzH1, 304
.set nb213_dxH2, 320
.set nb213_dyH2, 336
.set nb213_dzH2, 352
.set nb213_dxM, 368
.set nb213_dyM, 384
.set nb213_dzM, 400
.set nb213_qqH, 416
.set nb213_qqM, 432
.set nb213_c6, 448
.set nb213_c12, 464
.set nb213_six, 480
.set nb213_twelve, 496
.set nb213_vctot, 512
.set nb213_Vvdwtot, 528
.set nb213_fixO, 544
.set nb213_fiyO, 560
.set nb213_fizO, 576
.set nb213_fixH1, 592
.set nb213_fiyH1, 608
.set nb213_fizH1, 624
.set nb213_fixH2, 640
.set nb213_fiyH2, 656
.set nb213_fizH2, 672
.set nb213_fixM, 688
.set nb213_fiyM, 704
.set nb213_fizM, 720
.set nb213_fjx, 736
.set nb213_fjy, 752
.set nb213_fjz, 768
.set nb213_half, 784
.set nb213_three, 800
.set nb213_two, 816
.set nb213_rinvH1, 832
.set nb213_rinvH2, 848
.set nb213_rinvM, 864
.set nb213_krsqH1, 880
.set nb213_krsqH2, 896
.set nb213_krsqM, 912
.set nb213_krf, 928
.set nb213_crf, 944
.set nb213_is3, 960
.set nb213_ii3, 964
.set nb213_nri, 968
.set nb213_iinr, 976
.set nb213_jindex, 984
.set nb213_jjnr, 992
.set nb213_shift, 1000
.set nb213_shiftvec, 1008
.set nb213_facel, 1016
.set nb213_innerjjnr, 1024
.set nb213_ntia, 1032
.set nb213_innerk, 1036
.set nb213_n, 1040
.set nb213_nn1, 1044
.set nb213_nouter, 1048
.set nb213_ninner, 1052
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1064,%rsp         ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb213_nouter(%rsp)
        movl %eax,nb213_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb213_nri(%rsp)
        movq %rsi,nb213_iinr(%rsp)
        movq %rdx,nb213_jindex(%rsp)
        movq %rcx,nb213_jjnr(%rsp)
        movq %r8,nb213_shift(%rsp)
        movq %r9,nb213_shiftvec(%rsp)
        movq nb213_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb213_facel(%rsp)

        movq nb213_argkrf(%rbp),%rsi
        movq nb213_argcrf(%rbp),%rdi
        movsd (%rsi),%xmm1
        movsd (%rdi),%xmm2
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        movapd %xmm1,nb213_krf(%rsp)
        movapd %xmm2,nb213_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb213_half(%rsp)
        movl %ebx,nb213_half+4(%rsp)
        movsd nb213_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm3,%xmm4
        addpd  %xmm4,%xmm4      ## six
        movapd %xmm4,%xmm5
        addpd  %xmm5,%xmm5      ## twelve
        movapd %xmm1,nb213_half(%rsp)
        movapd %xmm2,nb213_two(%rsp)
        movapd %xmm3,nb213_three(%rsp)
        movapd %xmm4,nb213_six(%rsp)
        movapd %xmm5,nb213_twelve(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb213_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb213_charge(%rbp),%rdx
        movsd 8(%rdx,%rbx,8),%xmm3
        movsd 24(%rdx,%rbx,8),%xmm4

        movsd nb213_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb213_iqH(%rsp)
        movapd %xmm4,nb213_iqM(%rsp)

        movq  nb213_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb213_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb213_ntia(%rsp)
_nb_kernel213_x86_64_sse2.nb213_threadloop: 
        movq  nb213_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel213_x86_64_sse2.nb213_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel213_x86_64_sse2.nb213_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb213_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb213_n(%rsp)
        movl %ebx,nb213_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel213_x86_64_sse2.nb213_outerstart
        jmp _nb_kernel213_x86_64_sse2.nb213_end

_nb_kernel213_x86_64_sse2.nb213_outerstart: 
        ## ebx contains number of outer iterations
        addl nb213_nouter(%rsp),%ebx
        movl %ebx,nb213_nouter(%rsp)

_nb_kernel213_x86_64_sse2.nb213_outer: 
        movq  nb213_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb213_is3(%rsp)      ## store is3 

        movq  nb213_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb213_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb213_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb213_ii3(%rsp)

        addsd (%rax,%rbx,8),%xmm3       ## ox
        addsd 8(%rax,%rbx,8),%xmm4      ## oy
        addsd 16(%rax,%rbx,8),%xmm5     ## oz   
        addsd 24(%rax,%rbx,8),%xmm6     ## h1x
        addsd 32(%rax,%rbx,8),%xmm7     ## h1y
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm7,%xmm7
        movapd %xmm3,nb213_ixO(%rsp)
        movapd %xmm4,nb213_iyO(%rsp)
        movapd %xmm5,nb213_izO(%rsp)
        movapd %xmm6,nb213_ixH1(%rsp)
        movapd %xmm7,nb213_iyH1(%rsp)

        movsd %xmm2,%xmm6
        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 40(%rax,%rbx,8),%xmm6    ## h1z
        addsd 48(%rax,%rbx,8),%xmm0    ## h2x
        addsd 56(%rax,%rbx,8),%xmm1    ## h2y
        addsd 64(%rax,%rbx,8),%xmm2    ## h2z
        addsd 72(%rax,%rbx,8),%xmm3    ## mx
        addsd 80(%rax,%rbx,8),%xmm4    ## my
        addsd 88(%rax,%rbx,8),%xmm5    ## mz

        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm6,nb213_izH1(%rsp)
        movapd %xmm0,nb213_ixH2(%rsp)
        movapd %xmm1,nb213_iyH2(%rsp)
        movapd %xmm2,nb213_izH2(%rsp)
        movapd %xmm3,nb213_ixM(%rsp)
        movapd %xmm4,nb213_iyM(%rsp)
        movapd %xmm5,nb213_izM(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb213_vctot(%rsp)
        movapd %xmm4,nb213_Vvdwtot(%rsp)
        movapd %xmm4,nb213_fixO(%rsp)
        movapd %xmm4,nb213_fiyO(%rsp)
        movapd %xmm4,nb213_fizO(%rsp)
        movapd %xmm4,nb213_fixH1(%rsp)
        movapd %xmm4,nb213_fiyH1(%rsp)
        movapd %xmm4,nb213_fizH1(%rsp)
        movapd %xmm4,nb213_fixH2(%rsp)
        movapd %xmm4,nb213_fiyH2(%rsp)
        movapd %xmm4,nb213_fizH2(%rsp)
        movapd %xmm4,nb213_fixM(%rsp)
        movapd %xmm4,nb213_fiyM(%rsp)
        movapd %xmm4,nb213_fizM(%rsp)

        movq  nb213_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb213_pos(%rbp),%rsi
        movq  nb213_faction(%rbp),%rdi
        movq  nb213_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb213_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb213_ninner(%rsp),%ecx
        movl  %ecx,nb213_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb213_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel213_x86_64_sse2.nb213_unroll_loop
        jmp   _nb_kernel213_x86_64_sse2.nb213_checksingle
_nb_kernel213_x86_64_sse2.nb213_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb213_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb213_innerjjnr(%rsp)                   ## advance pointer (unrolled 2) 

        movq nb213_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb213_iqM(%rsp),%xmm3
        mulpd  nb213_iqH(%rsp),%xmm4

        movapd  %xmm3,nb213_qqM(%rsp)
        movapd  %xmm4,nb213_qqH(%rsp)

        movq nb213_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movl (%rsi,%rbx,4),%r9d
        movq nb213_vdwparam(%rbp),%rsi
        shll %r8d
        shll %r9d
        movl nb213_ntia(%rsp),%edi
        addl %edi,%r8d
        addl %edi,%r9d

        movlpd (%rsi,%r8,8),%xmm6       ## c6a
        movlpd (%rsi,%r9,8),%xmm7       ## c6b
        movhpd 8(%rsi,%r8,8),%xmm6      ## c6a c12a 
        movhpd 8(%rsi,%r9,8),%xmm7      ## c6b c12b 
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movapd %xmm4,nb213_c6(%rsp)
        movapd %xmm6,nb213_c12(%rsp)

        movq nb213_pos(%rbp),%rsi        ## base of pos[] 

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

    ## O interaction
    ## copy to xmm3-xmm5
    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5

    subpd nb213_ixO(%rsp),%xmm3
    subpd nb213_iyO(%rsp),%xmm4
    subpd nb213_izO(%rsp),%xmm5

    movapd %xmm3,%xmm13
    movapd %xmm4,%xmm14
    movapd %xmm5,%xmm15

        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5

        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3

    ## calc 1/rsq
    cvtpd2ps %xmm3,%xmm6
    rcpps %xmm6,%xmm6
    cvtps2pd %xmm6,%xmm6    ## lu in low xmm6 

    ## 1/x lookup seed in xmm6 
    movapd nb213_two(%rsp),%xmm4
    movapd %xmm3,%xmm5
    mulpd %xmm6,%xmm3       ## lu*rsq 
    subpd %xmm3,%xmm4       ## 2-lu*rsq 
    mulpd %xmm4,%xmm6       ## (new lu) 

    movapd nb213_two(%rsp),%xmm4
    mulpd %xmm6,%xmm5       ## lu*rsq 
    subpd %xmm5,%xmm4       ## 2-lu*rsq 
    mulpd %xmm6,%xmm4       ## xmm4=rinvsq 

    movapd %xmm4,%xmm3      ## rinvsq
    mulpd  %xmm4,%xmm4      ## rinv4
    mulpd  %xmm3,%xmm4      ## rinv6
    movapd %xmm4,%xmm5
    mulpd  %xmm5,%xmm5      ## rinv12
    mulpd  nb213_c6(%rsp),%xmm4
    mulpd  nb213_c12(%rsp),%xmm5
    movapd %xmm5,%xmm6
    subpd  %xmm4,%xmm6 ## Vvdw=vvdw12-vvdw6
    mulpd  nb213_six(%rsp),%xmm4
    mulpd  nb213_twelve(%rsp),%xmm5
    subpd  %xmm4,%xmm5
    mulpd  %xmm5,%xmm3  ## fscal

    addpd  nb213_Vvdwtot(%rsp),%xmm6
    movapd %xmm6,nb213_Vvdwtot(%rsp)

    mulpd  %xmm3,%xmm13 ## fx
    mulpd  %xmm3,%xmm14 ## fy
    mulpd  %xmm3,%xmm15 ## fz

    ## save j force temporarily
    movapd %xmm13,nb213_fjx(%rsp)
    movapd %xmm14,nb213_fjy(%rsp)
    movapd %xmm15,nb213_fjz(%rsp)

    ## increment i O force
    addpd nb213_fixO(%rsp),%xmm13
    addpd nb213_fiyO(%rsp),%xmm14
    addpd nb213_fizO(%rsp),%xmm15
    movapd %xmm13,nb213_fixO(%rsp)
    movapd %xmm14,nb213_fiyO(%rsp)
    movapd %xmm15,nb213_fizO(%rsp)
    ## finished O LJ interaction.


    ## do H1, H2, and M interactions in parallel.
    ## xmm0-xmm2 still contain j coordinates.                
    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb213_ixH1(%rsp),%xmm0
    subpd nb213_iyH1(%rsp),%xmm1
    subpd nb213_izH1(%rsp),%xmm2
    subpd nb213_ixH2(%rsp),%xmm3
    subpd nb213_iyH2(%rsp),%xmm4
    subpd nb213_izH2(%rsp),%xmm5
    subpd nb213_ixM(%rsp),%xmm6
    subpd nb213_iyM(%rsp),%xmm7
    subpd nb213_izM(%rsp),%xmm8

        movapd %xmm0,nb213_dxH1(%rsp)
        movapd %xmm1,nb213_dyH1(%rsp)
        movapd %xmm2,nb213_dzH1(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb213_dxH2(%rsp)
        movapd %xmm4,nb213_dyH2(%rsp)
        movapd %xmm5,nb213_dzH2(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb213_dxM(%rsp)
        movapd %xmm7,nb213_dyM(%rsp)
        movapd %xmm8,nb213_dzM(%rsp)
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

        movapd  nb213_three(%rsp),%xmm9
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

        movapd  nb213_half(%rsp),%xmm15
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

        movapd  nb213_three(%rsp),%xmm1
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

        movapd  nb213_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvH1
        mulpd   %xmm15,%xmm10 ##   rinvH2
    mulpd   %xmm15,%xmm11 ##   rinvM

        ## interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd %xmm9,%xmm1 ## copy of rinv
    movapd %xmm10,%xmm4
    movapd %xmm11,%xmm7
    movapd nb213_krf(%rsp),%xmm2
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
    movapd nb213_crf(%rsp),%xmm14
    subpd  %xmm14,%xmm2  ## rinv+krsq-crf
    subpd  %xmm14,%xmm5
    subpd  %xmm14,%xmm8
    movapd nb213_qqH(%rsp),%xmm12
    movapd nb213_qqM(%rsp),%xmm13
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
    addpd  nb213_vctot(%rsp),%xmm2
    addpd  %xmm8,%xmm5
    addpd  %xmm5,%xmm2
    movapd %xmm2,nb213_vctot(%rsp)

    mulpd  %xmm1,%xmm9  ## fscal
    mulpd  %xmm4,%xmm10
    mulpd  %xmm7,%xmm11

    ## move j forces to xmm0-xmm2
    movq nb213_faction(%rbp),%rdi
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

    ## add forces from O interaction
    addpd nb213_fjx(%rsp),%xmm0
    addpd nb213_fjy(%rsp),%xmm1
    addpd nb213_fjz(%rsp),%xmm2

        mulpd nb213_dxH1(%rsp),%xmm7
        mulpd nb213_dyH1(%rsp),%xmm8
        mulpd nb213_dzH1(%rsp),%xmm9
        mulpd nb213_dxH2(%rsp),%xmm10
        mulpd nb213_dyH2(%rsp),%xmm11
        mulpd nb213_dzH2(%rsp),%xmm12
        mulpd nb213_dxM(%rsp),%xmm13
        mulpd nb213_dyM(%rsp),%xmm14
        mulpd nb213_dzM(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb213_fixH1(%rsp),%xmm7
    addpd nb213_fiyH1(%rsp),%xmm8
    addpd nb213_fizH1(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb213_fixH2(%rsp),%xmm10
    addpd nb213_fiyH2(%rsp),%xmm11
    addpd nb213_fizH2(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb213_fixM(%rsp),%xmm13
    addpd nb213_fiyM(%rsp),%xmm14
    addpd nb213_fizM(%rsp),%xmm15

    movapd %xmm7,nb213_fixH1(%rsp)
    movapd %xmm8,nb213_fiyH1(%rsp)
    movapd %xmm9,nb213_fizH1(%rsp)
    movapd %xmm10,nb213_fixH2(%rsp)
    movapd %xmm11,nb213_fiyH2(%rsp)
    movapd %xmm12,nb213_fizH2(%rsp)
    movapd %xmm13,nb213_fixM(%rsp)
    movapd %xmm14,nb213_fiyM(%rsp)
    movapd %xmm15,nb213_fizM(%rsp)

    ## store back j forces from xmm0-xmm2
        movlpd %xmm0,(%rdi,%rax,8)
        movlpd %xmm1,8(%rdi,%rax,8)
        movlpd %xmm2,16(%rdi,%rax,8)
        movhpd %xmm0,(%rdi,%rbx,8)
        movhpd %xmm1,8(%rdi,%rbx,8)
        movhpd %xmm2,16(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb213_innerk(%rsp)
        jl   _nb_kernel213_x86_64_sse2.nb213_checksingle
        jmp  _nb_kernel213_x86_64_sse2.nb213_unroll_loop
_nb_kernel213_x86_64_sse2.nb213_checksingle: 
        movl  nb213_innerk(%rsp),%edx
        andl  $1,%edx
        jnz  _nb_kernel213_x86_64_sse2.nb213_dosingle
        jmp  _nb_kernel213_x86_64_sse2.nb213_updateouterdata
_nb_kernel213_x86_64_sse2.nb213_dosingle: 
        movq  nb213_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb213_innerjjnr(%rsp)

        movq nb213_charge(%rbp),%rsi     ## base of charge[] 

        xorpd %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
        mulsd  nb213_iqM(%rsp),%xmm3
        mulsd  nb213_iqH(%rsp),%xmm4

        movapd  %xmm3,nb213_qqM(%rsp)
        movapd  %xmm4,nb213_qqH(%rsp)

        movq nb213_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movq nb213_vdwparam(%rbp),%rsi
        shll %r8d
        movl nb213_ntia(%rsp),%edi
        addl %edi,%r8d

        movlpd (%rsi,%r8,8),%xmm6       ## c6a
        movhpd 8(%rsi,%r8,8),%xmm6      ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movapd %xmm4,nb213_c6(%rsp)
        movapd %xmm6,nb213_c12(%rsp)

        movq nb213_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2  and xmm4-xmm6
        movlpd (%rsi,%rax,8),%xmm4
        movlpd 8(%rsi,%rax,8),%xmm5
        movlpd 16(%rsi,%rax,8),%xmm6
    movapd %xmm4,%xmm0
    movapd %xmm5,%xmm1
    movapd %xmm6,%xmm2

        ## calc dr 
        subsd nb213_ixO(%rsp),%xmm4
        subsd nb213_iyO(%rsp),%xmm5
        subsd nb213_izO(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb213_dxO(%rsp)
        movapd %xmm5,nb213_dyO(%rsp)
        movapd %xmm6,nb213_dzO(%rsp)
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
        subsd nb213_ixH1(%rsp),%xmm4
        subsd nb213_iyH1(%rsp),%xmm5
        subsd nb213_izH1(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb213_dxH1(%rsp)
        movapd %xmm5,nb213_dyH1(%rsp)
        movapd %xmm6,nb213_dzH1(%rsp)
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
        subsd nb213_ixH2(%rsp),%xmm3
        subsd nb213_iyH2(%rsp),%xmm4
        subsd nb213_izH2(%rsp),%xmm5

        ## store dr 
        movapd %xmm3,nb213_dxH2(%rsp)
        movapd %xmm4,nb213_dyH2(%rsp)
        movapd %xmm5,nb213_dzH2(%rsp)
        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5

        ## move j coords to xmm4-xmm2
        movapd %xmm0,%xmm4
        movapd %xmm1,%xmm3
    ## xmm2 already contains z

        ## calc dr 
        subsd nb213_ixM(%rsp),%xmm4
        subsd nb213_iyM(%rsp),%xmm3
        subsd nb213_izM(%rsp),%xmm2

        ## store dr 
        movapd %xmm4,nb213_dxM(%rsp)
        movapd %xmm3,nb213_dyM(%rsp)
        movapd %xmm2,nb213_dzM(%rsp)

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## calculate krsq
        movsd nb213_krf(%rsp),%xmm0
        movsd %xmm0,%xmm1
        movsd %xmm0,%xmm2
        mulsd %xmm4,%xmm0
        mulsd %xmm5,%xmm1
        mulsd %xmm6,%xmm2
        movsd %xmm0,nb213_krsqM(%rsp)
        movsd %xmm1,nb213_krsqH2(%rsp)
        movsd %xmm2,nb213_krsqH1(%rsp)

        ## start with rsqH1 - put seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb213_three(%rsp),%xmm1
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb213_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb213_three(%rsp),%xmm1
        subsd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb213_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb213_rinvH1(%rsp)

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb213_three(%rsp),%xmm1
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb213_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb213_three(%rsp),%xmm1
        subsd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb213_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb213_rinvH2(%rsp)

        ## rsqM - seed in xmm2 
        cvtsd2ss %xmm4,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb213_three(%rsp),%xmm1
        mulsd   %xmm4,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb213_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb213_three(%rsp),%xmm1
        subsd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb213_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb213_rinvM(%rsp)

        ## do O interactions directly. xmm7=rsq
        cvtsd2ss %xmm7,%xmm2
        movapd   %xmm7,%xmm6
        rcpps    %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2
        movapd   nb213_two(%rsp),%xmm1
        movapd   %xmm1,%xmm0
        mulsd   %xmm2,%xmm7
        subsd   %xmm7,%xmm1
        mulsd   %xmm1,%xmm2 ## iter1 
        mulsd   %xmm2,%xmm6
        subsd   %xmm6,%xmm0
        mulsd   %xmm2,%xmm0 ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulsd   %xmm1,%xmm1 ## rinv4
        mulsd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulsd   %xmm2,%xmm2 ## rinvtwelve
        mulsd  nb213_c6(%rsp),%xmm1
        mulsd  nb213_c12(%rsp),%xmm2
        movapd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb213_Vvdwtot(%rsp),%xmm3
        mulsd  nb213_six(%rsp),%xmm1
        mulsd  nb213_twelve(%rsp),%xmm2
        subsd  %xmm1,%xmm2
        mulsd  %xmm0,%xmm2
        movapd %xmm2,%xmm4 ## total fsO 
        movsd %xmm3,nb213_Vvdwtot(%rsp)

        movapd nb213_dxO(%rsp),%xmm0
        movapd nb213_dyO(%rsp),%xmm1
        movapd nb213_dzO(%rsp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update O forces 
        movapd nb213_fixO(%rsp),%xmm3
        movapd nb213_fiyO(%rsp),%xmm4
        movapd nb213_fizO(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb213_fixO(%rsp)
        movsd %xmm4,nb213_fiyO(%rsp)
        movsd %xmm7,nb213_fizO(%rsp)
        ## update j forces with water O 
        movsd %xmm0,nb213_fjx(%rsp)
        movsd %xmm1,nb213_fjy(%rsp)
        movsd %xmm2,nb213_fjz(%rsp)

        ## H1 interactions
        movsd  nb213_rinvH1(%rsp),%xmm6
        movsd  %xmm6,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movsd  %xmm6,%xmm7
        movsd  nb213_krsqH1(%rsp),%xmm0
        addsd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        mulsd   nb213_two(%rsp),%xmm0
        subsd   nb213_crf(%rsp),%xmm6
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb213_qqH(%rsp),%xmm6   ## vcoul 
        mulsd   nb213_qqH(%rsp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH1 in xmm4 

        addsd  nb213_vctot(%rsp),%xmm6

        movapd nb213_dxH1(%rsp),%xmm0
        movapd nb213_dyH1(%rsp),%xmm1
        movapd nb213_dzH1(%rsp),%xmm2
        movsd %xmm6,nb213_vctot(%rsp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb213_fixH1(%rsp),%xmm3
        movapd nb213_fiyH1(%rsp),%xmm4
        movapd nb213_fizH1(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb213_fixH1(%rsp)
        movsd %xmm4,nb213_fiyH1(%rsp)
        movsd %xmm7,nb213_fizH1(%rsp)
        ## update j forces with water H1 
        addsd  nb213_fjx(%rsp),%xmm0
        addsd  nb213_fjy(%rsp),%xmm1
        addsd  nb213_fjz(%rsp),%xmm2
        movsd %xmm0,nb213_fjx(%rsp)
        movsd %xmm1,nb213_fjy(%rsp)
        movsd %xmm2,nb213_fjz(%rsp)

        ## H2 interactions 
        movsd  nb213_rinvH2(%rsp),%xmm5
        movsd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movsd  %xmm5,%xmm7
        movsd  nb213_krsqH2(%rsp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        mulsd   nb213_two(%rsp),%xmm0
        subsd   nb213_crf(%rsp),%xmm5
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb213_qqH(%rsp),%xmm5   ## vcoul 
        mulsd   nb213_qqH(%rsp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addsd  nb213_vctot(%rsp),%xmm5

        movapd nb213_dxH2(%rsp),%xmm0
        movapd nb213_dyH2(%rsp),%xmm1
        movapd nb213_dzH2(%rsp),%xmm2
        movsd %xmm5,nb213_vctot(%rsp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb213_fixH2(%rsp),%xmm3
        movapd nb213_fiyH2(%rsp),%xmm4
        movapd nb213_fizH2(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb213_fixH2(%rsp)
        movsd %xmm4,nb213_fiyH2(%rsp)
        movsd %xmm7,nb213_fizH2(%rsp)
        ## update j forces with water H2 
        addsd  nb213_fjx(%rsp),%xmm0
        addsd  nb213_fjy(%rsp),%xmm1
        addsd  nb213_fjz(%rsp),%xmm2
        movsd %xmm0,nb213_fjx(%rsp)
        movsd %xmm1,nb213_fjy(%rsp)
        movsd %xmm2,nb213_fjz(%rsp)

        ## M interactions 
        movsd  nb213_rinvM(%rsp),%xmm5
        movsd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movsd  %xmm5,%xmm7
        movsd  nb213_krsqM(%rsp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        mulsd   nb213_two(%rsp),%xmm0
        subsd   nb213_crf(%rsp),%xmm5
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb213_qqM(%rsp),%xmm5   ## vcoul 
        mulsd   nb213_qqM(%rsp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addsd  nb213_vctot(%rsp),%xmm5

        movapd nb213_dxM(%rsp),%xmm0
        movapd nb213_dyM(%rsp),%xmm1
        movapd nb213_dzM(%rsp),%xmm2
        movsd %xmm5,nb213_vctot(%rsp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update M forces 
        movapd nb213_fixM(%rsp),%xmm3
        movapd nb213_fiyM(%rsp),%xmm4
        movapd nb213_fizM(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb213_fixM(%rsp)
        movsd %xmm4,nb213_fiyM(%rsp)
        movsd %xmm7,nb213_fizM(%rsp)

        movq nb213_faction(%rbp),%rdi
        ## update j forces 
        addsd  nb213_fjx(%rsp),%xmm0
        addsd  nb213_fjy(%rsp),%xmm1
        addsd  nb213_fjz(%rsp),%xmm2
        movlpd (%rdi,%rax,8),%xmm3
        movlpd 8(%rdi,%rax,8),%xmm4
        movlpd 16(%rdi,%rax,8),%xmm5
        addsd %xmm0,%xmm3
        addsd %xmm1,%xmm4
        addsd %xmm2,%xmm5
        movlpd %xmm3,(%rdi,%rax,8)
        movlpd %xmm4,8(%rdi,%rax,8)
        movlpd %xmm5,16(%rdi,%rax,8)

_nb_kernel213_x86_64_sse2.nb213_updateouterdata: 
        movl  nb213_ii3(%rsp),%ecx
        movq  nb213_faction(%rbp),%rdi
        movq  nb213_fshift(%rbp),%rsi
        movl  nb213_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb213_fixO(%rsp),%xmm0
        movapd nb213_fiyO(%rsp),%xmm1
        movapd nb213_fizO(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

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
        movapd nb213_fixH1(%rsp),%xmm0
        movapd nb213_fiyH1(%rsp),%xmm1
        movapd nb213_fizH1(%rsp),%xmm2

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
        movapd nb213_fixH2(%rsp),%xmm0
        movapd nb213_fiyH2(%rsp),%xmm1
        movapd nb213_fizH2(%rsp),%xmm2

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

        ## accumulate Mi forces in xmm0, xmm1, xmm2 
        movapd nb213_fixM(%rsp),%xmm0
        movapd nb213_fiyM(%rsp),%xmm1
        movapd nb213_fizM(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

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
        movl nb213_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb213_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb213_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb213_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb213_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb213_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

       ## finish if last 
        movl nb213_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel213_x86_64_sse2.nb213_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb213_n(%rsp)
        jmp _nb_kernel213_x86_64_sse2.nb213_outer
_nb_kernel213_x86_64_sse2.nb213_outerend: 
        ## check if more outer neighborlists remain
        movl  nb213_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel213_x86_64_sse2.nb213_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel213_x86_64_sse2.nb213_threadloop
_nb_kernel213_x86_64_sse2.nb213_end: 
        movl nb213_nouter(%rsp),%eax
        movl nb213_ninner(%rsp),%ebx
        movq nb213_outeriter(%rbp),%rcx
        movq nb213_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1064,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret





.globl nb_kernel213nf_x86_64_sse2
.globl _nb_kernel213nf_x86_64_sse2
nb_kernel213nf_x86_64_sse2:     
_nb_kernel213nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb213nf_fshift, 16
.set nb213nf_gid, 24
.set nb213nf_pos, 32
.set nb213nf_faction, 40
.set nb213nf_charge, 48
.set nb213nf_p_facel, 56
.set nb213nf_argkrf, 64
.set nb213nf_argcrf, 72
.set nb213nf_Vc, 80
.set nb213nf_type, 88
.set nb213nf_p_ntype, 96
.set nb213nf_vdwparam, 104
.set nb213nf_Vvdw, 112
.set nb213nf_p_tabscale, 120
.set nb213nf_VFtab, 128
.set nb213nf_invsqrta, 136
.set nb213nf_dvda, 144
.set nb213nf_p_gbtabscale, 152
.set nb213nf_GBtab, 160
.set nb213nf_p_nthreads, 168
.set nb213nf_count, 176
.set nb213nf_mtx, 184
.set nb213nf_outeriter, 192
.set nb213nf_inneriter, 200
.set nb213nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb213nf_ixO, 0
.set nb213nf_iyO, 16
.set nb213nf_izO, 32
.set nb213nf_ixH1, 48
.set nb213nf_iyH1, 64
.set nb213nf_izH1, 80
.set nb213nf_ixH2, 96
.set nb213nf_iyH2, 112
.set nb213nf_izH2, 128
.set nb213nf_ixM, 144
.set nb213nf_iyM, 160
.set nb213nf_izM, 176
.set nb213nf_iqH, 192
.set nb213nf_iqM, 208
.set nb213nf_qqH, 224
.set nb213nf_qqM, 240
.set nb213nf_c6, 256
.set nb213nf_c12, 272
.set nb213nf_vctot, 288
.set nb213nf_Vvdwtot, 304
.set nb213nf_half, 320
.set nb213nf_three, 336
.set nb213nf_two, 352
.set nb213nf_rinvH1, 368
.set nb213nf_rinvH2, 384
.set nb213nf_rinvM, 400
.set nb213nf_krsqH1, 416
.set nb213nf_krsqH2, 432
.set nb213nf_krsqM, 448
.set nb213nf_krf, 464
.set nb213nf_crf, 480
.set nb213nf_nri, 496
.set nb213nf_iinr, 504
.set nb213nf_jindex, 512
.set nb213nf_jjnr, 520
.set nb213nf_shift, 528
.set nb213nf_shiftvec, 536
.set nb213nf_facel, 544
.set nb213nf_innerjjnr, 552
.set nb213nf_is3, 560
.set nb213nf_ii3, 564
.set nb213nf_ntia, 568
.set nb213nf_innerk, 572
.set nb213nf_n, 576
.set nb213nf_nn1, 580
.set nb213nf_nouter, 584
.set nb213nf_ninner, 588
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms


        push %r12
        push %r13
        push %r14
        push %r15

        subq $600,%rsp          ## local variable stack space (n*16+8)
        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb213nf_nouter(%rsp)
        movl %eax,nb213nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb213nf_nri(%rsp)
        movq %rsi,nb213nf_iinr(%rsp)
        movq %rdx,nb213nf_jindex(%rsp)
        movq %rcx,nb213nf_jjnr(%rsp)
        movq %r8,nb213nf_shift(%rsp)
        movq %r9,nb213nf_shiftvec(%rsp)
        movq nb213nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb213nf_facel(%rsp)

        movq nb213nf_argkrf(%rbp),%rsi
        movq nb213nf_argcrf(%rbp),%rdi
        movsd (%rsi),%xmm1
        movsd (%rdi),%xmm2
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        movapd %xmm1,nb213nf_krf(%rsp)
        movapd %xmm2,nb213nf_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb213nf_half(%rsp)
        movl %ebx,nb213nf_half+4(%rsp)
        movsd nb213nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb213nf_half(%rsp)
        movapd %xmm2,nb213nf_two(%rsp)
        movapd %xmm3,nb213nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb213nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb213nf_charge(%rbp),%rdx
        movsd 8(%rdx,%rbx,8),%xmm3
        movsd 24(%rdx,%rbx,8),%xmm4

        movsd nb213nf_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb213nf_iqH(%rsp)
        movapd %xmm4,nb213nf_iqM(%rsp)

        movq  nb213nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb213nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb213nf_ntia(%rsp)
_nb_kernel213nf_x86_64_sse2.nb213nf_threadloop: 
        movq  nb213nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel213nf_x86_64_sse2.nb213nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel213nf_x86_64_sse2.nb213nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb213nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb213nf_n(%rsp)
        movl %ebx,nb213nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel213nf_x86_64_sse2.nb213nf_outerstart
        jmp _nb_kernel213nf_x86_64_sse2.nb213nf_end

_nb_kernel213nf_x86_64_sse2.nb213nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb213nf_nouter(%rsp),%ebx
        movl %ebx,nb213nf_nouter(%rsp)

_nb_kernel213nf_x86_64_sse2.nb213nf_outer: 
        movq  nb213nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb213nf_is3(%rsp)            ## store is3 

        movq  nb213nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb213nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb213nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb213nf_ii3(%rsp)

        addsd (%rax,%rbx,8),%xmm3       ## ox
        addsd 8(%rax,%rbx,8),%xmm4      ## oy
        addsd 16(%rax,%rbx,8),%xmm5     ## oz   
        addsd 24(%rax,%rbx,8),%xmm6     ## h1x
        addsd 32(%rax,%rbx,8),%xmm7     ## h1y
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm7,%xmm7
        movapd %xmm3,nb213nf_ixO(%rsp)
        movapd %xmm4,nb213nf_iyO(%rsp)
        movapd %xmm5,nb213nf_izO(%rsp)
        movapd %xmm6,nb213nf_ixH1(%rsp)
        movapd %xmm7,nb213nf_iyH1(%rsp)

        movsd %xmm2,%xmm6
        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 40(%rax,%rbx,8),%xmm6    ## h1z
        addsd 48(%rax,%rbx,8),%xmm0    ## h2x
        addsd 56(%rax,%rbx,8),%xmm1    ## h2y
        addsd 64(%rax,%rbx,8),%xmm2    ## h2z
        addsd 72(%rax,%rbx,8),%xmm3    ## mx
        addsd 80(%rax,%rbx,8),%xmm4    ## my
        addsd 88(%rax,%rbx,8),%xmm5    ## mz

        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm6,nb213nf_izH1(%rsp)
        movapd %xmm0,nb213nf_ixH2(%rsp)
        movapd %xmm1,nb213nf_iyH2(%rsp)
        movapd %xmm2,nb213nf_izH2(%rsp)
        movapd %xmm3,nb213nf_ixM(%rsp)
        movapd %xmm4,nb213nf_iyM(%rsp)
        movapd %xmm5,nb213nf_izM(%rsp)

        ## clear vctot
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb213nf_vctot(%rsp)
        movapd %xmm4,nb213nf_Vvdwtot(%rsp)

        movq  nb213nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb213nf_pos(%rbp),%rsi
        movq  nb213nf_faction(%rbp),%rdi
        movq  nb213nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb213nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb213nf_ninner(%rsp),%ecx
        movl  %ecx,nb213nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb213nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel213nf_x86_64_sse2.nb213nf_unroll_loop
        jmp   _nb_kernel213nf_x86_64_sse2.nb213nf_checksingle
_nb_kernel213nf_x86_64_sse2.nb213nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb213nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb213nf_innerjjnr(%rsp)
        ## advance pointer (unrolled 2) 

        movq nb213nf_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb213nf_iqM(%rsp),%xmm3
        mulpd  nb213nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movapd  %xmm3,nb213nf_qqM(%rsp)
        movapd  %xmm4,nb213nf_qqH(%rsp)

        movq nb213nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movq nb213nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        movl nb213nf_ntia(%rsp),%edi
        addl %edi,%eax
        addl %edi,%ebx

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movlpd (%rsi,%rbx,8),%xmm7      ## c6b
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 
        movhpd 8(%rsi,%rbx,8),%xmm7     ## c6b c12b 
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb213nf_c6(%rsp)
        movapd %xmm6,nb213nf_c12(%rsp)

        movq nb213nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movapd nb213nf_ixO(%rsp),%xmm4
        movapd nb213nf_iyO(%rsp),%xmm5
        movapd nb213nf_izO(%rsp),%xmm6

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
        movapd nb213nf_ixH1(%rsp),%xmm4
        movapd nb213nf_iyH1(%rsp),%xmm5
        movapd nb213nf_izH1(%rsp),%xmm6

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
        movapd nb213nf_ixH2(%rsp),%xmm3
        movapd nb213nf_iyH2(%rsp),%xmm4
        movapd nb213nf_izH2(%rsp),%xmm5

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

        ## move ixM-izM to xmm2-xmm4  
        movapd nb213nf_iyM(%rsp),%xmm3
        movapd nb213nf_izM(%rsp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb213nf_ixM(%rsp),%xmm2
        subpd  %xmm0,%xmm2

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## calculate krsq
        movapd nb213nf_krf(%rsp),%xmm0
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        mulpd %xmm4,%xmm0
        mulpd %xmm5,%xmm1
        mulpd %xmm6,%xmm2
        movapd %xmm0,nb213nf_krsqM(%rsp)
        movapd %xmm1,nb213nf_krsqH2(%rsp)
        movapd %xmm2,nb213nf_krsqH1(%rsp)

        ## start with rsqH1 - put seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb213nf_three(%rsp),%xmm1
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb213nf_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb213nf_three(%rsp),%xmm1
        subpd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb213nf_half(%rsp),%xmm1   ## rinv 
        movapd  %xmm1,nb213nf_rinvH1(%rsp)

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb213nf_three(%rsp),%xmm1
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb213nf_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb213nf_three(%rsp),%xmm1
        subpd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb213nf_half(%rsp),%xmm1   ## rinv 
        movapd  %xmm1,nb213nf_rinvH2(%rsp)

        ## rsqM - seed in xmm2 
        cvtpd2ps %xmm4,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb213nf_three(%rsp),%xmm1
        mulpd   %xmm4,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb213nf_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb213nf_three(%rsp),%xmm1
        subpd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb213nf_half(%rsp),%xmm1   ## rinv 
        movapd  %xmm1,nb213nf_rinvM(%rsp)

        ## do O interactions directly - rsqO is in xmm7
        cvtpd2ps %xmm7,%xmm2
        movapd   %xmm7,%xmm6
        rcpps    %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2
        movapd   nb213nf_two(%rsp),%xmm1
        movapd   %xmm1,%xmm0
        mulpd   %xmm2,%xmm7
        subpd   %xmm7,%xmm1
        mulpd   %xmm1,%xmm2 ## iter1 
        mulpd   %xmm2,%xmm6
        subpd   %xmm6,%xmm0
        mulpd   %xmm2,%xmm0 ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulpd   %xmm1,%xmm1 ## rinv4
        mulpd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulpd   %xmm2,%xmm2 ## rinvtwelve
        mulpd  nb213nf_c6(%rsp),%xmm1
        mulpd  nb213nf_c12(%rsp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addpd  nb213nf_Vvdwtot(%rsp),%xmm3
        movapd %xmm3,nb213nf_Vvdwtot(%rsp)

        ## H1 interactions 
        movapd  nb213nf_rinvH1(%rsp),%xmm6
        addpd   nb213nf_krsqH1(%rsp),%xmm6
        subpd   nb213nf_crf(%rsp),%xmm6
        mulpd   nb213nf_qqH(%rsp),%xmm6   ## vcoul 
        addpd   nb213nf_vctot(%rsp),%xmm6
        movapd %xmm6,nb213nf_vctot(%rsp)

        ## H2 interactions 
        movapd  nb213nf_rinvH2(%rsp),%xmm6
        addpd   nb213nf_krsqH2(%rsp),%xmm6
        subpd   nb213nf_crf(%rsp),%xmm6
        mulpd   nb213nf_qqH(%rsp),%xmm6   ## vcoul 
        addpd   nb213nf_vctot(%rsp),%xmm6
        movapd %xmm6,nb213nf_vctot(%rsp)

        ## M interactions 
        movapd  nb213nf_rinvM(%rsp),%xmm6
        addpd   nb213nf_krsqM(%rsp),%xmm6
        subpd   nb213nf_crf(%rsp),%xmm6
        mulpd   nb213nf_qqM(%rsp),%xmm6   ## vcoul 
        addpd   nb213nf_vctot(%rsp),%xmm6
        movapd %xmm6,nb213nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb213nf_innerk(%rsp)
        jl   _nb_kernel213nf_x86_64_sse2.nb213nf_checksingle
        jmp  _nb_kernel213nf_x86_64_sse2.nb213nf_unroll_loop
_nb_kernel213nf_x86_64_sse2.nb213nf_checksingle: 
        movl  nb213nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz  _nb_kernel213nf_x86_64_sse2.nb213nf_dosingle
        jmp  _nb_kernel213nf_x86_64_sse2.nb213nf_updateouterdata
_nb_kernel213nf_x86_64_sse2.nb213nf_dosingle: 
        movq  nb213nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb213nf_innerjjnr(%rsp)

        movq nb213nf_charge(%rbp),%rsi     ## base of charge[] 

        xorpd %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
        mulsd  nb213nf_iqM(%rsp),%xmm3
        mulsd  nb213nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movapd  %xmm3,nb213nf_qqM(%rsp)
        movapd  %xmm4,nb213nf_qqH(%rsp)

        movq nb213nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movq nb213nf_vdwparam(%rbp),%rsi
        shll %eax
        movl nb213nf_ntia(%rsp),%edi
        addl %edi,%eax

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb213nf_c6(%rsp)
        movapd %xmm6,nb213nf_c12(%rsp)

        movq nb213nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb213nf_ixO(%rsp),%xmm4
        movapd nb213nf_iyO(%rsp),%xmm5
        movapd nb213nf_izO(%rsp),%xmm6

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
        movapd nb213nf_ixH1(%rsp),%xmm4
        movapd nb213nf_iyH1(%rsp),%xmm5
        movapd nb213nf_izH1(%rsp),%xmm6

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
        movapd nb213nf_ixH2(%rsp),%xmm3
        movapd nb213nf_iyH2(%rsp),%xmm4
        movapd nb213nf_izH2(%rsp),%xmm5

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
        ## move ixM-izM to xmm2-xmm4  
        movapd nb213nf_iyM(%rsp),%xmm3
        movapd nb213nf_izM(%rsp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb213nf_ixM(%rsp),%xmm2
        subpd  %xmm0,%xmm2

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## calculate krsq
        movsd nb213nf_krf(%rsp),%xmm0
        movsd %xmm0,%xmm1
        movsd %xmm0,%xmm2
        mulsd %xmm4,%xmm0
        mulsd %xmm5,%xmm1
        mulsd %xmm6,%xmm2
        movsd %xmm0,nb213nf_krsqM(%rsp)
        movsd %xmm1,nb213nf_krsqH2(%rsp)
        movsd %xmm2,nb213nf_krsqH1(%rsp)

        ## start with rsqH1 - put seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb213nf_three(%rsp),%xmm1
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb213nf_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb213nf_three(%rsp),%xmm1
        subsd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb213nf_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb213nf_rinvH1(%rsp)

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb213nf_three(%rsp),%xmm1
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb213nf_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb213nf_three(%rsp),%xmm1
        subsd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb213nf_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb213nf_rinvH2(%rsp)

        ## rsqM - seed in xmm2 
        cvtsd2ss %xmm4,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb213nf_three(%rsp),%xmm1
        mulsd   %xmm4,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb213nf_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb213nf_three(%rsp),%xmm1
        subsd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb213nf_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb213nf_rinvM(%rsp)

        ## do O interactions directly. xmm7=rsq
        cvtsd2ss %xmm7,%xmm2
        movapd   %xmm7,%xmm6
        rcpps    %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2
        movapd   nb213nf_two(%rsp),%xmm1
        movapd   %xmm1,%xmm0
        mulsd   %xmm2,%xmm7
        subsd   %xmm7,%xmm1
        mulsd   %xmm1,%xmm2 ## iter1 
        mulsd   %xmm2,%xmm6
        subsd   %xmm6,%xmm0
        mulsd   %xmm2,%xmm0 ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulsd   %xmm1,%xmm1 ## rinv4
        mulsd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulsd   %xmm2,%xmm2 ## rinvtwelve
        mulsd  nb213nf_c6(%rsp),%xmm1
        mulsd  nb213nf_c12(%rsp),%xmm2
        movapd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb213nf_Vvdwtot(%rsp),%xmm3
        movsd %xmm3,nb213nf_Vvdwtot(%rsp)

        ## H1 interactions 
        movsd  nb213nf_rinvH1(%rsp),%xmm6
        addsd   nb213nf_krsqH1(%rsp),%xmm6
        subsd   nb213nf_crf(%rsp),%xmm6
        mulsd   nb213nf_qqH(%rsp),%xmm6   ## vcoul 
        addsd   nb213nf_vctot(%rsp),%xmm6
        movsd %xmm6,nb213nf_vctot(%rsp)

        ## H2 interactions 
        movsd  nb213nf_rinvH2(%rsp),%xmm6
        addsd   nb213nf_krsqH2(%rsp),%xmm6
        subsd   nb213nf_crf(%rsp),%xmm6
        mulsd   nb213nf_qqH(%rsp),%xmm6   ## vcoul 
        addsd   nb213nf_vctot(%rsp),%xmm6
        movsd %xmm6,nb213nf_vctot(%rsp)

        ## M interactions 
        movsd  nb213nf_rinvM(%rsp),%xmm6
        addsd   nb213nf_krsqM(%rsp),%xmm6
        subsd   nb213nf_crf(%rsp),%xmm6
        mulsd   nb213nf_qqM(%rsp),%xmm6   ## vcoul 
        addsd   nb213nf_vctot(%rsp),%xmm6
        movsd %xmm6,nb213nf_vctot(%rsp)

_nb_kernel213nf_x86_64_sse2.nb213nf_updateouterdata: 
        ## get n from stack
        movl nb213nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb213nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb213nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb213nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb213nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb213nf_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

       ## finish if last 
        movl nb213nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel213nf_x86_64_sse2.nb213nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb213nf_n(%rsp)
        jmp _nb_kernel213nf_x86_64_sse2.nb213nf_outer
_nb_kernel213nf_x86_64_sse2.nb213nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb213nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel213nf_x86_64_sse2.nb213nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel213nf_x86_64_sse2.nb213nf_threadloop
_nb_kernel213nf_x86_64_sse2.nb213nf_end: 
        movl nb213nf_nouter(%rsp),%eax
        movl nb213nf_ninner(%rsp),%ebx
        movq nb213nf_outeriter(%rbp),%rcx
        movq nb213nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $600,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret


