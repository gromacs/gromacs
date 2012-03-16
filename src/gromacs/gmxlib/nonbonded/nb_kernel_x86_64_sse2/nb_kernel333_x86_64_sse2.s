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








.globl nb_kernel333_x86_64_sse2
.globl _nb_kernel333_x86_64_sse2
nb_kernel333_x86_64_sse2:       
_nb_kernel333_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb333_fshift, 16
.set nb333_gid, 24
.set nb333_pos, 32
.set nb333_faction, 40
.set nb333_charge, 48
.set nb333_p_facel, 56
.set nb333_argkrf, 64
.set nb333_argcrf, 72
.set nb333_Vc, 80
.set nb333_type, 88
.set nb333_p_ntype, 96
.set nb333_vdwparam, 104
.set nb333_Vvdw, 112
.set nb333_p_tabscale, 120
.set nb333_VFtab, 128
.set nb333_invsqrta, 136
.set nb333_dvda, 144
.set nb333_p_gbtabscale, 152
.set nb333_GBtab, 160
.set nb333_p_nthreads, 168
.set nb333_count, 176
.set nb333_mtx, 184
.set nb333_outeriter, 192
.set nb333_inneriter, 200
.set nb333_work, 208
        ## stack offsets for local variables 
        ## bottom of stack is cache-aligned for sse2 use 
.set nb333_ixO, 0
.set nb333_iyO, 16
.set nb333_izO, 32
.set nb333_ixH1, 48
.set nb333_iyH1, 64
.set nb333_izH1, 80
.set nb333_ixH2, 96
.set nb333_iyH2, 112
.set nb333_izH2, 128
.set nb333_ixM, 144
.set nb333_iyM, 160
.set nb333_izM, 176
.set nb333_iqM, 192
.set nb333_iqH, 208
.set nb333_dxO, 224
.set nb333_dyO, 240
.set nb333_dzO, 256
.set nb333_dxH1, 272
.set nb333_dyH1, 288
.set nb333_dzH1, 304
.set nb333_dxH2, 320
.set nb333_dyH2, 336
.set nb333_dzH2, 352
.set nb333_dxM, 368
.set nb333_dyM, 384
.set nb333_dzM, 400
.set nb333_qqM, 416
.set nb333_qqH, 432
.set nb333_rinvO, 448
.set nb333_rinvH1, 464
.set nb333_rinvH2, 480
.set nb333_rinvM, 496
.set nb333_rO, 512
.set nb333_rH1, 528
.set nb333_rH2, 544
.set nb333_rM, 560
.set nb333_tsc, 576
.set nb333_two, 592
.set nb333_c6, 608
.set nb333_c12, 624
.set nb333_vctot, 640
.set nb333_Vvdwtot, 656
.set nb333_fixO, 672
.set nb333_fiyO, 688
.set nb333_fizO, 704
.set nb333_fixH1, 720
.set nb333_fiyH1, 736
.set nb333_fizH1, 752
.set nb333_fixH2, 768
.set nb333_fiyH2, 784
.set nb333_fizH2, 800
.set nb333_fixM, 816
.set nb333_fiyM, 832
.set nb333_fizM, 848
.set nb333_fjx, 864
.set nb333_fjy, 880
.set nb333_fjz, 896
.set nb333_half, 912
.set nb333_three, 928
.set nb333_is3, 944
.set nb333_ii3, 948
.set nb333_nri, 952
.set nb333_iinr, 960
.set nb333_jindex, 968
.set nb333_jjnr, 976
.set nb333_shift, 984
.set nb333_shiftvec, 992
.set nb333_facel, 1000
.set nb333_innerjjnr, 1008
.set nb333_ntia, 1016
.set nb333_innerk, 1020
.set nb333_n, 1024
.set nb333_nn1, 1028
.set nb333_nouter, 1032
.set nb333_ninner, 1036
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1048,%rsp         ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb333_nouter(%rsp)
        movl %eax,nb333_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb333_nri(%rsp)
        movq %rsi,nb333_iinr(%rsp)
        movq %rdx,nb333_jindex(%rsp)
        movq %rcx,nb333_jjnr(%rsp)
        movq %r8,nb333_shift(%rsp)
        movq %r9,nb333_shiftvec(%rsp)
        movq nb333_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb333_facel(%rsp)

        movq nb333_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb333_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb333_half(%rsp)
        movl %ebx,nb333_half+4(%rsp)
        movsd nb333_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb333_half(%rsp)
        movapd %xmm2,nb333_two(%rsp)
        movapd %xmm3,nb333_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb333_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb333_charge(%rbp),%rdx
        movsd 8(%rdx,%rbx,8),%xmm3
        movsd 24(%rdx,%rbx,8),%xmm4
        movq nb333_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb333_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb333_iqH(%rsp)
        movapd %xmm4,nb333_iqM(%rsp)

        movq  nb333_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb333_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb333_ntia(%rsp)
_nb_kernel333_x86_64_sse2.nb333_threadloop: 
        movq  nb333_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel333_x86_64_sse2.nb333_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel333_x86_64_sse2.nb333_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb333_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb333_n(%rsp)
        movl %ebx,nb333_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel333_x86_64_sse2.nb333_outerstart
        jmp _nb_kernel333_x86_64_sse2.nb333_end

_nb_kernel333_x86_64_sse2.nb333_outerstart: 
        ## ebx contains number of outer iterations
        addl nb333_nouter(%rsp),%ebx
        movl %ebx,nb333_nouter(%rsp)

_nb_kernel333_x86_64_sse2.nb333_outer: 
        movq  nb333_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb333_is3(%rsp)      ## store is3 

        movq  nb333_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb333_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb333_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb333_ii3(%rsp)

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
        movapd %xmm3,nb333_ixO(%rsp)
        movapd %xmm4,nb333_iyO(%rsp)
        movapd %xmm5,nb333_izO(%rsp)
        movapd %xmm6,nb333_ixH1(%rsp)
        movapd %xmm7,nb333_iyH1(%rsp)

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
        movapd %xmm6,nb333_izH1(%rsp)
        movapd %xmm0,nb333_ixH2(%rsp)
        movapd %xmm1,nb333_iyH2(%rsp)
        movapd %xmm2,nb333_izH2(%rsp)
        movapd %xmm3,nb333_ixM(%rsp)
        movapd %xmm4,nb333_iyM(%rsp)
        movapd %xmm5,nb333_izM(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb333_vctot(%rsp)
        movapd %xmm4,nb333_Vvdwtot(%rsp)
        movapd %xmm4,nb333_fixO(%rsp)
        movapd %xmm4,nb333_fiyO(%rsp)
        movapd %xmm4,nb333_fizO(%rsp)
        movapd %xmm4,nb333_fixH1(%rsp)
        movapd %xmm4,nb333_fiyH1(%rsp)
        movapd %xmm4,nb333_fizH1(%rsp)
        movapd %xmm4,nb333_fixH2(%rsp)
        movapd %xmm4,nb333_fiyH2(%rsp)
        movapd %xmm4,nb333_fizH2(%rsp)
        movapd %xmm4,nb333_fixM(%rsp)
        movapd %xmm4,nb333_fiyM(%rsp)
        movapd %xmm4,nb333_fizM(%rsp)

        movq  nb333_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx     ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb333_pos(%rbp),%rsi
        movq  nb333_faction(%rbp),%rdi
        movq  nb333_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb333_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb333_ninner(%rsp),%ecx
        movl  %ecx,nb333_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb333_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel333_x86_64_sse2.nb333_unroll_loop
        jmp   _nb_kernel333_x86_64_sse2.nb333_checksingle
_nb_kernel333_x86_64_sse2.nb333_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb333_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb333_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb333_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb333_iqM(%rsp),%xmm3
        mulpd  nb333_iqH(%rsp),%xmm4
        movapd  %xmm3,nb333_qqM(%rsp)
        movapd  %xmm4,nb333_qqH(%rsp)

        movq nb333_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movl (%rsi,%rbx,4),%r9d
        movq nb333_vdwparam(%rbp),%rsi
        shll %r8d
        shll %r9d
        movl nb333_ntia(%rsp),%edi
        addl %edi,%r8d
        addl %edi,%r9d

        movlpd (%rsi,%r8,8),%xmm6       ## c6a
        movlpd (%rsi,%r9,8),%xmm7       ## c6b
        movhpd 8(%rsi,%r8,8),%xmm6      ## c6a c12a 
        movhpd 8(%rsi,%r9,8),%xmm7      ## c6b c12b 
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movapd %xmm4,nb333_c6(%rsp)
        movapd %xmm6,nb333_c12(%rsp)

        movq nb333_pos(%rbp),%rsi        ## base of pos[] 

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

    ## O interaction
    ## copy to xmm3-xmm5
    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5

    subpd nb333_ixO(%rsp),%xmm3
    subpd nb333_iyO(%rsp),%xmm4
    subpd nb333_izO(%rsp),%xmm5

    movapd %xmm3,nb333_dxO(%rsp)
    movapd %xmm4,nb333_dyO(%rsp)
    movapd %xmm5,nb333_dzO(%rsp)

        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5

        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
    ## xmm3=rsq

    cvtpd2ps %xmm3,%xmm5
    rsqrtps %xmm5,%xmm5
    cvtps2pd %xmm5,%xmm15    ## lu in low xmm2 

    ## lookup seed in xmm2 
    movapd %xmm15,%xmm5      ## copy of lu 
    mulpd %xmm15,%xmm15       ## lu*lu 
    movapd nb333_three(%rsp),%xmm7
    mulpd %xmm3,%xmm15       ## rsq*lu*lu                    
    movapd nb333_half(%rsp),%xmm6
    subpd %xmm15,%xmm7       ## 30-rsq*lu*lu 
    mulpd %xmm5,%xmm7
    mulpd %xmm6,%xmm7       ## xmm0=iter1 of rinv (new lu) 

    movapd %xmm7,%xmm5      ## copy of lu 
    mulpd %xmm7,%xmm7       ## lu*lu 
    movapd nb333_three(%rsp),%xmm15
    mulpd %xmm3,%xmm7       ## rsq*lu*lu                    
    movapd nb333_half(%rsp),%xmm6
    subpd %xmm7,%xmm15       ## 30-rsq*lu*lu 
    mulpd %xmm5,%xmm15
    mulpd %xmm6,%xmm15       ## xmm15=rinv

    mulpd %xmm15,%xmm3       ## xmm3=r 

    ## xmm15=rinv
    ## xmm3=r

    mulpd nb333_tsc(%rsp),%xmm3   ## rtab

    ## truncate and convert to integers
    cvttpd2pi %xmm3,%mm6

    ## convert back to float
    cvtpi2pd  %mm6,%xmm4

    ## multiply by 4
    pslld   $2,%mm6

    ## calculate eps
    subpd     %xmm4,%xmm3   ## xmm3=eps

    ## move to integer registers
    movd %mm6,%r10d
    psrlq $32,%mm6
    movd %mm6,%r11d

    ## multiply by 3
    lea  (%r10,%r10,2),%r10
        lea  (%r11,%r11,2),%r11

    ## xmm3=eps
    ## xmm15=rinv

        movq nb333_VFtab(%rbp),%rsi
    ## indices in r10, r11. Load dispersion and repulsion tables in parallel.
    movapd 32(%rsi,%r10,8),%xmm4        ## Y1d F1d  
    movapd 32(%rsi,%r11,8),%xmm12       ## Y2d F2d 
    movapd 64(%rsi,%r10,8),%xmm8        ## Y1r F1r  
    movapd 64(%rsi,%r11,8),%xmm13       ## Y2r F2r 
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm9
    unpcklpd %xmm12,%xmm4   ## Y1d Y2d 
    unpckhpd %xmm12,%xmm5   ## F1d F2d 
    unpcklpd %xmm13,%xmm8   ## Y1r Y2r 
    unpckhpd %xmm13,%xmm9   ## F1r F2r 

    movapd 48(%rsi,%r10,8),%xmm6        ## G1d H1d  
    movapd 48(%rsi,%r11,8),%xmm12       ## G2d H2d 
    movapd 80(%rsi,%r10,8),%xmm10       ## G1r H1r      
    movapd 80(%rsi,%r11,8),%xmm13       ## G2r H2r 
    movapd %xmm6,%xmm7
    movapd %xmm10,%xmm11
    unpcklpd %xmm12,%xmm6   ## G1d G2d 
    unpckhpd %xmm12,%xmm7   ## H1d H2d 
    unpcklpd %xmm13,%xmm10  ## G1r G2r 
    unpckhpd %xmm13,%xmm11  ## H1r H2r 
    ## dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11

    mulpd  %xmm3,%xmm7   ## Heps
    mulpd  %xmm3,%xmm11
    mulpd  %xmm3,%xmm6  ## Geps
    mulpd  %xmm3,%xmm10
    mulpd  %xmm3,%xmm7  ## Heps2
    mulpd  %xmm3,%xmm11
    addpd  %xmm6,%xmm5 ## F+Geps
    addpd  %xmm10,%xmm9
    addpd  %xmm7,%xmm5  ## F+Geps+Heps2 = Fp
    addpd  %xmm11,%xmm9
    addpd  %xmm7,%xmm7   ## 2*Heps2
    addpd  %xmm11,%xmm11
    addpd  %xmm6,%xmm7  ## 2*Heps2+Geps
    addpd  %xmm10,%xmm11

    addpd  %xmm5,%xmm7 ## FF = Fp + 2*Heps2 + Geps
    addpd  %xmm9,%xmm11
    mulpd  %xmm3,%xmm5 ## eps*Fp
    mulpd  %xmm3,%xmm9
    movapd nb333_c6(%rsp),%xmm12
    movapd nb333_c12(%rsp),%xmm13
    addpd  %xmm4,%xmm5 ## VV
    addpd  %xmm8,%xmm9

    mulpd  %xmm12,%xmm5 ## VV*c6 = vnb6
    mulpd  %xmm13,%xmm9 ## VV*c12 = vnb12
    addpd  %xmm9,%xmm5
    addpd  nb333_Vvdwtot(%rsp),%xmm5
    movapd %xmm5,nb333_Vvdwtot(%rsp)

    mulpd  %xmm12,%xmm7  ## FF*c6 = fnb6
    mulpd  %xmm13,%xmm11  ## FF*c12  = fnb12
    addpd  %xmm11,%xmm7

    mulpd  nb333_tsc(%rsp),%xmm7
    mulpd  %xmm15,%xmm7  ## -fscal
    xorpd  %xmm9,%xmm9

    subpd  %xmm7,%xmm9    ## fscal
    movapd %xmm9,%xmm10
    movapd %xmm9,%xmm11

    mulpd  nb333_dxO(%rsp),%xmm9    ## fx/fy/fz
    mulpd  nb333_dyO(%rsp),%xmm10
    mulpd  nb333_dzO(%rsp),%xmm11

    ## save j force temporarily
    movapd %xmm9,nb333_fjx(%rsp)
    movapd %xmm10,nb333_fjy(%rsp)
    movapd %xmm11,nb333_fjz(%rsp)

    ## increment i O force
    addpd nb333_fixO(%rsp),%xmm9
    addpd nb333_fiyO(%rsp),%xmm10
    addpd nb333_fizO(%rsp),%xmm11
    movapd %xmm9,nb333_fixO(%rsp)
    movapd %xmm10,nb333_fiyO(%rsp)
    movapd %xmm11,nb333_fizO(%rsp)
    ## finished O LJ interaction.


    ## do H1, H2, and M interactions in parallel.
    ## xmm0-xmm2 still contain j coordinates.                
    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb333_ixH1(%rsp),%xmm0
    subpd nb333_iyH1(%rsp),%xmm1
    subpd nb333_izH1(%rsp),%xmm2
    subpd nb333_ixH2(%rsp),%xmm3
    subpd nb333_iyH2(%rsp),%xmm4
    subpd nb333_izH2(%rsp),%xmm5
    subpd nb333_ixM(%rsp),%xmm6
    subpd nb333_iyM(%rsp),%xmm7
    subpd nb333_izM(%rsp),%xmm8

        movapd %xmm0,nb333_dxH1(%rsp)
        movapd %xmm1,nb333_dyH1(%rsp)
        movapd %xmm2,nb333_dzH1(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb333_dxH2(%rsp)
        movapd %xmm4,nb333_dyH2(%rsp)
        movapd %xmm5,nb333_dzH2(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb333_dxM(%rsp)
        movapd %xmm7,nb333_dyM(%rsp)
        movapd %xmm8,nb333_dzM(%rsp)
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

        movapd  nb333_three(%rsp),%xmm9
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

        movapd  nb333_half(%rsp),%xmm15
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

        movapd  nb333_three(%rsp),%xmm1
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

        movapd  nb333_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvH1
        mulpd   %xmm15,%xmm10 ##   rinvH2
    mulpd   %xmm15,%xmm11 ##   rinvM

        movapd  %xmm9,nb333_rinvH1(%rsp)
        movapd  %xmm10,nb333_rinvH2(%rsp)
        movapd  %xmm11,nb333_rinvM(%rsp)

        ## interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd nb333_tsc(%rsp),%xmm1
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

    ## multiply by three (copy, mult. by two, add back)
    movapd  %xmm1,%xmm10
    movapd  %xmm4,%xmm11
    movapd  %xmm7,%xmm12
    pslld   $1,%xmm1
    pslld   $1,%xmm4
    pslld   $1,%xmm7
    paddd   %xmm10,%xmm1
    paddd   %xmm11,%xmm4
    paddd   %xmm12,%xmm7

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

    movq nb333_VFtab(%rbp),%rsi

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
    movapd nb333_qqH(%rsp),%xmm12
    movapd nb333_qqM(%rsp),%xmm13
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
    addpd  nb333_vctot(%rsp),%xmm1
    addpd  %xmm9,%xmm5
    addpd  %xmm5,%xmm1
    movapd %xmm1,nb333_vctot(%rsp)

    movapd nb333_tsc(%rsp),%xmm10
    mulpd  %xmm10,%xmm3 ## fscal
    mulpd  %xmm10,%xmm7
    mulpd  %xmm11,%xmm10

    xorpd  %xmm4,%xmm4
    xorpd  %xmm8,%xmm8
    xorpd  %xmm11,%xmm11

    subpd  %xmm3,%xmm4
    subpd  %xmm7,%xmm8
    subpd  %xmm10,%xmm11

    mulpd  nb333_rinvH1(%rsp),%xmm4
    mulpd  nb333_rinvH2(%rsp),%xmm8
    mulpd  nb333_rinvM(%rsp),%xmm11

    ## move j forces to xmm0-xmm2
    movq nb333_faction(%rbp),%rdi
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

    ## add forces from O interaction
    addpd nb333_fjx(%rsp),%xmm0
    addpd nb333_fjy(%rsp),%xmm1
    addpd nb333_fjz(%rsp),%xmm2

        mulpd nb333_dxH1(%rsp),%xmm3
        mulpd nb333_dyH1(%rsp),%xmm4
        mulpd nb333_dzH1(%rsp),%xmm5
        mulpd nb333_dxH2(%rsp),%xmm7
        mulpd nb333_dyH2(%rsp),%xmm8
        mulpd nb333_dzH2(%rsp),%xmm9
        mulpd nb333_dxM(%rsp),%xmm10
        mulpd nb333_dyM(%rsp),%xmm11
        mulpd nb333_dzM(%rsp),%xmm12

    addpd %xmm3,%xmm0
    addpd %xmm4,%xmm1
    addpd %xmm5,%xmm2
    addpd nb333_fixH1(%rsp),%xmm3
    addpd nb333_fiyH1(%rsp),%xmm4
    addpd nb333_fizH1(%rsp),%xmm5

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb333_fixH2(%rsp),%xmm7
    addpd nb333_fiyH2(%rsp),%xmm8
    addpd nb333_fizH2(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb333_fixM(%rsp),%xmm10
    addpd nb333_fiyM(%rsp),%xmm11
    addpd nb333_fizM(%rsp),%xmm12

    movapd %xmm3,nb333_fixH1(%rsp)
    movapd %xmm4,nb333_fiyH1(%rsp)
    movapd %xmm5,nb333_fizH1(%rsp)
    movapd %xmm7,nb333_fixH2(%rsp)
    movapd %xmm8,nb333_fiyH2(%rsp)
    movapd %xmm9,nb333_fizH2(%rsp)
    movapd %xmm10,nb333_fixM(%rsp)
    movapd %xmm11,nb333_fiyM(%rsp)
    movapd %xmm12,nb333_fizM(%rsp)

    ## store back j forces from xmm0-xmm2
        movlpd %xmm0,(%rdi,%rax,8)
        movlpd %xmm1,8(%rdi,%rax,8)
        movlpd %xmm2,16(%rdi,%rax,8)
        movhpd %xmm0,(%rdi,%rbx,8)
        movhpd %xmm1,8(%rdi,%rbx,8)
        movhpd %xmm2,16(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb333_innerk(%rsp)
        jl    _nb_kernel333_x86_64_sse2.nb333_checksingle
        jmp   _nb_kernel333_x86_64_sse2.nb333_unroll_loop
_nb_kernel333_x86_64_sse2.nb333_checksingle: 
        movl  nb333_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel333_x86_64_sse2.nb333_dosingle
        jmp   _nb_kernel333_x86_64_sse2.nb333_updateouterdata
_nb_kernel333_x86_64_sse2.nb333_dosingle: 
        movq  nb333_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb333_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb333_charge(%rbp),%rsi     ## base of charge[] 

        movsd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
        mulsd  nb333_iqM(%rsp),%xmm3
        mulsd  nb333_iqH(%rsp),%xmm4
        movapd  %xmm3,nb333_qqM(%rsp)
        movapd  %xmm4,nb333_qqH(%rsp)

        movq nb333_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movq nb333_vdwparam(%rbp),%rsi
        shll %r8d
        movl nb333_ntia(%rsp),%edi
        addl %edi,%r8d

        movlpd (%rsi,%r8,8),%xmm6       ## c6a
        movlpd 8(%rsi,%r8,8),%xmm7      ## c12a

        movapd %xmm6,nb333_c6(%rsp)
        movapd %xmm7,nb333_c12(%rsp)

        movq nb333_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move j coordinates to local temp variables 
    movsd (%rsi,%rax,8),%xmm0
    movsd 8(%rsi,%rax,8),%xmm1
    movsd 16(%rsi,%rax,8),%xmm2

    ## xmm0 = jx
    ## xmm1 = jy
    ## xmm2 = jz

    ## O interaction
    ## copy to xmm3-xmm5
    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5

    subsd nb333_ixO(%rsp),%xmm3
    subsd nb333_iyO(%rsp),%xmm4
    subsd nb333_izO(%rsp),%xmm5

    movapd %xmm3,nb333_dxO(%rsp)
    movapd %xmm4,nb333_dyO(%rsp)
    movapd %xmm5,nb333_dzO(%rsp)

        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5

        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
    ## xmm3=rsq

    cvtsd2ss %xmm3,%xmm5
    rsqrtss %xmm5,%xmm5
    cvtss2sd %xmm5,%xmm15    ## lu in low xmm2 

    ## lookup seed in xmm2 
    movapd %xmm15,%xmm5      ## copy of lu 
    mulsd %xmm15,%xmm15       ## lu*lu 
    movapd nb333_three(%rsp),%xmm7
    mulsd %xmm3,%xmm15       ## rsq*lu*lu                    
    movapd nb333_half(%rsp),%xmm6
    subsd %xmm15,%xmm7       ## 30-rsq*lu*lu 
    mulsd %xmm5,%xmm7
    mulsd %xmm6,%xmm7       ## xmm0=iter1 of rinv (new lu) 

    movapd %xmm7,%xmm5      ## copy of lu 
    mulsd %xmm7,%xmm7       ## lu*lu 
    movapd nb333_three(%rsp),%xmm15
    mulsd %xmm3,%xmm7       ## rsq*lu*lu                    
    movapd nb333_half(%rsp),%xmm6
    subsd %xmm7,%xmm15       ## 30-rsq*lu*lu 
    mulsd %xmm5,%xmm15
    mulsd %xmm6,%xmm15       ## xmm15=rinv

    mulsd %xmm15,%xmm3       ## xmm3=r 

    ## xmm15=rinv
    ## xmm3=r

    mulsd nb333_tsc(%rsp),%xmm3   ## rtab

    ## truncate and convert to integers
    cvttsd2si %xmm3,%r10d

    ## convert back to float
    cvtsi2sd  %r10d,%xmm4

    ## multiply by 4
    shll   $2,%r10d

    ## calculate eps
    subsd     %xmm4,%xmm3   ## xmm3=eps

    ## multiply by 3
    lea  (%r10,%r10,2),%r10

    ## xmm3=eps
    ## xmm15=rinv

        movq nb333_VFtab(%rbp),%rsi
    movsd  32(%rsi,%r10,8),%xmm4
    movsd  40(%rsi,%r10,8),%xmm5
    movsd  48(%rsi,%r10,8),%xmm6
    movsd  56(%rsi,%r10,8),%xmm7
    movsd  64(%rsi,%r10,8),%xmm8
    movsd  72(%rsi,%r10,8),%xmm9
    movsd  80(%rsi,%r10,8),%xmm10
    movsd  88(%rsi,%r10,8),%xmm11
    ## dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11

    mulsd  %xmm3,%xmm7   ## Heps
    mulsd  %xmm3,%xmm11
    mulsd  %xmm3,%xmm6  ## Geps
    mulsd  %xmm3,%xmm10
    mulsd  %xmm3,%xmm7  ## Heps2
    mulsd  %xmm3,%xmm11
    addsd  %xmm6,%xmm5 ## F+Geps
    addsd  %xmm10,%xmm9
    addsd  %xmm7,%xmm5  ## F+Geps+Heps2 = Fp
    addsd  %xmm11,%xmm9
    addsd  %xmm7,%xmm7   ## 2*Heps2
    addsd  %xmm11,%xmm11
    addsd  %xmm6,%xmm7  ## 2*Heps2+Geps
    addsd  %xmm10,%xmm11

    addsd  %xmm5,%xmm7 ## FF = Fp + 2*Heps2 + Geps
    addsd  %xmm9,%xmm11
    mulsd  %xmm3,%xmm5 ## eps*Fp
    mulsd  %xmm3,%xmm9
    movapd nb333_c6(%rsp),%xmm12
    movapd nb333_c12(%rsp),%xmm13
    addsd  %xmm4,%xmm5 ## VV
    addsd  %xmm8,%xmm9

    mulsd  %xmm12,%xmm5 ## VV*c6 = vnb6
    mulsd  %xmm13,%xmm9 ## VV*c12 = vnb12
    addsd  %xmm9,%xmm5
    addsd  nb333_Vvdwtot(%rsp),%xmm5
    movsd %xmm5,nb333_Vvdwtot(%rsp)

    mulsd  %xmm12,%xmm7  ## FF*c6 = fnb6
    mulsd  %xmm13,%xmm11  ## FF*c12  = fnb12
    addsd  %xmm11,%xmm7

    mulsd  nb333_tsc(%rsp),%xmm7
    mulsd  %xmm15,%xmm7  ## -fscal
    xorpd  %xmm9,%xmm9

    subsd  %xmm7,%xmm9    ## fscal
    movapd %xmm9,%xmm10
    movapd %xmm9,%xmm11

    mulsd  nb333_dxO(%rsp),%xmm9    ## fx/fy/fz
    mulsd  nb333_dyO(%rsp),%xmm10
    mulsd  nb333_dzO(%rsp),%xmm11

    ## save j force temporarily
    movapd %xmm9,nb333_fjx(%rsp)
    movapd %xmm10,nb333_fjy(%rsp)
    movapd %xmm11,nb333_fjz(%rsp)

    ## increment i O force
    addsd nb333_fixO(%rsp),%xmm9
    addsd nb333_fiyO(%rsp),%xmm10
    addsd nb333_fizO(%rsp),%xmm11
    movsd %xmm9,nb333_fixO(%rsp)
    movsd %xmm10,nb333_fiyO(%rsp)
    movsd %xmm11,nb333_fizO(%rsp)
    ## finished O LJ interaction.


    ## do H1, H2, and M interactions in parallel.
    ## xmm0-xmm2 still contain j coordinates.                
    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subsd nb333_ixH1(%rsp),%xmm0
    subsd nb333_iyH1(%rsp),%xmm1
    subsd nb333_izH1(%rsp),%xmm2
    subsd nb333_ixH2(%rsp),%xmm3
    subsd nb333_iyH2(%rsp),%xmm4
    subsd nb333_izH2(%rsp),%xmm5
    subsd nb333_ixM(%rsp),%xmm6
    subsd nb333_iyM(%rsp),%xmm7
    subsd nb333_izM(%rsp),%xmm8

        movapd %xmm0,nb333_dxH1(%rsp)
        movapd %xmm1,nb333_dyH1(%rsp)
        movapd %xmm2,nb333_dzH1(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb333_dxH2(%rsp)
        movapd %xmm4,nb333_dyH2(%rsp)
        movapd %xmm5,nb333_dzH2(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movapd %xmm6,nb333_dxM(%rsp)
        movapd %xmm7,nb333_dyM(%rsp)
        movapd %xmm8,nb333_dzM(%rsp)
        mulsd  %xmm6,%xmm6
        mulsd  %xmm7,%xmm7
        mulsd  %xmm8,%xmm8
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
    addsd  %xmm7,%xmm6
    addsd  %xmm8,%xmm6

        ## start doing invsqrt for j atoms
    cvtsd2ss %xmm0,%xmm1
    cvtsd2ss %xmm3,%xmm4
    cvtsd2ss %xmm6,%xmm7
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm4,%xmm4
    rsqrtss %xmm7,%xmm7
    cvtss2sd %xmm1,%xmm1
    cvtss2sd %xmm4,%xmm4
    cvtss2sd %xmm7,%xmm7

        movapd  %xmm1,%xmm2
        movapd  %xmm4,%xmm5
    movapd  %xmm7,%xmm8

        mulsd   %xmm1,%xmm1 ## lu*lu
        mulsd   %xmm4,%xmm4 ## lu*lu
    mulsd   %xmm7,%xmm7 ## lu*lu

        movapd  nb333_three(%rsp),%xmm9
        movapd  %xmm9,%xmm10
    movapd  %xmm9,%xmm11

        mulsd   %xmm0,%xmm1 ## rsq*lu*lu
        mulsd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm7 ## rsq*lu*lu

        subsd   %xmm1,%xmm9
        subsd   %xmm4,%xmm10
    subsd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulsd   %xmm2,%xmm9
        mulsd   %xmm5,%xmm10
    mulsd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb333_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ## first iteration for rinvH1
        mulsd   %xmm15,%xmm10 ## first iteration for rinvH2
    mulsd   %xmm15,%xmm11 ## first iteration for rinvM

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulsd   %xmm2,%xmm2 ## lu*lu
        mulsd   %xmm5,%xmm5 ## lu*lu
    mulsd   %xmm8,%xmm8 ## lu*lu

        movapd  nb333_three(%rsp),%xmm1
        movapd  %xmm1,%xmm4
    movapd  %xmm1,%xmm7

        mulsd   %xmm0,%xmm2 ## rsq*lu*lu
        mulsd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm8 ## rsq*lu*lu

        subsd   %xmm2,%xmm1
        subsd   %xmm5,%xmm4
    subsd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulsd   %xmm1,%xmm9
        mulsd   %xmm4,%xmm10
    mulsd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb333_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvH1
        mulsd   %xmm15,%xmm10 ##   rinvH2
    mulsd   %xmm15,%xmm11 ##   rinvM

        movapd  %xmm9,nb333_rinvH1(%rsp)
        movapd  %xmm10,nb333_rinvH2(%rsp)
        movapd  %xmm11,nb333_rinvM(%rsp)

        ## interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd nb333_tsc(%rsp),%xmm1
    mulsd  %xmm9,%xmm0 ## r
    mulsd  %xmm10,%xmm3
    mulsd  %xmm11,%xmm6
    mulsd  %xmm1,%xmm0 ## rtab
    mulsd  %xmm1,%xmm3
    mulsd  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttsd2si %xmm0,%r8d
    cvttsd2si %xmm3,%r10d
    cvttsd2si %xmm6,%r12d

    ## convert back to float
    cvtsi2sd  %r8d,%xmm2
    cvtsi2sd  %r10d,%xmm5
    cvtsi2sd  %r12d,%xmm8

    ## multiply by 4
    shll  $2,%r8d
    shll  $2,%r10d
    shll  $2,%r12d

    movq nb333_VFtab(%rbp),%rsi

        lea  (%r8,%r8,2),%r8
    lea  (%r10,%r10,2),%r10
    lea  (%r12,%r12,2),%r12

    ## calculate eps
    subsd     %xmm2,%xmm0
    subsd     %xmm5,%xmm3
    subsd     %xmm8,%xmm6

    movapd    %xmm0,%xmm12 ## epsH1
    movapd    %xmm3,%xmm13 ## epsH2
    movapd    %xmm6,%xmm14 ## epsM

    ## Load LOTS of table data
    movsd (%rsi,%r8,8),%xmm0
    movsd 8(%rsi,%r8,8),%xmm1
    movsd 16(%rsi,%r8,8),%xmm2
    movsd 24(%rsi,%r8,8),%xmm3
    movsd (%rsi,%r10,8),%xmm4
    movsd 8(%rsi,%r10,8),%xmm5
    movsd 16(%rsi,%r10,8),%xmm6
    movsd 24(%rsi,%r10,8),%xmm7
    movsd (%rsi,%r12,8),%xmm8
    movsd 8(%rsi,%r12,8),%xmm9
    movsd 16(%rsi,%r12,8),%xmm10
    movsd 24(%rsi,%r12,8),%xmm11
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    mulsd  %xmm12,%xmm3  ## Heps
    mulsd  %xmm13,%xmm7
    mulsd  %xmm14,%xmm11
    mulsd  %xmm12,%xmm2  ## Geps
    mulsd  %xmm13,%xmm6
    mulsd  %xmm14,%xmm10
    mulsd  %xmm12,%xmm3  ## Heps2
    mulsd  %xmm13,%xmm7
    mulsd  %xmm14,%xmm11

    addsd  %xmm2,%xmm1  ## F+Geps
    addsd  %xmm6,%xmm5
    addsd  %xmm10,%xmm9
    addsd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addsd  %xmm7,%xmm5
    addsd  %xmm11,%xmm9
    addsd  %xmm3,%xmm3   ## 2*Heps2
    addsd  %xmm7,%xmm7
    addsd  %xmm11,%xmm11
    addsd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addsd  %xmm6,%xmm7
    addsd  %xmm10,%xmm11
    addsd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addsd  %xmm5,%xmm7
    addsd  %xmm9,%xmm11
    mulsd  %xmm12,%xmm1  ## eps*Fp
    mulsd  %xmm13,%xmm5
    mulsd  %xmm14,%xmm9
    movapd nb333_qqH(%rsp),%xmm12
    movapd nb333_qqM(%rsp),%xmm13
    addsd  %xmm0,%xmm1    ## VV
    addsd  %xmm4,%xmm5
    addsd  %xmm8,%xmm9
    mulsd  %xmm12,%xmm1  ## VV*qq = vcoul
    mulsd  %xmm12,%xmm5
    mulsd  %xmm13,%xmm9
    mulsd  %xmm12,%xmm3   ## FF*qq = fij
    mulsd  %xmm12,%xmm7
    mulsd  %xmm13,%xmm11

    ## accumulate vctot
    addsd  nb333_vctot(%rsp),%xmm1
    addsd  %xmm9,%xmm5
    addsd  %xmm5,%xmm1
    movsd %xmm1,nb333_vctot(%rsp)

    movapd nb333_tsc(%rsp),%xmm10
    mulsd  %xmm10,%xmm3 ## fscal
    mulsd  %xmm10,%xmm7
    mulsd  %xmm11,%xmm10

    xorpd  %xmm4,%xmm4
    xorpd  %xmm8,%xmm8
    xorpd  %xmm11,%xmm11

    subsd  %xmm3,%xmm4
    subsd  %xmm7,%xmm8
    subsd  %xmm10,%xmm11

    mulsd  nb333_rinvH1(%rsp),%xmm4
    mulsd  nb333_rinvH2(%rsp),%xmm8
    mulsd  nb333_rinvM(%rsp),%xmm11

    ## move j forces to xmm0-xmm2
    movq nb333_faction(%rbp),%rdi
        movsd (%rdi,%rax,8),%xmm0
        movsd 8(%rdi,%rax,8),%xmm1
        movsd 16(%rdi,%rax,8),%xmm2

    movapd %xmm4,%xmm3
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm7
    movapd %xmm8,%xmm9
    movapd %xmm11,%xmm10
    movapd %xmm11,%xmm12

    ## add forces from O interaction
    addsd nb333_fjx(%rsp),%xmm0
    addsd nb333_fjy(%rsp),%xmm1
    addsd nb333_fjz(%rsp),%xmm2

        mulsd nb333_dxH1(%rsp),%xmm3
        mulsd nb333_dyH1(%rsp),%xmm4
        mulsd nb333_dzH1(%rsp),%xmm5
        mulsd nb333_dxH2(%rsp),%xmm7
        mulsd nb333_dyH2(%rsp),%xmm8
        mulsd nb333_dzH2(%rsp),%xmm9
        mulsd nb333_dxM(%rsp),%xmm10
        mulsd nb333_dyM(%rsp),%xmm11
        mulsd nb333_dzM(%rsp),%xmm12

    addsd %xmm3,%xmm0
    addsd %xmm4,%xmm1
    addsd %xmm5,%xmm2
    addsd nb333_fixH1(%rsp),%xmm3
    addsd nb333_fiyH1(%rsp),%xmm4
    addsd nb333_fizH1(%rsp),%xmm5

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb333_fixH2(%rsp),%xmm7
    addsd nb333_fiyH2(%rsp),%xmm8
    addsd nb333_fizH2(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb333_fixM(%rsp),%xmm10
    addsd nb333_fiyM(%rsp),%xmm11
    addsd nb333_fizM(%rsp),%xmm12

    movsd %xmm3,nb333_fixH1(%rsp)
    movsd %xmm4,nb333_fiyH1(%rsp)
    movsd %xmm5,nb333_fizH1(%rsp)
    movsd %xmm7,nb333_fixH2(%rsp)
    movsd %xmm8,nb333_fiyH2(%rsp)
    movsd %xmm9,nb333_fizH2(%rsp)
    movsd %xmm10,nb333_fixM(%rsp)
    movsd %xmm11,nb333_fiyM(%rsp)
    movsd %xmm12,nb333_fizM(%rsp)

    ## store back j forces from xmm0-xmm2
        movsd %xmm0,(%rdi,%rax,8)
        movsd %xmm1,8(%rdi,%rax,8)
        movsd %xmm2,16(%rdi,%rax,8)

_nb_kernel333_x86_64_sse2.nb333_updateouterdata: 
        movl  nb333_ii3(%rsp),%ecx
        movq  nb333_faction(%rbp),%rdi
        movq  nb333_fshift(%rbp),%rsi
        movl  nb333_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb333_fixO(%rsp),%xmm0
        movapd nb333_fiyO(%rsp),%xmm1
        movapd nb333_fizO(%rsp),%xmm2

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
        movapd nb333_fixH1(%rsp),%xmm0
        movapd nb333_fiyH1(%rsp),%xmm1
        movapd nb333_fizH1(%rsp),%xmm2

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
        movapd nb333_fixH2(%rsp),%xmm0
        movapd nb333_fiyH2(%rsp),%xmm1
        movapd nb333_fizH2(%rsp),%xmm2

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

        ## accumulate Mi forces in xmm0, xmm1, xmm2 
        movapd nb333_fixM(%rsp),%xmm0
        movapd nb333_fiyM(%rsp),%xmm1
        movapd nb333_fizM(%rsp),%xmm2

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
        movl nb333_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb333_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb333_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb333_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb333_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb333_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb333_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel333_x86_64_sse2.nb333_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb333_n(%rsp)
        jmp _nb_kernel333_x86_64_sse2.nb333_outer
_nb_kernel333_x86_64_sse2.nb333_outerend: 
        ## check if more outer neighborlists remain
        movl  nb333_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel333_x86_64_sse2.nb333_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel333_x86_64_sse2.nb333_threadloop
_nb_kernel333_x86_64_sse2.nb333_end: 
        movl nb333_nouter(%rsp),%eax
        movl nb333_ninner(%rsp),%ebx
        movq nb333_outeriter(%rbp),%rcx
        movq nb333_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1048,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret







.globl nb_kernel333nf_x86_64_sse2
.globl _nb_kernel333nf_x86_64_sse2
nb_kernel333nf_x86_64_sse2:     
_nb_kernel333nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb333nf_fshift, 16
.set nb333nf_gid, 24
.set nb333nf_pos, 32
.set nb333nf_faction, 40
.set nb333nf_charge, 48
.set nb333nf_p_facel, 56
.set nb333nf_argkrf, 64
.set nb333nf_argcrf, 72
.set nb333nf_Vc, 80
.set nb333nf_type, 88
.set nb333nf_p_ntype, 96
.set nb333nf_vdwparam, 104
.set nb333nf_Vvdw, 112
.set nb333nf_p_tabscale, 120
.set nb333nf_VFtab, 128
.set nb333nf_invsqrta, 136
.set nb333nf_dvda, 144
.set nb333nf_p_gbtabscale, 152
.set nb333nf_GBtab, 160
.set nb333nf_p_nthreads, 168
.set nb333nf_count, 176
.set nb333nf_mtx, 184
.set nb333nf_outeriter, 192
.set nb333nf_inneriter, 200
.set nb333nf_work, 208
        ## stack offsets for local variables 
        ## bottom of stack is cache-aligned for sse2 use 
.set nb333nf_ixO, 0
.set nb333nf_iyO, 16
.set nb333nf_izO, 32
.set nb333nf_ixH1, 48
.set nb333nf_iyH1, 64
.set nb333nf_izH1, 80
.set nb333nf_ixH2, 96
.set nb333nf_iyH2, 112
.set nb333nf_izH2, 128
.set nb333nf_ixM, 144
.set nb333nf_iyM, 160
.set nb333nf_izM, 176
.set nb333nf_iqM, 192
.set nb333nf_iqH, 208
.set nb333nf_qqM, 224
.set nb333nf_qqH, 240
.set nb333nf_rO, 256
.set nb333nf_rH1, 272
.set nb333nf_rH2, 288
.set nb333nf_rM, 304
.set nb333nf_tsc, 320
.set nb333nf_c6, 336
.set nb333nf_c12, 352
.set nb333nf_vctot, 368
.set nb333nf_Vvdwtot, 384
.set nb333nf_half, 400
.set nb333nf_three, 416
.set nb333nf_is3, 432
.set nb333nf_ii3, 436
.set nb333nf_nri, 440
.set nb333nf_iinr, 448
.set nb333nf_jindex, 456
.set nb333nf_jjnr, 464
.set nb333nf_shift, 472
.set nb333nf_shiftvec, 480
.set nb333nf_facel, 488
.set nb333nf_innerjjnr, 496
.set nb333nf_ntia, 504
.set nb333nf_innerk, 508
.set nb333nf_n, 512
.set nb333nf_nn1, 516
.set nb333nf_nouter, 520
.set nb333nf_ninner, 524
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $536,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb333nf_nouter(%rsp)
        movl %eax,nb333nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb333nf_nri(%rsp)
        movq %rsi,nb333nf_iinr(%rsp)
        movq %rdx,nb333nf_jindex(%rsp)
        movq %rcx,nb333nf_jjnr(%rsp)
        movq %r8,nb333nf_shift(%rsp)
        movq %r9,nb333nf_shiftvec(%rsp)
        movq nb333nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb333nf_facel(%rsp)

        movq nb333nf_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb333nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb333nf_half(%rsp)
        movl %ebx,nb333nf_half+4(%rsp)
        movsd nb333nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb333nf_half(%rsp)
        movapd %xmm3,nb333nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb333nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb333nf_charge(%rbp),%rdx
        movsd 8(%rdx,%rbx,8),%xmm3
        movsd 24(%rdx,%rbx,8),%xmm4
        movq nb333nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb333nf_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb333nf_iqH(%rsp)
        movapd %xmm4,nb333nf_iqM(%rsp)

        movq  nb333nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb333nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb333nf_ntia(%rsp)
_nb_kernel333nf_x86_64_sse2.nb333nf_threadloop: 
        movq  nb333nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel333nf_x86_64_sse2.nb333nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel333nf_x86_64_sse2.nb333nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb333nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb333nf_n(%rsp)
        movl %ebx,nb333nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel333nf_x86_64_sse2.nb333nf_outerstart
        jmp _nb_kernel333nf_x86_64_sse2.nb333nf_end

_nb_kernel333nf_x86_64_sse2.nb333nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb333nf_nouter(%rsp),%ebx
        movl %ebx,nb333nf_nouter(%rsp)

_nb_kernel333nf_x86_64_sse2.nb333nf_outer: 
        movq  nb333nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb333nf_is3(%rsp)            ## store is3 

        movq  nb333nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb333nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb333nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb333nf_ii3(%rsp)

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
        movapd %xmm3,nb333nf_ixO(%rsp)
        movapd %xmm4,nb333nf_iyO(%rsp)
        movapd %xmm5,nb333nf_izO(%rsp)
        movapd %xmm6,nb333nf_ixH1(%rsp)
        movapd %xmm7,nb333nf_iyH1(%rsp)

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
        movapd %xmm6,nb333nf_izH1(%rsp)
        movapd %xmm0,nb333nf_ixH2(%rsp)
        movapd %xmm1,nb333nf_iyH2(%rsp)
        movapd %xmm2,nb333nf_izH2(%rsp)
        movapd %xmm3,nb333nf_ixM(%rsp)
        movapd %xmm4,nb333nf_iyM(%rsp)
        movapd %xmm5,nb333nf_izM(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb333nf_vctot(%rsp)
        movapd %xmm4,nb333nf_Vvdwtot(%rsp)

        movq  nb333nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx     ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb333nf_pos(%rbp),%rsi
        movq  nb333nf_faction(%rbp),%rdi
        movq  nb333nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb333nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb333nf_ninner(%rsp),%ecx
        movl  %ecx,nb333nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb333nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel333nf_x86_64_sse2.nb333nf_unroll_loop
        jmp   _nb_kernel333nf_x86_64_sse2.nb333nf_checksingle
_nb_kernel333nf_x86_64_sse2.nb333nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb333nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb333nf_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb333nf_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb333nf_iqM(%rsp),%xmm3
        mulpd  nb333nf_iqH(%rsp),%xmm4
        movapd  %xmm3,nb333nf_qqM(%rsp)
        movapd  %xmm4,nb333nf_qqH(%rsp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movq nb333nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movq nb333nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        movl nb333nf_ntia(%rsp),%edi
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
        movapd %xmm4,nb333nf_c6(%rsp)
        movapd %xmm6,nb333nf_c12(%rsp)

        movq nb333nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movapd nb333nf_ixO(%rsp),%xmm4
        movapd nb333nf_iyO(%rsp),%xmm5
        movapd nb333nf_izO(%rsp),%xmm6

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
        movapd nb333nf_ixH1(%rsp),%xmm4
        movapd nb333nf_iyH1(%rsp),%xmm5
        movapd nb333nf_izH1(%rsp),%xmm6

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
        movapd nb333nf_ixH2(%rsp),%xmm3
        movapd nb333nf_iyH2(%rsp),%xmm4
        movapd nb333nf_izH2(%rsp),%xmm5

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
        movapd nb333nf_iyM(%rsp),%xmm3
        movapd nb333nf_izM(%rsp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb333nf_ixM(%rsp),%xmm2
        subpd  %xmm0,%xmm2

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## rsqO - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb333nf_three(%rsp),%xmm0
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulpd   nb333nf_half(%rsp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm7,%xmm2
        movapd %xmm0,%xmm3
        mulpd %xmm0,%xmm0       ## lu*lu 
        mulpd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333nf_three(%rsp),%xmm0
        subpd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulpd nb333nf_half(%rsp),%xmm0   ## rinv 
        mulpd   %xmm0,%xmm7
        movapd  %xmm7,nb333nf_rO(%rsp)          ## r in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb333nf_three(%rsp),%xmm0
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulpd   nb333nf_half(%rsp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm0,%xmm3
        mulpd %xmm0,%xmm0       ## lu*lu 
        mulpd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333nf_three(%rsp),%xmm0
        subpd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulpd nb333nf_half(%rsp),%xmm0   ## rinv 
        mulpd  %xmm0,%xmm6
        movapd %xmm6,nb333nf_rH1(%rsp)          ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb333nf_three(%rsp),%xmm0
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulpd   nb333nf_half(%rsp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm0,%xmm3
        mulpd %xmm0,%xmm0       ## lu*lu 
        mulpd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333nf_three(%rsp),%xmm0
        subpd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulpd nb333nf_half(%rsp),%xmm0   ## rinv 
        mulpd %xmm0,%xmm5
        movapd %xmm5,nb333nf_rH2(%rsp)   ## r 

        ## rsqM - seed in xmm2 
        cvtpd2ps %xmm4,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb333nf_three(%rsp),%xmm0
        mulpd   %xmm4,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulpd   nb333nf_half(%rsp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm4,%xmm2
        movapd %xmm0,%xmm3
        mulpd %xmm0,%xmm0       ## lu*lu 
        mulpd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333nf_three(%rsp),%xmm0
        subpd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulpd nb333nf_half(%rsp),%xmm0   ## rinv 
        mulpd %xmm0,%xmm4
        movapd %xmm4,nb333nf_rM(%rsp)   ## r 

        ## do O interactions 
        ## rO is still in xmm7 
        mulpd   nb333nf_tsc(%rsp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movd %eax,%mm0
        movd %ebx,%mm1
        movq nb333nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        lea  (%rax,%rax,2),%rax ## idx *= 3 (total *=12 now) 
        lea  (%rbx,%rbx,2),%rbx

        ## Dispersion 
        movapd 32(%rsi,%rax,8),%xmm4    ## Y1 F1        
        movapd 32(%rsi,%rbx,8),%xmm3    ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 48(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 48(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb333nf_c6(%rsp),%xmm4
        mulpd  %xmm4,%xmm5       ## Vvdw6               
        addpd  nb333nf_Vvdwtot(%rsp),%xmm5
        movapd %xmm5,nb333nf_Vvdwtot(%rsp)

        ## Repulsion 
        movapd 64(%rsi,%rax,8),%xmm4    ## Y1 F1        
        movapd 64(%rsi,%rbx,8),%xmm3    ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 80(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 80(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb333nf_c12(%rsp),%xmm4
        mulpd  %xmm4,%xmm5 ## Vvdw12 
        addpd  nb333nf_Vvdwtot(%rsp),%xmm5
        movapd %xmm5,nb333nf_Vvdwtot(%rsp)

        ## Done with O interactions - now H1! 
        movapd nb333nf_rH1(%rsp),%xmm7
        mulpd nb333nf_tsc(%rsp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb333nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        lea  (%rax,%rax,2),%rax ## idx *= 3 (total *=12 now)   
        lea  (%rbx,%rbx,2),%rbx

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
        movapd nb333nf_qqH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addpd  nb333nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb333nf_vctot(%rsp)

        ## H2 interactions 
        movapd nb333nf_rH2(%rsp),%xmm7
        mulpd   nb333nf_tsc(%rsp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb333nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        lea  (%rax,%rax,2),%rax ## idx *= 3 (total *=12 now)
        lea  (%rbx,%rbx,2),%rbx

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
        movapd nb333nf_qqH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addpd  nb333nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb333nf_vctot(%rsp)

        ## M interactions 
        movapd nb333nf_rM(%rsp),%xmm7
        mulpd   nb333nf_tsc(%rsp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb333nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        lea  (%rax,%rax,2),%rax ## idx *= 3 (total *=12 now)
        lea  (%rbx,%rbx,2),%rbx

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
        movapd nb333nf_qqM(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addpd  nb333nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb333nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb333nf_innerk(%rsp)
        jl    _nb_kernel333nf_x86_64_sse2.nb333nf_checksingle
        jmp   _nb_kernel333nf_x86_64_sse2.nb333nf_unroll_loop
_nb_kernel333nf_x86_64_sse2.nb333nf_checksingle: 
        movl  nb333nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel333nf_x86_64_sse2.nb333nf_dosingle
        jmp   _nb_kernel333nf_x86_64_sse2.nb333nf_updateouterdata
_nb_kernel333nf_x86_64_sse2.nb333nf_dosingle: 
        movq  nb333nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb333nf_charge(%rbp),%rsi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb333nf_iqM(%rsp),%xmm3
        mulpd  nb333nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movapd  %xmm3,nb333nf_qqM(%rsp)
        movapd  %xmm4,nb333nf_qqH(%rsp)

        movq nb333nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movq nb333nf_vdwparam(%rbp),%rsi
        shll %eax
        movl nb333nf_ntia(%rsp),%edi
        addl %edi,%eax

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb333nf_c6(%rsp)
        movapd %xmm6,nb333nf_c12(%rsp)

        movq nb333nf_pos(%rbp),%rsi        ## base of pos[] 
        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coords to xmm0-xmm2 
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb333nf_ixO(%rsp),%xmm4
        movapd nb333nf_iyO(%rsp),%xmm5
        movapd nb333nf_izO(%rsp),%xmm6

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
        movapd nb333nf_ixH1(%rsp),%xmm4
        movapd nb333nf_iyH1(%rsp),%xmm5
        movapd nb333nf_izH1(%rsp),%xmm6

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
        movapd nb333nf_ixH2(%rsp),%xmm3
        movapd nb333nf_iyH2(%rsp),%xmm4
        movapd nb333nf_izH2(%rsp),%xmm5

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
        movapd nb333nf_iyM(%rsp),%xmm3
        movapd nb333nf_izM(%rsp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb333nf_ixM(%rsp),%xmm2
        subpd  %xmm0,%xmm2

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb333nf_three(%rsp),%xmm0
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb333nf_half(%rsp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm7,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333nf_three(%rsp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb333nf_half(%rsp),%xmm0   ## rinv 
        mulsd   %xmm0,%xmm7
        movapd  %xmm7,nb333nf_rO(%rsp)          ## r in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb333nf_three(%rsp),%xmm0
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb333nf_half(%rsp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333nf_three(%rsp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb333nf_half(%rsp),%xmm0   ## rinv 
        mulsd  %xmm0,%xmm6
        movapd %xmm6,nb333nf_rH1(%rsp)          ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb333nf_three(%rsp),%xmm0
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb333nf_half(%rsp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333nf_three(%rsp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb333nf_half(%rsp),%xmm0   ## rinv 
        mulsd %xmm0,%xmm5
        movapd %xmm5,nb333nf_rH2(%rsp)   ## r 

        ## rsqM - seed in xmm2 
        cvtsd2ss %xmm4,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb333nf_three(%rsp),%xmm0
        mulsd   %xmm4,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb333nf_half(%rsp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm4,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb333nf_three(%rsp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb333nf_half(%rsp),%xmm0   ## rinv 
        mulsd %xmm0,%xmm4
        movapd %xmm4,nb333nf_rM(%rsp)   ## r 

        ## do O interactions 
        movd %eax,%mm0
        ## rO is still in xmm7 
        mulsd   nb333nf_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%eax    ## lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm7,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax
        movq nb333nf_VFtab(%rbp),%rsi
        lea  (%rax,%rax,2),%rax ## idx *= 3 (total *=12 now)   

        ## Dispersion 
        movapd 32(%rsi,%rax,8),%xmm4    ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 48(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb333nf_c6(%rsp),%xmm4
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        addsd  nb333nf_Vvdwtot(%rsp),%xmm5
        movsd %xmm5,nb333nf_Vvdwtot(%rsp)

        ## Repulsion 
        movapd 64(%rsi,%rax,8),%xmm4    ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 80(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb333nf_c12(%rsp),%xmm4
        mulsd  %xmm4,%xmm5 ## Vvdw12 
        addsd  nb333nf_Vvdwtot(%rsp),%xmm5
        movsd %xmm5,nb333nf_Vvdwtot(%rsp)

        ## Done with O interactions - now H1! 
        movapd nb333nf_rH1(%rsp),%xmm7
        mulpd nb333nf_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb333nf_VFtab(%rbp),%rsi
        lea  (%rax,%rax,2),%rax ## idx *= 3 (total *=12 now)   

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
        movapd nb333nf_qqH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addsd  nb333nf_vctot(%rsp),%xmm5
        movlpd %xmm5,nb333nf_vctot(%rsp)

        ##  H2 interactions 
        movapd nb333nf_rH2(%rsp),%xmm7
        mulsd   nb333nf_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb333nf_VFtab(%rbp),%rsi
        lea  (%rax,%rax,2),%rax ## idx *= 3 (total *=12 now)   

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
        movapd nb333nf_qqH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addsd  nb333nf_vctot(%rsp),%xmm5
        movlpd %xmm5,nb333nf_vctot(%rsp)

        ## M interactions 
        movapd nb333nf_rM(%rsp),%xmm7
        mulsd   nb333nf_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb333nf_VFtab(%rbp),%rsi
        lea  (%rax,%rax,2),%rax ## idx *= 3 (total *=12 now)   

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
        movapd nb333nf_qqM(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addsd  nb333nf_vctot(%rsp),%xmm5
        movlpd %xmm5,nb333nf_vctot(%rsp)

_nb_kernel333nf_x86_64_sse2.nb333nf_updateouterdata: 
        ## get n from stack
        movl nb333nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb333nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb333nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb333nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb333nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb333nf_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb333nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel333nf_x86_64_sse2.nb333nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb333nf_n(%rsp)
        jmp _nb_kernel333nf_x86_64_sse2.nb333nf_outer
_nb_kernel333nf_x86_64_sse2.nb333nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb333nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel333nf_x86_64_sse2.nb333nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel333nf_x86_64_sse2.nb333nf_threadloop
_nb_kernel333nf_x86_64_sse2.nb333nf_end: 
        movl nb333nf_nouter(%rsp),%eax
        movl nb333nf_ninner(%rsp),%ebx
        movq nb333nf_outeriter(%rbp),%rcx
        movq nb333nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $536,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret

