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





.globl nb_kernel234_x86_64_sse2
.globl _nb_kernel234_x86_64_sse2
nb_kernel234_x86_64_sse2:       
_nb_kernel234_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb234_fshift, 16
.set nb234_gid, 24
.set nb234_pos, 32
.set nb234_faction, 40
.set nb234_charge, 48
.set nb234_p_facel, 56
.set nb234_argkrf, 64
.set nb234_argcrf, 72
.set nb234_Vc, 80
.set nb234_type, 88
.set nb234_p_ntype, 96
.set nb234_vdwparam, 104
.set nb234_Vvdw, 112
.set nb234_p_tabscale, 120
.set nb234_VFtab, 128
.set nb234_invsqrta, 136
.set nb234_dvda, 144
.set nb234_p_gbtabscale, 152
.set nb234_GBtab, 160
.set nb234_p_nthreads, 168
.set nb234_count, 176
.set nb234_mtx, 184
.set nb234_outeriter, 192
.set nb234_inneriter, 200
.set nb234_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb234_ixO, 0
.set nb234_iyO, 16
.set nb234_izO, 32
.set nb234_ixH1, 48
.set nb234_iyH1, 64
.set nb234_izH1, 80
.set nb234_ixH2, 96
.set nb234_iyH2, 112
.set nb234_izH2, 128
.set nb234_ixM, 144
.set nb234_iyM, 160
.set nb234_izM, 176
.set nb234_jxO, 192
.set nb234_jyO, 208
.set nb234_jzO, 224
.set nb234_jxH1, 240
.set nb234_jyH1, 256
.set nb234_jzH1, 272
.set nb234_jxH2, 288
.set nb234_jyH2, 304
.set nb234_jzH2, 320
.set nb234_jxM, 336
.set nb234_jyM, 352
.set nb234_jzM, 368
.set nb234_dxOO, 384
.set nb234_dyOO, 400
.set nb234_dzOO, 416
.set nb234_dxH1H1, 432
.set nb234_dyH1H1, 448
.set nb234_dzH1H1, 464
.set nb234_dxH1H2, 480
.set nb234_dyH1H2, 496
.set nb234_dzH1H2, 512
.set nb234_dxH1M, 528
.set nb234_dyH1M, 544
.set nb234_dzH1M, 560
.set nb234_dxH2H1, 576
.set nb234_dyH2H1, 592
.set nb234_dzH2H1, 608
.set nb234_dxH2H2, 624
.set nb234_dyH2H2, 640
.set nb234_dzH2H2, 656
.set nb234_dxH2M, 672
.set nb234_dyH2M, 688
.set nb234_dzH2M, 704
.set nb234_dxMH1, 720
.set nb234_dyMH1, 736
.set nb234_dzMH1, 752
.set nb234_dxMH2, 768
.set nb234_dyMH2, 784
.set nb234_dzMH2, 800
.set nb234_dxMM, 816
.set nb234_dyMM, 832
.set nb234_dzMM, 848
.set nb234_qqMM, 864
.set nb234_qqMH, 880
.set nb234_qqHH, 896
.set nb234_two, 912
.set nb234_c6, 944
.set nb234_c12, 960
.set nb234_vctot, 976
.set nb234_Vvdwtot, 992
.set nb234_fixO, 1008
.set nb234_fiyO, 1024
.set nb234_fizO, 1040
.set nb234_fixH1, 1056
.set nb234_fiyH1, 1072
.set nb234_fizH1, 1088
.set nb234_fixH2, 1104
.set nb234_fiyH2, 1120
.set nb234_fizH2, 1136
.set nb234_fixM, 1152
.set nb234_fiyM, 1168
.set nb234_fizM, 1184
.set nb234_fjxO, 1200
.set nb234_fjyO, 1216
.set nb234_fjzO, 1232
.set nb234_fjxH1, 1248
.set nb234_fjyH1, 1264
.set nb234_fjzH1, 1280
.set nb234_fjxH2, 1296
.set nb234_fjyH2, 1312
.set nb234_fjzH2, 1328
.set nb234_fjxM, 1344
.set nb234_fjyM, 1360
.set nb234_fjzM, 1376
.set nb234_half, 1392
.set nb234_three, 1408
.set nb234_tsc, 1424
.set nb234_fstmp, 1440
.set nb234_rsqOO, 1456
.set nb234_rsqH1H1, 1472
.set nb234_rsqH1H2, 1488
.set nb234_rsqH1M, 1504
.set nb234_rsqH2H1, 1520
.set nb234_rsqH2H2, 1536
.set nb234_rsqH2M, 1552
.set nb234_rsqMH1, 1568
.set nb234_rsqMH2, 1584
.set nb234_rsqMM, 1600
.set nb234_rinvOO, 1616
.set nb234_rinvH1H1, 1632
.set nb234_rinvH1H2, 1648
.set nb234_rinvH1M, 1664
.set nb234_rinvH2H1, 1680
.set nb234_rinvH2H2, 1696
.set nb234_rinvH2M, 1712
.set nb234_rinvMH1, 1728
.set nb234_rinvMH2, 1744
.set nb234_rinvMM, 1760
.set nb234_krf, 1776
.set nb234_crf, 1792
.set nb234_is3, 1808
.set nb234_ii3, 1812
.set nb234_nri, 1816
.set nb234_iinr, 1824
.set nb234_jindex, 1832
.set nb234_jjnr, 1840
.set nb234_shift, 1848
.set nb234_shiftvec, 1856
.set nb234_facel, 1864
.set nb234_innerjjnr, 1872
.set nb234_innerk, 1880
.set nb234_n, 1884
.set nb234_nn1, 1888
.set nb234_nouter, 1892
.set nb234_ninner, 1896
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1912,%rsp         ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb234_nouter(%rsp)
        movl %eax,nb234_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb234_nri(%rsp)
        movq %rsi,nb234_iinr(%rsp)
        movq %rdx,nb234_jindex(%rsp)
        movq %rcx,nb234_jjnr(%rsp)
        movq %r8,nb234_shift(%rsp)
        movq %r9,nb234_shiftvec(%rsp)
        movq nb234_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb234_facel(%rsp)

        movq nb234_argkrf(%rbp),%rsi
        movq nb234_argcrf(%rbp),%rdi
        movsd (%rsi),%xmm1
        movsd (%rdi),%xmm2
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        movapd %xmm1,nb234_krf(%rsp)
        movapd %xmm2,nb234_crf(%rsp)

        movq nb234_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb234_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb234_half(%rsp)
        movl %ebx,nb234_half+4(%rsp)
        movsd nb234_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb234_half(%rsp)
        movapd %xmm2,nb234_two(%rsp)
        movapd %xmm3,nb234_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb234_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb234_charge(%rbp),%rdx
        movsd 24(%rdx,%rbx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%rdx,%rbx,8),%xmm5

        movsd nb234_facel(%rsp),%xmm6
        mulsd  %xmm3,%xmm3
        mulsd  %xmm5,%xmm4
        mulsd  %xmm5,%xmm5
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb234_qqMM(%rsp)
        movapd %xmm4,nb234_qqMH(%rsp)
        movapd %xmm5,nb234_qqHH(%rsp)

        xorpd %xmm0,%xmm0
        movq  nb234_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq nb234_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movq  nb234_vdwparam(%rbp),%rax
        movlpd (%rax,%rdx,8),%xmm0
        movlpd 8(%rax,%rdx,8),%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb234_c6(%rsp)
        movapd %xmm1,nb234_c12(%rsp)

_nb_kernel234_x86_64_sse2.nb234_threadloop: 
        movq  nb234_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel234_x86_64_sse2.nb234_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel234_x86_64_sse2.nb234_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb234_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb234_n(%rsp)
        movl %ebx,nb234_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel234_x86_64_sse2.nb234_outerstart
        jmp _nb_kernel234_x86_64_sse2.nb234_end

_nb_kernel234_x86_64_sse2.nb234_outerstart: 
        ## ebx contains number of outer iterations
        addl nb234_nouter(%rsp),%ebx
        movl %ebx,nb234_nouter(%rsp)

_nb_kernel234_x86_64_sse2.nb234_outer: 
        movq  nb234_shift(%rsp),%rax        ## eax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb234_is3(%rsp)      ## store is3 

        movq  nb234_shiftvec(%rsp),%rax     ## eax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb234_iinr(%rsp),%rcx         ## ecx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb234_pos(%rbp),%rax      ## eax = base of pos[]  
        movl  %ebx,nb234_ii3(%rsp)

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
        movapd %xmm3,nb234_ixO(%rsp)
        movapd %xmm4,nb234_iyO(%rsp)
        movapd %xmm5,nb234_izO(%rsp)
        movapd %xmm6,nb234_ixH1(%rsp)
        movapd %xmm7,nb234_iyH1(%rsp)

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
        movapd %xmm6,nb234_izH1(%rsp)
        movapd %xmm0,nb234_ixH2(%rsp)
        movapd %xmm1,nb234_iyH2(%rsp)
        movapd %xmm2,nb234_izH2(%rsp)
        movapd %xmm3,nb234_ixM(%rsp)
        movapd %xmm4,nb234_iyM(%rsp)
        movapd %xmm5,nb234_izM(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb234_vctot(%rsp)
        movapd %xmm4,nb234_Vvdwtot(%rsp)
        movapd %xmm4,nb234_fixO(%rsp)
        movapd %xmm4,nb234_fiyO(%rsp)
        movapd %xmm4,nb234_fizO(%rsp)
        movapd %xmm4,nb234_fixH1(%rsp)
        movapd %xmm4,nb234_fiyH1(%rsp)
        movapd %xmm4,nb234_fizH1(%rsp)
        movapd %xmm4,nb234_fixH2(%rsp)
        movapd %xmm4,nb234_fiyH2(%rsp)
        movapd %xmm4,nb234_fizH2(%rsp)
        movapd %xmm4,nb234_fixM(%rsp)
        movapd %xmm4,nb234_fiyM(%rsp)
        movapd %xmm4,nb234_fizM(%rsp)

        movq  nb234_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb234_pos(%rbp),%rsi
        movq  nb234_faction(%rbp),%rdi
        movq  nb234_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb234_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb234_ninner(%rsp),%ecx
        movl  %ecx,nb234_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb234_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel234_x86_64_sse2.nb234_unroll_loop
        jmp   _nb_kernel234_x86_64_sse2.nb234_checksingle
_nb_kernel234_x86_64_sse2.nb234_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb234_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb234_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb234_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## load j O coordinates
    movlpd (%rsi,%rax,8),%xmm4
    movlpd 8(%rsi,%rax,8),%xmm5
    movlpd 16(%rsi,%rax,8),%xmm6
    movhpd (%rsi,%rbx,8),%xmm4
    movhpd 8(%rsi,%rbx,8),%xmm5
    movhpd 16(%rsi,%rbx,8),%xmm6

    ## xmm4 = Ox
    ## xmm5 = Oy
    ## xmm6 = Oz

    subpd nb234_ixO(%rsp),%xmm4
    subpd nb234_iyO(%rsp),%xmm5
    subpd nb234_izO(%rsp),%xmm6

    ## store dx/dy/dz
    movapd %xmm4,%xmm13
    movapd %xmm5,%xmm14
    movapd %xmm6,%xmm15

    ## square it
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        mulpd  %xmm6,%xmm6

        addpd  %xmm5,%xmm4
        addpd  %xmm6,%xmm4
    ## rsq in xmm4

        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb234_three(%rsp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb234_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb234_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb234_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm0,%xmm2       ## xmm0=iter2 of rinv (new lu) 

        mulpd %xmm2,%xmm4       ## xmm4=r 
        mulpd nb234_tsc(%rsp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1
    ## xmm1=eps 
    ## xmm2=rinv
    movapd %xmm4,%xmm3  ## eps
        pslld $3,%mm6           ## idx *= 8 

        movq nb234_VFtab(%rbp),%rsi
        movd %mm6,%r10d
        psrlq $32,%mm6
        movd %mm6,%r11d

    ## indices in r10, r11. Load dispersion and repulsion tables in parallel.
        movapd (%rsi,%r10,8),%xmm4          ## Y1d F1d  
        movapd (%rsi,%r11,8),%xmm0         ## Y2d F2d 
        movapd 32(%rsi,%r10,8),%xmm8        ## Y1r F1r  
        movapd 32(%rsi,%r11,8),%xmm3    ## Y2r F2r 
        movapd %xmm4,%xmm5
        movapd %xmm8,%xmm9
        unpcklpd %xmm0,%xmm4    ## Y1d Y2d 
        unpckhpd %xmm0,%xmm5    ## F1d F2d 
        unpcklpd %xmm3,%xmm8    ## Y1r Y2r 
        unpckhpd %xmm3,%xmm9    ## F1r F2r 

        movapd 16(%rsi,%r10,8),%xmm6        ## G1d H1d  
        movapd 16(%rsi,%r11,8),%xmm0            ## G2d H2d 
        movapd 48(%rsi,%r10,8),%xmm10           ## G1r H1r      
        movapd 48(%rsi,%r11,8),%xmm3        ## G2r H2r 
        movapd %xmm6,%xmm7
        movapd %xmm10,%xmm11
        unpcklpd %xmm0,%xmm6    ## G1d G2d 
        unpckhpd %xmm0,%xmm7    ## H1d H2d 
        unpcklpd %xmm3,%xmm10   ## G1r G2r 
        unpckhpd %xmm3,%xmm11   ## H1r H2r 
        ## tables ready, in xmm4-xmm7 and xmm8-xmm11

    mulpd  %xmm1,%xmm7   ## Heps
    mulpd  %xmm1,%xmm11
    mulpd  %xmm1,%xmm6  ## Geps
    mulpd  %xmm1,%xmm10
    mulpd  %xmm1,%xmm7  ## Heps2
    mulpd  %xmm1,%xmm11
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
    mulpd  %xmm1,%xmm5 ## eps*Fp
    mulpd  %xmm1,%xmm9
    addpd  %xmm4,%xmm5 ## VV
    addpd  %xmm8,%xmm9

    mulpd  nb234_c6(%rsp),%xmm5    ## VV*c6 = vnb6
    mulpd  nb234_c12(%rsp),%xmm9    ## VV*c12 = vnb12
    addpd  %xmm9,%xmm5
    addpd  nb234_Vvdwtot(%rsp),%xmm5
    movapd %xmm5,nb234_Vvdwtot(%rsp)

    mulpd  nb234_c6(%rsp),%xmm7     ## FF*c6 = fnb6
    mulpd  nb234_c12(%rsp),%xmm11     ## FF*c12  = fnb12
    addpd  %xmm11,%xmm7

    mulpd  nb234_tsc(%rsp),%xmm7
    mulpd  %xmm2,%xmm7
    xorpd  %xmm9,%xmm9

    subpd %xmm7,%xmm9
    mulpd %xmm9,%xmm13
    mulpd %xmm9,%xmm14
    mulpd %xmm9,%xmm15

    movapd nb234_fixO(%rsp),%xmm0
    movapd nb234_fiyO(%rsp),%xmm1
    movapd nb234_fizO(%rsp),%xmm2

    ## accumulate i forces
    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    movapd %xmm0,nb234_fixO(%rsp)
    movapd %xmm1,nb234_fiyO(%rsp)
    movapd %xmm2,nb234_fizO(%rsp)

        ## the fj's - start by accumulating forces from memory 
        movlpd (%rdi,%rax,8),%xmm3
        movlpd 8(%rdi,%rax,8),%xmm4
        movlpd 16(%rdi,%rax,8),%xmm5
        movhpd (%rdi,%rbx,8),%xmm3
        movhpd 8(%rdi,%rbx,8),%xmm4
        movhpd 16(%rdi,%rbx,8),%xmm5
        addpd %xmm13,%xmm3
        addpd %xmm14,%xmm4
        addpd %xmm15,%xmm5
        movlpd %xmm3,(%rdi,%rax,8)
        movlpd %xmm4,8(%rdi,%rax,8)
        movlpd %xmm5,16(%rdi,%rax,8)
        movhpd %xmm3,(%rdi,%rbx,8)
        movhpd %xmm4,8(%rdi,%rbx,8)
        movhpd %xmm5,16(%rdi,%rbx,8)
    ## done with OO interaction

    ## move j H1 coordinates to local temp variables 
    movq nb234_pos(%rbp),%rsi
    movlpd 24(%rsi,%rax,8),%xmm0
    movlpd 32(%rsi,%rax,8),%xmm1
    movlpd 40(%rsi,%rax,8),%xmm2
    movhpd 24(%rsi,%rbx,8),%xmm0
    movhpd 32(%rsi,%rbx,8),%xmm1
    movhpd 40(%rsi,%rbx,8),%xmm2

    ## xmm0 = H1x
    ## xmm1 = H1y
    ## xmm2 = H1z

    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb234_ixH1(%rsp),%xmm0
    subpd nb234_iyH1(%rsp),%xmm1
    subpd nb234_izH1(%rsp),%xmm2
    subpd nb234_ixH2(%rsp),%xmm3
    subpd nb234_iyH2(%rsp),%xmm4
    subpd nb234_izH2(%rsp),%xmm5
    subpd nb234_ixM(%rsp),%xmm6
    subpd nb234_iyM(%rsp),%xmm7
    subpd nb234_izM(%rsp),%xmm8

        movapd %xmm0,nb234_dxH1H1(%rsp)
        movapd %xmm1,nb234_dyH1H1(%rsp)
        movapd %xmm2,nb234_dzH1H1(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb234_dxH2H1(%rsp)
        movapd %xmm4,nb234_dyH2H1(%rsp)
        movapd %xmm5,nb234_dzH2H1(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb234_dxMH1(%rsp)
        movapd %xmm7,nb234_dyMH1(%rsp)
        movapd %xmm8,nb234_dzMH1(%rsp)
        mulpd  %xmm6,%xmm6
        mulpd  %xmm7,%xmm7
        mulpd  %xmm8,%xmm8
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
    addpd  %xmm7,%xmm6
    addpd  %xmm8,%xmm6

        ## start doing invsqrt for jH1 atoms
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

        movapd  nb234_three(%rsp),%xmm9
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

        movapd  nb234_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvH1H1 
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH2H1
    mulpd   %xmm15,%xmm11 ## first iteration for rinvMH1 

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb234_three(%rsp),%xmm1
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

        movapd  nb234_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvH1H1 
        mulpd   %xmm15,%xmm10 ##   rinvH2H1
    mulpd   %xmm15,%xmm11 ##   rinvMH1

        ## H1 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd %xmm9,%xmm1 ## copy of rinv
    movapd %xmm10,%xmm4
    movapd %xmm11,%xmm7
    movapd nb234_krf(%rsp),%xmm2
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
    movapd nb234_crf(%rsp),%xmm14
    subpd  %xmm14,%xmm2  ## rinv+krsq-crf
    subpd  %xmm14,%xmm5
    subpd  %xmm14,%xmm8
    movapd nb234_qqHH(%rsp),%xmm12
    movapd nb234_qqMH(%rsp),%xmm13
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
    addpd  nb234_vctot(%rsp),%xmm2
    addpd  %xmm8,%xmm5
    addpd  %xmm5,%xmm2
    movapd %xmm2,nb234_vctot(%rsp)

    mulpd  %xmm1,%xmm9  ## fscal
    mulpd  %xmm4,%xmm10
    mulpd  %xmm7,%xmm11

    ## move j H1 forces to xmm0-xmm2
    movq nb234_faction(%rbp),%rdi
        movlpd 24(%rdi,%rax,8),%xmm0
        movlpd 32(%rdi,%rax,8),%xmm1
        movlpd 40(%rdi,%rax,8),%xmm2
        movhpd 24(%rdi,%rbx,8),%xmm0
        movhpd 32(%rdi,%rbx,8),%xmm1
        movhpd 40(%rdi,%rbx,8),%xmm2

    movapd %xmm9,%xmm7
    movapd %xmm9,%xmm8
    movapd %xmm11,%xmm13
    movapd %xmm11,%xmm14
    movapd %xmm11,%xmm15
    movapd %xmm10,%xmm11
    movapd %xmm10,%xmm12

        mulpd nb234_dxH1H1(%rsp),%xmm7
        mulpd nb234_dyH1H1(%rsp),%xmm8
        mulpd nb234_dzH1H1(%rsp),%xmm9
        mulpd nb234_dxH2H1(%rsp),%xmm10
        mulpd nb234_dyH2H1(%rsp),%xmm11
        mulpd nb234_dzH2H1(%rsp),%xmm12
        mulpd nb234_dxMH1(%rsp),%xmm13
        mulpd nb234_dyMH1(%rsp),%xmm14
        mulpd nb234_dzMH1(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb234_fixH1(%rsp),%xmm7
    addpd nb234_fiyH1(%rsp),%xmm8
    addpd nb234_fizH1(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb234_fixH2(%rsp),%xmm10
    addpd nb234_fiyH2(%rsp),%xmm11
    addpd nb234_fizH2(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb234_fixM(%rsp),%xmm13
    addpd nb234_fiyM(%rsp),%xmm14
    addpd nb234_fizM(%rsp),%xmm15

    movapd %xmm7,nb234_fixH1(%rsp)
    movapd %xmm8,nb234_fiyH1(%rsp)
    movapd %xmm9,nb234_fizH1(%rsp)
    movapd %xmm10,nb234_fixH2(%rsp)
    movapd %xmm11,nb234_fiyH2(%rsp)
    movapd %xmm12,nb234_fizH2(%rsp)
    movapd %xmm13,nb234_fixM(%rsp)
    movapd %xmm14,nb234_fiyM(%rsp)
    movapd %xmm15,nb234_fizM(%rsp)

    ## store back j H1 forces from xmm0-xmm2
        movlpd %xmm0,24(%rdi,%rax,8)
        movlpd %xmm1,32(%rdi,%rax,8)
        movlpd %xmm2,40(%rdi,%rax,8)
        movhpd %xmm0,24(%rdi,%rbx,8)
        movhpd %xmm1,32(%rdi,%rbx,8)
        movhpd %xmm2,40(%rdi,%rbx,8)

        ## move j H2 coordinates to local temp variables 
    movq nb234_pos(%rbp),%rsi
    movlpd 48(%rsi,%rax,8),%xmm0
    movlpd 56(%rsi,%rax,8),%xmm1
    movlpd 64(%rsi,%rax,8),%xmm2
    movhpd 48(%rsi,%rbx,8),%xmm0
    movhpd 56(%rsi,%rbx,8),%xmm1
    movhpd 64(%rsi,%rbx,8),%xmm2

    ## xmm0 = H2x
    ## xmm1 = H2y
    ## xmm2 = H2z

    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb234_ixH1(%rsp),%xmm0
    subpd nb234_iyH1(%rsp),%xmm1
    subpd nb234_izH1(%rsp),%xmm2
    subpd nb234_ixH2(%rsp),%xmm3
    subpd nb234_iyH2(%rsp),%xmm4
    subpd nb234_izH2(%rsp),%xmm5
    subpd nb234_ixM(%rsp),%xmm6
    subpd nb234_iyM(%rsp),%xmm7
    subpd nb234_izM(%rsp),%xmm8

        movapd %xmm0,nb234_dxH1H2(%rsp)
        movapd %xmm1,nb234_dyH1H2(%rsp)
        movapd %xmm2,nb234_dzH1H2(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb234_dxH2H2(%rsp)
        movapd %xmm4,nb234_dyH2H2(%rsp)
        movapd %xmm5,nb234_dzH2H2(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb234_dxMH2(%rsp)
        movapd %xmm7,nb234_dyMH2(%rsp)
        movapd %xmm8,nb234_dzMH2(%rsp)
        mulpd  %xmm6,%xmm6
        mulpd  %xmm7,%xmm7
        mulpd  %xmm8,%xmm8
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
    addpd  %xmm7,%xmm6
    addpd  %xmm8,%xmm6

        ## start doing invsqrt for jH2 atoms
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

        movapd  nb234_three(%rsp),%xmm9
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

        movapd  nb234_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvH1H2 
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH2H2
    mulpd   %xmm15,%xmm11 ## first iteration for rinvMH2

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb234_three(%rsp),%xmm1
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

        movapd  nb234_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvH1H2
        mulpd   %xmm15,%xmm10 ##   rinvH2H2
    mulpd   %xmm15,%xmm11 ##   rinvMH2


        ## H2 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd %xmm9,%xmm1 ## copy of rinv
    movapd %xmm10,%xmm4
    movapd %xmm11,%xmm7
    movapd nb234_krf(%rsp),%xmm2
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
    movapd nb234_crf(%rsp),%xmm14
    subpd  %xmm14,%xmm2  ## rinv+krsq-crf
    subpd  %xmm14,%xmm5
    subpd  %xmm14,%xmm8
    movapd nb234_qqHH(%rsp),%xmm12
    movapd nb234_qqMH(%rsp),%xmm13
    mulpd  %xmm12,%xmm2 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulpd  %xmm12,%xmm5 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulpd  %xmm13,%xmm8 ## xmm6=voul=qq*(rinv+ krsq-crf)
    addpd  %xmm0,%xmm0 ## 2*krsq
    addpd  %xmm3,%xmm3
    addpd  %xmm6,%xmm6
    subpd  %xmm0,%xmm1 ## rinv-2*krsq
    subpd  %xmm3,%xmm4
    subpd  %xmm6,%xmm7
    mulpd  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulpd  %xmm12,%xmm4
    mulpd  %xmm13,%xmm7
    addpd  nb234_vctot(%rsp),%xmm2
    addpd  %xmm8,%xmm5
    addpd  %xmm5,%xmm2
    movapd  %xmm2,nb234_vctot(%rsp)

    mulpd  %xmm1,%xmm9  ## fscal
    mulpd  %xmm4,%xmm10
    mulpd  %xmm7,%xmm11

    ## move j H2 forces to xmm0-xmm2
    movq nb234_faction(%rbp),%rdi
        movlpd 48(%rdi,%rax,8),%xmm0
        movlpd 56(%rdi,%rax,8),%xmm1
        movlpd 64(%rdi,%rax,8),%xmm2
        movhpd 48(%rdi,%rbx,8),%xmm0
        movhpd 56(%rdi,%rbx,8),%xmm1
        movhpd 64(%rdi,%rbx,8),%xmm2

    movapd %xmm9,%xmm7
    movapd %xmm9,%xmm8
    movapd %xmm11,%xmm13
    movapd %xmm11,%xmm14
    movapd %xmm11,%xmm15
    movapd %xmm10,%xmm11
    movapd %xmm10,%xmm12

        mulpd nb234_dxH1H2(%rsp),%xmm7
        mulpd nb234_dyH1H2(%rsp),%xmm8
        mulpd nb234_dzH1H2(%rsp),%xmm9
        mulpd nb234_dxH2H2(%rsp),%xmm10
        mulpd nb234_dyH2H2(%rsp),%xmm11
        mulpd nb234_dzH2H2(%rsp),%xmm12
        mulpd nb234_dxMH2(%rsp),%xmm13
        mulpd nb234_dyMH2(%rsp),%xmm14
        mulpd nb234_dzMH2(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb234_fixH1(%rsp),%xmm7
    addpd nb234_fiyH1(%rsp),%xmm8
    addpd nb234_fizH1(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb234_fixH2(%rsp),%xmm10
    addpd nb234_fiyH2(%rsp),%xmm11
    addpd nb234_fizH2(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb234_fixM(%rsp),%xmm13
    addpd nb234_fiyM(%rsp),%xmm14
    addpd nb234_fizM(%rsp),%xmm15

    movapd %xmm7,nb234_fixH1(%rsp)
    movapd %xmm8,nb234_fiyH1(%rsp)
    movapd %xmm9,nb234_fizH1(%rsp)
    movapd %xmm10,nb234_fixH2(%rsp)
    movapd %xmm11,nb234_fiyH2(%rsp)
    movapd %xmm12,nb234_fizH2(%rsp)
    movapd %xmm13,nb234_fixM(%rsp)
    movapd %xmm14,nb234_fiyM(%rsp)
    movapd %xmm15,nb234_fizM(%rsp)

    ## store back j H2 forces from xmm0-xmm2
        movlpd %xmm0,48(%rdi,%rax,8)
        movlpd %xmm1,56(%rdi,%rax,8)
        movlpd %xmm2,64(%rdi,%rax,8)
        movhpd %xmm0,48(%rdi,%rbx,8)
        movhpd %xmm1,56(%rdi,%rbx,8)
        movhpd %xmm2,64(%rdi,%rbx,8)

        ## move j M coordinates to local temp variables 
    movq nb234_pos(%rbp),%rsi
    movlpd 72(%rsi,%rax,8),%xmm0
    movlpd 80(%rsi,%rax,8),%xmm1
    movlpd 88(%rsi,%rax,8),%xmm2
    movhpd 72(%rsi,%rbx,8),%xmm0
    movhpd 80(%rsi,%rbx,8),%xmm1
    movhpd 88(%rsi,%rbx,8),%xmm2

    ## xmm0 = Mx
    ## xmm1 = My
    ## xmm2 = Mz

    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb234_ixH1(%rsp),%xmm0
    subpd nb234_iyH1(%rsp),%xmm1
    subpd nb234_izH1(%rsp),%xmm2
    subpd nb234_ixH2(%rsp),%xmm3
    subpd nb234_iyH2(%rsp),%xmm4
    subpd nb234_izH2(%rsp),%xmm5
    subpd nb234_ixM(%rsp),%xmm6
    subpd nb234_iyM(%rsp),%xmm7
    subpd nb234_izM(%rsp),%xmm8

        movapd %xmm0,nb234_dxH1M(%rsp)
        movapd %xmm1,nb234_dyH1M(%rsp)
        movapd %xmm2,nb234_dzH1M(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb234_dxH2M(%rsp)
        movapd %xmm4,nb234_dyH2M(%rsp)
        movapd %xmm5,nb234_dzH2M(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb234_dxMM(%rsp)
        movapd %xmm7,nb234_dyMM(%rsp)
        movapd %xmm8,nb234_dzMM(%rsp)
        mulpd  %xmm6,%xmm6
        mulpd  %xmm7,%xmm7
        mulpd  %xmm8,%xmm8
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
    addpd  %xmm7,%xmm6
    addpd  %xmm8,%xmm6

        ## start doing invsqrt for jM atoms
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

        movapd  nb234_three(%rsp),%xmm9
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

        movapd  nb234_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvH1M 
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH2M
    mulpd   %xmm15,%xmm11 ## first iteration for rinvMM

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb234_three(%rsp),%xmm1
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

        movapd  nb234_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvH1M
        mulpd   %xmm15,%xmm10 ##   rinvH2M
    mulpd   %xmm15,%xmm11 ##   rinvMM


        ## M interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd %xmm9,%xmm1 ## copy of rinv
    movapd %xmm10,%xmm4
    movapd %xmm11,%xmm7
    movapd nb234_krf(%rsp),%xmm2
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
    movapd nb234_crf(%rsp),%xmm14
    subpd  %xmm14,%xmm2  ## rinv+krsq-crf
    subpd  %xmm14,%xmm5
    subpd  %xmm14,%xmm8
    movapd nb234_qqMH(%rsp),%xmm12
    movapd nb234_qqMM(%rsp),%xmm13
    mulpd  %xmm12,%xmm2 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulpd  %xmm12,%xmm5 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulpd  %xmm13,%xmm8 ## xmm6=voul=qq*(rinv+ krsq-crf)
    addpd  %xmm0,%xmm0 ## 2*krsq
    addpd  %xmm3,%xmm3
    addpd  %xmm6,%xmm6
    subpd  %xmm0,%xmm1 ## rinv-2*krsq
    subpd  %xmm3,%xmm4
    subpd  %xmm6,%xmm7
    mulpd  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulpd  %xmm12,%xmm4
    mulpd  %xmm13,%xmm7
    addpd  nb234_vctot(%rsp),%xmm2
    addpd  %xmm8,%xmm5
    addpd  %xmm5,%xmm2
    movapd  %xmm2,nb234_vctot(%rsp)

    mulpd  %xmm1,%xmm9  ## fscal
    mulpd  %xmm4,%xmm10
    mulpd  %xmm7,%xmm11

    ## move j M forces to xmm0-xmm2
    movq nb234_faction(%rbp),%rdi
        movlpd 72(%rdi,%rax,8),%xmm0
        movlpd 80(%rdi,%rax,8),%xmm1
        movlpd 88(%rdi,%rax,8),%xmm2
        movhpd 72(%rdi,%rbx,8),%xmm0
        movhpd 80(%rdi,%rbx,8),%xmm1
        movhpd 88(%rdi,%rbx,8),%xmm2

    movapd %xmm9,%xmm7
    movapd %xmm9,%xmm8
    movapd %xmm11,%xmm13
    movapd %xmm11,%xmm14
    movapd %xmm11,%xmm15
    movapd %xmm10,%xmm11
    movapd %xmm10,%xmm12

        mulpd nb234_dxH1M(%rsp),%xmm7
        mulpd nb234_dyH1M(%rsp),%xmm8
        mulpd nb234_dzH1M(%rsp),%xmm9
        mulpd nb234_dxH2M(%rsp),%xmm10
        mulpd nb234_dyH2M(%rsp),%xmm11
        mulpd nb234_dzH2M(%rsp),%xmm12
        mulpd nb234_dxMM(%rsp),%xmm13
        mulpd nb234_dyMM(%rsp),%xmm14
        mulpd nb234_dzMM(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb234_fixH1(%rsp),%xmm7
    addpd nb234_fiyH1(%rsp),%xmm8
    addpd nb234_fizH1(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb234_fixH2(%rsp),%xmm10
    addpd nb234_fiyH2(%rsp),%xmm11
    addpd nb234_fizH2(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb234_fixM(%rsp),%xmm13
    addpd nb234_fiyM(%rsp),%xmm14
    addpd nb234_fizM(%rsp),%xmm15

    movapd %xmm7,nb234_fixH1(%rsp)
    movapd %xmm8,nb234_fiyH1(%rsp)
    movapd %xmm9,nb234_fizH1(%rsp)
    movapd %xmm10,nb234_fixH2(%rsp)
    movapd %xmm11,nb234_fiyH2(%rsp)
    movapd %xmm12,nb234_fizH2(%rsp)
    movapd %xmm13,nb234_fixM(%rsp)
    movapd %xmm14,nb234_fiyM(%rsp)
    movapd %xmm15,nb234_fizM(%rsp)

    ## store back j M forces from xmm0-xmm2
        movlpd %xmm0,72(%rdi,%rax,8)
        movlpd %xmm1,80(%rdi,%rax,8)
        movlpd %xmm2,88(%rdi,%rax,8)
        movhpd %xmm0,72(%rdi,%rbx,8)
        movhpd %xmm1,80(%rdi,%rbx,8)
        movhpd %xmm2,88(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb234_innerk(%rsp)
        jl    _nb_kernel234_x86_64_sse2.nb234_checksingle
        jmp   _nb_kernel234_x86_64_sse2.nb234_unroll_loop
_nb_kernel234_x86_64_sse2.nb234_checksingle: 
        movl  nb234_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel234_x86_64_sse2.nb234_dosingle
        jmp   _nb_kernel234_x86_64_sse2.nb234_updateouterdata
_nb_kernel234_x86_64_sse2.nb234_dosingle: 
        movq  nb234_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb234_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## load j O coordinates
    movsd (%rsi,%rax,8),%xmm4
    movsd 8(%rsi,%rax,8),%xmm5
    movsd 16(%rsi,%rax,8),%xmm6

    ## xmm4 = Ox
    ## xmm5 = Oy
    ## xmm6 = Oz

    subsd nb234_ixO(%rsp),%xmm4
    subsd nb234_iyO(%rsp),%xmm5
    subsd nb234_izO(%rsp),%xmm6

    ## store dx/dy/dz
    movapd %xmm4,%xmm13
    movapd %xmm5,%xmm14
    movapd %xmm6,%xmm15

    ## square it
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        mulsd  %xmm6,%xmm6

        addsd  %xmm5,%xmm4
        addsd  %xmm6,%xmm4
    ## rsq in xmm4

        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb234_three(%rsp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb234_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb234_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb234_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm0,%xmm2       ## xmm0=iter2 of rinv (new lu) 

        mulsd %xmm2,%xmm4       ## xmm4=r 
        mulsd nb234_tsc(%rsp),%xmm4

        cvttsd2si %xmm4,%r10d   ## mm6 = lu idx 
        cvtsi2sd %r10d,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1
    ## xmm1=eps 
    ## xmm2=rinv
    movapd %xmm4,%xmm3  ## eps
        shll $3,%r10d           ## idx *= 8 

        movq nb234_VFtab(%rbp),%rsi

    ## indices in r10, r11. Load dispersion and repulsion tables in parallel.
        movapd (%rsi,%r10,8),%xmm4          ## Y1d F1d  
        movapd 32(%rsi,%r10,8),%xmm8        ## Y1r F1r  
        movhlps %xmm4,%xmm5
        movhlps %xmm8,%xmm9

        movapd 16(%rsi,%r10,8),%xmm6        ## G1d H1d  
        movapd 48(%rsi,%r10,8),%xmm10           ## G1r H1r      
        movhlps %xmm6,%xmm7
        movhlps %xmm10,%xmm11
        ## tables ready, in xmm4-xmm7 and xmm8-xmm11

    mulsd  %xmm1,%xmm7   ## Heps
    mulsd  %xmm1,%xmm11
    mulsd  %xmm1,%xmm6  ## Geps
    mulsd  %xmm1,%xmm10
    mulsd  %xmm1,%xmm7  ## Heps2
    mulsd  %xmm1,%xmm11
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
    mulsd  %xmm1,%xmm5 ## eps*Fp
    mulsd  %xmm1,%xmm9
    addsd  %xmm4,%xmm5 ## VV
    addsd  %xmm8,%xmm9

    mulsd  nb234_c6(%rsp),%xmm5    ## VV*c6 = vnb6
    mulsd  nb234_c12(%rsp),%xmm9    ## VV*c12 = vnb12
    addsd  %xmm9,%xmm5
    addsd  nb234_Vvdwtot(%rsp),%xmm5
    movsd %xmm5,nb234_Vvdwtot(%rsp)

    mulsd  nb234_c6(%rsp),%xmm7     ## FF*c6 = fnb6
    mulsd  nb234_c12(%rsp),%xmm11     ## FF*c12  = fnb12
    addsd  %xmm11,%xmm7

    mulsd  nb234_tsc(%rsp),%xmm7
    mulsd  %xmm2,%xmm7
    xorpd  %xmm9,%xmm9

    subpd  %xmm7,%xmm9
    mulsd %xmm9,%xmm13
    mulsd %xmm9,%xmm14
    mulsd %xmm9,%xmm15

    movq nb234_faction(%rbp),%rdi

    movapd nb234_fixO(%rsp),%xmm0
    movapd nb234_fiyO(%rsp),%xmm1
    movapd nb234_fizO(%rsp),%xmm2

    ## accumulate i forces
    addsd %xmm13,%xmm0
    addsd %xmm14,%xmm1
    addsd %xmm15,%xmm2
    movsd %xmm0,nb234_fixO(%rsp)
    movsd %xmm1,nb234_fiyO(%rsp)
    movsd %xmm2,nb234_fizO(%rsp)

        ## the fj's - start by accumulating forces from memory 
        addsd (%rdi,%rax,8),%xmm13
        addsd 8(%rdi,%rax,8),%xmm14
        addsd 16(%rdi,%rax,8),%xmm15
        movsd %xmm13,(%rdi,%rax,8)
        movsd %xmm14,8(%rdi,%rax,8)
        movsd %xmm15,16(%rdi,%rax,8)
    ## done with OO interaction

        ## move j H1 coordinates to local temp variables 
    movq nb234_pos(%rbp),%rsi
    movsd 24(%rsi,%rax,8),%xmm0
    movsd 32(%rsi,%rax,8),%xmm1
    movsd 40(%rsi,%rax,8),%xmm2

    ## xmm0 = H1x
    ## xmm1 = H1y
    ## xmm2 = H1z

    movsd %xmm0,%xmm3
    movsd %xmm1,%xmm4
    movsd %xmm2,%xmm5
    movsd %xmm0,%xmm6
    movsd %xmm1,%xmm7
    movsd %xmm2,%xmm8

    subsd nb234_ixH1(%rsp),%xmm0
    subsd nb234_iyH1(%rsp),%xmm1
    subsd nb234_izH1(%rsp),%xmm2
    subsd nb234_ixH2(%rsp),%xmm3
    subsd nb234_iyH2(%rsp),%xmm4
    subsd nb234_izH2(%rsp),%xmm5
    subsd nb234_ixM(%rsp),%xmm6
    subsd nb234_iyM(%rsp),%xmm7
    subsd nb234_izM(%rsp),%xmm8

        movsd %xmm0,nb234_dxH1H1(%rsp)
        movsd %xmm1,nb234_dyH1H1(%rsp)
        movsd %xmm2,nb234_dzH1H1(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb234_dxH2H1(%rsp)
        movsd %xmm4,nb234_dyH2H1(%rsp)
        movsd %xmm5,nb234_dzH2H1(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb234_dxMH1(%rsp)
        movsd %xmm7,nb234_dyMH1(%rsp)
        movsd %xmm8,nb234_dzMH1(%rsp)
        mulsd  %xmm6,%xmm6
        mulsd  %xmm7,%xmm7
        mulsd  %xmm8,%xmm8
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
    addsd  %xmm7,%xmm6
    addsd  %xmm8,%xmm6

        ## start doing invsqrt for jH1 atoms
    cvtsd2ss %xmm0,%xmm1
    cvtsd2ss %xmm3,%xmm4
    cvtsd2ss %xmm6,%xmm7
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm4,%xmm4
    rsqrtss %xmm7,%xmm7
    cvtss2sd %xmm1,%xmm1
    cvtss2sd %xmm4,%xmm4
    cvtss2sd %xmm7,%xmm7

        movsd  %xmm1,%xmm2
        movsd  %xmm4,%xmm5
    movsd  %xmm7,%xmm8

        mulsd   %xmm1,%xmm1 ## lu*lu
        mulsd   %xmm4,%xmm4 ## lu*lu
    mulsd   %xmm7,%xmm7 ## lu*lu

        movsd  nb234_three(%rsp),%xmm9
        movsd  %xmm9,%xmm10
    movsd  %xmm9,%xmm11

        mulsd   %xmm0,%xmm1 ## rsq*lu*lu
        mulsd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm7 ## rsq*lu*lu

        subsd   %xmm1,%xmm9
        subsd   %xmm4,%xmm10
    subsd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulsd   %xmm2,%xmm9
        mulsd   %xmm5,%xmm10
    mulsd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb234_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ## first iteration for rinvH1H1 
        mulsd   %xmm15,%xmm10 ## first iteration for rinvH2H1
    mulsd   %xmm15,%xmm11 ## first iteration for rinvMH1 

    ## second iteration step    
        movsd  %xmm9,%xmm2
        movsd  %xmm10,%xmm5
    movsd  %xmm11,%xmm8

        mulsd   %xmm2,%xmm2 ## lu*lu
        mulsd   %xmm5,%xmm5 ## lu*lu
    mulsd   %xmm8,%xmm8 ## lu*lu

        movsd  nb234_three(%rsp),%xmm1
        movsd  %xmm1,%xmm4
    movsd  %xmm1,%xmm7

        mulsd   %xmm0,%xmm2 ## rsq*lu*lu
        mulsd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm8 ## rsq*lu*lu

        subsd   %xmm2,%xmm1
        subsd   %xmm5,%xmm4
    subsd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulsd   %xmm1,%xmm9
        mulsd   %xmm4,%xmm10
    mulsd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb234_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvH1H1 
        mulsd   %xmm15,%xmm10 ##   rinvH2H1
    mulsd   %xmm15,%xmm11 ##   rinvMH1

        ## H1 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movsd %xmm9,%xmm1 ## copy of rinv
    movsd %xmm10,%xmm4
    movsd %xmm11,%xmm7
    movsd nb234_krf(%rsp),%xmm2
    mulsd  %xmm9,%xmm9  ## rinvsq
    mulsd  %xmm10,%xmm10
    mulsd  %xmm11,%xmm11
    mulsd  %xmm2,%xmm0 ## k*rsq
    mulsd  %xmm2,%xmm3
    mulsd  %xmm2,%xmm6
    movsd  %xmm0,%xmm2 ## copy of k*rsq
    movsd  %xmm3,%xmm5
    movsd  %xmm6,%xmm8
    addsd  %xmm1,%xmm2 ## rinv+krsq
    addsd  %xmm4,%xmm5
    addsd  %xmm7,%xmm8
    movsd  nb234_crf(%rsp),%xmm14
    subsd  %xmm14,%xmm2  ## rinv+krsq-crf
    subsd  %xmm14,%xmm5
    subsd  %xmm14,%xmm8
    movsd  nb234_qqHH(%rsp),%xmm12
    movsd  nb234_qqMH(%rsp),%xmm13
    mulsd  %xmm12,%xmm2 ## voul=qq*(rinv+ krsq-crf)
    mulsd  %xmm12,%xmm5 ## voul=qq*(rinv+ krsq-crf)
    mulsd  %xmm13,%xmm8 ## voul=qq*(rinv+ krsq-crf)
    addsd  %xmm0,%xmm0 ## 2*krsq
    addsd  %xmm3,%xmm3
    addsd  %xmm6,%xmm6
    subsd  %xmm0,%xmm1 ## rinv-2*krsq
    subsd  %xmm3,%xmm4
    subsd  %xmm6,%xmm7
    mulsd  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulsd  %xmm12,%xmm4
    mulsd  %xmm13,%xmm7
    addsd  nb234_vctot(%rsp),%xmm2
    addsd  %xmm8,%xmm5
    addsd  %xmm5,%xmm2
    movsd  %xmm2,nb234_vctot(%rsp)

    mulsd  %xmm1,%xmm9  ## fscal
    mulsd  %xmm4,%xmm10
    mulsd  %xmm7,%xmm11

    ## move j H1 forces to xmm0-xmm2
    movq nb234_faction(%rbp),%rdi
        movsd 24(%rdi,%rax,8),%xmm0
        movsd 32(%rdi,%rax,8),%xmm1
        movsd 40(%rdi,%rax,8),%xmm2

    movsd %xmm9,%xmm7
    movsd %xmm9,%xmm8
    movsd %xmm11,%xmm13
    movsd %xmm11,%xmm14
    movsd %xmm11,%xmm15
    movsd %xmm10,%xmm11
    movsd %xmm10,%xmm12

        mulsd nb234_dxH1H1(%rsp),%xmm7
        mulsd nb234_dyH1H1(%rsp),%xmm8
        mulsd nb234_dzH1H1(%rsp),%xmm9
        mulsd nb234_dxH2H1(%rsp),%xmm10
        mulsd nb234_dyH2H1(%rsp),%xmm11
        mulsd nb234_dzH2H1(%rsp),%xmm12
        mulsd nb234_dxMH1(%rsp),%xmm13
        mulsd nb234_dyMH1(%rsp),%xmm14
        mulsd nb234_dzMH1(%rsp),%xmm15

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb234_fixH1(%rsp),%xmm7
    addsd nb234_fiyH1(%rsp),%xmm8
    addsd nb234_fizH1(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb234_fixH2(%rsp),%xmm10
    addsd nb234_fiyH2(%rsp),%xmm11
    addsd nb234_fizH2(%rsp),%xmm12

    addsd %xmm13,%xmm0
    addsd %xmm14,%xmm1
    addsd %xmm15,%xmm2
    addsd nb234_fixM(%rsp),%xmm13
    addsd nb234_fiyM(%rsp),%xmm14
    addsd nb234_fizM(%rsp),%xmm15

    movsd %xmm7,nb234_fixH1(%rsp)
    movsd %xmm8,nb234_fiyH1(%rsp)
    movsd %xmm9,nb234_fizH1(%rsp)
    movsd %xmm10,nb234_fixH2(%rsp)
    movsd %xmm11,nb234_fiyH2(%rsp)
    movsd %xmm12,nb234_fizH2(%rsp)
    movsd %xmm13,nb234_fixM(%rsp)
    movsd %xmm14,nb234_fiyM(%rsp)
    movsd %xmm15,nb234_fizM(%rsp)

    ## store back j H1 forces from xmm0-xmm2
        movsd %xmm0,24(%rdi,%rax,8)
        movsd %xmm1,32(%rdi,%rax,8)
        movsd %xmm2,40(%rdi,%rax,8)

        ## move j H2 coordinates to local temp variables 
    movq nb234_pos(%rbp),%rsi
    movsd 48(%rsi,%rax,8),%xmm0
    movsd 56(%rsi,%rax,8),%xmm1
    movsd 64(%rsi,%rax,8),%xmm2

    ## xmm0 = H2x
    ## xmm1 = H2y
    ## xmm2 = H2z

    movsd %xmm0,%xmm3
    movsd %xmm1,%xmm4
    movsd %xmm2,%xmm5
    movsd %xmm0,%xmm6
    movsd %xmm1,%xmm7
    movsd %xmm2,%xmm8

    subsd nb234_ixH1(%rsp),%xmm0
    subsd nb234_iyH1(%rsp),%xmm1
    subsd nb234_izH1(%rsp),%xmm2
    subsd nb234_ixH2(%rsp),%xmm3
    subsd nb234_iyH2(%rsp),%xmm4
    subsd nb234_izH2(%rsp),%xmm5
    subsd nb234_ixM(%rsp),%xmm6
    subsd nb234_iyM(%rsp),%xmm7
    subsd nb234_izM(%rsp),%xmm8

        movsd %xmm0,nb234_dxH1H2(%rsp)
        movsd %xmm1,nb234_dyH1H2(%rsp)
        movsd %xmm2,nb234_dzH1H2(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb234_dxH2H2(%rsp)
        movsd %xmm4,nb234_dyH2H2(%rsp)
        movsd %xmm5,nb234_dzH2H2(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb234_dxMH2(%rsp)
        movsd %xmm7,nb234_dyMH2(%rsp)
        movsd %xmm8,nb234_dzMH2(%rsp)
        mulsd  %xmm6,%xmm6
        mulsd  %xmm7,%xmm7
        mulsd  %xmm8,%xmm8
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
    addsd  %xmm7,%xmm6
    addsd  %xmm8,%xmm6

        ## start doing invsqrt for jH2 atoms
    cvtsd2ss %xmm0,%xmm1
    cvtsd2ss %xmm3,%xmm4
    cvtsd2ss %xmm6,%xmm7
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm4,%xmm4
    rsqrtss %xmm7,%xmm7
    cvtss2sd %xmm1,%xmm1
    cvtss2sd %xmm4,%xmm4
    cvtss2sd %xmm7,%xmm7

        movsd  %xmm1,%xmm2
        movsd  %xmm4,%xmm5
    movsd  %xmm7,%xmm8

        mulsd   %xmm1,%xmm1 ## lu*lu
        mulsd   %xmm4,%xmm4 ## lu*lu
    mulsd   %xmm7,%xmm7 ## lu*lu

        movsd  nb234_three(%rsp),%xmm9
        movsd  %xmm9,%xmm10
    movsd  %xmm9,%xmm11

        mulsd   %xmm0,%xmm1 ## rsq*lu*lu
        mulsd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm7 ## rsq*lu*lu

        subsd   %xmm1,%xmm9
        subsd   %xmm4,%xmm10
    subsd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulsd   %xmm2,%xmm9
        mulsd   %xmm5,%xmm10
    mulsd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb234_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ## first iteration for rinvH1H2 
        mulsd   %xmm15,%xmm10 ## first iteration for rinvH2H2
    mulsd   %xmm15,%xmm11 ## first iteration for rinvMH2

    ## second iteration step    
        movsd  %xmm9,%xmm2
        movsd  %xmm10,%xmm5
    movsd  %xmm11,%xmm8

        mulsd   %xmm2,%xmm2 ## lu*lu
        mulsd   %xmm5,%xmm5 ## lu*lu
    mulsd   %xmm8,%xmm8 ## lu*lu

        movsd  nb234_three(%rsp),%xmm1
        movsd  %xmm1,%xmm4
    movsd  %xmm1,%xmm7

        mulsd   %xmm0,%xmm2 ## rsq*lu*lu
        mulsd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm8 ## rsq*lu*lu

        subsd   %xmm2,%xmm1
        subsd   %xmm5,%xmm4
    subsd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulsd   %xmm1,%xmm9
        mulsd   %xmm4,%xmm10
    mulsd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb234_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvH1H2
        mulsd   %xmm15,%xmm10 ##   rinvH2H2
    mulsd   %xmm15,%xmm11 ##   rinvMH2

        ## H2 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movsd %xmm9,%xmm1 ## copy of rinv
    movsd %xmm10,%xmm4
    movsd %xmm11,%xmm7
    movsd nb234_krf(%rsp),%xmm2
    mulsd  %xmm9,%xmm9  ## rinvsq
    mulsd  %xmm10,%xmm10
    mulsd  %xmm11,%xmm11
    mulsd  %xmm2,%xmm0 ## k*rsq
    mulsd  %xmm2,%xmm3
    mulsd  %xmm2,%xmm6
    movsd %xmm0,%xmm2 ## copy of k*rsq
    movsd %xmm3,%xmm5
    movsd %xmm6,%xmm8
    addsd  %xmm1,%xmm2 ## rinv+krsq
    addsd  %xmm4,%xmm5
    addsd  %xmm7,%xmm8
    movsd nb234_crf(%rsp),%xmm14
    subsd  %xmm14,%xmm2  ## rinv+krsq-crf
    subsd  %xmm14,%xmm5
    subsd  %xmm14,%xmm8
    movsd nb234_qqHH(%rsp),%xmm12
    movsd nb234_qqMH(%rsp),%xmm13
    mulsd  %xmm12,%xmm2 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulsd  %xmm12,%xmm5 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulsd  %xmm13,%xmm8 ## xmm6=voul=qq*(rinv+ krsq-crf)
    addsd  %xmm0,%xmm0 ## 2*krsq
    addsd  %xmm3,%xmm3
    addsd  %xmm6,%xmm6
    subsd  %xmm0,%xmm1 ## rinv-2*krsq
    subsd  %xmm3,%xmm4
    subsd  %xmm6,%xmm7
    mulsd  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulsd  %xmm12,%xmm4
    mulsd  %xmm13,%xmm7
    addsd  nb234_vctot(%rsp),%xmm2
    addsd  %xmm8,%xmm5
    addsd  %xmm5,%xmm2
    movsd  %xmm2,nb234_vctot(%rsp)

    mulsd  %xmm1,%xmm9  ## fscal
    mulsd  %xmm4,%xmm10
    mulsd  %xmm7,%xmm11

    ## move j H2 forces to xmm0-xmm2
    movq nb234_faction(%rbp),%rdi
        movsd 48(%rdi,%rax,8),%xmm0
        movsd 56(%rdi,%rax,8),%xmm1
        movsd 64(%rdi,%rax,8),%xmm2

    movsd %xmm9,%xmm7
    movsd %xmm9,%xmm8
    movsd %xmm11,%xmm13
    movsd %xmm11,%xmm14
    movsd %xmm11,%xmm15
    movsd %xmm10,%xmm11
    movsd %xmm10,%xmm12

        mulsd nb234_dxH1H2(%rsp),%xmm7
        mulsd nb234_dyH1H2(%rsp),%xmm8
        mulsd nb234_dzH1H2(%rsp),%xmm9
        mulsd nb234_dxH2H2(%rsp),%xmm10
        mulsd nb234_dyH2H2(%rsp),%xmm11
        mulsd nb234_dzH2H2(%rsp),%xmm12
        mulsd nb234_dxMH2(%rsp),%xmm13
        mulsd nb234_dyMH2(%rsp),%xmm14
        mulsd nb234_dzMH2(%rsp),%xmm15

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb234_fixH1(%rsp),%xmm7
    addsd nb234_fiyH1(%rsp),%xmm8
    addsd nb234_fizH1(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb234_fixH2(%rsp),%xmm10
    addsd nb234_fiyH2(%rsp),%xmm11
    addsd nb234_fizH2(%rsp),%xmm12

    addsd %xmm13,%xmm0
    addsd %xmm14,%xmm1
    addsd %xmm15,%xmm2
    addsd nb234_fixM(%rsp),%xmm13
    addsd nb234_fiyM(%rsp),%xmm14
    addsd nb234_fizM(%rsp),%xmm15

    movsd %xmm7,nb234_fixH1(%rsp)
    movsd %xmm8,nb234_fiyH1(%rsp)
    movsd %xmm9,nb234_fizH1(%rsp)
    movsd %xmm10,nb234_fixH2(%rsp)
    movsd %xmm11,nb234_fiyH2(%rsp)
    movsd %xmm12,nb234_fizH2(%rsp)
    movsd %xmm13,nb234_fixM(%rsp)
    movsd %xmm14,nb234_fiyM(%rsp)
    movsd %xmm15,nb234_fizM(%rsp)

    ## store back j H2 forces from xmm0-xmm2
        movsd %xmm0,48(%rdi,%rax,8)
        movsd %xmm1,56(%rdi,%rax,8)
        movsd %xmm2,64(%rdi,%rax,8)

        ## move j M coordinates to local temp variables 
    movq nb234_pos(%rbp),%rsi
    movsd 72(%rsi,%rax,8),%xmm0
    movsd 80(%rsi,%rax,8),%xmm1
    movsd 88(%rsi,%rax,8),%xmm2

    ## xmm0 = Mx
    ## xmm1 = My
    ## xmm2 = Mz

    movsd %xmm0,%xmm3
    movsd %xmm1,%xmm4
    movsd %xmm2,%xmm5
    movsd %xmm0,%xmm6
    movsd %xmm1,%xmm7
    movsd %xmm2,%xmm8

    subsd nb234_ixH1(%rsp),%xmm0
    subsd nb234_iyH1(%rsp),%xmm1
    subsd nb234_izH1(%rsp),%xmm2
    subsd nb234_ixH2(%rsp),%xmm3
    subsd nb234_iyH2(%rsp),%xmm4
    subsd nb234_izH2(%rsp),%xmm5
    subsd nb234_ixM(%rsp),%xmm6
    subsd nb234_iyM(%rsp),%xmm7
    subsd nb234_izM(%rsp),%xmm8

        movsd %xmm0,nb234_dxH1M(%rsp)
        movsd %xmm1,nb234_dyH1M(%rsp)
        movsd %xmm2,nb234_dzH1M(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb234_dxH2M(%rsp)
        movsd %xmm4,nb234_dyH2M(%rsp)
        movsd %xmm5,nb234_dzH2M(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb234_dxMM(%rsp)
        movsd %xmm7,nb234_dyMM(%rsp)
        movsd %xmm8,nb234_dzMM(%rsp)
        mulsd  %xmm6,%xmm6
        mulsd  %xmm7,%xmm7
        mulsd  %xmm8,%xmm8
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
    addsd  %xmm7,%xmm6
    addsd  %xmm8,%xmm6

        ## start doing invsqrt for jM atoms
    cvtsd2ss %xmm0,%xmm1
    cvtsd2ss %xmm3,%xmm4
    cvtsd2ss %xmm6,%xmm7
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm4,%xmm4
    rsqrtss %xmm7,%xmm7
    cvtss2sd %xmm1,%xmm1
    cvtss2sd %xmm4,%xmm4
    cvtss2sd %xmm7,%xmm7

        movsd  %xmm1,%xmm2
        movsd  %xmm4,%xmm5
    movsd  %xmm7,%xmm8

        mulsd   %xmm1,%xmm1 ## lu*lu
        mulsd   %xmm4,%xmm4 ## lu*lu
    mulsd   %xmm7,%xmm7 ## lu*lu

        movsd  nb234_three(%rsp),%xmm9
        movsd  %xmm9,%xmm10
    movsd  %xmm9,%xmm11

        mulsd   %xmm0,%xmm1 ## rsq*lu*lu
        mulsd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm7 ## rsq*lu*lu

        subsd   %xmm1,%xmm9
        subsd   %xmm4,%xmm10
    subsd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulsd   %xmm2,%xmm9
        mulsd   %xmm5,%xmm10
    mulsd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb234_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ## first iteration for rinvH1M 
        mulsd   %xmm15,%xmm10 ## first iteration for rinvH2M
    mulsd   %xmm15,%xmm11 ## first iteration for rinvMM

    ## second iteration step    
        movsd  %xmm9,%xmm2
        movsd  %xmm10,%xmm5
    movsd  %xmm11,%xmm8

        mulsd   %xmm2,%xmm2 ## lu*lu
        mulsd   %xmm5,%xmm5 ## lu*lu
    mulsd   %xmm8,%xmm8 ## lu*lu

        movsd  nb234_three(%rsp),%xmm1
        movsd  %xmm1,%xmm4
    movsd  %xmm1,%xmm7

        mulsd   %xmm0,%xmm2 ## rsq*lu*lu
        mulsd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm8 ## rsq*lu*lu

        subsd   %xmm2,%xmm1
        subsd   %xmm5,%xmm4
    subsd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulsd   %xmm1,%xmm9
        mulsd   %xmm4,%xmm10
    mulsd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb234_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvH1M
        mulsd   %xmm15,%xmm10 ##   rinvH2M
    mulsd   %xmm15,%xmm11 ##   rinvMM

        ## M interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movsd %xmm9,%xmm1 ## copy of rinv
    movsd %xmm10,%xmm4
    movsd %xmm11,%xmm7
    movsd nb234_krf(%rsp),%xmm2
    mulsd  %xmm9,%xmm9  ## rinvsq
    mulsd  %xmm10,%xmm10
    mulsd  %xmm11,%xmm11
    mulsd  %xmm2,%xmm0 ## k*rsq
    mulsd  %xmm2,%xmm3
    mulsd  %xmm2,%xmm6
    movsd %xmm0,%xmm2 ## copy of k*rsq
    movsd %xmm3,%xmm5
    movsd %xmm6,%xmm8
    addsd  %xmm1,%xmm2 ## rinv+krsq
    addsd  %xmm4,%xmm5
    addsd  %xmm7,%xmm8
    movsd nb234_crf(%rsp),%xmm14
    subsd  %xmm14,%xmm2  ## rinv+krsq-crf
    subsd  %xmm14,%xmm5
    subsd  %xmm14,%xmm8
    movsd nb234_qqMH(%rsp),%xmm12
    movsd nb234_qqMM(%rsp),%xmm13
    mulsd  %xmm12,%xmm2 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulsd  %xmm12,%xmm5 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulsd  %xmm13,%xmm8 ## xmm6=voul=qq*(rinv+ krsq-crf)
    addsd  %xmm0,%xmm0 ## 2*krsq
    addsd  %xmm3,%xmm3
    addsd  %xmm6,%xmm6
    subsd  %xmm0,%xmm1 ## rinv-2*krsq
    subsd  %xmm3,%xmm4
    subsd  %xmm6,%xmm7
    mulsd  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulsd  %xmm12,%xmm4
    mulsd  %xmm13,%xmm7
    addsd  nb234_vctot(%rsp),%xmm2
    addsd  %xmm8,%xmm5
    addsd  %xmm5,%xmm2
    movsd  %xmm2,nb234_vctot(%rsp)

    mulsd  %xmm1,%xmm9  ## fscal
    mulsd  %xmm4,%xmm10
    mulsd  %xmm7,%xmm11

    ## move j M forces to xmm0-xmm2
    movq nb234_faction(%rbp),%rdi
        movsd 72(%rdi,%rax,8),%xmm0
        movsd 80(%rdi,%rax,8),%xmm1
        movsd 88(%rdi,%rax,8),%xmm2

    movsd %xmm9,%xmm7
    movsd %xmm9,%xmm8
    movsd %xmm11,%xmm13
    movsd %xmm11,%xmm14
    movsd %xmm11,%xmm15
    movsd %xmm10,%xmm11
    movsd %xmm10,%xmm12

        mulsd nb234_dxH1M(%rsp),%xmm7
        mulsd nb234_dyH1M(%rsp),%xmm8
        mulsd nb234_dzH1M(%rsp),%xmm9
        mulsd nb234_dxH2M(%rsp),%xmm10
        mulsd nb234_dyH2M(%rsp),%xmm11
        mulsd nb234_dzH2M(%rsp),%xmm12
        mulsd nb234_dxMM(%rsp),%xmm13
        mulsd nb234_dyMM(%rsp),%xmm14
        mulsd nb234_dzMM(%rsp),%xmm15

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb234_fixH1(%rsp),%xmm7
    addsd nb234_fiyH1(%rsp),%xmm8
    addsd nb234_fizH1(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb234_fixH2(%rsp),%xmm10
    addsd nb234_fiyH2(%rsp),%xmm11
    addsd nb234_fizH2(%rsp),%xmm12

    addsd %xmm13,%xmm0
    addsd %xmm14,%xmm1
    addsd %xmm15,%xmm2
    addsd nb234_fixM(%rsp),%xmm13
    addsd nb234_fiyM(%rsp),%xmm14
    addsd nb234_fizM(%rsp),%xmm15

    movsd %xmm7,nb234_fixH1(%rsp)
    movsd %xmm8,nb234_fiyH1(%rsp)
    movsd %xmm9,nb234_fizH1(%rsp)
    movsd %xmm10,nb234_fixH2(%rsp)
    movsd %xmm11,nb234_fiyH2(%rsp)
    movsd %xmm12,nb234_fizH2(%rsp)
    movsd %xmm13,nb234_fixM(%rsp)
    movsd %xmm14,nb234_fiyM(%rsp)
    movsd %xmm15,nb234_fizM(%rsp)

    ## store back j M forces from xmm0-xmm2
        movsd %xmm0,72(%rdi,%rax,8)
        movsd %xmm1,80(%rdi,%rax,8)
        movsd %xmm2,88(%rdi,%rax,8)

_nb_kernel234_x86_64_sse2.nb234_updateouterdata: 
        movl  nb234_ii3(%rsp),%ecx
        movq  nb234_faction(%rbp),%rdi
        movq  nb234_fshift(%rbp),%rsi
        movl  nb234_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb234_fixO(%rsp),%xmm0
        movapd nb234_fiyO(%rsp),%xmm1
        movapd nb234_fizO(%rsp),%xmm2

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
        movapd nb234_fixH1(%rsp),%xmm0
        movapd nb234_fiyH1(%rsp),%xmm1
        movapd nb234_fizH1(%rsp),%xmm2

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
        movapd nb234_fixH2(%rsp),%xmm0
        movapd nb234_fiyH2(%rsp),%xmm1
        movapd nb234_fizH2(%rsp),%xmm2

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
        movapd nb234_fixM(%rsp),%xmm0
        movapd nb234_fiyM(%rsp),%xmm1
        movapd nb234_fizM(%rsp),%xmm2

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
        movl nb234_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb234_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb234_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb234_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb234_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb234_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb234_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel234_x86_64_sse2.nb234_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb234_n(%rsp)
        jmp _nb_kernel234_x86_64_sse2.nb234_outer
_nb_kernel234_x86_64_sse2.nb234_outerend: 
        ## check if more outer neighborlists remain
        movl  nb234_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel234_x86_64_sse2.nb234_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel234_x86_64_sse2.nb234_threadloop
_nb_kernel234_x86_64_sse2.nb234_end: 
        movl nb234_nouter(%rsp),%eax
        movl nb234_ninner(%rsp),%ebx
        movq nb234_outeriter(%rbp),%rcx
        movq nb234_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1912,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret




.globl nb_kernel234nf_x86_64_sse2
.globl _nb_kernel234nf_x86_64_sse2
nb_kernel234nf_x86_64_sse2:     
_nb_kernel234nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb234nf_fshift, 16
.set nb234nf_gid, 24
.set nb234nf_pos, 32
.set nb234nf_faction, 40
.set nb234nf_charge, 48
.set nb234nf_p_facel, 56
.set nb234nf_argkrf, 64
.set nb234nf_argcrf, 72
.set nb234nf_Vc, 80
.set nb234nf_type, 88
.set nb234nf_p_ntype, 96
.set nb234nf_vdwparam, 104
.set nb234nf_Vvdw, 112
.set nb234nf_p_tabscale, 120
.set nb234nf_VFtab, 128
.set nb234nf_invsqrta, 136
.set nb234nf_dvda, 144
.set nb234nf_p_gbtabscale, 152
.set nb234nf_GBtab, 160
.set nb234nf_p_nthreads, 168
.set nb234nf_count, 176
.set nb234nf_mtx, 184
.set nb234nf_outeriter, 192
.set nb234nf_inneriter, 200
.set nb234nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb234nf_ixO, 0
.set nb234nf_iyO, 16
.set nb234nf_izO, 32
.set nb234nf_ixH1, 48
.set nb234nf_iyH1, 64
.set nb234nf_izH1, 80
.set nb234nf_ixH2, 96
.set nb234nf_iyH2, 112
.set nb234nf_izH2, 128
.set nb234nf_ixM, 144
.set nb234nf_iyM, 160
.set nb234nf_izM, 176
.set nb234nf_jxO, 192
.set nb234nf_jyO, 208
.set nb234nf_jzO, 224
.set nb234nf_jxH1, 240
.set nb234nf_jyH1, 256
.set nb234nf_jzH1, 272
.set nb234nf_jxH2, 288
.set nb234nf_jyH2, 304
.set nb234nf_jzH2, 320
.set nb234nf_jxM, 336
.set nb234nf_jyM, 352
.set nb234nf_jzM, 368
.set nb234nf_qqMM, 384
.set nb234nf_qqMH, 400
.set nb234nf_qqHH, 416
.set nb234nf_tsc, 432
.set nb234nf_c6, 448
.set nb234nf_c12, 464
.set nb234nf_vctot, 480
.set nb234nf_Vvdwtot, 496
.set nb234nf_half, 512
.set nb234nf_three, 528
.set nb234nf_rsqOO, 544
.set nb234nf_rsqH1H1, 560
.set nb234nf_rsqH1H2, 576
.set nb234nf_rsqH1M, 592
.set nb234nf_rsqH2H1, 608
.set nb234nf_rsqH2H2, 624
.set nb234nf_rsqH2M, 640
.set nb234nf_rsqMH1, 656
.set nb234nf_rsqMH2, 672
.set nb234nf_rsqMM, 688
.set nb234nf_rinvOO, 704
.set nb234nf_rinvH1H1, 720
.set nb234nf_rinvH1H2, 736
.set nb234nf_rinvH1M, 752
.set nb234nf_rinvH2H1, 768
.set nb234nf_rinvH2H2, 784
.set nb234nf_rinvH2M, 800
.set nb234nf_rinvMH1, 816
.set nb234nf_rinvMH2, 832
.set nb234nf_rinvMM, 848
.set nb234nf_krf, 864
.set nb234nf_crf, 880
.set nb234nf_is3, 896
.set nb234nf_ii3, 900
.set nb234nf_nri, 904
.set nb234nf_iinr, 912
.set nb234nf_jindex, 920
.set nb234nf_jjnr, 928
.set nb234nf_shift, 936
.set nb234nf_shiftvec, 944
.set nb234nf_facel, 952
.set nb234nf_innerjjnr, 960
.set nb234nf_innerk, 968
.set nb234nf_n, 972
.set nb234nf_nn1, 976
.set nb234nf_nouter, 980
.set nb234nf_ninner, 984
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1000,%rsp         ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb234nf_nouter(%rsp)
        movl %eax,nb234nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb234nf_nri(%rsp)
        movq %rsi,nb234nf_iinr(%rsp)
        movq %rdx,nb234nf_jindex(%rsp)
        movq %rcx,nb234nf_jjnr(%rsp)
        movq %r8,nb234nf_shift(%rsp)
        movq %r9,nb234nf_shiftvec(%rsp)
        movq nb234nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb234nf_facel(%rsp)

        movq nb234nf_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb234nf_tsc(%rsp)

        movq nb234nf_argkrf(%rbp),%rsi
        movq nb234nf_argcrf(%rbp),%rdi
        movsd (%rsi),%xmm1
        movsd (%rdi),%xmm2
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        movapd %xmm1,nb234nf_krf(%rsp)
        movapd %xmm2,nb234nf_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb234nf_half(%rsp)
        movl %ebx,nb234nf_half+4(%rsp)
        movsd nb234nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb234nf_half(%rsp)
        movapd %xmm3,nb234nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb234nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb234nf_charge(%rbp),%rdx
        movsd 24(%rdx,%rbx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%rdx,%rbx,8),%xmm5

        movsd nb234nf_facel(%rsp),%xmm6
        mulsd  %xmm3,%xmm3
        mulsd  %xmm5,%xmm4
        mulsd  %xmm5,%xmm5
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb234nf_qqMM(%rsp)
        movapd %xmm4,nb234nf_qqMH(%rsp)
        movapd %xmm5,nb234nf_qqHH(%rsp)

        xorpd %xmm0,%xmm0
        movq  nb234nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq nb234nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movq  nb234nf_vdwparam(%rbp),%rax
        movlpd (%rax,%rdx,8),%xmm0
        movlpd 8(%rax,%rdx,8),%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb234nf_c6(%rsp)
        movapd %xmm1,nb234nf_c12(%rsp)

_nb_kernel234nf_x86_64_sse2.nb234nf_threadloop: 
        movq  nb234nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel234nf_x86_64_sse2.nb234nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel234nf_x86_64_sse2.nb234nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb234nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb234nf_n(%rsp)
        movl %ebx,nb234nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel234nf_x86_64_sse2.nb234nf_outerstart
        jmp _nb_kernel234nf_x86_64_sse2.nb234nf_end

_nb_kernel234nf_x86_64_sse2.nb234nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb234nf_nouter(%rsp),%ebx
        movl %ebx,nb234nf_nouter(%rsp)

_nb_kernel234nf_x86_64_sse2.nb234nf_outer: 
        movq  nb234nf_shift(%rsp),%rax        ## eax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb234nf_is3(%rsp)            ## store is3 

        movq  nb234nf_shiftvec(%rsp),%rax     ## eax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb234nf_iinr(%rsp),%rcx         ## ecx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb234nf_pos(%rbp),%rax      ## eax = base of pos[]  
        movl  %ebx,nb234nf_ii3(%rsp)

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
        movapd %xmm3,nb234nf_ixO(%rsp)
        movapd %xmm4,nb234nf_iyO(%rsp)
        movapd %xmm5,nb234nf_izO(%rsp)
        movapd %xmm6,nb234nf_ixH1(%rsp)
        movapd %xmm7,nb234nf_iyH1(%rsp)

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
        movapd %xmm6,nb234nf_izH1(%rsp)
        movapd %xmm0,nb234nf_ixH2(%rsp)
        movapd %xmm1,nb234nf_iyH2(%rsp)
        movapd %xmm2,nb234nf_izH2(%rsp)
        movapd %xmm3,nb234nf_ixM(%rsp)
        movapd %xmm4,nb234nf_iyM(%rsp)
        movapd %xmm5,nb234nf_izM(%rsp)

        ## clear vctot
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb234nf_vctot(%rsp)
        movapd %xmm4,nb234nf_Vvdwtot(%rsp)

        movq  nb234nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb234nf_pos(%rbp),%rsi
        movq  nb234nf_faction(%rbp),%rdi
        movq  nb234nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb234nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb234nf_ninner(%rsp),%ecx
        movl  %ecx,nb234nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb234nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel234nf_x86_64_sse2.nb234nf_unroll_loop
        jmp   _nb_kernel234nf_x86_64_sse2.nb234nf_checksingle
_nb_kernel234nf_x86_64_sse2.nb234nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb234nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb234nf_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb234nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move j coordinates to local temp variables 
        ## load ox, oy, oz, h1x
        movlpd (%rsi,%rax,8),%xmm0
        movlpd (%rsi,%rbx,8),%xmm2
        movhpd 8(%rsi,%rax,8),%xmm0
        movhpd 8(%rsi,%rbx,8),%xmm2
        movlpd 16(%rsi,%rax,8),%xmm3
        movlpd 16(%rsi,%rbx,8),%xmm5
        movhpd 24(%rsi,%rax,8),%xmm3
        movhpd 24(%rsi,%rbx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## ox 
        unpckhpd %xmm2,%xmm1 ## oy
        unpcklpd %xmm5,%xmm3 ## ox 
        unpckhpd %xmm5,%xmm4 ## oy
        movapd  %xmm0,nb234nf_jxO(%rsp)
        movapd  %xmm1,nb234nf_jyO(%rsp)
        movapd  %xmm3,nb234nf_jzO(%rsp)
        movapd  %xmm4,nb234nf_jxH1(%rsp)

        ## load h1y, h1z, h2x, h2y 
        movlpd 32(%rsi,%rax,8),%xmm0
        movlpd 32(%rsi,%rbx,8),%xmm2
        movhpd 40(%rsi,%rax,8),%xmm0
        movhpd 40(%rsi,%rbx,8),%xmm2
        movlpd 48(%rsi,%rax,8),%xmm3
        movlpd 48(%rsi,%rbx,8),%xmm5
        movhpd 56(%rsi,%rax,8),%xmm3
        movhpd 56(%rsi,%rbx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## h1y
        unpckhpd %xmm2,%xmm1 ## h1z
        unpcklpd %xmm5,%xmm3 ## h2x
        unpckhpd %xmm5,%xmm4 ## h2y
        movapd  %xmm0,nb234nf_jyH1(%rsp)
        movapd  %xmm1,nb234nf_jzH1(%rsp)
        movapd  %xmm3,nb234nf_jxH2(%rsp)
        movapd  %xmm4,nb234nf_jyH2(%rsp)

        ## load h2z, mx, my, mz
        movlpd 64(%rsi,%rax,8),%xmm0
        movlpd 64(%rsi,%rbx,8),%xmm2
        movhpd 72(%rsi,%rax,8),%xmm0
        movhpd 72(%rsi,%rbx,8),%xmm2
        movlpd 80(%rsi,%rax,8),%xmm3
        movlpd 80(%rsi,%rbx,8),%xmm5
        movhpd 88(%rsi,%rax,8),%xmm3
        movhpd 88(%rsi,%rbx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## h2z
        unpckhpd %xmm2,%xmm1 ## mx
        unpcklpd %xmm5,%xmm3 ## my
        unpckhpd %xmm5,%xmm4 ## mz
        movapd  %xmm0,nb234nf_jzH2(%rsp)
        movapd  %xmm1,nb234nf_jxM(%rsp)
        movapd  %xmm3,nb234nf_jyM(%rsp)
        movapd  %xmm4,nb234nf_jzM(%rsp)

        ## start calculating pairwise distances
        movapd nb234nf_ixO(%rsp),%xmm0
        movapd nb234nf_iyO(%rsp),%xmm1
        movapd nb234nf_izO(%rsp),%xmm2
        movapd nb234nf_ixH1(%rsp),%xmm3
        movapd nb234nf_iyH1(%rsp),%xmm4
        movapd nb234nf_izH1(%rsp),%xmm5
        subpd  nb234nf_jxO(%rsp),%xmm0
        subpd  nb234nf_jyO(%rsp),%xmm1
        subpd  nb234nf_jzO(%rsp),%xmm2
        subpd  nb234nf_jxH1(%rsp),%xmm3
        subpd  nb234nf_jyH1(%rsp),%xmm4
        subpd  nb234nf_jzH1(%rsp),%xmm5
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb234nf_rsqOO(%rsp)
        movapd %xmm3,nb234nf_rsqH1H1(%rsp)

        movapd nb234nf_ixH1(%rsp),%xmm0
        movapd nb234nf_iyH1(%rsp),%xmm1
        movapd nb234nf_izH1(%rsp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb234nf_jxH2(%rsp),%xmm0
        subpd  nb234nf_jyH2(%rsp),%xmm1
        subpd  nb234nf_jzH2(%rsp),%xmm2
        subpd  nb234nf_jxM(%rsp),%xmm3
        subpd  nb234nf_jyM(%rsp),%xmm4
        subpd  nb234nf_jzM(%rsp),%xmm5
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb234nf_rsqH1H2(%rsp)
        movapd %xmm3,nb234nf_rsqH1M(%rsp)

        movapd nb234nf_ixH2(%rsp),%xmm0
        movapd nb234nf_iyH2(%rsp),%xmm1
        movapd nb234nf_izH2(%rsp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb234nf_jxH1(%rsp),%xmm0
        subpd  nb234nf_jyH1(%rsp),%xmm1
        subpd  nb234nf_jzH1(%rsp),%xmm2
        subpd  nb234nf_jxH2(%rsp),%xmm3
        subpd  nb234nf_jyH2(%rsp),%xmm4
        subpd  nb234nf_jzH2(%rsp),%xmm5
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb234nf_rsqH2H1(%rsp)
        movapd %xmm3,nb234nf_rsqH2H2(%rsp)

        movapd nb234nf_ixH2(%rsp),%xmm0
        movapd nb234nf_iyH2(%rsp),%xmm1
        movapd nb234nf_izH2(%rsp),%xmm2
        movapd nb234nf_ixM(%rsp),%xmm3
        movapd nb234nf_iyM(%rsp),%xmm4
        movapd nb234nf_izM(%rsp),%xmm5
        subpd  nb234nf_jxM(%rsp),%xmm0
        subpd  nb234nf_jyM(%rsp),%xmm1
        subpd  nb234nf_jzM(%rsp),%xmm2
        subpd  nb234nf_jxH1(%rsp),%xmm3
        subpd  nb234nf_jyH1(%rsp),%xmm4
        subpd  nb234nf_jzH1(%rsp),%xmm5
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb234nf_rsqH2M(%rsp)
        movapd %xmm4,nb234nf_rsqMH1(%rsp)

        movapd nb234nf_ixM(%rsp),%xmm0
        movapd nb234nf_iyM(%rsp),%xmm1
        movapd nb234nf_izM(%rsp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb234nf_jxH2(%rsp),%xmm0
        subpd  nb234nf_jyH2(%rsp),%xmm1
        subpd  nb234nf_jzH2(%rsp),%xmm2
        subpd  nb234nf_jxM(%rsp),%xmm3
        subpd  nb234nf_jyM(%rsp),%xmm4
        subpd  nb234nf_jzM(%rsp),%xmm5
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb234nf_rsqMH2(%rsp)
        movapd %xmm4,nb234nf_rsqMM(%rsp)

        ## Invsqrt form rsq M-H2 (rsq in xmm0) and MM (rsq in xmm4) 
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1  ## luA
        cvtps2pd %xmm5,%xmm5  ## luB

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb234nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb234nf_half(%rsp),%xmm3   ## iter1 
        mulpd   nb234nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb234nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb234nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb234nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb234nf_rinvMH2(%rsp)
        movapd %xmm5,nb234nf_rinvMM(%rsp)

        movapd nb234nf_rsqOO(%rsp),%xmm0
        movapd nb234nf_rsqH1H1(%rsp),%xmm4
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1
        cvtps2pd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb234nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb234nf_half(%rsp),%xmm3   ## iter1 of  
        mulpd   nb234nf_half(%rsp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb234nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb234nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb234nf_half(%rsp),%xmm5   ## rinv
        movapd %xmm1,nb234nf_rinvOO(%rsp)
        movapd %xmm5,nb234nf_rinvH1H1(%rsp)

        movapd nb234nf_rsqH1H2(%rsp),%xmm0
        movapd nb234nf_rsqH1M(%rsp),%xmm4
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1
        cvtps2pd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb234nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb234nf_half(%rsp),%xmm3   ## iter1 
        mulpd   nb234nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb234nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb234nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb234nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb234nf_rinvH1H2(%rsp)
        movapd %xmm5,nb234nf_rinvH1M(%rsp)

        movapd nb234nf_rsqH2H1(%rsp),%xmm0
        movapd nb234nf_rsqH2H2(%rsp),%xmm4
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1
        cvtps2pd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb234nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb234nf_half(%rsp),%xmm3   ## iter1a 
        mulpd   nb234nf_half(%rsp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb234nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb234nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb234nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb234nf_rinvH2H1(%rsp)
        movapd %xmm5,nb234nf_rinvH2H2(%rsp)

        movapd nb234nf_rsqMH1(%rsp),%xmm0
        movapd nb234nf_rsqH2M(%rsp),%xmm4
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1
        cvtps2pd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb234nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb234nf_half(%rsp),%xmm3   ## iter1a 
        mulpd   nb234nf_half(%rsp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb234nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb234nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb234nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb234nf_rinvMH1(%rsp)
        movapd %xmm5,nb234nf_rinvH2M(%rsp)

        ## start with OO interaction 
        movapd nb234nf_rinvOO(%rsp),%xmm0
        movapd nb234nf_rsqOO(%rsp),%xmm4

                mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb234nf_tsc(%rsp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movd %eax,%mm0
        movd %ebx,%mm1

        movq nb234nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx

        ## dispersion 
        movlpd (%rsi,%rax,8),%xmm4      ## Y1   
        movlpd (%rsi,%rbx,8),%xmm3      ## Y2 
        movhpd 8(%rsi,%rax,8),%xmm4     ## Y1 F1        
        movhpd 8(%rsi,%rbx,8),%xmm3     ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%rsi,%rax,8),%xmm6    ## G1
        movlpd 16(%rsi,%rbx,8),%xmm3    ## G2
        movhpd 24(%rsi,%rax,8),%xmm6    ## G1 H1        
        movhpd 24(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## dispersion table ready, in xmm4-xmm7         
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb234nf_c6(%rsp),%xmm4
        mulpd  %xmm4,%xmm5       ## Vvdw6 

        ## Update Vvdwtot directly 
        addpd  nb234nf_Vvdwtot(%rsp),%xmm5
        movapd %xmm5,nb234nf_Vvdwtot(%rsp)

        ## repulsion 
        movlpd 32(%rsi,%rax,8),%xmm4    ## Y1   
        movlpd 32(%rsi,%rbx,8),%xmm3    ## Y2 
        movhpd 40(%rsi,%rax,8),%xmm4    ## Y1 F1        
        movhpd 40(%rsi,%rbx,8),%xmm3    ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 48(%rsi,%rax,8),%xmm6    ## G1
        movlpd 48(%rsi,%rbx,8),%xmm3    ## G2
        movhpd 56(%rsi,%rax,8),%xmm6    ## G1 H1        
        movhpd 56(%rsi,%rbx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 

        ## table ready, in xmm4-xmm7    
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb234nf_c12(%rsp),%xmm4
        mulpd  %xmm4,%xmm5

        addpd  nb234nf_Vvdwtot(%rsp),%xmm5
        movapd %xmm5,nb234nf_Vvdwtot(%rsp)

        ## H1-H1 interaction 
        movapd nb234nf_rinvH1H1(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb234nf_krf(%rsp),%xmm5
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb234nf_rsqH1H1(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subpd  nb234nf_crf(%rsp),%xmm6
        mulpd  nb234nf_qqHH(%rsp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf)      
        addpd  nb234nf_vctot(%rsp),%xmm6   ## local vctot summation variable 

        ## H1-H2 interaction  
        movapd nb234nf_rinvH1H2(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb234nf_krf(%rsp),%xmm5
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb234nf_rsqH1H2(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm6=rinv+ krsq 
        subpd  nb234nf_crf(%rsp),%xmm4
        mulpd  nb234nf_qqHH(%rsp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf)      
        addpd  %xmm4,%xmm6

        ## H1-M interaction 
        movapd nb234nf_rinvH1M(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb234nf_krf(%rsp),%xmm5
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb234nf_rsqH1M(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm6=rinv+ krsq 
        subpd  nb234nf_crf(%rsp),%xmm4
        mulpd  nb234nf_qqMH(%rsp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm4,%xmm6

        ## H2-H1 interaction 
        movapd nb234nf_rinvH2H1(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb234nf_krf(%rsp),%xmm5
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb234nf_rsqH2H1(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm6=rinv+ krsq 
        subpd  nb234nf_crf(%rsp),%xmm4
        mulpd  nb234nf_qqHH(%rsp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm4,%xmm6

        ## H2-H2 interaction 
        movapd nb234nf_rinvH2H2(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb234nf_krf(%rsp),%xmm5
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb234nf_rsqH2H2(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm6=rinv+ krsq 
        subpd  nb234nf_crf(%rsp),%xmm4
        mulpd  nb234nf_qqHH(%rsp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm4,%xmm6

        ## H2-M interaction 
        movapd nb234nf_rinvH2M(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb234nf_krf(%rsp),%xmm5
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb234nf_rsqH2M(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm6=rinv+ krsq 
        subpd  nb234nf_crf(%rsp),%xmm4
        mulpd  nb234nf_qqMH(%rsp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm4,%xmm6

        ## M-H1 interaction 
        movapd nb234nf_rinvMH1(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb234nf_krf(%rsp),%xmm5
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb234nf_rsqMH1(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm6=rinv+ krsq 
        subpd  nb234nf_crf(%rsp),%xmm4
        mulpd  nb234nf_qqMH(%rsp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm4,%xmm6

        ## M-H2 interaction 
        movapd nb234nf_rinvMH2(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb234nf_krf(%rsp),%xmm5
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb234nf_rsqMH2(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm6=rinv+ krsq 
        subpd  nb234nf_crf(%rsp),%xmm4
        mulpd  nb234nf_qqMH(%rsp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm4,%xmm6

        ## M-M interaction 
        movapd nb234nf_rinvMM(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb234nf_krf(%rsp),%xmm5
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb234nf_rsqMM(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm6=rinv+ krsq 
        subpd  nb234nf_crf(%rsp),%xmm4
        mulpd  nb234nf_qqMM(%rsp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm4,%xmm6
        movapd %xmm6,nb234nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb234nf_innerk(%rsp)
        jl    _nb_kernel234nf_x86_64_sse2.nb234nf_checksingle
        jmp   _nb_kernel234nf_x86_64_sse2.nb234nf_unroll_loop
_nb_kernel234nf_x86_64_sse2.nb234nf_checksingle: 
        movl  nb234nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel234nf_x86_64_sse2.nb234nf_dosingle
        jmp   _nb_kernel234nf_x86_64_sse2.nb234nf_updateouterdata
_nb_kernel234nf_x86_64_sse2.nb234nf_dosingle: 
        movq  nb234nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb234nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move j coordinates to local temp variables 
        ## load ox, oy, oz, h1x
        movlpd (%rsi,%rax,8),%xmm0
        movhpd 8(%rsi,%rax,8),%xmm0
        movlpd 16(%rsi,%rax,8),%xmm1
        movhpd 24(%rsi,%rax,8),%xmm1
        movlpd 32(%rsi,%rax,8),%xmm2
        movhpd 40(%rsi,%rax,8),%xmm2
        movlpd 48(%rsi,%rax,8),%xmm3
        movhpd 56(%rsi,%rax,8),%xmm3
        movlpd 64(%rsi,%rax,8),%xmm4
        movhpd 72(%rsi,%rax,8),%xmm4
        movlpd 80(%rsi,%rax,8),%xmm5
        movhpd 88(%rsi,%rax,8),%xmm5
        movsd  %xmm0,nb234nf_jxO(%rsp)
        movsd  %xmm1,nb234nf_jzO(%rsp)
        movsd  %xmm2,nb234nf_jyH1(%rsp)
        movsd  %xmm3,nb234nf_jxH2(%rsp)
        movsd  %xmm4,nb234nf_jzH2(%rsp)
        movsd  %xmm5,nb234nf_jyM(%rsp)
        unpckhpd %xmm0,%xmm0
        unpckhpd %xmm1,%xmm1
        unpckhpd %xmm2,%xmm2
        unpckhpd %xmm3,%xmm3
        unpckhpd %xmm4,%xmm4
        unpckhpd %xmm5,%xmm5
        movsd  %xmm0,nb234nf_jyO(%rsp)
        movsd  %xmm1,nb234nf_jxH1(%rsp)
        movsd  %xmm2,nb234nf_jzH1(%rsp)
        movsd  %xmm3,nb234nf_jyH2(%rsp)
        movsd  %xmm4,nb234nf_jxM(%rsp)
        movsd  %xmm5,nb234nf_jzM(%rsp)

        ## start calculating pairwise distances
        movapd nb234nf_ixO(%rsp),%xmm0
        movapd nb234nf_iyO(%rsp),%xmm1
        movapd nb234nf_izO(%rsp),%xmm2
        movapd nb234nf_ixH1(%rsp),%xmm3
        movapd nb234nf_iyH1(%rsp),%xmm4
        movapd nb234nf_izH1(%rsp),%xmm5
        subsd  nb234nf_jxO(%rsp),%xmm0
        subsd  nb234nf_jyO(%rsp),%xmm1
        subsd  nb234nf_jzO(%rsp),%xmm2
        subsd  nb234nf_jxH1(%rsp),%xmm3
        subsd  nb234nf_jyH1(%rsp),%xmm4
        subsd  nb234nf_jzH1(%rsp),%xmm5
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb234nf_rsqOO(%rsp)
        movapd %xmm3,nb234nf_rsqH1H1(%rsp)

        movapd nb234nf_ixH1(%rsp),%xmm0
        movapd nb234nf_iyH1(%rsp),%xmm1
        movapd nb234nf_izH1(%rsp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb234nf_jxH2(%rsp),%xmm0
        subsd  nb234nf_jyH2(%rsp),%xmm1
        subsd  nb234nf_jzH2(%rsp),%xmm2
        subsd  nb234nf_jxM(%rsp),%xmm3
        subsd  nb234nf_jyM(%rsp),%xmm4
        subsd  nb234nf_jzM(%rsp),%xmm5
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb234nf_rsqH1H2(%rsp)
        movapd %xmm3,nb234nf_rsqH1M(%rsp)

        movapd nb234nf_ixH2(%rsp),%xmm0
        movapd nb234nf_iyH2(%rsp),%xmm1
        movapd nb234nf_izH2(%rsp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb234nf_jxH1(%rsp),%xmm0
        subsd  nb234nf_jyH1(%rsp),%xmm1
        subsd  nb234nf_jzH1(%rsp),%xmm2
        subsd  nb234nf_jxH2(%rsp),%xmm3
        subsd  nb234nf_jyH2(%rsp),%xmm4
        subsd  nb234nf_jzH2(%rsp),%xmm5
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb234nf_rsqH2H1(%rsp)
        movapd %xmm3,nb234nf_rsqH2H2(%rsp)

        movapd nb234nf_ixH2(%rsp),%xmm0
        movapd nb234nf_iyH2(%rsp),%xmm1
        movapd nb234nf_izH2(%rsp),%xmm2
        movapd nb234nf_ixM(%rsp),%xmm3
        movapd nb234nf_iyM(%rsp),%xmm4
        movapd nb234nf_izM(%rsp),%xmm5
        subsd  nb234nf_jxM(%rsp),%xmm0
        subsd  nb234nf_jyM(%rsp),%xmm1
        subsd  nb234nf_jzM(%rsp),%xmm2
        subsd  nb234nf_jxH1(%rsp),%xmm3
        subsd  nb234nf_jyH1(%rsp),%xmm4
        subsd  nb234nf_jzH1(%rsp),%xmm5
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb234nf_rsqH2M(%rsp)
        movapd %xmm4,nb234nf_rsqMH1(%rsp)

        movapd nb234nf_ixM(%rsp),%xmm0
        movapd nb234nf_iyM(%rsp),%xmm1
        movapd nb234nf_izM(%rsp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb234nf_jxH2(%rsp),%xmm0
        subsd  nb234nf_jyH2(%rsp),%xmm1
        subsd  nb234nf_jzH2(%rsp),%xmm2
        subsd  nb234nf_jxM(%rsp),%xmm3
        subsd  nb234nf_jyM(%rsp),%xmm4
        subsd  nb234nf_jzM(%rsp),%xmm5
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb234nf_rsqMH2(%rsp)
        movapd %xmm4,nb234nf_rsqMM(%rsp)

        ## Invsqrt form rsq M-H2 (rsq in xmm0) and MM (rsq in xmm4) 
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1  ## luA
        cvtss2sd %xmm5,%xmm5  ## luB

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb234nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb234nf_half(%rsp),%xmm3   ## iter1 
        mulsd   nb234nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb234nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb234nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb234nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb234nf_rinvMH2(%rsp)
        movapd %xmm5,nb234nf_rinvMM(%rsp)

        movapd nb234nf_rsqOO(%rsp),%xmm0
        movapd nb234nf_rsqH1H1(%rsp),%xmm4
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1
        cvtss2sd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb234nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb234nf_half(%rsp),%xmm3   ## iter1 of  
        mulsd   nb234nf_half(%rsp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb234nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb234nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb234nf_half(%rsp),%xmm5   ## rinv
        movapd %xmm1,nb234nf_rinvOO(%rsp)
        movapd %xmm5,nb234nf_rinvH1H1(%rsp)

        movapd nb234nf_rsqH1H2(%rsp),%xmm0
        movapd nb234nf_rsqH1M(%rsp),%xmm4
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1
        cvtss2sd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb234nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb234nf_half(%rsp),%xmm3   ## iter1 
        mulsd   nb234nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb234nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb234nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb234nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb234nf_rinvH1H2(%rsp)
        movapd %xmm5,nb234nf_rinvH1M(%rsp)

        movapd nb234nf_rsqH2H1(%rsp),%xmm0
        movapd nb234nf_rsqH2H2(%rsp),%xmm4
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1
        cvtss2sd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb234nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb234nf_half(%rsp),%xmm3   ## iter1a 
        mulsd   nb234nf_half(%rsp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb234nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb234nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb234nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb234nf_rinvH2H1(%rsp)
        movapd %xmm5,nb234nf_rinvH2H2(%rsp)

        movapd nb234nf_rsqMH1(%rsp),%xmm0
        movapd nb234nf_rsqH2M(%rsp),%xmm4
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1
        cvtss2sd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb234nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb234nf_half(%rsp),%xmm3   ## iter1a 
        mulsd   nb234nf_half(%rsp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb234nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb234nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb234nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb234nf_rinvMH1(%rsp)
        movapd %xmm5,nb234nf_rinvH2M(%rsp)

        ## start with OO interaction 
        movsd nb234nf_rinvOO(%rsp),%xmm0
        movsd nb234nf_rsqOO(%rsp),%xmm4

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb234nf_tsc(%rsp),%xmm4

        cvttsd2si %xmm4,%ebx    ## mm6 = lu idx 
        cvtsi2sd %ebx,%xmm5
        subsd %xmm5,%xmm4
        movsd %xmm4,%xmm1       ## xmm1=eps 
        movsd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $3,%ebx

        movq nb234nf_VFtab(%rbp),%rsi

        ## dispersion 
        movlpd (%rsi,%rbx,8),%xmm4      ## Y1   
        movhpd 8(%rsi,%rbx,8),%xmm4     ## Y1 F1        
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%rsi,%rbx,8),%xmm6    ## G1
        movhpd 24(%rsi,%rbx,8),%xmm6    ## G1 H1        
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## dispersion table ready, in xmm4-xmm7         
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb234nf_c6(%rsp),%xmm4
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        ## Update Vvdwtot directly 
        addsd  nb234nf_Vvdwtot(%rsp),%xmm5
        movsd %xmm5,nb234nf_Vvdwtot(%rsp)

        ## repulsion 
        movlpd 32(%rsi,%rbx,8),%xmm4    ## Y1   
        movhpd 40(%rsi,%rbx,8),%xmm4    ## Y1 F1        

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 48(%rsi,%rbx,8),%xmm6    ## G1
        movhpd 56(%rsi,%rbx,8),%xmm6    ## G1 H1        

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 

        ## table ready, in xmm4-xmm7    
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb234nf_c12(%rsp),%xmm4
        mulsd  %xmm4,%xmm5

        addsd  nb234nf_Vvdwtot(%rsp),%xmm5
        movsd %xmm5,nb234nf_Vvdwtot(%rsp)

        ## H1-H1 interaction 
        movsd nb234nf_rinvH1H1(%rsp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        movsd nb234nf_krf(%rsp),%xmm5
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb234nf_rsqH1H1(%rsp),%xmm5   ## xmm5=krsq 
        movsd %xmm5,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subsd  nb234nf_crf(%rsp),%xmm6
        mulsd  nb234nf_qqHH(%rsp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  nb234nf_vctot(%rsp),%xmm6   ## local vctot summation variable 

        ## H1-H2 interaction  
        movsd nb234nf_rinvH1H2(%rsp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        movsd nb234nf_krf(%rsp),%xmm5
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb234nf_rsqH1H2(%rsp),%xmm5   ## xmm5=krsq 
        movsd  %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm6=rinv+ krsq 
        subsd  nb234nf_crf(%rsp),%xmm4
        mulsd  nb234nf_qqHH(%rsp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm4,%xmm6

        ## H1-M interaction 
        movsd nb234nf_rinvH1M(%rsp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        movsd nb234nf_krf(%rsp),%xmm5
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb234nf_rsqH1M(%rsp),%xmm5   ## xmm5=krsq 
        movsd  %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm6=rinv+ krsq 
        subsd  nb234nf_crf(%rsp),%xmm4
        mulsd  nb234nf_qqMH(%rsp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm4,%xmm6

        ## H2-H1 interaction 
        movsd nb234nf_rinvH2H1(%rsp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        movsd nb234nf_krf(%rsp),%xmm5
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb234nf_rsqH2H1(%rsp),%xmm5   ## xmm5=krsq 
        movsd  %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm6=rinv+ krsq 
        subsd  nb234nf_crf(%rsp),%xmm4
        mulsd  nb234nf_qqHH(%rsp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm4,%xmm6

        ## H2-H2 interaction 
        movsd nb234nf_rinvH2H2(%rsp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        movsd nb234nf_krf(%rsp),%xmm5
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb234nf_rsqH2H2(%rsp),%xmm5   ## xmm5=krsq 
        movsd  %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm6=rinv+ krsq 
        subsd  nb234nf_crf(%rsp),%xmm4
        mulsd  nb234nf_qqHH(%rsp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm4,%xmm6

        ## H2-M interaction 
        movsd nb234nf_rinvH2M(%rsp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        movsd nb234nf_krf(%rsp),%xmm5
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb234nf_rsqH2M(%rsp),%xmm5   ## xmm5=krsq 
        movsd  %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm6=rinv+ krsq 
        subsd  nb234nf_crf(%rsp),%xmm4
        mulsd  nb234nf_qqMH(%rsp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm4,%xmm6

        ## M-H1 interaction 
        movsd nb234nf_rinvMH1(%rsp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        movsd nb234nf_krf(%rsp),%xmm5
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb234nf_rsqMH1(%rsp),%xmm5   ## xmm5=krsq 
        movsd  %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm6=rinv+ krsq 
        subsd  nb234nf_crf(%rsp),%xmm4
        mulsd  nb234nf_qqMH(%rsp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm4,%xmm6

        ## M-H2 interaction 
        movsd nb234nf_rinvMH2(%rsp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        movsd nb234nf_krf(%rsp),%xmm5
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb234nf_rsqMH2(%rsp),%xmm5   ## xmm5=krsq 
        movsd  %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm6=rinv+ krsq 
        subsd  nb234nf_crf(%rsp),%xmm4
        mulsd  nb234nf_qqMH(%rsp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm4,%xmm6

        ## M-M interaction 
        movsd nb234nf_rinvMM(%rsp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        movsd nb234nf_krf(%rsp),%xmm5
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb234nf_rsqMM(%rsp),%xmm5   ## xmm5=krsq 
        movsd  %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm6=rinv+ krsq 
        subsd  nb234nf_crf(%rsp),%xmm4
        mulsd  nb234nf_qqMM(%rsp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm4,%xmm6
        movsd %xmm6,nb234nf_vctot(%rsp)

_nb_kernel234nf_x86_64_sse2.nb234nf_updateouterdata: 
        ## get n from stack
        movl nb234nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb234nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb234nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb234nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb234nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb234nf_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb234nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel234nf_x86_64_sse2.nb234nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb234nf_n(%rsp)
        jmp _nb_kernel234nf_x86_64_sse2.nb234nf_outer
_nb_kernel234nf_x86_64_sse2.nb234nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb234nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel234nf_x86_64_sse2.nb234nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel234nf_x86_64_sse2.nb234nf_threadloop
_nb_kernel234nf_x86_64_sse2.nb234nf_end: 
        movl nb234nf_nouter(%rsp),%eax
        movl nb234nf_ninner(%rsp),%ebx
        movq nb234nf_outeriter(%rbp),%rcx
        movq nb234nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1000,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret


