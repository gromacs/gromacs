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


.globl nb_kernel134_ia32_sse2
.globl _nb_kernel134_ia32_sse2
nb_kernel134_ia32_sse2: 
_nb_kernel134_ia32_sse2:        
.set nb134_p_nri, 8
.set nb134_iinr, 12
.set nb134_jindex, 16
.set nb134_jjnr, 20
.set nb134_shift, 24
.set nb134_shiftvec, 28
.set nb134_fshift, 32
.set nb134_gid, 36
.set nb134_pos, 40
.set nb134_faction, 44
.set nb134_charge, 48
.set nb134_p_facel, 52
.set nb134_argkrf, 56
.set nb134_argcrf, 60
.set nb134_Vc, 64
.set nb134_type, 68
.set nb134_p_ntype, 72
.set nb134_vdwparam, 76
.set nb134_Vvdw, 80
.set nb134_p_tabscale, 84
.set nb134_VFtab, 88
.set nb134_invsqrta, 92
.set nb134_dvda, 96
.set nb134_p_gbtabscale, 100
.set nb134_GBtab, 104
.set nb134_p_nthreads, 108
.set nb134_count, 112
.set nb134_mtx, 116
.set nb134_outeriter, 120
.set nb134_inneriter, 124
.set nb134_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb134_ixO, 0
.set nb134_iyO, 16
.set nb134_izO, 32
.set nb134_ixH1, 48
.set nb134_iyH1, 64
.set nb134_izH1, 80
.set nb134_ixH2, 96
.set nb134_iyH2, 112
.set nb134_izH2, 128
.set nb134_ixM, 144
.set nb134_iyM, 160
.set nb134_izM, 176
.set nb134_jxO, 192
.set nb134_jyO, 208
.set nb134_jzO, 224
.set nb134_jxH1, 240
.set nb134_jyH1, 256
.set nb134_jzH1, 272
.set nb134_jxH2, 288
.set nb134_jyH2, 304
.set nb134_jzH2, 320
.set nb134_jxM, 336
.set nb134_jyM, 352
.set nb134_jzM, 368
.set nb134_dxOO, 384
.set nb134_dyOO, 400
.set nb134_dzOO, 416
.set nb134_dxH1H1, 432
.set nb134_dyH1H1, 448
.set nb134_dzH1H1, 464
.set nb134_dxH1H2, 480
.set nb134_dyH1H2, 496
.set nb134_dzH1H2, 512
.set nb134_dxH1M, 528
.set nb134_dyH1M, 544
.set nb134_dzH1M, 560
.set nb134_dxH2H1, 576
.set nb134_dyH2H1, 592
.set nb134_dzH2H1, 608
.set nb134_dxH2H2, 624
.set nb134_dyH2H2, 640
.set nb134_dzH2H2, 656
.set nb134_dxH2M, 672
.set nb134_dyH2M, 688
.set nb134_dzH2M, 704
.set nb134_dxMH1, 720
.set nb134_dyMH1, 736
.set nb134_dzMH1, 752
.set nb134_dxMH2, 768
.set nb134_dyMH2, 784
.set nb134_dzMH2, 800
.set nb134_dxMM, 816
.set nb134_dyMM, 832
.set nb134_dzMM, 848
.set nb134_qqMM, 864
.set nb134_qqMH, 880
.set nb134_qqHH, 896
.set nb134_two, 912
.set nb134_c6, 944
.set nb134_c12, 960
.set nb134_vctot, 976
.set nb134_Vvdwtot, 992
.set nb134_fixO, 1008
.set nb134_fiyO, 1024
.set nb134_fizO, 1040
.set nb134_fixH1, 1056
.set nb134_fiyH1, 1072
.set nb134_fizH1, 1088
.set nb134_fixH2, 1104
.set nb134_fiyH2, 1120
.set nb134_fizH2, 1136
.set nb134_fixM, 1152
.set nb134_fiyM, 1168
.set nb134_fizM, 1184
.set nb134_fjxO, 1200
.set nb134_fjyO, 1216
.set nb134_fjzO, 1232
.set nb134_fjxH1, 1248
.set nb134_fjyH1, 1264
.set nb134_fjzH1, 1280
.set nb134_fjxH2, 1296
.set nb134_fjyH2, 1312
.set nb134_fjzH2, 1328
.set nb134_fjxM, 1344
.set nb134_fjyM, 1360
.set nb134_fjzM, 1376
.set nb134_half, 1392
.set nb134_three, 1408
.set nb134_tsc, 1424
.set nb134_fstmp, 1440
.set nb134_rsqOO, 1456
.set nb134_rsqH1H1, 1472
.set nb134_rsqH1H2, 1488
.set nb134_rsqH1M, 1504
.set nb134_rsqH2H1, 1520
.set nb134_rsqH2H2, 1536
.set nb134_rsqH2M, 1552
.set nb134_rsqMH1, 1568
.set nb134_rsqMH2, 1584
.set nb134_rsqMM, 1600
.set nb134_rinvOO, 1616
.set nb134_rinvH1H1, 1632
.set nb134_rinvH1H2, 1648
.set nb134_rinvH1M, 1664
.set nb134_rinvH2H1, 1680
.set nb134_rinvH2H2, 1696
.set nb134_rinvH2M, 1712
.set nb134_rinvMH1, 1728
.set nb134_rinvMH2, 1744
.set nb134_rinvMM, 1760
.set nb134_is3, 1808
.set nb134_ii3, 1812
.set nb134_innerjjnr, 1816
.set nb134_innerk, 1820
.set nb134_n, 1824
.set nb134_nn1, 1828
.set nb134_nri, 1832
.set nb134_nouter, 1836
.set nb134_ninner, 1840
.set nb134_salign, 1844
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $1848,%esp         ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb134_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb134_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb134_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb134_nouter(%esp)
        movl %eax,nb134_ninner(%esp)

        movl nb134_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb134_tsc(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb134_half(%esp)
        movl %ebx,nb134_half+4(%esp)
        movsd nb134_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb134_half(%esp)
        movapd %xmm2,nb134_two(%esp)
        movapd %xmm3,nb134_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb134_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb134_charge(%ebp),%edx
        movsd 24(%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb134_p_facel(%ebp),%esi
        movsd (%esi),%xmm6
        mulsd  %xmm3,%xmm3
        mulsd  %xmm5,%xmm4
        mulsd  %xmm5,%xmm5
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb134_qqMM(%esp)
        movapd %xmm4,nb134_qqMH(%esp)
        movapd %xmm5,nb134_qqHH(%esp)

        xorpd %xmm0,%xmm0
        movl  nb134_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb134_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb134_vdwparam(%ebp),%eax
        movlpd (%eax,%edx,8),%xmm0
        movhpd 8(%eax,%edx,8),%xmm0
        movhlps %xmm0,%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb134_c6(%esp)
        movapd %xmm1,nb134_c12(%esp)

_nb_kernel134_ia32_sse2.nb134_threadloop: 
        movl  nb134_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel134_ia32_sse2.nb134_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel134_ia32_sse2.nb134_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb134_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb134_n(%esp)
        movl %ebx,nb134_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel134_ia32_sse2.nb134_outerstart
        jmp _nb_kernel134_ia32_sse2.nb134_end

_nb_kernel134_ia32_sse2.nb134_outerstart: 
        ## ebx contains number of outer iterations
        addl nb134_nouter(%esp),%ebx
        movl %ebx,nb134_nouter(%esp)

_nb_kernel134_ia32_sse2.nb134_outer: 
        movl  nb134_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb134_is3(%esp)      ## store is3 

        movl  nb134_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb134_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb134_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb134_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3       ## ox
        addsd 8(%eax,%ebx,8),%xmm4      ## oy
        addsd 16(%eax,%ebx,8),%xmm5     ## oz   
        addsd 24(%eax,%ebx,8),%xmm6     ## h1x
        addsd 32(%eax,%ebx,8),%xmm7     ## h1y
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm7,%xmm7
        movapd %xmm3,nb134_ixO(%esp)
        movapd %xmm4,nb134_iyO(%esp)
        movapd %xmm5,nb134_izO(%esp)
        movapd %xmm6,nb134_ixH1(%esp)
        movapd %xmm7,nb134_iyH1(%esp)

        movsd %xmm2,%xmm6
        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 40(%eax,%ebx,8),%xmm6    ## h1z
        addsd 48(%eax,%ebx,8),%xmm0    ## h2x
        addsd 56(%eax,%ebx,8),%xmm1    ## h2y
        addsd 64(%eax,%ebx,8),%xmm2    ## h2z
        addsd 72(%eax,%ebx,8),%xmm3    ## mx
        addsd 80(%eax,%ebx,8),%xmm4    ## my
        addsd 88(%eax,%ebx,8),%xmm5    ## mz

        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm6,nb134_izH1(%esp)
        movapd %xmm0,nb134_ixH2(%esp)
        movapd %xmm1,nb134_iyH2(%esp)
        movapd %xmm2,nb134_izH2(%esp)
        movapd %xmm3,nb134_ixM(%esp)
        movapd %xmm4,nb134_iyM(%esp)
        movapd %xmm5,nb134_izM(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb134_vctot(%esp)
        movapd %xmm4,nb134_Vvdwtot(%esp)
        movapd %xmm4,nb134_fixO(%esp)
        movapd %xmm4,nb134_fiyO(%esp)
        movapd %xmm4,nb134_fizO(%esp)
        movapd %xmm4,nb134_fixH1(%esp)
        movapd %xmm4,nb134_fiyH1(%esp)
        movapd %xmm4,nb134_fizH1(%esp)
        movapd %xmm4,nb134_fixH2(%esp)
        movapd %xmm4,nb134_fiyH2(%esp)
        movapd %xmm4,nb134_fizH2(%esp)
        movapd %xmm4,nb134_fixM(%esp)
        movapd %xmm4,nb134_fiyM(%esp)
        movapd %xmm4,nb134_fizM(%esp)

        movl  nb134_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb134_pos(%ebp),%esi
        movl  nb134_faction(%ebp),%edi
        movl  nb134_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb134_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb134_ninner(%esp),%ecx
        movl  %ecx,nb134_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb134_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel134_ia32_sse2.nb134_unroll_loop
        jmp   _nb_kernel134_ia32_sse2.nb134_checksingle
_nb_kernel134_ia32_sse2.nb134_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb134_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb134_innerjjnr(%esp)            ## advance pointer (unrolled 2) 

        movl nb134_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move j coordinates to local temp variables 
        ## load ox, oy, oz, h1x
        movlpd (%esi,%eax,8),%xmm0
        movlpd (%esi,%ebx,8),%xmm2
        movhpd 8(%esi,%eax,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm2
        movlpd 16(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%ebx,8),%xmm5
        movhpd 24(%esi,%eax,8),%xmm3
        movhpd 24(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## ox 
        unpckhpd %xmm2,%xmm1 ## oy
        unpcklpd %xmm5,%xmm3 ## ox 
        unpckhpd %xmm5,%xmm4 ## oy
        movapd  %xmm0,nb134_jxO(%esp)
        movapd  %xmm1,nb134_jyO(%esp)
        movapd  %xmm3,nb134_jzO(%esp)
        movapd  %xmm4,nb134_jxH1(%esp)

        ## load h1y, h1z, h2x, h2y 
        movlpd 32(%esi,%eax,8),%xmm0
        movlpd 32(%esi,%ebx,8),%xmm2
        movhpd 40(%esi,%eax,8),%xmm0
        movhpd 40(%esi,%ebx,8),%xmm2
        movlpd 48(%esi,%eax,8),%xmm3
        movlpd 48(%esi,%ebx,8),%xmm5
        movhpd 56(%esi,%eax,8),%xmm3
        movhpd 56(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## h1y
        unpckhpd %xmm2,%xmm1 ## h1z
        unpcklpd %xmm5,%xmm3 ## h2x
        unpckhpd %xmm5,%xmm4 ## h2y
        movapd  %xmm0,nb134_jyH1(%esp)
        movapd  %xmm1,nb134_jzH1(%esp)
        movapd  %xmm3,nb134_jxH2(%esp)
        movapd  %xmm4,nb134_jyH2(%esp)

        ## load h2z, mx, my, mz
        movlpd 64(%esi,%eax,8),%xmm0
        movlpd 64(%esi,%ebx,8),%xmm2
        movhpd 72(%esi,%eax,8),%xmm0
        movhpd 72(%esi,%ebx,8),%xmm2
        movlpd 80(%esi,%eax,8),%xmm3
        movlpd 80(%esi,%ebx,8),%xmm5
        movhpd 88(%esi,%eax,8),%xmm3
        movhpd 88(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## h2z
        unpckhpd %xmm2,%xmm1 ## mx
        unpcklpd %xmm5,%xmm3 ## my
        unpckhpd %xmm5,%xmm4 ## mz
        movapd  %xmm0,nb134_jzH2(%esp)
        movapd  %xmm1,nb134_jxM(%esp)
        movapd  %xmm3,nb134_jyM(%esp)
        movapd  %xmm4,nb134_jzM(%esp)

        ## start calculating pairwise distances
        movapd nb134_ixO(%esp),%xmm0
        movapd nb134_iyO(%esp),%xmm1
        movapd nb134_izO(%esp),%xmm2
        movapd nb134_ixH1(%esp),%xmm3
        movapd nb134_iyH1(%esp),%xmm4
        movapd nb134_izH1(%esp),%xmm5
        subpd  nb134_jxO(%esp),%xmm0
        subpd  nb134_jyO(%esp),%xmm1
        subpd  nb134_jzO(%esp),%xmm2
        subpd  nb134_jxH1(%esp),%xmm3
        subpd  nb134_jyH1(%esp),%xmm4
        subpd  nb134_jzH1(%esp),%xmm5
        movapd %xmm0,nb134_dxOO(%esp)
        movapd %xmm1,nb134_dyOO(%esp)
        movapd %xmm2,nb134_dzOO(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb134_dxH1H1(%esp)
        movapd %xmm4,nb134_dyH1H1(%esp)
        movapd %xmm5,nb134_dzH1H1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb134_rsqOO(%esp)
        movapd %xmm3,nb134_rsqH1H1(%esp)

        movapd nb134_ixH1(%esp),%xmm0
        movapd nb134_iyH1(%esp),%xmm1
        movapd nb134_izH1(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb134_jxH2(%esp),%xmm0
        subpd  nb134_jyH2(%esp),%xmm1
        subpd  nb134_jzH2(%esp),%xmm2
        subpd  nb134_jxM(%esp),%xmm3
        subpd  nb134_jyM(%esp),%xmm4
        subpd  nb134_jzM(%esp),%xmm5
        movapd %xmm0,nb134_dxH1H2(%esp)
        movapd %xmm1,nb134_dyH1H2(%esp)
        movapd %xmm2,nb134_dzH1H2(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb134_dxH1M(%esp)
        movapd %xmm4,nb134_dyH1M(%esp)
        movapd %xmm5,nb134_dzH1M(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb134_rsqH1H2(%esp)
        movapd %xmm3,nb134_rsqH1M(%esp)

        movapd nb134_ixH2(%esp),%xmm0
        movapd nb134_iyH2(%esp),%xmm1
        movapd nb134_izH2(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb134_jxH1(%esp),%xmm0
        subpd  nb134_jyH1(%esp),%xmm1
        subpd  nb134_jzH1(%esp),%xmm2
        subpd  nb134_jxH2(%esp),%xmm3
        subpd  nb134_jyH2(%esp),%xmm4
        subpd  nb134_jzH2(%esp),%xmm5
        movapd %xmm0,nb134_dxH2H1(%esp)
        movapd %xmm1,nb134_dyH2H1(%esp)
        movapd %xmm2,nb134_dzH2H1(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb134_dxH2H2(%esp)
        movapd %xmm4,nb134_dyH2H2(%esp)
        movapd %xmm5,nb134_dzH2H2(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb134_rsqH2H1(%esp)
        movapd %xmm3,nb134_rsqH2H2(%esp)

        movapd nb134_ixH2(%esp),%xmm0
        movapd nb134_iyH2(%esp),%xmm1
        movapd nb134_izH2(%esp),%xmm2
        movapd nb134_ixM(%esp),%xmm3
        movapd nb134_iyM(%esp),%xmm4
        movapd nb134_izM(%esp),%xmm5
        subpd  nb134_jxM(%esp),%xmm0
        subpd  nb134_jyM(%esp),%xmm1
        subpd  nb134_jzM(%esp),%xmm2
        subpd  nb134_jxH1(%esp),%xmm3
        subpd  nb134_jyH1(%esp),%xmm4
        subpd  nb134_jzH1(%esp),%xmm5
        movapd %xmm0,nb134_dxH2M(%esp)
        movapd %xmm1,nb134_dyH2M(%esp)
        movapd %xmm2,nb134_dzH2M(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb134_dxMH1(%esp)
        movapd %xmm4,nb134_dyMH1(%esp)
        movapd %xmm5,nb134_dzMH1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb134_rsqH2M(%esp)
        movapd %xmm4,nb134_rsqMH1(%esp)

        movapd nb134_ixM(%esp),%xmm0
        movapd nb134_iyM(%esp),%xmm1
        movapd nb134_izM(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb134_jxH2(%esp),%xmm0
        subpd  nb134_jyH2(%esp),%xmm1
        subpd  nb134_jzH2(%esp),%xmm2
        subpd  nb134_jxM(%esp),%xmm3
        subpd  nb134_jyM(%esp),%xmm4
        subpd  nb134_jzM(%esp),%xmm5
        movapd %xmm0,nb134_dxMH2(%esp)
        movapd %xmm1,nb134_dyMH2(%esp)
        movapd %xmm2,nb134_dzMH2(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb134_dxMM(%esp)
        movapd %xmm4,nb134_dyMM(%esp)
        movapd %xmm5,nb134_dzMM(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb134_rsqMH2(%esp)
        movapd %xmm4,nb134_rsqMM(%esp)

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
        movapd  nb134_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb134_half(%esp),%xmm3   ## iter1 
        mulpd   nb134_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb134_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb134_half(%esp),%xmm1   ## rinv 
        mulpd   nb134_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb134_rinvMH2(%esp)
        movapd %xmm5,nb134_rinvMM(%esp)

        movapd nb134_rsqOO(%esp),%xmm0
        movapd nb134_rsqH1H1(%esp),%xmm4
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
        movapd  nb134_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb134_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb134_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb134_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb134_half(%esp),%xmm1   ## rinv 
        mulpd   nb134_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb134_rinvOO(%esp)
        movapd %xmm5,nb134_rinvH1H1(%esp)

        movapd nb134_rsqH1H2(%esp),%xmm0
        movapd nb134_rsqH1M(%esp),%xmm4
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
        movapd  nb134_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb134_half(%esp),%xmm3   ## iter1 
        mulpd   nb134_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb134_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb134_half(%esp),%xmm1   ## rinv 
        mulpd   nb134_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb134_rinvH1H2(%esp)
        movapd %xmm5,nb134_rinvH1M(%esp)

        movapd nb134_rsqH2H1(%esp),%xmm0
        movapd nb134_rsqH2H2(%esp),%xmm4
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
        movapd  nb134_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb134_half(%esp),%xmm3   ## iter1a 
        mulpd   nb134_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb134_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb134_half(%esp),%xmm1   ## rinv 
        mulpd   nb134_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb134_rinvH2H1(%esp)
        movapd %xmm5,nb134_rinvH2H2(%esp)

        movapd nb134_rsqMH1(%esp),%xmm0
        movapd nb134_rsqH2M(%esp),%xmm4
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
        movapd  nb134_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb134_half(%esp),%xmm3   ## iter1a 
        mulpd   nb134_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb134_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb134_half(%esp),%xmm1   ## rinv 
        mulpd   nb134_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb134_rinvMH1(%esp)
        movapd %xmm5,nb134_rinvH2M(%esp)

        ## start with OO interaction 
        movapd nb134_rinvOO(%esp),%xmm0
        movapd nb134_rsqOO(%esp),%xmm4

                mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb134_tsc(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movd %eax,%mm0
        movd %ebx,%mm1

        movl nb134_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx

        ## dispersion 
        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2 
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## dispersion table ready, in xmm4-xmm7         
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb134_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb134_c6(%esp),%xmm4
        mulpd  %xmm4,%xmm7       ## fijD 
        mulpd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addpd  nb134_Vvdwtot(%esp),%xmm5
        xorpd  %xmm3,%xmm3
        mulpd  nb134_tsc(%esp),%xmm7
        subpd  %xmm7,%xmm3
        movapd %xmm3,nb134_fstmp(%esp)
        movapd %xmm5,nb134_Vvdwtot(%esp)

        ## repulsion 
        movlpd 32(%esi,%eax,8),%xmm4    ## Y1   
        movlpd 32(%esi,%ebx,8),%xmm3    ## Y2 
        movhpd 40(%esi,%eax,8),%xmm4    ## Y1 F1        
        movhpd 40(%esi,%ebx,8),%xmm3    ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 48(%esi,%eax,8),%xmm6    ## G1
        movlpd 48(%esi,%ebx,8),%xmm3    ## G2
        movhpd 56(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 56(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 

        ## table ready, in xmm4-xmm7    
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb134_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb134_c12(%esp),%xmm4
        mulpd  %xmm4,%xmm7
        mulpd  %xmm4,%xmm5

        addpd  nb134_Vvdwtot(%esp),%xmm5
        movapd nb134_fstmp(%esp),%xmm3
        mulpd  nb134_tsc(%esp),%xmm7
        subpd  %xmm7,%xmm3
        movapd %xmm5,nb134_Vvdwtot(%esp)

        mulpd  %xmm3,%xmm0
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb134_dxOO(%esp),%xmm0
        mulpd nb134_dyOO(%esp),%xmm1
        mulpd nb134_dzOO(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb134_fixO(%esp),%xmm0
        addpd nb134_fiyO(%esp),%xmm1
        addpd nb134_fizO(%esp),%xmm2
        movapd %xmm3,nb134_fjxO(%esp)
        movapd %xmm4,nb134_fjyO(%esp)
        movapd %xmm5,nb134_fjzO(%esp)
        movapd %xmm0,nb134_fixO(%esp)
        movapd %xmm1,nb134_fiyO(%esp)
        movapd %xmm2,nb134_fizO(%esp)

        ## H1-H1 interaction 
        movapd nb134_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm6      ## xmm7=rinv 
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb134_qqHH(%esp),%xmm6   ## vcoul
        mulpd %xmm6,%xmm0

        addpd  nb134_vctot(%esp),%xmm6   ## local vctot summation variable 

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb134_dxH1H1(%esp),%xmm0
        mulpd nb134_dyH1H1(%esp),%xmm1
        mulpd nb134_dzH1H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb134_fixH1(%esp),%xmm0
        addpd nb134_fiyH1(%esp),%xmm1
        addpd nb134_fizH1(%esp),%xmm2
        movapd %xmm3,nb134_fjxH1(%esp)
        movapd %xmm4,nb134_fjyH1(%esp)
        movapd %xmm5,nb134_fjzH1(%esp)
        movapd %xmm0,nb134_fixH1(%esp)
        movapd %xmm1,nb134_fiyH1(%esp)
        movapd %xmm2,nb134_fizH1(%esp)

        ## H1-H2 interaction  
        movapd nb134_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb134_qqHH(%esp),%xmm7   ## vcoul
        mulpd %xmm7,%xmm0
        addpd  %xmm7,%xmm6

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb134_dxH1H2(%esp),%xmm0
        mulpd nb134_dyH1H2(%esp),%xmm1
        mulpd nb134_dzH1H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb134_fixH1(%esp),%xmm0
        addpd nb134_fiyH1(%esp),%xmm1
        addpd nb134_fizH1(%esp),%xmm2
        movapd %xmm3,nb134_fjxH2(%esp)
        movapd %xmm4,nb134_fjyH2(%esp)
        movapd %xmm5,nb134_fjzH2(%esp)
        movapd %xmm0,nb134_fixH1(%esp)
        movapd %xmm1,nb134_fiyH1(%esp)
        movapd %xmm2,nb134_fizH1(%esp)

        ## H1-M interaction 
        movapd nb134_rinvH1M(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb134_qqMH(%esp),%xmm7   ## vcoul
        mulpd %xmm7,%xmm0
        addpd  %xmm7,%xmm6

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb134_dxH1M(%esp),%xmm0
        mulpd nb134_dyH1M(%esp),%xmm1
        mulpd nb134_dzH1M(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb134_fixH1(%esp),%xmm0
        addpd nb134_fiyH1(%esp),%xmm1
        addpd nb134_fizH1(%esp),%xmm2
        movapd %xmm3,nb134_fjxM(%esp)
        movapd %xmm4,nb134_fjyM(%esp)
        movapd %xmm5,nb134_fjzM(%esp)
        movapd %xmm0,nb134_fixH1(%esp)
        movapd %xmm1,nb134_fiyH1(%esp)
        movapd %xmm2,nb134_fizH1(%esp)

        ## H2-H1 interaction 
        movapd nb134_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb134_qqHH(%esp),%xmm7   ## vcoul
        mulpd %xmm7,%xmm0
        addpd  %xmm7,%xmm6

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb134_fjxH1(%esp),%xmm3
        movapd nb134_fjyH1(%esp),%xmm4
        movapd nb134_fjzH1(%esp),%xmm5
        mulpd nb134_dxH2H1(%esp),%xmm0
        mulpd nb134_dyH2H1(%esp),%xmm1
        mulpd nb134_dzH2H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb134_fixH2(%esp),%xmm0
        addpd nb134_fiyH2(%esp),%xmm1
        addpd nb134_fizH2(%esp),%xmm2
        movapd %xmm3,nb134_fjxH1(%esp)
        movapd %xmm4,nb134_fjyH1(%esp)
        movapd %xmm5,nb134_fjzH1(%esp)
        movapd %xmm0,nb134_fixH2(%esp)
        movapd %xmm1,nb134_fiyH2(%esp)
        movapd %xmm2,nb134_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb134_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb134_qqHH(%esp),%xmm7   ## vcoul
        mulpd %xmm7,%xmm0
        addpd  %xmm7,%xmm6

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb134_fjxH2(%esp),%xmm3
        movapd nb134_fjyH2(%esp),%xmm4
        movapd nb134_fjzH2(%esp),%xmm5
        mulpd nb134_dxH2H2(%esp),%xmm0
        mulpd nb134_dyH2H2(%esp),%xmm1
        mulpd nb134_dzH2H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb134_fixH2(%esp),%xmm0
        addpd nb134_fiyH2(%esp),%xmm1
        addpd nb134_fizH2(%esp),%xmm2
        movapd %xmm3,nb134_fjxH2(%esp)
        movapd %xmm4,nb134_fjyH2(%esp)
        movapd %xmm5,nb134_fjzH2(%esp)
        movapd %xmm0,nb134_fixH2(%esp)
        movapd %xmm1,nb134_fiyH2(%esp)
        movapd %xmm2,nb134_fizH2(%esp)

        ## H2-M interaction 
        movapd nb134_rinvH2M(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb134_qqMH(%esp),%xmm7   ## vcoul
        mulpd %xmm7,%xmm0
        addpd  %xmm7,%xmm6

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb134_fjxM(%esp),%xmm3
        movapd nb134_fjyM(%esp),%xmm4
        movapd nb134_fjzM(%esp),%xmm5
        mulpd nb134_dxH2M(%esp),%xmm0
        mulpd nb134_dyH2M(%esp),%xmm1
        mulpd nb134_dzH2M(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb134_fixH2(%esp),%xmm0
        addpd nb134_fiyH2(%esp),%xmm1
        addpd nb134_fizH2(%esp),%xmm2
        movapd %xmm3,nb134_fjxM(%esp)
        movapd %xmm4,nb134_fjyM(%esp)
        movapd %xmm5,nb134_fjzM(%esp)
        movapd %xmm0,nb134_fixH2(%esp)
        movapd %xmm1,nb134_fiyH2(%esp)
        movapd %xmm2,nb134_fizH2(%esp)

        ## M-H1 interaction 
        movapd nb134_rinvMH1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb134_qqMH(%esp),%xmm7   ## vcoul
        mulpd %xmm7,%xmm0
        addpd  %xmm7,%xmm6

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb134_fjxH1(%esp),%xmm3
        movapd nb134_fjyH1(%esp),%xmm4
        movapd nb134_fjzH1(%esp),%xmm5
        mulpd nb134_dxMH1(%esp),%xmm0
        mulpd nb134_dyMH1(%esp),%xmm1
        mulpd nb134_dzMH1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb134_fixM(%esp),%xmm0
        addpd nb134_fiyM(%esp),%xmm1
        addpd nb134_fizM(%esp),%xmm2
        movapd %xmm3,nb134_fjxH1(%esp)
        movapd %xmm4,nb134_fjyH1(%esp)
        movapd %xmm5,nb134_fjzH1(%esp)
        movapd %xmm0,nb134_fixM(%esp)
        movapd %xmm1,nb134_fiyM(%esp)
        movapd %xmm2,nb134_fizM(%esp)

        ## M-H2 interaction 
        movapd nb134_rinvMH2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb134_qqMH(%esp),%xmm7   ## vcoul
        mulpd %xmm7,%xmm0
        addpd  %xmm7,%xmm6

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb134_fjxH2(%esp),%xmm3
        movapd nb134_fjyH2(%esp),%xmm4
        movapd nb134_fjzH2(%esp),%xmm5
        mulpd nb134_dxMH2(%esp),%xmm0
        mulpd nb134_dyMH2(%esp),%xmm1
        mulpd nb134_dzMH2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb134_fixM(%esp),%xmm0
        addpd nb134_fiyM(%esp),%xmm1
        addpd nb134_fizM(%esp),%xmm2
        movapd %xmm3,nb134_fjxH2(%esp)
        movapd %xmm4,nb134_fjyH2(%esp)
        movapd %xmm5,nb134_fjzH2(%esp)
        movapd %xmm0,nb134_fixM(%esp)
        movapd %xmm1,nb134_fiyM(%esp)
        movapd %xmm2,nb134_fizM(%esp)

        ## M-M interaction 
        movapd nb134_rinvMM(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb134_qqMM(%esp),%xmm7   ## vcoul
        mulpd %xmm7,%xmm0
        addpd  %xmm7,%xmm6
        movapd %xmm6,nb134_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb134_fjxM(%esp),%xmm3
        movapd nb134_fjyM(%esp),%xmm4
        movapd nb134_fjzM(%esp),%xmm5
        mulpd nb134_dxMM(%esp),%xmm0
        mulpd nb134_dyMM(%esp),%xmm1
        mulpd nb134_dzMM(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb134_fixM(%esp),%xmm0
        addpd nb134_fiyM(%esp),%xmm1
        addpd nb134_fizM(%esp),%xmm2
        movapd %xmm3,nb134_fjxM(%esp)
        movapd %xmm4,nb134_fjyM(%esp)
        movapd %xmm5,nb134_fjzM(%esp)
        movapd %xmm0,nb134_fixM(%esp)
        movapd %xmm1,nb134_fiyM(%esp)
        movapd %xmm2,nb134_fizM(%esp)

        movl nb134_faction(%ebp),%edi

        ## Did all interactions - now update j forces 
        ## Step1 - transpose fjxO, fjyO and fjzO, fjxH1
        movapd nb134_fjxO(%esp),%xmm0
        movapd nb134_fjyO(%esp),%xmm1
        movapd nb134_fjzO(%esp),%xmm2
        movapd nb134_fjxH1(%esp),%xmm3
        movapd %xmm0,%xmm4
        movapd %xmm2,%xmm5
        unpcklpd %xmm1,%xmm0  ## fjOxA fjOyA
        unpckhpd %xmm1,%xmm4  ## fjOxB fjOyB
        unpcklpd %xmm3,%xmm2  ## fjOzA fjH1xA
        unpckhpd %xmm3,%xmm5  ## fjOzB fjH1xB
        movlpd (%edi,%eax,8),%xmm1
        movlpd (%edi,%ebx,8),%xmm3
        movhpd 8(%edi,%eax,8),%xmm1
        movhpd 8(%edi,%ebx,8),%xmm3
        movlpd 16(%edi,%eax,8),%xmm6
        movlpd 16(%edi,%ebx,8),%xmm7
        movhpd 24(%edi,%eax,8),%xmm6
        movhpd 24(%edi,%ebx,8),%xmm7
        addpd  %xmm0,%xmm1
        addpd  %xmm4,%xmm3
        addpd  %xmm2,%xmm6
        addpd  %xmm5,%xmm7
        movlpd %xmm1,(%edi,%eax,8)
        movlpd %xmm3,(%edi,%ebx,8)
        movhpd %xmm1,8(%edi,%eax,8)
        movhpd %xmm3,8(%edi,%ebx,8)
        movlpd %xmm6,16(%edi,%eax,8)
        movlpd %xmm7,16(%edi,%ebx,8)
        movhpd %xmm6,24(%edi,%eax,8)
        movhpd %xmm7,24(%edi,%ebx,8)

        ## Step2 - transpose fjyH1, fjzH1 and fjxH2, fjyH2
        movapd nb134_fjyH1(%esp),%xmm0
        movapd nb134_fjzH1(%esp),%xmm1
        movapd nb134_fjxH2(%esp),%xmm2
        movapd nb134_fjyH2(%esp),%xmm3
        movapd %xmm0,%xmm4
        movapd %xmm2,%xmm5
        unpcklpd %xmm1,%xmm0  ## fjOxA fjOyA
        unpckhpd %xmm1,%xmm4  ## fjOxB fjOyB
        unpcklpd %xmm3,%xmm2  ## fjOzA fjH1xA
        unpckhpd %xmm3,%xmm5  ## fjOzB fjH1xB
        movlpd 32(%edi,%eax,8),%xmm1
        movlpd 32(%edi,%ebx,8),%xmm3
        movhpd 40(%edi,%eax,8),%xmm1
        movhpd 40(%edi,%ebx,8),%xmm3
        movlpd 48(%edi,%eax,8),%xmm6
        movlpd 48(%edi,%ebx,8),%xmm7
        movhpd 56(%edi,%eax,8),%xmm6
        movhpd 56(%edi,%ebx,8),%xmm7
        addpd  %xmm0,%xmm1
        addpd  %xmm4,%xmm3
        addpd  %xmm2,%xmm6
        addpd  %xmm5,%xmm7
        movlpd %xmm1,32(%edi,%eax,8)
        movlpd %xmm3,32(%edi,%ebx,8)
        movhpd %xmm1,40(%edi,%eax,8)
        movhpd %xmm3,40(%edi,%ebx,8)
        movlpd %xmm6,48(%edi,%eax,8)
        movlpd %xmm7,48(%edi,%ebx,8)
        movhpd %xmm6,56(%edi,%eax,8)
        movhpd %xmm7,56(%edi,%ebx,8)

        ## Step3 - transpose fjzH2, fjxM and fjyM, fjzM
        movapd nb134_fjzH2(%esp),%xmm0
        movapd nb134_fjxM(%esp),%xmm1
        movapd nb134_fjyM(%esp),%xmm2
        movapd nb134_fjzM(%esp),%xmm3
        movapd %xmm0,%xmm4
        movapd %xmm2,%xmm5
        unpcklpd %xmm1,%xmm0  ## fjOxA fjOyA
        unpckhpd %xmm1,%xmm4  ## fjOxB fjOyB
        unpcklpd %xmm3,%xmm2  ## fjOzA fjH1xA
        unpckhpd %xmm3,%xmm5  ## fjOzB fjH1xB
        movlpd 64(%edi,%eax,8),%xmm1
        movlpd 64(%edi,%ebx,8),%xmm3
        movhpd 72(%edi,%eax,8),%xmm1
        movhpd 72(%edi,%ebx,8),%xmm3
        movlpd 80(%edi,%eax,8),%xmm6
        movlpd 80(%edi,%ebx,8),%xmm7
        movhpd 88(%edi,%eax,8),%xmm6
        movhpd 88(%edi,%ebx,8),%xmm7
        addpd  %xmm0,%xmm1
        addpd  %xmm4,%xmm3
        addpd  %xmm2,%xmm6
        addpd  %xmm5,%xmm7
        movlpd %xmm1,64(%edi,%eax,8)
        movlpd %xmm3,64(%edi,%ebx,8)
        movhpd %xmm1,72(%edi,%eax,8)
        movhpd %xmm3,72(%edi,%ebx,8)
        movlpd %xmm6,80(%edi,%eax,8)
        movlpd %xmm7,80(%edi,%ebx,8)
        movhpd %xmm6,88(%edi,%eax,8)
        movhpd %xmm7,88(%edi,%ebx,8)

        ## should we do one more iteration? 
        subl $2,nb134_innerk(%esp)
        jl    _nb_kernel134_ia32_sse2.nb134_checksingle
        jmp   _nb_kernel134_ia32_sse2.nb134_unroll_loop
_nb_kernel134_ia32_sse2.nb134_checksingle: 
        movl  nb134_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel134_ia32_sse2.nb134_dosingle
        jmp   _nb_kernel134_ia32_sse2.nb134_updateouterdata
_nb_kernel134_ia32_sse2.nb134_dosingle: 
        movl  nb134_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb134_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move j coordinates to local temp variables 
        ## load ox, oy, oz, h1x
        movlpd (%esi,%eax,8),%xmm0
        movhpd 8(%esi,%eax,8),%xmm0
        movlpd 16(%esi,%eax,8),%xmm1
        movhpd 24(%esi,%eax,8),%xmm1
        movlpd 32(%esi,%eax,8),%xmm2
        movhpd 40(%esi,%eax,8),%xmm2
        movlpd 48(%esi,%eax,8),%xmm3
        movhpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movhpd 72(%esi,%eax,8),%xmm4
        movlpd 80(%esi,%eax,8),%xmm5
        movhpd 88(%esi,%eax,8),%xmm5
        movsd  %xmm0,nb134_jxO(%esp)
        movsd  %xmm1,nb134_jzO(%esp)
        movsd  %xmm2,nb134_jyH1(%esp)
        movsd  %xmm3,nb134_jxH2(%esp)
        movsd  %xmm4,nb134_jzH2(%esp)
        movsd  %xmm5,nb134_jyM(%esp)
        unpckhpd %xmm0,%xmm0
        unpckhpd %xmm1,%xmm1
        unpckhpd %xmm2,%xmm2
        unpckhpd %xmm3,%xmm3
        unpckhpd %xmm4,%xmm4
        unpckhpd %xmm5,%xmm5
        movsd  %xmm0,nb134_jyO(%esp)
        movsd  %xmm1,nb134_jxH1(%esp)
        movsd  %xmm2,nb134_jzH1(%esp)
        movsd  %xmm3,nb134_jyH2(%esp)
        movsd  %xmm4,nb134_jxM(%esp)
        movsd  %xmm5,nb134_jzM(%esp)

        ## start calculating pairwise distances
        movapd nb134_ixO(%esp),%xmm0
        movapd nb134_iyO(%esp),%xmm1
        movapd nb134_izO(%esp),%xmm2
        movapd nb134_ixH1(%esp),%xmm3
        movapd nb134_iyH1(%esp),%xmm4
        movapd nb134_izH1(%esp),%xmm5
        subsd  nb134_jxO(%esp),%xmm0
        subsd  nb134_jyO(%esp),%xmm1
        subsd  nb134_jzO(%esp),%xmm2
        subsd  nb134_jxH1(%esp),%xmm3
        subsd  nb134_jyH1(%esp),%xmm4
        subsd  nb134_jzH1(%esp),%xmm5
        movapd %xmm0,nb134_dxOO(%esp)
        movapd %xmm1,nb134_dyOO(%esp)
        movapd %xmm2,nb134_dzOO(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb134_dxH1H1(%esp)
        movapd %xmm4,nb134_dyH1H1(%esp)
        movapd %xmm5,nb134_dzH1H1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb134_rsqOO(%esp)
        movapd %xmm3,nb134_rsqH1H1(%esp)

        movapd nb134_ixH1(%esp),%xmm0
        movapd nb134_iyH1(%esp),%xmm1
        movapd nb134_izH1(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb134_jxH2(%esp),%xmm0
        subsd  nb134_jyH2(%esp),%xmm1
        subsd  nb134_jzH2(%esp),%xmm2
        subsd  nb134_jxM(%esp),%xmm3
        subsd  nb134_jyM(%esp),%xmm4
        subsd  nb134_jzM(%esp),%xmm5
        movapd %xmm0,nb134_dxH1H2(%esp)
        movapd %xmm1,nb134_dyH1H2(%esp)
        movapd %xmm2,nb134_dzH1H2(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb134_dxH1M(%esp)
        movapd %xmm4,nb134_dyH1M(%esp)
        movapd %xmm5,nb134_dzH1M(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb134_rsqH1H2(%esp)
        movapd %xmm3,nb134_rsqH1M(%esp)

        movapd nb134_ixH2(%esp),%xmm0
        movapd nb134_iyH2(%esp),%xmm1
        movapd nb134_izH2(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb134_jxH1(%esp),%xmm0
        subsd  nb134_jyH1(%esp),%xmm1
        subsd  nb134_jzH1(%esp),%xmm2
        subsd  nb134_jxH2(%esp),%xmm3
        subsd  nb134_jyH2(%esp),%xmm4
        subsd  nb134_jzH2(%esp),%xmm5
        movapd %xmm0,nb134_dxH2H1(%esp)
        movapd %xmm1,nb134_dyH2H1(%esp)
        movapd %xmm2,nb134_dzH2H1(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb134_dxH2H2(%esp)
        movapd %xmm4,nb134_dyH2H2(%esp)
        movapd %xmm5,nb134_dzH2H2(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb134_rsqH2H1(%esp)
        movapd %xmm3,nb134_rsqH2H2(%esp)

        movapd nb134_ixH2(%esp),%xmm0
        movapd nb134_iyH2(%esp),%xmm1
        movapd nb134_izH2(%esp),%xmm2
        movapd nb134_ixM(%esp),%xmm3
        movapd nb134_iyM(%esp),%xmm4
        movapd nb134_izM(%esp),%xmm5
        subsd  nb134_jxM(%esp),%xmm0
        subsd  nb134_jyM(%esp),%xmm1
        subsd  nb134_jzM(%esp),%xmm2
        subsd  nb134_jxH1(%esp),%xmm3
        subsd  nb134_jyH1(%esp),%xmm4
        subsd  nb134_jzH1(%esp),%xmm5
        movapd %xmm0,nb134_dxH2M(%esp)
        movapd %xmm1,nb134_dyH2M(%esp)
        movapd %xmm2,nb134_dzH2M(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb134_dxMH1(%esp)
        movapd %xmm4,nb134_dyMH1(%esp)
        movapd %xmm5,nb134_dzMH1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb134_rsqH2M(%esp)
        movapd %xmm4,nb134_rsqMH1(%esp)

        movapd nb134_ixM(%esp),%xmm0
        movapd nb134_iyM(%esp),%xmm1
        movapd nb134_izM(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb134_jxH2(%esp),%xmm0
        subsd  nb134_jyH2(%esp),%xmm1
        subsd  nb134_jzH2(%esp),%xmm2
        subsd  nb134_jxM(%esp),%xmm3
        subsd  nb134_jyM(%esp),%xmm4
        subsd  nb134_jzM(%esp),%xmm5
        movapd %xmm0,nb134_dxMH2(%esp)
        movapd %xmm1,nb134_dyMH2(%esp)
        movapd %xmm2,nb134_dzMH2(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb134_dxMM(%esp)
        movapd %xmm4,nb134_dyMM(%esp)
        movapd %xmm5,nb134_dzMM(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb134_rsqMH2(%esp)
        movapd %xmm4,nb134_rsqMM(%esp)

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
        movapd  nb134_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb134_half(%esp),%xmm3   ## iter1 
        mulsd   nb134_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb134_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb134_half(%esp),%xmm1   ## rinv 
        mulsd   nb134_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb134_rinvMH2(%esp)
        movapd %xmm5,nb134_rinvMM(%esp)

        movapd nb134_rsqOO(%esp),%xmm0
        movapd nb134_rsqH1H1(%esp),%xmm4
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
        movapd  nb134_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb134_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb134_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb134_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb134_half(%esp),%xmm1   ## rinv 
        mulsd   nb134_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb134_rinvOO(%esp)
        movapd %xmm5,nb134_rinvH1H1(%esp)

        movapd nb134_rsqH1H2(%esp),%xmm0
        movapd nb134_rsqH1M(%esp),%xmm4
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
        movapd  nb134_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb134_half(%esp),%xmm3   ## iter1 
        mulsd   nb134_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb134_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb134_half(%esp),%xmm1   ## rinv 
        mulsd   nb134_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb134_rinvH1H2(%esp)
        movapd %xmm5,nb134_rinvH1M(%esp)

        movapd nb134_rsqH2H1(%esp),%xmm0
        movapd nb134_rsqH2H2(%esp),%xmm4
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
        movapd  nb134_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb134_half(%esp),%xmm3   ## iter1a 
        mulsd   nb134_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb134_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb134_half(%esp),%xmm1   ## rinv 
        mulsd   nb134_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb134_rinvH2H1(%esp)
        movapd %xmm5,nb134_rinvH2H2(%esp)

        movapd nb134_rsqMH1(%esp),%xmm0
        movapd nb134_rsqH2M(%esp),%xmm4
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
        movapd  nb134_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb134_half(%esp),%xmm3   ## iter1a 
        mulsd   nb134_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb134_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb134_half(%esp),%xmm1   ## rinv 
        mulsd   nb134_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb134_rinvMH1(%esp)
        movapd %xmm5,nb134_rinvH2M(%esp)

        ## start with OO interaction 
        movsd nb134_rinvOO(%esp),%xmm0
        movsd nb134_rsqOO(%esp),%xmm4

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb134_tsc(%esp),%xmm4

        cvttsd2si %xmm4,%ebx    ## mm6 = lu idx 
        cvtsi2sd %ebx,%xmm5
        subsd %xmm5,%xmm4
        movsd %xmm4,%xmm1       ## xmm1=eps 
        movsd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $3,%ebx

        movl nb134_VFtab(%ebp),%esi

        ## dispersion 
        movlpd (%esi,%ebx,8),%xmm4      ## Y1   
        movhpd 8(%esi,%ebx,8),%xmm4     ## Y1 F1        
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%ebx,8),%xmm6    ## G1
        movhpd 24(%esi,%ebx,8),%xmm6    ## G1 H1        
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## dispersion table ready, in xmm4-xmm7         
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb134_two(%esp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb134_c6(%esp),%xmm4
        mulsd  %xmm4,%xmm7       ## fijD 
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addsd  nb134_Vvdwtot(%esp),%xmm5
        xorpd  %xmm3,%xmm3
        mulsd  nb134_tsc(%esp),%xmm7
        subsd  %xmm7,%xmm3
        movsd %xmm3,nb134_fstmp(%esp)
        movsd %xmm5,nb134_Vvdwtot(%esp)

        ## repulsion 
        movlpd 32(%esi,%ebx,8),%xmm4    ## Y1   
        movhpd 40(%esi,%ebx,8),%xmm4    ## Y1 F1        

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 48(%esi,%ebx,8),%xmm6    ## G1
        movhpd 56(%esi,%ebx,8),%xmm6    ## G1 H1        

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 

        ## table ready, in xmm4-xmm7    
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb134_two(%esp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb134_c12(%esp),%xmm4
        mulsd  %xmm4,%xmm7
        mulsd  %xmm4,%xmm5

        addsd  nb134_Vvdwtot(%esp),%xmm5
        movsd nb134_fstmp(%esp),%xmm3
        mulsd  nb134_tsc(%esp),%xmm7
        subsd  %xmm7,%xmm3
        movsd %xmm5,nb134_Vvdwtot(%esp)

        mulsd  %xmm3,%xmm0
        movsd %xmm0,%xmm1
        movsd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb134_dxOO(%esp),%xmm0
        mulsd nb134_dyOO(%esp),%xmm1
        mulsd nb134_dzOO(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb134_fixO(%esp),%xmm0
        addsd nb134_fiyO(%esp),%xmm1
        addsd nb134_fizO(%esp),%xmm2
        movsd %xmm3,nb134_fjxO(%esp)
        movsd %xmm4,nb134_fjyO(%esp)
        movsd %xmm5,nb134_fjzO(%esp)
        movsd %xmm0,nb134_fixO(%esp)
        movsd %xmm1,nb134_fiyO(%esp)
        movsd %xmm2,nb134_fizO(%esp)

        ## H1-H1 interaction 
        movsd  nb134_rinvH1H1(%esp),%xmm0
        movsd  %xmm0,%xmm6      ## xmm6=rinv 
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb134_qqHH(%esp),%xmm6   ## vcoul
        mulsd %xmm6,%xmm0
        addsd  nb134_vctot(%esp),%xmm6

        movsd %xmm0,%xmm1
        movsd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb134_dxH1H1(%esp),%xmm0
        mulsd nb134_dyH1H1(%esp),%xmm1
        mulsd nb134_dzH1H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb134_fixH1(%esp),%xmm0
        addsd nb134_fiyH1(%esp),%xmm1
        addsd nb134_fizH1(%esp),%xmm2
        movsd %xmm3,nb134_fjxH1(%esp)
        movsd %xmm4,nb134_fjyH1(%esp)
        movsd %xmm5,nb134_fjzH1(%esp)
        movsd %xmm0,nb134_fixH1(%esp)
        movsd %xmm1,nb134_fiyH1(%esp)
        movsd %xmm2,nb134_fizH1(%esp)

        ## H1-H2 interaction  
        movsd  nb134_rinvH1H2(%esp),%xmm0
        movsd  %xmm0,%xmm4      ## xmm7=rinv 
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb134_qqHH(%esp),%xmm4   ## vcoul
        mulsd %xmm4,%xmm0
        addsd  %xmm4,%xmm6

        movsd %xmm0,%xmm1
        movsd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb134_dxH1H2(%esp),%xmm0
        mulsd nb134_dyH1H2(%esp),%xmm1
        mulsd nb134_dzH1H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb134_fixH1(%esp),%xmm0
        addsd nb134_fiyH1(%esp),%xmm1
        addsd nb134_fizH1(%esp),%xmm2
        movsd %xmm3,nb134_fjxH2(%esp)
        movsd %xmm4,nb134_fjyH2(%esp)
        movsd %xmm5,nb134_fjzH2(%esp)
        movsd %xmm0,nb134_fixH1(%esp)
        movsd %xmm1,nb134_fiyH1(%esp)
        movsd %xmm2,nb134_fizH1(%esp)

        ## H1-M interaction 
        movsd  nb134_rinvH1M(%esp),%xmm0
        movsd  %xmm0,%xmm4      ## xmm7=rinv 
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb134_qqMH(%esp),%xmm4   ## vcoul
        mulsd %xmm4,%xmm0
        addsd  %xmm4,%xmm6

        movsd %xmm0,%xmm1
        movsd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb134_dxH1M(%esp),%xmm0
        mulsd nb134_dyH1M(%esp),%xmm1
        mulsd nb134_dzH1M(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb134_fixH1(%esp),%xmm0
        addsd nb134_fiyH1(%esp),%xmm1
        addsd nb134_fizH1(%esp),%xmm2
        movsd %xmm3,nb134_fjxM(%esp)
        movsd %xmm4,nb134_fjyM(%esp)
        movsd %xmm5,nb134_fjzM(%esp)
        movsd %xmm0,nb134_fixH1(%esp)
        movsd %xmm1,nb134_fiyH1(%esp)
        movsd %xmm2,nb134_fizH1(%esp)

        ## H2-H1 interaction 
        movsd  nb134_rinvH2H1(%esp),%xmm0
        movsd  %xmm0,%xmm4      ## xmm7=rinv 
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb134_qqHH(%esp),%xmm4   ## vcoul
        mulsd %xmm4,%xmm0
        addsd  %xmm4,%xmm6

        movsd %xmm0,%xmm1
        movsd %xmm0,%xmm2

        movapd nb134_fjxH1(%esp),%xmm3
        movapd nb134_fjyH1(%esp),%xmm4
        movapd nb134_fjzH1(%esp),%xmm5
        mulsd nb134_dxH2H1(%esp),%xmm0
        mulsd nb134_dyH2H1(%esp),%xmm1
        mulsd nb134_dzH2H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb134_fixH2(%esp),%xmm0
        addsd nb134_fiyH2(%esp),%xmm1
        addsd nb134_fizH2(%esp),%xmm2
        movsd %xmm3,nb134_fjxH1(%esp)
        movsd %xmm4,nb134_fjyH1(%esp)
        movsd %xmm5,nb134_fjzH1(%esp)
        movsd %xmm0,nb134_fixH2(%esp)
        movsd %xmm1,nb134_fiyH2(%esp)
        movsd %xmm2,nb134_fizH2(%esp)

        ## H2-H2 interaction 
        movsd  nb134_rinvH2H2(%esp),%xmm0
        movsd  %xmm0,%xmm4      ## xmm7=rinv 
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb134_qqHH(%esp),%xmm4   ## vcoul
        mulsd %xmm4,%xmm0
        addsd  %xmm4,%xmm6

        movsd %xmm0,%xmm1
        movsd %xmm0,%xmm2

        movsd nb134_fjxH2(%esp),%xmm3
        movsd nb134_fjyH2(%esp),%xmm4
        movsd nb134_fjzH2(%esp),%xmm5
        mulsd nb134_dxH2H2(%esp),%xmm0
        mulsd nb134_dyH2H2(%esp),%xmm1
        mulsd nb134_dzH2H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb134_fixH2(%esp),%xmm0
        addsd nb134_fiyH2(%esp),%xmm1
        addsd nb134_fizH2(%esp),%xmm2
        movsd %xmm3,nb134_fjxH2(%esp)
        movsd %xmm4,nb134_fjyH2(%esp)
        movsd %xmm5,nb134_fjzH2(%esp)
        movsd %xmm0,nb134_fixH2(%esp)
        movsd %xmm1,nb134_fiyH2(%esp)
        movsd %xmm2,nb134_fizH2(%esp)

        ## H2-M interaction 
        movsd  nb134_rinvH2M(%esp),%xmm0
        movsd  %xmm0,%xmm4      ## xmm7=rinv 
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb134_qqMH(%esp),%xmm4   ## vcoul
        mulsd %xmm4,%xmm0
        addsd  %xmm4,%xmm6

        movsd %xmm0,%xmm1
        movsd %xmm0,%xmm2

        movapd nb134_fjxM(%esp),%xmm3
        movapd nb134_fjyM(%esp),%xmm4
        movapd nb134_fjzM(%esp),%xmm5
        mulsd nb134_dxH2M(%esp),%xmm0
        mulsd nb134_dyH2M(%esp),%xmm1
        mulsd nb134_dzH2M(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb134_fixH2(%esp),%xmm0
        addsd nb134_fiyH2(%esp),%xmm1
        addsd nb134_fizH2(%esp),%xmm2
        movsd %xmm3,nb134_fjxM(%esp)
        movsd %xmm4,nb134_fjyM(%esp)
        movsd %xmm5,nb134_fjzM(%esp)
        movsd %xmm0,nb134_fixH2(%esp)
        movsd %xmm1,nb134_fiyH2(%esp)
        movsd %xmm2,nb134_fizH2(%esp)

        ## M-H1 interaction 
        movsd  nb134_rinvMH1(%esp),%xmm0
        movsd  %xmm0,%xmm4      ## xmm7=rinv 
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb134_qqMH(%esp),%xmm4   ## vcoul
        mulsd %xmm4,%xmm0
        addsd  %xmm4,%xmm6

        movsd %xmm0,%xmm1
        movsd %xmm0,%xmm2

        movapd nb134_fjxH1(%esp),%xmm3
        movapd nb134_fjyH1(%esp),%xmm4
        movapd nb134_fjzH1(%esp),%xmm5
        mulsd nb134_dxMH1(%esp),%xmm0
        mulsd nb134_dyMH1(%esp),%xmm1
        mulsd nb134_dzMH1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb134_fixM(%esp),%xmm0
        addsd nb134_fiyM(%esp),%xmm1
        addsd nb134_fizM(%esp),%xmm2
        movsd %xmm3,nb134_fjxH1(%esp)
        movsd %xmm4,nb134_fjyH1(%esp)
        movsd %xmm5,nb134_fjzH1(%esp)
        movsd %xmm0,nb134_fixM(%esp)
        movsd %xmm1,nb134_fiyM(%esp)
        movsd %xmm2,nb134_fizM(%esp)

        ## M-H2 interaction 
        movsd  nb134_rinvMH2(%esp),%xmm0
        movsd  %xmm0,%xmm4      ## xmm7=rinv 
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb134_qqMH(%esp),%xmm4   ## vcoul
        mulsd %xmm4,%xmm0
        addsd  %xmm4,%xmm6

        movsd %xmm0,%xmm1
        movsd %xmm0,%xmm2

        movapd nb134_fjxH2(%esp),%xmm3
        movapd nb134_fjyH2(%esp),%xmm4
        movapd nb134_fjzH2(%esp),%xmm5
        mulsd nb134_dxMH2(%esp),%xmm0
        mulsd nb134_dyMH2(%esp),%xmm1
        mulsd nb134_dzMH2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb134_fixM(%esp),%xmm0
        addsd nb134_fiyM(%esp),%xmm1
        addsd nb134_fizM(%esp),%xmm2
        movsd %xmm3,nb134_fjxH2(%esp)
        movsd %xmm4,nb134_fjyH2(%esp)
        movsd %xmm5,nb134_fjzH2(%esp)
        movsd %xmm0,nb134_fixM(%esp)
        movsd %xmm1,nb134_fiyM(%esp)
        movsd %xmm2,nb134_fizM(%esp)

        ## M-M interaction 
        movsd  nb134_rinvMM(%esp),%xmm0
        movsd  %xmm0,%xmm4      ## xmm7=rinv 
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb134_qqMM(%esp),%xmm4   ## vcoul
        mulsd %xmm4,%xmm0
        addsd  %xmm4,%xmm6
        movsd %xmm6,nb134_vctot(%esp)
        movsd %xmm0,%xmm1
        movsd %xmm0,%xmm2

        movapd nb134_fjxM(%esp),%xmm3
        movapd nb134_fjyM(%esp),%xmm4
        movapd nb134_fjzM(%esp),%xmm5
        mulsd nb134_dxMM(%esp),%xmm0
        mulsd nb134_dyMM(%esp),%xmm1
        mulsd nb134_dzMM(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb134_fixM(%esp),%xmm0
        addsd nb134_fiyM(%esp),%xmm1
        addsd nb134_fizM(%esp),%xmm2
        movsd %xmm3,nb134_fjxM(%esp)
        movsd %xmm4,nb134_fjyM(%esp)
        movsd %xmm5,nb134_fjzM(%esp)
        movsd %xmm0,nb134_fixM(%esp)
        movsd %xmm1,nb134_fiyM(%esp)
        movsd %xmm2,nb134_fizM(%esp)

        movl nb134_faction(%ebp),%edi

        ## Did all interactions - now update j forces 
        ## Step1 - merge forces
        movlpd nb134_fjxO(%esp),%xmm0
        movlpd nb134_fjzO(%esp),%xmm1
        movlpd nb134_fjyH1(%esp),%xmm2
        movlpd nb134_fjxH2(%esp),%xmm3
        movlpd nb134_fjzH2(%esp),%xmm4
        movlpd nb134_fjyM(%esp),%xmm5

        movhpd nb134_fjyO(%esp),%xmm0
        movhpd nb134_fjxH1(%esp),%xmm1
        movhpd nb134_fjzH1(%esp),%xmm2
        movhpd nb134_fjyH2(%esp),%xmm3
        movhpd nb134_fjxM(%esp),%xmm4
        movhpd nb134_fjzM(%esp),%xmm5

        movlps (%edi,%eax,8),%xmm6
        movhps 8(%edi,%eax,8),%xmm6
        movlps 16(%edi,%eax,8),%xmm7
        movhps 24(%edi,%eax,8),%xmm7
        addpd  %xmm6,%xmm0
        addpd  %xmm7,%xmm1
        movlps 32(%edi,%eax,8),%xmm6
        movhps 40(%edi,%eax,8),%xmm6
        movlps 48(%edi,%eax,8),%xmm7
        movhps 56(%edi,%eax,8),%xmm7
        addpd  %xmm6,%xmm2
        addpd  %xmm7,%xmm3
        movlps 64(%edi,%eax,8),%xmm6
        movhps 72(%edi,%eax,8),%xmm6
        movlps 80(%edi,%eax,8),%xmm7
        movhps 88(%edi,%eax,8),%xmm7
        addpd  %xmm6,%xmm4
        addpd  %xmm7,%xmm5

        movlpd %xmm0,(%edi,%eax,8)
        movhpd %xmm0,8(%edi,%eax,8)
        movlpd %xmm1,16(%edi,%eax,8)
        movhpd %xmm1,24(%edi,%eax,8)
        movlpd %xmm2,32(%edi,%eax,8)
        movhpd %xmm2,40(%edi,%eax,8)
        movlpd %xmm3,48(%edi,%eax,8)
        movhpd %xmm3,56(%edi,%eax,8)
        movlpd %xmm4,64(%edi,%eax,8)
        movhpd %xmm4,72(%edi,%eax,8)
        movlpd %xmm5,80(%edi,%eax,8)
        movhpd %xmm5,88(%edi,%eax,8)

_nb_kernel134_ia32_sse2.nb134_updateouterdata: 
        movl  nb134_ii3(%esp),%ecx
        movl  nb134_faction(%ebp),%edi
        movl  nb134_fshift(%ebp),%esi
        movl  nb134_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb134_fixO(%esp),%xmm0
        movapd nb134_fiyO(%esp),%xmm1
        movapd nb134_fizO(%esp),%xmm2

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
        movsd  (%edi,%ecx,8),%xmm3
        movsd  8(%edi,%ecx,8),%xmm4
        movsd  16(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,(%edi,%ecx,8)
        movsd  %xmm4,8(%edi,%ecx,8)
        movsd  %xmm5,16(%edi,%ecx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        movapd %xmm0,%xmm6
        movsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm6

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movapd nb134_fixH1(%esp),%xmm0
        movapd nb134_fiyH1(%esp),%xmm1
        movapd nb134_fizH1(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  24(%edi,%ecx,8),%xmm3
        movsd  32(%edi,%ecx,8),%xmm4
        movsd  40(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,24(%edi,%ecx,8)
        movsd  %xmm4,32(%edi,%ecx,8)
        movsd  %xmm5,40(%edi,%ecx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movapd nb134_fixH2(%esp),%xmm0
        movapd nb134_fiyH2(%esp),%xmm1
        movapd nb134_fizH2(%esp),%xmm2

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
        movsd  48(%edi,%ecx,8),%xmm3
        movsd  56(%edi,%ecx,8),%xmm4
        movsd  64(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,48(%edi,%ecx,8)
        movsd  %xmm4,56(%edi,%ecx,8)
        movsd  %xmm5,64(%edi,%ecx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## accumulate Mi forces in xmm0, xmm1, xmm2 
        movapd nb134_fixM(%esp),%xmm0
        movapd nb134_fiyM(%esp),%xmm1
        movapd nb134_fizM(%esp),%xmm2

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
        movsd  72(%edi,%ecx,8),%xmm3
        movsd  80(%edi,%ecx,8),%xmm4
        movsd  88(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,72(%edi,%ecx,8)
        movsd  %xmm4,80(%edi,%ecx,8)
        movsd  %xmm5,88(%edi,%ecx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## increment fshift force 
        movlpd (%esi,%edx,8),%xmm3
        movhpd 8(%esi,%edx,8),%xmm3
        movsd  16(%esi,%edx,8),%xmm4
        addpd  %xmm6,%xmm3
        addsd  %xmm7,%xmm4
        movlpd %xmm3,(%esi,%edx,8)
        movhpd %xmm3,8(%esi,%edx,8)
        movsd  %xmm4,16(%esi,%edx,8)

        ## get n from stack
        movl nb134_n(%esp),%esi
        ## get group index for i particle 
        movl  nb134_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb134_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb134_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb134_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb134_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb134_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel134_ia32_sse2.nb134_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb134_n(%esp)
        jmp _nb_kernel134_ia32_sse2.nb134_outer
_nb_kernel134_ia32_sse2.nb134_outerend: 
        ## check if more outer neighborlists remain
        movl  nb134_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel134_ia32_sse2.nb134_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel134_ia32_sse2.nb134_threadloop
_nb_kernel134_ia32_sse2.nb134_end: 
        emms

        movl nb134_nouter(%esp),%eax
        movl nb134_ninner(%esp),%ebx
        movl nb134_outeriter(%ebp),%ecx
        movl nb134_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb134_salign(%esp),%eax
        addl %eax,%esp
        addl $1848,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret




.globl nb_kernel134nf_ia32_sse2
.globl _nb_kernel134nf_ia32_sse2
nb_kernel134nf_ia32_sse2:       
_nb_kernel134nf_ia32_sse2:      
.set nb134nf_p_nri, 8
.set nb134nf_iinr, 12
.set nb134nf_jindex, 16
.set nb134nf_jjnr, 20
.set nb134nf_shift, 24
.set nb134nf_shiftvec, 28
.set nb134nf_fshift, 32
.set nb134nf_gid, 36
.set nb134nf_pos, 40
.set nb134nf_faction, 44
.set nb134nf_charge, 48
.set nb134nf_p_facel, 52
.set nb134nf_argkrf, 56
.set nb134nf_argcrf, 60
.set nb134nf_Vc, 64
.set nb134nf_type, 68
.set nb134nf_p_ntype, 72
.set nb134nf_vdwparam, 76
.set nb134nf_Vvdw, 80
.set nb134nf_p_tabscale, 84
.set nb134nf_VFtab, 88
.set nb134nf_invsqrta, 92
.set nb134nf_dvda, 96
.set nb134nf_p_gbtabscale, 100
.set nb134nf_GBtab, 104
.set nb134nf_p_nthreads, 108
.set nb134nf_count, 112
.set nb134nf_mtx, 116
.set nb134nf_outeriter, 120
.set nb134nf_inneriter, 124
.set nb134nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb134nf_ixO, 0
.set nb134nf_iyO, 16
.set nb134nf_izO, 32
.set nb134nf_ixH1, 48
.set nb134nf_iyH1, 64
.set nb134nf_izH1, 80
.set nb134nf_ixH2, 96
.set nb134nf_iyH2, 112
.set nb134nf_izH2, 128
.set nb134nf_ixM, 144
.set nb134nf_iyM, 160
.set nb134nf_izM, 176
.set nb134nf_jxO, 192
.set nb134nf_jyO, 208
.set nb134nf_jzO, 224
.set nb134nf_jxH1, 240
.set nb134nf_jyH1, 256
.set nb134nf_jzH1, 272
.set nb134nf_jxH2, 288
.set nb134nf_jyH2, 304
.set nb134nf_jzH2, 320
.set nb134nf_jxM, 336
.set nb134nf_jyM, 352
.set nb134nf_jzM, 368
.set nb134nf_qqMM, 384
.set nb134nf_qqMH, 400
.set nb134nf_qqHH, 416
.set nb134nf_two, 432
.set nb134nf_c6, 448
.set nb134nf_c12, 464
.set nb134nf_vctot, 480
.set nb134nf_Vvdwtot, 496
.set nb134nf_half, 512
.set nb134nf_three, 528
.set nb134nf_tsc, 544
.set nb134nf_rsqOO, 560
.set nb134nf_rsqH1H1, 576
.set nb134nf_rsqH1H2, 592
.set nb134nf_rsqH1M, 608
.set nb134nf_rsqH2H1, 624
.set nb134nf_rsqH2H2, 640
.set nb134nf_rsqH2M, 656
.set nb134nf_rsqMH1, 672
.set nb134nf_rsqMH2, 688
.set nb134nf_rsqMM, 704
.set nb134nf_rinvOO, 720
.set nb134nf_rinvH1H1, 736
.set nb134nf_rinvH1H2, 752
.set nb134nf_rinvH1M, 768
.set nb134nf_rinvH2H1, 784
.set nb134nf_rinvH2H2, 800
.set nb134nf_rinvH2M, 816
.set nb134nf_rinvMH1, 832
.set nb134nf_rinvMH2, 848
.set nb134nf_rinvMM, 864
.set nb134nf_is3, 912
.set nb134nf_ii3, 916
.set nb134nf_innerjjnr, 920
.set nb134nf_innerk, 924
.set nb134nf_n, 928
.set nb134nf_nn1, 932
.set nb134nf_nri, 936
.set nb134nf_nouter, 940
.set nb134nf_ninner, 944
.set nb134nf_salign, 948
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $952,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb134nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb134nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb134nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb134nf_nouter(%esp)
        movl %eax,nb134nf_ninner(%esp)

        movl nb134nf_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb134nf_tsc(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb134nf_half(%esp)
        movl %ebx,nb134nf_half+4(%esp)
        movsd nb134nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb134nf_half(%esp)
        movapd %xmm2,nb134nf_two(%esp)
        movapd %xmm3,nb134nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb134nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb134nf_charge(%ebp),%edx
        movsd 24(%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb134nf_p_facel(%ebp),%esi
        movsd (%esi),%xmm6
        mulsd  %xmm3,%xmm3
        mulsd  %xmm5,%xmm4
        mulsd  %xmm5,%xmm5
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb134nf_qqMM(%esp)
        movapd %xmm4,nb134nf_qqMH(%esp)
        movapd %xmm5,nb134nf_qqHH(%esp)

        xorpd %xmm0,%xmm0
        movl  nb134nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb134nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb134nf_vdwparam(%ebp),%eax
        movlpd (%eax,%edx,8),%xmm0
        movhpd 8(%eax,%edx,8),%xmm0
        movhlps %xmm0,%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb134nf_c6(%esp)
        movapd %xmm1,nb134nf_c12(%esp)

_nb_kernel134nf_ia32_sse2.nb134nf_threadloop: 
        movl  nb134nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel134nf_ia32_sse2.nb134nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel134nf_ia32_sse2.nb134nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb134nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb134nf_n(%esp)
        movl %ebx,nb134nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel134nf_ia32_sse2.nb134nf_outerstart
        jmp _nb_kernel134nf_ia32_sse2.nb134nf_end

_nb_kernel134nf_ia32_sse2.nb134nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb134nf_nouter(%esp),%ebx
        movl %ebx,nb134nf_nouter(%esp)

_nb_kernel134nf_ia32_sse2.nb134nf_outer: 
        movl  nb134nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb134nf_is3(%esp)            ## store is3 

        movl  nb134nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb134nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb134nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb134nf_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3       ## ox
        addsd 8(%eax,%ebx,8),%xmm4      ## oy
        addsd 16(%eax,%ebx,8),%xmm5     ## oz   
        addsd 24(%eax,%ebx,8),%xmm6     ## h1x
        addsd 32(%eax,%ebx,8),%xmm7     ## h1y
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm7,%xmm7
        movapd %xmm3,nb134nf_ixO(%esp)
        movapd %xmm4,nb134nf_iyO(%esp)
        movapd %xmm5,nb134nf_izO(%esp)
        movapd %xmm6,nb134nf_ixH1(%esp)
        movapd %xmm7,nb134nf_iyH1(%esp)

        movsd %xmm2,%xmm6
        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 40(%eax,%ebx,8),%xmm6    ## h1z
        addsd 48(%eax,%ebx,8),%xmm0    ## h2x
        addsd 56(%eax,%ebx,8),%xmm1    ## h2y
        addsd 64(%eax,%ebx,8),%xmm2    ## h2z
        addsd 72(%eax,%ebx,8),%xmm3    ## mx
        addsd 80(%eax,%ebx,8),%xmm4    ## my
        addsd 88(%eax,%ebx,8),%xmm5    ## mz

        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm6,nb134nf_izH1(%esp)
        movapd %xmm0,nb134nf_ixH2(%esp)
        movapd %xmm1,nb134nf_iyH2(%esp)
        movapd %xmm2,nb134nf_izH2(%esp)
        movapd %xmm3,nb134nf_ixM(%esp)
        movapd %xmm4,nb134nf_iyM(%esp)
        movapd %xmm5,nb134nf_izM(%esp)

        ## clear vctot
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb134nf_vctot(%esp)
        movapd %xmm4,nb134nf_Vvdwtot(%esp)

        movl  nb134nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb134nf_pos(%ebp),%esi
        movl  nb134nf_faction(%ebp),%edi
        movl  nb134nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb134nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb134nf_ninner(%esp),%ecx
        movl  %ecx,nb134nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb134nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel134nf_ia32_sse2.nb134nf_unroll_loop
        jmp   _nb_kernel134nf_ia32_sse2.nb134nf_checksingle
_nb_kernel134nf_ia32_sse2.nb134nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb134nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb134nf_innerjjnr(%esp)            ## advance pointer (unrolled 2) 

        movl nb134nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move j coordinates to local temp variables 
        ## load ox, oy, oz, h1x
        movlpd (%esi,%eax,8),%xmm0
        movlpd (%esi,%ebx,8),%xmm2
        movhpd 8(%esi,%eax,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm2
        movlpd 16(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%ebx,8),%xmm5
        movhpd 24(%esi,%eax,8),%xmm3
        movhpd 24(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## ox 
        unpckhpd %xmm2,%xmm1 ## oy
        unpcklpd %xmm5,%xmm3 ## ox 
        unpckhpd %xmm5,%xmm4 ## oy
        movapd  %xmm0,nb134nf_jxO(%esp)
        movapd  %xmm1,nb134nf_jyO(%esp)
        movapd  %xmm3,nb134nf_jzO(%esp)
        movapd  %xmm4,nb134nf_jxH1(%esp)

        ## load h1y, h1z, h2x, h2y 
        movlpd 32(%esi,%eax,8),%xmm0
        movlpd 32(%esi,%ebx,8),%xmm2
        movhpd 40(%esi,%eax,8),%xmm0
        movhpd 40(%esi,%ebx,8),%xmm2
        movlpd 48(%esi,%eax,8),%xmm3
        movlpd 48(%esi,%ebx,8),%xmm5
        movhpd 56(%esi,%eax,8),%xmm3
        movhpd 56(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## h1y
        unpckhpd %xmm2,%xmm1 ## h1z
        unpcklpd %xmm5,%xmm3 ## h2x
        unpckhpd %xmm5,%xmm4 ## h2y
        movapd  %xmm0,nb134nf_jyH1(%esp)
        movapd  %xmm1,nb134nf_jzH1(%esp)
        movapd  %xmm3,nb134nf_jxH2(%esp)
        movapd  %xmm4,nb134nf_jyH2(%esp)

        ## load h2z, mx, my, mz
        movlpd 64(%esi,%eax,8),%xmm0
        movlpd 64(%esi,%ebx,8),%xmm2
        movhpd 72(%esi,%eax,8),%xmm0
        movhpd 72(%esi,%ebx,8),%xmm2
        movlpd 80(%esi,%eax,8),%xmm3
        movlpd 80(%esi,%ebx,8),%xmm5
        movhpd 88(%esi,%eax,8),%xmm3
        movhpd 88(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## h2z
        unpckhpd %xmm2,%xmm1 ## mx
        unpcklpd %xmm5,%xmm3 ## my
        unpckhpd %xmm5,%xmm4 ## mz
        movapd  %xmm0,nb134nf_jzH2(%esp)
        movapd  %xmm1,nb134nf_jxM(%esp)
        movapd  %xmm3,nb134nf_jyM(%esp)
        movapd  %xmm4,nb134nf_jzM(%esp)

        ## start calculating pairwise distances
        movapd nb134nf_ixO(%esp),%xmm0
        movapd nb134nf_iyO(%esp),%xmm1
        movapd nb134nf_izO(%esp),%xmm2
        movapd nb134nf_ixH1(%esp),%xmm3
        movapd nb134nf_iyH1(%esp),%xmm4
        movapd nb134nf_izH1(%esp),%xmm5
        subpd  nb134nf_jxO(%esp),%xmm0
        subpd  nb134nf_jyO(%esp),%xmm1
        subpd  nb134nf_jzO(%esp),%xmm2
        subpd  nb134nf_jxH1(%esp),%xmm3
        subpd  nb134nf_jyH1(%esp),%xmm4
        subpd  nb134nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb134nf_rsqOO(%esp)
        movapd %xmm3,nb134nf_rsqH1H1(%esp)

        movapd nb134nf_ixH1(%esp),%xmm0
        movapd nb134nf_iyH1(%esp),%xmm1
        movapd nb134nf_izH1(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb134nf_jxH2(%esp),%xmm0
        subpd  nb134nf_jyH2(%esp),%xmm1
        subpd  nb134nf_jzH2(%esp),%xmm2
        subpd  nb134nf_jxM(%esp),%xmm3
        subpd  nb134nf_jyM(%esp),%xmm4
        subpd  nb134nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb134nf_rsqH1H2(%esp)
        movapd %xmm3,nb134nf_rsqH1M(%esp)

        movapd nb134nf_ixH2(%esp),%xmm0
        movapd nb134nf_iyH2(%esp),%xmm1
        movapd nb134nf_izH2(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb134nf_jxH1(%esp),%xmm0
        subpd  nb134nf_jyH1(%esp),%xmm1
        subpd  nb134nf_jzH1(%esp),%xmm2
        subpd  nb134nf_jxH2(%esp),%xmm3
        subpd  nb134nf_jyH2(%esp),%xmm4
        subpd  nb134nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb134nf_rsqH2H1(%esp)
        movapd %xmm3,nb134nf_rsqH2H2(%esp)

        movapd nb134nf_ixH2(%esp),%xmm0
        movapd nb134nf_iyH2(%esp),%xmm1
        movapd nb134nf_izH2(%esp),%xmm2
        movapd nb134nf_ixM(%esp),%xmm3
        movapd nb134nf_iyM(%esp),%xmm4
        movapd nb134nf_izM(%esp),%xmm5
        subpd  nb134nf_jxM(%esp),%xmm0
        subpd  nb134nf_jyM(%esp),%xmm1
        subpd  nb134nf_jzM(%esp),%xmm2
        subpd  nb134nf_jxH1(%esp),%xmm3
        subpd  nb134nf_jyH1(%esp),%xmm4
        subpd  nb134nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb134nf_rsqH2M(%esp)
        movapd %xmm4,nb134nf_rsqMH1(%esp)

        movapd nb134nf_ixM(%esp),%xmm0
        movapd nb134nf_iyM(%esp),%xmm1
        movapd nb134nf_izM(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb134nf_jxH2(%esp),%xmm0
        subpd  nb134nf_jyH2(%esp),%xmm1
        subpd  nb134nf_jzH2(%esp),%xmm2
        subpd  nb134nf_jxM(%esp),%xmm3
        subpd  nb134nf_jyM(%esp),%xmm4
        subpd  nb134nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb134nf_rsqMH2(%esp)
        movapd %xmm4,nb134nf_rsqMM(%esp)

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
        movapd  nb134nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb134nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb134nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb134nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb134nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb134nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb134nf_rinvMH2(%esp)
        movapd %xmm5,nb134nf_rinvMM(%esp)

        movapd nb134nf_rsqOO(%esp),%xmm0
        movapd nb134nf_rsqH1H1(%esp),%xmm4
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
        movapd  nb134nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb134nf_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb134nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb134nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb134nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb134nf_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb134nf_rinvOO(%esp)
        movapd %xmm5,nb134nf_rinvH1H1(%esp)

        movapd nb134nf_rsqH1H2(%esp),%xmm0
        movapd nb134nf_rsqH1M(%esp),%xmm4
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
        movapd  nb134nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb134nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb134nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb134nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb134nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb134nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb134nf_rinvH1H2(%esp)
        movapd %xmm5,nb134nf_rinvH1M(%esp)

        movapd nb134nf_rsqH2H1(%esp),%xmm0
        movapd nb134nf_rsqH2H2(%esp),%xmm4
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
        movapd  nb134nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb134nf_half(%esp),%xmm3   ## iter1a 
        mulpd   nb134nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb134nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb134nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb134nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb134nf_rinvH2H1(%esp)
        movapd %xmm5,nb134nf_rinvH2H2(%esp)

        movapd nb134nf_rsqMH1(%esp),%xmm0
        movapd nb134nf_rsqH2M(%esp),%xmm4
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
        movapd  nb134nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb134nf_half(%esp),%xmm3   ## iter1a 
        mulpd   nb134nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb134nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb134nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb134nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb134nf_rinvMH1(%esp)
        movapd %xmm5,nb134nf_rinvH2M(%esp)

        ## start with OO interaction 
        movapd nb134nf_rinvOO(%esp),%xmm0
        movapd nb134nf_rsqOO(%esp),%xmm4

                mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb134nf_tsc(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movd %eax,%mm0
        movd %ebx,%mm1

        movl nb134nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx

        ## dispersion 
        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2 
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 
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

        movapd nb134nf_c6(%esp),%xmm4
        mulpd  %xmm4,%xmm5       ## Vvdw6 

        ## Update Vvdwtot directly 
        addpd  nb134nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb134nf_Vvdwtot(%esp)

        ## repulsion 
        movlpd 32(%esi,%eax,8),%xmm4    ## Y1   
        movlpd 32(%esi,%ebx,8),%xmm3    ## Y2 
        movhpd 40(%esi,%eax,8),%xmm4    ## Y1 F1        
        movhpd 40(%esi,%ebx,8),%xmm3    ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 48(%esi,%eax,8),%xmm6    ## G1
        movlpd 48(%esi,%ebx,8),%xmm3    ## G2
        movhpd 56(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 56(%esi,%ebx,8),%xmm3    ## G2 H2 

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

        movapd nb134nf_c12(%esp),%xmm4
        mulpd  %xmm4,%xmm5

        addpd  nb134nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb134nf_Vvdwtot(%esp)

        ## All Coulomb interactions
        movapd nb134nf_rinvH1H1(%esp),%xmm0
        movapd nb134nf_rinvH1M(%esp),%xmm1
        addpd  nb134nf_rinvH1H2(%esp),%xmm0
        addpd  nb134nf_rinvH2M(%esp),%xmm1
        addpd  nb134nf_rinvH2H1(%esp),%xmm0
        addpd  nb134nf_rinvMH1(%esp),%xmm1
        addpd  nb134nf_rinvH2H2(%esp),%xmm0
        addpd  nb134nf_rinvMH2(%esp),%xmm1
        movapd nb134nf_rinvMM(%esp),%xmm2

        mulpd  nb134nf_qqHH(%esp),%xmm0
        mulpd  nb134nf_qqMH(%esp),%xmm1
        mulpd  nb134nf_qqMM(%esp),%xmm2
        addpd  %xmm1,%xmm0
        addpd  nb134nf_vctot(%esp),%xmm2
        addpd  %xmm2,%xmm0
        movapd %xmm0,nb134nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb134nf_innerk(%esp)
        jl    _nb_kernel134nf_ia32_sse2.nb134nf_checksingle
        jmp   _nb_kernel134nf_ia32_sse2.nb134nf_unroll_loop
_nb_kernel134nf_ia32_sse2.nb134nf_checksingle: 
        movl  nb134nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel134nf_ia32_sse2.nb134nf_dosingle
        jmp   _nb_kernel134nf_ia32_sse2.nb134nf_updateouterdata
_nb_kernel134nf_ia32_sse2.nb134nf_dosingle: 
        movl  nb134nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb134nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move j coordinates to local temp variables 
        ## load ox, oy, oz, h1x
        movlpd (%esi,%eax,8),%xmm0
        movhpd 8(%esi,%eax,8),%xmm0
        movlpd 16(%esi,%eax,8),%xmm1
        movhpd 24(%esi,%eax,8),%xmm1
        movlpd 32(%esi,%eax,8),%xmm2
        movhpd 40(%esi,%eax,8),%xmm2
        movlpd 48(%esi,%eax,8),%xmm3
        movhpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movhpd 72(%esi,%eax,8),%xmm4
        movlpd 80(%esi,%eax,8),%xmm5
        movhpd 88(%esi,%eax,8),%xmm5
        movsd  %xmm0,nb134nf_jxO(%esp)
        movsd  %xmm1,nb134nf_jzO(%esp)
        movsd  %xmm2,nb134nf_jyH1(%esp)
        movsd  %xmm3,nb134nf_jxH2(%esp)
        movsd  %xmm4,nb134nf_jzH2(%esp)
        movsd  %xmm5,nb134nf_jyM(%esp)
        unpckhpd %xmm0,%xmm0
        unpckhpd %xmm1,%xmm1
        unpckhpd %xmm2,%xmm2
        unpckhpd %xmm3,%xmm3
        unpckhpd %xmm4,%xmm4
        unpckhpd %xmm5,%xmm5
        movsd  %xmm0,nb134nf_jyO(%esp)
        movsd  %xmm1,nb134nf_jxH1(%esp)
        movsd  %xmm2,nb134nf_jzH1(%esp)
        movsd  %xmm3,nb134nf_jyH2(%esp)
        movsd  %xmm4,nb134nf_jxM(%esp)
        movsd  %xmm5,nb134nf_jzM(%esp)

        ## start calculating pairwise distances
        movapd nb134nf_ixO(%esp),%xmm0
        movapd nb134nf_iyO(%esp),%xmm1
        movapd nb134nf_izO(%esp),%xmm2
        movapd nb134nf_ixH1(%esp),%xmm3
        movapd nb134nf_iyH1(%esp),%xmm4
        movapd nb134nf_izH1(%esp),%xmm5
        subsd  nb134nf_jxO(%esp),%xmm0
        subsd  nb134nf_jyO(%esp),%xmm1
        subsd  nb134nf_jzO(%esp),%xmm2
        subsd  nb134nf_jxH1(%esp),%xmm3
        subsd  nb134nf_jyH1(%esp),%xmm4
        subsd  nb134nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb134nf_rsqOO(%esp)
        movapd %xmm3,nb134nf_rsqH1H1(%esp)

        movapd nb134nf_ixH1(%esp),%xmm0
        movapd nb134nf_iyH1(%esp),%xmm1
        movapd nb134nf_izH1(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb134nf_jxH2(%esp),%xmm0
        subsd  nb134nf_jyH2(%esp),%xmm1
        subsd  nb134nf_jzH2(%esp),%xmm2
        subsd  nb134nf_jxM(%esp),%xmm3
        subsd  nb134nf_jyM(%esp),%xmm4
        subsd  nb134nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb134nf_rsqH1H2(%esp)
        movapd %xmm3,nb134nf_rsqH1M(%esp)

        movapd nb134nf_ixH2(%esp),%xmm0
        movapd nb134nf_iyH2(%esp),%xmm1
        movapd nb134nf_izH2(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb134nf_jxH1(%esp),%xmm0
        subsd  nb134nf_jyH1(%esp),%xmm1
        subsd  nb134nf_jzH1(%esp),%xmm2
        subsd  nb134nf_jxH2(%esp),%xmm3
        subsd  nb134nf_jyH2(%esp),%xmm4
        subsd  nb134nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb134nf_rsqH2H1(%esp)
        movapd %xmm3,nb134nf_rsqH2H2(%esp)

        movapd nb134nf_ixH2(%esp),%xmm0
        movapd nb134nf_iyH2(%esp),%xmm1
        movapd nb134nf_izH2(%esp),%xmm2
        movapd nb134nf_ixM(%esp),%xmm3
        movapd nb134nf_iyM(%esp),%xmm4
        movapd nb134nf_izM(%esp),%xmm5
        subsd  nb134nf_jxM(%esp),%xmm0
        subsd  nb134nf_jyM(%esp),%xmm1
        subsd  nb134nf_jzM(%esp),%xmm2
        subsd  nb134nf_jxH1(%esp),%xmm3
        subsd  nb134nf_jyH1(%esp),%xmm4
        subsd  nb134nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb134nf_rsqH2M(%esp)
        movapd %xmm4,nb134nf_rsqMH1(%esp)

        movapd nb134nf_ixM(%esp),%xmm0
        movapd nb134nf_iyM(%esp),%xmm1
        movapd nb134nf_izM(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb134nf_jxH2(%esp),%xmm0
        subsd  nb134nf_jyH2(%esp),%xmm1
        subsd  nb134nf_jzH2(%esp),%xmm2
        subsd  nb134nf_jxM(%esp),%xmm3
        subsd  nb134nf_jyM(%esp),%xmm4
        subsd  nb134nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb134nf_rsqMH2(%esp)
        movapd %xmm4,nb134nf_rsqMM(%esp)

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
        movapd  nb134nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb134nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb134nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb134nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb134nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb134nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb134nf_rinvMH2(%esp)
        movapd %xmm5,nb134nf_rinvMM(%esp)

        movapd nb134nf_rsqOO(%esp),%xmm0
        movapd nb134nf_rsqH1H1(%esp),%xmm4
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
        movapd  nb134nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb134nf_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb134nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb134nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb134nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb134nf_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb134nf_rinvOO(%esp)
        movapd %xmm5,nb134nf_rinvH1H1(%esp)

        movapd nb134nf_rsqH1H2(%esp),%xmm0
        movapd nb134nf_rsqH1M(%esp),%xmm4
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
        movapd  nb134nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb134nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb134nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb134nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb134nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb134nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb134nf_rinvH1H2(%esp)
        movapd %xmm5,nb134nf_rinvH1M(%esp)

        movapd nb134nf_rsqH2H1(%esp),%xmm0
        movapd nb134nf_rsqH2H2(%esp),%xmm4
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
        movapd  nb134nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb134nf_half(%esp),%xmm3   ## iter1a 
        mulsd   nb134nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb134nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb134nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb134nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb134nf_rinvH2H1(%esp)
        movapd %xmm5,nb134nf_rinvH2H2(%esp)

        movapd nb134nf_rsqMH1(%esp),%xmm0
        movapd nb134nf_rsqH2M(%esp),%xmm4
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
        movapd  nb134nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb134nf_half(%esp),%xmm3   ## iter1a 
        mulsd   nb134nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb134nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb134nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb134nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb134nf_rinvMH1(%esp)
        movapd %xmm5,nb134nf_rinvH2M(%esp)

        ## start with OO interaction 
        movsd nb134nf_rinvOO(%esp),%xmm0
        movsd nb134nf_rsqOO(%esp),%xmm4

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb134nf_tsc(%esp),%xmm4

        cvttsd2si %xmm4,%ebx    ## mm6 = lu idx 
        cvtsi2sd %ebx,%xmm5
        subsd %xmm5,%xmm4
        movsd %xmm4,%xmm1       ## xmm1=eps 
        movsd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $3,%ebx

        movl nb134nf_VFtab(%ebp),%esi

        ## dispersion 
        movlpd (%esi,%ebx,8),%xmm4      ## Y1   
        movhpd 8(%esi,%ebx,8),%xmm4     ## Y1 F1        
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%ebx,8),%xmm6    ## G1
        movhpd 24(%esi,%ebx,8),%xmm6    ## G1 H1        
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

        movsd nb134nf_c6(%esp),%xmm4
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        ## Update Vvdwtot directly 
        addsd  nb134nf_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb134nf_Vvdwtot(%esp)

        ## repulsion 
        movlpd 32(%esi,%ebx,8),%xmm4    ## Y1   
        movhpd 40(%esi,%ebx,8),%xmm4    ## Y1 F1        

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 48(%esi,%ebx,8),%xmm6    ## G1
        movhpd 56(%esi,%ebx,8),%xmm6    ## G1 H1        

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

        movsd nb134nf_c12(%esp),%xmm4
        mulsd  %xmm4,%xmm5

        addsd  nb134nf_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb134nf_Vvdwtot(%esp)

        ## All Coulomb interactions
        movsd nb134nf_rinvH1H1(%esp),%xmm0
        movsd nb134nf_rinvH1M(%esp),%xmm1
        addsd  nb134nf_rinvH1H2(%esp),%xmm0
        addsd  nb134nf_rinvH2M(%esp),%xmm1
        addsd  nb134nf_rinvH2H1(%esp),%xmm0
        addsd  nb134nf_rinvMH1(%esp),%xmm1
        addsd  nb134nf_rinvH2H2(%esp),%xmm0
        addsd  nb134nf_rinvMH2(%esp),%xmm1
        movsd nb134nf_rinvMM(%esp),%xmm2

        mulsd  nb134nf_qqHH(%esp),%xmm0
        mulsd  nb134nf_qqMH(%esp),%xmm1
        mulsd  nb134nf_qqMM(%esp),%xmm2
        addsd  %xmm1,%xmm0
        addsd  nb134nf_vctot(%esp),%xmm2
        addsd  %xmm2,%xmm0
        movsd %xmm0,nb134nf_vctot(%esp)

_nb_kernel134nf_ia32_sse2.nb134nf_updateouterdata: 
        ## get n from stack
        movl nb134nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb134nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb134nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb134nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb134nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb134nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb134nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel134nf_ia32_sse2.nb134nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb134nf_n(%esp)
        jmp _nb_kernel134nf_ia32_sse2.nb134nf_outer
_nb_kernel134nf_ia32_sse2.nb134nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb134nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel134nf_ia32_sse2.nb134nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel134nf_ia32_sse2.nb134nf_threadloop
_nb_kernel134nf_ia32_sse2.nb134nf_end: 
        emms

        movl nb134nf_nouter(%esp),%eax
        movl nb134nf_ninner(%esp),%ebx
        movl nb134nf_outeriter(%ebp),%ecx
        movl nb134nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb134nf_salign(%esp),%eax
        addl %eax,%esp
        addl $952,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



