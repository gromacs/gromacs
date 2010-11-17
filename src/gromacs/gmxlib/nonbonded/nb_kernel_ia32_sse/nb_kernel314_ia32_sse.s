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



.globl nb_kernel314_ia32_sse
.globl _nb_kernel314_ia32_sse
nb_kernel314_ia32_sse:  
_nb_kernel314_ia32_sse: 
.set nb314_p_nri, 8
.set nb314_iinr, 12
.set nb314_jindex, 16
.set nb314_jjnr, 20
.set nb314_shift, 24
.set nb314_shiftvec, 28
.set nb314_fshift, 32
.set nb314_gid, 36
.set nb314_pos, 40
.set nb314_faction, 44
.set nb314_charge, 48
.set nb314_p_facel, 52
.set nb314_argkrf, 56
.set nb314_argcrf, 60
.set nb314_Vc, 64
.set nb314_type, 68
.set nb314_p_ntype, 72
.set nb314_vdwparam, 76
.set nb314_Vvdw, 80
.set nb314_p_tabscale, 84
.set nb314_VFtab, 88
.set nb314_invsqrta, 92
.set nb314_dvda, 96
.set nb314_p_gbtabscale, 100
.set nb314_GBtab, 104
.set nb314_p_nthreads, 108
.set nb314_count, 112
.set nb314_mtx, 116
.set nb314_outeriter, 120
.set nb314_inneriter, 124
.set nb314_work, 128
## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb314_ixO, 0
.set nb314_iyO, 16
.set nb314_izO, 32
.set nb314_ixH1, 48
.set nb314_iyH1, 64
.set nb314_izH1, 80
.set nb314_ixH2, 96
.set nb314_iyH2, 112
.set nb314_izH2, 128
.set nb314_ixM, 144
.set nb314_iyM, 160
.set nb314_izM, 176
.set nb314_jxO, 192
.set nb314_jyO, 208
.set nb314_jzO, 224
.set nb314_jxH1, 240
.set nb314_jyH1, 256
.set nb314_jzH1, 272
.set nb314_jxH2, 288
.set nb314_jyH2, 304
.set nb314_jzH2, 320
.set nb314_jxM, 336
.set nb314_jyM, 352
.set nb314_jzM, 368
.set nb314_dxOO, 384
.set nb314_dyOO, 400
.set nb314_dzOO, 416
.set nb314_dxH1H1, 432
.set nb314_dyH1H1, 448
.set nb314_dzH1H1, 464
.set nb314_dxH1H2, 480
.set nb314_dyH1H2, 496
.set nb314_dzH1H2, 512
.set nb314_dxH1M, 528
.set nb314_dyH1M, 544
.set nb314_dzH1M, 560
.set nb314_dxH2H1, 576
.set nb314_dyH2H1, 592
.set nb314_dzH2H1, 608
.set nb314_dxH2H2, 624
.set nb314_dyH2H2, 640
.set nb314_dzH2H2, 656
.set nb314_dxH2M, 672
.set nb314_dyH2M, 688
.set nb314_dzH2M, 704
.set nb314_dxMH1, 720
.set nb314_dyMH1, 736
.set nb314_dzMH1, 752
.set nb314_dxMH2, 768
.set nb314_dyMH2, 784
.set nb314_dzMH2, 800
.set nb314_dxMM, 816
.set nb314_dyMM, 832
.set nb314_dzMM, 848
.set nb314_qqMM, 864
.set nb314_qqMH, 880
.set nb314_qqHH, 896
.set nb314_two, 912
.set nb314_tsc, 928
.set nb314_c6, 944
.set nb314_c12, 960
.set nb314_six, 976
.set nb314_twelve, 992
.set nb314_vctot, 1008
.set nb314_Vvdwtot, 1024
.set nb314_fixO, 1040
.set nb314_fiyO, 1056
.set nb314_fizO, 1072
.set nb314_fixH1, 1088
.set nb314_fiyH1, 1104
.set nb314_fizH1, 1120
.set nb314_fixH2, 1136
.set nb314_fiyH2, 1152
.set nb314_fizH2, 1168
.set nb314_fixM, 1184
.set nb314_fiyM, 1200
.set nb314_fizM, 1216
.set nb314_fjxO, 1232
.set nb314_fjyO, 1248
.set nb314_fjzO, 1264
.set nb314_fjxH1, 1280
.set nb314_fjyH1, 1296
.set nb314_fjzH1, 1312
.set nb314_fjxH2, 1328
.set nb314_fjyH2, 1344
.set nb314_fjzH2, 1360
.set nb314_fjxM, 1376
.set nb314_fjyM, 1392
.set nb314_fjzM, 1408
.set nb314_half, 1424
.set nb314_three, 1440
.set nb314_rsqOO, 1456
.set nb314_rsqH1H1, 1472
.set nb314_rsqH1H2, 1488
.set nb314_rsqH1M, 1504
.set nb314_rsqH2H1, 1520
.set nb314_rsqH2H2, 1536
.set nb314_rsqH2M, 1552
.set nb314_rsqMH1, 1568
.set nb314_rsqMH2, 1584
.set nb314_rsqMM, 1600
.set nb314_rinvsqOO, 1616
.set nb314_rinvH1H1, 1632
.set nb314_rinvH1H2, 1648
.set nb314_rinvH1M, 1664
.set nb314_rinvH2H1, 1680
.set nb314_rinvH2H2, 1696
.set nb314_rinvH2M, 1712
.set nb314_rinvMH1, 1728
.set nb314_rinvMH2, 1744
.set nb314_rinvMM, 1760
.set nb314_fstmp, 1776
.set nb314_is3, 1792
.set nb314_ii3, 1796
.set nb314_innerjjnr, 1800
.set nb314_innerk, 1804
.set nb314_n, 1808
.set nb314_nn1, 1812
.set nb314_nri, 1816
.set nb314_nouter, 1820
.set nb314_ninner, 1824
.set nb314_salign, 1828
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $1832,%esp         ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb314_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb314_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb314_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb314_nouter(%esp)
        movl %eax,nb314_ninner(%esp)


        movl nb314_p_tabscale(%ebp),%eax
        movss (%eax),%xmm5
        shufps $0,%xmm5,%xmm5
        movaps %xmm5,nb314_tsc(%esp)
        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb314_half(%esp)
        movss nb314_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm3,%xmm4
        addps  %xmm4,%xmm4      ## 6.0
        movaps %xmm4,%xmm5
        addps  %xmm5,%xmm5      ## constant 12.0
        movaps %xmm1,nb314_half(%esp)
        movaps %xmm2,nb314_two(%esp)
        movaps %xmm3,nb314_three(%esp)
        movaps %xmm4,nb314_six(%esp)
        movaps %xmm5,nb314_twelve(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb314_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx),%ebx               ## ebx =ii 

        movl  nb314_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm5
        movss 12(%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movl nb314_p_facel(%ebp),%esi
        movss (%esi),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb314_qqMM(%esp)
        movaps %xmm4,nb314_qqMH(%esp)
        movaps %xmm5,nb314_qqHH(%esp)

        xorps %xmm0,%xmm0
        movl  nb314_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb314_p_ntype(%ebp),%edi
        imull (%edi),%ecx       ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb314_vdwparam(%ebp),%eax
        movlps (%eax,%edx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $0x55,%xmm1,%xmm1
        movaps %xmm0,nb314_c6(%esp)
        movaps %xmm1,nb314_c12(%esp)

_nb_kernel314_ia32_sse.nb314_threadloop: 
        movl  nb314_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel314_ia32_sse.nb314_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel314_ia32_sse.nb314_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb314_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb314_n(%esp)
        movl %ebx,nb314_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel314_ia32_sse.nb314_outerstart
        jmp _nb_kernel314_ia32_sse.nb314_end

_nb_kernel314_ia32_sse.nb314_outerstart: 
        ## ebx contains number of outer iterations
        addl nb314_nouter(%esp),%ebx
        movl %ebx,nb314_nouter(%esp)

_nb_kernel314_ia32_sse.nb314_outer: 
        movl  nb314_shift(%ebp),%eax            ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb314_is3(%esp)      ## store is3 

        movl  nb314_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb314_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb314_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb314_ii3(%esp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        addss (%eax,%ebx,4),%xmm3       ## ox
        addss 4(%eax,%ebx,4),%xmm4     ## oy
        addss 8(%eax,%ebx,4),%xmm5     ## oz
        addss 12(%eax,%ebx,4),%xmm6    ## h1x
        addss 16(%eax,%ebx,4),%xmm7    ## h1y
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm7,%xmm7
        movaps %xmm3,nb314_ixO(%esp)
        movaps %xmm4,nb314_iyO(%esp)
        movaps %xmm5,nb314_izO(%esp)
        movaps %xmm6,nb314_ixH1(%esp)
        movaps %xmm7,nb314_iyH1(%esp)

        movss %xmm2,%xmm6
        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 20(%eax,%ebx,4),%xmm6    ## h1z
        addss 24(%eax,%ebx,4),%xmm0    ## h2x
        addss 28(%eax,%ebx,4),%xmm1    ## h2y
        addss 32(%eax,%ebx,4),%xmm2    ## h2z
        addss 36(%eax,%ebx,4),%xmm3    ## mx
        addss 40(%eax,%ebx,4),%xmm4    ## my
        addss 44(%eax,%ebx,4),%xmm5    ## mz

        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm6,nb314_izH1(%esp)
        movaps %xmm0,nb314_ixH2(%esp)
        movaps %xmm1,nb314_iyH2(%esp)
        movaps %xmm2,nb314_izH2(%esp)
        movaps %xmm3,nb314_ixM(%esp)
        movaps %xmm4,nb314_iyM(%esp)
        movaps %xmm5,nb314_izM(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb314_vctot(%esp)
        movaps %xmm4,nb314_Vvdwtot(%esp)
        movaps %xmm4,nb314_fixO(%esp)
        movaps %xmm4,nb314_fiyO(%esp)
        movaps %xmm4,nb314_fizO(%esp)
        movaps %xmm4,nb314_fixH1(%esp)
        movaps %xmm4,nb314_fiyH1(%esp)
        movaps %xmm4,nb314_fizH1(%esp)
        movaps %xmm4,nb314_fixH2(%esp)
        movaps %xmm4,nb314_fiyH2(%esp)
        movaps %xmm4,nb314_fizH2(%esp)
        movaps %xmm4,nb314_fixM(%esp)
        movaps %xmm4,nb314_fiyM(%esp)
        movaps %xmm4,nb314_fizM(%esp)

        movl  nb314_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb314_pos(%ebp),%esi
        movl  nb314_faction(%ebp),%edi
        movl  nb314_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb314_innerjjnr(%esp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb314_ninner(%esp),%ecx
        movl  %ecx,nb314_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb314_innerk(%esp)   ## number of innerloop atoms 
        jge   _nb_kernel314_ia32_sse.nb314_unroll_loop
        jmp   _nb_kernel314_ia32_sse.nb314_single_check
_nb_kernel314_ia32_sse.nb314_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb314_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb314_innerjjnr(%esp)             ## advance pointer (unroll 4) 

        movl nb314_pos(%ebp),%esi       ## base of pos[] 

        leal  (%eax,%eax,2),%eax        ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx        ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move j coordinates to local temp variables
        ## Load Ox, Oy, Oz, H1x 
        movlps (%esi,%eax,4),%xmm1      ##  Oxa   Oya    -    -
        movlps (%esi,%ecx,4),%xmm4      ##  Oxc   Oyc    -    -
        movhps (%esi,%ebx,4),%xmm1      ##  Oxa   Oya   Oxb   Oyb 
        movhps (%esi,%edx,4),%xmm4      ##  Oxc   Oyc   Oxd   Oyd 
        movaps %xmm1,%xmm0              ##  Oxa   Oya   Oxb   Oyb 
        shufps $0x88,%xmm4,%xmm0       ##  Oxa   Oxb   Oxc   Oxd
        shufps $0xDD,%xmm4,%xmm1       ##  Oya   Oyb   Oyc   Oyd
        movlps 8(%esi,%eax,4),%xmm3     ##  Oza  H1xa    -    -
        movlps 8(%esi,%ecx,4),%xmm5     ##  Ozc  H1xc    -    -
        movhps 8(%esi,%ebx,4),%xmm3     ##  Oza  H1xa   Ozb  H1xb 
        movhps 8(%esi,%edx,4),%xmm5     ##  Ozc  H1xc   Ozd  H1xd 
        movaps %xmm3,%xmm2              ##  Oza  H1xa   Ozb  H1xb 
        shufps $0x88,%xmm5,%xmm2       ##  Oza   Ozb   Ozc   Ozd
        shufps $0xDD,%xmm5,%xmm3       ## H1xa  H1xb  H1xc  H1xd
        ## coordinates in xmm0-xmm3     
        ## store
        movaps %xmm0,nb314_jxO(%esp)
        movaps %xmm1,nb314_jyO(%esp)
        movaps %xmm2,nb314_jzO(%esp)
        movaps %xmm3,nb314_jxH1(%esp)

        ## Load H1y H1z H2x H2y 
        movlps 16(%esi,%eax,4),%xmm1
        movlps 16(%esi,%ecx,4),%xmm4
        movhps 16(%esi,%ebx,4),%xmm1
        movhps 16(%esi,%edx,4),%xmm4
        movaps %xmm1,%xmm0
        shufps $0x88,%xmm4,%xmm0
        shufps $0xDD,%xmm4,%xmm1
        movlps 24(%esi,%eax,4),%xmm3
        movlps 24(%esi,%ecx,4),%xmm5
        movhps 24(%esi,%ebx,4),%xmm3
        movhps 24(%esi,%edx,4),%xmm5
        movaps %xmm3,%xmm2
        shufps $0x88,%xmm5,%xmm2
        shufps $0xDD,%xmm5,%xmm3
        ## coordinates in xmm0-xmm3     
        ## store
        movaps %xmm0,nb314_jyH1(%esp)
        movaps %xmm1,nb314_jzH1(%esp)
        movaps %xmm2,nb314_jxH2(%esp)
        movaps %xmm3,nb314_jyH2(%esp)

        ## Load H2z Mx My Mz 
        movlps 32(%esi,%eax,4),%xmm1
        movlps 32(%esi,%ecx,4),%xmm4
        movhps 32(%esi,%ebx,4),%xmm1
        movhps 32(%esi,%edx,4),%xmm4
        movaps %xmm1,%xmm0
        shufps $0x88,%xmm4,%xmm0
        shufps $0xDD,%xmm4,%xmm1
        movlps 40(%esi,%eax,4),%xmm3
        movlps 40(%esi,%ecx,4),%xmm5
        movhps 40(%esi,%ebx,4),%xmm3
        movhps 40(%esi,%edx,4),%xmm5
        movaps %xmm3,%xmm2
        shufps $0x88,%xmm5,%xmm2
        shufps $0xDD,%xmm5,%xmm3
        ## coordinates in xmm0-xmm3     
        ## store
        movaps %xmm0,nb314_jzH2(%esp)
        movaps %xmm1,nb314_jxM(%esp)
        movaps %xmm2,nb314_jyM(%esp)
        movaps %xmm3,nb314_jzM(%esp)

        ## start calculating pairwise distances
        movaps nb314_ixO(%esp),%xmm0
        movaps nb314_iyO(%esp),%xmm1
        movaps nb314_izO(%esp),%xmm2
        movaps nb314_ixH1(%esp),%xmm3
        movaps nb314_iyH1(%esp),%xmm4
        movaps nb314_izH1(%esp),%xmm5
        subps  nb314_jxO(%esp),%xmm0
        subps  nb314_jyO(%esp),%xmm1
        subps  nb314_jzO(%esp),%xmm2
        subps  nb314_jxH1(%esp),%xmm3
        subps  nb314_jyH1(%esp),%xmm4
        subps  nb314_jzH1(%esp),%xmm5
        movaps %xmm0,nb314_dxOO(%esp)
        movaps %xmm1,nb314_dyOO(%esp)
        movaps %xmm2,nb314_dzOO(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb314_dxH1H1(%esp)
        movaps %xmm4,nb314_dyH1H1(%esp)
        movaps %xmm5,nb314_dzH1H1(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb314_rsqOO(%esp)
        movaps %xmm3,nb314_rsqH1H1(%esp)

        movaps nb314_ixH1(%esp),%xmm0
        movaps nb314_iyH1(%esp),%xmm1
        movaps nb314_izH1(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb314_jxH2(%esp),%xmm0
        subps  nb314_jyH2(%esp),%xmm1
        subps  nb314_jzH2(%esp),%xmm2
        subps  nb314_jxM(%esp),%xmm3
        subps  nb314_jyM(%esp),%xmm4
        subps  nb314_jzM(%esp),%xmm5
        movaps %xmm0,nb314_dxH1H2(%esp)
        movaps %xmm1,nb314_dyH1H2(%esp)
        movaps %xmm2,nb314_dzH1H2(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb314_dxH1M(%esp)
        movaps %xmm4,nb314_dyH1M(%esp)
        movaps %xmm5,nb314_dzH1M(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb314_rsqH1H2(%esp)
        movaps %xmm3,nb314_rsqH1M(%esp)

        movaps nb314_ixH2(%esp),%xmm0
        movaps nb314_iyH2(%esp),%xmm1
        movaps nb314_izH2(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb314_jxH1(%esp),%xmm0
        subps  nb314_jyH1(%esp),%xmm1
        subps  nb314_jzH1(%esp),%xmm2
        subps  nb314_jxH2(%esp),%xmm3
        subps  nb314_jyH2(%esp),%xmm4
        subps  nb314_jzH2(%esp),%xmm5
        movaps %xmm0,nb314_dxH2H1(%esp)
        movaps %xmm1,nb314_dyH2H1(%esp)
        movaps %xmm2,nb314_dzH2H1(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb314_dxH2H2(%esp)
        movaps %xmm4,nb314_dyH2H2(%esp)
        movaps %xmm5,nb314_dzH2H2(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb314_rsqH2H1(%esp)
        movaps %xmm3,nb314_rsqH2H2(%esp)

        movaps nb314_ixH2(%esp),%xmm0
        movaps nb314_iyH2(%esp),%xmm1
        movaps nb314_izH2(%esp),%xmm2
        movaps nb314_ixM(%esp),%xmm3
        movaps nb314_iyM(%esp),%xmm4
        movaps nb314_izM(%esp),%xmm5
        subps  nb314_jxM(%esp),%xmm0
        subps  nb314_jyM(%esp),%xmm1
        subps  nb314_jzM(%esp),%xmm2
        subps  nb314_jxH1(%esp),%xmm3
        subps  nb314_jyH1(%esp),%xmm4
        subps  nb314_jzH1(%esp),%xmm5
        movaps %xmm0,nb314_dxH2M(%esp)
        movaps %xmm1,nb314_dyH2M(%esp)
        movaps %xmm2,nb314_dzH2M(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb314_dxMH1(%esp)
        movaps %xmm4,nb314_dyMH1(%esp)
        movaps %xmm5,nb314_dzMH1(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,nb314_rsqH2M(%esp)
        movaps %xmm4,nb314_rsqMH1(%esp)

        movaps nb314_ixM(%esp),%xmm0
        movaps nb314_iyM(%esp),%xmm1
        movaps nb314_izM(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb314_jxH2(%esp),%xmm0
        subps  nb314_jyH2(%esp),%xmm1
        subps  nb314_jzH2(%esp),%xmm2
        subps  nb314_jxM(%esp),%xmm3
        subps  nb314_jyM(%esp),%xmm4
        subps  nb314_jzM(%esp),%xmm5
        movaps %xmm0,nb314_dxMH2(%esp)
        movaps %xmm1,nb314_dyMH2(%esp)
        movaps %xmm2,nb314_dzMH2(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb314_dxMM(%esp)
        movaps %xmm4,nb314_dyMM(%esp)
        movaps %xmm5,nb314_dzMM(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,nb314_rsqMH2(%esp)
        movaps %xmm4,nb314_rsqMM(%esp)

        ## start by doing reciprocal for OO
        movaps  nb314_rsqOO(%esp),%xmm7
        rcpps   %xmm7,%xmm2
        movaps  nb314_two(%esp),%xmm1
        mulps   %xmm2,%xmm7
        subps   %xmm7,%xmm1
        mulps   %xmm1,%xmm2 ## rinvsq 
        movaps %xmm2,nb314_rinvsqOO(%esp)

        ## next step is invsqrt - do two at a time.
        rsqrtps nb314_rsqH1H1(%esp),%xmm1
        rsqrtps nb314_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb314_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb314_rsqH1H1(%esp),%xmm1
        mulps   nb314_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb314_half(%esp),%xmm3   ## rinvH1H1 
        mulps   nb314_half(%esp),%xmm7   ## rinvH1H2 
        movaps  %xmm3,nb314_rinvH1H1(%esp)
        movaps  %xmm7,nb314_rinvH1H2(%esp)

        rsqrtps nb314_rsqH1M(%esp),%xmm1
        rsqrtps nb314_rsqH2H1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb314_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb314_rsqH1M(%esp),%xmm1
        mulps   nb314_rsqH2H1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb314_half(%esp),%xmm3
        mulps   nb314_half(%esp),%xmm7
        movaps  %xmm3,nb314_rinvH1M(%esp)
        movaps  %xmm7,nb314_rinvH2H1(%esp)

        rsqrtps nb314_rsqH2H2(%esp),%xmm1
        rsqrtps nb314_rsqH2M(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb314_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb314_rsqH2H2(%esp),%xmm1
        mulps   nb314_rsqH2M(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb314_half(%esp),%xmm3
        mulps   nb314_half(%esp),%xmm7
        movaps  %xmm3,nb314_rinvH2H2(%esp)
        movaps  %xmm7,nb314_rinvH2M(%esp)

        rsqrtps nb314_rsqMH1(%esp),%xmm1
        rsqrtps nb314_rsqMH2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb314_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb314_rsqMH1(%esp),%xmm1
        mulps   nb314_rsqMH2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb314_half(%esp),%xmm3
        mulps   nb314_half(%esp),%xmm7
        movaps  %xmm3,nb314_rinvMH1(%esp)
        movaps  %xmm7,nb314_rinvMH2(%esp)

        rsqrtps nb314_rsqMM(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb314_three(%esp),%xmm3
        mulps   nb314_rsqMM(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb314_half(%esp),%xmm3
        movaps  %xmm3,nb314_rinvMM(%esp)

        ## start with OO LJ interaction
        movaps nb314_rinvsqOO(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  %xmm1,%xmm1      ## rinv4
        mulps  %xmm0,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  nb314_c6(%esp),%xmm1
        mulps  nb314_c12(%esp),%xmm2
        movaps %xmm2,%xmm4
        subps  %xmm1,%xmm4
        addps  nb314_Vvdwtot(%esp),%xmm4
        mulps  nb314_six(%esp),%xmm1
        mulps  nb314_twelve(%esp),%xmm2
        movaps %xmm4,nb314_Vvdwtot(%esp)
        subps  %xmm1,%xmm2
        mulps  %xmm2,%xmm0      ## fscal 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb314_dxOO(%esp),%xmm0
        mulps nb314_dyOO(%esp),%xmm1
        mulps nb314_dzOO(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb314_fixO(%esp),%xmm0
        addps nb314_fiyO(%esp),%xmm1
        addps nb314_fizO(%esp),%xmm2
        movaps %xmm3,nb314_fjxO(%esp)
        movaps %xmm4,nb314_fjyO(%esp)
        movaps %xmm5,nb314_fjzO(%esp)
        movaps %xmm0,nb314_fixO(%esp)
        movaps %xmm1,nb314_fiyO(%esp)
        movaps %xmm2,nb314_fizO(%esp)

        ## Coulomb interactions - first H1H1
        movaps nb314_rinvH1H1(%esp),%xmm0

        movaps %xmm0,%xmm1
        mulps  nb314_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulps  nb314_tsc(%esp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movl nb314_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb314_two(%esp),%xmm7            ## two*Heps2 
        movaps nb314_qqHH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 
        ## update vctot 
        addps  nb314_vctot(%esp),%xmm5
        movaps %xmm5,nb314_vctot(%esp)

        xorps  %xmm1,%xmm1
        mulps  nb314_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5

        mulps nb314_dxH1H1(%esp),%xmm0
        mulps nb314_dyH1H1(%esp),%xmm1
        mulps nb314_dzH1H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb314_fixH1(%esp),%xmm0
        addps nb314_fiyH1(%esp),%xmm1
        addps nb314_fizH1(%esp),%xmm2
        movaps %xmm3,nb314_fjxH1(%esp)
        movaps %xmm4,nb314_fjyH1(%esp)
        movaps %xmm5,nb314_fjzH1(%esp)
        movaps %xmm0,nb314_fixH1(%esp)
        movaps %xmm1,nb314_fiyH1(%esp)
        movaps %xmm2,nb314_fizH1(%esp)

        ## H1-H2 interaction 
        movaps nb314_rinvH1H2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulps  nb314_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb314_two(%esp),%xmm7            ## two*Heps2 
        movaps nb314_qqHH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb314_vctot(%esp),%xmm5
        movaps %xmm5,nb314_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb314_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb314_dxH1H2(%esp),%xmm0
        mulps nb314_dyH1H2(%esp),%xmm1
        mulps nb314_dzH1H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb314_fixH1(%esp),%xmm0
        addps nb314_fiyH1(%esp),%xmm1
        addps nb314_fizH1(%esp),%xmm2
        movaps %xmm3,nb314_fjxH2(%esp)
        movaps %xmm4,nb314_fjyH2(%esp)
        movaps %xmm5,nb314_fjzH2(%esp)
        movaps %xmm0,nb314_fixH1(%esp)
        movaps %xmm1,nb314_fiyH1(%esp)
        movaps %xmm2,nb314_fizH1(%esp)

        ## H1-M interaction  
        movaps nb314_rinvH1M(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314_rsqH1M(%esp),%xmm1   ## xmm1=r 
        mulps  nb314_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb314_two(%esp),%xmm7            ## two*Heps2 
        movaps nb314_qqMH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb314_vctot(%esp),%xmm5
        movaps %xmm5,nb314_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb314_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb314_dxH1M(%esp),%xmm0
        mulps nb314_dyH1M(%esp),%xmm1
        mulps nb314_dzH1M(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb314_fixH1(%esp),%xmm0
        addps nb314_fiyH1(%esp),%xmm1
        addps nb314_fizH1(%esp),%xmm2
        movaps %xmm3,nb314_fjxM(%esp)
        movaps %xmm4,nb314_fjyM(%esp)
        movaps %xmm5,nb314_fjzM(%esp)
        movaps %xmm0,nb314_fixH1(%esp)
        movaps %xmm1,nb314_fiyH1(%esp)
        movaps %xmm2,nb314_fizH1(%esp)

        ## H2-H1 interaction 
        movaps nb314_rinvH2H1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulps  nb314_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb314_two(%esp),%xmm7            ## two*Heps2 
        movaps nb314_qqHH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb314_vctot(%esp),%xmm5
        movaps %xmm5,nb314_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb314_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        movaps nb314_fjxH1(%esp),%xmm3
        movaps nb314_fjyH1(%esp),%xmm4
        movaps nb314_fjzH1(%esp),%xmm5
        mulps nb314_dxH2H1(%esp),%xmm0
        mulps nb314_dyH2H1(%esp),%xmm1
        mulps nb314_dzH2H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb314_fixH2(%esp),%xmm0
        addps nb314_fiyH2(%esp),%xmm1
        addps nb314_fizH2(%esp),%xmm2
        movaps %xmm3,nb314_fjxH1(%esp)
        movaps %xmm4,nb314_fjyH1(%esp)
        movaps %xmm5,nb314_fjzH1(%esp)
        movaps %xmm0,nb314_fixH2(%esp)
        movaps %xmm1,nb314_fiyH2(%esp)
        movaps %xmm2,nb314_fizH2(%esp)

        ## H2-H2 interaction 
        movaps nb314_rinvH2H2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulps  nb314_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb314_two(%esp),%xmm7            ## two*Heps2 
        movaps nb314_qqHH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb314_vctot(%esp),%xmm5
        movaps %xmm5,nb314_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb314_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        movaps nb314_fjxH2(%esp),%xmm3
        movaps nb314_fjyH2(%esp),%xmm4
        movaps nb314_fjzH2(%esp),%xmm5
        mulps nb314_dxH2H2(%esp),%xmm0
        mulps nb314_dyH2H2(%esp),%xmm1
        mulps nb314_dzH2H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb314_fixH2(%esp),%xmm0
        addps nb314_fiyH2(%esp),%xmm1
        addps nb314_fizH2(%esp),%xmm2
        movaps %xmm3,nb314_fjxH2(%esp)
        movaps %xmm4,nb314_fjyH2(%esp)
        movaps %xmm5,nb314_fjzH2(%esp)
        movaps %xmm0,nb314_fixH2(%esp)
        movaps %xmm1,nb314_fiyH2(%esp)
        movaps %xmm2,nb314_fizH2(%esp)

        ## H2-M interaction 
        movaps nb314_rinvH2M(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314_rsqH2M(%esp),%xmm1   ## xmm1=r 
        mulps  nb314_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb314_two(%esp),%xmm7            ## two*Heps2 
        movaps nb314_qqMH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb314_vctot(%esp),%xmm5
        movaps %xmm5,nb314_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb314_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        movaps nb314_fjxM(%esp),%xmm3
        movaps nb314_fjyM(%esp),%xmm4
        movaps nb314_fjzM(%esp),%xmm5
        mulps nb314_dxH2M(%esp),%xmm0
        mulps nb314_dyH2M(%esp),%xmm1
        mulps nb314_dzH2M(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb314_fixH2(%esp),%xmm0
        addps nb314_fiyH2(%esp),%xmm1
        addps nb314_fizH2(%esp),%xmm2
        movaps %xmm3,nb314_fjxM(%esp)
        movaps %xmm4,nb314_fjyM(%esp)
        movaps %xmm5,nb314_fjzM(%esp)
        movaps %xmm0,nb314_fixH2(%esp)
        movaps %xmm1,nb314_fiyH2(%esp)
        movaps %xmm2,nb314_fizH2(%esp)

        ## M-H1 interaction 
        movaps nb314_rinvMH1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314_rsqMH1(%esp),%xmm1   ## xmm1=r 
        mulps  nb314_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb314_two(%esp),%xmm7            ## two*Heps2 
        movaps nb314_qqMH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb314_vctot(%esp),%xmm5
        movaps %xmm5,nb314_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb314_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        movaps nb314_fjxH1(%esp),%xmm3
        movaps nb314_fjyH1(%esp),%xmm4
        movaps nb314_fjzH1(%esp),%xmm5
        mulps nb314_dxMH1(%esp),%xmm0
        mulps nb314_dyMH1(%esp),%xmm1
        mulps nb314_dzMH1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb314_fixM(%esp),%xmm0
        addps nb314_fiyM(%esp),%xmm1
        addps nb314_fizM(%esp),%xmm2
        movaps %xmm3,nb314_fjxH1(%esp)
        movaps %xmm4,nb314_fjyH1(%esp)
        movaps %xmm5,nb314_fjzH1(%esp)
        movaps %xmm0,nb314_fixM(%esp)
        movaps %xmm1,nb314_fiyM(%esp)
        movaps %xmm2,nb314_fizM(%esp)

        ## M-H2 interaction 
        movaps nb314_rinvMH2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314_rsqMH2(%esp),%xmm1   ## xmm1=r 
        mulps  nb314_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb314_two(%esp),%xmm7            ## two*Heps2 
        movaps nb314_qqMH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb314_vctot(%esp),%xmm5
        movaps %xmm5,nb314_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb314_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        movaps nb314_fjxH2(%esp),%xmm3
        movaps nb314_fjyH2(%esp),%xmm4
        movaps nb314_fjzH2(%esp),%xmm5
        mulps nb314_dxMH2(%esp),%xmm0
        mulps nb314_dyMH2(%esp),%xmm1
        mulps nb314_dzMH2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb314_fixM(%esp),%xmm0
        addps nb314_fiyM(%esp),%xmm1
        addps nb314_fizM(%esp),%xmm2
        movaps %xmm3,nb314_fjxH2(%esp)
        movaps %xmm4,nb314_fjyH2(%esp)
        movaps %xmm5,nb314_fjzH2(%esp)
        movaps %xmm0,nb314_fixM(%esp)
        movaps %xmm1,nb314_fiyM(%esp)
        movaps %xmm2,nb314_fizM(%esp)

        ## M-M interaction 
        movaps nb314_rinvMM(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314_rsqMM(%esp),%xmm1   ## xmm1=r 
        mulps  nb314_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb314_two(%esp),%xmm7            ## two*Heps2 
        movaps nb314_qqMM(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb314_vctot(%esp),%xmm5
        movaps %xmm5,nb314_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb314_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        movaps nb314_fjxM(%esp),%xmm3
        movaps nb314_fjyM(%esp),%xmm4
        movaps nb314_fjzM(%esp),%xmm5
        mulps nb314_dxMM(%esp),%xmm0
        mulps nb314_dyMM(%esp),%xmm1
        mulps nb314_dzMM(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb314_fixM(%esp),%xmm0
        addps nb314_fiyM(%esp),%xmm1
        addps nb314_fizM(%esp),%xmm2
        movaps %xmm3,nb314_fjxM(%esp)
        movaps %xmm4,nb314_fjyM(%esp)
        movaps %xmm5,nb314_fjzM(%esp)
        movaps %xmm0,nb314_fixM(%esp)
        movaps %xmm1,nb314_fiyM(%esp)
        movaps %xmm2,nb314_fizM(%esp)

        movl nb314_faction(%ebp),%edi

        movd %mm0,%eax
        movd %mm1,%ebx
        movd %mm2,%ecx
        movd %mm3,%edx

        ## Did all interactions - now update j forces 
        ## 4 j waters with four atoms each.
        ## step 1 : transpose fjxO, fjyO, fjzO, fjxH1
        movaps nb314_fjxO(%esp),%xmm0
        movaps nb314_fjyO(%esp),%xmm1
        movaps nb314_fjzO(%esp),%xmm2
        movaps nb314_fjxH1(%esp),%xmm3
        movaps %xmm0,%xmm4
        movaps %xmm1,%xmm5
        unpcklps %xmm2,%xmm4
        unpcklps %xmm3,%xmm5
        unpckhps %xmm2,%xmm0
        unpckhps %xmm3,%xmm1
        movaps %xmm4,%xmm2
        movaps %xmm0,%xmm3
        unpcklps %xmm5,%xmm4
        unpckhps %xmm5,%xmm2

        unpcklps %xmm1,%xmm0
        unpckhps %xmm1,%xmm3
        ## results are now in xmm4, xmm2, xmm0, xmm3
        ## load the corresponding j forces from memory
        movlps   (%edi,%eax,4),%xmm1
        movlps   (%edi,%ebx,4),%xmm5
        movlps   (%edi,%ecx,4),%xmm6
        movlps   (%edi,%edx,4),%xmm7
        movhps   8(%edi,%eax,4),%xmm1
        movhps   8(%edi,%ebx,4),%xmm5
        movhps   8(%edi,%ecx,4),%xmm6
        movhps   8(%edi,%edx,4),%xmm7
        ## add
        addps    %xmm4,%xmm1
        addps    %xmm2,%xmm5
        addps    %xmm0,%xmm6
        addps    %xmm3,%xmm7
        ## store back
        movlps   %xmm1,(%edi,%eax,4)
        movlps   %xmm5,(%edi,%ebx,4)
        movlps   %xmm6,(%edi,%ecx,4)
        movlps   %xmm7,(%edi,%edx,4)
        movhps   %xmm1,8(%edi,%eax,4)
        movhps   %xmm5,8(%edi,%ebx,4)
        movhps   %xmm6,8(%edi,%ecx,4)
        movhps   %xmm7,8(%edi,%edx,4)

        ## step 2 : transpose fjyH1, fjzH1, fjxH2, fjyH2
        movaps nb314_fjyH1(%esp),%xmm0
        movaps nb314_fjzH1(%esp),%xmm1
        movaps nb314_fjxH2(%esp),%xmm2
        movaps nb314_fjyH2(%esp),%xmm3
        movaps %xmm0,%xmm4
        movaps %xmm1,%xmm5
        unpcklps %xmm2,%xmm4
        unpcklps %xmm3,%xmm5
        unpckhps %xmm2,%xmm0
        unpckhps %xmm3,%xmm1
        movaps %xmm4,%xmm2
        movaps %xmm0,%xmm3
        unpcklps %xmm5,%xmm4
        unpckhps %xmm5,%xmm2

        unpcklps %xmm1,%xmm0
        unpckhps %xmm1,%xmm3
        ## results are now in xmm4, xmm2, xmm0, xmm3
        ## load the corresponding j forces from memory
        movlps   16(%edi,%eax,4),%xmm1
        movlps   16(%edi,%ebx,4),%xmm5
        movlps   16(%edi,%ecx,4),%xmm6
        movlps   16(%edi,%edx,4),%xmm7
        movhps   24(%edi,%eax,4),%xmm1
        movhps   24(%edi,%ebx,4),%xmm5
        movhps   24(%edi,%ecx,4),%xmm6
        movhps   24(%edi,%edx,4),%xmm7
        ## add
        addps    %xmm4,%xmm1
        addps    %xmm2,%xmm5
        addps    %xmm0,%xmm6
        addps    %xmm3,%xmm7
        ## store back
        movlps   %xmm1,16(%edi,%eax,4)
        movlps   %xmm5,16(%edi,%ebx,4)
        movlps   %xmm6,16(%edi,%ecx,4)
        movlps   %xmm7,16(%edi,%edx,4)
        movhps   %xmm1,24(%edi,%eax,4)
        movhps   %xmm5,24(%edi,%ebx,4)
        movhps   %xmm6,24(%edi,%ecx,4)
        movhps   %xmm7,24(%edi,%edx,4)

        ## step 3 : transpose fjzH2, fjxM, fjyM, fjzM. xmm4 is scratch
        movaps nb314_fjzH2(%esp),%xmm0
        movaps nb314_fjxM(%esp),%xmm1
        movaps nb314_fjyM(%esp),%xmm2
        movaps nb314_fjzM(%esp),%xmm3

        movaps %xmm0,%xmm4
        movaps %xmm1,%xmm5
        unpcklps %xmm2,%xmm4
        unpcklps %xmm3,%xmm5
        unpckhps %xmm2,%xmm0
        unpckhps %xmm3,%xmm1
        movaps %xmm4,%xmm2
        movaps %xmm0,%xmm3
        unpcklps %xmm5,%xmm4
        unpckhps %xmm5,%xmm2

        unpcklps %xmm1,%xmm0
        unpckhps %xmm1,%xmm3

        ## results are now in xmm0, xmm1, xmm2, xmm3
        ## load the corresponding j forces from memory
        movlps   32(%edi,%eax,4),%xmm1
        movlps   32(%edi,%ebx,4),%xmm5
        movlps   32(%edi,%ecx,4),%xmm6
        movlps   32(%edi,%edx,4),%xmm7
        movhps   40(%edi,%eax,4),%xmm1
        movhps   40(%edi,%ebx,4),%xmm5
        movhps   40(%edi,%ecx,4),%xmm6
        movhps   40(%edi,%edx,4),%xmm7
        ## add
        addps    %xmm4,%xmm1
        addps    %xmm2,%xmm5
        addps    %xmm0,%xmm6
        addps    %xmm3,%xmm7
        ## store back
        movlps   %xmm1,32(%edi,%eax,4)
        movlps   %xmm5,32(%edi,%ebx,4)
        movlps   %xmm6,32(%edi,%ecx,4)
        movlps   %xmm7,32(%edi,%edx,4)
        movhps   %xmm1,40(%edi,%eax,4)
        movhps   %xmm5,40(%edi,%ebx,4)
        movhps   %xmm6,40(%edi,%ecx,4)
        movhps   %xmm7,40(%edi,%edx,4)

        ## should we do one more iteration? 
        subl $4,nb314_innerk(%esp)
        jl    _nb_kernel314_ia32_sse.nb314_single_check
        jmp   _nb_kernel314_ia32_sse.nb314_unroll_loop
_nb_kernel314_ia32_sse.nb314_single_check: 
        addl $4,nb314_innerk(%esp)
        jnz   _nb_kernel314_ia32_sse.nb314_single_loop
        jmp   _nb_kernel314_ia32_sse.nb314_updateouterdata
_nb_kernel314_ia32_sse.nb314_single_loop: 
        movl  nb314_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb314_innerjjnr(%esp)

        movl nb314_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates
        movlps (%esi,%eax,4),%xmm3              ##  Ox  Oy  
        movlps 16(%esi,%eax,4),%xmm4            ## H1y H1z 
        movlps 32(%esi,%eax,4),%xmm5            ## H2z  Mx 
        movhps 8(%esi,%eax,4),%xmm3             ##  Ox  Oy  Oz H1x
        movhps 24(%esi,%eax,4),%xmm4            ## H1y H1z H2x H2y
        movhps 40(%esi,%eax,4),%xmm5            ## H2z  Mx  My  Mz
        ## transpose
        movaps %xmm4,%xmm0
        movaps %xmm3,%xmm1
        movaps %xmm4,%xmm2
        movaps %xmm3,%xmm6
        shufps $18,%xmm5,%xmm4 ## (00010010)  h2x - Mx  - 
        shufps $193,%xmm0,%xmm3 ## (11000001)  Oy  - H1y - 
        shufps $35,%xmm5,%xmm2 ## (00100011) H2y - My  - 
        shufps $18,%xmm0,%xmm1 ## (00010010)  Oz  - H1z - 
        ##  xmm6: Ox - - H1x   xmm5: H2z - - Mz 
        shufps $140,%xmm4,%xmm6 ## (10001100) Ox H1x H2x Mx 
        shufps $136,%xmm2,%xmm3 ## (10001000) Oy H1y H2y My 
        shufps $200,%xmm5,%xmm1 ## (11001000) Oz H1z H2z Mz

        ## store all j coordinates in jO  
        movaps %xmm6,nb314_jxO(%esp)
        movaps %xmm3,nb314_jyO(%esp)
        movaps %xmm1,nb314_jzO(%esp)

        ## do O and H1 in parallel
        movaps nb314_ixO(%esp),%xmm0
        movaps nb314_iyO(%esp),%xmm1
        movaps nb314_izO(%esp),%xmm2
        movaps nb314_ixH1(%esp),%xmm3
        movaps nb314_iyH1(%esp),%xmm4
        movaps nb314_izH1(%esp),%xmm5
        subps  nb314_jxO(%esp),%xmm0
        subps  nb314_jyO(%esp),%xmm1
        subps  nb314_jzO(%esp),%xmm2
        subps  nb314_jxO(%esp),%xmm3
        subps  nb314_jyO(%esp),%xmm4
        subps  nb314_jzO(%esp),%xmm5

        movaps %xmm0,nb314_dxOO(%esp)
        movaps %xmm1,nb314_dyOO(%esp)
        movaps %xmm2,nb314_dzOO(%esp)
        movaps %xmm3,nb314_dxH1H1(%esp)
        movaps %xmm4,nb314_dyH1H1(%esp)
        movaps %xmm5,nb314_dzH1H1(%esp)

        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm3,%xmm4
        addps %xmm5,%xmm4       ## have rsq in xmm4
        ## Save H1 data in H1H1 
        movaps %xmm4,nb314_rsqH1H1(%esp)

        ## do 1/x for O and 1/sqrt(x) for H1
        rcpss  %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movss  nb314_two(%esp),%xmm2
        movaps  %xmm5,%xmm6
        mulss  %xmm1,%xmm0
        mulps   %xmm5,%xmm5
        subss  %xmm0,%xmm2
        movaps  nb314_three(%esp),%xmm7
        mulss  %xmm1,%xmm2      ## constant 1/r2


        mulps   %xmm4,%xmm5
        movss  %xmm2,%xmm0
        subps   %xmm5,%xmm7
        mulss  %xmm2,%xmm2
        mulps   %xmm6,%xmm7
        mulss  %xmm0,%xmm2      ## constant 1/r6
        mulps   nb314_half(%esp),%xmm7   ## rinv iH1 - j water 
        movss  %xmm2,%xmm1
        movaps %xmm7,nb314_rinvH1H1(%esp)

        mulss  %xmm2,%xmm2      ## constant 1/r12
        mulss  nb314_c6(%esp),%xmm1
        mulss  nb314_c12(%esp),%xmm2
        movss  %xmm2,%xmm3
        subss  %xmm1,%xmm3
        addss  nb314_Vvdwtot(%esp),%xmm3
        movss  %xmm3,nb314_Vvdwtot(%esp)
        mulss  nb314_six(%esp),%xmm1
        mulss  nb314_twelve(%esp),%xmm2
        subss  %xmm1,%xmm2
        mulss  %xmm2,%xmm0      ## fscal
        movss  %xmm0,%xmm1
        movss  %xmm0,%xmm2
        mulss  nb314_dxOO(%esp),%xmm0
        mulss  nb314_dyOO(%esp),%xmm1
        mulss  nb314_dzOO(%esp),%xmm2
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        subss   %xmm0,%xmm3
        subss   %xmm1,%xmm4
        subss   %xmm2,%xmm5
        movaps  %xmm3,nb314_fjxO(%esp)
        movaps  %xmm4,nb314_fjyO(%esp)
        movaps  %xmm5,nb314_fjzO(%esp)
        addss   nb314_fixO(%esp),%xmm0
        addss   nb314_fiyO(%esp),%xmm1
        addss   nb314_fizO(%esp),%xmm2
        movss  %xmm0,nb314_fixO(%esp)
        movss  %xmm1,nb314_fiyO(%esp)
        movss  %xmm2,nb314_fizO(%esp)

        ## do  H1 coulomb interaction
        movaps nb314_rinvH1H1(%esp),%xmm0   ## rinv 
        movaps %xmm0,%xmm1
        mulps  nb314_rsqH1H1(%esp),%xmm1        ## r
        mulps nb314_tsc(%esp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        movl nb314_VFtab(%ebp),%esi

        movlps (%esi,%ebx,4),%xmm4
        movlps (%esi,%ecx,4),%xmm3
        movlps (%esi,%edx,4),%xmm7
        movhps 8(%esi,%ebx,4),%xmm4
        movhps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%edx,4),%xmm7
        movaps %xmm3,%xmm6
        unpcklps %xmm7,%xmm6
        unpckhps %xmm7,%xmm3
        movaps %xmm4,%xmm5
        movaps %xmm4,%xmm7
        shufps $0x40,%xmm6,%xmm4
        shufps $0xE4,%xmm6,%xmm5
        movaps %xmm7,%xmm6
        shufps $0x48,%xmm3,%xmm6
        shufps $0xEC,%xmm3,%xmm7
        ## coulomb table ready, in xmm4-xmm7

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb314_two(%esp),%xmm7            ## two*Heps2 

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb314_qqHH(%esp),%xmm3
        movhps  nb314_qqMH(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## constant 11000001 

        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC 

        addps  nb314_vctot(%esp),%xmm5
        movaps %xmm5,nb314_vctot(%esp)

        mulps  nb314_tsc(%esp),%xmm3
        xorps  %xmm2,%xmm2
        subps  %xmm3,%xmm2
        mulps  %xmm2,%xmm0
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        mulps   nb314_dxH1H1(%esp),%xmm0
        mulps   nb314_dyH1H1(%esp),%xmm1
        mulps   nb314_dzH1H1(%esp),%xmm2
        ## update forces H1 - j water 
        movaps  nb314_fjxO(%esp),%xmm3
        movaps  nb314_fjyO(%esp),%xmm4
        movaps  nb314_fjzO(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb314_fjxO(%esp)
        movaps  %xmm4,nb314_fjyO(%esp)
        movaps  %xmm5,nb314_fjzO(%esp)
        addps   nb314_fixH1(%esp),%xmm0
        addps   nb314_fiyH1(%esp),%xmm1
        addps   nb314_fizH1(%esp),%xmm2
        movaps  %xmm0,nb314_fixH1(%esp)
        movaps  %xmm1,nb314_fiyH1(%esp)
        movaps  %xmm2,nb314_fizH1(%esp)

        ## i H2 & M simultaneously first get i particle coords: 
        movaps  nb314_ixH2(%esp),%xmm0
        movaps  nb314_iyH2(%esp),%xmm1
        movaps  nb314_izH2(%esp),%xmm2
        movaps  nb314_ixM(%esp),%xmm3
        movaps  nb314_iyM(%esp),%xmm4
        movaps  nb314_izM(%esp),%xmm5
        subps   nb314_jxO(%esp),%xmm0
        subps   nb314_jyO(%esp),%xmm1
        subps   nb314_jzO(%esp),%xmm2
        subps   nb314_jxO(%esp),%xmm3
        subps   nb314_jyO(%esp),%xmm4
        subps   nb314_jzO(%esp),%xmm5
        movaps %xmm0,nb314_dxH2H2(%esp)
        movaps %xmm1,nb314_dyH2H2(%esp)
        movaps %xmm2,nb314_dzH2H2(%esp)
        movaps %xmm3,nb314_dxMM(%esp)
        movaps %xmm4,nb314_dyMM(%esp)
        movaps %xmm5,nb314_dzMM(%esp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm1,%xmm0
        addps %xmm3,%xmm4
        addps %xmm2,%xmm0       ## have rsqH2 in xmm0 
        addps %xmm5,%xmm4       ## have rsqM in xmm4 

        ## start with H2, save data 
        movaps %xmm0,nb314_rsqH2H2(%esp)
        movaps %xmm4,nb314_rsqMM(%esp)
        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb314_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb314_half(%esp),%xmm3   ## rinv H2 - j water 
        mulps   nb314_half(%esp),%xmm7   ## rinv M - j water  

        movaps %xmm3,nb314_rinvH2H2(%esp)
        movaps %xmm7,nb314_rinvMM(%esp)

        movaps %xmm3,%xmm1
        mulps  nb314_rsqH2H2(%esp),%xmm1        ## xmm1=r 
        movaps %xmm3,%xmm0      ## xmm0=rinv 
        mulps  nb314_tsc(%esp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        movlps (%esi,%ebx,4),%xmm4
        movlps (%esi,%ecx,4),%xmm3
        movlps (%esi,%edx,4),%xmm7
        movhps 8(%esi,%ebx,4),%xmm4
        movhps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%edx,4),%xmm7
        movaps %xmm3,%xmm6
        unpcklps %xmm7,%xmm6
        unpckhps %xmm7,%xmm3
        movaps %xmm4,%xmm5
        movaps %xmm4,%xmm7
        shufps $0x40,%xmm6,%xmm4
        shufps $0xE4,%xmm6,%xmm5
        movaps %xmm7,%xmm6
        shufps $0x48,%xmm3,%xmm6
        shufps $0xEC,%xmm3,%xmm7
        ## coulomb table ready, in xmm4-xmm7
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb314_two(%esp),%xmm7            ## two*Heps2 

        xorps  %xmm3,%xmm3

        ## fetch charges to xmm3 (temporary) 
        movss   nb314_qqHH(%esp),%xmm3
        movhps  nb314_qqMH(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## constant 11000001

        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC 
        addps  nb314_vctot(%esp),%xmm5
        movaps %xmm5,nb314_vctot(%esp)

        xorps  %xmm1,%xmm1

        mulps nb314_tsc(%esp),%xmm3
        mulps %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps  %xmm1,%xmm0
        movaps  %xmm1,%xmm2
        mulps   nb314_dxH2H2(%esp),%xmm0
        mulps   nb314_dyH2H2(%esp),%xmm1
        mulps   nb314_dzH2H2(%esp),%xmm2
        ## update forces H1 - j water 
        movaps  nb314_fjxO(%esp),%xmm3
        movaps  nb314_fjyO(%esp),%xmm4
        movaps  nb314_fjzO(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb314_fjxO(%esp)
        movaps  %xmm4,nb314_fjyO(%esp)
        movaps  %xmm5,nb314_fjzO(%esp)
        addps   nb314_fixH2(%esp),%xmm0
        addps   nb314_fiyH2(%esp),%xmm1
        addps   nb314_fizH2(%esp),%xmm2
        movaps  %xmm0,nb314_fixH2(%esp)
        movaps  %xmm1,nb314_fiyH2(%esp)
        movaps  %xmm2,nb314_fizH2(%esp)

        ## do table for i M - j water interaction 
        movaps nb314_rinvMM(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314_rsqMM(%esp),%xmm1          ## xmm0=rinv, xmm1=r 
        mulps  nb314_tsc(%esp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        movlps (%esi,%ebx,4),%xmm4
        movlps (%esi,%ecx,4),%xmm3
        movlps (%esi,%edx,4),%xmm7
        movhps 8(%esi,%ebx,4),%xmm4
        movhps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%edx,4),%xmm7
        movaps %xmm3,%xmm6
        unpcklps %xmm7,%xmm6
        unpckhps %xmm7,%xmm3
        movaps %xmm4,%xmm5
        movaps %xmm4,%xmm7
        shufps $0x40,%xmm6,%xmm4
        shufps $0xE4,%xmm6,%xmm5
        movaps %xmm7,%xmm6
        shufps $0x48,%xmm3,%xmm6
        shufps $0xEC,%xmm3,%xmm7
        ## coulomb table ready, in xmm4-xmm7
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb314_two(%esp),%xmm7            ## two*Heps2 

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb314_qqMH(%esp),%xmm3
        movhps  nb314_qqMM(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## constant 11000001

        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC 
        addps  nb314_vctot(%esp),%xmm5
        movaps %xmm5,nb314_vctot(%esp)

        xorps  %xmm1,%xmm1

        mulps nb314_tsc(%esp),%xmm3
        mulps %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps  %xmm1,%xmm0
        movaps  %xmm1,%xmm2

        mulps   nb314_dxMM(%esp),%xmm0
        mulps   nb314_dyMM(%esp),%xmm1
        mulps   nb314_dzMM(%esp),%xmm2
        movaps  nb314_fjxO(%esp),%xmm3
        movaps  nb314_fjyO(%esp),%xmm4
        movaps  nb314_fjzO(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movl    nb314_faction(%ebp),%esi
        movaps  %xmm3,nb314_fjxO(%esp)
        movaps  %xmm4,nb314_fjyO(%esp)
        movaps  %xmm5,nb314_fjzO(%esp)
        addps   nb314_fixM(%esp),%xmm0
        addps   nb314_fiyM(%esp),%xmm1
        addps   nb314_fizM(%esp),%xmm2
        movaps  %xmm0,nb314_fixM(%esp)
        movaps  %xmm1,nb314_fiyM(%esp)
        movaps  %xmm2,nb314_fizM(%esp)

        ## update j water forces from local variables.
        ## transpose back first
        movaps  nb314_fjxO(%esp),%xmm0   ## Ox H1x H2x Mx 
        movaps  nb314_fjyO(%esp),%xmm1   ## Oy H1y H2y My
        movaps  nb314_fjzO(%esp),%xmm2   ## Oz H1z H2z Mz

        movaps  %xmm0,%xmm3
        movaps  %xmm0,%xmm4
        unpcklps %xmm1,%xmm3            ## Ox Oy - -
        shufps $0x1,%xmm2,%xmm4        ## h1x - Oz -
        movaps  %xmm1,%xmm5
        movaps  %xmm0,%xmm6
        unpcklps %xmm2,%xmm5            ## - - H1y H1z
        unpckhps %xmm1,%xmm6            ## h2x h2y - - 
        unpckhps %xmm2,%xmm1            ## - - My Mz

        shufps  $0x32,%xmm0,%xmm2 ## (00110010) h2z - Mx -
        shufps  $36,%xmm4,%xmm3 ## constant 00100100 ;# Ox Oy Oz H1x 
        shufps  $78,%xmm6,%xmm5 ## constant 01001110 ;# h1y h1z h2x h2y
        shufps  $232,%xmm1,%xmm2 ## constant 11101000 ;# h2z mx my mz

        movlps  (%esi,%eax,4),%xmm0
        movlps  16(%esi,%eax,4),%xmm1
        movlps  32(%esi,%eax,4),%xmm4
        movhps  8(%esi,%eax,4),%xmm0
        movhps  24(%esi,%eax,4),%xmm1
        movhps  40(%esi,%eax,4),%xmm4
        addps   %xmm3,%xmm0
        addps   %xmm5,%xmm1
        addps   %xmm2,%xmm4
        movlps   %xmm0,(%esi,%eax,4)
        movlps   %xmm1,16(%esi,%eax,4)
        movlps   %xmm4,32(%esi,%eax,4)
        movhps   %xmm0,8(%esi,%eax,4)
        movhps   %xmm1,24(%esi,%eax,4)
        movhps   %xmm4,40(%esi,%eax,4)

        decl nb314_innerk(%esp)
        jz    _nb_kernel314_ia32_sse.nb314_updateouterdata
        jmp   _nb_kernel314_ia32_sse.nb314_single_loop
_nb_kernel314_ia32_sse.nb314_updateouterdata: 
        movl  nb314_ii3(%esp),%ecx
        movl  nb314_faction(%ebp),%edi
        movl  nb314_fshift(%ebp),%esi
        movl  nb314_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb314_fixO(%esp),%xmm0
        movaps nb314_fiyO(%esp),%xmm1
        movaps nb314_fizO(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 in xmm0-xmm2 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        shufps $1,%xmm3,%xmm3
        shufps $1,%xmm4,%xmm4
        shufps $1,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2      ## xmm0-xmm2 has single force in pos0 

        ## increment i force 
        movss  (%edi,%ecx,4),%xmm3
        movss  4(%edi,%ecx,4),%xmm4
        movss  8(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,(%edi,%ecx,4)
        movss  %xmm4,4(%edi,%ecx,4)
        movss  %xmm5,8(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        movaps %xmm0,%xmm6
        movss %xmm2,%xmm7
        movlhps %xmm1,%xmm6
        shufps $8,%xmm6,%xmm6 ## constant 00001000      

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movaps nb314_fixH1(%esp),%xmm0
        movaps nb314_fiyH1(%esp),%xmm1
        movaps nb314_fizH1(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 in xmm0-xmm2 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        shufps $1,%xmm3,%xmm3
        shufps $1,%xmm4,%xmm4
        shufps $1,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2      ## xmm0-xmm2 has single force in pos0 

        ## increment i force 
        movss  12(%edi,%ecx,4),%xmm3
        movss  16(%edi,%ecx,4),%xmm4
        movss  20(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,12(%edi,%ecx,4)
        movss  %xmm4,16(%edi,%ecx,4)
        movss  %xmm5,20(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## constant 00001000      
        addps   %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movaps nb314_fixH2(%esp),%xmm0
        movaps nb314_fiyH2(%esp),%xmm1
        movaps nb314_fizH2(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 in xmm0-xmm2 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        shufps $1,%xmm3,%xmm3
        shufps $1,%xmm4,%xmm4
        shufps $1,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2      ## xmm0-xmm2 has single force in pos0 

        ## increment i force 
        movss  24(%edi,%ecx,4),%xmm3
        movss  28(%edi,%ecx,4),%xmm4
        movss  32(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,24(%edi,%ecx,4)
        movss  %xmm4,28(%edi,%ecx,4)
        movss  %xmm5,32(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## constant 00001000      
        addps   %xmm0,%xmm6

        ## accumulate Mi forces in xmm0, xmm1, xmm2 
        movaps nb314_fixM(%esp),%xmm0
        movaps nb314_fiyM(%esp),%xmm1
        movaps nb314_fizM(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 in xmm0-xmm2 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        shufps $1,%xmm3,%xmm3
        shufps $1,%xmm4,%xmm4
        shufps $1,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2      ## xmm0-xmm2 has single force in pos0 

        ## increment i force 
        movss  36(%edi,%ecx,4),%xmm3
        movss  40(%edi,%ecx,4),%xmm4
        movss  44(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,36(%edi,%ecx,4)
        movss  %xmm4,40(%edi,%ecx,4)
        movss  %xmm5,44(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## constant 00001000      
        addps   %xmm0,%xmm6

        ## increment fshift force  
        movlps  (%esi,%edx,4),%xmm3
        movss  8(%esi,%edx,4),%xmm4
        addps  %xmm6,%xmm3
        addss  %xmm7,%xmm4
        movlps  %xmm3,(%esi,%edx,4)
        movss  %xmm4,8(%esi,%edx,4)

        ## get n from stack
        movl nb314_n(%esp),%esi
        ## get group index for i particle 
        movl  nb314_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb314_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb314_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb314_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb314_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb314_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel314_ia32_sse.nb314_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb314_n(%esp)
        jmp _nb_kernel314_ia32_sse.nb314_outer
_nb_kernel314_ia32_sse.nb314_outerend: 
        ## check if more outer neighborlists remain
        movl  nb314_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel314_ia32_sse.nb314_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel314_ia32_sse.nb314_threadloop
_nb_kernel314_ia32_sse.nb314_end: 
        emms

        movl nb314_nouter(%esp),%eax
        movl nb314_ninner(%esp),%ebx
        movl nb314_outeriter(%ebp),%ecx
        movl nb314_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb314_salign(%esp),%eax
        addl %eax,%esp
        addl $1832,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret




.globl nb_kernel314nf_ia32_sse
.globl _nb_kernel314nf_ia32_sse
nb_kernel314nf_ia32_sse:        
_nb_kernel314nf_ia32_sse:       
.set nb314nf_p_nri, 8
.set nb314nf_iinr, 12
.set nb314nf_jindex, 16
.set nb314nf_jjnr, 20
.set nb314nf_shift, 24
.set nb314nf_shiftvec, 28
.set nb314nf_fshift, 32
.set nb314nf_gid, 36
.set nb314nf_pos, 40
.set nb314nf_faction, 44
.set nb314nf_charge, 48
.set nb314nf_p_facel, 52
.set nb314nf_argkrf, 56
.set nb314nf_argcrf, 60
.set nb314nf_Vc, 64
.set nb314nf_type, 68
.set nb314nf_p_ntype, 72
.set nb314nf_vdwparam, 76
.set nb314nf_Vvdw, 80
.set nb314nf_p_tabscale, 84
.set nb314nf_VFtab, 88
.set nb314nf_invsqrta, 92
.set nb314nf_dvda, 96
.set nb314nf_p_gbtabscale, 100
.set nb314nf_GBtab, 104
.set nb314nf_p_nthreads, 108
.set nb314nf_count, 112
.set nb314nf_mtx, 116
.set nb314nf_outeriter, 120
.set nb314nf_inneriter, 124
.set nb314nf_work, 128
## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb314nf_ixO, 0
.set nb314nf_iyO, 16
.set nb314nf_izO, 32
.set nb314nf_ixH1, 48
.set nb314nf_iyH1, 64
.set nb314nf_izH1, 80
.set nb314nf_ixH2, 96
.set nb314nf_iyH2, 112
.set nb314nf_izH2, 128
.set nb314nf_ixM, 144
.set nb314nf_iyM, 160
.set nb314nf_izM, 176
.set nb314nf_jxO, 192
.set nb314nf_jyO, 208
.set nb314nf_jzO, 224
.set nb314nf_jxH1, 240
.set nb314nf_jyH1, 256
.set nb314nf_jzH1, 272
.set nb314nf_jxH2, 288
.set nb314nf_jyH2, 304
.set nb314nf_jzH2, 320
.set nb314nf_jxM, 336
.set nb314nf_jyM, 352
.set nb314nf_jzM, 368
.set nb314nf_qqMM, 384
.set nb314nf_qqMH, 400
.set nb314nf_qqHH, 416
.set nb314nf_two, 432
.set nb314nf_tsc, 448
.set nb314nf_c6, 464
.set nb314nf_c12, 480
.set nb314nf_vctot, 496
.set nb314nf_Vvdwtot, 512
.set nb314nf_half, 528
.set nb314nf_three, 544
.set nb314nf_rsqOO, 560
.set nb314nf_rsqH1H1, 576
.set nb314nf_rsqH1H2, 592
.set nb314nf_rsqH1M, 608
.set nb314nf_rsqH2H1, 624
.set nb314nf_rsqH2H2, 640
.set nb314nf_rsqH2M, 656
.set nb314nf_rsqMH1, 672
.set nb314nf_rsqMH2, 688
.set nb314nf_rsqMM, 704
.set nb314nf_rinvsqOO, 720
.set nb314nf_rinvH1H1, 736
.set nb314nf_rinvH1H2, 752
.set nb314nf_rinvH1M, 768
.set nb314nf_rinvH2H1, 784
.set nb314nf_rinvH2H2, 800
.set nb314nf_rinvH2M, 816
.set nb314nf_rinvMH1, 832
.set nb314nf_rinvMH2, 848
.set nb314nf_rinvMM, 864
.set nb314nf_is3, 880
.set nb314nf_ii3, 884
.set nb314nf_innerjjnr, 888
.set nb314nf_innerk, 892
.set nb314nf_n, 896
.set nb314nf_nn1, 900
.set nb314nf_nri, 904
.set nb314nf_nouter, 908
.set nb314nf_ninner, 912
.set nb314nf_salign, 916
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $920,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb314nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb314nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb314nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb314nf_nouter(%esp)
        movl %eax,nb314nf_ninner(%esp)


        movl nb314nf_p_tabscale(%ebp),%eax
        movss (%eax),%xmm5
        shufps $0,%xmm5,%xmm5
        movaps %xmm5,nb314nf_tsc(%esp)
        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb314nf_half(%esp)
        movss nb314nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb314nf_half(%esp)
        movaps %xmm2,nb314nf_two(%esp)
        movaps %xmm3,nb314nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb314nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx),%ebx               ## ebx =ii 

        movl  nb314nf_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm5
        movss 12(%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movl nb314nf_p_facel(%ebp),%esi
        movss (%esi),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb314nf_qqMM(%esp)
        movaps %xmm4,nb314nf_qqMH(%esp)
        movaps %xmm5,nb314nf_qqHH(%esp)

        xorps %xmm0,%xmm0
        movl  nb314nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb314nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx       ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb314nf_vdwparam(%ebp),%eax
        movlps (%eax,%edx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $0x55,%xmm1,%xmm1
        movaps %xmm0,nb314nf_c6(%esp)
        movaps %xmm1,nb314nf_c12(%esp)

_nb_kernel314nf_ia32_sse.nb314nf_threadloop: 
        movl  nb314nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel314nf_ia32_sse.nb314nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel314nf_ia32_sse.nb314nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb314nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb314nf_n(%esp)
        movl %ebx,nb314nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel314nf_ia32_sse.nb314nf_outerstart
        jmp _nb_kernel314nf_ia32_sse.nb314nf_end

_nb_kernel314nf_ia32_sse.nb314nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb314nf_nouter(%esp),%ebx
        movl %ebx,nb314nf_nouter(%esp)

_nb_kernel314nf_ia32_sse.nb314nf_outer: 
        movl  nb314nf_shift(%ebp),%eax          ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb314nf_is3(%esp)            ## store is3 

        movl  nb314nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb314nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb314nf_pos(%ebp),%eax    ## eax = base of pos[]  
        movl  %ebx,nb314nf_ii3(%esp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        addss (%eax,%ebx,4),%xmm3       ## ox
        addss 4(%eax,%ebx,4),%xmm4     ## oy
        addss 8(%eax,%ebx,4),%xmm5     ## oz
        addss 12(%eax,%ebx,4),%xmm6    ## h1x
        addss 16(%eax,%ebx,4),%xmm7    ## h1y
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm7,%xmm7
        movaps %xmm3,nb314nf_ixO(%esp)
        movaps %xmm4,nb314nf_iyO(%esp)
        movaps %xmm5,nb314nf_izO(%esp)
        movaps %xmm6,nb314nf_ixH1(%esp)
        movaps %xmm7,nb314nf_iyH1(%esp)

        movss %xmm2,%xmm6
        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 20(%eax,%ebx,4),%xmm6    ## h1z
        addss 24(%eax,%ebx,4),%xmm0    ## h2x
        addss 28(%eax,%ebx,4),%xmm1    ## h2y
        addss 32(%eax,%ebx,4),%xmm2    ## h2z
        addss 36(%eax,%ebx,4),%xmm3    ## mx
        addss 40(%eax,%ebx,4),%xmm4    ## my
        addss 44(%eax,%ebx,4),%xmm5    ## mz

        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm6,nb314nf_izH1(%esp)
        movaps %xmm0,nb314nf_ixH2(%esp)
        movaps %xmm1,nb314nf_iyH2(%esp)
        movaps %xmm2,nb314nf_izH2(%esp)
        movaps %xmm3,nb314nf_ixM(%esp)
        movaps %xmm4,nb314nf_iyM(%esp)
        movaps %xmm5,nb314nf_izM(%esp)

        ## clear vctot  
        xorps %xmm4,%xmm4
        movaps %xmm4,nb314nf_vctot(%esp)
        movaps %xmm4,nb314nf_Vvdwtot(%esp)

        movl  nb314nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb314nf_pos(%ebp),%esi
        movl  nb314nf_faction(%ebp),%edi
        movl  nb314nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb314nf_innerjjnr(%esp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb314nf_ninner(%esp),%ecx
        movl  %ecx,nb314nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb314nf_innerk(%esp)         ## number of innerloop atoms 
        jge   _nb_kernel314nf_ia32_sse.nb314nf_unroll_loop
        jmp   _nb_kernel314nf_ia32_sse.nb314nf_single_check
_nb_kernel314nf_ia32_sse.nb314nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb314nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb314nf_innerjjnr(%esp)             ## advance pointer (unroll 4) 

        movl nb314nf_pos(%ebp),%esi     ## base of pos[] 

        leal  (%eax,%eax,2),%eax        ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx        ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move j coordinates to local temp variables
        ## Load Ox, Oy, Oz, H1x 
        movlps (%esi,%eax,4),%xmm1      ##  Oxa   Oya    -    -
        movlps (%esi,%ecx,4),%xmm4      ##  Oxc   Oyc    -    -
        movhps (%esi,%ebx,4),%xmm1      ##  Oxa   Oya   Oxb   Oyb 
        movhps (%esi,%edx,4),%xmm4      ##  Oxc   Oyc   Oxd   Oyd 
        movaps %xmm1,%xmm0              ##  Oxa   Oya   Oxb   Oyb 
        shufps $0x88,%xmm4,%xmm0       ##  Oxa   Oxb   Oxc   Oxd
        shufps $0xDD,%xmm4,%xmm1       ##  Oya   Oyb   Oyc   Oyd
        movlps 8(%esi,%eax,4),%xmm3     ##  Oza  H1xa    -    -
        movlps 8(%esi,%ecx,4),%xmm5     ##  Ozc  H1xc    -    -
        movhps 8(%esi,%ebx,4),%xmm3     ##  Oza  H1xa   Ozb  H1xb 
        movhps 8(%esi,%edx,4),%xmm5     ##  Ozc  H1xc   Ozd  H1xd 
        movaps %xmm3,%xmm2              ##  Oza  H1xa   Ozb  H1xb 
        shufps $0x88,%xmm5,%xmm2       ##  Oza   Ozb   Ozc   Ozd
        shufps $0xDD,%xmm5,%xmm3       ## H1xa  H1xb  H1xc  H1xd
        ## coordinates in xmm0-xmm3     
        ## store
        movaps %xmm0,nb314nf_jxO(%esp)
        movaps %xmm1,nb314nf_jyO(%esp)
        movaps %xmm2,nb314nf_jzO(%esp)
        movaps %xmm3,nb314nf_jxH1(%esp)

        ## Load H1y H1z H2x H2y 
        movlps 16(%esi,%eax,4),%xmm1
        movlps 16(%esi,%ecx,4),%xmm4
        movhps 16(%esi,%ebx,4),%xmm1
        movhps 16(%esi,%edx,4),%xmm4
        movaps %xmm1,%xmm0
        shufps $0x88,%xmm4,%xmm0
        shufps $0xDD,%xmm4,%xmm1
        movlps 24(%esi,%eax,4),%xmm3
        movlps 24(%esi,%ecx,4),%xmm5
        movhps 24(%esi,%ebx,4),%xmm3
        movhps 24(%esi,%edx,4),%xmm5
        movaps %xmm3,%xmm2
        shufps $0x88,%xmm5,%xmm2
        shufps $0xDD,%xmm5,%xmm3
        ## coordinates in xmm0-xmm3     
        ## store
        movaps %xmm0,nb314nf_jyH1(%esp)
        movaps %xmm1,nb314nf_jzH1(%esp)
        movaps %xmm2,nb314nf_jxH2(%esp)
        movaps %xmm3,nb314nf_jyH2(%esp)

        ## Load H2z Mx My Mz 
        movlps 32(%esi,%eax,4),%xmm1
        movlps 32(%esi,%ecx,4),%xmm4
        movhps 32(%esi,%ebx,4),%xmm1
        movhps 32(%esi,%edx,4),%xmm4
        movaps %xmm1,%xmm0
        shufps $0x88,%xmm4,%xmm0
        shufps $0xDD,%xmm4,%xmm1
        movlps 40(%esi,%eax,4),%xmm3
        movlps 40(%esi,%ecx,4),%xmm5
        movhps 40(%esi,%ebx,4),%xmm3
        movhps 40(%esi,%edx,4),%xmm5
        movaps %xmm3,%xmm2
        shufps $0x88,%xmm5,%xmm2
        shufps $0xDD,%xmm5,%xmm3
        ## coordinates in xmm0-xmm3     
        ## store
        movaps %xmm0,nb314nf_jzH2(%esp)
        movaps %xmm1,nb314nf_jxM(%esp)
        movaps %xmm2,nb314nf_jyM(%esp)
        movaps %xmm3,nb314nf_jzM(%esp)

        ## start calculating pairwise distances
        movaps nb314nf_ixO(%esp),%xmm0
        movaps nb314nf_iyO(%esp),%xmm1
        movaps nb314nf_izO(%esp),%xmm2
        movaps nb314nf_ixH1(%esp),%xmm3
        movaps nb314nf_iyH1(%esp),%xmm4
        movaps nb314nf_izH1(%esp),%xmm5
        subps  nb314nf_jxO(%esp),%xmm0
        subps  nb314nf_jyO(%esp),%xmm1
        subps  nb314nf_jzO(%esp),%xmm2
        subps  nb314nf_jxH1(%esp),%xmm3
        subps  nb314nf_jyH1(%esp),%xmm4
        subps  nb314nf_jzH1(%esp),%xmm5
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb314nf_rsqOO(%esp)
        movaps %xmm3,nb314nf_rsqH1H1(%esp)

        movaps nb314nf_ixH1(%esp),%xmm0
        movaps nb314nf_iyH1(%esp),%xmm1
        movaps nb314nf_izH1(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb314nf_jxH2(%esp),%xmm0
        subps  nb314nf_jyH2(%esp),%xmm1
        subps  nb314nf_jzH2(%esp),%xmm2
        subps  nb314nf_jxM(%esp),%xmm3
        subps  nb314nf_jyM(%esp),%xmm4
        subps  nb314nf_jzM(%esp),%xmm5
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb314nf_rsqH1H2(%esp)
        movaps %xmm3,nb314nf_rsqH1M(%esp)

        movaps nb314nf_ixH2(%esp),%xmm0
        movaps nb314nf_iyH2(%esp),%xmm1
        movaps nb314nf_izH2(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb314nf_jxH1(%esp),%xmm0
        subps  nb314nf_jyH1(%esp),%xmm1
        subps  nb314nf_jzH1(%esp),%xmm2
        subps  nb314nf_jxH2(%esp),%xmm3
        subps  nb314nf_jyH2(%esp),%xmm4
        subps  nb314nf_jzH2(%esp),%xmm5
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb314nf_rsqH2H1(%esp)
        movaps %xmm3,nb314nf_rsqH2H2(%esp)

        movaps nb314nf_ixH2(%esp),%xmm0
        movaps nb314nf_iyH2(%esp),%xmm1
        movaps nb314nf_izH2(%esp),%xmm2
        movaps nb314nf_ixM(%esp),%xmm3
        movaps nb314nf_iyM(%esp),%xmm4
        movaps nb314nf_izM(%esp),%xmm5
        subps  nb314nf_jxM(%esp),%xmm0
        subps  nb314nf_jyM(%esp),%xmm1
        subps  nb314nf_jzM(%esp),%xmm2
        subps  nb314nf_jxH1(%esp),%xmm3
        subps  nb314nf_jyH1(%esp),%xmm4
        subps  nb314nf_jzH1(%esp),%xmm5
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,nb314nf_rsqH2M(%esp)
        movaps %xmm4,nb314nf_rsqMH1(%esp)

        movaps nb314nf_ixM(%esp),%xmm0
        movaps nb314nf_iyM(%esp),%xmm1
        movaps nb314nf_izM(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb314nf_jxH2(%esp),%xmm0
        subps  nb314nf_jyH2(%esp),%xmm1
        subps  nb314nf_jzH2(%esp),%xmm2
        subps  nb314nf_jxM(%esp),%xmm3
        subps  nb314nf_jyM(%esp),%xmm4
        subps  nb314nf_jzM(%esp),%xmm5
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,nb314nf_rsqMH2(%esp)
        movaps %xmm4,nb314nf_rsqMM(%esp)

        ## start by doing reciprocal for OO
        movaps  nb314nf_rsqOO(%esp),%xmm7
        rcpps   %xmm7,%xmm2
        movaps  nb314nf_two(%esp),%xmm1
        mulps   %xmm2,%xmm7
        subps   %xmm7,%xmm1
        mulps   %xmm1,%xmm2 ## rinvsq 
        movaps %xmm2,nb314nf_rinvsqOO(%esp)

        ## next step is invsqrt - do two at a time.
        rsqrtps nb314nf_rsqH1H1(%esp),%xmm1
        rsqrtps nb314nf_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb314nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb314nf_rsqH1H1(%esp),%xmm1
        mulps   nb314nf_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb314nf_half(%esp),%xmm3   ## rinvH1H1 
        mulps   nb314nf_half(%esp),%xmm7   ## rinvH1H2 
        movaps  %xmm3,nb314nf_rinvH1H1(%esp)
        movaps  %xmm7,nb314nf_rinvH1H2(%esp)

        rsqrtps nb314nf_rsqH1M(%esp),%xmm1
        rsqrtps nb314nf_rsqH2H1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb314nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb314nf_rsqH1M(%esp),%xmm1
        mulps   nb314nf_rsqH2H1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb314nf_half(%esp),%xmm3
        mulps   nb314nf_half(%esp),%xmm7
        movaps  %xmm3,nb314nf_rinvH1M(%esp)
        movaps  %xmm7,nb314nf_rinvH2H1(%esp)

        rsqrtps nb314nf_rsqH2H2(%esp),%xmm1
        rsqrtps nb314nf_rsqH2M(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb314nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb314nf_rsqH2H2(%esp),%xmm1
        mulps   nb314nf_rsqH2M(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb314nf_half(%esp),%xmm3
        mulps   nb314nf_half(%esp),%xmm7
        movaps  %xmm3,nb314nf_rinvH2H2(%esp)
        movaps  %xmm7,nb314nf_rinvH2M(%esp)

        rsqrtps nb314nf_rsqMH1(%esp),%xmm1
        rsqrtps nb314nf_rsqMH2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb314nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb314nf_rsqMH1(%esp),%xmm1
        mulps   nb314nf_rsqMH2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb314nf_half(%esp),%xmm3
        mulps   nb314nf_half(%esp),%xmm7
        movaps  %xmm3,nb314nf_rinvMH1(%esp)
        movaps  %xmm7,nb314nf_rinvMH2(%esp)

        rsqrtps nb314nf_rsqMM(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb314nf_three(%esp),%xmm3
        mulps   nb314nf_rsqMM(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb314nf_half(%esp),%xmm3
        movaps  %xmm3,nb314nf_rinvMM(%esp)

        ## start with OO LJ interaction
        movaps nb314nf_rinvsqOO(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  %xmm1,%xmm1      ## rinv4
        mulps  %xmm0,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  nb314nf_c6(%esp),%xmm1
        mulps  nb314nf_c12(%esp),%xmm2
        movaps %xmm2,%xmm4
        subps  %xmm1,%xmm4
        addps  nb314nf_Vvdwtot(%esp),%xmm4
        movaps %xmm4,nb314nf_Vvdwtot(%esp)

        ## Coulomb interactions - first H1H1
        movaps nb314nf_rinvH1H1(%esp),%xmm0

        movaps %xmm0,%xmm1
        mulps  nb314nf_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulps  nb314nf_tsc(%esp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movl nb314nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb314nf_qqHH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## update vctot 
        addps  nb314nf_vctot(%esp),%xmm5
        movaps %xmm5,nb314nf_vctot(%esp)

        ## H1-H2 interaction 
        movaps nb314nf_rinvH1H2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314nf_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulps  nb314nf_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb314nf_qqHH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb314nf_vctot(%esp),%xmm5
        movaps %xmm5,nb314nf_vctot(%esp)

        ## H1-M interaction  
        movaps nb314nf_rinvH1M(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314nf_rsqH1M(%esp),%xmm1   ## xmm1=r 
        mulps  nb314nf_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb314nf_qqMH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb314nf_vctot(%esp),%xmm5
        movaps %xmm5,nb314nf_vctot(%esp)

        ## H2-H1 interaction 
        movaps nb314nf_rinvH2H1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314nf_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulps  nb314nf_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb314nf_qqHH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb314nf_vctot(%esp),%xmm5
        movaps %xmm5,nb314nf_vctot(%esp)

        ## H2-H2 interaction 
        movaps nb314nf_rinvH2H2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314nf_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulps  nb314nf_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb314nf_qqHH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb314nf_vctot(%esp),%xmm5
        movaps %xmm5,nb314nf_vctot(%esp)

        ## H2-M interaction 
        movaps nb314nf_rinvH2M(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314nf_rsqH2M(%esp),%xmm1   ## xmm1=r 
        mulps  nb314nf_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb314nf_qqMH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb314nf_vctot(%esp),%xmm5
        movaps %xmm5,nb314nf_vctot(%esp)

        ## M-H1 interaction 
        movaps nb314nf_rinvMH1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314nf_rsqMH1(%esp),%xmm1   ## xmm1=r 
        mulps  nb314nf_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb314nf_qqMH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb314nf_vctot(%esp),%xmm5
        movaps %xmm5,nb314nf_vctot(%esp)

        ## M-H2 interaction 
        movaps nb314nf_rinvMH2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314nf_rsqMH2(%esp),%xmm1   ## xmm1=r 
        mulps  nb314nf_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb314nf_qqMH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb314nf_vctot(%esp),%xmm5
        movaps %xmm5,nb314nf_vctot(%esp)

        ## M-M interaction 
        movaps nb314nf_rinvMM(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314nf_rsqMM(%esp),%xmm1   ## xmm1=r 
        mulps  nb314nf_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb314nf_qqMM(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb314nf_vctot(%esp),%xmm5
        movaps %xmm5,nb314nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb314nf_innerk(%esp)
        jl    _nb_kernel314nf_ia32_sse.nb314nf_single_check
        jmp   _nb_kernel314nf_ia32_sse.nb314nf_unroll_loop
_nb_kernel314nf_ia32_sse.nb314nf_single_check: 
        addl $4,nb314nf_innerk(%esp)
        jnz   _nb_kernel314nf_ia32_sse.nb314nf_single_loop
        jmp   _nb_kernel314nf_ia32_sse.nb314nf_updateouterdata
_nb_kernel314nf_ia32_sse.nb314nf_single_loop: 
        movl  nb314nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb314nf_innerjjnr(%esp)

        movl nb314nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates
        movlps (%esi,%eax,4),%xmm3              ##  Ox  Oy  
        movlps 16(%esi,%eax,4),%xmm4            ## H1y H1z 
        movlps 32(%esi,%eax,4),%xmm5            ## H2z  Mx 
        movhps 8(%esi,%eax,4),%xmm3             ##  Ox  Oy  Oz H1x
        movhps 24(%esi,%eax,4),%xmm4            ## H1y H1z H2x H2y
        movhps 40(%esi,%eax,4),%xmm5            ## H2z  Mx  My  Mz
        ## transpose
        movaps %xmm4,%xmm0
        movaps %xmm3,%xmm1
        movaps %xmm4,%xmm2
        movaps %xmm3,%xmm6
        shufps $18,%xmm5,%xmm4 ## (00010010)  h2x - Mx  - 
        shufps $193,%xmm0,%xmm3 ## (11000001)  Oy  - H1y - 
        shufps $35,%xmm5,%xmm2 ## (00100011) H2y - My  - 
        shufps $18,%xmm0,%xmm1 ## (00010010)  Oz  - H1z - 
        ##  xmm6: Ox - - H1x   xmm5: H2z - - Mz 
        shufps $140,%xmm4,%xmm6 ## (10001100) Ox H1x H2x Mx 
        shufps $136,%xmm2,%xmm3 ## (10001000) Oy H1y H2y My 
        shufps $200,%xmm5,%xmm1 ## (11001000) Oz H1z H2z Mz

        ## store all j coordinates in jO  
        movaps %xmm6,nb314nf_jxO(%esp)
        movaps %xmm3,nb314nf_jyO(%esp)
        movaps %xmm1,nb314nf_jzO(%esp)

        ## do O and H1 in parallel
        movaps nb314nf_ixO(%esp),%xmm0
        movaps nb314nf_iyO(%esp),%xmm1
        movaps nb314nf_izO(%esp),%xmm2
        movaps nb314nf_ixH1(%esp),%xmm3
        movaps nb314nf_iyH1(%esp),%xmm4
        movaps nb314nf_izH1(%esp),%xmm5
        subps  nb314nf_jxO(%esp),%xmm0
        subps  nb314nf_jyO(%esp),%xmm1
        subps  nb314nf_jzO(%esp),%xmm2
        subps  nb314nf_jxO(%esp),%xmm3
        subps  nb314nf_jyO(%esp),%xmm4
        subps  nb314nf_jzO(%esp),%xmm5

        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm3,%xmm4
        addps %xmm5,%xmm4       ## have rsq in xmm4
        ## Save H1 data in H1H1 
        movaps %xmm4,nb314nf_rsqH1H1(%esp)

        ## do 1/x for O and 1/sqrt(x) for H1
        rcpss  %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movss  nb314nf_two(%esp),%xmm2
        movaps  %xmm5,%xmm6
        mulss  %xmm1,%xmm0
        mulps   %xmm5,%xmm5
        subss  %xmm0,%xmm2
        movaps  nb314nf_three(%esp),%xmm7
        mulss  %xmm1,%xmm2      ## constant 1/r2


        mulps   %xmm4,%xmm5
        movss  %xmm2,%xmm0
        subps   %xmm5,%xmm7
        mulss  %xmm2,%xmm2
        mulps   %xmm6,%xmm7
        mulss  %xmm0,%xmm2      ## constant 1/r6
        mulps   nb314nf_half(%esp),%xmm7   ## rinv iH1 - j water 
        movss  %xmm2,%xmm1
        movaps %xmm7,nb314nf_rinvH1H1(%esp)

        mulss  %xmm2,%xmm2      ## constant 1/r12
        mulss  nb314nf_c6(%esp),%xmm1
        mulss  nb314nf_c12(%esp),%xmm2
        movss  %xmm2,%xmm3
        subss  %xmm1,%xmm3
        addss  nb314nf_Vvdwtot(%esp),%xmm3
        movss  %xmm3,nb314nf_Vvdwtot(%esp)

        ## do  H1 coulomb interaction
        movaps nb314nf_rinvH1H1(%esp),%xmm0   ## rinv 
        movaps %xmm0,%xmm1
        mulps  nb314nf_rsqH1H1(%esp),%xmm1      ## r
        mulps nb314nf_tsc(%esp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        movl nb314nf_VFtab(%ebp),%esi

        movlps (%esi,%ebx,4),%xmm4
        movlps (%esi,%ecx,4),%xmm3
        movlps (%esi,%edx,4),%xmm7
        movhps 8(%esi,%ebx,4),%xmm4
        movhps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%edx,4),%xmm7
        movaps %xmm3,%xmm6
        unpcklps %xmm7,%xmm6
        unpckhps %xmm7,%xmm3
        movaps %xmm4,%xmm5
        movaps %xmm4,%xmm7
        shufps $0x40,%xmm6,%xmm4
        shufps $0xE4,%xmm6,%xmm5
        movaps %xmm7,%xmm6
        shufps $0x48,%xmm3,%xmm6
        shufps $0xEC,%xmm3,%xmm7
        ## coulomb table ready, in xmm4-xmm7

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb314nf_qqHH(%esp),%xmm3
        movhps  nb314nf_qqMH(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## constant 11000001 

        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point xmm5 contains vcoul 

        addps  nb314nf_vctot(%esp),%xmm5
        movaps %xmm5,nb314nf_vctot(%esp)

        ## i H2 & M simultaneously first get i particle coords: 
        movaps  nb314nf_ixH2(%esp),%xmm0
        movaps  nb314nf_iyH2(%esp),%xmm1
        movaps  nb314nf_izH2(%esp),%xmm2
        movaps  nb314nf_ixM(%esp),%xmm3
        movaps  nb314nf_iyM(%esp),%xmm4
        movaps  nb314nf_izM(%esp),%xmm5
        subps   nb314nf_jxO(%esp),%xmm0
        subps   nb314nf_jyO(%esp),%xmm1
        subps   nb314nf_jzO(%esp),%xmm2
        subps   nb314nf_jxO(%esp),%xmm3
        subps   nb314nf_jyO(%esp),%xmm4
        subps   nb314nf_jzO(%esp),%xmm5
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm1,%xmm0
        addps %xmm3,%xmm4
        addps %xmm2,%xmm0       ## have rsqH2 in xmm0 
        addps %xmm5,%xmm4       ## have rsqM in xmm4 

        ## start with H2, save data 
        movaps %xmm0,nb314nf_rsqH2H2(%esp)
        movaps %xmm4,nb314nf_rsqMM(%esp)
        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb314nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb314nf_half(%esp),%xmm3   ## rinv H2 - j water 
        mulps   nb314nf_half(%esp),%xmm7   ## rinv M - j water  

        movaps %xmm3,nb314nf_rinvH2H2(%esp)
        movaps %xmm7,nb314nf_rinvMM(%esp)

        movaps %xmm3,%xmm1
        mulps  nb314nf_rsqH2H2(%esp),%xmm1      ## xmm1=r 
        movaps %xmm3,%xmm0      ## xmm0=rinv 
        mulps  nb314nf_tsc(%esp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        movlps (%esi,%ebx,4),%xmm4
        movlps (%esi,%ecx,4),%xmm3
        movlps (%esi,%edx,4),%xmm7
        movhps 8(%esi,%ebx,4),%xmm4
        movhps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%edx,4),%xmm7
        movaps %xmm3,%xmm6
        unpcklps %xmm7,%xmm6
        unpckhps %xmm7,%xmm3
        movaps %xmm4,%xmm5
        movaps %xmm4,%xmm7
        shufps $0x40,%xmm6,%xmm4
        shufps $0xE4,%xmm6,%xmm5
        movaps %xmm7,%xmm6
        shufps $0x48,%xmm3,%xmm6
        shufps $0xEC,%xmm3,%xmm7
        ## coulomb table ready, in xmm4-xmm7
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 

        xorps  %xmm3,%xmm3

        ## fetch charges to xmm3 (temporary) 
        movss   nb314nf_qqHH(%esp),%xmm3
        movhps  nb314nf_qqMH(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## constant 11000001

        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point xmm5 contains vcoul 
        addps  nb314nf_vctot(%esp),%xmm5
        movaps %xmm5,nb314nf_vctot(%esp)

        ## do table for i M - j water interaction 
        movaps nb314nf_rinvMM(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314nf_rsqMM(%esp),%xmm1        ## xmm0=rinv, xmm1=r 
        mulps  nb314nf_tsc(%esp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        movlps (%esi,%ebx,4),%xmm4
        movlps (%esi,%ecx,4),%xmm3
        movlps (%esi,%edx,4),%xmm7
        movhps 8(%esi,%ebx,4),%xmm4
        movhps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%edx,4),%xmm7
        movaps %xmm3,%xmm6
        unpcklps %xmm7,%xmm6
        unpckhps %xmm7,%xmm3
        movaps %xmm4,%xmm5
        movaps %xmm4,%xmm7
        shufps $0x40,%xmm6,%xmm4
        shufps $0xE4,%xmm6,%xmm5
        movaps %xmm7,%xmm6
        shufps $0x48,%xmm3,%xmm6
        shufps $0xEC,%xmm3,%xmm7
        ## coulomb table ready, in xmm4-xmm7
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb314nf_qqMH(%esp),%xmm3
        movhps  nb314nf_qqMM(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## constant 11000001

        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point xmm5 contains vcoul 
        addps  nb314nf_vctot(%esp),%xmm5
        movaps %xmm5,nb314nf_vctot(%esp)

        decl nb314nf_innerk(%esp)
        jz    _nb_kernel314nf_ia32_sse.nb314nf_updateouterdata
        jmp   _nb_kernel314nf_ia32_sse.nb314nf_single_loop
_nb_kernel314nf_ia32_sse.nb314nf_updateouterdata: 
        ## get n from stack
        movl nb314nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb314nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb314nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb314nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb314nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb314nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb314nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel314nf_ia32_sse.nb314nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb314nf_n(%esp)
        jmp _nb_kernel314nf_ia32_sse.nb314nf_outer
_nb_kernel314nf_ia32_sse.nb314nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb314nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel314nf_ia32_sse.nb314nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel314nf_ia32_sse.nb314nf_threadloop
_nb_kernel314nf_ia32_sse.nb314nf_end: 
        emms

        movl nb314nf_nouter(%esp),%eax
        movl nb314nf_ninner(%esp),%ebx
        movl nb314nf_outeriter(%ebp),%ecx
        movl nb314nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb314nf_salign(%esp),%eax
        addl %eax,%esp
        addl $920,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret

