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



.globl nb_kernel134_ia32_sse
.globl _nb_kernel134_ia32_sse
nb_kernel134_ia32_sse:  
_nb_kernel134_ia32_sse: 
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
        ## bottom of stack is cache-aligned for sse use 
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
.set nb134_c6, 928
.set nb134_c12, 944
.set nb134_tsc, 960
.set nb134_fscal, 976
.set nb134_vctot, 992
.set nb134_Vvdwtot, 1008
.set nb134_fixO, 1024
.set nb134_fiyO, 1040
.set nb134_fizO, 1056
.set nb134_fixH1, 1072
.set nb134_fiyH1, 1088
.set nb134_fizH1, 1104
.set nb134_fixH2, 1120
.set nb134_fiyH2, 1136
.set nb134_fizH2, 1152
.set nb134_fixM, 1168
.set nb134_fiyM, 1184
.set nb134_fizM, 1200
.set nb134_fjxO, 1216
.set nb134_fjyO, 1232
.set nb134_fjzO, 1248
.set nb134_fjxH1, 1264
.set nb134_fjyH1, 1280
.set nb134_fjzH1, 1296
.set nb134_fjxH2, 1312
.set nb134_fjyH2, 1328
.set nb134_fjzH2, 1344
.set nb134_fjxM, 1360
.set nb134_fjyM, 1376
.set nb134_fjzM, 1392
.set nb134_half, 1408
.set nb134_three, 1424
.set nb134_rsqOO, 1440
.set nb134_rsqH1H1, 1456
.set nb134_rsqH1H2, 1472
.set nb134_rsqH1M, 1488
.set nb134_rsqH2H1, 1504
.set nb134_rsqH2H2, 1520
.set nb134_rsqH2M, 1536
.set nb134_rsqMH1, 1552
.set nb134_rsqMH2, 1568
.set nb134_rsqMM, 1584
.set nb134_rinvOO, 1600
.set nb134_rinvH1H1, 1616
.set nb134_rinvH1H2, 1632
.set nb134_rinvH1M, 1648
.set nb134_rinvH2H1, 1664
.set nb134_rinvH2H2, 1680
.set nb134_rinvH2M, 1696
.set nb134_rinvMH1, 1712
.set nb134_rinvMH2, 1728
.set nb134_rinvMM, 1744
.set nb134_fstmp, 1760
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
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb134_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb134_half(%esp)
        movss nb134_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb134_half(%esp)
        movaps %xmm2,nb134_two(%esp)
        movaps %xmm3,nb134_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb134_iinr(%ebp),%ecx     ## ecx = pointer into iinr[]
        movl  (%ecx),%ebx               ## ebx =ii 

        movl  nb134_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm5
        movss 12(%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movl nb134_p_facel(%ebp),%esi
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
        movaps %xmm3,nb134_qqMM(%esp)
        movaps %xmm4,nb134_qqMH(%esp)
        movaps %xmm5,nb134_qqHH(%esp)

        xorps %xmm0,%xmm0
        movl  nb134_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb134_p_ntype(%ebp),%edi
        imull (%edi),%ecx ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb134_vdwparam(%ebp),%eax
        movlps (%eax,%edx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $0x55,%xmm1,%xmm1
        movaps %xmm0,nb134_c6(%esp)
        movaps %xmm1,nb134_c12(%esp)

_nb_kernel134_ia32_sse.nb134_threadloop: 
        movl  nb134_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel134_ia32_sse.nb134_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel134_ia32_sse.nb134_spinlock

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
        jg  _nb_kernel134_ia32_sse.nb134_outerstart
        jmp _nb_kernel134_ia32_sse.nb134_end

_nb_kernel134_ia32_sse.nb134_outerstart: 
        ## ebx contains number of outer iterations
        addl nb134_nouter(%esp),%ebx
        movl %ebx,nb134_nouter(%esp)

_nb_kernel134_ia32_sse.nb134_outer: 
        movl  nb134_shift(%ebp),%eax            ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb134_is3(%esp)      ## store is3 

        movl  nb134_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb134_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb134_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb134_ii3(%esp)

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
        movaps %xmm3,nb134_ixO(%esp)
        movaps %xmm4,nb134_iyO(%esp)
        movaps %xmm5,nb134_izO(%esp)
        movaps %xmm6,nb134_ixH1(%esp)
        movaps %xmm7,nb134_iyH1(%esp)

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
        movaps %xmm6,nb134_izH1(%esp)
        movaps %xmm0,nb134_ixH2(%esp)
        movaps %xmm1,nb134_iyH2(%esp)
        movaps %xmm2,nb134_izH2(%esp)
        movaps %xmm3,nb134_ixM(%esp)
        movaps %xmm4,nb134_iyM(%esp)
        movaps %xmm5,nb134_izM(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb134_vctot(%esp)
        movaps %xmm4,nb134_Vvdwtot(%esp)
        movaps %xmm4,nb134_fixO(%esp)
        movaps %xmm4,nb134_fiyO(%esp)
        movaps %xmm4,nb134_fizO(%esp)
        movaps %xmm4,nb134_fixH1(%esp)
        movaps %xmm4,nb134_fiyH1(%esp)
        movaps %xmm4,nb134_fizH1(%esp)
        movaps %xmm4,nb134_fixH2(%esp)
        movaps %xmm4,nb134_fiyH2(%esp)
        movaps %xmm4,nb134_fizH2(%esp)
        movaps %xmm4,nb134_fixM(%esp)
        movaps %xmm4,nb134_fiyM(%esp)
        movaps %xmm4,nb134_fizM(%esp)

        movl  nb134_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb134_pos(%ebp),%esi
        movl  nb134_faction(%ebp),%edi
        movl  nb134_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb134_innerjjnr(%esp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb134_ninner(%esp),%ecx
        movl  %ecx,nb134_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb134_innerk(%esp)   ## number of innerloop atoms 
        jge   _nb_kernel134_ia32_sse.nb134_unroll_loop
        jmp   _nb_kernel134_ia32_sse.nb134_single_check
_nb_kernel134_ia32_sse.nb134_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb134_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb134_innerjjnr(%esp)             ## advance pointer (unroll 4) 

        movl nb134_pos(%ebp),%esi       ## base of pos[] 

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
        movaps %xmm0,nb134_jxO(%esp)
        movaps %xmm1,nb134_jyO(%esp)
        movaps %xmm2,nb134_jzO(%esp)
        movaps %xmm3,nb134_jxH1(%esp)

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
        movaps %xmm0,nb134_jyH1(%esp)
        movaps %xmm1,nb134_jzH1(%esp)
        movaps %xmm2,nb134_jxH2(%esp)
        movaps %xmm3,nb134_jyH2(%esp)

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
        movaps %xmm0,nb134_jzH2(%esp)
        movaps %xmm1,nb134_jxM(%esp)
        movaps %xmm2,nb134_jyM(%esp)
        movaps %xmm3,nb134_jzM(%esp)

        ## start calculating pairwise distances
        movaps nb134_ixO(%esp),%xmm0
        movaps nb134_iyO(%esp),%xmm1
        movaps nb134_izO(%esp),%xmm2
        movaps nb134_ixH1(%esp),%xmm3
        movaps nb134_iyH1(%esp),%xmm4
        movaps nb134_izH1(%esp),%xmm5
        subps  nb134_jxO(%esp),%xmm0
        subps  nb134_jyO(%esp),%xmm1
        subps  nb134_jzO(%esp),%xmm2
        subps  nb134_jxH1(%esp),%xmm3
        subps  nb134_jyH1(%esp),%xmm4
        subps  nb134_jzH1(%esp),%xmm5
        movaps %xmm0,nb134_dxOO(%esp)
        movaps %xmm1,nb134_dyOO(%esp)
        movaps %xmm2,nb134_dzOO(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb134_dxH1H1(%esp)
        movaps %xmm4,nb134_dyH1H1(%esp)
        movaps %xmm5,nb134_dzH1H1(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb134_rsqOO(%esp)
        movaps %xmm3,nb134_rsqH1H1(%esp)

        movaps nb134_ixH1(%esp),%xmm0
        movaps nb134_iyH1(%esp),%xmm1
        movaps nb134_izH1(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb134_jxH2(%esp),%xmm0
        subps  nb134_jyH2(%esp),%xmm1
        subps  nb134_jzH2(%esp),%xmm2
        subps  nb134_jxM(%esp),%xmm3
        subps  nb134_jyM(%esp),%xmm4
        subps  nb134_jzM(%esp),%xmm5
        movaps %xmm0,nb134_dxH1H2(%esp)
        movaps %xmm1,nb134_dyH1H2(%esp)
        movaps %xmm2,nb134_dzH1H2(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb134_dxH1M(%esp)
        movaps %xmm4,nb134_dyH1M(%esp)
        movaps %xmm5,nb134_dzH1M(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb134_rsqH1H2(%esp)
        movaps %xmm3,nb134_rsqH1M(%esp)

        movaps nb134_ixH2(%esp),%xmm0
        movaps nb134_iyH2(%esp),%xmm1
        movaps nb134_izH2(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb134_jxH1(%esp),%xmm0
        subps  nb134_jyH1(%esp),%xmm1
        subps  nb134_jzH1(%esp),%xmm2
        subps  nb134_jxH2(%esp),%xmm3
        subps  nb134_jyH2(%esp),%xmm4
        subps  nb134_jzH2(%esp),%xmm5
        movaps %xmm0,nb134_dxH2H1(%esp)
        movaps %xmm1,nb134_dyH2H1(%esp)
        movaps %xmm2,nb134_dzH2H1(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb134_dxH2H2(%esp)
        movaps %xmm4,nb134_dyH2H2(%esp)
        movaps %xmm5,nb134_dzH2H2(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb134_rsqH2H1(%esp)
        movaps %xmm3,nb134_rsqH2H2(%esp)

        movaps nb134_ixH2(%esp),%xmm0
        movaps nb134_iyH2(%esp),%xmm1
        movaps nb134_izH2(%esp),%xmm2
        movaps nb134_ixM(%esp),%xmm3
        movaps nb134_iyM(%esp),%xmm4
        movaps nb134_izM(%esp),%xmm5
        subps  nb134_jxM(%esp),%xmm0
        subps  nb134_jyM(%esp),%xmm1
        subps  nb134_jzM(%esp),%xmm2
        subps  nb134_jxH1(%esp),%xmm3
        subps  nb134_jyH1(%esp),%xmm4
        subps  nb134_jzH1(%esp),%xmm5
        movaps %xmm0,nb134_dxH2M(%esp)
        movaps %xmm1,nb134_dyH2M(%esp)
        movaps %xmm2,nb134_dzH2M(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb134_dxMH1(%esp)
        movaps %xmm4,nb134_dyMH1(%esp)
        movaps %xmm5,nb134_dzMH1(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,nb134_rsqH2M(%esp)
        movaps %xmm4,nb134_rsqMH1(%esp)

        movaps nb134_ixM(%esp),%xmm0
        movaps nb134_iyM(%esp),%xmm1
        movaps nb134_izM(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb134_jxH2(%esp),%xmm0
        subps  nb134_jyH2(%esp),%xmm1
        subps  nb134_jzH2(%esp),%xmm2
        subps  nb134_jxM(%esp),%xmm3
        subps  nb134_jyM(%esp),%xmm4
        subps  nb134_jzM(%esp),%xmm5
        movaps %xmm0,nb134_dxMH2(%esp)
        movaps %xmm1,nb134_dyMH2(%esp)
        movaps %xmm2,nb134_dzMH2(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb134_dxMM(%esp)
        movaps %xmm4,nb134_dyMM(%esp)
        movaps %xmm5,nb134_dzMM(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,nb134_rsqMH2(%esp)
        movaps %xmm4,nb134_rsqMM(%esp)

        ## start by doing invsqrt for OO
        rsqrtps nb134_rsqOO(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb134_three(%esp),%xmm3
        mulps   nb134_rsqOO(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb134_half(%esp),%xmm3
        movaps  %xmm3,nb134_rinvOO(%esp)

        ## more invsqrt ops - do two at a time.
        rsqrtps nb134_rsqH1H1(%esp),%xmm1
        rsqrtps nb134_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb134_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb134_rsqH1H1(%esp),%xmm1
        mulps   nb134_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb134_half(%esp),%xmm3   ## rinvH1H1 
        mulps   nb134_half(%esp),%xmm7   ## rinvH1H2 
        movaps  %xmm3,nb134_rinvH1H1(%esp)
        movaps  %xmm7,nb134_rinvH1H2(%esp)

        rsqrtps nb134_rsqH1M(%esp),%xmm1
        rsqrtps nb134_rsqH2H1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb134_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb134_rsqH1M(%esp),%xmm1
        mulps   nb134_rsqH2H1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb134_half(%esp),%xmm3
        mulps   nb134_half(%esp),%xmm7
        movaps  %xmm3,nb134_rinvH1M(%esp)
        movaps  %xmm7,nb134_rinvH2H1(%esp)

        rsqrtps nb134_rsqH2H2(%esp),%xmm1
        rsqrtps nb134_rsqH2M(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb134_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb134_rsqH2H2(%esp),%xmm1
        mulps   nb134_rsqH2M(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb134_half(%esp),%xmm3
        mulps   nb134_half(%esp),%xmm7
        movaps  %xmm3,nb134_rinvH2H2(%esp)
        movaps  %xmm7,nb134_rinvH2M(%esp)

        rsqrtps nb134_rsqMH1(%esp),%xmm1
        rsqrtps nb134_rsqMH2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb134_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb134_rsqMH1(%esp),%xmm1
        mulps   nb134_rsqMH2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb134_half(%esp),%xmm3
        mulps   nb134_half(%esp),%xmm7
        movaps  %xmm3,nb134_rinvMH1(%esp)
        movaps  %xmm7,nb134_rinvMH2(%esp)

        rsqrtps nb134_rsqMM(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb134_three(%esp),%xmm3
        mulps   nb134_rsqMM(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb134_half(%esp),%xmm3
        movaps  %xmm3,nb134_rinvMM(%esp)

        ## start with OO LJ interaction
        movaps nb134_rinvOO(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb134_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulps  nb134_tsc(%esp),%xmm1

        movhlps %xmm1,%xmm2
    cvttps2pi %xmm1,%mm6
    cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
    movaps %xmm1,%xmm2
    mulps  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $3,%mm6
    pslld $3,%mm7

    movd %eax,%mm0
    movd %ebx,%mm1
    movd %ecx,%mm2
    movd %edx,%mm3

    movd %mm6,%eax
    psrlq $32,%mm6
    movd %mm7,%ecx
    psrlq $32,%mm7
    movd %mm6,%ebx
    movd %mm7,%edx

    movl nb134_VFtab(%ebp),%esi

    ## dispersion 
    movlps (%esi,%eax,4),%xmm5
    movlps (%esi,%ecx,4),%xmm7
    movhps (%esi,%ebx,4),%xmm5
    movhps (%esi,%edx,4),%xmm7 ## got half table 

    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 8(%esi,%eax,4),%xmm7
    movlps 8(%esi,%ecx,4),%xmm3
    movhps 8(%esi,%ebx,4),%xmm7
    movhps 8(%esi,%edx,4),%xmm3    ## other half of table  
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## dispersion table ready, in xmm4-xmm7 
    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp 
    mulps  nb134_two(%esp),%xmm7         ## two*Heps2 
    addps  %xmm6,%xmm7
    addps  %xmm5,%xmm7 ## xmm7=FF 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 

    movaps nb134_c6(%esp),%xmm4
    mulps  %xmm4,%xmm7   ## fijD 
    mulps  %xmm4,%xmm5   ## Vvdw6 
    movaps  %xmm7,nb134_fstmp(%esp)

    addps  nb134_Vvdwtot(%esp),%xmm5
    movaps %xmm5,nb134_Vvdwtot(%esp)

    ## repulsion 
    movlps 16(%esi,%eax,4),%xmm5
    movlps 16(%esi,%ecx,4),%xmm7
    movhps 16(%esi,%ebx,4),%xmm5
    movhps 16(%esi,%edx,4),%xmm7    ## got half table 

    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 24(%esi,%eax,4),%xmm7
    movlps 24(%esi,%ecx,4),%xmm3
    movhps 24(%esi,%ebx,4),%xmm7
    movhps 24(%esi,%edx,4),%xmm3    ## other half of table  
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## repulsion table ready, in xmm4-xmm7 
    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp 
    mulps  nb134_two(%esp),%xmm7         ## two*Heps2 
    addps  %xmm6,%xmm7
    addps  %xmm5,%xmm7 ## xmm7=FF 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 

    movaps nb134_c12(%esp),%xmm4
    mulps  %xmm4,%xmm7 ## fijR 
    mulps  %xmm4,%xmm5 ## Vvdw12 
    addps  nb134_fstmp(%esp),%xmm7
    mulps nb134_tsc(%esp),%xmm7

    addps  nb134_Vvdwtot(%esp),%xmm5
    movaps %xmm5,nb134_Vvdwtot(%esp)

        xorps %xmm2,%xmm2
        mulps  %xmm0,%xmm7
        subps %xmm7,%xmm2

        movaps %xmm2,%xmm0
        movaps %xmm2,%xmm1

    movd  %mm0,%eax
    movd  %mm1,%ebx
    movd  %mm2,%ecx
    movd  %mm3,%edx

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb134_dxOO(%esp),%xmm0
        mulps nb134_dyOO(%esp),%xmm1
        mulps nb134_dzOO(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb134_fixO(%esp),%xmm0
        addps nb134_fiyO(%esp),%xmm1
        addps nb134_fizO(%esp),%xmm2
        movaps %xmm3,nb134_fjxO(%esp)
        movaps %xmm4,nb134_fjyO(%esp)
        movaps %xmm5,nb134_fjzO(%esp)
        movaps %xmm0,nb134_fixO(%esp)
        movaps %xmm1,nb134_fiyO(%esp)
        movaps %xmm2,nb134_fizO(%esp)

        ## Coulomb interactions 
        ## start with H1-H1 interaction 
        movaps  nb134_rinvH1H1(%esp),%xmm7
        movaps  %xmm7,%xmm0
        mulps   %xmm0,%xmm0     ## xmm7=rinv, xmm0=rinvsq
        mulps  nb134_qqHH(%esp),%xmm7
        mulps  %xmm7,%xmm0      ## total fs H1 in xmm0 

        movaps %xmm7,%xmm6 ## local vctot summation variable

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb134_dxH1H1(%esp),%xmm0
        mulps nb134_dyH1H1(%esp),%xmm1
        mulps nb134_dzH1H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb134_fixH1(%esp),%xmm0
        addps nb134_fiyH1(%esp),%xmm1
        addps nb134_fizH1(%esp),%xmm2
        movaps %xmm3,nb134_fjxH1(%esp)
        movaps %xmm4,nb134_fjyH1(%esp)
        movaps %xmm5,nb134_fjzH1(%esp)
        movaps %xmm0,nb134_fixH1(%esp)
        movaps %xmm1,nb134_fiyH1(%esp)
        movaps %xmm2,nb134_fizH1(%esp)

        ## H1-H2 interaction 
        movaps  nb134_rinvH1H2(%esp),%xmm7
        movaps  %xmm7,%xmm0
        mulps   %xmm0,%xmm0     ## xmm7=rinv, xmm0=rinvsq
        mulps  nb134_qqHH(%esp),%xmm7
        mulps  %xmm7,%xmm0      ## total fs H1 in xmm0 

        addps  %xmm7,%xmm6 ## local vctot summation variable

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb134_dxH1H2(%esp),%xmm0
        mulps nb134_dyH1H2(%esp),%xmm1
        mulps nb134_dzH1H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb134_fixH1(%esp),%xmm0
        addps nb134_fiyH1(%esp),%xmm1
        addps nb134_fizH1(%esp),%xmm2
        movaps %xmm3,nb134_fjxH2(%esp)
        movaps %xmm4,nb134_fjyH2(%esp)
        movaps %xmm5,nb134_fjzH2(%esp)
        movaps %xmm0,nb134_fixH1(%esp)
        movaps %xmm1,nb134_fiyH1(%esp)
        movaps %xmm2,nb134_fizH1(%esp)

        ## H1-M interaction  
        movaps  nb134_rinvH1M(%esp),%xmm7
        movaps  %xmm7,%xmm0
        mulps   %xmm0,%xmm0     ## xmm7=rinv, xmm0=rinvsq
        mulps  nb134_qqMH(%esp),%xmm7
        mulps  %xmm7,%xmm0      ## total fs H1 in xmm0 

        addps  %xmm7,%xmm6 ## local vctot summation variable

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb134_dxH1M(%esp),%xmm0
        mulps nb134_dyH1M(%esp),%xmm1
        mulps nb134_dzH1M(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb134_fixH1(%esp),%xmm0
        addps nb134_fiyH1(%esp),%xmm1
        addps nb134_fizH1(%esp),%xmm2
        movaps %xmm3,nb134_fjxM(%esp)
        movaps %xmm4,nb134_fjyM(%esp)
        movaps %xmm5,nb134_fjzM(%esp)
        movaps %xmm0,nb134_fixH1(%esp)
        movaps %xmm1,nb134_fiyH1(%esp)
        movaps %xmm2,nb134_fizH1(%esp)

        ## H2-H1 interaction 
        movaps  nb134_rinvH2H1(%esp),%xmm7
        movaps  %xmm7,%xmm0
        mulps   %xmm0,%xmm0     ## xmm7=rinv, xmm0=rinvsq
        mulps  nb134_qqHH(%esp),%xmm7
        mulps  %xmm7,%xmm0      ## total fs H1 in xmm0 

        addps  %xmm7,%xmm6 ## local vctot summation variable

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb134_fjxH1(%esp),%xmm3
        movaps nb134_fjyH1(%esp),%xmm4
        movaps nb134_fjzH1(%esp),%xmm5
        mulps nb134_dxH2H1(%esp),%xmm0
        mulps nb134_dyH2H1(%esp),%xmm1
        mulps nb134_dzH2H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb134_fixH2(%esp),%xmm0
        addps nb134_fiyH2(%esp),%xmm1
        addps nb134_fizH2(%esp),%xmm2
        movaps %xmm3,nb134_fjxH1(%esp)
        movaps %xmm4,nb134_fjyH1(%esp)
        movaps %xmm5,nb134_fjzH1(%esp)
        movaps %xmm0,nb134_fixH2(%esp)
        movaps %xmm1,nb134_fiyH2(%esp)
        movaps %xmm2,nb134_fizH2(%esp)

        ## H2-H2 interaction 
        movaps  nb134_rinvH2H2(%esp),%xmm7
        movaps  %xmm7,%xmm0
        mulps   %xmm0,%xmm0     ## xmm7=rinv, xmm0=rinvsq
        mulps  nb134_qqHH(%esp),%xmm7
        mulps  %xmm7,%xmm0      ## total fs H1 in xmm0 

        addps  %xmm7,%xmm6 ## local vctot summation variable

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb134_fjxH2(%esp),%xmm3
        movaps nb134_fjyH2(%esp),%xmm4
        movaps nb134_fjzH2(%esp),%xmm5
        mulps nb134_dxH2H2(%esp),%xmm0
        mulps nb134_dyH2H2(%esp),%xmm1
        mulps nb134_dzH2H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb134_fixH2(%esp),%xmm0
        addps nb134_fiyH2(%esp),%xmm1
        addps nb134_fizH2(%esp),%xmm2
        movaps %xmm3,nb134_fjxH2(%esp)
        movaps %xmm4,nb134_fjyH2(%esp)
        movaps %xmm5,nb134_fjzH2(%esp)
        movaps %xmm0,nb134_fixH2(%esp)
        movaps %xmm1,nb134_fiyH2(%esp)
        movaps %xmm2,nb134_fizH2(%esp)

        ## H2-M interaction 
        movaps  nb134_rinvH2M(%esp),%xmm7
        movaps  %xmm7,%xmm0
        mulps   %xmm0,%xmm0     ## xmm7=rinv, xmm0=rinvsq
        mulps  nb134_qqMH(%esp),%xmm7
        mulps  %xmm7,%xmm0      ## total fs H1 in xmm0 

        addps  %xmm7,%xmm6 ## local vctot summation variable

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb134_fjxM(%esp),%xmm3
        movaps nb134_fjyM(%esp),%xmm4
        movaps nb134_fjzM(%esp),%xmm5
        mulps nb134_dxH2M(%esp),%xmm0
        mulps nb134_dyH2M(%esp),%xmm1
        mulps nb134_dzH2M(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb134_fixH2(%esp),%xmm0
        addps nb134_fiyH2(%esp),%xmm1
        addps nb134_fizH2(%esp),%xmm2
        movaps %xmm3,nb134_fjxM(%esp)
        movaps %xmm4,nb134_fjyM(%esp)
        movaps %xmm5,nb134_fjzM(%esp)
        movaps %xmm0,nb134_fixH2(%esp)
        movaps %xmm1,nb134_fiyH2(%esp)
        movaps %xmm2,nb134_fizH2(%esp)

        ## M-H1 interaction 
        movaps  nb134_rinvMH1(%esp),%xmm7
        movaps  %xmm7,%xmm0
        mulps   %xmm0,%xmm0     ## xmm7=rinv, xmm0=rinvsq
        mulps  nb134_qqMH(%esp),%xmm7
        mulps  %xmm7,%xmm0      ## total fs in xmm0 

        addps  %xmm7,%xmm6 ## local vctot summation variable

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb134_fjxH1(%esp),%xmm3
        movaps nb134_fjyH1(%esp),%xmm4
        movaps nb134_fjzH1(%esp),%xmm5
        mulps nb134_dxMH1(%esp),%xmm0
        mulps nb134_dyMH1(%esp),%xmm1
        mulps nb134_dzMH1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb134_fixM(%esp),%xmm0
        addps nb134_fiyM(%esp),%xmm1
        addps nb134_fizM(%esp),%xmm2
        movaps %xmm3,nb134_fjxH1(%esp)
        movaps %xmm4,nb134_fjyH1(%esp)
        movaps %xmm5,nb134_fjzH1(%esp)
        movaps %xmm0,nb134_fixM(%esp)
        movaps %xmm1,nb134_fiyM(%esp)
        movaps %xmm2,nb134_fizM(%esp)

        ## M-H2 interaction 
        movaps  nb134_rinvMH2(%esp),%xmm7
        movaps  %xmm7,%xmm0
        mulps   %xmm0,%xmm0     ## xmm7=rinv, xmm0=rinvsq
        mulps  nb134_qqMH(%esp),%xmm7
        mulps  %xmm7,%xmm0      ## total fs in xmm0 

        addps  %xmm7,%xmm6 ## local vctot summation variable

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb134_fjxH2(%esp),%xmm3
        movaps nb134_fjyH2(%esp),%xmm4
        movaps nb134_fjzH2(%esp),%xmm5
        mulps nb134_dxMH2(%esp),%xmm0
        mulps nb134_dyMH2(%esp),%xmm1
        mulps nb134_dzMH2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb134_fixM(%esp),%xmm0
        addps nb134_fiyM(%esp),%xmm1
        addps nb134_fizM(%esp),%xmm2
        movaps %xmm3,nb134_fjxH2(%esp)
        movaps %xmm4,nb134_fjyH2(%esp)
        movaps %xmm5,nb134_fjzH2(%esp)
        movaps %xmm0,nb134_fixM(%esp)
        movaps %xmm1,nb134_fiyM(%esp)
        movaps %xmm2,nb134_fizM(%esp)

        ## M-M interaction 
        movaps  nb134_rinvMM(%esp),%xmm7
        movaps  %xmm7,%xmm0
        mulps   %xmm0,%xmm0     ## xmm7=rinv, xmm0=rinvsq
        mulps  nb134_qqMM(%esp),%xmm7
        mulps  %xmm7,%xmm0      ## total fs in xmm0 

        addps  %xmm7,%xmm6 ## local vctot summation variable

        movaps %xmm0,%xmm1
        addps  nb134_vctot(%esp),%xmm6
        movaps %xmm6,nb134_vctot(%esp)
        movaps %xmm0,%xmm2

        movaps nb134_fjxM(%esp),%xmm3
        movaps nb134_fjyM(%esp),%xmm4
        movaps nb134_fjzM(%esp),%xmm5
        mulps nb134_dxMM(%esp),%xmm0
        mulps nb134_dyMM(%esp),%xmm1
        mulps nb134_dzMM(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb134_fixM(%esp),%xmm0
        addps nb134_fiyM(%esp),%xmm1
        addps nb134_fizM(%esp),%xmm2
        movaps %xmm3,nb134_fjxM(%esp)
        movaps %xmm4,nb134_fjyM(%esp)
        movaps %xmm5,nb134_fjzM(%esp)
        movaps %xmm0,nb134_fixM(%esp)
        movaps %xmm1,nb134_fiyM(%esp)
        movaps %xmm2,nb134_fizM(%esp)

        movl nb134_faction(%ebp),%edi
        ## update j forces 
        ## 4 j waters with four atoms each.
        ## step 1 : transpose fjxO, fjyO, fjzO, fjxH1
        movaps nb134_fjxO(%esp),%xmm0
        movaps nb134_fjyO(%esp),%xmm1
        movaps nb134_fjzO(%esp),%xmm2
        movaps nb134_fjxH1(%esp),%xmm3
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
        movaps nb134_fjyH1(%esp),%xmm0
        movaps nb134_fjzH1(%esp),%xmm1
        movaps nb134_fjxH2(%esp),%xmm2
        movaps nb134_fjyH2(%esp),%xmm3
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
        movaps nb134_fjzH2(%esp),%xmm0
        movaps nb134_fjxM(%esp),%xmm1
        movaps nb134_fjyM(%esp),%xmm2
        movaps nb134_fjzM(%esp),%xmm3

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
        subl $4,nb134_innerk(%esp)
        jl    _nb_kernel134_ia32_sse.nb134_single_check
        jmp   _nb_kernel134_ia32_sse.nb134_unroll_loop
_nb_kernel134_ia32_sse.nb134_single_check: 
        addl $4,nb134_innerk(%esp)
        jnz   _nb_kernel134_ia32_sse.nb134_single_loop
        jmp   _nb_kernel134_ia32_sse.nb134_updateouterdata
_nb_kernel134_ia32_sse.nb134_single_loop: 
        movl  nb134_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb134_innerjjnr(%esp)

        movl nb134_pos(%ebp),%esi
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
        movaps %xmm6,nb134_jxO(%esp)
        movaps %xmm3,nb134_jyO(%esp)
        movaps %xmm1,nb134_jzO(%esp)

        ## do O and M in parallel
        movaps nb134_ixO(%esp),%xmm0
        movaps nb134_iyO(%esp),%xmm1
        movaps nb134_izO(%esp),%xmm2
        movaps nb134_ixM(%esp),%xmm3
        movaps nb134_iyM(%esp),%xmm4
        movaps nb134_izM(%esp),%xmm5
        subps  nb134_jxO(%esp),%xmm0
        subps  nb134_jyO(%esp),%xmm1
        subps  nb134_jzO(%esp),%xmm2
        subps  nb134_jxO(%esp),%xmm3
        subps  nb134_jyO(%esp),%xmm4
        subps  nb134_jzO(%esp),%xmm5

        movaps %xmm0,nb134_dxOO(%esp)
        movaps %xmm1,nb134_dyOO(%esp)
        movaps %xmm2,nb134_dzOO(%esp)
        movaps %xmm3,nb134_dxMM(%esp)
        movaps %xmm4,nb134_dyMM(%esp)
        movaps %xmm5,nb134_dzMM(%esp)

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
        ## Save data 
        movaps %xmm0,nb134_rsqOO(%esp)
        movaps %xmm4,nb134_rsqMM(%esp)

        ## do 1/x for O
        rsqrtps %xmm0,%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb134_three(%esp),%xmm3
        mulps   %xmm0,%xmm1     ## rsq*lu*lu 
        subps   %xmm1,%xmm3     ## constant 30-rsq*lu*lu 
        mulps   %xmm2,%xmm3     ## lu*(3-rsq*lu*lu) 
        mulps   nb134_half(%esp),%xmm3
        movaps  %xmm3,nb134_rinvOO(%esp)        ## rinvH2 

        ## 1/sqrt(x) for M
        rsqrtps %xmm4,%xmm5
        movaps  %xmm5,%xmm6
        mulps   %xmm5,%xmm5
        movaps  nb134_three(%esp),%xmm7
        mulps   %xmm4,%xmm5
        subps   %xmm5,%xmm7
        mulps   %xmm6,%xmm7
        mulps   nb134_half(%esp),%xmm7   ## rinv iH1 - j water 
        movaps %xmm7,nb134_rinvMM(%esp)


        ## LJ table interaction
        movaps nb134_rinvOO(%esp),%xmm0
        movss %xmm0,%xmm1
        mulss  nb134_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulss  nb134_tsc(%esp),%xmm1

    cvttps2pi %xmm1,%mm6
    cvtpi2ps %mm6,%xmm3
        subss    %xmm3,%xmm1    ## xmm1=eps 
    movss %xmm1,%xmm2
    mulss  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $3,%mm6

    movd %eax,%mm0

    movl nb134_VFtab(%ebp),%esi
    movd %mm6,%eax

    ## dispersion 
    movlps (%esi,%eax,4),%xmm5
    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 8(%esi,%eax,4),%xmm7
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## dispersion table ready, in xmm4-xmm7 
    mulss  %xmm1,%xmm6      ## xmm6=Geps 
    mulss  %xmm2,%xmm7      ## xmm7=Heps2 
    addss  %xmm6,%xmm5
    addss  %xmm7,%xmm5      ## xmm5=Fp 
    mulss  nb134_two(%esp),%xmm7         ## two*Heps2 
    addss  %xmm6,%xmm7
    addss  %xmm5,%xmm7 ## xmm7=FF 
    mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addss  %xmm4,%xmm5 ## xmm5=VV 

    movss nb134_c6(%esp),%xmm4
    mulss  %xmm4,%xmm7   ## fijD 
    mulss  %xmm4,%xmm5   ## Vvdw6 
        xorps  %xmm3,%xmm3

        mulps  nb134_tsc(%esp),%xmm7
        subss  %xmm7,%xmm3
        movss  %xmm3,nb134_fstmp(%esp)

    addss  nb134_Vvdwtot(%esp),%xmm5
    movss %xmm5,nb134_Vvdwtot(%esp)

    ## repulsion 
    movlps 16(%esi,%eax,4),%xmm5
    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 24(%esi,%eax,4),%xmm7
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## table ready, in xmm4-xmm7 
    mulss  %xmm1,%xmm6      ## xmm6=Geps 
    mulss  %xmm2,%xmm7      ## xmm7=Heps2 
    addss  %xmm6,%xmm5
    addss  %xmm7,%xmm5      ## xmm5=Fp 
    mulss  nb134_two(%esp),%xmm7         ## two*Heps2 
    addss  %xmm6,%xmm7
    addss  %xmm5,%xmm7 ## xmm7=FF 
    mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addss  %xmm4,%xmm5 ## xmm5=VV 

    movss nb134_c12(%esp),%xmm4
    mulss  %xmm4,%xmm7 ## fijR 
    mulss  %xmm4,%xmm5 ## Vvdw12 
        movaps nb134_fstmp(%esp),%xmm3
        mulss  nb134_tsc(%esp),%xmm7
        subss  %xmm7,%xmm3

    addss  nb134_Vvdwtot(%esp),%xmm5
    movss %xmm5,nb134_Vvdwtot(%esp)

        mulss %xmm3,%xmm0
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movd  %mm0,%eax

        mulss  nb134_dxOO(%esp),%xmm0
        mulss  nb134_dyOO(%esp),%xmm1
        mulss  nb134_dzOO(%esp),%xmm2
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        subss   %xmm0,%xmm3
        subss   %xmm1,%xmm4
        subss   %xmm2,%xmm5
        movaps  %xmm3,nb134_fjxO(%esp)
        movaps  %xmm4,nb134_fjyO(%esp)
        movaps  %xmm5,nb134_fjzO(%esp)
        addss   nb134_fixO(%esp),%xmm0
        addss   nb134_fiyO(%esp),%xmm1
        addss   nb134_fizO(%esp),%xmm2
        movss  %xmm0,nb134_fixO(%esp)
        movss  %xmm1,nb134_fiyO(%esp)
        movss  %xmm2,nb134_fizO(%esp)

        ## do  M coulomb interaction
        movaps nb134_rinvMM(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        mulps  %xmm0,%xmm0      ## xmm0=rinvsq 

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb134_qqMH(%esp),%xmm3
        movhps  nb134_qqMM(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## constant 11000001 

        mulps  %xmm3,%xmm7 ## vcoul=rinv*qq
        movaps %xmm7,%xmm6
        addps  nb134_vctot(%esp),%xmm7
        movaps %xmm7,nb134_vctot(%esp)
        mulps  %xmm6,%xmm0

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        mulps   nb134_dxMM(%esp),%xmm0
        mulps   nb134_dyMM(%esp),%xmm1
        mulps   nb134_dzMM(%esp),%xmm2
        ## update forces M - j water 
        movaps  nb134_fjxO(%esp),%xmm3
        movaps  nb134_fjyO(%esp),%xmm4
        movaps  nb134_fjzO(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb134_fjxO(%esp)
        movaps  %xmm4,nb134_fjyO(%esp)
        movaps  %xmm5,nb134_fjzO(%esp)
        addps   nb134_fixM(%esp),%xmm0
        addps   nb134_fiyM(%esp),%xmm1
        addps   nb134_fizM(%esp),%xmm2
        movaps  %xmm0,nb134_fixM(%esp)
        movaps  %xmm1,nb134_fiyM(%esp)
        movaps  %xmm2,nb134_fizM(%esp)

        ## i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb134_ixH1(%esp),%xmm0
        movaps  nb134_iyH1(%esp),%xmm1
        movaps  nb134_izH1(%esp),%xmm2
        movaps  nb134_ixH2(%esp),%xmm3
        movaps  nb134_iyH2(%esp),%xmm4
        movaps  nb134_izH2(%esp),%xmm5
        subps   nb134_jxO(%esp),%xmm0
        subps   nb134_jyO(%esp),%xmm1
        subps   nb134_jzO(%esp),%xmm2
        subps   nb134_jxO(%esp),%xmm3
        subps   nb134_jyO(%esp),%xmm4
        subps   nb134_jzO(%esp),%xmm5
        movaps %xmm0,nb134_dxH1H1(%esp)
        movaps %xmm1,nb134_dyH1H1(%esp)
        movaps %xmm2,nb134_dzH1H1(%esp)
        movaps %xmm3,nb134_dxH2H2(%esp)
        movaps %xmm4,nb134_dyH2H2(%esp)
        movaps %xmm5,nb134_dzH2H2(%esp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm1,%xmm0
        addps %xmm3,%xmm4
        addps %xmm2,%xmm0       ## have rsqH1 in xmm0 
        addps %xmm5,%xmm4       ## have rsqH2 in xmm4 
        movaps  %xmm0,nb134_rsqH1H1(%esp)
        movaps  %xmm4,nb134_rsqH2H2(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb134_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb134_half(%esp),%xmm3   ## rinvH1H1
        mulps   nb134_half(%esp),%xmm7   ## rinvH2H2
        movaps  %xmm3,nb134_rinvH1H1(%esp)
        movaps  %xmm7,nb134_rinvH2H2(%esp)

        ## Do H1 coulomb interaction
        movaps %xmm3,%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        mulps  %xmm0,%xmm0      ## xmm0=rinvsq 

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb134_qqHH(%esp),%xmm3
        movhps  nb134_qqMH(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## constant 11000001 

        mulps  %xmm3,%xmm7 ##vcoul 
        mulps  %xmm7,%xmm0

        addps  nb134_vctot(%esp),%xmm7
        movaps %xmm7,nb134_vctot(%esp)

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        mulps   nb134_dxH1H1(%esp),%xmm0
        mulps   nb134_dyH1H1(%esp),%xmm1
        mulps   nb134_dzH1H1(%esp),%xmm2
        ## update forces H1 - j water 
        movaps  nb134_fjxO(%esp),%xmm3
        movaps  nb134_fjyO(%esp),%xmm4
        movaps  nb134_fjzO(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb134_fjxO(%esp)
        movaps  %xmm4,nb134_fjyO(%esp)
        movaps  %xmm5,nb134_fjzO(%esp)
        addps   nb134_fixH1(%esp),%xmm0
        addps   nb134_fiyH1(%esp),%xmm1
        addps   nb134_fizH1(%esp),%xmm2
        movaps  %xmm0,nb134_fixH1(%esp)
        movaps  %xmm1,nb134_fiyH1(%esp)
        movaps  %xmm2,nb134_fizH1(%esp)

        ## H2 Coulomb
        movaps nb134_rinvH2H2(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        mulps  %xmm0,%xmm0      ## xmm0=rinvsq 

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb134_qqHH(%esp),%xmm3
        movhps  nb134_qqMH(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## constant 11000001 

        mulps  %xmm3,%xmm7
        mulps  %xmm7,%xmm0

        addps  nb134_vctot(%esp),%xmm7   ## local vctot summation variable
        movaps %xmm7,nb134_vctot(%esp)

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        mulps   nb134_dxH2H2(%esp),%xmm0
        mulps   nb134_dyH2H2(%esp),%xmm1
        mulps   nb134_dzH2H2(%esp),%xmm2
        ## update forces H2 - j water 
        movaps  nb134_fjxO(%esp),%xmm3
        movaps  nb134_fjyO(%esp),%xmm4
        movaps  nb134_fjzO(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb134_fjxO(%esp)
        movaps  %xmm4,nb134_fjyO(%esp)
        movaps  %xmm5,nb134_fjzO(%esp)
        addps   nb134_fixH2(%esp),%xmm0
        addps   nb134_fiyH2(%esp),%xmm1
        addps   nb134_fizH2(%esp),%xmm2
        movaps  %xmm0,nb134_fixH2(%esp)
        movaps  %xmm1,nb134_fiyH2(%esp)
        movaps  %xmm2,nb134_fizH2(%esp)

        movl    nb134_faction(%ebp),%esi
        ## update j water forces from local variables.
        ## transpose back first
        movaps  nb134_fjxO(%esp),%xmm0   ## Ox H1x H2x Mx 
        movaps  nb134_fjyO(%esp),%xmm1   ## Oy H1y H2y My
        movaps  nb134_fjzO(%esp),%xmm2   ## Oz H1z H2z Mz

        movaps  %xmm0,%xmm3
        movaps  %xmm0,%xmm4
        unpcklps %xmm1,%xmm3            ## Ox Oy - -
        shufps $0x1,%xmm2,%xmm4        ## h1x - Oz -
        movaps  %xmm1,%xmm5
        movaps  %xmm0,%xmm6
        unpcklps %xmm2,%xmm5            ## - - H1y H1z
        unpckhps %xmm1,%xmm6            ## h2x h2y - - 
        unpckhps %xmm2,%xmm1            ## - - My Mz

        shufps  $50,%xmm0,%xmm2 ## (00110010) h2z - Mx -
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

        decl nb134_innerk(%esp)
        jz    _nb_kernel134_ia32_sse.nb134_updateouterdata
        jmp   _nb_kernel134_ia32_sse.nb134_single_loop
_nb_kernel134_ia32_sse.nb134_updateouterdata: 
        movl  nb134_ii3(%esp),%ecx
        movl  nb134_faction(%ebp),%edi
        movl  nb134_fshift(%ebp),%esi
        movl  nb134_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb134_fixO(%esp),%xmm0
        movaps nb134_fiyO(%esp),%xmm1
        movaps nb134_fizO(%esp),%xmm2

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
        movaps nb134_fixH1(%esp),%xmm0
        movaps nb134_fiyH1(%esp),%xmm1
        movaps nb134_fizH1(%esp),%xmm2

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
        movaps nb134_fixH2(%esp),%xmm0
        movaps nb134_fiyH2(%esp),%xmm1
        movaps nb134_fizH2(%esp),%xmm2

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
        movaps nb134_fixM(%esp),%xmm0
        movaps nb134_fiyM(%esp),%xmm1
        movaps nb134_fizM(%esp),%xmm2

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
        movl nb134_n(%esp),%esi
        ## get group index for i particle 
        movl  nb134_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb134_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb134_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb134_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb134_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb134_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel134_ia32_sse.nb134_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb134_n(%esp)
        jmp _nb_kernel134_ia32_sse.nb134_outer
_nb_kernel134_ia32_sse.nb134_outerend: 
        ## check if more outer neighborlists remain
        movl  nb134_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel134_ia32_sse.nb134_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel134_ia32_sse.nb134_threadloop
_nb_kernel134_ia32_sse.nb134_end: 
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



.globl nb_kernel134nf_ia32_sse
.globl _nb_kernel134nf_ia32_sse
nb_kernel134nf_ia32_sse:        
_nb_kernel134nf_ia32_sse:       
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
        ## bottom of stack is cache-aligned for sse use 
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
.set nb134nf_c6, 432
.set nb134nf_c12, 448
.set nb134nf_tsc, 464
.set nb134nf_vctot, 480
.set nb134nf_Vvdwtot, 496
.set nb134nf_three, 512
.set nb134nf_rsqOO, 528
.set nb134nf_rsqH1H1, 544
.set nb134nf_rsqH1H2, 560
.set nb134nf_rsqH1M, 576
.set nb134nf_rsqH2H1, 592
.set nb134nf_rsqH2H2, 608
.set nb134nf_rsqH2M, 624
.set nb134nf_rsqMH1, 640
.set nb134nf_rsqMH2, 656
.set nb134nf_rsqMM, 672
.set nb134nf_rinvOO, 688
.set nb134nf_rinvH1H1, 704
.set nb134nf_rinvH1H2, 720
.set nb134nf_rinvH1M, 736
.set nb134nf_rinvH2H1, 752
.set nb134nf_rinvH2H2, 768
.set nb134nf_rinvH2M, 784
.set nb134nf_rinvMH1, 800
.set nb134nf_rinvMH2, 816
.set nb134nf_rinvMM, 832
.set nb134nf_half, 880
.set nb134nf_is3, 896
.set nb134nf_ii3, 900
.set nb134nf_innerjjnr, 904
.set nb134nf_innerk, 908
.set nb134nf_n, 912
.set nb134nf_nn1, 916
.set nb134nf_nri, 920
.set nb134nf_nouter, 924
.set nb134nf_ninner, 928
.set nb134nf_salign, 932
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $936,%esp          ## local stack space 
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
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb134nf_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb134nf_half(%esp)
        movss nb134nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb134nf_half(%esp)
        movaps %xmm3,nb134nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb134nf_iinr(%ebp),%ecx     ## ecx = pointer into iinr[]
        movl  (%ecx),%ebx               ## ebx =ii 

        movl  nb134nf_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm5
        movss 12(%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movl nb134nf_p_facel(%ebp),%esi
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
        movaps %xmm3,nb134nf_qqMM(%esp)
        movaps %xmm4,nb134nf_qqMH(%esp)
        movaps %xmm5,nb134nf_qqHH(%esp)

        xorps %xmm0,%xmm0
        movl  nb134nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb134nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb134nf_vdwparam(%ebp),%eax
        movlps (%eax,%edx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $0x55,%xmm1,%xmm1
        movaps %xmm0,nb134nf_c6(%esp)
        movaps %xmm1,nb134nf_c12(%esp)

_nb_kernel134nf_ia32_sse.nb134nf_threadloop: 
        movl  nb134nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel134nf_ia32_sse.nb134nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel134nf_ia32_sse.nb134nf_spinlock

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
        jg  _nb_kernel134nf_ia32_sse.nb134nf_outerstart
        jmp _nb_kernel134nf_ia32_sse.nb134nf_end

_nb_kernel134nf_ia32_sse.nb134nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb134nf_nouter(%esp),%ebx
        movl %ebx,nb134nf_nouter(%esp)

_nb_kernel134nf_ia32_sse.nb134nf_outer: 
        movl  nb134nf_shift(%ebp),%eax          ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb134nf_is3(%esp)            ## store is3 

        movl  nb134nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb134nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb134nf_pos(%ebp),%eax    ## eax = base of pos[]  
        movl  %ebx,nb134nf_ii3(%esp)

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
        movaps %xmm3,nb134nf_ixO(%esp)
        movaps %xmm4,nb134nf_iyO(%esp)
        movaps %xmm5,nb134nf_izO(%esp)
        movaps %xmm6,nb134nf_ixH1(%esp)
        movaps %xmm7,nb134nf_iyH1(%esp)

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
        movaps %xmm6,nb134nf_izH1(%esp)
        movaps %xmm0,nb134nf_ixH2(%esp)
        movaps %xmm1,nb134nf_iyH2(%esp)
        movaps %xmm2,nb134nf_izH2(%esp)
        movaps %xmm3,nb134nf_ixM(%esp)
        movaps %xmm4,nb134nf_iyM(%esp)
        movaps %xmm5,nb134nf_izM(%esp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb134nf_vctot(%esp)
        movaps %xmm4,nb134nf_Vvdwtot(%esp)

        movl  nb134nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb134nf_pos(%ebp),%esi
        movl  nb134nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb134nf_innerjjnr(%esp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb134nf_ninner(%esp),%ecx
        movl  %ecx,nb134nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb134nf_innerk(%esp)         ## number of innerloop atoms 
        jge   _nb_kernel134nf_ia32_sse.nb134nf_unroll_loop
        jmp   _nb_kernel134nf_ia32_sse.nb134nf_single_check
_nb_kernel134nf_ia32_sse.nb134nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb134nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb134nf_innerjjnr(%esp)             ## advance pointer (unroll 4) 

        movl nb134nf_pos(%ebp),%esi     ## base of pos[] 

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
        movaps %xmm0,nb134nf_jxO(%esp)
        movaps %xmm1,nb134nf_jyO(%esp)
        movaps %xmm2,nb134nf_jzO(%esp)
        movaps %xmm3,nb134nf_jxH1(%esp)

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
        movaps %xmm0,nb134nf_jyH1(%esp)
        movaps %xmm1,nb134nf_jzH1(%esp)
        movaps %xmm2,nb134nf_jxH2(%esp)
        movaps %xmm3,nb134nf_jyH2(%esp)

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
        movaps %xmm0,nb134nf_jzH2(%esp)
        movaps %xmm1,nb134nf_jxM(%esp)
        movaps %xmm2,nb134nf_jyM(%esp)
        movaps %xmm3,nb134nf_jzM(%esp)

        ## start calculating pairwise distances
        movaps nb134nf_ixO(%esp),%xmm0
        movaps nb134nf_iyO(%esp),%xmm1
        movaps nb134nf_izO(%esp),%xmm2
        movaps nb134nf_ixH1(%esp),%xmm3
        movaps nb134nf_iyH1(%esp),%xmm4
        movaps nb134nf_izH1(%esp),%xmm5
        subps  nb134nf_jxO(%esp),%xmm0
        subps  nb134nf_jyO(%esp),%xmm1
        subps  nb134nf_jzO(%esp),%xmm2
        subps  nb134nf_jxH1(%esp),%xmm3
        subps  nb134nf_jyH1(%esp),%xmm4
        subps  nb134nf_jzH1(%esp),%xmm5
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
        movaps %xmm0,nb134nf_rsqOO(%esp)
        movaps %xmm3,nb134nf_rsqH1H1(%esp)

        movaps nb134nf_ixH1(%esp),%xmm0
        movaps nb134nf_iyH1(%esp),%xmm1
        movaps nb134nf_izH1(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb134nf_jxH2(%esp),%xmm0
        subps  nb134nf_jyH2(%esp),%xmm1
        subps  nb134nf_jzH2(%esp),%xmm2
        subps  nb134nf_jxM(%esp),%xmm3
        subps  nb134nf_jyM(%esp),%xmm4
        subps  nb134nf_jzM(%esp),%xmm5
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
        movaps %xmm0,nb134nf_rsqH1H2(%esp)
        movaps %xmm3,nb134nf_rsqH1M(%esp)

        movaps nb134nf_ixH2(%esp),%xmm0
        movaps nb134nf_iyH2(%esp),%xmm1
        movaps nb134nf_izH2(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb134nf_jxH1(%esp),%xmm0
        subps  nb134nf_jyH1(%esp),%xmm1
        subps  nb134nf_jzH1(%esp),%xmm2
        subps  nb134nf_jxH2(%esp),%xmm3
        subps  nb134nf_jyH2(%esp),%xmm4
        subps  nb134nf_jzH2(%esp),%xmm5
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
        movaps %xmm0,nb134nf_rsqH2H1(%esp)
        movaps %xmm3,nb134nf_rsqH2H2(%esp)

        movaps nb134nf_ixH2(%esp),%xmm0
        movaps nb134nf_iyH2(%esp),%xmm1
        movaps nb134nf_izH2(%esp),%xmm2
        movaps nb134nf_ixM(%esp),%xmm3
        movaps nb134nf_iyM(%esp),%xmm4
        movaps nb134nf_izM(%esp),%xmm5
        subps  nb134nf_jxM(%esp),%xmm0
        subps  nb134nf_jyM(%esp),%xmm1
        subps  nb134nf_jzM(%esp),%xmm2
        subps  nb134nf_jxH1(%esp),%xmm3
        subps  nb134nf_jyH1(%esp),%xmm4
        subps  nb134nf_jzH1(%esp),%xmm5
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
        movaps %xmm0,nb134nf_rsqH2M(%esp)
        movaps %xmm4,nb134nf_rsqMH1(%esp)

        movaps nb134nf_ixM(%esp),%xmm0
        movaps nb134nf_iyM(%esp),%xmm1
        movaps nb134nf_izM(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb134nf_jxH2(%esp),%xmm0
        subps  nb134nf_jyH2(%esp),%xmm1
        subps  nb134nf_jzH2(%esp),%xmm2
        subps  nb134nf_jxM(%esp),%xmm3
        subps  nb134nf_jyM(%esp),%xmm4
        subps  nb134nf_jzM(%esp),%xmm5
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
        movaps %xmm0,nb134nf_rsqMH2(%esp)
        movaps %xmm4,nb134nf_rsqMM(%esp)

        ## start by doing invsqrt for OO
        rsqrtps nb134nf_rsqOO(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb134nf_three(%esp),%xmm3
        mulps   nb134nf_rsqOO(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb134nf_half(%esp),%xmm3
        movaps  %xmm3,nb134nf_rinvOO(%esp)

        ## more invsqrt ops - do two at a time.
        rsqrtps nb134nf_rsqH1H1(%esp),%xmm1
        rsqrtps nb134nf_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb134nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb134nf_rsqH1H1(%esp),%xmm1
        mulps   nb134nf_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb134nf_half(%esp),%xmm3   ## rinvH1H1 
        mulps   nb134nf_half(%esp),%xmm7   ## rinvH1H2 
        movaps  %xmm3,nb134nf_rinvH1H1(%esp)
        movaps  %xmm7,nb134nf_rinvH1H2(%esp)

        rsqrtps nb134nf_rsqH1M(%esp),%xmm1
        rsqrtps nb134nf_rsqH2H1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb134nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb134nf_rsqH1M(%esp),%xmm1
        mulps   nb134nf_rsqH2H1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb134nf_half(%esp),%xmm3
        mulps   nb134nf_half(%esp),%xmm7
        movaps  %xmm3,nb134nf_rinvH1M(%esp)
        movaps  %xmm7,nb134nf_rinvH2H1(%esp)

        rsqrtps nb134nf_rsqH2H2(%esp),%xmm1
        rsqrtps nb134nf_rsqH2M(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb134nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb134nf_rsqH2H2(%esp),%xmm1
        mulps   nb134nf_rsqH2M(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb134nf_half(%esp),%xmm3
        mulps   nb134nf_half(%esp),%xmm7
        movaps  %xmm3,nb134nf_rinvH2H2(%esp)
        movaps  %xmm7,nb134nf_rinvH2M(%esp)

        rsqrtps nb134nf_rsqMH1(%esp),%xmm1
        rsqrtps nb134nf_rsqMH2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb134nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb134nf_rsqMH1(%esp),%xmm1
        mulps   nb134nf_rsqMH2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb134nf_half(%esp),%xmm3
        mulps   nb134nf_half(%esp),%xmm7
        movaps  %xmm3,nb134nf_rinvMH1(%esp)
        movaps  %xmm7,nb134nf_rinvMH2(%esp)

        rsqrtps nb134nf_rsqMM(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb134nf_three(%esp),%xmm3
        mulps   nb134nf_rsqMM(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb134nf_half(%esp),%xmm3
        movaps  %xmm3,nb134nf_rinvMM(%esp)

        ## start with OO LJ interaction
        movaps nb134nf_rinvOO(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb134nf_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulps  nb134nf_tsc(%esp),%xmm1

        movhlps %xmm1,%xmm2
    cvttps2pi %xmm1,%mm6
    cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
    movaps %xmm1,%xmm2
    mulps  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $3,%mm6
    pslld $3,%mm7

    movd %mm6,%eax
    psrlq $32,%mm6
    movd %mm7,%ecx
    psrlq $32,%mm7
    movd %mm6,%ebx
    movd %mm7,%edx

    movl nb134nf_VFtab(%ebp),%esi

    ## dispersion 
    movlps (%esi,%eax,4),%xmm5
    movlps (%esi,%ecx,4),%xmm7
    movhps (%esi,%ebx,4),%xmm5
    movhps (%esi,%edx,4),%xmm7 ## got half table 

    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 8(%esi,%eax,4),%xmm7
    movlps 8(%esi,%ecx,4),%xmm3
    movhps 8(%esi,%ebx,4),%xmm7
    movhps 8(%esi,%edx,4),%xmm3    ## other half of table  
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## dispersion table ready, in xmm4-xmm7 
    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 

    movaps nb134nf_c6(%esp),%xmm4
    mulps  %xmm4,%xmm5   ## Vvdw6 

    addps  nb134nf_Vvdwtot(%esp),%xmm5
    movaps %xmm5,nb134nf_Vvdwtot(%esp)

    ## repulsion 
    movlps 16(%esi,%eax,4),%xmm5
    movlps 16(%esi,%ecx,4),%xmm7
    movhps 16(%esi,%ebx,4),%xmm5
    movhps 16(%esi,%edx,4),%xmm7    ## got half table 

    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 24(%esi,%eax,4),%xmm7
    movlps 24(%esi,%ecx,4),%xmm3
    movhps 24(%esi,%ebx,4),%xmm7
    movhps 24(%esi,%edx,4),%xmm3    ## other half of table  
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## repulsion table ready, in xmm4-xmm7 
    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 

    movaps nb134nf_c12(%esp),%xmm4
    mulps  %xmm4,%xmm5 ## Vvdw12 

    addps  nb134nf_Vvdwtot(%esp),%xmm5
    movaps %xmm5,nb134nf_Vvdwtot(%esp)

        ## Coulomb interactions 
        movaps nb134nf_rinvH1H1(%esp),%xmm0
        movaps nb134nf_rinvH1M(%esp),%xmm1
        movaps nb134nf_rinvMM(%esp),%xmm2
        addps  nb134nf_rinvH1H2(%esp),%xmm0
        addps  nb134nf_rinvH2M(%esp),%xmm1
        addps  nb134nf_rinvH2H1(%esp),%xmm0
        addps  nb134nf_rinvMH1(%esp),%xmm1
        addps  nb134nf_rinvH2H2(%esp),%xmm0
        addps  nb134nf_rinvMH2(%esp),%xmm1

        mulps  nb134nf_qqHH(%esp),%xmm0
        mulps  nb134nf_qqMH(%esp),%xmm1
        mulps  nb134nf_qqMM(%esp),%xmm2

        addps  %xmm1,%xmm0
        addps  nb134nf_vctot(%esp),%xmm2
        addps  %xmm2,%xmm0
        movaps %xmm0,nb134nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb134nf_innerk(%esp)
        jl    _nb_kernel134nf_ia32_sse.nb134nf_single_check
        jmp   _nb_kernel134nf_ia32_sse.nb134nf_unroll_loop
_nb_kernel134nf_ia32_sse.nb134nf_single_check: 
        addl $4,nb134nf_innerk(%esp)
        jnz   _nb_kernel134nf_ia32_sse.nb134nf_single_loop
        jmp   _nb_kernel134nf_ia32_sse.nb134nf_updateouterdata
_nb_kernel134nf_ia32_sse.nb134nf_single_loop: 
        movl  nb134nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb134nf_innerjjnr(%esp)

        movl nb134nf_pos(%ebp),%esi
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
        movaps %xmm6,nb134nf_jxO(%esp)
        movaps %xmm3,nb134nf_jyO(%esp)
        movaps %xmm1,nb134nf_jzO(%esp)

        ## do O and M in parallel
        movaps nb134nf_ixO(%esp),%xmm0
        movaps nb134nf_iyO(%esp),%xmm1
        movaps nb134nf_izO(%esp),%xmm2
        movaps nb134nf_ixM(%esp),%xmm3
        movaps nb134nf_iyM(%esp),%xmm4
        movaps nb134nf_izM(%esp),%xmm5
        subps  nb134nf_jxO(%esp),%xmm0
        subps  nb134nf_jyO(%esp),%xmm1
        subps  nb134nf_jzO(%esp),%xmm2
        subps  nb134nf_jxO(%esp),%xmm3
        subps  nb134nf_jyO(%esp),%xmm4
        subps  nb134nf_jzO(%esp),%xmm5

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
        ## Save data 
        movaps %xmm0,nb134nf_rsqOO(%esp)
        movaps %xmm4,nb134nf_rsqMM(%esp)

        ## do 1/x for O
        rsqrtps %xmm0,%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb134nf_three(%esp),%xmm3
        mulps   %xmm0,%xmm1     ## rsq*lu*lu 
        subps   %xmm1,%xmm3     ## constant 30-rsq*lu*lu 
        mulps   %xmm2,%xmm3     ## lu*(3-rsq*lu*lu) 
        mulps   nb134nf_half(%esp),%xmm3
        movaps  %xmm3,nb134nf_rinvOO(%esp)      ## rinvH2 

        ## 1/sqrt(x) for M
        rsqrtps %xmm4,%xmm5
        movaps  %xmm5,%xmm6
        mulps   %xmm5,%xmm5
        movaps  nb134nf_three(%esp),%xmm7
        mulps   %xmm4,%xmm5
        subps   %xmm5,%xmm7
        mulps   %xmm6,%xmm7
        mulps   nb134nf_half(%esp),%xmm7   ## rinv iH1 - j water 
        movaps %xmm7,nb134nf_rinvMM(%esp)


        ## LJ table interaction
        movaps nb134nf_rinvOO(%esp),%xmm0
        movss %xmm0,%xmm1
        mulss  nb134nf_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulss  nb134nf_tsc(%esp),%xmm1

    cvttps2pi %xmm1,%mm6
    cvtpi2ps %mm6,%xmm3
        subss    %xmm3,%xmm1    ## xmm1=eps 
    movss %xmm1,%xmm2
    mulss  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $3,%mm6

    movl nb134nf_VFtab(%ebp),%esi
    movd %mm6,%eax

    ## dispersion 
    movlps (%esi,%eax,4),%xmm5
    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 8(%esi,%eax,4),%xmm7
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## dispersion table ready, in xmm4-xmm7 
    mulss  %xmm1,%xmm6      ## xmm6=Geps 
    mulss  %xmm2,%xmm7      ## xmm7=Heps2 
    addss  %xmm6,%xmm5
    addss  %xmm7,%xmm5      ## xmm5=Fp 
    mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addss  %xmm4,%xmm5 ## xmm5=VV 

    movss nb134nf_c6(%esp),%xmm4
    mulss  %xmm4,%xmm5   ## Vvdw6 

    addss  nb134nf_Vvdwtot(%esp),%xmm5
    movss %xmm5,nb134nf_Vvdwtot(%esp)

    ## repulsion 
    movlps 16(%esi,%eax,4),%xmm5
    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 24(%esi,%eax,4),%xmm7
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## table ready, in xmm4-xmm7 
    mulss  %xmm1,%xmm6      ## xmm6=Geps 
    mulss  %xmm2,%xmm7      ## xmm7=Heps2 
    addss  %xmm6,%xmm5
    addss  %xmm7,%xmm5      ## xmm5=Fp 
    mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addss  %xmm4,%xmm5 ## xmm5=VV 

    movss nb134nf_c12(%esp),%xmm4
    mulss  %xmm4,%xmm5 ## Vvdw12 
    addss  nb134nf_Vvdwtot(%esp),%xmm5
    movss %xmm5,nb134nf_Vvdwtot(%esp)

        ## do  M coulomb interaction
        movaps nb134nf_rinvMM(%esp),%xmm0

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb134nf_qqMH(%esp),%xmm3
        movhps  nb134nf_qqMM(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## constant 11000001 

        mulps  %xmm3,%xmm0
        addps  nb134nf_vctot(%esp),%xmm0
        movaps %xmm0,nb134nf_vctot(%esp)

        ## i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb134nf_ixH1(%esp),%xmm0
        movaps  nb134nf_iyH1(%esp),%xmm1
        movaps  nb134nf_izH1(%esp),%xmm2
        movaps  nb134nf_ixH2(%esp),%xmm3
        movaps  nb134nf_iyH2(%esp),%xmm4
        movaps  nb134nf_izH2(%esp),%xmm5
        subps   nb134nf_jxO(%esp),%xmm0
        subps   nb134nf_jyO(%esp),%xmm1
        subps   nb134nf_jzO(%esp),%xmm2
        subps   nb134nf_jxO(%esp),%xmm3
        subps   nb134nf_jyO(%esp),%xmm4
        subps   nb134nf_jzO(%esp),%xmm5
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm1,%xmm0
        addps %xmm3,%xmm4
        addps %xmm2,%xmm0       ## have rsqH1 in xmm0 
        addps %xmm5,%xmm4       ## have rsqH2 in xmm4 
        movaps  %xmm0,nb134nf_rsqH1H1(%esp)
        movaps  %xmm4,nb134nf_rsqH2H2(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb134nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb134nf_half(%esp),%xmm3   ## rinvH1H1
        mulps   nb134nf_half(%esp),%xmm7   ## rinvH2H2
        movaps  %xmm3,nb134nf_rinvH1H1(%esp)
        movaps  %xmm7,nb134nf_rinvH2H2(%esp)


        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb134nf_qqHH(%esp),%xmm3
        movhps  nb134nf_qqMH(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## constant 11000001 

        ## Do coulomb interactions      
        movaps nb134nf_rinvH1H1(%esp),%xmm0
        addps  nb134nf_rinvH2H2(%esp),%xmm0
        mulps  %xmm3,%xmm0
        addps  nb134nf_vctot(%esp),%xmm0
        movaps %xmm0,nb134nf_vctot(%esp)

        decl nb134nf_innerk(%esp)
        jz    _nb_kernel134nf_ia32_sse.nb134nf_updateouterdata
        jmp   _nb_kernel134nf_ia32_sse.nb134nf_single_loop
_nb_kernel134nf_ia32_sse.nb134nf_updateouterdata: 
        ## get n from stack
        movl nb134nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb134nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb134nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb134nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb134nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb134nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb134nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel134nf_ia32_sse.nb134nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb134nf_n(%esp)
        jmp _nb_kernel134nf_ia32_sse.nb134nf_outer
_nb_kernel134nf_ia32_sse.nb134nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb134nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel134nf_ia32_sse.nb134nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel134nf_ia32_sse.nb134nf_threadloop
_nb_kernel134nf_ia32_sse.nb134nf_end: 
        emms

        movl nb134nf_nouter(%esp),%eax
        movl nb134nf_ninner(%esp),%ebx
        movl nb134nf_outeriter(%ebp),%ecx
        movl nb134nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb134nf_salign(%esp),%eax
        addl %eax,%esp
        addl $936,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


