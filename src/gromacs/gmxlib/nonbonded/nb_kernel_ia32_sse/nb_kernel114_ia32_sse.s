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



.globl nb_kernel114_ia32_sse
.globl _nb_kernel114_ia32_sse
nb_kernel114_ia32_sse:  
_nb_kernel114_ia32_sse: 
.set nb114_p_nri, 8
.set nb114_iinr, 12
.set nb114_jindex, 16
.set nb114_jjnr, 20
.set nb114_shift, 24
.set nb114_shiftvec, 28
.set nb114_fshift, 32
.set nb114_gid, 36
.set nb114_pos, 40
.set nb114_faction, 44
.set nb114_charge, 48
.set nb114_p_facel, 52
.set nb114_p_krf, 56
.set nb114_p_crf, 60
.set nb114_Vc, 64
.set nb114_type, 68
.set nb114_p_ntype, 72
.set nb114_vdwparam, 76
.set nb114_Vvdw, 80
.set nb114_p_tabscale, 84
.set nb114_VFtab, 88
.set nb114_invsqrta, 92
.set nb114_dvda, 96
.set nb114_p_gbtabscale, 100
.set nb114_GBtab, 104
.set nb114_p_nthreads, 108
.set nb114_count, 112
.set nb114_mtx, 116
.set nb114_outeriter, 120
.set nb114_inneriter, 124
.set nb114_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb114_ixO, 0
.set nb114_iyO, 16
.set nb114_izO, 32
.set nb114_ixH1, 48
.set nb114_iyH1, 64
.set nb114_izH1, 80
.set nb114_ixH2, 96
.set nb114_iyH2, 112
.set nb114_izH2, 128
.set nb114_ixM, 144
.set nb114_iyM, 160
.set nb114_izM, 176
.set nb114_jxO, 192
.set nb114_jyO, 208
.set nb114_jzO, 224
.set nb114_jxH1, 240
.set nb114_jyH1, 256
.set nb114_jzH1, 272
.set nb114_jxH2, 288
.set nb114_jyH2, 304
.set nb114_jzH2, 320
.set nb114_jxM, 336
.set nb114_jyM, 352
.set nb114_jzM, 368
.set nb114_dxOO, 384
.set nb114_dyOO, 400
.set nb114_dzOO, 416
.set nb114_dxH1H1, 432
.set nb114_dyH1H1, 448
.set nb114_dzH1H1, 464
.set nb114_dxH1H2, 480
.set nb114_dyH1H2, 496
.set nb114_dzH1H2, 512
.set nb114_dxH1M, 528
.set nb114_dyH1M, 544
.set nb114_dzH1M, 560
.set nb114_dxH2H1, 576
.set nb114_dyH2H1, 592
.set nb114_dzH2H1, 608
.set nb114_dxH2H2, 624
.set nb114_dyH2H2, 640
.set nb114_dzH2H2, 656
.set nb114_dxH2M, 672
.set nb114_dyH2M, 688
.set nb114_dzH2M, 704
.set nb114_dxMH1, 720
.set nb114_dyMH1, 736
.set nb114_dzMH1, 752
.set nb114_dxMH2, 768
.set nb114_dyMH2, 784
.set nb114_dzMH2, 800
.set nb114_dxMM, 816
.set nb114_dyMM, 832
.set nb114_dzMM, 848
.set nb114_qqMM, 864
.set nb114_qqMH, 880
.set nb114_qqHH, 896
.set nb114_two, 912
.set nb114_c6, 928
.set nb114_c12, 944
.set nb114_six, 960
.set nb114_twelve, 976
.set nb114_vctot, 992
.set nb114_Vvdwtot, 1008
.set nb114_fixO, 1024
.set nb114_fiyO, 1040
.set nb114_fizO, 1056
.set nb114_fixH1, 1072
.set nb114_fiyH1, 1088
.set nb114_fizH1, 1104
.set nb114_fixH2, 1120
.set nb114_fiyH2, 1136
.set nb114_fizH2, 1152
.set nb114_fixM, 1168
.set nb114_fiyM, 1184
.set nb114_fizM, 1200
.set nb114_fjxO, 1216
.set nb114_fjyO, 1232
.set nb114_fjzO, 1248
.set nb114_fjxH1, 1264
.set nb114_fjyH1, 1280
.set nb114_fjzH1, 1296
.set nb114_fjxH2, 1312
.set nb114_fjyH2, 1328
.set nb114_fjzH2, 1344
.set nb114_fjxM, 1360
.set nb114_fjyM, 1376
.set nb114_fjzM, 1392
.set nb114_half, 1408
.set nb114_three, 1424
.set nb114_rsqOO, 1440
.set nb114_rsqH1H1, 1456
.set nb114_rsqH1H2, 1472
.set nb114_rsqH1M, 1488
.set nb114_rsqH2H1, 1504
.set nb114_rsqH2H2, 1520
.set nb114_rsqH2M, 1536
.set nb114_rsqMH1, 1552
.set nb114_rsqMH2, 1568
.set nb114_rsqMM, 1584
.set nb114_rinvsqOO, 1600
.set nb114_rinvH1H1, 1616
.set nb114_rinvH1H2, 1632
.set nb114_rinvH1M, 1648
.set nb114_rinvH2H1, 1664
.set nb114_rinvH2H2, 1680
.set nb114_rinvH2M, 1696
.set nb114_rinvMH1, 1712
.set nb114_rinvMH2, 1728
.set nb114_rinvMM, 1744
.set nb114_fstmp, 1760
.set nb114_is3, 1776
.set nb114_ii3, 1780
.set nb114_innerjjnr, 1784
.set nb114_innerk, 1788
.set nb114_n, 1792
.set nb114_nn1, 1796
.set nb114_nri, 1800
.set nb114_nouter, 1804
.set nb114_ninner, 1808
.set nb114_salign, 1812
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $1816,%esp         ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb114_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb114_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb114_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb114_nouter(%esp)
        movl %eax,nb114_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb114_half(%esp)
        movss nb114_half(%esp),%xmm1
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
        movaps %xmm1,nb114_half(%esp)
        movaps %xmm2,nb114_two(%esp)
        movaps %xmm3,nb114_three(%esp)
        movaps %xmm4,nb114_six(%esp)
        movaps %xmm5,nb114_twelve(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb114_iinr(%ebp),%ecx     ## ecx = pointer into iinr[]
        movl  (%ecx),%ebx               ## ebx =ii 

        movl  nb114_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm5
        movss 12(%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movl nb114_p_facel(%ebp),%esi
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
        movaps %xmm3,nb114_qqMM(%esp)
        movaps %xmm4,nb114_qqMH(%esp)
        movaps %xmm5,nb114_qqHH(%esp)

        xorps %xmm0,%xmm0
        movl  nb114_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb114_p_ntype(%ebp),%edi
        imull (%edi),%ecx ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb114_vdwparam(%ebp),%eax
        movlps (%eax,%edx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $0x55,%xmm1,%xmm1
        movaps %xmm0,nb114_c6(%esp)
        movaps %xmm1,nb114_c12(%esp)

_nb_kernel114_ia32_sse.nb114_threadloop: 
        movl  nb114_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel114_ia32_sse.nb114_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel114_ia32_sse.nb114_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb114_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb114_n(%esp)
        movl %ebx,nb114_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel114_ia32_sse.nb114_outerstart
        jmp _nb_kernel114_ia32_sse.nb114_end

_nb_kernel114_ia32_sse.nb114_outerstart: 
        ## ebx contains number of outer iterations
        addl nb114_nouter(%esp),%ebx
        movl %ebx,nb114_nouter(%esp)

_nb_kernel114_ia32_sse.nb114_outer: 
        movl  nb114_shift(%ebp),%eax            ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb114_is3(%esp)      ## store is3 

        movl  nb114_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb114_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb114_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb114_ii3(%esp)

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
        movaps %xmm3,nb114_ixO(%esp)
        movaps %xmm4,nb114_iyO(%esp)
        movaps %xmm5,nb114_izO(%esp)
        movaps %xmm6,nb114_ixH1(%esp)
        movaps %xmm7,nb114_iyH1(%esp)

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
        movaps %xmm6,nb114_izH1(%esp)
        movaps %xmm0,nb114_ixH2(%esp)
        movaps %xmm1,nb114_iyH2(%esp)
        movaps %xmm2,nb114_izH2(%esp)
        movaps %xmm3,nb114_ixM(%esp)
        movaps %xmm4,nb114_iyM(%esp)
        movaps %xmm5,nb114_izM(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb114_vctot(%esp)
        movaps %xmm4,nb114_Vvdwtot(%esp)
        movaps %xmm4,nb114_fixO(%esp)
        movaps %xmm4,nb114_fiyO(%esp)
        movaps %xmm4,nb114_fizO(%esp)
        movaps %xmm4,nb114_fixH1(%esp)
        movaps %xmm4,nb114_fiyH1(%esp)
        movaps %xmm4,nb114_fizH1(%esp)
        movaps %xmm4,nb114_fixH2(%esp)
        movaps %xmm4,nb114_fiyH2(%esp)
        movaps %xmm4,nb114_fizH2(%esp)
        movaps %xmm4,nb114_fixM(%esp)
        movaps %xmm4,nb114_fiyM(%esp)
        movaps %xmm4,nb114_fizM(%esp)

        movl  nb114_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb114_pos(%ebp),%esi
        movl  nb114_faction(%ebp),%edi
        movl  nb114_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb114_innerjjnr(%esp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb114_ninner(%esp),%ecx
        movl  %ecx,nb114_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb114_innerk(%esp)   ## number of innerloop atoms 
        jge   _nb_kernel114_ia32_sse.nb114_unroll_loop
        jmp   _nb_kernel114_ia32_sse.nb114_single_check
_nb_kernel114_ia32_sse.nb114_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb114_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb114_innerjjnr(%esp)             ## advance pointer (unroll 4) 

        movl nb114_pos(%ebp),%esi       ## base of pos[] 

        leal  (%eax,%eax,2),%eax        ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx        ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move j coordinates to local temp variables.
        ## Load Ox Oy Oz H1x
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
        movaps %xmm0,nb114_jxO(%esp)
        movaps %xmm1,nb114_jyO(%esp)
        movaps %xmm2,nb114_jzO(%esp)
        movaps %xmm3,nb114_jxH1(%esp)

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
        movaps %xmm0,nb114_jyH1(%esp)
        movaps %xmm1,nb114_jzH1(%esp)
        movaps %xmm2,nb114_jxH2(%esp)
        movaps %xmm3,nb114_jyH2(%esp)

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
        movaps %xmm0,nb114_jzH2(%esp)
        movaps %xmm1,nb114_jxM(%esp)
        movaps %xmm2,nb114_jyM(%esp)
        movaps %xmm3,nb114_jzM(%esp)

        ## start calculating pairwise distances
        movaps nb114_ixO(%esp),%xmm0
        movaps nb114_iyO(%esp),%xmm1
        movaps nb114_izO(%esp),%xmm2
        movaps nb114_ixH1(%esp),%xmm3
        movaps nb114_iyH1(%esp),%xmm4
        movaps nb114_izH1(%esp),%xmm5
        subps  nb114_jxO(%esp),%xmm0
        subps  nb114_jyO(%esp),%xmm1
        subps  nb114_jzO(%esp),%xmm2
        subps  nb114_jxH1(%esp),%xmm3
        subps  nb114_jyH1(%esp),%xmm4
        subps  nb114_jzH1(%esp),%xmm5
        movaps %xmm0,nb114_dxOO(%esp)
        movaps %xmm1,nb114_dyOO(%esp)
        movaps %xmm2,nb114_dzOO(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb114_dxH1H1(%esp)
        movaps %xmm4,nb114_dyH1H1(%esp)
        movaps %xmm5,nb114_dzH1H1(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb114_rsqOO(%esp)
        movaps %xmm3,nb114_rsqH1H1(%esp)

        movaps nb114_ixH1(%esp),%xmm0
        movaps nb114_iyH1(%esp),%xmm1
        movaps nb114_izH1(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb114_jxH2(%esp),%xmm0
        subps  nb114_jyH2(%esp),%xmm1
        subps  nb114_jzH2(%esp),%xmm2
        subps  nb114_jxM(%esp),%xmm3
        subps  nb114_jyM(%esp),%xmm4
        subps  nb114_jzM(%esp),%xmm5
        movaps %xmm0,nb114_dxH1H2(%esp)
        movaps %xmm1,nb114_dyH1H2(%esp)
        movaps %xmm2,nb114_dzH1H2(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb114_dxH1M(%esp)
        movaps %xmm4,nb114_dyH1M(%esp)
        movaps %xmm5,nb114_dzH1M(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb114_rsqH1H2(%esp)
        movaps %xmm3,nb114_rsqH1M(%esp)

        movaps nb114_ixH2(%esp),%xmm0
        movaps nb114_iyH2(%esp),%xmm1
        movaps nb114_izH2(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb114_jxH1(%esp),%xmm0
        subps  nb114_jyH1(%esp),%xmm1
        subps  nb114_jzH1(%esp),%xmm2
        subps  nb114_jxH2(%esp),%xmm3
        subps  nb114_jyH2(%esp),%xmm4
        subps  nb114_jzH2(%esp),%xmm5
        movaps %xmm0,nb114_dxH2H1(%esp)
        movaps %xmm1,nb114_dyH2H1(%esp)
        movaps %xmm2,nb114_dzH2H1(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb114_dxH2H2(%esp)
        movaps %xmm4,nb114_dyH2H2(%esp)
        movaps %xmm5,nb114_dzH2H2(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb114_rsqH2H1(%esp)
        movaps %xmm3,nb114_rsqH2H2(%esp)

        movaps nb114_ixH2(%esp),%xmm0
        movaps nb114_iyH2(%esp),%xmm1
        movaps nb114_izH2(%esp),%xmm2
        movaps nb114_ixM(%esp),%xmm3
        movaps nb114_iyM(%esp),%xmm4
        movaps nb114_izM(%esp),%xmm5
        subps  nb114_jxM(%esp),%xmm0
        subps  nb114_jyM(%esp),%xmm1
        subps  nb114_jzM(%esp),%xmm2
        subps  nb114_jxH1(%esp),%xmm3
        subps  nb114_jyH1(%esp),%xmm4
        subps  nb114_jzH1(%esp),%xmm5
        movaps %xmm0,nb114_dxH2M(%esp)
        movaps %xmm1,nb114_dyH2M(%esp)
        movaps %xmm2,nb114_dzH2M(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb114_dxMH1(%esp)
        movaps %xmm4,nb114_dyMH1(%esp)
        movaps %xmm5,nb114_dzMH1(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,nb114_rsqH2M(%esp)
        movaps %xmm4,nb114_rsqMH1(%esp)

        movaps nb114_ixM(%esp),%xmm0
        movaps nb114_iyM(%esp),%xmm1
        movaps nb114_izM(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb114_jxH2(%esp),%xmm0
        subps  nb114_jyH2(%esp),%xmm1
        subps  nb114_jzH2(%esp),%xmm2
        subps  nb114_jxM(%esp),%xmm3
        subps  nb114_jyM(%esp),%xmm4
        subps  nb114_jzM(%esp),%xmm5
        movaps %xmm0,nb114_dxMH2(%esp)
        movaps %xmm1,nb114_dyMH2(%esp)
        movaps %xmm2,nb114_dzMH2(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb114_dxMM(%esp)
        movaps %xmm4,nb114_dyMM(%esp)
        movaps %xmm5,nb114_dzMM(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,nb114_rsqMH2(%esp)
        movaps %xmm4,nb114_rsqMM(%esp)

        ## start by doing reciprocal for OO
        movaps  nb114_rsqOO(%esp),%xmm7
        rcpps   %xmm7,%xmm2
        movaps  nb114_two(%esp),%xmm1
        mulps   %xmm2,%xmm7
        subps   %xmm7,%xmm1
        mulps   %xmm1,%xmm2 ## rinvsq 
        movaps %xmm2,nb114_rinvsqOO(%esp)

        ## next step is invsqrt - do two at a time.
        rsqrtps nb114_rsqH1H1(%esp),%xmm1
        rsqrtps nb114_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb114_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb114_rsqH1H1(%esp),%xmm1
        mulps   nb114_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb114_half(%esp),%xmm3   ## rinvH1H1 
        mulps   nb114_half(%esp),%xmm7   ## rinvH1H2 
        movaps  %xmm3,nb114_rinvH1H1(%esp)
        movaps  %xmm7,nb114_rinvH1H2(%esp)

        rsqrtps nb114_rsqH1M(%esp),%xmm1
        rsqrtps nb114_rsqH2H1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb114_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb114_rsqH1M(%esp),%xmm1
        mulps   nb114_rsqH2H1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb114_half(%esp),%xmm3
        mulps   nb114_half(%esp),%xmm7
        movaps  %xmm3,nb114_rinvH1M(%esp)
        movaps  %xmm7,nb114_rinvH2H1(%esp)

        rsqrtps nb114_rsqH2H2(%esp),%xmm1
        rsqrtps nb114_rsqH2M(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb114_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb114_rsqH2H2(%esp),%xmm1
        mulps   nb114_rsqH2M(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb114_half(%esp),%xmm3
        mulps   nb114_half(%esp),%xmm7
        movaps  %xmm3,nb114_rinvH2H2(%esp)
        movaps  %xmm7,nb114_rinvH2M(%esp)

        rsqrtps nb114_rsqMH1(%esp),%xmm1
        rsqrtps nb114_rsqMH2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb114_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb114_rsqMH1(%esp),%xmm1
        mulps   nb114_rsqMH2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb114_half(%esp),%xmm3
        mulps   nb114_half(%esp),%xmm7
        movaps  %xmm3,nb114_rinvMH1(%esp)
        movaps  %xmm7,nb114_rinvMH2(%esp)

        rsqrtps nb114_rsqMM(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb114_three(%esp),%xmm3
        mulps   nb114_rsqMM(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb114_half(%esp),%xmm3
        movaps  %xmm3,nb114_rinvMM(%esp)

        ## start with OO LJ interaction
        movaps nb114_rinvsqOO(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  %xmm1,%xmm1      ## rinv4
        mulps  %xmm0,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  nb114_c6(%esp),%xmm1
        mulps  nb114_c12(%esp),%xmm2
        movaps %xmm2,%xmm4
        subps  %xmm1,%xmm4
        addps  nb114_Vvdwtot(%esp),%xmm4
        mulps  nb114_six(%esp),%xmm1
        mulps  nb114_twelve(%esp),%xmm2
        movaps %xmm4,nb114_Vvdwtot(%esp)
        subps  %xmm1,%xmm2
        mulps  %xmm2,%xmm0      ## fscal 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb114_dxOO(%esp),%xmm0
        mulps nb114_dyOO(%esp),%xmm1
        mulps nb114_dzOO(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb114_fixO(%esp),%xmm0
        addps nb114_fiyO(%esp),%xmm1
        addps nb114_fizO(%esp),%xmm2
        movaps %xmm3,nb114_fjxO(%esp)
        movaps %xmm4,nb114_fjyO(%esp)
        movaps %xmm5,nb114_fjzO(%esp)
        movaps %xmm0,nb114_fixO(%esp)
        movaps %xmm1,nb114_fiyO(%esp)
        movaps %xmm2,nb114_fizO(%esp)

        ## Coulomb interactions 
        ## start with H1-H1 interaction 
        movaps nb114_rinvH1H1(%esp),%xmm0
        movaps %xmm0,%xmm7
        mulps  %xmm0,%xmm0
        mulps  nb114_qqHH(%esp),%xmm7
        mulps  %xmm7,%xmm0
        addps  nb114_vctot(%esp),%xmm7
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb114_dxH1H1(%esp),%xmm0
        mulps nb114_dyH1H1(%esp),%xmm1
        mulps nb114_dzH1H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb114_fixH1(%esp),%xmm0
        addps nb114_fiyH1(%esp),%xmm1
        addps nb114_fizH1(%esp),%xmm2
        movaps %xmm3,nb114_fjxH1(%esp)
        movaps %xmm4,nb114_fjyH1(%esp)
        movaps %xmm5,nb114_fjzH1(%esp)
        movaps %xmm0,nb114_fixH1(%esp)
        movaps %xmm1,nb114_fiyH1(%esp)
        movaps %xmm2,nb114_fizH1(%esp)

        ## H1-H2 interaction 
        movaps nb114_rinvH1H2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb114_qqHH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fs H1-H2  
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb114_dxH1H2(%esp),%xmm0
        mulps nb114_dyH1H2(%esp),%xmm1
        mulps nb114_dzH1H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb114_fixH1(%esp),%xmm0
        addps nb114_fiyH1(%esp),%xmm1
        addps nb114_fizH1(%esp),%xmm2
        movaps %xmm3,nb114_fjxH2(%esp)
        movaps %xmm4,nb114_fjyH2(%esp)
        movaps %xmm5,nb114_fjzH2(%esp)
        movaps %xmm0,nb114_fixH1(%esp)
        movaps %xmm1,nb114_fiyH1(%esp)
        movaps %xmm2,nb114_fizH1(%esp)

        ## H1-M interaction  
        movaps nb114_rinvH1M(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb114_qqMH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fs H1-M  
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb114_dxH1M(%esp),%xmm0
        mulps nb114_dyH1M(%esp),%xmm1
        mulps nb114_dzH1M(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb114_fixH1(%esp),%xmm0
        addps nb114_fiyH1(%esp),%xmm1
        addps nb114_fizH1(%esp),%xmm2
        movaps %xmm3,nb114_fjxM(%esp)
        movaps %xmm4,nb114_fjyM(%esp)
        movaps %xmm5,nb114_fjzM(%esp)
        movaps %xmm0,nb114_fixH1(%esp)
        movaps %xmm1,nb114_fiyH1(%esp)
        movaps %xmm2,nb114_fizH1(%esp)

        ## H2-H1 interaction 
        movaps nb114_rinvH2H1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb114_qqHH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fs H2-H1 
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2
        movaps nb114_fjxH1(%esp),%xmm3
        movaps nb114_fjyH1(%esp),%xmm4
        movaps nb114_fjzH1(%esp),%xmm5
        mulps nb114_dxH2H1(%esp),%xmm0
        mulps nb114_dyH2H1(%esp),%xmm1
        mulps nb114_dzH2H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb114_fixH2(%esp),%xmm0
        addps nb114_fiyH2(%esp),%xmm1
        addps nb114_fizH2(%esp),%xmm2
        movaps %xmm3,nb114_fjxH1(%esp)
        movaps %xmm4,nb114_fjyH1(%esp)
        movaps %xmm5,nb114_fjzH1(%esp)
        movaps %xmm0,nb114_fixH2(%esp)
        movaps %xmm1,nb114_fiyH2(%esp)
        movaps %xmm2,nb114_fizH2(%esp)

        ## H2-H2 interaction 
        movaps nb114_rinvH2H2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb114_qqHH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fsH2H2
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2
        movaps nb114_fjxH2(%esp),%xmm3
        movaps nb114_fjyH2(%esp),%xmm4
        movaps nb114_fjzH2(%esp),%xmm5
        mulps nb114_dxH2H2(%esp),%xmm0
        mulps nb114_dyH2H2(%esp),%xmm1
        mulps nb114_dzH2H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb114_fixH2(%esp),%xmm0
        addps nb114_fiyH2(%esp),%xmm1
        addps nb114_fizH2(%esp),%xmm2
        movaps %xmm3,nb114_fjxH2(%esp)
        movaps %xmm4,nb114_fjyH2(%esp)
        movaps %xmm5,nb114_fjzH2(%esp)
        movaps %xmm0,nb114_fixH2(%esp)
        movaps %xmm1,nb114_fiyH2(%esp)
        movaps %xmm2,nb114_fizH2(%esp)

        ## H2-M interaction 
        movaps nb114_rinvH2M(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb114_qqMH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fs H2-M  
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2
        movaps nb114_fjxM(%esp),%xmm3
        movaps nb114_fjyM(%esp),%xmm4
        movaps nb114_fjzM(%esp),%xmm5
        mulps nb114_dxH2M(%esp),%xmm0
        mulps nb114_dyH2M(%esp),%xmm1
        mulps nb114_dzH2M(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb114_fixH2(%esp),%xmm0
        addps nb114_fiyH2(%esp),%xmm1
        addps nb114_fizH2(%esp),%xmm2
        movaps %xmm3,nb114_fjxM(%esp)
        movaps %xmm4,nb114_fjyM(%esp)
        movaps %xmm5,nb114_fjzM(%esp)
        movaps %xmm0,nb114_fixH2(%esp)
        movaps %xmm1,nb114_fiyH2(%esp)
        movaps %xmm2,nb114_fizH2(%esp)

        ## M-H1 interaction 
        movaps nb114_rinvMH1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb114_qqMH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fs M-H1 
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2
        movaps nb114_fjxH1(%esp),%xmm3
        movaps nb114_fjyH1(%esp),%xmm4
        movaps nb114_fjzH1(%esp),%xmm5
        mulps nb114_dxMH1(%esp),%xmm0
        mulps nb114_dyMH1(%esp),%xmm1
        mulps nb114_dzMH1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb114_fixM(%esp),%xmm0
        addps nb114_fiyM(%esp),%xmm1
        addps nb114_fizM(%esp),%xmm2
        movaps %xmm3,nb114_fjxH1(%esp)
        movaps %xmm4,nb114_fjyH1(%esp)
        movaps %xmm5,nb114_fjzH1(%esp)
        movaps %xmm0,nb114_fixM(%esp)
        movaps %xmm1,nb114_fiyM(%esp)
        movaps %xmm2,nb114_fizM(%esp)

        ## M-H2 interaction 
        movaps nb114_rinvMH2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb114_qqMH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fs M-H2 
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2
        movaps nb114_fjxH2(%esp),%xmm3
        movaps nb114_fjyH2(%esp),%xmm4
        movaps nb114_fjzH2(%esp),%xmm5
        mulps nb114_dxMH2(%esp),%xmm0
        mulps nb114_dyMH2(%esp),%xmm1
        mulps nb114_dzMH2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb114_fixM(%esp),%xmm0
        addps nb114_fiyM(%esp),%xmm1
        addps nb114_fizM(%esp),%xmm2
        movaps %xmm3,nb114_fjxH2(%esp)
        movaps %xmm4,nb114_fjyH2(%esp)
        movaps %xmm5,nb114_fjzH2(%esp)
        movaps %xmm0,nb114_fixM(%esp)
        movaps %xmm1,nb114_fiyM(%esp)
        movaps %xmm2,nb114_fizM(%esp)

        ## M-M interaction 
        movaps nb114_rinvMM(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb114_qqMM(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fs M-M 
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm7,nb114_vctot(%esp)
        movaps %xmm0,%xmm2
        movaps nb114_fjxM(%esp),%xmm3
        movaps nb114_fjyM(%esp),%xmm4
        movaps nb114_fjzM(%esp),%xmm5
        mulps nb114_dxMM(%esp),%xmm0
        mulps nb114_dyMM(%esp),%xmm1
        mulps nb114_dzMM(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb114_fixM(%esp),%xmm0
        addps nb114_fiyM(%esp),%xmm1
        addps nb114_fizM(%esp),%xmm2
        movaps %xmm3,nb114_fjxM(%esp)
        movaps %xmm4,nb114_fjyM(%esp)
        movaps %xmm5,nb114_fjzM(%esp)
        movaps %xmm0,nb114_fixM(%esp)
        movaps %xmm1,nb114_fiyM(%esp)
        movaps %xmm2,nb114_fizM(%esp)

        ## Did all interactions - now update j forces 
        movl nb114_faction(%ebp),%edi
        ## 4 j waters with four atoms each.
        ## step 1 : transpose fjxO, fjyO, fjzO, fjxH1
        movaps nb114_fjxO(%esp),%xmm0
        movaps nb114_fjyO(%esp),%xmm1
        movaps nb114_fjzO(%esp),%xmm2
        movaps nb114_fjxH1(%esp),%xmm3
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
        movaps nb114_fjyH1(%esp),%xmm0
        movaps nb114_fjzH1(%esp),%xmm1
        movaps nb114_fjxH2(%esp),%xmm2
        movaps nb114_fjyH2(%esp),%xmm3
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
        movaps nb114_fjzH2(%esp),%xmm0
        movaps nb114_fjxM(%esp),%xmm1
        movaps nb114_fjyM(%esp),%xmm2
        movaps nb114_fjzM(%esp),%xmm3

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
        subl $4,nb114_innerk(%esp)
        jl    _nb_kernel114_ia32_sse.nb114_single_check
        jmp   _nb_kernel114_ia32_sse.nb114_unroll_loop
_nb_kernel114_ia32_sse.nb114_single_check: 
        addl $4,nb114_innerk(%esp)
        jnz   _nb_kernel114_ia32_sse.nb114_single_loop
        jmp   _nb_kernel114_ia32_sse.nb114_updateouterdata
_nb_kernel114_ia32_sse.nb114_single_loop: 
        movl  nb114_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb114_innerjjnr(%esp)

        movl nb114_pos(%ebp),%esi
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
        movaps %xmm6,nb114_jxO(%esp)
        movaps %xmm3,nb114_jyO(%esp)
        movaps %xmm1,nb114_jzO(%esp)

        ## do O and M in parallel
        movaps nb114_ixO(%esp),%xmm0
        movaps nb114_iyO(%esp),%xmm1
        movaps nb114_izO(%esp),%xmm2
        movaps nb114_ixM(%esp),%xmm3
        movaps nb114_iyM(%esp),%xmm4
        movaps nb114_izM(%esp),%xmm5
        subps  nb114_jxO(%esp),%xmm0
        subps  nb114_jyO(%esp),%xmm1
        subps  nb114_jzO(%esp),%xmm2
        subps  nb114_jxO(%esp),%xmm3
        subps  nb114_jyO(%esp),%xmm4
        subps  nb114_jzO(%esp),%xmm5

        movaps %xmm0,nb114_dxOO(%esp)
        movaps %xmm1,nb114_dyOO(%esp)
        movaps %xmm2,nb114_dzOO(%esp)
        movaps %xmm3,nb114_dxMM(%esp)
        movaps %xmm4,nb114_dyMM(%esp)
        movaps %xmm5,nb114_dzMM(%esp)

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
        ## Save M data 
        movaps %xmm4,nb114_rsqMM(%esp)

        ## do 1/x for O and 1/sqrt(x) for M
        rcpss  %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movss  nb114_two(%esp),%xmm2
        movaps  %xmm5,%xmm6
        mulss  %xmm1,%xmm0
        mulps   %xmm5,%xmm5
        subss  %xmm0,%xmm2
        movaps  nb114_three(%esp),%xmm7
        mulss  %xmm1,%xmm2      ## constant 1/r2

        mulps   %xmm4,%xmm5
        movss  %xmm2,%xmm0
        subps   %xmm5,%xmm7
        mulss  %xmm2,%xmm2
        mulps   %xmm6,%xmm7
        mulss  %xmm0,%xmm2      ## constant 1/r6
        mulps   nb114_half(%esp),%xmm7   ## rinv iH1 - j water 
        movss  %xmm2,%xmm1
        movaps %xmm7,nb114_rinvMM(%esp)

        mulss  %xmm2,%xmm2      ## constant 1/r12
        mulss  nb114_c6(%esp),%xmm1
        mulss  nb114_c12(%esp),%xmm2
        movss  %xmm2,%xmm3
        subss  %xmm1,%xmm3
        addss  nb114_Vvdwtot(%esp),%xmm3
        movss  %xmm3,nb114_Vvdwtot(%esp)
        mulss  nb114_six(%esp),%xmm1
        mulss  nb114_twelve(%esp),%xmm2
        subss  %xmm1,%xmm2
        mulss  %xmm2,%xmm0      ## fscal
        movss  %xmm0,%xmm1
        movss  %xmm0,%xmm2
        mulss  nb114_dxOO(%esp),%xmm0
        mulss  nb114_dyOO(%esp),%xmm1
        mulss  nb114_dzOO(%esp),%xmm2
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        subss   %xmm0,%xmm3
        subss   %xmm1,%xmm4
        subss   %xmm2,%xmm5
        movaps  %xmm3,nb114_fjxO(%esp)
        movaps  %xmm4,nb114_fjyO(%esp)
        movaps  %xmm5,nb114_fjzO(%esp)
        addss   nb114_fixO(%esp),%xmm0
        addss   nb114_fiyO(%esp),%xmm1
        addss   nb114_fizO(%esp),%xmm2
        movss  %xmm0,nb114_fixO(%esp)
        movss  %xmm1,nb114_fiyO(%esp)
        movss  %xmm2,nb114_fizO(%esp)

        ## do  M coulomb interaction
        movaps nb114_rinvMM(%esp),%xmm0
        movaps %xmm0,%xmm4      ## xmm4=rinv
        mulps  %xmm0,%xmm0      ## xmm0=rinvsq 

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb114_qqMH(%esp),%xmm3
        movhps  nb114_qqMM(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## constant 11000001 

        mulps  %xmm3,%xmm4      ## xmm4=voul    
        mulps  %xmm4,%xmm0
        addps  nb114_vctot(%esp),%xmm4
        movaps %xmm4,nb114_vctot(%esp)

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        mulps   nb114_dxMM(%esp),%xmm0
        mulps   nb114_dyMM(%esp),%xmm1
        mulps   nb114_dzMM(%esp),%xmm2
        ## update forces M - j water 
        movaps  nb114_fjxO(%esp),%xmm3
        movaps  nb114_fjyO(%esp),%xmm4
        movaps  nb114_fjzO(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb114_fjxO(%esp)
        movaps  %xmm4,nb114_fjyO(%esp)
        movaps  %xmm5,nb114_fjzO(%esp)
        addps   nb114_fixM(%esp),%xmm0
        addps   nb114_fiyM(%esp),%xmm1
        addps   nb114_fizM(%esp),%xmm2
        movaps  %xmm0,nb114_fixM(%esp)
        movaps  %xmm1,nb114_fiyM(%esp)
        movaps  %xmm2,nb114_fizM(%esp)

        ## i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb114_ixH1(%esp),%xmm0
        movaps  nb114_iyH1(%esp),%xmm1
        movaps  nb114_izH1(%esp),%xmm2
        movaps  nb114_ixH2(%esp),%xmm3
        movaps  nb114_iyH2(%esp),%xmm4
        movaps  nb114_izH2(%esp),%xmm5
        subps   nb114_jxO(%esp),%xmm0
        subps   nb114_jyO(%esp),%xmm1
        subps   nb114_jzO(%esp),%xmm2
        subps   nb114_jxO(%esp),%xmm3
        subps   nb114_jyO(%esp),%xmm4
        subps   nb114_jzO(%esp),%xmm5
        movaps %xmm0,nb114_dxH1H1(%esp)
        movaps %xmm1,nb114_dyH1H1(%esp)
        movaps %xmm2,nb114_dzH1H1(%esp)
        movaps %xmm3,nb114_dxH2H2(%esp)
        movaps %xmm4,nb114_dyH2H2(%esp)
        movaps %xmm5,nb114_dzH2H2(%esp)
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
        movaps  %xmm0,nb114_rsqH1H1(%esp)
        movaps  %xmm4,nb114_rsqH2H2(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb114_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb114_half(%esp),%xmm3   ## rinvH1H1
        mulps   nb114_half(%esp),%xmm7   ## rinvH2H2
        movaps  %xmm3,nb114_rinvH1H1(%esp)
        movaps  %xmm7,nb114_rinvH2H2(%esp)

        ## Do H1 coulomb interaction
        movaps nb114_rinvH1H1(%esp),%xmm0
        movaps %xmm0,%xmm4      ## xmm4=rinv 
        mulps  %xmm0,%xmm0      ## xmm0=rinvsq 

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb114_qqHH(%esp),%xmm3
        movhps  nb114_qqMH(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## constant 11000001 

        mulps  %xmm3,%xmm4      ## xmm4=voul
        mulps  %xmm4,%xmm0
        addps  nb114_vctot(%esp),%xmm4
        movaps %xmm4,nb114_vctot(%esp)


        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        mulps   nb114_dxH1H1(%esp),%xmm0
        mulps   nb114_dyH1H1(%esp),%xmm1
        mulps   nb114_dzH1H1(%esp),%xmm2
        ## update forces H1 - j water 
        movaps  nb114_fjxO(%esp),%xmm3
        movaps  nb114_fjyO(%esp),%xmm4
        movaps  nb114_fjzO(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb114_fjxO(%esp)
        movaps  %xmm4,nb114_fjyO(%esp)
        movaps  %xmm5,nb114_fjzO(%esp)
        addps   nb114_fixH1(%esp),%xmm0
        addps   nb114_fiyH1(%esp),%xmm1
        addps   nb114_fizH1(%esp),%xmm2
        movaps  %xmm0,nb114_fixH1(%esp)
        movaps  %xmm1,nb114_fiyH1(%esp)
        movaps  %xmm2,nb114_fizH1(%esp)

        ## H2 Coulomb
        movaps nb114_rinvH2H2(%esp),%xmm0
        movaps %xmm0,%xmm4      ## xmm4=rinv 
        mulps  %xmm0,%xmm0      ## xmm0=rinvsq 

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb114_qqHH(%esp),%xmm3
        movhps  nb114_qqMH(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## constant 11000001 

        mulps  %xmm3,%xmm4      ## xmm4=voul
        mulps  %xmm4,%xmm0
        addps  nb114_vctot(%esp),%xmm4   ## local vctot summation variable
        movaps %xmm4,nb114_vctot(%esp)

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        mulps   nb114_dxH2H2(%esp),%xmm0
        mulps   nb114_dyH2H2(%esp),%xmm1
        mulps   nb114_dzH2H2(%esp),%xmm2
        ## update forces H2 - j water 
        movaps  nb114_fjxO(%esp),%xmm3
        movaps  nb114_fjyO(%esp),%xmm4
        movaps  nb114_fjzO(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb114_fjxO(%esp)
        movaps  %xmm4,nb114_fjyO(%esp)
        movaps  %xmm5,nb114_fjzO(%esp)
        addps   nb114_fixH2(%esp),%xmm0
        addps   nb114_fiyH2(%esp),%xmm1
        addps   nb114_fizH2(%esp),%xmm2
        movaps  %xmm0,nb114_fixH2(%esp)
        movaps  %xmm1,nb114_fiyH2(%esp)
        movaps  %xmm2,nb114_fizH2(%esp)

        movl    nb114_faction(%ebp),%esi
        ## update j water forces from local variables.
        ## transpose back first
        movaps  nb114_fjxO(%esp),%xmm0   ## Ox H1x H2x Mx 
        movaps  nb114_fjyO(%esp),%xmm1   ## Oy H1y H2y My
        movaps  nb114_fjzO(%esp),%xmm2   ## Oz H1z H2z Mz

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

        decl nb114_innerk(%esp)
        jz    _nb_kernel114_ia32_sse.nb114_updateouterdata
        jmp   _nb_kernel114_ia32_sse.nb114_single_loop
_nb_kernel114_ia32_sse.nb114_updateouterdata: 
        movl  nb114_ii3(%esp),%ecx
        movl  nb114_faction(%ebp),%edi
        movl  nb114_fshift(%ebp),%esi
        movl  nb114_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb114_fixO(%esp),%xmm0
        movaps nb114_fiyO(%esp),%xmm1
        movaps nb114_fizO(%esp),%xmm2

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
        movaps nb114_fixH1(%esp),%xmm0
        movaps nb114_fiyH1(%esp),%xmm1
        movaps nb114_fizH1(%esp),%xmm2

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
        movaps nb114_fixH2(%esp),%xmm0
        movaps nb114_fiyH2(%esp),%xmm1
        movaps nb114_fizH2(%esp),%xmm2

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
        movaps nb114_fixM(%esp),%xmm0
        movaps nb114_fiyM(%esp),%xmm1
        movaps nb114_fizM(%esp),%xmm2

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
        movl nb114_n(%esp),%esi
        ## get group index for i particle 
        movl  nb114_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb114_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb114_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb114_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb114_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb114_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel114_ia32_sse.nb114_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb114_n(%esp)
        jmp _nb_kernel114_ia32_sse.nb114_outer
_nb_kernel114_ia32_sse.nb114_outerend: 
        ## check if more outer neighborlists remain
        movl  nb114_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel114_ia32_sse.nb114_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel114_ia32_sse.nb114_threadloop
_nb_kernel114_ia32_sse.nb114_end: 
        emms

        movl nb114_nouter(%esp),%eax
        movl nb114_ninner(%esp),%ebx
        movl nb114_outeriter(%ebp),%ecx
        movl nb114_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb114_salign(%esp),%eax
        addl %eax,%esp
        addl $1816,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret





.globl nb_kernel114nf_ia32_sse
.globl _nb_kernel114nf_ia32_sse
nb_kernel114nf_ia32_sse:        
_nb_kernel114nf_ia32_sse:       
.set nb114nf_p_nri, 8
.set nb114nf_iinr, 12
.set nb114nf_jindex, 16
.set nb114nf_jjnr, 20
.set nb114nf_shift, 24
.set nb114nf_shiftvec, 28
.set nb114nf_fshift, 32
.set nb114nf_gid, 36
.set nb114nf_pos, 40
.set nb114nf_faction, 44
.set nb114nf_charge, 48
.set nb114nf_p_facel, 52
.set nb114nf_p_krf, 56
.set nb114nf_p_crf, 60
.set nb114nf_Vc, 64
.set nb114nf_type, 68
.set nb114nf_p_ntype, 72
.set nb114nf_vdwparam, 76
.set nb114nf_Vvdw, 80
.set nb114nf_p_tabscale, 84
.set nb114nf_VFtab, 88
.set nb114nf_invsqrta, 92
.set nb114nf_dvda, 96
.set nb114nf_p_gbtabscale, 100
.set nb114nf_GBtab, 104
.set nb114nf_p_nthreads, 108
.set nb114nf_count, 112
.set nb114nf_mtx, 116
.set nb114nf_outeriter, 120
.set nb114nf_inneriter, 124
.set nb114nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb114nf_ixO, 0
.set nb114nf_iyO, 16
.set nb114nf_izO, 32
.set nb114nf_ixH1, 48
.set nb114nf_iyH1, 64
.set nb114nf_izH1, 80
.set nb114nf_ixH2, 96
.set nb114nf_iyH2, 112
.set nb114nf_izH2, 128
.set nb114nf_ixM, 144
.set nb114nf_iyM, 160
.set nb114nf_izM, 176
.set nb114nf_jxO, 192
.set nb114nf_jyO, 208
.set nb114nf_jzO, 224
.set nb114nf_jxH1, 240
.set nb114nf_jyH1, 256
.set nb114nf_jzH1, 272
.set nb114nf_jxH2, 288
.set nb114nf_jyH2, 304
.set nb114nf_jzH2, 320
.set nb114nf_jxM, 336
.set nb114nf_jyM, 352
.set nb114nf_jzM, 368
.set nb114nf_qqMM, 384
.set nb114nf_qqMH, 400
.set nb114nf_qqHH, 416
.set nb114nf_two, 432
.set nb114nf_c6, 448
.set nb114nf_c12, 464
.set nb114nf_vctot, 480
.set nb114nf_Vvdwtot, 496
.set nb114nf_half, 512
.set nb114nf_three, 528
.set nb114nf_rsqOO, 544
.set nb114nf_rsqH1H1, 560
.set nb114nf_rsqH1H2, 576
.set nb114nf_rsqH1M, 592
.set nb114nf_rsqH2H1, 608
.set nb114nf_rsqH2H2, 624
.set nb114nf_rsqH2M, 640
.set nb114nf_rsqMH1, 656
.set nb114nf_rsqMH2, 672
.set nb114nf_rsqMM, 688
.set nb114nf_rinvsqOO, 704
.set nb114nf_rinvH1H1, 720
.set nb114nf_rinvH1H2, 736
.set nb114nf_rinvH1M, 752
.set nb114nf_rinvH2H1, 768
.set nb114nf_rinvH2H2, 784
.set nb114nf_rinvH2M, 800
.set nb114nf_rinvMH1, 816
.set nb114nf_rinvMH2, 832
.set nb114nf_rinvMM, 848
.set nb114nf_is3, 864
.set nb114nf_ii3, 868
.set nb114nf_innerjjnr, 872
.set nb114nf_innerk, 876
.set nb114nf_n, 880
.set nb114nf_nn1, 884
.set nb114nf_nri, 888
.set nb114nf_nouter, 892
.set nb114nf_ninner, 896
.set nb114nf_salign, 900
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $904,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb114nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb114nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb114nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb114nf_nouter(%esp)
        movl %eax,nb114nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb114nf_half(%esp)
        movss nb114nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb114nf_half(%esp)
        movaps %xmm2,nb114nf_two(%esp)
        movaps %xmm3,nb114nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb114nf_iinr(%ebp),%ecx     ## ecx = pointer into iinr[]
        movl  (%ecx),%ebx               ## ebx =ii 

        movl  nb114nf_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm5
        movss 12(%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movl nb114nf_p_facel(%ebp),%esi
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
        movaps %xmm3,nb114nf_qqMM(%esp)
        movaps %xmm4,nb114nf_qqMH(%esp)
        movaps %xmm5,nb114nf_qqHH(%esp)

        xorps %xmm0,%xmm0
        movl  nb114nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb114nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb114nf_vdwparam(%ebp),%eax
        movlps (%eax,%edx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $0x55,%xmm1,%xmm1
        movaps %xmm0,nb114nf_c6(%esp)
        movaps %xmm1,nb114nf_c12(%esp)

_nb_kernel114nf_ia32_sse.nb114nf_threadloop: 
        movl  nb114nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel114nf_ia32_sse.nb114nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel114nf_ia32_sse.nb114nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb114nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb114nf_n(%esp)
        movl %ebx,nb114nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel114nf_ia32_sse.nb114nf_outerstart
        jmp _nb_kernel114nf_ia32_sse.nb114nf_end

_nb_kernel114nf_ia32_sse.nb114nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb114nf_nouter(%esp),%ebx
        movl %ebx,nb114nf_nouter(%esp)

_nb_kernel114nf_ia32_sse.nb114nf_outer: 
        movl  nb114nf_shift(%ebp),%eax          ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb114nf_is3(%esp)            ## store is3 

        movl  nb114nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb114nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb114nf_pos(%ebp),%eax    ## eax = base of pos[]  
        movl  %ebx,nb114nf_ii3(%esp)

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
        movaps %xmm3,nb114nf_ixO(%esp)
        movaps %xmm4,nb114nf_iyO(%esp)
        movaps %xmm5,nb114nf_izO(%esp)
        movaps %xmm6,nb114nf_ixH1(%esp)
        movaps %xmm7,nb114nf_iyH1(%esp)

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
        movaps %xmm6,nb114nf_izH1(%esp)
        movaps %xmm0,nb114nf_ixH2(%esp)
        movaps %xmm1,nb114nf_iyH2(%esp)
        movaps %xmm2,nb114nf_izH2(%esp)
        movaps %xmm3,nb114nf_ixM(%esp)
        movaps %xmm4,nb114nf_iyM(%esp)
        movaps %xmm5,nb114nf_izM(%esp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb114nf_vctot(%esp)
        movaps %xmm4,nb114nf_Vvdwtot(%esp)

        movl  nb114nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb114nf_pos(%ebp),%esi
        movl  nb114nf_faction(%ebp),%edi
        movl  nb114nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb114nf_innerjjnr(%esp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb114nf_ninner(%esp),%ecx
        movl  %ecx,nb114nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb114nf_innerk(%esp)         ## number of innerloop atoms 
        jge   _nb_kernel114nf_ia32_sse.nb114nf_unroll_loop
        jmp   _nb_kernel114nf_ia32_sse.nb114nf_single_check
_nb_kernel114nf_ia32_sse.nb114nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb114nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb114nf_innerjjnr(%esp)             ## advance pointer (unroll 4) 

        movl nb114nf_pos(%ebp),%esi     ## base of pos[] 

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
        movaps %xmm0,nb114nf_jxO(%esp)
        movaps %xmm1,nb114nf_jyO(%esp)
        movaps %xmm2,nb114nf_jzO(%esp)
        movaps %xmm3,nb114nf_jxH1(%esp)

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
        movaps %xmm0,nb114nf_jyH1(%esp)
        movaps %xmm1,nb114nf_jzH1(%esp)
        movaps %xmm2,nb114nf_jxH2(%esp)
        movaps %xmm3,nb114nf_jyH2(%esp)

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
        movaps %xmm0,nb114nf_jzH2(%esp)
        movaps %xmm1,nb114nf_jxM(%esp)
        movaps %xmm2,nb114nf_jyM(%esp)
        movaps %xmm3,nb114nf_jzM(%esp)

        ## start calculating pairwise distances
        movaps nb114nf_ixO(%esp),%xmm0
        movaps nb114nf_iyO(%esp),%xmm1
        movaps nb114nf_izO(%esp),%xmm2
        movaps nb114nf_ixH1(%esp),%xmm3
        movaps nb114nf_iyH1(%esp),%xmm4
        movaps nb114nf_izH1(%esp),%xmm5
        subps  nb114nf_jxO(%esp),%xmm0
        subps  nb114nf_jyO(%esp),%xmm1
        subps  nb114nf_jzO(%esp),%xmm2
        subps  nb114nf_jxH1(%esp),%xmm3
        subps  nb114nf_jyH1(%esp),%xmm4
        subps  nb114nf_jzH1(%esp),%xmm5
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
        movaps %xmm0,nb114nf_rsqOO(%esp)
        movaps %xmm3,nb114nf_rsqH1H1(%esp)

        movaps nb114nf_ixH1(%esp),%xmm0
        movaps nb114nf_iyH1(%esp),%xmm1
        movaps nb114nf_izH1(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb114nf_jxH2(%esp),%xmm0
        subps  nb114nf_jyH2(%esp),%xmm1
        subps  nb114nf_jzH2(%esp),%xmm2
        subps  nb114nf_jxM(%esp),%xmm3
        subps  nb114nf_jyM(%esp),%xmm4
        subps  nb114nf_jzM(%esp),%xmm5
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
        movaps %xmm0,nb114nf_rsqH1H2(%esp)
        movaps %xmm3,nb114nf_rsqH1M(%esp)

        movaps nb114nf_ixH2(%esp),%xmm0
        movaps nb114nf_iyH2(%esp),%xmm1
        movaps nb114nf_izH2(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb114nf_jxH1(%esp),%xmm0
        subps  nb114nf_jyH1(%esp),%xmm1
        subps  nb114nf_jzH1(%esp),%xmm2
        subps  nb114nf_jxH2(%esp),%xmm3
        subps  nb114nf_jyH2(%esp),%xmm4
        subps  nb114nf_jzH2(%esp),%xmm5
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
        movaps %xmm0,nb114nf_rsqH2H1(%esp)
        movaps %xmm3,nb114nf_rsqH2H2(%esp)

        movaps nb114nf_ixH2(%esp),%xmm0
        movaps nb114nf_iyH2(%esp),%xmm1
        movaps nb114nf_izH2(%esp),%xmm2
        movaps nb114nf_ixM(%esp),%xmm3
        movaps nb114nf_iyM(%esp),%xmm4
        movaps nb114nf_izM(%esp),%xmm5
        subps  nb114nf_jxM(%esp),%xmm0
        subps  nb114nf_jyM(%esp),%xmm1
        subps  nb114nf_jzM(%esp),%xmm2
        subps  nb114nf_jxH1(%esp),%xmm3
        subps  nb114nf_jyH1(%esp),%xmm4
        subps  nb114nf_jzH1(%esp),%xmm5
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
        movaps %xmm0,nb114nf_rsqH2M(%esp)
        movaps %xmm4,nb114nf_rsqMH1(%esp)

        movaps nb114nf_ixM(%esp),%xmm0
        movaps nb114nf_iyM(%esp),%xmm1
        movaps nb114nf_izM(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb114nf_jxH2(%esp),%xmm0
        subps  nb114nf_jyH2(%esp),%xmm1
        subps  nb114nf_jzH2(%esp),%xmm2
        subps  nb114nf_jxM(%esp),%xmm3
        subps  nb114nf_jyM(%esp),%xmm4
        subps  nb114nf_jzM(%esp),%xmm5
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
        movaps %xmm0,nb114nf_rsqMH2(%esp)
        movaps %xmm4,nb114nf_rsqMM(%esp)

        ## start by doing reciprocal for OO
        movaps  nb114nf_rsqOO(%esp),%xmm7
        rcpps   %xmm7,%xmm2
        movaps  nb114nf_two(%esp),%xmm1
        mulps   %xmm2,%xmm7
        subps   %xmm7,%xmm1
        mulps   %xmm1,%xmm2 ## rinvsq 
        movaps %xmm2,nb114nf_rinvsqOO(%esp)

        ## next step is invsqrt - do two at a time.
        rsqrtps nb114nf_rsqH1H1(%esp),%xmm1
        rsqrtps nb114nf_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb114nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb114nf_rsqH1H1(%esp),%xmm1
        mulps   nb114nf_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb114nf_half(%esp),%xmm3   ## rinvH1H1 
        mulps   nb114nf_half(%esp),%xmm7   ## rinvH1H2 
        movaps  %xmm3,nb114nf_rinvH1H1(%esp)
        movaps  %xmm7,nb114nf_rinvH1H2(%esp)

        rsqrtps nb114nf_rsqH1M(%esp),%xmm1
        rsqrtps nb114nf_rsqH2H1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb114nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb114nf_rsqH1M(%esp),%xmm1
        mulps   nb114nf_rsqH2H1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb114nf_half(%esp),%xmm3
        mulps   nb114nf_half(%esp),%xmm7
        movaps  %xmm3,nb114nf_rinvH1M(%esp)
        movaps  %xmm7,nb114nf_rinvH2H1(%esp)

        rsqrtps nb114nf_rsqH2H2(%esp),%xmm1
        rsqrtps nb114nf_rsqH2M(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb114nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb114nf_rsqH2H2(%esp),%xmm1
        mulps   nb114nf_rsqH2M(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb114nf_half(%esp),%xmm3
        mulps   nb114nf_half(%esp),%xmm7
        movaps  %xmm3,nb114nf_rinvH2H2(%esp)
        movaps  %xmm7,nb114nf_rinvH2M(%esp)

        rsqrtps nb114nf_rsqMH1(%esp),%xmm1
        rsqrtps nb114nf_rsqMH2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb114nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb114nf_rsqMH1(%esp),%xmm1
        mulps   nb114nf_rsqMH2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb114nf_half(%esp),%xmm3
        mulps   nb114nf_half(%esp),%xmm7
        movaps  %xmm3,nb114nf_rinvMH1(%esp)
        movaps  %xmm7,nb114nf_rinvMH2(%esp)

        rsqrtps nb114nf_rsqMM(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb114nf_three(%esp),%xmm3
        mulps   nb114nf_rsqMM(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb114nf_half(%esp),%xmm3
        movaps  %xmm3,nb114nf_rinvMM(%esp)

        ## start with OO LJ interaction
        movaps nb114nf_rinvsqOO(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  %xmm1,%xmm1      ## rinv4
        mulps  %xmm0,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  nb114nf_c6(%esp),%xmm1
        mulps  nb114nf_c12(%esp),%xmm2
        movaps %xmm2,%xmm4
        subps  %xmm1,%xmm4
        addps  nb114nf_Vvdwtot(%esp),%xmm4
        movaps %xmm4,nb114nf_Vvdwtot(%esp)

        ## Coulomb interactions 
        ## all H-H interactions
        movaps nb114nf_rinvH1H1(%esp),%xmm0
        addps  nb114nf_rinvH1H2(%esp),%xmm0
        addps  nb114nf_rinvH2H1(%esp),%xmm0
        addps  nb114nf_rinvH2H2(%esp),%xmm0
        mulps  nb114nf_qqHH(%esp),%xmm0
        ## all M-H interactions
        movaps nb114nf_rinvH1M(%esp),%xmm1
        addps  nb114nf_rinvH2M(%esp),%xmm1
        addps  nb114nf_rinvMH1(%esp),%xmm1
        addps  nb114nf_rinvMH2(%esp),%xmm1
        mulps  nb114nf_qqMH(%esp),%xmm1
        ## The M-M interaction
        movaps nb114nf_rinvMM(%esp),%xmm2
        mulps  nb114nf_qqMM(%esp),%xmm2
        addps  %xmm1,%xmm0
        addps  nb114nf_vctot(%esp),%xmm2
        addps  %xmm2,%xmm0
        movaps %xmm0,nb114nf_vctot(%esp)
        ## should we do one more iteration? 
        subl $4,nb114nf_innerk(%esp)
        jl    _nb_kernel114nf_ia32_sse.nb114nf_single_check
        jmp   _nb_kernel114nf_ia32_sse.nb114nf_unroll_loop
_nb_kernel114nf_ia32_sse.nb114nf_single_check: 
        addl $4,nb114nf_innerk(%esp)
        jnz   _nb_kernel114nf_ia32_sse.nb114nf_single_loop
        jmp   _nb_kernel114nf_ia32_sse.nb114nf_updateouterdata
_nb_kernel114nf_ia32_sse.nb114nf_single_loop: 
        movl  nb114nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb114nf_innerjjnr(%esp)

        movl nb114nf_pos(%ebp),%esi
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
        movaps %xmm6,nb114nf_jxO(%esp)
        movaps %xmm3,nb114nf_jyO(%esp)
        movaps %xmm1,nb114nf_jzO(%esp)

        ## do O and M in parallel
        movaps nb114nf_ixO(%esp),%xmm0
        movaps nb114nf_iyO(%esp),%xmm1
        movaps nb114nf_izO(%esp),%xmm2
        movaps nb114nf_ixM(%esp),%xmm3
        movaps nb114nf_iyM(%esp),%xmm4
        movaps nb114nf_izM(%esp),%xmm5
        subps  nb114nf_jxO(%esp),%xmm0
        subps  nb114nf_jyO(%esp),%xmm1
        subps  nb114nf_jzO(%esp),%xmm2
        subps  nb114nf_jxO(%esp),%xmm3
        subps  nb114nf_jyO(%esp),%xmm4
        subps  nb114nf_jzO(%esp),%xmm5

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
        ## Save M data 
        movaps %xmm4,nb114nf_rsqMM(%esp)

        ## do 1/x for O and 1/sqrt(x) for M
        rcpss  %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movss  nb114nf_two(%esp),%xmm2
        movaps  %xmm5,%xmm6
        mulss  %xmm1,%xmm0
        mulps   %xmm5,%xmm5
        subss  %xmm0,%xmm2
        movaps  nb114nf_three(%esp),%xmm7
        mulss  %xmm1,%xmm2      ## constant 1/r2

        mulps   %xmm4,%xmm5
        movss  %xmm2,%xmm0
        subps   %xmm5,%xmm7
        mulss  %xmm2,%xmm2
        mulps   %xmm6,%xmm7
        mulss  %xmm0,%xmm2      ## constant 1/r6
        mulps   nb114nf_half(%esp),%xmm7   ## rinv iH1 - j water 
        movss  %xmm2,%xmm1
        movaps %xmm7,nb114nf_rinvMM(%esp)

        mulss  %xmm2,%xmm2      ## constant 1/r12
        mulss  nb114nf_c6(%esp),%xmm1
        mulss  nb114nf_c12(%esp),%xmm2
        movss  %xmm2,%xmm3
        subss  %xmm1,%xmm3
        addss  nb114nf_Vvdwtot(%esp),%xmm3
        movss  %xmm3,nb114nf_Vvdwtot(%esp)

        ## do  M coulomb interaction
        movaps nb114nf_rinvMM(%esp),%xmm0
        movaps %xmm0,%xmm4      ## xmm4=rinv

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb114nf_qqMH(%esp),%xmm3
        movhps  nb114nf_qqMM(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## constant 11000001 

        mulps  %xmm3,%xmm4      ## xmm4=voul    
        addps  nb114nf_vctot(%esp),%xmm4
        movaps %xmm4,nb114nf_vctot(%esp)

        ## i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb114nf_ixH1(%esp),%xmm0
        movaps  nb114nf_iyH1(%esp),%xmm1
        movaps  nb114nf_izH1(%esp),%xmm2
        movaps  nb114nf_ixH2(%esp),%xmm3
        movaps  nb114nf_iyH2(%esp),%xmm4
        movaps  nb114nf_izH2(%esp),%xmm5
        subps   nb114nf_jxO(%esp),%xmm0
        subps   nb114nf_jyO(%esp),%xmm1
        subps   nb114nf_jzO(%esp),%xmm2
        subps   nb114nf_jxO(%esp),%xmm3
        subps   nb114nf_jyO(%esp),%xmm4
        subps   nb114nf_jzO(%esp),%xmm5
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
        movaps  %xmm0,nb114nf_rsqH1H1(%esp)
        movaps  %xmm4,nb114nf_rsqH2H2(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb114nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb114nf_half(%esp),%xmm3   ## rinvH1H1
        mulps   nb114nf_half(%esp),%xmm7   ## rinvH2H2
        movaps  %xmm3,nb114nf_rinvH1H1(%esp)
        movaps  %xmm7,nb114nf_rinvH2H2(%esp)

        ## Do H1 & H2 coulomb interaction
        movaps nb114nf_rinvH1H1(%esp),%xmm0
        addps  nb114nf_rinvH2H2(%esp),%xmm0

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb114nf_qqHH(%esp),%xmm3
        movhps  nb114nf_qqMH(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## constant 11000001 

        mulps  %xmm3,%xmm0      ## xmm4=voul
        addps  nb114nf_vctot(%esp),%xmm0
        movaps %xmm0,nb114nf_vctot(%esp)

        decl nb114nf_innerk(%esp)
        jz    _nb_kernel114nf_ia32_sse.nb114nf_updateouterdata
        jmp   _nb_kernel114nf_ia32_sse.nb114nf_single_loop
_nb_kernel114nf_ia32_sse.nb114nf_updateouterdata: 
        ## get n from stack
        movl nb114nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb114nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb114nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb114nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb114nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb114nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb114nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel114nf_ia32_sse.nb114nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb114nf_n(%esp)
        jmp _nb_kernel114nf_ia32_sse.nb114nf_outer
_nb_kernel114nf_ia32_sse.nb114nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb114nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel114nf_ia32_sse.nb114nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel114nf_ia32_sse.nb114nf_threadloop
_nb_kernel114nf_ia32_sse.nb114nf_end: 
        emms

        movl nb114nf_nouter(%esp),%eax
        movl nb114nf_ninner(%esp),%ebx
        movl nb114nf_outeriter(%ebp),%ecx
        movl nb114nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb114nf_salign(%esp),%eax
        addl %eax,%esp
        addl $904,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


