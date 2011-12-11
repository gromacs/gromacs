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



.globl nb_kernel104_ia32_sse
.globl _nb_kernel104_ia32_sse
nb_kernel104_ia32_sse:  
_nb_kernel104_ia32_sse: 
.set nb104_p_nri, 8
.set nb104_iinr, 12
.set nb104_jindex, 16
.set nb104_jjnr, 20
.set nb104_shift, 24
.set nb104_shiftvec, 28
.set nb104_fshift, 32
.set nb104_gid, 36
.set nb104_pos, 40
.set nb104_faction, 44
.set nb104_charge, 48
.set nb104_p_facel, 52
.set nb104_p_krf, 56
.set nb104_p_crf, 60
.set nb104_Vc, 64
.set nb104_type, 68
.set nb104_p_ntype, 72
.set nb104_vdwparam, 76
.set nb104_Vvdw, 80
.set nb104_p_tabscale, 84
.set nb104_VFtab, 88
.set nb104_invsqrta, 92
.set nb104_dvda, 96
.set nb104_p_gbtabscale, 100
.set nb104_GBtab, 104
.set nb104_p_nthreads, 108
.set nb104_count, 112
.set nb104_mtx, 116
.set nb104_outeriter, 120
.set nb104_inneriter, 124
.set nb104_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use         
.set nb104_ixH1, 0
.set nb104_iyH1, 16
.set nb104_izH1, 32
.set nb104_ixH2, 48
.set nb104_iyH2, 64
.set nb104_izH2, 80
.set nb104_ixM, 96
.set nb104_iyM, 112
.set nb104_izM, 128
.set nb104_jxH1, 144
.set nb104_jyH1, 160
.set nb104_jzH1, 176
.set nb104_jxH2, 192
.set nb104_jyH2, 208
.set nb104_jzH2, 224
.set nb104_jxM, 240
.set nb104_jyM, 256
.set nb104_jzM, 272
.set nb104_dxH1H1, 288
.set nb104_dyH1H1, 304
.set nb104_dzH1H1, 320
.set nb104_dxH1H2, 336
.set nb104_dyH1H2, 352
.set nb104_dzH1H2, 368
.set nb104_dxH1M, 384
.set nb104_dyH1M, 400
.set nb104_dzH1M, 416
.set nb104_dxH2H1, 432
.set nb104_dyH2H1, 448
.set nb104_dzH2H1, 464
.set nb104_dxH2H2, 480
.set nb104_dyH2H2, 496
.set nb104_dzH2H2, 512
.set nb104_dxH2M, 528
.set nb104_dyH2M, 544
.set nb104_dzH2M, 560
.set nb104_dxMH1, 576
.set nb104_dyMH1, 592
.set nb104_dzMH1, 608
.set nb104_dxMH2, 624
.set nb104_dyMH2, 640
.set nb104_dzMH2, 656
.set nb104_dxMM, 672
.set nb104_dyMM, 688
.set nb104_dzMM, 704
.set nb104_qqHH, 720
.set nb104_qqMH, 736
.set nb104_qqMM, 752
.set nb104_vctot, 768
.set nb104_fixH1, 784
.set nb104_fiyH1, 800
.set nb104_fizH1, 816
.set nb104_fixH2, 832
.set nb104_fiyH2, 848
.set nb104_fizH2, 864
.set nb104_fixM, 880
.set nb104_fiyM, 896
.set nb104_fizM, 912
.set nb104_fjxH1, 928
.set nb104_fjyH1, 944
.set nb104_fjzH1, 960
.set nb104_fjxH2, 976
.set nb104_fjyH2, 992
.set nb104_fjzH2, 1008
.set nb104_fjxM, 1024
.set nb104_fjyM, 1040
.set nb104_fjzM, 1056
.set nb104_fjzMb, 1060
.set nb104_fjzMc, 1064
.set nb104_fjzMd, 1068
.set nb104_half, 1072
.set nb104_three, 1088
.set nb104_rsqH1H1, 1104
.set nb104_rsqH1H2, 1120
.set nb104_rsqH1M, 1136
.set nb104_rsqH2H1, 1152
.set nb104_rsqH2H2, 1168
.set nb104_rsqH2M, 1184
.set nb104_rsqMH1, 1200
.set nb104_rsqMH2, 1216
.set nb104_rsqMM, 1232
.set nb104_rinvH1H1, 1248
.set nb104_rinvH1H2, 1264
.set nb104_rinvH1M, 1280
.set nb104_rinvH2H1, 1296
.set nb104_rinvH2H2, 1312
.set nb104_rinvH2M, 1328
.set nb104_rinvMH1, 1344
.set nb104_rinvMH2, 1360
.set nb104_rinvMM, 1376
.set nb104_is3, 1392
.set nb104_ii3, 1396
.set nb104_innerjjnr, 1400
.set nb104_innerk, 1404
.set nb104_n, 1408
.set nb104_nn1, 1412
.set nb104_nri, 1416
.set nb104_nouter, 1420
.set nb104_ninner, 1424
.set nb104_salign, 1428
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $1432,%esp         ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb104_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb104_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb104_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb104_nouter(%esp)
        movl %eax,nb104_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb104_half(%esp)
        movss nb104_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb104_half(%esp)
        movaps %xmm3,nb104_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb104_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx),%ebx               ## ebx =ii 

        movl  nb104_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movss 12(%edx,%ebx,4),%xmm5
        movl nb104_p_facel(%ebp),%esi
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
        movaps %xmm3,nb104_qqHH(%esp)
        movaps %xmm4,nb104_qqMH(%esp)
        movaps %xmm5,nb104_qqMM(%esp)

_nb_kernel104_ia32_sse.nb104_threadloop: 
        movl  nb104_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel104_ia32_sse.nb104_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel104_ia32_sse.nb104_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb104_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb104_n(%esp)
        movl %ebx,nb104_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel104_ia32_sse.nb104_outerstart
        jmp _nb_kernel104_ia32_sse.nb104_end

_nb_kernel104_ia32_sse.nb104_outerstart: 
        ## ebx contains number of outer iterations
        addl nb104_nouter(%esp),%ebx
        movl %ebx,nb104_nouter(%esp)

_nb_kernel104_ia32_sse.nb104_outer: 
        movl  nb104_shift(%ebp),%eax            ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb104_is3(%esp)      ## store is3 

        movl  nb104_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb104_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb104_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb104_ii3(%esp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss 12(%eax,%ebx,4),%xmm3
        addss 16(%eax,%ebx,4),%xmm4
        addss 20(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb104_ixH1(%esp)
        movaps %xmm4,nb104_iyH1(%esp)
        movaps %xmm5,nb104_izH1(%esp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 24(%eax,%ebx,4),%xmm0
        addss 28(%eax,%ebx,4),%xmm1
        addss 32(%eax,%ebx,4),%xmm2
        addss 36(%eax,%ebx,4),%xmm3
        addss 40(%eax,%ebx,4),%xmm4
        addss 44(%eax,%ebx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb104_ixH2(%esp)
        movaps %xmm1,nb104_iyH2(%esp)
        movaps %xmm2,nb104_izH2(%esp)
        movaps %xmm3,nb104_ixM(%esp)
        movaps %xmm4,nb104_iyM(%esp)
        movaps %xmm5,nb104_izM(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb104_vctot(%esp)
        movaps %xmm4,nb104_fixH1(%esp)
        movaps %xmm4,nb104_fiyH1(%esp)
        movaps %xmm4,nb104_fizH1(%esp)
        movaps %xmm4,nb104_fixH2(%esp)
        movaps %xmm4,nb104_fiyH2(%esp)
        movaps %xmm4,nb104_fizH2(%esp)
        movaps %xmm4,nb104_fixM(%esp)
        movaps %xmm4,nb104_fiyM(%esp)
        movaps %xmm4,nb104_fizM(%esp)

        movl  nb104_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb104_pos(%ebp),%esi
        movl  nb104_faction(%ebp),%edi
        movl  nb104_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb104_innerjjnr(%esp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb104_ninner(%esp),%ecx
        movl  %ecx,nb104_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb104_innerk(%esp)   ## number of innerloop atoms 
        jge   _nb_kernel104_ia32_sse.nb104_unroll_loop
        jmp   _nb_kernel104_ia32_sse.nb104_single_check
_nb_kernel104_ia32_sse.nb104_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb104_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb104_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb104_pos(%ebp),%esi       ## base of pos[] 

        leal  (%eax,%eax,2),%eax        ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx        ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move j coordinates to local temp variables 
        movlps 12(%esi,%eax,4),%xmm2
        movlps 24(%esi,%eax,4),%xmm3
        movlps 36(%esi,%eax,4),%xmm4

        movlps 12(%esi,%ebx,4),%xmm5
        movlps 24(%esi,%ebx,4),%xmm6
        movlps 36(%esi,%ebx,4),%xmm7

        movhps 12(%esi,%ecx,4),%xmm2
        movhps 24(%esi,%ecx,4),%xmm3
        movhps 36(%esi,%ecx,4),%xmm4

        movhps 12(%esi,%edx,4),%xmm5
        movhps 24(%esi,%edx,4),%xmm6
        movhps 36(%esi,%edx,4),%xmm7

        ## current state:       
        ## xmm2= jxh1a  jyH1a  jxH1c  jyH1c 
        ## xmm3= jxH2a jyH2a jxH2c jyH2c 
        ## xmm4= jxMa jyMa jxMc jyMc 
        ## xmm5= jxH1b  jyH1b  jxH1d  jyH1d 
        ## xmm6= jxH2b jyH2b jxH2d jyH2d 
        ## xmm7= jxMb jyMb jxMd jyMd 

        movaps %xmm2,%xmm0
        movaps %xmm3,%xmm1
        unpcklps %xmm5,%xmm0    ## xmm0= jxH1a  jxH1b  jyH1a  jyH1b 
        unpcklps %xmm6,%xmm1    ## xmm1= jxH2a jxH2b jyH2a jyH2b 
        unpckhps %xmm5,%xmm2    ## xmm2= jxH1c  jxH1d  jyH1c  jyH1d 
        unpckhps %xmm6,%xmm3    ## xmm3= jxH2c jxH2d jyH2c jyH2d  
        movaps %xmm4,%xmm5
        movaps   %xmm0,%xmm6
        unpcklps %xmm7,%xmm4    ## xmm4= jxMa jxMb jyMa jyMb            
        unpckhps %xmm7,%xmm5    ## xmm5= jxMc jxMd jyMc jyMd     
        movaps   %xmm1,%xmm7
        movlhps  %xmm2,%xmm0    ## xmm0= jxH1a  jxH1b  jxH1c  jxH1d  
        movaps %xmm0,nb104_jxH1(%esp)
        movhlps  %xmm6,%xmm2    ## xmm2= jyH1a  jyH1b  jyH1c  jyH1d 
        movaps %xmm2,nb104_jyH1(%esp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,nb104_jxH2(%esp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,nb104_jyH2(%esp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,nb104_jxM(%esp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,nb104_jyM(%esp)

        movss  20(%esi,%eax,4),%xmm0
        movss  32(%esi,%eax,4),%xmm1
        movss  44(%esi,%eax,4),%xmm2

        movss  20(%esi,%ecx,4),%xmm3
        movss  32(%esi,%ecx,4),%xmm4
        movss  44(%esi,%ecx,4),%xmm5

        movhps 16(%esi,%ebx,4),%xmm0
        movhps 28(%esi,%ebx,4),%xmm1
        movhps 40(%esi,%ebx,4),%xmm2

        movhps 16(%esi,%edx,4),%xmm3
        movhps 28(%esi,%edx,4),%xmm4
        movhps 40(%esi,%edx,4),%xmm5

        shufps $204,%xmm3,%xmm0 ## constant 11001100
        shufps $204,%xmm4,%xmm1 ## constant 11001100
        shufps $204,%xmm5,%xmm2 ## constant 11001100
        movaps %xmm0,nb104_jzH1(%esp)
        movaps %xmm1,nb104_jzH2(%esp)
        movaps %xmm2,nb104_jzM(%esp)

        movaps nb104_ixH1(%esp),%xmm0
        movaps nb104_iyH1(%esp),%xmm1
        movaps nb104_izH1(%esp),%xmm2
        movaps nb104_ixH1(%esp),%xmm3
        movaps nb104_iyH1(%esp),%xmm4
        movaps nb104_izH1(%esp),%xmm5
        subps  nb104_jxH1(%esp),%xmm0
        subps  nb104_jyH1(%esp),%xmm1
        subps  nb104_jzH1(%esp),%xmm2
        subps  nb104_jxH2(%esp),%xmm3
        subps  nb104_jyH2(%esp),%xmm4
        subps  nb104_jzH2(%esp),%xmm5
        movaps %xmm0,nb104_dxH1H1(%esp)
        movaps %xmm1,nb104_dyH1H1(%esp)
        movaps %xmm2,nb104_dzH1H1(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb104_dxH1H2(%esp)
        movaps %xmm4,nb104_dyH1H2(%esp)
        movaps %xmm5,nb104_dzH1H2(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb104_rsqH1H1(%esp)
        movaps %xmm3,nb104_rsqH1H2(%esp)

        movaps nb104_ixH1(%esp),%xmm0
        movaps nb104_iyH1(%esp),%xmm1
        movaps nb104_izH1(%esp),%xmm2
        movaps nb104_ixH2(%esp),%xmm3
        movaps nb104_iyH2(%esp),%xmm4
        movaps nb104_izH2(%esp),%xmm5
        subps  nb104_jxM(%esp),%xmm0
        subps  nb104_jyM(%esp),%xmm1
        subps  nb104_jzM(%esp),%xmm2
        subps  nb104_jxH1(%esp),%xmm3
        subps  nb104_jyH1(%esp),%xmm4
        subps  nb104_jzH1(%esp),%xmm5
        movaps %xmm0,nb104_dxH1M(%esp)
        movaps %xmm1,nb104_dyH1M(%esp)
        movaps %xmm2,nb104_dzH1M(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb104_dxH2H1(%esp)
        movaps %xmm4,nb104_dyH2H1(%esp)
        movaps %xmm5,nb104_dzH2H1(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb104_rsqH1M(%esp)
        movaps %xmm3,nb104_rsqH2H1(%esp)

        movaps nb104_ixH2(%esp),%xmm0
        movaps nb104_iyH2(%esp),%xmm1
        movaps nb104_izH2(%esp),%xmm2
        movaps nb104_ixH2(%esp),%xmm3
        movaps nb104_iyH2(%esp),%xmm4
        movaps nb104_izH2(%esp),%xmm5
        subps  nb104_jxH2(%esp),%xmm0
        subps  nb104_jyH2(%esp),%xmm1
        subps  nb104_jzH2(%esp),%xmm2
        subps  nb104_jxM(%esp),%xmm3
        subps  nb104_jyM(%esp),%xmm4
        subps  nb104_jzM(%esp),%xmm5
        movaps %xmm0,nb104_dxH2H2(%esp)
        movaps %xmm1,nb104_dyH2H2(%esp)
        movaps %xmm2,nb104_dzH2H2(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb104_dxH2M(%esp)
        movaps %xmm4,nb104_dyH2M(%esp)
        movaps %xmm5,nb104_dzH2M(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb104_rsqH2H2(%esp)
        movaps %xmm3,nb104_rsqH2M(%esp)

        movaps nb104_ixM(%esp),%xmm0
        movaps nb104_iyM(%esp),%xmm1
        movaps nb104_izM(%esp),%xmm2
        movaps nb104_ixM(%esp),%xmm3
        movaps nb104_iyM(%esp),%xmm4
        movaps nb104_izM(%esp),%xmm5
        subps  nb104_jxH1(%esp),%xmm0
        subps  nb104_jyH1(%esp),%xmm1
        subps  nb104_jzH1(%esp),%xmm2
        subps  nb104_jxH2(%esp),%xmm3
        subps  nb104_jyH2(%esp),%xmm4
        subps  nb104_jzH2(%esp),%xmm5
        movaps %xmm0,nb104_dxMH1(%esp)
        movaps %xmm1,nb104_dyMH1(%esp)
        movaps %xmm2,nb104_dzMH1(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb104_dxMH2(%esp)
        movaps %xmm4,nb104_dyMH2(%esp)
        movaps %xmm5,nb104_dzMH2(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,nb104_rsqMH1(%esp)
        movaps %xmm4,nb104_rsqMH2(%esp)

        movaps nb104_ixM(%esp),%xmm0
        movaps nb104_iyM(%esp),%xmm1
        movaps nb104_izM(%esp),%xmm2
        subps  nb104_jxM(%esp),%xmm0
        subps  nb104_jyM(%esp),%xmm1
        subps  nb104_jzM(%esp),%xmm2
        movaps %xmm0,nb104_dxMM(%esp)
        movaps %xmm1,nb104_dyMM(%esp)
        movaps %xmm2,nb104_dzMM(%esp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,nb104_rsqMM(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb104_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb104_half(%esp),%xmm3   ## rinvMM 
        mulps   nb104_half(%esp),%xmm7   ## rinvMH2 
        movaps  %xmm3,nb104_rinvMM(%esp)
        movaps  %xmm7,nb104_rinvMH2(%esp)

        rsqrtps nb104_rsqH1H1(%esp),%xmm1
        rsqrtps nb104_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb104_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb104_rsqH1H1(%esp),%xmm1
        mulps   nb104_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb104_half(%esp),%xmm3
        mulps   nb104_half(%esp),%xmm7
        movaps  %xmm3,nb104_rinvH1H1(%esp)
        movaps  %xmm7,nb104_rinvH1H2(%esp)

        rsqrtps nb104_rsqH1M(%esp),%xmm1
        rsqrtps nb104_rsqH2H1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb104_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb104_rsqH1M(%esp),%xmm1
        mulps   nb104_rsqH2H1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb104_half(%esp),%xmm3
        mulps   nb104_half(%esp),%xmm7
        movaps  %xmm3,nb104_rinvH1M(%esp)
        movaps  %xmm7,nb104_rinvH2H1(%esp)

        rsqrtps nb104_rsqH2H2(%esp),%xmm1
        rsqrtps nb104_rsqH2M(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb104_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb104_rsqH2H2(%esp),%xmm1
        mulps   nb104_rsqH2M(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb104_half(%esp),%xmm3
        mulps   nb104_half(%esp),%xmm7
        movaps  %xmm3,nb104_rinvH2H2(%esp)
        movaps  %xmm7,nb104_rinvH2M(%esp)

        rsqrtps nb104_rsqMH1(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb104_three(%esp),%xmm3
        mulps   nb104_rsqMH1(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb104_half(%esp),%xmm3
        movaps  %xmm3,nb104_rinvMH1(%esp)

        ## start with H1-H1 interaction 
        movaps nb104_rinvH1H1(%esp),%xmm0
        movaps %xmm0,%xmm7
        mulps  %xmm0,%xmm0
        mulps  nb104_qqHH(%esp),%xmm7
        mulps  %xmm7,%xmm0
        addps  nb104_vctot(%esp),%xmm7
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb104_dxH1H1(%esp),%xmm0
        mulps nb104_dyH1H1(%esp),%xmm1
        mulps nb104_dzH1H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb104_fixH1(%esp),%xmm0
        addps nb104_fiyH1(%esp),%xmm1
        addps nb104_fizH1(%esp),%xmm2
        movaps %xmm3,nb104_fjxH1(%esp)
        movaps %xmm4,nb104_fjyH1(%esp)
        movaps %xmm5,nb104_fjzH1(%esp)
        movaps %xmm0,nb104_fixH1(%esp)
        movaps %xmm1,nb104_fiyH1(%esp)
        movaps %xmm2,nb104_fizH1(%esp)

        ## H1-H2 interaction 
        movaps nb104_rinvH1H2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb104_qqHH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fs H1-H2  
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb104_dxH1H2(%esp),%xmm0
        mulps nb104_dyH1H2(%esp),%xmm1
        mulps nb104_dzH1H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb104_fixH1(%esp),%xmm0
        addps nb104_fiyH1(%esp),%xmm1
        addps nb104_fizH1(%esp),%xmm2
        movaps %xmm3,nb104_fjxH2(%esp)
        movaps %xmm4,nb104_fjyH2(%esp)
        movaps %xmm5,nb104_fjzH2(%esp)
        movaps %xmm0,nb104_fixH1(%esp)
        movaps %xmm1,nb104_fiyH1(%esp)
        movaps %xmm2,nb104_fizH1(%esp)

        ## H1-M interaction  
        movaps nb104_rinvH1M(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb104_qqMH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fs H1-M  
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb104_dxH1M(%esp),%xmm0
        mulps nb104_dyH1M(%esp),%xmm1
        mulps nb104_dzH1M(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb104_fixH1(%esp),%xmm0
        addps nb104_fiyH1(%esp),%xmm1
        addps nb104_fizH1(%esp),%xmm2
        movaps %xmm3,nb104_fjxM(%esp)
        movaps %xmm4,nb104_fjyM(%esp)
        movaps %xmm5,nb104_fjzM(%esp)
        movaps %xmm0,nb104_fixH1(%esp)
        movaps %xmm1,nb104_fiyH1(%esp)
        movaps %xmm2,nb104_fizH1(%esp)

        ## H2-H1 interaction 
        movaps nb104_rinvH2H1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb104_qqHH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fs H2-H1 
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2
        movaps nb104_fjxH1(%esp),%xmm3
        movaps nb104_fjyH1(%esp),%xmm4
        movaps nb104_fjzH1(%esp),%xmm5
        mulps nb104_dxH2H1(%esp),%xmm0
        mulps nb104_dyH2H1(%esp),%xmm1
        mulps nb104_dzH2H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb104_fixH2(%esp),%xmm0
        addps nb104_fiyH2(%esp),%xmm1
        addps nb104_fizH2(%esp),%xmm2
        movaps %xmm3,nb104_fjxH1(%esp)
        movaps %xmm4,nb104_fjyH1(%esp)
        movaps %xmm5,nb104_fjzH1(%esp)
        movaps %xmm0,nb104_fixH2(%esp)
        movaps %xmm1,nb104_fiyH2(%esp)
        movaps %xmm2,nb104_fizH2(%esp)

        ## H2-H2 interaction 
        movaps nb104_rinvH2H2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb104_qqHH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fsH2H2
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2
        movaps nb104_fjxH2(%esp),%xmm3
        movaps nb104_fjyH2(%esp),%xmm4
        movaps nb104_fjzH2(%esp),%xmm5
        mulps nb104_dxH2H2(%esp),%xmm0
        mulps nb104_dyH2H2(%esp),%xmm1
        mulps nb104_dzH2H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb104_fixH2(%esp),%xmm0
        addps nb104_fiyH2(%esp),%xmm1
        addps nb104_fizH2(%esp),%xmm2
        movaps %xmm3,nb104_fjxH2(%esp)
        movaps %xmm4,nb104_fjyH2(%esp)
        movaps %xmm5,nb104_fjzH2(%esp)
        movaps %xmm0,nb104_fixH2(%esp)
        movaps %xmm1,nb104_fiyH2(%esp)
        movaps %xmm2,nb104_fizH2(%esp)

        ## H2-M interaction 
        movaps nb104_rinvH2M(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb104_qqMH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fs H2-M  
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2
        movaps nb104_fjxM(%esp),%xmm3
        movaps nb104_fjyM(%esp),%xmm4
        movaps nb104_fjzM(%esp),%xmm5
        mulps nb104_dxH2M(%esp),%xmm0
        mulps nb104_dyH2M(%esp),%xmm1
        mulps nb104_dzH2M(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb104_fixH2(%esp),%xmm0
        addps nb104_fiyH2(%esp),%xmm1
        addps nb104_fizH2(%esp),%xmm2
        movaps %xmm3,nb104_fjxM(%esp)
        movaps %xmm4,nb104_fjyM(%esp)
        movaps %xmm5,nb104_fjzM(%esp)
        movaps %xmm0,nb104_fixH2(%esp)
        movaps %xmm1,nb104_fiyH2(%esp)
        movaps %xmm2,nb104_fizH2(%esp)

        ## M-H1 interaction 
        movaps nb104_rinvMH1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb104_qqMH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fs M-H1 
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2
        movaps nb104_fjxH1(%esp),%xmm3
        movaps nb104_fjyH1(%esp),%xmm4
        movaps nb104_fjzH1(%esp),%xmm5
        mulps nb104_dxMH1(%esp),%xmm0
        mulps nb104_dyMH1(%esp),%xmm1
        mulps nb104_dzMH1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb104_fixM(%esp),%xmm0
        addps nb104_fiyM(%esp),%xmm1
        addps nb104_fizM(%esp),%xmm2
        movaps %xmm3,nb104_fjxH1(%esp)
        movaps %xmm4,nb104_fjyH1(%esp)
        movaps %xmm5,nb104_fjzH1(%esp)
        movaps %xmm0,nb104_fixM(%esp)
        movaps %xmm1,nb104_fiyM(%esp)
        movaps %xmm2,nb104_fizM(%esp)

        ## M-H2 interaction 
        movaps nb104_rinvMH2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb104_qqMH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fs M-H2 
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2
        movaps nb104_fjxH2(%esp),%xmm3
        movaps nb104_fjyH2(%esp),%xmm4
        movaps nb104_fjzH2(%esp),%xmm5
        mulps nb104_dxMH2(%esp),%xmm0
        mulps nb104_dyMH2(%esp),%xmm1
        mulps nb104_dzMH2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb104_fixM(%esp),%xmm0
        addps nb104_fiyM(%esp),%xmm1
        addps nb104_fizM(%esp),%xmm2
        movaps %xmm3,nb104_fjxH2(%esp)
        movaps %xmm4,nb104_fjyH2(%esp)
        movaps %xmm5,nb104_fjzH2(%esp)
        movaps %xmm0,nb104_fixM(%esp)
        movaps %xmm1,nb104_fiyM(%esp)
        movaps %xmm2,nb104_fizM(%esp)

        ## M-M interaction 
        movaps nb104_rinvMM(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb104_qqMM(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fs M-M 
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm7,nb104_vctot(%esp)
        movaps %xmm0,%xmm2
        movaps nb104_fjxM(%esp),%xmm3
        movaps nb104_fjyM(%esp),%xmm4
        movaps nb104_fjzM(%esp),%xmm5
        mulps nb104_dxMM(%esp),%xmm0
        mulps nb104_dyMM(%esp),%xmm1
        mulps nb104_dzMM(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb104_fixM(%esp),%xmm0
        addps nb104_fiyM(%esp),%xmm1
        addps nb104_fizM(%esp),%xmm2
        movaps %xmm3,nb104_fjxM(%esp)
        movaps %xmm4,nb104_fjyM(%esp)
        movaps %xmm5,nb104_fjzM(%esp)
        movaps %xmm0,nb104_fixM(%esp)
        movaps %xmm1,nb104_fiyM(%esp)
        movaps %xmm2,nb104_fizM(%esp)

        movl nb104_faction(%ebp),%edi

        ## Did all interactions - now update j forces 
        ## At this stage forces are still on the stack, in positions:
        ## fjxH1, fjyH1, fjzH1, ... , fjzM.
        ## Each position is a quadruplet of forces for the four 
        ## corresponding j waters, so we need to transpose them before
        ## adding to the memory positions.
        ## 
        ## This _used_ to be a simple transpose, but the resulting high number
        ## of unaligned 128-bit load/stores might trigger a possible hardware 
        ## bug on Athlon and Opteron chips, so I have worked around it
        ## to use 64-bit load/stores instead. The performance hit should be
        ## very modest, since the 128-bit unaligned memory instructions were
        ## slow anyway. 

        ## 4 j waters with three atoms each - first do 1st Hydrogen X & Y forces for 4 j particles 
        movaps nb104_fjxH1(%esp),%xmm0   ## xmm0= fjxH1a  fjxH1b  fjxH1c  fjxH1d 
        movaps nb104_fjyH1(%esp),%xmm2   ## xmm1= fjyH1a  fjyH1b  fjyH1c  fjyH1d
        movlps 12(%edi,%eax,4),%xmm3
        movlps 12(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0       ## xmm0= fjxH1a  fjyH1a  fjxH1b  fjyH1b
        unpckhps %xmm2,%xmm1       ## xmm1= fjxH1c  fjyH1c  fjxH1d  fjyH1d
        movhps 12(%edi,%ebx,4),%xmm3
        movhps 12(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,12(%edi,%eax,4)
        movlps %xmm4,12(%edi,%ecx,4)
        movhps %xmm3,12(%edi,%ebx,4)
        movhps %xmm4,12(%edi,%edx,4)

        ## constant 1st Hydrogen Z & 2nd hydrogen X forces for 4 j particles 
        movaps nb104_fjzH1(%esp),%xmm0    ## xmm0= fjzH1a   fjzH1b   fjzH1c   fjzH1d 
        movaps nb104_fjxH2(%esp),%xmm2   ## xmm1= fjxH2a  fjxH2b  fjxH2c  fjxH2d
        movlps 20(%edi,%eax,4),%xmm3
        movlps 20(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0       ## xmm0= fjzH1a  fjxH2a  fjzH1b  fjxH2b
        unpckhps %xmm2,%xmm1       ## xmm1= fjzH1c  fjxH2c  fjzH1d  fjxH2d
        movhps 20(%edi,%ebx,4),%xmm3
        movhps 20(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,20(%edi,%eax,4)
        movlps %xmm4,20(%edi,%ecx,4)
        movhps %xmm3,20(%edi,%ebx,4)
        movhps %xmm4,20(%edi,%edx,4)

        ## constant 2nd hydrogen Y & Z forces for 4 j particles 
        movaps nb104_fjyH2(%esp),%xmm0   ## xmm0= fjyH2a  fjyH2b  fjyH2c  fjyH2d 
        movaps nb104_fjzH2(%esp),%xmm2   ## xmm1= fjzH2a  fjzH2b  fjzH2c  fjzH2d
        movlps 28(%edi,%eax,4),%xmm3
        movlps 28(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0            ## xmm0= fjyH2a  fjzH2a  fjyH2b  fjzH2b
        unpckhps %xmm2,%xmm1            ## xmm1= fjyH2c  fjzH2c  fjyH2d  fjzH2d
        movhps 28(%edi,%ebx,4),%xmm3
        movhps 28(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,28(%edi,%eax,4)
        movlps %xmm4,28(%edi,%ecx,4)
        movhps %xmm3,28(%edi,%ebx,4)
        movhps %xmm4,28(%edi,%edx,4)

        ## Dummy (M) X & Y forces for 4 j particles 
        movaps nb104_fjxM(%esp),%xmm0   ## xmm0= fjxMa  fjxMb  fjxMc  fjxMd 
        movaps nb104_fjyM(%esp),%xmm2   ## xmm1= fjyMa  fjyMb  fjyMc  fjyMd
        movlps 36(%edi,%eax,4),%xmm3
        movlps 36(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0            ## xmm0= fjxMa  fjyMa  fjxMb  fjyMb
        unpckhps %xmm2,%xmm1            ## xmm1= fjxMc  fjyMc  fjxMd  fjyMd
        movhps 36(%edi,%ebx,4),%xmm3
        movhps 36(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,36(%edi,%eax,4)
        movlps %xmm4,36(%edi,%ecx,4)
        movhps %xmm3,36(%edi,%ebx,4)
        movhps %xmm4,36(%edi,%edx,4)


        ## Dummy (M) Z forces for 4 j particles 
        ## Just load the four Z coords into one reg. each
        movss 44(%edi,%eax,4),%xmm4
        movss 44(%edi,%ebx,4),%xmm5
        movss 44(%edi,%ecx,4),%xmm6
        movss 44(%edi,%edx,4),%xmm7
        ## add what we have on the stack
        addss nb104_fjzM(%esp),%xmm4
        addss nb104_fjzMb(%esp),%xmm5
        addss nb104_fjzMc(%esp),%xmm6
        addss nb104_fjzMd(%esp),%xmm7
        ## store back
        movss %xmm4,44(%edi,%eax,4)
        movss %xmm5,44(%edi,%ebx,4)
        movss %xmm6,44(%edi,%ecx,4)
        movss %xmm7,44(%edi,%edx,4)

        ## should we do one more iteration? 
        subl $4,nb104_innerk(%esp)
        jl    _nb_kernel104_ia32_sse.nb104_single_check
        jmp   _nb_kernel104_ia32_sse.nb104_unroll_loop
_nb_kernel104_ia32_sse.nb104_single_check: 
        addl $4,nb104_innerk(%esp)
        jnz   _nb_kernel104_ia32_sse.nb104_single_loop
        jmp   _nb_kernel104_ia32_sse.nb104_updateouterdata
_nb_kernel104_ia32_sse.nb104_single_loop: 
        movl  nb104_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb104_innerjjnr(%esp)

        movl nb104_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        xorps %xmm3,%xmm3
        xorps %xmm4,%xmm4
        xorps %xmm5,%xmm5
        movss 36(%esi,%eax,4),%xmm3             ## jxM  -  -  -
        movss 40(%esi,%eax,4),%xmm4             ## jyM  -  -  -
        movss 44(%esi,%eax,4),%xmm5             ## jzM  -  -  -  

        movlps 12(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1   -    -
        movss  20(%esi,%eax,4),%xmm7            ## xmm7 = jzH1   -    -    - 
        movhps 24(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1 jxH2 jyH2
        movss  32(%esi,%eax,4),%xmm2            ## xmm2 = jzH2   -    -    -

        ## have all coords, time for some shuffling.

        shufps $216,%xmm6,%xmm6 ## constant 11011000     ;# xmm6 = jxH1 jxH2 jyH1 jyH2 
        unpcklps %xmm2,%xmm7                    ## xmm7 = jzH1 jzH2   -    -
        movaps  nb104_ixM(%esp),%xmm0
        movaps  nb104_iyM(%esp),%xmm1
        movaps  nb104_izM(%esp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxM   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## constant 11100100    ;# xmm4 = jyM   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## constant 01000100    ;# xmm5 = jzM   0   jzH1 jzH2

        ## store all j coordinates in jM 
        movaps %xmm3,nb104_jxM(%esp)
        movaps %xmm4,nb104_jyM(%esp)
        movaps %xmm5,nb104_jzM(%esp)
        subps  %xmm3,%xmm0
        subps  %xmm4,%xmm1
        subps  %xmm5,%xmm2
        movaps %xmm0,nb104_dxMM(%esp)
        movaps %xmm1,nb104_dyMM(%esp)
        movaps %xmm2,nb104_dzMM(%esp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb104_three(%esp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb104_half(%esp),%xmm3   ## rinv iM- j water 

        xorps   %xmm1,%xmm1
        movaps  %xmm3,%xmm0
        xorps   %xmm4,%xmm4
        mulps   %xmm0,%xmm0     ## xmm0=rinvsq

        ## fetch charges to xmm4
        movss   nb104_qqMM(%esp),%xmm4
        movhps  nb104_qqMH(%esp),%xmm4

        mulps   %xmm4,%xmm3     ## xmm3=vcoul 
        mulps   %xmm3,%xmm0     ## total fscal 
        addps   nb104_vctot(%esp),%xmm3
        movaps  %xmm3,nb104_vctot(%esp)

        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2
        mulps   nb104_dxMM(%esp),%xmm0
        mulps   nb104_dyMM(%esp),%xmm1
        mulps   nb104_dzMM(%esp),%xmm2
        ## initial update for j forces 
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb104_fjxM(%esp)
        movaps  %xmm4,nb104_fjyM(%esp)
        movaps  %xmm5,nb104_fjzM(%esp)
        addps   nb104_fixM(%esp),%xmm0
        addps   nb104_fiyM(%esp),%xmm1
        addps   nb104_fizM(%esp),%xmm2
        movaps  %xmm0,nb104_fixM(%esp)
        movaps  %xmm1,nb104_fiyM(%esp)
        movaps  %xmm2,nb104_fizM(%esp)

        ## done with i M Now do i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb104_ixH1(%esp),%xmm0
        movaps  nb104_iyH1(%esp),%xmm1
        movaps  nb104_izH1(%esp),%xmm2
        movaps  nb104_ixH2(%esp),%xmm3
        movaps  nb104_iyH2(%esp),%xmm4
        movaps  nb104_izH2(%esp),%xmm5
        subps   nb104_jxM(%esp),%xmm0
        subps   nb104_jyM(%esp),%xmm1
        subps   nb104_jzM(%esp),%xmm2
        subps   nb104_jxM(%esp),%xmm3
        subps   nb104_jyM(%esp),%xmm4
        subps   nb104_jzM(%esp),%xmm5
        movaps %xmm0,nb104_dxH1M(%esp)
        movaps %xmm1,nb104_dyH1M(%esp)
        movaps %xmm2,nb104_dzH1M(%esp)
        movaps %xmm3,nb104_dxH2M(%esp)
        movaps %xmm4,nb104_dyH2M(%esp)
        movaps %xmm5,nb104_dzH2M(%esp)
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

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb104_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb104_half(%esp),%xmm3   ## rinv H1 - j water 
        mulps   nb104_half(%esp),%xmm7   ## rinv H2 - j water  

        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        movss   nb104_qqMH(%esp),%xmm6
        movhps  nb104_qqHH(%esp),%xmm6

        ## do coulomb interaction 
        movaps  %xmm3,%xmm0
        movaps  %xmm7,%xmm4
        mulps   %xmm0,%xmm0     ## rinvsq 
        mulps   %xmm4,%xmm4     ## rinvsq 
        mulps   %xmm6,%xmm3     ## vcoul 
        mulps   %xmm6,%xmm7     ## vcoul 
        movaps  %xmm3,%xmm2
        addps   %xmm7,%xmm2     ## total vcoul 
        mulps   %xmm3,%xmm0     ## fscal 

        addps   nb104_vctot(%esp),%xmm2
        mulps   %xmm4,%xmm7     ## fscal 
        movaps  %xmm2,nb104_vctot(%esp)
        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2
        mulps   nb104_dxH1M(%esp),%xmm0
        mulps   nb104_dyH1M(%esp),%xmm1
        mulps   nb104_dzH1M(%esp),%xmm2
        ## update forces H1 - j water 
        movaps  nb104_fjxM(%esp),%xmm3
        movaps  nb104_fjyM(%esp),%xmm4
        movaps  nb104_fjzM(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb104_fjxM(%esp)
        movaps  %xmm4,nb104_fjyM(%esp)
        movaps  %xmm5,nb104_fjzM(%esp)
        addps   nb104_fixH1(%esp),%xmm0
        addps   nb104_fiyH1(%esp),%xmm1
        addps   nb104_fizH1(%esp),%xmm2
        movaps  %xmm0,nb104_fixH1(%esp)
        movaps  %xmm1,nb104_fiyH1(%esp)
        movaps  %xmm2,nb104_fizH1(%esp)
        ## do forces H2 - j water 
        movaps %xmm7,%xmm0
        movaps %xmm7,%xmm1
        movaps %xmm7,%xmm2
        mulps   nb104_dxH2M(%esp),%xmm0
        mulps   nb104_dyH2M(%esp),%xmm1
        mulps   nb104_dzH2M(%esp),%xmm2
        movaps  nb104_fjxM(%esp),%xmm3
        movaps  nb104_fjyM(%esp),%xmm4
        movaps  nb104_fjzM(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movl    nb104_faction(%ebp),%esi
        movaps  %xmm3,nb104_fjxM(%esp)
        movaps  %xmm4,nb104_fjyM(%esp)
        movaps  %xmm5,nb104_fjzM(%esp)
        addps   nb104_fixH2(%esp),%xmm0
        addps   nb104_fiyH2(%esp),%xmm1
        addps   nb104_fizH2(%esp),%xmm2
        movaps  %xmm0,nb104_fixH2(%esp)
        movaps  %xmm1,nb104_fiyH2(%esp)
        movaps  %xmm2,nb104_fizH2(%esp)

        ## update j water forces from local variables 
        movlps  36(%esi,%eax,4),%xmm0
        movlps  12(%esi,%eax,4),%xmm1
        movhps  24(%esi,%eax,4),%xmm1
        movaps  nb104_fjxM(%esp),%xmm3
        movaps  nb104_fjyM(%esp),%xmm4
        movaps  nb104_fjzM(%esp),%xmm5
        movaps  %xmm5,%xmm6
        movaps  %xmm5,%xmm7
        shufps $2,%xmm6,%xmm6 ## constant 00000010
        shufps $3,%xmm7,%xmm7 ## constant 00000011
        addss   44(%esi,%eax,4),%xmm5
        addss   20(%esi,%eax,4),%xmm6
        addss   32(%esi,%eax,4),%xmm7
        movss   %xmm5,44(%esi,%eax,4)
        movss   %xmm6,20(%esi,%eax,4)
        movss   %xmm7,32(%esi,%eax,4)
        movaps   %xmm3,%xmm5
        unpcklps %xmm4,%xmm3
        unpckhps %xmm4,%xmm5
        addps    %xmm3,%xmm0
        addps    %xmm5,%xmm1
        movlps  %xmm0,36(%esi,%eax,4)
        movlps  %xmm1,12(%esi,%eax,4)
        movhps  %xmm1,24(%esi,%eax,4)

        decl  nb104_innerk(%esp)
        jz    _nb_kernel104_ia32_sse.nb104_updateouterdata
        jmp   _nb_kernel104_ia32_sse.nb104_single_loop
_nb_kernel104_ia32_sse.nb104_updateouterdata: 
        movl  nb104_ii3(%esp),%ecx
        movl  nb104_faction(%ebp),%edi
        movl  nb104_fshift(%ebp),%esi
        movl  nb104_is3(%esp),%edx

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movaps nb104_fixH1(%esp),%xmm0
        movaps nb104_fiyH1(%esp),%xmm1
        movaps nb104_fizH1(%esp),%xmm2

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
        movaps %xmm0,%xmm6
        movss %xmm2,%xmm7
        movlhps %xmm1,%xmm6
        shufps $8,%xmm6,%xmm6 ## constant 00001000      

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movaps nb104_fixH2(%esp),%xmm0
        movaps nb104_fiyH2(%esp),%xmm1
        movaps nb104_fizH2(%esp),%xmm2

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

        ## accumulate M i forces in xmm0, xmm1, xmm2 
        movaps nb104_fixM(%esp),%xmm0
        movaps nb104_fiyM(%esp),%xmm1
        movaps nb104_fizM(%esp),%xmm2

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
        movl nb104_n(%esp),%esi
        ## get group index for i particle 
        movl  nb104_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb104_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb104_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb104_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel104_ia32_sse.nb104_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb104_n(%esp)
        jmp _nb_kernel104_ia32_sse.nb104_outer
_nb_kernel104_ia32_sse.nb104_outerend: 
        ## check if more outer neighborlists remain
        movl  nb104_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel104_ia32_sse.nb104_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel104_ia32_sse.nb104_threadloop
_nb_kernel104_ia32_sse.nb104_end: 
        emms

        movl nb104_nouter(%esp),%eax
        movl nb104_ninner(%esp),%ebx
        movl nb104_outeriter(%ebp),%ecx
        movl nb104_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb104_salign(%esp),%eax
        addl %eax,%esp
        addl $1432,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret




.globl nb_kernel104nf_ia32_sse
.globl _nb_kernel104nf_ia32_sse
nb_kernel104nf_ia32_sse:        
_nb_kernel104nf_ia32_sse:       
.set nb104nf_p_nri, 8
.set nb104nf_iinr, 12
.set nb104nf_jindex, 16
.set nb104nf_jjnr, 20
.set nb104nf_shift, 24
.set nb104nf_shiftvec, 28
.set nb104nf_fshift, 32
.set nb104nf_gid, 36
.set nb104nf_pos, 40
.set nb104nf_faction, 44
.set nb104nf_charge, 48
.set nb104nf_p_facel, 52
.set nb104nf_p_krf, 56
.set nb104nf_p_crf, 60
.set nb104nf_Vc, 64
.set nb104nf_type, 68
.set nb104nf_p_ntype, 72
.set nb104nf_vdwparam, 76
.set nb104nf_Vvdw, 80
.set nb104nf_p_tabscale, 84
.set nb104nf_VFtab, 88
.set nb104nf_invsqrta, 92
.set nb104nf_dvda, 96
.set nb104nf_p_gbtabscale, 100
.set nb104nf_GBtab, 104
.set nb104nf_p_nthreads, 108
.set nb104nf_count, 112
.set nb104nf_mtx, 116
.set nb104nf_outeriter, 120
.set nb104nf_inneriter, 124
.set nb104nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use         
.set nb104nf_ixH1, 0
.set nb104nf_iyH1, 16
.set nb104nf_izH1, 32
.set nb104nf_ixH2, 48
.set nb104nf_iyH2, 64
.set nb104nf_izH2, 80
.set nb104nf_ixM, 96
.set nb104nf_iyM, 112
.set nb104nf_izM, 128
.set nb104nf_jxH1, 144
.set nb104nf_jyH1, 160
.set nb104nf_jzH1, 176
.set nb104nf_jxH2, 192
.set nb104nf_jyH2, 208
.set nb104nf_jzH2, 224
.set nb104nf_jxM, 240
.set nb104nf_jyM, 256
.set nb104nf_jzM, 272
.set nb104nf_dxMM, 288
.set nb104nf_dyMM, 304
.set nb104nf_dzMM, 320
.set nb104nf_qqHH, 336
.set nb104nf_qqMH, 352
.set nb104nf_qqMM, 368
.set nb104nf_vctot, 384
.set nb104nf_half, 400
.set nb104nf_three, 416
.set nb104nf_rsqH1H1, 432
.set nb104nf_rsqH1H2, 448
.set nb104nf_rsqH1M, 464
.set nb104nf_rsqH2H1, 480
.set nb104nf_rsqH2H2, 496
.set nb104nf_rsqH2M, 512
.set nb104nf_rsqMH1, 528
.set nb104nf_rsqMH2, 544
.set nb104nf_rsqMM, 560
.set nb104nf_rinvH1H1, 576
.set nb104nf_rinvH1H2, 592
.set nb104nf_rinvH1M, 608
.set nb104nf_rinvH2H1, 624
.set nb104nf_rinvH2H2, 640
.set nb104nf_rinvH2M, 656
.set nb104nf_rinvMH1, 672
.set nb104nf_rinvMH2, 688
.set nb104nf_rinvMM, 704
.set nb104nf_is3, 720
.set nb104nf_ii3, 724
.set nb104nf_innerjjnr, 728
.set nb104nf_innerk, 732
.set nb104nf_n, 736
.set nb104nf_nn1, 740
.set nb104nf_nri, 744
.set nb104nf_nouter, 748
.set nb104nf_ninner, 752
.set nb104nf_salign, 756
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $760,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb104nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb104nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb104nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb104nf_nouter(%esp)
        movl %eax,nb104nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb104nf_half(%esp)
        movss nb104nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb104nf_half(%esp)
        movaps %xmm3,nb104nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb104nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx),%ebx               ## ebx =ii 

        movl  nb104nf_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movss 12(%edx,%ebx,4),%xmm5
        movl nb104nf_p_facel(%ebp),%esi
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
        movaps %xmm3,nb104nf_qqHH(%esp)
        movaps %xmm4,nb104nf_qqMH(%esp)
        movaps %xmm5,nb104nf_qqMM(%esp)

_nb_kernel104nf_ia32_sse.nb104nf_threadloop: 
        movl  nb104nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel104nf_ia32_sse.nb104nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel104nf_ia32_sse.nb104nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb104nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb104nf_n(%esp)
        movl %ebx,nb104nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel104nf_ia32_sse.nb104nf_outerstart
        jmp _nb_kernel104nf_ia32_sse.nb104nf_end

_nb_kernel104nf_ia32_sse.nb104nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb104nf_nouter(%esp),%ebx
        movl %ebx,nb104nf_nouter(%esp)

_nb_kernel104nf_ia32_sse.nb104nf_outer: 
        movl  nb104nf_shift(%ebp),%eax          ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb104nf_is3(%esp)            ## store is3 

        movl  nb104nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb104nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb104nf_pos(%ebp),%eax    ## eax = base of pos[]  
        movl  %ebx,nb104nf_ii3(%esp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss 12(%eax,%ebx,4),%xmm3
        addss 16(%eax,%ebx,4),%xmm4
        addss 20(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb104nf_ixH1(%esp)
        movaps %xmm4,nb104nf_iyH1(%esp)
        movaps %xmm5,nb104nf_izH1(%esp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 24(%eax,%ebx,4),%xmm0
        addss 28(%eax,%ebx,4),%xmm1
        addss 32(%eax,%ebx,4),%xmm2
        addss 36(%eax,%ebx,4),%xmm3
        addss 40(%eax,%ebx,4),%xmm4
        addss 44(%eax,%ebx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb104nf_ixH2(%esp)
        movaps %xmm1,nb104nf_iyH2(%esp)
        movaps %xmm2,nb104nf_izH2(%esp)
        movaps %xmm3,nb104nf_ixM(%esp)
        movaps %xmm4,nb104nf_iyM(%esp)
        movaps %xmm5,nb104nf_izM(%esp)

        ## clear vctot
        xorps %xmm4,%xmm4
        movaps %xmm4,nb104nf_vctot(%esp)

        movl  nb104nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb104nf_pos(%ebp),%esi
        movl  nb104nf_faction(%ebp),%edi
        movl  nb104nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb104nf_innerjjnr(%esp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb104nf_ninner(%esp),%ecx
        movl  %ecx,nb104nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb104nf_innerk(%esp)         ## number of innerloop atoms 
        jge   _nb_kernel104nf_ia32_sse.nb104nf_unroll_loop
        jmp   _nb_kernel104nf_ia32_sse.nb104nf_single_check
_nb_kernel104nf_ia32_sse.nb104nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb104nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb104nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb104nf_pos(%ebp),%esi     ## base of pos[] 

        leal  (%eax,%eax,2),%eax        ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx        ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move j coordinates to local temp variables 
        movlps 12(%esi,%eax,4),%xmm2
        movlps 24(%esi,%eax,4),%xmm3
        movlps 36(%esi,%eax,4),%xmm4

        movlps 12(%esi,%ebx,4),%xmm5
        movlps 24(%esi,%ebx,4),%xmm6
        movlps 36(%esi,%ebx,4),%xmm7

        movhps 12(%esi,%ecx,4),%xmm2
        movhps 24(%esi,%ecx,4),%xmm3
        movhps 36(%esi,%ecx,4),%xmm4

        movhps 12(%esi,%edx,4),%xmm5
        movhps 24(%esi,%edx,4),%xmm6
        movhps 36(%esi,%edx,4),%xmm7

        ## current state:       
        ## xmm2= jxh1a  jyH1a  jxH1c  jyH1c 
        ## xmm3= jxH2a jyH2a jxH2c jyH2c 
        ## xmm4= jxMa jyMa jxMc jyMc 
        ## xmm5= jxH1b  jyH1b  jxH1d  jyH1d 
        ## xmm6= jxH2b jyH2b jxH2d jyH2d 
        ## xmm7= jxMb jyMb jxMd jyMd 

        movaps %xmm2,%xmm0
        movaps %xmm3,%xmm1
        unpcklps %xmm5,%xmm0    ## xmm0= jxH1a  jxH1b  jyH1a  jyH1b 
        unpcklps %xmm6,%xmm1    ## xmm1= jxH2a jxH2b jyH2a jyH2b 
        unpckhps %xmm5,%xmm2    ## xmm2= jxH1c  jxH1d  jyH1c  jyH1d 
        unpckhps %xmm6,%xmm3    ## xmm3= jxH2c jxH2d jyH2c jyH2d  
        movaps %xmm4,%xmm5
        movaps   %xmm0,%xmm6
        unpcklps %xmm7,%xmm4    ## xmm4= jxMa jxMb jyMa jyMb            
        unpckhps %xmm7,%xmm5    ## xmm5= jxMc jxMd jyMc jyMd     
        movaps   %xmm1,%xmm7
        movlhps  %xmm2,%xmm0    ## xmm0= jxH1a  jxH1b  jxH1c  jxH1d  
        movaps %xmm0,nb104nf_jxH1(%esp)
        movhlps  %xmm6,%xmm2    ## xmm2= jyH1a  jyH1b  jyH1c  jyH1d 
        movaps %xmm2,nb104nf_jyH1(%esp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,nb104nf_jxH2(%esp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,nb104nf_jyH2(%esp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,nb104nf_jxM(%esp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,nb104nf_jyM(%esp)

        movss  20(%esi,%eax,4),%xmm0
        movss  32(%esi,%eax,4),%xmm1
        movss  44(%esi,%eax,4),%xmm2

        movss  20(%esi,%ecx,4),%xmm3
        movss  32(%esi,%ecx,4),%xmm4
        movss  44(%esi,%ecx,4),%xmm5

        movhps 16(%esi,%ebx,4),%xmm0
        movhps 28(%esi,%ebx,4),%xmm1
        movhps 40(%esi,%ebx,4),%xmm2

        movhps 16(%esi,%edx,4),%xmm3
        movhps 28(%esi,%edx,4),%xmm4
        movhps 40(%esi,%edx,4),%xmm5

        shufps $204,%xmm3,%xmm0 ## constant 11001100
        shufps $204,%xmm4,%xmm1 ## constant 11001100
        shufps $204,%xmm5,%xmm2 ## constant 11001100
        movaps %xmm0,nb104nf_jzH1(%esp)
        movaps %xmm1,nb104nf_jzH2(%esp)
        movaps %xmm2,nb104nf_jzM(%esp)

        movaps nb104nf_ixH1(%esp),%xmm0
        movaps nb104nf_iyH1(%esp),%xmm1
        movaps nb104nf_izH1(%esp),%xmm2
        movaps nb104nf_ixH1(%esp),%xmm3
        movaps nb104nf_iyH1(%esp),%xmm4
        movaps nb104nf_izH1(%esp),%xmm5
        subps  nb104nf_jxH1(%esp),%xmm0
        subps  nb104nf_jyH1(%esp),%xmm1
        subps  nb104nf_jzH1(%esp),%xmm2
        subps  nb104nf_jxH2(%esp),%xmm3
        subps  nb104nf_jyH2(%esp),%xmm4
        subps  nb104nf_jzH2(%esp),%xmm5
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
        movaps %xmm0,nb104nf_rsqH1H1(%esp)
        movaps %xmm3,nb104nf_rsqH1H2(%esp)

        movaps nb104nf_ixH1(%esp),%xmm0
        movaps nb104nf_iyH1(%esp),%xmm1
        movaps nb104nf_izH1(%esp),%xmm2
        movaps nb104nf_ixH2(%esp),%xmm3
        movaps nb104nf_iyH2(%esp),%xmm4
        movaps nb104nf_izH2(%esp),%xmm5
        subps  nb104nf_jxM(%esp),%xmm0
        subps  nb104nf_jyM(%esp),%xmm1
        subps  nb104nf_jzM(%esp),%xmm2
        subps  nb104nf_jxH1(%esp),%xmm3
        subps  nb104nf_jyH1(%esp),%xmm4
        subps  nb104nf_jzH1(%esp),%xmm5
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
        movaps %xmm0,nb104nf_rsqH1M(%esp)
        movaps %xmm3,nb104nf_rsqH2H1(%esp)

        movaps nb104nf_ixH2(%esp),%xmm0
        movaps nb104nf_iyH2(%esp),%xmm1
        movaps nb104nf_izH2(%esp),%xmm2
        movaps nb104nf_ixH2(%esp),%xmm3
        movaps nb104nf_iyH2(%esp),%xmm4
        movaps nb104nf_izH2(%esp),%xmm5
        subps  nb104nf_jxH2(%esp),%xmm0
        subps  nb104nf_jyH2(%esp),%xmm1
        subps  nb104nf_jzH2(%esp),%xmm2
        subps  nb104nf_jxM(%esp),%xmm3
        subps  nb104nf_jyM(%esp),%xmm4
        subps  nb104nf_jzM(%esp),%xmm5
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
        movaps %xmm0,nb104nf_rsqH2H2(%esp)
        movaps %xmm3,nb104nf_rsqH2M(%esp)

        movaps nb104nf_ixM(%esp),%xmm0
        movaps nb104nf_iyM(%esp),%xmm1
        movaps nb104nf_izM(%esp),%xmm2
        movaps nb104nf_ixM(%esp),%xmm3
        movaps nb104nf_iyM(%esp),%xmm4
        movaps nb104nf_izM(%esp),%xmm5
        subps  nb104nf_jxH1(%esp),%xmm0
        subps  nb104nf_jyH1(%esp),%xmm1
        subps  nb104nf_jzH1(%esp),%xmm2
        subps  nb104nf_jxH2(%esp),%xmm3
        subps  nb104nf_jyH2(%esp),%xmm4
        subps  nb104nf_jzH2(%esp),%xmm5
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
        movaps %xmm0,nb104nf_rsqMH1(%esp)
        movaps %xmm4,nb104nf_rsqMH2(%esp)

        movaps nb104nf_ixM(%esp),%xmm0
        movaps nb104nf_iyM(%esp),%xmm1
        movaps nb104nf_izM(%esp),%xmm2
        subps  nb104nf_jxM(%esp),%xmm0
        subps  nb104nf_jyM(%esp),%xmm1
        subps  nb104nf_jzM(%esp),%xmm2
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,nb104nf_rsqMM(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb104nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb104nf_half(%esp),%xmm3   ## rinvMM 
        mulps   nb104nf_half(%esp),%xmm7   ## rinvMH2 
        movaps  %xmm3,nb104nf_rinvMM(%esp)
        movaps  %xmm7,nb104nf_rinvMH2(%esp)

        rsqrtps nb104nf_rsqH1H1(%esp),%xmm1
        rsqrtps nb104nf_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb104nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb104nf_rsqH1H1(%esp),%xmm1
        mulps   nb104nf_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb104nf_half(%esp),%xmm3
        mulps   nb104nf_half(%esp),%xmm7
        movaps  %xmm3,nb104nf_rinvH1H1(%esp)
        movaps  %xmm7,nb104nf_rinvH1H2(%esp)

        rsqrtps nb104nf_rsqH1M(%esp),%xmm1
        rsqrtps nb104nf_rsqH2H1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb104nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb104nf_rsqH1M(%esp),%xmm1
        mulps   nb104nf_rsqH2H1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb104nf_half(%esp),%xmm3
        mulps   nb104nf_half(%esp),%xmm7
        movaps  %xmm3,nb104nf_rinvH1M(%esp)
        movaps  %xmm7,nb104nf_rinvH2H1(%esp)

        rsqrtps nb104nf_rsqH2H2(%esp),%xmm1
        rsqrtps nb104nf_rsqH2M(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb104nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb104nf_rsqH2H2(%esp),%xmm1
        mulps   nb104nf_rsqH2M(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb104nf_half(%esp),%xmm3
        mulps   nb104nf_half(%esp),%xmm7
        movaps  %xmm3,nb104nf_rinvH2H2(%esp)
        movaps  %xmm7,nb104nf_rinvH2M(%esp)

        rsqrtps nb104nf_rsqMH1(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb104nf_three(%esp),%xmm3
        mulps   nb104nf_rsqMH1(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb104nf_half(%esp),%xmm3
        movaps  %xmm3,nb104nf_rinvMH1(%esp)

        ## all H-H interactions
        movaps nb104nf_rinvH1H1(%esp),%xmm0
        addps  nb104nf_rinvH1H2(%esp),%xmm0
        addps  nb104nf_rinvH2H1(%esp),%xmm0
        addps  nb104nf_rinvH2H2(%esp),%xmm0
        mulps  nb104nf_qqHH(%esp),%xmm0
        ## all M-H interactions
        movaps nb104nf_rinvH1M(%esp),%xmm1
        addps  nb104nf_rinvH2M(%esp),%xmm1
        addps  nb104nf_rinvMH1(%esp),%xmm1
        addps  nb104nf_rinvMH2(%esp),%xmm1
        mulps  nb104nf_qqMH(%esp),%xmm1
        ## The M-M interaction
        movaps nb104nf_rinvMM(%esp),%xmm2
        mulps  nb104nf_qqMM(%esp),%xmm2
        addps  %xmm1,%xmm0
        addps  nb104nf_vctot(%esp),%xmm2
        addps  %xmm2,%xmm0
        movaps %xmm0,nb104nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb104nf_innerk(%esp)
        jl    _nb_kernel104nf_ia32_sse.nb104nf_single_check
        jmp   _nb_kernel104nf_ia32_sse.nb104nf_unroll_loop
_nb_kernel104nf_ia32_sse.nb104nf_single_check: 
        addl $4,nb104nf_innerk(%esp)
        jnz   _nb_kernel104nf_ia32_sse.nb104nf_single_loop
        jmp   _nb_kernel104nf_ia32_sse.nb104nf_updateouterdata
_nb_kernel104nf_ia32_sse.nb104nf_single_loop: 
        movl  nb104nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb104nf_innerjjnr(%esp)

        movl nb104nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        xorps %xmm3,%xmm3
        xorps %xmm4,%xmm4
        xorps %xmm5,%xmm5
        movss 36(%esi,%eax,4),%xmm3             ## jxM  -  -  -
        movss 40(%esi,%eax,4),%xmm4             ## jyM  -  -  -
        movss 44(%esi,%eax,4),%xmm5             ## jzM  -  -  -  

        movlps 12(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1   -    -
        movss  20(%esi,%eax,4),%xmm7            ## xmm7 = jzH1   -    -    - 
        movhps 24(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1 jxH2 jyH2
        movss  32(%esi,%eax,4),%xmm2            ## xmm2 = jzH2   -    -    -

        ## have all coords, time for some shuffling.

        shufps $216,%xmm6,%xmm6 ## constant 11011000     ;# xmm6 = jxH1 jxH2 jyH1 jyH2 
        unpcklps %xmm2,%xmm7                    ## xmm7 = jzH1 jzH2   -    -
        movaps  nb104nf_ixM(%esp),%xmm0
        movaps  nb104nf_iyM(%esp),%xmm1
        movaps  nb104nf_izM(%esp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxM   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## constant 11100100    ;# xmm4 = jyM   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## constant 01000100    ;# xmm5 = jzM   0   jzH1 jzH2

        ## store all j coordinates in jM 
        movaps %xmm3,nb104nf_jxM(%esp)
        movaps %xmm4,nb104nf_jyM(%esp)
        movaps %xmm5,nb104nf_jzM(%esp)
        subps  %xmm3,%xmm0
        subps  %xmm4,%xmm1
        subps  %xmm5,%xmm2
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb104nf_three(%esp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb104nf_half(%esp),%xmm3   ## rinv iM- j water 

        xorps   %xmm1,%xmm1
        movaps  %xmm3,%xmm0
        xorps   %xmm4,%xmm4
        mulps   %xmm0,%xmm0     ## xmm0=rinvsq

        ## fetch charges to xmm4
        movss   nb104nf_qqMM(%esp),%xmm4
        movhps  nb104nf_qqMH(%esp),%xmm4

        mulps   %xmm4,%xmm3     ## xmm3=vcoul 
        addps   nb104nf_vctot(%esp),%xmm3
        movaps  %xmm3,nb104nf_vctot(%esp)

        ## done with i M Now do i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb104nf_ixH1(%esp),%xmm0
        movaps  nb104nf_iyH1(%esp),%xmm1
        movaps  nb104nf_izH1(%esp),%xmm2
        movaps  nb104nf_ixH2(%esp),%xmm3
        movaps  nb104nf_iyH2(%esp),%xmm4
        movaps  nb104nf_izH2(%esp),%xmm5
        subps   nb104nf_jxM(%esp),%xmm0
        subps   nb104nf_jyM(%esp),%xmm1
        subps   nb104nf_jzM(%esp),%xmm2
        subps   nb104nf_jxM(%esp),%xmm3
        subps   nb104nf_jyM(%esp),%xmm4
        subps   nb104nf_jzM(%esp),%xmm5
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

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb104nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb104nf_half(%esp),%xmm3   ## rinv H1 - j water 
        mulps   nb104nf_half(%esp),%xmm7   ## rinv H2 - j water  

        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        movss   nb104nf_qqMH(%esp),%xmm6
        movhps  nb104nf_qqHH(%esp),%xmm6

        ## do coulomb interaction 
        mulps   %xmm6,%xmm3     ## vcoul 
        mulps   %xmm6,%xmm7     ## vcoul 
        addps   %xmm7,%xmm3     ## total vcoul 
        addps   nb104nf_vctot(%esp),%xmm3
        movaps  %xmm3,nb104nf_vctot(%esp)

        decl  nb104nf_innerk(%esp)
        jz    _nb_kernel104nf_ia32_sse.nb104nf_updateouterdata
        jmp   _nb_kernel104nf_ia32_sse.nb104nf_single_loop
_nb_kernel104nf_ia32_sse.nb104nf_updateouterdata: 
        ## get n from stack
        movl nb104nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb104nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb104nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb104nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb104nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel104nf_ia32_sse.nb104nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb104nf_n(%esp)
        jmp _nb_kernel104nf_ia32_sse.nb104nf_outer
_nb_kernel104nf_ia32_sse.nb104nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb104nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel104nf_ia32_sse.nb104nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel104nf_ia32_sse.nb104nf_threadloop
_nb_kernel104nf_ia32_sse.nb104nf_end: 
        emms

        movl nb104nf_nouter(%esp),%eax
        movl nb104nf_ninner(%esp),%ebx
        movl nb104nf_outeriter(%ebp),%ecx
        movl nb104nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb104nf_salign(%esp),%eax
        addl %eax,%esp
        addl $760,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret




