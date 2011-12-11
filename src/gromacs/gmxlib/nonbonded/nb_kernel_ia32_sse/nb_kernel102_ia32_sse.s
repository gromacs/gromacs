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



.globl nb_kernel102_ia32_sse
.globl _nb_kernel102_ia32_sse
nb_kernel102_ia32_sse:  
_nb_kernel102_ia32_sse: 
.set nb102_p_nri, 8
.set nb102_iinr, 12
.set nb102_jindex, 16
.set nb102_jjnr, 20
.set nb102_shift, 24
.set nb102_shiftvec, 28
.set nb102_fshift, 32
.set nb102_gid, 36
.set nb102_pos, 40
.set nb102_faction, 44
.set nb102_charge, 48
.set nb102_p_facel, 52
.set nb102_p_krf, 56
.set nb102_p_crf, 60
.set nb102_Vc, 64
.set nb102_type, 68
.set nb102_p_ntype, 72
.set nb102_vdwparam, 76
.set nb102_Vvdw, 80
.set nb102_p_tabscale, 84
.set nb102_VFtab, 88
.set nb102_invsqrta, 92
.set nb102_dvda, 96
.set nb102_p_gbtabscale, 100
.set nb102_GBtab, 104
.set nb102_p_nthreads, 108
.set nb102_count, 112
.set nb102_mtx, 116
.set nb102_outeriter, 120
.set nb102_inneriter, 124
.set nb102_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use         
.set nb102_ixO, 0
.set nb102_iyO, 16
.set nb102_izO, 32
.set nb102_ixH1, 48
.set nb102_iyH1, 64
.set nb102_izH1, 80
.set nb102_ixH2, 96
.set nb102_iyH2, 112
.set nb102_izH2, 128
.set nb102_jxO, 144
.set nb102_jyO, 160
.set nb102_jzO, 176
.set nb102_jxH1, 192
.set nb102_jyH1, 208
.set nb102_jzH1, 224
.set nb102_jxH2, 240
.set nb102_jyH2, 256
.set nb102_jzH2, 272
.set nb102_dxOO, 288
.set nb102_dyOO, 304
.set nb102_dzOO, 320
.set nb102_dxOH1, 336
.set nb102_dyOH1, 352
.set nb102_dzOH1, 368
.set nb102_dxOH2, 384
.set nb102_dyOH2, 400
.set nb102_dzOH2, 416
.set nb102_dxH1O, 432
.set nb102_dyH1O, 448
.set nb102_dzH1O, 464
.set nb102_dxH1H1, 480
.set nb102_dyH1H1, 496
.set nb102_dzH1H1, 512
.set nb102_dxH1H2, 528
.set nb102_dyH1H2, 544
.set nb102_dzH1H2, 560
.set nb102_dxH2O, 576
.set nb102_dyH2O, 592
.set nb102_dzH2O, 608
.set nb102_dxH2H1, 624
.set nb102_dyH2H1, 640
.set nb102_dzH2H1, 656
.set nb102_dxH2H2, 672
.set nb102_dyH2H2, 688
.set nb102_dzH2H2, 704
.set nb102_qqOO, 720
.set nb102_qqOH, 736
.set nb102_qqHH, 752
.set nb102_vctot, 768
.set nb102_fixO, 784
.set nb102_fiyO, 800
.set nb102_fizO, 816
.set nb102_fixH1, 832
.set nb102_fiyH1, 848
.set nb102_fizH1, 864
.set nb102_fixH2, 880
.set nb102_fiyH2, 896
.set nb102_fizH2, 912
.set nb102_fjxO, 928
.set nb102_fjyO, 944
.set nb102_fjzO, 960
.set nb102_fjxH1, 976
.set nb102_fjyH1, 992
.set nb102_fjzH1, 1008
.set nb102_fjxH2, 1024
.set nb102_fjyH2, 1040
.set nb102_fjzH2, 1056
.set nb102_fjzH2b, 1060
.set nb102_fjzH2c, 1064
.set nb102_fjzH2d, 1068
.set nb102_half, 1072
.set nb102_three, 1088
.set nb102_rsqOO, 1104
.set nb102_rsqOH1, 1120
.set nb102_rsqOH2, 1136
.set nb102_rsqH1O, 1152
.set nb102_rsqH1H1, 1168
.set nb102_rsqH1H2, 1184
.set nb102_rsqH2O, 1200
.set nb102_rsqH2H1, 1216
.set nb102_rsqH2H2, 1232
.set nb102_rinvOO, 1248
.set nb102_rinvOH1, 1264
.set nb102_rinvOH2, 1280
.set nb102_rinvH1O, 1296
.set nb102_rinvH1H1, 1312
.set nb102_rinvH1H2, 1328
.set nb102_rinvH2O, 1344
.set nb102_rinvH2H1, 1360
.set nb102_rinvH2H2, 1376
.set nb102_is3, 1392
.set nb102_ii3, 1396
.set nb102_innerjjnr, 1400
.set nb102_innerk, 1404
.set nb102_n, 1408
.set nb102_nn1, 1412
.set nb102_nri, 1416
.set nb102_nouter, 1420
.set nb102_ninner, 1424
.set nb102_salign, 1428
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
        movl %eax,nb102_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb102_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb102_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb102_nouter(%esp)
        movl %eax,nb102_ninner(%esp)


        ## assume we have at least one i particle - start directly 
        movl  nb102_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb102_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movss 4(%edx,%ebx,4),%xmm5
        movl nb102_p_facel(%ebp),%esi
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
        movaps %xmm3,nb102_qqOO(%esp)
        movaps %xmm4,nb102_qqOH(%esp)
        movaps %xmm5,nb102_qqHH(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb102_half(%esp)
        movss nb102_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb102_half(%esp)
        movaps %xmm3,nb102_three(%esp)

_nb_kernel102_ia32_sse.nb102_threadloop: 
        movl  nb102_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel102_ia32_sse.nb102_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel102_ia32_sse.nb102_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb102_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb102_n(%esp)
        movl %ebx,nb102_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel102_ia32_sse.nb102_outerstart
        jmp _nb_kernel102_ia32_sse.nb102_end

_nb_kernel102_ia32_sse.nb102_outerstart: 
        ## ebx contains number of outer iterations
        addl nb102_nouter(%esp),%ebx
        movl %ebx,nb102_nouter(%esp)

_nb_kernel102_ia32_sse.nb102_outer: 
        movl  nb102_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb102_is3(%esp)      ## store is3 

        movl  nb102_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb102_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb102_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb102_ii3(%esp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss (%eax,%ebx,4),%xmm3
        addss 4(%eax,%ebx,4),%xmm4
        addss 8(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb102_ixO(%esp)
        movaps %xmm4,nb102_iyO(%esp)
        movaps %xmm5,nb102_izO(%esp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 12(%eax,%ebx,4),%xmm0
        addss 16(%eax,%ebx,4),%xmm1
        addss 20(%eax,%ebx,4),%xmm2
        addss 24(%eax,%ebx,4),%xmm3
        addss 28(%eax,%ebx,4),%xmm4
        addss 32(%eax,%ebx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb102_ixH1(%esp)
        movaps %xmm1,nb102_iyH1(%esp)
        movaps %xmm2,nb102_izH1(%esp)
        movaps %xmm3,nb102_ixH2(%esp)
        movaps %xmm4,nb102_iyH2(%esp)
        movaps %xmm5,nb102_izH2(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb102_vctot(%esp)
        movaps %xmm4,nb102_fixO(%esp)
        movaps %xmm4,nb102_fiyO(%esp)
        movaps %xmm4,nb102_fizO(%esp)
        movaps %xmm4,nb102_fixH1(%esp)
        movaps %xmm4,nb102_fiyH1(%esp)
        movaps %xmm4,nb102_fizH1(%esp)
        movaps %xmm4,nb102_fixH2(%esp)
        movaps %xmm4,nb102_fiyH2(%esp)
        movaps %xmm4,nb102_fizH2(%esp)

        movl  nb102_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                         ## number of innerloop atoms 

        movl  nb102_pos(%ebp),%esi
        movl  nb102_faction(%ebp),%edi
        movl  nb102_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb102_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb102_ninner(%esp),%ecx
        movl  %ecx,nb102_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb102_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel102_ia32_sse.nb102_unroll_loop
        jmp   _nb_kernel102_ia32_sse.nb102_single_check
_nb_kernel102_ia32_sse.nb102_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb102_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 

        addl $16,nb102_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb102_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx     ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move j coordinates to local temp variables 
        movlps (%esi,%eax,4),%xmm2
        movlps 12(%esi,%eax,4),%xmm3
        movlps 24(%esi,%eax,4),%xmm4

        movlps (%esi,%ebx,4),%xmm5
        movlps 12(%esi,%ebx,4),%xmm6
        movlps 24(%esi,%ebx,4),%xmm7

        movhps (%esi,%ecx,4),%xmm2
        movhps 12(%esi,%ecx,4),%xmm3
        movhps 24(%esi,%ecx,4),%xmm4

        movhps (%esi,%edx,4),%xmm5
        movhps 12(%esi,%edx,4),%xmm6
        movhps 24(%esi,%edx,4),%xmm7

        ## current state:       
        ## xmm2= jxOa  jyOa  jxOc  jyOc 
        ## xmm3= jxH1a jyH1a jxH1c jyH1c 
        ## xmm4= jxH2a jyH2a jxH2c jyH2c 
        ## xmm5= jxOb  jyOb  jxOd  jyOd 
        ## xmm6= jxH1b jyH1b jxH1d jyH1d 
        ## xmm7= jxH2b jyH2b jxH2d jyH2d 

        movaps %xmm2,%xmm0
        movaps %xmm3,%xmm1
        unpcklps %xmm5,%xmm0    ## xmm0= jxOa  jxOb  jyOa  jyOb 
        unpcklps %xmm6,%xmm1    ## xmm1= jxH1a jxH1b jyH1a jyH1b 
        unpckhps %xmm5,%xmm2    ## xmm2= jxOc  jxOd  jyOc  jyOd 
        unpckhps %xmm6,%xmm3    ## xmm3= jxH1c jxH1d jyH1c jyH1d  
        movaps %xmm4,%xmm5
        movaps   %xmm0,%xmm6
        unpcklps %xmm7,%xmm4    ## xmm4= jxH2a jxH2b jyH2a jyH2b                
        unpckhps %xmm7,%xmm5    ## xmm5= jxH2c jxH2d jyH2c jyH2d 
        movaps   %xmm1,%xmm7
        movlhps  %xmm2,%xmm0    ## xmm0= jxOa  jxOb  jxOc  jxOd  
        movaps %xmm0,nb102_jxO(%esp)
        movhlps  %xmm6,%xmm2    ## xmm2= jyOa  jyOb  jyOc  jyOd 
        movaps %xmm2,nb102_jyO(%esp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,nb102_jxH1(%esp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,nb102_jyH1(%esp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,nb102_jxH2(%esp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,nb102_jyH2(%esp)

        movss  8(%esi,%eax,4),%xmm0
        movss  20(%esi,%eax,4),%xmm1
        movss  32(%esi,%eax,4),%xmm2

        movss  8(%esi,%ecx,4),%xmm3
        movss  20(%esi,%ecx,4),%xmm4
        movss  32(%esi,%ecx,4),%xmm5

        movhps 4(%esi,%ebx,4),%xmm0
        movhps 16(%esi,%ebx,4),%xmm1
        movhps 28(%esi,%ebx,4),%xmm2

        movhps 4(%esi,%edx,4),%xmm3
        movhps 16(%esi,%edx,4),%xmm4
        movhps 28(%esi,%edx,4),%xmm5

        shufps $204,%xmm3,%xmm0 ## constant 11001100
        shufps $204,%xmm4,%xmm1 ## constant 11001100
        shufps $204,%xmm5,%xmm2 ## constant 11001100
        movaps %xmm0,nb102_jzO(%esp)
        movaps %xmm1,nb102_jzH1(%esp)
        movaps %xmm2,nb102_jzH2(%esp)

        movaps nb102_ixO(%esp),%xmm0
        movaps nb102_iyO(%esp),%xmm1
        movaps nb102_izO(%esp),%xmm2
        movaps nb102_ixO(%esp),%xmm3
        movaps nb102_iyO(%esp),%xmm4
        movaps nb102_izO(%esp),%xmm5
        subps  nb102_jxO(%esp),%xmm0
        subps  nb102_jyO(%esp),%xmm1
        subps  nb102_jzO(%esp),%xmm2
        subps  nb102_jxH1(%esp),%xmm3
        subps  nb102_jyH1(%esp),%xmm4
        subps  nb102_jzH1(%esp),%xmm5
        movaps %xmm0,nb102_dxOO(%esp)
        movaps %xmm1,nb102_dyOO(%esp)
        movaps %xmm2,nb102_dzOO(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb102_dxOH1(%esp)
        movaps %xmm4,nb102_dyOH1(%esp)
        movaps %xmm5,nb102_dzOH1(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb102_rsqOO(%esp)
        movaps %xmm3,nb102_rsqOH1(%esp)

        movaps nb102_ixO(%esp),%xmm0
        movaps nb102_iyO(%esp),%xmm1
        movaps nb102_izO(%esp),%xmm2
        movaps nb102_ixH1(%esp),%xmm3
        movaps nb102_iyH1(%esp),%xmm4
        movaps nb102_izH1(%esp),%xmm5
        subps  nb102_jxH2(%esp),%xmm0
        subps  nb102_jyH2(%esp),%xmm1
        subps  nb102_jzH2(%esp),%xmm2
        subps  nb102_jxO(%esp),%xmm3
        subps  nb102_jyO(%esp),%xmm4
        subps  nb102_jzO(%esp),%xmm5
        movaps %xmm0,nb102_dxOH2(%esp)
        movaps %xmm1,nb102_dyOH2(%esp)
        movaps %xmm2,nb102_dzOH2(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb102_dxH1O(%esp)
        movaps %xmm4,nb102_dyH1O(%esp)
        movaps %xmm5,nb102_dzH1O(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb102_rsqOH2(%esp)
        movaps %xmm3,nb102_rsqH1O(%esp)

        movaps nb102_ixH1(%esp),%xmm0
        movaps nb102_iyH1(%esp),%xmm1
        movaps nb102_izH1(%esp),%xmm2
        movaps nb102_ixH1(%esp),%xmm3
        movaps nb102_iyH1(%esp),%xmm4
        movaps nb102_izH1(%esp),%xmm5
        subps  nb102_jxH1(%esp),%xmm0
        subps  nb102_jyH1(%esp),%xmm1
        subps  nb102_jzH1(%esp),%xmm2
        subps  nb102_jxH2(%esp),%xmm3
        subps  nb102_jyH2(%esp),%xmm4
        subps  nb102_jzH2(%esp),%xmm5
        movaps %xmm0,nb102_dxH1H1(%esp)
        movaps %xmm1,nb102_dyH1H1(%esp)
        movaps %xmm2,nb102_dzH1H1(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb102_dxH1H2(%esp)
        movaps %xmm4,nb102_dyH1H2(%esp)
        movaps %xmm5,nb102_dzH1H2(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb102_rsqH1H1(%esp)
        movaps %xmm3,nb102_rsqH1H2(%esp)

        movaps nb102_ixH2(%esp),%xmm0
        movaps nb102_iyH2(%esp),%xmm1
        movaps nb102_izH2(%esp),%xmm2
        movaps nb102_ixH2(%esp),%xmm3
        movaps nb102_iyH2(%esp),%xmm4
        movaps nb102_izH2(%esp),%xmm5
        subps  nb102_jxO(%esp),%xmm0
        subps  nb102_jyO(%esp),%xmm1
        subps  nb102_jzO(%esp),%xmm2
        subps  nb102_jxH1(%esp),%xmm3
        subps  nb102_jyH1(%esp),%xmm4
        subps  nb102_jzH1(%esp),%xmm5
        movaps %xmm0,nb102_dxH2O(%esp)
        movaps %xmm1,nb102_dyH2O(%esp)
        movaps %xmm2,nb102_dzH2O(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb102_dxH2H1(%esp)
        movaps %xmm4,nb102_dyH2H1(%esp)
        movaps %xmm5,nb102_dzH2H1(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,nb102_rsqH2O(%esp)
        movaps %xmm4,nb102_rsqH2H1(%esp)

        movaps nb102_ixH2(%esp),%xmm0
        movaps nb102_iyH2(%esp),%xmm1
        movaps nb102_izH2(%esp),%xmm2
        subps  nb102_jxH2(%esp),%xmm0
        subps  nb102_jyH2(%esp),%xmm1
        subps  nb102_jzH2(%esp),%xmm2
        movaps %xmm0,nb102_dxH2H2(%esp)
        movaps %xmm1,nb102_dyH2H2(%esp)
        movaps %xmm2,nb102_dzH2H2(%esp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,nb102_rsqH2H2(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb102_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb102_half(%esp),%xmm3   ## rinvH2H2 
        mulps   nb102_half(%esp),%xmm7   ## rinvH2H1 
        movaps  %xmm3,nb102_rinvH2H2(%esp)
        movaps  %xmm7,nb102_rinvH2H1(%esp)

        rsqrtps nb102_rsqOO(%esp),%xmm1
        rsqrtps nb102_rsqOH1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb102_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb102_rsqOO(%esp),%xmm1
        mulps   nb102_rsqOH1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb102_half(%esp),%xmm3
        mulps   nb102_half(%esp),%xmm7
        movaps  %xmm3,nb102_rinvOO(%esp)
        movaps  %xmm7,nb102_rinvOH1(%esp)

        rsqrtps nb102_rsqOH2(%esp),%xmm1
        rsqrtps nb102_rsqH1O(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb102_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb102_rsqOH2(%esp),%xmm1
        mulps   nb102_rsqH1O(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb102_half(%esp),%xmm3
        mulps   nb102_half(%esp),%xmm7
        movaps  %xmm3,nb102_rinvOH2(%esp)
        movaps  %xmm7,nb102_rinvH1O(%esp)

        rsqrtps nb102_rsqH1H1(%esp),%xmm1
        rsqrtps nb102_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb102_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb102_rsqH1H1(%esp),%xmm1
        mulps   nb102_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb102_half(%esp),%xmm3
        mulps   nb102_half(%esp),%xmm7
        movaps  %xmm3,nb102_rinvH1H1(%esp)
        movaps  %xmm7,nb102_rinvH1H2(%esp)

        rsqrtps nb102_rsqH2O(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb102_three(%esp),%xmm3
        mulps   nb102_rsqH2O(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb102_half(%esp),%xmm3
        movaps  %xmm3,nb102_rinvH2O(%esp)

        ## start with OO interaction 
        movaps nb102_rinvOO(%esp),%xmm0
        movaps %xmm0,%xmm7
        mulps  %xmm0,%xmm0
        mulps  nb102_qqOO(%esp),%xmm7
        mulps  %xmm7,%xmm0
        addps  nb102_vctot(%esp),%xmm7
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb102_dxOO(%esp),%xmm0
        mulps nb102_dyOO(%esp),%xmm1
        mulps nb102_dzOO(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb102_fixO(%esp),%xmm0
        addps nb102_fiyO(%esp),%xmm1
        addps nb102_fizO(%esp),%xmm2
        movaps %xmm3,nb102_fjxO(%esp)
        movaps %xmm4,nb102_fjyO(%esp)
        movaps %xmm5,nb102_fjzO(%esp)
        movaps %xmm0,nb102_fixO(%esp)
        movaps %xmm1,nb102_fiyO(%esp)
        movaps %xmm2,nb102_fizO(%esp)

        ## O-H1 interaction 
        movaps nb102_rinvOH1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb102_qqOH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fsOH1  
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb102_dxOH1(%esp),%xmm0
        mulps nb102_dyOH1(%esp),%xmm1
        mulps nb102_dzOH1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb102_fixO(%esp),%xmm0
        addps nb102_fiyO(%esp),%xmm1
        addps nb102_fizO(%esp),%xmm2
        movaps %xmm3,nb102_fjxH1(%esp)
        movaps %xmm4,nb102_fjyH1(%esp)
        movaps %xmm5,nb102_fjzH1(%esp)
        movaps %xmm0,nb102_fixO(%esp)
        movaps %xmm1,nb102_fiyO(%esp)
        movaps %xmm2,nb102_fizO(%esp)

        ## O-H2 interaction  
        movaps nb102_rinvOH2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb102_qqOH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fsOH2  
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb102_dxOH2(%esp),%xmm0
        mulps nb102_dyOH2(%esp),%xmm1
        mulps nb102_dzOH2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb102_fixO(%esp),%xmm0
        addps nb102_fiyO(%esp),%xmm1
        addps nb102_fizO(%esp),%xmm2
        movaps %xmm3,nb102_fjxH2(%esp)
        movaps %xmm4,nb102_fjyH2(%esp)
        movaps %xmm5,nb102_fjzH2(%esp)
        movaps %xmm0,nb102_fixO(%esp)
        movaps %xmm1,nb102_fiyO(%esp)
        movaps %xmm2,nb102_fizO(%esp)

        ## H1-O interaction 
        movaps nb102_rinvH1O(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb102_qqOH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fsH1O 
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2
        movaps nb102_fjxO(%esp),%xmm3
        movaps nb102_fjyO(%esp),%xmm4
        movaps nb102_fjzO(%esp),%xmm5
        mulps nb102_dxH1O(%esp),%xmm0
        mulps nb102_dyH1O(%esp),%xmm1
        mulps nb102_dzH1O(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb102_fixH1(%esp),%xmm0
        addps nb102_fiyH1(%esp),%xmm1
        addps nb102_fizH1(%esp),%xmm2
        movaps %xmm3,nb102_fjxO(%esp)
        movaps %xmm4,nb102_fjyO(%esp)
        movaps %xmm5,nb102_fjzO(%esp)
        movaps %xmm0,nb102_fixH1(%esp)
        movaps %xmm1,nb102_fiyH1(%esp)
        movaps %xmm2,nb102_fizH1(%esp)

        ## H1-H1 interaction 
        movaps nb102_rinvH1H1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb102_qqHH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fsH1H1 
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2
        movaps nb102_fjxH1(%esp),%xmm3
        movaps nb102_fjyH1(%esp),%xmm4
        movaps nb102_fjzH1(%esp),%xmm5
        mulps nb102_dxH1H1(%esp),%xmm0
        mulps nb102_dyH1H1(%esp),%xmm1
        mulps nb102_dzH1H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb102_fixH1(%esp),%xmm0
        addps nb102_fiyH1(%esp),%xmm1
        addps nb102_fizH1(%esp),%xmm2
        movaps %xmm3,nb102_fjxH1(%esp)
        movaps %xmm4,nb102_fjyH1(%esp)
        movaps %xmm5,nb102_fjzH1(%esp)
        movaps %xmm0,nb102_fixH1(%esp)
        movaps %xmm1,nb102_fiyH1(%esp)
        movaps %xmm2,nb102_fizH1(%esp)

        ## H1-H2 interaction 
        movaps nb102_rinvH1H2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb102_qqHH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fsOH2  
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2
        movaps nb102_fjxH2(%esp),%xmm3
        movaps nb102_fjyH2(%esp),%xmm4
        movaps nb102_fjzH2(%esp),%xmm5
        mulps nb102_dxH1H2(%esp),%xmm0
        mulps nb102_dyH1H2(%esp),%xmm1
        mulps nb102_dzH1H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb102_fixH1(%esp),%xmm0
        addps nb102_fiyH1(%esp),%xmm1
        addps nb102_fizH1(%esp),%xmm2
        movaps %xmm3,nb102_fjxH2(%esp)
        movaps %xmm4,nb102_fjyH2(%esp)
        movaps %xmm5,nb102_fjzH2(%esp)
        movaps %xmm0,nb102_fixH1(%esp)
        movaps %xmm1,nb102_fiyH1(%esp)
        movaps %xmm2,nb102_fizH1(%esp)

        ## H2-O interaction 
        movaps nb102_rinvH2O(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb102_qqOH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fsH2O 
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2
        movaps nb102_fjxO(%esp),%xmm3
        movaps nb102_fjyO(%esp),%xmm4
        movaps nb102_fjzO(%esp),%xmm5
        mulps nb102_dxH2O(%esp),%xmm0
        mulps nb102_dyH2O(%esp),%xmm1
        mulps nb102_dzH2O(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb102_fixH2(%esp),%xmm0
        addps nb102_fiyH2(%esp),%xmm1
        addps nb102_fizH2(%esp),%xmm2
        movaps %xmm3,nb102_fjxO(%esp)
        movaps %xmm4,nb102_fjyO(%esp)
        movaps %xmm5,nb102_fjzO(%esp)
        movaps %xmm0,nb102_fixH2(%esp)
        movaps %xmm1,nb102_fiyH2(%esp)
        movaps %xmm2,nb102_fizH2(%esp)

        ## H2-H1 interaction 
        movaps nb102_rinvH2H1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb102_qqHH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fsH2H1 
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2
        movaps nb102_fjxH1(%esp),%xmm3
        movaps nb102_fjyH1(%esp),%xmm4
        movaps nb102_fjzH1(%esp),%xmm5
        mulps nb102_dxH2H1(%esp),%xmm0
        mulps nb102_dyH2H1(%esp),%xmm1
        mulps nb102_dzH2H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb102_fixH2(%esp),%xmm0
        addps nb102_fiyH2(%esp),%xmm1
        addps nb102_fizH2(%esp),%xmm2
        movaps %xmm3,nb102_fjxH1(%esp)
        movaps %xmm4,nb102_fjyH1(%esp)
        movaps %xmm5,nb102_fjzH1(%esp)
        movaps %xmm0,nb102_fixH2(%esp)
        movaps %xmm1,nb102_fiyH2(%esp)
        movaps %xmm2,nb102_fizH2(%esp)

        ## H2-H2 interaction 
        movaps nb102_rinvH2H2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps nb102_qqHH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fsH2H2 
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm7,nb102_vctot(%esp)
        movaps %xmm0,%xmm2
        movaps nb102_fjxH2(%esp),%xmm3
        movaps nb102_fjyH2(%esp),%xmm4
        movaps nb102_fjzH2(%esp),%xmm5
        mulps nb102_dxH2H2(%esp),%xmm0
        mulps nb102_dyH2H2(%esp),%xmm1
        mulps nb102_dzH2H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb102_fixH2(%esp),%xmm0
        addps nb102_fiyH2(%esp),%xmm1
        addps nb102_fizH2(%esp),%xmm2
        movaps %xmm3,nb102_fjxH2(%esp)
        movaps %xmm4,nb102_fjyH2(%esp)
        movaps %xmm5,nb102_fjzH2(%esp)
        movaps %xmm0,nb102_fixH2(%esp)
        movaps %xmm1,nb102_fiyH2(%esp)
        movaps %xmm2,nb102_fizH2(%esp)

        movl nb102_faction(%ebp),%edi

        ## Did all interactions - now update j forces 
        ## At this stage forces are still on the stack, in positions:
        ## fjxO, fjyO, fjzO, ... , fjzH2.
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


        ## 4 j waters with three atoms each - first do Oxygen X & Y forces for 4 j particles 
        movaps nb102_fjxO(%esp),%xmm0   ## xmm0= fjxOa  fjxOb  fjxOc  fjxOd 
        movaps nb102_fjyO(%esp),%xmm2   ## xmm1= fjyOa  fjyOb  fjyOc  fjyOd
        movlps (%edi,%eax,4),%xmm3
        movlps (%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0       ## xmm0= fjxOa  fjyOa  fjxOb  fjyOb
        unpckhps %xmm2,%xmm1       ## xmm1= fjxOc  fjyOc  fjxOd  fjyOd
        movhps (%edi,%ebx,4),%xmm3
        movhps (%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,(%edi,%eax,4)
        movlps %xmm4,(%edi,%ecx,4)
        movhps %xmm3,(%edi,%ebx,4)
        movhps %xmm4,(%edi,%edx,4)

        ## Oxygen Z & first hydrogen X forces for 4 j particles 
        movaps nb102_fjzO(%esp),%xmm0    ## xmm0= fjzOa   fjzOb   fjzOc   fjzOd 
        movaps nb102_fjxH1(%esp),%xmm2   ## xmm1= fjxH1a  fjxH1b  fjxH1c  fjxH1d
        movlps 8(%edi,%eax,4),%xmm3
        movlps 8(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0       ## xmm0= fjzOa  fjxH1a  fjzOb  fjxH1b
        unpckhps %xmm2,%xmm1       ## xmm1= fjzOc  fjxH1c  fjzOd  fjxH1d
        movhps 8(%edi,%ebx,4),%xmm3
        movhps 8(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,8(%edi,%eax,4)
        movlps %xmm4,8(%edi,%ecx,4)
        movhps %xmm3,8(%edi,%ebx,4)
        movhps %xmm4,8(%edi,%edx,4)


        ## First hydrogen Y & Z forces for 4 j particles 
        movaps nb102_fjyH1(%esp),%xmm0    ## xmm0= fjyH1a  fjyH1b  fjyH1c  fjyH1d 
        movaps nb102_fjzH1(%esp),%xmm2   ## xmm1= fjzH1a  fjzH1b  fjzH1c  fjzH1d
        movlps 16(%edi,%eax,4),%xmm3
        movlps 16(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0            ## xmm0= fjyH1a  fjzH1a  fjyH1b  fjzH1b
        unpckhps %xmm2,%xmm1            ## xmm1= fjyH1c  fjzH1c  fjyH1d  fjzH1d
        movhps 16(%edi,%ebx,4),%xmm3
        movhps 16(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,16(%edi,%eax,4)
        movlps %xmm4,16(%edi,%ecx,4)
        movhps %xmm3,16(%edi,%ebx,4)
        movhps %xmm4,16(%edi,%edx,4)


        ## Second hydrogen X & Y forces for 4 j particles 
        movaps nb102_fjxH2(%esp),%xmm0    ## xmm0= fjxH2a  fjxH2b  fjxH2c  fjxH2d 
        movaps nb102_fjyH2(%esp),%xmm2   ## xmm1= fjyH2a  fjyH2b  fjyH2c  fjyH2d
        movlps 24(%edi,%eax,4),%xmm3
        movlps 24(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0            ## xmm0= fjxH2a  fjyH2a  fjxH2b  fjyH2b
        unpckhps %xmm2,%xmm1            ## xmm1= fjxH2c  fjyH2c  fjxH2d  fjyH2d
        movhps 24(%edi,%ebx,4),%xmm3
        movhps 24(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,24(%edi,%eax,4)
        movlps %xmm4,24(%edi,%ecx,4)
        movhps %xmm3,24(%edi,%ebx,4)
        movhps %xmm4,24(%edi,%edx,4)


        ## Second hydrogen Z forces for 4 j particles 
        ## Just load the four Z coords into one reg. each
        movss 32(%edi,%eax,4),%xmm4
        movss 32(%edi,%ebx,4),%xmm5
        movss 32(%edi,%ecx,4),%xmm6
        movss 32(%edi,%edx,4),%xmm7
        ## add what we have on the stack
        addss nb102_fjzH2(%esp),%xmm4
        addss nb102_fjzH2b(%esp),%xmm5
        addss nb102_fjzH2c(%esp),%xmm6
        addss nb102_fjzH2d(%esp),%xmm7
        ## store back
        movss %xmm4,32(%edi,%eax,4)
        movss %xmm5,32(%edi,%ebx,4)
        movss %xmm6,32(%edi,%ecx,4)
        movss %xmm7,32(%edi,%edx,4)

        ## should we do one more iteration? 
        subl $4,nb102_innerk(%esp)
        jl    _nb_kernel102_ia32_sse.nb102_single_check
        jmp   _nb_kernel102_ia32_sse.nb102_unroll_loop
_nb_kernel102_ia32_sse.nb102_single_check: 
        addl $4,nb102_innerk(%esp)
        jnz   _nb_kernel102_ia32_sse.nb102_single_loop
        jmp   _nb_kernel102_ia32_sse.nb102_updateouterdata
_nb_kernel102_ia32_sse.nb102_single_loop: 
        movl  nb102_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb102_innerjjnr(%esp)

        movl nb102_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        xorps %xmm3,%xmm3
        xorps %xmm4,%xmm4
        xorps %xmm5,%xmm5

        movss (%esi,%eax,4),%xmm3               ## jxO  -  -  -
        movss 4(%esi,%eax,4),%xmm4              ## jyO  -  -  -
        movss 8(%esi,%eax,4),%xmm5              ## jzO  -  -  -  

        movlps 12(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1   -    -
        movss  20(%esi,%eax,4),%xmm7            ## xmm7 = jzH1   -    -    - 
        movhps 24(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1 jxH2 jyH2
        movss  32(%esi,%eax,4),%xmm2            ## xmm2 = jzH2   -    -    -

        ## have all coords, time for some shuffling.

        shufps $216,%xmm6,%xmm6 ## constant 11011000     ;# xmm6 = jxH1 jxH2 jyH1 jyH2 
        unpcklps %xmm2,%xmm7                    ## xmm7 = jzH1 jzH2   -    -
        movaps  nb102_ixO(%esp),%xmm0
        movaps  nb102_iyO(%esp),%xmm1
        movaps  nb102_izO(%esp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxO   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## constant 11100100    ;# xmm4 = jyO   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## constant 01000100    ;# xmm5 = jzO   0   jzH1 jzH2

        ## store all j coordinates in jO  
        movaps %xmm3,nb102_jxO(%esp)
        movaps %xmm4,nb102_jyO(%esp)
        movaps %xmm5,nb102_jzO(%esp)
        subps  %xmm3,%xmm0
        subps  %xmm4,%xmm1
        subps  %xmm5,%xmm2
        movaps %xmm0,nb102_dxOO(%esp)
        movaps %xmm1,nb102_dyOO(%esp)
        movaps %xmm2,nb102_dzOO(%esp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb102_three(%esp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb102_half(%esp),%xmm3   ## rinv iO - j water 

        xorps   %xmm1,%xmm1
        movaps  %xmm3,%xmm0
        xorps   %xmm4,%xmm4
        mulps   %xmm0,%xmm0     ## xmm0=rinvsq 
        ## fetch charges to xmm4 (temporary) 
        movss   nb102_qqOO(%esp),%xmm4

        movhps  nb102_qqOH(%esp),%xmm4

        mulps   %xmm4,%xmm3     ## xmm3=vcoul 
        mulps   %xmm3,%xmm0     ## total fscal 
        addps   nb102_vctot(%esp),%xmm3
        movaps  %xmm3,nb102_vctot(%esp)

        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2
        mulps   nb102_dxOO(%esp),%xmm0
        mulps   nb102_dyOO(%esp),%xmm1
        mulps   nb102_dzOO(%esp),%xmm2
        ## initial update for j forces 
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb102_fjxO(%esp)
        movaps  %xmm4,nb102_fjyO(%esp)
        movaps  %xmm5,nb102_fjzO(%esp)
        addps   nb102_fixO(%esp),%xmm0
        addps   nb102_fiyO(%esp),%xmm1
        addps   nb102_fizO(%esp),%xmm2
        movaps  %xmm0,nb102_fixO(%esp)
        movaps  %xmm1,nb102_fiyO(%esp)
        movaps  %xmm2,nb102_fizO(%esp)


        ## done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb102_ixH1(%esp),%xmm0
        movaps  nb102_iyH1(%esp),%xmm1
        movaps  nb102_izH1(%esp),%xmm2
        movaps  nb102_ixH2(%esp),%xmm3
        movaps  nb102_iyH2(%esp),%xmm4
        movaps  nb102_izH2(%esp),%xmm5
        subps   nb102_jxO(%esp),%xmm0
        subps   nb102_jyO(%esp),%xmm1
        subps   nb102_jzO(%esp),%xmm2
        subps   nb102_jxO(%esp),%xmm3
        subps   nb102_jyO(%esp),%xmm4
        subps   nb102_jzO(%esp),%xmm5
        movaps %xmm0,nb102_dxH1O(%esp)
        movaps %xmm1,nb102_dyH1O(%esp)
        movaps %xmm2,nb102_dzH1O(%esp)
        movaps %xmm3,nb102_dxH2O(%esp)
        movaps %xmm4,nb102_dyH2O(%esp)
        movaps %xmm5,nb102_dzH2O(%esp)
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
        movaps  nb102_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb102_half(%esp),%xmm3   ## rinv H1 - j water 
        mulps   nb102_half(%esp),%xmm7   ## rinv H2 - j water  

        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        ## do coulomb interaction 
        movaps  %xmm3,%xmm0
        movss   nb102_qqOH(%esp),%xmm6
        movaps  %xmm7,%xmm4
        movhps  nb102_qqHH(%esp),%xmm6
        mulps   %xmm0,%xmm0     ## rinvsq 
        mulps   %xmm4,%xmm4     ## rinvsq 
        mulps   %xmm6,%xmm3     ## vcoul 
        mulps   %xmm6,%xmm7     ## vcoul 
        movaps  %xmm3,%xmm2
        addps   %xmm7,%xmm2     ## total vcoul 
        mulps   %xmm3,%xmm0     ## fscal 

        addps   nb102_vctot(%esp),%xmm2
        mulps   %xmm4,%xmm7     ## fscal 
        movaps  %xmm2,nb102_vctot(%esp)
        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2
        mulps   nb102_dxH1O(%esp),%xmm0
        mulps   nb102_dyH1O(%esp),%xmm1
        mulps   nb102_dzH1O(%esp),%xmm2
        ## update forces H1 - j water 
        movaps  nb102_fjxO(%esp),%xmm3
        movaps  nb102_fjyO(%esp),%xmm4
        movaps  nb102_fjzO(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb102_fjxO(%esp)
        movaps  %xmm4,nb102_fjyO(%esp)
        movaps  %xmm5,nb102_fjzO(%esp)
        addps   nb102_fixH1(%esp),%xmm0
        addps   nb102_fiyH1(%esp),%xmm1
        addps   nb102_fizH1(%esp),%xmm2
        movaps  %xmm0,nb102_fixH1(%esp)
        movaps  %xmm1,nb102_fiyH1(%esp)
        movaps  %xmm2,nb102_fizH1(%esp)
        ## do forces H2 - j water 
        movaps %xmm7,%xmm0
        movaps %xmm7,%xmm1
        movaps %xmm7,%xmm2
        mulps   nb102_dxH2O(%esp),%xmm0
        mulps   nb102_dyH2O(%esp),%xmm1
        mulps   nb102_dzH2O(%esp),%xmm2
        movaps  nb102_fjxO(%esp),%xmm3
        movaps  nb102_fjyO(%esp),%xmm4
        movaps  nb102_fjzO(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movl    nb102_faction(%ebp),%esi
        movaps  %xmm3,nb102_fjxO(%esp)
        movaps  %xmm4,nb102_fjyO(%esp)
        movaps  %xmm5,nb102_fjzO(%esp)
        addps   nb102_fixH2(%esp),%xmm0
        addps   nb102_fiyH2(%esp),%xmm1
        addps   nb102_fizH2(%esp),%xmm2
        movaps  %xmm0,nb102_fixH2(%esp)
        movaps  %xmm1,nb102_fiyH2(%esp)
        movaps  %xmm2,nb102_fizH2(%esp)

        ## update j water forces from local variables 
        movlps  (%esi,%eax,4),%xmm0
        movlps  12(%esi,%eax,4),%xmm1
        movhps  24(%esi,%eax,4),%xmm1
        movaps  nb102_fjxO(%esp),%xmm3
        movaps  nb102_fjyO(%esp),%xmm4
        movaps  nb102_fjzO(%esp),%xmm5
        movaps  %xmm5,%xmm6
        movaps  %xmm5,%xmm7
        shufps $2,%xmm6,%xmm6 ## constant 00000010
        shufps $3,%xmm7,%xmm7 ## constant 00000011
        addss   8(%esi,%eax,4),%xmm5
        addss   20(%esi,%eax,4),%xmm6
        addss   32(%esi,%eax,4),%xmm7
        movss   %xmm5,8(%esi,%eax,4)
        movss   %xmm6,20(%esi,%eax,4)
        movss   %xmm7,32(%esi,%eax,4)
        movaps   %xmm3,%xmm5
        unpcklps %xmm4,%xmm3
        unpckhps %xmm4,%xmm5
        addps    %xmm3,%xmm0
        addps    %xmm5,%xmm1
        movlps  %xmm0,(%esi,%eax,4)
        movlps  %xmm1,12(%esi,%eax,4)
        movhps  %xmm1,24(%esi,%eax,4)

        decl  nb102_innerk(%esp)
        jz    _nb_kernel102_ia32_sse.nb102_updateouterdata
        jmp   _nb_kernel102_ia32_sse.nb102_single_loop
_nb_kernel102_ia32_sse.nb102_updateouterdata: 
        movl  nb102_ii3(%esp),%ecx
        movl  nb102_faction(%ebp),%edi
        movl  nb102_fshift(%ebp),%esi
        movl  nb102_is3(%esp),%edx

        ## accumulate Oi forces in xmm0, xmm1, xmm2 
        movaps nb102_fixO(%esp),%xmm0
        movaps nb102_fiyO(%esp),%xmm1
        movaps nb102_fizO(%esp),%xmm2

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
        movaps nb102_fixH1(%esp),%xmm0
        movaps nb102_fiyH1(%esp),%xmm1
        movaps nb102_fizH1(%esp),%xmm2

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
        movaps nb102_fixH2(%esp),%xmm0
        movaps nb102_fiyH2(%esp),%xmm1
        movaps nb102_fizH2(%esp),%xmm2

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

        ## increment fshift force  
        movlps  (%esi,%edx,4),%xmm3
        movss  8(%esi,%edx,4),%xmm4
        addps  %xmm6,%xmm3
        addss  %xmm7,%xmm4
        movlps  %xmm3,(%esi,%edx,4)
        movss  %xmm4,8(%esi,%edx,4)

        ## get n from stack
        movl nb102_n(%esp),%esi
        ## get group index for i particle 
        movl  nb102_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb102_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb102_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb102_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel102_ia32_sse.nb102_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb102_n(%esp)
        jmp _nb_kernel102_ia32_sse.nb102_outer
_nb_kernel102_ia32_sse.nb102_outerend: 
        ## check if more outer neighborlists remain
        movl  nb102_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel102_ia32_sse.nb102_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel102_ia32_sse.nb102_threadloop
_nb_kernel102_ia32_sse.nb102_end: 
        emms

        movl nb102_nouter(%esp),%eax
        movl nb102_ninner(%esp),%ebx
        movl nb102_outeriter(%ebp),%ecx
        movl nb102_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb102_salign(%esp),%eax
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





.globl nb_kernel102nf_ia32_sse
.globl _nb_kernel102nf_ia32_sse
nb_kernel102nf_ia32_sse:        
_nb_kernel102nf_ia32_sse:       
.set nb102nf_p_nri, 8
.set nb102nf_iinr, 12
.set nb102nf_jindex, 16
.set nb102nf_jjnr, 20
.set nb102nf_shift, 24
.set nb102nf_shiftvec, 28
.set nb102nf_fshift, 32
.set nb102nf_gid, 36
.set nb102nf_pos, 40
.set nb102nf_faction, 44
.set nb102nf_charge, 48
.set nb102nf_p_facel, 52
.set nb102nf_p_krf, 56
.set nb102nf_p_crf, 60
.set nb102nf_Vc, 64
.set nb102nf_type, 68
.set nb102nf_p_ntype, 72
.set nb102nf_vdwparam, 76
.set nb102nf_Vvdw, 80
.set nb102nf_p_tabscale, 84
.set nb102nf_VFtab, 88
.set nb102nf_invsqrta, 92
.set nb102nf_dvda, 96
.set nb102nf_p_gbtabscale, 100
.set nb102nf_GBtab, 104
.set nb102nf_p_nthreads, 108
.set nb102nf_count, 112
.set nb102nf_mtx, 116
.set nb102nf_outeriter, 120
.set nb102nf_inneriter, 124
.set nb102nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use         
.set nb102nf_ixO, 0
.set nb102nf_iyO, 16
.set nb102nf_izO, 32
.set nb102nf_ixH1, 48
.set nb102nf_iyH1, 64
.set nb102nf_izH1, 80
.set nb102nf_ixH2, 96
.set nb102nf_iyH2, 112
.set nb102nf_izH2, 128
.set nb102nf_jxO, 144
.set nb102nf_jyO, 160
.set nb102nf_jzO, 176
.set nb102nf_jxH1, 192
.set nb102nf_jyH1, 208
.set nb102nf_jzH1, 224
.set nb102nf_jxH2, 240
.set nb102nf_jyH2, 256
.set nb102nf_jzH2, 272
.set nb102nf_qqOO, 288
.set nb102nf_qqOH, 304
.set nb102nf_qqHH, 320
.set nb102nf_vctot, 336
.set nb102nf_half, 352
.set nb102nf_three, 368
.set nb102nf_rsqOO, 384
.set nb102nf_rsqOH1, 400
.set nb102nf_rsqOH2, 416
.set nb102nf_rsqH1O, 432
.set nb102nf_rsqH1H1, 448
.set nb102nf_rsqH1H2, 464
.set nb102nf_rsqH2O, 480
.set nb102nf_rsqH2H1, 496
.set nb102nf_rsqH2H2, 512
.set nb102nf_rinvOO, 528
.set nb102nf_rinvOH1, 544
.set nb102nf_rinvOH2, 560
.set nb102nf_rinvH1O, 576
.set nb102nf_rinvH1H1, 592
.set nb102nf_rinvH1H2, 608
.set nb102nf_rinvH2O, 624
.set nb102nf_rinvH2H1, 640
.set nb102nf_rinvH2H2, 656
.set nb102nf_is3, 672
.set nb102nf_ii3, 676
.set nb102nf_innerjjnr, 680
.set nb102nf_innerk, 684
.set nb102nf_n, 688
.set nb102nf_nn1, 692
.set nb102nf_nri, 696
.set nb102nf_nouter, 700
.set nb102nf_ninner, 704
.set nb102nf_salign, 708
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $712,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb102nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb102nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb102nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb102nf_nouter(%esp)
        movl %eax,nb102nf_ninner(%esp)


        ## assume we have at least one i particle - start directly 
        movl  nb102nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb102nf_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movss 4(%edx,%ebx,4),%xmm5
        movl nb102nf_p_facel(%ebp),%esi
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
        movaps %xmm3,nb102nf_qqOO(%esp)
        movaps %xmm4,nb102nf_qqOH(%esp)
        movaps %xmm5,nb102nf_qqHH(%esp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb102nf_half(%esp)
        movss nb102nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb102nf_half(%esp)
        movaps %xmm3,nb102nf_three(%esp)

_nb_kernel102nf_ia32_sse.nb102nf_threadloop: 
        movl  nb102nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel102nf_ia32_sse.nb102nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock  
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel102nf_ia32_sse.nb102nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb102nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb102nf_n(%esp)
        movl %ebx,nb102nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel102nf_ia32_sse.nb102nf_outerstart
        jmp _nb_kernel102nf_ia32_sse.nb102nf_end

_nb_kernel102nf_ia32_sse.nb102nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb102nf_nouter(%esp),%ebx
        movl %ebx,nb102nf_nouter(%esp)

_nb_kernel102nf_ia32_sse.nb102nf_outer: 
        movl  nb102nf_shift(%ebp),%eax          ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb102nf_is3(%esp)            ## store is3 

        movl  nb102nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb102nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb102nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb102nf_ii3(%esp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss (%eax,%ebx,4),%xmm3
        addss 4(%eax,%ebx,4),%xmm4
        addss 8(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb102nf_ixO(%esp)
        movaps %xmm4,nb102nf_iyO(%esp)
        movaps %xmm5,nb102nf_izO(%esp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 12(%eax,%ebx,4),%xmm0
        addss 16(%eax,%ebx,4),%xmm1
        addss 20(%eax,%ebx,4),%xmm2
        addss 24(%eax,%ebx,4),%xmm3
        addss 28(%eax,%ebx,4),%xmm4
        addss 32(%eax,%ebx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb102nf_ixH1(%esp)
        movaps %xmm1,nb102nf_iyH1(%esp)
        movaps %xmm2,nb102nf_izH1(%esp)
        movaps %xmm3,nb102nf_ixH2(%esp)
        movaps %xmm4,nb102nf_iyH2(%esp)
        movaps %xmm5,nb102nf_izH2(%esp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb102nf_vctot(%esp)

        movl  nb102nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                         ## number of innerloop atoms 

        movl  nb102nf_pos(%ebp),%esi
        movl  nb102nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb102nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb102nf_ninner(%esp),%ecx
        movl  %ecx,nb102nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb102nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel102nf_ia32_sse.nb102nf_unroll_loop
        jmp   _nb_kernel102nf_ia32_sse.nb102nf_single_check
_nb_kernel102nf_ia32_sse.nb102nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb102nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 

        addl $16,nb102nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb102nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx     ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move j coordinates to local temp variables 
        movlps (%esi,%eax,4),%xmm2
        movlps 12(%esi,%eax,4),%xmm3
        movlps 24(%esi,%eax,4),%xmm4

        movlps (%esi,%ebx,4),%xmm5
        movlps 12(%esi,%ebx,4),%xmm6
        movlps 24(%esi,%ebx,4),%xmm7

        movhps (%esi,%ecx,4),%xmm2
        movhps 12(%esi,%ecx,4),%xmm3
        movhps 24(%esi,%ecx,4),%xmm4

        movhps (%esi,%edx,4),%xmm5
        movhps 12(%esi,%edx,4),%xmm6
        movhps 24(%esi,%edx,4),%xmm7

        ## current state:       
        ## xmm2= jxOa  jyOa  jxOc  jyOc 
        ## xmm3= jxH1a jyH1a jxH1c jyH1c 
        ## xmm4= jxH2a jyH2a jxH2c jyH2c 
        ## xmm5= jxOb  jyOb  jxOd  jyOd 
        ## xmm6= jxH1b jyH1b jxH1d jyH1d 
        ## xmm7= jxH2b jyH2b jxH2d jyH2d 

        movaps %xmm2,%xmm0
        movaps %xmm3,%xmm1
        unpcklps %xmm5,%xmm0    ## xmm0= jxOa  jxOb  jyOa  jyOb 
        unpcklps %xmm6,%xmm1    ## xmm1= jxH1a jxH1b jyH1a jyH1b 
        unpckhps %xmm5,%xmm2    ## xmm2= jxOc  jxOd  jyOc  jyOd 
        unpckhps %xmm6,%xmm3    ## xmm3= jxH1c jxH1d jyH1c jyH1d  
        movaps %xmm4,%xmm5
        movaps   %xmm0,%xmm6
        unpcklps %xmm7,%xmm4    ## xmm4= jxH2a jxH2b jyH2a jyH2b                
        unpckhps %xmm7,%xmm5    ## xmm5= jxH2c jxH2d jyH2c jyH2d 
        movaps   %xmm1,%xmm7
        movlhps  %xmm2,%xmm0    ## xmm0= jxOa  jxOb  jxOc  jxOd  
        movaps %xmm0,nb102nf_jxO(%esp)
        movhlps  %xmm6,%xmm2    ## xmm2= jyOa  jyOb  jyOc  jyOd 
        movaps %xmm2,nb102nf_jyO(%esp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,nb102nf_jxH1(%esp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,nb102nf_jyH1(%esp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,nb102nf_jxH2(%esp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,nb102nf_jyH2(%esp)

        movss  8(%esi,%eax,4),%xmm0
        movss  20(%esi,%eax,4),%xmm1
        movss  32(%esi,%eax,4),%xmm2

        movss  8(%esi,%ecx,4),%xmm3
        movss  20(%esi,%ecx,4),%xmm4
        movss  32(%esi,%ecx,4),%xmm5

        movhps 4(%esi,%ebx,4),%xmm0
        movhps 16(%esi,%ebx,4),%xmm1
        movhps 28(%esi,%ebx,4),%xmm2

        movhps 4(%esi,%edx,4),%xmm3
        movhps 16(%esi,%edx,4),%xmm4
        movhps 28(%esi,%edx,4),%xmm5

        shufps $204,%xmm3,%xmm0 ## constant 11001100
        shufps $204,%xmm4,%xmm1 ## constant 11001100
        shufps $204,%xmm5,%xmm2 ## constant 11001100
        movaps %xmm0,nb102nf_jzO(%esp)
        movaps %xmm1,nb102nf_jzH1(%esp)
        movaps %xmm2,nb102nf_jzH2(%esp)

        movaps nb102nf_ixO(%esp),%xmm0
        movaps nb102nf_iyO(%esp),%xmm1
        movaps nb102nf_izO(%esp),%xmm2
        movaps nb102nf_ixO(%esp),%xmm3
        movaps nb102nf_iyO(%esp),%xmm4
        movaps nb102nf_izO(%esp),%xmm5
        subps  nb102nf_jxO(%esp),%xmm0
        subps  nb102nf_jyO(%esp),%xmm1
        subps  nb102nf_jzO(%esp),%xmm2
        subps  nb102nf_jxH1(%esp),%xmm3
        subps  nb102nf_jyH1(%esp),%xmm4
        subps  nb102nf_jzH1(%esp),%xmm5

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
        movaps %xmm0,nb102nf_rsqOO(%esp)
        movaps %xmm3,nb102nf_rsqOH1(%esp)

        movaps nb102nf_ixO(%esp),%xmm0
        movaps nb102nf_iyO(%esp),%xmm1
        movaps nb102nf_izO(%esp),%xmm2
        movaps nb102nf_ixH1(%esp),%xmm3
        movaps nb102nf_iyH1(%esp),%xmm4
        movaps nb102nf_izH1(%esp),%xmm5
        subps  nb102nf_jxH2(%esp),%xmm0
        subps  nb102nf_jyH2(%esp),%xmm1
        subps  nb102nf_jzH2(%esp),%xmm2
        subps  nb102nf_jxO(%esp),%xmm3
        subps  nb102nf_jyO(%esp),%xmm4
        subps  nb102nf_jzO(%esp),%xmm5

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
        movaps %xmm0,nb102nf_rsqOH2(%esp)
        movaps %xmm3,nb102nf_rsqH1O(%esp)

        movaps nb102nf_ixH1(%esp),%xmm0
        movaps nb102nf_iyH1(%esp),%xmm1
        movaps nb102nf_izH1(%esp),%xmm2
        movaps nb102nf_ixH1(%esp),%xmm3
        movaps nb102nf_iyH1(%esp),%xmm4
        movaps nb102nf_izH1(%esp),%xmm5
        subps  nb102nf_jxH1(%esp),%xmm0
        subps  nb102nf_jyH1(%esp),%xmm1
        subps  nb102nf_jzH1(%esp),%xmm2
        subps  nb102nf_jxH2(%esp),%xmm3
        subps  nb102nf_jyH2(%esp),%xmm4
        subps  nb102nf_jzH2(%esp),%xmm5

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
        movaps %xmm0,nb102nf_rsqH1H1(%esp)
        movaps %xmm3,nb102nf_rsqH1H2(%esp)

        movaps nb102nf_ixH2(%esp),%xmm0
        movaps nb102nf_iyH2(%esp),%xmm1
        movaps nb102nf_izH2(%esp),%xmm2
        movaps nb102nf_ixH2(%esp),%xmm3
        movaps nb102nf_iyH2(%esp),%xmm4
        movaps nb102nf_izH2(%esp),%xmm5
        subps  nb102nf_jxO(%esp),%xmm0
        subps  nb102nf_jyO(%esp),%xmm1
        subps  nb102nf_jzO(%esp),%xmm2
        subps  nb102nf_jxH1(%esp),%xmm3
        subps  nb102nf_jyH1(%esp),%xmm4
        subps  nb102nf_jzH1(%esp),%xmm5

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
        movaps %xmm0,nb102nf_rsqH2O(%esp)
        movaps %xmm4,nb102nf_rsqH2H1(%esp)

        movaps nb102nf_ixH2(%esp),%xmm0
        movaps nb102nf_iyH2(%esp),%xmm1
        movaps nb102nf_izH2(%esp),%xmm2
        subps  nb102nf_jxH2(%esp),%xmm0
        subps  nb102nf_jyH2(%esp),%xmm1
        subps  nb102nf_jzH2(%esp),%xmm2

        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,nb102nf_rsqH2H2(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb102nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb102nf_half(%esp),%xmm3   ## rinvH2H2 
        mulps   nb102nf_half(%esp),%xmm7   ## rinvH2H1 
        movaps  %xmm3,nb102nf_rinvH2H2(%esp)
        movaps  %xmm7,nb102nf_rinvH2H1(%esp)

        rsqrtps nb102nf_rsqOO(%esp),%xmm1
        rsqrtps nb102nf_rsqOH1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb102nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb102nf_rsqOO(%esp),%xmm1
        mulps   nb102nf_rsqOH1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb102nf_half(%esp),%xmm3
        mulps   nb102nf_half(%esp),%xmm7
        movaps  %xmm3,nb102nf_rinvOO(%esp)
        movaps  %xmm7,nb102nf_rinvOH1(%esp)

        rsqrtps nb102nf_rsqOH2(%esp),%xmm1
        rsqrtps nb102nf_rsqH1O(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb102nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb102nf_rsqOH2(%esp),%xmm1
        mulps   nb102nf_rsqH1O(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb102nf_half(%esp),%xmm3
        mulps   nb102nf_half(%esp),%xmm7
        movaps  %xmm3,nb102nf_rinvOH2(%esp)
        movaps  %xmm7,nb102nf_rinvH1O(%esp)

        rsqrtps nb102nf_rsqH1H1(%esp),%xmm1
        rsqrtps nb102nf_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb102nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb102nf_rsqH1H1(%esp),%xmm1
        mulps   nb102nf_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb102nf_half(%esp),%xmm3
        mulps   nb102nf_half(%esp),%xmm7
        movaps  %xmm3,nb102nf_rinvH1H1(%esp)
        movaps  %xmm7,nb102nf_rinvH1H2(%esp)

        rsqrtps nb102nf_rsqH2O(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb102nf_three(%esp),%xmm3
        mulps   nb102nf_rsqH2O(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb102nf_half(%esp),%xmm3
        movaps  %xmm3,nb102nf_rinvH2O(%esp)

        ## sum OO pot in xmm0, OH in xmm1 HH in xmm2 
        movaps nb102nf_rinvOO(%esp),%xmm0
        movaps nb102nf_rinvOH1(%esp),%xmm1
        movaps nb102nf_rinvH1H1(%esp),%xmm2
        addps  nb102nf_rinvOH2(%esp),%xmm1
        addps  nb102nf_rinvH1H2(%esp),%xmm2
        addps  nb102nf_rinvH1O(%esp),%xmm1
        addps  nb102nf_rinvH2H1(%esp),%xmm2
        addps  nb102nf_rinvH2O(%esp),%xmm1
        addps  nb102nf_rinvH2H2(%esp),%xmm2

        mulps  nb102nf_qqOO(%esp),%xmm0
        mulps  nb102nf_qqOH(%esp),%xmm1
        mulps  nb102nf_qqHH(%esp),%xmm2
        addps  nb102nf_vctot(%esp),%xmm0
        addps  %xmm2,%xmm1
        addps  %xmm1,%xmm0
        movaps  %xmm0,nb102nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb102nf_innerk(%esp)
        jl    _nb_kernel102nf_ia32_sse.nb102nf_single_check
        jmp   _nb_kernel102nf_ia32_sse.nb102nf_unroll_loop
_nb_kernel102nf_ia32_sse.nb102nf_single_check: 
        addl $4,nb102nf_innerk(%esp)
        jnz   _nb_kernel102nf_ia32_sse.nb102nf_single_loop
        jmp   _nb_kernel102nf_ia32_sse.nb102nf_updateouterdata
_nb_kernel102nf_ia32_sse.nb102nf_single_loop: 
        movl  nb102nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb102nf_innerjjnr(%esp)

        movl nb102nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        xorps %xmm3,%xmm3
        xorps %xmm4,%xmm4
        xorps %xmm5,%xmm5

        movss (%esi,%eax,4),%xmm3               ## jxO  -  -  -
        movss 4(%esi,%eax,4),%xmm4              ## jyO  -  -  -
        movss 8(%esi,%eax,4),%xmm5              ## jzO  -  -  -  

        movlps 12(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1   -    -
        movss  20(%esi,%eax,4),%xmm7            ## xmm7 = jzH1   -    -    - 
        movhps 24(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1 jxH2 jyH2
        movss  32(%esi,%eax,4),%xmm2            ## xmm2 = jzH2   -    -    -

        ## have all coords, time for some shuffling.

        shufps $216,%xmm6,%xmm6 ## constant 11011000     ;# xmm6 = jxH1 jxH2 jyH1 jyH2 
        unpcklps %xmm2,%xmm7                    ## xmm7 = jzH1 jzH2   -    -
        movaps  nb102nf_ixO(%esp),%xmm0
        movaps  nb102nf_iyO(%esp),%xmm1
        movaps  nb102nf_izO(%esp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxO   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## constant 11100100    ;# xmm4 = jyO   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## constant 01000100    ;# xmm5 = jzO   0   jzH1 jzH2

        ## store all j coordinates in jO  
        movaps %xmm3,nb102nf_jxO(%esp)
        movaps %xmm4,nb102nf_jyO(%esp)
        movaps %xmm5,nb102nf_jzO(%esp)
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
        movaps  nb102nf_three(%esp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb102nf_half(%esp),%xmm3   ## rinv iO - j water 

        xorps   %xmm1,%xmm1

        xorps   %xmm4,%xmm4

        ## fetch charges to xmm4 (temporary) 
        movss   nb102nf_qqOO(%esp),%xmm4

        movhps  nb102nf_qqOH(%esp),%xmm4

        mulps   %xmm4,%xmm3     ## xmm3=vcoul 

        addps   nb102nf_vctot(%esp),%xmm3
        movaps  %xmm3,nb102nf_vctot(%esp)

        ## done with i O Now do i H1 & H2 simultaneously: 
        movaps  nb102nf_ixH1(%esp),%xmm0
        movaps  nb102nf_iyH1(%esp),%xmm1
        movaps  nb102nf_izH1(%esp),%xmm2
        movaps  nb102nf_ixH2(%esp),%xmm3
        movaps  nb102nf_iyH2(%esp),%xmm4
        movaps  nb102nf_izH2(%esp),%xmm5
        subps   nb102nf_jxO(%esp),%xmm0
        subps   nb102nf_jyO(%esp),%xmm1
        subps   nb102nf_jzO(%esp),%xmm2
        subps   nb102nf_jxO(%esp),%xmm3
        subps   nb102nf_jyO(%esp),%xmm4
        subps   nb102nf_jzO(%esp),%xmm5
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
        movaps  nb102nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb102nf_half(%esp),%xmm3   ## rinv H1 - j water 
        mulps   nb102nf_half(%esp),%xmm7   ## rinv H2 - j water  

        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        ## do coulomb interaction 
        movaps  %xmm3,%xmm0
        movss   nb102nf_qqOH(%esp),%xmm6
        movaps  %xmm7,%xmm4
        movhps  nb102nf_qqHH(%esp),%xmm6
        mulps   %xmm6,%xmm3     ## vcoul 
        mulps   %xmm6,%xmm7     ## vcoul 
        addps   %xmm7,%xmm3     ## total vcoul 
        addps   nb102nf_vctot(%esp),%xmm3
        movaps  %xmm3,nb102nf_vctot(%esp)

        decl  nb102nf_innerk(%esp)
        jz    _nb_kernel102nf_ia32_sse.nb102nf_updateouterdata
        jmp   _nb_kernel102nf_ia32_sse.nb102nf_single_loop
_nb_kernel102nf_ia32_sse.nb102nf_updateouterdata: 
        ## get n from stack
        movl nb102nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb102nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb102nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb102nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb102nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel102nf_ia32_sse.nb102nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb102nf_n(%esp)
        jmp _nb_kernel102nf_ia32_sse.nb102nf_outer
_nb_kernel102nf_ia32_sse.nb102nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb102nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel102nf_ia32_sse.nb102nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel102nf_ia32_sse.nb102nf_threadloop
_nb_kernel102nf_ia32_sse.nb102nf_end: 
        emms

        movl nb102nf_nouter(%esp),%eax
        movl nb102nf_ninner(%esp),%ebx
        movl nb102nf_outeriter(%ebp),%ecx
        movl nb102nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb102nf_salign(%esp),%eax
        addl %eax,%esp
        addl $712,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


